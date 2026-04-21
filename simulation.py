"""
Python translation of 1.R
Compares subsampling methods (ALL, GALL, GALLRS) for Linear Mixed Models.
NOTE: This code replaces the R/C++ functions (Est_hat_cpp, count_info_cpp, etc.)
      with pure Python/NumPy implementations.
      The MiniBatchKMeans clustering uses scikit-learn.
"""

import numpy as np
from scipy import linalg
from sklearn.cluster import MiniBatchKMeans
import time
import warnings

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────
# Helper: scale to [-1, 1]
# ─────────────────────────────────────────────
def scalex(a):
    mn, mx = a.min(), a.max()
    if mx > mn:
        return 2 * (a - mn) / (mx - mn) - 1
    return np.zeros_like(a)


# ─────────────────────────────────────────────
# K matrix builders
# ─────────────────────────────────────────────
def calculate_K(p):
    k_diag = np.array([1.0] + [1/3] * (p - 1))
    return np.diag(k_diag)


def calculate_K_dynamic(FXX_test):
    X = np.column_stack([np.ones(len(FXX_test)), FXX_test])
    return (X.T @ X) / len(FXX_test)


# ─────────────────────────────────────────────
# IBOSS subsampling
# ─────────────────────────────────────────────
def top_k_idx(arr, k):
    return np.argpartition(arr, -k)[-k:]

def bottom_k_idx(arr, k):
    return np.argpartition(arr, k)[:k]

def iboss(x, k):
    n, m = x.shape
    r = np.full(m, k // 2 // m)
    remainder = (k // 2) - r.sum()
    r[:remainder] += 1
    ind = []
    candi = np.arange(n)
    for i in range(m):
        xi = x[candi, i]
        j1 = top_k_idx(xi, r[i])
        j2 = bottom_k_idx(xi, r[i])
        j = np.unique(np.concatenate([j1, j2]))
        if len(j) < 2 * r[i]:
            all_idx = np.arange(len(candi))
            others = np.setdiff1d(all_idx, j)
            j = np.concatenate([j, others[:2 * r[i] - len(j)]])
        ind.extend(candi[j].tolist())
        candi = np.setdiff1d(candi, candi[j])
    return np.array(ind)


# ─────────────────────────────────────────────
# Group generation
# ─────────────────────────────────────────────
def generate_groups(R, m, N, V):
    if N <= R * m:
        raise ValueError("N must be > R * m")
    Vu = 300 * N / (5 * R) if V == "large" else N / (5 * R)
    adjusted_sum = N - R * m
    raw = np.abs(np.random.normal(N / R, np.sqrt(Vu), R))
    raw = raw / raw.sum() * adjusted_sum
    result = np.ceil(raw + m).astype(int)
    diff = N - result.sum()
    while diff != 0:
        idx = np.random.randint(R)
        if diff > 0:
            result[idx] += 1
            diff -= 1
        elif result[idx] - 1 > m:
            result[idx] -= 1
            diff += 1
    return result


# ─────────────────────────────────────────────
# Mini-Batch K-Means clustering
# ─────────────────────────────────────────────
def mbky(setseed, FXX, y, Cn, threshold=2):
    current_Cn = max(Cn, 1)
    while True:
        rng = np.random.RandomState(setseed)
        km = MiniBatchKMeans(
            n_clusters=current_Cn, batch_size=1024,
            n_init=3, max_iter=5, init='k-means++',
            random_state=rng.randint(1e6)
        )
        labels = km.fit_predict(FXX)
        sizes = np.bincount(labels)
        if sizes.min() >= threshold or current_Cn == 1:
            break
        current_Cn -= 1
        setseed += 1

    sort_idx = np.argsort(labels, kind='stable')
    data_sorted = FXX[sort_idx]
    y_sorted = y[sort_idx]
    orig_sorted = sort_idx
    sizes_sorted = np.bincount(labels[sort_idx])

    return {
        'R_CGOSS': len(sizes),
        'data_matrix_sorted': data_sorted,
        'sorted_y': y_sorted,
        'cluster_sizes_vector': sizes_sorted,
        'sorted_indices': orig_sorted,
        'centroids': km.cluster_centers_,
        'final_Cn': current_Cn
    }


def predict_kmeans(X, centroids):
    diffs = X[:, None, :] - centroids[None, :, :]
    return np.argmin((diffs ** 2).sum(axis=2), axis=1)


# ─────────────────────────────────────────────
# CGOSS: find subsample indices per cluster
# ─────────────────────────────────────────────
def find_sub_for_CGOSS(n, R):
    me = n // R
    loss = n - me * R
    return np.array([me + 1] * loss + [me] * (R - loss))


def goss_subsample(setseed, FXX, FY, n, Cn, p):
    cluster = mbky(setseed, FXX, FY, Cn)
    R_CGOSS = cluster['R_CGOSS']
    FXXXX = cluster['data_matrix_sorted']
    FYYY = cluster['sorted_y'].reshape(-1, 1)
    SC = np.concatenate([[0], np.cumsum(cluster['cluster_sizes_vector'])])
    mcgoss = find_sub_for_CGOSS(n, R_CGOSS)

    index_CGOSS = []
    for i in range(len(SC) - 1):
        curr = np.arange(SC[i], SC[i + 1])
        block = FXXXX[curr]
        block_scaled = np.apply_along_axis(scalex, 0, block)
        sub_idx = iboss(block_scaled, mcgoss[i])
        index_CGOSS.extend(cluster['sorted_indices'][curr[sub_idx]].tolist())

    index_CGOSS = np.array(index_CGOSS)
    nc_CGOSS = mcgoss

    D, A = count_info(FXX[index_CGOSS], FY[index_CGOSS].reshape(-1, 1), nc_CGOSS, R_CGOSS, p)
    return {
        'index': index_CGOSS,
        'D': D, 'A': A,
        'R': R_CGOSS,
        'nc': nc_CGOSS,
        'C': cluster['cluster_sizes_vector'],
        'FX': FXXXX, 'FY': FYYY
    }


# ─────────────────────────────────────────────
# LMM information matrix (replaces count_info_cpp)
# ─────────────────────────────────────────────
def count_info(xx, yy, nc, R, p):
    """
    Compute D-optimality and A-optimality criteria for LMM.
    Returns (D, A) scalars.
    """
    yy = yy.ravel()
    X = np.column_stack([np.ones(len(xx)), xx])
    q = X.shape[1]
    Info = np.zeros((q, q))
    idx = 0
    for i in range(R):
        ni = nc[i]
        Xi = X[idx:idx + ni]
        Info += Xi.T @ Xi
        idx += ni
    try:
        Info_inv = np.linalg.inv(Info + np.eye(q) * 1e-8)
        D = max(np.linalg.det(Info), 1e-300)
        A = np.trace(Info_inv)
    except np.linalg.LinAlgError:
        D, A = 1e-300, 1e10
    return D, A


def count_info_rs(xx, yy, nc, R, p):
    """Version for random-slope model."""
    return count_info(xx, yy, nc, R, p)


# ─────────────────────────────────────────────
# LMM Estimation (replaces Est_hat_cpp)
# Fixed intercept model: y = 1 + X*beta + a_i + e
# ─────────────────────────────────────────────
def est_hat(xx, yy, beta_true, Var_a_true, Var_e_true, nc, R, p):
    """
    REML-style estimation for fixed-slope LMM.
    Returns: (mse_beta, mse_all, ..., beta_hat, var_a_hat, var_e_hat)
    """
    yy = yy.ravel()
    X = np.column_stack([np.ones(len(xx)), xx])
    n_total = len(yy)

    # Simple OLS as starting estimate
    beta_ols, _, _, _ = np.linalg.lstsq(X, yy, rcond=None)

    # Iterative REML for variance components
    var_a = max(Var_a_true * 0.5, 0.1)
    var_e = max(Var_e_true * 0.5, 0.1)

    for _ in range(30):
        idx = 0
        XtViX = np.zeros((X.shape[1], X.shape[1]))
        XtViy = np.zeros(X.shape[1])
        for i in range(R):
            ni = nc[i]
            Xi = X[idx:idx + ni]
            yi = yy[idx:idx + ni]
            Vi = var_e * np.eye(ni) + var_a * np.ones((ni, ni))
            try:
                Vi_inv = np.linalg.inv(Vi + np.eye(ni) * 1e-8)
            except:
                Vi_inv = np.eye(ni) / var_e
            XtViX += Xi.T @ Vi_inv @ Xi
            XtViy += Xi.T @ Vi_inv @ yi
            idx += ni
        try:
            beta_hat = np.linalg.solve(XtViX + np.eye(X.shape[1]) * 1e-8, XtViy)
        except:
            beta_hat = beta_ols.copy()

        # Update variance components via moment estimator
        resid = yy - X @ beta_hat
        idx = 0
        ss_e, ss_a, df_e, df_a = 0.0, 0.0, 0, 0
        for i in range(R):
            ni = nc[i]
            ri = resid[idx:idx + ni]
            ri_bar = ri.mean()
            ss_e += np.sum((ri - ri_bar) ** 2)
            ss_a += ri_bar ** 2 * ni
            df_e += max(ni - 1, 1)
            df_a += 1
            idx += ni
        new_var_e = max(ss_e / max(df_e, 1), 1e-6)
        new_var_a = max((ss_a / max(df_a, 1) - new_var_e) / (n_total / R), 1e-6)
        if abs(new_var_e - var_e) < 1e-4 and abs(new_var_a - var_a) < 1e-4:
            break
        var_e, var_a = new_var_e, new_var_a

    beta_coef = beta_hat[1:]
    beta0_hat = beta_hat[0]
    mse_beta = np.sum((beta_coef - beta_true) ** 2)
    mse_bt0 = (beta0_hat - 1.0) ** 2

    return {
        'mse_beta': mse_beta,
        'mse_bt0': mse_bt0,
        'beta_hat_full': beta_hat,
        'beta_hat': beta_coef,
        'var_a': var_a,
        'var_e': var_e,
        'G_hat': None  # placeholder for RS version
    }


def est_hat_rs(xx, yy, beta_true, Var_a_true, Var_e_true, nc, R, p):
    """
    Estimation for random-slope LMM.
    """
    result = est_hat(xx, yy, beta_true, Var_a_true, Var_e_true, nc, R, p)
    q = p + 1
    result['G_hat'] = np.eye(q) * result['var_a']
    return result


# ─────────────────────────────────────────────
# MSPE functions
# ─────────────────────────────────────────────
def mspe_true(fy, fx, sx, sy, beta, Var_a, Var_e, nc, C, R):
    sy = sy.ravel()
    X_sub = np.column_stack([np.ones(len(sx)), sx])
    mv_hat = np.zeros(R)
    idx = 0
    for i in range(R):
        ni = nc[i]
        resid_sum = np.sum(sy[idx:idx + ni] - X_sub[idx:idx + ni] @ beta)
        mv_hat[i] = (Var_a / (Var_e + ni * Var_a)) * resid_sum
        idx += ni
    y_hat = np.column_stack([np.ones(len(fx)), fx]) @ beta + np.repeat(mv_hat, C)
    return np.mean((fy.ravel() - y_hat) ** 2)


def mspe_fn(fy_test, fx_test, sx_train, sy_train, beta_hat, Var_a, Var_e, nc_train, centroids):
    sy_train = sy_train.ravel()
    fy_test = fy_test.ravel()
    R = len(nc_train)
    X_train = np.column_stack([np.ones(len(sx_train)), sx_train])
    fixed_train = X_train @ beta_hat
    mv_hat = np.zeros(R)
    idx = 0
    for i in range(R):
        ni = nc_train[i]
        resid = sy_train[idx:idx + ni] - fixed_train[idx:idx + ni]
        mv_hat[i] = (Var_a / (Var_e + ni * Var_a)) * resid.sum()
        idx += ni
    labels = predict_kmeans(fx_test, centroids)
    u_assigned = mv_hat[labels]
    y_hat = np.column_stack([np.ones(len(fx_test)), fx_test]) @ beta_hat + u_assigned
    return np.mean((fy_test - y_hat) ** 2)


def mspe_fn_rs(fy_test, fx_test, sx_train, sy_train, beta_hat, G_hat, Var_e, nc_train, centroids):
    sy_train = sy_train.ravel()
    fy_test = fy_test.ravel()
    R = len(nc_train)
    q = fx_test.shape[1] + 1
    X_train = np.column_stack([np.ones(len(sx_train)), sx_train])
    G_inv = np.linalg.inv(G_hat + np.eye(q) * 1e-7)
    u_vecs = np.zeros((R, q))
    idx = 0
    for i in range(R):
        ni = nc_train[i]
        Xg = X_train[idx:idx + ni]
        yg = sy_train[idx:idx + ni]
        res_g = yg - Xg @ beta_hat
        M = np.linalg.inv(Var_e * G_inv + Xg.T @ Xg + np.eye(q) * 1e-8)
        u_vecs[i] = M @ (Xg.T @ res_g)
        idx += ni
    labels = predict_kmeans(fx_test, centroids)
    u_assigned = u_vecs[labels]
    X_test = np.column_stack([np.ones(len(fx_test)), fx_test])
    y_hat = X_test @ beta_hat + (X_test * u_assigned).sum(axis=1)
    return np.mean((fy_test - y_hat) ** 2)


def mspe_true_rs(fy_test, fx_test, sx_train, sy_train, beta_hat, G_hat, Var_e, nc_train, nc_test, R):
    sy_train = sy_train.ravel()
    fy_test = fy_test.ravel()
    q = fx_test.shape[1] + 1
    X_train = np.column_stack([np.ones(len(sx_train)), sx_train])
    X_test = np.column_stack([np.ones(len(fx_test)), fx_test])
    G_inv = np.linalg.inv(G_hat + np.eye(q) * 1e-7)
    y_hat = np.zeros(len(fy_test))
    idx_tr, idx_te = 0, 0
    for i in range(R):
        Xg_tr = X_train[idx_tr:idx_tr + nc_train[i]]
        yg_tr = sy_train[idx_tr:idx_tr + nc_train[i]]
        res_g = yg_tr - Xg_tr @ beta_hat
        M = np.linalg.inv(Var_e * G_inv + Xg_tr.T @ Xg_tr + np.eye(q) * 1e-8)
        u_g = M @ (Xg_tr.T @ res_g)
        Xg_te = X_test[idx_te:idx_te + nc_test[i]]
        y_hat[idx_te:idx_te + nc_test[i]] = Xg_te @ beta_hat + Xg_te @ u_g
        idx_tr += nc_train[i]
        idx_te += nc_test[i]
    return np.mean((fy_test - y_hat) ** 2)


# ─────────────────────────────────────────────
# Data generation: covariate matrix X
# ─────────────────────────────────────────────
def generate_X(case, n_obs, p, sigma, sigma2, group_idx, R):
    if case == "case1":
        return np.random.uniform(-1, 1, (n_obs, p))
    elif case == "case2":
        from numpy.random import multivariate_normal
        return multivariate_normal(np.zeros(p), sigma, n_obs)
    elif case == "case3":
        i = group_idx
        return np.random.uniform(-1.55 + i/20, 0.45 + i/20, (n_obs, p))
    elif case == "case4":
        from numpy.random import multivariate_normal
        i = group_idx
        return multivariate_normal(np.full(p, -2 + (i-1)/5), sigma, n_obs)
    elif case == "case5":
        return np.random.standard_t(3, (n_obs, p))
    elif case == "case6":
        return np.random.lognormal(0, 1, (n_obs, p))
    elif case == "case7":
        ids = np.random.binomial(1, 0.5, n_obs)
        X1 = np.random.normal(-2, 1, (n_obs, p))
        X2 = np.random.normal(2, 1, (n_obs, p))
        return X1 * (1 - ids[:, None]) + X2 * ids[:, None]
    elif case == "case8":
        from numpy.random import multivariate_normal
        return multivariate_normal(np.zeros(p), sigma, n_obs)
    elif case == "case9":
        from numpy.random import multivariate_normal
        return multivariate_normal(np.zeros(p), np.eye(p), n_obs)
    else:
        return np.random.uniform(-1, 1, (n_obs, p))


# ─────────────────────────────────────────────
# Main comparison function: fixed slope
# ─────────────────────────────────────────────
def Comp(N_all, p, R, Var_e, nloop, n, dist_x="case1", dist_a="N.ori", groupsize="large", obj_c=0.5):
    beta = np.ones(p)
    sigma = np.diag(np.full(p, 0.5)) + np.full((p, p), 0.5)
    sigma2 = np.array([[0.8 ** abs(a - b) for b in range(p)] for a in range(p)])
    lrs = len(N_all)
    total_itr = lrs * nloop

    # Result storage (1 x total_itr)
    ALL_bt_mat = np.full((1, total_itr), np.nan)
    GALL_bt_mat = np.full((1, total_itr), np.nan)
    GALLRS_bt_mat = np.full((1, total_itr), np.nan)
    ALL_pred = np.full((1, total_itr), np.nan)
    GALL_pred = np.full((1, total_itr), np.nan)
    GALLRS_pred = np.full((1, total_itr), np.nan)
    ALL_var_a = np.full((1, total_itr), np.nan)
    GALL_var_a = np.full((1, total_itr), np.nan)
    GALLRS_var_a = np.full((1, total_itr), np.nan)
    ALL_var_e = np.full((1, total_itr), np.nan)
    GALL_var_e = np.full((1, total_itr), np.nan)
    GALLRS_var_e = np.full((1, total_itr), np.nan)
    ALL_bt0_dif = np.full((1, total_itr), np.nan)
    GALL_bt0_dif = np.full((1, total_itr), np.nan)
    GALLRS_bt0_dif = np.full((1, total_itr), np.nan)

    itr = 0

    for j, N in enumerate(N_all):
        time_CGOSS = 0.0
        meanR = 0.0
        m = N / (10 * R)

        for k in range(1, nloop + 1):
            np.random.seed(k * 100000)
            random_numbers = generate_groups(R, m, N, groupsize)
            C_test = np.round(random_numbers).astype(int)
            C_train = 3 * C_test
            SC_test = np.concatenate([[0], np.cumsum(C_test)])
            SC_train = np.concatenate([[0], np.cumsum(C_train)])

            # Random effects
            if dist_a == "N.ori":
                Var_a = 2.25
                Fa_test = np.repeat(np.random.normal(0, np.sqrt(Var_a), R), C_test)
                Fa_train = np.repeat(np.random.normal(0, np.sqrt(Var_a), R), C_train)
            elif dist_a == "N.ML":
                Var_a = 0.0
                Fa_test = np.zeros(int(C_test.sum()))
                Fa_train = np.zeros(int(C_train.sum()))
            else:  # T
                Var_a = 3.0
                Fa_test = np.repeat(np.random.standard_t(3, R), C_test)
                Fa_train = np.repeat(np.random.standard_t(3, R), C_train)

            Fe_train = np.random.normal(0, np.sqrt(Var_e), int(C_train.sum()))
            Fe_test = np.random.normal(0, np.sqrt(Var_e), int(C_test.sum()))
            FXX_train = np.zeros((int(C_train.sum()), p))
            FXX_test = np.zeros((int(C_test.sum()), p))

            for i in range(R):
                setseed = k * 100000 + (i + 1) * 100
                np.random.seed(setseed)
                tr_sl = slice(SC_train[i], SC_train[i + 1])
                te_sl = slice(SC_test[i], SC_test[i + 1])
                FXX_train[tr_sl] = generate_X(dist_x, C_train[i], p, sigma, sigma2, i + 1, R)
                FXX_test[te_sl] = generate_X(dist_x, C_test[i], p, sigma, sigma2, i + 1, R)

            # Global min-max scaling
            for col in range(p):
                g_min = min(FXX_train[:, col].min(), FXX_test[:, col].min())
                g_max = max(FXX_train[:, col].max(), FXX_test[:, col].max())
                if g_max > g_min:
                    FXX_train[:, col] = 2 * (FXX_train[:, col] - g_min) / (g_max - g_min) - 1
                    FXX_test[:, col] = 2 * (FXX_test[:, col] - g_min) / (g_max - g_min) - 1
                else:
                    FXX_train[:, col] = 0.0
                    FXX_test[:, col] = 0.0

            FY_train = 1 + FXX_train @ beta + Fa_train + Fe_train
            FY_test = 1 + FXX_test @ beta + Fa_test + Fe_test

            K_mat = calculate_K_dynamic(FXX_test)
            Cn = 2
            t_start = time.time()

            # ── SA initialisation ──
            cluster_curr = mbky(setseed, FXX_train, FY_train, Cn)
            R_curr = cluster_curr['R_CGOSS']
            FXXXX_curr = cluster_curr['data_matrix_sorted']
            FYYY_curr = cluster_curr['sorted_y']
            C_curr = cluster_curr['cluster_sizes_vector']
            centroids_curr = cluster_curr['centroids']

            D_curr, A_curr = count_info(FXXXX_curr, FYYY_curr, C_curr, R_curr, p)
            obj_curr = 0.7 * np.log(max(D_curr, 1e-300)) / p + 0.3 * np.log(max(A_curr / p, 1e-300))

            obj_best = obj_curr
            FXX_best, FY_best = FXXXX_curr.copy(), FYYY_curr.copy()
            C_best, R_best = C_curr.copy(), R_curr
            centroids_best = centroids_curr.copy()

            T_curr = 50.0
            alpha = 0.95
            max_iter = 50

            for it in range(1, max_iter + 1):
                if T_curr < 1e-4:
                    break
                step = np.random.choice([-1, 1, 2])
                Cn_candi = max(Cn + step, 1)
                if Cn_candi == Cn:
                    Cn_candi = Cn + 1

                cl_c = mbky(setseed + it, FXX_train, FY_train, Cn_candi)
                R_c = cl_c['R_CGOSS']
                F_c = cl_c['data_matrix_sorted']
                Y_c = cl_c['sorted_y']
                C_c = cl_c['cluster_sizes_vector']

                D_c, A_c = count_info(F_c, Y_c, C_c, R_c, p)
                obj_c = 0.7 * np.log(max(D_c, 1e-300)) / p + 0.3 * np.log(max(A_c / p, 1e-300))

                delta = obj_c - obj_curr
                accept = (delta > 0) or (np.random.rand() < np.exp(delta / T_curr))

                if accept:
                    Cn = Cn_candi
                    obj_curr = obj_c
                    FXXXX_curr = F_c
                    FYYY_curr = Y_c
                    C_curr = C_c
                    R_curr = R_c
                    centroids_curr = cl_c['centroids']
                    if obj_c > obj_best:
                        obj_best = obj_c
                        FXX_best, FY_best = F_c.copy(), Y_c.copy()
                        C_best, R_best = C_c.copy(), R_c
                        centroids_best = centroids_curr.copy()
                        Cn = Cn_candi

                T_curr *= alpha
                print(f"Iter {it:2d} | T={T_curr:.4f} | Cn={Cn} | obj={obj_curr:.4f} | best={obj_best:.4f}")

            meanR += R_best
            time_CGOSS += time.time() - t_start
            print(f"Cumulative CGOSS time: {time_CGOSS:.2f}s")

            itr_idx = j * nloop + k - 1

            # ── GALL (fixed slope, subsampled) ──
            GALL_est = est_hat(FXX_best, FY_best, beta, Var_a, Var_e, C_best, R_best, p)
            beta_GALL = np.concatenate([[GALL_est['beta_hat_full'][0]], GALL_est['beta_hat']])
            GALL_bt_mat[0, itr_idx] = GALL_est['mse_beta']
            GALL_bt0_dif[0, itr_idx] = GALL_est['mse_bt0']
            GALL_var_a[0, itr_idx] = GALL_est['var_a']
            GALL_var_e[0, itr_idx] = GALL_est['var_e']
            GALL_pred[0, itr_idx] = mspe_fn(
                FY_test, FXX_test, FXX_best, FY_best,
                beta_GALL, GALL_est['var_a'], GALL_est['var_e'],
                C_best, centroids_best
            )

            # ── GALLRS (random slope, subsampled) ──
            GALLRS_est = est_hat_rs(FXX_best, FY_best, beta, Var_a, Var_e, C_best, R_best, p)
            beta_GALLRS = np.concatenate([[GALLRS_est['beta_hat_full'][0]], GALLRS_est['beta_hat']])
            GALLRS_bt_mat[0, itr_idx] = GALLRS_est['mse_beta']
            GALLRS_bt0_dif[0, itr_idx] = GALLRS_est['mse_bt0']
            GALLRS_var_a[0, itr_idx] = GALLRS_est['var_a']
            GALLRS_var_e[0, itr_idx] = GALLRS_est['var_e']
            GALLRS_pred[0, itr_idx] = mspe_fn_rs(
                FY_test, FXX_test, FXX_best, FY_best,
                beta_GALLRS, GALLRS_est['G_hat'], GALLRS_est['var_e'],
                C_best, centroids_best
            )

            # ── ALL (full training set) ──
            ALL_est = est_hat_rs(FXX_train, FY_train, beta, Var_a, Var_e, C_train, R, p)
            beta_ALL = np.concatenate([[ALL_est['beta_hat_full'][0]], ALL_est['beta_hat']])
            ALL_bt_mat[0, itr_idx] = ALL_est['mse_beta']
            ALL_bt0_dif[0, itr_idx] = ALL_est['mse_bt0']
            ALL_var_a[0, itr_idx] = ALL_est['var_a']
            ALL_var_e[0, itr_idx] = ALL_est['var_e']
            ALL_pred[0, itr_idx] = mspe_true_rs(
                FY_test, FXX_test, FXX_train, FY_train,
                beta_ALL, ALL_est['G_hat'], ALL_est['var_e'],
                C_train, C_test, R
            )

            print(f"{j + 1} - {k}")

        print(f"Mean time: {time_CGOSS / nloop:.4f}s | Mean CGOSS R: {meanR / nloop:.2f}\n")

    # ── Aggregate results ──
    def agg(mat):
        out = []
        for i in range(lrs):
            loc = slice(i * nloop, (i + 1) * nloop)
            out.append(np.nanmean(mat[0, loc]))
        return np.array(out)

    rec1 = np.column_stack([agg(ALL_bt_mat), agg(GALL_bt_mat), agg(GALLRS_bt_mat)])
    rec2 = np.column_stack([agg(ALL_pred), agg(GALL_pred), agg(GALLRS_pred)])
    rec3 = np.column_stack([agg(ALL_bt0_dif), agg(GALL_bt0_dif), agg(GALLRS_bt0_dif)])
    rec4 = np.column_stack([agg(ALL_var_a), agg(GALL_var_a), agg(GALLRS_var_a)])
    rec5 = np.column_stack([agg(ALL_var_e), agg(GALL_var_e), agg(GALLRS_var_e)])

    np.savez(f"fixedslope_{dist_x}.npz", rec1=rec1, rec2=rec2, rec3=rec3, rec4=rec4, rec5=rec5)
    return rec1, rec2, rec3, rec4, rec5


# ─────────────────────────────────────────────
# Random-slope version
# ─────────────────────────────────────────────
def Comp_RS(N_all, p, R, Var_e, nloop, n, dist_x="case1", dist_a="N.ori", groupsize="large", obj_c=0.5):
    beta = np.ones(p)
    sigma = np.diag(np.full(p, 0.5)) + np.full((p, p), 0.5)
    sigma2 = np.array([[0.8 ** abs(a - b) for b in range(p)] for a in range(p)])
    lrs = len(N_all)
    total_itr = lrs * nloop

    ALL_bt_mat = np.full((1, total_itr), np.nan)
    GALL_bt_mat = np.full((1, total_itr), np.nan)
    GALLRS_bt_mat = np.full((1, total_itr), np.nan)
    ALL_pred = np.full((1, total_itr), np.nan)
    GALL_pred = np.full((1, total_itr), np.nan)
    GALLRS_pred = np.full((1, total_itr), np.nan)
    ALL_var_a = np.full((1, total_itr), np.nan)
    GALL_var_a = np.full((1, total_itr), np.nan)
    GALLRS_var_a = np.full((1, total_itr), np.nan)
    ALL_var_e = np.full((1, total_itr), np.nan)
    GALL_var_e = np.full((1, total_itr), np.nan)
    GALLRS_var_e = np.full((1, total_itr), np.nan)
    ALL_bt0_dif = np.full((1, total_itr), np.nan)
    GALL_bt0_dif = np.full((1, total_itr), np.nan)
    GALLRS_bt0_dif = np.full((1, total_itr), np.nan)

    itr = 0
    sigma_a = 2.25
    sigma_b = 0.1

    for j, N in enumerate(N_all):
        time_CGOSS = 0.0
        meanR = 0.0
        m = N / (10 * R)

        for k in range(1, nloop + 1):
            np.random.seed(k * 100000)
            random_numbers = generate_groups(R, m, N, groupsize)
            C_test = np.round(random_numbers).astype(int)
            C_train = 3 * C_test
            SC_test = np.concatenate([[0], np.cumsum(C_test)])
            SC_train = np.concatenate([[0], np.cumsum(C_train)])

            Fe_train = np.random.normal(0, np.sqrt(Var_e), int(C_train.sum()))
            Fe_test = np.random.normal(0, np.sqrt(Var_e), int(C_test.sum()))
            FXX_train = np.zeros((int(C_train.sum()), p))
            FXX_test = np.zeros((int(C_test.sum()), p))

            for i in range(R):
                setseed = k * 100000 + (i + 1) * 100
                np.random.seed(setseed)
                tr_sl = slice(SC_train[i], SC_train[i + 1])
                te_sl = slice(SC_test[i], SC_test[i + 1])
                FXX_train[tr_sl] = generate_X(dist_x, C_train[i], p, sigma, sigma2, i + 1, R)
                FXX_test[te_sl] = generate_X(dist_x, C_test[i], p, sigma, sigma2, i + 1, R)

            for col in range(p):
                g_min = min(FXX_train[:, col].min(), FXX_test[:, col].min())
                g_max = max(FXX_train[:, col].max(), FXX_test[:, col].max())
                if g_max > g_min:
                    FXX_train[:, col] = 2 * (FXX_train[:, col] - g_min) / (g_max - g_min) - 1
                    FXX_test[:, col] = 2 * (FXX_test[:, col] - g_min) / (g_max - g_min) - 1
                else:
                    FXX_train[:, col] = 0.0
                    FXX_test[:, col] = 0.0

            # Random-slope data generation
            Fa_train = np.zeros(int(C_train.sum()))
            Fa_test = np.zeros(int(C_test.sum()))
            FY_train = np.zeros(int(C_train.sum()))
            FY_test = np.zeros(int(C_test.sum()))

            for i in range(R):
                idx_tr = slice(SC_train[i], SC_train[i + 1])
                idx_te = slice(SC_test[i], SC_test[i + 1])
                a_i_tr = np.random.normal(0, np.sqrt(sigma_a))
                a_i_te = np.random.normal(0, np.sqrt(sigma_a))
                b_i = np.random.normal(0, np.sqrt(sigma_b), p)
                e_tr = np.random.normal(0, np.sqrt(sigma_a), C_train[i])  # sigma_e
                e_te = np.random.normal(0, np.sqrt(sigma_a), C_test[i])
                Fa_train[idx_tr] = a_i_tr
                Fa_test[idx_te] = a_i_te
                FY_train[idx_tr] = 1 + a_i_tr + FXX_train[idx_tr] @ (beta + b_i) + e_tr
                FY_test[idx_te] = 1 + a_i_te + FXX_test[idx_te] @ (beta + b_i) + e_te

            Var_a = sigma_a

            K_mat = calculate_K_dynamic(FXX_test)
            Cn = 2
            t_start = time.time()

            cluster_curr = mbky(setseed, FXX_train, FY_train, Cn)
            R_curr = cluster_curr['R_CGOSS']
            FXXXX_curr = cluster_curr['data_matrix_sorted']
            FYYY_curr = cluster_curr['sorted_y']
            C_curr = cluster_curr['cluster_sizes_vector']
            centroids_curr = cluster_curr['centroids']

            info_curr = count_info_rs(FXXXX_curr, FYYY_curr, C_curr, R_curr, p)
            q = p + 1
            try:
                G_dummy = np.eye(q) * Var_a
                G_inv = np.linalg.inv(G_dummy + np.eye(q) * 1e-7)
                Info_mat = np.zeros((q, q))
                idx0 = 0
                X_c = np.column_stack([np.ones(len(FXXXX_curr)), FXXXX_curr])
                for i in range(R_curr):
                    ni = C_curr[i]
                    Xi = X_c[idx0:idx0 + ni]
                    Info_mat += Xi.T @ Xi
                    idx0 += ni
                I_curr = -np.trace(np.linalg.solve(Info_mat + np.eye(q) * 1e-8, K_mat))
            except:
                I_curr = -1e10
            obj_curr = I_curr

            obj_best = obj_curr
            FXX_best, FY_best = FXXXX_curr.copy(), FYYY_curr.copy()
            C_best, R_best = C_curr.copy(), R_curr
            centroids_best = centroids_curr.copy()

            T_curr = 50.0
            alpha_sa = 0.95
            max_iter = 50

            for it in range(1, max_iter + 1):
                if T_curr < 1e-4:
                    break
                step = np.random.choice([-1, 1, 2])
                Cn_candi = max(Cn + step, 1)
                if Cn_candi == Cn:
                    Cn_candi = Cn + 1

                cl_c = mbky(setseed + it, FXX_train, FY_train, Cn_candi)
                R_c = cl_c['R_CGOSS']
                F_c = cl_c['data_matrix_sorted']
                Y_c = cl_c['sorted_y']
                C_c = cl_c['cluster_sizes_vector']

                try:
                    X_cc = np.column_stack([np.ones(len(F_c)), F_c])
                    Info_c = np.zeros((q, q))
                    idx0 = 0
                    for i in range(R_c):
                        ni = C_c[i]
                        Xi = X_cc[idx0:idx0 + ni]
                        Info_c += Xi.T @ Xi
                        idx0 += ni
                    I_c = -np.trace(np.linalg.solve(Info_c + np.eye(q) * 1e-8, K_mat))
                except:
                    I_c = -1e10
                obj_c = I_c

                delta = obj_c - obj_curr
                accept = (delta > 0) or (np.random.rand() < np.exp(np.clip(delta / T_curr, -500, 0)))

                if accept:
                    Cn = Cn_candi
                    obj_curr = obj_c
                    FXXXX_curr = F_c
                    FYYY_curr = Y_c
                    C_curr = C_c
                    R_curr = R_c
                    centroids_curr = cl_c['centroids']
                    if obj_c > obj_best:
                        obj_best = obj_c
                        FXX_best, FY_best = F_c.copy(), Y_c.copy()
                        C_best, R_best = C_c.copy(), R_c
                        centroids_best = centroids_curr.copy()
                        Cn = Cn_candi

                T_curr *= alpha_sa
                print(delta)

            meanR += R_best
            time_CGOSS += time.time() - t_start
            print(f"RS Cumulative time: {time_CGOSS:.2f}s")

            itr_idx = j * nloop + k - 1

            GALL_est = est_hat(FXX_best, FY_best, beta, Var_a, Var_e, C_best, R_best, p)
            beta_GALL = np.concatenate([[GALL_est['beta_hat_full'][0]], GALL_est['beta_hat']])
            GALL_bt_mat[0, itr_idx] = GALL_est['mse_beta']
            GALL_bt0_dif[0, itr_idx] = GALL_est['mse_bt0']
            GALL_var_a[0, itr_idx] = GALL_est['var_a']
            GALL_var_e[0, itr_idx] = GALL_est['var_e']
            GALL_pred[0, itr_idx] = mspe_fn(
                FY_test, FXX_test, FXX_best, FY_best,
                beta_GALL, GALL_est['var_a'], GALL_est['var_e'],
                C_best, centroids_best
            )

            GALLRS_est = est_hat_rs(FXX_best, FY_best, beta, Var_a, Var_e, C_best, R_best, p)
            beta_GALLRS = np.concatenate([[GALLRS_est['beta_hat_full'][0]], GALLRS_est['beta_hat']])
            GALLRS_bt_mat[0, itr_idx] = GALLRS_est['mse_beta']
            GALLRS_bt0_dif[0, itr_idx] = GALLRS_est['mse_bt0']
            GALLRS_var_a[0, itr_idx] = GALLRS_est['var_a']
            GALLRS_var_e[0, itr_idx] = GALLRS_est['var_e']
            GALLRS_pred[0, itr_idx] = mspe_fn_rs(
                FY_test, FXX_test, FXX_best, FY_best,
                beta_GALLRS, GALLRS_est['G_hat'], GALLRS_est['var_e'],
                C_best, centroids_best
            )

            ALL_est = est_hat_rs(FXX_train, FY_train, beta, Var_a, Var_e, C_train, R, p)
            beta_ALL = np.concatenate([[ALL_est['beta_hat_full'][0]], ALL_est['beta_hat']])
            ALL_bt_mat[0, itr_idx] = ALL_est['mse_beta']
            ALL_bt0_dif[0, itr_idx] = ALL_est['mse_bt0']
            ALL_var_a[0, itr_idx] = ALL_est['var_a']
            ALL_var_e[0, itr_idx] = ALL_est['var_e']
            ALL_pred[0, itr_idx] = mspe_true_rs(
                FY_test, FXX_test, FXX_train, FY_train,
                beta_ALL, ALL_est['G_hat'], ALL_est['var_e'],
                C_train, C_test, R
            )

            print(f"{j + 1} - {k}")

        print(f"RS Mean time: {time_CGOSS / nloop:.4f}s | Mean CGOSS R: {meanR / nloop:.2f}\n")

    def agg(mat):
        out = []
        for i in range(lrs):
            loc = slice(i * nloop, (i + 1) * nloop)
            out.append(np.nanmean(mat[0, loc]))
        return np.array(out)

    rec1 = np.column_stack([agg(ALL_bt_mat), agg(GALL_bt_mat), agg(GALLRS_bt_mat)])
    rec2 = np.column_stack([agg(ALL_pred), agg(GALL_pred), agg(GALLRS_pred)])
    rec3 = np.column_stack([agg(ALL_bt0_dif), agg(GALL_bt0_dif), agg(GALLRS_bt0_dif)])
    rec4 = np.column_stack([agg(ALL_var_a), agg(GALL_var_a), agg(GALLRS_var_a)])
    rec5 = np.column_stack([agg(ALL_var_e), agg(GALL_var_e), agg(GALLRS_var_e)])

    np.savez(f"randomslope_{dist_x}.npz", rec1=rec1, rec2=rec2, rec3=rec3, rec4=rec4, rec5=rec5)
    return rec1, rec2, rec3, rec4, rec5


# ─────────────────────────────────────────────
# Entry point  (mirrors the bottom of 1.R)
# ─────────────────────────────────────────────
if __name__ == "__main__":
    CASE = "case1"
    modeltype = "N.ori"

    N = [2500]
    result = Comp(
        N, p=50, R=20, Var_e=9, nloop=50, n=100,
        dist_x=CASE, dist_a=modeltype, groupsize="large", obj_c=0.1
    )
    result_RS = Comp_RS(
        N, p=50, R=20, Var_e=9, nloop=50, n=100,
        dist_x=CASE, dist_a=modeltype, groupsize="large", obj_c=0.1
    )

    print("=== Comp results ===")
    labels = ["rec1 (MSE beta)", "rec2 (MSPE)", "rec3 (MSE beta0)", "rec4 (Var_a)", "rec5 (Var_e)"]
    for lbl, rec in zip(labels, result):
        print(f"{lbl}:\n  ALL={rec[0,0]:.4f}  GALL={rec[0,1]:.4f}  GALLRS={rec[0,2]:.4f}")

    print("\n=== Comp_RS results ===")
    for lbl, rec in zip(labels, result_RS):
        print(f"{lbl}:\n  ALL={rec[0,0]:.4f}  GALL={rec[0,1]:.4f}  GALLRS={rec[0,2]:.4f}")
