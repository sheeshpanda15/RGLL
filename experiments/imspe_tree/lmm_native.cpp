// ML-EM on per-leaf sufficient statistics. See lmm.py and RGLL/lmm_fast_aligned.cpp.
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#if defined(_WIN32)
#define API extern "C" __declspec(dllexport)
#else
#define API extern "C" __attribute__((visibility("default")))
#endif

using Chol = void (*)(char*, int*, double*, int*, int*);
using Solve = void (*)(char*, int*, int*, double*, int*, double*, int*, int*);
using Gemm = void (*)(char*, char*, int*, int*, int*, double*, double*, int*,
                      double*, int*, double*, double*, int*);

constexpr double epsilon = 1e-8;

double dot(const double* a, const double* b, int n) {
    return std::inner_product(a, a + n, b, 0.0);
}

bool inverse_spd(double* matrix, int n, Chol factor, Chol invert) {
    char lower = 'L';
    int info = 0;
    factor(&lower, &n, matrix, &n, &info);
    if (info != 0) return false;
    invert(&lower, &n, matrix, &n, &info);
    if (info != 0) return false;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) matrix[i * n + j] = matrix[j * n + i];
    }
    return true;
}

API int lmm_fit(const double* gram, const double* weights, const double* beta0, double ve0,
                int q, int n, int k, int rank, double tolerance, int max_iterations,
                double* beta, double* effects, double* output,
                Chol factor, Chol invert, Solve solve, Gemm gemm) {
    const int stride = q + 1;
    const int block = stride * stride;
    double ve = ve0;
    double va = std::max(0.1, 0.1 * ve);
    double vb = rank > 1 ? std::max(0.01, 0.01 * ve) : 0.0;
    std::copy(beta0, beta0 + q, beta);
    std::vector<double> h_all(k * rank * rank), correction_all(k * q * rank);
    std::vector<double> info_matrix(q * q), solution(q), normal(q * q), residual_cov(rank);

    // H=(Z'Z + ve D^-1)^-1, correction=F'Z H, normal=ve M.
    auto posterior_and_information = [&](bool rhs) {
        std::fill(normal.begin(), normal.end(), 0.0);
        if (rhs) std::fill(solution.begin(), solution.end(), 0.0);
        for (int g = 0; g < k; ++g) {
            const double* a = gram + g * block;
            double* h = h_all.data() + g * rank * rank;
            double* correction = correction_all.data() + g * q * rank;
            if (rank == 1) {
                h[0] = 1.0 / (a[0] + ve / va);
                for (int i = 0; i < q; ++i) correction[i] = a[i * stride] * h[0];
                for (int i = 0; i < q; ++i) {
                    for (int j = 0; j < q; ++j) {
                        normal[i * q + j] += a[i * stride + j] - correction[i] * a[j * stride];
                    }
                    if (rhs) solution[i] += a[i * stride + q] - correction[i] * a[q];
                }
            } else {
                const double intercept_precision = ve / va;
                const double slope_precision = ve / vb;
                for (int i = 0; i < rank; ++i) {
                    std::copy_n(a + i * stride, rank, h + i * rank);
                    h[i * rank + i] += i == 0 ? intercept_precision : slope_precision;
                }
                if (!inverse_spd(h, rank, factor, invert)) return false;
                for (int i = 0; i < q; ++i) {
                    const double di = i == 0 ? intercept_precision : slope_precision;
                    for (int j = 0; j < q; ++j) {
                        const double dj = j == 0 ? intercept_precision : slope_precision;
                        // Z=F. A H = I - (ve D^-1) H and A - A H A = A H (ve D^-1).
                        const double coefficient = (i == j ? 1.0 : 0.0) - di * h[i * q + j];
                        correction[i * q + j] = coefficient;
                        normal[i * q + j] += coefficient * dj;
                    }
                    if (rhs) {
                        solution[i] += di * dot(h + i * q, a + q * stride, q);
                    }
                }
            }
        }
        for (int i = 0; i < q; ++i) {
            for (int j = 0; j < i; ++j) {
                const double value = 0.5 * (normal[i * q + j] + normal[j * q + i]);
                normal[i * q + j] = normal[j * q + i] = value;
            }
        }
        return true;
    };

    if (!posterior_and_information(false)) return 1;
    int used = 0;
    bool converged = false;
    for (; used < max_iterations; ) {
        double intercept_sum = 0.0, slope_sum = 0.0, residual_sum = 0.0;
        for (int g = 0; g < k; ++g) {
            const double* a = gram + g * block;
            const double* h = h_all.data() + g * rank * rank;
            double* u = effects + g * rank;
            double rss = a[block - 1];
            for (int i = 0; i < q; ++i) {
                const double xb = dot(a + i * stride, beta, q);
                rss += beta[i] * (xb - 2.0 * a[i * stride + q]);
                if (i < rank) residual_cov[i] = a[i * stride + q] - xb;
            }
            for (int i = 0; i < rank; ++i) u[i] = dot(h + i * rank, residual_cov.data(), rank);
            for (int i = 0; i < rank; ++i) {
                const double second = u[i] * u[i] + ve * h[i * rank + i];
                if (i == 0) intercept_sum += second;
                else slope_sum += second;
                rss += u[i] * (dot(a + i * stride, u, rank) - 2.0 * residual_cov[i]);
                rss += ve * dot(a + i * stride, h + i * rank, rank);
            }
            residual_sum += rss;
        }
        const double new_va = std::max(epsilon, intercept_sum / k);
        const double new_vb = rank > 1 ? std::max(epsilon, slope_sum / (k * (rank - 1))) : 0.0;
        const double new_ve = std::max(epsilon, residual_sum / n);
        const double variance_delta = std::max({std::abs(new_va - va) / (1 + va),
            std::abs(new_vb - vb) / (1 + vb), std::abs(new_ve - ve) / (1 + ve)});
        va = new_va;
        vb = new_vb;
        ve = new_ve;
        if (!posterior_and_information(true)) return 2;
        std::copy(normal.begin(), normal.end(), info_matrix.begin());
        char lower = 'L';
        int status = 0, one = 1;
        factor(&lower, &q, info_matrix.data(), &q, &status);
        if (status != 0) return 3;
        solve(&lower, &q, &one, info_matrix.data(), &q, solution.data(), &q, &status);
        if (status != 0) return 4;
        double beta_squared_change = 0.0;
        for (int i = 0; i < q; ++i) {
            const double difference = solution[i] - beta[i];
            beta_squared_change += difference * difference;
        }
        const double beta_delta = std::sqrt(beta_squared_change) / (1 + std::sqrt(dot(beta, beta, q)));
        std::copy(solution.begin(), solution.end(), beta);
        ++used;
        if (std::max(variance_delta, beta_delta) < tolerance) {
            converged = true;
            break;
        }
    }

    if (!inverse_spd(normal.data(), q, factor, invert)) return 5;
    std::vector<double> map(q * q), product(q * q), full(q * q);
    double score = ve;
    for (int g = 0; g < k; ++g) {
        const double* a = gram + g * block;
        const double* w = weights + g * q * q;
        const double* h = h_all.data() + g * rank * rank;
        const double* correction = correction_all.data() + g * q * rank;
        double* u = effects + g * rank;
        for (int i = 0; i < rank; ++i) residual_cov[i] = a[i * stride + q] - dot(a + i * stride, beta, q);
        for (int i = 0; i < rank; ++i) {
            u[i] = dot(h + i * rank, residual_cov.data(), rank);
            score += ve * dot(h + i * rank, w + i * q, rank);
        }
        std::fill(map.begin(), map.end(), 0.0);
        for (int i = 0; i < q; ++i) {
            for (int j = 0; j < rank; ++j) map[i * q + j] = -correction[i * rank + j];
            map[i * q + i] += 1.0;
        }
        char t = 'T', no = 'N';
        double one = 1.0, zero = 0.0;
        // In column-major BLAS, the row-major map buffer represents T'.
        // Thus these products form (T W) T' with W symmetric.
        gemm(&t, &no, &q, &q, &q, &one, map.data(), &q, const_cast<double*>(w), &q,
             &zero, product.data(), &q);
        gemm(&no, &no, &q, &q, &q, &one, product.data(), &q, map.data(), &q,
             &zero, full.data(), &q);
        score += ve * dot(normal.data(), full.data(), q * q);
    }
    if (!std::isfinite(score) || score <= 0 || !std::all_of(beta, beta + q, [](double x) { return std::isfinite(x); }) ||
        !std::all_of(effects, effects + k * rank, [](double x) { return std::isfinite(x); })) return 6;
    output[0] = score;
    output[1] = va;
    output[2] = vb;
    output[3] = ve;
    output[4] = used;
    output[5] = converged;
    return 0;
}
