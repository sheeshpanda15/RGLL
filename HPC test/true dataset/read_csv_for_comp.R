# ── 读取数据 ───────────────────────────────────────────────
dat <- read.csv("2004.csv", header = TRUE)

# 负数全部替换为 NA
dat[dat < 0] <- NA

# 删除缺失率过高的列（家庭净资产，缺失23.8%）
dat$R8378703 <- NULL

# 只保留 y 不缺失的行
dat <- dat[!is.na(dat$R8495200), ]

# 删除有任何 NA 的行
dat <- dat[complete.cases(dat), ]
cat("完整无缺失样本量:", nrow(dat), "\n")

# ── 提取 y 和 X ────────────────────────────────────────────
y <- as.numeric(dat$R8495200)

exclude <- c("R8495200",  # y
             "R0000100",  # CASEID
             "R0173600")  # SAMPLE_ID
X_cols <- setdiff(colnames(dat), exclude)
X <- as.matrix(dat[, X_cols])
cat("X 列数 p:", ncol(X), "\n")
cat("X 列名:", paste(colnames(X), collapse=", "), "\n\n")

# 全局 min-max 标准化 X 到 [-1, 1]
for (j in 1:ncol(X)) {
  lo <- min(X[, j]); hi <- max(X[, j])
  X[, j] <- if (hi > lo) 2 * (X[, j] - lo) / (hi - lo) - 1 else 0
}

# ── 分组 + 70/30 分割函数 ──────────────────────────────────
make_group_data <- function(group_col, min_size = 2, train_ratio = 0.7, seed = 42) {
  set.seed(seed)
  g <- as.integer(as.factor(dat[[group_col]]))

  valid_grps <- as.integer(names(table(g)[table(g) >= min_size]))
  keep <- g %in% valid_grps
  if (sum(!keep) > 0) cat("[", group_col, "] 剔除小组样本:", sum(!keep), "条\n")

  g   <- as.integer(as.factor(g[keep]))
  X_k <- X[keep, ]
  y_k <- y[keep]

  ord    <- order(g)
  X_s    <- X_k[ord, ]; y_s <- y_k[ord]; g_s <- g[ord]
  groups <- unique(g_s)

  train_idx <- c(); test_idx <- c()
  for (grp in groups) {
    idx  <- which(g_s == grp)
    n_tr <- min(max(1, round(length(idx) * train_ratio)), length(idx) - 1)
    sel  <- sample(idx, n_tr)
    train_idx <- c(train_idx, sort(sel))
    test_idx  <- c(test_idx,  sort(setdiff(idx, sel)))
  }

  list(
    FXX.train = X_s[train_idx, ],
    FY.train  = as.matrix(y_s[train_idx]),
    nc.train  = as.integer(table(factor(g_s[train_idx], levels = groups))),
    FXX.test  = X_s[test_idx, ],
    FY.test   = as.matrix(y_s[test_idx]),
    nc.test   = as.integer(table(factor(g_s[test_idx],  levels = groups))),
    R = length(groups),
    p = ncol(X)
  )
}

# ── 六种分组方案 ───────────────────────────────────────────
d_race   <- make_group_data("R0214700")               # 种族,   3组
d_sex    <- make_group_data("R0214800")               # 性别,   2组
d_region <- make_group_data("R8496500")               # 地区,   4组
d_worker <- make_group_data("R7898500")               # 雇主,   5组
d_urban  <- make_group_data("R8498600")               # 城乡,   3组
d_edu    <- make_group_data("R8497000", min_size = 30) # 教育,  ~15组

cat("=== 分组信息 ===\n")
for (name in c("d_race","d_sex","d_region","d_worker","d_urban","d_edu")) {
  d <- get(name)
  cat(sprintf("%-10s R=%-3d p=%-3d Train=%-5d Test=%-5d 组大小[%d-%d]\n",
              name, d$R, d$p, sum(d$nc.train), sum(d$nc.test),
              min(d$nc.train), max(d$nc.train)))
}
