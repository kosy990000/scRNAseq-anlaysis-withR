

#---------- 유틸 함수 ----------  #

compute_qc_limits <- function(df, use_quantile = FALSE, q = c(0.01, 0.99)) {
  rng <- function(x) {
    if (use_quantile) range(quantile(x, q, na.rm = TRUE)) else range(x, na.rm = TRUE)
  }
  list(
    total_counts      = rng(df$total_counts),        # for x/y log scales
    n_genes_by_counts = rng(df$n_genes_by_counts),   # for y log scale
    pct_counts_mt     = rng(df$pct_counts_mt)        # for color & y (linear)
  )
}

# 입력 데이터를 받아 qc_before와 qc_After를 보는 함수
qc_plots <- function(df) {
  # 로그 스케일 안전: 0/NA 제거용
  df_log <- subset(
    df,
    is.finite(total_counts) & is.finite(n_genes_by_counts) &
      total_counts > 0 & n_genes_by_counts > 0
  )
  
  # 1) total_counts 바이올린 (log Y)
  p1 <- ggplot(df_log, aes(x = condition, y = total_counts, fill = condition)) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.2, size = 0.3, alpha = 0.4) +
    scale_y_log10(labels = scales::comma) +
    labs(x = NULL, y = "n_total_counts") +
    theme_bw() +
    theme(legend.position = "none")
  
  # 2) pct_counts_mt 바이올린 (linear Y)
  p2 <- ggplot(df, aes(x = condition, y = pct_counts_mt, fill = condition)) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.2, size = 0.3, alpha = 0.4) +
    labs(x = NULL, y = "pct_counts_mt") +
    theme_bw() +
    theme(legend.position = "none")
  
  # 3) total_counts 히스토그램 (facet, 축 공유)
  p3 <- ggplot(df, aes(x = total_counts)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "black") +
    scale_x_log10(labels = scales::comma) +
    labs(x = "total_counts", y = "Frequency") +
    theme_bw() +
    facet_wrap(~condition, nrow = 1, scales = "fixed")
  
  # 4) total_counts vs n_genes_by_counts (facet, log-log, 축 공유)
  df_log2 <- subset(
    df,
    is.finite(total_counts) & is.finite(n_genes_by_counts) &
      total_counts > 0 & n_genes_by_counts > 0
  )
  p4 <- ggplot(df_log2, aes(x = total_counts, y = n_genes_by_counts)) +
    geom_point(alpha = 0.5, size = 0.6, color = "steelblue") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(x = "total_counts", y = "n_genes_by_counts") +
    theme_bw() +
    facet_wrap(~condition, nrow = 1, scales = "fixed")
  
  # 5) total_counts vs mt_counts (facet, log-log, 축 공유)
  df$mt_counts <- df$pct_counts_mt * df$total_counts
  df_log_mt <- subset(
    df,
    is.finite(total_counts) & is.finite(mt_counts) &
      total_counts > 0 & mt_counts > 0
  )
  p5 <- ggplot(df_log_mt, aes(x = total_counts, y = mt_counts)) +
    geom_point(alpha = 0.5, size = 0.6, color = "steelblue") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(x = "total_counts", y = "mt_counts") +
    theme_bw() +
    facet_wrap(~condition, nrow = 1, scales = "fixed") +   # ✅ 여기 괄호 닫기
    geom_abline(slope = 1, intercept = log10(0.1), 
                color = "red", linetype = "dashed") +
    geom_abline(slope = 1, intercept = log10(0.2), 
                color = "blue", linetype = "dashed") +
    geom_abline(slope = 1, intercept = log10(0.3), 
                color = "green", linetype = "dashed")
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5)
}


is_outlier <- function(x, nmads = 5, verbose = FALSE, na_to_false = TRUE) {
  # 길이 0 또는 전부 NA 방지
  if (length(x) == 0 || all(is.na(x))) {
    return(rep(FALSE, length(x)))
  }
  med <- median(x, na.rm = TRUE)
  dev <- mad(x, center = med, constant = 1, na.rm = TRUE)
  
  # mad가 0이거나 NA인 경우(변동 거의 없음) → 모두 FALSE
  if (!is.finite(med) || !is.finite(dev) || dev == 0) {
    res <- rep(FALSE, length(x))
    if (verbose) cat("Lower/Upper 미정 (dev=0 또는 NA). 모두 FALSE 처리.\n")
    return(res)
  }
  
  lower <- med - nmads * dev
  upper <- med + nmads * dev
  res <- (x < lower) | (x > upper)
  
  # NA 처리
  if (na_to_false) res[is.na(res)] <- FALSE
  
  if (verbose) cat("Lower:", lower, "Upper:", upper, "\n")
  return(res)
}

pct_counts_in_top_n_genes <- function(counts_mat, n = 20) {
  total_counts <- Matrix::colSums(counts_mat)
  top_counts <- vapply(
    seq_len(ncol(counts_mat)),
    function(j) {
      col <- counts_mat[, j]
      if (sum(col) == 0) return(0)
      sum(sort(col, decreasing = TRUE)[seq_len(min(n, length(col)))])
    },
    numeric(1)
  )
  ifelse(total_counts > 0, top_counts / total_counts, NA_real_)  # ← 0–1 스케일 유지
}
