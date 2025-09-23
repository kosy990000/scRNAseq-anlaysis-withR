


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(stringr)
  library(scran)
  library(patchwork)
  library(ggpubr)
})


path <- getwd()
folder_path = paste0(path,"/data/GSE")

# 또는 (권장: 운영체제 호환 안전)
df <- read.csv(file.path(folder_path, "ALL_files_QC_matrix_with_doublet.csv"))

# 컬럼 추가
df <- df %>%
  mutate(
    group = case_when(
      str_detect(file, "Q") ~ "Quercetin",
      str_detect(file, "V") ~ "Vector",
      TRUE ~ "Other"
    )
  )

#데이터 분리
# 리스트로 나누기
df_split <- split(df, df$group)

# 리스트 이름 확인
names(df_split)

# set group number
group_A = "Quercetin"
group_B = "Vector"

df$mt_counts <- df$pct_counts_mt * df$total_counts

# 2) 무한/NA 제외한 한계값 계산 (변수명: limits_ 로 충돌 회피)
limits_ <- list(
  total_counts = range(df$total_counts[is.finite(df$total_counts)], na.rm = TRUE),
  mt_counts    = range(df$mt_counts[is.finite(df$mt_counts)],    na.rm = TRUE)
)




# 3) 직접 df로 “직빵” 플롯 
conv_raw <- 
  ggplot(df, aes(x = total_counts, y = mt_counts, color = group)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "total_counts", y = "mt_counts", color = "group") +
  theme_bw() + theme(legend.position = "top") +
  coord_cartesian(
    xlim = c(0, 75000),   # total_counts 범위
    ylim = c(0, 5000)    # mt_counts 범위
  )+
  theme(strip.placement = "outside", strip.background = element_blank()) +
  stat_function(fun = function(x) 0.1 * x, color = "black", linetype = "dashed") +
  stat_function(fun = function(x) 0.2 * x, color = "blue",  linetype = "dashed") +
  stat_function(fun = function(x) 0.3 * x, color = "red",   linetype = "dashed") +
  scale_color_manual(
    values = c(
      "Quercetin" = "red",       # Treatment A
      "Vector" = "black" # Control
    )
  )

conv_raw

# -------------------------------------------------------------- #
# 3) 직접 df로 “직빵” 플롯 
conv <- df %>%
  filter(outlier_with_doublet == "FALSE") %>%
  ggplot(aes(x = total_counts, y = mt_counts, color = group)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "total_counts", y = "mt_counts", color = "group") +
  theme_bw() + theme(legend.position = "top") +
  coord_cartesian(
    xlim = c(0, 75000),   # total_counts 범위
    ylim = c(0, 5000)    # mt_counts 범위
  )+
  theme(strip.placement = "outside", strip.background = element_blank()) +
  stat_function(fun = function(x) 0.1 * x, color = "black", linetype = "dashed") +
  stat_function(fun = function(x) 0.2 * x, color = "blue",  linetype = "dashed") +
  stat_function(fun = function(x) 0.3 * x, color = "red",   linetype = "dashed") +
  scale_color_manual(
    values = c(
      "Quercetin" = "red",       # Treatment A
      "Vector" = "black" # Control
    )
  )

conv

# -----------------------------------------------------------------------------#
pval_df_1 <- df %>%
  summarise(
    test = list(wilcox.test(pct_counts_mt ~ group)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p = test$p.value,
    signif = case_when(
      p > 0.05 ~ "ns",
      p <= 0.05 & p > 0.01 ~ "*",
      p <= 0.01 & p > 0.001 ~ "**",
      p <= 0.001 & p > 0.0001 ~ "***",
      p <= 0.0001 ~ "****"
    ),
    label = paste0("p = ", signif)
  ) %>%
  select(p, signif, label)


conv2 <- ggplot(df, aes(x = group, y = pct_counts_mt)) +
  geom_violin(trim = TRUE, width = 1, fill = "red", alpha = 0.55) +
  # 박스플롯 오버레이(윤곽선만)
  geom_boxplot(width = 0.18, outlier.shape = NA,
               fill = "skyblue", color = "black", alpha = 0.8) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "pct_counts_mt") +
  theme_bw() +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",   # ns, *, **, ***
                     comparisons = list(c("Quercetin","Vector")))  # 그룹 이름 매칭

conv2
#------------------------------------------------------------------------------#


conv3 <- df %>%
  filter(outlier_with_doublet == "FALSE") %>%
  ggplot(aes(x = group, y = pct_counts_mt)) +
  geom_violin(trim = TRUE, width = 1, fill = "red", alpha = 0.55) +
  # 박스플롯 오버레이(윤곽선만)
  geom_boxplot(width = 0.18, outlier.shape = NA,
               fill = "skyblue", color = "black", alpha = 0.8) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "pct_counts_mt") +
  theme_bw() +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",   # ns, *, **, ***
                     comparisons = list(c("Quercetin","Vector")))  # 그룹 이름 매칭


conv3
