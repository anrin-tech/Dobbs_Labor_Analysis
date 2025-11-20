# ------------------------------------------------------------------------------
# 0. 準備 (Setup)
# ------------------------------------------------------------------------------

# 必要なパッケージのみインストール・読み込み
if (!require("ipumsr")) install.packages("ipumsr")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("fixest")) install.packages("fixest")
if (!require("broom")) install.packages("broom")

library(ipumsr)
library(tidyverse)
library(fixest)
library(broom)
library(lubridate)

# ファイルパス (環境に合わせて変更)
ddi_file  <- "/Users/LLENN/Desktop/Research/Data/cps_00001.xml"
data_file <- "/Users/LLENN/Desktop/Research/Data/cps_00001.dat"

# ------------------------------------------------------------------------------
# 1. データ構築 (Data)
# ------------------------------------------------------------------------------
print("データを読み込んでいます...")
cps_raw <- read_ipums_micro(ddi_file, data_file = data_file)

ban_states <- c(1, 5, 16, 21, 22, 28, 29, 40, 46, 47, 48, 54, 55)
protected_states <- c(6, 36, 53, 41, 17, 25, 9, 24)

print("データを構築中...")

df_all <- cps_raw %>%
  filter(AGE >= 18 & AGE <= 24) %>%
  filter(STATEFIP %in% c(ban_states, protected_states)) %>%
  mutate(
    Date = as.Date(paste(YEAR, MONTH, "01", sep = "-")),
    Treat = if_else(STATEFIP %in% ban_states, 1, 0),
    Post = if_else(Date >= as.Date("2022-07-01"), 1, 0),
    Female = if_else(SEX == 2, 1, 0),
    In_Labor_Force = if_else(LABFORCE == 2, 1, 0),
    In_School = if_else(SCHLCOLL %in% c(1, 2, 3, 4), 1, 0),
    rel_month = (year(Date) - 2022) * 12 + (month(Date) - 6)
  )

df_women <- df_all %>% filter(Female == 1)
df_men   <- df_all %>% filter(Female == 0)

# ------------------------------------------------------------------------------
# 2. 分析 (Analysis - Standard Clustered SE)
# ------------------------------------------------------------------------------
print("▼▼▼ 分析実行 (Standard Cluster SE) ▼▼▼")

# (1) Women DiD (LFP)
mod_women <- feols(In_Labor_Force ~ Treat * Post | STATEFIP + Date,
                   data = df_women, weights = ~ WTFINL, cluster = ~ STATEFIP)

# (2) Men Placebo (LFP)
mod_men <- feols(In_Labor_Force ~ Treat * Post | STATEFIP + Date,
                 data = df_men, weights = ~ WTFINL, cluster = ~ STATEFIP)

# (3) DDD (LFP)
mod_ddd <- feols(In_Labor_Force ~ Treat * Post * Female | STATEFIP^Female + Date^Female,
                 data = df_all, weights = ~ WTFINL, cluster = ~ STATEFIP)

# (4) Women DiD (School)
mod_school <- feols(In_School ~ Treat * Post | STATEFIP + Date,
                    data = df_women, weights = ~ WTFINL, cluster = ~ STATEFIP)

# ------------------------------------------------------------------------------
# 3. 結果出力
# ------------------------------------------------------------------------------
print_res <- function(label, model, param) {
  res <- coeftable(model)
  est <- res[param, "Estimate"]
  se  <- res[param, "Std. Error"]
  pv  <- res[param, "Pr(>|t|)"]
  stars <- ifelse(pv < 0.01, "***", ifelse(pv < 0.05, "**", ifelse(pv < 0.1, "*", "")))
  
  cat(paste0("\n--- ", label, " ---\n"))
  cat(paste0("Est: ", round(est, 4), " / SE: ", round(se, 4), "\n"))
  cat(paste0("P-val: ", round(pv, 4), " ", stars, "\n"))
}

print_res("Women LFP", mod_women, "Treat:Post")
print_res("Men Placebo", mod_men, "Treat:Post")
print_res("DDD Check", mod_ddd, "Treat:Post:Female")
print_res("Schooling", mod_school, "Treat:Post")

# ------------------------------------------------------------------------------
# 4. 図の作成 (Event Study)
# ------------------------------------------------------------------------------
print("\n▼▼▼ 図を作成中... ▼▼▼")

df_es <- df_women %>% filter(rel_month >= -18, rel_month <= 18)

es_model <- feols(In_Labor_Force ~ i(rel_month, Treat, ref = -1) | STATEFIP + Date,
                  data = df_es, weights = ~ WTFINL, cluster = ~ STATEFIP)

es_coef <- tidy(es_model) %>%
  filter(str_detect(term, "rel_month::Treat")) %>%
  mutate(
    rel_month = as.integer(str_extract(term, "-?\\d+")),
    conf_low  = estimate - 1.96 * std.error,
    conf_high = estimate + 1.96 * std.error
  )

gg_es <- ggplot(es_coef, aes(x = rel_month, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "red") +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.2, fill = "#0072B2") +
  geom_point(size = 2, color = "#0072B2") +
  geom_line(color = "#0072B2") +
  labs(x = "Months relative to Dobbs", y = "LFP Impact (pp)", 
       title = "Event-Study: Young Women's LFP") +
  theme_minimal()

ggsave("lfp_event_study_qtr.png", gg_es, width = 7, height = 4.5)
print("完了。")