# ------------------------------------------------------------------------------
# 0. 準備 (Setup)
# ------------------------------------------------------------------------------
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
# 1. データ構築 (Data Construction)
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
    
    # Outcome 1: LFP (Main)
    In_Labor_Force = if_else(LABFORCE == 2, 1, 0),
    
    # Outcome 2: Enrollment (Secondary)
    In_School = if_else(SCHLCOLL %in% c(1, 2, 3, 4), 1, 0),
    
    # Outcome 3: Employment (Robustness)
    # EMPSTAT: 10=At work, 12=Has job, not at work -> Employed
    # Note: Check if EMPSTAT exists. If not, assume LABFORCE=2 & EMPSTAT!=20s
    Is_Employed = if_else(EMPSTAT %in% c(10, 12), 1, 0),
    
    # Time Trend variable for Robustness
    Time_Index = as.numeric(Date),
    
    # Event Study Time
    rel_month = (year(Date) - 2022) * 12 + (month(Date) - 6)
  )

# サブセット作成
df_women <- df_all %>% filter(Female == 1)
df_men   <- df_all %>% filter(Female == 0)

# ------------------------------------------------------------------------------
# 2. メイン結果 (Main Results: Table 1 Replication)
# ------------------------------------------------------------------------------
print("▼▼▼ Main Results (Table 1) ▼▼▼")

# (1) Women DiD (LFP)
mod_main_women <- feols(In_Labor_Force ~ Treat * Post | STATEFIP + Date,
                        data = df_women, weights = ~ WTFINL, cluster = ~ STATEFIP)

# (2) Men Placebo (LFP)
mod_main_men <- feols(In_Labor_Force ~ Treat * Post | STATEFIP + Date,
                      data = df_men, weights = ~ WTFINL, cluster = ~ STATEFIP)

# (3) DDD (LFP)
mod_main_ddd <- feols(In_Labor_Force ~ Treat * Post * Female | STATEFIP^Female + Date^Female,
                      data = df_all, weights = ~ WTFINL, cluster = ~ STATEFIP)

# 結果表示用関数
print_res <- function(label, model, param) {
  res <- coeftable(model)
  if(param %in% rownames(res)){
    est <- res[param, "Estimate"]
    se  <- res[param, "Std. Error"]
    pv  <- res[param, "Pr(>|t|)"]
    stars <- ifelse(pv < 0.01, "***", ifelse(pv < 0.05, "**", ifelse(pv < 0.1, "*", "")))
    cat(paste0(str_pad(label, 25, "right"), ": Est=", round(est, 4), " / SE=", round(se, 4), " / P=", round(pv, 3), " ", stars, "\n"))
  } else {
    cat(paste0(label, ": Parameter not found\n"))
  }
}

print_res("Women LFP (Main)", mod_main_women, "Treat:Post")
print_res("Men LFP (Placebo)", mod_main_men, "Treat:Post")
print_res("DDD (Main)", mod_main_ddd, "Treat:Post:Female")

# ------------------------------------------------------------------------------
# 3. [修正版] DDD Event Study (Pre-trend check)
# ------------------------------------------------------------------------------
print("\n▼▼▼ Robustness 1: DDD Event Study (修正版・再開) ▼▼▼")

# 分析用データセット（イベント期間のみ）
df_es_ddd <- df_all %>% 
  filter(rel_month >= -18, rel_month <= 18) %>%
  # 【重要】ここで交差項変数を明示的に作ってしまう
  mutate(Treat_Female = Treat * Female)

# DDD Event Study Specification
# Treat_Female を使うことでエラーを回避
mod_es_ddd <- feols(In_Labor_Force ~ i(rel_month, Treat_Female, ref = -1) + 
                      i(rel_month, Treat, ref = -1) + 
                      i(rel_month, Female, ref = -1) 
                    | STATEFIP^Female + Date^Female,
                    data = df_es_ddd, weights = ~ WTFINL, cluster = ~ STATEFIP)

# 係数抽出 (変数名が変わったので正規表現も微調整)
es_ddd_coef <- tidy(mod_es_ddd) %>%
  filter(str_detect(term, "rel_month::.*:Treat_Female")) %>% 
  mutate(
    rel_month = as.integer(str_extract(term, "-?\\d+")),
    conf_low  = estimate - 1.96 * std.error,
    conf_high = estimate + 1.96 * std.error
  )

# Plot DDD Pre-trends
gg_es_ddd <- ggplot(es_ddd_coef, aes(x = rel_month, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "red") +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.2, color = "black") +
  geom_point(size = 2, color = "black") +
  labs(x = "Months relative to Dobbs", y = "Triple-Diff Estimate (pp)", 
       title = "DDD Event Study: Impact on Women relative to Men",
       subtitle = "Check: Should be flat around 0 before Dobbs") +
  theme_minimal()

ggsave("ddd_event_study.png", gg_es_ddd, width = 7, height = 4.5)
print(" -> ddd_event_study.png を保存しました。")

# ------------------------------------------------------------------------------
# 4. [必須] サブサンプル分析 (Composition check) - ここから残り全部
# ------------------------------------------------------------------------------
print("\n▼▼▼ Robustness 2: Subsample Analysis ▼▼▼")

# (A) Age 20-24 Only (高校生除外)
print("--- (A) Age 20-24 Only ---")
df_2024 <- df_women %>% filter(AGE >= 20)
mod_2024 <- feols(In_Labor_Force ~ Treat * Post | STATEFIP + Date,
                  data = df_2024, weights = ~ WTFINL, cluster = ~ STATEFIP)
print_res("Women 20-24 Only", mod_2024, "Treat:Post")

# (B) Non-Students Only (学生バイト除外)
print("--- (B) Non-Students Only ---")
df_nonstudent <- df_women %>% filter(In_School == 0)
mod_nonstudent <- feols(In_Labor_Force ~ Treat * Post | STATEFIP + Date,
                        data = df_nonstudent, weights = ~ WTFINL, cluster = ~ STATEFIP)
print_res("Non-Students Only", mod_nonstudent, "Treat:Post")

# ------------------------------------------------------------------------------
# 5. [必須] アウトカム分解 (Employment vs Unemployment)
# ------------------------------------------------------------------------------
print("\n▼▼▼ Robustness 3: Outcome Decomposition (Employment) ▼▼▼")

# (C) Employment Rate (EPOP)
mod_emp <- feols(Is_Employed ~ Treat * Post | STATEFIP + Date,
                 data = df_women, weights = ~ WTFINL, cluster = ~ STATEFIP)
print_res("Employment Rate", mod_emp, "Treat:Post")

# ------------------------------------------------------------------------------
# 6. [推奨] ショートウィンドウ & トレンド
# ------------------------------------------------------------------------------
print("\n▼▼▼ Robustness 4: Short Window & Trends ▼▼▼")

# (D) Short Window (2021-2022 End)
df_short <- df_women %>% filter(Date <= as.Date("2022-12-31"))
mod_short <- feols(In_Labor_Force ~ Treat * Post | STATEFIP + Date,
                   data = df_short, weights = ~ WTFINL, cluster = ~ STATEFIP)
print_res("Short Window (to 2022)", mod_short, "Treat:Post")

# (E) State-specific Linear Trends
mod_trend <- feols(In_Labor_Force ~ Treat * Post | STATEFIP + Date + STATEFIP[Time_Index],
                   data = df_women, weights = ~ WTFINL, cluster = ~ STATEFIP)
print_res("With State Trends", mod_trend, "Treat:Post")

print("\n完了しました。これらの数値を確認してください。")