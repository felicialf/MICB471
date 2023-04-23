# load packages
library(tidyverse)
library(rstatix)

# Growth Curve ------------------------------------------------------------
df <- read.csv("~/undergrad/amos files/20230217_MICB471.csv") 

p <- df %>% 
  mutate(time = seq(from = 0, to = nrow(df)*0.5 - 0.5, by = 0.5)) %>%
  select(-Time, -Temp)
p <- pivot_longer(p, p %>% select(-time) %>% colnames()) 

p <- p %>%
  mutate(row = str_sub(.$name, start=1L, end=1L)) %>%
  mutate(column = as.numeric(str_sub(.$name, start=2L, end=-1L))) %>%
  filter(column <= 7) %>%
  filter(!row %in% c("F", "G", "H"))

replicate <- c()
for(i in seq_along(p$row)){
  if (p[i, "row"] %in% c("A")){
    replicate = c(replicate, 1)}
  else if (p[i, "row"] %in% c("B")){
    replicate = c(replicate, 2)}
  else if (p[i, "row"] %in% c("C")){
    replicate = c(replicate, 3)}
  else if (p[i, "row"] %in% c("D")){
    replicate = c(replicate, 4)}
  else if (p[i, "row"] %in% c("E")){
    replicate = c(replicate, 5)}
}

treatment <- c()
for(i in seq_along(p$column)){
  if (p[i,"column"] == 1){
    treatment = c(treatment, "0% EtOH")}
  else if(p[i,"column"] == 2){
    treatment = c(treatment, "1% EtOH")}
  else if(p[i,"column"] == 3){
    treatment = c(treatment, "2% EtOH")}
  else if(p[i,"column"] == 4){
    treatment = c(treatment, "3% EtOH")}
  else if(p[i,"column"] == 5){
    treatment = c(treatment, "4% EtOH")}
  else if(p[i,"column"] == 6){
    treatment = c(treatment, "5% EtOH")}
  else if(p[i,"column"] == 7){
    treatment = c(treatment, "0% EtOH + Amp")}
  else{treatment = c(treatment, NA)}
}

concentration <- c()
for(i in seq_along(p$column)){
  if (p[i,"column"] == 1){
    concentration = c(concentration, 0)}
  else if(p[i,"column"] == 2){
    concentration = c(concentration, 1)}
  else if(p[i,"column"] == 3){
    concentration = c(concentration, 2)}
  else if(p[i,"column"] == 4){
    concentration = c(concentration, 3)}
  else if(p[i,"column"] == 5){
    concentration = c(concentration, 4)}
  else if(p[i,"column"] == 6){
    concentration = c(concentration, 5)}
  else if(p[i,"column"] == 7){
    concentration = c(concentration, 0)}
  else{concentration = c(concentration, NA)}
}

antibiotic <- c()
for(i in seq_along(p$column)){
  if (p[i,"column"] <= 6){
    antibiotic = c(antibiotic, F)}
  else if(p[i,"column"] == 7){
    antibiotic = c(antibiotic, T)}
  else{concentration = c(concentration, NA)}
}

p <- p %>%
  mutate(treatment = treatment) %>%
  mutate(replicate = replicate) %>%
  mutate(concentration = as.character(concentration)) %>%
  mutate(antibiotic = antibiotic) %>%
  filter(!is.na(treatment)) %>%
  group_by(treatment, time)
p$treatment <- factor(p$treatment, levels = c("0% EtOH","1% EtOH","2% EtOH","3% EtOH","4% EtOH","5% EtOH","0% EtOH + Amp"))

summary_stats <- p %>%
  group_by(treatment, time) %>%
  get_summary_stats(value)

p <- p %>% left_join(summary_stats, by = c("time", "treatment")) %>% ungroup()

my_times <- c(6, 12, 18, 24)
all <- p %>% filter(time != 0 & !str_detect(treatment, "Amp"))
fig <- p %>% filter(time %in% my_times & !str_detect(treatment, "Amp"))

# anova test
anova <- fig %>%
  group_by(time) %>%
  anova_test(value~treatment, dv = treatment)

# multiple comparisons test
multiple_comparisons <- fig %>%
  group_by(time) %>%
  tukey_hsd(value~treatment, p.adjust.method = "holm", detailed = TRUE) %>%
  add_xy_position()

# anova test
anova_all <- all %>%
  group_by(time) %>%
  anova_test(value~treatment, dv = treatment)

# multiple comparisons test
multiple_comparisons_all <- all %>%
  group_by(time) %>%
  tukey_hsd(value~treatment, p.adjust.method = "holm", detailed = TRUE) %>%
  add_xy_position()

# plot data
ggplot(p, aes(x=time, y=mean, group=treatment, color=concentration)) +
  geom_line(aes(linetype = antibiotic), alpha = 0.5) +
  geom_point(size = 1, alpha=0.5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, width=0.25), alpha = 0.5) +
  labs(y="Average OD600", x="Time [hours]", color = "% EtOH", linetype = "Ampicillin")

ggsave("20230217_MICB471_ethanol_growth_curve.png", height = 5, width = 8)

ggplot(fig %>% filter(time %in% my_times), aes(x=concentration, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = concentration),position = position_jitter(width=0.25), size = 2.5, alpha = 0.5) +
  ggpubr::stat_pvalue_manual(multiple_comparisons %>% filter(group1 == "0% EtOH"), label = "p = {scales::pvalue(p.adj)}") + 
  labs(x = "% EtOH", y = "OD600", color = "% EtOH") +
  facet_wrap(~time, labeller = label_both, ncol = 2)

ggsave("20230217_MICB471_ethanol_growth_curve_endpoints.png", height = 7.5, width = 10)

ggplot(all, aes(x=concentration, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = concentration),position = position_jitter(width=0.25), size = 2.5, alpha = 0.5) +
  ggpubr::stat_pvalue_manual(multiple_comparisons_all %>% filter(group1 == "0% EtOH"), label = "p = {scales::pvalue(p.adj)}") + 
  labs(x = "% EtOH", y = "OD600", color = "% EtOH") +
  facet_wrap(~time, labeller = label_both, ncol = 8)

ggsave("20230217_MICB471_ethanol_growth_curve_endpoints_all.png", height = 20, width = 30)

# calculate derivatives ---------------------------------------------------
deriv0 <- p %>%
  filter(treatment == "0% EtOH") %>%
  select(antibiotic, concentration, mean, time, treatment) %>%
  unique() %>%
  mutate(slope = (mean - lag(mean))/(time - lag(time)))

deriv1 <- p %>%
  filter(treatment == "1% EtOH") %>%
  select(antibiotic, concentration, mean, time, treatment) %>%
  unique() %>%
  mutate(slope = (mean - lag(mean))/(time - lag(time)))

deriv2 <- p %>%
  filter(treatment == "2% EtOH") %>%
  select(antibiotic, concentration, mean, time, treatment) %>%
  unique() %>%
  mutate(slope = (mean - lag(mean))/(time - lag(time)))

deriv3 <- p %>%
  filter(treatment == "3% EtOH") %>%
  select(antibiotic, concentration, mean, time, treatment) %>%
  unique() %>%
  mutate(slope = (mean - lag(mean))/(time - lag(time)))

deriv4 <- p %>%
  filter(treatment == "4% EtOH") %>%
  select(antibiotic, concentration, mean, time, treatment) %>%
  unique() %>%
  mutate(slope = (mean - lag(mean))/(time - lag(time)))

deriv5 <- p %>%
  filter(treatment == "5% EtOH") %>%
  select(antibiotic, concentration, mean, time, treatment) %>%
  unique() %>%
  mutate(slope = (mean - lag(mean))/(time - lag(time)))

deriv0amp <- p %>%
  filter(treatment == "0% EtOH + Amp") %>%
  select(antibiotic, concentration, mean, time, treatment) %>%
  unique() %>%
  mutate(slope = (mean - lag(mean))/(time - lag(time)))

# compile derivatives -----------------------------------------------------
deriv <- bind_rows(deriv0, deriv1, deriv2, deriv3, deriv4, deriv5, deriv0amp)

# plot data ---------------------------------------------------------------
ggplot(deriv, aes(x = time, y = slope, color = concentration)) + 
  geom_line(aes(linetype = antibiotic), alpha = 0.5) +
  geom_point(size = 1, alpha = 0.5) +
  labs(y="Growth Rate [OD600/hours]", x="Time [hours]", color = "% EtOH", linetype = "Ampicillin")

ggsave("20230217_MICB471_ethanol_growth_curve_derivative.png", height = 5, width = 8)



