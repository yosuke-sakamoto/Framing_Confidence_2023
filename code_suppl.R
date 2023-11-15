# Import library
library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(cowplot)
library(ggeffects)
library(rstatix)
source("anovakun_489.txt")

# Import data (variables are standardized)
df <- read.csv("data_behavior.csv", header = TRUE)
df$LogDiff <- scale(df$LogDiff)
df$LogTot <- scale(df$LogTot)
df$LogChVal <- scale(df$LogChVal)
df$LogUnChVal <- scale(df$LogUnChVal)
df$Conf <- scale(df$Conf)
df$LogChoiceRT <- scale(log(df$ChoiceRT))

df <- mutate(df, DiffRL = RVal - LVal, LogDiffRL = LogRVal - LogLVal)
df <- mutate(df, Baseline = if_else(TotVal <= 110, 0.1, if_else(TotVal <= 220, 0.2, 0.3)))
df_more <- filter(df, BlockCond == "MORE")
df_less <- filter(df, BlockCond == "LESS")
################################################################################
# Plot mean accuracy across baseline (Figure S1)
df %>% 
  group_by(id, Baseline, BlockCond) %>% 
  summarise(Corr = sum(Corr) / n()) %>% 
  group_by(Baseline, BlockCond) %>% 
  summarise(M = mean(Corr), SD = sd(Corr)) %>% 
  ggplot(aes(x = Baseline, y= M, color = BlockCond)) + geom_point(size = 4) + geom_line() + 
  geom_errorbar(aes(ymin = M - SD, ymax = M + SD), width = 0.01, position = position_dodge()) + 
  scale_x_continuous(breaks = seq(0.1, 0.3, 0.1), limits = c(0.095, 0.31), labels = c(50, 100, 150)) + 
  scale_y_continuous(breaks =seq(0.6, 0.9, 0.1), limits = c(0.6, 0.85)) + 
  xlab("Baseline level") + ylab("Mean accuracy") + 
  theme_classic(base_size = 16, base_line_size = 1) + scale_color_grey() +
  scale_color_manual(values = c("#080808", "#e0e0e0"), labels = c("More", "Less")) + 
  guides(color = guide_legend(title = "Frame", order = 1)) + 
  theme(legend.position = c(0.65, 0.8), 
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12), 
        axis.text = element_text(face = "bold")) -> ps1
ggsave("FigureS1.png", width = 600, height = 480, units = "px", scale = 2.4)
ggsave("FigureS1.pdf", width = 600, height = 480, units = "px", scale = 2.4)

## Mean accuracy across baseline (ANOVA, Table S1)
df %>% 
  group_by(id, Baseline, BlockCond) %>% 
  summarise(n = n(), Corr = sum(Corr) / n) %>% 
  as.data.frame() %>% 
  select(id,Baseline, BlockCond, Corr) %>% 
  pivot_wider(names_from = Baseline, values_from = Corr) -> df_aov
anovakun(df_aov[, 2:5], "AsB", BlockCond = c("LESS", "MORE"), Baseline = c(50, 100, 150))
p.adjust(c(0.0006, 0.0021, 0.1890), method = "bonferroni")

# Analyses on RT
fit_rt <- lmer(LogChoiceRT ~ Baseline * BlockCond + (Baseline + BlockCond|id), data = df, 
                control = lmerControl(optimizer = "bobyqa"))
summary(fit_rt)
Anova(fit_rt)
coef_rt <- as.data.frame(summary(fit_rt)$coef[-1, ])
p.adjust(coef_rt$Pr, method = "bonferroni")

## Plot RT across baseline (Figure S2a)
group_by(df, Baseline, BlockCond) %>% 
  summarise(M = mean(LogChoiceRT)) %>% 
  as.data.frame() -> df_meanRT
pred_rt <- ggpredict(fit_rt, terms = c("Baseline", "BlockCond"))
pred_rt$group <- factor(pred_rt$group, levels = c("MORE", "LESS"))
ggplot(pred_rt, aes(x, predicted, color = group)) + 
  geom_line(linewidth = 1.4) + geom_point(data = df_meanRT, aes(x = Baseline, y = M, color = BlockCond), shape = 19, size = 3.5) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, color = group), alpha = 0.1) + 
  theme_classic(base_size = 22, base_line_size = 1) + scale_color_grey() +
  xlab("Baseline level") + ylab("Estimated mean logRT") + 
  scale_color_manual(values = c("#080808", "#e0e0e0"), limits = c("MORE", "LESS"), labels = c(MORE = "More", LESS = "Less")) + 
  scale_x_continuous(breaks = seq(0.1, 0.3, 0.1), labels = c(50, 100, 150)) + 
  scale_y_continuous(breaks = seq(-0.5, 0.5, 0.25), limits = c(-0.55, 0.5)) + 
  guides(color = guide_legend("Frame")) + 
  theme(legend.position = c(0.2, 0.99), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16), 
        axis.text = element_text(face = "bold"), 
        plot.margin = unit(c(20, 10, 10, 10), "mm")) -> ps2a

# Model comparison for RT
## More frame
### No random slopes
m1_randi <- lmer(LogChoiceRT ~ LogDiff + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m2_randi <- lmer(LogChoiceRT ~ LogTot + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m3_randi <- lmer(LogChoiceRT ~ LogDiff + LogTot + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m4_randi <- lmer(LogChoiceRT ~ LogDiff * LogTot + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m5_randi <- lmer(LogChoiceRT ~ LogChVal + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m6_randi <- lmer(LogChoiceRT ~ LogUnChVal + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m7_randi <- lmer(LogChoiceRT ~ LogChVal + LogUnChVal + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m8_randi <- lmer(LogChoiceRT ~ LogChVal * LogUnChVal + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))

### Random slopes
m1_rands <- lmer(LogChoiceRT ~ LogDiff + (1 + LogDiff|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m2_rands <- lmer(LogChoiceRT ~ LogTot + (1 + LogTot|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m3_rands <- lmer(LogChoiceRT ~ LogDiff + LogTot + (1 + LogDiff + LogTot|id), 
                 data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m4_rands <- lmer(LogChoiceRT ~ LogDiff * LogTot + (1 + LogDiff + LogTot|id), 
                 data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m5_rands <- lmer(LogChoiceRT ~ LogChVal + (1 + LogChVal|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m6_rands <- lmer(LogChoiceRT ~ LogUnChVal + (1 + LogUnChVal|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m7_rands <- lmer(LogChoiceRT ~ LogChVal + LogUnChVal + (1 + LogChVal + LogUnChVal|id), 
                  data = df_more, control = lmerControl(optimizer = "nlminbwrap"))
m8_rands <- lmer(LogChoiceRT ~ LogChVal * LogUnChVal + (1 + LogChVal + LogUnChVal|id), 
                  data = df_more, control = lmerControl(optimizer = "nlminbwrap"))

### Calculate AIC
AIC_more <- round(AIC(m1_randi, m2_randi, m3_randi, m4_randi, 
                      m5_randi, m6_randi, m7_randi, m8_randi, 
                      m1_rands, m2_rands, m3_rands, m4_rands, 
                      m5_rands, m6_rands, m7_rands, m8_rands), 2)
bestFit_more <- rownames(AIC_more)[which(AIC_more$AIC == min(AIC_more$AIC))]

### Calculate BIC
BIC_more <- round(BIC(m1_randi, m2_randi, m3_randi, m4_randi, 
                      m5_randi, m6_randi, m7_randi, m8_randi, 
                      m1_rands, m2_rands, m3_rands, m4_rands, 
                      m5_rands, m6_rands, m7_rands, m8_rands), 2)
rownames(BIC_more)[which(BIC_more$BIC == min(BIC_more$BIC))]

## Less frame
### No random slopes
l1_randi <- lmer(LogChoiceRT ~ LogDiff + (1|id), data = df_less, 
                  control = lmerControl(optimizer = "nlminbwrap"))
l2_randi <- lmer(LogChoiceRT ~ LogTot + (1|id), data = df_less, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l3_randi <- lmer(LogChoiceRT ~ LogDiff + LogTot + (1|id), data = df_less, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l4_randi <- lmer(LogChoiceRT ~ LogDiff * LogTot + (1|id), data = df_less, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l5_randi <- lmer(LogChoiceRT ~ LogChVal + (1|id), data = df_less, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l6_randi <- lmer(LogChoiceRT ~ LogUnChVal + (1|id), data = df_less, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l7_randi <- lmer(LogChoiceRT ~ LogChVal + LogUnChVal + (1|id), data = df_less, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l8_randi <- lmer(LogChoiceRT ~ LogChVal * LogUnChVal + (1|id), data = df_less, 
                 control = lmerControl(optimizer = "nlminbwrap"))

### Random slopes
l1_rands <- lmer(LogChoiceRT ~ LogDiff + (1 + LogDiff|id), data = df_less, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l2_rands <- lmer(LogChoiceRT ~ LogTot + (1 + LogTot|id), data = df_less, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l3_rands <- lmer(LogChoiceRT ~ LogDiff + LogTot + (1 + LogDiff + LogTot|id), 
                  data = df_less, control = lmerControl(optimizer = "nlminbwrap"))
l4_rands <- lmer(LogChoiceRT ~ LogDiff * LogTot + (1 +  LogDiff + LogTot|id), 
                  data = df_less, control = lmerControl(optimizer = "nlminbwrap"))
l5_rands <- lmer(LogChoiceRT ~ LogChVal + (1 + LogChVal|id), data = df_less, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l6_rands <- lmer(LogChoiceRT ~ LogUnChVal + (1 + LogUnChVal|id), data = df_less, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l7_rands <- lmer(LogChoiceRT ~ LogChVal + LogUnChVal + (1 + LogChVal + LogUnChVal|id), 
                  data = df_less, control = lmerControl(optimizer = "nlminbwrap"))
l8_rands <- lmer(LogChoiceRT ~ LogChVal * LogUnChVal + (1 + LogChVal + LogUnChVal|id), 
                  data = df_less, control = lmerControl(optimizer = "nlminbwrap"))

### Calculate AIC
AIC_less <- round(AIC(l1_randi, l2_randi, l3_randi, l4_randi, 
                      l5_randi, l6_randi, l7_randi, l8_randi,
                      l1_rands, l2_rands, l3_rands, l4_rands,
                      l5_rands, l6_rands, l7_rands, l8_rands), 2)
bestFit_less <- rownames(AIC_less)[which(AIC_less$AIC == min(AIC_less$AIC))]

### Calculate BIC
BIC_less <- round(BIC(l1_randi, l2_randi, l3_randi, l4_randi, 
                      l5_randi, l6_randi, l7_randi, l8_randi,
                      l1_rands, l2_rands, l3_rands, l4_rands,
                      l5_rands, l6_rands, l7_rands, l8_rands), 2)
rownames(BIC_less)[which(BIC_less$BIC == min(BIC_less$BIC))]

# RT prediction of winning models
## `m7_rands` for more frame and and `l7_rands` for less frame
bestFit_more <- eval(parse(text = bestFit_more))
bestFit_less <- eval(parse(text = bestFit_less))
summary(bestFit_more)
coef_bestFit_more <- as.data.frame(summary(bestFit_more)$coef[-1, ])
p.adjust(coef_bestFit_more$Pr, method = "bonferroni")
summary(bestFit_less)
coef_bestFit_less <- as.data.frame(summary(bestFit_less)$coef[-1, ])
p.adjust(coef_bestFit_less$Pr, method = "bonferroni")

## Plot RT prediction (Figure S2b)
data.frame(Coefficient = c(summary(bestFit_more)$coef[1:3, 1], 
                           summary(bestFit_less)$coef[1:3, 1]), 
           SEM = c(summary(bestFit_more)$coef[1:3, 2], 
                   summary(bestFit_less)$coef[1:3, 2])) %>% 
  mutate(Variable = rep(c("Intercept", "Chosen", "Unchosen"), 2), 
         Condition = c(rep("MORE", 3), rep("LESS", 3))) %>% 
  ggplot(aes(x = Variable, y = Coefficient, fill = Condition)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Coefficient - SEM, ymax = Coefficient + SEM), 
                linewidth = 1.1, width = 0.2, position = position_dodge(0.9), 
                color = "#7d7d7d") + 
  scale_x_discrete(limit = c("Chosen", "Unchosen", "Intercept")) + 
  scale_y_continuous(breaks = seq(-1.5, 2, 0.5), limits =  c(-1.5, 1.8)) + 
  scale_fill_manual(values = c("#080808", "#e0e0e0"), 
                    limits = c("MORE", "LESS"), 
                    labels = c(MORE = "More", LESS = "Less")) + 
  theme_classic(base_size = 22, base_line_size = 1) + 
  guides(fill = guide_legend(title = "Frame")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(face = "bold"),
        legend.position = c(0.4, 0.8),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.margin = unit(c(10, 10, 10, 10), "mm")) -> ps2b

## Plot aggregated figures (Figure S2)
plot_grid(ps2a, ps2b, align = "h", labels = c("a", "b"), label_size = 32, nrow = 1)
ggsave("FigureS2.png", width = 1440, height = 840, units = "px", scale = 3.2)
ggsave("FigureS2.pdf", width = 1440, height = 840, units = "px", scale = 3.2)

# Plot mean confidence across baseline (Figure S3) 
df$BlockCond <- factor(df$BlockCond, levels = c("MORE", "LESS"))
df %>% 
  group_by(Baseline, BlockCond) %>% 
  summarise(C = mean(Conf), SD = sd(Conf)) %>% 
  ggplot(aes(x = Baseline, y = C, color = BlockCond, group = BlockCond)) + geom_line() +
  geom_point(size = 4, position = position_dodge()) + 
  geom_errorbar(aes(ymin = C - SD, ymax = C + SD), width = 0.01, position = position_dodge()) + 
  scale_x_continuous(breaks = seq(0.1, 0.3, 0.1), limits = c(0.095, 0.31), labels = c(50, 100, 150)) + 
  xlab("Baseline level") + ylab("Mean confidence") + 
  theme_classic(base_size = 16, base_line_size = 1) + scale_color_grey() +
  scale_color_manual(values = c("#080808", "#e0e0e0"), labels = c("More", "Less")) + 
  guides(color = guide_legend(title = "Frame", order = 1)) + 
  theme(legend.position = c(0.65, 0.8), 
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12), 
        axis.text = element_text(face = "bold")) -> ps3
ggsave("FigureS3.png", width = 600, height = 480, units = "px", scale = 2.4)
ggsave("FigureS3.pdf", width = 600, height = 480, units = "px", scale = 2.4)

# Regrssion analysis using "diff-(chosen-unchosen)" rule
## More frame
more_dch <- lmer(Conf ~ LogDiff + LogChVal + LogUnChVal + 
                   LogDiff:LogChVal + LogDiff:LogUnChVal + LogChVal:LogUnChVal + 
                   (1 + LogDiff + LogChVal + LogUnChVal|id), 
                   data = df_more, REML = F, 
                   control = lmerControl(optimizer = "bobyqa"))
summary(more_dch)
Anova(more_dch)
coef_bestFit2_more <- as.data.frame(summary(more_dch)$coef[-1, ])
p.adjust(Anova(more_dch)$Pr, method = "bonferroni")

## Less frame
less_dch <- lmer(Conf ~ LogDiff + LogChVal + LogUnChVal + 
                   LogDiff:LogChVal + LogDiff:LogUnChVal + LogChVal:LogUnChVal + 
                   (1 + LogDiff + LogChVal + LogUnChVal|id), 
                   data = df_less, REML = F, 
                   control = lmerControl(optimizer = "bobyqa"))
summary(less_dch)
Anova(less_dch)
coef_bestFit2_less <- as.data.frame(summary(less_dch)$coef[-1, ])
p.adjust(Anova(less_dch)$Pr, method = "bonferroni")

## Plot estimated parameters (Figure S4)
data.frame(Coefficient = c(summary(more_dch)$coef[, 1],
                           summary(less_dch)$coef[, 1]), 
           SEM = c(summary(more_dch)$coef[, 2], 
                   summary(less_dch)$coef[, 2])) %>% 
  mutate(Variable = rep(c("Intercept", "Diff", "Chosen", "Unchosen", "Diff*Chosen", "Diff*Unchosen", "Chosen*Unchosen"), 2), 
         Condition = c(rep("MORE", 7), rep("LESS", 7))) %>% 
  ggplot(aes(x = Variable, y = Coefficient, fill = Condition)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Coefficient - SEM, ymax = Coefficient + SEM), 
                linewidth = 1.1, width = 0.2, position = position_dodge(0.9), 
                color = "#7d7d7d") + 
  scale_x_discrete(limit = c("Diff", "Chosen", "Unchosen", "Diff*Chosen", "Diff*Unchosen", "Chosen*Unchosen", "Intercept")) + 
  scale_y_continuous(breaks = seq(-1.5, 1.5, 0.5), limits =  c(-1.2, 1.2)) + 
  scale_fill_manual(values = c("#080808", "#e0e0e0"), 
                    limits = c("MORE", "LESS"), 
                    labels = c(MORE = "More", LESS = "Less")) + 
  theme_classic(base_size = 22, base_line_size = 1) + 
  guides(fill = guide_legend(title = "Frame")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(face = "bold"),
        legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.margin = unit(c(10, 10, 10, 10), "mm")) -> ps4d
ggsave("FigureS4.png", width = 600, height = 480, units = "px", scale = 4.8)
ggsave("FigureS4.pdf", width = 600, height = 480, units = "px", scale = 4.8)

# Best fit models under the diff-sum rule (Figure S5)
best_ds_more <- lmer(Conf ~ LogDiff * LogTot + (1 +  LogDiff + LogTot|id), 
                     data = df_more, REML = F, 
                     control = lmerControl(optimizer = "nlminbwrap"))
best_ds_less <- lmer(Conf ~ LogDiff * LogTot + (1 +  LogDiff + LogTot|id), 
                     data = df_less, REML = F, 
                     control = lmerControl(optimizer = "nlminbwrap"))
data.frame(Coefficient = c(summary(best_ds_more)$coef[1:4, 1], 
                           summary(best_ds_less)$coef[1:4, 1]), 
           SEM = c(summary(best_ds_more)$coef[1:4, 2], 
                   summary(best_ds_less)$coef[1:4, 2])) %>% 
  mutate(Variable = rep(c("Intercept", "Diff", "Sum", "Diff*Sum"), 2), 
         Condition = c(rep("MORE", 4), rep("LESS", 4))) %>% 
  ggplot(aes(x = Variable, y = Coefficient, fill = Condition)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Coefficient - SEM, ymax = Coefficient + SEM), 
                linewidth = 1.1, width = 0.2, position = position_dodge(0.9), 
                color = "#7d7d7d") + 
  scale_x_discrete(limit = c("Diff", "Sum", "Diff*Sum", "Intercept")) + 
  scale_y_continuous(breaks = seq(-1.5, 1.5, 0.5), limits =  c(-1.5, 1.5)) + 
  scale_fill_manual(values = c("#080808", "#e0e0e0"), 
                    limits = c("MORE", "LESS"), 
                    labels = c(MORE = "More", LESS = "Less")) + 
  theme_classic(base_size = 22, base_line_size = 1) + 
  guides(fill = guide_legend(title = "Frame")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(face = "bold"),
        legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.margin = unit(c(10, 10, 10, 10), "mm")) -> ps5
ggsave("FigureS5.png", width = 600, height = 480, units = "px", scale = 3.6)
ggsave("FigureS5.pdf", width = 600, height = 480, units = "px", scale = 3.6)