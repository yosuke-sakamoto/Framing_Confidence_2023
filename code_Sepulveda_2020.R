# Import library
library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(ggeffects)

# Import data
df <- read.csv("DataFoodFramingNotebook.csv", header = TRUE)
df <- select(df, LValue, RValue, ValCh, ValUnCh, ChosenITM, Correct, ChoiceRT, Conf, ConfRT, BlockCond, Part)
colnames(df) <- c("LVal", "RVal", "DCE",  "DICE", "ChosenITM", "Corr", "ChoiceRT", "Conf", "ConfRT", "BlockCond", "id")
df$BlockCond <- if_else(df$BlockCond == 1, "like", "dislike")
df$DiffVal <- abs(df$LVal - df$RVal)
df$TotVal <- df$LVal + df$RVal
df$TotVal <- round(df$TotVal)
df$zDiffVal <- scale(df$DiffVal)
df$zTotVal <- scale(df$TotVal)
df$zDCE <- scale(df$DCE)
df$zDICE <- scale(df$DICE)
df$zChoiceRT <- scale(df$ChoiceRT)
df$zConf <- scale(df$Conf)

df$zTotVal <- split_quantile(df$zTotVal, 2)

df_like <- filter(df, BlockCond == "like")
df_dislike <- filter(df, BlockCond == "dislike")
################################################################################
# Confidence
fit_conf <- lmer(zConf ~ zTotVal + zDiffVal + BlockCond +  
                   zTotVal:zDiffVal + zDiffVal:BlockCond + BlockCond:zTotVal + 
                   (zTotVal + zDiffVal + BlockCond|id), 
                 data = df, REML = F, control = lmerControl(optimizer = "bobyqa"))
summary(fit_conf)
Anova(fit_conf)
plot(ggpredict(fit_conf, terms= c("zTotVal", "zDiffVal", "BlockCond")))
################################################################################
# Model comparison
## Like frame
### No random slopes
m1_randi <- lmer(zConf ~ zDiffVal + (1|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m2_randi <- lmer(zConf ~ zTotVal + (1|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m3_randi <- lmer(zConf ~ zDiffVal + zTotVal + (1|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m4_randi <- lmer(zConf ~ zDiffVal * zTotVal + (1|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m5_randi <- lmer(zConf ~ zDCE + (1|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m6_randi <- lmer(zConf ~ zDICE + (1|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m7_randi <- lmer(zConf ~ zDCE + zDICE + (1|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m8_randi <- lmer(zConf ~ zDCE * zDICE + (1|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))

### Random slopes
m1_rands <- lmer(zConf ~ zDiffVal + (1 + zDiffVal|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m2_rands <- lmer(zConf ~ zTotVal + (1 + zTotVal|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m3_rands <- lmer(zConf ~ zDiffVal + zTotVal + (1 + zDiffVal + zTotVal|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m4_rands <- lmer(zConf ~ zDiffVal * zTotVal + (1 +  zDiffVal + zTotVal|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m5_rands <- lmer(zConf ~ zDCE + (1 + zDCE|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m6_rands <- lmer(zConf ~ zDICE + (1 + zDICE|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m7_rands <- lmer(zConf ~ zDCE + zDICE + (1 + zDCE + zDICE|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))
m8_rands <- lmer(zConf ~ zDCE * zDICE + (1 + zDCE + zDICE|id), data = df_like, REML = F, control = lmerControl(optimizer = "bobyqa"))

### Calculate AIC
AIC_like <- round(AIC(m1_randi, m2_randi, m3_randi, m4_randi, 
                      m5_randi, m6_randi, m7_randi, m8_randi, 
                      m1_rands, m2_rands, m3_rands, m4_rands, 
                      m5_rands, m6_rands, m7_rands, m8_rands), 2)
bestFit_like <- rownames(AIC_like)[which(AIC_like$AIC == min(AIC_like$AIC))]

### Calculate BIC
BIC_like <- round(BIC(m1_randi, m2_randi, m3_randi, m4_randi, 
                      m5_randi, m6_randi, m7_randi, m8_randi, 
                      m1_rands, m2_rands, m3_rands, m4_rands, 
                      m5_rands, m6_rands, m7_rands, m8_rands), 2)
rownames(BIC_like)[which(BIC_like$BIC == min(BIC_like$BIC))]

## Dislike frame
### No random slopes
l1_randi <- lmer(zConf ~ zDiffVal + (1|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l2_randi <- lmer(zConf ~ zTotVal + (1|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l3_randi <- lmer(zConf ~ zDiffVal + zTotVal + (1|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l4_randi <- lmer(zConf ~ zDiffVal * zTotVal + (1|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l5_randi <- lmer(zConf ~ zDCE + (1|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l6_randi <- lmer(zConf ~ zDICE + (1|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l7_randi <- lmer(zConf ~ zDCE + zDICE + (1|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l8_randi <- lmer(zConf ~ zDCE * zDICE + (1|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))

### Random slopes
l1_rands <- lmer(zConf ~ zDiffVal + (1 + zDiffVal|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l2_rands <- lmer(zConf ~ zTotVal + (1 + zTotVal|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l3_rands <- lmer(zConf ~ zDiffVal + zTotVal + (1 + zDiffVal + zTotVal|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l4_rands <- lmer(zConf ~ zDiffVal * zTotVal + (1 +  zDiffVal + zTotVal|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l5_rands <- lmer(zConf ~ zDCE + (1 + zDCE|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l6_rands <- lmer(zConf ~ zDICE + (1 + zDICE|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l7_rands <- lmer(zConf ~ zDCE + zDICE + (1 + zDCE + zDICE|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))
l8_rands <- lmer(zConf ~ zDCE * zDICE + (1 + zDCE + zDICE|id), data = df_dislike, REML = F, control = lmerControl(optimizer = "bobyqa"))

### Calculate AIC
AIC_dislike <- round(AIC(l1_randi, l2_randi, l3_randi, l4_randi, 
                      l5_randi, l6_randi, l7_randi, l8_randi, 
                      l1_rands, l2_rands, l3_rands, l4_rands, 
                      l5_rands, l6_rands, l7_rands, l8_rands), 2)
bestFit_dislike <- rownames(AIC_dislike)[which(AIC_dislike$AIC == min(AIC_dislike$AIC))]

### Calculate BIC
BIC_dislike <- round(BIC(l1_randi, l2_randi, l3_randi, l4_randi, 
                      l5_randi, l6_randi, l7_randi, l8_randi, 
                      l1_rands, l2_rands, l3_rands, l4_rands, 
                      l5_rands, l6_rands, l7_rands, l8_rands), 2)
rownames(BIC_dislike)[which(BIC_dislike$BIC == min(BIC_dislike$BIC))]

# Confidence prediction of winning models
## `m8_rands` for more frame and and `l7_rands` for less frame
bestFit_like <- eval(parse(text = bestFit_like))
bestFit_dislike <- eval(parse(text = bestFit_dislike))
summary(bestFit_like)
summary(bestFit_dislike)
################################################################################
# Plot figures
## Confidence
pred_conf <- ggpredict(fit_conf, terms = c("zTotVal", "zDiffVal", "BlockCond"))
pred_conf$facet <- factor(pred_conf$facet, levels = c("like", "dislike"))
ggplot(pred_conf, aes(x, predicted, color = group)) + 
  geom_point(position = position_jitterdodge(seed = 1), size = 4) + 
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), linewidth = 1.5, position = position_jitterdodge(seed = 1)) + 
  facet_wrap(~ facet, labeller = labeller(facet = c(like = "Like", dislike = "Dislike"))) + 
  theme_classic(base_size = 22, base_line_size = 1) + scale_color_grey() +
  xlab("Summed value") + ylab("Predicted confidence") + 
  scale_x_discrete(labels = c("Low", "High")) + 
  guides(color = guide_legend("Value difference")) + 
  theme(legend.position = c(0.4, 0.2), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16), 
        axis.text = element_text(face = "bold"), 
        plot.margin = unit(c(20, 10, 10, 10), "mm")) -> p1

## AIC
AIC_diff_sum_like <- round(AIC(m1_randi, m2_randi, m3_randi, m4_randi, 
                               m1_rands, m2_rands, m3_rands, m4_rands), 2)
AIC_diff_sum_dislike <- round(AIC(l1_randi, l2_randi, l3_randi, l4_randi, 
                                  l1_rands, l2_rands, l3_rands, l4_rands), 2)
AIC_chosen_unchosen_like <- round(AIC(m5_randi, m6_randi, m7_randi, m8_randi, 
                                      m5_rands, m6_rands, m7_rands, m8_rands), 2)
AIC_chosen_unchosen_dislike <- round(AIC(l5_randi, l6_randi, l7_randi, l8_randi,
                                         l5_rands, l6_rands, l7_rands, l8_rands), 2)

rbind(AIC_diff_sum_like[which(AIC_diff_sum_like$AIC == min(AIC_diff_sum_like$AIC)), ], 
      AIC_diff_sum_dislike[which(AIC_diff_sum_dislike$AIC == min(AIC_diff_sum_dislike$AIC)), ], 
      AIC_chosen_unchosen_like[which(AIC_chosen_unchosen_like$AIC == min(AIC_chosen_unchosen_like$AIC)), ], 
      AIC_chosen_unchosen_dislike[which(AIC_chosen_unchosen_dislike$AIC == min(AIC_chosen_unchosen_dislike$AIC)), ]) %>% 
  mutate(frame = rep(c("like", "dislike"), 2), 
         rule = c(rep("ds", 2), rep("ch", 2))) %>% 
  ggplot(aes(x = frame, y = AIC)) + 
  geom_point(aes(shape = rule), size = 5) + 
  scale_x_discrete(limits = c("like", "dislike"), labels = c(like = "Like", dislike = "Dislike")) + 
  scale_y_continuous(breaks = seq(9200, 10200, 200), limits = c(9200, 10200)) + 
  scale_shape_manual(values = c(5, 8), labels = c(ch = "Chosen-unchosen", ds = "Diff-sum")) + 
  guides(shape = guide_legend(title = "Rule")) + 
  xlab("Frame") + ylab("AIC") + 
  theme_classic(base_size = 22, base_line_size = 1) + 
  theme(axis.text = element_text(face = "bold"), 
        legend.position = c(0.3, 0.85),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.margin = unit(c(30, 10, 10, 10), "mm")) -> p2

## Hierarchical linear regressions
data.frame(Coefficient = c(summary(bestFit_like)$coef[1:4, 1], 
                           summary(bestFit_dislike)$coef[1:3, 1], 0), 
           SEM = c(summary(bestFit_like)$coef[1:4, 2], 
                   summary(bestFit_dislike)$coef[1:3, 2], 0)) %>% 
  mutate(Variable = rep(c("Intercept", "Chosen", "Unchosen", "Interaction"), 2), 
         Condition = c(rep("like", 4), rep("dislike", 4))) -> df_coef
ggplot(df_coef, aes(x = Variable, y = Coefficient, fill = Condition)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Coefficient - SEM, ymax = Coefficient + SEM), 
                linewidth = 1.1, width = 0.2, position = position_dodge(0.9), 
                color = "#7d7d7d") + 
  scale_x_discrete(limit = c("Chosen", "Unchosen", "Interaction", "Intercept")) + 
  scale_y_continuous(breaks = seq(-1.0, 1.0, 0.5), limits =  c(-1.0, 1.0)) + 
  scale_fill_manual(values = c("#080808", "#e0e0e0"), 
                    limits = c("like", "dislike"), 
                    labels = c(like = "Like", dislike = "Dislike")) + 
  theme_classic(base_size = 22, base_line_size = 1) + 
  guides(fill = guide_legend(title = "Frame")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(face = "bold"),
        legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.margin = unit(c(10, 80, 10, 60), "mm")) -> p3

plot_grid(plot_grid(p1, p2, labels = c("a", "b"), label_size = 32), p3, labels = c("", "c"), label_size = 32, nrow = 2)
ggsave("figure4.png", width = 1440, height = 1200, units = "px", scale = 3.2)
