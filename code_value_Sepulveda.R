# Import library
library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(cowplot)
library(ggeffects)

# Import data (variables are standardized)
df <- read.csv("DataFoodFramingNotebook.csv", header = TRUE)
df <- select(df, LValue, RValue, ValCh, ValUnCh, ChosenITM, Correct, ChoiceRT, Conf, ConfRT, BlockCond, Part)
colnames(df) <- c("LVal", "RVal", "ChVal",  "UnChVal", "ChosenITM", "Corr", "ChoiceRT", "Conf", "ConfRT", "BlockCond", "id")
df$BlockCond <- if_else(df$BlockCond == 1, "like", "dislike")
df <- mutate(df, DiffRL = RVal - LVal)
df$DiffVal <- abs(df$LVal - df$RVal)
df <- subset(df, DiffVal != 0)
df$TotVal <- df$LVal + df$RVal
df$DiffVal <- scale(df$DiffVal)
df$TotVal <- scale(df$TotVal)
df$ChVal <- scale(df$ChVal)
df$UnChVal <- scale(df$UnChVal)
df$Conf <- scale(df$Conf)

df_like <- filter(df, BlockCond == "like")
df_dislike <- filter(df, BlockCond == "dislike")
################################################################################
# Choice probability
fit_choice <- glmer(ChosenITM ~ DiffRL * BlockCond + 
                      (DiffRL + BlockCond|id), data = df, 
                    family = "binomial", 
                    control = glmerControl(optimizer = "bobyqa"))
summary(fit_choice)
Anova(fit_choice)
coef_choice <- as.data.frame(summary(fit_choice)$coef[-1, ])
p.adjust(coef_choice$Pr, method = "bonferroni")

# Accuracy
fit_corr <- glmer(Corr ~ DiffVal + TotVal + BlockCond + 
                    DiffVal:TotVal + TotVal:BlockCond + BlockCond:DiffVal + 
                    (DiffVal + TotVal + BlockCond|id), data = df, 
                  family = "binomial", 
                  control = glmerControl(optimizer = "bobyqa"))
summary(fit_corr)
Anova(fit_corr)
p.adjust(Anova(fit_corr)$Pr, method = "bonferroni")

# Confidence
fit_conf <- lmer(Conf ~ TotVal * BlockCond +(TotVal + BlockCond|id), 
                 data = df, REML = F, control = lmerControl(optimizer = "bobyqa"))
summary(fit_conf)
Anova(fit_conf)
coef_conf <- as.data.frame(summary(fit_conf)$coef[-1, ])
p.adjust(coef_conf$Pr, method = "bonferroni")

# Model comparison
## Like frame
### No random slopes
m1_randi <- lmer(Conf ~ DiffVal + (1|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m2_randi <- lmer(Conf ~ TotVal + (1|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m3_randi <- lmer(Conf ~ DiffVal + TotVal + (1|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m4_randi <- lmer(Conf ~ DiffVal * TotVal + (1|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m5_randi <- lmer(Conf ~ ChVal + (1|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m6_randi <- lmer(Conf ~ UnChVal + (1|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m7_randi <- lmer(Conf ~ ChVal + UnChVal + (1|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m8_randi <- lmer(Conf ~ ChVal * UnChVal + (1|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))

### Random slopes
m1_rands <- lmer(Conf ~ DiffVal + (1 + DiffVal|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m2_rands <- lmer(Conf ~ TotVal + (1 + TotVal|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m3_rands <- lmer(Conf ~ DiffVal + TotVal + (1 + DiffVal + TotVal|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m4_rands <- lmer(Conf ~ DiffVal * TotVal + (1 +  DiffVal + TotVal|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m5_rands <- lmer(Conf ~ ChVal + (1 + ChVal|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m6_rands <- lmer(Conf ~ UnChVal + (1 + UnChVal|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m7_rands <- lmer(Conf ~ ChVal + UnChVal + (1 + ChVal + UnChVal|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
m8_rands <- lmer(Conf ~ ChVal * UnChVal + (1 + ChVal + UnChVal|id), data = df_like, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))

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
l1_randi <- lmer(Conf ~ DiffVal + (1|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l2_randi <- lmer(Conf ~ TotVal + (1|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l3_randi <- lmer(Conf ~ DiffVal + TotVal + (1|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l4_randi <- lmer(Conf ~ DiffVal * TotVal + (1|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l5_randi <- lmer(Conf ~ ChVal + (1|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l6_randi <- lmer(Conf ~ UnChVal + (1|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l7_randi <- lmer(Conf ~ ChVal + UnChVal + (1|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l8_randi <- lmer(Conf ~ ChVal * UnChVal + (1|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))

### Random slopes
l1_rands <- lmer(Conf ~ DiffVal + (1 + DiffVal|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l2_rands <- lmer(Conf ~ TotVal + (1 + TotVal|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l3_rands <- lmer(Conf ~ DiffVal + TotVal + (1 + DiffVal + TotVal|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l4_rands <- lmer(Conf ~ DiffVal * TotVal + (1 +  DiffVal + TotVal|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l5_rands <- lmer(Conf ~ ChVal + (1 + ChVal|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l6_rands <- lmer(Conf ~ UnChVal + (1 + UnChVal|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l7_rands <- lmer(Conf ~ ChVal + UnChVal + (1 + ChVal + UnChVal|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))
l8_rands <- lmer(Conf ~ ChVal * UnChVal + (1 + ChVal + UnChVal|id), data = df_dislike, REML = F, 
                 control = lmerControl(optimizer = "bobyqa"))

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
## `m8_rands` for like frame and and `l7_rands` for dislike frame
bestFit_like <- eval(parse(text = bestFit_like))
bestFit_dislike <- eval(parse(text = bestFit_dislike))
summary(bestFit_like)
coef_bestFit_like <- as.data.frame(summary(bestFit_like)$coef[-1, ])
p.adjust(coef_bestFit_like$Pr, method = "bonferroni")
summary(bestFit_dislike)
coef_bestFit_dislike <- as.data.frame(summary(bestFit_dislike)$coef[-1, ])
p.adjust(coef_bestFit_dislike$Pr, method = "bonferroni")
################################################################################
# Plot figures
## Choice probability (Figure 4a)
coef <- summary(fit_choice)$coef[, 1]
choice_pred <- function(x, c){
  return(1 / (1 + exp(-(coef[1] + coef[2] * x + coef[3] * c + coef[4] * x * c))))
}
df_choice_pred <- data.frame(x = runif(10000, -4.5, 4.5), 
                             c = c(rep(1, 5000), rep(0, 5000)))
df_choice_pred <- mutate(df_choice_pred, pred = choice_pred(x, c), 
                         BlockCond = ifelse(c == 0, "Like", "Dislike"))

df %>%  
  mutate(DiffRL = round(DiffRL, 4)) %>% 
  select(DiffRL, ChosenITM, BlockCond) %>% 
  group_by(DiffRL, BlockCond) %>% 
  summarise(p = sum(ChosenITM) / n(), n = n()) %>% 
  as.data.frame() -> df_choice_mean
df_choice_mean <- filter(df_choice_mean, n > 5)

ggplot(df_choice_pred, aes(x = x, y = pred, group = BlockCond)) + 
  geom_line(aes(linetype = BlockCond), linewidth = 1.3) + 
  geom_point(df_choice_mean, mapping = aes(x = DiffRL, y = p, color = BlockCond), size = 2.5) + 
  scale_x_continuous(breaks = seq(-4, 4, 1), limits = c(-4, 4), 
                     expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) + 
  scale_linetype_discrete(labels = c("Like", "Dislike")) + 
  scale_color_manual(values = c("#080808", "#e0e0e0"), limits = c("like", "dislike"), labels = c(like = "Like", dislike = "Dislike")) + 
  xlab("Value difference (right - left)") + 
  ylab("Probability of choosing right option") + 
  theme_classic(base_size = 22, base_line_size = 1) + 
  guides(linetype = guide_legend(title = "Frame", order = 1), color = guide_legend(title = "", order = 2), alpha = "none") + 
  theme(legend.position = c(0.95, 0.5), 
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16), 
        axis.text = element_text(face = "bold"), 
        plot.margin = unit(c(20, 10, 10, 10), "mm")) -> p1

## Confidence (Figure 4b)
pred_conf <- ggpredict(fit_conf, terms = c("TotVal", "BlockCond"))
pred_conf$group <- factor(pred_conf$group, levels = c("like", "dislike"))
ggplot(pred_conf, aes(x, predicted, color = group)) + 
  geom_line(linewidth = 1.4) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, color = group), alpha = 0.1) + 
  theme_classic(base_size = 22, base_line_size = 1) + scale_color_grey() +
  xlab("Summed value") + ylab("Estimated mean confidence") + 
  scale_color_manual(values = c("#080808", "#e0e0e0"), limits = c("like", "dislike"), labels = c(like = "Like", dislike = "Dislike")) + 
  scale_x_continuous(breaks = seq(-2, 3, 1), limits = c(-2, 3)) + 
  guides(color = guide_legend("Frame")) + 
  theme(legend.position = c(0.4, 0.2), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16), 
        axis.text = element_text(face = "bold"), 
        plot.margin = unit(c(20, 10, 10, 10), "mm")) -> p2

## AIC (Figure 4c)
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
  scale_y_continuous(breaks = seq(800, 9000, 200), limits = c(8000, 9000)) + 
  scale_shape_manual(values = c(5, 8), labels = c(ch = "Chosen-unchosen", ds = "Diff-sum")) + 
  guides(shape = guide_legend(title = "Rule")) + 
  xlab("Frame") + ylab("AIC") + 
  theme_classic(base_size = 22, base_line_size = 1) + 
  theme(axis.text = element_text(face = "bold"), 
        legend.position = c(0.3, 0.85),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.margin = unit(c(10, 10, 10, 10), "mm")) -> p3

## Hierarchical linear regressions (Figure 4d)
data.frame(Coefficient = c(summary(bestFit_like)$coef[1:4, 1], 
                           summary(bestFit_dislike)$coef[1:3, 1], 0), 
           SEM = c(summary(bestFit_like)$coef[1:4, 2], 
                   summary(bestFit_dislike)$coef[1:3, 2], 0)) %>% 
  mutate(Variable = rep(c("Intercept", "Chosen", "Unchosen", "Chosen*Unchosen"), 2), 
         Condition = c(rep("like", 4), rep("dislike", 4))) %>% 
  ggplot(aes(x = Variable, y = Coefficient, fill = Condition)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Coefficient - SEM, ymax = Coefficient + SEM), 
                linewidth = 1.1, width = 0.2, position = position_dodge(0.9), 
                color = "#7d7d7d") + 
  scale_x_discrete(limit = c("Chosen", "Unchosen", "Chosen*Unchosen", "Intercept")) + 
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
        plot.margin = unit(c(10, 10, 10, 10), "mm")) -> p4

## Aggregate multiple plots (Figure 4)
plot_grid(plot_grid(p1, p2, labels = c("a", "b"), label_size = 32), 
          plot_grid(p3, p4, align = "h", labels = c("c", "d"), label_size = 32), nrow = 2)
ggsave("Figure4.png", width = 1440, height = 1200, units = "px", scale = 3.2)
ggsave("Figure4.pdf", width = 1440, height = 1200, units = "px", scale = 3.2)

## Accuracy (Appendix C)
pred_corr <- ggpredict(fit_corr, terms = c("DiffVal", "BlockCond", "TotVal"))
pred_corr$group <- factor(pred_corr$group, levels = c("like", "dislike"))
ggplot(pred_corr, aes(x = x, y = predicted, group = group)) + 
  geom_line(linewidth = 1.4, aes(linetype = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, color = group), alpha = 0.1) + 
  facet_wrap(.~facet) + 
  theme_classic(base_size = 22, base_line_size = 1) + scale_color_grey() +
  xlab("Value difference") + ylab("") + 
  scale_y_continuous(breaks = seq(0.4, 1, 0.2), limits = c(0.3, 1)) + 
  scale_linetype_manual(name = "", values = c(1, 2), labels = c("Like", "Dislike")) + 
  guides(linetype = guide_legend("Frame"), color = "none") + 
  theme(legend.position = c(0.5, 0.2), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16), 
        axis.text = element_text(face = "bold")) -> p_appx2