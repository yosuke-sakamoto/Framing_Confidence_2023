# Import library
library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(cowplot)
library(ggeffects)

# Import data
df <- read.csv("data_behavior.csv", header = TRUE)
df$LogDiff <- scale(df$LogDiff)
df$LogTot <- scale(df$LogTot)
df$LogDCE <- scale(df$LogDCE)
df$LogDICE <- scale(df$LogDICE)
df$Conf <- scale(df$Conf)
df <- mutate(df, DiffRL = RVal - LVal, LogDiffRL = LogRVal - LogLVal, 
             Baseline = if_else(LVal == 50 | RVal == 50, 1, if_else(LVal == 100 | RVal == 100, 2, 3)))
df$Baseline <- as.factor(df$Baseline)
df_more <- filter(df, BlockCond == "MORE")
df_less <- filter(df, BlockCond == "LESS")
################################################################################
# Choice probability
fit_choice <- glmer(ChosenITM ~ LogDiffRL * BlockCond + 
                      (LogDiffRL + BlockCond|id), data = df, 
                    family = "binomial", 
                    control = glmerControl(optimizer = "bobyqa"))
summary(fit_choice)
Anova(fit_choice)

# Accuracy between conditions
df %>% 
  group_by(id, BlockCond) %>% 
  summarise(n = n(), prop = sum(Corr) / n) -> acc
t.test(filter(acc, BlockCond == "MORE")$prop, 
       filter(acc, BlockCond == "LESS")$prop, var.equal = T)
################################################################################
# Confidence
fit_conf <- lmer(Conf ~ Baseline + LogDiff + BlockCond +  
                   Baseline:LogDiff + LogDiff:BlockCond + BlockCond:Baseline + 
                   (Baseline + LogDiff + BlockCond|id), 
                 data = df, REML = F, control = lmerControl(optimizer = "bobyqa"))
summary(fit_conf)
Anova(fit_conf)
################################################################################
# Model comparison
## More frame
### No random slopes
m1_randi <- lmer(Conf ~ LogDiff + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m2_randi <- lmer(Conf ~ LogTot + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m3_randi <- lmer(Conf ~ LogDiff + LogTot + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m4_randi <- lmer(Conf ~ LogDiff * LogTot + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m5_randi <- lmer(Conf ~ LogDCE + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m6_randi <- lmer(Conf ~ LogDICE + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m7_randi <- lmer(Conf ~ LogDCE + LogDICE + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m8_randi <- lmer(Conf ~ LogDCE * LogDICE + (1|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))

### Random slopes
m1_rands <- lmer(Conf ~ LogDiff + (1 + LogDiff|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m2_rands <- lmer(Conf ~ LogTot + (1 + LogTot|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m3_rands <- lmer(Conf ~ LogDiff + LogTot + (1 + LogDiff + LogTot|id), 
                 data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m4_rands <- lmer(Conf ~ LogDiff * LogTot + (1 +  LogDiff + LogTot|id), 
                 data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m5_rands <- lmer(Conf ~ LogDCE + (1 + LogDCE|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m6_rands <- lmer(Conf ~ LogDICE + (1 + LogDICE|id), data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m7_rands <- lmer(Conf ~ LogDCE + LogDICE + (1 + LogDCE + LogDICE|id), 
                 data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
m8_rands <- lmer(Conf ~ LogDCE * LogDICE + (1 + LogDCE + LogDICE|id), 
                 data = df_more, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))

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
l1_randi <- lmer(Conf ~ LogDiff + (1|id), data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l2_randi <- lmer(Conf ~ LogTot + (1|id), data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l3_randi <- lmer(Conf ~ LogDiff + LogTot + (1|id), data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l4_randi <- lmer(Conf ~ LogDiff * LogTot + (1|id), data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l5_randi <- lmer(Conf ~ LogDCE + (1|id), data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l6_randi <- lmer(Conf ~ LogDICE + (1|id), data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l7_randi <- lmer(Conf ~ LogDCE + LogDICE + (1|id), data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l8_randi <- lmer(Conf ~ LogDCE * LogDICE + (1|id), data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))

### Random slopes
l1_rands <- lmer(Conf ~ LogDiff + (1 + LogDiff|id), data = df_less, 
                 REML = F, control = lmerControl(optimizer = "nlminbwrap"))
l2_rands <- lmer(Conf ~ LogTot + (1 + LogTot|id), data = df_less, 
                 REML = F, control = lmerControl(optimizer = "nlminbwrap"))
l3_rands <- lmer(Conf ~ LogDiff + LogTot + (1 + LogDiff + LogTot|id), 
                 data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l4_rands <- lmer(Conf ~ LogDiff * LogTot + (1 +  LogDiff + LogTot|id), 
                 data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l5_rands <- lmer(Conf ~ LogDCE + (1 + LogDCE|id), data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l6_rands <- lmer(Conf ~ LogDICE + (1 + LogDICE|id), data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l7_rands <- lmer(Conf ~ LogDCE + LogDICE + (1 + LogDCE + LogDICE|id), 
                 data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))
l8_rands <- lmer(Conf ~ LogDCE * LogDICE + (1 + LogDCE + LogDICE|id), 
                 data = df_less, REML = F, 
                 control = lmerControl(optimizer = "nlminbwrap"))

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

# Confidence prediction of winning models
## `m8_rands` for more frame and and `l8_rands` for less frame
bestFit_more <- eval(parse(text = bestFit_more))
bestFit_less <- eval(parse(text = bestFit_less))
summary(bestFit_more)
summary(bestFit_less)
################################################################################
# Plot figures
## Choice probability
coef <- summary(fit_choice)$coef[, 1]
choice_pred <- function(x, c){
  return(1 / (1 + exp(-(coef[1] + coef[2] * x + coef[3] * c + coef[4] * x * c))))
}
df_choice_pred <- data.frame(x = runif(10000, -0.4, 0.4), 
                             c = c(rep(1, 5000), rep(0, 5000)))
df_choice_pred <- mutate(df_choice_pred, pred = choice_pred(x, c), 
                         BlockCond = ifelse(c == 0, "More", "Less"))

df %>% 
  mutate(LogDiffRL = round(LogDiffRL, 4)) %>% 
  select(LogDiffRL, ChosenITM, BlockCond) %>% 
  group_by(LogDiffRL, BlockCond) %>% 
  summarise(p = sum(ChosenITM) / n()) %>% 
  as.data.frame() -> df_choice_mean

ggplot(df_choice_pred, aes(x = x, y = pred, group = BlockCond)) + 
  geom_line(aes(linetype = BlockCond), linewidth = 1.3) + 
  geom_point(df_choice_mean, mapping = aes(x = LogDiffRL, y = p, color = BlockCond), size = 2.5) + 
  scale_x_continuous(breaks = seq(-0.2, 0.2, 0.1), limits = c(-0.25, 0.25), 
                     expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  scale_linetype_discrete(labels = c("More", "Less")) + 
  scale_color_manual(values = c("#080808", "#e0e0e0"), limits = c("MORE", "LESS"), labels = c(MORE = "More", LESS = "Less")) + 
  xlab("Log-transformed value difference (right - left)") + 
  ylab("Probability of choosing right option") + 
  theme_classic(base_size = 22, base_line_size = 1) + 
  guides(linetype = guide_legend(title = "Frame", order = 1), color = guide_legend(title = "", order = 2), alpha = "none") + 
  theme(legend.position = c(0.95, 0.5), 
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16), 
        axis.text = element_text(face = "bold"), 
        plot.margin = unit(c(20, 10, 10, 10), "mm")) -> p1

## Confidence
pred_conf <- ggpredict(fit_conf, terms = c("Baseline", "LogDiff", "BlockCond"))
pred_conf$facet <- factor(pred_conf$facet, levels = c("MORE", "LESS"))
ggplot(pred_conf, aes(x, predicted, color = group)) + 
  geom_point(position = position_jitterdodge(seed = 1), size = 4) + 
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), linewidth = 1.5, position = position_jitterdodge(seed = 1)) + 
  facet_wrap(~ facet, labeller = labeller(facet = c(MORE = "More", LESS = "Less"))) + 
  theme_classic(base_size = 22, base_line_size = 1) + scale_color_grey() +
  xlab("Baseline level") + ylab("Predicted confidence") + 
  scale_x_discrete(labels = c("50", "100" ,"150")) + 
  guides(color = guide_legend("Value difference")) + 
  theme(legend.position = c(0.4, 0.2), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16), 
        axis.text = element_text(face = "bold"), 
        plot.margin = unit(c(20, 10, 10, 10), "mm")) -> p2

## AIC
AIC_diff_sum_more <- round(AIC(m1_randi, m2_randi, m3_randi, m4_randi, 
                               m1_rands, m2_rands, m3_rands, m4_rands), 2)
AIC_diff_sum_less <- round(AIC(l1_randi, l2_randi, l3_randi, l4_randi, 
                               l1_rands, l2_rands, l3_rands, l4_rands), 2)
AIC_chosen_unchosen_more <- round(AIC(m5_randi, m6_randi, m7_randi, m8_randi, 
                                      m5_rands, m6_rands, m7_rands, m8_rands), 2)
AIC_chosen_unchosen_less <- round(AIC(l5_randi, l6_randi, l7_randi, l8_randi,
                                      l5_rands, l6_rands, l7_rands, l8_rands), 2)

rbind(AIC_diff_sum_more[which(AIC_diff_sum_more$AIC == min(AIC_diff_sum_more$AIC)), ], 
      AIC_diff_sum_less[which(AIC_diff_sum_less$AIC == min(AIC_diff_sum_less$AIC)), ], 
      AIC_chosen_unchosen_more[which(AIC_chosen_unchosen_more$AIC == min(AIC_chosen_unchosen_more$AIC)), ], 
      AIC_chosen_unchosen_less[which(AIC_chosen_unchosen_less$AIC == min(AIC_chosen_unchosen_less$AIC)), ]) %>% 
  mutate(frame = rep(c("more", "less"), 2), 
         rule = c(rep("ds", 2), rep("ch", 2))) %>% 
  ggplot(aes(x = frame, y = AIC)) + 
  geom_point(aes(shape = rule), size = 5) + 
  scale_x_discrete(limits = c("more", "less"), labels = c(more = "More", less = "Less")) + 
  scale_shape_manual(values = c(5, 8), labels = c(ch = "Chosen-unchosen", ds = "Diff-sum")) + 
  guides(shape = guide_legend(title = "Rule")) + 
  xlab("Frame") + ylab("AIC") + 
  theme_classic(base_size = 22, base_line_size = 1) + 
  theme(axis.text = element_text(face = "bold"), 
        legend.position = c(0.3, 0.85),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.margin = unit(c(10, 10, 10, 10), "mm")) -> p3

## LMMs
data.frame(Coefficient = c(summary(bestFit_more)$coef[1:4, 1], 
                           summary(bestFit_less)$coef[1:4, 1]), 
           SEM = c(summary(bestFit_more)$coef[1:4, 2], 
                   summary(bestFit_less)$coef[1:4, 2])) %>% 
  mutate(Variable = rep(c("Intercept", "Chosen", "Unchosen", "Interaction"), 2), 
         Condition = c(rep("MORE", 4), rep("LESS", 4))) -> df_coef
ggplot(df_coef, aes(x = Variable, y = Coefficient, fill = Condition)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Coefficient - SEM, ymax = Coefficient + SEM), 
                linewidth = 1.1, width = 0.2, position = position_dodge(0.9), 
                color = "#7d7d7d") + 
  scale_x_discrete(limit = c("Chosen", "Unchosen", "Interaction", "Intercept")) + 
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
        plot.margin = unit(c(10, 10, 10, 10), "mm")) -> p4

plot_grid(plot_grid(p1, p2, labels = c("a", "b"), label_size = 32), plot_grid(p3, p4, align = "h", labels = c("c", "d"), label_size = 32), nrow = 2)
ggsave("figure3.png", width = 1440, height = 1200, units = "px", scale = 3.2)