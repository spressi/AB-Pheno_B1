library(tidyverse)

#behavior.valid = read_rds("behavior.valid.rds" %>% paste0(path.rds, .))

behavior.aov = behavior.valid %>% summarize(.by = c(subject, paradigm, SOA, congruency),
                                            rt = mean(rt)) %>% 
  mutate(SOA = SOA %>% as_factor())

# full analysis
behavior.aov %>%
  ez::ezANOVA(dv = rt, wid = subject,
              within = .(SOA, congruency),
              between = paradigm,
              type=2, detailed=T) %>% apa::anova_apa(force_sph_corr = T)

behavior.aov %>% summarize(.by = c(paradigm, SOA, congruency),
                           rt.se = se(rt), rt = mean(rt)) %>%
  ggplot(aes(x = SOA, y = rt, color = congruency)) + facet_wrap(~paradigm) +
  geom_errorbar(aes(ymin=rt-rt.se*1.96, ymax=rt+rt.se*1.96), size=2, position = dodge) +
  geom_point(size=6, position = dodge) +
  labs(y = "RT (ms)") +
  myGgTheme

# Dot Probe only
behavior.aov.dot = behavior.aov %>% filter(paradigm == "Dot Probe")

behavior.aov.dot %>% 
  ez::ezANOVA(dv = rt, wid = subject,
              within = .(SOA, congruency),
              type=2, detailed=T) %>% apa::anova_apa(force_sph_corr = T)

plot.rt.dot_probe = behavior.aov.dot %>% summarize(.by = c(SOA, congruency),
                               rt.se = se(rt), rt = mean(rt)) %>% 
  ggplot(aes(x = SOA, y = rt, color = congruency)) +
  geom_errorbar(aes(ymin=rt-rt.se*1.96, ymax=rt+rt.se*1.96), linewidth=2, position = dodge) +
  geom_point(size=6, position = dodge) +
  labs(y = "RT (ms)", x = "SOA (ms)") +
  myGgTheme
ggsave("plots/RT Dot Probe.png", plot=plot.rt.dot_probe, scale=1, device="png", dpi=300, units="px", width = 1920, height = 1080)

plot.rtbias.dot_probe =
  behavior.aov.dot %>% summarize(.by = c(subject, SOA, congruency),
                                 rt = mean(rt)) %>% 
  pivot_wider(names_from = congruency, values_from = rt) %>% 
  mutate(bias = angry - neutral) %>% 
  summarize(.by = c(SOA),
            rtbias = mean(bias),
            rtbias.se = se(bias)) %>% 
  ggplot(aes(x = SOA, y = rtbias, color = SOA)) +
  geom_hline(yintercept = 0) +
  #geom_col(color = "black") +
  geom_point(size=6) +
  geom_errorbar(aes(ymin=rtbias-rtbias.se*1.96, ymax=rtbias+rtbias.se*1.96), linewidth=2, width=.5) +
  labs(y = "RT-Bias (ms)", x = "SOA (ms)") +
  scale_color_manual(values = c("#FFC000", "red"), guide="none") + #relate to color coding of hypervigilence-avoidance slide
  myGgTheme
ggsave("plots/RT-Bias Dot Probe.png", plot=plot.rtbias.dot_probe, scale=1, device="png", dpi=300, units="px", width = 1920, height = 1080)

with(behavior.aov.dot %>% filter(SOA=="100") %>% pivot_wider(names_from = congruency, values_from = rt), 
     t.test(neutral, angry, paired=T)) %>% apa::t_apa(es_ci=T)
with(behavior.aov.dot %>% filter(SOA=="500") %>% pivot_wider(names_from = congruency, values_from = rt), 
     t.test(neutral, angry, paired=T)) %>% apa::t_apa(es_ci=T)

behavior.aov.dot %>% summarize(.by = c(SOA, congruency),
                               rt.se = se(rt), rt = mean(rt))

# Dual Probe only
behavior.aov.dual = behavior.aov %>% filter(paradigm == "Dual Probe")

behavior.aov.dual %>% 
  ez::ezANOVA(dv = rt, wid = subject,
              within = .(SOA, congruency),
              type=2, detailed=T) %>% apa::anova_apa(force_sph_corr = T)

behavior.aov.dual %>% summarize(.by = c(SOA, congruency),
                                rt.se = se(rt), rt = mean(rt)) %>% 
  ggplot(aes(x = SOA, y = rt, color = congruency)) +
  geom_errorbar(aes(ymin=rt-rt.se*1.96, ymax=rt+rt.se*1.96), size=2, position = dodge) +
  geom_point(size=6, position = dodge) +
  labs(y = "RT (ms)") +
  myGgTheme

behavior.aov.dual %>% summarize(.by = c(SOA),
                                rt.se = se(rt), rt = mean(rt))


# Analysis: Dual Responses ------------------------------------------------
behavior.valid %>% filter(paradigm == "Dual Probe") %>% 
  mutate(SOA = SOA %>% as_factor()) %>% 
  summarize(.by = c(subject, SOA),
            response_angry = mean(congruency == "angry")-.5) %>% 
  ez::ezANOVA(dv = response_angry, wid = subject,
              within = .(SOA),
              type=2, detailed=T) %>% apa::anova_apa(force_sph_corr = T)


# Reliability -------------------------------------------------------------
behavior.reliability = behavior.valid %>% 
  select(subject, paradigm, trial, SOA, congruency, rt) %>% 
  mutate(.by = c(subject, paradigm, SOA, congruency), #everything but trial and rt
         trial_within = 1:n()) %>% 
  mutate(split = if_else(trial_within %% 2 == 0, "even", "odd")) %>% 
  summarize(.by = c(subject, paradigm, SOA, congruency, split),
            rt = mean(rt))

#RT-differences
behavior.reliability %>% pivot_wider(names_from = c(congruency, split), values_from = rt) %>% 
  summarize(.by = c(paradigm, SOA),
            reliability = cor.test(neutral_odd - angry_odd, neutral_even - angry_even) %>% apa::cor_apa(r_ci=T, print=F),
            rel_sb = cor(neutral_odd - angry_odd, neutral_even - angry_even) %>% spearmanBrown()) %>% 
  arrange(rel_sb)


#RTs
behavior.reliability %>% pivot_wider(names_from = split, values_from = rt) %>% 
  summarize(.by = c(paradigm, SOA, congruency),
            reliability = cor.test(odd, even) %>% apa::cor_apa(r_ci=T, print=F),
            rel_sb = cor(odd, even) %>% spearmanBrown()) %>% 
  arrange(rel_sb)
