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

behavior.aov.dot %>% summarize(.by = c(SOA, congruency),
                               rt.se = se(rt), rt = mean(rt)) %>% 
  ggplot(aes(x = SOA, y = rt, color = congruency)) +
  geom_errorbar(aes(ymin=rt-rt.se*1.96, ymax=rt+rt.se*1.96), linewidth=2, position = dodge) +
  geom_point(size=6, position = dodge) +
  labs(y = "RT (ms)") +
  myGgTheme

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
            rel_sb = cor(neutral_odd - angry_odd, neutral_even - angry_even) %>% spearmanBrown())


#RTs
behavior.reliability %>% pivot_wider(names_from = split, values_from = rt) %>% 
  summarize(.by = c(paradigm, SOA, congruency),
            reliability = cor.test(odd, even) %>% apa::cor_apa(r_ci=T, print=F),
            rel_sb = cor(odd, even) %>% spearmanBrown())
