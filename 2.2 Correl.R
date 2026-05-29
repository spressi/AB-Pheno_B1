library(tidyverse)

eeg = read_rds("eeg.rds" %>% paste0(path.rds, .))
behavior.valid = read_rds("behavior.valid.rds" %>% paste0(path.rds, .))

data = behavior.valid %>% summarize(.by = c(subject, paradigm, SOA, congruency),
                             rt = mean(rt)) %>% 
  filter(paradigm == "Dot Probe", SOA == 100) %>% 
  pivot_wider(names_from = congruency, values_from = rt) %>% 
  transmute(subject = subject, RTBias = angry - neutral) %>% 
  left_join(eeg %>% select(subject, N2pc))

with(data, cor.test(RTBias, N2pc, alternative = "greater")) %>% apa::cor_apa(r_ci = T)


data.explore = behavior.valid %>% summarize(.by = c(subject, paradigm, SOA, congruency),
                                    rt = mean(rt)) %>% 
  #filter(paradigm == "Dot Probe", SOA == 100) %>% 
  pivot_wider(names_from = congruency, values_from = rt) %>% 
  mutate(RTBias = angry - neutral) %>% 
  left_join(eeg %>% select(subject, N2pc))

data.explore %>% summarize(.by = c(paradigm, SOA),
                           cor.test = cor.test(RTBias, N2pc) %>% apa::cor_apa(r_ci=T, print=F),
                           cor = cor(RTBias, N2pc, use = "pairwise.complete.obs")) %>% arrange(desc(abs(cor))) %>% select(-cor)
