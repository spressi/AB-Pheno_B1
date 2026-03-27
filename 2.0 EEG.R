library(tidyverse)

eeg = 
  "N2pc_BrainVision.csv" %>%
  #"N2pc_BrainVision_late.csv" %>% 
  #"N2pc_BrainVision_cherrypicked.csv" %>% 
  paste0(path.eeg.raw, "export/", .) %>% read_csv2() %>% 
  rename(subject = File) %>% 
  rename_all(\(x) x %>% str_replace("contraP78-ipsiP78-Diff_P78", "N2pc")) %>% 
  mutate(paradigm = case_when(subject %>% str_starts("a") ~ "Dot Probe",
                              subject %>% str_starts("b") ~ "Dual Probe",
                              T ~ NA) %>% as_factor(),
         N2pc = (N2pc_AN+N2pc_NA)/2,
         N2pc_Odd = (N2pc_AN_Odd+N2pc_NA_Odd)/2,
         N2pc_Even = (N2pc_AN_Even+N2pc_NA_Even)/2) %>% 
  relocate(paradigm, N2pc, N2pc_Odd, N2pc_Even, .before = N2pc_AN)

with(eeg, t.test(N2pc, alternative="less")) %>% apa::t_apa(es_ci=T)

eeg.long = eeg %>% select(subject, paradigm, AN = N2pc_AN, `NA` = N2pc_NA) %>% 
  pivot_longer(c(AN, `NA`), names_to = "side", values_to = "N2pc")
  
eeg.long %>% ez::ezANOVA(dv = N2pc,
                         wid = subject,
                         within = side,
                         between = paradigm,
                         type = 2, detailed = T) %>% apa::anova_apa()

eeg %>% summarize(N2pc.m = mean(N2pc, na.rm=T),
                  N2pc.se = se(N2pc, na.rm=T))

eeg.long %>% summarize(.by = side,
                       N2pc.m = mean(N2pc, na.rm=T),
                       N2pc.se = se(N2pc, na.rm=T))

eeg.long %>% summarize(.by = c(paradigm, side),
                       N2pc.se = se(N2pc, na.rm=T),
                       N2pc = mean(N2pc, na.rm=T)) %>% 
  ggplot((aes(y = N2pc, x = side, color = paradigm))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = N2pc - N2pc.se, ymax = N2pc + N2pc.se)) +
  geom_point() + 
  myGgTheme

# Reliabilities -----------------------------------------------------------
#TODO spearman brown correction
with(eeg, cor.test(N2pc_Even, N2pc_Odd)) %>% apa::cor_apa()
with(eeg, cor.test(N2pc_AN_Even, N2pc_AN_Odd)) %>% apa::cor_apa()
with(eeg, cor.test(N2pc_NA_Even, N2pc_NA_Odd)) %>% apa::cor_apa()
