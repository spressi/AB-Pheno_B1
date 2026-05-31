library(tidyverse)

eegLabImport = function(filePath) {
  filePath %>% read_csv(show_col_types = F) %>% 
    select(-counter) %>% #running subject number - not needed after subject code reconstruction
    mutate(subjectNum = subjectNum %>% str_pad(width = 2, pad = "0"), #add leading 0s
           paradigm = paradigm %>% recode_values(1 ~ "a", 2 ~ "b"),
           subject = paste0(paradigm, subjectNum),
           paradigm = paradigm %>% recode_values("a" ~ "Dot Probe", "b" ~ "Dual Probe") %>% as_factor(),
           angry = angry %>% recode_values(1 ~ "NA", 2 ~ "AN") %>% as_factor(),
           oddEven = if_else(trial %% 2 == 0, "Even", "Odd")) %>% 
    relocate(subject) %>% select(-subjectNum) %>% 
    mutate(N2pc = if_else(angry == "left", P8 - P7, P7 - P8)) %>% relocate(N2pc, .after = angry)
}

eegLabSummary = function(eegLabImport) {
  eegLabImport %>% summarize(.by = c(subject, paradigm), N2pc = mean(N2pc)) %>% 
    full_join(eegLabImport %>% summarize(.by = c(subject, paradigm, oddEven), N2pc = mean(N2pc)) %>% 
                pivot_wider(names_prefix = "N2pc_", names_from = oddEven, values_from = N2pc), 
              by = join_by(subject, paradigm)) %>% 
    full_join(eegLabImport %>% summarize(.by = c(subject, paradigm, angry), N2pc = mean(N2pc)) %>% 
                pivot_wider(names_prefix = "N2pc_", names_from = angry, values_from = N2pc), 
              by = join_by(subject, paradigm)) %>% 
    full_join(eegLabImport %>% summarize(.by = c(subject, paradigm, angry, oddEven), N2pc = mean(N2pc)) %>% 
                pivot_wider(names_prefix = "N2pc_", names_from = c(angry, oddEven), values_from = N2pc), 
              by = join_by(subject, paradigm)) %>% 
    
    #comparability to old processing (Reutter et al., 2017; 2019): 
    #weigh AN and NA evenly (even if they are comprised of unequal number of trials)
    #(barely has an effect, cf. number of trials)
    mutate(N2pc = (N2pc_AN+N2pc_NA)/2,
           N2pc_Odd = (N2pc_AN_Odd+N2pc_NA_Odd)/2,
           N2pc_Even = (N2pc_AN_Even+N2pc_NA_Even)/2)
}

eeg.trial = "EEGlab_180-300_prereg.csv" %>% paste0(path.eeg, .) %>% eegLabImport()

#number of trials
eeg.trial %>% count(subject, angry) %>% filter(n < trials.min) #no subject needs to be excluded according to preregistration
eeg.trial %>% count(subject) %>% summarize(valid.m = mean(n), valid.sd = sd(n), valid.min = min(n), valid.max = max(n)) #%>% mutate(across(everything(), \(x) x / trials.N))
eeg.trial %>% count(subject) %>% mutate(p = n / trials.N) %>% arrange(n)
eeg.trial %>% count(subject) %>% ggplot(aes(x = n)) + geom_histogram(color = "black") + geom_vline(xintercept = trials.min, color = "red", linetype = "dashed", linewidth = 2) + xlab("Valid EEG Trials") + myGgTheme


#subject-level aggregates
eeg = eeg.trial %>% eegLabSummary()


# Old: BrainVision Quick & Dirty
# eeg = 
#   "N2pc_BrainVision.csv" %>%
#   #"N2pc_BrainVision_late.csv" %>% 
#   #"N2pc_BrainVision_cherrypicked.csv" %>% 
#   paste0(path.eeg.raw, "export/", .) %>% read_csv2() %>% 
#   rename(subject = File) %>% 
#   rename_all(\(x) x %>% str_replace("contraP78-ipsiP78-Diff_P78", "N2pc")) %>% 
#   mutate(paradigm = case_when(subject %>% str_starts("a") ~ "Dot Probe",
#                               subject %>% str_starts("b") ~ "Dual Probe",
#                               T ~ NA) %>% as_factor(),
#          N2pc = (N2pc_AN+N2pc_NA)/2,
#          N2pc_Odd = (N2pc_AN_Odd+N2pc_NA_Odd)/2,
#          N2pc_Even = (N2pc_AN_Even+N2pc_NA_Even)/2) %>% 
#   relocate(paradigm, N2pc, N2pc_Odd, N2pc_Even, .before = N2pc_AN)

symdiff(eeg %>% pull(subject), behavior %>% pull(subject)) #TODO: check missing b12

eeg %>% write_rds("eeg.rds" %>% paste0(path.rds, .))

#TODO move upper part to separate preprocessing script (1.2.2 EEG?)

#eeg = read_rds("eeg.rds" %>% paste0(path.rds, .))

eeg = "EEGlab_150-220.csv" %>% paste0(path.eeg, .) %>% eegLabImport() %>% eegLabSummary() #cherry picked by grand average

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
  ggplot(aes(y = N2pc, x = side, color = paradigm)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = N2pc - N2pc.se, ymax = N2pc + N2pc.se), linewidth=2, position = dodge) +
  geom_point(size = 6, position = dodge) + 
  myGgTheme

# Reliabilities -----------------------------------------------------------
eeg %>% summarize(rel = cor.test(N2pc_Even, N2pc_Odd) %>% apa::cor_apa(r_ci=T, print=F),
                  rel_sb = cor(N2pc_Even, N2pc_Odd) %>% spearmanBrown())
eeg %>% summarize(rel = cor.test(N2pc_AN_Even, N2pc_AN_Odd) %>% apa::cor_apa(r_ci=T, print=F),
                  rel_sb = cor(N2pc_AN_Even, N2pc_AN_Odd) %>% spearmanBrown())
eeg %>% summarize(rel = cor.test(N2pc_NA_Even, N2pc_NA_Odd) %>% apa::cor_apa(r_ci=T, print=F),
                  rel_sb = cor(N2pc_NA_Even, N2pc_NA_Odd) %>% spearmanBrown())

eeg %>% summarize(.by = paradigm,
                  rel = cor.test(N2pc_Even, N2pc_Odd) %>% apa::cor_apa(r_ci=T, print=F),
                  rel_sb = cor(N2pc_Even, N2pc_Odd) %>% spearmanBrown())
