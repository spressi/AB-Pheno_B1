library(tidyverse)
#source("0 General.R")

files.behav = list.files(path.behav, pattern = ".log", full.names = T) %>% 
  Filter(\(x) x %>% grepl("_0", .) == F, .) %>% #get rid of training logs
  Filter(\(x) x %>% grepl("Pre.log", .) == F, .) #get rid of restarts

behavior = files.behav %>% lapply(\(x) x %>% read_delim(delim="\t", skip=3, show_col_types=F, name_repair="minimal")) %>% bind_rows() %>% 
  rename(event = `Event Type`) %>% select(Subject, event, Code, Time) %>% 
  separate(Subject, c("subject", "block")) %>% 
  filter(event == "Response" & Code != "1" |
           Code == "targets" |
           Code %>% grepl("distractors", .))
behavior = behavior %>% filter(Code %>% grepl("distractors", .)) %>% 
  mutate(.by = subject, trial = 1:n()) %>% 
  full_join(behavior, .) %>% 
  group_by(subject, block) %>% fill(trial, .direction = "down") %>% ungroup()
  
behavior %>% summarize(.by = subject, trial = max(trial)) %>% filter(trial!=576)

behavior %>% filter(event == "Response") %>% pull(Code) %>% as.integer() %>% summary()
#TODO remap response keys to correct numbers (depending on paradigm)
