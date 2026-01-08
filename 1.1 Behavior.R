library(tidyverse)
#source("0 General.R")

files.behav = list.files(path.behav, pattern = ".log", full.names = T) %>% 
  Filter(\(x) x %>% grepl("_0", .) == F, .) %>% #get rid of training logs
  Filter(\(x) x %>% grepl("Pre.log", .) == F, .) #get rid of restarts

behavior.overview = tibble(name = files.behav %>% pathToCode()) %>% separate(name, sep="_|-", into=c("subject", "block", "experiment", "noET"))
behavior.overview %>% count(subject) %>% filter(n != 2) #a07 no 2nd block
behavior.overview %>% filter(noET %>% is.na() == F) #a23: no eye-tracking 2nd block; b18, b22: no eye-tracking

behavior = files.behav %>% lapply(\(x) x %>% read_delim(delim="\t", skip=3, show_col_types=F, name_repair="minimal")) %>% bind_rows() %>% 
  rename(event = `Event Type`) %>% rename_with(tolower) %>% select(subject, event, code, time) %>% 
  separate(subject, c("subject", "block")) %>% 
  filter(event == "Response" & code != "1" | #get rid of 1 and 2 responses (space bar & experimenter key)
           code == "targets" |
           code %>% grepl("distractors", .)) %>% 
  mutate(paradigm = if_else(subject %>% str_starts("a"), "Dot Probe", "Dual Probe") %>% as_factor())
behavior = behavior %>% 
  filter(code %>% grepl("distractors", .)) %>% mutate(.by = subject, trial = 1:n()) %>% #use distractors to count trials
  full_join(behavior, .) %>% #join result back to original data.frame (will have NAs)
  group_by(subject, block) %>% fill(trial, .direction = "down") %>% ungroup() %>% #fill NAs
  filter(trial %>% is.na() == F) %>% #get rid of responses before first trial (i.e., responses to instructions)
  filter(event != lag(event) | event != "Response") #only keep first response in each trial (i.e., detect where multiple responses succeed each other without a Picture event)
#behavior %>% filter(event == lag(event), event == "Response") #check multiple responses
#behavior %>% count(subject, trial) %>% filter(n < 3) %>% left_join(behavior) #check missing responses

behavior %>% summarize(.by = subject, trial = max(trial)) %>% filter(trial!=576) #a07 no 2nd block (cf. behavior.overview)


# Remap Responses ---------------------------------------------------------
#behavior %>% filter(event == "Response") %>% summarize(.by = paradigm, summary = code %>% as.integer() %>% summary() %>% bind_rows())
#behavior %>% filter(event == "Response") %>% mutate(code = code %>% str_pad(2, "left", "0")) %>% ggplot(aes(x = code, fill = paradigm)) + geom_histogram(stat = "count", color = "black")

behavior = behavior %>% 
  mutate(event = case_when(code == "targets" ~ "target", #change "Picture" to "target"
                           event == "Picture" ~ "distractors", #change remaining "Picture"s to "distractors"
                           T ~ event %>% tolower())) %>% #keep response (but lower case)
  pivot_wider(names_from = event, values_from = c(code, time)) %>% 
  select(-code_target) %>% rename(code = code_distractors, response = code_response) %>% 
  
  #conditions
  mutate(code = code %>% gsub("distractors/", "", .)) %>% 
  separate(code, into = c("distractor_left", "target_left", NA, "distractor_right", "target_right"), sep = " ") %>% 
  
         #distractors
  mutate(angry = case_when(distractor_left %>% grepl("category1", .) ~ "left",
                           distractor_right %>% grepl("category1", .) ~ "right",
                           T ~ NA) %>% as_factor(),
         across(starts_with("distractor_"), \(x) x %>% gsub("category\\d+/", "", .) %>% gsub(".jpg", "", .)),
         
         #targets: Dot Probe
         across(starts_with("target_"), \(x) if_else(x == "NA", NA, x)),
         target_dotprobe = case_when(target_right %>% is.na() ~ "left",
                                     target_left %>% is.na() ~ "right",
                                     T ~ NA) %>% as_factor(),
         congruency_dotprobe = if_else(angry == target_dotprobe, "angry", "neutral") %>% as_factor(),
         targetKind_dotprobe = case_when(target_right %>% is.na() ~ target_left,
                                         target_left %>% is.na() ~ target_right,
                                         T ~ NA) %>% as_factor(),
         
         #targets: Dual Probe
         target_left = case_when(paradigm == "Dot Probe" ~ NA,
                                 T ~ target_left %>% as.integer()),
         target_right = case_when(paradigm == "Dot Probe" ~ NA,
                                  T ~ target_right %>% as.integer())
  ) %>% 
  
  #responses
  mutate(response = response %>% as.integer(),
         response = case_when(paradigm == "Dot Probe" ~ response - 1,
                              paradigm == "Dual Probe" ~ response - 3),
         response = if_else(response >= 5, response + 1, response), #5 key does not exist in Dual Probe => add 1
         keyAssign = case_when(paradigm == "Dual Probe" ~ NA, #no key assignment for dual probe
                               #only for dot probe:
                               subject %>% gsub("\\D", "", .) %>% as.integer() %>% {. %% 2} == 1 ~ 1, #subject number odd: 1
                               T ~ 2), #subject number even: 2
         response_dotprobe = case_when(paradigm == "Dual Probe" ~ NA,
                                 keyAssign == 1 ~ if_else(response==1, ":", "..") == targetKind_dotprobe,
                                 T ~ if_else(response==1, "..", ":") == targetKind_dotprobe),
         response_dual = case_when(paradigm == "Dot Probe" ~ NA,
                                   response %>% is.na() ~ NA,
                                   response == target_left ~ "left",
                                   response == target_right ~ "right",
                                   T ~ "incorrect") %>% as_factor(), #wrong response (target not present)
         congruency_dual = if_else(response_dual %>% as.character() == angry, "angry", "neutral") %>% as_factor() #note: in dual probe, this is the "response-congruency" (not a design parameter of the trial)
  ) %>%
  mutate(congruency = if_else(paradigm == "Dot Probe", congruency_dotprobe, congruency_dual)) %>% select(-congruency_dotprobe, -congruency_dual) %>% 
  
  #times
  mutate(rt = (time_response - time_target)/10, #response time in ms
         expositionCheck = (time_target - time_distractors)/10, #exposition time in ms (can be NA if premature response before target onset was given)
         SOA = if_else(expositionCheck < 500, 100, 500) %>% as_factor()) %>% #SOA = factorized exposition time
  select(subject, paradigm, block, trial, SOA, congruency, angry, response, rt, starts_with("distractor"), starts_with("target"), contains("dotprobe"), contains("dual"), expositionCheck, starts_with("time_"))


# Quality Checks ----------------------------------------------------------

#NAs
behavior %>% summarize(.by = c(subject, paradigm), NA_n = sum(response %>% is.na()), `NA_%` = mean(response %>% is.na())) %>% summarize(.by = paradigm, across(starts_with("NA_"), list(m = mean, sd = sd)))
#behavior %>% summarize(.by = c(subject, paradigm), NA_n = sum(response %>% is.na()), `NA_%` = mean(response %>% is.na())) %>% filter(NA_n > 0)
#behavior %>% filter(response %>% is.na())

#premature responses (i.e., before target onset)
#behavior %>% summarize(.by = c(subject, paradigm), pre_n = sum(expositionCheck %>% is.na()), `pre_%` = mean(expositionCheck %>% is.na())) %>% summarize(.by = paradigm, across(starts_with("pre_"), list(m = mean, sd = sd)))
behavior %>% summarize(.by = c(subject, paradigm), pre_n = sum(expositionCheck %>% is.na()), `pre_%` = mean(expositionCheck %>% is.na())) %>% filter(pre_n > 0)
#behavior %>% filter(expositionCheck %>% is.na()) %>% select(subject:trial, rt, contains("time")) %>% mutate(rt2 = (time_response - time_distractors)/10)

#wrong responses
behavior %>% filter(paradigm == "Dot Probe") %>% summarize(.by = c(subject, paradigm), wrong_n = sum(response_dotprobe == FALSE, na.rm=T), `wrong_%` = mean(response_dotprobe == FALSE, na.rm=T)) %>% summarize(.by = paradigm, across(starts_with("wrong_"), list(m = mean, sd = sd))) %>% 
  bind_rows(behavior %>% filter(paradigm == "Dual Probe") %>% summarize(.by = c(subject, paradigm), wrong_n = sum(response_dual == "incorrect", na.rm=T), `wrong_%` = mean(response_dual == "incorrect", na.rm=T)) %>% summarize(.by = paradigm, across(starts_with("wrong_"), list(m = mean, sd = sd))))
#behavior %>% filter(paradigm == "Dual Probe") %>% summarize(.by = c(subject, paradigm), wrong_n = sum(response_dual == "incorrect", na.rm=T), `wrong_%` = mean(response_dual == "incorrect", na.rm=T)) %>% filter(wrong_n > 0)

#=> unusable trials
behavior %>% summarize(.by = c(subject, paradigm), 
                       missing_n = sum(response %>% is.na() | expositionCheck %>% is.na() | (paradigm == "Dot Probe" & response_dotprobe == F) | (paradigm == "Dual Probe" & response_dual == "incorrect")),
                       `missing_%` = mean(response %>% is.na() | expositionCheck %>% is.na() | (paradigm == "Dot Probe" & response_dotprobe == F) | (paradigm == "Dual Probe" & response_dual == "incorrect"))) %>% 
  summarize(.by = paradigm, across(starts_with("missing_"), list(m = mean, sd = sd)))


#design
behavior %>% summarize(.by = paradigm,
                       SOA = mean(SOA == "100", na.rm=T),
                       congruency = mean(congruency == "angry", na.rm=T),
                       angry = mean(angry == "left", na.rm=T),
                       #only for dot probe:
                       target_dp = mean(target_dotprobe == "left", na.rm=T),
                       targetKind_dp = mean(targetKind_dotprobe == ":", na.rm=T))

#stimulus presentation
behavior %>% count(expositionCheck) #rarely, 1 frame was skipped (+ minimally unstable frame rate)
#behavior %>% pull(expositionCheck) %>% hist()
#behavior %>% pull(expositionCheck) %>% unique() %>% sort()

#response coding
#behavior %>% ggplot(aes(x = response, fill = paradigm)) + geom_histogram(stat = "count", color = "black") + facet_wrap(~paradigm, scales = "free")

#Dual Probe responses
behavior %>% summarize(response_dual_left = mean(response_dual == "left", na.rm=T),
                       response_dual_right = mean(response_dual == "right", na.rm=T),
                       response_dual_wrong = mean(response_dual == "incorrect", na.rm=T)) #%>% rowSums()


# Analysis ----------------------------------------------------------------
behavior.valid = behavior %>% filter(response %>% is.na() == F,
                                     expositionCheck %>% is.na() == F,
                                     response_dotprobe != F | response_dual != "incorrect") #kick out incorrect responses
