library(tidyverse)
#source("0 General.R")


# Impedances --------------------------------------------------------------
files.eeg.headers = list.files(path.eeg.raw, pattern = ".vhdr", full.names = T)
eeg.impedances.list = list()
for (file in files.eeg.headers) {
  #file = files.eeg.headers %>% sample(1) #for testing
  
  skip = 112 #start with skipping 112 lines (b06 for whatever reason)
  repeat {
    checkFile = file %>% 
      read_table(skip = skip, col_names = c("electrode", "impedance"), show_col_types = F, na = "???")
    if (checkFile %>% pull(electrode) %>% grepl("Fp1", .) %>% any()) { #check if Fp1 is contained (first electrode)
      break
    } else {
      skip = skip - 1
    }
  }
  
  eeg.impedances.list[[pathToCode(file)]] = checkFile
}
#tidy up
eeg.impedances = eeg.impedances.list %>% bind_rows(.id = "subject") %>% 
  mutate(electrode = electrode %>% gsub(":", "", .),
         time = if_else(subject %>% grepl("_2", .), "after", "before") %>% as_factor(),
         subject_session = subject,
         subject = subject %>% gsub("_2", "", .))


#sanity checks
eeg.impedances %>% count(electrode) %>% count(n) %>% rename(subjects = n, electrodes = nn)
eeg.impedances %>% filter(impedance %>% is.na()) %>% pull(subject_session) %>% unique() #TODO check NAs

##impedance change before/after
# with(eeg.impedances %>% summarize(.by = c(subject, time),
#                              impedance = mean(impedance, na.rm=T)) %>% 
#   pivot_wider(id_cols = subject, names_from = time, values_from = impedance),
#   t.test(before, after, paired=T)
# ) %>% apa::t_apa(es_ci=T)
eeg.impedances %>% summarize(.by = c(subject, time),
                             impedance = mean(impedance, na.rm=T))
#just saving another file does not work to measure impedances after experiment :(
# => need to click on impedance measurement again; otherwise previous values will just get carried over



# Markers -----------------------------------------------------------------
files.eeg.markers = list.files(path.eeg.raw, pattern = ".mrk", full.names = T) %>% 
  Filter(\(x) x %>% grepl("_2", .) == F, .) #get rid of everything containing "_2"
eeg.markers.list = list()
for (file in files.eeg.markers) {
  #file = files.eeg.markers %>% sample(1) #for testing
  
  eeg.markers.list[[pathToCode(file)]] = file %>% read.csv(skip=11, header=F, col.names = c("marker", "value", "sample", "size", "channel"))
}
#tidy up
eeg.markers = eeg.markers.list %>% bind_rows(.id = "subject") %>% tibble() %>% 
  filter(marker %>% grepl("Stimulus", .)) %>% 
  mutate(value = value %>% gsub("S\\s*", "", .) %>% as.integer(),
         paradigm = if_else(subject %>% grepl("a", .), "Dot Probe", "Dual Probe"))
eeg.markers %>% filter(paradigm == "Dot Probe") %>% count(value) %>% arrange()
eeg.markers %>% filter(paradigm == "Dual Probe") %>% count(value) %>% arrange(n)
