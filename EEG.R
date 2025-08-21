library(tidyverse)
#source("0 General.R")

writeCorrectedMarkers = F #rewrite marker file for subjects with inverted markers (low voltage = signal instead of high)

# Markers -----------------------------------------------------------------
files.eeg.markers = list.files(path.eeg.raw, pattern = ".mrk", full.names = T) %>% 
  Filter(\(x) x %>% grepl("_2", .) == F, .) #get rid of second file
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

eeg.markers %>% count(subject) %>% filter(n != 1152) %>% arrange(n)

#correct for inverted markers
invertedMarkers = eeg.markers %>% count(subject) %>% filter(n > 1152) %>% pull(subject) #some subjects had inverted markers (low voltage = signal instead of high). Every reset of the pins thus resulted in an additional marker
#eeg.markers %>% filter(subject %in% invertedMarkers) #255 is marker stop => recode as 255 - value
eeg.markers = eeg.markers %>% 
  mutate(value = if_else(subject %in% invertedMarkers, 255 - value, value)) %>% 
  filter(value != 0) %>% #get rid of stop markers
  #check assumption: first marker is always "Mk2" (last is "Mk2304")
  #summarize(.by = subject, marker.min = first(marker), marker.max = last(marker)) %>% filter(marker.min != "Mk2=Stimulus" | marker.max != "Mk2304=Stimulus", subject %in% invertedMarkers)
  mutate(.by = subject, trial = 1:n()) %>% 
  #check assumption: only subjects with inverted markers are affected by wrong marker order
  #mutate(markerCheck = paste0("Mk", trial+1, "=Stimulus")) %>% filter(markerCheck != marker) %>% pull(subject) %>% unique() %>% {. == invertedMarkers} %>% all()
  mutate(marker = paste0("Mk", trial+1, "=Stimulus"), #this also overwrites correct markers but result has been asserted to be the same
         output = paste(marker, paste0("S", value), sample, size, channel, sep = ","))
#eeg.markers %>% filter(subject %in% invertedMarkers)
eeg.markers %>% count(subject) %>% filter(n != 1152) %>% arrange(n)
  
if (writeCorrectedMarkers) {
  for (s in invertedMarkers) {
    #s = sample(invertedMarkers, 1) #for testing
    filename = files.eeg.markers %>% Filter(\(x) x %>% grepl(s, .), .)
    filename.copy = filename %>% gsub(".vmrk", "_original.vmrk", ., fixed = T)
    if (file.exists(filename.copy)==F) #careful! there is no option that prevents file.rename from overwriting an existing file => check yourself
      file.rename(filename, filename.copy) #this retains "last modified" as the original file creation date
    
    file = readLines(filename.copy)
    #file[12] #assert that line 12 is the last line that should remain unmodified
    file = c(file[1:12],
             eeg.markers %>% filter(subject == s) %>% pull(output))
    writeLines(file, filename)
  }
}

# Impedances --------------------------------------------------------------
files.eeg.headers = list.files(path.eeg.raw, pattern = ".vhdr", full.names = T)
eeg.impedances.list = list()
for (file in files.eeg.headers) {
  #file = files.eeg.headers %>% sample(1) #for testing
  
  skip = 112 #start with skipping 112 lines (b06 for whatever reason)
  repeat { #do-while loop
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
