library(tidyverse)


# Paths -------------------------------------------------------------------
path <- "C:/Data/AB_B1/Data/" #@work
#path = path %>% gsub("C:/Data", "D:/Arbeit", .) #@home

path.behav = "log/" %>% paste0(path, .)
path.seq = paste0(path, "../sequences/")

path.eye.raw = "Eye/" %>% paste0(path, .)
path.eye = "Summary/" %>% paste0(path.eye.raw, .)

path.eeg.raw = "EEG/" %>% paste0(path, .)
path.eeg = "" #TODO preprocessed files


# Functions ---------------------------------------------------------------
pathToCode = function(path, path.sep="/", file.ext="\\.") {
  first = path %>% gregexpr(path.sep, .) %>% lapply(max) %>% unlist() %>% {. + 1} %>% 
    pmax(1, .) #if first not found, set it to start of string
  last = path %>% gregexpr(file.ext, .) %>% lapply(max) %>% unlist() %>% {. - 1}
  last = ifelse(last < 1, sapply(path, str_length), last) #if last cannot be found, set it to end of string
  return(path %>% substring(first, last))
}
