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

se = function(x, na.rm = FALSE) {
  sd(x, na.rm) / sqrt(if(!na.rm) length(x) else sum(!is.na(x)))
}

correlation_out = function(coroutput) {
  names = coroutput$data.name %>% strsplit(" and ") %>% unlist()
  cat(paste0("r(", names[1], ", ", names[2], "): ", coroutput %>% apa::cor_apa(print=F)), "\n")
}

F_apa = function(x) {
  cat("F(", paste(x$parameter, collapse=", "), ") = ", 
      round(x$statistic, 2), 
      ", p ",
      ifelse(x$p.value < .001, "< .001",
             paste0("= ", round(x$p.value, 2))),
      "\n", sep="")
}

ez.ci = function(ez, conf.level = .95, sph.corr=T) {
  for (effect in ez$ANOVA$Effect) {
    index.sphericity = which(ez$`Sphericity Corrections`$Effect == effect)
    GGe = ifelse(sph.corr==F || length(index.sphericity)==0, 1, ez$`Sphericity Corrections`$GGe[index.sphericity])
    
    index.effect = which(ez$ANOVA$Effect == effect)
    ez$ANOVA %>% with(apaTables::get.ci.partial.eta.squared(F[index.effect], DFn[index.effect]*GGe, DFd[index.effect]*GGe, conf.level = conf.level)) %>% 
      sapply(round, digits=2) %>% 
      paste0(collapse=", ") %>% paste0(effect, ": ", round(conf.level*100), "% CI [", ., "]\n") %>% cat()
  }
}

lmer.ci = function(lmer, conf.level = .95, twotailed=T) {
  values = lmer %>% summary() %>% .$coefficients %>% .[, c("Estimate", "df")] %>% data.frame()
  effects = values %>% rownames()
  
  for (i in seq(effects)) {
    psych::r.con(values$Estimate[i], values$df[i], p = conf.level, twotailed=twotailed) %>% 
      round(digits=2) %>% 
      paste0(collapse=", ") %>% paste0(effects[i], ": ", round(conf.level*100), "% CI [", ., "]\n") %>% cat()
  }
}

dodge.width = .6
dodge = position_dodge(width=dodge.width)

#ggplot general theme
theme_set(myGgTheme <- theme_bw() + theme(
  #aspect.ratio = 1,
  plot.title = element_text(hjust = 0.5),
  panel.background = element_rect(fill="white", color="white"),
  legend.background = element_rect(fill="white", color="grey"),
  legend.key=element_rect(fill='white'),
  axis.text = element_text(color="black"),
  axis.ticks.x = element_line(color="black"),
  axis.line.x = element_line(color="black"),
  axis.line.y = element_line(color="black"),
  legend.text = element_text(size=14, color="black"),
  legend.title = element_text(size=14, color="black"),
  strip.text.x = element_text(size=12, color="black"),
  axis.text.x = element_text(size=16, color="black"),
  axis.text.y = element_text(size=16, color="black"),
  axis.title = element_text(size=16, color="black"))
)
