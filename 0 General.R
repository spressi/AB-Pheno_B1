library(tidyverse)
library(DescTools) #for Winsorize

exclusions = c()

trials.N = 576
trials.N.block = trials.N/2
hz.eeg = 500
hz.et = 250

trials.min = 20 #minimum amount of trials per subject

prestim = 1000
poststim = c(100, 500) #SOA
targetTime = 1500

expoID = "*" #all messages are beginning of trial

distance = 650 #screen center to subject eyes in mm
screen.width.mm = 520 #screen width in mm
screen.width.px = 1920 #screen width in pixels
screen.height.px = 1200 #screen height in pixels (only for plot)
pixsize = screen.width.mm / screen.width.px  # Pixel size in mm (for scan path length)
distance = 650 #screen center to subject eyes in mm


# Paths -------------------------------------------------------------------
path = "C:/Data/AB_B1/Data/" #@work
#path = path %>% gsub("C:/Data", "D:/Arbeit", .) #@home

path.que = "questionnaires/data_ab-pheno_b1.csv" %>% paste0(path, .)
path.behav = "log/" %>% paste0(path, .)
path.seq = "../sequences/" %>% paste0(path, .)

path.eye.raw = "Eye/" %>% paste0(path, .)
path.eye = "Summary/" %>% paste0(path.eye.raw, .)

path.eeg.raw = "EEG/" %>% paste0(path, .)
path.eeg = "export/EEGlab.csv" %>% paste0(path.eeg.raw, .)

path.rds = "" #project root directory


# Functions ---------------------------------------------------------------
#preprocessing
pathToCode = function(path, path.sep="/", file.ext="\\.") {
  first = path %>% gregexpr(path.sep, .) %>% lapply(max) %>% unlist() %>% {. + 1} %>% 
    pmax(1, .) #if first not found, set it to start of string
  last = path %>% gregexpr(file.ext, .) %>% lapply(max) %>% unlist() %>% {. - 1}
  last = ifelse(last < 1, sapply(path, str_length), last) #if last cannot be found, set it to end of string
  return(path %>% substring(first, last))
}

Winsorize.z = function(x, z = c(-2, 2), ...) {
  low = mean(x, na.rm=T) + sd(x, na.rm=T) * min(z)
  high = mean(x, na.rm=T) + sd(x, na.rm=T) * max(z)
  DescTools::Winsorize(x, val = c(low, high), ...)
}

#descriptives
checkContent = function(df, col, p.denominator=NA, print=T) {
  #symbol handling
  if (suppressWarnings(is.na(rlang::enexpr(p.denominator)) == F) && #p.denominator=NA
      exists(rlang::enexpr(p.denominator)) == F && #not a variable in global environment
      rlang::enexpr(p.denominator) %>% rlang::is_symbol()) { #column name passed without quotation
    p.denominator = rlang::ensym(p.denominator) %>% as.character() #cast column name to character for further evaluation
  }
  
  #type handling
  if (p.denominator %>% is.na() == F) { #not NA (for sum(n))
    if (p.denominator %>% is.numeric() == F) { #not numeric => must be a column name
      if (p.denominator %>% match(df %>% colnames()) %>% is.na()) {
        stop(paste(p.denominator, ": Column not found in data frame"))
      } else {
        p.denominator = df %>% pull(!!p.denominator) %>% unique() %>% length()
      }
    }
  }
  
  #if (p.denominator %>% is.na() == F && p.denominator %>% is.numeric() == F) warning("p.denominator not numeric. Using sum(n).")
  result = df %>% count(!!rlang::ensym(col), .drop=F) %>% 
    arrange(desc(n)) %>% 
    mutate(p = n / if_else(p.denominator %>% is.numeric(), p.denominator, sum(n)))
  if (print) {
    result %>% print(n = nrow(.))
    return(invisible(result))
  }
  return(result)
}

descriptives.list = list(m = mean, sd = sd, min = min, max = max)

se = function(x, na.rm = FALSE) {
  sd(x, na.rm) / sqrt(if(!na.rm) length(x) else sum(!is.na(x)))
}

#statistical analysis
correlation_out = function(coroutput) {
  names = coroutput$data.name %>% strsplit(" and ") %>% unlist()
  cat(paste0("r(", names[1], ", ", names[2], "): ", coroutput %>% apa::cor_apa(print=F)), "\n")
}

spearmanBrown = function(r, n=2) {
  n * r / (1 + (n-1) * r)
}

F_apa = function(x) {
  cat("F(", paste(x$parameter, collapse=", "), ") = ", 
      round(x$statistic, 2), 
      ", p ",
      ifelse(x$p.value < .001, "< .001",
             paste0("= ", round(x$p.value, 2))),
      "\n", sep="")
}

ez.ci = function(ez, conf.level = .9, sph.corr=T) { # 90% CIs are recommended for partial eta² (https://daniellakens.blogspot.com/2014/06/calculating-confidence-intervals-for.html#:~:text=Why%20should%20you%20report%2090%25%20CI%20for%20eta%2Dsquared%3F)
  for (effect in ez$ANOVA$Effect) {
    index.sphericity = which(ez$`Sphericity Corrections`$Effect == effect)
    GGe = ifelse(sph.corr==F || length(index.sphericity)==0, 1, ez$`Sphericity Corrections`$GGe[index.sphericity])
    
    index.effect = which(ez$ANOVA$Effect == effect)
    ez$ANOVA %>% with(apaTables::get.ci.partial.eta.squared(F[index.effect], DFn[index.effect]*GGe, DFd[index.effect]*GGe, conf.level = conf.level)) %>% 
      sapply(round, digits=2) %>% 
      paste0(collapse=", ") %>% paste0(effect, ": ", round(conf.level*100), "% CI [", ., "]\n") %>% cat()
  }
}

peta.ci <- #https://spressi.github.io/precision_workshop/2.2_CIs.html#/ci-around-partial-eta%C2%B2-1
  function(anova_table, conf.level = .9) { # 90% CIs are recommended for partial eta² (https://daniellakens.blogspot.com/2014/06/calculating-confidence-intervals-for.html#:~:text=Why%20should%20you%20report%2090%25%20CI%20for%20eta%2Dsquared%3F)
    
    result <- 
      apply(anova_table, 1, function(x) {
        ci <- 
          apaTables::get.ci.partial.eta.squared(
            F.value = x["F"], df1 = x["num Df"], df2 = x["den Df"], conf.level = conf.level
          )
        
        return(setNames(c(ci$LL, ci$UL), c("LL", "UL")))
      }) %>% 
      t() %>% 
      as.data.frame()
    
    result$conf.level <- conf.level
    
    return(result)
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

#plotting
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
