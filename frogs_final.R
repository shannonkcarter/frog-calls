#############################################################
###-----------------PREP WORKSPACE------------------------###
#############################################################

## Set working directory and clear memory
setwd("C:\\Users\\Shannon\\Dropbox\\Phenology field data_Dan\\east TX frog call sampling\\shannons_script")
rm(list = ls(all = T))

## Load required packages
library(car)
library(tidyverse)
library(multcomp)
library(reshape2)
library(sfsmisc)
library(RColorBrewer)
library(ggExtra)
library(lattice)
library(gridExtra)
library(lme4)
library(lmerTest)
library(cowplot)

## Load universal plotting elements
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.minor = element_blank(), 
                 panel.grid.major = element_blank(),
                 axis.text  = element_text(size = rel(1.7),colour = "black"),
                 axis.title = element_text(size = rel(2.0)),
                 axis.line  = element_line(colour = "black"))

#############################################################
###----------DATA PROCESSING AND ANALYSIS-----------------###
#############################################################

###---PROCESS AND CLEAN RAW DATA ------------------------------

## Load raw data for each pond
p1 <- read.csv("CallingDataPond1.csv")
p2 <- read.csv("CallingDataPond2.csv")
p3 <- read.csv("CallingDataPond3.csv")
p4 <- read.csv("CallingDataPond4.csv")
p5 <- read.csv("CallingDataPond5.csv")
p6 <- read.csv("CallingDataPond6.csv")
p7 <- read.csv("CallingDataPond7.csv")
p8 <- read.csv("CallingDataPond8.csv")

## Stack raw data into one df
p18 <- rbind(p1, p2, p3, p4, p5, p6, p7, p8)

## Change to long format and rename columns
mp18 <- melt(p18, id.vars = c("POND", "YEAR", "FROGTIME","FROGDATE"),
             measure.vars = c("HV", "HC", "BV", "BW", "RCA", "RCL", "RS", "GC", "PC", "PT", "AC", "RP"))
names(mp18) <- c("pond", "year", "time", "date", "sp", "calls")

## Adjust time variable
time  <- paste(mp18$date, mp18$time)
time1 <- strptime(time,"%m/%d/%Y %H:%M")
time2 <- as.POSIXlt(time1)
time3 <- as.Date(time2, "%m/%d/%Y %H:%M")
mp18$time <- time3
mp18 <- subset(mp18, select = c("pond", "year", "time", "sp", "calls"))

###---SUMMARIZE BY DAY AND ADD CUMULATIVE COLUMN ------------

## Sum the 6 daily observations (takes several hours)
sums <- ddply(mp18, .(pond, year, time, sp), summarize, 
              dailysum = sum(calls))

## Make a function that treats NA days like 0
cumsumxNA <- function(x) {
  x[which(is.na(x))] = 0
  return(cumsum(x))
}

## Run daily data through modified cumsum function
daily <- ddply(sums, .(pond, sp, year), transform, 
               cumsum  = cumsumxNA(dailysum))

## Clean it up and write it
daily <- daily[daily$year != 2000 ,]
daily$doy <- as.numeric(strftime(daily$time, format = "%j"))
write.csv(daily, "dailycalls.csv")

###---ADD PREDICTOR METRICS (FIRST, MEAN, LAST) -------------

## Define first, middle, and last calls and return corresponding date
callperiods <- ddply(daily, .(pond, sp, year), summarise,
                     firstdate = doy[min(which(cumsum >= (min(cumsum) + 5)))],
                     meddate   = doy[min(which(cumsum >= (max(cumsum)/2)))],
                     lastdate  = doy[max(which(cumsum <= (max(cumsum) - 5)))],
                     firstcall = (min(cumsum) + 5),
                     medcall   = (max(cumsum)/2) + 0.1,
                     lastcall  = (max(cumsum) - 5))

## Correct so median call date doesn't return day 1 and write
callperiods$meddate[callperiods$firstdate %in% NA] <- NA
write.csv(callperiods, "callperiods.csv")

###---OVERLAP BY CALLING PERIOD FOR ALL SP PAIRS ------------

## Build function to make matrix of overlap vs. predictors for all y-p-s

overlapmatrix <- function(y, p, s1, s2) {
  
  ### GETTING OVERLAP
  
  ## Subset each species from daily calls data
  of1 <- subset(daily, subset = (sp == s1 & year == y & pond == p))
  of2 <- subset(daily, subset = (sp == s2 & year == y & pond == p))
  
  ## Lowess function for each
  lf1 <- lowess(of1$dailysum, f = 1/50, iter = 3, delta = 4)
  lf2 <- lowess(of2$dailysum, f = 1/50, iter = 3, delta = 4)
  
  ## Calculate overlap
  d     <- data.frame(day = lf1$x, frog1 = lf1$y, frog2 = lf2$y)
  d$min <- pmin(d$frog1, d$frog2)
  inter <- integrate.xy(d$day, d$min)
  prop1 <- inter/integrate.xy(d$day, d$frog1)
  prop2 <- inter/integrate.xy(d$day, d$frog2)
  
  ### GETTING PREDICTORS
  
  ## Subset each species from predictor data
  pf1 <- subset(periods, sp == s1 & year == y & pond == p, select = c("firstdate", "meddate", "lastdate"))
  pf2 <- subset(periods, sp == s2 & year == y & pond == p, select = c("firstdate", "meddate", "lastdate"))
  colnames(pf1) <- c("f1start", "f1med", "f1last")
  colnames(pf2) <- c("f2start", "f2med", "f2last")
  
  ## Designate and return output
  out <- cbind(prop1, prop2, pf1, pf2)
  return(out)
}

## Make a matrix of all year-pond combos and apply function to matrix
all <- as.matrix(expand.grid(unique(daily$year), unique(daily$pond), unique(daily$sp), unique(daily$sp)))
overmat <- as.data.frame(t(mapply(overlapmatrix, y = all[,1], p = all[,2], s1 = all[,3], s2 = all[,4])))

## Bind to all combos
overmat <- as.data.frame(cbind(all, overmat))
colnames(overmat) <- c("year", "pond", "sp1", "sp2", "prop1", "prop2", "f1start", "f1med", "f1last", "f2start", "f2med", "f2last")

## Unlist all columns and write csv
overmat <- cbind(overmat[!sapply(overmat, is.list)], 
                 (t(apply(overmat[sapply(overmat, is.list)], 1, unlist))))
write.csv(overmat, "rawovermat.csv")

## Calculate RELATIVE start and med diff from dates.
overmat$startdiff <- abs((overmat$f1start - overmat$f2start)/(overmat$f1last - overmat$f1start))
overmat$meddiff   <- abs((overmat$f1med   - overmat$f2med)/(overmat$f1last - overmat$f1start))

## Clear like sp pairs
overmat <- subset(overmat, startdiff != 0)            # filter out HV-HV pairs
overmat <- subset(overmat, prop1 >= 0 & prop1 <= 1)   # overlap should be between 0 and 1
overmat <- subset(overmat, prop2 >= 0 & prop2 <= 1)

## Clean it up and write it
all <- subset(overmat, select = c(year, pond, sp1, sp2, prop1, prop2, startdiff, meddiff))
colnames(all) <- c("year", "pond", "sp1", "sp2", "overlap", "overlapf2", "startdiff", "meddiff")
all <- na.omit(all)
all <- data.frame(lapply(all, as.character), stringsAsFactors = F)

all$overlap <- as.numeric(all$overlap)
all$spsp <- paste(all$sp1, all$sp2)

write.csv(all, "overmat.csv")

###---FILTERING OUT SP PAIRS WITH LOW OVERLAP -----------------

## Identify sp-pairs with >= 15 y-p's with >= 0.15 overlap (42/66)
filter <- subset(all, subset = (overlap >= 0.15))
freq   <- data.frame(table(filter$spsp))

## Eliminate pairs where one sp did not call
all$startdiff[which(all$startdiff == Inf)] = NA

## Select only species pairs that meet criteria
fall <- subset(all, spsp %in% c("AC BW", "AC HV", "AC RCL", "BV AC", "BV BW","BV GC", "BV HC", 
                                "BV HV", "BV RCL","BV RS", "BW AC", "BW BV", "BW GC", "BW HC", 
                                "BW HV", "BW RCL", "BW RS", "GC AC", "GC BV", "GC BW", "GC HC",
                                "GC HV", "GC RCL", "HC AC", "HC BV", "HC BW", "HC GC", "HC HV", 
                                "HC RCL", "HV AC","HV BV", "HV BW", "HV GC", "HV HC", "HV RCL", 
                                "HV RS", "PC RS", "RCA BW", "RCA HV", "RCA RCL", "RCL AC", 
                                "RCL HV", "RCL RS", "RP RS", "RS BW", "RS HV", "RS PC", "RS RCL"))
## Rename sp sp
fall$spsp <- factor(fall$spsp,
                       levels = c("AC BW", "AC HV", "AC RCL", "BV AC", "BV BW","BV GC", 
                                  "BV HC", "BV HV", "BV RCL","BV RS", "BW AC", "BW BV", 
                                  "BW GC", "BW HC", "BW HV", "BW RCL", "BW RS", "GC AC", 
                                  "GC BV", "GC BW", "GC HC", "GC HV", "GC RCL", "HC AC", 
                                  "HC BV", "HC BW", "HC GC", "HC HV", "HC RCL", "HV AC",
                                  "HV BV", "HV BW", "HV GC", "HV HC", "HV RCL", "HV RS", 
                                  "PC RS", "RCA BW", "RCA HV", "RCA RCL", "RCL AC", "RCL HV",
                                  "RCL RS", "RP RS", "RS BW", "RS HV", "RS PC", "RS RCL"),
                       
                       labels = c("ACr : BWo", "ACr : HVe", "ACr : RCl", "BVa : ACr", "BVa : BWo", "BVa : GCa", 
                                  "BVa : HCi", "BVa : HVe", "BVa : RCl", "BVa : RSp", "BWo : ACr", "BWo : BVa", 
                                  "BWo : GCa", "BWo : HCi", "BWo : HVe", "BWo : RCl", "BWo : RSp", "GCa : ACr", 
                                  "GCa : BVa", "GCa : BWo", "GCa : HCi", "GCa : HVe", "GCa : RCl", "HCi : ACr", 
                                  "HCi : BVa", "HCi : BWo", "HCi : GCa", "HCi : HVe", "HCi : RCl", "HVe : ACr", 
                                  "HVe : BVa", "HVe : BWo", "HVe : GCa", "HVe : HCi", "HVe : RCl", "HVe : RSp", 
                                  "PCr : RSp", "RCa : BWo", "RCa : HVe", "RCa : RCl", "RCl : ACr", "RCl : HVe", 
                                  "RCl : RSp", "RPa : RSp", "RSp : BWo", "RSp : HVe", "RSp : PCr", "RSp : RCl"))


fall$pond <- as.factor(fall$pond)
write.csv(fall, "filtered_all.csv")

###---START   ~ OVERLAP: FIND LINEAR SLOPE FOR EVERY SPSP-P ------

## Function to find regression coefficients for sp-p
regression <- function(sp, p) {
  
  ## Subset each species-pond
  s <- subset(fall, spsp == sp & pond == p)
  
  ## Linear model
  l  <- lm(overlap ~ startdiff, data = s)
  ls <- summary(l) 
  
  ## Extract regression coefficient
  c  <- round(coef(l)[2], 3)             # regression slope
  r  <- round(ls$r.squared, 3)           # r-squared
  pv <- round(ls$coefficients[2, 4], 3)  # p-value
  
  ## Log model
  #l1  <- lm(overlap ~ log(startdiff), data = s)
  #ls1 <- summary(l1)
  
  ## Extract regression coefficient
  #c1  <- round(coef(l1)[2], 3)            # regression slope
  #r1  <- round(ls1$r.squared, 3)          # r-squared
  #pv1 <- round(ls1$coefficients[2, 4], 3) # p-value
  
  ## Output
  out <- cbind(sp, p, c, r, pv)#, c1, r1, pv1)
  return(out)
}

## Make a table to show number of years for each pond-spsp and choose only ones with >=4
e  <- as.data.frame(with(fall, table(spsp, pond)))
e1 <- subset(e, subset = e$Freq >= 4)

## Applying function to data and paste to sp-y
regmat <- as.data.frame(t(mapply(regression, sp = e1$spsp, p = e1$pond)))
regmat <- as.data.frame(cbind(e1, regmat))

## Clean it up and write csv
regmat <- subset(regmat, select = c("spsp", "pond", "V3", "V4", "V5"))#, "V6", "V7", "V8"))
colnames(regmat) <- c("spsp", "pond", "coef", "r2", "pval")#, "coef_log", "r2_log", "pval_log")
regmat$coef <- as.numeric(as.character(regmat$coef))
regmat$r2   <- as.numeric(as.character(regmat$r2))
regmat$pval <- as.numeric(as.character(regmat$pval))
#regmat$coef_log <- as.numeric(as.character(regmat$coef_log))
#regmat$r2_log   <- as.numeric(as.character(regmat$r2_log))
#regmat$pval_log <- as.numeric(as.character(regmat$pval_log))
regmat$pond  <- as.factor(regmat$pond)
write.csv(regmat, "start_over.csv")

###---MEDIAN  ~ OVERLAP: FIND LINEAR SLOPE FOR EVERY SPSP-P -------

## Filter out infinity values
fall$meddiff[which(fall$meddiff == Inf)] = NA

## Function to find regression coefficients for sp-p
reg_med <- function(sp, p) {
  
  ## Subset each species-pond
  s <- subset(fall, spsp == sp & pond == p)
  
  ## Linear model
  l  <- lm(overlap ~ meddiff, data = s)
  ls <- summary(l) 
  
  ## Extract regression coefficient
  c  <- round(coef(l)[2], 3)             # regression slope
  r  <- round(ls$r.squared, 3)           # r-squared
  pv <- round(ls$coefficients[2, 4], 3)  # p-value
  
  ## Output
  out <- cbind(sp, p, c, r, pv)
  return(out)
}

e  <- as.data.frame(with(fall, table(spsp, pond)))
e1 <- subset(e, subset = e$Freq >= 4)

## Applying function to data and paste to sp-y
regmat_med <- as.data.frame(t(mapply(reg_med, sp = e1$spsp, p = e1$pond)))
regmat_med <- as.data.frame(cbind(e1, regmat_med))

## Clean it up and write csv
regmat_med <- subset(regmat_med, select = c("spsp", "pond", "V3", "V4", "V5"))
colnames(regmat_med) <- c("spsp", "pond", "coef", "r2", "pval")
regmat_med$coef <- as.numeric(as.character(regmat_med$coef))
regmat_med$r2   <- as.numeric(as.character(regmat_med$r2))
regmat_med$pval <- as.numeric(as.character(regmat_med$pval))
regmat_med$pond <- as.factor(regmat_med$pond)
write.csv(regmat_med, "med_over.csv")

###---START   ~ OVERLAP: FIND LINEAR SLOPE FOR EVERY SPSP, POND AS RE ----------

## Pond has to be a factor
fall$pond <- as.factor(fall$pond)

## Finding regression coefficients with pond as a RE
reg_re <- function(sp) {
  
  ## Subset each species-pond
  s <- subset(fall, spsp == sp)
  
  ## Linear model
  l  <- lmer(overlap ~ startdiff + (1|pond), data = s) #need lme4 library make sure pond is a factor
  ls <- summary(l) 
  la <- anova(l)
  
  ## Extract regression coefficient
  c  <- summary(l)$coefficients[2, 1]   # regression slope
  se <- ls$coefficients[2,2]            # standard error
  p  <- la$`Pr(>F)`                     # p-value
  
  ## Output
  out <- cbind(sp, c, se, p)
  return(out)
}

reg_re("HCi : HVe")

## Applying function to data and paste to sp
reg_remat <- as.data.frame(t(mapply(reg_re, sp = fall$spsp)))
reg_remat <- as.data.frame(cbind(fall$spsp, reg_remat))

## Clean it up and write csv
reg_remat <- subset(reg_remat, select = c("fall$spsp", "V2", "V3", "V4"))
colnames(reg_remat) <- c("spsp", "coef", "se", "pval")
reg_remat$coef <- as.numeric(as.character(reg_remat$coef))
reg_remat$se <- as.numeric(as.character(reg_remat$se))
reg_remat$pval <- as.numeric(as.character(reg_remat$pval))
reg_remat <- unique(reg_remat)
write.csv(reg_remat, "start_over_re.csv")

###---MEDIAN  ~ OVERLAP: FIND LINEAR SLOPE FOR EVERY SPSP, POND AS RE ---------

## Pond has to be a factor
fall$pond <- as.factor(fall$pond)
fall$meddiff[is.infinite(fall$meddiff) | is.nan(fall$meddiff) ] <- NA 

## Finding regression coefficients with pond as a RE
reg_med_re <- function(sp) {
  
  ## Subset each species-pond
  s <- subset(fall, spsp == sp)
  
  ## Linear model
  l  <- lmer(overlap ~ meddiff + (1|pond), data = s) #need lme4 library make sure pond is a factor
  ls <- summary(l) 
  la <- anova(l)
  
  ## Extract regression coefficient
  c  <- summary(l)$coefficients[2, 1]   # regression slope
  se <- ls$coefficients[2,2]            # standard error
  p  <- la$`Pr(>F)`                     # p-value
  
  ## Output
  out <- cbind(sp, c, se, p)
  return(out)
}

reg_med_re("BVa : HVe")

## Applying function to data and paste to sp
reg_med_remat <- as.data.frame(t(mapply(reg_med_re, sp = fall$spsp)))
reg_med_remat <- as.data.frame(cbind(fall$spsp, reg_med_remat))

## Clean it up and write csv
reg_med_remat <- subset(reg_med_remat, select = c("fall$spsp", "V2", "V3", "V4"))
colnames(reg_med_remat) <- c("spsp", "coef", "se", "pval")
reg_med_remat$coef <- as.numeric(as.character(reg_med_remat$coef))
reg_med_remat$se <- as.numeric(as.character(reg_med_remat$se))
reg_med_remat$pval <- as.numeric(as.character(reg_med_remat$pval))
reg_med_remat <- unique(reg_med_remat)
write.csv(reg_med_remat, "med_over_re.csv")

###---START   ~ YEAR: FIND LINEAR SLOPE FOR EVERY SPSP-P ----------------
starttrend <- function(sp, p) {
  
  ## Subset each species-pond
  s <- subset(fall, spsp == sp & pond == p)
  
  ## Linear model
  l  <- lm(startdiff ~ year, data = s)
  ls <- summary(l) 
  
  ## Extract regression coefficient
  c  <- round(coef(l)[2], 3)             # regression slope
  r  <- round(ls$r.squared, 3)           # r-squared
  se <- ls$coefficients[2,2]             # standard error
  pv <- round(ls$coefficients[2, 4], 3)  # p-value
  
  ## Log model
  l1  <- lm(log(startdiff) ~ year, data = s)
  ls1 <- summary(l1)
  
  ## Extract coefficients
  c1  <- round(coef(l1)[2], 3)            # regression slope
  r1  <- round(ls1$r.squared, 3)          # r-squared
  pv1 <- round(ls1$coefficients[2, 4], 3) # p-value
  
  ## Output
  out <- cbind(sp, p, c, r, se, pv, c1, r1, pv1)
  return(out)
}

## Make a table to show number of years for each pond-spsp and choose only ones with >=4
e  <- as.data.frame(with(fall, table(spsp, pond)))
e1 <- subset(e, subset = e$Freq >= 4)

## Applying function tof data and paste to sp-y
reg_starttrend <- as.data.frame(t(mapply(starttrend, sp = e1$spsp, p = e1$pond)))
reg_starttrend <- as.data.frame(cbind(e1, reg_starttrend))

## Clean it up and write csv
reg_starttrend <- subset(reg_starttrend, select = c("spsp", "pond", "V3", "V4", "V5", "V6", "V7", "V8", "V9"))
colnames(reg_starttrend) <- c("spsp", "pond", "coef", "r2", "se", "pval", "coef_log", "r2_log", "pval_log")
reg_starttrend$coef <- as.numeric(as.character(reg_starttrend$coef))
reg_starttrend$r2   <- as.numeric(as.character(reg_starttrend$r2))
reg_starttrend$se <- as.numeric(as.character(reg_starttrend$se))
reg_starttrend$pval <- as.numeric(as.character(reg_starttrend$pval))
reg_starttrend$coef_log <- as.numeric(as.character(reg_starttrend$coef_log))
reg_starttrend$r2_log   <- as.numeric(as.character(reg_starttrend$r2_log))
reg_starttrend$pval_log <- as.numeric(as.character(reg_starttrend$pval_log))
reg_starttrend$pond <- as.factor(reg_starttrend$pond)
write.csv(reg_starttrend, "start_year.csv")


###---OVERLAP ~ YEAR: FIND LINEAR SLOPE FOR EVERY SPSP-P ---------------
overlaptrend <- function(sp, p) {
  
  ## Subset each species-pond
  s <- subset(fall, spsp == sp & pond == p)
  
  ## Linear model
  l  <- lm(overlap ~ year, data = s)
  ls <- summary(l) 
  
  ## Extract regression coefficient
  c  <- round(coef(l)[2], 3)             # regression slope
  r  <- round(ls$r.squared, 3)           # r-squared
  pv <- round(ls$coefficients[2, 4], 3)  # p-value
  
  ## Output
  out <- cbind(sp, p, c, r, pv)
  return(out)
}

## Make a table to show number of years for each pond-spsp and choose only ones with >=4
e  <- as.data.frame(with(fall, table(spsp, pond)))
e1 <- subset(e, subset = e$Freq >= 4)

## Applying function tof data and paste to sp-y
reg_overlaptrend <- as.data.frame(t(mapply(overlaptrend, sp = e1$spsp, p = e1$pond)))
reg_overlaptrend <- as.data.frame(cbind(e1, reg_overlaptrend))

## Clean it up and write csv
reg_overlaptrend <- subset(reg_overlaptrend, select = c("spsp", "pond", "V3", "V4", "V5"))
colnames(reg_overlaptrend) <- c("spsp", "pond", "coef", "r2", "pval")
reg_overlaptrend$coef <- as.numeric(as.character(reg_overlaptrend$coef))
reg_overlaptrend$r2   <- as.numeric(as.character(reg_overlaptrend$r2))
reg_overlaptrend$pval <- as.numeric(as.character(reg_overlaptrend$pval))
reg_overlaptrend$pond <- as.factor(reg_overlaptrend$pond)
write.csv(reg_overlaptrend, "over_year.csv")


###---START   ~ YEAR: FIND LINEAR SLOPE FOR EVERY SPSP, POND AS RE --------------

## Finding regression coefficients with pond as a RE
starttrend_re <- function(sp) {
  
  ## Subset each species-pond
  s <- subset(fall, spsp == sp)
  
  ## Linear model
  l  <- lmer(startdiff ~ year + (1|pond), data = s)
  ls <- summary(l) 
  la <- anova(l)
  
  ## Extract regression coefficient
  c  <- summary(l)$coefficients[2, 1]   # regression slope
  se <- ls$coefficients[2,2]            # standard error
  p  <- la$`Pr(>F)`                     # p-value
  
  ## Log model
  l1  <- lmer(log(startdiff) ~ year + (1|pond), data = s)
  ls1 <- summary(l1)
  la1 <- anova(l1)
  
  ## Extract regression coefficients
  c1 <- summary(l1)$coefficients[2, 1] # regression slope
  se1 <- ls1$coefficients[2,2]         # standard error
  p1 <- la1$`Pr(>F)`                   # p-value
  
  ## Output
  out <- cbind(sp, c, se, p, c1, se1, p1)
  return(out)
}

## Applying function to data and paste to sp
reg_starttrend_re <- as.data.frame(t(mapply(starttrend_re, sp = fall$spsp)))
reg_starttrend_re <- as.data.frame(cbind(fall$spsp, reg_starttrend_re))

## Clean it up and write csv
reg_starttrend_re <- subset(reg_starttrend_re, select = c("fall$spsp", "V2", "V3", "V4", "V5", "V6", "V7"))
colnames(reg_starttrend_re) <- c("spsp", "coef", "se", "pval", "coef_log", "se_log", "pval_log")
reg_starttrend_re$coef <- as.numeric(as.character(reg_starttrend_re$coef))
reg_starttrend_re$se   <- as.numeric(as.character(reg_starttrend_re$se))
reg_starttrend_re$pval <- as.numeric(as.character(reg_starttrend_re$pval))
reg_starttrend_re$coef_log <- as.numeric(as.character(reg_starttrend_re$coef_log))
reg_starttrend_re$se_log <- as.numeric(as.character(reg_starttrend_re$se_log))
reg_starttrend_re$pval_log <- as.numeric(as.character(reg_starttrend_re$pval_log))
reg_starttrend_re <- unique(reg_starttrend_re)
write.csv(reg_starttrend_re, "start_year_re.csv")

###---MEDIAN  ~ YEAR: FIND LINEAR SLOPE FOR EVERY SPSP, POND AS RE --------------

## Pond has to be a factor
fall$meddiff[is.infinite(fall$meddiff) | is.nan(fall$meddiff) ] <- NA 
fall$meddiff[fall$meddiff == 0] <- NA

## Finding regression coefficients with pond as a RE
medtrend_re <- function(sp) {
  
  ## Subset each species-pond
  s <- subset(fall, spsp == sp)
  
  ## Linear model
  l  <- lmer(meddiff ~ year + (1|pond), data = s)
  ls <- summary(l) 
  la <- anova(l)
  
  ## Extract regression coefficient
  c  <- summary(l)$coefficients[2, 1]   # regression slope
  se <- ls$coefficients[2,2]            # standard error
  p  <- la$`Pr(>F)`                     # p-value
  
  ## Log model
  l1  <- lmer(log(meddiff) ~ year + (1|pond), data = s)
  ls1 <- summary(l1)
  la1 <- anova(l1)
  
  ## Extract regression coefficients
  c1 <- summary(l1)$coefficients[2, 1] # regression slope
  se1 <- ls1$coefficients[2,2]         # standard error
  p1 <- la1$`Pr(>F)`                   # p-value
  
  ## Output
  out <- cbind(sp, c, se, p, c1, se1, p1)
  return(out)
}

## Applying function to data and paste to sp
reg_medtrend_re <- as.data.frame(t(mapply(medtrend_re, sp = fall$spsp)))
reg_medtrend_re <- as.data.frame(cbind(fall$spsp, reg_medtrend_re))

## Clean it up and write csv
reg_medtrend_re <- subset(reg_medtrend_re, select = c("fall$spsp", "V2", "V3", "V4", "V5", "V6", "V7"))
colnames(reg_medtrend_re) <- c("spsp", "coef", "se", "pval", "coef_log", "se_log", "pval_log")
reg_medtrend_re$coef <- as.numeric(as.character(reg_medtrend_re$coef))
reg_medtrend_re$se   <- as.numeric(as.character(reg_medtrend_re$se))
reg_medtrend_re$pval <- as.numeric(as.character(reg_medtrend_re$pval))
reg_medtrend_re$coef_log <- as.numeric(as.character(reg_medtrend_re$coef_log))
reg_medtrend_re$se_log <- as.numeric(as.character(reg_medtrend_re$se_log))
reg_medtrend_re$pval_log <- as.numeric(as.character(reg_medtrend_re$pval_log))
reg_medtrend_re <- unique(reg_medtrend_re)
write.csv(reg_medtrend_re, "med_year_re.csv")

###---OVERLAP ~ YEAR: FIND LINEAR SLOPE FOR EVERY SPSP, POND AS RE -----------

## Finding regression coefficients with pond as a RE
overlaptrend_re <- function(sp) {
  
  ## Subset each species-pond
  s <- subset(fall, spsp == sp)
  
  ## Linear model
  l  <- lmer(overlap ~ year + (1|pond), data = s)
  ls <- summary(l) 
  la <- anova(l)
  
  ## Extract regression coefficient
  c  <- summary(l)$coefficients[2, 1]   # regression slope
  se <- ls$coefficients[2, 2]           # standard error
  p  <- la$`Pr(>F)`                     # p-value
  
  ## Output
  out <- cbind(sp, c, se, p)
  return(out)
}
overlaptrend_re("HCi : HVe")

## Testing function
h <- subset(fall, spsp == "HCi : HVe")
m <- lmer(overlap ~ year + (1|pond), data = h)
summary(m)

## Applying function to data and paste to sp
reg_overlaptrend_re <- as.data.frame(t(mapply(overlaptrend_re, sp = fall$spsp)))
reg_overlaptrend_re <- as.data.frame(cbind(fall$spsp, reg_overlaptrend_re))

## Clean it up and write csv
reg_overlaptrend_re <- subset(reg_overlaptrend_re, select = c("fall$spsp", "V2", "V3", "V4"))
colnames(reg_overlaptrend_re) <- c("spsp", "coef", "se", "pval")
reg_overlaptrend_re$coef <- as.numeric(as.character(reg_overlaptrend_re$coef))
reg_overlaptrend_re$se <- as.numeric(as.character(reg_overlaptrend_re$se))
reg_overlaptrend_re$pval <- as.numeric(as.character(reg_overlaptrend_re$pval))
reg_overlaptrend_re <- unique(reg_overlaptrend_re)
write.csv(reg_overlaptrend_re, "over_year_re.csv")

###---WEATHER DATA ----------------------------------------------------

## Read and clean raw weather data
weather <- read.csv("WeatherData_long.csv", header = T)
weather <- na.omit(weather)
weather$date <- as.Date(weather$date, "%m/%d/%Y")
weather$doy  <- as.numeric(strftime(weather$date, format = "%j"))
weather$temp <- as.numeric(as.character(weather$temp))
weather$year <- as.factor(weather$year)

## Write a custom function that treats NA like 0
cumrainxNA <- function(x) {
  x[which(is.na(x))] = 0
  return(cumsum(x))
}

## Run daily data through modified cumsum function
weather <- ddply(weather, .(year, site), transform, 
               cumrain = cumrainxNA(rain),
               avgtemp = mean(temp))

write.csv(weather, "weather.csv")


###---MERGING WEATHER AND OVERLAP DATA ----------------------------------

## Subset relevant columns of overlap data
overlap <- subset(fall, select = c(year, pond, overlap, spsp))
overlap$pond <- as.numeric(overlap$pond)

## Add a site column
overlap$site <- ifelse(overlap$pond <= 4, "sf",
                       ifelse(overlap$pond >= 5, "dc",
                              NA))

## Average overlap values for 4 ponds at SF & DC
overlap <- ddply(overlap, .(year, spsp, site), transform,
                 avgover = mean(overlap))

## Subset relevant columns and delete duplicate rows
overlap <- subset(overlap, select = c(year, spsp, site, avgover))
overlap <- unique(overlap)

## Merging weather summary and overlap data
weather_over <- merge(overlap, weather_sum)
write.csv(weather_over, "weather_over.csv")


###---WEATHER ~ OVERLAP: FIND LINEAR SLOPE FOR EVERY SPSP-S ------------------

## Function to find regression coefficients for sp-p
regression <- function(sp, si) {
  
  ## Subset each species-pond
  s <- subset(weather_over, subset = (spsp == sp & site == si))
  
  ## Rain model
  l  <- lm(avgover ~ cumrain, data = s)
  ls <- summary(l) 
  
  ## Extract regression coefficient
  cr <- round(coef(l)[2], 6)              # regression slope
  rr <- round(ls$r.squared, 3)            # r-squared
  pr <- round(ls$coefficients[2, 4], 3)   # p-value
  
  ## Temp model
  l1  <- lm(avgover ~ avgtemp, data = s)
  ls1 <- summary(l1)
  
  ## Extract regression coefficient
  ct  <- round(coef(l1)[2], 6)             # regression slope
  rt  <- round(ls1$r.squared, 3)           # r-squared
  pt  <- round(ls1$coefficients[2, 4], 3)  # p-value
  
  ## Output
  out <- cbind(sp, si, cr, rr, pr, ct, rt, pt)
  return(out)
}

## Make a table to show number of years for each pond-spsp and choose only ones with >=4
e  <- as.data.frame(with(weather_over, table(spsp, site)))
e1 <- subset(e, subset = e$Freq >= 6)

## Applying function to data and paste to sp-y
regmat <- as.data.frame(t(mapply(regression, sp = e1$spsp, s = e1$site)))
regmat <- as.data.frame(cbind(e1, regmat))

## Clean it up and write csv
regmat <- subset(regmat, select = c("spsp", "site", "V3", "V4", "V5", "V6", "V7", "V8"))
colnames(regmat) <- c("spsp", "site", "coef_rain", "r2_rain", "pval_rain", "coef_temp", "r2_temp", "pval_temp")
regmat$coef_rain <- as.numeric(as.character(regmat$coef_rain))
regmat$r2_rain   <- as.numeric(as.character(regmat$r2_rain))
regmat$pval_rain <- as.numeric(as.character(regmat$pval_rain))
regmat$coef_temp <- as.numeric(as.character(regmat$coef_temp))
regmat$r2_temp   <- as.numeric(as.character(regmat$r2_temp))
regmat$pval_temp <- as.numeric(as.character(regmat$pval_temp))
write.csv(regmat, "weather_over_coef.csv")


###---HELPER FUNCITON TO GET MEANS AND SES ----------------

summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = F,
                      conf.interval = .95, .drop = T) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm = FALSE) {
    if (na.rm) sum(!is.na(x)) 
    else       length(x)
  }
  
  # The summary. For each group's df, return a vector with N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Run start ~ overlap coef through helper
start_over_err <- summarySE(start_over, measurevar = "coef", groupvars = c("spsp"))
write.csv(start_over_err,  "start_over_err.csv")

## Run med ~ overlap coef through helper
med_over_err <- summarySE(med_over, measurevar = "coef", groupvars = c("spsp"))
write.csv(med_over_err,  "med_over_err.csv")

## Run start ~ year coef through helper
start_year_err <- summarySE(start_year, measurevar = "coef_log", groupvars = c("spsp"))
write.csv(start_year_err,  "start_year_err.csv")

## Run overlap ~ year coef through helper
over_year_err  <- summarySE(over_year, measurevar = "coef", groupvars = c("spsp"))
write.csv(over_year_err,  "over_year_err.csv")

## Get mean and SE for each phenological period
p <- na.omit(periods)
p <- subset(p, subset = (duration > 0))
med_var <- summarySE (p, measurevar = "meddate", groupvars = "sp")
period_var <- summarySE(p, measurevar = "duration", groupvars = "sp")
start_var <- summarySE(p, measurevar = "firstdate", groupvars = "sp")



#############################################################
###---------------------LOAD DATA-------------------------###
#############################################################

## DAILY SUM AND YEARLY CUMULATIVE CALL FOR ALL Y-P-S
daily <- read.csv("dailycalls.csv", header = T)

## FIRST/MIDDLE/LAST/DURATION DATES FOR EACH Y-P-S
periods <- read.csv("callperiods.csv", header = T)
periods <- subset(periods, select = c("pond", "sp", "year", "firstdate", "meddate", "lastdate"))
periods$duration <- as.numeric(periods$lastdate - periods$firstdate)

## UNFILTERED ALL SP OVERLAP MATRIX
all <- read.csv("overmat.csv")

## FILTERED SP OVERLAP MATRIX
fall <- read.csv("filtered_all.csv", header = T)
fall$pond <- as.factor(fall$pond)

## METRIC ~ OVERLAP: ALL REGRESSION MODELS
start_over <- read.csv("start_over.csv", header = T)         ## START ~ OVERLAP
med_over   <- read.csv("med_over.csv", header = T)           ## MED ~ OVERLAP
start_over_re <- read.csv("start_over_re.csv", header = T)   ## START ~ OVERLAP; POND AS RE
med_over_re   <- read.csv("med_over_re.csv", header = T)     ## MED ~ OVERLAP; POND AS RE

## METRIC ~ OVERLAP: REGRESSION MODELS WITH MEANS/SE
start_over_err <- read.csv("start_over_err.csv",  header = T) ## START ~ OVERLAP
med_over_err   <- read.csv("med_over_err.csv", header = T)    ## MED ~ OVERLAP

## PHENOLOGY ~ YEAR: ALL REGRESSION MODELS
start_year <- read.csv("start_year.csv", header = T)       ## START ~ YEAR
over_year  <- read.csv("over_year.csv", header = T)        ## OVERLAP ~ YEAR
start_year_re <- read.csv("start_year_re.csv", header = T) ## START ~ YEAR; POND AS RE
over_year_re  <- read.csv("over_year_re.csv", header = T)  ## OVERLAP ~ YEAR; POND AS RE
med_year_re   <- read.csv("med_year_re.csv", header = T)   ## MEDIAN ~ YEAR; POND AS RE

## PHENOLOGY ~ YEAR: REGRESSION MODELS WITH MEANS/SE
start_year_err <- read.csv("start_year_err.csv",  header = T) ## START ~ YEAR
over_year_err  <- read.csv("over_year_err.csv", header = T)   ## OVERLAP ~ YEAR

## RAW WEATHER DATA
weather <- read.csv("weather.csv", header = T)         ## WEATHER DATA WITH CUMULATIVE RAIN
weather_sum <- read.csv("weather_sum.csv", header = T) ## TOTAL YEARLY RAIN FOR EACH SITE

## WEATHER X OVERLAP DATA
weather_over <- read.csv("weather_over.csv", header = T)           ## YEARLY WEATHER DATA WITH SPSP OVERLAP
weather_over_coef <- read.csv("weather_over_coef.csv", header = T) ## OVERLAP ~ WEATHER COEFFICIENTS

#############################################################
#############################################################
###---------------------PLOTS-----------------------------###
#############################################################
#############################################################

#############################################################
###-----------VISUALIZING THE RAW DISTRIBUTIONS-----------###
#############################################################

###---FUNCTION FOR PLOTTING Y-P-SP COMBOS -------------------

### Lowess function- tight fitting and used for overlap proportion value
overlapplot <- function(s1, s2, y, p) {
  
  ## Subset two species of interest
  s1 <- subset(daily, subset = (sp == s1 & year == y & pond == p))
  s2 <- subset(daily, subset = (sp == s2 & year == y & pond == p))
  
  ## Run lowess functions
  l1 <- as.data.frame(lowess(s1$dailysum, f = 1/20, iter = 3, delta = 4))
  l2 <- as.data.frame(lowess(s2$dailysum, f = 1/20, iter = 3, delta = 4))
  l  <- data.frame(day = l1$x, s1 = l1$y, s2 = l2$y)
  
  ## Plot
  p <- ggplot(l, aes(x = day, y = max(s1, s2))) + mytheme +
    xlab("Day of year") + ylab("Number of calling frogs") +
    geom_line(size = 2.5, colour = "darkorange", aes(x = day, y = s1, stat = "identity", position = "dodge")) +
    geom_line(size = 2.5, colour = "blue4",      aes(x = day, y = s2, stat = "identity", position = "dodge")) +
    geom_ribbon(aes(x = day, ymax = s1), ymin = 0, fill = "orangered", alpha = 0.5) +
    geom_ribbon(aes(x = day, ymax = s2), ymin = 0, fill = "navy",     alpha = 0.5) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
          axis.text.y = element_text(size = 40),
          axis.text.x = element_text(size = 40),
          axis.title.y = element_text(size = 40),
          axis.title.x = element_text(size = 40))
  
  ## Calculate areal overlap of lowess functions
  d     <- data.frame(day = l2$x, s1 = l1$y, s2 = l2$y)
  d$min <- pmin(d$s1, d$s2)
  inter <- integrate.xy(d$day, d$min)
  
  ## Print intersection and plot
  prop1 <- print(inter/integrate.xy(d$day, d$s1))
  prop2 <- print(inter/integrate.xy(d$day, d$s2))
  print(p)
}

overlapplot("RCL", "HV", 2007, 4) 

## Smoother plots, but y axis not indicative of abundance, so misleading. could be a way to id peaks.
o <- function(s1, s2, y, p) {
  
  ## Subset two sepcies of interest
  s_1 <- subset(daily, select = c("pond", "year", "sp", "dailysum", "doy"))
  s_2 <- subset(daily, select = c("pond", "year", "sp", "dailysum", "doy"))
  
  ## Remove zero entries for geom_density
  s_1 <- subset(daily, subset = (sp == s1 & year == y & pond == p & dailysum != 0))
  s_2 <- subset(daily, subset = (sp == s2 & year == y & pond == p & dailysum != 0))
  
  ## Plot
  p <- ggplot(s_1, aes(x = doy)) + mytheme + xlim(0, 365) +
    geom_density(color = "darkorange", fill = "orangered", alpha = 0.5, size = 2, adjust = 1/3) +
    geom_density(data = s_2, aes(x = doy), 
                 color = "blue4", fill = "navy", alpha = 0.5, size = 2, adjust = 1/3)
  
  ## Print plot
  print(p)
}

o("RCL", "HV", 2015, 4) 


###---FUNCTIONS FOR PLOTTING SINGLE POPULATION --------------

## Wrapped by year, choose pond and species
years <- function(s, p) {
  s <- subset(daily, sp == s & pond == p)
  
  p <- ggplot(s, aes(x = doy, y = dailysum)) + mytheme +
    xlab(NULL) + ylab (NULL) +
    #geom_point() +
    stat_smooth(colour = "blue4", se = F, span = 0.2, size = 1.5) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks  = element_blank(),
          strip.text.x = element_text(size = 15)) +
    facet_wrap(~year) #+ ylim(0, 4)
  
  return(p)
}

years("HV", 5)

## All ponds shown, wrapped by year, choose species
yearp <- function(s) {
  s <- subset(daily, sp == s)
  
  p <- ggplot(s, aes(x = doy, y = dailysum, colour = factor(pond))) + mytheme +
    xlab(NULL) + ylab (NULL) +
    #geom_point() +
    stat_smooth(span = 0.1, size = 1.5, se = F, alpha = 0.5) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks  = element_blank(),
          strip.text.x = element_text(size = 15)) +
    facet_wrap(~year)
  
  return(p)
}

yearp("HV")

## All years shown, wrapped by pond, choose sp
yeary <- function(s) {
  s <- subset(daily, sp == s)
  
  p <- ggplot(s, aes(x = doy, y = dailysum, colour = factor(year))) + mytheme +
    xlab(NULL) + ylab (NULL) +
    #geom_point() +
    stat_smooth(span = 0.1, size = 1.5, se = F, alpha = 0.5) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks  = element_blank(),
          strip.text.x = element_text(size = 15)) +
    facet_wrap(~pond)
  
  return(p)
}

yeary("RCL")

## All sp shown, wrapped by year, choose pond
yearsp <- function(p) {
  s <- subset(daily, pond == p & sp != "RCL") # taking out RCL since it's so much more abundant
  
  p <- ggplot(s, aes(x = doy, y = dailysum, color = factor(sp))) + mytheme +
    xlab(NULL) + ylab (NULL) +
    stat_smooth(span = 0.1, size = 1.5, se = F, alpha = 0.5) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks  = element_blank(),
          strip.text.x = element_text(size = 15)) +
    facet_wrap(~year)
  
  return(p)
}

yearsp(5)

## Not informative to do 2 sp here, because abundance/overlap disappears. only tells you about the shape







#############################################################
###------------OVERLAP ~ METRIC SCATTER PLOTS-------------###
#############################################################

###---POOLING ALL DATA ----------------------------------

## Start by overlap- pooling all data
so <- ggplot(fall, aes(x = log10(startdiff), y = overlap)) + mytheme +
  geom_point(size = 4, color = "black", alpha = 0.25) +
  #stat_smooth(method = lm, size = 3, se = F) +
  stat_smooth(method = lm, size = 2, se = F, aes(color = factor(spsp))) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        legend.position = "none") + ylim(0, 1) +
  ylab("Overlap in calling phenology") + xlab("log(Difference in first call date (days))")
so

## Median by overlap- pooling all data
mo <- ggplot(fall, aes(x = log(meddiff), y = overlap)) + mytheme +
  geom_point(size = 4, alpha = 0.25) + 
  #stat_smooth(method = lm, size = 3, se = F) + 
  stat_smooth(method = lm, size = 2, se = F, aes(color = factor(spsp))) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
  xlab("log(Difference in median call date (days))") + ylab("Overlap in calling phenology") +
  ylim(0, 1)
mo

## Isolating a well-represented pair for presentations
so_hvrcl <- ggplot(subset(fall, spsp == "HVe : RCl"), 
                aes(x = log10(startdiff), y = overlap, colour = pond)) + mytheme +
  geom_point(size = 4, alpha = 0.75) +
  stat_smooth(method = lm, size = 3, se = F) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        legend.background = element_rect(size = 1, colour = "black"),
        legend.position = c(0.9, 0.7),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15)) +
  ylab("Overlap in calling phenology") + xlab("log(Difference in first call date (days))") 

so_hvrcl


###---ALL DATA, WRAPPED BY FACTORS --------------------------

## Start ~ overlap, wrapped by spsp
so_sp <- ggplot(fall, aes(x = log(startdiff), y = overlap)) + #,color = factor(pond))) +
  geom_point(size = 2, aes(color = factor(pond))) + 
  stat_smooth(method = lm, size = 2, se = T, alpha = 0.5) + #, 
              #aes(color = factor(pond))) +
  facet_wrap(~ spsp) + mytheme +
  xlab("Difference in start date") + ylab("Overlap") +
  ylim(0, 1)
so_sp

## Median ~ overlap, wrapped by spsp
mo_sp <- ggplot(fall, aes(x = log(meddiff), y = overlap, color = factor(pond))) +
  geom_point(size = 2) + 
  stat_smooth(method = lm, size = 1, se = F) +
  facet_wrap(c("spsp")) + mytheme +
  #theme(legend.position = "none") +
  xlab("Difference in median date") + ylab("Overlap") +
  ylim(0, 1)
mo_sp

mo <- ggplot(fall, aes(x = log(meddiff), y = overlap)) +
  geom_point(size = 2) + 
  stat_smooth(method = lm, size = 1, se = F) +
  facet_wrap(c("spsp")) + mytheme +
  #theme(legend.position = "none") +
  xlab("Difference in median date") + ylab("Overlap") +
  ylim(0, 1)
mo





#############################################################
###------------OVERLAP ~ METRIC SUMMARY PLOTS-------------###
#############################################################

###---FOR START DIFF ---------------------------------------

## Subset significant regressions so we can highlight in plots
start_over_sig <- subset(start_over_re, subset = start_over_re$pval <= 0.05)

## Reverse levels of factor spsp so it appears in abc order
start_over$spsp = with(start_over, factor(spsp, levels = rev(levels(spsp))))
start_over$pond <- as.factor(start_over$pond)

## To export to ms, save as image 500 x 750
so <- ggplot(start_over, aes(x = coef, y = spsp)) + mytheme +
  # all data
  geom_point(size = 3, alpha = 0.5, aes(color = pond)) + 
  # averages
  geom_point(shape = 18, alpha = 1, size = 4,
             data = start_over_err, aes(x = coef, y = spsp)) +
  # error bars on averages
  geom_errorbarh(data = start_over_err, size = 1,
                 aes(xmin = coef - se, xmax = coef + se, height = 0)) +
  # significant regressions
  #geom_point(shape = 17, size = 6, colour = "red",
  #data = start_over_sig, aes(x = coef, y = spsp)) +
  # pond as RE
  #geom_point(size = 6, colour = "blue", alpha = 0.25,
  #          data = start_over_re, aes(x = coef, y = spsp)) +
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono", size = 12),
        legend.position = "NULL", #c(0.9, 0.7),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 15)) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) + xlim(-1.75, 1.2)+
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  xlab("Regression slope \n (numerical overlap ~ difference in start date)") + ylab("Species pair")

so <- ggMarginal(so, size = 7, margins = "x", type = "histogram",
                 col = "black", fill = "black", alpha = 0.75)
so

so1 <- ggplot(start_over_err, aes(x = coef, y = reorder(spsp, coef))) + mytheme +
  geom_point(size = 4, colour = "black", shape = 18) + 
  ## SE bars on averages
  geom_errorbarh(data = start_over_err, size = 1, aes(xmin = coef - se, xmax = coef + se, height = 0)) +
  ## Ponds shown independently
  geom_point(data = start_over, aes(x = coef, y = spsp, color = pond), 
             size = 3, alpha = 0.5) + 
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono", size = 12),
        legend.position = "NULL", #c(0.9, 0.7),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 15)) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) + xlim(-1.75, 1.2)+
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  xlab("Regression slope \n (numerical overlap ~ difference in start date)") + ylab("Species pair")
  

so1 <- ggMarginal(so1, size = 7, margins = "x", type = "histogram",
                  col = "black", fill = "black", alpha = 0.75)
so1


###---FOR MEDIAN DIFF ---------------------------------------

## Subset significant regressions so we can highlight in plots
med_over_sig <- subset(med_over, subset = med_over$pval <= 0.05)

## reverse levels of factor spsp so it appears in abc order
med_over$spsp = with(med_over, factor(spsp, levels = rev(levels(spsp))))
med_over$pond <- as.factor(med_over$pond)

## To export to ms, save as image 500 x 750
mo <- ggplot(med_over, aes(x = coef, y = spsp)) + mytheme +
  # all data
  geom_point(size = 3, alpha = 0.5, aes(color = pond)) + 
  # averages
  geom_point(shape = 18, alpha = 1, size = 4,
             data = med_over_err, aes(x = coef, y = spsp)) +
  # error bars on averages
  geom_errorbarh(data = med_over_err, size = 1,
                 aes(xmin = coef - se, xmax = coef + se, height = 0)) +
#     # significant regressions
#     geom_point(shape = 17, size = 6, colour = "red",
#               data = med_over_sig, aes(x = coef, y = spsp)) +
#     # pond as RE
#     geom_point(size = 6, colour = "blue", 
#                data = med_over_re, aes(x = coef, y = spsp)) +
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono", size = 12),
        legend.position = c(0.9, 0.75),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 15)) +
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) + xlim(-1.75, 1.2)+
  xlab("Regression slope \n (numerical overlap ~ difference in median date)") + ylab("Species pair")

mo <- ggMarginal(mo, size = 7, margins = "x", type = "histogram",
                 col = "black", fill = "black", alpha = 0.75)
mo

med_over_err$spsp = factor(med_over_err$spsp, levels = c("ACr : BWo", "ACr : HVe", "RCa : HVe", "HCi : HVe", "GCa : HVe", "HCi : ACr", 
                                                             "HVe : BWo", "ACr : RCl", "RPa : RSp", "GCa : ACr", "HVe : GCa", "BVa : ACr", 
                                                             "RCa : BWo", "RCl : HVe", "BVa : GCa", "GCa : RCl", "RCa : RCl", "BWo : RSp", 
                                                             "GCa : BVa", "BVa : BWo", "HCi : RCl", "BVa : RSp", "HVe : RSp", "RCl : RSp", 
                                                             "PCr : RSp", "BWo : HVe", "RSp : RCl", "GCa : BWo", "BVa : RCl", "HVe : BVa",
                                                             "HCi : BWo", "RSp : BWo", "HVe : ACr", "RSp : HVe", "HVe : HCi", "HVe : RCl", 
                                                             "BWo : HCi", "BWo : GCa", "GCa : HCi", "BWo : RCl", "HCi : GCa", "BWo : BVa", 
                                                             "RSp : PCr", "HCi : BVa", "BVa : HVe", "RCl : ACr", "BVa : HCi", "BWo : ACr"
))
med_over_err$spsp = with(med_over_err, factor(spsp, levels = rev(levels(spsp))))

mo1 <- ggplot(med_over_err, aes(x = coef, y = spsp)) + mytheme +
  # averages
  geom_point(shape = 18, alpha = 1, size = 4) +
  # error bars on averages
  geom_errorbarh(data = med_over_err, size = 1,
                 aes(xmin = coef - se, xmax = coef + se, height = 0)) +
  # all data
  geom_point(size = 3, alpha = 0.5, 
             data = start_over, aes(color = pond)) + 
  # style elements
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono", size = 12),
        legend.position = c(0.9, 0.75),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 15)) +
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) + xlim(-1.75, 1.2)+
  xlab("Regression slope \n (numerical overlap ~ difference in median date)") + ylab("Species pair")

mo1 <- ggMarginal(mo1, size = 7, margins = "x", type = "histogram",
                 col = "black", fill = "black", alpha = 0.75)
mo1
so1
###---START & MEDIAN DIFF AVERAGES TOGETHER -----------------

smo <- ggplot(med_over, aes(x = coef, y = spsp)) + mytheme +
  # averages for median
  geom_point(shape = 18, alpha = 0.75, fill = "black", size = 7,
             data = med_over_err, aes(x = coef, y = spsp)) +
  # error bars for median
  geom_errorbarh(data = med_over_err, size = 1, alpha = 0.75,
                 aes(xmin = coef - se, xmax = coef + se, height = 0)) +
  # averages for start
  geom_point(alpha = 0.75, colour = "darkgray", size = 7,
             data = start_over_err, aes(x = coef, y = spsp)) +
  # error bars for start
  geom_errorbarh(data = start_over_err, size = 1, colour = "darkgray", alpha = 0.75,
                 aes(xmin = coef - se, xmax = coef + se, height = 0)) +
  # style elements
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono"),
        axis.title.x = element_text(size = 17)) +
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  xlab("Regression coefficient\n(numerical overlap ~ difference in single metric)") + ylab("Species pair")
smo


###---R2 COEFFICIENT OF CORRELATION SUMMARY PLOTS -----------

## Start ~ overlap
rs <- ggplot(start_over, aes(x = r2, y = spsp)) + mytheme + 
  # all data
  geom_point(size = 5, alpha = 0.25) + 
#   # averages
#   geom_point(shape = 18, alpha = 1, fill = "black", size = 8,
#              data = regmater, aes(x = r2, y = spsp)) +
#   # error bars on averages
#   geom_errorbarh(data = regmater, size = 1, aes(xmin = r2 - se, xmax = r2 + se, height = 0)) +
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono"),
        axis.title.x = element_text(size = 17)) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  xlab("R-squared correlation coefficient (overlap ~ startdiff)") + ylab("Species pair")

ggMarginal(rs, size = 7, margins = "x", type = "histogram",
           col = "black", fill = "black", alpha = 0.5)

## Median ~ overlap
rm <- ggplot(med_over, aes(x = r2, y = spsp)) + mytheme + 
  # all data
  geom_point(size = 5, alpha = 0.25) + 
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono"),
        axis.title.x = element_text(size = 17)) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  xlab("R-squared correlation coefficient (overlap ~ startdiff)") + ylab("Species pair")

ggMarginal(rm, size = 7, margins = "x", type = "histogram",
           col = "black", fill = "black", alpha = 0.5)




#############################################################
###------------PHENOLOGY ~ YEAR SCATTER PLOTS-------------###
#############################################################

###---POOLING ALL DATA--------------------------------------

## Overlap by year, pooling all data
oy <- ggplot(fall, aes(x = year, y = overlap)) + mytheme +
  geom_point(size = 4, alpha = 0.25) +
  stat_smooth(method = lm, size = 3, se = F, aes(color = factor(spsp))) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
  ylab("Overlap in calling phenology") + xlab("Year")
oy

## Start by year, pooling all data
sy <- ggplot(fall, aes(x = year, y = log(startdiff))) + mytheme +
  geom_point(size = 4, alpha = 0.25, colour = "black") + #aes(color = factor(spsp))) +
  stat_smooth(method = lm, size = 3, se = F, aes(color = factor(spsp))) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
  ylab("Difference in start date") + xlab("Year")
sy

## Overlap by year, isolating well-represented pair
oy_hvrcl <- ggplot(subset(fall, spsp == "HVe : RCl"), 
                 aes(x = year, y = overlap, colour = pond)) + mytheme +
  geom_point(size = 4, alpha = 0.75) +
  stat_smooth(method = lm, size = 3, se = F) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        legend.background = element_rect(size = 1, colour = "black"),
        legend.position = c(0.9, 0.8),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15)) +
  ylab("Overlap in calling phenology") + xlab("Year")
oy_hvrcl

## Start by year, isolating well-represented pair
sy_hvrcl <- ggplot(subset(fall, spsp == "GCa : HVe"), 
                 aes(x = year, y = (startdiff), colour = pond)) + mytheme +
  geom_point(size = 4, alpha = 0.75) +
  stat_smooth(method = lm, size = 1, se = F) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        legend.background = element_rect(size = 1, colour = "black"),
        legend.position = "NULL", #c(0.9, 0.8),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15)) +
  ylab("log(Difference in first call date)") + xlab("Year")
sy_hvrcl


###---ALL DATA, WRAPPED BY FACTORS --------------------------

sy_sp <- ggplot(fall, aes(x = year, y = log(startdiff))) + mytheme +
  geom_point(size = 2, aes(color = factor(pond))) + 
  stat_smooth(method = lm, size = 2, se = T) +
  facet_wrap(~ spsp) +
  theme(legend.position = "none") +
  xlab("Year") + ylab("Difference in Start Date")
#ylim(0, 1)
sy_sp

sy_sp <- ggplot(fall, aes(x = year, y = overlap)) + mytheme +
  geom_point(size = 2, aes(color = factor(pond))) + 
  stat_smooth(method = lm, size = 2, se = T) +
  facet_wrap(~ spsp) +
  theme(legend.position = "none") +
  xlab("Year") + ylab("Overlap")
#ylim(0, 1)
sy_sp







#############################################################
###------------PHENOLOGY ~ YEAR SUMMARY PLOTS-------------###
#############################################################

###---FOR OVERLAP -------------------------------------------

over_year$spsp = with(over_year, factor(spsp, levels = rev(levels(spsp))))

## Subset significant regressions so we can highlight in plots
over_year_sig <- subset(over_year_re, subset = over_year_re$pval <= 0.05)
over_year$pond <- as.factor(over_year$pond)

## Overlap trend over years for each spsp-p, with RE model
oy <- ggplot(over_year, aes(x = coef, y = spsp)) + mytheme +
  ## Ponds shown independently
  geom_point(size = 3, alpha = 0.5, aes(color = pond)) + 
  ## Pond as RE model
  #   geom_point(size = 9, colour = "black", shape = 18,
  #              data = over_year_re, aes(x = coef, y = spsp)) +
  ## Ponds averaged
  geom_point(size = 4, colour = "black", shape = 18,
             data = over_year_err, aes(x = coef, y = spsp)) + 
  ## SE bars on averages
  geom_errorbarh(data = over_year_err, size = 1, aes(xmin = coef - se, xmax = coef + se, height = 0)) +
  ## Significant regressions
  #geom_point(shape = 8, size = 6,
  #          data = over_year_sig, aes(x = coef, y = spsp)) +
  ## Style elements
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono", size = 12),
        axis.text.x = element_text(size = 15),
        legend.position = c(0.9, 0.75),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 15)) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  xlab("Regression slope \n(numerical overlap ~ year)") + ylab("Species pair")

oy <- ggMarginal(oy, size = 7, margins = "x", type = "histogram",
                 col = "black", fill = "black", alpha = 0.75)
oy

## ORDERED BY MAGNITUDE
oy1 <- ggplot(over_year_err, aes(x = coef, y = reorder(spsp, coef))) + mytheme +
  geom_point(size = 4, colour = "black", shape = 18) + 
  ## SE bars on averages
  geom_errorbarh(data = over_year_err, size = 1, aes(xmin = coef - se, xmax = coef + se, height = 0)) +
  ## Ponds shown independently
  geom_point(data = over_year, aes(x = coef, y = spsp, color = pond), 
             size = 3, alpha = 0.5) + 
  ## Style elements
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono", size = 10),
        axis.text.x = element_text(size = 15),
        legend.position = "none",
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 15)) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  xlab("Regression slope \n(numerical overlap ~ year)") + ylab("Species pair")

oy1 <- ggMarginal(oy1, size = 7, margins = "x", type = "histogram",
                 col = "black", fill = "black", alpha = 0.75)
oy1

## Bufo woodhouseii
bwo <- subset(daily, subset = (sp == "BW" & doy > 364))
plot(bwo$cumsum ~ bwo$year)
ggplot(bwo, aes(x = year, y = cumsum, color = factor(pond))) +
  geom_point() + mytheme 

###---FOR START DIFFERENCE ----------------------------------

## Put spsp in alpha order
start_year$spsp = with(start_year, factor(spsp, levels = rev(levels(spsp))))

## Subset significant regressions so we can highlight in plots
start_year_sig <- subset(start_year_re, subset = start_year_re$pval <= 0.05)
start_year$pond <- as.factor(start_year$pond)

## Overlap trend over years for each spsp-p
## To export to ms, save as image with dimensions 500 x 750
sy <- ggplot(start_year, aes(x = coef_log, y = spsp)) + mytheme +
  ## Ponds shown independently
  geom_point(size = 3, alpha = 0.5, aes(color = pond)) + 
  #     ## Pond as RE model
  #     geom_point(size = 9, colour = "black", shape = 18,
  #                data = start_year_re, aes(x = coef_log, y = spsp)) +
  ## Ponds averaged
  geom_point(size = 4, colour = "black", shape = 18,
             data = start_year_err, aes(x = coef_log, y = spsp)) + 
  ## SE bars on averages
  geom_errorbarh(data = start_year_err, size = 1, aes(xmin = coef_log - se, xmax = coef_log + se, height = 0)) +
  ## Significant regressions
  #geom_point(shape = 8, size = 6,
  #            data = start_year_sig, aes(x = coef, y = spsp)) +
  ## Style elements
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono", size = 10),
        axis.text.x = element_text(size = 15),
        legend.position = "NULL",
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 15)) +
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  xlab("Regression slope \n(log(difference in start date) ~ year)") + ylab("Species pair") + xlim(-1, 1) 

sy <- ggMarginal(sy, size = 7, margins = "x", type = "histogram",
                 col = "black", fill = "black", alpha = 0.75)
sy

## Put spsp in the order of effect size for overlap trend
start_year_err$spsp = factor(start_year_err$spsp, levels = c("GCa : HCi", "GCa : BVa", "BVa : RCl", "RCa : BWo", "BVa : HCi", "BVa : GCa",
                                  "GCa : BWo", "BVa : RSp", "GCa : RCl", "ACr : BWo", "HCi : BWo", "HCi : GCa",
                                  "RCl : ACr", "RCa : HVe", "ACr : HVe", "BWo : HCi", "BWo : GCa", "GCa : ACr",
                                  "RSp : RCl", "BWo : BVa", "HVe : BWo", "BVa : ACr", "BWo : RCl", "HVe : GCa",
                                  "HVe : RCl", "HCi : BVa", "BWo : HVe", "RSp : HVe", "BWo : ACr", "HVe : BVa",
                                  "HVe : HCi", "RSp : BWo", "BVa : HVe", "BWo : RSp", "RCl : RSp", "HVe : RSp", 
                                  "RCa : RCl", "PCr : RSp", "RCl : HVe", "BVa : BWo", "ACr : RCl", "HCi : ACr", 
                                  "RPa : RSp", "RSp : PCr", "GCa : HVe", "HCi : HVe", "HCi : RCl", "HVe : ACr"))
start_year_err$spsp = with(start_year_err, factor(spsp, levels = rev(levels(spsp))))
start_year$pond <- as.factor(start_year$pond)

sy1 <- ggplot(start_year_err, aes(x = coef_log, y = spsp)) + mytheme +
  geom_point(size = 4, colour = "black", shape = 18) + 
  ## SE bars on averages
  geom_errorbarh(data = start_year_err, size = 1, aes(xmin = coef_log - se, xmax = coef_log + se, height = 0)) +
  ## Ponds shown independently
  geom_point(data = start_year, aes(x = coef_log, y = spsp, color = pond), 
             size = 3, alpha = 0.5) + 
  ## Style elements
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono", size = 10),
        axis.text.x = element_text(size = 15),
        legend.position = c(0.9, 0.75),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 15)) +
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  xlab("Regression slope \n(log(difference in start date) ~ year)") + ylab("Species pair") + xlim(-1, 1)
sy1 <- ggMarginal(sy1, size = 7, margins = "x", type = "histogram",
                 col = "black", fill = "black", alpha = 0.75)
sy1
oy1 # these aren't in the same order, so will have to manually put sy1 in the right order... :(

###---R2 COEFFICIENT OF CORRELATION SUMMARY PLOTS -----------

## Overlap ~ year
oyr <- ggplot(over_year, aes(x = r2, y = spsp)) + mytheme +
  geom_point(size = 5, alpha = 0.25, color = "black") + 
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono")) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  xlab("Regression coefficient (overlap ~ year)") + ylab("Species pair") 

ggMarginal(oyr, size = 7, margins = "x", type = "histogram",
           col = "black", fill = "black", alpha = 0.75)

## Start ~ year
syr <- ggplot(start_year, aes(x = r2, y = spsp)) + mytheme +
  geom_point(size = 5, alpha = 0.25, color = "black") + 
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono")) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  xlab("Regression coefficient (Startdiff ~ year)") + ylab("Species pair") 

ggMarginal(syr, size = 7, margins = "x", type = "histogram",
           col = "black", fill = "black", alpha = 0.75)





#############################################################
###---------------OTHER COOL SUMMARY PLOTS----------------###
#############################################################

###---COMMUNITY LEVEL PHENOLOGY GRAMS---------------------------------------------------

## Eliminate 0 observations to utilize geom_density
d <- subset(daily, select = c("pond", "year", "sp", "dailysum", "doy"))
d <- subset(d, subset = (dailysum != 0))

## Add a column with full species names
d$species <- factor(d$sp,
                    levels = c("AC", "BV", "BW", "GC", "HC", "HV", "PC", "PT", "RCA", "RCL", "RP", "RS"),
                    labels = c("A. crepitans", "B. valliceps", "B. woodhouseii", "G. carolinensis", "H. cineria",
                               "H. versicolor", "P. crucifer", "P. triseriata", "R. catesbeiana", "R. clamatans",
                               "R. palustris", "R. sphenocephala"))

## Stacked ~ sp
cp1 <- ggplot(d, aes(doy)) + mytheme +
  geom_density(size = 1.25, color = "darkslategray", fill = "darkslategray", alpha = 0.75, adjust = 1/2) + 
  facet_grid(species ~., scales = "free_y", switch = "y") + 
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        strip.text.y = element_text(size = 15, angle = 180, face = "italic")) +
  ylab(NULL) + xlab("Day of the year") 
cp1

## Lava plot (bubble)
cpbw <- ggplot(subset(d,sp == "BW"), aes(doy)) + mytheme +
  stat_density(size = 1.25, alpha = 0.75, adjust = 1/5,
               aes(ymax = ..density..,  ymin = -..density..),
               fill = "grey37", colour = "grey27",
               geom = "ribbon", position = "identity") +
  facet_grid(year ~.  , switch = "y") +
  scale_x_continuous(breaks = seq(0, 350, 50)) +
  theme(axis.ticks   = element_blank(), 
        axis.text.y  = element_blank(),
        axis.line.x  = element_line(size = 1),
        strip.text.y = element_text(size = 15, angle = 180, face = "italic")) +
  ylab(NULL) + xlab("Day of the year")
cpbw

cp2 <- ggplot(d, aes(doy)) + mytheme +
  stat_density(size = 1.5, alpha = 0.75, adjust = 1/5,
               aes(ymax = ..density..,  ymin = -..density..),
               fill = "darkslategray", colour = "darkslategray",
               geom = "ribbon", position = "identity") +
  facet_grid(species ~ ., switch = "y", scales = 'free') +
  scale_x_continuous(breaks = seq(0, 350, 50)) +
  theme(axis.ticks   = element_blank(), 
        axis.text.y  = element_blank(),
        axis.line.x  = element_line(size = 1),
        axis.text.x  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        strip.text.y = element_text(size = 12, angle = 180, face = "italic")) +
  ylab(NULL) + xlab("Day of year") + xlim(0, 365)
cp2

## whole community pooled-- visual trend for convergence
wholecomm <- ggplot(d, aes(doy)) + mytheme +
  stat_density(size = 1.5, alpha = 0.75, adjust = 1/5,
               aes(ymax = ..density..,  ymin = -..density..),
               fill = "darkslategray", colour = "darkslategray",
               geom = "ribbon", position = "identity") +
  facet_grid(year ~.  , switch = "y") +
  scale_x_continuous(breaks = seq(0, 350, 50)) +
  theme(axis.ticks   = element_blank(), 
        axis.text.y  = element_blank(),
        axis.line.x  = element_line(size = 1),
        axis.text.x  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        strip.text.y = element_text(size = 14, angle = 180, face = "italic")) +
  ylab(NULL) + xlab("Day of year") + xlim(0, 365)
wholecomm


sf <- subset(d, subset = (pond < 5))
dc <- subset(d, subset = (pond > 4))

cpsites <- ggplot(dc, aes(doy)) + mytheme +
  stat_density(size = 1.25, alpha = 0.75, adjust = 1/5,
               aes(ymax = ..density..,  ymin = 0),
               fill = "tomato1", colour = "tomato2",
               geom = "ribbon", position = "identity") +
  facet_grid(species ~. , switch = "y", scales = 'free') +
  scale_x_continuous(breaks = seq(0, 350, 50)) +
  theme(axis.ticks   = element_blank(), 
        axis.text.y  = element_blank(),
        axis.line.x  = element_line(size = 1),
        strip.text.y = element_text(size = 15, angle = 180, face = "italic")) +
  ylab(NULL) + xlab("Day of the year") +
  stat_density(data = sf, size = 1.25, alpha = 0.75, adjust = 1/5,
               aes(doy, ymax = 0, ymin = -..density..),
               fill = "navyblue", colour = "navyblue",
               geom = "ribbon", position = "identity")
cpsites

early <- subset(d, subset = (year < 2006))
late  <- subset(d, subset = (year > 2010))

cptime <- ggplot(early, aes(doy)) + mytheme +
  stat_density(size = 1.25, alpha = 0.75, adjust = 1/5,
               aes(ymax = ..density..,  ymin = 0),
               fill = "tomato1", colour = "tomato2",
               geom = "ribbon", position = "identity") +
  facet_grid(species ~. , switch = "y", scales = 'free') +
  scale_x_continuous(breaks = seq(0, 350, 50)) +
  theme(axis.ticks   = element_blank(), 
        axis.text.y  = element_blank(),
        axis.line.x  = element_line(size = 1),
        strip.text.y = element_text(size = 15, angle = 180, face = "italic")) +
  ylab(NULL) + xlab("Day of the year") +
  stat_density(data = late, size = 1.25, alpha = 0.75, adjust = 1/5,
               aes(doy, ymax = 0, ymin = -..density..),
               fill = "navyblue", colour = "navyblue",
               geom = "ribbon", position = "identity")
cptime


###---OVERLAP HEAT MAP------------------------------------------------------

## Make color palette
col <- colorRampPalette(brewer.pal(n = 9, "BuPu"))
psz <- 256
pal <- col(psz)

## Heat map of overlap for sp pairs, wrapped by year, averaged by pond
oy <- ggplot(all, aes(x = sp1, y = sp2, fill = overlap)) +
  geom_tile() +
  facet_wrap(~ year) +
  scale_fill_gradient2(low  = pal[1],
                       mid  = pal[psz/2],
                       high = pal[psz], 
                       midpoint = (max(all$overlap) +
                                     min(all$overlap))/2,
                       name = "Overlap") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = rel(1.5))) +
  xlab("Focal Species") + ylab("Secondary Species")
oy


###---BOXPLOTS SHOWING VARIATION IN DISTRIBUTIONS------------

## Remove non-calling years and write full sp name
p <- na.omit(periods)
p <- subset(p, subset = (duration > 0))
p$species <- factor(p$sp,
                    levels = c("AC", "BV", "BW", "GC", "HC", "HV", "PC", "PT", "RCA", "RCL", "RP", "RS"),
                    labels = c("A. crepitans", "B. valliceps", "B. woodhouseii", "G. carolinensis", "H. cineria",
                               "H. versicolor", "P. crucifer", "P. triseriata", "R. catesbeiana", "R. clamatans",
                               "R. palustris", "R. sphenocephala"))
p$species <- factor(p$species, levels=rev(levels(p$species)))

## Variation in duration
duration <- ggplot(p, aes(x = species, y = duration)) + mytheme +
  geom_boxplot(fill = "darkgreen", alpha = 0.5, width = 0.5) +
  coord_flip() +
  theme(axis.ticks   = element_blank(), 
        axis.line.x  = element_line(size = 1),
        axis.text.y  = element_text(size = 12, face = "italic"),
        axis.text.x  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  ylab("Duration of Calling Period") + xlab("Species")
duration

## Variation in median date (very similar to lava plot)
median <- ggplot(p, aes(x = species, y = meddate)) + mytheme +
  geom_boxplot(fill = "darkslategray", alpha = 0.5, width = 0.5) +
  geom_boxplot(data = p, aes(x = species, y = firstdate),
               color = 'darkgreen', fill = 'darkgreen', alpha = 0.5, width = 0.5) +
  geom_boxplot(aes(x = species, y = lastdate),
               color = 'darkred', fill = 'darkred', alpha = 0.5, width = 0.5, position = 'dodge') +
  coord_flip()+
  theme(axis.ticks   = element_blank(), 
        axis.line.x  = element_line(size = 1),
        axis.text.y  = element_text(size = 15, face = "italic")) +
  ylab("Median Date of Calling") + xlab("Species")
median

## Formatted for poster
bp <- ggplot(p, aes(x = species, y = duration)) + mytheme +
  geom_boxplot(fill = "darkgreen", alpha = 0.5, width = 0.5) +
  theme(axis.ticks   = element_blank(), 
        axis.line.x  = element_line(size = 1),
        axis.text.x  = element_text(size = 15, face = "italic", angle = 60, hjust = 1),
        axis.text.y  = element_text(size = 18),
        axis.title.y  = element_text(size = 18)) +
  ylab("Duration of calling period (days)") + xlab(NULL)
bp

bpm <- ggplot(p, aes(x = species, y = meddate)) + mytheme +
  geom_boxplot(fill = "darkgreen", alpha = 0.5, width = 0.5) +
  theme(axis.ticks   = element_blank(), 
        axis.line.x  = element_line(size = 1),
        axis.text.x  = element_text(size = 15, face = "italic", angle = 60, hjust = 1),
        axis.text.y  = element_text(size = 18),
        axis.title.y  = element_text(size = 18)) +
  ylab("Median call date (DOY)") + xlab(NULL)
bpm


#############################################################
###---------------------WEATHER PLOTS---------------------###
#############################################################

###---RAIN-----------------------------------------

## Raw daily data
rain <- ggplot(weather, aes(x = doy, y = rain, color = site)) + mytheme +
  geom_point(size = 2, alpha = 0.5) +
  # stat_smooth(size = 2, alpha = 0.75, se = F) +
  facet_wrap(~ year)
rain

## Cumulative rain
cumrain <- ggplot(subset(weather, year != 2003 & year != 2004 & year != 2005 & year != 2015), 
                  aes(x = doy, y = cumrain, color = site)) + mytheme +
  geom_point(size = 2, alpha = 0.15) +
  stat_smooth(size = 2, alpha = 0.5, se = T) +
  facet_wrap(~ year) + xlab("Day of year") + ylab("Cumulative Rainfall (cm)") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
cumrain

## Species rain response
ro_sp <- ggplot(weather_over, aes(x = cumrain, y = avgover, color = site)) +
  stat_smooth(method = lm, size = 2, se = T) +
  geom_point(size = 2) + 
  facet_wrap(~ spsp) +
  xlab("cumulative yearly rainfall (cm)") + ylab("overlap") +
  ylim(0, 1)
ro_sp

###---TEMPERATURE--------------------------------------

## Raw daily temp
temp <- ggplot(weather, aes(x = doy, y = temp, color = site)) + mytheme +
  geom_point(size = 2, alpha = 0.5) +
  stat_smooth(se = F, size = 2) +
  facet_wrap(~ year)
temp

## Annual average temp
avgtemp <- ggplot(subset(weather_sum, year < 2015 & year > 2005),
                  aes(x = (year), y = avgtemp, color = site)) + mytheme +
  #geom_bar(stat = "identity", position = "dodge", aes (fill = site)) +
  ylab("Average temperature (C) at 9:00pm") + xlab("Year") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  geom_point(size = 5) + geom_line(size = 3, alpha = 0.75) #+ 
  #stat_smooth(method = lm, size = 3)
avgtemp

## Species temperature response
to_sp <- ggplot(weather_over, aes(x = avgtemp, y = avgover, color = site)) +
  stat_smooth(method = lm, size = 2, se = T) +
  geom_point(size = 2) + 
  facet_wrap(~ spsp) +
  xlab("Average temperature") + ylab("overlap") +
  ylim(0, 1)
to_sp

###---SPECIES OVERLAP BY RAIN--------------------------------

## Subset significant regressions so we can highlight in plots
rain_over_sig <- subset(weather_over_coef, subset = (weather_over_coef$pval_rain <= 0.05))

## Overlap trend over years for each spsp-p, with RE model
ro <- ggplot(weather_over_coef, aes(x = coef_rain, y = spsp)) + mytheme +
  ## Ponds shown independently
  geom_point(size = 3, alpha = 0.5, aes(color = site)) + 
  ## Significant regressions
  geom_point(shape = 8, size = 6,
           data = rain_over_sig, aes(x = coef_rain, y = spsp)) +
  ## Style elements
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono", size = 10),
        legend.position = c(0.9, 0.8)) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  xlab("Regression slope \n(Numerical Overlap ~ Cumulative Yearly Rain)") + ylab("Species pair") 
  #xlim(-0.05, 0.05)

ro <- ggMarginal(ro, size = 7, margins = "x", type = "histogram",
                 col = "black", fill = "black", alpha = 0.75)
ro

###---SPECIES OVERLAP BY TEMPERATURE---------------------------

## Subset significant regressions so we can highlight in plots
temp_over_sig <- subset(weather_over_coef, subset = (weather_over_coef$pval_temp <= 0.05))

## Overlap trend over years for each spsp-p, with RE model
to <- ggplot(weather_over_coef, aes(x = coef_temp, y = spsp)) + mytheme +
  ## Ponds shown independently
  geom_point(size = 3, alpha = 0.5, aes(color = site)) + 
  ## Significant regressions
  geom_point(shape = 8, size = 6,
             data = temp_over_sig, aes(x = coef_temp, y = spsp)) +
  ## Style elements
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono", size = 10),
        legend.position = c(0.9, 0.8)) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  xlab("Regression slope \n(Numerical Overlap ~ Average Annual Temperature)") + ylab("Species pair") #+ xlim(-0.5, 0.5)

to <- ggMarginal(to, size = 7, margins = "x", type = "histogram",
                 col = "black", fill = "black", alpha = 0.75)
to



#############################################################
###-------------INDIVIDUAL SPECIES TRENDS-----------------###
#############################################################

###---START DATE TREND--------------------------------------

## All species
ggplot(periods, aes(x = year, y = firstdate)) +
  geom_point(size = 4, alpha = 0.5, aes(color = factor(pond))) +
  stat_smooth(method = lm, size = 2, se = T) + #,
              #aes(color = factor(pond))) +
  facet_wrap(~ sp) + mytheme

## GCa is only significant shifter
gc1 <- subset(periods, subset = (sp == 'GC'))
ggplot(gc1, aes(x = year, y = firstdate)) +
  geom_point(size = 4, aes(color = factor(pond))) +
  stat_smooth(method = lm, size = 2, se = T) #+ mytheme

###---MODEL FOR SPECIES TRENDS IN FIRST AND MEDIAN------------

## could not get function to work so go through sp manually 
f  <- summary(lmer(firstdate ~ year + (1|pond), data = s))
m  <- summary(lmer(meddate ~ year + (1|pond), data = s))



#############################################################
###-----------------FOR REVIEWER 2------------------------###
#############################################################

###---SOME EXTRA STUFF---------------------------------------
ggplot(all, aes(x = sp1, y = overlap)) + 
  geom_boxplot() +
  geom_point(alpha = 0.5, position = 'jitter', aes(color = sp2)) + 
  ylim(0.01, 0.99) + mytheme
ggplot(all, aes(x = sp2, y = overlapf2)) + 
  geom_boxplot() + 
  geom_point(alpha = 0.5, position = 'jitter', aes(color = sp1)) + 
  ylim(0.01, 0.99) +mytheme
ggplot(all, aes(x = ))
ggplot(all, aes(x = overlap, y = overlapf2)) +
  geom_point() + mytheme

ggplot(daily, aes(x = sp, y = cumsum, color = year)) +
  geom_point()

x <- ggplot(daily, aes(x = sp, y = cumsum, color = as.factor(year))) + geom_bar()
x

## end of year totals
sum <- subset(daily, subset = (doy == 364))
sum$pond <- as.factor(sum$pond)
x <- ggplot (sum, aes(x = year, y = cumsum)) + 
  geom_point(aes(color = pond)) + 
  #geom_smooth(size = 1.5, se = T) + 
  geom_smooth(size = 1.5, se = F) +
  facet_wrap(~ sp, scales = 'free_y') +
  ylab("abundance") 
x


v <- ggplot(sum, aes(x = reorder(sp, cumsum), y = cumsum)) +
  geom_boxplot() + xlab("species") + ylab("total abundance") + 
  ylim(0, 2000) #+ mytheme
v

## hard to interpret this because there's lots of missing data in earlier years
## so that's probably what's driving most minimal increases, with maybe HV actually increasing

m <- lmer(cumsum ~ year + (1 | pond), data = subset(sum, subset = (sp == "GC")))

## total number of calls, across years and sites


totals <- sum %>%
  group_by(sp) %>%
  summarize(totalcalls = sum(cumsum))
View(totals)  
  
y <- ggplot(totals, aes(x = sp, y = totalcalls)) +
  geom_point()
y



#############################################################
###-----------------FIGURES FOR PAPER---------------------###
#############################################################

###---FIGURE 3: TRENDS OVER TIME-----------------------------
## A. OVERLAP ~ YEAR
oy2 <- ggplot(over_year_err, aes(x = coef, y = reorder(spsp, coef))) + mytheme +
  geom_point(size = 4, colour = "black", shape = 18) + 
  ## SE bars on averages
  geom_errorbarh(data = over_year_err, size = 1, aes(xmin = coef - se, xmax = coef + se, height = 0)) +
  ## Ponds shown independently
  geom_point(data = over_year, aes(x = coef, y = spsp, color = pond), 
             size = 3, alpha = 0.5) + 
  ## Style elements
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono", size = 10),
        axis.text.x = element_text(size = 15),
        legend.position = "none",
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 15)) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  xlab("Regression slope \n(numerical overlap ~ year)") + ylab("Species pair")

oy2 <- ggMarginal(oy2, data = over_year, x = coef, size = 7, 
                  type = 'histogram', margins = "x", 
                  col = "black", fill = "black", alpha = 0.75)
oy2

## B. STARTDATE ~ YEAR
sy2 <- ggplot(start_year_err, aes(x = coef_log, y = spsp)) + mytheme +
  geom_point(size = 4, colour = "black", shape = 18) + 
  ## SE bars on averages
  geom_errorbarh(data = start_year_err, size = 1, aes(xmin = coef_log - se, xmax = coef_log + se, height = 0)) +
  ## Ponds shown independently
  geom_point(data = start_year, aes(x = coef_log, y = spsp, color = pond), 
             size = 3, alpha = 0.5) + 
  ## Style elements
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15),
        legend.position = c(0.85, 0.75),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_blank()) +
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) +
  xlab("Regression slope \n(log(difference in start date) ~ year)") + ylab("Species pair") + xlim(-1, 1)
sy2 <- ggMarginal(sy2, data = start_year, x = coef_log, size = 7, margins = "x", type = "histogram",
                  col = "black", fill = "black", alpha = 0.75)
sy2

## TOGETHER
fig3 <- plot_grid(oy2, sy2, rel_widths = c(1, 0.85), labels = "AUTO")
tiff("figure3.tiff", height = 19, width = 21, units = 'cm', res = 600)
plot(fig3)
dev.off()

###---FIGURE 4: BOXPLOT--------------------------------------
duration <- ggplot(p, aes(x = species, y = duration)) + mytheme +
  geom_boxplot(fill = "darkgreen", alpha = 0.5, width = 0.5) +
  coord_flip() +
  theme(axis.ticks   = element_blank(), 
        axis.line.x  = element_line(size = 1),
        axis.text.y  = element_text(size = 12, face = "italic"),
        axis.text.x  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  ylab("Duration of Calling Period") + xlab("Species")
duration

tiff("figure4.tiff", height = 13, width = 14, units = 'cm', res = 600)
plot(duration)
dev.off()

###---FIGURE 5: SEASONAL CALLING PATTERNS--------------------
cp2 <- ggplot(d, aes(doy)) + mytheme +
  stat_density(size = 1.5, alpha = 0.75, adjust = 1/5,
               aes(ymax = ..density..,  ymin = -..density..),
               fill = "darkslategray", colour = "darkslategray",
               geom = "ribbon", position = "identity") +
  facet_grid(species ~ ., switch = "y", scales = 'free') +
  scale_x_continuous(breaks = seq(0, 350, 50)) +
  theme(axis.ticks   = element_blank(), 
        axis.text.y  = element_blank(),
        axis.line.x  = element_line(size = 1),
        axis.text.x  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        strip.text.y = element_text(size = 12, angle = 180, face = "italic")) +
  ylab(NULL) + xlab("Day of year") + xlim(0, 365)
cp2

tiff("figure5.tiff", height = 17, width = 18, units = 'cm', res = 600)
plot(cp2)
dev.off()
###---FIGURE 6: OVERLAP ~ METRICS----------------------------
## A. OVERLAP ~ START
so2 <- ggplot(start_over_err, aes(x = coef, y = reorder(spsp, coef))) + mytheme +
  geom_point(size = 4, colour = "black", shape = 18) + 
  ## SE bars on averages
  geom_errorbarh(data = start_over_err, size = 1, aes(xmin = coef - se, xmax = coef + se, height = 0)) +
  ## Ponds shown independently
  geom_point(data = start_over, aes(x = coef, y = spsp, color = pond), 
             size = 3, alpha = 0.5) + 
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text("mono", size = 12),
        legend.position = "NULL", #c(0.9, 0.7),
        axis.title.x = element_text(size = 11.5),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 14)) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) + xlim(-1.75, 1.2)+
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  xlab("Regression slope \n (numerical overlap ~ difference in start date)") + ylab("Species pair")

so2 <- ggMarginal(so2, size = 7, margins = "x", type = "histogram",
                  col = "black", fill = "black", alpha = 0.75)
so2

## B. OVERLAP ~ MEDIAN
mo2 <- ggplot(med_over_err, aes(x = coef, y = spsp)) + mytheme +
  # averages
  geom_point(shape = 18, alpha = 1, size = 4) +
  # error bars on averages
  geom_errorbarh(data = med_over_err, size = 1,
                 aes(xmin = coef - se, xmax = coef + se, height = 0)) +
  # all data
  geom_point(size = 3, alpha = 0.5, 
             data = start_over, aes(color = pond)) + 
  # style elements
  theme(legend.background = element_rect(size = 1, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_blank(),
        legend.position = c(0.85, 0.75),
        axis.title.x = element_text(size = 11.5),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14)) +
  geom_vline(size = 1.5, linetype = 2,  xintercept = 0) +
  geom_hline(yintercept = seq(0.5, 48.5, 1), color = "lightgray", linetype = 2) + xlim(-1.75, 1.5)+
  xlab("Regression slope \n (numerical overlap ~ difference in median date)") + ylab("Species pair")

mo2 <- ggMarginal(mo2, size = 7, margins = "x", type = "histogram",
                  col = "black", fill = "black", alpha = 0.75)
#mo2

## TOGETHER
fig6 <- plot_grid(so2, mo2, rel_widths = c(1, 0.85), labels = "AUTO")
fig6
tiff("figure6.tiff", height = 19, width = 21, units = 'cm', res = 600)
plot(fig6)
dev.off()
