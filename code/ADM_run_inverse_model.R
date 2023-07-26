
#### Setup ####

# This code should be executed inside the R project that is provided in the SOM alongside
# this script to ensure that relative file paths work correctly

# load libraries
library("dplyr")
library("Rfast")
library("ggplot2")
library("progress")
library("matrixStats")
library("tictoc")
library("beepr")


# load inverse model functions

source("code/ADM_inverse_functions.R")

# specify a project name for printing on  model reports

projectname <- "Abri du Maras"  

#### read in isotope data ####

# read data input for all teeth
d18O.data <- read.csv("ADM_cleaned_d18O_model_input.csv")

# filter for correct tooth id
tooth.name <- "MARAS_687" # specify the sample name of the tooth

tooth.data <- d18O.data %>%
  filter(Tooth.ID == tooth.name)

#### define some empty objects ####

numtrials <- 100 # number of model solutions
# some empty objects for code to work
Edist <- rep(0,numtrials)
allTrials <- rep(0,numtrials)

#### running the Emeas function ####

## input variables:

# numtrials - number of trials to run iteration
# numsam - number of isotope samples/values
# length - distance between each sample and the next sample in mm*10, occlusal sample to go first
# dMeas - d18O values, occlusal sample to go first
# r1 - (real) isotope 1 sigma reproducibility
# r2 - (integer) sample length 1 sigma reproducibility (in mm * 10)
# r3 - (integer) sample depth 1-sigma reproducibility (in mm * 10)
# la - (integer) length of apposition (in mm * 10) --> get species specific from literature

r1 <- 0.2
r2 <- 1
r3 <- 5
la <- 60

# transform sample length data into 10*mm units for model input

tooth.data <- tooth.data %>%
  mutate(Passey.length = length.raw*10)

## run the Emeas function to generate Emeas estimates

Emeasout <- Emeasfun(numtrials = numtrials, numsam = length(tooth.data$mean.d18O), length = tooth.data$Passey.length, 
         dMeas = round(tooth.data$mean.d18O, 1), r1 = r1, r2 = r2, r3 = r3, la = la)

# write input parameters into a list to print in model reports
Emeas.params <- list(numtrials = numtrials, r1 = r1, r2 = r2, r3 = r3, la = la)

#### plot output ####  
ggplot(Emeasout, aes(x = totallength, y = allTrials))+
  theme_classic()+
  xlim(500,0)+
  geom_line()+
  geom_point(shape = 3)+
  geom_point(aes(x = totallength, y = dMeas), colour = "black", shape = 0)+
  geom_line(aes(x = totallength, y = dMeas), colour = "black")+
  geom_line(aes(x = totallength, y = dMeasError), colour = "black")+
  geom_point(aes(x = totallength, y = dMeasError), colour = "black", shape = 1)

#### evaluate Emeas to adjust df in mSolv code ####
# the distribution of Emeas values is stored in an object called 'Edist'
# values of Emeas/Edist should be close to DPE values from mSolv code
# adjust df stepwise to match Emeas/Edist to DPE

hist(Edist)
mean.edist <- mean(Edist)
mean.edist


#### running the mSolv1_1 function ######

#### input paramaters ####

nsolxns <- 100 # number of solutions to be computed
dMeas <- round(tooth.data$mean.d18O, 1) # isotope data input
numsam <- length(dMeas) # number of samples
openindx <- 1 # degree of openendedness, less than lm, openended (profile mature) --> index = 1; 
# close ended (enamel immature) index = lm
avelength <- round(mean(tooth.data$length.raw)*10, digits = 0) # average sample length
# length of maturation
lm <- 280 #in mm*10
#length of apposition in mm
lamm <- la/10
# enter sample depth as fraction of la in mm*10
# (given in original data as fraction of enamel thickness)
depth <- round(tooth.data$depth*lamm, 1)*10
finit <- 0.25 # initial mineral %weight

# some derived numbers used during modelling
numbefore <- ceiling(la/avelength)
numafter <- ceiling((lm-openindx)/avelength)+1
# empty object for DPE to enable writing to global env from inside function
MEST = matrix(0, nrow = (numsam+numbefore+numafter-1), ncol = nsolxns)
DPE = rep(0,nsolxns)
S = matrix(0, nrow = numsam, ncol = nsolxns)

# input parameters of the reference vector

# max sample length of the reference vector
maxlength = 45
# min sample length of the reference vector
minlength = 10
# min depth of the reference vector
mindepth = 2
# maximum value of the reference vector
maxratio = 16.8
# minimum value of the reference vector
minratio = 17.8
# sd for random draws that produce reference vector
stdev = 0.5
# damping factor
df = 0.01


## after adjusting df, the model needs to be tested for sensitivity to the reference vector
## --> check whether more oscillating (max and min further apart) generates drastically different solution
# df - damping factor. Needs to be chosen to minimize difference between the estimated measurement error
# (E~meas~) and the prediction error (E~pred~)

# write parameter input into a list for printing the model report
mSolvparams <- list(nsolxns = nsolxns, openindx = openindx, lm = lm, la = la, finit = finit, maxlength = maxlength, minlength = minlength, 
                    mindepth = mindepth,r1 = r1, r2 = r2, r3 = r3, maxratio = maxratio, minratio = minratio, stdev = stdev, df = df)

### run mSolv function ###

# this part of the script takes a long time to run, so I have designed it to work with a progress bar
# and a beep sound when it is finished
# everything needs to be selected from tic() until beep() and executed together
tic()
lfsolv <- mSolv_fun(nsolxns = nsolxns, numsam = numsam, finit = finit, la = la, lm = lm, openindx = 1, avelength = avelength, 
                    maxlength = maxlength, minlength = minlength, mindepth = mindepth, length = tooth.data$Passey.length,
                    dMeas = dMeas, r1 = r1, r2 = r2, r3 = r3, maxratio = maxratio, minratio = minratio, 
                    stdev = stdev, df = df, depth = depth)
toc()
beep(sound = 2)



# plot some example solutions of mSolv

ggplot(solvout, aes(x = totallength, y = dMeasd))+
  theme_classic()+
  geom_point(colour = "black")+
  geom_line(colour = "black")+
  geom_line(aes(x = totallength, y = trial1), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial2), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial3), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial4), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial5), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial6), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial7), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial8), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial9), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial10), colour = "lightblue")
  

#### evaluate DPE (= Epred, prediction error) vs Edist (= Emeas, estimated measurement error) ####

# DPE should be close to Edist
# adjust df to achieve this
# increase df to increase DPE

hist(DPE, breaks = 10)
mean(DPE)
# compare to measurement error calculated in the first part of the model
mean.edist

#### generate report of model parameters ####

rmarkdown::render("make_model_report.Rmd", params = list(
  project = projectname, 
  tooth = tooth.name, 
  Emeasparams = Emeas.params, 
  Edist = Edist, 
  mSolvparams = mSolvparams, 
  DPE = DPE
), output_file = paste0(tooth.name, "_inverse_model_report.html"))

#### extract output with mean and CI and write to csv ####

# extract only trial output

tdata <- solvout %>%
  select(-totallength)

# transpose
tdatatrans <- t(tdata)
tdatatrans

nl <- 1:(numsam+numbefore+numafter-1)
lengths <- paste("l",nl,sep="")

colnames(tdatatrans) <- lengths

# calculate mean and confidence interval (95%)

tdata_lci <- apply(tdatatrans,MARGIN=2,quantile,prob=0.025,na.rm=T) # calculate lower boundary of Conf Interval 
tdata_uci <- apply(tdatatrans,MARGIN=2,quantile,prob=0.975,na.rm=T) # calculate upper boundary of Conf Interval 
tdata_mean <- apply(tdatatrans,MARGIN=2,mean) # calculate mean (essentially the most likely model solution)

# bind into data frame and clean up names
tdataci <- cbind(solvout$totallength,tdata_lci,tdata_uci,tdata_mean)
colnames(tdataci) <- c("ci.length","lower.CI","upper.CI","mean")

# convert to data frame
tdataci.d <- as.data.frame(tdataci)

# convert lengths back into mm
tdataci.d$ci.length <- tdataci.d$ci.length/10

# combine all output
all.out <- cbind(solvout, tdataci.d)

all.out <- all.out %>%
  mutate(tooth = tooth.name) %>%
  mutate(Layer = tooth.data$Layer[1])

# write output to csv
outputfile <- paste(tooth.name, "_inverse_model_output.csv", sep = "")
write.csv(all.out, file = outputfile)



















