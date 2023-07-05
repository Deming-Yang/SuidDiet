library(scales)
library(viridisLite)
library(dplyr)
library(tidyr)
library(stats)

setwd("C:/Users/ydmag/Google Drive/Dissertation data/Extant suid seasonal diet/SuidDiet")

extant <- read.csv("data/suid_d13C_d18O.csv") #data from this study

ID.list <- levels(factor(extant$ID))
print(ID.list)

# [1] "NKU245" "NKU257" "NKU265" "P01"    "P03"    "P04"    "P05"    "P06"    "P07"    "SAM01"  "SAM02" 
# [12] "SAM03"  "SBL01"  "SBL02"  "SBL03"

#################data parsing################
P01 <- extant %>% filter(ID == ID.list[4])
P03 <- extant %>% filter(ID == ID.list[5])
P05 <- extant %>% filter(ID == ID.list[7])
NKU245 <- extant %>% filter(ID == ID.list[1])
NKU257 <- extant %>% filter(ID == ID.list[2])
NKU265 <- extant %>% filter(ID == ID.list[3])
P04 <- extant %>% filter(ID == ID.list[6])
P06 <- extant %>% filter(ID == ID.list[8])
SAM01 <- extant %>% filter(ID == ID.list[10])
SAM02 <- extant %>% filter(ID == ID.list[11])
SAM03 <- extant %>% filter(ID == ID.list[12])
SBL01 <- extant %>% filter(ID == ID.list[13])
SBL02 <- extant %>% filter(ID == ID.list[14])
SBL03 <- extant %>% filter(ID == ID.list[15])

################correction for molar-canine d13C spacing###############
P01 <- P01 %>% mutate( X.13C.corr = X.13C - 2)

P03 <- P03 %>% mutate( X.13C.corr = X.13C - 2)

P05 <- P05 %>% mutate( X.13C.corr = X.13C - 2)

################ calculate end member dietary values ################

#Cerling and Harris 1999

d.C4 <- -12.8

d.C3 <- -27.4 

epsilon <- 13.3 #Passy 2005, diet to molar enamel fractionation in domestic pigs

#nominal pure C4 diet, linear mixing
Enamel.C4 <- ((1 + d.C4/1e3)*(1 + epsilon/1e3) - 1) * 1e3 #0.3 per mill

#nominal pure C3 diet, linear mixing
Enamel.C3 <- ((1 + d.C3/1e3)*(1 + epsilon/1e3) - 1) * 1e3 #-14.5 per mill

#mixed diets, linear mixing
Enamel.75C4 <- ((1 + d.C4/1e3 * 0.75 + d.C3/1e3 * 0.25)*(1 + epsilon/1e3) - 1) * 1e3 #-3.4 per mill

Enamel.50C4 <- ((1 + d.C4/1e3 * 0.5 + d.C3/1e3 * 0.5)*(1 + epsilon/1e3) - 1) * 1e3 #-7.1 per mill

Enamel.25C4 <- ((1 + d.C4/1e3 * 0.25 + d.C3/1e3 * 0.75)*(1 + epsilon/1e3) - 1) * 1e3 #-10.8 per mill
