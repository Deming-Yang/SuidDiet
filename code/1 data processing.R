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
P07 <- extant %>% filter(ID == ID.list[9])
SAM01 <- extant %>% filter(ID == ID.list[10])
SAM02 <- extant %>% filter(ID == ID.list[11])
SAM03 <- extant %>% filter(ID == ID.list[12])
SBL01 <- extant %>% filter(ID == ID.list[13])
SBL02 <- extant %>% filter(ID == ID.list[14])
SBL03 <- extant %>% filter(ID == ID.list[15])

#################order by dist to ERJ#################
P01 <- P01[order(P01$dist, decreasing = T),]
P03 <- P03[order(P03$dist, decreasing = T),]
P05 <- P05[order(P05$dist, decreasing = T),]
NKU245 <- NKU245[order(NKU245$dist, decreasing = T),]
NKU257 <- NKU257[order(NKU257$dist, decreasing = T),]
NKU265 <- NKU265[order(NKU265$dist, decreasing = T),]
P04 <- P04[order(P04$dist, decreasing = T),]
P06 <- P06[order(P06$dist, decreasing = T),]
P07 <- P07[order(P07$dist, decreasing = T),]
SAM01 <- SAM01[order(SAM01$dist, decreasing = T),]
SAM02 <- SAM02[order(SAM02$dist, decreasing = T),]
SAM03 <- SAM03[order(SAM03$dist, decreasing = T),]
SBL01 <- SBL01[order(SBL01$dist, decreasing = T),]
SBL02 <- SBL02[order(SBL02$dist, decreasing = T),]
SBL03 <- SBL03[order(SBL03$dist, decreasing = T),]

################correction for molar-canine d13C spacing###############
P01 <- P01 %>% mutate( X.13C1750.corr = X.13C1750 - 2)

P03 <- P03 %>% mutate( X.13C1750.corr = X.13C1750 - 2)

P05 <- P05 %>% mutate( X.13C1750.corr = X.13C1750 - 2)

################ calculate end member dietary values after correction for the suess effect################

#Cerling et al 2015

d.C4 <- -10

d.C3 <- -26.6 

#closed canopy
d.C3.CC <- -32.6

epsilon <- 13.3 #Passy 2005, diet to molar enamel fractionation in domestic pigs

#nominal pure C4 diet, linear mixing
Enamel.C4 <- ((1 + d.C4/1e3)*(1 + epsilon/1e3) - 1) * 1e3 #3.2 per mill

#nominal pure C3 diet, linear mixing
Enamel.C3 <- ((1 + d.C3/1e3)*(1 + epsilon/1e3) - 1) * 1e3 #-13.6 per mill

#nominal closed canopy C3 diet
Enamel.C3.CC <- ((1 + d.C3.CC/1e3)*(1 + epsilon/1e3) - 1) * 1e3 #-19.7 per mill

#mixed diets, linear mixing
Enamel.75C4 <- ((1 + d.C4/1e3 * 0.75 + d.C3/1e3 * 0.25)*(1 + epsilon/1e3) - 1) * 1e3 #-1 per mill

Enamel.50C4 <- ((1 + d.C4/1e3 * 0.5 + d.C3/1e3 * 0.5)*(1 + epsilon/1e3) - 1) * 1e3 #-5.2 per mill

Enamel.25C4 <- ((1 + d.C4/1e3 * 0.25 + d.C3/1e3 * 0.75)*(1 + epsilon/1e3) - 1) * 1e3 #-9.4 per mill
