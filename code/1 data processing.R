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
