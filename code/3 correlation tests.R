library(scales)
library(viridisLite)
library(dplyr)
library(tidyr)
library(stats)

#######import data from Yang et al. 2020, Phacochoerus from Mpala, KE###### 
#####load in previously published data####
Yang.2020 <- read.csv("data/Yang et al 2020.csv") #data from this study

Yang.2020.delim <- Yang.2020 %>% separate(Original.ID,into=c("ID","serial"), 
                                          sep="-", convert = TRUE, extra = "merge")

Yang.2020.ID.list <- levels(factor(Yang.2020.delim$ID))
Yang.2020.ID.list

MPL1C <- Yang.2020.delim %>% filter(ID == Yang.2020.ID.list[1])

MPL1M <- Yang.2020.delim %>% filter(ID == Yang.2020.ID.list[2])

MPL2C <- Yang.2020.delim %>% filter(ID == Yang.2020.ID.list[3])

MPL2M <- Yang.2020.delim %>% filter(ID == Yang.2020.ID.list[4])

#######import data from Reid et al. 2019, Phacochoerus from Naivasha, KE######
#####load in previously published data####
Reid.2019 <- read.csv("data/Reid et al 2019.csv") #data from this study

Reid.2019.ID.list <- levels(factor(Reid.2019$Sample.ID))

B119 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[1])

B33 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[2])

B384 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[3])

B56.1 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[4])

B58.2 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[5])

###################potamochoerus canines#################
par(mfrow=c(1,3))

ccf(P01$X.13C.corr, P01$X.18O, main = "Cross-correlation: P01")

ccf(P03$X.13C.corr, P03$X.18O, main = "Cross-correlation: P03")

ccf(P05$X.13C.corr, P05$X.18O, main = "Cross-correlation: P05")

###################Phacochoerus molars#################
par(mfrow=c(4,3))

ccf(NKU245$X.13C, NKU245$X.18O, main = "Cross-correlation: NKU245")

ccf(NKU257$X.13C, NKU257$X.18O, main = "Cross-correlation: NKU257")

ccf(NKU265$X.13C, NKU265$X.18O, main = "Cross-correlation: NKU265")

ccf(P04$X.13C, P04$X.18O, main = "Cross-correlation: P04")

ccf(P06$X.13C, P06$X.18O, main = "Cross-correlation: P06")

ccf(P07$X.13C, P07$X.18O, main = "Cross-correlation: P07")

ccf(SAM01$X.13C, SAM01$X.18O, main = "Cross-correlation: SAM01")

ccf(SAM02$X.13C, SAM02$X.18O, main = "Cross-correlation: SAM02")

ccf(SAM03$X.13C, SAM03$X.18O, main = "Cross-correlation: SAM03")

ccf(SBL01$X.13C, SBL01$X.18O, main = "Cross-correlation: SBL01")

ccf(SBL02$X.13C, SBL02$X.18O, main = "Cross-correlation: SBL02")

ccf(SBL03$X.13C, SBL03$X.18O, main = "Cross-correlation: SBL03")

par(mfrow=c(3,3))

ccf(MPL1C$X.13C, MPL1C$X.18O, main = "Cross-correlation: MPL1C")

ccf(MPL1M$X.13C, MPL1M$X.18O, main = "Cross-correlation: MPL1M")

ccf(MPL2C$X.13C, MPL2C$X.18O, main = "Cross-correlation: MPL2C")

ccf(MPL2M$X.13C, MPL2M$X.18O, main = "Cross-correlation: MPL2M")

ccf(B119$X.13C, B119$X.18O, main = "Cross-correlation: B119")

ccf(B33$X.13C, B33$X.18O, main = "Cross-correlation: B33")

ccf(B384$X.13C, B384$X.18O, main = "Cross-correlation: B384")

ccf(B56.1$X.13C, B56.1$X.18O, main = "Cross-correlation: B56.1")

ccf(B58.2$X.13C, B58.2$X.18O, main = "Cross-correlation: B58.2")

