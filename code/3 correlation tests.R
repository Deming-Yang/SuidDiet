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

MPL1C <- MPL1C[order(MPL1C$Distance, decreasing = T),]
MPL1M <- MPL1M[order(MPL1M$Distance, decreasing = T),]
MPL2C <- MPL2C[order(MPL1C$Distance, decreasing = T),]
MPL2M <- MPL2M[order(MPL1M$Distance, decreasing = T),]

#######import data from Reid et al. 2019, Phacochoerus from Naivasha, KE######
#####load in previously published data####
Reid.2019 <- read.csv("data/Reid et al 2019.csv") #data from this study

Reid.2019.ID.list <- levels(factor(Reid.2019$Sample.ID))

B119 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[1])

B33 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[2])

B384 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[3])

B56.1 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[4])

B58.2 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[5])

B119 <- B119[order(B119$Distance, decreasing = T),]
B33 <- B33[order(B33$Distance, decreasing = T),]
B384 <- B384[order(B384$Distance, decreasing = T),]
B56.1 <- B56.1[order(B56.1$Distance, decreasing = T),]
B58.2 <- B58.2[order(B58.2$Distance, decreasing = T),]


###################potamochoerus canines#################

ccf.P01 <- ccf(P01$X.13C1750.corr, P01$X.18O, plot = F)

ccf.P03 <- ccf(P03$X.13C1750.corr, P03$X.18O, plot = F)

ccf.P05 <- ccf(P05$X.13C1750.corr, P05$X.18O, plot = F)

###################Phacochoerus molars#################

ccf.NKU245 <- ccf(NKU245$X.13C1750, NKU245$X.18O, plot = F)

ccf.NKU257 <- ccf(NKU257$X.13C1750, NKU257$X.18O, plot = F)

ccf.NKU265 <- ccf(NKU265$X.13C1750, NKU265$X.18O, plot = F)

ccf.P04 <- ccf(P04$X.13C1750, P04$X.18O, plot = F)

ccf.P06 <- ccf(P06$X.13C1750, P06$X.18O, plot = F)

ccf.P07 <- ccf(P07$X.13C1750, P07$X.18O, plot = F)

ccf.SAM01 <- ccf(SAM01$X.13C1750, SAM01$X.18O, plot = F)

ccf.SAM02 <- ccf(SAM02$X.13C1750, SAM02$X.18O, plot = F)

ccf.SAM03 <- ccf(SAM03$X.13C1750, SAM03$X.18O, plot = F)

ccf.SBL01 <- ccf(SBL01$X.13C1750, SBL01$X.18O, plot = F)

ccf.SBL02 <- ccf(SBL02$X.13C1750, SBL02$X.18O, plot = F)

ccf.SBL03 <- ccf(SBL03$X.13C1750, SBL03$X.18O, plot = F)

ccf.MPL1C <- ccf(MPL1C$X.13C, MPL1C$X.18O, plot = F)

ccf.MPL1M <- ccf(MPL1M$X.13C, MPL1M$X.18O, plot = F)

ccf.MPL2C <- ccf(MPL2C$X.13C, MPL2C$X.18O, plot = F)

ccf.MPL2M <- ccf(MPL2M$X.13C, MPL2M$X.18O, plot = F)

ccf.B119 <- ccf(B119$X.13C, B119$X.18O, plot = F)

ccf.B33 <- ccf(B33$X.13C, B33$X.18O, plot = F)

ccf.B384 <- ccf(B384$X.13C, B384$X.18O, plot = F)

ccf.B56.1 <- ccf(B56.1$X.13C, B56.1$X.18O, plot = F)

ccf.B58.2 <- ccf(B58.2$X.13C, B58.2$X.18O, plot = F)

