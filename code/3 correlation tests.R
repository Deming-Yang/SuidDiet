library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stats)

####helper functions, must run before analysis####
####Fn: processing data to account for time lags in isotope profiles####
t.lag <- function(x.series, y.series, x.lag) {
  #by default, the time lag is applided to x.series
  if(length(x.series) != length(y.series)){
    warning("Error! x and y lengths differ")
    break
  }else{
    n = length(x.series)
    if (x.lag > 0) {#shifting x to match y + x.lag
      new.x <- x.series[(1 + x.lag):n]
      new.y <- y.series[1 :(n - x.lag)]
    }else if(x.lag < 0) {#shifting y to match y - x.lag
      new.x <- x.series[1 :(n - abs(x.lag))]
      new.y <- y.series[(1 + abs(x.lag)):n]
      
    }else{#no shift
      new.x <- x.series
      new.y <- y.series
    }
    
    df <- data.frame(new.x, new.y)
    colnames(df) <- c("X.13C","X.18O")
    return(df)
  }
  
}

###################potamochoerus canines#################
P01.cor.k <- cor.test(P01$X.13C.corr, P01$X.18O, method = "kendall")
P01.cor.s <- cor.test(P01$X.13C.corr, P01$X.18O, method = "spearman")
P01.cor.k #NS
P01.cor.s #NS

P01.lag <- t.lag(P01$X.13C.corr, P01$X.18O, -2)#shifting x back 2 positions

P01.lag.cor.k <- cor.test(P01.lag$X.13C.corr, P01.lag$X.18O, method = "kendall")
P01.lag.cor.s <- cor.test(P01.lag$X.13C.corr, P01.lag$X.18O, method = "spearman")

P01.lag.cor.k #NS
P01.lag.cor.s #NS



P03.cor.k <- cor.test(P03$X.13C.corr, P03$X.18O, method = "kendall")
P03.cor.s <- cor.test(P03$X.13C.corr, P03$X.18O, method = "spearman")
P03.cor.k #NS
P03.cor.s #NS

P03.lag <- t.lag(P03$X.13C.corr, P03$X.18O, -2)#shifting x back 2 positions

P03.lag.cor.k <- cor.test(P03.lag$X.13C.corr, P03.lag$X.18O, method = "kendall")
P03.lag.cor.s <- cor.test(P03.lag$X.13C.corr, P03.lag$X.18O, method = "spearman")

P03.lag.cor.k # Significant ***
P03.lag.cor.s # Significant ***



P05.cor.k <- cor.test(P05$X.13C.corr, P05$X.18O, method = "kendall")
P05.cor.s <- cor.test(P05$X.13C.corr, P05$X.18O, method = "spearman")
P05.cor.k #NS
P05.cor.s #NS

P05.lag <- t.lag(P05$X.13C.corr, P05$X.18O, -2)#shifting x back 2 positions

P05.lag.cor.k <- cor.test(P05.lag$X.13C.corr, P05.lag$X.18O, method = "kendall")
P05.lag.cor.s <- cor.test(P05.lag$X.13C.corr, P05.lag$X.18O, method = "spearman")

P05.lag.cor.k # NS
P05.lag.cor.s # NS

###################Phacochoerus molars#################


NKU245.cor.k <- cor.test(NKU245$X.13C, NKU245$X.18O, method = "kendall")
NKU245.cor.s <- cor.test(NKU245$X.13C, NKU245$X.18O, method = "spearman")
NKU245.cor.k #significant *
NKU245.cor.s #significant **

NKU245.lag <- t.lag(NKU245$X.13C, NKU245$X.18O, -2)#shifting x back 2 positions

NKU245.lag.cor.k <- cor.test(NKU245.lag$X.13C, NKU245.lag$X.18O, method = "kendall")
NKU245.lag.cor.s <- cor.test(NKU245.lag$X.13C, NKU245.lag$X.18O, method = "spearman")

NKU245.lag.cor.k #NS
NKU245.lag.cor.s #NS



NKU257.cor.k <- cor.test(NKU257$X.13C, NKU257$X.18O, method = "kendall")
NKU257.cor.s <- cor.test(NKU257$X.13C, NKU257$X.18O, method = "spearman")
NKU257.cor.k #significant *
NKU257.cor.s #significant *

NKU257.lag <- t.lag(NKU257$X.13C, NKU257$X.18O, -2)#shifting x back 2 positions

NKU257.lag.cor.k <- cor.test(NKU257.lag$X.13C, NKU257.lag$X.18O, method = "kendall")
NKU257.lag.cor.s <- cor.test(NKU257.lag$X.13C, NKU257.lag$X.18O, method = "spearman")

NKU257.lag.cor.k #NS
NKU257.lag.cor.s #NS



NKU265.cor.k <- cor.test(NKU265$X.13C, NKU265$X.18O, method = "kendall")
NKU265.cor.s <- cor.test(NKU265$X.13C, NKU265$X.18O, method = "spearman")
NKU265.cor.k #significant *
NKU265.cor.s #significant *

NKU265.lag <- t.lag(NKU265$X.13C, NKU265$X.18O, -2)#shifting x back 2 positions

NKU265.lag.cor.k <- cor.test(NKU265.lag$X.13C, NKU265.lag$X.18O, method = "kendall")
NKU265.lag.cor.s <- cor.test(NKU265.lag$X.13C, NKU265.lag$X.18O, method = "spearman")

NKU265.lag.cor.k #NS, but close
NKU265.lag.cor.s #NS, but close





P04.cor.k <- cor.test(P04$X.13C, P04$X.18O, method = "kendall")
P04.cor.s <- cor.test(P04$X.13C, P04$X.18O, method = "spearman")
P04.cor.k #NS
P04.cor.s #NS

P04.lag <- t.lag(P04$X.13C, P04$X.18O, -2)#shifting x back 2 positions

P04.lag.cor.k <- cor.test(P04.lag$X.13C, P04.lag$X.18O, method = "kendall")
P04.lag.cor.s <- cor.test(P04.lag$X.13C, P04.lag$X.18O, method = "spearman")

P04.lag.cor.k # Significant **
P04.lag.cor.s # Significant **



P06.cor.k <- cor.test(P06$X.13C, P06$X.18O, method = "kendall")
P06.cor.s <- cor.test(P06$X.13C, P06$X.18O, method = "spearman")
P06.cor.k #NS
P06.cor.s #NS

P06.lag <- t.lag(P06$X.13C, P06$X.18O, -2)#shifting x back 2 positions

P06.lag.cor.k <- cor.test(P06.lag$X.13C, P06.lag$X.18O, method = "kendall")
P06.lag.cor.s <- cor.test(P06.lag$X.13C, P06.lag$X.18O, method = "spearman")

P06.lag.cor.k # NS
P06.lag.cor.s # NS

P07 <- extant %>% filter(ID == ID.list[9])

P07.cor.k <- cor.test(P07$X.13C, P07$X.18O, method = "kendall")
P07.cor.s <- cor.test(P07$X.13C, P07$X.18O, method = "spearman")
P07.cor.k #NS
P07.cor.s #NS

P07.lag <- t.lag(P07$X.13C, P07$X.18O, -2)#shifting x back 2 positions

P07.lag.cor.k <- cor.test(P07.lag$X.13C, P07.lag$X.18O, method = "kendall")
P07.lag.cor.s <- cor.test(P07.lag$X.13C, P07.lag$X.18O, method = "spearman")

P07.lag.cor.k # NS
P07.lag.cor.s # NS



SAM01.cor.k <- cor.test(SAM01$X.13C, SAM01$X.18O, method = "kendall")
SAM01.cor.s <- cor.test(SAM01$X.13C, SAM01$X.18O, method = "spearman")
SAM01.cor.k # Significant *
SAM01.cor.s # Significant **

SAM01.lag <- t.lag(SAM01$X.13C, SAM01$X.18O, -2)#shifting x back 2 positions

SAM01.lag.cor.k <- cor.test(SAM01.lag$X.13C, SAM01.lag$X.18O, method = "kendall")
SAM01.lag.cor.s <- cor.test(SAM01.lag$X.13C, SAM01.lag$X.18O, method = "spearman")

SAM01.lag.cor.k # NS
SAM01.lag.cor.s # NS



SAM02.cor.k <- cor.test(SAM02$X.13C, SAM02$X.18O, method = "kendall")
SAM02.cor.s <- cor.test(SAM02$X.13C, SAM02$X.18O, method = "spearman")
SAM02.cor.k # NS
SAM02.cor.s # NS

SAM02.lag <- t.lag(SAM02$X.13C, SAM02$X.18O, -2)#shifting x back 2 positions

SAM02.lag.cor.k <- cor.test(SAM02.lag$X.13C, SAM02.lag$X.18O, method = "kendall")
SAM02.lag.cor.s <- cor.test(SAM02.lag$X.13C, SAM02.lag$X.18O, method = "spearman")

SAM02.lag.cor.k # NS
SAM02.lag.cor.s # NS



SAM03.cor.k <- cor.test(SAM03$X.13C, SAM03$X.18O, method = "kendall")
SAM03.cor.s <- cor.test(SAM03$X.13C, SAM03$X.18O, method = "spearman")
SAM03.cor.k # NS
SAM03.cor.s # NS

SAM03.lag <- t.lag(SAM03$X.13C, SAM03$X.18O, -2)#shifting x back 2 positions

SAM03.lag.cor.k <- cor.test(SAM03.lag$X.13C, SAM03.lag$X.18O, method = "kendall")
SAM03.lag.cor.s <- cor.test(SAM03.lag$X.13C, SAM03.lag$X.18O, method = "spearman")

SAM03.lag.cor.k # Significant *
SAM03.lag.cor.s # Significant *



SBL01.cor.k <- cor.test(SBL01$X.13C, SBL01$X.18O, method = "kendall")
SBL01.cor.s <- cor.test(SBL01$X.13C, SBL01$X.18O, method = "spearman")
SBL01.cor.k # NS
SBL01.cor.s # NS

SBL01.lag <- t.lag(SBL01$X.13C, SBL01$X.18O, -2)#shifting x back 2 positions

SBL01.lag.cor.k <- cor.test(SBL01.lag$X.13C, SBL01.lag$X.18O, method = "kendall")
SBL01.lag.cor.s <- cor.test(SBL01.lag$X.13C, SBL01.lag$X.18O, method = "spearman")

SBL01.lag.cor.k # Significant *
SBL01.lag.cor.s # Significant *



SBL02.cor.k <- cor.test(SBL02$X.13C, SBL02$X.18O, method = "kendall")
SBL02.cor.s <- cor.test(SBL02$X.13C, SBL02$X.18O, method = "spearman")
SBL02.cor.k # NS
SBL02.cor.s # NS

SBL02.lag <- t.lag(SBL02$X.13C, SBL02$X.18O, -2)#shifting x back 2 positions

SBL02.lag.cor.k <- cor.test(SBL02.lag$X.13C, SBL02.lag$X.18O, method = "kendall")
SBL02.lag.cor.s <- cor.test(SBL02.lag$X.13C, SBL02.lag$X.18O, method = "spearman")

SBL02.lag.cor.k # Significant *
SBL02.lag.cor.s # Significant *



SBL03.cor.k <- cor.test(SBL03$X.13C, SBL03$X.18O, method = "kendall")
SBL03.cor.s <- cor.test(SBL03$X.13C, SBL03$X.18O, method = "spearman")
SBL03.cor.k # NS
SBL03.cor.s # NS

SBL03.lag <- t.lag(SBL03$X.13C, SBL03$X.18O, -2)#shifting x back 2 positions

SBL03.lag.cor.k <- cor.test(SBL03.lag$X.13C, SBL03.lag$X.18O, method = "kendall")
SBL03.lag.cor.s <- cor.test(SBL03.lag$X.13C, SBL03.lag$X.18O, method = "spearman")

SBL03.lag.cor.k # Significant *
SBL03.lag.cor.s # Significant *

#######import data from Yang et al. 2020, Phacochoerus from Mpala, KE###### 
#####load in previously published data####
Yang.2020 <- read.csv("data/Yang et al 2020.csv") #data from this study

Yang.2020.delim <- Yang.2020 %>% separate(Original.ID,into=c("ID","serial"), 
                                          sep="-", convert = TRUE, extra = "merge")

Yang.2020.ID.list <- levels(factor(Yang.2020.delim$ID))
Yang.2020.ID.list

MPL1C <- Yang.2020.delim %>% filter(ID == Yang.2020.ID.list[1])

MPL1C.cor.k <- cor.test(MPL1C$X.13C, MPL1C$X.18O, method = "kendall")
MPL1C.cor.s <- cor.test(MPL1C$X.13C, MPL1C$X.18O, method = "spearman")
MPL1C.cor.k # Significant *
MPL1C.cor.s # Significant *

MPL1C.lag <- t.lag(MPL1C$X.13C, MPL1C$X.18O, -2)#shifting x back 2 positions

MPL1C.lag.cor.k <- cor.test(MPL1C.lag$X.13C, MPL1C.lag$X.18O, method = "kendall")
MPL1C.lag.cor.s <- cor.test(MPL1C.lag$X.13C, MPL1C.lag$X.18O, method = "spearman")

MPL1C.lag.cor.k # Significant *
MPL1C.lag.cor.s # Significant *

MPL1M <- Yang.2020.delim %>% filter(ID == Yang.2020.ID.list[2])

MPL1M.cor.k <- cor.test(MPL1M$X.13C, MPL1M$X.18O, method = "kendall")
MPL1M.cor.s <- cor.test(MPL1M$X.13C, MPL1M$X.18O, method = "spearman")
MPL1M.cor.k # NS
MPL1M.cor.s # NS

MPL1M.lag <- t.lag(MPL1M$X.13C, MPL1M$X.18O, -2)#shifting x back 2 positions

MPL1M.lag.cor.k <- cor.test(MPL1M.lag$X.13C, MPL1M.lag$X.18O, method = "kendall")
MPL1M.lag.cor.s <- cor.test(MPL1M.lag$X.13C, MPL1M.lag$X.18O, method = "spearman")

MPL1M.lag.cor.k # Significant *
MPL1M.lag.cor.s # Significant *

MPL2C <- Yang.2020.delim %>% filter(ID == Yang.2020.ID.list[3])

MPL2C.cor.k <- cor.test(MPL2C$X.13C, MPL2C$X.18O, method = "kendall")
MPL2C.cor.s <- cor.test(MPL2C$X.13C, MPL2C$X.18O, method = "spearman")
MPL2C.cor.k # Significant *
MPL2C.cor.s # Significant *

MPL2C.lag <- t.lag(MPL2C$X.13C, MPL2C$X.18O, -2)#shifting x back 2 positions

MPL2C.lag.cor.k <- cor.test(MPL2C.lag$X.13C, MPL2C.lag$X.18O, method = "kendall")
MPL2C.lag.cor.s <- cor.test(MPL2C.lag$X.13C, MPL2C.lag$X.18O, method = "spearman")

MPL2C.lag.cor.k # Significant *
MPL2C.lag.cor.s # Significant *

MPL2M <- Yang.2020.delim %>% filter(ID == Yang.2020.ID.list[4])

MPL2M.cor.k <- cor.test(MPL2M$X.13C, MPL2M$X.18O, method = "kendall")
MPL2M.cor.s <- cor.test(MPL2M$X.13C, MPL2M$X.18O, method = "spearman")
MPL2M.cor.k # NS
MPL2M.cor.s # NS

MPL2M.lag <- t.lag(MPL2M$X.13C, MPL2M$X.18O, -2)#shifting x back 2 positions

MPL2M.lag.cor.k <- cor.test(MPL2M.lag$X.13C, MPL2M.lag$X.18O, method = "kendall")
MPL2M.lag.cor.s <- cor.test(MPL2M.lag$X.13C, MPL2M.lag$X.18O, method = "spearman")

MPL2M.lag.cor.k # Significant *
MPL2M.lag.cor.s # Significant *

#######import data from Reid et al. 2019, Phacochoerus from Naivasha, KE######
#####load in previously published data####
Reid.2019 <- read.csv("data/Reid et al 2019.csv") #data from this study

Reid.2019.ID.list <- levels(factor(Reid.2019$Sample.ID))

B119 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[1])

B119.cor.k <- cor.test(B119$X.13C, B119$X.18O, method = "kendall")
B119.cor.s <- cor.test(B119$X.13C, B119$X.18O, method = "spearman")
B119.cor.k # NS
B119.cor.s # NS

B119.lag <- t.lag(B119$X.13C, B119$X.18O, -2)#shifting x back 2 positions

B119.lag.cor.k <- cor.test(B119.lag$X.13C, B119.lag$X.18O, method = "kendall")
B119.lag.cor.s <- cor.test(B119.lag$X.13C, B119.lag$X.18O, method = "spearman")

B119.lag.cor.k # NS
B119.lag.cor.s # NS

B33 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[2])

B33.cor.k <- cor.test(B33$X.13C, B33$X.18O, method = "kendall")
B33.cor.s <- cor.test(B33$X.13C, B33$X.18O, method = "spearman")
B33.cor.k # NS
B33.cor.s # NS

B33.lag <- t.lag(B33$X.13C, B33$X.18O, -2)#shifting x back 2 positions

B33.lag.cor.k <- cor.test(B33.lag$X.13C, B33.lag$X.18O, method = "kendall")
B33.lag.cor.s <- cor.test(B33.lag$X.13C, B33.lag$X.18O, method = "spearman")

B33.lag.cor.k # NS
B33.lag.cor.s # NS

B384 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[3])

B384.cor.k <- cor.test(B384$X.13C, B384$X.18O, method = "kendall")
B384.cor.s <- cor.test(B384$X.13C, B384$X.18O, method = "spearman")
B384.cor.k # NS
B384.cor.s # NS

B384.lag <- t.lag(B384$X.13C, B384$X.18O, -2)#shifting x back 2 positions

B384.lag.cor.k <- cor.test(B384.lag$X.13C, B384.lag$X.18O, method = "kendall")
B384.lag.cor.s <- cor.test(B384.lag$X.13C, B384.lag$X.18O, method = "spearman")

B384.lag.cor.k # NS
B384.lag.cor.s # NS

B56.1 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[4])

B56.1.cor.k <- cor.test(B56.1$X.13C, B56.1$X.18O, method = "kendall")
B56.1.cor.s <- cor.test(B56.1$X.13C, B56.1$X.18O, method = "spearman")
B56.1.cor.k # NS
B56.1.cor.s # NS

B56.1.lag <- t.lag(B56.1$X.13C, B56.1$X.18O, -2)#shifting x back 2 positions

B56.1.lag.cor.k <- cor.test(B56.1.lag$X.13C, B56.1.lag$X.18O, method = "kendall")
B56.1.lag.cor.s <- cor.test(B56.1.lag$X.13C, B56.1.lag$X.18O, method = "spearman")

B56.1.lag.cor.k # NS
B56.1.lag.cor.s # NS

B58.2 <- Reid.2019 %>% filter(Sample.ID == Reid.2019.ID.list[5])

B58.2.cor.k <- cor.test(B58.2$X.13C, B58.2$X.18O, method = "kendall")
B58.2.cor.s <- cor.test(B58.2$X.13C, B58.2$X.18O, method = "spearman")
B58.2.cor.k # NS
B58.2.cor.s # NS

B58.2.lag <- t.lag(B58.2$X.13C, B58.2$X.18O, -2)#shifting x back 2 positions

B58.2.lag.cor.k <- cor.test(B58.2.lag$X.13C, B58.2.lag$X.18O, method = "kendall")
B58.2.lag.cor.s <- cor.test(B58.2.lag$X.13C, B58.2.lag$X.18O, method = "spearman")

B58.2.lag.cor.k # NS
B58.2.lag.cor.s # NS