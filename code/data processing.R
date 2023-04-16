library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(dplyr)
library(ggplot2)
library(CircStats)
library(stats)

####helper functions, must run before analysis####

####Fn 1: plotting arrows####
Plot.arrows <- function(x.corr, y.corr, col) {
  if(length(x.corr) != length(y.corr)){
    warning("Error! x and y lengths differ")
  }else{
    for(i in 2:length(x.corr)){
      arrows(x.corr[i-1], y.corr[i-1], x.corr[i], y.corr[i],angle=30, length = 0.15, col = col)
    }
  }
}

####Fn 2: processing data to account for time lags in isotope profiles####
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


setwd("C:/Users/ydmag/Google Drive/Dissertation data/Extant suid seasonal diet/SuidDiet")

extant <- read.csv("data/suid_d13C_d18O.csv") #data from this study

ID.list <- levels(factor(extant$ID))
print(ID.list)

# [1] "NKU245" "NKU257" "NKU265" "P01"    "P03"    "P04"    "P05"    "P06"    "P07"    "SAM01"  "SAM02" 
# [12] "SAM03"  "SBL01"  "SBL02"  "SBL03"

NKU245 <- extant %>% filter(ID == ID.list[1])

NKU245.lag <- t.lag(NKU245$X.13C, NKU245$X.18O, -2)#shifting x back 2 positions

NKU245.lag.cor.k <- cor.test(NKU245.lag$X.13C, NKU245.lag$X.18O, method = "kendall")
NKU245.lag.cor.s <- cor.test(NKU245.lag$X.13C, NKU245.lag$X.18O, method = "spearman")

NKU245.lag.cor.k #NS
NKU245.lag.cor.s #NS

NKU257 <- extant %>% filter(ID == ID.list[2])

NKU257.lag <- t.lag(NKU257$X.13C, NKU257$X.18O, -2)#shifting x back 2 positions

NKU257.lag.cor.k <- cor.test(NKU257.lag$X.13C, NKU257.lag$X.18O, method = "kendall")
NKU257.lag.cor.s <- cor.test(NKU257.lag$X.13C, NKU257.lag$X.18O, method = "spearman")

NKU257.lag.cor.k #NS
NKU257.lag.cor.s #NS

NKU265 <- extant %>% filter(ID == ID.list[3])

NKU265.lag <- t.lag(NKU265$X.13C, NKU265$X.18O, -2)#shifting x back 2 positions

NKU265.lag.cor.k <- cor.test(NKU265.lag$X.13C, NKU265.lag$X.18O, method = "kendall")
NKU265.lag.cor.s <- cor.test(NKU265.lag$X.13C, NKU265.lag$X.18O, method = "spearman")

NKU265.lag.cor.k #NS, but close
NKU265.lag.cor.s #NS, but close

P01 <- extant %>% filter(ID == ID.list[4])

P01.lag <- t.lag(P01$X.13C, P01$X.18O, -2)#shifting x back 2 positions

P01.lag.cor.k <- cor.test(P01.lag$X.13C, P01.lag$X.18O, method = "kendall")
P01.lag.cor.s <- cor.test(P01.lag$X.13C, P01.lag$X.18O, method = "spearman")

P01.lag.cor.k #NS
P01.lag.cor.s #NS

P03 <- extant %>% filter(ID == ID.list[5])

P03.lag <- t.lag(P03$X.13C, P03$X.18O, -2)#shifting x back 2 positions

P03.lag.cor.k <- cor.test(P03.lag$X.13C, P03.lag$X.18O, method = "kendall")
P03.lag.cor.s <- cor.test(P03.lag$X.13C, P03.lag$X.18O, method = "spearman")

P03.lag.cor.k # Significant ***
P03.lag.cor.s # Significant ***

P04 <- extant %>% filter(ID == ID.list[6])

P04.lag <- t.lag(P04$X.13C, P04$X.18O, -2)#shifting x back 2 positions

P04.lag.cor.k <- cor.test(P04.lag$X.13C, P04.lag$X.18O, method = "kendall")
P04.lag.cor.s <- cor.test(P04.lag$X.13C, P04.lag$X.18O, method = "spearman")

P04.lag.cor.k # Significant **
P04.lag.cor.s # Significant **
 



#consider canine-molar gap:


#in a combined plot, use regular plot function
plot(P03$X.13C, P03$X.18O, main = "d13C - d18O corss plot",
     xlim = c(-12,1),ylim = c(-3,6),
     pch = 16, col = "blue")
Plot.arrows(P03$X.13C, P03$X.18O, col = "blue")
points(P04$X.13C, P04$X.18O,
     pch = 16, col = "red")
Plot.arrows(P04$X.13C, P04$X.18O, col = "red")


P05 <- extant %>% filter(ID == ID.list[7])

#visualize vectors
ggplot(P05, aes(x = X.13C, y = X.18O)) +
  geom_segment(aes(xend = c(tail(X.13C, n = -1), NA), 
                   yend = c(tail(X.18O, n = -1), NA)),
               arrow = arrow(length = unit(0.4, "cm")),
               color = 4) +
  geom_point(size = 2, color = 4) +
  geom_text(aes(label = serial, x = X.13C + 0.3, y = X.18O -0.3))

#visualize vectors
ggplot(P04, aes(x = X.13C, y = X.18O)) +
  geom_segment(aes(xend = c(tail(X.13C, n = -1), NA), 
                   yend = c(tail(X.18O, n = -1), NA)),
               arrow = arrow(length = unit(0.4, "cm")),
               color = 4) +
  geom_point(size = 2, color = 4) +
  geom_text(aes(label = serial, x = X.13C + 0.3, y = X.18O -0.3))


atan2(diff(P03$X.18O,1), diff(P03$X.13C,1))

#calculate slopes of vectors
P03.vtr.ang <- atan2(diff(P03$X.18O,1), diff(P03$X.13C,1)) #convert to 0 to 360 degrees

rose.diag(P03.vtr.ang,bins=12)

hist.P03.vtr <- hist(P03.vtr.ang,breaks=seq(-180,180,30))
plot.P03.vtr <- data.frame(hist.P03.vtr$mids,hist.P03.vtr$density)
colnames(plot.P03.vtr) <- c("mid","density")

ggplot(plot.P03.vtr, aes(x=mid, y=density)) +
  geom_bar(stat="identity") +
  coord_polar(theta = "x", start = -pi/45)+
  scale_x_continuous(breaks = seq(0, 360, 60))


frequency()

ggplot(wind, aes(x = DirCat, fill = SpeedCat)) +
  geom_histogram(binwidth = 15, boundary = -7.5) +
  coord_polar() +
  scale_x_continuous(limits = c(0,360))

#visualize vectors
ggplot(P03, aes(x = X.13C, y = X.18O)) +
  geom_segment(aes(xend = c(tail(X.13C, n = -1), NA), 
                   yend = c(tail(X.18O, n = -1), NA)),
               arrow = arrow(length = unit(0.4, "cm")),
               color = 4) +
  geom_point(size = 2, color = 4) +
  geom_text(aes(label = serial, x = X.13C + 0.3, y = X.18O -0.3))

plot(coord_polar())

#######import data from Yang et al. 2020, Phacochoerus from Mpala, KE###### 

#######import data from Reid et al. 2019, Phacochoerus from Naivasha, KE######
