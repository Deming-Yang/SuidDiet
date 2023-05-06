library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(dplyr)
library(tidyr)
library(ggplot2)
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


#####simple plot function######
#in a combined plot, use regular plot function
plot(P03$X.13C, P03$X.18O, main = "d13C - d18O corss plot",
     xlim = c(-12,1),ylim = c(-3,6),
     pch = 16, col = "blue")
Plot.arrows(P03$X.13C, P03$X.18O, col = "blue")
points(P04$X.13C, P04$X.18O,
     pch = 16, col = "red")
Plot.arrows(P04$X.13C, P04$X.18O, col = "red")

#####advanced plot function######
#arbitrarilly define seasons using d18O data
#version 1: use diff(P03$X.18O)



ind.P03.dry <- which(diff(P03$X.18O) > 0)

ind.P03.rain <- which(diff(P03$X.18O) <= 0)

#determine color bins
bin.wd <- 0.3

col.bins.P03.18O <- ceiling((max(diff(P03$X.18O))-min(diff(P03$X.18O)))/bin.wd) + 1

col.P03.18O <- viridis(col.bins.P03.18O)

P03.18O.plot.col <- c(1,col.P03.18O[ceiling((diff(P03$X.18O) + abs(min(diff(P03$X.18O))))/bin.wd) + 1])

plot(P03$X.13C, P03$X.18O, main = "d13C - d18O corss plot",
     xlim = c(-12,1),ylim = c(-3,6),
     pch = 16, col = "black")
Plot.arrows(P03$X.13C, P03$X.18O, col = P03.18O.plot.col)

#version 1: use diff(P03$X.18O)



plot(P03$X.13C, P03$X.18O, main = "d13C - d18O corss plot",
     xlim = c(-14,1),ylim = c(-6,6),
     pch = 16, col = "black")
Plot.arrows(P03$X.13C, P03$X.18O, col = P03.18O.plot.col)

#version 1: use diff(P03$X.18O)

ind.P04.dry <- which(diff(P04$X.18O) > 0)

ind.P04.rain <- which(diff(P04$X.18O) <= 0)

#initialize vector
P04.18O.plot.col <- rep(1,length(P04$X.18O))
P04.18O.plot.col[ind.P04.dry + 1] <- col.dry
P04.18O.plot.col[ind.P04.rain + 1] <- col.rain
P04.18O.plot.col

points(P04$X.13C, P04$X.18O,
     pch = 16, col = "black")
Plot.arrows(P04$X.13C, P04$X.18O, col = P04.18O.plot.col)

#####alternative: positino o differently#####
#version 1: use diff(P03$X.18O)

ind.P03.dry <- which(diff(P03$X.18O) >= 0)

ind.P03.rain <- which(diff(P03$X.18O) < 0)

col.w.d <- viridis(6)

col.dry <- col.w.d[5]

col.rain <- col.w.d[3]

#initialize vector
P03.18O.plot.col <- rep(1,length(P03$X.18O))
P03.18O.plot.col[ind.P03.dry + 1] <- col.dry
P03.18O.plot.col[ind.P03.rain + 1] <- col.rain
P03.18O.plot.col





#version 1: use diff(P03$X.18O)

ind.P04.dry <- which(diff(P04$X.18O) >= 0)

ind.P04.rain <- which(diff(P04$X.18O) < 0)

#initialize vector
P04.18O.plot.col <- rep(1,length(P04$X.18O))
P04.18O.plot.col[ind.P04.dry + 1] <- col.dry
P04.18O.plot.col[ind.P04.rain + 1] <- col.rain
P04.18O.plot.col

points(P04$X.13C, P04$X.18O,
       pch = 16, col = "black")
Plot.arrows(P04$X.13C, P04$X.18O, col = P04.18O.plot.col)

###############
ind.P05.dry <- which(diff(P05$X.18O) > 0)

ind.P05.rain <- which(diff(P05$X.18O) <= 0)

#initialize vector
P05.18O.plot.col <- rep(1,length(P05$X.18O))
P05.18O.plot.col[ind.P05.dry + 1] <- col.dry
P05.18O.plot.col[ind.P05.rain + 1] <- col.rain
P05.18O.plot.col

points(P05$X.13C, P05$X.18O,
       pch = 16, col = "black")
Plot.arrows(P05$X.13C, P05$X.18O, col = P05.18O.plot.col)


#####test synpatric individuals####
#version 1: use diff(P03$X.18O)

ind.P05.dry <- which(diff(P05$X.18O) > 0)

ind.P05.rain <- which(diff(P05$X.18O) <= 0)

col.w.d <- viridis(5)

col.dry <- col.w.d[4]

col.rain <- col.w.d[2]

#initialize vector
P05.18O.plot.col <- rep(1,length(P05$X.18O))
P05.18O.plot.col[ind.P05.dry + 1] <- col.dry
P05.18O.plot.col[ind.P05.rain + 1] <- col.rain
P05.18O.plot.col

plot(P05$X.13C, P05$X.18O, main = "d13C - d18O corss plot",
     xlim = c(-14,1),ylim = c(-6,6),
     pch = 16, col = "black")
Plot.arrows(P05$X.13C, P05$X.18O, col = P05.18O.plot.col)

#version 1: use diff(P03$X.18O)

ind.P06.dry <- which(diff(P06$X.18O) > 0)

ind.P06.rain <- which(diff(P06$X.18O) <= 0)

#initialize vector
P06.18O.plot.col <- rep(1,length(P06$X.18O))
P06.18O.plot.col[ind.P06.dry + 1] <- col.dry
P06.18O.plot.col[ind.P06.rain + 1] <- col.rain
P06.18O.plot.col

points(P06$X.13C, P06$X.18O,
       pch = 16, col = "black")
Plot.arrows(P06$X.13C, P06$X.18O, col = P06.18O.plot.col)

#########################################################################
#tests#
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
