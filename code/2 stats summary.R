library(stats)
library(dplyr)
library(tidyr)

#use customized function to streamline the code
stats.summ <- function(array){
  n <- length(array)
  mean <- mean(array)
  sd <- sd(array)
  median <- median(array)
  min <- min(array)
  max <- max(array)
  range <- max(array) - min(array)
  res <- data.frame(n, mean, sd, median, min, max, range)
  colnames(res)<- c("N","Mean", "sd", "Median","Minimum", "Maximum", "Range")
  print(res, digits=3)
}

########stats summary of d13C#######
###################potamochoerus canines#################
stats.summ(P01$X.13C.corr)

stats.summ(P03$X.13C.corr)

stats.summ(P05$X.13C.corr)

###################Phacochoerus molars#################
stats.summ(NKU245$X.13C)

stats.summ(NKU257$X.13C)

stats.summ(NKU265$X.13C)

stats.summ(P04$X.13C)

stats.summ(P06$X.13C)

stats.summ(P07$X.13C)

stats.summ(SAM01$X.13C)

stats.summ(SAM02$X.13C)

stats.summ(SAM03$X.13C)

stats.summ(SBL01$X.13C)

stats.summ(SBL02$X.13C)

stats.summ(SBL03$X.13C)

########stats summary of d18O#######

stats.summ(P01$X.18O)

stats.summ(P03$X.18O)

stats.summ(P05$X.18O)

stats.summ(NKU245$X.18O)

stats.summ(NKU257$X.18O)

stats.summ(NKU265$X.18O)

stats.summ(P04$X.18O)

stats.summ(P06$X.18O)

stats.summ(P07$X.18O)

stats.summ(SAM01$X.18O)

stats.summ(SAM02$X.18O)

stats.summ(SAM03$X.18O)

stats.summ(SBL01$X.18O)

stats.summ(SBL02$X.18O)

stats.summ(SBL03$X.18O)
