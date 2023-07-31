library(scales)
library(viridisLite)

####Fn 1: plotting arrows####
Plot.arrows <- function(x.corr, y.corr, col, tr = NA) {
  require(scales)
  if(length(x.corr) != length(y.corr)){
    warning("Error! x and y lengths differ")
    break
  }else if(is.na(tr)){
    tr <- 0.5 #default transparancy
  }
  for(i in 2:length(x.corr)){
    arrows(x.corr[i-1], y.corr[i-1], x.corr[i], y.corr[i],lwd=3,
           angle=30, length = 0.12, col = alpha(col[i], tr))
  }
}

####Fn 2: use diff(d18O) to assign plot colors####
O.plot.col <- function(data) {
  
  ind.data1 <- which(diff(data) >= 0) #greater than 0
  
  ind.data0 <- which(diff(data) < 0) #smaller than 0
  
  require(viridisLite)
  
  col <- viridis(5)
  
  col.1 <- col[4]
  
  col.0 <- col[2]
  
  #initialize vector
  data.plot.col <- rep(1,length(data))
  data.plot.col[ind.data1 + 1] <- col.1
  data.plot.col[ind.data0 + 1] <- col.0
  
  return(data.plot.col)
  
}

#Fn 3
# Adds time series to plot w 2 prob density envelopes
tsdens = function(d, base = "black"){
  require(scales)
  #Check dimensions of d
  if(ncol(d) != 4){stop("d cols should be should be time, 2.5%, 50%, 97.5% CI")}
  
  base.rgb = col2rgb(base)
  cols = c(rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255, alpha = 0.25), 
           rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255, alpha = 1))
  
  polygon(c(d[, 1], rev(d[, 1])), c(d[, 2], rev(d[, 4])), 
          col = cols[1], border = NA) #95% CI
  lines(d[, 1], d[, 3], col = cols[3], lwd = 2)
}


###############plotting quasi-sympatric datasets###############
###############set 1: DR Congo ###############
par(mfrow=c(3,1))
#P01 Potamochoerus
plot(P01$X.13C1750.corr, P01$X.18O, main = "d13C - d18O biplot: ",
     xlim = c(-15,2),ylim = c(-5,7), xlab="d13C", ylab="d18O",
     pch = 15, col = "black",cex = 1.5)
Plot.arrows(P01$X.13C1750.corr, P01$X.18O, col = O.plot.col(P01$X.18O))

#P04 Phacochoerus
points(P04$X.13C1750, P04$X.18O,
     pch = 17, col = "black",cex = 1.5)
Plot.arrows(P04$X.13C1750, P04$X.18O, col = O.plot.col(P04$X.18O))
legend(-15,6,c("Potamochoerus","Phacochoerus"),pch=c(15,17))

###############set 2: Kenya ###############

#P05 Potamochoerus
plot(P05$X.13C1750.corr, P05$X.18O, main = "d13C - d18O biplot",
     xlim = c(-15,2),ylim = c(-5,7), xlab="d13C", ylab="d18O",
     pch = 15, col = "black",cex = 1.5)
Plot.arrows(P05$X.13C1750.corr, P05$X.18O, col = O.plot.col(P05$X.18O))

#P06 Phacochoerus
points(P06$X.13C1750, P06$X.18O,
       pch = 17, col = "black",cex = 1.5)
Plot.arrows(P06$X.13C1750, P06$X.18O, col = O.plot.col(P06$X.18O))

###############set 3: Malawi/Mozambique ###############

#P03 Potamochoerus
plot(P03$X.13C1750.corr, P03$X.18O, main = "d13C - d18O biplot",
     xlim = c(-15,2),ylim = c(-5,7), xlab="d13C", ylab="d18O",
     pch = 15, col = "black",cex = 1.5)
Plot.arrows(P03$X.13C1750.corr, P03$X.18O, col = O.plot.col(P03$X.18O))

#P07 Phacochoerus
points(P07$X.13C1750, P07$X.18O,
       pch = 17, col = "black",cex = 1.5)
Plot.arrows(P07$X.13C1750, P07$X.18O, col = O.plot.col(P07$X.18O))

#############supplementary figures#######################
######Supplementary Fig S1 cross-correlation analysis#######
par(mfrow=c(1,3))

ccf(P01$X.18O, P01$X.13C1750.corr, main = "Cross-correlation: P01")

ccf(P03$X.18O, P03$X.13C1750.corr, main = "Cross-correlation: P03")

ccf(P05$X.18O, P05$X.13C1750.corr, main = "Cross-correlation: P05")

###################Phacochoerus molars#################
######Supplementary Fig S2 cross-correlation analysis#######
par(mfrow=c(4,3))

ccf(NKU245$X.18O, NKU245$X.13C1750, main = "Cross-correlation: NKU245")

ccf(NKU257$X.18O, NKU257$X.13C1750, main = "Cross-correlation: NKU257")

ccf(NKU265$X.18O, NKU265$X.13C1750, main = "Cross-correlation: NKU265")

ccf(P04$X.18O, P04$X.13C1750, main = "Cross-correlation: P04")

ccf(P06$X.18O, P06$X.13C1750, main = "Cross-correlation: P06")

ccf(P07$X.18O, P07$X.13C1750, main = "Cross-correlation: P07")

ccf(SAM01$X.18O, SAM01$X.13C1750, main = "Cross-correlation: SAM01")

ccf(SAM02$X.18O, SAM02$X.13C1750, main = "Cross-correlation: SAM02")

ccf(SAM03$X.18O, SAM03$X.13C1750, main = "Cross-correlation: SAM03")

ccf(SBL01$X.18O, SBL01$X.13C1750, main = "Cross-correlation: SBL01")

ccf(SBL02$X.18O, SBL02$X.13C1750, main = "Cross-correlation: SBL02")

ccf(SBL03$X.18O, SBL03$X.13C1750, main = "Cross-correlation: SBL03")

######Supplementary Fig S3 cross-correlation analysis#######
par(mfrow=c(3,3))

ccf(MPL1C$X.18O, MPL1C$X.13C, main = "Cross-correlation: MPL1C")

ccf(MPL1M$X.18O, MPL1M$X.13C, main = "Cross-correlation: MPL1M")

ccf(MPL2C$X.18O, MPL2C$X.13C, main = "Cross-correlation: MPL2C")

ccf(MPL2M$X.18O, MPL2M$X.13C, main = "Cross-correlation: MPL2M")

ccf(B119$X.18O, B119$X.13C, main = "Cross-correlation: B119")

ccf(B33$X.18O, B33$X.13C, main = "Cross-correlation: B33")

ccf(B384$X.18O, B384$X.13C, main = "Cross-correlation: B384")

ccf(B56.1$X.18O, B56.1$X.13C, main = "Cross-correlation: B56.1")

ccf(B58.2$X.18O, B58.2$X.13C, main = "Cross-correlation: B58.2")

###### Supplementary Fig S4 results of inverse modeling #######
par(mfrow=c(2,2))
######### P03 d13C inverse output + 95% CI #####
plot(P03.all.out.13C$ci.length, P03.all.out.13C$mean, type = "l", lwd = 2, ylim = c(-15,1),
     main = "P03 d13C inverse output", xlab = "length (mm)", ylab = "d13C1750")
tsdens(P03.13C.CI) # add 95% CI as gray shading
points(max(P03$dist)-P03$dist, P03$X.13C1750.corr, col = "orange3", pch = 16, cex = 1.5) #measurements
lines(max(P03$dist)-P03$dist, P03$X.13C1750.corr, col = "orange3", lwd = 2, lty = 2)

######### P04 d13C inverse output + 95% CI #####
######### plot out 95% CI #####
plot(P04.all.out.13C$ci.length, P04.all.out.13C$mean, type = "l", lwd = 2, ylim = c(-10,5),
     main = "P04 d13C inverse output", xlab = "length (mm)", ylab = "d13C1750")
tsdens(P04.13C.CI) # add 95% CI as gray shading
points(max(P04$dist)-P04$dist, P04$X.13C, col = "orange3", pch = 16, cex = 1.5) #measurements
lines(max(P04$dist)-P04$dist, P04$X.13C, col = "orange3", lwd = 2, lty = 2)

######### P03 d18O inverse output + 95% CI #####
plot(P03.all.out.18O$ci.length, P03.all.out.18O$mean, type = "l", lwd = 2, ylim = c(-5,3),
     main = "P03 d18O inverse output", xlab = "length (mm)", ylab = "d18O")
tsdens(P03.18O.CI) # add 95% CI as gray shading
points(max(P03$dist)-P03$dist, P03$X.18O, col = "cyan4", pch = 16, cex = 1.5) #measurements
lines(max(P03$dist)-P03$dist, P03$X.18O, col = "cyan4", lwd = 2, lty = 2)

######### P04 d18O inverse output + 95% CI #####
plot(P04.all.out.18O$ci.length, P04.all.out.18O$mean, type = "l", lwd = 2, ylim = c(-1,8),
     main = "P04 d18O inverse output", xlab = "length (mm)", ylab = "d18O")
tsdens(P04.18O.CI) # add 95% CI as gray shading
points(max(P04$dist)-P04$dist, P04$X.18O, col = "cyan4", pch = 16, cex = 1.5) #measurements
lines(max(P04$dist)-P04$dist, P04$X.18O, col = "cyan4", lwd = 2, lty = 2)

######Supplementary Fig S4 cross-correlation analysis#######
#use mean output values of the inverse iterations

#remove the data points before and after the data series
P03.inv.13C.tr <- P03.all.out.13C$mean[2:21]

P03.inv.18O.tr <- P03.all.out.18O$mean[2:21]

#remove the data points before and after the data series
P04.inv.13C.tr <- P04.all.out.13C$mean[5:19]

P04.inv.18O.tr <- P04.all.out.18O$mean[5:19]

par(mfrow=c(2,2))

ccf(P03$X.18O, P03$X.13C1750.corr, main = "Cross-correlation: P03")

ccf(P04$X.18O, P04$X.13C1750, main = "Cross-correlation: P04")

ccf(P03.inv.18O.tr, P03.inv.13C.tr, main = "Cross-correlation: P03 inverse")

ccf(P04.inv.18O.tr, P04.inv.13C.tr, main = "Cross-correlation: P04 inverse")


######Supplementary Fig S5 reinterpretation of seasons based on dd18O/dt#######
#P03 Potamochoerus measured
par(mfrow=c(2,2))
plot(P03$X.13C1750.corr, P03$X.18O, main = "d13C - d18O biplot",
     xlim = c(-15,-1),ylim = c(-5,4), xlab="d13C", ylab="d18O",
     pch = 15, col = "black")
Plot.arrows(P03$X.13C1750.corr, P03$X.18O, col = O.plot.col(P03$X.18O))

#P04 Phacochoerus measured
plot(P04$X.13C1750, P04$X.18O, main = "d13C - d18O biplot",
     xlim = c(-9,4),ylim = c(0,7), xlab="d13C", ylab="d18O",
     pch = 15, col = "black")
Plot.arrows(P04$X.13C1750, P04$X.18O, col = O.plot.col(P04$X.18O))

plot(P03.inv.13C.tr, P03.inv.18O.tr, main = "d13C - d18O inverse output",
     xlim = c(-15,-1),ylim = c(-5,4), xlab="d13C", ylab="d18O",
     pch = 15, col = "black")
Plot.arrows(P03.inv.13C.tr, P03.inv.18O.tr, col = O.plot.col(P03.inv.18O.tr))

#P04 Phacochoerus inverse
plot(P04.inv.13C.tr, P04.inv.18O.tr, main = "d13C - d18O inverse output",
     xlim = c(-9,4),ylim = c(0,7), xlab="d13C", ylab="d18O",
     pch = 15, col = "black")
Plot.arrows(P04.inv.13C.tr, P04.inv.18O.tr, col = O.plot.col(P04.inv.18O.tr))