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
           angle=30, length = 0.15, col = alpha(col[i], tr))
  }
}

####Fn 2: use diff(d18O) to assign plot colors####
O.plot.col <- function(data) {
  
  ind.data1 <- which(diff(data) > 0) #greater than 0
  
  ind.data0 <- which(diff(data) <= 0) #smaller than 0
  
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


###############plotting quasi-sympatric datasets###############
###############set 1: DR Congo ###############

#P01 Potamochoerus
plot(P01$X.13C.corr, P01$X.18O, main = "d13C - d18O corss plot: ",
     xlim = c(-16,2),ylim = c(-5,7), xlab="d13C", ylab="d18O",
     pch = 15, col = "black")
Plot.arrows(P01$X.13C.corr, P01$X.18O, col = O.plot.col(P01$X.18O))

#P04 Phacochoerus
points(P04$X.13C, P04$X.18O,
     pch = 17, col = "black")
Plot.arrows(P04$X.13C, P04$X.18O, col = O.plot.col(P04$X.18O))
legend(-16,6,c("Potamochoerus","Phacochoerus"),pch=c(15,17))

###############set 2: Kenya ###############

#P05 Potamochoerus
plot(P05$X.13C.corr, P05$X.18O, main = "d13C - d18O corss plot",
     xlim = c(-16,2),ylim = c(-5,7), xlab="d13C", ylab="d18O",
     pch = 15, col = "black")
Plot.arrows(P05$X.13C.corr, P05$X.18O, col = O.plot.col(P05$X.18O))

#P06 Phacochoerus
points(P06$X.13C, P06$X.18O,
       pch = 17, col = "black")
Plot.arrows(P06$X.13C, P06$X.18O, col = O.plot.col(P06$X.18O))

###############set 3: Malawi/Mozambique ###############

#P03 Potamochoerus
plot(P03$X.13C.corr, P03$X.18O, main = "d13C - d18O corss plot",
     xlim = c(-16,2),ylim = c(-5,7), xlab="d13C", ylab="d18O",
     pch = 15, col = "black")
Plot.arrows(P03$X.13C.corr, P03$X.18O, col = O.plot.col(P03$X.18O))

#P07 Phacochoerus
points(P07$X.13C, P07$X.18O,
       pch = 17, col = "black")
Plot.arrows(P07$X.13C, P07$X.18O, col = O.plot.col(P07$X.18O))