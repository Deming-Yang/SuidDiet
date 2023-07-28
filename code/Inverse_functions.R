
#### Emeas function ####


Emeasfun <- function(numtrials, numsam, length, dMeas, r1, r2, r3, la){
  require("matrixStats")
  
  # Edist <- array(data = NA, dim =c(numsam,numtrials))
  # 
  # allTrials <- array(data = NA, dim =c(numsam,numtrials))
  
  Edist <- rep(0,numtrials)
  allTrials <- array(data=NA, dim=c(numtrials, numsam,3))
  
  for (j in 1:numtrials) {
    #vector on ones
    Isoslope <- rep(1,numsam)
    
    for (n in 1:numsam) {
      if(n == 1) Isoslope[n] <- (((abs(dMeas[n+1]-dMeas[n]))/(0.5*(length[n+1]+length[n]))))
      else
        if(n == numsam) Isoslope[n] <- (((abs(dMeas[n]-dMeas[n-1]))/(0.5*(length[n]+length[n-1]))))
        else
          Isoslope[n]=((((( (abs(dMeas[n+1] - dMeas[n])) / ( 0.5*(length[n+1]+length[n]) ))   
                          +   ( (abs(dMeas[n] - dMeas[n-1])) /(0.5* (length[n]+length[n-1]))))/2)))
    }
    
    RElength = r2*rnorm(numsam,mean = 0, sd = 1);		#generates vector of random length errors from normal distribution
    LengthError = Isoslope*RElength 	#calculates isotope error caused by length errors (eq **.**)  
    LEdMeas = LengthError + dMeas   	#adds LengthError to the measured data vector    
    
    #The following loops check to see if the 'error' isotope value is above or below 
    #the highest or lowest adjacent values in dMeas. 
    #as this would be a physically impossible result of sample length error
    
    for (n in 1:numsam) {
      
      if(n == 1){
        if(LEdMeas[n] < min(dMeas[n:n+1])){
          LEdMeas[n] <- (min(dMeas[n:n+1]))
        } else if(LEdMeas[n] > max(dMeas[n:n+1])){
          LEdMeas[n] <- (max(dMeas[n:n+1]))
        } else LEdMeas[n] <- LEdMeas[n]
      } else if(n == numsam){
        if(LEdMeas[n] > max(dMeas[n:n-1]))
          LEdMeas[n] <- min(dMeas[n:n-1])
        else if(LEdMeas[n] > max(dMeas[n:n-1])){
          LEdMeas[n] <- max(dMeas[n:n-1])
        } else LEdMeas[n] <- LEdMeas[n]
      } else
        if(LEdMeas[n] < min(dMeas[n-1:n+1])){
          LEdMeas[n] <- min(dMeas[n-1:n+1])
        } else if(LEdMeas[n] > max(dMeas[n-1:n+1])){
          LEdMeas[n] <- max(dMeas[n-1:n+1])
        } else LEdMeas[n] <- LEdMeas[n]
        
    }
    
    # The following section calculates the depth-dependent isotope error.  It fits a 
    # cubic spline to the the measured data dMeas, 
    # and creates a vector of interpolated delta values for each unit length on the x-axis. 
    # A new vector is created that is shifted la 
    # units, and the two vectors are subtracted to give DELTA-delta values reflecting the 
    # difference between the isotope ratio at the outside 
    # enamel surface and the enamel-dentine junction.  These values are then multiplied by 
    # the sampling depth uncertainty to give a per mil uncertainty.  
    
    totallength = rep(1,numsam)
    totallength[1] = length[1]
    
    for (n in 2:numsam) {
      totallength[n] <- totallength[n-1] + length[n]
    }
    
    addbefore = dMeas[1]*rep(1,la)
    
    addafter = dMeas[numsam]*rep(1,la)
    
    xx <- 0:totallength[numsam]
    
    dInterp = spline(totallength,dMeas,xout = xx)
    dInterp <- dInterp$y
    dInterp
    
    dInterpShift = c(addbefore,dInterp)  
    
    dInterp = c(dInterp, addafter)  
    
    Dd = dInterpShift-dInterp  
    
    Deltadelta = rep(1, numsam) 
    
    for (n in 1:numsam) {
      Deltadelta[n] = Dd[totallength[n]]
    }
    
    
    REdepth = (r3*rnorm(numsam, mean = 0, sd = 1))/la          #generates vector of random depth errors
    
    DepthError = Deltadelta*REdepth
    
    
    #end of depth error section   
    
    AnalysisError = r1*rnorm(numsam, mean = 0, sd = 1) 
    length(AnalysisError)
    
    SqSumError = (LengthError^2 + DepthError^2 + AnalysisError^2)^0.5  
    length(SqSumError)
    E = SqSumError %*% SqSumError #these are not even matrices
    
    Edist[j] <- E  
    
    # allTrials[j] = SqSumError  #declear allTrials as an array
    
    dMeasError = dMeas + SqSumError  #vector
    
    allTrials[j,,] <- cbind(totallength, dMeas, dMeasError)
    
  }   
  
  
  # totallength = rep(1,numsam) 
  # totallength[1] = length[1]  
  # 
  # for (n in 2:numsam) {
  #   totallength[n] = totallength[n-1]+length[n]
  # }
  
  # output <- list(cbind(totallength, dMeas, dMeasError), allTrials)
  
  output <- list(allTrials, Edist)
  
  return(output)
}

#### mSolv function ####

mSolv_fun <- function(nsolxns, numsam, finit, la, lm, openindx, avelength, maxlength, minlength, mindepth, length, dMeas, r1, r2, r3,
                      maxratio, minratio, stdev, df, depth){
  

  
  require("matrixStats")
  
  pb <- progress_bar$new(total = nsolxns, show_after = 0)
  
  for (n in 1:nsolxns) {
    pb$tick()
    Sys.sleep(1 / nsolxns)
    
    # addition of random error and definitions
    
    ###############################
    
    # add random error to length vector, by matrix multiplying the length precision
    # with a vector of the length numsam of random numbers between 0 and 2
    rlength = r2*rnorm(numsam)
    length = round(rlength, digits=0) + length #round to integers and add error
    
    #check that lengths are between min and max threshold specified previously
    for (a in 1:numsam) {
      if (length[a]>maxlength){
        length[a]=maxlength
      }
      if(length[a]<minlength){
        length[a]=minlength
      }
      else {
        length[a]=length[a]
      }
    }
    
    rdepth=r3*rnorm(numsam)
    rdepth = round(rdepth, digits = 0)
    depth1 = round(rdepth,digits=0) + depth
    
    
    for (z in 1:numsam) {
      if (depth1[z]>la){
        depth1[z]=la
      }
      if(depth1[z]<mindepth){
        depth1[z]=mindepth
      }
      else {
        depth1[z]=depth1[z]
      }
    }
    
    depth=depth1
    
    #Determines the number of m's distal and proximal to those that directly correspond with d's
    #numbefore reflects the m's that are sampled into at the beginning of the profile because of sampling depth
    #numafter reflects the m's that contribute to the isotope values of the
    #final samples in open-ended cases.
    
    #round to nearest higher integer
    numbefore=ceiling(la/avelength)
    numafter=ceiling((lm-openindx)/avelength)+1
    
    lengthbefore=avelength*rep(1,numbefore)
    lengthafter=avelength*rep(1,numafter)
    
    #more definitions
    fmat = 1 - finit
    numcol <<- numbefore+numsam+numafter
    length=c(lengthbefore,length,lengthafter,lm)
    depth=c(depth,lengthbefore)
    
    ###############################
    
    # constructing the averaging matrix
    
    # 1) binary matrix
    
    ##############################
    
    B=matrix(0,avelength,numcol)#make matrix (fill,nrows,ncols) full of zeroes
    B[,1]=1# fill first column with ones
    for (m in 2:numcol) {
      Fill=matrix(0,length[m],numcol)
      Fill[,m]=1
      B=rbind(B,Fill)
    }
    
    #expand B matrix with additional non contributing rows
    newnrow = sum(length) - nrow(B)
    addB = matrix(0, newnrow, numcol)
    B = rbind(B, addB)
    
    ncrow1 = (sum(length)-openindx)+1 # set non contributing rows
    ncrow2 = sum(length)
    
    B[ncrow1:ncrow2,] = 0 # fill noncontributing rows with 0
    
    #remove last column
    B <- B[,-numcol]
    
    ##################### maturation average of binary matrix ###################
    
    o = 1
    
    AB = (finit*(colmeans(B[o:(o+lm-1),])) + fmat*(colmeans(B[o:(o+lm-1),])))
    AB = AB/(sum(B[o,])) # make seed row for AB matrix, will be deleted later
    
    for (o in 1:(sum(length)-lm)) {
      p = (finit * (B[o,]) + fmat * (colmeans(B[o:(o+lm-1),])))
      if (sum(p) == 0){
        p = p
      } else p = (p/sum(p))#makes all rows sum to one
      AB = rbind(AB,p)
    }
    
    AB <- AB[-1,]#remove seed row
    
    ######################## cumulative length vector ##########################
    
    clength = (length[1])
    
    for (q in 2:numcol) {
      cl = length[q] + clength[q-1]
      clength = c(clength, cl)
    }
    
    ######################## final calculation of A #############################
    
    A = AB[1,]
    
    for (k in (numbefore + 1):(numsam + numbefore)) {
      E = AB[1,]
      for (j in 0:(depth[k - numbefore] - 1)) {
        e = colmeans(AB[(clength[k-1] - j+1):(clength[k] - j),])
        E = rbind(E, e)
      }
      E = E[-1,]
      meanE = colmeans(E)
      A = rbind(A, meanE)
    }
    
    A = A[-1,]
    
    ############################# inversion ######################################
    
    I = diag(numsam)
    
    dMeasr = dMeas + r1 * rnorm(numsam, mean = 0, sd = 1)
    
    NoB = numbefore
    
    NoA = numafter - 1
    
    mm = rep(1, (numsam+numbefore+numafter-1))
    
    for (x in 1:(numsam+numbefore+numafter-1)) {
      
      mm[x] = ((maxratio + minratio)/2 +  ((maxratio - minratio)/2) * cos( 2*pi*(x)/(34*1.2))) + stdev * rnorm(1, mean = 0, sd = 0.8)
    }
    
    # end of reference vector determination
    
    mEst <<- mm + (t(A) %*% (solve(A %*% t(A) + df * I))) %*% (dMeasr - A %*% mm)
    
    MEST[,n] <<- mEst
    dPred = A %*% mEst
    dpe = dPred - dMeas
    DPE[n] <<- t(dpe) %*% dpe
    S[,n] = dpe
    
  }
  
  dim = dim(MEST)
  dim = dim[1]
  trimMEST = MEST
  trimMEST = trimMEST[-(dim-NoA+1:dim),]
  trimMEST = trimMEST[-(1:NoB),]
  
  avlM = colmeans(t(MEST))    #average of all input signals, including preceding and open-ended m's.
  stlM = colSds(t(MEST))    #standard deviation of all input signals, including preceding and open-ended m's.
  
  avsM = colmeans(t(trimMEST)) #average of all input signals, with m's not spatially corresponding to d's excluded.
  stsM = colSds(t(trimMEST))    #standard deviation of all input signals, with m's not spatially corresponding to d's excluded.
  
  
  ############# plotting ###############
  
  #sets up x-axis
  
  totallength = rep(1, numsam)
  totallength[1] = length[1]
  
  xx = avelength*rep(1, numbefore)
  zz = avelength*rep(1, numafter-1)
  totallength = c(xx,totallength,zz)
  
  vec1 = dMeasr[1]*rep(1, numbefore)
  vec2 = dMeasr[numsam]*rep(1, numafter-1)
  dMeasd = c(vec1,dMeasr,vec2)
  
  for (n in 2:(numsam+numbefore+numafter-1)) {
    totallength[n] = (totallength[n-1]+length[n])
  }
  
  callength <- totallength - (numbefore + 1) * avelength #correct for data before
  
  solvout <- cbind(callength, dMeasd, MEST)
  
  ntrials=1:nsolxns
  tcols=paste("trial",ntrials,sep="")
  
  colnames(solvout) = c("callength", "dMeasd", tcols)
  
  solvout <<- as.data.frame(solvout)
  
  # return(solvout)
  
}