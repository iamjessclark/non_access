#source("Setup3.R")
#source("objects.R") # this creates new populations

for(i in 1:length(pops)){ 
  for(j in 1:nrow(Values)){
    exclusion <- Values[j,1]
    compliance <- Values[j,2]
    
    Clone <- pops[[i]]$clone(deep=T)
    
    month_total <- 0
    month_index <- 0

    HALT <- FALSE
    
    Clone$assignexclusion(exclusion)
    
    while(HALT==FALSE){
      
      Clone$runTimestep()
      month_total <-  month_total + 1
      month_index <- month_index + 1
      
      if(month_total==12|month_total==24|month_total==36|month_total==48|month_total==60|month_total==72){
        Clone$runMDA(pop=1,towns=NA, coverage=0.65, drug="IA", compliance=compliance)
      }
      
      MWs <- which(Clone$WM>0)
      FWs <- which(Clone$WF>0)
      worms <- c(FWs, setdiff(MWs, FWs)) # updated revision 2 
      
      PrevAll[[i]][j,month_index] <- ((length(intersect(which(Clone$age>=5*12), worms))/length(which(Clone$age>=5*12)))*100)*0.931 # 0.931 is the 93.1% sensitivity #(length(which((Clone$age/12)>=5 & Clone$Mf>0))/length(which((Clone$age/12)>=5)))*100
      PrevAccessOnly[[i]][j,month_index] <- ((length(intersect(which(Clone$sysExProb==0 & Clone$age>=5*12), worms))/length(which(Clone$sysExProb==0 & Clone$age>=5*12)))*100)*0.931 #(length(which(Clone$sysExProb==0 & (Clone$age/12)>=5 & Clone$Mf>0))/(length(which(Clone$sysExProb==0 & (Clone$age/12)>=5))))*100
      #Larvae[[i]][j,month_index] <- Clone$larvae
      PrevKids[[i]][j,month_index] <- ((length(intersect(which(Clone$sysExProb==0 & Clone$age>=12*6 & Clone$age<12*8), worms))/(length(which(Clone$sysExProb==0 & Clone$age>=12*6 & Clone$age<12*8))))*100)*0.931
      Prev20[[i]][j,month_index] <- (length(which(Clone$sysExProb==0 & Clone$age>=20*12 & Clone$Mf>1))/(length(which(Clone$sysExProb==0 & Clone$age>=20*12))))*100 # using MF not worms 
      
      if(month_total==12){
        People[[i]][j,1] <- length(which(Clone$treat==T))
      }
      
      if(month_total==24){
        People[[i]][j,2] <- length(which(Clone$treat==T))
      }
      
      if(month_total==36){
        People[[i]][j,3] <- length(which(Clone$treat==T))
      }
      
      if(month_total==48){
        People[[i]][j,4] <- length(which(Clone$treat==T))
      }
      
      if(month_total==60){
        People[[i]][j,5] <- length(which(Clone$treat==T))
      }
      
      if(month_total==144){
        HALT <- TRUE
      }
      
    }
    
    PrevAccessOnly[[i]] <- as.data.frame(PrevAccessOnly[[i]])
    PrevAll[[i]] <- as.data.frame(PrevAll[[i]])
    #Larvae[[i]] <- as.data.frame(Larvae[[i]])
    PrevKids[[i]] <- as.data.frame(PrevKids[[i]])
    Prev20[[i]] <- as.data.frame(Prev20[[i]])
    People[[i]] <- as.data.frame(People[[i]])
  }
  
  PrevAccessOnly[[i]] <- cbind(PrevAccessOnly[[i]], Values)
  PrevAll[[i]] <- cbind(PrevAll[[i]], Values)
  #Larvae[[i]] <- cbind(Larvae[[i]], Values)
  PrevKids[[i]] <- cbind(PrevKids[[i]], Values)
  Prev20[[i]] <- cbind(Prev20[[i]], Values)
  People[[i]] <- cbind(People[[i]], Values)
  
}

# how many rounds to reach elimination - here we can then say if you are at 10 rounds but with a starting prevalence of x then it is more likely you have a real prevalence of xhat. 
MDAroundstotal1 <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=2)
MDArounds_All <- Rounds_Func(measuredpop="All", 
                             threshold = 1, 
                             pops = pops, 
                             Values = Values, 
                             trt = "IA", 
                             objectout = MDAroundstotal1)

MDAroundstotal2 <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=2)
MDArounds_Access <- Rounds_Func(measuredpop="AccessOnly", 
                             threshold = 1, 
                             pops = pops, 
                             Values = Values, 
                             trt = "IA", 
                             objectout = MDAroundstotal2)

MDAroundstotal3 <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=2)
MDArounds_kids <- Rounds_Func(measuredpop="TASkids", 
                                threshold = 1, 
                                pops = pops, 
                                Values = Values, 
                                trt = "IA", 
                                objectout = MDAroundstotal3)

MDAroundstotal4 <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=2)
MDArounds_mf20 <- Rounds_Func(measuredpop="mf20", 
                                threshold = 2, 
                                pops = pops, 
                                Values = Values, 
                                trt = "IA", 
                                objectout = MDAroundstotal4)


#save.image(file='myEnvironment.RData')
