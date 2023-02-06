require(tidyverse)

#### running model to find how many rounds of trt needed to reach threshold under different observation patterns ####

Rounds_Func <- function(measuredpop, threshold, pops, Values, trt, objectout){
  for(i in 1:length(pops)){
    for(j in 1:nrow(Values)){
      exclusion <- Values[j,1]
      compliance <- 0.2
      
      Clone <- pops[[i]]$clone(deep=T)
      
      HALT <- FALSE
      month_ref <- 11 # start counting months just before MDA (MDA on month 12 in this example)
      month_total <- 1 # overall month count
      maxTime <- 20 # maximum number of years to run for
      MDArounds <- 0 # number of MDA rounds completed
      treatment = trt # Ivermection and albendazole
      MDAcov <- 0.65 # 65% coverage
      
      # run with MDA until mf < 1%, storing monthly prevalence
      while(HALT == FALSE){
        
        Clone$runTimestep() # run a time step
        month_ref = month_ref + 1
        
        if(month_ref == 12){ # if month 12 then do MDA
          Clone$runMDA(pop = 1, towns = NA, coverage=MDAcov, drug=treatment, compliance = compliance )
          month_ref <- 0 # reset month reference to 0 (start of new year)
          MDArounds <- MDArounds +1 # record MDA round
        }
        
        MWs <- which(Clone$WM>0)
        FWs <- which(Clone$WF>0)
        worms <- intersect(MWs, FWs)
        
        PrevAll_index <- ((length(intersect(which(Clone$age>=5*12), worms))/length(which(Clone$age>=5*12)))*100)*0.931
        PrevAccessOnly_index <- ((length(intersect(which(Clone$sysExProb==0 & Clone$age>=5*12), worms))/length(which(Clone$sysExProb==0 & Clone$age>=5*12)))*100)*0.931
        PrevKids_index <- ((length(intersect(which(Clone$sysExProb==0 & Clone$age>=12*6 & Clone$age<12*8), worms))/(length(which(Clone$sysExProb==0 & Clone$age>=12*6 & Clone$age<12*8))))*100)*0.931
        Prev20_index <- (length(which(Clone$sysExProb==0 & (Clone$age/12)>=20 & Clone$Mf>0))/(length(which(Clone$sysExProb==0 & (Clone$age/12)>=20))))*100 # using MF not worms 
        
        month_ref <- month_ref + 1
        month_total <- month_total +1
        
        if(measuredpop == "mf20"){
          if(Prev20_index <= threshold & month_ref==11){ # If <1% just before next round then stop
            EPHP <- TRUE # Assumed EPHP achieved
            HALT <- TRUE
          }
        }
        
        if(measuredpop == "All"){
          if(PrevAll_index <= threshold & month_ref==11){ # If <1% just before next round then stop
            EPHP <- TRUE # Assumed EPHP achieved
            HALT <- TRUE
          }
        }
        
        if(measuredpop == "AccessOnly"){
          if(PrevAccessOnly_index <= threshold & month_ref==11){ # If <1% just before next round then stop
            EPHP <- TRUE # Assumed EPHP achieved
            HALT <- TRUE
          }
        }
        
        if(measuredpop == "TASkids"){
          if(PrevKids_index <= threshold & month_ref==11){ # If <1% just before next round then stop
            EPHP <- TRUE # Assumed EPHP achieved
            HALT <- TRUE
          }
        }
        
        if(month_total >= maxTime*12){ # If exceed maxTime then stop
          EPHP <- FALSE
          HALT <- TRUE
        }
      }
      objectout[[i]][j,1] <- MDArounds
      objectout[[i]][j,2] <- EPHP
    }
  }
  
  return(objectout)
}

# for simulating indian populations 

#### sample populations from the initial runs that match your distribution of prevalence ####
drawUniformIDs2 <- function(Prevalence,minPrev,maxPrev,N,Nbins){
  # N is number of draws, Nbins in number of bins (categories)
  
  unifMin <- minPrev
  unifMax <- maxPrev
  
  RefPrev <- seq(unifMin,unifMax,length.out = Nbins)
  NeachRef <- ceiling(N/(Nbins-1)) #we generate a few more than N
  PrevInter <- cut(Prevalence,RefPrev)
  
  levels <- as.vector(levels(PrevInter))
  
  ids <- as.vector(sapply(levels(PrevInter),function(x) base::sample(which(PrevInter==x),size=NeachRef,replace=T)))
  ids <- sample(ids,N,replace = F) # keep only N
  return(ids)
}


#### data manip ####
Manip2Plot <- function(outputlist){
  outputflat <- flatten_dfr(outputlist)
  outputflat$time <- as.factor(outputflat$time)
  Q <- outputflat %>%
    group_by(time, compliance, exclusion) %>%
    summarise(prev_quant=quantile(mf, c(0.025, 0.5, 0.975)), quant=c(0.025, 0.5, 0.975))%>%
    pivot_wider(names_from = quant, values_from=prev_quant)
  colnames(Q) <- c("time", "compliance", "exclusion", "2.5%","50%", "97.5%")
  Q$time <- as.numeric(as.character(Q$time))
  return(Q)
}

#### plot function ####

PlotFunc <- function(data){
  plot <- ggplot()+ 
    geom_line(data=data, aes(x=time,y=`50%`), colour="purple") +
    geom_ribbon(data=data, aes(x=time, ymin=`2.5%`, ymax=`97.5%`), alpha=0.2, fill="purple")+
    labs(x="time (years)",y="mf prevalence (%)") +
    theme_minimal()+
    facet_grid(exclusion~compliance)+
    geom_hline(yintercept=1)
  return(plot)
}

#### make empty list ####
MakeList <- function(SimPopsName, ExVecName, AdVecName ){
  list.rep <- rep(list(list()),length(SimPopsName)) # this is inner level
  listout <- rep(list(list.rep), length(ExVecName)) # this is middle level
  listout <- rep(list(listout), length(AdVecName)) # this is outer level
  return(listout)
}


runModSysEx <- function(NoYears, exclusion, coverage, treatment, compliance, pop){
  
  HALT <- FALSE
  month_ref <- 11 # start counting months just before MDA
  month_total <- 1 # overall month counter, starting point
  
  prevs <- vector()
  prevs <- length(which(pop$sysExProb==0 & pop$Mf>0))/length(which(pop$sysExProb==0)) # intial value of prevalence is the baseline pre-treatment post burn-in
  maxTime <- NoYears # number of years
  MDArounds <- 0 # starting number of MDA rounds at the beginning of the session
  
  while(HALT==FALSE){
    
    pop$runTimestep() # run a timestep
    
    if(month_ref==12){
      pop$runMDA(pop=1,towns=NA, exclusion=exclusion, coverage=coverage, drug=treatment, compliance=compliance)
      month_ref <- 0
      MDArounds <- MDArounds + 1
    }
    
    mfprev <- length(which(pop$sysExProb==0 & pop$Mf>0))/length(which(pop$sysExProb==0))
    month_ref <- month_ref + 1
    month_total <- month_total + 1
    prevs[month_total] <- mfprev
    
    if(mfprev <= 0.01 & month_ref==11){
      EPHP <- TRUE
      HALT <- TRUE
    }
    
    if(month_total>= maxTime*12){
      EPHP <- FALSE
      HALT <- TRUE
    }
  }
  pre_df <- tibble(time=(1:month_total)/12,mf=100*prevs)
  return(pre_df)
}


GatesPlotFunc <- function(data){
  plot <- ggplot()+ 
    geom_line(data=data, aes(x=time,y=`50%`), colour="purple") +
    geom_ribbon(data=data, aes(x=time, ymin=`2.5%`, ymax=`97.5%`), alpha=0.2, fill="purple")+
    labs(x="time (years)",y="mf prevalence (%)") +
    theme_minimal()+
    geom_hline(yintercept=1)
  return(plot)
}

Manip2PlotGates <- function(outputlist){
  outputflat <- do.call(rbind,outputlist)
  outputflat$time <- as.factor(outputflat$time)
  Q <- outputflat %>%
    group_by(time) %>%
    summarise(prev_quant=quantile(mf, c(0.025, 0.5, 0.975)), quant=c(0.025, 0.5, 0.975))%>%
    pivot_wider(names_from = quant, values_from=prev_quant)
  colnames(Q) <- c("time", "2.5%","50%", "97.5%")
  Q$time <- as.numeric(as.character(Q$time))
  return(Q)
}

#### make communities ####
MakingMultiPops <- function(totalpops, populationsize){
  newPop <- list()
  for(i in 1:totalpops){
    # Set up population characteristics
    # In this analysis considering only one population so popData will only have one row
    nPops <- 1 #number of towns/settlements (set to 1 if not using movement functionality)
    popSize <- populationsize #mean population size
    popData <- generatePops(nPops,popSize)
    popData$prev <- NA
    popData$nation <- 1 #number of countries = 1
    dMatrix <- matrix(0,1,1) # only one population, so just a zero matrix
    
    # Up population characteristics to generate a new population
    newPop[[i]] <- Population$new(popData,dMatrix)
    
    # Burn in period defaults to 100 years (1200 months) to allow prevalence to equilibriate
    newPop[[i]]$burnin()
    print(i)
  }
  return(newPop)
}
#### function to look at the prevs ####

PrintPrevs <- function(poplist){
  newpopprevs <- vector()
  for(i in 1:length(poplist)){
    newpopprevs[i] <- (length(which((poplist[[i]]$age/12)>=5 & poplist[[i]]$Mf>0))/(poplist[[i]]$nHosts-length(which((poplist[[i]]$age/12)<5))))*100
  }
  return(newpopprevs)
}

#### reset the populations ####
reset <- function(populationlist, populationprevs){
  popindex <- which(populationprevs>=1)
  ClonePops <- list()
  for(i in 1:length(popindex)){
    ClonePops[[i]] <- populationlist[[popindex[i]]]$clone(deep=T)
  }
  return(ClonePops)
}

#### how many kids####
AgeStructure <- function(population){
  newpopages <- vector()
  for(i in 1:length(population)){
    newpopages[i] <- length(which(population[[i]]$age>=6 & population[[i]]$age<=7))
  }
  return(newpopages) 
}

#### Larvae data manipulation #####

LarvaeData <- function(LarvePop1, LarvePop2){
  LarvePop1 <- do.call(rbind, LarvePop1)
  colnames(LarvePop1) <- c("larvae", "exclusion", "compliance", "time")
  LarvePop2 <- do.call(rbind, LarvePop2)
  colnames(LarvePop2) <- c("larvae", "exclusion","compliance", "time")
  LarvePop1$population <- "low prevalence"
  LarvePop2$population <- "high prevalence"
  Larvae <- bind_rows(LarvePop1, LarvePop2)
  Larvae$population <- as.factor(Larvae$population)
  Larvae$exclusion[which(Larvae$exclusion==0)] <- "0%"
  Larvae$exclusion[which(Larvae$exclusion==0.1)] <- "10%"
  Larvae$exclusion[which(Larvae$exclusion==0.2)] <- "20%"
  Larvae$exclusion[which(Larvae$exclusion==0.3)] <- "30%"
  Larvae$exclusion[which(Larvae$exclusion==0.4)] <- "40%"
  Larvae$exclusion[which(Larvae$exclusion==0.5)] <- "50%"
  Larvae$exclusion <- as.factor(Larvae$exclusion)
  return(Larvae)
}

#### make pops ####
MakeHigh <- function(totalpops, populationsize){
  HighR0 <- list()
  for(i in 1:totalpops){
    # Set up population characteristics
    # In this analysis considering only one population so popData will only have one row
    nPops <- 1 #number of towns/settlements (set to 1 if not using movement functionality)
    popSize <- populationsize #mean population size
    popData <- generatePops(nPops,popSize)
    popData$prev <- NA
    popData$nation <- 1 #number of countries = 1
    dMatrix <- matrix(0,1,1) # only one population, so just a zero matrix
    
    # Up population characteristics to generate a new population
    HighR0[[i]] <- Population$new(popData,dMatrix)
    HighR0[[i]]$AssignHighVectorVars()   
    # Burn in period defaults to 100 years (1200 months) to allow prevalence to equilibriate
    HighR0[[i]] $burnin()
  }
  return(HighR0)
}

MakeLow <- function(totalpops, populationsize){
  lowR0 <- list()
  for(i in 1:totalpops){
    # Set up population characteristics
    # In this analysis considering only one population so popData will only have one row
    nPops <- 1 #number of towns/settlements (set to 1 if not using movement functionality)
    popSize <- populationsize #mean population size
    popData <- generatePops(nPops,popSize)
    popData$prev <- NA
    popData$nation <- 1 #number of countries = 1
    dMatrix <- matrix(0,1,1) # only one population, so just a zero matrix
    
    # Up population characteristics to generate a new population
    MadeLowPop <- Population$new(popData,dMatrix)
    MadeLowPop$AssignLowVectorVars()   
    # Burn in period defaults to 100 years (1200 months) to allow prevalence to equilibriate
    MadeLowPop$burnin()
    lowR0[[i]] <- MadeLowPop
  }
  return(lowR0)
} 


#### dealing with output lists to make figures for revisions  ####
prob_clear_list_func <- function(listname, clear_limit, Values, prevs){
  for(i in 1:length(listname)){
    #listname[[i]] <- cbind(listname[[i]], Values)
    listname[[i]]$population <- i
    listname[[i]]$starting <- prevs[i]
  }
  
  dataframe <- do.call(rbind, listname)
  
  shape <-  dataframe %>%
    rename(exclusion = Var1, adherence = Var2) %>%
    dplyr::select(exclusion, adherence, population, starting)
  
  dataframe <- dataframe %>%
    dplyr::select(-Var1, -Var2, -population, -starting) 
  
  renamecols <- seq(1:144)/12
  colnames(dataframe) <- renamecols
  
  dataframe <- bind_cols(dataframe, shape) %>%
    mutate(starting_prev = case_when(starting >= 70  ~ "70-75%",
                                     starting >= 65  ~ "65-69%",
                                     starting >= 60  ~ "60-64%",
                                     starting >= 55  ~ "55-59%",
                                     starting >= 50  ~ "50-54%",
                                     starting >= 45  ~ "45-49%",
                                     starting >= 40  ~ "40-44%",
                                     starting >= 35  ~ "35-39%",
                                     starting >= 30  ~ "30-34%",
                                     starting >= 25  ~ "25-29%",
                                     starting >= 20  ~ "20-24%",
                                     starting >= 15  ~ "15-19%",
                                     starting >= 10  ~ "10-14%",
                                     starting >= 5   ~ "05-9%", 
                                     starting > 0     ~ "0-4%")) %>%
    pivot_longer(-c(exclusion, adherence, population, starting, starting_prev), names_to = "Time", values_to = "Prevalence") 
  
  # take out NA's - these are populations that had no prevalence to start with
  dataframe <- dataframe[-which(is.na(dataframe$starting_prev)),]
  
  year5 <- dataframe %>% 
    filter(Time == 5) %>% 
    group_by(exclusion, adherence, starting_prev) %>% 
    summarise(n = n())
  
  number_clear = dataframe %>% 
    filter(Time == 5) %>% 
    group_by(exclusion, adherence, starting_prev) %>% 
    summarise(clear = sum(Prevalence <= clear_limit)) %>%
    merge(year5) %>%
    mutate(prob_clear = clear/n)
  
  return(number_clear)
}

#### dynamics when reaching threshold #### 

prev_dynnamics_func <- function(listname, Values, threshold, prevs){
  renamecols <- seq(1:144)/12

  for(i in 1:length(listname)){  
    listname[[i]] <- listname[[i]][,-c(145:146)]
    colnames(listname[[i]]) <- renamecols
    listname[[i]]$population <- i
    listname[[i]] <- cbind(listname[[i]], Values)
    
    listname[[i]] <- listname[[i]] %>% 
      rename(exclusion = Var1, adherence = Var2) 
    
    listname[[i]] <- listname[[i]] %>%
      pivot_longer(cols = -c(exclusion, adherence, population), names_to = "time", values_to = "prevalence")
  }   
  
    dataframe <- do.call(rbind, listname)
    IDs <- which(prevs<5)
    dataframe <- dataframe[!dataframe$population %in% IDs,]
    dataframe <- dataframe %>%
      filter(adherence == 0.2)
    
    EPHPops <- dataframe %>% 
      mutate(as_factor(population)) %>% 
      group_by(population) %>% 
      filter(time==5 & prevalence <= threshold)
    
    IDs <- unique(EPHPops$population)
    
    dataframe <- dataframe[dataframe$population %in% IDs,]
    
    return(dataframe)
  
}

#### Treatmetn rounds ####

Treatment_Rounds <- function(listname, observation, prevs, Values){
  for(i in 1:length(listname)){
    listname[[i]] <- as.data.frame(cbind(listname[[i]] , Values))
    listname[[i]]$population <- i
    listname[[i]]$starting <- prevs[i]
    listname[[i]]$observed <- observation
  }
  
  listnamedf <- do.call(rbind, listname)
  colnames(MDArounds_Alldf) <- c("rounds", "EPHP", "exclusion", "adherence", "population", "starting", "observed")
  listnamedf <- listnamedf %>%
    mutate(starting_prev = case_when(starting >= 70  ~ "70-75%",
                                     starting >= 65  ~ "65-69%",
                                     starting >= 60  ~ "60-64%",
                                     starting >= 55  ~ "55-59%",
                                     starting >= 50  ~ "50-54%",
                                     starting >= 45  ~ "45-49%",
                                     starting >= 40  ~ "40-44%",
                                     starting >= 35  ~ "35-39%",
                                     starting >= 30  ~ "30-34%",
                                     starting >= 25  ~ "25-29%",
                                     starting >= 20  ~ "20-24%",
                                     starting >= 15  ~ "15-19%",
                                     starting >= 10  ~ "10-14%",
                                     starting >= 5   ~ "05-9%", 
                                     starting > 0     ~ "0-4%"))
  return(listnamedf)
}

#### larvae ####

larvaefunc <- function(larvaelist){
  
  for(i in 1:length(larvaelist)){
    
    larvaelist[[i]]$population <- i
    larvaelist[[i]]$starting <- prevs[i]
  }
  
  Larvaedf <- do.call(rbind, larvaelist)
  shape <-  Larvaedf %>%
    rename(exclusion = Var1, adherence = Var2) %>%
    dplyr::select(exclusion, adherence, population, starting)
  
  Larvaedf <- Larvaedf %>%
    dplyr::select(-Var1, -Var2, -population, -starting) 
  
  renamecols <- seq(1:144)/12
  colnames(Larvaedf) <- renamecols
  
  dataframe <- bind_cols(Larvaedf, shape) %>%
    mutate(starting_prev = case_when(starting >= 70  ~ "70-75%",
                                     starting >= 65  ~ "65-69%",
                                     starting >= 60  ~ "60-64%",
                                     starting >= 55  ~ "55-59%",
                                     starting >= 50  ~ "50-54%",
                                     starting >= 45  ~ "45-49%",
                                     starting >= 40  ~ "40-44%",
                                     starting >= 35  ~ "35-39%",
                                     starting >= 30  ~ "30-34%",
                                     starting >= 25  ~ "25-29%",
                                     starting >= 20  ~ "20-24%",
                                     starting >= 15  ~ "15-19%",
                                     starting >= 10  ~ "10-14%",
                                     starting >= 5   ~ "05-9%", 
                                     starting > 0     ~ "0-4%")) %>%
    filter(adherence == 0.2) 
dataframe <-   dataframe %>%
    pivot_longer(-c(exclusion, adherence, population, starting, starting_prev), names_to = "time", values_to = "l3density") 
  
  return(dataframe)
  
}
