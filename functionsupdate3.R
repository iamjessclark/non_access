###############################
### Functions Spatial model ###
###############################
library(R6)

generatePops <- function(nPops=10,sizePop=100){
  x <- runif(nPops,0,sizePop)
  y <- runif(nPops,0,sizePop)
  name <- LETTERS[1:nPops] #Only nPops < 27
  pop <- round(rnorm(nPops,sizePop,sd=10))

  return(data.frame(name,pop,x,y))
}

findBorder <- function(data,nPops){
  #define border as midway between two central populations
  #(i.e. same number of populations on each side)
  midpoint_x <- (data$x[order(data$x)][c(ceiling((nPops+1)/2))]+data$x[order(data$x)][c(floor((nPops+1)/2))])/2
  return(midpoint_x)
}

splitBorder <- function(data,nPops,midpoint){
  data$nation <- ifelse(data$x<midpoint,1,2) #classify Nation 1 or 2
  data <- data[order(data$x),]
  row.names(data)<-c()
  data$name <- LETTERS[1:nPops]
  return(data)
}

Parameters <- R6Class(
  "Parameters",
  portable = F,
  public = list(
    ## Host ##
    tau = 0.00167, #death rate hosts
    #mImpRate = .0006, #mean importation rate
    #maxImpRate = .001, #max importation rate
    commRate = 0.5, #commuting rate per timestep (month)
    mink = .01, #min k (overdispersion gamma bite risk)
    maxk = .1, #max k (overdispersion gamma bite risk)

    ## Vector ##
    species = 1, #0 culex, 1 anopheles
    initialL3 = 5, #initial larval density
    lbda = 10, #bite rate mosq/monthx
    infecMosq = 0.37, #Proportion of mosquitos that get infected
    kappas1 = 4.395, #L3 uptake and development
    r1 = 0.055, #L3 uptake and development
    deathMosq = 5, #death rate of mosquitos
    minvth = 0, #min V/H ratio = 0
    maxvth = 10, #max V/H ratio = 70

    ## Worm ##
    nu = 0, #poly-monogamy parameter
    alpha = 1, #mf birth rate per FW
    psi1 = .414, #Proportion L3 leaving mosquito per bite
    psi2 = .32, #Proportion of L3 leaving mosquito that enter host
    s2 = .00275, #Proportion of L3 entering host that develop into adults
    mu = .0104, #death rate worms
    gamma = .1, #mf death rate

    ## Intervention ##
    sysCompN = .99, #systematic non-compliance bednets
    rhoBU = 0, #correlation non-compliance MDA & bite risk
    rhoCN = 0, #correlation non-compliance MDA & bednets
    sysComp = 0.35, #systematic non-compliance - the correlation between rounds 
    fecreds = c(6,9,9),#6, #months of reduction in fecundity due to MDA (DA) - developmental period in the host is 6-12 months
    sysEx = 0, # set this using the model running function, this is the proportion of the pop that are excluded
    #IDA: 100% mf kill and 100% remove adult worms (to represent sterlisation of adult worms)
    #IA: 99% mf kill, 35% adult kill, 9 month sterilisation
    mfKillMDAs =c(.95,.99,1), #DA, IA, IDA
    wKillMDAs = c(.55,.35,1), #DA, IA, IDA
    drugRef = c('DA','IA','IDA') #diethylcarbamazine&Abdnzl, Ivrmctn&Abdnzl, or all three? 
  )
)

Population <- R6Class(
  "Population",
  portable = F,
  public = list (
    p = Parameters$new(), #Load parameters
    nHost = NULL, #Total number host in pop
    nHosts = NULL, #Vector with number of hosts in each pop/town
    nPops = NULL, #Number of populations/towns
    nationIndex = NULL, #Index of which population/town is in which nation
    treatCov = 0, #Treatment coverage
    vectorSpec = 0, #0 for culex, 1 for anopheles
    larvae = 50, #average mean free stage larvae
    k = NULL, #host overdispersion parameter
    VtH = NULL, #Vector to host ratio
    town = array(),
    nation = array(),
    WF = array(), #Female number worms
    WM = array(), #Male number worms
    Mf = array(), #Number microfilarae
    age = array(), #Age host (in months)
    weights = matrix(), #weights according to population commuter model chosen
    borderWeights = matrix(), #border adjusted weights
    Tij = matrix(),  #Average number of individuals moving from i to j in one timestep
    treatprob = array(), # the probability someone will be treated given the coverage and degree of correlation between treatment rounds. 
    sysExProb = array(), # the probablity of a person being someone who is systematicaly excluded from surveillance
    treat = F, #Treatment received
    ttreat = array(),
    bednet = F, #Bednet used
    bnCov = 0, #Bednet coverage
    biteRisk = NULL, #bi - bite risk (heterogeneity)
    biteRate = NULL, #h(a) - age dependent biting rate
    #border = 0.5, #border strength 0=never cross border, 1=always cross border, 0.5=no effect
    dMatrix = matrix(), #distances between towns
    moveMatrix = matrix(), #who's moving, where to and where from
    
    
    # Initialize function - specify individual ages/bite risks and sample population parameters (e.g. k, VtH)#
    initialize = function(data,mydMatrix,nationalVtH=F){
      self$nHosts <- data$pop
      self$nationIndex <- data$nation
      self$dMatrix <- mydMatrix
      nHost <<- sum(nHosts)
      nPops <<- length(nHosts)
      larvae <<- rep(larvae,nPops)
      Tij <<- matrix(0,nPops,nPops)
      town <<- unlist(sapply(seq_along(nHosts),function(x){rep(LETTERS[x],nHosts[x])}), use.names=F)
      WF <<- rep(0,nHost) #Number of female adult worms
      WM <<- rep(0,nHost) #Number of male adult worms
      Mf <<- rep(0,nHost) #Number of microfilaria
      age <<- runif(nHost,0,100)*12 #in months
      treat <<- rep(F,nHost) #Boolean treated or not
      ttreat <<- rep(12*15,nHost) #Time since treatment (default to very high - 15 years, so that there is no influence)
      #treatprob <<- rbeta(nHost-(nHost*p$sysEx),shape1=treatCov*(1-p$sysComp)/p$sysComp,shape2=(1-treatCov)*(1-p$sysComp)/p$sysComp) #Prob each individual received MDA
      sysExProb <<- rbinom(nHost, 1, p$sysEx)
      treatprob[which(sysExProb==0 & age>=5)] <<- rbeta(length(which(sysExProb==0 & age>=5)),shape1=treatCov*(1-p$sysComp)/p$sysComp,shape2=(1-treatCov)*(1-p$sysComp)/p$sysComp) # for those who are NOT part of that group this is their treatment probability 
      treatprob[which(sysExProb==1 | age<5)] <<- 0 # for those who are part of that group, their treatment probability is 0 - they will never receive treatment 
      whichDrug <<- NA
      wKillMDA <<- NA
      mfKillMDA <<- NA
      drugID <<- NA
      #bednet <<- rep(F,nHost) #Boolean, using bednets or not
      bednet <<- rbinom(nHost,1,bnCov)
      if(nationalVtH) {
        VtH <<- runif(max(self$nationIndex),p$minvth,p$maxvth)
        } # sample one VtH for each nation from a uniform
      if(!nationalVtH) { VtH <<- rep(runif(1,p$minvth,p$maxvth),length(nationIndex)) }
      k <<- runif(1,min = p$mink, max = p$maxk) # sample one k for the population from a uniform
      biteRisk <<- rgamma(nHost,shape = k, rate = k) #individual bite risk shape is k, mean=1 so rate=shape
      biteRate <<- pmin(age/(9*12),1) #age in months - linear increase max out at 9 years - age scale factor for bite risk
    },
    ## Do worm dynamics
    wormDynamics = function(){
      names <- LETTERS[1:nPops] #list of population names
      for (i in seq_along(names)){
        inds <- which(town==names[i]) #indexes of who is in population i
        townHost <- length(inds)
        bnRed <- ifelse(bednet[inds],.03,1) #bednet reduction, if individual has bednets reduce by 97%, if false no reduction
        WF[inds] <<- pmax(WF[inds] - rpois(townHost,p$mu*WF[inds]) +
                      rpois(townHost,.5 * p$lbda * biteRisk[inds] * VtH[nationIndex[i]] * p$psi1 * p$psi2 * p$s2 * biteRate[inds] * larvae[i] * bnRed),0)
        WM[inds] <<- pmax(WM[inds] - rpois(townHost,p$mu*WM[inds]) +
                      rpois(townHost,.5 * p$lbda * biteRisk[inds] * VtH[nationIndex[i]] * p$psi1 * p$psi2 * p$s2 * biteRate[inds] * larvae[i] * bnRed),0)
      }
    },
    ## Do mf dynamics
    mfDynamics = function(){
      I = ifelse(WM>0,1,0) #Need males for reproduction
      Tr = ifelse(treat==T,1,0) #If treatment, no new Mf
      Mf <<- pmax(Mf - p$gamma * Mf + p$alpha * WF * I * (1-Tr),0)
    },
    ## Do larval dynamics
    larvaeDynamics = function(){
      names <- LETTERS[1:nPops]
      for(i in seq_along(names)){
        inds <- which(town==names[i]) #indexes of who is in population i
        if(p$species==0){ #Culex
          L3 <- p$kappas1 * (1 - exp(-p$r1 * Mf[inds] / p$kappas1))
        } else { #Anopheles
          L3 <- p$kappas1 * (1 - exp(-p$r1 * Mf[inds] / p$kappas1))^2
        }
        L3 <- L3 * biteRisk[inds]
        bnRed <- ifelse(bnCov==0,1,1-bnCov*(1-.03)) #bednet reduction
        meanlarvae <- sum(L3)/sum(biteRisk[inds]) * bnRed #Sum across town population rather than whole population
        larvae[i] <<- pmax(p$lbda * p$infecMosq * meanlarvae / (p$deathMosq + p$lbda * p$psi1),0)
      }
    },
    ## runs one time-step
    runTimestep = function(burnin=F){
      if (burnin==F){
        runCommute()
      }
      wormDynamics()
      mfDynamics()
                   #rngImportation(x=which(rnorm(nHost)<(1-exp(-.0006))))
      larvaeDynamics()
      aging()
      if(burnin==F){
        returnHome()
        ttreat <<- ttreat + 1 #update how long it's been since treatment
        }
      treat[which(ttreat >= rbinom(length(ttreat),2*p$fecreds[drugID],0.5))] <<- F #stochastic change
      #treat[which(ttreat == p$fecred)] <<- F #step change

    },
    
    assignexclusion = function(exclusion){
      p$sysEx <<- exclusion
      sysExProb <<- rbinom(nHost, 1, p$sysEx)
    },
    
    AssignHighVectorVars = function(){
      p$lbda <<- 611
      p$minvth <<- 60
      p$maxvth <<- 150
    },
    
    AssignLowVectorVars = function(){
      p$lbda <<- 3
      p$minvth <<- 20
      p$maxvth <<- 40
    },

    runMDA = function(pop=1,towns=NA,coverage,drug, compliance){
      if(coverage != treatCov | compliance!= p$sysComp ){ #check syst non-compliance probabilities been set for this coverage level already, if not then set them
        #| exclusion!= p$sysEx
        treatCov <<- coverage
        p$sysComp <<- compliance 
        #p$sysEx <<- exclusion
        #sysExProb <<- rbinom(nHost, 1, p$sysEx)
        treatprob[which(sysExProb==0 & age>=5)] <<- rbeta(length(which(sysExProb==0 & age>=5)),shape1=treatCov*(1-p$sysComp)/p$sysComp,shape2=(1-treatCov)*(1-p$sysComp)/p$sysComp) # for those who are NOT part of that group this is their treatment probability 
        treatprob[which(sysExProb==1 | age<5)] <<- 0 # for those who are part of that group, their treatment probability is 0 - they will never receive treatment 
        #treatprob <<- rbeta(nHost-(nHost*p$sysEx),shape1=treatCov*(1-p$sysComp)/p$sysComp,shape2=(1-treatCov)*(1-p$sysComp)/p$sysComp)
      }
      
      drugID <<- which(p$drugRef == drug)
      wKillMDA <<- p$wKillMDAs[drugID]
      mfKillMDA <<- p$mfKillMDAs[drugID]
      if(is.na(towns)){
        names = LETTERS[1:nPops]
        towns = names[which(nationIndex == pop)]
      }
      
      treatPops = which(town %in% towns) # this selects the whole population associated with the right town
      treated = treatPops[which(rbinom(length(treatPops),1,treatprob)==1)] # this then randomly assigns threatment to each individual in the population with a probablity equal to louise's term 
      # insert a check to make sure that the treated group is the right proportion of the population 
      WM[treated] <<- rbinom(length(treated),WM[treated],1-wKillMDA)
      WF[treated] <<- rbinom(length(treated),WF[treated],1-wKillMDA)
      Mf[treated] <<- Mf[treated]*rbinom(Mf[treated],1,1-mfKillMDA)
      treat[treated] <<- T #this will cause no mf production
      ttreat[treated] <<- rep(0,length(treated))
      whichDrug <- drugID
      
    },
    ## Random importation
    #rngImportation = function(x,xi=10 * 9.2 * .414 * .32 * .00275,mu=.0104){
    #   WF[x] <<- rep(round(.5 * xi * biteRate * 10/mu),length=length(x))
    #   WM[x] <<- rep(round(.5 * xi * biteRate * 10/mu),length=length(x))
    # },
    
    ##Aging & births
    aging = function(){
      age <<- age+1
      resetHosts(which(runif(nHost) < 1-exp(-p$tau)| age>1200)) #simulates deaths/births
      biteRate <<- pmin(age/(9*12),1) #Update biteRate - age in months - linear increase max out at 9 years
    },
    ## Do burnin period - lb length burnin (in months)
    burnin = function(lb=1200){
      i=0
      repeat{
        runTimestep(burnin=T)
        i = i+1
        if (i==lb){ break}
      }
    },
    ## reset hosts
    resetHosts = function(ID){
      sysExProb[ID] <<- rbinom(1, 1, p$sysEx)
      age[ID] <<- 0 #birth age set to 0
      WF[ID] <<- 0 #female worms set to 0
      WM[ID] <<- 0 #male worms set to 0
      Mf[ID] <<- 0 #Mf set to 0
      treat[ID] <<- F #no treatment
      ttreat[ID] <<- 15*12
      biteRisk[ID] <<- rgamma(length(ID),shape = k, rate = k) #Get new biterisk - shape is k, mean=1 so rate=shape
      biteRate[ID] <<- pmin(age[ID]/(9*12),1)
      bednet[ID] <<- F #no bednet
    },
    #Gravity model of trade, mass*mass/distance (all variables to the power of alphas close to 1)
    weightsTrade = function(){
      weights_Grav1 <- matrix(0,nPops,nPops)
      alpha <- rnorm(3,1,0.01) #powers close to 1
      i = 1
      repeat{
        weights_Grav1[i,] <- (nHosts[i]^alpha[1])*(nHosts)^alpha[2]/(dMatrix[i,]^alpha[3])
        if(i==nPops){break}
        i = i+1
      }
      weights_Grav1[which(weights_Grav1==Inf)] <- NA
      weights <<- weights_Grav1
    },
    #Einstein's gravity model, mass*mass/distance^2
    weightsGravity = function(){
      weights_Grav2 <- matrix(0,nPops,nPops)
      i = 1
      repeat{
        weights_Grav2[i,] <- nHosts[i]*nHosts/(dMatrix[i,]^2)
        if(i==nPops){break}
        i = i+1
      }
      weights_Grav2[which(weights_Grav2==Inf)] <- NA
      weights <<- weights_Grav2
    },
    #Radiation model
    #https://en.wikipedia.org/wiki/Radiation_law_for_human_mobility
    weightsRadiation = function(){
      weights_Rad <- matrix(0,nPops,nPops)
      s <- matrix(0,nPops,nPops)
      #check <- matrix(0,nPops,nPops)
      m <- nHosts
      i = 1
      repeat{
        j = 1
        repeat{
          #check[i,j] <- check[i,j]+1
          dist <- dMatrix[i,j]
          s[i,j] <- sum(nHosts[which((dMatrix[i,]<dist) & (dMatrix[i,]!=0))])
          weights_Rad[i,j] <- ifelse(i==j,NA,(m[i]*m[j])/((m[i]+s[i,j])*(m[i]+m[j]+s[i,j])))
          if(j==nPops){break}
          j = j+1
        }
        if(i==nPops){break}
        i = i+1
      }
      weights <<- weights_Rad
    },
    # Modify weights to incorporate border effect
    # border = some fraction that defines relative probability of crossing border during migration
    # 0.5 = no effect
    # 1 = always cross border
    # 0 = never cross border
    commuterCalc = function(border=0){
      borderWeights <- matrix(NA,nPops,nPops)
      i = 1
      repeat{
        j = 1
        repeat{
          if (i!=j) {
            borderWeights[i,j] <- ifelse(nationIndex[i]==nationIndex[j],weights[i,j]*(1-border),weights[i,j]*border)
          }
          if(j==nPops){break}
          j = j+1
        }
        borderWeights[i,] <- borderWeights[i,]/sum(borderWeights[i,],na.rm=T)
        if(i==nPops){break}
        i = i+1
      }
      borderWeights <<- borderWeights

      commuters <- rbinom(nPops,nHosts,p$commRate)
      Tij <<- matrix(0,nPops,nPops)
      Tij <<- commuters*borderWeights #Average number of individuals moving from i to j in one timestep
      if(nPops==1) {Tij <<- matrix(0,1,1)}

    },
    runCommute = function(){ #runs one time step commute, simultaneously across all towns
      names <- LETTERS[1:nPops]
      commuters <- matrix(rpois(matrix(1,nPops,nPops),Tij),nPops,nPops)

      totalMoves = data.frame(whoMove=c(),whereTo=c(),whereFrom=c(),stringsAsFactors = F)
      i=1
      repeat {
        numMove <- sum(commuters[i,],na.rm=T)
        numMove <- min(numMove,length(which(town==names[i])))
        if(numMove > 0){
          whoMove <- sample(which(town==names[i]),numMove)
          probs <- replace(borderWeights[i,],which(is.na(borderWeights[i,])),0)
          whereTo <- sample(names,length(whoMove),replace=T,probs)
          whereFrom <- rep(names[i],numMove)
          moves <- data.frame(whoMove,whereTo,whereFrom,stringsAsFactors = F)
          totalMoves <- rbind(totalMoves,moves,stringsAsFactors=F)
        }
        if(i==nPops){break}
        i = i+1
      }

      moveMatrix <<- totalMoves
      town[moveMatrix$whoMove] <<- moveMatrix$whereTo
      nHosts <<- as.vector(table(town)) #recalculate population sizes
      #Needed if permanently moving, but if just commuting then home, weights etc. don't change:
        #weightsTrade() #recalculate commuting weights
        #commuterCalc() #recalculate commuter numbers/border weights
    },
    returnHome = function(){
      town[moveMatrix$whoMove] <<- moveMatrix$whereFrom
      nHosts <<- as.vector(table(town)) #recalculate population sizes
      #weightsTrade() #recalculate commuting weights
      #commuterCalc() #recalculate commuter numbers/border weights
    }
  )
)



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

# make communities 
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
# function to look at the prevs

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

# how many kids 
AgeStructure <- function(population){
  newpopages <- vector()
  for(i in 1:length(population)){
    newpopages[i] <- length(which(population[[i]]$age>=6 & population[[i]]$age<=7))
  }
 return(newpopages) 
}

# Larvae data manipulation #

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

# make pops
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






