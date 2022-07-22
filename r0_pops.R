source("Setup3.R")

# make multiple high and low R0 pops ----

HighPops <- MakeHigh(totalpops = 100, populationsize = 1500) # min/max values from Irvine et al 2015
LowPops <- MakeLow(totalpops = 100, populationsize = 1500) # min/max values from Irvine et al 2015

highprevs <- PrintPrevs(HighPops)
lowprevs <- PrintPrevs(LowPops)

HighPrevPop <- HighPops[[25]]
LowPrevPop <- LowPops[[92]]

# run with low R0

exclusion_vec <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
compliance_vec <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
Values <- expand.grid(exclusion_vec, compliance_vec)

Peoplelowpop <- vector("list", nrow(Values))
PrevAlllowpop <- vector("list", nrow(Values))
PrevAccessOnlylowpop <- vector("list", nrow(Values))
Larvelowpop <- vector("list", nrow(Values))
PrevKidslowpop <- vector("list", nrow(Values))
Prev20lowpop <- vector("list", nrow(Values))

for(i in 1:nrow(Values)){
  exclusion <- Values[i,1]
  compliance <- Values[i,2]
  
  LowClone <- LowPrevPop$clone(deep=T)
  
  month_total <- 0
  treat_round <- 0
  
  HALT <- FALSE
  
  LowClone$assignexclusion(exclusion)
  
  while(HALT==FALSE){
    
    LowClone$runTimestep()
    month_total <-  month_total+1
    
    PrevAlllowpop[[i]][month_total] <- (length(which((LowClone$age/12)>=5 & LowClone$Mf>0))/length(which((LowClone$age/12)>=5)))*100
    PrevAccessOnlylowpop[[i]][month_total] <- (length(which(LowClone$sysExProb==0 & (LowClone$age/12)>=5 & LowClone$Mf>0))/(length(which(LowClone$sysExProb==0 & (LowClone$age/12)>=5))))*100
    Larvelowpop[[i]][month_total] <- LowClone$larvae
    #PrevKidslowpop[[i]][month_total] <- (length(which(LowClone$sysExProb==0 & (LowClone$age/12)>=6 & (LowClone$age/12)<8 & LowClone$Mf>0))/(length(which(LowClone$sysExProb==0 & (LowClone$age/12)>=6 & (LowClone$age/12)<8))))*100
    Prev20lowpop[[i]][month_total] <- (length(which(LowClone$sysExProb==0 & (LowClone$age/12)>=20 & LowClone$Mf>0))/(length(which(LowClone$sysExProb==0 & (LowClone$age/12)>=20))))*100
    
    if(month_total==12|month_total==24|month_total==36|month_total==48|month_total==60|month_total==72){
      LowClone$runMDA(pop=1,towns=NA, coverage=0.85, drug="IA", compliance=compliance)
    }
    
    
    if(month_total==12){
      Peoplelowpop[[i]] <- as.matrix(LowClone$treat,nrow=length(LowClone$treat), ncol=1)
    }
    
    if(month_total==24){
      Peoplelowpop[[i]] <- as.data.frame(cbind(Peoplelowpop[[i]], LowClone$treat))
    }
    
    if(month_total==36){
      Peoplelowpop[[i]] <- as.data.frame(cbind(Peoplelowpop[[i]], LowClone$treat))
    }
    
    if(month_total==48){
      Peoplelowpop[[i]] <- as.data.frame(cbind(Peoplelowpop[[i]], LowClone$treat))
    }
    
    if(month_total==60){
      Peoplelowpop[[i]] <- as.data.frame(cbind(Peoplelowpop[[i]], LowClone$treat))
    }
    
    if(month_total==144){
      HALT <- TRUE
    }
    
  }
  
  PrevAccessOnlylowpop[[i]] <- as.data.frame(PrevAccessOnlylowpop[[i]])
  PrevAccessOnlylowpop[[i]]$exclusion <- exclusion
  PrevAccessOnlylowpop[[i]]$compliance <- compliance
  PrevAccessOnlylowpop[[i]]$time <- 1:nrow(PrevAccessOnlylowpop[[i]])/12
  
  PrevAlllowpop[[i]] <- as.data.frame(PrevAlllowpop[[i]])
  PrevAlllowpop[[i]]$exclusion <- exclusion
  PrevAlllowpop[[i]]$compliance <- compliance
  PrevAlllowpop[[i]]$time <- 1:nrow(PrevAlllowpop[[i]])/12
  
  Larvelowpop[[i]] <- as.data.frame(Larvelowpop[[i]])
  Larvelowpop[[i]]$exclusion <- exclusion
  Larvelowpop[[i]]$compliance <- compliance
  Larvelowpop[[i]]$time <- 1:nrow(Larvelowpop[[i]])/12
  
  #PrevKidslowpop[[i]] <- as.data.frame(PrevKidslowpop[[i]])
  #PrevKidslowpop[[i]]$exclusion <- exclusion
  #PrevKidslowpop[[i]]$compliance <- compliance
  #PrevKidslowpop[[i]]$time <- 1:nrow(PrevKidslowpop[[i]])/12
  
  Prev20lowpop[[i]] <- as.data.frame(Prev20lowpop[[i]])
  Prev20lowpop[[i]]$exclusion <- exclusion
  Prev20lowpop[[i]]$compliance <- compliance
  Prev20lowpop[[i]]$time <- 1:nrow(Prev20lowpop[[i]])/12
  
  Peoplelowpop[[i]]$exclusion <- exclusion
  Peoplelowpop[[i]]$compliance <- compliance
  Peoplelowpop[[i]]$time <- 1:nrow(Peoplelowpop[[i]])/12
}

# Run with high R0

Peoplehighpop <- vector("list", nrow(Values))
PrevAllhighpop <- vector("list", nrow(Values))
PrevAccessOnlyhighpop <- vector("list", nrow(Values))
Larvehighpop <- vector("list", nrow(Values))
#PrevKidshighpop <- vector("list", nrow(Values))
Prev20highpop <- vector("list", nrow(Values))

for(i in 1:nrow(Values)){
  exclusion <- Values[i,1]
  compliance <- Values[i,2]
  
  HighClone <- HighPrevPop$clone(deep=T)
  
  month_total <- 0
  treat_round <- 0
  
  HALT <- FALSE
  
  HighClone$assignexclusion(exclusion)
  
  while(HALT==FALSE){
    
    HighClone$runTimestep()
    month_total <-  month_total+1
    
    PrevAllhighpop[[i]][month_total] <- (length(which((HighClone$age/12)>=5 & HighClone$Mf>0))/length(which((HighClone$age/12)>=5)))*100
    PrevAccessOnlyhighpop[[i]][month_total] <- (length(which(HighClone$sysExProb==0 & (HighClone$age/12)>=5 & HighClone$Mf>0))/(length(which(HighClone$sysExProb==0 & (HighClone$age/12)>=5))))*100
    Larvehighpop[[i]][month_total] <- HighClone$larvae
    #PrevKidshighpop[[i]][month_total] <- (length(which(HighClone$sysExProb==0 & (HighClone$age/12)>=6 & (HighClone$age/12)<8 & HighClone$Mf>0))/(length(which(HighClone$sysExProb==0 & (HighClone$age/12)>=6 & (HighClone$age/12)<8))))*100
    Prev20highpop[[i]][month_total] <- (length(which(HighClone$sysExProb==0 & (HighClone$age/12)>=20 & HighClone$Mf>0))/(length(which(HighClone$sysExProb==0 & (HighClone$age/12)>=20))))*100
    
    if(month_total==12|month_total==24|month_total==36|month_total==48|month_total==60|month_total==72|month_total==84|month_total==96|month_total==108|month_total==120|month_total==132|month_total==144){
      HighClone$runMDA(pop=1,towns=NA, coverage=0.85, drug="IA", compliance=compliance)
    }
    
    
    if(month_total==12){
      Peoplehighpop[[i]] <- as.matrix(HighClone$treat,nrow=length(HighClone$treat), ncol=1)
    }
    
    if(month_total==24){
      Peoplehighpop[[i]] <- as.data.frame(cbind(Peoplehighpop[[i]], HighClone$treat))
    }
    
    if(month_total==36){
      Peoplehighpop[[i]] <- as.data.frame(cbind(Peoplehighpop[[i]], HighClone$treat))
    }
    
    if(month_total==48){
      Peoplehighpop[[i]] <- as.data.frame(cbind(Peoplehighpop[[i]], HighClone$treat))
    }
    
    if(month_total==60){
      Peoplehighpop[[i]] <- as.data.frame(cbind(Peoplehighpop[[i]], HighClone$treat))
    }
    
    if(month_total==144){
      HALT <- TRUE
    }
    
  }
  
  PrevAccessOnlyhighpop[[i]] <- as.data.frame(PrevAccessOnlyhighpop[[i]])
  PrevAccessOnlyhighpop[[i]]$exclusion <- exclusion
  PrevAccessOnlyhighpop[[i]]$compliance <- compliance
  PrevAccessOnlyhighpop[[i]]$time <- 1:nrow(PrevAccessOnlyhighpop[[i]])/12
  
  PrevAllhighpop[[i]] <- as.data.frame(PrevAllhighpop[[i]])
  PrevAllhighpop[[i]]$exclusion <- exclusion
  PrevAllhighpop[[i]]$compliance <- compliance
  PrevAllhighpop[[i]]$time <- 1:nrow(PrevAllhighpop[[i]])/12
  
  Larvehighpop[[i]] <- as.data.frame(Larvehighpop[[i]])
  Larvehighpop[[i]]$exclusion <- exclusion
  Larvehighpop[[i]]$compliance <- compliance
  Larvehighpop[[i]]$time <- 1:nrow(Larvehighpop[[i]])/12
  
  #PrevKidshighpop[[i]] <- as.data.frame(PrevKidshighpop[[i]])
  #PrevKidshighpop[[i]]$exclusion <- exclusion
  #PrevKidshighpop[[i]]$compliance <- compliance
  #PrevKidshighpop[[i]]$time <- 1:nrow(PrevKidshighpop[[i]])/12
  
  Prev20highpop[[i]] <- as.data.frame(Prev20highpop[[i]])
  Prev20highpop[[i]]$exclusion <- exclusion
  Prev20highpop[[i]]$compliance <- compliance
  Prev20highpop[[i]]$time <- 1:nrow(Prev20highpop[[i]])/12
  
  Peoplehighpop[[i]]$exclusion <- exclusion
  Peoplehighpop[[i]]$compliance <- compliance
  Peoplehighpop[[i]]$time <- 1:nrow(Peoplehighpop[[i]])/12
}

#### Figures ####

for(i in 1:length(PrevAlllowpop)){
  PrevAlllowpop[[i]] <- as.data.frame(PrevAlllowpop[[i]])
  PrevAlllowpop[[i]]$exclusion <- Values[i,1]
  PrevAlllowpop[[i]]$compliance <- Values[i,2]
  PrevAlllowpop[[i]]$month <- 1:nrow(PrevAlllowpop[[i]])
  PrevAlllowpop[[i]]$population <- "Whole Population"
}

PrevAlllowpop <- do.call(rbind, PrevAlllowpop)
colnames(PrevAlllowpop)[1] <- "prevalence"

for(i in 1:length(PrevAllhighpop)){
  PrevAllhighpop[[i]] <- as.data.frame(PrevAllhighpop[[i]])
  PrevAllhighpop[[i]]$exclusion <- Values[i,1]
  PrevAllhighpop[[i]]$compliance <- Values[i,2]
  PrevAllhighpop[[i]]$month <- 1:nrow(PrevAllhighpop[[i]])
  PrevAllhighpop[[i]]$population <- "Whole Population"
}

PrevAllhighpop <- do.call(rbind, PrevAllhighpop)
colnames(PrevAllhighpop)[1] <- "prevalence"

for(i in 1:length(PrevAccessOnlylowpop)){
  PrevAccessOnlylowpop[[i]] <- as.data.frame(PrevAccessOnlylowpop[[i]])
  PrevAccessOnlylowpop[[i]]$exclusion <- Values[i,1]
  PrevAccessOnlylowpop[[i]]$compliance <- Values[i,2]
  PrevAccessOnlylowpop[[i]]$month <- 1:nrow(PrevAccessOnlylowpop[[i]])
  PrevAccessOnlylowpop[[i]]$population <- "Access Only"
}

PrevAccessOnlylowpop <- do.call(rbind, PrevAccessOnlylowpop)
colnames(PrevAccessOnlylowpop)[1] <- "prevalence"

for(i in 1:length(PrevAccessOnlyhighpop)){
  PrevAccessOnlyhighpop[[i]] <- as.data.frame(PrevAccessOnlyhighpop[[i]])
  PrevAccessOnlyhighpop[[i]]$exclusion <- Values[i,1]
  PrevAccessOnlyhighpop[[i]]$compliance <- Values[i,2]
  PrevAccessOnlyhighpop[[i]]$month <- 1:nrow(PrevAccessOnlyhighpop[[i]])
  PrevAccessOnlyhighpop[[i]]$population <- "Access Only"
}

PrevAccessOnlyhighpop <- do.call(rbind, PrevAccessOnlyhighpop)
colnames(PrevAccessOnlyhighpop)[1] <- "prevalence"

#for(i in 1:length(PrevKidslowpop)){
  #PrevKidslowpop[[i]] <- as.data.frame(PrevKidslowpop[[i]])
  #PrevKidslowpop[[i]]$exclusion <- Values[i,1]
  #PrevKidslowpop[[i]]$compliance <- Values[i,2]
  #PrevKidslowpop[[i]]$month <- 1:nrow(PrevKidslowpop[[i]])
  #PrevKidslowpop[[i]]$population <- "6 & 7 yearolds"
#}

#PrevKidslowpop <- do.call(rbind, PrevKidslowpop)
#colnames(PrevKidslowpop)[1] <- "prevalence"

for(i in 1:length(Prev20lowpop)){
  Prev20lowpop[[i]] <- as.data.frame(Prev20lowpop[[i]])
  Prev20lowpop[[i]]$exclusion <- Values[i,1]
  Prev20lowpop[[i]]$compliance <- Values[i,2]
  Prev20lowpop[[i]]$month <- 1:nrow(Prev20lowpop[[i]])
  Prev20lowpop[[i]]$population <- ">20 yearolds"
}

Prev20lowpop <- do.call(rbind, Prev20lowpop)
colnames(Prev20lowpop)[1] <- "prevalence"

#for(i in 1:length(PrevKidshighpop)){
  #PrevKidshighpop[[i]] <- as.data.frame(PrevKidshighpop[[i]])
  #PrevKidshighpop[[i]]$exclusion <- Values[i,1]
  #PrevKidshighpop[[i]]$compliance <- Values[i,2]
  #PrevKidshighpop[[i]]$month <- 1:nrow(PrevKidshighpop[[i]])
  #PrevKidshighpop[[i]]$population <- "6 & 7 yearolds"
#}

#PrevKidshighpop <- do.call(rbind, PrevKidshighpop)
#colnames(PrevKidshighpop)[1] <- "prevalence"

for(i in 1:length(Prev20highpop)){
  Prev20highpop[[i]] <- as.data.frame(Prev20highpop[[i]])
  Prev20highpop[[i]]$exclusion <- Values[i,1]
  Prev20highpop[[i]]$compliance <- Values[i,2]
  Prev20highpop[[i]]$month <- 1:nrow(Prev20highpop[[i]])
  Prev20highpop[[i]]$population <- ">20 yearolds"
}

Prev20highpop <- do.call(rbind, Prev20highpop)
colnames(Prev20highpop)[1] <- "prevalence"

prev67lowpop <- bind_rows(PrevAlllowpop, PrevAccessOnlylowpop, Prev20lowpop) #PrevKidslowpop)
prev67lowpop$exclusion[which(prev67lowpop$exclusion==0.0)] <- "0%"
prev67lowpop$exclusion[which(prev67lowpop$exclusion==0.1)] <- "10%"
prev67lowpop$exclusion[which(prev67lowpop$exclusion==0.5)] <- "50%"

prev67highpop <- bind_rows(PrevAllhighpop, PrevAccessOnlyhighpop, Prev20highpop) #PrevKidshighpop)
prev67highpop$exclusion[which(prev67highpop$exclusion==0.0)] <- "0%"
prev67highpop$exclusion[which(prev67highpop$exclusion==0.1)] <- "10%"
prev67highpop$exclusion[which(prev67highpop$exclusion==0.5)] <- "50%"

coul <- viridis(100)

lowpopplot <- prev67lowpop %>% 
  filter(compliance==0.2 & exclusion=="0%" |compliance==0.2 & exclusion=="10%" | compliance==0.2 & exclusion=="50%")%>%
  ggplot()+
  geom_line(aes(x=month/12, y=prevalence, group=population, colour=population))+
  facet_grid(exclusion~.)+
  theme_bw()+
  geom_hline(yintercept = 1, colour="red", alpha=0.3)+
  annotate("rect", xmin=5.5, xmax=12, ymin=0, ymax=7, alpha=0.3)+
  xlab("years")+
  scale_colour_manual(values=coul[c(1, 40, 84)])+
  scale_x_continuous(breaks = c(0:12))+
  scale_y_continuous(limits=c(0,7), breaks=seq(0,7, by=2))+
  theme(legend.position = "top", legend.title = element_text(colour="white"), text=element_text(size=14))+
  guides(color = guide_legend(override.aes = list(size = 3) ) )


highpopplot <- prev67highpop %>% 
  filter(compliance==0.2 & exclusion=="0%" |compliance==0.2 & exclusion=="10%" | compliance==0.2 & exclusion=="50%")%>%
  ggplot()+
  geom_line(aes(x=month/12, y=prevalence, group=population, colour=population))+
  facet_grid(exclusion~.)+
  theme_bw()+
  geom_hline(yintercept = 1, colour="red", alpha=0.3)+
  #annotate("rect", xmin=5.5, xmax=8, ymin=0, ymax=20, alpha=0.3)+
  xlab("years")+
  scale_colour_manual(values=coul[c(1, 40, 84)])+
  scale_x_continuous(breaks = c(0:12))+
  scale_y_continuous(limits=c(0,25),breaks=seq(0,26, by=5))+
  theme(legend.position = "none", text=element_text(size=14))+
  guides(color = guide_legend(override.aes = list(size = 3) ) )

prevplots <- plot_grid(highpopplot, lowpopplot, ncol=1, labels=c("A", "B"))

# larvae plots
Larvae <- LarvaeData(Larvelowpop, Larvehighpop)

coul <- magma(100)

Larvae %>%
  filter(compliance==0.2 &exclusion=="0%" & time<=5  |compliance==0.2& exclusion=="10%"& time<=5  |compliance==0.2&exclusion=="50%" & time<=5 ) %>%
  mutate(population=fct_recode(population, `high transmission`="high prevalence", `low transmission`="low prevalence"))%>%
  ggplot()+
  geom_line(aes(x=time, y=larvae, group=exclusion, colour=exclusion, linetype=exclusion), size=0.8)+
  facet_wrap(.~population, ncol=1, nrow=2)+
  xlab("years")+
  theme_bw()+
  theme(strip.text = element_text(size = 16), 
        axis.text = element_text( size = 14), 
        axis.title = element_text( size = 14))+
  scale_linetype_manual(values=c("solid", "dotted", 'dashed'))+
  scale_color_manual(values=c("#0A6777","#D8661B",  "#07FFBD"))+
  #coul[c(50,87, 35)]
  ylab("larval density")
ggsave("twopoplarvae.pdf")

# people treated plots

for(i in 1:length(Peoplelowpop)){
  colnames(Peoplelowpop[[i]]) <- c("1", "2", "3", "4", "5", "exclusion", "compliance", "time")
  Peoplelowpop[[i]]$NoTrts <- as.factor(rowSums(Peoplelowpop[[i]][1:5]==TRUE))
}

TrtRndslowpop <- list()
for(i in 1:length(Peoplelowpop)){
  colnames(Peoplelowpop[[i]]) <- c("1", "2", "3", "4", "5", "exclusion", "compliance", "time", "NoTrts")
  TrtRndslowpop[[i]] <- Peoplelowpop[[i]] %>%
    group_by(NoTrts)%>%
    summarise(counts=100*(n()/nrow(Peoplelowpop[[i]])))%>%
    pivot_wider(names_from = NoTrts, values_from = counts)
  TrtRndslowpop[[i]] <- as.data.frame(TrtRndslowpop[[i]])
}

NoRndslowpop <- rbindlist(list(TrtRndslowpop[[1]],TrtRndslowpop[[2]]), fill = TRUE)
for(i in 3:length(TrtRndslowpop)){
  NoRndslowpop <- rbindlist(list(NoRndslowpop, TrtRndslowpop[[i]]), fill = TRUE)
}

NoRndslowpop <- bind_cols(Values, NoRndslowpop)
colnames(NoRndslowpop) <- c("exclusion", "compliance", "0", "1", "2", "3", "4", "5")

NoRndslowpop <- NoRndslowpop %>%
  pivot_longer(cols=c(`0`:`5`), names_to = "Treatments", values_to = "people")

for(i in 1:length(Peoplehighpop)){
  colnames(Peoplehighpop[[i]]) <- c("1", "2", "3", "4", "5", "exclusion", "compliance", "time")
  Peoplehighpop[[i]]$NoTrts <- as.factor(rowSums(Peoplehighpop[[i]][1:5]==TRUE))
}

TrtRndshighpop <- list()
for(i in 1:length(Peoplehighpop)){
  colnames(Peoplehighpop[[i]]) <- c("1", "2", "3", "4", "5", "exclusion", "compliance", "time", "NoTrts")
  TrtRndshighpop[[i]] <- Peoplehighpop[[i]] %>%
    group_by(NoTrts)%>%
    summarise(counts=100*(n()/nrow(Peoplehighpop[[i]])))%>%
    pivot_wider(names_from = NoTrts, values_from = counts)
  TrtRndshighpop[[i]] <- as.data.frame(TrtRndshighpop[[i]])
}

NoRndshighpop <- rbindlist(list(TrtRndshighpop[[1]],TrtRndshighpop[[2]]), fill = TRUE)
for(i in 3:length(TrtRndshighpop)){
  NoRndshighpop <- rbindlist(list(NoRndshighpop, TrtRndshighpop[[i]]), fill = TRUE)
}

NoRndshighpop <- bind_cols(Values, NoRndshighpop)
colnames(NoRndshighpop) <- c("exclusion", "compliance", "0", "1", "2", "3", "4", "5")

NoRndshighpop <- NoRndshighpop %>%
  pivot_longer(cols=c(`0`:`5`), names_to = "Treatments", values_to = "people")

NoRndslowpop$exclusion[which(NoRndslowpop$exclusion==0)] <- "0%"
NoRndslowpop$exclusion[which(NoRndslowpop$exclusion==0.1)] <- "10%"
NoRndslowpop$exclusion[which(NoRndslowpop$exclusion==0.2)] <- "20%"
NoRndslowpop$exclusion[which(NoRndslowpop$exclusion==0.3)] <- "30%"
NoRndslowpop$exclusion[which(NoRndslowpop$exclusion==0.4)] <- "40%"
NoRndslowpop$exclusion[which(NoRndslowpop$exclusion==0.5)] <- "50%"

NoRndslowpop %>% 
  filter(compliance==0 & exclusion=="0%"|compliance==0 & exclusion=="10%"|compliance==0 & exclusion=="50%"|
           compliance==0.4 & exclusion=="0%"|compliance==0.4 & exclusion=="10%"|compliance==0.4 & exclusion=="50%"|
           compliance==1 & exclusion=="0%"|compliance==1 & exclusion=="10%"|compliance==1 & exclusion=="50%")%>%
  ggplot()+
  geom_col(aes(x=Treatments, y=people), fill="white", colour="black")+
  facet_grid(exclusion~compliance)+
  ylab("percentage treated")+
  theme_bw()
ggsave("peopletreatedlowpop.pdf")

NoRndshighpop$exclusion[which(NoRndshighpop$exclusion==0)] <- "0%"
NoRndshighpop$exclusion[which(NoRndshighpop$exclusion==0.1)] <- "10%"
NoRndshighpop$exclusion[which(NoRndshighpop$exclusion==0.2)] <- "20%"
NoRndshighpop$exclusion[which(NoRndshighpop$exclusion==0.3)] <- "30%"
NoRndshighpop$exclusion[which(NoRndshighpop$exclusion==0.4)] <- "40%"
NoRndshighpop$exclusion[which(NoRndshighpop$exclusion==0.5)] <- "50%"

NoRndshighpop %>%
  filter(compliance==0 & exclusion=="0%"|compliance==0 & exclusion=="10%"|compliance==0 & exclusion=="50%"|
           compliance==0.4 & exclusion=="0%"|compliance==0.4 & exclusion=="10%"|compliance==0.4 & exclusion=="50%"|
           compliance==1 & exclusion=="0%"|compliance==1 & exclusion=="10%"|compliance==1 & exclusion=="50%")%>%
  ggplot()+
  geom_col(aes(x=Treatments, y=people), fill="white", colour="black")+
  facet_grid(exclusion~compliance)+
  ylab("percentage treated")+
  theme_bw()
ggsave("peopletreatedhighpop.pdf")


#### access x adherence ####
# Population 1

exclusion_vec <- seq(from=0, to = 0.5, length.out=100)
compliance_vec <- seq(from = 0, to = 1, length.out=100)

datalowpop <- expand.grid(exclusion_vec, compliance_vec)
colnames(datalowpop) <- c("exclusion", "compliance")
datalowpop$prev_observed <- NA
datalowpop$prev_all <- NA

for(i in 1:nrow(datalowpop)){
  exclusion <- datalowpop[i,1]
  compliance <- datalowpop[i,2]
  
  LowClone <- LowPrevPop$clone(deep=T)
  
  mfprevPreTAS <- NA
  mfprevPreTAS_no_exclusion <- NA
  
  month_total <- 0
  treat_round <- 0
  
  mfprevPreTAS_no_exclusion <- NA
  mfprevTAS_notexcluded <- NA
  
  HALT <- FALSE
  
  LowClone$assignexclusion(exclusion)
  
  while(HALT==FALSE){
    
    LowClone$runTimestep()
    month_total <-  month_total+1
    
    if(month_total==12|month_total==24|month_total==36|month_total==48|month_total==60){
      LowClone$runMDA(pop=1,towns=NA, coverage=0.85, drug="IA", compliance=compliance)
    }
    
    
    if(month_total==66){
      datalowpop[i,3] <- (length(which(LowClone$sysExProb==0 & (LowClone$age/12)>=5 & LowClone$Mf>0))/(length(which(LowClone$sysExProb==0 & (LowClone$age/12)>=5))))*100
      datalowpop[i,4] <- (length(which((LowClone$age/12)>=5 & LowClone$Mf>0))/length(which((LowClone$age/12)>=5)))*100
    }                  
    
    if(month_total==66){
      HALT <- TRUE
    }
    
  }
}

Differencelowpop <- datalowpop %>%
  mutate(delta_prev=prev_all-prev_observed, 
         population="low prevalence")

# Population 2

datahighpop <- expand.grid(exclusion_vec, compliance_vec)
colnames(datahighpop) <- c("exclusion", "compliance")
datahighpop$prev_observed <- NA
datahighpop$prev_all <- NA

for(i in 1:nrow(datahighpop)){
  exclusion <- datahighpop[i,1]
  compliance <- datahighpop[i,2]
  
  HighClone <- HighPrevPop$clone(deep=T)
  
  mfprevPreTAS <- NA
  mfprevPreTAS_no_exclusion <- NA
  
  month_total <- 0
  treat_round <- 0
  
  mfprevPreTAS_no_exclusion <- NA
  mfprevTAS_notexcluded <- NA
  
  HALT <- FALSE
  
  HighClone$assignexclusion(exclusion)
  
  while(HALT==FALSE){
    
    HighClone$runTimestep()
    month_total <-  month_total+1
    
    if(month_total==12|month_total==24|month_total==36|month_total==48|month_total==60){
      HighClone$runMDA(pop=1,towns=NA, coverage=0.85, drug="IA", compliance=compliance)
    }
    
    
    if(month_total==66){
      datahighpop[i,3] <- (length(which(HighClone$sysExProb==0 & (HighClone$age/12)>=5 & HighClone$Mf>0))/(length(which(HighClone$sysExProb==0 & (HighClone$age/12)>=5))))*100
      datahighpop[i,4] <- (length(which((HighClone$age/12)>=5 & HighClone$Mf>0))/length(which((HighClone$age/12)>=5)))*100
    }
    
    if(month_total==66){
      HALT <- TRUE
    }
    
  }
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken  


### observ vs true prev Figures ####

Differencehighpop <- datahighpop %>%
  mutate(delta_prev=prev_all-prev_observed, 
         population="high prevalence")
coul <- magma(100)

deltaplot <- bind_rows(Differencelowpop, Differencehighpop)%>%
  #pivot_longer(cols=c(prev_all,prev_observed), values_to = "prevalence", names_to="obvs_pop")%>%
  ggplot(aes(x=compliance, y=exclusion*100, fill=delta_prev))+
  geom_tile()+
  ylab("Exclusion percentage")+
  facet_grid(population~.)+
  theme_cowplot()+
  scale_fill_viridis(option = "magma")+
  xlab("Adherence")+ ylab("Non-access percentage")+
  #labs(fill="\u0394 \n prevalence")+
  theme(legend.title = element_text(colour="white"))

figure2 <- plot_grid(prevplots, deltaplot, ncol=2)
ggsave("figure2panel.pdf")

save.image(file='myEnvironment.RData')

Larvae %>%
  filter(population=="low prevalence") %>%
  bind_cols(PrevAlllowpop)%>%
  ggplot()+
  geom_line(aes(x=larvae, y=prevalence))+
  facet_grid(`compliance...3` ~ `exclusion...2`)







