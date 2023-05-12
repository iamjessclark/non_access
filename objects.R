
# Objects for storing the different prevalences 

# make list of populations that have a prevalence >1%
pops <- makerevisionpops()
#save.image(file='myEnvironment.RData')
# get their prevalences - anyone over 5 with mf >1 
popsprevs <- PrintPrevs(pops)
# check prevalence spread
hist(popsprevs)
# reduce over representation of certain prevalences
which(popsprevs >=20 & popsprevs <= 25)
# randomly sample some to remove
sample(which(popsprevs >=20 & popsprevs <= 25) ,size = 1300, replace = FALSE)
# make a dummy list to ensure no mistakes - takes a while to make pops
pops2 <- pops
# remove a subset of the populations with over represented prevalence 
pops2 <- pops2[-sample(which(popsprevs >=20 & popsprevs <= 25) ,size = 1300,replace = TRUE)]
# check prevalence
popsprevs2 <- PrintPrevs(pops2)
hist(popsprevs2)
# save
pops <- pops2
popsprevs <- PrintPrevs(pops)
rm(pops2)
rm(popsprevs2)

# save ks and VtH 
ksandVtH <- ksVtH(pops)
ksandVtH$prevs <- popsprevs

# objects for simulations 
total_months <- 144
exclusion_vec <- c(0, 0.1, 0.5)
compliance_vec <- c(0, 0.2, 1)
Values <- expand.grid(exclusion_vec, compliance_vec)

# make objects to store model runs 
People <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=5)
PrevAll <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=total_months)
PrevAccessOnly <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=total_months)
#Larvae <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=total_months)
PrevKids <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=total_months)
Prev20 <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=total_months)

#--------------------------------------------------------------------------------------------------------
# used pre-revision 2 
#make multiple high and low R0 pops ----
# burn in done by these functions too
#HighPops <- MakeHigh(totalpops = 5000, populationsize = 1500) # min/max values from Irvine et al 2015 
#LowPops <- MakeLow(totalpops = 5000, populationsize = 1500) # min/max values from Irvine et al 2015 
#highprevs <- PrintPrevs(HighPops)
#lowprevs <- PrintPrevs(LowPops)

#prevs <- c(highprevs, lowprevs)
#prevs <- c(group1prevs, group2prevs, group3prevs, group4prevs, group5prevs)

#pops <- c(HighPops, LowPops)
#pops <- c(Group1, Group2, Group3, Group4, Group5)

#pops <- pops[-which(prevs==0)]
#prevs <- prevs[-which(prevs==0)]
