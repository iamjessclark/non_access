
# Objects for storing the different prevalences 

# make multiple high and low R0 pops ----
# burn in done by these functions too
HighPops <- MakeHigh(totalpops = 2000, populationsize = 1500) # min/max values from Irvine et al 2015 
LowPops <- MakeLow(totalpops = 3000, populationsize = 1500) # min/max values from Irvine et al 2015 

highprevs <- PrintPrevs(HighPops)
lowprevs <- PrintPrevs(LowPops)

prevs <- c(highprevs, lowprevs)

pops <- c(HighPops, LowPops)

total_months <- 144
exclusion_vec <- c(0, 0.1, 0.5)
compliance_vec <- c(0, 0.2, 1)
Values <- expand.grid(exclusion_vec, compliance_vec)

People <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=5)
PrevAll <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=total_months)
PrevAccessOnly <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=total_months)
Larvae <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=total_months)
PrevKids <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=total_months)
Prev20 <- lapply(1:length(pops), matrix, data= NA, nrow=nrow(Values), ncol=total_months)
