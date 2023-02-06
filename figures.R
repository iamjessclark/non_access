
require(ggsci)
source("Access MS Revisions.R")

#### probability clearance ####

ProbAll_clear <- prob_clear_list_func(PrevAll, 2, Values, prevs) %>%
  mutate(observations = "true prevalence")
ProbAccess_clear <- prob_clear_list_func(PrevAccessOnly, 2, Values, prevs)%>%
  mutate(observations = "community-wide with access")
ProbKids_clear <- prob_clear_list_func(PrevKids, 2, Values, prevs)%>%
  mutate(observations = "TAS")
Prob20_clear <- prob_clear_list_func(Prev20, 1, Values, prevs) %>%
  mutate(observations = "mf TAS (>20-year-olds)")

clear_df <- bind_rows(ProbAll_clear, ProbAccess_clear, ProbKids_clear, Prob20_clear)

clear_df %>%
  filter(adherence == 0.2, 
         starting_prev == "0-4%" | starting_prev == "05-9%" | starting_prev == "10-14%"|
         starting_prev == "15-19%" | starting_prev == "20-24%" | starting_prev == "25-29%"|
           starting_prev == "30-34%") %>%
  mutate(observations = factor(observations, levels = c("true prevalence" , "community-wide with access", "TAS", "mf TAS (>20-year-olds)")),
         exclusion = as_factor(exclusion), 
         exclusion = fct_recode(exclusion, "0%" = "0", "10%" = "0.1", "50%" = "0.5")) %>%
  ggplot()+
  geom_point(aes(x=starting_prev, y = prob_clear, colour = exclusion), size = 3)+
  geom_line(aes(x=starting_prev, y = prob_clear, group = exclusion))+
  facet_grid(observations ~ .) +
  theme_bw()+
  xlab("baseline prevalence") + ylab("Probability achieving threshold")+
  scale_colour_futurama()
ggsave("threshold probability.pdf")  


#### dynamics pre and post treatment ####

PrevAll_dynamics <- prev_dynnamics_func(PrevAll, Values, 2, prevs)
PrevAll_dynamics$observed <- "true prevalence"

PrevKids_dynamics <- prev_dynnamics_func(PrevKids, Values, 2, prevs) 
PrevKids_dynamics$observed <-  "TAS"

PrevAccess_dynamics <- prev_dynnamics_func(PrevAccessOnly, Values, 2, prevs)
PrevAccess_dynamics$observed <- "community-wide with access"

Prev20_dynamics <- prev_dynnamics_func(Prev20, Values, 1, prevs)
Prev20_dynamics$observed <- "mf TAS (>20-year-olds)"

dynamicsdf <- bind_rows(PrevAll_dynamics, PrevKids_dynamics, PrevAccess_dynamics, Prev20_dynamics)

  

dynamicsdf %>%
  mutate(exclusion = as_factor(exclusion), 
         exclusion = fct_recode(exclusion, "0%" = "0", "10%" = "0.1", "50%" = "0.5"),
         observed = factor(observed, levels = c("true prevalence" , "community-wide with access", "TAS", "mf TAS (>20-year-olds)"))) %>%
  group_by(exclusion, time, observed) %>%
  summarise(mean = mean(prevalence), sd=sd(prevalence), 
            ribbonmin = mean-sd, ribbonmax = mean + sd, ribbonmin = if_else(ribbonmin < 0, 0, ribbonmin)) %>% 
  ggplot()+
  geom_line(aes(x = as.integer(as.character(time)), y = mean, colour = exclusion)) + 
  geom_ribbon(aes(x = as.integer(as.character(time)), ymin=ribbonmin, ymax = ribbonmax, fill = exclusion), alpha = 0.1)+
  scale_y_continuous(breaks = seq(0, 10, by = 2))+
  facet_grid(observed~.,  scales="free")+
  scale_x_continuous(breaks=c(1:12))+
  theme_bw() +
  scale_colour_futurama()+
  ylab("prevalence") + xlab("years")
  ggsave("dynamics.pdf")
  
#### Number of treatment rounds ####
    
MDArounds_Alldf <- Treatment_Rounds(MDArounds_All, observation = "true prevalence", prevs = prevs, Values = Values)
  
MDArounds_accessdf <- Treatment_Rounds(MDArounds_Access, observation = "community-wide with access", prevs = prevs, Values = Values)
MDArounds_kidsdf <- Treatment_Rounds(MDArounds_kids, observation = "TAS", prevs = prevs, Values = Values)
MDArounds_mf20df <- Treatment_Rounds(MDArounds_mf20, observation = "mf TAS (>20-year-olds)", prevs = prevs, Values = Values)

roundsdf <- bind_rows(MDArounds_Alldf, MDArounds_accessdf, MDArounds_kidsdf, MDArounds_mf20df)

roundsdf %>%
  mutate(exclusion = as_factor(exclusion), 
         exclusion = fct_recode(exclusion, "0%" = "0", "10%" = "0.1", "50%" = "0.5"), 
         observed = factor(observed, levels = c("true prevalence" , "community-wide with access", "TAS", "mf TAS (>20-year-olds)"))) %>%
  filter(adherence == 0.2, 
         !is.na(starting_prev), 
         starting_prev == "0-4%" | starting_prev == "05-9%" | starting_prev == "10-14%"|
         starting_prev == "15-19%" | starting_prev == "20-24%" | starting_prev == "25-29%"|
         starting_prev == "30-34%") %>%
  ggplot()+
  geom_violin(aes(x = starting_prev, y = rounds, fill = exclusion, colour = exclusion), alpha = 0.8)+
  geom_jitter(aes(x = starting_prev, y = rounds, colour = exclusion),  alpha = 0.1) +
  facet_grid(observed~exclusion)+
  scale_fill_futurama()+
  scale_colour_futurama()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1.01, hjust=1))+
  xlab("baseline prevalence")+ylab("treatment rounds")
ggsave("treatrounds.pdf")


roundsdf %>% 
  mutate(exclusion = as_factor(exclusion), 
         exclusion = fct_recode(exclusion, "0%" = "0", "10%" = "0.1", "50%" = "0.5"), 
         observed = factor(observed, levels = c("true prevalence" , "community-wide with access", "TAS", "mf TAS (>20-year-olds)"))) %>%
  filter(adherence == 0.2, rounds <= 20, starting >= 5 & starting < 20,
         !is.na(starting_prev), 
         starting_prev == "0-4%" | starting_prev == "05-9%" | starting_prev == "10-14%"|
           starting_prev == "15-19%" | starting_prev == "20-24%" | starting_prev == "25-29%"|
           starting_prev == "30-34%") %>%
  ggplot()+
  geom_point(aes(x = starting, y = rounds, colour = exclusion), alpha = 0.45)+
  facet_grid(observed~exclusion)+
  scale_colour_futurama()+
  theme_bw()+
  xlab("baseline prevalence") + ylab("treatment rounds")
ggsave("treatroundsscatter.pdf")

#### larvae ####

Larvaedf <- larvaefunc(Larvae) 
Larvaedf %>% 
  mutate(exclusion = as_factor(exclusion), 
         exclusion = fct_recode(exclusion, "0%" = "0", "10%" = "0.1", "50%" = "0.5"))%>%
  group_by(exclusion, Time) %>%
  #summarise(mean = mean(l3density), sd= sd(l3density), ymin = mean-sd, ymax = mean + sd)%>%
  summarise(l3density = quantile(l3density, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75)) %>%
  pivot_wider(names_from = q, values_from = l3density) %>%
  ggplot()+
  geom_line(aes(x=as.numeric(as.character(Time)), y=`0.5`, group=exclusion, colour=exclusion, linetype=exclusion), size=0.8)+
  geom_ribbon(aes(x = as.numeric(as.character(Time)), ymin = `0.25`, ymax = `0.75`, fill = exclusion), alpha = 0.2)+
  facet_grid(exclusion ~.) +
  scale_x_continuous(breaks = c(1:12))+
  scale_colour_futurama()+
  scale_fill_futurama()+
  theme_bw()+
  xlab("years")+ylab("average mean free stage larvae")
ggsave("larvae_rev1.pdf")
