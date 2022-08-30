#HMS####
hm=hms%>%
 pivot_longer(cols=6:30, names_to = "Species", values_to="Count")%>%
  filter(!is.na(Count))%>%
  mutate(date=make_date(YR, MO, DAY))

hm1=hm%>%group_by(YR, Species)%>%mutate(cum_tot=cumsum(Count), #get cumulative sum for each year
                                  cum_pt=cum_tot/sum(Count),#calculate % of daily count/annual count
                                  MPD=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR, Species)%>%filter(MPD==1, !(YR<1983))%>%slice_head(n=1)  #select rows where cum_pt is >50%

hm1=hm1%>%filter(!Species%in%c("UNID.ACCIPITER", "UNID.BUTEO","UNID.EAGLE","UNID.FALCON","UNID.RAPTOR"))%>%
  mutate(troph=case_when(Species=="AMKE"~"mammal",
                         Species=="BAEA"~"mammal",
                         Species=="BLVU"~"scavenger",
                         Species=="BWHA"~"mammal",
                         Species=="COHA"~"avian",
                         Species=="GOEA"~"mammal",
                         Species=="MERL"~"avian",
                         Species=="NOGO"~"avian",
                         Species=="NOHA"~"mammal",
                         Species=="OSPR"~"mammal",
                         Species=="PEFA"~"avian",
                         Species=="RTHA"~"mammal",
                         Species=="RSHA"~"mammal",
                         Species=="RLHA"~"mammal",
                         Species=="SSHA"~"avian",
                         Species=="TUVU"~"scavenger"),
         size=case_when(Species=="AMKE"~"small",
                        Species=="BAEA"~"large",
                        Species=="BLVU"~"large",
                        Species=="BWHA"~"small",
                        Species=="COHA"~"small",
                        Species=="GOEA"~"large",
                        Species=="MERL"~"small",
                        Species=="NOGO"~"mid",
                        Species=="NOHA"~"small",
                        Species=="OSPR"~"large",
                        Species=="PEFA"~"mid",
                        Species=="RTHA"~"mid",
                        Species=="RSHA"~"mid",
                        Species=="RLHA"~"mid",
                        Species=="SSHA"~"small",
                        Species=="TUVU"~"large"))%>%na.omit()
ggplot(hm1, aes(y=Julian, x=YR, col=Species))+geom_point()+facet_wrap(~size+troph)+stat_smooth(method="lm")+ggtitle("Hawk Mt (east)")

#Hawk Ridge
