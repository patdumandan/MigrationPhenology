require(tidyr)
require(ggridges)

mig_window1=gos%>%group_by(Species, year)%>%
  filter(Count>0)%>%
  filter(Date==min(Date))%>%
  rename("first_obs"="Date")%>%
  filter(!Species%in%c("UA", "UB", "UF", "UE", "UR"))

mig_window2=gos%>%group_by(Species, year)%>%
  filter(Count>0)%>%
  filter(Date==max(Date))%>%
  rename("last_obs"="Date")%>%
  filter(!Species%in%c("UA", "UB", "UF", "UE", "UR"))

mig_dat=left_join(mig_window1,mig_window2, by=c("Species", "year"))
mig_dat$first_JD=as.POSIXlt(mig_dat$first_obs)$yday 
mig_dat$last_JD=as.POSIXlt(mig_dat$last_obs)$yday 
mig_dat=mig_dat%>%
  mutate(mig_window=last_JD-first_JD)

mig_dat=mig_dat%>%mutate(
                 diet=case_when(Species=="AK" ~"bird",
                                Species=="BE" ~"other",
                                Species=="CH" ~"bird",
                                Species=="BV" ~ "other",
                                Species=="TV" ~"other",
                                Species=="OS" ~"other",
                                Species== "NH" ~"mammal",
                                Species=="SS" ~"bird",
                                Species=="NG"~ "bird",
                                Species=="RS" ~"mammal",
                                Species=="BW"~"mammal",
                                Species=="RT" ~"mammal",
                                Species=="RL" ~"mammal",
                                Species=="SW"~"mammal",
                                Species=="FH"~"mammal",
                                Species=="GE"~"mammal",
                                Species=="ML"~"bird",
                                Species=="PG"~"bird",
                                Species=="PR"~"mammal",
                                Species=="MK"~"other"),
                 mammals=case_when(Species=="AK" ~"34.6",
                                   Species=="BE" ~"16",
                                   Species=="CH" ~"23.3",
                                   Species=="BV" ~ "92.5",
                                   Species=="TV" ~"79.7",
                                   Species=="OS" ~"0",
                                   Species== "NH" ~"48.2",
                                   Species=="SS" ~"2.72",
                                   Species=="NG"~ "24.7",
                                   Species=="RS" ~"39.5",
                                   Species=="BW"~"42.2",
                                   Species=="RT" ~"73.8",
                                   Species=="RL" ~"85.3",
                                   Species=="SW"~"48",
                                   Species=="FH"~"81.9",
                                   Species=="GE"~"80.6",
                                   Species=="ML"~"1.28",
                                   Species=="PG"~"5.11",
                                   Species=="PR"~"42.3",
                                   Species=="MK"~"2.22"),
                 birds=case_when(Species=="AK" ~"13.5",
                                 Species=="BE" ~"23.4",
                                 Species=="CH" ~"71.4",
                                 Species=="BV" ~ "7.5",
                                 Species=="TV" ~"6.42",
                                 Species=="OS" ~"0",
                                 Species== "NH" ~"38.6",
                                 Species=="SS" ~"89.4",
                                 Species=="NG"~ "74.2",
                                 Species=="RS" ~"6.43",
                                 Species=="BW"~"16.5",
                                 Species=="RT" ~"17",
                                 Species=="RL" ~"9.02",
                                 Species=="SW"~"10.6",
                                 Species=="FH"~"9.89",
                                 Species=="GE"~"17.8",
                                 Species=="ML"~"90.3",
                                 Species=="PG"~"93.5",
                                 Species=="PR"~"45",
                                 Species=="MK"~"0"))%>%
  filter(!Species%in%c("UA", "UB", "UF", "UE", "UR"))

bird_dat=mig_dat%>%filter(diet=="bird")
mam_dat=mig_dat%>%filter(diet=="mammal")

ggplot(mig_dat, aes(y=as.factor(Species), x=first_JD, col=Species))+geom_density_ridges2()+facet_wrap(~diet)
ggplot(mig_dat, aes(y=as.factor(Species), x=last_JD, col=Species))+geom_density_ridges2()+facet_wrap(~diet)
ggplot(mig_dat, aes(y=as.factor(Species), x=mig_window, col=Species))+geom_density_ridges2()+facet_wrap(~diet)

ggplot(mig_dat, aes(x=year, y=first_JD, col=Species))+geom_point()+facet_wrap(~diet)+stat_smooth(method="lm")+
  ylab("first observation date(Julian)")

ggplot(mig_dat, aes(x=year, y=last_JD, col=Species))+geom_point()+facet_wrap(~diet)+stat_smooth(method="lm")+
  ylab("last observation date(Julian)")

t1=gos%>%filter(Species=="AK")%>%
  filter(Count>0)%>%
  group_by(year)%>%
  filter(Date==min(Date))%>%
  rename("first_obs"="Date")

t2=gos%>%filter(Species=="AK")%>%
  filter(Count>0)%>%
  group_by(year)%>%
  filter(Date==max(Date))%>%
  rename("last_obs"="Date")

ak_dat=left_join(t1,t2, by=c("Species", "year"))
ak_dat$first_JD=as.POSIXlt(ak_dat$first_obs)$yday 
ak_dat$last_JD=as.POSIXlt(ak_dat$last_obs)$yday 
ak_dat=ak_dat%>%
  mutate(mig_window=last_JD-first_JD)

rt1=gos%>%filter(Species=="RT")
rt1$JD=as.POSIXlt(rt1$Date)$yday 

AK1=gos%>%filter(Species=="AK")
AK1$JD=as.POSIXlt(AK1$Date)$yday 
