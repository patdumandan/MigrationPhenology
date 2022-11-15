HawRd=read.csv("HawkRid.csv")
HawRd[, 3:24][is.na(HawRd[, 3:24])] <- 0

#without unknowns
hrtot=HawRd%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:19], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!year>2018) #get annual totals across all spp.

hrtot$Julian=as.POSIXlt(hrtot$Date)$yday 

#species-level####

hr_sp=hrtot%>%
  pivot_longer(cols=3:24, names_to = "Species", values_to="Count")%>%
  filter(!is.na(Count))

hr_sp=hr_sp%>%mutate(
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

hr_sp$mammals=as.numeric(hr_sp$mammals)
hr_sp$birds=as.numeric(hr_sp$birds)

hr_sp_max=hr_sp%>%group_by(Species,year,diet, birds, mammals)%>%
  filter(Count==max(Count))%>%slice_head(n=1)%>%select(-yr_tot)

hr_spp=hr_sp%>%group_by(Species,year,diet, birds, mammals)%>%summarise(yr_tot=sum(Count))

#categorical
ggplot(hr_sp_max, aes(x=year, y=Julian, col=Species))+geom_point()+facet_wrap(~diet)+
  stat_smooth(method="lm")+geom_hline(yintercept=257, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()

ggplot(hr_spp, aes(x=year, y=log(yr_tot), col=Species))+geom_point()+facet_wrap(~diet)+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")

#continuous
hr1=ggplot(hr_spp, aes(x=birds, y=log(yr_tot)))+geom_point()+xlab("% bird")+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")+ggtitle("Hawk Ridge")

hr2=ggplot(hr_spp, aes(x=mammals, y=log(yr_tot)))+geom_point()+xlab("% mammal")+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")+ggtitle("Hawk Ridge")

hr4=ggplot(hr_sp_max, aes(x=mammals, y=Julian))+geom_point()+
  stat_smooth(method="lm")+geom_hline(yintercept=257, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()+ggtitle("Hawk Ridge")+xlab("% mammal")

hr3=ggplot(hr_sp_max, aes(x=birds, y=Julian))+geom_point()+
  stat_smooth(method="lm")+geom_hline(yintercept=257, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()+ggtitle("Hawk Ridge")+xlab("% bird")

#mammals
hr_sp_mammd=hr_sp_max%>%filter(diet=="mammal")
hr_sp_mammc=hr_spp%>%filter(diet=="mammal")

hr5=ggplot(hr_sp_mammd, aes(x=year, y=Julian, col=mammals))+geom_point()+
  stat_smooth(method="lm")+geom_hline(yintercept=268, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()+ggtitle("mammal diet")

hr6=ggplot(hr_sp_mammc, aes(x=year, y=log(yr_tot), col=mammals))+geom_point()+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")+ggtitle("mammal diet")

ggarrange(hr2,hr4, hr6, hr5)

#birds
hr_sp_aved=hr_sp_max%>%filter(diet=="bird")
hr_sp_avec=hr_spp%>%filter(diet=="bird")

hr7=ggplot(hr_sp_aved, aes(x=year, y=Julian, col=birds))+geom_point()+
  stat_smooth(method="lm")+geom_hline(yintercept=257, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()+ggtitle("bird diet")

hr8=ggplot(hr_sp_avec, aes(x=year, y=log(yr_tot), col=birds))+geom_point()+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")+ggtitle("bird diet")

ggarrange(hr1,hr3, hr8, hr7)

summary(lm(hr_sp_max$Julian~hr_sp_max$birds))
