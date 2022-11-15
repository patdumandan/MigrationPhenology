require(dplyr)
require(tidyr)
require(lubridate)
require(ggplot2)
require(vegan)
#data manipulation####
GosMts=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/GosMts.csv")
GosMts[, 3:27][is.na(GosMts[, 3:27])] <- 0
GosMts$Julian=as.POSIXlt(GosMts$Date)$yday 
GosMts=GosMts%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:27], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!year>2018, !Julian>308)

#select only 1983-2018 (Aug 15-Nov 5)

#trend in max migration dates
gostotmax=GosMts%>%group_by(year)%>%filter(yr_tot==max(yr_tot))

ggplot(gostotmax, aes(x=year, y=Julian))+geom_point()+geom_hline(yintercept=268, linetype=2)+ggtitle("Goshutes")+ylab("peak migration date(Julian)")+
  stat_smooth(method="gam")+theme_classic()

#functional/species-level patterns
gos=GosMts%>%select(-yr_tot,)%>%
  pivot_longer(cols=3:27, names_to = "Species", values_to="Count")

gos=gos%>%mutate(
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

gos$mammals=as.numeric(gos$mammals)
gos$birds=as.numeric(gos$birds)
str(gos)
#annual totals per species
gos_sp=gos%>%
  group_by(year, Species, diet, mammals, birds)%>%
  summarise(yr_tot=sum(Count))%>%
  mutate(spcode=as.integer(as.character(Species)))%>%
  arrange(year)
  
#daily totals per species
gos_sp_day=gos%>%
  group_by(Species,diet, year, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(year)

#peak dates of migration for each species
gos_sp_day_max=gos_sp_day%>%
  group_by(diet, year, Species)%>%
  filter(total==max(total))%>%
  slice_head(n=1)%>%
  arrange(year)

#plot annual totals trend per species
ggplot(gos_sp, aes(x=year, y=log(yr_tot), col=Species))+geom_point()+
  stat_smooth(method="lm")+theme_classic()+facet_wrap(~Species)

#plot peak dates trend per species
ggplot(gos_sp_day_max, aes(x=year, y=Julian, col=Species))+geom_point()+
  stat_smooth(method="lm")+theme_classic()+facet_wrap(~Species)

#annual totals per functional group

gos_grp=gos%>%
  group_by(year, diet)%>%
  summarise(yr_tot=sum(Count))
  
#get daily totals by functional group
gos_grp_day=gos%>%
  group_by(year,Julian, diet)%>%
  summarise(total=sum(Count))%>%
  arrange(year)

#get peak dates per functional group
gos_grp_day_max=gos_grp_day%>%
  group_by(diet, year)%>%
  filter(total==max(total))%>%
  slice_head(n=1)%>%
  arrange(year)

#plot peak dates trend 
ggplot(gos_grp_day_max, aes(x=year, y=Julian, col=diet))+geom_point()+
  stat_smooth(method="gam")+theme_classic()

#plot annual trends per functional group
ggplot(gos_grp, aes(x=year, y=yr_tot, col=diet))+geom_point()+
  stat_smooth(method="gam")+theme_classic()

#categorical
gos_sp_max=gos%>%mutate(Julian=as.POSIXlt(Date)$yday)%>%
  group_by(year, Species, diet, mammals, birds)%>%
  filter(Count==max(Count), !Species%in%c("UA", "UB", "UE", "UF", "UR"))%>%slice_head(n=1)%>%
  arrange(year)

ggplot(gos_sp_max, aes(x=year, y=Julian, col=Species))+geom_point()+facet_wrap(~diet)+
  stat_smooth(method="lm")+geom_hline(yintercept=268, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()

ggplot(gos_sp, aes(x=year, y=log(yr_tot), col=Species))+geom_point()+facet_wrap(~diet)+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")


#continuous
p1=ggplot(gos_sp, aes(x=birds, y=log(yr_tot)))+geom_point()+xlab("% bird")+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")+ggtitle("Goshutes")

g2=ggplot(gos_sp, aes(x=mammals, y=log(yr_tot)))+geom_point()+xlab("% mammal")+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")+ggtitle("Goshutes")

gos_sp_max=gos%>%mutate(Julian=as.POSIXlt(Date)$yday)%>%
  group_by(year, Species, diet, mammals, birds)%>%
  filter(Count==max(Count), !Species%in%c("UA", "UB", "UE", "UF", "UR"))%>%slice_head(n=1)%>%
  arrange(year)

g4=ggplot(gos_sp_max, aes(x=mammals, y=Julian))+geom_point()+
  stat_smooth(method="lm")+geom_hline(yintercept=268, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()+ggtitle("Goshutes")+xlab("% mammal")

g3=ggplot(gos_sp_max, aes(x=birds, y=Julian))+geom_point()+
  stat_smooth(method="lm")+geom_hline(yintercept=268, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()+ggtitle("Goshutes")+xlab("% bird")

ggarrange(g1,g2,g3,g4)

#mammals
gos_sp_mammd=gos_sp_max%>%filter(diet=="mammal")
gos_sp_mammc=gos_sp%>%filter(diet=="mammal")

g5=ggplot(gos_sp_mammd, aes(x=year, y=Julian, col=mammals))+geom_point()+
  stat_smooth(method="lm")+geom_hline(yintercept=268, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()+ggtitle("mammal diet")

g6=ggplot(gos_sp_mammc, aes(x=year, y=log(yr_tot), col=mammals))+geom_point()+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")+ggtitle("mammal diet")

ggarrange(g2,g4, g6, g5)

#birds
gos_sp_aved=gos_sp_max%>%filter(diet=="bird")
gos_sp_avec=gos_sp%>%filter(diet=="bird")

g7=ggplot(gos_sp_aved, aes(x=year, y=Julian, col=birds))+geom_point()+
  stat_smooth(method="lm")+geom_hline(yintercept=268, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()+ggtitle("bird diet")

g8=ggplot(gos_sp_avec, aes(x=year, y=log(yr_tot), col=birds))+geom_point()+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")+ggtitle("bird diet")

ggarrange(g1,g3, g8, g7)

summary(lm(gos_sp_max$Julian~gos_sp_max$mammals))

#diversity values model####
g1=gsw%>%group_by(year)%>%filter(sd==max(sd))

t1=left_join(gos_sp_day_max, g1, by="year")%>%select(Species, diet, year, Julian.x, total, sd, Julian.y)%>%
  rename("sp_peak"="Julian.x", "div_index"="sd", "div_peak"="Julian.y")

group_peak=gos_grp_day_max%>%select(-total)%>%pivot_wider(names_from = diet, values_from = Julian )
group_numbers=gos_grp_day_max%>%select(-Julian)%>%pivot_wider(names_from = diet, values_from = total )%>%
  rename(bird_ct=bird, mam_ct=mammal, other_ct=other)

div_dat=g1%>%select(year, sd, sr, Julian, yr_tot)
div_dat1=left_join(div_dat, group_peak)
div_dat2=left_join(div_dat1, group_numbers)

t1$years = (t1$year - mean(t1$year))/(2 *sd(t1$year))

