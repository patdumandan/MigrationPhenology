require(dplyr)
require(tidyr)
require(lubridate)
require(ggplot2)
require(vegan)
GosMts=read.csv("GosMts.csv")

#site-level patterns
gostot=GosMts%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:27], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!year>2018)
gostot$Julian=as.POSIXlt(gostot$Date)$yday 

#trend in max migration dates
gostotmax=gostot%>%group_by(year)%>%filter(yr_tot==max(yr_tot))

ggplot(gostotmax, aes(x=year, y=Julian))+geom_point()+geom_hline(yintercept=268, linetype=2)+ggtitle("Goshutes")+ylab("peak migration date(Julian)")+
  stat_smooth(method="lm")+theme_classic()

#community metrics####
GosMts[, 3:27][is.na(GosMts[, 3:27])] <- 0
gsw=gostot%>%select( -Obs)%>%filter(!(year<1983))

require(vegan)
#calculate diversity and richness
gsw=gsw%>%mutate(sd=diversity(gsw[2:26], index="shannon"), sr=specnumber(gsw[2:26]))
gsw[, 2:26][is.na(gsw[, 2:26])] <- 0

#find date when SD nd SR were highest

g1=gsw%>%group_by(year)%>%filter(sd==max(sd))
g2=gsw%>%group_by(year)%>%filter(sr==max(sr))
g3=gsw%>%group_by(year)%>%filter(yr_tot==max(yr_tot))

#Goshutes
#par(mfrow=c(2,1))

plot(g1$Julian~g1$year, type="l", col="black", main="species diversity", ylab="peak diversity date", xlab="year", ylim=c(230,308))
abline(h=254, lty=2, col="black")
rect(xleft=1995, ybottom=0, ytop=308, xright=2005, col = rgb(0,0,0.5, 1/4)) #control

plot(g1$sd~g1$year, type="l")
abline(h=2.04, lty=2)
rect(xleft=1995, ybottom=0, ytop=308, xright=2005, col = rgb(0,0,0.5, 1/4)) #control

#SD
plot(g1$Julian~g1$year, type="l", col="red")
median(g1$Julian)
abline(h=254, lty=2, col="red")

#daily highest
lines(g3$Julian~g3$year, type="l", col="black")
median(g3$Julian)
abline(h=268, lty=2, col="black")

legend( x= "topright",inset=0,legend=c("Species Richness", "Species Diversity", "Highest Daily Count"),
        col=c("red", "blue", "black"), lty=1, cex=0.8, xpd=F)



#functional/species-level patterns
gos=GosMts%>%select(-TOTAL)%>%
  pivot_longer(cols=3:27, names_to = "Species", values_to="Count")%>%
  mutate(year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!is.na(Count))

gos=gos%>%mutate(Julian=as.POSIXlt(Date)$yday,
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

gos_sp=gos%>%mutate(Julian=as.POSIXlt(Date)$yday)%>%
  group_by(year, Species, diet, mammals, birds)%>%
  summarise(yr_tot=sum(Count))%>%filter(!Species%in%c("UA", "UB", "UE", "UF", "UR"))%>%
  arrange(year)

gos_sp_day=gos%>%mutate(Julian=as.POSIXlt(Date)$yday)%>%
  filter(!Julian<226, !Julian>308)%>%
  group_by(Species,diet, year, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(year)

gos_sp_day_max=gos_sp_day%>%
  group_by(diet, year, Species)%>%
  filter(total==max(total), !year>2018)%>%
  slice_head(n=1)%>%
  arrange(year)

ggplot(gos_sp_day_max, aes(x=year, y=log(total), col=Species))+geom_point()+
  stat_smooth(method="lm")+theme_classic()+facet_wrap(~diet)+
  geom_vline(xintercept=2005, lty=2)+
  geom_vline(xintercept=1996, lty=2)

ggplot(gos_sp_day_max, aes(x=year, y=Julian, col=Species))+geom_point()+
  stat_smooth(method="lm")+theme_classic()+facet_wrap(~diet)+
  geom_vline(xintercept=2005, lty=2)+
  geom_vline(xintercept=1996, lty=2)

gos_grp_day=gos%>%mutate(Julian=as.POSIXlt(Date)$yday)%>%
  filter(!Julian<226, !Julian>308)%>%
  group_by(year,Julian, diet)%>%
  summarise(total=sum(Count))%>%
  arrange(year)

gos_grp_day_max=gos_grp_day%>%
  group_by(diet, year)%>%
  filter(total==max(total))%>%
  slice_head(n=1)%>%
  arrange(year)

ggplot(gos_grp_day_max, aes(x=year, y=Julian, col=diet))+geom_point()+
  stat_smooth(method="gam")+theme_classic()+
  geom_vline(xintercept=2005, lty=2)+
  geom_vline(xintercept=1996, lty=2)

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

