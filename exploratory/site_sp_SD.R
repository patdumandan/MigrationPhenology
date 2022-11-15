#Goshutes sample####
require(dplyr)
require(lubridate)
require(tidyr)
require(vegan)
require(ggplot2)

#per species
#GOSHUTES####
GosMts=read.csv("GosMts.csv")

gos=GosMts%>%select(-TOTAL)%>%
  pivot_longer(cols=3:27, names_to = "Species", values_to="Count")%>%
  mutate(year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!is.na(Count))

GosMts[, 3:27][is.na(GosMts[, 3:27])] <- 0

#DAILY TOTALS
gostot=GosMts%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:27], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!year>2018)
gostot$Julian=as.POSIXlt(gostot$Date)$yday 

#GET MPD
#gosmpd=GosMts%>%select(-TOTAL)%>%
#  mutate(yr_tot=rowSums(.[3:27], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>% #get annual totals across all spp.
#  select(Date,year, month, date, Obs, yr_tot) #get only necessary columns

gos=gostot%>%group_by(year)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                     cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                     MPD=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(year)%>%filter(MPD==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#community metrics####
GosMts[, 3:27][is.na(GosMts[, 3:27])] <- 0
gsw=gostot%>%select( -Obs, -UA, -UR, -UB, -UE, -UF)%>%filter(!(year<1983))

gsw=gsw%>%mutate(sd=diversity(gsw[2:21], index="shannon"), sr=specnumber(gsw[2:21]))

#find date when SD nd SR were highest

g1=gsw%>%group_by(year)%>%filter(sd==max(sd))

#species-level####

gos_sp=GosMts%>%select(-TOTAL)%>%
  pivot_longer(cols=3:27, names_to = "Species", values_to="Count")%>%
  mutate(year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!is.na(Count))

spmpd=gos_sp%>%group_by(Species,Date)%>%summarise(sptot=sum(Count))%>%
  filter(!Species%in%c("UA", "UR", "UB", "UE", "UF"))
spmpd$Julian=as.POSIXlt(spmpd$Date)$yday 

indextab=gsw%>%select(Date,year, month, date, Julian, sd, sr)#all time sr and sd
indextab_maxsd=g1%>%select(Date,year, month, date, Julian, sd, sr)#all time sr and sd
indtab1=left_join(indextab_maxsd, spmpd)

#plot

ggplot(indtab1, aes(x=Julian, y=sd, col=year))+geom_line()+
  annotate("rect", xmin=250,  xmax=270, alpha=0.1, fill="blue",ymin=1.75, ymax=2.4)+
  ggtitle("Goshutes")+theme_classic()+geom_hline(yintercept=1.91, linetype=2)+
  ylab("H' index")+xlab("date of highest SD (Julian)")

ggplot(indtab1, aes(x=Julian, y=sptot, col=Species))+
  geom_line()+ggtitle("Goshutes")+ylab("Count")+xlab("date of highest SD (Julian)")+
  theme_classic()+annotate("rect", xmin=250, xmax=270, alpha=0.1, fill="blue", 
                           ymin=0, ymax=550)+geom_label(label=indtab1$Species)

#Hawk Mt####

hms=read.csv("HMSDay.csv")
#combine baea and goea data
hmstot=hms%>%
  mutate(BAEA=rowSums(.[9:11]), GOEA=rowSums(.[21:23]))%>%
  select(-OBSERVERS, -GOEA.I, -GOEA.A, -GOEA.U, -BAEA.I, -BAEA.U, -BAEA.A, -RAPTORS)%>%
  mutate(yr_tot=rowSums(.[6:30], na.rm=T), date=make_date(YR, MO, DAY))

hmstot$Julian=as.POSIXlt(hmstot$date)$yday 

hm=hmstot%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                  cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                  MPD=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(MPD==1, !(YR<1983))%>%slice_head(n=1)  #select rows where cum_pt is >50%

hm$Julian=as.POSIXlt(hm$date)$yday 
hm=hm%>%mutate(date=make_date(YR, MO, DAY))

#community metrics####
hsw=hms%>%select( -OBS,-HRS, -UNID.ACCIPITER,
                  -UNID.BUTEO, -UNID.EAGLE, -UNID.FALCON, -UNID.RAPTOR)%>%filter(!(YR<1983))

hsw[, 4:23][is.na(hsw[, 4:23])] <- 0

hsw=hsw%>%mutate(sd=diversity(hsw[4:23], index="shannon"), sr=specnumber(hsw[4:23]), 
                 Date=make_date(YR, MO, DAY), total=rowSums(.[4:23]))

hsw$Julian=as.POSIXlt(hsw$Date)$yday 

#find date when SD nd SR were highest

f1=hsw%>%group_by(YR)%>%filter(sd==max(sd))

#species-level####

hms_sp=hms%>%select(-UNID.ACCIPITER, -UNID.BUTEO,-UNID.EAGLE, -UNID.FALCON, -UNID.RAPTOR)%>%
  pivot_longer(cols=6:30, names_to = "Species", values_to="Count")%>%
  filter(!is.na(Count), !YR<1983, !Species%in%c("UNID.ACCIPITER",
                                                "UNID.BUTEO", 
                                                "UNID.EAGLE",
                                                "UNID.FALCON",
                                                "UNID.RAPTOR"))%>%
  mutate(Date=make_date(YR, MO, DAY))

hmsmpd=hms_sp%>%group_by(Species,Date)%>%summarise(sptot=sum(Count))
hmsmpd$Julian=as.POSIXlt(hmsmpd$Date)$yday 

indextabh=hsw%>%select(Date,YR, MO, DAY, Julian, sd, sr)#all time sr and sd
indextabh_maxsd=f1%>%select(Date,YR, MO, DAY, Julian, sd, sr)#all time sr and sd
indtab2=left_join(indextabh_maxsd, hmsmpd)

#plot

ggplot(indtab2, aes(x=Julian, y=sd, col=YR))+geom_line()+
  annotate("rect", xmin=230,  xmax=250, alpha=0.1, fill="blue",ymin=1.55, ymax=2.2)+
  annotate("rect", xmin=330,  xmax=350, alpha=0.1, fill="blue",ymin=1.55, ymax=2.2)+
  ggtitle("Hawk Mountain")+theme_classic()+geom_hline(yintercept=1.92, linetype=2)+
  ylab("H' index")+xlab("date of highest SD (Julian)")

ggplot(indtab2, aes(x=Julian, y=sptot, col=Species))+geom_label(label=indtab2$Species)+
  geom_line()+ggtitle("Hawk Mountain")+ylab("Count")+xlab("date of highest SD (Julian)")+
  theme_classic()+annotate("rect", xmin=230, xmax=250, alpha=0.1, fill="blue", 
                           ymin=0, ymax=80)+
  annotate("rect", xmin=330, xmax=350, alpha=0.1, fill="blue", 
           ymin=0, ymax=80)

#Hawk Ridge####
HawRd=read.csv("HawkRid.csv")
HawRd[, 3:24][is.na(HawRd[, 3:24])] <- 0

hrtot=HawRd%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:24], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!year>2015) #get annual totals across all spp.

hrtot$Julian=as.POSIXlt(hrtot$Date)$yday 

#community metrics####
rsw=hrtot%>%select( -ObsHrs, -UA, -UB, -UF, -UE, -UR)%>%filter(!(year<1983))

rsw=rsw%>%mutate(sd=diversity(rsw[2:18], index="shannon"), sr=specnumber(rsw[2:18]))

r1=rsw%>%group_by(year)%>%filter(sd==max(sd))

#species-level####

hr_sp=rsw%>%
  pivot_longer(cols=2:18, names_to = "Species", values_to="Count")%>%
  filter(!is.na(Count))

hrmpd=hr_sp%>%group_by(Species,Date)%>%summarise(sptot=sum(Count))
hrmpd$Julian=as.POSIXlt(hrmpd$Date)$yday 

indextabhr=rsw%>%select(Date,year, month, date, Julian, sd, sr)#all time sr and sd
indextabhr_maxsd=r1%>%select(Date,year, month, date, Julian, sd, sr)#all time sr and sd
indtab3=left_join(indextabhr_maxsd, hrmpd)

#plot

ggplot(indtab3, aes(x=Julian, y=sd, col=year))+geom_line()+
  annotate("rect", xmin=225,  xmax=240, alpha=0.1, fill="blue",ymin=1.55, ymax=2.1)+
  annotate("rect", xmin=255,  xmax=265, alpha=0.1, fill="blue",ymin=1.55, ymax=2.1)+
  ggtitle("Hawk Ridge")+theme_classic()+geom_hline(yintercept=1.78, linetype=2)+
  ylab("H' index")+xlab("date of highest SD (Julian)")

ggplot(indtab3, aes(x=Julian, y=sptot, col=Species))+geom_label(label=indtab3$Species)+
  geom_line()+ggtitle("Hawk Ridge")+ylab("Count")+xlab("date of highest SD (Julian)")+
  theme_classic()+annotate("rect", xmin=225, xmax=240, alpha=0.1, fill="blue", 
                           ymin=0, ymax=350)+
  annotate("rect", xmin=255, xmax=265, alpha=0.1, fill="blue", 
           ymin=0, ymax=350)

####DECOMPOSING DIVERSITY####

#GOSHUTES####
gos_sp_max=gos%>%mutate(Julian=as.POSIXlt(Date)$yday)%>%
  group_by(year, Species)%>%
  filter(Count==max(Count), !Species%in%c("UA", "UB", "UE", "UF", "UR"))%>%slice_head(n=1)%>%
  arrange(year)

gos_sp=gos%>%mutate(Julian=as.POSIXlt(Date)$yday)%>%
  group_by(year, Species)%>%
  summarise(yr_tot=sum(Count))%>%filter(!Species%in%c("UA", "UB", "UE", "UF", "UR"))%>%
  arrange(year)

ggplot(gos_sp_max, aes(x=year, y=Julian))+geom_point(col="blue")+
  facet_wrap(~Species)+geom_hline(yintercept=268, linetype=2)

ggplot(gos_sp, aes(x=year, y=log(yr_tot)))+geom_point()+facet_wrap(~Species)+
  stat_smooth(method="lm")

ggplot(gos_sp_max, aes(x=year, y=Julian))+geom_point()+facet_wrap(~Species)+
  stat_smooth(method="lm")+geom_hline(yintercept=268, linetype=2)


#DAILY TOTALS
gostot=GosMts%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:27], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!year>2018)
gostot$Julian=as.POSIXlt(gostot$Date)$yday 

gostotmax=gostot%>%group_by(year)%>%filter(yr_tot==max(yr_tot))

ggplot(gostotmax, aes(x=year, y=Julian))+geom_point()+geom_hline(yintercept=268, linetype=2)+ggtitle("Goshutes")+
  stat_smooth(method="lm")  

#hAWK mT####
#Hawk Mt####

hms=read.csv("HMSDay.csv")
#combine baea and goea data
hmstot=hms%>%
  mutate(BAEA=rowSums(.[9:11]), GOEA=rowSums(.[21:23]))%>%
  select(-OBSERVERS, -GOEA.I, -GOEA.A, -GOEA.U, -BAEA.I, -BAEA.U, -BAEA.A, -RAPTORS)%>%
  mutate(yr_tot=rowSums(.[6:30], na.rm=T), date=make_date(YR, MO, DAY))%>%filter(!YR<1983, !YR>2018)

hmstot$Julian=as.POSIXlt(hmstot$date)$yday 

#species-level####

hms_sp=hmstot%>%select(-yr_tot)%>%
  pivot_longer(cols=6:30, names_to = "Species", values_to="Count")%>%
  filter(!is.na(Count), !YR<1983, !Species%in%c("UNID.ACCIPITER",
                                                "UNID.BUTEO", 
                                                "UNID.EAGLE",
                                                "UNID.FALCON",
                                                "UNID.RAPTOR"))%>%
  mutate(Date=make_date(YR, MO, DAY))

#PEAK DAYS PER SPECIES
hms_sp_max=hms_sp%>%mutate(Julian=as.POSIXlt(Date)$yday)%>%
  group_by(YR, Species)%>%
  filter(Count==max(Count))%>%slice_head(n=1)%>%
  arrange(YR)

#PEAK DAYS FOR SITE
hmsmax=hmstot%>%
  select(-UNID.ACCIPITER, -UNID.BUTEO, -UNID.EAGLE, -UNID.FALCON, -UNID.RAPTOR)%>%
  mutate(Julian=as.POSIXlt(date)$yday)%>%
  group_by(YR)%>%
  filter(yr_tot==max(yr_tot))%>%
  arrange(YR)

ggplot(hmsmax, aes(x=YR, y=Julian))+geom_point()+stat_smooth(method="lm")+theme_classic()+
  geom_hline(yintercept=260, linetype=2)+ylab("peak migration date(Julian)")

ggplot(hms_sp_max, aes(x=YR, y=Julian))+geom_point()+stat_smooth(method="lm")+theme_classic()+
  facet_wrap(~Species)+geom_hline(yintercept=260, linetype=2)+ylab("peak migration date(Julian)")

ggplot(hms_sp_max, aes(x=YR, y=Count))+geom_point()+stat_smooth(method="lm")+theme_classic()+
  facet_wrap(~Species)+geom_hline(yintercept=260, linetype=2)+ylab("max counts (log)")

#Hawk Ridge####
HawRd=read.csv("HawkRid.csv")
HawRd[, 3:24][is.na(HawRd[, 3:24])] <- 0

#without unknowns
hrtot=HawRd%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:19], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!year>2018) #get annual totals across all spp.

hrtot$Julian=as.POSIXlt(hrtot$Date)$yday 

#species-level####

hr_sp=rsw%>%
  pivot_longer(cols=2:18, names_to = "Species", values_to="Count")%>%
  filter(!is.na(Count))

hr_sp_max=hr_sp%>%group_by(Species,year)%>%filter(Count==max(Count))%>%slice_head(n=1)%>%select(-yr_tot)
hrmax=hrtot%>%group_by(year)%>%filter(yr_tot==max(yr_tot))%>%slice_head(n=1)%>%select(-UA, -UB,-UE,-UF,-UR)%>%
  mutate(Julian=as.POSIXlt(Date)$yday )

#plot

ggplot(hrmax, aes(x=year, y=Julian))+geom_point()+stat_smooth(method="lm")+theme_classic()+
  geom_hline(yintercept=257, linetype=2)+ylab("peak migration date(Julian)")

ggplot(hr_sp_max, aes(x=year, y=Julian))+geom_point()+stat_smooth(method="lm")+theme_classic()+
  facet_wrap(~Species)+geom_hline(yintercept=257, linetype=2)+ylab("peak migration date(Julian)")

ggplot(hr_sp_max, aes(x=year, y=Count))+geom_point()+stat_smooth(method="lm")+theme_classic()+
  facet_wrap(~Species)+geom_hline(yintercept=257, linetype=2)+ylab("max counts (log)")

