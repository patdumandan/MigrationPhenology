require(dplyr)
require(lubridate)
require(tidyr)
require(vegan)
require(ggplot2)
GosMts=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/GosMts.csv")
GosMts[, 3:27][is.na(GosMts[, 3:27])] <- 0
GosMts$Julian=as.POSIXlt(GosMts$Date)$yday
GosMts=GosMts%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:27], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!year>2018, !Julian>308)
hms=read.csv("HMSDay.csv")
hms=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/HMSDay.csv")
#combine baea and goea data
HMS=hms%>%
  mutate(BAEA=rowSums(.[9:11]), GOEA=rowSums(.[21:23]))%>%
  select(-OBSERVERS, -GOEA.I, -GOEA.A, -GOEA.U, -BAEA.I, -BAEA.U, -BAEA.A, -RAPTORS)%>%
  mutate(yr_tot=rowSums(.[6:30], na.rm=T), date=make_date(YR, MO, DAY))

HMS$Julian=as.POSIXlt(HMS$date)$yday
#community metrics###
hsw=HMS%>%select( -OBS,-HRS, -UNID.ACCIPITER,
                  -UNID.BUTEO, -UNID.EAGLE, -UNID.FALCON, -UNID.RAPTOR)%>%filter(!(YR<1983))

hsw[, 4:23][is.na(hsw[, 4:23])] <- 0
hsw=hsw%>%mutate(sd=diversity(hsw[4:23], index="shannon"), sr=specnumber(hsw[4:23]),
                 Date=make_date(YR, MO, DAY), total=rowSums(.[4:23]))

GosMts[, 3:27][is.na(GosMts[, 3:27])] <- 0
GosMts$Julian=as.POSIXlt(GosMts$Date)$yday
GosMts=GosMts%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:27], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!year>2018, !Julian>308)
GosMts=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/GosMts.csv")
GosMts[, 3:27][is.na(GosMts[, 3:27])] <- 0
GosMts$Julian=as.POSIXlt(GosMts$Date)$yday
GosMts=GosMts%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:27], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!year>2018, !Julian>308)
gsw=GosMts%>%select(-UA, -UR, -UB, -UE, -UF)%>%
  mutate(specdiv=diversity(GosMts[3:22], index="shannon"), specric=specnumber(GosMts[3:22]), community="all")

hsw=hsw%>%mutate(sd=diversity(hsw[4:23], index="shannon"), sr=specnumber(hsw[4:23]),
                 Date=make_date(YR, MO, DAY), total=rowSums(.[4:23]), assemblage="Hawk Mountain")
gsw=GosMts%>%select(-UA, -UR, -UB, -UE, -UF)%>%
  mutate(specdiv=diversity(GosMts[3:22], index="shannon"),
         specric=specnumber(GosMts[3:22]), assemblage="Goshutes")
HawRd=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/HawkRid.csv")
HawRd=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/HawkRd.csv")
HawRd=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/HawRd.csv")
HawRd[, 3:24][is.na(HawRd[, 3:24])] <- 0

HawRd$Julian=as.POSIXlt(HawRd$Date)$yday
HawRd$Date=as.Date(HawRd$Date)
str(HawRd)
HawRd$Date=as.Date.character(HawRd$Date)
HawRd$Date=as.Date(HawRd$Date, "%m/%d/%Y")
HawRd$Julian=as.POSIXlt(HawRd$Date)$yday
#community metrics###
rsw=HawRd%>%select( -ObsHrs, -UA, -UB, -UF, -UE, -UR)%>%filter(!(year<1983))
#community metrics###
rsw=HawRd%>%select( -ObsHrs, -UA, -UB, -UF, -UE, -UR)%>%filter(!(year<1983))%>%
  mutate(year=year(Date), month=month(Date), date=day(Date))
#community metrics###
rsw=HawRd%>%select( -ObsHrs, -UA, -UB, -UF, -UE, -UR)%>%
  mutate(year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!(year<1983))

rsw=rsw%>%mutate(sd=diversity(rsw[2:18], index="shannon"), sr=specnumber(rsw[2:18]), assemblage="Hawk Ridge")

h1=gsw%>%select(Date, Julian, year, specdiv, community)
h1=gsw%>%select(Date, Julian, year, specdiv, assemblage)
h2=hsw%>%select(Date, Julian, year, specdiv, assemblage)
h3=rsw%>%select(Date, Julian, year, specdiv, assemblage)
hsw=hsw%>%mutate(specdiv=diversity(hsw[4:23], index="shannon"), sr=specnumber(hsw[4:23]),
                 Date=make_date(YR, MO, DAY), total=rowSums(.[4:23]), assemblage="Hawk Mountain")
#community metrics###
rsw=HawRd%>%select( -ObsHrs, -UA, -UB, -UF, -UE, -UR)%>%
  mutate(year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!(year<1983))
rsw=rsw%>%mutate(specdiv=diversity(rsw[2:18], index="shannon"), sr=specnumber(rsw[2:18]), assemblage="Hawk Ridge")
h1=gsw%>%select(Date, Julian, year, specdiv, assemblage)
h2=hsw%>%select(Date, Julian, year, specdiv, assemblage)
#community metrics###
hsw=HMS%>%select( -OBS,-HRS, -UNID.ACCIPITER,
                  -UNID.BUTEO, -UNID.EAGLE, -UNID.FALCON, -UNID.RAPTOR)%>%filter(!(YR<1983))%>%
  rename("year"=="YR")
#community metrics###
hsw=HMS%>%select( -OBS,-HRS, -UNID.ACCIPITER,
                  -UNID.BUTEO, -UNID.EAGLE, -UNID.FALCON, -UNID.RAPTOR)%>%filter(!(YR<1983))%>%
  rename("year"="YR")
hsw[, 4:23][is.na(hsw[, 4:23])] <- 0
hsw=hsw%>%mutate(specdiv=diversity(hsw[4:23], index="shannon"), sr=specnumber(hsw[4:23]),
                 Date=make_date(YR, MO, DAY), total=rowSums(.[4:23]), assemblage="Hawk Mountain")
h2=hsw%>%select(Date, Julian, year, specdiv, assemblage)

#community metrics###
hsw=HMS%>%select( -OBS,-HRS, -UNID.ACCIPITER,
                  -UNID.BUTEO, -UNID.EAGLE, -UNID.FALCON, -UNID.RAPTOR)%>%filter(!(YR<1983))%>%
  rename("year"="YR")
hsw=hsw%>%mutate(specdiv=diversity(hsw[4:23], index="shannon"), sr=specnumber(hsw[4:23]),
                 Date=make_date(YR, MO, DAY), total=rowSums(.[4:23]), assemblage="Hawk Mountain")%>%
  rename("date"="Date")
hsw=hsw%>%mutate(specdiv=diversity(hsw[4:23], index="shannon"), sr=specnumber(hsw[4:23]),
                 Date=make_date(year, MO, DAY), total=rowSums(.[4:23]), assemblage="Hawk Mountain")%>%
  rename("date"="Date")
hsw=hsw%>%mutate(specdiv=diversity(hsw[4:23], index="shannon"), sr=specnumber(hsw[4:23]),
                 Date=make_date(year, MO, DAY), total=rowSums(.[4:23]), assemblage="Hawk Mountain")
hsw=HMS%>%select( -OBS,-HRS, -UNID.ACCIPITER,
                  -UNID.BUTEO, -UNID.EAGLE, -UNID.FALCON, -UNID.RAPTOR)%>%filter(!(YR<1983))%>%
  rename("year"="YR")
hsw[, 4:23][is.na(hsw[, 4:23])] <- 0
hsw=hsw%>%mutate(specdiv=diversity(hsw[4:23], index="shannon"), sr=specnumber(hsw[4:23]),
                 Date=make_date(year, MO, DAY), total=rowSums(.[4:23]), assemblage="Hawk Mountain")
#H
h2=hsw%>%select(Date, Julian, year, specdiv, assemblage)
h3=rsw%>%select(Date, Julian, year, specdiv, assemblage)
sample_sites=rbind(h1,h2,h3)

require(ggplot2)
p3=ggplot(full_gos11, aes(x=Julian, y=specdiv, col=assemblage))+geom_point()+ylab("H' index")+
  ggtitle("sample assemblages")+theme_classic()+stat_smooth(method="gam")
ggplot(sample_sites, aes(x=Julian, y=specdiv, col=assemblage))+geom_point()+ylab("H' index")+
  ggtitle("sample assemblages")+theme_classic()+stat_smooth(method="gam")
ggplot(sample_sites, aes(x=Julian, y=specdiv, col=assemblage))+geom_point()+ylab("H' index")+
  ggtitle("assemblage-level phenology")+theme_classic()+stat_smooth(method="gam")

#sample
