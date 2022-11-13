require(dplyr)
require(lubridate)
require(tidyr)
require(vegan)
require(ggpubr)
require(ggplot2)

#per species
#GOSHUTES####
GosMts=read.csv("GosMts.csv")

#gos=GosMts%>%select(-TOTAL)%>%
#  pivot_longer(cols=3:27, names_to = "Species", values_to="Count")%>%
#  mutate(year=year(Date), month=month(Date), date=day(Date))%>%
#  filter(!is.na(Count))

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
gsw=gostot%>%select( -Obs)%>%filter(!(year<1983))

gsw=gsw%>%mutate(sd=diversity(gsw[2:26], index="shannon"), sr=specnumber(gsw[2:26]))

#find date when SD nd SR were highest

g1=gsw%>%group_by(year)%>%filter(sd==max(sd))
g2=gsw%>%group_by(year)%>%filter(sr==max(sr))
g3=gsw%>%group_by(year)%>%filter(yr_tot==max(yr_tot))

#Hawk Ridge###
HawRd=read.csv("HawkRid.csv")
HawRd[, 3:24][is.na(HawRd[, 3:24])] <- 0

hrtot=HawRd%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:24], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!year>2018) #get annual totals across all spp.

hrtot$Julian=as.POSIXlt(hrtot$Date)$yday 

hr=hrtot%>%select(Date,year, month, date,  yr_tot)%>%group_by(year)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                     cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                     MPD=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(year)%>%filter(MPD==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#community metrics####
HawRd[, 3:24][is.na(HawRd[, 3:24])] <- 0
hsw=hrtot%>%select( -ObsHrs)%>%filter(!(year<1983))

hsw=hsw%>%mutate(sd=diversity(hsw[3:24], index="shannon"), sr=specnumber(hsw[3:24]))

#Hawk Mt####

hms=read.csv("HMSDay.csv")
#combine baea and goea data
hms=hms%>%
  mutate(BAEA=rowSums(.[9:11]), GOEA=rowSums(.[21:23]))%>%
  select(-OBSERVERS, -GOEA.I, -GOEA.A, -GOEA.U, -BAEA.I, -BAEA.U, -BAEA.A, -RAPTORS)

hmstot=hms%>%
  mutate(yr_tot=rowSums(.[6:30], na.rm=T), date=make_date(YR, MO, DAY))%>%
  select(date, DAY, MO, YR, HRS, yr_tot)

hmstot$Julian=as.POSIXlt(hmstot$date)$yday 

hm=hmstot%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                   cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                   MPD=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(MPD==1, !(YR<1983))%>%slice_head(n=1)  #select rows where cum_pt is >50%

hm$Julian=as.POSIXlt(hm$date)$yday 
hm=hm%>%mutate(date=make_date(YR, MO, DAY))

#plot####
par(mfrow=c(1,3))

plot(gos$Julian~gos$year, type="l", col="blue", ylab="mean passage date(Julian)", xlab="year", main="Goshute Mts. (west)")
median(gos$Julian)
abline(h=268, lty=2)

plot(hr$Julian~hr$year, type="l", col="blue", ylab="mean passage date(Julian)", xlab="year", main="Hawk Ridge (central)")
median(hr$Julian)
abline(h=261, lty=2)

#community metrics####
dsw=hms%>%select( -OBS,-HRS)%>%filter(!(YR<1983))

dsw[, 4:28][is.na(dsw[, 4:28])] <- 0

dsw=dsw%>%mutate(sd=diversity(dsw[1:25], index="shannon"), sr=specnumber(dsw[1:25]), 
                 date=make_date(YR, MO, DAY), total=rowSums(.[4:28]))

dsw$Julian=as.POSIXlt(dsw$date)$yday 

#find date when SD nd SR were highest

f1=dsw%>%group_by(YR)%>%filter(sd==max(sd))
f2=dsw%>%group_by(YR)%>%filter(sr==max(sr))
f3=dsw%>%group_by(YR)%>%filter(total==max(total))

h1=hsw%>%group_by(year)%>%filter(sd==max(sd))
h2=hsw%>%group_by(year)%>%filter(sr==max(sr))
h3=hsw%>%group_by(year)%>%filter(yr_tot==max(yr_tot))

#PLOTS####
par(mfrow=c(3,1))
#target dates
#SR
plot(f2$Julian~f2$YR, type="l", col="blue", main="Hawk Mt(east)", ylab="Julian date", xlab="year")
median(f2$Julian)
abline(h=283, lty=2, col="blue")

#SD
lines(f1$Julian~f1$YR, type="l", col="red")
median(f1$Julian)
abline(h=265, lty=2, col="red")

#daily highest
lines(f3$Julian~f3$YR, type="l", col="black")
median(f3$Julian)
abline(h=260, lty=2, col="black")

#target dates
#SR

plot(h2$Julian~h2$year, type="l", col="blue", main="Hawk Ridge(central)", ylab="Julian date", xlab="year", ylim=c(230,305))
median(h2$Julian)
abline(h=269, lty=2, col="blue")

#SD
lines(h1$Julian~h1$year, type="l", col="red")
median(h1$Julian)
abline(h=240, lty=2, col="red")

#daily highest
lines(h3$Julian~h3$year, type="l", col="black")
median(h3$Julian)
abline(h=257, lty=2, col="black")

#Goshutes
plot(g2$Julian~g2$year, type="l", col="blue", main="Goshute Mts(west)", ylab="Julian date", xlab="year", ylim=c(230,320))
median(g2$Julian)
abline(h=264, lty=2, col="blue")

#SD
lines(g1$Julian~g1$year, type="l", col="red")
median(g1$Julian)
abline(h=254, lty=2, col="red")

#daily highest
lines(g3$Julian~g3$year, type="l", col="black")
median(g3$Julian)
abline(h=268, lty=2, col="black")

legend( x= "topright",inset=0,legend=c("Species Richness", "Species Diversity", "Highest Daily Count"),
       col=c("red", "blue", "black"), lty=1, cex=0.8, xpd=F)



