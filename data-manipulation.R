require(dplyr)
require(lubridate)
require(tidyr)

#per species
GosMts=read.csv("GosMts.csv")
gos=GosMts%>%select(-TOTAL)%>%
  pivot_longer(cols=3:27, names_to = "Species", values_to="Count")%>%
  mutate(year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!is.na(Count))

#Goshute Mts####
#per year: what I want
gostot=GosMts%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:27], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>% #get annual totals across all spp.
  select(Date,year, month, date, Obs, yr_tot) #get only necessary columns

gostot$Julian=as.POSIXlt(gostot$Date)$yday 

gos=gostot%>%group_by(year)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                     cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                     MPD=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
            group_by(year)%>%filter(MPD==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#Hawk Ridge###
HawRd=read.csv("HawkRid.csv")
hrtot=HawRd%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:24], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>% #get annual totals across all spp.
  select(Date,year, month, date,  yr_tot) #get only necessary columns

hrtot$Julian=as.POSIXlt(hrtot$Date)$yday 

hr=hrtot%>%group_by(year)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                     cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                     MPD=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(year)%>%filter(MPD==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#Hawk Mt####

hms=read.csv("HMSDaily.csv")
#combine baea and goea data
hms=hms%>%
  mutate(BAEA=rowSums(.[9:11]), GOEA=rowSums(.[21:23]))%>%
  select(-OBSERVERS, -GOEA.I, -GOEA.A, -GOEA.U, -BAEA.I, -BAEA.U, -BAEA.A, -RAPTORS)

hmstot=hms%>%
  mutate(yr_tot=rowSums(.[6:30], na.rm=T), date=make_date(YR, MO, DAY))%>%
  select(date, DAY, MO, YR, HRS, yr_tot)
  
hm=hmstot%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                   cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                   MPD=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(MPD==1, !(YR<1983))%>%slice_head(n=1)  #select rows where cum_pt is >50%

hm$Julian=as.POSIXlt(hm$date)$yday 

#plot####
par(mfrow=c(1,3))

plot(gos$Julian~gos$year, type="l", col="blue", ylab="mean passage date(Julian)", xlab="year", main="Goshute Mts. (west)")
median(gos$Julian)
abline(h=268, lty=2)

plot(hr$Julian~hr$year, type="l", col="blue", ylab="mean passage date(Julian)", xlab="year", main="Hawk Ridge (central)")
median(hr$Julian)
abline(h=261, lty=2)

plot(hm$Julian~hm$YR, type="l", col="blue", ylab="mean passage date(Julian)", xlab="year", main="Hawk Mt (east)")
median(hm$Julian)
abline(h=357, lty=2)

