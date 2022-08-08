require(dplyr)
require(lubridate)
require(tidyr)

#per species
gos=GosMts%>%select(-TOTAL)%>%
  pivot_longer(cols=3:27, names_to = "Species", values_to="Count")%>%
  mutate(year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!is.na(Count))

#Goshute Mts
#per year: what I want
gostot=GosMts%>%select(-TOTAL)%>%
  mutate(yr_tot=rowSums(.[3:27], na.rm=T), year=year(Date), month=month(Date), date=day(Date))%>% #get annual totals across all spp.
  select(Date,year, month, date, Obs, yr_tot) #get only necessary columns

gostot$Julian=as.POSIXlt(gostot$Date)$yday 

gos=gostot%>%group_by(year)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                     cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                     MPD=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
            group_by(year)%>%filter(MPD==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

plot(gos$Julian~gos$year, type="l", col="blue", ylab="mean passage date(Julian)", xlab="year", main="Goshute Mts.")
abline(h=268, lty=2)

