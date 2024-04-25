
#data manipulation####
qr=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/QuakerRidge.csv")

qr$Date=lubridate::mdy(qr$Date)

qr=qr%>%
  mutate(year=year(Date), month=month(Date), day=day(Date),
         AMKE=rowSums(.[4:7], na.rm=TRUE),
         BAEA=rowSums(.[8:13], na.rm=TRUE),
         BLVU=rowSums(.[14:15], na.rm=TRUE),
         BWHA=rowSums(.[16:19], na.rm=TRUE),
         COHA=rowSums(.[20:23], na.rm=TRUE),
         GOEA=rowSums(.[24:27], na.rm=TRUE),
         MERL=rowSums(.[28:31], na.rm=TRUE),
         NOGO=rowSums(.[33:35], na.rm=TRUE),
         NOHA=rowSums(.[36:40], na.rm=TRUE),
         OSPR=rowSums(.[41:42], na.rm=TRUE),
         PEFA=rowSums(.[43:46], na.rm=TRUE),
         RSHA=rowSums(.[48:51], na.rm=TRUE),
         RTHA=rowSums(.[52:55], na.rm=TRUE),
         SSHA=rowSums(.[57:60], na.rm=TRUE),
         TUVU=rowSums(.[62:64], na.rm=TRUE))%>%
  select(Date,year, month, day, Duration, Observer, 
         AMKE, BAEA, BLVU, BWHA, COHA, GOEA, MERL, MK, NOGO, NOHA, OSPR, PEFA,
         RL, RSHA, RTHA, SSHA, SW, SE, TUVU)%>%
  mutate(TOTAL=rowSums(.[7:24], na.rm=T))

qr$Julian=as.POSIXlt(qr$Date)$yday

QR=qr%>%filter(!year<1985, !Julian<231, !Julian >323)


#species-level counts####
qr_sp=QR%>%select(-TOTAL)%>%
  pivot_longer(cols=7:25, names_to = "Species", values_to="Count")

#annual totals per species
qr_sp_tot=qr_sp%>%
  group_by(year, Species)%>%
  summarise(sp_tot=sum(Count))

#determine uncommon species (based on non-detects annually)
qr_sp_min=qr_sp_tot%>%group_by(Species)%>%filter(sp_tot==min(sp_tot))

#daily totals per species
qr_sp_day=qr_sp%>%
  group_by(Species,year, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(year)

#10% passage date: site
qr_10pd_site=QR%>%group_by(year)%>%mutate(cum_tot=cumsum(TOTAL), #get cumulative sum for each year
                                          cum_pt=cum_tot/sum(TOTAL),#calculate % of daily count/annual count
                                          PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >10%
  group_by(year)%>%filter(PD_10==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(year, Julian, cum_pt)

#50% passage date

qr_50pd_site=QR%>%group_by(year)%>%mutate(cum_tot=cumsum(TOTAL), #get cumulative sum for each year
                                               cum_pt=cum_tot/sum(TOTAL),#calculate % of daily count/annual count
                                               PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >10%
  group_by(year)%>%filter(PD_50==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(year, Julian, cum_pt)

#90% passage date

qr_90pd_site=QR%>%group_by(year)%>%mutate(cum_tot=cumsum(TOTAL), #get cumulative sum for each year
                                               cum_pt=cum_tot/sum(TOTAL),#calculate % of daily count/annual count
                                               PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >10%
  group_by(year)%>%filter(PD_90==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(year, Julian, cum_pt)

ggplot(qr_10pd_site, aes(x=year, y=Julian))+geom_line()+theme_classic()+ylab("Julian date")+xlab("year")+ggtitle("10% passage date")+geom_hline(yintercept=256,linetype=2)

ggplot(qr_50pd_site, aes(x=year, y=Julian))+geom_line()+theme_classic()+ylab("Julian date")+xlab("year")+ggtitle("50% passage date")+geom_hline(yintercept=260,linetype=2)

ggplot(qr_90pd_site, aes(x=year, y=Julian))+geom_line()+theme_classic()+ylab("Julian date")+xlab("year")+ggtitle("90% passage date")+geom_hline(yintercept=276,linetype=2)
