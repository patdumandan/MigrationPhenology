require(dplyr)
require(tidyr)
require(lubridate)
require(ggplot2)
require(broom)
require(brms)

#data manipulation####
hms=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/HMSDay.csv")

#combine baea and goea data
#select 11 years after first count

HMS=hms%>%
  mutate(BAEA=rowSums(.[9:11]), GOEA=rowSums(.[21:23]))%>%
  select(-OBSERVERS, -GOEA.I, -GOEA.A, -GOEA.U, -BAEA.I, -BAEA.U, -BAEA.A, -RAPTORS, -GYRF, -MIKI, -STKI,-SWHA,
         -UNID.ACCIPITER,-UNID.BUTEO, -UNID.EAGLE, -UNID.FALCON, -UNID.RAPTOR)%>%
  mutate(yr_tot=rowSums(.[6:21], na.rm=T), date=make_date(YR, MO, DAY))%>%filter(!(YR<1945))

HMS$Julian=as.POSIXlt(HMS$date)$yday
HMS[, 6:21][is.na(HMS[, 6:21])] <- 0

#select only 1946-2018 (Aug 15-Dec 15: 226-348)

HMS=HMS%>%filter(!Julian<226, !Julian>348, !YR<1979)

#species-level counts
hms_sp=HMS%>%
  pivot_longer(cols=6:21, names_to = "Species", values_to="Count")

#annual totals per species
hms_sp_tot=hms_sp%>%
  group_by(YR, Species)%>%
  summarise(sp_tot=sum(Count))

#species relative abundances

hms_sp_ra=hms_sp_tot%>%
  group_by(Species)%>%
  summarise(ave_ct=mean(sp_tot))%>%
  mutate(rel_abund=ave_ct/sum(ave_ct))

#determine minimum year of non-0 detects for all species
hms_sp_min=hms_sp_tot%>%group_by(Species)%>%filter(sp_tot>0)%>%summarise(min_yr=min(YR))

#daily totals per species
hms_sp_day=hms_sp%>%
  group_by(Species,YR, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(YR)%>%
  filter(!YR<1983)

#site-level####

##10%PD####

#10% passage date: species
hms_10pd_sp=hms_sp_day%>%filter(!YR<1983)%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                        cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                        PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_10==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#10% passage date: site
hms_10pd_site=HMS%>%filter(!YR<1983)%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                          cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                          PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_10==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

#site-level
priors = c(prior(normal(0, 10), class = Intercept),
            prior(normal(0,10), class = "sds", coef="s(YR)"))

mod10=brm(bf(Julian~ 1 + s(YR)), 
       data=hms_10pd_site, iter=2500, family=poisson(link="log"), 
       prior = priors,
       control=list(adapt_delta=0.99))

mod10pred=data.frame(predict(mod10))
colnames(mod10pred)[1]="pred_JD"

m1pred=cbind(hms_10pd_site, mod10pred)

mean(mod10pred[32:36,]$pred_JD)-mean(mod10pred[1:5,]$pred_JD) #0.73
median(mod10pred[32:36,]$pred_JD)-median(mod10pred[1:5,]$pred_JD) #-0.80

##50% PD####
#50% passage date: species
hms_50pd_sp=hms_sp_day%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                        cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                        PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_50==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#50% passage date: site
hms_50pd_site=HMS%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                          cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                          PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_50==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

#visualize data

plot(hms_50pd_site$Julian~hms_50pd_site$YR, type="l", main="Hawk Mountain", ylab="50% passage date", xlab="year")

segments(x0=1979, x1=1983,y0=median(hms_50pd_site[1:5,]$Julian), y1=median(hms_50pd_site[1:5,]$Julian), lwd=2, col="red")
segments(x0=2014, x1=2018,y0=median(hms_50pd_site[35:39,]$Julian), y1=median(hms_50pd_site[35:39,]$Julian), lwd=2, col="red")

segments(x0=1979, x1=1983,y0=mean(hms_50pd_site[1:5,]$Julian), y1=mean(hms_50pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2014, x1=2018,y0=mean(hms_50pd_site[35:39,]$Julian), y1=mean(hms_50pd_site[35:39,]$Julian), lwd=2, col="orange")


abline(v=1979, lty=2, col="grey", lwd=2)
abline(v=1983, lty=2, col="grey", lwd=2)
abline(v=2014, lty=2, col="grey", lwd=2)
abline(v=2018, lty=2, col="grey", lwd=2)


mod50=brm(bf(Julian~ 1+s(YR)), 
          data=hms_50pd_site, iter=2500, family=poisson(), control=list(adapt_delta=0.99))

mod50pred=data.frame(predict(mod50))
colnames(mod50pred)[1]="pred_JD"

m5pred=cbind(hms_50pd_site, mod50pred)

mean(mod50pred[36:40,]$pred_JD)-mean(mod50pred[1:5,]$pred_JD) #-7.9
median(mod50pred[36:40,]$pred_JD)-median(mod50pred[1:5,]$pred_JD) #-7.9

##90% PD####
#90% passage date: species
hms_90pd_sp=hms_sp_day%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                        cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                        PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_90==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

##90% passage date: site
hms_90pd_site=HMS%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                          cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                          PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_90==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

#visualize data

plot(hms_90pd_site$Julian~hms_90pd_site$YR, type="l", main="Hawk Mountain", ylab="90% passage date", xlab="year")

segments(x0=1979, x1=1983,y0=median(hms_90pd_site[1:5,]$Julian), y1=median(hms_90pd_site[1:5,]$Julian), lwd=2, col="red")
segments(x0=2014, x1=2018,y0=median(hms_90pd_site[35:39,]$Julian), y1=median(hms_90pd_site[35:39,]$Julian), lwd=2, col="red")

segments(x0=1979, x1=1983,y0=mean(hms_90pd_site[1:5,]$Julian), y1=mean(hms_90pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2014, x1=2018,y0=mean(hms_90pd_site[35:39,]$Julian), y1=mean(hms_90pd_site[35:39,]$Julian), lwd=2, col="orange")


abline(v=1979, lty=2, col="grey", lwd=2)
abline(v=1983, lty=2, col="grey", lwd=2)
abline(v=2014, lty=2, col="grey", lwd=2)
abline(v=2018, lty=2, col="grey", lwd=2)

#site-level
mod90=brm(bf(Julian~ 1+s(YR)), 
          data=hms_90pd_site, iter=2500, family=poisson(), control=list(adapt_delta=0.99))

mod90pred=data.frame(predict(mod90))
colnames(mod90pred)[1]="pred_JD"

m9pred=cbind(hms_90pd_site, mod90pred)

mean(mod90pred[36:40,]$pred_JD)-mean(mod90pred[1:5,]$pred_JD) #-7.7
median(mod90pred[36:40,]$pred_JD)-median(mod90pred[1:5,]$pred_JD) 


#species-level####

##10% PD####
sp_df_mod10=hms_10pd_sp%>%
  select(Species, YR, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(YR)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

sp_df_res10=sp_df_mod10%>%select(-data, -model)
sp_df_res10preds=as.data.frame(sp_df_res10$mean_pred[,1])
colnames(sp_df_res10preds)[1]="mean_preds"

main_dat10=hms_10pd_sp%>%select(Species, YR, Julian)

sp_df_all=cbind(sp_df_res10preds, main_dat10)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1983:1987) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))

colnames(sp_df_all)[1]="Species"

sp_df_diff=sp_df_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=median(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
sp_diff=left_join(sp_df_diff, hms_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)


#visualize data####

plot(hms_10pd_site$Julian~hms_10pd_site$YR, type="l", main="Hawk Mountain", ylab="10% passage date", xlab="year")

lines(m1pred$pred_JD~m1pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1983, x1=1987,y0=median(hms_10pd_site[1:5,]$Julian), y1=median(hms_10pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2014, x1=2018,y0=median(hms_10pd_site[32:36,]$Julian), y1=median(hms_10pd_site[32:36,]$Julian), lwd=2, col="orange")

segments(x0=1983, x1=1987,y0=median(m1pred[1:5,]$pred_JD), y1=median(m1pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2014, x1=2018,y0=median(m1pred[32:36,]$pred_JD), y1=median(m1pred[32:36,]$pred_JD), lwd=2, col="red")

