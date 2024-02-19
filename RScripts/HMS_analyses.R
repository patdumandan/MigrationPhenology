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
           prior(normal(0, 10), class = "sds", coef="s(YR)"),
           prior(cauchy(0, 5), class=b, coef="sYR_1"))

mod10=brm(bf(Julian~ 1 + s(YR)), 
          data=hms_10pd_site, iter=2500, family=poisson(link="log"), 
          prior = priors,
          control=list(adapt_delta=0.99))

mod10pred=data.frame(predict(mod10))
colnames(mod10pred)[1]="pred_JD"

m1pred=cbind(hms_10pd_site, mod10pred)

mean(mod10pred[32:36,]$pred_JD)-mean(mod10pred[1:5,]$pred_JD) 
mean(hms_10pd_site[32:36,]$Julian)-mean(hms_10pd_site[1:5,]$Julian) 

##50% PD####
#50% passage date: species
hms_50pd_sp=hms_sp_day%>%filter(!YR<1983)%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                                           cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                                           PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_50==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#50% passage date: site
hms_50pd_site=HMS%>%filter(!YR<1983)%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                                             cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                                             PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_50==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

mod50=brm(bf(Julian~ 1+s(YR)), 
          data=hms_50pd_site, iter=2500, prior=priors,
          family=poisson(), control=list(adapt_delta=0.99))

mod50pred=data.frame(predict(mod50))
colnames(mod50pred)[1]="pred_JD"

m5pred=cbind(hms_50pd_site, mod50pred)

mean(mod50pred[32:36,]$pred_JD)-mean(mod50pred[1:5,]$pred_JD) #-7.9
mean(hms_50pd_site[32:36,]$Julian)-mean(hms_50pd_site[1:5,]$Julian) #-7.9

##90% PD####
#90% passage date: species
hms_90pd_sp=hms_sp_day%>%filter(!YR<1983)%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                                           cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                                           PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_90==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

##90% passage date: site
hms_90pd_site=HMS%>%filter(!YR<1983)%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                                             cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                                             PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_90==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

#site-level
mod90=brm(bf(Julian~ 1+s(YR)), 
          data=hms_90pd_site, iter=2500, 
          family=poisson(), 
          prior=priors,
          control=list(adapt_delta=0.99))

mod90pred=data.frame(predict(mod90))
colnames(mod90pred)[1]="pred_JD"

m9pred=cbind(hms_90pd_site, mod90pred)

mean(mod90pred[32:36,]$pred_JD)-mean(mod90pred[1:5,]$pred_JD) #-7.9
mean(hms_90pd_site[32:36,]$Julian)-mean(hms_90pd_site[1:5,]$Julian) #-7.9

#species-level####

##10% PD####
sp_df_mod10=hms_10pd_sp%>%
  select(Species, YR, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(YR)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

sp_df_res10=sp_df_mod10%>%select(-data, -model)
sp_df_res10preds=as.data.frame(sp_df_res10$preds[,1])
colnames(sp_df_res10preds)[1]="mean_preds"

main_dat10=hms_10pd_sp%>%select(Species, YR, Julian)

sp_df_all=cbind(sp_df_res10preds, main_dat10)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1983:1987) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))
sp_df_diff=sp_df_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
sp_diff10=left_join(sp_df_diff, hms_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

sp10diff=sum(sp_diff10$wgt_shift) #abundance-weighted shift
p10diff=mean(mod10pred[32:36,]$pred_JD)-mean(mod10pred[1:5,]$pred_JD)

sp10diff/p10diff #proportion of phenological shift by species phenology
1-(sp10diff/p10diff) #proportion of phenological shift by composition
p10diff*(1-(sp10diff/p10diff)) #phenological shift (composition)

##50% PD####
sp_df_mod50=hms_50pd_sp%>%
  select(Species, YR, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(YR)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

sp_df_res50=sp_df_mod50%>%select(-data, -model)
sp_df_res50preds=as.data.frame(sp_df_res50$preds[,1])
colnames(sp_df_res50preds)[1]="mean_preds"

main_dat50=hms_50pd_sp%>%select(Species, YR, Julian)

sp_df50_all=cbind(sp_df_res50preds, main_dat50)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1983:1987) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))
sp_df50_diff=sp_df50_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
sp_diff50=left_join(sp_df50_diff, hms_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

sp50diff=sum(sp_diff50$wgt_shift)
p50diff=mean(mod50pred[32:36,]$pred_JD)-mean(mod50pred[1:5,]$pred_JD)

sp50diff/p50diff
1-(sp50diff/p50diff)
p50diff*(1-(sp50diff/p50diff)) #proportion of phenological shift by composition

##90% PD####
sp_df_mod90=hms_90pd_sp%>%
  select(Species, YR, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(YR)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

sp_df_res90=sp_df_mod90%>%select(-data, -model)
sp_df_res90preds=as.data.frame(sp_df_res90$preds[,1])
colnames(sp_df_res90preds)[1]="mean_preds"

main_dat90=hms_90pd_sp%>%select(Species, YR, Julian)

sp_df90_all=cbind(sp_df_res90preds, main_dat90)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1983:1987) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))
sp_df90_diff=sp_df90_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
sp_diff90=left_join(sp_df90_diff, hms_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

sp90diff=sum(sp_diff90$wgt_shift)
p90diff=mean(mod90pred[32:36,]$pred_JD)-mean(mod90pred[1:5,]$pred_JD)

sp90diff/p90diff
sp90diff/p90diff #proportion of phenological shift by species phenology
1-(sp90diff/p90diff) #proportion of phenological shift by composition
p90diff*(1-(sp90diff/p90diff)) #phenological shift (composition)

#visualize data####

##shifts####
##10% PD####
plot(hms_10pd_site$Julian~hms_10pd_site$YR, type="l", main="Hawk Mountain", ylab="10% passage date", xlab="year")

lines(m1pred$pred_JD~m1pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1983, x1=1987,y0=mean(hms_10pd_site[1:5,]$Julian), y1=mean(hms_10pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2014, x1=2018,y0=mean(hms_10pd_site[32:36,]$Julian), y1=mean(hms_10pd_site[32:36,]$Julian), lwd=2, col="orange")

segments(x0=1983, x1=1987,y0=mean(m1pred[1:5,]$pred_JD), y1=mean(m1pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2014, x1=2018,y0=mean(m1pred[32:36,]$pred_JD), y1=mean(m1pred[32:36,]$pred_JD), lwd=2, col="red")

##50% PD####
plot(hms_50pd_site$Julian~hms_50pd_site$YR, type="l", ylab="50% passage date", xlab="year")

lines(m5pred$pred_JD~m5pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1983, x1=1987,y0=mean(hms_50pd_site[1:5,]$Julian), y1=mean(hms_50pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2014, x1=2018,y0=mean(hms_50pd_site[32:36,]$Julian), y1=mean(hms_50pd_site[32:36,]$Julian), lwd=2, col="orange")

segments(x0=1983, x1=1987,y0=mean(m5pred[1:5,]$pred_JD), y1=mean(m5pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2014, x1=2018,y0=mean(m5pred[32:36,]$pred_JD), y1=mean(m5pred[32:36,]$pred_JD), lwd=2, col="red")

##90% PD####
plot(hms_90pd_site$Julian~hms_90pd_site$YR, type="l", ylab="90% passage date", xlab="year")

lines(m9pred$pred_JD~m9pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1983, x1=1987,y0=mean(hms_90pd_site[1:5,]$Julian), y1=mean(hms_90pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2014, x1=2018,y0=mean(hms_90pd_site[32:36,]$Julian), y1=mean(hms_90pd_site[32:36,]$Julian), lwd=2, col="orange")

segments(x0=1983, x1=1987,y0=mean(m9pred[1:5,]$pred_JD), y1=mean(m9pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2014, x1=2018,y0=mean(m9pred[32:36,]$pred_JD), y1=mean(m9pred[32:36,]$pred_JD), lwd=2, col="red")


#contributions####
datphen=matrix(c(10,50,90,
                 -0.44,-6.81,-5.28, 
                 0.22,3.36,2.50,
                 -0.66, -10.17, -7.79), nrow=3, ncol=4, byrow=FALSE, 
               dimnames= list(c("1", "2", "3"),
                              c("metric","overall", "species phenology", "composition")))
dat_phen=as.data.frame(datphen)%>%
  pivot_longer(cols=2:4, names_to = "shift", values_to="value" )

ggplot(dat_phen)+geom_col(aes(x=shift, y=value))+ggtitle("Hawk Mountain")+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_x_discrete(limits=c("overall", "species phenology", "composition"))+
  facet_wrap(~metric)+theme_classic()+abline(h=0)+ylab("shift (days)")
