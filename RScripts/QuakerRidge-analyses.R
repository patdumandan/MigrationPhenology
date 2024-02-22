require(dplyr)
require(tidyr)
require(lubridate)
require(ggplot2)
require(broom)
require(brms)

#data manipulation####
qr=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/QuakerRidge.csv")

qr$Date=lubridate::mdy(qr$Date)

qk=qr%>%
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
         AMKE, BLVU, BWHA, COHA, GOEA, MERL, NOHA, PEFA,
          RTHA)%>%
  mutate(TOTAL=rowSums(.[7:15], na.rm=T))

qk$Julian=as.POSIXlt(qk$Date)$yday

QR=qk%>%filter(!year<2000, !Julian<231, !Julian >323)

#species-level counts#
qr_sp=QR%>%select(-TOTAL)%>%
  pivot_longer(cols=7:15, names_to = "Species", values_to="Count")

#annual totals per species
qr_sp_tot=qr_sp%>%
  group_by(year, Species)%>%
  summarise(sp_tot=sum(Count))

ggplot(qr_sp_tot)+geom_line(aes(x=year, y=log(sp_tot)))+facet_wrap(~Species)+theme_classic()

#species relative abundances
qr_sp_ra=qr_sp_tot%>%
  group_by(Species)%>%
  summarise(ave_ct=mean(sp_tot))%>%
  mutate(rel_abund=ave_ct/sum(ave_ct))%>%
  filter(!Species%in%c("SW", "RL", "SE", "MK"))
 
#determine uncommon species (based on non-detects annually)
qr_sp_min=qr_sp_tot%>%group_by(Species)%>%filter(sp_tot==min(sp_tot))

#daily totals per species
qr_sp_day=qr_sp%>%
  group_by(Species,year, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(year)

#site-level####

#10% passage date: site
qr_10pd_site=QR%>%group_by(year)%>%mutate(cum_tot=cumsum(TOTAL), #get cumulative sum for each year
                                          cum_pt=cum_tot/sum(TOTAL),#calculate % of daily count/annual count
                                          PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >10%
  group_by(year)%>%filter(PD_10==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(year, Julian, cum_pt)

#10% passage date: species
qr_10pd_sp=qr_sp_day%>%group_by(Species, year)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                      cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                      PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, year)%>%filter(PD_10==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

##10% PD####
priors = c(prior(normal(0, 10), class = Intercept),
           prior(normal(0, 10), class = "sds", coef="s(year)"),
           prior(cauchy(0, 5), class=b, coef="syear_1"))

qr_mod10=brm(bf(Julian~ 1 + s(year)), 
             data=qr_10pd_site, iter=2500, family=poisson(link="log"), 
             prior = priors,
             control=list(adapt_delta=0.99))

qr_mod10pred=data.frame(predict(qr_mod10))
colnames(qr_mod10pred)[1]="pred_JD"

qr_m1pred=cbind(qr_10pd_site, qr_mod10pred)

mean(qr_mod10pred[19:23,]$pred_JD)-mean(qr_mod10pred[1:5,]$pred_JD) 
mean(qr_10pd_site[19:23,]$Julian)-mean(qr_10pd_site[1:5,]$Julian) 

#50% passage date
qr_50pd_site=QR%>%group_by(year)%>%mutate(cum_tot=cumsum(TOTAL), #get cumulative sum for each year
                                          cum_pt=cum_tot/sum(TOTAL),#calculate % of daily count/annual count
                                          PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >10%
  group_by(year)%>%filter(PD_50==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(year, Julian, cum_pt)

#50% passage date: species
qr_50pd_sp=qr_sp_day%>%group_by(Species, year)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                        cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                        PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, year)%>%filter(PD_50==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

##50% PD####
qr_mod50=brm(bf(Julian~ 1 + s(year)), 
             data=qr_50pd_site, iter=2500, family=poisson(link="log"), 
             prior = priors,
             control=list(adapt_delta=0.99))

qr_mod50pred=data.frame(predict(qr_mod50))
colnames(qr_mod50pred)[1]="pred_JD"

qr_m5pred=cbind(qr_50pd_site, qr_mod50pred)

mean(qr_mod50pred[19:23,]$pred_JD)-mean(qr_mod50pred[1:5,]$pred_JD) 
mean(qr_50pd_site[19:23,]$Julian)-mean(qr_50pd_site[1:5,]$Julian) 

#90% passage date
qr_90pd_site=QR%>%group_by(year)%>%mutate(cum_tot=cumsum(TOTAL), #get cumulative sum for each year
                                          cum_pt=cum_tot/sum(TOTAL),#calculate % of daily count/annual count
                                          PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >10%
  group_by(year)%>%filter(PD_90==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(year, Julian, cum_pt)

#90% passage date: species
qr_90pd_sp=qr_sp_day%>%group_by(Species, year)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                        cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                        PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, year)%>%filter(PD_90==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

##90% PD####
qr_mod90=brm(bf(Julian~ 1 + s(year)), 
             data=qr_90pd_site, iter=2500, family=poisson(link="log"), 
             prior = priors,
             control=list(adapt_delta=0.99))

qr_mod90pred=data.frame(predict(qr_mod90))
colnames(qr_mod90pred)[1]="pred_JD"

qr_m9pred=cbind(qr_90pd_site, qr_mod90pred)

mean(qr_mod90pred[19:23,]$pred_JD)-mean(qr_mod90pred[1:5,]$pred_JD) 
mean(qr_90pd_site[19:23,]$Julian)-mean(qr_90pd_site[1:5,]$Julian) 

#species-level####

##10% PD####
qr_sp_df_mod10=qr_10pd_sp%>%
  select(Species, year, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(year)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

qr_sp_df_res10=qr_sp_df_mod10%>%select(-data, -model)
qr_sp_df_res10preds=as.data.frame(qr_sp_df_res10$preds[,1])
colnames(qr_sp_df_res10preds)[1]="mean_preds"

qr_main_dat10=qr_10pd_sp%>%select(Species, year, Julian)

qr_sp_df_all=cbind(qr_sp_df_res10preds, qr_main_dat10)%>%as.data.frame()%>%
  mutate(period=case_when(year%in%c(2000:2004) ~ "T1",
                          year%in%c(2018:2022) ~ "T2"))
qr_sp_df_diff=qr_sp_df_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
qr_sp_diff10=left_join(qr_sp_df_diff, qr_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

qr_sp10diff=sum(qr_sp_diff10$wgt_shift) #abundance-weighted shift
qr_p10diff=mean(qr_mod10pred[19:23,]$pred_JD)-mean(qr_mod10pred[1:5,]$pred_JD)

qr_sp10diff/qr_p10diff #proportion of phenological shift by species phenology
1-(qr_sp10diff/qr_p10diff) #proportion of phenological shift by composition
qr_p10diff*(1-(qr_sp10diff/qr_p10diff)) #phenological shift (composition)

##50% PD####
qr_sp_df_mod50=qr_50pd_sp%>%
  select(Species, year, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(year)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

qr_sp_df_res50=qr_sp_df_mod50%>%select(-data, -model)
qr_sp_df_res50preds=as.data.frame(qr_sp_df_res50$preds[,1])
colnames(qr_sp_df_res50preds)[1]="mean_preds"

qr_main_dat50=qr_50pd_sp%>%select(Species, year, Julian)

qr_sp_df_all=cbind(qr_sp_df_res50preds, qr_main_dat50)%>%as.data.frame()%>%
  mutate(period=case_when(year%in%c(2000:2004) ~ "T1",
                          year%in%c(2018:2022) ~ "T2"))
qr_sp_df_diff=qr_sp_df_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
qr_sp_diff50=left_join(qr_sp_df_diff, qr_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

qr_sp50diff=sum(qr_sp_diff50$wgt_shift) #abundance-weighted shift
qr_p50diff=mean(qr_mod50pred[19:23,]$pred_JD)-mean(qr_mod50pred[1:5,]$pred_JD)

qr_sp50diff/qr_p50diff #proportion of phenological shift by species phenology
1-(qr_sp50diff/qr_p50diff) #proportion of phenological shift by composition
qr_p50diff*(1-(qr_sp50diff/qr_p50diff)) #phenological shift (composition)

##90% PD####
qr_sp_df_mod90=qr_90pd_sp%>%
  select(Species, year, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(year)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

qr_sp_df_res90=qr_sp_df_mod90%>%select(-data, -model)
qr_sp_df_res90preds=as.data.frame(qr_sp_df_res90$preds[,1])
colnames(qr_sp_df_res90preds)[1]="mean_preds"

qr_main_dat90=qr_90pd_sp%>%select(Species, year, Julian)

qr_sp_df_all=cbind(qr_sp_df_res90preds, qr_main_dat90)%>%as.data.frame()%>%
  mutate(period=case_when(year%in%c(2000:2004) ~ "T1",
                          year%in%c(2018:2022) ~ "T2"))
qr_sp_df_diff=qr_sp_df_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
qr_sp_diff90=left_join(qr_sp_df_diff, qr_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

qr_sp90diff=sum(qr_sp_diff90$wgt_shift) #abundance-weighted shift
qr_p90diff=mean(qr_mod90pred[19:23,]$pred_JD)-mean(qr_mod90pred[1:5,]$pred_JD)

qr_sp90diff/qr_p90diff #proportion of phenological shift by species phenology
1-(qr_sp90diff/qr_p90diff) #proportion of phenological shift by composition
qr_p90diff*(1-(qr_sp90diff/qr_p90diff)) #phenological shift (composition)

#visualize data####
##shifts####

##10% PD####
plot(qr_10pd_site$Julian~qr_10pd_site$year, type="l", main="Quaker Ridge", ylab="10% passage date", xlab="year")
lines(qr_m1pred$pred_JD~qr_m1pred$year, lty=2, col="grey", lwd=2)

segments(x0=2000, x1=2004,y0=mean(qr_10pd_site[1:5,]$Julian), y1=mean(qr_10pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2018, x1=2022,y0=mean(qr_10pd_site[19:23,]$Julian), y1=mean(qr_10pd_site[19:23,]$Julian), lwd=2, col="orange")

segments(x0=2000, x1=2004,y0=mean(qr_mod10pred[1:5,]$pred_JD), y1=mean(qr_mod10pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2018, x1=2022,y0=mean(qr_mod10pred[19:23,]$pred_JD), y1=mean(qr_mod10pred[19:23,]$pred_JD), lwd=2, col="red")

##50% PD####
plot(qr_50pd_site$Julian~qr_50pd_site$year, type="l", main="Quaker Ridge", ylab="50% passage date", xlab="year")
lines(qr_m5pred$pred_JD~qr_m5pred$year, lty=2, col="grey", lwd=2)

segments(x0=2000, x1=2004,y0=mean(qr_50pd_site[1:5,]$Julian), y1=mean(qr_50pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2018, x1=2022,y0=mean(qr_50pd_site[19:23,]$Julian), y1=mean(qr_50pd_site[19:23,]$Julian), lwd=2, col="orange")

segments(x0=2000, x1=2004,y0=mean(qr_mod50pred[1:5,]$pred_JD), y1=mean(qr_mod50pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2018, x1=2022,y0=mean(qr_mod50pred[19:23,]$pred_JD), y1=mean(qr_mod50pred[19:23,]$pred_JD), lwd=2, col="red")

##90% PD####
plot(qr_90pd_site$Julian~qr_90pd_site$year, type="l", main="Quaker Ridge", ylab="90% passage date", xlab="year")
lines(qr_m9pred$pred_JD~qr_m9pred$year, lty=2, col="grey", lwd=2)

segments(x0=2000, x1=2004,y0=mean(qr_90pd_site[1:5,]$Julian), y1=mean(qr_90pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2018, x1=2022,y0=mean(qr_90pd_site[19:23,]$Julian), y1=mean(qr_90pd_site[19:23,]$Julian), lwd=2, col="orange")

segments(x0=2000, x1=2004,y0=mean(qr_mod90pred[1:5,]$pred_JD), y1=mean(qr_mod90pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2018, x1=2022,y0=mean(qr_mod90pred[19:23,]$pred_JD), y1=mean(qr_mod90pred[19:23,]$pred_JD), lwd=2, col="red")

#contributions####
qr_datphen=matrix(c(10,50,90,
                    1.25,-1.28,14.70, 
                    2.69,-1.73,-1.91,
                    -1.43,0.45 ,16.60 ), nrow=3, ncol=4, byrow=FALSE, 
                  dimnames= list(c("1", "2", "3"),
                                 c("metric","overall", "species phenology", "composition")))
qr_dat_phen=as.data.frame(qr_datphen)%>%
  pivot_longer(cols=2:4, names_to = "shift", values_to="value" )

ggplot(qr_dat_phen)+geom_col(aes(x=shift, y=value))+ggtitle("Quaker Ridge")+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_x_discrete(limits=c("overall", "species phenology", "composition"))+
  facet_wrap(~metric)+theme_classic()+abline(h=0)+ylab("shift (days)")

#compositional shifts####
qr_sp_tot=qr_sp%>%
  group_by(year, Species)%>%
  summarise(sp_tot=sum(Count))%>%
  group_by(year)%>%
  mutate(rel_abund=sp_tot/sum(sp_tot))

qr_sp_abundance2=qr_sp_tot%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(rel_abund~ 1+s(year)), 
                         data=data, 
                         iter=2500, 
                         family=zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit"),
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

qr_sp_relab=qr_sp_abundance2%>%select(-data, -model)
qr_sp_rela_preds=as.data.frame(qr_sp_relab$preds[,1])
colnames(qr_sp_rela_preds)[1]="mean_rela"

qr_sp_df_rela=cbind(qr_sp_rela_preds, qr_sp_tot)%>%as.data.frame()%>%
  mutate(period=case_when(year%in%c(2002:2004) ~ "T1",
                          year%in%c(2018:2022) ~ "T2"))

qr_sp_df_rela_diff=qr_sp_df_rela%>%
  select(Species, mean_rela,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_rela=mean(mean_rela))%>%
  pivot_wider(names_from = period, values_from=mean_rela, names_sep = ".")%>%
  mutate(shift= T2-T1)

#determine passage order of species (based on 50% passage date)
qr_ord=qr_50pd_sp%>%group_by(Species)%>%
  summarise(min_time=min(Julian))%>%arrange(min_time)

ggplot(qr_sp_df_rela_diff)+geom_col(aes(x=Species, y=shift))+theme_classic()+
  geom_hline(yintercept = 0, linetype = 2)+ggtitle("Quaker Ridge")+
  ylab("relative abundance shifts")+
  scale_x_discrete(limits=c("COHA","BWHA","PEFA", "BLVU", "AMKE", "MERL",
                            "NOHA", "GOEA", "RTHA"))

