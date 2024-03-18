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

HMS=HMS%>%filter(!Julian<226, !Julian>348, !YR<1990)

#species-level counts
hms_sp=HMS%>%
  pivot_longer(cols=6:21, names_to = "Species", values_to="Count")%>%
  filter(!YR<1990, !Species=="RLHA")

#annual totals per species
hms_sp_tot=hms_sp%>%
  group_by(YR, Species)%>%
  summarise(sp_tot=sum(Count))%>%
  filter(!YR<1990, !Species=="RLHA")%>%
  group_by(YR)%>%
  mutate(rel_abund=sp_tot/sum(sp_tot))

#species relative abundances
hms_sp_ra=hms_sp_tot%>%
  group_by(Species)%>%
  summarise(ave_ct=mean(sp_tot))%>%
  mutate(rel_abund=ave_ct/sum(ave_ct))

#determine minimum year of non-0 detects for all species
hms_sp_min=hms_sp_tot%>%group_by(Species)%>%filter(sp_tot>0)%>%
  summarise(min_yr=min(YR))

#daily totals per species
hms_sp_day=hms_sp%>%
  group_by(Species,YR, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(YR)%>%
  filter(!YR<1990, !Species=="RLHA")

#site-level####

##10%PD####

#10% passage date: species
hms_10pd_sp=hms_sp_day%>%filter(!YR<1990)%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                                           cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                                           PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_10==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#10% passage date: site
hms_10pd_site=HMS%>%filter(!YR<1990)%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                                             cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                                             PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_10==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

#site-level
priors = c(prior(normal(0, 10), class = Intercept),
           prior(normal(0, 10), class = "sds", coef="s(YR)"),
           prior(cauchy(0, 5), class=b, coef="sYR_1"))

hm_mod10=brm(bf(Julian~ 1 + s(YR)), 
          data=hms_10pd_site, iter=2500, family=poisson(link="log"), 
          prior = priors,
          control=list(adapt_delta=0.99))

#saveRDS(hm_mod10, "hms_mod10.RDS")

mod10pred=data.frame(predict(hm_mod10))
colnames(mod10pred)[1]="pred_JD"

m1pred=cbind(hms_10pd_site, mod10pred)

mean(mod10pred[25:29,]$pred_JD)-mean(mod10pred[1:5,]$pred_JD) 
mean(hms_10pd_site[25:29,]$Julian)-mean(hms_10pd_site[1:5,]$Julian) 

##50% PD####
#50% passage date: species
hms_50pd_sp=hms_sp_day%>%filter(!YR<1990)%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                                           cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                                           PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_50==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#50% passage date: site
hms_50pd_site=HMS%>%filter(!YR<1990)%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                                             cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                                             PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_50==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

hm_mod50=brm(bf(Julian~ 1+s(YR)), 
          data=hms_50pd_site, iter=2500, prior=priors,
          family=poisson(), control=list(adapt_delta=0.99))

#saveRDS(hm_mod50, "hms_mod50.RDS")

mod50pred=data.frame(predict(hm_mod50))
colnames(mod50pred)[1]="pred_JD"

m5pred=cbind(hms_50pd_site, mod50pred)

mean(mod50pred[25:29,]$pred_JD)-mean(mod50pred[1:5,]$pred_JD) #-7.9
mean(hms_50pd_site[25:29,]$Julian)-mean(hms_50pd_site[1:5,]$Julian) #-7.9

##90% PD####
#90% passage date: species
hms_90pd_sp=hms_sp_day%>%filter(!YR<1990)%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                                           cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                                           PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_90==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

##90% passage date: site
hms_90pd_site=HMS%>%filter(!YR<1990)%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                                             cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                                             PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_90==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

#site-level
hm_mod90=brm(bf(Julian~ 1+s(YR)), 
             data=hms_90pd_site, iter=2500, prior=priors,
             family=poisson(), control=list(adapt_delta=0.99))

#saveRDS(hm_mod90, "hms_mod90.RDS")
mod90pred=data.frame(predict(hm_mod90))
colnames(mod90pred)[1]="pred_JD"

m9pred=cbind(hms_90pd_site, mod90pred)

mean(mod90pred[25:29,]$pred_JD)-mean(mod90pred[1:5,]$pred_JD) #-7.9
mean(hms_90pd_site[25:29,]$Julian)-mean(hms_90pd_site[1:5,]$Julian) #-7.9

#species-level####

##10% PD####
sp_df_mod10=hms_10pd_sp%>%
  select(Species, YR, Julian)%>%group_by(Species)%>%
  nest()%>%
  mutate(model= purrr::map(data,~ brm(bf(Julian~ 1+s(YR)), 
                         data=., 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))))

#saveRDS(sp_df_mod10, "hms_sp_mod10.RDS")
#only get mean predictions
sp_df_mod10$preds=lapply(sp_df_mod10$model, predict)

sp_df_res10=sp_df_mod10%>%select(-data, -model)%>%
  unnest(cols = c(preds))%>%group_by(Species)%>%
  as.data.frame(sp_df_res10$preds)

hm_mod10_sp_meanpreds=as.data.frame(sp_df_res10$preds[,1])

colnames(hm_mod10_sp_meanpreds)[1]="mean_preds"

main_dat10=hms_10pd_sp%>%select(Species, YR, Julian)

hm_sp_df_all10=cbind(hm_mod10_sp_meanpreds, main_dat10)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))
hm_sp_df_diff10=hm_sp_df_all10%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)%>%
  mutate(diet=case_when(Species=="AMKE" ~ "insect",
                        Species=="BWHA" ~ "mammal",
                        Species=="OSPR" ~ "fish",
                        Species=="MERL" ~ "bird",
                        Species=="BAEA" ~ "fish",
                        Species=="PEFA" ~ "bird",
                        Species=="SSHA" ~ "bird",
                        Species=="COHA" ~ "mammal",
                        Species=="NOHA" ~ "mammal",
                        Species=="BLVU" ~ "carrion",
                        Species=="TUVU" ~ "carrion",
                        Species=="NOGO" ~ "bird",
                        Species=="GOEA" ~ "mammal",
                        Species=="RSHA" ~ "mammal",
                        Species=="RTHA" ~ "mammal"))
  

#relative abundance-weighted shift
hm_sp_diff10=left_join(hm_sp_df_diff10, hms_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

sp10diff=sum(hm_sp_diff10$wgt_shift) #abundance-weighted shift
hm_p10diff=mean(mod10pred[25:29,]$pred_JD)-mean(mod10pred[1:5,]$pred_JD)

sp10diff/hm_p10diff #proportion of phenological shift by species phenology
1-(sp10diff/hm_p10diff) #proportion of phenological shift by composition
hm_p10diff*(1-(sp10diff/hm_p10diff)) #phenological shift (composition)

##50% PD####
sp_df_mod50=hms_50pd_sp%>%
  select(Species, YR, Julian)%>%group_by(Species)%>%
  nest()%>%
  mutate(model= purrr::map(data,~ brm(bf(Julian~ 1+s(YR)), 
                                      data=., 
                                      iter=2500, family=poisson(), 
                                      prior=priors,
                                      control=list(adapt_delta=0.99))))

sp_df_mod50$preds=lapply(sp_df_mod50$model, predict)

#saveRDS(sp_df_mod50, "hms_sp_mod50.RDS")

sp_df_res50=sp_df_mod50%>%select(-data, -model)%>%
  unnest(cols = c(preds))%>%group_by(Species)%>%
  as.data.frame(sp_df_res50$preds)

hm_mod50_sp_meanpreds=as.data.frame(sp_df_res50$preds[,1])

colnames(hm_mod50_sp_meanpreds)[1]="mean_preds"

main_dat50=hms_50pd_sp%>%select(Species, YR, Julian)

hm_sp_df_all50=cbind(hm_mod50_sp_meanpreds, main_dat50)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))

hm_sp_df_diff50=hm_sp_df_all50%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)%>%
  mutate(diet=case_when(Species=="AMKE" ~ "insect",
                        Species=="BWHA" ~ "mammal",
                        Species=="OSPR" ~ "fish",
                        Species=="MERL" ~ "bird",
                        Species=="BAEA" ~ "fish",
                        Species=="PEFA" ~ "bird",
                        Species=="SSHA" ~ "bird",
                        Species=="COHA" ~ "mammal",
                        Species=="NOHA" ~ "mammal",
                        Species=="BLVU" ~ "carrion",
                        Species=="TUVU" ~ "carrion",
                        Species=="NOGO" ~ "bird",
                        Species=="GOEA" ~ "mammal",
                        Species=="RSHA" ~ "mammal",
                        Species=="RTHA" ~ "mammal"))

#relative abundance-weighted shift
hm_sp_diff50=left_join(hm_sp_df_diff50, hms_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

sp50diff=sum(hm_sp_diff50$wgt_shift) #abundance-weighted shift
hm_p50diff=mean(mod50pred[25:29,]$pred_JD)-mean(mod50pred[1:5,]$pred_JD)

sp50diff/hm_p50diff #proportion of phenological shift by species phenology
1-(sp50diff/hm_p50diff) #proportion of phenological shift by composition
hm_p50diff*(1-(sp50diff/hm_p50diff)) #phenological shift (composition)

##90% PD####
sp_df_mod90=hms_90pd_sp%>%
  select(Species, YR, Julian)%>%group_by(Species)%>%
  nest()%>%
  mutate(model= purrr::map(data,~ brm(bf(Julian~ 1+s(YR)), 
                                      data=., 
                                      iter=2500, family=poisson(), 
                                      prior=priors,
                                      control=list(adapt_delta=0.99))))

sp_df_mod90$preds=lapply(sp_df_mod90$model, predict)
#saveRDS(sp_df_mod90, "hms_sp_mod90.RDS")

sp_df_res90=sp_df_mod90%>%select(-data, -model)%>%
  unnest(cols = c(preds))%>%group_by(Species)%>%
  as.data.frame(sp_df_res90$preds)

hm_mod90_sp_meanpreds=as.data.frame(sp_df_res90$preds[,1])

colnames(hm_mod90_sp_meanpreds)[1]="mean_preds"

main_dat90=hms_90pd_sp%>%select(Species, YR, Julian)

hm_sp_df_all90=cbind(hm_mod90_sp_meanpreds, main_dat90)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))

hm_sp_df_diff90=hm_sp_df_all90%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)%>%
  mutate(diet=case_when(Species=="AMKE" ~ "insect",
                        Species=="BWHA" ~ "mammal",
                        Species=="OSPR" ~ "fish",
                        Species=="MERL" ~ "bird",
                        Species=="BAEA" ~ "fish",
                        Species=="PEFA" ~ "bird",
                        Species=="SSHA" ~ "bird",
                        Species=="COHA" ~ "mammal",
                        Species=="NOHA" ~ "mammal",
                        Species=="BLVU" ~ "carrion",
                        Species=="TUVU" ~ "carrion",
                        Species=="NOGO" ~ "bird",
                        Species=="GOEA" ~ "mammal",
                        Species=="RSHA" ~ "mammal",
                        Species=="RTHA" ~ "mammal"))


#relative abundance-weighted shift
hm_sp_diff90=left_join(hm_sp_df_diff90, hms_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

sp90diff=sum(hm_sp_diff90$wgt_shift) #abundance-weighted shift
hm_p90diff=mean(mod90pred[25:29,]$pred_JD)-mean(mod90pred[1:5,]$pred_JD)

sp90diff/hm_p90diff #proportion of phenological shift by species phenology
1-(sp90diff/hm_p90diff) #proportion of phenological shift by composition
hm_p90diff*(1-(sp90diff/hm_p90diff)) #phenological shift (composition)

#compositional shifts####
hm_sp_abundance=hms_sp_tot%>%
  group_by(Species)%>%nest()%>%
  mutate(model= purrr::map(data, ~brm(bf(rel_abund~ 1+s(YR)),
                                      data= .,
                                      iter=2500, 
                                      family=Beta(link = "logit", link_phi = "log"),
                                      prior=priors,
                                      control=list(adapt_delta=0.99))))

#saveRDS(hm_sp_abundance,"HMS_composition.RDS")

hm_sp_abundance$preds=lapply(hm_sp_abundance$model,predict)

hm_sp_abundance_preds=hm_sp_abundance%>%select(-data, -model)%>%
  unnest(cols = c(preds))%>%group_by(Species)%>%
  as.data.frame(preds)

YR=rep(1990:2018, each=1, times=15)

hm_sp_abundance_meanpreds=data.frame(hm_sp_abundance_preds$Species,hm_sp_abundance_preds$preds[,1])%>%
  cbind(YR)

colnames(hm_sp_abundance_meanpreds)[2]="mean_preds"
colnames(hm_sp_abundance_meanpreds)[1]="Species"

main_dat=hms_sp_tot%>%select(Species, YR, rel_abund)

hm_sp_abundance_diffs=inner_join(hm_sp_abundance_meanpreds, main_dat, by=c("Species", "YR"))%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))%>%
  filter(!is.na(period))%>%
  mutate(diet=case_when(Species=="AMKE" ~ "insect",
                        Species=="BWHA" ~ "mammal",
                        Species=="OSPR" ~ "fish",
                        Species=="MERL" ~ "bird",
                        Species=="BAEA" ~ "fish",
                        Species=="PEFA" ~ "bird",
                        Species=="SSHA" ~ "bird",
                        Species=="COHA" ~ "mammal",
                        Species=="NOHA" ~ "mammal",
                        Species=="BLVU" ~ "carrion",
                        Species=="TUVU" ~ "carrion",
                        Species=="NOGO" ~ "bird",
                        Species=="GOEA" ~ "mammal",
                        Species=="RSHA" ~ "mammal",
                        Species=="RTHA" ~ "mammal"))

hm_sp_abundance_diffs=hm_sp_abundance_diffs%>%
  select(Species, mean_preds,period, diet)%>%
  filter(!is.na(period))%>%
  group_by(Species, diet,period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)


#posterior predictions####
##site:10% PD####
hpreds10=as.data.frame(posterior_predict(hm_mod10))%>%
  mutate(iter=seq(1:5000))%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

quantile(hpreds10$shift, probs=c(0.025,0.975))
mean(hpreds10$shift)

length(which(hpreds10$shift<0))/length(hpreds10$shift)

##site:50% PD####
hpreds50=as.data.frame(posterior_predict(hm_mod50))%>%
  mutate(iter=seq(1:5000))%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

quantile(hpreds50$shift, probs=c(0.025,0.975))
mean(hpreds50$shift)

length(which(hpreds50$shift<0))/length(hpreds50$shift)

##site:90% PD####
hpreds90=as.data.frame(posterior_predict(hm_mod90))%>%
  mutate(iter=seq(1:5000))%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

quantile(hpreds90$shift, probs=c(0.025,0.975))
mean(hpreds90$shift)

length(which(hpreds90$shift<0))/length(hpreds90$shift)

##species: 10% PD####
spph=sp_df_mod10%>%
  mutate(preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

hms_splist=spph%>%select(Species, iter)

spph2=as.data.frame(spph$preds[,])

spph3=cbind(spph2, hms_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

spph_quant=spph3%>%
  group_by(Species)%>%
  reframe(
    q0=quantile(shift, 0.5),
    q1=quantile(shift, 0.025),
    q2=quantile(shift, 0.975))

spph_prob=spph3%>%
  group_by(Species)%>%
  reframe(
    prob_shift_advance=length(which(shift<0))/ length(shift),
    prob_shift_delay=length(which(shift>0))/ length(shift))

##species: 50% PD####
spph50=sp_df_mod50%>%
  mutate(preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

spph50_df=as.data.frame(spph50$preds[,])

spph50_dat=cbind(spph50_df, hms_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

spph50_quant=spph50_dat%>%
  group_by(Species)%>%
  reframe(
    q50=quantile(shift, 0.5),
    qLB=quantile(shift, 0.025),
    qUB=quantile(shift, 0.975))

spph50_prob=spph50_dat%>%
  group_by(Species)%>%
  reframe(
    prob_shift_advance=length(which(shift<0))/ length(shift),
    prob_shift_delay=length(which(shift>0))/ length(shift))

##species: 90% PD####
spph90=sp_df_mod90%>%
  mutate(preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

spph90_df=as.data.frame(spph90$preds[,])

spph90_dat=cbind(spph90_df, hms_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

spph90_quant=spph90_dat%>%
  group_by(Species)%>%
  reframe(
    q50=quantile(shift, 0.5),
    qLB=quantile(shift, 0.025),
    qUB=quantile(shift, 0.975))

spph90_prob=spph90_dat%>%
  group_by(Species)%>%
  reframe(
    prob_shift_advance=length(which(shift<0))/ length(shift),
    prob_shift_delay=length(which(shift>0))/ length(shift))

##composition####
hm_abund=hm_sp_abundance%>%
  mutate(preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()%>%select(-data, -model)%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))


hm_abund2=as.data.frame(hm_abund$preds[,])

hm_abund3=cbind(hm_abund2, hms_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

hm_abund_quant=hm_abund3%>%
  group_by(Species)%>%
  reframe(
    q0=quantile(shift, 0.5),
    q1=quantile(shift, 0.025),
    q2=quantile(shift, 0.975))

hm_abund_prob=hm_abund3%>%
  group_by(Species)%>%
  reframe(
    prob_shift_increase=length(which(shift<0))/ length(shift),
    prob_shift_decrease=length(which(shift>0))/ length(shift))
