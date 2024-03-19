require(dplyr)
require(tidyr)
require(lubridate)
require(ggplot2)
require(broom)
require(brms)

#data manipulation####
cape=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/CapeMay.csv")

CapeMay=cape%>%
  select(-Duration, -Observer,-RL,-GO,-MK,-SE,-SW,-UA, -UB, -UF,  -ZT, -TOTAL)%>%
  mutate(yr_tot=rowSums(.[2:16], na.rm=T),
         date=as.Date(Date, format="%m/%d/%Y"),
         YR=year(date), MO=month(date), DAY=day(date))%>%
  filter(!YR<1990, !YR>2018)

CapeMay$Julian=as.POSIXlt(CapeMay$date)$yday
CapeMay[, 2:16][is.na(CapeMay[, 2:16])] <- 0
CapeMay=CapeMay%>%filter(!Julian<244, !Julian>333)

#species-level counts
cm_sp=CapeMay%>%select(-yr_tot)%>%
  pivot_longer(cols=2:16, names_to = "Species", values_to="Count")

#annual totals per species
cm_sp_tot=cm_sp%>%
  group_by(YR, Species)%>%
  summarise(sp_tot=sum(Count))%>%
   group_by(YR)%>%
  mutate(rel_abund=sp_tot/sum(sp_tot))

ggplot(cm_sp_tot)+geom_line(aes(x=YR, y=log(sp_tot)))+facet_wrap(~Species)+theme_classic()

#species relative abundances
cm_sp_ra=cm_sp_tot%>%
  group_by(Species)%>%
  summarise(ave_ct=mean(sp_tot))%>%
  mutate(rel_abund=ave_ct/sum(ave_ct))

#determine minimum year of non-0 detects for all species
#use to determine which species to exclude
cm_sp_min=cm_sp_tot%>%group_by(Species)%>%filter(sp_tot>0)%>%summarise(n_yr=length(YR))

#daily totals per species
cm_sp_day=cm_sp%>%
  group_by(Species,YR, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(YR)

cm_sp_freq=cm_sp_day%>% 
  group_by(Julian)%>%
  summarise(day_freq=n())%>%
  arrange(Julian)

#site-level####

##10%PD####

#10% passage date: species
cm_10pd_sp=cm_sp_day%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                      cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                      PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_10==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#10% passage date: site
cm_10pd_site=CapeMay%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                            cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                            PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_10==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

#site-level
priors = c(prior(normal(0, 10), class = Intercept),
           prior(normal(0, 10), class = "sds", coef="s(YR)"),
           prior(cauchy(0, 5), class=b, coef="sYR_1"))

cm_mod10=brm(bf(Julian~ 1 + s(YR)), 
             data=cm_10pd_site, iter=2500, family=poisson(link="log"), 
             prior = priors,
             control=list(adapt_delta=0.99))

#saveRDS(cm_mod10, "CapeMay_mod10.RDS")
cm_mod10pred=data.frame(predict(cm_mod10))
colnames(cm_mod10pred)[1]="pred_JD"

cm_m1pred=cbind(cm_10pd_site, cm_mod10pred)

mean(cm_mod10pred[25:29,]$pred_JD)-mean(cm_mod10pred[1:5,]$pred_JD) 
mean(cm_10pd_site[25:29,]$Julian)-mean(cm_10pd_site[1:5,]$Julian) 

##50% PD####
#50% passage date: species
cm_50pd_sp=cm_sp_day%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                      cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                      PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_50==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#50% passage date: site
cm_50pd_site=CapeMay%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                            cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                            PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_50==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

cm_mod50=brm(bf(Julian~ 1+s(YR)), 
             data=cm_50pd_site, iter=2500, prior=priors,
             family=poisson(), control=list(adapt_delta=0.99))

#saveRDS(cm_mod50, "CM_mod50.RDS")
cm_mod50pred=data.frame(predict(cm_mod50))
colnames(cm_mod50pred)[1]="pred_JD"

cm_m5pred=cbind(cm_50pd_site, cm_mod50pred)

mean(cm_mod50pred[25:29,]$pred_JD)-mean(cm_mod50pred[1:5,]$pred_JD) #-7.9
mean(cm_50pd_site[25:29,]$Julian)-mean(cm_50pd_site[1:5,]$Julian) #-7.9

##90% PD####
#90% passage date: species
cm_90pd_sp=cm_sp_day%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                      cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                      PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_90==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

##90% passage date: site
cm_90pd_site=CapeMay%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                            cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                            PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_90==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

#site-level
cm_mod90=brm(bf(Julian~ 1+s(YR)), 
             data=cm_90pd_site, iter=2500, 
             family=poisson(), 
             prior=priors,
             control=list(adapt_delta=0.99))

#saveRDS(cm_mod90, "CM_mod90.RDS")

cm_mod90pred=data.frame(predict(cm_mod90))
colnames(cm_mod90pred)[1]="pred_JD"

cm_m9pred=cbind(cm_90pd_site, cm_mod90pred)

mean(cm_mod90pred[25:29,]$pred_JD)-mean(cm_mod90pred[1:5,]$pred_JD) #-7.9
mean(cm_90pd_site[25:29,]$Julian)-mean(cm_90pd_site[1:5,]$Julian) #-7.9

#species-level####

##10% PD####
cm_sp_df_mod10=cm_10pd_sp%>%
  select(Species, YR, Julian)%>%group_by(Species)%>%
  nest()%>%
  mutate(model= purrr::map(data,~ brm(bf(Julian~ 1+s(YR)), 
                                      data=., 
                                      iter=2500, family=poisson(), 
                                      prior=priors,
                                      control=list(adapt_delta=0.99))))

#saveRDS(cm_sp_df_mod10, "CM_sp_mod10.RDS")

cm_sp_df_mod10$preds=lapply(cm_sp_df_mod10$model, predict)

cm_sp_df_res10=cm_sp_df_mod10%>%select(-data, -model)%>%
  unnest(cols = c(preds))%>%group_by(Species)%>%filter(!Species=="NG")%>%
  as.data.frame(cm_sp_df_res10$preds)

cm_mod10_sp_meanpreds=as.data.frame(cm_sp_df_res10$preds[,1])

colnames(cm_mod10_sp_meanpreds)[1]="mean_preds"

cmng_sp_df_res10=cm_sp_df_mod10%>%select(-data, -model)%>%
  unnest(cols = c(preds))%>%group_by(Species)%>%filter(Species=="NG")%>%
  as.data.frame(cm_sp_df_res10$preds)

cmng_mod10_sp_meanpreds=as.data.frame(cmng_sp_df_res10$preds[,1])

colnames(cmng_mod10_sp_meanpreds)[1]="mean_preds"

cm_main_dat=cm_10pd_sp%>%select(Species, YR, Julian)%>%filter(!Species=="NG")
cmng_main_dat=cm_10pd_sp%>%select(Species, YR, Julian)%>%filter(Species=="NG")

cm_sp_df_all10=cbind(cm_mod10_sp_meanpreds, cm_main_dat)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))
cm_sp_df_diff10=cm_sp_df_all10%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)
 
cmng_sp_df_all10=cbind(cmng_mod10_sp_meanpreds, cmng_main_dat)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))%>%
  mutate(diet=case_when(Species=="NG" ~ "bird"))

cmng_sp_df_diff10=cmng_sp_df_all10%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
cmng_sp_diff10=left_join(cmng_sp_df_diff10, cm_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

cm_sp_diff10=left_join(cm_sp_df_diff10, cm_sp_ra, by="Species")%>%
  rbind(cmng_sp_diff10)%>%mutate(wgt_shift=shift*rel_abund)%>%
  mutate(diet=case_when(Species=="AK" ~ "insect",
                          Species=="BW" ~ "mammal",
                          Species=="OS" ~ "fish",
                          Species=="ML" ~ "bird",
                          Species=="BE" ~ "fish",
                          Species=="PG" ~ "bird",
                          Species=="SS" ~ "bird",
                          Species=="CH" ~ "mammal",
                          Species=="NH" ~ "mammal",
                          Species=="BV" ~ "carrion",
                          Species=="TV" ~ "carrion",
                          Species=="NG" ~ "bird",
                          Species=="GE" ~ "mammal",
                          Species=="RS" ~ "mammal",
                          Species=="RT" ~ "mammal"))

cm_sp10diff=sum(cm_sp_diff10$wgt_shift) #abundance-weighted shift
cm_p10diff=mean(cm_mod10pred[24:29,]$pred_JD)-mean(cm_mod10pred[1:5,]$pred_JD)

cm_sp10diff/cm_p10diff #proportion of phenological shift by species phenology
1-(cm_sp10diff/cm_p10diff) #proportion of phenological shift by composition
cm_p10diff*(1-(cm_sp10diff/cm_p10diff)) #phenological shift (composition)

##50% PD####
cm_sp_df_mod50=cm_50pd_sp%>%
  select(Species, YR, Julian)%>%group_by(Species)%>%
  nest()%>%
  mutate(model= purrr::map(data,~ brm(bf(Julian~ 1+s(YR)), 
                                      data=., 
                                      iter=2500, family=poisson(), 
                                      prior=priors,
                                      control=list(adapt_delta=0.99))))

cm_sp_df_mod50$preds=lapply(cm_sp_df_mod50$model, predict)
#saveRDS(cm_sp_df_mod50, "CM_sp_mod50.RDS")

cm_sp_df_res50=cm_sp_df_mod50%>%select(-data, -model)%>%
  unnest(cols = c(preds))%>%group_by(Species)%>%filter(!Species=="NG")%>%
  as.data.frame(cm_sp_df_res50$preds)

cm_mod50_sp_meanpreds=as.data.frame(cm_sp_df_res50$preds[,1])

colnames(cm_mod50_sp_meanpreds)[1]="mean_preds"

cmng_sp_df_res50=cm_sp_df_mod50%>%select(-data, -model)%>%
  unnest(cols = c(preds))%>%group_by(Species)%>%filter(Species=="NG")%>%
  as.data.frame(cm_sp_df_res50$preds)

cmng_mod50_sp_meanpreds=as.data.frame(cmng_sp_df_res50$preds[,1])

colnames(cmng_mod50_sp_meanpreds)[1]="mean_preds"

cm_main_dat=cm_10pd_sp%>%select(Species, YR, Julian)%>%filter(!Species=="NG")
cmng_main_dat=cm_10pd_sp%>%select(Species, YR, Julian)%>%filter(Species=="NG")

cm_sp_df_all50=cbind(cm_mod50_sp_meanpreds, cm_main_dat)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))

cm_sp_df_diff50=cm_sp_df_all50%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

cmng_sp_df_all50=cbind(cmng_mod50_sp_meanpreds, cmng_main_dat)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))%>%
  mutate(diet=case_when(Species=="NG" ~ "bird"))

cmng_sp_df_diff50=cmng_sp_df_all50%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
cmng_sp_diff50=left_join(cmng_sp_df_diff50, cm_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

cm_sp_diff50=left_join(cm_sp_df_diff50, cm_sp_ra, by="Species")%>%
  rbind(cmng_sp_diff50)%>%mutate(wgt_shift=shift*rel_abund)%>%
  mutate(diet=case_when(Species=="AK" ~ "insect",
                        Species=="BW" ~ "mammal",
                        Species=="OS" ~ "fish",
                        Species=="ML" ~ "bird",
                        Species=="BE" ~ "fish",
                        Species=="PG" ~ "bird",
                        Species=="SS" ~ "bird",
                        Species=="CH" ~ "mammal",
                        Species=="NH" ~ "mammal",
                        Species=="BV" ~ "carrion",
                        Species=="TV" ~ "carrion",
                        Species=="NG" ~ "bird",
                        Species=="GE" ~ "mammal",
                        Species=="RS" ~ "mammal",
                        Species=="RT" ~ "mammal"))

cm_sp50diff=sum(cm_sp_diff50$wgt_shift) #abundance-weighted shift
cm_p50diff=mean(cm_mod50pred[24:29,]$pred_JD)-mean(cm_mod50pred[1:5,]$pred_JD)

cm_sp50diff/cm_p50diff #proportion of phenological shift by species phenology
1-(cm_sp50diff/cm_p50diff) #proportion of phenological shift by composition
cm_p50diff*(1-(cm_sp50diff/cm_p50diff)) #phenological shift (composition)

##90% PD####
cm_sp_df_mod90=cm_90pd_sp%>%
  select(Species, YR, Julian)%>%group_by(Species)%>%
  nest()%>%
  mutate(model= purrr::map(data,~ brm(bf(Julian~ 1+s(YR)), 
                                      data=., 
                                      iter=2500, family=poisson(), 
                                      prior=priors,
                                      control=list(adapt_delta=0.99))))

cm_sp_df_mod90$preds=lapply(cm_sp_df_mod90$model, predict)
#saveRDS(cm_sp_df_mod90, "CM_sp_mod90.RDS")

cm_sp_df_res90=cm_sp_df_mod90%>%select(-data, -model)%>%
  unnest(cols = c(preds))%>%group_by(Species)%>%filter(!Species=="NG")%>%
  as.data.frame(cm_sp_df_res50$preds)

cm_mod90_sp_meanpreds=as.data.frame(cm_sp_df_res90$preds[,1])

colnames(cm_mod90_sp_meanpreds)[1]="mean_preds"

cmng_sp_df_res90=cm_sp_df_mod90%>%select(-data, -model)%>%
  unnest(cols = c(preds))%>%group_by(Species)%>%filter(Species=="NG")%>%
  as.data.frame(cm_sp_df_res90$preds)

cmng_mod90_sp_meanpreds=as.data.frame(cmng_sp_df_res90$preds[,1])

colnames(cmng_mod90_sp_meanpreds)[1]="mean_preds"

cm_main_dat=cm_10pd_sp%>%select(Species, YR, Julian)%>%filter(!Species=="NG")
cmng_main_dat=cm_10pd_sp%>%select(Species, YR, Julian)%>%filter(Species=="NG")

cm_sp_df_all90=cbind(cm_mod90_sp_meanpreds, cm_main_dat)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))

cm_sp_df_diff90=cm_sp_df_all90%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

cmng_sp_df_all90=cbind(cmng_mod90_sp_meanpreds, cmng_main_dat)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))%>%
  mutate(diet=case_when(Species=="NG" ~ "bird"))

cmng_sp_df_diff90=cmng_sp_df_all90%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
cmng_sp_diff90=left_join(cmng_sp_df_diff90, cm_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

cm_sp_diff90=left_join(cm_sp_df_diff90, cm_sp_ra, by="Species")%>%
  rbind(cmng_sp_diff90)%>%mutate(wgt_shift=shift*rel_abund)%>%
  mutate(diet=case_when(Species=="AK" ~ "insect",
                        Species=="BW" ~ "mammal",
                        Species=="OS" ~ "fish",
                        Species=="ML" ~ "bird",
                        Species=="BE" ~ "fish",
                        Species=="PG" ~ "bird",
                        Species=="SS" ~ "bird",
                        Species=="CH" ~ "mammal",
                        Species=="NH" ~ "mammal",
                        Species=="BV" ~ "carrion",
                        Species=="TV" ~ "carrion",
                        Species=="NG" ~ "bird",
                        Species=="GE" ~ "mammal",
                        Species=="RS" ~ "mammal",
                        Species=="RT" ~ "mammal"))

cm_sp90diff=sum(cm_sp_diff90$wgt_shift) #abundance-weighted shift
cm_p90diff=mean(cm_mod90pred[24:29,]$pred_JD)-mean(cm_mod90pred[1:5,]$pred_JD)

cm_sp90diff/cm_p90diff #proportion of phenological shift by species phenology
1-(cm_sp90diff/cm_p90diff) #proportion of phenological shift by composition
cm_p90diff*(1-(cm_sp90diff/cm_p90diff)) #phenological shift (composition)

#compositional shifts####
cm_sp_tot=cm_sp%>%
  group_by(YR, Species)%>%
  summarise(sp_tot=sum(Count))%>%
  group_by(YR)%>%
  mutate(rel_abund=sp_tot/sum(sp_tot))

cm_sp_abundance=cm_sp_tot%>%
  group_by(Species)%>%filter(!Species=="NG")%>%nest()%>%
  mutate(model= purrr::map(data, ~brm(bf(rel_abund~ 1+s(YR)),
                                      data= .,
                                      iter=2500, 
                                      family=Beta(link = "logit", link_phi = "log"),
                                      prior=priors,
                                      control=list(adapt_delta=0.99))))

cmng_sp_abundance=cm_sp_tot%>%
  group_by(Species)%>%filter(Species=="NG", !(YR=="2011"))%>%nest()%>%
  mutate(model= purrr::map(data, ~brm(bf(rel_abund~ 1+s(YR)),
                                      data= .,
                                      iter=2500, 
                                      family=Beta(link = "logit", link_phi = "log"),
                                      prior=priors,
                                      control=list(adapt_delta=0.99))))

cm_sp_abundance$preds=lapply(cm_sp_abundance$model,predict)
cmng_sp_abundance$preds=lapply(cmng_sp_abundance$model,predict)

cm_allsp_abundance=rbind(cm_sp_abundance, cmng_sp_abundance)

#saveRDS(cm_allsp_abundance,"CM_composition.RDS")

cm_sp_abundance_preds=cm_allsp_abundance%>%select(-data, -model)%>%
  unnest(cols = c(preds))%>%group_by(Species)%>%filter(!Species=="NG")%>%
  as.data.frame(preds)

cmng_sp_abundance_preds=cm_allsp_abundance%>%select(-data, -model)%>%
  unnest(cols = c(preds))%>%group_by(Species)%>%filter(Species=="NG")%>%
  as.data.frame(preds)

YR=rep(1990:2018, each=1, times=14)
YRng=rep(1990:2018, each=1, times=1)
YRng=YRng[!YRng=="2011"]

cm_allsp_abundance_preds=rbind(cm_sp_abundance_preds, cmng_sp_abundance_preds)

cm_sp_abundance_meanpreds=data.frame(cm_sp_abundance_preds$Species,cm_sp_abundance_preds$preds[,1])%>%
  cbind(YR)

colnames(cm_sp_abundance_meanpreds)[2]="mean_preds"
colnames(cm_sp_abundance_meanpreds)[1]="Species"

cmng_sp_abundance_meanpreds=data.frame(cmng_sp_abundance_preds$Species,cmng_sp_abundance_preds$preds[,1])%>%
  cbind(YRng)%>%rename(YR="YRng")

colnames(cmng_sp_abundance_meanpreds)[2]="mean_preds"
colnames(cmng_sp_abundance_meanpreds)[1]="Species"

cm_allsp_abundance_meanpreds=rbind(cm_sp_abundance_meanpreds, cmng_sp_abundance_meanpreds)

main_dat=cm_sp_tot%>%select(Species, YR, rel_abund)
main_datng=cm_sp_tot%>%select(Species, YR, rel_abund)%>%filter(Species=="NG")

cm_sp_abundance_diffs=inner_join(cm_allsp_abundance_meanpreds, main_dat, by=c("Species", "YR"))%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))%>%
  filter(!is.na(period))%>%
  mutate(diet=case_when(Species=="AK" ~ "insect",
                        Species=="BW" ~ "mammal",
                        Species=="OS" ~ "fish",
                        Species=="ML" ~ "bird",
                        Species=="BE" ~ "fish",
                        Species=="PG" ~ "bird",
                        Species=="SS" ~ "bird",
                        Species=="CH" ~ "mammal",
                        Species=="NH" ~ "mammal",
                        Species=="BV" ~ "carrion",
                        Species=="TV" ~ "carrion",
                        Species=="NG" ~ "bird",
                        Species=="GE" ~ "mammal",
                        Species=="RS" ~ "mammal",
                        Species=="RT" ~ "mammal"))

cm_sp_abundance_diffs=cm_sp_abundance_diffs%>%
  select(Species, mean_preds,period, diet)%>%
  filter(!is.na(period))%>%
  group_by(Species, diet,period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#posterior predictions####

##site: 10% PD####
cpreds10=as.data.frame(posterior_predict(cm_mod10))%>%
  mutate(iter=seq(1:5000))%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

quantile(cpreds10$shift, probs=c(0.025,0.975))
mean(cpreds10$shift)

length(which(cpreds10$shift<0))/length(cpreds10$shift)

##site: 50% PD####
cpreds50=as.data.frame(posterior_predict(cm_mod50))%>%
  mutate(iter=seq(1:5000))%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

quantile(cpreds50$shift, probs=c(0.025,0.975))
mean(cpreds50$shift)

length(which(cpreds50$shift<0))/length(cpreds10$shift)

##site: 90% PD####
cpreds90=as.data.frame(posterior_predict(cm_mod90))%>%
  mutate(iter=seq(1:5000))%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

quantile(cpreds90$shift, probs=c(0.025,0.975))
mean(cpreds90$shift)

length(which(cpreds90$shift<0))/length(cpreds90$shift)

##species: 10% PD####
cm_spph=cm_sp_df_mod10%>%filter(!Species=="NG")%>%
  mutate(preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

cmng_spph=cm_sp_df_mod10%>%filter(Species=="NG")%>%
  mutate(preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

cm_splist=cm_spph%>%select(Species, iter)
cmng_splist=cmng_spph%>%select(Species, iter)

cm_spph2=as.data.frame(cm_spph$preds[,])
cmng_spph2=as.data.frame(cmng_spph$preds[,])

cm_spph3=cbind(cm_spph2, cm_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

cmng_spph3=cbind(cmng_spph2, cmng_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:28)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V24", "V25", "V26", "V27", "V28") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

cmall_spph3=rbind(cm_spph3, cmng_spph3)

cm10_spph_quant=cmall_spph3%>%
  group_by(Species)%>%
  reframe(
    q0=quantile(shift, 0.5),
    q1=quantile(shift, 0.025),
    q2=quantile(shift, 0.975))

cm10_spph_prob=cmall_spph3%>%
  group_by(Species)%>%
  reframe(
    prob_shift_advance=length(which(shift<0))/ length(shift),
    prob_shift_delay=length(which(shift>0))/ length(shift))

##species: 50% PD####
cm50_spph=cm_sp_df_mod50%>%filter(!Species=="NG")%>%
  mutate(preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

cmng50_spph=cm_sp_df_mod50%>%filter(Species=="NG")%>%
  mutate(preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

cm_splist=cm_spph%>%select(Species, iter)
cmng_splist=cmng_spph%>%select(Species, iter)

cm50_spph2=as.data.frame(cm50_spph$preds[,])
cmng50_spph2=as.data.frame(cmng50_spph$preds[,])

cm50_spph3=cbind(cm50_spph2, cm_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

cmng50_spph3=cbind(cmng50_spph2, cmng_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:28)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V24", "V25", "V26", "V27", "V28") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

cmall50_spph3=rbind(cm50_spph3, cmng50_spph3)

cm50_spph_quant=cmall50_spph3%>%
  group_by(Species)%>%
  reframe(
    q0=quantile(shift, 0.5),
    q1=quantile(shift, 0.025),
    q2=quantile(shift, 0.975))

cm50_spph_prob=cmall50_spph3%>%
  group_by(Species)%>%
  reframe(
    prob_shift_advance=length(which(shift<0))/ length(shift),
    prob_shift_delay=length(which(shift>0))/ length(shift))

##species: 90% PD####
cm90_spph=cm_sp_df_mod90%>%filter(!Species=="NG")%>%
  mutate(preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

cmng90_spph=cm_sp_df_mod90%>%filter(Species=="NG")%>%
  mutate(preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

cm_splist=cm_spph%>%select(Species, iter)
cmng_splist=cmng_spph%>%select(Species, iter)

cm90_spph2=as.data.frame(cm90_spph$preds[,])
cmng90_spph2=as.data.frame(cmng90_spph$preds[,])

cm90_spph3=cbind(cm90_spph2, cm_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

cmng90_spph3=cbind(cmng90_spph2, cmng_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:28)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V24", "V25", "V26", "V27", "V28") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

cmall90_spph3=rbind(cm90_spph3, cmng90_spph3)

cm90_spph_quant=cmall90_spph3%>%
  group_by(Species)%>%
  reframe(
    q0=quantile(shift, 0.5),
    q1=quantile(shift, 0.025),
    q2=quantile(shift, 0.975))

cm90_spph_prob=cmall90_spph3%>%
  group_by(Species)%>%
  reframe(
    prob_shift_advance=length(which(shift<0))/ length(shift),
    prob_shift_delay=length(which(shift>0))/ length(shift))

##composition####
cm_abund=cm_sp_abundance%>%filter(!Species=="NG")%>%
  mutate(preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()%>%select(-data, -model)%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

cmng_abund=cmng_sp_abundance%>%
  mutate(preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()%>%select(-data, -model)%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

cm_abund2=as.data.frame(cm_abund$preds[,])
cmng_abund2=as.data.frame(cmng_abund$preds[,])

cm_abund3=cbind(cm_abund2, cm_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

cmng_abund3=cbind(cmng_abund2, cmng_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:28)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V24","V25", "V26", "V27", "V28") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

cmall_abund3=rbind(cm_abund3, cmng_abund3)

cm_abund_quant=cmall_abund3%>%
  group_by(Species)%>%
  reframe(
    q0=quantile(shift, 0.5),
    q1=quantile(shift, 0.025),
    q2=quantile(shift, 0.975))

cm_abund_prob=cmall_abund3%>%
  group_by(Species)%>%
  reframe(
    prob_shift_increase=length(which(shift>0))/ length(shift),
    prob_shift_decrease=length(which(shift<0))/ length(shift))
