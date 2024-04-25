require(dplyr)
require(tidyr)
require(lubridate)
require(ggplot2)
require(broom)
require(brms)

#data manipulation####
wach=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/WachusettMountain.csv")

WachMt=wach%>%
  select(-Duration, -Observer, -UA, -UB, -UF, -UE, -X, -BV, -RL, -MK, -SW, -GE, -NG)%>%
  mutate(yr_tot=rowSums(.[2:13], na.rm=T),
         date=as.Date(Date, format="%m/%d/%Y"),
         YR=year(date), MO=month(date), DAY=day(date))

WachMt$Julian=as.POSIXlt(WachMt$date)$yday
WachMt[, 2:13][is.na(WachMt[, 2:13])] <- 0

WachMt=WachMt%>%filter(!Julian<235, !Julian>320)

#species-level counts
wa_sp=WachMt%>%select(-yr_tot)%>%
  pivot_longer(cols=2:13, names_to = "Species", values_to="Count")

#annual totals per species
wa_sp_tot=wa_sp%>%
  group_by(YR, Species)%>%
  summarise(sp_tot=sum(Count))

#species relative abundances
wa_sp_ra=wa_sp_tot%>%
  group_by(Species)%>%
  summarise(ave_ct=mean(sp_tot))%>%
  mutate(rel_abund=ave_ct/sum(ave_ct))%>%
  filter(!Species%in%c("BV", "MK","RL", "SW", "GE", "NG"))

#determine minimum year of non-0 detects for all species
#use to determine which species to exclude
was_sp_min=wa_sp_tot%>%group_by(Species)%>%filter(sp_tot>0)%>%summarise(min_yr=min(YR))

#daily totals per species
wa_sp_day=wa_sp%>%
  group_by(Species,YR, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(YR)

#site-level####

##10%PD####

#10% passage date: species
wa_10pd_sp=wa_sp_day%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                             cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                             PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_10==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#10% passage date: site
wa_10pd_site=WachMt%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                       cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                      PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_10==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

#site-level
priors = c(prior(normal(0, 10), class = Intercept),
           prior(normal(0, 10), class = "sds", coef="s(YR)"),
           prior(cauchy(0, 5), class=b, coef="sYR_1"))

wa_mod10=brm(bf(Julian~ 1 + s(YR)), 
          data=wa_10pd_site, iter=2500, family=poisson(link="log"), 
          prior = priors,
          control=list(adapt_delta=0.99))

wa_mod10pred=data.frame(predict(wa_mod10))
colnames(wa_mod10pred)[1]="pred_JD"

wa_m1pred=cbind(wa_10pd_site, wa_mod10pred)

mean(wa_mod10pred[17:21,]$pred_JD)-mean(wa_mod10pred[1:5,]$pred_JD) 
mean(wa_10pd_site[17:21,]$Julian)-mean(wa_10pd_site[1:5,]$Julian) 

##50% PD####
#50% passage date: species
wa_50pd_sp=wa_sp_day%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                            cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                           PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_50==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#50% passage date: site
wa_50pd_site=WachMt%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                       cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                       PD_50=if_else(cum_pt>=0.5,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_50==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

wa_mod50=brm(bf(Julian~ 1+s(YR)), 
          data=wa_50pd_site, iter=2500, prior=priors,
          family=poisson(), control=list(adapt_delta=0.99))

wa_mod50pred=data.frame(predict(wa_mod50))
colnames(wa_mod50pred)[1]="pred_JD"

wa_m5pred=cbind(wa_50pd_site, wa_mod50pred)

mean(wa_mod50pred[17:21,]$pred_JD)-mean(wa_mod50pred[1:5,]$pred_JD) #-7.9
mean(wa_50pd_site[17:21,]$Julian)-mean(wa_50pd_site[1:5,]$Julian) #-7.9

##90% PD####
#90% passage date: species
wa_90pd_sp=wa_sp_day%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                            cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                            PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_90==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

##90% passage date: site
wa_90pd_site=WachMt%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                                             cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                                             PD_90=if_else(cum_pt>=0.9,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_90==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

#site-level
wa_mod90=brm(bf(Julian~ 1+s(YR)), 
          data=wa_90pd_site, iter=2500, 
          family=poisson(), 
          prior=priors,
          control=list(adapt_delta=0.99))

wa_mod90pred=data.frame(predict(wa_mod90))
colnames(wa_mod90pred)[1]="pred_JD"

wa_m9pred=cbind(wa_90pd_site, wa_mod90pred)

mean(wa_mod90pred[17:21,]$pred_JD)-mean(wa_mod90pred[1:5,]$pred_JD) #-7.9
mean(wa_90pd_site[17:21,]$Julian)-mean(wa_90pd_site[1:5,]$Julian) #-7.9

#phenological shifts####
wpreds10=as.data.frame(posterior_predict(wa_mod10))%>%
  mutate(iter=seq(1:5000))%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:21)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V17", "V18", "V19", "V20", "V21") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

quantile(wpreds10$shift, probs=c(0.025,0.975))
mean(wpreds10$shift)
length(which(wpreds10$shift>0))/length(wpreds10$shift)

wpreds50=as.data.frame(posterior_predict(wa_mod50))%>%
  mutate(iter=seq(1:5000))%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:21)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V17", "V18", "V19", "V20", "V21") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

quantile(wpreds50$shift, probs=c(0.025,0.975))
mean(wpreds50$shift)
length(which(wpreds50$shift>0))/length(wpreds50$shift)

wpreds90=as.data.frame(posterior_predict(wa_mod90))%>%
  mutate(iter=seq(1:5000))%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:21)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V17", "V18", "V19", "V20", "V21") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

quantile(wpreds90$shift, probs=c(0.025,0.975))
mean(wpreds90$shift)
length(which(wpreds90$shift>0))/length(wpreds90$shift)

#species-level####

##10% PD####
wa_sp_df_mod10=wa_10pd_sp%>%
  select(Species, YR, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(YR)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()

wa_sp_df_res10=wa_sp_df_mod10%>%select(-data, -model)
wa_sp_df_res10preds=as.data.frame(wa_sp_df_res10$preds[,1])
colnames(wa_sp_df_res10preds)[1]="mean_preds"

wa_main_dat10=wa_10pd_sp%>%select(Species, YR, Julian)

wa_sp_df_all=cbind(wa_sp_df_res10preds, wa_main_dat10)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(2002:2006) ~ "T1",
                          YR%in%c(2018:2022) ~ "T2"))
wa_sp_df_diff=wa_sp_df_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
wa_sp_diff10=left_join(wa_sp_df_diff, wa_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

wa_sp10diff=sum(wa_sp_diff10$wgt_shift) #abundance-weighted shift
wa_p10diff=mean(wa_mod10pred[17:21,]$pred_JD)-mean(wa_mod10pred[1:5,]$pred_JD)

wa_sp10diff/wa_p10diff #proportion of phenological shift by species phenology
1-(wa_sp10diff/wa_p10diff) #proportion of phenological shift by composition
wa_p10diff*(1-(wa_sp10diff/wa_p10diff)) #phenological shift (composition)

##50% PD####
wa_sp_df_mod50=wa_50pd_sp%>%
  select(Species, YR, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(YR)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

wa_sp_df_res50=wa_sp_df_mod50%>%select(-data, -model)
wa_sp_df_res50preds=as.data.frame(wa_sp_df_res50$preds[,1])
colnames(wa_sp_df_res50preds)[1]="mean_preds"

wa_main_dat50=wa_50pd_sp%>%select(Species, YR, Julian)

wa_sp_df50_all=cbind(wa_sp_df_res50preds, wa_main_dat50)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(2002:2006) ~ "T1",
                          YR%in%c(2018:2022) ~ "T2"))
wa_sp_df50_diff=wa_sp_df50_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
wa_sp_diff50=left_join(wa_sp_df50_diff, wa_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

wa_sp50diff=sum(wa_sp_diff50$wgt_shift)
wa_p50diff=mean(wa_mod50pred[17:21,]$pred_JD)-mean(wa_mod50pred[1:5,]$pred_JD)

wa_sp50diff/wa_p50diff #abundance-weighted shift
1-(wa_sp50diff/wa_p50diff) #proportion of phenological shift by composition
wa_p50diff*(1-(wa_sp50diff/wa_p50diff)) #phenological shift by composition

##90% PD####
wa_sp_df_mod90=wa_90pd_sp%>%
  select(Species, YR, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(YR)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

wa_sp_df_res90=wa_sp_df_mod90%>%select(-data, -model)
wa_sp_df_res90preds=as.data.frame(wa_sp_df_res90$preds[,1])
colnames(wa_sp_df_res90preds)[1]="mean_preds"

wa_main_dat90=wa_90pd_sp%>%select(Species, YR, Julian)

wa_sp_df90_all=cbind(wa_sp_df_res90preds, wa_main_dat90)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(2002:2006) ~ "T1",
                          YR%in%c(2018:2022) ~ "T2"))
wa_sp_df90_diff=wa_sp_df90_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
wa_sp_diff90=left_join(wa_sp_df90_diff, wa_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

wa_sp90diff=sum(wa_sp_diff90$wgt_shift)
wa_p90diff=mean(wa_mod90pred[17:21,]$pred_JD)-mean(wa_mod90pred[1:5,]$pred_JD)

wa_sp90diff/wa_p90diff #proportion of phenological shift by species phenology
1-(wa_sp90diff/wa_p90diff) #proportion of phenological shift by composition
wa_p90diff*(1-(wa_sp90diff/wa_p90diff)) #phenological shift (composition)

#visualize data####

##shifts####
##10% PD####
plot(wa_10pd_site$Julian~wa_10pd_site$YR, type="l", main="Wachusett Mountain", ylab="10% passage date", xlab="year")

lines(wa_m1pred$pred_JD~wa_m1pred$YR, lty=2, col="grey", lwd=2)

segments(x0=2002, x1=2006,y0=mean(wa_10pd_site[1:5,]$Julian), y1=mean(wa_10pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2018, x1=2022,y0=mean(wa_10pd_site[17:21,]$Julian), y1=mean(wa_10pd_site[17:21,]$Julian), lwd=2, col="orange")

segments(x0=2002, x1=2006,y0=mean(wa_mod10pred[1:5,]$pred_JD), y1=mean(wa_mod10pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2018, x1=2022,y0=mean(wa_mod10pred[17:21,]$pred_JD), y1=mean(wa_mod10pred[17:21,]$pred_JD), lwd=2, col="red")

##50% PD####
plot(wa_50pd_site$Julian~wa_50pd_site$YR, type="l", main="Wachusett Mountain", ylab="50% passage date", xlab="year")

lines(wa_m5pred$pred_JD~wa_m5pred$YR, lty=2, col="grey", lwd=2)

segments(x0=2002, x1=2006,y0=mean(wa_50pd_site[1:5,]$Julian), y1=mean(wa_50pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2018, x1=2022,y0=mean(wa_50pd_site[17:21,]$Julian), y1=mean(wa_50pd_site[17:21,]$Julian), lwd=2, col="orange")

segments(x0=2002, x1=2006,y0=mean(wa_mod50pred[1:5,]$pred_JD), y1=mean(wa_mod50pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2018, x1=2022,y0=mean(wa_mod50pred[17:21,]$pred_JD), y1=mean(wa_mod50pred[17:21,]$pred_JD), lwd=2, col="red")

##90% PD####
plot(wa_90pd_site$Julian~wa_90pd_site$YR, type="l", main="Wachusett Mountain", ylab="90% passage date", xlab="year")

lines(wa_m9pred$pred_JD~wa_m9pred$YR, lty=2, col="grey", lwd=2)

segments(x0=2002, x1=2006,y0=mean(wa_90pd_site[1:5,]$Julian), y1=mean(wa_90pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2018, x1=2022,y0=mean(wa_90pd_site[17:21,]$Julian), y1=mean(wa_90pd_site[17:21,]$Julian), lwd=2, col="orange")

segments(x0=2002, x1=2006,y0=mean(wa_mod90pred[1:5,]$pred_JD), y1=mean(wa_mod90pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2018, x1=2022,y0=mean(wa_mod90pred[17:21,]$pred_JD), y1=mean(wa_mod90pred[17:21,]$pred_JD), lwd=2, col="red")


#contributions####
wa_datphen=matrix(c(10,50,90,
                 0.62,-1.15,14.37, 
                 1.13,-0.34,5.05,
                 -0.82, -0.77, 9.32), nrow=3, ncol=4, byrow=FALSE, 
               dimnames= list(c("1", "2", "3"),
                              c("metric","overall", "species phenology", "composition")))
wa_dat_phen=as.data.frame(wa_datphen)%>%
  pivot_longer(cols=2:4, names_to = "shift", values_to="value" )

ggplot(wa_dat_phen)+geom_col(aes(x=shift, y=value))+ggtitle("Wachusett Mountain")+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_x_discrete(limits=c("overall", "species phenology", "composition"))+
  facet_wrap(~metric)+theme_classic()+abline(h=0)+ylab("shift (days)")

#compositional shifts####
wa_sp_tot=wa_sp%>%
  group_by(YR, Species)%>%
  summarise(sp_tot=sum(Count))%>%
  group_by(YR)%>%
  mutate(rel_abund=sp_tot/sum(sp_tot))

wa_sp_abundance2=wa_sp_tot%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(rel_abund~ 1+s(YR)), 
                         data=data, 
                         iter=2500, 
                         family=zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit"),
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

wa_sp_relab=wa_sp_abundance2%>%select(-data, -model)
wa_sp_rela_preds=as.data.frame(wa_sp_relab$preds[,1])
colnames(wa_sp_rela_preds)[1]="mean_rela"

wa_sp_df_rela=cbind(wa_sp_rela_preds, wa_sp_tot)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(2002:2006) ~ "T1",
                          YR%in%c(2018:2022) ~ "T2"))

wa_sp_df_rela_diff=wa_sp_df_rela%>%
  select(Species, mean_rela,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_rela=mean(mean_rela))%>%
  pivot_wider(names_from = period, values_from=mean_rela, names_sep = ".")%>%
  mutate(shift= T2-T1)

#determine passage order of species (based on 50% passage date)
wa_ord=wa_50pd_sp%>%group_by(Species)%>%
  summarise(min_time=min(Julian))%>%arrange(min_time)

ggplot(wa_sp_df_rela_diff)+geom_col(aes(x=Species, y=shift))+theme_classic()+
  geom_hline(yintercept = 0, linetype = 2)+ggtitle("Wachusett Mountain")+
  ylab("relative abundance shifts")+
  scale_x_discrete(limits=c("BE","RS","RT", "BW", "CH", "ML",
                            "PG", "NH", "OS", "SS", "AK", "TV"),
                   labels=c("BAEA", "RSHA", "RTHA", "BWHA", "COHA","MERL",
                            "PEFA", "NOHA", "OSPR", "SSHA", "AMKE", "TUVU"))

wa_sp_diff_rel=wa_sp_df_rela_diff%>%
  select(Species, shift)%>%
  rename(rela_shift=shift)

wa_sp_diff10_time=wa_sp_diff10%>%
  mutate(metric="10")%>%
  rename(time_shift=shift)%>%
  select(Species, time_shift, wgt_shift, metric)

wa_sp_diff50_time=wa_sp_diff50%>%
  mutate(metric="50")%>%
  rename(time_shift=shift)%>%
  select(Species, time_shift, wgt_shift, metric)

wa_sp_diff90_time=wa_sp_diff90%>%
  mutate(metric="90")%>%
  rename(time_shift=shift)%>%
  select(Species, time_shift, wgt_shift, metric)

wa_diff_all=do.call("rbind", list(wa_sp_diff10_time, wa_sp_diff50_time, wa_sp_diff90_time))%>%
  left_join(wa_sp_diff_rel, by="Species")

wa_tot_ave=wa_sp_tot%>%
  group_by(Species)%>%
  summarise(ave_tot=mean(rel_abund))

wa_diffs_all=left_join(wa_diff_all, wa_tot_ave, by="Species")%>%
  mutate(species=case_when(Species=="AK" ~ "AMKE",
                           Species=="BE" ~ "BAEA",
                           Species=="BW" ~ "BWHA",
                           Species=="CH" ~ "COHA",
                           Species=="ML" ~ "MERL",
                           Species=="NH" ~ "NOHA",
                           Species=="OS" ~ "OSPR",
                           Species=="PG" ~ "PEFA",
                           Species=="RS" ~ "RSHA",
                           Species=="RT" ~ "RTHA",
                           Species=="SS" ~ "SSHA",
                           Species=="TV" ~ "TUVU"))


waplot=ggplot(wa_diffs_all)+
  geom_point(aes(x=time_shift, y=rela_shift, col=species, size=ave_tot))+
  theme_classic()+geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)+ggtitle("Wachusett Mountain")+
  ylab("relative abundance shifts")+facet_wrap(~metric)+xlab("phenological shift")+
  scale_color_viridis_d()

waplot+guides(size=FALSE)

#POSTERIOR COMPOSITION####
wa_sp_abundance2=wa_sp_tot%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(rel_abund~ 1+s(YR)), 
                         data=data, 
                         iter=2500, 
                         family=zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit"),
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()


wa_spab=wa_sp_abundance2%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

wa_splist=wa_spab%>%select(Species, iter)

wa_spab2=as.data.frame(wa_spab$preds[,])

wa_spab3=cbind(wa_spab2, wa_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:21)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V17", "V18", "V19", "V20", "V21") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

wa_sp_quant=wa_spab3%>%
  group_by(Species)%>%
  reframe(
    q0=quantile(shift, 0.5),
    q1=quantile(shift, 0.025),
    q2=quantile(shift, 0.975))
