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
         -UNID.ACCIPITER,-UNID.BUTEO, -UNID.EAGLE, -UNID.FALCON, -UNID.RAPTOR, -TUVU, -BLVU)%>%
  mutate(yr_tot=rowSums(.[6:19], na.rm=T), date=make_date(YR, MO, DAY))%>%filter(!(YR<1956), !(YR>1990))

HMS$Julian=as.POSIXlt(HMS$date)$yday
HMS[, 6:19][is.na(HMS[, 6:19])] <- 0

#select only 1956-1990 (Aug 15-Dec 15: 226-348) and without vultures

HMS=HMS%>%filter(!Julian<226, !Julian>348, !YR<1956, !YR>1990)

#species-level counts
hms_sp=HMS%>%
  pivot_longer(cols=6:19, names_to = "Species", values_to="Count")

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
hms_sp_min=hms_sp_tot%>%group_by(Species)%>%filter(sp_tot>0)%>%
  summarise(n_yr=length(YR))

#daily totals per species
hms_sp_day=hms_sp%>%
  group_by(Species,YR, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(YR)

#site-level####

##10%PD####

#10% passage date: species
hms_10pd_sp=hms_sp_day%>%group_by(Species, YR)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                                                           cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                                                           PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, YR)%>%filter(PD_10==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%

#10% passage date: site
hms_10pd_site=HMS%>%group_by(YR)%>%mutate(cum_tot=cumsum(yr_tot), #get cumulative sum for each year
                                                             cum_pt=cum_tot/sum(yr_tot),#calculate % of daily count/annual count
                                                             PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(YR)%>%filter(PD_10==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(YR, Julian, cum_pt)

#site-level
priors = c(prior(normal(0, 10), class = Intercept),
           prior(normal(0, 10), class = "sds", coef="s(YR)"),
           prior(cauchy(0, 5), class=b, coef="sYR_1"))

hm1956_mod10=brm(bf(Julian~ 1 + s(YR)), 
             data=hms_10pd_site, iter=2500, family=poisson(link="log"), 
             prior = priors,
             control=list(adapt_delta=0.99))

mod10pred=data.frame(predict(hm1956_mod10))
colnames(mod10pred)[1]="pred_JD"

m1pred=cbind(hms_10pd_site, mod10pred)

mean(mod10pred[31:35,]$pred_JD)-mean(mod10pred[1:5,]$pred_JD) 
mean(hms_10pd_site[31:35,]$Julian)-mean(hms_10pd_site[1:5,]$Julian) 

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

hms1956_mod50=brm(bf(Julian~ 1+s(YR)), 
          data=hms_50pd_site, iter=2500, prior=priors,
          family=poisson(), control=list(adapt_delta=0.99))

mod50pred=data.frame(predict(hms1956_mod50))
colnames(mod50pred)[1]="pred_JD"

m5pred=cbind(hms_50pd_site, mod50pred)

mean(mod50pred[31:35,]$pred_JD)-mean(mod50pred[1:5,]$pred_JD) #-7.9
mean(hms_50pd_site[31:35,]$Julian)-mean(hms_50pd_site[1:5,]$Julian) #-7.9

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

#site-level
hm1956_mod90=brm(bf(Julian~ 1+s(YR)), 
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
         preds=purrr::map(model, posterior_predict))

sp_df_res10=sp_df_mod10%>%select(-data, -model)
sp_df_res10preds=as.data.frame(sp_df_res10$preds[,1])
colnames(sp_df_res10preds)[1]="mean_preds"

main_dat10=hms_10pd_sp%>%select(Species, YR, Julian)

sp_df_all=cbind(sp_df_res10preds, main_dat10)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1987) ~ "T1",
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
         preds=purrr::map(model, posterior_predict))

sp_df_res50=sp_df_mod50%>%select(-data, -model)
sp_df_res50preds=as.data.frame(sp_df_res50$preds[,1])
colnames(sp_df_res50preds)[1]="mean_preds"

main_dat50=hms_50pd_sp%>%select(Species, YR, Julian)

sp_df50_all=cbind(sp_df_res50preds, main_dat50)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1987) ~ "T1",
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
         preds=purrr::map(model, posterior_predict))

sp_df_res90=sp_df_mod90%>%select(-data, -model)
sp_df_res90preds=as.data.frame(sp_df_res90$preds[,1])
colnames(sp_df_res90preds)[1]="mean_preds"

main_dat90=hms_90pd_sp%>%select(Species, YR, Julian)

sp_df90_all=cbind(sp_df_res90preds, main_dat90)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1987) ~ "T1",
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

sp90diff/p90diff #proportion of phenological shift by species phenology
1-(sp90diff/p90diff) #proportion of phenological shift by composition
p90diff*(1-(sp90diff/p90diff)) #phenological shift (composition)

#visualize data####

##shifts####
##10% PD####
plot(hms_10pd_site$Julian~hms_10pd_site$YR, type="l", main="Hawk Mountain", ylab="10% passage date", xlab="year")

lines(m1pred$pred_JD~m1pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1990, x1=1987,y0=mean(hms_10pd_site[1:5,]$Julian), y1=mean(hms_10pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2014, x1=2018,y0=mean(hms_10pd_site[32:36,]$Julian), y1=mean(hms_10pd_site[32:36,]$Julian), lwd=2, col="orange")

segments(x0=1990, x1=1987,y0=mean(m1pred[1:5,]$pred_JD), y1=mean(m1pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2014, x1=2018,y0=mean(m1pred[32:36,]$pred_JD), y1=mean(m1pred[32:36,]$pred_JD), lwd=2, col="red")

##50% PD####
plot(hms_50pd_site$Julian~hms_50pd_site$YR, type="l", ylab="50% passage date", xlab="year")

lines(m5pred$pred_JD~m5pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1990, x1=1987,y0=mean(hms_50pd_site[1:5,]$Julian), y1=mean(hms_50pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2014, x1=2018,y0=mean(hms_50pd_site[32:36,]$Julian), y1=mean(hms_50pd_site[32:36,]$Julian), lwd=2, col="orange")

segments(x0=1990, x1=1987,y0=mean(m5pred[1:5,]$pred_JD), y1=mean(m5pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2014, x1=2018,y0=mean(m5pred[32:36,]$pred_JD), y1=mean(m5pred[32:36,]$pred_JD), lwd=2, col="red")

##90% PD####
plot(hms_90pd_site$Julian~hms_90pd_site$YR, type="l", ylab="90% passage date", xlab="year")

lines(m9pred$pred_JD~m9pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1990, x1=1987,y0=mean(hms_90pd_site[1:5,]$Julian), y1=mean(hms_90pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2014, x1=2018,y0=mean(hms_90pd_site[32:36,]$Julian), y1=mean(hms_90pd_site[32:36,]$Julian), lwd=2, col="orange")

segments(x0=1990, x1=1987,y0=mean(m9pred[1:5,]$pred_JD), y1=mean(m9pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2014, x1=2018,y0=mean(m9pred[32:36,]$pred_JD), y1=mean(m9pred[32:36,]$pred_JD), lwd=2, col="red")


#contributions####
datphen=matrix(c(10,50,90,
                 -0.80,-6.68,-5.21, 
                 0.12,3.23,2.41,
                 -0.92, -9.91, -7.62), nrow=3, ncol=4, byrow=FALSE, 
               dimnames= list(c("1", "2", "3"),
                              c("metric","overall", "species phenology", "composition")))
dat_phen=as.data.frame(datphen)%>%
  pivot_longer(cols=2:4, names_to = "shift", values_to="value" )

ggplot(dat_phen)+geom_col(aes(x=shift, y=value))+ggtitle("Hawk Mountain")+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_x_discrete(limits=c("overall", "species phenology", "composition"))+
  facet_wrap(~metric)+theme_classic()+abline(h=0)+ylab("shift (days)")

#compositional shifts####
hms_sp_tot=hms_sp%>%
  group_by(YR, Species)%>%
  summarise(sp_tot=sum(Count))%>%filter(!YR<1990,!Species=="RLHA")%>%
  group_by(YR)%>%
  mutate(rel_abund=sp_tot/sum(sp_tot))

sp_abundance3=hms_sp_tot%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(rel_abund~ 1+s(YR)), 
                         data=data, 
                         iter=2500, 
                         family=zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit"),
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, posterior_predict))%>%
  unnest(preds)%>%as.data.frame()

spab=sp_abundance3%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

hms_splist=spab%>%select(Species, iter)

spab2=as.data.frame(spab$preds[,])

spab3=cbind(spab2, hms_splist)%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:29)%>%
  mutate(period=case_when(year%in%c("V1", "V2", "V3", "V4", "V5") ~ "T1",
                          year%in%c("V25", "V26", "V27", "V28", "V29") ~ "T2"))%>%
  filter(!is.na(period))%>%select(-year)%>%group_by(Species,iter,period)%>%
  summarise(mean_est=mean(estimates))%>%
  pivot_wider(names_from=period, values_from = mean_est, names_sep = ".")%>%
  mutate(shift=T2-T1)

sp_quant=spab3%>%
  group_by(Species)%>%
  reframe(
    q50=quantile(shift, 0.5),
    qLB=quantile(shift, 0.025),
    qUB=quantile(shift, 0.975))

spab_prob=spab3%>%
  group_by(Species)%>%
  reframe(
    prob_shift_decrease=length(which(shift<0))/ length(shift),
    prob_shift_increase=length(which(shift>0))/ length(shift))

sp_relab=sp_abundance3%>%select(-data, -model)
sp_rela_preds=as.data.frame(sp_relab$preds[,1])
colnames(sp_rela_preds)[1]="mean_rela"

sp_df_rela=cbind(sp_rela_preds, hms_sp_tot)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))

sp_df_rela_diff=sp_df_rela%>%
  select(Species, mean_rela,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_rela=mean(mean_rela))%>%
  pivot_wider(names_from = period, values_from=mean_rela, names_sep = ".")%>%
  mutate(shift= T2-T1)

#determine passage order of species (based on 50% passage date)
hms_ord=hms_50pd_sp%>%group_by(Species)%>%
  summarise(min_time=min(Julian))%>%arrange(min_time)

ggplot(sp_df_rela_diff)+geom_col(aes(x=Species, y=shift))+theme_classic()+
  geom_hline(yintercept = 0, linetype = 2)+ggtitle("Hawk Mountain")+
  ylab("relative abundance shifts")+
  scale_x_discrete(limits=c("RLHA", "BAEA", "BWHA", "AMKE", "OSPR", "PEFA","SSHA", "MERL",
                            "NOHA", "COHA", "BLVU", "TUVU", "RSHA", "RTHA", "GOEA", "NOGO"))

#posterior predictions####
##site:10% PD####
hpreds10=as.data.frame(posterior_predict(mod10))%>%
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
hpreds50=as.data.frame(posterior_predict(mod50))%>%
  mutate(iter=seq(1:5000))%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:36)%>%
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
hpreds90=as.data.frame(posterior_predict(mod90))%>%
  mutate(iter=seq(1:5000))%>%
  pivot_longer(names_to="year", values_to="estimates", cols=1:36)%>%
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
spph=sp_df_mod10%>%unnest(preds)%>%as.data.frame()%>%
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
spph50=sp_df_mod50%>%unnest(preds)%>%as.data.frame()%>%
  group_by(Species)%>%
  mutate(iter=seq(1:5000))

hms_splist=spab%>%select(Species, iter)

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