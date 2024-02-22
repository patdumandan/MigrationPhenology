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
  filter(!YR<1995)

CapeMay$Julian=as.POSIXlt(CapeMay$date)$yday
CapeMay[, 2:16][is.na(CapeMay[, 2:16])] <- 0

#species-level counts
cm_sp=CapeMay%>%select(-yr_tot)%>%
  pivot_longer(cols=2:16, names_to = "Species", values_to="Count")

#annual totals per species
cm_sp_tot=cm_sp%>%
  group_by(YR, Species)%>%
  summarise(sp_tot=sum(Count))

ggplot(cm_sp_tot)+geom_line(aes(x=YR, y=log(sp_tot)))+facet_wrap(~Species)+theme_classic()

#species relative abundances
cm_sp_ra=cm_sp_tot%>%
  group_by(Species)%>%
  summarise(ave_ct=mean(sp_tot))%>%
  mutate(rel_abund=ave_ct/sum(ave_ct))

#determine minimum year of non-0 detects for all species
#use to determine which species to exclude
cm_sp_min=cm_sp_tot%>%group_by(Species)%>%filter(sp_tot>0)%>%summarise(min_yr=min(YR))

#daily totals per species
cm_sp_day=cm_sp%>%
  group_by(Species,YR, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(YR)

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

cm_mod10pred=data.frame(predict(cm_mod10))
colnames(cm_mod10pred)[1]="pred_JD"

cm_m1pred=cbind(cm_10pd_site, cm_mod10pred)

mean(cm_mod10pred[14:28,]$pred_JD)-mean(cm_mod10pred[1:5,]$pred_JD) 
mean(cm_10pd_site[14:28,]$Julian)-mean(cm_10pd_site[1:5,]$Julian) 

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

cm_mod50pred=data.frame(predict(cm_mod50))
colnames(cm_mod50pred)[1]="pred_JD"

cm_m5pred=cbind(cm_50pd_site, cm_mod50pred)

mean(cm_mod50pred[24:28,]$pred_JD)-mean(cm_mod50pred[1:5,]$pred_JD) #-7.9
mean(cm_50pd_site[24:28,]$Julian)-mean(cm_50pd_site[1:5,]$Julian) #-7.9

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

cm_mod90pred=data.frame(predict(cm_mod90))
colnames(cm_mod90pred)[1]="pred_JD"

cm_m9pred=cbind(cm_90pd_site, cm_mod90pred)

mean(cm_mod90pred[24:28,]$pred_JD)-mean(cm_mod90pred[1:5,]$pred_JD) #-7.9
mean(cm_90pd_site[24:28,]$Julian)-mean(cm_90pd_site[1:5,]$Julian) #-7.9

#species-level####

##10% PD####
cm_sp_df_mod10=cm_10pd_sp%>%
  select(Species, YR, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(YR)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

cm_sp_df_res10=cm_sp_df_mod10%>%select(-data, -model)
cm_sp_df_res10preds=as.data.frame(cm_sp_df_res10$preds[,1])
colnames(cm_sp_df_res10preds)[1]="mean_preds"

cm_main_dat10=cm_10pd_sp%>%select(Species, YR, Julian)

cm_sp_df_all=cbind(cm_sp_df_res10preds, cm_main_dat10)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1995:1999) ~ "T1",
                          YR%in%c(2018:2022) ~ "T2"))
cm_sp_df_diff=cm_sp_df_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
cm_sp_diff10=left_join(cm_sp_df_diff, cm_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

cm_sp10diff=sum(cm_sp_diff10$wgt_shift) #abundance-weighted shift
cm_p10diff=mean(cm_mod10pred[24:28,]$pred_JD)-mean(cm_mod10pred[1:5,]$pred_JD)

cm_sp10diff/cm_p10diff #proportion of phenological shift by species phenology
1-(cm_sp10diff/cm_p10diff) #proportion of phenological shift by composition
cm_p10diff*(1-(cm_sp10diff/cm_p10diff)) #phenological shift (composition)

##50% PD####
cm_sp_df_mod50=cm_50pd_sp%>%
  select(Species, YR, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(YR)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

cm_sp_df_res50=cm_sp_df_mod50%>%select(-data, -model)
cm_sp_df_res50preds=as.data.frame(cm_sp_df_res50$preds[,1])
colnames(cm_sp_df_res50preds)[1]="mean_preds"

cm_main_dat50=cm_50pd_sp%>%select(Species, YR, Julian)

cm_sp_df_all=cbind(cm_sp_df_res50preds, cm_main_dat50)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1995:1999) ~ "T1",
                          YR%in%c(2018:2022) ~ "T2"))
cm_sp_df_diff=cm_sp_df_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
cm_sp_diff50=left_join(cm_sp_df_diff, cm_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

cm_sp50diff=sum(cm_sp_diff50$wgt_shift) #abundance-weighted shift
cm_p50diff=mean(cm_mod50pred[24:28,]$pred_JD)-mean(cm_mod50pred[1:5,]$pred_JD)

cm_sp50diff/cm_p50diff #proportion of phenological shift by species phenology
1-(cm_sp50diff/cm_p50diff) #proportion of phenological shift by composition
cm_p50diff*(1-(cm_sp50diff/cm_p50diff)) #phenological shift (composition)

##90% PD####
cm_sp_df_mod90=cm_90pd_sp%>%
  select(Species, YR, Julian)%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(Julian~ 1+s(YR)), 
                         data=data, 
                         iter=2500, family=poisson(), 
                         prior=priors,
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

cm_sp_df_res90=cm_sp_df_mod90%>%select(-data, -model)
cm_sp_df_res90preds=as.data.frame(cm_sp_df_res90$preds[,1])
colnames(cm_sp_df_res90preds)[1]="mean_preds"

cm_main_dat90=cm_90pd_sp%>%select(Species, YR, Julian)

cm_sp_df_all=cbind(cm_sp_df_res90preds, cm_main_dat90)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1995:1999) ~ "T1",
                          YR%in%c(2018:2022) ~ "T2"))
cm_sp_df_diff=cm_sp_df_all%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_pd=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_pd, names_sep = ".")%>%
  mutate(shift= T2-T1)

#relative abundance-weighted shift
cm_sp_diff90=left_join(cm_sp_df_diff, cm_sp_ra, by="Species")%>%
  mutate(wgt_shift=shift*rel_abund)

cm_sp90diff=sum(cm_sp_diff90$wgt_shift) #abundance-weighted shift
cm_p90diff=mean(cm_mod90pred[24:28,]$pred_JD)-mean(cm_mod90pred[1:5,]$pred_JD)

cm_sp90diff/cm_p90diff #proportion of phenological shift by species phenology
1-(cm_sp90diff/cm_p90diff) #proportion of phenological shift by composition
cm_p90diff*(1-(cm_sp90diff/cm_p90diff)) #phenological shift (composition)

#visualize data####

##shifts####
##10% PD####
plot(cm_10pd_site$Julian~cm_10pd_site$YR, type="l", main="Cape May", ylab="10% passage date", xlab="year")

lines(cm_m1pred$pred_JD~cm_m1pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1995, x1=1999,y0=mean(cm_10pd_site[1:5,]$Julian), y1=mean(cm_10pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2018, x1=2022,y0=mean(cm_10pd_site[24:28,]$Julian), y1=mean(cm_10pd_site[24:28,]$Julian), lwd=2, col="orange")

segments(x0=1995, x1=1999,y0=mean(cm_mod10pred[1:5,]$pred_JD), y1=mean(cm_mod10pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2018, x1=2022,y0=mean(cm_mod10pred[24:28,]$pred_JD), y1=mean(cm_mod10pred[24:28,]$pred_JD), lwd=2, col="red")

##50% PD####
plot(cm_50pd_site$Julian~cm_50pd_site$YR, type="l", main="Cape May", ylab="50% passage date", xlab="year")

lines(cm_m5pred$pred_JD~cm_m5pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1995, x1=1999,y0=mean(cm_50pd_site[1:5,]$Julian), y1=mean(cm_50pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2018, x1=2022,y0=mean(cm_50pd_site[24:28,]$Julian), y1=mean(cm_50pd_site[24:28,]$Julian), lwd=2, col="orange")

segments(x0=1995, x1=1999,y0=mean(cm_mod50pred[1:5,]$pred_JD), y1=mean(cm_mod50pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2018, x1=2022,y0=mean(cm_mod50pred[24:28,]$pred_JD), y1=mean(cm_mod50pred[24:28,]$pred_JD), lwd=2, col="red")

##90% PD####
plot(cm_90pd_site$Julian~cm_90pd_site$YR, type="l", main="Cape May", ylab="90% passage date", xlab="year")

lines(cm_m9pred$pred_JD~cm_m9pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1995, x1=1999,y0=mean(cm_90pd_site[1:5,]$Julian), y1=mean(cm_90pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2018, x1=2022,y0=mean(cm_90pd_site[24:28,]$Julian), y1=mean(cm_90pd_site[24:28,]$Julian), lwd=2, col="orange")

segments(x0=1995, x1=1999,y0=mean(cm_mod90pred[1:5,]$pred_JD), y1=mean(cm_mod90pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2018, x1=2022,y0=mean(cm_mod90pred[24:28,]$pred_JD), y1=mean(cm_mod90pred[24:28,]$pred_JD), lwd=2, col="red")

#contributions####
cm_datphen=matrix(c(10,50,90,
                    1.16, 0.78,0.09, 
                    2.89, 1.83,0.63,
                    -1.73, -1.05, -0.57), nrow=3, ncol=4, byrow=FALSE, 
                  dimnames= list(c("1", "2", "3"),
                                 c("metric","overall", "species phenology", "composition")))
cm_dat_phen=as.data.frame(cm_datphen)%>%
  pivot_longer(cols=2:4, names_to = "shift", values_to="value" )

ggplot(cm_dat_phen)+geom_col(aes(x=shift, y=value))+ggtitle("Cape May")+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_x_discrete(limits=c("overall", "species phenology", "composition"))+
  facet_wrap(~metric)+theme_classic()+abline(h=0)+ylab("shift (days)")

#compositional shifts####
cm_sp_tot=cm_sp%>%
  group_by(YR, Species)%>%
  summarise(sp_tot=sum(Count))%>%
  group_by(YR)%>%
  mutate(rel_abund=sp_tot/sum(sp_tot))

cm_sp_abundance2=cm_sp_tot%>%
  nest(-Species)%>%
  mutate(model= list(brm(bf(rel_abund~ 1+s(YR)), 
                         data=data, 
                         iter=2500, 
                         family=zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit"),
                         control=list(adapt_delta=0.99))),
         preds=purrr::map(model, predict))%>%
  unnest(preds)%>%as.data.frame()

cm_sp_relab=cm_sp_abundance2%>%select(-data, -model)
cm_sp_rela_preds=as.data.frame(cm_sp_relab$preds[,1])
colnames(cm_sp_rela_preds)[1]="mean_rela"

cm_sp_df_rela=cbind(cm_sp_rela_preds, cm_sp_tot)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1995:1999) ~ "T1",
                          YR%in%c(2018:2022) ~ "T2"))

cm_sp_df_rela_diff=cm_sp_df_rela%>%
  select(Species, mean_rela,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_rela=mean(mean_rela))%>%
  pivot_wider(names_from = period, values_from=mean_rela, names_sep = ".")%>%
  mutate(shift= T2-T1)

#determine passage order of species (based on 50% passage date)
cm_ord=cm_50pd_sp%>%group_by(Species)%>%
  summarise(min_time=min(Julian))%>%arrange(min_time)

ggplot(cm_sp_df_rela_diff)+geom_col(aes(x=Species, y=shift))+theme_classic()+
  geom_hline(yintercept = 0, linetype = 2)+ggtitle("Cape May")+
  ylab("relative abundance shifts")+
  scale_x_discrete(limits=c("BW","OS","AK", "ML", "BE", "PG",
                            "SS", "CH", "NH", "BV", "TV", "NG",
                            "GE", "RS", "RT"),
                   labels=c("BWHA", "OSPR", "AMKE", "MERL", "BAEA",
                            "PEFA", "SSHA", "COHA", "NOHA", "BLVU",
                            "TUVU", "NOGO", "GOEA", "RSHA", "RTHA"))
