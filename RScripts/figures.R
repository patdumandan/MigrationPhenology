#Figure 1####
par(mfrow=c(3,1))

par(mar = c(2.25, 2.5, 1, 2.5))
x = seq(150, 600, by=1) 

y = dnorm(x, mean=220, sd=20) 
sy= dnorm(x, mean=340, sd=20) 


plot(x, y, type="l",lty="dashed",lwd=2, xlab="", xaxt="n", 
     ylab="Abundance",  xlim=c(150,410), ylim=c(0,0.025), yaxt="n", col="grey") 
lines(x, sy, lwd=2) 
abline(v=220, lty="dashed", lwd=2, col="grey")
abline(v=340, lty="dashed", lwd=2, col="grey")

arrows(x0=230, y0=0.02,x1=330, y1=0.02, lwd=2)
text(x=180, y=0.024, "assemblage shift")

x = seq(150, 600, by=1) 

ay = dnorm(x, mean=210, sd=20) 
by = dnorm(x, mean=230, sd=20)

aay = dnorm(x, mean=330, sd=20) 
bby = dnorm(x, mean=350, sd=20)


plot(x, ay, type="l", lty="dashed", lwd=2, xlab="",xaxt="n", 
     ylab="Abundance",  xlim=c(150,410)
     , ylim=c(0,0.025), yaxt="n", col="red") 
lines(x, by, lwd=2, col="blue", lty="dashed") 
lines(x, bby, lwd=2, col="blue") 
lines(x, aay, lwd=2, col="red") 
abline(v=220, lty="dashed", lwd=2, col="grey")
abline(v=340, lty="dashed", lwd=2, col="grey")

arrows(x0=222, y0=0.022,x1=337, y1=0.022, col="black", lwd=2)
text(x=180, y=0.022, "species phenology")

x = seq(150, 600, by=1) 

ay = dnorm(x, mean=210, sd=20) 
by = dnorm(x, mean=350, sd=20)

aay = dnorm(x, mean=210, sd=30) 
byy = dnorm(x, mean=350, sd=10) 


plot(x, ay, type="l", lty="dashed", lwd=2, xlab="time", 
     ylab="Abundance",  xlim=c(150,410),yaxt="n",
     ylim=c(0,0.05), col="red") 
lines(x, aay, lwd=2, col="red") 
lines(x, by, lwd=2, col="blue", lty="dashed") 
lines(x, byy, lwd=2, col="blue") 


abline(v=220, lty="dashed", lwd=2, col="grey")
abline(v=340, lty="dashed", lwd=2, col="grey")

arrows(x0=222, y0=0.038,x1=337, y1=0.038, col="black", lwd=2)
text(x=180, y=0.038, "composition")
arrows(x0=180, y0=0.02,x1=180, y1=0.01, col="black", lwd=2)
arrows(x0=320, y0=0.020,x1=320, y1=0.03, col="black", lwd=2)


#Figure 2####
cape=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/CapeMay.csv")

CapeMay=cape%>%
  select(-Duration, -Observer,-RL,-GO,-MK,-SE,-SW,-UA, -UB, -UF,  -ZT, -TOTAL)%>%
  mutate(yr_tot=rowSums(.[2:16], na.rm=T),
         date=as.Date(Date, format="%m/%d/%Y"),
         YR=year(date), MO=month(date), DAY=day(date))%>%
  filter(!YR<1990, !YR>2018)

CapeMay$Julian=as.POSIXlt(CapeMay$date)$yday
CapeMay[, 2:16][is.na(CapeMay[, 2:16])] <- 0

cm_sp=CapeMay%>%select(-yr_tot)%>%
  pivot_longer(cols=2:16, names_to = "Species", values_to="Count")

cm_sp_day_mean=cm_sp%>%
  group_by(Species,YR, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(YR)%>%
  group_by(Julian)%>%
  summarise(ave_ct=mean(total))%>%
  mutate(site="Cape May (CM)")

hms=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/HMSDay.csv")

#combine baea and goea data
#select 11 years after first count

HMS=hms%>%
  mutate(BAEA=rowSums(.[9:11]), GOEA=rowSums(.[21:23]))%>%
  select(-OBSERVERS, -GOEA.I, -GOEA.A, -GOEA.U, -BAEA.I, -BAEA.U, -BAEA.A, -RAPTORS, -GYRF, -MIKI, -STKI,-SWHA,
         -UNID.ACCIPITER,-UNID.BUTEO, -UNID.EAGLE, -UNID.FALCON, -UNID.RAPTOR)%>%
  mutate(yr_tot=rowSums(.[6:21], na.rm=T), date=make_date(YR, MO, DAY))%>%filter(!(YR<1990))

HMS$Julian=as.POSIXlt(HMS$date)$yday
HMS[, 6:21][is.na(HMS[, 6:21])] <- 0

hms_sp=HMS%>%
  pivot_longer(cols=6:21, names_to = "Species", values_to="Count")%>%
  filter(!YR<1990, !Species=="RLHA")

hms_sp_day_mean=hms_sp%>%
  group_by(Species,YR, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(YR)%>%
  group_by(Julian)%>%
  summarise(ave_ct=mean(total))%>%
  mutate(site="Hawk Mountain (HMS)")

sp_day_all=rbind(hms_sp_day_mean, cm_sp_day_mean)


fig1=ggplot()+geom_line(sp_day_all, mapping=aes(x=Julian, y=ave_ct, col=site))+
# scale_color_viridis(discrete=TRUE)+
ylab("average daily count")+
  xlab("Julian date")+
  theme_classic()+
 # annotate(geom="vline", x=c(226, 348), xintercept=c(226, 348), linetype=c("dashed", "dashed"))+
  #annotate(geom="text", x=c(221, 343), y=c(50, 50), label=c("Hawk Mountain start", "Hawk Mountain end"), angle=90)+
  #annotate(geom="vline", x=c(244, 333), xintercept=c(244, 333), linetype=c("dashed", "dashed"))+
  #annotate(geom="text", x=c(239, 328), y=c(50, 50), label=c("Cape May start", "Cape May end"), angle=90)+
  ggtitle("community phenology")+
  geom_bracket(
    xmin = 226 , xmax = 348, y.position = 85,
    label = paste("HMS"), type="expression"
  )+
  geom_vline(xintercept=271, col="gold", lty=2)+
  geom_vline(xintercept=281, col="darkgreen", lty=2)+
  geom_bracket(
    xmin = 244 , xmax = 333, y.position = 76,
    label = paste("CM"), type="expression"
  )+scale_color_manual(values=c("darkgreen", "gold"))

#Figure 3####

hm_datphen=matrix(c(10,50,90,
                    -0.39,-8.9,-3.21, 
                    1.27,3.14,1.65,
                    -3.21, -12.03, -4.86), nrow=3, ncol=4, byrow=FALSE, 
                  dimnames= list(c("1", "2", "3"),
                                 c("metric","overall", "species phenology", "composition")))
hm_dat_phen=as.data.frame(hm_datphen)%>%
  pivot_longer(cols=2:4, names_to = "shift", values_to="value" )

ggplot(hm_dat_phen)+geom_col(aes(x=shift, y=value))+ggtitle("Hawk Mountain")+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_x_discrete(limits=c("overall", "species phenology", "composition"))+
  facet_wrap(~metric)+theme_classic()+ylab("shift (days)")

cm_datphen=matrix(c(10,50,90,
                    3.37, 4.24,4.11, 
                    3.93, 2.62,2.45,
                    -0.55, 1.63, 1.66), nrow=3, ncol=4, byrow=FALSE, 
                  dimnames= list(c("1", "2", "3"),
                                 c("metric","overall", "species phenology", "composition")))
cm_dat_phen=as.data.frame(cm_datphen)%>%
  pivot_longer(cols=2:4, names_to = "shift", values_to="value" )

ggplot(cm_dat_phen)+geom_col(aes(x=shift, y=value))+ggtitle("Cape May")+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_x_discrete(limits=c("overall", "species phenology", "composition"))+
  facet_wrap(~metric)+theme_classic()+ylab("shift (days)")

#Figure 4####
#3a

hm_sp_df_rela=left_join(hm_sp_abundance_meanpreds, hms_sp_tot)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))

hm_sp_df_rela_diff=hm_sp_df_rela%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_rela=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_rela, names_sep = ".")%>%
  mutate(shift= T2-T1)%>%
  mutate(diet=case_when(Species=="AMKE" ~ "insect",
                        Species=="BWHA" ~ "mammal",
                        Species=="OSPR" ~ "fish",
                        Species=="MERL" ~ "bird",
                        Species=="BAEA" ~ "fish",
                        Species=="PEFA" ~ "bird",
                        Species=="SSHA" ~ "bird",
                        Species=="COHA" ~ "bird",
                        Species=="NOHA" ~ "mammal",
                        Species=="BLVU" ~ "carrion",
                        Species=="TUVU" ~ "carrion",
                        Species=="NOGO" ~ "bird",
                        Species=="GOEA" ~ "mammal",
                        Species=="RSHA" ~ "mammal",
                        Species=="RTHA" ~ "mammal"))

hm_sp_diff_rel=hm_sp_df_rela_diff%>%
  select(Species, shift, diet)%>%
  rename(rela_shift=shift)

hm_sp_diff10_time=hm_sp_diff10%>%
  mutate(metric="10")%>%
  rename(time_shift=shift)%>%
  select(Species, time_shift, wgt_shift, metric)

hm_sp_diff50_time=hm_sp_diff50%>%
  mutate(metric="50")%>%
  rename(time_shift=shift)%>%
  select(Species, time_shift, wgt_shift, metric)

hm_sp_diff90_time=hm_sp_diff90%>%
  mutate(metric="90")%>%
  rename(time_shift=shift)%>%
  select(Species, time_shift, wgt_shift, metric)

hm_diff_all=do.call("rbind", list(hm_sp_diff10_time, hm_sp_diff50_time, hm_sp_diff90_time))%>%
  left_join(hm_sp_diff_rel, by="Species")

hm_tot_ave=hms_sp_tot%>%
  group_by(Species)%>%
  summarise(ave_tot=mean(rel_abund))

hm_diffs_all=left_join(hm_diff_all, hm_tot_ave, by="Species")%>%
  mutate(col_grp=case_when(ave_tot > 0.1 ~ "a"))%>%
  mutate(data_labels = ifelse(col_grp == "a", Species, NA))

hmplot=ggplot(hm_diffs_all, aes(x=time_shift, y=rela_shift, col=diet, size=ave_tot))+
  geom_point()+facet_wrap(~metric)+
  geom_text(label=hm_diffs_all$data_labels, size=3, col="black", hjust = 1.8)+
  theme_classic()+geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)+ggtitle("Hawk Mountain")+
  ylab("relative abundance shifts")+xlab("phenological shifts")+
  scale_color_viridis_d()

hmplot+guides(size=FALSE)

cm_sp_df_rela=right_join(cm_sp_abundance_meanpreds, cm_sp_tot)%>%as.data.frame()%>%
  mutate(period=case_when(YR%in%c(1990:1994) ~ "T1",
                          YR%in%c(2014:2018) ~ "T2"))

cm_sp_df_rela_diff=cm_sp_df_rela%>%
  select(Species, mean_preds,period)%>%
  filter(!is.na(period))%>%
  group_by(Species, period)%>%
  summarise(mean_rela=mean(mean_preds))%>%
  pivot_wider(names_from = period, values_from=mean_rela, names_sep = ".")%>%
  mutate(shift= T2-T1)%>%
  mutate(diet=case_when(Species=="AK" ~ "insect",
                        Species=="BW" ~ "mammal",
                        Species=="OS" ~ "fish",
                        Species=="ML" ~ "bird",
                        Species=="BE" ~ "fish",
                        Species=="PG" ~ "bird",
                        Species=="SS" ~ "bird",
                        Species=="CH" ~ "bird",
                        Species=="NH" ~ "mammal",
                        Species=="BV" ~ "carrion",
                        Species=="TV" ~ "carrion",
                        Species=="NG" ~ "bird",
                        Species=="GE" ~ "mammal",
                        Species=="RS" ~ "mammal",
                        Species=="RT" ~ "mammal"))%>%
  mutate(spcode=case_when(Species=="AK" ~ "AMKE",
                          Species=="BW" ~ "BWHA",
                          Species=="OS" ~ "OSPR",
                          Species=="ML" ~ "MERL",
                          Species=="BE" ~ "BAEA",
                          Species=="PG" ~ "PEFA",
                          Species=="SS" ~ "SSHA",
                          Species=="CH" ~ "COHA",
                          Species=="NH" ~ "NOHA",
                          Species=="BV" ~ "BLVU",
                          Species=="TV" ~ "TUVU",
                          Species=="NG" ~ "NOGO",
                          Species=="GE" ~ "GOEA",
                          Species=="RS" ~ "RSHA",
                          Species=="RT" ~ "RTHA"))

cm_sp_diff_rel=cm_sp_df_rela_diff%>%
  select(Species, spcode, shift, diet)%>%
  rename(rela_shift=shift)

cm_sp_diff10_time=cm_sp_diff10%>%
  mutate(metric="10")%>%
  rename(time_shift=shift)%>%
  select(Species, time_shift, wgt_shift, metric)

cm_sp_diff50_time=cm_sp_diff50%>%
  mutate(metric="50")%>%
  rename(time_shift=shift)%>%
  select(Species, time_shift, wgt_shift, metric)

cm_sp_diff90_time=cm_sp_diff90%>%
  mutate(metric="90")%>%
  rename(time_shift=shift)%>%
  select(Species, time_shift, wgt_shift, metric)

cm_diff_all=do.call("rbind", list(cm_sp_diff10_time, cm_sp_diff50_time, cm_sp_diff90_time))%>%
  inner_join(cm_sp_diff_rel, by="Species")

cm_tot_ave=cm_sp_tot%>%
  group_by(Species)%>%
  summarise(ave_tot=mean(rel_abund))

cm_diffs_all=left_join(cm_diff_all, cm_tot_ave, by="Species")%>%
  mutate(col_grp=case_when(ave_tot > 0.1 ~ "a"))%>%
  mutate(data_labels = ifelse(col_grp == "a", spcode, NA))

cmplot=ggplot(cm_diffs_all, aes(x=time_shift, y=rela_shift, col=diet, size=ave_tot))+
  geom_point()+facet_wrap(~metric)+
  geom_text(label=cm_diffs_all$data_labels, size=3, col="black", hjust = 1.8)+
  theme_classic()+geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)+ggtitle("Cape May")+
  ylab("relative abundance shifts")+xlab("phenological shifts")+
  scale_color_viridis_d()

cmplot+guides(size=FALSE)

#Appendix####
#Fig 1
#HMS

hms_10pd_site=cbind(hms_10pd_site, mod10pred)%>%mutate(metric="onset")
hm10=ggplot(hms_10pd_site, aes(y=hms_10pd_site$Julian, x=hms_10pd_site$YR))+geom_point()+
  annotate("rect", ymin=min(hms_10pd_site$Julian), ymax=max(hms_10pd_site$Julian),
           xmin=1990, xmax=1994, alpha=0.2, fill="red")+xlab("Year")+ylab("10% passage date (Julian)")+
  annotate("rect", ymin=min(hms_10pd_site$Julian), ymax=max(hms_10pd_site$Julian),
           xmin=2014, xmax=max(hms_10pd_site$YR), alpha=0.2, fill="red")+
  theme_classic()+
  geom_line(aes(x=YR, y=pred_JD), linetype="dashed")

hms_50pd_site=cbind(hms_50pd_site, mod50pred)%>%mutate(metric="mean")
hm50=ggplot(hms_50pd_site, aes(y=hms_50pd_site$Julian, x=hms_50pd_site$YR))+geom_point()+
  annotate("rect", ymin=min(hms_50pd_site$Julian), ymax=max(hms_50pd_site$Julian),
           xmin=1990, xmax=1994, alpha=0.2, fill="red")+xlab("Year")+ylab("50% passage date (Julian)")+
  annotate("rect", ymin=min(hms_50pd_site$Julian), ymax=max(hms_50pd_site$Julian),
           xmin=2014, xmax=max(hms_50pd_site$YR), alpha=0.2, fill="red")+
  theme_classic()+
  geom_line(aes(x=YR, y=pred_JD), linetype="dashed")

hms_90pd_site=cbind(hms_90pd_site, mod90pred)%>%mutate(metric="end")
hm90=ggplot(hms_90pd_site, aes(y=hms_90pd_site$Julian, x=hms_90pd_site$YR))+geom_point()+
  annotate("rect", ymin=min(hms_90pd_site$Julian), ymax=max(hms_90pd_site$Julian),
           xmin=1990, xmax=1994, alpha=0.2, fill="red")+xlab("Year")+ylab("90% passage date (Julian)")+
  annotate("rect", ymin=min(hms_90pd_site$Julian), ymax=max(hms_90pd_site$Julian),
           xmin=2014, xmax=max(hms_90pd_site$YR), alpha=0.2, fill="red")+
  theme_classic()+
  geom_line(aes(x=YR, y=pred_JD), linetype="dashed")

hm1=ggarrange(hm10, hm50, hm90, nrow=3, ncol=1)
annotate_figure(hm1, top=text_grob("Hawk Mountain Sanctuary", face="bold"))

#Cape May
cm_10pd_site=cbind(cm_10pd_site, cm_mod10pred)%>%mutate(metric="onset")
cm10=ggplot(cm_10pd_site, aes(y=cm_10pd_site$Julian, x=cm_10pd_site$YR))+geom_point()+
  annotate("rect", ymin=min(cm_10pd_site$Julian), ymax=max(cm_10pd_site$Julian),
           xmin=1990, xmax=1994, alpha=0.2, fill="red")+xlab("Year")+ylab("10% passage date (Julian)")+
  annotate("rect", ymin=min(cm_10pd_site$Julian), ymax=max(cm_10pd_site$Julian),
           xmin=2014, xmax=max(cm_10pd_site$YR), alpha=0.2, fill="red")+
  theme_classic()+
  geom_line(aes(x=YR, y=pred_JD), linetype="dashed")

cm_50pd_site=cbind(cm_50pd_site, cm_mod50pred)%>%mutate(metric="mean")
cm50=ggplot(cm_50pd_site, aes(y=cm_50pd_site$Julian, x=cm_50pd_site$YR))+geom_point()+
  annotate("rect", ymin=min(cm_50pd_site$Julian), ymax=max(cm_50pd_site$Julian),
           xmin=1990, xmax=1994, alpha=0.2, fill="red")+xlab("Year")+ylab("50% passage date (Julian)")+
  annotate("rect", ymin=min(cm_50pd_site$Julian), ymax=max(cm_50pd_site$Julian),
           xmin=2014, xmax=max(cm_50pd_site$YR), alpha=0.2, fill="red")+
  theme_classic()+
  geom_line(aes(x=YR, y=pred_JD), linetype="dashed")

cm_90pd_site=cbind(cm_90pd_site, cm_mod90pred)%>%mutate(metric="end")
cm90=ggplot(cm_90pd_site, aes(y=cm_90pd_site$Julian, x=cm_90pd_site$YR))+geom_point()+
  annotate("rect", ymin=min(cm_90pd_site$Julian), ymax=max(cm_90pd_site$Julian),
           xmin=1990, xmax=1994, alpha=0.2, fill="red")+xlab("Year")+ylab("90% passage date (Julian)")+
  annotate("rect", ymin=min(cm_90pd_site$Julian), ymax=max(cm_90pd_site$Julian),
           xmin=2014, xmax=max(cm_90pd_site$YR), alpha=0.2, fill="red")+
  theme_classic()+
  geom_line(aes(x=YR, y=pred_JD), linetype="dashed")

cm1=ggarrange(cm10, cm50, cm90, nrow=3, ncol=1)
annotate_figure(cm1, top=text_grob("Cape May", face="bold"))
