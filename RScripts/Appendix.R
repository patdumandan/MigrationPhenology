#appendix####
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

#determine passage order of species (based on 10% passage date)
hms_ord10=hms_10pd_sp%>%group_by(Species)%>%
  summarise(min_time=min(Julian))%>%arrange(min_time)

hm10_abund=ggplot(hm_sp_abundance_diffs)+geom_col(aes(x=Species, y=shift, fill=diet))+theme_classic()+
  geom_hline(yintercept = 0, linetype = 2)+ggtitle("Hawk Mountain (10% passage date)")+
  ylab("relative abundance shifts")+
  scale_x_discrete(limits=c("BLVU", "BAEA", "AMKE", "OSPR","BWHA", "MERL",
                            "NOHA","PEFA", "COHA","SSHA", "TUVU", "NOGO",
                            "RSHA", "RTHA", "GOEA"))

hm10_phen=ggplot(hm_sp_df_diff10)+geom_col(aes(x=Species, y=shift, fill=diet))+theme_classic()+
  geom_hline(yintercept = 0, linetype = 2)+
  ylab("phenological shifts")+
  scale_x_discrete(limits=c("BLVU", "BAEA", "AMKE", "OSPR","BWHA", "MERL",
                            "NOHA","PEFA", "COHA","SSHA", "TUVU", "NOGO",
                            "RSHA", "RTHA", "GOEA"))

ggarrange(hm10_abund, hm10_phen, nrow=2, common.legend = TRUE, legend = "right")

#appendix####
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
                        Species=="COHA" ~ "mammal",
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
  geom_point()+
  geom_text(label=hm_diffs_all$data_labels, size=3, col="black", hjust = 1.8)+
  theme_classic()+geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)+ggtitle("Hawk Mountain")+
  ylab("relative abundance shifts")+facet_wrap(~metric)+xlab("phenological shifts")+
  scale_color_viridis_d()

hmplot+guides(size=FALSE)

