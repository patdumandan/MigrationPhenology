#Figure 1####
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

#Figure 2####

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

#Figure 3####
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
  geom_point()+facet_wrap(~metric)+
  geom_text(label=hm_diffs_all$data_labels, size=3, col="black", hjust = 1.8)+
  theme_classic()+geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)+ggtitle("Hawk Mountain")+
  ylab("relative abundance shifts")+xlab("phenological shifts")+
  scale_color_viridis_d()

hmplot+guides(size=FALSE)
