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
  geom_bracket(
    xmin = 244 , xmax = 333, y.position = 76,
    label = paste("CM"), type="expression"
  )+scale_color_manual(values=c("darkgreen", "gold"))

#Figure 2####

hm_datphen=matrix(c(10,50,90,
                    -0.54,-8.9,-3.74, 
                    1.41,3.14,1.62,
                    -1.95, -12.03, -5.35), nrow=3, ncol=4, byrow=FALSE, 
                  dimnames= list(c("1", "2", "3"),
                                 c("metric","overall", "species phenology", "composition")))
hm_dat_phen=as.data.frame(hm_datphen)%>%
  pivot_longer(cols=2:4, names_to = "shift", values_to="value" )

ggplot(hm_dat_phen)+geom_col(aes(x=shift, y=value))+ggtitle("Hawk Mountain")+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_x_discrete(limits=c("overall", "species phenology", "composition"))+
  facet_wrap(~metric)+theme_classic()+ylab("shift (days)")

#Figure 3####
#3a
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

#3b
hms_ord50=hms_50pd_sp%>%group_by(Species)%>%
  summarise(min_time=min(Julian))%>%arrange(min_time)

hm50_abund=ggplot(hm_sp_abundance_diffs)+geom_col(aes(x=Species, y=shift, fill=diet))+theme_classic()+
  geom_hline(yintercept = 0, linetype = 2)+ggtitle("Hawk Mountain (50% passage date)")+
  ylab("relative abundance shifts")+
  scale_x_discrete(limits=c("BWHA", "BAEA", "OSPR","AMKE","PEFA", "MERL",
                            "COHA","SSHA", "NOHA","BLVU", "TUVU", "RSHA",
                            "RTHA", "GOEA", "NOGO"))

hm50_phen=ggplot(hm_sp_df_diff50)+geom_col(aes(x=Species, y=shift, fill=diet))+theme_classic()+
  geom_hline(yintercept = 0, linetype = 2)+
  ylab("phenological shifts")+
  scale_x_discrete(limits=c("BWHA", "BAEA", "OSPR","AMKE","PEFA", "MERL",
                            "COHA","SSHA", "NOHA","BLVU", "TUVU", "RSHA",
                            "RTHA", "GOEA", "NOGO"))


ggarrange(hm50_abund, hm50_phen, nrow=2, common.legend = TRUE, legend = "right")

#3c
hms_ord90=hms_90pd_sp%>%group_by(Species)%>%
  summarise(min_time=min(Julian))%>%arrange(min_time)

hm90_abund=ggplot(hm_sp_abundance_diffs)+geom_col(aes(x=Species, y=shift, fill=diet))+theme_classic()+
  geom_hline(yintercept = 0, linetype = 2)+ggtitle("Hawk Mountain (90% passage date)")+
  ylab("relative abundance shifts")+
  scale_x_discrete(limits=c("BWHA", "OSPR","AMKE","PEFA", "MERL",
                            "SSHA", "COHA", "TUVU", "NOHA", "BAEA","BLVU",
                            "RSHA", "RTHA", "NOGO", "GOEA"))

hm90_phen=ggplot(hm_sp_df_diff90)+geom_col(aes(x=Species, y=shift, fill=diet))+theme_classic()+
  geom_hline(yintercept = 0, linetype = 2)+
  ylab("phenological shifts")+
  scale_x_discrete(limits=c("BWHA", "OSPR","AMKE","PEFA", "MERL",
                            "SSHA", "COHA", "TUVU", "NOHA", "BAEA","BLVU",
                            "RSHA", "RTHA", "NOGO", "GOEA"))


ggarrange(hm90_abund, hm90_phen, nrow=2, common.legend = TRUE, legend = "right")

