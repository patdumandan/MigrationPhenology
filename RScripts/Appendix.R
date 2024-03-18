#appendix####
##shifts####
##10% PD####

#CAPE MAY####
#visualize data####

##shifts####
##10% PD####
plot(cm_10pd_site$Julian~cm_10pd_site$YR, type="l", main="Cape May", ylab="10% passage date", xlab="year")

lines(cm_m1pred$pred_JD~cm_m1pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1990, x1=1994,y0=mean(cm_10pd_site[1:5,]$Julian), y1=mean(cm_10pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2014, x1=2018,y0=mean(cm_10pd_site[25:29,]$Julian), y1=mean(cm_10pd_site[25:29,]$Julian), lwd=2, col="orange")

segments(x0=1990, x1=1994,y0=mean(cm_mod10pred[1:5,]$pred_JD), y1=mean(cm_mod10pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2014, x1=2018,y0=mean(cm_mod10pred[25:29,]$pred_JD), y1=mean(cm_mod10pred[25:29,]$pred_JD), lwd=2, col="red")

##50% PD####
plot(cm_50pd_site$Julian~cm_50pd_site$YR, type="l", main="Cape May", ylab="50% passage date", xlab="year")

lines(cm_m5pred$pred_JD~cm_m5pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1990, x1=1994,y0=mean(cm_50pd_site[1:5,]$Julian), y1=mean(cm_50pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2014, x1=2018,y0=mean(cm_50pd_site[25:29,]$Julian), y1=mean(cm_50pd_site[25:29,]$Julian), lwd=2, col="orange")

segments(x0=1990, x1=1994,y0=mean(cm_mod50pred[1:5,]$pred_JD), y1=mean(cm_mod50pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2014, x1=2018,y0=mean(cm_mod50pred[25:29,]$pred_JD), y1=mean(cm_mod50pred[25:29,]$pred_JD), lwd=2, col="red")

##90% PD####
plot(cm_90pd_site$Julian~cm_90pd_site$YR, type="l", main="Cape May", ylab="90% passage date", xlab="year")

lines(cm_m9pred$pred_JD~cm_m9pred$YR, lty=2, col="grey", lwd=2)

segments(x0=1990, x1=1994,y0=mean(cm_90pd_site[1:5,]$Julian), y1=mean(cm_90pd_site[1:5,]$Julian), lwd=2, col="orange")
segments(x0=2014, x1=2018,y0=mean(cm_90pd_site[25:29,]$Julian), y1=mean(cm_90pd_site[25:29,]$Julian), lwd=2, col="orange")

segments(x0=1990, x1=1994,y0=mean(cm_mod90pred[1:5,]$pred_JD), y1=mean(cm_mod90pred[1:5,]$pred_JD), lwd=2, col="red")
segments(x0=2014, x1=2018,y0=mean(cm_mod90pred[25:29,]$pred_JD), y1=mean(cm_mod90pred[25:29,]$pred_JD), lwd=2, col="red")

#HMS####
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
  summarise(mean_time=mean(Julian))%>%arrange(mean_time)

hm50_abund=ggplot(hm_sp_abundance_diffs)+geom_col(aes(x=Species, y=shift, fill=diet))+theme_classic()+
  geom_hline(yintercept = 0, linetype = 2)+
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


fig3c=ggarrange(hm50_abund, hm90_phen, nrow=2, common.legend = TRUE, legend = "right")
annotate_figure(fig3c, top = text_grob("Hawk Mountain (90% passage date)", 
                                       face = "bold", size = 14))

fig3b=ggarrange(hm50_abund, hm50_phen, nrow=2, common.legend = TRUE, legend = "right")
annotate_figure(fig3b, top = text_grob("Hawk Mountain (50% passage date)", 
                                       face = "bold", size = 14))

fig3a=ggarrange(hm50_abund, hm10_phen, nrow=2, common.legend = TRUE, legend = "right")
annotate_figure(fig3a, top = text_grob("Hawk Mountain (10% passage date)", 
                                       face = "bold", size = 14))



