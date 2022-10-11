gos=gos%>%mutate(Julian=as.POSIXlt(Date)$yday,
                 diet=case_when(Species=="AK" ~"bird",
                                Species=="BE" ~"other",
                                Species=="CH" ~"bird",
                                Species=="BV" ~ "other",
                                Species=="TV" ~"other",
                                Species=="OS" ~"other",
                                Species== "NH" ~"mammal",
                                Species=="SS" ~"bird",
                                Species=="NG"~ "bird",
                                Species=="RS" ~"mammal",
                                Species=="BW"~"mammal",
                                Species=="RT" ~"mammal",
                                Species=="RL" ~"mammal",
                                Species=="SW"~"mammal",
                                Species=="FH"~"mammal",
                                Species=="GE"~"mammal",
                                Species=="ML"~"bird",
                                Species=="PG"~"bird",
                                Species=="PR"~"mammal",
                                Species=="MK"~"other"),
                 mammals=case_when(Species=="AK" ~"34.6",
                                      Species=="BE" ~"16",
                                      Species=="CH" ~"23.3",
                                      Species=="BV" ~ "92.5",
                                      Species=="TV" ~"79.7",
                                      Species=="OS" ~"0",
                                      Species== "NH" ~"48.2",
                                      Species=="SS" ~"2.72",
                                      Species=="NG"~ "24.7",
                                      Species=="RS" ~"39.5",
                                      Species=="BW"~"42.2",
                                      Species=="RT" ~"73.8",
                                      Species=="RL" ~"85.3",
                                      Species=="SW"~"48",
                                      Species=="FH"~"81.9",
                                      Species=="GE"~"80.6",
                                      Species=="ML"~"1.28",
                                      Species=="PG"~"5.11",
                                      Species=="PR"~"42.3",
                                      Species=="MK"~"2.22"),
                 birds=case_when(Species=="AK" ~"13.5",
                                  Species=="BE" ~"23.4",
                                  Species=="CH" ~"71.4",
                                  Species=="BV" ~ "7.5",
                                  Species=="TV" ~"6.42",
                                  Species=="OS" ~"0",
                                  Species== "NH" ~"38.6",
                                  Species=="SS" ~"89.4",
                                  Species=="NG"~ "74.2",
                                  Species=="RS" ~"6.43",
                                  Species=="BW"~"16.5",
                                  Species=="RT" ~"17",
                                  Species=="RL" ~"9.02",
                                  Species=="SW"~"10.6",
                                  Species=="FH"~"9.89",
                                  Species=="GE"~"17.8",
                                  Species=="ML"~"90.3",
                                  Species=="PG"~"93.5",
                                  Species=="PR"~"45",
                                  Species=="MK"~"0"))
                 %>%
  filter(!Species%in%c("UA", "UB", "UF", "UE", "UR"))

gos$mammals=as.numeric(gos$mammals)
gos$birds=as.numeric(gos$birds)
str(gos)

gos_sp=gos%>%mutate(Julian=as.POSIXlt(Date)$yday)%>%
  group_by(year, Species, diet, mammals, birds)%>%
  summarise(yr_tot=sum(Count))%>%filter(!Species%in%c("UA", "UB", "UE", "UF", "UR"))%>%
  arrange(year)

g1=ggplot(gos_sp, aes(x=birds, y=log(yr_tot)))+geom_point()+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")+ggtitle("Goshutes")

g2=ggplot(gos_sp, aes(x=mammals, y=log(yr_tot)))+geom_point()+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")+ggtitle("Goshutes")

gos_sp_max=gos%>%mutate(Julian=as.POSIXlt(Date)$yday)%>%
  group_by(year, Species, diet, mammals, birds)%>%
  filter(Count==max(Count), !Species%in%c("UA", "UB", "UE", "UF", "UR"))%>%slice_head(n=1)%>%
  arrange(year)

g4=ggplot(gos_sp_max, aes(x=mammals, y=Julian))+geom_point()+
  stat_smooth(method="lm")+geom_hline(yintercept=268, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()+ggtitle("Goshutes")

g3=ggplot(gos_sp_max, aes(x=birds, y=Julian))+geom_point()+
  stat_smooth(method="lm")+geom_hline(yintercept=268, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()+ggtitle("Goshutes")

ggarrange(g1,g2,g3,g4)

gos_sp_mammd=gos_sp_max%>%filter(diet=="mammal")
gos_sp_mammc=gos_sp%>%filter(diet=="mammal")

ggplot(gos_sp_mammd, aes(x=year, y=Julian, col=mammals))+geom_point()+
  stat_smooth(method="lm")+geom_hline(yintercept=268, linetype=2)+ylab("peak migration date(Julian)")+
  theme_classic()+ggtitle("mammal diet")

ggplot(gos_sp_mammc, aes(x=year, y=log(yr_tot), col=mammals))+geom_point()+
  stat_smooth(method="lm")+theme_classic()+ylab ("total counts (log)")+ggtitle("mammal diet")
