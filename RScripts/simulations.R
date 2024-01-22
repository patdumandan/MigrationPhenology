#simulation 1####
# assume that sp 1 is a late migrant, and is delaying migration
# assume that sp 2, 3, and 4 are not changing migration timing and abundances
# assume that all species have similar abundances over time

set.seed(10)
nsp=4
nrep=20

sp_dat=data.frame(
  year=rep(1980:2000, each=nsp, length.out=nrep*nsp),
  species=as.character(rep(1:4, length.out=nrep*nsp)))

sp_dat1=sp_dat%>%
  mutate(
    period=case_when(year%in%c(1980:1984) ~"T1",
                     year%in%c(1995:1999) ~ "T2",
                     year%in%c(1985:1994)~ "TN"),
    pd_10= case_when(species=="1" & period== "T1" ~ 222,
                     species=="1" & period== "T2" ~ 230,
                     species%in%c("2", "3", "4") & period== "T1"  ~ 198,
                     species%in%c("2", "3", "4")& period== "T2" ~ 198,
                     species%in%c("1", "2", "3", "4") & period=="TN"~ rpois(year%in%c(1985:1994), 202)))

#summarize unit-level values
sp_sum1=sp_dat1%>%group_by(species, period)%>%
  summarise(u_pd=mean(pd_10))%>%filter(!period=="TN")%>%
  mutate(abundance=50, 
         rel_abund=abundance/sum(abundance))

s1=sp_dat1%>%filter(year%in%c(1980:1984))
s2=sp_dat1%>%filter(year%in%c(1995:1999))

mean(s1$pd_10)
mean(s2$pd_10)

#create community data

comm_dat1=sp_dat1%>%group_by(year, period)%>%summarise(pd_10=mean(pd_10))

c1=comm_dat1%>%filter(year%in%c(1980:1984))
c2=comm_dat1%>%filter(year%in%c(1995:1999))

mean(c1$pd_10)
mean(c2$pd_10)

phen_comm=mean(c2$pd_10)-mean(c1$pd_10)

#NOT from sp_table
phen1_com=rep(212,100)
phen2_com=rep(215,100)

delta_P=phen2_com-phen1_com
set.seed(10)
hist(rpois(100,delta_P), xlab="shift (days)", main="aggregate-level shift")

mean(delta_P) #3 days delay
#visualize data####

ggplot(comm_dat1, aes(x=year, y=pd_10))+geom_line()+
  geom_hline(yintercept=mean(comm_dat1$pd_10), lty=2)+
  annotate("rect",xmin=1980, xmax=1984, ymin=min(comm_dat1$pd_10), ymax=max(comm_dat1$pd_10), 
           color="red", alpha=0.1, fill="red")+theme_classic()+
  annotate("rect",xmin=1999, xmax=1995, ymin=min(comm_dat1$pd_10), ymax=max(comm_dat1$pd_10), 
           color="red", alpha=0.1, fill="red")+theme_classic()+
  ggtitle("aggregate-level")+ylab("10% passage date")

ggplot(sp_dat1, aes(x=year, y=pd_10))+geom_line()+
  facet_wrap(~species)+ geom_hline(yintercept=mean(comm_dat1$pd_10), lty=2)+
  annotate("rect",xmin=1980, xmax=1984, ymin=min(sp_dat1$pd_10), ymax=240, 
           color="red", alpha=0.1, fill="red")+theme_classic()+
  annotate("rect",xmin=1999, xmax=1995, ymin=min(sp_dat1$pd_10), ymax=240, 
           color="red", alpha=0.1, fill="red")+theme_classic()+
  ggtitle("unit-level")+ylab("10% passage date")

# metrics for phenological shift#### 

##mean species-level shift#### 

#calculate phenological shift (species level)
phen_sp=mean(s2$pd_10)-mean(s1$pd_10)

#calculate phenological shift (community level)
phen_comm=mean(c2$pd_10)-mean(c1$pd_10)

##relative abundance-weighted phenological shift####

sp11=sp_sum1%>%
  pivot_wider(names_from = period, values_from=u_pd, names_sep = ".")%>%
  mutate(shift= T2-T1,
         wgt_shift=rel_abund*shift)

#contribution to community-level shift####

#relative abundance weighted

phen_sp_cont=sum(sp11$wgt_shift)/phen_comm #2 or 200% delay in species-level phenology is contributing to the community-level shift

phen_
mean(sp11$wgt_shift)
mean(sp11$shift)

#simulation 2####
# assume that sp 1 is a late migrant
# assume that sp 2,3, and 4 are not changing migration timing and abundances
# assume that sp 1 is the most abundant

sp_dat2=sp_dat%>%
  mutate(
    period=case_when(year%in%c(1980:1984) ~"T1",
                     year%in%c(1995:1999) ~ "T2",
                     year%in%c(1985:1994)~ "TN"),
    pd_10= case_when(species=="1" & period== "T1" ~ 222,
                     species=="1" & period== "T2" ~ 222,
                     species%in%c("2", "3","4") & period== "T1"  ~ 198,
                     species%in%c("2", "3", "4")& period== "T2" ~ 198,
                     species%in%c("1", "2", "3", "4") & period=="TN"~ rpois(year%in%c(1985:1994), 202)),
    rel_abund= case_when(species=="1"~ 0.31,
                         species%in%c("2", "3","4")~ 0.23))

#summarize unit-level values
s1=sp_dat2%>%filter(year%in%c(1980:1984))
s2=sp_dat2%>%filter(year%in%c(1995:1999))

mean(s1$pd_10)
mean(s2$pd_10)

#create community data

comm_dat2=sp_dat2%>%group_by(year, period)%>%summarise(pd_10=mean(pd_10))

c1=comm_dat2%>%filter(year%in%c(1980:1984))
c2=comm_dat2%>%filter(year%in%c(1995:1999))

mean(c1$pd_10)
mean(c2$pd_10)

#calculate phenological shift
phen_comm=mean(c2$pd_10)-mean(c1$pd_10)

#visualize data####

ggplot(comm_dat2, aes(x=year, y=pd_10))+geom_line()+
  geom_hline(yintercept=mean(comm_dat2$pd_10), lty=2)+
  annotate("rect",xmin=1980, xmax=1984, ymin=min(comm_dat2$pd_10), ymax=max(comm_dat2$pd_10), 
           color="red", alpha=0.1, fill="red")+theme_classic()+
  annotate("rect",xmin=1999, xmax=1995, ymin=min(comm_dat2$pd_10), ymax=max(comm_dat2$pd_10), 
           color="red", alpha=0.1, fill="red")+theme_classic()+
  ggtitle("aggregate-level")+ylab("10% passage date")


ggplot(sp_dat2, aes(x=year, y=pd_10))+geom_line()+
  facet_wrap(~species)+
  annotate("rect",xmin=1980, xmax=1984, ymin=min(sp_dat2$pd_10), ymax=max(sp_dat2$pd_10), 
           color="red", alpha=0.1, fill="red")+theme_classic()+
  annotate("rect",xmin=1999, xmax=1995, ymin=min(sp_dat2$pd_10), ymax=max(sp_dat2$pd_10), 
           color="red", alpha=0.1, fill="red")+theme_classic()+
  ggtitle("unit-level")+ylab("10% passage date")

# metrics for phenological shift#### 

##mean species-level shift#### 

#calculate phenological shift (species level)
phen_sp=mean(s2$pd_10)-mean(s1$pd_10)

#calculate phenological shift (community level)
phen_comm=mean(c2$pd_10)-mean(c1$pd_10)
#NOT from sp_table
phen1_com=rep(212,100)
phen2_com=rep(215,100)

delta_P=phen2_com-phen1_com
set.seed(10)
hist(rpois(100,delta_P), xlab="shift (days)", main="aggregate-level shift")

mean(delta_P) #3 days delay

##relative abundance-weighted phenological shift####

sp_dat21=sp_dat2%>%select(species, rel_abund)%>%
  unique(n=1)

sp_sum2=sp_dat2%>%group_by(species, period)%>%
  summarise(u_pd=mean(pd_10))%>%filter(!period=="TN")

sp22=sp_sum2%>%
  pivot_wider(names_from = period, values_from=u_pd, names_sep = ".")%>%
  mutate(shift=T2-T1)%>%
  left_join(sp_dat21, by="species")%>%
  mutate(wgt_shift=rel_abund*shift)

#contribution to community-level shift####

#relative abundance weighted

phen_sp_cont=sum(sp22$wgt_shift)/phen_comm 
phen_sp_cont2=sum(sp22$wgt_shift)/mean(delta_P)

mean(sp22$wgt_shift)
mean(sp22$shift)
