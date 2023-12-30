set.seed(1234)

nrep=20
nsp=3

#create species data frame
sp_dat=data.frame(
  year=rep(1980:2000, each=nsp, length.out=nrep*nsp),
  species=as.character(rep(1:3, length.out=nrep*nsp)))

sp_dat1=sp_dat%>%
  mutate(
    period=case_when(year%in%c(1980:1984) ~"T1",
                     year%in%c(1995:1999) ~ "T2",
                     year%in%c(1985:1994)~ "TN"),
    pd_10= case_when(species=="1" & period== "T1" ~ rpois(length(year%in%c(1980:1984)), 222),
                     species=="1" & period== "T2" ~ rpois(length(year%in%c(1995:1999)), 225),
                     species%in%c("2", "3") & period== "T1"  ~ 198,
                     species%in%c("2", "3")& period== "T2" ~ 198,
                     species%in%c("1", "2", "3") & period=="TN"~ rpois(year%in%c(1985:1994), 202)))

#summarize unit-level values
sp_sum1=sp_dat1%>%group_by(species, period)%>%summarise(u_pd=mean(pd_10))%>%filter(!period=="TN")

s1=sp_dat1%>%filter(year%in%c(1980:1984))
s2=sp_dat1%>%filter(year%in%c(1995:1999))

mean(s1$pd_10)
mean(s2$pd_10)

#calculate phenological shift 
phen_sp=mean(s2$pd_10)-mean(s1$pd_10)

#create community data

comm_dat=sp_dat1%>%group_by(year, period)%>%summarise(pd_10=mean(pd_10))

c1=comm_dat%>%filter(year%in%c(1980:1984))
c2=comm_dat%>%filter(year%in%c(1995:1999))

mean(c1$pd_10)
mean(c2$pd_10)

#calculate phenological shift
phen_comm=mean(c2$pd_10)-mean(c1$pd_10)

#plot

ggplot(comm_dat, aes(x=year, y=pd_10))+geom_line()+
  geom_hline(yintercept=mean(comm_dat$pd_10), lty=2)+
  annotate("rect",xmin=1980, xmax=1984, ymin=min(comm_dat$pd_10), ymax=max(comm_dat$pd_10), 
           color="red", alpha=0.1, fill="red")+theme_classic()+
  annotate("rect",xmin=1999, xmax=1995, ymin=min(comm_dat$pd_10), ymax=max(comm_dat$pd_10), 
           color="red", alpha=0.1, fill="red")+theme_classic()+
  ggtitle("aggregate-level")+ylab("10% passage date")


ggplot(sp_dat1, aes(x=year, y=pd_10))+geom_line()+
  facet_wrap(~species)+
  annotate("rect",xmin=1980, xmax=1984, ymin=min(sp_dat$pd_10), ymax=max(sp_dat$pd_10), 
           color="red", alpha=0.1, fill="red")+theme_classic()+
  annotate("rect",xmin=1999, xmax=1995, ymin=min(sp_dat$pd_10), ymax=max(sp_dat$pd_10), 
           color="red", alpha=0.1, fill="red")+theme_classic()+
  ggtitle("unit-level")+ylab("10% passage date")

