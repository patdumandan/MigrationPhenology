require(dplyr)
require(tidyr)
require(lubridate)
require(ggplot2)

#data manipulation####
hms=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/HMSDay.csv")

#combine baea and goea data
HMS=hms%>%
  mutate(BAEA=rowSums(.[9:11]), GOEA=rowSums(.[21:23]))%>%
  select(-OBSERVERS, -GOEA.I, -GOEA.A, -GOEA.U, -BAEA.I, -BAEA.U, -BAEA.A, -RAPTORS, -GYRF, -MIKI, -STKI,-SWHA,
         -UNID.ACCIPITER,-UNID.BUTEO, -UNID.EAGLE, -UNID.FALCON, -UNID.RAPTOR)%>%
  mutate(yr_tot=rowSums(.[6:21], na.rm=T), date=make_date(YR, MO, DAY))%>%filter(!(YR<1983))

HMS$Julian=as.POSIXlt(HMS$date)$yday
HMS[, 6:21][is.na(HMS[, 6:21])] <- 0

#select only 1983-2018 (Aug 15-dEC 15: 226-348)

HMS=HMS%>%filter(!Julian<226, !Julian>348)

#species-level 
hms_sp=HMS%>%
  pivot_longer(cols=6:21, names_to = "Species", values_to="Count")

#annual totals per species
hms_sp_tot=hms_sp%>%
  group_by(YR, Species)%>%
  summarise(sp_tot=sum(Count))

#determine uncommon species (based on non-detects annually)
hms_sp_min=hms_sp_tot%>%group_by(Species)%>%filter(sp_tot==min(sp_tot))

#daily totals per species
hms_sp_day=hms_sp%>%
  group_by(Species,YR, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(YR)

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

#create cubic cyclic splines in stan (like GAM by mgcv but manual and bayesian)

#LOAD PACKAGES
require(splines)
require(rstan)
require(brms)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#MANIPULATE DATA
hms_10pd_site$logJD=log(hms_10pd_site$Julian)
hms_10pd_site$years= (hms_10pd_site$YR - mean(hms_10pd_site$YR))/(2 *sd(hms_10pd_site$YR)) 


#FOR CREATING SPLINES
B1 = t(bs(hms_10pd_site$years, df=NULL, knots=NULL, degree=3, intercept = FALSE)) # creating the B-splines, degree=3 for cubic spline
num_basis1 = nrow(B1)

#FOR STAN MODEL
dat_list2=list(N=length(hms_10pd_site$years),
               pd_10=hms_10pd_site$Julian,
               years=hms_10pd_site$years,
               B1=B1,
               num_basis1=num_basis1)

mod_hms=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/exploratory/splines_10pd.stan", 
          data=dat_list2,
          iter=15000, chains=3)   

mod_hms_pois=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/exploratory/splines_poisson.stan", 
             data=dat_list2,
             iter=500, chains=3)   

#process results####
post=rstan::extract(mod_hms_pois)$Y_hat

#obtain mean for first and last 5 years of data
post=post%>%summarise_all(.funs = c(mean="mean"))

#plotting splines####
ff<-extract(mod_hms_pois)

Y_hat_med <- array(NA, length(hms_10pd_site$YR)) #median estimate
Y_hat_ub <- array(NA, length(hms_10pd_site$YR)) #upper boundary
Y_hat_lb <- array(NA, length(hms_10pd_site$YR)) #lower boundary

for (i in 1:length(hms_10pd_site$YR)) {
  Y_hat_med[i] <- median(ff$Y_hat[,i]);
  Y_hat_lb[i] <- quantile(ff$Y_hat[,i],probs = 0.025)
  Y_hat_ub[i] <- quantile(ff$Y_hat[,i],probs = 0.975)
}

plot(hms_10pd_site$Julian~hms_10pd_site$YR, pch=19, ylab="10% passage date", xlab="year", main="Hawk Mt. phenological trend",
     ylim=c(245, 276)) #plot raw data
lines(smooth.spline(hms_10pd_site$YR, Y_hat_med), col="blue")
lines(smooth.spline(hms_10pd_site$YR, Y_hat_lb), lty=2, col="red") #0.025
lines(smooth.spline(hms_10pd_site$YR, Y_hat_ub), lty=2, col="red") #0.975

exp_post2=as.data.frame(exp_post2)
exp_post=as.data.frame(exp(post))

diff1=exp_post[,1:10]
diff2=exp_post[,27:36]

diff1$rowmean_1=apply(diff1[,1:10],1,mean)
diff2$rowmean_2=apply(diff2[,1:10],1,mean)

diffs=cbind(diff1, diff2)%>%select(rowmean_1, rowmean_2)%>%mutate(med_diff=rowmean_2-rowmean_1)

quantile(diffs$med_diff, probs = c(0.025, 0.5, 0.975), na.rm = FALSE) #delay by 0.32 days [-1.51, 2.11]
hist(diffs$med_diff, xlab="shift in timing(d)", main="community phenological shift")
abline(v=0, lwd=3,lty=2)
length(which(diffs$med_diff>0))/length(diffs$med_diff) #64% probability of delayed migration
