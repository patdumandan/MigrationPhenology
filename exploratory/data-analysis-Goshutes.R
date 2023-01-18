require(dplyr)
require(tidyr)
require(lubridate)
require(ggplot2)
require(vegan)

#I. splines model####
#data manipulation####
GosMts=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/GosMts.csv")
GosMts[, 3:27][is.na(GosMts[, 3:27])] <- 0
GosMts$Julian=as.POSIXlt(GosMts$Date)$yday 
GosMts=GosMts%>%select(-TOTAL, -RS, -MK, -BV, -BE, -UA, -UB, -UF, -UR, -UE)%>%
  mutate(year=year(Date), month=month(Date), date=day(Date), total=rowSums(.[3:18]))%>%
  filter(!year>2018,!year<1984, !Julian>308)

#remove unknowns and uncommon species (BAEA??)
#select only 1983-2018 (Aug 15-Nov 5)

#species-level 
gos=GosMts%>%select(-total)%>%
  pivot_longer(cols=3:18, names_to = "Species", values_to="Count")

#annual totals per species
gos_sp_tot=gos%>%
  group_by(year, Species)%>%
  summarise(sp_tot=sum(Count))

#determine uncommon species (based on non-detects annually)
gos_sp_min=gos_sp_tot%>%group_by(Species)%>%filter(sp_tot==min(sp_tot))

#daily totals per species
gos_sp_day=gos%>%
  group_by(Species,year, Julian)%>%
  summarise(total=sum(Count))%>%
  arrange(year)


#10% passage date: species
gos_10pd_sp=gos_sp_day%>%group_by(Species, year)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                  cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                  PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(Species, year)%>%filter(PD_10==1)%>%slice_head(n=1)  #select rows where cum_pt is >50%


#10% passage date: site
gos_10pd_site=GosMts%>%group_by(year)%>%mutate(cum_tot=cumsum(total), #get cumulative sum for each year
                                               cum_pt=cum_tot/sum(total),#calculate % of daily count/annual count
                                               PD_10=if_else(cum_pt>=0.1,"1", "0"))%>% #flag rows with cumulative % >50%
  group_by(year)%>%filter(PD_10==1)%>%slice_head(n=1)%>%#select rows where cum_pt is >10%
  select(year, Julian, cum_pt)

#create cubic cyclic splines in stan (like GAM by mgcv but manual and bayesian)

#LOAD PACKAGES
require(splines)
require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#MANIPULATE DATA
gos_10pd_site$logJD=log(gos_10pd_site$Julian)
gos_10pd_site$years= (gos_10pd_site$year - mean(gos_10pd_site$year))/(2 *sd(gos_10pd_site$year)) 

#FOR CREATING SPLINES
B1 = t(bs(gos_10pd_site$years, df=NULL, knots=NULL, degree=3, intercept = FALSE)) # creating the B-splines, degree=3 for cubic spline
num_basis1 = nrow(B1)

#FOR STAN MODEL
dat_list2=list(N=length(gos_10pd_site$years),
               pd_10=gos_10pd_site$logJD,
               years=gos_10pd_site$years,
               B1=B1,
               num_basis1=num_basis1)

mod2=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/exploratory/splines_10pd.stan", 
          data=dat_list2,
          iter=300, chains=3)   

#process results####
post=rstan::extract(mod2)$Y_hat
post2=rstan::extract(mod2)$y_pred

#obtain mean for first and last 5 years of data
post=post%>%summarise_all(.funs = c(mean="mean"))

#plotting splines####
ff<-extract(mod2)

Y_hat_med <- array(NA, length(gos_10pd_site$year)) #median estimate
Y_hat_ub <- array(NA, length(gos_10pd_site$year)) #upper boundary
Y_hat_lb <- array(NA, length(gos_10pd_site$year)) #lower boundary

for (i in 1:length(gos_10pd_site$year)) {
  Y_hat_med[i] <- median(exp(ff$Y_hat[,i]));
  Y_hat_lb[i] <- quantile(exp(ff$Y_hat[,i]),probs = 0.025)
  Y_hat_ub[i] <- quantile(exp(ff$Y_hat[,i]),probs = 0.975)
}

plot(gos_10pd_site$Julian~gos_10pd_site$year, pch=19, ylab="10% passage date", xlab="year", main="Goshute Mts. phenological trend") #plot raw data
lines(smooth.spline(gos_10pd_site$year, Y_hat_med), col="blue")
lines(smooth.spline(gos_10pd_site$year, Y_hat_lb), lty=2, col="red") #0.025
lines(smooth.spline(gos_10pd_site$year, Y_hat_ub), lty=2, col="red") #0.975

exp_post2=as.data.frame(exp_post2)
exp_post=as.data.frame(exp(post))

diff1=exp_post[,1:10]
diff2=exp_post[,26:35]

diff1$rowmean_1=apply(diff1[,1:10],1,mean)
diff2$rowmean_2=apply(diff2[,1:10],1,mean)

diffs=cbind(diff1, diff2)%>%select(rowmean_1, rowmean_2)%>%mutate(med_diff=rowmean_2-rowmean_1)

quantile(diffs$med_diff, probs = c(0.025, 0.5, 0.975), na.rm = FALSE) #delay by 0.32 days [-1.51, 2.11]
hist(diffs$med_diff, xlab="shift in timing(d)", main="community phenological shift")
abline(v=0, lwd=3,lty=2)
length(which(diffs$med_diff>0))/length(diffs$med_diff) #64% probability of delayed migration

#II. breakpoint model####
#FOR STAN MODEL
dat_list2=list(N=length(gos_10pd_site$years),
               pd_10=gos_10pd_site$logJD,
               years=gos_10pd_site$years)

mod_bkpt=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/exploratory/period_bkpt.stan", 
          data=dat_list2,
          iter=300, chains=3)  
