require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

gos_dat=list(N=length(gostotmax$year),
              years=gostotmax$year,
              peak_date=gostotmax$Julian)
mod_gos=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/MigrationPhenology/bird_timing.stan", 
              data=gos_dat,
              iter=300, chains=3)                 
print(mod_gos, pars=c("bkpt", "slope1", "slope2", "alpha"))

posterior1=extract(mod_gos)$y_mu

mean_post1=apply(posterior1, 2, mean)

matplot(t(posterior1), col="grey", type="l", xaxt="n", ylab="peak migration date(Julian)")
lines(mean_post1~c(1:37), col="blue") 
points(gostotmax$Julian~c(1:37), pch=20)
axis(1,c(1983:2019),at=c(1:37))
mtext("Goshutes", side=3)

#breakpoint model: peak migration dates for bird eaters#####
mod_bird=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/MigrationPhenology/bird_timing.stan", 
              data=bird_dat,
              iter=300, chains=3)                 

bird_dat=list(N=length(birdmax$year),
              years=birdmax$year,
              peak_date=birdmax$Julian)

print(mod_bird, pars=c("bkpt", "slope1", "slope2", "alpha"))

posterior3=extract(mod_bird)$y_mu

mean_post3=apply(posterior3, 2, mean)
mean_post31=unique(mean_post3)

matplot(t(posterior3), col="grey", type="l", xaxt="n", ylab="peak migration date(Julian)")
lines(mean_post31~c(1:37), col="blue") 
points(birdmax$Julian~c(1:37), pch=20)
axis(1,c(1983:2019),at=c(1:37))
mtext("bird eaters", side=3)

#breakpoint model: peak migration dates for mammal eaters#####
mod_mam=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/MigrationPhenology/bird_timing.stan", 
             data=mam_dat,
             iter=300, chains=3)                 

mam_dat=list(N=length(mammax$year),
             years=mammax$year,
             peak_date=mammax$Julian)

print(mod_mam, pars=c("bkpt", "slope1", "slope2", "alpha"))

posterior4=extract(mod_mam)$y_mu

mean_post4=apply(posterior4, 2, mean)
mean_post41=unique(mean_post4)

lines(mean_post41~c(1:37), col="blue") 
matplot(t(posterior4), col="grey", type="l", xaxt="n", ylab="peak migration date(Julian)")
points(mammax$Julian~c(1:37), pch=20)
axis(1,c(1983:2019),at=c(1:37))
mtext("mammal eaters", side=3)

#count_trend####
gos_sp$years = (gos_sp$year - mean(gos_sp$year))/(2 *sd(gos_sp$year)) 

count_dat=list(N=length(gos_sp$year),
               Nsp=length(unique(gos_sp$Species)),
                years=gos_sp$years,
               count=gos_sp$yr_tot,
               spcode=gos_sp$spcode,
               zeros2=zeros2)

zeros2=rep(0,2)

mod_count=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/exploratory/count_trend_nbmod.stan", 
             data=count_dat,
             iter=300, chains=3)                 

print(mod_count, pars=c("beta[2]", "alpha"))

#max date trend####

gos_sp_day_max$years = (gos_sp_day_max$year - mean(gos_sp_day_max$year))/(2 *sd(gos_sp_day_max$year)) 
gos_sp_day_max$dates = (gos_sp_day_max$Julian - mean(gos_sp_day_max$Julian))/(2 *sd(gos_sp_day_max$Julian)) 

time_dat=list(N=length(gos_sp_day_max$year),
               Nsp=length(unique(gos_sp_day_max$Species)),
               years=gos_sp_day_max$years,
               count=gos_sp_day_max$total,
               spcode=gos_sp_day_max$spcode,
               zeros3=zeros3,
               Julian=gos_sp_day_max$dates)

zeros3=rep(0,3)

mod_time=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/exploratory/maxdate_trend_nbmod.stan", 
               data=time_dat,
               iter=300, chains=3)                 

print(mod_time, pars=c("beta[3]", "alpha"))

#diversity####
g1$years = (g1$year - mean(g1$year))/(2 *sd(g1$year)) 
g1$dates = (g1$Julian - mean(g1$Julian))/(2 *sd(g1$Julian)) 

dat_list1=list(N=length(g1$years),
               div_index=g1$specdiv,
               Julian=g1$dates,
               years=g1$years)

#simgple GLM####
mod1=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/exploratory/diversity_simple.stan", 
          data=dat_list1,
          iter=300, chains=3)   

print(mod1)

#splines model####
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
B1 = t(bs(gos_10pd_site$year, df=NULL, knots=NULL, degree=3, intercept = FALSE)) # creating the B-splines, degree=3 for cubic spline
num_basis1 = nrow(B1)

#FOR STAN MODEL
dat_list2=list(N=length(gos_10pd_site$year),
               pd_10=gos_10pd_site$logJD,
               years=gos_10pd_site$years,
               B1=B1,
               num_basis1=num_basis1)

mod2=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/exploratory/splines_10pd.stan", 
          data=dat_list2,
          iter=3000, chains=3)   

#MODEL RESULTS PROCESSING
post=rstan::extract(mod2)$Y_hat

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

plot(gos_10pd_site$Julian~gos_10pd_site$year, pch=19, ylab="10% passage date", xlab="year") #plot raw data
lines(smooth.spline(gos_10pd_site$year, Y_hat_med), col="blue")
lines(smooth.spline(gos_10pd_site$year, Y_hat_lb), lty=2, col="red") #0.025
lines(smooth.spline(gos_10pd_site$year, Y_hat_ub), lty=2, col="red") #0.975

