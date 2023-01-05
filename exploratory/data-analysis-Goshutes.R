require(dplyr)
require(tidyr)
require(lubridate)
require(ggplot2)
require(vegan)

#data manipulation####
GosMts=read.csv("D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/data/GosMts.csv")
GosMts[, 3:27][is.na(GosMts[, 3:27])] <- 0
GosMts$Julian=as.POSIXlt(GosMts$Date)$yday 
GosMts=GosMts%>%select(-TOTAL, -RS, -MK, -BV, -BE, -UA, -UB, -UF, -UR, -UE)%>%
  mutate(year=year(Date), month=month(Date), date=day(Date), total=rowSums(.[3:18]))%>%
  filter(!year>2018,!year<1984, !Julian>308)

#remove unknowns and uncommon species (BAEA??)
#select only 1983-2018 (Aug 15-Nov 5)

#functional/species-level patterns
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

#MODEL RESULTS PROCESSING####
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

exp_post=as.data.frame(exp(post))
diff1=exp_post[,1:10]
diff2=exp_post[,26:35]

diff1$rowmedian_1=apply(diff1[,1:10],1,median)
diff2$rowmedian_2=apply(diff2[,1:10],1,median)

diffs=cbind(diff1, diff2)%>%select(rowmedian_1, rowmedian_2)%>%mutate(med_diff=rowmedian_2-rowmedian_1)

quantile(diffs$med_diff, probs = c(0.025, 0.5, 0.975), na.rm = FALSE)



