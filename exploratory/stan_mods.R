require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

dat_list1=list(N=length(t1$Species),
               Nsp=length(unique(t1$Species)),
               div_index=t1$div_index,
               div_peak=t1$div_peak,
               years=t1$years,
               sp_time=t1$sp_peak,
               sp_count=t1$total, 
               spcode=t1$spcode)

#simgple GLM####
mod1=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/MigrationPhenology/diversity_simple.stan", 
           data=dat_list1,
           iter=300, chains=3)                 
                 

posterior=extract(mod1)

plot(t1$div_index ~ t1$years, pch = 20)
abline( mean(posterior$alpha), mean(posterior$b0), col = 6, lw = 2)
for (i in 1:300) {
  abline(posterior$alpha[i], posterior$b0[i], col = "gray", lty = 1)
}

#breakpoint model: index values#####
mod2=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/MigrationPhenology/diversity_bkpt.stan", 
          data=dat_list1,
          iter=300, chains=3)                 

print(mod2, pars=c("bkpt", "slope1", "slope2"))

posterior=extract(mod2)$y_mu

mean_post=apply(posterior, 2, mean)
mean_post1=unique(mean_post)
plot(t1$div_index~t1$year, pch=20, xaxt="n")

plot(mean_post1~c(0:35),xaxt="n", type="l", xlab="year", ylim=c(1.75,2.4), col="red", xaxt="n",ylab="diversity index")
points(unique(t1$div_index)~c(0:35), pch=20)
axis(1,c(1983:2018),at=c(0:35))
abline(v=14, lty=2)

#breakpoint model: peak diversity datess#####
mod3=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/MigrationPhenology/diversitydate_bkpt.stan", 
          data=dat_list1,
          iter=300, chains=3)                 

print(mod3, pars=c("bkpt", "slope1", "slope2"))

posterior3=extract(mod3)$y_mu

mean_post3=apply(posterior3, 2, mean)
mean_post31=unique(mean_post3)
plot(t1$div_peak~t1$year, pch=20, xaxt="n")

plot(mean_post31~c(0:35),xaxt="n", type="l", xlab="year", col="red",ylab="peak diversity date")
points(g1$sd~g1$year, pch=20, xaxt="n")
axis(1,c(1983:2018),at=c(0:35))
abline(v=14, lty=2)

dat_list2=list(N=length(div_dat2$year),
               div_index=div_dat2$sd,
               years=div_dat2$years,
               bird_time=div_dat2$bird,
               bird_count=div_dat2$bird_ct
               )

mod4=stan(file="D:/Dropbox (UFL)/PhD-stuff/MigrationPhenology/MigrationPhenology/bird-effect.stan", 
          data=dat_list2,
          iter=300, chains=3)                 

print(mod4, pars=c("b0", "time_eff", "count_eff"))


posterior=extract(mod4)$y_mu

summary(glm(sd~bird_ct+bird, data=div_dat2))

