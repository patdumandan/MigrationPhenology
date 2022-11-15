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
