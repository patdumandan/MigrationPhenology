# -toy model for calculating the effect of species' phenology and abundance on aggregate measures

#scenario: community-level delayed migration


nrep=10 #no. of replicates
nsp= 3 #no. of species
dchange= 3 #phenological shift (no of days delay)


#create data frame

# iterations are for when different time-series data lengths are tested

# -phen_shift= days of change 

set.seed(1234)

for (j in 1:nrep){
  
  # species data 
  
  sp_dat=data.frame(
    TS_dat=rep(1:nrep, 1),
    species=as.character(rep(1:3, each=nrep, length.out=nrep*nsp)),
    sp_shift=rep(rnorm(n=nsp*nrep, mean=dchange, sd=1))
  )
  
  #community data
  
  comm_dat=data.frame(
    TS_dat=rep(1:nrep, 1),
    comm_shift=rep(rnorm(n=nrep, mean=dchange, sd=1))
    
  )
}

#add relative abundances (assumed no change)

sp_dat=sp_dat%>%
  mutate(rel_abund=case_when(species=="1" ~ 0.8,
                             species=="2" ~ 0.1,
                             species=="3" ~ 0.1))
               
all_dat=right_join(comm_dat, sp_dat, by="TS_dat")


T1=all_dat%>%filter(TS_dat==1)
