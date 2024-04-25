# -toy model for calculating the effect of species' phenology and abundance on aggregate measures of phenological shift

#scenario: community-level delayed migration

phen1_com=rep(212,100)
phen2_com=rep(215,100)

delta_P=phen2_com-phen1_com
set.seed(10)
hist(rpois(100,delta_P), xlab="shift (days)", main="aggregate-level shift")

mean(delta_P) #3 days delay

#driver scenario 1#####
# species A (common species) delays migration

#community at T1

phen1_spA=rep(212, 100)
phen1_spB=rep(209, 100)
phen1_spC=rep(209, 100)
phen1_spD=rep(209, 100)

#community at T2

phen2_spA=rep(215, 100)
phen2_spB=rep(209, 100)
phen2_spC=rep(209, 100)
phen2_spD=rep(209, 100)

#abundances of species (for weighted shifts calculation)

relA_spA=50
relA_spB=50
relA_spC=50
relA_spD=50

#mean taxon-level shift

phen_spA=rpois(100,phen2_spA-phen1_spA)
phen_spB=rpois(100,phen2_spB-phen1_spB)
phen_spC=rpois(100,phen2_spC-phen1_spC)
phen_spD=rpois(100,phen2_spD-phen1_spD)

phen_sp=c(phen_spA+phen_spB+ phen_spC+phen_spD)

mean(phen_sp)

rand_prob_phen_sp= rpois(100, prob_phen_sp)

mean(rand_prob_phen_sp)
hist(rand_prob_phen_sp)

#relative abundance-weighted shift

phen_spA_ra=(phen2_spA-phen1_spA)*relA_spA
phen_spB_ra=(phen2_spB-phen1_spB)*relA_spB
phen_spC_ra=(phen2_spC-phen1_spC)*relA_spC

phen_sp_ra=c(phen_spA_ra+ phen_spB_ra+ phen_spC_ra)

prob_phen_sp_ra= mean(phen_sp_ra)/mean(delta_P) #80
set.seed(10)
rand_prob_phen_sp_ra= rpois(100, prob_phen_sp_ra)

mean(rand_prob_phen_sp_ra)
mean(prob_phen_sp_ra)

#compare the 2 metrics


#driver scenario 2#####
# species A (late migrant) increases abundance

#community at T1

phen1_spA=rep(212, 100)
phen1_spB=rep(209, 100)
phen1_spC=rep(209, 100)

#community at T2

phen2_spA=rep(212, 100)
phen2_spB=rep(209, 100)
phen2_spC=rep(209, 100)

#relative abundances of species (for weighted shifts calculation)

relA1_spA=0.5
relA1_spB=0.3
relA1_spC=0.2

relA2_spA=0.8
relA2_spB=0.1
relA2_spC=0.1

#mean taxon-level shift

phen_spA=(phen2_spA-phen1_spA)
phen_spB=(phen2_spB-phen1_spB)
phen_spC=(phen2_spC-phen1_spC)

phen_sp=c(phen_spA+ phen_spB+ phen_spC)

prob_phen_sp= phen_sp/mean(delta_P) 
set.seed(10)
rand_prob_phen_sp= rpois(100, prob_phen_sp)

mean(rand_prob_phen_sp)
mean(prob_phen_sp)

#relative abundance-weighted shift

phen_spA_ra=(phen2_spA*relA2_spA)-(phen1_spA*relA1_spA)
phen_spB_ra=(phen2_spB*relA2_spB)-(phen1_spB*relA1_spB)
phen_spC_ra=(phen2_spC*relA2_spC)-(phen1_spC*relA1_spC)

phen_sp_ra=c(phen_spA_ra+ phen_spB_ra+ phen_spC_ra)

prob_phen_sp_ra= mean(phen_sp_ra)/mean(delta_P) #80#
set.seed(10)
rand_prob_phen_sp_ra= rpois(100, prob_phen_sp_ra)

mean(rand_prob_phen_sp_ra)
mean(prob_phen_sp_ra)

#compare the 2 metrics
par(mfrow=c(1,2))
hist(rand_prob_phen_sp, xlab="probability", main="mean species phenology shift")
hist(rand_prob_phen_sp_ra, xlab="probability", main="abundance-weighted species phenology shift")


