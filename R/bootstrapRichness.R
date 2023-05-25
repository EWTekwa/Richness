bootstrapRichness <- function(Community, numBoot = 100){

# bootRichnessEsts.m
# Ed Tekwa Feb 8, 2022
# bootRichnessEsts.R
# Matt Whalen rewrote original matlab script for R - Mar 27, 2022
# function several estimates of richness:
#  - raw richness
#  - newly proposed method
#  - Chao1
#  - Chao2
#  - Abundance-based coverage estimator (ACE)
#  - Jackknife abundance estimator
#  - Jackknife incidence estimator
# Point estimates are calculated using script "RichnessEsts.R"
# This script calculates bootstrapped confidence intervals
# the spatial Community data: rows = transects, columns = species,
#                                     values = individual counts
# Whalen update 9 April 2023
  # - remove Clustering calculations
  # - add estimators Gamma Poisson and edited Omega estimates
# Whalen update 24 May 2023
  # - fixing minor errors in preparation for paper release


# Get point estimates from original dataset
Community = Community[ ,colSums(Community) > 0 ] # remove empty species columns
# run function to calculate all point estimates - ignore correction terms states and approximated detection probabilities
pointests = RichnessEsts( Community )[[1]]


# define variance function with normalizatin of n rather than n-1
varn <- function(x) mean((x-mean(x))^2)


# store bootstrapped estimates for raw richness and each estimator
expectedRichness_raw   = rep(0,numBoot)  # raw richness
expectedChao1 = rep(0,numBoot)           # Chao1
expectedGP    = rep(0,numBoot)           # Gamma-Poisson
expectedChao2 = rep(0,numBoot)           # Chao2
expectedACE   = rep(0,numBoot)           # ACE
expectedJK_a  = rep(0,numBoot)           # Jackknife (abundance)
expectedJK_i  = rep(0,numBoot)           # Jackknife (incidence)
expectedOmega = rep(0,numBoot)           # approximate Omega richness
expectedOmega_taylor   = rep(0,numBoot)  # approximated Omega richness using Taylor expansion
expectedOmega_taylor_0 = rep(0,numBoot)  # approximated Omega richness using the first correction term



# write observational process model
Di <- expression( 1-(1-((1-exp(-nm))*P))^k )
# second-order derivatives
d2Di_dnm2 <- D(D( Di, "nm" ), "nm")
d2Di_dP2  <- D(D( Di, "P" ), "P")
d2Di_dnmP  <- D(D( Di, "nm" ), "P")


# get the number of sampling units (e.g., transects or quadrats)
numSamplingUnit = nrow(Community)


# bootstrapping - scramble sampling units and resample
for( resample in 1:numBoot ){
  #  scramble sampling units (e.g., transects) first
    sampleSet_init = Community[ sample(1:numSamplingUnit,numSamplingUnit,replace = T), ]

    # using unit-scrambled data, next scramble the species
    sampleSetRaw   = sampleSet_init[ ,sample(1:pointests$Richness_raw,pointests$Richness_raw,replace = T) ]
    sampleSetChao1 = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$Chao1),replace = T) ]
    sampleSetGP    = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$GP),replace = T) ]
    sampleSetChao2 = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$Chao2),replace = T) ]
    sampleSetACE   = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$ACE),replace = T) ]
    sampleSetJK_a  = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$JK_a),replace = T) ]
    sampleSetJK_i  = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$JK_i),replace = T) ]
    sampleSetOmega          = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$Omega),replace = T) ]
    sampleSetOmega_taylor   = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$Omega_taylor),replace = T) ]
    sampleSetOmega_taylor_0 = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$Omega_taylor_0),replace = T) ]

    # take out empty species columns
    sampleSetRaw   = sampleSetRaw[ , colSums(sampleSetRaw) > 0 ]
    sampleSetChao1 = sampleSetChao1[ , colSums(sampleSetChao1) > 0 ]
    sampleSetGP    = sampleSetGP[ , colSums(sampleSetGP) > 0 ]
    sampleSetChao2 = sampleSetChao2[ , colSums(sampleSetChao2) > 0 ]
    sampleSetACE   = sampleSetACE[ , colSums(sampleSetACE) > 0 ]
    sampleSetJK_a  = sampleSetJK_a[ , colSums(sampleSetJK_a) > 0 ]
    sampleSetJK_i  = sampleSetJK_i[ , colSums(sampleSetJK_i) > 0 ]
    sampleSetOmega          = sampleSetOmega[ , colSums(sampleSetOmega) > 0 ]
    sampleSetOmega_taylor   = sampleSetOmega_taylor[ , colSums(sampleSetOmega_taylor) > 0 ]
    sampleSetOmega_taylor_0 = sampleSetOmega_taylor_0[ , colSums(sampleSetOmega_taylor_0) > 0 ]

    # compute raw estimate
    Richness_raw_boot = sum(colSums(sampleSetRaw)>0)

    # compute Chao1 estimate
    sampleSetChao1 = ceiling(sampleSetChao1) # convert data to discrete count
    raw = sum(colSums(sampleSetChao1) > 0)
    f1 = sum(colSums(sampleSetChao1) == 1) # number of singleton species
    f2 = sum(colSums(sampleSetChao1) == 2) # number of doubleton species
    Chao1_boot = raw+f1*(f1-1)/(2*(f2+1)) # Chao1 richness estimator

    # compute Gamma-Poisson mixture estimate (Chiu 2023, peerJ)
    f1 = sum(colSums(sampleSetGP) == 1) # number of singleton species
    f2 = sum(colSums(sampleSetGP) == 2) # number of doubleton species
    f3 = sum(colSums(sampleSetGP)==3) # number of tripleton species
    if( f3 == 0 ){
      f3c = 1
    } else {
      f3c = f3
    }
    if( f1 == 0 ){
      f1c = 1
    } else {
      f1c = f1
    }
    A = 2-(2*f2^2)/(3*f1c*f3c)
    if( f2 > 0 ){
      f0Chao1 = f1c^2/(2*f2)
    } else {
      f0Chao1 = f1c*(f1c-1)/2
    }
    if( A < 1 ){
      GP_boot = Richness_raw_boot + f0Chao1*max(c(.5, A))
    } else {
      GP_boot = Richness_raw_boot + f0Chao1
    }

    # compute Chao2 estimate
    raw = sum(colSums(sampleSetChao2) > 0)
    q1 = sum(colSums(sampleSetChao2 > 0) == 1) # number of species occurring in one transect only
    q2 = sum(colSums(sampleSetChao2 > 0) == 2) # number of species occurring in two transect only
    m = sum(colSums(sampleSetChao2 > 0))
    Chao2_boot = raw + ((m-1)/m)*q1*(q1-1)/(2*(q2+1)) # Chao2 richness estimator

    # compute abundance-based coverage estimator (ACE)
    sampleSetACE = ceiling(sampleSetACE) # convert data to discrete count
    raw = sum(colSums(sampleSetACE) > 0)
    f1 = sum(colSums(sampleSetACE) == 1) # number of singleton species
    f2 = sum(colSums(sampleSetACE) == 2) # number of doubleton species
    S_rare = sum(colSums(sampleSetACE) <= 10) # number of rare species (< = 10 individuals)
    S_abund = sum(colSums(sampleSetACE) > 10) # number of rare species (>10 individuals)
    n_rare = sum(sampleSetACE[,colSums(sampleSetACE) <= 10]) # total number of individuals in rare species
    C_ACE = 1 - f1/n_rare # sample coverage estimate
    wsumfa = 0
    for( a in 1:10 ){
      wsumfa = wsumfa + a*(a-1)*sum(colSums(Community) == a)
    }
    V2 = max(((S_rare/C_ACE)*wsumfa/(n_rare*(n_rare - 1)) - 1),0) # coefficient of variation
    if( C_ACE > 0 ){
      ACE_boot = S_abund + S_rare/C_ACE + (f1/C_ACE)*V2
    } else {
      ACE_boot = Chao1
    }

    # compute jackknife abundance estimator
    sampleSetACE = ceiling(sampleSetACE) # convert data to discrete count
    raw = sum(colSums(sampleSetJK_a) > 0)
    f1 = sum(colSums(sampleSetJK_a) == 1) # number of singleton species
    f2 = sum(colSums(sampleSetJK_a) == 2) # number of doubleton species
    JK_a_boot = raw + 2*f1-f2

    # compute jackknife incidence estimator
    raw = sum(colSums(sampleSetJK_i) > 0)
    q1 = sum(colSums(sampleSetJK_i > 0) == 1) # number of species occurring in one transect only
    q2 = sum(colSums(sampleSetJK_i > 0) == 2) # number of species occurring in two transect only
    m = sum(colSums(sampleSetJK_i > 0))
    if( m == 1){
      JK_i_boot = rawli
    } else {
      JK_i_boot = raw + (q1*(2*m-3)/m-q2*((m-2)^2)/(m*(m-1)))
    }

    # Omega_taylor
    # compute correction terms for proposed approximation method
    raw = sum(colSums(sampleSetOmega_taylor) > 0)
    P_detected = rep(0,length(raw)) # array of zeros to record occupancy for each species
    for( species in 1:ncol(sampleSetOmega_taylor) ) {
      P_detected[species] = sum(sampleSetOmega_taylor[,species] > 0)/numSamplingUnit # occupancy as number of transects occupied divided by number of transects
    }
    P_detected[P_detected==0] = NA
    # compute on full dataset - means, variances, and covariances
    mean_P_detected = mean(P_detected[P_detected>0], na.rm = TRUE )
    var_P_detected  = var(P_detected[P_detected>0], na.rm = TRUE )
    n_m_detected = colMeans(sampleSetOmega_taylor[,colSums(sampleSetOmega_taylor)>0])
    mean_n_m_detected = mean(n_m_detected)
    var_n_m_detected  = var(colMeans(sampleSetOmega_taylor[,colSums(sampleSetOmega_taylor)>0]))
    cov_nm_P_detected = cov( data.frame(x = colMeans(sampleSetOmega_taylor[,colSums(sampleSetOmega_taylor)>0]), y = na.omit(P_detected)), use = "complete.obs")
    if( length(cov_nm_P_detected[,2])>1 ){
      cov_nm_P_detected = cov_nm_P_detected[1,2]
    } else {
      cov_nm_P_detected = 0
    }

    k = numSamplingUnit
    nm = mean_n_m_detected
    P = mean_P_detected
    # add "eval" before "(" below if using ke in place of k
    Taylor_detectP_terms = (c( eval(Di),
                            eval(d2Di_dnm2)*var_n_m_detected/2,
                            eval(d2Di_dP2)*var_P_detected/2,
                            eval(d2Di_dnmP)*cov_nm_P_detected )) # Approximated detection probability in community

    # check for correction term relative to some threshold
    if( sum(Taylor_detectP_terms, na.rm = T) > 0.1 ){ # if sum of correction terms is positive and greater than a threshold
      D_omega = sum(Taylor_detectP_terms, na.rm = T) # use full correction
    } else {
      D_omega = Taylor_detectP_terms[1] # else, use 0th order correction
    }
    if( sum(Taylor_detectP_terms, na.rm = T) > 1 ){
      D_omega = 1
    }
    Omega_taylor_boot = raw/D_omega

    # Omega_talor_0 using first detection term
    raw = sum(colSums(sampleSetOmega_taylor_0) > 0)
    P_detected = rep(0,length(raw)) # array of zeros to record occupancy for each species
    for( species in 1:ncol(sampleSetOmega_taylor_0) ) {
      P_detected[species] = sum(sampleSetOmega_taylor_0[,species] > 0)/numSamplingUnit # occupancy as number of transects occupied divided by number of transects
    }
    P_detected[P_detected==0] = NA
    # compute on full dataset - means, variances, and covariances
    mean_P_detected = mean(P_detected[P_detected>0], na.rm = TRUE )
    n_m_detected = colMeans(sampleSetOmega_taylor_0[,colSums(sampleSetOmega_taylor_0)>0])
    mean_n_m_detected = mean(n_m_detected)

    k = numSamplingUnit
    nm = mean_n_m_detected
    P = mean_P_detected
    Omega_taylor_0_boot = raw/eval(Di)

    # compute Omega
    nm = n_m_detected
    P = na.omit(P_detected)
    D_means = eval(Di)
    D_mean = mean(D_means, na.rm = T)
    Omega_boot = raw/D_mean

    expectedRichness_raw[resample]   = Richness_raw_boot
    expectedChao1[resample] = Chao1_boot
    expectedGP[resample]    = GP_boot
    expectedChao2[resample] = Chao2_boot
    expectedACE[resample]   = ACE_boot
    expectedJK_a[resample] = JK_a_boot
    expectedJK_i[resample] = JK_i_boot
    expectedOmega[resample] = Omega_boot
    expectedOmega_taylor[resample]   = Omega_taylor_boot
    expectedOmega_taylor_0[resample] = Omega_taylor_0_boot
}




return(data.frame(Richness_raw = expectedRichness_raw,
            Omega = expectedOmega,
            Omega_T = expectedOmega_taylor,
            Omega_T0 = expectedOmega_taylor_0,
            Chao1 = expectedChao1, Chao2 = expectedChao2,
            ACE = expectedACE, JK_a = expectedJK_a, JK_i = expectedJK_i))
}
