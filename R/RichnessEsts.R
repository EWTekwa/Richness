RichnessEsts <- function( Community ){

# RichnessEsts.R is based on RichnessEsts.m
  # Matt Whalen rewrote matlab script for R - started Mar 27, 2022
# function calculates several estimates of richness:
  # - raw richness
  # - newly proposed method
  # - Chao1
  # - Gamma-Poisson
  # - Chao2
  # - Abundance-based coverage estimator (ACE)
  # - Jackknife abundance estimator
  # - Jackknife incidence estimator
  # - three versions of the estimator introduced in Tekwa et al. 2023
    # - Omega - exact estimator, recommended version
    # - Omega_T - second-order Taylor approximated estimator, useful for quantifying biases from means and variances in occupancy and observed abundance across species
    # - Omega_not - zeroth-order approximated estimator

# The new method calculates richness based on spatial Community data:
  # rows=transects, columns=species, values=individual counts





## write observational process model
# Dix = 1-exp(-nm/P) # local detection probability of species i at sampled site x
# Dix = Dixo*P
# Di = 1-(1-Dix)^k # detection probability of species i across all sampled sites x in community s
Di = expression( 1-(1-((1-exp(-nm/P))*P))^k )

## second-order derivatives
# d2Di_dnm2 = D(D( Di, "nm" ), "nm")
# d2Di_dP2  = D(D( Di, "P" ), "P")
# d2Di_dnmP = D(D( Di, "nm" ), "P")

# hard code second-order derivatives above for speed
d2Di_dnm2 = expression( -((1 - ((1 - exp(-nm/P)) * P))^(k - 1) *
                            (k * (exp(-nm/P) * (1/P) * (1/P) * P)) +
                            (1 - ((1 - exp(-nm/P)) * P))^((k - 1) - 1) *
                            ((k - 1) * (exp(-nm/P) * (1/P) * P)) *
                            (k * (exp(-nm/P) * (1/P) * P))) )
d2Di_dP2  = expression( -((1 - ((1 - exp(-nm/P)) * P))^(k - 1) *
                            (k * (exp(-nm/P) * (nm/P^2) +
                            ((exp(-nm/P) * (nm/P^2) * (nm/P^2) - exp(-nm/P) *
                            (nm * (2 * P)/(P^2)^2)) * P + exp(-nm/P) *
                            (nm/P^2)))) + (1 - ((1 - exp(-nm/P)) * P))^((k - 1) - 1) *
                            ((k - 1) * ((1 - exp(-nm/P)) - exp(-nm/P) * (nm/P^2) * P)) *
                            (k * ((1 - exp(-nm/P)) - exp(-nm/P) * (nm/P^2) * P))))
d2Di_dnmP = expression( (1 - ((1 - exp(-nm/P)) * P))^(k - 1) *
                            (k * ((exp(-nm/P) * (nm/P^2) * (1/P) -
                                   exp(-nm/P) * (1/P^2)) * P + exp(-nm/P) * (1/P))) -
                            (1 - ((1 - exp(-nm/P)) * P))^((k - 1) - 1) *
                            ((k - 1) * ((1 - exp(-nm/P)) - exp(-nm/P) * (nm/P^2) * P)) *
                            (k * (exp(-nm/P) * (1/P) * P)) )

numSamplingUnit =  nrow(Community) # get number of transects
Richness_raw = sum(colSums(Community)>0) # get raw richness
P_detected = rep(0,ncol(Community)) # array of zeros to record occupancy for each species
for( species in 1:ncol(Community) ) {
  P_detected[species] = sum(Community[,species] > 0)/numSamplingUnit # occupancy as number of transects occupied divided by number of transects
}
P_detected[P_detected==0] = NA

# compute on full dataset - means, variances, and covariances
mean_P_detected = mean(P_detected[P_detected>0], na.rm = TRUE )
var_P_detected  = var(P_detected[P_detected>0], na.rm = TRUE )
n_m_detected = colMeans(as.data.frame(Community[,colSums(Community)>0]))
mean_n_m_detected = mean(n_m_detected)
var_n_m_detected  = var(colMeans(as.data.frame(Community[,colSums(Community)>0])))
cov_nm_P_detected = cov( data.frame(x = colMeans(as.data.frame(Community[,colSums(Community)>0])), y = na.omit(P_detected)), use = "complete.obs")
if( length(cov_nm_P_detected[,2])>1 ){
  cov_nm_P_detected = cov_nm_P_detected[1,2]
} else {
  cov_nm_P_detected = 0
}

# calculate singletons, doubletons, and the Chao1 richness estimator
Community = ceiling(Community)
f1 = sum(colSums(Community) == 1) # number of singleton species
f2 = sum(colSums(Community) == 2) # number of doubleton species
Chao1 = Richness_raw+f1*(f1-1)/(2*(f2+1)) # Chao1 richness estimator

# compute GP (Gamma-Poisson mixture) estimate (Chiu 2023, peerJ)
f3 = sum(colSums(Community)==3) # number of tripleton species
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
  GP = Richness_raw + f0Chao1*max(c(.5, A))
} else {
  GP = Richness_raw + f0Chao1
}

# compute Chao2 estimate
q1 = sum(colSums(Community >= 1) == 1) # number of species occurring in one sample only
q2 = sum(colSums(Community >= 1) == 2) # number of species occurring in two samples only
m = sum(colSums(Community >= 1)) # total number of samples
Chao2 = Richness_raw + ((m-1)/m)*q1*(q1-1)/(2*(q2+1)) # Chao2 richness estimator

# compute abundance-based coverage estimator (ACE)
S_rare = sum(colSums(Community) <= 10 & colSums(Community) > 0 ) # number of rare species (<=10 individuals)
S_abund = sum(colSums(Community) > 10) # number of abundant species (>10 individuals)
n_rare = sum(Community[,colSums(Community) <= 10 & colSums(Community) > 0]) # total number of individuals in rare species
if( n_rare > 0 ){
  C_ACE = 1-f1/n_rare # sample coverage estimate
} else {
  C_ACE = 1
}
wsumfa = 0
for( a in 1:10 ){
  wsumfa = wsumfa + a*(a-1)*sum(colSums(Community) == a)
}
V2 = max(((S_rare/C_ACE)*wsumfa/(n_rare*(n_rare - 1)) - 1),0) # coefficient of variation
if( C_ACE > 0 ){
  ACE = S_abund + S_rare/C_ACE + (f1/C_ACE)*V2
} else {
  ACE = Chao1
}

# compute jackknife abundance estimator
JK_a = Richness_raw+2*f1-f2

# compute jackknife incidence estimator
if( m == 1){
  JK_i = Richness_raw
} else {
  JK_i = Richness_raw + (q1*(2*m-3)/m-q2*((m-2)^2)/(m*(m-1)))
}


# compute correction terms for Omega_T
k  = numSamplingUnit
nm = mean_n_m_detected
P  = mean_P_detected
Omega_detectP_terms = c( eval(Di),
  eval(d2Di_dnm2)*var_n_m_detected/2,
  eval(d2Di_dP2)*var_P_detected/2,
  eval(d2Di_dnmP)*cov_nm_P_detected ) # Approximated detection probability in community

# check for correction term relative to some threshold
if( sum(Omega_detectP_terms, na.rm = T) > 0.1 ){ # if sum of correction terms is positive and greater than a threshold
  D_omega = sum(Omega_detectP_terms, na.rm = T) # use full correction
} else {
  D_omega = Omega_detectP_terms[1] # else, use 0th order correction
}
if( sum(Omega_detectP_terms, na.rm = T) > 1 ){
  D_omega = 1
}
if( sum(Omega_detectP_terms, na.rm = T) < 0.1 ){
  D_omega = 0.1
}
Omega_taylor   = Richness_raw/D_omega
Omega_taylor_0 = Richness_raw/Omega_detectP_terms[1]


# assign the mean states of correction terms to return with the function
meanStates = c(mean_n_m_detected,mean_P_detected,
                var_n_m_detected,var_P_detected,
                cov_nm_P_detected )


# compute Omega
nm = n_m_detected
P = na.omit(P_detected)
D_means = eval(Di) # observation process model for each detected species
D_mean = mean(D_means)
Omega = Richness_raw/D_mean


if( Richness_raw == 0 ){
  Richness_raw = NA
  Richness_taylor = NA
  Richness_taylor_0 = NA
  Richness_omega = NA
  Chao1 = NA
  GP = NA
  Chao2 = NA
  ACE = NA
  JK_a = NA
  JK_i = NA
}




return(list(data.frame(Richness_raw = Richness_raw,
                  Chao1 = Chao1, GP = GP, Chao2 = Chao2,
                  ACE = ACE, JK_a = JK_a, JK_i = JK_i,
                  Omega = Omega, Omega_taylor = Omega_taylor, Omega_taylor_0 = Omega_taylor_0),
            meanStates = meanStates,
            Omega_detectP_terms = Omega_detectP_terms))
}
