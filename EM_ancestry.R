## 2.1. ANCESTRY BASED ON GENOTYPE LIKELIHOODS for ASW NA19818: model implementation ##

# Number of copies of each allele per site (SNP)
A <- LHs_ASW_NA19818$allele1
B <- LHs_ASW_NA19818$allele2

# Likelihoods for each genotype
gAA <- LHs_ASW_NA19818$Ind0.AA
gAB <- LHs_ASW_NA19818$Ind0.AB
gBB <- LHs_ASW_NA19818$Ind0.BB

like <- rbind(gAA, gAB, gBB)

# Starting parameter for the EM algorithm
pk <- c(0.5)

# Numeric optimization EM under the ALTERNATIVE model
EM_alternative <- function(pk, like) {
  
  # Likelihood * admixture proportions * genotype frequencies (assuming HW equilibrium: alleles frequencies independent)
  w0 <- (like[1, ] * admixture_prop[1, 1] * (1 - (0.6104 * pk))^2) / 3          # Homozygous 1
  w1 <- (like[2, ] * admixture_prop[1, 2] * 2 * pk * (1 - (0.6104 * pk))) / 3   # Heterozygous
  w2 <- (like[3, ] * admixture_prop[1, 3] * pk^2) / 3                           # Homozygous 2
  
  # Posterior calculation (based on homozygous 2):
  # Expectation (of the homozygous 2 genotype divided by the sum of the individual LHs)
  pkNew <- mean((w1 + 2*w2) / (2 * (w0 + w1 + w2)))
  
  return(pkNew)
}

average_pk <- ((0.6104 * pk) + (0.3896 * pk) + (0.0045 * pk)) / 3

# Numeric optimization EM under the NULL model
EM_null <- function(pk, like) {
  
  # Likelihood * admixture proportions * genotype frequencies (assuming HW equilibrium: alleles frequencies independent)
  w0 <- (like[1, ] * admixture_prop[1, 1] * (1 - average_pk)^2) / 3          # Homozygous 1
  w1 <- (like[2, ] * admixture_prop[1, 2] * 2 * pk * (1 - average_pk)) / 3   # Heterozygous
  w2 <- (like[3, ] * admixture_prop[1, 3] * pk^2) / 3                        # Homozygous 2
  
  # Posterior calculation (based on BB homozygous):
  # Expectation (of the BB genotype divided by the sum of the individual LHs)
  pkNew <- mean((w1 + 2*w2) / (2 * (w0 + w1 + w2)))
  
  return(pkNew)
}

# 20 EM iterations updating the parameter: find the optimal one
for(i in 1:20)
  print(pk <- EM_alternative(pk, like))  # Alternative model converging with pk = 0.5793603

for(i in 1:20)
  print(pk <- EM_null(pk, like))         # Null model converging with pk = 0.5704053


## 2.2. LHRT ##

LHRT <- 2 * (log(0.5793603) - log(0.5704053))

pvalue <- 1 - pchisq(LHRT, 1)

# As seen by the p-value obtained (0.86 with 1 degree of freedom, or 1 parameter of difference between the two
# alternative and the null model), the ASW NA19818 individual does not have a significant amount of Chinese ancestry


## 2.1. ANCESTRY BASED ON GENOTYPE LIKELIHOODS for ASW NA19818: model implementation ##

# Number of copies of each allele per site (SNP)
A <- LHs_ASW_NA19818$allele1
B <- LHs_ASW_NA19818$allele2

# Likelihoods for each genotype
gAA <- LHs_ASW_NA19818$Ind0.AA
gAB <- LHs_ASW_NA19818$Ind0.AB
gBB <- LHs_ASW_NA19818$Ind0.BB

like <- rbind(gAA, gAB, gBB)

# Starting parameter for the EM algorithm
pk <- c(0.5)

# Numeric optimization EM under the ALTERNATIVE model
EM_alternative <- function(pk, like) {
  
  # Likelihood * admixture proportions * genotype frequencies (assuming HW equilibrium: alleles frequencies independent)
  w0 <- (like[1, ] * admixture_prop[1, 1] * (1 - (0.6104 * pk))^2) / 3          # Homozygous 1
  w1 <- (like[2, ] * admixture_prop[1, 2] * 2 * pk * (1 - (0.6104 * pk))) / 3   # Heterozygous
  w2 <- (like[3, ] * admixture_prop[1, 3] * pk^2) / 3                           # Homozygous 2
  
  # Posterior calculation (based on homozygous 2):
  # Expectation (of the homozygous 2 genotype divided by the sum of the individual LHs)
  pkNew <- mean((w1 + 2*w2) / (2 * (w0 + w1 + w2)))
  
  return(pkNew)
}


