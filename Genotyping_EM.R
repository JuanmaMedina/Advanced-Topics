
# Assumed error rate
e <- 0.01

# Number of copies of each allele per individual
A <- c(7, 25, 5, 4, 0, 0)
B <- c(0, 1, 3, 4, 2, 4)

# Likelihoods for each genotype, assuming the error is the same for all the reads, 
# (according the GATK method). I.e: (1 - e) if the observed base matches true genotype
# and (e / 3) if it does not
gAA <- (1 - e)^A * (e / 3)^B
gAB <- (0.5 * (1 - e + e/3) )^(A+B)
gBB <- (1 - e)^B * (e / 3)^A

like <- rbind(gAA, gAB, gBB)

# Numeric optimization EM (takes PARAMETER and initial pre-calculated
# LIKELIHOOD of the genotypes as input)
EMgeno <- function(pk, like) {
  
  # Likelihood * priors (assuming HW equilibrium: alleles frequencies independent)
  # Q(Z), for Z = 0, 1, 2: distribution over hidden states (genotypes)
  w0 <- like[1, ] * (1 - pk)^2          # AA homozygous
  w1 <- like[2, ] * 2 * pk * (1 - pk)   # AB heterozygous
  w2 <- like[3, ] * pk^2                # BB homozygous
  
  # Posterior calculation (based on BB homozygous): 
  # Expectation (of the BB genotype divided by the sum of the individual LHs)
  pkNew <- mean((w1 + 2*w2) / (2 * (w0 + w1 + w2)))
  
  return(pkNew)
}

# Starting parameter for the EM algorithm
pk <- c(0.2)

# 10 EM iterations updating the parameter: find the optimal one
for(i in 1:10)
  print(pk <- EMgeno(pk, like))

# Convergence with parameter = 0.46
