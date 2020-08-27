########################################
########################################
#
# explicit simulation
#
########################################
########################################

##########
# two loci, unlinked
##########
# population represented by 2N * n_loci matrix
make_pop <- function(n_individuals, n_loci, init_freq, random_draw = F){
  pop <- array(NA, dim = c(2*n_individuals, n_loci))
  
  if (random_draw){
    pop <- array(sample(0:1, size = 2*n_individuals*n_loci, replace = T), dim = c(2*n_individuals, n_loci))
  }
  else{
    n_0 <- ceiling(n_individuals * init_freq[1] * 2)
    n_1 <- 2 * n_individuals - n_0
    for (locus in 1:n_loci){
      alleles <- c(rep(0, n_0), rep(1, n_1))
      pop[,locus] <- sample(alleles)
    }
  }

  return (pop)
}


# take population --> select 2 parents --> recombine to produce one offspring
mate <- function(pop, c=.5,  replace = T, fn = X_over_2){
  N <- dim(pop)[1] / 2
  loci <- dim(pop)[2]
  parental_chr_1 <- pop[sample(1:(N*2), size = 1),]
  parental_chr_2 <- pop[sample(1:(N*2), size = 1),]
  return(fn(parental_chr_1, parental_chr_2))
}

X_over_1 <- function(parental_1, parental_2){
  loci <- length(parental_1)
  offspring <- array(NA, c(2, loci))
  
  for (l in 1:loci){
    offspring[,l] <- sample(c(parental_1[l], parental_2[l]), size = 2, replace = F)
  }
  return (offspring)
}

X_over_2 <- function(parental_1, parental_2, c = .5){
  loci <- length(parental_1)
  offspring_chr_1 <- parental_1
  offspring_chr_2 <- parental_2
  X_over_RNG <- runif(n=loci-1)
  for (l in 2:loci){
    if (X_over_RNG[l-1] < c){
      chr_1_region <- offspring_chr_1[l:loci]
      chr_2_region <- offspring_chr_2[l:loci]
      offspring_chr_1[l:loci] <- chr_2_region
      offspring_chr_2[l:loci] <- chr_1_region
    }
  }
  return (rbind(offspring_chr_1, offspring_chr_2))
}

run_gen <- function(pop, N = "CONSTANT"){
  if (N == "CONSTANT") N = dim(pop)[1] / 2
  next_pop <- array(NA, dim = c(2*N, dim(pop)[2]))
  for (i in 1:N){
    next_pop[(2*i-1):(2*i),] <- mate(pop)
  }
  return(next_pop)
}


run_gens <- function(pop, gens, S_gens,S, N = "CONSTANT"){
  if (N == "CONSTANT") N = dim(pop)[1] / 2
  
  pop_list <- list()
  samp_list <- list()
  
  next_pop <- pop
  
  for (i in 1:(gens+1)){
    if (i %in% S_gens){
      
      next_pop <- run_gen(next_pop, N+S)
      samp_list[[i]] <- next_pop[1:S,]
      pop_list[[i]] <- next_pop[S+1:S+N,]
      
    }
    else{
      next_pop <-  run_gen(next_pop, N)
      pop_list[[i]] <- next_pop
      
    }
    
  }
  return(list("pops" = pop_list, "samples" = samp_list))
}


to_haplotype_freqs <- function(pop){
  n_chr <- dim(pop)[1]
  pAB <- 0
  pAb <- 0
  paB <- 0
  pab <- 0
  for (r in 1:n_chr){
    
    if (all(pop[r,] == c(1,1))){
      pAB <- pAB + 1
    }
    if (all(pop[r,] == c(1,0))){
      pAb <- pAb + 1
    }
    if (all(pop[r,] == c(0,1))){
      paB <- paB + 1
    }
    if (all(pop[r,] == c(0,0))){
      pab <- pab + 1
    }
  }
  return(list(pAB = pAB/n_chr, pAb = pAb/n_chr, paB = paB/n_chr, pab = pab/n_chr))
}
 

freqs_to_r2 <- function(pAB, pAb, paB, pab){
  pA <- pAB + pAb
  pB <- pAB + paB
  r <- (pAB*pab-pAb*paB)**2 / (pA*(1-pA)* pB*(1-pB))
  return(r)
}


r2_drift <- function(r_vec, S){
  mean_r <- mean(r_vec)
  r_drift <- (mean_r -  1/(2*S) ) / (1 - 1/(2*S))
  return(r_drift)
}


get_LD_est <- function(pAB, pAb, paB, pab, S){
  r_vec <- freqs_to_r2(pAB, pAb, paB, pab)
  r_drift <- r2_drift(r_vec, S)
  Ne_est <- 1/(3*r_drift)
  return(Ne_est)
}



# both with and without sampling

## our sample is a sample WITH REPLACEMENT from the gamete pool

reps <- 1000
r2_vec_full <- rep(NA, reps)
r2_vec_sample <- rep(NA, reps)
S = 50
N = 100
for (r in 1:reps){
  print(r)
  pop <- make_pop(n_individuals = N, n_loci = 2, init_freq = .5)
  pop_list <- run_gens(pop = pop, gens = 5, S_gens = c(), S = 0)
  final_pop <- pop_list$pops[[6]]
  sample_pop <- final_pop[sample(1:(2*N), size = 2*S, replace = T),]
  full_hap_freqs <- to_haplotype_freqs(pop = final_pop)
  sample_hap_freqs <- to_haplotype_freqs(pop = sample_pop)
  full_r2 <- freqs_to_r2(full_hap_freqs$pAB, full_hap_freqs$pAb, full_hap_freqs$paB, full_hap_freqs$pab)
  sample_r2 <- freqs_to_r2(sample_hap_freqs$pAB, sample_hap_freqs$pAb, sample_hap_freqs$paB, sample_hap_freqs$pab)
  r2_vec_full[r] <- full_r2
  r2_vec_sample[r] <- sample_r2
}      

full_mean_r2 <- mean(r2_vec_full)
sample_mean_r2 <- mean(r2_vec_sample)
sample_drift_r2 <- (sample_mean_r2 -  1/(2*S) ) / (1 - 1/(2*S))
full_est <- 1/(3*full_mean_r2)
full_est
sample_est <- 1/(3*sample_drift_r2)
sample_est



#############
# more than 2 loci, potentially linked
#############

# for more than 2 loci, need to change the way that r2 is calculated - ie calculate matrix for all pairs

# make_pop(), mate remain the same

get_r2_matrix <- function(pop){
  loci = dim(pop)[2]
  r2_matrix = array(NA, dim = c(loci, loci))
  for (locus_1 in 1:loci){
    genotype_1 <- pop[,locus_1]
    for (locus_2 in 1:(locus_1-1)){
      if (locus_1>locus_2){
        genotype_2 <- pop[, locus_2]
        if (length(unique(genotype_1)) == 2 & length(unique(genotype_1)) == 2){
          n = length(genotype_1)
          pA = sum(genotype_1) / n
          pB = sum(genotype_2) / n
          pAB= sum(genotype_1 & genotype_2) / n
          pab= sum( (1-genotype_1) & (1-genotype_2) ) / n
          pAb= sum(genotype_1 & (1-genotype_2) ) / n
          paB= sum( (1-genotype_1) & genotype_2 ) / n
          
          r2_phased   = (pAB * pab - pAb * paB)**2 / (pA * (1-pA) * pB * (1-pB))
          r2_matrix[locus_1, locus_2] <- r2_phased
        }
      }
    }
  }
  
  return (r2_matrix)
}

reps <- 1000
r2_vec_full <- rep(NA, reps)
r2_vec_sample <- rep(NA, reps)
S = 50
N = 100
n_loci = 5
for (r in 1:reps){
  print(r)
  pop <- make_pop(n_individuals = N, n_loci = n_loci, init_freq = .5)
  pop_list <- run_gens(pop = pop, gens = 5, S_gens = c(), S = 0)
  final_pop <- pop_list$pops[[6]]
  sample_pop <- final_pop[sample(1:(2*N), size = 2*S, replace = T),]
  full_r2_matrix <- get_r2_matrix(final_pop)
  sample_r2_matrix <- get_r2_matrix(sample_pop)
  r2_vec_full[r] <- mean(full_r2_matrix, na.rm=T)
  r2_vec_sample[r] <- mean(sample_r2_matrix, na.rm=T)
}      

full_mean_r2 <- mean(r2_vec_full)
sample_mean_r2 <- mean(r2_vec_sample)
sample_drift_r2 <- (sample_mean_r2 -  1/(2*S) ) / (1 - 1/(2*S))
full_est <- 1/(3*full_mean_r2)
full_est
sample_est <- 1/(3*sample_drift_r2)
sample_est









#############
# with migration
#############



# produce next_gen, and also a sample of next_gen . produced with a separate binomial draw from the previosus generation.
run_gen_sample <- function(pop,S, N = "CONSTANT"){
  if (N == "CONSTANT") N = dim(pop)[1] / 2
  next_pop <- array(NA, dim = c(2*N, dim(pop)[2]))
  for (i in 1:N){
    next_pop[(2*i-1):(2*i),] <- mate(pop)
  }
  sample_pop <- array(NA, dim = c(2*S, dim(pop)[2]))
  for (i in 1:S){
    sample_pop[(2*i-1):(2*i),] <- mate(pop)
  } 
  return(list("pop" = next_pop, "sample" = sample_pop))
}

# runs simulation, returning a list of the intermediate populations and a list of the requested samples
run_gens <- function(pop, gens, S_gens,S, N = "CONSTANT"){
  if (N == "CONSTANT") N = dim(pop)[1] / 2
  
  pop_list <- list()
  samp_list <- list()
  
  next_pop <- pop
  
  for (i in 1:(gens+1)){
    if (i %in% S_gens){
      
      next_pop <- run_gen(next_pop, N+S)
      samp_list[[i]] <- next_pop[1:S,]
      pop_list[[i]] <- next_pop[S+1:S+N,]
      
    }
    else{
      next_pop <-  run_gen(next_pop, N)
      pop_list[[i]] <- next_pop
      
    }

  }
  return(list("pops" = pop_list, "samples" = samp_list))
}

# nsplits populatin - ot used
# take_sample <- function(pop, S){
#   N = dim(pop)[1] / 2
#   S_idxs <- sample(1:N, size = S, replace = F)
#   S_idxs_genomes <- c(rbind((S_idxs*2-1), (S_idxs*2)))
#   sample_pop <- pop[S_idxs_genomes,]
#   new_pop <- pop[-S_idxs_genomes,]
#   return(list('sample' = sample_pop, 'pop' = new_pop))
# }

# calculate all pairwise r2 values for a population
calc_r2 <- function(pop){
  r2_matrix = array(NA, dim = rep(dim(pop)[2],2))
  for (locus_1 in 1:dim(pop)[2]){
    genotype_1 <- pop[,locus_1]
    for (locus_2 in 1:(locus_1-1)){
      if (locus_1>locus_2){
        genotype_2 <- pop[, locus_2]
        if (length(unique(genotype_1)) == 2 & length(unique(genotype_1)) == 2){
          n = length(genotype_1)
          pA = sum(genotype_1) / n
          pB = sum(genotype_2) / n
          pAB= sum(genotype_1 & genotype_2) / n
          pab= sum( (1-genotype_1) & (1-genotype_2) ) / n
          pAb= sum(genotype_1 & (1-genotype_2) ) / n
          paB= sum( (1-genotype_1) & genotype_2 ) / n
          
          r2_phased   = (pAB * pab - pAb * paB)**2 / (pA * (1-pA) * pB * (1-pB))
          r2_matrix[locus_1, locus_2] <- r2_phased
        }
      }
    }
  }

  return (r2_matrix)
}

# calculates mean r2 value for a population
get_mean_r2 <- function(pop){
  r2_matrix <- calc_r2(pop)
  return (mean(r2_matrix, na.rm = T))
}

# calculates r2 attributable to drift for a population
get_r2_drift <- function(samp){
  S <- dim(samp)[1] / 2
  r2_mean <- get_mean_r2(samp)
  r2_drift <- (r2_mean -  1/(2*S) ) / (1 - 1/(2*S))
  return (r2_drift)
}

# returns the LD estimate of Ne
estimate_LD <- function(samp){
  r2_drift <- get_r2_drift(samp)
  Ne_est = 1/(3*r2_drift)
  return (Ne_est)
}

to_haplotype_freqs <- function(pop){
  n_chr <- dim(pop)[1]
  pAB <- 0
  pAb <- 0
  paB <- 0
  pab <- 0
  for (r in 1:n_chr){

    if (all(pop[r,] == c(1,1))){
      pAB <- pAB + 1
    }
    if (all(pop[r,] == c(1,0))){
      pAb <- pAb + 1
    }
    if (all(pop[r,] == c(0,1))){
      paB <- paB + 1
    }
    if (all(pop[r,] == c(0,0))){
      pab <- pab + 1
    }
  }
  return(list(pAB = pAB/n_chr, pAb = pAb/n_chr, paB = paB/n_chr, pab = pab/n_chr))
}

# returns the mean of estiamtes of Ne from a pair of samples
estimate_LD_pair <- function(samp1, samp2){
  r2_drift_1 <- get_r2_drift(samp1)
  r2_drift_2 <- get_r2_drift(samp2)
  Ne_est = 1/(3*mean(r2_drift_1, r2_drift_2))
  
  return (Ne_est)
}

freq_to_fc <- function( A_t1, A_t2, B_t1, B_t2){
  return( ( ( (A_t1 - A_t2)**2 / ((A_t1 + A_t2)/2 - A_t1 * A_t2) ) + ( (B_t1-B_t2)**2 / ((B_t1+B_t2)/2 - B_t1*B_t2) ) ) / 2 )
}

# calculate the F_c value for genotype samples 1 and 2
get_FC_locus <- function(genotype_t1, genotype_t2){
  
  n <- length(genotype_t1)
  stopifnot( n == length(genotype_t2) )
  A_t1 <- sum(genotype_t1) / n
  A_t2 <- sum(genotype_t2) / n
  B_t1 <- 1-A_t1
  B_t2 <- 1-A_t2
  
  return (freq_to_fc( A_t1, A_t2, B_t1, B_t2))
}

# calculate mean FC value between 2 samples
get_mean_FC <- function(samp1, samp2){
  n_loci <- dim(samp1)[2]
  FC_total <- 0
  count <- 0
  for (i in 1:n_loci){
    FC_total <- FC_total + get_FC_locus(samp1[,i], samp2[,i])
    count <- count + 1
  }
  return (FC_total / count)
}

# calcualte estimate of Ne from 2 samples and the number of generations between them
estimate_AF <- function(samp1, samp2, t){
  S0 <- dim(samp1)[1] / 2
  St <- dim(samp2)[1] / 2
  mean_FC <- get_mean_FC(samp1, samp2)
  INV_est <- (mean_FC - 1/(2*S0) - 1/(2*St) ) / t
  Ne_est <-  t / (2*(mean_FC - 1/(2*S0) - 1/(2*St)))
  return (INV_est)
}

# get allele_frequencies from a population
get_freqs <- function(pop){
  return (colMeans(pop))
}

#######################################
if (!interactive()){
  

# TESTING AF ESTIAMTOR

repeats = 100
N = 100
S = 50
n_loci = 2
S_gens <- c(3)
estimates_AF <- c()
estimates_LD <- c()
FS_pop <- c()
FS_sample <- c()
FS_drift <- c()
FS_r2_arr <- array(NA, dim = c(repeats, max(S_gens+1)))
pAB_m <- array(NA, dim = c(repeats, max(S_gens+1)))

source("TY_SIM.R")

for (r in 1:repeats){
  print(r)
  print(paste0("\t", "Simulating"))
  gene_pool <- make_pop(n_individuals = N, n_loci =  n_loci, init_freq = c(.5,.5), random_draw = F)
  out <- run_gens(gene_pool, max(S_gens), S_gens = S_gens, S = S)
  samples <- out$samples
  
  popn_freqs <- to_haplotype_freqs(out$pops[[3]])
  samp_freqs <- to_haplotype_freqs(samples[[S_gens[1]]])
  
  popn_r2 <- freqs_to_r2(popn_freqs$pAB,popn_freqs$pAb,popn_freqs$paB,popn_freqs$pab )
  samp_r2 <- freqs_to_r2(samp_freqs$pAB,samp_freqs$pAb,samp_freqs$paB,samp_freqs$pab )
  
  FS_pop <- c(FS_pop, popn_r2)
  FS_sample <- c(FS_sample, samp_r2)
  
  # popn_est <- get_LD_est(popn_freqs$pAB,popn_freqs$pAb,popn_freqs$paB,popn_freqs$pab, 100000000000)
  # samp_est <- get_LD_est(samp_freqs$pAB,samp_freqs$pAb,samp_freqs$paB,samp_freqs$pab, 50)
  # 
  # estimates_LD <- c(estimates_LD, popn_est)
  # estimates_AF <- c(estimates_AF, samp_est)  
  
  # FS_pop <- c(FS_pop, get_mean_r2(out$pops[[3]]))
  # FS_sample<-c(FS_sample, get_mean_r2(samples[[S_gens[1]]]))
  # FS_drift<-c(FS_drift, get_r2_drift(samples[[S_gens[1]]]))
  
  # pAB_temp <- c()
  # for (g in 1:11){
  #   pAB_temp <- c(pAB_temp, mean(rowSums(out$pops[[g]]) == 2))
  #   print(pAB_temp)
  # }
  # print(length(pAB_temp))
  # print(length(pAB_v))
  # pAB_m <- rbind(pAB_m,  pAB_temp)
  # 
  # for (pop_idx in 1:length(out$pops)){
  #   FS_r2_arr[r,pop_idx] <- get_mean_r2(out$pops[[pop_idx]])
  # }
  # 
  print(paste0("\t", "AF est"))
  #est_AF <- estimate_AF(samples[[S_gens[1]]], samples[[S_gens[2]]], diff(S_gens))
  print(paste0("\t", "LD est"))
  #est_LD <- estimate_LD_pair(samples[[S_gens[1]]], samples[[S_gens[2]]])
  est_LD <- estimate_LD(samples[[S_gens[1]]])
  
  #estimates_AF <- c(estimates_AF, est_AF)
  estimates_LD <- c(estimates_LD, est_LD)
  
}


hist(pAB_m[,5], xlim = c(0,.5), breaks = 20)
hist(o$pAB[,5], xlim = c(0.,.5), breaks = 20)



pAB_mean <- pAB_v / repeats
plot(1:11, pAB_mean, type = 'l', ylim = c(0.23,0.27))
lines(1:11, colMeans(o$pAB),type = 'l')


plot(1:8, colMeans(FS_r2_arr)-1/(2*N), type = 'l', ylim = c(0,0.002))
lines(1:8, colMeans(o$r)-1/(2*N), type = 'l', ylim = c(0,0.01), col = 'red')
abline(h = exp_pop-1/(2*N), lty = 2)
abline(h = (exp_pop-1/(2*N))/2, lty = 2)


pAB_time <- c()
for (g in 1:8){
  pAB_time <- c(pAB_time, mean(rowSums(out$pops[[g]]) == 1))
}

plot(1:8, pAB_time, type = 'l', ylim = c(.4,.6))
lines(1:8, 1-colMeans(o$pab)-colMeans(o$pAB))

pAB_time <- c()
for (g in 1:8){
  pAB_time <- c(pAB_time, mean(rowSums(out$pops[[g]]) == 2))
}

plot(1:8, pAB_time, type = 'l', ylim = c(0,.5))
lines(1:8, colMeans(o$pAB))

par(mfrow = c(2,1), mai = rep(.5,4))

hist(1/(2*estimates_AF), xlim = c(50,200), breaks = round(diff(range(1/(2*estimates_AF)))/10) )
abline(v = mean(1/(2*estimates_AF), col = "blue"))
abline(v = N, col = 'red')

hist(estimates_LD, xlim = c(50,200), breaks = round(diff(range(estimates_LD))/10) )
abline(v = mean(estimates_LD), col = "blue")
abline(v = N, col = 'red')




hist(estimates_AF, breaks = 20)#, breaks = diff(range(estimates_AF))/20)
abline(v = mean(estimates_AF, col = "blue"))
abline(v = 1/(2*N), col = 'red')


}

