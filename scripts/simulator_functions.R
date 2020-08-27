library("matrixStats")
library(Rcpp)

sourceCpp(file = "propegate_new.cpp")


#
# EXAMPLE SIMULATORS
#

# one locus
sim.drift<-function(p0=0.5, N=100, t=10)
{
  p<-rep(NA, t+1)
  p[1]<-p0
  for (i in 1:t)
  {p[i+1]<-rbinom(1, size=2*N, prob=p[i])/(2*N)}
  return(p)			
}

# one locus, multiple loci
sim.drift2<-function(p0=0.5, loci=8, N=100, t=10)
{
  p<-matrix(nr=loci, nc=t+1)
  p[,1]<-p0
  for (i in 1:t)
  {p[,i+1]<-rbinom(loci, size=2*N, prob=p[,i])/(2*N)}
  return(p)
}

# one locus, multiple loci, 2 subpopulations
sim.drift3<-function(p0=0.5, loci=8, N_A=100, N_B=100, m_A=0.1, m_B=0.1, t=10)
{
  p_A<-matrix(nr=loci, nc=t+1)
  p_B<-matrix(nr=loci, nc=t+1)
  p_A[,1]<-p0
  p_B[,1]<-p0
  for (i in 1:t)
  {
    p_A[,i+1]<-rbinom(loci, size=2*N_A, prob=(1-m_A)*p_A[,i]+m_A*p_B[,i])/(2*N_A)
    p_B[,i+1]<-rbinom(loci, size=2*N_B, prob=(1-m_B)*p_B[,i]+m_B*p_A[,i])/(2*N_B)
  }
  return(list(p_A=p_A, p_B=p_B))
}


# one locus, multiple loci, multiple subpopulations

sim.drift.migration<-function(p0=0.5, loci=8, N_v = c(100,100), M = array(c(.9,.1,.1,.9), dim = c(2,2)), t=10)
{
  if (!check_M(M)) return("Bad M")
  p_pre <- array(data = NA, dim = c(loci, t+1, length(N_v)))
  p_post <- array(data = NA, dim = c(loci, t+1, length(N_v)))
  p_pre[,1,] <- p0
  for (i in 1:t)
  {
    
    p_post[,i,] <- p_pre[,i,] %*% M
    samp <- rbinom(n= length(p_post[,i,]), size = 2*N_v, prob = p_post[,i,])
    samp_arr <- array(samp, dim = c(length(N_v),loci))/(2*N_v)
    p_pre[,i+1,] <- samp_arr
    
  }
  return(list(pre = p_pre, post = p_post))
}

# one pair of loci (four halpotypes)
sim.drift4<-function(p0=c(0.25, 0.25, 0.25, 0.25), N=100, c=0.5, t=10)
{
  pAB<-rep(NA, t+1)
  pAb<-rep(NA, t+1)
  paB<-rep(NA, t+1)
  pab<-rep(NA, t+1)
  pA<-rep(NA, t+1)
  pB<-rep(NA, t+1)
  D<-rep(NA, t+1)
  r<-rep(NA, t+1)
  # FILL THESE VECTORS WITH INITIAL VALUES
  pAB[1]<-p0[1]
  pAb[1]<-p0[2]
  paB[1]<-p0[3]
  pab[1]<-p0[4]
  pA[1] <- pAB[1] + pAb[1]
  pB[1] <- pAB[1] + paB[1]
  D[1]<-(pAB[1]*pab[1]-pAb[1]*paB[1])
  r[1]<-(pAB[1]*pab[1]-pAb[1]*paB[1])**2 / (pA[1]*(1-pA[1])* pB[1]*(1-pB[1]))
  # PROPAGATION. NEED TO CONSIDER BOTH DRIFT AND RECOMBINATION
  for (i in 1:t)
  {
    # CALCULATE THE GAMETIC FREQ IN THE GAMETE POOL
    # WILL BE ERASED AND RE-WRITTEN IN EACH TIME STEP
    gameteAB<-pAB[i]-c*D[i]
    gameteAb<-pAb[i]+c*D[i]
    gameteaB<-paB[i]+c*D[i]
    gameteab<-pab[i]-c*D[i]
    # DRIFT, RANDOM SAMPLING FROM THE GAMETE POOL
    temp<-rmultinom(1, size=2*N, prob=c(gameteAB, gameteAb, gameteaB, gameteab))
    pAB[i+1]<-temp[1]/(2*N)
    pAb[i+1]<-temp[2]/(2*N)
    paB[i+1]<-temp[3]/(2*N)
    pab[i+1]<-temp[4]/(2*N)
    # FINALLY CALCULATE THE D IN THE NEXT GENERATION
    pA[i+1] <- pAB[i+1] + pAb[i+1]
    pB[i+1] <- pAB[i+1] + paB[i+1]
    D[i+1]<-(pAB[i+1]*pab[i+1]-pAb[i+1]*paB[i+1])
    #if ((pA[i+1]*(1-pA[i+1])* pB[i+1]*(1-pB[i+1])) == 0) print(pA[i+1])
    r[i+1]<-(pAB[i+1]*pab[i+1]-pAb[i+1]*paB[i+1])**2 / (pA[i+1]*(1-pA[i+1])* pB[i+1]*(1-pB[i+1]))
  }
  # RETURN A LIST OF THINGS
  return(list(pAB=pAB, pAb=pAb, paB=paB, pab=pab, D=D, pA=pA, pB=pB, r=r))
}


# many pairs of loci (four halpotypes)
sim.drift5<-function(p0=c(0.25, 0.25, 0.25, 0.25), N=100, c=0.5, t=10, loci = 1)
{
  pAB<-array(NA, dim = c(loci, t+1))
  pAb<-array(NA, dim = c(loci, t+1))
  paB<-array(NA, dim = c(loci, t+1))
  pab<-array(NA, dim = c(loci, t+1))
  pA<-array(NA, dim = c(loci, t+1))
  pB<-array(NA, dim = c(loci, t+1))
  D<-array(NA, dim = c(loci, t+1))
  r<-array(NA, dim = c(loci, t+1))
  # FILL THESE VECTORS WITH INITIAL VALUES
  pAB[,1]<-p0[1]
  pAb[,1]<-p0[2]
  paB[,1]<-p0[3]
  pab[,1]<-p0[4]
  pA[,1] <- pAB[,1] + pAb[,1]
  pB[,1] <- pAB[,1] + paB[,1]
  D[,1]<-(pAB[,1]*pab[,1]-pAb[,1]*paB[,1])
  r[,1]<-(pAB[,1]*pab[,1]-pAb[,1]*paB[,1])**2 / (pA[,1]*(1-pA[,1])* pB[,1]*(1-pB[,1]))
  # PROPAGATION. NEED TO CONSIDER BOTH DRIFT AND RECOMBINATION
  for (i in 1:t)
  {
    # RECOMBINATION
    gameteAB<-pAB[,i]-c*D[,i]
    gameteAb<-pAb[,i]+c*D[,i]
    gameteaB<-paB[,i]+c*D[,i]
    gameteab<-pab[,i]-c*D[,i]
    # DRIFT, RANDOM SAMPLING FROM THE GAMETE POOL
    temp<-rmultinomial(n=loci, size=2*N, prob=cbind(gameteAB, gameteAb, gameteaB, gameteab))
    pAB[,i+1]<-temp[,1]/(2*N)
    pAb[,i+1]<-temp[,2]/(2*N)
    paB[,i+1]<-temp[,3]/(2*N)
    pab[,i+1]<-temp[,4]/(2*N)
    # FINALLY CALCULATE THE D IN THE NEXT GENERATION
    pA[,i+1] <- pAB[,i+1] + pAb[,i+1]
    pB[,i+1] <- pAB[,i+1] + paB[,i+1]
    D[,i+1]<-(pAB[,i+1]*pab[,i+1]-pAb[,i+1]*paB[,i+1])
    r[,i+1]<-(pAB[,i+1]*pab[,i+1]-pAb[,i+1]*paB[,i+1])**2 / (pA[,i+1]*(1-pA[,i+1])*pB[,i+1]*(1-pB[,i+1]))
    
  }
  # RETURN A LIST OF THINGS
  return(list(pAB=pAB, pAb=pAb, paB=paB, pab=pab, D=D, pA=pA, pB=pB, r=r))
}


# many pairs of loci (four halpotypes), many subpopulations
sim.drift5_mig <- function(p0=c(0.25, 0.25, 0.25, 0.25), N_v=c(100,100), M = array(c(.9,.1,.1,.9), dim = c(2,2)), c=0.5, t=10, loci = 1, burn_in = 0)
{
  if (!check_M(M)) {
    print("Bad M")
    return("Bad M")
  }
  # pre migration statistics
  pAB_pre<-array(NA, dim = c(loci, t+1, length(N_v)))
  pAb_pre<-array(NA, dim = c(loci, t+1, length(N_v)))
  paB_pre<-array(NA, dim = c(loci, t+1, length(N_v)))
  pab_pre<-array(NA, dim = c(loci, t+1, length(N_v)))
  pA_pre<-array(NA, dim = c(loci, t+1, length(N_v)))
  pB_pre<-array(NA, dim = c(loci, t+1, length(N_v)))
  D_pre<-array(NA, dim = c(loci, t+1, length(N_v)))
  r_pre<-array(NA, dim = c(loci, t+1, length(N_v)))
  
  
  # post migration statistics
  pAB_post<-array(NA, dim = c(loci, t+1, length(N_v)))
  pAb_post<-array(NA, dim = c(loci, t+1, length(N_v)))
  paB_post<-array(NA, dim = c(loci, t+1, length(N_v)))
  pab_post<-array(NA, dim = c(loci, t+1, length(N_v)))
  pA_post<-array(NA, dim = c(loci, t+1, length(N_v)))
  pB_post<-array(NA, dim = c(loci, t+1, length(N_v)))
  D_post<-array(NA, dim = c(loci, t+1, length(N_v)))
  r_post<-array(NA, dim = c(loci, t+1, length(N_v)))
  
  # FILL THESE VECTORS WITH INITIAL VALUES
  if (is.vector(p0)){
    pAB_pre[,1,]<-p0[1]
    pAb_pre[,1,]<-p0[2]
    paB_pre[,1,]<-p0[3]
    pab_pre[,1,]<-p0[4]
  }
  if (is.array(p0)){
    if (dim(p0)[2] != length(N_v)) print("????????????????????????")
    for (sp in 1:length(N_v)){
      pAB_pre[,1,sp]<-p0[1,sp]
      pAb_pre[,1,sp]<-p0[2,sp]
      paB_pre[,1,sp]<-p0[3,sp]
      pab_pre[,1,sp]<-p0[4,sp]
    }
    
  }
  pA_pre[,1,] <- pAB_pre[,1,] + pAb_pre[,1,]
  pB_pre[,1,] <- pAB_pre[,1,] + paB_pre[,1,]
  D_pre[,1,]<-(pAB_pre[,1,]*pab_pre[,1,]-pAb_pre[,1,]*paB_pre[,1,])
  r_pre[,1,]<-D_pre[,1,]**2 / (pA_pre[,1,]*(1-pA_pre[,1,])* pB_pre[,1,]*(1-pB_pre[,1,]))
  
  
  for (i in 1:t)
  {
    
    # Migration
    pAB_post[,i,] <- pAB_pre[,i,] %*% M
    pAb_post[,i,] <- pAb_pre[,i,] %*% M
    paB_post[,i,] <- paB_pre[,i,] %*% M
    pab_post[,i,] <- pab_pre[,i,] %*% M
    
    # Calculate statistics

    pA_post[,i,] <- pAB_post[,i,] + pAb_post[,i,]
    pB_post[,i,] <- pAB_post[,i,] + paB_post[,i,]
    D_post[,i,] <- (pAB_post[,i,] * pab_post[,i,] - pAb_post[,i,] * paB_post[,i,])
    r_post[,i,] <- D_post[,i,]**2 / (pA_post[,i,]*(1-pA_post[,i,])*pB_post[,i,]*(1-pB_post[,i,]))
    
    # Recombine - make gamete pool
    gameteAB<-pAB_post[,i,]-c*D_post[,i,]
    gameteAb<-pAb_post[,i,]+c*D_post[,i,]
    gameteaB<-paB_post[,i,]+c*D_post[,i,]
    gameteab<-pab_post[,i,]-c*D_post[,i,]
    
    # Random sampling from gamete pool
    for (subpop in 1:length(N_v)){
      N_temp <- N_v[subpop]
      temp<-rmultinomial(n=loci, size=2*N_temp, prob=cbind(gameteAB[,subpop], gameteAb[,subpop], gameteaB[,subpop], gameteab[,subpop]))
      pAB_pre[,i+1,subpop]<-temp[,1]/(2*N_temp)
      pAb_pre[,i+1,subpop]<-temp[,2]/(2*N_temp)
      paB_pre[,i+1,subpop]<-temp[,3]/(2*N_temp)
      pab_pre[,i+1,subpop]<-temp[,4]/(2*N_temp)
    }
    
    # Calculate statistics
    pA_pre[,i+1,] <- pAB_pre[,i+1,] + pAb_pre[,i+1,]
    pB_pre[,i+1,] <- pAB_pre[,i+1,] + paB_pre[,i+1,]
    D_pre[,i+1,] <- (pAB_pre[,i+1,] * pab_pre[,i+1,] - pAb_pre[,i+1,] * paB_pre[,i+1,])
    r_pre[,i+1,] <- D_pre[,i+1,]**2 / (pA_pre[,i+1,]*(1-pA_pre[,i+1,])*pB_pre[,i+1,]*(1-pB_pre[,i+1,]))

  }
  # RETURN A LIST OF THINGS
  return(list(pAB_pre=pAB_pre, pAb_pre=pAb_pre, paB_pre=paB_pre, pab_pre=pab_pre, pAB_post=pAB_post, pAb_post=pAb_post, paB_post=paB_post, pab_post=pab_post, D_pre=D_pre, D_post=D_post, pA_pre=pA_pre, pA_post=pA_post, pB_pre=pB_pre, pB_post=pB_post, r_pre=r_pre, r_post = r_post))
}

propegate_population_old <- function(population_list, N_v = c(1000, 1000), M = array(c(.9,.1,.1,.9), dim = c(2,2)), c=0.5, t_add = 100){
  
  if (!check_M(M)) {
    print("Bad M")
    return("Bad M")
  }
  
  loci = population_list$loci
  last_pre = population_list$t_curr
  n_subpops = population_list$n_subpops
  
  haplotype_add <- array(NA, dim = c(loci, n_subpops, t_add))
  pAB_pre = abind(population_list$pAB_pre, haplotype_add, along = 3)
  pAb_pre = abind(population_list$pAb_pre, haplotype_add, along = 3)
  paB_pre = abind(population_list$paB_pre, haplotype_add, along = 3)
  pab_pre = abind(population_list$pab_pre, haplotype_add, along = 3)
  
  pAB_post = abind(population_list$pAB_post, haplotype_add, along = 3)
  pAb_post = abind(population_list$pAb_post, haplotype_add, along = 3)
  paB_post = abind(population_list$paB_post, haplotype_add, along = 3)
  pab_post = abind(population_list$pab_post, haplotype_add, along = 3)
  
  D_pre = abind(population_list$D_pre, haplotype_add, along = 3)
  r_pre = abind(population_list$r_pre, haplotype_add, along = 3)
  pA_pre = abind(population_list$pA_pre, haplotype_add, along = 3)
  pB_pre = abind(population_list$pB_pre, haplotype_add, along = 3)
  
  D_post = abind(population_list$D_post, haplotype_add, along = 3)
  r_post = abind(population_list$r_post, haplotype_add, along = 3)
  pA_post = abind(population_list$pA_post, haplotype_add, along = 3)
  pB_post = abind(population_list$pB_post, haplotype_add, along = 3)
  
  
  
  
  for (i in (last_pre):(last_pre+t_add-1))
  {
    
    # cat('\n\n\n')
    # print('i')
    # print(i)
    # cat('\n\n\n')
    # print('pAb_pre')
    # print(pAb_pre)
    
    ##############################################################
    ##############################################################
    ##############################################################
    ############## I am here in swapping dimensions ##############
    ##############################################################
    ##############################################################
    ##############################################################
    
    # Migration
    pAB_post[,,i] <- pAB_pre[,,i] %*% M
    pAb_post[,,i] <- pAb_pre[,,i] %*% M
    paB_post[,,i] <- paB_pre[,,i] %*% M
    pab_post[,,i] <- pab_pre[,,i] %*% M
    # print('M')
    # print(M)
    # print('pAb_post')
    # print(pAb_post)
    # 
    # # Calculate statistics
    # print(i)
    # print(dim(pA_post))
    # print(dim(pB_post))
    # print(dim(D_post))
    # print(dim(r_post))
    pA_post[,,i] <- pAB_post[,,i] + pAb_post[,,i]
    pB_post[,,i] <- pAB_post[,,i] + paB_post[,,i]
    D_post[,,i] <- (pAB_post[,,i] * pab_post[,,i] - pAb_post[,,i] * paB_post[,,i])
    r_post[,,i] <- D_post[,,i]**2 / (pA_post[,,i]*(1-pA_post[,,i])*pB_post[,,i]*(1-pB_post[,,i]))
    
    # Recombine - make gamete pool
    gameteAB<-pAB_post[,,i]-c*D_post[,,i]
    gameteAb<-pAb_post[,,i]+c*D_post[,,i]
    gameteaB<-paB_post[,,i]+c*D_post[,,i]
    gameteab<-pab_post[,,i]-c*D_post[,,i]
    # print('D_post')
    # print(D_post)
    # print('gameteAB')
    # print(gameteAb)
    
    # Random sampling from gamete pool
    for (subpop in 1:n_subpops){
      N_temp <- N_v[subpop]
      temp<-rmultinomial(n=loci, size=2*N_temp, prob=cbind(gameteAB[,subpop], gameteAb[,subpop], gameteaB[,subpop], gameteab[,subpop]))
      pAB_pre[,subpop, i+1]<-temp[,1]/(2*N_temp)
      pAb_pre[,subpop, i+1]<-temp[,2]/(2*N_temp)
      paB_pre[,subpop, i+1]<-temp[,3]/(2*N_temp)
      pab_pre[,subpop, i+1]<-temp[,4]/(2*N_temp)
    }
    
    # Calculate statistics
    #(dim(pA_pre))
    pA_pre[,,i+1] <- pAB_pre[,,i+1] + pAb_pre[,,i+1]
    pB_pre[,,i+1] <- pAB_pre[,,i+1] + paB_pre[,,i+1]
    D_pre[,,i+1] <- (pAB_pre[,,i+1] * pab_pre[,,i+1] - pAb_pre[,,i+1] * paB_pre[,,i+1])
    r_pre[,,i+1] <- D_pre[,,i+1]**2 / (pA_pre[,,i+1]*(1-pA_pre[,,i+1])*pB_pre[,,i+1]*(1-pB_pre[,,i+1]))
  }
  return(list(pAB_pre=pAB_pre, pAb_pre=pAb_pre, paB_pre=paB_pre, pab_pre=pab_pre, pAB_post=pAB_post, pAb_post=pAb_post, paB_post=paB_post, pab_post=pab_post, D_pre=D_pre, D_post=D_post, pA_pre=pA_pre, pA_post=pA_post, pB_pre=pB_pre, pB_post=pB_post, r_pre=r_pre, r_post = r_post, loci = loci, t_tot = NA, t_curr = last_pre+t_add, n_subpops = n_subpops))
}


### separated initialising and propegating ##

initialise_population <- function(p0=c(0.25, 0.25, 0.25, 0.25), n_subpops = 2, loci = 2, t_alloc = 1){
  # pre migration statistics
  pAB_pre<-array(NA, dim = c(loci, n_subpops, t_alloc))
  pAb_pre<-array(NA, dim = c(loci, n_subpops, t_alloc))
  paB_pre<-array(NA, dim = c(loci, n_subpops, t_alloc))
  pab_pre<-array(NA, dim = c(loci, n_subpops, t_alloc))
  pA_pre<-array(NA, dim = c(loci, n_subpops, t_alloc))
  pB_pre<-array(NA, dim = c(loci, n_subpops, t_alloc))
  D_pre<-array(NA, dim = c(loci, n_subpops, t_alloc))
  r_pre<-array(NA, dim = c(loci, n_subpops, t_alloc))
  
  
  # post migration statistics
  pAB_post<-array(NA, dim = c(loci, n_subpops, t_alloc-1))
  pAb_post<-array(NA, dim = c(loci, n_subpops, t_alloc-1))
  paB_post<-array(NA, dim = c(loci, n_subpops, t_alloc-1))
  pab_post<-array(NA, dim = c(loci, n_subpops, t_alloc-1))
  pA_post<-array(NA, dim = c(loci, n_subpops, t_alloc-1))
  pB_post<-array(NA, dim = c(loci, n_subpops, t_alloc-1))
  D_post<-array(NA, dim = c(loci, n_subpops, t_alloc-1))
  r_post<-array(NA, dim = c(loci, n_subpops, t_alloc-1))
  
  # FILL THESE VECTORS WITH INITIAL VALUES
  if (is.vector(p0)){
    pAB_pre[,,1]<-p0[1]
    pAb_pre[,,1]<-p0[2]
    paB_pre[,,1]<-p0[3]
    pab_pre[,,1]<-p0[4]
  }
  if (is.array(p0)){
    if (dim(p0)[2] != length(N_v)) print("????????????????????????")
    for (sp in 1:length(N_v)){
      pAB_pre[,sp,1]<-p0[1,sp]
      pAb_pre[,sp,1]<-p0[2,sp]
      paB_pre[,sp,1]<-p0[3,sp]
      pab_pre[,sp,1]<-p0[4,sp]
    }
    
  }
  pA_pre[,,1] <- pAB_pre[,,1] + pAb_pre[,,1]
  pB_pre[,,1] <- pAB_pre[,,1] + paB_pre[,,1]
  D_pre[,,1]<-(pAB_pre[,,1]*pab_pre[,,1]-pAb_pre[,,1]*paB_pre[,,1])
  r_pre[,,1]<-D_pre[,,1]**2 / (pA_pre[,,1]*(1-pA_pre[,,1])* pB_pre[,,1]*(1-pB_pre[,,1]))
  
  return(list(pAB_pre=pAB_pre, pAb_pre=pAb_pre, paB_pre=paB_pre, pab_pre=pab_pre, pAB_post=pAB_post, pAb_post=pAb_post, paB_post=paB_post, pab_post=pab_post, D_pre=D_pre, D_post=D_post, pA_pre=pA_pre, pA_post=pA_post, pB_pre=pB_pre, pB_post=pB_post, r_pre=r_pre, r_post = r_post, loci = loci, t_tot = NA, t_curr = 1, n_subpops = n_subpops, age = 0))
}

propegate_population_inner <- function(pAB_pre_first, pAb_pre_first, paB_pre_first, pab_pre_first, N_v = c(1000, 1000), M = array(c(.9,.1,.1,.9), dim = c(2,2)), c=0.5,loci, n_subpops,t_curr, t_add = 100){
  
  if (!check_M(M)) {
    print("Bad M")
    return("Bad M")
  }
  
  
  
  
  
  out_pre <- array(NA, dim = c(loci,  n_subpops, t_add+1) )
  pAB_pre = out_pre
  pAb_pre = out_pre
  paB_pre = out_pre
  pab_pre = out_pre
  
  out_post <- array(NA, dim = c(loci,  n_subpops, t_add) )
  pAB_post = out_post
  pAb_post = out_post
  paB_post = out_post
  pab_post = out_post
  
  r_post = out_post
  D_post = out_post
  pA_post = out_post
  pB_post = out_post
  
  
  
  pAB_pre[,,1] = pAB_pre_first
  pAb_pre[,,1] = pAb_pre_first
  paB_pre[,,1] = paB_pre_first
  pab_pre[,,1] = pab_pre_first
  
  
  
  for (i in (1):(t_add))
  {
    # Migration
    pAB_post[,,i] <- pAB_pre[,,i] %*% M
    pAb_post[,,i] <- pAb_pre[,,i] %*% M
    paB_post[,,i] <- paB_pre[,,i] %*% M
    pab_post[,,i] <- pab_pre[,,i] %*% M
    # 
    # # Calculate statistics
    # print(i)
    # print(dim(pA_post))
    # print(dim(pB_post))
    # print(dim(D_post))
    # print(dim(r_post))
    pA_post[,,i] <- pAB_post[,,i] + pAb_post[,,i]
    pB_post[,,i] <- pAB_post[,,i] + paB_post[,,i]
    D_post[,,i] <- (pAB_post[,,i] * pab_post[,,i] - pAb_post[,,i] * paB_post[,,i])
    r_post[,,i] <- D_post[,,i]**2 / (pA_post[,,i]*(1-pA_post[,,i])*pB_post[,,i]*(1-pB_post[,,i]))
    
    # Recombine - make gamete pool
    #print(D_post[,,])
    gameteAB<-pAB_post[,,i]-c*D_post[,,i]
    gameteAb<-pAb_post[,,i]+c*D_post[,,i]
    gameteaB<-paB_post[,,i]+c*D_post[,,i]
    gameteab<-pab_post[,,i]-c*D_post[,,i]
    
    # Random sampling from gamete pool
    for (subpop in 1:n_subpops){
      N_temp <- N_v[subpop]
      temp<-rmultinomial(n=loci, size=2*N_temp, prob=cbind(gameteAB[,subpop], gameteAb[,subpop], gameteaB[,subpop], gameteab[,subpop]))
      pAB_pre[,subpop,i+1]<-temp[,1]/(2*N_temp)
      pAb_pre[,subpop,i+1]<-temp[,2]/(2*N_temp)
      paB_pre[,subpop,i+1]<-temp[,3]/(2*N_temp)
      pab_pre[,subpop,i+1]<-temp[,4]/(2*N_temp)
    }
    
    # Calculate statistics
    #pA_pre[,i+1,] <- pAB_pre[,i+1,] + pAb_pre[,i+1,]
    #pB_pre[,i+1,] <- pAB_pre[,i+1,] + paB_pre[,i+1,]
    # D_pre[,i+1,] <- (pAB_pre[,i+1,] * pab_pre[,i+1,] - pAb_pre[,i+1,] * paB_pre[,i+1,])
    # r_pre[,i+1,] <- D_pre[,i+1,]**2 / (pA_pre[,i+1,]*(1-pA_pre[,i+1,])*pB_pre[,i+1,]*(1-pB_pre[,i+1,]))
  }
  return(list(pAB_pre=pAB_pre, pAb_pre=pAb_pre, paB_pre=paB_pre, pab_pre=pab_pre, pAB_post=pAB_post, pAb_post=pAb_post, paB_post=paB_post, pab_post=pab_post, D_post=D_post,r_post = r_post, pA_post=pA_post, pB_post=pB_post,loci = loci, t_tot = NA, t_curr = 1+t_add, n_subpops = n_subpops))
}

propegate_population_decomposed <- function(population_list, N_v = c(1000, 1000), M = array(c(.9,.1,.1,.9), dim = c(2,2)), c=0.5, t_add = 100){
  req <- get_req_info(population_list)
  return(propegate_population_inner(req$pAB_pre, req$pAb_pre, req$paB_pre, req$pab_pre, N_v, M, c, req$loci, req$n_subpops,req$t_curr,  t_add, req$age ))
}

get_req_info = function(pop){
  t_curr = pop$t_curr
  
  return( list(pAB_pre = pop$pAB_pre[,,t_curr], pAb_pre = pop$pAb_pre[,,t_curr], paB_pre = pop$paB_pre[,,t_curr],pab_pre = pop$pab_pre[,,t_curr], t_curr = pop$t_curr, n_subpops = as.integer(pop$n_subpops), loci = pop$loci, age = pop$age) )
}


### getting population to equilibrium ###

check_eqm <- function(pop, N_v, fst_v = c(), t_change = 1, window_size = 40, direction = FALSE){
  #print('1a')
  t_final <- pop$t_curr
  # print('2a')
  
  t_past  <- length(fst_v)
  # print('3a')
  
  if (t_final-t_change < window_size) return(0)
  #  print('4a')
  
  t_change = t_change + t_past
  #  print('5a')
  
  # print(fst_v)
  # print('pA')
  # print(pop$pA_post[,,t_final-1])
  # print('pB')
  # print(pop$pB_post[,,t_final-1])
  # 
  # print('fst raw')
  # print( get_Fst_vec(pop$pA_post[,,1:t_final-1], N_v, t_final-1) )
  
  pA_means = apply(pop$pA_post, c(1,3), mean)
  pA_filt = apply(pA_means, 1, function(r){any(r %in% c(0,1))} )
  pA_final = pop$pA_post[!pA_filt,,]
  
  #print(pop$pB_post)
  pB_means = apply(pop$pB_post, c(1,3), mean)
  #print(pB_means)
  pB_filt = apply(pB_means, 1, function(r) any(r %in% c(0,1)))
  #print(pB_filt)
  pB_final = pop$pB_post[!pB_filt,,]
  #print(pB_final)
  #print(1:t_final-1)
  
  fst_v = c(fst_v, colMeans(rbind(get_Fst_vec(pA_final[,,1:t_final-1], N_v, t_final-1), get_Fst_vec(pB_final[,,1:t_final-1], N_v, t_final-1))) )
  # print('6a')
  
  # print(fst_v)
  diff_fst = diff(fst_v)
  # print('7a')
  
  ma_diff_fst <- ma(diff_fst, n = window_size)
  #print('8a')
  
  # par(mfrow = c(3,1))
  # plot(1:length(fst_v), fst_v, type = 'l')
  # plot(1:length(diff_fst), diff_fst, type = 'l')
  # plot(1:length(ma_diff_fst), ma_diff_fst, type = 'l', ylim = c(-mean(ma_diff_fst, na.rm=T), max(ma_diff_fst, na.rm = T))); abline(h=0)
  # 
  #print(mean(diff_fst[t_change:(t_change+10)], na.rm = T))
  #  print(direction)
  # print(fst_v)
  # print(diff_fst)
  if (Inf %in% fst_v){
    #  print(pop$pA_post)
    #   print(pA_means)
    # print(pop$pB_post)
    # print(pB_means)
    
    #  print(pA_filt)
    #  print(pB_filt)
    
    #  print(pA_final)
    # print(pB_final)
    
  }
  #  print(mean(diff_fst[t_change:(t_change+10)]))
  
  if (direction == "none"){
    if (mean(diff_fst[t_change:(t_change+10)], na.rm = T) > 0) {direction = 0.5}
    else if (mean(diff_fst[t_change:(t_change+10)], na.rm = T) < 0) {direction = -0.5}
    else{print("how..."); print(fst_v); print(diff_fst)}
  }
  # print('9a')
  # print(direction)
  # print(ma_diff_fst[1:(length(fst_v)-1-window_size)])
  if ( all(direction * ma_diff_fst[1:(length(fst_v)-1-window_size)] > 0)) {return(c(direction, fst_v))}
  else {
    eqm = max(which.max(direction * ma_diff_fst[1:length(ma_diff_fst)] < 0) - t_past, 1)
    return (c(eqm, NA))
  }
}

do_burn_in <- function(pop, N_v, M, recom_rate, step_size = 50, window_size = 40){
  fst_v <- c()
  #print(111)
  pop <- propegate_population(pop, N_v = N_v, M = M, c = recom_rate, t_add = step_size)
  # print(222)
  # print(pop$pA_post[,,1])
  o = check_eqm(pop, N_v = N_v, direction = "none", window_size = window_size)
  # print(333)
  direction_or_eqm = o[1]
  fst_v = o[2:(step_size+1)]
  # print(444)
  
  if (direction_or_eqm > 0.7) return (pop)
  #print(555)
  
  cont = T
  while(cont){
    #print(666)
    
    pop = propegate_population(pop, N_v = N_v, M = M, c = recom_rate, t_add = step_size)
    
    #print(777)
    
    o = check_eqm(pop, N_v = N_v, direction = direction_or_eqm, fst_v = fst_v, window_size = window_size)
    #print(888)
    
    direction_or_eqm = o[1]
    #print(999)
    
    if (direction_or_eqm > 0.7){
      if (pop$t_curr - direction_or_eqm < 43){
        #print(1010101)
        
        pop <- propegate_population(pop, N_v = N_v, M = M, c = recom_rate, t_add = pop$t_curr - direction_or_eqm )
      }
      
      #print(111111)
      return (pop)
    } 
    
    fst_v <- o[(step_size+2) : (2*step_size+1)]
  }
}

check_eqm_r2 <- function(pop, N_v, r2v = c(), t_change = 1, window_size = 40, direction = FALSE){
  
  t_final <- pop$t_curr
  t_past  <- length(r2v)
  
  if (t_final-t_change < window_size) return(0)
  t_change = t_change + t_past
  
  
  
  r2v = c(r2v, apply(pop$r_post, 3, mean, na.rm = T) )
  diff_r2 = diff(r2v)
  ma_diff_r2 <- ma(diff_r2, n = window_size)
  # par(mfrow = c(3,1))
  # plot(1:length(r2v), r2v, type = 'l')
  # plot(1:length(diff_r2), diff_r2, type = 'l')
  # plot(1:length(ma_diff_r2), ma_diff_r2, type = 'l', ylim = c(-mean(ma_diff_r2, na.rm=T), max(ma_diff_r2, na.rm = T))); abline(h=0)
  #  
  #print(mean(diff_fst[t_change:(t_change+10)], na.rm = T))
  
  
  if (direction == "none"){
    if (mean(diff_r2[t_change:(t_change+3)], na.rm = T) > 0) {direction = 0.5}
    if (mean(diff_r2[t_change:(t_change+3)], na.rm = T) < 0) {direction = -0.5}
  }
  if ( all(direction * ma_diff_r2[1:(length(r2v)-1-window_size)] > 0)) {return(c(direction, r2v))}
  else {
    eqm = max(which.max(direction * ma_diff_r2[1:length(ma_diff_r2)] < 0) - t_past, 1)
    return (c(eqm, NA))
  }
}

do_burn_in_r2 <- function(pop, N_v, M, recom_rate, step_size = 8, window_size = 5){
  
  return(propegate_population(pop, N_v = N_v, M = M, c = recom_rate, t_add = 40))
  
  r2v <- c()
  pop <- propegate_population(pop, N_v = N_v, M = M, c = recom_rate, t_add = step_size)
  o = check_eqm_r2(pop, N_v = N_v, direction = "none", window_size = window_size)
  direction_or_eqm = o[1]
  r2v = o[2:(step_size+1)]
  
  
  if (direction_or_eqm > 0.7) return (pop)
  
  cont = T
  while(cont){
    pop = propegate_population(pop, N_v = N_v, M = M, c = recom_rate, t_add = step_size)
    o = check_eqm_r2(pop, N_v = N_v, direction = direction_or_eqm, r2v = r2v, window_size = window_size)
    direction_or_eqm = o[1]
    
    
    if (direction_or_eqm > 0.7) return (pop)
    
    r2v <- o[(step_size+2) : (2*step_size+1)]
    # print(fst_v)
    # return()
    # print(length(fst_v))
    
  }
}


#### get specified summary statistics ###
get_SS <- function(pop, N, m, Mfn, sampling_intervals, recom_rates, S_vec, grid_size, sample_subpops, max_loci, min_loci, n_loci, p =T){
  sample_gens = c(0, sampling_intervals)
  
  if (p) cat('N = ', N, ', m = ', m, '\n')
  
  N_v = rep(N, grid_size);# N_v = c(N, 1e7)
  M = Mfn(m, grid_size)
  
  Fc_mat <- array(NA, dim = c( length(recom_rates), length(S_vec), length(sampling_intervals)) )
  r2_mat <- array(NA, dim = c( length(recom_rates), length(S_vec), length(sampling_intervals) + 1) )
  first_burn = T
  
  idx_c = 1; 
  # choose c
  for (recom_rate in recom_rates){
    idx_S = 1;
    
    if (p) cat('\n\tc = ', recom_rate, '\n')
    
    # burn in
    if (m != 0 & first_burn) {cat('\t\tfst burn\n'); pop <- do_burn_in(pop, N_v, M, recom_rate); first_burn = F}
    else if (m != 0 & !first_burn) {cat('\t\tr2 burn\n');pop <- do_burn_in_r2(pop, N_v, M, recom_rate)}
    else if (m == 0){ pop <- propegate_population(pop, N_v, M, recom_rate, t_add = 30) }
    
    if (length(get_loci_remaining(pop)) < min_loci) {
      if (p) cat('\tReplacing population\n')
      pop <- initialise_population(n_subpops = grid_size, loci = max_loci, t_alloc = 1)
      pop <- do_burn_in(pop, N_v, M, recom_rate)
    }
    
    if (p) cat(paste0('\tAge: ', as.character(pop$age)))
    if (p) cat(',  Loci remaining: ', length(get_loci_remaining(pop)), '\n')
    
    sample_loci <- sample(get_loci_remaining(pop), n_loci, replace = F)
    
    # take samples & calculate SS
    final_gen = pop$t_curr - 1
    sp = sample_subpops[1]
    
    cat('\tSampling and calculating\n')
    
    #print(0)
    for (S in S_vec){
      #print(1)
      
      idx_g_Fc = 1; idx_g_r2 = 1
      #print(2)
      sMax <- take_sample_r(pop$pAB_post[sample_loci,sp, final_gen], pop$pAb_post[sample_loci,sp, final_gen], pop$paB_post[sample_loci,sp, final_gen], pop$pab_post[sample_loci,sp, final_gen], S)
      #print(3)
      
      
      r2_mat[idx_c, idx_S, idx_g_r2] = get_r2_drift(sMax$pAB, sMax$pAb, sMax$paB, sMax$pab, S = S, c = c); idx_g_r2 = idx_g_r2 + 1
      #print(4)
      for (g in sampling_intervals){
        #print(5)
        sx <- take_sample_r(pop$pAB_post[sample_loci,sp,final_gen - g], pop$pAb_post[sample_loci,sp,final_gen - g], pop$paB_post[sample_loci,sp, final_gen - g], pop$pab_post[sample_loci,sp,final_gen - g], S)
        
        #print(6)
        Fc_i <- (get_mean_fc(sMax$pAB+sMax$pAb, sx$pAB+sx$pAb, S, S, g) + get_mean_fc(sMax$pAB+sMax$pAb, sx$pAB+sx$pAb, S, S, g)) / 2
        #print(7)
        r2_mat[idx_c, idx_S, idx_g_r2] = get_r2_drift(sx$pAB, sx$pAb, sx$paB, sx$pab, S = S, c = c)
        #print(8)
        Fc_mat[idx_c, idx_S, idx_g_Fc] =  Fc_i
        #print(9)
        idx_g_Fc = idx_g_Fc + 1; idx_g_r2 = idx_g_r2 + 1
      }
      idx_S = idx_S + 1
    }
    idx_c = idx_c + 1
  }
  #rint(10)
  
  dimnames(Fc_mat)[[1]] = paste0('c', recom_rates); dimnames(r2_mat)[[1]] = paste0('c', recom_rates)
  dimnames(Fc_mat)[[2]] = paste0('S', S_vec); dimnames(r2_mat)[[2]] = paste0('S', S_vec)
  dimnames(Fc_mat)[[3]] = paste0('g', sampling_intervals); dimnames(r2_mat)[[3]] = paste0('g', sample_gens)
  
  #print(11)
  return(list(pop = pop, Fc = Fc_mat, r2 = r2_mat))
}

get_loci_remaining <- function(pop, threshold = 0.05){
  t_curr = pop$t_curr
  A=(rowMeans(pop$pA_post[,,t_curr-1]) - .5)**2 < (0.5 - threshold)**2
  B=(rowMeans(pop$pB_post[,,t_curr-1]) - .5)**2 < (0.5 - threshold)**2
  return((1:length(A))[A&B])
}

# get all possible time intervals
get_all_intervals = function(sample_gens){
  all_intervals = c()
  
  for (i in sample_gens){ for (j in sample_gens){all_intervals = c(all_intervals, i-j)} };
  all_intervals_clean = sort(unique(abs(all_intervals)))
  counts = table((all_intervals[all_intervals>=0]))
  
  return(all_intervals_clean[-1])
}

# get the number of times each time interval occurs
get_intervals_counts = function(sample_gens){
  all_intervals = c()
  
  for (i in sample_gens){ for (j in sample_gens){all_intervals = c(all_intervals, i-j)} };
  all_intervals_clean = sort(unique(abs(all_intervals)))
  counts = table((all_intervals[all_intervals>=0]))
  
  return(counts)
}

# get adjacencies from a list of subpopulations and the edgelength fo the sampled square
gen_adj <- function(ss, el){
  coords = array(NA, dim = c(length(ss), 3))
  coords[,1] = ss
  i=1
  for (subpop in ss){
    coords[i,2:3] = c(subpop %% el, (subpop%/%el)+1)
    i=i+1
  }
  
  adj_l = list(); for (d in 1:((el-1)*2)) {adj_l[[d]] = array(NA, dim = c(0,2))}
  for (i in 1:(length(ss)-1)){
    coord = coords[i,]
    for ( j in min((i+1), length(ss)):length(ss)){
      d = sum(abs(coords[i,2:3]-coords[j,2:3]))
      
      adj_l[[d]] <- rbind(adj_l[[d]], c(ss[i], ss[j]))
    }
  }
  adj_l_trim = list()
  for (d in 1:length(adj_l)){
    if (length(adj_l[[d]]) == 0){
      break
    }
    adj_l_trim[[d]] = adj_l[[d]]
  }
  return(adj_l_trim)
}

gen_1_adj = function(ss){
  idxs = 1:length(ss)
  edge_len = sqrt(length(ss))
  adj = array(NA, dim = c(0,2))
  for (idx in idxs){
    if ((idx - 1) %% edge_len != 0){
      print(idx); print(edge_len)
      adj = rbind(adj, c(ss[idx], ss[idx - 1]))
    }
    if ((idx) %% edge_len != 0){
      adj = rbind(adj, c(ss[idx], ss[idx + 1]))
    }
    if ((idx + edge_len) <= length(sample_subpops)){
      adj = rbind(adj, c(ss[idx], ss[idx + edge_len]))
    }
    if ((idx - edge_len) >= 1){
      adj = rbind(adj, c(ss[idx], ss[idx - edge_len]))
    }
  }
  
  return (adj)
  
}

pair_to_idx = function(n1, n2, nmax){
  return(n1+n2*nmax)
}


### get summary statistics when sampling from multiple subpopulations ###
get_SS_multi_pops <- function(pop, N, m, Mfn, sample_gens, recom_rates, S_list, grid_size, sample_subpops, max_loci, min_loci, n_loci){
  sample_gens
  all_intervals = get_all_intervals(sample_gens)
  
  adj = gen_adj(1:length(sample_subpops))
  print(adj)
  riji_dists = length(gen_adj(sample_subpops))
  
  cat('N = ', N, ', m = ', m, '\n')
  
  N_v = rep(N, grid_size);# N_v = c(N, 1e7)
  M = Mfn(m, grid_size)
  
  Fc_mat <- array(NA, dim = c( length(recom_rates), length(S_vec), length(all_intervals), length(sample_subpops)) )
  r2_mat <- array(NA, dim = c( length(recom_rates), length(S_vec), length(sampling_intervals) + 1, length(sample_subpops)) )
  riji_mat <- array(NA, dim = c( length(recom_rates), length(S_vec), length(sampling_intervals) + 1, riji_dists) )
  FcCor_mat <-  array(NA, dim = c( length(recom_rates), length(S_vec), length(all_intervals), riji_dists) )
  first_burn = T
  
  idx_c = 1; 
  # choose c
  for (recom_rate in recom_rates){
    idx_S = 1;
    
    cat('\n\tc = ', recom_rate, '\n')
    
    # burn in
    if (m != 0 & first_burn){
      cat('\t\tfst burn\n'); pop <- do_burn_in(pop, N_v, M, recom_rate); first_burn = F
    } else if (m != 0 & !first_burn) {
      cat('\t\tr2 burn\n');pop <- do_burn_in_r2(pop, N_v, M, recom_rate)
    } else if (m == 0){
      pop <- propegate_population(pop, N_v, M, recom_rate, t_add = 30)
    }
    
    if (length(get_loci_remaining(pop)) < min_loci) {
      cat('\tReplacing population\n')
      pop <- initialise_population(n_subpops = grid_size, loci = max_loci, t_alloc = 1)
      pop <- do_burn_in(pop, N_v, M, recom_rate)
    }
    
    cat(paste0('\tAge: ', as.character(pop$age)))
    cat(',  Loci remaining: ', length(get_loci_remaining(pop)), '\n')
    
    sample_loci <- sample(get_loci_remaining(pop), n_loci, replace = F)
    
    # take samples & calculate SS
    final_gen = pop$t_curr - 1
    sp = sample_subpops[1]
    
    cat('\tSampling and calculating\n')
    
    #print(0)
    for (S in S_vec){
      samples_by_t = list()
      idx_g = 1
      for (g in sample_gens){
        samples_by_pop = list()
        idx_sp = 1
        for (sp in sample_subpops){
          samples_by_pop[[idx_sp]] = take_sample_r(pop$pAB_post[sample_loci,sp, final_gen - g], pop$pAb_post[sample_loci,sp, final_gen - g], pop$paB_post[sample_loci,sp, final_gen - g], pop$pab_post[sample_loci,sp, final_gen - g], S)
          idx_sp = idx_sp + 1
        }
        samples_by_t[[idx_g]] = samples_by_pop
        idx_g = idx_g + 1
      }
      
      
      
      fcs = list()
      
      
      for (g_idx_1 in 1:length(sample_gens)){
        g1 = sample_gens[g_idx_1]
        
        s1s = samples_by_t[[g_idx_1]]
        
        r2s <- lapply(s1s, get_r2_drift_sample, S = S, c = c)
        r2_mat[idx_c, idx_S, g_idx_1, ] = unlist(r2s)
        
        
        for (d in 1:riji_dists){
          adj_d = adj[[d]]
          riji_mat[idx_c, idx_S, g_idx_1,d] = 0
          for (r in 1:dim(adj_d)[1]){
            r = adj_d[r,]
            riji_mat[idx_c, idx_S, g_idx_1,d] = riji_mat[idx_c, idx_S, g_idx_1,d] + mean(get_riji(s1s[[r[1]]], s1s[[r[2]]]))
          }
          riji_mat[idx_c, idx_S, g_idx_1,d] / dim(adj_d)[1]
        }
        
        
        
        
        for (g_idx_2 in (min(g_idx_1+1, length(sample_gens))): length(sample_gens)){
          g2 = sample_gens[g_idx_2]
          s2s = samples_by_t[[g_idx_2]]
          
          
          
          for (i in 1:length(s2s)) Fc_mat[idx_c, idx_S, which(all_intervals == g2 - g1),i ] = get_mean_fc_sample(s1s[[i]], s2s[[i]], S1 = S, S2 = S)
          
          
          for(d in 1:riji_dists){
            adj_d = adj[[d]]
            FcCor_mat[idx_c, idx_S,  which(all_intervals == g2 - g1),d] = 0
            for (r in 1:dim(adj_d)[1]){
              r = adj_d[r,]
              FcCor_mat[idx_c, idx_S,  which(all_intervals == g2 - g1),d] = FcCor_mat[idx_c, idx_S,  which(all_intervals == g2 - g1),d] + get_FcCor(s1s[[r[1]]], s2s[[r[1]]], s1s[[r[2]]], s2s[[r[2]]] , S_size = S)
            }
            FcCor_mat[idx_c, idx_S,  which(all_intervals == g_idx_2 - g_idx_1),d] = FcCor_mat[idx_c, idx_S,  which(all_intervals == g_idx_2 - g_idx_1),d] / dim(adj_d)[1]
          }
          
        }
        
      }
      idx_S = idx_S+1
    }
    idx_c = idx_c+1
  }
  
  
  
  print(all_intervals)
  #rint(10)
  
  dimnames(Fc_mat)[[1]] = paste0('c', recom_rates);         dimnames(r2_mat)[[1]] = paste0('c', recom_rates)
  dimnames(Fc_mat)[[2]] = paste0('S', S_vec);               dimnames(r2_mat)[[2]] = paste0('S', S_vec)
  dimnames(Fc_mat)[[3]] = paste0('g', all_intervals);       dimnames(r2_mat)[[3]] = paste0('g', sample_gens)
  dimnames(Fc_mat)[[4]] = paste0('sp', sample_subpops);     dimnames(r2_mat)[[4]] = paste0('sp', sample_subpops)
  
  
  dimnames(FcCor_mat)[[1]] = paste0('c', recom_rates);        dimnames(riji_mat)[[1]] = paste0('c', recom_rates)
  dimnames(FcCor_mat)[[2]] = paste0('S', S_vec);              dimnames(riji_mat)[[2]] = paste0('S', S_vec)
  dimnames(FcCor_mat)[[3]] = paste0('g', all_intervals);      dimnames(riji_mat)[[3]] = paste0('g', sample_gens)
  dimnames(FcCor_mat)[[4]] = paste0('sp', 1:riji_dists);        dimnames(riji_mat)[[4]] = paste0('sp', 1:riji_dists)
  
  #print(11)
  return(list(pop = pop, Fc = Fc_mat, r2 = r2_mat, riji = riji_mat, FcCor = FcCor_mat))
}

get_SS_multi_pops2 <- function(pop, N, m, Mfn, sample_gens, recom_rates, S_vec, grid_size, sample_subpops, max_loci, min_loci, n_loci, p = F){
  if (p) cat('N = ', N, ', m = ', m, '\n');
  
  N_v = rep(N, grid_size);# N_v = c(N, 1e7)
  M = Mfn(m, grid_size)
  Fc_mat <- array(NA, dim = c( length(S_vec), length(sample_gens), length(sample_gens), length(sample_subpops), length(recom_rates)) )
  r2_mat <- array(NA, dim = c( length(S_vec), length(sample_gens), length(recom_rates), length(sample_subpops)) )
  FcCor_mat <- array(NA, dim = c( length(S_vec), length(sample_gens), length(sample_gens), length(sample_subpops), length(sample_subpops), length(recom_rates)) )
  rirj_mat <- array(NA, dim = c( length(S_vec), length(sample_gens), length(recom_rates), length(sample_subpops), length(sample_subpops)) )
  first_burn = T
  idx_c = 1
  
  for (recom_rate in recom_rates){
    gc()
    if (p) cat('\n\tc = ', recom_rate, '\n')
    
    # burn in
    if (m != 0 & first_burn){
      if (p) cat('\t\tfst burn\n')
      pop <- do_burn_in(pop, N_v, M, recom_rate); first_burn = F
    } else if (m != 0 & !first_burn) {
      if (p) cat('\t\tr2 burn\n')
      pop <- do_burn_in_r2(pop, N_v, M, recom_rate)
    } else if (m == 0){
      pop <- propegate_population(pop, N_v, M, recom_rate, t_add = 40)
    }
    
    if (length(get_loci_remaining(pop)) < min_loci) {
      if (p) cat('\tReplacing population\n')
      pop <- initialise_population(n_subpops = grid_size, loci = max_loci, t_alloc = 1)
      pop <- do_burn_in(pop, N_v, M, recom_rate)
    }
    
    
    gc()
    
    if (p) cat(paste0('\tAge: ', as.character(pop$age)))
    if (p) cat(',  Loci remaining: ', length(get_loci_remaining(pop)), '\n')
    
    sample_loci <- sample(get_loci_remaining(pop), n_loci, replace = F)
    
    # take samples & calculate SS
    final_gen = pop$t_curr - 1
    sp = sample_subpops[1]
    
    if (p) cat('\tSampling and calculating\n')
    idx_S = 1;
    for (S in S_vec){
      samples_by_t = list()
      idx_g = 1
      for (g in sample_gens){
        samples_by_pop = list()
        idx_sp = 1
        for (sp in sample_subpops){
          samples_by_pop[[idx_sp]] = take_sample_r(pop$pAB_post[sample_loci,sp, final_gen - g], pop$pAb_post[sample_loci,sp, final_gen - g], pop$paB_post[sample_loci,sp, final_gen - g], pop$pab_post[sample_loci,sp, final_gen - g], S)
          idx_sp = idx_sp + 1
        }
        samples_by_t[[idx_g]] = samples_by_pop
        idx_g = idx_g + 1
      }
      
      for (idx_g_1 in 1:length(sample_gens)){
        g1 = sample_gens[idx_g_1]
        
        s1s = samples_by_t[[idx_g_1]]
        
        r2s <- lapply(s1s, get_r2_drift_sample, S = S, c = c)
        r2_mat[idx_S, idx_g_1, idx_c, ] = unlist(r2s)
        
        for (idx_subpop_1 in 1:length(sample_subpops)){
          for (idx_subpop_2 in 1:length(sample_subpops)){
            sp1 = sample_subpops[idx_subpop_1]; sp2 = sample_subpops[idx_subpop_2] 
            rirj_mat[idx_S, idx_g_1, idx_c, idx_subpop_1, idx_subpop_2] = mean(get_riji(s1s[[idx_subpop_1]], s1s[[idx_subpop_2]]))
          }
        }
        
        
        for (idx_g_2 in (min(idx_g_1+1, length(sample_gens))):length(sample_gens)){
          g2 = sample_gens[idx_g_2]
          s2s = samples_by_t[[idx_g_2]]
          
          for (idx_sp in 1:length(s2s)){
            Fc_mat[idx_S, idx_g_1, idx_g_2, idx_sp, idx_c ] = get_mean_fc_sample(s1s[[idx_sp]], s2s[[idx_sp]], S1 = S, S2 = S)
          }
          
          for (idx_subpop_1 in 1:length(sample_subpops)){
            for (idx_subpop_2 in 1:length(sample_subpops)){
              sp1 = sample_subpops[idx_subpop_1]; sp2 = sample_subpops[idx_subpop_2] 
              FcCor_mat[idx_S, idx_g_1, idx_g_2, idx_subpop_1, idx_subpop_2, idx_c] = get_FcCor( s1s[[idx_subpop_1]], s2s[[idx_subpop_1]], s1s[[idx_subpop_2]], s2s[[idx_subpop_2]] , S_size = S)
            }
          }
        }
      }
      idx_S = idx_S+1
    }
    idx_c = idx_c+1
  }
  
  
  
  dimnames(Fc_mat)[[1]] = paste0('S', S_vec);
  dimnames(Fc_mat)[[2]] = paste0('g', sample_gens);
  dimnames(Fc_mat)[[3]] = paste0('g', sample_gens);
  dimnames(Fc_mat)[[4]] = paste0('sp', sample_subpops);
  dimnames(Fc_mat)[[5]] = paste0('c', recom_rates);
  
  dimnames(r2_mat)[[1]] = paste0('S', S_vec);
  dimnames(r2_mat)[[2]] = paste0('g', sample_gens);
  dimnames(r2_mat)[[3]] = paste0('c', recom_rates);
  dimnames(r2_mat)[[4]] = paste0('sp', sample_subpops);
  
  dimnames(FcCor_mat)[[1]] = paste0('S', S_vec);
  dimnames(FcCor_mat)[[2]] = paste0('g', sample_gens);
  dimnames(FcCor_mat)[[3]] = paste0('g', sample_gens);
  dimnames(FcCor_mat)[[4]] = paste0('sp', sample_subpops);
  dimnames(FcCor_mat)[[5]] = paste0('sp', sample_subpops);
  dimnames(FcCor_mat)[[6]] = paste0('c', recom_rates);
  
  dimnames(rirj_mat)[[1]] = paste0('S', S_vec);
  dimnames(rirj_mat)[[2]] = paste0('g', sample_gens);
  dimnames(rirj_mat)[[3]] = paste0('c', recom_rates);
  dimnames(rirj_mat)[[4]] = paste0('sp', sample_subpops);
  dimnames(rirj_mat)[[5]] = paste0('sp', sample_subpops);
  
  
  #print(11)
  return(list(pop = pop, Fc = Fc_mat, r2 = r2_mat, rirj = rirj_mat, FcCor = FcCor_mat))
}

### ONLY TAKE CORRELATIONS WITH FOCAL POPULATION ! ### 
get_SS_multi_pops3 <- function(pop, N, m, Mfn, sample_gens, recom_rates, S_vec, grid_size, sample_subpops, max_loci, min_loci, n_loci, focal_pop, p = F){
  if (p) cat('N = ', N, ', m = ', m, '\n');
  
  N_v = rep(N, grid_size);# N_v = c(N, 1e7)
  M = Mfn(m, grid_size)
  Fc_mat <- array(NA, dim = c( length(S_vec), length(sample_gens), length(sample_gens), length(sample_subpops), length(recom_rates)) )
  r2_mat <- array(NA, dim = c( length(S_vec), length(sample_gens), length(recom_rates), length(sample_subpops)) )
  FcCor_mat <- array(NA, dim = c( length(S_vec), length(sample_gens), length(sample_gens), length(sample_subpops), length(sample_subpops), length(recom_rates)) )
  rirj_mat <- array(NA, dim = c( length(S_vec), length(sample_gens), length(recom_rates), length(sample_subpops), length(sample_subpops)) )
  first_burn = T
  idx_c = 1
  
  for (recom_rate in recom_rates){
    gc()
    if (p) cat('\n\tc = ', recom_rate, '\n')
    
    # burn in
    if (m != 0 & first_burn){
      if (p) cat('\t\tfst burn\n')
      pop <- do_burn_in(pop, N_v, M, recom_rate); first_burn = F
    } else if (m != 0 & !first_burn) {
      if (p) cat('\t\tr2 burn\n')
      pop <- do_burn_in_r2(pop, N_v, M, recom_rate)
    } else if (m == 0){
      pop <- propegate_population(pop, N_v, M, recom_rate, t_add = 30)
    }
    
    if (length(get_loci_remaining(pop)) < min_loci) {
      if (p) cat('\tReplacing population\n')
      pop <- initialise_population(n_subpops = grid_size, loci = max_loci, t_alloc = 1)
      pop <- do_burn_in(pop, N_v, M, recom_rate)
    }
    
    if (p) cat(paste0('\tAge: ', as.character(pop$age)))
    if (p) cat(',  Loci remaining: ', length(get_loci_remaining(pop)), '\n')
    
    sample_loci <- sample(get_loci_remaining(pop), n_loci, replace = F)
    
    # take samples & calculate SS
    final_gen = pop$t_curr - 1
    sp = sample_subpops[1]
    
    if (p) cat('\tSampling and calculating\n')
    idx_S = 1;
    for (S in S_vec){
      samples_by_t = list()
      idx_g = 1
      for (g in sample_gens){
        samples_by_pop = list()
        idx_sp = 1
        for (sp in sample_subpops){
          samples_by_pop[[idx_sp]] = take_sample_r(pop$pAB_post[sample_loci,sp, final_gen - g], pop$pAb_post[sample_loci,sp, final_gen - g], pop$paB_post[sample_loci,sp, final_gen - g], pop$pab_post[sample_loci,sp, final_gen - g], S)
          idx_sp = idx_sp + 1
        }
        samples_by_t[[idx_g]] = samples_by_pop
        idx_g = idx_g + 1
      }
      
      for (idx_g_1 in 1:length(sample_gens)){
        g1 = sample_gens[idx_g_1]
        
        s1s = samples_by_t[[idx_g_1]]
        
        r2s <- lapply(s1s, get_r2_drift_sample, S = S, c = c)
        r2_mat[idx_S, idx_g_1, idx_c, ] = unlist(r2s)
        
        for (idx_subpop_1 in 1:length(sample_subpops)){
          for (idx_subpop_2 in 1:length(sample_subpops)){
            sp1 = sample_subpops[idx_subpop_1]; sp2 = sample_subpops[idx_subpop_2] 
            if (sp1 == focal_pop | sp2 == focal_pop){
              rirj_mat[idx_S, idx_g_1, idx_c, idx_subpop_1, idx_subpop_2] = mean(get_riji(s1s[[idx_subpop_1]], s1s[[idx_subpop_2]]))
            }
          }
        }
        
        
        for (idx_g_2 in (min(idx_g_1+1, length(sample_gens))):length(sample_gens)){
          g2 = sample_gens[idx_g_2]
          s2s = samples_by_t[[idx_g_2]]
          
          for (idx_sp in 1:length(s2s)){
            Fc_mat[idx_S, idx_g_1, idx_g_2, idx_sp, idx_c ] = get_mean_fc_sample(s1s[[idx_sp]], s2s[[idx_sp]], S1 = S, S2 = S)
          }
          
          for (idx_subpop_1 in 1:length(sample_subpops)){
            for (idx_subpop_2 in 1:length(sample_subpops)){
              sp1 = sample_subpops[idx_subpop_1]; sp2 = sample_subpops[idx_subpop_2] 
              if (sp1 == focal_pop | sp2 == focal_pop){
                FcCor_mat[idx_S, idx_g_1, idx_g_2, idx_subpop_1, idx_subpop_2, idx_c] = get_FcCor( s1s[[idx_subpop_1]], s2s[[idx_subpop_1]], s1s[[idx_subpop_2]], s2s[[idx_subpop_2]] , S_size = S)
              }
            }
          }
        }
      }
      idx_S = idx_S+1
    }
    idx_c = idx_c+1
  }
  
  
  
  dimnames(Fc_mat)[[1]] = paste0('S', S_vec);
  dimnames(Fc_mat)[[2]] = paste0('g', sample_gens);
  dimnames(Fc_mat)[[3]] = paste0('g', sample_gens);
  dimnames(Fc_mat)[[4]] = paste0('sp', sample_subpops);
  dimnames(Fc_mat)[[5]] = paste0('c', recom_rates);
  
  dimnames(r2_mat)[[1]] = paste0('S', S_vec);
  dimnames(r2_mat)[[2]] = paste0('g', sample_gens);
  dimnames(r2_mat)[[3]] = paste0('c', recom_rates);
  dimnames(r2_mat)[[4]] = paste0('sp', sample_subpops);
  
  dimnames(FcCor_mat)[[1]] = paste0('S', S_vec);
  dimnames(FcCor_mat)[[2]] = paste0('g', sample_gens);
  dimnames(FcCor_mat)[[3]] = paste0('g', sample_gens);
  dimnames(FcCor_mat)[[4]] = paste0('sp', sample_subpops);
  dimnames(FcCor_mat)[[5]] = paste0('sp', sample_subpops);
  dimnames(FcCor_mat)[[6]] = paste0('c', recom_rates);
  
  dimnames(rirj_mat)[[1]] = paste0('S', S_vec);
  dimnames(rirj_mat)[[2]] = paste0('g', sample_gens);
  dimnames(rirj_mat)[[3]] = paste0('c', recom_rates);
  dimnames(rirj_mat)[[4]] = paste0('sp', sample_subpops);
  dimnames(rirj_mat)[[5]] = paste0('sp', sample_subpops);
  
  
  #print(11)
  return(list(pop = pop, Fc = Fc_mat, r2 = r2_mat, rirj = rirj_mat, FcCor = FcCor_mat))
}





### generating simulations in parallel###

# returns a model function f(N,m) --> summary stats
# takes model specification  parameters - n_demes, Mfn, S, sample_subpops, max_loci, min_loci
# and also takes paramteres specifying summary statistics to be calculated: sampling inteval, recom_rates
# the population used for the modelling is stored in the get_model_fn scope, which encloses model_fn

### only sample from focal population
get_model_fn = function(sampling_intervals, recom_rate, Mfn, S_vec, grid_size, sample_subpops, max_loci, min_loci, n_loci ){
  print(S_vec)
  pop = initialise_population(n_subpops = grid_size, loci = max_loci)
  
  model_fn = function(N,m, f = get_SS){
    o <- f(pop, N, m, Mfn, sampling_intervals, recom_rates, S_vec, grid_size, sample_subpops, max_loci, min_loci, n_loci)
    pop <<- o$pop
    return(o)
  }
}

# sample from multiple populations, also do correlations
### all correlations
get_model_fn_multi_pops = function(sample_gens, recom_rates, Mfn, S_vec, grid_size, sample_subpops, max_loci, min_loci, n_loci ){
  
  pop = initialise_population(n_subpops = grid_size, loci = max_loci)
  
  model_fn = function(N,m, f = get_SS_multi_pops2){
    o <- f(pop, N, m, Mfn, sample_gens, recom_rates, S_vec, grid_size, sample_subpops, max_loci, min_loci, n_loci)
    pop <<- o$pop
    return(o)
  }
}

### only correlations with focal population 
get_model_fn_multi_pops_focal = function(sample_gens, recom_rates, Mfn, S_vec, grid_size, sample_subpops, max_loci, min_loci, n_loci, focal_pop ){
  
  pop = initialise_population(n_subpops = grid_size, loci = max_loci)
  
  model_fn = function(N,m, f = get_SS_multi_pops3){
    o <- f(pop, N, m, Mfn, sample_gens, recom_rates, S_vec, grid_size, sample_subpops, max_loci, min_loci, n_loci, focal_pop)
    pop <<- o$pop
    return(o)
  }
}


### run in parallel
generate_ABC_sim_parallel = function(loci, sampling_intervals, recom_rates, Nm_list, model_fn_list, n_par, S_vec, Mfn, grid_size, sample_subpops){
  
  max_loci = (loci+5) * 1.1
  min_loci  = loci
  n_loci = loci
  
  cl = parallel::makeForkCluster(n_par, outfile = 'out.txt')
  doParallel::registerDoParallel(cl)
  
  
  SS_sim <- foreach (i = 0:(n_par-1), .noexport = c('propegate_population_inner')) %dopar%{
    Fc_out <- array(NA, dim = c(length(recom_rates), length(S_vec), length(sampling_intervals), 0))
    r2_out <- array(NA, dim = c(length(recom_rates), length(S_vec), length(sampling_intervals)+1, 0))
    
    
    Ns_loc = Nm_list[[i+1]][,1]
    ms_loc = Nm_list[[i+1]][,2]

    for (j in 1:length(Ns_loc)){

      o <- model_fn_list[[i+1]](Ns_loc[j], ms_loc[j], f = get_SS)

      Fc_out <- abind(Fc_out, o$Fc, along = 4)
      r2_out <- abind(r2_out, o$r2, along = 4)
    }

    
    list(Fc = Fc_out, r2 = r2_out)
  }
  stopCluster(cl)
  Nm_values = do.call(rbind, Nm_list)
  r2_out_glob = array(NA, dim = c(length(recom_rates), length(S_vec), length(sampling_intervals) + 1, 0))
  Fc_out_glob = array(NA, dim = c(length(recom_rates), length(S_vec), length(sampling_intervals)   , 0))

  for (threeD in SS_sim){

    Fc_out_glob = abind(Fc_out_glob, threeD$Fc, along = 4)
    r2_out_glob = abind(r2_out_glob, threeD$r2, along = 4)
  }
  return(list(Nm_values = Nm_values, r2_sim = r2_out_glob, Fc_sim = Fc_out_glob))
  
}

generate_ABC_sim_parallel_multi = function(loci, sample_gens, recom_rates, Nm_list, model_fn_list, n_par, S_vec, Mfn, grid_size, sample_subpops){
  

  max_loci = (loci+5) * 1.1
  min_loci  = loci
  n_loci = loci
  
  cl = parallel::makeForkCluster(n_par, outfile = 'out.txt')
  doParallel::registerDoParallel(cl)
  
  
  SS_sim <- foreach (i = 0:(n_par-1), .noexport = c('propegate_population_inner')) %dopar%{
    if (i == 0) s_inner = proc.time()
    
    Fc_out <- array(NA, dim = c( length(S_vec), length(sample_gens), length(sample_gens), length(sample_subpops), length(recom_rates), 0) )
    r2_out <- array(NA, dim = c( length(S_vec), length(sample_gens), length(recom_rates), length(sample_subpops), 0) )
    FcCor_out <- array(NA, dim = c( length(S_vec), length(sample_gens), length(sample_gens), length(sample_subpops), length(sample_subpops), length(recom_rates), 0) )
    rirj_out <- array(NA, dim = c( length(S_vec), length(sample_gens), length(recom_rates), length(sample_subpops), length(sample_subpops), 0) )
    first_burn = T
    Ns_loc = Nm_list[[i+1]][,1];
    ms_loc = Nm_list[[i+1]][,2];

    for (j in 1:length(Ns_loc)){
      
      gc()

      if (i == 0) print('')
      if (i == 0) print(paste(Ns_loc[j], ms_loc[j]))
      if (i == 0) s_loc = proc.time()
      
      o <- model_fn_list[[i+1]](Ns_loc[j], ms_loc[j], f = get_SS_multi_pops2)

      if (i == 0) print( paste0('that calc : ', (proc.time() - s_loc)[3] )); e_loc = proc.time()
      
      Fc_out <- abind(Fc_out, o$Fc, along = 6)
      r2_out <- abind(r2_out, o$r2, along = 5)
      FcCor_out <- abind(FcCor_out, o$FcCor, along = 7)
      rirj_out <- abind(rirj_out, o$rirj, along = 6)
      
      if (i == 0) print(paste0('other : ', (proc.time()-e_loc)[3]))
      if (i == 0) print(paste0('total : ', (proc.time()-s_inner)[3] ))
      
      
    }
    
    
    list(Fc = Fc_out, r2 = r2_out, FcCor = FcCor_out, rirj = rirj_out)
  }
  stopCluster(cl)
  Nm_values = do.call(rbind, Nm_list)
  
  Fc_out_pool = array(NA, dim = c( length(S_vec), length(sample_gens), length(sample_gens), length(sample_subpops), length(recom_rates), 0) )
  r2_out_pool = array(NA, dim = c( length(S_vec), length(sample_gens), length(recom_rates), length(sample_subpops), 0) )
  FcCor_out_pool = array(NA, dim = c( length(S_vec), length(sample_gens), length(sample_gens), length(sample_subpops), length(sample_subpops), length(recom_rates), 0) )
  rirj_out_pool = array(NA, dim = c( length(S_vec), length(sample_gens), length(recom_rates), length(sample_subpops), length(sample_subpops), 0) )
  
  for (Nm_subset_SS in SS_sim){
    
    Fc_out_pool = abind(Fc_out_pool, Nm_subset_SS$Fc, along = 6)
    r2_out_pool = abind(r2_out_pool, Nm_subset_SS$r2, along = 5)
    FcCor_out_pool = abind(FcCor_out_pool, Nm_subset_SS$FcCor, along = 7)
    rirj_out_pool = abind(rirj_out_pool, Nm_subset_SS$rirj, along = 6)
    }
  return(list(Nm_values = Nm_values, Fc_sim = Fc_out_pool, r2_sim = r2_out_pool, FcCor_sim = FcCor_out_pool, rirj_sim = rirj_out_pool))
  
}

get_square_around_centre = function(grid_size, square_size){
  if (grid_size%%2 == 1){
  centre = (grid_size / 2 + .5)
  centre_row  = (centre - (square_size-1)/2):(centre + (square_size-1)/2)
  edge_length = sqrt(grid_size)
  final = c()
  (square_size - 1)/2
  for (i in -((square_size - 1)/2):((square_size - 1)/2)){
  final = c(final, centre_row + i*edge_length)
  }
  }
  
  else{
    LL = grid_size / 2 - sqrt(grid_size)/2
    final = c()
    for (i in 1:square_size-1){
      final = c(final, (LL+(1:square_size)-1)+ i*sqrt(grid_size))
    }
  }
  return(final)
}



### processing summary statistics ###

# choose which summary statistics to keep
format_SS <- function(Fc_arr, r2_arr, cs = 'all', Ss = 'all', gs = 'all'){
  if (cs == 'all') c_labs = dimnames(Fc_arr)[[1]]
  else {
    c_labs <- paste0('c', cs)
  }
  
  if (Ss == 'all') S_labs = dimnames(Fc_arr)[[2]]
  else{
    S_labs <- paste0('S', Ss)
  }
  
  if (gs == 'all') {g_labs_Fc = dimnames(Fc_arr)[[3]]; g_labs_r2 = dimnames(r2_arr)[[3]]}
  else{
    g_labs_Fc <- paste0('g', gs); g_labs_r2 <- paste0('g', c(0,gs))
  }
  
  Fc_sub = Fc_arr[c_labs, S_labs, g_labs_Fc, , drop = F]
  Fc_out = apply(Fc_sub, c(2,3,4), mean)
  
  r2_sub = r2_arr[c_labs, S_labs, g_labs_r2, ,drop = F]
  r2_out = aperm(apply(r2_sub, c(1,2,4), mean), c(2,1,3))
  return(list(Fc = Fc_out, r2 = r2_out))
}

# reshape arrays
make_SS_vecs <- function(Fc3D, r23D){
  r2 = matrix(r23D[,,], nr = dim(r23D)[3], byrow = T)
  Fc = matrix(Fc3D[,,], nr = dim(Fc3D)[3], byrow = T)
  
  return(list(Fc = Fc, r2 = r2))
}

# convert to Ne estimates if necessary
do_construction_LD <- function(SS_LD_n_col, to_Ne = F, cs = NULL){
  out = array(NA, dim = dim(SS_LD_n_col))
  
  if (to_Ne){
    for (col in 1:length(cs)){
      out[,col] = ( 1 - SS_LD_n_col[,col]) / ( 2*cs[col]*(2-cs[col])*SS_LD_n_col[,col] )
    }
    return(out)
  }
  
  out[,1] = SS_LD_n_col[,1]
  for (col in 2:dim(out)[2]){
    out[,col] = SS_LD_n_col[,col]
  }
  return(out)
}

# convert to Ne estimates if necessary
do_construction_AF <-function(SS_AF_n_col, to_Ne = F, ts = NULL, S = NULL){
  out = array(NA, dim = dim(SS_AF_n_col))
  
  if (to_Ne){
    for (col in 1:length(ts)){
      out[,col] = .5/ ((SS_AF_n_col[,col] - 1/(2*S) - 1/(2*S) ) / ts[col])
    }
    return(out)
  }
  
  
  
  
  if (dim(SS_AF_n_col)[2] == 1) return(SS_AF_n_col)
  out = array(NA, dim = dim(SS_AF_n_col))
  out[,1] = SS_AF_n_col[,1]
  for (col in 2:dim(out)[2]){
    out[,col] = SS_AF_n_col[,col]
  }
  return(out)
}

diff_range = function(v){
  return(diff(range(v)))
}
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


### old way ###

# Varying migration rate
get_estimates_df <- function(N_v, ms, M_fn, p0, loci, burn_in, t_intervals, rest_period, sampling_repeats, sample_subpops='all', S=Inf, cs = c(.5) ){
  
  n_subpops <- length(N_v)
  sampling_period <- reset_period + max(t_intervals)
  if (sample_subpops == "all"){
    sample_subpops <- 1:length(N_v)
  }
  sample_subpops <- c(sample_subpops)
  
  
  v=rep(NA,length(cs) * length(ms) * length(sample_subpops) * (length(t_intervals) + 2*(length(t_intervals) - 1) ) * sampling_repeats )
  estimates = data.frame(m = v, sp = v, method = v, i = v, est = v)
  row = 1
  
  for (m in ms){
    cat(paste("m = ", as.character(m), '\n'))
    M = M_fn(m = m, n_subpops = n_subpops)
    
    t = burn_in + sampling_repeats * (max(t_intervals) + reset_period)
    
    for (c in cs){
      cat(paste("    c =", c, '\n'))
      cat("        Simulating\n")
      o = sim.drift5_mig(N_v = N_v, t = t,  M=M, loci = loci,  p0 = p0, c = c, burn_in = burn_in)
      
      cat("        Estimating\n")
      for (sp in sample_subpops){
        print(paste("            Sampling from", sp))
        for (i in 0:(sampling_repeats-1)){
          sample_gens = burn_in + t_intervals + sampling_period*i
          
          # LD ests
          v = c()
          for (gen in sample_gens){
            # take generation, subpopulation sample (potentially infinite)
            method_name = paste0("LD", as.character(c))
            s <- take_sample_r(o$pAB_post[,gen,sp], o$pAb_post[,gen,sp], o$paB_post[,gen,sp], o$pab_post[,gen,sp], S)
            estimates[row,] = c(m,sp,method_name, i, get_LD_est(s$pAB, s$pAb, s$paB, s$pab, S = S, c = c))
            row = row + 1
          }
          
          # AF ests
          s1 <- take_sample_r(o$pAB_post[,sample_gens[1],sp], o$pAb_post[,sample_gens[1],sp], o$paB_post[,sample_gens[1],sp], o$pab_post[,sample_gens[1],sp], S)
          for (gen in sample_gens[-1]){
            method_name = paste0("AF", gen - sample_gens[1])
            s2 <- take_sample_r(o$pAB_post[,gen,sp], o$pAb_post[,gen,sp], o$paB_post[,gen,sp], o$pab_post[,gen,sp], S)
            
            estimates[row,] =  c(m,sp, method_name, i, .5/get_AF_est(s1$pAB+s1$pAb, s2$pAB+s2$pAb, S, S, gen-sample_gens[1]))
            row = row + 1
            estimates[row,] =  c(m,sp, method_name, i, .5/get_AF_est(s1$pAB+s1$paB, s2$pAB+s2$paB, S, S, gen-sample_gens[1]))
            row = row + 1
          }
        }
      }
    }
  }
  return(estimates)
}

get_SS_df <- function(N_v, ms, M_fn, p0, loci, burn_in, t_intervals, rest_period, sampling_repeats, sample_subpops='all', S=Inf, cs = c(.5) ){
  
  n_subpops <- length(N_v)
  sampling_period <- reset_period + max(t_intervals)
  if (sample_subpops == "all"){
    sample_subpops <- 1:length(N_v)
  }
  sample_subpops <- c(sample_subpops)
  
  
  v=rep(NA,length(cs) * length(ms) * length(sample_subpops) * (length(t_intervals) + 2*(length(t_intervals) - 1) ) * sampling_repeats )
  estimates = data.frame(m = v, sp = v, method = v, i = v, est = v)
  row = 1
  
  for (m in ms){
    cat(paste("m = ", as.character(m), '\n'))
    M = M_fn(m = m, n_subpops = n_subpops)
    
    t = burn_in + sampling_repeats * (max(t_intervals) + reset_period)
    
    for (c in cs){
      cat(paste("    c =", c, '\n'))
      cat("        Simulating\n")
      o = sim.drift5_mig(N_v = N_v, t = t,  M=M, loci = loci,  p0 = p0, c = c, burn_in = burn_in)
      
      cat("        Estimating\n")
      for (sp in sample_subpops){
        print(paste("            Sampling from", sp))
        for (i in 0:(sampling_repeats-1)){
          sample_gens = burn_in + t_intervals + sampling_period*i
          
          # LD ests
          v = c()
          for (gen in sample_gens){
            # take generation, subpopulation sample (potentially infinite)
            method_name = paste0("r", as.character(c))
            s <- take_sample_r(o$pAB_post[,gen,sp], o$pAb_post[,gen,sp], o$paB_post[,gen,sp], o$pab_post[,gen,sp], S)
            estimates[row,] = c(m,sp,method_name, i, get_r2_drift(s$pAB, s$pAb, s$paB, s$pab, S = S, c = c))
            row = row + 1
          }
          
          # AF ests
          s1 <- take_sample_r(o$pAB_post[,sample_gens[1],sp], o$pAb_post[,sample_gens[1],sp], o$paB_post[,sample_gens[1],sp], o$pab_post[,sample_gens[1],sp], S)
          for (gen in sample_gens[-1]){
            method_name = paste0("AF", gen - sample_gens[1])
            s2 <- take_sample_r(o$pAB_post[,gen,sp], o$pAb_post[,gen,sp], o$paB_post[,gen,sp], o$pab_post[,gen,sp], S)
            
            estimates[row,] =  c(m,sp, method_name, i, get_mean_fc(s1$pAB+s1$pAb, s2$pAB+s2$pAb, S, S, gen-sample_gens[1]))
            row = row + 1
            estimates[row,] =  c(m,sp, method_name, i, get_mean_fc(s1$pAB+s1$paB, s2$pAB+s2$paB, S, S, gen-sample_gens[1]))
            row = row + 1
          }
        }
      }
    }
  }
  return(estimates)
}

get_means_df <- function(estimates_df){
  ms <- unique(estimates_df$m)
  methods <-unique(estimates_df$method)
  
  v = rep(NA, length(ms)*length(methods))
  mean_estimates_df = data.frame(m = v, method = v, est = v)
  row = 1
  print(methods)
  print(ms)
  for (m_i in ms){
    for (method_i in methods){
      d = as.numeric(subset(estimates_df, m == m_i & method == method_i)$est)
      mean_estimates_df[row,] = c(m_i, method_i, mean( d , na.rm=T) )
      row=row+1
      print(length(d)) # how many values are used to calculate mean
    }
  }
  
  mean_estimates_df$m = as.numeric(mean_estimates_df$m)
  mean_estimates_df$est = as.numeric(mean_estimates_df$est)
  
  return(mean_estimates_df)
}

trim_df <- function(means_df, methods = "all", max_m = 1){
  if (methods == "all") {methods_trimmed <- means_df}
  else {methods_trimmed <- means_df[means_df$method %in% methods,]}
  final <- methods_trimmed[methods_trimmed$m <= max_m,]
  return (final)
}


###
#
# sampling and estimating functions
#
###

### AF ###

# take single locus sample for each locus in frequencies vector
take_sample <- function(frequencies_vector, S){
  return(rbinom(length(frequencies_vector), 2*S, frequencies_vector) / (2*S))
}

# get Fc from frequencies of alleles A & a at times 1 & 2 (can be vectors)
freq_to_fc <- function( A_t1, A_t2, a_t1, a_t2){
  return( ( ( (A_t1 - A_t2)**2 / ((A_t1 + A_t2)/2 - A_t1 * A_t2) ) + ( (a_t1-a_t2)**2 / ((a_t1+a_t2)/2 - a_t1*a_t2) ) ) / 2 )
}
freq_to_fc = freq_to_fc_C

freq_to_fc_sqrt <- function( A_t1, A_t2, a_t1, a_t2){
  return( ( ( (A_t1 - A_t2) / sqrt((A_t1 + A_t2)/2 - A_t1 * A_t2) ) ) )
}
freq_to_fc_sqrt = freq_to_fc_sqrt_C

# get vector of Fc values from frequency of A at time 1 & 2
get_fc_vec <- function(A_t1, A_t2){
  B_t1 <- 1 - A_t1
  B_t2 <- 1 - A_t2
  
  FC <- freq_to_fc(A_t1, A_t2, B_t1, B_t2)
  
  return(FC)
}

get_fc_vec_sqrt <- function(A_t1, A_t2){
  B_t1 <- 1 - A_t1
  B_t2 <- 1 - A_t2
  
  FC <- freq_to_fc_sqrt(A_t1, A_t2, B_t1, B_t2)
  
  return(FC)
}

# get AF estiamte from frequencies A(t1), A(t2), sample sizes, and t2-t1
get_AF_est <- function(A_t1, A_t2, S1, S2, t){
  o <- filter_AF(A_t1, A_t2)
  A_t1 = o$t1
  A_t2 = o$t2
  mean_fc <- mean(get_fc_vec(A_t1, A_t2), na.rm = T)
  INV_est <- (mean_fc - 1/(2*S1) - 1/(2*S2) ) / t
  return(INV_est)
}

get_mean_fc <- function(A_t1, A_t2, S1, S2, t){
  o <- filter_AF(A_t1, A_t2)
  A_t1 = o$t1
  A_t2 = o$t2
  mean_fc <- mean(get_fc_vec(A_t1, A_t2), na.rm = T)
  return(mean_fc)
}

get_mean_fc_sample <- function(samp_t1, samp_t2, S1, S2, t){
  #print(t)
  A_t1 = samp_t1$pAB + samp_t1$pAb
  A_t2 = samp_t2$pAB + samp_t2$pAb
  
  B_t1 = samp_t1$pAB + samp_t2$paB
  B_t2 = samp_t2$pAB + samp_t2$paB
  
  return((get_mean_fc(A_t1, A_t2, S1, S2, t) + get_mean_fc(B_t1, B_t2, S1, S2, t)) / 2 )
}

get_fcv_sample <- function(samp_t1, samp_t2, S1, S2, t){
  A_t1 = samp_t1$pAB + samp_t1$pAb
  A_t2 = samp_t2$pAB + samp_t2$pAb
  
  B_t1 = samp_t1$pAB + samp_t2$paB
  B_t2 = samp_t2$pAB + samp_t2$paB
  
  return(list(get_fc_vec(A_t1, A_t2, S1, S2, t), get_fc_vec(B_t1, B_t2, S1, S2, t)) )
}

get_fcv_sample_freqs <- function(samp_t1, samp_t2, S1, S2, t){
  A_t1 = samp_t1$pA
  A_t2 = samp_t2$pA
  
  B_t1 = samp_t1$pB
  B_t2 = samp_t2$pB
  
  return(list(get_fc_vec(A_t1, A_t2), get_fc_vec(B_t1, B_t2)) )
}

get_fcv_sample_freqs_sqrt <- function(samp_t1, samp_t2, S1, S2, t){
  A_t1 = samp_t1$pA
  A_t2 = samp_t2$pA
  
  B_t1 = samp_t1$pB
  B_t2 = samp_t2$pB
  
  return(list(get_fc_vec_sqrt(A_t1, A_t2), get_fc_vec_sqrt(B_t1, B_t2)) )
}

# filter out alleles that start at frequency >0.95 or <0.05
filter_AF <- function(A_t1, A_t2){
  filt =  (A_t1>0.05 & A_t1<0.95) 
  return(list(t1 = A_t1[filt], t2 = A_t2[filt]))
}

filter_AF_pair = function(A_t1_s1, A_t2_s1, A_t1_s2, A_t2_s2){
  filt1 =  (A_t1_s1>0.05 & A_t1_s1<0.95) 
  filt2 =  (A_t1_s2>0.05 & A_t1_s2<0.95) 
  filt = (filt1&filt2)
  
  return(list( A_t1_s1[filt], A_t2_s1[filt], A_t1_s2[filt], A_t2_s2[filt]) )
  
}
filter_AF_pair = filter_AF_pair_C

get_FcCor <- function(S1_t1, S1_t2, S2_t1, S2_t2, S_size){
  As = filter_AF_pair(S1_t1$pAB+S1_t1$pAb, S1_t2$pAB+S1_t2$pAb, S2_t1$pAB+S2_t1$pAb, S2_t2$pAB+S2_t2$pAb)
  Bs = filter_AF_pair(S1_t1$pAB+S1_t1$paB, S1_t2$pAB+S1_t2$paB, S2_t1$pAB+S2_t1$paB, S2_t2$pAB+S2_t2$paB)
  
  S1_t1 = list(pA = As[[1]], pB = Bs[[1]])
  S1_t2 = list(pA = As[[2]], pB = Bs[[2]])
  S2_t1 = list(pA = As[[3]], pB = Bs[[3]])
  S2_t2 = list(pA = As[[4]], pB = Bs[[4]])
  
  fcv1 = get_fcv_sample_freqs_sqrt(S1_t1, S1_t2, S1 = S_size, S2 = S_size, t = NA)
  fcv2 = get_fcv_sample_freqs_sqrt(S2_t1, S2_t2, S1 = S_size, S2 = S_size, t = NA)


  
  cA = cor(fcv1[[1]], fcv2[[1]]); cB = cor(fcv1[[2]], fcv2[[2]]);
  
  if (is.na(cA)) return(NA)
  
  return(mean(cA, cB))
}


### LD ###

# take multinomial sample of haplotype frequencies
take_sample_r <- function(pAB, pAb, paB, pab, S){
  if (is.finite(S)){
    temp = gsl_mmm(2*S, cbind(pAB, pAb, paB, pab)) / (2*S)
    return(list(pAB = temp[,1], pAb = temp[,2], paB = temp[,3], pab = temp[,4], pA = temp[,1]+temp[,2], pB = temp[,1]+temp[,3]))
  }
  else{
    return(list(pAB=pAB, pAb=pAb, paB=paB, pab=pab))
  }

}

# calculate r2 vector from vectors of haplotype frequencies
freqs_to_r2 <- function(pAB, pAb, paB, pab){
  pA <- pAB + pAb
  pB <- pAB + paB
  r <- (pAB*pab-pAb*paB)**2 / (pA*(1-pA)* pB*(1-pB))
  return(r)
}

# adjust r2 vector to correct for sampling
r2_drift <- function(r_vec, S){
  mean_r <- mean(r_vec, na.rm=T)
  r_drift <- (mean_r -  1/(2*S) ) / (1 - 1/(2*S))
  return(r_drift)
}

# get LD estiamtes from hamplotype frequencies and sample size
get_LD_est <- function(pAB, pAb, paB, pab, S, c = .5){
  o <- filter_LD(pAB, pAb, paB, pab)
  pAB <- o$pAB; pAb <- o$pAb; paB <- o$paB; pab <- o$pab;
  r_vec <- freqs_to_r2(pAB, pAb, paB, pab)
  r_drift <- r2_drift(r_vec, S)
  Ne_est <- ( (1-c)**2 + c**2 ) / ( 2*c*(2-c)*r_drift )
  Ne_est <- ( 1 - r_drift) / ( 2*c*(2-c)*r_drift )
  
  return(Ne_est)
}

get_r2_drift<- function(pAB, pAb, paB, pab, S, c = .5){
  
  o <- filter_LD(pAB, pAb, paB, pab)
  pAB <- o$pAB; pAb <- o$pAb; paB <- o$paB; pab <- o$pab;
  r_vec <- freqs_to_r2(pAB, pAb, paB, pab)
  r_drift <- r2_drift(r_vec, S)
  
  return(r_drift)
}

get_r2_drift_sample <- function(samp, S, c = .5){
  pAB = samp$pAB
  pAb = samp$pAb
  paB = samp$paB
  pab = samp$pab
  
  return(get_r2_drift(pAB, pAb, paB, pab, S, C))
}

filter_LD <- function(pAB, pAb, paB, pab){
  pA <- pAB + pAb
  pB <- pAB + paB
  filt <- (pA>0.05 & pA<0.95 & pB>0.05 & pB<0.95)

  return(list(pAB = pAB[filt], pAb = pAb[filt], paB = paB[filt], pab = pab[filt]))
}

filter_LD_pair<- function(S1,S2){
  pA1 <- S1$pAB + S1$pAb
  pB1 <- S1$pAB + S1$paB
  filt1 <- (pA1>0.05 & pA1<0.95 & pB1>0.05 & pB1<0.95)
  
  
  pA2 <- S2$pAB + S2$pAb
  pB2 <- S2$pAB + S2$paB
  filt2 <- (pA2>0.05 & pA2<0.95 & pB2>0.05 & pB2<0.95)
  
  filt = (filt1 & filt2)
  
  return(list(S1 = list(pAB = S1$pAB[filt], pAb = S1$pAb[filt], paB = S1$paB[filt], pab = S1$pab[filt], pA = S1$pAB[filt]+S1$pAb[filt], pB = S1$pAB[filt]+S1$paB[filt]),
              S2 = list(pAB = S2$pAB[filt], pAb = S2$pAb[filt], paB = S2$paB[filt], pab = S2$pab[filt], pA = S2$pAB[filt]+S2$pAb[filt], pB = S2$pAB[filt]+S2$paB[filt]))
         )
}

get_riji = function(S1, S2){
  o = filter_LD_pair(S1,S2)
  S1 = o$S1; S2 = o$S2
  
  r1 = (S1$pAB*S1$pab-S1$pAb*S1$paB) / sqrt(S1$pA*(1-S1$pA)* S1$pB*(1-S1$pB))
  r2 = (S2$pAB*S2$pab-S2$pAb*S2$paB) / sqrt(S2$pA*(1-S2$pA)* S2$pB*(1-S2$pB))
  return(r1*r2/(sd(r1)*sd(r2)))
}


### Fst ###

# calculate Fst for each locus between two populations, given frequency vectors and population sizes
get_Fst <- function(freqs_1, freqs_2, N1, N2){
  freqs_mean <- (N1*freqs_1 + N2*freqs_2) / (N1+N2)
  H_t <- 2*freqs_mean*(1-freqs_mean)
  
  
  H_1 <- 2*freqs_1*(1-freqs_1)
  H_2 <- 2*freqs_2*(1-freqs_2)
  H_s <- (H_1 + H_2) / 2
  print(H_t)
  Fst <- (H_t - H_s) / H_t
  return (Fst)
}

# calculate vector of Fst values over generations
get_Fst_vec <- function(freqs_arr_1, freqs_arr_2, N1, N2, g_max){
  
  Fst_v <- c()
  for (g in 1:g_max){
    Fst_v <- c(Fst_v, mean(get_Fst(freqs_arr_1[,g], freqs_arr_2[,g], N1, N2) ))
  }
 return(Fst_v) 
}

# get Fst vector for many subpopulations
get_Fst_many <- function(freqs, N_v){
  freqs_means <- rowSums(t(t(freqs)*N_v)) / sum(N_v)
  
  num <- rowVars(freqs)
  den <- 2 * freqs_means * (1-freqs_means)

  return(num/den)
}

get_Fst_vec <- function(freqs_arr, N_v, g_max){
  v <- rep(NA, g_max)
  for (g in 1:g_max){
    v[g] <- mean(get_Fst_many(freqs_arr[,,g], N_v), na.rm=T)
  }
  return (v)
}

burn_in_checker <- function(ms, M_fn, cs, burn_in, N_v, loci, p0){
  M <- M_fn( m = min(ms), n_subpops = length(N_v))
  c <- min(cs)
  o <- sim.drift5_mig(N_v = N_v, t = burn_in, M=M, loci = loci, p0=p0, c = c)
  
  par(mfrow = c(2,2), mai = c(1,1,0.5,0.5))
  Fst <- get_Fst_vec(o$pA_pre, N_v, burn_in)
  plot(1:(burn_in), ma(Fst, n =1), type = "l", xlab = "generation")
  
  r2 <- colMeans(o$r_post[,,1], na.rm = T)
  plot(1:(burn_in+1), r2, type = 'l', xlab = "generation")
 
  
  ave_over_sp <- apply(o$pA_post, c(1,2), mean)
  fixation <- apply(ave_over_sp, 2, f)
  plot(1:(burn_in+1),fixation/loci, type='l')
  
  hist(ave_over_sp[,burn_in], breaks =seq(-0.05,1.05, 0.1))
  return (o)
}

f <- function(v){
  return(sum(abs(v-1) < 1e-10 | v < 1e-10))
}

get_methods <- function(t_intervals, cs){
  methods <- c()
  for (t in t_intervals[-1]){
    methods <- c(methods, paste0("AF", t-min(t_intervals)))
  }
  for (c in cs){
    methods <- c(methods, paste0("LD", c))
  }
  return(methods)
}


# colouring in
ma <- function(x, n = 5){ c(stats::filter(x, rep(1 / n, n), sides = 1)[n:length(x)], rep(NA, n-1)) }

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


###
#
# Migration matrices
#
###

# supplied m is the proortion of a subpopulation's gamete pool which comes from other subpopulations each generation.
# so panmixia is at m = (n_subpops-1)/n_subpops, NOT m = 1


# check that a matrix is a valid migration matrix. A valid M has all columns summing to 1:
# [ mAA, mAB]
# [ mBA, mBB]
# where mBA is the proportion of A's gamete pool which came from subpopulation B (ie migration from B to A)
check_M <- function(M){
  return(all(colSums(M) - 1 < 1e-14))
}

# symmetrical migration amongst n_subpops
construct_M_max_1 <- function(m, size_param){
  n_subpops = size_param
  M <- array(rep(m/n_subpops, n_subpops**2), dim = c(n_subpops, n_subpops))
  diag(M) <- 1- (m/n_subpops)*(n_subpops-1)
  return ( M )
}

# symmetrical migration amongst n_subpops, where m is parameterised correctly! ie as in Waples 2011
# size param is the number of subpopulations
construct_M_waples <- function(m, n_subpops){
  M <- array(rep(m/(n_subpops-1), n_subpops**2), dim = c(n_subpops, n_subpops))
  diag(M) <- 1 - m
  return ( M )
}

# size param is the edge_length
construct_M_4_step <- function(m, n_subpops){
  
  stopifnot(sqrt(n_subpops)%%1 == 0)
  m_each = m/4
  edge_length = sqrt(n_subpops)
  M <- array(data = NA, dim = c(n_subpops, n_subpops))
  for (reciever in 1:n_subpops){
    to_subpops <- c()
    if (reciever %% edge_length != 1){
      to_subpops = c(to_subpops, reciever-1)
    }
    if (reciever %% edge_length != 0){
      to_subpops = c(to_subpops, reciever+1)
    }
    if (reciever + edge_length < n_subpops){
      to_subpops = c(to_subpops, reciever + edge_length)
    }
    if (reciever - edge_length > 0){
      to_subpops = c(to_subpops, reciever - edge_length)
    }
    col_vals = rep(0,n_subpops)
    col_vals[to_subpops] = m_each
    col_vals[reciever] = 1 - sum(col_vals)
    M[,reciever] = col_vals
  }
  return(M)
}

construct_M_hex_step <- function(size_param, m){
  
}

# unidirectional migration
# M is of form:
#  [1-m , 0 ]
#  [ m  , 1 ]
construct_M_inf <- function(m, n_subpops = "irrelevant"){
  return (array(data = c(1-m, m, 0, 1), dim = c(2,2)))
}


###
#
# ABC 
#
###

m_prior_single_draw = function(){
  beta_out = 1
  while(beta_out > 0.8){
  beta_out = rbeta(1,.7, 1)
  beta_out[beta_out<0.01] = 0
  }
  return(beta_out)
}

m_prior_draw = function(n){
  v = rep(NA, n)
  for (i in 1:n){
    v[i] = m_prior_single_draw()
  }
  return(v)
}

N_prior_draw = function(n=1){
  return(runif(n, 200,25000))
}

Nm_draw <- function(size){
  Nm_unordered <- array(c(N_prior_draw(size), m_prior_draw(size)), dim = c(size, 2))
  
  dist_matrix = as.matrix(dist(scale(Nm_unordered), upper = T))
  diag(dist_matrix) = NA
  
  Nm_ordered <- array(NA, dim = c(size, 2))
  curr = which.max(Nm_unordered[,2])
  for (sample in 1:(size-1)){
    
    Nm_ordered[sample, ] = Nm_unordered[curr,]
    order = c(order, curr)
    r = dist_matrix[curr,]
    dist_matrix = dist_matrix[, colnames(dist_matrix) != as.character(curr)]
    curr = as.numeric(labels(which.min(r)))
  }
  Nm_ordered[size,] = Nm_unordered[curr,]
  order = c(order, curr)
  return(list(Nm_ordered, Nm_unordered))
}

sort_Nm_draw = function(Nm_unordered){
  size = dim(Nm_unordered)[1]
  dist_matrix = as.matrix(dist(scale(Nm_unordered), upper = T))
  diag(dist_matrix) = NA
  
  Nm_ordered <- array(NA, dim = c(size, 2))
  curr = which.max(Nm_unordered[,2])
  for (sample in 1:(size-1)){
    
    Nm_ordered[sample, ] = Nm_unordered[curr,]
    order = c(order, curr)
    r = dist_matrix[curr,]
    dist_matrix = dist_matrix[, colnames(dist_matrix) != as.character(curr)]
    curr = as.numeric(labels(which.min(r)))
  }
  Nm_ordered[size,] = Nm_unordered[curr,]
  order = c(order, curr)
  return(Nm_ordered)
}

Nm_draw_large <- function(size){
  Nm_u <- array(c(N_prior_draw(size), m_prior_draw(size)), dim = c(size, 2))
  
  split_and_sort = function(Nm_u){
    if (dim(Nm_u)[1] < 1000){
      return( sort_Nm_draw(Nm_u))
    }
    else {
      median_N = median(Nm_u[,1])
      median_m = median(Nm_u[,2])
      split1 = subset(Nm_u, Nm_u[,1] < median_N & Nm_u[,2] < median_m)
      split2 = subset(Nm_u, Nm_u[,1] > median_N & Nm_u[,2] < median_m)
      split3 = subset(Nm_u, Nm_u[,1] > median_N & Nm_u[,2] > median_m)
      split4 = subset(Nm_u, Nm_u[,1] < median_N & Nm_u[,2] > median_m)
      return (rbind(split_and_sort(split1), split_and_sort(split2), split_and_sort(split3), split_and_sort(split4)))
    }
  }
  return (split_and_sort(Nm_u))
}

do_abc = function(SS_sim, SS_real, Nm_values, N = NA, m = NA, draw = F, method = 'neuralnet', tol = .1, pca = F, transformation = c('none', 'none', 'none')){
  
  ### construct some posteriors:
  
  N_draws = dim(SS_sim)[1]
  
  # scale data
  scaled <-scale(rbind(SS_sim, SS_real) )
  stat.sim.scaled <- scaled[-(N_draws+1),]
  stat.obs.scaled <- scaled[N_draws+1,]
  
  if (draw){
    
    par(mfrow = c(3,3), mai = c(0.4,0.3,0.1,0.1))
    options(repr.plot.width=15, repr.plot.height=15)
    for (i in 1:dim(stat.sim.scaled)[2]){
      plot(Nm_values[,1], stat.sim.scaled[,i], main = paste('N vs stat', i))
    }
    
    par(mfrow = c(3,3), mai = c(0.4,0.3,0.1,0.1))
    for (i in 1:dim(stat.sim.scaled)[2]){
      plot(Nm_values[,2], stat.sim.scaled[,i], main = paste('m vs stat', i))
    }
    
    
    par(mfrow = c(3,3), mai = c(0.4,0.3,0.1,0.1))
    options(repr.plot.width=15, repr.plot.height=15)
    for (i in 1:dim(stat.sim.scaled)[2]){
      hist(stat.sim.scaled[,i])
      abline(v = stat.obs.scaled[i], col = "red", lty = 2, main = paste('stat', i))
    }
    
  }
  
  # construct PCs
  pc <- prcomp(stat.sim.scaled, scale = F, center= F)
  pc
  
  if (pca){
    pc.use <- pca
    pc.sim <- (pc$x[,pca])
    
    
    
    stat.obs.pc <- c()
    for (p in pc.use){
      stat.obs.pc <- c(stat.obs.pc, sum(stat.obs.scaled * pc$rotation[,p]))
    }
    
    if (draw){
      
      par(mfrow = c(3,3), mai = c(0.4,0.3,0.1,0.1))
      for (i in 1:dim(stat.sim.scaled)[2]){
        plot(Nm_values[,1], pc$x[,i], main = paste0("N vs PC",as.character(i)), pch = 16)
      }
      
      par(mfrow = c(3,3), mai = c(0.4,0.3,0.1,0.1))
      for (i in 1:dim(stat.sim.scaled)[2]){
        plot(Nm_values[,2], pc$x[,i], main = paste0("m vs PC",as.character(i)), pch = 16)
      }
      
      par(mfrow = c(3,3), mai = c(0.4,0.3,0.1,0.1))
      options(repr.plot.width=15, repr.plot.height=15)
      for (i in 1:dim(stat.sim.scaled)[2]){
        hist(pc$x[,i])
        abline(v = stat.obs.pc[i], col = "red", lty = 2, main = paste('PC', i))
      }
      
      
    }
    
    
    
  }
  if(pca){
    rej <- abc(target=stat.obs.pc, param=Nm_values, sumstat=pc.sim, tol=tol, method = method, transf = transformation) 
  }
  else{
    rej <- abc(target=stat.obs.scaled, param=Nm_values, sumstat=stat.sim.scaled, tol=tol, method = method,  transf = transformation) 
  }
  if (draw){
    par(mar=c(5.5,5.5,1,1))
    layout(matrix(c(1,1,1,2,3,4,5,5,5,6,7,8),ncol=3, byrow = T),heights=c(1,3,1,3))
    plot.new()
    text(0.5,0.5,"Adjusted values",cex=2,font=2)
    hist(rej$adj.values[,1], main = 'Posterior N', xlab = 'N', xlim = c(1000,21000), breaks=seq(-10000,250000,1000)); abline(v=N, lty = 2, col = 'red', cex =2);abline(v=summary(rej)[4,1], lty = 2, col = 'green', cex =2)
    hist(rej$adj.values[,2], main = 'Posterior m', xlab = 'm', xlim = c(0,1), breaks = seq(-1,2,.05)); abline(v=m, lty = 2, col = 'red', cex =2);abline(v=summary(rej)[4,2], lty = 2, col = 'green', cex =2)
    hist(rej$adj.values[,3], main = 'Posterior Nm', xlab = 'Nm', xlim = c(0,20000), breaks = seq(-100000,300000,500)); abline(v=N*m, lty = 2, col = 'red', cex =2);abline(v=(summary(rej)[4,3]), lty = 2, col = 'green', cex =2)
    
    plot.new()
    text(0.5,0.5,"Unadjusted values",cex=2,font=2)
    hist(rej$unadj.values[,1], main = 'Posterior N', xlab = 'N', xlim = c(1000,21000), breaks=seq(-10000,250000,1000)); abline(v=N, lty = 2, col = 'red', cex =2);abline(v=summary(rej)[4,1], lty = 2, col = 'green', cex =2)
    hist(rej$unadj.values[,2], main = 'Posterior m', xlab = 'm', xlim = c(0,1), breaks = seq(-1,2,.05)); abline(v=m, lty = 2, col = 'red', cex =2);abline(v=summary(rej)[4,2], lty = 2, col = 'green', cex =2)
    hist(rej$unadj.values[,3], main = 'Posterior Nm', xlab = 'Nm', xlim = c(0,20000), breaks = seq(-100000,300000,500)); abline(v=N*m, lty = 2, col = 'red', cex =2);abline(v=(summary(rej)[4,3]), lty = 2, col = 'green', cex =2)
    
  }
  
  return(rej)
}

get_dif = function(n){
  return(mean(abs(runif(n,2000,20000)-runif(n,2000,20000))))
}


# functions for processing and analysing multiple subpopualtion sim files. 

  
make_SS = function(sim_out, real_out_list, S_select, gens_select, subpops_select_sim, subpops_select_real, assume_equal = T, LD_info = T, focal_pop_sim = NULL, focal_pop_real = NULL, edge_length_sim, edge_length_real){
  if (assume_equal){
    # make SS_sim
    Fc = summarise_Fc(sim_out$Fc_sim, S_select, subpops_select_sim, gens_select)
    r2 = summarise_r2(sim_out$r2_sim, S_select, subpops_select_sim, gens_select)
    rirj = summarise_rirj(sim_out$rirj_sim, S_select, subpops_select_sim, gens_select, edge_length_sim)
    FcCor = summarise_FcCor(sim_out$FcCor_sim, S_select, subpops_select_sim, gens_select, edge_length_sim)
    
    print(dim(FcCor))
    print(dim(rirj))
    
    FcCor2d = make_2d(FcCor)
    rirj2d = make_2d(rirj)
    
    if (LD_info) SS_sim = cbind(t(Fc), t(r2), FcCor2d[,], rirj2d[,] )
    if (!LD_info) SS_sim = cbind(t(Fc), FcCor2d[,])
    
    ## make SS_real
    SS_real_list = list()
    for (i in 1:length(real_out_list)){
      real_out = real_out_list[[i]]
      Fc_r = summarise_Fc(real_out$Fc_sim, S_select, subpops_select_real, gens_select)
      r2_r = summarise_r2(real_out$r2_sim, S_select, subpops_select_real, gens_select)
      rirj_r = summarise_rirj(real_out$rirj_sim, S_select, subpops_select_real, gens_select, edge_length_real)
      FcCor_r = summarise_FcCor(real_out$FcCor_sim, S_select, subpops_select_real, gens_select, edge_length_real)
      FcCor2d_r = make_2d(FcCor_r)
      rirj2d_r = make_2d(rirj_r)
      
      SS_real = cbind(t(Fc_r), t(r2_r), FcCor2d_r[,], rirj2d_r[,] )
      
      if (LD_info) SS_real = cbind(t(Fc_r), t(r2_r), FcCor2d_r[,], rirj2d_r[,] )
      if (!LD_info)SS_real = cbind(t(Fc_r), FcCor2d_r[,])
      
      SS_real_list[[i]] = SS_real
    }
  }
  
  else{
    focal_pop_sim_str = c(paste0('sp', focal_pop_sim))
    focal_pop_real_str = c(paste0('sp', focal_pop_real))
    # make SS_sim
    Fc = summarise_Fc(sim_out$Fc_sim, S_select, focal_pop_sim, gens_select)
    r2 = summarise_r2(sim_out$r2_sim, S_select, focal_pop_sim, gens_select)
    rirj = summarise_rirj(sim_out$rirj_sim, S_select, subpops_select_sim, gens_select, focal = focal_pop_sim, edge_length = edge_length_sim)
    FcCor = summarise_FcCor(sim_out$FcCor_sim, S_select, subpops_select_sim, gens_select, focal = focal_pop_sim, edge_length = edge_length_sim)
    
    print(dim(FcCor))
    print(dim(rirj))
    
    FcCor2d = make_2d(FcCor)
    rirj2d = make_2d(rirj)
    
    if (LD_info) SS_sim = cbind(t(Fc), t(r2), FcCor2d[,], rirj2d[,] )
    if (!LD_info) SS_sim = cbind(t(Fc), FcCor2d[,])
    
    ## make SS_real
    SS_real_list = list()
    for (i in 1:length(real_out_list)){
      real_out = real_out_list[[i]]
      Fc_r = summarise_Fc(real_out$Fc_sim, S_select, focal_pop_real, gens_select)
      r2_r = summarise_r2(real_out$r2_sim, S_select, focal_pop_real, gens_select)
      rirj_r = summarise_rirj(real_out$rirj_sim, S_select, subpops_select_real, gens_select, focal = focal_pop_real, edge_length = edge_length_real)
      FcCor_r = summarise_FcCor(real_out$FcCor_sim, S_select, subpops_select_real, gens_select, focal = focal_pop_real, edge_length = edge_length_real)
      
      FcCor2d_r = make_2d(FcCor_r)
      rirj2d_r = make_2d(rirj_r)
      
      SS_real = cbind(t(Fc_r), t(r2_r), FcCor2d_r[,], rirj2d_r[,] )
      
      if (LD_info) SS_real = cbind(t(Fc_r), t(r2_r), FcCor2d_r[,], rirj2d_r[,] )
      if (!LD_info)SS_real = cbind(t(Fc_r), FcCor2d_r[,])
      
      SS_real_list[[i]] = SS_real
    }
  }
  return(list(sim = SS_sim, real = SS_real_list))
}


ceil = function(val, step = 1){
  if( val == 0 ) return (val)

  else {sign = -1; val = -val}
  
  temp = 0
  while (temp < val){
    temp = temp + step
  }
  return (temp*sign)
  
}


ceil_comp = function(val, step = 1){
  if( val == 0 ) return (val)
  if (val < 0){
    temp = 0
    while (temp>val){
      temp = temp - step
    }
    return(temp+step)
  }
  if (val>0){
    temp = 0
    while (temp<val){
      temp = temp + step
    }
    return(temp)
  }
}

floor_comp = function(val, step = 1){
  if( val == 0 ) return (val)
  if (val < 0){
    temp = 0
    while (temp>val){
      temp = temp - step
    }
    return(temp)
  }
  if (val>0){
    temp = 0
    while (temp<val){
      temp = temp + step
    }
    return(temp-step)
  }
}

ceil = function(nums, step){
  return(c(floor_comp(nums[1], step), ceil_comp(nums[2], step)))
}

ceil(c(-0.73,1.43), 0.2 )

ceil_comp(0.73, 0.2)


list_range = function(l, idx2,idx3,prior = F){
  if (prior){
    r1 = list_range_sub(l, idx2,idx3)
    r2 = list_range_sub(l[-length(l)], idx2,idx3)
    
    low = r2[1]
    high = min(r2[2]*1.1, r1[2]*1.02 )
    
    if (r1[2] > r2[2]*1.5) return(r2)
    else ( return(r1) )
    return(r2)
  }
  
  return(list_range_sub(l, idx2,idx3))
}

list_range_sub = function(l, idx2,idx3){
  v = c()
  for (e in l){
    
    v = c(v, e[,idx2,idx3])
  }
  
  x = range(v)
  
  return(x)
}

lines_from_list = function(l,x2,x3,y2,y3,ltys,cols, div = F, d2=NA,d3=NA){
  ltys = rep(ltys,100)
  cols = rep(cols,100)
  i=1
  
  if (div == F){
    for (e in l){
      
      
      lines(e[,x2,x3], e[,y2,y3], lty = ltys[[i]], col = cols[[i]])
      i=i+1
    }
  }
  
  if (div == T){
    for (e in l){
      
      
      lines(e[,x2,x3], e[,y2,y3]/e[,d2,d3], lty = ltys[[i]], col = cols[[i]])
      i=i+1
    }
  }
  
}

plot_compare_N = function(plot_list, legend, cols = 'nope', ltys = 'nope'){
  if (cols == 'nope'){ print('default cols'); cols = c('black')}
  if (ltys == 'nope') {print('default ltys');ltys = 1:length(plot_list)}
  
  # posterior means
  par(mfrow = c(2,1), mai = c(1,1,.3,.3))
  plot(NA, xlim = list_range(plot_list,1,1), ylim = list_range(plot_list, 2,1), xlab = 'N', ylab = expression(hat(N)))
  lines_from_list(plot_list, 1,1,2,1,ltys = 1:5, cols =c('black'))
  
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  plot(NA, xlim = list_range(plot_list,1,1), ylim = list_range(plot_list, 2,2), xlab = 'N', ylab = expression(hat(m)))
  lines_from_list(plot_list, 1,1,2,2,ltys = 1:5, cols =c('black'))
  
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  
  # posterior means
  plot_list_div = plot_list
  for ( i in 1:length(plot_list)){
    #plot_list_div[[i]] = abind(plot_list_div[[i]], (plot_list_div[[i]][,2,]-plot_list_div[[i]][,1,])/plot_list_div[[i]][,1,], along = 2)
    plot_list_div[[i]] = abind(plot_list_div[[i]], log2(plot_list_div[[i]][,2,]/plot_list_div[[i]][,1,]), along = 2)
  }
  
  par(mfrow = c(2,1), mai = c(1,1,.3,.3))
  plot(NA, xlim = list_range(plot_list_div,1,1), ylim = list_range(plot_list_div, 4,1), xlab = 'N', ylab = expression(hat(N)/N)); abline(h=0, col = 'red')
  lines_from_list(plot_list_div, 1,1,4,1,ltys = 1:5, cols =c('black'))
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  plot(NA, xlim = list_range(plot_list_div,1,1), ylim = list_range(plot_list_div, 4,2), xlab = 'N', ylab = expression(hat(m)/m)); abline(h=0, col = 'red')
  lines_from_list(plot_list_div, 1,1,4,2,ltys = 1:5, cols =c('black'))
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  
  # posterior variance
  par(mfrow = c(2,1), mai = c(1,1,.3,.3))
  plot(NA, xlim = list_range(plot_list,1,1), ylim = list_range(plot_list, 3,1), xlab = 'N', ylab = expression(Var(hat(N))))
  lines_from_list(plot_list, 1,1,3,1,ltys = 1:5, cols =c('black'))
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  plot(NA, xlim = list_range(plot_list,1,1), ylim = list_range(plot_list, 3,2), xlab = 'N', ylab = expression(Var(hat(m))))
  lines_from_list(plot_list, 1,1,3,2,ltys = 1:5, cols =c('black'))
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  
  
  plot_list_CV = plot_list
  for ( i in 1:length(plot_list)){
    plot_list_CV[[i]] = abind(plot_list_CV[[i]], sqrt(plot_list_CV[[i]][,3,])/plot_list_div[[i]][,2,], along = 2)
  }
  
  par(mfrow = c(2,1), mai = c(1,1,.3,.3))
  plot(NA, xlim = list_range(plot_list_CV,1,1), ylim = list_range(plot_list_CV, 4,1), xlab = 'N', ylab = expression(CV(hat(N))))
  lines_from_list(plot_list_CV, 1,1,4,1,ltys = 1:5, cols =c('black'))
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  plot(NA, xlim = list_range(plot_list_CV,1,1), ylim = list_range(plot_list_CV, 4,2), xlab = 'N', ylab = expression(CV(hat(m))))
  lines_from_list(plot_list_CV, 1,1,4,2,ltys = 1:5, cols =c('black'))
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  return()
  
}




plot_compare_N_neat = function(plot_list, legend, cols = 'nope', ltys = 'nope'){
  if (cols == 'nope'){ print('default cols'); cols = c('black')}
  if (ltys == 'nope') {print('default ltys');ltys = 1:length(plot_list)}
  
  # if (cols == 'nope'){ print('default cols'); cols = colorRampPalette(colors = c('red', 'blue'))(length(plot_list)) }
  # else {   cols = colorRampPalette(colors = cols)(length(plot_list)) }
  # 
  # if (ltys == 'nope') {print('default ltys'); ltys = 1}
  
  # posterior means
  plot_list_div = plot_list
  for ( i in 1:length(plot_list)){
    #plot_list_div[[i]] = abind(plot_list_div[[i]], (plot_list_div[[i]][,2,]-plot_list_div[[i]][,1,])/plot_list_div[[i]][,1,], along = 2)
    plot_list_div[[i]] = abind(plot_list_div[[i]], log2(plot_list_div[[i]][,2,]/plot_list_div[[i]][,1,]), along = 2)
  }
  
  
  plot_list_CV = plot_list
  for ( i in 1:length(plot_list)){
    plot_list_CV[[i]] = abind(plot_list_CV[[i]], sqrt(plot_list_CV[[i]][,3,])/plot_list_div[[i]][,2,], along = 2)
  }
  
  
  #pdf('../test.pdf', width = 10, height = 10)
  par(mfrow = c(2,2), mai = c(.5,.5,.2,.2), oma = c(4, 1, 0, 0))
  
  plot(NA, bty = 'l', xaxs = 'i', yaxs = 'i', xlim = list_range(plot_list_div,1,1), ylim = ceil(list_range(plot_list_div, 4,1),.5), xlab = 'N', ylab = expression(log[2](hat(N)/N)), xaxt = 'n'); abline(h=0, col = 'black', lw = 1); axis(1, at = c(400,seq(5000,20000,5000)))
  lines_from_list(plot_list_div, 1,1,4,1,ltys = ltys, cols = cols)
  fig_label('A', cex = 2)
  
  
  plot(NA,bty = 'l', xaxs = 'i', yaxs = 'i', xlim = list_range(plot_list_div,1,1), ylim = ceil(list_range(plot_list_div, 4,2),.5), xlab = 'N', ylab = expression(log[2](hat(m)/m)), xaxt = 'n'); abline(h=0, col = 'black', lw = 1); axis(1, at = c(400,seq(5000,20000,5000)))
  lines_from_list(plot_list_div, 1,1,4,2,ltys = ltys, cols = cols)
  fig_label('B', cex = 2)
  

  
  plot(NA, bty = 'l', xaxs = 'i', yaxs = 'i',xlim = list_range(plot_list_CV,1,1), ylim = ceil(c(0,list_range(plot_list_CV, 4,1)[2]),.1), xlab = 'N', ylab = expression(CV(hat(N))), xaxt = 'n'); axis(1, at = c(400,seq(5000,20000,5000)))
  lines_from_list(plot_list_CV, 1,1,4,1,ltys =ltys, cols =cols)
  fig_label('C', cex = 2)
  
  plot(NA,bty = 'l', xaxs = 'i', yaxs = 'i', xlim = list_range(plot_list_CV,1,1), ylim = ceil(c(0,list_range(plot_list_CV, 4,2)[2]),.1), xlab = 'N', ylab = expression(CV(hat(m))), xaxt = 'n'); axis(1, at = c(400,seq(5000,20000,5000)))
  lines_from_list(plot_list_CV, 1,1,4,2,ltys =ltys, cols =cols)
  fig_label('D', cex = 2)
  
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend('bottom', legend = legend, lty = ltys, col = cols,horiz = TRUE, inset = rep(0.02,4))
  #graphics.off()
  
  
  
  
  
  
  
  return()
  
}


plot_compare_m = function(plot_list, legend, cols = 'nope', ltys = 'nope'){
  if (cols == 'nope'){ print('default cols'); cols = c('black')}
  if (ltys == 'nope') {print('default ltys');ltys = 1:length(plot_list)}

  # if (cols == 'nope'){ print('default cols'); cols = colorRampPalette(colors = c('red', 'blue'))(length(plot_list)) }
  # else {   cols = colorRampPalette(colors = cols)(length(plot_list)) }
  # 
  # if (ltys == 'nope') {print('default ltys'); ltys = 1}
  
  
  
  # posterior means
  par(mfrow = c(2,1), mai = c(1,1,.3,.3))
  plot(NA, xlim = list_range(plot_list,1,2), ylim = list_range(plot_list, 2,1), xlab = 'm', ylab = expression(hat(N)))
  lines_from_list(plot_list, 1,2,2,1,ltys = ltys, cols = cols)
  
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  plot(NA, xlim = list_range(plot_list,1,2), ylim = list_range(plot_list, 2,2), xlab = 'm', ylab = expression(hat(m)))
  lines_from_list(plot_list, 1,2,2,2,ltys = ltys, cols = cols)
  
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  
  # posterior means
  plot_list_div = plot_list
  for ( i in 1:length(plot_list)){
    #plot_list_div[[i]] = abind(plot_list_div[[i]], (plot_list_div[[i]][,2,]-plot_list_div[[i]][,1,])/plot_list_div[[i]][,1,], along = 2)
    plot_list_div[[i]] = abind(plot_list_div[[i]], log2(plot_list_div[[i]][,2,]/plot_list_div[[i]][,1,]), along = 2)
    
  }
  par(mfrow = c(2,1), mai = c(1,1,.3,.3))
  plot(NA, xlim = list_range(plot_list_div,1,2), ylim = list_range(plot_list_div, 4,1), xlab = 'm', ylab = expression(hat(N)/N)); abline(h=0, col = 'black')
  lines_from_list(plot_list_div, 1,2,4,1,ltys = ltys, cols = cols)
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  plot(NA, xlim = list_range(plot_list_div,1,2), ylim = list_range(plot_list_div, 4,2), xlab = 'm', ylab = expression(hat(m)/m)); abline(h=0, col = 'black')
  lines_from_list(plot_list_div, 1,2,4,2,ltys = ltys, cols = cols)
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  
  # posterior variance
  par(mfrow = c(2,1), mai = c(1,1,.3,.3))
  plot(NA, xlim = list_range(plot_list,1,2), ylim = list_range(plot_list, 3,1), xlab = 'm', ylab = expression(Var(hat(N))))
  lines_from_list(plot_list, 1,2,3,1,ltys = ltys, cols = cols)
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  plot(NA, xlim = list_range(plot_list,1,2), ylim = list_range(plot_list, 3,2), xlab = 'm', ylab = expression(Var(hat(m))))
  #plot(NA, xlim = list_range(plot_list,1,2), ylim = c(0,0.08), xlab = 'm', ylab = expression(Var(hat(m))))
  
  lines_from_list(plot_list, 1,2,3,2,ltys = ltys, cols = cols)
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  
  
  plot_list_CV = plot_list
  for ( i in 1:length(plot_list)){
    plot_list_CV[[i]] = abind(plot_list_CV[[i]], sqrt(plot_list_CV[[i]][,3,])/plot_list_div[[i]][,2,], along = 2)
  }
  
  par(mfrow = c(2,1), mai = c(1,1,.3,.3))
  plot(NA, xlim = list_range(plot_list_CV,1,2), ylim = list_range(plot_list_CV, 4,1), xlab = 'm', ylab = expression(CV(hat(N))))
  lines_from_list(plot_list_CV, 1,2,4,1,ltys = ltys, cols = cols)
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  plot(NA, xlim = list_range(plot_list_CV,1,2), ylim = list_range(plot_list_CV, 4,2), xlab = 'm', ylab = expression(CV(hat(m))))
  lines_from_list(plot_list_CV, 1,2,4,2,ltys = ltys, cols = cols)
  legend('topright', legend = legend, lty = ltys, col = cols)
  
  
  
  
  
  
  
  
  return()
  
}

plot_compare_m_neat = function(plot_list, legend, cols = 'nope', ltys = 'nope'){
  if (cols == 'nope'){ print('default cols'); cols = c('black')}
  if (ltys == 'nope') {print('default ltys');ltys = 1:length(plot_list)}
  
  # if (cols == 'nope'){ print('default cols'); cols = colorRampPalette(colors = c('red', 'blue'))(length(plot_list)) }
  # else {   cols = colorRampPalette(colors = cols)(length(plot_list)) }
  # 
  # if (ltys == 'nope') {print('default ltys'); ltys = 1}
  
  

  
  # posterior means
  plot_list_div = plot_list
  for ( i in 1:length(plot_list)){
    #plot_list_div[[i]] = abind(plot_list_div[[i]], (plot_list_div[[i]][,2,]-plot_list_div[[i]][,1,])/plot_list_div[[i]][,1,], along = 2)
    plot_list_div[[i]] = abind(plot_list_div[[i]], log2(plot_list_div[[i]][,2,]/plot_list_div[[i]][,1,]), along = 2)
    
  }
  
  
  plot_list_CV = plot_list
  for ( i in 1:length(plot_list)){
    plot_list_CV[[i]] = abind(plot_list_CV[[i]], sqrt(plot_list_CV[[i]][,3,])/plot_list_div[[i]][,2,], along = 2)
  }
  print('1')
  #pdf('../test.pdf', width = 10, height = 10)
  par(mfrow = c(2,2), mai = c(.7,.7,.2,.2), oma = c(4, 1, 0, 0))
  
  
  plot(NA,bty = 'l',log = 'x', xaxs = 'i', yaxs = 'i', xlim = c(0.01,list_range(plot_list_div,1,2)[2]), ylim = ceil(list_range(plot_list_div, 4,1), .5), xlab = 'm', ylab = expression(log[2](hat(N)/N)), xaxt = 'n'); abline(h=0, col = 'black'); axis(1, at = c(0.01,0.02,0.04,0.1,0.2,0.4,0.8))
  lines_from_list(plot_list_div, 1,2,4,1,ltys = ltys, cols = cols)
  fig_label_logx('A', cex = 2)
  
  plot(NA,bty = 'l',log = 'x', xaxs = 'i', yaxs = 'i', xlim = c(0.01,list_range(plot_list_div,1,2)[2]), ylim = ceil(list_range(plot_list_div, 4,2),.5), xlab = 'm', ylab = expression(log[2](hat(m)/m)), xaxt = 'n'); abline(h=0, col = 'black'); axis(1, at = c(0.01,0.02,0.04,0.1,0.2,0.4,0.8))
  lines_from_list(plot_list_div, 1,2,4,2,ltys = ltys, cols = cols)
  fig_label_logx('B', cex = 2)
  


  
  plot(NA,bty = 'l',log = 'x', xaxs = 'i', yaxs = 'i', xlim = c(0.01,list_range(plot_list_CV,1,2)[2]), ylim = ceil(c(0,list_range(plot_list_CV, 4,1)[2]),.1), xlab = 'm', ylab = expression(CV(hat(N))), xaxt = 'n'); axis(1, at = c(0.01,0.02,0.04,0.1,0.2,0.4,0.8))
  lines_from_list(plot_list_CV, 1,2,4,1,ltys = ltys, cols = cols)
  fig_label_logx('C', cex = 2)
  
  plot(NA,bty = 'l',log = 'x', xaxs = 'i', yaxs = 'i', xlim = c(0.01,list_range(plot_list_CV,1,2)[2]), ylim = ceil(c(0,list_range(plot_list_CV, 4,2)[2]),.1), xlab = 'm', ylab = expression(CV(hat(m))), xaxt = 'n'); axis(1, at = c(0.01,0.02,0.04,0.1,0.2,0.4,0.8))
  lines_from_list(plot_list_CV, 1,2,4,2,ltys = ltys, cols = cols)
  fig_label_logx('D', cex = 2)
  
  
  
  
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend('bottom', legend = legend, lty = ltys, col = cols,horiz = TRUE, inset = rep(0.02,4))
  #graphics.off()

  
  return()
  
}


plot_compare_m_neat_Nm = function(plot_list, legend, cols = 'nope', ltys = 'nope'){
  if (cols == 'nope'){ print('default cols'); cols = rep('black',5)}
  if (ltys == 'nope') {print('default ltys');ltys = 1:5}
  
  cx = 1; cy = 1

  # posterior means
  plot_list_div = plot_list
  for ( i in 1:length(plot_list)){
    #plot_list_div[[i]] = abind(plot_list_div[[i]], (plot_list_div[[i]][,2,]-plot_list_div[[i]][,1,])/plot_list_div[[i]][,1,], along = 2)
    #print(i)
    #print(plot_list_div[[i]][,,1])

    plot_list_div[[i]] = cbind(plot_list_div[[i]][,1],log2(plot_list[[i]][,2]/plot_list[[i]][,1]))
  }
  plot_list_CV = plot_list
  for ( i in 1:length(plot_list)){

    plot_list_CV[[i]] = cbind(plot_list_CV[[i]][,1], sqrt(plot_list[[i]][,3])/plot_list[[i]][,1])
  }
  #pdf('../test.pdf', width = 10, height = 10)
  par(mfrow = c(1,2), mai = c(.9,1.3,.2,.2), oma = c(1, 0, 0, 0))
  
  plot(NA,bty = 'l',log = 'x', xaxs = 'i', yaxs = 'i', xlim = c(0.01,0.8), ylim = c(.5,-.5), xlab = 'm', ylab = bquote('log'['2']*'('*frac(hat('N'['e']*'m'),'N'['e']*'m')*')'), xaxt = 'n', yaxt = 'n'); abline(h=0, col = 'black'); axis(1, at = c(0.01,0.02,0.04,0.1,0.2,0.4,0.8))
  axis(2, at = c( seq(-1,1,.5)), cex.axis = cy)
  #axis(2, at = c(seq(yran[1], yran[2], 1)))
  rug(seq(-1,1, .1), ticksize = -0.02, side = 2)
  
  i = 1
  for (l in plot_list_div){
    #print(l)

    lines(l[,1]/2000, l[,2], lty = ltys[i], col = cols[i])
    print(i); print(l[,2])
    i = i+1
  }
  fig_label_logx('A', cex = 2)
  
  plot(NA,bty = 'l',log = 'x', xaxs = 'i', yaxs = 'i', xlim = c(0.01,0.8), ylim = c(-0,.8), xlab = 'm', ylab = bquote('CV'(hat('N'['e']*'m'))), xaxt = 'n'); abline(h=0, col = 'black'); axis(1, at = c(0.01,0.02,0.04,0.1,0.2,0.4,0.8))
  i = 1
  for (l in plot_list_CV){
    lines(l[,1]/2000, l[,2], lty = ltys[i], col = cols[i])
    i = i+1
  }
  fig_label_logx('B', cex = 2)
  
  
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend('bottom', legend = legend, lty = ltys, col = cols,horiz = TRUE, inset = rep(0.02,4))
  #graphics.off()
  
  
  return()
  
}


invert_N = function(results){
  results[,1:2,1,] = 2/results[,1:2,1,]
  return(results)
}


logtix = c(2,5,10); logtix = logtix * rep(c(1e-2, 1e-1, 1e0, 1e1), each = 3); logtix.char = as.character(logtix)

plot_compare_N_m_neat = function(plot_list_N, plot_list_m, legend, cols, ltys, prior = F, fname = '../test.pdf'){
  
  last = 1e10
  if (prior) last = length(plot_list_N)
  
  ##### N ##############
  
  # posterior means
  plot_list_div = plot_list_N
  for ( i in 1:length(plot_list_N)){
    #plot_list_div[[i]] = abind(plot_list_div[[i]], (plot_list_div[[i]][,2,]-plot_list_div[[i]][,1,])/plot_list_div[[i]][,1,], along = 2)
    plot_list_div[[i]] = abind(plot_list_div[[i]], log2(plot_list_div[[i]][,2,]/plot_list_div[[i]][,1,]), along = 2)
  }
  
  
  plot_list_CV = plot_list_N
  for ( i in 1:length(plot_list_N)){
    plot_list_CV[[i]] = abind(plot_list_CV[[i]], sqrt(plot_list_CV[[i]][,3,])/plot_list_CV[[i]][,1,], along = 2)
  }
  
  
  # plot_list_CV = plot_list_N
  # for ( i in 1:length(plot_list_N)){
  #   plot_list_CV[[i]] = abind(plot_list_CV[[i]], abs(plot_list_CV[[i]][,1,]-plot_list_CV[[i]][,2,])/plot_list_CV[[i]][,1,], along = 2)
  # }
  
  
  pdf(fname, width = 7, height = 8.55, pointsize = 9)
  par(mfrow = c(4,2), mai = c(.5,.5,.05,.2), oma = c(4, .2, .2, .3))
  
  cx = .95; cy = .90
  # N, Nhat bias
  plot(NA, bty = 'l', xaxs = 'i', yaxs = 'i', xlim = list_range(plot_list_div,1,1), ylim = ceil(list_range(plot_list_div, 4,1, prior),.5), xlab = expression(N[e]), ylab = expression(log[2](hat(N[e])/N[e])), xaxt = 'n', yaxt = 'n'); abline(h=0, col = 'black', lw = 1)
  axis(1, at = c(400,seq(5000,20000,5000)), cex.axis=cx); rug(seq(1000,20000,1000), ticksize = -0.02, side = 1)
  yran = ceil(list_range(plot_list_div, 4,1),.5)
  axis(2, at = c(seq(yran[1], yran[2], .5)), cex.axis = cy); rug(seq(yran[1], yran[2], .1), ticksize = -0.02, side = 2)
  
  lines_from_list(plot_list_div, 1,1,4,1, ltys = ltys, cols = cols)
  fig_label('A', cex = 2, region = 'figure')
  
  
  # N, mhat bias
  plot(NA, bty = 'l', xaxs = 'i', yaxs = 'i', xlim = list_range(plot_list_div,1,1), ylim = ceil(list_range(plot_list_div, 4,2, prior),.5), xlab = expression(N[e]), ylab = expression(log[2](hat(m)/m)), xaxt = 'n', yaxt = 'n'); abline(h=0, col = 'black', lw = 1)
  lines_from_list(plot_list_div, 1,1,4,2,ltys = ltys, cols = cols)
  fig_label('B', cex = 2, region = 'figure')
  axis(1, at = c(400,seq(5000,20000,5000)), cex.axis = cy); rug(seq(1000,20000,1000), ticksize = -0.02, side = 1)
  yran = ceil(list_range(plot_list_div, 4,2, prior),1)
  #axis(2, at = c( seq(0,3,.5)), cex.axis = cy)
  axis(2, at = c( seq(-1,5,.5)), cex.axis = cy)
  #axis(2, at = c(seq(yran[1], yran[2], 1)))
  rug(seq(yran[1], yran[2], .1), ticksize = -0.02, side = 2)
  
  

  # N, Nhat CV
  plot(NA,  bty = 'l', xaxs = 'i', yaxs = 'i', log = 'y',xlim = list_range(plot_list_CV,1,1), ylim = c(0.04,1 + 1e-100*ceil(list_range(plot_list_CV, 4,1, prior),.1)[2]), xlab = expression(N[e]), ylab = expression(CV(hat(N[e]))), xaxt = 'n', yaxt = 'n')
  axis(1, at = c(400,seq(5000,20000,5000)), cex.axis = cx); rug(seq(1000,20000,1000), ticksize = -0.02, side = 1)
  yran = ceil(list_range(plot_list_CV, 4,1, prior),.5)
  #axis(2, at = c(seq(yran[1], yran[2], 1)))
  axis(2, at=logtix, labels = logtix.char, cex.axis = cy)
  
  lines_from_list(plot_list_CV, 1,1,4,1,ltys =ltys, cols =cols)
  fig_label_logy('C', cex = 2, region = 'figure')
  
  
  # N, mhat CV
  plot(NA, bty = 'l', xaxs = 'i', yaxs = 'i', log = 'y',xlim = list_range(plot_list_CV,1,1), ylim  = c(0.04, 5+1e-100*ceil(list_range(plot_list_CV, 4,2, prior),.1)[2]), xlab = expression(N[e]), ylab = expression(CV(hat(m))), xaxt = 'n', yaxt = 'n')
  
  axis(1, at = c(400,seq(5000,20000,5000)), cex.axis=cx); rug(seq(1000,20000,1000), ticksize = -0.02, side = 1)
  axis(2, at=logtix, labels = logtix.char, cex.axis = cy)
  
  lines_from_list(plot_list_CV, 1,1,4,2,ltys =ltys, cols =cols)
  fig_label_logy('D', cex = 2, region = 'figure')
  
  
  
  
  ##### m ##############
  
  
  plot_list_div = plot_list_m
  for ( i in 1:length(plot_list_m)){
    #plot_list_div[[i]] = abind(plot_list_div[[i]], (plot_list_div[[i]][,2,]-plot_list_div[[i]][,1,])/plot_list_div[[i]][,1,], along = 2)
    plot_list_div[[i]] = abind(plot_list_div[[i]], log2(plot_list_div[[i]][,2,]/plot_list_div[[i]][,1,]), along = 2)
    
  }
  
  
  plot_list_CV = plot_list_m
  for ( i in 1:length(plot_list_m)){
    plot_list_CV[[i]] = abind(plot_list_CV[[i]], sqrt(plot_list_CV[[i]][,3,])/plot_list_CV[[i]][,1,], along = 2)
  }
  
  # plot_list_CV = plot_list_m
  # for ( i in 1:length(plot_list_m)){
  #   plot_list_CV[[i]] = abind(plot_list_CV[[i]], abs(plot_list_CV[[i]][,1,]-plot_list_CV[[i]][,2,])/plot_list_CV[[i]][,1,], along = 2)
  # }
  
  plot(NA, bty = 'l', log = 'x', xaxs = 'i', yaxs = 'i', xlim = c(0.01,list_range(plot_list_div,1,2)[2]), ylim = ceil(list_range(plot_list_div, 4,1, prior), .5), xlab = 'm', ylab = expression(log[2](hat(N[e])/N[e])), xaxt = 'n', yaxt = 'n'); abline(h=0, col = 'black')
  axis(1, at = c(0.01,0.02,0.04,0.1,0.2,0.4,0.8), cex.axis = cx); rug(diff(c(0.01,0.02,0.04,0.1,0.2,0.4,0.8))/2 + c(0.01,0.02,0.04,0.1,0.2,0.4,0.8)[-7], ticksize = -0.02, side = 1)
  yran = ceil(list_range(plot_list_div, 4,1, prior),1)
  axis(2, at = c(seq(yran[1], yran[2], 1)), cex.axis = cy); rug(seq(yran[1], yran[2], .1), ticksize = -0.02, side = 2)
  
  lines_from_list(plot_list_div, 1,2,4,1,ltys = ltys, cols = cols)
  fig_label_logx('E', cex = 2, region = 'figure')
  
  plot(NA, bty = 'l',log = 'x', xaxs = 'i', yaxs = 'i', xlim = c(0.01,list_range(plot_list_div,1,2)[2]), ylim = ceil(list_range(plot_list_div, 4,2, prior),.5), xlab = 'm', ylab = expression(log[2](hat(m)/m)), xaxt = 'n', yaxt='n'); abline(h=0, col = 'black')
  axis(1, at = c(0.01,0.02,0.04,0.1,0.2,0.4,0.8), cex.axis = cx); rug(diff(c(0.01,0.02,0.04,0.1,0.2,0.4,0.8))/2 + c(0.01,0.02,0.04,0.1,0.2,0.4,0.8)[-7], ticksize = -0.02, side = 1)
  
  yran = ceil(list_range(plot_list_div, 4,2, prior),.5)
  axis(2, seq(-2,2), cex.axis = cy); rug(seq(yran[1], yran[2], .1), ticksize = -0.02, side = 2)
  lines_from_list(plot_list_div, 1,2,4,2,ltys = ltys, cols = cols)
  fig_label_logx('F', cex = 2, region = 'figure')
  
  
  
  #
  plot(NA, bty = 'l',log = 'xy', xaxs = 'i', yaxs = 'i', xlim = c(0.01,list_range(plot_list_CV,1,2)[2]), ylim = c(0.02, 5+1e-100*ceil(list_range(plot_list_CV, 4,1, prior),.1)[2]), xlab = 'm', ylab = expression(CV(hat(N[e]))), xaxt = 'n', yaxt = 'n')
  axis(1, at = c(0.01,0.02,0.04,0.1,0.2,0.4,0.8), cex.axis = cx); rug(diff(c(0.01,0.02,0.04,0.1,0.2,0.4,0.8))/2 + c(0.01,0.02,0.04,0.1,0.2,0.4,0.8)[-7], ticksize = -0.02, side = 1)
  axis(2, at=logtix, labels = logtix.char, cex.axis = cy-0.02)
  
  lines_from_list(plot_list_CV, 1,2,4,1,ltys = ltys, cols = cols)
  fig_label_logxy('G', cex = 2, region = 'figure')
  
  plot(NA, bty = 'l',log = 'xy',xaxs = 'i', yaxs = 'i', xlim = c(0.01,list_range(plot_list_CV,1,2)[2]), ylim = c(0.04, 5+1e-100*ceil(list_range(plot_list_CV, 4,2, prior),.1)[2]), xlab = 'm', ylab = expression(CV(hat(m))), xaxt = 'n', yaxt = 'n')
  axis(1, at = c(0.01,0.02,0.04,0.1,0.2,0.4,0.8), cex.axis = cx); rug(diff(c(0.01,0.02,0.04,0.1,0.2,0.4,0.8))/2 + c(0.01,0.02,0.04,0.1,0.2,0.4,0.8)[-7], ticksize = -0.02, side = 1)
  axis(2, at=logtix, labels = logtix.char, cex.axis = cy)
  
  lines_from_list(plot_list_CV, 1,2,4,2,ltys = ltys, cols = cols)
  fig_label_logxy('H', cex = 2, region = 'figure')
  
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend('bottom', legend = legend, lty = ltys, col = cols,horiz = TRUE, inset = rep(0.02,4))
  graphics.off()
  
  return()
  
}


summarise_Fc = function(Fc_out, S, subpops, timepoints){
  if (length(timepoints) == 1){
    return(array(NA, dim = c(0, dim(Fc_out)[6] )))
  }
  
  S_filter = paste0('S', S)
  subpop_filter = paste0('sp', subpops)
  timepoint_filter = paste0('g', timepoints)
  gens = as.numeric(substr(dimnames(Fc_out)[[2]], 2,100)); intervals = get_all_intervals(timepoints); interval_counts = get_intervals_counts(timepoints)
  print(subpop_filter)
  
  Fc_filtered = Fc_out[S_filter, timepoint_filter, timepoint_filter, subpop_filter, , ,drop=F]
  
  Fc_final = array(0, dim = c(length(intervals), dim(Fc_out)[6] ))
  
  
  dimnames(Fc_final)[[1]] = as.list(paste0('sep', intervals))
  
  idx_gens_passed = which(gens %in% timepoints)
  for (idx_g_1 in idx_gens_passed){
    g1 = gens[idx_g_1]; g1str = paste0('g', g1)
    for (idx_g_2 in idx_gens_passed){
      if (idx_g_2 > idx_g_1){
        g2 = gens[idx_g_2];  g2str = paste0('g', g2);
        Fc_final[paste0('sep', g2-g1),] <- apply(Fc_filtered[S_filter, g1str, g2str, subpop_filter, , ,drop=F],6,mean) + Fc_final[paste0('sep', g2-g1),]
        
      }
    }
  }
  
  for (sep in intervals){
    Fc_final[paste0('sep', sep),] = Fc_final[paste0('sep', sep),] / interval_counts[1+which(intervals == sep)]
  }
  return(Fc_final)
}

summarise_r2 = function(r2_out, S, subpops, timepoints){
  S_filter = paste0('S', S)
  subpop_filter = paste0('sp', subpops)
  timepoint_filter = paste0('g', timepoints)
  gens = as.numeric(substr(dimnames(r2_out)[[2]], 2,100))
  recom_rates = as.numeric(substr(dimnames(r2_out)[[3]], 2,100)) 
  print(subpop_filter)
  
  r2_filtered = r2_out[S_filter, timepoint_filter, , subpop_filter, ,drop=F]
  
  r2_final = array(0, dim = c(length(recom_rates), dim(r2_out)[5] ))
  dimnames(r2_final)[[1]] =as.list(paste0('c', recom_rates))
  
  idx_gens_passed = which(gens %in% timepoints)
  for (recom_rate in recom_rates){
    recom_str = paste0('c', recom_rate)
    print(dim(r2_filtered[S_filter, timepoint_filter,recom_str, subpop_filter, ,drop=F]))
    r2_final[recom_str,] = apply(r2_filtered[S_filter, timepoint_filter,recom_str, subpop_filter, ,drop=F],5,mean)
  }
  
  return(r2_final)
}

summarise_FcCor = function(FcCor_out, S, subpops, timepoints,  edge_length, focal = F){
  
  
  S_filter = paste0('S', S)
  subpop_filter = paste0('sp', subpops)
  timepoint_filter = paste0('g', timepoints)
  gens = as.numeric(substr(dimnames(FcCor_out)[[2]], 2,100)); intervals = get_all_intervals(timepoints); interval_counts = get_intervals_counts(timepoints)
  
  if (length(timepoints) == 1 | length(subpops) == 1){
    arr = array(NA, dim = c(length(intervals), length(subpops)-1, dim(FcCor_out)[7]))
    return(arr)
  }
  
  
  adj = gen_adj(subpops, edge_length)
  if (focal){
    for (i in 1:length(adj)) {
      if (dim(adj[[i]])[1]>1) {adj[[i]] = adj[[i]][apply(adj[[i]][,], 1, function(x){focal %in% x}),, drop = F]}
      else if (dim(adj[[i]])[1]==1 & focal %in% adj[[i]]) {adj[[i]] = adj[[i]]}
      else {adj[[i]] = array(NA, dim = c(0,2))}
    }
    adj_temp = list()
    for (i in 1:length(adj)) {
      if (length(adj[[i]] > 0)) adj_temp[[i]] = adj[[i]]
    }
    adj = adj_temp
  }
  
  
  
  FcCor_filtered = FcCor_out[S_filter, timepoint_filter, timepoint_filter, subpop_filter, subpop_filter, , ,drop=F]
  
  FcCor_final = array(0, dim = c(length(intervals), length(adj), dim(FcCor_out)[7] ))
  dimnames(FcCor_final)[[1]] = as.list(paste0('sep', intervals))
  dimnames(FcCor_final)[[2]] = as.list(paste0('dist', 1:length(adj)))
  
  idx_gens_passed = which(gens %in% timepoints)
  for (dist in 1:length(adj)){
    
    for (pair_idx in 1: dim(adj[[dist]])[1]){
      for (idx_g_1 in idx_gens_passed){
        g1 = gens[idx_g_1]; g1str = paste0('g', g1)
        for (idx_g_2 in idx_gens_passed){
          if (idx_g_2 > idx_g_1){
            g2 = gens[idx_g_2];  g2str = paste0('g', g2);
            
            
            FcCor_final[paste0('sep', g2-g1), dist,] = FcCor_final[paste0('sep', g2-g1), dist,] + apply(FcCor_filtered[S_filter, g1str, g2str, paste0('sp',adj[[dist]][pair_idx,1]), paste0('sp',adj[[dist]][pair_idx,2]), , ,drop=F], 7, mean)
            
          }
        }
      }
    }
    FcCor_final[, dist,] = FcCor_final[, dist,] / dim(adj[[dist]])[1]
  }
  
  for (sep in intervals){
    FcCor_final[paste0('sep', sep),,] = FcCor_final[paste0('sep', sep),,] / interval_counts[1+which(intervals == sep)]
  }
  
  
  return(FcCor_final)
}

summarise_rirj = function(rirj_out, S, subpops, timepoints, edge_length, focal = F){
  S_filter = paste0('S', S)
  subpop_filter = paste0('sp', subpops)
  timepoint_filter = paste0('g', timepoints)
  gens = as.numeric(substr(dimnames(rirj_out)[[2]], 2,100))
  recom_rates = as.numeric(substr(dimnames(rirj_out)[[3]], 2,100))
  
  
  if (length(subpops) == 1){
    arr = array(NA, dim = c(length(recom_rates), 0, dim(rirj_out)[6]))
    dimnames(arr)[[1]] = as.list(paste0('c', recom_rates))
    return(arr)
  }
  
  adj = gen_adj(subpops, edge_length)
  if (focal){
    for (i in 1:length(adj)) {
      if (dim(adj[[i]])[1]>1) {adj[[i]] = adj[[i]][apply(adj[[i]][,], 1, function(x){focal %in% x}),, drop = F]}
      else if (dim(adj[[i]])[1]==1 & focal %in% adj[[i]]) {adj[[i]] = adj[[i]]}
      else {adj[[i]] = array(NA, dim = c(0,2))}
    }
    adj_temp = list()
    for (i in 1:length(adj)) {
      if (length(adj[[i]] > 0)) adj_temp[[i]] = adj[[i]]
    }
    adj = adj_temp
  }
  
  
  rirj_filtered = rirj_out[S_filter, timepoint_filter, , subpop_filter, subpop_filter, ,drop=F]
  
  rirj_final = array(0, dim = c(length(recom_rates), length(adj), dim(rirj_out)[6] ))
  dimnames(rirj_final)[[1]] = as.list(paste0('c', recom_rates))
  dimnames(rirj_final)[[2]] = as.list(paste0('dist', 1:length(adj)))
  
  
  for (d in 1:length(adj)){
    for (i in 1: dim(adj[[d]])[1]){
      rirj_final[,d,] = rirj_final[,d,] + apply(rirj_filtered[S_filter, timepoint_filter, , paste0('sp',adj[[d]][i,1]), paste0('sp',adj[[d]][i,2]), ,drop=F], c(3,6), mean)
    }
    rirj_final[,d,] = rirj_final[,d,] / dim(adj[[d]])[1]
  }
  return(rirj_final)
}

make_2d = function(three_d){
  if(dim(three_d)[2] == 0){
    return(array(NA, dim = c(dim(three_d)[3], 0)))
  }  
  
  if (dim(three_d)[1] == 0){
    return(array(NA, dim = c(dim(three_d)[3], 0)))
  }  
  
  d = dim(three_d)
  out = array(NA, dim = c(d[3], 0))
  for (i in 1:d[1]){
    for (j in 1:d[2]){
      out  = cbind(out, three_d[i,j,])
    }
  }
  return(out)
}

do_abc_PLS = function(SS_sim, SS_real, PLS_use, Nm_values, method = 'neuralnet', tol = .1, numnet = 10, sizenet = 5){
  
  ### construct some posteriors:
  
  N_draws = dim(SS_sim)[1]
  
  # scale data
  scaled <-scale(rbind(SS_sim, SS_real) )
  stat.sim.scaled <- scaled[-(N_draws+1),]
  stat.obs.scaled <- scaled[N_draws+1,]
  
  
  # do PLS
  
  ### simulated
  mod = pls::plsr(stat.sim.scaled[,1:2] ~ stat.sim.scaled, ncomp = PLS_use, center = T)
  #print(summary(mod))
  
  mod = pls::plsr(Nm_values ~ stat.sim.scaled, ncomp = PLS_use, center = T)
  
  #print(summary(mod))
  
  stat.sim.scaled_PLS = mod$scores
  
  ### real
  stat.obs.scaled_PLS = (stat.obs.scaled - colMeans(stat.sim.scaled)) %*%  mod$projection
  
  
  ### PLS abc
  out_PLS_full <- abc(target=stat.obs.scaled_PLS, param=Nm_values, sumstat=stat.sim.scaled_PLS, tol=tol, method = method, numnet = numnet, sizenet = sizenet)
  
  if (method == 'rejection'){
    out_PLS_full$adj.values = out_PLS_full$unadj.values
  }
  
  
  out_PLS_full = remove_outlier(out_PLS_full)
  
  
  
  return(out_PLS_full)
  
}

do_abc_PLS_Nm = function(SS_sim, SS_real, PLS_use, Nm_values, method = 'neuralnet', tol = .1, numnet = 10, sizenet = 5){
  
  ### construct some posteriors:
  
  N_draws = dim(SS_sim)[1]
  
  # scale data
  scaled <-scale(rbind(SS_sim, SS_real) )
  stat.sim.scaled <- scaled[-(N_draws+1),]
  stat.obs.scaled <- scaled[N_draws+1,]
  
  
  # do PLS
  
  ### simulated
  mod = pls::plsr(stat.sim.scaled[,1:2] ~ stat.sim.scaled, ncomp = PLS_use, center = T)
  #print(summary(mod))
  
  mod = pls::plsr(Nm_values ~ stat.sim.scaled, ncomp = PLS_use, center = T)
  
  #print(summary(mod))
  
  stat.sim.scaled_PLS = mod$scores
  
  ### real
  stat.obs.scaled_PLS = (stat.obs.scaled - colMeans(stat.sim.scaled)) %*%  mod$projection
  
  
  ### PLS abc
  out_PLS_full <- abc(target=stat.obs.scaled_PLS, param=Nm_values, sumstat=stat.sim.scaled_PLS, tol=tol, method = method, numnet = numnet, sizenet = sizenet)
  
  if (method == 'rejection'){
    out_PLS_full$adj.values = out_PLS_full$unadj.values
  }
  
  

  out_PLS_full = remove_outlier_Nm(out_PLS_full)
  
  return(out_PLS_full)
  
}



test_specific_Nm = function(fn, SS_sim, SS_target, n_comps, Nm_values, Nm_target, idx, method = 'neuralnet', tol = .2, sizenet = 4, numnet = 1){
  
  n_sims = dim(SS_sim)[1]
  n_real = diff(idx) + 1
  reps = length(SS_target)
  
  results = array(NA, dim = c(n_real,3,1, reps)) # n_sims, true/pred, N/m, reps
  dimnames(results) = list(paste0('sim', 1:n_real), c('True', 'Pred', 'Var'), c('N'))
  j=1
  for(i in idx[1]:idx[2]){
    print(i)
    for (r in 1:reps){
      print(r)
      
      out = fn(SS_sim, SS_target[[r]][i,], n_comps, Nm_values, method, tol, sizenet = sizenet, numnet = numnet)
      print('done')
      
      x = try(summary(out, print = F))
      if (class(x) == 'try-error') print(out$adj.values)
      
      results[j,1,,r] = Nm_target[i]
      results[j,2,,r] = summary(out, print = F)[4,1]
      results[j,3,,r] = colVars(out$adj.values)
    }
    j=j+1
    
  }
  
  return(results)
}


do_abc_trunc = function(SS_sim, SS_real, PLS_use, Nm_values,  method = 'neuralnet', tol = .1, numnet = 10, sizenet = 5){
  
  
  ### construct some posteriors:
  
  N_draws = dim(SS_sim)[1]
  
  # scale data
  scaled <-scale(rbind(SS_sim, SS_real) )
  stat.sim.scaled <- scaled[-(N_draws+1),]
  stat.obs.scaled <- scaled[N_draws+1,]
  
  
  # do PLS
  
  
  ### simulated
  mod = pls::plsr(stat.sim.scaled[,1:2] ~ stat.sim.scaled, ncomp = PLS_use, center = T)
  stat.sim.scaled_PLS = mod$scores
  
  ### real
  stat.obs.scaled_PLS = (stat.obs.scaled - colMeans(stat.sim.scaled)) %*%  mod$projection
  
  
  ### PLS abc
  out_PLS_full <- abc(target=stat.obs.scaled_PLS, param=Nm_values, sumstat=stat.sim.scaled_PLS, tol=tol, method = method, sizenet = sizenet, numnet = numnet)
  out_PLS = summary(out_PLS_full, print = F)
  
  
  N_min = out_PLS[1,1]
  N_max = out_PLS[7,1]
  m_min = out_PLS[1,2]
  m_max = out_PLS[7,2]
  
  filt = rep(FALSE, dim(Nm_values)[1])
  # truncation
  while (sum(filt) < 40){
    filt_N = Nm_values[,1] > N_min & Nm_values[,1] < N_max
    filt_m = Nm_values[,2] > m_min & Nm_values[,2] < m_max
    filt = filt_N & filt_m
    
    N_min = .9*N_min-1; N_max = N_max * 1.1+1
    m_min = .9*m_min-0.01; m_max = m_max * 1.1+0.01
    
  }
  
  
  stat.sim.scaled_trunc = stat.sim.scaled[filt,]
  Nm_values_trunc = Nm_values[filt,]
  
  
  out_linreg = abctools::saABC(theta = Nm_values_trunc[,1:2], stat.sim.scaled_trunc, plot = F)
  
  SS_sim_linreg = stat.sim.scaled_trunc %*% t(out_linreg$B)
  SS_real_linreg = stat.obs.scaled %*% t(out_linreg$B)
  
  
  out_trunc_full <- abc(target=SS_real_linreg, param=Nm_values_trunc, sumstat=SS_sim_linreg, tol=tol, method = method, sizenet = sizenet, numnet = numnet)
  
  if (method == 'rejection'){
    out_trunc_full$adj.values = out_trunc_full$unadj.values
  }
  
  
  out_trunc_full = remove_outlier(out_trunc_full)
  
  
  return(out_trunc_full)
}




do_abc_PCA = function(SS_sim, SS_real, PLS_use, Nm_values, method = 'neuralnet', tol = .1, numnet = 10, sizenet = 5){
  
  ### construct some posteriors:
  
  N_draws = dim(SS_sim)[1]
  
  # scale data
  scaled <-scale(rbind(SS_sim, SS_real) )
  stat.sim.scaled <- scaled[-(N_draws+1),]
  stat.obs.scaled <- scaled[N_draws+1,]
  
  
  # do PCA
  
  pc <- prcomp(stat.sim.scaled, scale = F, center= F)
  pc.sim <- (pc$x[,1:PLS_use])
  
  stat.obs.pc <- c()
  for (p in 1:PLS_use){
    stat.obs.pc <- c(stat.obs.pc, sum(stat.obs.scaled * pc$rotation[,p]))
  }
  
  out_PCA_full <- abc(target=stat.obs.pc, param=Nm_values, sumstat=pc.sim, tol=tol, method = method,numnet = numnet, sizenet = sizenet)
  ### PLS abc

  if (method == 'rejection'){
    out_PCA_full$adj.values = out_PCA_full$unadj.values
  }
  
  
  out_PCA_full = remove_outlier(out_PCA_full)
  
  
  return(out_PCA_full)
  
}

remove_outlier = function(ABC_o, N_range = c(1e3,1e6), m_range = c(0,0.8)){
  filt_N = in_range(ABC_o$adj.values[,1], N_range)
  filt_m = in_range(ABC_o$adj.values[,2], m_range)
  
  m_low = ABC_o$adj.values[,2] < m_range[1]
  m_high = ABC_o$adj.values[,2] > m_range[2]
  
  N_low = ABC_o$adj.values[,1] < N_range[1]
  N_high = ABC_o$adj.values[,1] > N_range[2]
  
  
  ABC_o$adj.values[m_low,2] = m_range[1]
  ABC_o$adj.values[m_high,2] = m_range[2]
  
  
  ABC_o$adj.values[N_low,1] = N_range[1]
  ABC_o$adj.values[N_high,1] = N_range[2]
  
  
  return(ABC_o)
}

remove_outlier = function(ABC_o, N_range = c(200,25000), m_range = c(0,0.8)){
  filt_N = in_range(ABC_o$adj.values[,1], N_range)
  filt_m = in_range(ABC_o$adj.values[,2], m_range)
  
  m_low = ABC_o$adj.values[,2] < m_range[1]
  m_high = ABC_o$adj.values[,2] > m_range[2]
  
  N_low = ABC_o$adj.values[,1] < N_range[1]
  N_high = ABC_o$adj.values[,1] > N_range[2]

  
  ABC_o$adj.values[m_low,2] = m_range[1]
  ABC_o$adj.values[m_high,2] = m_range[2]
  
  
  ABC_o$adj.values[N_low,1] = N_range[1]
  ABC_o$adj.values[N_high,1] = N_range[2]
  
  
  return(ABC_o)
}




  

remove_outlier_Nm = function(ABC_o, Nm_range = c(0,25000)){
  filt_Nm = in_range(ABC_o$adj.values[,1], Nm_range)

  low = ABC_o$adj.values[,1] < Nm_range[1]
  high = ABC_o$adj.values[,1] > Nm_range[2]
  

  
  
  ABC_o$adj.values[low,1] = Nm_range[1]
  ABC_o$adj.values[high,1] = Nm_range[2]
  

  
  
  return(ABC_o)
}



PLS_examine = function(SS_sim,Nm, PLS_use = 8 ){
  
  ### construct some posteriors:
  
  N_draws = dim(SS_sim)[1]
  
  # scale data
  scaled <-scale(rbind(SS_sim, SS_real) )
  stat.sim.scaled <- scaled[-(N_draws+1),]
  stat.obs.scaled <- scaled[N_draws+1,]
  
  
  # do PLS
  
  ### simulated
  print(dim(Nm))
  print(dim(stat.sim.scaled))
  mod = pls::plsr(Nm ~ stat.sim.scaled, ncomp = PLS_use, center = T)
  return(mod)
  
}



do_abc_PLS_inv = function(SS_sim, SS_real, PLS_use, Nm_values, method = 'neuralnet', tol = .1, numnet = 10, sizenet = 5){
  
  ### construct some posteriors:
  
  
  Nm_values[,1] = 1/Nm_values[,1]
  
  N_draws = dim(SS_sim)[1]
  
  # scale data
  scaled <-scale(rbind(SS_sim, SS_real) )
  stat.sim.scaled <- scaled[-(N_draws+1),]
  stat.obs.scaled <- scaled[N_draws+1,]
  
  
  # do PLS
  
  ### simulated
  mod = pls::plsr(stat.sim.scaled[,1:2] ~ stat.sim.scaled, ncomp = PLS_use, center = T)
  stat.sim.scaled_PLS = mod$scores
  
  ### real
  stat.obs.scaled_PLS = (stat.obs.scaled - colMeans(stat.sim.scaled)) %*%  mod$projection
  
  
  ### PLS abc
  out_PLS_full <- abc(target=stat.obs.scaled_PLS, param=Nm_values, sumstat=stat.sim.scaled_PLS, tol=tol, method = method, numnet = numnet, sizenet = sizenet)
  
  return(out_PLS_full)
  
}



cross_validate = function(fn, SS_sim, n_cross_or_selected, n_comps, Nm_values, method = 'neuralnet', tol = .2, numnet = 10, sizenet = 5){
  
  n_sims = dim(SS_sim)[1]
  
  if (length(n_cross_or_selected) == 1){
    selected = sample(1:n_sims, n_cross_or_selected, replace = F)
  } else{
    selected = n_cross_or_selected
  }
  
  results = array(NA, dim = c(length(selected),2,2)) # n_sims, true/pred, N/m
  dimnames(results) = list(paste0('sim', selected), c('True', 'Pred'), c('N', 'm'))
  
  for(i in 1:length(selected)){
    print(i)
    sel = selected[i]
    out = fn(SS_sim[-sel,], SS_sim[sel,], n_comps, Nm_values[-sel,], method, tol, numnet = numnet, sizenet = sizenet)
    results[i,1,] = Nm_values[sel,1:2]
    results[i,2,] = summary(out, print = F)[4,1:2]
  }
  
  return(results)
}


in_range = function(int, ran){
  ran = sort(ran)
  return(int > ran[1] & int < ran[2])
}

in_range(2.8, c(3,6))

check_CI = function(fn, SS_sim, n_cross_or_selected, n_comps, Nm_values, method = 'neuralnet', tol = .2, sizenet = sizenet, numnet = numnet){
  
  n_sims = dim(SS_sim)[1]
  
  if (length(n_cross_or_selected) == 1){
    selected = sample(1:n_sims, n_cross_or_selected, replace = F)
  } else{
    selected = n_cross_or_selected
  }
  
  results = array(NA, dim = c(length(selected),2,2)) # n_sims, true/pred, N/m
  dimnames(results) = list(paste0('sim', selected), c('True', 'Contained'), c('N', 'm'))
  
  for(i in 1:length(selected)){
    print(i)
    sel = selected[i]
    
    out = fn(SS_sim[-sel,], SS_sim[sel,], n_comps, Nm_values[-sel,], method, tol=tol,sizenet = sizenet, numnet = numnet)
    so = summary(out, print = F)
    
    N_range = so[c(2,6),1]; N_true = Nm_values[sel,1]
    m_range = so[c(2,6),2]; m_true = Nm_values[sel,2]
    
    

    results[i,1,] = c(N_true, m_true)
    results[i,2,] = c(in_range(N_true, N_range), in_range(m_true, m_range))
  }
  
  return(results)
}


check_CI_width = function(fn, SS_sim, n_cross_or_selected, n_comps, Nm_values, method = 'neuralnet', tol = .2, sizenet=sizenet, numnet = numnet){
  
  n_sims = dim(SS_sim)[1]
  
  if (length(n_cross_or_selected) == 1){
    selected = sample(1:n_sims, n_cross_or_selected, replace = F)
  } else{
    selected = n_cross_or_selected
  }
  
  results = array(NA, dim = c(length(selected),3,2)) # n_sims, true/pred, N/m
  dimnames(results) = list(paste0('sim', selected), c('True', 'Min', 'Max'), c('N', 'm'))
  
  for(i in 1:length(selected)){
    print(i)
    sel = selected[i]
    
    out = fn(SS_sim[-sel,], SS_sim[sel,], n_comps, Nm_values[-sel,], method, tol = tol, sizenet = sizenet, numnet = numnet)
    so = summary(out, print = F)
    
    N_range = so[c(2,6),1]; N_true = Nm_values[sel,1]
    m_range = so[c(2,6),2]; m_true = Nm_values[sel,2]
    
    
    
    results[i,1,] = c(N_true, m_true)
    results[i,2:3,1] = N_range
    results[i,2:3,2] = m_range
  }
  
  return(results)
}



test_specific = function(fn, SS_sim, SS_target, n_comps, Nm_values, Nm_target, idx, method = 'neuralnet', tol = .2, sizenet = 4, numnet = 1){
  
  n_sims = dim(SS_sim)[1]
  n_real = diff(idx) + 1
  reps = length(SS_target)
  
  results = array(NA, dim = c(n_real,3,2, reps)) # n_sims, true/pred, N/m, reps
  dimnames(results) = list(paste0('sim', 1:n_real), c('True', 'Pred', 'Var'), c('N', 'm'))
  j=1
  for(i in idx[1]:idx[2]){
    print(i)
    for (r in 1:reps){
      print(r)
      
      out = fn(SS_sim, SS_target[[r]][i,], n_comps, Nm_values, method, tol, sizenet = sizenet, numnet = numnet)
      print('done')
      
      x = try(summary(out, print = F))
      if (class(x) == 'try-error') print(out$adj.values)
      
      results[j,1,,r] = Nm_target[i,1:2]
      results[j,2,,r] = summary(out, print = F)[4,1:2]
      results[j,3,,r] = colVars(out$adj.values)
    }
    j=j+1
    
  }
  
  return(results)
}


var_m = 0.05752135 #var(m_prior_draw(1e7))
var_N = 51233441 #var(N_prior_draw(1e7))



loocv_err = function(results){
  N_err = (results[,1,1] - results[,2,1])**2 / var_N
  m_err = (results[,1,2] - results[,2,2])**2 / var_m
  
  
  
  return(cbind(N_err, m_err))
}

MAE_err = function(results){
  N_err = abs(results[,1,1] - results[,2,1]) / results[,1,1]
  m_err = abs(results[,1,2] - results[,2,2]) / results[,1,2]
  
  return(cbind(N_err, m_err))
}

RMS_err = function(results){
  N_err = abs(results[,1,1] - results[,2,1])
  m_err = abs(results[,1,2] - results[,2,2])
  
  return(cbind(N_err, m_err))
}


RMS_err = function(results){
  N_err = abs(log2((results[,1,1]/ results[,2,1])))
  m_mod = results[,2,2]
  m_mod[m_mod == 0] = 0.01
  m_err = abs(log2((results[,1,2]/ m_mod)))
  
  return(cbind(N_err, m_err))
}

RMS_err = function(results){
  N_err = abs(log2((results[,1,1]/ results[,2,1])))
  m_err = abs(log2((results[,1,2]/ results[,2,2])))
  Nc = cor(results[,1,1],results[,2,1])
  mc = cor(results[,1,2],results[,2,2])
  return(cbind(Nc,mc))
}


RMS_err = function(results){
  N_err = abs(results[,1,1] - results[,2,1])/results[,1,1]
  m_err = abs(results[,1,2] - results[,2,2])/results[,1,2]
  
  
  return(cbind(N_err, m_err))
}

RMS_err = function(results){
  N_err = abs(results[,1,1] - results[,2,1])/results[,1,1]
  m_err = abs(results[,1,2] - results[,2,2])/results[,1,2]
  
  
  n = dim(results)[1]
  N_draw = matrix(rep(N_prior_draw(100),n), nr = n, byrow = T)
  m_draw = matrix(rep(m_prior_draw_no0(100),n), nr = n, byrow = T)
  

  
  N_prior_err = rowMeans(abs(results[,1,1] - N_draw)/results[,1,1])
  m_prior_err = rowMeans(abs(results[,1,2] - m_draw)/results[,1,2])
  

  return(cbind(N_err, m_err))
}

m_prior_draw_no0 = function(n){
  v = c()
  while (length(v) < n){
    temp = m_prior_draw(1)
    if (temp != 0){
      v = c (v, temp)
    }
  }
  return(v)
}


get_sd_pool = function(n_pool, cov_pool, fr = 0.01){
  return(sqrt(0.5*fr*(1-fr) / (n_pool) * (1 + (2*n_pool-1)/cov_pool)))
}

get_sd_dip = function(n_dip, cov_dip, fr = 0.01){
  return(sqrt(fr*(1-fr) * (1+cov_dip) / (2*cov_dip*n_dip)))
}

# function to find how many diploid individuals would need to be sequenced to get same sd as pooled sequencing
# params are 1) size of pool 2) total pool coverage 3) individual seq coverage
get_diploid_equivalent = function(n_pool, cov_pool, cov_dip, fr = 0.01){
  sd_pool = get_sd_pool(n_pool, cov_pool, fr)
  curried_sd_dip = function(n_pool){
    return(get_sd_dip(n_pool, cov_dip, fr))
  }
  return( uniroot(f = function(n_dip){curried_sd_dip(n_dip) - sd_pool}, interval = c(0,1e5)  )$root )
}
get_diploid_equivalent_vec = Vectorize(get_diploid_equivalent)


# test single sample varying sample sizes
## across N, m



trim_mean = function(v, n_remove = 1){

  sorted_v = sort(v, decreasing = T)
  if (n_remove){
    acc = (n_remove+1):(length(v)-n_remove)
    v_trim = sorted_v[acc]
    return(mean(v_trim))
  }
    
  return(mean(v))
}




sur = function(v){
  return(sort(unique(round(v))))
}
get_S_vec = function(base, muls, AF = T, add_Inf = T){
  S_vec = c()
  S_vec_list = list()
  i=1
  for (mul in muls){
    S_vec = c(S_vec, sur(base*mul))
    S_vec_list[[i]] = sur(base*mul)
    i=i+1
  }
  
  if (AF){
    for (mul in muls){
      S_vec = c(S_vec, sur(get_diploid_equivalent_vec(1e100, sur(base*mul)*30, 30) ))
      S_vec_list[[i]] = sur(get_diploid_equivalent_vec(1e100, sur(base*mul)*30, 30) )
      i=i+1
    }
  }
  if (add_Inf){
    S_vec = c(S_vec, Inf)
    S_vec_list[i] = c(Inf)
  }
  return(list(sur(S_vec), S_vec_list))
}

select_bounded = function(n, Nm_values, N_range = c(-Inf, Inf), m_range = c(-Inf, Inf)){
  filt_N = in_range(Nm_values[,1], N_range)
  filt_m = in_range(Nm_values[,2], m_range)
  print(Nm_values[filt_N & filt_m,])
  print(which(filt_N & filt_m))
  return(sample(which(filt_N & filt_m), n ))
}
mean_new = function(d, ind){
  return(mean(d[ind]))
}
do_boot = function(d, n = 1e3){
  x = boot(d,mean_new, n)
  return(boot.ci(x)$normal[2:3])
}


combine = function(M1, M2, axis, M1_index, M2_slice){
  M1_pre = asub(M1,idx = 1:M1_index, dims = axis, drop = F)
  M1_post = asub(M1,idx = (M1_index+1):dim(M1)[axis], dims = axis, drop = F)
  M2_slice = asub(M2, idx = M2_slice, dims = axis, drop = F)
  return( abind(M1_pre, M2_slice, M1_post, along= axis))
}

compare_CV = function(LD, AF, LD_ci, AF_ci, ymax = 1){
  #par(mfrow=c(1,1), mai = c(1,1,0.2,0.2))
  
  plot(NA, yaxs = 'i' ,  bty = 'l', type = 'p', pch = 21, bg = 'white', cex = .8, ylab = 'Relative Error', xlab = 'Number of Sampled Subpopulations', xaxt = 'n', ylim = c(min(0,LD_ci[,,], AF_ci[,,]) ,max(LD_ci[,,]/100, AF_ci[,,]/100,ymax )) , xlim = c(1,9)) ; axis(1, at = c(1,2,3,5,9), labels = c(1,2,3,5,9))
  #plot(NA, yaxs = 'i' , frame = F, type = 'p', pch = 21, bg = 'white', cex = .8, ylab = 'Relative Error', xlab = 'Number of Sampled Subpopulations', xaxt = 'n', ylim = c(min(0,LD_ci[,,], AF_ci[,,]) ,1) , xlim = c(1,9)) ; axis(1, at = c(1,2,3,5,9), labels = c(1,2,3,5,9))
  
  for (i in 1:5){lines(rep(c(1,2,3,5,9)[i],2),LD_ci[i,,1])}
  lines(c(1,2,3,5,9), LD[,1], lw = .5)
  points(c(1,2,3,5,9), LD[,1], pch = 21, bg = 'white', cex = 1)
  
  
  for (i in 1:5){lines(rep(c(1,2,3,5,9)[i],2),AF_ci[i,,1], col = 'red')}
  lines(c(1,2,3,5,9), AF[,1], lw = .5, col  = 'red')
  points(c(1,2,3,5,9), AF[,1],pch = 21, bg = 'white', cex = 1, col = 'red')
  
  for (i in 1:5){lines(rep(c(1,2,3,5,9)[i],2),LD_ci[i,,2])}
  lines(c(1,2,3,5,9), LD[,2], lw = .5)
  points(c(1,2,3,5,9), LD[,2],frame = F, type = 'p', pch = 25, bg = 'white', cex = 1)
  
  for (i in 1:5){lines(rep(c(1,2,3,5,9)[i],2),AF_ci[i,,2], col = 'red')}
  lines(c(1,2,3,5,9), AF[,2], lw = .5, col  = 'red')
  points(c(1,2,3,5,9), AF[,2],pch = 25, bg = 'white', cex = 1, col = 'red')
  
  #legend('topleft',cex = .9,inset = 0.02, legend = c('N - LD', 'N - AF', 'm - LD', 'm - AF' ), pch = c(1,1,6, 6), col = c('black', 'red', 'black', 'red'))
  
}

compare_CV_tp = function(LD, AF, LD_ci, AF_ci){
  #par(mfrow=c(1,1), mai = c(1,1,0.2,0.2))
  
  plot(NA, yaxs = 'i' , bty = 'l' , type = 'p', pch = 21, bg = 'white', cex = .8, ylab = 'Relative Error', xlab = 'Number of Sampled Subpopulations', xaxt = 'n', ylim = c(min(0,LD_ci[,,], AF_ci[,,]) ,max(LD_ci[,,], AF_ci[,,] )) , xlim = c(1,4)) ; axis(1, at = c(1,2,3,4), labels = c(1,2,3,4))
  #plot(NA, yaxs = 'i' , frame = F, type = 'p', pch = 21, bg = 'white', cex = .8, ylab = 'Relative Error', xlab = 'Number of Sampled Subpopulations', xaxt = 'n', ylim = c(min(0,LD_ci[,,], AF_ci[,,]) ,1) , xlim = c(1,9)) ; axis(1, at = c(1,2,3,5,9), labels = c(1,2,3,5,9))
  
  for (i in 1:4){lines(rep(c(1,2,3,4)[i],2),LD_ci[i,,1])}
  lines(c(1,2,3,4), LD[,1], lw = .5)
  points(c(1,2,3,4), LD[,1], pch = 21, bg = 'white', cex = 1)
  
  
  for (i in 1:3){lines(rep(c(2,3,4)[i],2),AF_ci[i,,1], col = 'red')}
  lines(c(2,3,4), AF[,1], lw = .5, col  = 'red')
  points(c(2,3,4), AF[,1],pch = 21, bg = 'white', cex = 1, col = 'red')
  
  for (i in 1:4){lines(rep(c(1,2,3,4)[i],2),LD_ci[i,,2])}
  lines(c(1,2,3,4), LD[,2], lw = .5)
  points(c(1,2,3,4), LD[,2],frame = F, type = 'p', pch = 25, bg = 'white', cex = 1)
  
  for (i in 1:3){lines(rep(c(2,3,4)[i],2),AF_ci[i,,2], col = 'red')}
  lines(c(2,3,4), AF[,2], lw = .5, col  = 'red')
  points(c(2,3,4), AF[,2],pch = 25, bg = 'white', cex = 1, col = 'red')
  
  #legend('topleft',cex = .9,inset = 0.02, legend = c('N - LD', 'N - AF', 'm - LD', 'm - AF' ), pch = c(1,1,6, 6), col = c('black', 'red', 'black', 'red'))
  
}



frac_log = function(x,y,frac){
  lx = log2(x); ly = log2(y)
  exp_val = frac*ly + (1-frac)*lx
  return(2^exp_val)
}
frac_log = Vectorize(frac_log)


fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    print(ds)
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    #text(x[1], y[1], text, cex=cex, ...)
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      print(text)
      print(y)
      print(fig)
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  print(y)
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}


fig_label_logy <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    print(ds)
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    #text(x[1], y[1], text, cex=cex, ...)
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]

      y <- frac_log(y[1], y[2], fig[3:4])
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 0.4/100
  
  print(y)
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}


fig_label_logx <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    

    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig");
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- frac_log(x[1], x[2], fig[1:2])
      y <- y[1] + dy * fig[3:4]
      print(x)
    }
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 0.4/100
  sh <- strheight(text, cex=cex) * 60/100 
  

  
  x1 <- switch(pos,
               topleft     =log2(2^(x[1]) + sw), 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  #x1 = x1 + abs(0.4*x1
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}




fig_label_logxy <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig");
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- frac_log(x[1], x[2], fig[1:2])
      y <- frac_log(y[1], y[2], fig[3:4])
      print(x)
    }
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 0.4/100
  sh <- strheight(text, cex=cex) * 0.4/100 
  
  
  
  x1 <- switch(pos,
               topleft     =log2(2^(x[1]) + sw), 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  #x1 = x1 + abs(0.4*x1
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}





just_guess = function(Nm_target, idxs){
  n_target = diff(idxs)+1
  reps = 1000
  results = array(NA, dim = c(n_target,3,2, reps)) # n_sims, true/pred, N/m, reps
  dimnames(results) = list(paste0('sim', 1:n_target), c('True', 'Pred', 'Var'), c('N', 'm'))
  
  results[,1,,] = Nm_target[idxs[1]:idxs[2],]
  
  guess_m = m_prior_draw_no0(reps)
  guess_N = N_prior_draw(reps)
  
  results[,2,1,] = matrix(rep(guess_N,n_target),nc=reps, byrow = T)
  results[,2,2,] = matrix(rep(guess_m,n_target),nc=reps, byrow = T)
  
  results[,3,1,] = var(guess_N)
  results[,3,2,] = var(guess_m)
  
  return(results)
}
