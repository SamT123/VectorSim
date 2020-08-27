library(abc)

# generate "true data" from large stepping stone model
 
## choose  (N, m) to estimate
N = 1000
m = 0.1

n_subpops = 11**2
S = Inf

## Summary statstics to gather
sampling_intervals = c(2,4,6,8,10,20)
recom_rates = c(0.05, 0.2, 0.35, 0.5)


## simulation parameters
loci = 1e2
burn_in = 300 

reset_period = 10
sampling_repeats = 20

## some derived things
sampling_period = reset_period + max(t_intervals)
t_intervals = c(0,sampling_intervals)
ms = c(m)
M_fn = construct_M_4_step
p0 = rep(.25, 4)
N_v = rep(N, n_subpops)



o <- burn_in_checker(ms, M_fn, cs, burn_in, N_v, loci, p0)

SS_real <- get_SS_df(N_v = N_v,
                                  ms = ms,
                                  M_fn = M_fn,
                                  p0 = p0,
                                  loci = loci,
                                  burn_in = burn_in,
                                  t_intervals = t_intervals,
                                  rest_period = rest_period,
                                  sampling_repeats = sampling_repeats,
                                  sample_subpops = n_subpops/2+.5,
                                  S = S,
                                  cs = recom_rates)



SS_real_single = subset(SS_real, i==0)






library(coda)

# assuming the data are 10 samples of a normal distribution
# with mean 5.3 and sd 2.7
real_data <- SS_real
  
  # we want to use ABC to infer the parameters that were used. 
  # we sample from the same model and use mean and variance
  # as summary statstitics. We return true for ABC acceptance when
  # the difference to the data is smaller than a certain threshold
  
  meandata <- mean(data)
standarddeviationdata <- sd(data)

ABC_acceptance <- function(par){
  
  # prior to avoid negative standard deviation
  if (par[2] <= 0) return(F) 
  
  # stochastic model generates a sample for given par
  ## this will be a vector of r2 values and Fc values
  samples <- rnorm(10, mean =par[1], sd = par[2])
  
  # comparison with the observed summary statistics
  diffmean <- abs(mean(samples) - meandata)
  diffsd <- abs(sd(samples) - standarddeviationdata)
  if((diffmean < 0.1) & (diffsd < 0.2)) return(T) else return(F)
}


# we plug this in in a standard metropolis hastings MCMC, 
# with the metropolis acceptance exchanged for the ABC acceptance

run_MCMC_ABC <- function(startvalue, iterations){
  
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  
  for (i in 1:iterations){
    
    # proposalfunction for (N,m)
    proposal = rnorm(2,mean = chain[i,], sd= c(0.7,0.7))
    
    if(ABC_acceptance(proposal)){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(mcmc(chain))
}

posterior <- run_MCMC_ABC(c(4,2.3),300000)
plot(posterior)



















###############
## create a large dataset of N, m, AFx,y,z, LDx,y,z

n_subpops = 11**2
S = Inf
## Summary statstics to gather
sampling_intervals = c(2,4,6,8,10,20)
recom_rates = c(0.05, 0.2, 0.35, 0.5)

## simulation parameters
loci = 2
burn_in = 300 
reset_period = 10
sampling_repeats = 20

# sample_subpops
extension = 1
centre = n_subpops/2+.5
centre_line = c((centre-extension):(centre+extension))
sample_subpops = c(centre_line, centre_line+sqrt(n_subpops), centre_line-sqrt(n_subpops))
## some derived things
sampling_period = reset_period + max(t_intervals)
t_intervals = c(0,sampling_intervals)
ms = c(m)
M_fn = construct_M_4_step
p0 = rep(.25, 4)
N_v = rep(N, n_subpops)


m_prior_draw = function(n=1){
  beta_out = rbeta(n, 1, 3)
  adj = beta_out
  adj[adj < 0.01] = 0
  return(adj)
}

N_prior_draw = function(n=1){
  return(runif(n, 2000,20000))
}

n_draws = 100
i=1
N_sim = N_prior_draw(n_draws)
SS_sim = data.frame()
for (N in N_sim){
  m_sim = m_prior_draw(n_draws)
  print(i)
  i=i+1
  N_v = rep(N, n_subpops)
  SS_sim_Ni <- get_SS_df(N_v = N_v,
                    ms = m_sim,
                    M_fn = M_fn,
                    p0 = p0,
                    loci = loci,
                    burn_in = burn_in,
                    t_intervals = t_intervals,
                    rest_period = rest_period,
                    sampling_repeats = sampling_repeats,
                    sample_subpops = ,
                    S = S,
                    cs = recom_rates)
  SS_sim_Ni$N = N
  SS_sim = rbind(SS_sim, SS_sim_Ni)
}




  runif(0,)
hist(runif(1e5, 2000,20000), breaks = 100)



# estimate


# calculate PCs

# estimate with PCs

