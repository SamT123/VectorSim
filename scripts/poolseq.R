fr = 0.01

cov_dip = 3
n_dip = 10000

cov_pool = cov_dip*n_dip
n_pool = 5000

sd_dip = sqrt(fr*(1-fr) * (1+cov_dip) / (2*cov_dip*n_dip))
sd_dip

sd_pool= sqrt(0.5*fr*(1-fr) / (n_pool) * (1 + (2*n_pool-1)/cov_pool))
sd_pool
sd_pool / sd_dip


get_sd_pool = function(n_pool, cov_pool, fr = 0.01){
  return(sqrt(0.5*fr*(1-fr) / (n_pool) * (1 + (2*n_pool-1)/cov_pool)))
}

get_sd_dip = function(n_dip, cov_dip, fr = 0.01){
  return(sqrt(fr*(1-fr) * (1+cov_dip) / (2*cov_dip*n_dip)))
}


get_sd_dip(1000,30)
get_sd_dip(1000,1e100)

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






# function to find how many diploid individuals would need to be sequenced to get same sd as pooled sequencing
# params are 1) size of pool 2) total pool coverage 3) individual seq coverage
get_inf_equivalent = function(n_dip, cov_dip, fr = 0.01){
  sd_finite = get_sd_dip(n_dip, cov_dip, fr)
  curried_sd_dip = function(n_dip){
    return(get_sd_dip(n_dip, 1e100, fr))
  }
  return( uniroot(f = function(n_dip){curried_sd_dip(n_dip) - sd_finite}, interval = c(0,1e5)  )$root )
}
get_inf_equivalent_vec = Vectorize(get_diploid_equivalent)


get_inf_equivalent(100,30)
# how to use:

## I have a sample size S from full diploid sampling, and I want to know whether that amount of sequencing resources would ave been better used with poolseq
## so eg diploid sample size 500, coverage 30 - what is the equivalent for poolseq
S_base = 1/(  matrix(c(1,2,3,5,9)) %*% c(1,2,3,4,5) )
S_base = sort(unique(c(S_base)))


S_vec = S_base * 200; S_vec

S_vec_AF = c()
for (S in S_vec){
  S_vec_AF = c(S_vec_AF, get_diploid_equivalent(S*, S*20, 20))
}

get_diploid_equivalent_vec(1e3, S_vec*20, 20)/S_vec

pool_size = seq(10,5000,10)
seq_resource = 55 * 30
effective_size  = get_diploid_equivalent_vec(pool_size, seq_resource, 30)

plot(pool_size, effective_size, type = 'l', ylab = 'Equivalent size if using individual sequencing', asp = 1, ylim = c(0,5000), xlim = c(0,5000));axis(1);axis(2); abline(1,1, lty = 2, col = 'red')

get_diploid_equivalent_vec(10000, seq_resource, 20) / get_diploid_equivalent_vec(1e100, seq_resource, 20)


round(get_diploid_equivalent_vec(1e16, S_vec[-38]*30, 30)) == round(get_diploid_equivalent_vec(1e14, S_vec[-38]*30, 30))


round(get_diploid_equivalent_vec(1e100, S_vec[-38]*30, 30)) == round(get_diploid_equivalent_vec(1e30, S_vec[-38]*30, 30))


get_diploid_equivalent(1000, 500/2/3*30, 30)

get_sd_pool(n_pool = 1e3, cov_pool = 500/3/3*30)/get_sd_pool(n_pool = 1e100, cov_pool = 500/3/3*30)
get_sd_pool(n_pool = 1e3, cov_pool = 500/3/3*30)/get_sd_dip(n_dip =  500/3/3*30, cov_dip = 30)
S_vec_AF
S_vec

S_vec_AF / S_vec

get_diploid_equivalent(1e100, 500/6*30, 30)


get_diploid_equivalent(1e100, 500*30, 30)

plot(pool_size, get_diploid_equivalent_vec(pool_size, seq_resource, 30) / get_diploid_equivalent_vec(1e100, seq_resource, 30))




plot(pool_size, get_sd_pool(pool_size, seq_resource)/ get_sd_pool(1e100, seq_resource), ylim = c(0,10) ); abline(h = get_sd_dip(seq_resource/30, 30)/ get_sd_pool(1e100, seq_resource))



seq_resource = 83 * 30

get_sd_dip(seq_resource/30, 30)
get_sd_pool(500, seq_resource)
get_sd_pool(1e100, seq_resource)

pop = initialise_pop(pairs_per_recom_rate = 1000)

pop = propegate_population(pop, N, m, generations)

sample1 = take_sample(pop, size)

r2_1 = calculate_r2(sample1)

pop = propegate_population(pop, N, m, generations)

sample2 = take_sample(pop, size)

r2_2 = calculate_r2(sample2)

Fc = calculate_Fc(sample1, sample2, temporal_separation)



