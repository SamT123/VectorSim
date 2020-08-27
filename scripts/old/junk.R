


pop <- initialise_population(n_subpops = ns, loci = 200, t = 1e3)

N = 5000
m = 0.1

M <- construct_M_4_step(m, ns)
N_v = rep(N, ns)


fst_v <- get_Fst_vec(pop$pA_post, N_v = rep(1000,20), 1000)
r2 <- apply(pop$r_post, MARGIN = c(2), mean)

plot(1:(1000), (ma(fst_v, n=1)), type = "l" ); abline(v=c(eqm1, start2, eqm2), col = 'red')
plot(1:(1001), ma(r2, n=1), type = "l" ); abline(v=c(eqm1, start2, eqm2), col = 'red')

#### ABC MCMC

# have target SS


# initialise population
init_pop <- initialise_population(n_subpops = ns, loci = 200, t = 1e)

# choose N,m
N = 5000
m = 0.1

M <- construct_M_waples(m, ns)
N_v = rep(N, ns)


t_intervals = c(0,2,3,4,5,7,9,15)
reps = 10
gap = 5
one_period = max(t_intervals) + gap
total_period = one_period*reps
c=0.5

sample_subpops = 1:ns

# generate SS

o <- get_sample(pop,N_v, m, construct_M_waples, c, t_intervals, reps, gap, 1:10, min_loci = 100)

SS = o$SS
pop = o$pop
# return SS, population 

# accept/reject N,m

# choose N,m




#### ABC plain ####
ns = 5
m_prior_draw = function(n=1){
  beta_out = rbeta(n, 1, 3)
  adj = beta_out
  adj[adj < 0.01] = 0
  return(adj)
}

N_prior_draw = function(n=1){
  return(runif(n, 2000,20000))
}

i=1
n_draws = 3
N_sim = N_prior_draw(n_draws)
SS_sim = array(NA, dim = c(0, 8))
for (N in N_sim){
  cat(i)
  cat("\n")
  cat(N)
  cat("\n")
  i = i+1
  m_sim = m_prior_draw(n_draws)
  cat("Initialising Population\n")
  çc <- initialise_population(n_subpops = ns, loci = 200, t = 1e4)
  N_v = rep(N, ns)
  
  j=1
  for (m in rev(sort(m_sim))){
    cat(j)
    cat("\n")
    cat(m)
    cat("\n")
    
    o <- get_sample(pop, N_v, m, construct_M_waples, c, t_intervals, reps, gap, 1:ns, min_loci = 100)
    pop = o$pop
    SS = o$SS
    j = j+1
    SS_sim = rbind(SS_sim, SS)
  }
  
  
  
  
  
  
}

















o = get_popn_to_eqm(init_pop, M1, N_v, .5, 1)
pop <- o$popn
eqm1 <- o$eqm
pop = propegate_population(pop, N_v = N_v, M = M1, c=0.5, t_add = buffer)

start2 = pop$t_curr
o = get_popn_to_eqm(pop, M2, N_v, .5, pop$t_curr)
pop <- o$popn
eqm2 <- o$eqm
pop = propegate_population(pop, N_v = N_v, M = M2, c=0.5, t_add = buffer)








ns=10

## 
M1 = construct_M_waples(0.3, ns)
M2 = construct_M_waples(0.25, ns)
N_v = rep(1000, ns)


## samples needed
t_intervals = c(0,2,3,4,5,7,9,15)
reps = 10
gap = 5
one_period = max(t_intervals) + gap
total_period = one_period*reps
c=0.5

sample_subpops = 1:ns


init_pop <- initialise_population(n_subpops = ns, loci = 200, t = 1000)
# get to eqm
o = get_popn_to_eqm(init_pop, M1, N_v, c, 1)
pop <- o$popn
eqm <- o$eqm

# propegate
pop = propegate_population(pop, N_v = N_v, M = M1, c=c, t_add = total_period)

estimates = data.frame( 'subpop' = c(NA), 'method' = c(NA), 'sample_rep' = c(NA), 'est' = c(NA))


for (sp in sample_subpops){
  for (sample_rep in 0:(reps-1)){
    base = eqm + sample_rep * one_period
    sample_gens = base + t_intervals
    
    for (gen in sample_gens){
      s <- take_sample_r(pop$pAB_post[,gen,sp], pop$pAb_post[,gen,sp], pop$paB_post[,gen,sp], pop$pab_post[,gen,sp], S)
      estimates = rbind(estimates, c(sp,paste0("LD",c), sample_rep, get_r2_drift(s$pAB, s$pAb, s$paB, s$pab, S = S, c = c)))
    }
    
    s1 <- take_sample_r(pop$pAB_post[,sample_gens[1],sp], pop$pAb_post[,sample_gens[1],sp], pop$paB_post[,sample_gens[1],sp], pop$pab_post[,sample_gens[1],sp], S)
    for (gen in sample_gens[-1]){
      interval = gen - sample_gens[1]
      method_name = paste0("AF", interval)
      s2 <- take_sample_r(pop$pAB_post[,gen,sp], pop$pAb_post[,gen,sp], pop$paB_post[,gen,sp], pop$pab_post[,gen,sp], S)
      
      estimates =  rbind(estimates, c(sp, method_name, sample_rep, get_mean_fc(s1$pAB+s1$pAb, s2$pAB+s2$pAb, S, S, interval)))
      estimates =  rbind(estimates, c(sp, method_name, sample_rep, get_mean_fc(s1$pAB+s1$pAb, s2$pAB+s2$pAb, S, S, interval)))
      
    }
  }
}

estimates = estimates[-1,]

methods = unique(estimates$method)
mean_ests = list()

for(method_i in methods){
  print(method_i)
  mean_ests[[method_i]] = mean(as.numeric(subset(estimates, method == method_i)$est))
}


if (check_allele_distr(popn) < min_loci){
  pop = initialise_population(n_subpops = ns, loci = 200, t = 1000)
}

return(list(SS = mean_ests, pop = pop))








fst_v <- get_Fst_vec(pop$pA_post, N_v = rep(1000,20), 1000)
r2 <- apply(pop$r_post, MARGIN = c(2), mean)

plot(1:(1000), (ma(fst_v, n=1)), type = "l" ); abline(v=c(eqm1, start2, eqm2), col = 'red')
plot(1:(1001), ma(r2, n=1), type = "l" ); abline(v=c(eqm1, start2, eqm2), col = 'red')







plot(1:(sum(ts)), ma(diff(r2), n=20), type = "l" );abline(h=0); abline(v=cumsum(ts), col = 'red')
plot(2:(sum(ts)), ma(diff(fst_v), n=15), type = "l" );abline(h=0); abline(v=cumsum(ts), col = 'red')


ms = c(0.1, 0.05, 0.01)
ms = c(0.03,0.02,0.01)

Ms = list()
i=1
for (m in ms){
  Ms[[i]] <- construct_M_waples(m,20)
  i=i+1
}

ts = c(200,200,200)

init_pop <- initialise_population(n_subpops = 20, loci = 200, t = sum(ts))
prop_pop <- propegate_population(init_pop, N_v = rep(1000, 20),M_list = Ms, t_v = ts, c=0.05)


fst_v <- get_Fst_vec(prop_pop$pB_post, N_v = rep(1000,20), sum(ts))
r2 <- apply(prop_pop$r_post, MARGIN = c(2), mean)

plot(1:(sum(ts)), ma(fst_v, n=1), type = "l" ); abline(v=ts, col = 'red')
plot(0:(sum(ts)), ma(r2, n=1), type = "l" ); abline(v=cumsum(ts), col = 'red')

plot(1:(sum(ts)), ma(diff(r2), n=20), type = "l" );abline(h=0); abline(v=cumsum(ts), col = 'red')
plot(2:(sum(ts)), ma(diff(fst_v), n=15), type = "l" );abline(h=0); abline(v=cumsum(ts), col = 'red')

plot(2:(300), ma(diff(fst_v),n=2), type = "l" )

