# taking samples from multiple subpopulations
setwd('Desktop/VectorSim/scripts')

# saving and loading computed results - not simulated data
rm('real_out_list', 'sim_out')

#save.image(file = "../Data/cross_val_computed.RData")
load( "../Data/cross_val_computed.RData")

# loading packages etc
{
graphics.off()
rm(list = ls())
library(gridExtra)
library(boot)
library(abc)
library(abind)
library(Rcpp)
library(RcppArmadillo)
library(mc2d)
library(rdist)
library(pls)
library(bestNormalize)
library(scales)
library(foreach)
library(parallel)
library(doParallel)
source("simulator_functions.R")
propegate_population_inner = propegate_population_inner_C
propegate_population = propegate_population_decomposed
gc()
}

gc()

#  FIGURES  ####################
#
#  Bias
#
#
#  Across N,m (m = 0.05, N = 2000)
#  
#
#  CV error - sample size pooled sequencing - (ns, nt) = (2,3)
#
#
#  Credible intervals table
#
#
#  Optimised over space and time (m < 0.1, m > 0.1)
#
#
#  PCA colinearity of statistics
#
#
#
#  SI  ########################
#
#  Across N,m (m = 0.2, N = 10000)
#
#
#  Across N,m (m = 0,05, N = 2000) - all AF estimators
#
#
#  CV error - number of sampled subpopulations
#
#
#  Product Ne.m estimation 
#
#
#  Across N,m  (m = 0.2, N = 10000) - decreasing sample size AF estimator
#
#
#  Across m, single subpopulation infinite size
#
#
#  Accepted values plot
# 
#
#  Posterior histograms
#
#

####################################################
#   
#   1) Single sample, single subpopulation, N = 2,000
#   2) Multiple subpopulations, N = 2,000
#   3) Single subpopulation pool-seq - decreasing sample size
#   4) 3 subpopulations pool-seq - decreased sample size


#####################################
#####################################
######## Simulations ################
#####################################
#####################################

binom.test.one = function(successes, trials = 1000, p = 0.95){
  return (choose(trials, successes)*(p)^successes*(1-p)^(trials - successes))
}

binom.test = function(successes, trials, p){
  if (successes/trials == p) {return(.5)}
  
  if (successes/trials < p){
    tot_prob = 0
    for (s in 1:successes){
      tot_prob = tot_prob + binom.test.one(s, trials, p)
    }
    return(tot_prob)
  }
  
  if (successes/trials > p){
    tot_prob = 0
    for (s in successes:trials){
      tot_prob = tot_prob + binom.test.one(s, trials, p)
    }
    return(tot_prob)
  }
}


binom.test(960,1000, 0.95)

binom.test(964,1000, 0.95)*2
binom.test(935,1000, 0.95)*2

{

#info
{  loci = 1000; max_loci = ( loci + 5 ) * 1.1; min_loci = loci; n_loci = loci
sample_gens = c(0,5,10,20)
S_base = 1/(  matrix(c(1,2,3,5,9)) %*% c(1:4) )
S_base = c(S_base)

S_vec = get_S_vec(S_base, c(1000,200))[[1]]

S_vec_list = get_S_vec(S_base, c(1000,200))[[2]]

recom_rates = c(.50, .35, .20, .05)
Mfn = construct_M_4_step
grid_size_sim = 81
sample_subpops_sim = get_square_around_centre(81,3)
focal_pop = 41

grid_size_real = 121

sample_subpops_real = get_square_around_centre(121,3)

}


# testing a single N, m value
{

loci = 1000; max_loci = ( loci + 5 ) * 1.1; min_loci = loci; n_loci = loci
sample_gens = c(0,5,10,15,20)
recom_rates = c(.50, .35, .20, .05)
Mfn = construct_M_4_step
grid_size = 81
sample_subpops = get_square_around_centre(81,3)
focal_pop = 41
S_vec = c(100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 4500, 8000, 18000, Inf)
S_vec = c(S_vec, 2*S_vec)
pop = initialise_population(n_subpops = grid_size, loci = loci)

m = 0.9
N = 7000
Mfn = construct_M_4_step

model_fn <- get_model_fn_multi_pops_focal(sample_gens, recom_rates, Mfn, S_vec, grid_size_sim, sample_subpops_sim, max_loci, min_loci, n_loci, focal_pop)
Nm_list[[i]] =  Nm_draw_large(n_samp)
print('started')
Rprof()
o <- model_fn(N, m, f = get_SS_multi_pops3)
Rprof(NULL)
print(summaryRprof())

}

########### make data ############  
  
  
# generate simulated data
{

  
  n_par = 8
  n_samp = 250
  
  # constructing some things
  {
    model_fn_list = list()
    Nm_list = list()
    for (i in 1:n_par){
      model_fn_list[[i]] = get_model_fn_multi_pops(sample_gens, recom_rates, Mfn, S_vec, grid_size_sim, sample_subpops_sim, max_loci, min_loci, n_loci)
      Nm_list[[i]] =  Nm_draw_large(n_samp)
    }
  }
  
  
  s = proc.time()
  sim_out <- generate_ABC_sim_parallel_multi(loci, sample_gens, recom_rates, Nm_list, model_fn_list, n_par, S_vec, Mfn, grid_size_sim, sample_subpops_sim)
  e = proc.time(); print(e-s)
  
  saveRDS(sim_out, '../data/sim_out_wider_7.rds')
  sim_out <- readRDS('../data/sim_out_wider_7.rds')
  
  Nm_values_sim = sim_out$Nm_values
  Nm_values_sim = cbind(Nm_values_sim, apply(Nm_values_sim,1,prod))
}
sim_out = 0; gc()

# generate sim data - new S values
{
get_sim_info = function(sim_list){
  s = dimnames(sim_list$r2_sim)
  sizes =  as.numeric(substr(s[[1]], 2, nchar(s[[1]])))
  gens = as.numeric(substr(s[[2]], 2, nchar(s[[2]])))
  recoms = as.numeric(substr(s[[3]], 2, nchar(s[[3]])))
  subpops = as.numeric(substr(s[[4]], 3, nchar(s[[4]])))
  
  return(list(S = sizes, g = gens, c = recoms, sp = subpops))
}

trim_sim_out = function(sim_out, n){
  o = list()
  o[['Nm_values']] = sim_out[[1]][1:n,]
  for (i in 2:5){
    o[[i]] = asub(sim_out[[i]], idx = 1:n, dims = length(dim(sim_out[[i]])))
  }
  names(o) = names(sim_out)
  return(o)
}

add_S_fname = function(sim_out_original_fname, loci, grid_size_sim, S_vec_new){
  soi = readRDS(sim_out_original_fname)
  info = get_sim_info(soi)
  Nm_values_sim_original = soi$Nm_values
  
  soi = NULL; gc()
  
  max_loci = ( loci + 5 ) ; min_loci = loci; n_loci = loci
  sample_gens = info$g
  recom_rates = info$c
  Mfn = construct_M_4_step
  sample_subpops_sim = info$sp
  
  print(sample_gens)
  
  n_par = 8
  
  {
    Ns_l = split(Nm_values_sim_original[,1],ceiling(seq_along(Nm_values_sim_original[,1]) / (dim(Nm_values_sim_original)[1]/8) ))
    ms_l = split(Nm_values_sim_original[,2],ceiling(seq_along(Nm_values_sim_original[,2]) / (dim(Nm_values_sim_original)[1]/8) ))
    
    
    model_fn_list = list()
    Nm_list = list()
    for (i in 1:n_par){
      model_fn_list[[i]] = get_model_fn_multi_pops(sample_gens, recom_rates, Mfn, S_vec_new, grid_size_sim, sample_subpops_sim, max_loci, min_loci, n_loci)
      Nm_list[[i]] = matrix(c(Ns_l[[i]], ms_l[[i]]), nc = 2)
    }
  }
  
  
  s = proc.time()
  sim_out_copy <- generate_ABC_sim_parallel_multi(loci, sample_gens, recom_rates, Nm_list, model_fn_list, n_par, S_vec_new, Mfn, grid_size_sim, sample_subpops_sim)
  e = proc.time(); print(e-s)
  
  
  
  # combine old and new S values
  sim_out_original = readRDS(sim_out_original_fname)
  
  sim_out_combined = list()
  sim_out_combined[[1]] = sim_out_original[[1]]
  for (i in 2:5){
    sim_out_combined[[i]] = abind(sim_out_original[[i]], sim_out_copy[[i]], along = 1)
    names(sim_out_combined) = c()
  }
  
  names(sim_out_combined) = names(sim_out_original)
  
  return(sim_out_combined)
}

add_S = function(sim_out_original, loci, grid_size, S_vec_new){
  info = get_sim_info(sim_out_original)
  Nm_values_sim_original = sim_out_original$Nm_values
  
  soi = NULL; gc()
  
  max_loci = ( loci + 5 ) ; min_loci = loci; n_loci = loci
  sample_gens = info$g
  recom_rates = info$c
  Mfn = construct_M_4_step
  sample_subpops = info$sp
  
  n_par = 8
  
  {
    Ns_l = split(Nm_values_sim_original[,1],ceiling(seq_along(Nm_values_sim_original[,1]) / (dim(Nm_values_sim_original)[1]/8) ))
    ms_l = split(Nm_values_sim_original[,2],ceiling(seq_along(Nm_values_sim_original[,2]) / (dim(Nm_values_sim_original)[1]/8) ))
    
    
    model_fn_list = list()
    Nm_list = list()
    for (i in 1:n_par){
      model_fn_list[[i]] = get_model_fn_multi_pops(sample_gens, recom_rates, Mfn, S_vec_new, grid_size, sample_subpops, max_loci, min_loci, n_loci)
      Nm_list[[i]] = matrix(c(Ns_l[[i]], ms_l[[i]]), nc = 2)
    }
  }
  
  
  s = proc.time()
  sim_out_copy <- generate_ABC_sim_parallel_multi(loci, sample_gens, recom_rates, Nm_list, model_fn_list, n_par, S_vec_new, Mfn, grid_size, sample_subpops)
  e = proc.time(); print(e-s)
  
  
  
  # combine old and new S values

  sim_out_combined = list()
  sim_out_combined[[1]] = sim_out_original[[1]]
  for (i in 2:5){
    sim_out_combined[[i]] = abind(sim_out_original[[i]], sim_out_copy[[i]], along = 1)
    names(sim_out_combined) = c()
  }
  
  names(sim_out_combined) = names(sim_out_original)
  
  return(sim_out_combined)
}


sim_out_original_fname = '../data/sim_out_81_newS.rds'; new_fname= paste0(substr(sim_original_fname,1, nchar(sim_original_fname)-4), '_combined.rds')
S_vec_new = c(500000, 100000, 200, 1000/9, 100, 500/9)
sim_out_combined = add_S_fname(sim_out_original_fname, loci = 1000, grid_size_sim = 81, S_vec_new)
saveRDS(sim_out_combined, new_fname)
sim_out_combined = NULL; gc()
}

rm('real_out_list', 'sim_out')

# generate real data
{
  reps  = 9
  n_par = 8
  focal_pop = 61

  # have a component with varying N
  
  N1 = c(seq(400,1800,200),seq(2000, 20000, 1000))

  m1 = rep(0.05, 27)
  

  # have a component with varying m in [0.01, 1.0], N = 10,000
  
  N2 = rep(10000, 19)
  m2 = rev(seq(0.05, 0.95, 0.05))
  
  # have a component with varying m in [0.01, 1.0], N = 2,000
  
  N3 = rep(2000, 19)
  m3 = rev(seq(0.05, 0.95, 0.05))
  
  # have a component with varying m in [0.01, 1.0], N = 2,000
  
  N4 = rep(2000, 19)
  m4 = rev(seq(0.015, 0.195, 0.01))
  

  
  # have a component with varying m in [0, 0.05], much denser, N = 10,000
  
  #N3 = rep(10000, 26)
  #m3 = rev(seq(0, 0.05, 0.002))
  
  Ns = c(N1, N2, N3,N4)
  ms = c(m1, m2, m3,m4)
  
  Ns = c(N4)
  ms = c(m4)
  
  Ns = Nm_values_real[,1]
  ms = Nm_values_real[,2]
  
  # constructing some things
  {
    Ns_l = split(Ns,ceiling(seq_along(Ns)/15))
    ms_l = split(ms,ceiling(seq_along(ms)/15))
    
    model_fn_list = list()
    Nm_list = list()
    for (i in 1:n_par){
      model_fn_list[[i]] = get_model_fn_multi_pops_focal(sample_gens, recom_rates, Mfn, S_vec, grid_size_real, sample_subpops_real, max_loci, min_loci, n_loci, focal_pop)
      Nm_list[[i]] = matrix(c(Ns_l[[i]], ms_l[[i]]), nc = 2)
    }
  }
  
  real_out_list = list()
  s = proc.time()
  for (i in 1:reps){
    print(i)
    real_out_list[[i]] = generate_ABC_sim_parallel_multi(loci, sample_gens, recom_rates, Nm_list, model_fn_list, n_par, S_vec, Mfn, grid_size_real, sample_subpops_real)
  }  
  e = proc.time(); print(e-s)
  
  
  saveRDS(real_out_list, '../data/real_out_full_long_3.rds')
  real_out_list <- readRDS('../data/real_out_full_long_3.rds')
  
  Nm_values_real = real_out_list[[1]]$Nm_values
  Nm_values_real = cbind(Nm_values_real, apply(Nm_values_real,1,prod))
  
}

# generate real data - new S values
{
real_out_original_fname = '../data/real_out_121_newS.rds'; new_fname = paste0(substr(real_out_original_fname,1, nchar(real_out_original_fname)-4), '_combined.rds');
real_out_list = readRDS(real_out_original_fname)
real_out_list_combined = list()
for (i in 1:length(real_out_list)){
  real_out_list_combined[[i]] = add_S(real_out_list[[i]], loci = 5, grid_size_real, S_vec_new)
}
saveRDS(real_out_list_combined, new_fname)
}

########### load data #############

# combine  sim_files 
{


## size 100
sim_files = c('../data/sim_out_100.rds', '../data/sim_out_100_2.rds')

## size 121
sim_files = c('../data/sim_out_121_newS.rds')

## size 81
sim_files = c('../data/sim_out_81.rds', '../data/sim_out_81_2.rds')
sim_files = c('../data/sim_out_81_newS.rds', '../data/sim_out_81_newS_2.rds','../data/sim_out_81_newS_3.rds','../data/sim_out_81_newS_4.rds')
sim_files = c('../data/sim_out_full_S_1.rds', '../data/sim_out_full_S_2.rds','../data/sim_out_full_S_3.rds')
sim_files = c('../data/sim_out_final_1.rds', '../data/sim_out_final_2.rds')
sim_files = c('../data/sim_out_all_pops_1.rds')
sim_files = c('../data/sim_out_more_S_1.rds')
sim_files = c('../data/sim_out_wider_1.rds','../data/sim_out_wider_2.rds','../data/sim_out_wider_3.rds','../data/sim_out_wider_4.rds','../data/sim_out_wider_5.rds','../data/sim_out_wider_6.rds', '../data/sim_out_wider_7.rds')

### none yet

sim_out = list(); first = T
for (f in sim_files){
  if (first) {sim_out = readRDS(f); first = F}
  else{
    sim_i = readRDS(f)
    gc()
    sim_out$Nm_values = abind( sim_out$Nm_values,  sim_i$Nm_values, along = 1)
    for (i in 2:5){
      sim_out[[i]] = abind(sim_out[[i]], sim_i[[i]], along = length(dim(sim_out[[i]])))
      gc()
    }
  }
}

Nm_values_sim = sim_out$Nm_values
gc()
}
sim_i = NULL; gc()

# combine parallel real_list
{
  
### size 81
real_files = c('../data/real_out_81_newS_1.rds')
  
## size 100
real_files = c('../data/real_out_100.rds')

## size 121
real_files = c('../data/real_out_121_newS.rds','../data/real_out_121_newS_2.rds','../data/real_out_121_newS_3.rds','../data/real_out_121_newS_4.rds','../data/real_out_121_newS_5.rds')
real_files = c('../data/real_out_full_S_1.rds','../data/real_out_full_S_2.rds')
real_files = c('../data/real_out_final_1.rds')
real_files = c('../data/real_out_more_S_1.rds','../data/real_out_more_S_2.rds','../data/real_out_more_S_3.rds','../data/real_out_more_S_4.rds')

# full length
real_files = c('../data/real_out_full_long.rds','../data/real_out_full_long_2.rds', '../data/real_out_full_long_3.rds' )

real_out_list <- list()
for (f in real_files){
  real_out_list = c(real_out_list, readRDS(f))
}

Nm_values_real = real_out_list[[1]]$Nm_values

}

gc()
# add new real values
{
real_out_list1 = real_out_list
real_out_list1 = readRDS('../data/real_out_more_S_combined.rds')
real_out_list1 = real_out_list2
real_out_list2 = real_out_list
real_out_list_combined = list()

for (i in 1:min(length(real_out_list1), length(real_out_list2))){
  combined_sub = list()
  real_out_list_combined[[i]] = list()
  real_out_list_combined[[i]][[1]] = combine(real_out_list1[[i]][[1]], real_out_list2[[i]][[1]], axis = 1, M1_index = 95, M2_slice = 1:19)
  for (j in 2:length(real_out_list1[[i]])){
    real_out_list_combined[[i]][[j]] = combine(real_out_list1[[i]][[j]], real_out_list2[[i]][[j]], axis = length(dim(real_out_list[[i]][[j]])), M1_index = 95, M2_slice = 1:19 )
  }
  names(real_out_list_combined[[i]]) = names(real_out_list1[[i]])
}


for (i in 1:length(real_out_list_combined)){
  names(real_out_list_combined[[i]]) = names(real_out_list1[[i]])
}

saveRDS(real_out_list_combined, '../data/real_out_full_long.rds')
real_out_list = readRDS('../data/real_out_full_long.rds')
Nm_values_real = real_out_list[[1]]$Nm_values
rm('real_out_list_combined', 'real_out_list1', 'real_out_list2')
gc()

}



# filtering
{
# filter m > 0.8
{
  sim_out_filt = list()
  filt = Nm_values_sim[,2] < 0.8
  sim_out_filt[['Nm_values']] = sim_out[[1]][filt,]
  for (i in 2:length(sim_out)){
    sim_out_filt[[i]] = asub(sim_out[[i]], which(filt), dims = length(dim(sim_out[[i]])))
    sim_out[[i]] = 'null'
    gc()
  }
  gc()
  
  for (i in 1:length(sim_out)){
    names(sim_out_filt)[i] = names(sim_out)[i]
  }
  
  Nm_values_sim_filt = sim_out_filt$Nm_values
  
  gc()
}

# filter m > 0.8
{
real_out_list_filt = list(); i = 1
for (l in real_out_list){
  Nms = l$Nm_values
  filt = Nms[,2] < .8
  ll = list()
  ll[[1]] = l[[1]][filt,]
  for (j in 2:length(l)){
    tr = asub(l[[j]], which(filt), dims = length(dim(l[[j]])))
    ll[[j]] =  tr
  }
  real_out_list_filt[[i]] = ll
  i=i+1
}

for (i in 1:length(real_out_list_filt)){
  for (j in 1:length(real_out_list_filt[[i]])){
    names(real_out_list_filt[[i]])[j] = names(real_out_list[[i]])[j]
  }
}
  
Nm_values_real_filt = real_out_list_filt[[1]]$Nm_values
  

gc()
}

real_out_list = NULL; gc()

#sim_i = NULL; sim_out = NULL; sim_out_filt = NULL; real_out_list= NULL
gc()
}
}


##########################################################################
###################
###################   ABC and results
###################
##########################################################################


#####################################
#####################################
######## Best ABC method ############
#####################################
#####################################
{
  
# PLS vs prior truncation  vs PCA vs no reduction
  
# testing with 3 gens, 3 subpops
SS_3.3 = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100,30*round(1000/3/3),30)), gens_select = c(0,10,20), subpops_select_sim = sample_subpops_sim, subpops_select_real = 60:62, assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_3.3= SS_3.3$sim; SS_real_3.3 = SS_3.3$real

SS_3.3.ld = make_SS(sim_out, real_out_list, S_select = round(1000/3/3), gens_select = c(0,10,20), subpops_select_sim = 40:42, subpops_select_real = 60:62, assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_3.3 = SS_3.3.ld$sim; SS_real_3.3.ld = SS_3.3$real

selected = select_bounded(50, Nm_values_sim, m_range = c(0,0.8), N_range = c(400,22500))

o1 = cross_validate(do_abc_PLS, SS_sim_3.3, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1, sizenet = 4, numnet = 1)
o2 = cross_validate(do_abc_PCA, SS_sim_3.3, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1, sizenet = 4, numnet = 1)
o3 = cross_validate(do_abc_trunc, SS_sim_3.3, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1, sizenet = 4, numnet = 1)


o1 = cross_validate(do_abc_PLS, SS_sim_3.3, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1, sizenet = 4, numnet = 1)
o2 = cross_validate(do_abc_PLS, SS_sim_3.3, selected, n_comps =  5, Nm_values_sim, method = 'loclinear', tol = .1, sizenet = 4, numnet = 1)
o3 = cross_validate(do_abc_PLS, SS_sim_3.3, selected, n_comps =  5, Nm_values_sim, method = 'rejection', tol = .1, sizenet = 4, numnet = 1)



o1 = cross_validate(do_abc_PLS, SS_sim_3.3, selected, n_comps =  3, Nm_values_sim, method = 'neuralnet', tol = .1, sizenet = 4, numnet = 1)
o2 = cross_validate(do_abc_PLS, SS_sim_3.3, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1, sizenet = 4, numnet = 1)
o3 = cross_validate(do_abc_PLS, SS_sim_3.3, selected, n_comps =  7, Nm_values_sim, method = 'neuralnet', tol = .1, sizenet = 4, numnet = 1)



RMS_cv1 = RMS_err(o1); colMeans(RMS_cv1)
RMS_cv2 = RMS_err(o2); colMeans(RMS_cv2)
RMS_cv3 = RMS_err(o3); colMeans(RMS_cv3)


SS_unequal = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100,round(1000/3)*30,30)), gens_select = c(0,10,20), subpops_select_sim = 41, subpops_select_real = 61, assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_unequal= SS_unequal$sim; SS_real_unequal = SS_unequal$real

selected = sample(1:dim(Nm_values_sim)[1], 50)

o1 = cross_validate(do_abc_PLS, SS_sim_equal, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1)

o2 = cross_validate(do_abc_PLS, SS_sim_unequal, selected, n_comps =  2, Nm_values_sim, method = 'loclinear', tol = .05)
o3 = cross_validate(do_abc_PLS, SS_sim_unequal, selected, n_comps =  2, Nm_values_sim, method = 'neuralnet', tol = .05)

o4 = cross_validate(do_abc_PLS, SS_sim_unequal, selected, n_comps =  4, Nm_values_sim, method = 'neuralnet', tol = .2)
o5 = cross_validate(do_abc_PLS, SS_sim_equal, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .2)
o5.ll = cross_validate(do_abc_trunc, SS_sim_equal, selected, n_comps =  5, Nm_values_sim, method = 'loclinear', tol = .1)

o5.t = cross_validate(do_abc_trunc, SS_sim_equal, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .2)

o10 = cross_validate(do_abc_PLS, SS_sim_equal, selected, n_comps =  10, Nm_values_sim, method = 'neuralnet', tol = .2)
o10.trunc = cross_validate(do_abc_trunc, SS_sim_equal, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .2)
o10.15 = cross_validate(do_abc_PLS, SS_sim_equal, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .2, sizenet = 10)



o4 = cross_validate(do_abc_PLS, SS_sim_equal, selected, n_comps =  4, Nm_values_sim, method = 'neuralnet', tol = .2)





RMS_cv1 = RMS_err(o1); colMeans(RMS_cv1)
RMS_cv2 = RMS_err(o2); colMeans(RMS_cv2)
RMS_cv3 = RMS_err(o3); colMeans(RMS_cv3)
RMS_cv4 = RMS_err(o4); colMeans(RMS_cv4)
RMS_cv5 = RMS_err(o5); colMeans(RMS_cv5)
RMS_cv5 = RMS_err(o5.t); colMeans(RMS_cv5)
RMS_cv5 = RMS_err(o5.ll); colMeans(RMS_cv5)

RMS_cv10.15 = RMS_err(o10.15); colMeans(RMS_cv10.15)

RMS_cv10 = RMS_err(o10); colMeans(RMS_cv10)


o2 = cross_validate(do_abc_PLS, SS_sim_unequal, selected, n_comps =  3, Nm_values_sim, method = 'loclinear', tol = .2)





o7 = cross_validate(do_abc_trunc, SS_sim_two, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1)
o8 = cross_validate(do_abc_trunc, SS_sim_two, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .2)

oA = cross_validate(do_abc_trunc, SS_sim_two, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1)
oB = cross_validate(do_abc_PLS, SS_sim_two, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1)

RMS_cvA = RMS_err(oA); colMeans(RMS_cvA)
RMS_cvB = RMS_err(oB); colMeans(RMS_cvB)

RMS_cv1 = RMS_err(o1); colMeans(RMS_cv1)
RMS_cv2 = RMS_err(o2); colMeans(RMS_cv2)
RMS_cv3 = RMS_err(o3); colMeans(RMS_cv3)
RMS_cv4 = RMS_err(o4); colMeans(RMS_cv4)

RMS_cv5 = RMS_err(o5); colMeans(RMS_cv5)
RMS_cv6 = RMS_err(o6); colMeans(RMS_cv6)
RMS_cv7 = RMS_err(o7); colMeans(RMS_cv7)
RMS_cv8 = RMS_err(o8); colMeans(RMS_cv8)
RMS_cv9 = RMS_err(o9); colMeans(RMS_cv9)
}


method = 'neuralnet'
nc = 5
tol = .07
sizenet = 4
numnet = 1

t_end = Sys.time(); print(Sys.time())
#####################################
#####################################
######## Single sample ##############
#####################################
#####################################
{
  
# across m, N = 2000
{

  # do ABC
  {
    o_single_2000 = list()
    o_single_2000_full = list()
      
    i=1
    for (S in c(100, 500,1000,1550,3100,7750,15500,Inf)){
      SS = make_SS(sim_out, real_out_list, S_select = S, gens_select = c(0), subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
    
      one_m_single = test_specific(do_abc_PLS, SS_sim, SS_real, n_comps = min(nc, dim(SS_sim)[2]), Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol, method = method, sizenet = sizenet, numnet = numnet )
    
      o_single_2000_full[[i]] = one_m_single
      i=i+1
    }
    
    saveRDS(o_single_2000_full, '../data/o_single_2000_full.RDS')
  }
  
  # make plot
  {
  i=1
  for (l in o_single_2000_full){
    o_single_2000[[i]] = apply(l, c(1,2,3), mean) 
    i=i+1
  }
  
  pdf('../final_figures/single_2000.pdf', width = 13, height = 10)
  #plot_compare_m(o_single_2000, legend = paste('S =', c(500,1000,7750,15500)), cols = c('black', 'black'), ltys = 1:100)
  plot_compare_m_neat(o_single_2000[c(2,3,5,6)], legend = paste('S =', c(500,1000,7750,15500), ' '), cols = c('black', 'black'), ltys = 1:100)
  
  graphics.off()
  }
}

# across m,  N = 10000
{
  # do ABC
  {
    o_single_10000 = list()
    o_single_10000_full = list()
    
    i=1
    for (S in c(100, 500,1000,1550,3100,7750,15500,Inf)){
      SS = make_SS(sim_out, real_out_list, S_select = S, gens_select = c(0), subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = F, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      
      one_m_single = test_specific(do_abc_PLS, SS_sim, SS_real, n_comps = min(nc, dim(SS_sim)[2]), Nm_values_sim, Nm_values_real, idx = c(31,62), tol = tol, method = method,sizenet = sizenet, numnet = numnet )
      
      o_single_10000_full[[i]] = one_m_single
      i=i+1
    }
    saveRDS(o_single_10000_full, '../data/o_single_10000_full.RDS')
    
  }
  
  # make plot
  {
  i=1
  for (l in o_single_10000_full){
    o_single_10000[[i]] = apply(l, c(1,2,3), mean) 
    i=i+1
  }

  pdf('../final_figures/single_10000.pdf', width = 7, height = 7)
  idxs = c(2,4,5,8)
  plot_compare_m_neat(o_single_10000[idxs], legend = paste('S =', c(100, 500,1000,1550,3100,7750,15500,Inf)[idxs]), cols = c('black', 'black'), ltys = 1:100)
  graphics.off()
  }
}


# across N, m = 0.05
{
  # do ABC
  {
    o_single_N = list()
    o_single_N_full = list()
    
    i=1
    for (S in c(100, 500,1000,1550, 3100,7750,15500, Inf)){
      SS = make_SS(sim_out, real_out_list, S_select = S, gens_select = c(0), subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      
      one_m_single = test_specific(do_abc_PLS, SS_sim, SS_real, n_comps = min(nc, dim(SS_sim)[2]), Nm_values_sim, Nm_values_real, idx = c(1,27), tol = tol, method = method, sizenet = sizenet, numnet = numnet )
      
      o_single_N_full[[i]] = one_m_single
      i=i+1
    }
    
    saveRDS(o_single_N_full, '../data/o_single_N_full.RDS')
    
  }
  
  # make plot
  {
  i=1
  for (l in o_single_N_full){
    o_single_N[[i]] = apply(l, c(1,2,3), mean)
    i=i+1
  }
  pdf('../final_figures/single_vs_N.pdf', height = 7, width = 7, pointsize = 7)
  plot_compare_N_neat(o_single_N[c(2,3,5,7)], legend = paste('S =', c(500,1000,7750,Inf)), cols = c('black', 'black'), ltys = 1:100)
  graphics.off()
  }
}
  

# across N, m = 0.2
{
    # do ABC
    {
      o_single_N020 = list()
      o_single_N_full020 = list()
      
      i=1
      for (S in c(100, 500,1000,1550, 3100,7750,15500, Inf)){
        SS = make_SS(sim_out, real_out_list, S_select = S, gens_select = c(0), subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
        SS_sim = SS$sim; SS_real = SS$real
        
        one_m_single = test_specific(do_abc_PLS, SS_sim, SS_real, n_comps = min(nc, dim(SS_sim)[2]), Nm_values_sim, Nm_values_real, idx = c(63,76), tol = tol, method = method, sizenet = sizenet, numnet = numnet )
        
        o_single_N_full020[[i]] = one_m_single
        i=i+1
      }
      
      saveRDS(o_single_N_full020, '../data/o_single_N_full_020.RDS')
      
    }
    
    # make plot
    {
      i=1
      for (l in o_single_N_full020){
        o_single_N020[[i]] = apply(l, c(1,2,3), mean)
        i=i+1
      }
      pdf('../final_figures/single_vs_N020.pdf', height = 7, width = 7, pointsize = 7)
      plot_compare_N_neat(o_single_N020[c(2,3,5,7)], legend = paste('S =', c(500,1000,7750,Inf)), cols = c('black', 'black'), ltys = 1:100)
      graphics.off()
    }
  }

#combined plot  
{

  legends = c(paste('S =', c(500, 1550, 3100, Inf)), 'Prior')
  ltys = c(1,2,3,4,1)
  cols = c(rep('darkslategrey',4), 'green')
  
  N_plot = o_single_N[c(2,3,4,6)]; N_plot[[5]] =  guess_N_mean
  m_plot = o_single_2000[c(2,3,4,6)]; m_plot[[5]] =  guess_m_mean_N2000

  plot_compare_N_m_neat(N_plot, m_plot, legends, ltys = ltys, cols = cols, fname = '../final_figures/single_N_and_m.pdf', prior = T)
}
  
# cross validate
{
  # do CV
  {
  S_total = 1000

  selected_small = select_bounded(100, Nm_values_sim, m_range = c(0,0.1), N_range = c(400,22500))
  RMSE_out_small = array(NA, dim = c(0,2 ))
  RMSE_out_list_small = list(); i = 1
  
  sel = c(seq(1,25,2), 26:35)
  sel = c(seq(1,25,6), 26, 30, 34, 35)
  for (S in S_vec[S_vec>=100][sel]){
  
      SS = make_SS(sim_out, real_out_list, S_select = S, gens_select = c(0), subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      cv = cross_validate(do_abc_PLS, SS_sim, selected_small, n_comps =  min(nc, dim(SS_sim)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
      RMS_cv = RMS_err(cv); means = colMedians(RMS_cv)
      if (NaN %in% RMS_cv){
        break
      }
      RMSE_out_list_small[[i]] = RMS_cv; i = i+1
      RMSE_out_small = rbind(RMSE_out_small,means)
    
  }
  
  
  saveRDS(RMSE_out_list_small, '../data/RMSE_out_list_single_small.RDS')
  RMSE_out_list_small = readRDS('../data/RMSE_out_list_single_small.RDS')
  }
  
  {
    S_total = 1000
    
    selected_large = select_bounded(100, Nm_values_sim, m_range = c(0.1,0.8), N_range = c(400,22500))
    RMSE_out_large = array(NA, dim = c(0,2 ))
    RMSE_out_list_large = list(); i = 1
    
    sel = c(seq(1,25,2), 26:35)
    sel = c(seq(1,25,6), 26, 30, 34, 35)
    for (S in S_vec[S_vec>=100][sel]){
      
      SS = make_SS(sim_out, real_out_list, S_select = S, gens_select = c(0), subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      cv = cross_validate(do_abc_PLS, SS_sim, selected_large, n_comps =  min(nc, dim(SS_sim)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
      RMS_cv = RMS_err(cv); means = colMedians(RMS_cv)
      if (NaN %in% RMS_cv){
        break
      }
      RMSE_out_list_large[[i]] = RMS_cv; i = i+1
      RMSE_out_large = rbind(RMSE_out_large,means)
      
    }
    
    
    saveRDS(RMSE_out_list_large, '../data/RMSE_out_list_single_large.RDS')
    RMSE_out_list_large = readRDS('../data/RMSE_out_list_single_large.RDS')
  }
  
  # make plot
  
  pdf('../final_figures/CV_single.pdf', height = 7, width = 7, pointsize = 7)
  par(mfrow = c(2,1))
  {
  RMSE_out_large = matrix(unlist(lapply(RMSE_out_list_large,colMeans)), nc=2, byrow = T)


  
  max_y = 0
  i = 1
  for (x in RMSE_out_list_large){
      ci = do_boot(x[,2])
      max_y = max(max_y, ci)
    i = i+1
  }
  range(max_y,RMSE_out_large[,],1)
  plot(NA,log = 'x', xlab = expression(S), ylab = 'Relative Error', ylim = c(0,1), xlim = c(100,16000)); abline(h = RMSE_out_large[length(sel),1], lty = 5); axis(1, c(100,200,500,1000,2000,5000,10000,16000), labels = T)
  
  i = 1
  for (x in RMSE_out_list_large){
    ci = do_boot(x[,1])
    lines(rep((S_vec[S_vec>=100][sel][i]),2), ci, lw = .5)
    i = i+1
  }
  
  
  
  points(S_vec[S_vec>=100][sel], RMSE_out_large[,1],pch = 21,cex = 1, bg='white')
  
  
  i = 1
  for (x in RMSE_out_list_large){
    ci = do_boot(x[,2])
    lines(rep(S_vec[S_vec>=100][sel][i],2), ci, lw = .5, col = 'black')
    i = i+1
  }
  
  points(S_vec[S_vec>=100][sel], RMSE_out_large[,2],  pch = 25, cex = 1, col = 'black', bg = 'white'); abline(h = RMSE_out_large[length(sel),2], lty = 3, col = 'black')#, ylab = 'RMSE (m)', type = 'p',pch = 16,cex = .5, ylim = range(RMSE_out[,2])); abline(h = RMSE_out[1+sum(S_vec>100),2], lty = 2)

  legend('topright',cex = .9,inset = 0.02, legend = c('N', 'N (S = Inf)', 'm', 'm (S = Inf)' ), lty = c(NA, 5,NA,3), pch = c(1,NA,6, NA))
  }
  {
    RMSE_out_small = matrix(unlist(lapply(RMSE_out_list_small,colMeans)), nc=2, byrow = T)
    
    
    
    max_y = 0
    i = 1
    for (x in RMSE_out_list_small){
      ci = do_boot(x[,2])
      max_y = max(max_y, ci)
      i = i+1
    }
    range(max_y,RMSE_out_small[,],1)
    plot(NA,log = 'x', xlab = expression(S), ylab = 'Relative Error', ylim = c(0,1), xlim = c(100,16000)); abline(h = RMSE_out_small[length(sel),1], lty = 5); axis(1, c(100,200,500,1000,2000,5000,10000,16000), labels = T)
    
    i = 1
    for (x in RMSE_out_list_small){
      ci = do_boot(x[,1])
      lines(rep((S_vec[S_vec>=100][sel][i]),2), ci, lw = .5)
      i = i+1
    }
    
    
    
    points(S_vec[S_vec>=100][sel], RMSE_out_small[,1],pch = 21,cex = 1, bg='white')
    
    
    i = 1
    for (x in RMSE_out_list_small){
      ci = do_boot(x[,2])
      lines(rep(S_vec[S_vec>=100][sel][i],2), ci, lw = .5, col = 'black')
      i = i+1
    }
    
    points(S_vec[S_vec>=100][sel], RMSE_out_small[,2],  pch = 25, cex = 1, col = 'black', bg = 'white'); abline(h = RMSE_out_small[length(sel),2], lty = 3, col = 'black')#, ylab = 'RMSE (m)', type = 'p',pch = 16,cex = .5, ylim = range(RMSE_out[,2])); abline(h = RMSE_out[1+sum(S_vec>100),2], lty = 2)
    
    legend('topright',cex = .9,inset = 0.02, legend = c('N', 'N (S = Inf)', 'm', 'm (S = Inf)' ), lty = c(NA, 5,NA,3), pch = c(1,NA,6, NA))
  }
  graphics.off()
  
}

}

################################################
################################################
######## NUM SUBPOPS ###########################
################################################
################################################

{
# select an S per timepoint - try a couple
# investigate for a couple of 'over-time' sampling stategies
# look at how varying population size affects inference
# look at how varying migration rate affects inference
# look at what happens below the minimum simulated m 

## want to look for 'better with smaller sample size' type information - plain best strategy comes later


gens_select = c(0,10,20)
S_per_timepoint = 1000/3

S1 = round(S_per_timepoint)
S2 = round(S_per_timepoint / 2)
S3 = round(S_per_timepoint / 3)
S5 = round(S_per_timepoint / 5)
S9 = round(S_per_timepoint / 9)


# making statistics
{
v = c(1:2)
round(get_diploid_equivalent(1e100, S1*30, 30))
  
SS_one = make_SS(sim_out, real_out_list, S_select = S1, gens_select = gens_select, subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = T, LD_info = T,  focal_pop_sim = 61, focal_pop_real = 41, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_one = SS_one$sim[,]; SS_real_one = SS_one$real

SS_one_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S1*30, 30)), gens_select = gens_select, subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = T, LD_info = F,  focal_pop_sim = 61, focal_pop_real = 41, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_one_AF = SS_one_AF$sim[,]; SS_real_one_AF = SS_one_AF$real


SS_two = make_SS(sim_out, real_out_list, S_select = S2, gens_select = gens_select, subpops_select_sim = c(41,42), subpops_select_real = c(61,62), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_two = SS_two$sim[,]; SS_real_two = SS_two$real

SS_two_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S2*30, 30)), gens_select = gens_select, subpops_select_sim = c(41,42), subpops_select_real = c(61,62), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_two_AF = SS_two_AF$sim[,]; SS_real_two_AF = SS_two_AF$real


SS_three = make_SS(sim_out, real_out_list, S_select = S3, gens_select = gens_select, subpops_select_sim = c(40,41,42), subpops_select_real = c(60,61,62), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_three = SS_three$sim[,]; SS_real_three = SS_three$real

SS_three_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S3*30, 30)), gens_select = gens_select, subpops_select_sim = c(40,41,42), subpops_select_real = c(60,61,62), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_three_AF = SS_three_AF$sim[,]; SS_real_three_AF = SS_three_AF$real


SS_five = make_SS(sim_out, real_out_list, S_select = S5, gens_select = gens_select, subpops_select_sim = c(40,41,42,49,50), subpops_select_real = c(60,61,62,71,72), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_five = SS_five$sim[,]; SS_real_five = SS_five$real

SS_five_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S5*30, 30)), gens_select = gens_select, subpops_select_sim = c(40,41,42,49,50), subpops_select_real = c(60,61,62,71,72), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_five_AF = SS_five_AF$sim[,]; SS_real_five_AF = SS_five_AF$real


SS_nine = make_SS(sim_out, real_out_list, S_select = S9, gens_select = gens_select, subpops_select_sim = sample_subpops_sim, subpops_select_real = sample_subpops_real, assume_equal = T, LD_info = T, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_nine = SS_nine$sim[,]; SS_real_nine = SS_nine$real

SS_nine_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S9*30, 30)), gens_select = gens_select, subpops_select_sim = sample_subpops_sim, subpops_select_real = sample_subpops_real, assume_equal = T, LD_info = F, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_nine_AF = SS_nine_AF$sim[,]; SS_real_nine_AF = SS_nine_AF$real


}



# vary N, m = 0.05
{
  
  # run ABC
  {

    
  one_N = test_specific(do_abc_PLS, SS_sim_one, SS_real_one, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27), tol = tol, method = method, sizenet = sizenet, numnet = numnet)
  two_N = test_specific(do_abc_PLS, SS_sim_two, SS_real_two, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27), tol = tol , method = method, sizenet = sizenet, numnet = numnet)
  three_N = test_specific(do_abc_PLS, SS_sim_three, SS_real_three, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27), tol = tol , method = method, sizenet = sizenet, numnet = numnet)
  five_N = test_specific(do_abc_PLS, SS_sim_five, SS_real_five, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27), tol = tol, method = method, sizenet = sizenet, numnet = numnet)
  nine_N = test_specific(do_abc_PLS, SS_sim_nine, SS_real_nine, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27), tol = tol, method = method, sizenet = sizenet, numnet = numnet)
  
  one_N_AF = test_specific(do_abc_PLS, SS_sim_one_AF, SS_real_one_AF, n_comps = min(nc, dim(SS_sim_one_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(1,27), tol = tol, method = method, sizenet = sizenet, numnet = numnet)
  two_N_AF = test_specific(do_abc_PLS, SS_sim_two_AF, SS_real_two_AF, n_comps = min(nc, dim(SS_sim_two_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(1,27), tol = tol , method = method, sizenet = sizenet, numnet = numnet)
  three_N_AF = test_specific(do_abc_PLS, SS_sim_three_AF, SS_real_three_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27), tol = tol , method = method, sizenet = sizenet, numnet = numnet)
  five_N_AF = test_specific(do_abc_PLS, SS_sim_five_AF, SS_real_five_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
  nine_N_AF = test_specific(do_abc_PLS, SS_sim_nine_AF, SS_real_nine_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
  
  guess = just_guess(Nm_values_real, c(1,27))
  varyN_m0.05 = list(one_N, two_N, three_N, five_N, nine_N, one_N_AF, two_N_AF, three_N_AF, five_N_AF, nine_N_AF, guess)
  saveRDS(varyN_m0.05, '../data/varyN_m0.05')
  varyN_m0.05 = readRDS('../data/varyN_m0.05')
  
  }
  
  # processing
  varyN_m0.05_means = lapply(varyN_m0.05, apply, c(1,2,3), trim_mean,0)
    

  
  # make figure
  {
  pdf('../final_figures/n_sp_vs_N_m005.pdf', height = 7, width = 7)
  plot_compare_N_neat(list(one_N_mean,three_N_mean, nine_N_mean,one_N_mean_AF,three_N_mean_AF, nine_N_mean_AF), legend = c('1 (LD)  ', '3 (LD)  ', '9 (LD)  ', '1 (AF)  ', '3 (AF)  ', '9 (AF)  '), ltys = c(1,2,3,1,2,3), cols = c(rep('darkslategrey',3), rep('red',3)))
  graphics.off()
  }
  
  #old
  {
  plot_compare_N_neat(list(one_N_mean,two_N_mean,three_N_mean, five_N_mean, nine_N_mean), legend = c('1 sub', '2 sub', '3 sub', '5 sub', '9 sub'))
  plot_compare_N_neat(list(one_N_mean_AF,two_N_mean_AF,three_N_mean_AF, five_N_mean_AF, nine_N_mean_AF), legend = c('1 sub', '2 sub', '3 sub', '5 sub', '9 sub'))
  }
}

# vary N, m = 0.2
{
  
  # run ABC
  {
    one_N_mig = test_specific(do_abc_PLS, SS_sim_one, SS_real_one, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76), tol = tol, method = method, sizenet = sizenet, numnet = numnet)
    two_N_mig = test_specific(do_abc_PLS, SS_sim_two, SS_real_two, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76), tol = tol , method = method, sizenet = sizenet, numnet = numnet)
    three_N_mig = test_specific(do_abc_PLS, SS_sim_three, SS_real_three, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76), tol = tol , method = method, sizenet = sizenet, numnet = numnet)
    five_N_mig = test_specific(do_abc_PLS, SS_sim_five, SS_real_five, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
    nine_N_mig = test_specific(do_abc_PLS, SS_sim_nine, SS_real_nine, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76), tol = tol, method = method, sizenet = sizenet, numnet = numnet)
    
    one_N_AF_mig = test_specific(do_abc_PLS, SS_sim_one_AF, SS_real_one_AF, n_comps = min(nc, dim(SS_sim_one_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(63,76), tol = tol, method = method, sizenet = sizenet, numnet = numnet)
    two_N_AF_mig = test_specific(do_abc_PLS, SS_sim_two_AF, SS_real_two_AF, n_comps = min(nc, dim(SS_sim_two_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(63,76), tol = tol , method = method, sizenet = sizenet, numnet = numnet)
    three_N_AF_mig = test_specific(do_abc_PLS, SS_sim_three_AF, SS_real_three_AF, n_comps = 2, Nm_values_sim, Nm_values_real, idx = c(63,76), tol = tol , method = method, sizenet = sizenet, numnet = numnet)
    five_N_AF_mig = test_specific(do_abc_PLS, SS_sim_five_AF, SS_real_five_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
    nine_N_AF_mig = test_specific(do_abc_PLS, SS_sim_nine_AF, SS_real_nine_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
    
    guess = just_guess(Nm_values_real, c(63,76))
    
    varyN_m0.20 = list(one_N_mig, two_N_mig, three_N_mig, five_N_mig, nine_N_mig, one_N_AF_mig, two_N_AF_mig, three_N_AF_mig, five_N_AF_mig, nine_N_AF_mig, guess)
    saveRDS(varyN_m0.20, '../data/varyN_m0.20.RDS')
    varyN_m0.20 = readRDS('../data/varyN_m0.20.RDS')
    
    
  }
  
  # processing
  varyN_m0.20_means = lapply(varyN_m0.20, apply, c(1,2,3), mean)

  
  # make figure
  {
    pdf('../final_figures/n_sp_vs_N_m020.pdf', height = 7, width = 7, pointsize = 7)
    plot_compare_N_neat(varyN_m0.20_means[c(1,2,6,7,11)], legend = c('1 (Ind)  ', '2 (Ind)  ', '1 (Pool)  ', '2 (Pool)  ', 'Prior  '), ltys = c(1,2,1,2,1), cols = c(rep('darkslategrey',2), rep('red',2), '#4daf4a'))
    graphics.off()
  }
  
}

# vary m, N = 10000
{
  # run ABC
  {
  one_m = test_specific(do_abc_PLS, SS_sim_one, SS_real_one, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(31,62), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
  two_m = test_specific(do_abc_PLS, SS_sim_two, SS_real_two, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(31,62), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
  three_m = test_specific(do_abc_PLS, SS_sim_three, SS_real_three, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(31,62), tol = tol , method = method , sizenet = sizenet, numnet = numnet)
  five_m = test_specific(do_abc_PLS, SS_sim_five, SS_real_five, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(31,62), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
  nine_m = test_specific(do_abc_PLS, SS_sim_nine, SS_real_nine, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(31,62), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
  
  one_m_AF = test_specific(do_abc_PLS, SS_sim_one_AF, SS_real_one_AF, n_comps =  min(nc, dim(SS_sim_one_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(31,62), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
  two_m_AF = test_specific(do_abc_PLS, SS_sim_two_AF, SS_real_two_AF, n_comps =  min(nc, dim(SS_sim_two_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(31,62), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
  three_m_AF = test_specific(do_abc_PLS, SS_sim_three_AF, SS_real_three_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(31,62), tol = tol , method = method , sizenet = sizenet, numnet = numnet)
  five_m_AF = test_specific(do_abc_PLS, SS_sim_five_AF, SS_real_five_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(31,62), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
  nine_m_AF = test_specific(do_abc_PLS, SS_sim_nine_AF, SS_real_nine_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(31,62), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
  
  guess = just_guess(Nm_values_real, idxs = c(31,62))
  
  varym_N10000 = list(one_m, two_m, three_m, five_m, nine_m, one_m_AF, two_m_AF, three_m_AF, five_m_AF, nine_m_AF, guess)
  saveRDS(varym_N10000, '../data/varym_N10000.RDS')
}
  
  # processing
  varym_N10000_means = lapply(varym_N10000, apply, c(1,2,3), mean)
  
  # make figures
  {
  pdf('../final_figures/n_sp_vs_m_10000.pdf', height = 7, width = 7, pointsize = 10.5)
  plot_compare_m_neat(varym_N10000_means[c(1,2,6,7,11)], legend = c('1 (Ind)  ', '2 (Ind)  ', '1 (Pool)  ', '2 (Pool)  ', 'Prior'), ltys = c(1,5,1,5,1), cols = c(rep('darkslategrey',2), rep('red',2), '#4daf4a'))
  graphics.off()
  }
}

# vary m, N = 2000
{
  # run ABC
  {
    one_m_N2000 = test_specific(do_abc_PLS, SS_sim_one, SS_real_one, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
    two_m_N2000 = test_specific(do_abc_PLS, SS_sim_two, SS_real_two, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
    three_m_N2000= test_specific(do_abc_PLS, SS_sim_three, SS_real_three, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol , method = method , sizenet = sizenet, numnet = numnet)
    five_m_N2000 = test_specific(do_abc_PLS, SS_sim_five, SS_real_five, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
    nine_m_N2000= test_specific(do_abc_PLS, SS_sim_nine, SS_real_nine, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
    
    one_m_N2000_AF = test_specific(do_abc_PLS, SS_sim_one_AF, SS_real_one_AF, n_comps =  min(nc, dim(SS_sim_one_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
    two_m_N2000_AF = test_specific(do_abc_PLS, SS_sim_two_AF, SS_real_two_AF, n_comps =  min(nc, dim(SS_sim_two_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
    three_m_N2000_AF = test_specific(do_abc_PLS, SS_sim_three_AF, SS_real_three_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol , method = method , sizenet = sizenet, numnet = numnet)
    five_m_N2000_AF = test_specific(do_abc_PLS, SS_sim_five_AF, SS_real_five_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
    nine_m_N2000_AF = test_specific(do_abc_PLS, SS_sim_nine_AF, SS_real_nine_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol, method = method , sizenet = sizenet, numnet = numnet)
    
    guess = just_guess(Nm_values_real, idxs = c(83,114))
    
    
    varym_N2000 = list(one_m_N2000, two_m_N2000, three_m_N2000, five_m_N2000, nine_m_N2000, one_m_N2000_AF, two_m_N2000_AF, three_m_N2000_AF, five_m_N2000_AF, nine_m_N2000_AF, guess)
    saveRDS(varym_N2000, '../data/varym_N2000.RDS')
  }
  
  # processing
  varym_N2000_means = lapply(varym_N2000, apply, c(1,2,3), mean)
  
  {
    one_m_mean_N2000 = apply(one_m_N2000, c(1,2,3), mean)
    two_m_mean_N2000 = apply(two_m_N2000, c(1,2,3), mean)
    three_m_mean_N2000 = apply(three_m_N2000, c(1,2,3), mean)
    five_m_mean_N2000 = apply(five_m_N2000, c(1,2,3), mean)
    nine_m_mean_N2000 = apply(nine_m_N2000, c(1,2,3), mean)
    
    one_m_mean_N2000_AF = apply(one_m_N2000_AF, c(1,2,3), mean)
    two_m_mean_N2000_AF = apply(two_m_N2000_AF, c(1,2,3), mean)
    three_m_mean_N2000_AF = apply(three_m_N2000_AF, c(1,2,3), mean)
    five_m_mean_N2000_AF = apply(five_m_N2000_AF, c(1,2,3), mean)
    nine_m_mean_N2000_AF = apply(nine_m_N2000_AF, c(1,2,3), mean)
    
    guess_m_mean_N2000 = apply(guess, c(1,2,3), mean)
    
  }
  
  # make figures
  {
    pdf('../final_figures/n_sp_vs_m2', height = 7, width = 7)
    plot_compare_m_neat(list(one_m_mean_N2000,three_m_mean_N2000, nine_m_mean_N2000,one_m_mean_N2000_AF,three_m_mean_N2000AF, nine_m_mean_N2000_AF), legend = c('1 (LD)  ', '3 (LD)  ', '9 (LD)  ', '1 (AF)  ', '3 (AF)  ', '9 (AF)  '), ltys = c(1,2,3,1,2,3), cols = c(rep('darkslategrey',3), rep('red',3)))
    graphics.off()
  }
  

}




# make figure - m = 0.05, N = 2000
{
  N_plot = varyN_m0.05_means[c(1,3,6,7,11)]
  m_plot = varym_N2000_means[c(1,3,6,7,11)]
  
  legends = c('1 (Ind.)  ', '2 (Ind.)  ', '1 (Pool)  ', '2 (Pool)  ', 'Prior')
  ltys = c(1,5,1,5,1)
  cols = c(rep('darkslategrey',2), rep('red',2), '#4daf4a')
  
  plot_compare_N_m_neat(N_plot, m_plot, legends, ltys = ltys, cols = cols, fname = '../final_figures/n_sp_N_and_m.pdf', prior = T)
}


# make figure - m = 0.2, N = 10000
{
  N_plot = varyN_m0.20_means[c(1,2,6,7,11)]
  m_plot = varym_N10000_means[c(1,2,6,7,11)]
  
  legends = c('1 (Ind)  ', '2 (Ind)  ', '1 (Pool)  ', '2 (Pool)  ', 'Prior')
  ltys = c(1,5,1,5,1)
  cols = c(rep('darkslategrey',2), rep('red',2), '#4daf4a')
  
  plot_compare_N_m_neat(N_plot, m_plot, legends, ltys = ltys, cols = cols, fname = '../final_figures/n_sp_N_and_m_20_10000.pdf', prior = T)
}

# make figure - all poolseq
{
  N_plot = varyN_m0.05_means[c(6:11)]
  m_plot = varym_N2000_means[c(6:11)]
  
  legends = c('1 (AF)  ','2 (AF)  ', '3 (AF)  ','5 (AF)  ', '9 (AF)  ', 'Prior')
  ltys = c(1,1,2,3,4,1)
  colwheel = colorRampPalette(c('red', 'blue')); cols = colwheel(5)
  cols = c( c('dark red', 'red', 'red', 'red', 'red'), '#4daf4a')
  
  plot_compare_N_m_neat(N_plot, m_plot, legends, ltys = ltys, cols = cols, fname = '../final_figures/n_sp_N_and_m_many.pdf', prior = T)
}




# cross validation - small m
{
  selected_small = select_bounded(1000, Nm_values_sim, m_range = c(0,0.10), N_range = c(400,20000))

  # do ABC
  {
  cv1_ns_m_small = cross_validate(do_abc_PLS, SS_sim_one, selected_small, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2_ns_m_small = cross_validate(do_abc_PLS, SS_sim_two, selected_small, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv3_ns_m_small = cross_validate(do_abc_PLS, SS_sim_three, selected_small, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv5_ns_m_small = cross_validate(do_abc_PLS, SS_sim_five, selected_small, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv9_ns_m_small = cross_validate(do_abc_PLS, SS_sim_nine, selected_small, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  cv1_ns_AF_m_small = cross_validate(do_abc_PLS, SS_sim_one_AF, selected_small, n_comps =  min(nc, dim(SS_sim_one_AF)[2]), Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2_ns_AF_m_small = cross_validate(do_abc_PLS, SS_sim_two_AF, selected_small, n_comps =  min(nc, dim(SS_sim_two_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv3_ns_AF_m_small = cross_validate(do_abc_PLS, SS_sim_three_AF, selected_small, n_comps =  min(nc, dim(SS_sim_three_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv5_ns_AF_m_small = cross_validate(do_abc_PLS, SS_sim_five_AF, selected_small, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv9_ns_AF_m_small = cross_validate(do_abc_PLS, SS_sim_nine_AF, selected_small, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  spCVsmall = list(list(cv1_ns_m_small, cv2_ns_m_small, cv3_ns_m_small, cv5_ns_m_small, cv9_ns_m_small), list(cv1_ns_AF_m_small, cv2_ns_AF_m_small, cv5_ns_AF_m_small, cv5_ns_AF_m_small, cv9_ns_AF_m_small))
  saveRDS(spCVsmall, '../data/spCVsmall.RDS')
  spCVsmall = readRDS('../data/spCVsmall.RDS')
  
  }
  
  # processing
  {
  o_n_sp_m_small = array(NA, dim = c(0,2))
  for (cv in spCVsmall[[1]]) {RMS = RMS_err(cv); o_n_sp_m_small = rbind(o_n_sp_m_small,colMeans(RMS))}
  
  
  CIsLD_small = array(NA, dim = c(0,2,2))
  for (cv in spCVsmall[[1]]) {RMS = RMS_err(cv); CIsLD_small = abind(CIsLD_small, cbind(do_boot(RMS[,1]), do_boot(RMS[,2])), along = 1)}

  
  o_n_sp_AF_m_small = array(NA, dim = c(0,2))
  for (cv in spCVsmall[[2]]) {RMS = RMS_err(cv); o_n_sp_AF_m_small = rbind(o_n_sp_AF_m_small,colMeans(RMS))}
  
  
  CIsAF_small = array(NA, dim = c(0,2,2))
  for (cv in spCVsmall[[2]]) {RMS = RMS_err(cv); CIsAF_small = abind(CIsAF_small, cbind(do_boot(RMS[,1]), do_boot(RMS[,2])), along = 1)}

  }
  

}

# cross validation - large m
{
  selected_large = select_bounded(1000, Nm_values_sim, m_range = c(0.1,0.8), N_range = c(400,20000))
  
  # do ABC
  {
  cv1_ns_m_large = cross_validate(do_abc_PLS, SS_sim_one, selected_large, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2_ns_m_large = cross_validate(do_abc_PLS, SS_sim_two, selected_large, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv3_ns_m_large = cross_validate(do_abc_PLS, SS_sim_three, selected_large, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv5_ns_m_large = cross_validate(do_abc_PLS, SS_sim_five, selected_large, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv9_ns_m_large = cross_validate(do_abc_PLS, SS_sim_nine, selected_large, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  cv1_ns_AF_m_large = cross_validate(do_abc_PLS, SS_sim_one_AF, selected_large, n_comps =  min(nc, dim(SS_sim_one_AF)[2]), Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2_ns_AF_m_large = cross_validate(do_abc_PLS, SS_sim_two_AF, selected_large, n_comps =  min(nc, dim(SS_sim_two_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv3_ns_AF_m_large = cross_validate(do_abc_PLS, SS_sim_three_AF, selected_large, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv5_ns_AF_m_large = cross_validate(do_abc_PLS, SS_sim_five_AF, selected_large, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv9_ns_AF_m_large = cross_validate(do_abc_PLS, SS_sim_nine_AF, selected_large, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  spCVlarge = list(list(cv1_ns_m_large, cv2_ns_m_large, cv3_ns_m_large, cv5_ns_m_large, cv9_ns_m_large), list(cv1_ns_AF_m_large, cv2_ns_AF_m_large, cv3_ns_AF_m_large, cv5_ns_AF_m_large, cv9_ns_AF_m_large))
  saveRDS(spCVlarge, '../data/spCVlarge.RDS')
  spCVlarge = readRDS('../data/spCVlarge.RDS')
  }
  
  # processing
  {
    o_n_sp_m_large = array(NA, dim = c(0,2))
    for (cv in spCVlarge[[1]]) {RMS = RMS_err(cv); o_n_sp_m_large = rbind(o_n_sp_m_large,colMeans(RMS))}
    
    
    CIsLD_large = array(NA, dim = c(0,2,2))
    for (cv in spCVlarge[[1]]) {RMS = RMS_err(cv); CIsLD_large = abind(CIsLD_large, cbind(do_boot(RMS[,1]), do_boot(RMS[,2])), along = 1)}
    
    
    o_n_sp_AF_m_large = array(NA, dim = c(0,2))
    for (cv in spCVlarge[[2]]) {RMS = RMS_err(cv); o_n_sp_AF_m_large = rbind(o_n_sp_AF_m_large,colMeans(RMS))}
    
    
    CIsAF_large = array(NA, dim = c(0,2,2))
    for (cv in spCVlarge[[2]]) {RMS = RMS_err(cv); CIsAF_large = abind(CIsAF_large, cbind(do_boot(RMS[,1]), do_boot(RMS[,2])), along = 1)}
    
  }
  
}

# cross val plot
{
pdf(file = '../final_figures/CV_n_sp.pdf', height = 4.2, width = 7, pointsize = 7)
par(mfrow = c(1,2), mai = c(.7,.7,.4,.4), oma = c(3.2,0,0,0))
  #compare_CV(LD = o_n_sp, AF = o_n_sp_AF, LD_ci = CIsLD, AF_ci = CIsAF)
  
  compare_CV(LD = o_n_sp_m_small, AF = o_n_sp_AF_m_small, LD_ci = CIsLD_small, AF_ci = CIsAF_small,.4)
  fig_label('A', cex = 2)
  compare_CV(LD = o_n_sp_m_large, AF = o_n_sp_AF_m_large, LD_ci = CIsLD_large, AF_ci = CIsAF_large,.8)
  fig_label('B', cex = 2)
  
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('bottom',cex = .9, legend = c('N - LD      ', 'N - AF      ', 'm - LD      ', 'm - AF      ' ), pch = c(1,1,6, 6), col = c('black', 'red', 'black', 'red'), horiz = T, inset = rep(0.04,4))
graphics.off()
}





### not used

# cross validation
{
  
  # do ABC
  {
    cv1_ns = cross_validate(do_abc_PLS, SS_sim_one, selected, n_comps = nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv2_ns = cross_validate(do_abc_PLS, SS_sim_two, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv3_ns = cross_validate(do_abc_PLS, SS_sim_three, selected, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv5_ns = cross_validate(do_abc_PLS, SS_sim_five, selected, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv9_ns = cross_validate(do_abc_PLS, SS_sim_nine, selected, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    
    cv1_ns_AF = cross_validate(do_abc_PLS, SS_sim_one_AF, selected, n_comps = dim(SS_sim_one_AF)[2] , Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv2_ns_AF = cross_validate(do_abc_PLS, SS_sim_two_AF, selected, n_comps =  dim(SS_sim_two_AF)[2], Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv3_ns_AF = cross_validate(do_abc_PLS, SS_sim_three_AF, selected, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv5_ns_AF = cross_validate(do_abc_PLS, SS_sim_five_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv9_ns_AF = cross_validate(do_abc_PLS, SS_sim_nine_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    
    spCV = list(cv1_ns, cv2_ns, cv3_ns, cv5_ns, cv9_ns, cv1_ns_AF, cv2_ns_AF, cv3_ns_AF, cv5_ns_AF, cv9_ns_AF)
    saveRDS(spCV, '../data/spCV.RDS')
    spCV = readRDS(spCV, '../data/spCV.RDS')
    
  }
  
  
  
  
  # processing
  {
    o_n_sp = array(NA, dim = c(0,2))
    RMS_cv1 = RMS_err(cv1_ns); o_n_sp = rbind(o_n_sp,colMeans(RMS_cv1));
    RMS_cv2 = RMS_err(cv2_ns); o_n_sp = rbind(o_n_sp,colMeans(RMS_cv2));
    RMS_cv3 = RMS_err(cv3_ns); o_n_sp = rbind(o_n_sp,colMeans(RMS_cv3))
    RMS_cv5 = RMS_err(cv5_ns); o_n_sp = rbind(o_n_sp,colMeans(RMS_cv5))
    RMS_cv9 = RMS_err(cv9_ns); o_n_sp = rbind(o_n_sp,colMeans(RMS_cv9))
    
    CIsLD = array(NA, dim = c(5,2,2))
    CIsLD[1,,] =  cbind(do_boot(RMS_cv1[,1]), do_boot(RMS_cv1[,2]))
    CIsLD[2,,] =  cbind(do_boot(RMS_cv2[,1]), do_boot(RMS_cv2[,2]))
    CIsLD[3,,] =  cbind(do_boot(RMS_cv3[,1]), do_boot(RMS_cv3[,2]))
    CIsLD[4,,] =  cbind(do_boot(RMS_cv5[,1]), do_boot(RMS_cv5[,2]))
    CIsLD[5,,] =  cbind(do_boot(RMS_cv9[,1]), do_boot(RMS_cv9[,2]))
    
    
    
    o_n_sp_AF = array(NA, dim = c(0,2))
    RMS_cv1 = RMS_err(cv1_ns_AF); o_n_sp_AF = rbind(o_n_sp_AF,colMeans(RMS_cv1));
    RMS_cv2 = RMS_err(cv2_ns_AF); o_n_sp_AF = rbind(o_n_sp_AF,colMeans(RMS_cv2));
    RMS_cv3 = RMS_err(cv3_ns_AF); o_n_sp_AF = rbind(o_n_sp_AF,colMeans(RMS_cv3))
    RMS_cv5 = RMS_err(cv5_ns_AF); o_n_sp_AF = rbind(o_n_sp_AF,colMeans(RMS_cv5))
    RMS_cv9 = RMS_err(cv9_ns_AF); o_n_sp_AF = rbind(o_n_sp_AF,colMeans(RMS_cv9))
    
    CIsAF = array(NA, dim = c(5,2,2))
    CIsAF[1,,] =  cbind(do_boot(RMS_cv1[,1]), do_boot(RMS_cv1[,2]))
    CIsAF[2,,] =  cbind(do_boot(RMS_cv2[,1]), do_boot(RMS_cv2[,2]))
    CIsAF[3,,] =  cbind(do_boot(RMS_cv3[,1]), do_boot(RMS_cv3[,2]))
    CIsAF[4,,] =  cbind(do_boot(RMS_cv5[,1]), do_boot(RMS_cv5[,2]))
    CIsAF[5,,] =  cbind(do_boot(RMS_cv9[,1]), do_boot(RMS_cv9[,2]))
  }
  
  
}

# cross validation - tiny m
{
  selected_tiny = select_bounded(1000, Nm_values_sim, m_range = c(0,0.05), N_range = c(400,20000))
  
  # do ABC
  {
    cv1_ns_m_tiny = cross_validate(do_abc_PLS, SS_sim_one, selected_tiny, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv2_ns_m_tiny = cross_validate(do_abc_PLS, SS_sim_two, selected_tiny, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv3_ns_m_tiny = cross_validate(do_abc_PLS, SS_sim_three, selected_tiny, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv5_ns_m_tiny = cross_validate(do_abc_PLS, SS_sim_five, selected_tiny, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv9_ns_m_tiny = cross_validate(do_abc_PLS, SS_sim_nine, selected_tiny, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    
    cv1_ns_AF_m_tiny = cross_validate(do_abc_PLS, SS_sim_one_AF, selected_tiny, n_comps =  min(nc, dim(SS_sim_one_AF)[2]), Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv2_ns_AF_m_tiny = cross_validate(do_abc_PLS, SS_sim_two_AF, selected_tiny, n_comps =  min(nc, dim(SS_sim_two_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv3_ns_AF_m_tiny = cross_validate(do_abc_PLS, SS_sim_three_AF, selected_tiny, n_comps =  min(nc, dim(SS_sim_three_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv5_ns_AF_m_tiny = cross_validate(do_abc_PLS, SS_sim_five_AF, selected_tiny, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    cv9_ns_AF_m_tiny= cross_validate(do_abc_PLS, SS_sim_nine_AF, selected_tiny, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    
    spCVtiny = list(cv1_ns_m_tiny, cv2_ns_m_tiny, cv3_ns_m_tiny, cv5_ns_m_tiny, cv9_ns_m_tiny, cv1_ns_AF_m_tiny, cv2_ns_AF_m_tiny, cv3_ns_AF_m_tiny, cv5_ns_AF_m_tiny, cv9_ns_AF_m_tiny)
    saveRDS(spCVtiny, '../data/spCVtiny.RDS')
    
  }
  
  # processing
  {
    o_n_sp_m_tiny = array(NA, dim = c(0,2))
    RMS_cv1 = RMS_err(cv1_ns_m_tiny); o_n_sp_m_tiny = rbind(o_n_sp_m_tiny,colMeans(RMS_cv1));
    RMS_cv2 = RMS_err(cv2_ns_m_tiny); o_n_sp_m_tiny = rbind(o_n_sp_m_tiny,colMeans(RMS_cv2));
    RMS_cv3 = RMS_err(cv3_ns_m_tiny); o_n_sp_m_tiny = rbind(o_n_sp_m_tiny,colMeans(RMS_cv3))
    RMS_cv5 = RMS_err(cv5_ns_m_tiny); o_n_sp_m_tiny = rbind(o_n_sp_m_tiny,colMeans(RMS_cv5))
    RMS_cv9 = RMS_err(cv9_ns_m_tiny); o_n_sp_m_tiny = rbind(o_n_sp_m_tiny,colMeans(RMS_cv9))
    
    
    CIsLD_tiny = array(NA, dim = c(5,2,2))
    CIsLD_tiny[1,,] =  cbind(do_boot(RMS_cv1[,1]), do_boot(RMS_cv1[,2]))
    CIsLD_tiny[2,,] =  cbind(do_boot(RMS_cv2[,1]), do_boot(RMS_cv2[,2]))
    CIsLD_tiny[3,,] =  cbind(do_boot(RMS_cv3[,1]), do_boot(RMS_cv3[,2]))
    CIsLD_tiny[4,,] =  cbind(do_boot(RMS_cv5[,1]), do_boot(RMS_cv5[,2]))
    CIsLD_tiny[5,,] =  cbind(do_boot(RMS_cv9[,1]), do_boot(RMS_cv9[,2]))
    
    
    o_n_sp_AF_m_tiny = array(NA, dim = c(0,2))
    RMS_cv1 = RMS_err(cv1_ns_AF_m_tiny); o_n_sp_AF_m_tiny = rbind(o_n_sp_AF_m_tiny,colMeans(RMS_cv1));
    RMS_cv2 = RMS_err(cv2_ns_AF_m_tiny); o_n_sp_AF_m_tiny = rbind(o_n_sp_AF_m_tiny,colMeans(RMS_cv2));
    RMS_cv3 = RMS_err(cv3_ns_AF_m_tiny); o_n_sp_AF_m_tiny = rbind(o_n_sp_AF_m_tiny,colMeans(RMS_cv3))
    RMS_cv5 = RMS_err(cv5_ns_AF_m_tiny); o_n_sp_AF_m_tiny = rbind(o_n_sp_AF_m_tiny,colMeans(RMS_cv5))
    RMS_cv9 = RMS_err(cv9_ns_AF_m_tiny); o_n_sp_AF_m_tiny = rbind(o_n_sp_AF_m_tiny,colMeans(RMS_cv9))
    
    CIsAF_tiny = array(NA, dim = c(5,2,2))
    CIsAF_tiny[1,,] =  cbind(do_boot(RMS_cv1[,1]), do_boot(RMS_cv1[,2]))
    CIsAF_tiny[2,,] =  cbind(do_boot(RMS_cv2[,1]), do_boot(RMS_cv2[,2]))
    CIsAF_tiny[3,,] =  cbind(do_boot(RMS_cv3[,1]), do_boot(RMS_cv3[,2]))
    CIsAF_tiny[4,,] =  cbind(do_boot(RMS_cv5[,1]), do_boot(RMS_cv5[,2]))
    CIsAF_tiny[5,,] =  cbind(do_boot(RMS_cv9[,1]), do_boot(RMS_cv9[,2]))
  }
  
  # old
  {
    plot(c(1,2,3,5,9), o_n_sp_m[,1], type = 'p', pch = 16, ylab = 'RMSE (N)', xlab = 'Number of sampled subpopulations', xaxt = 'n', ylim = c(min(CIsA[,,1], CIsB[,,1]),max(CIsA[,,1], CIsB[,,1] ) ) ); axis(1, at = c(1,2,3,5,9), labels = c(1,2,3,5,9))
    for (i in 1:5){lines(rep(c(1,2,3,5,9)[i],2),CIsA[i,,1])}
    points(c(1,2,3,5,9), o_n_sp_AF_m[,1],pch = 16, col = 'red')
    for (i in 1:5){lines(rep(c(1,2,3,5,9)[i],2),CIsB[i,,1], col = 'red')}
    
    
    plot(c(1,2,3,5,9), o_n_sp_m[,2], type = 'p', pch = 16, ylab = 'RMSE (m)', xlab = 'Number of sampled subpopulations', xaxt = 'n', ylim = c(min(CIsA[,,2], CIsB[,,2]),max(CIsA[,,2], CIsB[,,2] ) ) ); axis(1, at = c(1,2,3,5,9), labels = c(1,2,3,5,9))
    for (i in 1:5){lines(rep(c(1,2,3,5,9)[i],2),CIsA[i,,2])}
    points(c(1,2,3,5,9), o_n_sp_AF_m[,2],pch = 16, col = 'red')
    for (i in 1:5){lines(rep(c(1,2,3,5,9)[i],2),CIsB[i,,2], col = 'red')}
  }
  
}


# cross validation - small N
{
  selected = select_bounded(50, Nm_values_sim, m_range = c(0,0.8), N_range = c(400,2000))
  
  cv1_ns_N = cross_validate(do_abc_PLS, SS_sim_one, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv2_ns_N = cross_validate(do_abc_PLS, SS_sim_two, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv3_ns_N = cross_validate(do_abc_PLS, SS_sim_three, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv5_ns_N = cross_validate(do_abc_PLS, SS_sim_five, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv9_ns_N = cross_validate(do_abc_PLS, SS_sim_nine, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  
  cv1_ns_AF_N = cross_validate(do_abc_PLS, SS_sim_one_AF, selected, n_comps =  nc, Nm_values_sim,method = method, tol = tol)
  cv2_ns_AF_N = cross_validate(do_abc_PLS, SS_sim_two_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv3_ns_AF_N = cross_validate(do_abc_PLS, SS_sim_three_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv5_ns_AF_N = cross_validate(do_abc_PLS, SS_sim_five_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv9_ns_AF_N = cross_validate(do_abc_PLS, SS_sim_nine_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  #cv3_one_far = cross_validate(do_abc_PLS, SS_sim_three_far, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1)
  #cv4 = cross_validate(do_abc_PLS, SS_sim_four, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1)
  #cv4_one_far = cross_validate(do_abc_PLS, SS_sim_four_one_far, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1)
  #cv5 = cross_validate(do_abc_PLS, SS_sim_five, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1)
  #cv5_one_far = cross_validate(do_abc_PLS, SS_sim_five_one_far, selected, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .1)
  
  o_n_sp_N = array(NA, dim = c(0,2))
  RMS_cv1 = RMS_err(cv1_ns_N); o_n_sp_N = rbind(o_n_sp_N,colMeans(RMS_cv1));
  RMS_cv2 = RMS_err(cv2_ns_N); o_n_sp_N = rbind(o_n_sp_N,colMeans(RMS_cv2));
  RMS_cv3 = RMS_err(cv3_ns_N); o_n_sp_N = rbind(o_n_sp_N,colMeans(RMS_cv3))
  RMS_cv5 = RMS_err(cv5_ns_N); o_n_sp_N = rbind(o_n_sp_N,colMeans(RMS_cv5))
  RMS_cv9 = RMS_err(cv9_ns_N); o_n_sp_N = rbind(o_n_sp_N,colMeans(RMS_cv9))
  
  
  o_n_sp_AF_N = array(NA, dim = c(0,2))
  RMS_cv1 = RMS_err(cv1_ns_AF_N); o_n_sp_AF_N = rbind(o_n_sp_AF_N,colMeans(RMS_cv1));
  RMS_cv2 = RMS_err(cv2_ns_AF_N); o_n_sp_AF_N = rbind(o_n_sp_AF_N,colMeans(RMS_cv2));
  RMS_cv3 = RMS_err(cv3_ns_AF_N); o_n_sp_AF_N = rbind(o_n_sp_AF_N,colMeans(RMS_cv3))
  RMS_cv5 = RMS_err(cv5_ns_AF_N); o_n_sp_AF_N = rbind(o_n_sp_AF_N,colMeans(RMS_cv5))
  RMS_cv9 = RMS_err(cv9_ns_AF_N); o_n_sp_AF_N = rbind(o_n_sp_AF_N,colMeans(RMS_cv9))
  
  
  plot(c(1,2,3,5,9), o_n_sp_N[,1], type = 'l', ylab = 'RMSE (N)', xlab = 'Number of sampled subpopulations', xaxt = 'n', ylim = c(min(o_n_sp_N[,1],o_n_sp_AF_N[,1]),max(o_n_sp_N[,1],o_n_sp_AF_N[,1] ) ) ); axis(1, at = c(1,2,3,5,9), labels = c(1,2,3,5,9))
  lines(c(1,2,3,5,9), o_n_sp_AF_N[,1], col = 'red')
  plot(c(1,2,3,5,9), o_n_sp_N[,2], type = 'l', ylab = 'RMSE (m)', xlab = 'Number of sampled subpopulations', xaxt = 'n', ylim = c(min(o_n_sp_N[,2],o_n_sp_AF_N[,2]),max(o_n_sp_N[,2],o_n_sp_AF_N[,2] ) ) ); axis(1, at = c(1,2,3,5,9), labels = c(1,2,3,5,9))
  lines(c(1,2,3,5,9), o_n_sp_AF_N[,2], col = 'red')
  
  #legend('topright', col = c('black'), lty = c(1,2), legend = c('adjacent popuations only', 'including distant population'))
  
}
}

#nearby subpopulations
{
  SS_five_near = make_SS(sim_out, real_out_list, S_select = S5, gens_select = gens_select, subpops_select_sim = c(40,41,42,50,32), subpops_select_real = c(60,61,62,71,51), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_five_near = SS_five_near$sim[,]; SS_real_five_near = SS_five_near$real
  
  SS_five_AF_near = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S5*30, 30)), gens_select = gens_select, subpops_select_sim = c(40,41,42,50,32), subpops_select_real = c(60,61,62,71,51), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_five_AF_near = SS_five_AF_near$sim[,]; SS_real_five_AF_near = SS_five_AF_near$real
  
  
  SS_five_one = make_SS(sim_out, real_out_list, S_select = S5, gens_select = gens_select, subpops_select_sim = c(40,41,42,49,50), subpops_select_real = c(60,61,62,71,72), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_five_one = SS_five_one$sim[,]; SS_real_five_one = SS_five_one$real
  
  SS_five_AF_one = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S5*30, 30)), gens_select = gens_select, subpops_select_sim = c(40,41,42,49,50), subpops_select_real = c(60,61,62,71,72), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_five_AF_one = SS_five_AF_one$sim[,]; SS_real_five_AF_one = SS_five_AF_one$real
  
  
  SS_five_two = make_SS(sim_out, real_out_list, S_select = S5, gens_select = gens_select, subpops_select_sim = c(40,41,42,33,51), subpops_select_real = c(60,61,62,52,72), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_five_two = SS_five_two$sim[,]; SS_real_five_two = SS_five_two$real
  
  SS_five_AF_two = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S5*30, 30)), gens_select = gens_select, subpops_select_sim = c(40,41,42,33,51), subpops_select_real = c(60,61,62,52,72), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_five_AF_two = SS_five_AF_two$sim[,]; SS_real_five_AF_two = SS_five_AF_two$real
}

##################################
##################################
######## Across time #############
##################################
##################################
{
# for a given across-space sampling strategy, how does estiamtion behave when varying across-time strategy
## across N, m

## want to look for 'better at smaller migration rate' type information - plain best strategy comes later

S_total = 1000
S_per_pop = 1000 / 3
gens_1 = c(0); S1 = round(S_per_pop)
gens_2 = c(0,20); S2 = round(S_per_pop / 2)
gens_3 = c(0,10, 20); S3 = round(S_per_pop / 3 )
gens_4 = c(0,5,10,20); S4 = round(S_per_pop / 4 )

gens_2.5 = c(0,5); gens_2.10 = c(0,10); gens_2.5 = c(0,5)

# make SS
{
SS_1 = make_SS(sim_out, real_out_list, S_select = S1, gens_select = gens_1, subpops_select_sim = c(40:42), subpops_select_real = c(60:62), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_1 = SS_1$sim; SS_real_1 = SS_1$real

SS_2 = make_SS(sim_out, real_out_list, S_select = S2, gens_select = gens_2, subpops_select_sim = c(40:42), subpops_select_real = c(60:62), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_2 = SS_2$sim; SS_real_2 = SS_2$real

SS_2.10 = make_SS(sim_out, real_out_list, S_select = S2, gens_select = gens_2.10, subpops_select_sim = c(40:42), subpops_select_real = c(60:62), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_2.10 = SS_2.10$sim; SS_real_2.10 = SS_2.10$real

SS_2.5 = make_SS(sim_out, real_out_list, S_select = S2, gens_select = gens_2.5, subpops_select_sim = c(40:42), subpops_select_real = c(60:62), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_2.5 = SS_2.5$sim; SS_real_2.5 = SS_2.5$real

SS_3 = make_SS(sim_out, real_out_list, S_select = S3, gens_select = gens_3, subpops_select_sim = c(40:42), subpops_select_real = c(60:62), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_3 = SS_3$sim; SS_real_3 = SS_3$real

SS_4 = make_SS(sim_out, real_out_list, S_select = S4, gens_select = gens_4, subpops_select_sim = c(40:42), subpops_select_real = c(60:62), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_4 = SS_4$sim; SS_real_4 = SS_4$real


SS_1_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S1*30, 30)), gens_select = gens_1, subpops_select_sim = c(40:42), subpops_select_real = c(60:62), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_1_AF = SS_1_AF$sim; SS_real_1_AF = SS_1_AF$real

SS_2_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S2*30, 30)), gens_select = gens_2, subpops_select_sim = c(40:42), subpops_select_real = c(60:62), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_2_AF = SS_2_AF$sim; SS_real_2_AF = SS_2_AF$real

SS_2.10_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S2*30, 30)), gens_select = gens_2.10, subpops_select_sim = c(40:42), subpops_select_real = c(60:62), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_2.10_AF = SS_2.10_AF$sim; SS_real_2.10_AF = SS_2.10_AF$real

SS_2.5_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S2*30, 30)), gens_select = gens_2.5, subpops_select_sim = c(40:42), subpops_select_real = c(60:62), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_2.5_AF = SS_2.5_AF$sim; SS_real_2.5_AF = SS_2.5_AF$real

SS_3_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S3*30, 30)), gens_select = gens_3, subpops_select_sim = c(40:42), subpops_select_real = c(60:62), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_3_AF = SS_3_AF$sim; SS_real_3_AF = SS_3_AF$real

SS_4_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S4*30, 30)), gens_select = gens_4, subpops_select_sim = c(40:42), subpops_select_real = c(60:62), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_4_AF = SS_4_AF$sim; SS_real_4_AF = SS_4_AF$real
}

# across N
{
  # do ABC
  {
  one_N_time = test_specific(do_abc_PLS, SS_sim_1, SS_real_1, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27),method = method, tol = tol)
  two_N_time = test_specific(do_abc_PLS, SS_sim_2, SS_real_2, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27),method = method, tol = tol)
  three_N_time = test_specific(do_abc_PLS, SS_sim_3, SS_real_3, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27), method = method, tol = tol)
  four_N_time = test_specific(do_abc_PLS, SS_sim_4, SS_real_4, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27), method = method, tol = tol)
  
  two_N.10_time = test_specific(do_abc_PLS, SS_sim_2.10, SS_real_2.10, n_comps = min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(1,27), method = method, tol = tol)
  two_N.5_time = test_specific(do_abc_PLS, SS_sim_2.5, SS_real_2.5, n_comps = min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(1,27), method = method, tol = tol)
  
  
  #one_N_time_AF = test_specific(do_abc_PLS, SS_sim_1_AF, SS_real_1_AF, n_comps = min(nc, dim(SS_sim_1_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(1,27), method = method, tol = tol)
  two_N_time_AF = test_specific(do_abc_PLS, SS_sim_2_AF, SS_real_2_AF, n_comps = min(nc, dim(SS_sim_2_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(1,27), method = method, tol = tol)
  three_N_time_AF = test_specific(do_abc_PLS, SS_sim_3_AF, SS_real_3_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27), method = method, tol = tol)
  four_N_time_AF = test_specific(do_abc_PLS, SS_sim_4_AF, SS_real_4_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(1,27), method = method, tol = tol)
  
  two_N.10_time_AF = test_specific(do_abc_PLS, SS_sim_2.10_AF, SS_real_2.10, n_comps = min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(1,27), method = method, tol = tol)
  two_N.5_time_AF = test_specific(do_abc_PLS, SS_sim_2.5_AF, SS_real_2.5, n_comps = min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(1,27), method = method, tol = tol)
  
  
  varyN_m0.05_tp = list(one_N_time, two_N_time, three_N_time, four_N_time, two_N.10_time, two_N.5_time, two_N_time_AF, three_N_time_AF, four_N_time_AF, two_N.10_time_AF,two_N.5_time_AF)
  saveRDS(varyN_m0.05_tp, '../data/varyN_m0.05_tp.RDS')
  }
  
  
  # processing
  {
  one_N_mean_time = apply(one_N_time, c(1,2,3), mean)
  two_N_mean_time = apply(two_N_time, c(1,2,3), mean)
  three_N_mean_time = apply(three_N_time, c(1,2,3), mean)
  four_N_mean_time = apply(four_N_time, c(1,2,3), mean)
  
  two_N_mean.10_time = apply(two_N.10_time, c(1,2,3), mean)
  two_N_mean.5_time = apply(two_N.5_time, c(1,2,3), mean)
  
  
  #one_N_mean_time_AF = apply(one_N_time, c(1,2,3), mean)
  two_N_mean_time_AF = apply(two_N_time_AF, c(1,2,3), mean)
  three_N_mean_time_AF = apply(three_N_time_AF, c(1,2,3), mean)
  four_N_mean_time_AF = apply(four_N_time_AF, c(1,2,3), mean)
  
  two_N_mean.10_time_AF = apply(two_N.10_time_AF, c(1,2,3), trim_mean, 0)
  two_N_mean.5_time_AF = apply(two_N.5_time_AF, c(1,2,3), trim_mean, 0)
  }
  
  
  # plotting
  {
    plot_compare_N_neat(...)
  }
  
  # old
  {
  plot_compare_N(list(one_N_mean_time,two_N_mean_time, three_N_mean_time, four_N_mean_time), legend = c('1 timepoint', '2 timepoints', '3 timepoints', '4 timepoints'), cols = c('black', 'black'), ltys = 1:100)
  plot_compare_N(list(two_m_mean_time, three_m_mean_time, four_m_mean_time), legend = c( '2 timepoints', '3 timepoints', '4 timepoints'), cols = c('black', 'black'), ltys = 1:100)
  
  plot_compare_m(list(two_m_mean_time,two_m_mean.10_time, two_m_mean.20_time), legend = c('5 sep', '10 sep', '20 sep'))
  }
}

# across N, m larger
{
  # do ABC
  {
    one_N_time_mig = test_specific(do_abc_PLS, SS_sim_1, SS_real_1, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76),method = method, tol = tol)
    two_N_time_mig = test_specific(do_abc_PLS, SS_sim_2, SS_real_2, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76),method = method, tol = tol)
    three_N_time_mig = test_specific(do_abc_PLS, SS_sim_3, SS_real_3, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76), method = method, tol = tol)
    four_N_time_mig = test_specific(do_abc_PLS, SS_sim_4, SS_real_4, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76), method = method, tol = tol)
    
    two_N.10_time_mig = test_specific(do_abc_PLS, SS_sim_2.10, SS_real_2.10, n_comps = min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(63,76), method = method, tol = tol)
    two_N.5_time_mig = test_specific(do_abc_PLS, SS_sim_2.5, SS_real_2.5, n_comps = min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(63,76), method = method, tol = tol)
    
    
    #one_N_time_AF = test_specific(do_abc_PLS, SS_sim_1_AF, SS_real_1_AF, n_comps = min(nc, dim(SS_sim_1_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(1,27), method = method, tol = tol)
    two_N_time_AF_mig = test_specific(do_abc_PLS, SS_sim_2_AF, SS_real_2_AF, n_comps = min(nc, dim(SS_sim_2_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(63,76), method = method, tol = tol)
    three_N_time_AF_mig = test_specific(do_abc_PLS, SS_sim_3_AF, SS_real_3_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76), method = method, tol = tol)
    four_N_time_AF_mig = test_specific(do_abc_PLS, SS_sim_4_AF, SS_real_4_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(63,76), method = method, tol = tol)
    
    two_N.10_time_AF_mig = test_specific(do_abc_PLS, SS_sim_2.10_AF, SS_real_2.10, n_comps = min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(63,76), method = method, tol = tol)
    two_N.5_time_AF_mig = test_specific(do_abc_PLS, SS_sim_2.5_AF, SS_real_2.5, n_comps = min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(63,76), method = method, tol = tol)
    
    
    varyN_m0.20_tp = list(one_N_time_mig, two_N_time_mig, three_N_time_mig, four_N_time_mig, two_N.10_time_mig, two_N.5_time_mig, two_N_time_AF_mig, three_N_time_AF_mig, four_N_time_AF_mig, two_N.10_time_AF_mig,two_N.5_time_AF_mig)
    saveRDS(varyN_m0.20_tp, '../data/varyN_m0.20_tp.RDS')
  }
  
  
  # processing
  {
    one_N_mean_time_mig = apply(one_N_time_mig, c(1,2,3), mean)
    two_N_mean_time_mig = apply(two_N_time_mig, c(1,2,3), mean)
    three_N_mean_time_mig = apply(three_N_time_mig, c(1,2,3), mean)
    four_N_mean_time_mig = apply(four_N_time_mig, c(1,2,3), mean)
    
    two_N_mean.10_time_mig = apply(two_N.10_time_mig, c(1,2,3), mean)
    two_N_mean.5_time_mig = apply(two_N.5_time_mig, c(1,2,3), mean)
    
    
    #one_N_mean_time_AF = apply(one_N_time, c(1,2,3), mean)
    two_N_mean_time_AF_mig = apply(two_N_time_AF_mig, c(1,2,3), mean)
    three_N_mean_time_AF_mig = apply(three_N_time_AF_mig, c(1,2,3), mean)
    four_N_mean_time_AF_mig = apply(four_N_time_AF_mig, c(1,2,3), mean)
    
    two_N_mean.10_time_AF_mig = apply(two_N.10_time_AF_mig, c(1,2,3), trim_mean, 0)
    two_N_mean.5_time_AF_mig = apply(two_N.5_time_AF_mig, c(1,2,3), trim_mean, 0)
  }
  
  
  # plotting
  {
    plot_compare_N_neat(...)
  }
  
  # old
  {
    plot_compare_N(list(one_N_mean_time,two_N_mean_time, three_N_mean_time, four_N_mean_time), legend = c('1 timepoint', '2 timepoints', '3 timepoints', '4 timepoints'), cols = c('black', 'black'), ltys = 1:100)
    plot_compare_N(list(two_m_mean_time, three_m_mean_time, four_m_mean_time), legend = c( '2 timepoints', '3 timepoints', '4 timepoints'), cols = c('black', 'black'), ltys = 1:100)
    
    plot_compare_m(list(two_m_mean_time,two_m_mean.10_time, two_m_mean.20_time), legend = c('5 sep', '10 sep', '20 sep'))
  }
}

# across m
{
  
  # do_ABC
  {
  one_m_time = test_specific(do_abc_PLS, SS_sim_1, SS_real_1, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
  two_m_time = test_specific(do_abc_PLS, SS_sim_2, SS_real_2, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
  three_m_time = test_specific(do_abc_PLS, SS_sim_3, SS_real_3, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
  four_m_time = test_specific(do_abc_PLS, SS_sim_4, SS_real_4, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
  
  two_m.10_time = test_specific(do_abc_PLS, SS_sim_2.10, SS_real_2.10, n_comps = min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
  two_m.5_time = test_specific(do_abc_PLS, SS_sim_2.5, SS_real_2.5, n_comps = min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
  
  #one_m_time_AF = test_specific(do_abc_PLS, SS_sim_1_AF, SS_real_1_AF, n_comps = min(nc, dim(SS_sim_1_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(28,46), method = method, tol = tol)
  two_m_time_AF = test_specific(do_abc_PLS, SS_sim_2_AF, SS_real_2_AF, n_comps = min(nc, dim(SS_sim_2_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
  three_m_time_AF = test_specific(do_abc_PLS, SS_sim_3_AF, SS_real_3_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
  four_m_time_AF = test_specific(do_abc_PLS, SS_sim_4_AF, SS_real_4_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
  
  two_m.10_time_AF = test_specific(do_abc_PLS, SS_sim_2.10_AF, SS_real_2.10, n_comps = min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
  two_m.5_time_AF = test_specific(do_abc_PLS, SS_sim_2.5_AF, SS_real_2.5, n_comps = min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
  
  varym_tp = list(one_m_time, two_m_time, three_m_time, four_m_time, two_m.10_time, two_m.5_time, two_m_time_AF, three_m_time_AF, four_m_time_AF, two_m.10_time_AF,two_m.5_time_AF)
  saveRDS(varym_tp, '../data/varym_tp.RDS')
  
  }
  
  # processing
  {
  one_m_mean_time = apply(one_m_time, c(1,2,3), mean)
  two_m_mean_time = apply(two_m_time, c(1,2,3), mean)
  three_m_mean_time = apply(three_m_time, c(1,2,3), mean)
  four_m_mean_time = apply(four_m_time, c(1,2,3), mean)
  
  two_m_mean.10_time = apply(two_m.10_time, c(1,2,3), mean)
  two_m_mean.5_time = apply(two_m.5_time, c(1,2,3), mean)
  
  #one_m_mean_time_AF = apply(one_m_time_AF, c(1,2,3), mean)
  two_m_mean_time_AF = apply(two_m_time_AF, c(1,2,3), mean)
  three_m_mean_time_AF = apply(three_m_time_AF, c(1,2,3), mean)
  four_m_mean_time_AF = apply(four_m_time_AF, c(1,2,3), mean)
  
  two_m_mean.10_time_AF = apply(two_m.10_time_AF, c(1,2,3), mean)
  two_m_mean.5_time_AF = apply(two_m.5_time_AF, c(1,2,3), mean)
  }
  
  # plotting
  {
    plot_compare_m_neat(...)
  }
  
  # old
  {
  plot_compare_m(list(one_m_mean_time,two_m_mean_time, three_m_mean_time, four_m_mean_time), legend = c('1 timepoint', '2 timepoints', '3 timepoints', '4 timepoints'), cols = c('black', 'black'), ltys = 1:100)
  plot_compare_m(list(two_m_mean_time, three_m_mean_time, four_m_mean_time), legend = c( '2 timepoints', '3 timepoints', '4 timepoints'), cols = c('black', 'black'), ltys = 1:100)
  
  plot_compare_m(list(two_m_mean_time,two_m_mean.10_time, two_m_mean.20_time), legend = c('5 sep', '10 sep', '20 sep'))
  }
}

# across m, N = 2000
{
  
  # do_ABC
  {
    one_m_time_N2000 = test_specific(do_abc_PLS, SS_sim_1, SS_real_1, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
    two_m_time_N2000 = test_specific(do_abc_PLS, SS_sim_2, SS_real_2, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
    three_m_time_N2000 = test_specific(do_abc_PLS, SS_sim_3, SS_real_3, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
    four_m_time_N2000 = test_specific(do_abc_PLS, SS_sim_4, SS_real_4, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
    
    two_m.10_time_N2000 = test_specific(do_abc_PLS, SS_sim_2.10, SS_real_2.10, n_comps = min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
    two_m.5_time_N2000 = test_specific(do_abc_PLS, SS_sim_2.5, SS_real_2.5, n_comps = min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
    
    #one_m_time_AF = test_specific(do_abc_PLS, SS_sim_1_AF, SS_real_1_AF, n_comps = min(nc, dim(SS_sim_1_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(28,46), method = method, tol = tol)
    two_m_time_AF_N2000 = test_specific(do_abc_PLS, SS_sim_2_AF, SS_real_2_AF, n_comps = min(nc, dim(SS_sim_2_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
    three_m_time_AF_N2000 = test_specific(do_abc_PLS, SS_sim_3_AF, SS_real_3_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
    four_m_time_AF_N2000 = test_specific(do_abc_PLS, SS_sim_4_AF, SS_real_4_AF, n_comps = nc, Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
    
    two_m.10_time_AF_N2000 = test_specific(do_abc_PLS, SS_sim_2.10_AF, SS_real_2.10, n_comps = min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
    two_m.5_time_AF_N2000 = test_specific(do_abc_PLS, SS_sim_2.5_AF, SS_real_2.5, n_comps = min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(28,62), method = method, tol = tol)
    
    varym_tp_N2000 = list(one_m_time_N2000, two_m_time_N2000, three_m_time_N2000, four_m_time_N2000, two_m.10_time_N2000, two_m.5_time_N2000, two_m_time_AF_N2000, three_m_time_AF_N2000, four_m_time_AF_N2000, two_m.10_time_AF_N2000,two_m.5_time_AF_N2000)
    saveRDS(varym_tp_N2000, '../data/varym_tp_N2000.RDS')
    
  }
  
  # processing
  {
    one_m_mean_time_N2000 = apply(one_m_time_N2000, c(1,2,3), mean)
    two_m_mean_time_N2000 = apply(two_m_time_N2000, c(1,2,3), mean)
    three_m_mean_time_N2000 = apply(three_m_time_N2000, c(1,2,3), mean)
    four_m_mean_time_N2000 = apply(four_m_time_N2000, c(1,2,3), mean)
    
    two_m_mean.10_time_N2000 = apply(two_m.10_time_N2000, c(1,2,3), mean)
    two_m_mean.5_time_N2000 = apply(two_m.5_time_N2000, c(1,2,3), mean)
    
    #one_m_mean_time_AF = apply(one_m_time_AF, c(1,2,3), mean)
    two_m_mean_time_AF_N2000 = apply(two_m_time_AF_N2000, c(1,2,3), mean)
    three_m_mean_time_AF_N2000 = apply(three_m_time_AF_N2000, c(1,2,3), mean)
    four_m_mean_time_AF_N2000 = apply(four_m_time_AF_N2000, c(1,2,3), mean)
    
    two_m_mean.10_time_AF_N2000 = apply(two_m.10_time_AF_N2000, c(1,2,3), mean)
    two_m_mean.5_time_AF_N2000 = apply(two_m.5_time_AF_N2000, c(1,2,3), mean)
  }
  
  # plotting
  {
    plot_compare_m_neat(...)
  }
  
  # old
  {
    plot_compare_m(list(one_m_mean_time,two_m_mean_time, three_m_mean_time, four_m_mean_time), legend = c('1 timepoint', '2 timepoints', '3 timepoints', '4 timepoints'), cols = c('black', 'black'), ltys = 1:100)
    plot_compare_m(list(two_m_mean_time, three_m_mean_time, four_m_mean_time), legend = c( '2 timepoints', '3 timepoints', '4 timepoints'), cols = c('black', 'black'), ltys = 1:100)
    
    plot_compare_m(list(two_m_mean_time,two_m_mean.10_time, two_m_mean.20_time), legend = c('5 sep', '10 sep', '20 sep'))
  }
}

# combined plot
{
m_plot = list(one_m_mean_time,two_m_mean_time, three_m_mean_time, four_m_mean_time, two_m_mean_time_AF, three_m_mean_time_AF, four_m_mean_time_AF)
N_plot = list(one_N_mean_time,two_N_mean_time, three_N_mean_time, four_N_mean_time, two_N_mean_time_AF, three_N_mean_time_AF, four_N_mean_time_AF)

legends = c('1 (LD)  ', '2 (LD)  ', '3 (LD)  ', '4 (LD)  ', '2 (AF)  ', '3 (AF)  ', '4 (AF)  ')
ltys = c(1,2,3,4,1,2,3)
cols = c(rep('darkslategrey',4), rep('red',3))

plot_compare_N_m_neat(N_plot, m_plot, legends, ltys = ltys, cols = cols, fname = '../final_figures/n_tp_N_and_m.pdf')
}

# larger m plot
{
  pdf('../final_figures/n_tp_vs_N_m020.pdf', height = 7, width = 7, pointsize = 7)
  plot_compare_N_neat(list(one_N_mean_time_mig, two_N_mean_time_mig, three_N_mean_time_mig, four_N_mean_time_mig, two_N_mean_time_AF_mig,three_N_mean_time_AF_mig, four_N_mean_time_AF_mig), legend = c('1 (LD)  ', '2 (LD)  ', '3 (LD)  ','4 (LD)  ', '2 (AF)  ', '3 (AF)  ', '4 (AF)  '), ltys = c(1,2,3,4,1,2,3), cols = c(rep('darkslategrey',4), rep('red',3)))
  graphics.off()
}

# cross validation
{
  selected
  
  # do ABC
  {
  cv1_time = cross_validate(do_abc_PLS, SS_sim_1, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2_time = cross_validate(do_abc_PLS, SS_sim_2, selected, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv3_time = cross_validate(do_abc_PLS, SS_sim_3, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv4_time = cross_validate(do_abc_PLS, SS_sim_4, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  cv2.10_time = cross_validate(do_abc_PLS, SS_sim_2.10, selected, n_comps =  min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2.5_time = cross_validate(do_abc_PLS, SS_sim_2.5, selected, n_comps =  min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  #cv1_time_AF = cross_validate(do_abc_PLS, SS_sim_1_AF, selected, n_comps =  min(nc, dim(SS_sim_1_AF)[2]), Nm_values_sim, method = method, tol = tol)
  cv2_time_AF = cross_validate(do_abc_PLS, SS_sim_2_AF, selected, n_comps =  min(nc, dim(SS_sim_2_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv3_time_AF = cross_validate(do_abc_PLS, SS_sim_3_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv4_time_AF = cross_validate(do_abc_PLS, SS_sim_4_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  cv2.10_time_AF = cross_validate(do_abc_PLS, SS_sim_2.10_AF, selected, n_comps =  min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2.5_time_AF = cross_validate(do_abc_PLS, SS_sim_2.5_AF, selected, n_comps =  min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  
  cv_tp = list(cv1_time, cv2_time, cv3_time, cv4_time, cv2.10_time, cv2.5_time, cv2_time_AF, cv3_time_AF, cv4_time_AF, cv2.10_time_AF,cv2.5_time_AF)
  saveRDS(cv_tp, '../data/cv_tp.RDS')
  
  }


  # processing
  {
  o_n_tp = array(NA, dim = c(0,2))
  RMS_cv1 = RMS_err(cv1_time); o_n_tp = rbind(o_n_tp,colMeans(RMS_cv1));
  RMS_cv2 = RMS_err(cv2_time); o_n_tp = rbind(o_n_tp,colMeans(RMS_cv2));
  RMS_cv3 = RMS_err(cv3_time); o_n_tp = rbind(o_n_tp,colMeans(RMS_cv3))
  RMS_cv4 = RMS_err(cv4_time); o_n_tp = rbind(o_n_tp,colMeans(RMS_cv4))
  
  
  CIsLD_tp = array(NA, dim = c(4,2,2))
  CIsLD_tp[1,,] =  cbind(do_boot(RMS_cv1[,1]), do_boot(RMS_cv1[,2]))
  CIsLD_tp[2,,] =  cbind(do_boot(RMS_cv2[,1]), do_boot(RMS_cv2[,2]))
  CIsLD_tp[3,,] =  cbind(do_boot(RMS_cv3[,1]), do_boot(RMS_cv3[,2]))
  CIsLD_tp[4,,] =  cbind(do_boot(RMS_cv4[,1]), do_boot(RMS_cv4[,2]))

  
  o_n_tp_AF = array(NA, dim = c(0,2))
  #RMS_cv1 = RMS_err(cv1_time_AF); o_n_tp_AF = rbind(o_n_tp_AF,colMeans(RMS_cv1));
  RMS_cv2 = RMS_err(cv2_time_AF); o_n_tp_AF = rbind(o_n_tp_AF,colMeans(RMS_cv2));
  RMS_cv3 = RMS_err(cv3_time_AF); o_n_tp_AF = rbind(o_n_tp_AF,colMeans(RMS_cv3))
  RMS_cv4 = RMS_err(cv4_time_AF); o_n_tp_AF = rbind(o_n_tp_AF,colMeans(RMS_cv4))
  
  
  CIsAF_tp = array(NA, dim = c(3,2,2))
  #CIsAF_tp[1,,] =  cbind(do_boot(RMS_cv1[,1]), do_boot(RMS_cv2[,2]))
  CIsAF_tp[1,,] =  cbind(do_boot(RMS_cv2[,1]), do_boot(RMS_cv2[,2]))
  CIsAF_tp[2,,] =  cbind(do_boot(RMS_cv3[,1]), do_boot(RMS_cv3[,2]))
  CIsAF_tp[3,,] =  cbind(do_boot(RMS_cv4[,1]), do_boot(RMS_cv4[,2]))
  }
  
  # old
  {
  
  plot(c(1,2,3,4), o_n_tp[,1], type = 'p', pch = 16, ylab = 'RMSE (N)', xlab = 'Number of sampled subpopulations', xaxt = 'n', ylim = c(min(CIsA[,,1], CIsB[,,1]),max(CIsA[,,1], CIsB[,,1] ) ) ); axis(1, at = c(1,2,3,5,9), labels = c(1,2,3,5,9))
  for (i in 1:4){lines(rep(c(1,2,3,5,9)[i],2),CIsA[i,,1])}
  points(c(2,3,4), o_n_tp_AF[,1],pch = 16, col = 'red')
  for (i in 1:3){lines(rep(c(2,3,5,9)[i],2),CIsB[i,,1], col = 'red')}
  
  
  plot(c(1,2,3,4), o_n_tp[,2], type = 'p', pch = 16, ylab = 'RMSE (m)', xlab = 'Number of sampled subpopulations', xaxt = 'n', ylim = c(min(CIsA[,,2], CIsB[,,2]),max(CIsA[,,2], CIsB[,,2] ) ) ); axis(1, at = c(1,2,3,5,9), labels = c(1,2,3,5,9))
  for (i in 1:4){lines(rep(c(1,2,3,5,9)[i],2),CIsA[i,,2])}
  points(c(2,3,4), o_n_tp_AF[,2],pch = 16, col = 'red')
  for (i in 1:3){lines(rep(c(2,3,5,9)[i],2),CIsB[i,,2], col = 'red')}
  


  plot(c(2,3,4), o_time[,1], type = 'l', ylab = 'RMSE (N)', xlab = 'Number of sampled timepoints', xaxt = 'n', ylim = c(min(o_time[,1]),max(o_time[,1]) ) ) ; axis(1, at = c(1,2,3,4,5), labels = c(1,2,3,4,5))
  plot(c(2,3,4), o_time[,2], type = 'l', ylab = 'RMSE (m)', xlab = 'Number of sampled timepoints', xaxt = 'n', ylim = c(min(o_time[,2]),max(o_time[,2]) ) ); axis(1, at = c(1,2,3,4,5), labels = c(1,2,3,4,5))
  
  plot(c(5,10,20), o_sep[,1], type = 'l', ylab = 'RMSE (N)', xlab = 'tsep', xaxt = 'n', ylim = c(min(o_sep[,1]),max(o_sep[,1]) ) ) ; axis(1, at = c(5,10,20), labels = c(5,10,20))
  plot(c(5,10,20), o_sep[,2], type = 'l', ylab = 'RMSE (m)', xlab = 'tsep', xaxt = 'n', ylim = c(min(o_sep[,2]),max(o_sep[,2]) ) ); axis(1, at = c(5,10,20), labels = c(5,10,20))
  }
}

# cross validation - small m
{
  selected_small
  
  # do ABC
  {
  cv1_time_m_small = cross_validate(do_abc_PLS, SS_sim_1, selected_small, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2_time_m_small = cross_validate(do_abc_PLS, SS_sim_2, selected_small, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv3_time_m_small = cross_validate(do_abc_PLS, SS_sim_3, selected_small, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv4_time_m_small = cross_validate(do_abc_PLS, SS_sim_4, selected_small, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  cv2.10_time_m_small = cross_validate(do_abc_PLS, SS_sim_2.10, selected_small, n_comps =  min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2.5_time_m_small = cross_validate(do_abc_PLS, SS_sim_2.5, selected_small, n_comps =  min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  #cv1_time_AF_m = cross_validate(do_abc_PLS, SS_sim_1_AF, selected, n_comps =  min(nc, dim(SS_sim_1_AF)[2]), Nm_values_sim, method = method, tol = tol)
  cv2_time_AF_m_small = cross_validate(do_abc_PLS, SS_sim_2_AF, selected_small, n_comps =  min(nc, dim(SS_sim_2_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv3_time_AF_m_small = cross_validate(do_abc_PLS, SS_sim_3_AF, selected_small, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv4_time_AF_m_small = cross_validate(do_abc_PLS, SS_sim_4_AF, selected_small, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  cv2.10_time_AF_m_small = cross_validate(do_abc_PLS, SS_sim_2.10_AF, selected_small, n_comps =  min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2.5_time_AF_m_small = cross_validate(do_abc_PLS, SS_sim_2.5_AF, selected_small, n_comps =  min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  }
  
  
  # processing
  {
  o_n_tp_small = array(NA, dim = c(0,2))
  RMS_cv1 = RMS_err(cv1_time_m_small); o_n_tp_small = rbind(o_n_tp_small,colMeans(RMS_cv1));
  RMS_cv2 = RMS_err(cv2_time_m_small); o_n_tp_small = rbind(o_n_tp_small,colMeans(RMS_cv2));
  RMS_cv3 = RMS_err(cv3_time_m_small); o_n_tp_small = rbind(o_n_tp_small,colMeans(RMS_cv3))
  RMS_cv4 = RMS_err(cv4_time_m_small); o_n_tp_small = rbind(o_n_tp_small,colMeans(RMS_cv4))
  
  
  CIsLD_tp_small = array(NA, dim = c(4,2,2))
  CIsLD_tp_small[1,,] =  cbind(do_boot(RMS_cv1[,1]), do_boot(RMS_cv1[,2]))
  CIsLD_tp_small[2,,] =  cbind(do_boot(RMS_cv2[,1]), do_boot(RMS_cv2[,2]))
  CIsLD_tp_small[3,,] =  cbind(do_boot(RMS_cv3[,1]), do_boot(RMS_cv3[,2]))
  CIsLD_tp_small[4,,] =  cbind(do_boot(RMS_cv4[,1]), do_boot(RMS_cv4[,2]))

  
  o_n_tp_AF_small = array(NA, dim = c(0,2))
  #RMS_cv1 = RMS_err(cv1_time_AF_m_small); o_n_tp_AF_small = rbind(o_n_tp_AF_small,colMeans(RMS_cv1));
  RMS_cv2 = RMS_err(cv2_time_AF_m_small); o_n_tp_AF_small = rbind(o_n_tp_AF_small,colMeans(RMS_cv2));
  RMS_cv3 = RMS_err(cv3_time_AF_m_small); o_n_tp_AF_small = rbind(o_n_tp_AF_small,colMeans(RMS_cv3))
  RMS_cv4 = RMS_err(cv4_time_AF_m_small); o_n_tp_AF_small = rbind(o_n_tp_AF_small,colMeans(RMS_cv4))
  
  
  CIsAF_tp_small = array(NA, dim = c(3,2,2))
  #CIsAF_tp_small[1,,] =  cbind(do_boot(RMS_cv2[,1]), do_boot(RMS_cv2[,2]))
  CIsAF_tp_small[1,,] =  cbind(do_boot(RMS_cv2[,1]), do_boot(RMS_cv2[,2]))
  CIsAF_tp_small[2,,] =  cbind(do_boot(RMS_cv3[,1]), do_boot(RMS_cv3[,2]))
  CIsAF_tp_small[3,,] =  cbind(do_boot(RMS_cv4[,1]), do_boot(RMS_cv4[,2]))
  }
  
  
  # old
  {
  plot(c(1,2,3,4), o_n_tp[,1], type = 'p', pch = 16, ylab = 'RMSE (N)', xlab = 'Number of sampled subpopulations', xaxt = 'n', ylim = c(min(CIsA[,,1], CIsB[,,1]),max(CIsA[,,1], CIsB[,,1] ) ) ); axis(1, at = c(1,2,3,5,9), labels = c(1,2,3,5,9))
  for (i in 1:4){lines(rep(c(1,2,3,5,9)[i],2),CIsA[i,,1])}
  points(c(2,3,4), o_n_tp_AF[,1],pch = 16, col = 'red')
  for (i in 1:3){lines(rep(c(2,3,5,9)[i],2),CIsB[i,,1], col = 'red')}
  
  
  plot(c(1,2,3,4), o_n_tp[,2], type = 'p', pch = 16, ylab = 'RMSE (m)', xlab = 'Number of sampled subpopulations', xaxt = 'n', ylim = c(min(CIsA[,,2], CIsB[,,2]),max(CIsA[,,2], CIsB[,,2] ) ) ); axis(1, at = c(1,2,3,5,9), labels = c(1,2,3,5,9))
  for (i in 1:4){lines(rep(c(1,2,3,5,9)[i],2),CIsA[i,,2])}
  points(c(2,3,4), o_n_tp_AF[,2],pch = 16, col = 'red')
  for (i in 1:3){lines(rep(c(2,3,5,9)[i],2),CIsB[i,,2], col = 'red')}
  
  

  
  plot(c(1,2,3,4), o_n_tp[,1], type = 'l', ylab = 'RMSE (N)', xlab = 'Number of sampled tiempoints', xaxt = 'n', ylim = c(min(o_n_tp[,1],o_n_tp_AF[,1]),max(o_n_tp[,1], o_n_tp_AF[,1] ) ) ); axis(1, at = c(1,2,3,4), labels = c(1,2,3,4))
  lines(c(2,3,4), o_n_tp_AF[,1], col = 'red')
  plot(c(1,2,3,4), o_n_tp[,2], type = 'l', ylab = 'RMSE (m)', xlab = 'Number of sampled timepoints', xaxt = 'n', ylim = c(min(o_n_tp[,2],o_n_tp_AF[,2]),max(o_n_tp[,2], o_n_tp_AF[,2] ) ) ); axis(1, at = c(1,2,3,4), labels = c(1,2,3,4))
  lines(c(2,3,4), o_n_tp_AF[,2], col = 'red')
  

  plot(c(5,10,20), o_sep[,1], type = 'l', ylab = 'RMSE (N)', xlab = 'tsep', xaxt = 'n', ylim = c(min(o_sep[,1]),max(o_sep[,1]) ) ) ; axis(1, at = c(5,10,20), labels = c(5,10,20))
  plot(c(5,10,20), o_sep[,2], type = 'l', ylab = 'RMSE (m)', xlab = 'tsep', xaxt = 'n', ylim = c(min(o_sep[,2]),max(o_sep[,2]) ) ); axis(1, at = c(5,10,20), labels = c(5,10,20))
  }
}

# cross validation - large m
{
  selected_large
  
  # do_ABC
  {
  cv1_time_m_large  = cross_validate(do_abc_PLS, SS_sim_1, selected_large, n_comps =  nc, Nm_values_sim,method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2_time_m_large  = cross_validate(do_abc_PLS, SS_sim_2, selected_large, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv3_time_m_large  = cross_validate(do_abc_PLS, SS_sim_3, selected_large, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv4_time_m_large  = cross_validate(do_abc_PLS, SS_sim_4, selected_large, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  cv2.10_time_m_large  = cross_validate(do_abc_PLS, SS_sim_2.10, selected_large, n_comps =  min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2.5_time_m_large  = cross_validate(do_abc_PLS, SS_sim_2.5, selected_large, n_comps =  min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  #cv1_time_AF_m_large  = cross_validate(do_abc_PLS, SS_sim_1_AF, selected_large, n_comps =  min(nc, dim(SS_sim_1_AF)[2]), Nm_values_sim, method = method, tol = tol)
  cv2_time_AF_m_large  = cross_validate(do_abc_PLS, SS_sim_2_AF, selected_large, n_comps =  min(nc, dim(SS_sim_2_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv3_time_AF_m_large  = cross_validate(do_abc_PLS, SS_sim_3_AF, selected_large, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv4_time_AF_m_large  = cross_validate(do_abc_PLS, SS_sim_4_AF, selected_large, n_comps =  nc, Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  
  cv2.10_time_AF_m_large  = cross_validate(do_abc_PLS, SS_sim_2.10_AF, selected_large, n_comps =  min(nc, dim(SS_sim_2.10_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  cv2.5_time_AF_m_large  = cross_validate(do_abc_PLS, SS_sim_2.5_AF, selected_large, n_comps =  min(nc, dim(SS_sim_2.5_AF)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
  }
  
  # processing
  {
  o_n_tp_large = array(NA, dim = c(0,2))
  RMS_cv1 = RMS_err(cv1_time_m_large ); o_n_tp_large = rbind(o_n_tp_large,colMeans(RMS_cv1));
  RMS_cv2 = RMS_err(cv2_time_m_large ); o_n_tp_large = rbind(o_n_tp_large,colMeans(RMS_cv2));
  RMS_cv3 = RMS_err(cv3_time_m_large ); o_n_tp_large = rbind(o_n_tp_large,colMeans(RMS_cv3))
  RMS_cv4 = RMS_err(cv4_time_m_large ); o_n_tp_large = rbind(o_n_tp_large,colMeans(RMS_cv4))
  
  
  
  CIsLD_tp_large = array(NA, dim = c(4,2,2))
  CIsLD_tp_large[1,,] =  cbind(do_boot(RMS_cv1[,1]), do_boot(RMS_cv1[,2]))
  CIsLD_tp_large[2,,] =  cbind(do_boot(RMS_cv2[,1]), do_boot(RMS_cv2[,2]))
  CIsLD_tp_large[3,,] =  cbind(do_boot(RMS_cv3[,1]), do_boot(RMS_cv3[,2]))
  CIsLD_tp_large[4,,] =  cbind(do_boot(RMS_cv4[,1]), do_boot(RMS_cv4[,2]))

  
  o_n_tp_AF_large = array(NA, dim = c(0,2))
  #RMS_cv1 = RMS_err(cv1_time_AF_m_large ); o_n_tp_AF_large = rbind(o_n_tp_AF_large,colMeans(RMS_cv1));
  RMS_cv2 = RMS_err(cv2_time_AF_m_large ); o_n_tp_AF_large = rbind(o_n_tp_AF_large,colMeans(RMS_cv2));
  RMS_cv3 = RMS_err(cv3_time_AF_m_large ); o_n_tp_AF_large = rbind(o_n_tp_AF_large,colMeans(RMS_cv3))
  RMS_cv4 = RMS_err(cv4_time_AF_m_large ); o_n_tp_AF_large = rbind(o_n_tp_AF_large,colMeans(RMS_cv4))
  
  
  CIsAF_tp_large = array(NA, dim = c(3,2,2))
  #CIsAF_tp_large[1,,] =  cbind(do_boot(RMS_cv2[,1]), do_boot(RMS_cv2[,2]))
  CIsAF_tp_large[1,,] =  cbind(do_boot(RMS_cv2[,1]), do_boot(RMS_cv2[,2]))
  CIsAF_tp_large[2,,] =  cbind(do_boot(RMS_cv3[,1]), do_boot(RMS_cv3[,2]))
  CIsAF_tp_large[3,,] =  cbind(do_boot(RMS_cv4[,1]), do_boot(RMS_cv4[,2]))
  }
  
  
  # old
  {
  plot(c(1,2,3,4), o_n_tp[,1], type = 'p', pch = 16, ylab = 'RMSE (N)', xlab = 'Number of sampled subpopulations', xaxt = 'n', ylim = c(min(CIsA[,,1], CIsB[,,1]),max(CIsA[,,1], CIsB[,,1] ) ) ); axis(1, at = c(1,2,3,4), labels = c(1,2,3,4))
  for (i in 1:4){lines(rep(c(1,2,3,5,9)[i],2),CIsA[i,,1])}
  points(c(2,3,4), o_n_tp_AF[,1],pch = 16, col = 'red')
  for (i in 1:3){lines(rep(c(2,3,5,9)[i],2),CIsB[i,,1], col = 'red')}
  
  
  plot(c(1,2,3,4), o_n_tp[,2], type = 'p', pch = 16, ylab = 'RMSE (m)', xlab = 'Number of sampled subpopulations', xaxt = 'n', ylim = c(min(CIsA[,,2], CIsB[,,2]),max(CIsA[,,2], CIsB[,,2] ) ) ); axis(1, at = c(1,2,3,4), labels = c(1,2,3,4))
  for (i in 1:4){lines(rep(c(1,2,3,5,9)[i],2),CIsA[i,,2])}
  points(c(2,3,4), o_n_tp_AF[,2],pch = 16, col = 'red')
  for (i in 1:3){lines(rep(c(2,3,5,9)[i],2),CIsB[i,,2], col = 'red')}

  
  plot(c(1,2,3,4), o_n_tp[,1], type = 'l', ylab = 'RMSE (N)', xlab = 'Number of sampled timepoints', xaxt = 'n', ylim = c(min(o_n_tp[,1],o_n_tp_AF[,1]),max(o_n_tp[,1], o_n_tp_AF[,1] ) ) ); axis(1, at = c(1,2,3,4), labels = c(1,2,3,4))
  lines(c(2,3,4), o_n_tp_AF[,1], col = 'red')
  plot(c(1,2,3,4), o_n_tp[,2], type = 'l', ylab = 'RMSE (m)', xlab = 'Number of sampled timepoints', xaxt = 'n', ylim = c(min(o_n_tp[,2],o_n_tp_AF[,2]),max(o_n_tp[,2], o_n_tp_AF[,2] ) ) ); axis(1, at = c(1,2,3,4), labels = c(1,2,3,4))
  lines(c(2,3,4), o_n_tp_AF[,2], col = 'red')
  
  
  plot(c(5,10,20), o_sep[,1], type = 'l', ylab = 'RMSE (N)', xlab = 'tsep', xaxt = 'n', ylim = c(min(o_sep[,1]),max(o_sep[,1]) ) ) ; axis(1, at = c(5,10,20), labels = c(5,10,20))
  plot(c(5,10,20), o_sep[,2], type = 'l', ylab = 'RMSE (m)', xlab = 'tsep', xaxt = 'n', ylim = c(min(o_sep[,2]),max(o_sep[,2]) ) ); axis(1, at = c(5,10,20), labels = c(5,10,20))
  
  }
}



# cross val plot
pdf(file = '../final_figures/CV_n_tp.pdf', height = 7, width = 7, pointsize = 7)
par(mfrow = c(2,1), mai = c(.5,.6,.4,.4), oma = c(5,1,0,0))
#compare_CV(LD = o_n_sp, AF = o_n_sp_AF, LD_ci = CIsLD, AF_ci = CIsAF)

compare_CV_tp(LD = o_n_tp_small, AF = o_n_tp_AF_small, LD_ci = CIsLD_tp_small, AF_ci = CIsAF_tp_small)
fig_label('A', cex = 2)
compare_CV_tp(LD = o_n_tp_large, AF = o_n_tp_AF_large, LD_ci = CIsLD_tp_large, AF_ci = CIsAF_tp_large)
fig_label('B', cex = 2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('bottom',cex = .9, legend = c('N - LD      ', 'N - AF      ', 'm - LD      ', 'm - AF      ' ), pch = c(1,1,6, 6), col = c('black', 'red', 'black', 'red'), horiz = T, inset = rep(0.04,4))
graphics.off()

# cross validation - small N
{
  selected = select_bounded(100, Nm_values_sim, m_range = c(0,0.8), N_range = c(400,2000))
  
  cv1_time_N = cross_validate(do_abc_PLS, SS_sim_1, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv2_time_N = cross_validate(do_abc_PLS, SS_sim_2, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv3_time_N = cross_validate(do_abc_PLS, SS_sim_3, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv4_time_N = cross_validate(do_abc_PLS, SS_sim_4, selected, n_comps =  nc, Nm_values_sim,method = method, tol = tol)
  
  cv2.10_time_N = cross_validate(do_abc_PLS, SS_sim_2.10, selected, n_comps =  nc, Nm_values_sim,method = method, tol = tol)
  cv2.20_time_N = cross_validate(do_abc_PLS, SS_sim_2.20, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  
  cv1_time_AF_N = cross_validate(do_abc_PLS, SS_sim_1_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv2_time_AF_N = cross_validate(do_abc_PLS, SS_sim_2_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv3_time_AF_N = cross_validate(do_abc_PLS, SS_sim_3_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv4_time_AF_N = cross_validate(do_abc_PLS, SS_sim_4_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  
  cv2.10_time_N = cross_validate(do_abc_PLS, SS_sim_2.10_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  cv2.20_time_N = cross_validate(do_abc_PLS, SS_sim_2.20_AF, selected, n_comps =  nc, Nm_values_sim, method = method, tol = tol)
  
  
  
  
  o_n_tp = array(NA, dim = c(0,2))
  RMS_cv1 = RMS_err(cv1_time_N); o_n_tp = rbind(o_n_tp,colMeans(RMS_cv1));
  RMS_cv2 = RMS_err(cv2_time_N); o_n_tp = rbind(o_n_tp,colMeans(RMS_cv2));
  RMS_cv3 = RMS_err(cv3_time_N); o_n_tp = rbind(o_n_tp,colMeans(RMS_cv3))
  RMS_cv4 = RMS_err(cv4_time_N); o_n_tp = rbind(o_n_tp,colMeans(RMS_cv4))
  
  
  
  o_n_tp_AF = array(NA, dim = c(0,2))
  RMS_cv1 = RMS_err(cv1_time_AF_N); o_n_tp_AF = rbind(o_n_tp_AF,colMeans(RMS_cv1));
  RMS_cv2 = RMS_err(cv2_time_AF_N); o_n_tp_AF = rbind(o_n_tp_AF,colMeans(RMS_cv2));
  RMS_cv3 = RMS_err(cv3_time_AF_N); o_n_tp_AF = rbind(o_n_tp_AF,colMeans(RMS_cv3))
  RMS_cv4 = RMS_err(cv4_time_AF_N); o_n_tp_AF = rbind(o_n_tp_AF,colMeans(RMS_cv4))
  
  
  
  
  plot(c(1,2,3,4), o_n_tp[,1], type = 'l', ylab = 'RMSE (N)', xlab = 'Number of sampled subpopulations', xaxt = 'n', ylim = c(min(o_n_tp[,1],o_n_tp_AF[,1]),max(o_n_tp[,1], o_n_tp_AF[,1] ) ) ); axis(1, at = c(1,2,3,4), labels = c(1,2,3,4))
  lines(c(1,2,3,4), o_n_tp_AF[,1], col = 'red')
  plot(c(1,2,3,4), o_n_tp[,2], type = 'l', ylab = 'RMSE (m)', xlab = 'Number of sampled subpopulations', xaxt = 'n', ylim = c(min(o_n_tp[,2],o_n_tp_AF[,2]),max(o_n_tp[,2], o_n_tp_AF[,2] ) ) ); axis(1, at = c(1,2,3,4), labels = c(1,2,3,4))
  lines(c(1,2,3,4), o_n_tp_AF[,2], col = 'red')
  
  
  
  
  plot(c(2,3,4), o_time[,1], type = 'l', ylab = 'RMSE (N)', xlab = 'Number of sampled timepoints', xaxt = 'n', ylim = c(min(o_time[,1]),max(o_time[,1]) ) ) ; axis(1, at = c(1,2,3,4,5), labels = c(1,2,3,4,5))
  plot(c(2,3,4), o_time[,2], type = 'l', ylab = 'RMSE (m)', xlab = 'Number of sampled timepoints', xaxt = 'n', ylim = c(min(o_time[,2]),max(o_time[,2]) ) ); axis(1, at = c(1,2,3,4,5), labels = c(1,2,3,4,5))
  
  plot(c(5,10,20), o_sep[,1], type = 'l', ylab = 'RMSE (N)', xlab = 'tsep', xaxt = 'n', ylim = c(min(o_sep[,1]),max(o_sep[,1]) ) ) ; axis(1, at = c(5,10,20), labels = c(5,10,20))
  plot(c(5,10,20), o_sep[,2], type = 'l', ylab = 'RMSE (m)', xlab = 'tsep', xaxt = 'n', ylim = c(min(o_sep[,2]),max(o_sep[,2]) ) ); axis(1, at = c(5,10,20), labels = c(5,10,20))
}


# make figure
# make figure

{
  pdf(file = '../final_figures/n_tp_CV.pdf', height = 7, width = 7, pointsize = 7)
  par(mfrow = c(1,1), mai = c(1,.8,0.4,0.4))
  compare_CV_tp(LD = o_n_tp, AF = o_n_tp_AF, LD_ci = CIsLD_tp, AF_ci = CIsAF_tp)
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend('bottom',cex = .9, legend = c('N - LD      ', 'N - AF      ', 'm - LD      ', 'm - AF      ' ), pch = c(1,1,6, 6), col = c('black', 'red', 'black', 'red'), horiz = T, inset = rep(0.02,4))
  graphics.off()
}


}



##################################
##################################
######## Across space + time #######
##################################
##################################

# 3D surface plot of LOOCV RMSE with different allocations across space and time
{
S_total = 1000


{
space_strategies_sim = list(c(41), c(41,42),c(40,41,42), c(40,41,42,50,51), sample_subpops_sim)
space_strategies_real = list(c(61), c(61,62),c(60,61,62), c(60,61,62,72,73), sample_subpops_real)

time_strategies = list(c(0), c(0,20), c(0,10,20), c(0,5,10,20))

selected 
RMSE_out_s_and_t = array(NA, dim = c(length(space_strategies_sim), length(time_strategies),2))

for (space_strategy_idx in 1:5){
  for (time_strategy_idx in 1:4){
    subpops_sim = space_strategies_sim[[space_strategy_idx]]; subpops_real = space_strategies_real[[space_strategy_idx]]
    timepoints = time_strategies[[time_strategy_idx]]
    S_loc = round(S_total / length(subpops_sim) / length(timepoints))
    SS = make_SS(sim_out, real_out_list, S_select = S_loc, gens_select = timepoints, subpops_select_sim = subpops_sim, subpops_select_real = subpops_real, assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
    SS_sim = SS$sim; SS_real = SS$real
    
    cv = cross_validate(do_abc_PLS, SS_sim, selected, n_comps =  min(nc, dim(SS_sim)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    RMS_cv = RMS_err(cv); means = colMeans(RMS_cv)
    RMSE_out_s_and_t[space_strategy_idx, time_strategy_idx,] = means
  }
}


RMSE_out_s_and_t_AF = array(NA, dim = c(length(space_strategies_sim), length(time_strategies),2))

for (space_strategy_idx in 1:5){
  for (time_strategy_idx in 2:4){
    if (!(time_strategy_idx == 2 & space_strategy_idx == 1 )){
    subpops_sim = space_strategies_sim[[space_strategy_idx]]; subpops_real = space_strategies_real[[space_strategy_idx]]
    timepoints = time_strategies[[time_strategy_idx]]
    S_loc = round(get_diploid_equivalent(1e100,30*round(S_total / length(subpops_sim) / length(timepoints)),30))
    SS = make_SS(sim_out, real_out_list, S_select = S_loc, gens_select = timepoints, subpops_select_sim = subpops_sim, subpops_select_real = subpops_real, assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
    SS_sim = SS$sim; SS_real = SS$real
    
    cv = cross_validate(do_abc_PLS, SS_sim, selected, n_comps =  min(nc, dim(SS_sim)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
    RMS_cv = RMS_err(cv); means = colMeans(RMS_cv)
    RMSE_out_s_and_t_AF[space_strategy_idx, time_strategy_idx,] = means
    }
  }
}

saveRDS(RMSE_out_s_and_t, '../data/RMSE_out_s_and_t_scaled.RDS') ### 1000 full genomes
saveRDS(RMSE_out_s_and_t_AF, '../data/RMSE_out_s_and_t_AF_scaled.RDS') ### 1000 full genomes

}

# small m
{
  space_strategies_sim = list(c(41), c(41,42),c(40,41,42), c(40,41,42,50,51), sample_subpops_sim)
  space_strategies_real = list(c(61), c(61,62),c(60,61,62), c(60,61,62,72,73), sample_subpops_real)
  
  time_strategies = list(c(0), c(0,20), c(0,10,20), c(0,5,10,20))
  
  selected_small = select_bounded(1000, Nm_values_sim, N_range = c(1000,20000), m_range = c(0,0.1))
  RMSE_out_s_and_t_m = array(NA, dim = c(length(space_strategies_sim), length(time_strategies),2))
  
  for (space_strategy_idx in 1:5){
    for (time_strategy_idx in 1:4){
      subpops_sim = space_strategies_sim[[space_strategy_idx]]; subpops_real = space_strategies_real[[space_strategy_idx]]
      timepoints = time_strategies[[time_strategy_idx]]
      S_loc = round(S_total / length(subpops_sim) / length(timepoints))
      SS = make_SS(sim_out, real_out_list, S_select = S_loc, gens_select = timepoints, subpops_select_sim = subpops_sim, subpops_select_real = subpops_real, assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      
      cv = cross_validate(do_abc_PLS, SS_sim, selected_small, n_comps =  min(nc, dim(SS_sim)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
      RMS_cv = RMS_err(cv); means = colMeans(RMS_cv)
      RMSE_out_s_and_t_m[space_strategy_idx, time_strategy_idx,] = means
    }
  }
  
  
  RMSE_out_s_and_t_AF_m = array(NA, dim = c(length(space_strategies_sim), length(time_strategies),2))
  
  for (space_strategy_idx in 1:5){
    for (time_strategy_idx in 2:4){
      if (time_strategy_idx == 2 & space_strategy_idx == 1 ) {}
      else{
      subpops_sim = space_strategies_sim[[space_strategy_idx]]; subpops_real = space_strategies_real[[space_strategy_idx]]
      timepoints = time_strategies[[time_strategy_idx]]
      S_loc = round(get_diploid_equivalent(1e100,30*round(S_total / length(subpops_sim) / length(timepoints)),30))
      SS = make_SS(sim_out, real_out_list, S_select = S_loc, gens_select = timepoints, subpops_select_sim = subpops_sim, subpops_select_real = subpops_real, assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      
      cv = cross_validate(do_abc_PLS, SS_sim, selected_small, n_comps =  min(nc, dim(SS_sim)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
      RMS_cv = RMS_err(cv); means = colMeans(RMS_cv)
      RMSE_out_s_and_t_AF_m[space_strategy_idx, time_strategy_idx,] = means
      }
    }
  }
  
  saveRDS(RMSE_out_s_and_t_m, '../data/RMSE_out_s_and_t_small_scaled.RDS') ### 1000 full genomes
  saveRDS(RMSE_out_s_and_t_AF_m, '../data/RMSE_out_s_and_t_AF_small_scaled.RDS')  ### 1000 genomes, allocated to poolseq
}


# large m
{
  space_strategies_sim = list(c(41), c(41,42),c(40,41,42), c(40,41,42,50,51), sample_subpops_sim)
  space_strategies_real = list(c(61), c(61,62),c(60,61,62), c(60,61,62,72,73), sample_subpops_real)
  
  time_strategies = list(c(0), c(0,20), c(0,10,20), c(0,5,10,20))
  selected_large = select_bounded(1000, Nm_values_sim, N_range = c(1000,20000), m_range = c(0,0.1))
  
  RMSE_out_s_and_t_m_large = array(NA, dim = c(length(space_strategies_sim), length(time_strategies),2))
  
  for (space_strategy_idx in 1:5){
    for (time_strategy_idx in 1:4){
      subpops_sim = space_strategies_sim[[space_strategy_idx]]; subpops_real = space_strategies_real[[space_strategy_idx]]
      timepoints = time_strategies[[time_strategy_idx]]
      S_loc = round(S_total / length(subpops_sim) / length(timepoints))
      SS = make_SS(sim_out, real_out_list, S_select = S_loc, gens_select = timepoints, subpops_select_sim = subpops_sim, subpops_select_real = subpops_real, assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      
      cv = cross_validate(do_abc_PLS, SS_sim, selected_large, n_comps =  min(nc, dim(SS_sim)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
      RMS_cv = RMS_err(cv); means = colMeans(RMS_cv)
      RMSE_out_s_and_t_m_large[space_strategy_idx, time_strategy_idx,] = means
    }
  }
  
  
  RMSE_out_s_and_t_AF_m_large = array(NA, dim = c(length(space_strategies_sim), length(time_strategies),2))
  
  for (space_strategy_idx in 1:5){
    for (time_strategy_idx in 2:4){
      if (time_strategy_idx == 2 & space_strategy_idx == 1 ) {}
      else{
      subpops_sim = space_strategies_sim[[space_strategy_idx]]; subpops_real = space_strategies_real[[space_strategy_idx]]
      timepoints = time_strategies[[time_strategy_idx]]
      S_loc = round(get_diploid_equivalent(1e100,30*round(S_total / length(subpops_sim) / length(timepoints)),30))
      SS = make_SS(sim_out, real_out_list, S_select = S_loc, gens_select = timepoints, subpops_select_sim = subpops_sim, subpops_select_real = subpops_real, assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      
      cv = cross_validate(do_abc_PLS, SS_sim, selected, n_comps =  min(nc, dim(SS_sim)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
      RMS_cv = RMS_err(cv); means = colMeans(RMS_cv)
      RMSE_out_s_and_t_AF_m_large[space_strategy_idx, time_strategy_idx,] = means
      }
    }
  }
  
  saveRDS(RMSE_out_s_and_t_m_large, '../data/RMSE_out_s_and_t_large_scaled.RDS') ### 1000 full genomes
  saveRDS(RMSE_out_s_and_t_AF_m_large, '../data/RMSE_out_s_and_t_AF_large_scaled.RDS')  ### 1000 genomes, allocated to poolseq
}



{

dat2 <- utils::stack(data.frame(RMSE_out_s_and_t_AF_m[,-1,1]),)
dat2$n_sp <- as.character(c(1,2,3,5,9))
dat2$ind <- rep(c(2,3,4),each = 5)

g1 = ggplot(mapping = aes(x=ind, y=n_sp), dat2)
g1 = g1 +  geom_tile(aes(fill = values), dat2, show.legend = F)
g1 = g1 + geom_text(aes(label = round(values, 3)), dat2, size = 4) 
g1 = g1+ xlab(bquote(n[t]))+ylab(bquote(n[s]))
g1 = g1 +   scale_y_discrete( labels = c(1,2,3,5,9)) 
#g = g +  scale_x_discrete(labels = c(2, 3, 4), position = "top") 
g1 = g1 +  scale_fill_gradient(low="green", high="red")
g1 = g1 + theme_minimal()   +  theme(text = element_text(size = 11))+ labs(tag = "A")


dat2 <- utils::stack(data.frame(RMSE_out_s_and_t_AF_m[,-1,2]),)
dat2$n_sp <- as.character(c(1,2,3,5,9))
dat2$ind <- rep(c(2,3,4),each = 5)

g2 = ggplot(mapping = aes(x=ind, y=n_sp), dat2)
g2 = g2 +  geom_tile(aes(fill = values), dat2, show.legend = F)
g2 = g2 + geom_text(aes(label = round(values, 3)), dat2, size = 4) 
g2 = g2 + xlab(bquote(n[t]))+ylab(bquote(n[s]))
g2 = g2 +   scale_y_discrete( labels = c(1,2,3,5,9)) 
#g = g +  scale_x_discrete(labels = c(2, 3, 4), position = "top") 
g2 = g2 +  scale_fill_gradient(low="green", high="red")
g2 = g2 + theme_minimal()   +  theme(text = element_text(size = 11))+ labs(tag = "B")
g2


dat2 <- utils::stack(data.frame(RMSE_out_s_and_t_AF_m_large[,-1,1]),)
dat2$n_sp <- as.character(c(1,2,3,5,9))
dat2$ind <- rep(c(2,3,4),each = 5)

g3 = ggplot(mapping = aes(x=ind, y=n_sp), dat2)
g3 = g3 +  geom_tile(aes(fill = values), dat2, show.legend = F)
g3 = g3 + geom_text(aes(label = round(values, 3)), dat2, size = 4) 
g3 = g3 + xlab(bquote(n[t]))+ylab(bquote(n[s]))
g3 = g3 +   scale_y_discrete( labels = c(1,2,3,5,9)) 
#g = g +  scale_x_discrete(labels = c(2, 3, 4), position = "top") 
g3 = g3 +  scale_fill_gradient(low="green", high="red")
g3 = g3 + theme_minimal()   +  theme(text = element_text(size = 11))+ labs(tag = "C")
g3


dat2 <- utils::stack(data.frame(RMSE_out_s_and_t_AF_m_large[,-1,2]),)
dat2$n_sp <- as.character(c(1,2,3,5,9))
dat2$ind <- rep(c(2,3,4),each = 5)

g4 = ggplot(mapping = aes(x=ind, y=n_sp), dat2)
g4 = g4 +  geom_tile(aes(fill = values), dat2, show.legend = F)
g4 = g4 + geom_text(aes(label = round(values, 3)), dat2,size = 4) 
g4 = g4 + xlab(bquote(n[t]))+ylab(bquote(n[s]))
g4 = g4 +   scale_y_discrete( labels = c(1,2,3,5,9) )
g4 = g4 +  scale_fill_gradient(low="green", high="red")
g4 = g4 + theme_minimal()   +  theme(text = element_text(size = 11))+ labs(tag = "D")
g4

pdf('../final_figures/heatmap.pdf', height = 2.8, width = 5)
grid.arrange( grobs = lapply(list(g1,g2,g3,g4), "+", theme(plot.margin=margin(.05,.05,.05,.05)))  , ncol = 2, )
graphics.off()
}




#small N
{
  space_strategies_sim = list(c(41), c(41,42),c(40,41,42), c(40,41,42,50,51), sample_subpops_sim)
  space_strategies_real = list(c(61), c(61,62),c(60,61,62), c(60,61,62,72,73), sample_subpops_real)
  
  time_strategies = list(c(0), c(0,20), c(0,10,20), c(0,5,10,20))
  
  selected = select_bounded(100, Nm_values_sim, m_range = c(0,0.8), N_range = c(400,2000))
  RMSE_out_s_and_t_N = array(NA, dim = c(length(space_strategies_sim), length(time_strategies),2))
  
  for (space_strategy_idx in 1:5){
    for (time_strategy_idx in 1:4){
      subpops_sim = space_strategies_sim[[space_strategy_idx]]; subpops_real = space_strategies_real[[space_strategy_idx]]
      timepoints = time_strategies[[time_strategy_idx]]
      S_loc = round(S_total / length(subpops_sim) / length(timepoints))
      SS = make_SS(sim_out, real_out_list, S_select = S_loc, gens_select = timepoints, subpops_select_sim = subpops_sim, subpops_select_real = subpops_real, assume_equal = F, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      
      cv = cross_validate(do_abc_PLS, SS_sim, selected, n_comps =  min(nc, dim(SS_sim)[2]), Nm_values_sim, method = method, tol = tol)
      RMS_cv = RMS_err(cv); means = colMeans(RMS_cv)
      RMSE_out_s_and_t_N[space_strategy_idx, time_strategy_idx,] = means
    }
  }
  
  
  RMSE_out_s_and_t_AF_N = array(NA, dim = c(length(space_strategies_sim), length(time_strategies),2))
  
  for (space_strategy_idx in 1:5){
    for (time_strategy_idx in 2:4){
      if (time_strategy_idx == 2 & space_strategy_idx == 1 ) break
      subpops_sim = space_strategies_sim[[space_strategy_idx]]; subpops_real = space_strategies_real[[space_strategy_idx]]
      timepoints = time_strategies[[time_strategy_idx]]
      S_loc = round(get_diploid_equivalent(1e100,30*round(S_total / length(subpops_sim) / length(timepoints)),30))
      SS = make_SS(sim_out, real_out_list, S_select = S_loc, gens_select = timepoints, subpops_select_sim = subpops_sim, subpops_select_real = subpops_real, assume_equal = F, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      
      cv = cross_validate(do_abc_PLS, SS_sim, selected, n_comps =  min(nc, dim(SS_sim)[2]), Nm_values_sim, method = method, tol = tol)
      RMS_cv = RMS_err(cv); means = colMeans(RMS_cv)
      RMSE_out_s_and_t_AF_N[space_strategy_idx, time_strategy_idx,] = means
    }
  }
  
  RMSE_out_s_and_t_N ### 1000 full genomes
  RMSE_out_s_and_t_AF_N ### 1000 genomes, allocated to poolseq
}

## idk
{
subpops_sim = space_strategies_sim[[2]]; subpops_real = space_strategies_real[[2]]
timepoints = time_strategies[[2]]
S_loc = round(S_total / length(subpops_sim) / length(timepoints))
SS = make_SS(sim_out, real_out_list, S_select = S_loc, gens_select = timepoints, subpops_select_sim = subpops_sim, subpops_select_real = subpops_real, assume_equal = F, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim = SS$sim; SS_real = SS$real

cv = cross_validate(do_abc_PLS, SS_sim, selected, n_comps =  min(5, dim(SS_sim)[2]), Nm_values_sim, method = 'neuralnet', tol = .1)
RMS_cv = RMS_err(cv); means = colMeans(RMS_cv)
}

}



##################################
##################################
######## Smaller sample size #######
##################################
##################################


gens_select = c(0,10,20)
S_per_timepoint = 1000/3

S1 = round(S_per_timepoint)
S2 = round(S_per_timepoint / 2)
S3 = round(S_per_timepoint / 3)
S5 = round(S_per_timepoint / 5)
S9 = round(S_per_timepoint / 9)


# single subpopulation
{
  for (S in c(S1, S2, S3, S5, S9)){
    SS_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S*30, 30)), gens_select = gens_select, subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = T, LD_info = F,  focal_pop_sim = 61, focal_pop_real = 41, edge_length_sim = 9, edge_length_real = 11 )
    SS_sim_AF = SS_AF$sim[,]; SS_real_AF = SS_AF$real
    
    out_res = test_specific(do_abc_PLS, SS_sim_AF, SS_real_AF, n_comps = min(nc, dim(SS_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol, method = method, sizenet = sizenet, numnet = numnet )
    
    o_AF_sample_size[[i]] = out_res
    i=i+1
  }
  
  saveRDS(o_AF_sample_size, '../data/o_AF_sample_size.RDS')
  
  {
    o_AF_sample_size_means = list()
    i=1
    for (l in o_AF_sample_size){
      o_AF_sample_size_means[[i]] = apply(l, c(1,2,3), mean)
      i=i+1
    }
    pdf('../final_figures/AF_sample_size.pdf', height = 7, width = 7, pointsize = 7)
    plot_compare_N_neat(o_AF_sample_size_means, legend = paste('S =', c(1000, 333)), cols = c('black', 'black'), ltys = 1:100)
    graphics.off()
  }
    
}
  
# multiple subpopultions
{
  ### vary m, N = 2000
  {
  o_AF_S_m_N2000 = list()
  i = 1
  for (S in c(1000/6, 500/6, 200/6, 100/6)){
    S = round(S)
    SS_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S*30, 30)), gens_select = gens_select, subpops_select_sim = c(40:41), subpops_select_real = c(60:61), assume_equal = T, LD_info = F,  focal_pop_sim = 61, focal_pop_real = 41, edge_length_sim = 9, edge_length_real = 11 )
    SS_sim_AF = SS_AF$sim[,]; SS_real_AF = SS_AF$real
    
    out_res_3 = test_specific(do_abc_PLS, SS_sim_AF, SS_real_AF, n_comps = min(nc, dim(SS_sim_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(83,114), tol = tol, method = method, sizenet = sizenet, numnet = numnet )
    
    o_AF_S_m_N2000[[i]] = out_res_3
    i=i+1
  }
  guess = just_guess(Nm_values_real, c(83,114))
  o_AF_S_m_N2000[[i]] = guess
  
  saveRDS(o_AF_S_m_N2000, '../data/o_AF_S_m_N2000.RDS')
  o_AF_S_m_N2000 = readRDS('../data/o_AF_S_m_N2000.RDS')
  
  
  
  # make figure
  {
    o_AF_S_m_N2000_means = list()
    i=1
    for (l in o_AF_S_m_N2000){
      o_AF_S_m_N2000_means[[i]] = apply(l, c(1,2,3), mean)
      i=i+1
    }
    pdf('../final_figures/AF_S_m_N2000.pdf', height = 7, width = 7, pointsize = 7)
    plot_compare_m_neat(o_AF_S_m_N2000_means, legend = paste('S =', round(c(1000, 500, 200, 100))), cols = c('black', 'black'), ltys = 1:5)
    graphics.off()
  }
  }
  
  ### vary m, N = 10000
  {
    o_AF_S_m_N10000 = list()
    i = 1
    for (S in c(1000/9, 500/9, 200/9, 100/9, 50/9)){
      S = round(S)
      SS_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S*30, 30)), gens_select = gens_select, subpops_select_sim = c(40,41,42), subpops_select_real = c(60,61,62), assume_equal = T, LD_info = F,  focal_pop_sim = 61, focal_pop_real = 41, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim_AF = SS_AF$sim[,]; SS_real_AF = SS_AF$real
      
      out_res_3 = test_specific(do_abc_PLS, SS_sim_AF, SS_real_AF, n_comps = min(nc, dim(SS_sim_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(31,62), tol = tol, method = method, sizenet = sizenet, numnet = numnet )
      
      o_AF_S_m_N10000[[i]] = out_res_3
      i=i+1
    }
    guess = just_guess(Nm_values_real, c(31,62))
    o_AF_S_m_N10000[[i]] = guess
    saveRDS(o_AF_S_m_N10000, '../data/o_AF_S_m_N10000.RDS')
    o_AF_S_m_N10000 = readRDS('../data/o_AF_S_m_N10000.RDS')
    
    
    
    # make figure
    {
      o_AF_S_m_N10000_means = list()
      i=1
      for (l in o_AF_S_m_N10000){
        o_AF_S_m_N10000_means[[i]] = apply(l, c(1,2,3), mean)
        i=i+1
      }
      pdf('../final_figures/AF_S_m_N_10000.pdf', height = 7, width = 7, pointsize = 7)
      plot_compare_m_neat(o_AF_S_m_N10000_means[-5], legend = paste('S =', round(c(1000, 500, 200, 100))), cols = c('black', 'black'), ltys = 1:5)
      graphics.off()
    }
  }
  
  
  ### vary N, m = 0.05
  {
  o_AF_S_N_m05 = list()
  i = 1
  for (S in c(1000/6, 500/6, 200/6, 100/6, 50/6)){
    S = round(S)
    SS_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S*30, 30)), gens_select = gens_select, subpops_select_sim = c(40:41), subpops_select_real = c(60:61), assume_equal = T, LD_info = F,  focal_pop_sim = 61, focal_pop_real = 41, edge_length_sim = 9, edge_length_real = 11 )
    SS_sim_AF = SS_AF$sim[,]; SS_real_AF = SS_AF$real
    
    out_res_3 = test_specific(do_abc_PLS, SS_sim_AF, SS_real_AF, n_comps = min(nc, dim(SS_sim_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(1,27), tol = tol, method = method, sizenet = sizenet, numnet = numnet )
    
    o_AF_S_N_m05[[i]] = out_res_3
    i=i+1
  }
  guess = just_guess(Nm_values_real, c(1,27))
  o_AF_S_N_m05[[i]] = guess
  saveRDS(o_AF_S_N_m05, '../data/o_AF_S_N_m05.RDS')
  o_AF_S_N_m05 = readRDS('../data/o_AF_S_N_m05.RDS')
  
  # make fig
  {
    o_AF_S_N_m05_means = list()
    i=1
    for (l in o_AF_S_N_m05){
      o_AF_S_N_m05_means[[i]] = apply(l, c(1,2,3), mean)
      i=i+1
    }
    pdf('../final_figures/AF_S_N_m05.pdf', height = 7, width = 7, pointsize = 7)
    plot_compare_N_neat(o_AF_S_N_m05_means[-5], legend = paste('S =', round(c(1000, 500, 200, 100))), cols = c('black', 'black'), ltys = 1:5)
    graphics.off()
  }
  }
  
  ### vary N, m = 0.20
  {
    o_AF_S_N_m20 = list()
    i = 1
    for (S in c(1000/9, 500/9, 200/9, 100/9, 50/9)){
      S = round(S)
      SS_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, S*30, 30)), gens_select = gens_select, subpops_select_sim = c(40,41,42), subpops_select_real = c(60,61,62), assume_equal = T, LD_info = F,  focal_pop_sim = 61, focal_pop_real = 41, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim_AF = SS_AF$sim[,]; SS_real_AF = SS_AF$real
      
      out_res_3 = test_specific(do_abc_PLS, SS_sim_AF, SS_real_AF, n_comps = min(nc, dim(SS_sim_AF)[2]), Nm_values_sim, Nm_values_real, idx = c(63,76), tol = tol, method = method, sizenet = sizenet, numnet = numnet )
      
      o_AF_S_N_m20[[i]] = out_res_3
      i=i+1
    }
    guess = just_guess(Nm_values_real, c(63,76))
    o_AF_S_N_m20[[i]] = guess
    saveRDS(o_AF_S_N_m20, '../data/o_AF_S_N_m20.RDS')
    o_AF_S_N_m20 = readRDS('../data/o_AF_S_N_m20.RDS')
    
    # make fig
    {
      o_AF_S_N_m20_means = list()
      i=1
      for (l in o_AF_S_N_m20){
        o_AF_S_N_m20_means[[i]] = apply(l, c(1,2,3), mean)
        i=i+1
      }
      pdf('../final_figures/AF_S_N_m20.pdf', height = 7, width = 7, pointsize = 7)
      plot_compare_N_neat(o_AF_S_N_m20_means[-5], legend = paste('S =', round(c(1000, 500, 200, 100))), cols = c('black', 'black'), ltys = 1:5)
      graphics.off()
    }
  }
}

  
# combined fig
{
  N_plot = o_AF_S_N_m05_means[-5]
  m_plot = o_AF_S_m_N2000_means[-5]
  
  legends = c(paste('S =', round(c(1000, 500, 200, 100))), 'Prior')
  ltys = c(1,2,3,4,1)
  cols = c(rep('red',4), '#4daf4a')
  
  plot_compare_N_m_neat(N_plot, m_plot, legends, ltys = ltys, cols = cols, fname = '../final_figures/AF_small_N_and_m.pdf', prior = T)
}



#multiple subpopulations CV
# cross validate
{
  # do CV
  {

    selected = select_bounded(1000, Nm_values_sim, m_range = c(0,0.8), N_range = c(500,22500))
    cv_out_AF_S = array(NA, dim = c(0,2 ))
    cv_out_AF_S_list = list(); i = 1
    S_loc = round(c(1000, 500, 200, 100))
    for (S in S_loc){
      S = round(S/6)
      SS = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100,S*30,30)), gens_select = c(0,10,20), subpops_select_sim = c(40:41), subpops_select_real = c(60:61), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      cv = cross_validate(do_abc_PLS, SS_sim, selected, n_comps =  min(nc, dim(SS_sim)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
      RMS_cv = RMS_err(cv); means = colMedians(RMS_cv)

      
      cv_out_AF_S_list[[i]] = RMS_cv; i = i+1
      cv_out_AF_S = rbind(cv_out_AF_S,means)
      
    }
    
    
    saveRDS(cv_out_AF_S_list, '../data/cv_out_AF_S_list.RDS')
    cv_out_AF_S_list = readRDS('../data/cv_out_AF_S_list.RDS')
  }
  
  # make plot
  {
    cv_out_AF_S = matrix(unlist(lapply(cv_out_AF_S_list,colMeans)), nc=2, byrow = T)
    
    pdf('../final_figures/CV_AF.pdf', height = 2.7, width = 5, pointsize = 10.5)
    
    
    max_y = 0
    i = 1
    for (x in cv_out_AF_S_list){
      ci = do_boot(x[,2])
      max_y = max(max_y, ci)
      i = i+1
    }
    
    par(mfrow = c(1,1), mai = c(.6,.7,.05,.2))
    plot(NA, xlab = expression(S), ylab = 'Relative Error', ylim = range(0,max_y), xlim = c(100,1000), xaxt = 'n'); axis(side = 1 , at = S_loc)
    
    i = 1
    for (x in cv_out_AF_S_list){
      ci = do_boot(x[,1])
      lines(-5+rep(S_loc[i],2), ci, lw = .5)
      lines(-5+S_loc[i]+c(-4,4), rep(ci[1],2), lw = .5)
      lines(-5+S_loc[i]+c(-4,4), rep(ci[2],2), lw = .5)
      
      i = i+1
    }
    
    points(-5+S_loc, cv_out_AF_S[,1],pch = 21,cex = .8, bg='black')
    
    
    i = 1
    for (x in cv_out_AF_S_list){
      ci = do_boot(x[,2])
      lines(5+rep(S_loc[i],2), ci, lw = .5, col = 'black')
      lines(5+S_loc[i]+c(-4,4), rep(ci[1],2), lw = .5)
      lines(5+S_loc[i]+c(-4,4), rep(ci[2],2), lw = .5)
      i = i+1
    }
    
    points(5+S_loc, cv_out_AF_S[,2],  pch = 21, cex = .8, col = 'black', bg = 'white');
    
    legend('topright',cex = .9,inset = 0.02, legend = c(expression(N[e]), 'm' ), pch = c(21,21), pt.bg = c('black', 'white'))
    graphics.off()
  }

  
  
  
  ## small N
  {
    S_total = 1000
    
    selected = select_bounded(100, Nm_values_sim, m_range = c(0,0.8), N_range = c(500,10000))
    RMSE_out_AF_small = array(NA, dim = c(0,2 ))
    RMSE_out_list_AF_small = list(); i = 1
    S_loc = round(c(1000, 500, 200, 100))
    for (S in S_loc){
      S = round(S/6)
      SS = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100,S*30,30)), gens_select = c(0,10,20), subpops_select_sim = c(40:41), subpops_select_real = c(60:61), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      cv = cross_validate(do_abc_PLS, SS_sim, selected, n_comps =  min(nc, dim(SS_sim)[2]), Nm_values_sim, method = method, tol = tol, sizenet = sizenet, numnet = numnet)
      RMS_cv = RMS_err(cv); means = colMedians(RMS_cv)
      
      
      RMSE_out_list_AF_small[[i]] = RMS_cv; i = i+1
      RMSE_out_AF_small = rbind(RMSE_out_AF_small,means)
      
    }
    
    
    saveRDS(RMSE_out_list_AF_small, '../data/RMSE_out_list_AF_small.RDS')
    RMSE_out_list_AF_small = readRDS('../data/RMSE_out_list_AF_small.RDS')
  }
  
  {
    RMSE_out_AF_small = matrix(unlist(lapply(RMSE_out_list_AF_small,colMeans)), nc=2, byrow = T)
    
    pdf('../final_figures/CV_AF.pdf', height = 4, width = 7, pointsize = 10.5)
    
    
    max_y = 0
    i = 1
    for (x in RMSE_out_list_AF_small){
      ci = do_boot(x[,2])
      max_y = max(max_y, ci)
      ci = do_boot(x[,1])
      max_y = max(max_y, ci)
      i = i+1
    }
    
    par(mfrow = c(1,1))
    plot(NA, xlab = expression(S), ylab = 'Relative Error', ylim = range(0,max_y,RMSE_out_AF_small[,]), xlim = c(100,1000), xaxt = 'n'); axis(side = 1 , at = S_loc)
    
    i = 1
    for (x in RMSE_out_list_AF_small){
      ci = do_boot(x[,1])
      lines(-3+rep(S_loc[i],2), ci, lw = .5)
      i = i+1
    }
    
    points(-3+S_loc, RMSE_out_AF_small[,1],pch = 21,cex = 1, bg='white')
    
    
    i = 1
    for (x in RMSE_out_list_AF_small){
      ci = do_boot(x[,2])
      lines(3+rep(S_loc[i],2), ci, lw = .5, col = 'black')
      i = i+1
    }
    
    points(3+S_loc, RMSE_out_AF_small[,2],  pch = 25, cex = 1, col = 'black', bg = 'white');
    
    legend('topright',cex = .9,inset = 0.02, legend = c('N', 'N (S = Inf)', 'm', 'm (S = Inf)' ), lty = c(NA, 5,NA,3), pch = c(1,NA,6, NA))
    graphics.off()
  }
  
}
##################################
##################################
######## CIs #####################
##################################
##################################

# LOOCV - how often are the true values within the 97.5% CI (coverage)
# width of confidence intervals - LOOCV variance
# make a few plots with confidence limits on them

{

  {
    SS_AF1000 = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, round(1000/2/3)*30, 30)), gens_select = c(0,10,20), subpops_select_sim = 40:42, subpops_select_real = 60:62, assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
    SS_sim_AF1000 = SS_AF1000$sim; SS_real_AF1000 = SS_AF1000$real
    
    SS_AF500 = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, round(500/2/3)*30, 30)), gens_select = c(0,10,20), subpops_select_sim = 40:41, subpops_select_real = 60:61, assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
    SS_sim_AF500 = SS_AF500$sim; SS_real_AF500 = SS_AF500$real
    
    SS_AF200 = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, round(200/2/3)*30, 30)), gens_select = c(0,10,20), subpops_select_sim = 40:42, subpops_select_real = 60:62, assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
    SS_sim_AF200 = SS_AF200$sim; SS_real_AF200 = SS_AF200$real
  }
  
  
  selectedCI = select_bounded(1000, Nm_values_sim, N_range = c(400, 20000), m_range = c(0.01, 0.8))
  CI_results_2 = check_CI_width(do_abc_PLS, SS_sim_AF1000, selectedCI, n_comps =  4, Nm_values_sim, method = 'neuralnet', tol = .5/10, sizenet = 4, numnet = 1)
  
  N_contained = CI_results_2[,1,1] > CI_results_2[,2,1] & CI_results_2[,1,1] < CI_results_2[,3,1]
  m_contained = CI_results_2[f,1,2] > CI_results_2[f,2,2] & CI_results_2[f,1,2] < CI_results_2[f,3,2]

  mean(N_contained)
  mean(m_contained)
  
  #plot
  {
  plot(NA, NA, xlim = range(CI_results_2[,1,1]), ylim = c(0,4))
  for (i in 1:length(selectedCI)){
    lines(rep(CI_results_2[i,1,1],2), c(CI_results_2[i,2,1]/CI_results_2[i,1,1], CI_results_2[i,3,1]/CI_results_2[i,1,1]))
  }
  abline(h=2)
  
  plot(NA, NA, xlim = range(CI_results_2[,1,2]), ylim = c(0,5))
  for (i in 1:length(selectedCI)){
    lines(rep(CI_results_2[i,1,2],2), c(CI_results_2[i,2,2]/CI_results_2[i,1,2], CI_results_2[i,3,2]/CI_results_2[i,1,2]))
  }
  abline (v = .025)
  }
  
  f = CI_results_2[,1,1] < 10000 & CI_results_2[,1,1] > 1000 #& CI_results_2[,1,2]  > 0.03 & CI_results_2[,2,2] < .4

  x=1
  # fully contained
  mean(CI_results_2[f,3,x] < 2 * CI_results_2[f,1,x] & CI_results_2[f,2,x] > .5 * CI_results_2[f,1,x])
  # contained low; m > 0.15
  mean(CI_results_2[f,2,x] > (.5 * CI_results_2[f,1,x]))
  # contained high
  mean(CI_results_2[f,3,x] < 2 * CI_results_2[f,1,x] )
  # mean width
  mean(log2(CI_results_2[f,3,1] / CI_results_2[f,2,1]))
  
  
  
  # check: all < 0.3, all < 0.03, interval [0.03,0.3], interval [0.015,0.025]
  # discard min 0 when checking for width only
  
  f =   CI_results_2[,1,2]  < 0.3  & CI_results_2[,2,2] !=0  #
  f =   CI_results_2[,1,2]  < 0.03  & CI_results_2[,2,2] !=0  #
  f =   CI_results_2[,1,2]  > 0.03 & CI_results_2[,1,2] < 0.3  & CI_results_2[,2,2] !=0  #
  f =   CI_results_2[,1,2]  > 0.015 & CI_results_2[,1,2] < 0.025 # & CI_results_2[,2,2] !=0  #
  colMeans(CI_results_2[f,,2])
  x=2
  # fully contained
  mean(CI_results_2[f,3,x] < 2 * CI_results_2[f,1,x] & CI_results_2[f,2,x] > .5 * CI_results_2[f,1,x])
  # contained low; m > 0.15
  mean(CI_results_2[f,2,x] > (.5 * CI_results_2[f,1,x]))
  # contained high
  mean(CI_results_2[f,3,x] < 2 * CI_results_2[f,1,x] )
  #width
  mean(log2(CI_results_2[f,3,x] / CI_results_2[f,2,x]))
}



##################################
##################################
### accepted + hist ##############
##################################
##################################

{
  
  # SS
  {
  SS_AF1000 = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, round(1000/2/3)*30, 30)), gens_select = c(0,10,20), subpops_select_sim = 40:42, subpops_select_real = 60:62, assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_AF1000 = SS_AF1000$sim; SS_real_AF1000 = SS_AF1000$real
  
  SS_AF1000.single = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, round(1000/3)*30, 30)), gens_select = c(0,10,20), subpops_select_sim = 41, subpops_select_real = 61, assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_AF1000.single = SS_AF1000.single$sim; SS_real_AF1000.single = SS_AF1000.single$real

  SS_AF500 = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, round(500/2/3)*30, 30)), gens_select = c(0,10,20), subpops_select_sim = 40:41, subpops_select_real = 60:61, assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_AF500 = SS_AF500$sim; SS_real_AF500 = SS_AF500$real
  
  SS_AF200 = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, round(200/2/3)*30, 30)), gens_select = c(0,10,20), subpops_select_sim = 40:42, subpops_select_real = 60:62, assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_AF200 = SS_AF200$sim; SS_real_AF200 = SS_AF200$real

  SS_single = make_SS(sim_out, real_out_list, S_select = Inf, gens_select = c(0,10,20), subpops_select_sim = 41, subpops_select_real = 61, assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_single = SS_single$sim; SS_real_single = SS_single$real
  
  SS_double = make_SS(sim_out, real_out_list, S_select = Inf, gens_select = c(0,10,20), subpops_select_sim = 41:42, subpops_select_real = 61:62, assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_double = SS_double$sim; SS_real_double = SS_double$real
  
  SS_ind = make_SS(sim_out, real_out_list, S_select = round(1000/9), gens_select = c(0,10,20), subpops_select_sim = 40:42, subpops_select_real = 60:62, assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_ind = SS_ind$sim; SS_real_ind = SS_ind$real
  }
  
    
  n1 = 58; print(round(Nm_values_real[n1,], 4))
  n2 = 91; print(round(Nm_values_real[n2,], 4))
  
  o1.AF <- do_abc_PLS(SS_sim_AF1000, SS_real_AF1000[[4]][n1,], PLS_use = 4, Nm_values_sim, method = 'neuralnet', tol = .1, numnet = 1, sizenet = 4)
  o2.AF <- do_abc_PLS(SS_sim_AF1000, SS_real_AF1000[[18]][n2,], PLS_use = 4, Nm_values_sim, method = 'neuralnet', tol = .1, numnet = 1, sizenet = 4)

  
  
  o1.AFsingle <- do_abc_PLS(SS_sim_AF1000.single, SS_real_AF1000.single[[4]][n1,], PLS_use = 2, Nm_values_sim, method = 'neuralnet', tol = .1, numnet = 1, sizenet = 4)
  o2.AFsingle <- do_abc_PLS(SS_sim_AF1000.single, SS_real_AF1000.single[[4]][n2,], PLS_use = 2, Nm_values_sim, method = 'neuralnet', tol = .1, numnet = 1, sizenet = 4)
  
  o1.single <- do_abc_PLS(SS_sim_ind, SS_real_ind[[21]][n1,], PLS_use = 3, Nm_values_sim, method = 'neuralnet', tol = .1, numnet = 1, sizenet = 4)
  o2.single <- do_abc_PLS(SS_sim_ind, SS_real_ind[[8]][n2,], PLS_use = 3, Nm_values_sim, method = 'neuralnet', tol = .1, numnet = 1, sizenet = 4)
  summary(o1.single)
  summary(o2.single)
  
  
  
  # colienarity
  {
  pdf('../final_figures/colinearity.pdf', height = 5, width = 7, pointsize=10.5)
  par(mai = c(1,1,.1,.1))
  plot(Nm_values_sim[,2], rowSds(SS_sim_single)/ rowMeans(SS_sim_single), pch = 16, cex = 0.4, ylab = bquote('CV('*bold(S)*')'), xlab = 'm')
  graphics.off()
  }
  
  # posteriors
  {
  pdf('../final_figures/posterior_hists.pdf', height = 6, width = 7, pointsize = 10.5)
  par(mfrow = c(3,2), mai = c(.6, .7,0.11,0.11), oma = c(1.5,0,0,0))
  hist(o2.AF$adj.values[,1],ylim = 2*c(0,250),xaxs = 'i', yaxs = 'i', bty = 'l', xlim = c(200,20000),breaks = seq(000,25000,300), xlab = expression(N[e]), col = 'coral', border = 'white', main = NULL, xaxt = 'n');axis(1, c(200,5000,10000,15000,20000)); abline(v =  Nm_values_real[n2,1], col = 'black', lw = 1, lty = 2);# abline(v = mean(o$adj.values[,1]), col = 'black')
  hist(o1.AF$adj.values[,1], xlim = c(200,20000),breaks = seq(000,25000,300), col = 'cyan3', border = 'white', add = T); abline(v =  Nm_values_real[n1,1], col = 'black', lw = 1, lty = 2);# abline(v = mean(o$adj.values[,1]), col = 'black')
  axis(1, c(200,5000,10000,15000,20000)); axis(2)
  fig_label('A', cex = 2)
  
  hist(o2.AF$adj.values[,2], ylim = 3*c(0,160), xaxs = 'i', yaxs = 'i', bty = 'l', xlim = c(0,.8),breaks = seq(0,.8,0.0125), xlab = expression(m), col = 'coral', border = 'white', main = NULL, xaxt = 'n'); abline(v =  Nm_values_real[n2,2], col = 'black', lw = 1, lty = 2);# abline(v = mean(o$adj.values[,1]), col = 'black')
  hist(o1.AF$adj.values[,2], xlim = c(0,.8),breaks = seq(0,.8,0.0125), col = 'cyan3', border = 'white', add = T); abline(v =  Nm_values_real[n1,2], col = 'black', lw = 1, lty = 2);# abline(v = mean(o$adj.values[,1]), col = 'black')
  axis(1, seq(0,0.8,0.1)); axis(2)
  fig_label('B', cex = 2)
  
  
  hist(o2.single$adj.values[,1], ylim = 2*c(0,100),xaxs = 'i', yaxs = 'i', bty = 'l', xlim = c(200,20000),breaks = seq(0,20000,250), xlab = expression(N[e]), col = 'coral', border = 'white', main = NULL, xaxt = 'n'); abline(v =  Nm_values_real[n2,1], col = 'black', lw = 1, lty = 2);# abline(v = mean(o$adj.values[,1]), col = 'black')
  hist(o1.single$adj.values[,1][in_range(o1.single$adj.values[,1], c(400,20000))], xlim = c(200,20000),breaks = seq(0,20000,250), col = 'cyan3', border = 'white', add = T); abline(v =  Nm_values_real[n1,1], col = 'black', lw = 1, lty = 2);# abline(v = mean(o$adj.values[,1]), col = 'black')
  axis(1, c(200,5000,10000,15000,20000)); axis(2)
  fig_label('C', cex = 2)
  
  
  hist(o2.single$adj.values[,2], ylim = 2*c(0,100),xaxs = 'i', yaxs = 'i', bty = 'l', xlim = c(0,.8),breaks = seq(0,.8,0.0125), xlab = expression(m), col = 'coral', border = 'white', main = NULL, xaxt = 'n');axis(1, seq(0,0.8,0.1)); abline(v =  Nm_values_real[n2,2], col = 'black', lw = 1, lty = 2);# abline(v = mean(o$adj.values[,1]), col = 'black')
  hist(o1.single$adj.values[,2], xlim = c(0,.8),breaks = seq(0,.8,0.0125), col = 'cyan3', border = 'white', add = T); abline(v =  Nm_values_real[n1,2], col = 'black', lw = 1, lty = 2);# abline(v = mean(o$adj.values[,1]), col = 'black')
  axis(1, seq(0,0.8,0.1)); axis(2)
  fig_label('D', cex = 2)
  
  
  hist(o2.AFsingle$adj.values[,1], ylim = 2*c(0,80),xaxs = 'i', yaxs = 'i', bty = 'l', xlim = c(200,20000),breaks = seq(0,20000,250), xlab = expression(N[e]), col = 'coral', border = 'white', main = NULL, xaxt = 'n'); abline(v =  Nm_values_real[n2,1], col = 'black', lw = 1, lty = 2);# abline(v = mean(o$adj.values[,1]), col = 'black')
  hist(o1.AFsingle$adj.values[,1], xlim = c(200,20000),breaks = seq(0,20000,250), col = 'cyan3', border = 'white', add = T); abline(v =  Nm_values_real[n1,1], col = 'black', lw = 1, lty = 2);# abline(v = mean(o$adj.values[,1]), col = 'black')
  axis(1, c(200,5000,10000,15000,20000)); axis(2)
  fig_label('E', cex = 2)
  
  
  hist(o2.AFsingle$adj.values[,2], ylim = 2*c(0,200),xaxs = 'i', yaxs = 'i', bty = 'l', xlim = c(0,.8),breaks = seq(0,.8,0.0125), xlab = expression(m), col = 'coral', border = 'white', main = NULL, xaxt = 'n');axis(1, seq(0,0.8,0.1)); abline(v =  Nm_values_real[n2,2], col = 'black', lw = 1, lty = 2);# abline(v = mean(o$adj.values[,1]), col = 'black')
  hist(o1.AFsingle$adj.values[,2], xlim = c(0,.8),breaks = seq(0,.8,0.0125), col = 'cyan3', border = 'white', add = T); abline(v =  Nm_values_real[n1,2], col = 'black', lw = 1, lty = 2);# abline(v = mean(o$adj.values[,1]), col = 'black')
  axis(1, seq(0,0.8,0.1)); axis(2)
  fig_label('F', cex = 2)
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", c(as.expression(bquote('N'['e']*' = 10000; m = 0.055' )), as.expression(bquote('N'['e']*' = 2000;   m = 0.4' ))), bty = "n", fill=c("cyan1", "coral"))
  graphics.off()
  
  }
  
  
  # Nm isolcine - estimation 
  {
  pdf('../final_figures/Nm_isocline.pdf', width = 7, height = 4, pointsize = 10.5)
    
  n = 89; print(round(Nm_values_real[n,], 4))
  o1 <- do_abc_PLS(SS_sim_single, SS_real_single[[9]][n,], PLS_use = 3, Nm_values_sim, method = 'neuralnet', tol = .02, numnet = 10, sizenet = 10)
  n = 117; print(round(Nm_values_real[n,], 4))
  o2 <- do_abc_PLS(SS_sim_single, SS_real_single[[8]][n,], PLS_use = 3, Nm_values_sim, method = 'neuralnet', tol = .02, numnet = 10, sizenet = 10)
  
  par(mfrow = c(1,2), mai = c(.7,.7,.2,.2), oma = c(4,.5,.5,0))
  
  # plot o1
  {
    plot(Nm_values_sim[,1], Nm_values_sim[,2], xlab = as.expression(bquote('N'['e'])), ylab = 'm', pch = 16, cex = .35, col = 'darkslategrey', xlim = c(200,20000), ylim = c(0,0.8), xaxs = 'i', yaxs = 'i')#, pch = 16)
    fig_label('A', cex = 2)
    
    points(o1$unadj.values[,1], o1$unadj.values[,2], pch = 16, col = 'red', cex = .7)
    points(o1$adj.values[,1], o1$adj.values[,2], pch = 16, col = 'deepskyblue', cex = .7)
    points(Nm_values_real[n,1],Nm_values_real[n,2], pch = 8, cex = 2, col = 'black', lw = 1.5)

    for (x in seq(200,10000,1000)){
      lines(seq(1,20000), x/seq(1,20000), col = 'darkslategrey', lw = 1, lty = 1)
    }

  }
  
  # plot 02
  {
    plot(Nm_values_sim[,1], Nm_values_sim[,2], xlab = as.expression(bquote('N'['e'])), ylab = 'm', pch = 16, cex = .35, col = 'darkslategrey', xlim = c(200,20000), ylim = c(0,0.8), xaxs = 'i', yaxs = 'i')#, pch = 16)
    fig_label('B', cex = 2)
    points(o2$unadj.values[,1], o2$unadj.values[,2], pch = 16, col = 'red', cex = .7, lw = 2)
    points(o2$adj.values[,1], o2$adj.values[,2], pch = 16, col = 'deepskyblue1', cex = .7, lw = 2)
    points(Nm_values_real[n,1],Nm_values_real[n,2], pch =8, cex = 2, col = 'black', lw = 1.5)

    
    for (x in seq(200,10000,1000)){
      lines(seq(1,20000), x/seq(1,20000), col = 'darkslategrey', lw = 1, lty = 1)
    }

  }
  
  # legend 
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend('bottom',inset = 0.04, col = c('grey','red', 'blue', 'black', 'black'), pch = c(16, 16, 16, 8, NA), lty = c(NA, NA, NA, NA, 1), legend = c('rejected','accepted (unadj.)', 'accepted (adj.)', 'target', as.expression(bquote('N'['e']*'m isocline' )) ),  nc = 3, cex = .9)
  
  graphics.off()
}
  
  
  # PCA
  
  {
    
    
  accept_idx1 = abs(rowProds(Nm_values_sim)-5000) < 200 & Nm_values_sim[,2]>0.33
  accept_idx2_small = abs(rowProds(Nm_values_sim)-1000) < 100 & Nm_values_sim[,2]<0.33
  accept_idx2_large = abs(rowProds(Nm_values_sim)-1000) < 100 & Nm_values_sim[,2]>0.33
  
  # colours
  {
  # Define colour pallete
  pal = colorRampPalette(c("blue", "red"))
  # Rank variable for colour assignment
  ord = findInterval(Nm_values_sim[accept_idx2_small,2], sort(Nm_values_sim[accept_idx2_small,2]))
  col_arr = pal(length(ord))[ord]; len = length(col_arr)
  
  pal2 = colorRampPalette(c("green", "dark green"))
  ord2 = findInterval(Nm_values_sim[accept_idx2_large,2], sort(Nm_values_sim[accept_idx2_large,2]))
  col_arr2 = pal2(length(ord2))[ord2]
  }
  
  {
    pdf('../final_figures/Nm_PCA.pdf', width = 7, height = 4, pointsize = 11)
    layout(matrix(c(1,2,4,1,3,4), nc=3, byrow = T), width = c(2,2,.5),height = c(1,1,1))
    par(mai = c(.6,.5,.3,.3))
    
    
    # plot sim values
    {
    plot(Nm_values_sim[,1], Nm_values_sim[,2], xlab = as.expression(bquote('N'['e'])), ylab = 'm', pch = 16, cex = .35, col = 'grey', xlim = c(200,20000), ylim = c(0,0.8), xaxs = 'i', yaxs = 'i', pty = 'm', xaxt = 'n'); axis(1, c(1000,5000,10000,150000,20000))#, pch = 16)
      for (x in c(100,seq(1000,10000,1000))){
        lines(seq(1,20000), x/seq(1,20000), col = 'darkslategrey', lw = 1, lty = 1)
      }
    points(Nm_values_sim[accept_idx2_small,1],Nm_values_sim[accept_idx2_small,2], pch = 16, col = col_arr, cex = 1)
    points(Nm_values_sim[accept_idx2_large,1],Nm_values_sim[accept_idx2_large,2], pch = 16, col = col_arr2, cex = 1)
    

    fig_label('A',cex = 2)
    }
    
    # remove negative correlations
    SS_sim_double_pc = SS_sim_double
    for (i in 1:dim(SS_sim_double_pc)[2]){SS_sim_double_pc[SS_sim_double_pc[,i]<0,i] = mean(SS_sim_double_pc[,i])}
    
    # pca single subpopulation
    pc <- prcomp(log(SS_sim_single[,]), scale = T, center= T)
    
    par(mai = c(0.2,0.2,0.1,0.1))
    plot(pc$x[,1], pc$x[,2], xlab = '', ylab = '', pch = 16, cex = .3, pty = 'm', xaxt = 'n', yaxt = 'n', col = 'grey', xlim = range(pc$x[,1])-c(0,6),  ylim =range(pc$x[,2]))#, xaxs = 'i', yaxs = 'i')#, xlim = c(-6,12), ylim = c(-.8,.8))
    title(ylab="PC2", line=0.1, cex.lab=1.2)
    title(xlab="PC1", line=0.1, cex.lab=1.2)
    text(-3.8,.58, bquote(n[s]*"= 1"), border = T)
    points(pc$x[accept_idx2_small,1], pc$x[accept_idx2_small,2], col = col_arr, pch = 16, cex = .8)
    points(pc$x[accept_idx2_large,1], pc$x[accept_idx2_large,2], col = col_arr2, pch = 16, cex = .8)
    fig_label('B',cex = 2)
    
    
    # pca double subpopualtions
    pc <- prcomp(log(SS_sim_double_pc), scale = T, center= T)

    par(mai = c(0.2,0.2,0.1,0.1))
    plot(pc$x[,1], pc$x[,2], xlab = '', ylab = '', pch = 16, cex = .3, pty = 'm', xaxt = 'n', yaxt = 'n', col = 'grey', xlim = range(pc$x[,1])-c(0,6), ylim =range(pc$x[,2]))#, xaxs = 'i', yaxs = 'i')#, xlim = c(-6,12), ylim = c(-.8,.8))
    title(ylab="PC2", line=0, cex.lab=1.2)
    title(xlab="PC1", line=0, cex.lab=1.2)
    text(-4,5.5, bquote(n[s]*'= 2'))
    points(pc$x[accept_idx2_small,1], pc$x[accept_idx2_small,2], col = col_arr, pch = 16, cex = .8)
    points(pc$x[accept_idx2_large,1], pc$x[accept_idx2_large,2], col = col_arr2, pch = 16, cex = .8)
    fig_label('C', cex = 2)
  
    
    #legend
    {
    legend_image <- as.raster(matrix(rev(c(pal(34),'white', pal2(60))), ncol=1),)
    par(mai = c(0,0,0,0))
    plot(NA,NA, xlim = c(0,5), ylim = c(0,5), bty = 'n', xaxt = 'n', yaxt = 'n', ylab = NA, xlab = NA)
  
    text(x=2, y = c(seq(1,3.95, l =9)), labels = c(paste0('    ',seq(0,0.8,l=9))))
    text(x=1.2, y = c(4.2), labels = bquote(bold('m')))
    
    rasterImage(legend_image, -1, 1, 1,4, interpolate = F)
    }
    graphics.off()
    
    }
    pc$sdev^2/sum(pc$sdev^2) * 100
  }

}



##################################
##################################
### Nm and Inf ##############
##################################
##################################
  
# Estimating Nm
{
  
  # do ABC
  {
    o_Nm = list()
    o_Nm_full = list()
    
    i=1
    for (S in c(100, 200, 500)){
      SS = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, 30*round(S / 3), 30)), gens_select = c(0,10,20), subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
      SS_sim = SS$sim; SS_real = SS$real
      
      oot = test_specific_Nm(do_abc_PLS_Nm, SS_sim, SS_real, n_comps = min(nc, dim(SS_sim)[2]), rowProds(Nm_values_sim), rowProds(Nm_values_real), idx = c(83,114), tol = tol, method = method, sizenet = sizenet, numnet = numnet )
      
      o_Nm_full[[i]] = oot
      i=i+1
    }
    
    saveRDS(o_Nm_full, '../data/o_Nm_full.RDS')
    o_Nm_full = readRDS('../data/o_Nm_full.RDS')
    
  }
  
  # make plot
  {
    i=1
    for (l in o_Nm_full){
      o_Nm[[i]] = apply(l[,,,], c(1,2), trim_mean,2) 
      i=i+1
    }
    pdf('../final_figures/Nm_estimation.pdf', width = 13, height = 6, pointsize = 10.5)
    plot_compare_m_neat_Nm(o_Nm, legend = c(100,200,500))
    graphics.off()
  }
}


random_fib = function(n){
  vals = rep(NA, n)
  vals[1:2] = c(1,1)
  
  for (i in 3:n){
    sign = sample(c(1,-1),1)
    vals[i] = vals[i-1] + sign*vals[i-2]
  }
  return(vals)
}


outs = array(NA, dim = c(1000,1000))
for (i in 1:1000){
  outs[i,] = random_fib(1000)
}


hist(abs(outs[,1000])^(1/1000), breaks = 100); abline(v=1.1319882487943, col = 'red')
mean(abs(outs[,1000])^(1/1000))

plot(1:100,colMeans(abs(outs)))
lines(1:100, 1.1319882487943^(1:100))

fibseq = random_fib(1e6)
plot(1:1000,fibseq)

# inf
{
  o_inf = list()
  o_inf_full = list()
  
  i=1
  for (S in c(Inf)){
    SS = make_SS(sim_out, real_out_list, S_select = 7750, gens_select = c(0,5,10,20), subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
    SS_sim = SS$sim; SS_real = SS$real
    
    oot = test_specific(do_abc_PLS, SS_sim, SS_real, n_comps = min(nc, dim(SS_sim)[2]), (Nm_values_sim), (Nm_values_real), idx = c(83,114), tol = tol, method = method, sizenet = sizenet, numnet = numnet )
    
    o_inf_full[[i]] = oot
    i=i+1
  }
  
  saveRDS(o_inf_full, '../data/o_inf_full.RDS')
  o_inf_full = readRDS('../data/o_inf_full.RDS')
  
  o_inf[[1]] = apply(o_inf_full[[1]], c(1,2,3), mean)
  
  pdf('../final_figures/S_inf_single_subpop.pdf', width = 7, height = 7, pointsize = 10.5)
  plot_compare_m_neat(o_inf, legend = c(NULL))
  graphics.off()
  saveRDS(o_Nm_full, '../data/o_single_2000_full.RDS')
}



##### How good would random guessing from the prior be #####
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

m_prior_draw_.05 = function(n){
  v = c()
  while (length(v) < n){
    temp = m_prior_draw(1)
    if (in_range(temp, c(0,0.05))){
      v = c (v, temp)
    }
  }
  return(v)
}

m_prior_draw(10000)
{
  ms_ests = cbind(m_prior_draw_no0(1e4), (m_prior_draw_no0(1e4)))
  Ns_ests = cbind(N_prior_draw(1e4), (N_prior_draw(1e4)))
  
  results_guess = abind(Ns_ests, ms_ests, along = 3)
  
  err_guess = RMS_err(results_guess)
  colMeans(err_guess)
}




##################################
##################################
###   Estiamting Nm ##############
##################################
##################################




####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################



#####################################
#####################################
######## LD vs AF (old) #############
#####################################
#####################################

{

# RMSE with varying S for some specific across space, time sampling strategies (eg we know that LD info is important when sampling from multiple subpopulations)
{
plot(Nm_values_sim[,2], summarise_rirj(sim_out$rirj_sim,S = 5000, subpops = c(40,41,42), timepoints = sample_gens, edge_length = 9,focal = 41)[1,,])

# select data - LD vs AF Inf sample size
SS_LD = make_SS(sim_out, real_out_list, S_select = Inf, gens_select = gens_select, subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = F, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
SS_sim_LD = SS_LD$sim; SS_real_LD = SS_LD$real

SS_AF = make_SS(sim_out, real_out_list, S_select = Inf, gens_select = gens_select, subpops_select_sim = c(41), subpops_select_real = c(61), assume_equal = F, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
SS_sim_AF = SS_AF$sim; SS_real_AF = SS_AF$real



SS_LD = make_SS(sim_out, real_out_list, S_select = Inf, gens_select = gens_select, subpops_select_sim = c(32,40,41,42,50), subpops_select_real = c(50,60,61,62,72), assume_equal = F, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_LD = SS_LD$sim; SS_real_LD = SS_LD$real

SS_AF = make_SS(sim_out, real_out_list, S_select = Inf, gens_select = gens_select, subpops_select_sim = c(32,40,41,42,50), subpops_select_real = c(50,60,61,62,72), assume_equal = F, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_AF = SS_AF$sim; SS_real_AF = SS_AF$real


# vary N

a_N = test_specific(do_abc_PLS, SS_sim_LD, SS_real_LD, n_comps = 5, Nm_values_sim, Nm_values_real, idx = c(1,37) )
b_N = test_specific(do_abc_PLS, SS_sim_AF, SS_real_AF, n_comps = 5, Nm_values_sim, Nm_values_real, idx = c(1,37) )

a_N_mean = apply(a_N, c(1,2,3), mean)
b_N_mean = apply(b_N, c(1,2,3), mean)

plot_compare_N(list(a_N_mean,b_N_mean), legend = c('LD + AF', 'AF only'))

# vary m

a_m = test_specific(do_abc_PLS, SS_sim_LD, SS_real_LD, n_comps = 3, Nm_values_sim, Nm_values_real, idx = c(38,52), tol = .03 )
b_m = test_specific(do_abc_PLS, SS_sim_AF, SS_real_AF, n_comps = 3, Nm_values_sim, Nm_values_real, idx = c(38,52), tol = .03 )

a_m_mean = apply(a_m, c(1,2,3), trim_mean)
b_m_mean = apply(b_m, c(1,2,3), trim_mean)

plot_compare_m(list(a_m_mean,b_m_mean), legend = c('LD + AF', 'AF only'))

# small m

a_small_m = test_specific(do_abc_PLS, SS_sim_LD, SS_real_LD, n_comps = 5, Nm_values_sim, Nm_values_real, idx = c(95,119) )
b_small_m = test_specific(do_abc_PLS, SS_sim_AF, SS_real_AF, n_comps = 5, Nm_values_sim, Nm_values_real, idx = c(95,119) )

a_small_m_mean = apply(a_small_m, c(1,2,3), mean)
b_small_m_mean = apply(b_small_m, c(1,2,3), mean)

plot_compare_m(list(a_small_m,b_small_m), legend = c('1 sub', '4 sub'))
}


# LD, AF all sample sizes
{
SS_sim_list = list()
SS_real_list = list()
i=1
for (S in sort(S_vec)){
  SS_LD = make_SS(sim_out, real_out_list, S_select = S, gens_select = c(0,20), subpops_select_sim = c(40,41,42), subpops_select_real = c(60,61,62), assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
  #SS_LD = make_SS(sim_out, real_out_list, S_select = S, gens_select = gens_select, subpops_select_sim = c(31,32,33,40,41,42,49,50,51), subpops_select_real = c(49,50,51,60,61,62,71,72,73), assume_equal = F, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
  #SS_LD = make_SS(sim_out, real_out_list, S_select = S, gens_select = gens_select, subpops_select_sim = c(41,42), subpops_select_real = c(61,62), assume_equal = F, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
  
  SS_sim_LD = SS_LD$sim[,]; SS_real_LD = SS_LD$real #lapply(SS_LD$real, asub, 5:8, 2)
  SS_sim_list[[i]] = SS_sim_LD
  SS_real_list[[i]] = SS_real_LD
  i=i+1
}
for (S in sort(S_vec)){
  SS_AF = make_SS(sim_out, real_out_list, S_select = S, gens_select = c(0,20), subpops_select_sim = c(40,41,42), subpops_select_real = c(60,61,62), assume_equal = T, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
  #SS_AF = make_SS(sim_out, real_out_list, S_select = S, gens_select = gens_select, subpops_select_sim = c(31,32,33,40,41,42,49,50,51), subpops_select_real = c(49,50,51,60,61,62,71,72,73), assume_equal = F, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
  #SS_AF = make_SS(sim_out, real_out_list, S_select = S, gens_select = gens_select, subpops_select_sim = c(41,42), subpops_select_real = c(61,62), assume_equal = F, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
  
  SS_sim_AF = SS_AF$sim; SS_real_AF = SS_AF$real
  SS_sim_list[[i]] = SS_sim_AF
  SS_real_list[[i]] = SS_real_AF
  i=i+1
}

}


# cross validation
{
  selected = select_bounded(100, Nm_values_sim, m_range = c(0,0.8), N_range = c(400,22500))
  results_list_cv = list(); i=1
  for (i in 1:length(SS_sim_list)){
    SS_sim_i = SS_sim_list[[i]]
    res = cross_validate(do_abc_PLS, SS_sim_i, selected, n_comps =  min(3, dim(SS_sim_i)[2]), Nm_values_sim, method = 'neuralnet', tol = .1)
    results_list_cv[[i]] =  res
    
  }
  
  
  mytest = cross_validate(do_abc_trunc, SS_sim_list[[1]], selected, n_comps =  4, Nm_values_sim, method = 'neuralnet', tol = .1)
  colMeans(RMS_err(mytest))
  colMeans(RMS_err(results_list_cv[[13]]))
  
  RMSE_N_LD = c()
  RMSE_m_LD = c()
  RMSE_N_AF = c()
  RMSE_m_AF = c()
  
  for (i in 1:length(S_vec)){
    o <- colMeans(RMS_err(results_list_cv[[i]]))
    RMSE_N_LD = c(RMSE_N_LD, o[1])
    RMSE_m_LD = c(RMSE_m_LD, o[2])
  }
  
  for (i in (length(S_vec)+1):(2*length(S_vec))){
    o <- colMeans(RMS_err(results_list_cv[[i]]))
    RMSE_N_AF = c(RMSE_N_AF, o[1])
    RMSE_m_AF = c(RMSE_m_AF, o[2])
  }
  
  par(mfrow = c(1,2))
  plot((sort(S_vec)),RMSE_N_AF , type = 'l', ylim = range(c(RMSE_N_AF, RMSE_N_LD)), xlab = expression(log[10](S)), ylab = 'RMSE (N)', xaxt = "n"); abline( h = RMSE_N_AF[length(RMSE_N_AF)], lty = 2)
  axis(1, at= 1:10/2, labels =1:10/2)
  lines((sort(S_vec)),RMSE_N_LD, col = 'red' ); abline( h = RMSE_N_LD[length(RMSE_N_LD)], lty = 2 , col = 'red')
  
  #par(mfrow = c(1,1))
  plot((sort(S_vec)),RMSE_m_AF , type = 'l', ylim = range(c(RMSE_m_AF, RMSE_m_LD)), xlab = expression(log[10](S)), ylab = 'RMSE (m)', xaxt = "n"); abline( h = RMSE_m_AF[length(RMSE_m_AF)], lty = 2)
  axis(1, at= 1:10/2, labels =1:10/2)
  lines((sort(S_vec)),RMSE_m_LD, col = 'red' ); abline( h = RMSE_m_LD[length(RMSE_m_LD)], lty = 2 , col = 'red')
  legend('topright', legend = c('AF only', 'AF + LD', 'AF, S = Inf', 'AF + LD, S = Inf'), col = c('black', 'red', 'black', 'red'), lty = c(1,1,2,2))
  
  
  par(mfrow = c(1,2))
  plot((sort(S_vec)[-11]),RMSE_N_AF[-11] , type = 'l', ylim = range(c(RMSE_N_AF, RMSE_N_LD)), xlab = expression(log[10](S)), ylab = 'RMSE (N)', xaxt = "n"); abline( h = RMSE_N_AF[length(RMSE_N_AF)], lty = 2)
  axis(1, at= 1:10/2, labels =1:10/2)
  lines((sort(S_vec)[-11]),RMSE_N_LD[-11], col = 'red' ); abline( h = RMSE_N_LD[length(RMSE_N_LD)], lty = 2 , col = 'red')
  
  #par(mfrow = c(1,1))
  plot((sort(S_vec)[-11]),RMSE_m_AF[-11] , type = 'l', ylim = range(c(RMSE_m_AF, RMSE_m_LD)), xlab = expression(log[10](S)), ylab = 'RMSE (m)', xaxt = "n"); abline( h = RMSE_m_AF[length(RMSE_m_AF)], lty = 2)
  axis(1, at= 1:10/2, labels =1:10/2)
  lines((sort(S_vec)[-11]),RMSE_m_LD[-11], col = 'red' ); abline( h = RMSE_m_LD[length(RMSE_m_LD)], lty = 2 , col = 'red')
  legend('topright', legend = c('AF only', 'AF + LD', 'AF, S = Inf', 'AF + LD, S = Inf'), col = c('black', 'red', 'black', 'red'), lty = c(1,1,2,2))
}


{
  SS_sim_list = list()
  SS_real_list = list()
  c(S_vec_list[[1]], S_vec_list[[2]])
  i=1
  for (S in S_vec){
    SS_LD = make_SS(sim_out, real_out_list, S_select = S, gens_select = c(0,10,20), subpops_select_sim = c(41,42), subpops_select_real = c(61,62), assume_equal = F, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
    SS_AF = make_SS(sim_out, real_out_list, S_select = S, gens_select = c(0,10,20), subpops_select_sim = c(41,42), subpops_select_real = c(61,62), assume_equal = F, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
    
    #SS_LD = make_SS(sim_out, real_out_list, S_select = S, gens_select = gens_select, subpops_select_sim = c(31,32,33,40,41,42,49,50,51), subpops_select_real = c(49,50,51,60,61,62,71,72,73), assume_equal = F, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
    #SS_LD = make_SS(sim_out, real_out_list, S_select = S, gens_select = gens_select, subpops_select_sim = c(41,42), subpops_select_real = c(61,62), assume_equal = F, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
    
    SS_sim_LD = SS_LD$sim[,]; SS_real_LD = SS_LD$real #lapply(SS_LD$real, asub, 5:8, 2)
    SS_sim_list[[i]] = SS_sim_LD
    SS_real_list[[i]] = SS_real_LD
    i=i+1
  }
  for (S in sort(S_vec)){
    SS_AF = make_SS(sim_out, real_out_list, S_select = S, gens_select = c(0,20), subpops_select_sim = c(41,42), subpops_select_real = c(61,62), assume_equal = F, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
    #SS_AF = make_SS(sim_out, real_out_list, S_select = S, gens_select = gens_select, subpops_select_sim = c(31,32,33,40,41,42,49,50,51), subpops_select_real = c(49,50,51,60,61,62,71,72,73), assume_equal = F, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
    #SS_AF = make_SS(sim_out, real_out_list, S_select = S, gens_select = gens_select, subpops_select_sim = c(41,42), subpops_select_real = c(61,62), assume_equal = F, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11)
    
    SS_sim_AF = SS_AF$sim; SS_real_AF = SS_AF$real
    SS_sim_list[[i]] = SS_sim_AF
    SS_real_list[[i]] = SS_real_AF
    i=i+1
  }
}

# against m
{
results_list = list(); i=1
for (i in 1:length(SS_sim_list)){
  SS_sim_i = SS_sim_list[[i]]
  SS_real_i = SS_real_list[[i]]
  
  res_full = test_specific(do_abc_PLS, SS_sim_i, SS_real_i, n_comps = 3, Nm_values_sim, Nm_values_real, idx = c(38,52), tol = 0.03 )
  res = apply(res_full, c(1,2,3), mean)
  results_list[[i]] =  res
}


leg = c(paste0('LD', sort(S_vec)), paste0('AF', sort(S_vec)))
colfunc_AF <- colorRampPalette(c("red", "green"))
cols = c(rep('black', length(S_vec)), colfunc_AF(length(S_vec)))
ltys = c(1:length(S_vec), rep(1,length(S_vec)))
plot_compare_m(results_list, legend = leg, cols = cols, ltys = ltys)

}




x=log(sort(S_vec)[-c(12)])
y=RMSE_N_AF[-c(12)]
M = lm(y~poly(x,3))
plot(seq(6,12), predict(M,data.frame(x=seq(6,12) )), type = 'l')
points(x,y)
x2=log(sort(S_vec)[-c(12)])
y2=RMSE_N_LD[-c(12)]
M2 = lm(y2~poly(x2,3))
lines(seq(6,12), predict(M2,data.frame(x2=seq(6,12) )), type = 'l', col = 'red')
points(x2,y2, col = 'red')

AFo = predict(M,data.frame(x=seq(6,12,0.01) ))
LDo = predict(M2,data.frame(x2=seq(6,12,0.01) ))

plot(seq(6,12,0.01), AFo, type = 'l')
lines(seq(6,12,0.01), LDo, col = 'red')


inverse = function (f, lower = -1000, upper = 10000) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]$root
}



  AFfun = inverse(function(S){predict(M, data.frame(x=S))})
  LDfun = inverse(function(S){predict(M2, data.frame(x2=S))})
  
  x= c()
  for (err in seq(550,2500, 100)){
    uniroot(function(Sdiff){AFfun(err)-LDfun(err)-Sdiff}, interval = c(-5,5))
    x = c(x, uniroot(function(Sdiff){AFfun(err)-LDfun(err)-Sdiff}, interval = c(-5,5))$root)
  }
  
  plot(seq(550,2500, 100), exp(x))

}

# examine accepted points
abc()
####################################################################################################################################################
####################################################################################################################################################SS_nine = make_SS(sim_out, real_out_list, S_select = S9, gens_select = gens_select, subpops_select_sim = c(31,32,33,40,41,42,49,50,51), subpops_select_real = c(31,32,33,40,41,42,49,50,51), assume_equal = T, LD_info = T, edge_length_sim = 9, edge_length_real = 9 )
SS_test = make_SS(sim_out, real_out_list, S_select = round(1000/2/3), gens_select = c(0,10,20), subpops_select_sim = 41:42, subpops_select_real = 61:62, assume_equal = T, LD_info = T,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_test = SS_test$sim; SS_real_test = SS_test$real

n = 32; print(round(Nm_values_real[n,], 4))
par(mfrow = c(2,2), ask = F)

vals = c()
for (i in 1:100){
o <- do_abc_PLS(SS_sim_nine_AF, SS_real_nine_AF[[19]][n,], PLS_use = 5, Nm_values_sim, method = 'neuralnet', tol = tol)#, numnet = 5, sizenet = 50)
if (var(o$adj.values[,1]) > 1e12) break
if (var(o$adj.values[,1] > max(vals))) o_big = o
vals = c(vals, var(o$adj.values[,1]))
}

var(o$adj.values[,1])

max(vals)
min(vals)

{
  par(mfrow = c(1,1))
  plot(Nm_values_sim[,1], Nm_values_sim[,2], xlab = 'N', ylab = 'm', pch = 16, cex = .5, col = 'grey', xlim = c(200,20000), ylim = c(0,1), xaxs = 'i', yaxs = 'i')#, pch = 16)
  
  points(o$unadj.values[,1], o$unadj.values[,2], pch = 16, col = 'red', cex = .7)
  points(o$adj.values[,1], o$adj.values[,2], pch = 16, col = 'blue', cex = .7)
  points(Nm_values_real[n,1],Nm_values_real[n,2], pch = 8, cex = 2, col = 'black', lw = 3)
  points(summary(o)[4,1],summary(o)[4,2], pch = 8, cex = 2, col = 'dark blue', lw = 3)
  
  
  for (x in seq(500,10000,1000)){
    lines(seq(2000,20000), x/seq(2000,20000), col = 'black', lw = 2, lty = 3)
  }
  
  
  
  legend('topright', col = c('red', 'green'), pch = c(16, 17), legend = c('accepted points', 'target point'))
  legend('topright', col = c('grey','red', 'blue', 'black', 'blue', 'black'), pch = c(16, 16, 16, 8,8, NA), lty = c(NA, NA, NA, NA,NA, 3), legend = c('rejected points','accepted points (unadjusted)', 'accepted points (adjusted)', 'target point', 'posterior mean', 'Nm isocline'))
  #legend('topright', col = c('grey','red', 'black', 'black'), pch = c(16, 16, 8, NA), lty = c(NA, NA, NA, 3), legend = c('rejected points','accepted points (unadjusted)', 'target point', 'Nm isocline'))
}



o_mod = o

summary(o_mod)
(o_mod$residuals)




o_mod = remove_outlier(o_mod)


summary(o_mod$adj.values)

summary(o_mod)





mod = PLS_examine(SS_sim_test, SS_real_test[[1]][1,], Nm_values_sim,PLS_use = 6)
summary(mod)


SS_test = make_SS(sim_out, real_out_list, S_select = Inf, gens_select = c(0,10,20), subpops_select_sim = 41, subpops_select_real = 61, assume_equal = F, LD_info = F,  focal_pop_sim = 41, focal_pop_real = 61, edge_length_sim = 9, edge_length_real = 11 )
SS_sim_test = SS_test$sim; SS_real_test = SS_test$real

n = 39; print(round(Nm_values_real[n,], 4))
par(mfrow = c(2,2), ask = F)

o <- do_abc_PLS(SS_sim_one_AF, SS_real_one_AF[[6]][n,], PLS_use = 2, Nm_values_sim, method = 'neuralnet', tol = .1, sizenet = 10)#, numnet = 5, sizenet = 50)

{
  par(mfrow = c(1,1))
  plot(Nm_values_sim[,1], Nm_values_sim[,2], xlab = 'N', ylab = 'm', pch = 16, cex = .5, col = 'grey', xlim = c(200,20000), ylim = c(0,1), xaxs = 'i', yaxs = 'i')#, pch = 16)
  
  points(o$unadj.values[,1], o$unadj.values[,2], pch = 16, col = 'red', cex = .7)
  points(o$adj.values[,1], o$adj.values[,2], pch = 16, col = 'blue', cex = .7)
  points(Nm_values_real[n,1],Nm_values_real[n,2], pch = 8, cex = 2, col = 'black', lw = 3)
  points(summary(o)[4,1],summary(o)[4,2], pch = 8, cex = 2, col = 'dark blue', lw = 3)
  
  
  for (x in seq(500,10000,1000)){
    lines(seq(2000,20000), x/seq(2000,20000), col = 'black', lw = 2, lty = 3)
  }
  
  
  
  legend('topright', col = c('red', 'green'), pch = c(16, 17), legend = c('accepted points', 'target point'))
  legend('topright', col = c('grey','red', 'blue', 'black', 'blue', 'black'), pch = c(16, 16, 16, 8,8, NA), lty = c(NA, NA, NA, NA,NA, 3), legend = c('rejected points','accepted points (unadjusted)', 'accepted points (adjusted)', 'target point', 'posterior mean', 'Nm isocline'))
  #legend('topright', col = c('grey','red', 'black', 'black'), pch = c(16, 16, 8, NA), lty = c(NA, NA, NA, 3), legend = c('rejected points','accepted points (unadjusted)', 'target point', 'Nm isocline'))
}




filt = Nm_values_sim[,1]>1000 & Nm_values_sim[,1]>0.01


{
  upperLidx = (abs(Nm_values_sim[,1] -8600) < 500 & abs(Nm_values_sim[,2] - .78) < 0.03)
  lowerRidx = (abs(Nm_values_sim[,1] - 17500) < 500 & abs(Nm_values_sim[,2] - .4) < 0.03)
  upperL = Nm_values_sim[upperLidx,]
  lowerR = Nm_values_sim[lowerRidx,]
  
  
  #plot(Nm_values_sim[,1], Nm_values_sim[,2], xlab = 'N', ylab = 'm', pch = 16)
  # points(upperL[,1], upperL[,2], col = 'green', pch = 16)
  # points(lowerR[,1], lowerR[,2], col = 'green', pch = 16)
  
  
  
  upperLests =  colMeans(SS_sim[upperLidx,])
  lowerRests = colMeans(SS_sim[lowerRidx,])
  
  somd = cbind(round(upperLests,-2),round(lowerRests,-2))
  colnames(somd) = c('upper L', 'lower R')
  somd
  
  pc <- prcomp(log2(SS_sim_test[filt,]), scale = F, center= F)
  
  plot(pc$x[,1], pc$x[,2], xlab = 'PC1', ylab = 'PC2', cex = 1)
  points(pc$x[accept_idx,1], pc$x[accept_idx,2], col = 'red', pch = 16)
  points(pc$x[lowerRidx,1], pc$x[lowerRidx,2], col = 'blue', pch = 16)
  
  legend('topleft' , col = c('red', 'blue'), legend = c('upper L', 'lower R'), pch = 16)
  
  pal = colorRampPalette(c("blue", "red"))(length(pc$x[,1]))
  plot(pc$x[,1], pc$x[,2], xlab = 'PC1', ylab = 'PC2', col = pal[order(Nm_values_sim[,2])], cex = 2*range01(Nm_values_sim[,1]), pch = 1)
  
}




#################################################################
#######################
#######################                  END
#######################
#################################################################








# posterior means
par(mfrow = c(2,1))
plot(a[,1,2], a[,2,1], ylim = c(7000,13000))
plot(a[,1,2], a[,2,2], ylim = c(0,1))

par(mfrow = c(2,1))
plot(a[,1,2], a[,2,1]/a[,1,1], ylim = c(.5,1.5))
plot(a[,1,2], a[,2,2]/a[,1,2], ylim = c(.5,1.5))

# posterior variance
par(mfrow = c(2,1))
plot(a[,1,2], a[,3,1], ylim = c(0,1e5))
plot(a[,1,2], a[,3,2], ylim = c(0,1e-4))

# posterior CV
par(mfrow = c(2,1))
plot(a[,1,2], sqrt(a[,3,1])/a[,1,1], ylim = c(0,.1))
plot(a[,1,2], sqrt(a[,3,2])/a[,1,2], ylim = c(0,.1))

# small m

a = test_specific(do_abc_PLS, SS_sim, SS_real[95:120,], n_comps = 5, Nm_values_sim, Nm_values_real[95:120,] )
# posterior means
par(mfrow = c(2,1))
plot(a[,1,2], a[,2,1], ylim = c(9000,11000))
plot(a[,1,2], a[,2,2], ylim = c(0,0.05))

par(mfrow = c(2,1))
plot(a[,1,2], a[,2,1]/a[,1,1], ylim = c(.5,1.5))
plot(a[,1,2], a[,2,2]/a[,1,2], ylim = c(0,3))

# posterior variance
par(mfrow = c(2,1))
plot(a[,1,2], a[,3,1], ylim = c(0,1.5e5))
plot(a[,1,2], a[,3,2], ylim = c(0,1e-4))

# posterior CV
par(mfrow = c(2,1))
plot(a[,1,2], sqrt(a[,3,1])/a[,1,1], ylim = c(0,.1))
plot(a[,1,2], sqrt(a[,3,2])/a[,1,2], ylim = c(0,2.1))



# cross validation
# cross validation nnet with tol ~ 0.1 is best
{
  tols = c(0.01, 0.02,0.05,.1, 0.2)
  
  cv.nnet = cv4abc(Nm_values_sim, sumstat = SS_sim, nval = 100, statistic = 'mean', tols = c(0.01, 0.02,0.05,.1, 0.2), method = 'neuralnet')
  plot(cv.nnet)
  summary(cv.nnet)
  
  cv.lin = cv4abc(Nm_values_sim, sumstat = SS_sim, nval = 100, statistic = 'mean', tols = c(0.01, 0.02,0.05,.1, 0.2), method = 'loclinear')
  plot(cv.lin, ask = F)
  summary(cv.lin)
  
  cv.rej = cv4abc(Nm_values_sim, sumstat = SS_sim, nval = 100, statistic = 'mean', tols = c(0.01, 0.02,0.05,.1, 0.2), method = 'rejection')
  par(mfrow = c(2,1))
  plot(cv.rej)
  summary(cv.rej)
  
  
  par(mfrow = c(2,1))
  cols = c('black', 'red', 'blue'); i = 1
  plot(NA, xlim = c(0,0.2), ylim = c(0,1), xlab = 'tolerance', ylab = 'LOO-CV prediction error (N)')
  for (obj in list(cv.nnet, cv.lin,cv.rej)){
    points(tols, summary(obj)[,1], col = cols[i]); i = i+1
  }
  legend('topright', legend = c('nerual net', 'loclinear', 'rejection'), col = cols, lty=1)
  
  cols = c('black', 'red', 'blue'); i = 1
  plot(NA, xlim = c(0,0.2), ylim = c(0,1), xlab = 'tolerance', ylab = 'LOO-CV prediction error (m)')
  for (obj in list(cv.nnet, cv.lin,cv.rej)){
    points(tols, summary(obj)[,2], col = cols[i]); i = i+1
  }
  legend('topright', legend = c('nerual net', 'loclinear', 'rejection'), col = cols, lty=1)
  
  
}


# N plot
{
  
  idx_min = 1
  idx_max = 37
  
  Ns_values = cbind(Nm_values_real[idx_min:idx_max,], apply(Nm_values_real[idx_min:idx_max,], 1, prod))
  
  mean_v_N.nn = array(NA, dim = c(0,3)); colnames(mean_v_N.nn) = c('N' , 'm', 'Nm')
  var_v_N.nn = array(NA, dim = c(0,3)); colnames(var_v_N.nn) = c('N' , 'm', 'Nm')
  
  mean_v_N.rej= array(NA, dim = c(0,3)); colnames(mean_v_N.rej) = c('N' , 'm', 'Nm')
  var_v_N.rej = array(NA, dim = c(0,3)); colnames(var_v_N.rej) = c('N' , 'm', 'Nm')
  
  for (r in idx_min:idx_max){
    o.nn <- do_abc(SS_sim, SS_real_M[r,], Nm_values_sim, method = 'neuralnet', tol = .1, pca = F)
    summary_o.nn = summary(o.nn)
    mean_v_N.nn = rbind(mean_v_N.nn, c(summary_o.nn[4,], prod(summary_o.nn[4,])) )
    var_v_N.nn = rbind(var_v_N.nn, c(colVars(o.nn$adj.values), var(apply(o.nn$adj.values, 1, prod)) ) )
    
    o.rej <- do_abc(SS_sim, SS_real_M[r,], Nm_values_sim, method = 'rejection', tol = .02, pca = F)
    summary_o.rej = summary(o.rej)
    mean_v_N.rej = rbind(mean_v_N.rej, c(summary_o.rej[4,], prod(summary_o.rej[4,])) )
    var_v_N.rej = rbind(var_v_N.rej, c(colVars(o.rej$unadj.values), var(apply(o.rej$unadj.values, 1, prod)) ) )
    
  }
  
  CV_v_N.nn = sqrt(var_v_N.nn) / mean_v_N.nn
  
  par(mfrow = c(3,1), mai = c(1,1,.2,.2))
  plot(Ns_values[,1], mean_v_N.nn[,1], main = ' m = 0.1', xlab = 'true N', ylab = expression(hat(N))); abline(a = 0, b = 1, col = 'red')
  plot(Ns_values[,1], mean_v_N.nn[,2], main = ' N = 1e4', xlab = 'true N', ylab = expression(hat(m)), ylim = c(0.03,0.07)) ; abline( h = 0.05, col = 'red')
  plot(Ns_values[,1], mean_v_N.nn[,3], main = ' N = 1e4', xlab = 'true N', ylab = expression(hat(Nm)), ylim = c(0,2000)) ; lines(Ns_values[,1], Ns_values[,3], col = 'red')
  
  plot(Ns_values[,1], var_v_N.nn[,1], main = ' m = 0.1', xlab = 'true N', ylab = expression(Var(hat(N))))
  plot(Ns_values[,1], var_v_N.nn[,2], main = ' N = 1e4', xlab = 'true N', ylab = expression(Var(hat(m))))
  plot(Ns_values[,1], var_v_N.nn[,3], main = ' N = 1e4', xlab = 'true N', ylab = expression(Var(hat(Nm))))
  
  plot(Ns_values[,1], CV_v_N.nn[,1], main = ' m = 0.1', xlab = 'true N', ylab = expression(Var(hat(N))), ylim = c(0,1))
  plot(Ns_values[,1], CV_v_N.nn[,2], main = ' N = 1e4', xlab = 'true N', ylab = expression(Var(hat(m))), ylim = c(0,1))
  plot(Ns_values[,1], CV_v_N.nn[,3], main = ' N = 1e4', xlab = 'true N', ylab = expression(Var(hat(Nm))), ylim = c(0,1))
  
  if (F){
    CV_v_N.rej = sqrt(var_v_N.rej) / mean_v_N.rej
    
    par(mfrow = c(3,1), mai = c(1,1,.2,.2))
    plot(Ns_values[,1], mean_v_N.rej[,1], main = ' m = 0.1', xlab = 'true N', ylab = expression(hat(N))); abline(a = 0, b = 1, col = 'red')
    plot(Ns_values[,1], mean_v_N.rej[,2], main = ' N = 1e4', xlab = 'true N', ylab = expression(hat(m)), ylim = c(0.05,0.15)) ; abline( h = 0.1, col = 'red')
    plot(Ns_values[,1], mean_v_N.rej[,3], main = ' N = 1e4', xlab = 'true N', ylab = expression(hat(Nm)), ylim = c(0,2000)) ; lines(Ns_values[,1], Ns_values[,3], col = 'red')
    
    plot(Ns_values[,1], var_v_N.rej[,1], main = ' m = 0.1', xlab = 'true N', ylab = expression(Var(hat(N))))
    plot(Ns_values[,1], var_v_N.rej[,2], main = ' N = 1e4', xlab = 'true N', ylab = expression(Var(hat(m))))
    plot(Ns_values[,1], var_v_N.rej[,3], main = ' N = 1e4', xlab = 'true N', ylab = expression(Var(hat(Nm))))
    
    plot(Ns_values[,1], CV_v_N.rej[,1], main = ' m = 0.1', xlab = 'true N', ylab = expression(Var(hat(N))), ylim = c(0,1))
    plot(Ns_values[,1], CV_v_N.rej[,2], main = ' N = 1e4', xlab = 'true N', ylab = expression(Var(hat(m))), ylim = c(0,1))
    plot(Ns_values[,1], CV_v_N.rej[,3], main = ' N = 1e4', xlab = 'true N', ylab = expression(Var(hat(Nm))), ylim = c(0,1))
  }
  
}

# m plot
{
  
  
  idx_min = 38
  idx_max = 56
  ms_values = cbind(Nm_values_real[idx_min:idx_max,], apply(Nm_values_real[idx_min:idx_max,], 1, prod))
  
  mean_v_m.nn = array(NA, dim = c(0,3)); colnames(mean_v_m.nn) = c('N' , 'm', 'Nm')
  var_v_m.nn = array(NA, dim = c(0,3)); colnames(var_v_m.nn) = c('N' , 'm', 'Nm')
  
  mean_v_m.rej= array(NA, dim = c(0,3)); colnames(mean_v_m.rej) = c('N' , 'm', 'Nm')
  var_v_m.rej = array(NA, dim = c(0,3)); colnames(var_v_m.rej) = c('N' , 'm', 'Nm')
  
  for (r in idx_min:idx_max){
    o.nn <- do_abc(SS_sim, SS_real_M[r,], Nm_values_sim, method = 'neuralnet', tol = .1, pca = F)
    summary_o.nn = summary(o.nn)
    mean_v_m.nn = rbind(mean_v_m.nn, c(summary_o.nn[4,], prod(summary_o.nn[4,])) )
    var_v_m.nn = rbind(var_v_m.nn, c(colVars(o.nn$adj.values), var(apply(o.nn$adj.values, 1, prod)) ) )
    
    o.rej <- do_abc(SS_sim, SS_real_M[r,], Nm_values_sim, method = 'rejection', tol = .1, pca = F)
    summary_o.rej = summary(o.rej)
    mean_v_m.rej = rbind(mean_v_m.rej, c(summary_o.rej[4,], prod(summary_o.rej[4,])) )
    var_v_m.rej = rbind(var_v_m.rej, c(colVars(o.rej$unadj.values), var(apply(o.rej$unadj.values, 1, prod)) ) )
    
  }
  
  CV_v_m.nn = sqrt(var_v_m.nn) / mean_v_m.nn
  
  par(mfrow = c(3,1), mai = c(1,1,.2,.2))
  plot(ms_values[,2], mean_v_m.nn[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(m)), ylim = c(0,1)); abline(a=0,b=1, col = 'red')
  plot(ms_values[,2], mean_v_m.nn[,1], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(N))) ; abline(h=1e4, col = 'red')
  plot(ms_values[,2], mean_v_m.nn[,3], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(Nm)), ylim = c(0,10000)) ; lines(ms_values[,2], ms_values[,3], col = 'red')
  
  par(mfrow = c(3,1), mai = c(1,1,.2,.2))
  plot(ms_values[,2], log2(mean_v_m.nn[,2]/ms_values[,2]), main = ' N = 1e4', xlab = 'true m', ylab = expression(log2(hat(m)/m))); abline(h = 0, col = 'red')
  plot(ms_values[,2], log2(mean_v_m.nn[,1]/ms_values[,1]), main = ' N = 1e4', xlab = 'true m', ylab = expression(log2(hat(N)/N))) ; abline(h = 0, col = 'red')
  plot(ms_values[,2], log2(mean_v_m.nn[,3]/ms_values[,3]), main = ' N = 1e4', xlab = 'true m', ylab = expression(log2(hat(Nm)/Nm))) ; abline(h = 0, col = 'red')
  
  
  par(mfrow = c(3,1), mai = c(1,1,.2,.2))
  plot(ms_values[,2], var_v_m.nn[,2], main = ' m = 0.1', xlab = 'true m', ylab = expression(Var(hat(m))))
  plot(ms_values[,2], var_v_m.nn[,1], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(N))))
  plot(ms_values[,2], var_v_m.nn[,3], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(Nm))), ylim = c(0,4e6))
  
  par(mfrow = c(1,1))
  plot(ms_values[,2], CV_v_m.nn[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(CV), ylim = c(0,1), pch = 16)
  points(ms_values[,2], CV_v_m.nn[,1], col = 'blue', pch = 16)
  points(ms_values[,2], CV_v_m.nn[,3], col = 'red', pch = 16)
  legend('topright',legend = c('N', 'm', 'Nm'), col = c('black', 'blue', 'red'), pch = 16)
  
  
  
  
  if (F){
    CV_v_m.rej = sqrt(var_v_m.rej) / mean_v_m.rej
    
    par(mfrow = c(3,1), mai = c(1,1,.2,.2))
    plot(ms_values[,2], mean_v_m.rej[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(m))); abline(a = 0, b = 1, col = 'red')
    plot(ms_values[,2], mean_v_m.rej[,1], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(N))) ; abline( h = 1e4, col = 'red')
    plot(ms_values[,2], mean_v_m.rej[,3], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(Nm))) ; lines(ms_values[,2], ms_values[,3], col = 'red')
    
    plot(ms_values[,2], var_v_m.rej[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(m))))
    plot(ms_values[,2], var_v_m.rej[,1], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(N))))
    plot(ms_values[,2], var_v_m.rej[,3], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(Nm))))
    
    par(mfrow = c(1,1))
    plot(ms_values[,2], CV_v_m.rej[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(CV), ylim = c(0,1), pch = 16)
    points(ms_values[,2], CV_v_m.rej[,1], ylim = c(0,1), col = 'blue', pch = 16)
    points(ms_values[,2], CV_v_m.rej[,3], ylim = c(0,1), col = 'red', pch = 16)
    legend('topright',legend = c('N', 'm', 'Nm'), col = c('black', 'blue', 'red'), pch = 16)
    
  }
}

# m plot multi
{
  
  
  idx_mins = c(38, 57, 76)
  cols = c('black', 'red', 'blue')
  Nms = list()
  means = list()
  vars = list()
  i = 1
  
  for (idx_min in idx_mins){
    idx_max = idx_min + 18
    
    
    ms_values = cbind(Nm_values_real[idx_min:idx_max,], apply(Nm_values_real[idx_min:idx_max,], 1, prod))
    
    mean_v_m.nn = array(NA, dim = c(0,3)); colnames(mean_v_m.nn) = c('N' , 'm', 'Nm')
    var_v_m.nn = array(NA, dim = c(0,3)); colnames(var_v_m.nn) = c('N' , 'm', 'Nm')
    
    mean_v_m.rej= array(NA, dim = c(0,3)); colnames(mean_v_m.rej) = c('N' , 'm', 'Nm')
    var_v_m.rej = array(NA, dim = c(0,3)); colnames(var_v_m.rej) = c('N' , 'm', 'Nm')
    
    for (r in idx_min:idx_max){
      o.nn <- do_abc(SS_sim, SS_real_M[r,], Nm_values_sim, method = 'neuralnet', tol = .1, pca = F)
      summary_o.nn = summary(o.nn)
      mean_v_m.nn = rbind(mean_v_m.nn, c(summary_o.nn[4,], prod(summary_o.nn[4,])) )
      var_v_m.nn = rbind(var_v_m.nn, c(colVars(o.nn$adj.values), var(apply(o.nn$adj.values, 1, prod)) ) )
      
      o.rej <- do_abc(SS_sim, SS_real_M[r,], Nm_values_sim, method = 'rejection', tol = .1, pca = F)
      summary_o.rej = summary(o.rej)
      mean_v_m.rej = rbind(mean_v_m.rej, c(summary_o.rej[4,], prod(summary_o.rej[4,])) )
      var_v_m.rej = rbind(var_v_m.rej, c(colVars(o.rej$unadj.values), var(apply(o.rej$unadj.values, 1, prod)) ) )
      
    }
    
    Nms[[i]] = ms_values
    means[[i]] = mean_v_m.nn
    vars[[i]] = var_v_m.nn
    i = i+1
  }
  
  
  CV_v_m.nn = sqrt(var_v_m.nn) / mean_v_m.nn
  
  
  # means
  par(mfrow = c(3,1), mai = c(1,1,.2,.2))
  
  plot(NA, ylim = c(0,1),xlim = c(0,1), main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(m))); abline(a=0,b=1, col = 'red')
  for (i in 1:3) lines(Nms[[i]][,2], means[[i]][,2], col = cols[i])
  legend('topright',legend = c(10000, 5000, 15000), col = cols, pch = 1)
  
  plot(NA, ylim = c(0,20000),xlim = c(0,1), main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(N)))#; abline(, col = 'red')
  for (i in 1:3) lines(Nms[[i]][,2], means[[i]][,1], col = cols[i])
  legend('topright',legend = c(10000, 5000, 15000), col = cols, pch = 1)
  
  plot(NA, ylim = c(0,20000),xlim = c(0,1), main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(Nm)))#; abline(, col = 'red')
  for (i in 1:3) lines(Nms[[i]][,2], means[[i]][,3], col = cols[i])
  legend('topright',legend = c(10000, 5000, 15000), col = cols, pch = 1)
  
  
  
  # vars
  
  par(mfrow = c(3,1), mai = c(1,1,.2,.2))
  
  plot(NA, ylim = c(0,0.08),xlim = c(0,1), main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(m))))#; abline(a=0,b=1, col = 'red')
  for (i in 1:3) lines(Nms[[i]][,2], vars[[i]][,2], col = cols[i])
  legend('topright',legend = c(10000, 5000, 15000), col = cols, pch = 1)
  
  plot(NA, ylim = c(0,2e7),xlim = c(0,1), main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(N))))#; abline(, col = 'red')
  for (i in 1:3) lines(Nms[[i]][,2], vars[[i]][,1], col = cols[i])
  legend('topright',legend = c(10000, 5000, 15000), col = cols, pch = 1)
  
  plot(NA, ylim = c(0,5e6),xlim = c(0,1), main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(Nm))))#; abline(, col = 'red')
  for (i in 1:3) lines(Nms[[i]][,2], vars[[i]][,3], col = cols[i])
  legend('topright',legend = c(10000, 5000, 15000), col = cols, pch = 1)
  
  
  # log plots
  par(mfrow = c(3,1), mai = c(1,1,.2,.2))
  
  plot(NA, ylim = c(-1,1),xlim = c(0,1), main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(m))); abline(h=0, col = 'black')
  for (i in 1:3) lines(Nms[[i]][,2], log2(means[[i]][,2]/Nms[[i]][,2]), col = cols[i])
  legend('topright',legend = c(10000, 5000, 15000), col = cols, pch = 1)
  
  plot(NA, ylim = c(-1,1),xlim = c(0,1), main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(N))); abline(h=0, col = 'black')
  for (i in 1:3) lines(Nms[[i]][,2], log2(means[[i]][,1]/Nms[[i]][,1]), col = cols[i])
  legend('topright',legend = c(10000, 5000, 15000), col = cols, pch = 1)
  
  plot(NA, ylim = c(-1,1),xlim = c(0,1), main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(Nm))); abline(h=0, col = 'black')
  for (i in 1:3) lines(Nms[[i]][,2], log2(means[[i]][,3]/Nms[[i]][,3]), col = cols[i])
  legend('topright',legend = c(10000, 5000, 15000), col = cols, pch = 1)
  
  
  # par(mfrow = c(1,1))
  # plot(ms_values[,2], CV_v_m.nn[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(CV), ylim = c(0,1), pch = 16)
  # points(ms_values[,2], CV_v_m.nn[,1], col = 'blue', pch = 16)
  # points(ms_values[,2], CV_v_m.nn[,3], col = 'red', pch = 16)
  # legend('topright',legend = c('N', 'm', 'Nm'), col = c('black', 'blue', 'red'), pch = 16)
  # 
  # 
  
  
  if (F){
    CV_v_m.rej = sqrt(var_v_m.rej) / mean_v_m.rej
    
    par(mfrow = c(3,1), mai = c(1,1,.2,.2))
    plot(ms_values[,2], mean_v_m.rej[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(m))); abline(a = 0, b = 1, col = 'red')
    plot(ms_values[,2], mean_v_m.rej[,1], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(N))) ; abline( h = 1e4, col = 'red')
    plot(ms_values[,2], mean_v_m.rej[,3], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(Nm))) ; lines(ms_values[,2], ms_values[,3], col = 'red')
    
    plot(ms_values[,2], var_v_m.rej[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(m))))
    plot(ms_values[,2], var_v_m.rej[,1], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(N))))
    plot(ms_values[,2], var_v_m.rej[,3], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(Nm))))
    
    par(mfrow = c(1,1))
    plot(ms_values[,2], CV_v_m.rej[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(CV), ylim = c(0,1), pch = 16)
    points(ms_values[,2], CV_v_m.rej[,1], ylim = c(0,1), col = 'blue', pch = 16)
    points(ms_values[,2], CV_v_m.rej[,3], ylim = c(0,1), col = 'red', pch = 16)
    legend('topright',legend = c('N', 'm', 'Nm'), col = c('black', 'blue', 'red'), pch = 16)
    
  }
}

# small m
{
  
  idx_min = 95
  idx_max = 120
  ms_values = cbind(Nm_values_real[idx_min:idx_max,], apply(Nm_values_real[idx_min:idx_max,], 1, prod))
  
  mean_v_m.nn = array(NA, dim = c(0,3)); colnames(mean_v_m.nn) = c('N' , 'm', 'Nm')
  var_v_m.nn = array(NA, dim = c(0,3)); colnames(var_v_m.nn) = c('N' , 'm', 'Nm')
  
  mean_v_m.rej= array(NA, dim = c(0,3)); colnames(mean_v_m.rej) = c('N' , 'm', 'Nm')
  var_v_m.rej = array(NA, dim = c(0,3)); colnames(var_v_m.rej) = c('N' , 'm', 'Nm')
  
  for (r in idx_min:idx_max){
    o.nn <- do_abc(SS_sim, SS_real_M[r,], Nm_values_sim, method = 'neuralnet', tol = .1, pca = F)
    summary_o.nn = summary(o.nn)
    mean_v_m.nn = rbind(mean_v_m.nn, c(summary_o.nn[4,], prod(summary_o.nn[4,])) )
    var_v_m.nn = rbind(var_v_m.nn, c(colVars(o.nn$adj.values), var(apply(o.nn$adj.values, 1, prod)) ) )
    
    o.rej <- do_abc(SS_sim, SS_real_M[r,], Nm_values_sim, method = 'rejection', tol = .1, pca = F)
    summary_o.rej = summary(o.rej)
    mean_v_m.rej = rbind(mean_v_m.rej, c(summary_o.rej[4,], prod(summary_o.rej[4,])) )
    var_v_m.rej = rbind(var_v_m.rej, c(colVars(o.rej$unadj.values), var(apply(o.rej$unadj.values, 1, prod)) ) )
    
  }
  
  CV_v_m.nn = sqrt(var_v_m.nn) / mean_v_m.nn
  
  par(mfrow = c(3,1), mai = c(1,1,.2,.2))
  plot(ms_values[,2], mean_v_m.nn[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(m))); abline(a=0,b=1, col = 'red')
  plot(ms_values[,2], mean_v_m.nn[,1], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(N))) ; abline(h=1e4, col = 'red')
  plot(ms_values[,2], mean_v_m.nn[,3], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(Nm))) ; lines(ms_values[,2], ms_values[,3], col = 'red')
  
  par(mfrow = c(3,1), mai = c(1,1,.2,.2))
  plot(ms_values[,2], log2(mean_v_m.nn[,2]/ms_values[,2]), main = ' N = 1e4', xlab = 'true m', ylab = expression(log2(hat(m)/m))); abline(h = 0, col = 'red')
  plot(ms_values[,2], log2(mean_v_m.nn[,1]/ms_values[,1]), main = ' N = 1e4', xlab = 'true m', ylab = expression(log2(hat(N)/N))) ; abline(h = 0, col = 'red')
  plot(ms_values[,2], log2(mean_v_m.nn[,3]/ms_values[,3]), main = ' N = 1e4', xlab = 'true m', ylab = expression(log2(hat(Nm)/Nm))) ; abline(h = 0, col = 'red')
  
  
  par(mfrow = c(3,1), mai = c(1,1,.2,.2))
  plot(ms_values[,2], var_v_m.nn[,2], main = ' m = 0.1', xlab = 'true m', ylab = expression(Var(hat(m))))
  plot(ms_values[,2], var_v_m.nn[,1], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(N))))
  plot(ms_values[,2], var_v_m.nn[,3], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(Nm))), ylim = c(0,4e6))
  
  par(mfrow = c(1,1))
  plot(ms_values[,2], CV_v_m.nn[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(CV), ylim = c(0,1), pch = 16)
  points(ms_values[,2], CV_v_m.nn[,1], col = 'blue', pch = 16)
  points(ms_values[,2], CV_v_m.nn[,3], col = 'red', pch = 16)
  legend('topright',legend = c('N', 'm', 'Nm'), col = c('black', 'blue', 'red'), pch = 16)
  
  
  
  
  if (F){
    CV_v_m.rej = sqrt(var_v_m.rej) / mean_v_m.rej
    
    par(mfrow = c(3,1), mai = c(1,1,.2,.2))
    plot(ms_values[,2], mean_v_m.rej[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(m))); abline(a = 0, b = 1, col = 'red')
    plot(ms_values[,2], mean_v_m.rej[,1], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(N))) ; abline( h = 1e4, col = 'red')
    plot(ms_values[,2], mean_v_m.rej[,3], main = ' N = 1e4', xlab = 'true m', ylab = expression(hat(Nm))) ; lines(ms_values[,2], ms_values[,3], col = 'red')
    
    plot(ms_values[,2], var_v_m.rej[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(m))))
    plot(ms_values[,2], var_v_m.rej[,1], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(N))))
    plot(ms_values[,2], var_v_m.rej[,3], main = ' N = 1e4', xlab = 'true m', ylab = expression(Var(hat(Nm))))
    
    par(mfrow = c(1,1))
    plot(ms_values[,2], CV_v_m.rej[,2], main = ' N = 1e4', xlab = 'true m', ylab = expression(CV), ylim = c(0,1), pch = 16)
    points(ms_values[,2], CV_v_m.rej[,1], ylim = c(0,1), col = 'blue', pch = 16)
    points(ms_values[,2], CV_v_m.rej[,3], ylim = c(0,1), col = 'red', pch = 16)
    legend('topright',legend = c('N', 'm', 'Nm'), col = c('black', 'blue', 'red'), pch = 16)
    
  }
}

# cross validation
n=100
print(Nm_values_real[n,])
o1=do_abc_PLS(X, X_r[n,], PLS_use = 5, Nm_values_sim, method = 'neuralnet', tol = .1)
o2=do_abc_trunc(X, X_r[n,], PLS_use = 5, Nm_values_sim, method = 'neuralnet', tol = .1)


chosen = sample(1:dim(Nm_values_sim)[1], 200)
cv1 = cross_validate(do_abc_PLS, X, chosen, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .2)
cv2 = cross_validate(do_abc_trunc, X, chosen, n_comps =  5, Nm_values_sim, method = 'neuralnet', tol = .2)

loocv_cv1 = loocv_err(cv1); colMeans(loocv_cv1)
loocv_cv2 = loocv_err(cv2); colMeans(loocv_cv2)

par(mfrow=c(1,1))
plot(cv1[,1,1],cv1[,2,1]/cv1[,1,1])
plot(cv1[,1,2],cv1[,2,2]/cv1[,1,2])


plot(cv2[,1,1],cv2[,2,1]/cv2[,1,1])
plot(cv2[,1,2],cv2[,2,2]/cv2[,1,2])


plot(cv1[,2,2]/cv1[,1,2],cv1[,1,2])


mae_cv1 = MAE_err(cv1); colMeans(mae_cv1)
mae_cv2 = MAE_err(cv2); colMeans(mae_cv2)

## calculate PLS stats for real data

install.packages('abctools')
filt = Nm_values_sim[,1] > 9000 & Nm_values_sim[,1] < 11000
filt2= Nm_values_sim[,2] > .01 & Nm_values_sim[,2] <.1
f= filt&filt2
Nm_values_sim[f,]
a = abctools::saABC(theta = Nm_values_sim[f,1:2], X[f,])

SS_sim_linreg = X %*% t(a$B)
SS_real_linreg = X_r %*% t(a$B)


n = 335; print(round(Nm_values_sim[n,], 3))
par(mfrow = c(2,2), ask = F)

o <- do_abc(SS_sim_equal, SS_sim_equal[n,], Nm_values_sim, N = Nm_values_sim[n,1], m = Nm_values_sim[n,2], draw = F, method = 'neuralnet', tol = .05, pca = F, transformation = rep('none',3))
plot(o, Nm_values_sim, ask = F)

{
  par(mfrow = c(1,1))
  plot(Nm_values_sim[,1], Nm_values_sim[,2], xlab = 'N', ylab = 'm', pch = 16, cex = .5, col = 'grey', xlim = c(2000,20000), ylim = c(0,1), xaxs='i', yaxs = 'i')#, pch = 16)
  
  points(o$unadj.values[,1], o$unadj.values[,2], pch = 16, col = 'red', cex = .7)
  points(o$adj.values[,1], o$unadj.values[,2], pch = 16, col = 'blue', cex = .7)
  points(Nm_values_sim[n,1],Nm_values_sim[n,2], pch = 8, cex = 2, col = 'black', lw = 3)
  points(summary(o)[4,1],summary(o)[4,2], pch = 8, cex = 2, col = 'blue', lw = 3)
  
  
  for (x in seq(500,10000,1000)){
    lines(seq(2000,20000), x/seq(2000,20000), col = 'black', lw = 2, lty = 3)
  }
  
  
  
  legend('topright', col = c('red', 'green'), pch = c(16, 17), legend = c('accepted points', 'target point'))
  legend('topright', col = c('grey','red', 'blue', 'black', 'blue', 'black'), pch = c(16, 16, 16, 8,8, NA), lty = c(NA, NA, NA, NA,NA, 3), legend = c('rejected points','accepted points (unadjusted)', 'accepted points (adjusted)', 'target point', 'posterior mean', 'Nm isocline'))
  #legend('topright', col = c('grey','red', 'black', 'black'), pch = c(16, 16, 8, NA), lty = c(NA, NA, NA, 3), legend = c('rejected points','accepted points (unadjusted)', 'target point', 'Nm isocline'))
  
  upperLidx = (abs(Nm_values_sim[,1] -8600) < 500 & abs(Nm_values_sim[,2] - .78) < 0.03)
  lowerRidx = (abs(Nm_values_sim[,1] - 17500) < 500 & abs(Nm_values_sim[,2] - .4) < 0.03)
  upperL = Nm_values_sim[upperLidx,]
  lowerR = Nm_values_sim[lowerRidx,]
  
  
  #plot(Nm_values_sim[,1], Nm_values_sim[,2], xlab = 'N', ylab = 'm', pch = 16)
  # points(upperL[,1], upperL[,2], col = 'green', pch = 16)
  # points(lowerR[,1], lowerR[,2], col = 'green', pch = 16)
  
  
  
  upperLests =  colMeans(SS_sim_linreg[upperLidx,])
  lowerRests = colMeans(SS_sim_linreg[lowerRidx,])
  
  somd = cbind(round(upperLests,-2),round(lowerRests,-2))
  colnames(somd) = c('upper L', 'lower R')
  somd
  
  pc <- prcomp(SS_sim, scale = T, center= T)
  
  plot(pc$x[,1], pc$x[,2], xlab = 'PC1', ylab = 'PC2', cex = 1)
  points(pc$x[upperLidx,1], pc$x[upperLidx,2], col = 'red', pch = 16)
  points(pc$x[lowerRidx,1], pc$x[lowerRidx,2], col = 'blue', pch = 16)
  
  legend('topleft' , col = c('red', 'blue'), legend = c('upper L', 'lower R'), pch = 16)
  
  pal = colorRampPalette(c("blue", "red"))(length(pc$x[,1]))
  plot(pc$x[,1], pc$x[,2], xlab = 'PC1', ylab = 'PC2', col = pal[order(Nm_values_sim[,2])], cex = 2*range01(Nm_values_sim[,1]), pch = 1)
  
}









############ TODO ##############

# Implement truncation abc method
# Implement cross validation testing
### loocv prediction error (vs param value, mean)
### absolute error / param value (vs param value, mean)
### param value vs predicted

### check whether my prediction error is the same as the one implemented
### 

# Check between:
### 0) PLS (different n_use) vs PLS+truncation vs PCA - a couple of different sampling stategies
### 1) sampling stategies
### 2) LD vs AF
### 3)



































# how do fc, r2, fccor, riji, PLS vary with migration rate

# fc
plot(NA, xlim = c(0,1), ylim = c(0,6e-4))
for (row in 1:5){
  lines(Nm_values_real[38:56,2], Fc[row,38:56])
}


# r2
plot(NA, xlim = c(0,1), ylim = c(0,4e-4))
for (row in 1:4){
  lines(Nm_values_real[38:56,2], r2[row,38:56])
}


# rirj
par(mfrow = c(2,2), mai = c(0.3,0.3,0.2,0.1))

for (dist in 1:4){
  plot(NA, xlim = c(0,1), ylim = range( rirj[,dist,38:56]))
  text(x = .8,y=max(rirj[,dist,38:56]),paste('Dist =', dist))
  for (row in 1:4){
    lines(Nm_values_real[38:56,2], rirj[row,dist,38:56])
  }
}

#FcCor

par(mfrow = c(2,2), mai = c(0.3,0.3,0.2,0.1))

for (dist in 1:4){
  plot(NA, xlim = c(0,1), ylim = range( FcCor[,dist,38:56]))
  text(x = .8,y=max(FcCor[,dist,38:56]),paste('Dist =', dist))
  for (row in 1:4){
    lines(Nm_values_real[38:56,2], FcCor[row,dist,38:56])
  }
}


## do abc



## examine abc results

## compare SS, sampling strategies, value of LD information














a= FcCor
dim(a) = c(120,20)

Fc2d = make_2d(FcCor)
r22d = make_2d(rirj)

X = cbind(Fc2d[,], r22d[,] )

install.packages('pls')
library(pls)

mod = pls::plsr(Nm_values_real[,1:2] ~ X, ncomp=5)

plot(mod, ncomp = 3, asp = 1, line = TRUE)

summary(mod)

par(mfrow = c(1,1))
plot(Nm_values_real[,2], mod$scores[,1])

data(yarn)
data(oliveoil)
data(gasoline)

dim(gasoline)

# do ABC
## need to decide what correlation parameters to use
## possible use PCA, PLS etc to reduce dimensionality


# calculate estimates for many methods, many deme sizes, many m values
{
  loci = 1000; max_loci = loci * 1.1; min_loci  = loci; n_loci = loci
  N = 100000
  Mfn = construct_M_4_step
  ms <- rev(c(.01, .1, .2, .5, .8, .95))
  n_demes_v <- c(100)
  el = 2
  sampling_intervals = c(2,5,10)
  all_intervals = get_all_intervals(c(0, sampling_intervals))
  recom_rates = c(.5, .35, .2, .05)
  AF_ests = array(NA, dim = c(length(ms), length(n_demes_v), length(all_intervals)))
  LD_ests = array(NA, dim = c(length(ms), length(n_demes_v), length(recom_rates) ))
  
  AF_raw = array(NA, dim = c(length(ms), length(n_demes_v), length(all_intervals)))
  LD_raw = array(NA, dim = c(length(ms), length(n_demes_v), length(recom_rates) ))
  
  FcCor = array(NA, dim = c(length(ms), length(n_demes_v), length(all_intervals), 2*el-2))
  riji  = array(NA, dim = c(length(ms), length(n_demes_v), length(recom_rates), 2*el-2))
  
  
  i = 1; j = 1; kAF = 1; kLD = 1
  for (grid_size in n_demes_v){
    sample_subpops = get_square_around_centre(grid_size, el)

    model_fn <- get_model_fn_multi_pops(sampling_intervals, recom_rates, Mfn, S, grid_size, sample_subpops, max_loci, min_loci, n_loci)
    i = 1
    for (m in ms){
      o <- model_fn(N, m)
      
      Fc_v = apply(o$Fc[,1,,], 2, mean)
      r2_v = apply(o$r2[,1,,], 1, mean)
      
      AF_ests[i,j,] = sampling_intervals/(2*Fc_v)/N
      LD_ests[i,j,] = (1-r2_v)/(2*recom_rates*(2-recom_rates)*r2_v)/N
      
      AF_raw[i,j,] = Fc_v
      LD_raw[i,j,] = r2_v
      
      
      
      FcCor[i,j,,] = apply(o$FcCor[,1,,], c(2,3), mean)
      riji[i,j,,] = apply(o$riji[,1,,], c(1,3), mean)
      i = i+1
    }
    j = j+1
    gc()
  }
}


{
  size_idx = 1
  size_selected = n_demes_v[size_idx]
  
  AF_ests_size = AF_raw[,size_idx,]
  LD_ests_size = LD_raw[,size_idx,]
  
  #,,interval/recomrate, subpop separation
  plot(ms, FcCor[,1,2,2])
  plot(ms, riji[,1,4,4])
  
  
  
  ############################
  par(mfrow = c(5,2), mai = c(.3,.3,0,0))
  for (i in 1:5){
    for (j in 1:2){
      plot(ms, FcCor[,1,i,j])
    }
  }
  
  
  ##############################
  par(mfrow = c(4,2), mai = c(.3,.3,0,0))
  for (i in 1:4){
    for (j in 1:2){
      plot(ms, riji[,1,i,j])
    }
  }
  
  a = convert_to_2d(FcCor[,1,,], riji[,1,,])
  pc <- prcomp(a, scale = T, center= T)
  par(mfrow = c(3,3))
  for (n in 1:6){
  plot(ms, pc$x[,n])
  }
  ################################
  
  convert_to_2d=function(fccor, riji){
    d1 = dim(fccor)
    d2 = dim(riji)
    dim(fccor) = c(d1[1], prod(d1[2:3]))
    dim(riji) = c(d2[1], prod(d2[2:3]))
    return(cbind(fccor, riji))
  }
  
  
  plot(ms, AF_raw[,1,2])
  plot(ms, LD_raw[,1,2])
  
  par(mfrow = c(2,1), mai = c(1,1,.3,.3))
  plot(NA, type = 'l', log = 'x', ylim = range(AF_ests_size), xlim = range(ms), xlab = 'm', ylab = expression(bar(Fc)), main = paste('n_demes = ', size_selected)); abline(h = 1); axis(side=2, at=seq(1,20,1), labels = T)
  
  
  colfunc_AF <- colorRampPalette(c("red", "green"))
  cols = colfunc_AF(length(sampling_intervals))
  for (col in 1:length(sampling_intervals)){
    lines(ms, AF_ests_size[,col], col = cols[col], cex = 1.5)
  }
  
  
  leg = c(paste('AF', sampling_intervals))
  col = c(cols)
  lty = c(rep(1, length(sampling_intervals)))
  
  legend('topright', legend = leg, col = col, lty = lty, title = 'method')
  
  
  plot(NA, type = 'l', log = 'x', ylim = range(LD_ests_size), xlim = range(ms), xlab = 'm', ylab = expression(bar(r^2)), main = paste('n_demes = ', size_selected)); abline(h = 1); axis(side=2, at=seq(1,20,1), labels = T)
  
  for (col in 1:length(recom_rates)){
    lines(ms, LD_ests_size[,col], lty=col, cex = 1.5)
    
  }
  
  leg = c(paste('LD', recom_rates))
  col = c(rep('black', length(recom_rates)))
  lty = c(1:length(recom_rates))
  
  legend('topright', legend = leg, col = col, lty = lty, title = 'method')
  
}




filt = Nm_values_sim[,2] < 0.05
mod = PLS_examine(SS_sim = SS_sim_one_AF[filt,], SS_real = SS_real_one_AF[[1]][1,], Nm = Nm_values_sim[filt,],PLS_use = dim(SS_sim_one_AF)[2])
summary(mod)


filt = Nm_values_sim[,2] > 0.05

mod = PLS_examine(SS_sim = SS_sim_three_AF[filt,], SS_real = SS_real_three_AF[[1]][1,], Nm = Nm_values_sim[filt,],PLS_use = dim(SS_sim_three_AF)[2])
summary(mod)





  