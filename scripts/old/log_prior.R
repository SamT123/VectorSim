setwd('Desktop/VectorSim/scripts/')

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


sim_out <- readRDS(file = '../data/sim_out_log_4.rds')
dimnames(sim_out$Fc_sim)[[1]]

# redefine priors
N_prior_draw_log = function(n=1){
  logged = runif(n,3,6)
  return(10^logged)
}

Nm_draw_log <- function(size){
  Nm_unordered <- array(c(N_prior_draw_log(size), m_prior_draw(size)), dim = c(size, 2))
  
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

Nm_draw_large_log <- function(size){
  Nm_u <- array(c(N_prior_draw_log(size), m_prior_draw(size)), dim = c(size, 2))
  
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

log_N = Nm_draw_large_log(1000)

hist((log_N[,1]))

# simulation
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
    pop = initialise_population(n_subpops = grid_size_sim, loci = loci)
    
    m = 0.4
    N = 7000
    Mfn = construct_M_4_step
    
    model_fn <- get_model_fn_multi_pops_focal(sample_gens, recom_rates, Mfn, S_vec, grid_size_sim, sample_subpops_sim, max_loci, min_loci, n_loci, focal_pop)
    Rprof()
    o <- model_fn(N, m, f = get_SS_multi_pops2)
    Rprof(NULL)
    print(summaryRprof())
    
  }
  
  ########### make data ############  
  
  
  # generate simulated data
  {
    n_par = 8
    n_samp = 350
    
    # constructing some things
    {
      model_fn_list = list()
      Nm_list = list()
      for (i in 1:n_par){
        model_fn_list[[i]] = get_model_fn_multi_pops(sample_gens, recom_rates, Mfn, S_vec, grid_size_sim, sample_subpops_sim, max_loci, min_loci, n_loci)
        Nm_list[[i]] =  Nm_draw_large_log(n_samp)
      }
    }
    
    
    s = proc.time()
    sim_out <- generate_ABC_sim_parallel_multi(loci, sample_gens, recom_rates, Nm_list, model_fn_list, n_par, S_vec, Mfn, grid_size_sim, sample_subpops_sim)
    e = proc.time(); print(e-s)
    
    saveRDS(sim_out, '../data/sim_out_log_5.rds')
    sim_out <- readRDS('../data/sim_out_log_5.rds')
    
    Nm_values_sim = sim_out$Nm_values
    Nm_values_sim = cbind(Nm_values_sim, apply(Nm_values_sim,1,prod))
  }
  
  

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
    reps  = 3
    n_par = 8
    focal_pop = 61
    
    
    # vary N, m = 0.05
    
    N1 = 10^seq(3.1,5.9,0.2)
    m1 = rep(0.05, length(N1))
    m2 = rep(0.30, length(N1))
    
    
    Ns = c(N1, N1)
    ms = c(m1, m2)
    
    
    # constructing some things
    {
      Ns_l = split(Ns,ceiling(seq_along(Ns)/4))
      ms_l = split(ms,ceiling(seq_along(ms)/4))
      
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
    
    
    saveRDS(real_out_list, '../data/real_out_log_4.rds')
    real_out_list <- readRDS('../data/real_out_log_4.rds')
    
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
    sim_files = c('../data/sim_out_log_2.rds','../data/sim_out_log_3.rds','../data/sim_out_log_4.rds','../data/sim_out_log_5.rds')

    
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
    real_files = c('../data/real_out_log_1.rds', '../data/real_out_log_2.rds', '../data/real_out_log_3.rds', '../data/real_out_log_4.rds')

    
    real_out_list <- list()
    for (f in real_files){
      real_out_list = c(real_out_list, readRDS(f))
    }
    
    Nm_values_real = real_out_list[[1]]$Nm_values
  }
  
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





# making statistics
{
  
  SS_one_AF = make_SS(sim_out, real_out_list, S_select = round(get_diploid_equivalent(1e100, round(1000/3/2)*30, 30)), gens_select = c(0,10,20), subpops_select_sim = c(41,42), subpops_select_real = c(61,62), assume_equal = T, LD_info = F,  focal_pop_sim = 61, focal_pop_real = 41, edge_length_sim = 9, edge_length_real = 11 )
  SS_sim_one_AF = SS_one_AF$sim[,]; SS_real_one_AF = SS_one_AF$real
  
}

out = test_specific(do_abc_PLS, SS_sim_one_AF, SS_real_one_AF, n_comps = 4, Nm_values_sim, Nm_values_real, idx = c(1,15), tol = tol, method = method, sizenet = sizenet, numnet = numnet)
out = test_specific(do_abc_PLS, SS_sim_one_AF, SS_real_one_AF, n_comps = 4, Nm_values_sim, Nm_values_real, idx = c(15,30), tol = tol, method = method, sizenet = sizenet, numnet = numnet)



varyN_log_0.05 = lapply(list(out), apply, c(1,2,3), mean)[[1]]
varyN_log_0.05 = lapply(list(out), apply, c(1,2,3), mean)[[1]]


plot(varyN_log_0.05[,1,1],varyN_log_0.05[,2,1]/varyN_log_0.05[,1,1], log = 'xy')
plot(log10(varyN_log_0.05[,1,1]),sqrt(varyN_log_0.05[,3,1])/varyN_log_0.05[,1,1], log = 'xy')

n = 9; Nm_values_real[n,]
o1 <- do_abc_PLS(SS_sim_one_AF, SS_real_one_AF[[5]][n,], PLS_use = 3, Nm_values_sim, method = 'neuralnet', tol = .07, numnet = 1, sizenet = 4)

plot(log10(Nm_values_sim[,1]), Nm_values_sim[,2], xlab = as.expression(bquote('N'['e'])), ylab = 'm', pch = 16, cex = .35, col = 'grey', xlim = c(2,7), ylim = c(0,0.8), xaxs = 'i', yaxs = 'i', pty = 'm')#, pch = 16)
points(log10(o1$unadj.values[,1]), o1$unadj.values[,2], pch = 16)
points(log10(o1$adj.values[,1]), o1$adj.values[,2], col = 'red', pch = 16)
points(log10(Nm_values_real[n,1]),Nm_values_real[n,2], pch = 8, col = 'blue')

points(Nm_values_sim[accept_idx2_small,1],Nm_values_sim[accept_idx2_small,2], pch = 16, col = col_arr, cex = .7)
points(Nm_values_sim[accept_idx2_large,1],Nm_values_sim[accept_idx2_large,2], pch = 16, col = col_arr2, cex = .7)















