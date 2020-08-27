# setwd('Desktop/VectorSim/scripts')

graphics.off()
rm(list = ls())

source("simulator_functions.R")
library(abc)
library(abind)
library(Rcpp)
library(RcppArmadillo)
library(mc2d)
library(rdist)
library(bestNormalize)
sourceCpp(file = "propegate.cpp")

propegate_population_inner = propegate_population_inner_C
propegate_population = propegate_population_decomposed



############################
# Migration induced biases #
############################

# check estimates for m = 0 case
{
loci = 1000
pop = initialise_population(n_subpops = 2, loci = loci)
sampling_intervals = c(2,10,15)
recom_rates = c(.5, .35, .2, .05)
m = 0.02
Mfn = construct_M_4_step
N = 4000
grid_size = 121
sample_subpops = c(61)


max_loci = loci * 1.1
min_loci  = loci
n_loci = loci

S_vec = c(Inf, 1000)

model_fn <- get_model_fn(sampling_intervals, recom_rates, Mfn, S_vec, grid_size, sample_subpops, max_loci, min_loci, n_loci)
Rprof()
o <- model_fn(N, m)
Rprof(NULL)
summaryRprof()
sampling_intervals/(2*colMeans(o$Fc[,1,]))
(1 - rowMeans(o$r2[,1,]) ) / (2 * recom_rates * (2-recom_rates) * rowMeans(o$r2[,1,]))
}


# calculate estimates for many methods, many deme sizes, many m values
{
loci = 5000
N = 100000
Mfn = construct_M_4_step
ms <- rev(c(.01, .05, .1, .2, .3, .4, .5, .6, .7, .8))
n_demes_v <- c(9,25,49,81,121, 169)
sampling_intervals = c(2,3,5,10,15, 20)
recom_rates = c(.5, .35, .2, .05)
AF_ests = array(NA, dim = c(length(ms), length(n_demes_v), length(sampling_intervals)))
LD_ests = array(NA, dim = c(length(ms), length(n_demes_v), length(recom_rates) ))
AF_raw = array(NA, dim = c(length(ms), length(n_demes_v), length(sampling_intervals)))
LD_raw = array(NA, dim = c(length(ms), length(n_demes_v), length(recom_rates) ))


i = 1; j = 1; kAF = 1; kLD = 1
for (grid_size in n_demes_v){
  sample_subpops = c(grid_size/2 + .5)  
  model_fn <- get_model_fn(sampling_intervals, recom_rates, Mfn, S_vec = c(Inf), grid_size, sample_subpops, max_loci, min_loci, n_loci)
  i = 1
  for (m in ms){
    o <- model_fn(N, m)
    
    Fc_v = colMeans(o$Fc[,1,])
    r2_v = rowMeans(o$r2[,1,])
    
    AF_ests[i,j,] = sampling_intervals/(2*Fc_v)/N
    LD_ests[i,j,] = (1-r2_v)/(2*recom_rates*(2-recom_rates)*r2_v)/N
    
    AF_raw[i,j,] = Fc_v
    LD_raw[i,j,] = r2_v
    i = i+1
  }
  
  j = j+1
  gc()
}

saveRDS(AF_ests, '../data/AF_ests_demo.rds')
saveRDS(LD_ests, '../data/LD_ests_demo.rds')

}

AF_ests = readRDS('../data/AF_ests_demo.rds')
LD_ests = readRDS('../data/LD_ests_demo.rds')
## make line plots demonstrating bias
### plot a particular deme size
{
size_idx = 3
size_selected = n_demes_v[size_idx]

AF_ests_size = AF_ests[,size_idx,]
LD_ests_size = LD_ests[,size_idx,]

pdf('../final_figures/bias.pdf', height = 4, width = 6, pointsize = 8)

par(mfrow = c(1,1), mai = c(.7,.9,.4,.4))
plot(NA, bty = 'l', type = 'l', yaxs = 'i', xaxs = 'i', ylim = log2(range(1,cbind(AF_ests_size, LD_ests_size),20)), xlim = range(ms), xlab = 'm', ylab = expression(hat(N[e])/N[e]), yaxt = 'n', xaxt = 'n'); abline(h = 0); axis(side = 1, at=ms, labels = ms); axis(side=2, at=c(0,1,2,3,4), labels = c(1,2,4,8,16))
 

colfunc_AF <- colorRampPalette(c("red", "green"))
cols = colfunc_AF(length(sampling_intervals))
for (col in 1:length(sampling_intervals)){
  lines(ms, log2(AF_ests_size[,col]), col = cols[col], cex = 1.5, lw = .7)
}

for (col in 1:length(recom_rates)){
  lines(ms, log2(LD_ests_size[,col]), lty=col, lw = 1)
  
}


leg = c(paste('AF: T = ', sampling_intervals), paste('LD: c = ', recom_rates), '','')

col = c(cols, rep('black', length(recom_rates)), NA, NA)
lty = c(rep(1, length(sampling_intervals)), 1:length(recom_rates), NA, NA)

  legend(0.02, 4.2, legend = leg, col = col, lty = lty, cex = 0.85, nc = 2)
  
  graphics.off()
  
  
}



{
  size_idx = 3
  size_selected = n_demes_v[size_idx]
  
  AF_ests_size = AF_ests[,size_idx,]
  LD_ests_size = LD_ests[,size_idx,]
  
  pdf('../final_figures/bias2.pdf', height = 4.5, width = 7, pointsize = 8)
  
  par(mfrow = c(1,1), mai = c(.7,.9,.4,.4))
  plot(NA, bty = 'l', type = 'l', yaxs = 'i', xaxs = 'i', ylim = log2(range(1,cbind(AF_ests_size, LD_ests_size),20)), xlim = range(ms), xlab = 'm', ylab = bquote(log[2]*(frac(hat(N[e])^LD0.05, N[e])) - 'log'[2](frac(hat(N[e]), N[e] ) )), yaxt = 'n', xaxt = 'n'); abline(h = 0); axis(side = 1, at=ms, labels = ms); axis(side=2, at=c(0,1,2,3,4))
  
  
  colfunc_AF <- colorRampPalette(c("red", "green"))
  cols = colfunc_AF(length(sampling_intervals))
  for (col in 1:length(sampling_intervals)){
    lines(ms, log2(LD_ests_size[,4]) - log2(AF_ests_size[,col]), col = cols[col], cex = 1.5, lw = .7)
  }
  
  for (col in 1:length(recom_rates)){
    lines(ms, log2(LD_ests_size[,4]) - log2(LD_ests_size[,col]), lty=col, lw = 1)
    
  }
  
  leg = c(paste('AF: t =', sampling_intervals), paste('LD: c = ', recom_rates), '', '')
  col = c(cols, rep('black', length(recom_rates)), NA, NA)
  lty = c(rep(1, length(sampling_intervals)), 1:length(recom_rates), NA, NA)
  
  legend(0.02, 4.2, legend = leg, col = col, lty = lty, cex = 0.9, nc = 2)
  
  graphics.off()
  
  
}

### plot a particular method
{
method_idx_AF = 6
method_idx_LD = 4

AF_ests_method = AF_ests[,,method_idx_AF]
LD_ests_method = LD_ests[,,method_idx_LD]
pdf('../final_figures/finite_size.pdf',width = 5, height = 3.2, pointsize = 5 )
par(mfrow = c(1,2), mai = c(.6,.6,0.1,0.1), oma = c(.7,0,0,0))
plot(NA, type = 'l',xaxs = 'i', log = 'x', ylim = c(1,20 ), xlim = range(ms), xlab = '',xaxt = 'n', ylab = expression(hat(N[e]))); abline (h = 1); axis(1, at =c(0.01,0.05,0.2,0.4,0.8)); axis(2, at =c(1,5,10,15,20))
mtext('m', side = 1, line = 2.2)

fig_label_logx('A', cex = 2)
colfunc_AF <- colorRampPalette(c("red", "blue"))
cols = c(colfunc_AF(length(n_demes_v)-1), 'black')
for (col in 1:length(n_demes_v)){
  lines(ms, AF_ests_method[,col], col = cols[col], cex = 1.5)
}
leg = n_demes_v
col = c(cols)
lty = c(rep(1, length(n_demes_v)))
#legend('topleft', legend = leg, col = col, lty = lty, title = 'demes')


plot(NA, type = 'l',xaxs = 'i', log = 'x', ylim = c(1,20), xlim = range(ms), xlab = '',xaxt = 'n', ylab = expression(hat(N[e]))); abline(h = 1); axis(1, at =c(0.01,0.05,0.2,0.4,0.8))
mtext('m', side = 1, line = 2.2)
fig_label_logx('B', cex = 2)

for (col in 1:length(n_demes_v)){
  lines(ms, LD_ests_method[,col], col=cols[col], cex = 1.5)
  
}
leg = n_demes_v
col = c(cols)
lty = c(rep(1, length(n_demes_v)))

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('bottom',inset = 0.02, legend = leg, col = col, lty = lty, horiz = T)
graphics.off()



}
{
  method_idx_AF = 6
  method_idx_LD = 4
  
  AF_ests_method = AF_ests[,,method_idx_AF]
  LD_ests_method = LD_ests[,,method_idx_LD]
  
  par(mfrow = c(2,1), mai = c(1,1,0.3,0.1))
  plot(NA, type = 'l', log = 'x', ylim = c(-1,1 ), xlim = range(ms), xlab = 'm', ylab = expression(hat(N[e])), main = paste("AF, separation =", sampling_intervals[method_idx_AF])); abline (h = 1)
  
  colfunc_AF <- colorRampPalette(c("red", "green"))
  cols = colfunc_AF(length(n_demes_v))
  for (col in 1:length(n_demes_v)){
    lines(ms, AF_ests_method[,4]-AF_ests_method[,col], col = cols[col], cex = 1.5)
  }
  leg = n_demes_v
  col = c(cols)
  lty = c(rep(1, length(n_demes_v)))
  legend('topright', legend = leg, col = col, lty = lty, title = 'demes')
  
  
  plot(NA, type = 'l', log = 'x', ylim = c(-1,2), xlim = range(ms), xlab = 'm', ylab = expression(hat(N[e])), main = paste("LD, recom rate = ", recom_rates[method_idx_LD])); abline(h = 1)
  for (col in 1:length(n_demes_v)){
    lines(ms, LD_ests_method[,4]-LD_ests_method[,col], col=cols[col], cex = 1.5)
    
  }
  leg = n_demes_v
  col = c(cols)
  lty = c(rep(1, length(n_demes_v)))
  legend('topright', legend = leg, col = col, lty = lty, title = 'demes')
  
}

### plot a particular method, particular size
{
method_idx_AF = 5
method_idx_LD = 4
size_idx = 5

AF_ests_method = AF_ests[,size_idx,method_idx_AF]
LD_ests_method = LD_ests[,size_idx,method_idx_LD]
size_selected = n_demes_v[size_idx]

par(mfrow = c(2,1), mai = c(1,1,0.3,0.1))
plot(NA, type = 'l', log = 'x', ylim = c(.5,max(AF_ests_method) ), xlim = range(ms), xlab = 'm', ylab = expression(hat(N[e])), main = paste("AF, separation =", sampling_intervals[method_idx_AF], ', n_demes = ', size_selected ) ); abline (h = 1)
lines(ms, AF_ests_method, cex = 1.5)

plot(NA, type = 'l', log = 'x', ylim = c(.5,max(LD_ests_method)), xlim = range(ms), xlab = 'm', ylab = expression(hat(N[e])), main = paste("LD, recom rate = ", recom_rates[method_idx_LD], ', n_demes = ', size_selected)); abline(h = 1)
lines(ms, LD_ests_method, cex = 1.5)
}


## make line plots showing Fc and r2 values
### plot a particular deme size
{
size_idx = 2
size_selected = n_demes_v[size_idx]

AF_ests_size = AF_raw[,size_idx,]
LD_ests_size = LD_raw[,size_idx,]

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

### plot a particular method
{
method_idx_AF = 6
method_idx_LD = 4

AF_ests_method = AF_raw[,,method_idx_AF]
LD_ests_method = LD_raw[,,method_idx_LD]

par(mfrow = c(2,1), mai = c(1,1,0.3,0.1))
plot(NA, type = 'l', log = 'x', ylim = range(AF_ests_method) , xlim = range(ms), xlab = 'm', ylab = expression(hat(N[e])), main = paste("AF, separation =", sampling_intervals[method_idx_AF])); abline (h = 1)

colfunc_AF <- colorRampPalette(c("red", "green"))
cols = colfunc_AF(length(n_demes_v))
for (col in 1:length(n_demes_v)){
  lines(ms, AF_ests_method[,col], col = cols[col], cex = 1.5)
}
leg = n_demes_v
col = c(cols)
lty = c(rep(1, length(n_demes_v)))
legend('topleft', legend = leg, col = col, lty = lty, title = 'demes')


plot(NA, type = 'l', log = 'x', ylim = range((LD_ests_method)), xlim = range(ms), xlab = 'm', ylab = expression(hat(N[e])), main = paste("LD, recom rate = ", recom_rates[method_idx_LD])); abline(h = 1)
for (col in 1:length(n_demes_v)){
  lines(ms, LD_ests_method[,col], col=cols[col], cex = 1.5)
  
}
leg = n_demes_v
col = c(cols)
lty = c(rep(1, length(n_demes_v)))
legend('topleft', legend = leg, col = col, lty = lty, title = 'demes')
}

### plot a particular method, particular size,  # not implemented
{
method_idx_AF = 5
method_idx_LD = 4
size_idx = 5

AF_ests_method = AF_raw[,size_idx,method_idx_AF]
LD_ests_method = LD_raw[,size_idx,method_idx_LD]
size_selected = n_demes_v[size_idx]

par(mfrow = c(2,1), mai = c(1,1,0.3,0.1))
plot(NA, type = 'l', log = 'x', ylim = range(AF_ests_method), xlim = range(ms), xlab = 'm', ylab = expression(Fc), main = paste("AF, separation =", sampling_intervals[method_idx_AF], ', n_demes = ', size_selected ) );
lines(ms, AF_ests_method, cex = 1.5)

plot(NA, type = 'l', log = 'x', ylim = range((LD_ests_method)), xlim = range(ms), xlab = 'm', ylab = expression(r^2), main = paste("LD, recom rate = ", recom_rates[method_idx_LD], ', n_demes = ', size_selected)); abline(h = 1)
lines(ms, LD_ests_method, cex = 1.5)


colfunc_AF <- colorRampPalette(c("red", "green"))
cols = colfunc_AF(length(n_demes_v))
for (col in 1:length(n_demes_v)){
  
}
leg = n_demes_v
col = c(cols)
lty = c(rep(1, length(n_demes_v)))
legend('topleft', legend = leg, col = col, lty = lty, title = 'demes')


for (col in 1:length(n_demes_v)){
  lines(ms, LD_ests_method[,col], col=cols[col], cex = 1.5)
  
}
leg = n_demes_v
col = c(cols)
lty = c(rep(1, length(n_demes_v)))
legend('topleft', legend = leg, col = col, lty = lty, title = 'demes')
}


###################
# ABC experiments #
###################

library(foreach)
library(parallel)
library(doParallel)




# generate simualted data set
{
loci = 1000; max_loci = ( loci + 5 ) * 1.1; min_loci = loci; n_loci = loci
  sample_gens = c(0,5,10,15,20,25)
  recom_rates = c(.50, .35, .20, .05)

Mfn = construct_M_4_step
grid_size = 100
sample_subpops = c(41)

S_vec = c(Inf, 100, 500, 1000, 400, 2000, 4000, 900, 4500, 9000)
n_par = 6
n_samp =10

# constructing some things
{
  
model_fn_list = list()
Nm_list = list()
for (i in 1:n_par){
  model_fn_list[[i]] = get_model_fn(sampling_intervals, recom_rates, Mfn, S, grid_size, sample_subpops, max_loci, min_loci, n_loci)
  Nm_list[[i]] = Nm_draw_large(n_samp)
}
}


s = proc.time()
sim_out <- generate_ABC_sim_parallel(loci, sampling_intervals, recom_rates, Nm_list, model_fn_list, n_par, S_vec, Mfn, grid_size, sample_subpops)
e = proc.time(); print(e-s)

#saveRDS(sim_out, '../data/sim_out_mega.rds')
sim_out <- readRDS('../data/sim_out_mega.rds')

Nm_values_sim = sim_out$Nm_values
Nm_values_sim = cbind(Nm_values_sim, apply(Nm_values_sim,1,prod))
}




# generate 'real data' set
{
  
# have a component with varying N
  
N1 = seq(2000, 20000, 500)
m1 = rep(0.05, 37)
  
# have a component with varying m in [0.01, 1.0], N = 10,000

N2 = rep(10000, 19)
m2 = rev(seq(0.05, 0.95, 0.05))
  
# have a component with varying m in [0.01, 1.0], N = 5,000

N3 = rep(5000, 19)
m3 = rev(seq(0.05, 0.95, 0.05))


# have a component with varying m in [0.01, 1.0], N = 15,000

N4 = rep(15000, 19)
m4 = rev(seq(0.05, 0.95, 0.05))
  
# have a component with varying m in [0, 0.05], much denser, N = 10,000

N5 = rep(10000, 26)
m5 = rev(seq(0, 0.05, 0.002))

Ns = c(N1, N2, N3, N4, N5)
ms = c(m1, m2, m3, m4, m5)


grid_size = 81
sample_subpops = c(41)

# constructing some things
{
Ns_l = split(Ns,ceiling(seq_along(Ns)/20))
ms_l = split(ms,ceiling(seq_along(ms)/20))

model_fn_list = list()
Nm_list = list()
for (i in 1:n_par){
  model_fn_list[[i]] = get_model_fn(sampling_intervals, recom_rates, Mfn, S, grid_size, sample_subpops, max_loci, min_loci, n_loci)
  Nm_list[[i]] = matrix(c(Ns_l[[i]], ms_l[[i]]), nc = 2)
}
}

s = proc.time()
real_out <- generate_ABC_sim_parallel(loci, sampling_intervals, recom_rates, Nm_list, model_fn_list, n_par, S_vec, Mfn, grid_size, sample_subpops)
e = proc.time(); print(e-s)


saveRDS(real_out, '../data/real_out_more.rds')
real_out <- readRDS('../data/real_out_more.rds')

Nm_values_real = real_out$Nm_values
Nm_values_real = cbind(Nm_values_real, apply(Nm_values_real,1,prod))

}
# do ABC



# construct some posteriors

### construct summary statistics
{
### s arrays: Ne_hat, and mean_CV
S = Inf
SS_list_sim <- format_SS(sim_out$Fc_sim, sim_out$r2_sim, Ss = c(S), gs = c(20))

cs = as.numeric(sub('.','', dimnames(SS_list_sim$r2)[[2]]))
ts = as.numeric(sub('.','', dimnames(SS_list_sim$Fc)[[2]]))

cleaned_SS_list_sim = make_SS_vecs(SS_list_sim$Fc, SS_list_sim$r2)
cs_ests_sim = do_construction_LD(cleaned_SS_list_sim$r2, to_Ne = T, cs = cs)
ts_ests_sim = do_construction_AF(cleaned_SS_list_sim$Fc, to_Ne = T, ts = ts, S = S)

full_SS_sim = cbind(cs_ests_sim[,], ts_ests_sim[,])
full_SS_sim_mean_CV = cbind(rowMeans(full_SS_sim), sqrt(rowVars(full_SS_sim))/rowMeans(full_SS_sim) )



SS_list_real <- format_SS(real_out$Fc_sim, real_out$r2_sim, Ss = c(S), gs = c(20))

cs = as.numeric(sub('.','', dimnames(SS_list_real$r2)[[2]]))
ts = as.numeric(sub('.','', dimnames(SS_list_real$Fc)[[2]]))

cleaned_SS_list_real = make_SS_vecs(SS_list_real$Fc, SS_list_real$r2)
cs_ests_real = do_construction_LD(cleaned_SS_list_real$r2, to_Ne = T, cs = cs)
ts_ests_real = do_construction_AF(cleaned_SS_list_real$Fc, to_Ne = T, ts = ts, S = S)

full_SS_real = cbind(cs_ests_real[,], ts_ests_real[,])
full_SS_real_mean_CV = cbind(rowMeans(full_SS_real), sqrt(rowVars(full_SS_real))/rowMeans(full_SS_real) )



}

### Coefficient of variation  plots
{
### 
plot(Nm_values_sim[,2], full_SS_sim_mean_var[,2])

### Coefficient of variation for specific true N
d_filt = (Nm_values_sim[,1] > 10000 & Nm_values_sim[,1] < 10500)
d = cbind(Nm_values_sim[d_filt, ], full_SS_sim_mean_CV[d_filt, ])

par(mfrow = c(1,1))
plot(Nm_values_sim[,2], full_SS_sim_mean_var[,2], ylab = 'CV( estimates )', xlab = 'm', pch = 16)#, col = cols(dim(d)[1])[sort(d[,1], index.return = T )$ix] )# /rowMeans(d[,3:6]) )
plot(d[,2], d[,4], ylab = 'CV( estimates )', xlab = 'm')#, col = cols(dim(d)[1])[sort(d[,1], index.return = T )$ix] )# /rowMeans(d[,3:6]) )
}




### look at poisterior for a specific 'real' SS vector
n = 119; print(round(Nm_values_real[n,], 3))
par(mfrow = c(2,2), ask = F)

o <- do_abc(full_SS_sim, full_SS_real[n,], Nm_values_sim, N = Nm_values_real[n,1], m = Nm_values_real[n,2], draw = T, method = 'neuralnet', tol = .1, pca = F, transformation = rep('none',3))
plot(o, Nm_values_sim, ask = F)


# cross validation nnet with tol ~ 0.1 is best
{
  tols = c(0.01, 0.02,0.05,.1, 0.2)
  
cv.nnet = cv4abc(Nm_values_sim, sumstat = full_SS_sim, nval = 200, statistic = 'mean', tols = c(0.01, 0.02,0.05,.1, 0.2), method = 'neuralnet')
plot(cv.nnet)
summary(cv.nnet)

cv.lin = cv4abc(Nm_values_sim, sumstat = full_SS_sim, nval = 200, statistic = 'mean', tols = c(0.01, 0.02,0.05,.1, 0.2), method = 'loclinear')
plot(cv.lin, ask = F)
summary(cv.lin)

cv.rej = cv4abc(Nm_values_sim, sumstat = full_SS_sim, nval = 200, statistic = 'mean', tols = c(0.01, 0.02,0.05,.1, 0.2), method = 'rejection')
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

# examine accepted and rejected points
{
par(mfrow = c(1,1))
plot(Nm_values_sim[,1], Nm_values_sim[,2], xlab = 'N', ylab = 'm', pch = 16, cex = .5, col = 'grey', xlim = c(2000,20000), ylim = c(0,1), xaxs='i', yaxs = 'i')#, pch = 16)

points(o$unadj.values[,1], o$unadj.values[,2], pch = 16, col = 'red', cex = .7)
points(o$adj.values[,1], o$unadj.values[,2], pch = 16, col = 'blue', cex = .7)
points(Nm_values_real[n,1],Nm_values_real[n,2], pch = 8, cex = 2, col = 'black', lw = 3)
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



upperLests =  colMeans(full_SS_sim[upperLidx,])
lowerRests = colMeans(full_SS_sim[lowerRidx,])

somd = cbind(round(upperLests,-2),round(lowerRests,-2))
colnames(somd) = c('upper L', 'lower R')
somd

pc <- prcomp(full_SS_sim, scale = T, center= T)

plot(pc$x[,1], pc$x[,2], xlab = 'PC1', ylab = 'PC2', cex = 1)
points(pc$x[upperLidx,1], pc$x[upperLidx,2], col = 'red', pch = 16)
points(pc$x[lowerRidx,1], pc$x[lowerRidx,2], col = 'blue', pch = 16)

legend('topleft' , col = c('red', 'blue'), legend = c('upper L', 'lower R'), pch = 16)

pal = colorRampPalette(c("blue", "red"))(length(pc$x[,1]))
plot(pc$x[,1], pc$x[,2], xlab = 'PC1', ylab = 'PC2', col = pal[order(Nm_values_sim[,2])], cex = 2*range01(Nm_values_sim[,1]), pch = 1)



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
  o.nn <- do_abc(full_SS_sim, full_SS_real[r,], Nm_values_sim, method = 'neuralnet', tol = .1, pca = F)
  summary_o.nn = summary(o.nn)
  mean_v_N.nn = rbind(mean_v_N.nn, c(summary_o.nn[4,], prod(summary_o.nn[4,])) )
  var_v_N.nn = rbind(var_v_N.nn, c(colVars(o.nn$adj.values), var(apply(o.nn$adj.values, 1, prod)) ) )
  
  o.rej <- do_abc(full_SS_sim, full_SS_real[r,], Nm_values_sim, method = 'rejection', tol = .02, pca = F)
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
  o.nn <- do_abc(full_SS_sim, full_SS_real[r,], Nm_values_sim, method = 'neuralnet', tol = .1, pca = F)
  summary_o.nn = summary(o.nn)
  mean_v_m.nn = rbind(mean_v_m.nn, c(summary_o.nn[4,], prod(summary_o.nn[4,])) )
  var_v_m.nn = rbind(var_v_m.nn, c(colVars(o.nn$adj.values), var(apply(o.nn$adj.values, 1, prod)) ) )
  
  o.rej <- do_abc(full_SS_sim, full_SS_real[r,], Nm_values_sim, method = 'rejection', tol = .1, pca = F)
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
      o.nn <- do_abc(full_SS_sim, full_SS_real[r,], Nm_values_sim, method = 'neuralnet', tol = .1, pca = F)
      summary_o.nn = summary(o.nn)
      mean_v_m.nn = rbind(mean_v_m.nn, c(summary_o.nn[4,], prod(summary_o.nn[4,])) )
      var_v_m.nn = rbind(var_v_m.nn, c(colVars(o.nn$adj.values), var(apply(o.nn$adj.values, 1, prod)) ) )
      
      o.rej <- do_abc(full_SS_sim, full_SS_real[r,], Nm_values_sim, method = 'rejection', tol = .1, pca = F)
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
  o.nn <- do_abc(full_SS_sim, full_SS_real[r,], Nm_values_sim, method = 'neuralnet', tol = .1, pca = F)
  summary_o.nn = summary(o.nn)
  mean_v_m.nn = rbind(mean_v_m.nn, c(summary_o.nn[4,], prod(summary_o.nn[4,])) )
  var_v_m.nn = rbind(var_v_m.nn, c(colVars(o.nn$adj.values), var(apply(o.nn$adj.values, 1, prod)) ) )
  
  o.rej <- do_abc(full_SS_sim, full_SS_real[r,], Nm_values_sim, method = 'rejection', tol = .1, pca = F)
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



###########

# want to investigate:

### different sample sizes - trade off between sample size and whether we have LD data

### samples from multipl epopulations
###### what summary statistics would be useful?
###### how to investigate whether they are useful easily? # perhaps correlation with parameter values of interest # correlation in Fc for example

### simulating a section of chromosome with simuPOP
###### can plug in
