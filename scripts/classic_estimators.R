library("mc2d")
library("ggplot2")
source("simulator_functions.R")

# WAPLES SINGLE LOCUS RESULTS

# no migration
fcmean = c()
ests = c()
Ns = seq(500, 5000, 100)
S1 = S2 = 50
t = 99
gs = c(1,4)
loci = 1e4


for (N in Ns){
  print(N)
o <- sim.drift2(p0 = 0.5, loci = loci, N = N, t = t)

s1 <- take_sample(o[,gs[1]], S1)
s2 <- take_sample(o[,gs[2]], S2)
fcv <- get_fc_vec(s1,s2)
fcmean = c(fcmean, mean(fcv))
AF <- get_AF_est(s1, s2, S1, S2, diff(gs))

AF_full <- get_AF_est(o[,gs[1]], o[,gs[2]], Inf, Inf, 3)

.5/AF
ests = c(ests,.5/AF_full)
}

plot(Ns, ests, xlab = 'true N', ylab = expression(hat(N[AF]))); abline(0,1, col = 'red')
plot(Ns, fcmean, xlab = 'true N', ylab = expression(hat(N[AF]))); abline(0,1, col = 'red')



####### some plotting
par(mfrow = c(1,1))
plot(NA, xlim = c(1,100), ylim = c(.25,.75), xlab = 't', ylab = 'Allele frequency')
for (i in 1:5){
  lines(1:100, o[i,], type = 'l')
}


N = 5000
S1 = S2 = 50
t = 99
gs = c(1,4)
loci = 5
o <- sim.drift2(p0 = 0.5, loci = loci, N = N, t = t)

for (i in 1:5){
  lines(1:100, o[i,], type = 'l', col = 'red')
}

N = 50000
S1 = S2 = 50
t = 99
gs = c(1,4)
loci = 5

o <- sim.drift2(p0 = 0.5, loci = loci, N = N, t = t)
for (i in 1:5){
  lines(1:100, o[i,], type = 'l', col = 'blue')
}


leg = c('500', '5000', '50000')
col = c('black', 'red', 'blue')
legend('topleft', legend = leg, col = col, lty = 1, cex = 1.3, title = expression(N[e]))










##########



# with migration

n_subpops = 10
m=.2

N_v = rep(N, n_subpops)
M <- construct_M(m = m, n_subpops = n_subpops)

o <- sim.drift.migration(p0=0.5, loci=loci, N = N_v, M = M, t=t+1)

# pre migration samples
pre_s1 <- take_sample(o$pre[,gs[1],1], S1)
pre_s2 <- take_sample(o$pre[,gs[2],1], S2)
# post migration samples
post_s1 <- take_sample(o$post[,gs[1],1], S1)
post_s2 <- take_sample(o$post[,gs[2],1], S2)


AF_pre <- get_AF_est(pre_s1, pre_s2, S1, S2, diff(gs))
AF_post <- get_AF_est(post_s1, post_s2, S1, S2, diff(gs))

AF_pre_full <- get_AF_est(o$pre[,gs[1],1],o$pre[,gs[2],1], Inf, Inf, 3)
AF_post_full <- get_AF_est(o$post[,gs[1],1],o$post[,gs[2],1], Inf, Inf, 3)


.5/AF_pre
.5/AF_post
.5/AF_pre_full
.5/AF_post_full

###############################
# SIMULATE TWO LOCI
###############################

### HAPLOTYPES - NO MIGRATION ###
ests= c()
Ns = seq(500, 5000,100)
S = 50
S_gens = c(497,500)
S_gens = c(97,100)
loci = 1000

pA_init = 0.5
pB_init = 0.5
p0 = c(pA_init*pB_init, pA_init*(1-pB_init), (1-pA_init)*pB_init, (1-pA_init)*(1-pB_init))
c = .4

for (N in Ns){
  print(N)
o <- sim.drift5(N=N, t=S_gens[2], loci = loci, p0=p0, c = c)


plot(1:(max(S_gens)+1), colMeans(o$r, na.rm=T)); abline(v=min(S_gens))

# LD estimate
s1 <- take_sample_r(o$pAB[,S_gens[1]], o$pAb[,S_gens[1]], o$paB[,S_gens[1]], o$pab[,S_gens[1]], S)
LD_est_sample <- get_LD_est(s1$pAB, s1$pAb, s1$paB, s1$pab, S, c)
LD_est_full <- get_LD_est(o$pAB, o$pAb, o$paB, o$pab, Inf, c)

LD_est_sample
ests = c(ests,LD_est_full)
}

plot(Ns, ests, xlab = 'true N', ylab = expression(hat(N[LD]))); abline(0,1, col = 'red')

# AF estimate
s2 <- take_sample_r(o$pAB[,S_gens[2]], o$pAb[,S_gens[2]], o$paB[,S_gens[2]], o$pab[,S_gens[2]], S)
AF_est <- get_AF_est(s1$pA, s2$pA, S, S, diff(S_gens))
.5/AF_est


########### some plotting

par(mfrow = c(1,1))
plot(NA, xlim = c(1,101), ylim = c(-0.07,0.07), xlab = 't', ylab = 'D')
for (i in 1:loci){
  lines(1:101, o$D[i,])
}

N = 5000
o <- sim.drift5(N=N, t=S_gens[2], loci = loci, p0=p0, c = c)


for (i in 1:loci){
  lines(1:101, o$D[i,], col = 'red')
}

N = 100000
o <- sim.drift5(N=N, t=S_gens[2], loci = loci, p0=p0, c = c)
for (i in 1:loci){
  lines(1:101, o$D[i,], col = 'blue')
}



leg = c('500', '5000', '100000')
col = c('black', 'red', 'blue')
legend('topleft', legend = leg, col = col, lty = 1, cex = 1.3, title = expression(N[e]) )






########### some plotting












##########


# varying c

N = 100
S = 50
gen = 100
loci = 1e2

pA_init = 0.5
pB_init = 0.5
p0 = c(pA_init*pB_init, pA_init*(1-pB_init), (1-pA_init)*pB_init, (1-pA_init)*(1-pB_init))
cs = c(.05,.1,.2,.3,.4,.5)

ests <- c()
ests_f <- c()

for (c in cs){
  print(c)
  o <- sim.drift5(N=N, t=gen, loci = loci, p0=p0, c = c)
  
  # LD estimate
  s1 <- take_sample_r(o$pAB[,gen], o$pAb[,gen], o$paB[,gen], o$pab[,gen], S)
  LD_est_sample <- get_LD_est(s1$pAB, s1$pAb, s1$paB, s1$pab, S, c)
  LD_est_full <- get_LD_est(o$pAB[,gen], o$pAb[,gen], o$paB[,gen], o$pab[,gen], Inf, c)
  
  ests <- c(ests,LD_est_sample)
  ests_f <- c(ests_f,LD_est_full)

}

plot(cs, ests, type = 'l', ylim = c(40,110), xlim = c(0,0.5))
lines(cs, ests_f)
abline(h = c(N/2,N), lty = 2)



# varying number of loci


N = 100
S = Inf
gen = 100
loci_v = c(10,30,50,100,200)

pA_init = 0.5
pB_init = 0.5
p0 = c(pA_init*pB_init, pA_init*(1-pB_init), (1-pA_init)*pB_init, (1-pA_init)*(1-pB_init))
c = .5

reps = 1000
ests_LDm <- array(NA, dim = c(length(0, reps)))
ests_AFm <- c()
for (loci in loci_v){
  print(loci)
  
  ests_LD <- c()
  ests_AF <- c()
  for (r in 1:reps){
  o <- sim.drift5(N=N, t=gen, loci = loci, p0=p0, c = c)
  
  # LD estimate
  s1 <- take_sample_r(o$pAB[,gen], o$pAb[,gen], o$paB[,gen], o$pab[,gen], S)
  LD_est_sample <- get_LD_est(s1$pAB, s1$pAb, s1$paB, s1$pab, S, c=c)
  LD_est_full <- get_LD_est(o$pAB[,gen], o$pAb[,gen], o$paB[,gen], o$pab[,gen], Inf, c)
  
  ests_LD <- c(ests_LD,LD_est_sample)
  
  s2 <- take_sample_r(o$pAB[,gen-3], o$pAb[,gen-3], o$paB[,gen-3], o$pab[,gen-3], S)
  
  ests_AF <- c(ests_AF, get_AF_est(s1$pAB + s1$pAb, s2$pAB + s2$pAb, S, S, 3))
  }
  ests_LDm <- rbind(ests_LDm, ests_LD)
  ests_AFm <- rbind(ests_AFm, ests_AF)
  
}


plot(loci_v, rowMeans(ests_LDm), type = 'l', ylim = c(50,200), xlim = c(0,210))
lines(loci_v, (.5/rowMeans(ests_AFm)))
abline(h = c(N/2,N), lty = 2)

hist(.5/ests_AFm[1,], breaks = 100); abline(v=100)
hist(ests_LDm[1,], breaks = 100); abline(v=100)

### HAPLOTYPES- MIGRATION ###

# check with single example

n_subpops <- 10
m = .3

N_v <- rep(N, n_subpops)
M <- construct_M(m, n_subpops)

pA_init = 0.5
pB_init = 0.5
p0 = c(pA_init*pB_init, pA_init*(1-pB_init), (1-pA_init)*pB_init, (1-pA_init)*(1-pB_init))

# simulate both haplotypes and single locus - make sure we get the sampl estimates
o <- sim.drift5_mig(N=N_v, M=M, c=.5, t=S_gens[2]+1, loci=loci, p0=p0)
o_single <- sim.drift.migration(p0=pA_init, loci=loci, N = N_v, M = M, t=S_gens[2])



# PRE MIGRATION SAMPLES

# LD estimate
s1 <- take_sample_r(o$pAB_pre[,S_gens[1],1], o$pAb_pre[,S_gens[1],1], o$paB_pre[,S_gens[1],1], o$pab_pre[,S_gens[1],1], S)
LD_est <- get_LD_est(s1$pAB, s1$pAb, s1$paB, s1$pab, S)
2*LD_est

# AF estimate
s2 <- take_sample_r(o$pAB_pre[,S_gens[2],1], o$pAb_pre[,S_gens[2],1], o$paB_pre[,S_gens[2],1], o$pab_pre[,S_gens[2],1], S)
AF_est <- get_AF_est(s1$pA, s2$pA, S, S, 3)
.5/AF_est

# # with single locus simulation

s1_single <- take_sample(o_single$pre[,S_gens[1],1],S)
s2_single <- take_sample(o_single$pre[,S_gens[2],1],S)
AF_est_single <- get_AF_est(s1_single, s2_single, S, S, diff(S_gens))
.5/AF_est_single

# POST MIGRATION SAMPLES

# LD estimate
s1 <- take_sample_r(o$pAB_post[,S_gens[1],1], o$pAb_post[,S_gens[1],1], o$paB_post[,S_gens[1],1], o$pab_post[,S_gens[1],1], S)
LD_est <- get_LD_est(s1$pAB, s1$pAb, s1$paB, s1$pab, S)
2*LD_est

# AF estimate
s2 <- take_sample_r(o$pAB_post[,S_gens[2],1], o$pAb_post[,S_gens[2],1], o$paB_post[,S_gens[2],1], o$pab_post[,S_gens[2],1], S)
AF_est <- get_AF_est(s1$pA, s2$pA, S, S, diff(S_gens))
.5/AF_est

# # with single locus simulation

s1_single <- take_sample(o_single$post[,S_gens[1],1],S)
s2_single <- take_sample(o_single$post[,S_gens[2],1],S)
AF_est_single <- get_AF_est(s1_single, s2_single, S, S, diff(S_gens))
.5/AF_est_single


#################################################
#################################################
#
#    Symmetrical migration
#    &
#    Isalnd model
#
#################################################
#################################################

# population structure
N = 1000
S = Inf
n_subpops = 10
ms = c(0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.9, 1)
M_fn = construct_M_waples

p0 = rep(.25, 4)

# estimation
t_intervals = c(0,2,3,4,5,7,9)
cs = c(0.05, 0.2, 0.35, 0.5)


N_v = rep(N, n_subpops)


# simulation parameters
loci = 1e1
burn_in = 300 

reset_period = 10
sampling_repeats = 20
sampling_period = reset_period + max(t_intervals)

o <- burn_in_checker(ms, M_fn, cs, burn_in, N_v, loci, p0)

estimates_sym <- get_estimates_df(N_v = N_v,
                                  ms = ms,
                                  M_fn = M_fn,
                                  p0 = p0,
                                  loci = loci,
                                  burn_in = burn_in,
                                  t_intervals = t_intervals,
                                  rest_period = rest_period,
                                  sampling_repeats = sampling_repeats,
                                  sample_subpops = 10,
                                  S = S,
                                  cs = cs)



means_sym <- get_means_df(estimates_sym)
max_m = 1
means_sym_plot <- trim_df(means_sym, max_m = max_m)

colfunc_AF <- colorRampPalette(c("red", "orange"))
ltys <- c(2,4,5,6)

# plot
ggplot(data=means_sym_plot,
       aes(x=m, y=est, colour=method, linetype = method)) +
  scale_color_manual( values=c(colfunc_AF(length(t_intervals)-1), rep("black",length(cs))) ) +
  scale_linetype_manual( values=c(rep(1, length(t_intervals)-1), ltys[1:(length(cs))] ) ) +
  geom_line() + theme_bw() + scale_x_continuous(trans = "log10", breaks = ms) + 
  scale_y_continuous(expression(hat(Ne)), breaks =seq(0,1.05*max(means_sym_plot$est),1000)) + ggtitle(paste0("Symmetric migration amongst ", n_subpops, " populations")) #+ ylim(c(900,10100))




### repeat for immigration from large pop

# N, S, ms, cs, t_intervals, loci, and all simulation parameters are held the same

# population structure

N_large  = 1e9
N_v = c(N, N_large)

S = 1000
#ms = c(0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.9)
M_fn = construct_M_inf
#p0 = rep(.25, 4)

# estimation
#t_intervals = c(0,2,3,4,5,7,9)
#cs = c(0.05, 0.2, 0.35, 0.5)


# simulation parameters
#loci = 1e2
#burn_in = 500 

#reset_period = 10
#sampling_repeats = 20
#sampling_period = reset_period + max(t_intervals)


loci = 1e2


o <- burn_in_checker(ms, M_fn, cs, burn_in, N_v, loci, p0)

estimates_inf <- get_estimates_df(N_v = N_v,
                                  ms = ms,
                                  M_fn = M_fn,
                                  p0 = p0,
                                  loci = loci,
                                  burn_in = burn_in,
                                  t_intervals = t_intervals,
                                  rest_period = rest_period,
                                  sampling_repeats = sampling_repeats,
                                  sample_subpops = 1,
                                  S = S,
                                  c = cs)
                     
means_inf <- get_means_df(estimates_inf)

max_m <- 0.3
means_inf_plot <- trim_df(means_inf, max_m = max_m)


ggplot(data=means_inf_plot,
       aes(x=m, y=est, colour=method, linetype = method)) +
  scale_color_manual( values=c(colfunc_AF(length(t_intervals)-1), rep("black",length(cs))) ) +
  scale_linetype_manual( values=c(rep(1, length(t_intervals)-1), ltys[1:(length(cs))] ) ) +
  geom_line() + theme_bw() + scale_x_continuous(trans = "log10", breaks = ms) + 
  scale_y_continuous(expression(hat(Ne)), breaks =seq(0,1.05*max(means_inf_plot$est),1000)) + ggtitle("Immigration from population of size 10^8") #+ ylim(c(900,10100))










########################################################
########################################################
########################################################
#
#
#    STEPPING STONE MODELS 
#
#
########################################################
########################################################
########################################################


edge_lengths <- c(3, 7, 11, 15, 19)
edge_lengths = c(11)
means_list <- list()

for (edge_length in edge_lengths){
print(paste0("Edge Length = ", edge_length))

# population structure
N = 1000
S = Inf
n_subpops = edge_length**2
sample_subpops <- c(n_subpops/2 + .5)

ms = c(0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.9, 1)
M_fn = construct_M_4_step
p0 = rep(.25, 4)

# estimation
t_intervals = c(0,2,3,4,5,7,9,15,30)
cs = c(0.05, 0.2, 0.35, 0.5)


N_v = rep(N, n_subpops)


# simulation parameters
loci = 1e1
burn_in = 500

reset_period = 10
sampling_repeats = 20
sampling_period = reset_period + max(t_intervals)

o <- burn_in_checker(ms, M_fn, cs, burn_in, N_v, loci, p0)


estimates_sym <- get_estimates_df(N_v = N_v,
                                  ms = ms,
                                  M_fn = M_fn,
                                  p0 = p0,
                                  loci = loci,
                                  burn_in = burn_in,
                                  t_intervals = t_intervals,
                                  rest_period = rest_period,
                                  sampling_repeats = sampling_repeats,
                                  sample_subpops = sample_subpops,
                                  S = S,
                                  cs = cs)


means_sym <- get_means_df(estimates_sym)
s <- paste0('EL', edge_length)
means_list[[s]] <- means_sym
}

######## 
saveRDS(means_list, file = "../data/step_data2.rds")
means_list <- readRDS("../data/step_data.rds")

# PLOT

# A particular deme size

means_examine = means_list$EL11
max_m = 1
means_examine_plot <- trim_df(means_examine, max_m = max_m)


colfunc_AF <- colorRampPalette(c("red", "orange"))
ltys <- c(2,4,5,6)

ggplot(data=means_examine_plot,
       aes(x=m, y=est, colour=method, linetype = method)) +
  scale_color_manual( values=c(colfunc_AF(length(t_intervals)-1), rep("black",length(cs))) ) +
  scale_linetype_manual( values=c(rep(1, length(t_intervals)-1), ltys[1:(length(cs))] ) ) +
  geom_line() + theme_bw() + scale_x_continuous(trans = "log10", breaks = ms) + 
  scale_y_continuous(expression(hat(Ne))) + ggtitle(paste0("stepping stone with  ", n_subpops, " populations")) + geom_hline(yintercept = 1000)#+ ylim(c(900,10100))


# A particular estimation method
means_full = bind_rows(means_list, .id = "demes")
means_full$demes = as.factor(as.character(as.numeric((substr(means_full$demes, 3,20)))**2))
means_full$demes = factor(means_full$demes, levels = edge_lengths**2)

max_m = 1
method = "AF9"
means_examine_plot <- trim_df(means_full, max_m = max_m, methods = c(method))


ggplot(data=means_examine_plot,
       aes(x=m, y=est, colour=demes)) +
  scale_color_manual( values=colfunc_AF(length(edge_lengths))) +
  geom_line() + theme_bw() + scale_x_continuous(trans = "log10", breaks = ms) + 
  scale_y_continuous(expression(hat(Ne))) + ggtitle(paste0("stepping stone with  ", n_subpops, " populations")) + geom_hline(yintercept = 1000)#+ ylim(c(900,10100))


#######
# differences
#######

## successive

methods <- get_methods(t_intervals, cs)

v = rep(NA, length(ms)*length(methods))
diff_list <- list()

for (i in 2:length(edge_lengths)){
row = 1
means_diff <- data.frame(m = v, method = v, est = v)
  
first_df  = means_list[[i-1]]
second_df = means_list[[i]]
for (m_i in ms){
  for (method_i in methods){
    first_est <- first_df[ which(first_df$m==m_i & first_df$method == method_i), 'est']
    second_est <- second_df[ which(second_df$m==m_i & second_df$method == method_i), 'est']
    means_diff[row,] <- c(m_i, method_i, (second_est - first_est) / second_est)
    row = row+1
  }
}
means_diff$m = as.numeric(means_diff$m)
means_diff$est = as.numeric(means_diff$est)
diff_list[[i-1]] = means_diff
}


library(gridExtra)
plot_list <- list()
i=1
max_m=1
for (df in diff_list){
  df <- trim_df(df, max_m = max_m)#, methods = c("AF2"))
  x=edge_lengths[[i+1]]**2
  y=edge_lengths[[i]]**2
  
  title_str <- paste0("(Est with ", x, " subpops - with ", y, " subpops) / ", x , " subpops")
  
  plot_list[[i]] <- ggplot(data=df,
         aes(x=m, y=est, colour=method, linetype = method)) +
    scale_color_manual( values=c(colfunc_AF(length(t_intervals)-1), rep("black",length(cs))) ) +
    scale_linetype_manual( values=c(rep(1, length(t_intervals)-1), ltys[1:(length(cs))] ) ) +
    geom_line() + theme_bw() + scale_x_continuous(trans = "log10", breaks = ms) + 
    scale_y_continuous(expression(delta)) + ggtitle(title_str) + geom_hline(yintercept = 0) + ylim(c(-0.5, 0.5))
  i = i+1
}
plot_list[['ncol']] <- 1
do.call(grid.arrange, plot_list)




## to largest simulation


methods <- get_methods(t_intervals, cs)

v = rep(NA, length(ms)*length(methods))
diff_list <- list()

for (i in 2:length(edge_lengths)){
  row = 1
  means_diff <- data.frame(m = v, method = v, est = v)
  
  first_df  = means_list[[i-1]]
  second_df = means_list[[length(edge_lengths)]]
  for (m_i in ms){
    for (method_i in methods){
      first_est <- first_df[ which(first_df$m==m_i & first_df$method == method_i), 'est']
      second_est <- second_df[ which(second_df$m==m_i & second_df$method == method_i), 'est']
      means_diff[row,] <- c(m_i, method_i, (second_est - first_est) / second_est)
      row = row+1
    }
  }
  means_diff$m = as.numeric(means_diff$m)
  means_diff$est = as.numeric(means_diff$est)
  diff_list[[i-1]] = means_diff
}


library(gridExtra)
plot_list <- list()
i=1
max_m=1
for (df in diff_list){
  df <- trim_df(df, max_m = max_m)#, methods = c("AF2"))
  x=edge_lengths[[length(edge_lengths)]]**2
  y=edge_lengths[[i]]**2
  
  title_str <- paste0("(Est with ", x, " subpops - with ", y, " subpops) / ", x , " subpops")
  
  plot_list[[i]] <- ggplot(data=df,
                           aes(x=m, y=est, colour=method, linetype = method)) +
    scale_color_manual( values=c(colfunc_AF(length(t_intervals)-1), rep("black",length(cs))) ) +
    scale_linetype_manual( values=c(rep(1, length(t_intervals)-1), ltys[1:(length(cs))] ) ) +
    geom_line() + theme_bw() + scale_x_continuous(trans = "log10", breaks = ms) + 
    scale_y_continuous(expression(delta)) + ggtitle(title_str) + geom_hline(yintercept = 0) + ylim(c(-0.5, 0.5))
  i = i+1
}
plot_list[['ncol']] <- 1
do.call(grid.arrange, plot_list)





max_m <- .3
means_diff_plot <- trim_df(means_diff, max_m = max_m)

ggplot(data=means_diff_plot,
       aes(x=m, y=est, colour=method, linetype = method)) +
  scale_color_manual( values=c(colfunc_AF(length(t_intervals)-1), rep("black",length(cs))) ) +
  scale_linetype_manual( values=c(rep(1, length(t_intervals)-1), ltys[1:(length(cs))] ) ) +
  geom_line() + theme_bw() + scale_x_continuous(trans = "log10", breaks = ms) + 
  scale_y_continuous(expression(delta)) + ggtitle("Difference between 100 subpop vs one large subpop estimates") + ylim(c(-0.1, 0.1))



########################################################################
########################################################################





N = 1000
S = Inf
ms = c(0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.9, 1)
M_fn = construct_M_inf
p0 = rep(.25, 4)
N_large  = 1e9
N_v = c(N, N_large)
#S = Inf
ms = c(0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.9)
M_fn = construct_M_inf
m = .1
M <- M_fn(m)



burn_in_checker(ms=c(.1), M_fn, cs=.5, burn_in=500, N_v=N_v, loci, p0)


o <- sim.drift5_mig(p0 = p0, N_v = N_v, M = M, c = .5, t = 100, loci = 1e2, burn_in = 0)

plot(1:100, get_Fst_vec(o$pA_post, N_v, 100)); abline(h=(1-2*m)/(4*N*m + 1 - 2*m))

hist((o$pA_post[,100,1]-0.5)**2, breaks = 100); abline(v = (1-2*m)/(4*N*m+1-2*m), col = 'red'); abline(v =2* (1-2*m)/(4*N*m+1-2*m), col = 'blue')




methods_to_plot <- c("LD", "AF2", "AF5", "AF9")

means_sym_plot <- trim_df(means_sym, methods_to_plot, max_m)

g <- ggplot(data=means_inf_plot,
      aes(x=m, y=est, colour=method)) + scale_color_manual(values=c(gg_color_hue(length(plot_lines)-1), "black") )+
      geom_line() + theme_bw() + scale_x_continuous(trans = "log10", breaks = ms) #+ scale_y_continuous(expression(hat(Ne)), breaks =seq(1000,10000,1000), limits=c(900,10100)) #+ ylim(c(900,10100))

g

# and add the symetrical migration plot
g <- g + geom_line(data=means_sym_plot,
                aes(x=m, y=est), linetype = 2)# + scale_color_manual(values=c(gg_color_hue(length(plot_lines)-1), "black"))
g



plot(1:5, 1:5, lty = 7, type = 'l')
