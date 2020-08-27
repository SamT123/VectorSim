source("simulator_functions.R")
library(dplyr)


############################
#
# SINGLE LINE
#
############################

####
#
# fully connected subpopulations
#
###


# population structure
N = 1000
S = Inf
n_subpops = 5
ms = c(0, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.9, 1)
M_fn = construct_M_waples

p0 = rep(.25, 4)

# estimation
t_intervals = c(0,20)
cs = c(0.5)

N_v = rep(N, n_subpops)

# simulation parameters
loci = 1e2
burn_in = 150

reset_period = 10
sampling_repeats = 20
sampling_period = reset_period + max(t_intervals)

o <- burn_in_checker(ms, M_fn, cs, burn_in, N_v, loci, p0)


estimates_sym_single <- get_estimates_df(N_v = N_v,
                                  ms = ms,
                                  M_fn = M_fn,
                                  p0 = p0,
                                  loci = loci,
                                  burn_in = burn_in,
                                  t_intervals = t_intervals,
                                  rest_period = rest_period,
                                  sampling_repeats = sampling_repeats,
                                  sample_subpops = 'all',
                                  S = S,
                                  cs = cs)



means_sym_single <- get_means_df(estimates_sym_single)
max_m = 1
means_sym_single_plot <- trim_df(means_sym_single, max_m = max_m, methods = c("AF20"))
means_sym_single_plot[means_sym_single_plot$m == 0,]$m = 0.001
means_sym_single_plot$est = means_sym_single_plot$est / N

label_locs = ms
label_locs[1] = 0.001

colfunc_AF <- colorRampPalette(c("red", "orange"))
ltys <- c(2,4,5,6)

# plot

png("../results/TY_single_connected.png")
ggplot(data=means_sym_single_plot, aes(x=m, y=est)) +
  geom_line() + theme_bw() + 
  scale_y_continuous(expression(hat(N) / N), breaks =seq(0, N*n_subpops)) + 
  scale_x_continuous("Migration rate", breaks = label_locs, labels = ms, trans = 'log10') +
  ggtitle(paste0("Fully connected migration amongst ", n_subpops, " populations")) #+ ylim(c(900,10100))
dev.off()

ggplot(data=means_sym_single_plot, aes(x=m, y=est)) +
  geom_line() + theme_bw() + 
  scale_y_continuous(expression(hat(N) / N), breaks =seq(0, N*n_subpops)) + 
  scale_x_continuous("Migration rate", breaks = label_locs, labels = ms) +
  ggtitle(paste0("Fully connected migration amongst ", n_subpops, " populations")) #+ ylim(c(900,10100))
  
###
#
# stepping stone migration
#
###
  


# population structure
N = 1000
S = Inf
n_subpops = 100
ms = c(0, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.9, 1)
M_fn = construct_M_4_step

p0 = rep(.25, 4)

# estimation
t_intervals = c(0,20)
cs = c(0.5)

N_v = rep(N, n_subpops)

# simulation parameters
loci = 1e2
burn_in = 150

reset_period = 10
sampling_repeats = 20
sampling_period = reset_period + max(t_intervals)

o <- burn_in_checker(ms, M_fn, cs, burn_in, N_v, loci, p0)


estimates_step_single <- get_estimates_df(N_v = N_v,
                                  ms = ms,
                                  M_fn = M_fn,
                                  p0 = p0,
                                  loci = loci,
                                  burn_in = burn_in,
                                  t_intervals = t_intervals,
                                  rest_period = rest_period,
                                  sampling_repeats = sampling_repeats,
                                  sample_subpops = 50,
                                  S = S,
                                  cs = cs)


means_step_single <- get_means_df(estimates_step_single)
max_m = 1
means_step_single_plot <- trim_df(means_step_single, max_m = max_m, methods = c("AF20"))
means_step_single_plot[means_step_single_plot$m == 0,]$m = 0.001
means_step_single_plot$est = means_step_single_plot$est / N

label_locs = ms
label_locs[1] = 0.001



# plot
png("../results/TY_single_step.png")

ggplot(data=means_step_single_plot, aes(x=m, y=est)) +
  geom_line() + theme_bw() + 
  scale_y_continuous(expression(hat(N) / N), breaks =seq(0, N*n_subpops)) + 
  scale_x_continuous(breaks = label_locs, labels = ms, trans = 'log10') +
  ggtitle(paste0("Stepping stone model with ", n_subpops, " populations")) #+ ylim(c(900,10100))
dev.off()

savefig("../results/TY_single_step.png")

############################
#
# VARYING N_SUBPOPS
#
############################

###
#
# fully connected
#
###

# population structure
N = 1000
S = Inf
n_subpops_sym = c(5,7,9,11,15)
ms = c(0, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.9, 1)
M_fn = construct_M_waples

p0 = rep(.25, 4)

# estimation
t_intervals = c(0,20)
cs = c(0.5)

means_sym_list = list()

for (n_subpops in n_subpops_sym){

print(n_subpops)
  
N_v = rep(N, n_subpops)


# simulation parameters
loci = 1e2
burn_in = 150

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
                                  sample_subpops = 'all',
                                  S = S,
                                  cs = cs)



means_sym <- get_means_df(estimates_sym)

means_sym_list[[paste0("d", n_subpops)]] = means_sym
}



means_sym_full = bind_rows(means_sym_list, .id = "demes")
means_sym_full <- trim_df(means_sym_full, max_m = 1, methods = c("AF20"))
means_sym_full[means_sym_full$m == 0,]$m = 0.001
means_sym_full$est = means_sym_full$est / N
means_sym_full$demes = as.factor((substr(means_sym_full$demes, 2,20)))
means_sym_full$demes = factor(means_sym_full$demes, levels = n_subpops_list)

colfunc_AF <- colorRampPalette(c("yellow3", "#f03b20"))
ltys <- c(2,4,5,6)


png("../results/TY_multi_connected.png")

ggplot(data=means_sym_full, aes(x=m, y=est, colour=demes)) +
  scale_color_manual( values=colfunc_AF(length(n_subpops_sym))) +
  geom_line() + theme_bw() + 
  scale_y_continuous(expression(hat(N) / N), breaks =seq(0, N*n_subpops)) + 
  scale_x_continuous(breaks = label_locs, labels = ms, trans = 'log10') +
 ggtitle(paste0("Fully connected migration")) #+ ylim(c(900,10100))
 
dev.off()

png("../results/TY_multi_connected_3d.png")

wireframe(est~m*demes,data=means_sym_full,
          screen=list(  z = 30 , x=-70),
          xlab="Migration rate",
          ylab="Number of demes",
          zlab = expression(hat(N) / N),
          main="Fully connected migration" )

dev.off()

means_sym_log = means_sym_full
means_sym_log$m[means_sym_log$m == 0] = 0.001
means_sym_log$m = log(means_sym_log$m)

png("../results/TY_multi_connected_3d_log.png")

wireframe(est~m*demes,data=means_sym_log,
          screen=list(  z = 30 , x=-70),
          xlab="log Migration rate",
          ylab="Number of demes",
          zlab = expression(hat(N) / N))

dev.off()

###
#
# stepping stone migration
#
###

# population structure
N = 1000
S = Inf
n_subpops_step = c(9,25,49, 81)
ms = c(0, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.9, 1)
M_fn = construct_M_4_step

p0 = rep(.25, 4)

# estimation
t_intervals = c(0,20)
cs = c(0.5)

means_step_list = list()

for (n_subpops in n_subpops_step){
  
  print(n_subpops)
  
  N_v = rep(N, n_subpops)
  
  
  # simulation parameters
  loci = 1e2
  burn_in = 150
  
  reset_period = 10
  sampling_repeats = 20
  sampling_period = reset_period + max(t_intervals)
  
  o <- burn_in_checker(ms, M_fn, cs, burn_in, N_v, loci, p0)
  
  
  estimates_step <- get_estimates_df(N_v = N_v,
                                    ms = ms,
                                    M_fn = M_fn,
                                    p0 = p0,
                                    loci = loci,
                                    burn_in = burn_in,
                                    t_intervals = t_intervals,
                                    rest_period = rest_period,
                                    sampling_repeats = sampling_repeats,
                                    sample_subpops = c(n_subpops/2+.5),
                                    S = S,
                                    cs = cs)
  
  
  
  means_step <- get_means_df(estimates_step)
  
  means_step_list[[paste0("d", n_subpops)]] = means_step
}






means_step_full = bind_rows(means_step_list, .id = "demes")

means_step_full <- trim_df(means_step_full, max_m = 1, methods = c("AF20"))
means_step_full[means_step_full$m == 0,]$m = 0.001
means_step_full$est = means_step_full$est / N
means_step_full$demes = as.factor((substr(means_step_full$demes, 2,20)))
means_step_full$demes = factor(means_step_full$demes, levels = n_subpops_list)

png("../results/TY_multi_step.png")

ggplot(data=means_step_full, aes(x=m, y=est, colour=demes)) +
  scale_color_manual( values=colfunc_AF(length(n_subpops_step))) +
  geom_line() + theme_bw() + 
  scale_y_continuous(expression(hat(N) / N), breaks =seq(0, N*n_subpops)) + 
  scale_x_continuous(breaks = label_locs, labels = ms, trans = 'log10') +
  ggtitle(paste0("Stepping stone migration")) #+ ylim(c(900,10100))
dev.off()




png("../results/TY_multi_step_3d.png")

wireframe(est~m*demes,data=means_step_full,
                      screen=list(  z = 30 , x=-70),
                      xlab="Migration rate",
                      ylab="Number of demes",
                      zlab = expression(hat(N) / N) ,
          main="Stepping stone migration")

dev.off()



means_step_log = means_step_full
means_step_log$m[means_step_log$m == 0] = 0.001
means_step_log$m = log(means_step_log$m)

png("../results/TY_multi_step_3d_log.png")


wireframe(est~m*demes,data=means_step_log,
          screen=list(  z = 30 , x=-70),
          xlab="log Migration rate",
          ylab="Number of demes",
          zlab = expression(hat(N) / N) )

dev.off()
#########################################################################################################

i = 1
j = 1
arr <- array(NA, dim = c(length(ms), length(n_subpops_list) ) )
for ( m in label_locs ){
  j=1
  for (d in n_subpops_list){
    arr[[i, j]] = means_sym_full[means_sym_full$m == m & means_sym_full$demes == d, ]$est
    j = j+1
  }
  i = i+1
  
}

persp(ms, n_subpops_list, arr, theta = -40, phi = 20, zlim = c(min(arr), max(arr)), main = "Stepping stone migration", xlab = "Migration rate", ylab = "Number of demes", N)



library(plotly)
# df_sample could be your plot_data    
df<-apply((as.matrix(means_sym_full[,-3], rownames.force = NA)), 2, as.numeric)
plot_ly(z=df[,c(1,3,2)]) %>% add_surface()


plot_ly(x=~df[,'m'], y=~df[,'demes'], z=~df[,'est'], type="mesh3d", mode="markers")





volcano






means_sym_relative = means_sym_full
for (r in 1:dim(means_sym_relative)[1]){
  means_sym_relative[r,4] = means_sym_relative[r,4]/means_sym_relative[r,1]
}


ggplot(data=means_sym_relative, aes(x=m, y=est, colour=demes)) +
  scale_color_manual( values=colfunc_AF(length(n_subpops_list))) +
  geom_line() + theme_bw() + 
  scale_y_continuous(expression(hat(N) / N), breaks =seq(0, N*n_subpops)) + 
  scale_x_continuous(breaks = label_locs, labels = ms, trans = 'log10') +
  ggtitle(paste0("Symmetric migration amongst ", n_subpops, " populations")) #+ ylim(c(900,10100))
