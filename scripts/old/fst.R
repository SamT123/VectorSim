source("simulator_functions.R")
# set equilibrium value with 500 generation average
# so set t large enough so we have 500 generations clearly at equilibrium at the end
# for each migration rate: set eqm = average Fst in last 500 gens
# time to eqm = first generation that reaches this average value
# can come up with nicer metric if needed

# population parameters

N = 1000
S = 500
n_subpops = 100
ms = c(0.01, 0.03, 0.05, 0.1, 0.3)
p0<-rep(0.25,4)
loci = 1e2

N_v = rep(N, n_subpops)


t = 500

# check that t is large enough for minimum m value
M = construct_M_4_step(m=min(ms), n_subpops = n_subpops)
o <- sim.drift5_mig(N_v = N_v, t = t, M=M, loci = loci, p0=p0)
fst_v <- get_Fst_vec(o$pA_pre, N_v, t)
plot(1:t, ma(fst_v, n =1), type = "l"); abline(v = t-500, lty = 2, col = 'red')

# calculate and plot t_eqm times
os_post = list()
os_pre  = list()
for (m in ms){
  print(m)
  M = construct_M(m=m, n_subpops = n_subpops)
  o <- sim.drift5_mig(N_v = N_v, t = t, M=M, loci = loci, p0=p0)
  os_post[[as.character(m)]] <- o$pA_post
  os_pre[[as.character(m)]] <- o$pA_pre
  
}

# raw plots
par(mfrow = c(6,2), mai=rep(.3,4))
ma_n = 10

for (pA in os_pre){
  fst_v <- get_Fst_vec(pA, N_v, t)
  plot(1:t, ma(fst_v, n = ma_n), type = "l", xlim = c(0, t))
  plot(1:(t-1), diff(ma(fst_v, ma_n)), type = "l", xlim = c(0, t))
}

for (pA in os_post){
  fst_v <- get_Fst_vec(pA, N_v, t)
  plot(1:t, ma(fst_v, n = ma_n), type = "l", xlim = c(0, t))
  plot(1:(t-1), diff(ma(fst_v, ma_n)), type = "l", xlim = c(0, t))
}


# estimate eqm
t_eqm_v <- c()
par(mfrow = c(6,2), mai=rep(.3,4))
for (pA in os_post){
  fst_v <- get_Fst_vec(pA, N_v, t)
  target_1 <-  mean(fst_v[(t-1000):t])
  target_2 <-  min(fst_v[(t-1000):t])
  target_3 <-  0
  t_eqm_1 <- which.max(ma(fst_v) > target_1*.97)
  t_eqm_2 <- which.max(ma(fst_v) > target_2)
  t_eqm_3 <- which.max(ma(diff(fst_v)) < target_3)
  
  plot(1:t, ma(fst_v, n = 1), type = "l", xlim = c(0, min(t, 3*t_eqm_1))); abline(v = c(t_eqm_1, t_eqm_2, t_eqm_3), col = c('red', 'green', 'blue'))
  plot(1:(t-1), diff(ma(fst_v, 3)), type = "l", xlim = c(0, min(t, 3*t_eqm_1)))
  
  t_eqm_v <- c(t_eqm_v, t_eqm_2)
}
  


# fit model t_eqm = a / m
df <- data.frame(t_eqm = t_eqm_v, m_inv = ms**-1)
mod <- lm(t_eqm ~ m_inv, data = df)
m_vals = seq(0.001,0.3,0.001)
predicted_t_eqm <- predict(mod, list(m_inv = m_vals**-1))

# 
par(mfrow = c(1,1), mai = rep(1,4))
plot(ms, t_eqm_v, pch = 1, cex = 2)
lines(m_vals, predicted_t_eqm, lty = 2, col = 'red')





