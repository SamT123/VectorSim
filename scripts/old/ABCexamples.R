th



# Approximate Bayes
library(abc)

# estimating the parameters of a normal distribution


## make real data set
mu_true  = 3.0
sig_true = 2.5
sample_size = 10000
d <- sort(rnorm(sample_size, mu_true, sig_true))


## prior distributions for parameters
mu_prior <- function(n=1){
  rnorm(n, 0,10)
}

sig_prior <- function(n=1){
  runif(n, 0, 10)
}



## simulate
n_sims = 1000

mu_sim = mu_prior(n_sims)
sig_sim = sig_prior(n_sims)


d_sim = array(NA, dim = c(n_sims, sample_size))

## choose summary stats

### euclidean distance
dist_out = rep(NA, n_sims)
dist_real_sample = c(0)

### sample mean and varaince
mu_out = rep(NA, n_sims)
sig_out = rep(NA, n_sims)
mu_real_sample = mean(d)
sig_real_sample = sqrt(var(d))

### sample quantiles
q25_out = rep(NA, n_sims)
q75_out = rep(NA, n_sims)
qsep_out = rep(NA, n_sims)
q25_real_sample = quantile(d)[2]
q75_real_sample = quantile(d)[4]
qsep_real_sample = q75_real_sample - q25_real_sample

for (i in 1:n_sims){
  o = sort(rnorm(sample_size, mu_sim[i], sig_sim[i]))
  d_sim[i, ] = o
  
  ### calculate SS
  dist_out[i] = dist(rbind(d, o))
  mu_out[i] = mean(o)
  sig_out[i] = sqrt(var(o))
  q25_out[i] = quantile(o)[2]
  q75_out[i] = quantile(o)[4]
  qsep_out[i] = q75_out[i] - q25_out[i]
}


## check that siummary stats correlate with parameters of interest
par(mfrow = c(2,1))
plot(mu_sim, dist_out); plot(sig_sim, dist_out)

plot(mu_sim, mu_out); plot(mu_sim, sig_out)
plot(sig_sim, mu_out); plot(sig_sim, sig_out)

plot(mu_sim, mu_out); plot(mu_sim, qsep_out)
plot(sig_sim, mu_out); plot(sig_sim, qsep_out)




## when doing the inference, botht he choise of summary stat and the method matter:
# rejection sampling is not great with smallish sample sizes
o <- abc(target = c(dist_real_sample), param = cbind(mu_sim, sig_sim), sumstat = cbind(dist_out), tol = .05, method = 'rejection')
summary(o)
hist(o)

# adding some adjustment really helps
o <- abc(target = c(mu_real_sample, qsep_real_sample), param = cbind(mu_sim, sig_sim), sumstat = cbind(mu_out, qsep_out), tol = .05, method = 'loclinear')
summary(o)
hist(o)

# rejection sampling is not great with smallish sample sizes
o <- abc(target = c(mu_real_sample, qsep_real_sample), param = cbind(mu_sim, sig_sim), sumstat = cbind(mu_out, qsep_out), tol = .05, method = 'rejection')
summary(o)
hist(o)

# adding some adjustment really helps
o <- abc(target = c(mu_real_sample, qsep_real_sample), param = cbind(mu_sim, sig_sim), sumstat = cbind(mu_out, qsep_out), tol = .05, method = 'loclinear')
summary(o)
hist(o)




# SMC ABC


sequential_start_time <- proc.time()
## make real data set
mu_true  = 3.0
sig_true = 2.5
sample_size = 10000
dat <- sort(rnorm(sample_size, mu_true, sig_true))

# Number of SMC rounds
rounds <- 5
# Number of iterations per round
N <- 10000
# Fraction of samples to save per round
tau <- c(0.2,rep(0.25,rounds-1))
# Number of samples to save from first round
Q <- N*tau[1]
# Threshold value for simulated vectors to be saved
epsilon <- rep(NA,rounds)
# Storage vector for dist measure from simulated data
dist_vals <- rep(NA,N)
# Draw initial mean and variance samples for first round from prior
tmp_mean <- runif(N, -20, 20)
tmp_variance <- runif(N, 0, 50)


## -------- 1st round -------- ##
# Draw N data sets and compare
for(i in 1:N){
  tmp_dat <- rnorm(length(dat), tmp_mean[i], sqrt(tmp_variance[i]))
  tmp_dat <- sort(tmp_dat)
  dist_vals[i] <- dist(rbind(dat,tmp_dat))
}
# Sort dist values and save the top Q means and variances
dist_indexes <- sort(dist_vals, index.return=T)
save_indexes <- dist_indexes$ix[1:Q]
epsilon[1] <- dist_indexes$x[Q] # first epsilon is the max dist value
saved_means <- tmp_mean[save_indexes]
saved_variances <- tmp_variance[save_indexes]


## -------- 2nd through rth rounds -------- ##
for(r in 2:rounds){
  
  curr_num_saved <- 0
  dist_vals <- rep(NA,Q)
  tmp_saved_means <- rep(NA,Q)
  tmp_saved_variances <- rep(NA,Q)
  while(curr_num_saved < Q){
    curr_mean <- sample(saved_means,1)
    curr_mean <- curr_mean + runif(1,-0.05,0.05)
    curr_variance <- sample(saved_variances,1)
    curr_variance <- curr_variance + runif(1,0,0.05)
    curr_dat <- rnorm(30,curr_mean,sqrt(curr_variance))
    curr_dat <- sort(curr_dat)
    curr_dist <- dist(rbind(dat,curr_dat))
    
    # Only save value if it is better than the previously saved max
    if(curr_dist < epsilon[r-1]){
      
      curr_num_saved <- curr_num_saved + 1
      dist_vals[curr_num_saved] <- curr_dist
      tmp_saved_means[curr_num_saved] <- curr_mean
      tmp_saved_variances[curr_num_saved] <- curr_variance
      
    }
    
  }
  
  # Save mean and variance samples
  dist_indexes <- sort(dist_vals, index.return=T)
  save_indexes <- dist_indexes$ix[1:(tau[r]*Q)]
  epsilon[r] <- dist_indexes$x[(tau[r]*Q)]
  saved_means <- tmp_saved_means[save_indexes]
  saved_variances <- tmp_saved_variances[save_indexes]
}


sequential_total_time <- proc.time() - sequential_start_time
# Make a data frame and print as a knitr table
res <- data.frame(Parameter=c("Normal Mean","Normal Variance"),
                  Mean=c(mean(saved_means), mean(saved_variances)),
                  SD=c(sd(saved_means),sd(saved_variances)))
knitr::kable(res, digits=3)
# Plot the posterior distribution of the parameters
par(mfrow=c(1,2))
hist(saved_means, main="Mean", xlab=expression(hat(mu)), ylab="")
abline(v=mean(saved_means),col="blue",lty="dashed",lwd=2)
hist(saved_variances, main="Variance", xlab=expression(hat(sigma)^2), ylab="")
abline(v=mean(saved_variances),col="blue",lty="dashed",lwd=2)
```

```{r}
timing <- data.frame(Algorithm=c("ABC","ABC-SMC"),
                     Time=c(normal_total_time[3],sequential_total_time[3]),
                     Improvement=c(normal_total_time[3]/normal_total_time[3],
                                   normal_total_time[3]/sequential_total_time[3]))
knitr::kable(timing, digits=3)
```
