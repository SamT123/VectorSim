###################################
#
# Investigating burn-in times
#
###################################

source('simulator_functions.R')
source('abc_functions.R')



ns = 81

N <- 1000
m = 0.05
print(N)
print(m)
pop <- initialise_population(n_subpops = ns, loci = 100, t=1)

M <- construct_M_4_step(m, ns)
N_v= rep(N, ns)

pop <- propegate_population(pop, N_v, M, c=0.05, t_add = 50)



N <- 100000
m = 0.01
M <- construct_M_4_step(m, ns)
N_v= rep(N, ns)


pop <- propegate_population(pop, N_v, M, c=0.05, t_add = 50)


fst_v = colMeans(rbind(get_Fst_vec(pop$pA_post, N_v, 50), get_Fst_vec(pop$pB_post, N_v, 50)))
plot(1:50, fst_v, type='l') ; abline(v=500); abline(v=640)
plot(40:1038, ma(diff(fst_v),n=40), type='l');     abline(h=0); abline(v=500); abline(v=640)


r2_v = apply(pop$r_post, MARGIN = c(3), mean)
plot(1:50, r2_v[1:50], type = 'l')


###

ns = 100

pop <- initialise_population(n_subpops = ns, loci = 100, t = 501)

M <- construct_M_4_step(0.08, ns)
pop <- propegate_population(pop, N_v, M, t_add = 500)

fst_v = get_Fst_vec(pop$pA_post, N_v, 500)
plot(1:500, fst_v)

r2_v = apply(pop$r_post, MARGIN = c(2), mean)
plot(1:700, r2_v[1:700])


