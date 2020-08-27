library(Rcpp)
library(RcppArmadillo)
sourceCpp(file = "propegate_new.cpp")

l = list(1,2,3,4,5)
l = test_fn_list(l)