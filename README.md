# VectorSim

# **MSc Thesis Repository**
## About
Author: Sam Turner

This repository contains code used for my Computational Methods in Ecology and Evolution MSc thesis.

## scripts

**ABC_experiments.R**
- All code relating to estimation of $N_e$ and $m$ by Approximate Bayesian Computation.

**bias_demostration**
- Demonstration that the use of classical linkage disequilibrium and allelic fluctuation estimators results in biased estimates in the stepping stone model

**classic_estimators.R**
- Demonstration of classical linkage disequilibrium and allelic fluctuation estimators on single subpopulations

**poolseq.R**
- Finding the relationship between individual sequencing allele frequency estimates and pooled sequencing estimates.

**propegate_new.cpp**
- The internal propagation script for the allele frequency estimator implemented in `simulator_functions.R`

**simulator_functions.R**
- All major functions are contained in this file, including the haplotype frequency simulator, sampling protocol, full ABC estimation method, and helper functions for evaluating the performance of the method. 


