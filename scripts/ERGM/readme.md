# Estimation Strategy: Multilayer ERGMs

## Description

### 0

Data check files was used to check for necessary imputation in data obtained from replication files and own pre-processing. 


### 1 

Cross-sectional Multilayer ERGM estimation. First, models are calibrated with Stochastic Approximation. Final results serve as starting values for the MCMLE approach in the ERGM package. This 'two-stage' approach results in the best trace plots for the sampled statistics. 


### 2

For documentation only. Temporal extension of the Multilayer ERGM is currently not possible because the ``btergm`` package does not support appropriate constraints on a multiplex network as required by the implementation of Chen 2021. 






