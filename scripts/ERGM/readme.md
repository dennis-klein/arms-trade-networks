# Estimation Strategy: Multilayer ERGMs

## Description

### 0

Data check files was used to check for necessary imputation in data obtained from replication files and own pre-processing. 


### 1 

Cross-sectional Multilayer ERGM estimation.We follow the specification in 

Chen, Ted Hsuan Yun. 2021. “Statistical Inference for Multilayer Networks in Political Science.” Political Science Research and Methods 9 (2). Cambridge University Press: 380–97. doi:10.1017/psrm.2019.49.

First, models are calibrated with Stochastic Approximation. Final results serve as starting values for the MCMLE approach in the ERGM package. This 'two-stage' approach results in the best trace plots for the sampled statistics. 

Additional out-of-sample (oos) validation is provided.

### 2

For documentation only and not pursued or validated further. Temporal extension of the Multilayer ERGM is currently not possible because the ``btergm`` package does not support appropriate constraints on a multiplex network as required by the implementation of Chen 2021. 






