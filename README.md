# Functional Data Analysis for Longitudinal Data with Informative Observation Times

The R scripts in this repository contain the code necessary to reproduce the simulations of the manuscript "Functional Data Analysis for Longitudinal Data with Informative Observation Times" (C. Weaver, L. Xiao, and W. Lu). 

* src.R contains the code necessary to run mean and covariance estimation as described in the manuscript

* simdat.R contains the code necessary to generate data according to the simulation study design

An example of running this code is as follows:

```{r}
param = list(n=150, rho=0.5, dist='gamma', "C"=1, sigz=1, sige=1)
data = datf(param)[[1]]
fit.o = face.sparse.weighted(data, weight='OBS') # OBS weighting scheme
fit.s = face.sparse.weighted(data, weight='SUBJ') # SUBJ weighting scheme
```
