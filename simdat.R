# Simulate data
# param: a list of setting values:
#  n: number of subjects
#  rho: frailty parameter (0-1)
#  dist: distribution of frailty variable ('lognormal' or 'gamma')
#  sigz: variance of fraily variable
#  sige: variance of white noise
datf = function(param){
  grid = seq(0.025, 0.975, 0.025)
  n = param[['n']]
  rho = param[['rho']]
  dist = param[['dist']]
  C = param[['C']]
  sigz = param[['sigz']]
  sige = param[['sige']]
  if(sigz==1){
    shape=1;scale=1;
    logmn = log(1^2 / sqrt(1 + 1^2))
    logsd = sqrt(log(1 + (1/ 1^2)))
  }
  if(sigz==.5){
    shape=2;scale = .5;
    logmn <- log(1^2 / sqrt(.5 + 1^2))
    logsd <- sqrt(log(1 + (.5 / 1^2)))
  }
  if(C==0){
    C = runif(n)
  }   else if(C==1){
    C = rep(1, n)
  }    else if(C==2){
    u1 = runif(n)
    u2 = runif(n,.2,1)
    C = ifelse(u1 < .2, 1, u2)
  } 
  C = sort(C)
  if(sigz==0){
    Zstar = rep(1,n)
  }  else {
    if(dist=="gamma") Zstar = rgamma(n,shape, scale=scale)
    if(dist=="lognormal") Zstar = rlnorm(n,logmn,logsd)
  }
  Z = (1-rho) + rho*Zstar
  m = rpois(n,10*Z*C)
  m[m<2] = 2
  t = lapply(1:n, function(i) sort(runif(m[i], 0, C[i])))
  truemuf = function(t){ return(t^2 + 2/3)}
  truemus = lapply(t, truemuf)
  truemugrid = truemuf(grid)
  muf = function(ti, Zi) Zi * (ti^2 + 2/3)
  sigmuf = function(ts){
    mat = matrix(NA,length(ts),length(ts))
    for(i in 1:length(ts)){
      for(j in 1:length(ts)){
        s = ts[i]
        t = ts[j]
        mat[i,j] = sigz * rho^2 * (s^2 +2/3) * (t^2 + 2/3)
      }
    }
    return(mat)
  }
  sigxf = function(ts){
    mat = matrix(NA,length(ts),length(ts))
    for(i in 1:length(ts)){
      for(j in 1:length(ts)){
        s = ts[i]
        t = ts[j]
        mat[i,j] = (1 - abs(s-t)) * (.05)*exp(s^2) * (.05)*exp(t^2)
      }
    }
    return(mat)
  }
  sigx = lapply(1:n,function(i) sigxf(t[[i]]))
  truesigxgrid = sigxf(grid)
  truesigmugrid = sigmuf(grid)
  truecovgrid = truesigxgrid + truesigmugrid
  mu = lapply(1:n,function(i) muf(t[[i]],Z[i]))
  X = lapply(1:n, function(i) rmvnorm(n=1, sigma=sigx[[i]]))
  eps = lapply(1:n, function(i) rnorm(m[i], 0, sige))
  Y = lapply(1:n, function(i){
    mu[[i]] + X[[i]] + eps[[i]]
  })
  ts = unlist(t)
  ord = order(ts)
  mus = unlist(mu)
  ys = unlist(Y)
  Zs = Cs = id = c()
  for(i in 1:length(t)){id = c(id, rep(i, length(t[[i]])))}
  for(i in 1:length(t)){Cs = c(Cs, rep(C[i], length(t[[i]])))}
  for(i in 1:length(t)){Zs = c(Zs, rep(Z[i], length(t[[i]])))}
  data=data.frame(argvals=ts, subj=id, y=ys, C=Cs, Z=Zs)
  return(list(data=data, truemugrid=truemugrid, truecovgrid=truecovgrid))
}