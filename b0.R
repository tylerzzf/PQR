
## heteroskedasticity model
rm(list=ls())
setwd("F:\\ARCH研究\\code\\Mult-pqr\\新方法")
source("Mult-quantile.R")
source("nobreaksimfunc.R")
setwd("F:\\ARCH研究\\code\\Mult-pqr\\新方法\\hete error")
source("DGP.R")

nrep = 250
p=2
taus=c(0.25,0.5,0.75)
cpt.num = 0
cpt.eps = 1e-6
eps = 1e-6
step.max = 50
set.seed(12345)
seeds = sample(1:999999, nrep, replace=FALSE)

n=T=40
tem.44 = nobreak.sim.func(n,T,p,taus,cpt.num,nrep,eps,cpt.eps,step.max)
num.acc.144 = tem.44$num.acc.1/nrep
num.acc.244 = tem.44$num.acc.2/nrep
num.acc.344 = tem.44$num.acc.3/nrep
time.44 = tem.44$time

n=40 
T=80
tem.48 = nobreak.sim.func(n,T,p,taus,cpt.num,nrep,eps,cpt.eps,step.max)
num.acc.148 = tem.48$num.acc.1/nrep
num.acc.248 = tem.48$num.acc.2/nrep
num.acc.348 = tem.48$num.acc.3/nrep
time.48 = tem.48$time

n=80
T=40
tem.84 = nobreak.sim.func(n,T,p,taus,cpt.num,nrep,eps,cpt.eps,step.max)
num.acc.184 = tem.84$num.acc.1/nrep
num.acc.284 = tem.84$num.acc.2/nrep
num.acc.384 = tem.84$num.acc.3/nrep
time.84 = tem.84$time

n=80
T=80
tem.88 = nobreak.sim.func(n,T,p,taus,cpt.num,nrep,eps,cpt.eps,step.max)
num.acc.188 = tem.88$num.acc.1/nrep
num.acc.288 = tem.88$num.acc.2/nrep
num.acc.388 = tem.88$num.acc.3/nrep
time.88 = tem.88$time

save.image("b0.RData")