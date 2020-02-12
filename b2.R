
## heteroskedasticity model
rm(list=ls())
setwd("F:\\ARCH研究\\code\\Mult-pqr\\新方法")
source("Mult-quantile.R")
source("twobreaksimfunc.R")
source("hausdorff.r")
setwd("F:\\ARCH研究\\code\\Mult-pqr\\新方法\\hete error")
source("DGP.R")

nrep = 250
p=2
taus=c(0.25,0.5,0.75)
cpt.num = 2
cpt.eps = 1e-6
eps = 1e-6
step.max = 50
set.seed(12345)
seeds = sample(1:999999, nrep, replace=FALSE)

n=T=40
tem.44 = twobreak.sim.func(n,T,p,taus,cpt.num,nrep,eps,cpt.eps,step.max)
num.acc.144 = tem.44$num.acc.1/nrep
num.acc.244 = tem.44$num.acc.2/nrep
num.acc.344 = tem.44$num.acc.3/nrep
date.acc.144 = tem.44$date.acc.1
date.acc.244 = tem.44$date.acc.2
date.acc.344 = tem.44$date.acc.3
time.44 = tem.44$time

n=40 
T=80
tem.48 = twobreak.sim.func(n,T,p,taus,cpt.num,nrep,eps,cpt.eps,step.max)
num.acc.148 = tem.48$num.acc.1/nrep
num.acc.248 = tem.48$num.acc.2/nrep
num.acc.348 = tem.48$num.acc.3/nrep
date.acc.148 = tem.48$date.acc.1
date.acc.248 = tem.48$date.acc.2
date.acc.348 = tem.48$date.acc.3
time.48 = tem.48$time

n=80
T=40
tem.84 = twobreak.sim.func(n,T,p,taus,cpt.num,nrep,eps,cpt.eps,step.max)
num.acc.184 = tem.84$num.acc.1/nrep
num.acc.284 = tem.84$num.acc.2/nrep
num.acc.384 = tem.84$num.acc.3/nrep
date.acc.184 = tem.84$date.acc.1
date.acc.284 = tem.84$date.acc.2
date.acc.384 = tem.84$date.acc.3
time.84 = tem.84$time

n=80
T=80
tem.88 = twobreak.sim.func(n,T,p,taus,cpt.num,nrep,eps,cpt.eps,step.max)
num.acc.188 = tem.88$num.acc.1/nrep
num.acc.288 = tem.88$num.acc.2/nrep
num.acc.388 = tem.88$num.acc.3/nrep
date.acc.188 = tem.88$date.acc.1
date.acc.288 = tem.88$date.acc.2
date.acc.388 = tem.88$date.acc.3
time.88 = tem.88$time

save.image("b2.RData")