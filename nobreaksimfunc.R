nobreak.sim.func <- function(n,T,p,taus,cpt.num,nrep,eps,cpt.eps,step.max){
  ptm <- proc.time();
  
  num.acc.1 = num.acc.2 = num.acc.3 = 0
  cpts.1=cpts.2=cpts.3= vector("list",nrep)
  coefs.pre = coefs = vector("list",nrep)
  num.false.1 = num.false.2 = num.false.3 = NULL
  for( i in c(1:nrep)){
    
    data.tem = DGP.func(n,T,cpt.num,i)
    y = data.tem$y
    X = data.tem$X
    beta = data.tem$beta
    
    lambda.range = seq(0.01,3,length.out=20)
    tuning.tem = IC.func(X,y,taus,n,T,p,cpt.eps,eps,lambda.range,step.max = 50)
    lambda = tuning.tem$lambda
    IC = tuning.tem$IC
    tem = PQRcpt.func(X,y,taus,n,T,p,cpt.eps,eps,lambda,step.max = 50)
    cpt = tem$cpt
    ## divide into three tau quantile level
    cpt.idx = which(cpt==1)
    cpt.1 = cpt[1:(cpt.idx[2]-1)]
    cpt.2 = cpt[cpt.idx[2]:(cpt.idx[3]-1)]
    cpt.3 = cpt[cpt.idx[3]:length(cpt)]
    cpt.1 = cpt.1[-(which(cpt.1==1))]
    cpt.2 = cpt.2[-(which(cpt.2==1))]
    cpt.3 = cpt.3[-(which(cpt.3==1))]
    
    
    if(length(cpt.1) == (length(beta)-1)){
      num.acc.1 = num.acc.1+1
    }else{
      num.false.1 = c(num.false.1,i)
    }
    if(length(cpt.2) == (length(beta)-1)){
      num.acc.2 = num.acc.2+1
    }else{
      num.false.2 = c(num.false.2,i)
    }
    if(length(cpt.3) == (length(beta)-1)){
      num.acc.3 = num.acc.3+1
    }else{
      num.false.3 = c(num.false.3,i)
    }
    cpts.1[[i]] = as.vector(cpt.1)
    cpts.2[[i]] = as.vector(cpt.2)
    cpts.3[[i]] = as.vector(cpt.3)
    
    coefs.pre[[i]] = tem$coef.pre
    coefs[[i]] = tem$coef
    print(i)
  }
  time=proc.time()-ptm
  list(num.acc.1 = num.acc.1, num.acc.2 = num.acc.2, num.acc.3 = num.acc.3,
       cpts.1 = cpts.1,cpts.2 = cpts.2 , cpts.3 = cpts.3,
       num.false.1 = num.false.1, num.false.2 = num.false.2,num.false.3 = num.false.3,
       time = time,coefs = coefs,coefs.pre = coefs.pre)
}