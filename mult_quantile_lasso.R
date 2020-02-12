PQRcpt.func <-  function(X,y,taus,n,T,p,cpt.eps=1e-2,eps=1e-5,lambda,step.max = 50)
{
  require(quantreg)
  
  ntau = length(taus)
  Taus = rep(taus, each=n*T)
  
  ### Step 1, preliminary estimator for beta
  Y = rep(y,ntau)
  X1 = cbind(X, matrix(0, nrow=n*T, ncol=p*T)) 
  idx = sapply(c(1:n), function(r){ seq((r*T*(p*T+p)-p*T+1), (r*T*(p*T+p)), by=1) } )
  X1.s = matrix(as.vector(t(X1))[-idx], nrow=n*T, ncol=p*T, byrow=T)
  X1.com = diag(ntau)%x%X1.s
  s <- rep(1:n,each=T)
  n <- length(unique(s))    
  nT <- length(s)           
  Z <- model.matrix(s~as.factor(s)-1)
  Z.com = rep(1,ntau)%x%Z
  XZ.com = cbind(X1.com,Z.com)##design matrix
  
  tem_pre = rq.fit.lasso.composite(XZ.com,Y,taus=Taus,
                                   lambda = c(rep(0,ntau*p*T),rep(T,n)),ws=1)
  
  beta.pre = matrix(tem_pre$coefficients[1:(ntau*p*T)],ncol=p,byrow=T)
  fix.pre = matrix(tem_pre$coefficients[-(1:(ntau*p*T))],ncol=1)
  
  ## initial weights
  ws.com = matrix(0,ncol = T,nrow = ntau)
  for(r in 1:ntau){
    ws = abs(beta.pre[((r-1)*T+2):(r*T), ]-beta.pre[((r-1)*T+1):(r*T-1), ])
    ws.tem = c(0,1/apply(ws, 1, sum))
    ws.com[r,] = ws.tem
  }
  ws.com = as.vector(t(ws.com))
  
  ###step 2, variable substitution
  X2 = NULL
  idx.x2 = NULL
  X2.tem = NULL
  for (r in 1:(T-1))
  {
    idx.x2 = c(idx.x2, seq(r, n*T, by=T))
    X2.tem.p = X
    if(p==1){
      X2.tem.p[idx.x2] = 0
    }else{
      X2.tem.p[idx.x2,] = 0
    }
    X2.tem = cbind(X2.tem, X2.tem.p)
  }
  X2 = cbind(X, X2.tem)
  X2.com = diag(ntau)%x%X2
  XZ2.com = cbind(X2.com,Z.com)
  
  ### post-lasso
  tol = 1
  step = 0
  beta.post.tem = 0
  
  while (tol > eps && step < step.max){
    
    G2 = c(rep(ws.com*n*lambda,each=p),rep(T*lambda,n))
    tem_fuse = rq.fit.lasso.composite(XZ2.com,Y,taus=Taus,lambda = G2,ws=1)
    theta = matrix(tem_fuse$coefficients[1:(ntau*p*T)],ncol=p,byrow=T)
    L = sqrt(apply(theta*theta,1,sum))
  
    cpt.1 = c(1:(ntau*T))[L > cpt.eps]
    cpt.1 = c(cpt.1,ntau*T+1)
    cpt = c(cpt.1%%T)
    cpt[which(cpt==0)]=T
    cpt.idx = which(cpt==1)
    
    X2.post.s = NULL
    for(q in 1:ntau){
      X2.post = NULL
      cpt.tem = c(cpt[cpt.idx[q]:(cpt.idx[q+1]-1)],T+1)
      cpt.n = length(cpt.tem)-1
      for(r in 1:cpt.n){
        X2.post.tem = matrix(0,ncol=p,nrow=ntau*n*T)
        for(k in 0:(n-1)){
          X2.post.tem[c(((q-1)*n*T+k*T+cpt.tem[r]):((q-1)*n*T+k*T+cpt.tem[r+1]-1)),] = X[c((k*T+cpt.tem[r]):(k*T+cpt.tem[r+1]-1)),]
        }
        X2.post = cbind(X2.post,X2.post.tem)
      }
      X2.post.s = cbind(X2.post.s,X2.post)
    }
    
    XZ3.com = cbind(X2.post.s,Z.com)
    
    G2.post = c(rep(0,dim(X2.post.s)[2]),rep(T*lambda,dim(Z.com)[2]))
    tem_post = rq.fit.lasso.composite(XZ3.com,Y,taus=Taus,lambda = G2.post, ws=1)
    
    gamma.post = matrix(tem_post$coefficients[1:dim(X2.post.s)[2]],ncol=p,byrow=T)
    
    fix.post = matrix(tem_post$coefficients[-(1:dim(X2.post.s)[2])],byrow=T)
    
    beta.post = unlist(c(sapply(c(1:(length(cpt.1)-1)),function(r){
      beta.re = cpt.1[r+1]-cpt.1[r]
      rep(gamma.post[r,],beta.re)})))
    beta.post = matrix(beta.post,ncol=p,byrow=T)
    
    tol = sqrt(sum(beta.post-beta.post.tem)^2)
    
    beta.post.tem = beta.post
    
    ###adaptive weight
    ws.com = matrix(0,ncol = T,nrow = ntau)
    for(r in 1:ntau){
      ws = abs(beta.post[((r-1)*T+2):(r*T), ]-beta.post[((r-1)*T+1):(r*T-1), ])
      ws.tem = c(0,1/apply(ws, 1, sum))
      ws.com[r,] = ws.tem
    }
    ws.com = as.vector(t(ws.com))
    ws.com[ws.com==Inf] = 1e+6
    
    step = step+1
  }
  
  coef.final = c(as.vector(t(beta.post)),fix.post)
  cpt = cpt[-length(cpt)]
  error = Y-XZ.com%*%coef.final
  error[error>0] = error[error>0]*Taus[error>0]
  error[error<0] = error[error<0]*(Taus[error<0]-1)
  list(coef = coef.final,error = error,cpt = cpt, fix = fix.post,coef.pre = tem_pre$coefficients)
}

IC.func<-function(X,y,taus,n,T,p,cpt.eps=1e-6,eps=1e-6,lambda.range,step.max = 50){
  IC.s = Inf
  rho = 0.5*log(min(n,T))/min(n,T)
  for (lambda in lambda.range){
    tem = PQRcpt.func(X,y,taus,n,T,p,cpt.eps=1e-6,eps=1e-6,lambda,step.max = 50)
    cpt.num = length(tem$cpt)-length(taus)
    error = sum(tem$error)
    fix = tem$fix
    IC = log(error)+rho*p*(cpt.num+length(taus))
    if(IC < IC.s){
      IC.s = IC
      lambda.s = lambda
    }
  }
  list(IC  =IC.s ,lambda = lambda.s)
}
# ##DGP
# rm(list = ls())
# n=10
# T=10
# p=2
# step.max = 50
# lambda=100
# eps =cpt.eps=1e-6
# taus=c(0.25,0.5,0.75)
# set.seed(1)
# x <- rnorm(n*T,0,2)
# X <- cbind(1,x)
# e <- rnorm(n*T)
# alpha = rep(rnorm(n),each=T)
# cpt.1 =NULL
# cpt.2 =NULL
# beta = c(1,2)
# cpt.1 = floor(T/2)
# Bcpt = c(rep(beta[1], cpt.1), rep(beta[2], (T-cpt.1)))
# B = cbind(1, rep(Bcpt, n))
# y <- alpha + apply(X*B, 1, sum) + e