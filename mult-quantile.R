
rq.fit.lasso2 <- function (x, y, taus, lambda, ws, beta = 0.9995, eps = 1e-06) 
{
  n <- length(y)
  p <- ncol(x)
  if (n != nrow(x)) 
    stop("x and y don't match n")
  if (length(lambda) == 1) 
    lambda <- rep(lambda, p)
  if (length(lambda) != p) 
    stop(paste("lambda must be either of length ", p, " or length one"))
  lambda = lambda*ws
  if (any(lambda < 0)) 
    stop("negative lambdas disallowed")
  R <- diag(lambda, nrow = length(lambda))
  R <- R[which(lambda != 0), , drop = FALSE]
  r <- rep(0, nrow(R))
  
  X <- rbind(x, R)
  Y <- c(y, r)
  N <- length(Y)
  rhs <- t(x)%*%(1-taus) + 0.5 * apply(R, 2, sum)
  d <- rep(1, N)
  u <- rep(1, N)
  wn <- rep(0, 10 * N)
  wn[1:N] <- rep(0.5, N)
  z <- .Fortran("rqfnb", as.integer(N), as.integer(p), a = as.double(t(as.matrix(X))), 
                c = as.double(-Y), rhs = as.double(rhs), d = as.double(d), 
                as.double(u), beta = as.double(beta), eps = as.double(eps), 
                wn = as.double(wn), wp = double((p + 3) * p), aa = double(p * 
                                                                            p), it.count = integer(3), info = integer(1), PACKAGE = "quantreg")
  if (z$info != 0) 
    stop(paste("Error info = ", z$info, "in stepy2: singular design"))
  coefficients <- -z$wp[1:p]
  names(coefficients) <- dimnames(x)[[2]]
  residuals <- y - x %*% coefficients
  it.count <- z$it.count
  list(coefficients = coefficients, residuals = residuals, 
       taus = taus, lambda = lambda, it = it.count)
}

PQRcpt.func <-  function(X,y,taus,n,T,p,cpt.eps=1e-5,eps=1e-5,lambda,step.max = 50)
{
  require(quantreg)
  ntau = length(taus)
  Taus = rep(taus, each=n*T)
  ### Step 1, intial estimator for beta
  Y = rep(y,ntau)
  X1 = matrix(0, nrow=n*T, ncol=p*T)
  X1 = cbind(X, matrix(0, nrow=n*T, ncol=p*T)) 
  idx = sapply(c(1:n), function(r){ seq((r*T*(p*T+p)-p*T+1), (r*T*(p*T+p)), by=1) } )
  X1.s = matrix(as.vector(t(X1))[-idx], nrow=n*T, ncol=p*T, byrow=T)
  X1.ss = diag(ntau)%x%X1.s
  s <- rep(1:n,each=T)
  n <- length(unique(s))    
  nT <- length(s)           
  Z <- model.matrix(s~as.factor(s)-1)
  ZZ = rep(1,ntau)%x%Z
  X.ss = cbind(X1.ss,ZZ)##design matrix
  tmp = rq.fit.lasso2(X.ss,Y,taus=Taus,lambda = c(rep(0,dim(X1.ss)[2]),rep(T,dim(ZZ)[2])),ws=1) 
  beta.ini = matrix(tmp$coefficients[1:dim(X1.ss)[2]],ncol=p,byrow=T)
  fix.1 = matrix(tmp$coefficients[-(1:dim(X1.ss)[2])],ncol=1)
  ## initial weights
  ws=NULL
  for(r in 1:ntau){
    ws.1 = abs(beta.ini[((r-1)*T+2):(r*T), ]-beta.ini[((r-1)*T+1):(r*T-1), ])
    ws.2 = 1/apply(ws.1, 1, sum)
    ws = c(ws,0,ws.2)
  }
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
  X2.s = diag(ntau)%x%X2
  X2.ss = cbind(X2.s,ZZ)
  ### post-lasso
  tol = 1
  step = 0
  beta.post.tem = 0
  while (tol > eps && step < step.max){
    G2 = c(rep(ws*n*lambda,each=p),rep(T*lambda,n))
    tmp.2 = rq.fit.lasso2(X2.ss,Y,taus=Taus,lambda = G2,ws=1)
    beta.s = matrix(tmp.2$coefficients[1:dim(X2.s)[2]],ncol=p,byrow=T)
    L = apply(beta.s*beta.s,1,sum)
    LL = sqrt(L)
    cpt.1 = c(1:(ntau*T))[LL > cpt.eps]
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
    X2.post.ss = cbind(X2.post.s,ZZ)
    
    G2.post = c(rep(0,dim(X2.post.s)[2]),rep(T*lambda,dim(ZZ)[2]))
    tmp.post = rq.fit.lasso2(X2.post.ss,Y,taus=Taus,lambda = G2.post,ws=1)
    gamma.post = matrix(tmp.post$coefficients[1:dim(X2.post.s)[2]],ncol=p,byrow=T)
    fix.post = matrix(tmp.post$coefficients[-(1:dim(X2.post.s)[2])],byrow=T)
    beta.post = unlist(c(sapply(c(1:(length(cpt.1)-1)),function(r){
      beta.re = cpt.1[r+1]-cpt.1[r]
      rep(gamma.post[r,],beta.re)})))
    beta.post = matrix(beta.post,ncol=p,byrow=T)
    tol = sqrt( sum(beta.post-beta.post.tem)^2 )
    beta.post.tem = beta.post
    ###adaptive weight
    ws=NULL
    for(r in 1:ntau){
      ws.1 = abs(beta.post[((r-1)*T+2):(r*T), ]-beta.post[((r-1)*T+1):(r*T-1), ])
      ws.2 = 1/apply(ws.1, 1, sum)
      ws = c(ws,0,ws.2)
    }
    ws[ws==Inf] = 1e+6 ##?
    step = step+1
  }
  coef.final = c(as.vector(t(beta.post)),fix.post)
  cpt = cpt[-length(cpt)]
  error = Y-X.ss%*%coef.final
  error[error>0] = error[error>0]*Taus[error>0]
  error[error<0] = error[error<0]*(Taus[error<0]-1)
  list(coef = coef.final,error = error,cpt = cpt, fix = fix.post)
}

IC.func<-function(X,y,taus,n,T,p,cpt.eps=1e-6,eps=1e-6,lambda.range,step.max = 50){
  IC.s = Inf
  rho = 0.2*log(min(n,T))/min(n,T)
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
# n=10
# T=10
# p=2
# step.max = 50
# lambda=20
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