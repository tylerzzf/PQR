DGP.func<-function(n,T,cpt.num,i){
  set.seed(seeds[i])
  x <- rnorm(n*T,0,sqrt(2))
  X <- cbind(1,x)
  e <- rnorm(n*T,0,1)*(1+0.5*abs(x))
  alpha = rep(rnorm(n),each=T)
  cpt.1 =NULL
  cpt.2 =NULL
  if(cpt.num==0){
    beta = c(1)
    Bcpt = c(rep(beta[1], T))
    B = cbind(1, rep(Bcpt, n)) 
    y <- alpha + apply(X*B, 1, sum) + e ### X(11,12,..1t,21,...,nt)
  }
  if(cpt.num==1){
    beta = c(1,0)
    cpt.1 = floor(T/2)
    Bcpt = c(rep(beta[1], cpt.1), rep(beta[2], (T-cpt.1)))
    B = cbind(1, rep(Bcpt, n)) 
    y <- alpha + apply(X*B, 1, sum) + e ### X(11,12,..1t,21,...,nt) 
  }
  if(cpt.num==2){
    beta = c(1,0,1)
    cpt.1 = floor(T/3)+1
    cpt.2 = floor(2*T/3)+1
    Bcpt = c(rep(beta[1], (cpt.1-1)), rep(beta[2], (cpt.2-cpt.1)),rep(beta[3],(T-cpt.2+1)))
    B = cbind(1, rep(Bcpt, n)) 
    y <- alpha + apply(X*B, 1, sum) + e ### X(11,12,..1t,21,...,nt)
  }
  list(X=X,y=y,beta=beta,cpt.1 = cpt.1,cpt.2 = cpt.2)
}