hausdorff.func <- function(A,B){
  hAB = NULL
  hBA = NULL
  Asize = length(A)
  Bsize = length(B)
  AA = matrix(rep(A,Bsize),nrow = Bsize)
  BB = matrix(rep(B,Asize),ncol = Asize,byrow=FALSE)
  for(i in (1:Bsize)){
    hAB=c(hAB,min(abs(AA-BB)[i,]))
  }
  hAB = max(hAB)
  for(i in (1:Asize)){
    hBA=c(hBA,min(abs(t(BB-AA)[i,])))
  }
  hBA = max(hBA)
  distance = max(hAB,hBA)
}