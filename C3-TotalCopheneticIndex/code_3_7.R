library(Zseq)
library(gmp)
library(ape)
library(CollessLike)

big.factorial = function(n){
  if(n<2) return(1)
  return(Factorial(n+1)[n+1])
}

big.double.factorial = function(n){
  if(n<2) return(1)
  m = (n+2+n%%2)/2
  return(Factorial.Double(m,odd=(n%%2==1))[m])
}

big.binomial = function(n,k){
  return(big.factorial(n)/(big.factorial(k)*big.factorial(n-k)))
}

Cknk = function(k,n){
  return(big.binomial(n,k)*((big.double.factorial(2*k-3)*
          big.double.factorial(2*(n-k)-3))/(2*big.double.factorial(2*n-3))))
}

EUPhi = function(n){
    return(as.numeric((n*(n-1)/4)*(big.double.factorial(2*n-2)/
                                     big.double.factorial(2*n-3)-2)))
}

EUS = function(n){
    return(as.numeric(n*(big.double.factorial(2*n-2)/big.double.factorial(2*n-3)-1)))
}

term.cov = function(n){
  return(((13*n^2-9*n-2)*2^(n-5)*big.factorial(n))/(big.double.factorial(2*n-3))-
           (n*(n-1)/2)*(5*n-2)) 
} 

compute.EUcov = function(n.max=500){
  terms = lapply(2:n.max,term.cov)
  terms = c(0,terms)
  exp.values = list(0)
  for(n in 2:n.max){
    sums = 0
    if(n>2){ 
      for(k in 2:(n-1)){ 
        sums = sums + Cknk(k,n)*exp.values[[k]]
      }
      sums = 2*sums
    }
    sums = sums + terms[[n]]
    exp.values[[n]] = sums
    print(n)
  }
  exp.values = sapply(exp.values, as.numeric)
  write.table(exp.values,file=paste("E_U(SxPhi)_",n.max,".txt",sep = ""),
              row.names = F,col.names = F)
  return(exp.values)
}

compute.cov = function(exp.values,n.max = 500){
  cov.form = function(i)return(exp.values[i]-EUS(i)*EUPhi(i))
  cov.values = sapply(1:n.max, cov.form)
  write.table(cov.values,file = paste("cov_U(S_Phi)_",n.max,".txt",sep = ""),
              row.names = F,col.names = F)
  return(cov.values)
}

exp.values.cov = compute.EUcov(1000)
cov.values = compute.cov(exp.values.cov,1000)

exp.values.cov[3:20]
cov.values[3:20]

summary(lm(log(cov.values[900:1000])~log(900:1000)))

plot(log(1:1000),log(cov.values),xlab="Log of the number of leaves",
     ylab="Log of the variance")
reg.cov=lm(log(cov.values[500:1000])~log(500:1000))
abline(reg.cov,col="violet",lwd=2) 

##From other sections:
#trees = list()
#all.sackin.index = list() 
#all.cophen.index = list()
#for(n in 3:8){
#  trees[[n]] = read.tree(file = paste("trees_n",n,".txt",sep = ""))
#  all.sackin.index[[n]] = sapply(trees[[n]],sackin.index)
#  all.cophen.index[[n]] = sapply(trees[[n]],cophen.index)
#}

covariancesU = function(n){
  len = length(all.sackin.index[[n]])
  value = cov(all.sackin.index[[n]],all.cophen.index[[n]]*(len-1)/len)
  return(value)
}
real.cov.values = sapply(3:7,covariancesU)

pearson.cor = function(n.max){
  return(cov.values[4:n.max]/sqrt(var.values.S[4:n.max]*var.values.Phi[4:n.max]))
}
pearson.cor(20)