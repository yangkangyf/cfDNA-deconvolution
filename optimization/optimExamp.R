
N = 1000 #number individuals
M = 3 #number of parameters without intercept

X = matrix(rnorm(N*M),nrow=N,ncol=M)

Xp = cbind(rep(1,N),X)

betas = runif(M+1)

Y = Xp%*%betas + rnorm(N)


LRLogLik <- function(ibetas,iY,iX) {
    N = length(iY)
    ev = iX%*%ibetas
    sig2 = var(iY-ev)
    ret = -N/2*log(2*pi)-N/2*log(sig2) - 1/(2*sig2)*(sum((iY-ev)^2))
    return(-ret)
}


ores= optim(rep(0,4),LRLogLik,iY=Y,iX=Xp)


#methylation of tissues
Nt = 3 #number of tissues
Nm = 10 #number of methylation sites
Ni = 1 #number of individuals
Rd = 10 #read depth

#methylatio precentage of each tis at each CpG
R = matrix(runif(Nt*Nm),nrow=Nt,ncol=Nm,byrow=T) 

Mtemp = matrix(runif(Nt*Ni),nrow=Ni,ncol=Nt,byrow=T)

#percentage of tissue of each individual
M = matrix(0,nrow=Ni,ncol=Nt)
for (i in 1:Ni) {
    ci = Mtemp[i,]
    M[i,] = ci/sum(ci)
}

#get some noise
eps = matrix(rnorm(Ni*Nm,0,0.01),nrow=Ni,ncol=Nm,byrow=T)

#methylation of each CpG in each individual
O = M%*%R+eps

#make read depths
D = matrix(rpois(Nm*Ni,Rd),nrow=Ni,ncol=Nm,byrow=T)

D
#draw counts
C = matrix(0,nrow=Ni,ncol=Nm)
C
for (i in 1:Ni) {
  print(D[i,])
  print(O[i,])
  C[i,] = rbinom(Nm,D[i,],O[i,])
  print(C[i,])
}

#the goal is to get back M for each individual (at this point independently). FS we start with O[1,]
M1  = C[1,]
U1 = D[1,]-C[1,]

iguess = runif(Nt)
iguess = iguess/sum(iguess)

O1est = iguess%*%R
O1est 

#methylated sites, unmethylated sites, probability of methylation
logLikCounts <- function(m1,u1,mguess,R) {
  pm = mguess%*%R
  ll = sum(dbinom(m1, (m1+u1), pm, log = T))
  return(-ll)
}

logLikCounts(M1,U1,iguess,R)

#now optimizie
iguess
res = optim(iguess,logLikCounts,m1=M1,u1=U1, R=R, upper=1,lower=0,method="L-BFGS-B")

res
res$par
sum(res$par)
cbind(res,M[1,])
res

TCLik <- function(TCP,T,O) {
  TCP = TCP/sum(TCP)
  eO = t(TCP%*%T)
  N = length(TCP)      
  print(N)
  sig2 = var(O-eO) #this might be variance of O1, check later
  ret = -N/2*log(2*pi)-N/2*log(sig2) - 1/(2*sig2)*(sum((O-eO)^2))
  return(-ret)
}

length(O1)