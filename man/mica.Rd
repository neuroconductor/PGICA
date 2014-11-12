\name{mica}
\alias{mica}
\title{
This is the main function for Parallel ICA algorithm
}
\description{
This function is the implementation of a parallel likelihood-based ICA algorithm. It can run in parallel on clusters and multi-core personal computers.
}
\usage{
mica(W0=0,n.b=1000,maxit=200,maxN=75,l.b=-10,u.b=10,N0=19,epsilon=10^(-4),
hc=0,alpha=1,ind=500,nproc=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W0}{
W0 is the initial value for mixing matrix W. It is calculated from function StVal().
}
  \item{n.b}{
n.b is the number of bins when estimating densities of underlying independent components.
}
  \item{maxit}{
maxit is the maximum number of iterations for the algorithm.
}
  \item{maxN}{
maxN is the maximum number of Gaussian distributions used to estimate density.
}
  \item{l.b}{
l.b is a lower bound of the interval on which we will estimate density.
}
  \item{u.b}{
u.b is an upper bound of the interval on which we will estimate density.
}
  \item{N0}{
N0 is the starting number of Gaussian distributions for density estimation.
}
  \item{epsilon}{
epsilon is the tolerance for the difference of estimations between 2 successive iterations.
}
  \item{hc}{
hc is not used here.
}
  \item{alpha}{
alpha is the step size for Newton algorithm.
}
  \item{ind}{
ind is a tuning parameter to speed up big matrix multiplication.
}
  \item{nproc}{
nproc is the number of parallel processes to use.
}
}
\details{
The Independent Component Analysis (ICA) model can be written as: X=AS, where X is T by V matrix,
A is T by T square matrix and S is T by V matrix. We assume that the rows of S are independent. 
Denote the inverse of A as W, then with W and X we can calculate S. The goal of our parallel ICA
algorithm is to calculate A (thus S) from X with independence assumption. In neuroimaging analysis,
The rows of S can be intepreted as functional networks.
}
\value{
mica returns a list of estimations. The first component of the list is the value of W.
}
\references{
Ani Eloyan, Ciprian M.Crainiceanu and Brian S.Caffo: Likelihood Based Population Independent Component Analysis
}
\author{
Ani Eloyan, Shaojie Chen, Lei Huang and Huitong Qiu et.al.
}

\examples{
%%Example 1: simulation example
######################  Generating simulation data
n=2000
m=2
N.s=3

A.true1=matrix(c(0.75,0.25,0.5,-0.5), 2, 2, byrow=TRUE)
A.true2=matrix(c(1,0,0.5,-0.5), 2, 2, byrow=TRUE)
A.true3=matrix(c(1,0.5,0.75,1), 2, 2, byrow=TRUE)

W.true1=solve(A.true1)
W.true2=solve(A.true2)
W.true3=solve(A.true3)

S1=rgamma(n,shape=4,scale=0.25)
S2=rgamma(n,shape=2,scale=2)
S.true=cbind(S1,S2) 
tr=PGICA:::trans(S.true,W.true1)
S.true=tr[[1]]

X.true1=S.true\%*\%A.true1
X.true2=S.true\%*\%A.true2
X.true3=S.true\%*\%A.true3

require(fastICA)
f1=fastICA(X.true1,n.comp=m)
f2=fastICA(X.true2,n.comp=m)
f3=fastICA(X.true3,n.comp=m)

X=c(list(t(f1$X\%*\%f1$K)),list(t(f2$X\%*\%f2$K)),list(t(f3$X\%*\%f3$K)))
K=c(list(f1$K),list(f2$K),list(f3$K))

X.full=c()
for(i in 1:N.s){
X.full=cbind(X.full,t(X[[i]]))
}

f=fastICA(X.full,n.comp=m)
tr=PGICA:::trans(f$S,f$W)
S.f=tr[[1]]
W1=solve(f$A[,1:m])
W2=solve(f$A[,(m+1):(2*m)])
W3=solve(f$A[,(2*m+1):(3*m)])

XtX=t(X.full)\%*\%X.full
sv=svd(XtX)
Sigma=diag(sqrt(sv$d))
SigmaInv=diag(1/sqrt(sv$d))
U=sv$u
V=X.full\%*\%U\%*\%SigmaInv

V.l=V[,1:m]
Sigma.l=Sigma[1:m,1:m]
U.l=t(U)[1:m,]
X.app=V.l\%*\%Sigma.l\%*\%U.l

X.a=c()
W0=c()
for(i in 1:N.s){
	X.a[[i]]=t(X.app[,((i-1)*m+1):(i*m)])
	W0=c(W0,list(t(solve(Sigma[1:m,1:m]\%*\%U.l[,((i-1)*m+1):(i*m)]))))
}

f=fastICA(X.full,n.comp=m)
S.f=f$S

dir.create('./Sim')
setwd("./Sim")

for(i in 1:N.s){
	X=X.a[[i]]
	save(X,file=paste(paste("app",i,sep=""),".rda",sep=""))
	}

save(W0,A.true1,A.true2,A.true3,K,W1,W2,W3,file="W0.rda")
save(S.f,S.true,file="Sfiles.rda")
setwd("..")

####################  Use PGICA to analyze above generated data
maxit=100
maxN=50
N0=19
epsilon=10^-3
W0=0
hc=0
u.b=10
l.b=-10
alpha=0.5
require(fastICA)
	
m=2
n.b=1000
N.s=3
n=2000
ind=500
setwd("./Sim")
fileDir=getwd()
files = dir(fileDir, pattern = "*", full.names = TRUE)
files=files[c(3:5,1,2)]
		
load("./W0.rda")
load(files[1])
f=fastICA(t(X),n.comp=m)
tr=PGICA:::trans(f$S,t(f$A))
S.f=t(tr[[1]])
A.f=t(tr[[2]])
	
for(subj in 2:N.s){
	W.temp=solve(W0[[1]]\%*\%t(A.f))\%*\%W0[[subj]]
	W0[[subj]]=W.temp
	}
W0[[1]]=solve(t(A.f))
require(parallel)
res=mica(W0,n.b,maxit,maxN,l.b,u.b,N0,epsilon,hc,alpha,ind,nproc=2)
tmp=((solve(res[[1]][[1]]))\%*\%S.f)
#hist(tmp[1,])
#hist(tmp[2,])
#hist(S.true[,1])
#hist(S.true[,2])

%% Example 2: Real Data Example
#This real data example is time-consuming
#m=20;
#N.s=1;
#alpha=0.5;
#setwd("..")
#dir.create('./data')
#data(PC,package="PGICA")
#save(PC,file="./data/sample.rda")
#StVal("./data/",m=20,N.s=1,V=30000)
#fileDir=getwd()
#files = dir(fileDir, pattern = "*.rda", full.names = TRUE)
#nfiles = length(files)
#outfile="m20.rda"
#setwd(fileDir)
#load("W0.rda")
#n=30000

#maxit=80
#maxN=50
#N0=19
#ind=100
#epsilon=10^-3
#hc=0
#u.b=10
#l.b=-10
#n.b=1000
#files=files[c(2,1)]
#res=mica(W0,n.b,maxit,maxN,l.b,u.b,N0,epsilon,hc,alpha)
}

\keyword{ ICA }
\keyword{ parallel computing }
