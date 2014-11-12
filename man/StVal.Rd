\name{StVal}
\alias{StVal}
\title{
Computing Starting Values for PGICA
}
\description{
The main function uses a iterative algorithm. This function calculates the starting values for the algorithm.
}
\usage{
StVal(fileDir, m = 20, N.s = 150, V = 67749)
}
\arguments{
  \item{fileDir}{
fileDir is a file directory containing the original .rda files.
}
  \item{m}{
m is the number of independent components you want.
%%     ~~Describe \code{m} here~~
}
  \item{N.s}{
N.s is the number of subjects.
%%     ~~Describe \code{N.s} here~~
}
  \item{V}{
V is the number of observed signals. In the fMRI case, it is the number of voxels.
%%     ~~Describe \code{V} here~~
}
}
\details{
We used the PVD algorithm by C.M.Crainiceanu et.al to compute starting values for our algorithm.
}
\value{
This function has no returned values, instead it will write the output into .rda files. The W0.rda file contains starting values for W matrix and the appi.rda files are the starting values for S. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Ciprian M.Crainiceanu, Brian S.Caffo, Sheng Luo, Vadim Zipunikov
Population Value Decomposition, A Framework for the Analysis of Image Populations
\url{http://biostats.bepress.com/jhubiostat/paper220/}
%% ~put references to the literature/web site here ~
}
\author{
Ani Eloyan, Shaojie Chen, Lei Huang, Huitong Qiu
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
dir.create('./data')
data(PC,package="PGICA")
save(PC,file="./data/sample.rda")
StVal("./data/",m=20,N.s=1,V=30000)
# You will find files W0.rda and app1.rda under your current working directory}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Starting Value}

