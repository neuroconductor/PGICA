\name{PC}
\alias{PC}
\docType{data}
\title{
A sample fMRI data in .rda format.
}
\description{
This is a subsample of real fMRI data from the 1000 Connectome Project.
}
\usage{data(PC)}
\format{
PC is a 30000 by 20 matrix, where 30000 is the number of voxels in brain. 
}
\details{
This is not the original data, but the data after a dimension reduction(from 67749 by 100 to 67749 by 20) using PCA. Then a sample of 30000 voxels are picked out of 67749. The sample is achieved by setting the seed to be 0.
}
\source{
http://www.nitrc.org/projects/fcon_1000/
}
\references{
Ani Eloyan, Ciprian M.Crainiceanu and Brian S. Caffo: Likelihood Based Population Independent Component Analysis
}
\examples{
data(PC)
}
\keyword{datasets}
