\name{PGICA-package}
\alias{PGICA-package}
\alias{PGICA}
\docType{package}
\title{
Parallel Group ICA Algorithm
}
\description{
This package implements a Group ICA Alorithm which can run in parallel on clusters and multi-core personal computers. Unlike existing ICA algorithms, this parallel algorithm can analyze very big data. It can be used in many applications including signal processing and neuroimage data analysis.
}
\details{
\tabular{ll}{
Package: \tab PGICA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-12\cr
License: \tab GPL-2\cr
}
The two main functions in this package are StVal and mica. mica is the main ICA function and StVal calculates initial values for mica.
}
\author{
Ani Eloyan, Shaojie Chen, Lei Huang, Huitong Qiu

Maintainer:  <schen89@jhu.edu> Shaojie Chen
}
\references{
Ani Eloyan, Ciprian M. Crainiceanu and Brian S. Caffo: Likelihood Based Population Independent Component Analysis
}
\keyword{ICA}
\keyword{Parallel Computing}
