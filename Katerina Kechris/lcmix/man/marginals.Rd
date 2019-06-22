\name{marginals}

\alias{marginals}
\alias{marginals.mdmixmod}

\title{
Extract Marginal Models
}

\description{
Extract the marginal models from an object of class \code{mdmixmod}.
}

\usage{
\method{marginals}{mdmixmod}(x, ...)
}

\arguments{
\item{x}{an object of class \code{mdmixmod}.}
\item{...}{currently unused.}
}

\details{
\code{marginals} fits and returns the marginal models from \code{x}, that is, objects of class \code{mixmod} for each of the data sources used in the fitting of \code{x}, with the corresponding distribution.  
}

\value{A list.}

\author{
Daniel Dvorkin
}

\examples{\dontrun{
data(CiData)
joint <- mdmixmod(CiData, c(2,3,2), family=c("pvii", "norm", "norm"),
    topology="chained")
marginals(joint)
# $binding
# PVII mixture model ('pvii')
# Data 'binding' of length 10244 fitted to 2 components
# Model statistics:
#      iter      llik      qval       bic    iclbic 
#    353.00 -14409.37 -20330.50 -28883.38 -40725.65 
# 
# $expression
# Normal mixture model ('mvnorm')
# Data 'expression' of size 10244-by-4 fitted to 3 components
# Model statistics:
#       iter       llik       qval        bic     iclbic 
#      42.00  -47499.54  -50052.71  -95405.40 -100511.73 
# 
# $conservation
# PVII mixture model ('pvii')
# Data 'conservation' of length 10244 fitted to 2 components
# Model statistics:
#       iter      llik      qval       bic    iclbic 
#     897.00 -14259.84 -17797.79 -28584.32 -35660.23
}}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}}, \code{\link{posteriors}}
}
