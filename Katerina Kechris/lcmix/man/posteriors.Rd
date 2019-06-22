\name{posteriors}

\alias{posteriors}
\alias{posteriors.mdmixmod}

\title{
Joint And Marginal Model Posterior Probabilities
}

\description{
Extract the joint and marginal model posterior probabilities from an object of class \code{mdmixmod}.
}

\usage{
\method{posteriors}{mdmixmod}(x, ...)
}

\arguments{
\item{x}{an object of class \code{mdmixmod}.}
\item{...}{currently unused.}
}

\details{
\code{posteriors} returns the joint posterior probability for \code{x} and the marginal probabilities for each of the data sources used in the fitting of \code{x}.  
}

\value{A list.  The first element has the name \code{JOINT}, indicating the posterior probability for the joint model; the other elements have the names of the data sources used in the fitting of \code{x}.}

\author{
Daniel Dvorkin
}

\examples{\dontrun{
data(CiData)
joint <- mdmixmod(CiData, c(2,3,2), family=c("pvii", "norm", "norm"),
    topology="chained")
pj = posteriors(joint)
lapply(pj, function(p) p[1:3,])
# $JOINT
#           [,1]      [,2]
# [1,] 0.0500650 0.9499350
# [2,] 0.1592164 0.8407836
# [3,] 0.5812072 0.4187928
# 
# $binding
#           [,1]      [,2]
# [1,] 0.1924187 0.8075813
# [2,] 0.2221210 0.7778790
# [3,] 0.8393147 0.1606853
# 
# $expression
#           [,1]         [,2]       [,3]
# 1 0.0003780798 9.543479e-01 0.04527400
# 2 0.3904756217 2.536705e-18 0.60952438
# 3 0.0008795689 9.806697e-01 0.01845077
# 
# $conservation
#           [,1]       [,2]
# [1,] 0.5916767 0.40832326
# [2,] 0.8832577 0.11674231
# [3,] 0.9465522 0.05344779
}}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}}, \code{\link{marginals}}
}
