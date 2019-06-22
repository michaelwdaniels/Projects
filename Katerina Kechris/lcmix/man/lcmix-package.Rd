\name{lcmix-package}

\alias{lcmix-package}
\alias{lcmix}

\concept{parameter estimation}
\concept{maximum likelihood estimation}
\concept{maximum likelihood estimator}
\concept{mixture model estimation}
\concept{EM algorithm}
\concept{expectation maximization algorithm}

\docType{package}

\title{
Layered and Chained Mixture Models
}

\description{
Fit layered and chained mixture models, which are special cases of Bayesian belief networks, to lists of data.
}

\details{
In the layered model, a single categorical hidden variable generates a layer of categorical hidden variables, each of which in turn generates an observed variable.  The chained model is a hidden Markov model with a non-homogeneous alphabet, that is, each of the categorical hidden variables may have its own sample space. 

Hidden variables have a categorical distribution.  Observed variables may take on any distribution.  The current release of the package provides for univariate and multivariate normal (Gaussian), Pearson Type VII, exponential, gamma, and Weibull distributions for observed variables.

\strong{Running under 64-bit \R is strongly recommended!}

Developed at the Computational Bioscience Program of the University of Colorado Denver | Anschutz Medical Campus.  Work supported by National Library of Medicine (National Institutes of Health) grant T15 LM009451.
} % line breaks inside paragraphs prevent Rd2txt from truncating lines

\author{
Daniel Dvorkin <daniel.dvorkin@gmail.com>
}

\references{
McLachlan, G.J. and Thriyambakam, K.  (2008)  \emph{The EM Algorithm and Extensions}, John Wiley & Sons.
}

\keyword{package}
\keyword{models}
\keyword{cluster}
\keyword{classif}

\seealso{
\code{\link{mdmixmod}} for fitting layered and chained models to multiple data sources, \code{\link{mixmod}} for fitting mixture models to single data sources; packages \code{mixtools} and \code{mclust}.
}
