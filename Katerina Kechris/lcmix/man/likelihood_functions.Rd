
\name{likelihood functions}

\alias{likelihood}

\alias{logLik.mixmod}

\alias{qval}
\alias{qval.mixmod}
\alias{qfun}
\alias{qfun.mixmod}
\alias{qfun.mdmixmod}

\alias{aic}
\alias{aic.mixmod}
\alias{bic}
\alias{bic.mixmod}
\alias{entropy}
\alias{entropy.mixmod}
\alias{iclbic}
\alias{iclbic.mixmod}
\alias{siclbic}
\alias{siclbic.mdmixmod}

\title{
Likelihood-Based Mixture Model Statistics
}

\description{
See Details for a description of the individual functions.
}

\usage{
\method{logLik}{mixmod}(object, ...)

\method{qval}{mixmod}(x, map=FALSE, ...)
\method{qfun}{mixmod}(x, map=FALSE, ...)
\method{qfun}{mdmixmod}(x, map=FALSE, ...)

\method{aic}{mixmod}(x, ...)
\method{bic}{mixmod}(x, ...)
\method{entropy}{mixmod}(x, map=FALSE, ...)
\method{iclbic}{mixmod}(x, map=FALSE, ...)
\method{siclbic}{mdmixmod}(x, map=FALSE, ...)
}

\arguments{
\item{x, object}{an object of class \code{\link{mixmod}} or \code{\link{mdmixmod}}.}
\item{map}{logical; if \code{TRUE}, the maximum a posteriori (MAP) estimates rather than the posterior probabilities will be used when estimating expectations with respect to hidden data.}
\item{...}{currently unused.}
}

\details{
\code{logLik} calculates \eqn{\mathcal L_X(\hat\theta)}{L(theta|X)}, the log-likelihood of the estimated parameters \eqn{\theta}{theta} with respect to the observed data \eqn{X}{X}, while \code{qval} calculates the \dQuote{Q-value}, the expectation with respect to the hidden data of the log-likelihood with respect to the complete data:  \eqn{\mbox Q(\hat\theta) = \mbox E [ \mathcal L_{X, \mathcal Y}(\hat\theta) ]}{Q(theta) = E[L(theta|X,Y)]} for \code{mixmod} and \eqn{\mbox Q(\hat\theta) = \mbox E [ \mathcal L_{X, \mathcal Y, \mathcal Y_0}(\hat\theta) ]}{Q(theta) = E[L(theta|X,Y,Y0)]} for \code{mdmixmod}.  \code{qfun} returns the hidden and observed portions of the Q-value separately, as elements of a vector.

\code{aic}, \code{bic}, \code{entropy}, \code{iclbic}, and \code{siclbic} calculate various information criteria for model selection with mixture models of class \code{mixmod} and \code{mdmixmod}.  These criteria are Akaike's information criterion (AIC, Akaike, 1974), the Bayes information criterion (BIC, Schwarz, 1978), the classification entropy (Biernacki et al., 2000), the integrated complete likelihood BIC (ICL-BIC, Biernacki et al., 2000), and the simplified ICL-BIC (SICL-BIC) for objects of class \code{mdmixmod}, respectively.  They are defined as follows:
\tabular{rcl}{
	AIC \tab = \tab \eqn{2 \mathcal L_X(\hat\theta) - 2 |\Theta|}{2 L(theta|X) - 2 |Theta|} \cr
	BIC \tab = \tab \eqn{2 \mathcal L_X(\hat\theta) - |\Theta| \log N}{2 L(theta|X) - |Theta| log(N)} \cr
	entropy \tab = \tab \eqn{2 \mathcal L_X(\hat\theta) - 2 \mbox Q(\hat\theta)}{2 L(theta|X) - 2 Q(theta)} \cr
	ICL-BIC \tab = \tab \eqn{2 \mbox Q(\hat\theta) - |\Theta| \log N}{2 Q(theta) - |Theta| log(N)} \cr
	SICL-BIC \tab = \tab \eqn{2 \mbox E[\mathcal L_{X, \mathcal Y_0}(\hat{\theta})] - |\Theta| \log N}{2 E[L(theta|X,Y0)] - |Theta| log(N)} (\code{mdmixmod} only)
}
where \eqn{|\Theta|}{|Theta|} is the size of the parameter space and \eqn{N}{N} is the size of the data.  Generally, the model which provides the highest value of any information criterion should be selected.  Current testing indicates that ICL-BIC is preferred for \code{mixmod} and BIC for \code{mdmixmod}.
}

\references{
Akaike, H. (1974)  A new look at the statistical model identification, \emph{IEEE Transactions on Automatic Control} \bold{19(6)}, 716--723.

Biernacki, C. and Celeux, G. and Govaert, G.  (2000)  Assessing a mixture model for clustering with the integrated completed likelihood, \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence} \bold{22(7)}, 719--725.

McLachlan, G.J. and Thriyambakam, K.  (2008)  \emph{The EM Algorithm and Extensions}, John Wiley & Sons.
    
Schwarz, G. (1978)  Estimating the dimension of a model, \emph{The Annals of Statistics} \bold{6(2)}, 461--464.
}

\value{
A numeric vector for \code{qfun}, a numeric scalar for the other functions.
}

\author{
Daniel Dvorkin
}

\note{
Some authors define AIC, BIC, and ICL-BIC as the negative of the quantities given in Details. 
}

\seealso{
\code{\link{mixmod}} and \code{\link{mdmixmod}} for details of the hidden variable structure.
}

\examples{\dontrun{
data(CiData)
fit <- mixmod(CiData$expression, 3)
bic(fit)            # -95405.4
qval(fit)           # -50055.35
qval(fit, map=TRUE) # -49738.53
}}

\keyword{models}
\keyword{cluster}
\keyword{classif}