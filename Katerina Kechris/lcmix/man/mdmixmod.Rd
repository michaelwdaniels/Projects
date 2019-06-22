\name{mdmixmod}

\alias{mdmixmod}
\alias{print.mdmixmod}

\concept{parameter estimation}
\concept{maximum likelihood estimation}
\concept{maximum likelihood estimator}
\concept{mixture model estimation}
\concept{EM algorithm}
\concept{expectation maximization algorithm}

\title{
Multiple Data Mixture Models
}

\description{
Fits a layered or chained mixture model to a list representing multiple sources of data, using a choice of distributions and number of components for each data source.
}

\usage{
mdmixmod(X, K, K0=min(K), topology=LC_TOPOLOGY, family=NULL, prior=NULL, 
    prefit=TRUE, iter.max=LC_ITER_MAX, dname=deparse(substitute(X)))
\method{print}{mdmixmod}(x, ...)
}

\arguments{
\item{X}{a list of observed data sources; the elements must be numeric vectors, matrices, or data frames.  Each element of \code{X} must have the same data length, that is, length of vectors or number of rows of matrices and data frames; but each element of \code{X} may be of arbitrary width, that is, number of columns of matrices and data frames.}
\item{K}{the vector of numbers of mixture components for the hidden variables corresponding to each observed data source.  If \code{length(K) < length(X)}, the values in \code{K} will be recycled; thus \code{K} may be a scalar, in which case the same number of components will be used throughout.  If \code{length(K) > length(X)}, then \code{K} will be truncated and a warning will be given.}
\item{K0}{the number of mixture components for the top-level hidden variable.}
\item{topology}{one of the model topologies in \code{\link{LC_TOPOLOGY}}, either \code{"layered"} or \code{"chained"}.  By default, \code{"layered"} is used.  Partial matches are allowed.}
\item{family}{a vector of names of distribution families to be used in fitting the models for each observed data source; each element of \code{family} must be one of the distribution names in \code{\link{LC_FAMILY}}.  The usual recycling and truncation rules are followed.  By default, \code{"normal"} is used and recycled to \code{c("normal", \dots, "normal")}.  Partial matches are allowed.}
\item{prior}{prior probability distribution on \eqn{\mathcal Y_0}{Y_0}.  This feature is under development and its use is not currently recommended.}
\item{prefit}{logical; if \code{TRUE}, the marginal models will be fitted first and the resulting weights used for initialization.}
\item{iter.max}{the maximum number of iterations for the EM algorithm, by default equal to \code{\link{LC_ITER_MAX}}.}
\item{dname}{the name of the data.}
\item{x}{an object of class \code{mixmod}.}
\item{...}{further arguments to \code{print.default}.}
}

\details{
In the layered model, a top-level hidden categorical random variable \eqn{\mathcal{Y}_0}{Y_0}, which can take on values from 1 to some positive integer \eqn{K_0}, generates categorical hidden random variables \eqn{\mathcal{Y}_1, \dots, \mathcal{Y}_Z}{Y_1, \dots, Y_Z} for some positive integer \eqn{Z}.  For \eqn{z = 1,\dots,Z}, each \eqn{\mathcal{Y}_z}{Y_z} can take on values from 1 to some positive integer \eqn{K_z}.  In the chained model, \eqn{\mathcal{Y}_0}{Y_0} generates \eqn{\mathcal{Y}_1}{Y_1}, which in turn generates \eqn{\mathcal{Y}_2}{Y_2}, etc., up to \eqn{\mathcal{Y}_{Z-1}}{Y_{Z-1}}, which generates \eqn{\mathcal{Y}_Z}{Y_Z}.

In both models, the \eqn{\mathcal{Y}_z}{Y_z}'s generate the observed mixture random variables \eqn{\mathcal{X}_1, \dots, \mathcal{X}_Z}{X_1, \dots, X_Z}, from which the elements of the observed data \code{X} are assumed to be drawn.  (That is, \code{Z = length(X)}, the number of list elements in \code{X}.)  The relationship between each \eqn{\mathcal{Y}_z}{Y_z} and \eqn{\mathcal{X}_z}{X_z} is the same as the relationship between \eqn{\mathcal{Y}}{Y} and \eqn{\mathcal{X}}{X} in \code{\link{mixmod}}.

As in \code{\link{mixmod}}, the EM algorithm attempts to maximize the Q-value, that is, the expected complete data (hidden and observed variables) log-likelihood.
}

\value{
A list of class \code{mdmixmod}, a subclass of \code{mixmod}, having the following elements:
\item{N}{the length of the data, that is, \code{length(X[[1]])} if \code{X[[1]]} is a vector, or \code{nrow(X[[1]])} if \code{X[[1]]} is a matrix or data frame.}
\item{Z}{the size of the data, that is, \code{Z = length(X)}.}
\item{D}{the vector of widths of the data, that is, \code{D[z] = 1} if \code{X[[z]]} is a vector, or \code{D[z] = ncol(X[[z]])} if \code{X[[z]]} is a matrix or data frame.}
\item{K}{the vector of the numbers of components in the lower-level mixture models.}
\item{K0}{the number of components in the top-level mixture model, that is, \eqn{K_0}{K_0}.}
\item{X}{the original data, with data frames converted to matrices.  If the elements of \code{X} were not named, they will be named \code{"X1",\dots,"XZ"} here.}
\item{npar}{the total number of parameters in the model.}
\item{npar.hidden}{the number of parameters for the hidden component portion of the model.}
\item{npar.observed}{the number of parameters for the observed data portion of the model.}
\item{iter}{the number of iterations required to fit the model.}
\item{params}{the parameters estimated for the model.  This is a list with elements \code{hidden} and \code{observed}, corresponding to distribution for the hidden and observed portions of the model.  \code{hidden} has elements \code{prob0}, the vector of probabilities for the possible values for \eqn{\mathcal{Y}_0}{Y_0}, and \code{cprob}, the list of matrices of conditional probabilities for the possible values of the \eqn{\mathcal{Y}_z}{Y_z}'s.  In the layered model, these are \eqn{K_0 \times K_z}{K_0-by-K_z} matrices of which the \eqn{(k_0,k_z)}th element represents \eqn{P(\mathcal{Y}_z = k_z \mid \mathcal{Y}_0 = k_0)}{P(Y_z = k_z | Y_0 = k_0)}.  In the chained model, the \eqn{z}th element of \code{cprob} is a \eqn{K_{z-1} \times K_z}{K_{z-1}-by-K_z} matrix of which the \eqn{(k_{z-1},k_z)}th element represents \eqn{P(\mathcal{Y}_z = k_z \mid \mathcal{Y}_{z-1} = k_{z-1})}{P(Y_z = k_z | Y_{z-1} = k_{z-1})}.  The chained model \code{hidden} also has elements \code{probz}, a list of vectors representing the marginal probabilities for the \eqn{\mathcal{Y}_z}{Y_z}'s, and \code{rprob}, representing the conditional probabilities of the \eqn{\mathcal{Y}_z}{Y_z}'s given the values of the \eqn{\mathcal{Y}_{z+1}}{Y_{z+1}}'s.  The elements of \code{observed} depend on the distribution family chosen in fitting the model.}
\item{stats}{a vector with named elements corresponding to the number of iterations, log-likelihood, Q-value, and BIC for the estimated parameters.}
\item{weights}{a list with elements \code{U}, a \eqn{N \times K_0}{N-by-K_0} matrix; \code{V}, a \eqn{Z}-length list of \eqn{K_0}-length lists of \eqn{N \times K_z}{N-by-K_z} matrices in the layered model, or a \eqn{Z}-length list of \eqn{K_{z-1}}-length lists of \eqn{N \times K_z}{N-by-K_z} matrices in the layered model; and \code{W}, a \eqn{Z}-length list of \eqn{N \times K_z}{N-by-K_z} matrices; representing weights used in the M-step of the EM algorithm for estimating the final set of parameters for the observed data portion of the model.}
\item{pdfs}{a list with elements \code{G}, \code{alpha}, \code{beta}, and \code{gamma}, representing various estimated density functions for the data.  \code{gamma} represents the estimated density of the observed data across all data sources under the fitted mixture model.}
\item{posterior}{the \eqn{N \times K_0}{N-by-K_0} matrix of which the \eqn{(n,k_0)}th element is the estimated posterior probability that the \eqn{n}th observation (across all data sources) was generated by the \eqn{k_0}th component.  Equal to the \code{U} element of \code{weights}.}
\item{assignment}{the vector of length \eqn{N} of which the \eqn{n}th element is the most probable top-level component to have generated the \eqn{n}th observation.  In other words, \code{assignment[n] = which.max(posterior[n,])}.}
\item{iteration.params}{a list of length \code{iter} giving the estimated parameters at each iteration of the algorithm.}
\item{iteration.stats}{a data frame of \code{iter} rows giving iteration statistics, as in \code{stats}, at each iteration of the algorithm.}
\item{topology}{the topology of the model.}
\item{family}{the vector of names of the distribution families used in the model. See \code{\link{LC_FAMILY}}.}
\item{distn}{the vector of names of the actual distributions used in the model. See \code{\link{LC_FAMILY}}.}
\item{iter.max}{the maximum number of distributions allowed in model fitting.}
\item{dname}{the name of the data.}
\item{dattr}{attributes of the data, used by model likelihood functions to determine if the data have been scaled or otherwise transformed.}
\item{zvec}{the vector of names of \code{X}; if the elements of \code{X} are unnamed, names are assigned.}
\item{kvec}{a list of which the \eqn{z}{z}th element is a vector of integers from 1 to \eqn{K_z}{K_z}.}
\item{k0vec}{a vector of integers from 1 to \eqn{K_0}{K_0}.}
\item{prior}{the value of the \code{prior} parameter used in model fitting.  See Arguments.}
\item{marginals}{if \code{prefit} is \code{TRUE}, the marginal fits to the data, otherwise \code{NULL}.}
}

\references{
McLachlan, G.J. and Thriyambakam, K.  (2008)  \emph{The EM Algorithm and Extensions}, John Wiley & Sons.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{LC_FAMILY}} for distributions and families; \code{\link{mixmod}} for fitting single-data mixture models; \code{\link{reporting}} and \code{\link{likelihood}} for model reporting; \code{\link{rocinfo}} for performance evaluation; \code{convergencePlot} for behavior of the algorithm; \code{\link{simulation}} for simulating from the parameters of a model.
}

\examples{\dontrun{
data(CiData)
data(CiGene)
fit <- mdmixmod(CiData, c(2,3,2), topology="chained",
           family=c("pvii", "norm", "pvii"))
fit
# Chained (PVII, normal, PVII) mixture model ('pvii', 'mvnorm', 'pvii')
# Data 'CiData' of size 10244-by-(1,4,1) fitted to 2 (2,3,2) components
# Model statistics:
#       iter       llik       qval        bic     iclbic 
#     377.00  -75859.81  -87065.28 -152310.62 -174721.56 
margs <- marginals(fit)
allFits <- c(list(chained=fit), margs)
plot(multiroc(allFits, CiGene$target))
}}

\keyword{graphs}
\keyword{models}
\keyword{cluster}
\keyword{classif}
