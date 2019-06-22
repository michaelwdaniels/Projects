\name{ROC curves}

\alias{rocinfo}
\alias{rocinfo.default}
\alias{rocinfo.mixmod}
\alias{plot.rocinfo}
\alias{lines.rocinfo}

\alias{rocauc}
\alias{rocauc.default}
\alias{rocauc.mixmod}

\alias{multiroc}
\alias{plot.multiroc}

\concept{ROC curve}

\title{
Receiver Operating Characteristic (ROC) Curves
}

\description{
Functions to calculate and plot Receiver Operating Characteristic (ROC) curves for one or more sets of predictions, and for one or more objects of class \code{mixmod} or \code{mdmixmod}.
}

\usage{
\method{rocinfo}{default}(x, labels, quasi=FALSE, ...)
\method{rocinfo}{mixmod}(x, labels, quasi=FALSE, ...)
\method{plot}{rocinfo}(x, legend="x", cex.legend=1, auc=TRUE, dca=TRUE,
    col=1, lty=1, lwd=1, ylab="true positive rate",
    xlab=ifelse(x$quasi, "all positive rate", "false positive rate"),
    grid=FALSE, gridres=0.1, gridcol="lightgray", gridlty="dashed", ...)
\method{rocauc}{default}(x, labels, quasi=FALSE, ...)
\method{rocauc}{mixmod}(x, labels, quasi=FALSE, ...)

multiroc(x, labels, quasi=FALSE)
\method{plot}{multiroc}(x, legend=names(x), cex.legend=1, auc=TRUE, dca=FALSE, bw=FALSE,
    col=(if(bw) rep(1, length(x)) else 1:length(x)),
    lty=(if(bw) 1:length(x) else rep(1, length(x))),
    lwd=rep(1, length(x)), ylab="true positive rate",
    xlab=ifelse(x[[1]]$quasi, "all positive rate", "false positive rate"),
    grid=FALSE, gridres=0.1, gridcol="lightgray", gridlty="dashed", ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{for \code{rocinfo} and \code{rocauc}, a vector of probabilities, or an object of class \code{mixmod} or \code{mdmixmod}; for \code{multiroc}, a list of objects suitable to be passed as arguments to \code{rocinfo}; for plotting methods, an object of the appropriate class.}
\item{labels}{a vector of logical values of the same length as the probabilities, in which a value of \code{TRUE} indicates a true positive; for \code{multiroc}, optionally a list of such vectors of the same length, with the same element names, as \code{x}.}
\item{quasi}{logical; if \code{TRUE}, the x-axis of the ROC curve denotes the total positive rate rather than the false positive rate.}
\item{legend}{character or \code{FALSE}; if character, a legend will be drawn with the name(s) in \code{legend}; if \code{FALSE}, no legend will be drawn.}
\item{cex.legend}{magnification to be used for legend text.}
\item{auc}{logical; if \code{TRUE}, AUC, or Area Under the Curve, will be shown in the legend.  Larger AUC indicates better performance.  No effect if \code{legend} is \code{FALSE}.}
\item{dca}{logical; if \code{TRUE}, DCA, or Distance of Closest Approach of the curve to (0,1), will be shown in the legend.  Smaller DCA indicates better performance.  No effect if \code{legend} is \code{FALSE}.}
\item{bw}{logical; if \code{TRUE}, a black-and-white version of the plot will be drawn.  Convenience function for drawing multiple ROC curves on the same plot.}
\item{col, lty, lwd, ylab, xlab}{as the equivalent arguments to \code{plot.default} or \code{lines.default}.  For \code{plot.multiroc}, \code{col}, \code{lty}, and \code{lwd} should be vectors of the same length as \code{x}.}
\item{grid}{logical; if \code{TRUE}, draw a grid on the plot.}
\item{gridres, gridcol, gridlty}{resolution, color, and line type for the grid.  No effect if \code{grid} is \code{FALSE}.}
\item{...}{further arguments to \code{plot.default} or \code{lines.default} for the plotting methods, no effect for the other functions.}
}

\details{
\code{rocinfo} calculates the ROC curve x coordinates (false positive rate, or all positive rate if \code{quasi} is \code{TRUE}), y coordinates (true positive rate), AUC, and DCA for the given probabilities and labels.  \code{multiroc} performs the same function for multiple probabilities, and optionally multiple labels.
}

\value{
For \code{rocinfo}, a list of class \code{rocinfo} having the following elements:
\item{tpr}{the true positive rate in ascending order of the probabilities.}
\item{fpr}{the false positive rate in ascending order of the probabilities.}
\item{apr}{the all positive rate in ascending order of the probabilities.}
\item{auc}{the Area Under the Curve.}
\item{dca}{the Distance of Closest Approach.}
\item{quasi}{logical; the same as the argument to \code{rocinfo}.}

For \code{rocauc}, the numerical value of the AUC.

For \code{multiroc}, a list of class \code{multiroc} having elements of class \code{rocinfo}.

The \code{plot} and \code{lines} methods are used for their side effects.
}

\author{
Daniel Dvorkin
}

\examples{\dontrun{

### multiple ROC curves for a single set of labels

data(CiData)
data(CiGene)
CiJoint <- lapply(namedList("layered", "chained"), function(top)
    mdmixmod(CiData, c(2,3,2), family=c("pvii", "norm", "pvii"),
    topology=top))
CiMarginal <- marginals(CiJoint$layered) # same for $layered and $chained
CiFits <- c(CiJoint, CiMarginal)
sapply(CiFits, rocauc, labels=CiGene$target)
#     layered      chained      binding   expression conservation 
#   0.8507258    0.8913765    0.8596395    0.7943286    0.7479479
plot(multiroc(CiFits, CiGene$target), grid=TRUE)

### multiple ROC curves for multiple labels

DlFits <- mdmixmod(DlData, c(2,2,2), family=c("norm", "norm", "pvii"))
AllFits = list(Ci=CiJoint$chained, Dl=DlFits)
plot(multiroc(AllFits, list(Ci=CiGene$target, Dl=DlGene$target)))
}}

\keyword{hplot}
