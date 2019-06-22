### visible_utilities.R:  visible utility functions for the lcmix package

# Given a matrix x, return x divided by its sum (allDiv) or with its columns divided by its column sums (colDiv) or rows divided by its row sums (rowDiv).  These are common operations and these functions can be particularly when x is being calculated on the fly.  E.g. instead of writing "x = foo(...) ; x/rowSums(x)" you can simply write "rowDiv(foo(...))".  If log is TRUE, the logarithm of the answer will be returned; this is useful for precision, e.g.,  "exp(rowDiv(x, TRUE))" may be more precise than "rowDiv(x)".  TO DO:  look for places where these can be used to replace existing code, especially rowDiv in place of rowSums.
allDiv <- function(x, log=FALSE)
{
	if(log)
		log(x) - log(sum(x))
	else x / sum(x)
}
colDiv <- function(x, log=FALSE)
{
	if(log) {
		t(t(log(x)) - log(colSums(x)))
	} else t(t(x) / colSums(x))
}
rowDiv <- function(x, log=FALSE)
{
	if(log)
		log(x) - log(rowSums(x))
	else x / rowSums(x)
}

# Given a list of objects "...", bind them in the sense that vectors are concatenated, while matrices and vectors are bound column-wise (default) or row-wise.  The type of the first object in "..." determines the operation applied, so be careful!

anyBind <- function(..., byrow=FALSE)
{
	if(isMDF(list(...)[[1]]))
		if(byrow)
			rbind(...)
		else
			cbind(...)
	else
		c(...)
}

# Given a square matrix x, mirror the lower triangular portion to the upper triangular portion (by default) to create a symmetric matrix.
asSymmetric <- function(x, keep.lower=TRUE)
{
	dims = dim(x)
	if(dims[[1]] != dims[[2]])
		stop("'x' must be a square matrix")
	
	idx = if(keep.lower) upper.tri(x) else lower.tri(x)
	x[idx] = t(x)[idx]
	
	return(x)
}

# Return a printable "collapsed" version of a vector.  If "space" is TRUE, then the elements will be set off by spaces.  If quote is TRUE, then single quotes will be added around the elements of the vector. If as.expr is TRUE, the string will be rendered as an expression which can be parsed with evpat(), e.g. collapseVec(1:3, as.expr=TRUE) = "c(1, 2, 3)".
collapseVec <- function(x, space=TRUE, quote=FALSE, as.expr=FALSE)
{	
	collapse = ifelse(space, ", ", ",")
	if(quote) x = sprintf("'%s'", x)
	
	retn = sprintf("(%s)", paste(x, collapse=collapse))
	if(as.expr) retn = paste("c", retn, sep="")
	return(retn)
}

# Given a matrix x, return a matrix consisting of a sample of its columns or rows.  Arguments "size", "replace", and "prob" are as the arguments to sample().
colSample <- function(x, size, replace=FALSE, prob=NULL)
	x[, sample(1:ncol(x), size=size, replace=replace, prob=prob)]
rowSample <- function(x, size, replace=FALSE, prob=NULL)
	x[sample(1:nrow(x), size=size, replace=replace, prob=prob), ]

# Given a matrix x, return a matrix with the columns sorted by column sums or the rows sorted by row sums.  Arguments "decreasing" and "na.last" are as the arguments to sort().
colSort <- function(x, decreasing=FALSE, na.last=NA)
	x[, order(colSums(x), decreasing=decreasing, na.last=na.last)]
rowSort <- function(x, decreasing=FALSE, na.last=NA)
	x[order(rowSums(x), decreasing=decreasing, na.last=na.last), ]

# Return an N-length vector containing the [Manhattan | squared Euclidean] distances between all rows of the N-by-D matrix x and the D-length vector center.  This is by analogy to the built-in function mahalanobis() which returns the squared Mahalanobis distance.  The "..." argument exists solely to allow manhattan() and euclidean() to be called with additional named arguments, such as the "cov" and "inverted" arguments to mahalanobis(); these arguments will have no effect on the value returned by manhattan() or euclidean().	
euclidean <- function(x, center, ...)
{
	if(is.vector(x)) x = matrix(x, nrow=1)
	
	colSums((t(as.matrix(x)) - center)^2)
}
manhattan <- function(x, center, ...)
{
	if(is.vector(x)) x = matrix(x, nrow=1)
	
	colSums(abs(t(as.matrix(x)) - center))
}

# Evaluate parsable text, e.g., evpat("c(1, 2, 3)") + 1 = c(2, 3, 4)
evpat <- function(x) eval(parse(text=x))

# replace all values in x less than min with min, and all values greater than max with max
forceRange <- function(x, min=-Inf, max=Inf)
{
	x[x < min] = min
	x[x > max] = max
	
	return(x)
}

# Get indexed elements of a vector or rows of a matrix or data frame, or of these elements within a list
getIdx <- function(x, idx)
{
	if(is.list(x) && !is.data.frame(x))
		lapply(x, getIdx, idx=idx)
	else
		if(isMDF(x)) x[idx,] else x[idx]
}

# Return TRUE if "x" is a matrix or a data frame, otherwise FALSE
isMDF <- function(x) (is.matrix(x) || is.data.frame(x))

# Given an N-by-K matrix x, returns an N-by-K logical matrix of which the (n,k)th element is TRUE if x_{n,k} = [max|min](x_{n,1}, ..., x_{n,K}), otherwise FALSE.
isRowMax <- function(x)
{
	retn = (rowRanks(x) == ncol(x))
	dimnames(retn) = dimnames(x)
	return(retn)
}
isRowMin <- function(x) isRowMax(-x)

# Replacement for expand.grid so the last factor, rather than the first, varies fastest; x must be a list, or something convertible to a list.
listGrid <- function(x, KEEP.OUT.ATTRS=TRUE, stringsAsFactors=FALSE)
	rev(expand.grid(rev(as.list(x)), KEEP.OUT.ATTRS=KEEP.OUT.ATTRS, 
		stringsAsFactors=stringsAsFactors))

# Summary functions for numeric lists
listSum <- function(x) Reduce("+", x)
listMean <- function(x, w=NULL)
{
	if(is.null(w))
		listSum(x) / length(x)
	else if(!is.numeric(w))
		stop("'w' must be numeric or NULL")
	else {
		w = rep(w, length.out=length(x))
		listSum(mapply("*", w, x, SIMPLIFY=FALSE)) / sum(w)
	}
}

# Given a matrix x, returns the logarithm of absolute value of the determinant.  Like the built-in function det(), this is a convenience wrapper around determinant().
labsDet <- function(x)
{
	retn = determinant(x, logarithm=TRUE)$modulus
	attributes(retn) = NULL
	return(retn)
}

# Given an N-by-K matrix x, returns an N-by-K logical matrix of which the (n,k)th element is 1 if x_{n,k} = max(x_{n,1}, ..., x_{n,K}), otherwise 0.
matrix2design <- function(x) 1 * isRowMax(x)

# Return the row and column indices of the maximum (or minimum) element in a matrix; this is what it seems like which.max() (which.min()) ought to return for a matrix, but doesn't.  The as.vector() wrapping makes the function return a simple numeric vector rather than a one-row data frame.  Thanks to <http://r.789695.n4.nabble.com/returning-the-largest-element-in-an-array-matrix-td795214.html>
matrixMaxIdx <- function(x) as.vector(which(x == max(x), arr.ind=TRUE))
matrixMinIdx <- function(x) as.vector(which(x == min(x), arr.ind=TRUE))

# adjust p-values treating the entire matrix "x" as one set of p-values, then return the results with the appropriate shape and row/column names
matrixPadjust <- function(x, method=p.adjust.methods, n=length(x))
{
	if(!(is.matrix(x) && is.numeric(x)))
		stop("'x' must be a numeric matrix")
	
	retn = p.adjust(as.numeric(x), method=method, n=n)
	retn = matrix(retn, ncol=ncol(x))
	dimnames(retn) = dimnames(x)
	
	return(retn)
}

# Given a matrix x, return the rank of the matrix as determined by QR decomposition.  This is one of the methods available in Matrix::rankMatrix().
matrixRank <- function(x) qr(x, LAPACK=FALSE)$rank

# Given a square matrix x, return the trace of the matrix, that is, the sum of the elements along the diagonal.  This implementation is notably faster than the more obvious "sum(diag(x))".  Not called "trace" to avoid stepping on the built-in R function.  No checking is done, so it is the user's responsibility to ensure that x is square; the function will work with non-square matrices as well!
matrixTrace <- function(x) sum(x[row(x)==col(x)])

# Return the MLE of the covariance of x assuming that mu is the mean.  By default, if x is univariate, a scalar will be returned, but if simplify is FALSE, then a matrix will be returned regardless.
mleCov <- function(x, mean=colMeans(as.matrix(x)), simplify=TRUE)
{
	x = as.matrix(x)
	retn = .mleTcov(t(x), mean, nrow(x))
	if(simplify) drop(retn) else retn
}
mleVar <- mleCov

# Return the MLE of the standard deviation of x, or if x is a matrix or data frame, the vector of column standard deviations.
mleSd <- function(x, mean=colMeans(as.matrix(x)))
	drop(sqrt(diag(mleCov(x, mean, simplify=F))))

# Multivariate gamma and log-gamma functions.  The argument "x" is the scalar or vector value to be evaluated; "D" is the dimension of the function.
mgamma <- function(x, D) exp(lmgamma(x, D))
lmgamma <- function(x, D)
{
	if(length(x) > 1)
		D*(D-1)/4 * LC_LOGPI + rowSums(lgamma(outer(x, 0.5*(1-(1:D)), "+")))
	else
		D*(D-1)/4 * LC_LOGPI + sum(lgamma(x + 0.5*(1-(1:D))))
}
mdigamma <- function(x, D) 
{
	if(length(x) > 1)
		rowSums(digamma(outer(x, 0.5*(1-(1:D)), "+")))
	else
		sum(digamma(x + 0.5*(1-(1:D))))
}

# Construct a named list, data frame, or vector with the names of the objects being passed in; see <http://stackoverflow.com/questions/3057341/how-to-use-rs-ellipsis-feature-when-writing-your-own-function>.
namedList <- function(...)
{
	retn = list(...)
	names(retn) = as.list(substitute(list(...)))[-1L]
	
	return(retn)
}
namedDataFrame <- function(..., stringsAsFactors=FALSE) 
	as.data.frame(namedList(...), stringsAsFactors=stringsAsFactors)
namedVector <- function(...) unlist(namedList(...))

notification <- function(x, depth=0)
{ # Takes a message to print to the screen and a "depth" (number of leading 
  # tabs), and returns the result with the time and date.  Usually used with 
  # message(), e.g. message(time.notif(...)). Useful for keeping track of the 
  # progress of loops.
	
	lead = ifelse(depth, paste(rep("\t", depth), collapse=""), "")
	
	sprintf("%s%s @ %s", lead, x, format(Sys.time(), usetz=T))
		# I like the formatted Sys.time() better than date()
}

# extract numeric columns from data frame
numericols <- function(x, include.logical=TRUE)
{
	if(!is.data.frame(x)) stop("'x' must be a data frame")

	if(include.logical)
		x[sapply(x, function(col) (is.numeric(col) || is.logical(col)))]
	else
		x[sapply(x, is.numeric)]
}

# Given a vector x, scale the vector to lie within the range [0,1] or [a,b]
range01 <- function(x)
{
	xmax = max(x, na.rm=TRUE)
	xmin = min(x, na.rm=TRUE)
	
	if(xmin == xmax) # all elements of x are equal
		rep(0.5, length(x))
	else
		(x-xmin)/(xmax-xmin)
}
rangeab <- function(x, a, b)
{
	if(b <= a) stop("'b' must be strictly greater than 'a'")
	(range01(x) * (b-a)) + a
}

# Recursive list apply:  recursively apply a function 'f' to a list 'object' containing, at some level, objects of a class appearing in 'classes'.  This is a simpler alternative to 'base::rapply'.
rlapply <- function(object, f, classes="ANY", ...)
{
	if(is.list(object)) lapply(object, function(x)
		if("ANY" %in% classes)
			f(x, ...)
		else
			if(any(sapply(classes, function(cls) is(x, cls))))
				f(x, ...)
			else if(is.list(x))
				rlapply(x, f, classes, ...)
			else
				if(is.vector(x) && (length(x) > 1))
					rlapply(as.list(x), f, classes, ...)
				else
					NULL
	)
	else if(is.vector(object) && length(object) > 1)
		rlapply(as.list(object), f, classes, ...)
	else
		stop("'object' must be a list")
}

# Standardize data; if X is a vector, return the standardized vector; if X is a matrix or data frame, return with standardized columns  In the case of method="sd", this means dividing data by its standard deviation; in the case of method="mean", this means dividing data by its mean.  If mle is TRUE, then the maximum likelihood estimator of standard deviation rather than the unbiased estimator will be used.  (The MLE argument has no effect if method="mean".)  if center is TRUE, the data will also be centered to have mean (or column means) 0.  Thus "standardize(X, 'sd', FALSE, TRUE)" is equivalent to "scale(X)".  The return value will have attribute standardized, which is a list with elements $scale and $center, which have the same meanings as scaled:scale and scaled:center in the return value attributes of scale().
standardize <- function(X, method=c("sd", "mean"), mle=TRUE, center=FALSE)
{
	method = match.arg(method)
	
	asmat = is.matrix(X) || is.data.frame(X)
	if(asmat) tX = t(X) # we will use this repeatedly below
	
	mu = if(asmat) colMeans(X) else mean(X)
	if(method == "sd") sigma = if(mle) mleSd(X) else sd(X)
	
	if(center) {
		center = mu
		if(asmat) {
			tX = tX-center
			X = t(tX)
		} else X = X-center
	} else center = NULL
	
	scale = switch(method, sd=sigma, mean=mu)
	X = if(asmat) t(tX/scale) else X/scale
	
	retn = X
	attr(retn, "standardized") = list(scale=scale, center=center)
	return(retn)
}
standardise <- standardize

# Return k evenly spaced quantiles of X (or quantiles of the columns of X, if X is a matrix or data frame) which are by default ordered from greatest to least.
uniformQuantiles <- function(X, K, decreasing=TRUE)
{
	if(isMDF(X)) {
		D = ncol(X)
		decreasing = rep(decreasing, length.out=D)
			# allow for vector of "decreasing" values that vary by column of X
		retn = sapply(1:ncol(X), function(d)
			uniformQuantiles(X[,d], K, decreasing[[d]]))
		colnames(retn) = colnames(X)
		return(retn)
	} else {
		# probs = (if(decreasing) K:1 else 1:K) / (K+1)
			# unbiased uniform quantiles
		# probs = ((if(decreasing) K:1 else 1:K) - 0.5) / K
			# equal-area uniform quantiles
		# probs = ((if(decreasing) K:1 else 1:K) - 1) / (K-1)
			# uniform quantiles out to [0,1]
		if(K == 1) median(X) else {
			probs = ((if(decreasing) K:1 else 1:K) - 0.5) / K
			retn = quantile(X, probs, names=FALSE)
			if(var(retn) == 0) { # extreme distribution, so go to extremes
				probs = ((if(decreasing) K:1 else 1:K) - 1) / (K-1)
				quantile(X, probs, names=FALSE)
			} else return(retn)
		}
	}
}

# Given a string, change the first letter in the string (if any) to uppercase
upperFirst <- function(x)
{
	pos = regexpr("[A-za-z]", x)
		# check for upper-case already there so as not to capitalize the next letter (if any) after the first upper-case one
	if(pos > 0)
		substr(x, pos, pos) = toupper(substr(x, pos, pos))
	
	return(x)
}

### SANDBOX FUNCTIONS (TO DO:  we need these, so don't delete them!  But they can probably be tightened up quite a bit, and then documented and made visible.

.outerMin <- function(x,y)
{
	K1 = length(x)
	K2 = length(y)
	
	retn = matrix(nrow=K1, ncol=K2)
	for(k1 in 1:K1) for(k2 in 1:K2)
		retn[k1,k2] = min(x[k1], y[k2])
	
	return(retn)
}

.outerMax <- function(x,y)
{
	K1 = length(x)
	K2 = length(y)
	
	retn = matrix(nrow=K1, ncol=K2)
	for(k1 in 1:K1) for(k2 in 1:K2)
		retn[k1,k2] = max(x[k1], y[k2])
	
	return(retn)
}

.cdf2pdf <- function(X)
{
	if(is.matrix(X)) {
		K1 = nrow(X)
		K2 = ncol(X)
		if((K1 < 2) || (K2 < 2))
			stop("'X' must be at least a 2x2 matrix")
		X[2:K1,] = X[2:K1,] - X[1:(K1-1),]
		X[,2:K2] = X[,2:K2] - X[,1:(K2-1)]
	} else {
		K = length(X)
		if(K < 2)
			stop("'X' must be a vector of at least length 2")
		X[2:K] = X[2:K] - X[1:(K-1)]
	}
	
	return(X)
}

.discreteFrechetCopula <- function(p1, p2)
{
	P1 = cumsum(p1)
	P2 = cumsum(p2)
	
	.outerMin(P1, P2)
}

.discreteIndependenceCopula <- function(p1, p2)
{
	P1 = cumsum(p1)
	P2 = cumsum(p2)
	
	outer(P1, P2)
}

.dependentDiscretePdf <- function(p1, p2, dependence=1)
{
	if((dependence < 0) || (dependence > 1))
		stop("'dependence' must be between 0 and 1 (inclusive)")

	Ci = .discreteIndependenceCopula(p1, p2)
	Cf = .discreteFrechetCopula(p1, p2) # Frechet upper bound
	
	Cd = dependence*Cf + (1-dependence)*Ci
	
	.cdf2pdf(Cd)
}

.dependentDiscreteCpdf <- function(p1, p2, dependence=1)
	.dependentDiscretePdf(p1, p2, dependence) / p1
