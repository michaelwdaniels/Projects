### model_reporting.R:  model reporting functions for the lcmix package.  A number of functions here are set as methods using the function setMethodS3(), which requires the package "R.methodsS3".

### GENERAL MODEL EXTRACTION FUNCTIONS

# get component assignments from model
setMethodS3("assignment", "mixmod",
function(x, ...) apply(posterior(x), 1, which.max))

# get design matrix from model
setMethodS3("design", "mixmod",
function(x, ...) diag(x$K)[assignment(x),])

# get number of iterations from model
setMethodS3("iter", "mixmod",
function(x, ...) x$stats[["iter"]])

# get parameters from model
setMethodS3("params", "mixmod",
function(x, ...) x$params)

# get posterior probabilities from model
setMethodS3("posterior", "mixmod",
function(x, ...) x$posterior)

# get posterior probabilities from model
setMethodS3("stats", "mixmod",
function(x, ...) x$stats)

### PREDICTION METHODS

setMethodS3("predict", "mixmod",
function(x, newdata=NULL, use.train=FALSE,
	train=(if(use.train) x$train else NULL), weights.list=FALSE, ...)
{
	X = if(is.null(newdata)) x$X else newdata
	params = x$params
	
	if(use.train && is.null(x$train))
		stop("'x' has no training data, 'use.train=TRUE' invalid")
	
	if(is(x, "mdmixmod")) {
		
		estep <- switch(x$topology,
			layered = .estep_mdmixmod_layered(),
			chained = .estep_mdmixmod_chained())
		pstep <- switch(x$topology,
			layered = .pstep_mdmixmod_layered(x$distn, train),
			chained = .pstep_mdmixmod_chained(x$distn, train))
	
#*#		# <experimental>
#*#		
#*#		if(TRUE) if((!use.train) && (!is.null(x$train))) {
#*#		
#*#			Unt = posterior(x)[!x$train$istrain,]
#*#			hid = params$hidden
#*#			hid$prob0 = colMeans(Unt)
#*#			hid$lprob0 = log(hid$prob0)
#*#			
#*#			if(x$topology == "layered") {
#*#				
#*#				hid$prob = lapply(x$zvec, function(z)
#*#					drop(hid$prob0 %*% hid$cprob[[z]]))
#*#				hid$lprob = lapply(hid$prob, log)
#*#			
#*#			} else if(x$topology == "chained") {
#*#				
#*#				; # TO DO
#*#			
#*#			} else stop(sprintf("unrecognized topology '%s'", x$topology))
#*#			
#*#			params$hidden = hid
#*#		}
#*#		
#*#		# </experimental>
		
		pdfs = pstep(X, params, x$kvec, x$zvec)
		weights = estep(pdfs, params, x$zvec, x$kvec, x$k0vec)
		
		if(weights.list) weights else weights$U
	
	} else {
	
		estep <- .estep_mixmod()
		pstep <- .pstep_mixmod(x$distn, train)
		
#*#		# <experimental>
#*#		
#*#		if(FALSE) if((!use.train) && (!is.null(x$train))) {
#*#			
#*#			; # TO DO
#*#			
#*#		}
#*#		
#*#		# </experimental>
		
		pdfs = pstep(X, params, x$kvec)
		weights = estep(pdfs)
		
		if(weights.list) weights else weights$W
	}
})

### MULTIPLE-DATA MODEL EXTRACTION FUNCTIONS

# Get marginals for object of class "mdixmod".
setMethodS3("marginals", "mdmixmod",
function(x, ...)
{
	if(is.null(x$marginals)) {
		lapply(x$zvec, function(z) {
			res = mixmod(x$X[[z]], x$K[[z]], x$family[[z]], dname=x$zvec[[z]], 
				iter.max=x$iter.max)
			return(res)
		})
	} else x$marginals
})

# Get full and marginal posteriors for object of class "mdmixmod".
setMethodS3("posteriors", "mdmixmod",
function(x, ...) c(list(JOINT=posterior(x)), lapply(marginals(x), posterior)))

# Get full and marginal predictions for object of class "mdmixmod".
setMethodS3("predictions", "mdmixmod",
function(x, newdata=NULL, use.train=FALSE,
	train=(if(use.train) x$train else NULL), weights.list=FALSE, ...)
{
	JOINT = predict(x, newdata=newdata, use.train=use.train, train=train, 
		weights.list=weights.list, ...)
	
	mx = marginals(x)
	margs = if(is.null(newdata))
		lapply(mx, predict, use.train=use.train,
			train=train, weights.list=weights.list, ...)
	else
		lapply(seq(length(mx)), function(z)
			predict(mx[[z]], newdata=newdata[[z]], use.train=use.train, 
				train=train, weights.list=weights.list, ...))
	names(margs) = names(mx) # in case of the "else" case above
	
	c(namedList(JOINT), margs)
})

### GENERAL LIKELIHOOD, ETC. FUNCTIONS
# TO DO:  there are a bunch of BIC, ABIC, deviance, log-likelihood, mutual information, etc. type model reporting functions sitting around in flydev_all_chrom/{include.R, mi.R, lcmix.OLD.R}.  What we need to do is decide on a final set of performance measures and include those here, without all the cruft.

# AIC (Akaike information criterion) for an object of class "mixmod" or "mdmixmod".
setMethodS3("aic", "mixmod",
function(x, ...) 2*logLik(x) - 2*x$npar)

# BIC (Bayes information criterion) for an object of class "mixmod" or "mdmixmod".
setMethodS3("bic", "mixmod",
function(x, ...) 2*logLik(x) - x$npar*log(x$N))

# Classification entropy for an object of class "mixmod" or "mdmixmod".
setMethodS3("entropy", "mixmod",
function(x, map=FALSE, ...) 2 * (logLik(x) - qval(x, map=map)))

# ICL-BIC (Bayes information criterion using the integrated complete likelihood) for an object of class "mixmod" or "mdmixmod".  If map == FALSE, the expected entropy method will be used; if TRUE, the MAP (maximum a posteriori) method will be used.
setMethodS3("iclbic", "mixmod",
function(x, map=FALSE, ...) 2*qval(x, map=map) - x$npar*log(x$N))

# Log-likelihood for an object of class "mixmod" or "mdmixmod".  The argument "object" is used here rather than "x" for conformance to the R generic logLik() definition.
setMethodS3("logLik", "mixmod",
function(object, ...) sum(object$pdfs$lgamma))

# Q-function and Q-value for an object of class "mixmod" or "mdmixmod".  If map == FALSE, the expected entropy method will be used; if TRUE, the MAP (maximum a posteriori) method will be used.
setMethodS3("qfun", "mixmod",
function(x, map=FALSE, ...)
	.qfun_mixmod(x$weights, x$params, x$pdfs, x$train, map=map))
setMethodS3("qfun", "mdmixmod",
function(x, map=FALSE, ...)
	.qfun_mdmixmod(x$weights, x$params, x$pdfs, x$train, map=map))
setMethodS3("qval", "mixmod",
function(x, map=FALSE, ...) sum(qfun(x, map=map)))

# "Simple ICL-BIC" for an object of class "mdmixmod" (different from ICL-BIC in that it uses only the top-level classification entropy).
setMethodS3("siclbic", "mdmixmod",
function(x, map=FALSE, ...)
{
	entropy = -2 * ifelse(map,
		sum(isRowMax(x$posterior) * log(x$posterior)),
		sum(x$posterior * log(x$posterior)))
	
	bic(x) - entropy
		
})

### PRINTING FUNCTIONS

# Print a representation of an object of class "mdmixmod".
setMethodS3("print", "mdmixmod",
function(x, ...)
{
	headstrbase = "%s %s mixture model %s\n"
	family = collapseVec(sapply(x$family, .pretty_family))
	distn = collapseVec(x$distn, quote=TRUE)
	headstr = upperFirst(sprintf(headstrbase, x$topology, family, distn))
	cat(headstr)
	
	parstrbase = "Data '%s' of size %i-by-%s fitted to %i %s components\n"
	D = collapseVec(x$D, space=FALSE)
	K = collapseVec(x$K, space=FALSE)
	parstr = sprintf(parstrbase, x$dname, x$N, D, x$K0, K)
	cat(parstr)
	
	cat("Model statistics:\n")
	print(x$stats, ...)
})

# Print a representation of an object of class "mixmod".
setMethodS3("print", "mixmod",
function(x, ...)
{
	family = .pretty_family(x$family)
	
	headstrbase = "%s mixture model ('%s')\n"
	headstr = upperFirst(sprintf(headstrbase, family, x$distn))
	cat(headstr)
	
	parstrbase = "Data '%s' of %s fitted to %i components\n"
	moddim = ifelse(x$D > 1,
		sprintf("size %i-by-%i", x$N, x$D),
		sprintf("length %i", x$N))
	parstr = sprintf(parstrbase, x$dname, moddim, x$K)
	cat(parstr)
	
	cat("Model statistics:\n")
	print(x$stats, ...)
})

### PLOTTING FUNCTIONS

# Plot the convergence of the log-likelihood for an object of class "mixmod" or "mdmixmod".  If show.qval is TRUE, also show Q-value.
convergencePlot <- function(x, show.qval=FALSE, main="")
{
	# draw log-likelihood always
	ylim = range(x$iteration.stats$llik)
	ylab = "log-likelihood"
	plot(x$iteration.stats$llik, type="l", xlab="iteration", ylab=ylab,
		ylim=ylim, main=main, axes=FALSE)
	axis(1)
	axis(2, pretty(range(x$iteration.stats$llik)))
	
	# draw Q-value PRN
	if(show.qval) {
		ylim = range(x$iteration.stats$qval)
		ylab = "Q-value"
		par(new=TRUE)
		plot(x$iteration.stats$qval, col="red", type="l",
		     axes=FALSE, xlab="", ylab="")
		axis(2, pretty(range(x$iteration.stats$qval)),
			col.ticks="red", col.axis="red", padj=1)
		mtext("Q-value", 2, col="red", padj=-3.5)
	}
}

### ROC FUNCTIONS

# Calculate ROC curve coordinates, AUC, and DCA for list of probabilities and (optionally, a list of) labels.
multiroc <- function(x, labels, k=1, quasi=FALSE)
{
	# data munging and sanity check; the idea here is that x and labels should be lists of equal length, with the elements of x being numeric vectors and the elements of labels being Boolean (or integer 0/1) vectors
	if(is.matrix(x)) # most common non-list case
		x = as.list(as.data.frame(x))
	else if(!is.list(x)) { # fallback (e.g., single vector of numbers)
		predname = deparse(substitute(x))
		x = list(x)
		names(x) = predname
	}
	
	if(is.matrix(labels))
		labels = as.list(as.data.frame(labels))
	else if(!is.list(labels))
		labels = list(labels)
	if(length(labels) > length(x))
		stop("more labels than x")
	else if(length(labels) < length(x))
		labels = rep(labels, length.out=length(x))
	names(labels) = names(x)
	
	# get rocinfo for each set of x and labels
	retn = mapply(rocinfo, x, labels, k=k, quasi=quasi, SIMPLIFY=FALSE)
	
	# package it up, send it back
	class(retn) = "multiroc"
	return(retn)
}

setMethodS3("plot", "multiroc",
function(x, legend=names(x), cex.legend=1, auc=TRUE, dca=FALSE, bw=FALSE, 
	col=(if(bw) rep(1, length(x)) else 1:length(x)),
	lty=(if(bw) 1:length(x) else rep(1, length(x))),
	lwd=rep(1, length(x)), ylab="true positive rate",
	xlab=ifelse(x[[1]]$quasi, "all positive rate", "false positive rate"),
	grid=FALSE, gridres=0.1, gridcol="lightgray", gridlty="dashed", gridlwd=1, 
	cut=NULL, ...)
{
	# data munging; the idea here is that legend, col, lty, and lwd should all be vectors of the same length of x; note that legend may also be FALSE, in which case none of this happens
	xlen = length(x)
	if(is.character(legend)) {
		if(length(legend) < xlen) legend = rep(legend, length.out=xlen)
		if(length(col) < xlen) col = rep(col, length.out=xlen)
		if(length(lty) < xlen) lty = rep(lty, length.out=xlen)
		if(length(lwd) < xlen) lwd = rep(lwd, length.out=xlen)
	}
	
	# set up blank plot
	plot(0.5, 0.5, pch="", xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab, ...)
	
	# draw grid PRN
	if(grid) {
		gridseq = seq(0, 1, gridres)
		abline(h=gridseq, v=gridseq, col=gridcol, lty=gridlty, lwd=gridlwd)
	}
	
	# draw ROC curves
	for(i in 1:xlen) lines(x[[i]], col=col[[i]], lty=lty[[i]], lwd=lwd[[i]])
	
	# add cuts PRN
	if(!is.null(cut)) for(i in 1:xlen) for(ct in cut)
	{
		xcoord = if(x[[i]]$quasi) x[[i]]$apr else x[[i]]$fpr
		idx = which.min((x[[i]]$x-ct)^2)
		xidx = xcoord[idx]
		yidx = x[[i]]$tpr[idx]
		points(xidx, yidx, col=col[[i]], pch=15, cex=1)
		text(xidx, yidx+0.035, ct, col=col[[i]], cex=1)
	}
	
	# add legend PRN
	if(is.character(legend)) {
		if(auc) legend = paste(legend,
			sapply(x, function(xi) sprintf("AUC=%.3f", xi$auc)),
			sep=", ")
		if(dca) legend = paste(legend,
			sapply(x, function(xi) sprintf("DCA=%.3f", xi$dca)),
			sep=", ")
		legend("bottomright", legend=legend, lty=lty, col=col, lwd=lwd,
			cex=cex.legend)
	}

})

# Given a vector of probabilities and a vector of labels, calculate the ROC curve.  The argument "x" should be a vector of probabilities, or a matrix for which the first column contains the probability of assignment to component 1, and "labels" should be TRUE when the true component is 1, otherwise FALSE.  (See Cortes and Mohri (2004) "Confidence Intervals for the Area under the ROC Curve" for a different, Wilcox-test based method of calculation, but this method works well too, and allows quasi-ROC-curve AUC calculation.)
setMethodS3("rocauc", "default",
function(x, labels, k=1, quasi=FALSE, ...) rocinfo(x, labels, k, quasi)$auc)

# ROC AUC for mixture model probability of assignment to component 1.
setMethodS3("rocauc", "mixmod",
function(x, labels, k=1, quasi=FALSE, ...)
	rocauc(posterior(x)[,k], labels, k, quasi))

# Calculate ROC curve coordinates, AUC, and DCA for single set of probabilities and labels.
setMethodS3("rocinfo", "default",
function(x, labels, k=1, quasi=FALSE, ...)
{
	# data munging PRN
	if(isMDF(x)) x = as.vector(x[,k])
	
	# sanity check
	N = length(x)
	if(length(labels) != N) stop("x and labels have unequal lengths")
	
	# put labels in proper order
	predorder = order(x, decreasing=TRUE)
	isPos = labels[predorder]
	isNeg = !isPos
	
	# calculate true positive rate, false positive rate, all positive rate
	numPos = sum(isPos)
	tpr = cumsum(isPos) / sum(isPos)
	fpr = cumsum(isNeg) / sum(isNeg)
	apr = (1:N) / N
	
	# calculate AUC and DCA
	xcoord = if(quasi) apr else fpr
	auc = sum(tpr * diff(c(0, xcoord))) # "step function integration"
	dca = sqrt(min(euclidean(cbind(tpr, xcoord), c(1,0))))
	
	# package it up, send it back
	retn = namedList(tpr, fpr, apr, auc, dca, quasi)
	retn$x = x[predorder]
	class(retn) = "rocinfo"
	return(retn)
})

# ROC info for mixture model probability of assignment to component 1.
setMethodS3("rocinfo", "mixmod",
function(x, labels, k=1, quasi=FALSE, ...)
	rocinfo(posterior(x)[,k], labels, k, quasi))

setMethodS3("lines", "rocinfo",
function(x, col=1, lty=1, lwd=1, cut=NULL, ...)
{
	xcoord = if(x$quasi) x$apr else x$fpr
	lines(x$tpr ~ xcoord, type="l", col=col, lty=lty, lwd=lwd, ...)
	if(!is.null(cut)) for(ct in cut)
	{
		idx = which.min((x$x-ct)^2)
		xidx = xcoord[idx]
		yidx = x$tpr[idx]
		points(xidx, yidx, col=col, pch=15, cex=1)
		text(xidx, yidx+0.035, ct, col=col, cex=1)
	}
})

setMethodS3("plot", "rocinfo",
function(x, legend="x", cex.legend=1, auc=TRUE, dca=FALSE,
	col=1, lty=1, lwd=1, ylab="true positive rate",
	xlab=ifelse(x$quasi, "all positive rate", "false positive rate"),
	grid=FALSE, gridres=0.1, gridcol="lightgray", gridlty="dashed", gridlwd=1,
	cut=NULL, ...)
{
	# set up blank plot
	plot(0.5, 0.5, pch="", xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab, ...)
	
	# draw grid PRN
	if(grid) {
		gridseq = seq(0, 1, gridres)
		abline(h=gridseq, v=gridseq, col=gridcol, lty=gridlty, lwd=gridlwd)
	}

	# draw actual ROC curve
	xcoord = if(x$quasi) x$apr else x$fpr
	lines(x$tpr ~ xcoord, type="l", col=col, lty=lty, lwd=lwd)
	
	# draw cut points PRN
	if(!is.null(cut)) for(ct in cut)
	{
		idx = which.min((x$x-ct)^2)
		xidx = xcoord[idx]
		yidx = x$tpr[idx]
		points(xidx, yidx, col=col, pch=15, cex=1)
		text(xidx, yidx+0.035, ct, col=col, cex=1)
	}
	
	# add legend PRN
	if(is.character(legend)) {
		if(auc) legend = paste(legend, sprintf("AUC=%.3f", x$auc), sep=", ")
		if(dca) legend = paste(legend, sprintf("DCA=%.3f", x$dca), sep=", ")
		legend("bottomright", legend=legend, lty=lty, col=col, lwd=lwd,
			cex=cex.legend)
	}
})

# Means for lists of rocinfo and multiroc objects
rocinfoListMean <- function(x)
{
	calcFun <- function(name)
	{
		meas = sapply(x, function(y) y[[name]])
		
		if(is.matrix(meas)) {
			mean = rowMeans(meas)
			sd = apply(meas, 1, sd)
			se = sd / sqrt(ncol(meas))
		} else {
			mean = mean(meas)
			sd = sd(meas)
			se = sd / sqrt(length(meas))
		}
		
		namedList(mean, se, sd)
	}
	
	measRes = lapply(namedList("tpr", "fpr", "apr", "auc", "dca", "x"), calcFun)
	retn = lapply(measRes, function(mr) mr$mean)
	retn$quasi = x[[1]]$quasi
	class(retn) = "rocinfo"
	
	attr(retn, "sd") = lapply(measRes, function(mr) mr$sd)
	attr(retn, "se") = lapply(measRes, function(mr) mr$se)
	
	return(retn)
}
multirocListMean <- function(x)
{
	I = length(x)
	J = length(x[[1]])
	
	retn = lapply(1:J, function(j)
		rocinfoListMean(lapply(1:I, function(i) x[[i]][[j]])))
	names(retn) = names(x[[1]])
	class(retn) = "multiroc"
	return(retn)
}

### OTHER PERFORMANCE MEASURES

# Given a matrix or vector of probabilities and a vector of labels, report on various measures related to true and false positives and negatives, relative to the kth model component.  The argument "prob" should be a matrix for which the first column contains the probability of assignment to component 1, and "labels" should be TRUE when the true component is 1, otherwise FALSE.  TO DO:  update documentation to show that "prob" can be a vector (in the range 0 to 1, otherwise it makes no sense) if k=1, or a matrix for any value of k.
setMethodS3("perfmeas", "default",
function(prob, labels, k=1, lfdr=NULL, midrange=FALSE, quasi=FALSE, 
	aslist=FALSE, ...)
{
	# munge probabilities as necessary
	if(!isMDF(prob))
	{	
		if(k != 1) stop("'k' must be equal to 1 for vector 'prob'")
		prob = cbind(prob, 1-prob)
	}
	
	# get predictions from probabilities
	pred = if(!is.null(lfdr))
		(prob[,k] >= 1-lfdr)
	else if(midrange)
		(prob[,k] >= sum(range(prob[,k]))/2)
	else
		(apply(prob, 1, which.max) == k)
	
	# total number of predictions and labels
	N = length(pred)
	if(N != length(labels))
		stop("prediction and label lengths do not match")
	
	# prediction positives and negatives
	predpos = sum(pred)
	predneg = sum(!pred)
	
	# label positives and negatives
	labelpos = sum(labels)
	labelneg = sum(!labels)
	
	# contingency table
	truepos = sum(pred & labels)
	falsepos = sum(pred & !labels)
	trueneg = sum(!pred & !labels)
	falseneg = sum(!pred & labels)
	
	# 0-1 performance measures
	perffun <- function(numerator, denominator)
		ifelse((numerator==0) && (denominator==0), 1, 
			exp(log(numerator) - log(denominator)))
	oddsr = perffun(truepos*trueneg, falsepos*falseneg)
	loddsr = log(oddsr)
	sensit = perffun(truepos, truepos+falseneg)
	specif = perffun(trueneg, trueneg+falsepos)
	precis = perffun(truepos, truepos+falsepos)
	recall = sensit
	f1meas = ifelse((precis == 0) && (recall == 0), 0,
		perffun(2*precis*recall, precis+recall))
	mattcc = perffun(truepos*trueneg - falsepos*falseneg,
		exp(0.5 * (log(truepos+falseneg) + log(truepos+falsepos)
			 + log(trueneg+falsepos) + log(trueneg+falseneg))))
			 	# Matthews correlation coefficient, as recommended by Baldi et al. (2000), "Assessing the accuracy of prediction algorithms for classification: an overview", Bioinformatics 16(5):412-424, doi 10.1093/bioinformatics/16.5.412
	tfdr = 1 - precis
	
	# other performance measures
	rocauc = rocauc(prob, labels, k, quasi)
	
	# package it up
	retn = namedList(N, predpos, predneg, labelpos, labelneg, truepos, 	
		falsepos, trueneg, falseneg, oddsr, loddsr, sensit, specif, precis, recall, f1meas, mattcc, tfdr, rocauc)
	
	# send it back
	if(aslist) retn else unlist(retn)
})

# Given a mixture model and a vector of labels, report on various measures related to true and false positives and negatives, relative to the first model component.  The argument "labels" should be TRUE when the true component is 1, otherwise FALSE.  "..." are further arguments to predict().  TO DO:  update documentation to allow for passing arguments to predict().
setMethodS3("perfmeas", "mixmod",
function(mod, labels, k=1, fdr=NULL, quasi=FALSE, aslist=FALSE, ...)
	perfmeas(predict(mod, ...), labels, k, fdr, quasi, aslist))

### SUB-MIXTURES (MIXTURES OF DATA IN MIXTURE COMPONENTS)

setMethodS3("submix", "mixmod",
function(mod, assignment=c("map", "lfdr"), assigncut=0.1, ...)
{
	# get assignments, data, and training
	isassign = switch(match.arg(assignment),
		map = {
			lapply(mod$kvec, function(k) (mod$assignment==k))
		},
		lfdr = {
			pcut = 1-assigncut
			post = predict(mod)
			lapply(mod$kvec, function(k) 
				(post[,k] >= pcut))
		}
	)
	X = lapply(isassign, function(i)
			if(sum(i)) getIdx(mod$X, i) else NULL)
	if(is.null(mod$train)) {
		train = train.weight = NULL
	} else {
		train = lapply(isassign, function(i)
			if(sum(i)) getIdx(mod$train$train, i) else NULL)
		train.weight = mod$train$weight
	}
			
	
	X = switch(match.arg(assignment),
		map = lapply(mod$kvec, function(k)
			getIdx(mod$X, mod$assignment==k)),
		lfdr = {
			pcut = 1-assigncut
			post = predict(mod)
			lapply(mod$kvec, function(k)
			{
				isk = (post[,k] >= pcut)
				nisk = sum(isk)
				if(nisk) getAssignment(mod$X, isk) else NULL
			})
		}
	)
	
	# extract additional arguments to mixmod()
	dotargs = list(...)
	dotargnames = names(dotargs)
	for(argname in c("K", "family"))
		assign(argname,
			replicate(mod$K,
			          if(argname %in% dotargnames)
			          	dotargs[[argname]]
			          else
			          	mod[[argname]],
		              simplify=FALSE))
})