### simulation.R:  simulation and resampling functions for the lcmix package.

# resample from data to create a data set with the same size, distribution, and hidden parameters as the fitted model "x", an object of class "mixmod" or "mdmixmod"
setMethodS3("resampleFromFit", "mixmod",
function(x, n=x$N, replace=TRUE, hidpar=x$params$hidden, ...)
{
	# simulate hidden data
	Y = .simulate_hidden_mixdata(n, hidpar)
	
	# simulate observed data
	x$X = as.matrix(x$X)
	X = matrix(nrow=n, ncol=ncol(x$X))
	for(k in x$kvec)
	{
		isk = (Y == k)
		nisk = sum(isk)
		prob = x$weights$W[,k]
		X[isk,] = rowSample(x$X, nisk, replace, prob)
	}
	X = drop(X)
	
	# package it up, send it back
	namedList(X, Y)
})
setMethodS3("resampleFromFit", "mdmixmod",
function(x, n=x$N, replace=TRUE, hidpar=x$params$hidden,
	topology=x$topology, ...)
{
	# simulate hidden data
	hidden = .simulate_hidden_mdmixdata(n, hidpar, topology)
	Y = hidden$Y
	Y0 = hidden$Y0
	
	# simulate observed data
	X = lapply(x$zvec, function(z)
	{
		xXz = as.matrix(x$X[[z]])
		Xz = matrix(nrow=n, ncol=ncol(xXz))
		for(kz in x$kvec[[z]])
		{
			iskz = (Y[[z]] == kz)
			niskz = sum(iskz)
			prob = x$weights$W[[z]][,kz]
			Xz[iskz,] = rowSample(xXz, niskz, replace, prob)
		}
		drop(Xz)
	})
	
	# package it up, send it back
	namedList(X, Y, Y0)
})

# simulate data with similar performance to the Ci data
simulateCiTypeData <- function(n=10000, par=LC_SIMPAR, topology=LC_TOPOLOGY)
{
	topology = match.arg(topology)
	
	if(topology == "layered") {
		par$hidden$cprob = lapply(par$hidden$prob, function(p)
			.dependentDiscreteCpdf(par$hidden$prob0, p, dependence=0.60))
	} else { # topology == "chained"
		par$hidden$cprob = list(
			binding = .dependentDiscreteCpdf(par$hidden$prob0,
				par$hidden$prob$binding, dependence=0.85),
			expression = .dependentDiscreteCpdf(par$hidden$prob$binding,
				par$hidden$prob$expression, dependence=0.75),
			conservation = .dependentDiscreteCpdf(par$hidden$prob$expression,
				par$hidden$prob$conservation, dependence=0.75))
	}
	
	simulateMdmixdata(n, c("norm", "mvnorm", "norm"), par, topology)
}

# simulate a data set with the same size, distribution, and parameters as the fitted model "x", an object of class "mixmod" or "mdmixmod"
setMethodS3("simulateFromFit", "mixmod",
function(x, n=x$N, ...)
	simulateMixdata(n, x$distn, x$params))
setMethodS3("simulateFromFit", "mdmixmod",
function(x, n=x$N, ...)
	simulateMdmixdata(n, x$distn, x$params, x$topology))

	
# Simulate multiple mixture data of the given size, distributions, parameters, and topology; params should be as the $params element of the return value of mdmixmod().  Return a list with elements $X (a list of observed data), $Y (a list of hidden components corresponding to the observed data), and $Y0 (top-level hidden components), with attributes @topology (a string giving the topology, i.e. "layered" or "chained") and @params (the params used in the simulation.)
simulateMdmixdata <- function(n, distn, params, topology=LC_TOPOLOGY)
{
	## gather statistics, munge input PRN
		
	topology = match.arg(topology)
	
	Z = length(params$hidden$cprob)
	zvec = names(params$hidden$cprob)
	if(is.null(zvec))
		zvec = paste("X", 1:Z, sep="")
	names(params$hidden$cprob) = names(params$observed) = names(distn) = zvec
		# just to make sure everything is using the same names
	names(zvec) = zvec
		# so results of lapply-type functions have proper names
	
	D = sapply(zvec, function(z) switch(distn[[z]],
			norm     = 1,
			mvnorm   = ncol(params$observed[[z]]$mean),
			weisd    = 1,
			mvweisd  = ncol(params$observed[[z]]$shape),
			nbinom   = 1,
			mvnbinom = ncol(params$observed[[z]]$size),
			npois    = 1,
			mvpois   = stop("multivariate Poisson not supported; use negative binomial (family 'nbinom') instead"),
			bern     = 1,
			mvbern   = ncol(params$observed[[z]]$prob),
			gamma    = 1,
			mvgamma  = ncol(params$observed[[z]]$shape),
			pvii     = 1,
			mvpvii   = ncol(params$observed[[z]]$mean)))
			
	
	K0 = length(params$hidden$prob0)
	K = sapply(params$hidden$cprob, ncol)
	
	## simulate hidden data
	
	hidden = .simulate_hidden_mdmixdata(n, params$hidden, topology)
	Y = hidden$Y
	Y0 = hidden$Y0
	
	## simulate observed data
	
	X = lapply(zvec, function(z)
	{
		par = params$observed[[z]]
		
		if(D[[z]] == 1) {
			res = vector(mode="numeric", length=n)
		} else {
			res = matrix(0, nrow=n, ncol=D[[z]])
			if(distn[[z]] %in% c("mvnorm", "mvpvii"))
				colnames(res) = colnames(params$observed$mean)
			else if(distn[[z]] %in% c("mvweisd", "mvgamma"))
				colnames(res) = colnames(params$observed$shape)
			else if(distn[[z]] == "mvnbinom")
				colnames(res) = colnames(params$observed$size)
			else if(distn[[z]] == "mvpois")
				stop("multivariate Poisson not supported; use negative binomial (family 'nbinom') instead")
			else if(distn[[x]] == "mvbern")
				colnames(res) = colnames(params$observed$prob)
		}
		
		for(kz in 1:K[[z]])
		{
			iskz = (Y[[z]] == kz)
			niskz = sum(iskz)
			
			if(D[[z]] == 1) res[iskz] = switch(distn[[z]],
				norm = rnorm(niskz, par$mean[kz], par$sd[kz]),
				weisd = rweisd(niskz, par$shape[kz], par$decay[kz]),
				gamma = rgamma(niskz, par$shape[kz], par$rate[kz]),
				pvii = rpvii(niskz, par$mean[kz], par$scale[kz], 
					par$shape[kz]),
				weisd = rweisd(niskz, par$size[kz], par$prob[kz]),
				stop(sprintf("unknown distribution '%s'", distn)))
			else res[iskz,] = switch(distn[[z]],
				mvnorm = rmvnorm(niskz, par$mean[kz,], par$cov[[kz]]),
				mvweisd = rmvweisd(niskz, par$shape[kz,], par$decay[kz,], 
					par$corr[[kz]]),
				mvgamma = rmvgamma(niskz, par$shape[kz,], par$rate[kz,], 
					par$corr[[kz]]), # mvpvii, mvnbinom, mvpois, mvbern TO DO
				stop(sprintf("unknown distribution '%s'", distn)))
		}
		
		return(res)
	})
	
	## package and return
	
	retn = list(X=X, Y=Y, Y0=Y0)
	attr(retn, "n") = n
	attr(retn, "distn") = distn
	attr(retn, "params") = params
	attr(retn, "topology") = topology
	
	return(retn)
}
			
# Simulate mixture data of the given data size, distribution, and parameters; params should be as the $params element of the return value of mixmod(). Return a list with elements $X (observed data) and $Y (hidden component) and attribute @params (the parameters used in the simulation).
simulateMixdata <- function(n, distn, params)
{
	## simulate hidden data
	
	Y = .simulate_hidden_mixdata(n, params$hidden)
	
	## simulate observed data
	
	# gather statistics and set up observed data
	
	D = switch(distn,
			norm     = 1,
			mvnorm   = ncol(params$observed$mean),
			weisd    = 1,
			mvweisd  = ncol(params$observed$shape),
			gamma    = 1,
			mvgamma  = ncol(params$observed$shape),
			pvii     = 1,
			mvpvii   = ncol(params$observed$mean),
			nbinom   = 1,
			mvnbinom = ncol(params$observed$lambda),
			pois     = 1,
			mvpois   = stop("multivariate Poisson not supported; use negative binomial (family 'nbinom') instead"),
			bern     = 1,
			mvbern   = ncol(params$observed$prob),
			stop(sprintf("distn '%s' not yet supported", distn)))
	
	if(D == 1) {
		X = vector(mode="numeric", length=n)
	} else {
		X = matrix(0, nrow=n, ncol=D)
		if(distn %in% c("mvnorm", "mvpvii"))
			colnames(X) = colnames(params$observed$mean)
		else if(distn %in% c("mvweisd", "mvgamma"))
			colnames(X) = colnames(params$observed$shape)
		else if(distn == "mvnbinom")
			colnames(X) = colnames(params$observed$size)
		else if(distn == "mvpois")
			colnames(X) = stop("multivariate Poisson not supported; use negative binomial (family 'nbinom') instead")
		else if(distn == "mvbern")
			colnames(X) = colnames(params$observed$prob)
	}
	
	# simulate
	
	par = params$observed
	
	for(k in 1:length(params$hidden$prob))
	{
		isk = (Y == k)
		nisk = sum(isk)
		
		if(D == 1) X[isk] = { # TO DO:  replace these chains of "if-else-if" with "switch" statements a la simulateMdmixdata().
			if(distn == "norm")
				rnorm(nisk, par$mean[k], par$sd[k])
			else if(distn == "weisd")
				rweisd(nisk, par$shape[k], par$decay[k])
			else if(distn == "gamma")
				rgamma(nisk, par$shape[k], par$rate[k])
			else if(distn == "pvii")
				rpvii(nisk, par$mean[k], par$scale[k], par$shape[k])
			else if(distn == "nbinom")
				rnbinom(nisk, par$size[k], par$prob[k])
			else if(distn == "pois")
				rpois(nisk, par$lambda[k])
			else if(distn == "bern")
				rbern(nisk, par$prob[k])
		} else X[isk,] = {
			if(distn == "mvnorm")
				rmvnorm(nisk, par$mean[k,], par$cov[[k]])
			else if(distn == "mvweisd")
				rmvweisd(nisk, par$shape[k,], par$decay[k,], par$corr[[k]])
			else if(distn == "mvgamma")
				rmvgamma(nisk, par$shape[k,], par$rate[k,], par$corr[[k]])
			else if(distn == "mvpvii")
				rmvpvii(nisk, par$mean[k,], par$scale[[k]], par$shape[k])
			else if(distn == "mvnbinom") # TO DO
				stop("multivariate negative binomial not yet supported")
			else if(distn == "mvpois")
				stop("multivariate Poisson not supported; use negative binomial (family 'nbinom') instead")
			else if(distn == "mvbern")
				rmvbern(nisk, par$prob[k,], par$corr[[k]])
		}
	}
				
	## package and return
	
	retn = list(X=X, Y=Y)
	attr(retn, "n") = n
	attr(retn, "distn") = distn
	attr(retn, "params") = params
	return(retn)
}

# Cross-validation to determine appropriate training weight by maximizing (quasi-)ROC AUC; "..." are further arguments to mixmod() or mdmixmod().  The "applyfun" argument may be a drop-in replacement for lapply(), e.g. for parallelization.
train.xval <- function(X, K, train, k=1, 
	train.weights=c(0, 1, 5, 10, 20, 50, 100),
	type=c("mixmod", "mdmixmod"), quasi=TRUE, iter=30, train.prop=0.5, 
	applyfun=lapply, seed.coef=NULL, ...)
{
	type = match.arg(type)
	mixfun = get(type)
	
	N = switch(type,
		mixmod = ifelse(isMDF(X), nrow(X), length(X)),
		mdmixmod = ifelse(isMDF(X[[1]]), nrow(X[[1]]), length(X[[1]])))
	Nseq = seq(N)
	
	isLabeled = (train == k)
	numLabeled = sum(isLabeled)
	Ntrain = as.integer(round(numLabeled * train.prop))
	if((Ntrain == 0) || (numLabeled < (Ntrain+1)))
		stop("insufficient labeled examples for chosen training proportion")
	labeledIdx = which(isLabeled)
	
	unsup = (0 %in% train.weights)
	if(unsup) us = mixfun(X, K, ...)
		# we only need to do unsupervised once, so get it out of the way
	
	utw = sort(unique(train.weights))
	tws = setdiff(utw, 0)
	names(tws) = tws
	runs = do.call(rbind, applyfun(seq(iter), function(i)
	{
		if(!is.null(seed.coef)) set.seed(seed.coef * i)
		
		trainIdx = sample(labeledIdx, Ntrain)
		testIdx = setdiff(Nseq, trainIdx)
		trainLoc = train; trainLoc[testIdx] = 0
		
		fits = lapply(tws, function(tw)
			mixfun(X, K, train=trainLoc, train.weight=tw, ...))
		if(unsup) fits = c(list("0" = us), fits)
		
		labels = isLabeled[testIdx]
		sapply(fits, function(x) rocauc(predict(x)[testIdx, k], labels))
	}))
	
	means = colMeans(runs)
	
	retn = as.integer(utw[[which.max(means)]])
	attr(retn, "means") = means
	attr(retn, "runs") = runs
	class(retn) = c("train.xval", "integer")
	return(retn)
}
setMethodS3("print", "train.xval", function(x, ...) print(as.integer(x), ...))
setMethodS3("summary", "train.xval",
function(x, ...)
{
	retn = data.frame(mean=attr(x, "means"), sd=apply(attr(x, "runs"), 2, sd))
	retn$se = retn$sd / sqrt(nrow(attr(x, "runs")))
	return(retn)
})

# Another take on cross-validation, using the traditional "k-fold" method, except that we do _not_ call it that here since "K" and "k" already mean something different.
foldXval <- function(X, folds=10, K, train, k=1,
	train.weights=c(0, 1, 5, 10, 20, 50, 100), type=c("mixmod", "mdmixmod"), quasi=TRUE, applyfun=lapply, seed=NULL, ...)
{
	if(folds < 2)
		stop("number of folds must be at least 2")
	type = match.arg(type)
	mixfun = get(type)
	
	N = switch(type,
		mixmod = ifelse(isMDF(X), nrow(X), length(X)),
		mdmixmod = ifelse(isMDF(X[[1]]), nrow(X[[1]]), length(X[[1]])))
	Nseq = seq(N)
	
	isLabeled = (train == k)
	numLabeled = sum(isLabeled)
	if(folds > numLabeled)
		stop("number of folds greater than length of labeled data")
	labeledIdx = which(isLabeled)
	unlabeledIdx = which(!isLabeled)
	
	unsup = (0 %in% train.weights)
	if(unsup) us = mixfun(X, K, ...)
		# we only need to do unsupervised once, so get it out of the way
	
	if(!is.null(seed)) set.seed(seed)
	foldTestGroup = sample(rep(1:folds, length.out=N))
	foldTestLabeledIdx = split(labeledIdx, foldTestGroup) # MARK MARK
	
	utw = sort(unique(train.weights))
	tws = setdiff(utw, 0)
	names(tws) = tws
	runs = do.call(rbind, applyfun(1:folds, function(f)
	{
		if(!is.null(seed.coef)) set.seed(seed.coef * i)
		
		trainIdx = sample(labeledIdx, Ntrain)
		testIdx = setdiff(Nseq, trainIdx)
		trainLoc = train; trainLoc[testIdx] = 0
		
		fits = lapply(tws, function(tw)
			mixfun(X, K, train=trainLoc, train.weight=tw, ...))
		if(unsup) fits = c(list("0" = us), fits)
		
		labels = isLabeled[testIdx]
		sapply(fits, function(x) rocauc(predict(x)[testIdx, k], labels))
	}))
	
	means = colMeans(runs)
	
	retn = as.integer(utw[[which.max(means)]])
	attr(retn, "means") = means
	attr(retn, "runs") = runs
	class(retn) = c("train.xval", "integer")
	return(retn)
}
setMethodS3("print", "train.xval", function(x, ...) print(as.integer(x), ...))
setMethodS3("summary", "train.xval",
function(x, ...)
{
	retn = data.frame(mean=attr(x, "means"), sd=apply(attr(x, "runs"), 2, sd))
	retn$se = retn$sd / sqrt(nrow(attr(x, "runs")))
	return(retn)
})

