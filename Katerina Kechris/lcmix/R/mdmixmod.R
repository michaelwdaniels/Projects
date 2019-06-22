### mdmixmod.R:  define the mdmixmod (multiple data mixture models) function for the lcmix package.

mdmixmod <- function(X, K, K0=min(K), topology=LC_TOPOLOGY, 
	family=NULL, prior=NULL, train=NULL, train.weight=1, decreasing=TRUE, 
	prefit=TRUE, iter.max=LC_ITER_MAX, dname=deparse(substitute(X)))
{
	## munge data, gather preliminary information, and do sanity checks
		
	if(!is.list(X)) stop("'X' must be a list")
	
	dname = dname
		# keep the implicit call to deparse(substitute(X)) from getting screwed up by whatever manipulation we do below (this is a nasty side-effect of lazy evaluation, I suspect)
	dattr = lapply(X, attributes)
	
	X = lapply(X, function(Xz) drop(as.matrix(Xz)))
		# put all data into matrix or vector form
	
	Z = length(X)
	if(Z < 2) stop("single data source, use mixmod() instead")
	if(is.null(names(X)))
		names(X) = paste("X", 1:Z, sep="")
	zvec = names(X)
	Z = length(zvec)
	names(zvec) = zvec
		# so the results of lapply() and related functions will always have the names of the elements of X
	
	univariate = !sapply(X, is.matrix)
	ntmp = sapply(zvec, function(z)
		ifelse(univariate[[z]], length(X[[z]]), nrow(X[[z]])))
	if(min(ntmp) != max(ntmp))
		stop("unequal data sizes")
	else
		N = ntmp[[1]]
	
	D = sapply(zvec, function(z) ifelse(univariate[[z]], 1, ncol(X[[z]])))
	
	topology = match.arg(topology)
	
	if(is.null(family)) family = rep("normal", Z) else
	{
		family = match.arg(family, names(LC_FAMILY), several.ok=TRUE)
		if(length(family) > Z)
			warning("length(family) > Z, truncating to fit")
		if(length(family) != Z)
			family = rep(family, length.out=Z)
	}
	names(family) = zvec
	dummy = lapply(zvec, function(z) .sanity_check_data(X[[z]], family[[z]]))
	
	distn = sapply(zvec, function(z)
		ifelse(univariate[[z]], LC_FAMILY[[(family[[z]])]]$uni, 
			                    LC_FAMILY[[(family[[z]])]]$multi))
	
	if(length(K) > Z)
		warning("length(K) > Z, truncating to fit")
	if(length(K) != Z)
		K = rep(K, length.out=Z)
	names(K) = zvec
	kvec = lapply(zvec, function(z) 1:K[[z]])
	
	k0vec = 1:K0
	
	if(length(decreasing) > Z)
		warning("length(decreasing) > Z, truncating to fit")
	if(length(decreasing) != Z)
		decreasing = rep(decreasing, length.out=Z)
	names(decreasing) = zvec
	
	npar.hidden = switch(topology,
		layered = K0-1 + K0*sum(K-1),
		chained = K0-1 + sum(c(K0, K[1:(Z-1)]) * (K-1))
		#chained = sum(c(K0, K[1:(Z-1)]) * (K-1)) + K0-1 - ((K0-1)*K[[1]])
	)
	npar.observed = sum(sapply(zvec, function(z)
		.npar_observed(distn[[z]], K[[z]], D[[z]])))
	npar = npar.hidden + npar.observed
	
	## handle training PRN
	
	if(train.weight == 0) train = NULL
	if(train.weight < 0) stop("negative 'train.weight' is not valid")
	
	if(is.null(train)) {
		
		train.raw = NULL
		bicsub = npar * log(N)
	
	} else {
		
		# preserve old state for passing to mixmod() in initialization
		train.raw = train
		
		# sanity checks and component fixing
		if(length(train) != N)
			stop("'train' must be NULL or of same length as data")
		if(!all(train %in% (0:(K0+1))))
			stop("'train' values must be integers between 0 and K0+1 inclusive")
		train[(train==0)] = K0+1
		
		# expansion
		weight = train.weight
		Ktrain = K0+1
		Tprime = diag(Ktrain)[train,]
		istrain = (train <= K0)
		wTprime = Tprime
		wTprime[istrain,] = weight * wTprime[istrain,]
		ltprob = log(colSums(wTprime)) - log(sum(wTprime))
		tprob = exp(ltprob)
		rpbase = diag(K0)
		lrpbase = log(rpbase)
		Tpdf = t(t(Tprime) * tprob)
		
		# vector for multiplying labeled rows of weight matrices, effective data set size, and coefficient for dividing column means calculated from such training-weight-multiplied matrices to give proper probabilities
		twmult = rep.int(1, N) + (weight-1)*istrain
		Neff = sum(twmult)
		ltwcoef = log(Neff) - log(N)
		twcoef = exp(ltwcoef)
		
		# package it up
		train = namedList(train, weight, Ktrain, Tprime, istrain, tprob, ltprob, 
			rpbase, lrpbase, Tpdf, twmult, Neff, twcoef, ltwcoef)
		
		# additional calculations
		bicsub = npar * log(Neff)
	}
		
	## define internal functions
	
	estep <- switch(topology,
		layered = .estep_mdmixmod_layered(),
		chained = .estep_mdmixmod_chained())
	mstep <- switch(topology,
		layered = .mstep_mdmixmod_layered(distn, train),
		chained = .mstep_mdmixmod_chained(distn, train))
	pstep <- switch(topology,
		layered = .pstep_mdmixmod_layered(distn, train),
		chained = .pstep_mdmixmod_chained(distn, train))
	sstep <- .sstep_mdmixmod(bicsub, train)
	
	## initialize
	
	# build weights (quasi-E-step)
	
	if(prefit) {
	
		marginals = lapply(zvec, function(z)
			mixmod(X[[z]], K[[z]], family=family[[z]],
				train=train.raw, train.weight=train.weight,
				decreasing=decreasing[[z]], iter.max=iter.max, dname=zvec[[z]]))
		
		W = lapply(marginals, posterior)
		
		U = switch(topology,
			layered = listMean(lapply(zvec, function(z)
				if((K[[z]] == K0) && is.null(prior)) W[[z]] else
					posterior(mixmod(X[[z]], K0, family=family[[z]],
						prior=prior, iter.max=iter.max,
						train=train.raw, train.weight=train.weight)))),
			chained = if((K[[1]] == K0) && is.null(prior)) W[[1]] else
				posterior(mixmod(X[[1]], K0, family=family[[1]],
					prior=prior, iter.max=iter.max,
					train=train.raw, train.weight=train.weight)))
	
	} else {
		
		marginals = NULL # so we have something to return
			
		W = mapply(.init_weights, X, K, distn, SIMPLIFY=FALSE)
				
		U = switch(topology,
			layered = listMean(lapply(zvec, function(z)
				.init_weights(X[[z]], K0, distn[[z]], prior, 
					decreasing=decreasing[[z]]))),
			chained = .init_weights(X[[1]], K0, distn[[1]], prior, 
				decreasing=decreasing[[1]]))
	}
	
	V = switch(topology,
		layered = lapply(W, function(Wz) .vgenDep(U, Wz)),
		chained = c(list(.vgenDep(U, W[[1]])),
			lapply(2:Z, function(z) .vgenDep(W[[z-1]], W[[z]]))))
	names(V) = zvec # the "chained" option above produces an unnamed list
	
	weights = namedList(U, V, W)
	
	# initial M-, P-, and S-steps
	params = mstep(X, weights, K, kvec, K0, k0vec, zvec, prior=prior)
	pdfs = pstep(X, params, kvec, zvec)
	stats = sstep(0, weights, params, pdfs)
	
	# save initial information
	iteration.params = list(params)
	iteration.stats = list(stats)
	
	## iterate
	
	old_llik = stats[["llik"]] * 2 # dummy value to ensure iteration
	
	iter = 1
	
	while(abs(stats[["llik"]]/old_llik - 1) > LC_ITER_TOL && iter <= iter.max)
	{
		# preserve previous log-likelihood
		old_llik = stats[["llik"]]
		
		# E-step
		weights = estep(pdfs, params, zvec, kvec, k0vec)
						
		# M-step
		params = mstep(X, weights, K, kvec, K0, k0vec, zvec, params, 
			prior=prior)
				
		# post-EM calculations
		pdfs = pstep(X, params, kvec, zvec)
		stats = sstep(iter, weights, params, pdfs)
		
		# update iteration count and save results
		iter = iter + 1
		iteration.params[[iter]] = params
		iteration.stats[[iter]] = stats
			# increment iter before saving because iteration.params[[1]] is actually the parameters from iteration 0, etc., and the same for iteration.params
	}
	
	## final calculations
	
	weights = estep(pdfs, params, zvec, kvec, k0vec)
	
	posterior = weights$U
	rownames(posterior) = if(isMDF(X[[1]])) rownames(X[[1]]) else names(X[[1]])
	
	assignment = apply(posterior, 1, which.max)
	
	## package and return
	
	iteration.stats = Reduce(rbind, iteration.stats)
	rownames(iteration.stats) = NULL
	iteration.stats = as.data.frame(iteration.stats)
	
	retn = namedList(N, Z, D, K, K0, X, decreasing, 
		npar, npar.hidden, npar.observed, train, iter, params, stats, weights, pdfs, posterior, assignment, iteration.params, iteration.stats, topology, family, distn, prefit, iter.max, dname, dattr, zvec, kvec, k0vec, prior, marginals, estep, mstep, pstep, sstep)
	
	class(retn) = c("mdmixmod", "mixmod")
	
	return(retn)
}