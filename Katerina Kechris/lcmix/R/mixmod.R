### mixmod.R:  define the mixmod (single data mixture models) function for the lcmix package.

mixmod <- function(X, K, family=names(LC_FAMILY), prior=NULL, train=NULL, 
	train.weight=1, init=0, decreasing=TRUE, iter.max=LC_ITER_MAX, 
	dname=deparse(substitute(X)))
{
	## munge data, gather preliminary information, and do sanity checks
	
	dname = dname
		# keep the implicit call to deparse(substitute(X)) from getting screwed up by whatever manipulation we do below (this is a nasty side-effect of lazy evaluation, I suspect)
	dattr = attributes(X)
	
	X = drop(as.matrix(X)) # get matrix or vector form
	multivariate = is.matrix(X)
	if(multivariate) {
		N = nrow(X)
		D = ncol(X)
		univariate = FALSE
	} else {
		N = length(X)
		D = 1
		univariate = TRUE
	}
	
	family = match.arg(family)
	dummy = .sanity_check_data(X, family) # if we get through this, we're okay
	distn = ifelse(univariate, LC_FAMILY[[family]]$uni,
		LC_FAMILY[[family]]$multi)
		
	kvec = 1:K
	
	npar.hidden = K-1
	npar.observed = .npar_observed(distn, K, D)
	npar = npar.hidden + npar.observed
	
	## handle training PRN
	
	if(train.weight == 0) train = NULL
	if(train.weight < 0) stop("negative 'train.weight' is not valid")
	
	if(is.null(train)) bicsub = npar * log(N) else {
		
		# sanity checks and component fixing
		if(length(train) != N)
			stop("'train' must be NULL or of same length as data")
		if(!all(train %in% (0:(K+1))))
			stop("'train' values must be integers between 0 and K+1 inclusive")
		train[(train==0)] = K+1
		
		# expansion
		weight = train.weight
		Ktrain = K+1
		Tprime = diag(Ktrain)[train,]
		istrain = (train <= K)
		wTprime = Tprime
		wTprime[istrain,] = weight * wTprime[istrain,]
		ltprob = log(colSums(wTprime)) - log(sum(wTprime))
		tprob = exp(ltprob)
		rpbase = diag(K)
		lrpbase = log(rpbase)
		
		# vector for multiplying labeled rows of weight matrices, effective data set size, and coefficient for dividing column means calculated from such training-weight-multiplied matrices to give proper probabilities
		twmult = rep.int(1, N) + (weight-1)*istrain
		Neff = sum(twmult)
		ltwcoef = log(Neff) - log(N)
		twcoef = exp(ltwcoef)
		
		# package it up
		train = namedList(train, weight, Ktrain, Tprime, istrain, tprob, ltprob, 
			rpbase, lrpbase, twmult, Neff, twcoef, ltwcoef)
		
		# additional calculations
		bicsub = npar * log(Neff)
	}
	
	## define internal functions
	
	estep <- .estep_mixmod()
	mstep <- .mstep_mixmod(distn, train)
	pstep <- .pstep_mixmod(distn, train)
	sstep <- .sstep_mixmod(bicsub, train)
	
	## initialize
	
	# quasi E-step
	W = if(0 && init) {
		
		xr = apply(as.matrix(X), 2, rank, ties.method="random")
		xrr = rank(rowSums(xr), ties.method="random")
		
		res = lapply(1:init, function(i)
		{
			idx = sample(1:N, K)
			idx = idx[order(xrr[idx], decreasing=TRUE)]
			mu = as.matrix(xr[idx,])
						
			D = sapply(1:K, function(k) manhattan(xr, mu[k,])) # distance matrix
			S = 1/(D+LC_EPS) # similarity matrix
			if(!is.null(prior)) S = t(t(S) * prior)
			W = rowDiv(S)
print(cor(W, X))			
			weights = list(W=W)
			#for(j in 1:init) {
				params = mstep(X, weights, prior=prior)
				pdfs = pstep(X, params, kvec)
				weights = estep(pdfs)
			#}
			stats = sstep(0, weights, params, pdfs)
			list(W=weights$W, Q=stats[["qval"]])
		})
		
		res[[which.max(sapply(res, function(r) r$Q))]]$W
	
	} else if(1 && init) {
		theta = get(sprintf("thetahat.%s", distn))(X)
		res = replicate(init, simplify=FALSE, expr={
			idx = sample(1:N, K)
			mu = if(multivariate) X[idx,] else X[idx]
#if(multivariate) mu = apply(mu, 2, sort, decreasing=decreasing)
			G = sapply(1:K, function(k) switch(distn,
				norm = dnorm(X, mu[k], theta$sd),
				mvnorm = dmvnorm(X, mu[k,], theta$cov),
				pvii = dpvii(X, mu[k], theta$scale, theta$shape),
				mvpvii = dmvpvii(X, mu[k,], theta$scale, theta$shape)
			))
			for(i in 1:10)
			{
				weights = if(i == 1)
					list(W=.bad2zero(rowDiv(G)))
				else
					estep(pdfs)
				params = mstep(X, weights, prior=prior)
				pdfs = pstep(X, params, kvec)
			}
			stats = sstep(0, weights, params, pdfs)
			list(W=weights$W, val=stats[["llik"]])
		})
		W = res[[which.max(sapply(res, function(r) r$val))]]$W
		xr = apply(as.matrix(X), 2, rank)
		mu = drop((t(W) %*% xr) / colSums(W))
		mus = if(multivariate) rowSums(log(mu)) else mu
		W = W[,order(mus, decreasing=TRUE)]
print(drop((t(W) %*% as.matrix(X)) / colSums(W)))
		W
	} else .init_weights(X, K,
		distn=distn, prior=prior, train=train, decreasing=decreasing)
	weights = namedList(W)
	
	# M-step
	params = mstep(X, weights, prior=prior)
						
	# post-EM calculations
	pdfs = pstep(X, params, kvec)
	stats = sstep(0, weights, params, pdfs)
	
	# save results
	iteration.params = list(params)
	iteration.stats = list(stats)
	
	## iterate
	
	old_llik = stats[["llik"]] - 100 # dummy value to ensure iteration
		
	iter = 1
	
	while(abs(stats[["llik"]]/old_llik - 1) > LC_ITER_TOL && iter <= iter.max)
	{
		# preserve previous log-likelihood
		old_llik = stats[["llik"]]
		
		# E-step
		weights = estep(pdfs)
		
		# M-step
		params = mstep(X, weights, params, prior=prior)
				
		# post-EM calculations
		pdfs = pstep(X, params, kvec)
		stats = sstep(iter, weights, params, pdfs)
		
		# update iteration count and save results
		iter = iter + 1
		iteration.params[[iter]] = params
		iteration.stats[[iter]] = stats
			# increment iter before saving because iteration.params[[1]] is actually the parameters from iteration 0, etc., and the same for iteration.params	
	}
	
	## final calculations
	
	weights = estep(pdfs)
	
	posterior = weights$W
	rownames(posterior) = if(isMDF(X)) rownames(X) else names(X)
	
	assignment = apply(posterior, 1, which.max)
	
	## package and return
	
	iteration.stats = Reduce(rbind, iteration.stats)
	rownames(iteration.stats) = NULL
	iteration.stats = as.data.frame(iteration.stats)
	
	retn = namedList(N, D, K, X, decreasing, npar, npar.hidden, npar.observed, 
		train, iter, params, stats, weights, pdfs, posterior, assignment, iteration.params, iteration.stats, family, distn, prior, iter.max, dname, dattr, kvec, estep, mstep, pstep, sstep)
	
	class(retn) = "mixmod"
	
	return(retn)
}
