### invisible_utilities:  invisible utility functions for the lcmix package

# Replace bad numeric values (NA, NaN, Inf) with 0
.bad2zero <- function(x)
{
	x[is.na(x)] = 0
	x[is.nan(x)] = 0
	x[is.infinite(x)] = 0
	
	return(x)
}

# Model E-step functions:  given no arguments; return functions of parameters and PDF matrices in the single-data case; or functions of parameters, PDF matrices, vector of data source names, and vectors of possible component values in the multiple-data case; which calculate weights.
.estep_mixmod <- function()
{
	function(pdfs)
	{
		lW = pdfs$lbeta - pdfs$lgamma
		W = exp(lW)
		
		namedList(W, lW)
	}
}
.estep_mdmixmod_layered <- function()
{
	function(pdfs, params, zvec, kvec, k0vec)
	{	
		## put parameters and PDFs into mathematical notation
		
		lQ = params$hidden$lcprob
		
		lG = pdfs$lG
		lGprime = pdfs$lGprime
		lbeta = pdfs$lbeta
		lgamma = pdfs$lgamma
		
		## now do the actual calculations
		
		lU = lbeta - lgamma
		U = exp(lU)
		
		lV = lapply(zvec, function(z) lapply(k0vec, function(k0)
			lU[,k0] + t(t(lG[[z]]) + lQ[[z]][k0,]) - lGprime[[z]][,k0]))
		V = rlapply(lV, exp, "matrix")
		
		W = lapply(V, listSum)
		
		## send it back

		namedList(U, V, W)
	}
}
.estep_mdmixmod_chained <- function()
{
	function(pdfs, params, zvec, kvec, k0vec)
	{
		## put parameters and PDFs into mathematical notation
		
		lP = params$hidden$lprob
		lQ = params$hidden$lcprob
		
		lalpha0 = pdfs$lalpha0
		lalpha = pdfs$lalpha
		lbeta0 = pdfs$lbeta0
		lbeta = pdfs$lbeta
		lgamma = pdfs$lgamma
		
		## now do the actual calculations
		
		lU = lbeta0 - lgamma
		U = exp(lU)
		lV1 = lapply(k0vec, function(k0)
			t(t(lalpha0[,k0] + lbeta[[1]])
			  + lQ[[1]][k0,] - lP[[1]]) - lgamma)
		lV2plus = lapply(2:length(zvec), function(z)
			lapply(kvec[[z-1]], function(kzm1)
				t(t(lalpha[[z-1]][,kzm1] + lbeta[[z]])
			  	  + lQ[[z]][kzm1,] - lP[[z]]) - lgamma))
		#lV426:  stay away!
		lV = c(list(lV1), lV2plus)
		names(lV) = zvec
		V = rlapply(lV, exp, "matrix")
		W = lapply(V, listSum)
		
		## send it back
		
		namedList(U, V, W)
	}
}

# Return a fully formed "dvec" vector given a matrix or data frame
.extract_dvec <- function(X)
{ 
	retn = 1:ncol(X)
	names(retn) = colnames(X)
	
	return(retn)
}

# Quickly calculates the inverse of double (NOT INTEGER) valued matrix X; does no sanity checking, offers no options, and in general should only be used when you're sure this is what you want.  To iterate over a list, and really speed things up, you can pass in the identity matrix; e.g, if Y is a list of 3x3 matrices, you can write "I3 = diag(3) ; res = lapply(Y, .fast_inv, I3)" which will be considerably faster than just "res = lapply(Y, .fast_inv)".  The fastest way to convert X from integer to double is just to call ".fast_inv(X+0.0)".  The ".list_fast_inv" wrapper is for convenience and speed; applied to a list, it returns a list.  Finally, "*sym()" guarantees that the return value will be symmetric, because inversion can sometimes introduce slight asymmetries; these methods should of course be used only with matrices that are supposed to be symmetric, like covariance matrices.
#*# TO DO:  either find a way to write equally fast code using pure R (no calls to "La_dgesv" or other base package ".Call" functions) or link directly to the relevant LAPACK functions.  Per <http://stats.stackexchange.com/questions/14951/efficient-calculation-of-matrix-inverse-in-r>, "chol2inv(chol(.))" may be much faster than "solve(.)" (and guarantee symmetry of the results?)
#*# .fast_inv <- function(X, ID=diag(nrow(X)))
#*# 	.Call("La_dgesv", X, ID, LC_EPS, PACKAGE="base")
.fast_inv  <- function(X, ID=diag(nrow(X))) solve(X, ID)
.fast_inv_sym <- function(X, ID=diag(nrow(X)))
{
	res = .fast_inv(X, ID)
	(res + t(res)) / 2.0
}
.list_fast_inv <- function(L)
	lapply(L, .fast_inv, ID=diag(nrow(L[[1]])))
.list_fast_inv_sym <- function(L)
	lapply(L, .fast_inv_sym, ID=diag(nrow(L[[1]])))

# Quickly calculates the logarithm of the absolute value of the determinant of double (NOT INTEGER) valued matrix X; does no sanity checking, offers no options, and in general should only be used when you're sure this is what you want.  The ".vec_fast_labs_det" wrapper is for convenience; applied to a list, it returns a vector.
#*# TO DO:  either find a way to write equally fast code using pure R (no calls to "det_ge_real" or other base package ".Call" functions) or link directly to the relevant LAPACK functions.
#*# .fast_labs_det <- function(X)
#*# 	.Call("det_ge_real", X, TRUE, PACKAGE="base")$modulus
.fast_labs_det <- function(X) determinant(X, logarithm=TRUE)$modulus
.vec_fast_labs_det <- function(L) sapply(L, .fast_labs_det)

# Calculate the initial weight matrix ("W") for data "X" with "K" means, means being chosen appropriately. The matrix that is returned is N-by-K.
.init_weights <- function(X, K, distn, prior=NULL, train=NULL, decreasing=TRUE)
{
	X = as.matrix(X)
	X = switch(distn,
		bern = X,
		mvbern = X,
		apply(X, 2, rank) # default
	)
	mu = switch(distn,
		bern = as.matrix((if(decreasing) K:1 else 1:K) - 0.5) / K,
			# compare to uniformQuantiles()
		mvbern = replicate(ncol(X), ((if(decreasing) K:1 else 1:K) - 0.5) / K),
		uniformQuantiles(X, K, decreasing=decreasing) # default
	)
	
	if(!is.null(train))
	{
		Kseq = 1:K
	
		trainIdx = which((colSums(train$Tprime) > 0)[Kseq])
		mu.old = mu
		Tstar = train$Tprime[,trainIdx, drop=FALSE]
		mu[trainIdx,] = (t(Tstar) %*% X) / colSums(Tstar)
		
		nonTrainIdx = setdiff(Kseq, trainIdx)
		nnt = length(nonTrainIdx)
		if(nnt) { # fill in non-training means with previously calculated values
			distSums = sapply(nonTrainIdx, function(nti)
				sum(manhattan(mu[trainIdx,], mu.old[nti,])))
			distRanks = rank(-distSums, ties.method="first")
				# "-" because we want to get the means which are farthest away from the means calculated from the training data
			useIdx = which(distRanks <= nnt)
			mu[nonTrainIdx,] = mu.old[nonTrainIdx, drop=FALSE][useIdx, drop=FALSE]
		}
	}
	
	D = sapply(seq(K), function(k) manhattan(X, mu[k,])) # distance matrix
	S = 1/(D+LC_EPS) # similarity matrix
	if(!is.null(prior)) S = t(t(S) * prior)
	
	retn = rowDiv(S)
	if(!is.null(train))
		retn[train$istrain,] = train$Tprime[train$istrain,1:K]
	return(retn)
}

# TESTING

.inititer_norm <- function(x, mu)
{
	# initialization (partial M- and P-steps)
	sigma = sqrt(colMeans(outer(x, mu, "-")^2))
		# we do E[X - mu]^2 here to avoid negative variances since the "means" are randomly chosen; below in the M-step, we use the simpler E[X^2] - mu^2 method
	G = mapply(function(mk, sk) dnorm(x, mk, sk), mu, sigma)
	
	# E- and M-steps
	W = rowDiv(G)
	csw = colSums(W)
	p = colMeans(W)
	mu = drop(x %*% W) / csw
	sigma = sqrt(colSums(x^2 * W)/csw - mu^2)
	
	# P- and S-steps
	G = mapply(function(mk, sk) dnorm(x, mk, sk), mu, sigma)
	gamma = drop(G %*% p)
	lgamma = log(gamma)
	llik = sum(lgamma)
	
	# prevent inversion
	muorder = order(mu, decreasing=TRUE)
	W = W[,muorder]
	mu = mu[muorder]
	p = p[muorder]
	
	# entropy
	entr = -2 * sum(W*log(W))
	
	# package and return
	namedList(W, mu, p, llik, entr)
}

.inititer_mvnorm <- function(X, mu)
{
	# gather useful statistics and munge data
	N = nrow(X)
	K = nrow(mu)
	kvec = 1:K
	tX = t(X)
	
	# initialization (partial M- and P-steps)
	Sigma = lapply(kvec, function(k) tcrossprod(tX - mu[k,]) / N)
		# see note on the initial "sigma" in .inititer_norm()
	G = sapply(kvec, function(k) dmvnorm(X, mu[k,], Sigma[[k]]))
	
	# E- and M-steps
	W = rowDiv(G)
	csw = colSums(W)
	p = colMeans(W)
	mu = t(tX %*% W) / csw
	Sigma = lapply(kvec, function(k)
		crossprod(X, W[,k]*X)/csw[k] - tcrossprod(mu[k,]))
	
	# P- and S-steps
	G = sapply(kvec, function(k) dmvnorm(X, mu[k,], Sigma[[k]]))
	gamma = drop(G %*% p)
	lgamma = log(gamma)
	llik = sum(lgamma)
	
	# prevent inversion
	muorder = order(rowSums(mu), decreasing=TRUE)
	W = W[,muorder]
	mu = mu[muorder,]
	p = p[muorder]
	
	# entropy
	entr = -2 * sum(W*log(W))
	
	# package and return
	namedList(W, mu, p, llik, entr)
}

.init_weights2 <- function(X, K, distn="normal", nstart=100, seed=123) # nstart=100, seed=NULL
{
	X = drop(as.matrix(X))
	
	multivariate = is.matrix(X)
	iterfun <- if(multivariate) .inititer_mvnorm else .inititer_norm
	sampfun <- if(multivariate) function()
		rowSort(rowSample(X, K), decreasing=TRUE)
	else function()
		sort(sample(X, K), decreasing=TRUE)
	
	N = ifelse(multivariate, nrow(X), length(X))
	if(multivariate) tX = t(X)
	if(!is.null(seed)) set.seed(seed)
	
	res = lapply(1:nstart, function(.) iterfun(X, sampfun()))
	llik = sapply(res, function(r) r$llik)
	idx = which.max(llik)
	
	retn = res[[idx]]$W
	attr(retn, "mu") = res[[idx]]$mu
	attr(retn, "p") = res[[idx]]$p
	attr(retn, "llik") = res[[idx]]$llik
	return(retn)
}

# /TESTING

# Given weight matrices "A" and "B", with column means a and b, calculate the transition matrix C such that:
# (1) a C = b
# (2) the value sum((A C - B)^2) is minimized
# using .dependentDiscreteCpdf (implying positive dependence).
.least_squares_transition_matrix <- function(A, B)
{
	a = colMeans(A)
	b = colMeans(B)
	toptim <- function(d) sum(((A %*% .dependentDiscreteCpdf(a, b, d)) - B)^2)
	
	.dependentDiscreteCpdf(a, b, optimize(toptim, c(0,1))$minimum)
}

# Deal with outliers in log-pdf-by-component matrices so extreme their density is effectively 0 for all components; the value -300 is sufficient to prevent instability but not actually change the results in any meaningful way
.lGfilter <- function(lG)
{
	if(min(lG) < -300) {
		K = ncol(lG)
		lG[(rowSums(lG < -300) == K), ] = -300
	}
	return(lG)
}
	
# Workhorse function for (unweighted) MLE of covariance; takes a D-by-N matrix, a D-length vector of means, and N, and returns a MLE covariance matrix.
.mleTcov <- function(tX, mu, N) { tXadj = tX - mu ; tcrossprod(tXadj) / N }

# Model M-step functions:  given distribution(s), return functions of data and weights (and data source names in the multiple-data-model case) which calculate parameters.
.mstep_mixmod <- function(distn, train)
{
	mstep.observed <- .mstep_observed(distn)
	
	function(X, weights, theta=NULL, prior=NULL) {
		
		W = weights$W
		if(!is.null(train)) W = W * train$twmult
		
		lprob = if(is.null(prior)) log(colMeans(W)) else log(prior)
		if(!is.null(train)) lprob = lprob - train$ltwcoef
		prob = exp(lprob)
		
		hidden = namedList(prob, lprob)
		
		observed = mstep.observed(X, W, theta=theta$observed)
		
		if(!is.null(train)) {
			tprob = train$tprob
			ltprob = train$ltprob
			tau = drop(t(W) %*% train$Tprime[,train$Ktrain])
			lrprob = rbind(train$lrpbase, log(tau) - log(sum(tau)))
			rprob = exp(lrprob)
			lcprob = t(ltprob + lrprob) - lprob
			cprob = exp(lcprob)
			train = namedList(tprob, ltprob, cprob, lcprob, rprob, lrprob)
		}
		
		namedList(hidden, observed, train)
	}
}
.mstep_mdmixmod_layered <- function(distn, train)
{
	mstep.observed <- mapply(.mstep_observed, distn)
	
	function(X, weights, K, kvec, K0, k0vec, zvec, params=NULL, prior=NULL)
	{
		U = weights$U
		V = weights$V
		W = weights$W
		
		if(!is.null(train)) {
			U = U * train$twmult
			V = rlapply(V, function(x) x * train$twmult,
					classes=c("matrix", "data.frame"))
			W = lapply(W, function(x) x * train$twmult)
		}
		
		lprob0 = if(is.null(prior)) log(colMeans(U)) else log(prior)
		ljprob = lapply(V, function(Vz)
			t(sapply(k0vec, function(k0) log(colMeans(Vz[[k0]])))))
		lprob = lapply(W, function(Wz) log(colMeans(Wz)))
		if(!is.null(train)) {
			lprob0 = lprob0 - train$ltwcoef
			ljprob = lapply(ljprob, function(jz) jz - train$ltwcoef)
			lprob = lapply(lprob, function(pz) pz - train$ltwcoef)
		}
		lcprob = lapply(ljprob, function(ljpz) ljpz - lprob0)
		
		prob0 = exp(lprob0)
		cprob = lapply(lcprob, exp)
		prob = lapply(lprob, exp)
		
		hidden = namedList(prob0, prob, cprob, lprob0, lprob, lcprob)
		
		observed = lapply(zvec, function(z)
			mstep.observed[[z]](X[[z]], W[[z]], theta=params$observed[[z]]))
		
		if(!is.null(train)) {
			tprob = train$tprob
			ltprob = train$ltprob
			tau = drop(t(U) %*% train$Tprime[,train$Ktrain])
			lrprob = rbind(train$lrpbase, log(tau) - log(sum(tau)))
			rprob = exp(lrprob)
			lcprob = t(train$ltprob + lrprob) - lprob0
			cprob = exp(lcprob)
			train = namedList(tprob, ltprob, cprob, lcprob, rprob, lrprob)
		}
		
		namedList(hidden, observed, train)
	}
}
.mstep_mdmixmod_chained <- function(distn, train)
{
	mstep.observed <- mapply(.mstep_observed, distn)
	
	function(X, weights, K, kvec, K0, k0vec, zvec, params=NULL, prior=NULL)
	{
		U = weights$U
		V = weights$V
		W = weights$W
		
		if(!is.null(train)) {
			U = U * train$twmult
			V = rlapply(V, function(x) x * train$twmult,
					classes=c("matrix", "data.frame"))
			W = lapply(W, function(x) x * train$twmult)
		}
		
		z2plus = 2:length(zvec)
		
		lprob0 = if(is.null(prior)) log(colMeans(U)) else log(prior)
		ljprob1 = t(sapply(k0vec, function(k0) log(colMeans(V[[1]][[k0]]))))
		ljprob2plus = lapply(z2plus, function(z)
			t(sapply(kvec[[z-1]], function(kzm1)
				log(colMeans(V[[z]][[kzm1]])))))
		ljprob = c(list(ljprob1), ljprob2plus)
		names(ljprob) = zvec
		lprob = lapply(W, function(Wz) log(colMeans(Wz)))
		if(!is.null(train)) {
			lprob0 = lprob0 - train$ltwcoef
			ljprob = lapply(ljprob, function(jz) jz - train$ltwcoef)
			lprob = lapply(lprob, function(pz) pz - train$ltwcoef)
		}
		lcprob1 = ljprob[[1]] - lprob0
		lcprob2plus = lapply(z2plus, function(z) ljprob[[z]] - lprob[[z-1]])
		lcprob = c(list(lcprob1), lcprob2plus)
		names(lcprob) = zvec
		lrprob = lapply(zvec, function(z) t(ljprob[[z]]) - lprob[[z]])
		
		prob0 = exp(lprob0)
		cprob = lapply(lcprob, exp)
		rprob = lapply(lrprob, exp)
		prob = lapply(lprob, exp)
		
		hidden = namedList(prob0, prob, cprob, rprob,
		     lprob0, lprob, lcprob, lrprob)
				
		observed = lapply(zvec, function(z)
			mstep.observed[[z]](X[[z]], W[[z]], theta=params$observed[[z]]))
		
		if(!is.null(train)) {
			tprob = train$tprob
			ltprob = train$ltprob
			tau = drop(t(U) %*% train$Tprime[,train$Ktrain])
			lrprob = rbind(train$lrpbase, log(tau) - log(sum(tau)))
			rprob = exp(lrprob)
			lcprob = t(train$ltprob + lrprob) - lprob0
			cprob = exp(lcprob)
			train = namedList(tprob, ltprob, cprob, lcprob, rprob, lrprob)
		}
		
		namedList(hidden, observed, train)
	}
}

# Return a function which performs weighted observed-distribution parameter estimation for the single-data mixture model M-step (can also be used in the multiple-data mixture model if handled appropriately.)
.mstep_observed <- function(distn) get(sprintf(".thetahat_%s", distn))

# Combined E- and M-steps for univariate Pearson Type VII parameter estimation.  Takes data, transposed data, weights, and a list "theta" containing the sum of weights, log of sum of weights, and current parameter values.  Returns updated theta with an additional element:  the value of the objective function.
.mvpvii_naive_emstep <- function(X, tX, w, sw, theta)
{
	# E-step
	
	D = ncol(X)
	N = nrow(X)
	X.prime = t(t(X) - theta$mean)
	iSigma = .fast_inv_sym(theta$scale)
	I = diag(D)
	
	alpha.prime = theta$shape + 0.5
	Psi = lapply(1:N, function(n) tcrossprod(X.prime[n,]))
	Lambda.prime = lapply(Psi, function(Psin) I + 0.5*Psin%*%iSigma)
	
	iLambda.prime = .list_fast_inv_sym(Lambda.prime)
	T.prime = lapply(iLambda.prime, "*", y=alpha.prime)
	mgammaterm = sum(digamma(alpha.prime - (0:(D-1))/2))
	labsdetterm = .vec_fast_labs_det(Psi)
	t.star = mgammaterm - labsdetterm
	
	wtp = lapply(1:N, function(n) w[n] * T.prime[[n]])
	swtp = listSum(wtp)
	wtpx = lapply(1:N, function(n) wtp[[n]] %*% X[n,])
	swtpx = listSum(wtpx)
	wpsitp = lapply(1:N, function(n) w[n] * Psi[[n]] %*% T.prime[[n]])
	swpsitp = listSum(wpsitp)
	swts = sum(w * t.star)

print(swpsitp)
	
	# M-step
		
	mu = drop(.fast_inv_sym(swtp) %*% swtpx)
	Sigma = swpsitp / sw
	
	toptim <- function(alpha) -lmgamma(alpha, D)*sw - alpha*swts
	res = nlminb(theta$shape, toptim, lower=LC_EPS+(D-1)/2)
	alpha = res$par
	
	# send it back
	list(mean=mu, scale=Sigma, shape=alpha, obj=res$objective)
}
.mvpvii_emstep <- function(X, tX, w, sw, theta)
{
	# E-step
	
	shape.prime = theta$shape + 0.5*ncol(X)
	tX.prime = tX - theta$mean
	X.prime = t(tX.prime)
	delta = rowSums((X.prime %*% .fast_inv_sym(theta$scale)) * X.prime)
	rate.prime = 1 + 0.5*delta
	
	tprime  = shape.prime / rate.prime
	tstar   = digamma(shape.prime) - log(rate.prime)
	
	# M-step
	
	wtp = w * tprime
	mu = colSums(wtp * X) / sum(wtp)
	
	D = ncol(X)
	#Sigma = tX.prime %*% (wtp * X.prime) / sw
	Sigma = (diag(D) + tX.prime %*% (wtp * X.prime)) / (1 + sw)
	
	swts = sum(w * tstar)
	topt <- function(alpha) lgamma(alpha)*sw - alpha*swts
	grad <- function(alpha) digamma(alpha)*sw - swts
	hess <- function(alpha) matrix(trigamma(alpha)*sw)
	res = nlminb(theta$shape, topt, gradient=grad, hessian=hess, lower=1+LC_EPS)
	alpha = res$par
	
	
	# send it back
	list(mean=mu, scale=Sigma, shape=alpha, obj=res$objective)
}

# Initial parameter estimation for multivariate Pearson Type VII.  Takes data, weights.  Returns the estimated parameters.
.mvpvii_init_params <- function(X, w)
{
	theta.mvnorm = thetahat.mvnorm(X, w)

	mean = theta.mvnorm$mean
	shape = 2
	scale = theta.mvnorm$cov * shape
	
	namedList(mean, scale, shape)
}

# Calculate the number of parameters for distribution of observed data in a single-data mixture model
.npar_observed <- function(distn, K, D)
{
	K * switch(distn,
		norm     = 2,
		mvnorm   = D + D*(D+1)/2,
		weisd    = 2,
		mvweisd  = 2*D + D*(D-1)/2,
		gamma    = 2,
		mvgamma  = 2*D + D*(D-1)/2,
		exp      = 1,
		mvexp    = D + D*(D-1)/2,
		pvii     = 3,
		mvpvii   = D + D*(D+1)/2 + 1,
		nbinom   = 2,
		mvnbinom = 2*D + D*(D-1)/2,
		pois     = 1,
		mvpois   = stop("multivariate Poisson not supported; use negative binomial (family 'nbinom') instead"),
		bern     = 1,
		mvbern   = D + D*(D-1)/2,
		stop(sprintf("unknown distribution '%s'", distn))
	)
}

# Return a nicely formatted string version of a family name, suitable for printing
.pretty_family <- function(family)
{
	switch(family,
		weibull   = "Weibull",
		pvii      = "Pearson VII",
		nbinom    = "negative binomial",
		poisson   = "Poisson",
		bernoulli = "Bernoulli",
		family)
}

# Model P-step functions:  given distribution(s), return functions of data, number(s) of components, and parameters (and data source names in the multiple-data-model case) which calculate PDF matrices
.pstep_mixmod <- function(distn, train)
{
	pstep.observed <- .pstep_observed(distn)
	
	function(X, params, kvec)
	{
		lG = .lGfilter(pstep.observed(X, kvec, params$observed))
		G = exp(lG) # g_n,k = f(x_n | y=k)
				
		if(is.null(train)) {
			lbeta = t(params$hidden$lprob + t(lG))
			beta = exp(lbeta)
			gamma = drop(G %*% params$hidden$prob) # f(X)
				# "drop" so we have an N-length vector rather than a N-by-1 matrix
			lgamma = log(gamma)
		} else {
			Gtrain = train$Tprime %*% t(params$train$cprob)
			lGtrain = log(Gtrain)
			lGstar = lGtrain + lG
			Gstar = exp(lGstar)
			lbeta = t(params$hidden$lprob + t(lGstar))
			beta = exp(lbeta)
			gamma = drop(Gstar %*% params$hidden$prob)
			lgamma = log(gamma)
		}
		
		namedList(G, beta, gamma, lG, lbeta, lgamma)
	}
}
.pstep_mdmixmod_layered <- function(distn, train)
{
	pstep.observed <- mapply(.pstep_observed, distn)
	
	function(X, params, kvec, zvec)
	{
		lG = lapply(zvec, function(z) .lGfilter(pstep.observed[[z]](
			X[[z]], kvec[[z]], params$observed[[z]])))
		G = lapply(lG, exp)
		
		Gprime = lapply(zvec, function(z)
			G[[z]] %*% t(params$hidden$cprob[[z]]))
		lGprime = lapply(Gprime, log)
		
		if(is.null(train)) {
			lGstar = listSum(lGprime)
			Gstar = exp(lGstar)
		} else {
			Gtrain = train$Tprime %*% t(params$train$cprob)
			lGtrain = log(Gtrain)
			lGstar = lGtrain + listSum(lGprime)
			Gstar = exp(lGstar)
		}
		
		lbeta0 = t(params$hidden$lprob0 + t(lGstar))
		beta0 = exp(lbeta0)
		gamma = drop(Gstar %*% params$hidden$prob0)
		lgamma = log(gamma)
		
		namedList(G, Gprime, beta0, gamma, lG, lGprime, lbeta0, lgamma)
	}
}
.pstep_mdmixmod_chained <- function(distn, train)
{
	
	pstep.observed <- mapply(.pstep_observed, distn)
	
	function(X, params, kvec, zvec)
	{	
		## initial component-specific PDF calculations
		
		lG = lapply(zvec, function(z) .lGfilter(pstep.observed[[z]](
			X[[z]], kvec[[z]], params$observed[[z]])))
		
		G = lapply(lG, exp)
		
		## put parameters into mathematical notation
			
		p0 = params$hidden$prob0
		P = params$hidden$prob
		Q = params$hidden$cprob
		R = params$hidden$rprob
		
		lP = params$hidden$lprob
		
		## now do the actual calculations
		
		Z = length(zvec)
		z2plus = 2:Z
		zminus = (Z-1):1
		N = nrow(G[[1]])
		
		alpha0 = if(is.null(train))
			outer(rep.int(1, N), p0)
		else
			train$Tpdf %*% params$train$rprob
		
		lalpha0 = log(alpha0)
		
		alpha = vector(mode="list", length=Z)
		lalpha = vector(mode="list", length=Z)
		lalpha[[1]] = log(alpha0 %*% Q[[1]]) + lG[[1]]
		alpha[[1]] = exp(lalpha[[1]])
		for(z in z2plus)
		{
			lalpha[[z]] = log(alpha[[z-1]] %*% Q[[z]]) + lG[[z]]
			alpha[[z]] = exp(lalpha[[z]])
		}
		names(alpha) = zvec
		
		beta = vector(mode="list", length=Z)
		lbeta = vector(mode="list", length=Z)
		lbeta[[Z]] = t(lP[[Z]] + t(lG[[Z]]))
		beta[[Z]] = exp(lbeta[[Z]])
		
		for(z in zminus)
		{
			lbeta[[z]] = log(beta[[z+1]] %*% R[[z+1]]) + lG[[z]]
			beta[[z]] = exp(lbeta[[z]])
		}
		names(beta) = zvec
		
		if(is.null(train)) {
			beta0 = beta[[1]] %*% R[[1]]
			lbeta0 = log(beta0)
		} else {
			Gtrain = train$Tprime %*% t(params$train$cprob)
			lGtrain = log(Gtrain)
			lbeta0 = log(beta[[1]] %*% R[[1]]) + lGtrain
			beta0 = exp(lbeta0)
		}
			
		gamma = rowSums(beta0)
		lgamma = log(gamma)
		
		## send it back
		
		namedList(G, alpha0, alpha, beta0, beta, gamma,
			lG, lalpha0, lalpha, lbeta0, lbeta, lgamma)
	}
}

# Return a function which performs PDF estimation for the observed data given the components of a single-data mixture model (can also be used in multiple-data mixture models if handled appropriately.)
.pstep_observed <- function(distn)
{
	switch(distn,
		
		norm = function(X, kvec, obspar) 
			sapply(kvec, function(k)
				dnorm(X, obspar$mean[k], obspar$sd[k], log=TRUE)),
		mvnorm = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvnorm(X, obspar$mean[k,], obspar$cov[[k]], log=TRUE)),
		
		weisd = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dweisd(X, obspar$shape[k], obspar$decay[k], log=TRUE)),
		mvweisd = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvweisd(X, obspar$shape[k,], obspar$decay[k,],
					obspar$corr[[k]], log=TRUE)),
		
		gamma = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dgamma(X, obspar$shape[k], obspar$rate[k], log=TRUE)),
		mvgamma = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvgamma(X, obspar$shape[k,], obspar$rate[k,],
					obspar$corr[[k]], log=TRUE)),
		
		exp = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dexp(X, obspar$rate[k], log=TRUE)),
		mvexp = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvexp(X, obspar$rate[k,],
					obspar$corr[[k]], log=TRUE)),
		
		pvii = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dpvii(X, obspar$mean[k], obspar$scale[k],
					obspar$shape[k], log=TRUE)),
		mvpvii = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvpvii(X, obspar$mean[k,], obspar$scale[[k]],
					obspar$shape[k], log=TRUE)),
		
		nbinom = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dnbinom(X, obspar$size[k], obspar$prob[k], log=TRUE)),
		mvnbinom = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvnbinom(X, obspar$size[k,], obspar$prob[k,],
					obspar$corr[[k]], log=TRUE)),
		
		pois = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dpois(X, obspar$lambda[k], log=TRUE)),
		mvpois = stop("multivariate Poisson not supported; use negative binomial (family 'nbinom') instead"),
		
		bern = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dpois(X, obspar$prob[k], log=TRUE)),
		mvbern = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvbern(X, obspar$prob[k,], obspar$corr[[k]], log=TRUE))
	
	)
}

# Combined E- and M-steps for univariate Pearson Type VII parameter estimation.  Takes data, weights, and a list "theta" containing the sum of weights, log of sum of weights, and current parameter values.  Returns updated theta with an additional element:  the value of the objective function.
.pvii_emstep <- function(x, w, sw, theta)
{	
	# E-step
	
	shape.prime = theta$shape + 0.5
	Delta = (x - theta$mean)^2
	delta = Delta / theta$scale
	rate.prime  = 1 + 0.5*delta
	
	tprime  = shape.prime / rate.prime
	tstar   = digamma(shape.prime) - log(rate.prime)
	
	# M-step
	
	wtp = w * tprime
	mu = sum(wtp * x) / sum(wtp)
	
	sigmasq = sum(wtp * Delta) / sw
	
	swts = sum(w * tstar)
	toptim <- function(alpha) lgamma(alpha)*sw - alpha*swts
	res = nlminb(theta$shape, toptim, lower=1+LC_EPS)
	alpha = res$par
	
	# send it back
	
	list(mean=mu, scale=sigmasq, shape=alpha, obj=res$objective)
}

# Initial parameter estimation for univariate Pearson Type VII.  Takes data, weights.  Returns the estimated parameters.
.pvii_init_params <- function(x, w)
{
	theta.norm = thetahat.norm(x, w)
	
	mean = theta.norm$mean
	shape = 2
	scale = theta.norm$var * shape
	
	namedList(mean, scale, shape)
}

# "Q-function" (expected log-likelihood) for single-data (.qfun_mixmod()) and multiple-data (.qfun_mdmixmod()) mixture models. If map==FALSE, the expected value will be returned; if TRUE, the MAP (maximum a posteriori) estimate will be returned.  The .bad2zero() calls are to deal with what happens when we try to take the log of 0 probabilities.
.qfun_mixmod <- function(weights, params, pdfs, train, map=FALSE)
{
	W = weights$W
	if(map) W = isRowMax(W)
	if(!is.null(train)) W = W * train$twmult
	
	hidden = sum(W %*% .bad2zero(params$hidden$lprob))
	observed = sum(weights$W * .bad2zero(pdfs$lG))
	
	train = if(is.null(train)) 0 else
		sum(W * (train$Tprime %*% t(.bad2zero(params$train$lcprob))))
		
	namedVector(hidden, observed, train)
}
.qfun_mdmixmod <- function(weights, params, pdfs, train, map=FALSE)
{
	U = weights$U
	V = weights$V
	W = weights$W
	
	zvec = names(W)
	names(zvec) = zvec
	Z = length(zvec)
	
	if(map) {
		U = isRowMax(U)
		V = lapply(V, .v2map)
		W = lapply(W, isRowMax)
	}
	
	if(!is.null(train)) {
		U = U * train$twmult
		V = rlapply(V, function(x) x * train$twmult,
				classes=c("matrix", "data.frame"))
		W = lapply(W, function(x) x * train$twmult)
	}
	
	lp = params$hidden$lprob0
	lQ = params$hidden$lcprob
	lG = pdfs$lG
	
	term1 = sum(U %*% .bad2zero(lp))
	term2 = sum(sapply(zvec, function(z)
		sum(sapply(1:length(V[[z]]), function(kparens)
			sum(V[[z]][[kparens]] %*% .bad2zero(lQ[[z]][kparens,]))))))
	
	term3 = sum(sapply(zvec, function(z)
		sum(W[[z]] * .bad2zero(lG[[z]]))))
	
	term4 = if(is.null(train)) 0 else
		sum(U * (train$Tprime %*% t(.bad2zero(params$train$lcprob))))
	
	c(hidden=term1+term2, observed=term3, train=term4)
}

# Weighted multivariate standard normal (means all equal to 0, variances all equal to 1) correlation (equivalent to covariance) matrix estimation.  The columns of X are already assumed to have mean 0, so there isn no need to estimate means and center.  The current engine function, .rhohatxy, is not terribly efficient and finding a better way to do this is TO DO.  (We can find the zeros of the derivative of the log-likelihood by solving a cubic (or maybe quartic?) equation, but it may have imaginary zeros or lie outside [-1, 1] or all kinds of other messiness.  Optimization is much easier.)
.rhohatxy <- function(x, y, w)
{
	t1 = sum(w)
	t2 = sum(w * (x-y)^2)
	t3 = sum(w * x * y)
	
	# negative log-likelihood for bivariate normal with standard marginals
	toptim <- function(rho)
	{
		u1 = 1 - rho^2
		u2 = 1 - rho
		
		0.5*log(u1)*t1 + 0.5*t2/u1 + u2*t3/u1
	}
	
	optimize(toptim, c(-1, 1))$minimum
}
.rhohat <- function(X, w)
{
	D = ncol(X)
	retn = diag(D)
	
	if(D > 1) for(i in 1:(D-1)) for(j in (i+1):D)
		retn[i,j] = retn[j,i] = .rhohatxy(X[,i], X[,j], w)
	
	drop(retn)
}

# Ensure that input data to a single-data mixture model (or a list component in a multiple-data mixture model) meets the requirements of the given distribution family.
.sanity_check_data <- function(X, family)
{
	if((family %in% LC_NONNEGFAM) && min(X) <= 0)
	{
		family = .pretty_family(family)
		stop(sprintf("%s family requires non-negative data", family))
	}
	
	return(0) # normal result, got through the above
}

# Turn a weight matrix into a matrix of almost-0 and almost-1 values, where the maximum in each row gets assigned almost-1 and all other positions almost-0
.sharpen_weight_matrix <- function(W)
{
	res = diag(ncol(W))[apply(W, 1, which.max),] + LC_EPS
	
	res / rowSums(res)
}

# Hidden data simulation:  take a number of samples to simulate, a distribution, and hidden parameters (and a topology for .simulate_hidden_mdmixdata()) and return hidden data (Y for .simulate_hidden_mixdata(), a list with elements $Y and $Y0 for .simulate_hidden_mdmixdata()).  No sanity checking is done, so everything has to be properly prepared; e.g., see calls to .simulate_hidden_mdmixdata() from simulateMdmixdata() to see how arguments should be set up.
.simulate_hidden_mixdata <- function(n, par)
{
	rcategorical(n, par$prob)
}
.simulate_hidden_mdmixdata <- function(n, par, topology)
{
	# gather information
	K0 = length(par$prob0)
	K = sapply(par$cprob, ncol)
	Z = length(par$cprob)
	zvec = names(par$cprob)
	names(zvec) = zvec

	# simulate top-level
	Y0 = rcategorical(n, par$prob0)
	
	# simulate lower-level
	if(topology == "layered") {
	
		Y = lapply(zvec, function(z)
		{
			res = vector(mode="numeric", length=n)
			for(k0 in 1:K0)
			{
				isk0 = (Y0 == k0)
				nisk0 = sum(isk0)
				res[isk0] = rcategorical(nisk0, par$cprob[[z]][k0,])
			}
			return(res)
		})
	
	} else if(topology == "chained") {
	
		Y = vector(mode="list", length=Z)
		Y[[1]] = vector(mode="numeric", length=n)
		
		for(k0 in 1:K0)
		{
			isk0 = (Y0 == k0)
			nisk0 = sum(isk0)
			Y[[1]][isk0] = rcategorical(nisk0, 
				par$cprob[[1]][k0,])
		}
		
		if(Z > 1) for(z in 2:Z)
		{
			Y[[z]] = vector(mode="numeric", length=n)
			for(kzm1 in 1:K[[z-1]])
			{
				iskzm1 = (Y[[z-1]] == kzm1)
				niskzm1 = sum(iskzm1)
				Y[[z]][iskzm1] = rcategorical(niskzm1, 
					par$cprob[[z]][kzm1,])
			}
		}
		
		names(Y) = zvec
	
	} else stop(sprintf("unrecognized topology '%s'", topology))
	
	# package it up, send it back
	namedList(Y, Y0)
}

# Model S-step functions:  given no arguments, return functions of iteration, weights, parameters, PDFs, and BIC adjustment factor (size of parameter space times log of sample size) which calculate iteration statistics
.sstep_mixmod <- function(bicsub, train)
{
	function(iter, weights, params, pdfs)
	{
		llik = ifelse(is.null(train), sum(pdfs$lgamma),
			sum(pdfs$lgamma*train$twmult))
		qval = sum(.qfun_mixmod(weights, params, pdfs, train))
		bic = 2*llik - bicsub
		iclbic = 2*qval - bicsub
		
		namedVector(iter, llik, qval, bic, iclbic)
	}
}
.sstep_mdmixmod <- function(bicsub, train)
{
	function(iter, weights, params, pdfs)
	{
		llik = ifelse(is.null(train), sum(pdfs$lgamma),
			sum(pdfs$lgamma*train$twmult))
		qval = sum(.qfun_mdmixmod(weights, params, pdfs, train))
		bic = 2*llik - bicsub
		iclbic = 2*qval - bicsub
		
		namedVector(iter, llik, qval, bic, iclbic)
	}
}

# Return a list containing the outer products of the rows of matrix X.
.tcrossprodcols <- function(X) lapply(1:ncol(X), function(n) tcrossprod(X[,n]))
.tcrossprodrows <- function(X) lapply(1:nrow(X), function(n) tcrossprod(X[n,]))

# Weighted parameter estimation functions for mixtures.  If argument "theta" is not NULL, it should be a list of parameters with the same names as those to be estimated, reflecting previous iteration's values of those parameters.  (Note that it will have no effect in the case of distributions with closed form MLEs.)
.thetahat_bern <- function(x, W, aslist=TRUE, theta=NULL, ...)
{
	res = sapply(1:ncol(W), function(k)
		thetahat.bern(x, W[,k], aslist=FALSE))
	
	p = as.vector(res)
		# the kth element of p is the success probability for the kth component
	
	if(aslist)
		list(prob=p)
	else
		c(prob=p)
}
.thetahat_mvbern <- function(X, W, aslist=TRUE, theta=NULL, ...)
{# TEMP version without dependency; fixing this is TO DO
	
	kvec = 1:ncol(W)
		
	res = sapply(kvec, function(k)
		thetahat.mvbern(X, W[,k], aslist=TRUE), simplify=FALSE)
			
	prob = t(sapply(kvec, function(k) res[[k]]$prob))
#	corr = lapply(kvec, function(k) res[[k]]$corr)
		
	if(aslist)
		namedList(prob) # namedList(prob, corr)
	else
		namedVector(prob) # namedVector(prob, corr)
}
.thetahat_exp <- function(x, W, aslist=TRUE, theta=NULL, ...)
{
	res = sapply(1:ncol(W), function(k)
		thetahat.exp(x, W[,k], aslist=FALSE))
	
	lambda = as.vector(res)
		# the kth element of lambda is the rate for the kth component
	
	if(aslist)
		list(rate=lambda)
	else
		c(rate=lambda)
}
.thetahat_mvexp <- function(X, W, aslist=TRUE, theta=NULL, ...)
{
	kvec = 1:ncol(W)
	
	res = sapply(kvec, function(k)
		thetahat.mvexp(X, W[,k], aslist=TRUE), simplify=FALSE)
	
	lambda = t(sapply(kvec, function(k) res[[k]]$rate))
		# the kth row of lambda is the vector of the column rates of X for the kth component
	rho = lapply(kvec, function(k) res[[k]]$corr)
	
	if(aslist)
		list(rate=lambda, corr=rho)
	else
		c(rate=lambda, corr=rho)
}
.thetahat_gamma <- function(x, W, aslist=TRUE, theta=NULL, ...)
{
	res = sapply(1:ncol(W), function(k)
		thetahat.gamma(x, W[,k], aslist=FALSE, shape.guess=theta$shape[k]))
	
	sigma = res["shape",]
	lambda = res["rate",]
		# the kth element of sigma is the shape for the kth component, and similarly for lambda (rate)
	
	if(aslist)
		list(shape=sigma, rate=lambda)
	else
		c(shape=sigma, rate=lambda)
}
.thetahat_mvgamma <- function(X, W, aslist=TRUE, theta=NULL, ...)
{
	kvec = 1:ncol(W)
	
	res = sapply(kvec,
	             function(k)
	             	thetahat.mvgamma(X, W[,k], aslist=TRUE, 
	             		shape.guess=theta$shape[k,]),
		         simplify=FALSE)
	
	sigma = t(sapply(kvec, function(k) res[[k]]$shape))
		# the kth row of sigma is the vector of the column shapes of X for the kth component
	lambda = t(sapply(kvec, function(k) res[[k]]$rate))
		# the kth row of lambda is the vector of the column rates of X for the kth component
	rho = lapply(kvec, function(k) res[[k]]$corr)
	
	if(aslist)
		list(shape=sigma, rate=lambda, corr=rho)
	else
		c(shape=sigma, rate=lambda, corr=rho)
}
.thetahat_norm <- function(x, W, aslist=TRUE, theta=NULL, ...)
{
	res = sapply(1:ncol(W), function(k)
		thetahat.norm(x, W[,k], aslist=FALSE))
	
	mu = res["mean",]
	sigmasq = res["var",]
	sigma = res["sd",]
		# the kth element of mu is the mean for the kth component, and similarly for sigmasq (variance) and sigma (std. dev.)
	
	if(aslist)
		list(mean=mu, var=sigmasq, sd=sigma)
	else
		c(mean=mu, var=sigmasq, sd=sigma)
}
.thetahat_mvnorm <- function(X, W, aslist=TRUE, theta=NULL, ...)
{
	
	kvec = 1:ncol(W)
	
	res = sapply(kvec, function(k)
		thetahat.mvnorm(X, W[,k], aslist=TRUE), simplify=FALSE)
	
	mu = t(sapply(kvec, function(k) res[[k]]$mean))
		# the kth row of mu is the vector of the column means of X for the kth component
	Sigma = lapply(kvec, function(k) res[[k]]$cov)
		# the kth list element of Sigma is the covariance matrix for the kth component
		
	if(aslist)
		list(mean=mu, cov=Sigma)
	else
		c(mean=mu, cov=Sigma)
}
.thetahat_pvii <- function(x, W, aslist=TRUE, theta=NULL, iter.max=0, ...)
{
	res = sapply(1:ncol(W), function(k) {
		if(!is.null(theta)) theta = lapply(theta, function(par) par[k])
			# get component-specific parameters to pass to estimation fn.
		thetahat.pvii(x, W[,k], aslist=FALSE, iter.max=iter.max, theta=theta)
	})
	
	mean  = res["mean", ]
	scale = res["scale",]
	shape = res["shape",]
	
	if(aslist)
		namedList(mean, scale, shape)
	else
		namedVector(mean, scale, shape)
}
.thetahat_mvpvii <- function(X, W, aslist=TRUE, theta=NULL, iter.max=0, ...)
{
	kvec = 1:ncol(W)
	
	res = lapply(kvec, function(k) {
		if(!is.null(theta)) theta = lapply(theta, function(par)
				if(is.list(par)) par[[k]]
				else if(is.matrix(par)) par[k,]
				else par[k])
					# get component-specific parameters to pass on
		thetahat.mvpvii(X, W[,k], aslist=TRUE, iter.max=iter.max, theta=theta)
	})
		
	mean  = t(sapply(kvec, function(k) res[[k]]$mean))
	scale = lapply(kvec, function(k) res[[k]]$scale)
	shape = sapply(kvec, function(k) res[[k]]$shape)
		# see .thetahat_mvnorm() above to understand the formation of these data structures
	
	if(aslist)
		namedList(mean, scale, shape)
	else
		namedVector(mean, scale, shape)
}
.thetahat_weisd <- function(x, W, aslist=TRUE, theta=NULL, ...)
{
	res = sapply(1:ncol(W), function(k)
		thetahat.weisd(x, W[,k], aslist=FALSE, shape.guess=theta$shape[k]))
	
	sigma = res["shape",]
	delta = res["decay",]
		# the kth element of sigma is the shape for the kth component, and similarly for delta (decay)
	
	if(aslist)
		list(shape=sigma, decay=delta)
	else
		c(shape=sigma, decay=delta)
}
.thetahat_mvweisd <- function(X, W, aslist=TRUE, theta=NULL, ...)
{
	kvec = 1:ncol(W)
	
	res = sapply(kvec,
	             function(k)
	             	thetahat.mvweisd(X, W[,k], aslist=TRUE, 
	             		shape.guess=theta$shape[k,]),
		         simplify=FALSE)
	
	sigma = t(sapply(kvec, function(k) res[[k]]$shape))
		# the kth row of sigma is the vector of the column shapes of X for the kth component
	delta = t(sapply(kvec, function(k) res[[k]]$decay))
		# the kth row of delta is the vector of the column decays of X for the kth component
	rho = lapply(kvec, function(k) res[[k]]$corr)
	
	list(shape=sigma, decay=delta, corr=rho)
}
.thetahat_nbinom <- function(x, W, aslist=TRUE, theta=NULL, ...)
{
	res = sapply(1:ncol(W), function(k)
		thetahat.nbinom(x, W[,k], aslist=FALSE, size.guess=theta$size[k]))
	
	size = res["size",]
	prob = res["prob",]
	
	if(aslist)
		namedList(size, prob)
	else
		namedVector(size, prob)
}
.thetahat_mvnbinom <- function(X, W, aslist=TRUE, theta=NULL, ...)
{# TEMP version without dependency; fixing this is TO DO
	
	kvec = 1:ncol(W)
	
	res = sapply(kvec,
	             function(k)
	             	thetahat.mvnbinom(X, W[,k], aslist=TRUE, 
	             		size.guess=theta$size[k,]),
		         simplify=FALSE)
	
	size = t(sapply(kvec, function(k) res[[k]]$size))
	prob = t(sapply(kvec, function(k) res[[k]]$prob))
#	corr = lapply(kvec, function(k) res[[k]]$corr) TO DO
	
	namedList(size, prob) # namedList(size, prob, corr) TO DO
}
.thetahat_pois <- function(x, W, aslist=TRUE, theta=NULL, ...)
{
	res = sapply(1:ncol(W), function(k)
		thetahat.pois(x, W[,k], aslist=FALSE))
	
	lambda = as.vector(res)
		# the kth element of lambda is the lambda for the kth component
	
	if(aslist)
		namedList(lambda)
	else
		namedVector(lambda)
}
.thetahat_mvpois <- function(X, W, aslist=TRUE, theta=NULL, ...)
{
	stop("multivariate Poisson not supported; use negative binomial (family 'nbinom') instead")
}

# Build "Q" (transition) matrices under under local independence assumption (.qgen) or under linear dependence assumption (.qgenDep)
.qgen <- function(A, B)
{
	J = crossprod(A, B) # N f(A,B)
	
	J / rowSums(J) # f(B|A)
}
.qgenDep <- function(A, B)
{
	res = sapply(1:ncol(B), function(kb) coef(nnls(A, B[,kb]))) # N f(B|A)
	
	res / rowSums(res) # f(B|A)
}

# error function for unsupported multivariate distribution families
.usmv <- function(family)
	stop(sprintf("multivariate family '%s' not yet supported", family))

# Build "V" (joint weight) matrices under local independence assumption (.vgen) or under linear dependence assumption (.vgenDep)
.vgen <- function(A, B) lapply(1:ncol(A), function(ka) A[,ka] * B)
.vgenDep <- function(A, B)
{
	R = .qgenDep(B, A) # f(A|B)
	tB = t(B)
	
	lapply(1:ncol(A), function(ka) t(tB * R[,ka]))
}

# Convert "V" (joint weight) matrix into MAP estimate.  TO DO:  this is horribly slow.
.v2map <- function(V)
{
	N = nrow(V[[1]])
	Ka = length(V)
	
	map.raw = lapply(1:N, function(n) { # MAP(V) indexed by n instead of ka
		res = t(sapply(1:Ka, function(ka) V[[ka]][n,]))
		1 * (res == max(res))
	})
	
	lapply(1:Ka, function(ka) t(sapply(1:N, function(N) map.raw[[N]][ka,])))
}