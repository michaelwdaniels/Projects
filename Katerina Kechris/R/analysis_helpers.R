source('R/includes.R')
# --------------------------------- WRAPPERS --------------------------------- #

# Functions to wrap the modeling and plotting portions of the analyses, which
# are repeated under lots of conditions (i.e. training set size)

divideYeastData <- function(feature.table, date.table, labels, 
                            train.date = NULL, 
                            train.class = c('positive', 'negative', 'both'),
                            num.train.examples = NULL, seed, negatives = c('all','constant','match','match_test'), constant = constant){
  
  # Divide the yeast data set into a training and testing set based on desired
  # example class, date a positive example was discovered, or number of training
  # examples desired.
  
  set.seed(seed)
  
  # Format arguments
  train.class = train.class[1]
  if (!is.null(train.date) & train.class != 'positive'){ 
    stop('Dates only available when using positive training examples only')}
  if (train.class == 'both' & is.null(num.train.examples)){
    stop('If using both positive and negative training examples, must specify 
         the number of training examples to use')
  }
  
  # If using both pos and neg training examples, just split by the number of
  # training examples desired and we are done
  if (train.class == 'both'){
    is.train = row.names(feature.table) %in% sample(row.names(feature.table), 
                                                    num.train.examples)
    train.data = feature.table[is.train,]
    test.data = feature.table[!is.train,]
  } 
  
  # If using only pos training examples, first separate by label class
  if (train.class == 'positive'){
    train.data = feature.table[labels == 1,]
    test.data = feature.table[labels == 0,]
    if (negatives == 'all'){
      test.data <- test.data}
    if (negatives == 'match'){
      idx <- sample(1:nrow(test.data),train.size) # match number of negatives to current training size
      test.data <- test.data[idx,]}
    if (negatives == 'constant'){
      idx <- sample(1:nrow(test.data),constant) # match number of negatives to the total number of positive labels
      test.data <- test.data[idx,]}
    if (negatives == 'match_test'){
      inTrain <- sample(1:length(labels),size = floor(p[i]*length(labels)), replace = FALSE)
      
      training <- feature.table[inTrain,]
      
      idx <- sample(1:length(which(labels==1)),size = floor(p[i]*length(which(labels==1))), replace = FALSE)
      
      test.data <- test.data[idx,]}
      
    # Then separate by date
    if (!is.null(train.date)){
      good.date = unique(row.names(date.table[date.table$year < train.date,]))
      is.good = row.names(train.data) %in% good.date
      test.data = rbind(test.data, train.data[!is.good,])
      train.data = train.data[is.good,]
    }
    
    # Then separate by number of training examples
    if (!is.null(num.train.examples)){
      if (num.train.examples < nrow(train.data)){
        is.train = row.names(train.data) %in% sample(row.names(train.data), 
                                                        num.train.examples)
        test.data = rbind(test.data, train.data[!is.train,])
        train.data = train.data[is.train,]
      } else{
        stop('Number of training examples desired is larger than training set')
      }
    }
    is.train = row.names(feature.table) %in% row.names(train.data)
  }
  
  # If using only neg training examples, first separate by label class
  if (train.class == 'negative'){
    train.data = feature.table[labels == 0,]
    test.data = feature.table[labels == 1,]
    
    # The separate by number of training examples
    if (!is.null(num.train.examples)){
      if (num.train.examples < nrow(train.data)){
        is.train = row.names(train.data) %in% sample(row.names(train.data), 
                                                     num.train.examples)
        test.data = rbind(test.data, train.data[!is.train,])
        train.data = train.data[is.train,]
      } else{
        stop('Number of training examples desired is larger than training set')
      }
    }
    is.train = row.names(feature.table) %in% row.names(train.data)
  }

  # Return divided data
  results = list(feature.table = feature.table,
                 is.train = is.train,
                 train = train.data,
                 test = test.data)
  return(results)
}

divideYeastData.jack <- function(feature.table, date.table, labels, 
                                 train.date = NULL, 
                                 train.class = c('positive', 'negative', 'both'),
                                 num.train.examples = NULL, seed, k, negatives = c('all','constant','match','match_test'),constant){
  
  # Divide the yeast data set into a training and testing set based on desired
  # example class, date a positive example was discovered, or number of training
  # examples desired. k is the leave one out iteration
  
  set.seed(seed)
  
  # Format arguments
  train.class = train.class[1]
  if (!is.null(train.date) & train.class != 'positive'){ 
    stop('Dates only available when using positive training examples only')}
  if (train.class == 'both' & is.null(num.train.examples)){
    stop('If using both positive and negative training examples, must specify 
         the number of training examples to use')
  }
  
  # If using both pos and neg training examples, just split by the number of
  # training examples desired and we are done
  if (train.class == 'both'){
    is.train = row.names(feature.table) %in% sample(row.names(feature.table), 
                                                    num.train.examples)
    train.data = feature.table[is.train,]
    test.data = feature.table[!is.train,]
  } 
  
  # If using only pos training examples, first separate by label class
  if (train.class == 'positive'){
    train.data = feature.table[labels == 1,]
    test.data = feature.table[labels == 0,]
    if (negatives == 'all'){
      test.data <- test.data}
    if (negatives == 'match'){
      idx <- sample(1:nrow(test.data),train.size) # match number of negatives to current training size
      test.data <- test.data[idx,]}
    if (negatives == 'match_test'){
      idx <- sample(1:nrow(test.data),num.train.examples) # match number of negatives to number of positives in test set
      test.data <- test.data[idx,]}
    if (negatives == 'constant'){
       idx <- sample(1:nrow(test.data),constant) # match number of negatives to the total number of negatives
       test.data <- test.data[idx,]}
       
    # Then separate by date
    if (!is.null(train.date)){
      good.date = unique(row.names(date.table[date.table$year < train.date,]))
      is.good = row.names(train.data) %in% good.date
      test.data = rbind(test.data, train.data[!is.good,])
      train.data = train.data[is.good,]
    }
    
    # Then separate by number of training examples
    if (!is.null(num.train.examples)){
      if (num.train.examples < nrow(train.data)){
        is.train = row.names(train.data) %in% row.names(train.data)[-k]
        test.data = rbind(test.data, train.data[!is.train,])
        train.data = train.data[is.train,]
      } else{
        stop('Number of training examples desired is larger than training set')
      }
    }
    is.train = row.names(feature.table) %in% row.names(train.data)
  }
  
  # If using only neg training examples, first separate by label class
  if (train.class == 'negative'){
    train.data = feature.table[labels == 0,]
    test.data = feature.table[labels == 1,]
    
    # The separate by number of training examples
    if (!is.null(num.train.examples)){
      if (num.train.examples < nrow(train.data)){
        is.train = row.names(train.data) %in% sample(row.names(train.data), 
                                                     num.train.examples)
        test.data = rbind(test.data, train.data[!is.train,])
        train.data = train.data[is.train,]
      } else{
        stop('Number of training examples desired is larger than training set')
      }
    }
    is.train = row.names(feature.table) %in% row.names(train.data)
  }
  
  # Return divided data
  results = list(feature.table = feature.table,
                 is.train = is.train,
                 train = train.data,
                 test = test.data)
  return(results)
  }


getLCMIXmodels <- function(feature.table, is.train, test.labels, K0 = 2, TW = 10){
  
  BINARY_FEATURES2 = BINARY_FEATURES[BINARY_FEATURES %in% names(feature.table)]
  INTEGER_FEATURES2 = INTEGER_FEATURES[INTEGER_FEATURES %in% names(feature.table)]
  OPEN_FEATURES2 = OPEN_FEATURES[OPEN_FEATURES %in% names(feature.table)]
  CLOSED_FEATURES2 = CLOSED_FEATURES[CLOSED_FEATURES %in% names(feature.table)]
  
  # Make feature tables
  X = list(
    bern = do.call(cbind, feature.table[,BINARY_FEATURES2]),
    nbin = do.call(cbind, feature.table[,INTEGER_FEATURES2]),
    norm = cbind(
      do.call(cbind, feature.table[,OPEN_FEATURES2]),
      log(do.call(cbind, feature.table[,CLOSED_FEATURES2]))
    )
  )
  
  # Marginal model selection
  message(notification("marginal model selection", 1))
  margsel = msapply(namedList(2,3,4,5), function(K)
    mapply(X, names(X), FUN=function(x, fam)
      iclbic(mixmod(x, K=K, family=fam), map=T)))
  K = apply(margsel, 1, which.max) + 1
  
  # Joint model fitting
  message(notification("fitting joint models", 1))
  train = as.numeric(is.train)
  fm1 = mcparallel(mdmixmod(X, K=K, K0=K0, family=names(X),
                            train=train, train.weight=TW)) 
  fm2 = mcparallel(mdmixmod(X, K=K, K0=K0, family=names(X)))
  fits = mccollect(list(fm1, fm2))
  names(fits) = c("semisup", "unsup")
  
  # Make predictions for semisup and unsup models
  preds = lapply(fits, function(x) predict(x, type="prob")[!is.train, 1])
  
  # Report performance
  reports = c("precis", "recall", "f1meas")
  suppressWarnings(rept1 <- as.data.frame(sapply(preds, perfmeas, 
                                                 labels=test.labels)[reports,]))
  rept1$report = reports
  rept1$FDR = 0
  suppressWarnings(rept2 <- as.data.frame(sapply(preds, perfmeas, lfdr=.20, 
                                                 labels=test.labels)[reports,]))
  rept2$report = reports
  rept2$FDR = .2
  suppressWarnings(rept3 <- as.data.frame(sapply(preds, perfmeas, lfdr=.01, 
                                                 labels=test.labels)[reports,]))
  rept3$report = reports
  rept3$FDR = .01
  combined.report = rbind(melt(rept1, id.vars = c('report', 'FDR'), 
                               measure.vars = c('semisup', 'unsup')),
                          melt(rept2, id.vars = c('report', 'FDR'), 
                               measure.vars = c('semisup', 'unsup')),
                          melt(rept3, id.vars = c('report', 'FDR'), 
                               measure.vars = c('semisup', 'unsup')))
  
  # Return prediction probabilities and performance report
  return(list(preds = preds, performance = combined.report))
  
}


mdmixmod1 <- function (X, K, K0 = min(K), topology = LC_TOPOLOGY, family = NULL, 
          prior = NULL, train = NULL, train.weight = 1, decreasing = TRUE, 
          prefit = TRUE, iter.max = LC_ITER_MAX, dname = deparse(substitute(X))) 
{
  if (!is.list(X)) 
    stop("'X' must be a list")
  dname = dname
  dattr = lapply(X, attributes)
  X = lapply(X, function(Xz) drop(as.matrix(Xz)))
  Z = length(X)
  if (Z < 2) 
    stop("single data source, use mixmod() instead")
  if (is.null(names(X))) 
    names(X) = paste("X", 1:Z, sep = "")
  zvec = names(X)
  Z = length(zvec)
  names(zvec) = zvec
  univariate = !sapply(X, is.matrix)
  ntmp = sapply(zvec, function(z) ifelse(univariate[[z]], length(X[[z]]), 
                                         nrow(X[[z]])))
  if (min(ntmp) != max(ntmp)) {
    stop("unequal data sizes")} else N = ntmp[[1]]
  D = sapply(zvec, function(z) ifelse(univariate[[z]], 1, ncol(X[[z]])))
  topology = match.arg(topology)
  if (is.null(family)) {
    family = rep("normal", Z)} else {
    family = match.arg(family, names(LC_FAMILY), several.ok = TRUE)
    if (length(family) > Z) 
      warning("length(family) > Z, truncating to fit")
    if (length(family) != Z) 
      family = rep(family, length.out = Z)
  }
  names(family) = zvec
  dummy = lapply(zvec, function(z) .sanity_check_data(X[[z]], 
                                                      family[[z]]))
  distn = sapply(zvec, function(z) ifelse(univariate[[z]], 
                                          LC_FAMILY[[(family[[z]])]]$uni, LC_FAMILY[[(family[[z]])]]$multi))
  if (length(K) > Z) 
    warning("length(K) > Z, truncating to fit")
  if (length(K) != Z) 
    K = rep(K, length.out = Z)
  names(K) = zvec
  kvec = lapply(zvec, function(z) 1:K[[z]])
  k0vec = 1:K0
  if (length(decreasing) > Z) 
    warning("length(decreasing) > Z, truncating to fit")
  if (length(decreasing) != Z) 
    decreasing = rep(decreasing, length.out = Z)
  names(decreasing) = zvec
  npar.hidden = switch(topology, layered = K0 - 1 + K0 * sum(K - 1), chained = K0 - 1 + sum(c(K0, K[1:(Z - 1)]) * (K - 1)))
  npar.observed = sum(sapply(zvec, function(z) .npar_observed(distn[[z]], 
                                                              K[[z]], D[[z]])))
  npar = npar.hidden + npar.observed
  
  if (train.weight == 0) 
    train = NULL
  if (train.weight < 0) 
    stop("negative 'train.weight' is not valid")
  if (is.null(train)) {
    train.raw = NULL
    bicsub = npar * log(N)} else {
    train.raw = train
    if (length(train) != N) 
      stop("'train' must be NULL or of same length as data")
    if (!all(train %in% (0:(K0 + 1)))) 
      stop("'train' values must be integers between 0 and K0+1 inclusive")
    train[(train == 0)] = K0 + 1
    weight = train.weight
    Ktrain = K0 + 1
    Tprime = diag(Ktrain)[train, ]
    istrain = (train <= K0)
    wTprime = Tprime
    wTprime[istrain, ] = weight * wTprime[istrain, ]
    ltprob = log(colSums(wTprime)) - log(sum(wTprime))
    tprob = exp(ltprob)
    rpbase = diag(K0)
    lrpbase = log(rpbase)
    Tpdf = t(t(Tprime) * tprob)
    twmult = rep.int(1, N) + (weight - 1) * istrain
    Neff = sum(twmult)
    ltwcoef = log(Neff) - log(N)
    twcoef = exp(ltwcoef)
    train = namedList(train, weight, Ktrain, Tprime, istrain, 
                      tprob, ltprob, rpbase, lrpbase, Tpdf, twmult, Neff, 
                      twcoef, ltwcoef)
    bicsub = npar * log(Neff)
  }
  estep <- switch(topology, layered = .estep_mdmixmod_layered(), 
                  chained = .estep_mdmixmod_chained())
  mstep <- switch(topology, layered = .mstep_mdmixmod_layered(distn, 
                                                              train), chained = .mstep_mdmixmod_chained(distn, train))
  pstep <- switch(topology, layered = .pstep_mdmixmod_layered(distn, 
                                                              train), chained = .pstep_mdmixmod_chained(distn, train))
  sstep <- .sstep_mdmixmod(bicsub, train)
  if (prefit) {
    marginals = lapply(zvec, function(z) mixmod(X[[z]], K[[z]], 
                                                family = family[[z]], train = train.raw, train.weight = train.weight, 
                                                decreasing = decreasing[[z]], iter.max = iter.max, 
                                                dname = zvec[[z]]))
    W = lapply(marginals, posterior)
    U = switch(topology, layered = listMean(lapply(zvec, 
                                                   function(z) if ((K[[z]] == K0) && is.null(prior)) W[[z]] else posterior(mixmod(X[[z]], 
                                                                                                                                  K0, family = family[[z]], prior = prior, iter.max = iter.max, 
                                                                                                                                  train = train.raw, train.weight = train.weight)))), 
               chained = if ((K[[1]] == K0) && is.null(prior)) W[[1]] else posterior(mixmod(X[[1]], 
                                                                                            K0, family = family[[1]], prior = prior, iter.max = iter.max, 
                                                                                            train = train.raw, train.weight = train.weight)))
  }
  else {
    marginals = NULL
    W = mapply(.init_weights, X, K, distn, SIMPLIFY = FALSE)
    U = switch(topology, layered = listMean(lapply(zvec, 
                                                   function(z) .init_weights(X[[z]], K0, distn[[z]], 
                                                                             prior, decreasing = decreasing[[z]]))), chained = .init_weights(X[[1]], 
                                                                                                                                             K0, distn[[1]], prior, decreasing = decreasing[[1]]))
  }
  V = switch(topology, layered = lapply(W, function(Wz) .vgenDep(U, 
                                                                 Wz)), chained = c(list(.vgenDep(U, W[[1]])), lapply(2:Z, 
                                                                                                                     function(z) .vgenDep(W[[z - 1]], W[[z]]))))
  names(V) = zvec
  weights = namedList(U, V, W)
  params = mstep(X, weights, K, kvec, K0, k0vec, zvec, prior = prior)
  pdfs = pstep(X, params, kvec, zvec)
  stats = sstep(0, weights, params, pdfs)
  iteration.params = list(params)
  iteration.stats = list(stats)
  old_llik = stats[["llik"]] * 2
  iter = 1
  while (abs(stats[["llik"]]/old_llik - 1) > LC_ITER_TOL && 
         iter <= iter.max) {
    old_llik = stats[["llik"]]
    weights = estep(pdfs, params, zvec, kvec, k0vec)
    params = mstep(X, weights, K, kvec, K0, k0vec, zvec, 
                   params, prior = prior)
    pdfs = pstep(X, params, kvec, zvec)
    stats = sstep(iter, weights, params, pdfs)
    iter = iter + 1
    iteration.params[[iter]] = params
    iteration.stats[[iter]] = stats
  }
  weights = estep(pdfs, params, zvec, kvec, k0vec)
  posterior = weights$U
  rownames(posterior) = if (isMDF(X[[1]])) 
    rownames(X[[1]])
  else names(X[[1]])
  assignment = apply(posterior, 1, which.max)
  iteration.stats = Reduce(rbind, iteration.stats)
  rownames(iteration.stats) = NULL
  iteration.stats = as.data.frame(iteration.stats)
  retn = namedList(N, Z, D, K, K0, X, decreasing, npar, npar.hidden, 
                   npar.observed, train, iter, params, stats, weights, pdfs, 
                   posterior, assignment, iteration.params, iteration.stats, 
                   topology, family, distn, prefit, iter.max, dname, dattr, 
                   zvec, kvec, k0vec, prior, marginals, estep, mstep, pstep, 
                   sstep)
  class(retn) = c("mdmixmod", "mixmod")
  return(retn)
}

getLCMIXmodels.semi.only <- function(feature.table, is.train, test.labels, K0 = 2, TW = 10){
  
  BINARY_FEATURES2 = BINARY_FEATURES[BINARY_FEATURES %in% names(feature.table)]
  INTEGER_FEATURES2 = INTEGER_FEATURES[INTEGER_FEATURES %in% names(feature.table)]
  OPEN_FEATURES2 = OPEN_FEATURES[OPEN_FEATURES %in% names(feature.table)]
  CLOSED_FEATURES2 = CLOSED_FEATURES[CLOSED_FEATURES %in% names(feature.table)]

  # Make feature tables
  X = list(
    bern = do.call(cbind, feature.table[,BINARY_FEATURES2]),
    nbin = do.call(cbind, list(feature.table[,INTEGER_FEATURES2])),
    norm = cbind(
      do.call(cbind, feature.table[,OPEN_FEATURES2]),
      log(do.call(cbind, feature.table[,CLOSED_FEATURES2]))
    )
  )
  
  # Marginal model selection
  message(notification("marginal model selection", 1))
  margsel = msapply(namedList(2,3,4,5), function(K)
    mapply(X, names(X), FUN=function(x, fam)
      iclbic(mixmod(x, K=K, family=fam), map=T)))
  K = apply(margsel, 1, which.max) + 1
  
  # Joint model fitting
  message(notification("fitting joint models", 1))
  train = as.numeric(is.train)
  fm1 = mdmixmod(X, K=K, K0=K0, family=names(X),
                            train=train, train.weight=TW)
  fits = list(fm1)
  names(fits) = "semisup"
  
  # Make predictions for semisup and unsup models
  preds = lapply(fits, function(x) predict(x, type="prob")[!is.train, 1])
  preds.train = lapply(fits, function(x) predict(x, type="prob")[is.train, 1])
  preds.all = lapply(fits, function(x) predict(x, type="prob")[, 1])
  # Report performance
  reports = c("precis", "recall", "f1meas")
  suppressWarnings(rept1 <- as.data.frame(sapply(preds, perfmeas, 
                                                 labels=test.labels)[reports,]))
  rept1$report = reports
  rept1$FDR = 0
  suppressWarnings(rept2 <- as.data.frame(sapply(preds, perfmeas, lfdr=.20, 
                                                 labels=test.labels)[reports,]))
  rept2$report = reports
  rept2$FDR = .2
  suppressWarnings(rept3 <- as.data.frame(sapply(preds, perfmeas, lfdr=.01, 
                                                 labels=test.labels)[reports,]))
  rept3$report = reports
  rept3$FDR = .01
  combined.report = rbind(melt(rept1, id.vars = c('report', 'FDR')), 
                               
                          melt(rept2, id.vars = c('report', 'FDR')), 
                               
                          melt(rept3, id.vars = c('report', 'FDR') )
                               )
  
  # Return prediction probabilities and performance report
  return(list(preds = preds, performance = combined.report,fit=fm1,preds_train=preds.train,preds_all=preds.all))
  
}

