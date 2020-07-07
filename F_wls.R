## F_wls.R

####################################
##			 wls function 		  ##
####################################
## input:
## 		design matrix x
##		y
##		number of slices (same as Slice Inverse Regression)

## output:
##		weighted leverage score for every variable
##		higher the score, more important the variable is.
####################################
wlsFunc <- function(x, y, nslice=10, cn1=0.1, cn2=1, choose.dir=FALSE){
	## require package MASS and dr
	require(MASS)
	require(dr)

	## basic information about x and y
	n <- dim(x)[1]
	p <- dim(x)[2]
	h <- nslice
	x = scale(x,scale=FALSE)

	## slice
	if( is.factor(y) == 1 ){
		index <- as.numeric(factor(y))
		nh <- summary(factor(y))
		h <- length(nh)
	}
	if( is.factor(y) != 1 ){
		slice <- dr.slices.arc(y,h)
		index <- slice$slice.indicator		# Slice Index
		nh <- slice$slice.sizes				# Observations per Slice
	}
	ph <- nh/n

	####################################
	## calculate wls
	####################################
	wls <- c()								# Weighted Leverage Score

	####################################
	## scenario 1
	####################################
	svdx <- svd(x)
	u <- svdx$u
	d <- svdx$d
	v <- svdx$v

	if( choose.dir==FALSE) {
		dir = min(n, p)
	} else {
		## selection of d
		theta = d^2/(d[1])^2 + 1
		loglik <- penalty <- rep(0, length(d) )
		for( i in 1:length(d) ){
		  if(i < length(d)) {
				loglik[i] <- sum( log(theta[(i+1):length(d)]) + 1 - theta[(i+1):length(d)] )
			} else loglik[i] = 0
			penalty[i] <- i*cn1/sqrt(n)
		}
		BIC = -loglik + penalty
		(dir <- which.min(BIC))
		print(dir)
	}

	## calculate WLS
	w <- matrix(ncol = dir, nrow=length(nh))
	for(j in 1:dir){
  		for(i in 1:length(nh)) w[i,j] <- sum(u[,j] * (index == i))/nh[i]
	}

	uut <- array(dim = c(length(nh), dir, dir))			# UUT Array
	for(i in 1:length(nh)) uut[i,,] <- (nh[i]) * ( w[i,] %*% t(w[i,]) )

	## LEVERAGE SCORES
	for(j in 1:p) wls[j] <- t(v[j,1:dir]) %*% colSums(uut) %*% v[j,1:dir]


	####################################
	## BIC
	####################################
	wls.sort = sort(wls, decreasing = TRUE)

	loglik <- penalty <- rep(0, min(n,p))
	for(k in 1:min(n,p)){
		temp_loglik = sum(wls.sort[1:k])
		(loglik[k] = -log(temp_loglik))
		penalty[k] = (log(n) + cn2*log(p))*k/max(n,p)
	}
	BIC <- loglik + penalty
	(sel.k <- which.min(BIC))

	select <- order(wls, decreasing = T)[1:sel.k]
	return( list(wls=wls, select=select ) )
}
