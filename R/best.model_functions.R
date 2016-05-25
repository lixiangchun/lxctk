
stanfit2mcmc.list <- function(stanFit)
{
	codaObject <- coda::mcmc.list(lapply(1:ncol(stanFit), function(x) {coda::mcmc(as.array(stanFit)[,x,])}))
	codaObject
}

best.student_t <- function(y, stanDso, meanY=mean(y), sdY=sd(y), 
							unifLo=sdY/1000, unifHi=sdY*1000, 
							normalSigma=sdY*100,
							expLambda=1/abs(meanY-1.0), ...)
{
	if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
	dataList = list(y = y, Ntotal = length(y), meanY = meanY, sdY = sdY,
					unifLo=unifLo, unifHi=unifHi,normalSigma=normalSigma,expLambda=expLambda)
	parameters = c( "mu" , "sigma" , "nu" )     # The parameters to be monitored
	# Get MC sample of posterior:
	stanFit <- sampling(object=stanDso, data=dataList, pars=parameters, ...)
	stanFit
}

best.robust_t_test <- function(y, grp, stanDso, meanY=mean(y), sdY=sd(y),
								unifLo=sdY/1000, unifHi=sdY*1000, 
								normalSigma=sdY*100,
								expLambda=1/abs(meanY-1.0),
								parameters=c('mu', 'sigma', 'nu', 'mu_diff', 'sigma_diff', 'nu_diff'), ...)
{
	Ny <- length(y)
	Nx <- length(grp)
	if ( any( grp!=1 & grp!=2 ) ) { stop("All x values must be 1 or 2.") }
	if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
	if ( Nx != Ny ) { stop("y and grp must be same length.") }
	dataList = list(y = y, x = grp, Ntotal = Nx, meanY = meanY, sdY = sdY,
					unifLo=unifLo, unifHi=unifHi, normalSigma=normalSigma, expLambda=expLambda)
	##parameters = c( "mu" , "sigma" , "nu" )     # The parameters to be monitored
	# Get MC sample of posterior:
	stanFit <- sampling(object=stanDso, data=dataList, pars=parameters, ...)
	stanFit
}

#load('x.RData')
#options(mc.cores = parallel::detectCores())
#r = best.student_t(x, stanDso, unifLo=0, unifHi=.Machine$double.xmax, iter=20000)
#r = best.student_t(x, stanDso)#unifLo=0, unifHi=.Machine$double.xmax, iter=20000)


