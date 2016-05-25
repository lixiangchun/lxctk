
bayes.t.test <- function(vals, grp, stanDso, ...)
{
	Ny <- length(vals)
	Nx <- length(grp)

	if ( any( grp!=1 & grp!=2 ) ) { stop("All x values must be 1 or 2.") }
	if ( any( !is.finite(vals) ) ) { stop("All y values must be finite.") }
	if ( Nx != Ny ) { stop("vals and grp must be same length.") }

	dataList = list(y = vals, x = grp, Ntotal = Nx, meanY = mean(vals), sdY = sd(vals))

	parameters = c( "mu" , "sigma" , "nu" )     # The parameters to be monitored

	# Get MC sample of posterior:
	stanFit <- sampling(object=stanDso, data=dataList, pars=parameters, ...)
	#stanFit <- sampling( object=stanDso, data = dataList, 
	#                     pars = parameters, # optional
	#                     chains = nChains, iter = 4000, warmup = 2000, thin = thinSteps )
	#
	# For consistency with JAGS-oriented functions in DBDA2E collection, 
	# convert stan format to coda format:
	##codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
	##                                function(x) { mcmc(as.array(stanFit)[,x,]) } ) )

	stanFit
}

