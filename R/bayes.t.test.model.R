  # THE MODEL.
original.model <- function() {
  modelString = "
  data {
    int<lower=1> Ntotal ;
    int x[Ntotal] ;
    real y[Ntotal] ;
    real meanY ;
    real sdY ;
  }
  transformed data {
    real unifLo ;
    real unifHi ;
    real normalSigma ;
    real expLambda ;
    unifLo <- sdY/1000 ;
    unifHi <- sdY*1000 ;
    normalSigma <- sdY*100 ;
    //expLambda <- 1/29.0 ;
    expLambda <- 1/(meanY-1) ;
  }
  parameters {
    real<lower=0> nuMinusOne ; 
    real mu[2] ;               // 2 groups
    real<lower=0> sigma[2] ;   // 2 groups
  }
  transformed parameters {
    real<lower=0> nu ;         // actually lower=1
    nu <- nuMinusOne + 1 ;
  }
  model {
    sigma ~ uniform( unifLo , unifHi ) ; // vectorized
    mu ~ normal( meanY , normalSigma ) ; // vectorized
    nuMinusOne ~ exponential( expLambda ) ;
    for ( i in 1:Ntotal ) {
      y[i] ~ student_t( nu , mu[x[i]] , sigma[x[i]] ) ;
    }
  }
  " # close quote for modelString
  modelString
}

modified.model <- function() {
  modelString = "
  data {
    int<lower=1> Ntotal ;
    int x[Ntotal] ;
    real y[Ntotal] ;
    real meanY ;
    real sdY ;
  }
  transformed data {
    real unifLo ;
    real unifHi ;
    real normalSigma ;
    real expLambda ;
    unifLo <- sdY/1000 ;
    unifHi <- sdY*1000 ;
    normalSigma <- sdY*100 ;
    expLambda <- 1/(meanY-1) ;
  }
  parameters {
    real<lower=0> nuMinusOne[2] ; // 2 groups, added by LiXC
    real mu[2] ;               // 2 groups
    real<lower=0> sigma[2] ;   // 2 groups
  }
  transformed parameters {
    real<lower=0> nu[2] ;         // actually lower=1, 2 groups, added by LiXC
    nu[1] <- nuMinusOne[1] + 1 ;
    nu[2] <- nuMinusOne[2] + 1 ;
  }
  model {
    sigma ~ uniform( unifLo , unifHi ) ; // vectorized
    mu ~ normal( meanY , normalSigma ) ; // vectorized
    nuMinusOne ~ exponential( expLambda ) ;
    for ( i in 1:Ntotal ) {
      y[i] ~ student_t( nu[x[i]] , mu[x[i]] , sigma[x[i]] ) ;
    }
  }
  " # close quote for modelString
  modelString
}

student_t_distribution_model <- function() {
  modelString = "
  data {
    int<lower=1> Ntotal ;
    real y[Ntotal] ;
    real meanY ;
    real sdY ;
  }
  transformed data {
    real unifLo ;
    real unifHi ;
    real normalSigma ;
    real expLambda ;
    unifLo <- sdY/1000 ;
    unifHi <- sdY*1000 ;
    normalSigma <- sdY*100 ;
    expLambda <- 1/(meanY - 1.0) ;
  }
  parameters {
    real<lower=0> nuMinusOne ; 
    real mu ;
    real<lower=0> sigma ;
  }
  transformed parameters {
    real<lower=0> nu ;  // actually lower=1
    nu <- nuMinusOne + 1 ;
  }
  model {
    sigma ~ uniform( unifLo , unifHi ) ;
    mu ~ normal( meanY , normalSigma ) ; 
    nuMinusOne ~ exponential( expLambda ) ;
    y ~ student_t( nu , mu , sigma ) ;
  }
  " # close quote for modelString
modelString
}
