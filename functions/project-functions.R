require(mvtnorm)
require(numDeriv)
source("ghyp-functions.R")

#' generate from a truncated normal distribution with lower bound a
#' @param mu,sig mean and standard deviation of the original Normal
#' @param a lower bound
rtnorm <- function(n,a, mu, sig) {
  qnorm(runif(n,pnorm(a,mean=mu,sd=sig),1),mean=mu,sd=sig)
}


# log P(y, V|theta, X, W) = log P(y|V, theta, X, W) + log P(V|theta, X, W)
gh.loglik<-function(y, v, mu, alpha, sigma, psi, lambda){
  ll<-dgig(v, psi=psi, eta=1, lambda=lambda, log = TRUE)
  ll<-ll+dnorm(y, mean=mu, sd=sigma, log=TRUE)
  sum(ll)
}


#' Statistics to plot analytic posterior using a 2-D grid method
#' @param seq1 grid of var1 values
#' @param seq2 grid of var2 values
#' @param lpmat matrix of joint log posterior evaluations of var1 and var2
#' @return a list of grid values, densities, mean, and credible interval
grid.plot.stat <- function(seq1, seq2, lpmat, alpha = .95) {
  # calculate margins
  dens12 <- exp(lpmat-max(lpmat)) # unnormalized p(var1, var2 | v)
  dens1 <- rowSums(dens12) # p(var1 | v)
  d1 <- seq1[2]-seq1[1]
  dens1 <- dens1/sum(dens1)/d1 # normalize
  mean1 <- sum(seq1*dens1)*d1 # posterior mean
  # credible interval
  clim <- c((1-alpha)/2, 1 - (1-alpha)/2)
  L1 <- seq1[which.min(abs(cumsum(dens1)*d1 - clim[1]))]
  U1 <- seq1[which.min(abs(cumsum(dens1)*d1 - clim[2]))]
  dens2 <- colSums(dens12) # p(var2 | y)
  d2 <- seq2[2]-seq2[1]
  dens2 <- dens2/sum(dens2)/d2
  mean2 <- sum(seq2*dens2)*d2
  L2 <- seq2[which.min(abs(cumsum(dens2)*d2 - clim[1]))]
  U2 <- seq2[which.min(abs(cumsum(dens2)*d2 - clim[2]))]
  list(seq1 = seq1, dens1 = dens1, mean1 = mean1, pCI = c(L1, U1),
       seq2 = seq2, dens2 = dens2, mean2 = mean2, lCI = c(L2, U2),
       dens12 = dens12)
}


# log P(beta, alpha | y, X, W, V, gamma, psi, lambda)
ba.loglik<-function(y, X, W, V, beta, alpha, gamma){
  sigma2<-exp(W%*%gamma)*V
  del<-as.numeric(1/sigma2)
  X_tilde<-cbind(X, V)
  XDX<-crossprod(X_tilde, X_tilde*del)
  XDy<-crossprod(X_tilde, y*del)
  beta.hat<-solve(XDX, XDy)
  Sig.hat<-solve(XDX)
  dmvnorm(c(beta, alpha), mean=beta.hat, sigma=Sig.hat, log=TRUE)
}


# Gibbs sampler for c(alpha, beta) with a flat prior. 
# c(beta, alpha) ~ N(c(beta.hat, alpha.hat), sigma.hat^2), where c(beta.hat, alpha.hat) is the MLE of 
# the regression model yi | xi, wi, Vi ~ N(c(xi, Vi)%*%beta, exp(wi%*%gamma)*Vi) 
ba.gibbs<-function(y, X, W, V, gamma){
  delta<-(exp(W%*%gamma)*V)^(-1)
  M<-lm(y~cbind(X,V)-1,weights=delta)
  as.numeric(rmvnorm(1, mean=M$coefficients, sigma=vcov(M)))
}


# log P(gamma | y, X, W, V, beta, alpha, psi, lambda)
gamma.loglik<-function(y, X, W, V, beta, alpha, gamma){
  R<-(y-X%*%beta-V*alpha)^2/V
  -0.5*sum(R/exp(W%*%gamma)+W%*%gamma)
}


#' Returns the mean and covariance matrix of the MVN proposal matching the mode and quadrature of \code{gamma.loglik}
gam.mv<-function(y, X, W, V, beta, alpha){
  R<-(y-X%*%beta-V*alpha)^2/V
  gam.reg<-glm(R~W-1, family=Gamma("log"), control=list(maxit=50))
  prop.mode <- gam.reg$coefficients
  prop.v <- vcov(gam.reg)
  list(mode=prop.mode, v=prop.v)
}


#' MIID sampler of gamma. Proposal is a MVN matching the mode and quadrature of \code{gamma.loglik}
gamma.miid<-function(y, X, W, V, beta, alpha, old.gam){
  
  gam.prop.mv<-gam.mv(y, X, W, V, beta, alpha)
  gam.prop.mu<-gam.prop.mv$mode
  gam.prop.V<-gam.prop.mv$v
  new.gam<-t(rmvnorm(1,mean=gam.prop.mu,sigma=gam.prop.V))
  
  gam.lacc.old<-gamma.loglik(y, X, W, V, beta, alpha, old.gam)-
    dmvnorm(as.numeric(old.gam), mean=gam.prop.mu, sigma=gam.prop.V,log=TRUE)
  gam.lacc.new<-gamma.loglik(y, X, W, V, beta, alpha, new.gam)-
    dmvnorm(as.numeric(new.gam), mean=gam.prop.mu, sigma=gam.prop.V,log=TRUE)
  gam.lacc<-gam.lacc.new - gam.lacc.old
  
  if (gam.lacc>0 || runif(1) < exp(gam.lacc)) accept<-TRUE
  else {
    accept<-FALSE
    new.gam<-old.gam
  }
  
  list(gamma=new.gam, accept=accept)
}


# log P(V | y, X, W, beta, alpha, gamma, psi, lambda)
# Using a flat prior, the posterior is a GIG(psii.hat, etai.hat, lambdai.hat)
v.logpost<-function(y, X, W, V, beta, alpha, gamma, psi, lambda){
  Ai<-3/2-lambda
  Bi<-psi/2+alpha^2/(2*exp(W%*%gamma))
  Ci<-psi/2+(y-X%*%beta)^2/(2*exp(W%*%gamma))
  #-(Ai*log(V)+Bi*V+Ci/V)
  lambdai.hat<-1-Ai
  etai.hat<-sqrt(Ci/Bi)
  psii.hat<-sqrt(4*Bi*Ci)
  dgig(V, psii.hat, etai.hat, lambdai.hat, log = TRUE)
}


# Gibbs sampler for updating the missing V.
# P(V | y, X, W, beta, alpha, gamma, psi, lambda) = GIG(psii.hat, etai.hat, lambdai.hat)
v.gibbs<-function(y, X, W, beta, alpha, gamma, psi, lambda){
  n<-length(y)
  
  Ai<-3/2-lambda
  Bi<-psi/2+alpha^2/(2*exp(W%*%gamma))
  Ci<-psi/2+(y-X%*%beta)^2/(2*exp(W%*%gamma))
  lambdai.hat<-1-Ai
  etai.hat<-sqrt(Ci/Bi)
  psii.hat<-sqrt(4*Bi*Ci)
  rgig(n, psii.hat, etai.hat, lambdai.hat)
}


# log P(psi, lambda | y, X, W, V, beta, alpha, gamma)
psi.lam.loglik<-function(v, psi, lambda){
  n<-length(v)
  S<-sum(log(v))
  T<-sum(v+1/v)
  lambda*S-0.5*psi*T-n*log(2*besselK(psi, nu = lambda, expon.scaled = TRUE))+n*psi
}


# posterior draw from P(psi, lambda | y, X, W, V, beta, alpha, gamma)
psi.lam.mwg<-function(psi.lam.old, V, rwsd){
  accept<-c(FALSE, FALSE)
  names(accept)<-c("psi", "lambda")
  
  psl.lp.old<-psi.lam.loglik(V, psi.lam.old["psi"], psi.lam.old["lambda"])+dnorm(psi.lam.old["psi"],0,8,log=TRUE) #log(psi.lam.old["psi"])#

  for (name in c("psi", "lambda")){
    psi.lam.new<-psi.lam.old
    psi.lam.new[name]<-psi.lam.new[name]+rwsd[name]*rnorm(1)
    
    if (psi.lam.new["psi"]>0){
      psl.lp.new<-psi.lam.loglik(V, psi.lam.new["psi"], psi.lam.new["lambda"])+dnorm(psi.lam.new["psi"],0,8,log=TRUE)#log(psi.lam.new["psi"])#
      psl.lacc<-psl.lp.new-psl.lp.old
      
      if (psl.lacc>0 || runif(1)<exp(psl.lacc)){
        psi.lam.old<-psi.lam.new
        psl.lp.old<-psl.lp.new 
        accept[name]<-TRUE
      }
    }
  }
  
  list(psi.lam=psi.lam.old, accept=accept)
}

#' Gibbs sampler from the joint posterior
#' 
#' @param y vector of survival time
#' @param X mean covariate matrix
#' @param W log-variance covariate matrix
#' @param beta0 initial value of mean coefficients
#' @param gamma0 initial value of log-variance coefficients
#' @param psi0,lambda0 initial values of the parameters of V~GIG(psi,1,lambda)
#' @param V0 initial values of V sampled from rgig(888, psi0, 1, lambda0)
#' @param alpha0 initial value of the mixing coefficient alpha
#' @param censor.id a vector of ids of the censored cases
#' @param censored whether the data has censored cases. If true, in each mcmc cycle sample from a truncated normal to replace the censoring time
#' @param ba.fixed logical, whether or not beta and alpha are known
#' @param gam.fixed logical, whether or not gamma is known
#' @param psl.fixed logical, whether or not psi and lambda are known
#' @param v.fixed logical, whether or not V is known
#' @param v.out logical, whether to output V in every mcmc cycle. If true, add an nsamples*n matrix in the returned list
#' @param y.out logical, whether to output y in every mcmc cycle. If true, add an nsamples*n matrix in the returned list
#' @param rwsd jump size for psi and lambda in mwg, a named vector
#' @return a list of mcmc samples for each paramter and acceptance rate for gamma, psi and lambda, optionally a matrix of y and the missing V in each mcmc iteration
mcmc<-function(nsamples, y, X, W, beta0, alpha0, gamma0, V0, psi0, lambda0, censor.id=NULL, censored=FALSE,
                ba.fixed=FALSE, gam.fixed=FALSE, psl.fixed=FALSE, v.fixed=FALSE, v.out=FALSE, rwsd, y.out=FALSE) {
  burn<-min(floor(nsamples/10),1e3)
  paccept<-numeric(3)
  names(paccept)<-c("gamma","psi","lambda")
  
  n<-length(y)
  px<-ncol(X)
  if (is.null(ncol(X))) px<-1
  pw<-ncol(W)
  if (is.null(ncol(W))) pw<-1
  
  beta.alpha<-matrix(NA, nsamples, px+1)
  gamma<-matrix(NA, nsamples, pw)
  psi.lam<-matrix(NA, nsamples, 2)
  if (v.out) V<-matrix(NA, nsamples, n)
  if (y.out) Y<-matrix(NA, nsamples, n)
  
  beta.alpha.curr<-c(beta0, alpha0)
  gamma.curr<-as.matrix(gamma0)
  psi.lam.curr<-c("psi"=psi0, "lambda"=lambda0)
  V.curr<-V0
  y.curr<-y
  
  for (ii in (-burn+1):nsamples){
    
    if (censored){
      #update y for censored cases
      mu.curr<-(cbind(X,V.curr)%*%beta.alpha.curr)[,1]
      sig2.curr<-(exp(W%*%gamma.curr)*V.curr)[,1]
      y.curr[censor.id]<-sapply(censor.id, function(i){rtnorm(1, a=y[i], mu.curr[i], sig=sqrt(sig2.curr[i]))})
      if (y.out && ii>0) Y[ii,]<-y.curr
    }
  
    #update c(beta, alpha)
    if (!ba.fixed){
      beta.alpha.curr<-ba.gibbs(y.curr, X, W, V.curr, gamma.curr)
      if (ii>0) beta.alpha[ii, ]<-beta.alpha.curr
    }
    
    #update gamma
    if (!gam.fixed){
      gam.prop<-gamma.miid(y.curr, X, W, V.curr, as.matrix(beta.alpha.curr[1:px]), beta.alpha.curr[px+1], as.matrix(gamma.curr))
      if (gam.prop$accept){
        gamma.curr<-gam.prop$gamma
        paccept["gamma"]<-paccept["gamma"]+1
      }
      if (ii>0) gamma[ii,]<-gamma.curr
    }
    
    
    #update V
    if (!v.fixed){
      V.curr<-v.gibbs(y.curr, X, W, as.matrix(beta.alpha.curr[1:px]), beta.alpha.curr[px+1], gamma.curr, psi.lam.curr["psi"], psi.lam.curr["lambda"])
      if (v.out){
        if (ii>0) V[ii,]<-V.curr
      }
    }
    
    #update psi and lambda
    if (!psl.fixed){
      psi.lam.prop<-psi.lam.mwg(psi.lam.curr, V.curr, rwsd)
      psi.lam.curr<-psi.lam.prop$psi.lam
      paccept[c("psi", "lambda")]<-paccept[c("psi", "lambda")]+psi.lam.prop$accept
      if (ii>0) psi.lam[ii,]<-psi.lam.curr
    }
    
  }
  
  paccept <- paccept/(nsamples+burn) # acceptance rate
  if (!gam.fixed || !psl.fixed) message("MH Acceptance Rate:")
  if (!gam.fixed) message("gamma: ", round(paccept[1]*100), "%")
  if (!psl.fixed) {
    message("psi: ", round(paccept[2]*100), "%")
    message("lambda: ", round(paccept[3]*100), "%")
  }
  if (v.out && y.out){
    ans<-list(beta=beta.alpha[,1:px], alpha=beta.alpha[,px+1], gamma=gamma, psi=psi.lam[,1], lambda=psi.lam[,2], accept=paccept, V=V, Y=Y)
  } else if (y.out){
    ans<-list(beta=beta.alpha[,1:px], alpha=beta.alpha[,px+1], gamma=gamma, psi=psi.lam[,1], lambda=psi.lam[,2], accept=paccept, Y=Y)
  } else if (v.out){
    ans<-list(beta=beta.alpha[,1:px], alpha=beta.alpha[,px+1], gamma=gamma, psi=psi.lam[,1], lambda=psi.lam[,2], accept=paccept, V=V)
  }
  else ans<-list(beta=beta.alpha[,1:px], alpha=beta.alpha[,px+1], gamma=gamma, psi=psi.lam[,1], lambda=psi.lam[,2], accept=paccept)
  ans
}
