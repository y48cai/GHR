source("project-functions.R")
# experiment with simulated dataset
# simulate some data
n1<-500
X1<-rtnorm(n1, 10, 10, 5)
W1<-rtnorm(n1, 3, 3, 1)
beta1<-1
alpha1<-1
gamma1<-1
psi1<-20
lambda1<-20
V1<-rgig(n1, psi1, 1, lambda1) #V1 | X1, W1 ~ rgig(psi1, 1, lambda1)
mu1<-X1*beta1+alpha1*V1
sigma2<-exp(W1*gamma1)*V1
y1<-rnorm(n1, mu1, sqrt(sigma2)) # Y1 | V1, X1, W1 ~ N(mu1, sigma2)

#---------check likelihoods and analytic posteriors vs. mcmc samples----------
# ---------------------psi, lambda-------------------
# check unsimplified vs simplified log likelihoods
ntest<-10
Psi<-psi1+runif(ntest, 5, 10)
Lambda<-lambda1+runif(ntest, 5, 10)

llu<-sapply(1:ntest, function(ii){
  psi<-Psi[ii]
  lambda<-Lambda[ii]
  gh.loglik(y1, V1, mu1, alpha1, sqrt(sigma2), psi, lambda)
})

lls<-sapply(1:ntest, function(ii){
  psi.lam.loglik(V1, Psi[ii], Lambda[ii])
})

llu-lls

# check analytic posterior vs mcmc samples. Prior: psi~N(0,5)
# assuming beta, alpha, gamma, V are fixed
result.pl<-mcmc(20000, y1, X1, W1, beta1, alpha1, gamma1, V1, psi1, lambda1,
                ba.fixed=TRUE, gam.fixed=TRUE, psl.fixed=FALSE, v.fixed=TRUE, v.out=FALSE, rwsd=c("psi"=0.4, "lambda"=0.4))
colMeans(cbind(result.pl$psi, result.pl$lambda))

pseq<-seq(5,15, length.out=500)
lseq<-seq(5,15, length.out=500)
psl<-as.matrix(expand.grid(pseq, lseq))
lp.mat<-apply(psl, 1, function(theta){psi.lam.loglik(V1, theta[1], theta[2])}) #+dnorm(theta[1],0,5,log=TRUE)
lp.mat<-matrix(lp.mat, 500, 500)
psi.lam.stat<-grid.plot.stat(pseq, lseq, lp.mat, alpha = .95)

par(mfrow=c(1,2))
hist(result.pl$psi, breaks=100, prob=T, main=expression(p(psi*" | "*bold(y),bold(X),bold(W),alpha, beta, gamma)), xlab=expression(psi))
lines(psi.lam.stat$seq1, psi.lam.stat$dens1, type='l', col="red")
hist(result.pl$lambda, breaks=100, prob=T, main=expression(p(lambda*" | "*bold(y),bold(X),bold(W),alpha, beta, gamma)), xlab=expression(lambda))
lines(psi.lam.stat$seq2, psi.lam.stat$dens2, type='l', col="red")

#-------------------beta, alpha------------------
# check unsimplified vs simplified log likelihoods
Beta<-beta1+runif(ntest, 1, 2)
Alpha<-alpha1+runif(ntest, 1, 2)

llu<-sapply(1:ntest, function(ii){
  beta<-Beta[ii]
  alpha<-Alpha[ii]
  gh.loglik(y1, V1, beta*X1+alpha*V1, alpha, sqrt(sigma2), psi1, lambda1)
})

lls<-sapply(1:ntest, function(ii){
  ba.loglik(y1, X1, W1, V1, Beta[ii], Alpha[ii], as.matrix(gamma1))
})

llu-lls

# check analytic posterior vs mcmc samples, assuming V, gamma, psi, lambda are fixed
result.ba<-mcmc(10000, y1, X1, W1, beta1, alpha1, gamma1, V1, psi1, lambda1,
                ba.fixed=FALSE, gam.fixed=TRUE, psl.fixed=TRUE, v.fixed=TRUE, v.out=FALSE, rwsd=c("psi"=0.4, "lambda"=0.4))
colMeans(cbind(result.ba$beta, result.ba$alpha))

bseq<-seq(0.1,2, length.out=500)
aseq<-seq(-2,2, length.out=500)
ba<-as.matrix(expand.grid(bseq, aseq))
lp.mat<-apply(ba, 1, function(theta){ba.loglik(y1, X1, W1, V1, theta[1], theta[2], as.matrix(gamma1))})
lp.mat<-matrix(lp.mat, 500, 500)
ba.stat<-grid.plot.stat(bseq, aseq, lp.mat, alpha = .95)

par(mfrow=c(1,2))
hist(result.ba$beta, breaks=100, prob=T, main=expression(p(beta*" | "*bold(y),bold(X),bold(W),gamma,psi,lambda)), xlab=expression(beta))
lines(ba.stat$seq1, ba.stat$dens1, type='l', col="red")
hist(result.ba$alpha, breaks=100, prob=T, main=expression(p(alpha*" | "*bold(y),bold(X),bold(W),gamma,psi,lambda)), xlab=expression(alpha))
lines(ba.stat$seq2, ba.stat$dens2, type='l', col="red")

#---------------------gamma------------------------
# check unsimplified vs simplified log likelihoods
Gamma<-gamma1+runif(ntest,1,2)

llu<-sapply(1:ntest, function(ii){
  gamma<-Gamma[ii]
  sigma<-sqrt(exp(W1*gamma)*V1)
  gh.loglik(y1, V1, beta1*X1+alpha1*V1, alpha1, sigma, psi1, lambda1)
})

lls<-sapply(1:ntest, function(ii){
  gamma.loglik(y1, X1, W1, V1, as.matrix(beta1), alpha1, as.matrix(Gamma[ii]))
})

llu-lls

# compare target (true) to proposal (approximation)
par(mfrow=c(1,1))
gmv<-gam.mv(y1, X1, W1, V1, as.matrix(beta1), alpha1)
gseq<-seq(0.8,1.2, length.out=500)
lp.dens<-sapply(1:length(gseq), function(ii){gamma.loglik(y1, X1, W1, V1, as.matrix(beta1), alpha1, as.matrix(gseq[ii]))})
gdens<-exp(lp.dens-max(lp.dens))
dg <- gseq[2]-gseq[1]
gdens <- gdens/sum(gdens)/dg
plot(gseq, gdens, type='l', lwd=1.5, xlab=expression(gamma), ylab="Density")
curve(dnorm(x, gmv$mode, sqrt(gmv$v)), add = TRUE, col = "red", lwd=1.5)
legend("topleft", legend = c("True Posterior", "Normal Approx."), fill = c("black", "red"))

# check analytic posterior vs mcmc samples. Assuming V, alpha, beta, psi, lambda fixed
result.gam<-mcmc(10000, y1, X1, W1, beta1, alpha1, gamma1, V1, psi1, lambda1,
                 ba.fixed=TRUE, gam.fixed=FALSE, psl.fixed=TRUE, v.fixed=TRUE, v.out=FALSE, rwsd=c("psi"=0.4, "lambda"=0.4))
mean(result.gam$gamma)

hist(result.gam$gamma, breaks=100, prob=T, main=expression(p(gamma*" | "*bold(y),bold(X),bold(W),beta,alpha,psi,lambda)), xlab=expression(gamma))
lines(gseq, gdens, col="red")
#----------------------missing data (V)-----------------------------
# check unsimplified vs simplified log likelihoods
# only V varies, everything else fixed 
X.sim<-rep(X1[1], ntest)
W.sim<-rep(W1[1], ntest)
V.sim<-rgig(ntest, psi1, 1, lambda1)
mu.sim<-X.sim*beta1+V.sim*alpha1
sig2.sim<-exp(W.sim*gamma1)*V.sim
#y.sim<-rnorm(ntest, mean=mu.sim, sd=sqrt(sig2.sim))
y.sim<-rep(y1[1], ntest)

ldu<-sapply(1:ntest, function(ii){
  gh.loglik(y=y.sim[ii], V.sim[ii], mu=mu.sim[ii], alpha=alpha1, sigma=sqrt(sig2.sim)[ii], psi=psi1, lambda=lambda1)
})

lds<-sapply(1:ntest, function(ii){v.logpost(y.sim[ii], X.sim[ii], W.sim[ii], V.sim[ii], beta1, alpha1, gamma1, psi1, lambda1)})

ldu-lds


# pick a few Vi's to check mcmc samples matches analytic log posterior
result.v<-mcmc(10000, y1[1:20], X1[1:20], W1[1:20], beta1, alpha1, gamma1, V1, psi1, lambda1,
               ba.fixed=TRUE, gam.fixed=TRUE, psl.fixed=TRUE, v.fixed=FALSE, v.out=TRUE, rwsd=c("psi"=0.4, "lambda"=0.4))


Vind <- sort(sample(1:20, 6)) # pick a few Vi's at random
par(mfrow = c(2,3))
invisible(sapply(Vind, function(vind) {
  hist(result.v$V[,vind], breaks = 100, freq = FALSE,
       xlab = parse(text = paste0("v[", vind, "]")),
       main = parse(text = paste0("p(v[", vind,
                                  "]*\" | \"*theta,bold(D))")))
  vSeq <- range(result.v$V[,vind])
  vSeq <- seq(vSeq[1], vSeq[2], len = 100)
  vPdf <- v.logpost(y = y1[vind], X=X1[vind], W=W1[vind], V = vSeq, 
                    beta = beta1, alpha = alpha1, gamma = gamma1, psi = psi1, lambda=lambda1)
  vPdf <- exp(vPdf - max(vPdf))
  vPdf <- vPdf/sum(vPdf)/(vSeq[2]-vSeq[1])
  lines(vSeq, vPdf, col = "red")
}))

#-------------------sample from joint posterior-------------------------------
result1<-mcmc(10000, y1, X1, W1, beta0=5, alpha0=5, gamma0=gamma1, V0=V1, psi0=5, lambda0=5,
              ba.fixed=FALSE, gam.fixed=FALSE, psl.fixed=FALSE, v.fixed=FALSE, v.out=FALSE, rwsd=c("psi"=0.1, "lambda"=0.3))

# point estimates and credible intervals
post.mean<-c(mean(result1$beta), mean(result1$alpha), colMeans(result1$gamma), mean(result1$psi), mean(result1$lambda))
names(post.mean)<-c("beta", "alpha", "gamma", "psi", "lambda")
sapply(1:5, function(i){quantile(result1[[i]], c(0.025,0.975))})

# MCMC histogram, true parameter values vs posterior means
theta.true<-c(beta1, alpha1, gamma1, psi1, lambda1)
par(mfrow = c(2,3))
invisible(lapply(1:5, function(i){
  hist(result1[[i]], prob=T, breaks=100, xlab=parse(text=names(result1)[i]),
       main=parse(text = paste0("p(", names(result1)[i], ")")))
  abline(v=theta.true[i], col="blue")
  abline(v=post.mean[i], col="red")
}))
plot(-10:0, xlim=c(1,10), ylim=c(1,10), xlab="", ylab="", xaxt='n', yaxt='n', axes=FALSE)
legend("topleft", legend = c("True Parameter", "Posterior Mean"), fill = c("blue", "red"), cex=1.5)

# another MCMC run, this time with a prior imposed on psi.
result2<-mcmc(10000, y1, X1, W1, beta0=5, alpha0=5, gamma0=gamma1, V0=V1, psi0=5, lambda0=5,
              ba.fixed=FALSE, gam.fixed=FALSE, psl.fixed=FALSE, v.fixed=TRUE, v.out=FALSE, rwsd=c("psi"=0.4, "lambda"=0.4))

# compare the true pdf of y_i|x_i, w_i, theta0 vs the pdf of y_i|x_i, w_i, theta,
# where theta is the point estimates of the previous two mcmc runs
par(mfrow=c(2,3))
invisible(sapply(yind, function(i) {
  plot(x=yseq, y=dghyp(yseq, X1[i]*beta1, alpha1, sqrt(exp(W1[i]*gamma1)), psi1, lambda1),type='l', 
       xlab = parse(text = paste0("y", i)), ylab="Density")
  lines(yseq, dghyp(yseq, X1[i]*post.mean["beta"], post.mean["alpha"], sqrt(exp(W1[i]*post.mean["gamma"])),
                    post.mean["psi"], post.mean["lambda"]), type='l', ylab="Density", ylim=c(0,0.03),col="blue")
  lines(yseq, dghyp(yseq, X1[i]*post.mean2["beta"], post.mean2["alpha"], sqrt(exp(W1[i]*post.mean2["gamma"])),
                    post.mean2["psi"], post.mean2["lambda"]),col="red")
}))
plot(-10:0, xlim=c(1,10), ylim=c(1,10), xlab="", ylab="", xaxt='n', yaxt='n', axes=FALSE)
legend("topleft", legend = c(expression(p(y[i]*" | "*theta[0],x[i],W[i])), expression(p(y*" | "*theta[1],x[i],W[i])), expression(p(y*" | "*theta[2],x[i],w[i]))), 
       fill = c("black","blue", "red"), cex=1.5)

#-------------------------------------------------------------
require("survival")
source("ghyp-functions.R")
source("hlm-functions.R")

# data analysis with colon dataset
# preprocess dataset
colon1<-colon[colon$etype==2,-2]
rownames(colon1)<-seq(length=nrow(colon1))
colon1$sex<-as.factor(colon1$sex)
colon1$obstruct<-as.factor(colon1$obstruct)
colon1$perfor<-as.factor(colon1$perfor)
colon1$adhere<-as.factor(colon1$adhere)
colon1$status<-as.factor(colon1$status)
colon1$differ<-as.factor(colon1$differ)
colon1$extent<-as.factor(colon1$extent)
colon1$surg<-as.factor(colon1$surg)
colon1$node4<-as.factor(colon1$node4)
colon1$etype<-as.factor(colon1$etype)
colon1$nodes<-as.integer(colon1$nodes)

# identify observations with missing covariates
colon1$miss.nodes<-FALSE
colon1$miss.nodes[is.na(colon1$nodes)]<-TRUE
colon1$miss.differ<-FALSE
colon1$miss.differ[is.na(colon1$differ)]<-TRUE

miss.ids<-colon1[c(colon1$id[colon1$miss.nodes],colon1$id[colon1$miss.differ]),"id"]
non.miss.ids<-colon1$id[-miss.ids]
colon.non.miss<-colon1[non.miss.ids,]
censor.id<-colon.non.miss[colon.non.miss$status==0,]$id
censor.id<-as.character(censor.id)
noncensor.id<-colon.non.miss[colon.non.miss$status==1,]$id

summary(lm(time~rx+sex+age+obstruct+perfor+adhere+nodes+differ+extent+surg+node4-1,data=colon1)) 
M.lm<-lm(time~rx+sex+age+obstruct+perfor+adhere+nodes+differ+extent+surg+node4-1,data=colon1)
summary(glm((M.lm$residuals)^2~model.matrix(M.lm)-1,family=Gamma("log"))) 

# Determine significant covariates using an HLM setting, and fit an HLM as the 
# starting values for beta and hat in MCMC.
X<-model.matrix(lm(time~rx+age+obstruct+differ+node4, data=colon1[non.miss.ids,]))
W<-model.matrix(lm(time~rx+obstruct+node4, data=colon1[non.miss.ids,]))
y<-colon1$time[non.miss.ids]
names(y)<-non.miss.ids
M.hlm<-hlm.fit(y=y,X,W)
beta0<-M.hlm$beta
gamma0<-M.hlm$gamma
alpha0<-10
psi0<-10
lambda0<-10
V0<-rgig(nrow(X),psi=psi0, eta=1, lambda=lambda0)

# run MCMC
result<-mcmc(nsamples=10000, y=y, X=X, W=W, beta0=beta0, alpha0=alpha0, gamma0=gamma0, V0=V0, psi0=psi0, lambda0=lambda0, 
             censor.id=censor.id, censored=TRUE, rwsd=c("psi"=0.2, "lambda"=0.2))

# posterior sample means - point estimate for theta
beta.hat<-colMeans(result$beta)
alpha.hat<-mean(result$alpha)
gamma.hat<-colMeans(result$gamma)
psi.hat<-mean(result$psi)
lambda.hat<-mean(result$lambda) 

mu.hat<-X%*%beta.hat
sigma.hat<-sqrt(exp(W%*%gamma.hat))

# confidence intervals for parameters
beta.CI<-sapply(1:ncol(X), function(i) quantile(x=(result$beta)[,i], probs=c(0.025, 0.975)))
colnames(beta.CI)<-names(beta0)
beta.est<-rbind(beta.CI, beta.hat)

gamma.CI<-sapply(1:ncol(W), function(i) quantile(x=(result$gamma)[,i], probs=c(0.025, 0.975)))
colnames(gamma.CI)<-names(gamma0)
gamma.est<-rbind(gamma.CI, gamma.hat)

agl.CI<-cbind(quantile(result$alpha, probs=c(0.025, 0.975)), quantile(result$psi, probs=c(0.025, 0.975)), quantile(result$lambda, probs=c(0.025, 0.975)))
names(agl.CI)<-c("alpha", "psi", "lambda")
agl.est<-rbind(agl.CI, c(alpha.hat, psi.hat, lambda.hat))

# E(y) and var(y)
# E(y)=X%*%beta+alpha*E(V), only for non-censored cases
EV<-besselK(psi.hat, nu = lambda.hat+1, expon.scaled = TRUE)/besselK(psi.hat, nu = lambda.hat, expon.scaled = TRUE)
Ey<-X%*%beta.hat+alpha.hat*EV
names(Ey)<-non.miss.ids
names(mu.hat)<-non.miss.ids
names(sigma.hat)<-non.miss.ids

# var(y)=alpha^2*var(V)+exp(W%*%gamma)*E(V)
varV<-besselK(psi.hat, nu = lambda.hat+2, expon.scaled = TRUE)/besselK(psi.hat, nu = lambda.hat, expon.scaled = TRUE)-
  (besselK(psi.hat, nu = lambda.hat+1, expon.scaled = TRUE)/besselK(psi.hat, nu = lambda.hat, expon.scaled = TRUE))^2
vary<-alpha.hat^2*varV+exp(W[colon.non.miss$status==1, ]%*%gamma.hat)*EV
sqrt(vary)
mean(sqrt(vary))

# compare observed vs predicted life time 
cbind(colon.non.miss$time[colon.non.miss$status==1], Ey[colon.non.miss$status==1, ])
plot(colon.non.miss$time[colon.non.miss$status==1], Ey[colon.non.miss$status==1, ], xlab=expression(y[i]), 
     ylab=expression(E(y[i])), main="Expected Survival time vs. Observed Survival time")

# visual check: p(y_i|x_i, w_i, theta.hat)
yind<-sort(sample(noncensor.id, 9))
par(mfrow=c(2,5))
invisible(
  sapply(yind, function(i){
    yseq<-seq(1,5000,length.out = 2000)
    plot(yseq, dghyp(yseq, mu.hat[names(mu.hat)==i], alpha.hat, sigma.hat[names(sigma.hat)==i], psi.hat, lambda.hat),
         type='l', xlab=parse(text=paste0("y[",i,"]")), ylab="Probability Density", cex=1.5)
    abline(v=c(y[names(y)==i], Ey[names(Ey)==i]),col=c("blue","red"))
  }))
plot(-1, xlim=c(0,1), ylim=c(0,1),xlab="",ylab="",xaxt='n',yaxt='n', axes=FALSE)
par(cex=0.7)
legend("topleft", legend=c(expression(p(y[i]*" | "*hat(theta), x_i, w_i)), "Observed", "Estimated"), fill=c("black","blue","red"))

# check coverage of confidence interval (0, Upper bound)
require(ghyp)
CI<-sapply(noncensor.id, function(i){
  GH.hat<-ghyp(lambda=lambda.hat, chi=psi.hat, psi=psi.hat, mu=as.numeric(mu.hat[row.names(mu.hat)==i]), 
               sigma=as.numeric(sigma.hat[rownames(sigma.hat)==i]), gamma=alpha.hat)
  lower.p<-round(ghyp::pghyp(0, GH.hat),8)
  ghyp::qghyp(1-lower.p, GH.hat)
})