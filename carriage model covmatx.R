if(!is.null(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

library(rjags)
# library(coda)
# library(mcmcplots)

set.seed(123)

#load data
ds <- (read.csv("carriage.csv"))[,2:9]

#model formulation
code <-"
model{
for (t in 1:n){
N_carr[t,] ~ dmulti(prev[t,], N_swab[t])  
for (m in 1:k){
numerator[t,m] <- (exp(beta0.st[m] + beta1.st[m]*t1[t] + beta2.st[m]*t2[t]))
}

denominator[t] <- sum(numerator[t,])
for (m in 1:k){
prev[t,m] <- numerator[t,m]/denominator[t]
}}

beta0~dnorm(0, 0.0001) # mean for random intercepts
beta1~dnorm(0, 0.0001) # mean for random slopes
beta2~dnorm(0, 0.0001) # mean for random slopes
sigma0~dunif(0, 100) # SD of intercepts
sigma1~dunif(0, 100) # SD of slopes
sigma2~dunif(0, 100) # SD of slopes
tau.beta1<-1/sigma1^2

rho~dunif(-1, 1) # correlation between intercepts and slopes

#var-covar matrix for the random effects
Sigma[1, 1] <- sigma0^2  
Sigma[2, 2] <- sigma2^2
Sigma[1, 2] <- rho*sigma0*sigma2
Sigma[2, 1] <- rho*sigma0*sigma2

InvSigma[1:2, 1:2] <- inverse(Sigma[,])

for (i in 1:54) {
B.hat[i, 1] <- beta0
B.hat[i, 2] <- beta2
B[i, 1:2]~dmnorm(B.hat[i, ], InvSigma[,]) 
beta0.st[i] <- B[i, 1] 
beta2.st[i] <- B[i, 2]

beta1.st[i]~dnorm( beta1 , tau.beta1)
}


for(j in 1:55){
exp.beta0.st[j] <- exp(beta0.st[j])
exp.beta1.st[j] <- exp(beta1.st[j])
exp.beta2.st[j] <- exp(beta2.st[j])
prev.ratio.fin[j] <- prev[7,j]/prev[1,j]
log.prev.ratio.fin[j] <- log(prev.ratio.fin[j])
for(t in 1:7){
log.prevratio[j,t] <- log(prev.ratio[j,t])
prev.ratio[j,t] <- prev[t,j]/prev[1,j]}
}


beta0.st[k] <- 0
beta1.st[k] <- 0
beta2.st[k] <- 0


}"


# parameters for model input
N_carr <- t(ds[,2:8])     #carriage
N_swab <- rowSums(N_carr) #total swabs each time period
n      <- nrow(N_carr)    #7 time periods
k      <- ncol(N_carr)    #55 serotypes (including negative pneumo)

#categorical time
t0 <- c(1,0,0,0,0,0,0) #pre-vaccine
t1 <- c(0,1,1,1,1,1,1) #early post-vaccine
t2 <- c(0,0,0,1,1,1,1) #late post-vaccine

# run model
dat <- list(n=n, t1=t1, t2=t2, N_carr=N_carr,  N_swab=N_swab, k=k)
mod <- jags.model(textConnection(code), data=dat, n.chains=3)
adapt(object = mod, n.iter = 500000)
update(mod,10000)

post <- coda.samples(mod, thin=5, c('beta0.st', 'beta1.st', 'beta2.st', 
                                         'beta0', 'beta1', 'beta2', 
                                         'prev', 'prev.ratio', 'log.prevratio', 'prev.ratio.fin', 'log.prev.ratio.fin',
                                         'exp.beta0.st', 'exp.beta1.st','exp.beta2.st',
                                         'sigma0','sigma1','sigma2', 'rho'), n.iter=150000)


sum.post <- summary(post)
# View(sum.post$quantiles)

# write.csv(sum.post13$quantiles, "postquants.csv")
# prevsamps <- as.data.frame(jpos13[[1]][,773:1157])
# write.csv(prevsamps, "prevsamps.csv")

#plot(post, col=c("blue", "purple", "green"))

#########################################
#MCMC Diagnostics
#########################################
# test <- sort(effectiveSize(post))
summary(effectiveSize(post)[which(effectiveSize(post)!=0)])   

# gew <- geweke.diag(post)[[1]][[1]] 
hist(geweke.diag(post)[[1]][[1]])

# par(mfrow=c(3,2))
# autocorr.plot(as.mcmc(post[1]))
gelman <- (gelman.diag(post, multivariate = FALSE))

#############################################################################################################
#Inference
#############################################################################################################
par(mfrow=c(1,1))

sts <- (as.character(ds[,1]))

order.b0 <- sts[order(sum.post$quantiles[,3][2:56], decreasing = TRUE)]
order.b1 <- sts[order(sum.post$quantiles[,3][58:112], decreasing = TRUE)]
order.b2 <- sts[order(sum.post$quantiles[,3][113:167], decreasing = TRUE)]
order.prevratio <- sts[order(sum.post$quantiles[,3][1543:1597], decreasing = TRUE)]
order.logprevratio <- sts[order(sum.post$quantiles[,3][333:387], decreasing = TRUE)]



median(sum.post$quantiles[,3][1158:1542])
# [1] 1.303473


#beta0s
par(mfrow=c(1,1))
caterplot(post, parms = "beta0.st", labels = order.b0, style = c("plain"), quantiles = list(inner=c(0.025,0.975)), bty="n", lwd = c(1,1), xlab="Log(Relative Risk Ratio)", cex.labels = 1, col="grey40")
title(expression(paste("Model Estimates of ", beta [0][i])))
mtext(expression(bold('Log(Relative Risk Ratio at time 0)')), side = 1, line = 3)
mtext(expression(bold('Serotype')), side =2, line = 3)
#beta1s
caterplot(post, parms = "beta1.st", labels = order.b1, style = c("plain"), quantiles = list(inner=c(0.025,0.975)), bty="n", lwd = c(1,1), cex.labels = 1, col="grey40")
title(expression(paste("Model Estimates of ", beta [1][i])))
mtext(expression(bold('Log(Relative Risk Ratio per year)')), side = 1, line = 3)
mtext(expression(bold('Serotype')), side =2, line = 3)
#beta2s
caterplot(post, parms = "beta2.st", labels = order.b2, style = c("plain"), quantiles = list(inner=c(0.025,0.975)), bty="n", lwd = c(1,1), cex.labels = 1, col="grey40")
title(expression(paste("Model Estimates of ", beta [2][i])))
mtext(expression(bold('Log(Relative Risk Ratio per year)')), side = 1, line = 3)
mtext(expression(bold('Serotype')), side =2, line = 3)

caterplot(post, parms = "log.prev.ratio.fin", labels = order.logprevratio, style = c("plain"), quantiles = list(inner=c(0.025,0.975)), lwd = c(1,1), bty='n', cex.labels = 1, col="grey40")
title("Log(prevalence ratio) y6/y0")
caterpoints(runif(55, 0, 0), type="l", lty=2, col="gray60", lwd=3)
caterpoints(runif(55, log(1.479), log(1.479)), type="l", lty=4, col="black", lwd=3)
legend(x=.65, y=4.5,legend=c('log(prev. ratio): 0','Constant log(prev. ratio): 0.39'),lwd=2, lty=c(2,4), col=c('gray60','black'),bty = "n")
mtext(expression(bold('Log(prevalence ratio) over study period')), side = 1, line = 3)
mtext(expression(bold('Serotype')), side =2, line = 3)

