###############
#use a latent variable model to find V and S based on observed values of T_obs=V+S+error and S_obs=S+error
##############

###simulate data###
#actual values#
S<-rlnorm(100,5,.1)
V<-rlnorm(100,3,.5)
Tot<-S+V
#observed values#
T_obs<-Tot+rnorm(100,0,5)
S_obs<-S+rnorm(100,0,8)
#get informative priors for S and V based on observed values
#i'm going put a gamma prior on S and V. I want to get
#the parameters for those distributions
#use the method of moments: mean(x)=a/b; var(x)=a/b^2
#if I take log of both sides of each this is a system of linear equations

#do it for S_obs
S_par<-exp(solve(matrix(c(1,1,-1,-2),2,2),matrix(log(c(mean(S_obs),var(S_obs))),2,1)))
#do it for V_est (our ugly inital estimate of V where we ignore negatives)
V_est<-(T_obs-S_obs)[which((T_obs-S_obs)>0)]
V_par<-exp(solve(matrix(c(1,1,-1,-2),2,2),matrix(log(c(mean(V_est),var(V_est))),2,1)))

#let's see how we did
hist(S_obs,freq = F)
curve(dgamma(x,S_par[1],S_par[2]),from=0,to=300,add = T)
hist(V_est,freq = F)
curve(dgamma(x,V_par[1],V_par[2]),from=0,to=100,add=T)
#looks good to me let's use these

#We'll set up a model that mimics the experiment where we are intersted in V and S, but only
#get to see the noisy estimates of S_obs and T_obs
mod<-"model{
for(i in 1:n){
#set prior for unobserved vars

V[i]~dgamma(v_alpha,v_beta)
S[i]~dgamma(s_alpha,s_beta)
#set up observation model as acutal values plus observation error
Sobs[i]~dnorm(S[i],tau.to)
Tobs[i]~dnorm(V[i]+S[i],tau.so)
}
#set some priors
sig.v~dunif(0,1000)
tau.v<-pow(sig.v,-2)
sig.s~dunif(0,1000)
tau.s<-pow(sig.v,-2)
sig.so~dunif(0,1000)
sig.to~dunif(0,1000)
tau.so<-pow(sig.so,-2)
tau.to<-pow(sig.to,-2)
}"

#i'll use jags, which is an external program for mcmc
library(rjags)
#grab the data we want to feed it
jdata<-list('n'=100,'Sobs'=S_obs,"Tobs"=T_obs,"v_alpha"=V_par[1],"v_beta"=V_par[2],"s_alpha"=S_par[1],"s_beta"=S_par[2])
#initialize the model
tot.mod<-jags.model(textConnection(mod),data = jdata,n.chains = 5)
#burn-in
update(tot.mod,1000)
#get samples from it for what we are intested V, S and the obs errors
fit.mod<-coda.samples(tot.mod,c('V','S','sig.to','sig.so'),n.iter = 5000)
#get summary of our 5 chains
f<-summary(fit.mod)
#look at S values
f$statistics[grep("S",rownames(f$statistics)),]
f$quantiles[grep("S",rownames(f$quantiles)),]
#look at S values
f$statistics[grep("V",rownames(f$statistics)),]
f$quantiles[grep("V",rownames(f$quantiles)),]
#notice: no negative values for V and the Estimates of V and S are uncorrelated
cor(f$statistics[grep("S",rownames(f$statistics)),"Mean"],f$statistics[grep("V",rownames(f$statistics)),"Mean"])
#unlike they would be if we used T_obs-S_obs
cor(S_obs,T_obs-S_obs)

#compare out fitted values to the actual values
#check correlation for estimates for V and S only applies to the simulated data. 
# I simulated them to be uncorrelated and wanted to show the this estimation routine 
# preserved that while simple subtraction would induce a spurious correlation. 
#In the real system we would expect V and S to be correlated, i am guessing. 
# when there is more of one there is more of the other.
s.est<-f$statistics[1:100,'Mean']
v.est<-f$statistics[101:200,'Mean']
plot(S,s.est)
summary(lm(s.est~S))
abline(0,1)
plot(V,v.est)
summary(lm(v.est~V))
abline(0,1)
#looks good to me.


