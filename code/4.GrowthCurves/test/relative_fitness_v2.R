### Don's script to attempt to understand and partially re-create the results
### sent by Jay and Canan

####load the data####
#grab the clone names
clns<-c("anc","M26","M23","M54","M41","M21","M19","M17","M79","M13","M4")

##set http targets for data pull
#targ<-NULL
#for(i in 1:length(clns))targ[i]<-paste0("https://raw.githubusercontent.com/LennonLab/SporeMut/master/code/4.GrowthCurves/output/",clns[i],"_new.fit.parms.txt")

##grab data
#dat<-list()
#for(i in 1:length(clns))dat[[i]]<-read.table(targ[i],sep = ",",header = T)
#names(dat)<-clns
input_dat<-read.csv("comp_data.csv")
head(input_dat)
dat_order<-c(1:6,unlist(lapply(paste0("^",clns[-1],"_"),function(x)grep(x,input_dat$Curve))))
ord_list<-rep(clns,each=6)

#put this is form this code used previously
dat<-split(input_dat[dat_order,-1],f=factor(ord_list,levels = clns))



#simple back-of-the-envelope calculations first to see find relative fitness of M79
#find weighted mean and variance for both ancestor and m79
w_anc<-(1/dat[['anc']]$umax.se^2)/sum((1/dat[['anc']]$umax.se^2))
anc_mean<-sum(w_anc*dat[['anc']]$umax)
anc_var<-sum(w_anc*(dat[['anc']]$umax-anc_mean)^2)/(1-sum(w_anc^2))

w_m79<-(1/dat[['M79']]$umax.se^2)/sum((1/dat[['M79']]$umax.se^2))
m79_mean<-sum(w_m79*dat[['M79']]$umax)
m79_var<-sum(w_m79*(dat[['M79']]$umax-m79_mean)^2)/(1-sum(w_m79^2))

#what is relative fitness
rf<-m79_mean/anc_mean;rf
#what is the uncertainty
rf_se<-m79_mean/anc_mean*sqrt((sqrt(m79_var)/m79_mean)^2+(sqrt(anc_var)/anc_mean)^2)
#what are the (naive) 95% CI of this
#lower
rf-1.96*rf_se
#upper
rf+1.96*rf_se
# Hmm relative fitness of M79 quite different than the estimates in the "umax-clones.pdf" graph
# BUT the effect is the same, the lower confidence interval > 1 (so far so good)

#set up data for analyses
var.mat<-data.frame(sapply(dat,"[",7))^2
#make weights based on inverse variance within clones
wts_mat<-sweep(1/var.mat,MARGIN = 2,STATS =colSums(1/var.mat) ,FUN = "/")

#find weighted means
apply(as.data.frame(sapply(dat,"[",6))*as.matrix(wts_mat),2,sum)
####
#plot raw: these are not the real mean and se because each point has
#different amounts of uncertainty that we cannot just throw away.
naive_mean<-apply(as.data.frame(sapply(dat,"[",6)),2,mean)
naive_se<-apply(as.data.frame(sapply(dat,"[",6)),2,sd)/sqrt(6)
plt_naive<-t(rbind(naive_mean-1.96*naive_se,naive_mean,naive_mean+1.96*naive_se))
plot(plt_naive[11:2,2],1:10,ylab="clone",xlab='naive estimate of umax',
     xlim=c(0.01,.16),yaxt='n',pch=21,bg=16,bty='n',main=bquote(mu*'max'))
for(i in 11:2)points(dat[[i]]$umax,rep(12-i,6),pch=21,col="skyblue",bg="skyblue")
points(plt_naive[11:2,2],1:10,pch=21,bg=16)
arrows(plt_naive[11:2,2],1:10,plt_naive[11:2,1],1:10,length = 0)
arrows(plt_naive[11:2,2],1:10,plt_naive[11:2,3],1:10,length = 0)
abline(v=1,lty=2)
axis(side = 2,at=1:10,labels = clns[11:2],las=2)
####

#create a long-form data.frame for later use
umax_df<-data.frame(umax=stack(as.data.frame(sapply(dat,"[",6)))$values,
                            lvl=as.factor(sapply(strsplit(as.character(stack(as.data.frame(sapply(dat,"[",6)))$ind),"\\."),"[",1)))
                    

library(rjags)
# now let's do it bayes-style. from what I gather from the code, the sampling is simply done using the 6
# replicates for each clone and the ancestor. So, I will do that I will assume that the umax for each rep (for each clone)
# is sampling a population mean. in addition, will use the umax_se to define weights

#make a jags model that estimates a umax for the replicates of clones and the ancestor,
# using the weights the estimates. we will sample a log-normal dist for umax to make
# sure the estimates are all positive during sampling and assume each clone is sampling a dist
# with unique variance
mod<-"model{
     for(i in 1:6){
       for(j in 1:10){
       #here is the model for each clone
       M_mu[i,j]~dlnorm(mu_m[j],tau_m[j]*wts_m[i,j])
       #M_mu[i,j]~dlnorm(mu_m[j],tau_m[j])
       }
       #here is the model of the ancestor
      An_mu[i]~dlnorm(mu_a,tau_a*wts_a[i])
      #An_mu[i]~dlnorm(mu_a,tau_a)
     }
     #priors
     for(i in 1:10){
      mu_m[i]~dnorm(0,.001)
      tau_m[i]~dgamma(.01,.01)
     }
     mu_a~dnorm(0,.001)
     tau_a~dgamma(.01,.01)
     
     #calc in here just for fun. I will do the rest 'by hand'
     rf79<-exp(mu_m[8])/exp(mu_a)
     diff79<-exp(mu_m[8])-exp(mu_a)
}"     
#Fit it
#gather data to feed jags
jdat<-list("M_mu"=as.data.frame(sapply(dat,"[",6))[,-1],"An_mu"=dat[['anc']]$umax,"wts_m"=wts_mat[,2:11],"wts_a"=wts_mat[,1])
#build the model object
jmod<-jags.model(textConnection(mod),n.chains = 4,data = jdat,n.adapt = 1000)
#update it
update(jmod,1000)
#sample it
ans<-coda.samples(jmod,c("mu_m","mu_a"),10000)
#look at plot to make sure they converge
x11() #windows() if on a Windows machine
#par(ask=T)
plot(ans)
#the samples look great so sets peek at the summaries of the posterior dists
summ<-summary(ans)
#check out the median of the posterior estimates
# we exponentiate them because we used a log-normal
exp(summ$quantiles)

#lets compare these medians with the simple weighted means
exp(summ$quantiles[,"50%"])
apply(as.data.frame(sapply(dat,"[",6))*as.matrix(wts_mat),2,sum)
#nice. all makes sense and is working

##Let's calcuate the relative fitness across all chains

#combine all of the chains 
vals<-rbind(as.matrix(ans[[1]]),as.matrix(ans[[2]]),as.matrix(ans[[3]]),as.matrix(ans[[4]]))
#find the quartiles of the relative fitness
apply(vals,2,function(x)median(exp(x)))
#as another check on our code, these should duplicate what we saw in the summary step above
exp(summ$quantiles[,"50%"])
#bingo!

#first lets look at m79 and compare it to our back of envelope calcs
m79_rel<-quantile(exp(vals[,9])/exp(vals[,1]),c(.025,.5,.975))
#notice these are more variable than our back-of-envelope 
rbind(m79_rel,c(rf-1.96*rf_se,rf,rf+1.96*rf_se)) 
#close, but more variation here than back-of-envelope estimate (which was only approximate)

#calc relative fitness for rest of clones and plot them in a style
#similar to umax-clones.pdf
plt_qnts<-matrix(0,10,3)
rownames(plt_qnts)<-clns[-1]
colnames(plt_qnts)<-c("2.5%","50%","97.5%")
for(i in 1:10)plt_qnts[i,]<-quantile(exp(vals[,i+1])/exp(vals[,1]),c(.025,.5,.975))
#kill the old window
dev.off()
#make a new one
x11() #or windows()
#make the relative fitness vs close plots
plot(plt_qnts[10:1,2],1:10,ylab="clone",xlab='Relative fitness',
     xlim=c(0,2),yaxt='n',pch=21,bg=16,bty='n',main=bquote(mu*'max'))
arrows(plt_qnts[10:1,2],1:10,plt_qnts[10:1,1],1:10,length = 0)
arrows(plt_qnts[10:1,2],1:10,plt_qnts[10:1,3],1:10,length = 0)
abline(v=1,lty=2)
axis(side = 2,at=1:10,labels = clns[11:2],las=2)

# whoa...some very different results here M26, M23, M19, M17, M13, M4 are all quite
# different than the umax-clones-pdf plot results with respect to the relationship to
# the vertical line at 1 

#Do the same analysis but in the way that Lennon Lab did it
#Next I will do this with the long-form data

#set up data for one-hot coding
one_hot<-cbind(anc=rep(1,dim(umax_df)[1]),sapply(clns[-1],function(x)as.numeric(umax_df$lvl==x)))

#this model is using a Gaussian with a single variance and fitting the 
# ancestor as the intercept and each clone as a different beta coef. The same as 
# you would get from the R formula umax~clone where the ancestor is coded as the first factor in "clone"
# here we use a different version of the weights because we don't want them relative by clone 
mod_1<-"model{
     for(i in 1:n){
       y[i]~dnorm(mu[i],tau/wts[i])
       mu[i]<-inprod(X[i,],b)
     }
     #priors
     for(i in 1:11){
      b[i]~dnorm(0,.001)
     }
     tau~dgamma(.01,.01)
     
}"     
#Fit it

jdat<-list("y"=umax_df$umax,"X"=one_hot,"n"=dim(umax_df)[1],"wts"=1000*stack(var.mat)$value)
jmod_1<-jags.model(textConnection(mod_1),n.chains = 4,data = jdat,n.adapt = 1000)
update(jmod_1,1000)
#sample the betas. we will use these to reconstruct the previous graph and the umax-clones.pdf graph
umax_dif<-coda.samples(jmod_1,c("b"),10000)
summary(umax_dif)
# let's peek at the tau to make sure it is converging 
plot(coda.samples(jmod_1,c("tau"),10000))

#next, we will use the posterior chains to create relative fitness from coefs
rel_fit<-t(sapply(paste0("b[",2:11,"]"),function(x)quantile((as.matrix(umax_dif[[1]])[,"b[1]"]+as.matrix(umax_dif[[1]])[,x])/(as.matrix(umax_dif[[1]])[,"b[1]"]),c(.025,.5,.975))))
rownames(rel_fit)<-clns[-1]

#plot these. they should be very similar to the last graph
x11()

plot(rel_fit[10:1,2],1:10,ylab="clone",xlab='Fitness difference',xlim=c(0,2),yaxt='n',pch=21,bg=16,bty='n',main=bquote(mu*'max'~'Version 2'))
arrows(rel_fit[10:1,2],1:10,rel_fit[10:1,1],1:10,length = 0)
arrows(rel_fit[10:1,2],1:10,rel_fit[10:1,3],1:10,length = 0)
abline(v=1,lty=2)
axis(side = 2,at=1:10,labels = clns[11:2],las=2)

#compared to the previous,method this one slightly diff variation around estimates
# the estimates themselves are very close
cbind(rel_fit[,2],plt_qnts[,2]) #math is still working!

#NOW about the x-axis weirdness on umax-clones.pdf. 
#look at fitness DIFFERENCES not relative fitness (this is what is plotted in the figure I was sent
# but a 1 was added to the x-axis value)
umax_post<-summary(umax_dif)
umax_ans<-umax_post$quantiles[2:11,c("2.5%","50%","97.5%")]
rownames(umax_ans)<-clns[-1]

plot(umax_ans[10:1,2],1:10,ylab="clone",xlab='Fitness difference',xlim=c(-.1,.1),yaxt='n',pch=21,bg=16,bty='n')
arrows(umax_ans[10:1,2],1:10,umax_ans[10:1,1],1:10,length = 0)
arrows(umax_ans[10:1,2],1:10,umax_ans[10:1,3],1:10,length = 0)
abline(v=0,lty=2)
axis(side = 2,at=1:10,labels = clns[11:2],las=2)

#notice the x-axis. next I will remake this but change the zero before the decimals on the x-axis
# to ones...note this is not the same as adding one

plot(umax_ans[10:1,2],1:10,ylab="clone",xlab='Fitness difference',xlim=c(-.1,.1),yaxt='n',xaxt='n',pch=21,bg=16,bty='n')
arrows(umax_ans[10:1,2],1:10,umax_ans[10:1,1],1:10,length = 0)
arrows(umax_ans[10:1,2],1:10,umax_ans[10:1,3],1:10,length = 0)
abline(v=0,lty=2)
axis(side = 2,at=1:10,labels = clns[11:2],las=2)
axis(side=1, at=c(-.1,-.05,0,.05,.1),labels= c(NA,-1.05,1,1.05,NA))




