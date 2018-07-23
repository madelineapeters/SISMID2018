#Read in data
bsflu_data <- read.table("https://kingaa.github.io/sbied/stochsim/bsflu_data.txt")

#Create name vectors
statenames = c("S","I","R1")
paramnames = c("Beta","mu_I","mu_R1","rho")

#Evaluation of y_n given x_n and theta
dmeas = Csnippet("
  lik = dpois(B,rho*R1+1e-6,give_log);
")

#Draw of y_n given x_n and theta; add small value 1e-6 to avoid sampling from Poisson distribution 0
#Need small value added to have neglible effect on where parameters end up
rmeas = Csnippet("
  B = rpois(rho*R1+1e-6);
")

#Create object with equations
rproc = Csnippet("
  double N = 763;
  double t1 = rbinom(S,1-exp(-Beta*I/N*dt));
  double t2 = rbinom(I,1-exp(-mu_I*dt));
  double t3 = rbinom(R1,1-exp(-mu_R1*dt));
  S  -= t1;
  I  += t1 - t2;
  R1 += t2 - t3;
")

init = Csnippet("
   S = 762;
   I = 1;
   R1 = 0;
")

#These snippets implement parameter transformations that we'll want
#For dealing with constraints on parameters
#Exponentially transform parameters
fromEst = Csnippet("
    TBeta = exp(Beta);
    Tmu_I = exp(mu_I);
    Trho = expit(rho);
")

#Log transform parameters
toEst = Csnippet("
  TBeta = log(Beta);
  Tmu_I = log(mu_I);
  Trho = logit(rho);
")

#Build pomp object
bsflu = pomp(
  data=subset(bsflu_data,select=-C),
  times="day",t0=0,
  rmeasure=rmeas,dmeasure=dmeas,
  rprocess=euler.sim(rproc,delta.t=1/12),
  initializer=init,
  fromEstimationScale=fromEst,toEstimationScale=toEst,
  statenames=statenames,
  paramnames=paramnames
)

plot(bsflu,main="")

#Testing the codes (check rprocess and rmeasure codes)

#Set up test parameters
params = c(Beta=2,mu_I=1,rho=0.9,mu_R1=1/3,mu_R2=1/2)

#Run and plot some simulations
y = simulate(bsflu,params=params,nsim=10,as.data.frame=TRUE)

#Plot
library(ggplot2) 
theme_set(theme_bw()) 
library(reshape2)
ggplot(data=y,mapping=aes(x=time,y=B,group=sim))+
  geom_line()

#Check basic particle filter (pf depends on rprocess and dmeasure)
pf = pfilter(bsflu,params=params,Np=1000)
plot(pf)

#Setting up the estimation problem

#Treat mu_R1 and mu_R2 as known, and fix these at empirical means
fixed_params = with(bsflu_data,c(mu_R1=1/(sum(B)/512),mu_R2=1/(sum(C)/512)))

#Parallelize process
library(foreach)
library(doParallel)
registerDoParallel()

#Running particle filter
library(doRNG)
registerDoRNG(625904618)
pf = bake(file="pf.rds",{
  foreach(i=1:10,.packages='pomp',
          .export=c("bsflu","fixed_params")
  ) %dopar% {
    pfilter(bsflu,params=c(Beta=2,mu_I=1,rho=0.9,fixed_params),Np=10000)
  }
})
L_pf = logmeanexp(sapply(pf,logLik),se=TRUE)

#Building a picture of likelihood surface
results = as.data.frame(as.list(c(coef(pf[[1]]),loglik=L_pf[1],loglik=L_pf[2])))
write.csv(results,file="bsflu_params.csv",row.names=FALSE)

#A local search of the likelihood surface
registerDoRNG(482947940)
mifs_local = 
  foreach(i=1:20,
          .packages='pomp', #load package
          .combine=c, #combine using c operator
          .export=c("bsflu","fixed_params")
  ) %dopar% {
    mif2( bsflu,
          start=c(Beta=2,mu_I=1,rho=0.9,fixed_params),
          Np=2000, #number particles, can use less with iterative filtering because of improved fit with parameter perturbations
          Nmif=50, #numer of interative filters
          cooling.type="geometric",
          cooling.fraction.50=0.5, #after 50 iterations, perturbations decreased to half original size
          transform=TRUE,
          rw.sd=rw.sd(Beta=0.02,mu_I=0.02,rho=0.02) #initial sd of perturbations
    )
  }

#plot results
ggplot(data=melt(conv.rec(mifs_local)), #convergence record - plot of what happened as we did those iterations
       aes(x=iteration,y=value,group=L1,color=factor(L1)))+
  geom_line()+
  guides(color=FALSE)+
  facet_wrap(~variable,scales="free_y")+
  theme_bw()
#loglik not likelihood of model we're interested in, but likelihood of model with perturbations

#Final filtering iteration generates approximation to likelihood at the results point estimate, but not a good enough reliable inference:
  #1) Parameter perturbations are applied in the last filtering iteration, so likelihood shown here is not identical to the model of interest
  #2)mif2 is usually carried out with fewer particles than needed for a good likelihood evaluations: errors in mif2 average out over many iterations of filtering
#Evaluate the likelihood, together with SE, using replicated particles filters at each point estimate:
registerDoRNG(900242057)
results_local = bake(file="lik_local.rds",{
  foreach(mf=mifs_local,.packages='pomp',.combine=rbind) %dopar%
  {
    evals <- replicate(10, logLik(pfilter(mf,Np=20000))) #taking mf runs, running pf on it; need a lot more particles because don't have added variability
    ll <- logmeanexp(evals,se=TRUE)
    c(coef(mf),loglik=ll[1],loglik=ll[2])
  }
})
results_local = as.data.frame(results_local)

pairs(~loglik+Beta+mu_I+rho,data=results_local,pch=16)

#Add points to our database
results <- rbind(results,results_local[names(results)])
write.csv(results,file="bsflu_params.csv",row.names=FALSE)

#Global search of likelihood surface using randomised starting values
#Local search algorith: starts somewhere and looks nearby

params_box = rbind( #set based on final range of parameter estimates from original pf
  Beta=c(1,5),
  mu_I=c(0.5,3),
  rho=c(0.5,1)
)

registerDoRNG(1270401374)
guesses = as.data.frame(apply(params_box,1,function(x)runif(300,x[1],x[2])))
mf1 = mifs_local[[1]]
results_global = bake(file="box_search_global.rds",{
  foreach(guess=iter(guesses,"row"), #300 guesses from params_box
          .packages='pomp',
          .combine=rbind,
          .options.multicore=list(set.seed=TRUE),
          .export=c("mf1","fixed_params")
  ) %dopar% {
    mf = mif2(mf1,start=c(unlist(guess),fixed_params))
    mf = mif2(mf,Nmif=100) #100 iterations
    ll = replicate(10,logLik(pfilter(mf,Np=100000))) #compute loglik
    ll = logmeanexp(ll,se=TRUE) #mean and SE on loglik
    c(coef(mf),loglik=ll[1],loglik=ll[2])
  }
})
#cool more quickly than theory suggests, but run multiple times (heat and cool, heat and cool)

#store output to database
results_global <- as.data.frame(results_global)
results <- rbind(results,results_global[names(results)])
write.csv(results,file="bsflu_params.csv",row.names=FALSE)
#SE less than 1, can think of more meaningful quantity

######Exercise##########################################################
#Fitting SEIR^3 model
########################################################################

#Create name vectors
statenames = c("S","E","I","R1","R2")
paramnames = c("Beta","mu_E","mu_I","mu_R1","mu_R2","rho")

#Evaluation of y_n given x_n and theta
dmeas = Csnippet("
                 lik = dpois(B,rho*R1+1e-6,give_log);
                 ")

#Draw of y_n given x_n and theta; add small value 1e-6 to avoid sampling from Poisson distribution 0
#Need small value added to have neglible effect on where parameters end up
rmeas = Csnippet("
                 B = rpois(rho*R1+1e-6);
                 ")

#Create object with equations
rproc = Csnippet("
                    double N = 763;
                    double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
                    double dN_EI = rbinom(E,1-exp(-mu_E*dt));
                    double dN_IR = rbinom(I,1-exp(-mu_I*dt));
                    double dN_R1R2 = rbinom(R1,1-exp(-mu_R1*dt));
                    double dN_R2R3 = rbinom(R2,1-exp(-mu_R2*dt));
                    S -= dN_SE;
                    E += dN_SE - dN_EI;
                    I += dN_EI - dN_IR;
                    R1 += dN_IR - dN_R1R2;
                    R2 += dN_R1R2 - dN_R2R3;
                    ")

init = Csnippet("
                    S = 762;
                    E = 0;
                    I = 1;
                    R1 = 0;
                    R2 = 0;
                    ")

#These snippets implement parameter transformations that we'll want
#For dealing with constraints on parameters
#Exponentially transform parameters
fromEst = Csnippet("
                   TBeta = exp(Beta);
                   Tmu_I = exp(mu_I);
                   Tmu_E = exp(mu_E);
                   Tmu_R1 = exp(mu_R1);
                   Tmu_R2 = exp(mu_R2);
                   Trho = expit(rho);
                   ")

#Log transform parameters
toEst = Csnippet("
                 TBeta = log(Beta);
                 Tmu_I = log(mu_I);
                 Tmu_E = log(mu_E);
                 Tmu_R1 = log(mu_R1);
                 Tmu_R2 = log(mu_R2);
                 Trho = logit(rho);
                 ")

#Build pomp object
bsflu = pomp(
  data=subset(bsflu_data,select=-C),
  times="day",t0=0,
  rmeasure=rmeas,dmeasure=dmeas,
  rprocess=euler.sim(rproc,delta.t=1/12),
  initializer=init,
  fromEstimationScale=fromEst,toEstimationScale=toEst,
  statenames=statenames,
  paramnames=paramnames
)

plot(bsflu,main="")

#Testing the codes (check rprocess and rmeasure codes)

#Set up test parameters
params = c(Beta=2,mu_I=1,mu_E = 1,rho=0.9,mu_R1=1/3,mu_R2=1/2)

#Run and plot some simulations
y = simulate(bsflu,params=params,nsim=10,as.data.frame=TRUE)

#Plot
library(ggplot2) 
theme_set(theme_bw()) 
library(reshape2)
ggplot(data=y,mapping=aes(x=time,y=B,group=sim))+
  geom_line()

#Check basic particle filter (pf depends on rprocess and dmeasure)
pf = pfilter(bsflu,params=params,Np=1000)
plot(pf)

#Setting up the estimation problem

#Treat mu_R1 and mu_R2 as known, and fix these at empirical means
fixed_params = with(bsflu_data,c(mu_R1=1/(sum(B)/512),mu_R2=1/(sum(C)/512)))

#Parallelize process
library(foreach)
library(doParallel)
registerDoParallel()

#Running particle filter
library(doRNG)
registerDoRNG(625904618)
pf =
  foreach(i=1:10,.packages='pomp',
          .export=c("bsflu","fixed_params")
  ) %dopar% {
    pfilter(bsflu,params=c(Beta=2,mu_I=1,mu_E=1,rho=0.9,fixed_params),Np=10000)
  }

L_pf = logmeanexp(sapply(pf,logLik),se=TRUE)

#Building a picture of likelihood surface
results = as.data.frame(as.list(c(coef(pf[[1]]),loglik=L_pf[1],loglik=L_pf[2])))
write.csv(results,file="bsflu_params.csv",row.names=FALSE)

#A local search of the likelihood surface
registerDoRNG(482947940)
mifs_local = 
  foreach(i=1:20,
          .packages='pomp', #load package
          .combine=c, #combine using c operator
          .export=c("bsflu","fixed_params")
  ) %dopar% {
    mif2( bsflu,
          start=c(Beta=2,mu_I=1,mu_E=1,rho=0.9,fixed_params),
          Np=5000, #number particles, can use less with iterative filtering because of improved fit with parameter perturbations
          Nmif=50, #numer of interative filters
          cooling.type="geometric",
          cooling.fraction.50=0.5, #after 50 iterations, perturbations decreased to half original size
          transform=TRUE,
          rw.sd=rw.sd(Beta=0.02,mu_E=0.02,mu_I=0.02,rho=0.02) #initial sd of perturbations
    )
  }

#plot results
ggplot(data=melt(conv.rec(mifs_local)), #convergence record - plot of what happened as we did those iterations
       aes(x=iteration,y=value,group=L1,color=factor(L1)))+
  geom_line()+
  guides(color=FALSE)+
  facet_wrap(~variable,scales="free_y")+
  theme_bw()
#loglik not likelihood of model we're interested in, but likelihood of model with perturbations

#Final filtering iteration generates approximation to likelihood at the results point estimate, but not a good enough reliable inference:
#1) Parameter perturbations are applied in the last filtering iteration, so likelihood shown here is not identical to the model of interest
#2)mif2 is usually carried out with fewer particles than needed for a good likelihood evaluations: errors in mif2 average out over many iterations of filtering
#Evaluate the likelihood, together with SE, using replicated particles filters at each point estimate:
registerDoRNG(900242057)
results_local = 
  foreach(mf=mifs_local,.packages='pomp',.combine=rbind) %dopar%
  {
    evals <- replicate(10, logLik(pfilter(mf,Np=20000))) #taking mf runs, running pf on it; need a lot more particles because don't have added variability
    ll <- logmeanexp(evals,se=TRUE)
    c(coef(mf),loglik=ll[1],loglik=ll[2])
  }

results_local = as.data.frame(results_local)

pairs(~loglik+Beta+mu_E+mu_I+rho,data=results_local,pch=16)

#Add points to our database
results <- rbind(results,results_local[names(results)])
write.csv(results,file="bsflu_params.csv",row.names=FALSE)

#Global search of likelihood surface using randomised starting values
#Local search algorith: starts somewhere and looks nearby

params_box = rbind( #set based on final range of parameter estimates from original pf
  Beta=c(1,5),
  mu_E=c(1,7),
  mu_I=c(0.5,3),
  rho=c(0.5,1)
)

registerDoRNG(1270401374)
guesses = as.data.frame(apply(params_box,1,function(x)runif(10,x[1],x[2])))
mf1 = mifs_local[[1]]
results_global = 
  foreach(guess=iter(guesses,"row"), #10 guesses from params_box
          .packages='pomp',
          .combine=rbind,
          .options.multicore=list(set.seed=TRUE),
          .export=c("mf1","fixed_params")
  ) %dopar% {
    mf = mif2(mf1,start=c(unlist(guess),fixed_params))
    mf = mif2(mf,Nmif=10) #50 iterations
    ll = replicate(10,logLik(pfilter(mf,Np=500))) #compute loglik
    ll = logmeanexp(ll,se=TRUE) #mean and SE on loglik
    c(coef(mf),loglik=ll[1],loglik=ll[2])
  }
#cool more quickly than theory suggests, but run multiple times (heat and cool, heat and cool)

#store output to database
results_global <- as.data.frame(results_global)
results <- rbind(results,results_global[names(results)])
write.csv(results,file="bsflu_params.csv",row.names=FALSE)
#SE less than 1, can think of more meaningful quantity

library(plyr)
all <- ldply(list(guess=guesses, result=subset(results, loglik > max(loglik)-50)), .id="type") 
pairs(~loglik+Beta+mu_I+mu_E+rho, data=all, col=ifelse(all$type=="guess", grey(0.5), "red"), pch=16)
