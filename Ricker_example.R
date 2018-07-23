##########################
###Ricker example
##########################

library(ggplot2)
library(reshape2)
library(pomp)

pompExample(ricker)

#plot pomp object to get time series
plot(ricker)

#simulate to pomp object, try to simulate from that process, make new pomp object
x = simulate(ricker)

plot(x) 
#N and e are latent variables (e is normally distributed random variable with mean 0 and SD sigma)

y = as.data.frame(ricker)
y.sim = as.data.frame(x)

#can simulate and get dataframe directly
y.sim = simulate(ricker,as.data.frame=T,nsim=10)

#structure command, gives view of structure of object
str(y.sim)

#data viewed of single realisation of model, can plot against other simulations
x = simulate(ricker,nsim=9, as.data.frame=T,include.data=TRUE)
ggplot(data=x,aes(x=time,y=y))+geom_line()+facet_wrap(sim~.,ncol=2)

#coefficients
coef(ricker)
coef(ricker,'phi')
coef(ricker) = c(phi=20,c=1,r=44,sigma=0.3,N.0=10,e.0=0)
coef(ricker)

#help
methods?pomp

#examples
install.packages('pompExamples',repos='https://kingaa.github.io')
library(pompExamples)


#pfilter pomp object
pf = pfilter(ricker, Np=1000)
class(pf)

logLik(pf)

pf = pfilter(pf)
logLik(pf)

pf = pfilter(pf,Np=100)
logLik(pf)

