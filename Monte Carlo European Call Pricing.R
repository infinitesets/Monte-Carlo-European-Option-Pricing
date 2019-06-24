library(pracma)

tic() 
#define constants
S_0<-50
T<-1
K<-50
r0<-0.07
sigma<-0.13
n<-10^4 # sample size
d<-12   #time size
delta<-T/d
err <- 0.05

R<-matrix(numeric(n*(d+1)), nrow=n)

R[,1]<-r0
#generate n sample paths of Brownian Motion
for(i in 2:(d+1)){
  R[,i]<-R[,(i-1)]+0.18*(0.086-R[,(i-1)])*delta+0.02*sqrt(delta)*rnorm(n)
}

S<-matrix(numeric(n*(d+1)), nrow=n)
S[,1]<-log(S_0)

for(i in 2:(d+1)){
  S[,i] <- S[,(i-1)]+(R[,(i-1)]+log(1-0.15*sigma^2/2)/0.15)*delta+(sigma*sqrt(rgamma(n, shape=delta/0.15, scale = 0.15))*rnorm(n))
}
S1 <- exp(S)

ar <- delta*apply(R[,-1], 1, sum)

EuroCallPayoff<-pmax(S1[,d+1]-K,0)*exp(-ar*T)

req_sample_size <- ceiling((2.58*1.2*sd(EuroCallPayoff)/err)^2)


R<-matrix(numeric(req_sample_size*(d+1)), nrow=req_sample_size)

R[,1]<-r0

for(i in 2:(d+1)){
  R[,i]<-R[,(i-1)]+0.18*(0.086-R[,(i-1)])*delta+0.02*sqrt(delta)*rnorm(req_sample_size)
}

S<-matrix(numeric(req_sample_size*(d+1)), nrow=req_sample_size)
S[,1]<-log(S_0)

for(i in 2:(d+1)){
  S[,i] <- S[,(i-1)]+(R[,(i-1)]+log(1-((0.15*sigma^2)/2))/0.15)*delta+(sigma*sqrt(rgamma(req_sample_size, shape=(delta/0.15), scale = 0.15))*rnorm(req_sample_size))
}
S1 <- exp(S)

ar <- delta*apply(R[,-1], 1, sum)

EuroCallPayoff<-pmax(S1[,d+1]-K,0)*exp(-ar*T)
EuroCallPrice<-mean(EuroCallPayoff)
error_eurocall<-2.58*sd(EuroCallPayoff)/sqrt(req_sample_size)
error_eurocall
EuroCallPrice
req_sample_size
toc()

tic()
#define constants
S_0<-50
T<-1
K<-50
r0<-0.07
sigma<-0.13
n<-10^3 # sample size
d<-365   #time size
delta<-T/d
err <- 0.05

R<-matrix(numeric(n*(d+1)), nrow=n)

R[,1]<-r0
#generate n sample paths of Brownian Motion
for(i in 2:(d+1)){
  R[,i]<-R[,(i-1)]+0.18*(0.086-R[,(i-1)])*delta+0.02*sqrt(delta)*rnorm(n)
}

S<-matrix(numeric(n*(d+1)), nrow=n)
S[,1]<-log(S_0)

for(i in 2:(d+1)){
  S[,i] <- S[,(i-1)]+(R[,(i-1)]+log(1-0.15*sigma^2/2)/0.15)*delta+(sigma*sqrt(rgamma(n, shape=(delta/0.15), scale = 0.15))*rnorm(n))
}
S1 <- exp(S)

ar <- delta*apply(R[,-1], 1, sum)

EuroCallPayoff<-pmax(S1[,d+1]-K,0)*exp(-ar*T)

req_sample_size <- ceiling((2.58*1.2*sd(EuroCallPayoff)/err)^2)


R<-matrix(numeric(req_sample_size*(d+1)), nrow=req_sample_size)

R[,1]<-r0

for(i in 2:(d+1)){
  R[,i]<-R[,(i-1)]+0.18*(0.086-R[,(i-1)])*delta+0.02*sqrt(delta)*rnorm(req_sample_size)
}

S<-matrix(numeric(req_sample_size*(d+1)), nrow=req_sample_size)
S[,1]<-log(S_0)

for(i in 2:(d+1)){
  S[,i] <- S[,(i-1)]+(R[,(i-1)]+log(1-((0.15*sigma^2)/2))/0.15)*delta+(sigma*sqrt(rgamma(req_sample_size, shape=delta/0.15, scale = 0.15))*rnorm(req_sample_size))
}
S1 <- exp(S)

ar <- delta*apply(R[,-1], 1, sum)

EuroCallPayoff<-pmax(S1[,d+1]-K,0)*exp(-ar*T)
EuroCallPrice<-mean(EuroCallPayoff)
error_eurocall<-2.58*sd(EuroCallPayoff)/sqrt(req_sample_size)
error_eurocall
EuroCallPrice
req_sample_size
toc()


#price Asian geometric mean call option with European call option as control variate
#define constants
source('C:\\Users\\seanc\\Downloads\\Sample_GBM.R')
S_0<-50
T<-1
K<-50
r0<-0.07
sigma<-0.13
n<-10^4 # sample size
d<-52   #time size
delta<-T/d
err <- 0.05



x<-matrix(rnorm(n*d),nrow=n)

#generate n sample paths of Brownian Motion
BM<-sqrt(delta)*t(apply(x,1,cumsum))

#generate n sample paths of stock price
grid<-seq(delta,MT,length.out=d) #time grid
S<-S_0*exp(sweep(sigma*BM,MARGIN=2,(r-sigma^2/2)*grid,'+'))


S<- Sample_GBM(S_0,T,r0,sigma,d,n)

#Use the simple Monte Carlo Method to price European Call Options as control variate
EuroCallPayoff<-pmax(S[,d]-K,0)*exp(-r0*T)
EuroCallPrice<-mean(EuroCallPayoff)
ExactEuroCall<-S_0*pnorm((log(S_0/K)+(r0+sigma^2/2)*T)/(sigma*sqrt(T)))-K*exp(-r0*T)*pnorm((log(S_0/K)+(r0-sigma^2/2)*T)/(sigma*sqrt(T))) #exact price of European call option
EuroCallPrice
error <- abs(ExactEuroCall-EuroCallPrice)
error
#Use monte carlo method with control variate to price Asian Geometric mean call option
GeoMean<-apply(S^(1/d),1,prod) #You must use the sample paths for European option and Asian option to create correlation!!!!!
GeoCallPayoff<-pmax(GeoMean-K,0)*exp(-r0*T)
GeoCallSMPrice<-mean(GeoCallPayoff) #simple monte carlo estimation of geometric mean call option
hat_beta<-cov(GeoCallPayoff,EuroCallPayoff)/var(EuroCallPayoff)
GeoCallMCVPrice<-GeoCallSMPrice+hat_beta*(ExactEuroCall-EuroCallPrice) #estimate with control variate

#compute exact price of Asian geometric mean call option
BarT<-T*(1+1/d)/2
Barsigma<-sqrt(sigma^2*(2+1/d)/3)
Barr<-r0+(Barsigma^2-sigma^2)/2
ExactGeoCall<-(S_0*pnorm((log(S_0/K)+(Barr+Barsigma^2/2)*BarT)/(Barsigma*sqrt(BarT)))-K*exp(-Barr*BarT)*pnorm((log(S_0/K)+(Barr-Barsigma^2/2)*BarT)/(Barsigma*sqrt(BarT))))*exp(Barr*BarT-r0*T)

#true errors
error_sm<-abs(ExactGeoCall-GeoCallSMPrice)
error_mcv<-abs(ExactGeoCall-GeoCallMCVPrice)

cat("With S(0)=",S_0, ", r=",r0,", volatility=",sigma, ", T=", T, ", K=",K, "monitoring frequency=",d,  ", and sample size n=",n, ". The simple MC estimate of Geometric Mean Asian call option price is", GeoCallSMPrice, ", with absolute error as", error_sm, ". The MC estimate with European call option price as control variate of Geometric Mean Asian call option price is", GeoCallMCVPrice, ", with absolute error as", error_mcv, ". " )
