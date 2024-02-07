library(survival)
library(riskRegression)
library(cmprsk)
library(pec)
library(randomForestSRC)
library(ggplot2)
library(flexsurv)

#estimate result
mcdata1 = gen2(base='weibull',N=1000,px1=0.5,b1=c(0.7,-0.8,1),beta1=c(-0.8,1),beta2=c(1,-0.8),
               w1=c(0.8,1.2),w2=c(1,1.2),low=0,up=5)
dt=mcdata1[,c(1,2,3,4)]
Y=dt[,c(1,2)]
X=dt[,c(3,4)]

#estimate result
fit=crcure(Y,X,emmax = 50,eps = 1e-04,nboot = 100)
estb=fit$b
estbeta1=fit$beta1
estbeta2=fit$beta2
ests01=fit$s01

x=c(0,0.032);z=c(1,x) #predict x

#my pre
estpai=estb%*%z
esteb1=exp(estbeta1%*%x)
f1_pre=estpai*(1-ests01^esteb1)

#fine-gray
crrfit1 <- crr(dt$Time, dt$Status, X, failcode=1, cencode=0)
predict_crr <- predict.crr(crrfit1, data.frame('x1'=0,'x2'=0.032))
fg1=predict_crr[,2]

#cox
cox_model <- coxph(Surv(Time, Status==1) ~ x1+x2, data = dt)
coxpre <- survfit(cox_model, newdata = data.frame('x1'=0,'x2'=0.032))$cumhaz
coxpre1=coxpre[dt$Status==1]
coxpre <- 1 - exp(-coxpre1)

#crrsf
rfsrc_model <- rfsrc(Surv(Time, Status) ~ ., data = dt, ntree = 100, cause = 1)
rsfpre=predict(rfsrc_model,cause=1,newdata=data.frame('x1'=0,'x2'=0.032),times=t_values)
crrsf=rf[dt$Status==1]

#cs
cs <- CSC(Hist(Time, Status) ~ x1 + x2, data = dt,cause=1)
cspre=predict(cs,cause=1,newdata=data.frame('x1'=0,'x2'=0.032),times=t_values)
cs1=cspre[1]$absRisk

#mixture cure
library(flexsurv)
cure_model <- flexsurvreg(Surv(Time, Status==1) ~ x1 + x2, data = dt, dist = "weibull")
mc <- predict(cure_model, newdata = data.frame('x1'=0,'x2'=0.032),type = "cumhaz",times=t_values)
mc=mc[[1]][[1]][,2]
mcpre=1 - exp(-mc);mcpre=mcpre[1:a,1]