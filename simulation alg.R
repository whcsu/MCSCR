gen2=function(base,N,px1,b1,beta1,beta2,w1,w2,low,up) {
  x1=rbinom(N,1,px1)
  x2=rnorm(N,0,1)
  x=cbind (x1,x2)
  bX1=x%*%beta1;bX2=x%*%beta2
  # proportion of the uncured
  logit1 = cbind (1,x1,x2) %*% b1
  pi1 = exp(logit1)/(1+exp(logit1)) 
  y = rbinom (N,1,prob=pi1)
  u=runif (N,0,1);C = runif(N,low,up)
  if(base=='exp'){
    f2=function(t,i,w2,bX2,u){
      ({exp(-w2*t)}^(exp(bX2)[i]))-u[i]
    }
    f1=function(t,i,w1,bX1,w2,bX2,u){
      ({exp(-w1*t)}^(exp(bX1)[i]))*({exp(-w2*t)}^(exp(bX2)[i]))-u[i]
    }
    t1=t2=numeric(N)
    for(j in 1:N){
      t1[j] <- uniroot(f1,interval=c(0,200),i=j,w1=w1,bX1=bX1,w2=w2,bX2=bX2,u=u,tol=0.0001)$root
      t2[j] <- uniroot(f2,interval=c(0,200),i=j,w2=w2,bX2=bX2,u=u,tol=0.0001)$root
    }
    h1=exp(bX1)*w1
    h2=exp(bX2)*w2
  }
  if(base=='weibull'){
    f2=function(t,i,w2,bX2,u){
      ({exp(-(w2[1]*t)^(w2[2]))}^(exp(bX2)[i]))-u[i]
    }
    f1=function(t,i,w1,bX1,w2,bX2,u){
      ({exp(-(w1[1]*t)^(w1[2]))}^(exp(bX1)[i]))*({exp(-(w2[1]*t)^(w2[2]))}^(exp(bX2)[i]))-u[i]
    }
    t1=t2=numeric(N)
    for(j in 1:N){
      t1[j] <- uniroot(f1,interval=c(0,200),i=j,w1=w1,bX1=bX1,w2=w2,bX2=bX2,u=u,tol=0.0001)$root
      t2[j] <- uniroot(f2,interval=c(0,200),i=j,w2=w2,bX2=bX2,u=u,tol=0.0001)$root
    }
    h1=exp(bX1)*w1[1]*w1[2]*(w1[1]*t1)^(w1[2]-1)
    h2=exp(bX2)*w2[1]*w2[2]*(w2[1]*t1)^(w2[2]-1)
  }
  if(base=='lognorm'){
    f2=function(t,i,w2,bX2,u){
      ({1-pnorm((log(t)-w2[1])/w2[2])}^(exp(bX2)[i]))-u[i]
    }
    f1=function(t,i,w1,bX1,w2,bX2,u){
      ({1-pnorm((log(t)-w1[1])/w1[2])}^(exp(bX1)[i]))*({1-pnorm((log(t)-w2[1])/w2[2])}^(exp(bX2)[i]))-u[i]
    } 
    t1=t2=numeric(N)
    for(j in 1:N){
      t1[j] <- uniroot(f1,interval=c(0,200),i=j,w1=w1,bX1=bX1,w2=w2,bX2=bX2,u=u,tol=0.0001)$root
      t2[j] <- uniroot(f2,interval=c(0,200),i=j,w2=w2,bX2=bX2,u=u,tol=0.0001)$root
    }
    h1=exp(bX1)*exp(-(log(t1)-w1[1])^2/(2*w1[2]^2))/(w1[2]*t1*sqrt(2*pi))
    h2=exp(bX2)*exp(-(log(t1)-w2[1])^2/(2*w2[2]^2))/(w2[2]*t1*sqrt(2*pi))
  }
  p1=h1/(h1+h2)
  y1 = rbinom (N,1,prob=p1)
  Status1=ifelse(y1==1,1,2)
  Status=ifelse(y==0,2,Status1)
  Time=ifelse(y==0,t2,t1)
  deta=as.numeric(Time<=C)
  Time=ifelse(deta==0,C,Time)
  Status=ifelse(deta==0,0,Status)
  s=u
  mcdata = data.frame (Time,Status,x,y,s)
}
mcdata1 = gen2(base='weibull',N=1000,px1=0.5,b1=c(0.7,-0.8,1),beta1=c(-0.8,1),beta2=c(1,-0.8),
               w1=c(0.8,1.2),w2=c(1,1.2),low=0,up=5)
dt=mcdata1[,c(1,2,3,4,5)]