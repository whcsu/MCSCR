#breslow h0(t) S0(t)
crcure=function(Y,X,emmax = 50,eps = 1e-04,nboot = 10){
  #baseline function
  smsurv <-
    function(Time,Status,X,beta,w,k){    
      death_point <- sort(unique(subset(Time, Status==k)))
      coxexp <- exp((beta)%*%t(X))  
      lambda <- numeric()
      event <- numeric()
      for(i in 1: length(death_point)){
        event[i] <- sum(is.na(as.numeric(Time==death_point[i]&Status==k)))
        temp <- sum(as.numeric(Time>=death_point[i])*w*drop(coxexp))
        temp1 <- event[i]
        lambda[i] <- temp1/temp
      }
      HHazard <- numeric()
      for(i in 1:length(Time)){
        HHazard[i] <- sum(as.numeric(Time[i]>=death_point)*lambda)
        if(Time[i]>max(death_point))HHazard[i] <- Inf
        if(Time[i]<min(death_point))HHazard[i] <- 0
      }
      survival <- exp(-HHazard)
      list(survival=survival,lambda=lambda)
    }
  library(cmprsk)
  library(survival)
  #data prepare
  n = nrow(Y)
  colnames(Y)=c("Time","Status")
  Time = Y[,1]
  Status = Y[,2]
  DT=data.frame(Y,X)  
  DT1=DT[,-1] 
  DT1$Status=ifelse(DT$Status==1,1,0)   
  Z=as.matrix(cbind(rep(1,n),X))
  colnames(Z) = c("(Intercept)",colnames(X))
  bnm = colnames(Z)
  nb = ncol(Z)
  betanm = colnames(X)
  nbeta = ncol(X)
  
  #initial prameter
  b=glm(Status~.,data=DT1,family=quasibinomial(link="logit"))$coef
  beta1=crr(DT$Time, DT$Status, X, failcode=1, cencode=0)$coef
  beta2=crr(DT$Time, DT$Status, X, failcode=2, cencode=0)$coef
  
  #em
  em <-
    function(DT,X,Z,b,beta1,beta2,emmax,eps)
    { 
      I1=ifelse(DT$Status==1,1,0)
      I2=ifelse(DT$Status==2,1,0)
      deta=ifelse(DT$Status==0,0,1)  
      w=I1 
      n <- length(w)
      
      l2 <- function(beta1,X,lamda1,s01,I1,w) {
        ebX1=exp(drop(beta1 %*% t(X)))
        h1=lamda1*ebX1;s1=s01^ebX1
        log_likelihood1=I1*log(w*h1)+w*log(s1)
        return(-sum(log_likelihood1))
      }
      l3 <- function(beta2,X,I2,lamda2,s02,w) {
        ebX2=exp(drop(beta2 %*% t(X)))
        h2=lamda2*ebX2;s2=s02^ebX2
        log_likelihood2=I2*log(h2)+log(s2)
        return(-sum(log_likelihood2))
      }
      #em
      convergence=1000;i=1
      while (convergence > eps & i < emmax){  
        pai=matrix(exp((b)%*%t(Z))/(1+exp((b)%*%t(Z))),ncol=1)
        s01 <- smsurv(DT[,1],DT[,2],X,beta1,w,1)$survival 
        s02 <- smsurv(DT[,1],DT[,2],X,beta2,w,2)$survival
        lamda1=smsurv(DT[,1],DT[,2],X,beta1,w,1)$lambda
        lamda2=smsurv(DT[,1],DT[,2],X,beta2,w,2)$lambda
        s1=drop(s01^(exp((beta1)%*%t(X))))
        s2=drop(s02^(exp((beta2)%*%t(X))))
        
        # E step 
        w <- I1+((1-I1)*pai*s1)/(1-pai+pai*s1)
        DTb=data.frame(w,X)
        
        # M step
        logit=glm(w~.,data=DTb,family=quasibinomial(link="logit"))
        update_b <- logit$coef
        opbeta1=optim(par=c(rep(0.5,nbeta)),fn=l2,w=w,X=X,I1=I1,lamda1=lamda1,s01=s01)
        opbeta2=optim(par=c(rep(0.5,nbeta)),fn=l3,w=w,X=X,I2=I2,lamda2=lamda2,s02=s02)
        update_beta1=opbeta1$par;update_beta2=opbeta2$par
        update_s01 <-smsurv(DT[,1],DT[,2],X,beta1,w,1)$survival
        update_s02 <-smsurv(DT[,1],DT[,2],X,beta2,w,2)$survival
        
        convergence<-sum(c(update_b-b,beta1-update_beta1,update_beta2-beta2)^2)
        b <- update_b
        beta1 <- update_beta1 
        beta2 <- update_beta2
        s01<-update_s01
        s02<-update_s02
        i <- i+1
      }
      em <- list(logit=logit,fg1=fg1,s01=s01,s02=s02,fg2=fg2,b=b,beta1=beta1,beta2=beta2,s1=s1,s2=s2,tau=convergence)
    }
  
  #estimate result
  emfit <- em(DT,X,Z,b,beta1,beta2,50,1e-4)
  logit=emfit$logit
  b=emfit$b;beta1=emfit$beta1;beta2=emfit$beta2
  s1=emfit$s1;s2=emfit$s2
  s01=emfit$s01;s02=emfit$s02
  
  #bootstrap
  b_boot<-matrix(rep(0,nboot*nb), nrow=nboot)
  beta1_boot<-matrix(rep(0,nboot*(nbeta)), nrow=nboot)
  beta2_boot<-matrix(rep(0,nboot*(nbeta)), nrow=nboot)
  iter <- matrix(rep(0,nboot),ncol=1)
  data1<-subset(DT,Status==1);data0<-subset(DT,Status==0)
  data2=subset(DT,Status==2)
  n1<-nrow(data1);n0<-nrow(data0);n2=nrow(data2)
  i<-1
  while (i<=nboot){
    id1<-sample(1:n1,n1,replace=TRUE);id0<-sample(1:n0,n0,replace=TRUE)
    id2<-sample(1:n2,n2,replace=TRUE)
    bootdata<-rbind(data1[id1,],data0[id0,],data2[id2,])
    bootX <- bootdata[,betanm]
    bootZ=as.matrix(cbind(rep(1,n),bootX))
    bootfit <- em(bootdata,bootX,bootZ,b,beta1,beta2,100,1e-7)
    b_boot[i,] <- bootfit$b
    beta1_boot[i,] <- bootfit$beta1
    beta2_boot[i,] <- bootfit$beta2
    if (bootfit$tau<eps) i<-i+1}
  b_var <- apply(b_boot, 2, var)
  beta1_var <- apply(beta1_boot, 2, var)
  beta2_var <- apply(beta2_boot, 2, var)
  b_sd <- sqrt(b_var)
  beta1_sd <- sqrt(beta1_var)
  beta2_sd <- sqrt(beta2_var)
  fit<-list()
  fit$logit<- logit
  fit$fg1=fg1;fit$fg=fg2
  fit$b <- b
  fit$beta1 <- beta1;fit$beta2 <- beta2
  fit$b_var <- b_var
  fit$b_sd <- b_sd
  fit$b_zvalue <- fit$b/b_sd
  fit$b_pvalue <- (1-pnorm(abs(fit$b_zvalue)))*2
  fit$beta1_var <- beta1_var;fit$beta2_var <- beta2_var
  fit$beta1_sd <- beta1_sd;fit$beta2_sd <- beta2_sd
  fit$beta1_zvalue <- fit$beta1/beta1_sd;fit$beta2_zvalue <- fit$beta2/beta2_sd
  fit$beta1_pvalue <- (1-pnorm(abs(fit$beta1_zvalue)))*2	;fit$beta2_pvalue <- (1-pnorm(abs(fit$beta2_zvalue)))*2	
  fit$s1 <- s1;fit$s2 <- s2
  fit$Time <- Time
  fit$s01=s01;fit$s02=s02
  return(fit)
}