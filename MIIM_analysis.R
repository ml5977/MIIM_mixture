###R packages required
rpackage.list <- c("nlme", "MASS", "mgcv", "splines2", "survival")
lapply(rpackage.list, require, character.only = TRUE)

###simulated data loading
data1 <- read.csv(".../data_generation_continuous.csv") #continuous outcome
data2 <- read.csv(".../data_generation_binary.csv") #binary outcome
data3 <- read.csv(".../data_generation_survival.csv") #survival outcome

###Functions for MIIM
##MIIM for continuous and categorical outcomes
si <- function(theta, y, x, z, p1, p2, p3, weights, family, opt=TRUE, k=5, fx=FALSE){
  #data processing
  beta1 <- theta[c(1:p1)]/sqrt(sum(theta[c(1:p1)]^2))
  beta2 <- theta[c((p1+1):(p1+p2))]/sqrt(sum(theta[c((p1+1):(p1+p2))]^2))
  beta3 <- theta[c((p1+p2+1):(p1+p2+p3))]/sqrt(sum(theta[c((p1+p2+1):(p1+p2+p3))]^2))
  a1 <- as.matrix(x[,c(1:p1)])%*%beta1
  a2 <- as.matrix(x[,c((p1+1):(p1+p2))])%*%beta2
  a3 <- as.matrix(x[,c((p1+p2+1):(p1+p2+p3))])%*%beta3
  
  #GAM formula
  dat <- as.data.frame(cbind(y,a1, a2, a3, z));colnames(dat) <- c("y","a1","a2","a3",colnames(z))
  formula.t <- as.formula(paste('y~ti(a1,fx=',fx,',k=',k,')+ti(a2,fx=',fx,',k=',k,')+ti(a3,fx=',fx,',k=',k,')',
                                '+ti(a1,a2,fx=',fx,',k=',k,')',
                                '+ti(a1,a3,fx=',fx,',k=',k,')',
                                '+ti(a2,a3,fx=',fx,',k=',k,')',
                                "+",paste(colnames(z),collapse = '+')))
  
  #GAM function
  b <- gam(formula.t,data=dat,family=family,weights=weights,method="ML")
  if (opt) return(b$gcv.ubre) else {
    b$beta <- c(beta1, beta2, beta3)
    return(b)
  }
}

##MIIM-HVS for continuous and categorical outcomes (with Lasso penalties)
si.gl <- function(theta, y, x, z, p1, p2, p3, l1, l2, l3, weights, family, opt=TRUE, k=5, fx=FALSE){
  #data processing
  beta1 <- theta[c(1:p1)]/sqrt(sum(theta[c(1:p1)]^2))
  beta2 <- theta[c((p1+1):(p1+p2))]/sqrt(sum(theta[c((p1+1):(p1+p2))]^2))
  beta3 <- theta[c((p1+p2+1):(p1+p2+p3))]/sqrt(sum(theta[c((p1+p2+1):(p1+p2+p3))]^2))
  a1 <- as.matrix(x[,c(1:p1)])%*%beta1
  a2 <- as.matrix(x[,c((p1+1):(p1+p2))])%*%beta2
  a3 <- as.matrix(x[,c((p1+p2+1):(p1+p2+p3))])%*%beta3
  
  #GAM formula
  dat <- as.data.frame(cbind(y,a1, a2, a3, z));colnames(dat) <- c("y","a1","a2","a3",colnames(z))
  formula.t <- as.formula(paste('y~ti(a1,fx=',fx,',k=',k,')+ti(a2,fx=',fx,',k=',k,')+ti(a3,fx=',fx,',k=',k,')',
                                '+ti(a1,a2,fx=',fx,',k=',k,')',
                                '+ti(a1,a3,fx=',fx,',k=',k,')',
                                '+ti(a2,a3,fx=',fx,',k=',k,')',
                                "+",paste(colnames(z),collapse = '+')))
  
  #GAM function
  b <- gam(formula.t,data=dat,family=family,weights=weights,method="ML")
  if (opt) return(b$gcv.ubre+ l1*sum(abs(beta1)) + l2*sum(abs(beta2)) + l3*sum(abs(beta3))) else {
    b$beta <- c(beta1, beta2, beta3)
    return(b)
  }
}

##MIIM for survival outcome
si_surv <- function(theta,t,d,x,z,p1,p2,p3,opt=TRUE,k=5,fx=FALSE) {
  #data processing
  beta1 <- theta[c(1:p1)]/sqrt(sum(theta[c(1:p1)]^2))
  beta2 <- theta[c((p1+1):(p1+p2))]/sqrt(sum(theta[c((p1+1):(p1+p2))]^2))
  beta3 <- theta[c((p1+p2+1):(p1+p2+p3))]/sqrt(sum(theta[c((p1+p2+1):(p1+p2+p3))]^2))
  a1 <- as.matrix(x[,c(1:p1)])%*%beta1
  a2 <- as.matrix(x[,c((p1+1):(p1+p2))])%*%beta2
  a3 <- as.matrix(x[,c((p1+p2+1):(p1+p2+p3))])%*%beta3
  
  #GAM formula
  dat <- as.data.frame(cbind(t, d, a1, a2, a3, z));colnames(dat) <- c("t","d","a1","a2","a3",colnames(z))
  formula.t <- as.formula(paste('t~ti(a1,fx=',fx,',k=',k,')+ti(a2,fx=',fx,',k=',k,')+ti(a3,fx=',fx,',k=',k,')',
                                '+ti(a1,a2,fx=',fx,',k=',k,')',
                                '+ti(a1,a3,fx=',fx,',k=',k,')',
                                '+ti(a2,a3,fx=',fx,',k=',k,')',
                                "+",paste(colnames(z),collapse = '+')))
  
  #GAM function
  b <- gam(formula.t,family=cox.ph, weights=d, method="REML",data=dat)
  if (opt) return(b$gcv.ubre) else {
    b$beta <- c(beta1, beta2, beta3)
    return(b)
  }
}

###Analysis
p1 <- 5; p2 <-2; p3 <- 4
mixture.name <- c("x11","x12","x13","x14","x15","x21","x22","x31","x32","x33","x34")

###MIIM for continuous
#Step 1: initial value (\beta) from conventional linear regression
lm.fit <- lm(y ~ x11 + x12 + x13 + x14 + x15 + x21 + x22 + x31 + x32 + x33 + x34 + z1 + z2, data=data1)
temp.lm.coef <- lm.fit$coefficients[c(2:12)]
temp.lm.coef[c(1:p1)] <- temp.lm.coef[c(1:p1)]/sqrt(sum(temp.lm.coef[c(1:p1)]^2)) #group1
temp.lm.coef[c((p1+1):(p1+p2))] <- temp.lm.coef[c((p1+1):(p1+p2))]/sqrt(sum(temp.lm.coef[c((p1+1):(p1+p2))]^2)) #group2
temp.lm.coef[c((p1+p2+1):(p1+p2+p3))] <- temp.lm.coef[c((p1+p2+1):(p1+p2+p3))]/sqrt(sum(temp.lm.coef[c((p1+p2+1):(p1+p2+p3))]^2)) #group3
th0 <- temp.lm.coef
#Step 2: estimate beta
#get updated beta using no penalization to escape local optimum
f0 <- nlm(si,th0,y=data1$y,x=data1[,mixture.name],z=data1[,c("z1","z2")],p1=p1, p2=p2, p3=p3,
          weights=NULL,family="gaussian",fx=TRUE,k=5);print("f0 finished") 
#get final beta with smoothing parameter selection
f1 <- nlm(si,f0$estimate,y=data1$y,x=data1[,mixture.name],z=data1[,c("z1","z2")],p1=p1, p2=p2, p3=p3,
          weights=NULL,family="gaussian",k=5);print("f1 finished") 
#Step 3: estimate all other parameters
#get link function and alpha
b <- si(f1$estimate,y=data1$y,x=data1[,mixture.name],z=data1[,c("z1","z2")],p1=p1, p2=p2, p3=p3,
        weights=NULL,family="gaussian",opt=FALSE) 
#summary output
summary(b)
plot(b, pages=1)
b$beta


###MIIM-HVS for continuous
f.lambda1 <- f.lambda2 <- f.lambda3 <- 50 #tuning parameters
# get updated beta using no penalization
f0.gl <- nlm(si.gl,th0,y=data1$y,x=data1[,mixture.name],z=data1[,c("z1","z2")],p1=p1, p2=p2, p3=p3,
             l1=f.lambda1,l2=f.lambda2,l3=f.lambda3,weights=NULL,family="gaussian",fx=TRUE,k=5);print("f0.gl finished") 
# get final beta with smoothing parameter selection
f1.gl <- nlm(si.gl,f0.gl$estimate,y=data1$y,x=data1[,mixture.name],z=data1[,c("z1","z2")],p1=p1, p2=p2, p3=p3,
             l1=f.lambda1,l2=f.lambda2,l3=f.lambda3,weights=NULL,family="gaussian",k=5);print("f1.gl finished") 
# get link function and gamma
b.gl <- si.gl(f1.gl$estimate,y=data1$y,x=data1[,mixture.name],z=data1[,c("z1","z2")],p1=p1, p2=p2, p3=p3,
              l1=f.lambda1,l2=f.lambda2,l3=f.lambda3,weights=NULL,family="gaussian",opt=FALSE) 


###MIIM for binary outcome
#initial values from logistic regression
logit.fit <- glm(y ~ x11 + x12 + x13 + x14 + x15 + x21 + x22 + x31 + x32 + x33 + x34 + z1 + z2, family="binomial", data=data2)
temp.lm.coef <- logit.fit$coefficients[c(2:12)]
temp.lm.coef[c(1:p1)] <- temp.lm.coef[c(1:p1)]/sqrt(sum(temp.lm.coef[c(1:p1)]^2)) #group1
temp.lm.coef[c((p1+1):(p1+p2))] <- temp.lm.coef[c((p1+1):(p1+p2))]/sqrt(sum(temp.lm.coef[c((p1+1):(p1+p2))]^2)) #group2
temp.lm.coef[c((p1+p2+1):(p1+p2+p3))] <- temp.lm.coef[c((p1+p2+1):(p1+p2+p3))]/sqrt(sum(temp.lm.coef[c((p1+p2+1):(p1+p2+p3))]^2)) #group3
th0 <- temp.lm.coef
#get updated beta using no penalization to escape local optimum
f0 <- nlm(si,th0,y=data2$y,x=data2[,mixture.name],z=data2[,c("z1","z2")],p1=p1, p2=p2, p3=p3,
          weights=NULL,family="binomial",fx=TRUE,k=5);print("f0 finished") 
#get final beta with smoothing parameter selection
f1 <- nlm(si,f0$estimate,y=data2$y,x=data2[,mixture.name],z=data2[,c("z1","z2")],p1=p1, p2=p2, p3=p3,
          weights=NULL,family="binomial",k=5);print("f1 finished") 
#Step 3: estimate all other parameters
#get link function and alpha
b <- si(f1$estimate,y=data2$y,x=data2[,mixture.name],z=data2[,c("z1","z2")],p1=p1, p2=p2, p3=p3,
        weights=NULL,family="binomial",opt=FALSE) 


###MIIM for survival outcome
#initial values from logistic regression
cox.fit <- coxph(Surv(surv.t, surv.e) ~ x11 + x12 + x13 + x14 + x15 + x21 + x22 + x31 + x32 + x33 + x34 + z1 + z2, data=data3)
temp.lm.coef <- cox.fit$coefficients[c(1:11)]
temp.lm.coef[c(1:p1)] <- temp.lm.coef[c(1:p1)]/sqrt(sum(temp.lm.coef[c(1:p1)]^2)) #group1
temp.lm.coef[c((p1+1):(p1+p2))] <- temp.lm.coef[c((p1+1):(p1+p2))]/sqrt(sum(temp.lm.coef[c((p1+1):(p1+p2))]^2)) #group2
temp.lm.coef[c((p1+p2+1):(p1+p2+p3))] <- temp.lm.coef[c((p1+p2+1):(p1+p2+p3))]/sqrt(sum(temp.lm.coef[c((p1+p2+1):(p1+p2+p3))]^2)) #group3
th0 <- temp.lm.coef
# get updated beta using no penalization
miim.start.time <- Sys.time()
f0 <- nlm(si_surv,th0,t=data3$surv.t,d=data3$surv.e,x=data3[,mixture.name],z=data3[,c("z1","z2")],p1=p1, p2=p2, p3=p3,
          fx=TRUE,k=5);print("f0 finished") 
# get final beta with smoothing parameter selection
f1 <- nlm(si_surv,f0$estimate,t=data3$surv.t,d=data3$surv.e,x=data3[,mixture.name],z=data3[,c("z1","z2")],p1=p1, p2=p2, p3=p3,
          k=5);print("f1 finished") 
# get link function and gamma
b <- si_surv(f1$estimate,t=data3$surv.t,d=data3$surv.e,x=data3[,mixture.name],z=data3[,c("z1","z2")],p1=p1, p2=p2, p3=p3,
             opt=FALSE) 



