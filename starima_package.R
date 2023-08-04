
#the STARIMA program package aims to do space-time analysis and modelling. 
#It include space-time autocorrelation analysis(STACF), space-time partial autocorrelation analysis(STPACF),and STARIMA modelling (fitting and prediction).

starima_fit <- function (Z=NULL,W=NULL,p=0,d=0,q=0) 
{
  #the program to estimate STARIMA (p,d,q) model based on Hannan-Rissanen calibration 
  #algorithm. 
  #Author: 
  #Dr Jiaqiu Wang, Tao Cheng
  #Spatio-temporal Data Mining Group - University College London
  #London, Britain
  #
  #Version: 2.0  
  #January 2012

  if (d!=0)
  {	
	RAW<-Z;
	Z<-diff(Z,lag=d,differences=1);
  }
  else
  {
  	RAW<-Z;
  }


  Len_total<-nrow(Z);
  N <- ncol(Z);  
  
  # Estimation of STAR model
  
  if (p!=0 && q ==0)
  {
    Z_AR<-lagconvert(Z,W,p);
  
    model<-regress(Z_AR$Z,Z_AR$X);
  
    fit.results<-matrix(as.matrix(model$FIT),length(as.matrix(model$FIT))/N,N);
    fit.res<-matrix(as.matrix(model$RES),length(as.matrix(model$RES))/N,N);
    
  }
  
  # Estimation of STMA model
  
   if (p==0 && q !=0)
   {
     # Stage 1 pre-estimates the epsilon
     
     eps <- epsilon_estimation (Z,W);
     
     # Stage 2, estimating q by using epsilon(t-1)...epsilon(t-q) in stage 1
     
     zx <- preliminary_estimation_ma(Z,W,q,eps$n_len,eps$res);
     
     # Stage 3, calibrating preliminary parameters with iteratives in stage 2
     
     model <- calibrating_ma(W,q,zx$z,zx$x,N);
     
    fit.results<-model$FIT;
    fit.res<-model$RES;
     
   }
  
   # Estimation of STARMA model
  
   if (p!=0 && q !=0)
   {
  
    # Stage 1 pre-estimates the epsilon
  
    eps <- epsilon_estimation (Z,W);
    
          
    # Stage 2, estimating p and q by using epsilon(t-1)...epsilon(t-q) in stage 1

    zx <- preliminary_estimation(Z,W,p,q,eps$n_len,eps$res);

    # Stage 3, calibrating preliminary parameters with iteratives in stage 2

    model <- calibrating (W,p,q,zx$z,zx$x,N);
  
    # Get fitted results
  
    fit.results<-model$FIT;
    fit.res<-model$RES;
        
    }
  
    if (p==0 && q ==0)
    {
      warning("Error: p and q shouldn't be zero at the same time");
      return(0);
    }
    
 
    # Restore differenced series to RAW series
 	s<-recover(RAW,fit.results,d);

    # Evaluate NRMSE
    NRmse<-vector();

  
    for (i in 1:N)
    {
      NRmse[i]<-sqrt(t(s$RES[,i])%*%s$RES[,i]/length(s$RES[,i]))/sd(s$OBS[,i]);
    }
    
  
return(list(BETA=model$BETA,PVAL=model$PVAL,NRMSE=NRmse,OBS=s$OBS,FIT=s$RESULTS,RES=s$RES,P=p,D=d,Q=q,W=W));

}

starima_pre <- function (Z=NULL,model=NULL) 
{
  #the program to estimate STARIMA (p,d,q) model based on Hannan-Rissanen calibration 
  #algorithm. 
  #Author: 
  #Dr Jiaqiu Wang, Tao Cheng
  #Spatio-temporal Data Mining Group - University College London
  #London, Britain
  #
  #Version: 1.3  
  #January 2012

  d<-model$D;

  if (d!=0)
  {	
	RAW<-Z;
	Z<-diff(Z,lag=model$D,differences=1);
  }
  else
  {
  	RAW<-Z;
  }

  
  N <- ncol(Z);
 
  beta<-model$BETA;
  p<-model$P;
  q<-model$Q;
  W<-model$W;

  SpatialLag <- length(names(W));
  
  # STAR model
  if (p!=0 && q==0)
  {
   Z_AR<-lagconvert(Z,W,p);
   LEN<-dim(Z_AR$X);
   Regressor<-as.matrix(Z_AR$X);
   pre<- Regressor%*%beta;
   eps<- as.matrix(Z_AR$Z)-as.matrix(pre);
   
   pre.results<-matrix(pre,length(pre)/N,N);

  # Restore differenced series to RAW series
  s<-recover(RAW,pre.results,d);

  # Evaluate NRMSE
  NRmse<-vector();
  
  for (i in 1:N)
  {
    NRmse[i]<-sqrt(t(s$RES[,i])%*%s$RES[,i]/length(s$RES[,i]))/sd(s$OBS[,i]);
  }

   
   return(list(PRE=s$RESULTS,OBS=s$OBS,NRMSE=NRmse));

  }
  
  # STMA model
  if (p==0 && q!=0 )
  {
    for (i in 1:N)
    {
     if (i==1)
     {
      Z_RAW<-as.matrix(Z[(q+1):nrow(Z),i]);
        
     }
     else
     {
      Z_RAW<-rbind(as.matrix(Z_RAW),as.matrix(Z[(q+1):nrow(Z),i]));
     }
    }
    LEN<-length(Z_RAW);
    epsilon <- matrix(0,LEN,(SpatialLag+1)*q);
  }
  
  # STARIMA model
  if (p !=0 && q!=0)
  {
   Z_AR<-lagconvert(Z[(q+1):nrow(Z),],W,p);
   LEN<-dim(Z_AR$X);
   epsilon <- matrix(0,LEN,(SpatialLag+1)*q);
  }
  
  minimum_rmse<-1000;
  maxIterative<-100;
  
  for (k in 1:maxIterative)
  {
    if (p!=0)        
    {
    Regressor<- cbind(as.matrix(Z_AR$X),epsilon);
    }
    else
    {
    Regressor<- epsilon;
    }
    
    results<- Regressor%*%beta;
    
    if(p!=0)
    {
    eps<- as.matrix(Z_AR$Z)-as.matrix(results);
    }
    else{
    eps<- as.matrix(Z_RAW)-as.matrix(results);
    }
  
    rmse<-sqrt(t(eps)%*%eps/length(eps));
  
    old_epsilon<-epsilon;
  
    eps<-matrix(eps,length(eps)/N,N);
  
    Y_MA<-lagconvert(eps,W,q);
   
    epsilon <- rbind(old_epsilon[1:(q*N),],Y_MA$X);

   if (minimum_rmse>rmse)
   {
     pre.res<-eps;
     pre<-results;
   }
  }
  
  pre.results<-matrix(pre,length(pre)/N,N);
 
  # Restore differenced series to RAW series
  s<-recover(RAW,pre.results,d);

  # Evaluate NRMSE

  NRmse<-vector();
  
  for (i in 1:N)
  {
    NRmse[i]<-sqrt(t(s$RES[,i])%*%s$RES[,i]/length(s$RES[,i]))/sd(s$OBS[,i]);
  }

  res<-s$OBS-s$RESULTS;
   
  return(list(PRE=s$RESULTS,OBS=s$OBS,NRMSE=NRmse,RES=res));

}

stacf <- function(Z=NULL, W=NULL, nLags=20)
{

#the program to estimate space-time autocorrelation function (STACF)

T_len <- nrow(Z)
N <- ncol(Z)

for (i in 1:N)
{
	Z[,i]<-znorm(Z[,i]);
}


# Compute the ST-ACF according to Pfeifer and Deutsch.

stacf <- matrix(0,nLags+1, 1);
stacf[1,1] <- 1;

W<-t(W);

for (s in 1:nLags)
{
	Rlk <- 0;
	Rll <- 0;
	Rkk <- 0;

# Compute the space time autocovariance

	for (tt in 1:(T_len-s))
	{
		Rlk <- Rlk + Z[tt,]%*%t((Z[(tt+s),])%*%W);
	}

	for (tt in 1:T_len)
	{			
		Rll <- Rll + (Z[tt,]%*%W)%*%t((Z[tt,]%*%W));
		Rkk <- Rkk + t(Z[tt,])%*%Z[tt,];
	}

	stacf[s+1,1] <- (Rlk/sqrt(Rll*Rkk))*(T_len/(T_len-s));
}

# Calculate confidence intervals

nSE <- 2
sigmaQ <- sqrt(1/((T_len-nLags)*N));
ubound <- sigmaQ*nSE;
lbound <- sigmaQ*-nSE;

plot(stacf, main="Space-Time Autocorrelation Function",xlab="Time Lags", ylab="ST-Autocorrelation", type="h", col="red", ylim=c(min(stacf), max(stacf)));

#library(fields);

abline(ubound, 0, lty=2, col="blue");
abline(lbound, 0, lty=2, col="blue");
abline(0, 0);

return(list(stacf=stacf, ubound=ubound, lbound=lbound))

}

stpacf <- function (Z=NULL,W=NULL,nLags=20)
{
#the program to compute the space-time partial autocorrelation function (STPACF) 
#   of a univariate, stochastic time series. STPACF is computed by fitting
#   successive space-time autoregressive models of orders 1,2, ... by ordinary least
#   squares, retaining the last coefficient of each regression. When called
#   with no output arguments, STPACF plots the partial correlation sequence with
#   confidence bounds.
#   


T_len <- nrow(Z);

# Initialize partialSTACF

partialSTACF <- matrix(0,nLags+1, 1);
partialSTACF [1,1] <- 1;


W_fit <- list();
W_fit$W<-W;

for (order in 1:nLags)

{

model<-starima_fit(Z,W_fit,order,0,0);

len<-length(model$BETA);

partialSTACF[order+1,1]<-model$BETA[len];

}

# Calculate confidence intervals

nSE <- 2;
P<-0;
sigmaQ <- sqrt(T_len-P-1);
ubound <- nSE/sigmaQ;
lbound <- -nSE/sigmaQ;

plot(partialSTACF, main="Space-Time Partial Autocorrelation Function",xlab="Time Lags", ylab="ST Partial Autocorrelation", type="h", col="red", ylim=c(-0.25, 1));

#library(fields);

abline(ubound, 0, lty=2, col="blue");
abline(lbound, 0, lty=2, col="blue");
abline(0, 0);


return(list(stpacf=partialSTACF, ubound=ubound, lbound=lbound));

}

# Stage 1, implementation

epsilon_estimation <- function (Z,W)
{

  TN <- dim(Z); 
  L<-TN[1];
  N <- TN[2];
  
  if (L<30)
  {a<-trunc(L/10);}
  else
  {a <- 3;}
  
  c<-log(L^a);
  
  aic<-matrix(NA,c,2);
  
  
  for (n in 1:c) 
  {
 
   zx <- lagconvert(Z,W,n);
   
   para<-regress(zx$Z,zx$X);

   epsilon<-matrix(as.matrix(para$RES),length(para$RES)/N,N);
   
   sigma<-var(para$RES);
   
   aic[n,1]<-log(sigma)+(2*n)/length(para$RES);
   aic[n,2]<-n;
   
  }
  
  min_aic<-min(aic[,1]);
  pos<-which.min(aic[,1]);
  n=aic[pos,2];
  
  zx <- lagconvert(Z,W,n);
   
  para<-regress(zx$Z,zx$X);
   
  epsilon<-matrix(as.matrix(para$RES),length(para$RES)/N,N);
  
  fitval<-as.matrix(para$FIT);
  
  fit_val<-matrix(fitval,length(fitval)/N,N);
    
  return(list(n_len=n,fit.results=fit_val,res=epsilon));
  
}

# Stage 2, implementation

preliminary_estimation <- function (Z,W,p,q,n,epsilon)
{
  TN <- dim(Z);
  LEN<-TN[1];
  N<-TN[2];
  
  SpatialLag <- length(names(W));
  
  Y_MA <- lagconvert(epsilon,W,q);
  
  Z_AR <- lagconvert(Z[(n+1):LEN,],W,p);
  
  Z_AR_LEN<-(dim(Z_AR$X))[1];
  Y_MA_LEN<-(dim(Y_MA$X))[1];
  
  if (Z_AR_LEN>Y_MA_LEN)
  {
      start_point<-(Z_AR_LEN-Y_MA_LEN)/N;
      Z_AR <- lagconvert(Z[(n+start_point+1):LEN,],W,p);
  }

  if (Z_AR_LEN<Y_MA_LEN)
  {  
      start_point<-(Y_MA_LEN-Z_AR_LEN)/N;
      EPS_LEN<-(dim(epsilon))[1];
      Y_MA <- lagconvert(epsilon[(start_point+1):EPS_LEN,],W,q);
  }
  
  X<-cbind(as.matrix(Z_AR$X),as.matrix(Y_MA$X));
  
  LEN_X<-(dim(X))[1];
  
  Z_RAW<-lagconvert(Z[(q+1):LEN,],W,p);
  
  LEN_Z<-(dim(as.matrix(Z_RAW$Z)))[1];
  
  eps<-matrix(0,(LEN_Z-LEN_X),(SpatialLag+1)*q);
    
  #for (i in 1:(q*N))
  #{
  #  eps[,i]<-rnorm((LEN_Z-LEN_X));
  #}
  
  #x_ar<-rbind(as.matrix(Z_RAW$X)[(p+1):(LEN_Z-LEN_X),],as.matrix(Z_AR$X));
  x_ar<- as.matrix(Z_RAW$X);
  y_ma<-rbind(as.matrix(eps),as.matrix(Y_MA$X));
  X<-cbind(x_ar,y_ma);
    
  return(list(z=Z_RAW$Z,x=X));
}

preliminary_estimation_ma <- function (Z,W,q,n,epsilon)
{
  TN <- dim(Z);
  LEN<-TN[1];
  N<-TN[2];
  
  SpatialLag <- length(names(W));
  
  Y_MA <- lagconvert(epsilon,W,q);
  
  LEN_MA<-(dim(as.matrix(Y_MA$X)))[1];
  
  for (i in 1:N)
  {
    if (i==1)
    {
      Z_RAW<-as.matrix(Z[(q+1):LEN,i]);
        
    }else
    {
    Z_RAW<-rbind(as.matrix(Z_RAW),as.matrix(Z[(q+1):LEN,i]));
    }
  }
  
  LEN_Z<-length(Z_RAW);
  
  eps<-matrix(0,(LEN_Z-LEN_MA),(SpatialLag+1)*q);
    
  #for (i in 1:(q*N))
  #{
  #  eps[,i]<-rnorm((LEN_Z-LEN_MA));
  #}
  
  y_ma<-rbind(as.matrix(eps),as.matrix(Y_MA$X));
    
  return(list(z=Z_RAW,x=y_ma));
}

# Stage 3, implementation

calibrating <- function (W,p,q,z,x,N)
{
  
  maxIterative <- 100;
  LEN<-(dim(as.matrix(x)))[2];
    
  SpatialLag <- length(W);
  
  for (k in 1:maxIterative)
  {
    old_epsilon<- x[,((SpatialLag+1)*p+1):LEN];
    
    para<-regress(z,x);
    
    eps<-para$RES;
    
    rmse<-sqrt(t(eps)%*%eps/length(eps));
    
    eps<-matrix(as.matrix(eps),length(eps)/N,N);    
    
    if (k!=maxIterative)
    {
      Y_MA<-lagconvert(eps,W,q);
      epsilon<-Y_MA$X;
      epsilon<-rbind(as.matrix(old_epsilon[1:(q*N),]),as.matrix(epsilon));
      x[,(((SpatialLag+1)*p)+1):LEN]<-epsilon;
      
    }
  }
  
  fitval<-as.matrix(para$FIT);
  
  fit_val<-matrix(fitval,length(fitval)/N,N);
  
  return (list(BETA=para$BETA,PVAL=para$PVAL,FIT=fit_val,RES=eps));
  
}

calibrating_ma <- function (W,q,z,x,N)
{
  
  maxIterative <- 100;
  LEN<-(dim(as.matrix(x)))[2];
    
  SpatialLag <- length(W);
  
  for (k in 1:maxIterative)
  {
    old_epsilon<- x;
    
    para<-regress(z,x);
    
    eps<-para$RES;
    
    rmse<-sqrt(t(eps)%*%eps/length(eps));
    
    eps<-matrix(as.matrix(eps),length(eps)/N,N);
    
    if (k!=maxIterative)
    {
      Y_MA<-lagconvert(eps,W,q);
      epsilon<-Y_MA$X;
      epsilon<-rbind(as.matrix(old_epsilon[1:(q*N),]),as.matrix(epsilon));
      x<-epsilon;
      
    }
  }
  
  fitval<-as.matrix(para$FIT);
  
  fit_val<-matrix(fitval,length(fitval)/N,N);
  
  return (list(BETA=para$BETA,PVAL=para$PVAL,FIT=fit_val,RES=eps));
  
}


regress <- function (z,x)
{
  
  lm.model <- lm(z~x-1);
  lm.summary <- summary(lm.model);
 
  beta <- as.vector(lm.model$coef);
  fitted_values<-as.matrix(lm.model$fitted.values);
  pval <- as.vector(lm.summary$coef[,4]);
  eps <- lm.summary$res;
  
  return (list(BETA=beta,PVAL=pval,FIT=fitted_values,RES=eps));
  
}

lagconvert <- function (Z,W,order)
{
  
  rows<-(dim(Z)[1]);
  N<-(dim(Z)[2]);

  for (i in 1:N)
   {

     x_matrix <- lagmatrix(Z[,i],order);
     z_tmp <- x_matrix[(order+1):rows,1];
     x_tmp <- x_matrix[(order+1):rows,2:(order+1)];
 
     if (i==1)
     {
       z<-z_tmp;
       x<-x_tmp;
     }
     else
    {
       z<-c(z,z_tmp);
       x<-rbind(as.matrix(x),as.matrix(x_tmp));
    }
   
   }
   
   if (is.list(W))
    {SpatialLag <- length(W);}
   else
    {return (list(Z=z,X=x));}
     
  
   Z_Lag<-list();
   x_Lag<-list();
   
   for (s in 1:SpatialLag)
   {
    Z_Lag[[s]]<- Z%*%W[[s]];
    x_Lag[[s]]<-matrix(NA,(rows-order)*N,N);
    
    x_Lag_tmp<-matrix(NA,(rows-order),N);
    

    for (i in 1:N)
    {
     x_lag_matrix <- lagmatrix(Z_Lag[[s]][,i],order);
     x_Lag_tmp <- x_lag_matrix[(order+1):rows,2:(order+1)];
     
     if (i==1)
       {x_Lag[[s]]<-x_Lag_tmp;}
     else
       {x_Lag[[s]]<-rbind(as.matrix(x_Lag[[s]]),as.matrix(x_Lag_tmp));}
    }
    x<-cbind(as.matrix(x),as.matrix(x_Lag[[s]]));
   }
     
     
   return (list(Z=z,X=x));
}

recover<-function(Z.raw=NULL,Z.diff=NULL,d=0)
{
	if (d!=0)
	{
		z.len<-nrow(Z.raw);
		fit.len<-nrow(Z.diff);
		obs.len<-fit.len+d;
		obs<-Z.raw[(z.len-obs.len+1):z.len,];
		results<-obs[1:fit.len,] + Z.diff;
		obs.raw<-Z.raw[(z.len-fit.len+1):z.len,];
		res<-obs.raw-results;
	}
	else
	{
		z.len<-nrow(Z.raw);
		fit.len<-nrow(Z.diff);
		results<-Z.diff;
		obs.raw<-Z.raw[(z.len-fit.len+1):z.len,];
		res<-obs.raw-results;

	}
	
	return(list(OBS=obs.raw,RESULTS=results,RES=res));
		
}

lagmatrix <- function (Z,order)
{
  rows<-length(Z);
  cols<-order+1;
  x_matrix<-matrix(NA,rows,cols);
  
  for (i in 0:order)
  {
    x_matrix[(i+1):rows,(i+1)]<-Z[1:(rows-i)];
  }
  
  return(x_matrix);
}

# normalize vector

znorm <- function (Z=NULL)
{

Z<-as.vector(Z);

avg<-mean(Z);
std<-sd(Z);

Zscore<-(Z-avg)/std;

return(t(Zscore));

}

# calculate NRMSE index

NRMSE <- function(res=NULL, obs, pred, sd=TRUE, na.rm=F)
	{
	#Inputs:
	#res = vector of residuals
	#obs = vector of observed values
	#pred = vector of predicted values
	#sd = should standard deviation be used for normalisation? If FALSE, min max is used
	#na.rm = how should missing values (NA) be treated.

	if(is.null(res))
		{
		res <- obs-pred
		}
	if(sd==T)
		{
		NRMSE <- (sqrt(mean((res)^2, na.rm=na.rm)))/(sd(obs))
		}
	else
		{
		NRMSE <- (sqrt(mean((res)^2, na.rm=na.rm)))/(max(obs, na.rm=na.rm)-min(obs, na.rm=na.rm))
		}
	return(NRMSE=NRMSE)
	}