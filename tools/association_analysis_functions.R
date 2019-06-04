#' Compute linear association between a continuous variable y and another
#' variable x, which is transformed into a factor if it is discrete.
#' Z is a dataframe to be conditioned upon: compute y's residuals with it.
#' x,y,z, are all named vectors/matrices
#' @param x - vector (numeric or discrete)
#' @param y - numeric vector
#' @param z - optional - a data frame of variables to adjust for
linear_association_analysis<-function(x,y,z=NULL){
  if(!is.null(z)){
    df1 = data.frame(y=y,z)
    lm1 = lm(y~.,data=df1)
    y = lm1$residuals
  }
  inds = intersect(names(x),names(y))
  x = x[inds];y=y[inds]
  inds = !is.na(x) & !is.na(y)
  x = x[inds];y=y[inds]
  xcopy = as.numeric(as.character(x))
  if(sum(is.na(x))>length(x)/2){
    x = as.factor(x)
  }
  df = data.frame(y=y,x=x)
  lm_model1 = lm(y~x,data=df)
  lm_model0 = lm(y~1,data=df)
  lm_res = summary(lm_model1)
  lm_an = anova(lm_model1,lm_model0)
  pval = as.matrix(lm_an)[2,6]
  r2 = lm_res$r.squared
  rho = sqrt(r2)
  betas = lm_res$coefficients[-1,1]
  names(betas) = rownames(lm_res$coefficients)[-1]
  if(length(betas)==1 && betas[1]<0){rho=-rho}
  return(c(pval=pval,r2=r2,rho=rho,betas))
}

#' Go over the column pairs in m and run the function func on each one. 
#' Store the value from the field f in the function's output in the result
pairwise_eval<-function(m,func,f,...){
  n = ncol(m)
  res = matrix(NA,nrow=n,ncol=n)
  colnames(res) = colnames(m)
  rownames(res) = colnames(m)
  for(i in 1:n){
    for(j in 1:i){
      curr_output = func(m[,i],m[,j],...)
      res[i,j] = curr_output[f]
      res[j,i] = res[i,j]
    }
  }
  return(res)
}

# # test
# n=100
# x = rbinom(n,1,0.3)
# y = rnorm(n,sd = 0.2) + x
# z = matrix(rnorm(n*10),ncol=10,nrow=n)
# rownames(z) = 1:n
# names(x) = 1:n
# names(y) = 1:n
# cor(x,y)
# cor.test(x,y,method="spearman")
# lmtest_res = linear_association_analysis(x,y,z)
# abs(lmtest_res["rho"] - cor(x,y))
# corrs = pairwise_eval(z,linear_association_analysis,"rho")
# corrs2 = cor(z)
# hist(corrs - corrs2)
# sum(corrs - corrs2 > 0.00001)
