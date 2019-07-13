
regress_out<-function(x,covs,...){
  form = x~.
  d = data.frame(x=x,covs)
  m = lm(form,data=d,...)
  return(m$residuals)
}

# Simple pairwise
get_correlation_matrix<-function(x,y,z,...){
  if(!is.null(z)){
    x = t(apply(x,2,regress_out,covs=z))
    y = t(apply(y,2,regress_out,covs=z))
  }
  return(cor(x,y,...))
}

get_explained_variance_using_PCA<-function(x,y,ncomp=10,...){
  xpca = prcomp(x,center = T,scale. = F)
  ypca = prcomp(y,center = T,scale. = F)
  xvars = xpca$sdev^2; yvars = ypca$sdev^2
  xtotvar = sum(xvars); ytotvar = sum(yvars)
  xprop = xvars/xtotvar; yprop = yvars/ytotvar
  
  # For QCs the y explained by x below should fit the r2s
  # r2s = c()
  # for(j in 1:ncomp){
  #   if(is.null(z)){
  #     df = data.frame(y = ypca$x[,j],xpca$x[,1:ncomp])
  #   }
  #   else{
  #     df = data.frame(y = ypca$x[,j],xpca$x[,1:ncomp],z=z)
  #   }
  #   lm_res = summary(lm(y~.,data=df))
  #   r2s[j] = lm_res$r.squared
  # }
  
  pc_cross_cor = cor(xpca$x[,1:ncomp],ypca$x[,1:ncomp])
  pc_cross_r2 = pc_cross_cor^2
  y_explained_by_x = colSums(pc_cross_r2)
  explained_prop = sum(y_explained_by_x*yprop[1:ncomp])
  explained_prop = explained_prop / sum(yprop[1:ncomp])
  return(explained_prop)
}

apply_function_on_pairs<-function(l,f,...){
  n = length(l)
  m = matrix(NA,nrow=n,ncol=n,dimnames = list(names(l),names(l)))
  for(i in 1:n){
    for(j in 1:n){
      m[i,j] = f(l[[i]],l[[j]],...)
      m[j,i] = m[i,j]
      if(i==j){break}
    }
  }
  return(m)
}



