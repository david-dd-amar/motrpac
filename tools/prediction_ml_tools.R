library(limma)
library(randomForest)
library(CMA)
library(e1071)
library(tensorflow)
library(keras)
library(RMTL)
# install_keras(method="conda")
# install_tensorflow(version = "1.12",method = "conda")
# set training data

dl_multitask_regression1<-function(x,y,batch_size=1,epochs=20,...){
  
  if(is.null(colnames(x))){colnames(x)=as.character(1:ncol(x))}
  means = apply(x,2,mean)
  sds = apply(x,2,sd)
  cols = sds > 0
  
  x = x[,cols];means=means[cols];sds=sds[cols]
  names(means) = colnames(x)
  names(sds) = colnames(x)
  
  for(j in 1:ncol(x)){x[,j] = (x[,j]-means[j])/sds[j]}
  x_train = x
  input_shape = c(ncol(x),1)
  x_train <- array_reshape(x_train, c(nrow(x_train), ncol(x), 1))
  
  # set model
  model <- keras_model_sequential() %>%
    layer_conv_1d(filters = 1024, kernel_size = 3, activation = 'relu',
                  input_shape = input_shape) %>% 
    layer_dropout(rate = 0.25) %>% 
    layer_flatten() %>% 
    layer_dense(units = 128, activation = 'relu') %>% 
    layer_dropout(rate = 0.5) %>% 
    layer_dense(units = ncol(y), activation = 'softmax')
   
  
   # compile model
    model %>% compile(
      loss = "mean_squared_error", 
      optimizer = "adam"
    )
  
    # fit model
    model %>% fit(
      x = x_train, 
      y = y, 
      epochs = 100,
      batch_size = 1,
      verbose = 1
    )
    
    test_predictions <- model %>% predict(x_train)
    diag(cor(test_predictions,y))
    
    obj = list(model=model,xnames = colnames(x))
    class(obj) = "dl"
}

# Apply func on the columns of x before applying the ml method mlFunc
# If topK is specified, rank the features by their scores (decreasing order)
# and take the topK. If topK is null, use the threshold thr. If both are null,
# do not perform feature selection. 
# ... are parameters for mlFunc.
feature_selection_wrapper<-function(x,y,func,mlFunc,topK=100,thr=NULL,...){
  feature_scores = apply(x,2,func)
  scores_thr = min(feature_scores)
  if(!is.null(topK) && topK < ncol(x)){
    scores_thr = sort(feature_scores,decreasing = T)[topK]
  }
  if(is.null(topK) && !is.null(thr)){
    scores_thr = thr
  }
  selected_features = feature_scores >= scores_thr
  selected_features = colnames(x)[selected_features]
  x = x[,selected_features]
  model = mlFunc(x=x,y=y,...)
  obj = list(model=model,selected_features=selected_features)
  class(obj) = "feature_selection_wrapper"
  return(obj)
}

predict.feature_selection_wrapper<-function(obj,newdata){
  newdata = newdata[,obj$selected_features]
  return(predict(obj$model,newdata=newdata))
}

coeff_of_var<-function(x){
  mu = mean(x)
  if(mu==0){return(0)}
  cv = sd(x)/mu
}

MTL_wrapper<-function(x,y,...){
  xx =list()
  yy = list()
  for(j in 1:ncol(y)){
    xx[[j]] = x
    yy[[j]] = as.matrix(y[,j])
  }
  obj = list(model=MTL(xx,yy,...))
  class(obj) = "MTL_wrapper"
  return(obj)
}

predict.MTL_wrapper<-function(obj,newdata){
  te_x = list()
  for(j in 1:ncol(obj$model$W)){
    te_x[[j]] = as.matrix(newdata)
  }
  raw_preds = predict(obj$model,te_x)
  m = c()
  for(i in 1:length(raw_preds)){
    m = cbind(m,raw_preds[[i]])
  }
  return(m)
}



