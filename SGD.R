## load the data
train = read.csv('forestfiretrain.csv')
test = read.csv('forestfiretest.csv')

getRMSE <- function(pred, actual) {
  ## TODO: given a vector or predicted values and actual values, calculate RMSE
  return (sqrt(mean((pred-actual)^2,na.rm=TRUE)))
}

## Function to calculate R^2
getrsq<-function(predictedy, y)
{
  return(1-(sum((y-predictedy)^2))/(sum((y-mean(y))^2)))
}

addIntercept <- function(mat) {
  ## add intercept to the matrix
  allones= rep(1, nrow(mat))
  return(cbind(Intercept=allones, mat))
}

predictSamples <- function(beta, mat) {
    ## TODO: compute the predicted value using matrix multiplication
    ## Note that for a single row of mat, pred = sum_i (beta_i * feature_i)
    return (mat %*% beta)
}

MAX_EPOCH = 100

sgd <- function(learn.rate, train, test, epoch=MAX_EPOCH) {
  ## convert the train and test to matrix format
  train.mat = as.matrix(train) 
  test.mat = as.matrix(test)

  ## TODO: get the number of rows in the train matrix
  N = nrow(train.mat)
  ## TODO: get the number of dimensions (columns) in the train matrix
  d = ncol(train.mat)

  ## standardize the columns of both matrices
  for (i in 1:(d-1)){
    ## TODO: standardize the train and test matrices
    colAvg = mean(train.mat[ ,i])
    colSD = sd(train.mat[, i])
    train.mat[, i] = (train.mat[ ,i] - colAvg)/colSD
    test.mat[, i] = (test.mat[ ,i] - colAvg)/colSD
  }

  ## add a feature to represent the intercept
  tmat <- addIntercept(train.mat[, -d])
  testmat <- addIntercept(test.mat[, -d])

  ## initialize all the coefficients to be 0.5
  beta = rep(0.5,d)
  j = 1
  mse.df <- NULL
  # predict training residuals
  pred_train = predictSamples(beta, tmat)
  pred_test = predictSamples(beta, testmat)
  tMse = getRMSE(pred_train, train$area)
  testMSE = getRMSE(pred_test, test$area)
  tR2 = getrsq(pred_train, train$area)
  testR2 = getrsq(pred_test, test$area)
  mse.df <- rbind(mse.df, data.frame(epoch=j, trainMSE=tMse, testMSE=testMSE, trainR2=tR2, testR2=testR2))

  while(j < MAX_EPOCH){  
    j=j+1;
    # for each row in the training data
    for (n in seq(1:N)){
      ##TODO: update beta according to slide #6 in APA-reg2 
      beta = beta + learn.rate * tmat[n,] *((train$area[n] - tmat[n,] %*% beta) );
    }
    pred_train = predictSamples(beta, tmat)
    pred_test = predictSamples(beta, testmat)
    # tmp_test <- data.frame(pred=pred_test, actual=test$area, type="test")
    # tmp_train <- data.frame(pred=pred_train, actual=train$area, type="train")
    # tmp <- rbind(tmp_train, tmp_test)
    # ggplot(tmp, aes(x=pred, y=actual, color=type)) + theme_bw() + geom_point()

    tMse = getRMSE(pred_train, train$area)
    testMSE = getRMSE(pred_test, test$area)
    tR2 = getrsq(pred_train, train$area)
    testR2 = getrsq(pred_test, test$area)
    mse.df <- rbind(mse.df, data.frame(epoch=j, trainMSE=tMse, testMSE=testMSE, trainR2=tR2, testR2=testR2))
  } 
  return(mse.df)
}

## learning rate 1
l1.df <- sgd(0.00025, train, test)
## learning rate 2
l2.df <- sgd(0.00575, train, test)
## learning rate 3
l3.df <- sgd(0.0065, train, test)

## plot the data
library(reshape2)
library(ggplot2)
ml.df <- melt(l1.df, id=c("epoch"))
ggplot(subset(ml.df, variable %in% c("trainMSE", "testMSE")), aes(x=epoch, y=value, color=variable)) + geom_point() + geom_line() + 
  theme_bw() + scale_color_brewer(palette="Dark2") + xlab("Epoch") + ylab("RMSE") + 
  facet_grid(variable ~ ., scales="free") + theme(legend.position="none")
ggsave(file="Q4-L1.pdf", width=8, height=6)

ml.df <- melt(l2.df, id=c("epoch"))
ggplot(subset(ml.df, variable %in% c("trainMSE", "testMSE")), aes(x=epoch, y=value, color=variable)) + geom_point() + geom_line() + 
  theme_bw() + scale_color_brewer(palette="Dark2") + xlab("Epoch") + ylab("RMSE") + 
  facet_grid(variable ~ ., scales="free") + theme(legend.position="none")
ggsave(file="Q4-L2.pdf", width=8, height=6)

ml.df <- melt(l3.df, id=c("epoch"))
ggplot(subset(ml.df, variable %in% c("trainMSE", "testMSE")), aes(x=epoch, y=value, color=variable)) + geom_point() + geom_line() + 
  theme_bw() + scale_color_brewer(palette="Dark2") + xlab("Epoch") + scale_y_log10("RMSE") + 
  facet_grid(variable ~ ., scales="free") + theme(legend.position="none")
ggsave(file="Q4-L3.pdf", width=8, height=6)


## TODO: fit an MLR model to it
fit <- lm(area~FFMC+DMC+DC+ISI+temp+RH+wind+rain, data=train)
## calculate MSE for LM model
lm.mse = getRMSE(predict(fit, test), test$area)
lm.train.mse = getRMSE(predict(fit,train), train$area)

library(xtable)
results = data.frame(Model="LM", Train=lm.train.mse, Test=lm.mse)
results <- rbind(results, data.frame(Model="SGD L1", 
        Train=l1.df$train[l1.df$epoch==25], Test=l1.df$test[l1.df$epoch==25]))
results <- rbind(results, data.frame(Model="SGD L2", 
        Train=l2.df$train[l2.df$epoch==10], Test=l2.df$test[l2.df$epoch==10]))
results <- rbind(results, data.frame(Model="SGD L3", 
        Train=l3.df$train[l3.df$epoch==1], Test=l3.df$test[l3.df$epoch==1]))
print(xtable(results, digits=3), include.rownames=FALSE)
