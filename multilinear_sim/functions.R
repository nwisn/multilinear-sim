#--------------------------------------------------------------------------------
# modeling functions

# fit OLS linear model
run_lm <- function(title = "OLS",
                   weights = NULL,
                   predictors = NULL,
                   Xsub.training,
                   Xsub.testing,
                   Y.training){
  if(is.null(predictors)) predictors <- colnames(Xsub.training[,-1])
  df <- cbind(data.frame(Y=Y.training), as.data.frame(Xsub.training[,-1][,predictors]))
  if(is.null(weights)) weights <- rep(1, nrow(df))
  model <- lm(Y ~ ., data = df, weights = weights)
  Yhat.train <- predict(model, newdata=as.data.frame(Xsub.training[,-1]))
  Yhat.test <- predict.lm(model, newdata=as.data.frame(Xsub.testing[,-1]))
  ret <- list(model=model, Yhat.train=Yhat.train, Yhat.test=Yhat.test, title=title)
  return(ret)
}

# fit glmnet model
run_glmnet <- function(alpha = 0,
                       title = "ridge",
                       measure = "mse",
                       lambda = "lambda.1se",
                       weights = NULL,
                       dfmax = NULL,
                       pmax = NULL,
                       Xsub.training,
                       Xsub.testing,
                       Y.training
){
  require(glmnet)
  if(is.null(weights)) weights <- rep(1, nrow(Xsub.training))
  if(is.null(dfmax)) dfmax <- ncol(Xsub.training) + 1
  if(is.null(pmax)) pmax <- min(dfmax *2+20, ncol(Xsub.training))
  model <- cv.glmnet(Xsub.training[,-1], Y.training, alpha = alpha, type.measure=measure, weights=weights, dfmax=dfmax, pmax=pmax)
  Yhat.train <- predict(model, newx=Xsub.training[,-1], lambda=model[[lambda]])
  Yhat.test <- predict(model, newx=Xsub.testing[,-1], lambda=model[[lambda]])
  ret <- list(model=model, Yhat.train=as.numeric(Yhat.train), Yhat.test=as.numeric(Yhat.test), title=title, lambda=lambda)
  return(ret)
}

# get model results
get_model_results <- function(results,
                              alpha = 1,
                              Xsub.training,
                              Xsub.testing,
                              Y.training,
                              Y.testing){
  require(ggsci)
  require(gridExtra)
  df.y.training <- data.frame(Y=Y.training, Yhat=results$Yhat.train, error=results$Yhat.train-Y.training, subset=rep("training set", length(Y.training)), model=rep(results$title, length(Y.training)))
  df.y.testing <- data.frame(Y=Y.testing, Yhat=results$Yhat.test, error=results$Yhat.test-Y.testing, subset=rep("testing set", length(Y.testing)),  model=rep(results$title, length(Y.testing)))

  cortest.training <- cor.test(df.y.training$Y, df.y.training$Yhat, method = "p", exact = F)
  cortest.testing <- cor.test(df.y.testing$Y, df.y.testing$Yhat, method = "p", exact = F)
  r.training <- cortest.training$estimate
  r.lower.training <- cortest.training$conf.int[1]
  r.upper.training <- cortest.training$conf.int[2]
  r.testing <- cortest.testing$estimate
  r.lower.testing <- cortest.testing$conf.int[1]
  r.upper.testing <- cortest.testing$conf.int[2]
  r.training.p <- cortest.training$p.value
  r.testing.p <- cortest.testing$p.value

  rmse.training <- sqrt(mean(df.y.training$error^2))
  rmse.testing <- sqrt(mean(df.y.testing$error^2))
  rmse.training.ci <- quantile(replicate(3000, sqrt(mean(sample(df.y.training$error^2, replace = T)))), c(0.025, 0.975))
  rmse.testing.ci <- quantile(replicate(3000, sqrt(mean(sample(df.y.testing$error^2, replace = T)))), c(0.025, 0.975))
  rmse.training.lower <- rmse.training.ci[1]
  rmse.training.upper <- rmse.training.ci[2]
  rmse.testing.lower <- rmse.testing.ci[1]
  rmse.testing.upper <- rmse.testing.ci[2]

  ttest.training <- t.test(df.y.training$error)
  ttest.testing <- t.test(df.y.testing$error)
  bias.training.lower <- ttest.training$conf.int[1]
  bias.training.upper <- ttest.training$conf.int[2]
  bias.testing.lower <- ttest.testing$conf.int[1]
  bias.testing.upper <- ttest.testing$conf.int[2]

  tilttest.training <- cor.test(df.y.training$Yhat-df.y.training$Y, (df.y.training$Yhat+df.y.training$Y)/2)
  tilttest.testing <- cor.test(df.y.testing$Yhat-df.y.testing$Y, (df.y.testing$Yhat+df.y.testing$Y)/2)
  rtilt.training <- tilttest.training$estimate
  rtilt.testing <- tilttest.testing$estimate
  rtilt.training.lower <- tilttest.training$conf.int[1]
  rtilt.training.upper <- tilttest.training$conf.int[2]
  rtilt.testing.lower <- tilttest.testing$conf.int[1]
  rtilt.testing.upper <- tilttest.testing$conf.int[2]
  rtilt.training.p <- tilttest.training$p.value
  rtilt.testing.p <- tilttest.testing$p.value

  df.stats <- data.frame(subset = c("training set", "testing set"),
                         model = rep(results$title, 2),
                         r = c(r.training, r.testing),
                         r.lower = c(r.lower.training, r.lower.testing),
                         r.upper = c(r.upper.training, r.upper.testing),
                         r.p = c(r.training.p, r.testing.p),
                         RMSE = c(rmse.training, rmse.testing),
                         RMSE.lower = c(rmse.training.lower, rmse.testing.lower),
                         RMSE.upper = c(rmse.training.upper, rmse.testing.upper),
                         bias = c(ttest.training$estimate, ttest.testing$estimate),
                         bias.lower = c(bias.training.lower, bias.testing.lower),
                         bias.upper = c(bias.training.upper, bias.testing.upper),
                         bias.p = c(ttest.training$p.value, ttest.testing$p.value),
                         upper = c(ttest.training$estimate+1.96*sd(df.y.training$error,na.rm=T), ttest.testing$estimate+1.96*sd(df.y.testing$error,na.rm=T)),
                         lower = c(ttest.training$estimate-1.96*sd(df.y.training$error,na.rm=T), ttest.testing$estimate-1.96*sd(df.y.testing$error,na.rm=T)),
                         r.tilt = c(rtilt.training, rtilt.testing),
                         r.tilt.lower = c(rtilt.training.lower, rtilt.testing.lower),
                         r.tilt.upper = c(rtilt.training.upper, rtilt.testing.upper),
                         r.tilt.p = c(rtilt.training.p, rtilt.testing.p)
  )

  df.y <- rbind(df.y.training, df.y.testing)
  text.rmse <- paste0("RMSE=",signif(df.stats$RMSE,3))
  text.r <- paste0("R=",signif(df.stats$r,2))
  text.p <- paste0("p=",signif(df.stats$r.p,2))
  text.bias <- paste0("bias=",signif(df.stats$bias,2))
  text.biasp <- paste0("p=",signif(df.stats$bias.p,2))

  ret <- list(df.y=df.y, df.stats=df.stats)
  return(ret)
}
