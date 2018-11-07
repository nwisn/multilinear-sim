#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(ggsci)
library(MASS)
library(gridExtra)
library(caret)
library(glmnet)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {



  output$figure <- renderPlot({

    # dynamic inputs -----------------------------------------------------------------

    p <- input$p # number of predictors in complete
    frac <- input$frac # fraction of predictors in available data
    psub <- ceiling(frac*p) # number of predictors in available data
    n <- input$n*2 # number of total samples
    epsilon <- input$epsilon # noise
    alpha <- input$alpha # skew
    randomseed <- input$randomseed

    set.seed(randomseed)

    # Generate data -----------------------------------------------------------------
    # specify true model coefficients
    b <- as.matrix(data.frame(b = rnorm(p+1, mean = 0, sd = 1)))

    # simulate X and compute true Y under model
    #X <- cbind(rep(1,n), do.call(rbind, replicate(n, rnorm(p, mean = 0, sd = 1), simplify = F)))
    X <- cbind(rep(1,n), do.call(rbind, replicate(n, rsn(p, xi = 0, omega = 1, alpha = alpha), simplify = F)))
    colnames(X) <- paste0("X", 1:ncol(X))
    Y <- X %*% b + rnorm(n, mean=0, sd=1)*epsilon

    # take subset of predictors as 'available features'
    inTrain <- createDataPartition(y=Y,p=0.5,list=F)
    Xsub.training <- X[inTrain,c(1,2:(psub+1))]
    Xsub.testing <- X[-inTrain,c(1,2:(psub+1))]
    Y.training <- Y[inTrain]
    Y.testing <- Y[-inTrain]

    source("functions.R")

    #--------------------------------------------------------------------
    # compute

    results <- list()
    gg <- list()

    # fit lm on all available features
    results$lm <- run_lm(title="OLS",
                         Xsub.training = Xsub.training,
                         Xsub.testing = Xsub.testing,
                         Y.training = Y.training)
    gg$lm <- get_model_results(results$lm,
                               Xsub.training = Xsub.training,
                               Xsub.testing = Xsub.testing,
                               Y.training = Y.training,
                               Y.testing = Y.testing)

    # fit ridge on all available features
    results$ridge <- run_glmnet(alpha=0,title="ridge",
                                Xsub.training = Xsub.training,
                                Xsub.testing = Xsub.testing,
                                Y.training = Y.training)
    gg$ridge <- get_model_results(results$ridge,
                                  Xsub.training = Xsub.training,
                                  Xsub.testing = Xsub.testing,
                                  Y.training = Y.training,
                                  Y.testing = Y.testing)

    # fit lasso on all available features
    results$lasso <- run_glmnet(alpha=1,title="lasso",
                                Xsub.training = Xsub.training,
                                Xsub.testing = Xsub.testing,
                                Y.training = Y.training)
    gg$lasso <- get_model_results(results$lasso,
                                  Xsub.training = Xsub.training,
                                  Xsub.testing = Xsub.testing,
                                  Y.training = Y.training,
                                  Y.testing = Y.testing)
    # get features from lasso
    beta.lasso <- as.matrix(coef(results$lasso$model$glmnet.fit, s=results$lasso$model[[results$lasso$lambda]]))
    features.lasso <- names(beta.lasso[-1,1])[beta.lasso[-1,1] != 0]

    # fit lm using subset of features
    results$lm.subset <- run_lm(title="OLS-subset", predictors = features.lasso,
                                Xsub.training = Xsub.training,
                                Xsub.testing = Xsub.testing,
                                Y.training = Y.training)
    gg$lm.subset <- get_model_results(results$lm.subset,
                                      Xsub.training = Xsub.training,
                                      Xsub.testing = Xsub.testing,
                                      Y.training = Y.training,
                                      Y.testing = Y.testing)


    # compute inverse density weights
    fd <- approxfun(density(Y, adjust = 2)) # density for inverse weighting
    dens <- fd(Y.training)
    normdens <- dens/max(dens)
    weights <- 1/normdens

    # fit lasso with weights
    results$lasso.w <- run_glmnet(alpha=1,title="weighted lasso", weights = weights,
                                  Xsub.training = Xsub.training,
                                  Xsub.testing = Xsub.testing,
                                  Y.training = Y.training)
    gg$lasso.w <- get_model_results(results$lasso.w,
                                    Xsub.training = Xsub.training,
                                    Xsub.testing = Xsub.testing,
                                    Y.training = Y.training,
                                    Y.testing = Y.testing)
    # get features from lasso
    beta.lasso.w <- as.matrix(coef(results$lasso.w$model$glmnet.fit, s=results$lasso.w$model[[results$lasso.w$lambda]]))
    features.lasso.w <- names(beta.lasso.w[-1,1])[beta.lasso.w[-1,1] != 0]

    # fit lm using subset of features
    results$lm.subset.w <- run_lm(title="weighted OLS-subset", predictors = features.lasso.w, weights = weights,
                                  Xsub.training = Xsub.training,
                                  Xsub.testing = Xsub.testing,
                                  Y.training = Y.training)
    gg$lm.subset.w <- get_model_results(results$lm.subset.w,
                                        Xsub.training = Xsub.training,
                                        Xsub.testing = Xsub.testing,
                                        Y.training = Y.training,
                                        Y.testing = Y.testing)




    #-------------------------------------------------------------------
    # final plot

    df.y.final <- do.call(rbind, lapply(gg, function(this_gg) this_gg$df.y))
    df.stats.final <- do.call(rbind, lapply(gg, function(this_gg) this_gg$df.stats))

    gg.cor <- ggplot(subset(df.y.final, df.y.final$subset=="testing set")) + aes(x = Y, y=Yhat, color = model) +
      geom_point(alpha=min(1/n/2+0.5, 1), size = 1) +
      geom_abline(slope=0,intercept=1,lty=3)+
      geom_vline(xintercept=0,lty=3)+
      geom_abline(data=subset(df.stats.final,df.stats.final$subset=="testing set"),aes(intercept=bias,slope=1),lty=2)+
      geom_abline(data=subset(df.stats.final,df.stats.final$subset=="testing set"),aes(intercept=upper,slope=1),lty=2)+
      geom_abline(data=subset(df.stats.final,df.stats.final$subset=="testing set"),aes(intercept=lower,slope=1),lty=2)+
      facet_grid(~model)+scale_color_d3()+
      ggtitle("Prediction Correspondence in Test Set") +
      guides(color=F)

    gg.error <- ggplot(subset(df.y.final, df.y.final$subset=="testing set")) + aes(x = Y, y=(Yhat-Y), color = model) +
      geom_point(alpha=min(1/n/2+0.5, 1), size = 1) +
      geom_abline(slope=0,intercept=0,lty=3)+
      geom_vline(xintercept=0,lty=3)+
      geom_abline(data=subset(df.stats.final,df.stats.final$subset=="testing set"),aes(intercept=bias,slope=0),lty=2)+
      geom_abline(data=subset(df.stats.final,df.stats.final$subset=="testing set"),aes(intercept=upper,slope=0),lty=2)+
      geom_abline(data=subset(df.stats.final,df.stats.final$subset=="testing set"),aes(intercept=lower,slope=0),lty=2)+
      facet_grid(~model)+scale_color_d3()+
      ggtitle("Prediction Errors in Test Set") +
      guides(color=F)

    gg.meandiff <- ggplot(subset(df.y.final, df.y.final$subset=="testing set")) + aes(x = (Y+Yhat)/2, y=(Yhat-Y), color = model) +
      geom_point(alpha=min(1/n/2+0.5, 1), size = 1) +
      geom_abline(slope=0,intercept=0,lty=3)+
      geom_vline(xintercept=0,lty=3)+
      geom_abline(data=subset(df.stats.final,df.stats.final$subset=="testing set"),aes(intercept=bias,slope=0),lty=2)+
      geom_abline(data=subset(df.stats.final,df.stats.final$subset=="testing set"),aes(intercept=upper,slope=0),lty=2)+
      geom_abline(data=subset(df.stats.final,df.stats.final$subset=="testing set"),aes(intercept=lower,slope=0),lty=2)+
      facet_grid(~model)+scale_color_d3()+
      ggtitle("Mean-Difference (Bland-Altman) Plot") +
      guides(color=F)

    this_df <- subset(df.stats.final, df.stats.final$subset=="testing set")

    this_df$model.ordered <- factor(this_df$model, levels=this_df$model[order(this_df$r, decreasing = F)])
    gg.r <- ggplot(this_df) + aes(x=model.ordered, y=r, fill=model) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept=0,lty=1)+
      geom_errorbar( aes(x=model.ordered, y=r, ymin=r.lower, ymax=r.upper), width=0.4, color = "black", alpha=0.9, size=1.3)+
      scale_fill_d3() +
      coord_flip() +
      ggtitle("R") +
      guides(fill=F) +
      xlab("")

    this_df$model.ordered <- factor(this_df$model, levels=this_df$model[order(this_df$RMSE, decreasing = T)])
    gg.rmse <- ggplot(this_df) + aes(x=model.ordered, y=RMSE, fill=model) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept=0,lty=1)+
      geom_errorbar( aes(x=model.ordered, y=RMSE, ymin=RMSE.lower, ymax=RMSE.upper), width=0.4, color = "black", alpha=0.9, size=1.3)+
      scale_fill_d3() +
      coord_flip() +
      ggtitle("RMSE") +
      guides(fill=F) +
      xlab("")

    this_df$model.ordered <- factor(this_df$model, levels=this_df$model[order(abs(this_df$bias), decreasing = T)])
    gg.bias <- ggplot(this_df) + aes(x=model.ordered, y=bias, fill=model) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept=0,lty=1)+
      geom_errorbar( aes(x=model.ordered, y=bias, ymin=bias.lower, ymax=bias.upper), width=0.4, color = "black", alpha=0.9, size=1.3)+
      scale_fill_d3() +
      coord_flip() +
      ggtitle("Mean bias") +
      guides(fill=F) +
      xlab("")

    this_df$model.ordered <- factor(this_df$model, levels=this_df$model[order(abs(this_df$r.tilt), decreasing = T)])
    gg.tilt <- ggplot(this_df) + aes(x=model.ordered, y=r.tilt, fill=model) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept=0,lty=1)+
      geom_errorbar( aes(x=model.ordered, y=r.tilt, ymin=r.tilt.lower, ymax=r.tilt.upper), width=0.4, color = "black", alpha=0.9, size=1.3)+
      scale_fill_d3() +
      coord_flip() +
      ggtitle("Tilt") +
      guides(fill=F) +
      xlab("")

    lay <- rbind(c(7,7),
                 c(1,1),
                 #c(1,1,1,1),
                 c(2,2),
                 #c(2,2,2,2),
                 c(3,4),
                 #c(3,3,4,4),
                 c(5,6)
                 #c(5,5,6,6)
                 )

    gg.final <- grid.arrange(gg.error, gg.meandiff, gg.r, gg.rmse, gg.bias, gg.tilt, gg.cor, layout_matrix = lay)

    gg.final

  })

})
