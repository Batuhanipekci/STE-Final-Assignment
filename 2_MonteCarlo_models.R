models = function(y, dist = "norm", shape = NA)  {
  
  require(rugarch)
  require(forecast)
  require(e1071)
  require(quantmod)
  require(timeSeries)
  require(Metrics)
  require(RSNNS)
  require(smooth)
  require(Mcomp)
  
                              forecasts <- list()

                              arma_101 <- arima(y,order  = c(1, 0, 1))
                              u <- arma_101$residuals
                              
                              training <- y[1:(length(y)-60)]
                              
                              # Moving Average
                                u.training <- u[1:(length(u)-60)] 
                                forecasts[["movingaverage"]]  <- sma(u.training^2, h = 60, order = 5)$forecast
                              
                              # GARCH(1,1)
                              
                              spec_GARCH_11 <- ugarchspec(variance.model = list(model = "sGARCH",
                                                                                garchOrder = c(1,1),
                                                                                submodel = NULL,
                                                                                external.regressors = NULL,
                                                                                variance.targeting = FALSE),
                                                          distribution.model = dist, 
                                                          fixed.pars = list(shape = shape),
                                                          mean.model     = list(armaOrder = c(3,1),
                                                                                external.regressors = NULL) 
                              )
                              
                              garch_11 <- ugarchfit(spec = spec_GARCH_11 , data = training, solver.control = list(tol = 1e-12))
                              forecast_garch_11 <- ugarchforecast(garch_11, n.ahead = 60)
                              forecasts[["GARCH_11"]] <- forecast_garch_11@forecast$sigmaFor
                                 
                              # EGARCH(1,1)
                              
                              spec_EGARCH_11 <- ugarchspec(variance.model = list(model = "eGARCH",
                                                                                 garchOrder = c(1,1),
                                                                                 submodel = NULL,
                                                                                 external.regressors = NULL,
                                                                                 variance.targeting = FALSE),
                                                           distribution.model = dist,
                                                           fixed.pars = list(shape = shape),
                                                           mean.model     = list(armaOrder = c(3,1),
                                                                                 external.regressors = NULL) 
                              )
                              
                              egarch_11 <- ugarchfit(spec = spec_EGARCH_11 , data = training, solver.control = list(tol = 1e-12))
                              forecast_egarch_11 <- ugarchforecast(egarch_11, n.ahead = 60)
                              forecasts[["EGARCH_11"]] <- forecast_egarch_11@forecast$sigmaFor
                              
                              # GJRGARCH(1,1)
                              
                              spec_GJRGARCH_11 <- ugarchspec(variance.model = list(model = "gjrGARCH",
                                                                                   garchOrder = c(1,1),
                                                                                   submodel = NULL,
                                                                                   external.regressors = NULL,
                                                                                   variance.targeting = FALSE),
                                                             distribution.model = dist, 
                                                             fixed.pars = list(shape = shape),
                                                             mean.model     = list(armaOrder = c(3,1),
                                                                                   external.regressors = NULL)
                                                             )
                              
                              gjrgarch_11 <- ugarchfit(spec = spec_GJRGARCH_11 , data = training, solver.control = list(tol = 1e-12))
                              forecast_gjrgarch_11 <- ugarchforecast(gjrgarch_11, n.ahead = 60)

                              forecasts[["GJRGARCH_11"]] <- forecast_gjrgarch_11@forecast$sigmaFor
                              
                              # TGARCH (1,1)
                              spec_TGARCH_11 <- ugarchspec(variance.model = list(model = "fGARCH",
                                                                                 garchOrder = c(1,1),
                                                                                 submodel = "TGARCH",
                                                                                 external.regressors = NULL,
                                                                                 variance.targeting = FALSE),
                                                           distribution.model = dist,
                                                           fixed.pars = list(shape = shape),
                                                           mean.model     = list(armaOrder = c(3,1),
                                                                                 external.regressors = NULL)
                                                           )
                        
                              tgarch_11 <- ugarchfit(spec = spec_TGARCH_11 , data = training, solver.control = list(tol = 1e-12))
                              forecast_tgarch_11 <- ugarchforecast(tgarch_11, n.ahead = 60)
                              forecasts[["TGARCH_11"]] <- forecast_tgarch_11@forecast$sigmaFor
                              
                              sma <- sma(training, h = 60, order = 5)
                              forecasts[["moving_average"]] <- sma$forecast
                              
                              # SVM
                              
                                
                              x <- as.data.frame(na.omit(cbind(as.ts(y), lag(as.ts(y)))))
                              names(x) <- c("y_t", "y_t.1")
                              
                              c.mean <- svm(x = x[,-1], y = x[,1], scale = FALSE,
                                            type = "eps-regression", kernel = "radial",
                                            cost = 0.05, epsilon = 0.0001, gamma = 12.5 )
                             
                              
                              f <- list()
                              
                              for(n in 60:1){
                                
                                u.svm <- as.ts(c.mean$residuals)
                                u_tr.svm <- u.svm[(61-n):(length(u.svm) - n)]  #Subset to training set
                                u_tr.1.svm <- na.omit(lag(u_tr.svm))
                                u_tr_squared.svm <- as.data.frame(na.omit(cbind(u_tr.svm, u_tr.1.svm))^2)  #Final training set
                                names(u_tr_squared.svm) <- c("y_t", "y_t.1")
                                

                              
                                  svm     <- svm(x = as.numeric(u_tr_squared.svm[,-1]), y = as.numeric(u_tr_squared.svm[,1]), 
                                                  scale = FALSE,
                                                  type = "eps-regression", kernel = "radial",
                                                  cost = 10, epsilon = 0.00005, gamma = 50 )
                                
                                forecast_svmr <- predict(svm, as.data.frame(u_tr_squared.svm[,-1]))
                                
                                f[[61-n]] <- forecast_svmr[length(forecast_svmr)]
                              }
                              
                              forecasts[["ff_SVM"]] <- unlist(f)
                              
                              # RNN
                              x <- as.data.frame(na.omit(cbind(as.ts(y), lag(as.ts(y)))))
                              names(x) <- c("y_t", "y_t.1")
                              
                              rnntr <- jordan(x[,-1], x[,1],
                                              size = 4, 
                                              sigmoid = 'tanh',
                                              update_rule = 'sgd',
                                              learningrate = 0.05)
                              
                              fit.j = predict(rnntr, as.data.frame(x[,-1]))
                              
                              u.j = y[-1] - fit.j
                              u_tr.j <- u.j[(61-n):(length(u.j) - n)]  #Subset to training set
                              u_tr.1.j <- na.omit(lag(u_tr.j))
                              u_tr_squared.j <- as.data.frame(na.omit(cbind(u_tr.j, u_tr.1.j))^2)  #Final training set
                              names(u_tr_squared.j) <- c("y_t", "y_t.1")
                              
                              #Jordan 
                              
                              j <- list()
                              
                              for(n in 60:1){
                                
                                u.j = y[-1] - fit.j
                                u_tr.j <- u.j[(61-n):(length(u.j) - n)]  #Subset to training set
                                u_tr.1.j <- na.omit(lag(u_tr.j))
                                u_tr_squared.j <- as.data.frame(na.omit(cbind(u_tr.j, u_tr.1.j))^2)  #Final training set
                                names(u_tr_squared.j) <- c("y_t", "y_t.1")
                                
                                
                                rnntru <- jordan(u_tr_squared.j[,-1], u_tr_squared.j[,1], 
                                                 size = 4, 
                                                 sigmoid = 'tanh',
                                                 update_rule = 'sgd',
                                                 learningrate = 0.05)

                                forecast_j <- predict(rnntru, as.data.frame(u_tr_squared.j[,-1]))
                                
                                j[[61-n]] <- forecast_j[length(forecast_j)]
                              }
                              
                              forecasts[["Jordan"]] <- unlist(j)
                              
                              
                      return(forecasts)
}
                      
                      