setwd("/Users/apple/Documents/Ders dost/Master/2nd Semester/Special Topics in Econometrics/Paper/R/")

library(e1071)
library(timeSeries)
library(quantmod)
library(tseries)
library(rugarch)
library(aTSA)
library(fpp2)
library(forecast)
library(gridExtra)
library(tseries)
library(astsa)
library(FitAR)
library(MSGARCH)

# 1.Data

y <- readRDS('Data/differenced_series.RData')

y <- (y-mean(y))/sd(y)
acf(y)
pacf(y)

# 2.Information Criteria to choose among ARIMA models

AIC <- list()
BIC <- list()
  
            for(d in 0:2){
              for(q in 0:5){
                  for(p in 1:5){
                        temp <- Arima(y,order  = c(p, d, q))
                        AIC[[paste('p =', p,',','d =', d, ',' ,'q =', q)]] <- temp$aic
                        BIC[[paste('p =', p,',','d =', d,',' ,'q =', q)]] <- temp$bic
                       
                   }
            
               }
            }

scores <- as.data.frame(cbind(parameters = names(AIC), AIC = unlist(AIC), BIC = unlist(BIC)), stringsAsFactors = F)
row.names(scores) = 1:90
scores$AIC <- as.numeric(scores$AIC)
scores$BIC <- as.numeric(scores$BIC)

best <- scores[intersect(order(unique(scores$BIC))[1:15] , order(unique(scores$AIC))[1:15]), ]
best 
## best AIC : (1,0,5), (3,0,0), (5,0,1)
## best BIC : (3,0,0), (3,0,1), (4,0,0) 
## select: (3,0,0), (3,0,1)


# 3.Residual diagnostics

par(mar=c(1,1,1,1))
dev.off()
arima_300 <- arima(y, order  = c(3, 0, 0))
LBQPlot(arima_300$residuals, lag.max = 15)

arima_301 <- arima(y,order  = c(3, 0, 1))
LBQPlot(arima_301$residuals, lag.max = 15)
## It makes sense to select (3,0,1)
pacf(arima_301$residuals)

# 4.Exploration of standardized residuals from ARIMA(3,0,1)

std_res <- (arima_301$residuals - mean(arima_301$residuals)) / sd(arima_301$residuals)

plot(std_res, main = 'Standardized Residuals (3,0,1)')

qqnorm(y, main = "QQ-plot of the Standardized Residuals from ARIMA(3,0,1)")
qqline(y, col = "red")

Acf(std_res, main = 'ACF, Standardized Residuals from ARIMA(3,0,1)')
Pacf(std_res, main = 'PACF, Standardized Residuals from ARIMA(3,0,1)')

## Check for heteroscedasticity
Box.test(arima_301$residuals^2, 1) #0.00001
LBQPlot(arima_301$residuals^2, lag.max = 15)



# 5.Estimating Conditional Variance
y = wo_y
wo_y = (wo_y - mean(wo_y))/sd(wo_y)

training <- y[1:939]
arima_301 <- arima(y,order  = c(3, 0, 1))

test <- arima_301$residuals^2
test <- test[940:999] 

##################
## GARCH / ARCH ##
##################

results <- list()
temp <- list()

y_est_var <- (y[(length(y) - 59):length(y),] - mean(y))^2

forecasts[["movingaverage"]]  <- sma(y_est_var, h = 60, order = 5)$forecast

# ARCH(1)

spec_ARCH_1 <- ugarchspec(variance.model = list(model = "sGARCH",
                                                  garchOrder = c(1,0),
                                                  submodel = NULL,
                                                  external.regressors = NULL,
                                                  variance.targeting = FALSE),
                            
                            mean.model     = list(armaOrder = c(0,1),
                                                  external.regressors = NULL) 
)

arch_1 <- ugarchfit(spec = spec_ARCH_1 , data = training, solver.control = list(trace=0))
forecast_arch_1 <- ugarchforecast(arch_1, n.ahead = 60)


temp[["ARCH_1"]] <- forecast_arch_1
# ARCH(2)

spec_ARCH_2 <- ugarchspec(variance.model = list(model = "sGARCH",
                                                garchOrder = c(2,0),
                                                submodel = NULL,
                                                external.regressors = NULL,
                                                variance.targeting = FALSE),
                          
                          mean.model     = list(armaOrder = c(0,1),
                                                external.regressors = NULL) 
)

arch_2 <- ugarchfit(spec = spec_ARCH_2 , data = training, solver.control = list(trace=0))
forecast_arch_2 <- ugarchforecast(arch_2, n.ahead = 60)

temp[["ARCH_2"]] <- forecast_arch_2


# ARCH(3)

spec_ARCH_3 <- ugarchspec(variance.model = list(model = "sGARCH",
                                                garchOrder = c(3,0),
                                                submodel = NULL,
                                                external.regressors = NULL,
                                                variance.targeting = FALSE),
                          
                          mean.model     = list(armaOrder = c(0,1),
                                                external.regressors = NULL) 
)

arch_3 <- ugarchfit(spec = spec_ARCH_3 , data = training, solver.control = list(trace=0))
forecast_arch_3 <- ugarchforecast(arch_3, n.ahead = 60)

temp[["ARCH_3"]] <- forecast_arch_3


# ARCH(4)

spec_ARCH_4 <- ugarchspec(variance.model = list(model = "sGARCH",
                                                garchOrder = c(4,0),
                                                submodel = NULL,
                                                external.regressors = NULL,
                                                variance.targeting = FALSE),
                          
                          mean.model     = list(armaOrder = c(0,1),
                                                external.regressors = NULL) 
)

arch_4 <- ugarchfit(spec = spec_ARCH_4 , data = training, solver.control = list(trace=0))
forecast_arch_4 <- ugarchforecast(arch_4, n.ahead = 60)


temp[["ARCH_4"]] <- forecast_arch_4

# GARCH(1,1)

spec_GARCH_11 <- ugarchspec(variance.model = list(model = "sGARCH",
                                                       garchOrder = c(1,1),
                                                       submodel = NULL,
                                                       external.regressors = NULL,
                                                       variance.targeting = FALSE),
                               
                                  mean.model     = list(armaOrder = c(0,1),
                                                        external.regressors = NULL) 
                                 )

garch_11 <- ugarchfit(spec = spec_GARCH_11 , data = training, solver.control = list(trace=0))
forecast_garch_11 <- ugarchforecast(garch_11, n.ahead = 60)

temp[["GARCH_11"]] <- forecast_garch_11

# GARCH(2,1)

spec_GARCH_21 <- ugarchspec(variance.model = list(model = "sGARCH",
                                                  garchOrder = c(2,1),
                                                  submodel = NULL,
                                                  external.regressors = NULL,
                                                  variance.targeting = FALSE),
                            
                            mean.model     = list(armaOrder = c(0,1),
                                                  external.regressors = NULL) 
)

garch_21 <- ugarchfit(spec = spec_GARCH_21 , data = training, solver.control = list(trace=0))
forecast_garch_21 <- ugarchforecast(garch_21, n.ahead = 60)

temp[["GARCH_21"]] <- forecast_garch_21

# GARCH(2,2) # did not converge

# EGARCH(1,1)

spec_EGARCH_11 <- ugarchspec(variance.model = list(model = "eGARCH",
                                                  garchOrder = c(1,1),
                                                  submodel = NULL,
                                                  external.regressors = NULL,
                                                  variance.targeting = FALSE),
                            
                            mean.model     = list(armaOrder = c(0,1),
                                                  external.regressors = NULL) 
)

egarch_11 <- ugarchfit(spec = spec_EGARCH_11 , data = training, solver.control = list(trace=0))
forecast_egarch_11 <- ugarchforecast(egarch_11, n.ahead = 60)

temp[["EGARCH_11"]] <- forecast_egarch_11

# EGARCH(2,1)

spec_EGARCH_21 <- ugarchspec(variance.model = list(model = "eGARCH",
                                                   garchOrder = c(2,1),
                                                   submodel = NULL,
                                                   external.regressors = NULL,
                                                   variance.targeting = FALSE),
                             
                             mean.model     = list(armaOrder = c(0,1),
                                                   external.regressors = NULL) 
)

egarch_21 <- ugarchfit(spec = spec_EGARCH_21 , data = training, solver.control = list(trace=0))
forecast_egarch_21 <- ugarchforecast(egarch_21, n.ahead = 60)

temp[["EGARCH_21"]] <- forecast_egarch_21

# EGARCH(2,2)

spec_EGARCH_22 <- ugarchspec(variance.model = list(model = "eGARCH",
                                                   garchOrder = c(2,2),
                                                   submodel = NULL,
                                                   external.regressors = NULL,
                                                   variance.targeting = FALSE),
                             
                             mean.model     = list(armaOrder = c(0,1),
                                                   external.regressors = NULL) 
)

egarch_22 <- ugarchfit(spec = spec_EGARCH_22 , data = training, solver.control = list(trace=0))
forecast_egarch_22 <- ugarchforecast(egarch_22, n.ahead = 60)

temp[["EGARCH_22"]] <- forecast_egarch_22

# GJRGARCH(1,1)

spec_GJRGARCH_11 <- ugarchspec(variance.model = list(model = "gjrGARCH",
                                                   garchOrder = c(1,1),
                                                   submodel = NULL,
                                                   external.regressors = NULL,
                                                   variance.targeting = FALSE),
                             
                             mean.model     = list(armaOrder = c(0,1),
                                                   external.regressors = NULL) 
)

gjrgarch_11 <- ugarchfit(spec = spec_GJRGARCH_11 , data = training, solver.control = list(trace=0))
forecast_gjrgarch_11 <- ugarchforecast(gjrgarch_11, n.ahead = 60)

temp[["GJRGARCH_11"]] <- forecast_gjrgarch_11

# GJRGARCH(1,2)

spec_GJRGARCH_12 <- ugarchspec(variance.model = list(model = "gjrGARCH",
                                                     garchOrder = c(1,2),
                                                     submodel = NULL,
                                                     external.regressors = NULL,
                                                     variance.targeting = FALSE),
                               
                               mean.model     = list(armaOrder = c(0,1),
                                                     external.regressors = NULL) 
)

gjrgarch_12 <- ugarchfit(spec = spec_GJRGARCH_12 , data = training, solver.control = list(trace=0))
forecast_gjrgarch_12 <- ugarchforecast(gjrgarch_12, n.ahead = 60)

temp[["GJRGARCH_12"]] <- forecast_gjrgarch_12

# GJRGARCH(2,1)

spec_GJRGARCH_21 <- ugarchspec(variance.model = list(model = "gjrGARCH",
                                                     garchOrder = c(2,1),
                                                     submodel = NULL,
                                                     external.regressors = NULL,
                                                     variance.targeting = FALSE),
                               
                               mean.model     = list(armaOrder = c(0,1),
                                                     external.regressors = NULL) 
)

gjrgarch_21 <- ugarchfit(spec = spec_GJRGARCH_21 , data = training, solver.control = list(trace=0))
forecast_gjrgarch_21 <- ugarchforecast(gjrgarch_21, n.ahead = 60)

temp[["GJRGARCH_21"]] <- forecast_gjrgarch_21

# GJRGARCH (2,2)

spec_GJRGARCH_22 <- ugarchspec(variance.model = list(model = "gjrGARCH",
                                                     garchOrder = c(2,2),
                                                     submodel = NULL,
                                                     external.regressors = NULL,
                                                     variance.targeting = FALSE),
                               
                               mean.model     = list(armaOrder = c(0,1),
                                                     external.regressors = NULL) 
)

gjrgarch_22 <- ugarchfit(spec = spec_GJRGARCH_22 , data = training, solver.control = list(trace=0))
forecast_gjrgarch_22 <- ugarchforecast(gjrgarch_22, n.ahead = 60)

temp[["GJRGARCH_22"]] <- forecast_gjrgarch_22

# TGARCH (1,1)

spec_TGARCH_11 <- ugarchspec(variance.model = list(model = "fGARCH",
                                                     garchOrder = c(1,1),
                                                     submodel = "TGARCH",
                                                     external.regressors = NULL,
                                                     variance.targeting = FALSE),
                               
                               mean.model     = list(armaOrder = c(0,1),
                                                     external.regressors = NULL) 
)

tgarch_11 <- ugarchfit(spec = spec_TGARCH_11 , data = training, solver.control = list(trace=0))
forecast_tgarch_11 <- ugarchforecast(tgarch_11, n.ahead = 60)

temp[["TGARCH_11"]] <- forecast_tgarch_11

# TGARCH (1,2)

spec_TGARCH_12 <- ugarchspec(variance.model = list(model = "fGARCH",
                                                   garchOrder = c(1,2),
                                                   submodel = "TGARCH",
                                                   external.regressors = NULL,
                                                   variance.targeting = FALSE),
                             
                             mean.model     = list(armaOrder = c(0,1),
                                                   external.regressors = NULL) 
)

tgarch_12 <- ugarchfit(spec = spec_TGARCH_12 , data = training, solver.control = list(trace=0))
forecast_tgarch_12 <- ugarchforecast(tgarch_12, n.ahead = 60)

temp[["TGARCH_12"]] <- forecast_tgarch_12

# TGARCH (2,1)

spec_TGARCH_21 <- ugarchspec(variance.model = list(model = "fGARCH",
                                                   garchOrder = c(2,1),
                                                   submodel = "TGARCH",
                                                   external.regressors = NULL,
                                                   variance.targeting = FALSE),
                             
                             mean.model     = list(armaOrder = c(0,1),
                                                   external.regressors = NULL) 
)

tgarch_21 <- ugarchfit(spec = spec_TGARCH_21 , data = training, solver.control = list(trace=0))
forecast_tgarch_21 <- ugarchforecast(tgarch_21, n.ahead = 60)

temp[["TGARCH_21"]] <- forecast_tgarch_21


# TGARCH(2,2)

spec_TGARCH_22 <- ugarchspec(variance.model = list(model = "fGARCH",
                                                   garchOrder = c(2,2),
                                                   submodel = "TGARCH",
                                                   external.regressors = NULL,
                                                   variance.targeting = FALSE),
                             
                             mean.model     = list(armaOrder = c(0,1),
                                                   external.regressors = NULL) 
)

tgarch_22 <- ugarchfit(spec = spec_TGARCH_22 , data = training, solver.control = list(trace=0))
forecast_tgarch_22 <- ugarchforecast(tgarch_22, n.ahead = 60)


temp[["TGARCH_22"]] <- forecast_tgarch_22







#####################
### Recurrent SVM ###
#####################


y <- (y - mean(y))/sd(y)

x <- as.data.frame(na.omit(cbind(as.ts(y), lag(as.ts(y)))))
names(x) <- c("y_t", "y_t.1")

# SVM - tuning

cost1 <- c(0.0005, 0.001, 0.01, 0.05, 0.1, 1, 5, 20)
epsilon1 <- c(10^-5, 5*10^-5, seq(10^-4, 10^-3, 2*10^-4), 10^-3, 5*10^-3, 10^-2, 5*10^-2, seq(10^-3, 10^-2, 10^-3),0.1)
gamma1 <- c(0.0005, 0.005,0.001, 0.01, 0.05, 1, 5, 20, 100)
d1 = c(1,3)


tune <- tune(svm, y_t ~ y_t.1 , data = x, ranges = list(epsilon = epsilon1, cost = cost1, d = d1, kernel = "polynomial"))
tune$best.parameters

c.mean <- svm(x = x[,-1], y = x[,1], scale = FALSE, type = "eps-regression", kernel = "radial", 
              gamma = tune$best.parameters$gamma, cost =  tune$best.parameters$cost, epsilon =  tune$best.parameters$epsilon)

u.svm <- as.ts(c.mean$residuals)
n    <- 59
u_tr.svm <- u.svm[1:(length(u.svm) - n)]  
u_tr.1.svm <- na.omit(lag(u_tr.svm))
u_tr_squared.svm <- as.data.frame(cbind(as.ts(u_tr.svm), lag(as.ts(u_tr.svm)))^2)  #Final training set
names(u_tr_squared.svm) <- c("y_t", "y_t.1")

tune.u <- tune(svm, y_t ~ y_t.1 , data = u_tr_squared.svm, ranges = list(epsilon = epsilon1, cost = cost1, gamma = gamma1))

tune.u$best.parameters


# Recurrent-SVM-Gaussian Variance Estimation

f <- list()

for(n in 60:1){
  
  u.svm <- as.ts(c.mean$residuals)
  u_tr.svm <- u.svm[(61-n):(length(u.svm) - n)]  
  u_tr.1.svm <- na.omit(lag(u_tr.svm))
  u_tr_squared.svm <- as.data.frame(na.omit(cbind(u_tr.svm, u_tr.1.svm))^2)  
  
  names(u_tr_squared.svm) <- c("y_t", "y_t.1")
  
  
            while(p_count < 5){
              
              print(i)
              
              
              
              
              svmr     <- svm(y_t ~ y_t.1 , data = u_tr_squared.svm, scale = FALSE,
                              type = "eps-regression", kernel = "radial",
                              gamma = tune.u$best.parameters$gamma, cost = tune.u$best.parameters$cost, epsilon = tune.u$best.parameters$epsilon)
              
              
              test    <- Box.test(svmr$residuals, lag = 1, type = "Ljung-Box")
              p_val   <- test$p.value
              print(p_val)
              p_count <- ifelse(p_val > 0.1, p_count + 1, 0)
              
              w <- svmr$residuals
              u_tr_squared.svm[,3] <- w
              
              i <- i + 1
            }
            
  
  forecast_svmr <- predict(svmr, as.data.frame(u_tr_squared.svm[,-1]))
  
  f[[61-n]] <- forecast_svmr[length(forecast_svmr)]
}

temp[["SVM"]] <- unlist(f)






#######################
### Jordan networks ###
#######################

breaks <- 10
set.seed(222)
folds <- cut(1:nrow(x), breaks = breaks, labels = FALSE)
nnet.param <- expand.grid("learningrate" = seq(from = 0.01, to = 0.1, by = 0.01))

results <- as.data.frame(matrix(NA, ncol = nrow(nnet.param), nrow = breaks))


          for(n in 1:nrow(nnet.param)){
              for (i in 1:breaks) {
                  # Splitting the data into training and validation
                  
                  idx.val <- which(folds == i, arr.ind = TRUE)
                  cv.train <- x[-idx.val,]
                  cv.val <- x[idx.val,]
                  
                  
                  # Training the rnn model with a number of nodes n (size) and learning rate
                  neuralnet <- jordan(cv.train[,-1], cv.train[,1], size = 4, learningrate =  nnet.param$learningrate[n])
                  yhat.val <- predict(neuralnet, as.data.frame(cv.val[,-1], nrow = 100))
                  results[i, n] <- mse(as.vector(yhat.val),cv.val[,1])
              }
          }


nnet.param$avg_mse <- apply(results, 2, mean)
opt.mse <- which.min(nnet.param$avg_mse)

model_rnn <- jordan(x[,-1], x[,1],
size = 4,
sigmoid = 'tanh',
update_rule = 'sgd',
learningrate =  nnet.param[opt.mse,])

rnn <- predict(model_rnn, as.data.frame(x[,1]))

u <- x[,1] - rnn
j <- list()

          for(t in 60:1){
              print(t)
              
              u_tr <- u[(61-t):(length(u) - t)]
              u_tr.1 <- na.omit(lag(u_tr))
              u_squared <- na.omit(cbind(u_tr, u_tr.1))^2
              
              u_test <- u[(length(u) - n):length(u)]
              u_test_squared <- u_test^2
              
              ### Tune conditional variance
              
              
              breaks <- 10
              set.seed(222)
              folds <- cut(1:nrow(u_squared), breaks = breaks, labels = FALSE)
              nnet.param.u <- expand.grid("learningrate" = seq(from = 0.001, to = 0.1, by = 0.01))
              
              
              results.u <- as.data.frame(matrix(NA, ncol = nrow(nnet.param.u), nrow = breaks))
              
              for(n in 1:nrow(nnet.param.u)){
                  for (i in 1:breaks) {
                      # Splitting the data into training and validation
                      
                      idx.val <- which(folds == i, arr.ind = TRUE)
                      cv.train <- u_squared[-idx.val,]
                      cv.val <- u_squared[idx.val,]
                      
                      # Training the rnn model with a number of nodes n (size) and learning rate
                      neuralnet <- jordan(cv.train[,-1], cv.train[,1], size = 4, learningrate =  nnet.param.u$learningrate[n])
                      
                      yhat.val <- predict(neuralnet, as.data.frame(cv.val[,-1]))
                      results.u[i, n] <- mse(as.vector(yhat.val),cv.val[,1])
                  }
              }
              
              
              nnet.param.u$avg_mse <- apply(results.u, 2, mean)
              opt.mse <- which.min(nnet.param.u$avg_mse)
              
              ###
              
              model_rnn_u <- jordan(u_squared[,-1], u_squared[,1],
              size = 4,
              sigmoid = 'tanh',
              update_rule = 'sgd',
              learningrate = nnet.param[opt.mse,])
              
              ann_garch_forecast <- predict(model_rnn_u, as.data.frame(u_squared[-1]^2))
              
              j[[t]] <- ann_garch_forecast[length(ann_garch_forecast)]
          }


temp[["Jordan"]] <- unlist(j)


#saveRDS(temp, "forecasts.RData")
