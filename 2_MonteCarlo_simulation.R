setwd("/Users/apple/Documents/Ders dost/Master/2nd Semester/Special Topics in Econometrics/Paper/R")
library(rugarch)
library(MLmetrics)
library(Metrics)
library(RSNNS)
library(e1071)
library(forecast)
library(sandwich)

# ARMA(1,0) - GARCH(1,1) simulations following the paper

K            <- 0.0005
alpha        <- 0.1
delta        <- 0.8
theta        <- 0.5 
c            <- 0


MonteCarlo50 <-    function(dist, period, shape, K, alpha, delta, c, theta){
                        
                        sim <- data.frame(matrix(NA, nrow = 1000, ncol = 50))
                        names(sim) <- paste0('sim_', 1:50)
                        for(i in 1:50){
                         
                          sim.spec.norm   <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                                                        mean.model = list(armaOrder = c(1,0), include.mean = FALSE),
                                                        distribution.model = dist, 
                                                        fixed.pars         = list(shape = shape, omega = K, alpha1 = alpha, 
                                                                                  beta1 = delta, ar1 = theta, mu = c))
                                            
                                            path.sgarch <- ugarchpath(sim.spec.norm, n.sim = period, n.start = 1)
                                            sim[,i]    <- as.vector(fitted(path.sgarch))
                          }
                        return(sim)
                      }

sim1000_norm <- MonteCarlo50(dist = 'norm', 
                             period = 1000,
                             shape = NA,
                             K = K,
                             alpha = alpha,
                             delta = delta,
                             c = 0,
                             theta = theta)

sim500_norm <- MonteCarlo50(dist = 'norm', 
                             period = 500,
                             shape = NA,
                             K = K,
                             alpha = alpha,
                             delta = delta,
                             c = 0,
                             theta = theta)

sim1000_std <- MonteCarlo50(dist = 'std', 
                            period = 1000,
                            shape = 5.2,
                            K = K,
                            alpha = alpha,
                            delta = delta, 
                            c = 0,
                            theta = theta)

sim500_std <- MonteCarlo50(dist = 'std', 
                            period = 500,
                            shape = 5.2,
                            K = K,
                            alpha = alpha,
                            delta = delta,
                            c = 0,
                            theta = theta)


# Sourcing the models for Monte Carlo simulation
# Saving the forecasts and test sets (actual) of the simulated data
source("2_MonteCarlo_models.R")

forcst_s500n <- list()
actual_s500n <- list()

                    for( i in 1:length(sim500_norm)){ 
                      
                      forcst_s500n[[i]]  <- models(sim500_norm[,i])
                      
                      y <- sim500_norm[,i]
                      arma_100 <- arima(y,order  = c(1, 0, 0))
                      u <- arma_100$residuals
                      
                      u_test <- u[(length(u)-59):length(u)]
                      u_test_squared <- u_test^2
                      
                      actual_s500n[[i]] <- u_test^2
                    
                    }

forcst_s500s <- list()
actual_s500s <- list()

                    for( i in 1:length(sim500_std)){
                        
                    forcst_s500s[[i]]  <- models(sim500_std[,i], dist = "std", shape = 5.2)
                    
                    
                    y <- sim500_std[,i]
                    arma_100 <- arima(y,order  = c(1, 0, 0))
                    u <- arma_100$residuals
                    
                    
                    u_test <- u[(length(u)-59):length(u)]
                    u_test_squared <- u_test^2
                    
                    actual_s500s[[i]] <- u_test^2
                    }

forcst_s1000s <- list()
actual_s1000s <- list()

                    for( i in 1:length(sim1000_std)){
                        
                    forcst_s1000s[[i]] <- garchmodels(sim1000_std[,i], dist = "std", shape = 5.2)
                    
                    y <- sim1000_std[,i]
                    arma_100 <- arima(y,order  = c(1, 0, 0))
                    u <- arma_100$residuals
                    
                    
                    u_test <- u[(length(u)-59):length(u)]
                    u_test_squared <- u_test^2
                    
                    actual_s1000s[[i]] <- u_test^2
                    }

forcst_s1000n <- list()
actual_s1000n <- list()

                    for( i in 1:length(sim1000_norm)){
                        
                    forcst_s1000n[[i]] <- garchmodels(sim1000_norm[,i])
                    
                    y <- sim1000_norm[,i]
                    arma_100 <- arima(y,order  = c(1, 0, 0))
                    u <- arma_100$residuals
                    
                    
                    u_test <- u[(length(u)-59):length(u)]
                    u_test_squared <- u_test^2
                    
                    actual_s1000n[[i]] <- u_test^2
                    }

# Calculating MAE and DA

mae_s1000n <- list()
da_s1000n <- list()

                    for(i in 1:length(sim1000_norm)){
                      
                      mae_s1000n[[i]] <- list()
                      da_s1000n[[i]] <- list()
                      
                      for(j in 1:length(forcst_s1000n[[i]])){
                        
                        mae_s1000n[[i]][[j]] <- mae(actual_s1000n[[i]],forcst_s1000n[[i]][[j]])
                        names(mae_s1000n[[i]])[j] <- names(forcst_s1000n[[i]])[j]
                        da_s1000n[[i]][[j]] <- mean(sign(diff(actual_s1000n[[i]], lag = 1)) == sign(diff(forcst_s1000n[[i]][[j]], lag = 1)), na.rm = T)
                        names(da_s1000n[[i]])[j] <- names(forcst_s1000n[[i]])[j]
                        
                      }
                    }

mae_s500n <- list()
da_s500n <- list()

                  for(i in 1:length(sim500_norm)){
                    
                    mae_s500n[[i]] <- list()
                    da_s500n[[i]] <- list()
                    
                    for(j in 1:length(forcst_s500n[[i]])){
                      
                      mae_s500n[[i]][[j]] <- mae(actual_s500n[[i]],forcst_s500n[[i]][[j]])
                      names(mae_s500n[[i]])[j] <- names(forcst_s500n[[i]])[j]
                      da_s500n[[i]][[j]] <- mean(sign(diff(actual_s500n[[i]], lag = 1)) == sign(diff(forcst_s500n[[i]][[j]], lag = 1)), na.rm = T)
                      names(da_s500n[[i]])[j] <- names(forcst_s500n[[i]])[j]
                      
                    }
                  }

mae_s500s <- list()
da_s500s <- list()

                    for(i in 1:length(sim500_std)){
                      
                      mae_s500s[[i]] <- list()
                      da_s500s[[i]] <- list()
                      
                      for(j in 1:length(forcst_s500s[[i]])){
                        
                        mae_s500s[[i]][[j]] <- mae(actual_s500s[[i]],forcst_s500s[[i]][[j]])
                        names(mae_s500s[[i]])[j] <- names(forcst_s500s[[i]])[j]
                        da_s500s[[i]][[j]] <- mean(sign(diff(actual_s500s[[i]], lag = 1)) == sign(diff(forcst_s500s[[i]][[j]], lag = 1)), na.rm = T)
                        names(da_s500s[[i]])[j] <- names(forcst_s500s[[i]])[j]
                        
                      }
                    }

mae_s1000s <- list()
da_s1000s <- list()

                    for(i in 1:length(sim1000_std)){
                      
                      mae_s1000s[[i]] <- list()
                      da_s1000s[[i]] <- list()
                      
                      for(j in 1:length(forcst_s1000s[[i]])){
                        
                        mae_s1000s[[i]][[j]] <- mae(actual_s1000s[[i]],forcst_s1000s[[i]][[j]])
                        names(mae_s1000s[[i]])[j] <- names(forcst_s1000s[[i]])[j]
                        da_s1000s[[i]][[j]] <- mean(sign(diff(actual_s1000s[[i]], lag = 1)) == sign(diff(forcst_s1000s[[i]][[j]], lag = 1)), na.rm = T)
                        names(da_s1000s[[i]])[j] <- names(forcst_s1000s[[i]])[j]
                        
                      }
                    }

# Calculating residuals and DM test

res_s1000n <- list()

                    for(i in c("movingaverage", "Jordan", "ff_SVM","GARCH_11", "EGARCH_11", "TGARCH_11", "GJRGARCH_11")){
                      
                      for(j in 1:50){
                        
                        res_s1000n[[i]][[j]] <- as.numeric(forcst_s1000n[[j]][[i]] - actual_s1000n[[j]])
                        
                      }
                    }

res_s500n <- list()
                    for(i in c("movingaverage", "Jordan", "ff_SVM","GARCH_11", "EGARCH_11", "TGARCH_11", "GJRGARCH_11")){
                      
                      for(j in 1:50){
                        
                        res_s500n[[i]][[j]] <- as.numeric(forcst_s500n[[j]][[i]] - actual_s500n[[j]])
                        
                      }
                    }

res_s500s <- list()
                    for(i in c("movingaverage", "Jordan", "ff_SVM","GARCH_11", "EGARCH_11", "TGARCH_11", "GJRGARCH_11")){
                      
                      for(j in 1:50){
                        
                        res_s500s[[i]][[j]] <- as.numeric(forcst_s500s[[j]][[i]] - actual_s500s[[j]])
                        
                      }
                    }

res_s1000s <- list()

                    for(i in c("movingaverage", "Jordan", "ff_SVM","GARCH_11", "EGARCH_11", "TGARCH_11", "GJRGARCH_11")){
                      
                      for(j in 1:50){
                        
                        res_s1000s[[i]][[j]] <- as.numeric(forcst_s1000s[[j]][[i]] - actual_s1000s[[j]])
                        
                      }
                    }
                    



dm = list()

                    for(g in c("movingaverage","Jordan","ff_SVM", "GARCH_11","EGARCH_11", "TGARCH_11", "GJRGARCH_11")){
                      dm_ma <- list()
                      
                      for(k in c("movingaverage","Jordan","ff_SVM", "GARCH_11","EGARCH_11", "TGARCH_11", "GJRGARCH_11")){
                        
                        if(g == k){ 
                          dm[[g]] = list(0) 
                        } else {
                          dm[[g]][[k]] <- list()
                          
                          dm_ma_j_v1 <- c()
                          dm_ma_j_v2 <- c()
                          dm_ma_j_v3 <- c()
                          dm_ma_j_v4 <- c()
                          
                          for(i in 1:50){
                            
                            dm_ma_j_v1[i] <- as.numeric(dm.test(res_s500n[[g]][[i]], res_s500n[[k]][[i]], alternative = "less")$p.value)
                            dm_ma_j_v2[i] <- as.numeric(dm.test(res_s500s[[g]][[i]], res_s500s[[k]][[i]], alternative = "less")$p.value)
                            dm_ma_j_v3[i] <- as.numeric(dm.test(res_s1000s[[g]][[i]], res_s1000s[[k]][[i]], alternative = "less")$p.value)
                            dm_ma_j_v4[i] <- as.numeric(dm.test(res_s1000n[[g]][[i]], res_s1000n[[k]][[i]], alternative = "less")$p.value)
                            
                          }
                          dm[[g]][[k]][["s500n"]] <- mean(dm_ma_j_v1, na.rm = T)
                          dm[[g]][[k]][["s500s"]] <- mean(dm_ma_j_v2, na.rm = T)
                          dm[[g]][[k]][["s1000s"]] <- mean(dm_ma_j_v3, na.rm = T)
                          dm[[g]][[k]][["s1000n"]] <- mean(dm_ma_j_v4, na.rm = T)
                        }
                        
                      }
                    }

dmrev = list()

                    for(g in c("movingaverage","Jordan","ff_SVM", "GARCH_11","EGARCH_11", "TGARCH_11", "GJRGARCH_11")){
                      
                      
                      for(k in c("GJRGARCH_11","TGARCH_11","EGARCH_11", "GARCH_11","ff_SVM", "Jordan", "movingaverage")){
                        
                        if(g == k){ 
                          dmrev[[g]] = list(0) 
                        } else {
                          dmrev[[g]][[k]] <- list()
                          
                          dm_ma_j_v1 <- c()
                          dm_ma_j_v2 <- c()
                          dm_ma_j_v3 <- c()
                          dm_ma_j_v4 <- c()
                          
                          for(i in 1:50){
                            
                            dm_ma_j_v1[i] <- as.numeric(dm.test(res_s500n[[g]][[i]], res_s500n[[k]][[i]], alternative = "less")$p.value)
                            dm_ma_j_v2[i] <- as.numeric(dm.test(res_s500s[[g]][[i]], res_s500s[[k]][[i]], alternative = "less")$p.value)
                            dm_ma_j_v3[i] <- as.numeric(dm.test(res_s1000s[[g]][[i]], res_s1000s[[k]][[i]], alternative = "less")$p.value)
                            dm_ma_j_v4[i] <- as.numeric(dm.test(res_s1000n[[g]][[i]], res_s1000n[[k]][[i]], alternative = "less")$p.value)
                            
                          }
                          dmrev[[g]][[k]][["s500n"]] <- mean(dm_ma_j_v1, na.rm = T)
                          dmrev[[g]][[k]][["s500s"]] <- mean(dm_ma_j_v2, na.rm = T)
                          dmrev[[g]][[k]][["s1000s"]] <- mean(dm_ma_j_v3, na.rm = T)
                          dmrev[[g]][[k]][["s1000n"]] <- mean(dm_ma_j_v4, na.rm = T)
                        }
                        
                      }
                    }
