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

# Data cleaning

data <- read.csv("Data/TRY_USD.csv", stringsAsFactors = F)

data$newDate <- paste0(substr(data$Date, 9, 12), "-",substr(data$Date,1,3), "-" , substr(data$Date,5,6))
substr(data$newDate, 6, 8) <- ifelse(substr(data$newDate, 6, 8) == "Jan", "01",
                                     ifelse(substr(data$newDate, 6, 8) == "Feb", "02",
                                            ifelse(substr(data$newDate, 6, 8) == "Mar", "03",
                                                   ifelse(substr(data$newDate, 6, 8) == "Apr", "04",
                                                          ifelse(substr(data$newDate, 6, 8) == "May", "05",
                                                                 ifelse(substr(data$newDate, 6, 8) == "Jun", "06",
                                                                        ifelse(substr(data$newDate, 6, 8) == "Jul", "07",
                                                                               ifelse(substr(data$newDate, 6, 8) == "Aug", "08",
                                                                                      ifelse(substr(data$newDate, 6, 8) == "Sep", "09",
                                                                                             ifelse(substr(data$newDate, 6, 8) == "Oct", "10",
                                                                                                    ifelse(substr(data$newDate, 6, 8) == "Nov", "11",
                                                                                                           ifelse(substr(data$newDate, 6, 8) == "Dec", "12"," "))))))))))))
data$newDate <- gsub('[a-z]+', '', data$newDate)
data$newDate  <- as.Date(data$newDate)

data <- data[,c(7,2)]
data[,2] <- as.numeric(data[,2])
data <- data[order(data[,1]),]
data <- data[(nrow(data) - 999):nrow(data),]

# Exploration of the Original Series

I <- as.xts(data[,2], order.by = data[,1])
I <- log(I)

theme_set(theme_bw())

autoplot(I) +
  labs("TRY/USD Historical Data",
       y = "TRY/USD",
       x = 'Date') 

hist(I, col = "blue", breaks = 50, main = 'Histogram of TRY/USD (Original)')

dev.off()

stat_org_name = c('minimum', 'maximum','mean', 'variance', 'skewness', 'kurtosis', 'JB', 'KS', 'Q(6)', 'Q*(6)', 'ARCH(4)', 'ADF', 'KPSS')
stat_org_val = as.numeric(c(min(I),
                            max(I),
                            mean(I), 
                            var(I), 
                            skewness(I), 
                            kurtosis(I),
                            as.numeric(tseries::jarque.bera.test(I)$statistic,
                            ks.test(I, "pnorm", mean(I), sd(I))$statistic),
                            Box.test(I, 6)$statistic, 
                            Box.test(I^2, 6)$statistic,
                            arch.test(arima(I, order = c(3,1,0)))[1,4]),
                            tseries::adf.test(I)$statistic,
                            tseries::kpss.test(I)$statistic
                          )

stat_org_val <- round(as.numeric(stat_org_val), 5)
stat_org_pvalues = c('','','','','','',
                     jarque.bera.test(I)$p.value,
                     ks.test(I, "pnorm", mean(I), sd(I))$p.value,
                     Box.test(I, 6)$p.value,
                     Box.test(I^2, 6)$p.value,
                     arch.test(arima(I, order = c(3,1,0)))[1,5],
                     paste('>', tseries::adf.test(I)$p.value),
                     paste('<', tseries::kpss.test(I)$p.value)
                    )

stat_org_pvalues <- ifelse(stat_org_pvalues == '0','< 1.11e-16', stat_org_pvalues)
stat_org = cbind(stat_org_name, stat_org_val, stat_org_pvalues)
stat_org <- as.data.frame(stat_org)
colnames(stat_org) <- c('Original Series', 'Statistic', 'P-value')
rownames(stat_org) <- NULL
dev.off()

grid.table(stat_org, rows = NULL , theme = ttheme_minimal())

qqnorm(I,
       main = "QQ-plot of Original Series")
qqline(I,
       col = "red")
dev.off()


#Differenced series
y <- na.omit(diff(I))

autoplot(y) +
  labs("TRY/USD Historical Data, Differenced Series",
       y = "TRY/USD (differenced)",
       x = 'Date') 
y[which(y = max(y))]
min(y)

hist(I, col = "blue", breaks = 50, main = 'Histogram of TRY/USD (differenced)')

dev.off()

stat_diff_name = c('minimum','maximum','mean', 'variance', 'skewness', 'kurtosis', 'JB', 'KS', 'Q(6)', 'Q*(6)', 'ARCH(4)', 'ADF', 'KPSS')
stat_diff_val = as.numeric(c(min(y),
                             max(y),
                             mean(y), 
                             var(y), 
                             skewness(y), 
                             kurtosis(y),
                             as.numeric(tseries::jarque.bera.test(y)$statistic,
                                       ks.test(y, "pnorm", mean(y), sd(y))$statistic),
                             Box.test(y, 6)$statistic, 
                             Box.test(y^2, 6)$statistic,
                             arch.test(arima(y, order = c(3,0,0)))[1,4]),
                             tseries::adf.test(y)$statistic,
                             tseries::kpss.test(y)$statistic
)

stat_diff_val <- round(as.numeric(stat_diff_val), 5)
stat_diff_pvalues = c('','','','','','',
                     jarque.bera.test(y)$p.value,
                     ks.test(y, "pnorm", mean(y), sd(y))$p.value,
                     round(Box.test(y, 6)$p.value,2),
                     Box.test(y^2, 6)$p.value,
                     arch.test(arima(y, order = c(3,1,0)))[1,5],
                     paste('<', tseries::adf.test(y)$p.value),
                     paste('>', tseries::kpss.test(y)$p.value)
)


stat_diff_pvalues <- ifelse(stat_diff_pvalues == '0','< 1.11e-16', stat_diff_pvalues)
stat_diff = cbind(stat_diff_name, stat_diff_val, stat_diff_pvalues)
stat_diff <- as.data.frame(stat_diff)
colnames(stat_diff) <- c('Differenced Series', 'Statistic', 'P-value')
rownames(stat_diff) <- NULL
dev.off()

grid.table(stat_diff, rows = NULL , theme = ttheme_minimal())
dev.off()

qqnorm(y,
       main = "QQ-plot of the Differenced Series")
qqline(y,
       col = "red")

dev.off()


acf(y)
pacf(y)

acf(I)
pacf(I)


closer.look <- y[1:699]
plot(closer.look, type = "l")

#saveRDS(y, 'differenced_series.RData')



### Model results ###

source("1_Real_data_models.R")

real_actual <- (y[(length(y) - 59):length(y),] - mean(y))^2
real_forcst <- temp

# MAE and DA
for(i in 1:length(real_forcst)){
  final_mae[[i]] <- mae(real_actual, real_forcst[[i]])
  names(final_mae)[[i]] <- names(real_forcst)[i]
  
  final_da[[i]] <- mean(sign(diff(as.ts(real_actual), lag = 1)) == sign(diff(as.ts(real_forcst[[i]]), lag = 1)), na.rm = T)   
  names(final_da)[[i]] <- names(real_forcst)[i]
}

# Residuals and DM

residuals_actual <- list()

              for(i in 1:20){
                
                residuals_actual[[i]] <- real_actual - real_forcst[[i]]
                names(real_forcst)[i] <- names(real_forcst)[i]
                residuals_actual[[i]] <- as.numeric(residuals_actual[[i]])
                
              }


dm_actual_org = list()
        
              for(g in c("movingaverage","Jordan","SVM", "GARCH_11","EGARCH_11", "TGARCH_11", "GJRGARCH_11")){
                
                for(k in c("GJRGARCH_11","TGARCH_11","EGARCH_11", "GARCH_11","SVM", "Jordan", "movingaverage")){
                  
                  if(g == k){ 
                    
                    dm_actual_org[[g]] = list(0) 
                    
                  } else {
                    
                    dm_actual_org[[g]][[k]] <- as.numeric(dm.test(residuals_actual[[g]], residuals_actual[[k]], alternative = "less")$p.value)
                    
                  }
                }
              }



dm_actual_org_rev = list()
              
              for(g in c("GJRGARCH_11","TGARCH_11","EGARCH_11", "GARCH_11","SVM", "Jordan", "movingaverage")){
                
                for(k in c("movingaverage","Jordan","SVM", "GARCH_11","EGARCH_11", "TGARCH_11", "GJRGARCH_11")){
                  
                  if(g == k){ 
                    
                    dm_actual_org_rev[[g]] = list(0) 
                    
                  } else {
                    
                    dm_actual_org_rev[[g]][[k]] <- as.numeric(dm.test(residuals_actual[[g]], residuals_actual[[k]], alternative = "less")$p.value)
                  }
                }
              }

