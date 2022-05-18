# Necessary packages
library(tictoc)
library(parallel)
library(rugarch)
library(readxl)
library(tseries)
library(xts)
library(plot.matrix)
library(RColorBrewer)
library(binaryLogic)
library(TTR)
library(dplyr)
library(quantmod)

#--------------------------------------------------#
#                  Initialize VIX                  #
#--------------------------------------------------#

getSymbols(c("^VIX"),from = as.Date("2012-02-29"), to = as.Date("2021-02-27"))
data <- VIX$VIX.Close

# Plot data
plot(xts(data), main = "")

# Plot the logarithmic data
returns = log(data)
plot(xts(returns), main = "")

# Plot the differenced data
returns <- na.omit(diff(returns))
plot(xts(returns), main = "")

# Plot the ACF
acf(returns , main = "")

# Test for stationarity with Augmented Dickey-Fuller Test
adf.test(returns)


#--------------------------------------------------#
#                  Initialize RV                   #
#--------------------------------------------------#

data <- read_excel("data_rv_9year.xlsx")

# Plot the data
plot(xts(read.zoo(data)), main = "")

# Plot the logarithmic data
returns = log(read.zoo(data))
plot(xts(returns), main = "")

# Plot the ACF
acf(ts(returns), main = "")

# Test for stationarity with Augmented Dickey-Fuller Test
adf.test(returns)


#--------------------------------------------------#
#              Generate models                     #
#--------------------------------------------------#

# Generate models up to an ARFIMA(4,d,4)
ar = 4
ma = 4
tic()
cl = makePSOCKcluster(detectCores()-1, useXDR = FALSE)
fit = autoarfima(
  data = returns, ar.max = ar, ma.max = ma, cluster = cl, criterion = "BIC",
  method = "partial", arfima = TRUE, include.mean = TRUE, return.all = TRUE
)

stopCluster(cl)
toc()

#--------------------------------------------------#
#                Generate Heatmaps                 #
#--------------------------------------------------#

M <- fit$rank.matrix
M <- as.data.frame(M)

types_list = c(0,1)

xtick = c()
ytick = c()

for(i in 0:as.integer(ar)){
  xtick = c(xtick, paste("p =", i))
}
for(i in 0:as.integer(ma)){
  ytick = c(ytick, paste("q =", i))
}

par(mfrow = c(2,1))
colfunc <- colorRampPalette(c("white", "blue"))

# For loop creating heatmaps
for(elem in types_list) {
  heatmap <-
    dplyr::filter(M, ARFIMA == elem)[-c(length(M), length(M) - 2, length(M) - 3)]
  
  heatmap_matrix <- matrix(nrow = ar+1, ncol = ma+1)
  
  for (i in 1:nrow(heatmap)) {
    row <- heatmap[i, ]
    heatmap_matrix[row$AR+1, row$MA+1] <- row$BIC
  }
  
  heatmap_matrix <- (heatmap_matrix - mean(na.omit(heatmap_matrix))) / sd(na.omit(heatmap_matrix))
  
  heat_min = min(na.omit(heatmap_matrix))
  heat_max = max(na.omit(heatmap_matrix))
  
  N = 30
  
  if(elem == 0){
    k = -0.36 
    mainname = "ARMA(p,q)"
  } else{
    k=-0.2 
    mainname = "ARFIMA(p,d,q)"
  }
  
  if(k <= heat_min | k >= heat_max){
    k = (heat_min + heat_max)/2
  }

  par(mar = c(6, 6, 4, 5))
  res <- plot(
    heatmap_matrix,
    col = c(heat.colors(n = N, alpha = 0.9), colfunc(N)),
    breaks = c(seq(heat_min, k, length.out = N), seq(k, heat_max, length.out = N)),
    na.col = "green",
    digits = 3,
    cex = 0.8,
    main = mainname,
    xlab = '',
    ylab = '',
    axis.col = NULL,
    axis.row = NULL,
    key=list(side=4,  
             font=2, 
             cex.axis=0.75), 
    fmt.key="%.2f", 
    polygon.key=NULL, 
    axis.key=NULL, 
    border= NA,
    spacing.key=c(1,0.5,0) 
  )
  
  for (i in 1:length(res$cell.polygon)) {
    text <- res$cell.text[[i]]
    args <- res$cell.polygon[[i]]
    args$border     <- 1
    do.call("polygon", args)
    do.call("text", text)
  }
  
  axis(
    side = 2,
    at = seq(1, ar+1, 1),
    labels = rev(xtick),
    las = 2,
    cex.axis = 0.9
  )
  axis(
    side = 1,
    at = seq(1, ma+1, 1),
    labels = ytick,
    las = 2,
    cex.axis = 0.9
  )
 }

# Plot residuals for chosen models

# RV
fit_rv <- fit
par(mfrow = c(3,1))
plot(ts(fit_rv$fit[[8]]@fit$z), main = "ARMA(2, 1)", ylab="Residuals")
plot(ts(fit_rv$fit[[32]]@fit$z), main = "ARFIMA(1,0.371,1)", ylab="Residuals")
plot(ts(fit_rv$fit[[40]]@fit$z), main = "ARFIMA(4,0.419,2)", ylab="Residuals")

# VIX
fit_vix <- fit
par(mfrow = c(1,1))
plot(ts(fit_vix$fit[[2]]@fit$z), main = "ARMA(1, 0)", ylab="Residuals")


# QQ-plots og Ljung-Box Test
# First choose the model 
model = 2
# Get the standardized residuals of the model
ressta <- fit$fit[[model]]@fit$z
par(mfrow = c(1,1))
# Make the QQ-plot
qqnorm(ressta, pch=19, xlab = "Theoretical quantiles", cex.main=2,
       ylab = "Quantiles for the standardized residuals",main="ARFIMA(4, 0.419, 2)", 
       cex.lab=1.5);
qqline(ressta)

# Perform the Ljung-Box Test for the model for a chosen lag
Box.test(fit$fit[[model]]@fit$z, lag = 1, type = "Ljung")

#--------------------------------------------------#
#                  Forecast RV                     #
#--------------------------------------------------#
# Choose model to forecast
model = 32
T = length(returns)

# Initialize full dataset
data <- read_excel("data_rv_10year.xlsx")
data <- log(read.zoo(data))

# Make outsample
outsample <- data[index(data)>="2021-02-26"]

out_of_sample <- NROW(outsample)
nahead <- out_of_sample

# Rolling forecast
dates_out_of_sample <- index(outsample)
arma_spec = getspec(fit$fit[[model]])
arma_fit <- arfimafit(spec = arma_spec, data = data, out.sample = out_of_sample)

# Forecast log-returns along the whole out-of-sample
arma_fore <- arfimaforecast(arma_fit, n.roll = out_of_sample-1)
forecast_log_returns <- xts(arma_fore@forecast$seriesFor[1, ], dates_out_of_sample)

basket = cbind(xts(returns[index(returns)>"2020-02-26"]), 
               outsample,
               forecast_log_returns)
zoo.basket <- xts(as.zoo(basket))
plot(zoo.basket, screens = 1, col = c("black","palegreen4","orangered3"), main = "")
addLegend(legend.loc = "bottomleft", legend = c("Original", "Realized", "Forecasted"),
          lty = 1, lwd = 3, y.intersp=2,col = c("black","palegreen4","orangered3"))

# One-year forecast
nahead <- out_of_sample
arfimaforecast <- arfimaforecast(fit$fit[[model]], n.ahead = nahead)

# Apply dates to forecast
forecast_dates <- index(outsample)
combi <- data.frame(forecast_dates, arfimaforecast@forecast$seriesFor[,1])

# Generate confidence intervals
sigma <- tail(coef(fit$fit[[model]]), 1)
drift <- function(i){
  sigma_h <- sigma*sqrt(i*(1+i/T))
}

sigma_h <- sapply(1:nahead, drift)

conf95 <- cbind(c(arfimaforecast@forecast$seriesFor - 1.96*sigma_h), c(arfimaforecast@forecast$seriesFor + 1.96*sigma_h))
conf80 <- cbind(c(arfimaforecast@forecast$seriesFor - 1.28*sigma_h), c(arfimaforecast@forecast$seriesFor + 1.28*sigma_h))
combi2 <- data.frame(forecast_dates, conf95)
combi3 <- data.frame(forecast_dates, conf80)

# Plot realized values in conjunction with confidence intervals
basket = cbind(xts(returns[index(returns)>"2020-02-26"]), 
               xts(head(outsample, nahead)), 
               xts(read.zoo(combi)),
               xts(read.zoo(combi2)),
               xts(read.zoo(combi3)))
zoo.basket <- xts(as.zoo(basket))
plot.xts(zoo.basket, screens = 1, col = c("black", "palegreen4","orangered3", "skyblue4", "skyblue4", "skyblue3", "skyblue3"), main = "")
addLegend(legend.loc = "bottomleft", legend = c("Original", "Realized", "Forecasted", "Conf. int. 80","Conf. int. 95"),
          lty = 1, lwd = 3, y.intersp = 2, col = c("black","palegreen4","orangered3", "skyblue3", "skyblue4"))

# Color confidence intervals
shadearea <- function(dates, ts1, ts2, color){
  m <- matrix(c(ts1, ts2), ncol = 2)
  m.xts <- xts(m, order.by = dates)
  print(addPolygon(m.xts, on = 1, col = color))
}

shadearea(dates_out_of_sample,
          arfimaforecast@forecast$seriesFor + 1.96*sigma_h,
          arfimaforecast@forecast$seriesFor + 1.28*sigma_h,
          rgb(74/255,112/255,139/255,0.5))

shadearea(dates_out_of_sample,
          arfimaforecast@forecast$seriesFor + 1.28*sigma_h,
          arfimaforecast@forecast$seriesFor - 1.28*sigma_h,
          rgb(108/255,166/255,205/255,0.5))
shadearea(dates_out_of_sample,
          arfimaforecast@forecast$seriesFor - 1.28*sigma_h,
          arfimaforecast@forecast$seriesFor - 1.96*sigma_h,
          rgb(74/255,112/255,139/255,0.5))


