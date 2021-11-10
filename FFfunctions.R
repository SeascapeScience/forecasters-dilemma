
linearFF <- function(Y0 = 0.6, # Initial condition
                     beta = .5, # Parameter for linear autoregression model (slope)
                     alpha = 2, # Parameter for linear autoregression model (intercept)
                     epsilon = 0.5, # Error/uncertainty term
                     rho0 = 0, # Reaction effect (0 means no effect)
                     trusti = 0, # Trust effect (0 means no effect)
                     nt = 100, # Number of time steps
                     iters = 1, # Number of feedback iterations to compute forecast
                     memori = 1, # Number of time steps back the memory goes (for calculating rho and tau)
                     Zcast = NA, # Predetermined forecast vector (instead of forecasting dynamically)
                     maintitle = NA, # Main title
                     plotit = T) # Flag for turning on/off plot
{
  Y <- Y0
  if (is.na(Zcast[1])) {
    Z <- Y0
  } else {
    Z <- Zcast[1]
  }
  if (length(Zcast) == 1)
  {
    Zcast[1:nt] <- Zcast[1]
  }
  iters <- round(max(iters,1))
  for (i in 2:nt)
  {
    rho <- 0
    tau <- 0
    if (is.na(Zcast[1]))
    {
      for (j in 1:iters)
      {
        Z[i] <- beta*Y[i-1] + alpha + (runif(1)-.5)*epsilon - rho*tau
        zi <- sum(Z[max(1,i-memori+1):i])
        yim1 <- sum(Y[max(1,i-memori):(i-1)])
        zim1 <- sum(Z[max(1,i-memori):(i-1)])
        rho <- tanh(zi)*rho0
        tau <- exp(-trusti*abs(zim1-yim1)/yim1)
      } 
    } else {
      Z[i] <- Zcast[i]
      zi <- sum(Z[max(1,i-memori+1):i])
      yim1 <- sum(Y[max(1,i-memori):(i-1)])
      zim1 <- sum(Z[max(1,i-memori):(i-1)])
      rho <- tanh(zi)*rho0
      tau <- exp(-trusti*abs(zim1-yim1)/yim1)
    }
    Y[i] <- beta*Y[i-1] + alpha + (runif(1)-.5)*epsilon - rho*tau
  }
  if (plotit)
  {
    plot(Z,type='l',
         lty = 1, 
         col='blue',
         ylim = c(0,max(max(Z),max(Y))),
         xlab = 'Time step',
         ylab = 'Y',
         lwd = 3,
         main = maintitle)
    lines(Y, lwd = 2)
    legend(nt*.4,1.5,
           legend = c('actual','forecast'), 
           col = c('black','blue'), lty = 1:1, lwd = 2,
           cex = 1)
  }
  dataout <- list("Y" = Y, "Z" = Z)
  return(dataout)
}