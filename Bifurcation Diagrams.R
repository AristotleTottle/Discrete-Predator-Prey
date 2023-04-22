############################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Always run the three lines of code below #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
############################################

library(rgl)
a <- 1
b <- 1

#####################################
#   bifurcation diagram, k1 fixed   #
#| | | | | | | | | | | | | | | | | |#
#| | | | | | | | | | | | | | | | | |#
#V V V V V V V V V V V V V V V V V V#
#####################################

k1Fixed <- 2.5
k2Int <- 3.2
k2End <- 3.9
k2Seq <- seq(from = k2Int, to = k2End, by = (k2End-k2Int)/899)
xPoints <- data.frame(num <- 1:500)
yPoints <- data.frame(num <- 1:500)

for(t in 1:length(k2Seq))
{
  fx <- (k1Fixed*k2Seq[t]+a*k2Seq[t]-k2Seq[t]-a)/(k1Fixed*k2Seq[t]+a*b)
  fy <- 1-(1+b*fx)/k2Seq[t]
  x0 <- c(fx*runif(1, min = .8, max = .9))
  y0 <- c(fy*runif(1, min = .8, max = .9))
  # x0 <- c(fx*runif(1, min = .4, max = .5))
  # y0 <- c(fy*runif(1, min = .7, max = .8))
  for(z in 1:799)
  {
    x0 <- c(x0, k1Fixed*x0[z]*(1-x0[z])+a*x0[z]*y0[z])
    y0 <- c(y0, k2Seq[t]*y0[z]*(1-y0[z])-b*x0[z]*y0[z])
  }
  xPoints[paste(t)] <- x0[301:800]
  yPoints[paste(t)] <- y0[301:800]
}
xPoints <- xPoints[-c(1)]
yPoints <- yPoints[-c(1)]

str <- paste("k1 = ", k1Fixed, sep = "")
plot(x = NULL, xlim = c(k2Seq[1], k2Seq[length(k2Seq)]), ylim = c(0,1.1), xlab = "k2", ylab = "x", main = str)
for(t in 1:length(k2Seq))
{
  k2Temp <- seq(k2Seq[t], k2Seq[t], length = 500)
  points(x = k2Temp, y = xPoints[,t], cex = .001)
}

#####################################
#   bifurcation diagram, k2 fixed   #
#| | | | | | | | | | | | | | | | | |#
#| | | | | | | | | | | | | | | | | |#
#V V V V V V V V V V V V V V V V V V#
#####################################

k2Fixed <- 3.5
k1Int <- 1.8
k1End <- 3.45
k1Seq <- seq(from = k1Int, to = k1End, by = (k1End-k1Int)/599)
xPoints <- data.frame(num <- 1:500)
yPoints <- data.frame(num <- 1:500)

for(t in 1:length(k1Seq))
{
  fx <- (k1Seq[t]*k2Fixed+a*k2Fixed-k2Fixed-a)/(k1Seq[t]*k2Fixed+a*b)
  fy <- 1-(1+b*fx)/k2Fixed
  # x0 <- c(.4)
  # y0 <- c(.4)
  x0 <- c(fx*runif(1, min = .96, max = .97))
  y0 <- c(fy*runif(1, min = .96, max = .97))
  # x0 <- c(fx*runif(1, min = .2, max = .4))
  # y0 <- c(fy*runif(1, min = .2, max = .4))
  for(z in 1:599)
  {
    x0 <- c(x0, k1Seq[t]*x0[z]*(1-x0[z])+a*x0[z]*y0[z])
    y0 <- c(y0, k2Fixed*y0[z]*(1-y0[z])-b*x0[z]*y0[z])
  }
  xPoints[paste(t)] <- x0[301:800]
  yPoints[paste(t)] <- y0[301:800]
}

xPoints <- xPoints[-c(1)]
yPoints <- yPoints[-c(1)]

str <- paste("k2 = ", k2Fixed, sep = "")
plot(x = NULL, xlim = c(k1Seq[1], k1Seq[length(k1Seq)]), ylim = c(0,1.1), xlab = "k1", ylab = "y", main = str)
for(t in 1:length(k1Seq))
{
  k1Temp <- seq(k1Seq[t], k1Seq[t], length = 500)
  points(x = k1Temp, y = yPoints[,t], cex = .001)
}

######################################
#  3D bifurcation diagram, k1 fixed  #
#| | | | | | | | |  | | | | | | | | |#
#| | | | | | | | |  | | | | | | | | |#
#V V V V V V V V V  V V V V V V V V V#
######################################

k1Fixed <- 2.5
k2Int <- 3.2
k2End <- 3.8
k2Seq <- seq(from = k2Int, to = k2End, by = (k2End-k2Int)/1599)

xPoints <- data.frame(num <- 1:300)
yPoints <- data.frame(num <- 1:300)
for(t in 1:length(k2Seq))
{
  fx <- (k1Fixed*k2Seq[t]+a*k2Seq[t]-k2Seq[t]-a)/(k1Fixed*k2Seq[t]+a*b)
  fy <- 1-(1+b*fx)/k2Seq[t]
  # x0 <- c(fx*runif(1, min = .4, max = .5))
  # y0 <- c(fy*runif(1, min = .7, max = .8))
  x0 <- c(fx*runif(1, min = .98, max = .99))
  y0 <- c(fy*runif(1, min = .98, max = .99))
  for(z in 1:599)
  {
    x0 <- c(x0, k1Fixed*x0[z]*(1-x0[z])+a*x0[z]*y0[z])
    y0 <- c(y0, k2Seq[t]*y0[z]*(1-y0[z])-b*x0[z]*y0[z])
  }
  xPoints[paste(t)] <- x0[301:600]
  yPoints[paste(t)] <- y0[301:600]
}
xPoints <- xPoints[-c(1)]
yPoints <- yPoints[-c(1)]
plot3d(x = NULL, xlab = "x", ylab = "k_2", zlab = "y", scale = .5)
for(t in 1:length(k2Seq))
{
  k2Temp <- seq(k2Seq[t], k2Seq[t], length = length(xPoints[,t]))
  points3d(x = xPoints[,t], y = k2Temp, z = yPoints[,t])
}

######################################
#  3D bifurcation diagram, k2 fixed  #
#| | | | | | | | |  | | | | | | | | |#
#| | | | | | | | |  | | | | | | | | |#
#V V V V V V V V V  V V V V V V V V V#
######################################

k2Fixed <- 3.525
k1Int <- 2
k1End <- 3.39
k1Seq <- seq(from = k1Int, to = k1End, by = (k1End-k1Int)/1599)

xPoints <- data.frame(num <- 1:300)
yPoints <- data.frame(num <- 1:300)
for(t in 1:length(k1Seq))
{
  fx <- (k1Seq[t]*k2Fixed+a*k2Fixed-k2Fixed-a)/(k1Seq[t]*k2Fixed+a*b)
  fy <- 1-(1+b*fx)/k2Fixed
  x0 <- c(fx*runif(1, min = .96, max = .97))
  y0 <- c(fy*runif(1, min = .96, max = .97))
  # x0 <- c(fx*runif(1, min = .4, max = .6))
  # y0 <- c(fy*runif(1, min = .85, max = .9))
  for(z in 1:599)
  {
    x0 <- c(x0, k1Seq[t]*x0[z]*(1-x0[z])+a*x0[z]*y0[z])
    y0 <- c(y0, k2Fixed*y0[z]*(1-y0[z])-b*x0[z]*y0[z])
  }
  xPoints[paste(t)] <- x0[301:600]
  yPoints[paste(t)] <- y0[301:600]
}
xPoints <- xPoints[-c(1)]
yPoints <- yPoints[-c(1)]
plot3d(x = NULL, xlab = "x", ylab = "k_1", zlab = "y", scale = .5)
for(t in 1:length(k1Seq))
{
  k1Temp <- seq(k1Seq[t], k1Seq[t], length = length(xPoints[,t]))
  points3d(x = xPoints[,t], y = k1Temp, z = yPoints[,t])
}