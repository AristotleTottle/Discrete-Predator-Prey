library(rgl)
k1 <- 2.133306
k2 <- 3.5
a <- 1
b <- 1

### equilibrium points at origin and along x/y axes

origin <- c(0,0)
yAxis <- c(0, (k2-1)/k2)
xAxis <- c((k1-1)/k1, 0)

### fourth equilibrium point with positive reals

fourthx <- (k1*k2+a*k2-k2-a)/(k1*k2+a*b)
fourthy <- 1-(1+b*fourthx)/k2

fourth <- c(fourthx, fourthy)

equilx <- c(origin[1], yAxis[1], xAxis[1], fourth[1])
equily <- c(origin[2], yAxis[2], xAxis[2], fourth[2])

### Jacobian matrices & eigenvalues

originM <- matrix(data = c(k1-2*k1*origin[1]+a*origin[2], a*origin[1], -b*origin[2], k2 - 2*k2*origin[2]-b*origin[1]), byrow = TRUE, nrow = 2, ncol = 2)
originEigen <- c(eigen(originM)$values[1], eigen(originM)$values[2])

xAxisM <- matrix(data = c(k1-2*k1*xAxis[1]+a*xAxis[2], a*xAxis[1], -b*xAxis[2], k2 - 2*k2*xAxis[2]-b*xAxis[1]), byrow = TRUE, nrow = 2, ncol = 2)
xAxisEigen <- c(eigen(xAxisM)$values[1], eigen(xAxisM)$values[2])

yAxisM <- matrix(data = c(k1-2*k1*yAxis[1]+a*yAxis[2], a*yAxis[1], -b*yAxis[2], k2 - 2*k2*yAxis[2]-b*yAxis[1]), byrow = TRUE, nrow = 2, ncol = 2)
yAxisEigen <- c(eigen(yAxisM)$values[1], eigen(yAxisM)$values[2])

fourthM <- matrix(data = c(k1-2*k1*fourthx+a*fourthy, a*fourthx, -b*fourthy, 
                           k2 - 2*k2*fourthy-b*fourthx), byrow = TRUE, nrow = 2, ncol = 2)
fourthEigen <- c(eigen(fourthM)$values[1], eigen(fourthM)$values[2])

### initial conditions


x <- c(.4)
y <- c(.4)
# x <- c(fourthx*.9)
# y <- c(fourthy*.9)
# x <- c(fourthx*.55)
# y <- c(fourthy*.75)
# x <- c(0.4)
# y <- c(0.4)

timeStop <- 1000

for(t in 1:timeStop)
{
  x <- c(x, k1*x[t]*(1-x[t])+a*x[t]*y[t])
  y <- c(y, k2*y[t]*(1-y[t])-b*x[t]*y[t])
}

### plot x and y together

mainTitle <- paste("k1 = ", k1, ", k2 = ", k2, sep = "")
plot(x = x[300:length(x)], y = y[300:length(y)], xlim = c(0, 1.05), ylim = c(0, 1.05), 
     main = mainTitle, xlab = "x_n", ylab = "y_n", cex = .4, pch = 19)
### plot(x = x[0:401], y = y[0:401], xlim = c(min(equilx)*1.5, max(equilx)*1.5), ylim = c(min(equily)*1.5, max(equily)*1.5))

# plot equilibrium points as red
points(x = origin[1], y = origin[2], col = "red", pch = 19, cex = 1)
points(x = xAxis[1], y = xAxis[2], col = "red", pch = 19, cex = 1)
points(x = yAxis[1], y = yAxis[2], col = "red", pch = 19, cex = 1)
points(x = fourth[1], y = fourth[2], col = "red", pch = 19, cex = 1)

for(t in 1:timeStop)
{
  arrows(x0 = x[t], y0 = y[t], x1 = x[t+1], y1 = y[t+1], code = 2, length = 0.07)
}

### plot n and x_n

plot(x = 0, y = x[1], xlim = c(0, timeStop), ylim = c(min(x), max(x)), cex = 0.5, xlab = "time", ylab = "x_n")

for(t in 1:timeStop)
{
  points(x = t, y = x[t+1], cex = 0.5)
}

for(t in 1:timeStop)
{
  arrows(x0 = t-1, y0 = x[t], x1 = t, y1 = x[t+1], code = 2, length = 0.07)
}

### plot n and y_n

plot(x = 0, y = y[1], xlim = c(0, timeStop), ylim = c(min(y), max(y)), cex = 0.5)

for(t in 1:timeStop)
{
  points(x = t, y = y[t+1])
}

for(t in 1:timeStop)
{
  arrows(x0 = t-1, y0 = y[t], x1 = t, y1 = y[t+1], code = 2, length = 0.07)
}

### Bar chart (x)

coordinates <- data.frame(x = x[timeStop-100:timeStop], y = y[timeStop-100:timeStop])
sort.x <- coordinates[order(x[timeStop-100:timeStop]), ]

intx1 <- min(x[timeStop-100:timeStop])*.99
intx2 <- max(x[timeStop-100:timeStop])/.99


intIntervals <- 20 ## number of sub-intervals


intCounts <- c(1:intIntervals)*0 ## counts number of x-values falling into each sub-interval
intC <- 1 ## says which interval we're on for for-loop

intervalStep <- (intx2 - intx1) / intIntervals
interval <- c(intx1, intx1 + intervalStep) ## 1st interval, turns into 2nd interval, etc.

intervals <- matrix(data = interval, nrow = intIntervals, ncol = 2, byrow = TRUE)

for(r in 2:intIntervals)  ## start for getting names of intervals
{
  intervals[r, 1] <- interval[1] + (r-1)*intervalStep
  intervals[r, 2] <- interval[1] + r*intervalStep
}

intervalNames <- c(paste(intervals[1,1], "-", intervals[1,2]))

for(r in 2:intIntervals)
{
  intervalNames <- c(intervalNames, paste(intervals[r,1], "-", intervals[r,2]))
}  ## end for getting names of intervals

for(t in 0:100)
{
  if(interval[1] <= sort.x[t+1, 1] & sort.x[t+1, 1] < interval[2]) {
    intCounts[intC] <- intCounts[intC] + 1
  } else {
    intC <- intC + 1
    interval[1] <- interval[1] + intervalStep
    interval[2] <- interval[2] + intervalStep
    intCounts[intC] <- intCounts[intC] + 1
  }
}

barplot(intCounts, xlab = "Sub Intervals", ylab = "Occurences", names.arg = intervalNames, main = "Occurences for x-values in sub intervals")

### Bar Chart (y)

coordinates <- data.frame(x = x, y = y)
sort.y <- coordinates[order(y), ]

inty1 <- min(y)*.99
inty2 <- max(y)/.99


intIntervals <- 75 ## number of sub-intervals


intCounts <- c(1:intIntervals)*0 ## counts number of x-values falling into each sub-interval
intC <- 1 ## says which interval we're on for for-loop

intervalStep <- (inty2 - inty1) / intIntervals
interval <- c(inty1, inty1 + intervalStep) ## 1st interval, turns into 2nd interval, etc.

intervals <- matrix(data = interval, nrow = intIntervals, ncol = 2, byrow = TRUE)

for(r in 2:intIntervals)  ## start for getting names of intervals
{
  intervals[r, 1] <- interval[1] + (r-1)*intervalStep
  intervals[r, 2] <- interval[1] + r*intervalStep
}

intervalNames <- c(paste(intervals[1,1], "-", intervals[1,2]))

for(r in 2:intIntervals)
{
  intervalNames <- c(intervalNames, paste(intervals[r,1], "-", intervals[r,2]))
}  ## end for getting names of intervals

for(t in 0:timeStop)
{
  if(interval[1] <= sort.y[t+1, 2] & sort.y[t+1, 2] < interval[2]) {
    intCounts[intC] <- intCounts[intC] + 1
  } else {
    intC <- intC + 1
    interval[1] <- interval[1] + intervalStep
    interval[2] <- interval[2] + intervalStep
    intCounts[intC] <- intCounts[intC] + 1
  }
}

barplot(intCounts, xlab = "Sub Intervals", ylab = "Occurences", names.arg = intervalNames, main = "Occurences for y-values in sub intervals")

### Lyapunov exponent, one trajectory

eps <- 0.0001

distances <- data.frame(num <- 1:timeStop)
n <- 1  ### counter
k <- 60 ### x_j only goes to x_(timeStop-k)

for(i in 1:timeStop)
{
  for(j in (i + 1):(timeStop + 1))
  {
    d <- abs(sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2))
    if(d < eps && (timeStop + 1 - j) >= k)
    {
      dist <- c(d)
      for(t in 1:(timeStop + 1 - j))
      {
        dist <- c(dist, abs(sqrt((x[i+t] - x[j+t])^2 + (y[i+t] - y[j+t])^2)))
      }
      if((j-2) > 0)
      {
        for(v in 1:(j-2))
        {
          dist <- c(dist, 0)
        }
      }
      distances[paste(n)] <- dist
      n <- n + 1
    }
  }
}

distances <- distances[-c(1)]
lnDist <- log(distances)
lnMeanDist <- rowMeans(lnDist)

for(t in 1:length(lnMeanDist))
{
  if(lnMeanDist[t] == -Inf || lnMeanDist[t] == Inf || is.nan(lnMeanDist[t]))
  {
    break;
  }
}
kAdj <- t - 2
plot(x = 0:kAdj, y = lnMeanDist[1:(kAdj+1)], xlab = "time k", ylab = "log of average distance")

xV <- seq(from = 0, to = 60, by = 4)
xV <- c(0:60)
lmFit <- lm(formula = lnMeanDist[(xV+1)] ~ xV)
abline(lmFit)
lyapunov <- lmFit$coefficients[2]
lyapunov

### Lyapunov exponent, multiple trajectories

centerPair <- c(x[1], y[1])
eps <- 0.0000001 ## radius of circle containing points; for all distances between two points, |d_1 - d_2| <= 2*eps
numPairs <- 10
angle <- 0
xPairs <- c()
yPairs <- c()

for(t in 1:numPairs)
{
  xPairs <- c(xPairs, centerPair[1] + runif(1)*eps*cos(angle))
  yPairs <- c(yPairs, centerPair[2] + runif(1)*eps*sin(angle))
  angle <- angle + 2*pi/numPairs
}
xPairs <- c(xPairs, centerPair[1])
yPairs <- c(yPairs, centerPair[2])

pairs <- data.frame(x = xPairs, y = yPairs)

k <- 300 ### stop at t=k

distances <- data.frame(num <- 1:(k+1))
n <- 1

for(i in 1:(numPairs + 1))
{
  xM <- c(pairs[i, 1])
  yM <- c(pairs[i, 2])
  for(t in 1:k)
  {
    xM <- c(xM, k1*xM[t]*(1-xM[t])+a*xM[t]*yM[t])
    yM <- c(yM, k2*yM[t]*(1-yM[t])-b*xM[t]*yM[t])
    # xM <- c(xM, 3*xM[t])
    # yM <- c(yM, 2*yM[t])
  }
  if(i+1 <= numPairs + 1)
  {
    for(j in (i+1):(numPairs + 1))
    {
      xP <- c(pairs[j, 1])
      yP <- c(pairs[j, 2])
      dist <- c(sqrt((xM[1]-xP[1])^2+(yM[1]-yP[1])^2))
      for(t in 1:k)
      {
        xP <- c(xP, k1*xP[t]*(1-xP[t])+a*xP[t]*yP[t])
        yP <- c(yP, k2*yP[t]*(1-yP[t])-b*xP[t]*yP[t])
        # xP <- c(xP, 3*xP[t])
        # yP <- c(yP, 2*yP[t])
        dist <- c(dist, sqrt((xM[t+1]-xP[t+1])^2+(yM[t+1]-yP[t+1])^2))
      }
      distances[paste(n)] <- dist
      n <- n + 1
    }
  }
}

distances <- distances[-c(1)]
distancesDiv <- data.frame(num <- 1:k)
n <- 1

for(t in 1:length(distances[1,]))
{
  distDiv <- c()
  for(z in 2:length(distances[,1]))
  {
    distDiv <- c(distDiv, distances[z,t]/distances[z-1,t])
  }
  distancesDiv[paste(n)] <- distDiv
  n <- n + 1
}
distancesDiv <- distancesDiv[-c(1)]
lnMeanDist <- log(rowMeans(distancesDiv))

kAdj <- 1
for(t in 1:length(lnMeanDist))
{
  if(lnMeanDist[t] == -Inf || lnMeanDist[t] == Inf || is.nan(lnMeanDist[t]))
  {
    break;
  }
  kAdj <- kAdj + 1
}

plot(x = 0:(kAdj-1), y = lnMeanDist[1:kAdj], xlab = "time k", ylab = "lyapunov")

s <- floor(kAdj/2)
cycleStart <- lnMeanDist[s]
nCycle <- 1

for(t in (s+1):kAdj)
{
  if(abs(cycleStart - lnMeanDist[t]) < 0.001)
  {
    break;
  }
  nCycle <- nCycle + 1
}

lyapunov <- mean(lnMeanDist[s:(s + nCycle - 1)])
lyapunov

### bifurcation diagram, k1 fixed

k1Fixed <- 2.5
k2Int <- 3.2
k2End <- 3.8
k2Seq <- seq(from = k2Int, to = k2End, by = (k2End-k2Int)/899)
xPoints <- data.frame(num <- 1:500)
yPoints <- data.frame(num <- 1:500)

for(t in 1:length(k2Seq))
{
  fx <- (k1Fixed*k2Seq[t]+a*k2Seq[t]-k2Seq[t]-a)/(k1Fixed*k2Seq[t]+a*b)
  fy <- 1-(1+b*fx)/k2Seq[t]
  # x0 <- c(fx*runif(1, min = .8, max = .9))
  # y0 <- c(fy*runif(1, min = .8, max = .9))
  x0 <- c(fx*runif(1, min = .4, max = .5))
  y0 <- c(fy*runif(1, min = .7, max = .8))
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


### bifurcation diagram, k2 fixed

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
  x0 <- c(.4)
  y0 <- c(.4)
  # x0 <- c(fx*runif(1, min = .96, max = .97))
  # y0 <- c(fy*runif(1, min = .96, max = .97))
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


### 3D bifurcation diagram, across a rectangle

k1Int <- 3.35
k1End <- 3.50
k2Int <- 2.8
k2End <- 3
k1Seq <- seq(from = k1Int, to = k1End, by = (k1End-k1Int)/399)
k2Seq <- seq(from = k2Int, to = k2End, by = (k2End-k2Int)/399)
plot3d(x = NULL, xlab = "k_1", ylab = "k_2", zlab = "x")

for(t1 in 1:length(k1Seq))
{
  xPoints <- data.frame(num <- 1:300)
  yPoints <- data.frame(num <- 1:300)
  
  for(t2 in 1:length(k2Seq))
  {
    x0 <- c(runif(1, min = .5, max = .7))
    y0 <- c(runif(1, min = .5, max = .7))
    for(z in 1:599)
    {
      x0 <- c(x0, k1Seq[t1]*x0[z]*(1-x0[z])+a*x0[z]*y0[z])
      y0 <- c(y0, k2Seq[t2]*y0[z]*(1-y0[z])-b*x0[z]*y0[z])
    }
    xPoints[paste(t2)] <- x0[301:600]
    yPoints[paste(t2)] <- y0[301:600]
  }
  xPoints <- xPoints[-c(1)]
  yPoints <- yPoints[-c(1)]
  for(t in 1:length(k2Seq))
  {
    k1Temp <- seq(k1Seq[t1], k1Seq[t1], length = length(xPoints[,t]))
    k2Temp <- seq(k2Seq[t], k2Seq[t], length = length(xPoints[,t]))
    points3d(x = k1Temp, y = k2Temp, z = xPoints[,t])
  }
}

### 3D bifurcation diagram, single line
k1Int <- 1.5
k1End <- 3.4
k1Seq <- seq(from = k1Int, to = k1End, by = (k1End-k1Int)/899)

singLine <- function(k1){
  # -.1*k1+3.75
  # -.7*k1+5
  -0.7*k1+5.3
}
k2Seq <- singLine(k1Seq)

xPoints <- data.frame(num <- 1:300)
yPoints <- data.frame(num <- 1:300)
for(t in 1:length(k2Seq))
{
  fx <- (k1Seq[t]*k2Seq[t]+a*k2Seq[t]-k2Seq[t]-a)/(k1Seq[t]*k2Seq[t]+a*b)
  fy <- 1-(1+b*fx)/k2Seq[t]
  x0 <- c(fx*runif(1, min = .9, max = 1))
  y0 <- c(fy*runif(1, min = .9, max = 1))
  for(z in 1:599)
  {
    x0 <- c(x0, k1Seq[t]*x0[z]*(1-x0[z])+a*x0[z]*y0[z])
    y0 <- c(y0, k2Seq[t]*y0[z]*(1-y0[z])-b*x0[z]*y0[z])
  }
  xPoints[paste(t)] <- x0[301:600]
  yPoints[paste(t)] <- y0[301:600]
}
xPoints <- xPoints[-c(1)]
yPoints <- yPoints[-c(1)]
plot3d(x = NULL, xlab = "k_1", ylab = "k_2", zlab = "x", scale = .5)
for(t in 1:length(k2Seq))
{
  k1Temp <- seq(k1Seq[t], k1Seq[t], length = length(xPoints[,t]))
  k2Temp <- seq(k2Seq[t], k2Seq[t], length = length(xPoints[,t]))
  points3d(x = k1Temp, y = k2Temp, z = xPoints[,t])
}

### 3D bifurcation diagram, single line, only k1 is plotted
k1Int <- 1.5
k1End <- 3.4
k1Seq <- seq(from = k1Int, to = k1End, by = (k1End-k1Int)/899)

singLine <- function(k1){
  # -.1*k1+3.75
  # -.7*k1+5
  -0.7*k1+5.3
}
k2Seq <- singLine(k1Seq)

xPoints <- data.frame(num <- 1:300)
yPoints <- data.frame(num <- 1:300)
for(t in 1:length(k2Seq))
{
  fx <- (k1Seq[t]*k2Seq[t]+a*k2Seq[t]-k2Seq[t]-a)/(k1Seq[t]*k2Seq[t]+a*b)
  fy <- 1-(1+b*fx)/k2Seq[t]
  x0 <- c(fx*runif(1, min = .3, max = 2))
  y0 <- c(fy*runif(1, min = .3, max = 2))
  for(z in 1:599)
  {
    x0 <- c(x0, k1Seq[t]*x0[z]*(1-x0[z])+a*x0[z]*y0[z])
    y0 <- c(y0, k2Seq[t]*y0[z]*(1-y0[z])-b*x0[z]*y0[z])
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


### 3D bifurcation diagram, k1 fixed, x/y

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


### 3D bifurcation diagram, k2 fixed, x/y

k2Fixed <- 3.5
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

### lyapunov estimates over interval, k1 fixed

k1Fix <- 2.5
k2Int <- 3.2
k2End <- 3.6
k2Seq <- seq(from = k2Int, to = k2End, by = (k2End-k2Int)/499)
lyapunovPoints <- c()

for(z in 1:length(k2Seq))
{
  fthX <- (k1Fix*k2Seq[z]+a*k2Seq[z]-k2Seq[z]-a)/(k1Fix*k2Seq[z]+a*b)
  fthY <- 1-(1+b*fourthx)/k2Seq[z]
  # x0 <- c(fx*runif(1, min = .8, max = .9))
  # y0 <- c(fy*runif(1, min = .8, max = .9))
  # x0 <- c(fx*runif(1, min = .4, max = .5))
  # y0 <- c(fy*runif(1, min = .7, max = .8))
  # centerPair <- c(x[1], y[1])
  centerPair <- c(runif(1, min = 0.4, max = 0.5)*fthX, runif(1, min = 0.7, max = 0.8)*fthY)
  eps <- 0.0000001 ## radius of circle containing points; for all distances between two points, |d_1 - d_2| <= 2*eps
  numPairs <- 10
  angle <- 0
  xPairs <- c()
  yPairs <- c()
  
  for(t in 1:numPairs)
  {
    xPairs <- c(xPairs, centerPair[1] + runif(1)*eps*cos(angle))
    yPairs <- c(yPairs, centerPair[2] + runif(1)*eps*sin(angle))
    angle <- angle + 2*pi/numPairs
  }
  xPairs <- c(xPairs, centerPair[1])
  yPairs <- c(yPairs, centerPair[2])
  
  pairs <- data.frame(x = xPairs, y = yPairs)
  
  k <- 300 ### stop at t=k
  
  distances <- data.frame(num <- 1:(k+1))
  n <- 1
  
  for(i in 1:(numPairs + 1))
  {
    xM <- c(pairs[i, 1])
    yM <- c(pairs[i, 2])
    for(t in 1:k)
    {
      xM <- c(xM, k1Fix*xM[t]*(1-xM[t])+a*xM[t]*yM[t])
      yM <- c(yM, k2Seq[z]*yM[t]*(1-yM[t])-b*xM[t]*yM[t])
    }
    if(i+1 <= numPairs + 1)
    {
      for(j in (i+1):(numPairs + 1))
      {
        xP <- c(pairs[j, 1])
        yP <- c(pairs[j, 2])
        dist <- c(sqrt((xM[1]-xP[1])^2+(yM[1]-yP[1])^2))
        for(t in 1:k)
        {
          xP <- c(xP, k1Fix*xP[t]*(1-xP[t])+a*xP[t]*yP[t])
          yP <- c(yP, k2Seq[z]*yP[t]*(1-yP[t])-b*xP[t]*yP[t])
          dist <- c(dist, sqrt((xM[t+1]-xP[t+1])^2+(yM[t+1]-yP[t+1])^2))
        }
        distances[paste(n)] <- dist
        n <- n + 1
      }
    }
  }
  
  distances <- distances[-c(1)]
  distancesDiv <- data.frame(num <- 1:k)
  n <- 1
  
  for(t in 1:length(distances[1,]))
  {
    distDiv <- c()
    for(z in 2:length(distances[,1]))
    {
      distDiv <- c(distDiv, distances[z,t]/distances[z-1,t])
    }
    distancesDiv[paste(n)] <- distDiv
    n <- n + 1
  }
  distancesDiv <- distancesDiv[-c(1)]
  lnMeanDist <- log(rowMeans(distancesDiv))
  
  kAdj <- 1
  for(t in 1:length(lnMeanDist))
  {
    if(lnMeanDist[t] == -Inf || lnMeanDist[t] == Inf || is.nan(lnMeanDist[t]))
    {
      break;
    }
  }
  kAdj <- t - 1
  
  s <- floor(kAdj/2)
  cycleStart <- lnMeanDist[s]
  nCycle <- 1
  
  for(t in (s+1):kAdj)
  {
    if(abs(cycleStart - lnMeanDist[t]) < 0.001)
    {
      break;
    }
    nCycle <- nCycle + 1
  }
  
  lyapunovPoints <- c(lyapunovPoints, mean(lnMeanDist[s:(s + nCycle - 1)]))
}

str <- paste("k_1 = ", k1Fix, sep = "")
plot(x = k2Seq, y = lyapunovPoints, xlab = "k2", ylab = "Lyapunov estimates", 
     main = str, cex = 0.3, ylim = c(-1, 1))
abline(h = 0)

### lyapunov estimates over interval, k2 fixed

k2Fix <- 3.5
k1Int <- 1.8
k1End <- 3.45
k1Seq <- seq(from = k1Int, to = k1End, by = (k1End-k1Int)/699)
lyapunovPoints <- c()

for(z in 1:length(k1Seq))
{
  fthX <- (k1Fix*k2Seq[z]+a*k2Seq[z]-k2Seq[z]-a)/(k1Fix*k2Seq[z]+a*b)
  fthY <- 1-(1+b*fourthx)/k2Seq[z]
  # x0 <- c(fx*runif(1, min = .8, max = .9))
  # y0 <- c(fy*runif(1, min = .8, max = .9))
  # x0 <- c(fx*runif(1, min = .2, max = .4))
  # y0 <- c(fy*runif(1, min = .2, max = .4))
  # centerPair <- c(x[1], y[1])
  centerPair <- c(runif(1, min = 0.2, max = 0.4)*fthX, runif(1, min = 0.2, max = 0.4)*fthY)
  eps <- 0.0000001 ## radius of circle containing points; for all distances between two points, |d_1 - d_2| <= 2*eps
  numPairs <- 10
  angle <- 0
  xPairs <- c()
  yPairs <- c()
  
  for(t in 1:numPairs)
  {
    xPairs <- c(xPairs, centerPair[1] + runif(1)*eps*cos(angle))
    yPairs <- c(yPairs, centerPair[2] + runif(1)*eps*sin(angle))
    angle <- angle + 2*pi/numPairs
  }
  xPairs <- c(xPairs, centerPair[1])
  yPairs <- c(yPairs, centerPair[2])
  
  pairs <- data.frame(x = xPairs, y = yPairs)
  
  k <- 300 ### stop at t=k
  
  distances <- data.frame(num <- 1:(k+1))
  n <- 1
  
  for(i in 1:(numPairs + 1))
  {
    xM <- c(pairs[i, 1])
    yM <- c(pairs[i, 2])
    for(t in 1:k)
    {
      xM <- c(xM, k1Seq[z]*xM[t]*(1-xM[t])+a*xM[t]*yM[t])
      yM <- c(yM, k2Fix*yM[t]*(1-yM[t])-b*xM[t]*yM[t])
    }
    if(i+1 <= numPairs + 1)
    {
      for(j in (i+1):(numPairs + 1))
      {
        xP <- c(pairs[j, 1])
        yP <- c(pairs[j, 2])
        dist <- c(sqrt((xM[1]-xP[1])^2+(yM[1]-yP[1])^2))
        for(t in 1:k)
        {
          xP <- c(xP, k1Seq[z]*xP[t]*(1-xP[t])+a*xP[t]*yP[t])
          yP <- c(yP, k2Fix*yP[t]*(1-yP[t])-b*xP[t]*yP[t])
          dist <- c(dist, sqrt((xM[t+1]-xP[t+1])^2+(yM[t+1]-yP[t+1])^2))
        }
        distances[paste(n)] <- dist
        n <- n + 1
      }
    }
  }
  
  distances <- distances[-c(1)]
  distancesDiv <- data.frame(num <- 1:k)
  n <- 1
  
  for(t in 1:length(distances[1,]))
  {
    distDiv <- c()
    for(z in 2:length(distances[,1]))
    {
      distDiv <- c(distDiv, distances[z,t]/distances[z-1,t])
    }
    distancesDiv[paste(n)] <- distDiv
    n <- n + 1
  }
  distancesDiv <- distancesDiv[-c(1)]
  lnMeanDist <- log(rowMeans(distancesDiv))
  
  for(t in 1:length(lnMeanDist))
  {
    if(lnMeanDist[t] == -Inf || lnMeanDist[t] == Inf || is.nan(lnMeanDist[t]))
    {
      break;
    }
  }
  kAdj <- t - 1
  
  s <- floor(kAdj/2)
  cycleStart <- lnMeanDist[s]
  nCycle <- 1
  
  for(v in (s+1):kAdj)
  {
    diff <- abs(cycleStart - lnMeanDist[v])
    if(diff < 0.001)
    {
      break;
    }
    else
    {
      nCycle <- nCycle + 1
    }
  }
  
  lyapunovPoints <- c(lyapunovPoints, mean(lnMeanDist[s:(s + nCycle - 1)]))
}

str <- paste("k2 = ", k2Fix, sep = "")
plot(x = k1Seq, y = lyapunovPoints, xlab = "k1", ylab = "Lyapunov estimates", 
     main = str, cex = 0.3, ylim = c(-1,1))
abline(h = 0)

### real part over interval, k1 fixed

k1Fixed <- 2.5
k2Int <- 3.2
k2End <- 3.6
k2Seq <- seq(from = k2Int, to = k2End, by = (k2End-k2Int)/899)
plot(x = NULL, xlab = "k2", ylab = "Reals", xlim = c(-2,2), ylim = c(-2,2))
for(t in 1:length(k2Seq))
{
  fx <- (k1Fixed*k2Seq[t]+a*k2Seq[t]-k2Seq[t]-a)/(k1Fixed*k2Seq[t]+a*b)
  fy <- 1-(1+b*fx)/k2Seq[t]
  
  fM <- matrix(data = c(k1Fixed-2*k1Fixed*fx+a*fy, a*fx, -b*fy, 
                             k2Seq[t] - 2*k2Seq[t]*fy-b*fx), byrow = TRUE, nrow = 2, ncol = 2)
  if(Im(eigen(fM)$values[1]) == 0)
  {
    points(x = k2Seq[t], y = eigen(fM)$values[1])
    points(x = k2Seq[t], y = eigen(fM)$values[2])
  }
  else
  {
    points(x = k2Seq[t], y = Re(eigen(fM)$values[1]))
  }
}

plot(x = k2Seq, y = Reals, xlab = "k2", ylab = "Reals")

### real part over interval, k2 fixed

k2Fixed <- 3.4
k1Int <- 3.33
k1End <- 4
k1Seq <- seq(from = k1Int, to = k1End, by = (k1End-k1Int)/999)
plot(x = NULL, xlab = "Reals", ylab = "k1", ylim = c(k1Int,k1End), xlim = c(-3,-.5))
imReals <- c()
reReals <- data.matrix(num <- 1:2)
for(t in 1:length(k1Seq))
{
  fx <- (k1Seq[t]*k2Fixed+a*k2Fixed-k2Fixed-a)/(k1Seq[t]*k2Fixed+a*b)
  fy <- 1-(1+b*fx)/k2Fixed
  
  fM <- matrix(data = c(k1Seq[t]-2*k1Seq[t]*fx+a*fy, a*fx, -b*fy, 
                        k2Fixed - 2*k2Fixed*fy-b*fx), byrow = TRUE, nrow = 2, ncol = 2)
  if(Im(eigen(fM)$values[1]) == 0)
  {
    points(y = k1Seq[t], x = eigen(fM)$values[1], cex = 0.4)
    points(y = k1Seq[t], x = eigen(fM)$values[2], cex = 0.4)
    reReals <- cbind(reReals, c(eigen(fM)$values[1], eigen(fM)$values[2]))
  }
  else
  {
    points(y = k1Seq[t], x = Re(eigen(fM)$values[1]), cex = 0.4)
    imReals <- c(imReals, Re(eigen(fM)$values[1]))
  }
}
reReals <- reReals[,-c(1)]
data <- data.frame(k1Var = k1Seq[1:length(imReals)], iReal = imReals)
line <- lm(k1Var ~ iReal, data = data)
abline(line)

data2 <- data.frame(iReals = imReals, k1Vars = k1Seq[1:length(imReals)])
line2 <- lm(iReals ~ k1Vars, data = data2)
newdata2 <- data.frame(k1Vars = k1Seq[(length(imReals) + 1):length(k1Seq)])
predict2 <- predict(line2, newdata = newdata2)

k1Est <- 0
div <- 2
for(t in 1:(length(k1Seq)-length(imReals)))
{
  divTemp <- abs((predict2[t]-reReals[1,t])/(predict2[t])-reReals[2,t])
  if(divTemp < 1)
  {
    divTemp <- 1/divTemp
  }
  if(divTemp < div)
  {
    div <- divTemp
    k1Est <- k1Seq[t+length(imReals)]
  }
}

real1 <- reReals[2,]
R1 <- data.frame(res = reReals[1,])
qdata1 <- data.frame(reReals = real1, k1Vars = k1Seq[(length(imReals) + 1):length(k1Seq)])
quad1 <- lm(k1Vars ~ poly(reReals,2), data = qdata1)
predQuad1 <- predict(quad1, newdata = R1)
lines(x = real1, y = predQuad1, col = "red")

### imaginary part over interval, k1 fixed

k1Fixed <- 2.5
k2Int <- 3.2
k2End <- 3.6
k2Seq <- seq(from = k2Int, to = k2End, by = (k2End-k2Int)/899)
Imags <- c()

for(t in 1:length(k2Seq))
{
  fx <- (k1Fixed*k2Seq[t]+a*k2Seq[t]-k2Seq[t]-a)/(k1Fixed*k2Seq[t]+a*b)
  fy <- 1-(1+b*fx)/k2Seq[t]
  
  fM <- matrix(data = c(k1Fixed-2*k1Fixed*fx+a*fy, a*fx, -b*fy, 
                        k2Seq[t] - 2*k2Seq[t]*fy-b*fx), byrow = TRUE, nrow = 2, ncol = 2)
  Imags <- c(Imags, abs(Im(eigen(fM)$values[1])))
}

plot(x = k2Seq, y = Imags, xlab = "k2", ylab = "Imaginary")

### real and imaginary part over interval, k1 fixed

k1Fixed <- 2
k2Int <- 0.2
k2End <- 5
k2Seq <- seq(from = k2Int, to = k2End, by = (k2End-k2Int)/899)
Reals <- c()
Imags <- c()

for(t in 1:length(k2Seq))
{
  fx <- (k1Fixed*k2Seq[t]+a*k2Seq[t]-k2Seq[t]-a)/(k1Fixed*k2Seq[t]+a*b)
  fy <- 1-(1+b*fx)/k2Seq[t]
  
  fM <- matrix(data = c(k1Fixed-2*k1Fixed*fx+a*fy, a*fx, -b*fy, 
                        k2Seq[t] - 2*k2Seq[t]*fy-b*fx), byrow = TRUE, nrow = 2, ncol = 2)
  Reals <- c(Reals, Re(eigen(fM)$values[1]))
  Imags <- c(Imags, abs(Im(eigen(fM)$values[1])))
}
plot3d(x = Reals, y = k2Seq, z = Imags, xlab = "Real", 
       ylab = "k_2", zlab = "Imaginary", scale = .5)

### real and imaginary part over interval, k2 fixed

k2Fixed <- 3.4
k1Int <- 3.3
k1End <- 3.5
k1Seq <- seq(from = k1Int, to = k1End, by = (k1End-k1Int)/899)
plot3d(x = NULL, xlab = "Real", ylab = "k_1", zlab = "Imaginary", scale = .5)
for(t in 1:length(k1Seq))
{
  fx <- (k1Seq[t]*k2Fixed+a*k2Fixed-k2Fixed-a)/(k1Seq[t]*k2Fixed+a*b)
  fy <- 1-(1+b*fx)/k2Fixed
  
  fM <- matrix(data = c(k1Seq[t]-2*k1Seq[t]*fx+a*fy, a*fx, -b*fy, 
                        k2Fixed - 2*k2Fixed*fy-b*fx), byrow = TRUE, nrow = 2, ncol = 2)
  if(Im(eigen(fM)$values[1]) == 0)
  {
    points3d(x = eigen(fM)$values[1], y = k1Seq[t], z = 0)
    points3d(x = eigen(fM)$values[2], y = k1Seq[t], z = 0)
  }
  else
  {
    points3d(x = Re(eigen(fM)$values[1]), y = k1Seq[t], z = Im(eigen(fM)$values[1]))
    points3d(x = Re(eigen(fM)$values[2]), y = k1Seq[t], z = Im(eigen(fM)$values[2]))
  }
}

### imaginary over k1 and k2 variable

k1Int <- 0.1
k1End <- 6
k2Int <- 0.1
k2End <- 6
k1Seq <- seq(from = k1Int, to = k1End, by = (k1End-k1Int)/399)
k2Seq <- seq(from = k2Int, to = k2End, by = (k2End-k2Int)/399)
Imags <- data.frame(num <- 1:length(k2Seq))

plot3d(x = NULL, xlab = "k_1", ylab = "k_2", zlab = "Imaginary", scale = .5)
for(t1 in 1:length(k1Seq))
{
  ImagsTemp <- c()
  for(t2 in 1:length(k2Seq))
  {
    fx <- (k1Seq[t1]*k2Seq[t2]+a*k2Seq[t2]-k2Seq[t2]-a)/(k1Seq[t1]*k2Seq[t2]+a*b)
    fy <- 1-(1+b*fx)/k2Seq[t2]
    fM <- matrix(data = c(k1Seq[t1]-2*k1Seq[t1]*fx+a*fy, a*fx, -b*fy, 
                          k2Seq[t2] - 2*k2Seq[t2]*fy-b*fx), byrow = TRUE, nrow = 2, ncol = 2)
    ImEig <- abs(Im(eigen(fM)$values[1]))
    if(ImEig == 0)
    {
      ImagsTemp <- c(ImagsTemp, -1)
    }
    else
    {
      ImagsTemp <- c(ImagsTemp, ImEig)
    }
    
  }
  k1Temp <- seq(k1Seq[t1], k1Seq[t1], length = length(ImagsTemp))
  points3d(x = k1Temp, y = k2Seq, z = ImagsTemp)
}


### real over k1 and k2 variable

k1Int <- 0.1
k1End <- 6
k2Int <- 0.1
k2End <- 6
k1Seq <- seq(from = k1Int, to = k1End, by = (k1End-k1Int)/199)
k2Seq <- seq(from = k2Int, to = k2End, by = (k2End-k2Int)/199)
Imags <- data.frame(num <- 1:length(k2Seq))
plot3d(x = NULL, xlab = "k_1", ylab = "k_2", zlab = "Real", scale = .5)
for(t1 in 1:length(k1Seq))
{
  ImagsTemp <- c()
  for(t2 in 1:length(k2Seq))
  {
    fx <- (k1Seq[t1]*k2Seq[t2]+a*k2Seq[t2]-k2Seq[t2]-a)/(k1Seq[t1]*k2Seq[t2]+a*b)
    fy <- 1-(1+b*fx)/k2Seq[t2]
    fM <- matrix(data = c(k1Seq[t1]-2*k1Seq[t1]*fx+a*fy, a*fx, -b*fy, 
                          k2Seq[t2] - 2*k2Seq[t2]*fy-b*fx), byrow = TRUE, nrow = 2, ncol = 2)
    
    ImEig <- Re(eigen(fM)$values[1])
    ImagsTemp <- c(ImagsTemp, ImEig)
  }
  k1Temp <- seq(k1Seq[t1], k1Seq[t1], length = length(ImagsTemp))
  points3d(x = k1Temp, y = k2Seq, z = ImagsTemp)
}

### "R" over k1 and k2 variable

k1Int <- 3.25
k1End <- 3.5
k2Int <- 3
k2End <- 3.5
k1Seq <- seq(from = k1Int, to = k1End, by = (k1End-k1Int)/199)
k2Seq <- seq(from = k2Int, to = k2End, by = (k2End-k2Int)/199)
Rs <- data.frame(num <- 1:length(k2Seq))
plot3d(x = NULL, xlab = "k_1", ylab = "k_2", zlab = "Real", scale = .5)
for(t1 in 1:length(k1Seq))
{
  RsTemp <- c()
  for(t2 in 1:length(k2Seq))
  {
    fx <- (k1Seq[t1]*k2Seq[t2]+a*k2Seq[t2]-k2Seq[t2]-a)/(k1Seq[t1]*k2Seq[t2]+a*b)
    fy <- 1-(1+b*fx)/k2Seq[t2]
    fM <- matrix(data = c(k1Seq[t1]-2*k1Seq[t1]*fx+a*fy, a*fx, -b*fy, 
                          k2Seq[t2] - 2*k2Seq[t2]*fy-b*fx), byrow = TRUE, nrow = 2, ncol = 2)
    if(Im(eigen(fM)$values[1]) == 0)
    {
      RsTemp <- c(RsTemp, sqrt(eigen(fM)$values[1]^2 + eigen(fM)$values[2]^2))
    }
    else
    {
      RsTemp <- c(RsTemp, sqrt(Re(eigen(fM)$values[1])^2 + Im(eigen(fM)$values[1])^2))
    }
  }
  k1Temp <- seq(k1Seq[t1], k1Seq[t1], length = length(RsTemp))
  points3d(x = k1Temp, y = k2Seq, z = RsTemp)
}

### export png of x vs y plot, k1 fixed, ALWAYS SET X/Y LIMITS

k1Fixed <- 2
k2Int <- 3.5
k2End <- 4.6
k2Seq <- seq(from = k2Int, to = k2End, by = (k2End-k2Int)/1599)
setwd("C:/Users/bciro/Desktop/Discrete Project/Videos/PNG Files")
dir.create(paste("k1 = ", k1Fixed, sep = ""))
for(t in 1:length(k2Seq))
{
  fx <- (k1Fixed*k2Seq[t]+a*k2Seq[t]-k2Seq[t]-a)/(k1Fixed*k2Seq[t]+a*b)
  fy <- 1-(1+b*fx)/k2Seq[t]
  x0 <- c(fx*runif(1, min = .96, max = .97))
  y0 <- c(fy*runif(1, min = .96, max = .97))
  tStop <- 1300
  for(t1 in 1:tStop)
  {
    x0 <- c(x0, k1Fixed*x0[t1]*(1-x0[t1])+a*x0[t1]*y0[t1])
    y0 <- c(y0, k2Seq[t]*y0[t1]*(1-y0[t1])-b*x0[t1]*y0[t1])
  }
  
  
  mainTitle <- paste("k1 = ", k1Fixed, "/k1 = ", k1Fixed, ", k2 = ", sprintf("%.13f",k2Seq[t]), sep = "")
  
  png(paste(mainTitle, ".png", sep = ""), width = 1000, height = 1000)
  plot(x = x0[300:length(x0)], y = y0[300:length(y0)], xlim = c(.38, 1), 
       ylim = c(.3,.8), main = mainTitle, cex = .4, xlab = "x_n", ylab = "y_n", pch = 16)
  points(x = fx, y = fy, col = "red", pch = 19, cex = 1)
  dev.off()
}

### export pngs of x vs y plot, k2 fixed, ALWAYS SET X/Y LIMITS

k2Fixed <- 3.55
k1Int <- 2.8
k1End <- 3.38
k1Seq <- seq(from = k1Int, to = k1End, by = (k1End-k1Int)/1599)
setwd("C:/Users/bciro/Desktop/Discrete Project/Videos/PNG Files")
dir.create(paste("k2 = ", k2Fixed, sep = ""))
for(t in 1:length(k1Seq))
{
  fx <- (k1Seq[t]*k2Fixed+a*k2Fixed-k2Fixed-a)/(k1Seq[t]*k2Fixed+a*b)
  fy <- 1-(1+b*fx)/k2Fixed
  x0 <- c(fx*runif(1, min = .96, max = .97))
  y0 <- c(fy*runif(1, min = .96, max = .97))
  # x0 <- c(fx*runif(1, min = .4, max = .6))
  # y0 <- c(fy*runif(1, min = .85, max = .9))
  tStop <- 1300
  for(t1 in 1:tStop)
  {
    x0 <- c(x0, k1Seq[t]*x0[t1]*(1-x0[t1])+a*x0[t1]*y0[t1])
    y0 <- c(y0, k2Fixed*y0[t1]*(1-y0[t1])-b*x0[t1]*y0[t1])
  }
  
  
  mainTitle <- paste("k2 = ", k2Fixed, "/k1 = ", sprintf("%.13f",k1Seq[t]), ", k2 = ", k2Fixed, sep = "")
  
  png(paste(mainTitle, ".png", sep = ""), width = 1000, height = 1000)
  plot(x = x0[300:length(x0)], y = y0[300:length(y0)], xlim = c(.15, 1.1), 
    ylim = c(.05,.75), main = mainTitle, cex = .4, xlab = "x_n", ylab = "y_n", pch = 16)
  points(x = fx, y = fy, col = "red", pch = 19, cex = 1)
  dev.off()
}
