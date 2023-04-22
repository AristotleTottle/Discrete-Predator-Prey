library(lattice)
library(phaseR)
k <- c(1)

f = function(t, y, parameters) {
  list(c(parameters[1]*y-parameters[1]*y^2))
}

x <- seq(0, 3, .01)
plot(x, f(x), type = 'l', ylim = c(-3,3), xlab = "x", 
     ylab = "dx/dt", main = "k = 12")
abline(h = 0)
abline(v = 0)

xn <- c(.4)
equil1 <- 0
equil2 <- 1 - 1/k

timeStop <- 10

for(t in 1:timeStop)
{
  xn <- c(xn, k*xn[t]*(1-xn[t]))
}
xn1 <- xn[2:length(xn)]

plot(x = 0.9, y = 0, xlim = c(0,4), ylim = c(0,1), xlab = "k", ylab = "x_n", 
     col = "red", pch = 19)
points(x = c(0.9, 1.5, 3.1, 3.1, 3.5, 3.5, 3.5, 3.5), y = c(0, .333, .558, .764, 
                                .826, .500, .874, .382), col = "red", pch = 19)

###plot x_n vs x_{n+1}

plot(x = xn[1:timeStop], y = xn1, xlim = c(0, max(xn)), ylim = c(0, max(xn1)))
points(x = equil1, y = equil1, col = "red", pch = 19, cex = 1)
points(x = equil2, y = equil2, col = "red", pch = 19, cex = 1)

###plot n vs. x_n

plot(x = 0:timeStop, y = xn, ylim = c(0,1), xlab = "time n", ylab = "x_n", 
     main = paste("k = ", k, ", x_0 = ", xn[1], sep = ""), cex = .7)
for(t in 1:timeStop)
{
  arrows(x0 = t-1, y0 = xn[t], x1 = t, y1 = xn[t+1], code = 2, length = 0.0)
}

###plot x_n on one-dimensional line

stripplot(xn)

### orbit diagram

r1 <- 3.82
r2 <- 3.84
r <- seq(from = r1, to = r2, by = (r2-r1)/799)
rPoints <- data.frame(num <- 1:200)

for(t in 1:length(r))
{
  x0 <- c(runif(1))
  for(z in 1:599)
  {
    x0 <- c(x0, r[t]*x0[z]*(1-x0[z]))
  }
  rPoints[paste(t)] <- x0[401:600]
}

rPoints <- rPoints[-c(1)]

plot(x = NULL, xlim = c(r[1], r[length(r)]), ylim = c(0,1), xlab = "k", ylab = "x")

for(t in 1:length(r))
{
  rTemp <- seq(r[t], r[t], length = 200)
  points(x = rTemp, rPoints[,t], cex = .001)
}
points(x = c(0.9, 1.5, 3.1, 3.1, 3.5, 3.5, 3.5, 3.5), y = c(0, .333, .558, .764, 
                                .826, .500, .874, .382), col = "red", pch = 19)


abline(v = c(0.9, 1.5, 3.1, 3.5), col = "red")

### Lyapunov Exponent

rInt <- 0
rEnd <- 4
rSeq <- seq(from = rInt, to = rEnd, by = (rEnd - rInt)/1599)
lyapunov <- 0
plot(x = NULL, ylim = c(-4,1), xlim = c(rSeq[1], rSeq[length(rSeq)]),
     xlab = "k", ylab = "Lyapunov exponent")
abline(h = 0)
xTemp <- c()
for(t in 1:length(rSeq))
{
  xTemp <- c(runif(1))
  for(t1 in 1:99)
  {
    xTemp <- c(xTemp, rSeq[t]*xTemp[t1]*(1-xTemp[t1]))
  }
  lyapunov <- 0
  for(t1 in 1:length(xTemp))
  {
    lyapunov <- lyapunov + log(abs(rSeq[t]-2*rSeq[t]*xTemp[t1]))/length(xTemp)
  }
  points(x = rSeq[t], y = lyapunov, cex = 0.001)
}


### Phase Portrait, Flow Field

x0 <- c(0.1, 3, -0.1, 2)
str = paste("k = ", k, sep = "")
plot(x = NULL, xlim = c(0,3), ylim = c(-2,3), xlab = "t'", ylab = "x")
logisticFlow <- flowField(f, xlim = c(0,3), ylim = c(-2,3), parameters = k, system = "one.dim")
abline(h = 0, col = "red")
abline(h = 1, col = "red")
logisticTraj <- trajectory(f, y0 = x0, tlim = c(0,3), parameters = k, system = "one.dim")

logisticPhase <- phasePortrait(f, ylim = c(-1,2), parameters = k, xlab = "x", ylab = "dx/dt'")


