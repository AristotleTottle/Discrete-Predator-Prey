library(lattice)
k <- 1
xn <- c(.4)
equil1 <- 0
equil2 <- 1 - 1/k

timeStop <- 30

for(t in 1:timeStop)
{
  xn <- c(xn, k*xn[t]*(1-xn[t]))
}
xn1 <-c()
for(t in 1:timeStop)
{
  xn1 <-c(xn1, xn[t], xn[t])
}
xn1 <- c(xn1, xn[length(xn)])

xn2 <- c(0)
for(t in 1:timeStop)
{
  xn2 <-c(xn2, xn[t], xn[t])
}

###plot x_n vs n

str = paste("k = ", k, ", ", "x_0 = ", xn[1], sep = "")
plot(x = 0:timeStop, y = xn, xlim = c(0, timeStop), ylim = c(0, max(xn)), 
     xlab = "time n", ylab = "x_n", main = str)

for(t in 1:timeStop)
{
  arrows(x0 = t-1, y0 = xn[t], x1 = t, y1 = xn[t+1], code = 2, length = 0.01)
}
# abline(h = equil1, col = "red", lwd = 2)
# abline(h = equil2, col = "red", lwd = 2)

###plot x_(n+1) vs x_n

plot(x = xn1, y = xn2, xlim = c(0, max(xn1)), ylim = c(0, max(xn2)), 
     xlab = "x_n", ylab = "x_(n+1)", main = str)
points(x = equil1, y = equil1, col = "red", pch = 19, cex = 1)
points(x = equil2, y = equil2, col = "red", pch = 19, cex = 1)

for(t in 1:timeStop)
{
  arrows(x0 = xn1[t], y0 = xn2[t], x1 = xn1[t+1], y1 = xn2[t+1], code = 2, length = 0.07)
}

###plot x_n on one-dimensional line

stripplot(xn)

### orbit diagram

r1 <- 3.3
r2 <- 4
r <- seq(from = r1, to = r2, by = (r2-r1)/799)
rPoints <- data.frame(num <- 1:300)

for(t in 1:length(r))
{
  x0 <- c(runif(1))
  for(z in 1:799)
  {
    x0 <- c(x0, r[t]*x0[z]*(1-x0[z]))
  }
  rPoints[paste(t)] <- x0[501:800]
}

rPoints <- rPoints[-c(1)]

plot(x = NULL, xlim = c(r[1], r[length(r)]), ylim = c(0,1), xlab = "k", ylab = "x")

for(t in 1:length(r))
{
  rTemp <- seq(r[t], r[t], length = 300)
  points(x = rTemp, rPoints[,t], cex = .001)
}

### lyapunov exponent

lyapunov <- 0
n <- 100
for(t in 1:n)
{
  lyapunov <- lyapunov + log(abs(k-2*k*xn[t]))/n
}
lyapunov