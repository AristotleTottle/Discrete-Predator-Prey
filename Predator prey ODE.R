library(phaseR)
k1 <- 0.5
k2 <- -0.5
a <- 1
b <- 1

p0 <- c(k1, a, k2, b)
# xy <- matrix(runif(8, min = 0, max = 3), nrow = 4, ncol = 2, byrow = FALSE)
xy <- matrix(c(2,2.5, -1,1, -0.7,-0.7, 1,-1), nrow = 4, ncol = 2, byrow = TRUE)
predatorPrey <- function(t, y, parameters) {
  dx <- parameters[1]*y[1]*(1-y[1]) + parameters[2]*y[1]*y[2]
  dy <- parameters[3]*y[2]*(1-y[2]) - parameters[4]*y[1]*y[2]
  list(c(dx, dy))
}

xEquil <- c(0, 0, 1, (k1*k2+k2)/(k1*k2+1))
yEquil <- c(0, 1, 0, (k1*k2-k1)/(k1*k2+1))
origin <- c(0,0)
yAxis <- c(0, 1)
xAxis <- c(1, 0)
fourth <- c((k1*k2+k2)/(k1*k2+1), (k1*k2-k1)/(k1*k2+1))

fourthM <- matrix(data = c(k1-2*k1*fourth[1]+a*fourth[2], a*fourth[1], -b*fourth[2], 
                  k2 - 2*k2*fourth[2]-b*fourth[1]), byrow = TRUE, nrow = 2, ncol = 2)
fourthEigen <- c(eigen(fourthM)$values[1], eigen(fourthM)$values[2])

str = paste("k_1 = ", k1, ", k_2 = ", k2, sep = "")
plot(x = NULL, xlim = c(-2,3), ylim = c(-2,3), xlab = "x", ylab = "y", main = str)
points(x = xEquil, y = yEquil, col = "red", cex = 3, pch = 20)
predPreyFlow <- flowField(predatorPrey, xlim = c(-2,3), ylim = c(-2,3), 
                          parameters = p0, arrow.type = "equal", frac = .75, points = 15, arrow.head = .10)
predPreyNull <- nullclines(predatorPrey, xlim = c(-2,3), ylim = c(-2,3), parameters = p0, 
                           col = c("blue", "darkgreen"))
predPreyTraj <- trajectory(predatorPrey, y0 = xy, tlim = c(0,20), parameters = p0, cex = 3)