rm(list=ls()) ; gc() ; graphics.off()

library(DistatisR)
library(devtools)
load_all("~/git/sGSVD/")


data(DistAlgo)
res.d <- distatis(DistAlgo)
M <- res.d$res4Cmat$C

res.eig <- eigen(M)
res.speig <- sparseEIGEN(X = M, k = 4L, rds = rep(1.25, 4), init = res.eig$vectors[, c(1, 1, 1, 1)], )

apply(DistAlgo * array(rep(res.speig$vectors[, 2], each = 36), c(6, 6, 4)), 1:2, sum)

array(rep(1:4, each = 9), c(3, 3, 4)) * array(rep(c(1, 0, 1, 0), each = 9), c(3, 3, 4))



res.speig2 <- sparseEIGENpos(X = M, k = 4L, rds = rep(1.6, 4), init = "svd")

crossprod(res.speig2$vectors)

x <- abs(rnorm(10))
U <- projOrth_then_projPos_then_projL1L2(abs(rnorm(10)), rds = 2, OrthSpace = rep(0, 10), grp = 1:10, itermax = 1000, eps = 1e-10)$x
# U <- projL1L2(rnorm(10), rds = 2)$x
c(L1 = normL1(U), L2 = normL2(U))
res <- projOrth_then_projPos_then_projL1L2(vec = x, rds = 2, OrthSpace = U, grp = 1:10, itermax = 1000, eps = 1e-10)
crossprod(x, U)
round(cbind(x, U), 3)



projOrth(vec = c(sqrt(2), sqrt(2)), OrthSpace = projL2(c(0.1, 0.7))$x)
