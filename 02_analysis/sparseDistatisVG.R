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



res.speig2 <- sparseEIGENpos(X = M, k = 4L, rds = rep(1.25, 4), init = res.eig$vectors[, c(1, 1, 1, 1)], projPriority = "")

