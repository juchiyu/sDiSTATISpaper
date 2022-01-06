rm(list = ls()) ; gc() ; graphics.off()

library(DistatisR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

# library(sGSVD)
devtools::load_all("~/git/sGSVD/")

# Data
noise <- matrix(runif(100), 10, 10)
M <- toeplitz(10:1) + 1/2 * (noise + t(noise))
rownames(M) <- sprintf("row%02i", 1:10)
colnames(M) <- sprintf("col%02i", 1:10)

## Heatmap with ggplot
tibble(row = rownames(M), as_tibble(M)) %>% 
  pivot_longer(col01:col10, names_to = "col", values_to = "x") %>%
  ggplot(aes(col, row, fill = x)) + 
  geom_tile(color = "white") + 
  labs(x = "", y = "", fill = "") +
  scale_fill_viridis_c() + theme_bw() + coord_equal() +
  theme(axis.text.x = element_text(angle = 90))

## Sparse eigen-decomposition
res.speig <- sparseEIGEN(X = M, k = 5L, init = "svd", rds = rep(1.3, 5))

###### REAL DATA ########
gord <- read_csv("../data/gordon333_label.csv")
systemlabel <- read_csv("../data/systemlabel_G3.csv")
load("../data/cfDiSTATIS_GordonParcel.RData")


# res.eig <- sparseEIGEN(cubes$rcube[,,1], init = "svd", rds = 0.5 * sqrt(c(333, 333)))

## Sparse EIGEN on the RV coefficient matrix
## It's going to need too much work
ZeCubes <- cubes$zcube
ZeCubes[is.infinite(ZeCubes)] <- 0.99
res.distatis <- distatis(ZeCubes, Distance = FALSE)
## Easy : sparsify the compromise
rvmat <- res.distatis$res4Cmat$C
rownames(rvmat) <- sprintf("row%02i", 1:nrow(rvmat))
colnames(rvmat) <- sprintf("col%02i", 1:ncol(rvmat))

## Heatmap with ggplot
tibble(row = rownames(rvmat), as_tibble(rvmat)) %>% 
  pivot_longer(starts_with("col"), names_to = "col", values_to = "x") %>%
  ggplot(aes(col, row, fill = x)) + 
  geom_tile(color = "white") + 
  labs(x = "", y = "", fill = "") +
  scale_fill_viridis_c() + theme_bw() + coord_equal() +
  theme(axis.text.x = element_text(angle = 90))

## Eigen value decomposition
res.eig <- eigen(rvmat)
ev <- res.eig$values

## Sparse eigen-decomposition for one sparsity parameter
res.speig <- sparseEIGEN(rvmat, k = 10L, init = "svd", rds = sqrt(nrow(rvmat)) * rep(0.6, 10))

parz <- seq(sqrt(nrow(rvmat)), 1, length = 20)

res.list <- vector(mode = "list", length = length(parz))
init <- "svd"
for (i in seq_along(parz)) {
  res.list[[i]] <- sparseEIGEN(rvmat, k = 10L, init = init, rds = rep(parz[i], 10))
  init <- res.list[[i]]$vectors
}

nbzeros <- lapply(res.list,  \(el) cumsum(colSums(el$vectors == 0)))
zeroratio <- lapply(nbzeros, \(x) x / (1:10 * nrow(rvmat)))
fitratio <- lapply(res.list,  \(el) cumsum(el$values) / cumsum(ev))
# names(nbzeros) <- sprintf("Par%02i", seq_along(parz))

dat.si <- data.frame(
  k = rep(1:10, length(parz)),
  rds = rep(parz, each = 10),
  nbzeros = unlist(nbzeros),
  zeroratio = unlist(zeroratio),
  fitratio = unlist(fitratio)
)

dat.si %<>% mutate(SI = zeroratio * fitratio)

theta <- seq(pi, 3/2*pi, length.out = 150)

siplot <- dat.si %>% filter(k < 5) %>% 
  ggplot(aes(zeroratio, fitratio, color = factor(k))) + 
  geom_hline(yintercept = 1, color = "darkorchid4", alpha = 0.6, linetype = 3) + 
  geom_vline(xintercept = 1, color = "darkorchid4", alpha = 0.6, linetype = 3) + 
  lapply(seq(0.25, 1.25, by = 0.25), function(r) annotate("path",
                                                          x = 1 + r*cos(theta),
                                                          y = 1 + r*sin(theta),
                                                          color = "darkorchid4")) + 
  lapply(seq(0.125, 1.5, by = 0.25), function(r) annotate("path",
                                                          x = 1 + r*cos(theta),
                                                          y = 1 + r*sin(theta),
                                                          color = "darkorchid4", size = 0.2))  +
  geom_abline(intercept = 0, slope = 1, color = "darkorchid4", alpha = 0.6, linetype = 3) +
  geom_point(aes(color = factor(k), size = SI)) +
  geom_line(aes(color = factor(k))) +
  theme_bw() + 
  coord_equal(xlim = 0:1, ylim = 0:1) + 
  guides(fill = "none") +
  labs(x = "Zero ratio", y = "Fit ratio", size = "Sparsity\nIndex", color = "Number of\nDimensions") + 
  theme(
    axis.title = element_text(color = "#42376B"), 
    axis.text = element_text(color = "#42376B"), 
    title = element_text(color = "#42376B"),
    panel.border = element_rect(size = 1.5, color = "#42376B", fill = NA),
    # plot.margin = unit(c(0.1,0,0,0), unit = "cm"),
    panel.grid = element_blank())

siplot

rdsopt <- (dat.si %>% filter(k < 5) %>%
             slice(which.max(SI)))[, "rds"]
kopt <- (dat.si %>% filter(k < 5) %>%
           slice(which.max(SI)))[, "k"]

res.opt <- sparseEIGEN(rvmat, k = 10L, init = "svd", rds = rep(rdsopt, 10))

Fopt <- data.frame(
  ID = sprintf("Subj%02i", 1:10),
  res.opt$vectors %*% diag(sqrt(res.opt$values)))

tau <- 100 * res.opt$values / sum(ev)
labz <- sprintf("Dim. %i (%0.2f%%)", 1:10, tau)

ggplot(data.frame(Fopt), aes(X1, X2, label = ID)) + 
  geom_hline(yintercept = 0, color = "darkorchid4", alpha = 0.6) + 
  geom_vline(xintercept = 0, color = "darkorchid4", alpha = 0.6) +
  geom_point(size = 2, color = "#42376B") + 
  ggrepel::geom_text_repel(, color = "#42376B") +
  labs(x = labz[1], y = labz[2]) + 
  theme_light() +  
  theme(
    axis.title = element_text(color = "#42376B"), 
    axis.text = element_text(color = "#42376B"), 
    title = element_text(color = "#42376B"))
    
## Now let's decompose the correct matrix !
compmat <- res.distatis$res4Splus$Splus
rownames(compmat) <- sprintf("row%02i", 1:nrow(compmat))
colnames(compmat) <- sprintf("col%02i", 1:ncol(compmat))

## Heatmap with ggplot
tibble(row = rownames(compmat), as_tibble(compmat)) %>% 
  pivot_longer(starts_with("col"), names_to = "col", values_to = "x") %>%
  ggplot(aes(col, row, fill = x)) + 
  geom_tile() + 
  labs(x = "", y = "", fill = "") +
  scale_fill_viridis_c() + theme_bw() + coord_equal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank())

res.eig.comp <- eigen(compmat)

plot(res.eig.comp$values, type = "h")



