rm(list = ls()) ; gc() ; graphics.off()

library(DistatisR)
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyr)
library(readr)
library(data4PCCAR)
library(PTCA4CATA)

# library(sGSVD)
devtools::load_all("~/git/sGSVD/")

###### REAL DATA ########
gord <- read_csv("../data/gordon333_label.csv")
systemlabel <- read_csv("../data/systemlabel_G3.csv")
load("../data/cfDiSTATIS_GordonParcel.RData")


## DiSTATIS
ZeCubes <- cubes$zcube
ZeCubes[is.infinite(ZeCubes)] <- 0.99
res.distatis <- distatis(ZeCubes, Distance = FALSE)
## Sparsify the compromise
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


## Eigen value decomposition
res.eig <- eigen(compmat)
ev <- res.eig$values

plot(ev, type = "h")

## Sparse eigen-decomposition for one sparsity parameter
res.speig <- sparseEIGEN(compmat, k = 10L, init = "svd", rds = sqrt(nrow(compmat)) * rep(0.6, 10))

parz <- seq(sqrt(nrow(compmat)), 1, length = 10)

res.list <- vector(mode = "list", length = length(parz))
init <- "svd"
for (i in seq_along(parz)) {
  res.list[[i]] <- sparseEIGEN(compmat, k = 10L, init = init, rds = rep(parz[i], 10))
  init <- res.list[[i]]$vectors
}

save(res.list, file = "res.list.RData")

nbzeros <- lapply(res.list,  \(el) cumsum(colSums(el$vectors == 0)))
zeroratio <- lapply(nbzeros, \(x) x / (1:10 * nrow(compmat)))
fitratio <- lapply(res.list,  \(el) cumsum(el$values) / cumsum(ev[1:10]))
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

res.opt <- sparseEIGEN(compmat, k = 10L, init = "svd", rds = rep(rdsopt, 10))

Fopt <- data.frame(
  CommName,
  res.opt$vectors %*% diag(sqrt(res.opt$values)))

tau <- 100 * res.opt$values / sum(ev)
labz <- sprintf("Dim. %i (%0.2f%%)", 1:10, tau)

ggplot(Fopt, aes(X1, X2, label = Community)) + 
  geom_hline(yintercept = 0, color = "darkorchid4", alpha = 0.6) + 
  geom_vline(xintercept = 0, color = "darkorchid4", alpha = 0.6) +
  geom_point(aes(color = I(CommColor)), size = 2) + 
  ggrepel::geom_text_repel(aes(color = I(CommColor))) +
  labs(x = labz[1], y = labz[2]) + 
  theme_light() +  
  theme(
    axis.title = element_text(color = "#42376B"), 
    axis.text = element_text(color = "#42376B"), 
    title = element_text(color = "#42376B"))


## PLOT ONLY THE MEANS
Fopt %>% 
  group_by(Community, CommColor) %>% 
  summarize(X1 = mean(X1), X2 = mean(X2)) %>%
  ggplot(aes(X1, X2, label = Community)) + 
  geom_hline(yintercept = 0, color = "darkorchid4", alpha = 0.6) + 
  geom_vline(xintercept = 0, color = "darkorchid4", alpha = 0.6) +
  geom_point(aes(color = I(CommColor)), size = 2) + 
  ggrepel::geom_text_repel(aes(color = I(CommColor))) +
  labs(x = labz[1], y = labz[2]) + 
  theme_light() +  
  theme(
    axis.title = element_text(color = "#42376B"), 
    axis.text = element_text(color = "#42376B"), 
    title = element_text(color = "#42376B"))

### Bootstrap for DiSTATIS
fi <- res.distatis$res4Splus$F
fi[, 1] <- -fi[, 1] ; fi[, 2] <- -fi[, 2]
distatis.boot <- Boot4Mean(fi,
                      CommName$Community, niter = 1000)

mfi <- distatis.boot$GroupMeans
mfi.boot <- distatis.boot$BootCube
colnames(fi) <- colnames(mfi) <- colnames(mfi.boot) <- paste0("Dimension ", c(1:ncol(fi)))

distatis.fi.plot <- createFactorMap(
  fi,
  col.points = CommName$CommColor,
  col.labels = CommName$CommColor,
  alpha.points = 0.1, cex = 2, text.cex = 3,
  col.background = NULL, 
  col.axes = "#42376B", 
  width.axes = 1, alpha.axes = 0.5)


distatis.mfi <- createFactorMap(
  mfi,
  col.points = CommColor.gc[rownames(mfi), "CommColor"],
  col.labels = CommColor.gc[rownames(mfi), "CommColor"],
  alpha.points = 0.7, cex = 2, text.cex = 3, 
  col.background = NULL, col.axes = "#42376B",
  width.axes = 1, alpha.axes = 0.5,
  constraints = minmaxHelper(fi[,c(1:2)]))

distatis.mfi.ci <- MakeCIEllipses(
  mfi.boot,
  col = CommColor.gc[rownames(mfi.boot), "CommColor"],
  alpha.ellipse = 0.1, 
  line.size = 0.5, alpha.line = 0.2)

distatis.label <- createxyLabels.gen(
  lambda = ev, tau = 100*(ev/sum(ev)))

## get plots
fmap.distatis <- distatis.mfi$zeMap_background + distatis.fi.plot$zeMap_dots + distatis.mfi.ci + distatis.mfi$zeMap_dots + distatis.mfi$zeMap_text + distatis.label + ggtitle("DiSTATIS: Factor Scores")
fmap.distatis

sdistatis.fi <- res.opt$vectors %*% diag(sqrt(res.opt$values))
sdistatis.fi[, 2] <- -sdistatis.fi[, 2]

sdistatis.boot <- Boot4Mean(sdistatis.fi,
                       CommName$Community, niter = 1000)

smfi <- sdistatis.boot$GroupMeans
smfi.boot <- sdistatis.boot$BootCube
colnames(sdistatis.fi) <- colnames(smfi) <- colnames(smfi.boot) <- paste0("Dimension ", c(1:ncol(sdistatis.fi)))

sdistatis.fi.plot <- createFactorMap(
  sdistatis.fi,
  col.points = CommName$CommColor,
  col.labels = CommName$CommColor,
  alpha.points = 0.1, cex = 2, text.cex = 3,
  col.background = NULL, 
  col.axes = "#42376B", 
  width.axes = 1, alpha.axes = 0.5)

sdistatis.mfi <- createFactorMap(
  smfi,
  col.points = CommColor.gc[rownames(smfi), "CommColor"],
  col.labels = CommColor.gc[rownames(smfi), "CommColor"],
  alpha.points = 1, cex = 2, text.cex = 3, force = FALSE,
  col.background = NULL, col.axes = "#42376B",
  width.axes = 1, alpha.axes = 0.5,
  constraints = minmaxHelper(sdistatis.fi[,c(1:2)]))

sdistatis.mfi.ci <- MakeCIEllipses(
  smfi.boot,
  col = CommColor.gc[rownames(smfi.boot), "CommColor"],
  alpha.ellipse = 0.1, 
  line.size = 0.5, 
  alpha.line = 0.2)

sdistatis.label <- createxyLabels.gen(
  lambda = res.opt$values,
  tau = 100*(res.opt$values/sum(ev)))

## get plots
# options(ggrepel.max.overlaps = 0)
fmap.sparsedistatis <- sdistatis.mfi$zeMap_background + sdistatis.fi.plot$zeMap_dots + sdistatis.mfi.ci + sdistatis.mfi$zeMap_dots + sdistatis.mfi$zeMap_text + sdistatis.label+ ggtitle("Sparse DiSTATIS: Factor Scores")

fmap.sparsedistatis


## Save the plots
save(fmap.distatis, file = "fmap.distatis.RData")
save(fmap.sparsedistatis, file = "fmap.sparsedistatis.RData")
