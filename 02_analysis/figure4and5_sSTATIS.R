library(dplyr)
library(ExPosition)
library(PTCA4CATA)
library(ggplot2)
# library(ggrepel)
# library(tidyverse)
# library(plotly)
# library(RColorBrewer)
# library(kableExtra)
# library(ggpubr)
# library(ggplotify)
# library(SPAFAC)
# library(sGSVD)
# library(inlmisc)
# library(pheatmap)
# library(grid)

## Source functions
source("functions/PlotFactor.R")
source("functions/createxyLabels.gen.digit.R")
source("functions/createLabel.gen.digit.R")
source("functions/plot.smca.results.R")
source("functions/PlotMyScree.R")
source("functions/PrettyBarPlot.MuSu.R")
source("functions/sparsity_helpers.R")
source("functions/sparsity_helpers.R")
source("functions/statis.R")
source("functions/SUMPCAnormMat.R")

## Load the data
load("01_data/MFA_turkeysensory.rda")
rownames(grand.tab) <- c("HORNEADO", "SABORI", "DE FUD", "CAMPESTRE", "VIRGINIA", "NATURAL", "CLÃSICA", "ALPINO")

## design by raters
raters <- sapply(strsplit(colnames(grand.tab), "\\."), "[", 2)
raters.coldx <- data.frame(t(raters))
colnames(raters.coldx) <- colnames(grand.tab)

raters.fac <- factor(raters, levels = unique(raters))

## design by flavor
flavors <- sapply(strsplit(colnames(grand.tab), "\\."), "[", 1)
flavors.coldx <- data.frame(t(flavors))
colnames(flavors.coldx) <- colnames(grand.tab)

## design by flavor when grand.tab is ordered by flavors
flavors.ord <- rep(unique(flavors), each = 8)
flavors.ord.coldx <- data.frame(t(flavors.ord))

## get color ------------------
col.raters <- list()
col.raters$gc <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#E8C245','#a65628','#f781bf') %>% as.matrix
col.raters$oc <- rep(col.raters$gc, each = 12) %>% as.matrix
rownames(col.raters$gc) <- unique(raters)
rownames(col.raters$oc) <- colnames(grand.tab)

## when the columns are ordered by flavors (instead of raters)
col.flavor.idx <- c("#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac", "#c3bc3f", "#8cc2ca")
names(col.flavor.idx) <- unique(flavors)
col.flavors <- list()
col.flavors$gc <- as.matrix(col.flavor.idx)
col.flavors$oc <- as.matrix(recode(flavors, !!!col.flavor.idx))
rownames(col.flavors$oc) <- colnames(grand.tab)


## get theme -----------------
ggtheme <- theme(
  axis.title = element_text(size = 6, color = "#42376B"),
  axis.text = element_text(size = 6, color = "#42376B"),
  title = element_text(size = 6, color = "#42376B"),
  panel.border = element_rect(linewidth = 1.5, color = "#42376B", fill = NA))
ggtheme.noborder <- theme(
  axis.title = element_text(size = 6, color = "#42376B"),
  axis.text = element_text(size = 6, color = "#42376B"),
  title = element_text(size = 6, color = "#42376B"))

##### STUCK HERE: MEXPOSITION id DEPRECATED
## -------------------------------------------------------------------
# res.sts.byraters <- mpSTATIS(grand.tab, column.design = raters.coldx, statis.prepro.option = 'Plain_STATIS', DESIGN = raters, graph = FALSE)
grand.tab.list <- lapply(split(data.frame(t(grand.tab)), raters), t)

res.sts.byraters <- statis(
  LaGrandeTable = grand.tab.list,
  center = TRUE, 
  scale = "SS1",
  Norm = "SUMPCA")

# USE STATIS WITH SUMPCA NORMALZATION AND KEEP SS1


# res.sts.byflavor <- mpSTATIS(grand.tab, column.design = flavors.coldx, statis.prepro.option = 'Plain_STATIS', DESIGN = flavors, graph = FALSE)

# ## reverse components 1 and 2 from the compromise
# res.sts.byraters$mexPosition.Data$Table$fi[,c(1:2)] <- res.sts.byraters$mexPosition.Data$Table$fi[,c(1:2)]*-1
# res.sts.byraters$mexPosition.Data$Table$partial.fi.array[,c(1:2),] <- res.sts.byraters$mexPosition.Data$Table$partial.fi.array[,c(1:2),]*-1
# res.sts.byraters$mexPosition.Data$Table$Q[,c(1:2)] <- res.sts.byraters$mexPosition.Data$Table$Q[,c(1:2)]*-1

grand.tab.preproc <- lapply(lapply(grand.tab.list, function(x) expo.scale(x, center = TRUE, scale = "SS1")), SUMPCAnormMat)



## -------------------------------------------------------------------
eigres.sts <- data.frame(
  eig = res.sts.byraters$res4Cmat$eigValues,
  tau = res.sts.byraters$res4Cmat$tau)

rv.sts.scree <- 
  PlotMyScree(eigres.sts, lwd = 0.5, cex = 1.3, text.cex = 6) + 
  ggtitle("Scree plot for non-sparsed\nRV matrix from STATIS") + 
  ggtheme + 
  theme(
    plot.title = element_text(face = "plain"),
    panel.border = element_rect(linewidth = 1, color = "#42376B", fill = NA))


## -------------------------------------------------------------------
options(ggrepel.max.overlap = Inf)
sts.labels <- createxyLabels.gen.digit(
  1, 2,
  lambda = res.sts.byraters$res4Cmat$eigValues,
  tau = res.sts.byraters$res4Cmat$tau,
  axisName = "Component ")

## sts.constraints
sts.fi <- res.sts.byraters$res4Splus$Fi
sts.pfi <- res.sts.byraters$res4Splus$PartialFi
colnames(sts.fi) <- colnames(sts.pfi) <- paste0("Dimension ", c(1:ncol(sts.fi)))
rownames(sts.pfi) <- rownames(sts.fi)
dimnames(sts.pfi)[[3]] <- unique(raters)
sts.constraints <- minmaxHelper4Partial(sts.fi, sts.pfi)

## Splus.U1
### Fi
Splus.fi <- createFactorMap(sts.fi,
                            constraints = sts.constraints,
                            text.cex = 1.5,
                            cex = 1.5,
                            col.points = "#42376B",
                            col.labels = "#42376B",
                            col.background = NULL,
                            col.axes = "#42376B",
                            width.axes = 0.5,
                            label.axisName = "Component ",
                            alpha.axes = 0.5,
                            title = "STATIS: global and partial\nfactor scores of products (rows)")
Splus.pfi <- createPartialFactorScoresMap(sts.fi,
                                          sts.pfi,
                                          colors4Blocks = col.raters$gc,
                                          colors4Items = "#42376B",size.points = 0.5,
                                          alpha.lines = 0.5, alpha.points = 0.7, size.lines = 0.5)

sts.cmpfi.plot <- Splus.fi$zeMap_background + Splus.pfi$linesColByBlocks + Splus.pfi$pointsColByBlocks + Splus.fi$zeMap_dots + Splus.fi$zeMap_text + sts.labels + ggtheme

### Fj
Splus.fj <- createFactorMap(
  do.call("rbind", res.sts.byraters$res4Splus$Fj),
  col.labels = col.flavors$oc,
  col.points = col.flavors$oc,
  text.cex = 2,
  col.background = NULL,
  col.axes = "#42376B",
  width.axes = 0.5,
  label.axisName = "Component ",
  alpha.axes = 0.5,
  title = "STATIS: column factor scores of \nraters and flavors (columns)")
sts.cmpfj.plot <- Splus.fj$zeMap + sts.labels + ggtheme

### Barplots - fj
sts.cmpfj.dat <- data.frame(do.call("rbind", res.sts.byraters$res4Splus$Q)[,c(1,2)], 
                            flavors = flavors, raters = raters, 
                            ID = rownames(do.call("rbind", res.sts.byraters$res4Splus$Fj)), 
                            IDnum = c(1:length(flavors)))
colnames(sts.cmpfj.dat)[1:2] <- paste0("Dimension.", c(1:2))

(sts.cmpfj.bar1 <- PrettyBarPlot.MuSu(sts.cmpfj.dat$IDnum, sts.cmpfj.dat$Dimension.1, grp = raters, 
                                      color = drop(col.raters$gc[unique(sts.cmpfj.dat$raters),]), 
                                      label.bar = flavors,
                                      grp.label.size = 1.5,
                                      grp.label.position = 1,
                                      bar.label.size = 1.25,
                                      title = "STATIS: column factor scores",
                                      ylab = "Component 1") + ylim(c(-1,1)) + 
    ggtheme.noborder)

(sts.cmpfj.bar2 <- PrettyBarPlot.MuSu(sts.cmpfj.dat$IDnum, sts.cmpfj.dat$Dimension.2, grp = raters, 
                                      color = drop(col.raters$gc[unique(sts.cmpfj.dat$raters),]), 
                                      label.bar = flavors,
                                      grp.label.size = 1.5,
                                      grp.label.position = 1.22,
                                      bar.label.size = 1.25,
                                      title = "STATIS: column factor scores",
                                      ylab = "Component 2") + ylim(c(-1.25,1.22)) + 
    ggtheme.noborder)


## -------------------------------------------------------------------
RVmat <- res.sts.byraters$mexPosition.Data$InnerProduct$RVMatrix
rownames(RVmat) <- colnames(RVmat) <- unique(raters)

# pheatmap(res.sts.byraters$mexPosition.Data$InnerProduct$RVMatrix, clustering_method = 'ward.D2')

disRV2 <- dist(res.sts.byraters$mexPosition.Data$InnerProduct$fi[,c(1:2)])
hc <- hclust(disRV2, method = "ward.D2")
hc.3 <- cutree(hc, k = 3)
names(hc.3) <- sub(".", "", names(hc.3))

## hc.3 colors
hc.col <- list()
hc.col.idx <- c("1" = "#6388b4",
                "2" = "#ffae34",
                "3" = "#ef6f6a")
hc.col$oc <- recode(hc.3, !!!hc.col.idx) %>% as.matrix
rownames(hc.col$oc) <- names(hc.3)
hc.col$gc <- hc.col.idx %>% as.matrix

hc.dx <- hc.3 %>% factor %>% data.frame
colnames(hc.dx) <- "cluster"
hc.col.list <- list(cluster = drop(hc.col$gc))


## heatmap ??? NOT WORKING!!!
rv.heat <- pheatmap(RVmat, clustering_distance_rows = disRV2, clustering_distance_cols = disRV2, clustering_method = 'ward.D2', color = colorRampPalette(c("white", "darkorchid4"))(50),
                    annotation_row = hc.dx,
                    annotation_col = hc.dx,
                    annotation_colors = hc.col.list, annotation_legend = FALSE, fontsize = 5,
                    treeheight_col = 3, treeheight_row = 3, angle_col = "90", 
                    legend.width = 5, cellwidth = 5, cellheight = 5)
rv.heat$gtable$grobs[[5]]$gp <- gpar(col = col.raters$gc[hc$order,])
rv.heat$gtable$grobs[[4]]$gp <- gpar(col = col.raters$gc[hc$order,])
rv.heat$gtable$grobs[[10]]$children[[1]]$width <- rv.heat$gtable$grobs[[10]]$children[[1]]$width/2
rv.heat$gtable$grobs[[10]]$children[[1]]$height <- rv.heat$gtable$grobs[[10]]$children[[1]]$height*0.6
rv.heat$gtable$grobs[[10]]$children[[1]]$y <- rv.heat$gtable$grobs[[10]]$children[[1]]$y*0.6 + unit(0.3, units = "npc")
rv.heat$gtable$grobs[[10]]$children[[2]]$y <- rv.heat$gtable$grobs[[10]]$children[[2]]$y*0.6 + unit(0.3, units = "npc")
rv.heat$gtable$grobs[[10]]$children[[1]]$x <- rv.heat$gtable$grobs[[10]]$children[[1]]$x - unit(0.3, units = "npc")
rv.heat$gtable$grobs[[10]]$children[[2]]$x <- rv.heat$gtable$grobs[[10]]$children[[2]]$x - unit(0.6, units = "npc")
# rv.heat$gtable$grobs[[10]]$children[[1]]$vjust <- unit(1, units = "cm")
# rv.heat$gtable$grobs[[10]]$children[[2]]$vjust <- unit(10, units = "npc")
# rv.heat$gtable$grobs[[10]]$children[[1]]$hjust <- unit(0.8, units = "npc")
# rv.heat$gtable$grobs[[10]]$children[[2]]$hjust <- rep(unit(1.5, units = "npc"), 5)

## plot RV factor scores
rv.fimat <- res.sts.byraters$mexPosition.Data$InnerProduct$fi
rownames(rv.fimat) <- sub(".", "", rownames(rv.fimat))
rv.fimat[,1] <- rv.fimat[,1]*-1

rv.fi <- createFactorMap(rv.fimat,
                         text.cex = 1.5,
                         cex = 1,
                         col.points = col.raters$gc[rownames(rv.fimat),],
                         col.labels = col.raters$gc[rownames(rv.fimat),],
                         # col.points = hc.col$oc,
                         # col.labels = hc.col$oc,
                         col.background = NULL,
                         col.axes = "#42376B",
                         width.axes = 0.5,
                         label.axisName = "Component ",
                         alpha.axes = 0.5,
                         title = "STATIS: factor scores\nof the RV matrix")
rv.label <- createxyLabels.gen(lambda = res.sts.byraters$mexPosition.Data$InnerProduct$eigs,
                               tau = res.sts.byraters$mexPosition.Data$InnerProduct$t, 
                               axisName = "Component ")

rv.fi.plot <- rv.fi$zeMap + rv.label + ggtheme + theme(panel.border = element_rect(size = 1, color = "#42376B", fill = NA), plot.title = element_text(size = 6))

## plot barplot for alpha values
alpha2plot <- res.sts.byraters$mexPosition.Data$InnerProduct$alphaWeights %>% as.matrix
rownames(alpha2plot) <- sub(".", "", rownames(alpha2plot))
colnames(alpha2plot) <- "Weights for each tables (\alpha)"
rv.alpha <- PrettyBarPlot2(alpha2plot[,1], threshold = 0, color4bar = col.raters$gc[order(alpha2plot)], sortValues = TRUE, font.size = 2.5) +
  ggtitle("Table weights") +
  ylab("Weights") + ggtheme.noborder


## -------------------------------------------------------------------
## Subspace 1
I.c <- nrow(RVmat)

rv.eig <- eigen(RVmat)

K <- 2L
U0 <- abs(rv.eig$vectors[, 1:K])
res.speig <- NULL
params <- seq(sqrt(I.c), 1, length = 25)

for (i in seq_along(params)) {
  cat(sprintf("Radius: %0.2f\n", params[i]))
  res.speig[[i]] <- tryCatch(
    sparseEIGEN(X = RVmat, k = K, init = U0, rds = rep(params[i], K)),
    error = function(e) NA)
  if (!is.na(res.speig[[i]])){
    U0 <- res.speig[[i]]$vectors
  }
}


res.speig.sparseindex <- lapply(res.speig, sparseIndexEigen, eigenValues = rv.eig$values)

si.mat.tmp <- sapply(res.speig.sparseindex, 
                     \(.) c(SI = .$SI[1], 
                            fitRatio = .$fitRatio[1],
                            zeroRatio = .$zeroRatio[1])) 
si.dat <- data.frame(
  radius = params,
  t(si.mat.tmp)
)

si.dat$max <- "FALSE"
si.dat$max[which.max(si.dat$SI)] <- "MAX"

theta <- seq(pi, 3/2*pi, length.out = 150)
siplot1 <- ggplot(si.dat, aes(zeroRatio, fitRatio)) + 
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
  geom_point(aes(size = SI*0.1)) + geom_line() +
  guides(color = guide_legend(override.aes = list(shape = 19, color = "black"), nrow = 1, byrow = TRUE), shape = 19) +
  theme_bw() + 
  scale_size_continuous(limits = 0:1) + 
  coord_equal(xlim=0:1, ylim = 0:1) + 
  labs(x = "Zero ratio", y = "Fit", size = "Sparsity\nIndex") + 
  ggtitle("Subspace 1") +
  theme(
    axis.title = element_text(size = 6, color = "#42376B"), 
    axis.text = element_text(size = 6, color = "#42376B"), 
    title = element_text(size = 6, color = c("#42376B", "#FFFFFF", "#FFFFFF")),
    legend.key.size = unit(0.1, "cm"), legend.text = element_text(size = 6), legend.position = "",
    # legend.spacing.x = unit(0.05, "cm"),
    # legend.margin = margin(0,0,0,0), legend.box.margin = margin(-5,0,0,0),
    # legend.title.align = 0, legend.title = element_text(size = 6),
    panel.border = element_rect(size = 1.5, color = "#42376B", fill = NA),
    plot.margin = unit(c(0.1,0,0,0), unit = "cm"),
    panel.grid = element_blank()) + 
  with(si.dat[si.dat$max == "MAX",], annotate(geom = "point", x = zeroRatio, y = fitRatio, color = "#7F6CC7", size = 1)) +
  with(si.dat[si.dat$max == "MAX",], annotate(geom = "segment", x = zeroRatio + 0.25, y = fitRatio + 0.2, xend = zeroRatio+ 0.02, yend = fitRatio + 0.02, arrow = arrow(length = unit(0.05, "inches"), type = "closed"), color = "darkorchid4", size = 0.6)) + #0.05, 0.5
  # with(si.dat[si.dat$max == "MAX",], annotate(geom = "segment", x = zeroRatio + 0.025, y = fitRatio + 0.022, xend = zeroRatio+ 0.02, yend = fitRatio + 0.02, arrow = arrow(length = unit(0.04, "inches"), type = "closed"), color = "#7F6CC7", size = 0.01)) + #0.04, 0.01
  with(si.dat[si.dat$max == "MAX",], annotate(geom = "label", x = zeroRatio + 0.26, y = fitRatio + 0.18, color = "darkorchid4", label = substring(sprintf("%.3f", SI), 2), fill = "#FAFAFA", label.padding = unit(0.04, "in"), size = 2))

siplot1

rds.opt1 <- si.dat$radius[which.max(si.dat$SI)]

U0 <- abs(rv.eig$vectors[, 1:K])
rv.speig1 <- sparseEIGEN(X = RVmat, k = K, init = U0, rds = rep(rds.opt1, K))

## check
table(rv.speig1$vectors[,1] != 0, hc.3)

## Subspace 2
torem <- which(rv.speig1$vectors[, 1] != 0)
RVmat2 <- RVmat[-torem, -torem]
I.c2 <- nrow(RVmat2)
rv.eig2 <- eigen(RVmat2)

U0 <- abs(rv.eig2$vectors[, 1:K])
res.speig2 <- NULL
params2 <- seq(sqrt(I.c2), 1, length = 25)

for (i in seq_along(params2)) {
  cat(sprintf("Radius: %0.2f\n", params2[i]))
  res.speig2[[i]] <- tryCatch(
    sparseEIGEN(X = RVmat2, k = K, init = U0, rds = rep(params2[i], K)),
    error = function(e) NA)
  if (!is.na(res.speig2[[i]])){
    U0 <- res.speig2[[i]]$vectors
  }
}


res.speig.sparseindex2 <- lapply(res.speig2, sparseIndexEigen, eigenValues = rv.eig2$values)

si.mat.tmp2 <- sapply(res.speig.sparseindex2, 
                      \(.) c(SI = .$SI[1], 
                             fitRatio = .$fitRatio[1],
                             zeroRatio = .$zeroRatio[1])) 
si.dat2 <- data.frame(
  radius = params2,
  t(si.mat.tmp2)
)

si.dat2$max <- "FALSE"
si.dat2$max[which.max(si.dat$SI)] <- "MAX"

siplot2 <- ggplot(si.dat2, aes(zeroRatio, fitRatio)) + 
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
  geom_point(aes(size = SI*0.1)) + geom_line() +
  guides(color = guide_legend(override.aes = list(shape = 19, color = "black"), nrow = 1, byrow = TRUE), shape = 19) +
  theme_bw() + 
  scale_size_continuous(limits = 0:1) + 
  coord_equal(xlim=0:1, ylim = 0:1) + 
  labs(x = "Zero ratio", y = "Fit", size = "Sparsity\nIndex") + 
  ggtitle("Subspace 2") +
  theme(
    axis.title = element_text(size = 6, color = "#42376B"), 
    axis.text = element_text(size = 6, color = "#42376B"), 
    title = element_text(size = 6, color = c("#42376B", "#FFFFFF", "#FFFFFF")),
    legend.key.size = unit(0.1, "cm"), legend.text = element_text(size = 6), legend.position = "",
    # legend.spacing.x = unit(0.05, "cm"),
    # legend.margin = margin(0,0,0,0), legend.box.margin = margin(-5,0,0,0),
    # legend.title.align = 0, legend.title = element_text(size = 6),
    panel.border = element_rect(size = 1.5, color = "#42376B", fill = NA),
    plot.margin = unit(c(0.1,0,0,0), unit = "cm"),
    panel.grid = element_blank()) + 
  with(si.dat2[si.dat2$max == "MAX",], annotate(geom = "point", x = zeroRatio, y = fitRatio, color = "#7F6CC7", size = 1)) +
  with(si.dat2[si.dat2$max == "MAX",], annotate(geom = "segment", x = zeroRatio + 0.25, y = fitRatio + 0.2, xend = zeroRatio+ 0.02, yend = fitRatio + 0.02, arrow = arrow(length = unit(0.05, "inches"), type = "closed"), color = "darkorchid4", size = 0.6)) + #0.05, 0.5
  # with(si.dat2[si.dat2$max == "MAX",], annotate(geom = "segment", x = zeroRatio + 0.025, y = fitRatio + 0.022, xend = zeroRatio+ 0.02, yend = fitRatio + 0.02, arrow = arrow(length = unit(0.04, "inches"), type = "closed"), color = "#7F6CC7", size = 0.01)) + #0.04, 0.01
  with(si.dat2[si.dat2$max == "MAX",], annotate(geom = "label", x = zeroRatio + 0.26, y = fitRatio + 0.18, color = "darkorchid4", label = substring(sprintf("%.3f", SI), 2), fill = "#FAFAFA", label.padding = unit(0.04, "in"), size = 2))

siplot2

rds.opt2 <- si.dat2$radius[which.max(si.dat2$SI)]

U0 <- abs(rv.eig2$vectors[, 1:K])
rv.speig2 <- sparseEIGEN(X = RVmat2, k = K, init = U0, rds = rep(rds.opt2, K))

## Third subspace
torem2 <- which(rv.speig2$vectors[, 1] != 0)
RVmat3 <- RVmat2[-torem2, -torem2]
rv.eig3 <- eigen(RVmat3)

## Organize
U1 <- rv.speig1$vectors[, 1]
U3 <- U2 <- rep(0, I.c)
U2[-torem] <- rv.speig2$vectors[, 1]
U3[-torem][-torem2] <- abs(rv.eig3$vectors[, 1])

Usp <- data.frame(U1 = U1, U2 = U2, U3 = U3)
rownames(Usp) <- unique(raters)

## Scree plot
eig.srv <- c(sum(rv.eig$values), sum(rv.eig2$values), sum(rv.eig3$values))
eigres.srv <- data.frame(eig = eig.srv, tau = eig.srv/sum(eig.srv))
rv.ssts.scree <- PlotMyScree(eigres.srv, lwd = 0.5, cex = 1.3, title = "Scree plot for the sparsified\nRV space from sSTATIS", text.cex = 5) + theme(axis.title = element_text(size = 5, color = "#42376B"), axis.text = element_text(size = 5, color = "#42376B"), title = element_text(face = "plain", size = 5, color = "#42376B"), plot.title = element_text(face = "plain"), panel.border = element_rect(size = 0.5, color = "#42376B", fill = NA))

## check
class <- Usp %>% transmute(class = ifelse(U1 > 0, "U1", ifelse(U2 > 0, "U2", "U3")))
table(class$class, hc.3)


## -------------------------------------------------------------------
dat.U <- data.frame(
  ID = rownames(Usp),
  # IDnum = as.numeric(factor(rownames(Usp))),
  Usp, cluster = factor(hc.3))

dat.U4lol <- dat.U[hc$order,]
dat.U4lol$IDnum <- c(1:nrow(dat.U4lol))
colnames(dat.U4lol)[2:4] <- c(paste0("Subspace ", c(1:3)))

tryplot <- dat.U4lol %>% pivot_longer(cols = `Subspace 1`:`Subspace 3`,
             names_to = "vector", names_prefix = "",
             values_to = "Weights") %>% data.frame

sprv.loli <- ggplot(tryplot, aes(IDnum, Weights, color = cluster)) + 
  geom_point(size = 1) + 
  geom_segment(aes(xend = IDnum, y = 0, yend = Weights)) + 
  facet_wrap(~ vector, ncol = 1, strip.position = "top") +
  color_palette(palette = hc.col$gc) +
  theme_minimal() + labs(x = "") + 
  theme(axis.text.x = element_blank(), 
        legend.position = "NA",
        axis.text = element_text(size = 6, color = "#42376B"),
        axis.title = element_text(size = 6, color = "#42376B"),
        strip.text = element_text(size = 6, color = "#42376B"),
        plot.margin = unit(c(0.1,0,0,0), unit = "cm"),
        strip.background = element_rect(color = "transparent", fill = "transparent", size = 0.1))


# lol2 <- ggplot(dat.U4lol, aes(IDnum, U2, color = cluster)) + 
#   geom_point() + 
#   geom_segment(aes(xend = IDnum, y = 0, yend = U2)) + 
#   color_palette(palette = hc.col$gc) +
#   theme_minimal()  + labs(x = "")+ theme(axis.text.x = element_blank())
# 
# lol3 <- ggplot(dat.U4lol, aes(IDnum, U3, color = cluster)) + 
#   geom_point() + 
#   geom_segment(aes(xend = IDnum, y = 0, yend = U3)) + 
#   color_palette(palette = hc.col$gc) +
#   theme_minimal()  + labs(x = "")+ theme(axis.text.x = element_blank())
# 
# ggpubr::ggarrange(lol1, lol2, lol3, nrow = 3, common.legend = TRUE)



## -------------------------------------------------------------------
## alpha
dat.alpha <- t(t(dat.U[, c("U1", "U2", "U3")])/colSums(dat.U[, c("U1", "U2", "U3")])) %>% as.data.frame
alpha.mat <- dat.alpha[raters,]

which1 <- alpha.mat$U1 != 0
which2 <- alpha.mat$U2 != 0
which3 <- alpha.mat$U3 != 0

## Compute Splus
#!! the columns are reordered
# Splus.U1 <- SPAFAC:::ComputeSplus(grand.tab.listbyflavor, 1/sqrt(dat.alpha$U1), isCube = FALSE)
# Splus.U2 <- SPAFAC:::ComputeSplus(grand.tab.listbyflavor, 1/sqrt(dat.alpha$U2), isCube = FALSE)
# Splus.U3 <- SPAFAC:::ComputeSplus(grand.tab.listbyflavor, 1/sqrt(dat.alpha$U3), isCube = FALSE)

Splus1.svd <- GSVD::gsvd(grand.tab.preproc[,which1], LW = rep(1/nrow(grand.tab), nrow(grand.tab)), RW = alpha.mat$U1[which1])
Splus2.svd <- GSVD::gsvd(grand.tab.preproc[,which2], LW = rep(1/nrow(grand.tab), nrow(grand.tab)), RW = alpha.mat$U2[which2])
Splus3.svd <- GSVD::gsvd(grand.tab.preproc[,which3], LW = rep(1/nrow(grand.tab), nrow(grand.tab)), RW = alpha.mat$U3[which3])

colnames(Splus1.svd$u) <- colnames(Splus1.svd$v) <- colnames(Splus1.svd$p) <- colnames(Splus1.svd$q) <- paste0("Dimension ", c(1:ncol(Splus1.svd$u)))
colnames(Splus2.svd$u) <- colnames(Splus2.svd$v) <- colnames(Splus2.svd$p) <- colnames(Splus2.svd$q) <- paste0("Dimension ", c(1:ncol(Splus2.svd$u)))
colnames(Splus3.svd$u) <- colnames(Splus3.svd$v) <- colnames(Splus3.svd$p) <- colnames(Splus3.svd$q) <- paste0("Dimension ", c(1:ncol(Splus3.svd$u)))

Splus2.svd$u[,2] <- Splus2.svd$u[,2]*-1
Splus2.svd$v[,2] <- Splus2.svd$v[,2]*-1
Splus2.svd$p[,2] <- Splus2.svd$p[,2]*-1
Splus2.svd$q[,2] <- Splus2.svd$q[,2]*-1

## organize output
#================= FUNC =================
STATIS.Sout <- function(dattab, svdres, M, alpha, design){
  out <- list(svd = list(u = svdres$u, v = svdres$v, d = svdres$d),
              gsvd = list(p = svdres$p, q = svdres$q, d = svdres$d, M = M, W = alpha),
              alpha = alpha,
              eig = svdres$l,
              tau = svdres$l/sum(svdres$l),
              fi = t(t(svdres$p) * svdres$d),
              fj = t(t(svdres$q) * svdres$d),
              partial.fi = array(dim = c(dim(svdres$u), length(unique(design)))))
  design.fac <- factor(design, unique(design))
  dattab.list <- split(as.data.frame(t(dattab)), design.fac) %>% lapply(t)
  dattab.vlist <- split(as.data.frame(svdres$q), design.fac) %>% lapply(as.matrix)
  out$partial.fi <- mapply("%*%", dattab.list, dattab.vlist) %>% array(dim = c(dim(svdres$u), length(unique(design))), 
                                                                       dimnames = list(rownames(svdres$u), 
                                                                                       colnames(svdres$u), 
                                                                                       unique(design)))
  out$svd.list$X <- dattab.list
  out$svd.list$v <- dattab.vlist
  return(out)
}
#========================================

# get factorized design
Splus1.res <- STATIS.Sout(grand.tab.preproc[,which1], Splus1.svd, rep(1/nrow(grand.tab), nrow(grand.tab)), alpha.mat$U1[which1], raters[which1])
Splus2.res <- STATIS.Sout(grand.tab.preproc[,which2], Splus2.svd, rep(1/nrow(grand.tab), nrow(grand.tab)), alpha.mat$U2[which2], raters[which2])
Splus3.res <- STATIS.Sout(grand.tab.preproc[,which3], Splus3.svd, rep(1/nrow(grand.tab), nrow(grand.tab)), alpha.mat$U3[which3], raters[which3])

##*** CHECK ***
##*SPAFAC:::ComputeSplus(Splus1.res$partial.fi, unique(alpha.mat$U1[which1]), isCube = TRUE)
##*Splus1.res$fi


## -------------------------------------------------------------------
eigres.srv.cmp <- list()
eigres.srv.cmp$subsp1 <- data.frame(eig = Splus1.res$eig,
                                    tau = Splus1.res$tau)
eigres.srv.cmp$subsp2 <- data.frame(eig = Splus2.res$eig,
                                    tau = Splus2.res$tau)
eigres.srv.cmp$subsp3 <- data.frame(eig = Splus3.res$eig,
                                    tau = Splus3.res$tau)

# srv.cmp.scree1 <- PlotMyScree(eigres.srv.cmp$subsp1, lwd = 0.5, cex = 1.3, title = "Subspace 1") + ggtheme + theme(plot.title = element_text(size = 4, font = "plain"))
# srv.cmp.scree2 <- PlotMyScree(eigres.srv.cmp$subsp2, lwd = 0.5, cex = 1.3, title = "Subspace 2") + ggtheme + theme(plot.title = element_text(size = 4, font = "plain"))
# srv.cmp.scree3 <- PlotMyScree(eigres.srv.cmp$subsp3, lwd = 0.5, cex = 1.3, title = "Subspace 3") + ggtheme + theme(plot.title = element_text(size = 4, font = "plain"))

srv.cmp.long <- rbind(eigres.srv.cmp$subsp1, eigres.srv.cmp$subsp2, eigres.srv.cmp$subsp3)
srv.cmp.long$Subspace <- c(rep("Subspace 1", 7), rep("Subspace 2", 7), rep("Subspace 3", 7))
srv.cmp.long$x <- rep(c(1:7), 3)

sprv.cmp.scree.adj <- ggplot(srv.cmp.long, aes(x = x, y = tau)) +
      geom_line(color = "grey40", size = 0.5) +
      geom_point(color = "#42376B", size = 1.1) +
  facet_wrap(~Subspace, ncol = 1) +
      scale_x_continuous(breaks=c(1:7)) +
  scale_y_continuous(name = bquote(atop(bold(), paste('\n\n    Percentage of \nvariance explained (%)'))),
                       sec.axis = sec_axis(~.*(srv.cmp.long$eig[1]/srv.cmp.long$tau[1]), name = "Pseudo-eigenvalues")) +
    xlab("Components") +
    theme(text = element_text(size = 6, color = "#42376B"),
          legend.position = "none",
          strip.text = element_text(face = "plain", color = "#42376B", size = 6),
          axis.text.y.left = element_text(face = "plain", color = "#42376B", angle = 90, hjust = 0.5, size = 4),
          axis.text.y.right = element_text(face = "plain", color = "#42376B", angle = 270, hjust = 0.5, size = 3),
          axis.text.x = element_text(size = 6, color = "#42376B"),
          panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0.1,0,0,0), unit = "cm"),
          panel.border = element_rect(color = "#42376B", fill = "transparent"),
          strip.background = element_rect(fill = "transparent", size = 0.1))

# sprv.cmp.scree <- ggpubr::ggarrange(srv.cmp.scree1 + rremove("xlab") + rremove("ylab"), 
#                                     srv.cmp.scree2 + rremove("xlab") + rremove("ylab"), 
#                                     srv.cmp.scree3 + rremove("xlab") + rremove("ylab"), nrow = 3, common.legend = TRUE, legend = "bottom")

# sprv.cmp.scree.adj <- annotate_figure(sprv.cmp.scree, left = grid::textGrob(bquote(atop(bold(), paste('\n    Percentage of  variance explained (%)'))), rot = 90, vjust = 1, gp = grid::gpar(cex = 0.5)),
#                                       right = grid::textGrob("Pseudo-eigenvalues", rot = 270, vjust = 1, gp = grid::gpar(cex = 0.5)),
#                                       bottom = grid::textGrob("Components", vjust = -1, gp = grid::gpar(cex = 0.5)))
# 
sprv.cmp.scree.adj



## -------------------------------------------------------------------
## Labels
stsU1.labels <- createxyLabels.gen.digit(1,2,
                                   lambda = Splus1.res$eig,
                                   tau = Splus1.res$tau,
                                   digit4tau = 2,
                                   axisName = "Component ")

stsU2.labels <- createxyLabels.gen.digit(1,2,
                                   lambda = Splus2.res$eig,
                                   tau = Splus2.res$tau,
                                   digit4tau = 2,
                                   axisName = "Component ")

stsU3.labels <- createxyLabels.gen.digit(1,2,
                                   lambda = Splus3.res$eig,
                                   tau = Splus3.res$tau,
                                   digit4tau = 2,
                                   axisName = "Component ")

## sts.constraints
stsU1.constraints <- minmaxHelper4Partial(Splus1.res$fi, Splus1.res$partial.fi)
stsU3.constraints <- minmaxHelper4Partial(Splus3.res$fi, Splus3.res$partial.fi)

## Splus.U1
### Fi
Splus1.fi <- createFactorMap(Splus1.res$fi,
                             constraints = stsU1.constraints,
                             text.cex = 1.5,
                             cex = 1.5,
                             col.points = "#42376B",
                             col.labels = "#42376B",
                             col.background = NULL,
                             col.axes = "#42376B",
                             width.axes = 0.5,
                             label.axisName = "Component ",
                             alpha.axes = 0.5,
                             title = "Subspace 1: global and partial\nfactor scores of products\n(rows)")
Splus1.pfi <- createPartialFactorScoresMap(Splus1.res$fi,
                                           Splus1.res$partial.fi,
                                           colors4Blocks = col.raters$gc[dimnames(Splus1.res$partial.fi)[3][[1]],],
                                           colors4Items = "#42376B",
                                           size.points = 0.5,
                                           alpha.lines = 0.7, alpha.points = 0.7, size.lines = 0.5)

srv.cmpfi1.plot <- Splus1.fi$zeMap_background + Splus1.fi$zeMap_text + Splus1.pfi$linesColByBlocks + Splus1.pfi$pointsColByBlocks + Splus1.fi$zeMap_dots + stsU1.labels + ggtheme

### Fj
Splus1.fj <- createFactorMap(Splus1.res$fj,
                             col.labels = col.flavors$oc[which1],
                             col.points = col.flavors$oc[which1],
                             text.cex = 2,
                             col.background = NULL,
                             col.axes = "#42376B",
                             width.axes = 0.5,
                             label.axisName = "Component ",
                             alpha.axes = 0.5,
                             title = "STATIS: column factor scores of \nraters and flavors (columns)")
srv.cmpfj1.plot <- Splus1.fj$zeMap + stsU1.labels + ggtheme

### Barplots - fj
srv.cmpfj1.dat <- data.frame(Splus1.res$fj[,c(1,2)], flavors = flavors[which1], raters = raters[which1], 
                             ID = rownames(Splus1.res$fj), IDnum = c(1:length(flavors[which1])))


(srv.cmpfj1.bar1 <- PrettyBarPlot.MuSu(srv.cmpfj1.dat$IDnum, srv.cmpfj1.dat$Dimension.1, grp = raters[which1], 
                                       color = drop(col.raters$gc[unique(srv.cmpfj1.dat$raters),]), 
                                       label.bar = flavors[which1],
                                       grp.label.size = 1.5,
                                       grp.label.position = 0.22,
                                       bar.label.size = 1.25,
                                       title = "Subspace 1:\ncolumn factor scores",
                                       ylab = "Component 1") + ylim(c(-0.2,0.22)) + 
    ggtheme.noborder)

(srv.cmpfj1.bar2 <- PrettyBarPlot.MuSu(srv.cmpfj1.dat$IDnum, srv.cmpfj1.dat$Dimension.2, grp = raters[which1], 
                                       color = drop(col.raters$gc[unique(srv.cmpfj1.dat$raters),]), 
                                       label.bar = flavors[which1],
                                       grp.label.size = 1.5,
                                       grp.label.position = 0.18,
                                       bar.label.size = 1.25,
                                       title = "",
                                       ylab = "Component 2") + ylim(c(-0.14,0.18)) + 
    ggtheme.noborder + theme(plot.title = element_blank()))
## Splus.U2
### Fi
Splus2.fi2plot <- Splus2.res$fi[,1:2]
Splus2.fi2plot[,1:2] <- Splus2.res$fi[,1:2]*-1
Splus2.pfi2plot <- Splus2.res$partial.fi
Splus2.pfi2plot[,1:2,] <- Splus2.res$partial.fi[,1:2,]*-1
stsU2.constraints <- minmaxHelper4Partial(Splus2.fi2plot, Splus2.pfi2plot)

Splus2.fi <- createFactorMap(Splus2.fi2plot,
                             constraints = stsU2.constraints,
                             text.cex = 1.5,
                             cex = 1.5,
                             col.points = "#42376B",
                             col.labels = "#42376B",
                             col.background = NULL,
                             col.axes = "#42376B",
                             width.axes = 0.5,
                             label.axisName = "Component ",
                             alpha.axes = 0.5,
                             title = "Subspace 2: global and partial\nfactor scores of products\n(rows)")
Splus2.pfi <- createPartialFactorScoresMap(Splus2.fi2plot,
                                           Splus2.pfi2plot,
                                           colors4Blocks = col.raters$gc[dimnames(Splus2.pfi2plot)[3][[1]],],
                                           colors4Items = "#42376B",
                                           size.points = 0.5,
                                           alpha.lines = 0.7, alpha.points = 0.7, size.lines = 0.5)

srv.cmpfi2.plot <- Splus2.fi$zeMap_background + Splus2.fi$zeMap_text + Splus2.pfi$linesColByBlocks + Splus2.pfi$pointsColByBlocks+ Splus2.fi$zeMap_dots + stsU2.labels + ggtheme

### Fj
Splus2.fj <- createFactorMap(Splus2.res$fj,
                             col.labels = col.flavors$oc[which2],
                             col.points = col.flavors$oc[which2],
                             text.cex = 2,
                             col.background = NULL,
                             col.axes = "#42376B",
                             width.axes = 0.5,
                             label.axisName = "Component ",
                             alpha.axes = 0.5,
                             title = "STATIS: column factor scores of \nraters and flavors (columns)")
srv.cmpfj2.plot <- Splus2.fj$zeMap + stsU2.labels + ggtheme

### Barplots - fj
srv.cmpfj2.dat <- data.frame(Splus2.res$fj[,c(1,2)]*-1, flavors = flavors[which2], raters = raters[which2], 
                             ID = rownames(Splus2.res$fj), IDnum = c(1:length(flavors[which2])))


(srv.cmpfj2.bar1 <- PrettyBarPlot.MuSu(srv.cmpfj2.dat$IDnum, srv.cmpfj2.dat$Dimension.1, grp = raters[which2], 
                                       color = drop(col.raters$gc[unique(srv.cmpfj2.dat$raters),]), 
                                       label.bar = flavors[which2],
                                       grp.label.size = 1.5,
                                       grp.label.position = 0.18,
                                       bar.label.size = 1.25,
                                       title = "Subspace 2:\ncolumn factor scores",
                                       ylab = "Component 1") + ylim(c(-0.16,0.18)) + 
    ggtheme.noborder)

(srv.cmpfj2.bar2 <- PrettyBarPlot.MuSu(srv.cmpfj2.dat$IDnum, srv.cmpfj2.dat$Dimension.2, grp = raters[which2], 
                                       color = drop(col.raters$gc[unique(srv.cmpfj2.dat$raters),]),
                                       label.bar = flavors[which2],
                                       grp.label.size = 1.5,
                                       grp.label.position = 0.15,
                                       bar.label.size = 1.25,
                                       title = "",
                                       ylab = "Component 2") + ylim(c(-0.11,0.15)) + 
    ggtheme.noborder + theme(plot.title = element_blank()))

## Splus.U3
### Fi
Splus3.fi <- createFactorMap(Splus3.res$fi,
                             constraints = stsU3.constraints,
                             text.cex = 1.5,
                             cex = 1.5,
                             col.points = "#42376B",
                             col.labels = "#42376B",
                             col.background = NULL,
                             col.axes = "#42376B",
                             width.axes = 0.5,
                             label.axisName = "Component ",
                             alpha.axes = 0.5,
                             title = "sSTATIS: global and partial\nfactor scores of products\n(rows)")
Splus3.pfi <- createPartialFactorScoresMap(Splus3.res$fi,
                                           Splus3.res$partial.fi,
                                           colors4Blocks = col.raters$gc[dimnames(Splus3.res$partial.fi)[3][[1]],],
                                           colors4Items = "#42376B",
                                           size.points = 0.5,
                                           alpha.lines = 0.7, alpha.points = 0.7, size.lines = 0.5)

srv.cmpfi3.plot <- Splus3.fi$zeMap_background + Splus3.pfi$linesColByBlocks + Splus3.pfi$pointsColByBlocks + Splus3.fi$zeMap_text + Splus3.fi$zeMap_dots + stsU3.labels + ggtheme

### Fj
Splus3.fj <- createFactorMap(Splus3.res$fj,
                             col.labels = col.flavors$oc[which3],
                             col.points = col.flavors$oc[which3],
                             text.cex = 2,
                             col.background = NULL,
                             col.axes = "#42376B",
                             width.axes = 0.5,
                             label.axisName = "Component ",
                             alpha.axes = 0.5,
                             title = "Subspace 3: column factor scores of \nraters and flavors (columns)")
srv.cmpfj3.plot <- Splus3.fj$zeMap + stsU3.labels + ggtheme

### Barplots - fj
srv.cmpfj3.dat <- data.frame(Splus3.res$fj[,c(1,2)], flavors = flavors[which3], raters = raters[which3], 
                             ID = rownames(Splus3.res$fj), IDnum = c(1:length(flavors[which3])))


(srv.cmpfj3.bar1 <- PrettyBarPlot.MuSu(srv.cmpfj3.dat$IDnum, srv.cmpfj3.dat$Dimension.1, grp = raters[which3], 
                                       color = drop(col.raters$gc[unique(srv.cmpfj3.dat$raters),]), 
                                       label.bar = flavors[which3],
                                       grp.label.size = 1.5,
                                       grp.label.position = 0.2,
                                       bar.label.size = 1.25,
                                       title = "Subspace 3:\ncolumn factor scores",
                                       ylab = "Component 1") + ylim(c(-0.18,0.2)) + 
    ggtheme.noborder)

(srv.cmpfj3.bar2 <- PrettyBarPlot.MuSu(srv.cmpfj3.dat$IDnum, srv.cmpfj3.dat$Dimension.2, grp = raters[which3], 
                                       color = drop(col.raters$gc[unique(srv.cmpfj3.dat$raters),]), 
                                       label.bar = flavors[which3],
                                       grp.label.size = 1.5,
                                       grp.label.position = 0.15,
                                       bar.label.size = 1.25,
                                       title = "",
                                       ylab = "Component 2") + ylim(c(-0.12,0.15)) + 
    ggtheme.noborder + theme(plot.title = element_blank()))


## -------------------------------------------------------------------
RV.ns.alpha <- res.sts.byraters$mexPosition.Data$InnerProduct$alphaWeights
names(RV.ns.alpha) <- sub(".","",names(RV.ns.alpha))
alpha.mat.cp <- recode(raters, !!!RV.ns.alpha)

# Sparsity index - by raters
### Run iterations
K <- 4L

I <- 8
J <- ncol(grand.tab)
J.raters <- length(unique(raters))

parz <- expand.grid(
  rdsleft = sqrt(I),#seq(sqrt(I), 1, length = 15),
  rdsright = seq(sqrt(J.raters), 1, length = 15))

parz <- parz[order(rowSums(parz^2/c(I,J.raters)), decreasing = T), ]

iter <- 1
res.ssts.list <- NULL
U0 <- NULL
V0 <- NULL

for (rds.iter in 1:nrow(parz)) {
  cat(sprintf("Left radius: %0.2f - Right radius: %0.2f\n", parz[rds.iter, 1], parz[rds.iter, 2]))
  res.ssts.list[[iter]] <- sparseGSVD(grand.tab.preproc, 
                                      LW = rep(1/nrow(grand.tab), nrow(grand.tab)), 
                                      RW = alpha.mat.cp,
                                      k = K,
                                      init = "svd", seed = 1000, 
                                      rdsLeft = rep(parz[rds.iter, 1], K), 
                                      rdsRight = rep(parz[rds.iter, 2], K),
                                      grpRight = raters.fac)
  # U0 <- res.ssts.list[[iter]]$p
  # V0 <- res.ssts.list[[iter]]$q
  iter <- iter + 1
}

### record results
dat.si.byraters <- data.frame(
  parz = parz, 
  SI = unname(t(data.frame(SI = sapply(res.ssts.list, extract_si, "SI", TRUE)))))

dat.si.m.byraters <- dat.si.byraters %>% tidyr::pivot_longer(starts_with("SI"), names_to = "k", names_prefix = "SI\\.")

dat.fit.zeros.byraters <- data.frame(
  rdsleft = dat.si.m.byraters$parz.rdsleft, 
  rdsright = dat.si.m.byraters$parz.rdsright, 
  k = dat.si.m.byraters$k,
  fit = pivot_the_tab(sapply(res.ssts.list, extract_si, "r1", TRUE)),
  zeros = pivot_the_tab(sapply(res.ssts.list, extract_si, "r4", TRUE)),
  SI = dat.si.m.byraters$value)

theta <- seq(pi, 3/2*pi, length.out = 150)


## ----fitzeroratioplot-----------------------------------------------
dimcol <- GetTolColors(n = 7, scheme = "discrete rainbow")
names(dimcol) <- as.factor(c(7:1))
# dimcol <- dimcol[c("2","3","4","5","6","7")]
dimcol["2"] <- "#f7c856"

dat.fit.zeros.byraters <- na.omit(dat.fit.zeros.byraters)
# dat.fit.zeros <- na.omit(dat.fit.zeros %>% filter(k %in% c("2", "3", "4", "5","6")))
dat.fit.zeros.byraters$max <- ifelse(dat.fit.zeros.byraters$k > 1 & dat.fit.zeros.byraters$SI == max(dat.fit.zeros.byraters$SI[which(dat.fit.zeros.byraters$k > 1)]), "MAX", "NOTMAX")
dat.fit.zeros.byraters$MAX <- ifelse(dat.fit.zeros.byraters$k > 1 & dat.fit.zeros.byraters$SI == max(dat.fit.zeros.byraters$SI[which(dat.fit.zeros.byraters$k > 1)]), "darkorchid4", "white")
dat.fit.zeros.byraters$alpha <- ifelse(dat.fit.zeros.byraters$k > 1 & dat.fit.zeros.byraters$SI == max(dat.fit.zeros.byraters$SI[which(dat.fit.zeros.byraters$k > 1)]), 1, 0.2)

siplot <- ggplot(dat.fit.zeros.byraters, aes(zeros, fit)) + 
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
  geom_point(aes(fill = k, alpha = I(alpha)), color = dat.fit.zeros.byraters$MAX, shape = 21, stroke = 1) +
  guides(color = guide_legend(override.aes = list(shape = 19, color = dimcol), nrow = 1, byrow = TRUE), shape = 19) +
  theme_bw() + 
  scale_color_manual(breaks= as.factor(c(1:7)), values = dimcol) +
  scale_fill_manual(breaks= as.factor(c(1:7)), values = dimcol) +
  coord_equal(xlim=0:1, ylim = 0:1) + 
  labs(x = "Zero ratio", y = "Fit", size = "Sparsity\nIndex", fill = "Number of\nComponents") + 
  # ggtitle("Ratio of zeros as a\nfunction of the fit ratio") + 
  theme(
    axis.title = element_text(size = 6, color = "#42376B"), 
    axis.text = element_text(size = 6, color = "#42376B"), 
    title = element_text(size = 6, color = "#42376B"),
    legend.key.size = unit(0.1, "cm"), legend.text = element_text(size = 6), legend.position = "right",
    legend.spacing.x = unit(0.05, "cm"),
    legend.margin = margin(0,0,0,0), legend.box.margin = margin(-5,0,0,0),
    legend.title.align = 0, legend.title = element_text(size = 6),
    panel.border = element_rect(size = 1, color = "#42376B", fill = NA),
    plot.margin = unit(c(0.1,0,0,0), unit = "cm"),
    panel.grid = element_blank()) + 
  with(dat.fit.zeros.byraters[dat.fit.zeros.byraters$max == "MAX",], annotate(geom = "point", x = zeros, y = fit, alpha = alpha, color = dimcol[k], size = 1)) +
  with(dat.fit.zeros.byraters[dat.fit.zeros.byraters$max == "MAX",], annotate(geom = "segment", x = zeros + 0.45, y = fit + 0.4, xend = zeros+ 0.02, yend = fit + 0.02, arrow = arrow(length = unit(0.05, "inches"), type = "closed"), color = "darkorchid4", size = 0.6)) + #0.05, 0.5
  # with(dat.fit.zeros.byraters[dat.fit.zeros.byraters$max == "MAX",], annotate(geom = "segment", x = zeros + 0.025, y = fit + 0.022, xend = zeros+ 0.02, yend = fit + 0.02, arrow = arrow(length = unit(0.04, "inches"), type = "closed"), color = dimcol[k], size = 0.01)) + #0.04, 0.01
  with(dat.fit.zeros.byraters[dat.fit.zeros.byraters$max == "MAX",], annotate(geom = "label", x = zeros + 0.46, y = fit + 0.38, color = "darkorchid4", label = substring(sprintf("%.3f", SI), 2), fill = dimcol[k], size = 2, label.padding = unit(0.04, "in")))

siplot


## -------------------------------------------------------------------
optres <- dat.fit.zeros.byraters[dat.fit.zeros.byraters$max == "MAX",]
kopt <- optres$k
rdsleftopt <- optres$rdsleft
rdsrightopt <- optres$rdsright
SIopt <- optres$SI

kopt;rdsleftopt/sqrt(I);rdsrightopt/sqrt(J);SIopt


## -------------------------------------------------------------------
## run opt sparseSTATIS
ssts.svd<- sparseGSVD(grand.tab.preproc, 
                      LW = rep(1/nrow(grand.tab), nrow(grand.tab)), 
                      RW = alpha.mat.cp,
                      k = kopt,
                      init = "svd", seed = 1000, 
                      rdsLeft = rep(rdsleftopt, kopt), 
                      rdsRight = rep(rdsrightopt, kopt),
                      grpRight = raters.fac)
colnames(ssts.svd$u) <- colnames(ssts.svd$v) <- colnames(ssts.svd$p) <- colnames(ssts.svd$q) <- paste0("Dimension ", c(1:ncol(ssts.svd$u)))
ssts.res <- STATIS.Sout(grand.tab.preproc, ssts.svd, rep(1/nrow(grand.tab), nrow(grand.tab)), alpha.mat.cp, raters.fac)


## -------------------------------------------------------------------
eigres.scmp <- data.frame(eig = ssts.res$eig,
                          tau = ssts.res$tau)

(scmp.scree <- PlotMyScree(eigres.scmp, cex = 1.3, lwd = 0.5, title = "Scree plot for the sparsified grand table from sSTATIS", text.cex = 6) + theme(axis.title = element_text(size = 6, color = "#42376B"), axis.text = element_text(size = 6, color = "#42376B"), plot.title = element_text(face = "plain"),  title = element_text(size = 6, color = "#42376B"), panel.border = element_rect(size = 1, color = "#42376B", fill = NA), axis.text.y = element_text(size = 4)))


## -------------------------------------------------------------------
options(ggrepel.max.overlap = Inf)

ssts.labels <- createxyLabels.gen.digit(1,2,
                                  lambda = ssts.res$eig,
                                  tau = ssts.res$tau,
                                  digit4tau = 2,
                                  axisName = "Component ")

## sts.constraints
ssts.constraints <- minmaxHelper4Partial(ssts.res$fi, ssts.res$partial.fi)

### Fi
ssts.fi <- createFactorMap(ssts.res$fi,
                           constraints = ssts.constraints,
                           text.cex = 1.5,
                           cex = 1.5,
                           col.points = "#42376B",
                           col.labels = "#42376B",
                           col.background = NULL,
                           col.axes = "#42376B",
                           width.axes = 0.5,
                           label.axisName = "Component ",
                           alpha.axes = 0.5,
                           title = "sSTATIS: global and partial\nfactor scores of products (rows)")
ssts.pfi <- createPartialFactorScoresMap(ssts.res$fi,
                                         ssts.res$partial.fi,
                                         colors4Blocks = col.raters$gc[dimnames(ssts.res$partial.fi)[3][[1]],],
                                         colors4Items = "#42376B", size.points = 0.5,
                                          alpha.lines = 0.5, alpha.points = 0.7, size.lines = 0.5)

ssts.cmpfi.plot <- ssts.fi$zeMap_background + ssts.pfi$pointsColByBlocks + ssts.pfi$linesColByBlocks + ssts.fi$zeMap_dots + ssts.fi$zeMap_text + ssts.labels + ggtheme

### Fj
ssts.fj <- createFactorMap(ssts.res$fj,
                           col.labels = col.flavors$oc,
                           col.points = col.flavors$oc,
                           text.cex = 2,
                           col.background = NULL,
                           col.axes = "#42376B",
                           width.axes = 0.5,
                           label.axisName = "Component ",
                           alpha.axes = 0.5,
                           title = "sSTATIS: column factor scores of \nraters and flavors (columns)")
ssts.cmpfj.plot <- ssts.fj$zeMap + ssts.labels + ggtheme

### Barplots - fj
ssts.cmpfj.dat <- data.frame(ssts.res$fj[,c(1,2)], flavors = flavors, raters = raters, 
                             ID = rownames(ssts.res$fj), IDnum = c(1:length(flavors)))


(ssts.cmpfj.bar1 <- PrettyBarPlot.MuSu(ssts.cmpfj.dat$IDnum, ssts.cmpfj.dat$Dimension.1, grp = raters, 
                                       color = drop(col.raters$gc[unique(ssts.cmpfj.dat$raters),]), 
                                       label.bar = flavors,
                                       grp.label.size = 1.5,
                                       grp.label.position = 0.2,
                                       bar.label.size = 1.25,
                                       title = "sSTATIS: column factor scores",
                                       ylab = "Component 1") + ylim(c(-0.13,0.2)) + 
    ggtheme.noborder)
(ssts.cmpfj.bar2 <- PrettyBarPlot.MuSu(ssts.cmpfj.dat$IDnum, ssts.cmpfj.dat$Dimension.2, grp = raters, 
                                       color = drop(col.raters$gc[unique(ssts.cmpfj.dat$raters),]), 
                                       label.bar = flavors,
                                       grp.label.size = 1.5,
                                       grp.label.position = 0.2,
                                       bar.label.size = 1.25,
                                       title = "sSTATIS: column factor scores",
                                       ylab = "Component 2") + ylim(c(-0.13,0.2)) + 
    ggtheme.noborder)



## -------------------------------------------------------------------
# png(filename="Figure1-SSTATIS.png", width = 140, height = 200, units = "mm", bg = "white",res = 800)
# gridExtra::grid.arrange(grobs = list(rv.sts.scree, rv.fi.plot, rv.alpha,
#                                      siplot, scmp.scree,
#                                      sts.cmpfi.plot, ssts.cmpfi.plot,
#                                      sts.cmpfj.bar1, ssts.cmpfj.bar1,
#                                      sts.cmpfj.bar2, ssts.cmpfj.bar2),
#                         widths = c(0.165, 0.165, 0.165, 0.165, 0.165, 0.165),
#                         heights = c(0.18, 0.16,0.29,0.19, 0.18),
#                         layout_matrix = rbind(c(1,1,2,2,3,3),
#                                               c(4,4,4,5,5,5),
#                                               c(6,6,6,7,7,7),
#                                               c(8,8,8,9,9,9),
#                                               c(10,10,10,11,11,11))
# )
# dev.off()

png(filename="Figure1-SSTATIS(change).png", width = 140, height = 200, units = "mm", bg = "white",res = 800)
gridExtra::grid.arrange(grobs = list(rv.sts.scree, rv.fi.plot, rv.alpha,
                                     siplot, scmp.scree,
                                     sts.cmpfi.plot, ssts.cmpfi.plot,
                                     sts.cmpfj.bar1, ssts.cmpfj.bar1,
                                     sts.cmpfj.bar2, ssts.cmpfj.bar2),
                        widths = c(0.165, 0.166, 0.164, 0.165, 0.165, 0.165),
                        heights = c(0.18, 0.16, 0.165, 0.165, 0.165, 0.165),
                        layout_matrix = rbind(c(1,1,2,2,3,3),
                                              c(4,4,4,5,5,5),
                                              c(6,6,8,8,8,8),
                                              c(6,6,10,10,10,10),
                                              c(7,7,9,9,9,9),
                                              c(7,7,11,11,11,11))
)
dev.off()


## -------------------------------------------------------------------
png(filename="Figure2-SSTATIS.png", width = 140, height = 200, units = "mm", bg = "white",res = 800)
gridExtra::grid.arrange(grobs = list(rv.ssts.scree, as.grob(rv.heat),
                                     siplot1, siplot2, sprv.loli, sprv.cmp.scree.adj,
                                     srv.cmpfi1.plot, srv.cmpfi2.plot, srv.cmpfi3.plot,
                                     srv.cmpfj1.bar1, srv.cmpfj2.bar1, srv.cmpfj3.bar1,
                                     srv.cmpfj1.bar2, srv.cmpfj2.bar2, srv.cmpfj3.bar2),
                        widths = c(0.165, 0.165, 0.165, 0.155, 0.175, 0.165),
                        heights = c(0.03, 0.06, 0.06, 0.06, 0.06, 0.06, 0.04, 0.23,0.15, 0.15),
                        layout_matrix = rbind(c(NA,NA,1,1,2,2),
                                              c(3,3,1,1,2,2),
                                              c(3,3,1,1,2,2),
                                              c(3,3,5,5,6,6),
                                              c(4,4,5,5,6,6),
                                              c(4,4,5,5,6,6),
                                              c(4,4,5,5,6,6),
                                              c(7,7,8,8,9,9),
                                              c(10,10,11,11,12,12),
                                              c(13,13,14,14,15,15))
)
dev.off()

