## ----setup, include = FALSE-----------------------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
# detach("package:SPAFAC", unload = TRUE)
# devtools::install_github("juchiyu/SPAFAC")

library(PTCA4CATA)
library(DistatisR)
library(ExPosition)
library(ggplot2)
library(tidyr)
library(dplyr)
library(pheatmap)

devtools::load_all("~/git/sGSVD/")
devtools::load_all("~/git/SPAFAC/")

tool.path <- "functions/"

source(paste0(tool.path, "PlotMyScree.R"))
source(paste0(tool.path, "DiStatis.preproc.R"))
source(paste0(tool.path, "PlotFactor.R"))
source(paste0(tool.path, "createxyLabels.gen.digit.R"))
source(paste0(tool.path, "createLabel.gen.digit.R"))
source(paste0(tool.path, "plot.smca.results.R"))
source(paste0(tool.path, "PrettyBarPlot.MuSu.R"))

# other functions
extract_si <- function(sdisca, si = "SI", fromSGSVD = FALSE) {
  if ("list" %in% class(sdisca)) {
    if (fromSGSVD) {
      return(sdisca$SI[[si]])
    }else {
      return(sdisca$sparsity$SI[[si]])
    }
  } else if (is.na(sdisca)) {
    return(NA)
  } else {
    stop("Something is a foot!")
  }
}

pivot_the_tab <- function(dat) {
  df <- data.frame(V = unname(t(data.frame(dat))))
  df %>% 
    tidyr::pivot_longer(starts_with("V"), names_to = "k", names_prefix = "V\\.") %>%
    select(value) %>% purrr::as_vector() %>% unname()
}



## Our favorite ggplot themes -----------------
ggtheme <- theme(
  axis.title = element_text(size = 6, color = "#42376B"),
  axis.text = element_text(size = 6, color = "#42376B"),
  title = element_text(size = 6, color = "#42376B"),
  panel.border = element_rect(size = 1.5, color = "#42376B", fill = NA)
)

ggtheme.noborder <- theme(
    axis.title = element_text(size = 6, color = "#42376B"),
    axis.text = element_text(size = 6, color = "#42376B"),
    title = element_text(size = 6, color = "#42376B")
  )

options(ggrepel.max.overlap = Inf)



## ----load data------------------------------------------------------
data("sortingWines", package = "DistatisR")

wine.dx <- sortingWines$winesDescription
wine.col <- wine.dx$ColorCode

free.exp <- sortingWines$freeSortExperts
colnames(free.exp) <- paste0(colnames(sortingWines$freeSortExperts), ".freEx")
tern.exp <- sortingWines$ternarySortExperts
colnames(tern.exp) <- paste0(colnames(sortingWines$ternarySortExperts), ".terEx")
tern.nov <- sortingWines$ternarySortNovices
colnames(tern.nov) <- paste0(colnames(sortingWines$ternarySortNovices), ".terNv")

gtab <- cbind(free.exp, tern.exp, tern.nov)
exp.dx <- sub('.*\\.', "", colnames(gtab))

exp.col <- createColorVectorsByDesign(makeNominalData(as.matrix(exp.dx)))
rownames(exp.col$gc) <- sub('.*\\.', "", rownames(exp.col$gc))

DisCube <- DistanceFromSort(gtab)

nbgrp <- factor(ifelse(unname(apply(gtab, 2, max) <= 3), "LessThan3", "MoreThan3"))
nbgrp.col <- c("darkorchid4", "darkolivegreen")[nbgrp]


## ----discube-preproc------------------------------------------------
DisCube.proc <- apply(DisCube, 3, DiStatis.preproc, simplify = FALSE)


## ----distatis-------------------------------------------------------
distatis.res <- distatis(DisCube)
G12 <- data.frame(distatis.res$res4Cmat$G[, 1:2])
# res.hc <- hclust(dist(G12), method = "ward.D2")
res.hc <- hclust(dist(distatis.res$res4Cmat$C), method = "ward.D2")
cl <- paste0("cluster", cutree(res.hc, 3))


## ----preparation----------------------------------------------------
RVmat <- distatis.res$res4Cmat$C
rv.eig <- eigen(RVmat)


## ----fig1a rv.scree-------------------------------------------------
rv.eigres <- data.frame(
  eig = distatis.res$res4Cmat$eigValues,
  tau = distatis.res$res4Cmat$tau)

rv.scree <- PlotMyScree(
  rv.eigres, lwd = 0.25, cex = 0.5, text.cex = 6) + 
  ggtitle("Scree plot for non-sparse\nRV matrix from DiSTATIS") + 
  ggtheme + 
  theme(plot.title = element_text(face = "plain"),  
        panel.border = element_rect(size = 1, color = "#42376B", fill = NA))

(fig1a <- rv.scree + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))


## ----fig1b compromise.scree-----------------------------------------
cmp.eigres <- data.frame(
  eig = distatis.res$res4Splus$eigValues,
  tau = distatis.res$res4Splus$tau)

cmp.scree <- PlotMyScree(
  cmp.eigres, lwd = 0.25, cex = 0.5, text.cex = 6) + 
  ggtitle("Scree plot for non-sparse\nCompromise matrix from DiSTATIS") + 
  ggtheme + 
  theme(plot.title = element_text(face = "plain"),  
        panel.border = element_rect(size = 1, color = "#42376B", fill = NA))

(fig1b <- cmp.scree + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))


## ----fig1c cmp.sparsity, cache = TRUE-------------------------------
I.cmp <- nrow(distatis.res$res4Splus$Splus)
K <- 3L

CMPmat <- distatis.res$res4Splus$Splus

cmp.U0 <- distatis.res$res4Splus$eigVectors[, 1:K]
res.cmp.speig <- NULL
cmp.params <- seq(sqrt(I.cmp), 1, length = 20)

for (i in seq_along(cmp.params)) {
  # cat(sprintf("Radius: %0.2f\n", cmp.params[i]))
  res.cmp.speig[[i]] <- tryCatch(
    sparseEIGEN(X = CMPmat, k = K, init = cmp.U0, rds = rep(cmp.params[i], K)),
    error = function(e) NA)
  if (!is.na(res.cmp.speig[[i]])){
    U0 <- res.cmp.speig[[i]]$vectors
  }
}

res.cmp.speig.sparseindex <- lapply(res.cmp.speig, sparseIndexEigen, eigenValues = distatis.res$res4Splus$eigValues)


## -------------------------------------------------------------------
cmp.si.mat.tmp <- sapply(res.cmp.speig.sparseindex,
                     \(.) c(SI = .$SI,
                            fitRatio = .$fitRatio,
                            zeroRatio = .$zeroRatio))

cmp.si.dat <-  data.frame(
  radius = cmp.params,
  t(cmp.si.mat.tmp)) %>% 
  pivot_longer(
    -radius, 
    names_to = c(".value", "k"),
    names_pattern = "([a-zA-Z]+)([0-9]+)") %>%
  drop_na() %>%
  mutate(
    max = ifelse(k > 1 & radius > 1 & SI == max(SI[k > 1 & radius > 1 ]), "MAX", "NOTMAX"),
    MAX = ifelse(k > 1 & radius > 1  & SI == max(SI[k > 1 & radius > 1 ]), "darkorchid4", "white"),
    alpha = ifelse(k > 1 & radius > 1 & SI == max(SI[k > 1 & radius > 1 ]), 1, 0.2))

theta <- seq(pi, 3/2*pi, length.out = 150)

dimcol <- inlmisc::GetColors(n = 7, scheme = "discrete rainbow")
names(dimcol) <- c(7:1)
dimcol["2"] <- "#F7C856"


cmp.siplot <- cmp.si.dat %>% ggplot(aes(zeroRatio, fitRatio)) + 
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
  geom_point(aes(fill = k, alpha = I(alpha)), color = cmp.si.dat$MAX, shape = 21, stroke = 1) +
  # geom_line(aes(color = k)) +
  guides(color = guide_legend(override.aes = list(shape = 19, color = dimcol), nrow = 1, byrow = TRUE), shape = 19) +
  theme_bw() + 
  scale_color_manual(breaks = as.factor(c(1:K)), values = dimcol) +
  scale_fill_manual(breaks = as.factor(c(1:K)), values = dimcol) +
  coord_equal(xlim=0:1, ylim = 0:1) + 
  labs(x = "Zero ratio", y = "Fit", size = "Sparsity\nIndex", fill = "Number of\nComponents") + 
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
  with(cmp.si.dat[cmp.si.dat$max == "MAX",], annotate(geom = "point", x = zeroRatio, y = fitRatio, alpha = alpha, color = dimcol[k], size = 1)) +
  with(cmp.si.dat[cmp.si.dat$max == "MAX",], annotate(geom = "segment", x = zeroRatio - 0.5, y = fitRatio - 0.5, xend = zeroRatio-0.05, yend = fitRatio-0.05, arrow = arrow(length = unit(0.05, "inches"), type = "closed"), color = "darkorchid4", size = 0.6)) + 
  with(cmp.si.dat[cmp.si.dat$max == "MAX",], annotate(geom = "label", x = zeroRatio - 0.5, y = fitRatio - 0.5, color = "darkorchid4", label = substring(sprintf("%.3f", SI), 2), fill = dimcol[k], size = 2, label.padding = unit(0.04, "in")))

(fig1c <- cmp.siplot)


## ----fig2b heatmap--------------------------------------------------
distrow <- dist(G12[, 1:2])
distcol <- distrow

annot.groups <-  data.frame(Cluster = cl, row.names = rownames(RVmat))

fig2b <- pheatmap(RVmat, show_rownames = FALSE, show_colnames = FALSE, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2", color = colorRampPalette(c("white", "darkorchid3"))(50), gaps_col = cl, gaps_row = cl, annotation_row = annot.groups, annotation_col = annot.groups, , annotation_legend = FALSE, fontsize = 5,
                    treeheight_col = 7, treeheight_row = 7, angle_col = "90", 
                    legend.width = 5, cellwidth = 2, cellheight = 2, annotation_colors = list(Cluster = c(cluster1 = "#F8766D", cluster2 = "#00BA38", cluster3 = "#619CFF")))


fig2b$gtable$grobs[[8]]$children[[1]]$width <- fig2b$gtable$grobs[[8]]$children[[1]]$width/2
fig2b$gtable$grobs[[8]]$children[[1]]$height <- fig2b$gtable$grobs[[8]]$children[[1]]$height*0.6
fig2b$gtable$grobs[[8]]$children[[1]]$y <- fig2b$gtable$grobs[[8]]$children[[1]]$y*0.6 + unit(0.3, units = "npc")
fig2b$gtable$grobs[[8]]$children[[2]]$y <- fig2b$gtable$grobs[[8]]$children[[2]]$y*0.6 + unit(0.3, units = "npc")
fig2b$gtable$grobs[[8]]$children[[1]]$x <- fig2b$gtable$grobs[[8]]$children[[1]]$x #- unit(0.3, units = "npc")
fig2b$gtable$grobs[[8]]$children[[2]]$x <- fig2b$gtable$grobs[[8]]$children[[2]]$x - unit(0.3, units = "npc")



## ----fig1d rv.factorscores------------------------------------------
col4clusters <- c(cluster1 = "#F8766D", cluster2 = "#00BA38", cluster3 = "#619CFF")

rv.fi <- distatis.res$res4Cmat$G

## labels
rv.labels <- createxyLabels.gen(
  1, 2,
  lambda = distatis.res$res4Cmat$eigValues,
  tau = distatis.res$res4Cmat$tau,
  axisName = "Component ")

## constraints
colnames(rv.fi) <- paste0("Dimension ", c(1:ncol(rv.fi)))
rv.constraints <- minmaxHelper(rv.fi)

## Fi
RVmap <- createFactorMap(
  rv.fi,
  constraints = rv.constraints,
  text.cex = 1.5,
  cex = 1.5,
  col.background = NULL,
  col.axes = "#42376B",
  col.points = col4clusters[cl],#exp.col$oc,
  # col.labels = col4clusters[cl],#exp.col$oc,
  width.axes = 0.5,
  label.axisName = "Component ",
  alpha.axes = 0.5,
  title = "DiSTATIS: factor scores for RV"
)

(fig1d <- RVmap$zeMap_background + RVmap$zeMap_dots + rv.labels + ggtheme)


## ----fig1e alphas---------------------------------------------------
## plot barplot for alpha values
rv.alpha2plot <- distatis.res$res4Cmat$alpha  %>% as.matrix
rownames(rv.alpha2plot) <- sub(".", "", rownames(rv.alpha2plot))
colnames(rv.alpha2plot) <- "Weights for each tables (\alpha)"
rv.alpha <- PrettyBarPlot2(rv.alpha2plot[order(cl),1], threshold = 0, 
                           color4bar = col4clusters[sort(cl)],#exp.col$oc, 
                           font.size = 2, ) +
  ggtitle("Table weights") +
  ylab("Weights") + ggtheme.noborder

(fig1e <- rv.alpha + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line = element_line(color = "#42376B"))   )


## ----fig1f compromise.factorscores----------------------------------
cmp.fi <- distatis.res$res4Splus$F

## labels
cmp.labels <- createxyLabels.gen(
  1, 2,
  lambda = distatis.res$res4Splus$eigValues,
  tau = distatis.res$res4Splus$tau,
  axisName = "Component ")

## constraints
colnames(cmp.fi) <- paste0("Dimension ", c(1:ncol(cmp.fi)))
cmp.constraints <- minmaxHelper(cmp.fi)

## Fi
CMPmap <- createFactorMap(
  cmp.fi,
  constraints = cmp.constraints,
  text.cex = 1.5,
  cex = 1.5,
  col.background = NULL,
  col.axes = "#42376B",
  col.points = wine.col,
  col.labels = wine.col,
  width.axes = 0.5,
  label.axisName = "Component ",
  alpha.axes = 0.5,
  title = "DiSTATIS: factor scores for compromise"
)

(fig1f <- CMPmap$zeMap_background + CMPmap$zeMap_dots + CMPmap$zeMap_text + cmp.labels + ggtheme)



## ----fig1g compromise.factorscores----------------------------------
rds.opt.cmp <- cmp.si.dat$radius[cmp.si.dat$max == "MAX"]
kopt.cmp <- as.integer(max(as.integer(cmp.si.dat$k[cmp.si.dat$max == "MAX"]), 2L))

U0 <- distatis.res$res4Splus$eigVectors[, 1:kopt.cmp]
cmp.speig <- sparseEIGEN(X = CMPmat, k = kopt.cmp, init = U0, rds = rep(rds.opt.cmp, kopt.cmp))

dat.lol.cmp.sp <- cmp.speig$vectors[, 1:kopt.cmp] %>%
  data.frame(wine.dx) %>% 
  arrange(Wines) %>% 
  mutate(i = 1:nrow(.)) %>%
  relocate(i) %>%
  pivot_longer(starts_with("X")) %>%
  mutate(name = gsub("X", "Dimension ", name))

g.lol.cmp.sp <- dat.lol.cmp.sp %>%
  ggplot(aes(i, value, color = Type)) + 
  geom_hline(yintercept = 0, color = "grey") + 
  geom_point() + 
  geom_segment(aes(xend = i, y = 0, yend = value)) + 
  facet_grid(name~., scale = "free") + 
  scale_color_manual(values = setNames(unique(wine.dx$ColorCode), c("W", "R", "P")), labels = c("White", "Red", "Rosé")) +
  labs(x = "", y = "", color = "") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        axis.text = element_text(size = 6, color = "#42376B"),
        axis.title = element_text(size = 6, color = "#42376B"),
        strip.text = element_text(size = 6, color = "#42376B"),
        plot.margin = unit(c(0.1,0,0,0), unit = "cm"),
        strip.background = element_rect(color = "transparent", fill = "transparent", size = 0.1), 
        legend.position = "bottom")

(fig1g <- g.lol.cmp.sp + ggtheme.noborder)


## ----fig1h compromise.sparse.factorscores---------------------------
## Optimal factor scores
cmp.sp.fi <- cmp.speig$vectors %*% diag(sqrt(cmp.speig$values))
rownames(cmp.sp.fi) <- wine.dx$Wines

## labels
cmp.sp.labels <- createxyLabels.gen.digit(1,2,
                                  lambda = cmp.speig$values,
                                  tau = 100 * cmp.speig$values / sum(cmp.eigres$eig),
                                  digit4tau = 2,
                                  axisName = "Component ")

## sts.constraints
cmp.sp.constraints <- minmaxHelper(cmp.sp.fi)

### Fi
cmp.sp.fmap <- createFactorMap(cmp.sp.fi,
                           constraints = cmp.sp.constraints,
                           text.cex = 1.5,
                           cex = 1.5,
                           col.points = wine.col,
                           col.labels = wine.col,
                           col.background = NULL,
                           col.axes = "#42376B",
                           width.axes = 0.5,
                           label.axisName = "Component ",
                           alpha.axes = 0.5,
                           title = "sDiSTATIS: sparse factor\nscores of Compromise")


cmp.sp.fi.plot <- cmp.sp.fmap$zeMap_background + cmp.sp.fmap$zeMap_dots + cmp.sp.fmap$zeMap_text + cmp.sp.labels + ggtheme

(fig1h <- cmp.sp.fi.plot)


## ----figure 1, include = FALSE--------------------------------------
png(filename = "Figure1-SDiSTATIS.png", width = 140, height = 200, units = "mm", bg = "white", res = 600)
gridExtra::grid.arrange(
  grobs = list(fig1a, fig1b, fig1c,
               ggplotify::as.grob(fig2b), fig1d, fig1e,
               # fig1f, fig1g, fig1h),
               fig1f, fig1h),
  widths = c(1, 1, 1, 1, 1, 1),
  heights = c(1, 1, 1, 1),
  layout_matrix = rbind(c(1,1,2,2,3,3),
                        c(4,4,4,5,5,5),
                        c(4,4,4,6,6,6),
                        c(7,7,7,8,8,8)))
dev.off()


## ----first decomposition, cache = TRUE------------------------------
I.c <- nrow(distatis.res$res4Cmat$C)


K <- 2L
U0 <- abs(rv.eig$vectors[, 1:K])
res.speig <- NULL
params <- seq(sqrt(I.c), 1, length = 25)

for (i in seq_along(params)) {
  # cat(sprintf("Radius: %0.2f\n", params[i]))
  res.speig[[i]] <- tryCatch(
    sparseEIGEN(X = RVmat, k = K, init = U0, rds = rep(params[i], K)),
    error = function(e) NA)
  if (!is.na(res.speig[[i]])){
    U0 <- res.speig[[i]]$vectors
  }
}


res.speig.sparseindex <- lapply(res.speig[-(22:23)], sparseIndexEigen, eigenValues = rv.eig$values)


## ----optimum 1------------------------------------------------------

si.mat.tmp <- sapply(res.speig.sparseindex, 
                     \(.) c(SI = .$SI[1], 
                            fitRatio = .$fitRatio[1],
                            zeroRatio = .$zeroRatio[1])) 
si.dat <- data.frame(
  radius = params[-(22:23)],
  t(si.mat.tmp)
)

si.dat$max <- "FALSE"
si.dat$max[which.max(si.dat$SI)] <- "MAX"

theta <- seq(pi, 3/2*pi, length.out = 150)


rds.opt1 <- si.dat$radius[which.max(si.dat$SI)]
  
U0 <- abs(rv.eig$vectors[, 1:K])
rv.speig1 <- sparseEIGEN(X = RVmat, k = K, init = U0, rds = rep(rds.opt1, K))


## ----second decomposition, cache = TRUE-----------------------------
torem <- which(rv.speig1$vectors[, 1] != 0)
RVmat2 <- RVmat[-torem, -torem]
I.c2 <- nrow(RVmat2)
rv.eig2 <- eigen(RVmat2)

U0 <- abs(rv.eig2$vectors[, 1:K])
res.speig2 <- NULL
params2 <- seq(sqrt(I.c2), 1, length = 25)

for (i in seq_along(params2)) {
  # cat(sprintf("Radius: %0.2f\n", params2[i]))
  res.speig2[[i]] <- tryCatch(
    sparseEIGEN(X = RVmat2, k = K, init = U0, rds = rep(params2[i], K)),
    error = function(e) NA)
  if (!is.na(res.speig2[[i]])){
    U0 <- res.speig2[[i]]$vectors
  }
}


res.speig.sparseindex2 <- lapply(res.speig2, sparseIndexEigen, eigenValues = rv.eig2$values)


## ----optimum 2------------------------------------------------------
si.mat.tmp2 <- sapply(res.speig.sparseindex2, 
                     \(.) c(SI = .$SI[1], 
                            fitRatio = .$fitRatio[1],
                            zeroRatio = .$zeroRatio[1])) 
si.dat2 <- data.frame(
  radius = params2,
  t(si.mat.tmp2)
)

si.dat2$max <- "FALSE"
si.dat2$max[which.max(si.dat2$SI)] <- "MAX"

rds.opt2 <- si.dat2$radius[which.max(si.dat2$SI)]
  
U0 <- abs(rv.eig2$vectors[, 1:K])
rv.speig2 <- sparseEIGEN(X = RVmat2, k = K, init = U0, rds = rep(rds.opt2, K))


## ----last decomposition---------------------------------------------
torem2 <- which(rv.speig2$vectors[, 1] != 0)
RVmat3 <- RVmat2[-torem2, -torem2]
rv.eig3 <- eigen(RVmat3)


## ----combine all three sparse eigenvectors--------------------------
U1 <- rv.speig1$vectors[, 1]
U3 <- U2 <- rep(0, I.c)
U2[-torem] <- rv.speig2$vectors[, 1]
U3[-torem][-torem2] <- abs(rv.eig3$vectors[, 1])

Usp <- data.frame(U1 = U1, U2 = U2, U3 = U3)
dat.U <- data.frame(
  ID = rownames(Usp),
  IDnum = as.numeric(rownames(Usp)),
  Usp, exp.dx)


## ----create alphas--------------------------------------------------
DisCube.proc.arr <- simplify2array(DisCube.proc)

alphas <- data.frame(lapply(Usp, \(.) ./sum(.)))

Splus1 <- ComputeSplus(DisCube.proc.arr, alphas[, 1], isCube = TRUE)
Splus2 <- ComputeSplus(DisCube.proc.arr, alphas[, 2], isCube = TRUE)
Splus3 <- ComputeSplus(DisCube.proc.arr, alphas[, 3], isCube = TRUE)

eigen.Splus1 <- eigen(Splus1, symmetric = TRUE)
eigen.Splus2 <- eigen(Splus2, symmetric = TRUE)
eigen.Splus3 <- eigen(Splus3, symmetric = TRUE)


## ----fig2a----------------------------------------------------------
## Scree plot
eig.srv <- c(sum(rv.eig$values), sum(rv.eig2$values), sum(rv.eig3$values))
eigres.srv <- data.frame(eig = eig.srv, tau = eig.srv/sum(eig.srv))

(fig2a <- PlotMyScree(eigres.srv, lwd = 0.5, cex = 1.3, title = "Scree plot for the sparsified\nRV space from sDiSTATIS", text.cex = 5) + theme(axis.title = element_text(size = 5, color = "#42376B"), axis.text = element_text(size = 5, color = "#42376B"), title = element_text(face = "plain", size = 5, color = "#42376B"), plot.title = element_text(face = "plain"), panel.border = element_rect(size = 0.5, color = "#42376B", fill = NA)))


## ----fig2c----------------------------------------------------------
screedata1 <- data.frame(
  eig = eigen.Splus1$values,
  tau = 100 * eigen.Splus1$values / sum(eigen.Splus1$values),
  Subspace = "Subspace 1",
  components = seq_along(eigen.Splus1$values))
screedata2 <- data.frame(
  eig = eigen.Splus2$values,
  tau = 100 * eigen.Splus2$values / sum(eigen.Splus2$values),
  Subspace = "Subspace 2",
  components = seq_along(eigen.Splus2$values))
screedata3 <- data.frame(
  eig = eigen.Splus3$values,
  tau = 100 * eigen.Splus3$values / sum(eigen.Splus3$values),
  Subspace = "Subspace 3",
  components = seq_along(eigen.Splus3$values))

screedata <- do.call("rbind", list(screedata1, screedata2, screedata3))

screeplot <- screedata %>% 
  ggplot(aes(x = components, y = tau)) +
      geom_line(color = "grey40", size = 0.5) +
      geom_point(color = "#42376B", size = 1.1) +
    facet_wrap(~Subspace, ncol = 1) +
    scale_y_continuous(name = bquote(atop(bold(.(title)),paste('\n\n    Percentage of \nvariance explained (%)'))),
                       sec.axis = sec_axis(~.*(screedata$eig[1]/screedata$tau[1]), name = "Pseudo-eigenvalues")) +
    xlab("Components") +
    scale_x_continuous(breaks=unique(screedata$components)) +
    theme(text = element_text(size = 6, color = "#42376B"),
          legend.position = "none",
          strip.text = element_text(face = "plain", color = "#42376B", size = 6),
          axis.text.y.left = element_text(face = "plain", color = "#42376B", angle = 90, hjust = 0.5, size = 4),
          axis.text.y.right = element_text(face = "plain", color = "#42376B", angle = 270, hjust = 0.5, size = 3),
          axis.text.x = element_text(size = 6, color = "#42376B"),
          panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0.1,0,0,0), unit = "cm"),
          panel.border = element_rect(color = "#42376B", fill = "transparent"),
          strip.background = element_rect(fill = "transparent", size = 0.05), 
          panel.spacing =  unit(0, "cm"))



(fig2c <- screeplot)


## ----fig2d----------------------------------------------------------
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
  geom_point(color = "darkorchid4") + geom_line() +
  theme_bw() + scale_size_continuous(limits = 0:1) + 
  coord_equal(xlim=0:1, ylim = 0:1) + 
  labs(x = "Zero ratio", y = "Fit ratio", size = "Sparsity\nIndex") + 
  theme(
    panel.border = element_rect(size = 1.5, color = "#42376B", fill = NA),
    panel.grid = element_blank())  +
  with(si.dat[si.dat$max == "MAX",], annotate(geom = "point", x = zeroRatio, y = fitRatio, color = "#7F6CC7", size = 1)) +
  with(si.dat[si.dat$max == "MAX",], annotate(geom = "segment", x = zeroRatio + 0.25, y = fitRatio + 0.2, xend = zeroRatio+ 0.02, yend = fitRatio + 0.02, arrow = arrow(length = unit(0.05, "inches"), type = "closed"), color = "darkorchid4", size = 0.6)) + #0.05, 0.5
  with(si.dat[si.dat$max == "MAX",], annotate(geom = "label", x = zeroRatio + 0.26, y = fitRatio + 0.18, color = "darkorchid4", label = substring(sprintf("%.3f", SI), 2), fill = "#FAFAFA", label.padding = unit(0.04, "in"), size = 2))

(fig2d <- siplot1 + theme(text = element_text(size = 6, color = "#42376B")))


## ----fig2e----------------------------------------------------------
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
  geom_point(color = "darkorchid4") + geom_line() +
  theme_bw() + scale_size_continuous(limits = 0:1) + 
  coord_equal(xlim=0:1, ylim = 0:1) + 
  labs(x = "Zero ratio", y = "Fit ratio", size = "Sparsity\nIndex") + 
  theme(
    panel.border = element_rect(size = 1.5, color = "#42376B", fill = NA),
    panel.grid = element_blank()) +
  with(si.dat2[si.dat2$max == "MAX",], annotate(geom = "point", x = zeroRatio, y = fitRatio, color = "#7F6CC7", size = 1)) +
  with(si.dat2[si.dat2$max == "MAX",], annotate(geom = "segment", x = zeroRatio + 0.25, y = fitRatio + 0.2, xend = zeroRatio+ 0.02, yend = fitRatio + 0.02, arrow = arrow(length = unit(0.05, "inches"), type = "closed"), color = "darkorchid4", size = 0.6)) + #0.05, 0.5
  with(si.dat2[si.dat2$max == "MAX",], annotate(geom = "label", x = zeroRatio + 0.26, y = fitRatio + 0.18, color = "darkorchid4", label = substring(sprintf("%.3f", SI), 2), fill = "#FAFAFA", label.padding = unit(0.04, "in"), size = 2))

(fig2e <- siplot2 + theme(text = element_text(size = 6, color = "#42376B")))


## ----fig2f----------------------------------------------------------

Usp.o <- Usp %>% mutate(cl = cl)
Usp.o <- Usp.o[res.hc$order, ] %>% mutate(IDnum = 1:n())

lol1.cl <- ggplot(Usp.o, aes(x = IDnum, y = U1, color = cl)) +
  geom_point(size = 1) +
  geom_segment(aes(xend = IDnum, y = 0,  yend = U1), size = .5) + 
  coord_flip() +
  theme_minimal() + 
  labs(x = "", y = "", color = "Cluster", title = "Subspace 1") + 
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "NA",
    text = element_text(size = 6, color = "#42376B"),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(0.1, 0, 0, 0), unit = "cm"),
    strip.background = element_rect(color = "transparent", fill = "transparent", size = 0.1)
  )

lol2.cl <- ggplot(Usp.o, aes(x = IDnum, y = U2, color = cl)) + 
  geom_point(size = 1) + 
  geom_segment(aes(xend = IDnum, y = 0,  yend = U2), size = .5) + 
  coord_flip() +
  theme_minimal() + 
  labs(x = "", y = "", color = "Cluster", title = "Subspace 2") + 
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "NA",
    text = element_text(size = 6, color = "#42376B"),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(0.1, 0, 0, 0), unit = "cm"),
    strip.background = element_rect(color = "transparent", fill = "transparent", size = 0.1)
  )

lol3.cl <- ggplot(Usp.o, aes(x = IDnum, y = U3, color = cl)) + 
  geom_point(size = 1) + 
  geom_segment(aes(xend = IDnum, y = 0,  yend = U3), size = .5) + 
  coord_flip() +
  theme_minimal() + 
  labs(x = "", y = "", color = "Cluster", title = "Subspace 3") + 
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "NA",
    text = element_text(size = 6, color = "#42376B"),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(0.1, 0, 0, 0), unit = "cm"),
    strip.background = element_rect(color = "transparent", fill = "transparent", size = 0.1)
  )

(fig2f <- ggpubr::ggarrange(lol1.cl, lol2.cl, lol3.cl, ncol = 3, common.legend = TRUE, legend = "bottom") + ggtheme.noborder)



## ----compute all three eigen subspaces------------------------------
f1 <- t(t(eigen.Splus1$vectors) * sqrt(abs(eigen.Splus1$values)))
f2 <- t(t(eigen.Splus2$vectors) * sqrt(abs(eigen.Splus2$values)))
f3 <- t(t(eigen.Splus3$vectors) * sqrt(abs(eigen.Splus3$values)))

rownames(f1) <- rownames(f2) <- rownames(f3) <- dimnames(DisCube)[[1]]
colnames(f1) <- paste("Factor", 1:ncol(f1))
colnames(f2) <- paste("Factor", 1:ncol(f2))
colnames(f3) <- paste("Factor", 1:ncol(f3))


## ----fig2g----------------------------------------------------------
cmp.fi1 <- f1

## labels
cmp.labels1 <- createxyLabels.gen(
  1, 2,
  lambda = eigen.Splus1$values,
  tau = 100 * eigen.Splus1$values / sum(eigen.Splus1$values),
  axisName = "Component ")

## constraints
colnames(cmp.fi1) <- paste0("Dimension ", c(1:ncol(cmp.fi1)))
cmp.constraints1 <- minmaxHelper(cmp.fi1)

## Fi
CMPmap1 <- createFactorMap(
  cmp.fi1,
  constraints = cmp.constraints1,
  text.cex = 1.5,
  cex = 1.5,
  col.background = NULL,
  col.axes = "#42376B",
  col.points = wine.col,
  col.labels = wine.col,
  width.axes = 0.5,
  label.axisName = "Component ",
  alpha.axes = 0.5,
  title = "Subspace 1: global factor \nscores of wines"
)

(fig2g <- CMPmap1$zeMap_background + CMPmap1$zeMap_dots + CMPmap1$zeMap_text + cmp.labels1 + ggtheme)



## ----fig2h----------------------------------------------------------
cmp.fi2 <- f2

## labels
cmp.labels2 <- createxyLabels.gen(
  1, 2,
  lambda = eigen.Splus2$values,
  tau = 100 * eigen.Splus2$values / sum(eigen.Splus2$values),
  axisName = "Component ")

## constraints
colnames(cmp.fi2) <- paste0("Dimension ", c(1:ncol(cmp.fi2)))
cmp.constraints2 <- minmaxHelper(cmp.fi2)

## Fi
CMPmap2 <- createFactorMap(
  cmp.fi2,
  constraints = cmp.constraints2,
  text.cex = 1.5,
  cex = 1.5,
  col.background = NULL,
  col.axes = "#42376B",
  col.points = wine.col,
  col.labels = wine.col,
  width.axes = 0.5,
  label.axisName = "Component ",
  alpha.axes = 0.5,
  title = "Subspace 2: global factor \nscores of wines"
)


(fig2h <- CMPmap2$zeMap_background + CMPmap2$zeMap_dots + CMPmap2$zeMap_text + cmp.labels2 + ggtheme)


## ----fig2i----------------------------------------------------------
cmp.fi3 <- f3

## labels
cmp.labels3 <- createxyLabels.gen(
  1, 2,
  lambda = eigen.Splus3$values,
  tau = 100 * eigen.Splus3$values / sum(eigen.Splus3$values),
  axisName = "Component ")

## constraints
colnames(cmp.fi3) <- paste0("Dimension ", c(1:ncol(cmp.fi3)))
cmp.constraints3 <- minmaxHelper(cmp.fi3)

## Fi
CMPmap3 <- createFactorMap(
  cmp.fi3,
  constraints = cmp.constraints3,
  text.cex = 1.5,
  cex = 1.5,
  col.background = NULL,
  col.axes = "#42376B",
  col.points = wine.col,
  col.labels = wine.col,
  width.axes = 0.5,
  label.axisName = "Component ",
  alpha.axes = 0.5,
  title = "Subspace 3: global factor \nscores of wines"
)


(fig2i <- CMPmap3$zeMap_background + CMPmap3$zeMap_dots + CMPmap3$zeMap_text + cmp.labels3 + ggtheme)


## ----figure 2, include = FALSE--------------------------------------
png(filename = "Figure2-SDiSTATIS.png", width = 140, height = 200, units = "mm", bg = "white", res = 600)
gridExtra::grid.arrange(
  grobs = list(fig2a, fig2c, 
               fig2d, fig2e, fig2f,
               fig2g, fig2h, fig2i),
  widths = c(1, 1, 1, 1, 1, 1),
  heights = c(1, 1, 1, 1),
  layout_matrix = rbind(c(1,1,1,2,2,2),
                        c(3,3,5,5,5,5),
                        c(4,4,5,5,5,5),
                        c(6,6,7,7,8,8)))
dev.off()

