## ----setup, include=FALSE-------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
# detach("package:SPAFAC", unload = TRUE)
# devtools::install_github("juchiyu/SPAFAC")

library(FactoMineR)
library(factoextra)
library(PTCA4CATA)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(plotly)
library(RColorBrewer)
library(kableExtra)
library(ggpubr)
library(ggplotify)
library(SPAFAC)
library(sGSVD)
library(inlmisc)
# devtools:install_github("jokergoo/ComplexHeatmap")
# library(ComplexHeatmap)

load("../01_data/MFA_turkeysensory.rda")

source("../02_analysis_sparseDiSTATIS/functions/PlotFactor.R")
source("../02_analysis_sparseDiSTATIS/functions/createxyLabels.gen.digit.R")
source("../02_analysis_sparseDiSTATIS/functions/createLabel.gen.digit.R")
source("../02_analysis_sparseDiSTATIS/functions/plot.smca.results.R")
source("../02_analysis_sparseDiSTATIS/functions/PlotMyScree.R")
source("../02_analysis_sparseDiSTATIS/functions/PrettyBarPlot.MuSu.R")

# other functions
extract_si <- function(sdisca, si = "SI") {
  if (class(sdisca) == "list") {
    return(sdisca$sparsity$SI[[si]])
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

rownames(grand.tab) <- gsub("PECHUGA DE PAVO ", "", rownames(grand.tab))
rownames(grand.tab) <- c("HORNEADO", "SABORI", "DE FUD", "CAMPESTRE", "VIRGINIA", "NATURAL", "CLÁSICA", "ALPINO")
raters <- sapply(strsplit(colnames(grand.tab), "\\."), "[", 2)
flavors <- sapply(strsplit(colnames(grand.tab), "\\."), "[", 1)

## get color ------------------
col.raters <- list()
col.raters$gc <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#E8C245','#a65628','#f781bf') %>% as.matrix
col.raters$oc <- rep(col.raters$gc, each = 12) %>% as.matrix
rownames(col.raters$gc) <- unique(raters)
rownames(col.raters$oc) <- colnames(grand.tab)

## design by flavor
flavors.coldx <- data.frame(t(flavors))
colnames(flavors.coldx) <- colnames(grand.tab)

## when the columns are ordered by flavors (instead of raters)
col.flavor.idx <- c("#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac", "#c3bc3f", "#8cc2ca")
names(col.flavor.idx) <- unique(flavors)
col.flavors <- list()
col.flavors$gc <- col.flavor.idx %>% as.matrix()
col.flavors$oc <- recode(flavors, !!!col.flavor.idx) %>% as.matrix()
rownames(col.flavors$oc) <- colnames(grand.tab)

## get theme -----------------
ggtheme <- theme(axis.title = element_text(size = 6, color = "#42376B"), axis.text = element_text(size = 6, color = "#42376B"), title = element_text(size = 6, color = "#42376B"), panel.border = element_rect(size = 1.5, color = "#42376B", fill = NA))
ggtheme.noborder <- theme(axis.title = element_text(size = 6, color = "#42376B"), axis.text = element_text(size = 6, color = "#42376B"), title = element_text(size = 6, color = "#42376B"))


## -------------------------------------------------------------------
res.mfa <- MFA(grand.tab, group = tapply(raters, raters, length), name.group = unique(raters), graph = FALSE)
## reorganize partial fi
mfa.partial.fi <- lapply(split(res.mfa$ind$coord.partiel, 0:63 %/% 8), matrix, nrow = 8) %>% simplify2array %>% aperm(c(3,2,1))
dimnames(mfa.partial.fi) <- list(dimnames(res.mfa$ind$coord)[[1]], paste0("Dimension ", c(1:ncol(res.mfa$ind$coord.partiel))), rownames(res.mfa$group$coord))

colnames(res.mfa$global.pca$var$coord) <- colnames(res.mfa$ind$coord) <- paste0("Dimension ", c(1:ncol(res.mfa$ind$coord)))



## -------------------------------------------------------------------
eigres <- as.data.frame(res.mfa$eig[,c(1:2)])
colnames(eigres) <- c("eig", "tau")
PlotMyScree(eigres, cex = 2, text.cex = 8, lwd = 0.5)+ theme(axis.title = element_text(size = 6, color = "#42376B"), axis.text = element_text(size = 6, color = "#42376B"), panel.border = element_rect(size = 1.5, fill = NA, color = "#42376B"))


## -------------------------------------------------------------------
## Lables
mfa.label <- createxyLabels.gen(1,2,
                                lambda = eigres$eig,
                                tau = eigres$tau,
                                axisName = "Component "
)

## mfa constraints
mfa.constraints <- minmaxHelper4Partial(res.mfa$ind$coord,
                                        mfa.partial.fi)

## row factor scores
mfa.fi <- createFactorMap(res.mfa$ind$coord,
                          constraints = mfa.constraints,
                          text.cex = 2,
                          cex = 1.5,
                          col.points = "#42376B",
                          col.labels = "#42376B",
                          col.background = NULL,
                          col.axes = "#42376B",
                          width.axes = 0.5,
                          label.axisName = "Component ",
                          alpha.axes = 0.5,
                          title = "MFA: global and partial\nfactor scores of products\n(rows)")
# mfa.fi$zeMap + ggtheme + mfa.label

## partial factor scores
mfa.pfi <- createPartialFactorScoresMap(res.mfa$ind$coord,
                                        mfa.partial.fi,
                                        colors4Blocks = col.raters$gc,
                                        colors4Items = "#42376B",
                                        alpha.lines = 0.7, alpha.points = 0.7, size.lines = 0.5
)

mfa.fig1 <- mfa.fi$zeMap_background + mfa.pfi$linesColByBlocks + mfa.pfi$pointsColByBlocks + mfa.fi$zeMap_dots + mfa.fi$zeMap_text + mfa.label + ggtheme

# column factor scores
mfa.fj <- createFactorMap(res.mfa$global.pca$var$coord,
                          col.labels = col.flavors$oc,
                          col.points = col.flavors$oc,
                          text.cex = 2,
                          col.background = NULL,
                          col.axes = "#42376B",
                          width.axes = 0.5,
                          label.axisName = "Component ",
                          alpha.axes = 0.5,
                          title = "MFA: column factor scores of \nraters and flavors (columns)")
mfa.fig2 <- mfa.fj$zeMap + mfa.label + ggtheme

# mfa.fj.bar1 <- PrettyBarPlot2(res.mfa$global.pca$var$coord[,1],
#                color4bar = col.raters$oc,threshold = 0, main = "MFA: column factor scores - Component 1") + ggtheme.noborder
# mfa.fj.bar2 <- PrettyBarPlot2(res.mfa$global.pca$var$coord[,2],
#                color4bar = col.raters$oc,threshold = 0, main = "MFA: column factor scores - Component 2") + ggtheme.noborder



## -------------------------------------------------------------------
## bars from ggplot ===========

### rearrange the data
dat4bar.mfa <- data.frame(res.mfa$global.pca$var$coord[,c(1,2)], flavors = flavors, raters = raters, 
                          ID = rownames(res.mfa$global.pca$var$coord), IDnum = c(1:length(flavors)))

### new bars
(mfa.fj.bar1 <- PrettyBarPlot.MuSu(dat4bar.mfa$IDnum, dat4bar.mfa$Dimension.1, grp = raters, 
                                  color = drop(col.raters$gc), 
                                  label.bar = flavors,
                                  grp.label.size = 1.5,
                                  grp.label.position = 1.4,
                                  bar.label.size = 1.25,
                                  title = "MFA: column factor scores",
                                  ylab = "Component 1") + ylim(c(-1.2,1.4)) + ggtheme.noborder)

mfa.fj.bar2 <- PrettyBarPlot.MuSu(dat4bar.mfa$IDnum, dat4bar.mfa$Dimension.2, grp = raters, 
                                  color = drop(col.raters$gc), 
                                  label.bar = flavors,
                                  grp.label.size = 1.5,
                                  grp.label.position = 1.5,
                                  bar.label.size = 1.25,
                                  title = "MFA: column factor scores",
                                  ylab = "Component 2") + ylim(c(-1.5,1.5)) + ggtheme.noborder


## ----allparz, cache = TRUE, include = FALSE, message = FALSE, warning=FALSE----
### Run iterations
K <- 5L

I <- 8
J <- ncol(grand.tab)
J.raters <- length(unique(raters))
J.flavor <- length(unique(flavors))

parz <- expand.grid(
  rdsleft = sqrt(I), #seq(sqrt(I), 1, length = 15),
  rdsright = seq(sqrt(J.raters), 1, length = 15))

parz <- parz[order(rowSums(parz^2/c(I,J.raters)), decreasing = T), ]

iter <- 1
res.smfa.list.byraters <- NULL
U0 <- NULL
V0 <- NULL

for (rds.iter in 1:nrow(parz)) {
  cat(sprintf("Left radius: %0.2f - Right radius: %0.2f\n", parz[rds.iter, 1], parz[rds.iter, 2]))
  res.smfa.list.byraters[[iter]] <- tryCatch(
    sparseMFA(X = as.matrix(grand.tab), column.design = raters,
              components = K,
              sparseOption = "subtable",
              center = TRUE, scale = TRUE, mfa.scale = TRUE,
              orthogonality = "loadings",
              rdsLeft = rep(parz[rds.iter, 1], K), rdsRight = rep(parz[rds.iter, 2], K), 
              grpLeft = NULL), error = function(e) NA
  )
  if (!is.na(res.smfa.list.byraters[[iter]])){
    U0 <- res.smfa.list.byraters[[iter]]$gsvd$p
    V0 <- res.smfa.list.byraters[[iter]]$gsvd$q
  }
  iter <- iter + 1
}

### record results
dat.si.byraters <- data.frame(
  parz = parz, 
  SI = unname(t(data.frame(SI = sapply(res.smfa.list.byraters, extract_si, "SI")))))

dat.si.m.byraters <- dat.si.byraters %>% tidyr::pivot_longer(starts_with("SI"), names_to = "k", names_prefix = "SI\\.")

dat.fit.zeros.byraters <- data.frame(
  rdsleft = dat.si.m.byraters$parz.rdsleft, 
  rdsright = dat.si.m.byraters$parz.rdsright, 
  k = dat.si.m.byraters$k,
  fit = pivot_the_tab(sapply(res.smfa.list.byraters, extract_si, "r1")),
  zeros = pivot_the_tab(sapply(res.smfa.list.byraters, extract_si, "r4")),
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
    legend.key.size = unit(0.1, "cm"), legend.text = element_text(size = 6), 
    legend.position = "right",
    legend.spacing.x = unit(0.1, "cm"),
    legend.margin = margin(0,0,0,0), legend.box.margin = margin(-5,0,0,0),
    legend.title.align = 0, legend.title = element_text(size = 6),
    panel.border = element_rect(size = 1.5, color = "#42376B", fill = NA),
    plot.margin = unit(c(0.1,0,0,0), unit = "cm"),
    panel.grid = element_blank()) + 
  with(dat.fit.zeros.byraters[dat.fit.zeros.byraters$max == "MAX",], annotate(geom = "point", x = zeros, y = fit, alpha = alpha, color = dimcol[k], size = 1)) +
  with(dat.fit.zeros.byraters[dat.fit.zeros.byraters$max == "MAX",], annotate(geom = "segment", x = zeros + 0.25, y = fit + 0.2, xend = zeros+ 0.02, yend = fit + 0.02, arrow = arrow(length = unit(0.05, "inches"), type = "closed"), color = "darkorchid4", size = 0.6)) + #0.05, 0.5
  # with(dat.fit.zeros.byraters[dat.fit.zeros.byraters$max == "MAX",], annotate(geom = "segment", x = zeros + 0.025, y = fit + 0.025, xend = zeros+ 0.02, yend = fit + 0.02, arrow = arrow(length = unit(0.04, "inches"), type = "closed"), color = dimcol[k], size = 0.01)) + #0.04, 0.01
  with(dat.fit.zeros.byraters[dat.fit.zeros.byraters$max == "MAX",], annotate(geom = "label", x = zeros + 0.26, y = fit + 0.18, color = "darkorchid4", label = substring(sprintf("%.3f", SI), 2), fill = dimcol[k], size = 2, label.padding = unit(0.04, "in")))

siplot


## -------------------------------------------------------------------
optres <- dat.fit.zeros.byraters[dat.fit.zeros.byraters$max == "MAX",]
kopt <- optres$k
rdsleftopt <- optres$rdsleft
rdsrightopt <- optres$rdsright
SIopt <- optres$SI

kopt;rdsleftopt/sqrt(I);rdsrightopt/sqrt(J);SIopt


## -------------------------------------------------------------------
res.smfa <- sparseMFA(X = as.matrix(grand.tab), column.design = raters,
                      components = 7L, #kopt,
                      sparseOption = "subtable",
                      center = TRUE, scale = TRUE, mfa.scale = TRUE,
                      orthogonality = "loadings",
                      rdsLeft = rep(rdsleftopt, 7L), rdsRight = rep(rdsrightopt, 7L), 
                      grpLeft = NULL)
# res.smfa <- sparseMFA(X = as.matrix(grand.tab), column.design = raters,
#                       components = 7L, #kopt,
#                       sparseOption = "subtable",
#                       center = TRUE, scale = TRUE, mfa.scale = TRUE,
#                       orthogonality = "loadings",
#                       rdsLeft = rep(sqrt(I), 7L), rdsRight = rep(sqrt(J.raters), 7L), 
#                       grpLeft = NULL)
colnames(res.smfa$fi) <- colnames(res.smfa$fj) <- colnames(res.smfa$partial.fi) <- paste0("Dimension ", c(1:ncol(res.smfa$fi)))
dimnames(res.smfa$partial.fi)[[3]] <- unique(raters)



## -------------------------------------------------------------------
smfa.eigres <- data.frame(eig = res.smfa$eig,
                          tau = 100*res.smfa$eig/sum(res.mfa$eig[,1]))
smfa.scree <- PlotMyScree(smfa.eigres, cex = 2, text.cex = 8, lwd = 0.5, color.sig = c(rep("#42376B", 2), rep("grey80", 5)))+ theme(axis.title = element_text(size = 6, color = "#42376B"), axis.text = element_text(size = 6, color = "#42376B"), panel.border = element_rect(size = 1.5, fill = NA, color = "#42376B"))
smfa.scree


## -------------------------------------------------------------------
## Lables
smfa.label <- createxyLabels.gen.digit(1,2,
                                 lambda = res.smfa$eig,
                                 tau = 100*res.smfa$eig/sum(eigres$eig),
                                 digit4tau = 2,
                                 axisName = "Component "
)

## mfa constraints
smfa.constraints <- minmaxHelper4Partial(res.smfa$fi,
                                         res.smfa$partial.fi)

## row factor scores
smfa.fi <- createFactorMap(res.smfa$fi,
                           constraints = smfa.constraints,
                           text.cex = 2,
                           cex = 1.5,
                           col.points = "#42376B",
                           col.labels = "#42376B",
                           col.background = NULL,
                           col.axes = "#42376B",
                           width.axes = 0.5,
                           label.axisName = "Component ",
                           alpha.axes = 0.5,
                           title = "sMFA: global and partial\nfactor scores of products\n(rows)")
# mfa.fi$zeMap + ggtheme + mfa.label

## partial factor scores
smfa.pfi <- createPartialFactorScoresMap(res.smfa$fi,
                                         res.smfa$partial.fi,
                                         colors4Blocks = col.raters$gc,
                                         colors4Items = "#42376B",
                                         alpha.lines = 0.7, alpha.points = 0.7, size.lines = 0.5
)

smfa.fig1 <- smfa.fi$zeMap_background + smfa.pfi$linesColByBlocks + smfa.pfi$pointsColByBlocks + smfa.fi$zeMap_dots + smfa.fi$zeMap_text + smfa.label + ggtheme

# column factor scores
smfa.fj <- createFactorMap(res.smfa$fj,
                           col.labels = col.raters$oc,
                           col.points = col.raters$oc,
                           text.cex = 2,
                           col.background = NULL,
                           col.axes = "#42376B",
                           width.axes = 0.5,
                           label.axisName = "Component ",
                           alpha.axes = 0.5,
                           title = "sMFA: column factor scores of \nraters and flavors (columns)")
smfa.fig2 <- smfa.fj$zeMap + smfa.label + ggtheme
smfa.fj.bar1 <- PrettyBarPlot2(res.smfa$fj[,1],
                               color4bar = col.raters$oc, color4ns = "white", threshold = 0.001, line.col = "white", signifOnly = FALSE, main = "sMFA: column factor scores - Component 1") + ggtheme.noborder
smfa.fj.bar2 <- PrettyBarPlot2(res.smfa$fj[,2],
                               color4bar = col.raters$oc, color4ns = "white", threshold = 0.001, line.col = "white", signifOnly = FALSE, main = "sMFA: column factor scores - Component 1") + ggtheme.noborder


## -------------------------------------------------------------------
## bars from ggplot ===========

### rearrange the data
dat4bar.smfa <- data.frame(res.smfa$fj[,c(1,2)], flavors = flavors, raters = raters, 
                          ID = rownames(res.smfa$fj), IDnum = c(1:length(flavors)))
smfa.color.gc <- drop(col.raters$gc)

### new bars
(smfa.fj.bar1 <- PrettyBarPlot.MuSu(dat4bar.smfa$IDnum, dat4bar.smfa$Dimension.1, grp = raters, 
                                  color = drop(col.raters$gc), 
                                  label.bar = flavors,
                                  grp.label.size = 1.5,
                                  grp.label.position = 1.4,
                                  bar.label.size = 1.25,
                                  title = "sMFA: column factor scores",
                                  ylab = "Component 1") + ylim(c(-1,1.5))+ ggtheme.noborder)

(smfa.fj.bar2 <- PrettyBarPlot.MuSu(dat4bar.smfa$IDnum, dat4bar.smfa$Dimension.2, grp = raters, 
                                  color = drop(col.raters$gc), 
                                  label.bar = flavors,
                                  grp.label.size = 1.52,
                                  grp.label.position = 1.4,
                                  bar.label.size = 1.25,
                                  title = "sMFA: column factor scores",
                                  ylab = "Component 2") + ylim(c(-0.8,1.52)) + ggtheme.noborder)


## -------------------------------------------------------------------
# png(filename="Figure-SMFA.png", width = 140, height = 200, units = "mm", bg = "white",res = 800)
# gridExtra::grid.arrange(grobs = list(smfa.scree, siplot,
#                                      mfa.fig1, smfa.fig1,
#                                      mfa.fj.bar1, smfa.fj.bar1,
#                                      mfa.fj.bar2, smfa.fj.bar2),
#                         widths = c(0.5, 0.5),
#                         heights = c(0.2,0.40,0.2, 0.2),
#                         layout_matrix = rbind(c(1,2),
#                                               c(3,4),
#                                               c(5,6),
#                                               c(7,8))
# )
# dev.off()

png(filename="Figure-SMFA(change).png", width = 140, height = 200, units = "mm", bg = "white",res = 800)
gridExtra::grid.arrange(grobs = list(smfa.scree, siplot,
                                     mfa.fig1, smfa.fig1,
                                     mfa.fj.bar1, smfa.fj.bar1,
                                     mfa.fj.bar2, smfa.fj.bar2),
                        widths = c(0.33, 0.17, 0.5),
                        heights = c(0.2,0.2,0.2,0.2, 0.2),
                        layout_matrix = rbind(c(1,1,2),
                                              c(3,5,5),
                                              c(3,7,7),
                                              c(4,6,6),
                                              c(4,8,8))
)
dev.off()

