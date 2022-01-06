###############
## Run sDiCA ##
###############

## ---- check.fi ----
sdica.fi <- res.dica.test$grp$coord; sdica.fi %>% head

## ---- check.fj ----
sdica.fj <- res.dica.test$var$coord[1:6, 1:4]; sdica.fj

## ---- check.eig ----
sdica.eig <- res.dica.test$eig$eigenvalue; sdica.eig

## ---- plot.scree ----
PlotScree(res.dica.test$eig$eigenvalue)

## ---- ind.map ----
dica.test.ind <- PlotFactor(lambda = res.dica.test$eig$eigenvalue,
                            tau = res.dica.test$eig$percentageOfVariance,
                            f = res.dica.test$grp$coord,
                            design = NULL,
                            col.list = ind.col$gc[rownames(res.dica.test$grp$coord),],
                            alpha.indiv.point = 0.5,
                            pch.indiv.point = 19,
                            xaxis = x2plot,
                            yaxis = y2plot
)

dica.test.ind$f.map$zeMap + dica.test.ind$label

for (dim in c(x2plot, y2plot)){
  PrettyBarPlot2(res.dica.test$grp$coord[,dim], threshold = 0, #1/nrow(res.dica.test$ind$coord), 
                 font.size = 5,
                 color4bar = gplots::col2hex(ind.col$gc[names(res.dica.test$grp$coord[,dim]),])) %>% print
  
}

## ---- var.map ----
dica.test.var.cj <- getMeans(res.dica.test$var$eta2, factor = Demo2use4col$Overlapped.Gene, FUN = sum) %>% as.matrix
## multiply by the inverse of the eigenvalues
dica.test.var.cj.relative <- dica.test.var.cj %*% diag(res.dica.test$eig$eigenvalue^-1)
dica.test.var <- createFactorMap(dica.test.var.cj.relative,
                                 text.cex = 3,
                                 col.points = col.gene$gc[rownames(dica.test.var.cj),],
                                 col.labels = col.gene$gc[rownames(dica.test.var.cj),],
                                 pch = 17,
                                 alpha.points = 0.8,
                                 axis1 = x2plot,
                                 axis2 = y2plot)

dica.test.var$zeMap_background + dica.test.var$zeMap_dots + dica.test.var$zeMap_text +
  dica.test.ind$label

dica.test.barx <- PrettyBarPlot2(dica.test.var.cj[,x2plot], threshold = 0,# 1/nrow(res.dica.test$var$eta2), 
                                 color.letter = col.gene$gc[rownames(dica.test.var.cj),],
                                 color4bar = col.gene$gc[rownames(dica.test.var.cj),],
                                 color4ns = "grey90")
dica.test.bary <- PrettyBarPlot2(dica.test.var.cj[,y2plot], threshold = 0,# 1/nrow(res.dica.test$var$eta2), 
                                 color.letter = col.gene$gc[rownames(dica.test.var.cj),],
                                 color4bar = col.gene$gc[rownames(dica.test.var.cj),],
                                 color4ns = "grey90")

dica.test.barx
dica.test.bary

## ---- lev.map ----
dica.test.lev <- createFactorMap(res.dica.test$var$coord,
                                 col.points = col.coord$oc[rownames(res.dica.test$var$coord),],
                                 axis1 = x2plot,
                                 axis2 = y2plot)

dica.test.lev$zeMap_background + dica.test.lev$zeMap_dots +
  dica.test.ind$label

PrettyBarPlot2(res.dica.test$var$coord[,x2plot], threshold = 1/nrow(res.dica.test$var$coord),
               color.letter = col.coord$oc[rownames(res.dica.test$var$coord),],
               color4bar = col.coord$oc[rownames(res.dica.test$var$coord),],
               color4ns = "grey60",
               font.size = 2, signifOnly = TRUE)
PrettyBarPlot2(res.dica.test$var$coord[,y2plot], threshold = 1/nrow(res.dica.test$var$coord),
               color.letter = col.coord$oc[rownames(res.dica.test$var$coord),],
               color4bar = col.coord$oc[rownames(res.dica.test$var$coord),],
               color4ns = "grey90",font.size = 2, signifOnly = TRUE)
