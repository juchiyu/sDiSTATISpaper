plot.sdica.results <- function(sdica.res, dica.res, dim1 = 1, dim2 = 2, y, grp.col) {

  dat.sdicasupp <- data.frame(
    sdica.res$fii,
    cla = y)

  dat.sdica <- data.frame(
    sdica.res$fi,
    cla = unique(y))

  tau <- 100 * sdica.res$eig / sum(dica.res$TExPosition.Data$eigs)

  labz <- sapply(1:2, function(i) sprintf("Dim %i (%2.1f%%)", i, tau[i]))

  gind <- ggplot(dat.sdicasupp, aes(X1, X2, color = cla)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point() + stat_ellipse() +
    geom_label(data = dat.sdica, aes(label = cla), show.legend = FALSE) +
    theme_bw() +
    labs(x = labz[1], y = labz[2], color = "Group")

  print(gind)


  sdica.fj <- PTCA4CATA::createFactorMap(sdica.res$fj,
                                         col.points = grp.col,
                                         col.labels = grp.col,
                                         text.cex = 2, alpha.labels = 0.5,
                                         col.background = NULL,
                                         col.axes = "orchid4",
                                         alpha.axes = 0.5)
  gvar <- sdica.fj$zeMap_background + sdica.fj$zeMap_dots + sdica.fj$zeMap_text

  print(gvar)

}
