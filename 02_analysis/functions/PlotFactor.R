PlotFactor <- function(lambda = lambda,
                       tau = tau,
                       digit4lambda = 3, 
                       digit4tau = 0, 
                       f = f,
                       design = NULL,
                       col.list = col.list,
                       xaxis = 1,
                       yaxis = 2,
                       label.axisName = "Component ",
                       title = NULL,
                       alpha.indiv.point = 0.1,
                       pch.indiv.point = 19,
                       alpha.mean.point = 0.8,
                       pch.mean.point = 17,
                       cex.mean.point = 3,
                       cex.mean.text = 3,
                       mean.constraints = NULL,
                       indiv.constraints = NULL
){
  require(PTCA4CATA)
  map.label <- createxyLabels.gen.digit(lambda = lambda,
                                  tau = tau,
                                  digit4lambda = digit4lambda, 
                                  digit4tau = digit4tau, 
                                  x_axis = xaxis,
                                  y_axis = yaxis, 
                                  axisName = label.axisName)
  # Row factor scores for individuals - compute
  colnames(f) <- paste0("Dimension ", c(1:ncol(f)))
  if(!is.null(design)){  
    f.mean <- getMeans(f, design)
    f.boot <- Boot4Mean(f, design, niter = 1000)
    colnames(f.boot$BootCube) <- paste0("Dimension ", c(1:ncol(f.boot$BootCube)))
    
    col4ind <- col.list$oc
    col4grp <- col.list$gc
  }else{
    col4ind <- col.list
    col4grp <- col.list
  }
  
  if(is.null(mean.constraints)){
    mean.constraints = minmaxHelper4Brick(f.boot$BootCube[,c(xaxis,yaxis),])
  }
  if(is.null(indiv.constraints)){
    indiv.constraints = minmaxHelper(f, axis1 = xaxis, axis2 = yaxis)
  }
  
  # plotting
  f.map <- createFactorMap(f,
                           col.points = col4ind, 
                           col.labels = col4ind,
                           col.background = NULL,
                           col.axes = "#42376B",width.axes = 0.5,
                           title = title,
                           alpha.axes = 0.5,
                           alpha.points = alpha.indiv.point,
                           pch = pch.indiv.point,
                           axis1 = xaxis,
                           axis2 = yaxis,
                           constraints = indiv.constraints)
  
  # fi.map$zeMap_background + fi.map$zeMap_dots + map.label
  
  if (!is.null(design)){
  f.mean.plot <- createFactorMap(f.mean,
                                 col.points = col4grp[rownames(f.mean),],
                                 col.labels = col4grp[rownames(f.mean),],
                                 col.background = NULL,
                                 col.axes = "#42376B",
                                 alpha.axes = 0.5,
                                 alpha.points = alpha.mean.point,
                                 cex = cex.mean.point,
                                 text.cex = cex.mean.text,
                                 pch = pch.mean.point,
                                 axis1 = xaxis,
                                 axis2 = yaxis,
                                 constraints = mean.constraints)
  
  f.CI <- MakeCIEllipses(f.boot$BootCube[,c(xaxis,yaxis),],
                         names.of.factors = paste0("Dimension ", c(xaxis,yaxis)),
                         col = col4grp[rownames(f.boot$BootCube),],
                         alpha.ellipse = 0.1, 
                         line.size = 0.5, alpha.line = 0.2)
  f.TI <- MakeToleranceIntervals(f[,c(xaxis,yaxis)],
                                 design = design,
                                 col = col4grp[rownames(f.mean),],
                                 names.of.factors = paste0("Dimension ", c(xaxis, yaxis)),
                                 alpha.ellipse = 0.05)
  }
  if(!is.null(design)){
    return(list(f.map = f.map, mean.map = f.mean.plot, CI = f.CI, TI = f.TI, label = map.label))
  }else{
    return(list(f.map = f.map, label = map.label))
    }
  
  
  # example for output
  # f.map$zeMap_background + f.map$zeMap_dots +
  #   TI + CI +
  #   mean.plot$zeMap_dots + mean.plot$zeMap_text + label
}