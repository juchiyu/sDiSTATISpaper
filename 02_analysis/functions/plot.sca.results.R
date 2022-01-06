plot.sca.results <- function(res.sca, res.ca, dim1 = 1, dim2 = 2, color.row.col = c("black", "red")) {
  
  dat.caplot.rows <- data.frame(
    PC  = res.sca$fi,
    chapter = rownames(res.sca$fi))
  
  dat.caplot.cols <- data.frame(
    PC  = res.sca$fj,
    chapter = rownames(res.sca$fj))
  
  tau <- 100 * res.sca$eig / sum(res.ca$eig[,1])  
  
  
  ggplot() + 
    geom_hline(yintercept = 0, color = "grey") +
    geom_vline(xintercept = 0, color = "grey") +
    geom_point(data = dat.caplot.rows, mapping = aes_string(sprintf("PC.%i", dim1) , sprintf("PC.%i", dim2)), color = color.row.col[1], size = 2) + 
    geom_point(data = dat.caplot.cols, mapping = aes_string(sprintf("PC.%i", dim1) , sprintf("PC.%i", dim2)), color = color.row.col[2], size = 2, shape = 17) + 
    geom_text_repel(data = dat.caplot.rows, mapping = aes_string(sprintf("PC.%i", dim1) , sprintf("PC.%i", dim2), label = "chapter"), color = color.row.col[1]) + 
    geom_text_repel(data = dat.caplot.cols, mapping = aes_string(sprintf("PC.%i", dim1) , sprintf("PC.%i", dim2), label = "chapter"), color = color.row.col[2]) + 
    theme_minimal() + 
    labs(x = sprintf("Dim%i (%2.1f%%)", dim1, tau[dim1]),
         y = sprintf("Dim%i (%2.1f%%)", dim2, tau[dim2]))
  
}