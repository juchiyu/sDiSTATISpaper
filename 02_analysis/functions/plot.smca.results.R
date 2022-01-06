plot.smca.results <- function(res.smca, res.mca, dim1 = 1, dim2 = 2, color = c("#458b00", "#510e53")) {
  
  FacSco <- data.frame(
    PC  = res.smca$fj,
    item = rownames(res.smca$fj),
    question = gsub("\\.[a-z]$", "", rownames(res.smca$fj)))
  
  FacSco$questionType <- ifelse(as.numeric(gsub("^Q", "", FacSco$question)) < 22, "Indep", "Interdep")
  
  tau <- 100 * res.smca$eig / sum(res.mca$eig[,1])  
  
  
  ggplot(FacSco, aes_string(sprintf("PC.%i", dim1) , sprintf("PC.%i", dim2), color = "questionType", label = "item")) + 
    geom_hline(yintercept = 0, color = "grey") +
    geom_vline(xintercept = 0, color = "grey") +
    geom_point(size = 2, shape = 17) + 
    geom_text_repel() +
    scale_color_manual(values = color) +
    theme_minimal() + 
    labs(x = sprintf("Dim%i (%2.1f%%)", dim1, tau[dim1]),
         y = sprintf("Dim%i (%2.1f%%)", dim2, tau[dim2]))
  
}