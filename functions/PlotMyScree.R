PlotMyScree <- function(eigres, grouped = FALSE, group.col = NULL, color.sig =  "#42376B", color.ns = "grey60",
                        cex = 1.1, text.cex = 10, lwd = 1, title = NULL){
  data2plot <- eigres %>% as.data.frame
  if (grouped){
    if (is.null(group.col)){
      stop("I am not creating colors for you!")
    }
    data2plot %<>% group_by(group) %>% mutate(x = 1:n())
    p <- ggplot(data2plot, aes(x = x, y = tau, color = group)) +
      geom_line(size = lwd) +
      geom_point(size = cex) +
      scale_x_continuous(breaks=c(1:max(data2plot$x))) +
      scale_color_manual(values = group.col)
  }else{
    p <- ggplot(data2plot, aes(x = 1:length(eig), y = tau)) +
      geom_line(color = "grey40", size = lwd) +
      geom_point(color = color.sig, size = cex) +
      scale_x_continuous(breaks=c(1:nrow(eigres)))
  }
    # geom_hline(yintercept = 1, linetype = "dashed", color = "darkgreen", size = lwd) +
    p + scale_y_continuous(name = bquote(atop(bold(), paste('\n\n    Percentage of \nvariance explained (%)'))),
                       sec.axis = sec_axis(~.*(eigres$eig[1]/eigres$tau[1]), name = "Pseudo-eigenvalues")) +
      ggtitle(title) +
    xlab("Components") +
    theme(text = element_text(size = text.cex),
          legend.position = "none",
          plot.title = element_text(face = "bold", size = text.cex),
          axis.text.y.left = element_text(angle = 90, hjust = 0.5),
          axis.text.y.right = element_text(angle = 270, hjust = 0.5),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(color = "black", fill = "transparent"))
}
