PlotMyScree <- function(eigres, color.sig =  "#42376B", color.ns = "grey60",
                        cex = 1.1, text.cex = 10, lwd = 1, title = NULL){
  eigres %>% as.data.frame %>% ggplot(aes(x = 1:length(eig), y = tau)) +
    geom_line(color = "grey40", size = lwd) +
    geom_point(color = color.sig, size = cex) +
    # geom_hline(yintercept = 1, linetype = "dashed", color = "darkgreen", size = lwd) +
    scale_y_continuous(name = bquote(atop(bold(.(title)),paste('\n\n    Percentage of \nvariance explained (%)'))),
                       sec.axis = sec_axis(~.*(eigres$eig[1]/eigres$tau[1]), name = "Pseudo-eigenvalues")) +
    xlab("Components") +
    scale_x_continuous(breaks=c(1:nrow(eigres))) +
    theme(text = element_text(size = text.cex),
          legend.position = "none",
          axis.text.y.left = element_text(angle = 90),
          axis.text.y.right = element_text(angle = 270),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(color = "black", fill = "transparent"))
}
