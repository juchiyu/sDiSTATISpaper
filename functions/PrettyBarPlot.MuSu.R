PrettyBarPlot.MuSu <- function(x, y, grp, label.bar, color, bar.label.size = 1, grp.label.size = 1, grp.label.position = 1.2, title = NULL, ylab = NULL){
  theme4bar <- theme_bw() +
    theme(
      axis.line.y = element_line(colour = "black"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.position = "None")

  data.frame(x = x, y = y, grp = grp, label.bar = label.bar, IDnum = as.numeric(as.factor(x))) %>%
    ggplot(aes(x = x, y = y, fill = grp)) +
    geom_col(show.legend = FALSE) +
    geom_text(
      aes(
        y = y + sign(y) * 0.01,
        label = ifelse(near(y,0), "", label.bar),
        color = grp,
        hjust = +(y < 0),
      ),
      angle = 90,
      size = bar.label.size
    ) +
    geom_text(
      mapping = aes(
        x = meanID,
        y = grp.label.position,
        label = grp,
        color = ifelse(near(y,0), "ns", grp)
      ),
      size = grp.label.size,
      fontface = "bold",
      data = . %>% group_by(grp) %>% summarize(meanID = mean(IDnum), y = mean(abs(y)))
    ) +
    geom_hline(yintercept = 0) +
    ylab(ylab) +
    ggtitle(title) +
    scale_x_continuous(expand = c(0,0)) +
    scale_color_manual(values = c(color, "ns" = "grey80")) +
    scale_fill_manual(values = c(color, "ns" = "grey80")) +
    theme4bar
}
