theta <- seq(pi, 3/2*pi, length.out = 150)

dat.fit.zeros <- data.frame(
  zeros = rep(seq(0, 1, length = 10), 2),
  fit = c(seq(1, 0, length = 10)^(1/4), seq(1, 0, length = 10)^(1/2)),
  k = factor(rep(1:2, e = 10))
)

dat.fit.zeros$SI <- with(dat.fit.zeros, 1 - sqrt((1-zeros)^2 + (1-fit)^2))
dat.fit.zeros$max <- ifelse(dat.fit.zeros$SI == max(dat.fit.zeros$SI), "MAX", "NOTMAX")
dat.fit.zeros$MAX <- ifelse(dat.fit.zeros$SI == max(dat.fit.zeros$SI), "darkorchid4", "white")
dat.fit.zeros$alpha <- ifelse(dat.fit.zeros$SI == max(dat.fit.zeros$SI), 1, 0.3)



siplot <- ggplot(dat.fit.zeros, aes(zeros, fit)) + 
  geom_hline(yintercept = 1, color = "darkorchid4", alpha = 0.6, linetype = 3) + 
  geom_vline(xintercept = 1, color = "darkorchid4", alpha = 0.6, linetype = 3) + 
  geom_abline(intercept = 0, slope = 1, color = "darkorchid4", alpha = 0.6, linetype = 3) +
  geom_line(aes(color = k), alpha = 0.4) +
  geom_point(aes(size = SI, fill = k, alpha = I(alpha)), color = dat.fit.zeros$MAX, shape = 21, stroke = 1) +
  theme_bw() + 
  # scale_shape_manual(values = 21) + 
  # scale_color_manual(values = c("NOMAX" = 0, "MAX" = "darkorchid4")) + 
  coord_equal(xlim=0:1, ylim = 0:1) + 
  labs(x = "Zero ratio", y = "Fit ratio", color = "Dim.", size = "Sparsity\nIndex") + 
  ggtitle("Ratio of zeros as a function of the fit ratio") + 
  theme(
    panel.grid = element_blank()) + 
  lapply(seq(0.25, 1.25, by = 0.25), function(r) annotate("path",
                                                          x = 1 + r*cos(theta),
                                                          y = 1 + r*sin(theta),
                                                          color = "darkorchid4")) + 
  lapply(seq(0.125, 1.5, by = 0.25), function(r) annotate("path",
                                                          x = 1 + r*cos(theta),
                                                          y = 1 + r*sin(theta),
                                                          color = "darkorchid4", size = 0.2) )  + 
  with(dat.fit.zeros[which.max(dat.fit.zeros$SI),], annotate(geom = "segment", x = zeros + 0.1, y = fit + 0.1, xend = zeros+ 0.02, yend = fit + 0.02, arrow = arrow(length = unit(0.1, "inches"), type = "closed"), color = "darkorchid4", size = 1)) +
  with(dat.fit.zeros[which.max(dat.fit.zeros$SI),], annotate(geom = "segment", x = zeros + 0.022, y = fit + 0.022, xend = zeros+ 0.021, yend = fit + 0.021, arrow = arrow(length = unit(0.08, "inches"), type = "closed"), color = "#FFE5DB", size = 0.01)) + 
  with(dat.fit.zeros[which.max(dat.fit.zeros$SI),], annotate(geom = "label", x = zeros + 0.11, y = fit + 0.11, color = "darkorchid4", label = substring(sprintf("%.3f", SI), 2), fill = "#FFF1EC"))


siplot
