createxyLabels.gen.digit <- function(x_axis = 1, y_axis = 2, lambda, tau, digit4lambda = 3, digit4tau = 0, axisName = "Dimension "){
  xyLabels = labs(x = createLabel.gen.digit(zeAxis = x_axis, lambda = lambda[x_axis], 
                                      tau = tau[x_axis], digit4lambda = digit4lambda, digit4tau = digit4tau, axisName), 
                  y = createLabel.gen.digit(zeAxis = y_axis, lambda = lambda[y_axis], tau = tau[y_axis], digit4lambda = digit4lambda, digit4tau = digit4tau, axisName))
  return(xyLabels)
}