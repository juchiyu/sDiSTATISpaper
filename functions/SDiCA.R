#' Sparse Discriminant Correspondence Analysis
#'
#' @param Y A data frame or matrix with factors
#' @param design A vector that describes the row membership of Y
#' @param c1 The radius of l1 ball (constraint). It can take values from 1 to the square root of the number of row. The closer to 1 it is, the more sparsified are the observations.
#' @param c2 The radius of l2 ball (constraint). It can take values from 1 to the square root of the number of columns. The closer to 1 it is, the more sparsified are the variables.
#' @param n The number of components we want to obtain
#' @param meth Method to create the sparse singular vectors, by default "cgsvd" (can also be "gpmd")
#' @param init Method to initialize the singular vectors, by default "rand" (can also be "svd")
#' @param v.partition If true, creates a group constraint based on the ariables partition (by default False)
#' @param Gcol A group partition of the data columns (categories)
#' @param Grow A group partition of the groups
#' @param order If TRUE, the dimensions will be decreasingly ordered by their eigen value.
#' @param double.centering
#' @param row.w A vector containing the weight of the observations (default is a vector of 1). It must be of the same size as the number of observation.
#'
#' @return Returns the sparse MCA for the given constraints.
#' @export
#'
#' @examples
#'data(cheese)
#'SDiCA(cheese, c1 = sqrt(nrow(cheese))/2, c2 = sqrt(ncol(cheese))/2, n = 4)

SDiCA <- function(Y, 
                  design, 
                  c1, 
                  c2, 
                  n = 5, 
                  meth ='cgsvd', 
                  init = "svd", 
                  v.partition = F, 
                  Grow = NULL,
                  Gcol = NULL,
                  order = T, 
                  row.w = NULL,
                  double.centering = T,
                  orth.first = T,
                  itermax.pi = 1000,
                  itermax.pocs = 1000) {

  if(sum(is.na(Y))>0) stop("Error. Missing values.")

  X.ind <- tab_disjonctif(Y)

  if (length(design) != nrow(X.ind)) stop("Error. The length of the design vector must be equal to the number of rows of the data.")

  X <- aggregate(X.ind, by = list(design), sum)
  rownames(X) <- X$Group.1
  X <- X[,-1]

  if (n > min(ncol(X), nrow(X))) n <- min(ncol(X), nrow(X))

  if (length(c1) == 1) {
    c1 <- rep (c1, times = n)
  } else if (length(c1) != n)
    stop("Error, the length of c1 must be equal to the number of dimension.")

  if (length(c2) == 1) {
    c2 <- rep (c2, times = n)
  } else if (length(c2) != n)
    stop("Error, the length of c2 must be equal to the number of dimension.")

  if (meth =='cgsvd') {
    res <- cgsvd(Y = Y, X = X, c1 = c1, c2 = c2, R = n, init = init, v.partition = v.partition, Grow = Grow, Gcol = Gcol, row.w = row.w, double.centering = double.centering, orth.first = orth.first,itermax.pi = itermax.pi, itermax.pocs = itermax.pocs)
  }else if (meth =='gpmd') {
    warning("Caution, only the first value of c1 and c2 are considered in PMD.")
    res <- gpmd(Y = Y, X = X, K = n, sumabsu = c1[1], sumabsv = c2[1])
  }

  if ((sum(sort(res$D, d = T) != res$D) > 0) & order) {
    new.order <- order(res$D, decreasing = T)
    res$D <- sort(res$D, d = T)
    res$P <- res$P[, new.order]
    res$Q <- res$Q[, new.order]

  }

  if(is.null(Gcol)) {Gcol <- partition_variables(Y)}
  Itot <- 1/length(Gcol)*sum(sapply(1:length(Gcol), function(i) {length(Gcol[[i]])-1}))

  eig <- as.data.frame(cbind(dim = seq(from = 1, by = 1, length = length(res$D)),
                             eigenvalue = res$D,
                             percentageOfVariance = res$D / Itot * 100,
                             cumulatedPercentageOfVariance = cumsum(res$D / Itot * 100)
                             )
                       )

  # groups
  F <- diag(as.numeric(res$r)^(-1)) %*% res$P %*% diag(sqrt(res$D))[,1:n] #diag(as.numeric(res$r)^(-1/2)) %*% res$P %*% diag(sqrt(res$D))[,1:n]
  F2 <- F^2
  contrib <- diag(as.numeric(res$r)) %*% F2 %*% diag(1/res$D[1:n]) #(F2) %*% diag(1/colSums(F2)) * 100 #
  cos2 <- t(t(F2)%*%diag(1/rowSums(F2)))
  
  col <- paste("dim", seq(from = 1, by = 1, length = n))
  colnames(F) <- colnames(contrib) <- colnames (cos2) <- col
  rownames(F) <- rownames(contrib) <- rownames (cos2) <- rownames (X)

  grp <- list(coord = F, contrib = contrib, cos2 = cos2)
  
  #categories
  G <- diag(as.numeric(res$c)^(-1)) %*% res$Q %*% diag(sqrt(res$D))[,1:n] #diag(as.numeric(res$c)^(-1/2)) %*% res$Q %*% diag(sqrt(res$D))[,1:n]
  G2 <- G^2
  contrib <- diag(as.numeric(res$c)) %*% G2 #contribution absolue
  cos2 <- t(t(G2)%*%diag(1/rowSums(G2)))
  
  partition <- partition_variables(Y)
  
  eta2 <- c()
  for (j in 1:n) {
    eta2 <- cbind(eta2, sapply(1:ncol(Y), function(i) {sum(contrib[partition[[i]],j]) }))
  }
  
  colnames(G) <- colnames(contrib) <- colnames (cos2) <- colnames(eta2) <- col
  rownames(G) <- rownames(contrib) <- rownames (cos2) <- colnames (X)
  rownames(eta2) <- colnames(Y)
   
  var <- list(coord = G, contrib = contrib %*% diag(1/res$D[1:n]), cos2 = cos2, eta2 = eta2)
  
  
  return(list (cgsvd = res,
               eig = eig,
               var = var,
               grp = grp,
               other =list(Xinit = Y,
                           Xgroup = X,
                           Xdisj = X.ind,
                           Gcol = Gcol,
                           Grow = Grow)))
}
