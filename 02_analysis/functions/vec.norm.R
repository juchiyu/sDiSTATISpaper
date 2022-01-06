vec.norm <- function(x){
  # equivalent to ||x||, where x is a vector
  as.numeric(sqrt(crossprod(x)))
}