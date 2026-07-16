#######################################################################################
#INPUT: CDM(X);  d by n matrix X as X=(x_1,...,x_n),
#where n (> 3) is the sample size and d (> 1) is the dimension.

#OUTPUT: 
#values[j]; The estimator of the j-th eigenvalue by the CD method. 
#vectors[, j]; The estimator of the j-th eigenvector by the CD method.
#scores[, j];  The estimator of the j-th PC scores for x_1,...,x_n by the CD method.
#######################################################################################

CDM <- function(X, random='False'){
  d <- dim(X)[1]
  n <- dim(X)[2]
  n1 <- as.integer(ceiling(n/2))
  n2 <- n - n1
  n12 <- list(n1, n2)
  r <- min(n2-1, d)
  
  if (random=='False'){
    index <- c(1:n)
    Xcdm <- list(X[, 1:n1], X[, (n1+1):n])
  } else if (random=='True'){
    pi <- matrix(1/n, 1, n)
    index <- sample(c(1:n), n, replace=FALSE, pi)
    Xcdm <- list(X[, index[1:n1]], X[, index[(n1+1):n]])
  }
  Mean <- list(apply(Xcdm[[1]], 1, mean), apply(Xcdm[[2]], 1, mean))
  for (i in 1:2){
    Xcdm[[i]] <- sweep(Xcdm[[i]], 1, Mean[[i]], '-')
  }
  
  Sd <- t(Xcdm[[1]]) %*% Xcdm[[2]] / sqrt((n1 - 1) * (n2 - 1))
  eig <- svd(Sd)
  cdmval <- eig$d[1:r]
  crossvec <- list(eig$u, eig$v)
  
  cdmvec <- matrix(0, d, r)
  cdmscore <- matrix(0, n, r)
  
  for (i in 1:r){
    crossvec[[2]][, i] <- sign(as.numeric(crossvec[[1]][, i] %*% t(Xcdm[[1]]) %*% Xcdm[[2]] %*% crossvec[[2]][, i])) * crossvec[[2]][, i]
    
    #eigenvalue and eigenvector
    h <- list(0, 0)
    for (j in 1:2){
      h[[j]] <- Xcdm[[j]] %*% crossvec[[j]][, i] / sqrt(cdmval[i] * (n12[[j]] - 1))
    }
    
    cdmvec[, i] <- as.vector(h[[1]] + h[[2]]) / 2
    cdmvec[, i] <- cdmvec[, i] / norm(cdmvec[, i], type='2')
    
    
    #score
    r_score <- list(numeric(n1), numeric(n2))
    for (j in 1:2){
      for (k in 1:n12[[j]]){
        r_score[[j]][k] <- crossvec[[j]][k, i] * sqrt(n12[[j]] * cdmval[i])
      }
    }
    for (k in 1:n){
      cdmscore[index[k], i] <- c(r_score[[1]], r_score[[2]])[k]
    }
  }
  
  return(list(values=cdmval, vectors=cdmvec, scores=cdmscore))
}
