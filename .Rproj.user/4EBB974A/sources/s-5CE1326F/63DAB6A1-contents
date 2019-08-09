#' Title
#'
#' @param x Matrix.
#' @param k Vector or number.
#' @param type "Gram" or "predictor".
#' @param lambda Lagrange multiplier.
#' @param ncomp Number of PCs needed.
#' @param center. If TRUE, centering x.
#' @param scale. If TRUE, scaling x.
#' @param bess_tol Stop condition.
#' @param bess_maxiter Maximal iterate number.
#'
#' @return A dataframe.
#' @examples
#' \dontrun{
#' data("pitprops")
#' result <- bsPCA(pitprops,k=7,type="Gram")
#' }
#' @export
bsPCA <- function(x,k=ncol(x),type=c("Gram","predictor"),lambda=10000,ncomp=min(dim(x)),
                  center.=TRUE,scale.=FALSE,bess_tol=1e-3,bess_maxiter=100){
  if(length(k) == 1){
    k <- rep(k,ncomp)
  }else{
    ncomp <- length(k)
  } #Determine the number of PCs needed 'ncomp' and sparsity 'k'

  p <- ncol(x)
  X <- switch(type,
              predictor = {
                n<-dim(x)[1]
                p<-dim(x)[2]
                if (n/p>=100){
                  cat("You may wish to restart and use a more efficient way \n")
                  cat("let the argument x be the sample covariance/correlation matrix and set type=Gram \n")
                }
                X<-scale(x,center=center.,scale=scale.)
              },
              Gram = {
                n<-dim(x)[1]
                x_temp <- scale(x,center=center.,scale=scale.)
                X <- t(x_temp)%*%x_temp/(n-1)
              }
  ) #Deterimine X is data matrix or covariance matrix

  cen <- attr(X, "scaled:center")
  sc <- attr(X, "scaled:scale")

  Xp <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
  Xp <- X

  svdobj <- svd(X) # singualr value decomposition
  v <- svdobj$v # the PC of classical pca,If x is covariance matrix锛寁 is the same.
  totalvariance <- sum((svdobj$d)^2) # the totalvariance of x
  alpha <- as.matrix(v[,1:ncomp,drop=FALSE])

  W <- matrix(0,p,ncomp) # initialize PCs matrix
  sdev <- rep(0,ncomp)# additional explanation standard deviation

  ccs <- seq(ncomp)# Generate a vector of length ncomp
  for(cc in ccs){
    if(type=="predictor"){
      res <- spca1(Xp,sparsity=k[cc],lambda,bess_tol,bess_maxiter)
      w <- res$w
      W[,cc] <- w
      sdev[cc] <- t(Xp%*%w)%*%(Xp%*%w)
    }else{
      res <- spca2(Xp,sparsity=k[cc],lambda,bess_tol,bess_maxiter)
      w <- res$w
      W[,cc] <- w
    }

    # deflate the data matrix or covariance matrix
    if(type=="predictor"){
      Xp <- Xp-X%*%w%*%t(w)
    }else{
      temp <- (t(w)%*%Xp%*%w)[1]
      Xp <- Xp-temp*w%*%t(w)
    }

    if(type=="predictor"){
      if(cc < ncomp && all(abs(Xp)<1e-14)){
        W <- W[,1:cc,drop=FALSE]
        break
      }
    }else if(cc < ncomp && all(abs(Xp)<1e-14)){
      W <- W[,1:cc,drop=FALSE]
      break
    }
  }

  Z <- Xp%*%W
  qrZ <- qr(Z)
  RqrZ <- qr.R(qrZ)
  sdev <- diag(RqrZ)^2

  bspc <- list(sdev=sdev,rotation=W,X = Xp)
  return(bspc)
}




spca1 <- function(Xp,sparsity,lambda,bess_tol=1e-3,bess_maxiter=100){
  n <- nrow(Xp)
  p <- ncol(Xp)
  residual <- rep(1,p)#Define residual
  sacrifice <- rep(1,p)#Define sacrifice
  active <- rep(1,p)#Define active set and inactive set

  #Initialize
  w <- rep(1/p,p)
  for(j in 1:p){
    residual[j] <- (lambda*w[j]-t(Xp[,j])%*%Xp%*%w)/(n-lambda)
  }
  for(j in 1:p){
    sacrifice[j] <- (lambda-n)*(w[j]+residual[j])^2
  }
  m <- sort(sacrifice,decreasing=T)
  active <- sacrifice >= m[sparsity]
  ##update active set and inactive set
  w[!active] <- 0
  svd_act <- svd(Xp[,active,drop=FALSE])
  w[active] <- svd_act$v[,1]

  obj_old <- 0
  ii <- 0

  while(ii <= bess_maxiter){
    wt <- w
    res <- Exchange1(Xp,w,sparsity,sacrifice,residual,active,lambda)
    w <- res$w
    sacrifice <- res$sacrifice
    active <- res$active

    obj <- -t(w)%*%t(Xp)%*%Xp%*%w+lambda*(t(w)%*%w-1)
    if(abs((obj-obj_old)/(obj+0.01)) < bess_tol){
      break
    }else{
      ii <- ii + 1
    }
    obj_old <- obj
  }
  w <- w/sqrt(sum(w^2))#normalize
  return(list(w=w,obj=obj))
}





Exchange1 <- function(Xp,w,sparsity,sacrifice,residual,active,lambda){
  n <- nrow(Xp)
  p <- ncol(Xp)
  obj <- -t(w)%*%t(Xp)%*%Xp%*%w+lambda*(t(w)%*%w-1)
  w_new <- w
  c <- 1
  m <- sort(sacrifice,decreasing=T)
  while(c < min(sparsity,p-sparsity)){
    ##Determine the active subset and inactive subset for exchange
    Exch_active <- active & (sacrifice < m[sparsity-c])
    Exch_inactive <- !active & (sacrifice >= m[sparsity+c])
    ##Exchange the subset
    active_new <- active & !Exch_active | Exch_inactive
    ##Update the active sset and inactive set
    w_new[!active_new] <- 0
    x <- Xp[,active_new,drop=FALSE]
    svd_act <- svd(x)
    w_new[active_new] <- svd_act$v[,1]
    ##Compute the objective function value
    obj_new <- -t(w_new)%*%t(Xp)%*%Xp%*%w_new+lambda*(t(w_new)%*%w_new-1)

    if(obj > obj_new){
      obj <- obj_new
      active <- active_new
      w <- w_new
    }
    c <- c+1
  }

  ##Update residual
  for(j in 1:p){
    if(active[j]){
      residual[j] <- 0
    }else{
      residual[j] <- (lambda*w[j]-t(Xp[,j])%*%Xp%*%w)/(n-lambda)
    }
  }
  ##Update sacrifice
  for(j in 1:p){
    sacrifice[j] <- (lambda-n)*(w[j]+residual[j])^2
  }
  return(list(w=w,sacrifice=sacrifice,active=active,obj=obj))
}




spca2 <- function(Sp,sparsity,lambda,bess_tol=1e-3,
                  bess_maxiter=100){

  p <- dim(Sp)[2]
  residual <- rep(1,p)#Define residual
  sacrifice <- rep(1,p)#Define sacrifice
  active <- rep(1,p)#Define active set and inactive set

  w <- rep(1/p,p)
  for(j in 1:p){
    residual[j] <- (lambda*w[j]-t(Sp[,j])%*%w)/(Sp[j,j] - lambda)
  }
  for(j in 1:p){
    sacrifice[j] <- (lambda-Sp[j,j])*(w[j]+residual[j])^2
  }
  m <- sort(sacrifice,decreasing=T)
  active <- sacrifice >= m[sparsity]
  ##update active set and inactive set
  w[!active] <- 0
  svd_act <- svd(Sp[,active,drop=FALSE])
  w[active] <- svd_act$v[,1]

  obj_old <- 0
  ii <- 0

  while(ii <= bess_maxiter){
    wt <- w
    res <- Exchange2(Sp,w,sparsity,sacrifice,residual,active,lambda)
    w <- res$w
    sacrifice <- res$sacrifice
    active <- res$active

    obj <- -t(w)%*%Sp%*%w+lambda*(t(w)%*%w-1)
    if(abs((obj-obj_old)/(obj+0.01)) < bess_tol){
      break
    }else{
      ii <- ii + 1
    }
    obj_old <- obj
  }
  w <- w/sqrt(sum(w^2)) #normalize
  return(list(w=w,obj=obj))
}


Exchange2 <- function(Sp,w,sparsity,sacrifice,residual,active,lambda){
  p <- ncol(Sp)
  obj <- -t(w)%*%Sp%*%w+lambda*(t(w)%*%w-1)
  w_new <- w
  c <- 1
  m <- sort(sacrifice,decreasing=T)
  while(c < min(sparsity,p-sparsity)){
    ##Determine the active subset and inactive subset for exchange
    Exch_active <- active & (sacrifice < m[sparsity-c])
    Exch_inactive <- !active & (sacrifice >= m[sparsity+c])
    ##Exchange the subset
    active_new <- active & !Exch_active | Exch_inactive
    ##Update the active sset and inactive set
    w_new[!active_new] <- 0
    s <- Sp[,active_new,drop=FALSE]
    svd_act <- svd(s)
    w_new[active_new] <- svd_act$v[,1]
    ##Compute the objective function value
    obj_new <- -t(w_new)%*%Sp%*%w_new+lambda*(t(w_new)%*%w_new-1)

    if(obj > obj_new){
      obj <- obj_new
      active <- active_new
      w <- w_new
    }
    c <- c + 1
  }

  ##Update residual
  for(j in 1:p){
    if(active[j]){
      residual[j] <- 0
    }else{
      residual[j] <- (lambda*w[j]-t(Sp[,j])%*%w)/(Sp[j,j] - lambda)
    }
  }
  ##Update sacrifice
  for(j in 1:p){
    sacrifice[j] <- (lambda-Sp[j,j])*(w[j]+residual[j])^2
  }
  return(list(w=w,sacrifice=sacrifice,active=active,obj=obj))
}
