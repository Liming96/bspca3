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
bspca3 <- function(x,k=ncol(x),type=c("Gram","predictor"),lambda=10000,ncomp=min(dim(x)),
                   center.=TRUE,scale.=FALSE,bess_tol=1e-3,bess_maxiter=100){
  if(length(k) == 1){
    k <- rep(k,ncomp)
  }else{
    ncomp <- length(k)
  } #确定主成分个数ncomp和稀疏度k

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
  ) #确定x是数据矩阵还是协方差矩阵

  cen <- attr(X, "scaled:center")
  sc <- attr(X, "scaled:scale")

  Xp <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
  Xp <- X

  svdobj<-svd(X) # 奇异值分解
  v<-svdobj$v # classical pca主成分,若x是协方差矩阵，v不变
  totalvariance<-sum((svdobj$d)^2) # x的总方差
  alpha<-as.matrix(v[,1:ncomp,drop=FALSE])

  W <- matrix(0,p, ncomp) # 初始化主成分矩阵
  sdev <- rep(0,ncomp)#附加解释标准差

  ccs <- seq(ncomp)#产生一个长度为ncomp的向量
  for(cc in ccs){
    if(type=="predictor"){
      #res <- pre_spca(Xp,sparsity=k[cc],lambda,Q[,seq_len(cc-1),drop=FALSE],
      #bess_tol,bess_maxiter)#计算第cc个主成分
      res <- pre_spca3(Xp,sparsity=k[cc],lambda,bess_tol,bess_maxiter)#计算第cc个主成分
    }else{
      res <- cov_spca3(Xp,sparsity=k[cc],lambda,bess_tol,bess_maxiter)#计算第cc个主成分
    }
    w <- res$w
    W[,cc] <- w #第cc个主成分为w


    if(type=="predictor"){
      Xp <- Xp-X%*%w%*%t(w)#对数据矩阵进行deflate
    }else{
      temp <- (t(w)%*%Xp%*%w)[1]
      Sp <- Xp-temp*w%*%t(w)
    }#根据x是数据矩阵还是协方差矩阵进行deflate

    if(type=="predictor"){
      if(cc < ncomp && all(abs(Xp)<1e-14)){
        W <- W[,1:cc,drop=FALSE]
        break
      }
    }else if(cc < ncomp && all(abs(Sp)<1e-14)){
      W <- W[,1:cc,drop=FALSE]
      break
    }
  }

  Z <- Xp%*%W
  qrZ <- qr(Z)
  RqrZ <- qr.R(qrZ)
  sdev <- diag(RqrZ)^2

  bspc <- list(sdev=sdev,rotation=W,center = if(is.null(cen)) FALSE else cen,
               scale = if(is.null(sc)) FALSE else sc,
               X = Xp)
  class(bspc) <- c("bsprcomp", "prcomp")
  return(bspc)
}


pre_spca3 <- function(Xp,sparsity,lambda,bess_tol=1e-3,bess_maxiter=100){
  n <- nrow(Xp)
  p <- ncol(Xp)
  residual <- rep(1,p)#定义标准化残差residual
  sacrifice <- rep(1,p)#定义牺牲向量sacrifice
  active <- rep(1,p)#定义活跃集和非活跃集

  #初始化
  w <- rep(1/p,p)#随机产生长度为p的正态分布数据，对w进行初始化
  for(j in 1:p){
    residual[j] <- (lambda*w[j]-t(Xp[,j])%*%Xp%*%w)/(t(Xp[,j])%*%Xp[,j]-lambda)
  }#初始化残差向量
  for(j in 1:p){
    sacrifice[j] <- (lambda-t(Xp[,j])%*%Xp[,j])*(w[j]+residual[j])^2
  }#初始化residual
  m <- sort(sacrifice,decreasing=T)
  active <- sacrifice >= m[sparsity]#初始化活跃集和非活跃集
  #更新活跃集和非活跃集的载荷值
  w[!active] <- 0
  svd_act <- svd(Xp[,active,drop=FALSE])
  w[active] <- svd_act$v[,1]

  obj_old <- -Inf #最开始目标函数的值
  ii <- 0

  while(ii <= bess_maxiter){
    wt <- w
    res <- Exchange(Xp,w,sparsity,sacrifice,residual,active,lambda)
    w <- res$w
    sacrifice <- res$sacrifice
    active <- res$active

    #标准化
    w <- w/sqrt(sum(w^2))#对w进行标准化

    obj <- -t(w)%*%t(Xp)%*%Xp%*%w+lambda*(t(w)%*%w-1)
    if(abs((obj-obj_old)/obj) < bess_tol){
      break
    }else{
      ii <- ii + 1
    }
    obj_old <- obj
  }
  w <- w/sqrt(sum(w^2))#对w进行标准化
  return(list(w=w,obj=obj))
}

cov_spca3 <- function(Sp,sparsity,lambda,bess_tol=1e-3,bess_maxiter=100){

  p <- dim(Sp)[2]
  residual <- rep(1,p)#定义标准化残差residual
  sacrifice <- rep(1,p)#定义牺牲向量sacrifice
  active <- rep(1,p)#定义活跃集和非活跃集

  w <- rep(1/p,p)#随机产生长度为p的正态分布数据，对w进行初始化
  for(j in 1:p){
    residual[j] <- (lambda*w[j]-t(Sp[,j])%*%w)/(Sp[j,j] - lambda)
  }#初始化残差向量
  for(j in 1:p){
    sacrifice[j] <- (lambda-Sp[j,j])*(w[j]+residual[j])^2
  }#初始化sacrifice
  m <- sort(sacrifice,decreasing=T)
  active <- sacrifice >= m[sparsity]#确定活跃集和非活跃集
  ##更新活跃集和非活跃集的载荷值
  w[!active] <- 0
  svd_act <- svd(Sp[,active,drop=FALSE])
  w[active] <- svd_act$v[,1]

  obj_old <- -Inf #最开始目标函数的值
  ii <- 0

  while(ii <= bess_maxiter){
    wt <- w
    res <- Exchange_cov(Sp,w,sparsity,sacrifice,residual,active,lambda)
    w <- res$w
    sacrifice <- res$sacrifice
    active <- res$active

    #标准化
    w <- w/sqrt(sum(w^2))#对w进行标准化

    obj <- -t(w)%*%Sp%*%w+lambda*(t(w)%*%w-1)
    if(abs((obj-obj_old)/obj) < bess_tol){
      break
    }else{
      ii <- ii + 1
    }
    obj_old <- obj
  }
  w <- w/sqrt(sum(w^2))#对w进行标准化
  return(list(w=w,obj=obj))
}


Exchange <- function(Xp,w,sparsity,sacrifice,residual,active,lambda){
  p <- ncol(Xp)
  obj <- -t(w)%*%t(Xp)%*%Xp%*%w+lambda*(t(w)%*%w-1)#计算目标函数值
  w_new <- w
  c <- 1
  m <- sort(sacrifice,decreasing=T)
  while(c < min(sparsity,p-sparsity)){
    #确定用来交换的活跃子集和非活跃子集
    Exch_active <- active & (sacrifice < m[sparsity-c])
    Exch_inactive <- !active & (sacrifice >= m[sparsity+c])
    #交换子集
    active_new <- active & !Exch_active | Exch_inactive
    #更新活跃集和非活跃集的载荷值
    w_new[!active_new] <- 0
    x <- Xp[,active_new,drop=FALSE]
    svd_act <- svd(x)
    w_new[active_new] <- svd_act$v[,1]
    #计算目标函数值
    #browser()
    obj_new <- -t(w_new)%*%t(Xp)%*%Xp%*%w_new+lambda*(t(w_new)%*%w_new-1)

    if(obj > obj_new){
      obj <- obj_new
      active <- active_new
      w <- w_new
    }else{
      c <- c + 1
    }
  }

  #更新残差向量residual
  for(j in 1:p){
    if(active[j]){
      residual[j] <- 0
    }else{
      residual[j] <- (lambda*w[j]-t(Xp[,j])%*%Xp%*%w)/(t(Xp[,j])%*%Xp[,j]-lambda)
    }
  }
  #计算sacrifice
  for(j in 1:p){
    sacrifice[j] <- (lambda-t(Xp[,j])%*%Xp[,j])*(w[j]+residual[j])^2
  }
  return(list(w=w,sacrifice=sacrifice,active=active,obj=obj))
}



Exchange_cov <- function(Sp,w,sparsity,sacrifice,residual,active,lambda){
  p <- ncol(Sp)
  obj <- -t(w)%*%Sp%*%w+lambda*(t(w)%*%w-1)#计算目标函数值
  w_new <- w
  c <- 1
  m <- sort(sacrifice,decreasing=T)
  while(c < min(sparsity,p-sparsity)){
    #确定用来交换的活跃子集和非活跃子集
    Exch_active <- active & (sacrifice < m[sparsity-c])
    Exch_inactive <- !active & (sacrifice >= m[sparsity+c])
    #交换子集
    active_new <- active & !Exch_active | Exch_inactive
    #更新活跃集和非活跃集的载荷值
    w_new[!active_new] <- 0
    s <- Sp[,active_new,drop=FALSE]
    svd_act <- svd(s)
    w_new[active_new] <- svd_act$v[,1]
    #计算目标函数值
    #browser()
    obj_new <- -t(w_new)%*%Sp%*%w_new+lambda*(t(w_new)%*%w_new-1)

    if(obj > obj_new){
      obj <- obj_new
      active <- active_new
      w <- w_new
    }else{
      c <- c + 1
    }
  }

  #更新残差向量residual
  for(j in 1:p){
    if(active[j]){
      residual[j] <- 0
    }else{
      residual[j] <- (lambda*w[j]-t(Sp[,j])%*%w)/(Sp[j,j] - lambda)
    }
  }
  #计算sacrifice
  for(j in 1:p){
    sacrifice[j] <- (lambda-Sp[j,j])*(w[j]+residual[j])^2
  }
  return(list(w=w,sacrifice=sacrifice,active=active,obj=obj))
}

