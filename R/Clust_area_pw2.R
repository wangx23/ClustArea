#' clustered area level model pairwise

#' @return
#' @export
#' @import sae
Clust_area_pw2 <- function(indexy, y, x, index_d, vardir,
                       nu = 1, gam = 3, lam = 0.5,
                       lam0 = 0.001, maxiter= 500, tol = 1e-4,
                       seed = 1234,
                       method = "kmeans", k0 = 10)
{

  n0 <- length(y)
  uindex <- unique(indexy)
  nr <- length(uindex)
  ns <- as.numeric(table(indexy)) # number of replicates
  nbar <- mean(ns)

  ## sparse function
  D <- matrix(0,nr*(nr-1)/2,nr)
  for(j in 1:(nr-1))
  {
    indexj <- (2*nr-j)*(j-1)/2
    indexvj <- indexj + (1:(nr-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):nr)] <- -1
  }

  len_d <- length(index_d)
  Ip <- diag(1,len_d,len_d)
  AtA <- t(D)%*%D %x% Ip

  x_d <- x[,index_d, drop = FALSE]

  # for X expand
  Xm <- matrix(0, n0, nr*len_d)
  for(i in 1:nr)
  {
    Xm[indexy == uindex[i],(len_d*(i-1) + 1) : (len_d*i)] <- x_d[indexy == uindex[i],]
  }

  #### use sae to initial sig2 ###

  fit0 <- eblupFH(y~0+x, vardir = vardir)


  ##### two cases ###
  if(len_d == ncol(x))
  {
    # no same part initial of beta

    Xty <- t(Xm) %*% y
    Xinv <- inverseR(indexy, x_d, lam0)
    beta00 <-  matrix(Xinv %*% Xty,ncol= len_d, byrow = TRUE)

   # if(min(ns) ==1)
    #{
      res00_ada1  <- adaptive_start(beta00, k0 = k0, method = method)
      cluster00 <- res00_ada1$cluster00
      res00_ada2 <- lm_clust_reg(dom = indexy, y = y, x = x, index_d = index_d,
                                 cluster = cluster00)
      beta00_ada <- matrix(res00_ada2, ncol = len_d, byrow = TRUE)[cluster00,,drop= FALSE]
      beta00 <- beta00_ada
    #}



    ##### initial of others
    deltam <- D %*% beta00
    vm <-  matrix(0, nr*(nr-1)/2 , len_d)

    beta_cur <- c(t(beta00))

    resids0 <- y - Xm %*% beta_cur

    lossj <- function(xx)
    {
      lossfun2(xx, vardir, resids0)
    }

    sig2_cur <- optimize(lossj, interval = c(0,mean(resids0^2)))$minimum



    ### iterations
    for(j in 1:maxiter)
    {
      # update beta and eta
      Wm <- diag(1/(nbar*(vardir + sig2_cur)))
      Xty <- t(Xm) %*% Wm %*% y

      wm <- diag(Wm)
      Xinv <- inverseR(indexy, x_d*sqrt(wm), nu = nu)

      beta_new <-  Xinv %*% (Xty + nu*c(t(deltam -  vm/nu) %*% D))
      betam_new <- matrix(beta_new, nr, len_d, byrow = TRUE)

      # update sigma_v^2
      resids <- y - Xm %*% beta_new

      lossj <- function(xx)
      {
        lossfun2(xx, vardir, resids)
      }

      sigv2_new <- optimize(lossj, interval = c(0,mean(resids^2)))$minimum


      Dbeta <- D %*% betam_new
      psim <- Dbeta + vm/nu

      deltam_new <- t(sapply(1:nrow(psim),function(xx) mcp(psim[xx,],lam,nu,gam)))
      vm <- vm + nu * (Dbeta - deltam_new)

      diff_norm <- sqrt(sum((Dbeta- deltam_new)^2))

      diff2 <- max(abs(Dbeta- deltam_new))
      deltam <- deltam_new
      sig2_cur <- sigv2_new
      if(diff2 < tol)
      {break}


    }

    cluster_est <- getgroup(t(deltam), n = nr)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)

    muhat <- rowSums(x * betam_new[indexy,])
    ratio_est <- sig2_cur/(vardir + sig2_cur)
    area_est <- ratio_est * y + (1-ratio_est)*muhat

    out <- list(beta = betam_new, sig2 = sig2_cur,
                betac = beta_c, cluster = cluster_est,
                muhat = muhat, area_est = area_est,
                vardir = vardir,deltam = deltam, niters = j)

  }
  if(len_d < ncol(x))
  {
    # having same part
    # initial of beta

    z <- x[,-index_d, drop = FALSE]

    ztz <- solve(t(z) %*% z)%*%t(z)
    Qz <- diag(1, n0, n0) - z%*% ztz

    Xty <- t(Xm)%*% Qz %*% y

    Xinv <- inverseR2(indexy, z, x_d, lam0)

    beta00 <-  matrix(Xinv %*% Xty,ncol= len_d, byrow = TRUE)

   # if(min(ns) ==1)
    #{
      res00_ada1  <- adaptive_start(beta00, k0 = k0, method = method)
      cluster00 <- res00_ada1$cluster00
      res00_ada2 <- lm_clust_reg(dom = indexy, y = y, x = x, index_d = index_d,
                                 cluster = cluster00)
      beta00_ada <- matrix(res00_ada2[-(1:ncol(z))], ncol = len_d, byrow = TRUE)[cluster00,,drop= FALSE]
      beta00 <- beta00_ada
    #}

    eta00 <- ztz %*% (y - Xm%*%beta00)


    ##### initial of others
    deltam <- D %*% beta00
    vm <-  matrix(0, nr*(nr-1)/2 , len_d)

    beta_cur <- c(t(beta00))

    resids0 <- y - z %*% eta00 - Xm %*% beta_cur
    lossj <- function(xx)
    {
      lossfun2(xx, vardir, resids0)
    }

    sig2_cur <- optimize(lossj, interval = c(0,mean(resids0^2)))$minimum

    ### iterations
    for(j in 1:maxiter)
    {
      # update beta and eta
      Wm <- diag(1/(vardir + sig2_cur)*(1/nbar))
      ztz <-  solve(t(z) %*% Wm %*% z)%*%t(z) %*% Wm
      Qz <- Wm - Wm %*% z%*% ztz

      wm <- diag(Wm)
      Xty <- t(Xm) %*% Qz %*% y

      #Xinv <- solve(XtX + nu*AtA)
      Xinv <- inverseR2(indexy, z*sqrt(wm), x_d*sqrt(wm), nu)

      beta_new <-  Xinv %*% (Xty + nu*c(t(deltam -  vm/nu) %*% D))
      betam_new <- matrix(beta_new, nr, len_d, byrow = TRUE)
      eta_new <- ztz %*% (y - Xm%*%beta_new)
      # update sigma_v^2
      resids <- y - z %*% eta_new - Xm %*% beta_new

      lossj <- function(xx)
      {
        lossfun2(xx, vardir, resids)
      }


      sigv2_new <- optimize(lossj, interval = c(0,mean(resids^2)))$minimum


      Dbeta <- D %*% betam_new
      psim <- Dbeta + vm/nu

      deltam_new <- t(do.call("cbind",lapply(1:nrow(psim),function(xx) mcp(psim[xx,],lam,nu,gam))))
      vm <- vm + nu * (Dbeta - deltam_new)

     diff_norm <- sqrt(sum((Dbeta- deltam_new)^2))

      diff2 <- max(abs(Dbeta- deltam_new))

      deltam <- deltam_new
      sig2_cur <- sigv2_new
      if(diff2 < tol)
      {break}
    }

    cluster_est <- getgroup(t(deltam), n = nr)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)

    muhat <- z %*% eta_new + rowSums(x_d * betam_new[indexy,])
    ratio_est <- sig2_cur/(vardir + sig2_cur)
    area_est <- ratio_est * y + (1-ratio_est)*muhat


    out <- list(beta = betam_new, eta = eta_new,
                betac = beta_c, sig2 = sig2_cur,
                cluster = cluster_est,
                muhat = muhat, area_est = area_est, vardir = vardir,
                deltam = deltam, niters = j)

  }

  out$beta00 <- beta00

  out$loss <- lossfun2(sigv2_new, vardir, resids)
  out$diff_max <- diff2
  out$diff_norm <- diff_norm

  return(out)

}


#' @export
Clust_area_pw_bic2 <- function(indexy, y, x, index_d, vardir, lamvec,
                              nu = 1, gam = 3, lam0 = 0.001, maxiter= 500, tol = 1e-4,
                              seed = 1234)
{
  K_est <- rep(0, length(lamvec))
  loss_value0 <- loss_value <- rep(0, length(lamvec))
  for(j in 1:length(lamvec))
  {
    resj <- Clust_area_pw2(indexy, y, x, index_d, vardir,
                          nu = nu, gam = gam, lam = lamvec[j],
                          lam0 = lam0, maxiter= maxiter, tol = tol,
                          seed = seed)
   # fitj <- fit_clust_area(dom = indexy, y = y, x = x, index_d = index_d,
  #                 vardir = vardir, cluster = resj$cluster)
    loss_value0[j] <- resj$loss
    #loss_value[j] <- -2*fitj$fit$goodness[1]
    K_est[j] <- length(unique(resj$cluster))
  }

  res_df <- data.frame(lambda = lamvec,Khat = K_est,
                       loss = loss_value, loss0 = loss_value0)
  return(res_df)
}



## grad of function of sigma_v^2
# resids y - zeta - xbeta
gradfun2 <- function(sigv2, vardir, resids)
{
  value1 <- sum(1/(vardir + sigv2))
  value2 <- sum(resids^2/(vardir + sigv2)^2)
  value <- value1 - value2
  return(value)
}

lossfun2 <- function(sigv2, vardir, resids)
{
  value1 <-  sum(log(vardir + sigv2))
  value2 <- sum(resids^2/(vardir + sigv2))
  value <- value1 + value2
  return(value)
}


