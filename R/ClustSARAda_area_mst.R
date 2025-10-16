#### clustered area level model with SAR(1) model fixed tau and rho
#' @import sae
#' @importFrom dplyr filter select mutate union groups row_number
#' @importFrom igraph as_data_frame mst E graph_from_data_frame components
#' @importFrom igraph "E<-"
#' @importFrom magrittr %>%
#' @importClassesFrom Matrix dgCMatrix
#' @import Matrix
#' @import spdep
#' @return
#' @export


ClustSARAda_area_mst <- function(dom, y, x, index_d, vardir, graph0,
                              betam0, tauf, rhof,
                              adjW, nu = 1, gam = 3, lam = 0.5,
                              lam0 = 0.001, maxiter= 500, tol = 1e-4)
{

  n0 <- length(y)
  uindex <- unique(dom)
  m <- length(uindex)
  ns <- as.numeric(table(dom)) # number of replicates
  nbar <- mean(ns)

  len_d <- length(index_d)
  Ip <- diag(1,len_d,len_d)
  Im <- diag(m)

  ### mst and H matrix, the graph is based on nr###
  E(graph0)$weight <- apply(incidentGen(graph0, m) %*% as.matrix(betam0), 1, norm, "2")
  graph_mst <- mst(graph0)
  Hm <- incidentGen(graph_mst, m)
  #Hm <- rbind(matrix(1/sqrt(m),1,m), Hm)

  Am <- Hm %x% Ip
  AtA <- (t(Hm) %*% Hm) %x% Ip


  #### fixed design information ########

  x_d <- x[,index_d, drop = FALSE]
  Xm <- matrix(0, n0, m*len_d)
  for(i in 1:m)
  {
    Xm[dom == uindex[i],(len_d*(i-1) + 1) : (len_d*i)] <- x_d[dom == uindex[i],]
  }

  beta_cur <- c(t(betam0))

  ##### two cases ###
  ## case 1
  if(len_d == ncol(x))
  {
    ##### initial of others
    deltam <- Hm %*% betam0
    vm <-  matrix(0, (m-1) , len_d)


    Cm <- (Im - rhof* t(adjW)) %*% (Im - rhof*adjW)
    Sm <- diag(vardir) + tauf * solve(Cm)
    Sm_inv <- solve(Sm)
    Wm <- 1/m*(Sm_inv*(1/nbar))
    Xmt <- t(Xm) %*% Wm

    nu <- mean(diag(Sm_inv))
    #nu <- mean(1/vardir)
    ### iterations
    for(j in 1:maxiter)
    {
      # update beta and eta


      beta_new <- solve(Xmt %*% Xm + nu * AtA) %*% (Xmt %*% y + nu*c(as.matrix(t(deltam -  vm/nu) %*% Hm)))
      betam_new <- matrix(beta_new, m, len_d, byrow = TRUE)

      Dbeta <- Hm %*% betam_new
      psim <- Dbeta + vm/nu
      deltam_new <- t(do.call("cbind",lapply(1:nrow(psim),function(xx) mcp(psim[xx,],lam,nu,gam))))

      vm <- vm + nu * (Dbeta - deltam_new)

      diff_norm <- sqrt(sum((Dbeta- deltam_new)^2))
      diff2 <- max(abs(Dbeta- deltam_new))

      deltam <- deltam_new
      beta_cur <- beta_new

      if(diff2 < tol)
      {break}
    }

    resids_new <- as.numeric(y -Xm %*% beta_new)


    cluster_est <- get_sp_cluster2(graph_mst, m, len_d, betam_new)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)

    Cm <- (Im - rhof* t(adjW)) %*% (Im - rhof*adjW)
    Sm <- diag(vardir) + tauf * solve(Cm)

    muhat <- rowSums(x * beta_c[cluster_est,-1])
    ratio_est <- tauf * solve(Cm) %*% solve(Sm)
    area_est <- muhat + ratio_est %*%(y - muhat)

    out <- list(beta = betam_new,
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

    eta_cur <- ztz %*% (y - Xm %*% beta_cur)

    ##### initial of others
    deltam <- Am %*% beta_cur
    vm <-  matrix(0, (m-1) , len_d)

    Cm <- (Im - rhof* t(adjW)) %*% (Im - rhof*adjW)
    Sm <- diag(vardir) + tauf * solve(Cm)
    Sm_inv <- solve(Sm)
    Wm <- 1/m*(Sm_inv*(1/nbar))
    ztz <-  solve(t(z) %*% Wm %*% z)%*% t(z) %*% Wm
    Qz <- Wm - Wm %*% z%*% ztz

    nu <- mean(diag(Sm_inv))
    #nu <- mean(1/vardir)


    ### iterations
    for(j in 1:maxiter)
    {
      # update beta and eta

      Xinv <- solve(t(Xm) %*%Qz%*% Xm + nu * AtA)
      beta_new <- Xinv %*% (t(Xm) %*% Qz %*% y + nu*c(as.matrix(t(deltam -  vm/nu) %*% Hm)))

      betam_new <- matrix(beta_new, m, len_d, byrow = TRUE)
      eta_new <- as.numeric(ztz %*% (y - Xm%*%beta_new))


      Dbeta <- as.matrix(Hm %*% betam_new)
      psim <- Dbeta + vm/nu

      deltam_new <- t(do.call("cbind",lapply(1:nrow(psim),function(xx) mcp(psim[xx,],lam,nu,gam))))

      vm <- vm + nu * (Dbeta - deltam_new)

      diff_norm <- sqrt(sum((Dbeta - deltam_new)^2))
      diff2 <- max(abs(Dbeta- deltam_new))

      deltam <- deltam_new
      beta_cur <- beta_new

      if(diff2 < tol)
      {break}
    }

    resids_new <- as.numeric(y - z %*% eta_new - Xm %*% beta_new)

    cluster_est <- get_sp_cluster2(graph_mst, m, len_d, betam_new)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)


    Cm <- (Im - rhof* t(adjW)) %*% (Im - rhof*adjW)
    Sm <- diag(vardir) + tauf * solve(Cm)

    muhat <- z %*% eta_new +  rowSums(x_d * beta_c[cluster_est,-1])
    ratio_est <- tauf * solve(Cm) %*% solve(Sm)
    area_est <- muhat + ratio_est %*%(y - muhat)

    out <- list(beta = betam_new, eta = eta_new,
                betac = beta_c,
                cluster = cluster_est,
                muhat = muhat, area_est = area_est, vardir = vardir,
                deltam = deltam, niters = j)

  }

  out$beta00 <- betam0
  out$loss <- lossSARfun(tauf, rhof, vardir, resids_new, adjW)
  out$diff_max <- diff2


  return(out)

}

#' @export
ClustSARAda_area_mst_ada <- function(dom, y, x, index_d, vardir,
                                graph0, betam0, tauf, rhof, adjW,
                                nu = 1, gam = 3, lam = 0.5,
                                lam0 = 0.001, maxiter= 500,
                                tol = 1e-4)
{
  m <- length(unique(dom))
  res <- ClustSARAda_area_mst(dom, y, x, index_d, vardir,
                         graph0 = graph0, betam0 = betam0,
                         tauf = tauf, rhof = rhof,
                         adjW = adjW,
                          nu = nu, gam = gam, lam = lam,
                         lam0 = lam0, maxiter= maxiter,
                         tol = tol)

  graph_updated <- graph0
  E(graph_updated)$weight <- apply(incidentGen(graph0, m) %*% as.matrix(res$betac[res$cluster,-1]), 1, norm, "2")

  #betam0_up <- initialfun2(dom = dom, y = y, x = x, index_d, graph_updated)

  betam0_up <- as.matrix(res$betac[res$cluster,-1])
  res_updated <- ClustSARAda_area_mst(dom, y, x, index_d, vardir, graph0 = graph_updated,
                                 betam0 = betam0_up, tauf = tauf, rhof = rhof, adjW = adjW,
                                 nu = nu, gam = gam, lam = lam,
                                 lam0 = lam0, maxiter= maxiter,
                                 tol = tol)

  return(list(res = res_updated, graph_updated = graph_updated))
}



#' @export
ClustSARAda_area_mst_ada_bic <- function(dom, y, x, index_d, vardir,
                                    graph0, betam0, tauf, rhof,
                                    adjW, lamvec,
                                    nu = 1, gam = 3, lam0 = 0.001,
                                    maxiter= 500, tol = 1e-4)
{
  m <- length(unique(dom))
  nlam <- length(lamvec)
  len_d <- length(index_d)


  K_est <- rep(0, nlam)
  loss_value0 <- rep(0, nlam)
  beta_array <- array(0, dim = c(m, len_d, nlam))
  #SAR_est <- matrix(0, nlam, 2)

  for(j in 1:length(lamvec))
  {
    resj <- ClustSARAda_area_mst_ada(dom, y, x, index_d, vardir,
                                graph0 = graph0, betam0 = betam0,
                                tauf = tauf, rhof = rhof, adjW = adjW,
                                nu = nu, gam = gam, lam = lamvec[j],
                                lam0 = lam0, maxiter= maxiter,
                                tol = tol)
    resj_updated <- resj$res

    loss_value0[j] <- resj_updated$loss
    K_est[j] <- length(unique(resj_updated$cluster))
    beta_array[,,j] <- resj_updated$beta
   # SAR_est[j,] <- c(resj_updated$tauf, resj_updated$rhof)
  }



  res_df <- data.frame(lambda = lamvec,Khat = K_est,
                       loss0 = loss_value0)

  return(list(res_df = res_df, beta_array = beta_array))
}


## grad of function of tau and rho
# resids y - zeta - xbeta
# only for each dom one observation

lossSARfun <- function(tau, rho, vardir, resids, adjW)
{
  Im <- diag(1,length(vardir))
  adjW <- as.matrix(adjW)
  Cm <- (Im - rho* t(adjW)) %*% (Im - rho*adjW)
  Sm <- diag(vardir) + tau*solve(Cm)

  value1 <- determinant(Sm, logarithm = TRUE)$modulus
  value2 <- t(resids) %*% solve(Sm) %*% resids
  value <- value1 + value2
  return(value)
}


gradSARfun <- function(tau, rho, vardir, resids, adjW)
{
  Im <- diag(1,length(vardir))
  adjW <- as.matrix(adjW)
  adjWt <- t(adjW)
  Cm <- (Im - rho* adjWt) %*% (Im - rho*adjW)
  Sm <- diag(vardir) + tau*solve(Cm)
  Sm_inv <- solve(Sm)
  Cm_inv <- solve(Cm)

  value_tau1 <- sum(diag(Sm_inv %*% Cm_inv))
  value_tau2 <- t(resids) %*% Sm_inv %*% Cm_inv %*% Sm_inv %*% resids
  value_tau <- value_tau1 - value_tau2

  temp <- adjWt - 2*rho * adjWt %*% adjW + adjW
  value_rho1 <- sum(diag(tau* Sm_inv %*% Cm_inv %*% temp %*% Cm_inv))
  value_rho2 <- tau * t(resids) %*% Sm_inv %*% Cm_inv %*% temp %*% Cm_inv %*% Sm_inv %*% resids
  value_rho <- value_rho1 - value_rho2

  value <- c(value_tau, value_rho)
  return(value)
}



#' @export
Est_clustSARAda_area_mst <- function(dom, y, x, index_d, vardir,
                                  graph0, betam0, tauf = tauf, rhof = rhof, adjW,
                                  lamvec, nu = 1, gam = 3,
                                  lam0 = 0.001, maxiter= 500, tol = 1e-4)
{
  n0 <- length(y)
  len_d <- length(index_d)
  nz <- ncol(x) - len_d
  m <- length(unique(dom))

  res_bic <- ClustSARAda_area_mst_ada_bic(dom=dom, y, x, index_d, vardir,
                                       graph0= graph0, betam0 = betam0,
                                       tauf = tauf, rhof = rhof, adjW = adjW,
                                       lamvec = lamvec,
                                       nu = nu, gam = gam,
                                       lam0 = lam0, maxiter= maxiter,
                                       tol = tol)



  bic_df <- res_bic$res_df %>%
    mutate(bic = loss0 + log(n0)*(Khat*len_d + nz + 1),
           mbic = loss0 + log(m)*log(n0)*(Khat*len_d + nz +1))


  res_fm <- ClustSARAda_area_mst_ada(dom = dom, y, x, index_d, vardir,
                                  graph0 = graph0,betam0 = betam0,
                                  tauf = tauf, rhof = rhof, adjW = adjW,
                                  nu = nu, gam = gam, lam = lamvec[which.min(bic_df$mbic)],
                                  lam0 = lam0, maxiter= maxiter, tol = tol)
  res_fm <- res_fm$res
  area_estm <- res_fm$area_est

  fit_fm <- NULL
  tryCatch({
    fit_fm <- fit_clustSAR_area(dom = dom, y = y, x = x, index_d = index_d,
                                vardir = vardir, adjW = adjW, cluster = res_fm$cluster)
    area_estm <- fit_fm$area_est_ada
  }, error = function(e) {
    area_estm <- res_fm$area_est
  })

  out <- list(bic_df = bic_df,
              resfm = res_fm,
              cluster_estm = res_fm$cluster,
              area_estm = area_estm,
              refit = fit_fm,
              out_bic = res_bic)

  return(out)

}









