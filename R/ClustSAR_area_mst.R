#### clustered area level model with SAR(1) model
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


ClustSAR_area_mst <- function(dom, y, x, index_d, vardir, graph0,
                              betam0, tau0 = NULL, rho0 = NULL,
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

    if(is.null(tau0) | is.null(rho0))
    {
      resids0 <- y - Xm %*% beta_cur

      lossj <- function(xx)
      {
        lossSARfun(xx[1],xx[2],vardir, resids0, adjW)
      }

      gradj <- function(xx)
      {
        gradSARfun(xx[1],xx[2], vardir, resids0, adjW)
      }

      resj <- optim(par = c(mean(resids0^2),0), fn = lossj, gr = gradj,
                    method = "L-BFGS-B",
                    lower = c(0,-0.999), upper = c(mean(y^2),0.999))
      tau_cur <- resj$par[1]
      rho_cur <- resj$par[2]
    }else{
      tau_cur <- tau0
      rho_cur <- rho0
    }

    ### iterations
    for(j in 1:maxiter)
    {
      # update beta and eta

      Cm <- (Im - rho_cur* t(adjW)) %*% (Im - rho_cur*adjW)
      Sm <- diag(vardir) + tau_cur * solve(Cm)
      Sm_inv <- solve(Sm)
      Wm <- 1/m*(Sm_inv*(1/nbar))
      Xmt <- t(Xm) %*% Wm

     # nu <- mean(diag(Sm_inv))
      nu <- mean(1/vardir)

      beta_new <- solve(Xmt %*% Xm + nu * AtA) %*% (Xmt %*% y + nu*c(as.matrix(t(deltam -  vm/nu) %*% Hm)))
      betam_new <- matrix(beta_new, m, len_d, byrow = TRUE)

      # update tau and rho
      resids_new <- as.numeric(y -Xm %*% beta_new)

      lossj <- function(xx)
      {
        lossSARfun(xx[1],xx[2],vardir, resids_new, adjW)
      }

      gradj <- function(xx)
      {
        gradSARfun(xx[1],xx[2], vardir, resids_new, adjW)
      }

      resj <- optim(par = c(mean(resids_new^2),0), fn = lossj, gr = gradj,
                    method = "L-BFGS-B",
                    lower = c(0,-0.999), upper = c(mean(y^2),0.999),
                    control = list(maxit = 20))
      tau_new <- resj$par[1]
      rho_new <- resj$par[2]


      Dbeta <- Hm %*% betam_new
      psim <- Dbeta + vm/nu
      deltam_new <- t(do.call("cbind",lapply(1:nrow(psim),function(xx) mcp(psim[xx,],lam,nu,gam))))

      vm <- vm + nu * (Dbeta - deltam_new)

      diff_norm <- sqrt(sum((Dbeta- deltam_new)^2))
      diff2 <- max(abs(Dbeta- deltam_new))

      deltam <- deltam_new
      tau_cur <- tau_new
      rho_cur <- rho_new
      beta_cur <- beta_new

      if(diff2 < tol)
      {break}
    }


    cluster_est <- get_sp_cluster2(graph_mst, m, len_d, betam_new)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)

    Cm <- (Im - rho_cur* t(adjW)) %*% (Im - rho_cur*adjW)
    Sm <- diag(vardir) + tau_cur * solve(Cm)

    muhat <- rowSums(x * beta_c[cluster_est,-1])
    ratio_est <- tau_cur * solve(Cm) %*% solve(Sm)
    area_est <- muhat + ratio_est %*%(y - muhat)

    out <- list(beta = betam_new, tau = tau_cur, rho = rho_cur,
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

    if(is.null(tau0) | is.null(rho0))
    {
      resids0 <- y - z %*% eta_cur - Xm %*% beta_cur

      lossj <- function(xx)
      {
        lossSARfun(xx[1],xx[2],vardir, resids0, adjW)
      }

      gradj <- function(xx)
      {
        gradSARfun(xx[1],xx[2], vardir, resids0, adjW)
      }

      resj <- optim(par = c(mean(resids0^2),0), fn = lossj, gr = gradj,
            method = "L-BFGS-B",
            lower = c(0,-0.999), upper = c(mean(y^2),0.999))
      tau_cur <- resj$par[1]
      rho_cur <- resj$par[2]
    }else{
      tau_cur <- tau0
      rho_cur <- rho0
    }

    ### iterations
    for(j in 1:maxiter)
    {
      # update beta and eta
      Cm <- (Im - rho_cur* t(adjW)) %*% (Im - rho_cur*adjW)
      Sm <- diag(vardir) + tau_cur * solve(Cm)
      Sm_inv <- solve(Sm)
      Wm <- 1/m*(Sm_inv*(1/nbar))
      ztz <-  solve(t(z) %*% Wm %*% z)%*% t(z) %*% Wm
      Qz <- Wm - Wm %*% z%*% ztz

      #nu <- mean(diag(Sm_inv))
      nu <- mean(1/vardir)

      Xinv <- solve(t(Xm) %*%Qz%*% Xm + nu * AtA)
      beta_new <- Xinv %*% (t(Xm) %*% Qz %*% y + nu*c(as.matrix(t(deltam -  vm/nu) %*% Hm)))

      betam_new <- matrix(beta_new, m, len_d, byrow = TRUE)
      eta_new <- as.numeric(ztz %*% (y - Xm%*%beta_new))


      # update tau and rho
      resids_new <- as.numeric(y - z %*% eta_new - Xm %*% beta_new)

      lossj <- function(xx)
      {
        lossSARfun(xx[1],xx[2],vardir, resids_new, adjW)
      }

      gradj <- function(xx)
      {
        gradSARfun(xx[1],xx[2], vardir, resids_new, adjW)
      }

      resj <- optim(par = c(mean(resids_new^2),0), fn = lossj, gr = gradj,
                    method = "L-BFGS-B",
                    lower = c(0,-0.999), upper = c(mean(y^2),0.999),
                    control = list(maxit = 20))
      tau_new <- resj$par[1]
      rho_new <- resj$par[2]

      Dbeta <- as.matrix(Hm %*% betam_new)
      psim <- Dbeta + vm/nu

      deltam_new <- t(do.call("cbind",lapply(1:nrow(psim),function(xx) mcp(psim[xx,],lam,nu,gam))))

      vm <- vm + nu * (Dbeta - deltam_new)

      diff_norm <- sqrt(sum((Dbeta - deltam_new)^2))
      diff2 <- max(abs(Dbeta- deltam_new))

      deltam <- deltam_new
      beta_cur <- beta_new
      tau_cur <- tau_new
      rho_cur <- rho_new

      if(diff2 < tol)
      {break}
    }

    cluster_est <- get_sp_cluster2(graph_mst, m, len_d, betam_new)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)


    Cm <- (Im - rho_cur* t(adjW)) %*% (Im - rho_cur*adjW)
    Sm <- diag(vardir) + tau_cur * solve(Cm)

    muhat <- z %*% eta_new +  rowSums(x_d * beta_c[cluster_est,-1])
    ratio_est <- tau_cur * solve(Cm) %*% solve(Sm)
    area_est <- muhat + ratio_est %*%(y - muhat)

    out <- list(beta = betam_new, eta = eta_new,
                betac = beta_c,
                tau = tau_cur,
                rho = rho_cur,
                cluster = cluster_est,
                muhat = muhat, area_est = area_est, vardir = vardir,
                deltam = deltam, niters = j)

  }

  out$beta00 <- betam0
  out$loss <- lossSARfun(tau_cur, rho_cur, vardir, resids_new, adjW)
  out$diff_max <- diff2


  return(out)

}

#' @export
ClustSAR_area_mst_ada <- function(dom, y, x, index_d, vardir,
                                graph0, betam0, tau0 = NULL, rho0 = NULL, adjW,
                                nu = 1, gam = 3, lam = 0.5,
                                lam0 = 0.001, maxiter= 500,
                                tol = 1e-4)
{
  m <- length(unique(dom))
  res <- ClustSAR_area_mst(dom, y, x, index_d, vardir,
                         graph0 = graph0, betam0 = betam0,
                         tau0 = tau0, rho0 = rho0, adjW = adjW,
                          nu = nu, gam = gam, lam = lam,
                         lam0 = lam0, maxiter= maxiter,
                         tol = tol)

  graph_updated <- graph0
  E(graph_updated)$weight <- apply(incidentGen(graph0, m) %*% as.matrix(res$betac[res$cluster,-1]), 1, norm, "2")

  #betam0_up <- initialfun2(dom = dom, y = y, x = x, index_d, graph_updated)

  betam0_up <- as.matrix(res$betac[res$cluster,-1])
  res_updated <- ClustSAR_area_mst(dom, y, x, index_d, vardir, graph0 = graph_updated,
                                 betam0 = betam0_up, tau0 = tau0, rho0 = rho0, adjW = adjW,
                                 nu = nu, gam = gam, lam = lam,
                                 lam0 = lam0, maxiter= maxiter,
                                 tol = tol)

  return(list(res = res_updated, graph_updated = graph_updated))
}



#' @export
ClustSAR_area_mst_ada_bic <- function(dom, y, x, index_d, vardir,
                                    graph0, betam0, tau0 = NULL, rho0 = NULL,
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
  SAR_est <- matrix(0, nlam, 2)

  for(j in 1:length(lamvec))
  {
    resj <- ClustSAR_area_mst_ada(dom, y, x, index_d, vardir,
                                graph0 = graph0,
                                betam0 = betam0, tau0 = tau0, rho0 = rho0,
                                adjW = adjW,
                                nu = nu, gam = gam, lam = lamvec[j],
                                lam0 = lam0, maxiter= maxiter,
                                tol = tol)
    resj_updated <- resj$res

    loss_value0[j] <- resj_updated$loss
    K_est[j] <- length(unique(resj_updated$cluster))
    beta_array[,,j] <- resj_updated$beta
    SAR_est[j,] <- c(resj_updated$tau, resj_updated$rho)
  }



  res_df <- data.frame(lambda = lamvec,Khat = K_est,
                       loss0 = loss_value0)

  return(list(res_df = res_df, beta_array = beta_array, SAR_est = SAR_est))
}


## grad of function of tau and rho
# resids y - zeta - xbeta
# only for each dom one observation

#' @export
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


#' @export
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

#
# tempfun <- function(xx)
# {
#   lossSARfun(xx[1],xx[2],vardir, resids, adjW)
# }
#
#
#
# grad(func = tempfun, c(tau,rho)) -gradSARfun(tau, rho, vardir, resids, adjW)
#




# incidentGen <- function(graph, m) {
#   # return the incident matrix of a graph
#   linkMatrix <- igraph::as_data_frame(graph) %>%
#     dplyr::mutate(edge = row_number()) %>% dplyr::select(edge, from, to)
#   linkMatrix$from <- linkMatrix$from %>% as.integer()
#   linkMatrix$to <- linkMatrix$to %>% as.integer()
#   H <- as(matrix(0, dim(linkMatrix)[1], m), "dgCMatrix")
#   H[as.matrix(linkMatrix[, c(1, 2)])] = 1
#   H[as.matrix(linkMatrix[, c(1, 3)])] = -1
#   H
# }

# get_sp_cluster2 <- function(graph_mst, m, len_d, betam_new)
# {
#   tol <- 1e-3
#   E(graph_mst)$weight = apply(incidentGen(graph_mst, m) %*% betam_new, 1, norm, "2")
#   link_matrix <- igraph::as_data_frame(graph_mst) %>%
#     dplyr::mutate(edge = row_number()) %>% filter(weight < tol) %>%
#     dplyr::select(edge, from, to)
#
#   graph_updated <-  graph_from_data_frame(link_matrix[, 2:3], vertices = c(1:m), directed = F)
#
#   cluster_est <- components(graph_updated)$membership
#   return(cluster_est)
# }
#



#### initial values ####
#
# initialfun2 <- function(dom, y, x, index_d, graph0, lam0 = 0.001,
#                         method = "kmeans", k0 = 10)
# {
#   n0 <- length(y)
#   uindex <- unique(dom)
#   m <- length(uindex)
#
#   len_d <- length(index_d)
#   Ip <- diag(1,len_d,len_d)
#
#   Hm <- incidentGen(graph0, m)
#   #Hm <- rbind(matrix(1/sqrt(m),1,m), Hm)
#
#   AtA <- (t(Hm) %*% Hm) %x% Ip
#
#   x_d <- x[,index_d, drop = FALSE]
#   Xm <- matrix(0, n0, m*len_d)
#   for(i in 1:m)
#   {
#     Xm[dom == uindex[i],(len_d*(i-1) + 1) : (len_d*i)] <- x_d[dom == uindex[i],]
#   }
#
#
#   if(len_d == ncol(x))
#   {
#     beta00 <- matrix(solve(t(Xm) %*% Xm + lam0*AtA) %*% t(Xm) %*% y,
#                      ncol = len_d, byrow = TRUE)
#
#     # res00_ada1  <- adaptive_start(beta00, k0 = k0, method = method)
#     # cluster00 <- res00_ada1$cluster00
#     # res00_ada2 <- lm_clust_reg(dom = dom, y = y, x = x, index_d = index_d,
#     #                            cluster = cluster00)
#     # beta00_ada <- matrix(res00_ada2, ncol = len_d, byrow = TRUE)[cluster00,,drop= FALSE]
#   }
#
#   if(len_d < ncol(x))
#   {
#     z <- x[,-index_d, drop = FALSE]
#     ztz <- solve(t(z) %*% z)%*%t(z)
#     Qz <- diag(1, n0, n0) - z%*% ztz
#
#     beta00 <- matrix(solve(t(Xm) %*% Qz%*% Xm + lam0*AtA) %*% t(Xm)%*% Qz %*% y,
#                      ncol = len_d, byrow = TRUE)
#
#     # res00_ada1  <- adaptive_start(beta00, k0 = k0, method = method)
#     # cluster00 <- res00_ada1$cluster00
#     # res00_ada2 <- lm_clust_reg(dom = dom, y = y, x = x, index_d = index_d,
#     #                            cluster = cluster00)
#     # beta00_ada <- matrix(res00_ada2[-(1:ncol(z))], ncol = len_d, byrow = TRUE)[cluster00,,drop= FALSE]
#   }
#
#
#   return(beta00)
# }
#


