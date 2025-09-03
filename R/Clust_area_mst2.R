#' clustered area level model, new mst algorithm based on admm

#' @import sae
#' @importFrom dplyr filter select mutate union groups row_number
#' @importFrom igraph as_data_frame mst E graph_from_data_frame components
#' @importFrom igraph "E<-"
#' @importFrom magrittr %>%
#' @importClassesFrom Matrix dgCMatrix
#' @import Matrix
#' @return
#' @export
Clust_area_mst2 <- function(dom, y, x, index_d, vardir, graph0, betam0,
                           rho = 1, gam = 3, lam = 0.5,
                           lam0 = 0.001, maxiter= 500, tol = 1e-4)
{

  n0 <- length(y)
  uindex <- unique(dom)
  m <- length(uindex)
  ns <- as.numeric(table(dom)) # number of replicates
  nbar <- mean(ns)

  len_d <- length(index_d)
  Ip <- diag(1,len_d,len_d)

  ### mst and H matrix, the graph is based on nr###
  E(graph0)$weight <- apply(incidentGen(graph0, m) %*% as.matrix(betam0), 1, norm, "2")
  graph_mst <- igraph::mst(graph0)
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
    deltam <- Am %*% beta_cur
    vm <-  matrix(0, (m-1) , len_d)

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
      wm <- 1/(m*nbar*(vardir + sig2_cur))
      Xmt <- t(wm * Xm)

      rho <- mean(1/(vardir + sig2_cur))

      beta_new <- solve(Xmt %*% Xm + rho * AtA) %*% (Xmt %*% y + rho*c(as.matrix(t(deltam -  vm/rho) %*% Hm)))
      betam_new <- matrix(beta_new, m, len_d, byrow = TRUE)

      # update sigma_v^2
      resids <- y -Xm %*% beta_new

      lossj <- function(xx)
      {
        lossfun2(xx, vardir, resids)
      }

      sig2_new <- optimize(lossj, interval = c(0,mean(resids^2)))$minimum

      Dbeta <- Hm %*% betam_new
      psim <- Dbeta + vm/rho
      deltam_new <- t(do.call("cbind",lapply(1:nrow(psim),function(xx) mcp(psim[xx,],lam,rho,gam))))

      vm <- vm + rho * (Dbeta - deltam_new)

      diff_norm <- sqrt(sum((Dbeta- deltam_new)^2))
      diff2 <- max(abs(Dbeta- deltam_new))

      deltam <- deltam_new
      sig2_cur <- sig2_new

      if(diff2 < tol)
      {break}


    }


    cluster_est <- get_sp_cluster2(graph_mst, m, len_d, betam_new)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)

    muhat <- rowSums(x * betam_new[dom,])
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

    eta_cur <- ztz %*% (y - Xm%*%beta_cur)

    ##### initial of others
    deltam <- Am %*% beta_cur
    vm <-  matrix(0, (m-1) , len_d)

    resids0 <- y - z %*% eta_cur - Xm %*% beta_cur

    lossj <- function(xx)
    {
      lossfun2(xx, vardir, resids0)
    }

    sig2_cur <- optimize(lossj, interval = c(0,mean(resids0^2)))$minimum


    ### iterations
    for(j in 1:maxiter)
    {
      # update beta and eta
      wm <- 1/m*(1/(vardir + sig2_cur)*(1/nbar))
      Wm <- diag(wm)
      ztz <-  solve(t(wm*z) %*% z)%*%t(wm *z)
      Qz <- Wm - Wm %*% z%*% ztz

      rho <- mean(1/(vardir + sig2_cur))

      Xinv <- solve(t(Xm) %*%Qz%*% Xm + rho * AtA)
      beta_new <- Xinv %*% (t(Xm) %*% Qz %*% y + rho*c(as.matrix(t(deltam -  vm/rho) %*% Hm)))

      betam_new <- matrix(beta_new, m, len_d, byrow = TRUE)
      eta_new <- ztz %*% (y - Xm%*%beta_new)


      # update sigma_v^2
      resids_new <- y - z %*% eta_new - Xm %*% beta_new

      lossj <- function(xx)
      {
        lossfun2(xx, vardir, resids_new)
      }

      sig2_new <- optimize(lossj, interval = c(0,mean(resids_new^2)))$minimum

      Dbeta <- as.matrix(Hm %*% betam_new)
      psim <- Dbeta + vm/rho

      deltam_new <- t(do.call("cbind",lapply(1:nrow(psim),function(xx) mcp(psim[xx,],lam,rho,gam))))

      vm <- vm + rho * (Dbeta - deltam_new)

      diff_norm <- sqrt(sum((Dbeta - deltam_new)^2))
      diff2 <- max(abs(Dbeta- deltam_new))

      deltam <- deltam_new
      sig2_cur <- sig2_new
      beta_cur <- beta_new

      if(diff2 < tol)
      {break}
    }

    cluster_est <- get_sp_cluster2(graph_mst, m, len_d, betam_new)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)

    muhat <- z %*% eta_new +  rowSums(x_d * betam_new[dom,])
    ratio_est <- sig2_cur/(vardir + sig2_cur)
    area_est <- ratio_est * y + (1-ratio_est)*muhat

    out <- list(beta = betam_new, eta = eta_new,
                betac = beta_c, sig2 = sig2_cur,
                cluster = cluster_est,
                muhat = muhat, area_est = area_est, vardir = vardir,
                deltam = deltam, niters = j)

  }

  out$beta00 <- betam0
  out$loss <- lossfun2(sig2_new, vardir, resids_new)
  out$diff_max <- diff2


  return(out)

}

#' @export
Clust_area_mst_ada2 <- function(dom, y, x, index_d, vardir,
                                graph0, betam0,
                               rho = 1, gam = 3, lam = 0.5,
                               lam0 = 0.001, maxiter= 500,
                               tol = 1e-4)
{
  m <- length(unique(dom))
  res <- Clust_area_mst2(dom, y, x, index_d, vardir,
                         graph0 = graph0, betam0 = betam0,
                        rho = rho, gam = gam, lam = lam,
                        lam0 = lam0, maxiter= maxiter,
                        tol = tol)

  graph_updated <- graph0
  E(graph_updated)$weight <- apply(incidentGen(graph0, m) %*% as.matrix(res$betac[res$cluster,-1]), 1, norm, "2")

  #betam0_up <- initialfun2(dom = dom, y = y, x = x, index_d, graph_updated)

  betam0_up <- as.matrix(res$betac[res$cluster,-1])
  res_updated <- Clust_area_mst2(dom, y, x, index_d, vardir, graph0 = graph_updated,
                                betam0 = betam0_up,
                                rho = rho, gam = gam, lam = lam,
                                lam0 = lam0, maxiter= maxiter,
                                tol = tol)

  return(list(res = res_updated, graph_updated = graph_updated))
}



#' @export
Clust_area_mst_ada_bic2 <- function(dom, y, x, index_d, vardir,
                                   graph0, betam0, lamvec,
                                   rho = 1, gam = 3, lam0 = 0.001,
                                   maxiter= 500, tol = 1e-4)
{
  m <- length(unique(dom))
  nlam <- length(lamvec)
  len_d <- length(index_d)


  K_est <- rep(0, nlam)
  loss_value0 <- loss_value <- rep(0, nlam)
  beta_array <- array(0, dim = c(m, len_d, nlam))


  for(j in 1:length(lamvec))
  {
    resj <- Clust_area_mst_ada2(dom, y, x, index_d, vardir,
                               graph0 = graph0, betam0 = betam0,
                               rho = rho, gam = gam, lam = lamvec[j],
                               lam0 = lam0, maxiter= maxiter,
                               tol = tol)
    resj_updated <- resj$res

    loss_value0[j] <- resj_updated$loss
    K_est[j] <- length(unique(resj_updated$cluster))
    beta_array[,,j] <- resj_updated$beta
  }



  res_df <- data.frame(lambda = lamvec,Khat = K_est,
                       loss = loss_value, loss0 = loss_value0)

  return(list(res_df = res_df, beta_array = beta_array))
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



incidentGen <- function(graph, m) {
  # return the incident matrix of a graph
  linkMatrix <- igraph::as_data_frame(graph) %>%
    dplyr::mutate(edge = row_number()) %>% dplyr::select(edge, from, to)
  linkMatrix$from <- linkMatrix$from %>% as.integer()
  linkMatrix$to <- linkMatrix$to %>% as.integer()
  H <- as(matrix(0, dim(linkMatrix)[1], m), "dgCMatrix")
  H[as.matrix(linkMatrix[, c(1, 2)])] = 1
  H[as.matrix(linkMatrix[, c(1, 3)])] = -1
  H
}

#' @export
get_sp_cluster2 <- function(graph_mst, m, len_d, betam_new)
{
  tol <- 1e-3
  E(graph_mst)$weight = apply(incidentGen(graph_mst, m) %*% betam_new, 1, norm, "2")
  link_matrix <- igraph::as_data_frame(graph_mst) %>%
    dplyr::mutate(edge = row_number()) %>% filter(weight < tol) %>%
    dplyr::select(edge, from, to)

  graph_updated <-  graph_from_data_frame(link_matrix[, 2:3], vertices = c(1:m), directed = F)

  cluster_est <- components(graph_updated)$membership
  return(cluster_est)
}




#' @export

initialfun2 <- function(dom, y, x, index_d, graph0, lam0 = 0.001,
                       method = "kmeans", k0 = 10)
{
  n0 <- length(y)
  uindex <- unique(dom)
  m <- length(uindex)

  len_d <- length(index_d)
  Ip <- diag(1,len_d,len_d)

  Hm <- incidentGen(graph0, m)
  #Hm <- rbind(matrix(1/sqrt(m),1,m), Hm)

  AtA <- (Matrix::t(Hm) %*% Hm) %x% Ip

  x_d <- x[,index_d, drop = FALSE]
  Xm <- matrix(0, n0, m*len_d)
  for(i in 1:m)
  {
    Xm[dom == uindex[i],(len_d*(i-1) + 1) : (len_d*i)] <- x_d[dom == uindex[i],]
  }


  if(len_d == ncol(x))
  {
    beta00 <- matrix(solve(t(Xm) %*% Xm + lam0*AtA) %*% t(Xm) %*% y,
                     ncol = len_d, byrow = TRUE)

    # res00_ada1  <- adaptive_start(beta00, k0 = k0, method = method)
    # cluster00 <- res00_ada1$cluster00
    # res00_ada2 <- lm_clust_reg(dom = dom, y = y, x = x, index_d = index_d,
    #                            cluster = cluster00)
    # beta00_ada <- matrix(res00_ada2, ncol = len_d, byrow = TRUE)[cluster00,,drop= FALSE]
  }

  if(len_d < ncol(x))
  {
    z <- x[,-index_d, drop = FALSE]
    ztz <- solve(t(z) %*% z)%*%t(z)
    Qz <- diag(1, n0, n0) - z%*% ztz

    beta00 <- matrix(solve(t(Xm) %*% Qz%*% Xm + lam0*AtA) %*% t(Xm)%*% Qz %*% y,
                     ncol = len_d, byrow = TRUE)

    # res00_ada1  <- adaptive_start(beta00, k0 = k0, method = method)
    # cluster00 <- res00_ada1$cluster00
    # res00_ada2 <- lm_clust_reg(dom = dom, y = y, x = x, index_d = index_d,
    #                            cluster = cluster00)
    # beta00_ada <- matrix(res00_ada2[-(1:ncol(z))], ncol = len_d, byrow = TRUE)[cluster00,,drop= FALSE]
  }


  return(beta00)
}



