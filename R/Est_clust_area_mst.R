#' based on given data, give the area level model estimator
#' @import sae
#' @importFrom dplyr filter select mutate union groups row_number
#' @importFrom igraph as_data_frame mst E graph_from_data_frame components
#' @importFrom igraph "E<-"
#' @importFrom magrittr %>%
#' @importClassesFrom Matrix dgCMatrix
#' @return
#' @export
Est_clust_area_mst <- function(dom, y, x, index_d, vardir, graph0,
                               gamma0,
                               lamvec, rho = 1, gam = 3,
                               lam0 = 0.001, maxiter= 200, max_inner = 10, tol = 1e-4)
{
  n0 <- length(y)
  len_d <- length(index_d)
  nz <- ncol(x) - len_d
  m <- length(unique(dom))

  res_bic <- Clust_area_mst_ada_bic(dom=dom, y, x, index_d, vardir, graph0= graph0,
                                    lamvec, gamma0, rho = rho, gam = gam,
                                    lam0 = lam0, maxiter= maxiter,
                                    tol = tol)



  bic_df <- res_bic$res_df %>%
    mutate(bic = loss0 + log(n0)*(Khat*len_d + nz + 1),
           mbic = loss0 + log(m)*log(n0)*(Khat*len_d + nz +1))


  res_fm <- Clust_area_mst_ada(dom = dom, y, x, index_d, vardir, graph0 = graph0,
                               gamma0 = gamma0,
                               rho = rho, gam = gam, lam = lamvec[which.min(bic_df$mbic)],
                               lam0 = lam0, maxiter= maxiter, tol = tol)
  res_fm <- res_fm$res
  area_estm <- res_fm$area_est

  tryCatch({
    fit_fm <- fit_clust_area(dom = dom, y = y, x = x, index_d = index_d,
                             vardir = vardir, cluster = res_fm$cluster)
    area_estm <- fit_fm$area_est_ada
  }, error = function(e) {
    area_estm <- res_fm$area_est
  })




  out <- list(bic_df = bic_df,
              resfm = res_fm,
              cluster_estm = res_fm$cluster,
              area_estm = area_estm)

  return(out)

}


