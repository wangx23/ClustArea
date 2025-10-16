#' based on given data, give the area level model estimator with SAR
#' @import sae
#' @importFrom dplyr filter select mutate union groups row_number
#' @importFrom igraph as_data_frame mst E graph_from_data_frame components
#' @importFrom igraph "E<-"
#' @importFrom magrittr %>%
#' @importClassesFrom Matrix dgCMatrix
#' @import spdep
#' @return
#' @export
Est_clustSAR_area_mst <- function(dom, y, x, index_d, vardir,
                               graph0, betam0, tau0 = NULL, rho0 = NULL, adjW,
                               lamvec, nu = 1, gam = 3,
                               lam0 = 0.001, maxiter= 500, tol = 1e-4)
{
  n0 <- length(y)
  len_d <- length(index_d)
  nz <- ncol(x) - len_d
  m <- length(unique(dom))

  res_bic <- ClustSAR_area_mst_ada_bic(dom=dom, y, x, index_d, vardir,
                                     graph0= graph0, betam0 = betam0,
                                     tau0 = tau0, rho0 = rho0, adjW = adjW,
                                     lamvec = lamvec,
                                     nu = nu, gam = gam,
                                     lam0 = lam0, maxiter= maxiter,
                                     tol = tol)



  bic_df <- res_bic$res_df %>%
    mutate(bic = loss0 + log(n0)*(Khat*len_d + nz + 1),
           mbic = loss0 + log(m)*log(n0)*(Khat*len_d + nz +1))


  res_fm <- ClustSAR_area_mst_ada(dom = dom, y, x, index_d, vardir,
                                  graph0 = graph0,betam0 = betam0,
                                  tau0 = tau0, rho0 = rho0, adjW = adjW,
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


