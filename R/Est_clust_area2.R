#'  wrap up function pairwise
#'  @import sae
#' @export
Est_clust_area2 <- function(dom, y, x, index_d, vardir, lamvec,
                           nu = 1, gam = 3,
                           lam0 = 0.001, maxiter= 500, tol = 1e-4,
                           seed = 1234)
{
  n0 <- length(y)
  len_d <- length(index_d)
  nz <- ncol(x) - len_d
  res_bic <- Clust_area_pw_bic2(indexy=dom, y, x, index_d, vardir,
                               lamvec, nu = nu, gam = gam,
                               lam0 = lam0, maxiter= maxiter, tol = tol,
                               seed = seed)

  npp <- length(unique(dom))

  bic_df <- res_bic %>%
    mutate(bic = loss0 + log(n0)*(Khat*len_d + nz + 1),
           mbic = loss0 + log(npp)*log(n0)*(Khat*len_d + nz +1))

  ### fitted
  res_f <- Clust_area_pw2(indexy = dom, y, x, index_d, vardir,
                         nu = nu, gam = gam, lam = lamvec[which.min(bic_df$bic)],
                         lam0 = lam0, maxiter= maxiter, tol = tol,
                         seed = seed)
  area_est <- res_f$area_est

  tryCatch({
    fit_f <- fit_clust_area(dom = dom, y = y, x = x, index_d = index_d,
                            vardir = vardir, cluster = res_f$cluster)
    area_est <- fit_f$area_est_ada
  }, error = function(e) {
    area_est <- res_f$area_est
  })

  ##
  res_fm <- Clust_area_pw2(indexy = dom, y, x, index_d, vardir,
                          nu = nu, gam = gam, lam = lamvec[which.min(bic_df$mbic)],
                          lam0 = lam0, maxiter= maxiter, tol = tol,
                          seed = seed)
  area_estm <- res_fm$area_est

  tryCatch({
    fit_fm <- fit_clust_area(dom = dom, y = y, x = x, index_d = index_d,
                             vardir = vardir, cluster = res_fm$cluster)
    area_estm <- fit_fm$area_est_ada
  }, error = function(e) {
    area_estm <- res_fm$area_est
  })


  out <- list(resf = res_f,
              #fit = fit_f,
              bic_df = bic_df,
              cluster_est = res_f$cluster,
              area_est = area_est,
              resfm = res_fm,
              #fitm = fit_fm,
              cluster_estm = res_fm$cluster,
              area_estm = area_estm)

  return(out)

}



