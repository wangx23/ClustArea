#' fit the model based on a give cluster structure and given vardir ####

#' @export
#' @import sae
fit_clust_area <- function(dom, y, x, index_d, vardir, cluster)
{
  n0 <- length(y)
  uindex <- unique(dom)
  nr <- length(uindex)
  ns <- as.numeric(table(dom)) # number of replicates
  len_d <- length(index_d)

  x_d <- x[,index_d, drop = FALSE]

  ng <- length(unique(cluster))


  Ux <- matrix(0, n0, ng*len_d)
  for(j in 1:ng)
  {
    uindexj <- uindex[cluster==j]
    indexj <- dom %in% uindexj
    Ux[indexj,(len_d*(j-1) + 1) : (len_d*j)] <- x_d[indexj,]
  }

  names_d <- paste(rep(paste("x",1:len_d,sep=""),ng), rep(paste("c",1:ng,sep=""), each = len_d),sep="_")
  if(len_d == ncol(x))
  {
    #Ux <- Xm%*%W
    Ux <- Ux
    colnames(Ux) <- names_d
  }else{
    #Ux <- cbind(x[,-index_d, drop = FALSE], Xm%*%W)
    Ux <- cbind(x[,-index_d, drop = FALSE], Ux)
    colnames(Ux) <- c(paste("z",1:(ncol(x)-len_d),sep=""), names_d)
  }

  df_f <- as.data.frame(cbind(y, vardir, Ux))

  forFH <-as.formula(paste("y~0+",paste(colnames(Ux), collapse = "+")))
  res_mse <- mseFH(forFH, vardir = vardir, data = df_f, method ="REML", MAXITER = 500)
  res_fh <- res_mse$est
  mse <- res_mse$mse


  out <- list(area_est_ada =res_fh$eblup,
              fit =res_fh$fit,
              mse = mse,
              resids = y - Ux %*% res_fh$fit$estcoef[,1])
}
