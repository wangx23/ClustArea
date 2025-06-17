# some functions used in the algorithm
#### MCP and SCAD ####
sfun <- function(x, th)
{
  xn <- sqrt(sum(x^2))
  thval <- 1 - th/xn
  
  if(xn ==0)
  {
    value <- xn
  }else{
    value <- thval*((thval) >0)*x
  }
}
mcp <- function(x,lam,nu,gam)
{
  temp <- gam*lam
  xn <- sqrt(sum(x^2))
  if(xn <= temp)
  {
    z <- sfun(x,lam/nu) / (1 - 1/(gam*nu))
  }else{
    z <- x
  }
  z <- as.matrix(z)
  return(z)
}

### penalty function 
mcp_pen <- function(x, lam, gam)
{
  xn <- sqrt(sum(x^2))
  if (xn <=  gam * lam){
    lam * xn - xn ^ 2 / (2 * gam) 
  }else{
    gam * lam^2 / 2
  }
}

scad <- function(x,lam,nu,gam)
{
  temp1 <- lam/nu
  temp2 <- gam * lam
  xn <- sqrt(sum(x^2))
  
  if(xn <= lam + temp1)
  {
    z <- sfun(x, temp1)
  }else if(xn <= temp2 & xn >= lam + temp1)
  {
    z <- sfun(x, temp2/((gam-1)*nu))/(1- 1/((gam - 1 )*nu))
  }else{
    z <- x
  }
  
  z <- as.matrix(z)
  
  return(z)
  
}

adaptive_start <- function(beta00, method = "kmeans", k0 = 10, 
                           seed = 1234)
{

  nr <- nrow(beta00)
  np <- ncol(beta00)
  if(np > 1)
  {
    if(method =="kmeans")
    {
      set.seed(seed)
      
      kk <- max(k0, round(sqrt(nr)))
      reskm <- kmeans(beta00, centers = kk, nstart = 50)
      beta00 <- as.matrix(reskm$centers[reskm$cluster,])
      cluster00 <- reskm$cluster
    }
    
    if(method =="scc")
    {
      library(scclust) 
      ### initial based on scc ###
      dist_locs <- distances::distances(beta00)
      size_constraint <- min(round(nr/k0),round(sqrt(nr)))
      scc <- sc_clustering(dist_locs, size_constraint)
      #scc <- sc_clustering(dist_locs, size_constraint, seed_method = "exclusion_order")
      refined_clustering <- hierarchical_clustering(dist_locs,
                                                    size_constraint = size_constraint,
                                                    existing_clustering = scc)
      
      index_scc <- as.vector(refined_clustering) + 1
      
      beta00 <- Reduce("rbind",lapply(1:max(index_scc), function(xx){colMeans(beta00[index_scc==xx,,drop=FALSE])}))
      beta00 <- beta00[index_scc,,drop= FALSE]
      cluster00 <- index_scc
    }
  }
  
  if(np ==1)
  {
    ### just split equally 
    kk <- max(k0, round(sqrt(nr)))
    quantile_breaks <- quantile(beta00, probs = seq(0, 1, length.out = kk+1), na.rm = TRUE)
    x_factor <- cut(beta00, breaks = quantile_breaks, include.lowest = TRUE, labels = 1:(kk))
    
    cluster00 <- as.numeric(x_factor)
    beta00 <- Reduce("rbind",lapply(1:max(cluster00), function(xx){colMeans(beta00[cluster00==xx,,drop=FALSE])}))
    beta00 <- beta00[cluster00,,drop= FALSE]
    
  }

  return(list(beta00 = beta00, cluster00 = cluster00))
}


### linear regression model without vardir ###
lm_clust_reg <- function(dom, y, x, index_d, cluster)
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
  
  # df_f <- as.data.frame(cbind(y, Ux))
  
  coefu <- as.numeric(solve(t(Ux) %*% Ux + 0.0001*diag(ncol(Ux))) %*% t(Ux) %*% y)
  # forlm <-as.formula(paste("y~0+",paste(colnames(Ux), collapse = "+")))
  # res_lm <- lm(forlm, data = df_f)
  # 
  # return(coef(res_lm))
  return(coefu)
}




getgroup = function(deltam, n, tol = 1e-2)
{
  p = nrow(deltam)
  b2value =sqrt(colMeans(deltam^2))
  b2value[b2value <= tol] = 0
  
  d2 = matrix(0, n, n)
  for(j in 1:(n-1))
  {
    indexj1 = (2*n -j)*(j-1)/2 + 1
    indexj2 = indexj1 + n - j - 1
    d2[(n - indexj2 + indexj1):n,j] = b2value[indexj1:indexj2]
  }
  d2 = t(d2) + d2
  
  
  ngj = 1:n
  groupest = rep(0,n)
  j = 0
  
  while (length(ngj) >0)
  {
    j = j + 1
    gj = (1:n)[d2[ngj[1],] ==0]
    indexj = ngj %in% gj
    gj = ngj[indexj]
    ngj = ngj[!indexj]
    groupest[gj] = j * rep(1,length(gj))
  }
  
  
  return(groupest)
  
}



#### quick inverse ###
inverseR <- function(indexy, x,nu = 1)
{
  uindex <- unique(indexy)
  n <- length(uindex)
  p <- ncol(x)
  Ip <- diag(1,p,p)
  
  DB <- matrix(0,p,p)
  
  A0B <- matrix(0, n*p, p)
  matinv <- matrix(0, n*p, n*p)
  
  for(i in 1:n)
  {
    xmati <- x[indexy==uindex[i],,drop=FALSE]
    matj <- solve(t(xmati)%*%xmati + nu*n*Ip)
    DB <- DB  + matj
    indexi <- ((i-1)*p+1):(i*p)
    A0B[indexi,] <- matj
    matinv[indexi,indexi] <- matj
  }
  
  temp1 <- solve(1/nu*Ip - DB)
  
  A1 <- A0B %*% temp1 %*% t(A0B)
  
  
  matinv <- matinv + A1 
  return(matinv)
}





inverseR2 <- function(indexy,z, x,nu = 1)
{
  uindex <- unique(indexy)
  n <- length(uindex)
  p <- ncol(x)
  q <- ncol(z)
  Ip <- diag(1,p,p)
  
  # D^A(-1)B ###
  DB <- matrix(0,p,p)
  zA0z <- matrix(0, q, q)
  A0z <- matrix(0,p,q)
  
  A0B <- matrix(0, n*p, p)
  A0Z <- matrix(0, n*p, q)
  matinv <- matrix(0, n*p, n*p)
  
  for(i in 1:n)
  {
    xmati <- x[indexy==uindex[i],,drop=FALSE]
    zmati <- z[indexy==uindex[i],,drop=FALSE]
    matj <- solve(t(xmati)%*%(xmati) + nu*n*Ip)
    DB <- DB  + matj
    zA0z <- zA0z + t(zmati)%*%(zmati) -  t(zmati)%*%(xmati) %*% matj%*% t(xmati)%*%(zmati)
    A0z <- A0z + matj%*% t(xmati)%*%(zmati)
    indexi <- ((i-1)*p+1):(i*p)
    A0B[indexi,] <- matj
    A0Z[indexi,] <- matj%*%t(xmati)%*%(zmati)
    matinv[indexi,indexi] <- matj
    
  }
  
  temp1 <- solve(1/nu*Ip - DB - A0z %*% solve(zA0z) %*%t(A0z))
  
  A1 <- A0B %*% temp1 %*% t(A0B)
  A2 <- A0Z %*% solve(zA0z)%*% t(A0z)
  A3 <- A2 %*%temp1
  A22 <- A3 %*% t(A2)
  A12 <- A3 %*% t(A0B)
  
  matinv <- matinv+ A0Z %*% solve(zA0z) %*% t(A0Z) + A1 + A22 + A12 + t(A12)
  return(matinv)
}



#### fit the model with given cluster structure and output is the estimate####
# fit_clust_reg <- function(dom, y, x, index_d, vardir, cluster)
# {
#   n0 <- length(y)
#   uindex <- unique(dom)
#   nr <- length(uindex)
#   ns <- as.numeric(table(dom)) # number of replicates
#   len_d <- length(index_d)
#   
#   x_d <- x[,index_d, drop = FALSE]
#   
#   ng <- length(unique(cluster))
#   
#   
#   Ux <- matrix(0, n0, ng*len_d)
#   for(j in 1:ng)
#   {
#     uindexj <- uindex[cluster==j]
#     indexj <- dom %in% uindexj
#     Ux[indexj,(len_d*(j-1) + 1) : (len_d*j)] <- x_d[indexj,]
#   }
#   
#   names_d <- paste(rep(paste("x",1:len_d,sep=""),ng), rep(paste("c",1:ng,sep=""), each = len_d),sep="_")
#   if(len_d == ncol(x))
#   {
#     #Ux <- Xm%*%W
#     Ux <- Ux
#     colnames(Ux) <- names_d
#   }else{
#     #Ux <- cbind(x[,-index_d, drop = FALSE], Xm%*%W)
#     Ux <- cbind(x[,-index_d, drop = FALSE], Ux)
#     colnames(Ux) <- c(paste("z",1:(ncol(x)-len_d),sep=""), names_d)
#   }
#   
#   df_f <- as.data.frame(cbind(y, vardir, Ux))
#   
#   forFH <-as.formula(paste("y~0+",paste(colnames(Ux), collapse = "+")))
#   res_fh <- eblupFH(forFH, vardir = vardir, data = df_f, method="ML")
#   
#   return(res_fh$fit$estcoef)
# }