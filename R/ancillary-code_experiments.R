# experiment ancillary_code
quality_algo <- function(emb, data, ndim, ...) {
  if (inherits(emb, "error")) {
    quals <- data.table(c("auc_rnx2", "q_local", "q_global"), replicate(length(data), rep(NA, length(data))))
    colnames(quals) <- c("meas", names(data))
  } else {  
    points <- extract_points(emb, ndim = ndim) # extract embedding space coordinates in ndim dimensions
    assertMatrix(points)
    rm(emb)
    
    es_dists <- as.matrix(dist(points)) # embedding space distances
    cs_dists <- data # comparing space distances (param space or fs with differing metrics)
     
    quals <- vapply(cs_dists, function(x, y) quality(y, x), y = es_dists, FUN.VALUE = numeric(3))
    quals <- as.data.table(quals, keep.rownames = TRUE)
    colnames(quals) <- c("meas", names(data)) 
  } 
  quals
}  

embedding_algo <- function(dist_obj, method = c("isomap", "umap", "diffmap", "mds", "tsne"), ...) {
  method <- match.arg(method, c("isomap", "umap", "diffmap", "mds", "tsne"))
  
  # change to assert: accept only matrices?
  dist_mat <- if (is.matrix(dist_obj)) {dist_obj} else {as.matrix(dist_obj)}
  
  emb <- tryCatch(
    switch(
      method,
      "isomap" = vegan::isomap(dist_mat, ...),
      "umap" = umap::umap(dist_mat, input = "dist", ...),
      "diffmap" =  diffuse2(dist_mat, ...),
      "mds" = cmdscale(dist_mat, ...),
      "tsne" = Rtsne::Rtsne(dist_mat, is_distance = TRUE, ...)
    ),
    error = function(c) {
      emb <- NA
      class(emb) <- "error"
      emb
    }
  )
  
  if (method == "tsne") class(emb) <- "tsne"
  
  emb
}

nquality <- function(x, ...) {
  UseMethod("quality")
}

quality.default <- function(d1, d2) {
  auc <- auc_rnx(d1, d2, weight = "log10")
  q_local <- local_q(d1, d2)  
  q_global <- global_q(d1, d2)
  
  c(auc_rnx = auc, q_local = q_local, q_global = q_global)
}

quality.embedding <- function(embs, p_space) {
  auc <- auc_rnx(embs, p_space = p_space, ndim = 3, weight = "log10")
  q_local <- local_q(embs, p_space = p_space, ndim = 3)  
  q_global <- global_q(embs, p_space = p_space, ndim = 3)
  
  c(auc_rnx = auc, q_local = q_local, q_global = q_global)
}


R_nx2 <- function(x, ...) {
  UseMethod("R_nx2") 
}

R_nx2.default <- function(d1, d2) {
  assertMatrix(d1)
  assertMatrix(d2)
  
  Q <- coRanking::coranking(d1,
                            d2,
                            input_Xi = "dist")
  
  nQ <- nrow(Q)
  N <- nQ + 1
  
  Qnx <- diag(apply(apply(Q, 2, cumsum), 1, cumsum)) /
    seq_len(nQ) / N
  
  Rnx <- ((N - 1) * Qnx - seq_len(nQ)) /
    (N - 1 - seq_len(nQ))
  Rnx[-nQ]
}

R_nx2.embedding <- function(object, ndim = 2, p_space = FALSE) {
  # chckpkg("coRanking")
  # if (!object@has.org.data) stop("object requires original data")
  ld_dist <- as.matrix(dist(object$points[, 1:ndim]))
  
  if (p_space) {
    Q <- coRanking::coranking(object$p_dist,
                              ld_dist,
                              input_Xi = "dist")
  } else {
    Q <- coRanking::coranking(object$f_dist,
                              ld_dist,
                              input_Xi = "dist")
  }
  
  nQ <- nrow(Q)
  N <- nQ + 1
  
  Qnx <- diag(apply(apply(Q, 2, cumsum), 1, cumsum)) /
    seq_len(nQ) / N
  
  Rnx <- ((N - 1) * Qnx - seq_len(nQ)) /
    (N - 1 - seq_len(nQ))
  Rnx[-nQ]
}

auc_rnx <- function(x, ...) {
  UseMethod("auc_rnx")
}

auc_rnx.default <- function(d1, d2, weight = "inv") {
  rnx <- R_nx2(d1, d2)
  
  weight <- match.arg(weight, c("inv", "ln", "log", "log10"))
  switch(
    weight,
    inv   = auc_ln_k_inv(rnx),
    log   = auc_log_k(rnx),
    ln    = auc_log_k(rnx),
    log10 = auc_log10_k(rnx),
    stop("wrong parameter for weight")
  )
}

auc_rnx.embedding <- function(object, ndim = 2, p_space = FALSE, weight = "inv") {
  rnx <- R_nx2(object, ndim = ndim, p_space = p_space)
  
  weight <- match.arg(weight, c("inv", "ln", "log", "log10"))
  switch(
    weight,
    inv   = auc_ln_k_inv(rnx),
    log   = auc_log_k(rnx),
    ln    = auc_log_k(rnx),
    log10 = auc_log10_k(rnx),
    stop("wrong parameter for weight")
  )
}

auc_ln_k_inv <- function(rnx) {
  Ks <- seq_along(rnx)
  return(sum(rnx / Ks) / sum(1 / Ks))
}

auc_log_k <- function(rnx) {
  Ks <- seq_along(rnx)
  return(sum(rnx * log(Ks)) / sum(log(Ks)))
}

auc_log10_k <- function(rnx) {
  Ks <- seq_along(rnx)
  return(sum(rnx * log10(Ks)) / sum(log10(Ks)))
}

local_q <- function(x, ...) {
  UseMethod("local_q")
}

local_q.default <- function(d1, d2, ...) {
  assertMatrix(d1)
  assertMatrix(d2)
  
  Q <- coRanking::coranking(d1,
                            d2,
                            input_Xi = "dist")
  
  nQ <- nrow(Q)
  N <- nQ + 1
  
  Qnx <- diag(apply(apply(Q, 2, cumsum), 1, cumsum)) / seq_len(nQ) / N
  lcmc <- Qnx - seq_len(nQ) / nQ
  
  Kmax <- which.max(lcmc)
  
  Qlocal <- sum(lcmc[1:Kmax]) / Kmax
  return(as.vector(Qlocal))
}

local_q.embedding <- 
  function(object, ndim = 3, p_space) {
    ld_dist <- as.matrix(dist(object$points[, 1:ndim]))
    
    if (p_space) {
      Q <- coRanking::coranking(object$p_dist,
                                ld_dist,
                                input_Xi = "dist")
    } else {
      Q <- coRanking::coranking(object$f_dist,
                                ld_dist,
                                input_Xi = "dist")
    }
    
    nQ <- nrow(Q)
    N <- nQ + 1
    
    Qnx <- diag(apply(apply(Q, 2, cumsum), 1, cumsum)) / seq_len(nQ) / N
    lcmc <- Qnx - seq_len(nQ) / nQ
    
    Kmax <- which.max(lcmc)
    
    Qlocal <- sum(lcmc[1:Kmax]) / Kmax
    return(as.vector(Qlocal))
  }

global_q <- function(x, ...) {
  UseMethod("global_q")
}

global_q.default <- function(d1, d2) {
  assertMatrix(d1)
  assertMatrix(d2)
  
  Q <- coRanking::coranking(d1,
                            d2,
                            input_Xi = "dist")
  
  nQ <- nrow(Q)
  N <- nQ + 1
  
  Qnx <- diag(apply(apply(Q, 2, cumsum), 1, cumsum)) / seq_len(nQ) / N
  lcmc <- Qnx - seq_len(nQ) / nQ
  
  Kmax <- which.max(lcmc)
  
  Qglobal <- sum(lcmc[(Kmax + 1):nQ]) / (N - Kmax)
  return(as.vector(Qglobal))
}

global_q.embedding <- 
  function(object, ndim = 3, p_space){
    ld_dist <- as.matrix(dist(object$points[, 1:ndim]))
    
    if (p_space) {
      Q <- coRanking::coranking(object$p_dist,
                                ld_dist,
                                input_Xi = "dist")
    } else {
      Q <- coRanking::coranking(object$f_dist,
                                ld_dist,
                                input_Xi = "dist")
    }
    
    nQ <- nrow(Q)
    N <- nQ + 1
    
    Qnx <- diag(apply(apply(Q, 2, cumsum), 1, cumsum)) / seq_len(nQ) / N
    lcmc <- Qnx - seq_len(nQ) / nQ
    
    Kmax <- which.max(lcmc)
    
    Qglobal <- sum(lcmc[(Kmax + 1):nQ]) / (N - Kmax)
    return(as.vector(Qglobal))
  }


# help fun for plot embedding - S3 class to extract embedding coordinates
extract_points <- function(x, ...) {
  UseMethod("extract_points")
}

# S3 method for isomap
extract_points.isomap <- function(embedding, ndim = dim(embedding$points)[2]) {
  embedding$points[, 1:ndim]
}

# S3 method for umap
extract_points.umap <- function(embedding, ndim = dim(embedding$layout)[2]) {
  embedding$layout[, 1:ndim]  
}

# S3 method for diffusionMap
extract_points.diffuse <- function(embedding, ndim = dim(embedding$X)[2]) {
  embedding$X[, 1:ndim]
}

# S3 method for mds
extract_points.matrix <- function(embedding, ndim = dim(embedding)[2]) {
  embedding[, 1:ndim]
}

# S3 method for tsne
extract_points.tsne <- function(embedding, ndim = dim(embedding$Y)[2]) {
  embedding$Y[, 1:ndim]
}

ps_dist <- function(fn_dat, k) {
  as.matrix(
    isomapdist(
      dist(
        do.call(cbind, fn_dat[[2]])
      ),
      k = k
    )
  )
}

diffuse2 <- function(D, eps.val = epsilonCompute(D), neigen = NULL, t = 0, 
                     maxdim = 50, delta = 10^-5, maxiter = 100000) {
  start = proc.time()[3]
  D = as.matrix(D)
  n = dim(D)[1]
  K = exp(-D^2/(eps.val))
  v = sqrt(apply(K, 1, sum))
  A = K/(v %*% t(v))
  ind = which(A > delta, arr.ind = TRUE)
  Asp = sparseMatrix(i = ind[, 1], j = ind[, 2], x = A[ind], 
                     dims = c(n, n))
  f = function(x, A = NULL) {
    as.matrix(A %*% x)
  }
  cat("Performing eigendecomposition\n")
  if (is.null(neigen)) {
    neff = min(maxdim + 1, n)
  }
  else {
    neff = min(neigen + 1, n)
  }
  decomp = arpack(f, extra = Asp, sym = TRUE, options = list(which = "LA", maxiter = maxiter,
                                                             nev = neff, n = n, ncv = max(min(c(n, 4 * neff)))))
  psi = decomp$vectors/(decomp$vectors[, 1] %*% matrix(1, 1, 
                                                       neff))
  phi = decomp$vectors * (decomp$vectors[, 1] %*% matrix(1, 
                                                         1, neff))
  eigenvals = decomp$values
  cat("Computing Diffusion Coordinates\n")
  if (t <= 0) {
    lambda = eigenvals[-1]/(1 - eigenvals[-1])
    lambda = rep(1, n) %*% t(lambda)
    if (is.null(neigen)) {
      lam = lambda[1, ]/lambda[1, 1]
      neigen = min(which(lam < 0.05))
      neigen = min(neigen, maxdim)
      eigenvals = eigenvals[1:(neigen + 1)]
      cat("Used default value:", neigen, "dimensions\n")
    }
    X = psi[, 2:(neigen + 1)] * lambda[, 1:neigen]
  }
  else {
    lambda = eigenvals[-1]^t
    lambda = rep(1, n) %*% t(lambda)
    if (is.null(neigen)) {
      lam = lambda[1, ]/lambda[1, 1]
      neigen = min(which(lam < 0.05))
      neigen = min(neigen, maxdim)
      eigenvals = eigenvals[1:(neigen + 1)]
      cat("Used default value:", neigen, "dimensions\n")
    }
    X = psi[, 2:(neigen + 1)] * lambda[, 1:neigen]
  }
  cat("Elapsed time:", signif(proc.time()[3] - start, digits = 4), 
      "seconds\n")
  y = list(X = X, phi0 = phi[, 1], eigenvals = eigenvals[-1], 
           eigenmult = lambda[1, 1:neigen], psi = psi, phi = phi, 
           neigen = neigen, epsilon = eps.val)
  class(y) = "diffuse"
  return(y)
}