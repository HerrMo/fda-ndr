## computing optimal embeddings

# Script to compute the embeddings based on the optimal values obtained via tuning.
#
t1 <- Sys.time()
source("R/setup.R")


meas <- c("q_local", "auc_rnx")
metric <- c("ps_l2", "ps_geo", "fs_dist", "fs_geo") # fs_dist = fs_l2



                  ######################
                  ### SIMULATED DATA ###
                  ######################

l_embs <- l_vals <- vector("list", 11 * 5 * 4) # 11 sets, 5 meths, 4 space-metric-combos
full_names <- vector("character", 11 * 5 * 4)


# compute embeddings
count <- 1
set.seed(4)
for (met in meths) {
  for (dat in settings) {
    # get dist object of setting
    dsts <- dist_list1[[dat]]$dists

    # optimal entries by space for method and data
    opt_vals <- get_opt_vals(met, dat)

    # compute embedding for every space and metric combi
    for (i in 1:4) {
      # save parameter values etc.
      l_vals[[count]] <- opt_vals[i, ]

      # prepare param vals for do.call
      vals <-
        opt_vals %>%
        ungroup() %>%
        select(-data, -method, -space, -value)
      if (met == "mds") {
        # mds dimension parameter must be renamed correctly to k
        spc_vals <- c(list(dist_obj = dsts, method = met),
                      list(k = unname(vals[i, ] %>% pull(1))))
      } else {
        spc_vals <- c(list(dist_obj = dsts, method = met), as.list(vals[i, ]))
      }

      # embedding ---
      if (met != "dmap") {
        l_embs[[count]] <- do.call(embedding_algo, spc_vals)
      } else {
        spc_vals$method <- "diffmap"
        l_embs[[count]] <- do.call(embedding_algo, spc_vals)
      }

      # compute meta data as string and store for plotting as title ---
      ttl <- paste0(dat, ":", " ", met, ".", l_vals[[count]]$space, "-", l_vals[[count]]$value)
      full_names[count] <- ttl

      count <- count + 1
    }
  }
}

names(l_embs) <- full_names
l_embs_lin <- l_embs[str_detect(names(l_embs), "a1_l.l2|p1_l.l2|c1_l.l2|a2_l.l2|p2_l.l2|i2_l.l2")]
l_embs_nl <- l_embs[str_detect(names(l_embs), "a1_l.l2|p1_l.l2|c1_l.l2|a2_l.l2|p2_l.l2|i2_l.l2", negate = TRUE)]





                      #################
                      ### REAL DATA ###
                      #################

#####################                   #####################
#####################   Spectro data    #####################
#####################                   #####################

meths2 <- c("umap", "tsne", "diffmap", "isomap",  "mds")
# mea <- "q_local" # "auc_rnx" # first q_local then auc_rnx

l_embs_spectro <- list()
for (mea in meas) {
  set.seed(1)
  for (mtr in metric[3:4]) {
    if (mea == "auc_rnx") {
      ttl1 <- if (mtr == "fs_dist") ": $AUC^{euc}_{Rnx} -" else ": $AUC^{geo}_{Rnx} -"
    } else {
      ttl1 <- if (mtr == "fs_dist") ": $Q_{local}^{euc} -" else ": $Q_{local}^{geo} -"
    }
    ttl1 <- sapply(meths2, function(x) paste0(x, ttl1))
    opt_vals <- sapply(meths2, get_opt_val, mea = mea, space = mtr, data = res_sp)


    for (i in 1:5) {
      ttl1[i] <- paste0(ttl1[i], round(opt_vals[i], 2))   # TeX(paste0(ttl1[i], round(opt_vals[i], 2), "$"))
    }

    tt <- sapply(meths2, get_opt_entry, mea = mea, space = mtr, data = res_sp)
    ttt <- lapply(tt, function(x) do.call(embedding_algo,
                                          c(list(dist_obj = real_l2), x[,-c(1:3)])))
    names(ttt) <- ttl1
    l_embs_spectro <- c(l_embs_spectro, ttt)
    print(tt)
  }
}


#####################                   #####################
#####################  Earthquake data  #####################
#####################                   #####################

l_embs_earthquake <- list()
for (mea in meas) {
  set.seed(1)
  for (mtr in metric[3:4]) {
    if (mea == "auc_rnx") {
      ttl1 <- if (mtr == "fs_dist") ": $AUC^{euc}_{Rnx} -" else ": $AUC^{geo}_{Rnx} -"
    } else {
      ttl1 <- if (mtr == "fs_dist") ": $Q_{local}^{euc} -" else ": $Q_{local}^{geo} -"
    }
    ttl1 <- sapply(meths2, function(x) paste0(x, ttl1))
    opt_vals <- sapply(meths2, get_opt_val, mea = mea, space = mtr, data = res_eq)


    for (i in 1:5) {
      ttl1[i] <- paste0(ttl1[i], round(opt_vals[i], 2))   # TeX(paste0(ttl1[i], round(opt_vals[i], 2), "$"))
    }

    tt <- sapply(meths2, get_opt_entry, mea = mea, space = mtr, data = res_eq)
    ttt <- lapply(tt, function(x) do.call(embedding_algo,
                                          c(list(dist_obj = eq_l2), x[,-c(1:3)])))
    names(ttt) <- ttl1
    l_embs_earthquake <- c(l_embs_earthquake, ttt)
    print(tt)
  }
}


#####################                   #####################
#####################     COIL data     #####################
#####################                   #####################


l_embs_coil <- list()
set.seed(2)
for (mea in meas) {
  if (mea == "auc_rnx") set.seed(1) # effects become more obvious
  for (mtr in metric[3:4]) {
    if (mea == "auc_rnx") {
      meta_dat <- if (mtr == "fs_dist") ": $AUC^{euc}_{Rnx} -" else ": $AUC^{geo}_{Rnx} -"
    } else {
      meta_dat <- if (mtr == "fs_dist") ": $Q_{local}^{euc} -" else ": $Q_{local}^{geo} -"
    }
    meta_dat <- sapply(meths2, function(x) paste0(x, meta_dat))

    opt_perf_vals <- sapply(meths2, get_opt_val, mea = mea, space = mtr, data = res_coil)
    for (i in 1:5) {
      meta_dat[i] <- paste0(meta_dat[i], round(opt_perf_vals[i], 2))
    }

    opt_params <- sapply(meths2, get_opt_entry, mea = mea, space = mtr, data = res_coil)
    opt_embs <- lapply(opt_params, function(x) do.call(embedding_algo,
                                               c(list(dist_obj = coil_l2), x[,-c(1:3)])))


    names(opt_embs) <- meta_dat
    l_embs_coil <- c(l_embs_coil, opt_embs)
  }
}

Sys.time() - t1
