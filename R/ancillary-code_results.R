
# TODO
# arg dist_list: only if supplied (relevant for sim data only)

plot_grid3d <- function(embs_list, dist_list = NULL, color = NULL, ncol = 5, params = NULL, ylim2 = NULL, subtitle = TRUE,
                        bycol = FALSE, color_pal = rbPal, ttl_cex = 1.5, ...) {
  # define plot settings
  n_plts <- length(embs_list)
  if (bycol) {
    par(mfcol = c(n_plts/ncol, ncol = ncol),
        oma = c(0, 3, 3, 0))
  } else {
    par(mfrow = c(n_plts/ncol, ncol = ncol),
        oma = c(0, 3, 3, 0))
  }

  # specify plot arguments for plots individually
  if (!missing(params)) {
    argmts <- vector("list", length(params))
    names(argmts) <- names(params)
    for (par in seq_along(params)) {
      param <- names(params)[par]

      # default will be deployed to plots which should be left unchanged
      default <- formals(scatterplot3d)[[param]]

      argmts[[par]] <- set_param(n_plts,
                                 default,
                                 params[[param]]$vals,
                                 params[[param]]$pos)
    }
  }
  for (i in seq_along(embs_list)) {

    # extract embedding coordinates
    coords <- extract_points(embs_list[[i]])
    n_obs <- nrow(coords)

    # color code for plot
    col <- 1:n_obs
    if (!is.null(dist_list)) {
      dat <- str_split(names(embs_list[i]), ":")[[1]][1]
      col <- dist_list[[dat]]$data$param[[1]]
    } else if (!is.null(color)) {
      col <- color
    }

    cols <- color_pal(n_obs)[as.numeric(cut(col, breaks = n_obs))]

    # in case ndim > 3
    if (ncol(coords) > 3) coords <- coords[, 1:3]

    print(paste0(i, ": ", names(embs_list[i])))
    if (ncol(coords) == 2) {
      args2d <- list(x = coords,
                     col = cols,
                     pch = 20,
                     # mar = c(1, 0, 0, 1), # uncomment for mds
                     cex.main = 1.5,
                     cex.lab = 0.001,
                     xlab = TeX(""),
                     ylab = TeX(""),
                     xaxt = 'n',
                     yaxt = 'n'
      )
      if (!missing(ylim2)) {
        y_count <- 1
        if (i %in% ylim2$pos) {
          args2d <- c(args2d, list(ylim = c(ylim2$vals1[y_count],
                                            c(ylim2$vals2[y_count]))))
          y_count <- y_count + 1
        }
      }
      do.call(plot, args2d)
    } else {
      temp_args <- c(list(x = coords,
                          color = cols,
                          xlab = TeX(""),
                          ylab = TeX(""),
                          zlab = TeX(""),
                          mar = if (subtitle) c(4, 0, 0, 1) else c(1, 0, 0, 1), # bottom = 1 for mds only and lin sets, else 4
                          tick.marks = FALSE,
                          lab = c(3, 3, 0.5),
                          lab.z = c(3, 0.5)),
                     ...)
      if (!missing(params)) temp_args <- c(temp_args, lapply(argmts, function(x) x[i]))
      do.call(scatterplot3d, temp_args)
    }
    if (subtitle) {
      nice_ttl <- makeNiceTtl(names(embs_list[i]))
      mtext(nice_ttl, side = 1, 2, cex = ttl_cex)
    }
  }

  gridplt <- recordPlot()
  gridplt
}

# TODO: assign attributes to each embedding with metric, value, etc
#       rewrite makeNiceTtl based on attributes
makeNiceTtl <- function(name, setting = TRUE) {
  plt_ttl <- name
  meth <- str_extract(plt_ttl, "umap|isomap|tsne|dmap|mds")
  meth <- switch(meth,
                 "mds" = "MDS",
                 "umap" = "UMAP",
                 "isomap" = "ISOMAP",
                 "tsne" = "t-SNE",
                 "dmap" = "DIFFMAP")
  perf <- round(as.numeric(str_split(plt_ttl, "-")[[1]][2]), 2)
  if (setting) {
    regex_sets <- "a1_l|p1_l|c1_l|a2_l|p2_l|i2_l|a2_sr|a3_hx|a3_sr|a3_sc|a3_tp"
    setting <- str_extract(plt_ttl, regex_sets)
    setting <- str_replace(setting, "_", "-")
  }

  metric <- "dir"
  if (str_detect(plt_ttl, "geo")) metric <- "geo"

  if (str_detect(plt_ttl, "Q_")) {
    ttl <- TeX(paste0("$Q_{local}^{", metric, "} = $", perf))
  } else {
    ttl <- TeX(paste0("$AUC^{", metric, "}_{Rnx} = $", perf))
  }

  ttl
}

set_param <- function(n, default, vals, pos) {
  param <- rep(default, n)
  param[pos] <- vals
  param
}

get_opt_val <- function(meth, mea, space, data) {
  data[method == meth & meas == mea, max(get(space))]
}

get_opt_entry <- function(meth, mea, space, data) {
  temp <- data[method == meth & meas == mea]
  temp <- temp[which.max(temp[[space]])]
  Filter(function(x) !all(is.na(x)), temp)
}

get_opt_vals <- function(meth, dat, res = res_sim2, n = 1) {
  selc <- switch(
    meth,
    "tsne" = c("value", "perplexity", "dims", "theta", "max_iter",
               "eta", "exaggeration", "data", "method", "space"),
    "umap" = c("value", "n_neighbors", "n_components", "min_dist", "n_epochs", "init", "data", "method", "space"),
    "isomap" = c("value", "k", "ndim", "data", "method", "space"),
    "mds" = c("value", "k_m", "data", "method", "space"),
    "dmap" = c("value", "t", "neigen", "eps.val", "data", "method", "space"))

  vals <-
    res %>%
    filter(meas != "q_global") %>%
    filter(meas != "q_local") %>%
    filter(data == dat) %>%
    filter(method == meth) %>%
    group_by(space) %>%
    top_n(n, value) %>%
    select(one_of(selc)) %>%
    arrange(desc(value))

  if (nrow(vals) != (n * 4)) warning("Seems like there are no unique optimal value")
  vals
}

leftpad <- function(padVec, pad.len = max(nchar(padVec))) {
  vapply(padVec,
         function(x) {
           if (nchar(x) > 0) {
             # creating padding for each entry
             padding <- paste(rep(0, pad.len - nchar(x)), collapse = "")
             # pasting padding and corresponding entry
             paste0(padding, x, collapse = "")
           } else {  # if entry ""
             x
           }
         }, "", USE.NAMES = FALSE)
}

plot_funs_temp <- function(data, col = NULL, xscale = NULL) {
  n <- nrow(data)
  grid_len <- ncol(data)
  df_dat <- data.frame(
    args =  if (is.null(xscale)) {rep(1:grid_len, n)} else {rep(xscale, n)},
    vals = c(t(data)),
    id = as.factor(rep(1:n, each = grid_len))
  )

  if (!is.null(col)) df_dat$col <- as.factor(rep(col, each = grid_len))

  ggplot(df_dat) +
    geom_line(aes(x = args,
                  y = vals,
                  group = id,
                  colour = if (is.null(col)) {id} else {col})) +
    theme_classic() +
    theme(legend.position = "None",
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 10),
          plot.margin = margin(t = 1, 0, 0, 0, "cm"))
}
