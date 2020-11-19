source("R/setup.R")

# Fig 1
# fda ndr framework

cols <- rbPal(1000)[as.numeric(cut(dist_list1$a3_sc.l2$data$param[[1]], breaks = 1000))]
scatterplot3d(matrix(unlist(dist_list1$a3_sc.l2$data$param, recursive = TRUE),
                     ncol = 3),
              color = cols,
              pch = 20,
              angle = 130,
              box = FALSE,
              axis = FALSE)

plot(l_embs_nl$`a3_sc.l2: isomap.fs_geo-0.911131345906121`$points[, 1],
     l_embs_nl$`a3_sc.l2: isomap.fs_geo-0.911131345906121`$points[, 2],
     col = cols,
     pch = 20,
     xaxt = "n",
     yaxt = "n",
     ann = FALSE)

lrg_smpl <- sample(1:1000, 50)
plot_funs_temp(dist_list1$a3_sc.l2$data$dat$funs[lrg_smpl,],
               col = cols[lrg_smpl]) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_color_manual(values = cols[lrg_smpl])


# Fig 2
# example observations simulated data
thm <- theme(axis.title = element_text(size = rel(1.5)),
             axis.text = element_text(size = rel(1.1)),
             axis.text.x = element_text(angle = -45))

ex_labs <- list(xlab("t"), ylab("x(t)"))
emb_labs <- list(xlab(TeX("y_1")), ylab(TeX("y_2")))

smpl <- sample(1:1000, 10)

xscl <- seq(0, 1, length = 200)

ex_data_p2 <-
  plot_funs_temp(dist_list1$p2_l.l2$data$dat$funs[smpl, ],
             col = dist_list1$p2_l.l2$data$dat$param[[1]][smpl],
             xscale = xscl) +
  ex_labs + thm
ex_data_ap <-
  plot_funs_temp(dist_list1$c1_l.l2$data$dat$funs[smpl, ],
             col = dist_list1$c1_l.l2$data$dat$param[[1]][smpl],
             xscale = xscl) +
  ex_labs + thm
ex_data_sr2 <- plot_funs_temp(dist_list1$a2_sr.l2$data$dat$funs[smpl, ],
                          col = dist_list1$a2_sr.l2$data$param[[1]][smpl],
                          xscale = xscl) +
  ex_labs + thm
ex_data_hx3 <- plot_funs_temp(dist_list1$a3_hx.l2$data$dat$funs[smpl, ],
                          col = dist_list1$a3_hx.l2$data$param[[1]][smpl],
                          xscale = xscl) +
  ex_labs + thm

plot_grid(plotlist = list(ex_data_ap, ex_data_p2, ex_data_sr2, ex_data_hx3),
          labels = c("c1-l", "p1-l", "a2-sr", "a3-hx"),
          hjust = 0, vjust = 1, nrow = 1)


# Fig 8
### example obersvations real data

# eq data
smpl_eq <- sample(1:1558, size = 10)
plt_dat <- eq$funs[smpl_eq, ]
str(plt_dat)
plt_eq <- plot_funs_temp(plt_dat, hypo_dists[smpl_eq]) +
  xlab("Time in seconds") +
  ylab("Velocity")
plt_eq


# sp data
n_smpl <- 10
smpl <- sample(1:nrow(real), n_smpl)
funs <- real[smpl, ]

plt_sp <- plot_funs_temp(funs, lbls[smpl]) +
  xlab("Wavelength in nanometer") +
  ylab("Intensity")
plt_sp


plot_grid(plotlist = list(plt_eq, plt_sp), nrow = 1, labels = c("Eq", "Sp"), hjust = 0)

# coil data
dt_coil <- data.table::as.data.table(coil)

n_pxls <- 128
obs <- c(1, 11, 33, 50)
n_obs <- length(obs)
tt_coil <- dt_coil[obs, ]
colnames(tt_coil) <- paste0("px", 1:ncol(tt_coil))
tt_coil <- melt(tt_coil, measure.vars = colnames(tt_coil))
tt_coil[, id := rep(obs, n_pxls ^ 2)]

pixels <- expand.grid(0:(n_pxls - 1), 0:(n_pxls - 1))

tt_coil[, x := rep(pixels$Var1, each = n_obs)]
tt_coil[, y := rep(pixels$Var2, each = n_obs)]
setnames(tt_coil, "value", "intensity")

tt_coil$id_fac <- as.factor(tt_coil$id)
levels(tt_coil$id_fac) <- c("0°", "90°", "180°", "270°")

ggplot(data = tt_coil) +
  geom_tile(aes(x, y, fill = intensity)) +
  facet_wrap(~ id_fac, nrow = 6, ncol = 6)  +
  coord_flip() +
  scale_x_reverse() +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 15))
