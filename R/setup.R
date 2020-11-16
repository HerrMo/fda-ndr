library(ggplot2)
library(igraph) # needed for diffmap
library(Matrix) # -"-
library(funData) # needed for earthquake data
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(scatterplot3d)
library(latex2exp)
library(cowplot)

source("R/ancillary-code_results.R")
source("R/ancillary-code_experiments.R")


# data and methods
meths <- c("umap", "tsne", "dmap", "isomap",  "mds")
nic_meths <- c("MDS", "UMAP", "t-SNE", "ISOMAP", "DIFFMAP")
dat_sets <- c("a1_l", "p1_l", "c1_l", "a2_l", "p2_l", "i2_l", "a2_sr", "a3_hx", "a3_sr", "a3_sc", "a3_tp")


# parameter settings and performance values of sim data (tuning results)
load("data/res_tuning_sim_data.RData")
setnames(res_sim, "fs_dist", "fs_l2")
res_sim2 <- pivot_longer(res_sim, c("fs_l2", "fs_geo", "ps_l2", "ps_geo"), names_to = "space")

# parameter settings and performance values of real data (tunin results)
load("data/res_tuning_real_spectro.RData")
res_sp <- res
load("data/res_tuning_real_earthquake.RData")
load("data/res_tuning_real_coil.RData")


# distance matrices and raw simulated fun data
load("data/distance_matrices1.RData")
settings <- paste(rep(dat_sets, each = 1), rep("l2", 1), sep = ".")
nic_settings <- sapply(strsplit(settings, "\\."), function(x) x[1])
nic_settings <- sub("_", "-", nic_settings)
names(dist_list1) <- settings
# param 1 is constant, but ususally used for plotting, param 2 varies
dist_list1$p1_l.l2$data$param[1] <- dist_list1$p1_l.l2$data$param[2]


# preprocessing real data

# eq data

eq1 <- readRDS("data/seissol_tweedie_smooth.rds")
# extract funData object as in Happ et al.
ind <- which(eq1$hypo.dist < 40000)
x <- seq(0, 30, by = .5)
eq <- funData(argvals = x, X = as.matrix(log1p(eq1$Bodenbewegung[ind, ])))
# Convert to fun data object as used for this setting
eq <- list(funs = eq@X, grid = eq@argvals[[1]])
eq_l2 <- as.matrix(dist(eq$funs))
hypo_dists <- eq1$hypo.dist[ind]

# sp data

real_train <- read.table("data/EthanolLevel_TRAIN.txt")
real_test <- read.table("data/EthanolLevel_TEST.txt")
real <- as.matrix(rbind(real_train, real_test)[, -1])
lbls <-  rbind(real_train, real_test)[[1]]
real_l2 <- as.matrix(dist(real))

# coil data

load("data/coil20.RData")
coil <- as.matrix(coil20[coil20$Label == 1, -16385])

coil_ids <- attributes(coil)[[2]][[1]]
coil_ids <- str_split(coil_ids, "_", simplify = TRUE)
padded_ids <- leftpad(coil_ids[, 2])

coil_correct <- coil[order(padded_ids),]
coil_l2 <- as.matrix(dist(coil_correct))


# default color palette
rbPal <- colorRampPalette(c("yellow", "blue"))


### loading stored optimal embeddings

# data for linear and nonlinear embeddings splitted to be uploadable to github
load("data/embs_nonlinear_01.RData")
load("data/embs_nonlinear_02.RData")
load("data/embs_nonlinear_03.RData")
l_embs_nl <- c(l_embs_nl1, l_embs_nl2, l_embs_nl3)
rm(list =  c("l_embs_nl1", "l_embs_nl2", "l_embs_nl3"))

load("data/embs_linear_01.RData")
load("data/embs_linear_02.RData")
load("data/embs_linear_03.RData")
load("data/embs_linear_04.RData")
l_embs_lin <- c(l_embs_lin1, l_embs_lin2, l_embs_lin3, l_embs_lin4)
rm(list = c("l_embs_lin1", "l_embs_lin2", "l_embs_lin3", "l_embs_lin4"))

load("data/embs_spectro.RData")
load("data/embs_earthquake.RData")
load("data/embs_coil.RData")
