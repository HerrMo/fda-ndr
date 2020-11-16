## plots
source("R/setup.R")

# optimal embedding results for simulated and real data




## SIMULATED DATA

## linear settings ------------------------

mds_lin <- l_embs_lin[str_detect(names(l_embs_lin), "mds")][c(1, 5, 9, 13, 17, 21)]
umap_lin <- l_embs_lin[str_detect(names(l_embs_lin), "umap")][c(1, 5, 9, 13, 17, 21)]
tsne_lin <- l_embs_lin[str_detect(names(l_embs_lin), "tsne")][c(4, 6, 9, 13, 17, 21)] # 21 with 320° and 13 with 160°
imap_lin <- l_embs_lin[str_detect(names(l_embs_lin), "isomap")][c(1, 5, 9, 13, 17, 21)] # 17 with 50°
dmap_lin <- l_embs_lin[str_detect(names(l_embs_lin), "dmap")][c(1, 5, 11, 13, 19, 21)] # 19 with 280° and 21 with 260°
lin_embs <- c(mds_lin, umap_lin, tsne_lin, imap_lin, dmap_lin)

plot_grid3d(lin_embs,
            dist_list1,
            params = list(angle = list(vals = c(160, 320, 50, 280, 260),
                                       pos = c(16, 18, 23, 29, 30))),
            ylim2 = list(pos = c(1, 19, 20, 21),
                         vals1 = c(-45, -45, -30, -75),
                         vals2 = c(45, 45, 30, 75)),
            cex.lab = 0.75,
            pch = 20,
            ncol = 6,
            subtitle = FALSE)

nic_meths <- c("MDS", "UMAP", "t-SNE", "ISOMAP", "DIFFMAP")
adjs_y1 <- c(0.95, 0.7, 0.5, 0.28, 0.05)
for (i in 1:5) mtext(nic_meths[i], side = 2, font = 2, line = 1, outer = TRUE, adj = adjs_y1[i], cex = 1.5)
adjs_x1 <- c(0.03, 0.25, 0.4, 0.6, 0.75, 0.9)
for (i in 1:6) mtext(nic_settings[i], line = 1, font = 2, outer = TRUE, adj = adjs_x1[i], cex = 2)



## nonlinear settings ---------------------
embs_nl_tuid <- l_embs_nl[str_detect(names(l_embs_nl), "mds", negate = TRUE)]
embs_nl_mds <- l_embs_nl[str_detect(names(l_embs_nl), "mds")]

### mds embeddings for nonlinear settings
plot_grid3d(embs_nl_mds[c(1, 5, 10, 15, 20)], dist_list1, pch = 20, subtitle = FALSE)

mtext(c("MDS"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.5, cex = 2) # mds only
adjs_x2 <- c(0.07, 0.275, 0.5, 0.7, 0.9)
for (i in 1:5) mtext(nic_settings[-c(1:6)][i], line = 1, font = 2, outer = TRUE, adj = adjs_x2[i], cex = 2)


### other embeddings for nonlinear settings
# for outer margin annotations of methods and settings run 'annotation' (l 88)

# Fig 5 A
# 3d plots for nonlinear settings based on parameter space optimization via AUC-dir
plot_grid3d(embs_nl_tuid[str_detect(names(embs_nl_tuid), "ps_l2")],
            dist_list1,
            params = list(angle = list(vals = c(240, 160, 60), pos = c(1, 10, 3))),
            cex.lab = 0.75,
            pch = 20)

# Fig 5 B
# 3d plots for nonlinear settings based on parameter space optimization via AUC-geo
plot_grid3d(embs_nl_tuid[str_detect(names(embs_nl_tuid), "ps_geo")], dist_list1,
            ylim2 = list(pos = 11, vals1 = -1000, vals2 = 1000), cex.lab = 0.75, pch = 20)

# Fig 6 A
# 3d plots for nonlinear settings based on function space optimization via AUC-dir
plot_grid3d(embs_nl_tuid[str_detect(names(embs_nl_tuid), "fs_l2")], dist_list1, cex.lab = 0.75, pch = 16)

# Fig 6 B
# 3d plots for nonlinear settings based on function space optimization via AUC-geo
plot_grid3d(embs_nl_tuid[str_detect(names(embs_nl_tuid), "fs_geo")],
            dist_list1,
            cex.lab = 0.001,
            pch = 20,
            ylim2 = list(pos = c(11, 12), vals1 = c(-1000, -400), vals2 = c(1000, 400)),
            params = list(angle = list(vals = c(160, 200), pos = c(5, 13))))

# annotation
for (i in 1:5) mtext(nic_settings[-c(1:6)][i], line = 1, font = 2, outer = TRUE, adj = adjs_x2[i], cex = 2)
adjs_y2 <- c(0.95, 0.67, 0.42, 0.1)
for (i in 1:4) mtext(nic_meths[-1][i], side = 2, font = 2, line = 1, outer = TRUE, adj = adjs_y2[i], cex = 2)


# Fig. 7
# Dotplot
dpPal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


dp_dat <-
  res_sim2 %>%
  filter(!meas %in% c("q_global", "q_local")) %>%
  filter(data %in% settings[-c(1:6)]) %>%
  filter(method != "mds") %>%
  group_by(data, method, space) %>%
  top_n(1, value) %>%
  select(one_of(c("meas", "data", "method", "space", "value")))

dp_dat$data <- as.factor(dp_dat$data)
dp_dat$space <- as.factor(dp_dat$space)
dp_dat$method <- as.factor(dp_dat$method)
levels(dp_dat$data) <- c("a2-sr", "a3-hx", "a3-sc", "a3-sr", "a3-tp")
levels(dp_dat$space) <- c(TeX("$AUC^{geo}_{Rnx}$ (fs)"),
                      TeX("$AUC^{dir}_{Rnx}$ (fs)"),
                      TeX("$AUC^{geo}_{Rnx}$ (ps)"),
                      TeX("$AUC^{dir}_{Rnx}$ (ps)"))
levels(dp_dat$method) <- c("DIFFMAP", "ISOMAP", "t-SNE", "UMAP")

# geom_jitter: position on "yaxis" might be slightly different
ggplot(dp_dat, aes(x = value, y = data, color = method, shape = method)) +
  geom_jitter(size = 3) +
  facet_wrap(~ space,
             labeller = label_parsed) +
  scale_colour_manual(values = dpPal, name = "Method") +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        legend.position = "top",
        axis.line = element_blank(),
        strip.text.x = element_text(color = "black", size = 15),
        strip.background = element_rect(fill = "white"),
        axis.text.y = element_text(size = 15)) +
  scale_shape_discrete(name = "Method")

# Table 3
res_sim2 %>%
  filter(data %in% settings[-c(1:6)]) %>% # only nonlinear settings
  filter(data != "a3_sr.l2") %>% # remove because no good embedding possible
  filter(meas == "auc_rnx") %>%
  filter(!method %in% "mds") %>%
  group_by(data, method, space) %>%
  top_n(1, value) %>%
  select(value) %>%
  separate(space, c("space", "metric"), sep = "_") %>%
  group_by(data, method, metric) %>%
  summarize(diff = abs(diff(range(value)))) %>%
  group_by(method, metric) %>%
  summarize(md = mean(diff))


###### REAL DATA

# new order: mds, umap, tsne, imap, diffmap to match sim data ordering
reordering <- c(5, 1:2, 4, 3, 10, 6:7, 9, 8, 15, 11:12, 14, 13, 20, 16:17, 19, 18)

##### Earthquake data
plot_grid3d(l_embs_earthquake[reordering],
            color = hypo_dists,
            cex.lab = 0.75,
            pch = 20,
            subtitle = TRUE,
            params = list(angle = list(vals = c(rep(260, 4), # mds
                                                260, 260, 260, 260, # isomap
                                                100, 260, 300, 240, # umap
                                                300, 200, 120, 265 # tsne
            ),
            pos = c(1, 6, 11, 16,
                    4, 9, 14, 19,
                    2, 7, 12, 17,
                    3, 8, 13, 18
            ))),
            box = FALSE,
            ncol = 4,
            bycol = TRUE)

mtext(c("t-SNE"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.51, cex = 1.5)
mtext(c("DIFFMAP"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.05, cex = 1.5)
mtext(c("ISOMAP"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.29, cex = 1.5)
mtext(c("UMAP"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.73, cex = 1.5)
mtext(c("MDS"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.93, cex = 1.5)

##### COIL data
clPal <- colorRampPalette(rainbow(15))

plot_grid3d(l_embs_coil[reordering],
            color = 1:72,
            cex.lab = 0.75,
            pch = 20,
            subtitle = TRUE,
            params = list(angle = list(vals = c(rep(120, 4), # mds
                                                rep(120, 4), # isomap
                                                265, 265, # umap
                                                265, 265 # tsne
            ),
            pos = c(1, 6, 11, 16,
                    4, 9, 14, 19,
                    12, 17,
                    13, 18
            ))),
            box = FALSE,
            ncol = 4,
            bycol = TRUE,
            color_pal = clPal)

mtext(c("t-SNE"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.51, cex = 1.5)
mtext(c("DIFFMAP"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.05, cex = 1.5)
mtext(c("ISOMAP"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.29, cex = 1.5)
mtext(c("UMAP"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.73, cex = 1.5)
mtext(c("MDS"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.93, cex = 1.5)


##### Spectro data
spPal <- colorRampPalette(brewer.pal(12,"Paired"))
plot_grid3d(l_embs_spectro[reordering],
            color = 1:1004,
            cex.lab = 0.75,
            pch = 20,
            subtitle = TRUE,
            params = list(angle = list(vals = c(rep(270, 4), # mds
                                                270, 270, 260, 270, # isomap
                                                300, 260, 245, 240, # umap
                                                240, 240, 100, 265 # tsne
            ),
            pos = c(1, 6, 11, 16,
                    4, 9, 14, 19,
                    2, 7, 12, 17,
                    3, 8, 13, 18
            ))),
            box = FALSE,
            ncol = 4,
            bycol = TRUE,
            color_pal = spPal)

mtext(c("t-SNE"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.51, cex = 1.5)
mtext(c("DIFFMAP"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.05, cex = 1.5)
mtext(c("ISOMAP"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.29, cex = 1.5)
mtext(c("UMAP"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.73, cex = 1.5)
mtext(c("MDS"), side = 2, font = 2, line = 1, outer = TRUE, adj = 0.93, cex = 1.5)

