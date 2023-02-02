# Ancom2分析 https://github.com/FrederickHuangLin/ANCOMBC
# https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
# https://www.bioconductor.org/packages/release/bioc/html/ANCOMBC.html
library(ANCOMBC)
library(phyloseq)
library(tidyverse)

# 构建phyloseq数据结构
map <- read.table("data/metadata.txt", header=T, row.names=1, sep="\t", 
                 comment.char="", stringsAsFactors=F)
otu <- read.table("data/otutab.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F)
tax <- read.table("data/taxonomy.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F)
tse <- phyloseq(sample_data(map),otu_table(as.matrix(otu), taxa_are_rows=TRUE), tax_table(as.matrix(tax)))
tse$Group = factor(tse$Group, levels = c("obese", "overweight", "lean"))
# run ancom function

output = ancombc2(data = tse, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "Group", rand_formula = NULL,
                  p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "Group", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(-1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2, 2),
                                       solver = "ECOS",
                                       B = 100))

res_prim = output$res

# Original results
df = res_prim %>%
  dplyr::select(taxon, contains(group)) 
df_fig = df %>%
  filter(diff_GroupOE == 1 | diff_GroupWT == 1) %>%
  mutate(lfc_OE = ifelse(diff_GroupOE == 1, 
                                 lfc_GroupOE, 0),
         lfc_WT = ifelse(diff_GroupWT == 1, 
                           lfc_GroupWT, 0)) %>%
  transmute(taxon, 
            `OE vs. KO` = round(lfc_OE, 2),
            `WT vs. KO` = round(lfc_WT, 2)) %>%
  pivot_longer(cols = `OE vs. KO`:`WT vs. KO`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

lo = floor(min(df_fig$value))
up = ceiling(max(df_fig$value))
mid = (lo + up)/2
fig = df_fig %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to KO subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("", fig, )

# Pseudo-count sensitivity analysis
tab_sens <- output$pseudo_sens_tab

sens = tab_sens %>%
  transmute(taxon, sens_group = GroupOE) %>%
  left_join(df %>%
              transmute(taxon, diff_group = diff_GroupOE), 
            by = "taxon") %>%
  mutate(group = "OE vs. KO") %>%
  bind_rows(
    tab_sens %>%
      transmute(taxon, sens_group = GroupWT) %>%
      left_join(df %>%
                  transmute(taxon, diff_group = diff_GroupWT), 
                by = "taxon") %>%
      mutate(group = "WT vs. KO")
  )
sens$diff_group = recode(sens$diff_group * 1, 
                           `1` = "Significant",
                           `0` = "Nonsignificant")

fig_sens = sens %>%
  ggplot(aes(x = taxon, y = sens_group, color = diff_group)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  facet_grid(rows = vars(group), scales = "free") +
  labs(x = NULL, y = "Sensitivity Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
ggsave("", fig_sens, )

# global test
res_global = output$res_global

df = res_prim %>%
  dplyr::select(taxon, contains(group)) 
df_fig_global = df %>%
  left_join(res_global %>%
              transmute(taxon, diff_group = diff_abn)) %>%
  dplyr::filter(diff_group == 1) %>%
  mutate(lfc_OE = ifelse(diff_GroupOE == 1, 
                         lfc_GroupOE, 0),
         lfc_WT = ifelse(diff_GroupWT == 1, 
                         lfc_GroupWT, 0)) %>%
  transmute(taxon, 
            `OE vs. KO` = round(lfc_OE, 2),
            `WT vs. KO` = round(lfc_WT, 2)) %>%
  pivot_longer(cols = `OE vs. KO`:`WT vs. KO`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

lo = floor(min(df_fig_global$value))
up = ceiling(max(df_fig_global$value))
mid = (lo + up)/2
fig = df_fig_global %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes for globally significant taxa") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("", fig, )

# Pseudo-count sensitivity analysis for the global test
sens_global = tab_sens %>%
  transmute(taxon, sens_global = global) %>%
  left_join(res_global %>%
              transmute(taxon, diff_global = diff_abn * 1), 
            by = "taxon") 
sens_global$diff_global = recode(sens_global$diff_global, 
                                 `1` = "Significant",
                                 `0` = "Nonsignificant")

fig_sens_global = sens_global %>%
  ggplot(aes(x = taxon, y = sens_global, color = diff_global)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  labs(x = NULL, y = "Sensitivity Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
ggsave("", fig_sens_global, )

# pairwise directional test
res_pair = output$res_pair

df_fig_pair = res_pair %>%
  filter(diff_GroupOE == 1 | diff_GroupWT == 1 |
           diff_GroupWT_GroupOE == 1) %>%
  mutate(lfc_OE = ifelse(diff_GroupOE == 1, 
                           lfc_GroupOE, 0),
         lfc_WT = ifelse(diff_GroupWT == 1, 
                                 lfc_GroupWT, 0),
         lfc_WT_OE = ifelse(diff_GroupWT_GroupOE == 1, 
                                      lfc_GroupWT_GroupOE, 0)) %>%
  transmute(taxon, 
            `OE vs. KO` = round(lfc_OE, 2), 
            `WT vs. KO` = round(lfc_WT, 2),
            `WT vs. OE` = round(lfc_WT_OE, 2)
  ) %>%
  pivot_longer(cols = `OE vs. KO`:`WT vs. OE`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)
df_fig_pair$group = factor(df_fig_pair$group, 
                           levels = c("OE vs. KO",
                                      "WT vs. KO",
                                      "WT vs. OE"))

lo = floor(min(df_fig_pair$value))
up = ceiling(max(df_fig_pair$value))
mid = (lo + up)/2
fig_pair = df_fig_pair %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold change of pairwise comparisons") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Pseudo-count sensitivity analysis
sens_pair = tab_sens %>%
  transmute(taxon, sens_pair = `KO - OE`) %>%
  left_join(res_pair %>%
              transmute(taxon, 
                        diff_pair = diff_GroupOE), 
            by = "taxon") %>%
  mutate(group = "KO vs. OE") %>%
  bind_rows(
    tab_sens %>%
      transmute(taxon, sens_pair = `KO - WT`) %>%
      left_join(res_pair %>%
                  transmute(taxon, 
                            diff_pair = diff_GroupWT), 
                by = "taxon") %>%
      mutate(group = "KO vs. WT")
  ) %>%
  bind_rows(
    tab_sens %>%
      transmute(taxon, sens_pair = `OE - WT`) %>%
      left_join(res_pair %>%
                  transmute(taxon, 
                            diff_pair = diff_GroupWT_GroupOE), 
                by = "taxon") %>%
      mutate(group = "WT vs. OE")
  )
sens_pair$diff_pair = recode(sens_pair$diff_pair * 1, 
                             `1` = "Significant",
                             `0` = "Nonsignificant")
sens_pair$group = factor(sens_pair$group, 
                         levels = c("KO vs. OE",
                                    "KO vs. WT",
                                    "WT vs. OE"))

fig_sens_pair = sens_pair %>%
  ggplot(aes(x = taxon, y = sens_pair, color = diff_pair)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  facet_grid(rows = vars(group), scales = "free") +
  labs(x = NULL, y = "Sensitivity Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
fig_sens_pair

# Dunnett’s type of test,  is designed for making comparisons 
# of several experimental groups with the control or the reference group.

# Of note is that the ANCOM-BC2 primary results 
# do not control the mdFDR for the comparison of multiple groups.
res_dunn = output$res_dunn

df_fig_dunn = res_dunn %>%
  filter(diff_GroupOE == 1 | diff_GroupWT == 1) %>%
  mutate(lfc_OE = ifelse(diff_GroupOE == 1, 
                         lfc_GroupOE, 0),
         lfc_WT = ifelse(diff_GroupWT == 1, 
                         lfc_GroupWT, 0)) %>%
  transmute(taxon, 
            `OE vs. KO` = round(lfc_OE, 2),
            `WT vs. KO` = round(lfc_WT, 2)) %>%
  pivot_longer(cols = `OE vs. KO`:`WT vs. KO`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

lo = floor(min(df_fig_dunn$value))
up = ceiling(max(df_fig_dunn$value))
mid = (lo + up)/2
fig_dunn = df_fig_dunn %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_dunn

sens_dunn = tab_sens %>%
  transmute(taxon, sens_group = GroupOE) %>%
  left_join(res_dunn %>%
              transmute(taxon, diff_group = diff_GroupOE), 
            by = "taxon") %>%
  mutate(group = "OE vs. KO") %>%
  bind_rows(
    tab_sens %>%
      transmute(taxon, sens_group = GroupWT) %>%
      left_join(df %>%
                  transmute(taxon, diff_group = diff_GroupWT), 
                by = "taxon") %>%
      mutate(group = "WT vs. KO")
  )
sens_dunn$diff_group = recode(sens_dunn$diff_group * 1, 
                             `1` = "Significant",
                             `0` = "Nonsignificant")

fig_sens_dunn = sens_dunn %>%
  ggplot(aes(x = taxon, y = sens_group, color = diff_group)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  facet_grid(rows = vars(group), scales = "free") +
  labs(x = NULL, y = "Sensitivity Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
fig_sens_dunn

# trend test
res_trend = output$res_trend

df_fig_trend = res_trend %>%
  dplyr::filter(diff_abn) %>%
  transmute(taxon,
            lfc = lfc_GroupOE,
            se = se_GroupOE,
            q_val,
            group = "OE - KO") %>%
  bind_rows(
    res_trend %>%
      dplyr::filter(diff_abn) %>%
      transmute(taxon,
                lfc = lfc_GroupWT,
                se = se_GroupOE,
                q_val,
                group = "WT - KO")
  )

df_fig_trend$group = factor(df_fig_trend$group, 
                            levels = c("OE - KO", "WT - KO"))

fig_trend = df_fig_trend %>%
  ggplot(aes(x = group, y = lfc, fill = group)) + 
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = .2,
                position = position_dodge(.9)) +
  facet_wrap(vars(taxon), nrow = 2, scales = "free") +
  labs(x = NULL, y = NULL, title = "Log fold change as compared to obese subjects") +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.9, 0.2))
fig_trend

sens_trend = tab_sens %>%
  transmute(taxon, sens_group = GroupOE) %>%
  left_join(res_trend %>%
              transmute(taxon, diff_group = diff_abn), 
            by = "taxon") %>%
  mutate(group = "OE vs. KO") %>%
  bind_rows(
    tab_sens %>%
      transmute(taxon, sens_group = GroupWT) %>%
      left_join(res_trend %>%
                  transmute(taxon, diff_group = diff_abn), 
                by = "taxon") %>%
      mutate(group = "WT vs. KO")
  )
sens_trend$diff_group = recode(sens_trend$diff_group * 1, 
                              `1` = "Significant",
                              `0` = "Nonsignificant")

fig_sens_trend = sens_trend %>%
  ggplot(aes(x = taxon, y = sens_group, color = diff_group)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  facet_grid(rows = vars(group), scales = "free") +
  labs(x = NULL, y = "Sensitivity Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
fig_sens_trend
