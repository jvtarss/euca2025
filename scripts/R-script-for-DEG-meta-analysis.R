library(tximport)
library(DESeq2)
library(ggVennDiagram)
library(pheatmap)
library(tidyverse)
library(edgeR)
library(RColorBrewer)
library(metaRNASeq)
library(cowplot)

setwd("/mnt/main/euca/salmon-quantifications")
output_dir_name <- "jackknife_metaRNASeq_fishercomb_v5_direct_refactored"
dir.create(output_dir_name, showWarnings = FALSE)

alpha <- 0.05
lfc_threshold_meta <- 1
lfc_threshold_individual <- 1

samples <- data.frame(
  sample = c("WW_D22-1", "WW_D22-2", "WW_D22-3", "WS_D22-1", "WS_D22-2", "WS_D22-3",
             "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4",
             "W-K_1", "W-K_2", "W-K_3", "+W-K_1", "+W-K_2", "+W-K_3"),
  condition = factor(c(rep("Control", 3), rep("Drought", 3),
                       rep("Control", 4), rep("Drought", 4),
                       rep("Drought", 3), rep("Control", 3))),
  study = factor(c(rep("SRP406292", 6), rep("SRP458400", 8), rep("SRP178354", 6)))
)
samples$condition <- relevel(samples$condition, ref = "Control")

files <- file.path(".", samples$study, paste0(samples$sample, ".quant.sf"))
names(files) <- samples$sample
if (any(!file.exists(files))) { stop("Arquivos de quantificação faltando.") }

tx2gene <- read.delim("/mnt/main/euca/references/tx2gene.tsv")

txi_gene <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")
dge <- DGEList(counts = txi_gene$counts, samples = samples)

design_filter <- model.matrix(~study + condition, data = samples)
keep <- filterByExpr(dge, design = design_filter)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
all_counts_raw <- dge$counts
all_samples_info <- dge$samples

cpm_normalized <- cpm(dge, log = TRUE, prior.count = 2)

pdf(file.path(output_dir_name, "pcoa_plots.pdf"), width = 10, height = 8)
pcoa_all_data <- cmdscale(dist(t(cpm_normalized)), k = 2, eig = TRUE)
pcoa_all_df <- data.frame(pcoa_all_data$points, all_samples_info)
percent_explained <- pcoa_all_data$eig * 100 / sum(pcoa_all_data$eig)
print(ggplot(pcoa_all_df, aes(X1, X2, color = condition, shape = study)) +
  geom_point(size = 4, alpha = 0.8) + theme_bw(base_size = 14) +
  labs(x = paste0("PCoA1 (", round(percent_explained[1], 1), "%)"),
       y = paste0("PCoA2 (", round(percent_explained[2], 1), "%)"),
       title = "PCoA - Todas as Amostras") + scale_color_brewer(palette = "Set1"))

s_id <- "SRP406292"
study_samples_idx <- all_samples_info$study == s_id
pcoa_study_data <- cmdscale(dist(t(cpm_normalized[, study_samples_idx])), k = 2, eig = TRUE)
pcoa_study_df <- data.frame(pcoa_study_data$points, all_samples_info[study_samples_idx, ])
percent_explained_study <- pcoa_study_data$eig * 100 / sum(pcoa_study_data$eig)
print(ggplot(pcoa_study_df, aes(X1, X2, color = condition)) + geom_point(size = 4, alpha = 0.8) + theme_bw(base_size = 14) + labs(x = paste0("PCoA1 (", round(percent_explained_study[1], 1), "%)"), y = paste0("PCoA2 (", round(percent_explained_study[2], 1), "%)"), title = paste("PCoA - Estudo:", s_id)) + scale_color_brewer(palette = "Set1"))

s_id <- "SRP458400"
study_samples_idx <- all_samples_info$study == s_id
pcoa_study_data <- cmdscale(dist(t(cpm_normalized[, study_samples_idx])), k = 2, eig = TRUE)
pcoa_study_df <- data.frame(pcoa_study_data$points, all_samples_info[study_samples_idx, ])
percent_explained_study <- pcoa_study_data$eig * 100 / sum(pcoa_study_data$eig)
print(ggplot(pcoa_study_df, aes(X1, X2, color = condition)) + geom_point(size = 4, alpha = 0.8) + theme_bw(base_size = 14) + labs(x = paste0("PCoA1 (", round(percent_explained_study[1], 1), "%)"), y = paste0("PCoA2 (", round(percent_explained_study[2], 1), "%)"), title = paste("PCoA - Estudo:", s_id)) + scale_color_brewer(palette = "Set1"))

s_id <- "SRP178354"
study_samples_idx <- all_samples_info$study == s_id
pcoa_study_data <- cmdscale(dist(t(cpm_normalized[, study_samples_idx])), k = 2, eig = TRUE)
pcoa_study_df <- data.frame(pcoa_study_data$points, all_samples_info[study_samples_idx, ])
percent_explained_study <- pcoa_study_data$eig * 100 / sum(pcoa_study_data$eig)
print(ggplot(pcoa_study_df, aes(X1, X2, color = condition)) + geom_point(size = 4, alpha = 0.8) + theme_bw(base_size = 14) + labs(x = paste0("PCoA1 (", round(percent_explained_study[1], 1), "%)"), y = paste0("PCoA2 (", round(percent_explained_study[2], 1), "%)"), title = paste("PCoA - Estudo:", s_id)) + scale_color_brewer(palette = "Set1"))
dev.off()

analyze_study_deseq2 <- function(study_id, counts, info) {
  subset_info <- info[info$study == study_id, ]
  subset_counts <- round(counts[, subset_info$sample])
  if (length(unique(subset_info$condition)) < 2 || any(table(subset_info$condition) < 2)) return(NULL)
  dds <- DESeqDataSetFromMatrix(subset_counts, subset_info, ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "Drought", "Control"), alpha = alpha)
  res_shrunken <- lfcShrink(dds, contrast = c("condition", "Drought", "Control"), type = "ashr", res = res, quiet = TRUE)
  data.frame(gene = rownames(res), logFC = res_shrunken$log2FoldChange, P.Value = res$pvalue, adj.P.Val = res$padj, study = study_id, row.names = rownames(res))
}

run_meta_analysis <- function(study_ids, counts, info, iter_id) {
  results_list <- lapply(study_ids, analyze_study_deseq2, counts = counts, info = info)
  results_list <- Filter(Negate(is.null), results_list)
  names(results_list) <- sapply(results_list, function(df) unique(df$study))
  if (length(results_list) < 2) return(list(meta_res_full = NULL, significant_meta_DEGs = data.frame(), study_results = results_list))

  common_genes <- Reduce(intersect, lapply(results_list, rownames))
  if (length(common_genes) == 0) return(list(meta_res_full = NULL, significant_meta_DEGs = data.frame(), study_results = results_list))
  
  pval_list <- lapply(results_list, function(res) { pvals <- res[common_genes, "P.Value"]; pvals[is.na(pvals)] <- 1; pvals })
  fisher_res <- metaRNASeq::fishercomb(indpval = pval_list, BHth = alpha)
  
  logFC_matrix <- do.call(cbind, lapply(results_list, function(res) res[common_genes, "logFC"]))
  mean_logFC <- rowMeans(logFC_matrix, na.rm = TRUE)
  heterogeneity_var <- apply(logFC_matrix, 1, var, na.rm = TRUE)
  
  weights <- 1 / (do.call(cbind, pval_list) + 1e-10)
  weights_norm <- weights / rowSums(weights, na.rm = TRUE)
  weighted_logFC <- rowSums(logFC_matrix * weights_norm, na.rm = TRUE)

  meta_res <- data.frame(
    gene = common_genes, meta.logFC.mean = mean_logFC, meta.logFC.weighted = weighted_logFC,
    meta.P.Value = fisher_res$rawpval, meta.adj.P.Val = fisher_res$adjpval,
    heterogeneity.variance = heterogeneity_var, row.names = common_genes
  )
  sig_degs <- meta_res %>% filter(!is.na(meta.adj.P.Val), meta.adj.P.Val < alpha, abs(meta.logFC.weighted) >= lfc_threshold_meta) %>%
              mutate(direction = ifelse(meta.logFC.weighted > 0, "Up", "Down"))
  
  message(paste("Iter:", iter_id, "- DEGs:", nrow(sig_degs)))
  list(meta_res_full = meta_res, significant_meta_DEGs = sig_degs, study_results = results_list)
}

all_study_ids <- unique(all_samples_info$study)
full_analysis <- run_meta_analysis(all_study_ids, all_counts_raw, all_samples_info, "FULL")
DEGs_sig_full_meta <- full_analysis$significant_meta_DEGs

jackknife_study_sets <- lapply(all_study_ids, function(id) setdiff(all_study_ids, id))
names(jackknife_study_sets) <- paste0("Excl_", all_study_ids)
jackknife_study_sets <- Filter(function(x) length(x) >= 2, jackknife_study_sets)
jackknife_results <- lapply(names(jackknife_study_sets), function(name) {
  run_meta_analysis(jackknife_study_sets[[name]], all_counts_raw, all_samples_info, name)
})
names(jackknife_results) <- names(jackknife_study_sets)

if (nrow(DEGs_sig_full_meta) > 0) {
  jk_sig_gene_sets <- lapply(jackknife_results, function(res) res$significant_meta_DEGs$gene)
  jk_sig_gene_sets[["FULL_META"]] <- DEGs_sig_full_meta$gene
  robust_genes <- Reduce(intersect, Filter(function(x) !is.null(x) && length(x) > 0, jk_sig_gene_sets))
  DEGs_robustos_jackknife_df <- full_analysis$meta_res_full[robust_genes, ] %>%
    mutate(direction = ifelse(meta.logFC.weighted > 0, "Up", "Down"))
} else {
  DEGs_robustos_jackknife_df <- data.frame()
}

all_individual_results_df <- bind_rows(full_analysis$study_results)
common_degs_summary <- all_individual_results_df %>%
  filter(!is.na(adj.P.Val) & adj.P.Val < alpha, abs(logFC) >= lfc_threshold_individual) %>%
  mutate(direction = sign(logFC)) %>%
  group_by(gene) %>%
  summarise(n_studies = n(), n_directions = n_distinct(direction), .groups = 'drop') %>%
  filter(n_studies == length(all_study_ids), n_directions == 1)

DEGs_comuns_final_df <- DEGs_robustos_jackknife_df %>%
  filter(gene %in% common_degs_summary$gene)

write.csv(DEGs_robustos_jackknife_df, file.path(output_dir_name, "meta_analysis_ROBUST_DEGs.csv"))
write.csv(DEGs_comuns_final_df, file.path(output_dir_name, "meta_analysis_COMMON_DEGs.csv"))
write.csv(DEGs_sig_full_meta, file.path(output_dir_name, "meta_analysis_SIGNIFICANT_DEGs.csv"))

venn_list_studies <- lapply(full_analysis$study_results, function(res) {
  res %>% filter(!is.na(adj.P.Val), adj.P.Val < alpha, abs(logFC) >= lfc_threshold_individual) %>% pull(gene)
})
venn_list_studies[["MetaAnalise"]] <- DEGs_sig_full_meta$gene
if (sum(sapply(venn_list_studies, length)) > 0) {
  p_venn_s <- ggVennDiagram(venn_list_studies, label = "count") + scale_fill_distiller(palette = "Blues", direction = 1)
  ggsave(file.path(output_dir_name, "venn_studies_vs_meta.pdf"), plot = p_venn_s, width = 10, height = 8)
}

jaccard_data <- bind_rows(lapply(names(jackknife_results), function(iter_name) {
  jk_degs <- jackknife_results[[iter_name]]$significant_meta_DEGs$gene
  full_degs <- DEGs_sig_full_meta$gene
  intersect_val <- length(intersect(full_degs, jk_degs))
  union_val <- length(union(full_degs, jk_degs))
  jaccard_dissim <- if (union_val == 0) 0 else 1 - (intersect_val / union_val)
  data.frame(Iteration = iter_name, Jaccard_Dissimilarity = jaccard_dissim, 
             Genes_Lost = length(setdiff(full_degs, jk_degs)), 
             Genes_Gained = length(setdiff(jk_degs, full_degs)))
}))
write.csv(jaccard_data, file.path(output_dir_name, "jackknife_jaccard_dissimilarity.csv"), row.names = FALSE)

deg_summary_list <- c(
  lapply(names(full_analysis$study_results), function(s_name) {
    degs <- full_analysis$study_results[[s_name]] %>% filter(!is.na(adj.P.Val), adj.P.Val < alpha, abs(logFC) >= lfc_threshold_individual)
    data.frame(Analysis = s_name, Total = nrow(degs), Up = sum(degs$logFC > 0), Down = sum(degs$logFC < 0))
  }),
  list(
    FullMeta = data.frame(Analysis = "FullMeta", Total = nrow(DEGs_sig_full_meta), Up = sum(DEGs_sig_full_meta$direction == "Up"), Down = sum(DEGs_sig_full_meta$direction == "Down")),
    RobustJK = data.frame(Analysis = "RobustJK", Total = nrow(DEGs_robustos_jackknife_df), Up = sum(DEGs_robustos_jackknife_df$direction == "Up"), Down = sum(DEGs_robustos_jackknife_df$direction == "Down")),
    CommonAll = data.frame(Analysis = "CommonAll", Total = nrow(DEGs_comuns_final_df), Up = sum(DEGs_comuns_final_df$direction == "Up"), Down = sum(DEGs_comuns_final_df$direction == "Down"))
  )
)
deg_summary_final_df <- bind_rows(deg_summary_list)
write.csv(deg_summary_final_df, file.path(output_dir_name, "deg_summary_report.csv"), row.names = FALSE)

if (nrow(DEGs_robustos_jackknife_df) > 1) {
  sig_degs_presence_list <- c(list(FULL_META = DEGs_sig_full_meta$gene), lapply(jackknife_results, function(res) res$significant_meta_DEGs$gene))
  presence_matrix <- sapply(sig_degs_presence_list, function(gene_set) as.integer(DEGs_robustos_jackknife_df$gene %in% gene_set))
  rownames(presence_matrix) <- DEGs_robustos_jackknife_df$gene
  
  pheatmap(presence_matrix, 
           filename = file.path(output_dir_name, "heatmap_robust_degs_presence.pdf"),
           main = "Presence of Robust DEGs in Meta-Analyses",
           cluster_rows = (nrow(presence_matrix) > 1), 
           cluster_cols = (ncol(presence_matrix) > 1),
           show_rownames = (nrow(presence_matrix) <= 50),
           color = c("white", "darkblue"))
}