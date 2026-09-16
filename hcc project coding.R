library(Seurat)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(ggrepel)

seu <- readRDS("D:/hcc_project/data/processed/seu_umap.rds")
cat("تم تحميل seu بنجاح\n")



res_df <- res_df %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "أعلى في النقائل",
      padj < 0.05 & log2FoldChange < -1 ~ "أعلى في PT",
      TRUE ~ "غير معنوي"
    )
  )

table(res_df$significance)



fig_dir <- "D:/hcc_project/results/figures/"
existing_figs <- list.files(fig_dir, pattern = "^[0-9]+_.*\\.png$")
nums <- as.integer(sub("^([0-9]+)_.*", "\\1", existing_figs))
next_num <- max(nums, na.rm = TRUE) + 1
fig_num <- sprintf("%02d", next_num)
fig_path <- file.path(fig_dir, paste0(fig_num, "_volcano_DE_hepatocyte.png"))






res_df <- res_df %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "أعلى في النقائل",
      padj < 0.05 & log2FoldChange < -1 ~ "أعلى في PT",
      TRUE ~ "غير معنوي"
    )
  )

table(res_df$significance)

top_genes <- res_df %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(10)

fig_dir <- "D:/hcc_project/results/figures/"
existing_figs <- list.files(fig_dir, pattern = "^[0-9]+_.*\\.png$")
nums <- as.integer(sub("^([0-9]+)_.*", "\\1", existing_figs))
next_num <- max(nums, na.rm = TRUE) + 1
fig_num <- sprintf("%02d", next_num)
fig_path <- file.path(fig_dir, paste0(fig_num, "_volcano_DE_hepatocyte.png"))

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("أعلى في النقائل" = "red", 
                                "أعلى في PT" = "blue", 
                                "غير معنوي" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(data = top_genes, aes(label = gene), 
                           size = 3, color = "black", max.overlaps = 15) +
  labs(title = "DE Genes: Hepatocyte Cells (Metastasis vs PT)",
       x = "log2 Fold Change", y = "-log10(adjusted p-value)",
       color = "الدلالة") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(p_volcano)

ggsave(filename = fig_path, plot =

         
         
         log_path <- "D:/hcc_project/results/project_log.txt"
       log_entry <- paste0(
         "[", Sys.time(), "] ",
         "Created Volcano plot for DE genes (Hepatocyte, Metastasis vs PT). ",
         "All significant genes downregulated in metastasis (liver function genes). ",
         "Saved: ", fig_path, "\n"
       )
       cat(log_entry, file = log_path, append = TRUE)
       cat("تم تسجيل الخطوة في project_log.txt\n")

       
       
      
       fig_dir <- "D:/hcc_project/results/figures/"
       fig_path2 <- file.path(fig_dir, paste0(fig_num, "_volcano_DE_hepatocyte_v2.png"))
       
       p_volcano2 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
         geom_point(alpha = 0.7, size = 2) +
         scale_color_manual(values = c(
           "أعلى في النقائل" = "#D62728",   # أحمر واضح
           "أعلى في PT"      = "#1F77B4",   # أزرق غامق واضح
           "غير معنوي"       = "#C7C7C7"    # رمادي فاتح خفيف
         )) +
         geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.4) +
         geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
         ggrepel::geom_text_repel(
           data = top_genes, aes(label = gene),
           size = 3.5, color = "black", fontface = "bold",
           max.overlaps = 15, box.padding = 0.4
         ) +
         labs(
           title = "الجينات المختلفة التعبير: خلايا Hepatocyte (النقائل مقابل PT)",
           x = "log2 Fold Change",
           y = "-log10(adjusted p-value)",
           color = "الدلالة"
         ) +
         theme_minimal(base_size = 13) +
         theme(
           plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
           legend.position = "right",
           panel.grid.minor = element_blank()
         )
       
       print(p_volcano2)
       
       ggsave(filename = fig_path2, plot = p_volcano2, width = 9, height = 7, dpi = 300)
       cat("تم حفظ النسخة المحسّنة في:", fig_path2, "\n")
       
       
       
       
       # ============================================
       # رسمة Volcano بلوحة ألوان جديدة (أزرق بحري + روزي)
       # ============================================
       fig_dir <- "D:/hcc_project/results/figures/"
       fig_path3 <- file.path(fig_dir, paste0(fig_num, "_volcano_DE_hepatocyte_v3.png"))
       
       p_volcano3 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
         geom_point(alpha = 0.75, size = 2) +
         scale_color_manual(values = c(
           "أعلى في النقائل" = "#E4557A",   # روزي قوي
           "أعلى في PT"      = "#00637A",   # أزرق بحري (petrol)
           "غير معنوي"       = "#D9D9D9"    # رمادي فاتح
         )) +
         geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.4) +
         geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
         ggrepel::geom_text_repel(
           data = top_genes, aes(label = gene),
           size = 3.5, color = "black", fontface = "bold",
           max.overlaps = 15, box.padding = 0.4
         ) +
         labs(
           title = "الجينات المختلفة التعبير: خلايا Hepatocyte (النقائل مقابل PT)",
           x = "log2 Fold Change",
           y = "-log10(adjusted p-value)",
           color = "الدلالة"
         ) +
         theme_minimal(base_size = 13) +
         theme(
           plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
           legend.position = "right",
           panel.grid.minor = element_blank()
         )
       
       print(p_volcano3)
       
       ggsave(filename = fig_path3, plot = p_volcano3, width = 9, height = 7, dpi = 300)
       cat("تم حفظ النسخة بالألوان الجديدة في:", fig_path3, "\n")
       
       # ============================================
       # تسجيل الخطوة في اللوج
       # ============================================
       log_path <- "D:/hcc_project/results/project_log.txt"
       log_entry <- paste0(
         "[", Sys.time(), "] ",
         "Finalized Volcano plot color scheme (petrol blue +
         
         
         
         
         # ============================================
# رسمة Volcano بألوان باستيل (لافندر + مينتي) - بالإنجليزي بالكامل
# ============================================
fig_dir <- "D:/hcc_project/results/figures/"
fig_path4 <- file.path(fig_dir, paste0(fig_num, "_volcano_DE_hepatocyte_v4.png"))

# نعيد تسمية عمود significance بالإنجليزي
res_df <- res_df %>%
  mutate(
    significance_en = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up in Metastasis",
      padj < 0.05 & log2FoldChange < -1 ~ "Up in PT",
      TRUE ~ "Not Significant"
    )
  )

p_volcano4 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance_en)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_color_manual(values = c(
    "Up in Metastasis" = "#C9A0DC",   # لافندر
         "Up in PT"         = "#8FD9C4",   # مينتي فاتح
         "Not Significant"  = "#E0E0E0"    # رمادي فاتح جدًا
       )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = top_genes, aes(label = gene),
    size = 3.5, color = "black", fontface = "bold",
    max.overlaps = 15, box.padding = 0.4
  ) +
  labs(
    title = "Differentially Expressed Genes: Hepatocyte Cells (Metastasis vs PT)",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "Significance"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_volcano4)

ggsave(filename = fig_path4, plot = p_volcano4, width = 9, height = 7, dpi = 300)
cat("Saved pastel version to:", fig_path4, "\n")

# ============================================
# تسجيل الخطوة في اللوج
# ============================================
log_path <- "D:/hcc_project/results/project_log.txt"
log_entry <- paste0(
  "[", Sys.time(), "] ",
  "Finalized Volcano plot with pastel color scheme (lavender + mint) and English labels. ",
  "Saved: ", fig_path4, "\n"
)
cat(log_entry, file = log_path, append = TRUE)
cat("تم تسجيل الخطوة في project_log.txt\n")























fig_dir <- "D:/hcc_project/results/figures/"
fig_path4 <- file.path(fig_dir, "10_volcano_DE_hepatocyte_v4.png")

res_df <- res_df %>%
  mutate(
    significance_en = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up in Metastasis",
      padj < 0.05 & log2FoldChange < -1 ~ "Up in PT",
      TRUE ~ "Not Significant"
    )
  )

p_volcano4 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance_en)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_color_manual(values = c(
    "Up in Metastasis" = "#C9A0DC",
    "Up in PT"         = "#8FD9C4",
    "Not Significant"  = "#E0E0E0"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = top_genes, aes(label = gene),
    size = 3.5, color = "black", fontface = "bold",
    max.overlaps = 15, box.padding = 0.4
  ) +
  labs(
    title = "Differentially Expressed Genes: Hepatocyte Cells (Metastasis vs PT)",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "Significance"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_volcano4)

ggsave(filename = fig_path4, plot = p_volcano4, width = 9, height = 7, dpi = 300)
cat("Saved pastel version to:", fig_path4, "\n")




log_path <- "D:/hcc_project/results/project_log.txt"
log_entry <- paste0(
  "[", Sys.time(), "] ",
  "FINAL Volcano plot confirmed: pastel colors (lavender/mint/grey), English labels. ",
  "File: 10_volcano_DE_hepatocyte_v4.png — this is the approved version for the report.\n"
)
cat(log_entry, file = log_path, append = TRUE)
cat("تم تسجيل النسخة النهائية في project_log.txt\n")









p_volcano5 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance_en)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_color_manual(values = c(
    "Up in Metastasis" = "#C9A0DC",   # لافندر
    "Up in PT"         = "#8FD9C4",   # مينتي
    "Not Significant"  = "#A8C8E8"    # أزرق فاتح باستيل بدل الرمادي
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = top_genes, aes(label = gene),
    size = 3.5, color = "black", fontface = "bold",
    max.overlaps = 15, box.padding = 0.4
  ) +
  labs(
    title = "Differentially Expressed Genes: Hepatocyte Cells (Metastasis vs PT)",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "Significance"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_volcano5)

fig_dir <- "D:/hcc_project/results/figures/"
fig_path5 <- file.path(fig_dir, "11_volcano_DE_hepatocyte_final.png")
ggsave(filename = fig_path5, plot = p_volcano5, width = 9, height = 7, dpi = 300)
cat("Saved final version to:", fig_path5, "\n")




log_path <- "D:/hcc_project/results/project_log.txt"
log_entry <- paste0(
  "[", Sys.time(), "] ",
  "FINAL APPROVED Volcano plot: pastel colors (lavender/mint/light-blue), English labels. ",
  "File: 11_volcano_DE_hepatocyte_final.png — approved for final report.\n"
)
cat(log_entry, file = log_path, append = TRUE)
cat("تم تسجيل الاعتماد النهائي\n")analysis_script.R



length(list.files("D:/hcc_project/results/figures/"))
length(list.files("D:/hcc_project/results/tables/"))