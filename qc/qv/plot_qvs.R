#!/bin/Rscript

library(ggplot2)
library(data.table)
library("scales")

# Read input data
# Sample AsmGroup QV NumSequences TotalSize SeqGroup
# data <- fread("qv_plot_input.tsv", header=T, sep="\t")
data <- fread("qv_plot_input_fix_seqgroup.tsv", header=T, sep="\t")
data$AsmGroup <- factor(data$AsmGroup, levels = c("PriorBest", "T2T-Y"))
df_sorted <- data[order(data$Sample, data$AsmGroup), ]

# Jitter position for better visualization
jitter_pos <- position_jitter(seed = 42, width = 0.2, height = 0)

# Plot violin plots for QV on the Y axis and AsmGroup on the X axis; points colored by SeqGroup and line connecting same Samples
p <- ggplot(df_sorted, aes(x=AsmGroup, y = QV)) +
  theme_bw() +
  geom_line(aes(x=AsmGroup, y = QV, group = Sample), position = jitter_pos, linewidth = 0.2, color = "darkgray") +
  geom_point(aes(color = SeqGroup), position = jitter_pos, shape=21, stroke = 0.2, alpha = 0.8) +
  geom_violin(alpha=0.2, linewidth=0.2) +
  geom_boxplot(width=0.1, position=position_dodge(0.9), outlier.shape = NA) +
  scale_y_continuous(labels=comma) +
  labs(
     x="Assembly Group",
     y="QV (Phred)")

ggsave("qv_violin_plot.png", plot=p, width=3, height=3)
ggsave("qv_violin_plot.pdf", plot=p, width=3, height=3)

# Plot violin plots for num. sequences on the Y axis and AsmGroup on the X axis; points colored by SeqGroup and line connecting same Samples
# Invert the Y axis for better visualization
p <- ggplot(df_sorted, aes(x=AsmGroup, y = NumSequences)) +
  theme_bw() +
  geom_line(aes(x=AsmGroup, y = NumSequences, group = Sample), position = jitter_pos, linewidth = 0.2, color = "darkgray") +
  geom_point(aes(color = SeqGroup), position = jitter_pos, shape=21, stroke = 0.2, alpha = 0.8) +
  geom_violin(alpha=0.2, linewidth=0.2) +
  geom_boxplot(width=0.1, position=position_dodge(0.9), outlier.shape = NA) +
  scale_y_continuous(labels=comma, trans = "log10", breaks = scales::log_breaks(n = 12)) +
  labs(
     x="Assembly Group",
     y="Num. of sequences")
ggsave("numseq_violin_plot.png", plot=p, width=3, height=3)
ggsave("numseq_violin_plot.pdf", plot=p, width=3, height=3)

# Plot QV comparison across Y assemblies, colored by SeqGroup and sized by TotalSize
p <- ggplot(df_sorted, aes(x=AsmGroup, y = QV, group = Sample, color = SeqGroup)) +
  theme_bw() +
  geom_line(position = jitter_pos, linewidth = 0.3, color = "darkgray") +
  geom_point(position = jitter_pos, shape=21, aes(size = TotalSize), stroke = 1) +
  scale_size_continuous(labels=comma, range = c(1, 10), name = "Assembly Size") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") + 
  labs(title="QV comparison across Y assemblies",
	   x="Assembly Group",
	   y="QV (Phred)")

ggsave("qv_plot_by_asmsize.png", plot=p, width=8, height=6)
ggsave("qv_plot_by_asmsize.pdf", plot=p, width=8, height=6)

# Plot QV comparison across Y assemblies, colored by SeqGroup and sized by 1/NumSequences
p <- ggplot(df_sorted, aes(x=AsmGroup, y = QV, group = Sample, color = SeqGroup)) +
  theme_bw() +
  geom_line(position = jitter_pos, linewidth = 0.3, color = "darkgray") +
  geom_point(position = jitter_pos, shape=21, aes(size = 1/NumSequences), stroke = 1) +
  scale_size_continuous(labels=comma, range = c(1, 10), name = "1/Num.ofSequences") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") + 
  labs(title="QV comparison across Y assemblies",
	   x="Assembly Group",
	   y="QV (Phred)")

ggsave("qv_plot_by_numseq.png", plot=p, width=8, height=6)
ggsave("qv_plot_by_numseq.pdf", plot=p, width=8, height=6)

# Plot NumSequences across Y assemblies, colored by SeqGroup and sized by QV
p <- ggplot(df_sorted, aes(x=AsmGroup, y = NumSequences, group = Sample, color = SeqGroup)) +
  theme_bw() +
  geom_line(position = jitter_pos, linewidth = 0.3, color = "darkgray") +
  geom_point(position = jitter_pos, shape=21, aes(size = QV), stroke = 1) +
  scale_size_continuous(labels=comma, range = c(1, 10), name = "QV (Phred)") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") + 
  labs(title="Num. of sequences across Y assemblies",
     x="Assembly Group",
     y="Num. of sequences")

ggsave("numseq_plot_by_qv.png", plot=p, width=8, height=6)
ggsave("numseq_plot_by_qv.pdf", plot=p, width=8, height=6)

# summary(data, digits = 3)
# Summary by group
summary(data[data$AsmGroup=="PriorBest", ], digits = 3)
summary(data[data$AsmGroup=="T2T-Y", ], digits = 3)

# Plot PriorBest vs T2T-Y QV comparison
data= fread("qv_plot_by_asmset_input_fixseqgroup.tsv", header=T, sep="\t")
p <- ggplot(data, aes(x=PriorBestYQV, y = CuratedYQV, color = SeqGroup)) +
  theme_bw() +
  geom_abline(slope=1, intercept=0, linetype="dashed", color = "darkgray") +
  geom_point(shape=21, stroke = 1, alpha = 0.8, aes(size = CuratedYSize)) +
  scale_x_continuous(labels=comma, limits = c(20, 65)) +
  scale_y_continuous(labels=comma, limits = c(20, 65)) +
  scale_size_continuous(labels=comma, range = c(3, 7), name = "T2T-Y Size") +
  # scale_color_brewer(palette = "Set1") +
  # scale_fill_brewer(palette = "Set1") + 
  labs(title="QV comparison: PriorBestY vs T2T-Y",
     x="Prior Best Y QV (Phred)",
     y="T2T-Y QV (Phred)")

ggsave("qv_plot_by_asmgroup.png", plot=p, width=8, height=6)


