library(ggplot2)
library(patchwork)

issues <- read.table("issue_counts_lengths_by_seqclass.nonoverlap.txt", header=TRUE, sep="\t")

issues$SeqClass <- factor(issues$SeqClass, levels=c("PAR1", "XTR", "XDR", "AMPL", "CEN-DYZ3", "SAT", "DYZ19", "HET", "OTHER", "PAR2", "TEL", "UNASSIGNED"))
issues$IssueType <- factor(issues$IssueType, levels=c("ERRBASE", "ERRSTRUCT", "NGAP"))
issues$Sample <- factor(issues$Sample)
# summary(issues)

# shared fill colors across all plots
issue_colors <- c(
  "ERRBASE"  = rgb(204,153,255, maxColorValue=255),
  "ERRSTRUCT"= rgb(255,153,153, maxColorValue=255),
  "NGAP"     = rgb(0,0,0,       maxColorValue=255)
)

# set a common font theme for all plots
theme_set(theme_bw(base_size = 6, base_family = "Arial"))

# 1. Plot number of issues by sequence class, faceted by issue type
common_theme <- list(
  geom_boxplot(aes(fill=IssueType), outlier.size = 0, linewidth = 0.1, outlier.shape = NA),
  geom_point(aes(fill=IssueType), position = position_jitter(seed = 42, width = 0.2, height = 0), shape=21, stroke = 0.2, alpha = 0.5),
  labs(x="Sequence Class", y="Num of Issues"),
  theme_bw(),
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")
)

p_errbase <- ggplot(subset(issues, IssueType=="ERRBASE"), aes(x=SeqClass, y=NumIssues)) +
  common_theme +
  scale_y_continuous(breaks = scales::breaks_width(100)) +
  scale_fill_manual(values = issue_colors) +
  ggtitle("ERRBASE")

p_errstruct <- ggplot(subset(issues, IssueType=="ERRSTRUCT"), aes(x=SeqClass, y=NumIssues)) +
  common_theme +
  scale_y_continuous(breaks = scales::breaks_width(2)) +
  scale_fill_manual(values = issue_colors) +
  ggtitle("ERRSTRUCT")

p_ngap <- ggplot(subset(issues, IssueType=="NGAP"), aes(x=SeqClass, y=NumIssues)) +
  common_theme +
  scale_y_continuous(breaks = scales::breaks_width(2)) +
  scale_fill_manual(values = issue_colors) +
  ggtitle("NGAP")

p_errbase / p_errstruct / p_ngap

ggsave("num_issues_by_seqclass.png", width=8, height=12)

# 2. Plot number of samples with issues by sequence class, faceted by issue type
samples_with_issues <- aggregate(NumIssues ~ SeqClass + IssueType, data=issues,
                                 FUN=function(x) sum(x > 0))
names(samples_with_issues)[3] <- "NumSamplesWithIssues"

ggplot(samples_with_issues, aes(x=SeqClass, y=NumSamplesWithIssues, fill=IssueType)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=issue_colors) +
  scale_y_continuous(limits=c(0, 142), breaks=scales::breaks_width(10)) +
  labs(x="Sequence Class", y="Num of Samples with Issues") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave("num_samples_with_issues_by_seqclass.png", width=8, height=4)

# 3. Plot IssueLength by sequence class, scaled to ClassLength
issues$IssueLengthScaled <- 100 * issues$IssueLength / issues$ClassLength

ggplot(issues, aes(x=SeqClass, y=IssueLengthScaled, fill=IssueType)) +
  geom_boxplot(outlier.size = 0, linewidth = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(seed = 42), shape=21, stroke = 0.2, alpha = 0.2) +
  stat_summary(aes(colour=IssueType),
              fun.min=min, fun.max=max, fun=median,
              geom="errorbar", linewidth=0.3, width=0.4,
              position=position_dodge(width=0.8)) +
  stat_summary(aes(fill=IssueType, colour=IssueType),
              fun=median,
              geom="crossbar", linewidth=0.3, width=0.4, fatten=1,
              position=position_dodge(width=0.8)) +
  geom_boxplot(aes(fill=IssueType, colour=IssueType),
              outlier.shape = 21,
              outlier.size = 1, outlier.stroke = 0.2, outlier.alpha = 0.2,
              linewidth = 0.25) +
  scale_fill_manual(values=issue_colors) +
  scale_colour_manual(values=issue_colors) +
  scale_y_continuous(breaks = scales::breaks_width(10)) +
  labs(x="Sequence Class", y="% affected") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggsave("issue_length_scaled_by_seqclass.png", width=120, height=50, units="mm")
ggsave("issue_length_scaled_by_seqclass.pdf", width=120, height=50, units="mm")

# 3-1. Zoom in
ggplot(issues, aes(x=SeqClass, y=IssueLengthScaled, fill=IssueType)) +
  geom_boxplot(outlier.size = 0, linewidth = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(seed = 42), shape=21, stroke = 0.2, alpha = 0.2) +
  stat_summary(aes(colour=IssueType),
              fun.min=min, fun.max=max, fun=median,
              geom="errorbar", linewidth=0.3, width=0.4,
              position=position_dodge(width=0.8)) +
  stat_summary(aes(fill=IssueType, colour=IssueType),
              fun=median,
              geom="crossbar", linewidth=0.3, width=0.4, fatten=1,
              position=position_dodge(width=0.8)) +
  geom_boxplot(aes(fill=IssueType, colour="gray"),
              outlier.shape = 21,
              outlier.size = 1, outlier.stroke = 0.2, outlier.alpha = 0.2,
              linewidth = 0.25) +
  scale_fill_manual(values=issue_colors) +
  scale_colour_manual(values=issue_colors) +
  scale_y_continuous(breaks = scales::breaks_width(10)) +
  ylim(0, 10) +
  labs(x="Sequence Class", y="% affected") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggsave("issue_length_scaled_by_seqclass_zoomed.png", width=120, height=50, units="mm")
ggsave("issue_length_scaled_by_seqclass_zoomed.pdf", width=120, height=50, units="mm")

# 3-2. Plot IssueLength by sequence class, scaled to ClassLength, with no UNASSIGNED
issues_no_unassigned <- subset(issues, SeqClass != "UNASSIGNED")
ggplot(issues_no_unassigned, aes(x=SeqClass, y=IssueLengthScaled, fill=IssueType)) +
  geom_point(position = position_jitterdodge(seed = 42), shape=21, stroke = 0.2, alpha = 0.2) +
  stat_summary(aes(colour=IssueType),
               fun.min=min, fun.max=max, fun=median,
               geom="errorbar", linewidth=0.3, width=0.4,
               position=position_dodge(width=0.8)) +
  stat_summary(aes(fill=IssueType, colour=IssueType),
               fun=median,
               geom="crossbar", linewidth=0.3, width=0.4, fatten=1,
               position=position_dodge(width=0.8)) +
  # geom_boxplot(aes(fill=IssueType),
  #               outlier.shape = 21,
  #               outlier.size = 1, outlier.stroke = 0.2, outlier.alpha = 0.7,
  #               linewidth = 0.25) +
  scale_fill_manual(values=issue_colors) +
  scale_colour_manual(values=issue_colors) +
  scale_y_continuous(breaks = scales::breaks_width(2)) +
  labs(x="Sequence Class", y="% Class affected per sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave("issue_length_scaled_by_seqclass_no_unassigned.png", width=120, height=80, units="mm")

# 4. Plot stacked bar plots for each sample how many bases are affected by each issue type
ggplot(issues, aes(x=Sample, y=IssueLength, fill=IssueType)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=issue_colors) +
  labs(x="Sample", y="Issues affected (bp)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
ggsave("issues_bps_by_sample.png", width=12, height=6)
