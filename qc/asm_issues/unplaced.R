library(ggplot2)

# Use a gray color pallette for NumIssues
colors <- c(
  "0" = rgb(220,220,220, maxColorValue=255),
  "1" = rgb(150,150,150, maxColorValue=255),
  "2" = rgb(100,100,100, maxColorValue=255),
  "3" = rgb(50,50,50, maxColorValue=255),
  "4" = rgb(20,20,20, maxColorValue=255),
  ">4" = rgb(0,0,0, maxColorValue=255)
)

issues <- read.table("issue_counts_lengths_by_seqclass.nonoverlap.gap.txt", header=TRUE, sep="\t")

issues$SeqClass <- factor(issues$SeqClass, levels=c("PAR1",
  "XTR", "XDR", "AMPL", "CEN-DYZ3", "SAT", "DYZ19", "HET",
  "OTHER", "PAR2", "TEL", "UNASSIGNED"))

issues$IssueType <- factor(issues$IssueType, levels=c("NGAP"))
issues$Sample <- factor(issues$Sample)

# Reverse order of NumIssues for better stacking in the plot
issues$NumIssues <- factor(issues$NumIssues, levels=c("0", "1", "2", "3", "4", ">4"))

# Exclude samples with 0 NGAPs for the plot
# issues <- subset(issues, NumIssues != "0")
# --> include "0"

# Plot a stacked bar plot by number of NGAP issues by sequence class
ggplot(issues, aes(x=SeqClass, fill=NumIssues)) +
  geom_bar() +
  labs(x="sequence class", y="num. of samples with gaps") +
  scale_y_continuous(breaks = scales::breaks_width(5)) +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(values=colors, name="num. of gaps") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ggtitle("gaps by sequence classes in T2T-Y")

ggsave("ngap_issues_by_seqclass.png", width=8, height=4)

issues_prior <- read.table("gaps_in_prior_wi_0.txt", header=TRUE, sep="\t")

issues_prior$Sample <- factor(issues_prior$Sample)
issues_prior$SeqClass <- factor(issues_prior$SeqClass, levels=c("PAR1",
  "XTR", "XDR", "AMPL", "CEN-DYZ3", "SAT", "DYZ19", "HET",
  "OTHER", "PAR2", "TEL", "UNASSIGNED"))
# Reverse order of NumIssues for better stacking in the plot
issues_prior$NumIssues <- factor(issues_prior$NumIssues, levels=c("0", "1", "2", "3", "4", ">4"))

ggplot(issues_prior, aes(x=SeqClass, fill=NumIssues)) +
  geom_bar() +
  labs(x="sequence class", y="num. of samples with gaps") +
  scale_y_continuous(breaks = scales::breaks_width(5)) +
  scale_fill_manual(values=colors, name="num. of gaps") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ggtitle("gaps by sequence classes in prior Ys")

ggsave("ngap_issues_by_seqclass.prior.png", width=8, height=4)

# Compare the two plots side by side for each sequence class
common_cols <- c("Sample", "SeqClass", "NumIssues")
dat <- rbind(
  data.frame(issues[, common_cols], Label="T2T-Y"),
  data.frame(issues_prior[, common_cols], Label="Prior")
)

# Ensure all SeqClass x Label combinations exist so empty bars are shown
all_combinations <- expand.grid(
  SeqClass = levels(issues$SeqClass),
  Label = c("Prior", "T2T-Y"),
  stringsAsFactors = FALSE
)
dat <- merge(all_combinations, dat, by=c("SeqClass", "Label"), all.x=TRUE)
dat$SeqClass <- factor(dat$SeqClass, levels=levels(issues$SeqClass))
dat$Label <- factor(dat$Label, levels=c("Prior", "T2T-Y"))
dat$NumIssues <- factor(dat$NumIssues, levels=c("0", "1", "2", "3", "4", ">4"))
# keep NA rows so empty bars are reserved in facets; geom_bar will skip them

ggplot(dat, aes(x=Label, fill=NumIssues)) + 
    geom_blank() +
    geom_bar(data=subset(dat, !is.na(NumIssues))) +
    facet_grid(~ SeqClass, scales="free_x", space="free_x", switch="x") +
    labs(x="Sequence Class", y="Num. of samples with gaps") +
    scale_fill_manual(values=colors, name="Num. of gaps") +
    scale_x_discrete(drop=FALSE) +
    scale_y_continuous(breaks = scales::breaks_width(5)) +
    theme_bw() +
    theme(
      axis.title = element_text(size=6),
      axis.text = element_text(size=5),
      axis.text.x = element_text(angle=45, hjust=1),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(size=5),
      panel.spacing = unit(0.1, "lines"),
      legend.key.height = unit(0.3, "cm"),
      legend.key.width  = unit(0.3, "cm"),
      legend.text = element_text(size=5),
      legend.title = element_text(size=5)
    )

ggsave("ngap_issues_by_seqclass.prior_vs_new.png", width=6, height=3)
ggsave("ngap_issues_by_seqclass.prior_vs_new.pdf", width=100, height=60, units="mm")