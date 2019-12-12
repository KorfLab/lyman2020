library(ggplot2)
library(ggsignif)
library(reshape2)
library(RColorBrewer)
library(cowplot)

setwd("~/Documents/Programming/python/lyman2020/")
selector <- read.delim("selector_out_2019_12_12.txt")
s2 <- selector[which(selector$txs == 1), c(3,5,7:10)]
s3 <- melt(s2, id.vars = c("len", "exons", "exp"),
           variable.name = "level", value.name = "isoforms")
s3$introns <- s3$exons - 1
s3$intronNormExp <- s3$exp / s3$introns
s3$lenNormExp <- s3$exp / (s3$len / 1000)
s3$doubleNormExp <- s3$lenNormExp / s3$introns
expQuantile <- quantile(s3$exp)
s3$expFactor <- ifelse(s3$exp <= expQuantile[2], "low",
                       ifelse(s3$exp <= expQuantile[4], "medium", "high"))
s3$expFactor <- factor(s3$expFactor, levels = c("low", "medium", "high"))
# drop intron counts containing too few examples to be useful
s4 <- s3[which(s3$introns < 8),]
# drop intron counts == 1, since these genes will never have more
# than one estimated isoform according to Alex's method.
s5 <- s4[which(s4$introns > 1),]
# set more informative names for iso4, iso6 and iso8
s5$sensitivity <- ifelse(s5$level == "iso4", "common",
                         ifelse(s5$level == "iso6", "common + uncommon",
                                "common + uncommon + rare"))
s5$sensitivity <- factor(s5$sensitivity, levels = c("common",
                                                    "common + uncommon",
                                                    "common + uncommon + rare"))
# for box plots, use only most conservative isoform estimates
s6 <- s5[which(s5$sensitivity == "common"),]

# plot isoform estimates vs coverage
scatter_panel <-  ggplot(data = s5, aes(x = doubleNormExp, y = isoforms)) +
  geom_point() +
  geom_smooth(aes(color = sensitivity)) +
  facet_wrap(~introns, ncol = 6) +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values = c("#de2d26", "#ffbf00", "#32936f")) +
  labs(x = "coverage", y = "predicted isoforms") +
  theme_bw() +
  theme(legend.position = "none")

# with comparisons by ggsignif library (wilcox.test method)
logBox_signif_facet <- ggplot(data = s6, aes(x = expFactor,
                                       y = isoforms,
                                       fill = expFactor)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("medium", "high"))) +
  facet_wrap(~as.factor(introns), ncol = 6) +
  scale_y_log10() +
  labs(x = "expression", y = "predicted isoforms") +
  scale_fill_brewer(palette = "Reds") +
  theme_bw() +
  theme(legend.position = "none")

# with manual comparisons using ggsignif
s6med <- s6[which(s6$expFactor == "medium"), c("introns","isoforms")]
s6high <- s6[which(s6$expFactor == "high"), c("introns","isoforms")]

pvals <- c(wilcox.test(s6med[s6med$introns == 2, "isoforms"],
                       s6high[s6high$introns == 2, "isoforms"])$p.value,
           wilcox.test(s6med[s6med$introns == 3, "isoforms"],
                       s6high[s6high$introns == 3, "isoforms"])$p.value,
           wilcox.test(s6med[s6med$introns == 4, "isoforms"],
                       s6high[s6high$introns == 4, "isoforms"])$p.value,
           wilcox.test(s6med[s6med$introns == 5, "isoforms"],
                       s6high[s6high$introns == 5, "isoforms"])$p.value,
           wilcox.test(s6med[s6med$introns == 6, "isoforms"],
                       s6high[s6high$introns == 6, "isoforms"])$p.value,
           wilcox.test(s6med[s6med$introns == 7, "isoforms"],
                       s6high[s6high$introns == 7, "isoforms"])$p.value)

logBox_manual_single <- ggplot(data = s6,
                               aes(x = as.factor(introns),
                                   y = isoforms,
                                   fill = expFactor)) +
  geom_boxplot() +
  geom_signif(annotation = formatC(pvals[1], digits = 1),
              y_position = 1.78, xmin = 1, xmax = 1.25,
              tip_length = c(0.02, 0.38)) +
  geom_signif(annotation = formatC(pvals[2], digits = 1),
              y_position = 1.52, xmin = 2, xmax = 2.25,
              tip_length = c(0.02, 0.15)) +
  geom_signif(annotation = formatC(pvals[3], digits = 1),
              y_position = 1.85, xmin = 3, xmax = 3.25,
              tip_length = c(0.02, 0.14)) +
  geom_signif(annotation = formatC(pvals[4], digits = 1),
              y_position = 2.59, xmin = 4, xmax = 4.25,
              tip_length = c(0.02, 0.35)) +
  geom_signif(annotation = formatC(pvals[5], digits = 1),
              y_position = 2.12, xmin = 5, xmax = 5.25,
              tip_length = c(0.02, 0.15)) +
  geom_signif(annotation = formatC(pvals[6], digits = 1),
              y_position = 2.1, xmin = 6, xmax = 6.25,
              tip_length = c(0.02, 0.16)) +
  scale_y_log10() +
  labs(y = "predicted isoforms") +
  scale_fill_brewer(palette = "Reds") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

# exon count chart
bar_single <- ggplot(data = s6, aes(x = as.factor(introns))) +
  geom_bar(fill = "gray80") +
  labs(x = "annotated introns", y = "genes") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

bar_facet <- ggplot(data = s6, aes(x = expFactor)) +
  geom_bar(fill = "gray80") +
  facet_wrap(~introns, ncol = 6) +
  labs(x = "expression", y = "genes") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

# only scatterplot faceted
pdf("scatter_box_bar_3panel_vertical_v1_9x12.pdf", height = 9, width = 12)
plot_grid(scatter_panel, logBox_manual_single, bar_single,
          labels = c("a", "b", "c"),
          ncol = 1,
          rel_heights = c(2.5, 2.5, 1),
          align = "v",
          axis = "l")
dev.off()

# all panels faceted
# plot_grid(scatter_panel,
#           logBox_signif_facet + theme(axis.title.x = element_blank()),
#           bar_facet,
#           labels = c("a", "b", "c"),
#           ncol = 1,
#           rel_heights = c(2.5, 2.5, 1),
#           align = "v",
#           axis = "l")

# all panels faceted, fewer labels
pdf("scatter_box_bar_3panel_vertical_v2_9x12.pdf", height = 9, width = 12)
plot_grid(scatter_panel,
          logBox_signif_facet +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  strip.background.x = element_blank(),
                  strip.text.x = element_blank()),
          bar_facet +
            theme(strip.background.x = element_blank(),
                  strip.text.x = element_blank()),
          labels = c("a", "b", "c"),
          ncol = 1,
          rel_heights = c(3, 2.5, 1),
          align = "v",
          axis = "l")
dev.off()

# not faceted by intron count
scatter_single <- ggplot(data = s5, aes(x = doubleNormExp, y = isoforms)) +
  geom_point() +
  geom_vline(xintercept = expQuantile[2], color = "gray60") +
  geom_vline(xintercept = expQuantile[4], color = "gray60") +
  geom_smooth(aes(color = sensitivity)) +
  annotate("text", x = 100, y = 10000000, label = "low", color = "gray60") +
  annotate("text", x = 20000, y = 10000000, label = "medium", color = "gray60") +
  annotate("text", x = 5000000, y = 10000000, label = "high", color = "gray60") +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values = c("#de2d26", "#ffbf00", "#32936f")) +
  labs(x = "coverage", y = "predicted isoforms") +
  theme_bw() +
  theme(legend.position = "none")

s5$expIso <- factor(paste(s5$sensitivity, s5$expFactor),
                    c("common high",
                      "common medium",
                      "common low",
                      "common + uncommon high",
                      "common + uncommon medium",
                      "common + uncommon low",
                      "common + uncommon + rare high",
                      "common + uncommon + rare medium",
                      "common + uncommon + rare low"))


i4 <- s5[s5$level == "iso4",]
i6 <- s5[s5$level == "iso6",]
i8 <- s5[s5$level == "iso8",]

annoFrame <- data.frame(sensitivity = c("common",
                                        "common + uncommon",
                                        "common + uncommon + rare"),
                        start = rep("medium", 3),
                        end = rep("high", 3),
                        y = c(3, 6.2, 8.9),
                        label = c(wilcox.test(i4[i4$expFactor == "medium", "isoforms"],
                                              i4[i4$expFactor == "high", "isoforms"])$p.value,
                                  wilcox.test(i6[i6$expFactor == "medium", "isoforms"],
                                              i6[i6$expFactor == "high", "isoforms"])$p.value,
                                  wilcox.test(i8[i8$expFactor == "medium", "isoforms"],
                                              i8[i8$expFactor == "high", "isoforms"])$p.value))

box_s5 <- ggplot(data = s5, aes(x = expFactor, y = isoforms)) +
  geom_boxplot(fill = c("#fee0d2", "#fc9272", "#de2d26",
                        "#fff4d4", "#ffe79e", "#ffbf00",
                        "#daf7ec", "#69cfa9", "#32936f")) +
  geom_signif(data = annoFrame, aes(xmin=start,
                                    xmax=end,
                                    y_position=y,
                                    annotations=formatC(label, digits = 1)),
              manual = TRUE)+
  facet_wrap(~sensitivity) +
  scale_y_log10() +
  labs(x = "expression") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank())

pdf("scatter_box_2panel_horizontal_6x12.pdf", height = 6, width = 12)
plot_grid(scatter_single, box_s5, labels = c("a", "b"),
          align = "h", axis = "t")
dev.off()

