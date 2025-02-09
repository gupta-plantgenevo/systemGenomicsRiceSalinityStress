library(ggplot2)
library(forcats)
#plotdat = read.table("module_pathway_diff/GO_BP/Golev3_forPlot", header = T, sep="\t")
plotdat = read.table("ricesalinity/data/module_pathway_diff/GO_BP/Sraw/plotDat_golev4_med0.95.txt", header = T, sep="\t")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#col = gg_color_hue(2)
col = c('lightpink2','lightblue3')

f2_1 <- ggplot(plotdat, aes(y = fct_inorder(GO.Name), x = medS, color=Env)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = medS-ci, xmax = medS+ci), height = 0.25) +
  scale_color_manual(values = col) + 
  xlab("Selection Strength") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"), 
        legend.position = "none")
f2_1


w = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/swar_pc_plot.txt", header = T, sep = "\t")
tmp = w[,-c(3,4)]

val = c()
val = append(val, as.character(tmp[1,-c(1,6)]))
val = append(val, as.character(tmp[2,-c(1,6)]))
val = append(val, as.character(tmp[3,-c(1,6)]))
val = append(val, as.character(tmp[4,-c(1,6)]))
val=as.numeric(val)
t = data.frame(env=c("PCW","PCW","PCW","PCW","PCW","PCW","PCW","PCW","PCW","PCW","PCW","PCW","PCS","PCS","PCS","PCS"),
               PC=c("PC1","PC1","PC1","PC1","PC2","PC2","PC2","PC2","PC4","PC4","PC4","PC4","PC5","PC5","PC5","PC5"),
               ID=c("beta","dir","indir","tot","beta","dir","indir","tot","beta","dir","indir","tot","beta","dir","indir","tot"),
               Value=val)

PClabel = as_labeller(c(PC1="PC1[normal]", PC2="PC2[normal]", PC4="PC3[normal]", PC5="PC2[salt]"), 
                      default = label_parsed)
#col = c('#FFCB05','#00274C')
col = c('lightpink2','lightblue3')
f2_2 <- ggplot(t,aes(x = ID, y = Value, fill = env)) + 
  geom_bar(position="stack",  stat = "identity") + 
  scale_y_continuous(breaks = seq(-0.006,0.006,0.003),limits = c(-0.006,0.004)) + 
  scale_fill_manual(values = col) + 
  facet_grid(~PC,labeller = PClabel) +
  ylab("Selection") +
  theme_bw() + 
  scale_x_discrete(labels = c(expression(~beta), "D", "I", "T", expression(~beta), "D", "I", "T", expression(~beta), "D", "I", "T")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        strip.text.x = element_text(size = 14, color = "black"),
        axis.line.y = element_line(color = 'black'),
        axis.line.x = element_blank(),
        axis.title.x =element_blank(),
        axis.text = element_text(size=16, color="gray10"),
        axis.title =element_text(size=24, color="black"),
        axis.ticks.length.y = unit(0.3, "cm"), 
        axis.ticks.x = element_blank(), legend.position="none" )

library(cowplot)
fig2 = plot_grid(f2_1,
  plot_grid(NULL, f2_2, NULL,
            ncol = 1, nrow = 3, rel_heights = c(1,2,1)),
  nrow = 1, ncol = 2, labels = c('a.','b.'),
  label_size = 18, label_fontface = "bold")



