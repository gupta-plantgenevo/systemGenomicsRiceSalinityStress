library(ggplot2)
library(dplyr)

plotData = read.table("ricesalinity/data/Fitness/Fitness_byGeno_Env.txt",  header = T, sep="\t")
tmp1 = plotData[which(plotData$Fecundity_Salt < plotData$Fecundity_Wet),]
tmp2 = plotData[which(plotData$Fecundity_Salt > plotData$Fecundity_Wet),]
a1 = ggplot() + geom_point(aes(x = tmp1$Fecundity_Wet, y = tmp1$Fecundity_Salt), size = 3, color = "lightblue3") + 
  geom_point(aes(x = tmp2$Fecundity_Wet, y = tmp2$Fecundity_Salt), size = 3, color = "lightpink2") + 
  geom_abline(intercept = 0, slope = 1, color="gray30", linetype = 3) + 
  xlab("Fecundity in Normal field") + ylab("Fecundity in Saline field") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black', size = 0.65),
        axis.title = element_text(size=20, color="black"),
        axis.text = element_text(size=16), 
        axis.ticks.length = unit(0.3, "cm"))

varH2_an = read.table("ricesalinity/data/QuanGen/mixedANOVA/INDtrans_SS_MS_F_P_H2_wo0.5.txt", header = T, sep="\t")
varH2_an = varH2_an[which(varH2_an$fdr_g < 0.001),]
a2 = varH2_an %>%
  ggplot(aes(x=H2)) +
  geom_histogram(binwidth=0.05, fill="gray70", color="#e9ecef", alpha=0.9) +
  theme_bw() + scale_x_continuous(breaks = seq(0, 1, by=0.25), expand=c(0,0)) + 
  scale_y_continuous(breaks = seq(0, 1600, by=500), expand=c(0,0)) +
  xlab(bquote('Heritability '~(H^2))) + ylab("Count") + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.text = element_text(size=16),
        axis.title = element_text(size=20, color="black"),
        axis.line = element_line(color = 'black', linewidth = 0.65),
        axis.ticks.length=unit(0.3, "cm"))


sel_coeff_wet = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_wet.w0fecun.nonStd.txt", header=TRUE, sep = "\t")
sel_coeff_salt = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_salt.w0fecun.nonStd.txt", header=TRUE, sep = "\t")
a3 = ggplot() + 
  stat_qq(aes(sample = abs(sel_coeff_wet$S)), colour = "lightblue3") + 
  stat_qq(aes(sample = abs(sel_coeff_salt$S)), colour = "lightpink2") + 
  theme_bw() + ylab(expression(Selection~differential~italic("|S|"))) + xlab("Theoretical Quantiles") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_text(size=20, color="black"), 
        axis.text = element_text(size=16), 
        axis.line = element_line(color = 'black', size = 0.65),
        axis.ticks.length=unit(0.3, "cm"))


a4 = ggplot() + 
  stat_qq(aes(sample = (sel_coeff_wet$S)), colour = "lightblue3") + 
  stat_qq(aes(sample = (sel_coeff_salt$S)), colour = "lightpink2") + 
  theme_bw() + ylab(expression(Selection~differential~italic("S"))) + xlab("Theoretical Quantiles") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_text(size=20, color="black"), 
        axis.text = element_text(size=16), 
        axis.line = element_line(color = 'black', size = 0.65),
        axis.ticks.length=unit(0.3, "cm"))

a5 = ggplot() + 
  stat_qq(aes(sample = (sel_coeff_wet$C)), colour = "lightblue3") + 
  stat_qq(aes(sample = (sel_coeff_salt$C)), colour = "lightpink2") + 
  theme_bw() + ylab(expression(Selection~differential~italic("C"))) + xlab("Theoretical Quantiles") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_text(size=20, color="black"), 
        axis.text = element_text(size=16), 
        axis.line = element_line(color = 'black', size = 0.65),
        axis.ticks.length=unit(0.3, "cm"))


alldat = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/CNAP_all.infor.txt", header = T, sep="\t")
unique(alldat$Condition)
cn = alldat[which(alldat$Condition == "Conditionally Neutral"),]
samedit_sel = alldat[which(alldat$Condition == "Detrimental"| alldat$Condition == "Beneficial"),]
ap_wetben = alldat[which(alldat$Condition == "Antagonistically Pleiotropic (Wet Beneficial)"),]
ap_saltben = alldat[which(alldat$Condition == "Antagonistically Pleiotropic (Salt Beneficial)"),]
a6 = ggplot() + theme_bw() +
  geom_point(data=cn, aes(x=S_salt, y=S_wet), color='gray85', size=2.5) +  
  geom_point(data=samedit_sel, aes(x=S_salt, y=S_wet), color='black', size=2.5) + 
  geom_point(data=ap_wetben, aes(x=S_salt, y=S_wet), color='royalblue', size=2.5) +  
  geom_point(data=ap_saltben, aes(x=S_salt, y=S_wet), color='magenta3', size=2.5) +
  xlim(-0.4,0.5) + ylim(-0.4,0.5) +
  labs(x = expression(~italic(S)~ '(Saline field)'), y=expression(~italic(S)~ '(Normal field)')) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_text(size=20, color="black"), 
        axis.text = element_text(size=16), 
        axis.line = element_line(color = 'black', size = 0.65),
        axis.ticks.length=unit(0.3, "cm"))


library(cowplot)
fig1 = plot_grid(a1, a2, a3, a4, a5, a6, 
          labels = c('a.', 'b.', 'c.', 'd.', 'e.', 'f.'), 
          nrow = 2, ncol = 3, label_size = 18,
          label_fontface = "bold")

