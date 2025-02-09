##Plot
#From the selection results I just saved the plotting points in the format below (gr is group -- wet beta, wet gamma etc)
#Pheno	Env	SelType	gr	Selcoeff	SE
#LOP	Wet	beta	wb	0.062366308	0.03374911
phenoPlot = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/phenoSel_plotData.txt", header = T, sep = "\t")
phenoPlot$Pheno = as.character(c(1,1,1,1,2,2,2,2,3,3,3,3))
phenoPlot$gr = as.character(c(1,2,3,4,1,2,3,4,1,2,3,4))

Phenolab = c("LOP", "Chl_a", "FT")
names(Phenolab) <- c("1","2","3")
library(ggplot2)

#colors = c('#00274C',"lightsteelblue4",'#FFCB05', "#CFC096")
colors = c("lightblue3","lightsteelblue4","lightpink2","mistyrose3")
ggplot(phenoPlot, aes(x = gr, y = Selcoeff, fill = gr)) + 
  geom_bar(position="stack",  stat = "identity") +
  geom_errorbar(aes(ymin=Selcoeff-SE, ymax=Selcoeff+SE), width=.2,
                position=position_dodge(0.05)) + 
  geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colors) + 
  facet_grid(~Pheno,labeller = labeller(Pheno=Phenolab)) +
  ylab("Selection Gradient") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        strip.text.x = element_text(size = 16, color = "black"),
        axis.line.y = element_line(color = 'black'),
        axis.line.x = element_blank(),
        axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y =element_text(size=20, color="black"), 
        axis.text.y = element_text(size=16, color="gray10"), 
        axis.ticks.length.y = unit(0.3, "cm"), 
        axis.ticks.x = element_blank())

print(g)

