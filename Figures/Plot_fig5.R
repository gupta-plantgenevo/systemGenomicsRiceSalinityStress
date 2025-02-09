setwd("ricesalinity/data/eQTL/ind_data/fdr_trans1Mb/")

prop_wetsame = read.table("wet_common_samesign_cistrans.txt", header = T, sep="\t")
prop_wetdiff= read.table("wet_common_diffsign_cistrans.txt", header = T, sep="\t")
prop_saltsame = read.table("salt_common_samesign_cistrans.txt", header = T, sep="\t")
prop_saltdiff = read.table("salt_common_diffsign_cistrans.txt", header = T, sep="\t")

wetcomm = rbind(wetsamesign, wetdiffsign)
saltcomm = rbind(saltsamesign, saltdiffsign)

colors = c("#FF8A65", "#4DB6AC")
#colors = c('#FFCB05','#00274C')
p1 = ggplot(data=wetcomm, aes(beta.x, beta.y, colour = grp))+ 
  geom_point(alpha = 0.7, size = 2) +
  scale_colour_manual(values=colors) + 
  scale_y_continuous(breaks = c(-2.5, 0, 2.5)) + 
  labs(x=expression(~italic(cis)~"-eQTL effect size"), y=expression("mean "~italic(trans)~"-eQTL effect size")) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'),
        axis.title.x =element_text(size=16, color="black"),
        axis.text.x = element_text(size=12, color="gray10"),
        axis.title.y =element_text(size=16, color="black"),
        axis.text.y = element_text(size=12, color="gray10"), 
        axis.ticks.length.y = unit(0.3, "cm"), 
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.position = "none")

p2 = ggplot(data=saltcomm, aes(beta.x, beta.y, colour = grp))+ 
  geom_point(alpha = 0.7, size = 2) +
  scale_colour_manual(values=colors) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) + 
  labs(x=expression(~italic(cis)~"-eQTL effect size"), y=expression("mean "~italic(trans)~"-eQTL effect size")) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'),
        axis.title.x =element_text(size=16, color="black"),
        axis.text.x = element_text(size=12, color="gray10"),
        axis.title.y =element_text(size=16, color="black"),
        axis.text.y = element_text(size=12, color="gray10"), 
        axis.ticks.length.y = unit(0.3, "cm"), 
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.position = "none")

wetcov = read.table("ricesalinity/data/QuanGen/Ind_Wet_covariates.all.txt", header = T, sep="\t")
saltcov = read.table("ricesalinity/data/QuanGen/Ind_Salt_covariates.all.txt", header = T, sep="\t")

z = as.data.frame(wetcov[which(wetcov$geneID %in% wetsamesign$gene),]$ExpPolymorphismWet)
names(z)[1] = "ExpPol"
z$grp = rep("Reinforcing", times=nrow(z))
y = as.data.frame(wetcov[which(wetcov$geneID %in% wetdiffsign$gene),]$ExpPolymorphismWet)
names(y)[1] = "ExpPol"
y$grp = rep("Compensating", times=nrow(y))
x = rbind(z,y)
wilcox.test(z$ExpPol, y$ExpPol, alternative = "greater")

library(ggpubr)
colors = c("#4DB6AC", "#FF8A65")

b1 = ggplot(data=x, aes(grp, ExpPol, fill = grp))+ 
  geom_boxplot(alpha = 0.7) +
  scale_colour_manual(values=colors) + 
  labs(x="", y ="Inter-varietal Expression Variation") +
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"), 
                     label = "p.format",size = 5, label.x = 1.4) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'),
        axis.title.x =element_text(size=16, color="black"),
        axis.text.x = element_text(size=12, color="gray10"),
        axis.title.y =element_text(size=16, color="black"),
        axis.text.y = element_text(size=12, color="gray10"), 
        axis.ticks.length.y = unit(0.3, "cm"), 
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.position = "none")



z = as.data.frame(saltcov[which(saltcov$geneID %in% saltsamesign$gene),]$ExpPolymorphismSalt)
names(z)[1] = "ExpPol"
z$grp = rep("Reinforcing", times=nrow(z))
y = as.data.frame(saltcov[which(saltcov$geneID %in% saltdiffsign$gene),]$ExpPolymorphismSalt)
names(y)[1] = "ExpPol"
y$grp = rep("Compensating", times=nrow(y))
x = rbind(z,y)
wilcox.test(z$ExpPol, y$ExpPol, alternative = "greater")


b2 = ggboxplot(x, x = "grp", y = "ExpPol", fill = "grp", palette = c('#00274C', '#FFCB05'), alpha=0.7) +  
  labs(x="", y ="Inter-varietal Expression Variation") +
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"), 
                     label = "p.format",size = 5, label.x = 1.4, label.y = 6) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'),
        axis.title.x =element_text(size=16, color="black"),
        axis.text.x = element_text(size=12, color="gray10"),
        axis.title.y =element_text(size=16, color="black"),
        axis.text.y = element_text(size=12, color="gray10"), 
        axis.ticks.length.y = unit(0.3, "cm"), 
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.position = "none")





setwd("ricesalinity/data/eQTL/ind_data/ems_paper_analysis/")
library(dplyr)
library(ggplot2)
maf = read.table("SNPid_maf.txt", header = F, sep="\t")
names(maf)[1] = "SNP_id"; names(maf)[2] = "maf"

cis.salt.eqtl = read.table("../fdr_trans1Mb/cis-eQTLs.salt.leadSNPs.txt", header = T, sep="\t")
cis.wet.eqtl = read.table("../fdr_trans1Mb/cis-eQTLs.wet.leadSNPs.txt", header = T, sep="\t")
trans.salt.eqtl = read.table("../fdr_trans1Mb/trans-eQTLs.salt.1mb.leadSNPs.txt", header = T, sep="\t")
trans.wet.eqtl = read.table("../fdr_trans1Mb/trans-eQTLs.wet.1mb.leadSNPs.txt", header = T, sep="\t")
names(cis.salt.eqtl)[5] = "SNP_id"; names(cis.wet.eqtl)[5] = "SNP_id"; 
names(trans.salt.eqtl)[5] = "SNP_id"; names(trans.wet.eqtl)[5] = "SNP_id"; 

cis.wet = maf[which(maf$SNP_id %in% cis.wet.eqtl$SNP_id),]
cis.salt = maf[which(maf$SNP_id %in% cis.salt.eqtl$SNP_id),]
trans.wet = maf[which(maf$SNP_id %in% trans.wet.eqtl$SNP_id),]
trans.salt = maf[which(maf$SNP_id %in% trans.salt.eqtl$SNP_id),]
mafleft_wet = data.frame(setdiff(maf, rbind(cis.wet, trans.wet)))
mafleft_salt = data.frame(setdiff(maf, rbind(cis.salt, trans.salt)))

t = as.data.frame(cis.wet[,2])
names(t)[1] = "maf"
t$cond = rep("Cis-eQTL", times = nrow(t))

t1 = as.data.frame(trans.wet[,2])
names(t1)[1] = "maf"
t1$cond = rep("Trans-eQTL", times = nrow(t1))

t2 = data.frame(mafleft_wet[,2])
names(t2)[1] = "maf"
t2$cond = rep("Background", times = nrow(t2))
test = rbind(t2, t, t1)

br = seq(0,0.5, by=0.05)
wet = test %>% group_by(cond) %>% mutate(maf_bin = cut(maf, breaks=br)) %>% 
  group_by(cond, maf_bin) %>% summarise(sum_bin = n()) %>% 
  group_by(cond) %>% mutate(percent = sum_bin/sum(sum_bin)) %>% 
  as.data.frame() 
wet[c(1,11,21),]$sum_bin = c(0,0,0)
wet[c(1,11,21),]$percent = c(0,0,0)


t = as.data.frame(cis.salt[,2])
names(t)[1] = "maf"
t$cond = rep("Cis-eQTL", times = nrow(t))

t1 = as.data.frame(trans.salt[,2])
names(t1)[1] = "maf"
t1$cond = rep("Trans-eQTL", times = nrow(t1))

t2 = as.data.frame(mafleft_salt[,2])
names(t2)[1] = "maf"
t2$cond = rep("Background", times = nrow(t2))
test = rbind(t2, t, t1)

br = seq(0.0,0.5, by=0.05)
salt = test %>% group_by(cond) %>% mutate(maf_bin = cut(maf, breaks=br)) %>% 
  group_by(cond, maf_bin) %>% summarise(sum_bin = n()) %>% 
  group_by(cond) %>% mutate(percent = sum_bin/sum(sum_bin)) %>% 
  as.data.frame() 
salt[c(1,11,21),]$sum_bin = c(0,0,0)
salt[c(1,11,21),]$percent = c(0,0,0)

colors = c("gray50", "#00308F", "#FFCB05")
g1 = ggplot(wet, aes(x=maf_bin, y=percent*100, fill=cond)) + 
  geom_col(colour="black",width=0.7,position=position_dodge(0.7), alpha=0.7) +
  scale_fill_manual(values=colors, name = "Condition") + 
  labs(x="MAF", y = "Percent SNPs") +
  scale_x_discrete(labels = c("0.05","0.1", "0.15","0.2","0.25","0.3","0.35","0.4","0.45","0.5")) +
  guides(fill = guide_legend(title = "")) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'),
        axis.title.x =element_text(size=16, color="black"),
        axis.text.x = element_text(size=12, color="gray10"),
        axis.title.y =element_text(size=16, color="black"),
        axis.text.y = element_text(size=12, color="gray10"), 
        axis.ticks.length.y = unit(0.3, "cm"), 
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.position = c(0.85, 0.85), 
        legend.title = element_text(size=12,face="bold"),
        legend.text = element_text(size=12,face="bold"))

g2 = ggplot(salt, aes(x=maf_bin, y=percent*100, fill=cond)) + 
  geom_col(colour="black",width=0.7,position=position_dodge(0.7), alpha=0.7) +
  scale_fill_manual(values=colors, name = "Condition") + 
  labs(x="MAF", y = "Percent SNPs") +
  scale_x_discrete(labels = c("0.05","0.1", "0.15","0.2","0.25","0.3","0.35","0.4","0.45","0.5")) +
  guides(fill = guide_legend(title = "")) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'),
        axis.title.x =element_text(size=16, color="black"),
        axis.text.x = element_text(size=12, color="gray10"),
        axis.title.y =element_text(size=16, color="black"),
        axis.text.y = element_text(size=12, color="gray10"), 
        axis.ticks.length.y = unit(0.3, "cm"), 
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.position = c(0.85, 0.85), 
        legend.title = element_text(size=12,face="bold"),
        legend.text = element_text(size=12,face="bold"))

bottom_row <- plot_grid(p2, b2, labels = c('a.', 'b.'), label_size = 14)
plot_grid(bottom_row, g2, labels = c('', 'c.'), label_size = 14, ncol = 1)

bottom_row <- plot_grid(p1, b1, labels = c('a.', 'b.'), label_size = 14)
plot_grid(bottom_row, g1, labels = c('c.', ''), label_size = 14, ncol = 1)
