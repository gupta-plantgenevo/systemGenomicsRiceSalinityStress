setwd("ricesalinity/data/module_pathway_diff/decoherence/ind_salt_wet/revisedSel/")
#This will analyse and plot the results
library(data.table)
library(ggplot2)
library(dplyr)
results = fread("ricesalinity/data/module_pathway_diff/decoherence/ind_salt_wet/deco_trwith_above01S_lessthan50percentNA.txt", header = TRUE, sep="\t")
#results = fread("deco_tr_withSmean_above05.sortedQ.txt", nrows=32000, header = TRUE, sep="\t")

res_sig = results[which(results$qvalue<0.05),]
dim(res_sig)
tmp = res_sig[,c(4,5,1,2,3,6)]
write.table(tmp, "ricesalinity/data/module_pathway_diff/decoherence/ind_salt_wet/decoSIG_trwith_above01S_lessthan50percentNA.txt", row.names = F, col.names = T, quote = F, sep = "\t" )

library(sqldf)
res_sig_highwet = res_sig[which(abs(res_sig$cor_Wet) > abs(res_sig$cor_Salt)),]
res_sig_highsalt = res_sig[which(abs(res_sig$cor_Wet) < abs(res_sig$cor_Salt)),]
remain = sqldf('SELECT * FROM results EXCEPT SELECT * FROM res_sig')
dim(res_sig_highwet)
dim(res_sig_highsalt)

#ggplot(tmp, aes(x=Year, y=abs(Corr), fill=TrType)) + geom_boxplot() + xlab("Correlation by Year") + ylab("Difference in the magnitude of correlation in Transcripts") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot() + 
  geom_point(data=remain, aes(x=cor_Wet, y=cor_Salt), colour="grey90") + 
  geom_point(data=res_sig_highwet, aes(x=cor_Wet, y=cor_Salt), colour="lightblue3") + 
  geom_point(data=res_sig_highsalt, aes(x=cor_Wet, y=cor_Salt), colour="lightpink2") + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 0.75) + 
  xlab("CorrelationNormal") + ylab("CorrelationSalt") + 
  xlim(-0.6,1) + ylim(-0.6,1) + 
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.50), 
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.ticks.length=unit(0.3, "cm"))


t=res_sig_highsalt$tr1; t=append(t, res_sig_highsalt$tr2); 
x=as.data.frame(table(t))
t=res_sig_highwet$tr1; t=append(t, res_sig_highwet$tr2); 
y=as.data.frame(table(t))
t=res_sig$tr1; t=append(t, res_sig$tr2); 
z=as.data.frame(table(t))
View(x); View(y); View(z)
write.table(z, "uniqueTr.sig.txt", row.names = F, col.names = T, quote = F, sep = "\t")
median(z$Freq)
dim(z[which(z$Freq > 12),])
z = z[which(z$Freq > 12),]
write.table(z, "uniqueTr.sig.aboveMed12.txt", row.names = F, col.names = T, quote = F, sep = "\t")

t = x[which(x$t %in% (setdiff(x$t, intersect(x$t, y$t)))),]

# volcano plot
plot((results$cor_Wet) - (results$cor_Salt),-log10(results$p_value),bty='n',col=as.factor(results$qvalue<0.05),xlab='difference in Spearman correlation (Wet-Salt)',ylab='-log10, p-value')

t=res_sig$tr1; t=append(t, res_sig$tr2); 
length(unique(t))




# qq-plot
qqplot(-log10(runif(dim(pairs)[1])),-log10(plasticdeco$p_value),pch=20,xlab='p-value, uniform distribution',ylab='p-value, differential correlation test',bty='n',xlim=c(0,max(-log10(plasticdeco$p_value))),ylim=c(0,max(-log10(plasticdeco$p_value))))
x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)

# example PLot of a pair of transcripts
i<-which( results$p_value==min(results$p_value))
plot( e_all_norm[pairs$V1[i],which(predictor==0)] , e_all_norm[pairs$V2[i],which(predictor==0)] ,bty='n',xlab='Normalized expression of gene 1',ylab='Normalized expression of gene 2',main='Year 2016');
abline(lm ( e_all_norm[pairs$V2[i],which(predictor==0)] ~ e_all_norm[pairs$V1[i],which(predictor==0)] ),lty=2)

plot( e_all_norm[pairs$V1[i],which(predictor==1)] , e_all_norm[pairs$V2[i],which(predictor==1)] ,bty='n',xlab='Normalized expression of gene 1',ylab='Normalized expression of gene 2',main='Year 2017');
abline(lm ( e_all_norm[pairs$V2[i],which(predictor==1)] ~ e_all_norm[pairs$V1[i],which(predictor==1)] ),lty=2)


