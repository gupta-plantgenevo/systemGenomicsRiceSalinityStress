library(rstatix)
library(dplyr)
library(coin)
library(cowplot)
library(ggplot2)

sel_coeff_wet = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_normal.w0fecun.nonStd.txt", header=TRUE, sep = "\t")
sel_coeff_salt = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_salt.w0fecun.nonStd.txt", header=TRUE, sep = "\t")

###2. S comparison -- direction and magnitude
##2.1. Just getting the # tr and S values
#In wet
median(abs(sel_coeff_wet$S))                 #Median of S in wet
count(sel_coeff_wet[which(sel_coeff_wet$S > 0),])     #Number of transcripts with positive S
count(sel_coeff_wet[which(sel_coeff_wet$S < 0),])     #Number of transcripts with negative S
median(sel_coeff_wet[which(sel_coeff_wet$S > 0),]$S)  #Median of positive S values
median(sel_coeff_wet[which(sel_coeff_wet$S < 0),]$S)  #Median of negative S values

#In Salt
median(abs(sel_coeff_salt$S))                #Median of S in salt
count(sel_coeff_salt[which(sel_coeff_salt$S > 0),])
count(sel_coeff_salt[which(sel_coeff_salt$S < 0),])
median(sel_coeff_salt[which(sel_coeff_salt$S > 0),]$S)
median(sel_coeff_salt[which(sel_coeff_salt$S < 0),]$S)

##Testing whether S in wet is significantly different that S in Salt (overall S, not just the magnistude)
ks.test(sel_coeff_wet$S, sel_coeff_salt$S, alternative ="less") 
#Wet S is significantly less than salt S; D = 0.30955; p-val < 2.2e-16

##2.2 Test whether positive values are significantly different in magnitude than negative values within the same treatment
#Wet
Selwet_S = sel_coeff_wet[which(sel_coeff_wet$S > 0),]$S
type_S = rep("Pos",times=length(Selwet_S))
Selwet_S = append(Selwet_S, abs(sel_coeff_wet[which(sel_coeff_wet$S < 0),]$S))
type_S = append(type_S, rep("Neg",times=length(sel_coeff_wet[which(sel_coeff_wet$S < 0),]$S)))
type_S = as.factor(type_S)

#Test whether one group is significantly higher/lower than the other
#From the plot is looks like positive values might be significant greater than the negative values 
wilcox.test(abs(Selwet_S)~type_S, na.rm=TRUE, exact = FALSE, alternative = "less")


##This tells us that positive directional selection is stronger than negative seelction in Wet conditions(Mann Whitney test, z-value = -4.4949, P-val = 3.48e-06)

#Salt Condition
Selsalt_S = sel_coeff_salt[which(sel_coeff_salt$S > 0),]$S
type_S = rep("Pos",times=length(Selsalt_S))
Selsalt_S = append(Selsalt_S, abs(sel_coeff_salt[which(sel_coeff_salt$S < 0),]$S))
type_S = append(type_S, rep("Neg",times=length(sel_coeff_salt[which(sel_coeff_salt$S < 0),]$S)))
type_S = as.factor(type_S)

#Test whether one group is significantly higher/lower than the other
wilcox.test(abs(Selsalt_S)~type_S, na.rm=TRUE,exact = FALSE)
##This tells us that positive directional selection is stronger under salinity  (Mann Whitney test, z-value = -7.813, P-val = 2.791e-15)


##2.3 Test whether positive (and negative) selection is stronger in wet vs salt
#Positive Selection
tmp = sel_coeff_wet[which(sel_coeff_wet$S > 0),]$S
tmp2 = rep("Wet",times=length(tmp))
tmp = append(tmp, sel_coeff_salt[which(sel_coeff_salt$S > 0),]$S)
tmp2 = append(tmp2, rep("Salt",times=length(sel_coeff_salt[which(sel_coeff_salt$S > 0),]$S)))
tmp2 = as.factor(tmp2)

#Test whether pos. sel is significantly higher/lower in wet vs salt
wilcox.test(tmp~tmp2, na.rm=TRUE, paired=FALSE, exact = FALSE)
##This tells us that positive directional selection is stronger under salinity than under wet (Mann Whitney test, z-value = 5.9709, P-val = 1.179e-09)

#Negative Selection
tmp = sel_coeff_wet[which(sel_coeff_wet$S < 0),]$S
tmp2 = rep("Wet",times=length(tmp))
tmp = append(tmp, sel_coeff_salt[which(sel_coeff_salt$S < 0),]$S)
tmp2 = append(tmp2, rep("Salt",times=length(sel_coeff_salt[which(sel_coeff_salt$S < 0),]$S)))
tmp2 = as.factor(tmp2)

#Test whether neg. sel is significantly higher/lower in wet vs salt
wilcox.test(abs(tmp)~tmp2, na.rm=TRUE, exact = FALSE)
##This tells us that negative directional selection is stronger under salinity than under wet (Mann Whitney test, z-value = -2.3609, P-val = 0.009115)


##Looking into if and what transcripts are significant even after bonferroni correction
## S
##Wet
s_wet_correctedPval = p.adjust(sel_coeff_wet$Pval_S, method = "bonferroni")
sel_coeff_wet$Pval_S.corrected = s_wet_correctedPval
tmp = sel_coeff_wet[which(sel_coeff_wet$Pval_S.corrected < 0.05),]
Wet_S_sig = tmp
nrow(Wet_S_sig[which(Wet_S_sig$S < 0),])

##Salt
s_salt_correctedPval = p.adjust(sel_coeff_salt$Pval_S, method = "bonferroni")
sel_coeff_salt$Pval_S.corrected = s_salt_correctedPval
tmp = sel_coeff_salt[which(sel_coeff_salt$Pval_S.corrected < 0.05),]



##2. Testing whether the magnitude of C is significant different between conditions -- this can be done but it doesnt really inform much since the signs means two very different things 

median(abs(wet_salt_mag[which(wet_salt_mag$Env == 'Wet'),]$SelMag_C))
median(abs(wet_salt_mag[which(wet_salt_mag$Env == 'Salt'),]$SelMag_C))


#Two-sided test to test whether different
wilcox_test(SelMag_C~Env, data = wet_salt_mag, na.rm=TRUE, paired=FALSE, alternative = "two.sided", exact = FALSE)
##This tell us that the magnitude of selection varies significantly between the env (Two-sided Mann Whitney U test, P-val < 2.2e16)

#Test whether one group is significantly higher/lower than the other
#From the plot is looks like saline might be stronger
wilcox_test(SelMag_C~Env, data = wet_salt_mag, na.rm=TRUE, paired=FALSE, alternative = "greater", exact = FALSE)
##This tells us that selection under salinity is stronger (Mann Whitney test, z-value = 20.158, P-val < 2.2e-16)

###2. C comparison -- direction and magnitude
##2.1. Just getting the # tr and C values
#In wet
median((sel_coeff_wet$C))                 #Median of S in wet
count(sel_coeff_wet[which(sel_coeff_wet$C > 0),])     #Number of transcripts with positive S
count(sel_coeff_wet[which(sel_coeff_wet$C < 0),])     #Number of transcripts with negative S
median(sel_coeff_wet[which(sel_coeff_wet$C > 0),]$C)  #Median of positive S values
median(sel_coeff_wet[which(sel_coeff_wet$C < 0),]$C)  #Median of negative S values

#In Salt
median((sel_coeff_salt$C))                #Median of S in salt
count(sel_coeff_salt[which(sel_coeff_salt$C > 0),])
count(sel_coeff_salt[which(sel_coeff_salt$C < 0),])
median(sel_coeff_salt[which(sel_coeff_salt$C > 0),]$C)
median(sel_coeff_salt[which(sel_coeff_salt$C < 0),]$C)

##Testing whether S in wet is significantly different that S in Salt (overall S, not just the magnistude)
ks.test(sel_coeff_wet$C, sel_coeff_salt$C, alternative ="less") 
#Wet S is significantly less than salt S; D = 0.078192; p-val < 2.2e-16

##2.2 Test whether positive values are significantly different in magnitude than negative values within the same treatment
#Wet
Selwet_C = sel_coeff_wet[which(sel_coeff_wet$C > 0),]$C
type_C = rep("Pos",times=length(Selwet_C))
Selwet_C = append(Selwet_C, sel_coeff_wet[which(sel_coeff_wet$C < 0),]$C)
type_C = append(type_C, rep("Neg",times=length(sel_coeff_wet[which(sel_coeff_wet$C < 0),]$C)))
type_C = as.factor(type_C)

#Test whether one group is significantly higher/lower than the other
wilcox_test(abs(Selwet_C)~type_C, na.rm=TRUE, paired=FALSE, alternative = "greater", exact = FALSE)
##This tells us that negative C (stabilizing sel. is stronger than disruptive selection) is stronger in Wet conditions (Mann Whitney test, z-value = 10.772, P-val < 2.2e-16)

#Salt Condition
Selsalt_C = sel_coeff_salt[which(sel_coeff_salt$C > 0),]$C
type_C = rep("Pos",times=length(Selsalt_C))
Selsalt_C = append(Selsalt_C, sel_coeff_salt[which(sel_coeff_salt$C < 0),]$C)
type_C = append(type_C, rep("Neg",times=length(sel_coeff_salt[which(sel_coeff_salt$C < 0),]$C)))
type_C = as.factor(type_C)

#Test whether one group is significantly higher/lower than the other
wilcox_test(abs(Selsalt_C)~type_C, na.rm=TRUE, paired=FALSE, alternative = "two.sided", exact = FALSE)
##No significant difference between stabilizing or disruptive sel in the salt condition -- Mann-Whitney two-sided P-val = 0.6283


##2.3 Test whether positive (and negative) selection is stronger in wet vs salt
#Positive Selection
tmp = sel_coeff_wet[which(sel_coeff_wet$C > 0),]$C
tmp2 = rep("Wet",times=length(tmp))
tmp = append(tmp, sel_coeff_salt[which(sel_coeff_salt$C > 0),]$C)
tmp2 = append(tmp2, rep("Salt",times=length(sel_coeff_salt[which(sel_coeff_salt$C > 0),]$C)))
tmp2 = as.factor(tmp2)

#Test whether pos. sel is significantly higher/lower in wet vs salt
wilcox_test(tmp~tmp2, na.rm=TRUE, paired=FALSE, alternative = "greater", exact = FALSE)
##This tells us that positive C (disruptive selection) is stronger under salinity than under wet (Mann Whitney test, z-value = 17.947, P-val < 2.2e-16)

#Negative Selection
tmp = sel_coeff_wet[which(sel_coeff_wet$C < 0),]$C
tmp2 = rep("Wet",times=length(tmp))
tmp = append(tmp, sel_coeff_salt[which(sel_coeff_salt$C < 0),]$C)
tmp2 = append(tmp2, rep("Salt",times=length(sel_coeff_salt[which(sel_coeff_salt$C < 0),]$C)))
tmp2 = as.factor(tmp2)

#Test whether neg. sel is significantly higher/lower in wet vs salt
wilcox.test(abs(tmp)~tmp2, na.rm=TRUE, paired=FALSE,  exact = FALSE)
##This tells us that negative C (stabilizing selection) is stronger under salinity than under wet (Mann Whitney test, z-value = 11.308, P-val = < 2.2e-16)
##Looking into if and what transcripts are significant even after bonferroni correction
## C
##Wet
c_wet_correctedPval = p.adjust(sel_coeff_wet$Pval_C, method = "bonferroni")
sel_coeff_wet$Pval_C.corrected = c_wet_correctedPval
tmp = sel_coeff_wet[which(sel_coeff_wet$Pval_C.corrected < 0.05),]

##Salt
c_salt_correctedPval = p.adjust(sel_coeff_salt$Pval_C, method = "bonferroni")
sel_coeff_salt$Pval_C.corrected = c_salt_correctedPval
tmp = sel_coeff_salt[which(sel_coeff_salt$Pval_C.corrected < 0.05),]
