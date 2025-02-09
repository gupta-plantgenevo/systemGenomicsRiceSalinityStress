#######This prepares the input file for selection analysis on other physiological traits
library(lme4)
library(lmerTest)
library(dplyr)
library(stringr)
library(corrplot)

###Functions that will be used by the script
#Function to find the outliers that are beyond 3SD from the normalized data
findOutlier <- function(data, cutoff = 3) {
  ## Calculate the sd
  sds <- apply(data, 2, sd, na.rm = TRUE)
  avg <- apply(data, 2, mean, na.rm =TRUE)
  ## Identify the cells with value greater than cutoff * sd (column wise)
  result <- mapply(function(d, m, s) {
    which(abs((d - m)/s) > cutoff)
  }, data, avg, sds)
  result
}

#Function to filter out the outliers that have been identified using the above function (will add NA inplace of the outliers)
removeOutlier <- function(data, outliers) {
  result <- mapply(function(d, o) {
    res <- d
    res[o] <- NA
    return(res)
  }, data, outliers)
  return(as.data.frame(result))
}

##Seelction Function
selestimate_fecun <- function(fecun){
  S = c()
  C = c()
  Pval_S = c()
  Serr = c()
  Pval_C = c()
  tval_S = c()
  tval_C = c()
  Cerr = c()
  ##This loop will fit a linear and quadratic model and then save the S and C values (regression coeff) and the respective P values in separate vectors
  ##Takes a while to run
  for(i in 4:ncol(fecun)){
    tr = fecun[,i]
    linear_model = lmerTest::lmer(Fecundity ~ tr + (1|Block), data = fecun)
    Serr = append(Serr, summary(linear_model)$coefficients[2,2])
    S = append(S, summary(linear_model)$coefficients[2,1])
    tval_S = append(tval_S, summary(linear_model)$coefficients[2,4])
    Pval_S = append(Pval_S, summary(linear_model)$coefficients[2,5])
    tr2 = tr^2
    quad_model = lmer(Fecundity ~ tr + tr2 + (1|Block), data = fecun)
    Cerr = append(Cerr, summary(quad_model)$coefficients[3,2])
    C = append(C, 2*(summary(quad_model)$coefficients[3,1]))
    tval_C = append(tval_C, (summary(quad_model)$coefficients[3,4]))
    Pval_C = append(Pval_C, summary(quad_model)$coefficients[3,5])
  }
  output = data.frame(S,Serr,tval_S,Pval_S,C,Cerr,tval_C,Pval_C)
  return(output)
}

###Doing it for INDICA -- Saline Fecundity:
##Reading the files:
fitness_all = read.table("ricesalinity/data/Fitness/Indica_salt_fitness.txt",  header=TRUE, sep="\t", check.names = FALSE)
salt_phys = read.table("ricesalinity/data/Physio_otherMeasurements/salt_rawdata.txt", header=TRUE, sep="\t", check.names = FALSE)

salt_phys$chl_b_T2 = as.numeric(salt_phys$chl_b_T2)
salt_phys$Na1.K1 = as.numeric(salt_phys$Na_T1/salt_phys$K_T1)
salt_phys$Na2.K2 = as.numeric(salt_phys$Na_T2/salt_phys$K_T2)

ind_salt = salt_phys %>% filter(Plot %in% ind_salt_fitness_all$Plot)

ind_salt$Trt = rep("SS", times=nrow(ind_salt))
data_for_cov = na.omit(ind_salt)
corrmat = cor(data_for_cov[,c(6:20)])
corrplot(corrmat, method = 'number')
#So, due to high correlation between traits (>0.6) I can exclude some of them and focus on the ones that are not correlated.
#Also, probably should explore just taking Na/K ratios instead of Na and K individually
#Later will try removing T1 data points (since nonsig selection was observed in the later analysis on T1s)

##The first step is removing individuals with 0 fecundity individuals:
names(ind_salt_fitness_all)[12] = "Fecundity"
ind_salt_fitness_filtered = ind_salt_fitness_all %>% filter(Fecundity > 0 )
str(ind_salt_fitness_filtered)       ##removes 24 individuals with 0 fecundity of a total of 389 individuals; Now has a total of 365 individuals

#Standardize the fecundity fitness
ind_salt_fitness_filtered$Fecundity = ind_salt_fitness_filtered$Fecundity/(mean(ind_salt_fitness_filtered$Fecundity))
ind_salt_fitness_filtered = ind_salt_fitness_filtered[,c(8,12,10)]

##We then need to remove individuals with 0 fecundity from the RNASeq DataSet as well
ind_salt_sel = ind_salt %>% filter(Plot %in% ind_salt_fitness_filtered$Plot)
#str(ind_geno_salt_sel)

##Next we need to standardize the phenotypes (xi - mean(x)/sd(s)) 
#(I  think I will find and remove outliers for phenotypic measurements, seems like what is done in lit as well-- Confirm with Niels)
ind_salt_sel = ind_salt_sel[,c(1,6:18)]
tmp = merge(ind_salt_fitness_filtered, ind_salt_sel, by="Plot")
ind_salt_sel = tmp
ind_salt_std = as.data.frame(scale(ind_salt_sel[,-c(1:3)]))            ##Standardizing all columns. NOTE that scale removes NA by default
ind_salt_std$Plot = ind_salt_sel$Plot

##Removing from each column the relative (or standardized) values that are outside of 3SD limits
outliers = findOutlier(ind_salt_std)
ind_salt_standardized = removeOutlier(ind_salt_std, outliers)
ind_salt_std = ind_salt_standardized

##Now here I will merge the RNASeq data with fitness data
tmp = merge(ind_salt_fitness_filtered, ind_salt_std, by=c("Plot"))
#str(tmp)
total = ncol(tmp)

data_forSel_analysis_salt_IND = tmp
write.table(data_forSel_analysis_salt_IND, "ricesalinity/data/SelectionAnalysis/InputFiles/IND_salt_fecundity.physio.input.new.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

#Running Selection
fecundity = read.table("ricesalinity/data/SelectionAnalysis/InputFiles/IND_salt_fecundity.physio.input.new.txt", header=TRUE, sep="\t", check.names = FALSE)
fecundity = tmp
Pheno = colnames(fecundity)[-c(1:3)]
fecundity$Block = as.factor(fecundity$Block)
#tmp = sel_coeff
sel_coeff = selestimate_fecun(fecundity)
str(sel_coeff)
sel_coeff$Pheno = Pheno
sel_coeff = sel_coeff[,c(9,1:8)]
#str(sel_coeff)

##Beta
linear_model = lmer(Fecundity ~ V_LOP + Na_T1 + K_T1 + chl_a_T1 + Na_T2 + chl_a_T2 + Daysto50Flower + (1|Block), data = fecundity)
summary(linear_model)$coefficients

#gamma
t1 = (fecundity$V_LOP)^2; t2 = (fecundity$Na_T1)^2; t3 = (fecundity$K_T1)^2; t4 = (fecundity$chl_a_T1)^2; 
t5 = (fecundity$Na_T2)^2; t6 = (fecundity$chl_a_T2)^2; t7 = (fecundity$Daysto50Flower)^2
quad_model = lmer(Fecundity ~ V_LOP + t1 + Na_T1 + t2 + K_T1 + t3 + chl_a_T1 + t4 + Na_T2 + t5 + chl_a_T2 + t6 + Daysto50Flower + t7 + (1|Block), data = fecundity)
summary(quad_model)$coefficients


data_for_cov = ind_salt %>% filter(Plot %in% ind_salt_fitness_all$Plot)
data_for_cov = na.omit(data_for_cov)
data_for_cov = data_for_cov[,-c(1:5)]
covmat = cov(data_for_cov)    #Variance covariance mat 
test = solve(covmat)  #inverse 
test = as.matrix(test)
tmp = as.matrix(sel_coeff$S)
beta = test %*% tmp
#gamma = (test %*% tmp) %*% test 
sel_coeff$beta = as.vector(beta)
write.table(sel_coeff, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/Physio_SelCoeff_Indica_Salt.block.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
corrmat = cor(data_for_cov)

###Doing it for INDICA -- Wet Fecundity:
##Reading the files:
ind_wet_fitness_all = read.table("ricesalinity/data/Fitness/Indica_normal_fitness.txt",  header=TRUE, sep="\t", check.names = FALSE)
ind_wet = read.table("ricesalinity/data/Physio_otherMeasurements/wet_rawdata.txt", header=TRUE, sep="\t", check.names = FALSE)
ind_wet = ind_wet %>% filter(Plot %in% ind_wet_fitness_all$Plot)
ind_wet$Na1.K1 = as.numeric(ind_wet$Na_T1/ind_wet$K_T1)
ind_wet$Na2.K2 = as.numeric(ind_wet$Na_T2/ind_wet$K_T2)
ind_wet$Daysto50Flower = as.numeric(ind_wet$Daysto50Flower)

data_for_cov = na.omit(ind_wet)
corrmat = cor(data_for_cov[,c(6:18)])
corrplot(corrmat, method = 'number')

##The first step is removing individuals with 0 fecundity individuals:
names(ind_wet_fitness_all)[12] = "Fecundity"
ind_wet_fitness_filtered = ind_wet_fitness_all %>% filter(Fecundity > 0 )
str(ind_wet_fitness_filtered)       ##Now has a total of 384 individuals from 410

#Standardize the fecundity fitness
ind_wet_fitness_filtered$Fecundity = ind_wet_fitness_filtered$Fecundity/(mean(ind_wet_fitness_filtered$Fecundity))
ind_wet_fitness_filtered = ind_wet_fitness_filtered[,c(8,12,10)]

##We then need to remove individuals with 0 fecundity from the RNASeq DataSet as well
ind_wet_sel = ind_wet %>% filter(Plot %in% ind_wet_fitness_filtered$Plot)
#str(ind_geno_salt_sel)

##Next we need to standardize the phenotypes (xi - mean(x)/sd(s)) 
#(I dont think I will find and remove outliers for phenotypic measurements -- Confirm with Niels)
ind_wet_sel = ind_wet_sel[,c(1,6:18)]
tmp = merge(ind_wet_fitness_filtered, ind_wet_sel, by="Plot")
ind_wet_sel = tmp
ind_wet_std = as.data.frame(scale(ind_wet_sel[,-c(1:3)]))            ##Standardizing all columns. NOTE that scale removes NA by default
ind_wet_std$Plot = ind_wet_sel$Plot

##Removing from each column the relative (or standardized) values that are outside of 3SD limits
outliers = findOutlier(ind_wet_std)
ind_wet_standardized = removeOutlier(ind_wet_std, outliers)
ind_wet_std = ind_wet_standardized

##Now here I will merge the RNASeq data with fitness data
tmp = merge(ind_wet_fitness_filtered, ind_wet_std, by=c("Plot"))
#str(tmp)
total = ncol(tmp)

data_forSel_analysis_wet_IND = tmp
write.table(data_forSel_analysis_wet_IND, "ricesalinity/data/SelectionAnalysis/InputFiles/IND_wet_fecundity.physio.input.new.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

#Running Selection
fecundity = read.table("ricesalinity/data/SelectionAnalysis/InputFiles/IND_wet_fecundity.physio.input.new.txt", header=TRUE, sep="\t", check.names = FALSE)
fecundity = tmp
Pheno = colnames(fecundity)[-c(1:3)]
fecundity$Block = as.factor(fecundity$Block)
#tmp = sel_coeff
sel_coeff = selestimate_fecun(fecundity)
str(sel_coeff)
sel_coeff$Pheno = Pheno
sel_coeff = sel_coeff[,c(9,1:8)]
#str(sel_coeff)

##Beta
linear_model = lmer(Fecundity ~ V_LOP + Na_T1 + K_T1 + chl_a_T1 + Na_T2 + chl_a_T2 + Daysto50Flower + (1|Block), data = fecundity)
summary(linear_model)$coefficients

#gamma
t1 = (fecundity$V_LOP)^2; t2 = (fecundity$Na_T1)^2; t3 = (fecundity$K_T1)^2; t4 = (fecundity$chl_a_T1)^2; 
t5 = (fecundity$Na_T2)^2; t6 = (fecundity$chl_a_T2)^2; t7 = (fecundity$Daysto50Flower)^2
quad_model = lmer(Fecundity ~ V_LOP + t1 + Na_T1 + t2 + K_T1 + t3 + chl_a_T1 + t4 + Na_T2 + t5 + chl_a_T2 + t6 + Daysto50Flower + t7 + (1|Block), data = fecundity)
summary(quad_model)$coefficients


##I just printed and saved the results in a text file
#For S/C and beta/gamma.


