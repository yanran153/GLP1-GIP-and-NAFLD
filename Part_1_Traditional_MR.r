
################################ Part 1 Traditional Mendelian Randomization ################################

#This part aims to explore the relationship between circulating GIP/GLP1 and NAFLD as well as its related traits,
#And this parts contains three steps, MR analysis,sensitive analysis and result combination, among which we used inverse variance weighted to combine the results from discovery and validation datasets      
#In the subsequent analysis, we exclude 'GLP1_fasing' due to absense of significat instrument variables
#The code for leave-one-out analysis could be found in Part S1 leave-one-out.r

######################## Step 1 MR Analysis ########################

library(TwoSampleMR)
library(gtools)
library(mr.raps)
library(penalized)
library(data.table)

################ Define the Function we used ################

mr_circlulating_binary = function(data_exposure,exposure_name,data_outcome,outcome_name){
  #This function is used to perform Mendelian Randomization on binary outcome
  #parameter@data_exposure: Output from format_data()
  #parameter@exposure_name: Name of the exposure
  #parameter@data_outcome: Output from format_data()
  #parameter@outcome_name: Name of outcome
  #output: A dataframe including the result from IVW-R/weighted median/MR_RAPS  or from Wald Raios
  
  dat = harmonise_data(data_exposure,data_outcome,action=3)
  dat = subset(dat,mr_keep==TRUE)
if(dim(dat)[1]>2){
  res = mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_weighted_median","mr_egger_regression"))
  or = generate_odds_ratios(res)
  mr_frame = dat[,c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome")]
  names(mr_frame)[c(1:5)] = c("rsid","beta.x","se.x","beta.y","se.y")
  MR_RAPS  =  mr.raps.overdispersed.robust(mr_frame$beta.x, mr_frame$beta.y, mr_frame$se.x, mr_frame$se.y,loss.function = "huber", k = 1.345, 
  initialization = c("l2"), suppress.warning = FALSE, diagnosis = FALSE, niter = 20, tol = .Machine$double.eps^0.5)
  or$b = or$or
  or$lo_ci = or$or_lci95 
  or$up_ci = or$or_uci95 
  or = or[,-c(12,13,14)]
  
  or[5,]= NA
  or[5,c(1,2,3,4,6)]= or[2,c(1,2,3,4,6)]
  or[,"method"] = as.character(or[,"method"])
  or[5,"method"]="MR-RAPS"
  or[5,c("b","se","pval","lo_ci","up_ci")]=c(exp(MR_RAPS$beta.hat),MR_RAPS$beta.se,MR_RAPS$beta.p.value,
   											 exp(MR_RAPS$beta.hat-1.96*MR_RAPS$beta.se),exp(MR_RAPS$beta.hat+1.96*MR_RAPS$beta.se))
  or$outcome = outcome_name
  or$exposure = exposure_name  
  return(or)																	
}else{
  res = mr(dat)
  or = generate_odds_ratios(res)
  or$b = or$or
  or$lo_ci = or$or_lci95 
  or$up_ci = or$or_uci95 
  or = or[,-c(12,13,14)]
  or$outcome = outcome_name
  or$exposure = exposure_name 
  return(or)
}
}

mr_circlulating_continous = function(data_exposure,exposure_name,data_outcome,outcome_name){
  #This function is used to perform Mendelian Randomization on continuous outcome
  #parameter@data_exposure: Output from format_data()
  #parameter@exposure_name: Name of the exposure
  #parameter@data_outcome: Output from format_data()
  #parameter@outcome_name: Name of outcome
  #output: A dataframe including the result from IVW-R/weighted median/MR_RAPS or from Wald Raios
  
  dat = harmonise_data(data_exposure,data_outcome,action = 3)
  dat = subset(dat,mr_keep==TRUE)
if(dim(dat)[1]>=2){
  res = mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_weighted_median","mr_egger_regression"))
  mr_frame = dat[,c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome")]
  names(mr_frame)[c(1:5)] = c("rsid","beta.x","se.x","beta.y","se.y")
  MR_RAPS  =  mr.raps.overdispersed.robust(mr_frame$beta.x, mr_frame$beta.y, mr_frame$se.x, mr_frame$se.y,loss.function = "huber", k = 1.345, 
  initialization = c("l2"), suppress.warning = FALSE, diagnosis = FALSE, niter = 20, tol = .Machine$double.eps^0.5)
  
  res$lo_ci = res$b-1.96*res$se
  res$up_ci = res$b+1.96*res$se
  res[5,]= NA
  res[5,c(1,2,3,4,6)]= res[2,c(1,2,3,4,6)]
  res[,"method"] = as.character(res[,"method"])
  res[5,"method"]="MR-RAPS"
  res[5,c("b","se","pval","lo_ci","up_ci")]=c(MR_RAPS$beta.hat,MR_RAPS$beta.se,MR_RAPS$beta.p.value, (MR_RAPS$beta.hat-1.96*MR_RAPS$beta.se),
                                                                        (MR_RAPS$beta.hat+1.96*MR_RAPS$beta.se))
  res$outcome = outcome_name
  res$exposure = exposure_name  
  return(res)																	
}else{
  res = mr(dat)
  res$outcome = outcome_name
  res$exposure = exposure_name 
  res$lo_ci = res$b-1.96*res$se
  res$up_ci = res$b+1.96*res$se
  return(res)
}
}

################ NAFLD ################

exposure_analysis = c("GLP1_2h","GIP_2h","GIP_fasting")
outcome_file = c("NAFLD_Discovery_Processed.txt","NAFLD_Validation.txt") #This is the file name of preprocessed outcome data
Result_NAFLD = NULL
setwd("~/GLP1/Analysis/Data/Outcome/NAFLD") #This is the file path of preprocessed outcome data
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   outcome_origin = fread(outcome)
   for(j in 1:length(exposure_analysis)){
     #exposure = fread(paste0("~/GLP1/Analysis/Data/Exposure/Cirrculating_concentration/",exposure_analysis[j],"_reset.txt")) #This is the file path of preprocessed exposure data
	 exposure_name = exposure_analysis[j]
	 print(exposure_name)
	 #exposure = exposure[exposure$EAF>0.01 & exposure$EAF<0.99,]
     #exposure = subset(exposure,exposure$P<5e-08)
     #exposure = format_data(exposure,snp_col = "SNP",beta_col = "Beta",se_col = "SE",pval_col = "P",effect_allele_col = "EA",other_allele_col = "NEA",
                            #eaf_col = "EAF",chr_col = "Chr",pos_col = "BP")
     #exposure_dat = clump_data(exposure,clump_kb = 10000,clump_r2 = 0.01,clump_p1 = 1,clump_p2 = 1)
	 exposure_dat = fread(paste0("~/GLP1/Analysis/Data/Exposure/Cirrculating_concentration/",exposure_analysis[j],"_clump.txt"))
	 if(outcome=="NAFLD_Discovery_Processed.txt"){
	    outcome_name = "NAFLD_Discovery_Processed"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = mr_circlulating_binary(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD = rbind(Result_NAFLD,result)
	 }else if(outcome=="NAFLD_Validation.txt"){
	    outcome_name = "NAFLD_Validation"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,variant_id%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "variant_id",beta_col = "lnOR",se_col = "standard_error",effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",pval_col = "p_value")
		result = mr_circlulating_binary(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD = rbind(Result_NAFLD,result)
	 }
   }
}

################ NAFLD related Traits ################& exposure_name!="GIP_2h"

exposure_analysis = c("GLP1_2h","GIP_2h","GIP_fasting")
outcome_file = c("Cirrhosis_Processed.txt","PDFF_Processed.txt","GGT.tsv.gz","ALT.tsv.gz","ALP.tsv.gz")  #This is the file name of preprocessed outcome data
Result_NAFLD_Related = NULL
setwd("~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits") #This is the file path of preprocessed outcome data
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   outcome_origin = fread(outcome)
   for(j in 1:length(exposure_analysis)){
     #exposure = fread(paste0("~/GLP1/Analysis/Data/Exposure/Cirrculating_concentration/",exposure_analysis[j],"_reset.txt")) #This is the file path of preprocessed exposure data
	 exposure_name = exposure_analysis[j]
	 #exposure = exposure[exposure$EAF>0.01 & exposure$EAF<0.99,]
     #exposure<-subset(exposure,exposure$P<5e-08)
     #exposure = format_data(exposure,snp_col = "SNP",beta_col = "Beta",se_col = "SE",pval_col = "P",effect_allele_col = "EA",other_allele_col = "NEA",
                            #eaf_col = "EAF",chr_col = "Chr",pos_col = "BP")
     #exposure_dat = clump_data(exposure,clump_kb = 10000,clump_r2 = 0.01,clump_p1 = 1,clump_p2 = 1)
	 exposure_dat = fread(paste0("~/GLP1/Analysis/Data/Exposure/Cirrculating_concentration/",exposure_analysis[j],"_clump.txt"))
	 if(outcome=="Cirrhosis_Processed.txt"){
	    outcome_name = strsplit(outcome,"_Processed.txt")[[1]][1]
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = mr_circlulating_binary(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD_Related = rbind(Result_NAFLD_Related,result)
	 }else if(outcome=="PDFF_Processed.txt"){
	    outcome_name = "PDFF"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Effect",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = mr_circlulating_continous(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD_Related = rbind(Result_NAFLD_Related,result)
	 }else if(outcome=="GGT.tsv.gz"|outcome=="ALT.tsv.gz"|outcome=="ALP.tsv.gz"){
	    outcome_name = strsplit(outcome,".tsv.gz")[[1]][1]
		print(outcome_name)
		outcome_match = subset(outcome_origin,SNP%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",
		                         pval_col = "P",eaf_col="A1FREQ")
		result = mr_circlulating_continous(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD_Related = rbind(Result_NAFLD_Related,result)
	 }
	 
   }
}


######################## Step 2 Sensitive Analysis ########################

library(TwoSampleMR)
library(data.table)

################ Define the Function we used ################

mr_senstive = function(data_exposure,exposure_name,data_outcome,outcome_name){
  #This function is used to perform sensitive analysis of Mendelian Randomization
  #parameter@data_exposure: Output from format_data()
  #parameter@exposure_name: Name of the exposure
  #parameter@data_outcome: Output from format_data()
  #parameter@outcome_name: Name of outcome
  #output: A dataframe including the sensitivity result
  
  dat = harmonise_data(data_exposure,data_outcome,action=3)
  dat = subset(dat,mr_keep==TRUE)
  if(dim(dat)[1]>2){
    het = mr_heterogeneity(dat)
	het = subset(het, method=="Inverse variance weighted")
    plt = mr_pleiotropy_test(dat)
	sensitive_result = c(exposure_name,outcome_name,het$Q,het$Q_pval,plt$egger_intercept,plt$se,plt$pval)
	return(sensitive_result)
  }
}

################ Conduct the senstive analysis ################

######## Pleiotropy and Heterogeneity #######
exposure_analysis = c("GLP1_2h","GIP_2h","GIP_fasting")
outcome_file = c("NAFLD_Discovery_Processed.txt","NAFLD_Validation.txt","Cirrhosis_Processed.txt","GGT.tsv.gz","ALT.tsv.gz","ALP.tsv.gz")
Result_sensitive = NULL
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   if(outcome=="NAFLD_Discovery_Processed.txt"|outcome=="NAFLD_Validation.txt"){
    setwd("~/GLP1/Analysis/Data/Outcome/NAFLD")
   }else{
    setwd("~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits")
   }
   outcome_origin = fread(outcome)
   for(j in 1:length(exposure_analysis)){
     #exposure = fread(paste0("~/GLP1/Analysis/Data/Exposure/Cirrculating_concentration/",exposure_analysis[j],"_reset.txt")) #This is the file path of preprocessed exposure data
	 exposure_name = exposure_analysis[j]
	 #exposure = exposure[exposure$EAF>0.01 & exposure$EAF<0.99,]
     #exposure<-subset(exposure,exposure$P<5e-08)
     #exposure = format_data(exposure,snp_col = "SNP",beta_col = "Beta",se_col = "SE",pval_col = "P",effect_allele_col = "EA",other_allele_col = "NEA",
                            #eaf_col = "EAF",chr_col = "Chr",pos_col = "BP")
     #exposure_dat = clump_data(exposure,clump_kb = 10000,clump_r2 = 0.01,clump_p1 = 1,clump_p2 = 1)
	 exposure_dat = fread(paste0("~/GLP1/Analysis/Data/Exposure/Cirrculating_concentration/",exposure_analysis[j],"_clump.txt"))
	 if(outcome=="NAFLD_Discovery_Processed.txt"|outcome=="Cirrhosis_Processed.txt"){
	    outcome_name = strsplit(outcome,"_Processed")[[1]][1]
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = mr_senstive(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_sensitive = rbind(Result_sensitive,result)
	 }else if(outcome=="NAFLD_Validation.txt"){
	    outcome_name = strsplit(outcome,"_Processed")[[1]][1]
	    outcome_match = subset(outcome_origin,variant_id%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "variant_id",beta_col = "lnOR",se_col = "standard_error",effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",pval_col = "p_value")
		result =  mr_senstive(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_sensitive = rbind(Result_sensitive,result)
	 }else if(outcome=="GGT.tsv.gz"|outcome=="ALT.tsv.gz"|outcome=="ALP.tsv.gz"){
	    outcome_name = strsplit(outcome,".tsv.gz")[[1]][1]
		outcome_match = subset(outcome_origin,SNP%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",
		                         pval_col = "P",eaf_col="A1FREQ")
		result = mr_senstive(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_sensitive = rbind(Result_sensitive,result)
	 }
	 
   }
}
colnames(Result_sensitive) = c("Exposure","Outcome","Q","Q_pval","Egger_intercept","SE","Intercept_pval")


######################## Step 3 Result Combination ########################

library(data.table)

IVW_meta = function(data1,data2){
   data1_use = data1[,c("outcome","exposure","method","b","se")]
   data1_use$b = log(data1_use$b)
   data2_use = data2[,c("b","se")]
   data2_use$b = log(data2_use$b)
   data_use = cbind(data1_use,data2_use)
   data_use = subset(data_use,method=="Inverse variance weighted (multiplicative random effects)"|method=="MR Egger"|method=="Weighted median"|method=="MR-RAPS"|method=="Wald ratio")
   colnames(data_use)=c("outcome","exposure","method","b","se","b.1","se.1")
   data_use$w1 = 1/(data_use$se)^2
   data_use$w2 = 1/(data_use$se.1)^2
   b.meta.up = data_use$b*data_use$w1+data_use$b.1*data_use$w2
   b.meta.down = apply(data_use[,c("w1","w2")],1,sum)
   data_use$b.meta = b.meta.up/b.meta.down
   data_use$se.meta  = sqrt(1/b.meta.down)
   data_use$Z = data_use$b.meta/data_use$se.meta
   data_use$P.meta = 2 * pnorm(abs(data_use$Z), lower.tail=FALSE)
   data_use$OR = exp(data_use$b.meta)
   data_use$CIU = exp(data_use$b.meta+1.96*data_use$se.meta)
   data_use$CID = exp(data_use$b.meta-1.96*data_use$se.meta)
   return(data_use)
}

result_NAFLD = fread("~/GLP1/Analysis/Traditional_MR/Result_NAFLD.csv")
data1 = subset(result_NAFLD,outcome=="NAFLD_Discovery_Processed")
data2 = subset(result_NAFLD,outcome=="NAFLD_Validation")

NAFLD_meta = IVW_meta(data1,data2)