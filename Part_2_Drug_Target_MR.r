
################################ Part 2 Drug Target Mendelian Randomization ################################

#This part aims to explore the relationship between drug target gene and NAFLD as well as its related traits.
#This part includes three steps, drug target MR using drug repsonse downstream biomarker, drug target mr using gene expressions from whole blood tissue, drug target mr using gene expressions from highly expressed tissues
#In this part, positive control analyses were also conducted to valid the instruments of downstream biomarker and gene expressions, respectivly
#In step 3, we analyzed all tissues where GLP1R/GIPR is significant. For simplicity, we just provide the code for the pancreas here and the code for another tissues could be found in Part S2
#The code for leave-one-out analysis could be found in Part_S1_leave-one-out.r

library(TwoSampleMR)
library(data.table)
library(mr.raps)


######################## Step 1 Drug Target MR using drug repsonse downstream biomarker ########################

mr_drug_binary = function(data_exposure,exposure_name,data_outcome,outcome_name){
  #This function is used to perform Mendelian Randomization on binary outcome
  #parameter@data_exposure: Output from format_data()
  #parameter@exposure_name: Name of the exposure
  #parameter@data_outcome: Output from format_data()
  #parameter@outcome_name: Name of outcome
  #output: A dataframe including the result from IVW-R/Weighted_Median/MR_Egger/MR_RAPS or from Wald Ratio
  #Note that the estimates here represented the effect on outcomes per S.D. change in HbA1c lowering
  
  dat = harmonise_data(data_exposure,data_outcome,action=3)
  dat = subset(dat,mr_keep==TRUE)
if(dim(dat)[1]>2){
  res = mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_weighted_median","mr_egger_regression"))
  res$b = -res$b
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
  or[5,c("b","se","pval","lo_ci","up_ci")]=c(exp(-MR_RAPS$beta.hat),MR_RAPS$beta.se,MR_RAPS$beta.p.value,
   											 exp(-MR_RAPS$beta.hat-1.96*MR_RAPS$beta.se),exp(-MR_RAPS$beta.hat+1.96*MR_RAPS$beta.se))
  or$outcome = outcome_name
  or$exposure = exposure_name  
  return(or)																	
}else{
  res = mr(dat)
  res$b = -res$b
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

mr_drug_continous = function(data_exposure,exposure_name,data_outcome,outcome_name){
  #This function is used to perform Mendelian Randomization on continuous outcome
  #parameter@data_exposure: Output from format_data()
  #parameter@exposure_name: Name of the exposure
  #parameter@data_outcome: Output from format_data()
  #parameter@outcome_name: Name of outcome
  #output: A dataframe including the result from IVW-R/weighted median/MR_RAPS or from Wald Ratio
  #Note that the estimates here represented the effect on outcomes per S.D. change in HbA1c lowering
  
  dat = harmonise_data(data_exposure,data_outcome,action = 3)
  dat = subset(dat,mr_keep==TRUE)
if(dim(dat)[1]>2){
  res = mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_weighted_median","mr_egger_regression"))
  res$b = -res$b
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
  res[5,c("b","se","pval","lo_ci","up_ci")]=c(-MR_RAPS$beta.hat,MR_RAPS$beta.se,MR_RAPS$beta.p.value, (-MR_RAPS$beta.hat-1.96*MR_RAPS$beta.se),
                                                                        (-MR_RAPS$beta.hat+1.96*MR_RAPS$beta.se))
  res$outcome = outcome_name
  res$exposure = exposure_name  
  return(res)																	
}else{
  res = mr(dat)
  res$b = -res$b
  res$outcome = outcome_name
  res$exposure = exposure_name 
  res$lo_ci = res$b-1.96*res$se
  res$up_ci = res$b+1.96*res$se
  return(res)
}
}

################ NAFLD ################

exposure_analysis = c("GLP1R","GIPR","Combine")
outcome_file = c("NAFLD_Discovery_Processed.txt","NAFLD_Validation.txt") #This is the file name of preprocessed outcome data
Result_NAFLD = NULL
setwd("~/GLP1/Analysis/Data/Outcome/NAFLD") #This is the file path of preprocessed outcome data
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   outcome_origin = fread(outcome)
   for(j in 1:length(exposure_analysis)){
	   exposure_name = exposure_analysis[j]
	   exposure_dat = fread(paste0("~/GLP1/Analysis/Data/Exposure/Downstream/HbA1c_",exposure_analysis[j],"_clump.txt"))
	 if(outcome=="NAFLD_Discovery_Processed.txt"){
	    outcome_name = "NAFLD_Discovery_Processed"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = mr_drug_binary(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD = rbind(Result_NAFLD,result)
	 }else if(outcome=="NAFLD_Validation.txt"){
	    outcome_name = "NAFLD_Validation"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,variant_id%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "variant_id",beta_col = "lnOR",se_col = "standard_error",effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",pval_col = "p_value")
		result = mr_drug_binary(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD = rbind(Result_NAFLD,result)
	 }
   }
}

################ NAFLD Related Traits ################

exposure_analysis = c("GLP1R","GIPR","Combine")
outcome_file = c("Cirrhosis_Processed.txt","PDFF_Processed.txt","GGT.tsv.gz","ALT.tsv.gz","ALP.tsv.gz") #This is the file name of preprocessed outcome data
Result_NAFLD_Related = NULL
setwd("~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits") #This is the file path of preprocessed outcome data
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   outcome_origin = fread(outcome)
   for(j in 1:length(exposure_analysis)){
	 exposure_name = exposure_analysis[j]
	 exposure_dat = fread(paste0("~/GLP1/Analysis/Data/Exposure/Downstream/HbA1c_",exposure_analysis[j],"_clump.txt"))
	 if(outcome=="Cirrhosis_Processed.txt"){
	    outcome_name = strsplit(outcome,"_Processed.txt")[[1]][1]
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = mr_drug_binary(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD_Related = rbind(Result_NAFLD_Related,result)
	 }else if(outcome=="PDFF_Processed.txt"){
	    outcome_name = "PDFF"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Effect",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result =  mr_drug_continous(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD_Related = rbind(Result_NAFLD_Related,result)
	 }else if(outcome=="GGT.tsv.gz"|outcome=="ALT.tsv.gz"|outcome=="ALP.tsv.gz"){
	    outcome_name = strsplit(outcome,".tsv.gz")[[1]][1]
		print(outcome_name)
		outcome_match = subset(outcome_origin,SNP%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",
		                         pval_col = "P",eaf_col="A1FREQ")
		result =  mr_drug_continous(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD_Related = rbind(Result_NAFLD_Related,result)
	 }
	 
   }
}

################ Postitive Control ################

exposure_analysis = c("GLP1R","GIPR","Combine")
setwd("~/GLP1/Analysis/Data/Exposure/Downstream")
outcome_total = fread("~/GLP1/Analysis/Data/Outcome/Postitive_Control/T2D_GWAS.txt")
Result_T2D = NULL
for(i in 1:length(exposure)){
  exposure_dat = fread(paste0("HbA1c_",exposure[i],"_clump.txt"))
  exposure_dat$id.exposure = exposure[i]
  outcome_match = subset(outcome_total,rsID%in%exposure_dat$SNP)
  if(dim(outcome_match)[1]!=0){
    outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsID",beta_col = "Fixed-effects_beta",se_col = "Fixed-effects_SE",
	                          effect_allele_col = "effect_allele",other_allele_col = "other_allele",
							  pval_col = " Fixed-effects_p-value",eaf_col = "effect_allele_frequency")
    dat = harmonise_data(exposure_dat,outcome_dat,action=3)
    dat = subset(dat,mr_keep==TRUE)
    res = mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_weighted_median","mr_egger_regression"))
	res$b = -res$b
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
    or[5,c("b","se","pval","lo_ci","up_ci")]=c(exp(-MR_RAPS$beta.hat),MR_RAPS$beta.se,MR_RAPS$beta.p.value,
   											 exp(-MR_RAPS$beta.hat-1.96*MR_RAPS$beta.se),exp(-MR_RAPS$beta.hat+1.96*MR_RAPS$beta.se))
    or$outcome = "T2D"
	or$exposure = exposure[i]
    het = mr_heterogeneity(dat)
    plt = mr_pleiotropy_test(dat)
	
	if(subset(het,method=="Inverse variance weighted")$Q_pval<=0.05){
      print(paste0("Heterogeneity:",exposure[i]))
    }else if((dim(dat)[1]>2)&(plt$pval<=0.05)){
      print(paste0("Pleiotropy:",exposure[i]))
    }
	
	Result_T2D = rbind(Result_T2D,or)
  }else{
    print(paste0("No corresponding SNPs of ",exposure[i]," was found in GWAS"))
  }
}


######################## Step 2 Drug Target MR using target gene expressions from whole blood tissue ########################

mr_drug_binary = function(data_exposure,exposure_name,data_outcome,outcome_name){
  #This function is used to perform Mendelian Randomization on binary outcome
  #parameter@data_exposure: Output from format_data()
  #parameter@exposure_name: Name of the exposure
  #parameter@data_outcome: Output from format_data()
  #parameter@outcome_name: Name of outcome
  #output: A dataframe including the result from IVW-R/Weighted_Median/MR_RAPS or from Wald Ratio
  
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

mr_drug_continous = function(data_exposure,exposure_name,data_outcome,outcome_name){
  #This function is used to perform Mendelian Randomization on continuous outcome
  #parameter@data_exposure: Output from format_data()
  #parameter@exposure_name: Name of the exposure
  #parameter@data_outcome: Output from format_data()
  #parameter@outcome_name: Name of outcome
  #output: A dataframe including the result from IVW-R/Weighted_Median/MR_RAPS or from Wald Ratio
  
  dat = harmonise_data(data_exposure,data_outcome,action = 3)
  dat = subset(dat,mr_keep==TRUE)
if(dim(dat)[1]>2){
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

exposure_analysis = c("GLP1R","GIPR")
outcome_file = c("NAFLD_Discovery_Processed.txt","NAFLD_Validation.txt") #This is the file name of preprocessed outcome data
Result_NAFLD = NULL
setwd("~/GLP1/Analysis/Data/Outcome/NAFLD") #This is the file path of preprocessed outcome data
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   outcome_origin = fread(outcome)
   for(j in 1:length(exposure_analysis)){
	   exposure_name = exposure_analysis[j]
	   exposure_dat = fread(paste0("~/GLP1/Analysis/Data/Exposure/Gene_expression/",exposure_analysis[j],".txt"))
	 if(outcome=="NAFLD_Discovery_Processed.txt"){
	    outcome_name = "NAFLD_Discovery_Processed"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = mr_drug_binary(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD = rbind(Result_NAFLD,result)
	 }else if(outcome=="NAFLD_Validation.txt"){
	    outcome_name = "NAFLD_Validation"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,variant_id%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "variant_id",beta_col = "lnOR",se_col = "standard_error",effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",pval_col = "p_value")
		result = mr_drug_binary(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD = rbind(Result_NAFLD,result)
	 }
   }
}

################ NAFLD Related Traits ################

exposure_analysis = c("GLP1R","GIPR")
outcome_file = c("Cirrhosis_Processed.txt","PDFF_Processed.txt","GGT.tsv.gz","ALT.tsv.gz","ALP.tsv.gz") #This is the file name of preprocessed outcome data
Result_NAFLD_Related = NULL
setwd("~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits") #This is the file path of preprocessed outcome data
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   outcome_origin = fread(outcome)
   for(j in 1:length(exposure_analysis)){
	 exposure_name = exposure_analysis[j]
	 exposure_dat = fread(paste0("~/GLP1/Analysis/Data/Exposure/Gene_expression/",exposure_analysis[j],".txt"))
	 if(outcome=="Cirrhosis_Processed.txt"){
	    outcome_name = strsplit(outcome,"_Processed.txt")[[1]][1]
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = mr_drug_binary(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD_Related = rbind(Result_NAFLD_Related,result)
	 }else if(outcome=="PDFF_Processed.txt"){
	    outcome_name = "PDFF"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Effect",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result =  mr_drug_continous(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD_Related = rbind(Result_NAFLD_Related,result)
	 }else if(outcome=="GGT.tsv.gz"|outcome=="ALT.tsv.gz"|outcome=="ALP.tsv.gz"){
	    outcome_name = strsplit(outcome,".tsv.gz")[[1]][1]
		print(outcome_name)
		outcome_match = subset(outcome_origin,SNP%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",
		                         pval_col = "P",eaf_col="A1FREQ")
		result =  mr_drug_continous(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_NAFLD_Related = rbind(Result_NAFLD_Related,result)
	 }
	 
   }
}


######################## Step 3 Drug Targe MR in Specific Tissues ########################

library(TwoSampleMR)
library(data.table)
library(mr.raps)

exposure_analysis = "GLP1R"
outcome_file = c("NAFLD_Discovery_Processed.txt","NAFLD_Validation.txt","Cirrhosis_Processed.txt","PDFF_Processed.txt","GGT.tsv.gz","ALT.tsv.gz","ALP.tsv.gz")
Result_Pancreas = NULL
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   if(length(grep("NAFLD",outcome))!=0){
     setwd("~/GLP1/Analysis/Data/Outcome/NAFLD")
	 outcome_origin = fread(outcome)
   }else{
     setwd("~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits")
	 outcome_origin = fread(outcome)
   }
   exposure_name = exposure_analysis
   exposure_dat = fread(paste0("~/GLP1/Analysis/Data/Exposure/Gtex/",exposure_analysis,".txt"))
   if(outcome=="NAFLD_Discovery_Processed.txt"){
	    outcome_name = "NAFLD_Discovery_Processed"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = mr_drug_binary(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_Pancreas = rbind(Result_Pancreas,result)
	}else if(outcome=="NAFLD_Validation.txt"){
	    outcome_name = "NAFLD_Validation"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,variant_id%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "variant_id",beta_col = "lnOR",se_col = "standard_error",effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",pval_col = "p_value")
		result = mr_drug_binary(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_Pancreas = rbind(Result_Pancreas,result)
	}else if(outcome=="Cirrhosis_Processed.txt"){
	    outcome_name = strsplit(outcome,"_Processed.txt")[[1]][1]
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = mr_drug_binary(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_Pancreas = rbind(Result_Pancreas,result)
	 }else if(outcome=="PDFF_Processed.txt"){
	    outcome_name = "PDFF"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Effect",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result =  mr_drug_continous(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_Pancreas = rbind(Result_Pancreas,result)
	 }else if(outcome=="GGT.tsv.gz"|outcome=="ALT.tsv.gz"|outcome=="ALP.tsv.gz"){
	    outcome_name = strsplit(outcome,".tsv.gz")[[1]][1]
		print(outcome_name)
		outcome_match = subset(outcome_origin,SNP%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",
		                         pval_col = "P",eaf_col="A1FREQ")
		result =  mr_drug_continous(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_Pancreas = rbind(Result_Pancreas,result)
	 }
}

######################## Step S1 Sensitive Analysis ########################

mr_senstive = function(data_exposure,exposure_name,data_outcome,outcome_name){
  #This function is used to perform sensitive analysis of Mendelian Randomization
  #parameter@data_exposure: Output from format_data()
  #parameter@exposure_name: Name of the exposure
  #parameter@data_outcome: Output from format_data()
  #parameter@outcome_name: Name of outcome
  #output: A dataframe including results for heterogeneity and pleiotropy test
  
  dat = harmonise_data(data_exposure,data_outcome,action=3)
  dat = subset(dat,mr_keep==TRUE)
  if(dim(dat)[1]>2){
    het = mr_heterogeneity(dat)
	het = subset(het, method=="Inverse variance weighted")
    plt = mr_pleiotropy_test(dat)
	sensitive_result = data.frame(t(c(exposure_name,outcome_name,het$Q,het$Q_pval,plt$egger_intercept,plt$se,plt$pval)))
	colnames(sensitive_result) = c("Exposure","Outcome","Q","Q_pval","Egger_intercept","SE","Intercept_pval")
	return(sensitive_result)
	
  }
}

######## Pleiotropy and Heterogeneity #######

#exposure_analysis = c("GLP1R","GIPR")
exposure_analysis = c("GLP1R","GIPR","Combine")
outcome_file = c("NAFLD_Discovery_Processed.txt","NAFLD_Validation.txt","Cirrhosis_Processed.txt","PDFF_Processed.txt","GGT.tsv.gz","ALT.tsv.gz","ALP.tsv.gz")
Result_Sensitive = NULL
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   if(length(grep("NAFLD",outcome))!=0){
     setwd("~/GLP1/Analysis/Data/Outcome/NAFLD")
	 outcome_origin = fread(outcome)
   }else{
     setwd("~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits")
	 outcome_origin = fread(outcome)
   }
   for(j in 1:length(exposure_analysis)){
     exposure_name = exposure_analysis[j]
     exposure_dat = fread(paste0("~/GLP1/Analysis/Data/Exposure/Downstream/HbA1c_",exposure_analysis[j],"_clump.txt"))
	 #exposure_dat = fread(paste0("~/GLP1/Analysis/Data/Exposure/Gene_expression/",exposure_analysis[j],".txt"))
	  if(outcome=="NAFLD_Discovery_Processed.txt"){
	    outcome_name = "NAFLD_Discovery_Processed"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = mr_senstive(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_Sensitive = rbind(Result_Sensitive,result)
	}else if(outcome=="NAFLD_Validation.txt"){
	    outcome_name = "NAFLD_Validation"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,variant_id%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "variant_id",beta_col = "lnOR",se_col = "standard_error",effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",pval_col = "p_value")
		result = mr_senstive(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_Sensitive = rbind(Result_Sensitive,result)
	}else if(outcome=="Cirrhosis_Processed.txt"){
	    outcome_name = strsplit(outcome,"_Processed.txt")[[1]][1]
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = mr_senstive(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_Sensitive = rbind(Result_Sensitive,result)
	 }else if(outcome=="PDFF_Processed.txt"){
	    outcome_name = "PDFF"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Effect",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result =  mr_senstive(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_Sensitive = rbind(Result_Sensitive,result)
	 }else if(outcome=="GGT.tsv.gz"|outcome=="ALT.tsv.gz"|outcome=="ALP.tsv.gz"){
	    outcome_name = strsplit(outcome,".tsv.gz")[[1]][1]
		print(outcome_name)
		outcome_match = subset(outcome_origin,SNP%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",
		                         pval_col = "P",eaf_col="A1FREQ")
		result =  mr_senstive(exposure_dat,exposure_name,outcome_dat,outcome_name)
	    Result_Sensitive = rbind(Result_Sensitive,result)
	 }
	 
   }
    
}


######################## Step S2 Positive Control for Gene Expression########################

library(TwoSampleMR)
library(data.table)

################ T2D ################

exposure = c("GLP1R","GIPR")
setwd("~/GLP1/Analysis/Data/Exposure/Gene_expression")
outcome_total = fread("~/GLP1/Analysis/Data/Outcome/Postitive_Control/T2D_GWAS.txt")
Result_T2D = NULL
for(i in 1:length(exposure)){
  exposure_dat = fread(paste0(exposure[i],".txt"))
  outcome_match = subset(outcome_total,rsID%in%exposure_dat$SNP)
  if(dim(outcome_match)[1]!=0){
    outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsID",beta_col = "Fixed-effects_beta",se_col = "Fixed-effects_SE",
	                          effect_allele_col = "effect_allele",other_allele_col = "other_allele",
							  pval_col = "Fixed-effects_p-value",eaf_col = "effect_allele_frequency")
    dat = harmonise_data(exposure_dat,outcome_dat,action=3)
    dat = subset(dat,mr_keep==TRUE)
    res = mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_weighted_median","mr_egger_regression"))
	or = generate_odds_ratios(res)
    het = mr_heterogeneity(dat)
    plt = mr_pleiotropy_test(dat)
	
	or$outcome = "T2D"
	or$exposure = exposure[i]
	if(subset(het,method=="Inverse variance weighted")$Q_pval<=0.05){
      print(paste0("Heterogeneity:",exposure[i]))
    }else if((dim(dat)[1]>2)&(plt$pval<=0.05)){
      print(paste0("Pleiotropy:",exposure[i]))
    }
	
	Result_T2D = rbind(Result_T2D,or)
  }else{
    print(paste0("No corresponding SNPs of ",exposure[i]," was found in GWAS"))
  }
}

################ BMI ################

exposure = c("GLP1R","GIPR")
setwd("~/GLP1/Analysis/Data/Exposure/Gene_expression")
outcome_total = fread("~/GLP1/Analysis/Data/Outcome/Postitive_Control/BMI_GWAS.txt")
Result_BMI = NULL
for(i in 1:length(exposure)){
  exposure_dat = fread(paste0(exposure[i],".txt"))
  outcome_match = subset(outcome_total,SNP%in%exposure_dat$SNP)
  if(dim(outcome_match)[1]!=0){
    outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "Tested_Allele",
                                     other_allele_col = "Other_Allele",pval_col = "P",eaf_col = "Freq_Tested_Allele_in_HRS")
    dat = harmonise_data(exposure_dat,outcome_dat,action=3)
    dat = subset(dat,mr_keep==TRUE)
    res = mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_weighted_median","mr_egger_regression"))
    het = mr_heterogeneity(dat)
    plt = mr_pleiotropy_test(dat)
	
	res$outcome = "BMI"
	res$exposure = exposure[i]
	res$lo_ci = res$b-1.96*res$se
    res$up_ci = res$b+1.96*res$se
	if(subset(het,method=="Inverse variance weighted")$Q_pval<=0.05){
      print(paste0("Heterogeneity:",exposure[i]))
    }else if((dim(dat)[1]>2)&(plt$pval<=0.05)){
      print(paste0("Pleiotropy:",exposure[i]))
    }
	
	Result_BMI = rbind(Result_BMI,res)
  }else{
    print(paste0("No corresponding SNPs of ",exposure[i]," was found in GWAS"))
  }

}

############ Find the valid proxy SNP based on LD Matrix from 1000 Genome ############

#All the vcf files were downloaded from https://gwas.mrcieu.ac.uk/

######### BMI Origin #########

library(vcfR)
BMIdat = read.vcfR("~/GLP1/BMI_sup/Data/ieu-b-40.vcf.gz")
head(BMIdat@meta,12)
str(BMIdat)

Fix = as.data.frame(BMIdat@fix[,(1:5)])
Gt = as.data.frame(BMIdat@gt[,2])
colnames(Gt)="GT"

summary_dat = strsplit(as.character(Gt$GT),split=":")
beta = as.numeric(sapply(summary_dat,"[",1))
se = as.numeric(sapply(summary_dat,"[",2))
pval = as.numeric(sapply(summary_dat,"[",3))
pval = 10^(-pval)
eaf = as.numeric(sapply(summary_dat,"[",4))
N = as.numeric(sapply(summary_dat,"[",5))
ID = sapply(summary_dat,"[",6)
mydat = data.frame(SNP = ID, beta = beta,se = se,pval = pval, maf = eaf, N = N)
colnames(mydat)=c("ID","Beta","Se","Pval","MAF","N")
data_all = merge(mydat,Fix)

write.table(data_all,"~/GLP1/BMI_sup/Data/BMI.summats.txt",quote=FALSE,row.names=FALSE,sep="\t")

######### GLP1R Origin #########

library(vcfR)
GLP1R = read.vcfR("~/GLP1/BMI_sup/Data/eqtl-a-ENSG00000112164.vcf.gz")
head(GLP1R@meta,12)
str(GLP1R)

Fix = as.data.frame(GLP1R@fix[,(1:5)])
Gt = as.data.frame(GLP1R@gt[,2])
colnames(Gt)="GT"

summary_dat = strsplit(as.character(Gt$GT),split=":")
beta = as.numeric(sapply(summary_dat,"[",1))
se = as.numeric(sapply(summary_dat,"[",2))
pval = as.numeric(sapply(summary_dat,"[",3))
pval = 10^(-pval)
eaf = as.numeric(sapply(summary_dat,"[",4))

N = as.numeric(sapply(summary_dat,"[",5))
ID = sapply(summary_dat,"[",6)
mydat = data.frame(SNP = ID, beta = beta,se = se,pval = pval, maf = eaf, N = N)
colnames(mydat)=c("ID","Beta","Se","Pval","MAF","N")
data_all = merge(mydat,Fix)
write.table(data_all,"~/GLP1/BMI_sup/Data/GLP1R.summats.txt",quote=FALSE,row.names=FALSE,sep="\t")

######### GIPR Origin #########

library(vcfR)
GIPR = read.vcfR("~/GLP1/BMI_sup/Data/eqtl-a-ENSG00000010310.vcf.gz")
head(GIPR@meta,12)
str(GIPR)

Fix = as.data.frame(GIPR@fix[,(1:5)])
Gt = as.data.frame(GIPR@gt[,2])
colnames(Gt)="GT"

summary_dat = strsplit(as.character(Gt$GT),split=":")
beta = as.numeric(sapply(summary_dat,"[",1))
se = as.numeric(sapply(summary_dat,"[",2))
pval = as.numeric(sapply(summary_dat,"[",3))
pval = 10^(-pval)
eaf = as.numeric(sapply(summary_dat,"[",4))

N = as.numeric(sapply(summary_dat,"[",5))
ID = sapply(summary_dat,"[",6)
mydat = data.frame(SNP = ID, beta = beta,se = se,pval = pval, maf = eaf, N = N)
colnames(mydat)=c("ID","Beta","Se","Pval","MAF","N")
data_all = merge(mydat,Fix)
write.table(data_all,"~/GLP1/BMI_sup/Data/GIPR.summats.txt",quote=FALSE,row.names=FALSE,sep="\t")

######### SNP List For GLP1R LDMatrix ######### 

SNP_all = fread("~/GLP1/BMI_sup/Annotation/snnplist.avinput") #snnplist.avinput is a file including the rsID, position, effect allele, other allele of SNPs
SNP_chr = subset(SNP_all,CHR=="chr6")
SNP_gene = subset(SNP_chr,START>=39016557-100000&STOP<=39059079+100000) #hg19
SNP_use = SNP_gene$rsID

GLP1R = read.table("~/GLP1/BMI_sup/Data/GLP1R.summats.txt",head=TRUE)
GLP1R_sig = subset(GLP1R,Pval<=5e-8)
SNP_GLP1R = GLP1R_sig$ID

SNP_for_LD = intersect(SNP_use,SNP_GLP1R)
write.table(SNP_for_LD,"~/GLP1/BMI_sup/Data/GLP1R_SNPlist.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

######### SNP List For GIPR LDMatrix ######### 

SNP_all = fread("~/GLP1/BMI_sup/Annotation/snnplist.avinput")
SNP_chr = subset(SNP_all,CHR=="chr19")
SNP_gene = subset(SNP_chr,START>=46171479-100000&STOP<=46186980+100000)#hg19
SNP_use = SNP_gene$rsID

GIPR = read.table("~/GLP1/BMI_sup/Data/GIPR.summats.txt",head=TRUE)
GIPR_sig = subset(GIPR,Pval<=5e-8)
SNP_GIPR = GIPR_sig$ID

SNP_for_LD = intersect(SNP_use,SNP_GIPR)
write.table(SNP_for_LD,"~/GLP1/BMI_sup/Data/GIPR_SNPlist.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

######### Calculate LDMatrix #########
#In this step, GIPR could be alternatived by GLP1R to obtain necessary information 

setwd("/data_200t/kuazuxue/data/1000genome")
library(BEDMatrix)
Geno = BEDMatrix("Mergeplk",simple=TRUE) # There are the bfiles of 1000 genome
SNPlist = read.table("~/GLP1/BMI_sup/Data/GIPR_SNPlist.txt",stringsAsFactor=FALSE)[,1]
Genotype = Geno[,SNPlist]

library('pheatmap')
correlation = matrix(data = NA, nrow = dim(Genotype)[2], ncol = dim(Genotype)[2])
for(i in 1:dim(Genotype)[2]){
  for(j in 1:dim(Genotype)[2]){
    correlation[i,j]= cor(Genotype[,i],Genotype[,j])
  }
}
correlation = data.frame(correlation)

name = SNPlist
colnames(correlation) = name
rownames(correlation) = name
png("~/GLP1/BMI_sup/Correlation_GIPR.png",width=21.6,height=10.2,units = "in",res = 720) 
pheatmap(as.matrix(correlation),cluster_row = FALSE, cluster_col = FALSE,
         display_numbers = 'TRUE', number_format = "%.2f", number_color="black")
dev.off()		 

######### Check if the proxy is in BMI GWAS #########

BMI_gwas = fread("~/GLP1/BMI_sup/Data/BMI.summats.txt")
BMI_snp = BMI_gwas$ID

#For GIPR, the proxy SNP is rs2075142 or rs2334255
"rs2075142"%in%BMI_snp #TRUE
"rs2334255"%in%BMI_snp #TRUE
#For GLP1R, the proxy SNP is rs9380825
"rs9380825"%in%BMI_snp #TRUE

######### Conduct New MR with proxy SNP #########  

library(TwoSampleMR)
exposure = read.table("~/GLP1/BMI_sup/Data/GLP1R.summats.txt",head=TRUE)
exposure = format_data(exposure,snp_col = "ID",beta_col = "Beta",se_col = "Se",pval_col = "Pval",effect_allele_col = "ALT",
                       other_allele_col = "REF",eaf_col = "MAF",samplesize_col="N")
exposure_proxy = subset(exposure,SNP=="rs9380825") #rs2075142 for GIPR /rs9380825 for GLP1R 
outcome_dat=read_outcome_data(snps = exposure_proxy$SNP,filename = "~/GLP1/BMI_sup/Data/BMI.summats.txt",sep = "\t",
                               snp_col = "ID",beta_col = "Beta",se_col = "Se",effect_allele_col = "ALT",other_allele_col = "REF",
                               pval_col = "Pval",eaf_col = "MAF")
dat=harmonise_data(exposure_proxy,outcome_dat)
result = mr(dat)
result









