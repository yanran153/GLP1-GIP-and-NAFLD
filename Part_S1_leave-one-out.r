
################################ Part S1 Leave-one-out Analysis ################################

#This part aims to conduct leave-one-out analysis for all the mendeiian randomization

library(data.table)
library(TwoSampleMR)

######################## Define the Leave-one-out function ########################

leave_one_out = function(data_exposure,exposure_name,data_outcome,outcome_name,binary){
   dat = harmonise_data(data_exposure,data_outcome,action=3)
   dat = subset(dat,mr_keep==TRUE)
   IV_number = dim(dat)
   if(IV_number[1]>2){
     if(binary=="TRUE"){
        LOO_inner = NULL
		for(i in 1:IV_number[1]){
          dat_inner = dat[-i,]
		  IV_name = as.character(dat[i,"SNP"])
		  res_inner = mr(dat_inner,method_list=c("mr_ivw_mre"))
		  res_or = generate_odds_ratios(res_inner)
          res_output = res_or[,c("or","se","pval","or_lci95","or_uci95")]
          res_output$IV = IV_name 
          LOO_inner = rbind(LOO_inner,res_output)
        }
		res_all = mr(dat,method_list=c("mr_ivw_mre"))
		res_all_or = generate_odds_ratios(res_all)
		res_output_all = res_all_or[,c("or","se","pval","or_lci95","or_uci95")]
        res_output_all$IV = 'All' 
        LOO_inner = rbind(LOO_inner,res_output_all)
		LOO_inner$exposure = exposure_name
		LOO_inner$outcome = outcome_name
		LOO = LOO_inner[,c(7,8,6,1,2,3,4,5)]
        colnames(LOO) = c("Exposure","Outcome","SNP","Estimation","SE","P_value","CI_low","CI_up")
      }else{
        LOO_inner = NULL
		for(i in 1:IV_number[1]){
          dat_inner = dat[-i,]
		  IV_name = as.character(dat[i,"SNP"])
		  res_inner = mr(dat_inner,method_list=c("mr_ivw_mre"))
	      res_inner$CI_low = res_inner$b - 1.96*res_inner$se
		  res_inner$CI_up = res_inner$b + 1.96*res_inner$se
          res_output = res_inner[,c("b","se","pval","CI_low","CI_up")]
          res_output$IV = IV_name 
          LOO_inner = rbind(LOO_inner,res_output)
        }
		res_all = mr(dat,method_list=c("mr_ivw_mre"))
        res_all$CI_low = res_all$b - 1.96*res_all$se
		res_all$CI_up = res_all$b + 1.96*res_all$se
        res_output_all = res_all[,c("b","se","pval","CI_low","CI_up")]
        res_output_all$IV = 'All' 
        LOO_inner = rbind(LOO_inner,res_output_all)
		LOO_inner$exposure = exposure_name
		LOO_inner$outcome = outcome_name
		LOO = LOO_inner[,c(7,8,6,1,2,3,4,5)]
        colnames(LOO) = c("Exposure","Outcome","SNP","Estimation","SE","P_value","CI_low","CI_up")
      }
	  return(LOO)
   } 
 }
 
 
######################## Leave-one-out Analysis for Step 1 Traditional MR Analysis ########################

exposure_analysis = c("GLP1_2h","GIP_2h","GIP_fasting")
outcome_file = c("NAFLD_Discovery_Processed.txt","NAFLD_Validation.txt","Cirrhosis_Processed.txt","PDFF_Processed.txt","GGT.tsv.gz","ALT.tsv.gz","ALP.tsv.gz")
Result_Loo = NULL
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
		result = leave_one_out(exposure_dat,exposure_name,outcome_dat,outcome_name,TRUE)
	    Result_Loo = rbind(Result_Loo,result)
	 }else if(outcome=="NAFLD_Validation.txt"){
	    outcome_name = strsplit(outcome,"_Processed")[[1]][1]
	    outcome_match = subset(outcome_origin,variant_id%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "variant_id",beta_col = "lnOR",se_col = "standard_error",effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",pval_col = "p_value")
		result =  leave_one_out(exposure_dat,exposure_name,outcome_dat,outcome_name,TRUE)
	    Result_Loo = rbind(Result_Loo,result) 
	 }else if(outcome=="GGT.tsv.gz"|outcome=="ALT.tsv.gz"|outcome=="ALP.tsv.gz"){
	    outcome_name = strsplit(outcome,".tsv.gz")[[1]][1]
		outcome_match = subset(outcome_origin,SNP%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",
		                         pval_col = "P",eaf_col="A1FREQ")
		result = leave_one_out(exposure_dat,exposure_name,outcome_dat,outcome_name,FALSE)
	    Result_Loo = rbind(Result_Loo,result)
	 }else if (outcome == "PDFF_Processed.txt"){
	    outcome_name = "PDFF"
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
	    outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Effect",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = leave_one_out(exposure_dat,exposure_name,outcome_dat,outcome_name,FALSE)
		Result_Loo = rbind(Result_Loo,result)
	 }
	 
   }
}
write.csv(Result_Loo,"~/GLP1/Analysis/Traditional_MR/Result_LOO_IVWR.csv")


######################## Leave-one-out Analysis for Step 2 Drug Target MR Analysis with gene expressions ########################

exposure_analysis = c("GLP1R","GIPR")
outcome_file = c("NAFLD_Discovery_Processed.txt","NAFLD_Validation.txt","Cirrhosis_Processed.txt","PDFF_Processed.txt","GGT.tsv.gz","ALT.tsv.gz","ALP.tsv.gz")
Result_Loo = NULL
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   if(outcome=="NAFLD_Discovery_Processed.txt"|outcome=="NAFLD_Validation.txt"){
    setwd("~/GLP1/Analysis/Data/Outcome/NAFLD")
   }else{
    setwd("~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits")
   }
   outcome_origin = fread(outcome)
   for(j in 1:length(exposure_analysis)){
     exposure_name = exposure_analysis[j]
     exposure_dat = fread(paste0("~/GLP1/Analysis/Data/Exposure/Gene_expression/",exposure_analysis[j],".txt"))
	 if(outcome=="NAFLD_Discovery_Processed.txt"|outcome=="Cirrhosis_Processed.txt"){
	    outcome_name = strsplit(outcome,"_Processed")[[1]][1]
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = leave_one_out(exposure_dat,exposure_name,outcome_dat,outcome_name,TRUE)
	    Result_Loo = rbind(Result_Loo,result)
	 }else if(outcome=="NAFLD_Validation.txt"){
	    outcome_name = strsplit(outcome,"_Processed")[[1]][1]
	    outcome_match = subset(outcome_origin,variant_id%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "variant_id",beta_col = "lnOR",se_col = "standard_error",effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",pval_col = "p_value")
		result =  leave_one_out(exposure_dat,exposure_name,outcome_dat,outcome_name,TRUE)
	    Result_Loo = rbind(Result_Loo,result) 
	 }else if(outcome=="GGT.tsv.gz"|outcome=="ALT.tsv.gz"|outcome=="ALP.tsv.gz"){
	    outcome_name = strsplit(outcome,".tsv.gz")[[1]][1]
		outcome_match = subset(outcome_origin,SNP%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",
		                         pval_col = "P",eaf_col="A1FREQ")
		result = leave_one_out(exposure_dat,exposure_name,outcome_dat,outcome_name,FALSE)
	    Result_Loo = rbind(Result_Loo,result)
	 }else if (outcome == "PDFF_Processed.txt"){
	    outcome_name = "PDFF"
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
	    outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Effect",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = leave_one_out(exposure_dat,exposure_name,outcome_dat,outcome_name,FALSE)
		Result_Loo = rbind(Result_Loo,result)
	 }
	 
   }
}
write.csv(Result_Loo,"~/GLP1/Analysis/Drug_Target_MR/Result_LOO_IVWR.csv")


######################## Leave-one-out Analysis for Step 2 Drug Target MR Analysis with downstream biomarkers ########################

exposure_analysis = c("GLP1R","GIPR","Combine")
outcome_file = c("NAFLD_Discovery_Processed.txt","NAFLD_Validation.txt","Cirrhosis_Processed.txt","PDFF_Processed.txt","GGT.tsv.gz","ALT.tsv.gz","ALP.tsv.gz")
Result_Loo = NULL
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   if(outcome=="NAFLD_Discovery_Processed.txt"|outcome=="NAFLD_Validation.txt"){
    setwd("~/GLP1/Analysis/Data/Outcome/NAFLD")
   }else{
    setwd("~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits")
   }
   outcome_origin = fread(outcome)
   for(j in 1:length(exposure_analysis)){
     exposure_name = exposure_analysis[j]
     exposure_dat = fread(paste0("~/GLP1/Analysis/Data/Exposure/Downstream/HbA1c_",exposure_analysis[j],"_clump.txt"))
	 if(outcome=="NAFLD_Discovery_Processed.txt"|outcome=="Cirrhosis_Processed.txt"){
	    outcome_name = strsplit(outcome,"_Processed")[[1]][1]
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = leave_one_out(exposure_dat,exposure_name,outcome_dat,outcome_name,TRUE)
	    Result_Loo = rbind(Result_Loo,result)
	 }else if(outcome=="NAFLD_Validation.txt"){
	    outcome_name = strsplit(outcome,"_Processed")[[1]][1]
	    outcome_match = subset(outcome_origin,variant_id%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "variant_id",beta_col = "lnOR",se_col = "standard_error",effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",pval_col = "p_value")
		result =  leave_one_out(exposure_dat,exposure_name,outcome_dat,outcome_name,TRUE)
	    Result_Loo = rbind(Result_Loo,result) 
	 }else if(outcome=="GGT.tsv.gz"|outcome=="ALT.tsv.gz"|outcome=="ALP.tsv.gz"){
	    outcome_name = strsplit(outcome,".tsv.gz")[[1]][1]
		outcome_match = subset(outcome_origin,SNP%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",
		                         pval_col = "P",eaf_col="A1FREQ")
		result = leave_one_out(exposure_dat,exposure_name,outcome_dat,outcome_name,FALSE)
	    Result_Loo = rbind(Result_Loo,result)
	 }else if (outcome == "PDFF_Processed.txt"){
	    outcome_name = "PDFF"
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
	    outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Effect",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = leave_one_out(exposure_dat,exposure_name,outcome_dat,outcome_name,FALSE)
		Result_Loo = rbind(Result_Loo,result)
	 }
	 
   }
}
write.csv(Result_Loo,"~/GLP1/Analysis/Drug_Target_MR/Result_LOO_Biomarker_IVWR.csv")


######################## Leave-one-out Analysis for Step 2 Drug Target MR Analysis with Highly Expressed Tissues ########################

leave_one_out = function(data_exposure,exposure_name,data_outcome,outcome_name,binary,tissue_name){
   dat = harmonise_data(data_exposure,data_outcome,action=3)
   dat = subset(dat,mr_keep==TRUE)
   IV_number = dim(dat)
   if(IV_number[1]>2){
     if(binary=="TRUE"){
        LOO_inner = NULL
		for(i in 1:IV_number[1]){
          dat_inner = dat[-i,]
		  IV_name = as.character(dat[i,"SNP"])
		  res_inner = mr(dat_inner,method_list=c("mr_ivw_mre"))
		  res_or = generate_odds_ratios(res_inner)
          res_output = res_or[,c("or","se","pval","or_lci95","or_uci95")]
          res_output$IV = IV_name 
          LOO_inner = rbind(LOO_inner,res_output)
        }
		res_all = mr(dat,method_list=c("mr_ivw_mre"))
		res_all_or = generate_odds_ratios(res_all)
		res_output_all = res_all_or[,c("or","se","pval","or_lci95","or_uci95")]
        res_output_all$IV = 'All' 
        LOO_inner = rbind(LOO_inner,res_output_all)
		LOO_inner$exposure = exposure_name
		LOO_inner$outcome = outcome_name
		LOO_inner$tissue = tissue_name
		LOO = LOO_inner[,c(7,8,9,6,1,2,3,4,5)]
        colnames(LOO) = c("Exposure","Outcome","Tissue","SNP","Estimation","SE","P_value","CI_low","CI_up")
      }else{
        LOO_inner = NULL
		for(i in 1:IV_number[1]){
          dat_inner = dat[-i,]
		  IV_name = as.character(dat[i,"SNP"])
		  res_inner = mr(dat_inner,method_list=c("mr_ivw_mre"))
	      res_inner$CI_low = res_inner$b - 1.96*res_inner$se
		  res_inner$CI_up = res_inner$b + 1.96*res_inner$se
          res_output = res_inner[,c("b","se","pval","CI_low","CI_up")]
          res_output$IV = IV_name 
          LOO_inner = rbind(LOO_inner,res_output)
        }
		res_all = mr(dat,method_list=c("mr_ivw_mre"))
        res_all$CI_low = res_all$b - 1.96*res_all$se
		res_all$CI_up = res_all$b + 1.96*res_all$se
        res_output_all = res_all[,c("b","se","pval","CI_low","CI_up")]
        res_output_all$IV = 'All' 
        LOO_inner = rbind(LOO_inner,res_output_all)
		LOO_inner$exposure = exposure_name
		LOO_inner$outcome = outcome_name
		LOO_inner$tissue = tissue_name
		LOO = LOO_inner[,c(7,8,9,6,1,2,3,4,5)]
        colnames(LOO) = c("Exposure","Outcome","Tissue","SNP","Estimation","SE","P_value","CI_low","CI_up")
      }
	  return(LOO)
   } 
 }
 
exposure_analysis = list.files("~/GLP1/data/Exposure/Gtex/SNP_Transfer/SNP_exposure/clumped")
outcome_file = c("NAFLD_Discovery_Processed.txt","NAFLD_Validation.txt","Cirrhosis_Processed.txt","PDFF_Processed.txt","GGT.tsv.gz","ALT.tsv.gz","ALP.tsv.gz")
Result_Loo = NULL
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   if(outcome=="NAFLD_Discovery_Processed.txt"|outcome=="NAFLD_Validation.txt"){
    setwd("~/GLP1/Analysis/Data/Outcome/NAFLD")
   }else{
    setwd("~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits")
   }
   outcome_origin = fread(outcome)
   print(outcome)
   for(j in 1:length(exposure_analysis)){
       exposure_gene = strsplit(exposure_analysis[j],"_")[[1]][1]
	   exposure_tissue = strsplit(strsplit(exposure_analysis[j],"R_")[[1]][2],".txt")[[1]][1]
	   exposure_dat = fread(paste0("~/GLP1/data/Exposure/Gtex/SNP_Transfer/SNP_exposure/clumped/",exposure_analysis[j]))
	 if(outcome=="NAFLD_Discovery_Processed.txt"|outcome=="Cirrhosis_Processed.txt"){
	    outcome_name = strsplit(outcome,"_Processed")[[1]][1]
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = leave_one_out(exposure_dat,exposure_gene,outcome_dat,outcome_name,TRUE,exposure_tissue)
	    Result_Loo = rbind(Result_Loo,result)
	 }else if(outcome=="NAFLD_Validation.txt"& exposure_tissue!="Whole_Blood"){
	    outcome_name = strsplit(outcome,"_Processed")[[1]][1]
	    outcome_match = subset(outcome_origin,variant_id%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "variant_id",beta_col = "lnOR",se_col = "standard_error",effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",pval_col = "p_value")
		result =  leave_one_out(exposure_dat,exposure_gene,outcome_dat,outcome_name,TRUE,exposure_tissue)
	    Result_Loo = rbind(Result_Loo,result) 
	 }else if(outcome=="GGT.tsv.gz"|outcome=="ALT.tsv.gz"|outcome=="ALP.tsv.gz"){
	    outcome_name = strsplit(outcome,".tsv.gz")[[1]][1]
		outcome_match = subset(outcome_origin,SNP%in%exposure_dat$SNP)
		outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",
		                         pval_col = "P",eaf_col="A1FREQ")
		result = leave_one_out(exposure_dat,exposure_gene,outcome_dat,outcome_name,FALSE,exposure_tissue)
	    Result_Loo = rbind(Result_Loo,result)
	 }else if (outcome == "PDFF_Processed.txt"){
	    outcome_name = "PDFF"
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
	    outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Effect",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		result = leave_one_out(exposure_dat,exposure_gene,outcome_dat,outcome_name,FALSE,exposure_tissue)
		Result_Loo = rbind(Result_Loo,result)
	 }
	 
   }
}
write.csv(Result_Loo,"~/GLP1/Analysis/Drug_Target_MR/Result_LOO_Tissue_IVWR.csv")