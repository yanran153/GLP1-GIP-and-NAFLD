
################################ Part S2 Drug Target Analysis in Different Tissues ################################

#This part aims to explore specific causal effect of gene expressions in their highly expressed tissue and NAFLD as well as its clinicaly relevent traits
#The code for leave-one-out analysis could be found in Part_S1_leave_one_out.r

library(data.table)
library(TwoSampleMR)
library(gtools)

######################## Identifiy the targeted tissues ########################

setwd("~/GLP1/data/Exposure/Gtex/tissue")
file_all = list.files("~/GLP1/data/Exposure/Gtex/tissue")
file_tissue = file_all[grep(".v8.egenes.txt.gz",file_all)]
file_tissue = mixedsort(file_tissue)
gene = "GLP1R" ##GIPR

gene_tissue_list = NULL
for(i in 1:length(file_tissue)){
  data_tissue = fread(file_tissue[i])
  for(j in 1:length(gene)){
    data_gene = subset(data_tissue,gene_name==gene[j])
	if(dim(data_gene)[1]!=0){
	  if(data_gene$qval<=0.05){
	   tissue_name = strsplit(file_tissue[i],".v8.egenes.txt.gz")[[1]][1]
	   gene_tissue = c(gene[j],tissue_name)
	   gene_tissue_list = rbind(gene_tissue_list,gene_tissue)
	  }
	}
	
  }
}
colnames(gene_tissue_list) = c("Gene","Highly-expressed Tissue")


######################## Match cis-snp ########################

k1 = as.integer(args[1]) 
tissue_match = read.table("~/GLP1/data/Exposure/Gtex/SNP_Transfer/Gene_tissue_list.txt",head=TRUE,stringsAsFactor=FALSE)
tissue_name = tissue_match[k1,2]
gene = tissue_match[k1,1]

library(data.table)
library(TwoSampleMR)

setwd("~/GLP1/data/Exposure/Gtex/tissue")
data_tissue = fread(paste0("~/GLP1/data/Exposure/Gtex/tissue/",tissue_name,".v8.egenes.txt.gz"))
data_gene = subset(data_tissue,gene_name==gene)
gene_ens = subset(data_gene,gene_name==gene)$gene_id
gene_start = subset(data_gene,gene_name==gene)$gene_start 
gene_end = subset(data_gene,gene_name==gene)$gene_end

gene_tissue = fread(paste0("~/GLP1/data/Exposure/Gtex/tissue/",tissue_name,".v8.signif_variant_gene_pairs.txt.gz"))
gene_target = data.frame(subset(gene_tissue,gene_id==gene_ens))
gene_snp_information = strsplit(as.character(gene_target[,1]),"_")
gene_target$CHR = sapply(gene_snp_information,'[',1)
gene_target$POS = sapply(gene_snp_information,'[',2)
gene_target$NEA = sapply(gene_snp_information,'[',3)
gene_target$EA = sapply(gene_snp_information,'[',4)
gene_cis = subset(gene_target,POS>=gene_start-100000&POS<=gene_end+100000)
gene_cis_sig = subset(gene_cis,pval_nominal<=5e-8)

gene_for_annovar = gene_cis[,c("CHR","POS","POS","NEA","EA")]

#fwrite(gene_for_annovar,paste0("~/GLP1/data/Exposure/Gtex/SNP_Transfer/SNP_tmp/",gene,"_",tissue_name,"_cis.avinput"),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
#system(paste0("/home/yuanzhongshang/annovar/annotate_variation.pl ~/GLP1/data/Exposure/Gtex/SNP_Transfer/SNP_tmp/",gene,"_",tissue_name,"_cis.avinput ~/GLP1/BMI_sup/Annotation/ -filter -build hg38 -dbtype avsnp150"))

gene_annovar_match = fread(paste0("~/GLP1/data/Exposure/Gtex/SNP_Transfer/SNP_tmp/",gene,"_",tissue_name,"_cis.avinput.hg38_avsnp150_dropped"))
gene_for_gene = merge(gene_cis,gene_annovar_match,by.x=c("POS","NEA","EA"),by.y=c("V5","V6","V7"))
gene_dat_mr = format_data(gene_for_gene,snp_col = "V2",beta_col = "slope",se_col="slope_se",effect_allele_col = "EA",other_allele_col = "NEA",
		                pval_col = "pval_nominal",eaf_col="maf",pos_col="POS",chr_col="CHR")
gene_dat_clump = clump_data(gene_dat_mr,clump_kb = 10000,clump_r2 = 0.3)
write.table(gene_dat_clump,paste0("~/GLP1/data/Exposure/Gtex/SNP_Transfer/SNP_exposure/",gene,"_",tissue_name,".txt"),quote=FALSE,row.names=FALSE,sep="\t")


######################## Conduct Mendelian Randomization of NAFLD ########################

mr_drug_binary = function(data_exposure,exposure_name,data_outcome,outcome_name,tissue_name){
  #This function is used to perform Mendelian Randomization on binary outcome
  #parameter@data_exposure: Output from format_data()
  #parameter@exposure_name: Name of the exposure
  #parameter@data_outcome: Output from format_data()
  #parameter@outcome_name: Name of outcome
  #output: A dataframe including the result from IVW-R/weighted median/MR_RAPS  or from Wald Raios
  
  dat = harmonise_data(data_exposure,data_outcome,action=3)
  dat = subset(dat,mr_keep==TRUE)
if(dim(dat)[1]>=2){
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
  or$id.exposure = tissue_name
  return(or)																	
}else if(dim(dat)[1]==1){
  res = mr(dat)
  or = generate_odds_ratios(res)
  or$b = or$or
  or$lo_ci = or$or_lci95 
  or$up_ci = or$or_uci95 
  or = or[,-c(12,13,14)]
  or$outcome = outcome_name
  or$exposure = exposure_name
  or$id.exposure = tissue_name  
  return(or)
}
}

mr_drug_continous = function(data_exposure,exposure_name,data_outcome,outcome_name,tissue_name){
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
  res$id.exposure = tissue_name
  return(res)																	
}else if(dim(dat)[1]==1){
  res = mr(dat)
  res$outcome = outcome_name
  res$exposure = exposure_name 
  res$lo_ci = res$b-1.96*res$se
  res$up_ci = res$b+1.96*res$se
  res$id.exposure = tissue_name
  return(res)
}
}

exposure_analysis = list.files("~/GLP1/data/Exposure/Gtex/SNP_Transfer/SNP_exposure/clumped")
outcome_file = c("NAFLD_Discovery_Processed.txt","NAFLD_Validation.txt") #This is the file name of preprocessed outcome data
Result_NAFLD = NULL
setwd("~/GLP1/Analysis/Data/Outcome/NAFLD")
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   outcome_origin = fread(outcome)
   for(j in 1:length(exposure_analysis)){
	   exposure_gene = strsplit(exposure_analysis[j],"_")[[1]][1]
	   exposure_tissue = strsplit(strsplit(exposure_analysis[j],"R_")[[1]][2],".txt")[[1]][1]
	   exposure_dat = fread(paste0("~/GLP1/data/Exposure/Gtex/SNP_Transfer/SNP_exposure/clumped/",exposure_analysis[j]))
	 if(outcome=="NAFLD_Discovery_Processed.txt"){
	    outcome_name = "NAFLD_Discovery"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		if(dim(outcome_match)[1]>0){
		 outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		 result = mr_drug_binary(exposure_dat,exposure_gene,outcome_dat,outcome_name,exposure_tissue)
	     Result_NAFLD = rbind(Result_NAFLD,result)
		}
	 }else if(outcome=="NAFLD_Validation.txt"){
	    outcome_name = "NAFLD_Validation"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,variant_id%in%exposure_dat$SNP)
		if(dim(outcome_match)[1]>0){
		 outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "variant_id",beta_col = "lnOR",se_col = "standard_error",effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",pval_col = "p_value")
		 result = mr_drug_binary(exposure_dat,exposure_gene,outcome_dat,outcome_name,exposure_tissue)
	     Result_NAFLD = rbind(Result_NAFLD,result)
		}
	 }
   }
}


exposure_analysis = list.files("~/GLP1/data/Exposure/Gtex/SNP_Transfer/SNP_exposure/clumped")
outcome_file = c("Cirrhosis_Processed.txt","PDFF_Processed.txt","GGT.tsv.gz","ALT.tsv.gz","ALP.tsv.gz") #This is the file name of preprocessed outcome data
Result_NAFLD_Related = NULL
setwd("~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits") #This is the file path of preprocessed outcome data
for(i in 1:length(outcome_file)){
   outcome = outcome_file[i]
   outcome_origin = fread(outcome)
   for(j in 1:length(exposure_analysis)){
	   exposure_gene = strsplit(exposure_analysis[j],"_")[[1]][1]
	   exposure_tissue = strsplit(strsplit(exposure_analysis[j],"R_")[[1]][2],".txt")[[1]][1]
	   exposure_dat = fread(paste0("~/GLP1/data/Exposure/Gtex/SNP_Transfer/SNP_exposure/clumped/",exposure_analysis[j]))
	 if(outcome=="Cirrhosis_Processed.txt"){
	    outcome_name = strsplit(outcome,"_Processed.txt")[[1]][1]
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		if(dim(outcome_match)[1]>0){
		 outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		 result = mr_drug_binary(exposure_dat,exposure_gene,outcome_dat,outcome_name,exposure_tissue)
	     Result_NAFLD_Related = rbind(Result_NAFLD_Related,result)
		}
	 }else if(outcome=="PDFF_Processed.txt"){
	    outcome_name = "PDFF"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		if(dim(outcome_match)[1]>0){
		 outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Effect",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		 result =  mr_drug_continous(exposure_dat,exposure_gene,outcome_dat,outcome_name,exposure_tissue)
	     Result_NAFLD_Related = rbind(Result_NAFLD_Related,result)
		}
	 }else if(outcome=="GGT.tsv.gz"|outcome=="ALT.tsv.gz"|outcome=="ALP.tsv.gz"){
	    outcome_name = strsplit(outcome,".tsv.gz")[[1]][1]
		print(outcome_name)
		outcome_match = subset(outcome_origin,SNP%in%exposure_dat$SNP)
		if(dim(outcome_match)[1]>0){
		 outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",
		                         pval_col = "P",eaf_col="A1FREQ")
		 result =  mr_drug_continous(exposure_dat,exposure_gene,outcome_dat,outcome_name,exposure_tissue)
	     Result_NAFLD_Related = rbind(Result_NAFLD_Related,result)
		}
	 }
	 
   }
}


#################Sensitive Analysis ###########################
mr_senstive = function(data_exposure,exposure_name,data_outcome,outcome_name,tissue_name){
  #This function is used to perform sensitive analysis of Mendelian Randomization
  #parameter@data_exposure: Output from format_data()
  #parameter@exposure_name: Name of the exposure
  #parameter@data_outcome: Output from format_data()
  #parameter@outcome_name: Name of outcome
  #output: A dataframe including the result from IVW-R/weighted median/MR_RAPS  or from Wald Raios
  
  dat = harmonise_data(data_exposure,data_outcome,action=3)
  dat = subset(dat,mr_keep==TRUE)
  if(dim(dat)[1]>2){
    het = mr_heterogeneity(dat)
	het = subset(het, method=="Inverse variance weighted")
    plt = mr_pleiotropy_test(dat)
	sensitive_result = c(exposure_name,outcome_name,dim(dat)[1],tissue_name,het$Q,het$Q_pval,plt$egger_intercept,plt$se,plt$pval)
	return(sensitive_result)
  }
}

######## Pleiotropy and Heterogeneity #######

exposure_analysis = list.files("~/GLP1/data/Exposure/Gtex/SNP_Transfer/SNP_exposure/clumped")
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
       exposure_gene = strsplit(exposure_analysis[j],"_")[[1]][1]
	   exposure_tissue = strsplit(strsplit(exposure_analysis[j],"R_")[[1]][2],".txt")[[1]][1]
	   exposure_dat = fread(paste0("~/GLP1/data/Exposure/Gtex/SNP_Transfer/SNP_exposure/clumped/",exposure_analysis[j]))
	  if(outcome=="NAFLD_Discovery_Processed.txt"){
	    outcome_name = "NAFLD_Discovery"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		if(dim(outcome_match)[1]>0){
		 outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		 result = mr_senstive(exposure_dat,exposure_gene,outcome_dat,outcome_name,exposure_tissue)
	     Result_Sensitive = rbind(Result_Sensitive,result)
		}
	}else if(outcome=="NAFLD_Validation.txt"){
	    outcome_name = "NAFLD_Validation"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,variant_id%in%exposure_dat$SNP)
		if(dim(outcome_match)[1]>0){
		 outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "variant_id",beta_col = "lnOR",se_col = "standard_error",effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",pval_col = "p_value")
		 result = mr_senstive(exposure_dat,exposure_gene,outcome_dat,outcome_name,exposure_tissue)
	     Result_Sensitive = rbind(Result_Sensitive,result)
		}
	}else if(outcome=="Cirrhosis_Processed.txt"){
	    outcome_name = strsplit(outcome,"_Processed.txt")[[1]][1]
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		if(dim(outcome_match)[1]>0){
		 outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Beta",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		 result = mr_senstive(exposure_dat,exposure_gene,outcome_dat,outcome_name,exposure_tissue)
	     Result_Sensitive = rbind(Result_Sensitive,result)
		}
	 }else if(outcome=="PDFF_Processed.txt"){
	    outcome_name = "PDFF"
		print(outcome_name)
	    outcome_match = subset(outcome_origin,rsName%in%exposure_dat$SNP)
		if(dim(outcome_match)[1]>0){
		 outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "rsName",beta_col = "Effect",se_col="Se",effect_allele_col = "Amin",other_allele_col = "Amaj",
		                         pval_col = "Pval",eaf_col="MAF_PC")
		 result =  mr_senstive(exposure_dat,exposure_gene,outcome_dat,outcome_name,exposure_tissue)
	     Result_Sensitive = rbind(Result_Sensitive,result)
		}
	 }else if(outcome=="GGT.tsv.gz"|outcome=="ALT.tsv.gz"|outcome=="ALP.tsv.gz"){
	    outcome_name = strsplit(outcome,".tsv.gz")[[1]][1]
		print(outcome_name)
		outcome_match = subset(outcome_origin,SNP%in%exposure_dat$SNP)
		if(dim(outcome_match)[1]>0){
		 outcome_dat = format_data(outcome_match,type = "outcome",snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",
		                         pval_col = "P",eaf_col="A1FREQ")
		 result =  mr_senstive(exposure_dat,exposure_gene,outcome_dat,outcome_name,exposure_tissue)
	     Result_Sensitive = rbind(Result_Sensitive,result)
		}
	 }
	 
   }
    
}
colnames(Result_Sensitive) = c("Exposure","Outcome","SNPN","Tissue","Q","Q_pval","Egger_intercept","SE","Intercept_pval")
