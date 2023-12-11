
################################ Part 3 Mediation Mendelian Randomization ################################

#This part aims to explore the mediating roles of Type 2 Diabetes and Body mass index on the association between gene expression and NAFLD
#We just calculate the indirect effect of gene expressions on NAFLD and the effects of T2D/BMI on NALFD because the others have calcuated in the previous part 
#Mediated Effect is calculated on https://amplab.shinyapps.io/MEDCI/ 
#The code for mediation analysis with replication dataset could be found in Part_S3_Mediation_Analysis_Replication.r

library(data.table)
library(TwoSampleMR)
library(mr.raps)

######################## Step 1 T2D-->NAFLD ########################

exp<-fread("~/GLP1/Analysis/Data/Outcome/Postitive_Control/T2D_GWAS.txt")
colnames(exp) = c("CHR","POS","Info","SNP","EA","OA","EAF","BETA","SE","P")
exp$P = as.numeric(exp$P)
exp = exp[exp$EAF>0.01 & exp$EAF<0.99,]
exposure <-subset(exp,exp$P<5e-08)

exposure = format_data(exposure,snp_col = "SNP",beta_col = "BETA",se_col = "SE",pval_col = "P",
                       effect_allele_col = "EA",other_allele_col = "OA", eaf_col = "EAF",
                       chr_col = "CHR",pos_col = "POS")
exposure_dat<-clump_data(exposure,clump_kb = 10000,clump_r2 = 0.01,clump_p1 = 1,clump_p2 = 1, pop = "EUR")
exposure_dat = fread("~/GLP1/Analysis/Data/Outcome/Postitive_Control/T2D_clump.txt")
outcome_dat<-read_outcome_data(snps = exposure_dat$SNP,filename = "~/GLP1/Analysis/Data/Outcome/NAFLD/NAFLD_Discovery_Processed.txt",sep = "\t",
                               snp_col = "rsName",beta_col = "Beta",se_col = "Se",effect_allele_col = "Amin",other_allele_col = "Amaj",pval_col = "Pval",
							   eaf_col="MAF_PC")

dat<-harmonise_data(exposure_dat,outcome_dat,action=3)
dat = subset(dat,mr_keep==TRUE)
res = mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_weighted_median","mr_egger_regression"))
or<-generate_odds_ratios(res)

mr_frame = dat[,c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome")]
names(mr_frame)[c(1:5)] = c("rsid","beta.x","se.x","beta.y","se.y")
MR_RAPS  =  mr.raps.overdispersed.robust(mr_frame$beta.x, mr_frame$beta.y, mr_frame$se.x, mr_frame$se.y,loss.function = "huber", k = 1.345, 
            initialization = c("l2"), suppress.warning = FALSE, diagnosis = FALSE, niter = 20, tol = .Machine$double.eps^0.5)
or$b = or$or
or$lo_ci = or$or_lci95 
or$up_ci = or$or_uci95 
or = or[,-c(12,13,14)]
or[4,]= NA
or[4,c(1,2,3,4,6)]= or[2,c(1,2,3,4,6)]
or[,"method"] = as.character(or[,"method"])
or[4,"method"]="MR-RAPS"
or[4,c("b","se","pval","lo_ci","up_ci")]=c(exp(MR_RAPS$beta.hat),MR_RAPS$beta.se,MR_RAPS$beta.p.value,
   										   exp(MR_RAPS$beta.hat-1.96*MR_RAPS$beta.se),exp(MR_RAPS$beta.hat+1.96*MR_RAPS$beta.se))
or$outcome = "NAFLD"
or$exposure = "T2D" 

het = mr_heterogeneity(dat)
plt = mr_pleiotropy_test(dat)


######################## Step 2 BMI-->NAFLD ########################

exp<-fread("~/GLP1/Analysis/Data/Outcome/Postitive_Control/BMI_GWAS.txt")
colnames(exp) = c("CHR","POS","SNP","EA","OA","EAF","BETA","SE","P","N")
exp = exp[exp$EAF>0.01 &exp$EAF<0.99,]
exposure<-subset(exp,exp$P<5e-08)

exposure = format_data(exposure,snp_col = "SNP",beta_col = "BETA",se_col = "SE",pval_col = "P",
                       effect_allele_col = "EA",other_allele_col = "OA", eaf_col = "EAF",
                       chr_col = "CHR",pos_col = "POS")
exposure_dat<-clump_data(exposure,clump_kb = 10000,clump_r2 = 0.01,clump_p1 = 1,clump_p2 = 1,pop = "EUR")
exposure_dat = fread("~/GLP1/Analysis/Data/Outcome/Postitive_Control/BMI_clump.txt")
outcome_dat<-read_outcome_data(snps = exposure_dat$SNP,filename = "~/GLP1/Analysis/Data/Outcome/NAFLD/NAFLD_Discovery_Processed.txt",sep = "\t",
                               snp_col = "rsName",beta_col = "Beta",se_col = "Se",effect_allele_col = "Amin",other_allele_col = "Amaj",pval_col = "P",
							   eaf_col="MAF_PC")
							   
dat<-harmonise_data(exposure_dat,outcome_dat,action=3)
dat = subset(dat,mr_keep==TRUE)
res = mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_weighted_median"))
or<-generate_odds_ratios(res)

mr_frame = dat[,c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome")]
names(mr_frame)[c(1:5)] = c("rsid","beta.x","se.x","beta.y","se.y")
MR_RAPS  =  mr.raps.overdispersed.robust(mr_frame$beta.x, mr_frame$beta.y, mr_frame$se.x, mr_frame$se.y,loss.function = "huber", k = 1.345, 
            initialization = c("l2"), suppress.warning = FALSE, diagnosis = FALSE, niter = 20, tol = .Machine$double.eps^0.5)
or$b = or$or
or$lo_ci = or$or_lci95 
or$up_ci = or$or_uci95 
or = or[,-c(12,13,14)]
or[4,]= NA
or[4,c(1,2,3,4,6)]= or[2,c(1,2,3,4,6)]
or[,"method"] = as.character(or[,"method"])
or[4,"method"]="MR-RAPS"
or[4,c("b","se","pval","lo_ci","up_ci")]=c(exp(MR_RAPS$beta.hat),MR_RAPS$beta.se,MR_RAPS$beta.p.value,
   										   exp(MR_RAPS$beta.hat-1.96*MR_RAPS$beta.se),exp(MR_RAPS$beta.hat+1.96*MR_RAPS$beta.se))
or$outcome = "NAFLD"
or$exposure = "BMI" 

het = mr_heterogeneity(dat)
plt = mr_pleiotropy_test(dat)

######################## Step 3 Direct Effect ########################

library(TwoSampleMR)
library(data.table)

######### GLP1R->T2D->NAFLD #########
exposure_dat1_origin = fread("~/GLP1/BMI_sup/Data/GLP1R.summats.txt")
exposure_dat2_origin = fread("~/GLP1/Analysis/Data/Outcome/Postitive_Control/T2D_GWAS.txt")
colnames(exposure_dat2_origin) = c("CHR","POS","Info","SNP","EA","OA","EAF","BETA","SE","P")
exposure_dat1_origin$exposure = "GLP1R"
exposure_dat2_origin$exposure = "T2D"
exposure_dat2_origin$P = as.numeric(exposure_dat2_origin$P)
IV_candidate = intersect(exposure_dat1_origin$ID,exposure_dat2_origin$SNP)

exposure_dat1_mvmr = subset(exposure_dat1_origin,ID%in%IV_candidate)
exposure_dat2_mvmr = subset(exposure_dat2_origin,SNP%in%IV_candidate)


exposure_dat1_mvmr = format_data(exposure_dat1_mvmr,snp_col = "ID",beta_col = "Beta",se_col = "Se",pval_col = "Pval",effect_allele_col = "ALT",
                       other_allele_col = "REF",eaf_col = "MAF")
exposure_dat1_mvmr$exposure = "GLP1R"
exposure_dat2_mvmr = format_data(exposure_dat2_mvmr,snp_col = "SNP",beta_col = "BETA",se_col = "SE",pval_col = "P",effect_allele_col = "EA",
                       other_allele_col = "OA",eaf_col = "EAF")
exposure_dat2_mvmr$exposure = "T2D"
exposure_dat2_mvmr = exposure_dat2_mvmr[,colnames(exposure_dat1_mvmr)]

exposure_mvmr = rbind(exposure_dat1_mvmr,exposure_dat2_mvmr)
exposure_mvmr_end = NULL
for(i in 1:length(IV_candidate)){
   exposure_judge = subset(exposure_mvmr,SNP==IV_candidate[i])
   if(exposure_judge$pval.exposure[1]<=5e-8|exposure_judge$pval.exposure[2]<=5e-8){
    exposure_mvmr_end = rbind(exposure_mvmr_end,exposure_judge)   
   }
}
exposure_mvmr_clump = clump_data(exposure_mvmr_end,clump_kb = 10000,clump_r2 = 0.01)
exposure_mvmr_dat1 = subset(exposure_mvmr_clump,exposure=="GLP1R")
exposure_mvmr_dat2 = subset(exposure_mvmr_clump,exposure=="T2D")
mvmr_snp = union(exposure_mvmr_dat1$SNP,exposure_mvmr_dat2$SNP)
exposure_mvmr_dat = subset(xposure_mvmr_end,SNP%in%mvmr_snp)
outcome_dat = read_outcome_data(snps = exposure_mvmr_dat$SNP,filename = "~/GLP1/Analysis/Data/Outcome/NAFLD/NAFLD_Discovery_Processed.txt",sep = "\t",
                               snp_col = "rsName",beta_col = "Beta",se_col = "Se",effect_allele_col = "Amin",other_allele_col = "Amaj",pval_col = "Pval",
							   eaf_col="MAF_PC")
mvdat = mv_harmonise_data(exposure_mvmr_clump,outcome_dat) 
res = mv_multiple(mvdat) 
res


######### GIPR->T2D->NAFLD #########
exposure_dat1_origin = fread("~/GLP1/BMI_sup/Data/GIPR.summats.txt")
exposure_dat2_origin = fread("~/GLP1/Analysis/Data/Outcome/Postitive_Control/T2D_GWAS.txt")
colnames(exposure_dat2_origin) = c("CHR","POS","Info","SNP","EA","OA","EAF","BETA","SE","P")
exposure_dat1_origin$exposure = "GIPR"
exposure_dat2_origin$exposure = "T2D"
exposure_dat2_origin$P = as.numeric(exposure_dat2_origin$P)
IV_candidate = intersect(exposure_dat1_origin$ID,exposure_dat2_origin$SNP)

exposure_dat1_mvmr = subset(exposure_dat1_origin,ID%in%IV_candidate)
exposure_dat2_mvmr = subset(exposure_dat2_origin,SNP%in%IV_candidate)


exposure_dat1_mvmr = format_data(exposure_dat1_mvmr,snp_col = "ID",beta_col = "Beta",se_col = "Se",pval_col = "Pval",effect_allele_col = "ALT",
                       other_allele_col = "REF",eaf_col = "MAF")
exposure_dat1_mvmr$exposure = "GIPR"
exposure_dat2_mvmr = format_data(exposure_dat2_mvmr,snp_col = "SNP",beta_col = "BETA",se_col = "SE",pval_col = "P",effect_allele_col = "EA",
                       other_allele_col = "OA",eaf_col = "EAF")
exposure_dat2_mvmr$exposure = "T2D"
exposure_dat2_mvmr = exposure_dat2_mvmr[,colnames(exposure_dat1_mvmr)]

exposure_mvmr = rbind(exposure_dat1_mvmr,exposure_dat2_mvmr)
exposure_mvmr_end = NULL
for(i in 1:length(IV_candidate)){
   exposure_judge = subset(exposure_mvmr,SNP==IV_candidate[i])
   if(exposure_judge$pval.exposure[1]<=5e-8|exposure_judge$pval.exposure[2]<=5e-8){
    exposure_mvmr_end = rbind(exposure_mvmr_end,exposure_judge)   
   }
}
exposure_mvmr_clump = clump_data(exposure_mvmr_end,clump_kb = 10000,clump_r2 = 0.01)
exposure_mvmr_dat1 = subset(exposure_mvmr_clump,exposure=="GIPR")
exposure_mvmr_dat2 = subset(exposure_mvmr_clump,exposure=="T2D")
mvmr_snp = union(exposure_mvmr_dat1$SNP,exposure_mvmr_dat2$SNP)
exposure_mvmr_dat = subset(xposure_mvmr_end,SNP%in%mvmr_snp)
outcome_dat = read_outcome_data(snps = exposure_mvmr_dat$SNP,filename = "~/GLP1/Analysis/Data/Outcome/NAFLD/NAFLD_Discovery_Processed.txt",sep = "\t",
                               snp_col = "rsName",beta_col = "Beta",se_col = "Se",effect_allele_col = "Amin",other_allele_col = "Amaj",pval_col = "Pval",
							   eaf_col="MAF_PC")
mvdat = mv_harmonise_data(exposure_mvmr_clump,outcome_dat) 
res = mv_multiple(mvdat) 
res


######### GLP1R->BMI->NAFLD #########
library(TwoSampleMR)
exposure_dat1_origin = fread("~/GLP1/BMI_sup/Data/GLP1R.summats.txt")
exposure_dat2_origin = fread("~/GLP1/Analysis/Data/Outcome/Postitive_Control/BMI_GWAS.txt")
colnames(exposure_dat2_origin) = c("CHR","POS","SNP","EA","OA","EAF","BETA","SE","P","N")
exposure_dat1_origin$exposure = "GLP1R"
exposure_dat2_origin$exposure = "BMI"
exposure_dat2_origin$P = as.numeric(exposure_dat2_origin$P)
IV_candidate = intersect(exposure_dat1_origin$ID,exposure_dat2_origin$SNP)

exposure_dat1_mvmr = subset(exposure_dat1_origin,ID%in%IV_candidate)
exposure_dat2_mvmr = subset(exposure_dat2_origin,SNP%in%IV_candidate)


exposure_dat1_mvmr = format_data(exposure_dat1_mvmr,snp_col = "ID",beta_col = "Beta",se_col = "Se",pval_col = "Pval",effect_allele_col = "ALT",
                       other_allele_col = "REF",eaf_col = "MAF")
exposure_dat1_mvmr$exposure = "GLP1R"
exposure_dat2_mvmr = format_data(exposure_dat2_mvmr,snp_col = "SNP",beta_col = "BETA",se_col = "SE",pval_col = "P",effect_allele_col = "EA",
                       other_allele_col = "OA",eaf_col = "EAF")
exposure_dat2_mvmr$exposure = "BMI"
exposure_dat2_mvmr = exposure_dat2_mvmr[,colnames(exposure_dat1_mvmr)]

exposure_mvmr = rbind(exposure_dat1_mvmr,exposure_dat2_mvmr)
exposure_mvmr_end = NULL
for(i in 1:length(IV_candidate)){
   exposure_judge = subset(exposure_mvmr,SNP==IV_candidate[i])
   if(exposure_judge$pval.exposure[1]<=5e-8|exposure_judge$pval.exposure[2]<=5e-8){
    exposure_mvmr_end = rbind(exposure_mvmr_end,exposure_judge)   
   }
}
exposure_mvmr_clump = clump_data(exposure_mvmr_end,clump_kb = 10000,clump_r2 = 0.01)
exposure_mvmr_dat1 = subset(exposure_mvmr_clump,exposure=="GLP1R")
exposure_mvmr_dat2 = subset(exposure_mvmr_clump,exposure=="BMI")
mvmr_snp = union(exposure_mvmr_dat1$SNP,exposure_mvmr_dat2$SNP)
exposure_mvmr_dat = subset(exposure_mvmr_end,SNP%in%mvmr_snp)
outcome_dat = read_outcome_data(snps = exposure_mvmr_dat$SNP,filename = "~/GLP1/Analysis/Data/Outcome/NAFLD/NAFLD_Discovery_Processed.txt",sep = "\t",
                               snp_col = "rsName",beta_col = "Beta",se_col = "Se",effect_allele_col = "Amin",other_allele_col = "Amaj",pval_col = "Pval",
							   eaf_col="MAF_PC")
mvdat = mv_harmonise_data(exposure_mvmr_clump,outcome_dat) 
res = mv_multiple(mvdat) 
res

######### GIPR->BMI->NAFLD #########
library(TwoSampleMR)
library(data.table)
exposure_dat1_origin = fread("~/GLP1/BMI_sup/Data/GIPR.summats.txt")
exposure_dat2_origin = fread("~/GLP1/Analysis/Data/Outcome/Postitive_Control/BMI_GWAS.txt")
colnames(exposure_dat2_origin) = c("CHR","POS","SNP","EA","OA","EAF","BETA","SE","P","N")
exposure_dat1_origin$exposure = "GIPR"
exposure_dat2_origin$exposure = "BMI"
exposure_dat2_origin$P = as.numeric(exposure_dat2_origin$P)
IV_candidate = intersect(exposure_dat1_origin$ID,exposure_dat2_origin$SNP)

exposure_dat1_mvmr = subset(exposure_dat1_origin,ID%in%IV_candidate)
exposure_dat2_mvmr = subset(exposure_dat2_origin,SNP%in%IV_candidate)


exposure_dat1_mvmr = format_data(exposure_dat1_mvmr,snp_col = "ID",beta_col = "Beta",se_col = "Se",pval_col = "Pval",effect_allele_col = "ALT",
                       other_allele_col = "REF",eaf_col = "MAF")
exposure_dat1_mvmr$exposure = "GIPR"
exposure_dat2_mvmr = format_data(exposure_dat2_mvmr,snp_col = "SNP",beta_col = "BETA",se_col = "SE",pval_col = "P",effect_allele_col = "EA",
                       other_allele_col = "OA",eaf_col = "EAF")
exposure_dat2_mvmr$exposure = "BMI"
exposure_dat2_mvmr = exposure_dat2_mvmr[,colnames(exposure_dat1_mvmr)]

exposure_mvmr = rbind(exposure_dat1_mvmr,exposure_dat2_mvmr)
exposure_mvmr_end = NULL
for(i in 1:length(IV_candidate)){
   exposure_judge = subset(exposure_mvmr,SNP==IV_candidate[i])
   if(exposure_judge$pval.exposure[1]<=5e-8|exposure_judge$pval.exposure[2]<=5e-8){
    exposure_mvmr_end = rbind(exposure_mvmr_end,exposure_judge)   
   }
}
exposure_mvmr_clump = clump_data(exposure_mvmr_end,clump_kb = 10000,clump_r2 = 0.01)
exposure_mvmr_dat1 = subset(exposure_mvmr_clump,exposure=="GIPR")
exposure_mvmr_dat2 = subset(exposure_mvmr_clump,exposure=="BMI")
mvmr_snp = union(exposure_mvmr_dat1$SNP,exposure_mvmr_dat2$SNP)
exposure_mvmr_dat = subset(exposure_mvmr_end,SNP%in%mvmr_snp)
outcome_dat = read_outcome_data(snps = exposure_mvmr_dat$SNP,filename = "~/GLP1/Analysis/Data/Outcome/NAFLD/NAFLD_Discovery_Processed.txt",sep = "\t",
                               snp_col = "rsName",beta_col = "Beta",se_col = "Se",effect_allele_col = "Amin",other_allele_col = "Amaj",pval_col = "Pval",
							   eaf_col="MAF_PC")
mvdat = mv_harmonise_data(exposure_mvmr_clump,outcome_dat) 
res = mv_multiple(mvdat) 
res








