
################################ Part 0 Data Preprocessing ################################

#This part aims to perform data preprocessing in order to save computing time or to regularize the data format.
#For Gtex data, we choose ANNOVAR to match chr:pos and rsID for SNPs, and here we just take pancreas tissue as an example for simplicity. 


######################## Step 1 Outcome Data ########################

library(TwoSampleMR)
library(data.table)

################ NAFLD ################

setwd("~/GLP1/Analysis/Data/Outcome/NAFLD")
library(data.table)
Outcome_Discovery  = fread("NAFL_Discovery.txt")
Outcome_Discovery$Beta = log(Outcome_Discovery$Effect)
Outcome_Discovery$Se = abs(log(Outcome_Discovery$Effect)/qnorm(Outcome_Discovery$Pval/2))
Outcome_Discovery$MAF_PC = (Outcome_Discovery$MAF_PC)/100
write.table(Outcome_Discovery,"~/GLP1/Analysis/Data/Outcome/NAFLD/NAFLD_Discovery_Processed.txt",quote=FALSE,row.names=FALSE,sep="\t")

################ Cirrhosis ################

setwd("~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits")
library(data.table)
Cirrhosis  = fread("Cirrhosis.txt")
Cirrhosis$Beta = log(Cirrhosis$Effect)
Cirrhosis$Se = abs(log(Cirrhosis$Effect)/qnorm(Cirrhosis$Pval/2))
Cirrhosis$MAF_PC = (Cirrhosis$MAF_PC)/100
write.table(Cirrhosis,"~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits/Cirrhosis_Processed.txt",quote=FALSE,row.names=FALSE,sep="\t")

################ PDFF ################

setwd("~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits")
library(data.table)
PDFF = fread("PDFF.txt")
PDFF$Se = sqrt(((PDFF$Effect)^2)/qchisq(PDFF$Pval, 1, lower.tail = F))
PDFF$MAF_PC = (PDFF$MAF_PC)/100
write.table(PDFF,"~/GLP1/Analysis/Data/Outcome/NAFLD_Related_Traits/PDFF_Processed.txt",quote=FALSE,row.names=FALSE,sep="\t")


######################## Step 2 Exposure Data ########################

library(TwoSampleMR)
library(data.table)

################ Downstream Biomarker - HbA1c ################
data_reponse = fread("~/GLP1/Supplementary_MR/Data/HbA1c/D1.tsv.gz")
data_reponse_GLP1R = subset(data_reponse,CHR == 6 & BP >39016557-100000 & BP<39059079+100000)
data_reponse_GLP1R$P = as.numeric(data_reponse_GLP1R$P)
data_reponse_GLP1R = subset(data_reponse_GLP1R,P<5e-8)
data_reponse_GLP1R$Sample_Size = 438069
GLP1R_dat = format_data(data_reponse_GLP1R,snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",
                        other_allele_col = "ALLELE0",pval_col = "P",eaf_col="A1FREQ",chr_col = "chromosome",
                        samplesize_col = "Sample_Size")
GLP1R_dat$PVE = 2*(GLP1R_dat$beta.exposure)^2/(2*(GLP1R_dat$beta.exposure)^2+2*GLP1R_dat$samplesize.exposure*(GLP1R_dat$se.exposure)^2)
GLP1R_dat$F_Statistics = ((GLP1R_dat$samplesize.exposure-2)*GLP1R_dat$PVE)/(1-GLP1R_dat$PVE)
GLP1R_dat = subset(GLP1R_dat,F_Statistics>10)
GLP1R_clump<-clump_data(GLP1R_dat,clump_kb = 10000,clump_r2 = 0.3)
write.table(data_reponse_GLP1R,"~/GLP1/Analysis/Data/Exposure/Downstream/HbA1c_GLP1R_clump.txt",quote=FALSE,row.names=FALSE,sep="\t")

data_reponse_GIPR = subset(data_reponse,CHR == 19 & BP >46171479-100000 & BP<46186980+100000 )
data_reponse_GIPR$P = as.numeric(data_reponse_GIPR$P)
data_reponse_GIPR = subset(data_reponse_GIPR,P<5e-8)
data_reponse_GIPR$Sample_Size = 438069
GIPR_dat = format_data(data_reponse_GIPR,snp_col = "SNP",beta_col = "BETA",se_col="SE",effect_allele_col = "ALLELE1",
                       other_allele_col = "ALLELE0",pval_col = "P",eaf_col="A1FREQ",chr_col = "chromosome",
                       samplesize_col = "Sample_Size")
GIPR_dat$PVE = 2*(GIPR_dat$beta.exposure)^2/(2*(GIPR_dat$beta.exposure)^2+2*GIPR_dat$samplesize.exposure*(GIPR_dat$se.exposure)^2)
GIPR_dat$F_Statistics = ((GIPR_dat$samplesize.exposure-2)*GIPR_dat$PVE)/(1-GIPR_dat$PVE)
GIPR_dat = subset(GIPR_dat,F_Statistics>10)
GIPR_clump<-clump_data(GIPR_dat,clump_kb = 10000,clump_r2 = 0.3)
write.table(GIPR_clump,"~/GLP1/Analysis/Data/Exposure/Downstream/HbA1c_GIPR_clump.txt",quote=FALSE,row.names = FALSE,sep="\t")

combine = rbind(GLP1R_clump,GIPR_clump)
combine$id.exposure = "Combine"
write.table(combine,"~/GLP1/Analysis/Data/Exposure/Downstream/HbA1c_Combine_clump.txt",quote=FALSE,row.names = FALSE,sep="\t")

################ eQTLGene ################

exposure_extract = function(Gene_ENS,Gene_start,Gene_stop,Gene_chr,Gene_ref){
  #This function aims to extract valid instruments for Gene expression by standard we defined in the paper and conventional significant threshold 
  #parameter@Gene_ENS: Ensembl ID for the target gene
  #parameter@Gene_start: Transcription Start Site of target Gene(hg19/GRCh37)
  #parameter@Gene_stop: Transcription end site of target Gene(hg19/GRCh37)
  #parameter@Gene_chr: Chromosome of target gene
  #parameter@Gene_ref: eQTL from eQTLGene released in 2019
  #output: A dataframe including valid instrument variables

  eqtl_gwas<-extract_instruments(outcome=paste0("eqtl-a-",Gene_ENS),clump = FALSE)
  eqtl_gwasc<-clump_data(eqtl_gwas,clump_kb = 10000,clump_r2 = 0.3,pop = "EUR")
  gene_gwas<-subset(eqtl_gwasc,chr.exposure== Gene_chr & pos.exposure>Gene_start-100000 & pos.exposure<Gene_stop+100000)
  gene_gwas$PVE = 2*(gene_gwas$beta.exposure)^2/(2*(gene_gwas$beta.exposure)^2+2*gene_gwas$samplesize.exposure*(gene_gwas$se.exposure)^2)
  gene_gwas$F_Statistics = ((gene_gwas$samplesize.exposure-2)*gene_gwas$PVE)/(1-gene_gwas$PVE)
  gwas_exposure = subset(gene_gwas,F_Statistics>=10)
                                 
  gene_origin = subset(Gene_ref,(SNP%in%gwas_exposure$SNP)&(Gene==Gene_ENS))
  gene_gwas = merge(gwas_exposure,gene_origin,by.x=c("SNP","effect_allele.exposure","other_allele.exposure"),by.y=c("SNP","AssessedAllele","OtherAllele"))
  gene_end = gene_gwas[,1:17]
                                 
  return(gene_end)
}

Gene_information = read.table("~/GLP1/Analysis/Data/Exposure/Geneinformation.txt",head=TRUE) #This file includes all the information required of function 'exposure_extract'
Gene_2019 = fread("~/GLP1/Analysis/Data/Exposure/Gene_expression/eQTLGene_2019.txt")
for(i in 1:dim(Gene_information)[1]){
  Name_gene = Gene_information[i,1]
  ENS_gene = Gene_information[i,5]
  Chr_gene = Gene_information[i,4]
  Start_gene = Gene_information[i,2]
  Stop_gene = Gene_information[i,3]
  IV_for_gene = exposure_extract(ENS_gene,Start_gene,Stop_gene,Chr_gene,Gene_2019)
  write.table(IV_for_gene,paste0(Name_gene,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
}

################ Gtex ################

setwd("~/GLP1/data/Exposure/Gtex/tissue")
tissue_gene = fread("Pancreas.v8.egenes.txt.gz")
gene_ens = subset(tissue_gene,gene_name=="GLP1R")$gene_id
gene_start = subset(tissue_gene,gene_name=="GLP1R")$gene_start 
gene_end = subset(tissue_gene,gene_name=="GLP1R")$gene_end
gene_tissue = data.frame(subset(fread("Pancreas.v8.signif_variant_gene_pairs.txt.gz"),gene_id==gene_ens))

gene_snp_information = strsplit(as.character(gene_tissue[,1]),"_")
gene_tissue$CHR = sapply(gene_snp_information,'[',1)
gene_tissue$POS = sapply(gene_snp_information,'[',2)
gene_tissue$NEA = sapply(gene_snp_information,'[',3)
gene_tissue$EA = sapply(gene_snp_information,'[',4)
gene_cis = subset(gene_tissue,POS>=gene_start-100000&POS<=gene_end+100000)

gene_for_annovar = gene_cis[,c("CHR","POS","POS","NEA","EA")]
fwrite(gene_for_annovar,"~/GLP1/Analysis/Data/Exposure/Gtex/match_tmp/GLP1R_cis.avinput",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
system("~/annovar/annotate_variation.pl ~/GLP1/Analysis/Data/Exposure/Gtex/match_tmp/GLP1R_cis.avinput ~/GLP1/BMI_sup/Annotation/ -filter -build hg38 -dbtype avsnp150")

gene_annovar_match = fread("~/GLP1/Analysis/Data/Exposure/Gtex/match_tmp/GLP1R_cis.avinput.hg38_avsnp150_dropped")
gene_for_GLP1R = merge(gene_cis,gene_annovar_match,by.x=c("POS","NEA","EA"),by.y=c("V5","V6","V7"))
GLP1R_dat = format_data(gene_for_GLP1R,snp_col = "V2",beta_col = "slope",se_col="slope_se",effect_allele_col = "EA",other_allele_col = "NEA",
		                pval_col = "pval_nominal",eaf_col="maf",pos_col="POS",chr_col="CHR")
GLP1R_datc = clump_data(GLP1R_dat,clump_kb = 10000,clump_r2 = 0.3)
write.table(GLP1R_datc,"~/GLP1/Analysis/Data/Exposure/Gtex/GLP1R.txt",quote=FALSE,row.names=FALSE,sep="\t")