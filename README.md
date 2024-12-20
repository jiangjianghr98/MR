#1.1芬兰数据处理：
rm(list=ls())

library(vroom)
library(tidyr)
library(dplyr)
library(data.table)

biofsci <- vroom('finngen_R10_PAIN.gz', col_names = TRUE)
head(biofsci)
colnames(biofsci)
biofsci <- biofsci %>%
  rename(
    SNP = rsids,
    CHR = "#chrom",
    BP = pos,
    effect_allele = alt,
    other_allele = ref,
    P = pval,
    EAF = af_alt,
    BETA = beta,
    SE = sebeta,
  ) %>%
  select(SNP, CHR, BP, effect_allele, other_allele, P, EAF, BETA, SE) %>%
  mutate(P = as.numeric(P)) 

write.csv(biofsci, "pain.csv", row.names = FALSE)

data=fread("pain.csv")



#1.2GWAS数据处理
rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
library("pacman")
p_load(VariantAnnotation,gwasglue,TwoSampleMR,dplyr,tidyr,CMplot)
p_load(data.table,tidyr,dplyr) 


## 方法1: readVcf读取数据
infile="ieu-b-4874.vcf.gz" 
pvalfilter= 5e-08
expo_data_MR1 <- readVcf(infile) %>% 
  gwasvcf_to_TwoSampleMR(type = "exposure") # 使用readVcf()读取VCF文件,并使用gwasvcf_to_TwoSampleMR()转换格式

outcome_data_MR1 <- readVcf(infile) %>% 
  gwasvcf_to_TwoSampleMR(type = "outcome") # 使用readVcf()读取VCF文件,并使用gwasvcf_to_TwoSampleMR()转换格式
              
expo_data_MR1<-subset(expo_data_MR1, pval.exposure<pvalfilter) # 以5e-08为阈值,过滤expo_data_MR数据集
write.csv(expo_data_MR1, file="expo_data_MR1.csv", row.names=F) # 将过滤后的结果out_data_MR写入CSV文件
#write.csv(expo_data_MR1, file="expo_data_MR1.txt", row.names=F) # 将过滤后的结果out_data_MR写入CSV文件

write.csv(outcome_data_MR1, file="outcome_data_MR1.csv", row.names=F) # 将过滤后的结果out_data_MR写入CSV文件


#1.3结局数据处理
rm(list=ls())
library("pacman")
p_load(data.table,tidyr,dplyr) 
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("VariantAnnotation",force = TRUE)

samplesize = 384+294770
mydata = fread("finngen_R10_M13_FORESTIER.gz") 
#mydata = fread("bbj-a-159.vcf.gz",skip = "#CHROM")
head(mydata)
mydata$samplesize= samplesize
fwrite(mydata,"finngen_R10_M13_FORESTIER.csv")

#2.1孟德尔随机化
rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
if (!require("devtools")) {
  install.packages("devtools")
} else {}
if (!require("data.table")) {
  install.packages("data.table")
} else {}
if (!require("TwoSampleMR")) {
  devtools::install_github("MRCIEU/TwoSampleMR") 
} else {}
#remotes::install_github("MRCIEU/MRInstruments")

p_load(data.table,dplyr,tidyr,ggplot2)


#remotes::install_github("yulab-smu/yulab.utils")
#library(yulab.utils)
#install_zip_gh("MRCIEU/ieugwasr")

library("ieugwasr")
library(MRInstruments)
library(plyr)
library(dplyr)
library(data.table)


dataname1="./out/finngen_R10_PAIN.csv"

GWAS_1 <- fread(dataname1)
GWAS_1$PHENO<-"Pain" 
head(GWAS_1)
GWAS_1 <- as.data.frame(GWAS_1)

out_data <- format_data(
  GWAS_1,
  type="outcome",
  phenotype_col = "PHENO",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  eaf_col = "af_alt",
  pval_col = "pval",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  samplesize_col ="samplesize" ,
  chr_col = "#chrom",
  pos_col = "pos")

  rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
if (!require("devtools")) {
  install.packages("devtools")
} else {}

if (!require("data.table")) {
  install.packages("data.table")
} else {}

if (!require("TwoSampleMR")) {
  devtools::install_github("MRCIEU/TwoSampleMR") 
} else {}

#remotes::install_github("MRCIEU/MRInstruments")

p_load(data.table,dplyr,tidyr,ggplot2)


#remotes::install_github("yulab-smu/yulab.utils")
#library(yulab.utils)
#install_zip_gh("MRCIEU/ieugwasr")
library(ieugwasr)
library(MRInstruments)
library(plyr)
library(dplyr)
library(data.table)


pfilter=0.05

if(!dir.exists("ORdata05")){
  dir.create("ORdata05")
}

if(!dir.exists("Pleiotropydata05")){
  dir.create("Pleiotropydata05")
}

if(!dir.exists("Result05")){
  dir.create("Result05")
}

FileNames <- list.files(path="exprosure/", pattern = "*.csv")
FileNames
#FileNames=FileNames[1:2]

exp_dat <- list() 
ex_pore <- c() 
for(i in c(1:length(FileNames))){
  IV <- fread(paste0("exprosure","/",FileNames[i]))
  #IV <- fread(paste0("exprosure","/",FileNames[1]))
  head(IV)
  IV$PHENO <- FileNames[i] 
  #IV$PHENO <- FileNames[1] 
  IV=as.data.frame(IV)
  #head(IV)

  IV1<-format_data(IV,
                   type="exposure",
                   phenotype_col = "PHENO", 
                   snp_col = "rsid",
                   beta_col = "beta",
                   se_col = "standard_error",
                   eaf_col = "standard_error",
                   effect_allele_col = "effect_allele",
                   other_allele_col = "other_allele",
                   pval_col = "p_value",
                   samplesize_col = "n",
                   chr_col = "chromosome",
                   pos_col = "base_pair_location")

  exp_dat[[i]] <- IV1 
  
  ex_pore<-c(ex_pore,FileNames[i])
}
save.image("exposure.Rdata")
load("exposure.Rdata")

dataname1="./out/finngen_R10_PAIN.csv"
GWAS_1 <- fread(dataname1)
head(GWAS_1)
allSNP <- do.call(rbind, exp_dat)
GWAS_2 <- subset(GWAS_1,GWAS_1$rsids %in% allSNP$SNP) 
rm(GWAS_1)
GWAS_2$PHENO<-"Pain" #要修改
head(GWAS_2)
#head(GWAS_1)
GWAS_2 <- as.data.frame(GWAS_2)

out_data <- format_data(
  GWAS_2,
  type="outcome",
  phenotype_col = "PHENO",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  eaf_col = "af_alt",
  pval_col = "pval",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  samplesize_col ="samplesize" ,
  chr_col = "#chrom",
  pos_col = "pos")

out_dat <- list() 
out_dat[[1]] <- out_data
out_come<-c("Pain") 

save.image("outcome.Rdata")

load("outcome.Rdata")

results <- list()
#ex_pore=ex_pore[1]

for (i in c(1:length(ex_pore))){
  #for (i in ex_pore){
  for (j in c(1:length(out_come))){
    # 调和数据
    dat <- harmonise_data(
      exposure_dat = exp_dat[[i]],
      outcome_dat = out_dat[[j]],
      action = 2
    )
    #去除弱工具变量
    beta_squared <-(dat$beta.exposure)^2 
    se_squared <-(dat$se.exposure)^2 
    sample_adjusted <- dat$samplesize.exposure - 2 
    dat$R2 <- (2 * beta_squared) / 
      ((2 * beta_squared) + (2 * dat$samplesize.exposure * se_squared)) 
    dat$f <- dat$R2 * sample_adjusted / (1 - dat$R2) 
    dat$meanf<- mean(dat$f) 
    dat<-dat[dat$f>10,] 

    res <- mr(dat) 
    result_or <- generate_odds_ratios(res)
    
    if(result_or$pval[3]<pfilter){
    #filename <- basename(sub("\\.txt$","",i))
    filename <- basename(sub("\\.txt$","",ex_pore[i])) 
    filename2 <- sub("\\.csv$", "", paste0(filename))
    dir.create(paste0("./Result05/",filename2))
    write.table(dat, 
                file = paste0("./Result05/",filename2,"/harmonise.csv"),
                row.names = F, sep = ",", quote = F)
    write.table(result_or[,5:ncol(result_or)], 
                file = paste0("./Result05/",filename2,"/OR.csv"),
                row.names = F, sep = ",", quote = F)
    
    write.table(result_or[, 5:ncol(result_or)], 
                file = paste0("./ORdata05/",filename, "_OR.csv"), 
                row.names = FALSE, sep = ",", quote = F)          
    p1 <- mr_scatter_plot(res, dat)
    ggsave(p1[[1]], file=paste0("./Result05/",filename2,"/scatter.pdf"), 
           width=8, height=8)
    
    pleiotropy <- mr_pleiotropy_test(dat)
    write.table(pleiotropy, file = paste0("./Result05/",filename2,"/pleiotropy.csv"),
                sep = ",", quote = F, row.names=F)
    
    write.table(pleiotropy, 
                file = paste0("./Pleiotropydata05/",filename, "_pleiotropy.csv"), 
                sep = ",", quote = F, row.names=F)
    
    heterogeneity <- mr_heterogeneity(dat)
    write.table(heterogeneity, file = paste0("./Result05/",filename2,"/heterogeneity.csv"),
                sep = ",", quote = F, row.names=F)
    
    presso <- run_mr_presso(dat, NbDistribution = 1000)
    capture.output(presso, file = paste0("./Result05/",filename2,"/presso.csv"))
    
    singlesnp_res <- mr_singlesnp(dat)
    singlesnpOR <- generate_odds_ratios(singlesnp_res)
    write.table(singlesnpOR, file=paste0("./Result05/",filename2,"/singlesnpOR.csv"),
                row.names = F, sep = ",", quote = F)
    
    p2 <- mr_forest_plot(singlesnp_res)
    ggsave(p2[[1]], file=paste0("./Result05/",filename2,"/forest.pdf"), width=8, height=8)
    
    sen_res <- mr_leaveoneout(dat)
    p3 <- mr_leaveoneout_plot(sen_res)
    ggsave(p3[[1]], file=paste0("./Result05/",filename2,"/sensitivity-analysis.pdf"), 
           width=8, height=8)
    
    res_single <- mr_singlesnp(dat)
    p4 <- mr_funnel_plot(singlesnp_res)
    ggsave(p4[[1]], file=paste0("./Result05/",filename2,"/funnelplot.pdf"), width=8, height=8)

    res$exposure=ex_pore[i]
    res$outcome=out_come[j]
    
    print(paste0("------", ex_pore[i], " & ",out_come[j],"------"))
    print(generate_odds_ratios(res))
    
    results[[length(out_come)*(i-1)+j]] <- generate_odds_ratios(res)
    }
  }
}


results_allIV <- do.call(rbind, results) 
fwrite(results_allIV,"result06.csv")

#2.2绘图

rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
library("pacman")
p_load(forestploter,grid,ggplot2,data.table)
mrresultsfile <- "allOR_Inform06.csv"

#biofsci = read.table(mrresultsfile, header = T, sep = "\t")
biofsci = fread(mrresultsfile, header = T)
colnames(biofsci)

biofsci$pval = as.numeric(biofsci$pval)
biofsci$pval <- ifelse(biofsci$pval<0.05, "<0.05", sprintf("%.4f", biofsci$pval))

biofsci$estimate <- paste0(format(round(biofsci$or, 4), nsmall = 4), " (",  
                           format(round(biofsci$or_lci95, 4), nsmall = 4), "-",
                           format(round(biofsci$or_uci95, 4), nsmall = 4), ")")

biofsci$Trails <-  biofsci$reportedTrait

biofsci$Trails = ifelse(is.na(biofsci$Trails), "", biofsci$Trails)
biofsci$method = ifelse(is.na(biofsci$method), "", biofsci$method)
biofsci$nsnp = ifelse(is.na(biofsci$nsnp), "", biofsci$nsnp)  
biofsci$pval = ifelse(is.na(biofsci$pval), "", biofsci$pval)

colnames(biofsci)
biofsci$` ` <- paste(rep(" ", 15), collapse = " ")   
biofsci$`OR(95%CI)` <- biofsci$estimate

colnames(biofsci)
biofsci$Trails[duplicated(biofsci$Trails)] <- " "
biofsci1=biofsci[,c("Trails","method", "nsnp","pval","OR(95%CI)"," ")]

colnames(biofsci1)
                                     
tm <- forest_theme(base_size = 10,
                   refline_col = "red",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")
  
colnames(biofsci1)
pdf("allforestall07.pdf", width=10, height=9)
p=forest(biofsci1,  
         est = biofsci$or,
         lower = biofsci$or_lci95,
         upper = biofsci$or_uci95,
         #arrow_lab = c("Placebo Better", "Treatment Better"),
         sizes = 0.4,
         ci_column = 6,
         ref_line = 1,
         xlim = c(0.8, 1.2),
         ticks_at = c(0.8, 0.9, 1, 1.1, 1.2),
         footnote = "",
         theme = tm)
p
dev.off()

#2.3反向孟德尔

  library(ieugwasr)
  library(MRInstruments)
  library(plyr)
  library(dplyr)
  library(data.table)
  library(TwoSampleMR)
  library(ggplot2)
  library(foreach)
  library(purrr)
  

if(!dir.exists("Result10")){ 
  dir.create("Result10") 
}


#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)
ivwFile="ivw68-815.csv"    

biof=read.csv(ivwFile, header=T, sep=",", check.names=F)
outcomeid=biof$id
outcomeid
#outcomeid=outcomeid[1]

exposureID=unique(biof$id.outcome)
exposure_dat=fread("GCST90200815_buildGRCh38_clump.txt")
exposure_dat$samplesize= 4488
exposure_dat$PHENO="DESEASE"
fwrite(exposure_dat,"finn.csv")
head(exposure_dat)

exposure_dat <- read_exposure_data("finn.csv", 
                              sep = ",",       
                              #phenotype_col = "PHENO", 
                              snp_col = "variant_id",     
                              beta_col = "beta", 
                              se_col = "standard_error",   
                              effect_allele_col = "effect_allele", 
                              other_allele_col = "other_allele", 
                              pval_col = "p_value",
                              eaf_col = "effect_allele_frequency", 
                              samplesize_col = "samplesize", 
                              chr_col = "#chrom",
                              pos_col = "pos",
                              clump=F 
)  


workdir <- getwd() 
workdir <- paste0(getwd(),"/out")
workdir

outcomefile <- list.files(path = workdir, pattern = "^G.*\\.tsv", full.names = TRUE)

print(outcomefile)

result=c()
for(i in outcomefile){
  tryCatch({  

    #outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=i)
    outdata=fread(i)
    #outdata=fread("out/GCST90274770.tsv.gz")
    outdata$PHENO <- gsub("\\.tsv.gz$", "", basename(i))
    #outdata$PHENO <- outcomefile[i] 
    outdata=as.data.frame(outdata)
    #outdata$pval <- 10^(-outdata$neg_log_10_p_value)
    colnames(outdata)
    outc_data <- format_data(
      outdata,
      type="outcome",
      snps=exposure_dat$SNP,
      header=T,
      #phenotype_col="trait",
      phenotype_col="PHENO",
      snp_col = "rsid",
      beta_col = "beta",
      se_col = "standard_error", 
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "p_value",
      eaf_col = "effect_allele_frequency",
      chr_col = "chromosome",
      pos_col = "base_pair_location"
    )
 
      #expo_snps <- exposure_dat$SNP
      #outc_snps <- outc_data$SNP

    #outc_data$outcome="Colorectal cancer"
    dat <- harmonise_data(exposure_dat, outc_data)
    
    

    #presso=run_mr_presso(dat)
    #write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file=paste0(i, ".table.MR-PRESSO_Global.csv"))
    #write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file=paste0(i, ".table.MR-PRESSO_Outlier.csv"))
    

    mr_result=mr(dat)
    result_or=generate_odds_ratios(mr_result)
    
    #newid <- sub("_.tsv.gz", "", basename(i))
    newid <- gsub("\\.tsv.gz$", "", basename(i))

    filename <- paste0("Result10/", newid) 
    if (!dir.exists(filename)) { 
      dir.create(filename) 
    }
    
     
    filename2 <- paste0("Result11/OR/", newid) 
   
    fwrite(dat, file = paste0(filename,"/harmonise.csv"),sep = ",", quote = F)
    fwrite(result_or, file = paste0(filename,"/OR.csv"),sep = ",", quote = F,row.names=F)
    pleiotropy <- mr_pleiotropy_test(dat)
    write.table(pleiotropy, file = paste0(filename, "/pleiotropy.csv"),sep = ",", quote = F,row.names=F)
    
    heterogeneity <- mr_heterogeneity(dat)
    write.table(heterogeneity, file = paste0(filename, "/heterogeneity.csv"),sep = ",", quote = F,row.names=F)
    
    #presso <- run_mr_presso(dat, NbDistribution = 1000)
    #capture.output(presso, file = paste0(filename, "/presso.txt"),sep = ",", quote = F,row.names=F)
    
    singlesnp_res <- mr_singlesnp(dat)
    singlesnpOR <- generate_odds_ratios(singlesnp_res)
    write.table(singlesnpOR,file=paste0(filename,"/singlesnpOR.csv"),sep = "\t",quote = F, row.names=F)
    
    p1 <- mr_scatter_plot(mr_result, dat)
    ggsave(p1[[1]], file=paste0(filename,"/scatter.pdf"), width=8, height=8)
    singlesnp_res<- mr_singlesnp(dat)
    singlesnpOR=generate_odds_ratios(singlesnp_res)
    write.table(singlesnpOR,file=paste0(filename,"/singlesnpOR.csv"),sep = ",", quote = F, row.names=F)
    p2 <- mr_forest_plot(singlesnp_res)
    ggsave(p2[[1]], file=paste0(filename,"/forest.pdf"), width=8, height=8)
    sen_res<- mr_leaveoneout(dat)
    p3 <- mr_leaveoneout_plot(sen_res)
    ggsave(p3[[1]], file=paste0(filename,"/sensitivity-analysis.pdf"), width=8, height=8)
    res_single <- mr_singlesnp(dat)
    p4 <- mr_funnel_plot(singlesnp_res)
    ggsave(p4[[1]], file=paste0(filename,"/funnelplot.pdf"), width=8, height=8)
    list(result_or = result_or, singlesnp_res = singlesnp_res, sen_res = sen_res)
    
    result=c(result, mr_result$pval[3])

   # if (!is.na(result_or$pval[3]) && result_or$pval[3] < 0.05) {
      if (T) {

      outTab=dat[dat$mr_keep=="TRUE",]
    }
  }, error = function(e) {

    cat("Error occurred for file:", i, "\n")
    cat("Error message:", conditionMessage(e), "\n")
  })
}

#3.1计算中介效应值

p_load(data.table,dplyr,tidyr,ggplot2)


#读取OR17_beta1文件夹下所有OR_beta1.csv
beta1FileList=list.files(path = "OR17_beta1", pattern = "*.csv", full.names = T)
beta1FileList

#读取ORdata19_beta2文件夹下所有OR_beta2.csv
beta2FileList=list.files(path = "ORdata19_beta2", pattern = "*.csv", full.names = T)
beta2FileList

#读取ORdata20_betaAll文件夹下所有OR_betaAll.csv
betaAllFileList=list.files(path = "ORdata20_betaAll", pattern = "*.csv", full.names = T)
betaAllFileList

for (i in 1:length(beta1FileList)) {
  beta1File=beta1FileList[i]
  
  file_parts <- strsplit(basename(beta1File), split = "_")
  
  second_part <- file_parts[[1]][2]
  beta2File=paste0("ORdata19_beta2/",second_part,"_OR_beta2.csv")
  
  first_part <- file_parts[[1]][1]
  betaAllFile=paste0("ORdata20_betaAll/",first_part,"_OR_betaAll.csv")
  
  output_folder <- "medResult21"
  if (!file.exists(output_folder)) {
    dir.create(output_folder)
  }
  
  if (file.exists(beta2File) & file.exists(betaAllFile)) {
    
    #beta1File="ebi-a-GCST90001662_GCST90200850_OR_beta1.csv"   
    #beta2File="GCST90200850_OR_beta2.csv"                   
    #betaAllFile="ebi-a-GCST90001662_OR_betaAll.csv"           
    
    
    extract_data <- function(file_name) {
      data <- read.csv(file_name, header=T, sep=",", check.names=F)
      b <- data[3, "b"]
      se <- data[3, "se"]
      exposures <- unique(data[,"exposure"])
      return(list(b = b, se = se, exposures = exposures))
    }
    
    beta1 = extract_data(beta1File)
    beta2 = extract_data(beta2File)
    betaAllRT = read.csv(betaAllFile)
    bataAll=betaAllRT[3,"b"]
    
    outcomes = unique(betaAllRT[,"outcome"])
    calculate_effects = function(b1, se1, b2, se2, bAll) {
      beta12 = b1 * b2
      se12 = sqrt((b1 ^ 2) * (se1 ^ 2) + (b2 ^ 2) * (se2 ^ 2))
      ciLow= beta12 -1.96 * se12
      ciHigh = beta12 +1.96 * se12
      betaDirect = bAll - beta12
      beta12Ratio = beta12 / bAll
      ciLowRatio = ciLow / bAll
      ciHighRatio = ciHigh / bAll
      zValue = beta12 / se12
      pvalue = 2 * pnorm(q=abs(zValue), lower.tail = FALSE)
      
      return(list(beta12 = beta12, se12 = se12, ciLow = ciLow, ciHigh = ciHigh, betaDirect = betaDirect, beta12Ratio = beta12Ratio, ciLowRatio = ciLowRatio, ciHighRatio = ciHighRatio, pvalue = pvalue))
    }
    
    effects = calculate_effects(beta1$b, beta1$se, beta2$b, beta2$se, bataAll)
    result <- data.frame(
      "Gut Microbes" = beta1$exposures,
      "Immune cell" = beta2$exposures,
      "Outcome" = outcomes,
      "Mediated effect" = paste0(signif(effects$beta12, digits = 3), "(", signif(effects$ciLow, digits = 3), ", ", signif(effects$ciHigh, digits = 3), ")"),
      "Mediated proportion" = paste0(signif(effects$beta12Ratio * 100, digits = 3), "%", "(", signif(effects$ciLowRatio * 100, digits = 3), "%, ", signif(effects$ciHighRatio * 100, digits = 3), "%)")
     # "p value" = effects$pvalue
    )
    
    
    
    my_string <- basename(beta1File)
    new_string <- sub("_OR_beta1.csv", "", my_string)
    
    write.csv(result, file = paste0(output_folder, "/", sub("_OR_beta1.csv", "", basename(beta1File)), "_Med_mr_Effect21", ".csv"), row.names=F)
  } else {
    warning(paste("Files", beta2File, "and/or", betaAllFile, "do not exist. Skipping this file."))
    next
  }
}

#3.2绘图

library(grid)
library(readr)
library(forestploter)

biofsci <- do.call(rbind, lapply(grep(".csv$", dir(), value = TRUE),
                              function(i) read.csv(i, header = T, sep = ",", check.names = F)))

lineVec <- cumsum(c(1, table(biofsci[, c('exposureTrait', 'outcomeTrait')])))
biofsci$' ' <- paste(rep(" ", 10), collapse = " ")
biofsci$'OR(95% CI)' <- ifelse(is.na(biofsci$or), "", sprintf("%.3f (%.3f to %.3f)", biofsci$or, biofsci$or_lci95, biofsci$or_uci95))
biofsci$pval <- ifelse(biofsci$pval < 0.001, "<0.001", sprintf("%.3f", biofsci$pval))
biofsci$exposure <- ifelse(is.na(biofsci$exposure), "", biofsci$exposure)
biofsci$nsnp <- ifelse(is.na(biofsci$nsnp), "", biofsci$nsnp)
biofsci[duplicated(biofsci[, c('exposureTrait','outcomeTrait')]), c('exposureTrait', 'outcomeTrait')] <- ""

tm <- forestploter::forest_theme(base_size = 15, ci_pch = 15, ci_lty = 1, ci_lwd = 1.5,
                                 ci_col = "black", ci_Theight = 0.2,
                                 refline_lty = "dashed", refline_lwd = 1, 
                                 refline_col = "grey40", xaxis_cex = 0.8, 
                                 footnote_cex = 0.6, footnote_col = "darkblue") 

plot <- forestploter::forest(biofsci[, c("exposureTrait","outcomeTrait","nsnp","method","pval", " ","OR(95% CI)")],
                             est = biofsci$or, lower = biofsci$or_lci95, upper = biofsci$or_uci95, 
                             ci_column = 6, ref_line = 1, xlim = c(0.5, 1.5), theme = tm)

nature_colors <- c("#009E73", "#56B4E9", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

boxcolor <- nature_colors[as.numeric(as.factor(biofsci$method)) %% length(nature_colors) + 1]

for(i in 1:nrow(biofsci)) {
  plot <- edit_plot(plot, col = 6, row = i, which = "ci", gp = gpar(fill = boxcolor[i], fontsize=25))
}

for(i in which(as.numeric(gsub('<', "", biofsci$pval)) < 0.05)) {
  plot <- edit_plot(plot, col=5, row = i, which = "text", gp = gpar(fontface="bold"))
}

plot <- add_border(plot, part = "header", row = c(1, lineVec), gp = gpar(lwd = 1))
plot <- edit_plot(plot, col = 1:ncol(biofsci), row = 1:nrow(biofsci), which = "text", gp = gpar(fontsize=12), x = unit(0.5, "npc"), hjust = unit(0.5, "npc"))
plot <- edit_plot(plot, col = 1:ncol(biofsci), which = "text", hjust = unit(0.5, "npc"), x = unit(0.5, "npc"), part="header")

pdf("mr_forest.pdf", width = 28, height = 8)
print(plot)
dev.off()

