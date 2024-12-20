#1.1芬兰数据处理：
rm(list=ls())
# 载入所有需要的库
library(vroom)
library(tidyr)
library(dplyr)
library(data.table)

#读取文件
#芬兰数据网址 https://finngen.gitbook.io/documentation/data-description
biofsci <- vroom('finngen_R10_PAIN.gz', col_names = TRUE)
#查看列名
head(biofsci)
colnames(biofsci)
# 重命名列，并转换数据类型
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
  mutate(P = as.numeric(P)) # 转换P值为数值型

# 保存数据到新的csv文件
write.csv(biofsci, "pain.csv", row.names = FALSE)

data=fread("pain.csv")



#1.2GWAS数据处理
rm(list=ls())
# 判断是否已经安装了“pacman”包，如果没有就安装它
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
# 设置Bioconductor镜像地址为中国科技大学的镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
# 加载“pacman”包，用于方便加载其他的R包
library("pacman")
# 使用“p_load”函数加载所需的R包
# 在加载这些包之前，“pacman”会先检查是否已经安装，如果没有会自动安装
p_load(VariantAnnotation,gwasglue,TwoSampleMR,dplyr,tidyr,CMplot)
# 批量安装加载需要的包
p_load(data.table,tidyr,dplyr) 


## 方法1: readVcf读取数据
# 输入文件
infile="ieu-b-4874.vcf.gz" # 这里定义了输入的VCF格式文件
pvalfilter= 5e-08

# 作为暴露文件,并进行格式转换
expo_data_MR1 <- readVcf(infile) %>% 
  gwasvcf_to_TwoSampleMR(type = "exposure") # 使用readVcf()读取VCF文件,并使用gwasvcf_to_TwoSampleMR()转换格式

# 作为结局文件,并进行格式转换
outcome_data_MR1 <- readVcf(infile) %>% 
  gwasvcf_to_TwoSampleMR(type = "outcome") # 使用readVcf()读取VCF文件,并使用gwasvcf_to_TwoSampleMR()转换格式

# 根据p值阈值过滤                
expo_data_MR1<-subset(expo_data_MR1, pval.exposure<pvalfilter) # 以5e-08为阈值,过滤expo_data_MR数据集


# 将过滤后的数据写入CSV文件
write.csv(expo_data_MR1, file="expo_data_MR1.csv", row.names=F) # 将过滤后的结果out_data_MR写入CSV文件
#write.csv(expo_data_MR1, file="expo_data_MR1.txt", row.names=F) # 将过滤后的结果out_data_MR写入CSV文件

write.csv(outcome_data_MR1, file="outcome_data_MR1.csv", row.names=F) # 将过滤后的结果out_data_MR写入CSV文件


#1.3结局数据处理
rm(list=ls())
# 判断是否已经安装了“pacman”包，如果没有就安装它
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
# 设置Bioconductor镜像地址为中国科技大学的镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
# 加载“pacman”包，用于方便加载其他的R包
# 加载pacman包
library("pacman")
# 批量安装加载需要的包
p_load(data.table,tidyr,dplyr) 
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("VariantAnnotation",force = TRUE)


#样本量
samplesize = 384+294770
mydata = fread("finngen_R10_M13_FORESTIER.gz") #读取数据
#mydata = fread("bbj-a-159.vcf.gz",skip = "#CHROM")
head(mydata)
mydata$samplesize= samplesize
fwrite(mydata,"finngen_R10_M13_FORESTIER.csv")

#2.1孟德尔随机化
rm(list=ls())
# 判断是否已经安装了“pacman”包，如果没有就安装它
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
# 设置Bioconductor镜像地址为中国科技大学的镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
# 加载“pacman”包，用于方便加载其他的R包
# 如果没有 devtools 包,则安装 devtools 包
if (!require("devtools")) {
  install.packages("devtools")
} else {}

# 如果没有 data.table 包,则安装 data.table 包  
if (!require("data.table")) {
  install.packages("data.table")
} else {}

# 如果没有 TwoSampleMR 包,则从 GitHub 安装 TwoSampleMR 包
if (!require("TwoSampleMR")) {
  devtools::install_github("MRCIEU/TwoSampleMR") 
} else {}

# 从 GitHub 安装 MRInstruments 包
#remotes::install_github("MRCIEU/MRInstruments")

# 在加载这些包之前，“pacman”会先检查是否已经安装，如果没有会自动安装
p_load(data.table,dplyr,tidyr,ggplot2)

##  安装包
#remotes::install_github("yulab-smu/yulab.utils")
#library(yulab.utils)
#install_zip_gh("MRCIEU/ieugwasr")
# 导入ieugwasr包
library("ieugwasr")
library(MRInstruments)
library(plyr)
library(dplyr)
library(data.table)


dataname1="./out/finngen_R10_PAIN.csv"
# 读取GWAS数据 
GWAS_1 <- fread(dataname1)
GWAS_1$PHENO<-"Pain" #要修改
head(GWAS_1)
GWAS_1 <- as.data.frame(GWAS_1)
# 格式化为outcome数据，要修改
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
# 判断是否已经安装了“pacman”包，如果没有就安装它
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
# 设置Bioconductor镜像地址为中国科技大学的镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
# 加载“pacman”包，用于方便加载其他的R包
# 如果没有 devtools 包,则安装 devtools 包
if (!require("devtools")) {
  install.packages("devtools")
} else {}

# 如果没有 data.table 包,则安装 data.table 包  
if (!require("data.table")) {
  install.packages("data.table")
} else {}

# 如果没有 TwoSampleMR 包,则从 GitHub 安装 TwoSampleMR 包
if (!require("TwoSampleMR")) {
  devtools::install_github("MRCIEU/TwoSampleMR") 
} else {}

# 从 GitHub 安装 MRInstruments 包
#remotes::install_github("MRCIEU/MRInstruments")

# 在加载这些包之前，“pacman”会先检查是否已经安装，如果没有会自动安装
p_load(data.table,dplyr,tidyr,ggplot2)

##  安装包
#remotes::install_github("yulab-smu/yulab.utils")
#library(yulab.utils)
#install_zip_gh("MRCIEU/ieugwasr")
# 导入ieugwasr包
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

# 获取当前工作目录下所有后缀为.txt的文件,存入 FileNames 列表中
FileNames <- list.files(path="exprosure/", pattern = "*.csv")
FileNames
#FileNames=FileNames[1:2]
# ---暴露数据 -------------------------------------
# 创建空列表用于存储exposure数据
exp_dat <- list() 
# 创建空向量用于存储表型
ex_pore <- c()
# 遍历 FileNames 列表  
for(i in c(1:length(FileNames))){
  # 读取数据
  IV <- fread(paste0("exprosure","/",FileNames[i]))
  #IV <- fread(paste0("exprosure","/",FileNames[1]))
  # 添加表型列
  head(IV)
  IV$PHENO <- FileNames[i] 
  #IV$PHENO <- FileNames[1] 
  IV=as.data.frame(IV)
  #head(IV)
  # 格式化为 TwoSampleMR 所需格式
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
  # 将格式化后数据存入列表
  exp_dat[[i]] <- IV1 
  
  # 将表型加入表型向量
  ex_pore<-c(ex_pore,FileNames[i])
}
# 保存
save.image("exposure.Rdata")
# 读取
load("exposure.Rdata")

#####结局数据
# GWAS数据文件
dataname1="./out/finngen_R10_PAIN.csv"
# 读取GWAS数据 
GWAS_1 <- fread(dataname1)
head(GWAS_1)
# 绑定所有exposure数据
allSNP <- do.call(rbind, exp_dat)
# 保留GWAS数据中在exposure数据中的SNP，可能是GWAS_1$rsids或者GWAS_1$SNP
GWAS_2 <- subset(GWAS_1,GWAS_1$rsids %in% allSNP$SNP) 
# 删除GWAS_1释放内存
rm(GWAS_1)
# 添加表型列
GWAS_2$PHENO<-"Pain" #要修改
head(GWAS_2)
#head(GWAS_1)
GWAS_2 <- as.data.frame(GWAS_2)
# 格式化为outcome数据，要修改
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

# 将outcome数据存入列表
out_dat <- list() 
out_dat[[1]] <- out_data
# 定义outcome表型向量
out_come<-c("Pain") 
# 保存
save.image("outcome.Rdata")
# 读取
load("outcome.Rdata")
# 创建空列表存储结果
results <- list()
#ex_pore=ex_pore[1]


# 遍历exposure和outcome组合
for (i in c(1:length(ex_pore))){
  #for (i in ex_pore){
  for (j in c(1:length(out_come))){
    # 调和数据
    dat <- harmonise_data(
      exposure_dat = exp_dat[[i]],
      outcome_dat = out_dat[[j]],
      action = 2
    )
   # dat<-subset(dat,mr_keep==TRUE) #删除不符合要求的数据
    #去除弱工具变量
    beta_squared <-(dat$beta.exposure)^2 #beta的平方
    se_squared <-(dat$se.exposure)^2 #se的平方
    sample_adjusted <- dat$samplesize.exposure - 2 #样本量减2
    dat$R2 <- (2 * beta_squared) / 
      ((2 * beta_squared) + (2 * dat$samplesize.exposure * se_squared)) #计算R2
    dat$f <- dat$R2 * sample_adjusted / (1 - dat$R2) #计算f
    dat$meanf<- mean(dat$f) #计算f的均值
    dat<-dat[dat$f>10,] #去除f小于10的数据

    # MR分析
    res <- mr(dat)
    # 生成优比值  
    result_or <- generate_odds_ratios(res)
    
    # 保存结果
    if(result_or$pval[3]<pfilter){
    #filename <- basename(sub("\\.txt$","",i))
    filename <- basename(sub("\\.txt$","",ex_pore[i])) 
    filename2 <- sub("\\.csv$", "", paste0(filename))
    dir.create(paste0("./Result05/",filename2))
    # 输出结果
    write.table(dat, 
                file = paste0("./Result05/",filename2,"/harmonise.csv"),
                row.names = F, sep = ",", quote = F)
    write.table(result_or[,5:ncol(result_or)], 
                file = paste0("./Result05/",filename2,"/OR.csv"),
                row.names = F, sep = ",", quote = F)
    
    # 保存 odds ratio 数据，并使用 filename 作为前缀
    write.table(result_or[, 5:ncol(result_or)], 
                file = paste0("./ORdata05/",filename, "_OR.csv"), 
                row.names = FALSE, sep = ",", quote = F)
    # 绘制散点图             
    p1 <- mr_scatter_plot(res, dat)
    ggsave(p1[[1]], file=paste0("./Result05/",filename2,"/scatter.pdf"), 
           width=8, height=8)
    
    # 进行多重比较校正
    pleiotropy <- mr_pleiotropy_test(dat)
    write.table(pleiotropy, file = paste0("./Result05/",filename2,"/pleiotropy.csv"),
                sep = ",", quote = F, row.names=F)
    
    #单独输出多效性数据
    write.table(pleiotropy, 
                file = paste0("./Pleiotropydata05/",filename, "_pleiotropy.csv"), 
                sep = ",", quote = F, row.names=F)
    
    # 进行异质性检验
    heterogeneity <- mr_heterogeneity(dat)
    write.table(heterogeneity, file = paste0("./Result05/",filename2,"/heterogeneity.csv"),
                sep = ",", quote = F, row.names=F)
    
    # MR-PRESSO
    presso <- run_mr_presso(dat, NbDistribution = 1000)
    capture.output(presso, file = paste0("./Result05/",filename2,"/presso.csv"))
    
    # 单个SNP分析 
    singlesnp_res <- mr_singlesnp(dat)
    singlesnpOR <- generate_odds_ratios(singlesnp_res)
    write.table(singlesnpOR, file=paste0("./Result05/",filename2,"/singlesnpOR.csv"),
                row.names = F, sep = ",", quote = F)
    
    # 森林图
    p2 <- mr_forest_plot(singlesnp_res)
    ggsave(p2[[1]], file=paste0("./Result05/",filename2,"/forest.pdf"), width=8, height=8)
    
    # Leave-one-out分析
    sen_res <- mr_leaveoneout(dat)
    p3 <- mr_leaveoneout_plot(sen_res)
    ggsave(p3[[1]], file=paste0("./Result05/",filename2,"/sensitivity-analysis.pdf"), 
           width=8, height=8)
    
    # Funnel plot检验，如果不存在publication bias,那么数据点应该对称分布在两侧,呈漏斗形分布。
    res_single <- mr_singlesnp(dat)
    p4 <- mr_funnel_plot(singlesnp_res)
    ggsave(p4[[1]], file=paste0("./Result05/",filename2,"/funnelplot.pdf"), width=8, height=8)
    # 添加exposure和outcome名称
    res$exposure=ex_pore[i]
    res$outcome=out_come[j]
    
    # 打印分析信息
    print(paste0("------", ex_pore[i], " & ",out_come[j],"------"))
    print(generate_odds_ratios(res))
    
    # 存储结果
    results[[length(out_come)*(i-1)+j]] <- generate_odds_ratios(res)
    }
  }
}

# 绑定所有结果
results_allIV <- do.call(rbind, results) 
fwrite(results_allIV,"result06.csv")

#2.2绘图
#分析前准备
# 使用“ls()”函数获取当前环境中的所有变量、函数和对象的名称，然后传递给“list”参数
# 因此，“rm(list=ls())”表示删除当前环境中的所有变量、函数和对象，相当于清空当前环境
rm(list=ls())
# 判断是否已经安装了“pacman”包，如果没有就安装它
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
# 设置Bioconductor镜像地址为中国科技大学的镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# 加载“pacman”包，用于方便加载其他的R包
library("pacman")
# 使用“p_load”函数加载所需的R包
p_load(forestploter,grid,ggplot2,data.table)
mrresultsfile <- "allOR_Inform06.csv"
# 定义变量mrresultsfile,用于存储MR结果文件的文件名
# 使用read.table()函数读取MR结果文件
# 文件名为mrresultsfile变量的值"MRresults.txt"
# header = T表示文件第一行作为变量名
# sep = "\t" 表示列之间使用制表符分隔
#biofsci = read.table(mrresultsfile, header = T, sep = "\t")
biofsci = fread(mrresultsfile, header = T)
colnames(biofsci)
# 处理P值的显示格式
# 对biofsci数据集的pval字段进行处理
# 使用ifelse()进行条件判断
# 如果pval小于0.001显示"<0.001" 
# 否则使用sprintf()格式化显示P值到小数点后4位
biofsci$pval = as.numeric(biofsci$pval)
biofsci$pval <- ifelse(biofsci$pval<0.05, "<0.05", sprintf("%.4f", biofsci$pval))

# 格式化估计值列，这里可以调整小数点位数，比如把2 改为4
biofsci$estimate <- paste0(format(round(biofsci$or, 4), nsmall = 4), " (",  
                           format(round(biofsci$or_lci95, 4), nsmall = 4), "-",
                           format(round(biofsci$or_uci95, 4), nsmall = 4), ")")

biofsci$Trails <-  biofsci$reportedTrait

# 处理NA值,使用ifelse()对多个字段进行处理
# 如果值为NA,则替换为""

biofsci$Trails = ifelse(is.na(biofsci$Trails), "", biofsci$Trails)
biofsci$method = ifelse(is.na(biofsci$method), "", biofsci$method)
biofsci$nsnp = ifelse(is.na(biofsci$nsnp), "", biofsci$nsnp)  
biofsci$pval = ifelse(is.na(biofsci$pval), "", biofsci$pval)

# 添加空格用于格式调整
# 使用rep()函数重复生成20个空格字符
# 使用paste()用空格拼接成一个字符串
# 赋值给新生成的变量 biofsci` `
colnames(biofsci)
biofsci$` ` <- paste(rep(" ", 15), collapse = " ")   
biofsci$`OR(95%CI)` <- biofsci$estimate
# 提取class. 和 .id 之间的内容
# 仅提取class.和.id之间的内容

colnames(biofsci)
biofsci$Trails[duplicated(biofsci$Trails)] <- " "
biofsci1=biofsci[,c("Trails","method", "nsnp","pval","OR(95%CI)"," ")]
# 处理OR值的显示格式
# 使用ifelse()对OR字段进行处理  
# 如果OR是NA,显示"",否则显示格式化的OR值和95%可信区间
colnames(biofsci1)

# 设置森林图的主题,定义图形元素的显示格式  
# 使用forest_theme()函数
# 定义基础字体大小,参考线颜色,脚注格式等参数                                       
tm <- forest_theme(base_size = 10,
                   refline_col = "red",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")

# 绘制森林图                  
# 使用forest()函数,传递数据及参数
# biofsci[,c(1:5,9:10)]:选择需要的列数据
# est/lower/upper: OR值和置信区间     
# arrow_lab: 曲线两个方向的标签文本  
# sizes: 点的大小
# ci_column: 置信区间所在列
# ref_line: 添加参考线 
# xlim: x轴范围
# footnote: 脚注
# theme: 使用设置的主题对象    
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
rm(list=ls())
  # 判断是否已经安装了“pacman”包，如果没有就安装它
  if(!require("pacman")) install.packages("pacman",update = F,ask = F)
  # 设置Bioconductor镜像地址为中国科技大学的镜像
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
  
  # 如果没有 devtools 包,则安装 devtools 包
  if (!require("devtools")) {
    install.packages("devtools")
  } else {}
  
  # 如果没有 data.table 包,则安装 data.table 包  
  if (!require("data.table")) {
    install.packages("data.table")
  } else {}
  
  # 如果没有 TwoSampleMR 包,则从 GitHub 安装 TwoSampleMR 包
  if (!require("TwoSampleMR")) {
    devtools::install_github("MRCIEU/TwoSampleMR") 
  } else {}
  
  # 从 GitHub 安装 MRInstruments 包
  #remotes::install_github("MRCIEU/MRInstruments")
  # 加载“pacman”包，用于方便加载其他的R包
  library(pacman)
  # 在加载这些包之前，“pacman”会先检查是否已经安装，如果没有会自动安装
  p_load(data.table,dplyr,tidyr,ggplot2)
  
  ##  安装包
  #remotes::install_github("yulab-smu/yulab.utils")
  #library(yulab.utils)
  #install_zip_gh("MRCIEU/ieugwasr")
  # 导入ieugwasr包
  library(ieugwasr)
  library(MRInstruments)
  library(plyr)
  library(dplyr)
  library(data.table)
  library(TwoSampleMR)
  library(ggplot2)
  library(foreach)
  library(purrr)
  

if(!dir.exists("Result10")){ #判断文件夹是否存在
  dir.create("Result10") #创建文件夹
}


#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")

#引用包
library(TwoSampleMR)
ivwFile="ivw68-815.csv"     #IVW方法过滤的结果文件
#读取免疫细胞数据作为结局
biof=read.csv(ivwFile, header=T, sep=",", check.names=F)
outcomeid=biof$id
outcomeid
#outcomeid=outcomeid[1]
#提取暴露数据(疾病的数据)
exposureID=unique(biof$id.outcome)
exposure_dat=fread("GCST90200815_buildGRCh38_clump.txt")
exposure_dat$samplesize= 4488
exposure_dat$PHENO="DESEASE"
fwrite(exposure_dat,"finn.csv")
head(exposure_dat)
# 读取输入文件
exposure_dat <- read_exposure_data("finn.csv", # 读取文件名存储在MRFile中的CSV文件
                              sep = ",",       # 列分隔符为逗号
                              #phenotype_col = "PHENO", # 暴露表型列名
                              snp_col = "variant_id",     # SNP ID列名
                              beta_col = "beta", # 暴露效应大小列名
                              se_col = "standard_error",   # 暴露效应标准误列名
                              effect_allele_col = "effect_allele", # 暴露效应等位基因列名
                              other_allele_col = "other_allele", # 暴露其他等位基因列名
                              pval_col = "p_value",
                              eaf_col = "effect_allele_frequency", # 暴露等位基因频率列名
                              samplesize_col = "samplesize", # 暴露样本大小列名
                              chr_col = "#chrom",
                              pos_col = "pos",
                              clump=F # 进行连锁不平衡剔除
)  # 不进行连锁不平衡剔除



#对结局数据进行循环(免疫细胞)

workdir <- getwd() 
workdir <- paste0(getwd(),"/out")
workdir
# 列出所有以meta开头的txt文件
outcomefile <- list.files(path = workdir, pattern = "^G.*\\.tsv", full.names = TRUE)
# 打印文件全路径 
print(outcomefile)

result=c()
for(i in outcomefile){
  tryCatch({  
    # 开始捕获可能的错误
    ##提取结局数据
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
      # 提取两个数据框的snp列  
      #expo_snps <- exposure_dat$SNP
      #outc_snps <- outc_data$SNP
    #将暴露数据和结局数据合并
    #outc_data$outcome="Colorectal cancer"
    dat <- harmonise_data(exposure_dat, outc_data)
    
    
    #MR-PRESSO异常值检测(偏倚的SNP)
    #presso=run_mr_presso(dat)
    #write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file=paste0(i, ".table.MR-PRESSO_Global.csv"))
    #write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file=paste0(i, ".table.MR-PRESSO_Outlier.csv"))
    
    #孟德尔随机化分析
    mr_result=mr(dat)
    result_or=generate_odds_ratios(mr_result)
    
    #newid <- sub("_.tsv.gz", "", basename(i))
    newid <- gsub("\\.tsv.gz$", "", basename(i))

    filename <- paste0("Result10/", newid) #创建文件夹
    if (!dir.exists(filename)) { #判断文件夹是否存在
      dir.create(filename) #创建文件夹
    }
    
     
    filename2 <- paste0("Result11/OR/", newid) #创建文件夹
   
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
    
    #绘图
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
    
    #输出孟德尔随机化分析的结果
    result=c(result, mr_result$pval[3])
    #对IVW方法pvalue小于0.05的结果可视化
   # if (!is.na(result_or$pval[3]) && result_or$pval[3] < 0.05) {
      if (T) {
      #输出用于孟德尔随机化的工具变量
      outTab=dat[dat$mr_keep=="TRUE",]
    }
  }, error = function(e) {
    # 处理错误，这里打印错误信息
    cat("Error occurred for file:", i, "\n")
    cat("Error message:", conditionMessage(e), "\n")
  })
}

#3.1计算中介效应值
rm(list=ls())
# 判断是否已经安装了“pacman”包，如果没有就安装它
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
# 设置Bioconductor镜像地址为中国科技大学的镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
# 加载“pacman”包，用于方便加载其他的R包
library("pacman")
# 在加载这些包之前，“pacman”会先检查是否已经安装，如果没有会自动安装
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
  
  # 使用 strsplit 函数分割字符串
  file_parts <- strsplit(basename(beta1File), split = "_")
  
  # 获取第二部分beta2文件名
  second_part <- file_parts[[1]][2]
  beta2File=paste0("ORdata19_beta2/",second_part,"_OR_beta2.csv")
  
  # 获取第一部分betaall文件名
  first_part <- file_parts[[1]][1]
  betaAllFile=paste0("ORdata20_betaAll/",first_part,"_OR_betaAll.csv")
  
  # 创建输出文件夹，如果不存在的话
  output_folder <- "medResult21"
  if (!file.exists(output_folder)) {
    dir.create(output_folder)
  }
  
  if (file.exists(beta2File) & file.exists(betaAllFile)) {
    
    #beta1File="ebi-a-GCST90001662_GCST90200850_OR_beta1.csv"     #免疫细胞-代谢物MR分析结果OR_beta1文件
    #beta2File="GCST90200850_OR_beta2.csv"                   #代谢物-疾病MR分析结果OR_beta2文件
    #betaAllFile="ebi-a-GCST90001662_OR_betaAll.csv"              #免疫细胞-疾病MR分析结果OR_betaAll文件
    
    
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
    
    
    
    # 定义字符串
    my_string <- basename(beta1File)
    # 使用sub函数将_OR_beta1.csv替换为空
    new_string <- sub("_OR_beta1.csv", "", my_string)
    # 在最后，将结果写入新的CSV文件
    
    # 写结果到新的文件夹中
    write.csv(result, file = paste0(output_folder, "/", sub("_OR_beta1.csv", "", basename(beta1File)), "_Med_mr_Effect21", ".csv"), row.names=F)
  } else {
    # 如果文件不存在，就跳到下一次循环
    warning(paste("Files", beta2File, "and/or", betaAllFile, "do not exist. Skipping this file."))
    next
  }
}

#3.2绘图
rm(list=ls())
# 判断是否已经安装了“pacman”包，如果没有就安装它
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
# 设置Bioconductor镜像地址为中国科技大学的镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
# 加载“pacman”包，用于方便加载其他的R包
library("pacman")
# 在加载这些包之前，“pacman”会先检查是否已经安装，如果没有会自动安装
p_load(grid,readr,forestploter)

# 引用所需包
library(grid)
library(readr)
library(forestploter)

# 数据整理
biofsci <- do.call(rbind, lapply(grep(".csv$", dir(), value = TRUE),
                              function(i) read.csv(i, header = T, sep = ",", check.names = F)))

lineVec <- cumsum(c(1, table(biofsci[, c('exposureTrait', 'outcomeTrait')])))
biofsci$' ' <- paste(rep(" ", 10), collapse = " ")
biofsci$'OR(95% CI)' <- ifelse(is.na(biofsci$or), "", sprintf("%.3f (%.3f to %.3f)", biofsci$or, biofsci$or_lci95, biofsci$or_uci95))
biofsci$pval <- ifelse(biofsci$pval < 0.001, "<0.001", sprintf("%.3f", biofsci$pval))
biofsci$exposure <- ifelse(is.na(biofsci$exposure), "", biofsci$exposure)
biofsci$nsnp <- ifelse(is.na(biofsci$nsnp), "", biofsci$nsnp)
biofsci[duplicated(biofsci[, c('exposureTrait','outcomeTrait')]), c('exposureTrait', 'outcomeTrait')] <- ""

# 主题参数设置
tm <- forestploter::forest_theme(base_size = 15, ci_pch = 15, ci_lty = 1, ci_lwd = 1.5,
                                 ci_col = "black", ci_Theight = 0.2,
                                 refline_lty = "dashed", refline_lwd = 1, 
                                 refline_col = "grey40", xaxis_cex = 0.8, 
                                 footnote_cex = 0.6, footnote_col = "darkblue") 

# 生成并修改图形
plot <- forestploter::forest(biofsci[, c("exposureTrait","outcomeTrait","nsnp","method","pval", " ","OR(95% CI)")],
                             est = biofsci$or, lower = biofsci$or_lci95, upper = biofsci$or_uci95, 
                             ci_column = 6, ref_line = 1, xlim = c(0.5, 1.5), theme = tm)

# 定义颜色
nature_colors <- c("#009E73", "#56B4E9", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# 应用颜色到你的数据
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

# 输出图形
pdf("mr_forest.pdf", width = 28, height = 8)
print(plot)
dev.off()

