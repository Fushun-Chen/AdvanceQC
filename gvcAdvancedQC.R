# verion1
# author: chenfushun
# the purpose of program compare simp file with TCGA and COSMIC database. input: simp files, cancer name, depth num, code dir. output: comparison file and figure.
library(data.table)
library(ggplot2)
library(stats)
library(RColorBrewer)
library(VennDiagram)
library(Cairo)
# Function 
# Base Substitution 
BaseSubNum <- function(base.sub){
  # Computes the base substitution num
  # Args : base.sub is type of base substitution.
  # Returns : the number of base substitution
  base.sub.freq <- data.frame(table(base.sub))
  G_A_C_T =sum(base.sub.freq$Freq[base.sub.freq$base.sub %in% c("G>A","C>T")])
  G_T_C_A =sum(base.sub.freq$Freq[base.sub.freq$base.sub %in% c("G>T","C>A")])
  A_G_T_C =sum(base.sub.freq$Freq[base.sub.freq$base.sub %in% c("A>G","T>C")])
  G_C_C_G =sum(base.sub.freq$Freq[base.sub.freq$base.sub %in% c("G>C","C>G")])
  A_C_T_G =sum(base.sub.freq$Freq[base.sub.freq$base.sub %in% c("A>C","T>G")])
  A_T_T_A =sum(base.sub.freq$Freq[base.sub.freq$base.sub %in% c("A>T","T>A")])

  base.sub.type.num <- data.frame(base_type = c("G>A+C>T", "G>T+C>A", "A>G+T>C", "G>C+C>G", "A>C+T>G", "A>T+T>A"), base_type_num = c(G_A_C_T, G_T_C_A, A_G_T_C, G_C_C_G, A_C_T_G, A_T_T_A)) 
  return(base.sub.type.num) 
}
# pnorm 
Pnorm <- function(x, mean.value, sd.value){
  if(x > mean.value){
    pnorm.value <- 1 - pnorm(x, mean = mean.value, sd = sd.value)
  }else{
    pnorm.value <- pnorm(x, mean = mean.value, sd = sd.value)
  }
  return(pnorm.value)
}
# 
PatientCompareMC3tumor <- function(tumor.mc3.exonic.depth.value, patient.exonic.value, output.file){
  # mut_depth
  tumor.mc3.exonic.depth.mut.depth <- data.frame(Tumor_Sample_Barcode = tumor.mc3.exonic.depth.value$Tumor_Sample_Barcode, Depth = tumor.mc3.exonic.depth.value$t_depth)
  depth.pvalue <- format(as.numeric(ks.test(tumor.mc3.exonic.depth.mut.depth$Depth, patient.exonic.value$Depth)$p.value),digits=3)

  pdf(paste(output.file,".depth.pdf", sep = ""))
  p <- ggplot(tumor.mc3.exonic.depth.mut.depth, aes(x=Depth)) +  geom_density(size = 1.2, color =  brewer.pal(8,"Set1")[2]) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15))
  p1 <- p + geom_density(data = patient.exonic.value, size = 1.2, color =  brewer.pal(8,"Set1")[1])  + annotate("text", x=Inf, y = Inf, label = paste("pvalue:",depth.pvalue, sep = "") , vjust=2, hjust=1.5) + xlab("Depth")
  print(p1)
  dev.off()

  # mut_Rate
  tumor.mc3.exonic.depth.t.rate <- data.frame(Tumor_Sample_Barcode = tumor.mc3.exonic.depth.value$Tumor_Sample_Barcode, MutAF_tumor = (tumor.mc3.exonic.depth.value$t_alt_count/tumor.mc3.exonic.depth.value$t_depth))
  mut.rate.pvalue <- format(as.numeric(ks.test(tumor.mc3.exonic.depth.t.rate$MutAF_tumor, patient.exonic.value$MutAF_tumor)$p.value),digits=3)

  pdf(paste(output.file,".mut_rate.pdf", sep = ""))
  p <- ggplot(tumor.mc3.exonic.depth.t.rate, aes(x=MutAF_tumor)) +  geom_density(size = 1.2, color =  brewer.pal(8,"Set1")[2]) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15)) + xlab("Mut Rate")
  p1 <- p + geom_density(data = patient.exonic.value,  size = 1.2, color =  brewer.pal(8,"Set1")[1])  + annotate("text", x=Inf, y = Inf, label = paste("pvalue:",mut.rate.pvalue, sep = ""), vjust=2, hjust=1.5)
  print(p1)
  dev.off()

  # MutNum, 1KG, dbSNP rate, AS-ratio
  tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value <- tumor.mc3.exonic.depth.value[,.( Mut_num = length(t_depth), MC3_1KG_MAF_ALL_rate = c(length(`1KG_MAF_ALL`) - sum(`1KG_MAF_ALL` == "-")) / length(`1KG_MAF_ALL`) , rsID_rate = c(length(rsID) - sum(rsID == "-")) / length(rsID), nosy_sy_rate = sum(Mutation_type %in% c("nonsynonymous SNV", "stopgain")) / sum(Mutation_type == "synonymous SNV")),by = Tumor_Sample_Barcode]
  tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$nosy_sy_rate[is.infinite(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$nosy_sy_rate)] <- NA
  # MutNum
  patient.exonic.mut.num <- nrow(patient.exonic.value)

  patient.mut.greater.TCGA.rate <- sum(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$Mut_num > patient.exonic.mut.num) / nrow(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value)
  patient.mut.less.TCGA.rate <- sum(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$Mut_num <= patient.exonic.mut.num) / nrow(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value)
  patient.exonic.mut.num.pnorm <- Pnorm(patient.exonic.mut.num, mean(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$Mut_num, na.rm = TRUE), sd(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$Mut_num, na.rm = TRUE))
  pdf(paste(output.file,".mut_num_density.pdf", sep = ""))
  p <- ggplot(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value, aes(x=Mut_num)) +  geom_density(size = 1.2, color =  brewer.pal(8,"Set1")[2]) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15)) + xlab("Mutation Num")
  p1 <- p + geom_vline(xintercept = c(as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$Mut_num,0.05)), as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$Mut_num,0.95))), linetype="dotted", color = "red", size=1.2) + annotate("text", x = as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$Mut_num,0.05)), y = Inf, label = "5%", vjust = 2,hjust = 1.1) + annotate("text", x = as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$Mut_num,0.95)), y = Inf, label = "95%", vjust = 2,hjust = 1.1)
  p2 <- p1 + geom_vline(xintercept = patient.exonic.mut.num, size=1.2, colour=brewer.pal(8,"Set1")[1]) 
  print(p2)
  dev.off()

  # 1KG_MAF_ALL
  patient.1KG.MAF.ALL.Num.rate <- c(length(patient.exonic.value$"1KG_MAF_ALL") - sum(patient.exonic.value$"1KG_MAF_ALL" == "-")) / length(patient.exonic.value$"1KG_MAF_ALL")
  patient.1KG.MAF.ALL.mean.rate <-  mean(c(as.numeric(patient.exonic.value$"1KG_MAF_ALL"[patient.exonic.value$"1KG_MAF_ALL" != "-"])))

  patient.1KG.rate.greater.TCGA.rate <- sum(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$MC3_1KG_MAF_ALL_rate > patient.1KG.MAF.ALL.Num.rate) / nrow(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value)
  patient.1KG.rate.less.TCGA.rate <- sum(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$MC3_1KG_MAF_ALL_rate <= patient.1KG.MAF.ALL.Num.rate) / nrow(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value)
  
  patient.1KG.MAF.ALL.Num.rate.pnorm <- Pnorm(patient.1KG.MAF.ALL.Num.rate, mean(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$MC3_1KG_MAF_ALL_rate, na.rm = TRUE), sd(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$MC3_1KG_MAF_ALL_rate, na.rm = TRUE))
  pdf(paste(output.file,".1KG_rate_density.pdf", sep = ""))
  p <- ggplot(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value, aes(x=MC3_1KG_MAF_ALL_rate)) +  geom_density(size = 1.2, color =  brewer.pal(8,"Set1")[2]) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15)) + xlab("1KG MAF rate")
  p1 <- p + geom_vline(xintercept = c(as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$MC3_1KG_MAF_ALL_rate,0.05)), as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$MC3_1KG_MAF_ALL_rate,0.95))), linetype="dotted", color = "red", size=1.2) + annotate("text", x = as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$MC3_1KG_MAF_ALL_rate,0.05)), y = Inf, label = "5%", vjust = 2,hjust = 1.1) + annotate("text", x = as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$MC3_1KG_MAF_ALL_rate,0.95)), y = Inf, label = "95%", vjust = 2,hjust = 1.1) 
  p2 <- p1 + geom_vline(xintercept = patient.1KG.MAF.ALL.Num.rate, size=1.2, colour=brewer.pal(8,"Set1")[1]) 
  print(p2)
  dev.off()

  # rsID rate 
  patient.rsID.rate <- c(length(patient.exonic.value$rsID) - sum(patient.exonic.value$rsID == "-")) / length(patient.exonic.value$rsID)
  patient.rsID.rate.greater.TCGA.rate <- sum(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$rsID_rate > patient.rsID.rate) / nrow(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value)
  patient.rsID.rate.less.TCGA.rate <- sum(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$rsID_rate <= patient.rsID.rate) / nrow(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value)
  patient.rsID.rate.pnorm <- Pnorm(patient.rsID.rate, mean(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$rsID_rate, na.rm = TRUE), sd(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$rsID_rate, na.rm = TRUE))
  
  pdf(paste(output.file,".rsID_rate_density.pdf", sep = ""))
  p <- ggplot(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value, aes(x=rsID_rate)) +  geom_density(size = 1.2, color =  brewer.pal(8,"Set1")[2]) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15)) + xlab("rsID rate")
  p1 <- p + geom_vline(xintercept = c(as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$rsID_rate,0.05)), as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$rsID_rate,0.95))), linetype="dotted", color = "red", size=1.2) + annotate("text", x = as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$rsID_rate,0.05)), y = Inf, label = "5%", vjust = 2,hjust = 1.1) + annotate("text", x = as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$rsID_rate,0.95)), y = Inf, label = "95%", vjust = 2,hjust = 1.1) 
  p2 <- p1 + geom_vline(xintercept = patient.rsID.rate, size=1.2, colour=brewer.pal(8,"Set1")[1]) 
  print(p2)
  dev.off()

  # nosy_sy_rate 
  # patient NoSy_num Sy_num
  patient.nosy.num <- sum(patient.exonic.value$Mutation_type %in% c("nonsynonymous SNV", "stopgain")) 
  patient.sy.num <-  sum(patient.exonic.value$Mutation_type == "synonymous SNV")
  # patient nosy/sy in dbsnp
  patient.exonic.dbsnp.value <- patient.exonic.value[patient.exonic.value$rsID != "-",]
  patient.nosy.sy.indbsnp.rate <- sum(patient.exonic.dbsnp.value$Mutation_type %in% c("nonsynonymous SNV", "stopgain")) / sum(patient.exonic.dbsnp.value$Mutation_type == "synonymous SNV")
  # patient nosy/sy in nodbsnp
  patient.exonic.nodbsnp.value <- patient.exonic.value[patient.exonic.value$rsID == "-",]
  patient.nosy.sy.nodbsnp.rate <- sum(patient.exonic.nodbsnp.value$Mutation_type %in% c("nonsynonymous SNV", "stopgain")) / sum(patient.exonic.nodbsnp.value$Mutation_type == "synonymous SNV")

  # patient As_ratio
  patient.as.ratio <- patient.nosy.num / patient.sy.num

  patient.nosy.sy.rate.greater.TCGA.rate <- sum(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$nosy_sy_rate > patient.as.ratio, na.rm = TRUE) / nrow(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value) 
  patient.nosy.sy.rate.less.TCGA.rate <- sum(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$nosy_sy_rate <= patient.as.ratio, na.rm = TRUE) / nrow(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value) 

  patient.as.ratio.pnorm <- Pnorm(patient.as.ratio, mean(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$nosy_sy_rate, na.rm = TRUE), sd(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$nosy_sy_rate, na.rm = TRUE))
  
  pdf(paste(output.file,".nosy_sy_rate_density.pdf", sep = ""))
  p <- ggplot(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value, aes(x=nosy_sy_rate)) +  geom_density(size = 1.2, color =  brewer.pal(8,"Set1")[2]) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15)) + xlab("AS-Ratio")
  p1 <- p + geom_vline(xintercept = c(as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$nosy_sy_rate, 0.05, na.rm = T)), as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$nosy_sy_rate, 0.95, na.rm = T))), linetype="dotted", color = "red", size=1.2) + annotate("text", x = as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$nosy_sy_rate, 0.05, na.rm = T)), y = Inf, label = "5%", vjust = 2,hjust = 1.1) + annotate("text", x = as.numeric(quantile(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$nosy_sy_rate, 0.95, na.rm = T)), y = Inf, label = "95%", vjust = 2,hjust = 1.1) 
  p2 <- p1 + geom_vline(xintercept = patient.as.ratio, size=1.2, colour=brewer.pal(8,"Set1")[1]) 
  print(p2)
  dev.off()

  # Base Substitution 
  # tumor mc3 base substitution
  tumor.mc3.exonic.depth.value.base.substitution <- data.frame(base_substitution = paste(tumor.mc3.exonic.depth.value$Reference_Allele, tumor.mc3.exonic.depth.value$Tumor_Seq_Allele2, sep = ">"))
  tumor.mc3.exonic.depth.value.base.sub.num <- BaseSubNum(tumor.mc3.exonic.depth.value.base.substitution$base_substitution)
  # patient base sub 
  patient.exonic.value.base.substitution <- data.frame(base_substitution = paste(patient.exonic.value$Ref, patient.exonic.value$Mut, sep = ">"))
  patient.exonic.value.base.sub.num <- BaseSubNum(patient.exonic.value.base.substitution)

  # the correlation of tumor base sub and patient base sub
  tumor.mc3.exonic.depth.value.base.sub.num$type <- paste("MC3",cancer.name,sep = "-")
  tumor.mc3.exonic.depth.value.base.sub.num$base_type_rate <- tumor.mc3.exonic.depth.value.base.sub.num$base_type_num / sum(tumor.mc3.exonic.depth.value.base.sub.num$base_type_num)
  patient.exonic.value.base.sub.num$type <- "Patient"
  patient.exonic.value.base.sub.num$base_type_rate <- patient.exonic.value.base.sub.num$base_type_num / sum(patient.exonic.value.base.sub.num$base_type_num)
  mc3.patient.base.sub.cor <- format(cor(tumor.mc3.exonic.depth.value.base.sub.num$base_type_rate, patient.exonic.value.base.sub.num$base_type_rate, method = "kendall"), digits = 3)
  mc3.patient.base.sub.value <- rbind(tumor.mc3.exonic.depth.value.base.sub.num, patient.exonic.value.base.sub.num)
  mc3.patient.base.sub.value$base_type <- factor(mc3.patient.base.sub.value$base_type, levels = c("G>A+C>T", "G>T+C>A", "A>G+T>C", "G>C+C>G", "A>C+T>G", "A>T+T>A"))

  pdf(paste(output.file,".base_sub.pdf", sep = ""))
  p0 <- ggplot(mc3.patient.base.sub.value, aes(x = type, y = base_type_rate, fill = base_type)) + geom_bar(stat="identity") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15)) + scale_fill_manual(values=c('DarkRed','orange','yellow','SteelBlue1','blue','SteelBlue4')) + annotate("text", x = 2.3, y = 1.05, label = paste("cor:",mc3.patient.base.sub.cor, sep = ""),size = 5)
  print(p0)
  dev.off()
  # return type num
  df.result <- data.frame(TCGAvsPatient_depth_pvalue = depth.pvalue, TCGAvsPatient_mutRate_pvalue = mut.rate.pvalue, Patient_exonic_mut_num = patient.exonic.mut.num, Patient_exonic_mutNum_greater_TCGA_rate = patient.mut.greater.TCGA.rate, Patient_exonic_mutNum_less_TCGA_rate = patient.mut.less.TCGA.rate, Patient_exonic_mut_num_pnorm = patient.exonic.mut.num.pnorm, Patient_1KG_MAF_ALL_Num_rate = patient.1KG.MAF.ALL.Num.rate, Patient_1KG_MAF_ALL_mean_rate = patient.1KG.MAF.ALL.mean.rate, Patient_1KG_rate_greater_TCGA_rate = patient.1KG.rate.greater.TCGA.rate, Patient_1KG_rate_less_TCGA_rate = patient.1KG.rate.less.TCGA.rate, Patient_1KG_MAF_ALL_Num_rate_pnorm = patient.1KG.MAF.ALL.Num.rate.pnorm, Patient_rsID_rate = patient.rsID.rate, Patient_rsID_rate_greater_TCGA_rate = patient.rsID.rate.greater.TCGA.rate, Patient_rsID_rate_less_TCGA_rate = patient.rsID.rate.less.TCGA.rate, Patient_rsID_rate_pnorm = patient.rsID.rate.pnorm, Patient_NoSy_num = patient.nosy.num, Patient_Sy_num = patient.sy.num, Patient_nosy_sy_rate = patient.as.ratio, Patient_nosy_sy_rate_greater_TCGA_rate = patient.nosy.sy.rate.greater.TCGA.rate, Patient_nosy_sy_rate_less_TCGA_rate = patient.nosy.sy.rate.less.TCGA.rate, Patient_as_ratio_pnorm = patient.as.ratio.pnorm, Patient_nosy_sy_rate_dbsnp = patient.nosy.sy.indbsnp.rate, Patient_nosy_sy_rate_nodbsnp = patient.nosy.sy.nodbsnp.rate, MC3_patient_baseSub_cor = mc3.patient.base.sub.cor)
  return(df.result)

}

# the mut number of patient in CosMic and TCGA
PatientInCosMicMC3MutAANum <- function(patient.exonic.value, cosmic.value, tcga.value){

  # the number of patient chr_pos in CosMic and TCGA
  patient.mut.pos.in.cosmic.tcga.num <- sum(paste(patient.exonic.value$Chromosome,patient.exonic.value$Position,sep = "_") %in% c(paste(cosmic.value$Chr,cosmic.value$Pos, sep = "_"), paste(tcga.value$Chromosome,tcga.value$Start_Position, sep = "_")))
  # the number of patient chr_pos_type in CosMic and TCGA
  patient.mut.pos.type.in.cosmic.tcga.num <- sum(paste(patient.exonic.value$Chromosome,patient.exonic.value$Position, patient.exonic.value$Ref,patient.exonic.value$Mut, sep = "_") %in% c(paste(cosmic.value$Chr,cosmic.value$Pos, cosmic.value$Ref,cosmic.value$Mut, sep = "_"), paste(tcga.value$Chromosome,tcga.value$Start_Position, tcga.value$Reference_Allele,tcga.value$Tumor_Seq_Allele2, sep = "_")))
  
  # the number of patient gene_AApos in CosMic and TCGA
  patient.AA.pos.in.cosmic.tcga.num <- sum(paste(patient.exonic.value$Gene_name, substr(patient.exonic.value$AA_change,2,nchar(patient.exonic.value$AA_change) - 1), sep = "_") %in% c(paste(cosmic.value$Gene_name, substr(cosmic.value$AA_change,2, nchar(cosmic.value$AA_change) - 1), sep = "_"), paste(tcga.value$Hugo_Symbol, substr(tcga.value$AA_change,2, nchar(tcga.value$AA_change) - 1), sep = "_")))
  # the number of patient gene_AApos_type in CosMic and TCGA
  patient.AA.pos.type.in.cosmic.tcga.num <- sum(paste(patient.exonic.value$Gene_name, patient.exonic.value$AA_change,sep = "_") %in% c(paste(cosmic.value$Gene_name, cosmic.value$AA_change, sep = "_"), paste(tcga.value$Hugo_Symbol, tcga.value$AA_change, sep = "_")))
  
  #data frame result 
  df.result <- data.frame(patient_mut_pos_in_CosMic_TCGA_Num = patient.mut.pos.in.cosmic.tcga.num, patient_mut_pos_type_in_CosMic_TCGA_Num = patient.mut.pos.type.in.cosmic.tcga.num, patient_AA_pos_in_CosMic_TCGA_Num = patient.AA.pos.in.cosmic.tcga.num, patient_AA_pos_type_in_CosMic_TCGA_Num = patient.AA.pos.type.in.cosmic.tcga.num)
  return(df.result)
}

# the mut number of patient in tumor
PatientInTumorMutAANum <- function(patient.exonic.value, tumor.value){
  tumor.name <- paste(unique(tumor.value$Cancer_name), collapse = "_")
  # the number of patient chr_pos in tumor
  patient.mut.pos.in.tumor.num <- sum(paste(patient.exonic.value$Chromosome,patient.exonic.value$Position,sep = "_") %in% c(paste(tumor.value$Chromosome, tumor.value$Start_Position, sep = "_")))
  # the number of patient chr_pos_type in tumor
  patient.mut.pos.type.in.tumor.num <- sum(paste(patient.exonic.value$Chromosome,patient.exonic.value$Position, patient.exonic.value$Ref,patient.exonic.value$Mut, sep = "_") %in% c(paste(tumor.value$Chromosome,tumor.value$Start_Position, tumor.value$Reference_Allele,tumor.value$Tumor_Seq_Allele2, sep = "_")))
  
  # the number of patient gene_AApos in tumor
  patient.AA.pos.in.tumor.num <- sum(paste(patient.exonic.value$Gene_name, substr(patient.exonic.value$AA_change,2,nchar(patient.exonic.value$AA_change) - 1), sep = "_") %in% c(paste(tumor.value$Hugo_Symbol, substr(tumor.value$AA_change,2, nchar(tumor.value$AA_change) - 1), sep = "_")))
  # the number of patient gene_AApos_type in tumor
  patient.AA.pos.type.in.tumor.num <- sum(paste(patient.exonic.value$Gene_name, patient.exonic.value$AA_change,sep = "_") %in% c(paste(tumor.value$Hugo_Symbol, tumor.value$AA_change, sep = "_")))
  
  #data frame result 
  df.result <- data.frame(patient_mut_pos_in = patient.mut.pos.in.tumor.num, patient_mut_pos_type_in = patient.mut.pos.type.in.tumor.num, patient_AA_pos_in = patient.AA.pos.in.tumor.num, patient_AA_pos_type_in = patient.AA.pos.type.in.tumor.num)
  colnames(df.result) <- paste(colnames(df.result), tumor.name, "Num", sep = "_")
  return(df.result)
}

# the rate of Patient mut in DataBase(tcga or cosmic)
PatientMutDataBaseRate <- function(patient.exonic.value, mut.pos.gene.aa.value){

  colnames(mut.pos.gene.aa.value) <- c("Chr", "Pos", "Ref", "Mut", "Gene_name", "AA_change", "Barcode")
  # mut chr pos
  mut.chr.pos <- paste(mut.pos.gene.aa.value$Chr, mut.pos.gene.aa.value$Pos, sep = "_")[paste(mut.pos.gene.aa.value$Chr, mut.pos.gene.aa.value$Pos, sep = "_") %in% paste(patient.exonic.value$Chromosome, patient.exonic.value$Position, sep = "_")] 
  df.mut.chr.pos.freq <- unique(data.frame(table(mut.chr.pos)))
  mut.chr.pos.rate <- mean(df.mut.chr.pos.freq$Freq) / length(unique(mut.pos.gene.aa.value$Barcode))
  # mut chr pos type
  mut.chr.pos.type <- paste(mut.pos.gene.aa.value$Chr, mut.pos.gene.aa.value$Pos, mut.pos.gene.aa.value$Ref, mut.pos.gene.aa.value$Mut, sep = "_")[paste(mut.pos.gene.aa.value$Chr, mut.pos.gene.aa.value$Pos, mut.pos.gene.aa.value$Ref, mut.pos.gene.aa.value$Mut, sep = "_") %in% paste(patient.exonic.value$Chromosome, patient.exonic.value$Position, patient.exonic.value$Ref, patient.exonic.value$Mut, sep = "_")] 
  df.mut.chr.pos.type.freq <- unique(data.frame(table(mut.chr.pos.type)))
  mut.chr.pos.type.rate <- mean(df.mut.chr.pos.type.freq$Freq) / length(unique(mut.pos.gene.aa.value$Barcode))
  # mut gene AA 
  mut.gene.AA <- paste(mut.pos.gene.aa.value$Gene_name, substr(mut.pos.gene.aa.value$AA_change,2, nchar(mut.pos.gene.aa.value$AA_change) - 1), sep = "_")[paste(mut.pos.gene.aa.value$Gene_name, substr(mut.pos.gene.aa.value$AA_change,2, nchar(mut.pos.gene.aa.value$AA_change) - 1), sep = "_") %in% paste(patient.exonic.value$Gene_name, substr(patient.exonic.value$AA_change,2, nchar(patient.exonic.value$AA_change) - 1), sep = "_")] 
  df.mut.gene.AA.freq <- unique(data.frame(table(mut.gene.AA)))
  mut.gene.AA.rate <- mean(df.mut.gene.AA.freq$Freq) / length(unique(mut.pos.gene.aa.value$Barcode))
  # mut gene AA type
  mut.gene.AA.type <- paste(mut.pos.gene.aa.value$Gene_name, mut.pos.gene.aa.value$AA_change, sep = "_")[paste(mut.pos.gene.aa.value$Gene_name, mut.pos.gene.aa.value$AA_change, sep = "_") %in% paste(patient.exonic.value$Gene_name, patient.exonic.value$AA_change, sep = "_")] 
  df.mut.gene.AA.type.freq <- unique(data.frame(table(mut.gene.AA.type)))
  mut.gene.AA.type.rate <- mean(df.mut.gene.AA.type.freq$Freq) / length(unique(mut.pos.gene.aa.value$Barcode))
  #data frame result 
  df.result <- data.frame(patient_chr_pos_mean_rate = mut.chr.pos.rate, patient_chr_pos_type_mean_rate = mut.chr.pos.type.rate, patient_chr_AA_mean_rate = mut.gene.AA.rate, patient_chr_AA_type_mean_rate = mut.gene.AA.type.rate)
  return(df.result)

}

# Venn figure function
VennFigure <- function(setlist, sample.num, filename, output.dir){
  if(sample.num == 2){
    xname <- list(
      Set1=setlist[[1]]$A,
      Set2=setlist[[2]]$A
    )
    fil=c('red','blue')
    color=c('darkred','white','darkblue')
    darkfil=c('darkred','darkblue')
    filenames=filename
  }else if(sample.num == 3){
    xname <- list(
      Set1=setlist[[1]]$A,
      Set2=setlist[[2]]$A,
      Set3=setlist[[3]]$A)
    fil<-c("red", "blue", "green")
    color<-c("darkred", "white", "darkblue", "white",
             "white", "white", "darkgreen")
    darkfil<-c("darkred", "darkblue", "darkgreen")
    filenames=filename
  }else if(sample.num == 4){
    xname <- list(
      Set1=setlist[[1]]$A,
      Set2=setlist[[2]]$A,
      Set3=setlist[[3]]$A,
      Set4=setlist[[4]]$A)
    fil<-c("cornflowerblue", "green", "yellow", "darkorchid1")
    color<- c("orange", "white", "darkorchid4", "white",
              "white", "white", "white", "white", "darkblue", "white",
              "white", "white", "white", "darkgreen", "white")
    darkfil<- c("darkblue", "darkgreen", "orange", "darkorchid4")
    filenames=filename
  }else if(sample.num == 5){
    xname=list(
      Set1=setlist[[1]]$A,
      Set2=setlist[[2]]$A,
      Set3=setlist[[3]]$A,
      Set4=setlist[[4]]$A,
      Set5=setlist[[5]]$A
    )
    fil=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
    color='black'
    darkfil=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
    filenames=filename
  }

  venn.plot<-venn.diagram(
    x = xname,
    filename = NULL,
    category.names = filenames,
    col = "transparent",
    fill = fil,
    alpha = 0.5,
    label.col = color,
    cex = 1.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.default.pos = "text",
    cat.col = darkfil,
    cat.cex = 1.5,
    cat.fontfamily = "serif",
    cat.dist = 0.07,
    cat.pos = 0
  )
  CairoPNG(file = paste(output.dir, "/", "Venn.png", sep = ""))
  grid.draw(venn.plot)
  dev.off()
}


# args
args<-commandArgs(TRUE)
if(length(args) < 7){
  stop('intFile answerNum codeDir depthNum cancerName sampleName outDir')
}
# 

answer.num <- args[length(args) - 5]
code.dir <- args[length(args) - 4]
depth.num <- args[length(args) - 3]
cancer.name <- args[length(args) - 2]
sample.all.names <- args[length(args) - 1]
output.dir <- args[length(args)]

#
sample.names <- strsplit(sample.all.names,',',fixed = F)[[1]]
###
#code.dir <- "/disk/fushun/Project/GVCAdvancedQC_5_15/code/"
#patient.file <- "/disk/fushun/Project/GVC_train/GVCQC_beirui_5_5_simp/data/BR0504.empirical.Somatic.WES.snv.simp"
#depth.num <- 100
#cancer.name <- "LUAD"
#output.file <- "/disk/fushun/Project/GVCAdvancedQC_5_15/test/BR0504"
##

# reading tcga mc3 value, mc3 value in exonic , getting tumor value, getting mean depth 
mc3.value <- fread(paste(code.dir,"/","tcga_mc3_SNP_sort_ann_new_apart.txt", sep = ""))
mc3.exonic.value <- mc3.value[mc3.value$Mutation_type != ".",]
#
if(!(cancer.name %in% unique(mc3.value$Cancer_name))){
  stop(paste("TCGA has not the cancer name : ", cancer.name,sep = ""))

}

tumor.mc3.exonic.value <- mc3.exonic.value[mc3.exonic.value$Cancer_name %in% cancer.name,]
# tumor mean depth
tumor.meandepth <- tumor.mc3.exonic.value[,.(mean_depth = mean(t_depth), mutNum = length(t_depth)), by = Tumor_Sample_Barcode]
tumor.meandepth$mean_depth_format <- 0
tumor.meandepth$mean_depth_format[tumor.meandepth$mean_depth < 100] <- as.numeric(format(tumor.meandepth$mean_depth[tumor.meandepth$mean_depth < 100], scientific = T, digits = 1))
tumor.meandepth$mean_depth_format[tumor.meandepth$mean_depth >= 100] <- as.numeric(format(tumor.meandepth$mean_depth[tumor.meandepth$mean_depth >= 100], scientific = T, digits = 2))
# getting tumor value in depth num
depth.num.barcode <- tumor.meandepth$Tumor_Sample_Barcode[tumor.meandepth$mean_depth_format == depth.num]
tumor.mc3.exonic.depth.value <- tumor.mc3.exonic.value[tumor.mc3.exonic.value$Tumor_Sample_Barcode %in% depth.num.barcode,]
# MC3 synoy and Nosy 
nosy.mc3.value <- mc3.value[mc3.value$Mutation_type %in% c("nonsynonymous SNV", "stopgain"), ]
sy.mc3.value <- mc3.value[mc3.value$Mutation_type == "synonymous SNV", ]
# MC3 synoy and Nosy in Tumor
nosy.tumor.mc3.value <- nosy.mc3.value[nosy.mc3.value$Cancer_name == cancer.name,]
sy.tumor.mc3.value <- sy.mc3.value[sy.mc3.value$Cancer_name == cancer.name,]
# Cosmic value
cosmic.value <- fread(paste(code.dir,"/","chr_pos_cosmic_value_Substition_noheader_sort_ann_apart_peopleNum.txt", sep = ""))
# cosmic synoy and Nosy
nosy.cosmic.value <- cosmic.value[cosmic.value$Mutation_type %in% c("nonsynonymous SNV", "stopgain"),]
sy.cosmic.value <- cosmic.value[cosmic.value$Mutation_type == "synonymous SNV",]

# 
all.patient.result.list <- list()
for(i in 1:(length(args) - 6)){
  patient.file <- args[i]
  sample.name <- sample.names[i]
  output.file <- paste(output.dir, "/", sample.name, sep = "")
  # patient value
  patient.value <- fread(patient.file)
  patient.exonic.value <- patient.value[patient.value$Function_region == "exonic",]
  patient.mut.num <- nrow(patient.value)
  patient.nosy.value <- patient.exonic.value[patient.exonic.value$Mutation_type %in% c("nonsynonymous SNV", "stopgain"),]
  patient.sy.value <-  patient.exonic.value[patient.exonic.value$Mutation_type == "synonymous SNV",]
  # compare patient type with mc3 tumor
  patient.mc3tumor.exonic.mutnum.mutdepth.1KG.rsID.AsRatio.baseSub <- PatientCompareMC3tumor(tumor.mc3.exonic.depth.value, patient.exonic.value, output.file)
  ### patient pos in CosMiC, TCGA(MC3) , Tumor(MC3) 
  # the mut number of patient in CosMic and TCGA
  nosy.patient.cosmic.mc3.num <- PatientInCosMicMC3MutAANum(patient.nosy.value, nosy.cosmic.value,nosy.mc3.value)
  colnames(nosy.patient.cosmic.mc3.num) <- paste(colnames(nosy.patient.cosmic.mc3.num), "nosy", sep = "_")
  sy.patient.cosmic.mc3.num <- PatientInCosMicMC3MutAANum(patient.sy.value, sy.cosmic.value,sy.mc3.value)
  colnames(sy.patient.cosmic.mc3.num) <- paste(colnames(sy.patient.cosmic.mc3.num), "sy", sep = "_")
  # the mut number of patient in tumor
  nosy.patient.tumor.num <- PatientInTumorMutAANum(patient.nosy.value, nosy.tumor.mc3.value)
  colnames(nosy.patient.tumor.num) <- paste(colnames(nosy.patient.tumor.num), "nosy", sep = "_")
  sy.patient.tumor.num <- PatientInTumorMutAANum(patient.sy.value, sy.tumor.mc3.value)
  colnames(sy.patient.tumor.num) <- paste(colnames(sy.patient.tumor.num), "sy", sep = "_")
  # the rate of Patient mut in cosmic
  # nosy cosmic
  nosy.cosmic.apart.value <- nosy.cosmic.value[, c("Chr", "Pos", "Ref", "Mut", "Gene_name", "AA_change", "Sample_name")]
  patient.mut.AA.cosmic.nosy.rate <- PatientMutDataBaseRate(patient.nosy.value, nosy.cosmic.apart.value)
  colnames(patient.mut.AA.cosmic.nosy.rate) <- paste("Cosmic_nosy", colnames(patient.mut.AA.cosmic.nosy.rate), sep = "_")
  # sy cosmic
  sy.cosmic.apart.value <- sy.cosmic.value[, c("Chr", "Pos", "Ref", "Mut", "Gene_name", "AA_change", "Sample_name")]
  patient.mut.AA.cosmic.sy.rate <- PatientMutDataBaseRate(patient.sy.value, sy.cosmic.apart.value)
  colnames(patient.mut.AA.cosmic.sy.rate) <- paste("Cosmic_sy", colnames(patient.mut.AA.cosmic.sy.rate), sep = "_")
  # the rate of Patient mut in MC3
  # nosy mc3
  nosy.mc3.apart.value <- nosy.mc3.value[,c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Hugo_Symbol", "AA_change", "Tumor_Sample_Barcode")]
  patient.mut.AA.mc3.nosy.rate <- PatientMutDataBaseRate(patient.nosy.value, nosy.mc3.apart.value)
  colnames(patient.mut.AA.mc3.nosy.rate) <- paste("mc3_nosy", colnames(patient.mut.AA.mc3.nosy.rate), sep = "_")
  # sy mc3
  sy.mc3.apart.value <- sy.mc3.value[, c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Hugo_Symbol", "AA_change", "Tumor_Sample_Barcode")]
  patient.mut.AA.mc3.sy.rate <- PatientMutDataBaseRate(patient.sy.value, sy.mc3.apart.value)
  colnames(patient.mut.AA.mc3.sy.rate) <- paste("mc3_sy", colnames(patient.mut.AA.mc3.sy.rate), sep = "_")
  # nosy tumor
  nosy.tumor.mc3.apart.value <- nosy.tumor.mc3.value[,c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Hugo_Symbol", "AA_change", "Tumor_Sample_Barcode")]
  patient.mut.AA.tumor.mc3.nosy.rate <- PatientMutDataBaseRate(patient.nosy.value, nosy.tumor.mc3.apart.value)
  colnames(patient.mut.AA.tumor.mc3.nosy.rate) <- paste(cancer.name, "nosy", colnames(patient.mut.AA.tumor.mc3.nosy.rate), sep = "_")
  # sy tumor
  sy.tumor.mc3.apart.value <- sy.tumor.mc3.value[, c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Hugo_Symbol", "AA_change", "Tumor_Sample_Barcode")]
  patient.mut.AA.tumor.mc3.sy.rate <- PatientMutDataBaseRate(patient.sy.value, sy.tumor.mc3.apart.value)
  colnames(patient.mut.AA.tumor.mc3.sy.rate) <- paste(cancer.name, "sy", colnames(patient.mut.AA.tumor.mc3.sy.rate), sep = "_")
  # all result
  df.patient.mut.num <- data.frame(sample_name = sample.name, patient_mut_num = patient.mut.num)
  df.patient.tcga.mc3.result <- cbind(patient.mc3tumor.exonic.mutnum.mutdepth.1KG.rsID.AsRatio.baseSub, nosy.patient.cosmic.mc3.num, sy.patient.cosmic.mc3.num, nosy.patient.tumor.num, sy.patient.tumor.num, patient.mut.AA.cosmic.nosy.rate, patient.mut.AA.cosmic.sy.rate, patient.mut.AA.mc3.nosy.rate, patient.mut.AA.mc3.sy.rate, patient.mut.AA.tumor.mc3.nosy.rate, patient.mut.AA.tumor.mc3.sy.rate)
  df.patient.mut.tcga.mc3.result <- cbind(df.patient.mut.num, df.patient.tcga.mc3.result)

  all.patient.result.list[[sample.name]] <- df.patient.mut.tcga.mc3.result

  write.table(df.patient.mut.tcga.mc3.result, paste(output.dir,"/", sample.name,".result.value.txt",sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

}

all.patient.result.value <- do.call(rbind, all.patient.result.list)

# preformance
if(answer.num >= 2 & (length(args) - 6) > 1){
  input.chr.pos.value.list <- list()
  for(i in 1:(length(args) - 6)){
    
    patient.file <- args[i]
    sample.name <- sample.names[i]
    patient.value <- fread(patient.file)
    input.chr.pos.value.list[[sample.name]] <- data.frame(chr = patient.value$Chromosome, pos = patient.value$Position)
  }
  input.chr.pos.value <- do.call(rbind, input.chr.pos.value.list)
  input.chr.pos.freq <- data.frame(table(input.chr.pos.value))
  input.chr.pos.answer <- input.chr.pos.freq[input.chr.pos.freq$Freq >= answer.num,]

  result.preformance.list <- list()
  for(j in 1:length(input.chr.pos.value.list)){

    patient.chr.pos.value <- input.chr.pos.value.list[[j]]
    patient.answer.num <- sum(paste(patient.chr.pos.value$chr, patient.chr.pos.value$pos, sep = "_") %in% paste(input.chr.pos.answer$chr, input.chr.pos.answer$pos, sep = "_"))
    patient.percision <- patient.answer.num / nrow(patient.chr.pos.value)
    patient.sensitivity <- patient.answer.num / nrow(input.chr.pos.answer)
    patient.preformance.value <- data.frame(sample_name = names(input.chr.pos.value.list)[j], patient_answer_num = patient.answer.num, answer_num = nrow(input.chr.pos.answer), percision = patient.percision, sensitivity = patient.sensitivity, F1_score = c(2 * patient.percision * patient.sensitivity) / c(patient.percision + patient.sensitivity))
    result.preformance.list[[j]] <- patient.preformance.value
  }

  result.preformance.value <- do.call(rbind, result.preformance.list)

  all.patient.result.value <- cbind(all.patient.result.value, result.preformance.value) 
  write.table(all.patient.result.value, paste(output.dir,"/", "all.result.value.txt",sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
  # venn figure
  venn.chr.pos.list <- list()
  for(j in 1:length(input.chr.pos.value.list)){

    patient.chr.pos.value <- input.chr.pos.value.list[[j]]
    patient.chr.pos.value$A <- paste(patient.chr.pos.value$chr, patient.chr.pos.value$pos, sep = "_") 
    venn.chr.pos.list[[j]] <- patient.chr.pos.value
  }
  sample.num <- length(input.chr.pos.value.list)
  filename <- sample.names
  VennFigure(venn.chr.pos.list, sample.num, filename, output.dir)

}



