#!/usr/bin/env Rscript
# verion1
# author: chenfushun
# the purpose of program combine the result of many sample gvcAdvacedQC, and compare with TCGA cancer database. input: many sample gvcAdvancedQC result, mc3 database, cancer name, output dir. output: comparison figure.

# load library
suppressPackageStartupMessages(library("optparse", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("data.table", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("ggplot2", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("RColorBrewer", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))

##Specify desired options in a list

option_list <- list(
    make_option(c("-p","--patient-file"), help="many samples gvcAdvancedQC result"),
    make_option(c("-t","--tcga-file"), help="the tcga file"),
    make_option(c("-c","--cancer-name"), help="cancer name"),
    make_option(c("-s","--sample-name"), help="sample name"),
    make_option(c("-o","--output-dir"), help="the output dir, use to save the result")

)

# Function 
multi_figure_fun <- function(patient_value, mc3.value, cancer.name, SampleName, output_dir){
	
	# patient
	patient.exonic.value <- patient_value[,c("Patient_exonic_mut_num", "Patient_1KG_MAF_ALL_Num_rate", "Patient_rsID_rate", "Patient_NoSy_num", "Patient_nosy_sy_rate")]
	colnames(patient.exonic.value) <- c("Mut_num", "M_1KG_MAF_ALL_rate", "rsID_rate", "NoSy_num", "nosy_sy_rate")
	# mc3
	mc3.exonic.value <- mc3.value[mc3.value$Mutation_type != ".",]
	tumor.mc3.exonic.value <- mc3.exonic.value[mc3.exonic.value$Cancer_name %in% cancer.name,]
	tumor.mc3.exonic.depth.value <- tumor.mc3.exonic.value

	tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value <- tumor.mc3.exonic.depth.value[,.( Mut_num = length(t_depth), M_1KG_MAF_ALL_rate = c(length(`1KG_MAF_ALL`) - sum(`1KG_MAF_ALL` == "-")) / length(`1KG_MAF_ALL`) , rsID_rate = c(length(rsID) - sum(rsID == "-")) / length(rsID), nosy_sy_rate = sum(Mutation_type %in% c("nonsynonymous SNV", "stopgain")) / sum(Mutation_type == "synonymous SNV"), NoSy_num = sum(Mutation_type %in% c("nonsynonymous SNV", "stopgain"))),by = Tumor_Sample_Barcode]
	tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$nosy_sy_rate[is.infinite(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value$nosy_sy_rate)] <- NA

	# Mut_num exonic
	pdf(paste(output_dir, SampleName, ".", cancer.name,".mut_num_density.pdf", sep = ""))
	p <- ggplot(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value, aes(x=Mut_num)) +  geom_density(size = 1.2, color =  brewer.pal(8,"Set1")[2]) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15)) + xlab("Mutation Num") + xlim(0, 10000)
	p1 <- p + geom_density(data = patient.exonic.value,  size = 1.2, color =  brewer.pal(8,"Set1")[1]) 
	print(p1)
	dev.off()

	# 1KG_MAF_ALL_rate
	pdf(paste(output_dir, SampleName, ".", cancer.name,".1KG_MAF_ALL_rate.pdf", sep = ""))
	p <- ggplot(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value, aes(x=M_1KG_MAF_ALL_rate)) +  geom_density(size = 1.2, color =  brewer.pal(8,"Set1")[2]) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15)) + xlab("1KG_MAF_ALL_rate") 
	p1 <- p + geom_density(data = patient.exonic.value,  size = 1.2, color =  brewer.pal(8,"Set1")[1]) 
	print(p1)
	dev.off()

	# rsID_rate
	pdf(paste(output_dir, SampleName, ".", cancer.name,".rsID_rate.pdf", sep = ""))
	p <- ggplot(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value, aes(x=rsID_rate)) +  geom_density(size = 1.2, color =  brewer.pal(8,"Set1")[2]) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15)) + xlab("rsID_rate") 
	p1 <- p + geom_density(data = patient.exonic.value,  size = 1.2, color =  brewer.pal(8,"Set1")[1]) 
	print(p1)
	dev.off()

	# nosy_sy_rate
	pdf(paste(output_dir, SampleName, ".", cancer.name,".nosy_sy_rate.pdf", sep = ""))
	p <- ggplot(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value, aes(x=nosy_sy_rate)) +  geom_density(size = 1.2, color =  brewer.pal(8,"Set1")[2]) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15)) + xlab("nosy_sy_rate") 
	p1 <- p + geom_density(data = patient.exonic.value,  size = 1.2, color =  brewer.pal(8,"Set1")[1]) 
	print(p1)
	dev.off()

	# NoSy_num
	pdf(paste(output_dir, SampleName, ".", cancer.name,".NoSy_num.pdf", sep = ""))
	p <- ggplot(tumor.mc3.exonic.depth.mutNum.1KG.rsID.nosy.value, aes(x=NoSy_num)) +  geom_density(size = 1.2, color =  brewer.pal(8,"Set1")[2]) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black")) + theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=15)) + xlab("NoSy_num") + xlim(0, 5000)
	p1 <- p + geom_density(data = patient.exonic.value,  size = 1.2, color =  brewer.pal(8,"Set1")[1]) 
	print(p1)
	dev.off()

}

# Get command line options
arguments <- parse_args(OptionParser(usage = "%prog [options]", option_list = option_list), positional_arguments = 0)
opt <- arguments$options

PatientFile <- opt$`patient-file`
Mc3File <- opt$`tcga-file`
CancerName <- opt$`cancer-name`
SampleName <- opt$`sample-name`
OutputDir <- opt$`output-dir`

patient_value <- fread(PatientFile, header = T, sep = "\t")
mc3.value <- fread(Mc3File, header = T, sep = "\t")
multi_figure_fun(patient_value, mc3.value, CancerName, SampleName, OutputDir)
