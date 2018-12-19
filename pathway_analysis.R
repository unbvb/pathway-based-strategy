###################################################################################
##  Our strategy in discovering significant prognosis/diagnosis biomarkers takes four steps:
##  Step one: data processing,including normalization & survival/histology information processing
##  Step two: Transferring gene expression data into pathway score data
##  Step three: Picking up significant pathways
##  Step four: Finding all gene expression data and searching singnificant genes
##  Attention: gene expression data and characteristic data processing is omitted in our example 
###################################################################################

################################################################
## Relevant data information required 
## ex: Normalized gene expression data
## ex.kegg: KEGG pathway score data 
## kegg path: KEGG pathway ID and the corresponding gene symbols
## pheno: Patient characterist or clinic information
################################################################

#####
setwd()#set your destination folder
source("FAIME.R")
#######libraries#############
library(survival)

################Step one:data processing,including normalization & survival/histology information processing###########

################Step two:Transferring gene expression data into pathway score data###########
#ex:Normalized gene expression data matrix,rownames are probe ID, colnames are sample name 
#genewprobe:A vector of gene Symbol for rows of dat, the names of which is probeName
#geneset2genek:An one-to-one matrix with KEGG pathway ID, and the corresponding gene symbols

ex.kegg = runFAIME(ex, genewprobe, geneset2genek)
checkOrder(ex,pheno)
#orderDeal(ex,ex.kegg,pheno)
checkDiagData <- function(ex) #When conducting diagnosis analysis, all values in one probe can not be all same
################Step three:Picking up significant pathways###########
sig_prognosis_pathway = progPathway(os.surv,rs.surv,10e-3,10e-3)
sig_diagnosis_pathway = diagPathway(adc_pos,10e-45)

################Step four: Finding all gene expression data and searching singnificant genes################
sig_prognosis_gene = progGene(os.surv,rs.surv,ctf1=10e-7,ctf2=10e-7)
sig_diagnosis_gene = diagGeene(adc_pos,10e-45)


#smple
library("hgu133plus2.db")
library("annotate")

eset <- read.delim("~/data/gcrma_expr.txt")
OUT <- select(hgu133plus2.db, as.character(eset$X), c("SYMBOL", "ENTREZID", "GENENAME"))
ex = eset[,2:ncol(eset)]
rownames(ex) = eset[,1]
genewprobe = OUT$SYMBOL
names(genewprobe) = OUT$PROBEID
keggpath = read.delim("/home/smy/data/Genesets/kegg.txt")
geneset2genek = data.frame(keggpath$path_id, keggpath$symbol)
kegg = runFAIME(ex, genewprobe, geneset2genek)
pheno <- read.delim("F:/data/lung_cancer_13_2.pheno",header = TRUE)

checkOrder(ex,pheno)
orderDeal(ex,kegg,pheno)

################Step two################
#survival information in overall survival 
lung.death <- as.character(pheno$death)
lung.death[which(pheno$death %in% pheno$death[423])]="x: nothing"
death = unlist(strsplit(format(lung.death, justify="none"), ": "))[1:nrow(pheno)*2]
lung.survival = pheno$os_year
x <- death
d = which(death=="yes"|death=="dead"|death==1|death=="DEAD"|death=="deceased")
l = which(death=="no"|death=="alive"|death==0|death=="ALIVE")
x[d]=1
x[l]=0
x[c(-d,-l)]=NA
os.surv <- Surv(lung.survival,as.numeric(x))
  

#survival information in recurrence-free survival 
lung.recu <- as.character(pheno$recurrence)
lung.recu[which(pheno$recurrence %in% pheno$recurrence[423])]="x: nothing"
recu = unlist(strsplit(format(lung.recu, justify="none"), ": "))[1:nrow(pheno)*2]
recu.time = pheno$rfs_year
y <- recu
d = which(recu=="yes"|recu=="Y"|recu==1|recu=="relapsed"|recu=="U")
l = which(recu=="no"|recu=="N"|recu==0|recu=="not relapsed")
y[d]=1
y[l]=0
y[c(-d,-l)]=NA
y <- as.numeric(y)
rs.surv <- Surv(recu.time,y)

adc_pos <- which(pheno$histology2=="ADC")
scc_pos <- which(pheno$histology2=="SCC")
################Step three################
sig_prognosis_pathway = progPathway(os.surv,rs.surv,10e-3,10e-3)
sig_diagnosis_pathway = diagPathway(adc_pos,10e-45)
################Step four################
sig_prognosis_gene = progGene(os.surv,rs.surv,ctf1=10e-7,ctf2=10e-7)
sig_diagnosis_gene = diagGene(adc_pos,10e-45)