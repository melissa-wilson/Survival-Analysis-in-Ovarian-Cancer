library(GEOquery)
library(reshape2)
library(survival)
library(ggplot2)
library(GGally)
library(survMisc)
library(limma)
library(BatchQC)
library(survminer)
library(hgu133a.db)
library(annotate)
library(sva)
library(coxphw)
library(magrittr)
library(dplyr)
library(Greg)


########## getting data



# get denkert data
denkert.dat = getGEO("GSE14764", GSEMatrix =TRUE, getGPL=FALSE)
denkert.dat = denkert.dat[[1]]
denkert.exp = exprs(denkert.dat)
denkert.phenotype = pData(denkert.dat)

# get bonome data
bonome.dat = getGEO("GSE26712", GSEMatrix =TRUE, getGPL=FALSE)
bonome.dat = bonome.dat[[1]]
bonome.exp = exprs(bonome.dat)
bonome.phenotype = pData(bonome.dat)
# get ride of normal cells
bonome.phenotype = bonome.phenotype[11:195,]
bonome.exp = bonome.exp[,row.names(bonome.phenotype)]

# keep only neccessary phenotype data columns
denkert.phenotype = denkert.phenotype[,c(2,26,39,40,41)]
bonome.phenotype = bonome.phenotype[,c(2,24,38,40)]

# keep only those with serous ovca
denkert.phenotype = denkert.phenotype[denkert.phenotype[,3]=="serous ovca",]
denkert.exp = denkert.exp[,row.names(denkert.phenotype)]
# all data in bonome from serous tumors

# is there batch effect? takes like 10 minutes to load
batchQC(cbind(denkert.exp,bonome.exp),batch = c(rep(1,ncol(denkert.exp)),rep(2,ncol(bonome.exp)))) # yes
all.expr = ComBat(cbind(denkert.exp,bonome.exp),batch = c(rep(1,ncol(denkert.exp)),rep(2,ncol(bonome.exp))))

# turn bonome's survival years into months and make into same format
bonome.phenotype$survival.months = 12*as.numeric(bonome.phenotype$`survival years:ch1`)
denkert.phenotype$survival.months = as.numeric(denkert.phenotype$`overall survival time:ch1`)
denkert.phenotype$survival.event = as.integer(denkert.phenotype$`overall survival event:ch1`)
bonome.phenotype$survival.event = bonome.phenotype$`status:ch1`
bonome.phenotype[bonome.phenotype$`status:ch1`=="AWD (alive with disease)",6] = 0
bonome.phenotype[bonome.phenotype$`status:ch1`=="NED (no evidence of disease)",6] = 0
bonome.phenotype[bonome.phenotype$`status:ch1`=="DOD (dead of disease)",6] = 1

# combine phenotype data
phenotype = rbind(denkert.phenotype[,c(1,2,6,7)],bonome.phenotype[,c(1,2,5,6)])
phenotype$survival.event = as.integer(phenotype$survival.event)

# need to map probe id to gene symbol
x <- hgu133aSYMBOL
mapped.probes = mappedkeys(x)
genesym.probeid = as.data.frame(x[mapped.probes])

# find probe id of genes of interest
brca1 = genesym.probeid[genesym.probeid[,2]=="BRCA1",1]
he4 = genesym.probeid[genesym.probeid[,2]=="WFDC2",1]
p53 = genesym.probeid[genesym.probeid[,2]=="TP53",1]
rassf1a = genesym.probeid[genesym.probeid[,2]=="RASSF1",1]






############## survival analysis





# adds column to data.frame
add.gene = function(dat, gene){
  
  # get name of gene
  gene.name = deparse(substitute(gene))
  
  # get gene expression level
  gene.expr = all.expr[gene[1],]
  
  # add new column
  dat = cbind(dat,gene.expr)
  
  # rename last column
  colnames(dat)[ncol(dat)] = gene.name
  
  return(dat)
}


# add genes
ov = add.gene(phenotype, brca1)
ov = add.gene(ov, he4)
ov = add.gene(ov, p53)
ov = add.gene(ov, rassf1a)

# create the survival object
survival = Surv(ov$survival.months,ov$survival.event)

# create the cox regression model
fit = coxph(survival~rassf1a+brca1+he4+p53,data = ov)
summary(fit)

# test proportional hazards assumption. tests if scaled Schoenfeld residuals ind of time
test.ph = cox.zph(fit)
test.ph # brca1 significant
ggcoxzph(test.ph) # decrease in est for brca1 as time went on

# need to split times intervals so can have time and brca1 interaction
ov.timesplit = timeSplitter(data = ov, by = 4,
               event_var = "survival.event",
               event_start_status = 0,
               time_var = "survival.months")

# should be similar to previous model
interval.fit = coxph(Surv(Start_time, Stop_time, survival.event)~rassf1a+brca1+he4+p53,data = ov.timesplit)
summary(interval.fit)

# add time interaction
time.var.fit = coxph(Surv(Start_time, Stop_time, survival.event)~rassf1a+brca1+he4+p53+brca1:Start_time,data = ov.timesplit)
summary(time.var.fit)

# drop rassf1a
time.2 = coxph(Surv(Start_time, Stop_time, survival.event)~brca1+he4+p53+brca1:Start_time,data = ov.timesplit)
summary(time.2)

# drop p53
time.3 = coxph(Surv(Start_time, Stop_time, survival.event)~brca1+he4+brca1:Start_time,data = ov.timesplit)
summary(time.3)

# drop he4
time.4 = coxph(Surv(Start_time, Stop_time, survival.event)~brca1+brca1:Start_time,data = ov.timesplit)
summary(time.4)

# time 3 maximizes concordance, final model
ggforest(time.3,data = ov.timesplit) # left decreased death, right increased death


# split into different groups to graph. use original data
ov$brca1.group = 0
ov$he4.group = 0
ov[which(ov$brca1<median(ov$brca1)),"brca1.group"] = "lower"
ov[which(ov$brca1>=median(ov$brca1)),"brca1.group"] = "upper"
ov[which(ov$he4<median(ov$he4)),"he4.group"] = "lower"
ov[which(ov$he4>=median(ov$he4)),"he4.group"] = "upper"
ov$brca1.group=as.factor(ov$brca1.group)
ov$he4.group=as.factor(ov$he4.group)

fit.km = survfit(Surv(ov$survival.months,ov$survival.event)~brca1.group + he4.group, data = ov)
ggsurvplot(fit.km, legend = "bottom", legend.title = "Group",
           legend.labs = c("Low BRCA1, Low HR4", "Low BRCA1, High HE4", "High BRCA1, Low HE4", "High BRCA1, High HE4")) + guides(colour = guide_legend(nrow = 2))
