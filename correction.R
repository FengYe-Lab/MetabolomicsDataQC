#              RAW data 
#                                     id    batch seq1   Type group X2.Hydroxypyridine Pyruvic.acid L.Lactic.acid Glycolic.acid L.Alanine
# QC_Y17I06_QC1_1         QC_Y17I06_QC1_1       1    1     QC    QC               6693        57174       1462805         55643    743047
# MM1_Y17I06_MB01_1       MM1_Y17I06_MB01_1     1    2 Sample    MM               6272        21584       2222712         38390    257925
# MM1_Y17I06_MB02_1       MM1_Y17I06_MB02_1     1    3 Sample    MM               6434        34124       1090932         36995   1169152
# MM1_Y17I06_MB03_1       MM1_Y17I06_MB03_1     1    4 Sample    MM               5894        14243       1734526         24866    510929
# MM1_Y17I06_MB04_1       MM1_Y17I06_MB04_1     1    5 Sample    MM               6571        16520       1441516         52105    838068
# MM1_Y17I06_MB05_1       MM1_Y17I06_MB05_1     1    6 Sample    MM               4589        31903       2093716         31234    331844

##   qc-RLSC and QC-RFSC come form statTarget2.0
##  
Corr = function(RAW_data,Method=c("QC-RLSC","quantileComBat","eigenMS","vsn","batchNorm","qc-RLSC","QC-RFSC")){
###########################   QC-RLSC
  if(Method=="QC-RLSC"){
QC_RLSC <- function(tab, colv, or) {
  # create table of the same sizeas initial
  tab_corr <- tab
  # For each variable (columns) in the initial table
  for (i in 1:ncol(tab)) {
    # fit loess curve to the QCs
    ll <- stats::loess(tab[which(colv == 1), i] ~ or[which(colv == 1)])
    # approximate the curve for all the samples
    aa <- stats::approx(x = or[which(colv == 1)], y = ll$fitted, xout = or)
    
    # correct the variable according to the curve for all the samples
    tab_corr[,i] <- tab[,i] / aa$y
    
    # print which variable has been corrected in order to monitor the progress
    #print(i)
  }
  return(tab_corr)
}
############################ 
rownames(RAW_data)=RAW_data[,1]
QC_RLSC_Pre=RAW_data[,-c(1:5)]
type= RAW_data[,"Type"]
type= sub("QC",1,type)
type= sub("Sample",2,type)
type=as.numeric(type)
seq1=RAW_data[,"seq1"]
QC_RLSC_data<-QC_RLSC(tab=QC_RLSC_Pre,colv=type,or=seq1)
QC_RLSC_data=cbind(RAW_data[,c("id","batch","seq1","Type","group")],QC_RLSC_data) 
return(QC_RLSC_data)
}
if(Method=="quantileComBat"){
quantileComBatFUN<-function(RAW_data){
  if(!c(require(preprocessCore),require(sva),require(pcaMethods) ,require(stringr))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("preprocessCore");biocLite("sva");biocLite("pcaMethods");biocLite("stringr");biocLite("bladderbatch")
  }
 rownames(RAW_data)= RAW_data[,1]
  metabdata=RAW_data[,-c(1:5)]
  ### un-log2
  x<-metabdata
  ### get 1/2 minimum values of each metab in non-imputed data
  minval <-metabdata
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      if (is.na(x[i, j])) {
        minval[i, j] <- summary(x[, j])[1]/2
      }
      else {
        minval[i, j] <- x[i, j]
      }
    }
  }
  ### bpca impute
  bpca0<-pca(x, nPcs=3, method="bpca")
  bpca<-completeObs(bpca0)
  
  ### replace imputed values <0 with minval
  for (i in 1:nrow(bpca)) {
    for (j in 1:ncol(bpca)) {
      if (bpca[i, j] < 0) {
        bpca[i, j] <- minval[i, j]
      }
    }
  }
  ### quantile normalize and annotate
  norm<-log2(t(normalize.quantiles(t(bpca))))
  rownames(norm)<-rownames(bpca)
  colnames(norm)<-colnames(bpca)
  phenos=unique(RAW_data$group)
  # phenos<-unique(c(A.PHENOS, QC.PHENOS))
  
  ### ComBat prep - create pheno variable from original data
  grp<-as.character(rownames(norm))
  names(grp)<-grp
  
  # phenos<-unique(c(A.PHENOS, QC.PHENOS))
  
  for (i in 1:length(phenos)){
    
    x<-phenos[i]
    
    these<-grep(x, grp, value=TRUE)
    
    grp[these]<-i
    
  }
  pheno<-as.factor(grp)
  ### specify pheno vector
  mod<-model.matrix(~as.factor(pheno))
  
  ### specify batch vector
  batch<-RAW_data[, "batch"]
  
 cat("\n","run ComBat on data","\n")
  norm2<-as.data.frame(t(ComBat(t(norm), mod=mod, batch=batch)))
  ### add batch and seq data
  norm3<-cbind(RAW_data[, c("id", "batch", "seq1","Type","group")], norm2)
  rownames(norm3)<-norm3$id
  
  return(norm3)
  
}
quantileComBat_data<-quantileComBatFUN(RAW_data)
return(quantileComBat_data)
}
if(Method=="eigenMS"){
# if only one PHENO, PHENO must contain some identifier that applies to all sample IDs
# metabolites are removed if that metabolite has no data for an entire pheno group
eigenMSFUN<-function(RAW_data){
  source("I:/Bioinfo1/Projects/QCCorrInMetabolomics/code/EigenMS.R")

  rownames(RAW_data)= RAW_data[,1]
  metabdata=RAW_data[,-c(1:5)]
  data.t<-as.data.frame(t(metabdata))
  grp<-as.character(names(data.t))
  names(grp)<-grp
  phenos<-unique(RAW_data$group)
  for (i in 1:length(phenos)){
       x<-phenos[i]
    
       these<-grep(x, grp, value=TRUE)
    
       grp[these]<-i
    
    }

  trt<-as.factor(grp)
  
  prot.info<-cbind(rownames(data.t), rownames(data.t))
  
  part1<-eig_norm1(m=data.t, treatment=trt, prot.info=prot.info)
  
  norm<-eig_norm2(rv=part1)$norm_m 
  
  norm2<-as.data.frame(t(norm))
  
  norm3<-cbind(RAW_data[, c("id", "batch", "seq1","Type","group")], norm2)
  rownames(norm3)<-norm3$id
  
  return(norm3)
  
}

eigenMS_data<-eigenMSFUN(RAW_data)
return(eigenMS_data)
}   
if(Method=="vsn"){
vsnFUN<-function(RAW_data){
  if(!require(vsn)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("vsn")
  }
  library(vsn)
  rownames(RAW_data)= RAW_data[,1]
  metabdata=RAW_data[,-c(1:5)]
  
  d<-t(metabdata)
  
  norm<-as.data.frame(t(justvsn(d)))
  
  norm2<-cbind(RAW_data[, c("id", "batch", "seq1","Type","group")], norm)
  rownames(norm2)<-norm2$id
  
  return(norm2)
  
}

vsn_data<-vsnFUN(RAW_data)
return(vsn_data)
}
if(Method=="batchNorm"){
batchNormFUN<-function(RAW_data){
rownames(RAW_data)=RAW_data$id
  metabdata=RAW_data[,-c(1:5)]
  
  metabdata.t<-as.data.frame(t(metabdata))
  QC_RAW_data=RAW_data[RAW_data$Type=="QC",]
  sample_RAW_data=RAW_data[RAW_data$Type=="Sample",]
  
    # qc<-grep(z, rownames(RAW_data), value=TRUE)
    qc<-rownames(QC_RAW_data)
    # a<-grep(paste(A.PHENOS[which(names(A.PHENOS)==z)], collapse="|"), rownames(RAW_data), value=TRUE)
    a<-rownames(sample_RAW_data)
    ### eq1
    # get total abundance, TAq, for each control sample by adding all of the log2 peak intensities for each qc sample
    TA.q<-apply(metabdata.t[, qc], 2, function(x) sum(x, na.rm=TRUE))
    
    
    # get median total abundance across all qc and real data samples, median(TA*)
    sums<-apply(metabdata.t[, c(qc, a)], 2, function(x) sum(x, na.rm=TRUE))
    
    med.TA.<-median(sums, na.rm=TRUE)
    
    
    # get scaling factor, TA.scale.q, for each control sample
    TA.scale.q<-(med.TA./TA.q)
    
    
    # multiply control-sample-specific scaling factors by observed log2 peak areas for each metab within control sample for Yp,q
    scaled<-data.frame(mapply(`*`, metabdata.t[, qc], TA.scale.q, SIMPLIFY=FALSE))
    colnames(scaled)<-colnames(metabdata.t[, qc])
    rownames(scaled)<-rownames(metabdata.t[, qc])
    
    
    # using Yp,q as outcome variables, create separate linear regression models for each metabolite and each batch.
    # treat the injection order of the control sample as a linear predictor
    scaled.t<-as.data.frame(t(scaled))
    scaled.t$id<-rownames(scaled.t)
    
    batch.pos<-RAW_data[qc, c("id", "batch", "seq1")]
    
    merged<-merge(batch.pos, scaled.t, by="id")
    rownames(merged)<-merged$id
    merged$seq1<-as.numeric(merged$seq1)
    
    ### select metab cols
    pcs<-colnames(merged[,-c(1:3)])
    
    nbatches<-length(unique(RAW_data[, "batch"]))
    
    ### create dfs for the resultant values
    results.base<-as.data.frame(matrix(NA, nrow=length(pcs), ncol=nbatches))
    rownames(results.base)<-pcs
    colnames(results.base)<-paste(1:nbatches, sep="")
    results.R<-results.base
    
    for (j in 1:length(pcs)){
      
      myMetab<-pcs[j]
      
      for (i in 1:nbatches){
        
        subset<-merged[which(merged$batch==i), ]
        
        subset$seq1<-as.numeric(subset$seq1)
        
        formula<-as.formula(paste(myMetab, "~ seq1", sep=""))
        
        if (sum(is.na(subset[, myMetab]))>=2){
          
          results.base[myMetab, i]<-NA
          results.R[myMetab, i]<-NA
          
        } else {
          
          summary<-summary(lm(formula, data=subset, na.action=na.exclude))$coefficients
          
          base.pb<-summary["(Intercept)", "Estimate"]
          R.pb<-summary["seq1", "Estimate"]
          
          results.base[myMetab, i]<-base.pb
          results.R[myMetab, i]<-R.pb
          
        }
        
      }
      
    }  
    
    ### eq2
    # get total abundance, TA.s, for each real sample by adding up all the log2 peak intensities for each sample
    TA.s<-apply(metabdata.t[, c(qc, a)], 2, function(x) sum(x, na.rm=TRUE)) 
    # median(TA*) doesn't change
    # get scaling factor, TA.scale.s, for each real sample
    TA.scale.s<-(med.TA./TA.s)
    # get the median of peak p across all qc samples, median(C.p.qc.*)
    med.C.p.qc<-apply(metabdata.t[, qc], 1, function(x) median(x, na.rm=TRUE))
    # get vector of non-normalized abundances of peak p for all real samples for the numerator
    # C<-raw.t[, c(mQC, m01, m12)] ???
    # get part 1 of the numerator, Cp,s*TA.scale.s
    part1<-data.frame(mapply(`*`, metabdata.t[, c(qc, a)], TA.scale.s, SIMPLIFY=FALSE))
    colnames(part1)<-colnames(metabdata.t[, c(qc, a)])
    rownames(part1)<-rownames(metabdata.t[, c(qc, a)])
    # get full numerator, Cp,s*TA.scale.s*median(Cp,qc,*)
    part1.t<-as.data.frame(t(part1))
    numerator<-data.frame(mapply(`*`, part1.t, med.C.p.qc, SIMPLIFY=FALSE))
    colnames(numerator)<-colnames(part1.t)
    rownames(numerator)<-rownames(part1.t)
    numerator.t<-as.data.frame(t(numerator))
    # loop through numerators and calculate normalized values
    normalized<-as.data.frame(matrix(NA, nrow=dim(numerator.t)[1], ncol=dim(numerator.t)[2]))
    rownames(normalized)<-pcs
    colnames(normalized)<-colnames(numerator.t)
    samples<-colnames(numerator.t)
    batch.pos<-RAW_data[c(qc, a), c("id", "batch", "seq1")]
    
    for (j in 1:length(samples)){
      
      mySample<-samples[j]
      
      for (i in 1:length(pcs)){
        
        myMetab<-pcs[i]
        
        numerator.p.s<-numerator.t[myMetab, mySample]
        
        batch<-batch.pos[mySample, "batch"]
        
        seq1<-as.numeric(batch.pos[mySample, "seq1"])
        
        base<-results.base[myMetab, batch]
        
        R<-results.R[myMetab, batch]
        
        denominator<-(R*seq1)+base
        
        final<-numerator.p.s/denominator
        
        normalized[myMetab, mySample]<-final
        
      }
      
    }  
  norm<-as.data.frame(t(normalized))
  # norm1<-do.call("rbind", norm)
  norm2<-norm[rownames(RAW_data), ]
  
  norm3<-cbind(RAW_data[, c("id", "batch", "seq1","Type","group")], norm2)
  rownames(norm3)<-norm3$id
  
  return(norm3)
  
}
batchNormalizer_data<-batchNormFUN(RAW_data)
return(batchNormalizer_data)
}
if(Method=="qc-RLSC"){
  qcRLSC<-function(RAW_data){
 if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
if(!require(statTarget))
 BiocManager::install("statTarget", version = "devel")

    library(statTarget)
    Dr_path=getwd()
    name=colnames(RAW_data)[-c(1:5)]
    statTarget_RAW_data=t(RAW_data[,-c(2:5)])
    row.names(statTarget_RAW_data)[1]=""
    colnames(statTarget_RAW_data)=statTarget_RAW_data[1,]
    statTarget_RAW_data=data.frame(statTarget_RAW_data[-1,])
    statTarget_RAW_data=as.data.frame(cbind(name,statTarget_RAW_data))
    write.csv(statTarget_RAW_data,"statTarget_RAW_data.csv",row.names = F)
    
    statTarget_RAW_Pheno_data=RAW_data[,c("id","batch","group","seq1")]
    statTarget_RAW_Pheno_data$group=sub("QC","NA",statTarget_RAW_Pheno_data$group)
    names(statTarget_RAW_Pheno_data)=c("sample","batch","class","order")
    write.csv(statTarget_RAW_Pheno_data,"statTarget_RAW_Pheno_data.csv",row.names = F)
    
    
    samPeno <- paste(Dr_path,"statTarget_RAW_Pheno_data.csv", sep="/")
    samFile <- paste(Dr_path,"statTarget_RAW_data.csv", sep="/")
dir.create("qcRLSC(statTarget)")
cat("\n")
cat("The smoothing parameter for QC-RLSC which controls the bias-variance tradeoff in QC-RLSC method 
if the QCspan is set at '0', the generalised cross-validation will be performed to avoid overfitting 
    the observed data.","\n")

QCspan = readline("Please enter the parameters of QCspan:")
cat("\n","Lets you specify local constant regression (i.e., the Nadaraya-Watson estimator, degree=0),
  local linear regression (degree=1),or local polynomial fits (degree=2, the default) for QC-RLSC","\n")
degree = readline("Please enter the parameters of degree:")
cat("\n","The parameter for imputation method i.e., nearest neighbor averaging,'KNN'; minimum values, 'min'; Half of minimum values, 'minHalf'; median values, 'median'.","\n")
cat("KNN,","min,","minHalf,","median")
imputeM = readline("Please enter the parameters of imputeM:")
shiftCor(samPeno,samFile, Frule = 0.8, MLmethod = "QCRLSC", ntree = 500,
         QCspan , degree , imputeM, plot = FALSE)

qcRLSC_data=read.csv("qcRLSC(statTarget)/statTarget/shiftCor/After_shiftCor/shift_all_cor.csv")
rownames(qcRLSC_data)=qcRLSC_data$sample
qcRLSC_data=cbind(RAW_data[,c("id","batch","seq1","Type","group")],qcRLSC_data[,-c(1:3)])
  }
  qcRLSC_data=qcRLSC(RAW_data)
  return(qcRLSC_data)
}
  
if(Method=="QC-RFSC"){
  QCRFSC=function(RAW_data){
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    if(!require(statTarget))
      BiocManager::install("statTarget", version = "devel")
    
    library(statTarget)
    Dr_path=getwd()
    name=colnames(RAW_data)[-c(1:5)]
    statTarget_RAW_data=t(RAW_data[,-c(2:5)])
    row.names(statTarget_RAW_data)[1]=""
    colnames(statTarget_RAW_data)=statTarget_RAW_data[1,]
    statTarget_RAW_data=data.frame(statTarget_RAW_data[-1,])
    statTarget_RAW_data=as.data.frame(cbind(name,statTarget_RAW_data))
    write.csv(statTarget_RAW_data,"statTarget_RAW_data.csv",row.names = F)
    
    statTarget_RAW_Pheno_data=RAW_data[,c("id","batch","group","seq1")]
    statTarget_RAW_Pheno_data$group=sub("QC","NA",statTarget_RAW_Pheno_data$group)
    names(statTarget_RAW_Pheno_data)=c("sample","batch","class","order")
    write.csv(statTarget_RAW_Pheno_data,"statTarget_RAW_Pheno_data.csv",row.names = F)
    
    
    samPeno <- paste(Dr_path,"statTarget_RAW_Pheno_data.csv", sep="/")
    samFile <- paste(Dr_path,"statTarget_RAW_data.csv", sep="/")
    dir.create("QCRFSC(statTarget)")
    cat("\n","The parameter for imputation method i.e., nearest neighbor averaging,'KNN'; minimum values, 'min'; Half of minimum values, 'minHalf'; median values, 'median'.","\n")
    cat("Please choose the following method:","\n","KNN,","min,","minHalf,","median","\n")
    imputeM = readline("Please enter the parameters of imputeM:")
    shiftCor(samPeno,samFile, Frule = 0.8, MLmethod = "QCRFSC", ntree = 500,
             QCspan = 0, degree = 2, imputeM, plot = FALSE)
    QCRFSC_data=read.csv("QCRFSC(statTarget)/statTarget/shiftCor/After_shiftCor/shift_all_cor.csv")
    rownames(QCRFSC_data)=QCRFSC_data$sample
    QCRFSC_data=cbind(RAW_data[,c("id","batch","seq1","Type","group")],QCRFSC_data[,-c(1:3)])
  }
  QCRFSC_data=QCRFSC(RAW_data)
    return(QCRFSC_data)
}
}


# QCRLSC_data=Corr(RAW_data = RAW_data,Method = "QC-RLSC")
# quantileComBat_data=Corr(RAW_data = RAW_data,Method = "quantileComBat")
# eigenMS_data=Corr(RAW_data = RAW_data,Method = "eigenMS")
# vsn_data=Corr(RAW_data = RAW_data,Method = "vsn")
# batchNorm_data=Corr(RAW_data = RAW_data,Method = "batchNorm")
# qcRLSC_data=Corr(RAW_data = RAW_data,Method ="qc-RLSC" )
# QCRFSC_data=Corr(RAW_data = RAW_data,Method = "QC-RFSC")

