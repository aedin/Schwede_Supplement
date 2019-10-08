#######################
## Functions to predict Ov Molecular Subtypes
#########################
# To run source("./src/ProcSubtypes.R")


##########################################################
## Functions For processing Birrer and AOCS Annotations
########################################################

## Stroma fRMA normalized Data
procBirrerData<-function(birrer.frma) {
    rownames(birrer.frma) <- birrer.frma[,1]
    birrer.frma <- birrer.frma[,-1]
    colnames(birrer.frma) <- sub(".CEL.gz","",colnames(birrer.frma))
    
    ## Edits to colnames Birrer.frma
    vec <- sub("HGU133Plus2.0.CEL","",colnames(birrer.frma))
    vec <- sub("HGU133Plus2.CEL","",vec)
    vec <- sub("X090506birrer","",vec)
    vec <- sub("IVTf","",vec)
    vec <- sub("IVTF","",vec)
    vec <- sub("X091406birrer","",vec)
    vec <- sub("X091905birrer","",vec)
    vec <- sub("X0900606birrer","",vec)
    vec <- sub("X091405birrer","",vec)
    vec <- sub("X092905birrer","",vec)
    vec <- sub("SM","SM-",vec)
    vec <- sub("GSM-","GSM",vec)
    colnames(birrer.frma) <- vec
    return(birrer.frma)
}



procReadSig<-function(sigfile=file.path("data","Schwede_SupplTable 1 -TumorStromalGenes.csv")) {
  geneSig<-scan(sigfile, what="tx", sep="\n")
  geneSig<-sapply(geneSig, strsplit, ",")
  names(geneSig)<-unlist(as.character(lapply(geneSig, function(x) x[1])))
  geneSig<-lapply(geneSig, function(x) x[-1])
  geneSig<-lapply(geneSig, function(x) x[!x==""])
  return(geneSig)
}




procAnnot<-function(annot) {
  colnames(annot) <-apply(annot[1,], 2, as.character)
  annot<-annot[-1,]
  annot[1:2,]
  ## Subset to the set of stroma samples which are in GSE1850 or are normal stroma
  annot$In_GSE1850<-as.character(annot[,2])%in% pData(GSE18520_eset)[,1]
  annot<-annot[annot$"Sample Type"=="Normal Stroma" | annot$In_GSE1850,]
  
  annot<-annot[!rownames(annot)%in%"4",]
  #  > comparelists(annot[,2],  mm$ID)$Set.Diff
  #  [1] "386"  "412"  "1231" "1660" "N13"  "N28"  "N37"  "N41"  "N43"  "N47" 
  #  [11] "N5"   "N6"   "317"  "351"  "744"  "358"  "718"  "934"  "N14"  "N18"
  
  ## Merge matching column Names
  colnames(annot)<-tolower(sub(" ", "_" ,colnames(annot)))
  colnames(annot)[2] <-"alt_sample_name"
  # colnames(annot)[3] <-"sample_type"
  colnames(annot)<-sub(":", "", colnames(annot))
  colnames(annot)[4]<-"days_to_death"
  
  annot$array_id<-as.character(annot$array_id)
  ## Correct mis-matched IDs
  annot$array_id[grep("786", annot$array_id)]<-"SM-N786"
  annot$array_id[grep("692", annot$array_id)]<-"SM-N692"
  annot$array_id[grep("694", annot$array_id)] <-"SM-N694"
  annot$array_id[grep("656", annot$array_id)]<-"SM-N656"
  annot$array_id[grep("970", annot$array_id)]<-"SM-N970"
  
  ## Subset annotation for GSE18520 and cbind both sets
  
  tAnnot<-pData(GSE18520_eset)[pData(GSE18520_eset)[,1]%in% annot[,2],]  
  
  # Convert surival months to days
  annot[,4] <-as.numeric(annot[,4])/12*365
  newCols<-setdiff(colnames(annot), colnames(tAnnot))
  tAnnot<-cbind(tAnnot, data.frame(matrix(NA, ncol=length(newCols), nrow=nrow(tAnnot))))
  colnames(tAnnot)[c(ncol(tAnnot)-length(newCols)+1):ncol(tAnnot)]<-newCols
  newCols<-setdiff(colnames(tAnnot),colnames(annot))
  annot<-cbind(annot, data.frame(matrix(NA, ncol=length(newCols), nrow=nrow(annot))))
  colnames(annot)[c(ncol(annot)-length(newCols)+1):ncol(annot)]<-newCols
  rownames(annot)<-annot$array_id
  annotNew<-annot[,colnames(tAnnot)]
  annotNew<-rbind(tAnnot, annot)
  return(annotNew)
}




######################################
#Molecular Subtyvincepes classification 
#######################################

robustScaling<-function(data, std = c("none", "scale", "robust"), verbose=FALSE) {
  std <- match.arg(std)
  switch(std
      , scale = {
    data <- scale(data, center = TRUE, scale = TRUE)
    if (verbose) {
      message("standardization of the gene expressions")
    }
  }  , robust = {
    data <- apply(data, 2, function(x) {
      return((rescale(x, q = rescale.q, na.rm = TRUE) - 
                0.5) * 2)
    })
    if (verbose) {
      message("robust standardization of the gene expressions")
    }
  }  , none = {
    if (verbose) {
      message("no standardization of the gene expressions")
    }
  })
 return(data)
}



build.tcga<-function(eSet=TCGA_eset,  TCGAgenes=as.character(tcgaGenes[,1])) {
  MS<- levels(factor(TCGA_eset$SUBTYPE))
  tcga_scm<-sapply(MS, function(x) rowMeans(exprs(eSet[TCGAgenes,grep(x, TCGA_eset$SUBTYPE)])))
  #save(TCGA_eset, tcga_scm, TCGAgenes, file="TCGA_scm_Classifier_Verhaark_100genes.Rda") 
  return(tcga_scm)
}


build.tcga.mod2<-function(eSet=TCGA_eset, TCGAgenes=tcgaGenes){
  modelMat<-model.matrix(~-1+tcgaGenes$CLASS)
  colnames(modelMat) = sort(m)
  modelMat[1:2,]
  pData(TCGA_eset)$SUBTYPE[1:2]

  sapply(m, function(i) {
    gg<-as.character(tcgaGenes$gene)[modelMat[,i]==1]
    TCGA_eset[gg, grep(i, TCGA_eset$SUBTYPE)]
    tt<-lmFit(TCGA_eset[tcgaGenes[,1],!is.na(TCGA_eset$SUBTYPE)], model.matrix(~-1+TCGA_eset$SUBTYPE)) 
    })
}


build.tcga.mod3<-function(eSet=TCGA_eset, TCGAgenes=tcgaGenes){
  m= levels(factor(TCGA_eset$SUBTYPE))
  mS<-sub("IMM", "IMR", toupper(substr(m,1,3)))
  ## do a simple weighted sum over genes for each subtype
  pred<-t(apply(exprs(TCGA_eset[tcgaGenes[,1],] ), 2, function(x)(t(modelMat)%*%x)/apply(modelMat,2,sum)))
  pred<-data.frame(pred)
  colnames(pred) = colnames(modelMat)
  pred$Predict  = apply(pred,1,function(x)colnames(modelMat)[which.max(x)])
  return(pred)
  }


predict.tcga<-function(eSet=TCGA_eset,  tcga.model=tcga_scm, corThres=0.7){
  
 # MSgenes are the molecular subtype gene signatures
 # load(datadir,"TCGA_scm_Classifier_Verhaark_100genes.Rda")

  if (inherits(eSet,"ExpressionSet")) eSet<-exprs(eSet)
  eSet<-data.frame(eSet)
  indgenes<-intersect(rownames(tcga.model), rownames(eSet))
  eSet<-eSet[indgenes,]
  tcga.model<-tcga.model[indgenes,]
  #print(dim(eSet))
  #print(class(eSet))
 
    testRes<-data.frame(t(sapply(eSet, function(x) cor(x,tcga.model, method="spearman"))))
  #print(dim(testRes))
  colnames(testRes) = colnames(tcga.model)
  
  
  nprob <- t(apply(X = testRes[,1:4], MARGIN = 1, FUN = function(x) {
    return(abs(x)/sum(abs(x), na.rm = TRUE))
  }))
  
  colnames(nprob) = colnames(testRes)
  
  
  testRes$Predict<-apply(testRes[,1:4], 1, function(x) ifelse (max(x)>corThres,colnames(testRes)[which.max(x)], "Unclassified"))
  
 
  ## Add posterior probability (approach from genefu:::intrinsic.cluster.predict)
  
  testRes$PredictPosteriorProb<- sapply(seq_along(testRes$Predict), function(x) if(!testRes$Predict[x]=="Unclassified") return(nprob[x, testRes$Predict[x]]) else return(NA))
  
  
  return(testRes)
}
 






build.aocs<-function(aocs_eSet,MSgenes ) {
  aocs_scm<-sapply(as.character(paste("C", c(1,2,4,5), sep="")), function(x) rowMeans(exprs(aocs_eSet[MSgenes,grep(x,aocs_eSet$MS)])))
 # save(aocs_eSet, aocs_scm, MSgenes, file="AOCS_Classifier_TothillGenes.Rda") 
  return(aocs_scm)
}

predict.aocs<-function (eSet=aocs_eSet, aocs.model=aocs.scm, corThres=0.7) {
  
  # Load "allgenes" ,"aocs.scm" which are the list of genes and model respectively
  if (is.null(aocs.model)) load(file.path(datadir,"AOCS_Classifier_TothillGenes.Rda"))
  
  # Extract the gene subset from eSet
  if (inherits(eSet,"ExpressionSet")) eSet<-exprs(eSet)
  
  indgenes<-intersect(rownames(aocs.model), rownames(eSet))
#  print(length(indgenes))
  if (length(indgenes)< nrow(aocs.model)) print(paste(length(indgenes), "identifiers mapped to the", nrow(aocs.model),  "AOCS Affymeytrix probeset identifiers"))
  eSet<-data.frame(eSet[indgenes,])
 # print(dim(eSet))
  
  aocs.model<-aocs.model[indgenes,]
  #print(dim(aocs.model))
  
  # Classify each sample
  testRes<-data.frame(t(sapply(eSet, function(x) cor(x, aocs.model, method="spearman"))))
  colnames(testRes) = colnames(aocs.scm)
  
  nprob <- t(apply(testRes[,1:4], 1, function(x) abs(x)/sum(abs(x), na.rm = TRUE)))

  testRes$Predict<-apply(testRes[,1:4], 1, function(x) ifelse (max(x)>corThres,colnames(testRes)[which.max(x)], "Unclassified"))
  
  ## Add posterior probability (approach from genefu:::intrinsic.cluster.predict)
 
  testRes$PredictPosteriorProb<- sapply(seq_along(testRes$Predict), function(x) if(!testRes$Predict[x]=="Unclassified") return(nprob[x, testRes$Predict[x]]) else return(NA))
  
  return(testRes)
}


plot.predict<-function(x=testRes, corThres=0.7, colLab=NULL, reorderByClass=TRUE, xlab="Case", ...) { 
   
    aocs<-c("gray80", "grey20", "navy", "cyan")
    names(aocs)[1:4] =paste0("C",c(1,2,4,5))
  
      
    tcga<- c("cyan", "grey20", "gray80", "navy")
    names(tcga)[1:4] <-c( "Differentiated", "Immunoreactive", "Mesenchymal", "Proliferative") 
    
    if (colnames(x)[1]%in% names(aocs)) ssCols<-aocs
    if (colnames(x)[1]%in% names(tcga)) ssCols<-tcga
    
    #print(ssCols)
    ## Correlation Plot
    
    #ylim=round(range(x[,1:4]), 1)+c(-0.1, 0.1)
    ylim=c(0.3,1)
    if(reorderByClass) {
      ind<-order(x$Predict)
      x<-x[ind,]
      if (!is.null(colLab)) colLab = colLab[ind]
    }
    
    plot(1:nrow(x), x[,1], ylim=ylim, bg=ssCols[1],  col="black",ylab="predict", pch=21,xlab=xlab, xaxt="n")
    for (i in 2:4) points(1:nrow(x), x[,i], bg=ssCols[i], col="black", pch=21)
    if (!is.null(colLab)) axis(1,1:nrow(x), labels=colLab,  tck = 1, col = "grey", lty = "dotted", las=2)
    abline(h=corThres, col="red", lty=2)
    legend("bottomleft", legend=names(ssCols), fill=ssCols,horiz =TRUE, cex=0.6)
    #legend("bottomleft", legend=levels(factor(x$Predict)), fill=ssCols,horiz =TRUE, cex=0.8)
    
    
   # require(ggplot2)
  #  x2<-melt(birrerPredict[,1:4])
   # p<- ggplot(x2, aes(x=rep(seq_along(x[,1]),4),y=value, colour=variable))+geom_point(size=3)
  #  p+scale_colour_manual(values=ssCols)+ theme_bw()  
  #  + facet_wrap(~x$Predict)
    
    if(reorderByClass) { 
      abline(v=cumsum(table(x$Predict[order(x$Predict)]))+0.5, col=ssCols)
      return(ind)
    }
    }


distCor <- function(x, use = "pairwise.complete.obs") {
         co.x <- cor(x, use = use)
         dist.co.x <- 1 - co.x
         return(as.dist(dist.co.x))
       }


procAOCS<-function() {
  ## Other clinical annotation NO extra Info
  require(GEOquery)
  #aocsAnnot<-read.csv("http://www.ncbi.nlm.nih.gov/geosuppl/?acc=GSE9891&file=GSE9891%5Fclinical%5Fanns%2Ecsv")
  GSE9891<-getGEO("GSE9891")
  # this contain 2 series, the main AOCS series and 30 microdissected
  #Reanalyzed by: GSE40785:238  
  aocs_eSet<-GSE9891[[1]][,GSE9891[[1]]$characteristics_ch1.2=="Subtype : Ser/PapSer" & GSE9891[[1]]$characteristics_ch1.1=="Type : Malignant" ]
  summary(pData(aocs_eSet))
  
  
  ## The GSE does not come with MS calls.. get that from suppl data from Tothill et al., paper.
  aocs<-read.csv(file.path(datadir,"AOCSclinicaldata.csv"), as.is=TRUE, row.names=1)
  rownames(aocs) <-toupper(sub(".cel$", "", rownames(aocs)))
  
  aocs<-aocs[sampleNames(aocs_eSet),]
  aocs<-cbind(aocs,pData(aocs_eSet))
  
  # curate aocs check
  table(sub("Consolidated.Grade : ", "",aocs[, 35]) == aocs$Grade)
  aocs<-aocs[,!duplicated(colnames(aocs))]
  pData(aocs_eSet) = aocs
  ## Note, we are excluding C3, C6 and NC....
  aocs_eSet$MS<-paste("C", aocs_eSet$k.means.Group, sep="")
  aocs_eSet$MS<-sub("CNC", "", aocs_eSet$MS)
  aocs_eSet$MS<-sub("C3", "", aocs_eSet$MS)
  aocs_eSet$MS<-sub("C6", "", aocs_eSet$MS)
  aocs_eSet$MS[aocs_eSet$MS==""]<-NA
  save(aocs_eSet, file=file.path(datadir,"AOCS_eSet.RData"))
}


predStats<-function(rTrue, rPred){
  ## Expect 2 binary vectors with 1=TRUE, 0=FALSE
  ## A true positive (TP)
  ## B false positve (FP)
  ## C false negative (FN)
  ## D true negative (TN)
  ## Precision is the same as PPV
  ## See http://en.wikipedia.org/wiki/Specificity_(tests)
  
  
  rTab<- table(rTrue, rPred)
  #  print(rTab)
  A<-rTab["1","1"]
  B<-rTab["0","1"]
  C<-rTab["1","0"]
  D<-rTab["0","0"]
  #  print(c(A,B,C,D))
  accuracy = sum(A,D)/sum(A,B,C,D)
  precision = A/sum(A,B)    ## Same at the PPV
  sensitivity = A/sum(A,C)
  specificity  = D/sum(B,D)
  PPV = A/sum(A,B)
  NPV =  D/sum(C,D)
  statsRes<-c(accuracy=accuracy,sensitivity=sensitivity,specificity =specificity ,PPV=PPV,NPV=NPV)
  statsRes<-round(statsRes,3)
  
  return(statsRes)
}

fisherT<-function(vec1, vec2, universe=NULL, FisherAlternative="greater") {
       vec1<-as.character(vec1)
       vec2<-as.character(vec2)
       #print(comparelists(vec1, vec2))
       if (is.null(universe)) universe <-union(vec1,vec2)
       v1=factor(universe%in%vec1 , levels=c(TRUE, FALSE))
       v2=factor(universe%in%vec2, levels=c(TRUE, FALSE))
       contingency_table<-matrix(as.vector(table(universe%in%vec1,universe%in%vec2))[4:1],2,2)
       print("contingency_table")
       print(contingency_table)
       #res<-fisher.test(table(v1,v2), alternative=FisherAlternative)
       res<-fisher.test(contingency_table,alternative=FisherAlternative )
       return(res)
}



###################
## Barcode Functions
####################

myauc <- function(x, cov1) {
  require(ROC)
  # x is vector (gene) of values
  # cov1 is the phenoData (binary vector, 0,1)
  # dxrule.sca is a simple function object ifelse(x > thresh, 1, 0)
  # Return the AUC for a gene
  AUC(rocdemo.sca(truth = cov1, data = x, rule = dxrule.sca))
}





runResampAUC<-function(eSet=birrer_eSet, cov=as.numeric(birrer_eSet$sample_type%in%"Tumor"),nResamp=100){
  set.seed(123)
  out<- sapply(1:nResamp, function(i) {esApply(eSet[, sample(1:ncol(eSet), size=nrow(eSet), replace=TRUE)], 1, function(i) myauc(i, cov1= cov))})		
  colnames(out) = paste("N",1:nResamp, sep="")
  rownames(out) = featureNames(eSet)	
  return(out)
}

heatplot<-function (dataset, dend = c("both", "row", "column", "none"), 
                    cols.default = TRUE, lowcol = "green", highcol = "red", scale = "none", 
                    classvec = NULL, classvecCol = NULL, classvec2 = NULL, classvec2Col = NULL,
                    distfun = NULL, returnSampleTree = FALSE, method = "ave", dualScale = TRUE, 
                    zlim = c(-3, 3), scaleKey = TRUE, ...) 
{
  library(gplots)
  data <- array2ade4(dataset)
  data <- as.matrix(data)
  if (dualScale) {
    print(paste("Data (original) range: ", round(range(data), 
                                                 2)[1], round(range(data), 2)[2]), sep = "")
    data <- t(scale(t(data)))
    print(paste("Data (scale) range: ", round(range(data), 
                                              2)[1], round(range(data), 2)[2]), sep = "")
    data <- pmin(pmax(data, zlim[1]), zlim[2])
    print(paste("Data scaled to range: ", round(range(data), 
                                                2)[1], round(range(data), 2)[2]), sep = "")
  }
  distEisen <- function(x, use = "pairwise.complete.obs") {
    co.x <- cor(x, use = use)
    dist.co.x <- 1 - co.x
    return(as.dist(dist.co.x))
  }
  cols <- function(low = lowcol, high = highcol, ncolors = 123) {
    low <- col2rgb(low)/255
    if (is.character(high)) 
      high <- col2rgb(high)/255
    col <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
                                                        high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
    return(col)
  }
  cols.gentleman <- function() {
    library(RColorBrewer)
    hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
    return(rev(hmcol))
  }
  if (cols.default) 
    plotcols = cols.gentleman()
  else plotcols = cols()
  if (is.null(distfun)) 
    distf = distEisen
  else distf = function(x) dist(t(x), method = distfun)
  Colv <- FALSE
  Rowv <- FALSE
  dend <- match.arg(dend)
  dend = tolower(dend)
  if (dend %in% c("row", "r")) {
    print(dend)
    Rowv = as.dendrogram(hclust(distf(t(data)), method = method))
  }
  if (dend %in% c("column", "col", "c")) {
    print(dend)
    Colv = as.dendrogram(hclust(distf(data), method = method))
    dend = "column"
  }
  if (dend %in% c("both", "TRUE")) {
    Colv = as.dendrogram(hclust(distf(data), method = method))
    Rowv = as.dendrogram(hclust(distf(t(data)), method = method))
    dend = "both"
  }
  RSideColors = CSideColors = NULL
  if (any(!is.null(classvec), !is.null(classvec2))) {
    proc.classvec <- function(classvec,classvecCol) {
      classvec = as.factor(classvec)
      if (is.null(classvecCol)) 
        classvecCol = getcol(length(levels(classvec)))
      SideCols = factor(classvec, labels = classvecCol)
      print(cbind(Class = levels(classvec), Color = levels(SideCols)))
      SideCols = as.character(SideCols)
      nSC = length(SideCols)
      return(list(nSC, SideCols))
    }
    if (!is.null(classvec)) {
      out = proc.classvec(classvec,classvecCol)
      nSC = out[[1]]
      SideCols = out[[2]]
      if (!nSC %in% dim(data)) 
        print("Error: classvec length not equal to nrow or ncol in data")
      if (nSC == nrow(data)) 
        RSideColors = SideCols
      if (nSC == ncol(data)) 
        CSideColors = SideCols
    }
    if (!is.null(classvec2)) {
      out = proc.classvec(classvec2,classvec2Col)
      nSC = out[[1]]
      SideCols = out[[2]]
      if (!nSC %in% dim(data)) 
        print("Error: classvec2 length not equal to nrow or ncol in data")
      if (nSC == nrow(data)) 
        RSideColors = SideCols
      if (nSC == ncol(data)) 
        CSideColors = SideCols
    }
  }
  if (all(is.null(RSideColors), is.null(CSideColors))) 
    heatmap.2(data, Colv = Colv, Rowv = Rowv, col = plotcols, 
              scale = scale, trace = "none", density.info = "none", 
              zlim = zlim, dendrogram = dend, ...)
  if (all(!is.null(RSideColors), is.null(CSideColors))) 
    heatmap.2(data, Colv = Colv, Rowv = Rowv, col = plotcols, 
              scale = scale, trace = "none", density.info = "none", 
              RowSideColors = RSideColors, dendrogram = dend, zlim = zlim, 
              ...)
  if (all(is.null(RSideColors), !is.null(CSideColors))) 
    heatmap.2(data, Colv = Colv, Rowv = Rowv, col = plotcols, 
              scale = scale, trace = "none", density.info = "none", 
              ColSideColors = CSideColors, dendrogram = dend, zlim = zlim, 
              ...)
  if (all(!is.null(RSideColors), !is.null(CSideColors))) 
    heatmap.2(data, Colv = Colv, Rowv = Rowv, col = plotcols, 
              scale = scale, trace = "none", density.info = "none", 
              RowSideColors = RSideColors, ColSideColors = CSideColors, 
              zlim = zlim, dendrogram = dend, ...)
  if (returnSampleTree) 
    return(Colv)
}
