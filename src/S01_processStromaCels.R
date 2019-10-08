##########################
## Testing Stroma Tumor Birrer Dataset
## Tumor Data GSE18520
## Stroma cel files in matt's directory
## Aedin March 2013
##########################

require(curatedOvarianData)
source("/cbhomes/aedin/scripts/affy/make_eSet.R")
source("/home/aedin/Dropbox/Ovarian/TumorStromaSubtypes_AC/src/ProcSubtypes.R")

dataDir<- "/cbhomes/schwede/Birrer Stroma Cel Files"
## get stroma cell files (in Matt's Dir)
stroma <- read.csv(file.path(dataDir,"frma.csv"))

## Read original annotation that Mike gave us with cel files
annot<-read.csv(file.path(dataDir,"Stroma Sample Database 9-19-06.csv"), header=TRUE, as.is=TRUE, skip=1)
annot2<- read.csv(file.path(dataDir,"Birrer clinical data tumor.csv"))
rm(annot2)

## Output Dir
setwd("~/Dropbox/Ovarian/Matt")


## Mok data set GSE18520
require(curatedOvarianData)
#data(package="curatedOvarianData")
data("GSE18520_eset")


## Or load already read in: 
load("TumStroma.RData")
birrer.frma<-procBirrerData(stroma)
geneSig<-procReadSig()
annotNew<-procAnnot(annot)


## Matching of stroma and tumor samples
mm<-read.csv("Birrer_Mok_Matched_Epi_Stroma.csv", as.is=TRUE)
if (all(sub("_Epi", "", mm[,1])==sub("_Stroma", "", mm[,2]))) {
  mm$ID= sub("_Epi", "", mm[,1])
}

annotNew$in_matchedList<-annotNew$alt_sample_name%in%mm[mm$ID%in%annotNew$alt_sample_name,"ID"]

write.csv(annotNew, file="UpdatedAnnotation.csv")
## Check, manually do... remove annoted samples for which there are no cell files

problems<-annotNew[!rownames(annotNew)%in%colnames(birrer.frma),]
write.csv(problems, file="problemsMatchingCelFiles2Annot.csv")

## For the moment will exclude and just create new eSet
annotNew<-annotNew[rownames(annotNew)%in%colnames(birrer.frma),]

## After I exclude the ones for which there are no cell files, I see that there are no matched stroma:tumor pairs for the follow
## "312" "317" "351" "715" "718" "744" These are the one listed in problems of which two are listed in Birrer_Mok_Matched_Epi_Stroma as being there, however all have Failed in status.  Therefore will exclude these

names(which(table(annotNew$alt_sample_name)==1))
#  "312" "317" "351" "715" "718" "744" "N13" "N14" "N18" "N28" "N37" "N41" "N43" "N47" "N5"  "N6"
## The ones with N are the normal tissue
## for the following have tumor but no stroma
annotNew[annotNew$alt_sample_name%in%names(which(table(annotNew$alt_sample_name)==1))&!annotNew$sample_type=="Normal Stroma", ]

# Therefore exclude these unmatched pairs, leaving 86 samples which include 10 normal and 76 matched pairs (n=38)
annotNew<-annotNew[!rownames(annotNew)%in%rownames(annotNew[annotNew$alt_sample_name%in%names(which(table(annotNew$alt_sample_name)==1))&!annotNew$sample_type=="Normal Stroma", ]),]

## Subset the birrer data and make eSet
rownames(annotNew)%in%colnames(birrer.frma)
annotNew$ShortLabel=paste(as.character(factor(annotNew$sample_type, labels=c("", "S", "T"))), annotNew$alt_sample_name, sep="")

fData(birrer_eSet)<-makefData(birrer_eSet)
save(birrer_eSet, file="birrer_eSet.RData")

## Update Annotation
## From Rmd File

load(file=file.path(datadir,"birrer_eSet_complete.RData"))
birrer_eSet$sample_type<- factor(birrer_eSet$sample_type, levels=c("tumor", "Tumor Stroma", "Normal Stroma"), labels=c("Tumor", "Tumor_Stroma", "Normal_Stroma"))
birrer_eSet$ShortLabel<-factor(birrer_eSet$ShortLabel, levels=c("T", "S", "N"))
#birrer_eSet$ShortLabel <-as.character(factor(birrer_eSet$sample_type, labels=c("N","T", "S")))
#table(birrer_eSet$sample_type,birrer_eSet$ShortLabel)

# batch effect by date, can't remove as confounded
#heatplot(cor((exprs(birrer_eSet))), classvec=substr(birrer_eSet$batch,1,7))
#birrer_eSet$batch_yrsMonth= factor((substr(birrer_eSet$batch,1,7)))
#yrs= levels(birrer_eSet$batch_yrsMonth)
#legend("topleft", yrs,fill=getcol(length(yrs)))


makefData<-function(eSet, chip=hgu133plus2SYMBOL){
  fData<-toTable(chip)
  fData= rbind(fData, cbind(probe_id=featureNames(eSet)[!(featureNames(eSet)%in%fData[,1])], symbol=NA))
  rownames(fData) =fData[,1]
  fData= fData[featureNames(eSet),]
  return(fData)
}

fData(birrer_eSet)<-makefData(birrer_eSet)
save(birrer_eSet,file=file.path(datadir,"birrer_eSet_complete.RData"))

##


pdf(file="heatplot_SchwedeSig_BirrerData.pdf")
heatplot(birrer_eSet[c(geneSig[[1]],geneSig[[2]]),], classvec=birrer_eSet$sample_type, classvec2=rep(c("green", "yellow"), times=sapply(geneSig, length)[1:2]), labCol=birrer_eSet$ShortLabel, labRow="")

legend("topright", legend=c("Normal", "Tumor", "Stroma"), fill=c("red", "blue", "green"))
dev.off()


## Annot provides annotation for 79 samples
## mm lists 34 tumor and stroma pairs, of which the 34 are a subset 
## GSE18520  and annot$Sample (Annot[,2])
##  Interestingly enough 3 of the mm sample have "failed" status in annot

## Normal stroma
## annot[annot$"Sample Type"=="Normal Stroma",] 



###################
## Functions
####################

procAnnot<-function(annot) {
  colnames(annot) <-apply(annot[1,], 2, as.character)
  annot<-annot[-1,]
  annot[1:2,]
  ## Subset to the set of stroma samples which are in GSE1850 or are normal stroma
  annot$In_GSE1850<-as.character(annot[,2])%in% pData(GSE18520_eset)[,1]
  annot<-annot[annot$"Sample Type"=="Normal Stroma" | annot$In_GSE1850,]
  
  ## In mm - 34 matched set, 10 normal.. Should have at least 78
  ## GSE18520_eset has 63 samples, but only 44 have stroma
  ## birrer.frma has 135 samples, all stroma, tumor, normal
  ## comparelists(annot[,2], pData(GSE18520_eset)[,1]) gives intersect of 44.  Max samples is 44x2, +10 which is 98... 
  ## in annot there are 45. Stroma samples 332 has matching records, except that its Status is either Completed or Failed???. There is an array for it. Will remove the duplicate
  
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
  
####$$$$$$$$$$$$$$$$###
  ## NOTE####
   
  ##  Not sure if there should be cel files for 312 and 715
  ##  mm thinks there should be
  # Array.Names_Epi Array.Names_Stroma  ID
  # 3          312_Epi         312_Stroma 312
  # 7         715_Epi         715_Stroma 715
  
  ## tHERE ARE IN ANNOT BUT LISTED AS FAILED.
  # array_id alt_sample_name  sample_type days_to_death chemo_status
  # 57   SM-312             312 Tumor Stroma          48          Sen
  # 66   SM-715             715 Tumor Stroma          33          Sen
  
  # matching_epithelial array status definitions in_gse1850
  # 57                       Yes Failed                   TRUE
  # 66                       Yes Failed                   TRUE
  
  ## Can't find cel file name to match
  
  
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
