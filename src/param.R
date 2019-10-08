#r params
wkdir<-file.path("/Users/aculhane/Documents/GitHub/Schwede_Supplement")

inputdir<-file.path(wkdir,"data", "input")  # Raw Data
datadir<-file.path(wkdir, "data", "Processed")  # Processed, Working Data

srcdir<-file.path(wkdir,"src") # Scripts R code


#Output
figdir<-file.path(wkdir, "fig")
outdir<-file.path(wkdir, "output")
supplDir<-outdir
barcodefile<- file.path(supplDir, "barcodeRes.RData")
curatedOvarianDataEsetsfile=file.path(supplDir, "esets.rda") 

SuppTab <- file.path(supplDir, "Supplementary Tables All.xlsx")

source(file.path(srcdir,"ProcSubtypes.R"))