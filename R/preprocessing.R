# R package management
library(renv)
renv::activate()

suppressPackageStartupMessages(library(oligo))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(org.Hs.eg.db))

dir.create("./results")

sample.info<-read.delim("./data/cell_culture_sample_info.tab",row.names=1)


# lets install the annotation and probe set database (these are binary packages)
install.packages("./brain_array/hugene20sthsentrezg.db_19.0.0.tar.gz",repos=NULL,type="source")
install.packages("./brain_array/pd.hugene20st.hs.entrezg_19.0.0.tar.gz",repos=NULL,type="source")
library(hugene20sthsentrezg.db)
library(pd.hugene20st.hs.entrezg)

#read in cell files
raw.data <- read.celfiles(filenames = paste0("./cel_files/",sample.info$filename), sampleNames = rownames(sample.info), pkgname = "pd.hugene20st.hs.entrezg");



# perform RMA and extract values
rma.data <- rma(raw.data)
rma.expr <- exprs(rma.data)


# save normalized expression
save(sample.info, file = paste0("./results/","sample.info.rdata"))
save(rma.expr, file = paste0("./results/","rma.expr - rma expression matrix - BEAS2B (brainarray annotation).rdata"))

#lets build an annotation file
annotations<-data.frame(select(hugene20sthsentrezg.db,keys=keys(hugene20sthsentrezg.db),columns=c("ENTREZID","SYMBOL","GENENAME")))
rownames(annotations)<-annotations[,1]
save(annotations, file = "./results/annotations_for_cell_culture_v19.rdata")




