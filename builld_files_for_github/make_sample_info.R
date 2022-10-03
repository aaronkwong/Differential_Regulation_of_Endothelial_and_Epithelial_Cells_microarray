
#setwd("C:/Users/Aaron/Desktop/CEL")
setwd(gdpath("Masters/Project 5 - Epi and Endo IR/Raw Data/Cell Culture Pre Processing by Aaron/Cel files only name modified"))

#make a folder with only the cell files and set it as the working directory
filenames   <- dir(path="C:/Users/wonga/OneDrive - UHN/Masters/Project5 Epi and Endo IR/Raw Data/Cell Culture Pre Processing by Aaron/Cel files only name modified",pattern = "CEL$");
sample.info <- data.frame(
	row.names = filenames
	
)


sample.info           <- as.data.frame(cbind(filenames, do.call(rbind, strsplit(filenames, split = "_|-|\\."))))
sample.info$V8        <- rep(0,nrow(sample.info))
colnames(sample.info) <- c("filename", "array", "scan.date", "scan.id", "repl", "cat", "cell", "treat")

#cell type list
cell.type<-list()
for (i in 1:nrow(sample.info)){
	if (sample.info[i,"cat"]%in% list(1,2,3,4,5)){
		a<-"endo"
	}else if (sample.info[i,"cat"]%in% list(6,7,8,9,10)){
		a<-"epi"
	}
	cell.type[i]<-a
}
sample.info$cell<-as.character(cell.type)

#treatment list
treatment<-list()
for (i in 1:nrow(sample.info)){
	if (sample.info[i,"cat"]%in% list(1,6)){
		a<-"con"
	}else if (sample.info[i,"cat"]%in% list(2,7)){
		a<-"c6"
	}else if (sample.info[i,"cat"]%in% list(3,8)){
		a<-"c6r2"
	}else if (sample.info[i,"cat"]%in% list(4,9)){
		a<-"c18"
	}else if (sample.info[i,"cat"]%in% list(5,10)){
		a<-"c18r2"
	}
	treatment[i]<-a
}
sample.info$treat<-as.character(treatment)

sample.info$cell      <- factor(sample.info$cell);
sample.info$treat     <- factor(sample.info$treat)
sample.info           <- sample.info[order(sample.info$cat, sample.info$repl), ]
rownames(sample.info) <- paste(sample.info$cell, sample.info$treat, sample.info$repl, sep = ".")

write.table(
	sample.info,
	file=gdpath("Masters/Data for Others/Gaowa/Gaowa_compare_cell_types/github/data/cell_culture_sample_info.tab"),
	col.names=TRUE,
	sep="\t",
	row.names=TRUE
	)


