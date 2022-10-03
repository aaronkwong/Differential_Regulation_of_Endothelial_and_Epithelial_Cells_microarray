# This script was authord by both Ricardo Zamel and Aaron Wong

# limma analysis
library(limma)

# load rma expression
base::load("./results/rma.expr - rma expression matrix - BEAS2B (brainarray annotation).rdata")
# load sample info
base::load("./results/sample.info.rdata")
# load annotation
base::load("./results/annotations_for_cell_culture_v19.rdata")


# I'll use an as-one-factor limma parametization and extract contrasts of interest
treat      <- sample.info$treat
cell       <- sample.info$cell
treat.cell <- factor(paste0(treat, ".", cell))

#create design matrix
design <- model.matrix(~0+treat.cell)
colnames(design) <- levels(treat.cell)
fit <- lmFit(rma.expr, design)

# add contrats for all the comparisons we want to make
cont.matrix <- makeContrasts(
	con.epi.endo     = con.endo - con.epi,
	c6.epi.endo      = c6.endo - c6.epi,
	c18.epi.endo     = c18.endo - c18.epi,
	c6r2.epi.endo    = c6r2.endo - c6r2.epi,
	c18r2.epi.endo   = c18r2.endo - c18r2.epi,
	epi.con.c6       = c6.epi - con.epi,
	epi.con.c18      = c18.epi - con.epi,
	epi.con.c6r2     = c6r2.epi - con.epi,
	epi.con.c18r2    = c18r2.epi - con.epi,
	epi.c6.c18       = (c18.epi - con.epi) - (c6.epi - con.epi),
	epi.c6r2.c18r2   = (c18r2.epi - con.epi) - (c6r2.epi - con.epi),
	endo.con.c6      = c6.endo - con.endo,
	endo.con.c18     = c18.endo - con.endo,
	endo.con.c6r2    = c6r2.endo - con.endo,
	endo.con.c18r2   = c18r2.endo - con.endo,
	endo.c6.c18      = (c18.endo - con.endo) - (c6.endo - con.endo),
	endo.c6r2.c18r2  = (c18r2.endo - con.endo) - (c6r2.endo - con.endo),
	levels           = design
);

# make contrasts
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#extract p values and FDR correction
fit2.p       <- fit2$p.value
fit2.p05     <- fit2.p <= 0.05
fit2.q       <- apply(fit2.p, 2, p.adjust, method = "fdr")
fit2.q05     <- fit2.q <= 0.05

# summarize 
means <- as.data.frame(sapply(as.character(levels(treat.cell)), function (x) {
	y <- paste(strsplit(x, split = "\\.")[[1]], sep=".")
	rowMeans(rma.expr[, intersect(grep(y[1], colnames(rma.expr), fixed = TRUE), grep(y[2], colnames(rma.expr), fixed = TRUE))]);
}))

#calculate logFC
log2r <- data.frame(
	con.epi.endo     = means$con.endo   - means$con.epi,
	c6.epi.endo      = means$c6.endo    - means$c6.epi,
	c18.epi.endo     = means$c18.endo   - means$c18.epi,
	c6r2.epi.endo    = means$c6r2.endo  - means$c6r2.epi,
	c18r2.epi.endo   = means$c18r2.endo - means$c18r2.epi,
	epi.con.c6       = means$c6.epi     - means$con.epi,
	epi.con.c18      = means$c18.epi    - means$con.epi,
	epi.con.c6r2     = means$c6r2.epi   - means$con.epi,
	epi.con.c18r2    = means$c18r2.epi  - means$con.epi,
	epi.c6.c18       = means$c18.epi    - means$c6.epi,
	epi.c6r2.c18r2   = means$c18r2.epi  - means$c6r2.epi,
	endo.con.c6      = means$c6.endo    - means$con.endo,
	endo.con.c18     = means$c18.endo   - means$con.endo,
	endo.con.c6r2    = means$c6r2.endo  - means$con.endo,
	endo.con.c18r2   = means$c18r2.endo - means$con.endo,
	endo.c6.c18      = means$c18.endo   - means$c6.endo,
	endo.c6r2.c18r2  = means$c18r2.endo - means$c6r2.endo
)


#create a summary table
colnames(log2r) <- paste0(colnames(log2r), ".log2r")
colnames(fit2.p)     <- paste0(colnames(fit2.p), ".p")
colnames(fit2.q)     <- paste0(colnames(fit2.q), ".q")
limma.summ <- cbind(fit2.p, fit2.q, log2r)

# add annotation
limma.anno <- cbind(limma.summ, annotations[rownames(limma.summ), ])

# save summarized object
save(limma.anno, file = "./results/limma analysis all contrasts.rdata")

# write out a table for easy excel sorting
write.table(limma.anno, file = "./results/limma analysis all contrasts.tab",sep="\t")







