library(limma)
library(affy)
library(genefilter)
library(ggplot2)
library(heatmap3)
library(readxl)
library(topGO)
library(hgu133plus2.db)

#Loading dataset
Expression <- read_excel('D:/Scar formation/Expression.xlsx')
data <- as.data.frame(Expression)
data <- data.frame(data[,-1], row.names = data[,1])
datamx <- as.matrix(data)

#DEGs Testing - Pairwise Analysis
fac <- factor(c('C','S','C','S','C','S'))
patient <- factor((c(2,2,3,3,4,4)))
design <- model.matrix(~fac)
fit <- lmFit(datamx, design)
colnames(coef(fit))
fit <- eBayes(fit)
tt <- topTable(fit, coef =2)
tt.all<- topTable(fit, coef=2, number="ALL")
con <- makeContrasts(patient2 - patient3, levels = design)
con.fit<- contrasts.fit(fit, con)
con.fit<- eBayes(con.fit.)
colnames(con.fit)
tt.con <- topTable(con.fit, coef = "patient2 - patient3")
write.table(tt.con, "Patient2-Patient3.txt", sep = "\t", quote = FALSE)

#Top GO Analysis
genes <- tt.all[,5]
names(genes) <- rownames(tt.all)
topDiffGenes <- function(allScore) {return(allScore < 0.05)}
sum(topDiffGenes(genes))
sampleGOdata <- new("topGOdata",description = "Scar formation BP", ontology="BP",allGenes = genes, geneSel =topDiffGenes,nodeSize=10,annot=annFUN.db, affyLib = "hgu133plus2.db")
sampleGOdata
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic ="fisher")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic ="ks")
resultKS.elim <- runTest(sampleGOdata, algorithm= "elim", statistic ="ks")
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, classicKS = resultKS,elimKS = resultKS.elim, orderBy= "elimKS", ranksOf="classicFisher", topNodes = 20)
write.table(allRes, "Top All GO.txt", sep = "\t", quote = FALSE)
pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(upGOdata,names(pValue.classic))
gSize <- gstat$Annotated/max(gstat$Annotated)*4
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize, col = gCol)
#Visualizing
sel.go <- names(pValue.classic)[pValue.elim<pValue.classic]
cbind(termStat(sampleGOdata,sel.go), elim = pValue.elim[sel.go], classic = pValue.classic[sel.go])
showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')

#UPREGULATED GENES ONTOLOGY (SCORE BASED ON ADJUSTED P-VALUE < 0.05)
logFC <- tt.all[,1]
pvals <- tt.all[,5]

filter <- data.frame(which(logFC>0 & pvals <0.05))
a <- filter[,1]
upGOsample <- as.matrix(x=tt.all[a,])
genes.up <- upGOsample[,5]
names(genes.up) <- rownames(upGOsample)
upGOdata <- new("topGOdata",
                description = "Upregulated genes",
                ontology = "BP",
                allGenes = genes.up,
                geneSel = topDiffGenes,
                nodeSize = 5,
                annot = annFUN.db,
                affyLib = "hgu133plus2.db")
upGOdata

resultFisher.up <- runTest(upGOdata, algorithm = "classic", statistic = "fisher")
resultKS.up <- runTest(upGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim.up <- runTest(upGOdata, algorithm ="elim", statistic = "ks")
upRes <- GenTable(upGOdata, classicFisher = resultFisher.up, classicKS = resultKS.up, elimKS = resultKS.elim.up, orderBy = "elimKS", ranksOf = "classicFisher",topNodes = 100)
pValue.classic <- score(resultKS.up)
pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(upGOdata,names(pValue.classic))
sel.go <- names(pValue.classic)[pValue.elim<pValue.classic]
cbind(termStat(upGOdata,sel.go), elim = pValue.elim[sel.go], classic = pValue.classic[sel.go])
showSigOfNodes(upGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')

#DOWNREGULATED GENES ONTOLOGY
filter <- data.frame(which(logFC<0 & pvals <0.05))
a <- filter[,1]
downGOsample <- as.matrix(x=tt.all[a,])
genes.down <- downGOsample[,5]
names(genes.down) <- rownames(downGOsample)
downGOdata <- new("topGOdata",
                description = "Downregulated genes P<0.05",
                ontology = "BP",
                allGenes = genes.down,
                geneSel = topDiffGenes,
                nodeSize = 5,
                annot = annFUN.db,
                affyLib = "hgu133plus2.db")
downGOdata

resultFisher.down <- runTest(downGOdata, algorithm = "classic", statistic = "fisher")
resultKS.down <- runTest(downGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim.down <- runTest(downGOdata, algorithm ="elim", statistic = "ks")
downRes <- GenTable(downGOdata, classicFisher = resultFisher.down, classicKS = resultKS.down, elimKS = resultKS.elim.down, orderBy = "elimKS", ranksOf = "classicFisher",topNodes = 20)
pValue.classic.down <- score(resultKS.down)
pValue.elim.down <- score(resultKS.elim.down)[names(pValue.classic.down)]
gstat.down <- termStat(downGOdata,names(pValue.classic.down))
sel.go.down <- names(pValue.classic.down)[pValue.elim.down<pValue.classic.down]
cbind(termStat(downGOdata,sel.go.down), elim = pValue.elim.down[sel.go.down], classic = pValue.classic.down[sel.go.down])
showSigOfNodes(downGOdata, score(resultKS.elim.down), firstSigNodes = 5, useInfo = 'all')


require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
ontogoly <- getBM(mart=mart, attributes=c("affy_hg_u133_plus_2","hgnc_symbol","ensembl_gene_id", 
                                          "gene_biotype", "external_gene_name","go_id",
                                          "name_1006","namespace_1003"), 
                  filter="affy_hg_u133_plus_2", values=rownames(upGOsample), uniqueRows=TRUE)
ontology.matched <- match(rownames(upGOsample),ontogoly$affy_hg_u133_plus_2)
GOlist <- data.frame(rownames(upGOsample),ontogoly[ontology.matched,c("affy_hg_u133_plus_2","gene_biotype","hgnc_symbol")])
rownames(GOlist) <- GOlist[,1]
GOlist <- GOlist[,-1]
