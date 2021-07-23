library(SummarizedExperiment)
library(TCGAbiolinks)

qry.rna <- GDCquery(project = "NCICCR-DLBCL",
                    data.category= "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts")
GDCdownload(qry.rna)

dat <- qry.rna[[1]][[1]]
table(as.factor(dat$sample_type))

rnas.raw <- GDCprepare(qry.rna, summarizedExperiment = TRUE)

data <- assay(rnas.raw)
rownames(data) <- rowData(rnas.raw)$external_gene_name
head(rownames(data))

dataFilt <- TCGAanalyze_Filtering(tabDF = data,
                                  method = "quantile",
                                  qnt.cut = 0.25)
dim(dataFilt)

dataNorm <- TCGAanalyze_Normalization(tabDF = dataFilt, geneInfo = geneInfo)
dim(dataNorm)
