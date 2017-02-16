library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(dplyr)
library(GOstats)

# initializing biomaRt reference
mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart <- useDataset("hsapiens_gene_ensembl",mart)

# getting file names (each of which has a gene list of interest)
path = "~/Documents/git/unified_gene_expression/results/rpe_v_retina/"
file.names <- dir(path, pattern =".txt")

# background genes
background = read.delim(file = paste0(path, "background/GeneIDs-all.txt"))
names(background) = "Gene"

# getting entrez IDs for background genes
background.entrez = getBM(attributes=c("external_gene_name", "entrezgene"),filters="external_gene_name",values=background, mart=mart)$entrezgene

### DF to hold GO annotations ###
GO.enrich = data.frame(matrix(nrow = 0, ncol = 9))
colnames(GO.enrich) = c("module_color", "testDirection", "GOBPID", "Pvalue", "OddsRatio", "ExpCount", "Count", "Size", "Term")

for(geneList in file.names){
  genes = read.delim(file = paste0(path, geneList), header = F)
  names(genes) = "Gene"
  
  # getting entrez IDs for genes of interest
  entrez_ids = getBM(attributes=c("external_gene_name", "entrezgene"),filters="external_gene_name",values=genes, mart=mart)$entrezgene
  
  # Testing over-enrichment
  # genes given must be unique or GOstats won't work
  params.over <- new('GOHyperGParams',
                geneIds=unique(entrez_ids),
                universeGeneIds=unique(background.entrez),
                ontology='BP',
                pvalueCutoff=0.01,
                conditional=F,
                testDirection='over',
                annotation="org.Hs.eg.db"
  )
  hgOver <- hyperGTest(params.over)
  
  result.over <- summary(hgOver)
  
  new.result = bind_cols(data.frame(module_color = rep(gsub(".txt", "", gsub("GeneIDs-", "", geneList)), nrow(result.over)),
                                    testDirection = rep("over", nrow(result.over))), result.over)
  
  GO.enrich = bind_rows(GO.enrich, new.result)
  
}

save(GO.enrich, file = "~/Documents/git/unified_gene_expression/results/rpe_v_retina/string_networks/GO_enrichment.Rdata")





# do no p-value cut-off, see how many tests get returned, and use that for bonferonni

# p-value adjustment
p.adjust(p, method = "bonferonni", n = length(p))


