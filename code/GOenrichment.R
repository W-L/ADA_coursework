library(data.table)
library(topGO)
library(GOplot)

# parsing functions
tidy <- function(string){
  new <- strsplit(string, split = 'gene:')[[1]][2]
  return(new)
}

split_GOs <- function(GOstring){
  GO_list <- strsplit(GOstring, split=',')
  return(GO_list[[1]])
}

mod_name <- function(name){
  new_name <- strsplit(name, split="[.]")[[1]][2]
  end_name <- gsub('__mRNA', '', new_name)
}

add_genes_to_terms <- function(genes2GO, results){
  for (i in seq(1, nrow(results))){
    GOterm <- results$GO.ID[i]
    results[i, "genes"] <- paste(genes2GO[GOterm][[1]], collapse=', ')
  }
  return(results)
}


# load data and filt for significance
male <- fread("/Users/lweilguny/Desktop/report/DGE/tables/msvsmm.complete.txt")
female <- fread("/Users/lweilguny/Desktop/report/DGE/tables/fmvsfs.complete.txt")
gender <- fread("/Users/lweilguny/Desktop/report/DGE/tables/msvsfs.complete.txt")
male_sign <- subset(male, padj<0.01)
female_sign <- subset(female, padj<0.01)
gender_sign <- subset(gender, padj<0.01)

# store just the gene names
male_genes <- male_sign$Id
female_genes <- female_sign$Id
gender_genes <- gender_sign$Id

# select genes
maturation <- union(male_genes, female_genes) 
gender_effect <- intersect(union(male_genes, female_genes), gender_genes)
genes_of_interest <- maturation[!maturation %in% gender_effect]

# reduce to interesting columns
cand_male <- subset(male_sign, Id %in% genes_of_interest, select=c("Id", "baseMean", "log2FoldChange", "padj"))
cand_female <- subset(female_sign, Id %in% genes_of_interest, select=c("Id", "baseMean", "log2FoldChange", "padj"))
cands <- rbind(cand_male, cand_female)
cands_sort <- cands[order(cands[,'Id'],cands[,'padj']),]
candidates <-  cands_sort[!duplicated(cands_sort$Id),]
candidates$Id <- sapply(candidates$Id, tidy)

# load annotations from eggNOG
annotations <- fread("~/Desktop/report/smansoni_proteome.fasta.unwrapped.mod.emapper.annotations")
gene_ann <- annotations[,c(2,6)]
gene_ann$V6 <- sapply(gene_ann$V6, split_GOs)
gene_ann$V2 <- sapply(gene_ann$V2, mod_name)

# get GO terms into format for topGO
GOlist <- list()
for (i in seq(1,length(gene_ann$V6))){
  GOlist <- c(GOlist, gene_ann$V6[i])
  names(GOlist)[[i]] <- gene_ann$V2[i]
}
str(GOlist)

# define all annotations as gene universe
# and select genes of interest - named factor vector
universe <- gene_ann$V2
geneList <- factor(as.integer(universe %in% candidates$Id))
names(geneList) <- universe

# Go analysis
# initialize topGOdata object
GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, annotationFun = annFUN.gene2GO, gene2GO = GOlist)
GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = geneList, annotationFun = annFUN.gene2GO, gene2GO = GOlist)
GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = geneList, annotationFun = annFUN.gene2GO, gene2GO = GOlist)

# define stat object with Fisher's exact and run test
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "fisher")
res_BP <- getSigGroups(GOdata_BP, test.stat)
res_MF <- getSigGroups(GOdata_MF, test.stat)
res_CC <- getSigGroups(GOdata_CC, test.stat)
hist(score(res_BP), 50, xlab = "p-values", main="Distribution of p-values for GO 'Biological Process'")
hist(score(res_MF), 50, xlab = "p-values", main="Distribution of p-values for GO 'Molecular Function'")
hist(score(res_CC), 50, xlab = "p-values", main="Distribution of p-values for GO 'Cellular Compartment'")

res_table_BP <- GenTable(GOdata_BP, pval = res_BP,  orderBy = "pval", ranksOf = "pval", topNodes = 1000)
res_table_MF <- GenTable(GOdata_MF, pval = res_MF,  orderBy = "pval", ranksOf = "pval", topNodes = 1000)
res_table_CC <- GenTable(GOdata_CC, pval = res_CC,  orderBy = "pval", ranksOf = "pval", topNodes = 1000)

res_table_BP$pval <- as.numeric(res_table_BP$pval)
res_table_MF$pval <- as.numeric(res_table_MF$pval)
res_table_CC$pval <- as.numeric(res_table_CC$pval)

printGraph(GOdata_BP, res_BP, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)


# add annotated genes to results table
gt_BP <- genesInTerm(GOdata_BP)
gt_MF <- genesInTerm(GOdata_MF)
gt_CC <- genesInTerm(GOdata_CC)

results_BP <- add_genes_to_terms(genes2GO = gt_BP, results = res_table_BP)
results_MF <- add_genes_to_terms(genes2GO = gt_MF, results = res_table_MF)
results_CC <- add_genes_to_terms(genes2GO = gt_CC, results = res_table_CC)

results_BP$Category <- rep("BP", nrow(results_BP))
results_MF$Category <- rep("MF", nrow(results_BP))
results_CC$Category <- rep("CC", nrow(results_BP))

full_res <- rbind(results_BP, results_MF, results_CC)
full_res$genes <- as.factor(full_res$genes)

# bubble plot
# rename some columns to create plot object
names(full_res) <- c("ID", "term", "Annotated", "Significant", "Expected", "adj_pval", "genes", "category")
names(candidates) <- c("ID", "baseMean", "logFC", "adj_pval")

terms <- full_res
genes <- candidates

circle_dat_mod <- function(terms, genes){
  
  colnames(terms) <- tolower(colnames(terms))
  terms$genes <- toupper(terms$genes)
  genes$ID <- toupper(genes$ID)
  tgenes <- strsplit(as.vector(terms$genes), ', ')
  if (length(tgenes[[1]]) == 1) tgenes <- strsplit(as.vector(terms$genes), ',')
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  logFC <- sapply(unlist(tgenes), function(x) genes$logFC[match(x, genes$ID)])
  if(class(logFC) == 'factor'){
    logFC <- gsub(",", ".", gsub("\\.", "", logFC))
    logFC <- as.numeric(logFC)
  }
  s <- 1; zsc <- c()
  for (c in 1:length(count)){
    value <- 0
    e <- s + count[c] - 1
    value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1))
    value <- value[!is.na(value)]
    zsc <- c(zsc, sum(value) / sqrt(count[c]))
    s <- e + 1
  }
  if (is.null(terms$id)){
    df <- data.frame(category = rep(as.character(terms$category), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }else{
    df <- data.frame(category = rep(as.character(terms$category), count), ID = rep(as.character(terms$id), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
  return(df)
}

GO_plot <- circle_dat_mod(full_res, candidates)
GOBubble(data=GO_plot, labels = 3, display="single", ID = T, table.legend = T, table.col = T)



