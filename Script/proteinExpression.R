if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
## install hpar and biomaRt
BiocManager::install("hpar")
BiocManager::install("biomaRt")
## load libraries
library("hpar")
library("biomaRt")
require(tidyverse)

#####
# Supporting functions

data_summary <- function(data, varname, groupnames) {
  require(plyr)
  summary_func <- function(x, col) {
    c(mean = mean(x[[col]], na.rm = TRUE),
      sd = sd(x[[col]], na.rm = TRUE))
    }
  data_sum <- ddply(data, groupnames, .fun = summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#####
# Input for analysis

# An example id list for TPD52-like proteins
ids <- (c("ENSG00000076554",
          "ENSG00000111907",
          "ENSG00000101150",
          "ENSG00000170777"))

# Or using biomaRt to wrangle a list of ensembl ids
my_genes_df <- read.delim("data/rabs.tsv")
my_entrez_ids <- as.character(my_genes_df$Human)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl") 
my_ensembl_ids <- getBM(filters="entrezgene_id",
                        attributes=c("ensembl_gene_id","entrezgene_id"),
                        values=my_entrez_ids,
                        mart=ensembl)
ids <- my_ensembl_ids$ensembl_gene_id

# or retrieve Ensembl IDs for a GO Term
ids <- getBM(attributes = 'ensembl_gene_id', 
             filters = 'go', 
             values = 'GO:0016192', 
             mart = ensembl)

#####
# Process data

# retrieve RNA expression data for cell lines for all ids
gcl <- lapply(ids, getHpa, hpadata = "rnaGeneCellLine")
gcl <- Reduce(rbind,gcl)

# load look-up table for cell lines and match origins to cell lines
cellLineLUT <- read_delim("Data/CellLines.txt",delim = "\t", show_col_types = FALSE)
gcl <- merge(gcl,cellLineLUT)
gcl$Origin <- as.factor(gcl$Origin)

# there is a sorting issue with gene names:
# 1) sorting will place e.g. Rab11 before Rab2
# 2) AP3 needs to come before Rab1
geneNameList <- unique(gcl$Gene.name)
firstLetterList <- substr(geneNameList,1,2)
matchList <- regmatches(geneNameList, regexec('?[0-9]+', geneNameList))
geneNumberList <- unlist({matchList[lengths(matchList)==0] <- length(matchList); matchList})
df <- data.frame(Gene.name = geneNameList,
                 First.letter = firstLetterList,
                 Gene.number = as.numeric(geneNumberList))
df <- df[order(df[,2],df[,3],df[,1]),]
orderedGeneNameList <- as.factor(df$Gene.name)

# adjust plot heights depending on rough size of query
if(length(geneNameList) <= 10) {
  heightVar <- 128
} else if(length(geneNameList) <= 60) {
  heightVar <- 257
} else {
  heightVar <- 512
}

# make the large plot to compare expression of each gene in cell lines
p1 <-  gcl %>%
  mutate(across(Gene.name, factor, levels = orderedGeneNameList)) %>%
  mutate(ordering = as.numeric(Origin), Cell.line = fct_reorder(Cell.line,ordering)) %>%
  ggplot(aes(x = Cell.line, y = pTPM, fill = Origin)) +
  geom_col() +
  theme(text=element_text(size=8), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = "bottom") +
  labs(x = "Cell line", y = "pTPM") +
  facet_wrap(.~Gene.name)
p1
ggsave("Output/Plots/geneCellLines.png", plot = p1, width = 170, height = heightVar, units = "mm")

# now do individual plots of the expression in cell lines
plot_gcl_individual <- function(x,df) {
  iDF <- df[df$Gene.name==x,]
  if(nrow(iDF) == 0) {
    # if there is no data (nrow == 0) we will skip
    return()
  }
  iPlot <-  iDF %>%
    mutate(ordering = as.numeric(Origin), Cell.line = fct_reorder(Cell.line,ordering)) %>%
    ggplot(aes(x = Cell.line, y = pTPM, fill = Origin)) +
    geom_col(show.legend = FALSE) +
    theme(text=element_text(size=6), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),) +
    labs(x = "Cell line", y = "pTPM") +
    facet_wrap(.~Gene.name)
  plotName <- paste0("Output/Plots/geneCellLines_",x,".png")
  ggsave(plotName, plot = iPlot, width = 120, height = 70, units = "mm")
  return()
}
invisible(lapply(orderedGeneNameList, plot_gcl_individual, df = gcl))

# summary data
gcl_summary <- data_summary(gcl,"pTPM","Gene.name")
p1_summary <- ggplot(gcl_summary, aes(x = reorder(Gene.name, pTPM), y = pTPM)) +
  geom_point(fill = "gray") +
  geom_errorbar(aes(ymin = pTPM - sd, ymax = pTPM + sd), width=.2,
                position = position_dodge(.9)) +
  theme(text = element_text(size = 8)) +
  labs(x = "Gene", y = "Mean pTPM") +
  coord_flip()
p1_summary
ggsave("Output/Plots/geneCellLinesSummary.png", plot = p1_summary, width = 170, height = heightVar, units = "mm")


#####
# Protein expression in normal tissue

# retrieve data
nt <- lapply(ids, getHpa, hpadata = "hpaNormalTissue")
nt <- Reduce(rbind,nt)
# expression level needs to be factored and relevelled to make sense
nt$Level <- factor(nt$Level, order = TRUE, levels =c('Not detected', 'Low', 'Medium', 'High'))
nt$Tissue <- as.factor(nt$Tissue)

# make the large plot to show protein expression for all genes
p2 <- nt %>%
  mutate(across(Gene.name, factor, levels = orderedGeneNameList)) %>%
  mutate(ordering = as.numeric(Tissue), Cell.type = fct_reorder(Cell.type,ordering)) %>%
  ggplot(aes(x = Cell.type, y = Level, fill = Tissue)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(text=element_text(size=8), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = "bottom") +
  labs(x = "Tissue", y = "Protein expression") +
  facet_wrap(.~Gene.name)
p2
ggsave("Output/Plots/normalTissue.png", plot = p2, width = 170, height = heightVar, units = "mm")

# now plot individuals of protein expression
plot_nt_individual <- function(x,df) {
  iDF <- df[df$Gene.name==x,]
  if(nrow(iDF) <= 1) {
    # if there is no data (nrow == 0) or only the NA row (nrow == 1), we will skip
    return()
  }
  iPlot <-  df %>%
    mutate(ordering = as.numeric(Tissue), Cell.type = fct_reorder(Cell.type,ordering)) %>%
    filter(Gene.name==x) %>%
    ggplot(aes(x = Cell.type, y = Level, fill = Tissue)) +
    geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
    theme(text=element_text(size=6), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),) +
    labs(x = "Tissue", y = "Protein expression") +
    facet_wrap(.~Gene.name)
  plotName <- paste0("Output/Plots/normalTissue_",x,".png")
  ggsave(plotName, plot = iPlot, width = 120, height = 70, units = "mm")
  return()
}
invisible(lapply(orderedGeneNameList, plot_nt_individual, df = nt))

#####
# Gene expression in normal tissue

# retrieve data
gt <- lapply(ids, getHpa, hpadata = "rnaGeneTissue")
gt <- Reduce(rbind,gt)

# make large plot of gene expression in different tissues, uses pTPM (protein transcripts per million)
p3 <- gt %>%
  mutate(across(Gene.name, factor, levels = orderedGeneNameList)) %>%
  ggplot(aes(x = Tissue, y = pTPM, fill = Tissue)) +
  geom_col() +
  theme(text=element_text(size=8), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = "bottom") +
  labs(x = "Tissue", y = "pTPM") +
  facet_wrap(.~Gene.name)
p3
ggsave("Output/Plots/geneTissue.png", plot = p3, width = 170, height = heightVar, units = "mm")

# now plot individual expression across tissues
plot_gt_individual <- function(x,df) {
  iDF <- df[df$Gene.name==x,]
  if(nrow(iDF) <= 1) {
    # if there is no data (nrow == 0) or only the NA row (nrow == 1), we will skip
    return()
  }
  iPlot <-  iDF %>%
    mutate(across(Gene.name, factor, levels = orderedGeneNameList)) %>%
    ggplot(aes(x = Tissue, y = pTPM, fill = Tissue)) +
    geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
    theme(text=element_text(size=6), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),) +
    labs(x = "Tissue", y = "pTPM") +
    facet_wrap(.~Gene.name)
  plotName <- paste0("Output/Plots/geneTissue_",x,".png")
  ggsave(plotName, plot = iPlot, width = 120, height = 70, units = "mm")
  return()
}
invisible(lapply(orderedGeneNameList, plot_gt_individual, df = gt))

# summary data
gt_summary <- data_summary(gt,"pTPM","Gene.name")
p3_summary <- ggplot(gt_summary, aes(x = reorder(Gene.name, pTPM), y = pTPM)) +
  geom_point(fill = "gray") +
  geom_errorbar(aes(ymin = pTPM - sd, ymax = pTPM + sd), width=.2,
                position = position_dodge(.9)) +
  theme(text = element_text(size = 8)) +
  labs(x = "Gene", y = "Mean pTPM") +
  coord_flip()
p3_summary
ggsave("Output/Plots/geneTissueSummary.png", plot = p3_summary, width = 170, height = heightVar, units = "mm")