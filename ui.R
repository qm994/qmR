library(shiny)
library(VennDiagram)
library(stringr)
library(d3heatmap)
library(RColorBrewer)
library(DT)
library(KEGGREST)#bioconductor
library(dplyr)
library(rsconnect)

deSets <- vector(mode = "list")

all_C1 <- read.table("C1-MIA.allDE.tab", sep = "\t", header = TRUE)
all_C1$probeID <- NULL
all_C1$AveExpr <- NULL
all_C4 <- read.table("C4-MIA.allDE.tab", sep = "\t", header = TRUE)
all_C4$probeID <- NULL
all_C4$AveExpr <- NULL
temp_C1_up <- filter(all_C1, logFC > 0)
temp_C1_up <- temp_C1_up[complete.cases(temp_C1_up), ]
temp_C1_down <- filter(all_C1, logFC < 0)
temp_C1_down <- temp_C1_down[complete.cases(temp_C1_down), ]
temp_C4_up <- filter(all_C4, logFC > 0)
temp_C4_up <- temp_C4_up[complete.cases(temp_C4_up), ]
temp_C4_down <- filter(all_C4, logFC < 0)
temp_C4_down <- temp_C4_down[complete.cases(temp_C4_down), ]

deSets[["C1-MIA.up"]] <- temp_C1_up
deSets[["C1-MIA.down"]] <- temp_C1_down
deSets[["C4-MIA.up"]] <- temp_C4_up
deSets[["C4-MIA.down"]] <- temp_C4_down

Gem_all <- read.table("GEM_HuR-GEM_IgG.allDE.tab", sep = "\t", header = TRUE)
Gem_all$probeID <- NULL
Gem_all$AveExpr <- NULL
MIA_all <- read.table("MIA_HuR-MIA_IgG.allDE.tab", sep = "\t", header = TRUE)
MIA_all$probeID <- NULL
MIA_all$AveExpr <- NULL
OLAP_all <- read.table("OLAP_HuR-OLAP_IgG.allDE.tab", sep = "\t", header = TRUE)
OLAP_all$probeID <- NULL
OLAP_all$AveExpr <- NULL

Gem_up <- filter(Gem_all, logFC > 0)
Gem_up <- Gem_up[complete.cases(Gem_up), ]

Gem_down <- filter(Gem_all, logFC < 0)
Gem_down <- Gem_down[complete.cases(Gem_down), ]

MIA_up <- filter(MIA_all, logFC > 0)
MIA_up <- MIA_up[complete.cases(MIA_up), ]

MIA_down <- filter(MIA_all, logFC < 0)
MIA_down <- MIA_down[complete.cases(MIA_down), ]

OLAP_up <- filter(OLAP_all, logFC > 0)
OLAP_up <- OLAP_up[complete.cases(OLAP_up), ]

OLAP_down <- filter(OLAP_all, logFC < 0)
OLAP_down <- OLAP_down[complete.cases(OLAP_down), ]

deSets[["GEM_HuR-GEM_IgG.up"]] <- Gem_up
deSets[["GEM_HuR-GEM_IgG.up"]] <- deSets[["GEM_HuR-GEM_IgG.up"]][!(deSets[["GEM_HuR-GEM_IgG.up"]]$gene %in% "Y_RNA"), ]
deSets[["GEM_HuR-GEM_IgG.down"]] <- Gem_down
deSets[["MIA_HuR-MIA_IgG.up"]] <- MIA_up
deSets[["MIA_HuR-MIA_IgG.up"]] <- deSets[["MIA_HuR-MIA_IgG.up"]][!(deSets[["MIA_HuR-MIA_IgG.up"]]$gene %in% "Y_RNA"), ]
deSets[["MIA_HuR-MIA_IgG.down"]] <- MIA_down
deSets[["OLAP_HuR-OLAP_IgG.up"]] <- OLAP_up
deSets[["OLAP_HuR-OLAP_IgG.up"]] <- deSets[["OLAP_HuR-OLAP_IgG.up"]][!(deSets[["OLAP_HuR-OLAP_IgG.up"]]$gene %in% "Y_RNA"), ]
deSets[["OLAP_HuR-OLAP_IgG.down"]] <- OLAP_down
# RNA seq loading
rnaseq.mia.c.hvl <- read.table("rnaseq.MIA_C_HvL.sigDE.txt", sep = "\t", header = TRUE)
rnaseq.mia.c.hvl$ensembl <- NULL
names(rnaseq.mia.c.hvl)[1] <- "gene"

rnaseq.mia.c.hvl.up <- filter(rnaseq.mia.c.hvl, logFC > 0)
rnaseq.mia.c.hvl.up <- rnaseq.mia.c.hvl.up[complete.cases(rnaseq.mia.c.hvl.up), ]
rnaseq.mia.c.hvl.down <- filter(rnaseq.mia.c.hvl, logFC < 0)
rnaseq.mia.c.hvl.down <- rnaseq.mia.c.hvl.down[complete.cases(rnaseq.mia.c.hvl.down), ]

rnaseq.mia.k.hvl <- read.table("rnaseq.MIA_K_HvL.sigDE.txt", sep = "\t", header = TRUE)
rnaseq.mia.k.hvl$ensembl <- NULL
names(rnaseq.mia.k.hvl)[1] <- "gene"

rnaseq.mia.k.hvl.up <- filter(rnaseq.mia.k.hvl, logFC > 0)
rnaseq.mia.k.hvl.up <- rnaseq.mia.k.hvl.up[complete.cases(rnaseq.mia.k.hvl.up), ]
rnaseq.mia.k.hvl.down <- filter(rnaseq.mia.k.hvl, logFC < 0)
rnaseq.mia.k.hvl.down <- rnaseq.mia.k.hvl.down[complete.cases(rnaseq.mia.k.hvl.down), ]

rnaseq.mia.h.cvk <- read.table("rnaseq.MIA_H_CvK.sigDE.txt", sep = "\t", header = TRUE)
rnaseq.mia.h.cvk$ensembl <- NULL
names(rnaseq.mia.h.cvk)[1] <- "gene"

rnaseq.mia.h.cvk.up <- filter(rnaseq.mia.h.cvk, logFC > 0)
rnaseq.mia.h.cvk.up <- rnaseq.mia.h.cvk.up[complete.cases(rnaseq.mia.h.cvk.up), ]
rnaseq.mia.h.cvk.down <- filter(rnaseq.mia.h.cvk, logFC < 0)
rnaseq.mia.h.cvk.down <- rnaseq.mia.h.cvk.down[complete.cases(rnaseq.mia.h.cvk.down), ]

rnaseq.mia.l.cvk <- read.table("rnaseq.MIA_L_CvK.sigDE.txt", sep = "\t", header = TRUE)
rnaseq.mia.l.cvk$ensembl <- NULL
names(rnaseq.mia.l.cvk)[1] <- "gene"

rnaseq.mia.l.cvk.up <- filter(rnaseq.mia.l.cvk, logFC > 0)
rnaseq.mia.l.cvk.up <- rnaseq.mia.l.cvk.up[complete.cases(rnaseq.mia.l.cvk.up), ]
rnaseq.mia.l.cvk.down <- filter(rnaseq.mia.l.cvk, logFC < 0)
rnaseq.mia.l.cvk.down <- rnaseq.mia.l.cvk.down[complete.cases(rnaseq.mia.l.cvk.down), ]

deSets[["rnaseq.MIA.C.HvL.up"]] <- rnaseq.mia.c.hvl.up
deSets[["rnaseq.MIA.C.HvL.down"]] <- rnaseq.mia.c.hvl.down
deSets[["rnaseq.MIA.K.HvL.up"]] <- rnaseq.mia.k.hvl.up 
deSets[["rnaseq.MIA.K.HvL.down"]] <- rnaseq.mia.k.hvl.down
deSets[["rnaseq.MIA.H.CvK.up"]] <- rnaseq.mia.h.cvk.up
deSets[["rnaseq.MIA.H.CvK.down"]] <- rnaseq.mia.h.cvk.down
deSets[["rnaseq.MIA.L.CvK.up"]] <- rnaseq.mia.l.cvk.up
deSets[["rnaseq.MIA.L.CvK.down"]] <- rnaseq.mia.l.cvk.down

# the info prepare for heatmap
expressionTable <- read.table("resistantTotal.rawAnno.tab",sep="\t",header = TRUE)

expressionTable <- expressionTable[ , c("annotation","MIA.total.1","MIA.total.2","MIA.total.3","C1.total.1","C1.total.2","C1.total.3",
                                        "C4.total.1","C4.total.2","C4.total.3")]
expressionTable <- expressionTable[!duplicated(expressionTable$annotation), ]

pancRIP_expression <- read.table("pancRIP.rawAnno.tab",sep = "\t", header = TRUE)
pancRIP_expression <- pancRIP_expression[ ,c(2:14)]
pancRIP_expression <- pancRIP_expression[!duplicated(pancRIP_expression$annotation), ]

expression.rnaseq <- read.table("rnaSeqTPM.tab", sep = "\t", header = TRUE)
expression.rnaseq <- expression.rnaseq[ ,c(2:14)]
names(expression.rnaseq)[1] <- "gene"
expression.rnaseq <- expression.rnaseq[!duplicated(expression.rnaseq$gene), ]

# load the text file which contains the DNA repair genes
DNA_repair_genes <- read.csv("dsNonHomoRepair.txt", sep = ",", header = FALSE)
DNA_string_genes <- as.character(unlist(DNA_repair_genes))


#We need to eliminate the Y_RNA siince they only exist in up tables
deSets[["C1-MIA.up"]] <- deSets[["C1-MIA.up"]][!(deSets[["C1-MIA.up"]]$gene %in% "Y_RNA"), ]
deSets[["C4-MIA.up"]] <- deSets[["C4-MIA.up"]][!(deSets[["C4-MIA.up"]]$gene %in% "Y_RNA"), ]

#Extract the gene column as a string
C1UPstringenes <- str_c(unlist(deSets[["C1-MIA.up"]]$gene), sep = ",")
C4UPstringenes <- str_c(unlist(deSets[["C4-MIA.up"]]$gene), sep = ",")
C1DOWNstringenes <- str_c(unlist(deSets[["C1-MIA.down"]]$gene), sep = ",")
C4DOWNstringenes <- str_c(unlist(deSets[["C4-MIA.down"]]$gene), sep = ",")

GemUpstringenes <- str_c(unlist(deSets[["GEM_HuR-GEM_IgG.up"]]$gene), sep = ",")
GemDownstringenes <- str_c(unlist(deSets[["GEM_HuR-GEM_IgG.down"]]$gene), sep = ",")
MIAUpstringenes <- str_c(unlist(deSets[["MIA_HuR-MIA_IgG.up"]]$gene), sep = ",")
MIADownstringenes <- str_c(unlist(deSets[["MIA_HuR-MIA_IgG.down"]]$gene), sep = ",")
OLAPUpstringenes <- str_c(unlist(deSets[["OLAP_HuR-OLAP_IgG.up"]]$gene), sep = ",")
OLAPDownstringenes <- str_c(unlist(deSets[["OLAP_HuR-OLAP_IgG.down"]]$gene), sep = ",")

expressionTable_genes <- str_c(unlist(expressionTable$annotation), sep = ",")
pancRIP_expression_genes <- str_c(unlist(pancRIP_expression$annotation), sep = ",")
expression.rnaseq_genes <- str_c(unlist(expression.rnaseq$gene), sep = ",")

#decide the lowercase objects in gene column(True or False)
C1UPGetLowercase <- str_match(C1UPstringenes, "[a-z]*")
C4UPGetLowercase <- str_match(C4UPstringenes, "[a-z]*")
C1DOWNGetLowercase <- str_match(C1DOWNstringenes, "[a-z]*")
C4DOWNGetLowercase <- str_match(C4DOWNstringenes, "[a-z]*")

GemUpGetLowercase <- str_match(GemUpstringenes, "[a-z]*")
GemDownGetLowercase <- str_match(GemDownstringenes, "[a-z]*")
MIAUpGetLowercase <- str_match(MIAUpstringenes, "[a-z]*")
MIADownGetLowercase <- str_match(MIADownstringenes, "[a-z]*")
OLAPUpGetLowercase <- str_match(OLAPUpstringenes, "[a-z]*")
OLAPDownGetLowercase <- str_match(OLAPDownstringenes, "[a-z]*")


expressionTable_lowercase <- str_match(expressionTable_genes, "[a-z]*")
pancRIP_expression_lowercase <- str_match(pancRIP_expression_genes, "[a-z]*")
expression.rnaseq_lowercase <- str_match(expression.rnaseq_genes, "[a-z]*")

#Eliminate the lowercase gene names
deSets[["C1-MIA.up"]] <- deSets[["C1-MIA.up"]][!(deSets[["C1-MIA.up"]]$gene %in% C1UPGetLowercase), ]
deSets[["C4-MIA.up"]] <- deSets[["C4-MIA.up"]][!(deSets[["C4-MIA.up"]]$gene %in% C4UPGetLowercase), ]
deSets[["C1-MIA.dowm"]] <- deSets[["C1-MIA.down"]][!(deSets[["C1-MIA.down"]]$gene %in% C1DOWNGetLowercase), ]
deSets[["C4-MIA.down"]] <- deSets[["C4-MIA.down"]][!(deSets[["C4-MIA.down"]]$gene %in% C4DOWNGetLowercase), ]

deSets[["GEM_HuR-GEM_IgG.up"]] <- deSets[["GEM_HuR-GEM_IgG.up"]][!(deSets[["GEM_HuR-GEM_IgG.up"]]$gene %in% GemUpGetLowercase), ]
deSets[["GEM_HuR-GEM_IgG.down"]] <- deSets[["GEM_HuR-GEM_IgG.down"]][!(deSets[["GEM_HuR-GEM_IgG.down"]]$gene %in% GemDownGetLowercase), ]
deSets[["MIA_HuR-MIA_IgG.up"]] <- deSets[["MIA_HuR-MIA_IgG.up"]][!(deSets[["MIA_HuR-MIA_IgG.up"]]$gene %in% MIAUpGetLowercase), ]
deSets[["MIA_HuR-MIA_IgG.down"]] <- deSets[["MIA_HuR-MIA_IgG.down"]][!(deSets[["MIA_HuR-MIA_IgG.down"]]$gene %in% MIADownGetLowercase), ]
deSets[["OLAP_HuR-OLAP_IgG.up"]] <- deSets[["OLAP_HuR-OLAP_IgG.up"]][!(deSets[["OLAP_HuR-OLAP_IgG.up"]]$gene %in% OLAPUpGetLowercase), ]
deSets[["OLAP_HuR-OLAP_IgG.down"]] <- deSets[["OLAP_HuR-OLAP_IgG.down"]][!(deSets[["OLAP_HuR-OLAP_IgG.down"]]$gene %in% OLAPDownGetLowercase), ]

expressionTable <- expressionTable[!(expressionTable$annotation %in% expressionTable_lowercase), ]
pancRIP_expression <- pancRIP_expression[!(pancRIP_expression$annotation %in% pancRIP_expression_lowercase), ]
expression.rnaseq <- expression.rnaseq[!(expression.rnaseq$gene %in% expression.rnaseq_lowercase), ]

rownames(expressionTable) <- expressionTable$annotation
rownames(pancRIP_expression) <- pancRIP_expression$annotation
rownames(expression.rnaseq) <- expression.rnaseq$gene
expressionTable$annotation <- NULL
pancRIP_expression$annotation <- NULL
expression.rnaseq$gene <- NULL

testnames1 <- list("C1-MIA", "C4-MIA")
testnames2 <- list("GEM_HuR-GEM_IgG","MIA_HuR-MIA_IgG","OLAP_HuR-OLAP_IgG")
testnames3 <- list("rnaseq.MIA.C.HvL", "rnaseq.MIA.K.HvL", "rnaseq.MIA.H.CvK", "rnaseq.MIA.L.CvK")



# TCA KEGG PATHWAY
gene_list <- list()

genes <- as.list(keggGet("hsa00020")[[1]]$GENE)

a <- 1
for(x in 1:length(genes)){
  if(x%%2 == 0){
    gene_list[a] = genes[x]
    a = a+1
  }
}

gene_string <- unlist(gene_list)
gene_string1 <- unlist(lapply(gene_string, function(x) str_split(x, "; ", 2)))
TCA_gene_string <- as.character()

b <- 1
for(y in 1:length(gene_string1)){
  if(y%%2 != 0){
    TCA_gene_string[b] = gene_string1[y]
    b= b+1
  }
}

# Hedgehog signaling pathway
gene_list_1 <- list()
genes1 <- as.list(keggGet("hsa04340")[[1]]$GENE)

c <- 1
for(x in 1:length(genes1)){
  if(x%%2 == 0){
    gene_list_1[c] = genes1[x]
    c = c+1
  }
}

gene_string_1 <- unlist(gene_list_1)
gene_string2 <- unlist(lapply(gene_string_1, function(x) str_split(x, "; ", 2)))
hedgehog_gene_string <- as.character()
d <- 1

for(y in 1:length(gene_string2)){
  if(y%%2 != 0){
    hedgehog_gene_string[d] = gene_string2[y]
    d= d+1
  }
}

# Cell cycle pathway
original_cell_list <- list()
original_cell_genes <- as.list(keggGet("hsa04110")[[1]]$GENE)

a1 <- 1
for(x in 1:length(original_cell_genes)){
  if(x%%2 == 0){
    original_cell_list[a1] = original_cell_genes[x]
    a1 = a1+1
  }
}

original_cell_string <- unlist(original_cell_list)
new_cell_string <- unlist(lapply(original_cell_string, function(x) str_split(x, "; ", 2)))
CellCycle_gene_string <- as.character()
b1 <- 1

for(y in 1:length(new_cell_string)){
  if(y%%2 != 0){
    CellCycle_gene_string[b1] = new_cell_string[y]
    b1= b1+1
  }
}

#Apoptosis 
original_apop_list <- list()
original_apop_genes <- as.list(keggGet("hsa04210")[[1]]$GENE)

a2 <- 1
for(x in 1:length(original_apop_genes)){
  if(x%%2 == 0){
    original_apop_list[a2] = original_apop_genes[x]
    a2 = a2+1
  }
}

original_apop_string <- unlist(original_apop_list)
new_apop_string <- unlist(lapply(original_apop_string, function(x) str_split(x, "; ", 2)))
Apoptosis_gene_string <- as.character()
b2 <- 1

for(y in 1:length(new_apop_string)){
  if(y%%2 != 0){
    Apoptosis_gene_string[b2] = new_apop_string[y]
    b2= b2+1
  }
}

# Adherens Junction
original_junction_list <- list()
original_junction_genes <- as.list(keggGet("hsa04520")[[1]]$GENE)

a3 <- 1
for(x in 1:length(original_junction_genes)){
  if(x%%2 == 0){
    original_junction_list[a3] = original_junction_genes[x]
    a3 = a3+1
  }
}

original_junction_string <- unlist(original_junction_list)
new_junction_string <- unlist(lapply(original_junction_string, function(x) str_split(x, "; ", 2)))
Adherens_gene_string <- as.character()
b3 <- 1

for(y in 1:length(new_junction_string)){
  if(y%%2 != 0){
    Adherens_gene_string[b3] = new_junction_string[y]
    b3= b3+1
  }
}

# Pancreatic Cancer
original_pancreatic_list <- list()
original_pancreatic_genes <- as.list(keggGet("hsa05212")[[1]]$GENE)

a4 <- 1
for(x in 1:length(original_pancreatic_genes)){
  if(x%%2 == 0){
    original_pancreatic_list[a4] = original_pancreatic_genes[x]
    a4 = a4+1
  }
}

original_pancreatic_string <- unlist(original_pancreatic_list)
new_pancreatic_string <- unlist(lapply(original_pancreatic_string, function(x) str_split(x, "; ", 2)))
Pancreatic_gene_string <- as.character()
b4 <- 1

for(y in 1:length(new_pancreatic_string)){
  if(y%%2 != 0){
    Pancreatic_gene_string[b4] = new_pancreatic_string[y]
    b4= b4+1
  }
}
#ABC TRANSPORTERS
original_ABC_list <- list()
original_ABC_genes <- as.list(keggGet("hsa02010")[[1]]$GENE)

a5 <- 1
for(x in 1:length(original_ABC_genes)){
  if(x%%2 == 0){
    original_ABC_list[a5] = original_ABC_genes[x]
    a5 = a5+1
  }
}

original_ABC_string <- unlist(original_ABC_list)
new_ABC_string <- unlist(lapply(original_ABC_string, function(x) str_split(x, "; ", 2)))
ABC_gene_string <- as.character()
b5 <- 1

for(y in 1:length(new_ABC_string)){
  if(y%%2 != 0){
    ABC_gene_string[b5] = new_ABC_string[y]
    b5 = b5 + 1
  }
}


shinyUI(navbarPage( title = "Choosing the test",
                    tabsetPanel(
                      # tab 1
                      tabPanel("Microarray data", 
                               fluidRow(
                                 column(3, selectInput(inputId = "deTest1",
                                                       label = "Choose Test 1",
                                                       choices = c("C1-MIA","C4-MIA","GEM_HuR-GEM_IgG",
                                                                   "MIA_HuR-MIA_IgG","OLAP_HuR-OLAP_IgG",
                                                                   "rnaseq.MIA.C.HvL", "rnaseq.MIA.K.HvL",
                                                                   "rnaseq.MIA.H.CvK", "rnaseq.MIA.L.CvK"
                                                       ))
                                 ),
                                 column(3, radioButtons("direction1", "Change in expression",
                                                        c("Up"="up", "Down"="down"))
                                 ),
                                 column(3, selectInput("deTest2",
                                                       "Choose Test 2",
                                                       choices = c("C1-MIA", "C4-MIA","GEM_HuR-GEM_IgG",
                                                                   "MIA_HuR-MIA_IgG","OLAP_HuR-OLAP_IgG",
                                                                   "rnaseq.MIA.C.HvL", "rnaseq.MIA.K.HvL",
                                                                   "rnaseq.MIA.H.CvK", "rnaseq.MIA.L.CvK"
                                                       ))
                                 ),
                                 column(3, radioButtons("direction2", "Change in expression",
                                                        c("Up"="up", "Down"="down"))
                                 )
                                 
                               ),
                               
                               
                               fluidRow(
                                 column(3, numericInput("number1", "Choose the number", value = 5)
                                 ),
                                 column(3, tags$strong("download the gene list"), downloadButton("downloadData1", "Download")
                                 ),
                                 column(3, numericInput("number2", "Choose the number", value = 5)
                                 ),
                                 column(3, tags$strong("download the gene list"), downloadButton("downloadData2", "Download")
                                 )
                                 
                                 
                               ),
                               
                               br(),
                               
                               fluidRow(
                                 column(12, wellPanel(
                                   
                                   # when the tests1 are micrarray datasets
                                   conditionalPanel(condition = "['C1-MIA', 'C4-MIA'].indexOf(input.deTest1) != -1 || 
                                                    ['GEM_HuR-GEM_IgG','MIA_HuR-MIA_IgG','OLAP_HuR-OLAP_IgG'].indexOf(input.deTest1) != -1",
                                                    
                                                    actionButton("check_sort_1", "Sort (Yes/No)"),
                                                    
                                                    fluidRow(column(5, dataTableOutput("topDEgenes1")))
                                   ),
                                   
                                   # when the test1 are rna seq datasets
                                   conditionalPanel(condition = "['rnaseq.MIA.C.HvL', 'rnaseq.MIA.K.HvL',
                                                    'rnaseq.MIA.H.CvK', 'rnaseq.MIA.L.CvK'].indexOf(input.deTest1) != -1",
                                                    actionButton("check_sort_3", "Sort (Yes/No)"),
                                                    
                                                    fluidRow(column(5, dataTableOutput("topDEgenes3")))
                                   )
                                 )
                                 )
                               ),
                               fluidRow(
                                 column(12, wellPanel(
                                   
                                   # when test2 are microarray datasets
                                   conditionalPanel(condition = "['C1-MIA', 'C4-MIA'].indexOf(input.deTest2) != -1 || 
                                                    ['GEM_HuR-GEM_IgG','MIA_HuR-MIA_IgG','OLAP_HuR-OLAP_IgG'].indexOf(input.deTest2) != -1",
                                                    
                                                    actionButton("check_sort_2", "Sort (Yes/No)"),
                                                    
                                                    fluidRow(column(5, dataTableOutput("topDEgenes2")), offset = 1)
                                                    
                                   ),
                                   
                                   # when test2 are rna seq datasets
                                   conditionalPanel(condition = "['rnaseq.MIA.C.HvL', 'rnaseq.MIA.K.HvL',
                                                    'rnaseq.MIA.H.CvK', 'rnaseq.MIA.L.CvK'].indexOf(input.deTest2) != -1",
                                                    actionButton("check_sort_4", "Sort (Yes/No)"),
                                                    
                                                    fluidRow(column(5, dataTableOutput("topDEgenes4")), offset = 1)
                                   )
                                 )
                                 )# second column ends
                               ),# fluidrow ends 
                               
                               
                               fluidRow(
                                 column(5, tags$strong("Download the Overlap Genes from venn diagram"), downloadButton("downloadData3", "Download"))
                               ),
                               fluidRow(
                                 column(12, plotOutput("vennPlot")
                                 )
                               )
                               ),
                      
                      #tab2
                      tabPanel("c1-c4-mia",
                               
                               sidebarPanel(
                                 wellPanel(
                                   tags$strong("Updating the Heatmap"),
                                   actionButton("Button", label = "Update heatmap"),
                                   
                                   tags$hr(),
                                   
                                   checkboxGroupInput("heatmap.1", "Select samples to include:",
                                                      choices = colnames(expressionTable)[1:9],
                                                      selected = colnames(expressionTable)[1:9]
                                   ),
                                   
                                   tags$hr(),
                                   
                                   selectInput("check_gene", tags$h5("Choosing the genes to remove the row"),
                                               choices = "", selected = "", multiple = TRUE)
                                   
                                 )),
                               mainPanel(
                                 d3heatmapOutput("GeneHeatmap", height = "400px", width = "100%")
                                 
                               )
                               
                      ),
                      
                      #tab3 is the heatmap based on another expression table
                      tabPanel("GEM-OLAP-MIA",
                               
                               sidebarPanel(
                                 wellPanel(
                                   tags$strong("Updating the heatmap"),
                                   actionButton("Button2", "Updating the Heatmap"),
                                   
                                   checkboxGroupInput("heatmap.2", "Select samples to include:",
                                                      choices = colnames(pancRIP_expression)[1:12],
                                                      selected = colnames(pancRIP_expression)[1:12]
                                   ),
                                   selectInput("check_gene.2", tags$h5("Choosing the genes to remove the row"),
                                               choices = "", selected = "")
                                 )
                                 
                               ),
                               
                               mainPanel(
                                 d3heatmapOutput("GeneHeatmap2", height = "400px", width = "100%")
                               )
                      ),
                      
                      #tab4 is the heatmap based on another expression table
                      tabPanel("rna-mia",
                               
                               sidebarPanel(
                                 wellPanel(
                                   tags$strong("Updating the heatmap"),
                                   actionButton("Button3", "Updating the heatmap"),
                                   checkboxGroupInput("heatmap.3", "Select samples to include:",
                                                      choices = colnames(expression.rnaseq)[1:12],
                                                      selected = colnames(expression.rnaseq)[1:12]
                                   ),
                                   selectInput("check_gene.3", tags$h5("Choosing the genes to remove the row"),
                                               choices = "", selected = "", multiple = TRUE)
                                 )
                               ),
                               
                               mainPanel(
                                 d3heatmapOutput("GeneHeatmap3", height = "400px", width = "100%") 
                               )
                      ),
                      #tab5 combine all the pathways in one selectInput
                      tabPanel(
                        "Pathway Selection",
                        wellPanel(
                          fluidRow(
                            
                            column(5, 
                                   selectInput("pathways", "Choosing the pathway",
                                               choices = c("Apoptosis", "DNA Repair", "Cell Cycle",
                                                           "Hedgehog signaling", "TCA Cycle", "Adherens Junction",
                                                           "Pancreatic Cancer", "ABC Transporters")
                                   ),
                                   radioButtons("Tables", "Choosing the table", c("C1-MIA" = "C1-MIA",
                                                                                  "C4-MIA" = "C4-MIA",
                                                                                  "GEM_HuR-GEM_IgG","MIA_HuR-MIA_IgG",
                                                                                  "OLAP_HuR-OLAP_IgG","rnaseq.MIA.C.HvL",
                                                                                  "rnaseq.MIA.K.HvL","rnaseq.MIA.H.CvK",
                                                                                  "rnaseq.MIA.L.CvK")),
                                   radioButtons("Regulation", "Whether up or down", c("UP" = "up",
                                                                                      "DOWN"= "down"))
                                   
                            ), #column 1 ends
                            column(5,
                                   fileInput(
                                     "Upload.file", "Loading the plain gene list",
                                     accept = c("text/csv",
                                                "text/comma-separated-values,text/plain",
                                                ".csv"))
                                   
                                   
                            )# column 2 ends
                          ),
                          fluidRow(
                            offset = 1,
                            downloadButton("Download.tab5.1", "Dowloadin the table of pathway results"),
                            downloadButton("Download.tab5.2", "Downloading the table based on file input")
                            
                          )  
                          
                        ),#wellpanel ends
                        tags$hr(),
                        
                        
                        mainPanel(
                          fluidRow(
                            h3("Searching result"),
                            dataTableOutput("table1")
                          ),
                          fluidRow(
                            h3("The searching result for gene list"),
                            dataTableOutput("gene.text"),
                            verbatimTextOutput("text.message")
                            
                          )
                          
                        )#mainpanel ends
                        
                      )
                               )
                             )# navbarpage ends
)



        
        