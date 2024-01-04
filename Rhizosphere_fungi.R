devtools::install_github("GraceYoon/SPRING",force = T)
devtools::install_github("GraceYoon/SPRING")

devtools::install_github("stefpeschel/NetCoMi", 
                         dependencies = c("Depends", "Imports", "LinkingTo"),
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))

# Loading required packages
library(BiocManager)
library(WGCNA)
library(NetCoMi); packageVersion("NetCoMi")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2") # graphics
library(readxl)       # necessary to import the data from Excel file
library(dplyr)        # filter and reformat data frames
library(tibble)       # Needed for converting column to row names
library(qiime2R)
library(tidyverse)
library(microbiome)

setwd("C:/Users/User/OneDrive/Desktop/Taiwan/Rhizosphere")

otu_mat<- read_excel("Rhizosphere_fungi.xlsx", sheet = "OTU matrix")
tax_mat<- read_excel("Rhizosphere_fungi.xlsx", sheet = "Taxonomy table")
samples_df <- read_excel("Rhizosphere_fungi.xlsx", sheet = "Samples")

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu") 

tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")

samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

physeq201 <- phyloseq(OTU, TAX, samples)
physeq201

#aggregating taxa in the Genus level
complete_agg <- tax_glom(physeq201, 'Genus')
#complete_agg
tax_table(complete_agg)
complete_agg <- aggregate_taxa(complete_agg, 'Genus')
complete_agg

# keep only taxa that were ....................................................
minTotRelAbun <- 5e-5
sums <- taxa_sums(complete_agg)
#sums
keepTaxa <- taxa_names(complete_agg)[which((sums / sum(sums)) > minTotRelAbun)]
#keepTaxa
filt_complete_agg <- prune_taxa(keepTaxa, complete_agg)
filt_complete_agg



# select data of condition.
Healthy_physeq <- subset_samples(filt_complete_agg, Condition =="Healthy")
Disease_physeq <- subset_samples(filt_complete_agg, Condition =="Diseased")

# Saving the phyloseq objects for updated plant condition
saveRDS(Healthy_physeq, "condition/Healthy.rds")
saveRDS(Disease_physeq, "condition/Disease.rds")

#Sparcc#######################################################################################################

#Healthy********************************************************************************************************
#15 Samples
Healthy_data <- readRDS("condition/Healthy.rds")

top_Healthy <- prune_taxa(names(sort(taxa_sums(Healthy_data),TRUE)[1:50]), Healthy_data)
plot_heatmap(top_Healthy)
Healthy_data

# data(Healthy_data)

data <-
  Healthy_data %>%
  tax_glom("Genus") %>%
  transform_sample_counts(function(x)100* x / sum(x)) %>%
  psmelt() %>%
  as_tibble()

# highest abundance: all samples pooled together
data %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)

#SparCC
net_single_Healthy <- netConstruct(Healthy_data,
                                   verbose = 3,
                                   filtTax = "highestFreq",
                                   filtTaxPar = list(highestFreq = 100),
                                   filtSamp = "totalReads",
                                   filtSampPar = list(totalReads = 1000),
                                   zeroMethod = "none", normMethod = "none",
                                   measure = "sparcc",
                                   sparsMethod = "threshold", thresh = 0.4,
                                   adjust = "adaptBH",
                                   dissFunc = "signed", 
                                   seed = 123456)

saveRDS(net_single_Healthy, "condition/Healthy_condition/Healthy_Sparcc network.rds")

props_single_Healthy <- netAnalyze(net_single_Healthy, 
                                   clustMethod = "cluster_fast_greedy",
                                   hubPar = "eigenvector", 
                                   hubQuant = 0.95)

saveRDS(props_single_Healthy, "condition/Healthy_condition/Healthy_Sparcc analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Healthy)
net.summary

p_Healthy <- plot(props_single_Healthy,
                  shortenLabels = "none",
                  labelScale = FALSE,
                  rmSingles = "all",
                  nodeSize = "eigenvector",
                  nodeColor = "cluster",
                  hubBorderCol = "blue",
                  cexNodes = 1,
                  cexLabels = 0.5,
                  edgeWidth = 1,
                  highlightHubs = TRUE,
                  cexHubs = 1.5,
                  title1 = "Healthy condition Network on Genus with SparCC Method-Rhizosphere", 
                  showTitle = TRUE,
                  cexTitle = 1.5)

saveRDS(p_Healthy, "condition/Healthy_condition/Healthy_Sparcc plot.rds")


from <- p_Healthy$labels$labels1[p_Healthy$q1$Edgelist$from]
to <- p_Healthy$labels$labels1[p_Healthy$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Healthy$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Healthy$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Healthy$q1$Edgelist$from)))
edges$from <- p_Healthy$q1$Edgelist$from
edges$to <- p_Healthy$q1$Edgelist$to
edges$weight <- p_Healthy$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "condition/Healthy_condition/Healthy edge data.csv", row.names = FALSE)

hubs <- props_single_Healthy$hubs$hubs1
write(hubs, "condition/Healthy_condition/Healthy condition Hubs.txt")

node.lables <- p_Healthy$labels$labels1
clust <- props_single_Healthy$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Healthy$centralities$degree1[nodes$lable]
evs=props_single_Healthy$centralities$eigenv1[nodes$lable]
betweennesses=props_single_Healthy$centralities$between1[nodes$lable]
closenesses=props_single_Healthy$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "condition/Healthy_condition/Healthy node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "condition/Healthy_condition/Healthy hub data.csv", row.names = FALSE)

#Disease********************************************************************************************************
#15 Samples
Disease_data <- readRDS("condition/Disease.rds")

top_Disease <- prune_taxa(names(sort(taxa_sums(Disease_data),TRUE)[1:50]), Disease_data)
plot_heatmap(top_Disease)
Disease_data

data(Disease_data)

data <-
  Disease_data %>%
  tax_glom("Genus") %>%
  transform_sample_counts(function(x)100* x / sum(x)) %>%
  psmelt() %>%
  as_tibble()

# highest abundance: all samples pooled together
data %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)

#SparCC
net_single_Disease <- netConstruct(Disease_data,
                                   verbose = 3,
                                   filtTax = "highestFreq",
                                   filtTaxPar = list(highestFreq = 100),
                                   filtSamp = "totalReads",
                                   filtSampPar = list(totalReads = 1000),
                                   zeroMethod = "none", normMethod = "none",
                                   measure = "sparcc",
                                   sparsMethod = "threshold", thresh = 0.4,
                                   adjust = "adaptBH",
                                   dissFunc = "signed", 
                                   seed = 123456)

saveRDS(net_single_Disease, "condition/Disease_condition/Disease_Sparcc network.rds")

props_single_Disease <- netAnalyze(net_single_Disease, 
                                   clustMethod = "cluster_fast_greedy",
                                   hubPar = "eigenvector", 
                                   hubQuant = 0.95)

saveRDS(props_single_Disease, "condition/Disease_condition/Disease_Sparcc analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Disease)
net.summary

p_Disease <- plot(props_single_Disease,
                  shortenLabels = "none",
                  labelScale = FALSE,
                  rmSingles = "all",
                  nodeSize = "eigenvector",
                  nodeColor = "cluster",
                  hubBorderCol = "blue",
                  cexNodes = 1,
                  cexLabels = 0.5,
                  edgeWidth = 1,
                  highlightHubs = TRUE,
                  cexHubs = 1.5,
                  title1 = "Disease condition Network on Genus with SparCC Method-Rhizosphere", 
                  showTitle = TRUE,
                  cexTitle = 1.5)

saveRDS(p_Disease, "condition/Disease_condition/Disease_Sparcc plot.rds")


from <- p_Disease$labels$labels1[p_Disease$q1$Edgelist$from]
to <- p_Disease$labels$labels1[p_Disease$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Disease$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Disease$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Disease$q1$Edgelist$from)))
edges$from <- p_Disease$q1$Edgelist$from
edges$to <- p_Disease$q1$Edgelist$to
edges$weight <- p_Disease$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "condition/Disease_condition/Disease edge data.csv", row.names = FALSE)

hubs <- props_single_Disease$hubs$hubs1
write(hubs, "condition/Disease_condition/Disease condition Hubs.txt")

node.lables <- p_Disease$labels$labels1
clust <- props_single_Disease$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Disease$centralities$degree1[nodes$lable]
evs=props_single_Disease$centralities$eigenv1[nodes$lable]
betweennesses=props_single_Disease$centralities$between1[nodes$lable]
closenesses=props_single_Disease$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "condition/Disease_condition/Disease node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "condition/Disease_condition/Disease hub data.csv", row.names = FALSE)

# Split the phyloseq object into two groups
Healthy_physeq <- subset_samples(filt_complete_agg, Condition =="Healthy")
Disease_physeq <- subset_samples(filt_complete_agg, Condition =="Diseased")
Healthy_physeq
Disease_physeq

# Network construction
net_condition <- netConstruct(data = Healthy_physeq ,
                              data2 = Disease_physeq,
                              verbose = 3,
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              zeroMethod = "none", normMethod = "none",
                              measure = "sparcc",
                              sparsMethod = "threshold", thresh = 0.4,
                              adjust = "adaptBH",
                              dissFunc = "signed", 
                              seed = 123456)


props_condition <- netAnalyze(net_condition, 
                              clustMethod = "cluster_fast_greedy",
                              hubPar = "eigenvector", 
                              hubQuant = 0.95)
summary(props_condition,groupNames = c("Healthy", "Diseased"))

# Network plot
plot(props_condition, sameLayout = TRUE, layoutGroup = "union", 
     nodeSize = "clr", repulsion = 0.9, cexTitle = 3.7, cexNodes = 2,
     cexLabels = 2, groupNames = c("Healthy", "Disease"))



plot(props_condition,
     shortenLabels = "none",
     labelScale = FALSE,
     rmSingles = "inboth",
     nodeSize = "eigenvector",
     nodeColor = "cluster",
     hubBorderCol = "blue",
     sameLayout = TRUE, 
     layoutGroup = "union",
     cexNodes = 1,
     cexLabels = 1.0,
     edgeWidth = 1,
     highlightHubs = TRUE,
     cexHubs = 1.5,
     groupNames = c("Healthy", "Disease"), 
     showTitle = TRUE,
     cexTitle = 1.5)

# Network comparison
comp_condition_orig <- netCompare(props_condition, permTest = TRUE, nPerm = 1000, 
                                  storeAssoPerm = TRUE,
                                  fileStoreAssoPerm = "assoPerm_comp",
                                  storeCountsPerm = FALSE,
                                  adjust = "adaptBH",
                                  cores = 4,
                                  seed = 123456)
summary(comp_condition_orig,groupNames = c("Healthy", "Diseased"))

