# Hello
#replace here ancient withhigh
#replace here recent withlow

setwd("~/Downloads/functionalenrichmentkog_Marcin") # Here is your R script and analysis
path="./Module_gene_list_txts" # This is a folder name 

setwd("~/Downloads/preserved_P")
dir_original <- getwd()

# Load everything needed
library(qvalue)
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# objects names are the same and we will overwrite here! 
load(file="GeneTree_and_dynamicColors_USA_high.RData")
geneTree_high <- geneTree_USA
dynamicColors_high <- dynamicColors_USA 

load(file="GeneTree_and_dynamicColors_USA_low.RData")
geneTree_low <- geneTree_USA
dynamicColors_low <- dynamicColors_USA 


load(file="mergedColors_low.RData")
str(mergedColors_low)
load(file="mergedColors_high.RData")
str(mergedColors_high)

#import data low
# Part 1 Importing library, data we need setting up original directory-------------------------------
dir_original <- getwd()  # remember to bring it back after you introduce new directory

# Importing data, text = average intensities of DE genes - from contrasts (any), fdr less then 5 %  -------------------------------
# 8 replicates each, for low data
my_data <- read.table("normalized_intens_all_genes.txt", header=T)
my_data_low <- my_data[, c(1:13, 38:49)] 

dim(my_data_low) # 63265 25

# Add vector for adjusting the gene set we will use for analysis, e.g. DE genes only 
DE_genes_low <- read.table("AllContrasts_DEgenelists_q005_all_in_one.txt")
dim(DE_genes_low)  # 38316 1
DE_genes_low_unique <- unique(DE_genes_low)
dim(DE_genes_low_unique) # 12119
DE_genes_low_unique_vector <- DE_genes_low_unique[1:dim(DE_genes_low_unique)[1],] # 12119 DE genes - combining all 3 contrasts
length(DE_genes_low_unique_vector) # 12119
head(DE_genes_low_unique_vector, 10)

# Match genes of interest with whole data set
matched_geneIDs_low <- my_data_low[match(DE_genes_low_unique_vector, 
                                               my_data_low$gene), ]
dim(matched_geneIDs_low)  # 12119    25

low_data <- matched_geneIDs_low
dim(low_data) # 12119    25

# Remove subset [1:100] for computing  -------------------------------
#my_data_low <- low_data[1:500,] # adjust
my_data_low <- low_data

# Make-up -------------------------------
my_low_data <- my_data_low[,c(-1)] 
rownames(my_low_data) <- my_data_low[,1]
str(my_low_data) # must be data frame

# Transpose data for correlations -------------------------------
datExpr1_low <- as.data.frame(t(my_low_data))


# Import data high -------------------------------------------------------
# Part 1 Importing library, data we need setting up original directory-------------------------------
my_data <- read.table("normalized_intens_all_genes.txt", header=T)
my_data_high <- my_data[, c(1, 14:37)]
dim(my_data_high) # 63265 25

# Add vector for adjusting the gene set we will use for analysis, e.g. DE genes only 
DE_genes_high <- read.table("AllContrasts_DEgenelists_q005_all_in_one.txt")
dim(DE_genes_high)  # 38316 1
DE_genes_high_unique <- unique(DE_genes_high)
dim(DE_genes_high_unique) # 12119
DE_genes_high_unique_vector <- DE_genes_high_unique[1:dim(DE_genes_high_unique)[1],] # 12119 DE genes - combining all 3 contrasts
length(DE_genes_high_unique_vector) # 12119
head(DE_genes_high_unique_vector, 10)

# Match genes of interest with whole data set
matched_geneIDs_high <- my_data_high[match(DE_genes_high_unique_vector, 
                                                 my_data_high$gene), ]
dim(matched_geneIDs_high)  # 12119    25

high_data <- matched_geneIDs_high
dim(high_data) # 12119    25

# Remove subset [1:100] for computing  -------------------------------
#my_data_high <- high_data[1:500,] #adjust
my_data_high <- high_data

# Make-up -------------------------------
my_high_data <- my_data_high[,c(-1)] 
rownames(my_high_data) <- my_data_high[,1]
str(my_high_data) # must be data frame

# Transpose data for correlations -------------------------------
datExpr1_high <- as.data.frame(t(my_high_data))

# # Number of data sets that we work with ---------------------------------
nSets = 2  # this is 2 networks! 
# Object that will contain the expression data
multiExpr = list()
multiExpr[[1]] = list(data = datExpr1_low)   # low datExpr
multiExpr[[2]] = list(data = datExpr1_high)   # high datExpr

# Names for the two sets
setLabels = c("low", "high")
# Important: components of multiExpr must carry identificating names
names(multiExpr) <- setLabels
# Display the dimensions of the expression data (if you are confused by this construct, ignore it):
lapply(multiExpr, lapply, dim)  # they  need to be matching

# Loading of module labels ------------------------------------------------
colorlow <- mergedColors_low
### inNetwork <- colorlow != "grey" #change to grey or get new file

### Creating inNetwork logical vector
imc_high <- read.table("Intramodular connectivity_high.txt", header=T)
imc_low <- read.table("Intramodular connectivity_low.txt", header=T)

dim(imc_high) == dim(imc_low)
dim(imc_high) 

cutoff <- 0  # can be adjusted for highly connected genes

vect1 <- imc_high$kTotal > cutoff
table(vect1 <- imc_high$kTotal > cutoff)

vect2 <- imc_low$kTotal > cutoff
table(vect2 <- imc_low$kTotal > cutoff)

pdf("1202_2017_Histogram of kTotal for high and low with cut of zero.pdf")
hist(imc_high$kTotal[vect1])
hist(imc_low$kTotal[vect2])
dev.off()
#

inNetwork <- imc_high$kTotal > cutoff | imc_low$kTotal > cutoff
table(inNetwork)

pdf("1202_2017_Histogram of kTotal for high or low bigger then zero in one of the network.pdf")
par(mfrow=c(2,1))
hist(imc_high$kTotal[inNetwork], xlim=c(0, 700))
hist(imc_low$kTotal[inNetwork], xlim=c(0, 700))
dev.off()
#

unique(colorlow)
length(colorlow)

colorhigh <- mergedColors_high

length(colorhigh)
unique(colorhigh)

length(colorhigh) == length(colorlow)  # TRUE ok!

#colorList <- list(colorlow, colorhigh)
#names(colorList) <- setLabels


colorList_lowhigh <- list(colorlow, colorhigh)

names(colorList_lowhigh) <- setLabels

## Network 1 is low and Network 2 is high
# Match  identified modules to the reference labels
colorList_refhigh = list()
for (set in 1:nSets) {
  colorList_refhigh[[set]] = matchLabels(colorList_lowhigh[[set]], reference = colorhigh)
}
#

colorList_reflow = list()
for (set in 1:nSets) {
  colorList_reflow[[set]] = matchLabels(colorList_lowhigh[[set]], reference = colorlow)
}
#

# Load other stuff that were clalculated for 1 week! ----------------------
load(file="2001_Preserved_network.RData")
load(file = "2001_Low-High-clusterRepro_5000_2.RData")

# Match  identified modules to the reference labels

# extra nessesary stuff ---------------------------------------------------
dendrograms = list(geneTree_low, geneTree_high)

# Calculate the contingency table and p-values
overlap = overlapTable(unlist(colorList_refhigh[2])[inNetwork], unlist(colorList_refhigh[1])[inNetwork]) #changed 2001
# These both unlisted lists containing 10439 genes 
# The numMat will encode color. We use -log of the p value.
numMat = -log10(overlap$pTable)
numMat[numMat >50] = 50


# Additional information for the plot. These will be used shortly.
xLabels = paste("M", sort(unique(unlist(colorList_refhigh[1]))))
yLabels = paste("M", sort(unique(unlist(colorList_refhigh[2]))))
xSymbols = paste(sort(unique(unlist(colorList_refhigh[1]))), ": ", table(unlist(colorList_refhigh[1][inNetwork])), sep = "")
ySymbols = paste(sort(unique(unlist(colorList_refhigh[2]))), ": ", table(unlist(colorList_refhigh[2][inNetwork])), sep = "")




# Also genes loaded  ------------------------------------------------------
my_data <- read.table("normalized_intens_all_genes.txt", header=T)
dim(my_data) # 63265 49
#Add vector for adjusting the gene set we will use for analysis, e.g. DE genes only 
DE_genes_USA <- read.table("AllContrasts_DEgenelists_q005_all_in_one.txt")
dim(DE_genes_USA)  # 38316 1
DE_genes_USA_unique <- unique(DE_genes_USA)
dim(DE_genes_USA_unique) # 12119
DE_genes_USA_unique_vector <- DE_genes_USA_unique[1:10439,] # 12119 DE genes - combining all 3 contrasts

length(DE_genes_USA_unique_vector) # 12119
head(DE_genes_USA_unique_vector, 10)

# Match genes of interest with whole data set
matched_geneIDs <- my_data[match(DE_genes_USA_unique_vector, my_data$gene), ]
dim(matched_geneIDs)  # 12119  49

USA_data <- matched_geneIDs
head(USA_data)

Genes_all_name=matched_geneIDs[,1]
Genes = substring(Genes_all_name, 6, 20)

#These condition must be true
length(Genes) == length(unlist(colorList_refhigh[2])[inNetwork])


# Making a table with Genes ID and colours --------------------------------
mytable_left = data.frame(Genes, unlist(colorList_refhigh[2])[inNetwork])
colnames(mytable_left) <- c("Genes", "Left")

mytable_right = data.frame(Genes, unlist(colorList_refhigh[1])[inNetwork])
colnames(mytable_right) <- c("Genes", "Right")


# Colours needs to go alphabetically and I will use new vectors to assign these --------------------------------------
leftcol = sort(unique(unlist(colorList_refhigh[2])[inNetwork]))
rightcol = sort(unique(unlist(colorList_refhigh[1])[inNetwork]))


mycollist_left <- list()
for (i in 1:length(leftcol)) {
  print(paste("My left color", i, "is now taken"))
  print(paste("which is", leftcol[i]))
  mycollist_left[[i]] <- mytable_left[mytable_left$Left == leftcol[i], 1]
}
names(mycollist_left) <- paste("left",leftcol)


mycollist_right <- list()
for (i in 1:length(rightcol)) {
  print(paste("My right color", i, "is now taken"))
  print(paste("which is", rightcol[i]))
  mycollist_right[[i]] <- mytable_right[mytable_right$Right == rightcol[i], 1]
}
names(mycollist_right) <- paste("right",rightcol)


# Generating the final list  ----------------------------------------------
newlists_left_colors <- list()
for (i in 1:length(mycollist_left)) {
  print(paste("My left color", i, "is now taken which is", leftcol[i]))
  for (j in 1:length(mycollist_right)){
    print(paste("With right color", j, "which is", rightcol[j]))
    name <- paste("my_color_left_",leftcol[i],"_and_my_color_right_",rightcol[j], sep="")
    tmp <- intersect(mycollist_left[[i]], mycollist_right[[j]])  #removed list
    newlists_left_colors[[name]] <- tmp
  }
}
str(newlists_left_colors)

# Print all the lists into the files for further KOG annotation--------------------------------------
d1_preserved <- "DE_preserved"
if (!file.exists(d1_preserved)) dir.create(d1_preserved)
dir_temp_preserved <- paste0("./", d1_preserved, "")
setwd(dir_temp_preserved)

for (z in 1:length(newlists_left_colors)) {
  pre <- "Preserved_"
  myname <- names(newlists_left_colors[z])
  post <- ".txt"
  write.table(newlists_left_colors[[z]], file=paste0(pre, myname, post, sep=""), row.names = F, col.names = F) # not required as separate doc
}

setwd(dir_original)
