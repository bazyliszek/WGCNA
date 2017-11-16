# http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/ModulePreservation/Tutorials/simulation-PPInetworks.pdf
unlink(".RData")
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
my_data_low <- my_data[, c(1, 8:13, 20:25, 32:37, 44:49)] 
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
my_data_high <- my_data[, c(1:7, 14:19, 26:31, 38:43)] 
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
setLabels = c("Low", "High")
# Important: components of multiExpr must carry identificating names
names(multiExpr) <- setLabels
# Display the dimensions of the expression data (if you are confused by this construct, ignore it):
lapply(multiExpr, lapply, dim)  # they  need to be matching

# Loading of module labels ------------------------------------------------
colorLow <- mergedColors_low
### inNetwork <- colorLow != "grey" #change to grey or get new file

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

pdf("03112017_Histogram of kTotal for high and low with cut of zero.pdf")
hist(imc_high$kTotal[vect1])
hist(imc_low$kTotal[vect2])
dev.off()
#

inNetwork <- imc_high$kTotal > cutoff | imc_low$kTotal > cutoff
table(inNetwork)

pdf("03112017_Histogram of kTotal for high or low bigger then zero in one of the network.pdf")
par(mfrow=c(2,1))
hist(imc_high$kTotal[inNetwork], xlim=c(0, 700))
hist(imc_low$kTotal[inNetwork], xlim=c(0, 700))
dev.off()
#

unique(colorLow)
length(colorLow)

colorHigh <- mergedColors_high

length(colorHigh)
unique(colorHigh)

length(colorHigh) == length(colorLow)  # TRUE ok!

#colorList <- list(colorLow, colorHigh)
#names(colorList) <- setLabels


colorList_LowHigh <- list(colorLow, colorHigh)

names(colorList_LowHigh) <- setLabels

## Network 1 is low and Network 2 is high
# Match  identified modules to the reference labels
colorList_refhigh = list()
for (set in 1:nSets) {
  colorList_refhigh[[set]] = matchLabels(colorList_LowHigh[[set]], reference = colorHigh)
}
#

colorList_reflow = list()
for (set in 1:nSets) {
  colorList_reflow[[set]] = matchLabels(colorList_LowHigh[[set]], reference = colorLow)
}
#


# Calculation of module preservation statistics ---------------------------
#my_perm <- 1
my_perm <- 200   #  adjust to 200 in original script and 1 for tests
#mp <- modulePreservation(multiExpr, 
#                         colorList_refhigh, dataIsExpr = T, networkType = 'signed',
#referenceNetworks = c(1:2), loadPermutedStatistics = FALSE, 
#nPermutations = my_perm, verbose = 3,  calculateQvalue = T, calculateClusterCoeff = T, greyName = T, )
#referenceNetwork means that set 1 is a reference first and then set 2 is a reference 

mp <- modulePreservation(multiExpr,
                         colorList_refhigh, # once could test  colorList = c(colorList_refhigh, colorList_reflow)
                         dataIsExpr = T,
                         networkType = 'signed',
                         referenceNetworks = c(1:2),
                         loadPermutedStatistics = FALSE,
                         nPermutations = my_perm,
                         calculateQvalue = T,       # This should be set T
                         calculateClusterCoeff = T,   # This should be set T
                         greyName = T,
                         verbose = 3, maxModuleSize = 5000)


str(mp)
save(mp, file = "2001_Preserved network.RData")
#load("2001_Preserved network.RData")

str(mp)

# Calculation of In-Group proportion --------------------------------------
library(clusterRepro)
#install.packages("impute") # impute is not available for 3.1.0 
#and we do not have missing expression data so we do exclude from further analysis 

#source("http://bioconductor.org/biocLite.R")
#biocLite("impute")
library(impute)

# Impute missing data and calculate eigengenes
impExpr = list()
for (set in 1:nSets)
{
impExpr[[set]] = list(data =t(impute.knn(t(multiExpr[[set]]$data))$data))
}

eigengenes = list()
for (set in 1:nSets)
{
eigengenes[[set]] = multiSetMEs(impExpr, universalColors = colorList_refhigh[[set]], excludeGrey = TRUE) # #changed 2001
for (ss in 1:nSets)
{
rownames(eigengenes[[set]][[ss]]$data) = rownames(multiExpr[[ss]]$data)
}}

checkSets(impExpr, checkStructure = TRUE)

# Here comes the IGP calculation
library(clusterRepro)
cr = list()
set.seed(20)
for (ref in 1:nSets)
{
cr[[ref]] = list()
for (test in 1:nSets)
{
printFlush(system.time({
cr[[ref]][[test]] = clusterRepro(Centroids = as.matrix(eigengenes[[ref]][[test]]$data),
New.data = as.matrix(impExpr[[test]]$data),
Number.of.permutations = 10000)}))  # adjust 100 to 10000
collectGarbage()
}
}

# Save the results
save(cr, file = "2001_Low-High-clusterRepro_5000_2.RData")

# Stats -------------------------------------------------------------------
# Load the module preservation statistics
# load(file = "Low-High-clusterRepro_10000.RData")
# load(file = "Preserved network.RData")

ref <-  2 # Select the High data as reference   # changed 20012016
test <- 1 #  Select the Low data as test
# We could switch ref and test if we want 
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])

print(signif(statsZ[, "Zsummary.pres", drop = FALSE],2))
# Compare preservation to quality:
print(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

# Plots -------------------------------------------------------------------
ref <-  2   #changed 20012016
test <- 1
ind <-  1  # what is this for?
stats <-  mp$preservation$observed[[ref]][[test]]
labelsX <- rownames(stats)
labelsX[labelsX=="gold"] <- "grey99"  # if gold is there it will be replaces by grey99
modColors <- labelsX
plotMods <- !(modColors %in% c("grey", "grey99"))  
moduleSizes <- stats[plotMods, 1]
textLabels <- match(modColors, standardColors(length(labelsX)))[plotMods]  # DB: cols no from labelsX
colorLabels <- labelsX[plotMods]

# Further graphs ----------------------------------------------------------
nModules = sum(plotMods)
nPlots = 6
plotData = list()

# Fill up the plotData
plotData[[1]] = plotData[[2]] = matrix(0, nModules, nPlots)
plotData[[1]][, c(1:4)] = moduleSizes
plotData[[2]][, 1] = mp$preservation$Z[[ref]][[test]]$Zsummary.pres[plotMods]
plotData[[2]][, 2] = mp$preservation$observed[[ref]][[test]]$medianRank.pres[plotMods]

# Match the modulePreservation ordering of modules to that of clusterRepro
crLabels = sort(unique(colorLabels))
mp2cr = match(colorLabels, crLabels)

# Scatterplots of IGP and p-value vs. module size
plotData[[2]][, 3] = cr[[ref]][[test]]$Actual.IGP   ### CHECK - DID NOT WORK     DB: DID work for me
plotData[[2]][, 4] = -log10(cr[[ref]][[test]]$p.value + 1e-4)    ### CHECK - DID NOT WORK      DB: DID work for me

# Scatterplot of observed IGP vs. Zsummary and medianRank
plotData[[1]][, c(5,6)] = plotData[[2]][, c(1:2)]
plotData[[2]][, c(5,6)] = plotData[[2]][, 3]

# Plot annotation
xLabs = c(rep("Module size", 4), "Zsummary", "Median rank")
yLabs = c("Zsummary", "Median rank", "Observed IGP", "-log10(IGP perm p)", "Observed IGP", "Observed IGP")
mains = spaste(LETTERS[1:nPlots], ". ", #rep("Ref: Human, Test: Chimp\n", nPlots),
c(yLabs[1:4], paste(yLabs[5:6], "vs.", xLabs[5:6])),
c("", "", "", "", "\n", "\n"))

# Scatterplot options
verbose = c(rep(FALSE, 4), rep(TRUE, 2))
ablines = list(c(0, 2, 10), NA, NA, c(-log10(0.05), -log10(0.05/nModules)), NA, NA)
abColors = list(c("black", "blue", "darkgreen"), NA, NA, c("blue", "red"), NA, NA)
logs = c("x", "x", "x", "x", "", "")
invertY = c(FALSE, TRUE, rep(FALSE, 4))
verSP = function(...) { verboseScatterplot(..., abline = TRUE) }

cexLabels = 1.4
#sizeGrWindow(6,9)
pdf(file = "03112017_Low_Specific-NetworkAndIGPStatistics.pdf", w=6, h=9, onefile = FALSE)
par(mfrow = c(3,2))
par(mar = c(3.3, 3.3, 3.2, 0.5))
par(mgp = c(2, 0.6, 0))
for (p in 1:nPlots)
{
x = plotData[[1]][, p]
y = plotData[[2]][, p]
miny = min(y, ablines[[p]], na.rm = TRUE)
maxy = max(y, ablines[[p]], na.rm = TRUE)
miny = miny - (maxy-miny)*0.1
maxy = maxy + (maxy-miny)*0.1
(if (verbose[p]) verSP else plot ) (plotData[[1]][, p], plotData[[2]][, p],
main = mains[p],
xlab = xLabs[p],
ylab = yLabs[p],
cex.main = cexLabels, cex.lab = cexLabels, cex.axis = cexLabels,
bg = colorLabels,
col = colorLabels, cex = 2.2,
ylim = if (invertY[p]) c(maxy, miny) else c(miny, maxy),
pch = 21,
log = logs[p])
labelPoints(plotData[[1]][, p], plotData[[2]][, p], textLabels, cex = cexLabels, offs = 0.06)
if (!is.na(ablines[[p]][[1]]))
for (al in 1:length(ablines[[p]]))
abline(h = ablines[[p]][[al]], col = abColors[[p]][[al]], lty = 2)
}
dev.off()

# Plots of all statistics -------------------------------------------------
figNames = c("quality", "preservation")
useStats = list(c(1:8)+1, c(9:24)+1) # There are 24 statistics, CHECK in STATS!
sectioning = list(c(2,4), c(4,5))
dims = list(sectioning[[1]] * 1.8, sectioning[[2]]*1.8)

for (f in 1:2)  {
pdf(file=spaste("03112017_LowHigh-LowSpecific-modulePreservation-", figNames[f],"-%02d.pdf"),

w=dims[[f]][2], h=dims[[f]][1], onefile = FALSE)
for (ref in 1:nSets) for (test in 1:nSets) if (ref!=test)
{
stats = cbind(mp$quality$Z[[ref]][[test]],   # 19 25
mp$referenceSeparability$Z[[ref]][[test]][, -1, drop = FALSE],  # 19 1
mp$preservation$Z[[ref]][[test]][, -1],  # 19 13
mp$accuracy$Z[[ref]][[test]][, -1], # 19 3
mp$testSeparability$Z[[ref]][[test]][, -1, drop = FALSE]) # 19 1
labelsX = rownames(stats)
labelsX[labelsX=="gold"] = "grey99"
modColors = labelsX
moduleSizes = stats[, 1]
plotMods = !(modColors %in% c("grey", "grey99")) # 17
plotStats = useStats[[f]]
textLabels = match(modColors, standardColors(length(labelsX)))[plotMods]  # 17
par(mfrow = sectioning[[f]])
par(mar = c(3.0, 3.0, 4, 0.4))
par(mgp = c(1.7, 0.6, 0))

  for (s in plotStats)
{
min = min(stats[plotMods, s], na.rm = TRUE)
max = max(stats[plotMods, s], na.rm = TRUE)
minMS = min(moduleSizes[plotMods])
maxMS = max(moduleSizes[plotMods])
nms = colnames(stats)
if (max < 10) max = 10
if (min > -max/5) min = -max/5

plot(moduleSizes[plotMods], stats[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
main = paste("Ref.:", setLabels[ref], "\nTest:", setLabels[test], "\n", nms[s]),
cex = 1.2,
cex.main = 1.0,
ylab = nms[s], xlab = "Module size", log = "x",
#ylim = c(min, max + 0.1 * (max-min)),
#ylim = c(5, 1400 + 0.1 * (max-min)),
ylim = range(c(min, max + 0.1 * (max-min)), na.rm=TRUE),
xlim = c(minMS/(maxMS/minMS)^0.1, maxMS*(maxMS/minMS)^0.1)
)
labelPoints(moduleSizes[plotMods], stats[plotMods, s], textLabels, cex = 0.90, offs = 0.06)
abline(h=0)
abline(h=2, col = "blue", lty = 2)
abline(h=10, col = "darkgreen", lty = 2)
}
}
dev.off()
}

# Output of complete results  ---------------------------------------------
# This variable will contain the summary table
summaryTable = NULL
# Loop over all combinations of reference and tests sets
for (ref in 1:nSets) for (test in 1:nSets) if (ref!=test)
{
modules = rownames(mp$preservation$Z[[ref]][[test]]);
nMods = length(modules);
sizes = mp$preservation$Z[[ref]][[test]][, 1];
acc = matrix(NA, nMods, 3);
if (test!=4)
{
acc[match(rownames(mp$accuracy$observed[[ref]][[test]]), modules), ] =
mp$accuracy$observed[[ref]][[test]][, -1, drop = FALSE];
colnames(acc) = colnames(mp$accuracy$observed[[ref]][[test]])[-1];
accZ = mp$accuracy$Z[[ref]][[test]][, -1, drop = FALSE];
acc.log.p = mp$accuracy$log.p[[ref]][[test]][, -1, drop = FALSE];
acc.log.pBonf = mp$accuracy$log.pBonf[[ref]][[test]][, -1, drop = FALSE];
} else {
accZ = matrix(NA, nMods, 3);
acc.log.p = matrix(NA, nMods, 3);
acc.log.pBonf = matrix(NA, nMods, 3);
colnames(acc) = colnames(mp$accuracy$observed[[1]][[2]])[-1];
colnames(accZ) = colnames(mp$accuracy$Z[[1]][[2]])[-1];
colnames(acc.log.p) = colnames(mp$accuracy$log.p[[1]][[2]])[-1];
colnames(acc.log.pBonf) = colnames(mp$accuracy$log.pBonf[[1]][[2]])[-1];
}
# Table of results for this reference-test combination
tab = cbind(referenceSet = rep(setLabels[ref], nMods),
testSet = rep(setLabels[test], nMods),
moduleLabel = modules,
moduleSize = sizes,
mp$quality$observed[[ref]][[test]][, -1, drop = FALSE],
mp$preservation$observed[[ref]][[test]][, -1, drop = FALSE],
acc,
mp$referenceSeparability$observed[[ref]][[test]][, -1, drop = FALSE],
mp$testSeparability$observed[[ref]][[test]][, -1, drop = FALSE],
mp$quality$Z[[ref]][[test]][, -1, drop = FALSE],
mp$quality$log.p[[ref]][[test]][, -1, drop = FALSE],
mp$quality$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
mp$preservation$Z[[ref]][[test]][, -1, drop = FALSE],
mp$preservation$log.p[[ref]][[test]][, -1, drop = FALSE],
mp$preservation$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
accZ,
acc.log.p,
acc.log.pBonf,
mp$referenceSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
mp$referenceSeparability$log.p[[ref]][[test]][, -1, drop = FALSE],
mp$referenceSeparability$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
mp$testSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
mp$testSeparability$log.p[[ref]][[test]][, -1, drop = FALSE],
mp$testSeparability$log.pBonf[[ref]][[test]][, -1, drop = FALSE]
)
# Add the table to the main table.
if (is.null(summaryTable)) summaryTable = tab else summaryTable = rbind(summaryTable, tab);
}
# Save the table in csv format.
write.table(summaryTable, file = "LowHigh-lowSpecific-completeResults.csv", row.names = FALSE,
sep = ",", quote = FALSE);

#install.packages("flashClust")
library(flashClust)
nSets
#Marcin´s dendograms
dendrograms = list(geneTree_low, geneTree_high)

# Calculate the contingency table and p-values
overlap = overlapTable(unlist(colorList_refhigh[2])[inNetwork], unlist(colorList_refhigh[1])[inNetwork]) #changed 2001
# The numMat will encode color. We use -log of the p value.
numMat = -log10(overlap$pTable)
numMat[numMat > signif(max(numMat), 2)] = signif(max(numMat), 2) # changed from 50!
# Prepare for generating a color-coded plot of the overlap table. The text of the table will consist of
# counts and corresponding p-values.
textMat = paste(overlap$countTable, "\n", signif(overlap$pTable, 2))
dim(textMat) = dim(numMat)
# Additional information for the plot. These will be used shortly.
xLabels = paste("M", sort(unique(unlist(colorList_refhigh[1]))))
yLabels = paste("M", sort(unique(unlist(colorList_refhigh[2]))))
xSymbols = paste(sort(unique(unlist(colorList_refhigh[1]))), ": ", table(unlist(colorList_refhigh[1][inNetwork])), sep = "")
ySymbols = paste(sort(unique(unlist(colorList_refhigh[2]))), ": ", table(unlist(colorList_refhigh[2][inNetwork])), sep = "")

### Figure 3 publication
pdf("03112017_MotivationFigure-dendrosAndTable.pdf")
fp = TRUE
layout(matrix(c(1,2,5, 3,4,5), 3, 2),
       heights = c(3, 1, 5.5), widths = c(1, 1));
par(mgp = c(3, 1, 0))
plotDendroAndColors(dendrograms[[2]],
                    cbind(unlist(colorList_refhigh[2])[inNetwork], unlist(colorList_refhigh[1])[inNetwork]),
                    c("High P modules", "Low P modules"),
                    setLayout = FALSE,
                    marAll = c(1, 6, 2.7, 0.2),
                    addGuide = FALSE,
                    main = "A. High P dendrogram\nand module colors", cex.main = 1.2,
                    dendroLabels = FALSE, hang = 0.03, cex.colorLabels = 0.7, abHeight = 0.8, 
                    ylim = c(0.5,1))

par(mgp = c(3, 1, 0))
plotDendroAndColors(dendrograms[[1]],
                    cbind(unlist(colorList_reflow[1])[inNetwork], unlist(colorList_reflow[2])[inNetwork]), # [1] is low and [2] is high
                    c("Low P modules", "High P modules"),
                    setLayout = FALSE,
                    marAll = c(1, 6, 2.7, 0.2),
                    addGuide = FALSE,
                    main = "B. Low P dendrogram\nand module colors", cex.main = 1.2,
                    dendroLabels = FALSE, hang = 0.03, cex.colorLabels = 0.7, abHeight = 0.8, 
                    ylim = c(0.5,1))



# Plot the overlap table
fcex = 0.4
pcex = 0.4
fcexl = 1.00
pcexl = 1.00
par(mar = c(6, 7, 2, 1.0))
labeledHeatmap(Matrix = numMat,
               xLabels = xLabels, xSymbols = xSymbols,
               yLabels = yLabels, ySymbols = ySymbols,
               colorLabels = TRUE,
               colors = greenWhiteRed(100)[50:100],
               textMatrix = textMat, cex.text = if (fp) fcex else pcex, setStdMargins = FALSE,
               cex.lab = if (fp) fcexl else pcexl,
               xColorWidth = 0.08,
               main = "C. High modules (rows) vs. Low modules (columns)", cex.main = 1.2)
dev.off()

