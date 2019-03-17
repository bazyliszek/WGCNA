# 2103 MWDB This script makes analysis on all files 
# Load library  ------------------------------------------
library(qvalue)
library(plyr)
library(ggplot2) 
library(plyr)

# KOG analysis for subcategories ------------------------------------------
setwd("~/Downloads/functionalenrichmentkog_Marcin_P16042017") # Here is your R script and analysis
path="./DE_preserved_P" # This is a folder name 

dir_original <- getwd()

# check this server for updated version of  ftp://ftp.ncbi.nih.gov/pub/COG/
# Structure of KOG is following: 
# kog_id (KOG0001) ,
# kog_sub_cat (1-to-5 letters), and 
# KOG_cat (e.g. INFORMATION_STORAGE_AND_PROCESSING)

# check for category of "D" assignment between database and jgi

# This file has been taken from KOG database genome.jgi.doe.gov and it does not need to be updated
KOG <-read.table("KOGcatsubcatdesc.txt", sep="\t", header = T)
head(KOG)
dim(KOG)

# This file has been downlaoded by DB from EMBL
kog <-read.table("kogsubcatid.txt", sep="", header = T)
head(kog)
dim(kog)

# Reduced set of genes that are on array taking into account the problmes with probes
jgi_sub <-read.csv("Dappu1_JGI_Paralog_Tandem_DB_adjusted.csv", sep=",", header = T) ## Notice I removed GENE: for matching 
head(jgi_sub)
dim(jgi_sub)[1]  # 24679

# ddply(jgi_sub, ~ KOG_JGI,
#       summarize, counts = length(KOG_JGI))
# dim(ddply(jgi_sub, ~ KOG_JGI,
#       summarize, counts = length(KOG_JGI)))
# head(ddply(jgi_sub, ~ KOG_JGI,
#           summarize, counts = length(KOG_JGI)))

#replacing the headers
names(kog) <- c("kog_sub_cat", "KOG_JGI", "kog_id_desc")

# DB removed 0s and/or NAs here
jgikog <- na.omit(merge.data.frame(jgi_sub, kog, by="KOG_JGI", all.x=T))  #based on KOG_JGI in jgi and kog_id in kog
dim(jgikog)
head(jgikog) 

jgikog$kog_sub_cat <- as.character(jgikog$kog_sub_cat)

kog_id_count <- as.data.frame(table(jgikog$KOG_JGI))[-1,]
names(kog_id_count) <- c("kog_id", "Counts_array")


# loop for main
KOG_cat <- LETTERS[c(1:23, 25:26)]
KOG_cat_table <- as.data.frame(KOG_cat)

for (i in 1:length(KOG_cat)) {
  KOG_cat_table[i, 2] <- sum(grepl(KOG_cat[i], jgikog$kog_sub_cat))
}

names(KOG_cat_table) <- c("KOG_cat", "Counts_array")
# write.csv(KOG_cat_table, file="KOG_subcat_table_on-array.csv")  # not required as separate doc


# Import data for DEgenes for KOG enrichment! -----------------------------------------------------------------
# Creat the list for all the modules and give them name 

ff <- list.files(path="./DE_preserved_P", full.names=TRUE)
myfilelist <- lapply(ff, read.table, col.names = "Gene") # header = F
names(myfilelist) <- list.files(path="./DE_preserved_P", full.names=FALSE)
for (i in 1:length(myfilelist)) {
names(myfilelist[[i]]) <- "gene" 
print(dim(myfilelist[[i]]))
}

length(ff) #432

# Filtration step ---------------------------------------------------------
DE_genes_kog <- list()
for (i in 1:length(myfilelist)) {
  DE_genes_kog[[i]] <- na.omit(merge(myfilelist[[i]], jgikog, by="gene", all.x=T))
}

names(DE_genes_kog) <- list.files(path="./DE_preserved_P", full.names=FALSE)
str(DE_genes_kog)

# counting
DE_KOG_JGI_count <- list()
for (i in 1:length(myfilelist)) {
  DE_KOG_JGI_count[[i]] <- as.data.frame(table(DE_genes_kog[[i]]$kog_sub_cat))[-1,]  # KOG_JGI replaced kog_sub_cat 
}

names(DE_KOG_JGI_count) <- list.files(path="./DE_preserved_P", full.names=FALSE) # 266

for (i in 1:length(myfilelist)) {
  if (class(DE_KOG_JGI_count[[i]]) == "data.frame"){  # needed to add this
      names(DE_KOG_JGI_count[[i]]) <- c("DE_KOG_JGI_count", "Counts_DE")
  }}

# Loop for main -----------------------------------------------------------
DE_KOG_cat <- LETTERS[c(1:23, 25:26)]

DE_KOG_cat_list <-  list()
for (i in 1:length(myfilelist)) {
DE_KOG_cat_table <- as.data.frame(DE_KOG_cat) # this stay 
for (j in 1:length(DE_KOG_cat)) {
  DE_KOG_cat_table[j, 2] <- sum(grepl(DE_KOG_cat[j], DE_genes_kog[[i]]$kog_sub_cat))
}
names(DE_KOG_cat_table) <- c("DE_KOG_cat", "Counts_DE")
DE_KOG_cat_list[[i]] <- DE_KOG_cat_table
}
names(DE_KOG_cat_list) <- list.files(path="./DE_preserved_P", full.names=FALSE)


# Write cvs DE_KOG_cat_tables ---------------------------------------------
d1_DE_KOG_subcat_table <- "DE_KOG_subcat_table"
if (!file.exists(d1_DE_KOG_subcat_table)) dir.create(d1_DE_KOG_subcat_table)
dir_temp_DE_KOG_subcat_table <- paste0("./", d1_DE_KOG_subcat_table, "")
setwd(dir_temp_DE_KOG_subcat_table)

for (i in 1:length(myfilelist)) {
  pre <- "DE_KOG_subcat_table_125"
  myname <- names(DE_KOG_cat_list[i])
  post <- ".csv"
  write.csv(DE_KOG_cat_list[[i]], file=paste0(pre, myname, post)) # not required as separate doc
}
setwd(dir_original)

# Preparation for Fisher test ---------------------------------------------
DE_KOG_tes_list <-  list()   # DE_KOG_tes_list was tes previously
for (i in 1:length(myfilelist)) {
DE_KOG_tes_list[[i]] <- cbind(KOG_cat_table[,1:2], DE_KOG_cat_list[[i]][,2])  # check this last one with DBwhy they are equal?
colnames(DE_KOG_tes_list[[i]]) <- c("KOG_sub_cat", "Counts_array", "Counts_DE")
DE_KOG_tes_list[[i]]$DE_others <- sum(DE_KOG_tes_list[[i]]$Counts_DE) - DE_KOG_tes_list[[i]]$Counts_DE
DE_KOG_tes_list[[i]]$ns_in_sub_cat <- DE_KOG_tes_list[[i]]$Counts_array - DE_KOG_tes_list[[i]]$Counts_DE  #check this as they cancell each other
DE_KOG_tes_list[[i]]$ns_others <- (sum(DE_KOG_tes_list[[i]]$Counts_array) - sum(DE_KOG_tes_list[[i]]$Counts_DE)) - DE_KOG_tes_list[[i]]$ns_in_sub_cat
DE_KOG_tes_list[[i]]$sum <- DE_KOG_tes_list[[i]]$Counts_DE + DE_KOG_tes_list[[i]]$DE_others + DE_KOG_tes_list[[i]]$ns_in_sub_cat + DE_KOG_tes_list[[i]]$ns_others
DE_KOG_tes_list[[i]]
}
names(DE_KOG_tes_list) <- list.files(path="./DE_preserved_P", full.names=FALSE)

# Write csv subcat_table ---------------------------------------------------------------
d1_DE_KOG_subcat_table_152 <- "DE_KOG_subcat_table_152"
if (!file.exists(d1_DE_KOG_subcat_table_152)) dir.create(d1_DE_KOG_subcat_table_152)
dir_temp_DE_KOG_subcat_table_152 <- paste0("./", d1_DE_KOG_subcat_table_152, "")
setwd(dir_temp_DE_KOG_subcat_table_152)

for (i in 1:length(myfilelist)) {
  pre <- "DE_KOG_subcat_table_152"
  myname <- names(DE_KOG_tes_list[i])
  post <- ".csv"
  write.csv(DE_KOG_tes_list[[i]], file=paste0(pre, myname, post)) # not required as separate doc
}
setwd(dir_original)

# Test fisher -------------------------------------------------------------
test_fisher_tes_list <- list()
p.value_tes_list <- list()
odds.ratio_tes_list <- list()
q_tes_list <- list()
final_tes_list <- list()

for (i in 1:length(myfilelist)) {
test_fisher_tes_list[[i]] <- DE_KOG_tes_list[[i]][, 3:6]
p.value_tes_list[[i]] <- apply(test_fisher_tes_list[[i]], 1, function(x) fisher.test(matrix(x, nr=2))$p.value)
odds.ratio_tes_list[[i]] <- apply(test_fisher_tes_list[[i]], 1, function(x) fisher.test(matrix(x, nr=2))$estimate)

# multiple testing correction according to Benjamini-Hochberg, 1995
q_tes_list[[i]] <- p.adjust(p.value_tes_list[[i]], method="fdr")
final_tes_list[[i]] <- cbind(DE_KOG_tes_list[[i]][,1], test_fisher_tes_list[[i]], p.value_tes_list[[i]], odds.ratio_tes_list[[i]], q_tes_list[[i]])

colnames(final_tes_list[[i]]) <- c("KOG_sub_cat",  "Counts_DE", "DE_others", "ns_in_sub_cat",
  "ns_others", "p.value", "odds.ratio", "q" )
str(final_tes_list[[i]])
}
names(final_tes_list) <- list.files(path="./DE_preserved_P", full.names=FALSE)


end_tes_list <- list()
for (i in 1:length(myfilelist)) {
end_tes_list[[i]] <- arrange(final_tes_list[[i]], q_tes_list[[i]])
}
names(end_tes_list) <- list.files(path="./DE_preserved_P", full.names=FALSE)

# Cutoff q test -----------------------------------------------------------
cutoff_q_tes <- 0.05   # adjust cut-off if required
new_tes_list <- list()
new_matched_tes <- list()

for (i in 1:length(myfilelist)) {
new_tes_list[[i]] <- subset(end_tes_list[[i]], q < cutoff_q_tes)
new_matched_tes[[i]] <- merge(new_tes_list[[i]], KOG, by.x="KOG_sub_cat", by.y="KOG_sub_cat", all.x=T)
}
names(new_matched_tes) <- list.files(path="./DE_preserved_P", full.names=FALSE)

# KOG_sub_cat_enrichment ---------------------------------------------------------------------
d1_KOG_sub_cat_enrichment_qvalue005_206 <- "DE_KOG_sub_cat_enrichment_qvalue005_206"
if (!file.exists(d1_KOG_sub_cat_enrichment_qvalue005_206)) dir.create(d1_KOG_sub_cat_enrichment_qvalue005_206)
dir_temp_KOG_sub_cat_enrichment_qvalue005_206 <- paste0("./", d1_KOG_sub_cat_enrichment_qvalue005_206, "")
setwd(dir_temp_KOG_sub_cat_enrichment_qvalue005_206)

for (i in 1:length(myfilelist)) {
  pre <- "KOG_sub_cat_enrichment_qvalue005_206"
  myname <- names(new_matched_tes[i])
  post <- ".csv"
  write.csv(new_matched_tes[[i]], file=paste0(pre, myname, post)) # not required as separate doc
}
setwd(dir_original)


# testes_list("qplot_p.value_q_KOG_sub_cat.value.testes_list")
#par(mfrow=c(2,2)) # this does not work for ggplot2
#qplot(q, binwidth = 0.05)
#qplot(p.value, binwidth = 0.05)
# dev.off()


# Part 2 KOG analysis for main categories (KOG_cat is main) ------------------------------------------------------------------
testes_list <- list()

for (i in 1:length(myfilelist)) {
testes_list[[i]] <- data.frame(matrix(NA, nrow=4, ncol=7))
colnames(testes_list[[i]]) <- c("KOG_cat", "Counts_array", "Counts_DE", "DE_others", "ns_in_cat", "ns_others", "sum")
testes_list[[i]]$KOG_cat <- c("INFORMATION_STORAGE_AND_PROCESSING", 
  "CELLULAR_PROCESSES_AND_SIGNALING", "METABOLISM", "POORLY_CHARACTERIZED")

# A, B, J, K, L
inf <- colSums(DE_KOG_tes_list[[i]][c(1:2,10:12), 2:3], dims=1)   #### CHECK TES here with DB 

# M, N, O, T, U, V, W, Y, Z 
cel <- colSums(DE_KOG_tes_list[[i]][c(13:15,20:25), 2:3], dims=1)

# C, D, E, F, G, H, I, P, Q
met <- colSums(DE_KOG_tes_list[[i]][c(3:9, 16:17), 2:3], dims=1)

# R, S
poo <- colSums(DE_KOG_tes_list[[i]][c(18:19), 2:3], dims=1)

testes_list[[i]][1, 2:3] <- inf
testes_list[[i]][2, 2:3] <- cel
testes_list[[i]][3, 2:3] <- met
testes_list[[i]][4, 2:3] <- poo

testes_list[[i]]

### Preparation for Fisher testes
testes_list[[i]]$DE_others <- sum(DE_KOG_tes_list[[i]]$Counts_DE) - testes_list[[i]]$Counts_DE  ## DE_KOG_tes_list was tes before
testes_list[[i]]$ns_in_cat <- testes_list[[i]]$Counts_array - testes_list[[i]]$Counts_DE
testes_list[[i]]$ns_others <- (sum(testes_list[[i]]$Counts_array) - sum(testes_list[[i]]$Counts_DE)) - testes_list[[i]]$ns_in_cat
testes_list[[i]]$sum <- testes_list[[i]]$Counts_DE + testes_list[[i]]$DE_others + testes_list[[i]]$ns_in_cat + testes_list[[i]]$ns_others

}

testes_list[[i]]
names(testes_list) <- list.files(path="./DE_preserved_P", full.names=FALSE)

# Write the csv here ------------------------------------------------------
d1_DE_KOG_cat_table_268 <- "DE_KOG_cat_table_268"
if (!file.exists(d1_DE_KOG_cat_table_268)) dir.create(d1_DE_KOG_cat_table_268)
dir_temp_DE_KOG_cat_table_268 <- paste0("./", d1_DE_KOG_cat_table_268, "")
setwd(dir_temp_DE_KOG_cat_table_268)

for (i in 1:length(myfilelist)) {
  pre <- "DE_KOG_cat_table_267"
  myname <- names(testes_list[i])
  post <- ".csv"
  write.csv(testes_list[[i]], file=paste0(pre, myname, post)) # not required as separate doc
}
setwd(dir_original)

# empty lists 
test_fisher_testes_list <- list()
p.value_testes_list <- list()
odds.ratio_testes_list <- list()
q_testes_list <- list()
final_testes_list <- list()

for (i in 1:length(myfilelist)) {
test_fisher_testes_list[[i]] <- testes_list[[i]][, 3:6]
p.value_testes_list[[i]] <- apply(test_fisher_testes_list[[i]], 1, function(x) fisher.test(matrix(x, nr=2))$p.value)
odds.ratio_testes_list[[i]] <- apply(test_fisher_testes_list[[i]], 1, function(x) fisher.test(matrix(x, nr=2))$estimate)

# multiple testing correction according to Benjamini-Hochberg, 1995
q_testes_list[[i]] <- p.adjust(p.value_testes_list[[i]], method="fdr")
final_testes_list[[i]] <- cbind(testes_list[[i]][,1], test_fisher_testes_list[[i]], p.value_testes_list[[i]], odds.ratio_testes_list[[i]], q_testes_list[[i]])

colnames(final_testes_list[[i]]) <- c("KOG_cat",  "Counts_DE", "DE_others", "ns_in_cat", 
  "ns_others", "p.value", "odds.ratio", "q" )
str(final_testes_list[[i]])
}

names(final_testes_list) <- list.files(path="./DE_preserved_P", full.names=FALSE)

# Final list 
end_testes_list <- list()
for (i in 1:length(myfilelist)) {
end_testes_list[[i]] <- arrange(final_testes_list[[i]], q_testes_list[[i]])
}
names(end_testes_list) <- list.files(path="./DE_preserved_P", full.names=FALSE)

# Csv down the files  ---------------------------------------------------
d1_end_testes_list <- "end_testes_list"
if (!file.exists(d1_end_testes_list)) dir.create(d1_end_testes_list)
dir_temp_end_testes_list <- paste0("./", d1_end_testes_list, "")
setwd(dir_temp_end_testes_list)

for (i in 1:length(myfilelist)) {
  pre <- "DE_KOG_cat_table_267"
  myname <- names(end_testes_list[i])
  post <- ".csv"
  write.csv(end_testes_list[[i]], file=paste0(pre, myname, post)) # not required as separate doc
}
setwd(dir_original)

# Adding cutoff for q value
cutoff_q_testes_list <- 0.05
withcutoff <- list()
for (i in 1:length(myfilelist)) {
withcutoff[[i]] <- subset(end_testes_list[[i]], q < cutoff_q_testes_list)
}
names(withcutoff) <- list.files(path="./DE_preserved_P", full.names=FALSE)

# Write the cvs -----------------------------------------------------------
d1_KOG_cat_enrichment_with_cutoff <- "KOG_cat_enrichment_with_cutoff"
if (!file.exists(d1_KOG_cat_enrichment_with_cutoff)) dir.create(d1_KOG_cat_enrichment_with_cutoff)
dir_temp_KOG_cat_enrichment_with_cutoff <- paste0("./", d1_KOG_cat_enrichment_with_cutoff, "")
setwd(dir_temp_KOG_cat_enrichment_with_cutoff)

for (i in 1:length(myfilelist)) {
  pre <- "KOG_cat_enrichment_with_cutoff"
  myname <- names(withcutoff[i])
  post <- ".csv"
  write.csv(withcutoff[[i]], file=paste0(pre, myname, post)) # not required as separate doc
}

setwd(dir_original)
