# Load libraries ----------------------------------------------------------
library(data.table)

setwd("~/Dropbox/KOG_Oct2017/rdsprojects2016frischd-01welcome_to_our_lovely_shared_folderall_ancient_recentModule_gene_list_txts")
dir_original <- getwd()

annot_path="./annotated" # Here we will be printing Annotated_files
path_sup="./needed"  # Here you need 2 text files

setwd(path_sup)
# Load supportive files
kog <-read.table("kogsubcatid.txt", sep="", header = T)
head(kog)
dim(kog)

# Reduced set of genes that are on array taking into account the problmes with probes
jgi_sub <-read.csv("Dappu1_JGI_Paralog_Tandem_DB_adjusted.csv", sep=",", header = T) ## Notice I removed GENE: for matching 
head(jgi_sub)
dim(jgi_sub)[1]  # 24679

setwd(dir_original)
# 
# mymodule <-read.table("black_module_gene_list.txt", sep="", header = F)
# head(mymodule)
# firstmerging <- merge(mymodule, jgi_sub, by.x="V1", by.y="gene")
# dim(firstmerging) # 920
# secondmerging <- merge(firstmerging, kog, by.x="KOG_JGI", by.y="kog_id", all.x = TRUE) # adjust TRUE or FALSE depending what you need NB! agree DB
# write.csv(path=annot_path, secondmerging, file=paste0(annotated, XXX)) # n
# dim(secondmerging)
# 

# Set parameters ----------------------------------------------------------
home.folder <- dir_original
module.file.pattern <- 'gene_list.txt'
module.files <- list.files(path=".", pattern = module.file.pattern)
module.files
myfilelist <- lapply(module.files, read.table, col.names = "Gene") # giving name

#Print new files
setwd(annot_path)
for (i in 1:length(myfilelist)) {
  pre <- "Annotated_"
  myname <- module.files[i]
  mymodule <- myfilelist[i]
  firstmerging <- merge(mymodule, jgi_sub, by.x="Gene", by.y="gene")
  dim(firstmerging)
  secondmerging <- merge(firstmerging, kog, by.x="KOG_JGI", by.y="kog_id", all.x = TRUE) # adjust TRUE or FALSE depending what you need NB! agree DB
  post <- ".csv"
  write.csv(secondmerging, file=paste0(pre, myname, post)) # not required as separate doc
}

setwd(dir_original)

