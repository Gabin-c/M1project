#During all this analyse we gonna use function made in function_dds
source("function_dds.R")


### import data ----
#Firstly we import the 3 data table which we will use during this analysis. 
#The "counts_table" contains the counting of read mapped on each genes by samples,
#"airway_metadata" contains the information about the samples design 
#and "anno" contains the informations about the genes. 
counts_table <- read.csv("Data/airway_scaledcounts.csv") 
airway_metadata <- read.csv("Data/airway_metadata.csv")
anno <- read.csv("Data/annotables_grch38.csv")
anno <- anno %>% select(ensgene,symbol)

### Create the dds object ---
dds <- DESeqDataSetFromMatrix(counts_table,colData=airway_metadata,design = ~dex,tidy = TRUE)
# Set reference of experience, here "control"
colData(dds)$dex <- relevel(colData(dds)$dex , ref="control")
# To display experiment design.
colData(dds) 
# To display column which biological condition is set.
design(dds)

### Explorate data ----
#Transformation of dds counts table in data frame to use ggplot package
count_table_dds <- as.data.frame(counts(dds))
#Observation of count distribution 
par(mfrow=c(1,2))
for( i in 1:8){
  p <- ggplot(data=count_table_dds, aes(log(count_dds[,i]+1))) + geom_histogram(breaks=seq(0,14,1),col="black",fill="grey")+theme_light()+labs(title=colnames(count_table_dds)[i], x="Count value (number of read by genes) in log(count+1)",y="Count frequency") + theme_bw()
 plot(p)
}
#Number for Null for each sample
apply(count_table_dds, 2 ,FUN = function(x) sum(x==0))
#Observation of depth for each sample using deph.plot function
depth <- colSums(count_table_dds)
depth <- as.data.frame(depth)
depth$sample <- row.names(depth)
ggplot(depth, aes( x=sample ,y=depth))+ geom_bar(stat="identity",col="black", fill="white")+labs(title = "Depth of each sample", x="Sample", y="Depth")+theme_bw()

#Differential expression analysis and normalization
dds <- DESeq(dds)
#Extraction of the DE result
res <- results(dds)
# The scale factor resulting of normalization are stock in size factor of our design.
colData(dds) 
