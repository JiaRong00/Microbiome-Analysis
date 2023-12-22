#Load the packages
library(dada2); library(DECIPHER); library(phangorn); library(phyloseq)

#store directory path
pathF="./FWD"
pathR="./REV"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(pathF,  pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(pathR, pattern="_R2_001.fastq", full.names = TRUE))

#Check that all forward and reverse files are present
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

#Extract sample names, assuming file names have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.namesR<- sapply(strsplit(basename(fnRs), "_"), `[`, 1)

#Forward & Reverse Read quality 
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in a new folder named filtered
filtFs <- file.path(pathF, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(pathR, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Read filtering
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 

#Error rate model for forward & reverse read estimates.
set.seed(100)
errF <- learnErrors(filtFs,nbases = 1e8, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(filtRs, nbases = 1e8, multithread=TRUE, randomize=TRUE)


#Sample inference(De-noising)
getN <- function(x) sum(getUniques(x))
trackF <- c()
trackR <- c()
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sample in sample.names) {
  cat("Processing:", sample, "\n") #Show which sample is being merged
  dadaFs <- dada(filtFs[[sample]], err=errF, multithread=TRUE)
  dadaRs <- dada(filtRs[[sample]], err=errR, multithread=TRUE)
  trackF <- append(trackF, getN(dadaFs))
  trackR <- append(trackR, getN(dadaRs))
  merger <- mergePairs(dadaFs, filtFs[[sample]], dadaRs, filtRs[[sample]], verbose=TRUE)
  mergers[[sample]] <- merger
}

#Construct sequence Table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "seqtab.rds")

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) 
saveRDS(seqtab.nochim,"seqtab_final.rds") 
write.csv(seqtab.nochim,"seqtab.csv", row.names=TRUE)

#Check reads lost through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, trackF, trackR, sapply(mergers, getN), rowSums(seqtab.nochim), round(rowSums(seqtab.nochim)/out[,1]*100, 1))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "percentage_retained")
rownames(track) <- sample.names
saveRDS(track, "track_reads_final.rds")
write.csv(track, "track_reads_final.csv")

#DECIPHER for taxonomy
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("SILVA_SSU_r138_2019.RData") 
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest

# Assigning ASV IDs (ASV_1, ASV_2 etc)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

#FASTA file containing DNA seq of each ASV
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_sequence_all.fa")

#count table (number of reads for each ASV for each sample):
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_table_all.tsv", sep="\t", quote=F, col.names=NA)

#Taxonomy table
asv_tax <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

#Construct phylogenetic tree
ASVs.nochim = DNAStringSet(colnames(seqtab.nochim))
names(ASVs.nochim) = paste0("ASV_", 1:ncol(seqtab.nochim))
alignment = AlignSeqs(ASVs.nochim, anchor=NA, processors=30)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
save(fitGTR, file="fitGTR.RData")
detach("package:phangorn", unload=TRUE)
plot(treeNJ, main="NJ")
ape::write.tree(treeNJ, file='ASVs_tree.txt')

#Read metadata file 
samples.out <- rownames(seqtab.nochim)
samdf <- read.csv("metadata.csv")
rownames(samdf) <- samples.out
colnames(seqtab.nochim)<-paste0("ASV_", 1:ncol(seqtab.nochim)) 

#Combining into a phyloseq object. 
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(asv_tax),
               phy_tree(fitGTR$tree),
               refseq(ASVs.nochim))
ps
save(ps, file="phyloseq.RData")


