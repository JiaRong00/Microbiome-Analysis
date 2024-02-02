#Load the packages
library(dada2); library(DECIPHER); library(phangorn); library(phyloseq); library(ShortRead); library(Biostrings)

#store directory path. If multiple sequencing runs were carried out, process each run separately up till the sample infeference step.
pathF="./Run1/FWD"
pathR="./Run1/REV"

# Assuming forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(pathF,  pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(pathR, pattern="_R2_001.fastq", full.names = TRUE))

#Check that all forward and reverse files are present
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

#Primers
FWD <- "ACCTGCGGARGGATCA"  ## CHANGE ME to your forward primer sequence
REV <- "GAGATCCRTTGYTRAAAGTT"  ## CHANGE ME to your reverse primer sequence
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
        RevComp = Biostrings::reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#Pre-filter to remove Ns
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#Check for primers
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#If primers are of fixed length and only present at the start/end of the reads, we can use trimleft or trimright parameter in trimandfilter function of dada2 to remove

#Remove primers using cutadapt
cutadapt <- "/usr/local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine. If using windows single-file executable, remember to add '.exe'   
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#Sanity Check
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#Extract sample names, assuming file names have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(cutFs), "_"), `[`, 1)
sample.namesR<- sapply(strsplit(basename(cutRs), "_"), `[`, 1)

#Forward & Reverse Read quality. Usually the plots for samples in the same run should look similar, hence not necesaary to look at all plots.
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2]) #Normally the plot quality is worse for reverse reads. May choose to proceed with just forward reads only


# Place filtered files in a new folder named filtered
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#Read filtering. truncate at length at which read quality fall below 30, but ensure there is minimum of of 12bp for overlap.
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(240,160), #truncate at 240 for forward reads, 160 for reverse reads.
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) #multithread=FALSE for windows.
saveRDS(out, "./Run1/out.rds")
#May decrease maxEE to reduce computation time for sample inference step. 

#Error rate model for forward & reverse read estimates. 
set.seed(100)
errF <- learnErrors(filtFs,nbases = 1e8, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(filtRs, nbases = 1e8, multithread=TRUE, randomize=TRUE)
plotErrors(errF, nominalQ=TRUE) #black line should be a good fit to the data
#However, if binned quality reads are used then the error model will likely not be a good fit. In this case, refer to @JacobRPrice's solution on https://github.com/benjjneb/dada2/issues/1307
saveRDS(errF, "./Run1/errF.rds")
saveRDS(errR, "./Run1/errR.rds")
        
#Sample inference(De-noising)
names(filtFs) <- sample.names
names(filtRs) <- sample.names
getN <- function(x) sum(getUniques(x))
trackF <- c()
trackR <- c()
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sample in sample.names[1:]) { #Can change number in bracket to process eg. first 10 samples and continue next day if the sample inference step takes too long
  cat("Processing:", sample, "\n") #Show which sample is being merged
  dadaFs <- dada(filtFs[[sample]], err=errF, multithread=TRUE)
  dadaRs <- dada(filtRs[[sample]], err=errR, multithread=TRUE)
  trackF[sample] <- getN(dadaFs)
  trackR[sample]<- getN(dadaRs)
  mergers[[sample]] <- mergePairs(dadaFs, filtFs[[sample]], dadaRs, filtRs[[sample]], verbose=TRUE) #Change Maxmismatch to 1 or 2 from default 0 if needed
}

#Construct sequence Table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "./Run1/seqtab.rds")

#If multiple sequencing runs conducted, combine the sequence tables for each run before proceeding
track <- cbind(out, trackF, trackR, sapply(mergers, getN))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
saveRDS(track, "track_reads.rds")
run1 <-readRDS("./Run1/seqtab.rds")
run2 <-readRDS("./Run2/seqtab.rds")
st.all <- mergeSequenceTables(run1, run2)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) 
saveRDS(seqtab.nochim,"seqtab_final.rds") 
write.csv(seqtab.nochim,"seqtab.csv", row.names=TRUE)

#Check reads lost through pipeline
nonchim <- rowSums(seqtab.nochim)
percentage_retained <- round(rowSums(seqtab.nochim)/st.all[,1]*100, 1)
track <- cbind(st.all, nonchim, percentage_retained)
saveRDS(track, "track_reads_final.rds")
write.csv(track, "track_reads_final.csv")

#DECIPHER for taxonomy. Note that this method of classification is stricter than assignTaxonomy() in dada2
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("SILVA_SSU_r138_2019.RData") 
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest

#Alternatively, stick to assignTaxonomy in dada2
#taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
#taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v138.1.fa.gz")

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


