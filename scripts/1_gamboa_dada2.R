library(dada2); packageVersion("dada2")

path <- "panama_16s_filtered/"
list.files(path)

fnFs <- sort(list.files(path, pattern="_fastp_R1_001.fastq.gz", full.names = TRUE)) #import forward reads
fnRs <- sort(list.files(path, pattern="_fastp_R2_001.fastq.gz", full.names = TRUE)) #import reverse reads
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #extract sample names


#plotQualityProfile(fnFs[1:2]) #plot quality score of forward reads
#plotQualityProfile(fnRs[1:2]) #plot quality score of reverse reads

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) #create object holding filtered F reads.
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) #create object holding filtered R reads
names(filtFs) <- sample.names #name
names(filtRs) <- sample.names #name

#This following step will take a few minutes. Get a drink or take a small break!
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, #specify input+output
                     truncLen=c(240,150),#where we are cutting off our sequences
                     trimLeft=c(0,50), 
                     maxN=0, maxEE=c(1,2), truncQ=2, rm.phix=TRUE, #default parameters
                     compress=TRUE, multithread=FALSE) #my multitread is set to FALSE bc I run on Windows. Set yours to TRUE if you're on a MAC

head(out)

errF <- learnErrors(filtFs, multithread=FALSE, nbases = 1e+09) #my multitread is set to FALSE bc I run on Windows. Set yours to TRUE if you're on a MAC
errR <- learnErrors(filtRs, multithread=FALSE, nbases = 1e+09) #my multitread is set to FALSE bc I run on Windows. Set yours to TRUE if you're on a MAC
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=FALSE) #merge paired reads
head(mergers[[1]]) # inspect the merger data.frame from the first sample

seqtab <- makeSequenceTable(mergers) #create sequence table that is imput for chimera removing function
dim(seqtab) #Check number of ASV's before chimera removal

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=FALSE) #my multitread is set to FALSE bc I run on Windows. Set yours to TRUE if you're on a MAC

dim(seqtab.nochim)#check number of ASV's after chimera removal

sum(seqtab.nochim)/sum(seqtab) 


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#This step takes a few minutes.

taxa <- assignTaxonomy(seqtab.nochim, "tax/GTDB_bac-arc_ssu_r86.fa.gz", multithread=FALSE) #my multitread is set to FALSE bc I run on Windows. Set yours to TRUE if you're on a MAC
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

seqs <- colnames(seqtab.nochim)
headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  headers[i] <- paste(">ASV", i, sep="_")
}

#fasta file
fasta <- c(rbind(headers, seqs)) #make fasta file
write(fasta, "panama_dada2_GTDB.fa") #write out as .fasta file

#count table
count_tab <- t(seqtab.nochim) #make count table
row.names(count_tab) <- sub(">", "", headers) #replace rownames with something that is easy to follow
write.csv(count_tab, "panama_dada2_GTDB_counts.csv", row.names=TRUE) #write out as .csv file


# tax table:
tax <- taxa #make taxonomy table
row.names(tax) <- sub(">", "", headers) #replace rownames with something that is easy to follow
write.csv(tax, "panama_dada2_GTDB_taxonomy.csv", row.names=TRUE) #write out as .csv file
write.csv(track, "reads_tracked_panama_GTDB.csv", row.names=TRUE) #write out as .csv file


