#!/usr/bin/env Rscript


HAPGEN <- "~/Software/hapgen2/hapgen2"
HAPLO <- "~/Data/HM3/CEU.chr10.hap"
LEGEND <- "~/Data/HM3/hapmap3.r2.b36.chr10.legend"
MAP <- "~/Data/HM3/genetic_map_chr10_combined_b36.txt"
OUT <- "sim.out"
POPSIZE <- "11418"

# sample size, N/2 cases, N/2 controls
N <- 500

# risk ratio per allele dosage
alleleRR <- c(1.1, 1.2, 1.3, 1.4, 1.5, 2)

# reference allele
refAllele <- 0

################################################################################
# Split genotypes into blocks

blocksize <- 500

# read legend
legend <- read.table(LEGEND, sep="", header=TRUE, stringsAsFactors=FALSE)

# read haplotype
hapd <- scan(HAPLO, what=integer())
hap <- matrix(hapd, byrow=TRUE, nrow=nrow(legend))

hapg <- apply(hap, 1, function(h) {
   apply(matrix(h, byrow=TRUE, ncol=3), 1, paste, collapse="")
})
hapg <- t(hapg)

# read map
map <- read.table(MAP, header=TRUE, sep="")

blkid <- ceiling(1:nrow(legend) / blocksize)
nblk <- length(unique(blkid))
if(nblk > 1000)
   stop("too many blocks for nblk")

blkdir <- "blocks"
dir.create(blkdir)

for(i in 1:nblk)
{
   w <- blkid == i
   HAPLO2 <- tail(strsplit(HAPLO, "/")[[1]], 1)
   LEGEND2 <- tail(strsplit(LEGEND, "/")[[1]], 1)

   hapblknm <- sprintf("%s/%s_%04d", blkdir, HAPLO2, i)
   legblknm <- sprintf("%s/%s_%04d", blkdir, LEGEND2, i)

   write.table(hap[w, ], file=hapblknm, quote=FALSE, col.names=FALSE,
	 row.names=FALSE, sep=" ")

   write.table(legend[w, ], file=legblknm, quote=FALSE, col.names=TRUE,
	 row.names=FALSE, sep=" ")
}

# There are many ambiguous genotypes in this data (like 101 or 011), we ignore
# that, as we're only interested in strict monomorphism
monomorphic <- apply(hapg, 1, function(x) {
   length(table(x)) == 1
})


################################################################################
# Simulate each block using the same causal SNPs, one SNP per block

causal <- sapply(1:nblk, function(i) {
   sample(legend$position[blkid == i & !monomorphic], 1)
})

snpRR <- sample(alleleRR, length(causal), replace=TRUE)


for(i in 1:nblk)
{
   # risks: <allele BP> <risk allele> <hetero risk> <homo risk>
   risks <- sprintf("%s %s %s %s", causal[i], refAllele,
	 snpRR[i], snpRR[i]^2)

   hapblknm <- sprintf("%s/%s_%04d", blkdir, HAPLO2, i)
   legblknm <- sprintf("%s/%s_%04d", blkdir, LEGEND2, i)
   outblknm <- sprintf("%s_%04d", OUT, i)

   s <- sprintf("%s -h ../%s -l ../%s -m %s -dl %s -n %s %s -o %s -no_haps_output -Ne 11418",
         HAPGEN, hapblknm, legblknm, MAP, risks, N/2, N/2, outblknm, POPSIZE)

   dir <- sprintf("block_%04d", i)
   dir.create(dir)
   old.dir <- getwd()
   setwd(dir)
   system(s)
   system(sprintf("~/Data/HapGen/hapgen2tped.py %s", outblknm))
 
   write.table(as.integer(causal[i]),
	 file=sprintf("causal_%04d.txt", i), row.names=FALSE,
	 col.names=FALSE, quote=FALSE)
   setwd(old.dir)
}

################################################################################
# Combine blocks into one dataset using PLINK

out.tped <- sprintf("%s.tped", OUT)

files <- sapply(1:nblk, function(i) {
   sprintf(fmt="block_%04d/sim.out_%04d.tped", i, i)
})
s <- paste("cat", paste(files, collapse=" "), ">", out.tped)
system(s)

# Convert to BED file
# All FAM files are identical, so use the first
system(sprintf(
   "plink --noweb --tped %s --tfam %s --make-bed --out %s",
      out.tped, "block_0001/sim.out_0001.fam", OUT
))

p <- legend$position %in% causal
d <- data.frame(
   Index=which(p),
   RS=legend$rs[p],
   BP=causal,
   RiskRatio=snpRR
)

write.table(d, quote=FALSE,
   file="causal.txt", row.names=FALSE, col.names=TRUE)

bim <- read.table(sprintf("%s.bim", OUT), header=FALSE, sep="",
      stringsAsFactors=FALSE)
causal.bin <- as.integer(bim$V2 %in% d$RS)

write.table(causal.bin, file="causal_bin.txt", row.names=FALSE,
      col.names=FALSE, quote=FALSE)


