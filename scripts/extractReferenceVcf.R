# Script to extract references from the GRCH Reference VCF
library(data.table)
library(bedr)

REFERENCE <- '/mnt/humanSV/data/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.gz'
ASSEMBLY_MINI <- '/home/dkoppstein/ngschool/src/fantavier/10_survivor/assemblytics_sv.vcf'

vcm <- read.vcf(ASSEMBLY_MINI)
vcr <- read.vcf(REFERENCE)

# Convert to data.table
dtm <- as.data.table(vcm$vcf)
dtr <- as.data.table(vcr$vcf)

# Extract Info fields
dtm_info <- as.data.table(t(strsplit(dtm$INFO, ";")))

dtm_info <- dtm[1:nrow(dtm), tstrsplit(INFO, ";"),]
dtr_info <- dtr[1:nrow(dtr), tstrsplit(INFO, ";"),]

dtm_info_split <- dtm_info[, .(SVTYPE = tstrsplit(V2, '=' )[2], CHR2 = tstrsplit(V4, "=")[2], END = tstrsplit(V5, "="), SVLEN = tstrsplit(V6, "=")),]

dtr_info_split <- dtm_info[, .(SVTYPE = tstrsplit(V2, '=' )[2], CHR2 = tstrsplit(V4, "=")[2], END = tstrsplit(V5, "="), SVLEN = tstrsplit(V6, "=")),]

