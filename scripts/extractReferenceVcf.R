# Script to extract references from the GRCH Reference VCF
library(data.table)
library(bedr)

REFERENCE <- '/mnt/humanSV/data/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.gz'
ASSEMBLY_MINI <- '/home/dkoppstein/ngschool/src/fantavier/10_survivor/assemblytics_sv.vcf'

vcm <- read.vcf(ASSEMBLY_MINI)
vcr <- read.vcf(REFERENCE)


dtm <- as.data.table(vcm$vcf)
dtr <- as.data.table(vcr$vcf)



