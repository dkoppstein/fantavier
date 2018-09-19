# fantavier
Awesome SV analysis pipeline

# How to use
1. symlink the data directory in /mnt to the fantavier folder

```
cd fantavier
ln -s /mnt/humanSV/data .
```

2. Then run using e.g.

```
snakemake -j 2 -- 5_racon/assembly_consensus.fasta
```

Where the last file is file you want to create. 