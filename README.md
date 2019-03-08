# Pipeline for phasing experiments from GIAB data

Note that this Snakemake pipeline contains paths to the reference genome and to WhatsHap that you
might need to adapt to use it.

## Phasing from RTG and 10x data

The results of the following analysis can be found at
```
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_MPI_whatshap_08232018/
```

To produce a chromosome-length phasing of small variants for the Ashkenazim trio, we combined variant
calls from Real Time Genomics [1] with phased blocks produced by 10x Genomics [2,3,4]. The single sample
10x Genomics VCF files were combined into multi-sample VCF using bcftools and all VCFs were split by 
chromosome (to facilitate easy parallelization with Snakemake). Then, WhatsHap (version 0.15+14.ga105b78, [5])
was used in pedigree-aware mode [6] using the following command line:

```
whatshap phase --ped AJ.ped --indels --reference hg19.fasta rtg.vcf 10x-merged.vcf | bgzip > output.vcf
```

[1] ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/Rutgers_IlluminaHiSeq300X_rtg_11052015/rtg_allCallsV2.vcf.gz

[2] ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_09302016/NA24143_hg19/NA24143_hg19_phased_variants.vcf.gz

[3] ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_09302016/NA24149_hg19/NA24149_hg19_phased_variants.vcf.gz

[4] ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_09302016/NA24385_hg19/NA24385_hg19_phased_variants.vcf.gz

[5] M. Martin, M. Patterson, S. Garg, S. O. Fischer, N. Pisanti, G. W. Klau, A. Schönhuth, T. Marschall.
WhatsHap: fast and accurate read-based phasing. bioRxiv, 2016. DOI: 10.1101/085050

[6] S. Garg, M. Martin, T. Marschall. Read-Based Phasing of Related Individuals. Bioinformatics (Proceedings of ISMB), 32, pp. i234-i242, 2016. DOI: 10.1093/bioinformatics/btw276
