# ROH, Load, and Demography among _Ursus americanus_, _U. arctos_, and _U. maritimus_ 
Code for processing WGS data from bears for ROH calling, variant annotation, and identifying locations of variants with respect to ROHs.

Preprocessing:
   - Downloading _U. maritimus_ samples from NCBI; SRA info for others in table S1
   - Mapping reads to reference genomes
   - Marking duplicates, adding/replacing read groups using GATK
   - Indexing with SAMtools
   - Calling variants & filtering with GATK
   - Selecting scaffolds with PLINK (& removing X & Y chromosomes)

ROH analyses and heterozygosity estimates:
   - ROH calling using PLINK
   - HO estimated with VCFtools

Variant annotation and genomic location:
   - SnpEff annotation
   - Heterozygous or homozygous state of variants via fitlering in BCFtools
   - Location of variants with respect to ROHs using BEDTools and custom scripts
   - Cactus assembly alignment (for _U. americanus_)

