Step 1: Submit the script to inspect quality of fastqc file (01-fastqc.sh)
```
#!/bin/bash
#SBATCH -J fastqc
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p nocona
#SBATCH -a 1-4:1
#SBATCH --export=ALL

for name in $(cat names1.txt); do fastqc ${name}.fastq.gz; done
```
Step 2: Trimming the data based on quality and length filter. Use fastp for this. You could either
download it from conda or load the module if you hpc has it. (02-trim-S1.sh)

```
#!/bin/bash
#SBATCH -J s1_trimming
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p quanah
#SBATCH -a 1-4:1
#SBATCH --export=ALL


fastp \
  -i S1_F1_S1_R1_001.fastq \
  -I S1_F1_S1_R2_001.fastq \
  -o clean_l120_q30/cleaned_S1_F1_S1_R1_001.fastq \
  -O clean_l120_q30/cleaned_S1_F1_S1_R2_001.fastq \
  -q 20 \
  --length_required=120 \
```

step 3:Trimming for S4 pools.  (03-trim-s4.sh)
```
#!/bin/bash
#SBATCH -J s4_f
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p quanah
#SBATCH -a 1-4:1
#SBATCH --export=ALL


fastp \
  -i S4_F_S6_R1_001.fastq \
  -I S4_F_S6_R2_001.fastq \
  -o clean_l120_q30/cleaned_S4_F_S6_R1_001.fastq \
  -O clean_l120_q30/cleaned_S4_F_S6_R2_001.fastq \
  -q 20 \
  --length_required=120 \
```
Step 4:	Convert FASTQ file to for s1 and s4 uBAM files. I am submitting it as an array script and i have the name for fastq files of s1 and s4 pools and their parents saved as all_names.txt
(04-ubam.sh)
```
#!/bin/bash
#SBATCH -J ubam
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p nocona
#SBATCH -a 1-4:1
#SBATCH --export=ALL

prefixNum="$SLURM_ARRAY_TASK_ID"
prefixNumP="$prefixNum"p
prefix=`sed -n "$prefixNumP" all_names.txt`

gatk FastqToSam --java-options "-Xmx8G" \
    --FASTQ clean_l120_q30/"$prefix"_R1_001.fastq.gz \
    --FASTQ2 clean_l120_q30/"$prefix"_R2_001.fastq.gz \
    --OUTPUT uBAM/"$prefix".unmapped.bam \
    --READ_GROUP_NAME Silene_vulgaris \
    --SAMPLE_NAME "$prefix" \
    --LIBRARY_NAME NEBNext-"$prefix" \
    --PLATFORM_UNIT HT23LBBXX.6.1101 \
    --PLATFORM illumina \
    --SEQUENCING_CENTER OMRF \
    --RUN_DATE 2018-03-12T08:55:00+0500
```
Step 5:	Map the reads to reference genome for both s1 and s4 pools and parents of them. (05-abam.sh). We need gatk, samtools and bwamem for running this script. Also, please make sure
to index the reference genome (index generated from bwa, samtools and picard all are needed).
```
#!/bin/bash
#SBATCH -J abam_all
#SBATCH -N10
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p nocona
#SBATCH -a 1-4:1
#SBATCH --export=ALL

prefixNum="$SLURM_ARRAY_TASK_ID"
prefixNumP="$prefixNum"p
prefix=`sed -n "$prefixNumP" all_names.txt`

gatk --java-options "-Dsamjdk.compression_level=5 -Xms3000m" SamToFastq \
    --INPUT uBAM/"$prefix".unmapped.bam \
    --FASTQ /dev/stdout \
    --INTERLEAVE true \
    --NON_PF true \
    | \
    bwa mem -K 100000000 -p -v 3 -t 16 \
    -Y SileneGenomeV3/S_vulgaris_v3.fasta  \
    /dev/stdin - \
    | \
    samtools view -1 - > aBAM/"$prefix".aligned.bam
```
Step 6: Merge aligned and unaligned BAM files (06-mbam.sh)
```
#!/bin/bash
#SBATCH -J merge
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p nocona
#SBATCH -a 1-4:1
#SBATCH --export=ALL

prefixNum="$SLURM_ARRAY_TASK_ID"
prefixNumP="$prefixNum"p
prefix=`sed -n "$prefixNumP" all_names.txt`

gatk --java-options "-Dsamjdk.compression_level=5 -Xms3000m" MergeBamAlignment \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ALIGNED_BAM aBAM/"$prefix".aligned.bam \
    --UNMAPPED_BAM uBAM/"$prefix".unmapped.bam \
     --OUTPUT mBAM/"$prefix".merged.bam \
    --REFERENCE_SEQUENCE SileneGenomeV3/S_vulgaris_v3.fasta \
    --PAIRED_RUN true \
    --SORT_ORDER 'unsorted' \
    --IS_BISULFITE_SEQUENCE false \
    --ALIGNED_READS_ONLY false \
    --CLIP_ADAPTERS false \
    --MAX_RECORDS_IN_RAM 2000000 \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --PROGRAM_RECORD_ID 'bwamem' \
    --PROGRAM_GROUP_VERSION '0.7.17-r1188' \
    --PROGRAM_GROUP_COMMAND_LINE 'bwa mem -K 100000000 -p -v 3 -t 16 -Y SileneGenomeV3/S_vulgaris_v3.fasta /dev/stdin -' \
    --PROGRAM_GROUP_NAME 'bwamem' \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --UNMAP_CONTAMINANT_READS true
```
Step 7: Sort the alignment by coordinates and fix tags (07-sbam.sh)
```
#!/bin/bash
#SBATCH -J sbam
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p quanah
#SBATCH -a 1-4:1
#SBATCH --export=ALL

prefixNum="$SLURM_ARRAY_TASK_ID"
prefixNumP="$prefixNum"p
prefix=`sed -n "$prefixNumP" all_names.txt`

gatk --java-options "-Dsamjdk.compression_level=5 -Xms4000m" SortSam \
    --INPUT mBAM/"$prefix".merged.bam \
    --OUTPUT /dev/stdout \
    --SORT_ORDER 'coordinate' \
    --CREATE_INDEX false \
    --CREATE_MD5_FILE false \
    | \
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms500m" SetNmMdAndUqTags \
    --INPUT /dev/stdin \
    --OUTPUT sBAM/"$prefix".sorted.bam \
    --CREATE_INDEX true \
    --CREATE_MD5_FILE true \
    --REFERENCE_SEQUENCE SileneGenomeV3/S_vulgaris_v3.fasta
```
step 8: Remove the reads that are optical duplicates. The output of this file was the input file for estimating allele frequency difference (08-dbam.sh)

```
#!/bin/bash
#SBATCH -J dbam
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p quanah
#SBATCH -a 1-4:1
#SBATCH --export=ALL

prefixNum="$SLURM_ARRAY_TASK_ID"
prefixNumP="$prefixNum"p
prefix=`sed -n "$prefixNumP" all_names.txt`


gatk --java-options "-Dsamjdk.compression_level=5 -Xms4000m" MarkDuplicates \
    --INPUT sBAM/"$prefix".sorted.bam \
    --OUTPUT dBAM/"$prefix".dedup.bam \
    --METRICS_FILE dBAM/"$prefix".dedup.metrics \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --ASSUME_SORT_ORDER 'coordinate' \
    --CREATE_INDEX true \
    --CREATE_MD5_FILE true

```
step 9: Generate a pileup file. I am showing for s4 pools only, but I did the same thing for s1 as well. (09-pileup.sh)

```
#!/bin/bash
#SBATCH -J pileup_s4
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p nocona
#SBATCH --export=ALL

samtools mpileup -B dBAM/cleaned_S4_H_S4.dedup.bam dBAM/cleaned_S4_F_S6.dedup.bam -f SileneGenomeV3/S_vulgaris_v3.fasta > mpileup_ref/S4_female_herm.mpileup
```
step 10: Create synchronized file that contains the allele frequencies after filtering for base quality. I did it with same code for s1 too. (10-sync.sh)
```
#!/bin/bash
#SBATCH -J sync_s1
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p quanah
#SBATCH -a 1-43:1
#SBATCH --export=ALL

java -ea -Xmx7g -jar popoolation2_1201/mpileup2sync.jar --input mpileup/S4_herm_female.mpileup --output mpileup/S4_herm_female_java.sync --fastq-type sanger --min-qual 20 --threads 8

```
step 10: Estimate allele frequ3ency difference for s4 pools (10-sllele_freq_diff.sh)
```

#!/bin/bash
#SBATCH -J 15_cov
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p nocona
#SBATCH --export=ALL

perl ../popoolation2_1201/snp-frequency-diff.pl --input S4_herm_female_java.sync --output-prefix filter_min_cov_15/S4_herm_female --min-count 6 --min-coverage 15 --max-coverage 100
```
step 11:	Fishers exact test for identifying significant allele frequency difference for s4 (11-fisher_test.sh)
```
#!/bin/bash
#SBATCH -J cov_15
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p nocona
#SBATCH --export=ALL


perl ../../popoolation2_1201/fisher-test.pl --input S4_herm_female_java.sync --output filter_cov_15/filtered_s4.sync.fet --min-count 6 --min-coverage 15 --max-coverage 100 --suppress-noninformative
```

step 12:	Exclude the site that are homozygous for parent in both pools
```

awk 'NR == FNR { patterns[$1]; next } !($1 in patterns)' snp_to_filter_out_pattern.txt filtered_s1.sync.fet > filtered_s1.sync.fet_homozygous.txt

awk 'NR == FNR { patterns[$1]; next } !($1 in patterns)' snp_to_filter_out_pattern.txt filtered_s4.sync.fet > filtered_s4.sync.fet_homozygous.txt
```
step 13:
A) Fst estimation for s1. Input file is the sync file generated by popoolation pooling software. Poolsize.txt is the text file containing pool sizes for our samples. In this case it is,
S1_herm_female_java.1   12
S1_herm_female_java.2   6

#fst code used
```
grenedalf/bin/grenedalf fst --sync-path /lustre/scratch/askhanal/Silene/mpileup/S1_herm_female_java.sync --reference-genome-dict-file /lustre/scratch/askhanal/Silene/SileneGenomeV3/S_vulgaris_v3.dict --filter-sample-min-coverage 15 --filter-sample-max-coverage 120 --filter-total-only-biallelic-snps --window-type single --method unbiased-nei --pool-sizes poolsize.txt
```
B) Fst estimation for s4. Input file is the sync file generated by popoolation pooling software. Poolsize.txt is the text file containing pool sizes for our samples. In this case, it is,
S4_herm_female_java.1   12
S4_herm_female_java.2   6

 #fst code used
```
grenedalf/bin/grenedalf fst --sync-path /lustre/scratch/askhanal/Silene/mpileup/S4_herm_female_java.sync --reference-genome-dict-file /lustre/scratch/askhanal/Silene/SileneGenomeV3/S_vulgaris_v3.dict --filter-sample-min-coverage 15 --filter-sample-max-coverage 120 --filter-total-only-biallelic-snps --window-type single --method unbiased-nei --pool-sizes poolsize_s4.txt
```
j)	I removed the sites with fst value of “Nan” and negative fst values. Also, with the following code remove the site in fst that is homozygous for parent
```
awk 'NR == FNR { patterns[$1]; next } !($1 in patterns)' snp_to_filter_out_pattern.txt fst_s4_0_1_new_co_v2.txt > fst_s4_0_1_new_co_v2_filtered_homozygous.txt

awk 'NR == FNR { patterns[$1]; next } !($1 in patterns)' snp_to_filter_out_pattern.txt fst_s4_0_1_new_co_v2.txt > fst_s4_0_1_new_co_v2_filtered_homozygous.txt
```

