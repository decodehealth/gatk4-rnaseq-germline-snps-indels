# SNP Calling Pipeline
Author: Coleman Harris, Decode Health

### Inputs: 
- FASTQ R1 & R2 sample files
- Reference files
	- GTF 
	- FASTA
	- VCF files (Note: must be from same reference as GTF and FASTA)

### Outputs:
- SNP locations [.vcf]

# SETTING UP REFERENCES
Note: This only needs to be run once for a given set of reference files (e.g. GTF, FASTA, etc.)

## FASTA setup 
### Generate sequence dict
	java -jar /mnt/iquityazurefileshare1/IQuity/GATK_Testing/gatk-workflows/picard.jar \
		CreateSequenceDictionary \
		R=/mnt/iquityazurefileshare1/GTFs/SNPs/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		O=/mnt/iquityazurefileshare1/GTFs/SNPs/Homo_sapiens.GRCh38.dna.primary_assembly.dict

### Generate FASTA index file
	samtools faidx \
		Homo_sapiens.GRCh38.dna.primary_assembly.fa

------
## gtfToCallingIntervals  
### R script to read GTF
	Rscript --no-save -<<'RCODE'
		gtf = read.table("/mnt/iquityazurefileshare1/GTFs/Homo_sapiens.GRCh38.104.gtf", sep="\t")
		gtf = subset(gtf, V3 == "exon")
		write.table(data.frame(chrom=gtf[,'V1'], start=gtf[,'V4'], end=gtf[,'V5']), "exome.bed", quote = F, sep="\t", col.names = F, row.names = F)
	RCODE

### Copy into new file
	awk '{print $1 "\t" ($2 - 1) "\t" $3}' exome.bed > exome.fixed.bed

### Generate interval list
	/usr/lib/gatk/gatk \
		BedToIntervalList \
		-I exome.fixed.bed \
		-O /mnt/iquityazurefileshare1/GTFs/Homo_sapiens.GRCh38.104.gtf.exons.interval_list \
		-SD /mnt/iquityazurefileshare1/GTFs/SNPs/Homo_sapiens.GRCh38.dna.primary_assembly.dict	

------
## STAR generate references 
### Prepare for STAR
	mkdir /mnt/iquityazurefileshare1/GTFs/SNPs/STAR2_5/

### Run STAR
	STAR \
		--runMode genomeGenerate \
		--genomeDir /mnt/iquityazurefileshare1/GTFs/SNPs/STAR2_5/ \
		--genomeFastaFiles /mnt/iquityazurefileshare1/GTFs/SNPs/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		--sjdbGTFfile /mnt/iquityazurefileshare1/GTFs/Homo_sapiens.GRCh38.104.gtf \
		--sjdbOverhang 100 \
		--runThreadN 6

### Zip results [Optional]
	tar -zcvf star-HUMAN-refs.tar.gz mkdir /mnt/iquityazurefileshare1/GTFs/SNPs/STAR2_5/

------

# ANALYZING A SAMPLE 
------
## STAR alignment ✅
	STAR \
		--genomeDir /mnt/iquityazurefileshare1/GTFs/SNPs/STAR2_5/ \
		--runThreadN 6 \
		--readFilesIn /mnt/iquityazurefileshare1/IQuity/GATK_Testing/2343-TA-1_R1_sequence.fastq /mnt/iquityazurefileshare1/IQuity/GATK_Testing/2343-TA-1_R2_sequence.fastq \
		--sjdbOverhang 100 \
		--outSAMtype BAM SortedByCoordinate \
		--twopassMode Basic \
		--outFileNamePrefix "/mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing_2/2343-TA-1_"

------
## Mark duplicates ✅
	/usr/lib/gatk/gatk \
		MarkDuplicates \
		--INPUT /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_Aligned.sortedByCoord.out.bam \
		--OUTPUT /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_Aligned.dedupped.bam  \
		--CREATE_INDEX true \
		--VALIDATION_STRINGENCY SILENT \
		--METRICS_FILE /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_DUP.metrics		

JG time estimate: <5min
------
## Run SplitNCigar ✅
	/usr/lib/gatk/gatk \
		SplitNCigarReads \
		-R /mnt/iquityazurefileshare1/GTFs/SNPs/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		-I /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_Aligned.dedupped.bam  \
		-O /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_Aligned.dedupped.split.bam
	
JG time estimate: 5-10min
------
## Label read1 groups ✅
	java -jar /mnt/iquityazurefileshare1/IQuity/GATK_Testing/gatk-workflows/picard.jar \
		AddOrReplaceReadGroups \
		I= /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_Aligned.dedupped.split.bam \
		O= /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_Aligned.dedupped.split.labeled.bam \
		RGPL=illumina \
		RGPU=unit1 \
		RGLB=lib1 \
		RGSM=2343-TA-1 \
		RGID=1	

JG time estimate: 5-10min
------				
## Run BaseRecalibrator ✅
	/usr/lib/gatk/gatk \
		BaseRecalibrator \
		-R /mnt/iquityazurefileshare1/GTFs/SNPs/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		-I /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_Aligned.dedupped.split.labeled.bam \
		--use-original-qualities \
		-O /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_recal_data.csv \
		-known-sites /mnt/iquityazurefileshare1/GTFs/SNPs/homo_sapiens_clinically_associated.vcf \
		-known-sites /mnt/iquityazurefileshare1/GTFs/SNPs/homo_sapiens_phenotype_associated.vcf \
		-known-sites /mnt/iquityazurefileshare1/GTFs/SNPs/homo_sapiens_somatic_incl_consequences.vcf \
		-known-sites /mnt/iquityazurefileshare1/GTFs/SNPs/homo_sapiens_structural_variations.vcf


	## old files from GATK/Broad references
		-known-sites /mnt/iquityazurefileshare1/IQuity/GATK_Testing/gatk-workflows/inputs/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf \
    -known-sites /mnt/iquityazurefileshare1/IQuity/GATK_Testing/gatk-workflows/inputs/Homo_sapiens_assembly19_1000genomes_decoy/Mills_and_1000G_gold_standard.indels.b37.sites.vcf \
		-known-sites /mnt/iquityazurefileshare1/IQuity/GATK_Testing/gatk-workflows/inputs/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.known_indels.vcf \

JG time estimate: 5-10min
------
## Run BQSR ✅
	/usr/lib/gatk/gatk \
    ApplyBQSR \
    --add-output-sam-program-record \
    -R /mnt/iquityazurefileshare1/GTFs/SNPs/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -I /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_Aligned.dedupped.split.labeled.bam \
    --use-original-qualities \
    -O /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_Aligned.dedupped.split.labeled.recal.bam  \
    --bqsr-recal-file /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_recal_data.csv

JG time estimate: <5min
------
## HaplotypeCaller ✅
	/usr/lib/gatk/gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
		HaplotypeCaller \
		-R /mnt/iquityazurefileshare1/GTFs/SNPs/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		-I /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1_Aligned.dedupped.split.labeled.recal.bam \
		-L /mnt/iquityazurefileshare1/GTFs/Homo_sapiens.GRCh38.104.gtf.exons.interval_list \
		-O /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1.haplo.vcf.gz \
		-dont-use-soft-clipped-bases \
		--standard-min-confidence-threshold-for-calling 20 

	## optional --dbsnp flag
		--dbsnp /mnt/iquityazurefileshare1/IQuity/GATK_Testing/gatk-workflows/inputs/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf	

JG time estimate: 5-10min?
------
## VariantFiltration ✅
	/usr/lib/gatk/gatk \
		VariantFiltration \
		--R /mnt/iquityazurefileshare1/GTFs/SNPs/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		--V /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1.haplo.vcf.gz \
		--window 35 \
		--cluster 3 \
		--filter-name "FS" \
		--filter "FS>30.0" \
		--filter-name "QD" \
		--filter "QD<2.0" \
		-O /mnt/iquityazurefileshare1/IQuity/GATK_Testing/210830_Testing/2343-TA-1.haplo.variant_filtered.vcf.gz
		
**JG time estimate: 5-10min?**

