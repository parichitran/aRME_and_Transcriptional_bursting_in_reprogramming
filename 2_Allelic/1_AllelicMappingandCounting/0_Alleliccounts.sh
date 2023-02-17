#! /bin/bash

#Step 0: Make index for both genomes
STAR --runMode genomeGenerate --genomeDir /pathto/129S1_SvImJ_VS_CAST/Ref_insilico/129S1_SvImJ/index --genomeFastaFiles 129S1_SvImJ_pseudo.fa --runThreadN 10 --limitGenomeGenerateRAM 33524399488

STAR --runMode genomeGenerate --genomeDir /pathto/129S1_SvImJ_VS_CAST/Ref_insilico/CAST_EiJ/index --genomeFastaFiles CAST_EiJ_pseudo.fa --runThreadN 10 --limitGenomeGenerateRAM 33524399488

#Move all the reads to this folder and run the scripts
mv /media/gayenlab/GayenLab10TB1/02_github/aRME_and_Transcriptional_bursting_in_reprogramming/0_Reads/*.fastq.gz /1_AllelicMappingandCounting
#Step 1: Getting RNA-seq data raw data 


#Create table files from gtf which will be used in downstream

awk '$3=="exon" && ($1 ~ /^[1-9]/ || $1 == "X" || $1 == "Y")' /full/path/to/gtf/Mus_musculus.GRCm38.102.gtf | cut -f1,4,5,9 | awk -F $'\t' 'BEGIN {OFS = FS} {split($4, a, "gene_id \""); split(a[2], b, "\""); print $1, $2-1, $3, b[1]}' > /full/path/to/output/Mus_musculus.GRCm38.102.EXONS.bed


grep "^#CHROM" /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.exons.vcf | cut -f1-5 > /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.snp_table.txt

grep -v "^#" /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.exons.vcf | sort -V | cut -f1-5 >> /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.snp_table.txt

date

#Process every reads in loop to get transcript counts at SNP positions
for f in `ls -1 *.fastq.gz | sed 's/.fastq.gz//' `
 do
 echo  ""${f}
 #1.Aligne fastq to 129S1
 STAR --readFilesIn ${f}.fastq.gz\
     --outFileNamePrefix ${f}_mat. \
     --runThreadN 11 --outSAMtype SAM \
     --outSAMattrRGline ID:mat \
     --readFilesCommand zcat \
     --genomeDir "/pathto/129S1_SvImJ_VS_CAST/Ref_insilico/129S1_SvImJ/index" \
     --outFilterMultimapNmax 1 --sjdbGTFfile /full/path/to/gtf/Mus_musculus.GRCm38.102.gtf

 wait
 #Aligne fastq to CAST
 STAR --readFilesIn ${f}.fastq.gz \
     --outFileNamePrefix ${f}_pat. \
     --runThreadN 11 --outSAMtype SAM \
     --outSAMattrRGline ID:pat \
     --readFilesCommand zcat \
     --genomeDir  "/pathto/129S1_SvImJ_VS_CAST/Ref_insilico/CAST_EiJ/index" \
     --outFilterMultimapNmax 1 --sjdbGTFfile /full/path/to/gtf/Mus_musculus.GRCm38.102.gtf

 #Sorting Bam
 wait
 samtools sort -n -o ${f}_mat.Nsorted.sam -@ 4 ${f}_mat.Aligned.out.sam
 wait
 samtools sort -n -o ${f}_pat.Nsorted.sam -@ 4 ${f}_pat.Aligned.out.sam


 wait
 #merge mat and pat
 python3 1_alleleseq_merge_stream_v2.py \
       --mat_sam ${f}_mat.Nsorted.sam \
       --pat_sam ${f}_pat.Nsorted.sam \
       --o ${f}_merged.sam \
       --paired 0
 wait
 #sort and convert to bam
 samtools sort -o ${f}_merged_SEreads.sorted.bam ${f}_merged.sam

 #read_count using allelcount
 wait
 python2 2_allelecounter.py --vcf "/full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.exons.vcf.gz" \
       --bam ${f}_merged_SEreads.sorted.bam \
       --ref "/full/path/129S1_SvImJ_pseudo.fa" \
       --sample F1 --min_cov 0 --min_baseq 20 --min_mapq 10 \
       --o ${f}_merged_SEreads.stat_0.txt

 wait
 # remove last file sam
 rm *.sam
 wait
 #final output gene* count file 
 Rscript --vanilla 3_counts_to_snp_genes.R \
        -d . \
        -n ${f}_merged_SEreads \
        -r ${f} \
        -o . \
        -v "/full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.snp_table.txt" \
        -b "/full/path/to/Mus_musculus.GRCm38.102.EXONS.bed"
done

date


#Only 2snps atleat 3 reads containing genes are taken with average of the reads

python3 4_129S1_and_CAST_removing_biased_SNPs_Alleiccounts.py
