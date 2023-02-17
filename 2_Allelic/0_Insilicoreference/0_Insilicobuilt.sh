#! /bin/bash
#Downloading reference genome on which insilico genomes will be built

#
curl -C  - http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz -o Mus_musculus.GRCm38.102.gtf.gz

#
curl -C  - http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -o Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

#VCF file


#Making insilico allelic reference genomes
python3 prepare_reference.py --PSEUDOREF True --HETVCF True \
  --pseudoref_dir /full/path/to/dir/for/pseudo/ref/out/ \
  --vcf_dir /full/path/todir/for/vcf/outputs/ \
  --ref /full/path/to/ref/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
  --name_mat 129S1_SvImJ --name_pat CAST_EiJ \
  --vcf_joint /full/path/to/vcf/.. \
  --gtf /full/path/to/gtf/Mus_musculus.GRCm38.102.gtf



