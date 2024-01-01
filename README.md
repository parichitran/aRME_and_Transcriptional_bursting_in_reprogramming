# Autosomal Random Monoallelic Expressions and its association with Transcriptiona bursting during reprogramming.

Data is downloaded from GSE153847 using the download list from SRA(672files)




 
Fastqc(0.11.9) & Multiqc (1.12)reports were generated



 
                        
insilico reference genome for 129 and CAST were generated on mm10(GRC38.102) 






Mapped(STAR-2.7.10a) the reads to its corresponding allelic origin and  respective count files eventually



To exclude any false positive, we only considered those genes with at least 2 informative SNPs (at least 3 reads per SNP site).  we took an average of SNP wise reads to have the allelic read counts and biased SNP too removed???.And cells with atleast 200 reads(alleleA+alleleB>=200) are considered.



RPKM Normalisation is done.Consider only those genes with avg10RPKM and atleast 25 cells expression across cells at each daypoint.n=611 cells.Again cells with extreme library size( >1.3&<0.85)were removed




Performed genome-wide allele-specific burst kinetics analysis using 
SCALE (Jiang et al., 2017).




Gprofiler is used for getting GeneOntologies and Metascape for cell type profiling in various category of  genes


