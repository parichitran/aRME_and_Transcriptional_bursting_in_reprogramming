library(SCALE)


#Get input of both 129 & CAST alleles

alleleA=read.table('129S1_filtered_grt_10_f25cells_exp.tsv',sep = "\t",header=TRUE)
alleleB=read.table('CAST_filtered_grt_10_f25cells_exp.tsv', sep="\t",header=TRUE)

alleleA=data.frame(alleleA, row.names="ID")
alleleB=data.frame(alleleB, row.names="ID")

genename <- rownames(alleleA)
sampname <- colnames(alleleA)

N = apply(alleleA + alleleB, 2, sum)
lib.size = N/mean(N)


#spikein_read=spike_in
#remove cells with extreme library size factors

lib.size.filter=(lib.size>0.85 & lib.size < 1.3)
alleleA=alleleA[,lib.size.filter]; alleleB=alleleB[,lib.size.filter]
dim(alleleA)

#Number of cells are stored for unfiltered analysis
write.table(alleleA, file = 'cellforunfiltered', col.names = TRUE,row.names = FALSE, quote = FALSE, sep = '\t')


abkt = c(0, 1, 10000, 0)

# Gene classification
gene.class.obj <- gene_classify(alleleA = as.matrix(alleleA),alleleB = as.matrix(alleleB))
A.prop = gene.class.obj$A.prop
B.prop = gene.class.obj$B.prop
gene.category = gene.class.obj$gene.category
results.list = gene.class.obj$results.list


# Estimate bursting kinetics parameters
cellsize=rep(1, ncol(alleleA))

allelic.kinetics.obj = allelic_kinetics(alleleA = alleleA,alleleB = alleleB,abkt = abkt,gene.category = gene.category,cellsize = cellsize, pdf = TRUE)

bandwidth = allelic.kinetics.obj$bandwidth
konA = allelic.kinetics.obj$konA; konB = allelic.kinetics.obj$konB
koffA = allelic.kinetics.obj$koffA; koffB = allelic.kinetics.obj$koffB
sA = allelic.kinetics.obj$sA; sB = allelic.kinetics.obj$sB
sizeA = sA/koffA; sizeB = sB/koffB


# Nonparametric test on whether the two alleles share the same burst frequency and burst size.
diff.allelic.obj = diff_allelic_bursting(alleleA = alleleA,alleleB = alleleB,cellsize = cellsize,gene.category = gene.category,abkt = abkt,allelic.kinetics.obj = allelic.kinetics.obj,mode = 'corrected')
pval.kon = diff.allelic.obj$pval.kon; pval.size = diff.allelic.obj$pval.size
alleleA = as.matrix(alleleA)
alleleB = as.matrix(alleleB)


# Chi-square test on whether the two alleles fire independently.

non.ind.obj = non_ind_bursting(alleleA = alleleA, alleleB = alleleB,gene.category = gene.category,results.list = results.list)
pval.ind = non.ind.obj$pval.ind; non.ind.type = non.ind.obj$non.ind.type
#i=which(genename=='Btf3l4')
#j=which(genename=='Xist')
#allelic_plot(alleleA = alleleA, alleleB = alleleB,gene.class.obj = gene.class.obj,allelic.kinetics.obj = allelic.kinetics.obj,diff.allelic.obj = diff.allelic.obj,non.ind.obj = non.ind.obj, i= i)
source("/media/gayenlab/GayenLab/software/SCALE-master_edited/SCALE-master/package/R/output_table.R")
SCALE.output=output_table(alleleA=alleleA, alleleB=alleleB,gene.class.obj = gene.class.obj,allelic.kinetics.obj = allelic.kinetics.obj,diff.allelic.obj = diff.allelic.obj,non.ind.obj = non.ind.obj)

#allelic_plot(alleleA = alleleA, alleleB = alleleB,gene.class.obj = gene.class.obj,allelic.kinetics.obj = allelic.kinetics.obj,diff.allelic.obj = diff.allelic.obj,non.ind.obj = non.ind.obj, i= j)
write.table(SCALE.output, file = 'SCALE.output.txt', col.names = TRUE,row.names = FALSE, quote = FALSE, sep = '\t')


