library(data.table)
library(log4r)
library(RSQLite)
library(DBI)
library(ggplot2)
library(RSQLite)

options(width=155)
options(scipen=999) ## dont write number as scientific note

cfg = new.env()
cfg$assay = 'methylation'
cfg$species = 'GRCh37'
cwd = '/ifs/work/solitlab/huw/solit/study/hiseq/Project_06000_DD'
cfg$cwd = cwd

cfg$resDir = paste0(cwd, '/r_fang')
cfg$pre = 'scc'
cfg$assay = 'methylation'

delFiles = T

source('~/pipeline/configure_env.R')

setwd(cwd)
cwd

fread('target.tsv') -> target
target[, R12 := 'R1']
target[grep("R2_001", fastq), R12 := 'R2']
target[, samplename := sub(".*\\/(.*)_IGO.*", "\\1", fastq)]
target[, samplename := sub("_test", "", samplename)]
target[, samplename := paste0("pat-", samplename)]
target

target = dcast(target, samplename ~ R12, value.var = 'fastq')
target[, patientname := substr(samplename, 1, 7)]
target


target[, sampledir := paste0(cfg$initDir, '/', samplename)]
target[, {system(paste('mkdir -p ', sampledir))}, by=1:nrow(target)]
dir.exists(target$sampledir)

## fastqc
fastqc.jobname = paste0("fastqc")
fastq = paste(c(target$R1, target$R2), collapse=" ")
fastqc.cmd = bsub.head(fastqc.jobname, mem=10, cpu=20)
fastqc.cmd = paste0(fastqc.cmd, " \"", cfg$fastqc, " -o ", cfg$fastqcDir, " -j ", cfg$java, " -f fastq -t 20 ", fastq, " \"")
fastqc.cmd

exe.jobs(fastqc.cmd, logger)

## make genome for bismark
bismark.db.jobname = 'bismark.db'
bismark.db.cmd = bsub.head(bismark.db.jobname, mem=20, cpu=10)
bismark.db.cmd = paste0(bismark.db.cmd, " \"", cfg$bismark, '/bismark_genome_preparation --path_to_bowtie ', cfg$bowtie2dir, ' --verbose /ifs/work/solitlab/huw/study/db/hsa/grch37_bismark \"')
bismark.db.cmd

exe.jobs(bismark.db.cmd, cfg$logger)

## align
target[, bismark.jobname := paste0(samplename, '.bismark')]
target[, bamfile := paste0(cfg$sampledir, '/', cfg$samplename, '_pe.bam')] 
target[, bismark.cmd := bsub.head(bismark.jobname, mem=60, cpu=30, postdone = 'bismark.db', We = '30:00')]  
target[, bismark.cmd := paste0(bismark.cmd, " \"", cfg$bismark, "/bismark --path_to_bowtie ", cfg$bowtie2dir, " --quiet  -o ", sampledir)]
target[, bismark.cmd := paste0(bismark.cmd, " --samtools_path ", cfg$samtoolsdir, " -B ", samplename, " --nucleotide_coverage ")]
target[, bismark.cmd := paste0(bismark.cmd, " ", cfg$bismarkGRCh37, " -1 ", R1, " -2 ", R2, " \"")]
target[1, bismark.cmd]
target$bamfile
colnames(target)

exe.jobs(target$bismark.cmd, cfg$logger)

## deduplicate
#target[, sorted.bam.file := paste0(sampledir, '/', samplename, '_pe_sorted.bam')] 
target[, bamfile := paste0(sampledir, '/', samplename, '_pe.bam')] 
target[, dedup.bam.file := paste0(sampledir, '/', samplename, '.deduplicated_pe.bam')] 
target[, dedup.jobname := paste0(samplename, '.dedup')]
target[, dedup.cmd := bsub.head(dedup.jobname, mem=30, cpu=10, postdone = bismark.jobname, We = '30:00'), by=1:nrow(target)]  
target[, dedup.cmd := paste0(dedup.cmd, " \"", cfg$bismark, '/deduplicate_bismark -p --output_dir ', sampledir, ' --bam --samtools_path ', cfg$samtools, ' ', bamfile, "\"")]
target[2, dedup.cmd]

exe.jobs(target$dedup.cmd, cfg$logger)

## sort bam file
target[, sorted.bam.file := paste0(sampledir, '/', samplename, '_pe.deduplicated_sorted.bam')] 
target[, bam.sort.jobname := paste0(samplename, '.sort.bam')]
target[, bam.sort.cmd := bsub.head(bam.sort.jobname, mem=60, cpu=30, postdone = dedup.jobname, We = '30:00'), by=1:nrow(target)]  
target[, bam.sort.cmd := paste0(bam.sort.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx60g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' SortSam I=', dedup.bam.file, " O=", sorted.bam.file, " CREATE_INDEX=true SORT_ORDER=coordinate\"")]
target[2, bam.sort.cmd]

exe.jobs(target$bam.sort.cmd, cfg$logger)

## collect alignment summary metrics
target[, align.sum.file := paste0(sampledir, '/', samplename, '_align_sum.txt')]
target[, align.sum.jobname := paste0(samplename, '.align.sum')]
target[, align.sum.cmd := bsub.head(align.sum.jobname, mem=10, cpu=3, postdone = bam.sort.jobname, We = '10:00'), by = 1:nrow(target)]  
target[, align.sum.cmd := paste0(align.sum.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx10g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' CollectAlignmentSummaryMetrics R=', cfg$genomeGRCh37Fasta, " I=", sorted.bam.file, " O=", align.sum.file, " \"")]
target[2, align.sum.cmd]

exe.jobs(target$align.sum.cmd, cfg$logger)

## collect gc bias metrics
target[, gc.bias.metrics.file := paste0(sampledir, '/', samplename, '_gc_bias_metrics.txt')]
target[, gc.bias.metrics.pdf.file := paste0(sampledir, '/', samplename, '_gc_bias_metrics.pdf')]
target[, gc.bias.metrics.sum.file := paste0(sampledir, '/', samplename, '_gc_bias_metrics_sum.txt')]
target[, gc.bias.jobname := paste0(samplename, '.gc.bias')]
target[, gc.bias.cmd := bsub.head(gc.bias.jobname, mem=10, cpu=3, postdone = bam.sort.jobname, We = '10:00'), by = 1:nrow(target)]  
target[, gc.bias.cmd := paste0(gc.bias.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx10g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' CollectGcBiasMetrics R=', cfg$genomeGRCh37Fasta, " I=", sorted.bam.file, " O=", gc.bias.metrics.file, " CHART=", gc.bias.metrics.pdf.file, " S=", gc.bias.metrics.sum.file, "\"")]
target[1, gc.bias.cmd]

exe.jobs(target$gc.bias.cmd, cfg$logger)

## collect hs metrics
# mmuMethylInterval	= '/home/huw/program/truseq-methyl-capture-epic-manifest-file.bed'
target[, hs.metrics.file := paste0(sampledir, '/', samplename, '_hs_metrics.txt')]
target[, hs.jobname := paste0(samplename, '.hs')]
target[, hs.cmd := bsub.head(hs.jobname, mem=10, cpu=3, postdone = bam.sort.jobname, We = '10:00'), by = 1:nrow(target)]  
target[, hs.cmd := paste0(hs.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx10g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' CollectHsMetrics R=', cfg$genomeGRCh37Fasta, " I=", sorted.bam.file, " O=", hs.metrics.file, " TARGET_INTERVALS=", cfg$hsaMethylTargetInterval, " BAIT_INTERVALS=", cfg$hsaMethylBaitInterval, "\"")]
target[1, hs.cmd]
file.exists(cfg$hsaMethylTargetInterval)

exe.jobs(target$hs.cmd, cfg$logger)

## collect insert size
target[, insert.size.metrics.file := paste0(sampledir, '/', samplename, '_insert_size_metrics.txt')]
target[, insert.size.metrics.histgram.file := paste0(sampledir, '/', samplename, '_insert_size_metrics_histgram.pdf')]
target[, insert.size.jobname := paste0(samplename, '.insert.size')]
target[, insert.size.cmd := bsub.head(insert.size.jobname, mem=10, cpu=3, postdone = bam.sort.jobname, We = '10:00'), by = 1:nrow(target)]  
target[, insert.size.cmd := paste0(insert.size.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx10g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' CollectInsertSizeMetrics R=', cfg$genomeGRCh37Fasta, " I=", sorted.bam.file, " O=", insert.size.metrics.file, " H=", insert.size.metrics.histgram.file, " M=0.5\"")]
target[1, insert.size.cmd]

exe.jobs(target$insert.size.cmd, cfg$logger)

## collect targeted pcr metrics
target[, targeted.pcr.metrics.file := paste0(sampledir, '/', samplename, '_targeted.pcr_metrics.txt')]
target[, targeted.pcr.metrics.histgram.file := paste0(sampledir, '/', samplename, '_targeted.pcr_metrics_histgram.pdf')]
target[, targeted.pcr.jobname := paste0(samplename, '.targeted.pcr')]
target[, targeted.pcr.cmd := bsub.head(targeted.pcr.jobname, mem=10, cpu=3, postdone = bam.sort.jobname, We = '10:00'), by = 1:nrow(target)]  
target[, targeted.pcr.cmd := paste0(targeted.pcr.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx10g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' CollectInsertSizeMetrics R=', cfg$genomeGRCh37Fasta, " I=", sorted.bam.file, " O=", targeted.pcr.metrics.file, " H=", targeted.pcr.metrics.histgram.file, " M=0.5\"")]
target[1, targeted.pcr.cmd]

exe.jobs(target$targeted.pcr.cmd, cfg$logger)

## bismark_methylation_extractor
target[, bismark.ex.jobname := paste0(samplename, ".bismark.ex")]
target[, bismark.ex.cmd := bsub.head(bismark.ex.jobname, mem=30, We ='10:00', cpu=22, postdone = dedup.jobname), by=1:nrow(target)]  
target[, bismark.ex.cmd := paste0(bismark.ex.cmd, " \"", cfg$bismark, "/bismark_methylation_extractor --samtools_path ", cfg$samtoolsdir, " -o ", sampledir, " --parallel 5 --bedGraph")]
target[, bismark.ex.cmd := paste0(bismark.ex.cmd, " ", dedup.bam.file, " \"")]
target[1, bismark.ex.cmd]

exe.jobs(target$bismark.ex.cmd, cfg$logger)

## bismark2report
#target[, report.file := paste0(sampledir, '/', samplename, '_PE_report.txt')]
target[, dedup.report.file := paste0(sampledir, '/', samplename, '_pe.deduplicated_report.txt')]
target[, bismark.report.jobname := paste0(samplename, ".bismark.report")]
target[, bismark.report.cmd := bsub.head(bismark.report.jobname, mem=2, We ='10:00', cpu=1, postdone = bismark.jobname), by=1:nrow(target)]  
target[, bismark.report.cmd := paste0(bismark.report.cmd, " \"", cfg$bismark, "/bismark2report --dir ", sampledir, " --alignment_report ", dedup.report.file, "\"")]
target[1, bismark.report.cmd]

exe.jobs(target$bismark.report.cmd, cfg$logger)

## not run
## bismark2summary
target[, summary.file := paste0(sampledir, '/', samplename, '_PE_summary.txt')]
target[, bismark.summary.jobname := paste0(samplename, ".bismark.summary")]
target[, bismark.summary.cmd := bsub.head(bismark.summary.jobname, mem=5, We ='10:00', cpu=1, postdone = ), by=1:nrow(target)]  
target[, bismark.summary.cmd := paste0(bismark.summary.cmd, " \"", bismark, "/bismark2summary -o ", sampledir, "/", samplename, "_pe ", dedup.bam.file, "\"")]
target[1, bismark.summary.cmd]
target$dedup.bam.file
exe.jobs(target$bismark.summary.cmd, cfg$logger)
target$dedup.bam.file

## methylkit
library(methylKit)

target.interval = fread(cfg$hsaMethylBaitIntervalOrig, autostart=48)
setkey(target.interval, V1, V2, V3)

save.image()

## convert the bismark results for methylkit
## min coverage > 10
## methylkit.file for methylkit library
## methylkit.bed.file for calling C/G for strand information
## methylkit.base.file the actual C/G base data 
## revise the methylkit.file to add the strand information

target[, cov.file := paste0(sampledir, '/', samplename, '_pe.deduplicated.bismark.cov.gz')]
target[, methylkit.file := paste0(sampledir, '/', samplename, '_pe.deduplicated.bismark_methylkit.mincov10.txt')]
file.exists(target$cov.file)

target[, {
	met.cov = fread(paste0("zcat ", cov.file))
	met.cov = met.cov[(V5 + V6) > 10,] 
	met.cov[, V1 := paste0('chr', V1)]
	met.cov = foverlaps(met.cov, target.interval[, 1:3])
	head(met.cov)
	met.cov = met.cov[!is.na(V2),]
	met.cov[, coverage := V5 + V6]
	met.cov[, freqT := round(100 * V5 / coverage, 2)]
	met.cov[, freqC := round(100 * V6 / coverage, 2)]
	met.cov[, chrBase := paste0(V1, '.', i.V2)]
	met.cov[, strand := 'F']
	setnames(met.cov, 'i.V2', 'base')
	setnames(met.cov, 'V1', 'chr')
	met.cov = met.cov[, .(chrBase, chr, base, strand, coverage, freqC, freqT)]
	fwrite(met.cov, file=methylkit.file, sep="\t", quote=F)
}, by=1:nrow(target)]

## decide the strand for bismark coverage results
target$samplename
target[, methylkit.bed.file := paste0(sampledir, '/', samplename, '_pe.deduplicated.bismark_methylkit.mincov10.bed')]
target[, {
	cmd = paste0("awk '{if(NR==1){}else{print($2, $3-1, $3)}}' ", methylkit.file, " | sed \"s/ /\t/g\" | sort -k1 -k2 | sed \"s/chr//\" > ", methylkit.bed.file)
	cmd
	cat(samplename, " now ... \n\t", cmd)
	system(cmd)
}, by=1:nrow(target)]

target[, getCG.jobname := paste0(samplename, '.getCG')]
target[, base.gc.file := paste0(sampledir, '/', samplename, '_base_GC.txt')]
target[, getCG.cmd := bsub.head(getCG.jobname, mem=20, cpu=1, We='10:00')]
target[, getCG.cmd := paste0(getCG.cmd, " \"", cfg$bedtools, '/bedtools getfasta -fi ', cfg$genomeGRCh37Fasta, ' -bed ', methylkit.bed.file, ' -fo ', base.gc.file, ' -tab ', "\"")]
target[1, getCG.cmd]
exe.jobs(target$getCG.cmd, logger)

## revise the methylkit.file to add the strand information
target$samplename
target[, {
	#methylkit.file = target$methylkit.file[3]
	#base.gc.file = target$base.gc.file[3]
	#rm(methylkit.file, base.gc.file)
	cat("\n", samplename, ' now ...\n')
	met.cov = fread(methylkit.file, colClass=c(strand='character'))
	base.gc = fread(base.gc.file, header=F)
	base.gc[, chr := paste0('chr', sub(":.*", "", V1))]
	base.gc[, base :=  sub(".*-", "", V1)]
	base.gc[, tag := paste0(chr, '.', base)]
	merge(met.cov, base.gc[, .(tag, V2)], by.x = 'chrBase', by.y='tag') -> met.cov
	met.cov[, strand := 'X']
	met.cov[V2 == 'C', strand := 'F']
	met.cov[V2 == 'G', strand := 'R']
	head(met.cov)
	cat('\t total: ', nrow(met.cov), "\n")
	cat('\t strand F: ', nrow(met.cov[strand == 'F',]), "\n")
	cat('\t strand R: ', nrow(met.cov[strand == 'R',]), "\n")
	cat('\t strand X: ', nrow(met.cov[strand == 'X',]), "\n")
	if(nrow(met.cov[strand == 'X',]) > 0){
		cat(samplename, "\n\tsomething wrong here, row number ", which(met.cov$strand == 'X'))
	}
	met.cov[, V2 := NULL]
	fwrite(met.cov, file= methylkit.file, quote=F, sep="\t")
}, by = 1:nrow(target)]

## not run, segfault, memory not mapped
## methylkit::processBismarkAln.R
library(methylKit)
target[, group := substr(samplename, 9, 10)]
save(target, file='target.RData')
save(cfg, file='cfg.RData')
try{
    methyl.o = processBismarkAln(location = as.list(target$sorted.bam.file), sample.id=as.list(target$samplename), 
				 assembly='GRCh37', save.folder=cfg$methylDir, save.context=c("CpG"),
				 read.context=c("CpG"), nolap = F, mincov=10, minqual=20, phred64=F,
				 treatment=target$group, save.db=T)
}

try(
    {methyl.o = processBismarkAln(location = target$sorted.bam.file[1], sample.id=target$samplename[1], 
				  assembly='GRCh37', save.folder=cfg$methylDir, save.context=c("CpG"),
				  read.context=c("CpG"), nolap = F, mincov=10, minqual=20, phred64=F) }
)

## use the methylkit.file 
save.image()
save(target, cfg, file='apr15.RData')
library(methylKit)
target$samplename
methyl.o = methRead(as.list(target$methylkit.file), sample.id=as.list(target$samplename), assembly='GRCh37', treatment=factor(target$group), context='CpG')
## filter
methyl.fil.o = filterByCoverage(methyl.o, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=100)
rm(methyl.o)
save(methyl.fil.o, file=paste0(cfg$methylDir, '/methyl_filtered_obj.RData'))

## import this raw data into the raw sqlite
library(RSQLite)
library(DBI)

con = dbConnect(drv=SQLite(), dbname=paste0(cfg$methylDir, './sqlite/methyl.raw.sqlite'))
target[, {dt.m = fread(methylkit.file)
       dt.m[, samplename := samplename]
       dbWriteTable(con, 'methylRaw', as.data.frame(dt.m), append=T)
}, by=1:nrow(target)]

dbSendStatement(con, 'create index iindex on methylRaw (samplename, chr, base)')

dbDisconnet(con)


for(i in 1:nrow(target)){
	pdf(file=paste0(cfg$methylDir, '/', target$samplename[i], '.pdf'), width=12, height=6)
	getMethylationStats	(methyl.fil.o[[i]],	plot=T,		both.strands=T)
	getCoverageStats	(methyl.fil.o[[i]],	plot=TRUE,	both.strands=T)
	dev.off()
}

## combine, destrand for both strand when analyzing CpG
meth = unite(methyl.fil.o, destrand=T, min.per.group=1L)
fwrite(meth, file=paste0(cfg$methylDir, '/methylation_matrix.tsv'), sep="\t", quote=F)
save(meth, file=paste0(cfg$methylDir, '/meth_united.RData'))

pdf(file=paste0(cfg$methylDir, '/corr_pca_cluster.pdf'))
getCorrelation(meth, plot=TRUE)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
PCASamples(meth, screeplot=TRUE)
PCASamples(meth)
dev.off()

myDiff = calculateDiffMeth(meth, mc.cores=5)
myDiff0p    = getMethylDiff(myDiff, difference=0, qvalue=0.01)

pdf(file=paste0(cfg$methylDir, '/diff_methyl.pdf'))
diffMethPerChr(myDiff, plot=T, qvalue.cutoff=0.01, meth.cutoff=25)
dev.off()

library(genomation)
tx.file = system.file('extdata', 'hg19_gencode_genesymbol.bed.txt', package = 'methylKit')
gene.obj=readTranscriptFeatures(tx.file, up.flank=2000, down.flank=2000)
diffAnn=annotateWithGeneParts(as(myDiff0p, "GRanges"), gene.obj)
diffAnn
head(getAssociationWithTSS(diffAnn))
getTargetAnnotationStats(diffAnn, percentage=TRUE, precedence=TRUE)

cpg.file = system.file('extdata', 'hg19_cpg_island_tag_genesymbol.bed.txt', package = 'methylKit')
cpg.obj = readFeatureFlank(cpg.file, feature.flank.name=c("CpGi","shores"))
# convert methylDiff object to GRanges and annotate
diffCpGann = annotateWithFeatureFlank(as(myDiff0p,"GRanges"), cpg.obj$CpGi,cpg.obj$shores, feature.name = "CpGi",flank.name = "shores")
diffCpGann

pdf(file=paste0(cfg$methylDir, '/dmr_annotation.pdf'))
plotTargetAnnotation(diffAnn,precedence=TRUE, main="DMR")
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"), main="DMR by CpG island")
dev.off()

getFeatsWithTargetsStats(diffAnn,percentage=TRUE)
getFeatsWithTargetsStats(diffCpGann,percentage=TRUE)

## regional DMR promoter
gene.obj$promoters
promoters = regionCounts(methyl.fil.o, gene.obj$promoters)
promoters

promoter.meth = unite(promoters, destrand=T, min.per.group=1L)
promoter.meth
fwrite(promoter.meth, file=paste0(cfg$methylDir, '/methylation_promoter_region_matrix.tsv'), sep="\t", quote=F)
save(promoter.meth, file=paste0(cfg$methylDir, '/promoter_region_meth_united.RData'))

pdf(file=paste0(cfg$methylDir, '/promoter_region_corr_pca_cluster.pdf'))
getCorrelation(promoter.meth, plot=TRUE)
clusterSamples(promoter.meth, dist="correlation", method="ward", plot=TRUE)
PCASamples(promoter.meth, screeplot=TRUE)
PCASamples(promoter.meth)
dev.off()

myDiff.promoter.region = calculateDiffMeth(promoter.meth, mc.cores=3)
myDiff0p.promoter    = getMethylDiff(myDiff.promoter.region, difference=0, qvalue=0.01)

rm(myDiff25p.promoter.hi,myDiff25p.promoter.lo)

## regional DMR CpG island
cpgi.region = regionCounts(methyl.fil.o, cpg.obj$CpGi)
cpgi.shore.region = regionCounts(methyl.fil.o, cpg.obj$shores)

cpgi.region.unite  = unite(cpgi.region, destrand=T, min.per.group=1L)
cpgi.shore.unite   = unite(cpgi.shore.region, destrand=T, min.per.group=1L)
fwrite(cpgi.region.unite, file=paste0(cfg$methylDir, '/methylation_cpgi_region_unite_matrix.tsv'), sep="\t", quote=F)
fwrite(cpgi.shore.unite, file=paste0(cfg$methylDir, '/methylation_cpgi_shore_unite_matrix.tsv'), sep="\t", quote=F)
save(cpgi.region.unite, file=paste0(cfg$methylDir, '/methylation_cpgi_region_matrix.RData'))
save(cpgi.shore.unite, file=paste0(cfg$methylDir, '/methylation_cpgi_shore_matrix.RData'))

pdf(file=paste0(cfg$methylDir, '/cpgi_region_corr_pca_cluster.pdf'))
getCorrelation(cpgi.region.unite, plot=TRUE)
clusterSamples(cpgi.region.unite, dist="correlation", method="ward", plot=TRUE)
PCASamples(cpgi.region.unite, screeplot=TRUE)
PCASamples(cpgi.region.unite)
dev.off()

pdf(file=paste0(cfg$methylDir, '/cpgi_shore_region_corr_pca_cluster.pdf'))
getCorrelation(cpgi.shore.unite, plot=TRUE)
clusterSamples(cpgi.shore.unite, dist="correlation", method="ward", plot=TRUE)
PCASamples(cpgi.shore.unite, screeplot=TRUE)
PCASamples(cpgi.shore.unite)
dev.off()

myDiff.cpgi.region =  calculateDiffMeth(cpgi.region.unite, mc.cores=5)
myDiff.cpgi.shore  =  calculateDiffMeth(cpgi.shore.unite,  mc.cores=5)
myDiff0p.cpgi.region   = getMethylDiff(myDiff.cpgi.region, difference=0, qvalue=0.01)
myDiff0p.cpgi.shore    = getMethylDiff(myDiff.cpgi.shore, difference=0, qvalue=0.01)


## annotation all the variables
## gene regions:
tx.file = system.file('extdata', 'hg19_gencode_genesymbol.bed.txt', package = 'methylKit')
tx = fread(tx.file)
tx2 = tx[, 1:4]
colnames(tx2) = paste0('gene.', c('chr', 'start', 'end', 'tag'))
tx2
rm(tx)


## myDiff25p: annotate by single CpG DMR of gene parts
diffAnn = annotateWithGeneParts(as(myDiff0p, "GRanges"), gene.obj)
diffanno.tss = as.data.table(diffAnn@dist.to.TSS)
xx = as.data.table(as.data.frame(myDiff0p@.Data, col.names=colnames(myDiff0p)))
diffanno.tss = cbind(diffanno.tss, xx[diffanno.tss$target.row,]) 
diffanno.tss = merge(diffanno.tss, tx2, by.x = 'feature.name', by.y = 'gene.tag', all.x = T, all.y = F)
diffanno.tss = diffanno.tss[!is.na(meth.diff),]
diffanno.tss

getTargetAnnotationStats(diffAnn, percentage=TRUE, precedence=TRUE)

## annotate by DMR of CpG locations
diffCpGann = annotateWithFeatureFlank(as(myDiff0p,"GRanges"), cpg.obj$CpGi, cpg.obj$shores, feature.name = "CpGi", flank.name = "shores")
getFeatsWithTargetsStats(diffCpGann,percentage=TRUE)
diffCpGanno.v2 = cbind(xx, diffCpGann@members)
diffCpGanno.v2

## annotate myDiff25p.promoter
xx = as.data.table(as.data.frame(myDiff.promoter.region@.Data, col.names=colnames(myDiff.promoter.region)))
xx[, tag := paste0(chr, ':', start)]
tmp = tx2[, tag := paste0(gene.chr, ':', gene.start)]
myDiff.promoter = merge(xx, tmp, by.x='tag', by.y='tag')
dim(myDiff.promoter[pvalue < 0.01 & qvalue < 0.01,])

##  annotation of cpgi region
cpg.file = system.file('extdata', 'hg19_cpg_island_tag_genesymbol.bed.txt', package = 'methylKit')
cpg.file
cpg = fread(cpg.file)
head(cpg)
colnames(cpg) = paste0('cpg.', c('chr', 'start', 'end', 'tag'))
cpg[, tag := paste0(cpg.chr, ':', cpg.start)]

xx = as.data.table(as.data.frame(myDiff.cpgi.region@.Data, col.names=colnames(myDiff.cpgi.region)))
xx[, type := 'cpgi']
xx2 = as.data.table(as.data.frame(myDiff.cpgi.shore@.Data, col.names=colnames(myDiff.cpgi.shore)))
xx2[, type := 'shore']
rbind(xx, xx2) -> xx3
cpgi.region.shore = xx3
rm(xx, xx2, xx3)

xx = annotateWithGeneParts(as(cpgi.region.shore, "GRanges"), gene.obj)
as.data.table(xx@dist.to.TSS) -> xx.tss
cpgi.region.shore.anno = cbind(cpgi.region.shore[xx.tss$target.row,], xx.tss)
cpgi.region.shore.anno[pvalue < 0.01 & qvalue < 0.01,]
table(sign(cpgi.region.shore.anno[pvalue < 0.01 & qvalue < 0.01, meth.diff]))

diffanno.tss # single CpG associated with TSS
diffCpGanno.v2 # single CpG associated with CpG island or shore
myDiff.promoter # promoter region count, DMR, and associated with gene symbol
cpgi.region.shore.anno # cpgi region count, DMR, and 

options(scipen = 10)
diffanno.tss[grep('GATA3', feature.name),]

fwrite(diffanno.tss, file=paste0(cfg$methylDir, '/diff_cpg_genepart_anno_tss.tsv'), sep="\t", quote=F)
fwrite(diffCpGanno.v2, file=paste0(cfg$methylDir, '/diff_cpg_cpg_shore_anno.tsv'), sep="\t", quote=F)
fwrite(myDiff.promoter, file=paste0(cfg$methylDir, '/diff_promoter_region_anno.tsv'), sep="\t", quote=F)
fwrite(cpgi.region.shore.anno, file=paste0(cfg$methylDir, '/diff_cpg_island_region_shore_anno.tsv'), sep="\t", quote=F)

# diffanno.tss.gr = as(diffanno.tss, "GRanges")
# diffCpGanno.gr  = as(diffCpGanno.v2, "GRanges")
# myDiff.promoter.gr = as(myDiff.promoter, "GRanges")
# cpgi.region.shore.anno.gr  = as(cpgi.region.shore.anno, "GRanges")

tx = fread(tx.file)
tx3 = tx[, 1:8]; tx3
tx3[, symbol := sub(".*_", "", V4)] 
colnames(tx3) = c('chr', 'start', 'end', 'tag', 'score', 'strand', 'cds.start', 'cds.end', 'symbol')
tx3[symbol != 'NA' & symbol != '',]
tx.gr = as(tx3, 'GRanges')
tx.gr
rm(tx3, tx)

save(tx.gr, diffanno.tss, diffCpGanno.v2, myDiff.promoter, cpgi.region.shore.anno, file='last4.Rdata') ## 

# prepare for plotting specific region
con = dbConnect(drv=SQLite(), dbname=paste0(cfg$methylDir, '/sqlite/methyl.raw.sqlite'))

head(cpg)
plot.cpg = function(gene){
	g.cpgi = cpg[grep(gene, cpg.tag),]
	if(nrow(g.cpgi) < 1) { return(paste0('no cpgi found for gene ', gene))
	for(i in 
	g.bmin = g.cpgi$cpg.start - 1
	g.bmax = g.cpgi$cpg.end + 1
	g.chr  = g.cpgi$cpg.chr
	g.bmin;g.bmax;g.chr
	g.cpg.sel = dbGetQuery(con, paste0('select * from methylRaw where coverage > 10 and  chr="', g.chr, '" and base > ', g.bmin, ' and base < ', g.bmax))
	g.cpg.sel = as.data.table(g.cpg.sel)
	g.cpg.sel[, t12:= 'T1']
	g.cpg.sel[grep('T2', samplename), t12:= 'T2']
	g.cpg.sel[, tag := paste0(chr, ':', base)]
	g.cpg.sel
	sig.cpg = diffanno.tss[chr == g.chr & start > g.bmin & end < g.bmax & qvalue < 0.01 & abs(meth.diff) > 25, tag]
	sig.cpg
	g.cpg.sel[tag %in% sig.cpg,]
	gg = ggplot(g.cpg.sel, aes(x = base, y = samplename, color = freqC)) +
		geom_point(size=3) + xlab('') + scale_colour_gradient(high = "#132B43", low = "#56B1F7")
	ggsave(gg, file=paste0(cfg$methylDir, '/gene_gata3_CpG509.pdf'), width=29, height=2)
}

sync('r_fang/methyl')

## ma plot
tmp = diffanno.tss[abs(dist.to.feature) < 2000, ]
tmp[, Significance := 'N.S.']
tmp[qvalue < 0.001 & abs(meth.diff) > 50, Significance := 'Sig.']
table(tmp$Significance)
tmp
gg = ggplot(tmp, aes(x = meth.diff, y = -log(pvalue) + 3, color = Significance)) + 
       geom_point(size=1) + xlab("Methylation level (T2 vs T1)") + ylab("-log10(pvalue)") 
ggsave(gg, file=paste0(cfg$methylDir, '/maplot.png'), width=4, height=3)
sync('r_fang/methyl')



dbDisconnet(con)

ls()


xxa = xx[1:10,]
xxa
xxa[, sampleid := sub(".*(\\d+)$", "\\1", variable)]
value.var = colnames(meth.dt)[5:ncol(meth.dt)]
merge(meth.dt, value.var = colnames(meth.dt)[5:ncol(meth.dt)])
dcast(meth.dt)

## convert obj to db
cfg$methylDir
objDB = makeMethylDB(methyl.fil.o, 'r_fang/methyl/exMethDB')

library(GenomicRanges)
my.win=GRanges(seqnames="chr21", ranges=IRanges(start=seq(from=9764513,by=10000,length.out=20),width=5000) )
 
# selects the records that lie inside the regions
selectByOverlap(myobj[[1]], my.win)


perc.meth = percMethylation(meth)

library(methylKit)

pdf(file=paste0(cfg$methylDir, '/mix_model_cost.pdf'))
myMixdml = myDiff.to.mixmdl(myDiff, plot=T, main='SCC T2 vsT1')
plotCost(myMixmdl, main='cost function')
dev.off()

save.image()
## END
### https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html













## generate the alignment file
targets = fread(mapfile, header = F)
## group == library
targets = setNames(targets, c('group', 'sampleName', 'runID', 'fastqDir', 'seqType')) # 
targets.aln = copy(targets)
targets.aln[, fastq1 := '']
targets.aln[, fastq2 := '']
targets.aln[, fastq1.full := '']
targets.aln[, fastq2.full := '']
targets.aln = targets.aln[0,]
targets.aln

for(i in 1:nrow(targets)){
	tmp = copy(targets[i,])
	fastqs = dir(path=tmp$fastqDir, pattern="*.fastq.gz")
	fastq1 = fastqs[grep("_R1_", fastqs)]
	tmp = tmp[rep(1,length(fastq1)),]
	tmp[, fastq1 := fastq1]
	tmp[, fastq2 := sub("_R1_", "_R2_", fastq1)]
	tmp[, fastq1.full := paste0(fastqDir, '/', fastq1)]
	tmp[, fastq2.full := paste0(fastqDir, '/', fastq2)]
	if(! (all(file.exists(tmp$fastq1.full)) & all(file.exists(tmp$fastq2.full)))){
		cat("fastq file name problem:\n", file = logfile, append=T)
		write.table(tmp, logfile, append=T)
		break()
	}
	rbind(targets.aln, tmp) -> targets.aln
}
targets.aln

## align BWA
my $r1adaptor = 'AGATCGGAAGAGCACACGTCT';
my $r2adaptor = 'AGATCGGAAGAGCGTCGTGTA';
cwd = getwd()

## sample group runID
targets.aln[, sgrTag := paste0(sampleName, '_', group, '_', runID)]
targets.aln[, sgrDir := paste0(resDir, '/', sgrTag)]
targets.aln[, {system(paste0('mkdir -p ', sgrDir))}, by = 1:nrow(targets.aln)]

## zcat fastq
targets.aln[, fastq1.zcatfile := file.path(sgrDir, sub(".gz", "", fastq1))]
targets.aln[, fastq2.zcatfile := file.path(sgrDir, sub(".gz", "", fastq2))]

targets.aln[, zcat1.jobname := paste0("zcat.", sampleName, '.', fastq1)]
targets.aln[, zcat1.cmd := bsub.head(zcat1.jobname, mem=1, cpu=1, We='6:59', cwd)]
targets.aln[, zcat1.cmd:= paste0(zcat1.cmd, ' "/bin/zcat ', fastq1.full, ' > ', fastq1.zcatfile, '"')]

targets.aln[, zcat2.jobname := paste0("zcat.", sampleName, '.', fastq2)]
targets.aln[, zcat2.cmd := bsub.head(zcat2.jobname, mem=1, cpu=1, We='6:59', cwd)]
targets.aln[, zcat2.cmd := paste0(zcat2.cmd, ' "/bin/zcat ', fastq2.full, ' > ', fastq2.zcatfile, '"')]
targets.aln$zcat.cmd2[1]

exe.jobs(targets.aln$zcat1.cmd, cfg$logger)
exe.jobs(targets.aln$zcat2.cmd, cfg$logger)

## convert quality score
targets.aln[, fastq1.cqsfile := file.path(sgrDir, sub(".gz", ".cqs", fastq1))]
targets.aln[, fastq2.cqsfile := file.path(sgrDir, sub(".gz", ".cqs", fastq2))]
targets.aln$fastq1.cqsfile[1]

targets.aln[, cqs1.jobname := paste0("cqs1.", sampleName, '.', fastq1)]
targets.aln[, cqs1.cmd := bsub.head(cqs1.jobname, mem=1, cpu=1, We='2:59', cwd, zcat1.jobname), by=1:nrow(targets.aln)]
targets.aln[, cqs1.cmd := paste0(cqs1.cmd, ' "', ConvertQualityScore, ' --input ', fastq1.zcatfile, ' --output ', fastq1.cqsfile, '"')]
targets.aln$cqs1.cmd[1]

targets.aln[, cqs2.jobname := paste0("cqs2.", sampleName, '.', fastq1)]
targets.aln[, cqs2.cmd := bsub.head(cqs2.jobname, mem=1, cpu=1, We='2:59', cwd, zcat2.jobname), by = 1:nrow(targets.aln)]
targets.aln[, cqs2.cmd := paste0(cqs2.cmd, ' "', ConvertQualityScore, ' --input ', fastq2.zcatfile, ' --output ', fastq2.cqsfile, '"')]
targets.aln$cqs2.cmd[1]

exe.jobs(targets.aln$cqs2.cmd, cfg$logger)
exe.jobs(targets.aln$cqs2.cmd, cfg$logger)

## discard read less than half length
## find the read length
read.len = function(x){
	gzfile(x, open='r') -> fh
	readLines(fh, 2) -> h
	nchar(h[2]) }

targets.aln[, read.len := sapply(targets.aln$fastq1.full, read.len)]
targets.aln$read.len

targets.aln[, fastq1.cutadaptfile :=  file.path(sgrDir, sub(".gz", ".cutadapt", fastq1))]
targets.aln[, fastq2.cutadaptfile :=  file.path(sgrDir, sub(".gz", ".cutadapt", fastq2))]
targets.aln$fastq1.cutadaptfile[1]

targets.aln[, cutadapt.jobname := paste0("cutadapt.", sgrTag)]
targets.aln[, cutadapt.cmd := bsub.head(cutadapt.jobname, mem=1, cpu=1, We='2:59', cwd, postdone = paste(cqs1.jobname, cqs2.jobname, collapse=" ")), by = 1:nrow(targets.aln)]
targets.aln[, cutadapt.cmd := paste(cutadapt.cmd, ' "', PYTHON, CUTADAPT, '-f fastq -a', clipR1, '-A', clipR2, '--overlap 10 --minimum-length', floor(read.len/2))]
targets.aln[, cutadapt.cmd := paste(cutadapt.cmd, '--quality-cutoff', bqtrim, '-o', fastq1.cutadaptfile, '--paired-output', fastq2.cutadaptfile, fastq1.cqsfile, fastq2.cqsfile, ' "')]
targets.aln$cutadapt.cmd[1]

exe.jobs(targets.aln$cutadapt.cmd, logger)

## breakpoint
wait4jobs(targets.aln$cutadapt.jobname, logger)
testOutputFiles = c(targets.aln$fastq1.cutadapt, targets.aln$fastq2.cutadapt)

if(delFiles){
	del.files(targets.aln$fastq1.zcatfile, logger)
	del.files(targets.aln$fastq2.zcatfile, logger)
	del.files(targets.aln$fastq1.cqsfile,  logger)
	del.files(targets.aln$fastq2.cqsfile,  logger)
}

## alignment
## targets.aln[, sgrTag := paste0(sampleName, '_', group, '_', runID)]
targets.aln[, reference := B37_BWA_INDEX]

targets.aln[, readgroup := paste0("@RG\\\\tID:", sgrTag, '_', seqType, '\\\\tPL:Illumina\\\\tPU:', sgrTag, "\\\\tLB:", sampleName, '_', group, "\\\\tSM:", sampleName)]
targets.aln[, bamfile := paste0(sgrDir, '/', sgrTag, '_', fastq1, '.bam')] # 
targets.aln$bamfile[1]
file.exists(targets.aln$bamfile)
file.exists(targets.aln$fastq2.cutadaptfile)

targets.aln[, bwa.jobname := paste0("bwa.", sgrTag, fastq1)]
targets.aln[, bwa.cmd := bsub.head(bwa.jobname, mem=15, cpu=12, We='2:59', cwd)]
targets.aln[, bwa.cmd := paste0(bwa.cmd, ' "', bwa, ' mem -t 8 -T 0 -R \\\"', readgroup, '\\\" ', reference, ' ', fastq1.cutadaptfile, ' ', fastq2.cutadaptfile, ' | ', cfg$samtools, ' view -Shb -o ', bamfile, ' - "')]
targets.aln$bwa.cmd[1]

exe.jobs(targets.aln$bwa.cmd, logger)

## breakpoint
wait4jobs(targets.aln$bwa.jobname, logger)

## clean files
testInputFiles = c(targets.aln$fastq1.cutadapt, targets.aln$fastq2.cutadapt)
testOutputFiles = targets.aln$bamFile
if(delFiles & all(file.size(testOutputFiles) > 2000000)){
	del.files(targets.aln$fastq1.cutadapt, logger)
	del.files(targets.aln$fastq2.cutadapt, logger)
}

## merge bam file for the same library/group 
## may several r1/r2 fastq files
tmpdir='/scratch/huw/'
targets.bam = targets.aln[, .(bamFileList = paste(paste0(' I=', bamfile), collapse = ' '), bwaJobList = paste(bwa.jobname, collapse = ' ')), by = list(sampleName, group)]
targets.bam[, libMergedBamFile := paste0(resDir, '/', sampleName, '_library', group, '_merged.bam')]
targets.bam[, merge.jobname := paste0('mergebam.', sampleName, group)]
targets.bam[, merge.cmd := bsub.head(merge.jobname, mem=30, cpu=8, We='2:26', cwd=cwd, postdone = bwaJobList)]
targets.bam[, merge.cmd := paste0(merge.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, '/picard.jar MergeSamFiles ', bamFileList, ' O=', libMergedBamFile, ' SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=', tmpdir, ' CREATE_INDEX=TRUE USE_THREADING=FALSE MAX_RECORDS_IN_RAM=5000000 "')] 
targets.bam$merge.cmd[2]

exe.jobs(targets.bam$merge.cmd, logger)

## markDuplicates for reads from the same library/group
## merge bam file from the same library/group
## output in resDir (r_fang)
targets.bam[, mkdupedLibBamFile := paste0(resDir, '/', sampleName, '_library', group, '_merged_mkdup.bam')]
targets.bam[, libmkdupMetricsFile := paste0(matricsDir, '/', sampleName, '_library', group, '_merged_mkdup_metrics.txt')]
targets.bam[, mkdup.jobname := paste0('mkdup.', sampleName, '.', group)]
targets.bam[, mkdup.cmd := bsub.head(mkdup.jobname, mem=30, cpu=3, We='8:26', cwd=cwd, postdone = merge.jobname)]
targets.bam[, mkdup.cmd := paste0(mkdup.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', cfg$PICARD, '/picard.jar MarkDuplicates I=', libMergedBamFile, ' O=', mkdupedLibBamFile, ' METRICS_FILE=', libmkdupMetricsFile, ' VALIDATION_STRINGENCY=LENIENT TMP_DIR=', tmpdir, ' REMOVE_DUPLICATES=TRUE CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=5000000 "')] 
targets.bam$mkdup.cmd[3]

file.exists(targets.bam$libMergedBamFile)
exe.jobs(targets.bam$mkdup.cmd, logger)

## merge bam file for each sample
tmpdir='/scratch/huw/'
targets.smpBam = targets.bam[, .(BamFileList = paste(paste0(' I=', mkdupedLibBamFile), collapse = ' '),  mkdupJobList = paste(mkdup.jobname, collapse = ' '), bamFile = paste(mkdupedLibBamFile, collapse=' '), fileNumber = length(mkdupedLibBamFile)), by = list(sampleName)]
colnames(targets.smpBam)
targets.smpBam
targets.smpBam[, smpBamFile := paste0(resDir, '/', sampleName, '_merged_mkdup_smp.bam')]
targets.smpBam[, smp.jobname := paste0('smp.merge.', sampleName)]
targets.smpBam$smp.jobname
targets.smpBam[, smp.cmd := bsub.head(smp.jobname, mem=30, cpu=8, We='8:26', cwd=cwd, postdone = mkdupJobList)]
targets.smpBam[, smp.cmd := paste0(smp.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', cfg$PICARD, '/picard.jar MergeSamFiles ', BamFileList, ' TMP_DIR=', tmpdir, ' DIR=', tmpdir, ' O=', smpBamFile, ' SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=', tmpdir, ' CREATE_INDEX=TRUE USE_THREADING=FALSE MAX_RECORDS_IN_RAM=5000000 "')] 
targets.smpBam$smp.cmd[3] 

exe.jobs(targets.smpBam$smp.cmd, logger)

## breakpoint before collect metrics
wait4jobs(targets.smpBam$smp.jobname, logger)

## clean files
testInputFiles = targets.smpBam$bamFile
testOutputFiles = targets.smpBam$smpBamFile
if(delFiles & all(file.size(testOutputFiles) > 2000000)){
	del.files(targets.aln$fastq1.cutadapt, logger)
}
#else{stop("something might wrong here: the _merged_mkdup_smp.bam file for sample dont have the right size\n")}

##  calculate HS matrix  for the bam file for each sample
## the baits targets and the bam file
## variant_pipeline: 1402
tmpdir='/scratch/huw/'
targets.smpBam[, hmsFile := paste0(matricsDir, '/', 'hms_', sampleName, '_HsMatrix.txt')]
targets.smpBam[, hms.jobname := paste0('hms.', sampleName)]
targets.smpBam[, hms.cmd := bsub.head(hms.jobname, mem=10, cpu=1, We='8:26', cwd=cwd, postdone = smp.jobname)]
targets.smpBam[, hms.cmd := paste0(hms.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar CalculateHsMetrics I=', smpBamFile, ' TMP_DIR=', tmpdir, ' O=', hmsFile, ' REFERENCE_SEQUENCE=', genomeFasta, ' METRIC_ACCUMULATION_LEVEL=SAMPLE BAIT_INTERVALS=', baits_ilist, ' BAIT_SET_NAME=', wesImpactTarget, ' TARGET_INTERVALS=', targets_ilist, ' VALIDATION_STRINGENCY=LENIENT "')]
targets.smpBam$hms.cmd[3]

exe.jobs(targets.smpBam$hms.cmd, logger)

## for pair end, calculate the insert size
tmpdir='/scratch/huw/'
targets.smpBam[, insertSizeFile := paste0(matricsDir, '/', 'readlen_', sampleName, '_insertSizeFile.txt')]
targets.smpBam[, insertSizeHistgramFile := paste0(matricsDir, '/', 'readlen_', sampleName, '_insertSizeFile.txt')]
targets.smpBam[, readlen.jobname := paste0('readlen.', sampleName)]
targets.smpBam[, readlen.cmd := bsub.head(readlen.jobname, mem=10, cpu=1, We='8:26', cwd=cwd, postdone = smp.jobname)]
targets.smpBam[, readlen.cmd := bsub.head(readlen.jobname, mem=10, cpu=1, We='8:26', cwd=cwd)]
targets.smpBam[, readlen.cmd := paste0(readlen.cmd, ' "', JAVA, '/java -Xms256m -Xmx10g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar CollectInsertSizeMetrics I=', smpBamFile)]
targets.smpBam[, readlen.cmd := paste0(readlen.cmd, ' O=', insertSizeFile, ' REFERENCE_SEQUENCE=', genomeFasta, ' METRIC_ACCUMULATION_LEVEL=SAMPLE HISTOGRAM_FILE=', insertSizeHistgramFile)]
targets.smpBam[, readlen.cmd := paste0(readlen.cmd, ' VALIDATION_STRINGENCY=LENIENT TMP_DIR=', tmpdir, ' "')]
targets.smpBam$readlen.cmd[3]

exe.jobs(targets.smpBam$readlen.cmd, logger)

##  CollectAlignmentSummaryMetrics
tmpdir='/scratch/huw/'
targets.smpBam[, alnSummaryFile := paste0(matricsDir, '/', 'alnSum_', sampleName, '_alnSummaryFile.txt')]
targets.smpBam[, alnSum.jobname := paste0('alnSum.', sampleName)]
targets.smpBam[, alnSum.cmd := bsub.head(alnSum.jobname, mem=10, cpu=1, We='8:26', cwd=cwd, postdone = smp.jobname)]
#targets.smpBam[, alnSum.cmd := bsub.head(alnSum.jobname, mem=10, cpu=1, We='8:26', cwd=cwd)]
targets.smpBam[, alnSum.cmd := paste0(alnSum.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar  CollectAlignmentSummaryMetrics I=', smpBamFile)]
targets.smpBam[, alnSum.cmd := paste0(alnSum.cmd, ' O=', alnSummaryFile, ' REFERENCE_SEQUENCE=', genomeFasta, ' METRIC_ACCUMULATION_LEVEL=SAMPLE TMP_DIR=', tmpdir, ' "')]
targets.smpBam$alnSum.cmd[3]

exe.jobs(targets.smpBam$alnSum, logger)

## CollectOxoGMetrics
tmpdir='/scratch/huw/'
targets.smpBam[, OxoFile := paste0(matricsDir, '/', 'Oxo_', sampleName, '_OxoFile.txt')]
targets.smpBam[, Oxo.jobname := paste0('Oxo.', sampleName)]
targets.smpBam[, Oxo.cmd := bsub.head(Oxo.jobname, mem=4, cpu=1, We='8:26', cwd=cwd, postdone = smp.jobname)]
targets.smpBam[, Oxo.cmd := bsub.head(Oxo.jobname, mem=4, cpu=1, We='8:26', cwd=cwd)]
targets.smpBam[, Oxo.cmd := paste0(Oxo.cmd, ' "', JAVA, '/java -Xms256m -Xmx4g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar CollectOxoGMetrics  I=', smpBamFile)]
targets.smpBam[, Oxo.cmd := paste0(Oxo.cmd, ' REFERENCE_SEQUENCE=', genomeFasta, ' O=', OxoFile, ' DB_SNP=', DB_SNP, ' VALIDATION_STRINGENCY=LENIENT')];
targets.smpBam[, Oxo.cmd := paste0(Oxo.cmd,' TMP_DIR=', tmpdir, '  "')];
targets.smpBam$Oxo.cmd[3]

exe.jobs(targets.smpBam$Oxo.cmd, logger)

## depth of coverage
tmpdir='/scratch/huw/'
targets.smpBam[, depthCoverageFile := paste0(matricsDir, '/', 'depthCoverage_', sampleName, '_depthCoverageFile.txt')]
targets.smpBam[, depthCoverage.jobname := paste0('depthCoverage.', sampleName)]
targets.smpBam[, depthCoverage.cmd := bsub.head(depthCoverage.jobname, mem=4, cpu=1, We='8:26', cwd=cwd, postdone = smp.jobname)]
targets.smpBam[, depthCoverage.cmd := bsub.head(depthCoverage.jobname, mem=4, cpu=1, We='8:26', cwd=cwd)]
targets.smpBam[, depthCoverage.cmd := paste0(depthCoverage.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar -T DepthOfCoverage')]
targets.smpBam[, depthCoverage.cmd := paste0(depthCoverage.cmd, ' -I ', smpBamFile, ' -R ', genomeFasta, ' -o ', depthCoverageFile, ' -L ', FP_INT, ' -rf BadCigar')]
targets.smpBam[, depthCoverage.cmd := paste0(depthCoverage.cmd, ' -mmq 20 -mbq 0 -omitLocusTable -omitSampleSummary -baseCounts --includeRefNSites -omitIntervals "')]
targets.smpBam$depthCoverage.cmd[3]

exe.jobs(targets.smpBam$depthCoverage.cmd, logger)

## GC bias metrics
tmpdir = "/scratch/huw/"
targets.smpBam[, gcbiasFile := paste0(matricsDir, '/', 'gcbias_', sampleName, '_gcbiasFile.txt')]
targets.smpBam[, gcbiasSummaryFile := paste0(matricsDir, '/', 'gcbias_', sampleName, '_gcbiasSummaryFile.txt')]
targets.smpBam[, gcbiasChartFile := paste0(matricsDir, '/', 'gcbias_', sampleName, '_gcbiasChartFile.pdf')]
targets.smpBam[, gcbias.jobname := paste0('gcbias.', sampleName)]
targets.smpBam[, gcbias.cmd := bsub.head(gcbias.jobname, mem=4, cpu=1, We='8:26', cwd=cwd, postdone = smp.jobname)]
targets.smpBam[, gcbias.cmd := bsub.head(gcbias.jobname, mem=4, cpu=1, We='8:26', cwd=cwd)]
targets.smpBam[, gcbias.cmd := paste0(gcbias.cmd, ' "', JAVA, '/java -Xms256m -Xmx4g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar CollectGcBiasMetrics I=', smpBamFile)]
targets.smpBam[, gcbias.cmd := paste0(gcbias.cmd, ' O=', gcbiasFile, ' REFERENCE_SEQUENCE=', genomeFasta, ' SUMMARY_OUTPUT=', gcbiasSummaryFile)]
targets.smpBam[, gcbias.cmd := paste0(gcbias.cmd, ' CHART_OUTPUT=', gcbiasChartFile, ' TMP_DIR=', tmpdir, ' "')]
targets.smpBam$gcbias.cmd[3]

exe.jobs(targets.smpBam$gcbias.cmd, logger)

## summarize metrics file
## need add

## breakpoint after collect metrics
testJobname = c(targets.smpBam$gcbias.jobname, targets.smpBam$Oxo.jobname, targets.smpBam$depthCoverage.jobname, targets.smpBam$alnSum.jobname, targets.smpBam$readlen.jobname, targets.smpBam$hms.jobname)
wait4jobs(testJobname, logger)

## cleaning and realignment
## organize the targets
targets.smpGroup = fread(groupfile, header = F)
targets.smpGroup = setNames(targets.smpGroup, c('sampleName', 'smpGroup'))
targets.smpGroup
targets.smpGroup = merge(targets.smpBam[, .(sampleName, smpBamFile, smp.jobname)], targets.smpGroup, by = 'sampleName')
targets.smpGroup

targets.recal = targets.smpGroup[, .(bamFileList = paste(paste0('-I ', smpBamFile), collapse=' '), smpJobnameList = paste(smp.jobname, collapse = ', ')), by=smpGroup]
targets.recal

## REALIGNTARGETCREATOR 
tmpdir = "/scratch/huw/"
targets.recal[, realnGenTargetFile := paste0(matricsDir, '/', 'realn_',smpGroup, '_realnGenTargetFile.intervals')]
targets.recal[, realnGen.jobname := paste0('realnGen_', smpGroup)]
targets.recal[, realnGen.cmd := bsub.head(realnGen.jobname, mem=15, cpu=1, We='8:26', cwd=cwd, postdone = smpJobnameList)]
targets.recal[, realnGen.cmd := bsub.head(realnGen.jobname, mem=15, cpu=1, We='8:26', cwd=cwd)]
targets.recal[, realnGen.cmd := paste0(realnGen.cmd, ' "', JAVA, '/java -Xms256m -Xmx15g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK)]
targets.recal[, realnGen.cmd := paste0(realnGen.cmd, '/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ', genomeFasta, ' -known ', DB_SNP, ' -known ', MILLS_1000G, ' ', bamFileList)]
## keep the same the -known with the next step
targets.recal[, realnGen.cmd := paste0(realnGen.cmd, ' --out ', realnGenTargetFile, ' "')]
targets.recal$realnGen.cmd[1]

exe.jobs(targets.recal$realnGen.cmd, logger)
for(i in 1:nrow(targets.recal)){system(targets.recal$realnGen.cmd[i])}

targets.smpGroup$realnGenTargetFile = targets.recal[targets.smpGroup$smpGroup, realnGenTargetFile]
targets.smpGroup

## IndelRealigner
tmpdir = "/scratch/huw/"
targets.recal[, indelRealnOutFile := paste0(matricsDir, '/', 'realn_',smpGroup, '_indelRealn_OutFile.txt')]
targets.recal[, realnIndel.jobname := paste0('realnIndel.', smpGroup)] 
targets.recal[, realnIndel.cmd := bsub.head(realnIndel.jobname, mem=15, cpu=1, We='8:26', cwd=cwd, postdone = realnGen.jobname)]
targets.recal[, realnIndel.cmd := bsub.head(realnIndel.jobname, mem=15, cpu=1, We='8:26', cwd=cwd)]
targets.recal[, realnIndel.cmd := paste0(realnIndel.cmd, ' "', JAVA, '/java -Xms256m -Xmx15g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')]
targets.recal[, realnIndel.cmd := paste0(realnIndel.cmd, ' -T IndelRealigner -R ', genomeFasta, ' -known ', DB_SNP, ' -known ', MILLS_1000G, ' -targetIntervals ', realnGenTargetFile)]
targets.recal[, realnIndel.cmd := paste0(realnIndel.cmd, ' -L ', targets_bed)]
targets.recal[, realnIndel.cmd := paste0(realnIndel.cmd, ' --noOriginalAlignmentTags -S LENIENT --maxReadsForRealignment 500000 --maxReadsInMemory 3000000 --maxReadsForConsensuses 500000 -rf BadCigar ')]
targets.recal[, realnIndel.cmd := paste0(realnIndel.cmd, bamFileList, ' -nWayOut _realnIndel.bam "')]  ## filename change here from .bam to _realnindel.bam
# bamFileList => ' -I file1.bam -I file2.bam'
# the output bamfiles will be ' -I file1_realnIndel.bam -I file2_realnIndel.bam'
targets.recal$realnIndel.cmd[1]

#for(i in 1:nrow(targets.recal)){system(targets.recal$realnIndel.cmd[i])}
exe.jobs(targets.recal$realnIndel.cmd, logger)

# generated by the above realnIndel step nWayOut
# they are under the cwd !!!
targets.smpGroup[, smpRealnBamFile := sub(".bam", "_realnIndel.bam", smpBamFile)] 
targets.smpGroup[, smpRealnBaiFile := sub(".bam", "_realnIndel.bai", smpBamFile)] 
targets.smpGroup$realnIndel.jobname = targets.recal[targets.smpGroup$smpGroup, realnIndel.jobname]
targets.smpGroup

wait4jobs(targets.recal$realnIndel.jobname, logger)
for(i in 1:length(targets.smpGroup$smpRealnBamFile)){ system(paste0('mv ', targets.smpGroup$smpRealnBamFile[i], ' ', resDir))}
for(i in 1:length(targets.smpGroup$smpRealnBamFile)){ system(paste0('mv ', targets.smpGroup$smpRealnBaiFile[i], ' ', resDir))}

file.exists(targets.smpGroup$smpRealnBamFile)
file.exists(targets.smpGroup$smpBamFile)

## Realign use -known
## while recal use --knownSites

## baserecal 
targets.smpGroup[, baserecalOutFile := paste0(matricsDir, '/', 'realn_', sampleName, '_recal_OutFile.txt')]
targets.smpGroup[, baserecal.jobname := paste0('baserecal.', sampleName)]
#targets.smpGroup[, baserecal.cmd := bsub.head(baserecal.jobname, mem=40, cpu=12, We='8:26', cwd=cwd)] , postdone = realnIndel.jobname)]
targets.smpGroup[, baserecal.cmd := bsub.head(baserecal.jobname, mem=40, cpu=12, We='8:26', cwd=cwd)]
targets.smpGroup[, baserecal.cmd := paste0(baserecal.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')]
targets.smpGroup[, baserecal.cmd := paste0(baserecal.cmd, ' -T BaseRecalibrator -R ', genomeFasta, ' --knownSites ', DB_SNP, ' --knownSites ', MILLS_1000G)] 
targets.smpGroup[, baserecal.cmd := paste0(baserecal.cmd, ' --knownSites ', HAPMAP, ' --knownSites ', OMNI_1000G, ' --knownSites ', PHASE1_SNPS_1000G)]
targets.smpGroup[, baserecal.cmd := paste0(baserecal.cmd, ' --covariate ContextCovariate --covariate CycleCovariate --covariate QualityScoreCovariate --covariate ReadGroupCovariate -rf BadCigar --num_cpu_threads_per_data_thread 12 ')]
targets.smpGroup[, baserecal.cmd := paste0(baserecal.cmd, ' -S LENIENT -rf BadCigar --out ', baserecalOutFile, ' -I ', smpRealnBamFile, ' "')]
targets.smpGroup$baserecal.cmd[1]
targets.smpGroup$baserecal.cmd

exe.jobs(targets.smpGroup$baserecal.cmd, logger)

## printreads
targets.smpGroup[, printReadsBQSROutFile := paste0(matricsDir, '/', 'realn_', sampleName, '_bqsr_outfile.txt')]
targets.smpGroup[, printReadsBamOutFile  := sub('.bam', "_recal.bam", smpRealnBamFile)]
targets.smpGroup[, printReads.jobname := paste0('printReads.', sampleName)]
targets.smpGroup[, printReads.cmd := bsub.head(printReads.jobname, mem=30, cpu=6, We='8:26', cwd=cwd, postdone = baserecal.jobname)]
targets.smpGroup[, printReads.cmd := bsub.head(printReads.jobname, mem=30, cpu=6, We='8:26', cwd=cwd)]
targets.smpGroup[, printReads.cmd := paste0(printReads.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')]
targets.smpGroup[, printReads.cmd := paste0(printReads.cmd, ' -T PrintReads -R ', genomeFasta, ' --emit_original_quals -BQSR ', baserecalOutFile)]
targets.smpGroup[, printReads.cmd := paste0(printReads.cmd, ' --num_cpu_threads_per_data_thread 6 -rf BadCigar --downsampling_type NONE' )]
targets.smpGroup[, printReads.cmd := paste0(printReads.cmd, ' --out ', printReadsBamOutFile, ' -I ', smpRealnBamFile, ' "')]
targets.smpGroup$printReads.cmd[1]

exe.jobs(targets.smpGroup$printReads.cmd, logger)

## breakpoint before SNP calling
testJobname = c(targets.smpGroup$printReads.jobname, targets.smpGroup$baserecal.jobname, targets.smpGroup$realnindel.jobname)
wait4jobs(testJobname, logger)

testOutputFiles = c(targets.smpGroup$printReadsBamOutFile)
if(delFiles & file.size(testOutputFiles) > 20000000){
	del.files(targets.smpGroup$smpBamFile, logger)
	del.files(targets.smpRealnBamFile$smpBamFile, logger)
}
#else{stop("something might wrong here: _merged_mkdup_smp_realnIndel_recal.bam\n")}

## split targets.bed by chrom
bed = fread(targets5bp_bed)
chrs = unique(bed$V1)
for(i in chrs){
	fwrite(bed[V1 == i,], quote=F, sep="\t", col.names=F, row.names=F, file=paste0(matricsDir, '/targets5bp_', i, '.bed' ))
}

## build pair 
targets.pair = fread(pairfile, header=F)
targets.pair = setNames(targets.pair, c('normal', 'tumor'))
targets.pair[, normalBamFile := paste0(initDir, '/', normal, '/', normal, '_merged_mkdup_smp_recal_realnIndel.bam')]
targets.pair[, tumorBamFile := paste0(initDir, '/', tumor, '/', tumor, '_merged_mkdup_smp_recal_realnIndel.bam')]
file.exists(targets.pair$tumorBamFile)
file.exists(targets.pair$normalBamFile)

## build chrPair
targets.chrPair = targets.pair[rep(1:nrow(targets.pair), each=length(chrs)), ]
targets.chrPair[, chr := rep(chrs, times=nrow(targets.pair))]
targets.chrPair[, chr.bed := paste0(matricsDir, '/targets5bp_', chr, '.bed')]
targets.chrPair

## call haplocalltyper by chromosome
##  https://software.broadinstitute.org/gatk/best-practices/workflow?=id=11145
## 4.0 version
## haplotypeCaller
## ImportGenomicsDB
## genotypeGVCFs
## VariantRecalibrator
## ApplyRecalibration
targets.germline = targets.chrPair[!duplicated(paste0(normal, chr.bed)),]
dim(targets.germline)
targets.germline[, haplocallGVCFChrFile := paste0(initDir, '/', normal, '/haplocall_', normal, '_', chr, '.g.vcf.gz')]
targets.germline[, haplocall.jobname := paste0('haplocall.', normal, '.', chr)]
targets.germline[, haplocall.cmd := bsub.head(haplocall.jobname, mem=30, cpu=10), by=1:nrow(targets.germline)]
targets.germline[, haplocall.cmd := paste0(haplocall.cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "'  HaplotypeCaller  -R ", genomeFasta, " -L ", chr.bed, " -I ", normalBamFile, " -O ", haplocallGVCFChrFile, " -ERC GVCF --dbsnp ", DB_SNP, " --max-reads-per-alignment-start 0 \"")]
targets.germline$haplocall.cmd[1]

exe.jobs(targets.germline$haplocall.cmd, logger)

file.exists(targets.germline$haplocallGVCFChrFile)

## combine haplocall g.vcf per sample
targets.germline[, {
	chrgvcf = paste(paste(' --variant ', haplocallGVCFChrFile), collapse=' ')
	haplocallGVCFFile = paste0(initDir, '/', normal, '/haplo_', .BY, '.g.vcf.gz')
	haplocat.jobname = paste0('haplocat.', .BY)
	pd = paste(haplocall.jobname, collapse=' ')
	cmd = bsub.head(haplocat.jobname, mem=30, cpu=10, postdone = pd)
	cmd = paste0(cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "'  CombineGVCFs ")
	cmd = paste0(cmd, " -R ", genomeFasta,  chrgvcf, " -O ", haplocallGVCFFile, " \"")
	exe.jobs(cmd, logger)
	}, by=normal]

## submitted above

file.exists(targets.germline$haplocallGVCFFile)

## combine all samples' gvcf files into one file
sampleGVcf= paste(paste(' --variant ', targets.germline$haplocallGVCFFile), collapse=' ')
haplocallGVCFFile = paste0(haploDir, '/haplo_combined_all_sample', '.g.vcf.gz')
combine.gvcf.jobname = 'haplo.combine.gvcf'
cmd = bsub.head(combine.gvcf.jobname, mem=10, cpu=1)
cmd = paste0(cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "'  CombineGVCFs ")
cmd = paste0(cmd, " -R ", genomeFasta,  ' ', sampleGVcf, " -O ", haplocallGVCFFile, " \"")

file.exists(targets.germline$haplocallGVCFFile)

## call genotype
haplocallVCFFile = paste0(haploDir, '/haplo_combined_all_sample', '.vcf.gz')
haplocall.genotype.jobname = 'haplotype.genotype'
cmd = bsub.head(haplocall.genotype.jobname, mem=50, cpu=10)
cmd = paste0(cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "'  GenotypeGVCFs")
cmd = paste0(cmd, " -R ", genomeFasta,  ' -V ', haplocallGVCFFile, " -O ", haplocallVCFFile, " \"")

##  recalibrating haplocall SNPs
haplocall.vcf.recal.jobname = 'haplocall.vcf.recal' 
haplocallVCFRecalFile 	      = paste0(haploDir, 'haplo_combined_all_sample_recal')
haplocallVCFRecalTranchesFile = paste0(haploDir, 'haplo_combined_all_sample_recal_tranches')
haplocallVCFRecalRscriptFile  = paste0(haploDir, 'haplo_combined_all_sample_recal_rscript.r')
cmd = bsub.head(haplocall.vcf.recal.jobname, mem=50, cpu=10, postdone = haplocall.genotype.jobname)
cmd = paste0(cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' VariantRecalibrator ")
cmd = paste0(cmd, " -R ", genomeFasta,  ' -V ', haplocallVCFFile, " --resource hapmap,known=false,training=true,truth=true,prior=15.0:", HAPMAP)
cmd = paste0(cmd, " --resource omni,known=false,training=true,truth=false,prior=12.0:", OMNI_1000G)
cmd = paste0(cmd, " --resource 1000G,known=false,training=true,truth=false,prior=10.0:", PHASE1_SNPS_1000G)
cmd = paste0(cmd, " --resource dbsnp,known=true,training=false,truth=false,prior=2.0:", DB_SNP)
cmd = paste0(cmd, " -an QD -an MQ -an MQRankSUM -an ReadPosRankSum -an FS -an SOR")
cmd = paste0(cmd, " -mode BOTH -O ", haplocallVCFRecalFile)
cmd = paste0(cmd, " --tranches-file ", haplocallVCFRecalTranchesFile)
cmd = paste0(cmd, " --rscript-file ", haplocallVCFRecalRscriptFile, " \"")
cmd

exe.jobs(cmd, logger)

##  apply haplocall recal
## be noticed the haplocallRecalVCFFile is the vcf file
## be noticed the haplocallVCFRecalFile is the recal parameter file
haplocallRecalVCFFile = paste0(haploDir, '/haplo_combined_all_sample_recal', '.vcf.gz')
haplocall.vcf.recal.apply.jobname = 'haplocall.vcf.recal.apply'
cmd = bsub.head(haplocall.vcf.recal.apply.jobname, mem=50, cpu=10, postdone = haplocall.vcf.recal.jobname)
cmd = paste0(cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "'  ApplyVQSR ")
cmd = paste0(cmd, " -R ", genomeFasta,  ' -V ', haplocallVCFFile, ' -O ', haplocallRecalVCFFile)
cmd = paste0(cmd, " -ts_filter_level 99.0 --tranches-file ", haplocallVCFRecalTranchesFile)
cmd = paste0(cmd, " --recal-file ", haplocallVCFRecalFile, " -mode BOTH \" ")

exe.jobs(cmd, logger)

##  Filter haplocall VCF
haplocallFilterVCFFile = paste0(haploDir, '/', pre, '_haplo_combined_all_sample_recal_filter', '.vcf.gz')
haplocall.vcf.filter.jobname = 'haplocall.vcf.filter'
cmd = bsub.head(haplocall.vcf.filter.jobname, mem=50, cpu=10, postdone = haplocall.vcf.recal.apply.jobname)
cmd = paste0(cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "'  FilterVcf ")
cmd = paste0(cmd, " -I ", haplocallRecalVCFFile, ' -O ', haplocallFilterVCFFile, ' -R ', genomeFasta)
cmd = paste0(cmd, " --MIN_AB 0.1 ") ## minimal allele frequency for germline mutations, Vaf
cmd = paste0(cmd, " --MIN_DP 5 ") ## minimal depth for germline mutations

exe.jobs(cmd, logger)

## generate pon file, panel of normal for mutect2
haplocallVCFFile 
ponfile = paste0(initDir, '/pon.vcf.gz')
cmd = bsub.head(jobname ='pon', , mem=30, cpu = 5)
cmd = paste0(cmd, " \"", GATK4, " --java-options '-Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' CreateSomaticPanelOfNormals ")
cmd = paste0(cmd, " -vcfs ", haplocallFilterVCFFile)
cmd = paste0(cmd, " -O ", ponfile, "\"")
cmd

exe.jobs(cmd, logger)

## mutect2 by chromosome
targets.chrPair[, mutect2callVCFChrFile  := paste0(initDir, '/', tumor, '/mutect2_', tumor, '_', chr, '.vcf')]
targets.chrPair[, mutect2call.jobname := paste0('mutect2call.', chr, '.', tumor)]
targets.chrPair[, mutect2call.cmd := bsub.head(mutect2call.jobname, mem=30, cpu=2, We='8:26', cwd=cwd)]
targets.chrPair[, mutect2call.cmd := paste0(mutect2call.cmd, " \"", GATK4, " --java-options '-Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' Mutect2")]
targets.chrPair[, mutect2call.cmd := paste0(mutect2call.cmd, " --reference ", genomeFasta, " -I ", tumorBamFile, " -tumor ", tumor, " -I ", normalBamFile, " -normal ", normal)]
targets.chrPair[, mutect2call.cmd := paste0(mutect2call.cmd, " --germline-resource ", gnomad, " --af-of-alleles-not-in-resource 0.00003125 --panel-of-normals ", ponfile, " -O ", mutect2callVCFChrFile, " \"")]
targets.chrPair$mutect2call.cmd[1]

exe.jobs(targets.chrPair$mutect2call.cmd, logger)

## combine mutect2 VCF for each sample
targets.chrPair[, {
	chrvcf = paste(paste(' -I ', mutect2callVCFChrFile, sep=' '), collapse = ' ')
	mutect2callVCFFile = paste0(initDir, '/', tumor, '/mutect2_', tumor, '_all_chr.vcf')
	mutect2cat.jobname = paste0('mutect2cat.', .BY)
	pd = paste(mutect2call.jobname, collapse= ' ')
	cmd = bsub.head(mutect2cat.jobname, mem = 10,  cpu=1, postdone = pd)
	cmd = paste0(cmd, " \" ", GATK4, " --java-options  '-Xms256m, -Xmx10g -XX:UseGCOverheadingLimit -Djava.io.tmpdir=", tmpdir, " '  GatherVcfs ")
	cmd = paste0(cmd, ' -R ', genomeFasta, ' ', chrvcf, ' -o', mutect2callVCFFile, "\"")
	exe.jobs(cmd, logger)
	}, by=tumor]

targets.pair[, mutect2callVCFFile = paste0(initDir, '/', tumor, '/mutect2_', tumor, '_all_chr.vcf')]
targets.pair[, mutect2cat.jobname = paste0('mutect2cat.', tumor)]

## combine mutect2 vcf files for all samples
mutect2.combine.jobname = 'mutect2.combine'
mutect2SampleVCFs = paste(paste(' -I ', targets.chrPair$mutect2callVCFFile, sep=' '), collapse = ' ')
mutect2VCFFile = paste0(mutectDir, '/mutect2_all_sample.vcf')
pd = paste(targets.pair$mutect2cat.jobname, collapse= ' ')
cmd = bsub.head(mutect2.combine.jobname, mem = 10, cpu = 1, postdone = pd)
cmd = paste0(cmd, "\" ", GATK4, " --java-options '-Xms256m, -Xmx10g -XX:UseGCOverheadingLimit -Djava.io.tmpdir=", tmpdir, " ' GatherVcfs ")
cmd = paste0(cmd,  " -R ", genomeFasta, ' -o ', mutect2VCFFile,  "\"")

## estimate cross sample contamination: step 0:
mutect2biallele = paste0(mutectDir, '/biallele.vcf')
mutect2.selectbiallele.jobname = 'mutect2.biallele'
cmd = bsub.head(mutect2.selectbiallele.jobname, mem=50, cpu=10, postdone = mutect2.vcf.recal.apply.jobname)
cmd = paste0(cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' SelectVariants")
cmd = paste0(cmd, " --restrict-alleles-to BIALLELIC -R ", genomeFasta, ' -V ', mutect2VCFFile, " -O ", mutect2biallele, "\"")

exe.jobs(cmd, logger)

## estimate cross sample contamination: step 1:
mutect2biallelePileUpFile = paste0(mutectDir, '/mutect2_combined_all_sample_recal_filter', '.vcf.gz')
mutect2.getpileup.jobname = 'mutect2.pileup'
cmd = bsub.head(mutect2.getpileup.jobname, mem=50, cpu=10, postdone = mutect2.vcf.recal.apply.jobname)
cmd = paste0(cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' GetPileupSummaries")
cmd = paste0(cmd, " -I ", tumorBamFile, ' -V ', mutect2biallele, " -O ", mutect2biallelePileUpFile, "\"")

exe.jobs(cmd, logger)

## estimate cross sample contamination: step 2:
mutect2contaminationFile = paste0(mutectDir, '/mutect2_combined_all_sample_recal_filter', '.vcf.gz')
mutect2.cal.contamintation.jobname = 'mutect2.calContamination'
cmd = bsub.head(mutect2.cal.contamination.jobname, mem=50, cpu=10, postdone = mutect2.vcf.recal.apply.jobname)
cmd = paste0(cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' CalculateContamination")
cmd = paste0(cmd, " -I ", mutect2biallelePileUpFile, ' -O ', mutect2contaminationFile, "\"")

exe.jobs(cmd, logger)

## there is no recal vcf for mutect results
mutect2FilterVCFFile = paste0(mutectDir, '/mutect2_combined_all_sample_recal_filter', '.vcf.gz')
mutect2.vcf.filter.jobname = 'mutect2.vcf.filter'
cmd = bsub.head(mutect2.vcf.filter.jobname, mem=50, cpu=10, postdone = mutect2.vcf.recal.apply.jobname)
cmd = paste0(cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "'  FilterMutectCalls")
cmd = paste0(cmd, " -V ", mutect2VCFFile, ' -O ', mutect2FilterVCFFile)
cmd = paste0(cmd, " --contamination-table mutect2contaminationFile \"")

exe.jobs(cmd, logger)

## somaticsniper
targets.pair[, snipercallVCFOutFile  := paste0(sniperDir, '/', tumor, '.vcf')]
targets.pair[, snipercall.jobname := paste0('snipercall.', tumor)]
#targets.pair[, snipercall.cmd := bsub.head(snipercall.jobname, mem=90, cpu=30, We='8:26', cwd=cwd, postdone = baserecal.jobname)]
targets.pair[, snipercall.cmd := bsub.head(snipercall.jobname, mem=4, cpu=2, We='8:26', cwd=cwd)]
targets.pair[, snipercall.cmd := paste0(snipercall.cmd, ' "', SOMATICSNIPER, '/bam-somaticsniper -F vcf -f ', genomeFasta, ' -q 1 ', tumorBamFile, ' ', normalBamFile, ' ',  snipercallVCFOutFile, ' "')]
targets.pair$snipercall.cmd[1]

exe.jobs(targets.pair$snipercall.cmd, logger)

## wait for all the mutation calller finish.
wait4jobs()

## combine vcf from different calling
## now just from mutect2 and haplotypecall
## mutect2FilterVCFFile
## haplocallFilterVCFFile
callerCombinedVCFFile = paste0(varDir, '/', pre, '_caller_combined.vcf.gz')
merge.caller.vcf.jobname = 'merge.caller.vcf'
cmd = bsub.head(merge.caller.vcf.jobname, mem=50, cpu=10)
cmd = paste0(cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' MergeVcfs")
cmd = paste0(cmd, " -I ", mutect2FilterVCFFile, ' -I ', haplocallFilterVCFFile, " -O ", callerCombinedVCFFile, " \"")

exe.jobs(cmd, logger)

## split vcf for vcf2maf to add tumor-id
## working here
## callerCombinedVCFFile
targets.pair[, callerCombinedSampleVCFFile = paste0(initDir, '/', tumor, '/caller_combined_', tumor, '.vcf')] ; 
targets.pair[, vcfsplit.jobname := paste0('vcfsplit.', tumor)]
targets.pair[, vcfsplit.cmd := bsub.head(vcfsplit.jobname, mem=10, cpu=1, We='8:26', cwd=cwd, postdone = merge.caller.vcf.jobname)
targets.pair[, vcfsplit.cmd := paste0(vcfsplit.cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' SelectVaraints")
targets.pair[, vcfsplit.cmd := paste0(cmd, " -R ", genomeFasta, ' -V ', callerCombinedVCFFile, " --sample-name ", tumor)
targets.pair[, vcfsplit.cmd := paste0(cmd, " -O ", callerCombinedSampleVCFFile, " \"")]
targets.pair$vcfsplit.cmd[1]

exe.jobs(targets.pair$vcfsplit.cmd, logger) # 

## vcf2maf for each sample after combined vcf from different callers
ncbi = 'GRCh37'
targets.pair[, sampleMafFile := paste0(initDir, '/', tumor, '/caller_combined_', tumor, '.maf')] ; 
targets.pair[, vcf2maf.jobname := paste0('vcf2maf.', tumor)]
targets.pair[, vcf2maf.cmd := bsub.head(vcf2maf.jobname, mem=30, cpu=3, We='8:26', cwd=cwd, postdone = vcfsplit.jobname)]
targets.pair[, vcf2maf.cmd := paste0(vcf2maf.cmd, " \"", PERL, "/perl ", VCF2MAF, "/vcf2maf.pl --input ", callerCombinedSampleVCFFile)]
targets.pair[, vcf2maf.cmd := paste0(vcf2maf.cmd, " --output-maf ", sampleMafFile)]
targets.pair[, vcf2maf.cmd := paste0(vcf2maf.cmd, ' --ref-fasta ', genomeFasta, ' --tmp-dir ', tmpdir, '  --ncbi ', ncbi, ' --vep-path ', VEP, ' --vep-data ', VEP)]
targets.pair[, vcf2maf.cmd := paste0(vcf2maf.cmd, " --tumor-id ", tumor, " --normal-id ", normal, " --filter-mutect ", ExAC_VCF, " \"")]
targets.pair$vcf2maf.cmd[2]

exe.jobs(targets.pair$vcf2maf.cmd, logger) 

## combine maf files from diferent samples
combinedMafFile = paste0(varDir, '/', pre, '.maf')
sampleMafFiles = paste(targets.pair$sampleMafFile, collapse=' ')
combine.maf.jobname = 'combine.maf'
combineMaf.cmd := bsub.head(combine.maf.jobname, mem=10, cpu=1, We='8:26', cwd=cwd, postdone = paste(targets.pair$vc2maf.jobname, collapse=' ')]
combineMaf.cmd := paste0(combineMaf.cmd, " \"cat ",  sampleMafFiles, " > ", combinedMafFile, " \"")
combineMaf.cmd

exe.jobs(combineMaf.cmd, logger)

## facets to get copy number and CCF 
## MafAnnotate
MINCOV = 0
BASEQ = 20
MAPQ = 15
targets.pair[, facetsCountNOutFile  := paste0(facetsDir, '/sample/', normal, '_merged_mkdup_smp_realnIndel_recal.bam.dat')]
targets.pair[, facetsCountN.jobname := paste0('facetsCountN.', normal)]
targets.pair[, facetsCountN.cmd := bsub.head(facetsCountN.jobname, mem=8, cpu=4, We='8:26', cwd=cwd)]
targets.pair[, facetsCountN.cmd := paste0(facetsCountN.cmd, ' "', facets, '/bin/GetBaseCounts  --thread 4 --filter_improper_pair --sort_output --fasta ')]
targets.pair[, facetsCountN.cmd := paste0(facetsCountN.cmd, genomeFasta, ' --vcf ', FACETS_DB_SNP, ' --maq ', MAPQ, ' --baq ', BASEQ, ' --cov ', MINCOV)]
targets.pair[, facetsCountN.cmd := paste0(facetsCountN.cmd, ' --bam ', normalBamFile, ' --out ', facetsCountNOutFile, ' "')] 
targets.pair$facetsCountN.cmd[1]

exe.jobs(targets.pair$facetsCountN.cmd[!duplicated(targets.pair$normalBamFile)], logger)

targets.pair[, facetsCountTOutFile  := paste0(facetsDir, '/sample/', tumor, '_merged_mkdup_smp_realnIndel_recal.bam.dat')]
targets.pair[, facetsCountT.jobname := paste0('facetsCountT.', tumor)]
targets.pair[, facetsCountT.cmd := bsub.head(facetsCountT.jobname, mem=8, cpu=4, We='8:26', cwd=cwd)]
targets.pair[, facetsCountT.cmd := paste0(facetsCountT.cmd, ' "', facets, '/bin/GetBaseCounts  --thread 4 --filter_improper_pair --sort_output --fasta ')]
targets.pair[, facetsCountT.cmd := paste0(facetsCountT.cmd, genomeFasta, ' --vcf ', FACETS_DB_SNP, ' --maq ', MAPQ, ' --baq ', BASEQ, ' --cov ', MINCOV)]
targets.pair[, facetsCountT.cmd := paste0(facetsCountT.cmd, ' --bam ', tumorBamFile, ' --out ', facetsCountTOutFile, ' "')] 
targets.pair$facetsCountT.cmd[1]

exe.jobs(targets.pair$facetsCountT.cmd[!duplicated(targets.pair$tumorBamFile)], logger)

## facets mergeTN.R
targets.pair[, facetsMergeOutFile  := paste0(facetsDir, '/sample/count_merged_', tumor, '_', normal, '.dat.gz')]
targets.pair[, facetsMerge.jobname := paste0('facetsMerge.', tumor)]
targets.pair[, facetsMerge.cmd := bsub.head(facetsMerge.jobname, mem=18, cpu=4, We='8:26', cwd=cwd, c(facetsCountT.jobname, facetsCountN.jobname))]
#targets.pair[, facetsMerge.cmd := bsub.head(facetsMerge.jobname, mem=18, cpu=4, We='8:26', cwd=cwd)]
targets.pair[, facetsMerge.cmd := paste0(facetsMerge.cmd, ' "', facets, '/mergeTN.R ', facetsCountTOutFile, ' ', facetsCountNOutFile, ' ', facetsMergeOutFile, ' "')]
targets.pair$facetsMerge.cmd[1]

exe.jobs(targets.pair$facetsMerge.cmd[!duplicated(targets.pair$tumorBamFile)], logger)

## call facets
## facets.suite facets.lib dir tag file ggenome pc c
pcval = 200
cval = 100

# below is the same as the facets/facets_RUN.sh
targets.pair[, facetsRunOutDir  := paste0(facetsDir, '/', tumor)] ; 
for(i in 1:nrow(targets.pair)){system(paste0('mkdir -p ', targets.pair$facetsRunOutDir[i]))}
targets.pair[, facetsRun.jobname := paste0('facetsRun.', tumor)]
targets.pair[, facetsRun.cmd := bsub.head(facetsRun.jobname, mem=3, cpu=3, We='8:26', cwd=cwd, postdone = targets.pair$facetsRun.jobname)]
targets.pair[, facetsRun.cmd := paste0(facetsRun.cmd, ' "', FACETS_SUITE, '/doFacets.R', ' --cval ', cval, ' --purity_cval ', pcval, ' --genome hg19')]
targets.pair[, facetsRun.cmd := paste0(facetsRun.cmd, ' --counts_file ', facetsMergeOutFile, ' --TAG ', tumor, ' --directory ', facetsRunOutDir, ' --R_lib ', RLIB_PATH, ' "')]
targets.pair$facetsRun.cmd[2]

exe.jobs(targets.pair$facetsRun.cmd, logger)

wait4jobs(targets.pair$facetsRun.jobname)

## combine seg file
targets.pair[, facetsSegFile := paste0(facetsRunOutDir, '/', tumor, '_hisens.seg')]
targets.pair[, facetsSegFile]
seg.list = lapply(targets.pair$facetsSegFile, fread, sep="\t")
rbindlist(seg.list) -> seg.all
seg.all
fwrite(seg.all, file=paste0(facetsDir, '/seg_all_hisen.seg'), sep = '\t')
rm(seg.list, seg.all)

## gene level
source('~/program/facets-suite-1.0.1/geneLevel.R')
impact341_targets = '/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv3/picard_targets.interval_list'
IMPACT341_targets = suppressWarnings(fread(paste0('grep -v "^@" ', impact341_targets)))
setnames(IMPACT341_targets, c("chr", "start", "end", "strand", "name"))
setkey(IMPACT341_targets, chr, start, end)
head(IMPACT341_targets)

targets.pair[, facetsCncfFile := paste0(facetsDir, '/', tumor, '/', tumor, '_hisens.cncf.txt')]
targets.pair[, facetsGeneLevelFile := paste0(facetsDir, '/', tumor, '/', tumor, '_hisens.geneLevel.txt')]
for(i in 1:nrow(targets.pair)){
	cncf_file = targets.pair$facetsCncfFile[i]
	get_gene_level_calls(cncf_files = cncffile, gene_targets = IMPACT341_targets) -> tmp
	tmp[, Tumor_Sample_Barcode := targets.pair$tumor[i]]
	head(tmp)
	fwrite(tmp, file=targets.pair$facetsGeneLevelFile[i], sep="\t")
}

# combine gene level into one
gl.list = lapply(targets.pair$facetsGeneLevelFile, fread, sep="\t")
rbindlist(gl.list) -> gl.all
fwrite(gl.all, file=paste0(facetsDir, '/', pre, '_facets_gene_level'), sep = '\t')
rm(gl.list, gl.all)

facetsRunOutDir  
## get facets.files for mafAnno
targets.pair[, rdata.file := paste0(facetsRunOutDir, '/', tumor, '.Rdata')] ## need revise
tmp = targets.pair[, .(tumor, rdata.file)]
setNames(tmp, c('Tumor_Sample_Barcode', 'Rdata_filename'))
facets.rdata.file = paste0(facetsDir, '/rdata_file.txt')
fwrite(tmp, sep="\t", quote=F, col.names=T, row.names=F, file=rdata.file)

## mafAnno: add CNV to each haplotypecaller maf files
## no postdone ever since after facets
combinedMafAnnoMafFile = paste0(varDir, '/', pre, '_mafAnno.maf')
mafanno.jobname = 'combine.maf'
mafanno.cmd := bsub.head(mafanno.jobname, mem=10, cpu=3, We='8:26', cwd=cwd)
mafanno.cmd := paste0(mafanno.cmd, " \"", FACETS_SUITE, "/mafanno.R -m ", combinedMafFile, " -f ", facets.rdata.file, " -o ", combinedMafAnnoMafFile, " \"")
mafanno.cmd

exe.jobs(mafanno.cmd, logger) 

## mafAnno: add CNV to each mutect maf files
#targets.pair[, mutectmafAnnoOutFile  := paste0(mutectDir, '/', tumor, '_CNV.maf')] ; 
#targets.pair[, mutectmafAnno.jobname := paste0('mutectmafAnno.', tumor)]
#targets.pair[, mutectmafAnno.cmd := bsub.head(mutectmafAnno.jobname, mem=3, cpu=3, We='8:26', cwd=cwd)]
#targets.pair[, mutectmafAnno.cmd := paste0(mutectmafAnno.cmd, ' "', FACETS_SUITE, '/mafAnno_v2.R -m ', mutect2mafOutFile, ' -f ', facetsHisensRdata, ' -o ', mutectmafAnnoOutFile, ' "')]
#targets.pair$mutectmafAnno.cmd[2]
#
#exe.jobs(targets.pair$mutectmafAnno.cmd, logger) # 

# filter maf file: common mutations
combinedMafAnnoMafFile
ngsfilterCommonMafFile = paste0(varDir, '/', pre, '_mafAnno_filterCommon.maf')
ngsfilter.common.jobname = 'ngsfilter.common'
ngsfilter.common.cmd := bsub.head(ngsfilter.common.jobname, mem=10, cpu=3, We='8:26', cwd=cwd)
ngsfilter.common.cmd := paste0(ngsfilter.common.cmd, " \"", NGSFilter, "/filter_common_variants.R -m ", combinedMafAnnoMafFile, " -o ", combinedMafAnnoMafFile, "\"")
ngsfilter.common.cmd

exe.jobs(ngsfilter.common.cmd, logger) 

# filter maf file: low confidence mutations
ngsfilterLowConfidenceMafFile = paste0(varDir, '/', pre, '_mafAnno_filterCommon_lowConfidence.maf')
ngsfilter.low.conf.jobname = 'ngsfilter.low.conf'
ngsfilter.low.conf.cmd := bsub.head(ngsfilter.low.conf.jobname, mem=10, cpu=3, We='8:26', cwd=cwd)
ngsfilter.low.conf.cmd := paste0(ngsfilter.low.conf.cmd, " \"", NGSFilter, "/filter_low_conf.R -m ", ngsfilterCommonMafFile, " -o ", ngsfilterLowConfidenceMafFile, "\"")
ngsfilter.low.conf.cmd

exe.jobs(ngsfilter.common.cmd, logger) 

# ./filter_blacklist_regions.R -m input.maf -o output.maf
ngsfilterLowConfidenceMafFile
ngsfilterBlacklistMafFile = paste0(varDir, '/', pre, '_mafAnno_filterCommon_lowConfidence_blacklist.maf')
ngsfilter.blacklist.jobname = 'ngsfilter.blacklist'
ngsfilter.blacklist.cmd := bsub.head(ngsfilter.blacklist.jobname, mem=10, cpu=3, We='8:26', cwd=cwd)
ngsfilter.blacklist.cmd := paste0(ngsfilter.blacklist.cmd, " \"", NGSFilter, "/filter_blacklist.R -m ", ngsfilterLowConfidenceMafFile, " -o ", ngsfilterBlacklistMafFile, "\"")
ngsfilter.blacklist.cmd

exe.jobs(ngsfilter.blacklist.cmd, logger) 

# ./tag_hotspots.py -m input.maf -itxt hotspot-list-union-v1-v2.txt -o output.maf
ngsfilterBlacklistMafFile
ngsfilterHotspotMafFile = paste0(varDir, '/', pre, '_mafAnno_filterCommon_lowConfidence_blacklist_hotspot.maf')
ngsfilter.hotspot.jobname = 'ngsfilter.hotspot'
ngsfilter.hotspot.cmd := bsub.head(ngsfilter.hotspot.jobname, mem=10, cpu=3, We='8:26', cwd=cwd)
ngsfilter.hotspot.cmd := paste0(ngsfilter.hotspot.cmd, " \"", NGSFilter, "/tag_hotspots.py -m ", ngsfilterBlacklistMafFile, " -o ", ngsfilterHotspotMafFile, " -itxt hotspot-list-union-v1-v2.txt \"")
ngsfilter.hotspot.cmd

exe.jobs(ngsfilter.hotspot.cmd, logger) 

# hybrid genome

targets.pair[, bicNormalBam := paste0('Proj_07813_DF/r_001/alignments/Proj_07813_DF_indelRealigned_recal_s_', gsub("-", "_", normal), ".bam") ]
targets.pair[, bicTumorBam := paste0('Proj_07813_DF/r_001/alignments/Proj_07813_DF_indelRealigned_recal_s_', gsub("-", "_", tumor), ".bam") ]
file.exists(targets.pair$bicNormalBam)
file.exists(targets.pair$bicTumorBam)

## examine HLA types by polysolver
targets.pair[, hla.polysolverOutDir  := paste0(initDir, '/', normal, '/polysolver')] ; 
targets.pair[, {system(paste0('mkdir -p ', hla.polysolverOutDir))}, by=1:nrow(targets.pair)]

targets.pair[, hla.polysolver.jobname := paste0('hla.polysolver.', normal)]
targets.pair[, hla.polysolver.cmd := bsub.head(hla.polysolver.jobname, mem=30, cpu=3, We='8:26', cwd=cwd)] 
targets.pair[, hla.polysolver.cmd := paste0(hla.polysolver.cmd, " \" bash ", polysolver, "/scripts/shell_call_hla_type")]
targets.pair[, hla.polysolver.cmd := paste0(hla.polysolver.cmd, " ", bicNormalBamFile, " Unknown ", " 1 hg19 STDFQ")] 
targets.pair[, hla.polysolver.cmd := paste0(hla.polysolver.cmd, " 0 ", hla.polysolverOutDir, " \"")]
targets.pair$hla.polysolver.cmd[3]
cmd = targets.pair$hla.polysolver.cmd[!duplicated(targets.pair$normal)]
cmd[2]

exe.jobs(cmd, logger)

## polysolver step 2
targets.pair[, hla.polysolver2.jobname := paste0('hla.polysolver2.', normal)]
targets.pair[, hla.polysolver2.cmd := bsub.head(hla.polysolver2.jobname, mem=30, cpu=3, We='8:26', cwd=cwd, postdone=hla.polysolver.jobname), by=1:nrow(targets.pair)] 
targets.pair[, hla.polysolver2.cmd := paste0(hla.polysolver2.cmd, " \" bash ", polysolver, "/scripts/shell_call_hla_mutations_from_type")]
targets.pair[, hla.polysolver2.cmd := paste0(hla.polysolver2.cmd, " ", normalBamFile, " ", tumorBamFile, " ", hla.polysolverOutDir, "/winners.hla.txt hg19 STDFQ")] 
targets.pair[, hla.polysolver2.cmd := paste0(hla.polysolver2.cmd, " ", hla.polysolverOutDir, " \"")]
targets.pair$hla.polysolver2.cmd[3]
cmd = targets.pair$hla.polysolver2.cmd[!duplicated(targets.pair$normal)]
cmd[2]

exe.jobs(cmd, logger)

## polysolver step 3
targets.pair[, hla.polysolver3.jobname := paste0('hla.polysolver3.', normal)]
targets.pair[, hla.polysolver3.cmd := bsub.head(hla.polysolver3.jobname, mem=30, cpu=3, We='8:26', cwd=cwd, postdone=hla.polysolver2.jobname), by=1:nrow(targets.pair)] 
targets.pair[, hla.polysolver3.cmd := paste0(hla.polysolver3.cmd, " \" bash ", polysolver, "/scripts/shell_annotate_hla_mutations indiv ")]
targets.pair[, hla.polysolver3.cmd := paste0(hla.polysolver3.cmd, " ", hla.polysolverOutDir, " \"")]
targets.pair$hla.polysolver3.cmd[3]
cmd = targets.pair$hla.polysolver3.cmd[!duplicated(targets.pair$normal)]
cmd[2]

exe.jobs(cmd, logger)

## examine HLA by SOAP-HLA
polysolver
targets.pair[, soapHlaOutDir  := paste0(initDir, '/', normal, '/soapHla')] ; 
targets.pair[, {system(paste0('mkdir -p ', soapHlaOutDir))}, by = 1:nrow(targets.pair)]

targets.pair[, soaphla.jobname := paste0('soaphla.', normal)]
targets.pair[, soaphla.cmd := bsub.head(soaphla.jobname, mem=30, cpu=3, We='8:26', cwd=cwd)] 
targets.pair[, soaphla.cmd := paste0(soaphla.cmd, " \" ", PERL, '/perl ', SOAPHLA, "/MHC_autopipeline.pl -i")]
targets.pair[, soaphla.cmd := paste0(soaphla.cmd, " ", normalBamFile, " -od ", soapHlaOutDir, " -v hg19 \"")] 
targets.pair$soaphla.cmd[1]
cmd = targets.pair$soaphla.cmd[!duplicated(targets.pair$normal)]
cmd[2]

exe.jobs(cmd, logger)

## maf2vcf.pl to get vcf files
neoDir = paste0(resDir, '/neo')
if(!dir.exists(neoDir)) {dir.create(neoDir)}
maffile = 'scc_oncokb.maf'
maf2vcf.jobname = 'maf2vcf'
cmd = bsub.head(jobname =maf2vcf.jobname, mem=30, cpu = 5)
cmd = paste0(cmd, " \"", PERL, "/perl ",  VCF2MAF, "/maf2vcf.pl --output-dir ", neoDir, " --input-maf ", maffile, " \"")
cmd

exe.jobs(cmd, logger)

## annotate vcf file
##vcfFile = paste0(neoDir, '/', sub("\\.maf", ".vcf", maffile))
##annoVcfFile = paste0(neoDir, '/anno_', sub("\\.maf", ".vcf", maffile))
##annotate.vcf.jobname = 'annotate.vcf'
##cmd = bsub.head(annotate.vcf.jobname, mem=40, cpu=44, postdone=maf2vcf.jobname)
##cmd = paste0(cmd, " \"", VEPv88, " -i ", vcfFile, " -o ", annoVcfFile, " -offline --assembly GRCh37 --species human --format vcf --vcf --symbol --plugin Downstream --plugin Wildtype --terms SO --dir_plugins ", VEP_plugins, " -fork 40\"") 
##cmd = paste0(" ", VEPv88, " -i ", vcfFile, " -o ", annoVcfFile, " -offline --assembly GRCh37 --format vcf --vcf --symbol --plugin Downstream --plugin Wildtype --terms SO --dir_plugins ", VEP_plugins, " -fork 30") 
##cmd
##system(cmd)

#write.table(cmd, file='x')
#exe.jobs(cmd, logger)

## split vcf per sample 
targets.pair[, bicSampleName := paste0('s_', gsub("-", "_", tumor))]
targets.pair[, mafVcfFile := paste0(initDir, '/', tumor, '/neo_', tumor, '.vcf')]
targets.pair[, splitmafvcf.jobname := paste0("splitmafvcf.", tumor)]
targets.pair[, splitmafvcf.cmd := bsub.head(splitmafvcf.jobname, mem=10, We='1:11', cpu=1, postdone = maf2vcf.jobname), by=1:nrow(targets.pair)]
targets.pair[, splitmafvcf.cmd := paste0(splitmafvcf.cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' SelectVariants")]
targets.pair[, splitmafvcf.cmd := paste0(splitmafvcf.cmd, " -R ", genomeFasta, ' -V ', vcfFile, " --sample-name ", bicSampleName, " -O ", mafVcfFile, "\"")]
targets.pair[1, splitmafvcf.cmd]

exe.jobs(targets.pair$splitmafvcf.cmd, logger)

## annotate vcf per sample
targets.pair[, mafAnnoVcfFile := paste0(initDir, '/', tumor, '/neo_anno_', tumor, '.vcf')]
targets.pair[, annovcf.jobname := paste0("annovcf.", tumor)]
targets.pair[, annovcf.cmd := bsub.head(annovcf.jobname, mem=20, We='1:11', cpu=5, postdone = splitmafvcf.jobname), by=1:nrow(targets.pair)]
targets.pair[, annovcf.cmd := paste0(annovcf.cmd, " \"", VEPv88, " -i ", mafVcfFile, " -o ", mafAnnoVcfFile, " -offline --assembly GRCh37 --format vcf --vcf --symbol --plugin Downstream --plugin Wildtype --terms SO --dir_plugins ", VEP_plugins, " -fork 4 --force_overwrite\"") ]
targets.pair[1, annovcf.cmd]

exe.jobs(targets.pair$annovcf.cmd, logger)

## get protein fasta
targets.pair[, mutPepFile := paste0(initDir, '/', tumor, '/mut_', tumor, '.fasta')]
targets.pair[, mut.pep.jobname := paste0("mut.pep.", tumor)]
targets.pair[, mutpep.cmd := bsub.head(mut.pep.jobname, mem=10, We='1:11', cpu=1, postdone = annovcf.jobname), by = 1:nrow(targets.pair)]
targets.pair[, mutpep.cmd := paste0(mutpep.cmd, " \"", pvacseq, "  generate_protein_fasta ", mafAnnoVcfFile, " 21 ", mutPepFile, " \"")]
targets.pair[1, mutpep.cmd]

exe.jobs(targets.pair$mutpep.cmd, logger)

## run neoantigen prediction
#pvacseq run  DS-bla-185-T2 <- anno.vcf  Test  HLA-G*01:09,HLA-E*01:01,H2-IAb  NetMHC PickPocket NNalign output  -e 9,10
targets.pair[, hlaFile := paste0(initDir, '/', normal, '/winners.hla.nofreq.txt')]
targets.pair[, allele := {
	hlaFile = targets.pair$hlaFile[3];fread(hlaFile, header=F) -> hla; hla
	hlaFile = targets.pair$hlaFile[3];
	scan(hlaFile, character(), sep="\n") -> hla; 
	hla = strsplit(hla, "\t")
	names(hla) = lapply(hla, "[[", 1)
	hla = lapply(hla, function(x){x[2:length(x)]})
	hla2 = lapply(hla, function(x){strsplit(x, "_") -> xx; unlist(lapply(xx, function(xx){paste0(toupper(xx[1]), "-", toupper(xx[2]), "*", paste(xx[3:length(xx)], collapse=":"), collapse="")}))})
	#hla2 = lapply(hla, function(x){strsplit(x, "_") -> xx; unlist(lapply(xx, function(xx){paste0(toupper(xx[1]), "-", toupper(xx[2]), "*", xx[3], ':', xx[4], collapse="")}))})
	hla2 = paste(unlist(hla2), collapse=",")
	hla2 }, by=1:nrow(targets.pair)]
targets.pair$allele

targets.pair[, sampleNeoDir := paste0(initDir, '/', tumor, '/neoantigeon')]
targets.pair[, {system(paste0("mkdir -p ", sampleNeoDir))}, by = 1:nrow(targets.pair)]

targets.pair[, neo.jobname := paste0("neo.", tumor)]
targets.pair[, neo.cmd := bsub.head(neo.jobname, mem=20, We='1:11', cpu=1, postdone = mut.pep.jobname)]
targets.pair[, neo.cmd := paste0(neo.cmd, " \"", pvacseq, " run ", mafAnnoVcfFile, " ", tumor, " ", allele, " NetMHCpan ", sampleNeoDir, " -e 9,10 \"")]
targets.pair[1, neo.cmd]


exe.jobs(targets.pair$neo.cmd[1], logger)

## peronsal genome
maf = fread(ngsfilterHotspotMafFile)
vc.sel = c('Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'In_Frame_Ins', 'In_Frame_Del')
maf.neo = maf[Variant_Classification %in% vc.sel, .(Tumor_Sample_Barcode, HGVSp_Short, SWISSPROT, Transcript_ID, HGVSc, HGVSp, Variant_Classification, Hugo_Symbol)]
maf.neo[, loc := sub(".*?([0-9]+).*", "\\1", HGVSp_Short)]
maf.neo[, loc2 := {xx = sub("([0-9]+)", "loc", HGVSp_Short); xx = sub(".*?([0-9]+).*", "\\1", xx); xx}]
maf.neo = maf.neo[Variant_Classification == 'Missense_Mutation',]
maf.neo = copy(maf[Variant_Classification == 'Missense_Mutation', .(Tumor_Sample_Barcode, HGVSp_Short, SWISSPORT, Transcript_ID, HGVSc, HGVSp, Variant_Classification, Hugo_Symbol, loc, )])
maf.neo[, t.seq := '']
maf.neo[, n.seq := '']
maf.neo = maf.neo[loc != '',]
maf.neo[, id := paste0('id', 1:nrow(maf.neo))]
maf.neo[, tagN := paste0(">", id, 'N')]
maf.neo[, tagT := paste0(">", id, 'T')]

## get 25 aa, with the middle changed aa
for(i in 1:nrow(tmp)){
	pepi = pep.mmu[[tmp$ensembl_peptide_id[i]]]
	pepi
	if(is.null(pepi)) next;
	l = getLength(pepi); l
	loc = as.numeric(tmp[i, loc]); loc
	if(loc > l) {
		tmp[i, n.seq := 'stop']
		tmp[i, t.seq := 'stop']
		next;
	}
	loc.r = loc + 12; if(loc.r > l) loc.r = l; loc.r
	loc.l = loc - 12; if(loc.l < 1) loc.l = 1; loc.l
	if(tmp[i, Variant_Classification] == 'Missense_Mutation'){
		# get HGVSp
		hgvsp = tmp[i, HGVSp_Short]
		hgvsp = sub('p.', '', hgvsp); hgvsp
		a.l = substr(hgvsp, 1, 1); a.l
		a.r = substr(hgvsp, nchar(hgvsp), nchar(hgvsp)); a.r
		# get 25 aa for normal
		xx = seqinr::getSequence(getFrag(pepi, loc.l, loc.r), as.string=T);xx = xx[[1]];xx
		tmp[i, n.seq := xx]
		# get 25 aa for tumor
		xx.l = seqinr::getSequence(getFrag(pepi, loc.l, loc-1), as.string=T); xx.l = xx.l[[1]]; xx.l
		xx.r = seqinr::getSequence(getFrag(pepi, loc+1, loc.r), as.string=T); xx.r = xx.r[[1]]; xx.r
		xx.m = seqinr::getSequence(getFrag(pepi, loc, loc), as.string=T); xx.m = xx.m[[1]]; xx.m
		if(xx.m == a.l){
		tmp[i, n.seq := xx]
		tmp[i, t.seq := paste0(xx.l, a.r, xx.r)]
		}else{
			tmp[i, n.seq := 'stop']
			tmp[i, t.seq := 'stop']
			#stop(paste0('something wrong i=', i, "\n", paste(tmp[i, 1:9], collapse="\t\n"), "\nleft of HGVSp:", xx.m, "\n"))
		}
	}
}

## neoantigen


## 
import pyGeno.bootstrap as B
B.importGenome("Human.GRCh37.75.tar.gz")






