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

fread('target.tsv') ->target.1
target.1
target.1[, R12 := 'R1']
target.1[grep("R2_001", fastq), R12 := 'R2']
target.1[, samplename := sub(".*\\/(.*)_IGO.*", "\\1", fastq)]
target.1[, samplename := sub("_test", "", samplename)]
target.1[, samplename := paste0("pat-", samplename)]
target.1

target = dcast(target.1, samplename ~ R12, value.var = 'fastq')
target[, patientname := substr(samplename, 1, 7)]
target[, 1:5]


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

### summarize hs 
cmd = paste0('find ', cfg$initDir, ' -name "*metrics_ex*" > hs.flist')
cmd
system(cmd)
hs.flist = fread('hs.flist', header=F)
hs.flist = hs.flist$V1
hs.flist
pat.names = sub(".*(pat-...-..).*", "\\1", hs.flist)
pat.names
hs.dt = lapply(hs.flist, fread, sep="\t")
names(hs.dt) = pat.names
names(hs.dt)

hs.dt = rbindlist(hs.dt, use.names=T, fill=T)
hs.dt[, sample := pat.names]
hs.dt
hs.dt.s = melt(hs.dt, id.vars = 'sample')
hs.dt.s
unique(hs.dt.s$variable)
hs.dt.s.1 = hs.dt.s[variable %in% c('MEAN_TARGET_CONVERAGE', 'PCT_USABLE_BASES_ON_TARGET', 'FOLD_ENRICHMENT', 'PCT_TARGET_BASES_2X', 'PCT_TARGET_BASES_10X', 'PF_UNIQUE_READS'),]
hs.dt.s.1[, sample := factor(sample, levels=sort(unique(sample)))]
hs.dt.s.1$sample
hs.dt.s.1
gg = ggplot(hs.dt.s.1, aes(x=sample, y=value)) + geom_bar(stat='identity') + facet_wrap(~variable, scales='free')
gg = gg + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=.5)) + xlab('') + ylab('')
ggsave(gg, file=paste0(cfg$methylDir, '/summary.pdf'), width=9, height=6)
sync('r_fang/methyl')

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

## its actually mincov 3 as a base line filter
target[, cov.file := paste0(sampledir, '/', samplename, '_pe.deduplicated.bismark.cov.gz')]
target[, methylkit.file := paste0(sampledir, '/', samplename, '_pe.deduplicated.bismark_methylkit.mincov3.txt')]
file.exists(target$cov.file)

target[, {
	cat('now dealing with ', cov.file, ' \n')
	met.cov = fread(paste0("zcat ", cov.file))
	## its actually mincov 3 as a base line filter
	met.cov = met.cov[(V5 + V6) > 5,] 
	met.cov[, V1 := paste0('chr', V1)]
	met.cov = foverlaps(met.cov, target.interval[, 1:3])
	head(met.cov)
	met.cov = met.cov[!is.na(V2),]
	met.cov[, coverage := V5 + V6]
	met.cov[, freqC := round(100 * V5 / coverage, 2)]
	met.cov[, freqT := round(100 * V6 / coverage, 2)]
	met.cov[, chrBase := paste0(V1, '.', i.V2)]
	met.cov[, strand := 'F']
	setnames(met.cov, 'i.V2', 'base')
	setnames(met.cov, 'V1', 'chr')
	met.cov = met.cov[, .(chrBase, chr, base, strand, coverage, freqC, freqT)]
	fwrite(met.cov, file=methylkit.file, sep="\t", quote=F)
}, by=1:nrow(target)]

## decide the strand for bismark coverage results
target$samplename
target[, methylkit.bed.file := paste0(sampledir, '/', samplename, '_pe.deduplicated.bismark_methylkit.mincov3.bed')]
target[, {
	cat(samplename, " now ... \n\t")
	cmd = paste0("awk '{if(NR==1){}else{print($2, $3-1, $3)}}' ", methylkit.file, " | sed \"s/ /\t/g\" | sort -k1 -k2 | sed \"s/chr//\" > ", methylkit.bed.file)
	cat("\t\t\t', cmd, '\n")
	system(cmd)
}, by=1:nrow(target)]

target[, getCG.jobname := paste0(samplename, '.getCG')]
target[, base.gc.file := paste0(sampledir, '/', samplename, '_base_GC.txt')]
target[, getCG.cmd := bsub.head(getCG.jobname, mem=20, cpu=1, We='10:00')]
target[, getCG.cmd := paste0(getCG.cmd, " \"", cfg$bedtools, '/bedtools getfasta -fi ', cfg$genomeGRCh37Fasta, ' -bed ', methylkit.bed.file, ' -fo ', base.gc.file, ' -tab ', "\"")]
target[1, getCG.cmd]
exe.jobs(target$getCG.cmd, logger)

wait4jobs(target$getCG.jobname)

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
target[,.(samplename, group)]
save(target, file='target.RData')
save(cfg, file='cfg.RData')

## use the methylkit.file 

save.image()

library(methylKit)
mincov = 8
head(target$methylkit.file)
methyl.o = methRead(as.list(target$methylkit.file), sample.id=as.list(target$samplename), assembly='GRCh37', treatment=factor(target$group), context='CpG')
## filter
methyl.fil.o = filterByCoverage(methyl.o, lo.count=mincov, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
methyl.fil.o
#rm(methyl.o)
meth.raw = methyl.fil.o
meth.raw
meth.raw.file = paste0(cfg$methylDir, '/methyl_raw_fil_obj_mincov', mincov, '.RData')
meth.raw.file
save(meth.raw, file=meth.raw.file)
rm(methyl.o, methyl.fil.o)

## import this raw data into the raw sqlite
library(RSQLite)
library(DBI)

# raw.cov10.hiC99.
con = dbConnect(drv=RSQLite::SQLite(), dbname=paste0(cfg$methylDir, '/sqlite/methyl.db'))
dbSendStatement(con, 'drop table if exists methylRaw')
target[, {dt.m = fread(methylkit.file)
       dt.m = dt.m[coverage > mincov & freqC < 99.9,]
       dt.m[, samplename := samplename]
       dbWriteTable(con, 'methylRaw', as.data.frame(dt.m), append=T)
}, by=1:nrow(target)]

dbSendStatement(con, 'create index iidnex on methylRaw (samplename, chr, base)')
#dbSendStatement(con, 'drop table methylRaw')
#dbGetQuery(con, paste0('select * from methylRaw limit 10'))

#dbDisconnet(con)

## combine, destrand for both strand when analyzing CpG
# ======================================================
## based on the manifest annotation file from illumina
meth.raw.file
#save(meth.fil.o, file=methyl.raw.file)

meth.raw.file = meth.raw.file
region.file = cfg$hsaMethylAnno
gene.obj.file = system.file('extdata', 'hg19_gencode_genesymbol.bed.txt', package = 'methylKit')
output.file = paste0(cfg$cwd, '/meth_bait')
pdf.file = paste0(cfg$methylDir, '/cpg_bait')
dm.bait.jobname = 'meth.bait.diff'
dm.bait.cmd = bsub.head(cpu=10, mem=20, jobname=dm.bait.jobname, We='10:30')
dm.bait.cmd = paste0(dm.bait.cmd, " \" /home/huw/local/bin/Rscript --vanilla /home/huw/program/fun/my.calculate.diff.meth.r ") 
dm.bait.cmd = paste0(dm.bait.cmd, ' --meth-raw-file ', meth.raw.file, ' --region-file ', region.file, ' --gene-obj-file ', gene.obj.file, ' --output-file ', output.file, " --pdf-file ", pdf.file, "\"")
dm.bait.cmd

exe.jobs(dm.bait.cmd, logger)

# ======================================================
## based on the individual CpG, no region-file here
gene.obj.file = system.file('extdata', 'hg19_gencode_genesymbol.bed.txt', package = 'methylKit')
output.file = paste0(cfg$cwd, '/meth_cpg')
pdf.file = paste0(cfg$methylDir, '/cpg')
dm.jobname = 'meth.cpg.diff'
dm.cmd = bsub.head(cpu=10, mem=20, jobname=dm.jobname, We='10:30')
dm.cmd = paste0(dm.cmd, "\" /home/huw/local/bin/Rscript --vanilla /home/huw/program/fun/my.calculate.diff.meth.r ") 
dm.cmd = paste0(dm.cmd, ' --meth-raw-file ', meth.raw.file, ' --gene-obj-file ', gene.obj.file, ' --output-file ', output.file, " --pdf-file ", pdf.file, "\"")
dm.cmd
exe.jobs(dm.cmd, logger)

# ======================================================
## based on the promoter regions 
source('~/program/fun/granges_2_bed.r')
gene.obj=readTranscriptFeatures(gene.obj.file, up.flank=2000, down.flank=2000)
granges.2.bed(gene.obj$promoters, 'promoter.bed')
rm(gene.obj)

region.file = normalizePath('promoter.bed')
region.file
gene.obj.file = system.file('extdata', 'hg19_gencode_genesymbol.bed.txt', package = 'methylKit')
output.file = paste0(cfg$cwd, '/meth_tss')
pdf.file = paste0(cfg$methylDir, '/cpg_tss')
dm.tss.jobname = 'meth.tss.diff'
dm.tss.cmd = bsub.head(cpu=10, mem=20, jobname=dm.tss.jobname, We='10:00')
dm.tss.cmd = paste0(dm.tss.cmd, "\" /home/huw/local/bin/Rscript --vanilla /home/huw/program/fun/my.calculate.diff.meth.r ") 
dm.tss.cmd = paste0(dm.tss.cmd, ' --meth-raw-file ', meth.raw.file, ' --region-file ', region.file, ' --gene-obj-file ', gene.obj.file, ' --output-file ', output.file, " --pdf-file ", pdf.file, "\"")
dm.tss.cmd

exe.jobs(dm.bait.cmd, logger)
exe.jobs(dm.tss.cmd, logger)
exe.jobs(dm.cmd, logger)

# the following files will be created
# save(meth.meth, file = paste0(opt$output.file, '_meth_unite.RData')) # united counts
# save(meth.diff, file = paste0(opt$output.file, '_meth_diff.RData')) # DM
# save(meth.diff.anno, file = paste0(opt$output.file, '_meth_diff_anno.RData')) # annotated DM
# save(meth.diff.hyper.anno, file = paste0(opt$output.file, '_meth_diff_hyper_anno.RData')) # hyper DM 
# save(meth.diff.hypo.anno, file = paste0(opt$output.file, '_meth_diff_hypo_anno.RData')) # hypo DM 
# save(meth.diff.anno.df, file = paste0(opt$output.file, '_meth_diff_anno_df.RData')) # a combination of annotation and meth.diff
# pdf.file3 = paste0(opt$pdf.file, "_anno_summary.pdf")	# summary file
# pdf.file2 = paste0(opt$pdf.file, '_summary.pdf') # summary file
# pdf.file = paste0(opt$pdf.file, '_MAplot.pdf') # MA plot

## pca analysis
load('meth_cpg_meth_unite.RData')
pca = getData(meth.meth)
as.data.table(pca) -> pca
pca[, .(numCs1/coverage1, numCs2/coverage2, numCs3/coverage3, numCs4/coverage4, numCs5/coverage5, numCs6/coverage6)] -> pca.2
pca.2[rowSums(is.na(pca.2)) == 0,] -> pca.3
source('~/program/fun/rowSds.r')
rowSds(pca.3) -> tmp
pca.3[tmp > 0.1, ] -> pca.4
dim(pca.4)
prcomp(pca.4) -> pca.5
pca.5[[2]] -> pca.6
row.names(pca.6) = meth.meth@sample.ids
pca.7 = as.data.table(pca.6, keep.rownames=T)
pca.7[, t12 := as.numeric(substr(rn, 10, 10))]
gg = ggplot(pca.7, aes(x = PC1, y = PC2, label=rn, color=t12)) + geom_point()
gg = gg + geom_text_repel() + xlab(paste0('PC1 (', round(100*pca.5[[1]][1],2), '%)'))
gg = gg + ylab(paste0('PC2 (', round(100*pca.5[[1]][2],2), '%)'))
gg = gg + theme(legend.position = "none") 
ggsave(gg, file=paste0(cfg$methylDir, '/pca_ggplot.pdf'), width=4, height=4)
sync('r_fang/methyl')

## import CpG methy.diff to the database
con = dbConnect(drv=RSQLite::SQLite(), dbname=paste0(cfg$methylDir, '/sqlite/methyl.db'))
dbSendStatement(con, 'drop table if exists methylCpGDiff')
x = load(paste0(cfg$cwd, '/meth_cpg_meth_diff.RData'))
x
system(paste0(' ls -lh ', cfg$cwd, '/meth_cpg_meth_diff.RData'))
tmp = getData(meth.diff)
head(tmp)
dbWriteTable(con, 'methylCpGDiff', tmp, append=F)
dbSendStatement(con, 'create index iCpGDiff on methylCpGDiff (chr, start)')
#dbDisconnect(con)
dbGetQuery(con, paste0('select * from methylCpGDiff limit 10'))
dbGetQuery(con, paste0('select count(*) from methylCpGDiff limit 10'))

## import Illumina annotated CpG island methy.diff to the database
dbSendStatement(con, 'drop table if exists methylBaitDiff')
load('meth_bait_meth_diff.RData')
tmp = getData(meth.diff)
dbWriteTable(con, 'methylBaitDiff', tmp, append=F)
dbSendStatement(con, 'create index iBaitDiff on methylBaitDiff (chr, start)')
dbGetQuery(con, paste0('select * from methylBaitDiff limit 10'))
#dbSendStatement(con, 'drop table methylRaw')
#dbGetQuery(con, paste0('select * from methylRaw limit 10'))
#dbDisconnet(con)

# focused on CpG and illumina annotated CpG island
output.file = paste0(cfg$cwd, '/meth_bait')
xx = load(paste0(output.file, '_meth_diff_anno_df.RData'))
xx = load(paste0(output.file, '_meth_diff_anno.RData'))
xx = load(paste0(output.file, '_meth_diff.RData'))

## specific gene region plot, parameters
gene.obj.file = system.file('extdata', 'hg19_gencode_genesymbol.bed.txt', package = 'methylKit')
cpg.meth.Rfile = normalizePath(paste0('./meth_cpg_meth_diff.RData'))
system(paste0('ls -lh ',  cpg.meth.Rfile))
cpgi.meth.Rfile = normalizePath(paste0('./meth_bait_meth_diff.RData'))
output.dir = cfg$methylDir
cpgi.anno.file = cfg$hsaMethylAnno

gene.symbol = 'KRT6A'; chr='chr12'; chr.start=52875224; chr.end=52896391
cmd = paste0('Rscript ~/program/fun/plot.meth.multi.r --cpgi-anno-file ', cpgi.anno.file, ' --gene-anno-file ', gene.obj.file, ' --cpgi-meth-Rfile ', cpgi.meth.Rfile, ' --cpg-meth-Rfile ', cpg.meth.Rfile, ' --output-dir ', output.dir, ' --gene-symbol ', gene.symbol, ' --chr ', chr, '  --chr-start ', chr.start, ' --chr-end ', chr.end, ' --width ', 10, ' --height ', 11, ' --cpgi-ucsc-anno-file ', cfg$hg19.cpgi, ' --meth-db ', cfg$methylDir)
cmd
exe.jobs(cmd, logger)
sync('r_fang/methyl')

## GATA3: chr10 8091374 8098329
gene.symbol = 'GATA3'; chr = 'chr10'; chr.start = 8096000; chr.end = 8105329
cmd = paste0('Rscript ~/program/fun/plot.meth.multi.r --cpgi-anno-file ', cpgi.anno.file, ' --gene-anno-file ', gene.obj.file, ' --cpgi-meth-Rfile ', cpgi.meth.Rfile, ' --cpg-meth-Rfile ', cpg.meth.Rfile, ' --output-dir ', output.dir, ' --gene-symbol ', gene.symbol, ' --chr ', chr, '  --chr-start ', chr.start, ' --chr-end ', chr.end, ' --width ', 10, ' --height ', 11, ' --cpgi-ucsc-anno-file ', cfg$hg19.cpgi, ' --meth-db ', cfg$methylDir)
exe.jobs(cmd, logger)
sync('r_fang/methyl')

# FOXA 138,057,364-38,065,717 
gene.symbol = 'FOXA1'; chr = 'chr14'; chr.start = 38057000; chr.end = 38065768
cmd = paste0('Rscript ~/program/fun/plot.meth.multi.r --cpgi-anno-file ', cpgi.anno.file, ' --gene-anno-file ', gene.obj.file, ' --cpgi-meth-Rfile ', cpgi.meth.Rfile, ' --cpg-meth-Rfile ', cpg.meth.Rfile, ' --output-dir ', output.dir, ' --gene-symbol ', gene.symbol, ' --chr ', chr, '  --chr-start ', chr.start, ' --chr-end ', chr.end, ' --width ', 10, ' --height ', 11, ' --cpgi-ucsc-anno-file ', cfg$hg19.cpgi, ' --meth-db ', cfg$methylDir)
exe.jobs(cmd, logger)
sync('r_fang/methyl')

## PPARG
ss = 'chr3:12,310,146-12,375,710'
ss = gsub(",", "", ss); ss
chr = sub("(.*):.*-.*", "\\1", ss); chr
chr.start = sub(".*:(.*)-.*", "\\1", ss); chr.start
chr.end = sub(".*:.*-(.*)", "\\1", ss); chr.end
gene.symbol = 'PPARG'
cmd = paste0('Rscript ~/program/fun/plot.meth.multi.r --cpgi-anno-file ', cpgi.anno.file, ' --gene-anno-file ', gene.obj.file, ' --cpgi-meth-Rfile ', cpgi.meth.Rfile, ' --cpg-meth-Rfile ', cpg.meth.Rfile, ' --output-dir ', output.dir, ' --gene-symbol ', gene.symbol, ' --chr ', chr, '  --chr-start ', chr.start, ' --chr-end ', chr.end, ' --width ', 10, ' --height ', 11, ' --cpgi-ucsc-anno-file ', cfg$hg19.cpgi, ' --meth-db ', cfg$methylDir)
exe.jobs(cmd, logger)
sync('r_fang/methyl')

## 
ss = 'chr12:49,403,615-49,458,864'; 
gene.symbol = 'KMT2D'
ss = gsub(",", "", ss); ss
chr = sub("(.*):.*-.*", "\\1", ss); chr
chr.start = sub(".*:(.*)-.*", "\\1", ss); chr.start
chr.end = sub(".*:.*-(.*)", "\\1", ss); chr.end
cmd = paste0('Rscript ~/program/fun/plot.meth.multi.r --cpgi-anno-file ', cpgi.anno.file, ' --gene-anno-file ', gene.obj.file, ' --cpgi-meth-Rfile ', cpgi.meth.Rfile, ' --cpg-meth-Rfile ', cpg.meth.Rfile, ' --output-dir ', output.dir, ' --gene-symbol ', gene.symbol, ' --chr ', chr, '  --chr-start ', chr.start, ' --chr-end ', chr.end, ' --width ', 10, ' --height ', 11, ' --cpgi-ucsc-anno-file ', cfg$hg19.cpgi, ' --meth-db ', cfg$methylDir)
exe.jobs(cmd, logger)
sync('r_fang/methyl')


library(methylKit)

pdf(file=paste0(cfg$methylDir, '/mix_model_cost.pdf'))
myMixdml = myDiff.to.mixmdl(myDiff, plot=T, main='SCC T2 vsT1')
plotCost(myMixmdl, main='cost function')
dev.off()

save.image()

## summary hs metrics
target$hs.metrics.file
target[, hs.metrics.ex.file := sub(".txt", "_ex.txt", hs.metrics.file)]
target[, {
	cmd = paste0('sed -n 7,8p ', hs.metrics.file, ' > ', hs.metrics.ex.file)
	system(cmd)
}, by = 1:nrow(target)]

hs.dt = lapply(target$hs.metrics.ex.file, fread)
names(hs.dt) = target$samplename
hs.dt
dd = t(hs.dt[[1]][, 30:54])
for(i in 2:length(hs.dt)){
	dd = cbind(dd, t(hs.dt[[1]][, 30:54]))
}
colnames(dd) = names(hs.dt)
dd
dd = as.data.table(dd, keep.rownames=T)
dd = dd[rn == 'PCT_TARGET_BASES_10X',]
dd = melt(dd)
dd
gg = ggplot(dd, aes(x = variable, y = value)) + geom_bar(stat='identity') + 
	theme(axis.text.x = element_text(angle=90, hjust=0.5)) +
	ylim(0,1) + xlab('') + ylab('% 10X target coverage')
ggsave(gg, file=paste0(cfg$methylDir, '/pct_target_10x_coverage.pdf'), width=3, height=3)
sync('r_fang/methyl')

## import the RNA-Seq results
library(readxl)
res = as.data.table(read_excel('../blca_cmo_06155_2016/rnaseq/SCC.T2vsT1.DESeq2.DEG.xlsx'))
head(res)
res[, sig := 'not']
res[log2FoldChange > 1 & padj < 0.05, sig := 'up']
res[log2FoldChange < -1 & padj < 0.05, sig := 'dn']
head(res)
res[Gene == 'RNASE7',]
res[Gene == 'FOXA1',]
res[sig == 'up', Gene]
meth.diff.anno.df[symbol == 'GATA3',]

gene.obj.file = system.file('extdata', 'hg19_gencode_genesymbol.bed.txt', package = 'methylKit')
cpg.meth.Rfile = normalizePath('meth_cpg_meth_diff.RData')
cpgi.meth.Rfile = normalizePath('meth_bait_meth_diff.RData')
cpgi.anno.file = cfg$hsaMethylAnno
xx = load('meth_bait_meth_diff_anno_df.RData'); xx
meth.diff.anno.df[symbol %in% res[sig == 'up', Gene], gep := 'up']
meth.diff.anno.df[symbol %in% res[sig == 'dn', Gene], gep := 'dn']
tmp = meth.diff.anno.df[!is.na(gep) & abs(dist.to.feature) < 10000 & qvalue < 0.05 & abs(meth.diff) > 25,]
tmp[symbol == 'GATA3', ]
tmp[, start2 := start - abs(dist.to.feature) - 2000]
tmp[, end2 := end + abs(dist.to.feature) + 2000]
tmp[, start.old := start]
tmp[, end.old := end]
tmp[, start := start2]
tmp[, end := end2]
tmp[meth.diff > 25 & gep == 'dn', output := paste0(cfg$methylDir, '/hyperDM_DN')]
tmp[meth.diff > 25 & gep == 'up', output := paste0(cfg$methylDir, '/hyperDM_UP')]
tmp[meth.diff < -25 & gep == 'dn', output := paste0(cfg$methylDir, '/hypoDM_DN')]
tmp[meth.diff < -25 & gep == 'up', output := paste0(cfg$methylDir, '/hypoDM_UP')]
tmp[meth.diff < -25, DM := 'hypo']
tmp[meth.diff > 25, DM := 'hyper']
target.bed.file = paste0(cfg$methylDir, '/GEP_DM.bed')
fwrite(tmp, file=target.bed.file, sep="\t", quote=F)

for(i in unique(tmp$output)){
	unlink(i, force=T, recursive=T)
}

for(i in unique(tmp$output)){ dir.create(i) }

tmp
output.dir = paste0(cfg$methylDir, '/DM_GEP_plot'); output.dir
cmd = paste0('Rscript ~/program/fun/plot.meth.multi.r --cpgi-anno-file ', cpgi.anno.file, ' --gene-anno-file ', gene.obj.file, ' --cpgi-meth-Rfile ', cpgi.meth.Rfile, ' --cpg-meth-Rfile ', cpg.meth.Rfile, ' --output-dir ', output.dir, ' --width ', 10, ' --height ', 11, ' --cpgi-ucsc-anno-file ', cfg$hg19.cpgi, ' --meth-db ', cfg$methylDir, ' --target-bed ', target.bed.file)
cmd
exe.jobs(cmd, logger)
sync('r_fang/methyl')

meth.diff.anno.df[gep == 'up' & abs(dist.to.feature) < 4000 & qvalue < 0.05 & meth.diff < -25,]
meth.diff.anno.df[gep == 'dn' & abs(dist.to.feature) < 4000 & qvalue < 0.05 & meth.diff > 25,]

xx = load('meth_bait_meth_diff_anno_df.RData')
meth.diff.anno.df[, DM := 'N.S.']
meth.diff.anno.df[qvalue < 0.001 & meth.diff < -25, DM := 'Hypo']
meth.diff.anno.df[qvalue < 0.001 & meth.diff >  25, DM := 'Hyper']
tmp = table(meth.diff.anno.df$tss, meth.diff.anno.df$DM)
tmp = melt(tmp)
tmp
gg = ggplot(tmp, aes(x=Var2, y=value)) + geom_bar(stat='identity') + facet_wrap(~Var1) + 
	xlab('')  + ylab('Differetially Methylated CpG islands \nannotated in Illumina manifest file')
ggsave(gg, file=paste0(cfg$methylDir, '/DM_cat.pdf'), width=3, height=4)
sync('r_fang/methyl')

meth.diff.anno.df[meth.diff>25 & qvalue < 0.001, DM := 'hyper']
meth.diff.anno.df[meth.diff<-25 & qvalue < 0.001, DM := 'hypo']
table(meth.diff.anno.df[, .(DM, gep, tss)])

meth.diff.anno.df
xx = load('meth_bait_meth_unite.RData')
meth.meth
dm.gr = as(meth.diff.anno.df[abs(meth.diff > 25) & qvalue < 0.001, ], 'GRanges')
findOverlaps(as(meth.meth, 'GRanges'), dm.gr)

## END
### https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
