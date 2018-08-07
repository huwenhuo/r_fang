library(data.table)
library(log4r)

resDir = 'r_fang'
inPre = 'Proj_07813_DF'
pre = 'scc'
species = 'b37'
assay = 'dnaseq'

## copy configure.R configure_dir.R into resDir
source(paste0(resDir, '/configure.R'))

## to generate pon Panel of Normal
pon = fread("r_fang/bamfile", header=F)
colnames(pon) = 'bamfile'
pon[, base := sub(".*recal_", "", bamfile)]
pon[, base := sub(".bam", "", base)]
pon[, vcffile := paste0(resDir, '/pon/', base, '.vcf.gz')]
pon[, pon.jobname := paste0('pon.', base)]
pon[, cmd := bsub.head(jobname = pon.jobname, mem=30, cpu = 5)]
pon[, cmd := paste0(cmd, " \"", GATK4, " --java-options '-Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' Mutect2")]
pon[, cmd := paste0(cmd, " -R ", genomeFasta, " -I ", bamfile, " -tumor ", base, " --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -O ", vcffile, " \"")]
pon$cmd[1]

exe.jobs(pon$cmd, logger)

## generate pon file, panel of normal
ponfile = paste0(resDir, '/pon/pon56Sample.vcf.gz')
cmd = bsub.head(jobname ='pon', , mem=30, cpu = 5)
cmd = paste0(cmd, " \"", GATK4, " --java-options '-Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' CreateSomaticPanelOfNormals ")
cmd = paste0(cmd, paste(paste0(" -vcfs ", pon$vcffile), collapse = " "))
cmd = paste0(cmd, " -O ", ponfile, "\"")
cmd
write.table(cmd, file='pon.sh', quote=F, sep="\t")

exe.jobs(cmd, logger)

### start here
targets = fread(mapfile, header = F)
## group == library
targets = setNames(targets, c('group', 'sampleName', 'runID', 'fastqDir', 'seqType')) # 
targets
targets.aln = copy(targets)
targets.aln[, fastq1 := '']
targets.aln[, fastq2 := '']
targets.aln[, fastq1.full := '']
targets.aln[, fastq2.full := '']
targets.aln = targets.aln[0,]
targets.aln

for(i in 1:nrow(targets)){
	tmp = copy(targets[i,])
	tmp
	fastqs
	fastqs = dir(path=tmp$fastqDir, pattern="*.fastq.gz")
	fastq1 = fastqs[grep("_R1_", fastqs)]
	fastq1
	tmp = tmp[rep(1,length(fastq1)),]
	tmp
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
##my $r1adaptor = 'AGATCGGAAGAGCACACGTCT';
##my $r2adaptor = 'AGATCGGAAGAGCGTCGTGTA';
cwd = getwd()

## r_fang/sampleName_group_runID as file folder for each 
targets.aln[, sgrTag := paste0(sampleName, '_', group, '_', runID)]
targets.aln[, sgrDir := paste0(resDir, '/', sgrTag)]
targets.aln$sgrDir[1]
targets.aln$sgrDir

for(i in 1:nrow(targets.aln)){system(paste0('mkdir -p ', targets.aln$sgrDir[i]))}

## zcat fastq
targets.aln[, fastq1.zcatfile := file.path(sgrDir, sub(".gz", "", fastq1))]
targets.aln[, fastq2.zcatfile := file.path(sgrDir, sub(".gz", "", fastq2))]

targets.aln[, zcat.jobname := paste0("zcat.", sampleName, '.', fastq1)]
targets.aln[, zcat.cmd1 := bsub.head(zcat.jobname, mem=1, cpu=1, We='6:59', cwd)]
targets.aln[, zcat.cmd1:= paste0(zcat.cmd1, ' "/bin/zcat ', fastq1.full, ' > ', fastq1.zcatfile, '"')]
targets.aln$zcat.cmd1[1]

targets.aln[, zcat.jobname := paste0("zcat.", sampleName, '.', fastq2)]
targets.aln[, zcat.cmd2 := bsub.head(zcat.jobname, mem=1, cpu=1, We='6:59', cwd)]
targets.aln[, zcat.cmd2 := paste0(zcat.cmd2, ' "/bin/zcat ', fastq2.full, ' > ', fastq2.zcatfile, '"')]
targets.aln$zcat.cmd2[1]

exe.jobs(targets.aln$zcat.cmd1, logger)
exe.jobs(targets.aln$zcat.cmd2, logger)

## convert quality score
targets.aln[, fastq1.cqsfile := file.path(sgrDir, sub(".gz", ".cqs", fastq1))]
targets.aln[, fastq2.cqsfile := file.path(sgrDir, sub(".gz", ".cqs", fastq2))]
targets.aln$fastq1.cqsfile[1]

targets.aln[, cqs.jobname := paste0("cqs.", sampleName, '.', fastq1)]
targets.aln[, cqs.cmd1 := bsub.head(cqs.jobname, mem=1, cpu=1, We='2:59', cwd, zcat.jobname)]
targets.aln[, cqs.cmd1 := paste0(cqs.cmd1, ' "', ConvertQualityScore, ' --input ', fastq1.zcatfile, ' --output ', fastq1.cqsfile, '"')]
targets.aln$cqs.cmd1[1]

targets.aln[, cqs.jobname := paste0("cqs.", sampleName, '.', fastq1)]
targets.aln[, cqs.cmd2 := bsub.head(cqs.jobname, mem=1, cpu=1, We='2:59', cwd, zcat.jobname)]
targets.aln[, cqs.cmd2 := paste0(cqs.cmd2, ' "', ConvertQualityScore, ' --input ', fastq2.zcatfile, ' --output ', fastq2.cqsfile, '"')]
targets.aln$cqs.cmd2[1]

for(i in 1:nrow(targets.aln)){system(targets.aln$cqs.cmd1[i])}
for(i in 1:nrow(targets.aln)){system(targets.aln$cqs.cmd2[i])}

## discard read less than half length
read.len = function(x){
	gzfile(x, open='r') -> fh
	readLines(fh, 2) -> h
	nchar(h[2])
}
targets.aln[, read.len := sapply(targets.aln$fastq1.full, read.len)]
targets.aln$read.len

targets.aln[, fastq1.cutadaptfile :=  file.path(sgrDir, sub(".gz", ".cutadapt", fastq1))]
targets.aln[, fastq2.cutadaptfile :=  file.path(sgrDir, sub(".gz", ".cutadapt", fastq2))]
targets.aln$fastq1.cutadaptfile[1]

targets.aln[, cutadapt.jobname := paste0("cutadapt.", sgrTag)]
targets.aln[, cutadapt.cmd := bsub.head(cutadapt.jobname, mem=1, cpu=1, We='2:59', cwd, postdone = cqs.jobname)]
targets.aln[, cutadapt.cmd := paste(cutadapt.cmd, ' "', PYTHON, CUTADAPT, '-f fastq -a', clipR1, '-A', clipR2, '--overlap 10 --minimum-length', floor(read.len/2))]
targets.aln[, cutadapt.cmd := paste(cutadapt.cmd, '--quality-cutoff', bqtrim, '-o', fastq1.cutadaptfile, '--paired-output', fastq2.cutadaptfile, fastq1.cqsfile, fastq2.cqsfile, ' "')]
targets.aln$cutadapt.cmd[1]

#for(i in 1:nrow(targets.aln)){system(targets.aln$cutadapt.cmd[i])}
exe.jobs(targets.aln$cutadapt.cmd, logger)

## breakpoint
wait4jobs(targets.aln$cutadapt.jobname, logger)
testOutputFiles = c(targets.aln$fastq1.cutadapt, targets.aln$fastq2.cutadapt)
if(delFiles){
	del.files(targets.aln$fastq1.zcatfile, logger)
	del.files(targets.aln$fastq2.zcatfile, logger)
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
targets.aln[, bwa.cmd := paste0(bwa.cmd, ' "', bwa, ' mem -t 8 -T 0 -R \\\"', readgroup, '\\\" ', reference, ' ', fastq1.cutadaptfile, ' ', fastq2.cutadaptfile, ' | ', samtools, ' view -Shb -o ', bamfile, ' - "')]
targets.aln$bwa.cmd[1]

for(i in 1:nrow(targets.aln)){system(targets.aln$bwa.cmd[i])}
exe.jobs(targets.aln$bwa.cmd, logger)

## breakpoint
wait4jobs(targets.aln$cutadapt.jobname, logger)

## clean files
testInputFiles = c(targets.aln$fastq1.cutadapt, targets.aln$fastq2.cutadapt)
testOutputFiles = targets.aln$bamFile
if(delFiles & all(file.size(testOutputFiles) > 2000000)){
	del.files(targets.aln$fastq1.cutadapt, logger));
	del.files(targets.aln$fastq2.cutadapt, logger));
}

## merge bam file for the same library/group 
## may several r1/r2 fastq files
tmpdir='/scratch/huw/'
targets.bam = targets.aln[, .(bamFileList = paste(paste0(' I=', bamfile), collapse = ' '), bwaJobList = paste(bwa.jobname, collapse = ' ')), by = list(sampleName, group)]
targets.bam
targets.bam$bwaJobList
targets.bam[, libMergedBamFile := paste0(resDir, '/', sampleName, '_library', group, '_merged.bam')]
targets.bam[, merge.jobname := paste0('mergebam.', sampleName, group)]
targets.bam$merge.jobname
targets.bam[, merge.cmd := bsub.head(merge.jobname, mem=30, cpu=8, We='2:26', cwd=cwd, postdone = bwaJobList)]
targets.bam[, merge.cmd := paste0(merge.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar MergeSamFiles ', bamFileList, ' O=', libMergedBamFile, ' SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=', tmpdir, ' CREATE_INDEX=TRUE USE_THREADING=FALSE MAX_RECORDS_IN_RAM=5000000 "')] 
targets.bam$merge.cmd[2]

#for(i in 1:nrow(targets.bam)){system(targets.bam$merge.cmd[i])}
exe.jobs(targets.bam$merge.cmd, logger)

## markDuplicates for reads from the same library/group
## merge bam file from the same library/group
## output in resDir (r_fang)
tmpdir='/scratch/huw/'
targets.bam[, mkdupedLibBamFile := paste0(resDir, '/', sampleName, '_library', group, '_merged_mkdup.bam')]
targets.bam[, libmkdupMetricsFile := paste0(matricsDir, '/', sampleName, '_library', group, '_merged_mkdup_metrics.txt')]
targets.bam[, mkdup.jobname := paste0('mkdup.', sampleName, '.', group)]
targets.bam$mkdup.jobname
targets.bam[, mkdup.cmd := bsub.head(mkdup.jobname, mem=30, cpu=3, We='8:26', cwd=cwd, postdone = merge.jobname)]
targets.bam[, mkdup.cmd := paste0(mkdup.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar MarkDuplicates I=', libMergedBamFile, ' O=', mkdupedLibBamFile, ' METRICS_FILE=', libmkdupMetricsFile, ' VALIDATION_STRINGENCY=LENIENT TMP_DIR=', tmpdir, ' REMOVE_DUPLICATES=TRUE CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=5000000 "')] 
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
targets.smpBam[, smp.cmd := paste0(smp.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar MergeSamFiles ', BamFileList, ' TMP_DIR=', tmpdir, ' DIR=', tmpdir, ' O=', smpBamFile, ' SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=', tmpdir, ' CREATE_INDEX=TRUE USE_THREADING=FALSE MAX_RECORDS_IN_RAM=5000000 "')] 
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

targets.pair[, bicNormalBamFile := paste0('Proj_07813_DF/r_001/alignments/Proj_07813_DF_indelRealigned_recal_s_', gsub("-", "_", normal), ".bam") ]
targets.pair[, bicTumorBamFile := paste0('Proj_07813_DF/r_001/alignments/Proj_07813_DF_indelRealigned_recal_s_', gsub("-", "_", tumor), ".bam") ]
file.exists(targets.pair$bicNormalBamFile)
file.exists(targets.pair$bicTumorBamFile)

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
## examine the mutatiosn in HLA
targets.pair[, hla.polysolver2OutDir  := paste0(initDir, '/', tumor, '/polysolver')] ; 
targets.pair[, hla.polysolver2.jobname := paste0('hla.polysolver2.', tumor)]
targets.pair[, hla.polysolver2.cmd := bsub.head(hla.polysolver2.jobname, mem=30, cpu=3, We='8:26', cwd=cwd, postdone=hla.polysolver.jobname), by=1:nrow(targets.pair)] 
targets.pair[, hla.polysolver2.cmd := paste0(hla.polysolver2.cmd, " \" bash ", polysolver, "/scripts/shell_call_hla_mutations_from_type")]
targets.pair[, hla.polysolver2.cmd := paste0(hla.polysolver2.cmd, " ", bicNormalBamFile, " ", bicTumorBamFile, " ", hla.polysolverOutDir, "/winners.hla.txt hg19 STDFQ")] 
targets.pair[, hla.polysolver2.cmd := paste0(hla.polysolver2.cmd, " ", hla.polysolver2OutDir, " \"")]
targets.pair$hla.polysolver2.cmd[3]

exe.jobs(targets.pair$hla.polysolver2.cmd, logger)

## polysolver step 3
## annotate the mutations in the HLA
targets.pair[, hla.polysolver3.jobname := paste0('hla.polysolver3.', tumor)]
targets.pair[, hla.polysolver3.cmd := bsub.head(hla.polysolver3.jobname, mem=30, cpu=3, We='8:26', cwd=cwd, postdone=hla.polysolver2.jobname), by=1:nrow(targets.pair)] 
targets.pair[, hla.polysolver3.cmd := paste0(hla.polysolver3.cmd, " \" bash ", polysolver, "/scripts/shell_annotate_hla_mutations indiv ")]
targets.pair[, hla.polysolver3.cmd := paste0(hla.polysolver3.cmd, " ", hla.polysolver2OutDir, " \"")]
targets.pair$hla.polysolver3.cmd[3]

exe.jobs(targets.pair$hla.polysolver3.cmd, logger)

## ## examine HLA by SOAP-HLA
## polysolver
## targets.pair[, soapHlaOutDir  := paste0(initDir, '/', normal, '/soapHla')] ; 
## targets.pair[, {system(paste0('mkdir -p ', soapHlaOutDir))}, by = 1:nrow(targets.pair)]
## 
## targets.pair[, soaphla.jobname := paste0('soaphla.', normal)]
## targets.pair[, soaphla.cmd := bsub.head(soaphla.jobname, mem=30, cpu=3, We='8:26', cwd=cwd)] 
## targets.pair[, soaphla.cmd := paste0(soaphla.cmd, " \" ", PERL, '/perl ', SOAPHLA, "/MHC_autopipeline.pl -i")]
## targets.pair[, soaphla.cmd := paste0(soaphla.cmd, " ", normalBamFile, " -od ", soapHlaOutDir, " -v hg19 \"")] 
## targets.pair$soaphla.cmd[1]
## cmd = targets.pair$soaphla.cmd[!duplicated(targets.pair$normal)]
## cmd[2]
## 
## exe.jobs(cmd, logger)

## maf2vcf.pl to get vcf files
neoDir = paste0(resDir, '/neo')
if(!dir.exists(neoDir)) {dir.create(neoDir)}
maffile = 'scc_oncokb.maf'
maf2vcf.jobname = 'maf2vcf'
cmd = bsub.head(jobname =maf2vcf.jobname, mem=30, cpu = 5)
cmd = paste0(cmd, " \"", PERL, "/perl ",  VCF2MAF, "/maf2vcf.pl --output-dir ", neoDir, " --input-maf ", maffile, " \"")
cmd

exe.jobs(cmd, logger)

## split vcf per sample 
vcffile
targets.pair[, bicSampleName := paste0('s_', gsub("-", "_", tumor))]
targets.pair[, mafVcfFile := paste0(initDir, '/', tumor, '/neo_', tumor, '.vcf')]
targets.pair[, splitmafvcf.jobname := paste0("splitmafvcf.", tumor)]
targets.pair[, splitmafvcf.cmd := bsub.head(splitmafvcf.jobname, mem=10, We='1:11', cpu=1, postdone = maf2vcf.jobname), by=1:nrow(targets.pair)]
targets.pair[, splitmafvcf.cmd := paste0(splitmafvcf.cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' SelectVariants")]
targets.pair[, splitmafvcf.cmd := paste0(splitmafvcf.cmd, " -R ", genomeFasta, ' -V ', vcfFile, " --sample-name ", bicSampleName, " -O ", mafVcfFile, "\"")]
targets.pair[1, splitmafvcf.cmd]

exe.jobs(targets.pair$splitmafvcf.cmd, logger)

## use the vcf file from mutect from BIC
targets.pair[, bicMutectVcfFile := paste0('Proj_07813_DF/r_001/variants/snpsIndels/mutect/Proj_07813_DF_s_', gsub("-", "_", normal), "_s_", gsub("-", "_", tumor), '_mutect_calls.vcf')]
file.exists(targets.pair$bicMutectVcfFile)

## filter vcf per sample
targets.pair[, filterMutectVcfFile := paste0(initDir, '/', tumor, '/neo_filter_', tumor, '.vcf')]
targets.pair[, filtervcf.jobname := paste0("filtervcf.", tumor)]
targets.pair[, filtervcf.cmd := bsub.head(filtervcf.jobname, mem=20, We='1:11', cpu=5)]
targets.pair[, filtervcf.cmd := paste0(filtervcf.cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' VariantFiltration")]
targets.pair[, filtervcf.cmd := paste0(filtervcf.cmd, " --filter-expression 'AD > 3' --filter-name AD3 ")]
targets.pair[, filtervcf.cmd := paste0(filtervcf.cmd, " -R ", genomeFasta, " -O ", filterMutectVcfFile, " -V ", bicMutectVcfFile, " \"")]
targets.pair[1, filtervcf.cmd]

exe.jobs(targets.pair$filtervcf.cmd, logger)

## exclude filtered in vcf
targets.pair[, exFilterMutectVcfFile := paste0(initDir, '/', tumor, '/neo_ex_filter_', tumor, '.vcf')]
targets.pair[, exFiltervcf.jobname := paste0("exFiltervcf.", tumor)]
targets.pair[, exFiltervcf.cmd := bsub.head(exFiltervcf.jobname, mem=5, We='1:11', cpu=1, postdone = filtervcf.jobname), by=1:nrow(targets.pair)]
targets.pair[, exFiltervcf.cmd := paste0(exFiltervcf.cmd, " \"", GATK4, "  --java-options '-Xms256m -Xmx5g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=", tmpdir, "' SelectVariants")]
targets.pair[, exFiltervcf.cmd := paste0(exFiltervcf.cmd, " --exclude-filtered true -O ", exFilterMutectVcfFile, " -V ", filterMutectVcfFile, " --sample-name ", bicSampleName, " \"")]
targets.pair[1, exFiltervcf.cmd]

exe.jobs(targets.pair$exFiltervcf.cmd, logger)

## annotate vcf per sample
targets.pair[, annoVcfFile := paste0(initDir, '/', tumor, '/neo_anno_', tumor, '.vcf')]
targets.pair[, annovcf.jobname := paste0("annovcf.", tumor)]
targets.pair[, annovcf.cmd := bsub.head(annovcf.jobname, mem=20, We='1:11', cpu=5, postdone = exFiltervcf.jobname), by=1:nrow(targets.pair)]
targets.pair[, annovcf.cmd := paste0(annovcf.cmd, " \"", PERL, '/perl ',  VEPv88, " -i ", exFilterMutectVcfFile, " -o ", annoVcfFile, " -offline --assembly GRCh37 --format vcf --vcf --symbol --plugin Downstream --plugin Wildtype --terms SO --dir_plugins ", VEP_plugins, " -fork 4 --force_overwrite\"") ]
targets.pair[1, annovcf.cmd]

exe.jobs(targets.pair$annovcf.cmd, logger)

## get protein fasta
#targets.pair[, mutPepFile := paste0(initDir, '/', tumor, '/mut_', tumor, '.fasta')]
#targets.pair[, mut.pep.jobname := paste0("mut.pep.", tumor)]
#targets.pair[, mutpep.cmd := bsub.head(mut.pep.jobname, mem=10, We='1:11', cpu=1, postdone = annovcf.jobname), by = 1:nrow(targets.pair)]
#targets.pair[, mutpep.cmd := paste0(mutpep.cmd, " \"", pvacseq, "  generate_protein_fasta ", mafAnnoVcfFile, " 21 ", mutPepFile, " \"")]
#targets.pair[1, mutpep.cmd]

# exe.jobs(targets.pair$mutpep.cmd, logger)

## run neoantigen prediction
#pvacseq run  DS-bla-185-T2 <- anno.vcf  Test  HLA-G*01:09,HLA-E*01:01,H2-IAb  NetMHC PickPocket NNalign output  -e 9,10
targets.pair[, hlaFile := paste0(hla.polysolverOutDir, '/winners.hla.txt')]
targets.pair[, allele := {
	hlaFile = targets.pair$hlaFile[3];fread(hlaFile, header=F) -> hla; hla
	hlaFile = targets.pair$hlaFile[3];
	scan(hlaFile, character(), sep="\n") -> hla; 
	hla = strsplit(hla, "\t")
	names(hla) = lapply(hla, "[[", 1)
	hla = lapply(hla, function(x){x[2:length(x)]})
	hla2 = lapply(hla, function(x){strsplit(x, "_") -> xx; unlist(lapply(xx, function(xx){paste0(toupper(xx[1]), "-", toupper(xx[2]), "*", paste(xx[3:4], collapse=":"), collapse="")}))})
	#hla2 = lapply(hla, function(x){strsplit(x, "_") -> xx; unlist(lapply(xx, function(xx){paste0(toupper(xx[1]), "-", toupper(xx[2]), "*", xx[3], ':', xx[4], collapse="")}))})
	hla2 = paste(unlist(hla2), collapse=",")
	hla2 }, by=1:nrow(targets.pair)]
targets.pair$allele

targets.pair[, sampleNeoDir := paste0(initDir, '/', tumor, '/neoantigeon')]
targets.pair[, {if(!dir.exists(sampleNeoDir)){dir.create(sampleNeoDir)}else{system(paste0("rm -rf ", sampleNeoDir, "/*"))}}, by = 1:nrow(targets.pair)]

targets.pair[, neo.jobname := paste0("neo.", tumor)]
targets.pair[, neo.cmd := bsub.head(neo.jobname, mem=30, We='20:11', cpu=5, postdone = annovcf.jobname), by=1:nrow(targets.pair)]
targets.pair[, neo.cmd := paste0(neo.cmd, " \"", pvacseq, " run --iedb-install-directory  ", IEDB, " ", annoVcfFile, " ", tumor, " ", allele, " NetMHCpan ", sampleNeoDir, " -e 9,10 \"")]

targets.pair[1, neo.cmd]
file.exists(targets.pair$annoVcfFile)
file.size(targets.pair$annoVcfFile)

exe.jobs(targets.pair$neo.cmd, logger)

##
targets.pair[, neoFile := paste0(sampleNeoDir, '/MHC_Class_I/', tumor, '.combined.parsed.tsv')]
file.exists(targets.pair$neoFile)
file.size(targets.pair$neoFile)
lapply(targets.pair$neoFile, fread) -> neo.dt
names(neo.dt) = targets.pair$tumor
neo.dt = Map(cbind, neo.dt, bcr = names(neo.dt))
rbindlist(neo.dt) -> neo.dt
head(neo.dt)
dim(neo.dt)

fwrite(neo.dt, sep="\t", quote=F, row.names=F, file=paste0(neoDir, '/', pre, '_neoantigen.tsv'))
rm(neo.dt)



