pre = 'scc'

source('~/program/configure_dnaseq.R')

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
		msg = "fastq file name problem"
		fatal(logger, msg)
		stop(msg)
	}
	rbind(targets.aln, tmp) -> targets.aln
}
targets.aln

## r_fang/sampleName_group_runID as file folder for each 
targets.aln[, sgrTag := paste0(sampleName, '_', group, '_', runID)]
targets.aln[, sgrDir := paste0(initDir, '/', sgrTag)]
targets.aln[, {system(paste0('mkdir -p ', sgrDir))}, by=1:nrow(targets.aln)]

## zcat fastq
targets.aln[, fastq1.zcatfile := file.path(sgrDir, sub(".gz", "", fastq1))]
targets.aln[, zcat1.jobname := paste0(resDir, '.', "zcat1.", sampleName, '.', fastq1)]
targets.aln[, zcat1.cmd := bsub.head(zcat1.jobname, mem=1, cpu=1, We='6:59', cwd), by=1:nrow(targets.aln)]
targets.aln[, zcat1.cmd := paste0(zcat1.cmd, ' "/bin/zcat ', fastq1.full, ' > ', fastq1.zcatfile, '"')]

targets.aln[, fastq2.zcatfile := file.path(sgrDir, sub(".gz", "", fastq2))]
targets.aln[, zcat2.jobname := paste0(resDir, '.', "zcat2.", sampleName, '.', fastq2)]
targets.aln[, zcat2.cmd := bsub.head(zcat2.jobname, mem=1, cpu=1, We='6:59', cwd), by=1:nrow(targets.aln)]
targets.aln[, zcat2.cmd := paste0(zcat2.cmd, ' "/bin/zcat ', fastq2.full, ' > ', fastq2.zcatfile, '"')]

exe.jobs(targets.aln$zcat1.cmd, logger)
exe.jobs(targets.aln$zcat2.cmd, logger)

## convert quality score
targets.aln[, fastq1.cqsfile := file.path(sgrDir, sub(".gz", ".cqs", fastq1))]
targets.aln[, cqs1.jobname := paste0(resDir, '.', "cqs1.", sampleName, '.', fastq1)]
targets.aln[, cqs1.cmd := bsub.head(cqs1.jobname, mem=1, cpu=1, We='2:59', cwd, zcat1.jobname), by=1:nrow(targets.aln)]
targets.aln[, cqs1.cmd := paste0(cqs1.cmd, ' "', ConvertQualityScore, ' --input ', fastq1.zcatfile, ' --output ', fastq1.cqs1file, '"')]
targets.aln$cqs1.cmd[1]

targets.aln[, fastq2.cqsfile := file.path(sgrDir, sub(".gz", ".cqs", fastq2))]
targets.aln[, cqs2.jobname := paste0(resDir, '.', "cqs.", sampleName, '.', fastq2)]
targets.aln[, cqs2.cmd := bsub.head(cqs2.jobname, mem=1, cpu=1, We='2:59', cwd, zcat2.jobname), by=1:nrow(targets.aln)]
targets.aln[, cqs2.cmd := paste0(cqs2.cmd, ' "', ConvertQualityScore, ' --input ', fastq2.zcatfile, ' --output ', fastq2.cqsfile, '"')]
targets.aln$cqs2.cmd[1]

exe.jobs(targets.aln$cqs1.cmd, logger)
exe.jobs(targets.aln$cqs2.cmd, logger)

## discard read less than half length
read.len = function(x){
	gzfile(x, open='r') -> fh
	readLines(fh, 2) -> h
	nchar(h[2])
}
targets.aln[, read.len := sapply(targets.aln$fastq1.full, read.len)]
targets.aln$read.len


## cutadapt need both fastq file jobs done
wait4jobs(c(targets.aln$cqs1.jobname, targets.aln$cqs2.jobname), logger)

file.size(targets.aln$fastq1.zcatfile)
del.files(targets.aln$fastq1.zcatfile, logger)
del.files(targets.aln$fastq2.zcatfile, logger)

# cutadapt
targets.aln[, fastq1.cutadaptfile :=  file.path(sgrDir, sub(".gz", ".cutadapt", fastq1))]
targets.aln[, fastq2.cutadaptfile :=  file.path(sgrDir, sub(".gz", ".cutadapt", fastq2))]
targets.aln[, cutadapt.jobname := paste0(resDir, '.', "cutadapt.", sgrTag)]
targets.aln[, cutadapt.cmd := bsub.head(cutadapt.jobname, mem=1, cpu=1, We='2:59', cwd, postdone=c(fastq1.cqsfile, fastq2.cqsfile)), by=1:nrow(targets.aln)]
targets.aln[, cutadapt.cmd := paste0(cutadapt.cmd, ' "', PYTHON, ' ', CUTADAPT, ' -f fastq -a ', clipR1, ' -A ', clipR2, ' --overlap 10 --minimum-length ', floor(read.len/2))]
targets.aln[, cutadapt.cmd := paste0(cutadapt.cmd, ' --quality-cutoff ', bqtrim, ' -o ', fastq1.cutadaptfile, ' --paired-output ', fastq2.cutadaptfile, ' ', fastq1.cqsfile, ' ', fastq2.cqsfile, ' "')]
targets.aln$cutadapt.cmd[1]

exe.jobs(targets.aln$cutadapt.cmd, logger)

# align
targets.aln[, reference := B37_BWA_INDEX]
targets.aln[, readgroup := paste0("@RG\\\\tID:", sgrTag, '_', seqType, '\\\\tPL:Illumina\\\\tPU:', sgrTag, "\\\\tLB:", sampleName, '_', group, "\\\\tSM:", sampleName)]
targets.aln[, bamfile := paste0(sgrDir, '/', sgrTag, '_', fastq1, '.bam')] # 
targets.aln$bamfile[1]
file.exists(targets.aln$bamfile)
file.exists(targets.aln$fastq2.cutadaptfile)

targets.aln[, bwa.jobname := paste0(resDir, '.', "bwa.", sgrTag, '_', fastq1)]
targets.aln[, bwa.cmd := bsub.head(bwa.jobname, mem=15, cpu=12, We='6:59', cwd, postdone=cutadapt.jobname), by=1:nrow(targets.aln)]
targets.aln[, bwa.cmd := paste0(bwa.cmd, ' "', bwa, ' mem -t 8 -T 0 -R \\\"', readgroup, '\\\" ', reference, ' ', fastq1.cutadaptfile, ' ', fastq2.cutadaptfile, ' | ', samtools, ' view -Shb -o ', bamfile, ' - "')]

exe.jobs(targets.aln$bwa.cmd, logger)

## merge bam file for the same library/group 
## may several r1/r2 fastq files
tmpdir='/scratch/huw/'
targets.bam = targets.aln[, .(bamFileList = paste(paste0(' I=', bamfile), collapse = ' '), bwaJobList = paste(bwa.jobname, collapse = ' ')), by = list(sampleName, group)]
targets.bam
targets.bam$bwaJobList
targets.bam[, libMergedBamFile := paste0(resDir, '/', sampleName, '_library', group, '_merged.bam')]
targets.bam[, merge.jobname := paste0(resDir, '.', 'mergebam.', sampleName, group)]
targets.bam$merge.jobname
targets.bam[, merge.cmd := bsub.head(merge.jobname, mem=30, cpu=8, We='2:26', cwd=cwd, postdone = bwaJobList), by=1:nrow(targets.bam)]
targets.bam[, merge.cmd := paste0(merge.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar MergeSamFiles ', bamFileList, ' O=', libMergedBamFile, ' SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=', tmpdir, ' CREATE_INDEX=TRUE USE_THREADING=FALSE MAX_RECORDS_IN_RAM=5000000 "')] 
targets.bam$merge.cmd[2]

exe.jobs(targets.bam$merge.cmd, logger)

## markDuplicates for reads from the same library/group
## output in resDir (r_fang)
tmpdir='/scratch/huw/'
targets.bam[, mkdupedLibBamFile := paste0(resDir, '/', sampleName, '_library', group, '_merged_mkdup.bam')]
targets.bam[, libmkdupMetricsFile := paste0(matricsDir, '/', sampleName, '_library', group, '_merged_mkdup_metrics.txt')]
targets.bam[, mkdup.jobname := paste0(resDir, '.', 'mkdup.', sampleName, '.', group)]
targets.bam$mkdup.jobname
targets.bam[, mkdup.cmd := bsub.head(mkdup.jobname, mem=30, cpu=3, We='8:26', cwd=cwd, postdone = merge.jobname), by=1:nrow(targets.bam)]
targets.bam[, mkdup.cmd := paste0(mkdup.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar MarkDuplicates I=', libMergedBamFile, ' O=', mkdupedLibBamFile, ' METRICS_FILE=', libmkdupMetricsFile, ' VALIDATION_STRINGENCY=LENIENT TMP_DIR=', tmpdir, ' REMOVE_DUPLICATES=TRUE CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=5000000 "')] 
targets.bam$mkdup.cmd[3]

exe.jobs(targets.bam$mkdup.cmd, logger)

## merge bam file for each sample
tmpdir='/scratch/huw/'
targets.smpBam[, sampleDir := paste0(initDir, '/', sampleName), by=1:nrow(targets.smpBam)]
targets.smpBam[, {system(paste0('mkdir -p ', sampleDir))}, by=1:nrow(targets.smpBam)]
targets.smpBam = targets.bam[, .(BamFileList = paste(paste0(' I=', mkdupedLibBamFile), collapse = ' '),  mkdupJobList = paste(mkdup.jobname, collapse = ' '), bamFile = paste(mkdupedLibBamFile, collapse=' '), fileNumber = length(mkdupedLibBamFile)), by = list(sampleName)]
colnames(targets.smpBam)
targets.smpBam
targets.smpBam[, smpBamFile := paste0(sampleDir, '/', sampleName, '_merged_mkdup_smp.bam')]
targets.smpBam[, smp.jobname := paste0(resDir, '.', 'smp.merge.', sampleName)]
targets.smpBam[, smp.cmd := bsub.head(smp.jobname, mem=30, cpu=8, We='8:26', cwd=cwd, postdone = mkdupJobList), by=1:nrow(targets.smpBam)]
targets.smpBam[, smp.cmd := paste0(smp.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar MergeSamFiles ', BamFileList, ' TMP_DIR=', tmpdir, ' O=', smpBamFile, ' SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE USE_THREADING=FALSE MAX_RECORDS_IN_RAM=5000000 "')] 
targets.smpBam$smp.cmd[11] 

exe.jobs(targets.smpBam$smp.cmd, logger) # 20mb

## breakpoint before collect metrics
wait4jobs(targets.smpBam$smp.jobname, logger)

# test generated files
file.exists(targets.smpBam$smpBamFile)
file.size(targets.smpBam$smpBamFile)

## clean files


targets.aln[, deldir(sgrDir), by=1:nrow(targets.aln)]

targets.bam[, libMergedBaiFile := sub("bam", "bai", libMergedBamFile)]
targets.bam[, mkdupedLibBaiFile := sub("bam", "bai", mkdupedLibBamFile)]
targets.bam[, delfile.2(libMergedBaiFile, logger), by=1:nrow(targets.bam)]
targets.bam[, delfile.2(libMergedBamFile, logger), by=1:nrow(targets.bam)]
targets.bam[, delfile.2(mkdupedLibBaiFile, logger), by=1:nrow(targets.bam)]
targets.bam[, delfile.2(mkdupedLibBamFile, logger), by=1:nrow(targets.bam)]

##  calculate HS matrix  for the bam file for each sample
## the baits targets and the bam file
## variant_pipeline: 1402
tmpdir='/scratch/huw/'
targets.smpBam[, hmsFile := paste0(sampleDir, '/', 'hms_', sampleName, '_HsMatrix.txt')]
targets.smpBam[, hms.jobname := paste0(resDir, '.', 'hms.', sampleName)]
targets.smpBam[, hms.cmd := bsub.head(hms.jobname, mem=10, cpu=1, We='8:26', cwd=cwd), by=1:nrow(targets.smpBam)] # dont need postdone since there is wait4jobs before
targets.smpBam[, hms.cmd := paste0(hms.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar CalculateHsMetrics I=', )]
targets.smpBam[, hms.cmd := paste0(hms.cmd, smpBamFile, ' TMP_DIR=', tmpdir, ' O=', hmsFile, ' REFERENCE_SEQUENCE=', genomeFasta, ' METRIC_ACCUMULATION_LEVEL=SAMPLE BAIT_INTERVALS=')]
targets.smpBam[, hms.cmd := paste0(hms.cmd, smpBamFile, baits_ilist, ' BAIT_SET_NAME=', wesImpactTarget, ' TARGET_INTERVALS=', targets_ilist, ' VALIDATION_STRINGENCY=LENIENT "')]
targets.smpBam$hms.cmd[3]

for(i in 1:nrow(targets.smpBam)){system(targets.smpBam$hms.cmd[i])}

## for pair end, calculate the insert size
tmpdir='/scratch/huw/'
targets.smpBam[, insertSizeFile := paste0(sampleDir, '/', 'readlen_', sampleName, '_insertSizeFile.txt')]
targets.smpBam[, insertSizeHistgramFile := paste0(sampleDir, '/', 'readlen_', sampleName, '_insertSizeFile.txt')]
targets.smpBam[, readlen.jobname := paste0(resDir, '.', 'readlen.', sampleName)]
targets.smpBam[, readlen.cmd := bsub.head(readlen.jobname, mem=10, cpu=1, We='8:26', cwd=cwd, postdone=hms.jobname), by=1:nrow(targets.smpBam)]
targets.smpBam[, readlen.cmd := paste0(readlen.cmd, ' "', JAVA, '/java -Xms256m -Xmx10g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar CollectInsertSizeMetrics I=', smpBamFile)]
targets.smpBam[, readlen.cmd := paste0(readlen.cmd, ' O=', insertSizeFile, ' REFERENCE_SEQUENCE=', genomeFasta, ' METRIC_ACCUMULATION_LEVEL=SAMPLE HISTOGRAM_FILE=', insertSizeHistgramFile)]
targets.smpBam[, readlen.cmd := paste0(readlen.cmd, ' VALIDATION_STRINGENCY=LENIENT TMP_DIR=', tmpdir, ' "')]
targets.smpBam$readlen.cmd[3]

for(i in 1:nrow(targets.smpBam)){system(targets.smpBam$readlen.cmd[i])}

##  CollectAlignmentSummaryMetrics
tmpdir='/scratch/huw/'
targets.smpBam[, alnSummaryFile := paste0(sampleDir, '/', 'alnSum_', sampleName, '_alnSummaryFile.txt')]
targets.smpBam[, alnSum.jobname := paste0(resDir, '.', 'alnSum.', sampleName)]
targets.smpBam[, alnSum.cmd := bsub.head(alnSum.jobname, mem=10, cpu=1, We='8:26', cwd=cwd, postdone=readlen.jobname), by=1:nrow(targets.smpBam)]
targets.smpBam[, alnSum.cmd := paste0(alnSum.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar  CollectAlignmentSummaryMetrics I=', smpBamFile)]
targets.smpBam[, alnSum.cmd := paste0(alnSum.cmd, ' O=', alnSummaryFile, ' REFERENCE_SEQUENCE=', genomeFasta, ' METRIC_ACCUMULATION_LEVEL=SAMPLE TMP_DIR=', tmpdir, ' "')]
targets.smpBam$alnSum.cmd[3]

for(i in 1:nrow(targets.smpBam)){system(targets.smpBam$alnSum.cmd[i])}

## CollectOxoGMetrics
tmpdir='/scratch/huw/'
targets.smpBam[, OxoFile := paste0(sampleDir, '/', 'Oxo_', sampleName, '_OxoFile.txt')]
targets.smpBam[, Oxo.jobname := paste0(resDir, '.', 'Oxo.', sampleName)]
targets.smpBam[, Oxo.cmd := bsub.head(Oxo.jobname, mem=4, cpu=1, We='8:26', cwd=cwd, postdone = alnSum.jobname), by=1:nrow(targets.smpBam)]
targets.smpBam[, Oxo.cmd := paste0(Oxo.cmd, ' "', JAVA, '/java -Xms256m -Xmx4g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar CollectOxoGMetrics  I=', smpBamFile)]
targets.smpBam[, Oxo.cmd := paste0(Oxo.cmd, ' REFERENCE_SEQUENCE=', genomeFasta, ' O=', OxoFile, ' DB_SNP=', DB_SNP, ' VALIDATION_STRINGENCY=LENIENT')];
targets.smpBam[, Oxo.cmd := paste0(Oxo.cmd,' TMP_DIR=', tmpdir, '  "')];
targets.smpBam$Oxo.cmd[3]

for(i in 1:nrow(targets.smpBam)){system(targets.smpBam$Oxo.cmd[i])}

## depth of coverage
tmpdir='/scratch/huw/'
targets.smpBam[, depthCoverageFile := paste0(sampleDir, '/', 'depthCoverage_', sampleName, '_depthCoverageFile.txt')]
targets.smpBam[, depthCoverage.jobname := paste0(resDir, '.', 'depthCoverage.', sampleName)]
targets.smpBam[, depthCoverage.cmd := bsub.head(depthCoverage.jobname, mem=4, cpu=1, We='8:26', cwd=cwd, postdone=Oxo.jobname), by=1:nrow(targets.smpBam)]
targets.smpBam[, depthCoverage.cmd := paste0(depthCoverage.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar -T DepthOfCoverage')]
targets.smpBam[, depthCoverage.cmd := paste0(depthCoverage.cmd, ' -I ', smpBamFile, ' -R ', genomeFasta, ' -o ', depthCoverageFile, ' -L ', FP_INT, ' -rf BadCigar')]
targets.smpBam[, depthCoverage.cmd := paste0(depthCoverage.cmd, ' -mmq 20 -mbq 0 -omitLocusTable -omitSampleSummary -baseCounts --includeRefNSites -omitIntervals "')]
targets.smpBam$depthCoverage.cmd[3]

for(i in 1:nrow(targets.smpBam)){system(targets.smpBam$depthCoverage.cmd[i])}

## GC bias metrics
tmpdir = "/scratch/huw/"
targets.smpBam[, gcbiasFile := paste0(sampleDir, '/', 'gcbias_', sampleName, '_gcbiasFile.txt')]
targets.smpBam[, gcbiasSummaryFile := paste0(sampleDir, '/', 'gcbias_', sampleName, '_gcbiasSummaryFile.txt')]
targets.smpBam[, gcbiasChartFile := paste0(sampleDir, '/', 'gcbias_', sampleName, '_gcbiasChartFile.pdf')]
targets.smpBam[, gcbias.jobname := paste0('gcbias.', sampleName)]
targets.smpBam[, gcbias.cmd := bsub.head(gcbias.jobname, mem=4, cpu=1, We='8:26', cwd=cwd, postdone = depthCoverage.jobname), by=1:nrow(targets.smpBam)]
targets.smpBam[, gcbias.cmd := bsub.head(gcbias.jobname, mem=4, cpu=1, We='8:26', cwd=cwd), by=1:nrow(targets.smpBam)]
targets.smpBam[, gcbias.cmd := paste0(gcbias.cmd, ' "', JAVA, '/java -Xms256m -Xmx4g -Djava.io.tmpdir=', tmpdir, ' -jar ', PICARD, '/picard.jar CollectGcBiasMetrics I=', smpBamFile)]
targets.smpBam[, gcbias.cmd := paste0(gcbias.cmd, ' O=', gcbiasFile, ' REFERENCE_SEQUENCE=', genomeFasta, ' SUMMARY_OUTPUT=', gcbiasSummaryFile)]
targets.smpBam[, gcbias.cmd := paste0(gcbias.cmd, ' CHART_OUTPUT=', gcbiasChartFile, ' TMP_DIR=', tmpdir, ' "')]
targets.smpBam$gcbias.cmd[3]

for(i in 1:nrow(targets.smpBam)){system(targets.smpBam$gcbias.cmd[i])}

## summarize metrics file
## need add

## cleaning and realignment bam files
## organize the targets
targets.smpGroup = fread(groupfile, header = F)
targets.smpGroup = setNames(targets.smpGroup, c('sampleName', 'smpGroup'))
targets.smpGroup
targets.smpGroup = merge(targets.smpBam[, .(sampleName, smpBamFile, smp.jobname)], targets.smpGroup, by = 'sampleName')
targets.smpGroup

## perform baserecal before realignment
## mutect2 and haplotypercaller will run without realignment sice realignment is not necessary for both of these

## baserecal by sample wise: calculate the parameters:
targets.smpGroup[, baserecalOutFile := paste0(sampleDir, '/', sampleName, '_recal_OutFile.txt')]
targets.smpGroup[, baserecal.jobname := paste0(resDir, '.', 'baserecal.', sampleName)]
targets.smpGroup[, baserecal.cmd := bsub.head(baserecal.jobname, mem=40, cpu=12, We='8:26', cwd=cwd), by=1:nrow(targets.smpGroup)]
targets.smpGroup[, baserecal.cmd := paste0(baserecal.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')]
targets.smpGroup[, baserecal.cmd := paste0(baserecal.cmd, ' -T BaseRecalibrator -R ', genomeFasta, ' --knownSites ', DB_SNP, ' --knownSites ', MILLS_1000G)] 
targets.smpGroup[, baserecal.cmd := paste0(baserecal.cmd, ' --knownSites ', HAPMAP, ' --knownSites ', OMNI_1000G, ' --knownSites ', PHASE1_SNPS_1000G)]
targets.smpGroup[, baserecal.cmd := paste0(baserecal.cmd, ' --covariate ContextCovariate --covariate CycleCovariate --covariate QualityScoreCovariate ')]
targets.smpGroup[, baserecal.cmd := paste0(baserecal.cmd, ' --covariate ReadGroupCovariate -rf BadCigar --num_cpu_threads_per_data_thread 12 ')]
targets.smpGroup[, baserecal.cmd := paste0(baserecal.cmd, ' -L ', targets5bp_bed)]
targets.smpGroup[, baserecal.cmd := paste0(baserecal.cmd, ' -S LENIENT -rf BadCigar --out ', baserecalOutFile, ' -I ', smpBamFile, ' "')]
targets.smpGroup$baserecal.cmd[1]

exe.jobs(targets.smpGroup$baserecal.cmd, logger)

## printreads
targets.smpGroup[, printReadsBamOutFile  := sub('.bam', "_recal.bam", smpBamFile)]
targets.smpGroup[, printReads.jobname := paste0(resDir, '.', 'printReads.', sampleName)]
targets.smpGroup[, printReads.cmd := bsub.head(printReads.jobname, mem=30, cpu=6, We='8:26', cwd=cwd, postdone = baserecal.jobname), by=1:nrow(targets.smpGroup)]
targets.smpGroup[, printReads.cmd := paste0(printReads.cmd, ' "', JAVA, '/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')]
targets.smpGroup[, printReads.cmd := paste0(printReads.cmd, ' -T PrintReads -R ', genomeFasta, ' --emit_original_quals -BQSR ', baserecalOutFile)]
targets.smpGroup[, printReads.cmd := paste0(printReads.cmd, ' --num_cpu_threads_per_data_thread 6 -rf BadCigar --downsampling_type NONE' )]
targets.smpGroup[, printReads.cmd := paste0(printReads.cmd, ' --out ', printReadsBamOutFile, ' -I ', smpBamFile, ' "')]
targets.smpGroup$printReads.cmd[1]

exe.jobs(targets.smpGroup$printReads.cmd, logger)


## REALIGNTARGETCREATOR \ref{realign} 
## these two steps are not necessary for mutect2 and haplotypecaller and can run later
## https://gatkforums.broadinstitute.org/gatk/discussion/2362/best-approach-for-realignertargetcreator-and-indelrealigner
targets.realn = targets.smpGroup[, .(bamFileList = paste(paste0('-I ', printReadsBamOutFile), collapse=' '), smpJobnameList = paste(printReads.jobname, collapse = ' ')), by=smpGroup]
targets.realn
tmpdir = "/scratch/huw/"

targets.realn[, realnGenTargetFile := paste0(matricsDir, '/', 'realn_', smpGroup, '_realnGenTargetFile.intervals')]
targets.realn[, realnGen.jobname := paste0(resDir, '.', 'realnGen_', smpGroup)]
targets.realn[, realnGen.cmd := bsub.head(realnGen.jobname, mem=15, cpu=1, We='8:26', cwd=cwd), by=1:nrow(targets.realn)]
targets.realn[, realnGen.cmd := paste0(realnGen.cmd, ' "', JAVA, '/java -Xms256m -Xmx15g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK)]
targets.realn[, realnGen.cmd := paste0(realnGen.cmd, '/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ', genomeFasta, ' -known ', DB_SNP, ' -known ', MILLS_1000G, ' ', bamFileList)]
targets.realn[, realnGen.cmd := paste0(realnGen.cmd, ' --out ', realnGenTargetFile, ' -L ', targets5bp_bed, ' "')]
targets.realn$realnGen.cmd[1]

exe.jobs(targets.realn$realnGen.cmd, logger)

## IndelRealigner
tmpdir = "/scratch/huw/"
targets.realn[, indelRealnOutFile := paste0(matricsDir, '/', 'realn_',smpGroup, '_indelRealn_OutFile.txt')]
targets.realn[, nwayoutMapFile := paste0(matricsDir, '/', 'realn_',smpGroup, '_mapfile.map')]
targets.realn[, {
	cat(NULL, file=nwayoutMapFile, append=F)
	bamfiles = gsub("-I ", "", bamFileList)
	bamfiles = basename(unlist(strsplit(bamfiles, split=" ")))
	bamoutfiles = sub(".bam", "_realnIndel.bam", bamfiles)
	for(i in 1:length(bamfiles)){
		cat(bamfiles[i], "\t", 'r_fang/', bamoutfiles[i], "\n", sep="", file=nwayoutMapFile, append=T)
	}
	}, by=1:nrow(targets.realn)]
	
targets.realn[, realnIndel.jobname := paste0(resDir, '.', 'realnIndel.', smpGroup)] 
targets.realn[, realnIndel.cmd := bsub.head(realnIndel.jobname, mem=15, cpu=1, We='8:26', cwd=cwd, postdone = realnGen.jobname), by=1:nrow(targets.realn)]
targets.realn[, realnIndel.cmd := paste0(realnIndel.cmd, ' "', JAVA, '/java -Xms256m -Xmx15g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')]
targets.realn[, realnIndel.cmd := paste0(realnIndel.cmd, ' -T IndelRealigner -R ', genomeFasta, ' -known ', DB_SNP, ' -known ', MILLS_1000G, ' -targetIntervals ', realnGenTargetFile)]
targets.realn[, realnIndel.cmd := paste0(realnIndel.cmd, ' -L ', targets5bp_bed)]
targets.realn[, realnIndel.cmd := paste0(realnIndel.cmd, ' --noOriginalAlignmentTags -S LENIENT --maxReadsForRealignment 500000 --maxReadsInMemory 3000000 --maxReadsForConsensuses 500000 -rf BadCigar ')]
targets.realn[, realnIndel.cmd := paste0(realnIndel.cmd, bamFileList, ' -nWayOut ', nwayoutMapFile, ' "')]  ## filename change here from .bam to _realnindel.bam
# bamFileList => ' -I file1.bam -I file2.bam'
# the output bamfiles will be ' -I file1_realnIndel.bam -I file2_realnIndel.bam'
targets.realn$realnIndel.cmd[1]

exe.jobs(targets.realn$realnIndel.cmd, logger)

## Realign use -known
## while recal use --knownSites

waitJobname = c(targets.smpGroup$printReads.jobname, targets.smpGroup$baserecal.jobname, targets.smpGroup$realnindel.jobname)
wait4jobs(waitJobname, logger)

## breakpoint before SNP calling
del.files(targets.smpGroup$smpBamFile, logger, targets.smpGroup$smpRealnBamFile, 20000000) # delete if smpBamFile exists and > 20mb
del.files(targets.smpGroup$smpRealnBamFile, logger, targets.smpGroup$printReadsBamOutFile, 20000000) # delete if smpBamFile exists and > 20mb

## construct pair samples
targets.pair = fread(pairfile, header=F)
targets.pair = setNames(targets.pair, c('normal', 'tumor'))
targets.pair[, normalSmpBamFile := paste0(resDir, '/', normal,'/', normal, '_merged_mkdup_smp.bam')]
targets.pair[, tumorSmpBamFile := paste0(resDir, '/', tumor, '/', tumor ,'_merged_mkdup_smp.bam')]
targets.pair[, normalSmpRecalBamFile := sub(".bam", "_recal.bam", normalSmpBamFile)]
targets.pair[, tumorSmpRecalBamFile := sub(".bam", "_recal.bam", tumorSmpBamFile)]
targets.pair[, normalSmpRecalRealnBamFile := sub(".bam", "_recal_realnIndel.bam", normalSmpBamFile)]
targets.pair[, tumorSmpRecalRealnBamFile  := sub(".bam", "_recal_realnIndel.bam", tumorSmpBamFile)]

targets.pair[, printReads.jobnames := paste0(resDir, '.', 'printReads.',normal, ' ', resDir, '.', 'printReads.', tumor)]

## mutect2 calling
## which dont need realign
## https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php
targets.pair[, mutectcallVCFOutFile  := paste0(mutectDir, '/', tumor, '.vcf')]
targets.pair[, mutectcall.jobname := paste0(resDir, '.', 'mutectcall.', tumor)]
targets.pair[, mutectcall.cmd := bsub.head(mutectcall.jobname, mem=90, cpu=30, We='8:26', cwd=cwd, postdone=printReads.jobnames), by=1:nrow(targets.pair)]
targets.pair[, mutectcall.cmd := paste0(mutectcall.cmd, ' "', JAVA, '/java -Xms256m -Xmx90g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')]
targets.pair[, mutectcall.cmd := paste0(mutectcall.cmd, ' -T MuTect2 -R ', genomeFasta, ' --dbsnp ', DB_SNP, ' --cosmic ', COSMIC, ' -I:normal ', normalSmpRecalBamFile)]
targets.pair[, mutectcall.cmd := paste0(mutectcall.cmd, ' -I:tumor ', tumorSmpRecalBamFile, ' -L ', targets5bp_bed, ' -o ', mutectcallVCFOutFile, ' -rf BadCigar --num_cpu_threads_per_data_thread 30 "')]
targets.pair$mutectcall.cmd[20]

exe.jobs(targets.pair$mutectcall.cmd, logger)

## wait for the realign
## somaticsniper
targets.pair[, snipercallVCFOutFile  := paste0(sniperDir, '/', tumor, '.vcf')]
targets.pair[, snipercall.jobname := paste0(resDir, '.', 'snipercall.', tumor)]
targets.pair[, snipercall.cmd := bsub.head(snipercall.jobname, mem=90, cpu=30, We='8:26', cwd=cwd, postdone =printReads.jobname), by=1:nrow(targets.pair)]
targets.pair[, snipercall.cmd := paste0(snipercall.cmd, ' "', SOMATICSNIPER, '/bam-somaticsniper -F vcf -f ', genomeFasta, ' -q 1 ', tumorBamFile, ' ', normalBamFile, ' ',  snipercallVCFOutFile, ' "')]
targets.pair$snipercall.cmd[1]

exe.jobs(targets.pair$snipercall.cmd, logger)

## facets: GetBaseCounts for each of the normal and tumor samples: normal samples were used more than once
MINCOV = 0;BASEQ = 20;MAPQ = 15
targets.pair[, facetsCountNOutFile  := paste0(facetsDir, '/sample/', normal, '_merged_mkdup_smp_recal_realnIndel.bam.dat')]
targets.pair[, facetsCountN.jobname := paste0(resDir, '.', 'facetsCountN.', normal)]
targets.pair[, facetsCountN.cmd := bsub.head(facetsCountN.jobname, mem=8, cpu=4, We='8:26', cwd=cwd, postdone=printReads.jobname), by=1:nrow(targets.pair)]
targets.pair[, facetsCountN.cmd := paste0(facetsCountN.cmd, ' "', facets, '/bin/GetBaseCounts  --thread 4 --filter_improper_pair --sort_output --fasta ')]
targets.pair[, facetsCountN.cmd := paste0(facetsCountN.cmd, genomeFasta, ' --vcf ', FACETS_DB_SNP, ' --maq ', MAPQ, ' --baq ', BASEQ, ' --cov ', MINCOV)]
targets.pair[, facetsCountN.cmd := paste0(facetsCountN.cmd, ' --bam ', normalBamFile, ' --out ', facetsCountNOutFile, ' "')] 
targets.pair$facetsCountN.cmd[1]

exe.jobs(targets.pair$facetsCountN.cmd[!duplicated(targets.pair$normalBamFile)], logger)

targets.pair[, facetsCountTOutFile  := paste0(facetsDir, '/sample/', tumor, '_merged_mkdup_smp_recal_realnIndel.bam.dat')]
targets.pair[, facetsCountT.jobname := paste0(resDir, '.', 'facetsCountT.', tumor)]
targets.pair[, facetsCountT.cmd := bsub.head(facetsCountT.jobname, mem=8, cpu=4, We='8:26', cwd=cwd, postdone=printReads.jobname), by=1:nrow(targets.pair)]
targets.pair[, facetsCountT.cmd := paste0(facetsCountT.cmd, ' "', facets, '/bin/GetBaseCounts  --thread 4 --filter_improper_pair --sort_output --fasta ')]
targets.pair[, facetsCountT.cmd := paste0(facetsCountT.cmd, genomeFasta, ' --vcf ', FACETS_DB_SNP, ' --maq ', MAPQ, ' --baq ', BASEQ, ' --cov ', MINCOV)]
targets.pair[, facetsCountT.cmd := paste0(facetsCountT.cmd, ' --bam ', tumorBamFile, ' --out ', facetsCountTOutFile, ' "')] 
targets.pair$facetsCountT.cmd[1]

exe.jobs(targets.pair$facetsCountT.cmd[!duplicated(targets.pair$tumorBamFile)], logger)

## facets mergeTN.R
targets.pair[, facetsMergeOutFile  := paste0(facetsDir, '/sample/count_merged_', tumor, '_', normal, '.dat.gz')]
targets.pair[, facetsMerge.jobname := paste0(resDir, '.', 'facetsMerge.', tumor)]
targets.pair[, facetsMerge.cmd := bsub.head(facetsMerge.jobname, mem=18, cpu=4, We='8:26', cwd=cwd, c(facetsCountT.jobname, facetsCountN.jobname)), by=1:nrow(targets.pair)]
targets.pair[, facetsMerge.cmd := paste0(facetsMerge.cmd, ' "', facets, '/mergeTN.R ', facetsCountTOutFile, ' ', facetsCountNOutFile, ' ', facetsMergeOutFile, ' "')]
targets.pair$facetsMerge.cmd[1]

exe.jobs(targets.pair$facetsMerge.cmd[!duplicated(targets.pair$tumorBamFile)], logger)

## facets: run
## facets.suite facets.lib dir tag file ggenome pc c
pcval = 200
cval = 100
targets.pair[, facetsRunOutDir  := paste0(facetsDir, '/', tumor)] ; 
for(i in 1:nrow(targets.pair)){system(paste0('mkdir -p ', targets.pair$facetsRunOutDir[i]))}
targets.pair[, facetsRun.jobname := paste0(resDir, '.', 'facetsRun.', tumor)]
targets.pair[, facetsRun.cmd := bsub.head(facetsRun.jobname, mem=3, cpu=3, We='8:26', cwd=cwd, postdone = targets.pair$facetsRun.jobname), by=1:nrow(targets.pair)]
targets.pair[, facetsRun.cmd := paste0(facetsRun.cmd, ' "', FACETS_SUITE, '/doFacets.R', ' --cval ', cval, ' --purity_cval ', pcval, ' --genome hg19')]
targets.pair[, facetsRun.cmd := paste0(facetsRun.cmd, ' --counts_file ', facetsMergeOutFile, ' --TAG ', tumor, ' --directory ', facetsRunOutDir, ' --R_lib ', RLIB_PATH, ' "')]
targets.pair$facetsRun.cmd[2]

exe.jobs(targets.pair$facetsRun.cmd, logger)

## facets: gene level
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
# gl.list = lapply(targets.pair$facetsGeneLevelFile, fread, sep="\t")
# rbindlist(gl.list) -> gl.all
# fwrite(gl.all, file=paste0(facetsDir, '/', pre, '_facets_gene_level'), sep = '\t')
# rm(gl.list, gl.all)

## muTect: DMP_rescue.py 
#targets.pair[, DMPrescueOutFile  := paste0(mutectDir, '/', tumor, '_rescue.vcf')] ; 
#targets.pair[, DMPrescue.jobname := paste0(resDir, '.', 'DMPrescue.', tumor)]
#targets.pair[, DMPrescue.cmd := bsub.head(DMPrescue.jobname, mem=3, cpu=1, We='8:26', cwd=cwd), by=1:nrow(targets.pair)]
#targets.pair[, DMPrescue.cmd := paste0(DMPrescue.cmd, ' "', exeDir, '/rescue/DMP_rescue.py -v ', mutectcallVCFOutFile, ' --txt ', mutectcallTXTOutFile, ' -t ', tumor, ' -n ', normal, ' -o ', DMPrescueOutFile, ' "')]
#targets.pair$DMPrescue.cmd[2]
#
#exe.jobs(targets.pair$DMPrescue.cmd, logger)

#glbs.TNRatio=float(5)
#glbs.TotalDepth=float(0)
#glbs.AlleleDepth=float(3)
#glbs.VariantFreq=float(0.01)
#tum_coverage  = int(tum_alt_ad) + int(tum_ref_ad)
#tum_alt_freq = float(tum_alt_ad)/ tum_coverage
#norm_alt_freq = float(norm_alt_ad)/ (int(norm_alt_ad) + int(norm_ref_ad))
#
## Now see if this passes all the thresholds!
#if tum_coverage >= glbs.TotalDepth and tum_alt_ad >= glbs.AlleleDepth and tum_alt_freq >= glbs.VariantFreq and tum_alt_freq >= glbs.TNRatio * norm_alt_freq:


## vcf2maf
## muTect: VCF2MAF for each sample 
## DMP_rescue was not used here
ncbi = 'GRCh37'
targets.pair[, mutect2mafOutFile  := paste0(mutectDir, '/', tumor, '.maf')] ; 
targets.pair[, mutect2maf.jobname := paste0(resDir, '.', 'mutect2maf.', tumor)]
targets.pair[, mutect2maf.cmd := bsub.head(mutect2maf.jobname, mem=3, cpu=3, We='8:26', cwd=cwd, postdone = targets.pair$mutectcall.jobname), by=1:nrow(targets.pair)]
targets.pair[, mutect2maf.cmd := paste0(mutect2maf.cmd, ' "', PERL, '/perl ', VCF2MAF, '/vcf2maf.pl --input ', mutectcallVCFOutFile, ' --output-maf ',  mutect2mafOutFile)]
targets.pair[, mutect2maf.cmd := paste0(mutect2maf.cmd, ' --ref-fasta ', genomeFasta, ' --tmp-dir ', tmpdir, '  --ncbi ', ncbi, ' --vep-path ', VEP, ' --vep-data ', VEP)]
targets.pair[, mutect2maf.cmd := paste0(mutect2maf.cmd, ' --tumor-id ', tumor, ' --normal-id ', normal, ' --filter-mutect ', ExAC_VCF, ' "')] # --updown-length 10 ??
targets.pair$mutect2maf.cmd[2]

exe.jobs(targets.pair$mutect2maf.cmd, logger) # 

## mafAnno: add CNV to each mutect maf files
targets.pair[, facetsHisensRdata := paste0(facetsDir, '/', tumor, '/', tumor, '_hisens.Rdata')] ; 
targets.pair[, mutectmafAnnoOutFile  := paste0(mutectDir, '/', tumor, '_CNV.maf')] ; 
targets.pair[, mutectmafAnno.jobname := paste0(resDir, '.', 'mutectmafAnno.', tumor)]
targets.pair[, mutectmafAnno.cmd := bsub.head(mutectmafAnno.jobname, mem=3, cpu=3, We='8:26', cwd=cwd, postdone=paste(mutect2maf.jobname, facetsRun.jobname, sep=' ')), by=1:nrow(targets.pair)]
targets.pair[, mutectmafAnno.cmd := paste0(mutectmafAnno.cmd, ' "', FACETS_SUITE, '/mafAnnov2.R -m ', mutect2mafOutFile, ' -f ', facetsHisensRdata, ' -o ', mutectmafAnnoOutFile, ' -n ', tumor, ' "')]
targets.pair$mutectmafAnno.cmd[2]

exe.jobs(targets.pair$mutectmafAnno.cmd, logger) # 

## combine mutect results
maf.list = lapply(targets.pair$mutect2mafOutFile, fread, sep="\t")
rbindlist(maf.list) -> maf.all
fwrite(maf.all, file=paste0(varDir, '/', pre, '.maf'), sep = '\t')

# sync before is haplotypercaller
wait4jobs(targets.pair$mutectmafAnno.jobname, logger)

### call SNP by haplotypecaller
#UnifiedGenotyper unifiedGenotyper.vcf
#HaplotypeCaller haplotypecaller.vcf
#CombineVariants of above tow => HaplotypeCaller_RAW.vcf
#VariantRecalibrator input:HaplotypeCaller_RAW.vcf 
#ApplyRecalibration input:HaplotypeCaller_RAW.vcf output:HaplotypeCaller_SNP_vqsr.vcf
#CombineVariants: output:_UnifiedGenotyper_RAW.vcf  input:ugVars
#VariantRecalibrator: 
#ApplyRecalibration:
#VariantRecalibrator
#ApplyRecalibration

## split bam file by chromosomome
## does not seem necessary
#smp = unique(c(targets.pair$normal,targets.pair$tumor))
#fai = fread(genomeFAI)
#fai = setNames(fai, c('chrom', 'len', 'acc.len', 'V4', 'V5' ))
#setkey(fai, chrom)
#smpBamFile.rep = rep(smp, each=nrow(fai))
#targets.chrBam = data.table(sampleName = smpBamFile.rep, chrom = rep(fai$chrom, times = length(smpBamFile)))
#targets.chrBam
#length(smpBamFile)
#nrow(fai)
#nrow(targets.chrBam)
#targets.chrBam
#targets.chrBam[, smpBamFile := paste0(resDir, '/', sampleName, '_merged_mkdup_smp_realnIndel_recal.bam')]
#targets.chrBam[, chrBamFile := paste0(resDir, '/', sampleName, '/', sampleName, '_', chrom, '_merged_mkdup_smp_realnIndel_recal.bam')]
#targets.chrBam[, sampleDir:= paste0(resDir, '/', sampleName)]
#targets.chrBam[, {system(paste('mkdir -p', sampleDir))}, by=1:nrow(targets.chrBam)]
#
#targets.chrBam[, splitBam.jobname := paste0(resDir, '.', 'splitBam.', sampleName)]
#targets.chrBam[, splitBam.cmd := bsub.head(splitBam.jobname, mem=9, cpu=3, We='8:26', cwd=cwd), by=1:nrow(targets.chrBam)]
#targets.chrBam[, splitBam.cmd := paste0(splitBam.cmd, ' "', samtools, ' view -b -o ', chrBamFile, ' ', smpBamFile, ' ', chrom, ' "')]
#targets.chrBam$splitBam.cmd[1]
#
#exe.jobs(targets.chrBam$splitBam.cmd, logger)


## haplotypecaller
targets.pair[, haplocallVCFOutFile  := paste0(haploDir, '/', tumor, '.vcf')]
targets.pair[, haplocall.jobname := paste0(resDir, '.', 'haplocall.', normal)]
targets.pair[, haplocall.cmd := bsub.head(haplocall.jobname, mem=90, cpu=30, We='8:26', cwd=cwd), by=1:nrow(targets.pair)]
targets.pair[, haplocall.cmd := paste0(haplocall.cmd, ' "', JAVA, '/java -Xms256m -Xmx90g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')]
targets.pair[, haplocall.cmd := paste0(haplocall.cmd, ' -T HaplotypeCaller -R ', genomeFasta, ' -L ', targets5bp_bed, ' --dbsnp ', DB_SNP, ' --downsampling_type NONE --annotation AlleleBalanceBySample')]
targets.pair[, haplocall.cmd := paste0(haplocall.cmd, ' --annotation ClippingRankSumTest --read_filter BadCigar --num_cpu_threads_per_data_thread 30 --out ', haplocallVCFOutFile)]
targets.pair[, haplocall.cmd := paste0(haplocall.cmd, ' --max-alternate-alleles 2 -minPruning 4 -I ', normalSmpBamFile, ' "')]
targets.pair$haplocall.cmd[1]
targets.pair$haplocall.cmd[!duplicated(targets.pair$normal)]

exe.jobs(targets.pair$haplocall.cmd[!duplicated(targets.pair$normal)], logger) 

## combine haplotypercaller results
inputHaplocallVCFFile = paste(targets.pair$tumor, targets.pair$haplocallVCFOutFile, sep=" ")
inputHaplocallVCFFile = paste(inputHaplocallVCFFile, collapse=" --variant:")
inputHaplocallVCFFile = paste0("--variant:", inputHaplocallVCFFile)
inputHaplocallVCFFile

haplocallCombinedVCFOutFile  = paste0(haploDir, '/', 'combined_haplocall.vcf')
combineHaploVcf.jobname = paste0(resDir, '.', 'combineHaploVcf')
combineHaploVcf.cmd = bsub.head(combineHaploVcf.jobname, mem=3, cpu=1, We='826', cwd=cwd, postdone=paste(targets.pair$haplocall.jobname, collapse=" "))
combineHaploVcf.cmd = paste0(combineHaploVcf.cmd, ' "', JAVA, '/java -Xms256m -Xmx3g -XX-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')
combineHaploVcf.cmd = paste0(combineHaploVcf.cmd, ' -T CombineVariants -R ', genomeFasta, ' --assumeIdenticalSamples ') 
combineHaploVcf.cmd = paste0(combineHaploVcf.cmd, ' -o ', haplocallCombinedVCFOutFile)
combineHaploVcf.cmd = paste0(combineHaploVcf.cmd, ' -I ', inputHaplocallVCFFile, ' "')
combineHaploVcf.cmd

exe.jobs( combineHaploVcf.cmd, logger)

## haplotypecaller: VariantRecalibrator: SNP
varRecalSNPOutFile  = paste0(haploDir, '/', 'combined_haplocall_varRecal_SNP.recal')
varRecalSNPTranchesOutFile  = paste0(haploDir, '/', 'combined_haplocall_varRecal_tranches_SNP.tranches')
varRecalSNPPlotOutFile  = paste0(haploDir, '/', 'combined_haplocall_varRecal_plots_SNP.R')
varRecalSNP.jobname = paste0(resDir, '.', 'varRecalSNP')
varRecalSNP.cmd = bsub.head(varRecalSNP.jobname, mem=10, cpu=4, We='8:26', cwd=cwd, postdone = combineHaploVcf.jobname)
varRecalSNP.cmd = paste0(varRecalSNP.cmd, ' "', JAVA, '/java -Xms256m -Xmx10g -XX-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')
varRecalSNP.cmd = paste0(varRecalSNP.cmd, ' -T VariantRecalibrator -R ', genomeFasta, ' -input ', haplocallCombinedVCFOutFile, ' -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ', HAPMAP)
varRecalSNP.cmd = paste0(varRecalSNP.cmd, ' -resource:omni,known=false,training=true,truth=true,prior=12.0 ', OMNI_1000G) 
varRecalSNP.cmd = paste0(varRecalSNP.cmd, ' -resource:1000G,known=false,training=true,truth=false,prior=10.0 ', PHASE1_SNPS_1000G)
varRecalSNP.cmd = paste0(varRecalSNP.cmd, ' -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ', DB_SNP)
varRecalSNP.cmd = paste0(varRecalSNP.cmd, ' -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 ')
varRecalSNP.cmd = paste0(varRecalSNP.cmd, ' -recalFile ', varRecalSNPOutFile, ' -tranchesFile ', varRecalSNPTranchesOutFile, ' -rscriptFile ', varRecalSNPPlotOutFile, ' -nt 4 "')
varRecalSNP.cmd

exe.jobs(varRecalSNP.cmd, logger)

## haplotypecaller: apply VariantRecalibrator: SNP
varRecaledSNPFile  = paste0(haploDir, '/', 'combined_haplocall_varRecaled_SNP.vcf')
varRecalSNPApply.jobname = paste0(resDir, '.', 'varRecalSNP.Apply')
varRecalSNPApply.jobname
varRecalSNPApply.cmd = bsub.head(varRecalSNPApply.jobname, mem=3, cpu=1, We='8:26', cwd=cwd, postdone = varRecalSNP.jobname)
varRecalSNPApply.cmd = paste0(varRecalSNPApply.cmd, ' "', JAVA, '/java -Xms256m -Xmx3g -XX-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')
varRecalSNPApply.cmd = paste0(varRecalSNPApply.cmd, ' -T  ApppyRecalibration -R ', genomeFasta, ' -input ', haplocallCombinedVCFOutFile, ' --ts_filter_level 99.0 -mode SNP ')
varRecalSNPApply.cmd = paste0(varRecalSNPApply.cmd, ' -tranchesFile ', varRecalSNPTranchesOutFile)
varRecalSNPApply.cmd = paste0(varRecalSNPApply.cmd, ' -recalFile ', varRecalSNPOutFile)
varRecalSNPApply.cmd = paste0(varRecalSNPApply.cmd, ' -o ', varRecalSNPedFile)
varRecalSNPApply.cmd

exe.jobs(varRecalSNPApply.cmd, logger)

## haplotypecaller: VariantRecalibrator: INDEL
varRecalINDELOutFile  = paste0(haploDir, '/', 'combined_haplocall_varRecal_INDEL.recal')
varRecalINDELTranchesOutFile  = paste0(haploDir, '/', 'combined_haplocall_varRecal_tranches_INDEL.tranches')
varRecalINDELPlotOutFile  = paste0(haploDir, '/', 'combined_haplocall_varRecal_plots_INDEL.R')
varRecalINDEL.jobname = paste0(resDir, '.', 'varRecalINDEL')
varRecalINDEL.cmd = bsub.head(varRecalINDEL.jobname, mem=10, cpu=4, We='8:26', cwd=cwd, postdone = combineHaploVcf.jobname)
varRecalINDEL.cmd = paste0(varRecalINDEL.cmd, ' "', JAVA, '/java -Xms256m -Xmx10g -XX-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')
varRecalINDEL.cmd = paste0(varRecalINDEL.cmd, ' -T VariantRecalibrator -R ', genomeFasta, ' -input ', varRecalSNPOutFile, ' -resource:mills,known=false,training=true,truth=true,prior=12.0 ', MILL_1000G)
varRecalINDEL.cmd = paste0(varRecalINDEL.cmd, ' -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 ')
varRecalINDEL.cmd = paste0(varRecalINDEL.cmd, ' -recalFile ', varRecalINDELOutFile, ' -tranchesFile ', varRecalINDELTranchesOutFile, ' -rscriptFile ', varRecalINDELPlotOutFile, ' -nt 4 "')
varRecalINDEL.cmd

exe.jobs(varRecalINDEL.cmd, logger)

## haplotypecaller: apply VariantRecalibrator: INDEL
varRecaledINDELFile  = paste0(haploDir, '/', 'combined_haplocall_varRecaled_INDEL.vcf')
varRecalINDELApply.jobname = paste0(resDir, '.', 'varRecalINDEL.Apply')
varRecalINDELApply.cmd = bsub.head(varRecalINDELApply.jobname, mem=3, cpu=1, We='8:26', cwd=cwd, postdone = varRecalINDEL.jobname)
varRecalINDELApply.cmd = paste0(varRecalINDELApply.cmd, ' "', JAVA, '/java -Xms256m -Xmx3g -XX-UseGCOverheadLimit -Djava.io.tmpdir=', tmpdir, ' -jar ', GATK, '/GenomeAnalysisTK.jar')
varRecalINDELApply.cmd = paste0(varRecalINDELApply.cmd, ' -T  ApppyRecalibration -R ', genomeFasta, ' -input ', varRecalSNPOutFile, ' --ts_filter_level 99.0 -mode INDEL ')
varRecalINDELApply.cmd = paste0(varRecalINDELApply.cmd, ' -tranchesFile ', varRecalINDELTranchesOutFile)
varRecalINDELApply.cmd = paste0(varRecalINDELApply.cmd, ' -recalFile ', varRecalINDELOutFile)
varRecalINDELApply.cmd = paste0(varRecalINDELApply.cmd, ' -o ', varRecalINDELedFile)
varRecalINDELApply.cmd

exe.jobs(varRecalINDELApply.cmd, logger)
## haplotypecaller: VCF2MAF for each sample 
## DMP_rescue was not used here
#ncbi = 'GRCh37'
#targets.pair[, haplocall2mafOutFile  := paste0(haploDir, '/', tumor, '.maf')] ; 
#targets.pair[, haplocall2maf.jobname := paste0(resDir, '.', 'haplocall2maf.', tumor)]
#targets.pair[, haplocall2maf.cmd := bsub.head(haplocall2maf.jobname, mem=3, cpu=3, We='8:26', cwd=cwd, postdone=haplocall.jobname), by=1:nrow(targets.pair)] ## need add postdone
#targets.pair[, haplocall2maf.cmd := paste0(haplocall2maf.cmd, ' "', PERL, '/perl ', VCF2MAF, '/vcf2maf.pl --input ', haplocallVCFOutFile, ' --output-maf ',  haplocall2mafOutFile)]
#targets.pair[, haplocall2maf.cmd := paste0(haplocall2maf.cmd, ' --ref-fasta ', genomeFasta, ' --tmp-dir ', tmpdir, '  --ncbi ', ncbi, ' --vep-path ', VEP, ' --vep-data ', VEP)]
#targets.pair[, haplocall2maf.cmd := paste0(haplocall2maf.cmd, ' --tumor-id ', tumor, ' --normal-id ', normal, ' --filter-mutect ', ExAC_VCF, ' "')] # --updown-length 10 ??
#targets.pair$haplocall2maf.cmd[2]
#
#exe.jobs(targets.pair$haplocall2maf.cmd, logger) # 

wait4jobs(c(targets.pair$mutect2maf.jobname, targets.pair$haplocall2maf), logger)

## mafAnno: add CNV to each haplotypecaller maf files
#targets.pair[, haplomafAnnoOutFile  := paste0(haploDir, '/', tumor, '_CNV.maf')] ; 
#targets.pair[, haplomafAnno.jobname := paste0(resDir, '.', 'haplomafAnno.', tumor)]
#targets.pair[, haplomafAnno.cmd := bsub.head(haplomafAnno.jobname, mem=3, cpu=3, We='8:26', cwd=cwd, postdone=paste(haplocall2maf.jobanme, facetsRun.jobname, sep=" ")), by=1:nrow(targets.pair)]
#targets.pair[, haplomafAnno.cmd := paste0(haplomafAnno.cmd, ' "', FACETS_SUITE, '/mafAnnov2.R -m ', haplocall2mafOutFile, ' -f ', facetsHisensRdata, ' -o ', haplomafAnnoOutFile, ' -n ', tumor, ' "')]
#targets.pair$haplomafAnno.cmd[2]
#
#exe.jobs(targets.pair$haplomafAnno.cmd, logger) # 

## combine all the happlocalltype and mutect maf files
#maf.list = lapply(c(targets.pair$haplocall2mafOutFile, targets.pair$mutect2mafOutFile), fread, sep="\t")
#rbindlist(maf.list) -> maf.all
#fwrite(maf.all, file=paste0(varDir, '/', pre, '.maf'), sep = '\t')

# filter
ngsfilter
# exome cnv
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DMP_CNV", job_hold => "$ssfj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_DMP_CNV.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/exome_cnv.pl -pre $pre -result $output/variants/copyNumber/dmp_cnv -berger $target_design -bamlist $output/intFiles/$pre\_sv_bam_list.txt -patient $patient -std_covg $target_std_normals -genome $species -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group`;
# structural variation
# generate bam file list: sampleName bamfileName
targets.pair[, group := substr(tumor, 1, 10)]
bamlistFile = paste0(strvarDir, '/bamlistFile.txt')
targets.pair[, delly2DELOutFile := paste0(strvarDir, '/', tumor, '_belly2_DEL.bcf')]
targets.pair[, delly2DEL.jobname := paste0(resDir, '.', 'delly2DEL.', tumor)]
targets.pair[, delly2DEL.cmd := bsub.head(delly2DEL.jobname, mem=10, cpu=3, We='8:26', cwd=cwd), by=1:nrow(targets.pair)]
targets.pair[, delly2DEL.cmd := paste0(delly2DEL.cmd, ' "', belly2, ' call -t DEL -g ', genomeFasta, ' -o ', delly2DELOutFile, ' ', tumorBamFile, ' ', normalBamFile, '"')]
targets.pair$delly2DEL.cmd[2]
exe.jobs(targets.pair$delly2DEL.cmd, logger) # 

targets.pair[, group := substr(tumor, 1, 10)]
bamlistFile = paste0(strvarDir, '/bamlistFile.txt')
targets.pair[, delly2DUPOutFile := paste0(strvarDir, '/', tumor, '_belly2_DUP.bcf')]
targets.pair[, delly2DUP.jobname := paste0(resDir, '.', 'delly2DUP.', tumor)]
targets.pair[, delly2DUP.cmd := bsub.head(delly2DUP.jobname, mem=10, cpu=3, We='8:26', cwd=cwd, delly2DEL.jobname), by=1:nrow(targets.pair)]
targets.pair[, delly2DUP.cmd := paste0(delly2DUP.cmd, ' "', belly2, ' call -t DUP -g ', genomeFasta, ' -o ', delly2DUPOutFile, ' ', tumorBamFile, ' ', normalBamFile, '"')]
targets.pair$delly2DUP.cmd[2]
exe.jobs(targets.pair$delly2DUP.cmd, logger) # 

targets.pair[, delly2INVOutFile := paste0(strvarDir, '/', tumor, '_belly2_INV.bcf')]
targets.pair[, delly2INV.jobname := paste0(resDir, '.', 'delly2INV.', tumor)]
targets.pair[, delly2INV.cmd := bsub.head(delly2INV.jobname, mem=10, cpu=3, We='8:26', cwd=cwd, delly2DUP.jobname), by=1:nrow(targets.pair)]
targets.pair[, delly2INV.cmd := paste0(delly2INV.cmd, ' "', belly2, ' call -t INV -g ', genomeFasta, ' -o ', delly2INVOutFile, ' ', tumorBamFile, ' ', normalBamFile, '"')]
targets.pair$delly2INV.cmd[2]
exe.jobs(targets.pair$delly2INV.cmd, logger) # 

targets.pair[, delly2BNDOutFile := paste0(strvarDir, '/', tumor, '_belly2_BND.bcf')]
targets.pair[, delly2BND.jobname := paste0(resDir, '.', 'delly2BND.', tumor)]
targets.pair[, delly2BND.cmd := bsub.head(delly2BND.jobname, mem=10, cpu=3, We='8:26', cwd=cwd, delly2INV.jobname ), by=1:nrow(targets.pair)]
targets.pair[, delly2BND.cmd := paste0(delly2BND.cmd, ' "', belly2, ' call -t BND -g ', genomeFasta, ' -o ', delly2BNDOutFile, ' ', tumorBamFile, ' ', normalBamFile, '"')]
targets.pair$delly2BND.cmd[2]
exe.jobs(targets.pair$delly2BND.cmd, logger) # 

targets.pair[, delly2INSOutFile := paste0(strvarDir, '/', tumor, '_belly2_INS.bcf')]
targets.pair[, delly2INS.jobname := paste0(resDir, '.', 'delly2INS.', tumor)]
targets.pair[, delly2INS.cmd := bsub.head(delly2INS.jobname, mem=10, cpu=3, We='8:26', cwd=cwd, delly2BND.jobname), by=1:nrow(targets.pair)]
targets.pair[, delly2INS.cmd := paste0(delly2INS.cmd, ' "', belly2, ' call -t INS -g ', genomeFasta, ' -o ', delly2INSOutFile, ' ', tumorBamFile, ' ', normalBamFile, '"')]
targets.pair$delly2INS.cmd[2]
exe.jobs(targets.pair$delly2INS.cmd, logger) # 

## combine delly
targets.pair[, delly2INSOutFile := paste0(strvarDir, '/', tumor, '_belly2_INS.bcf')]
targets.pair[, delly2INS.jobname := paste0(resDir, '.', 'delly2INS.', tumor)]
targets.pair[, delly2INS.cmd := bsub.head(delly2INS.jobname, mem=10, cpu=3, We='8:26', cwd=cwd, delly2BND.jobname), by=1:nrow(targets.pair)]
targets.pair[, delly2INS.cmd := paste0(delly2INS.cmd, ' "', belly2, ' call -t INS -g ', genomeFasta, ' -o ', delly2INSOutFile, ' ', tumorBamFile, ' ', normalBamFile, '"')]
targets.pair$delly2INS.cmd[2]
exe.jobs(targets.pair$delly2INS.cmd, logger) # 
# cDNA contamination
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CDNA_CONTAM", job_hold => "$hold_value", cpu => "1", mem => "4", cluster_out => "$output/progress/$pre\_$uID\_CDNA_CONTAM.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/check_cDNA_contamination.py -s $output/variants/structVar/delly/$pre\_AllAnnotatedSVs.txt -o $output/metrics/$pre\_cDNA_contamination.txt`;
# hybrid genome





