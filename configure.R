library(log4r)
library(data.table)

cwd = getwd()

source(paste0('~/pipeline/configure_dir.R'))

if(assay == 'methylation'){
	source('~/program/fun/param.r')
	initDir = paste0(resDir, '/initFiles')
	progressDir = paste0(resDir, '/progress')
	matricsDir = paste0(resDir, '/matrics')
	fastqcDir = paste0(resDir, '/fastqc')
	statsDir = paste0(resDir, '/stats')
	methylDir = paste0(resDir, '/methylDir')
	system(paste0('mkdir -p ', resDir)) ## r_fang
	system(paste0('mkdir -p ', progressDir)) ## jobname.done
	system(paste0('mkdir -p ', matricsDir))  ## collected metrics
	system(paste0('mkdir -p ', initDir)) 
	system(paste0('mkdir -p ', fastqcDir)) 
	system(paste0('mkdir -p ', statsDir)) 
}

logfile = paste0(cwd, '/', '_', gsub("/", "_", resDir), '_log')
logfile
create.logger() -> logger
logfile(logger) = logfile
level(logger) = 'INFO' 

## this for set up @ENV
path			= '/home/huw/program/homer/bin:/home/huw/local/bin:/home/huw/perl5/CPAN/bin:/home/huw/program/cufflink221/:/home/huw/program/tophat213/:/home/huw/program/bowtie2/:/home/huw/program/bin:/home/huw/program/blat:/home/huw/program/weblogo:/home/huw/program/gs/bin:/home/huw/program/ngsplot:$PATH'
JAVA_HOME 		= '/home/huw/program/jdk7/bin/java'
JAVA_HOME 		= '/opt/common/CentOS_6-dev/java/jdk1.8.0_31/'
NGSPLOT 		= '/home/huw/program/ngsplot'
TMPDIR			= '/scratch/huw'

if(species == 'mm9') {
	genomeFasta = MM9_FASTA
	genomeBWA =  MM9_BWA_INDEX
	genomeFAI = MM9_FAI
	DB_SNP = "" 
	FP_INT = "" 
	FP_TG  = ""
}

if(species == 'mm10_custom') {
	genomeFasta = MM10_CUSTOM_FASTA
	genomeBWA =  MM10_BWA_INDEX
	DB_SNP = paste0(dataDir, "/mm10/mm10_snp142.vcf") 
	FP_INT = "" 
	FP_TG  = ""
}

if(species == 'mm10') {
	genomeFasta = MM10_FASTA
	genomeBWA =  MM10_BWA_INDEX
	DB_SNP = paste0(dataDir, "/mm10/mm10_snp142.vcf") 
	FP_INT = "" 
	FP_TG  = ""
}

if(species == 'hybrid') {
	genomeFasta = B37_MM10_HYBRID_FASTA
	#genomeBWA =  B37_BWA_INDEX
	DB_SNP = paste0(dataDir, "/b37/dbsnp_138.b37.vcf")
	FP_INT = paste0(dataDir, "/b37/Agilent51MBExome__b37__FP_intervals.list")
	FP_TG  = paste0(dataDir, "/b37/Agilent51MBExome__b37__FP_tiling_genotypes.txt")
}

if(grepl('hg19', species, ignore.case=T)){
	   genomeFasta = HG19_FASTA
	   #bismarkHg19
	   genomeFAI = HG19_FAI
}

## configure files
if(assay == 'dnaseq' & assay.step == 'pipeline'){
	source('~/program/fun/param.r')
	pairfile = paste0(resDir, '/', pre, '_sample_pairing.txt')
	mapfile = paste0(resDir, '/', pre, '_sample_mapping.txt')
	groupfile = paste0(resDir, '/', pre, '_sample_grouping.txt')

	statsDir = paste0(resDir, '/stats')
	system(paste0('mkdir -p ', statsDir)) ## err std files

	if(!file.exists(groupfile)){
		stop('groupfile not exist.')
	}

	if(!file.exists(mapfile)){
		stop('mapfile not exist.')
	}

	if(!exists('species')){
		cat('Species must provided\n\n')
		cat('Environment parameters are not imported!!\n\n')
		return
	}

	if(!exists('wesImpactTarget')){
		cat('wesImpactTarget must provided\n\n')
		cat('Environment parameters are not imported!!\n\n')
		return
	}

	wesImpactTarget = 'AgilentExon_51MB_b37_v3' # baits
	initDir = paste0(resDir, '/initFiles')
	progressDir = paste0(resDir, '/progress')
	matricsDir = paste0(resDir, '/matrics')
	varDir = paste0(resDir, '/variation')
	haploDir = paste0(varDir, '/haplotypecaller')
	mutectDir = paste0(varDir, '/mutect')
	sniperDir = paste0(varDir, '/somaticsniper')
	facetsDir = paste0(varDir, '/facets')
	strvarDir = paste0(varDir, '/strvar')
	system(paste0('mkdir -p ', resDir)) ## r_fang
	system(paste0('mkdir -p ', progressDir)) ## jobname.done
	system(paste0('mkdir -p ', matricsDir))  ## collected metrics
	system(paste0('mkdir -p ', varDir))  
	system(paste0('mkdir -p ', haploDir)) 
	system(paste0('mkdir -p ', mutectDir)) 
	system(paste0('mkdir -p ', sniperDir)) 
	system(paste0('mkdir -p ', facetsDir)) 
	system(paste0('mkdir -p ', facetsDir, '/sample')) 
	system(paste0('mkdir -p ', strvarDir)) 
	system(paste0('mkdir -p ', initDir)) 

	wesImpactTarget = 'AgilentExon_51MB_b37_v3' # baits
	baits_ilist = paste0(targetsDir, '/', wesImpactTarget, '/', wesImpactTarget, '_baits.ilist')
	targets_ilist = paste0(targetsDir, '/', wesImpactTarget, '/', wesImpactTarget, '_targets.ilist')
	targets_bed = paste0(targetsDir, '/', wesImpactTarget, '/', wesImpactTarget, '_targets.bed')
	targets5bp_ilist = paste0(targetsDir, '/', wesImpactTarget, '/', wesImpactTarget, '_targets_plus5bp.ilist')
	targets5bp_bed = paste0(targetsDir, '/', wesImpactTarget, '/', wesImpactTarget, '_targets_plus5bp.bed')
}

## this for set up @ENV
path			= '/home/huw/program/homer/bin:/home/huw/local/bin:/home/huw/perl5/CPAN/bin:/home/huw/program/cufflink221/:/home/huw/program/tophat213/:/home/huw/program/bowtie2/:/home/huw/program/bin:/home/huw/program/blat:/home/huw/program/weblogo:/home/huw/program/gs/bin:/home/huw/program/ngsplot:$PATH'
JAVA_HOME 		= '/home/huw/program/jdk7/bin/java'
JAVA_HOME 		= '/opt/common/CentOS_6-dev/java/jdk1.8.0_31/'
NGSPLOT 		= '/home/huw/program/ngsplot'
TMPDIR			= '/scratch/huw'

if(species == 'mm9') {
	genomeFasta = MM9_FASTA
	genomeBWA =  MM9_BWA_INDEX
	genomeFAI = MM9_FAI
	DB_SNP = "" 
	FP_INT = "" 
	FP_TG  = ""
}

if(species == 'mm10_custom') {
	genomeFasta = MM10_CUSTOM_FASTA
	genomeBWA =  MM10_BWA_INDEX
	DB_SNP = paste0(dataDir, "/mm10/mm10_snp142.vcf") 
	FP_INT = "" 
	FP_TG  = ""
}

if(species == 'mm10') {
	genomeFasta = MM10_FASTA
	genomeBWA =  MM10_BWA_INDEX
	DB_SNP = paste0(dataDir, "/mm10/mm10_snp142.vcf") 
	FP_INT = "" 
	FP_TG  = ""
}

if(species == 'hybrid') {
	genomeFasta = B37_MM10_HYBRID_FASTA
	#genomeBWA =  B37_BWA_INDEX
	DB_SNP = paste0(dataDir, "/b37/dbsnp_138.b37.vcf")
	FP_INT = paste0(dataDir, "/b37/Agilent51MBExome__b37__FP_intervals.list")
	FP_TG  = paste0(dataDir, "/b37/Agilent51MBExome__b37__FP_tiling_genotypes.txt")
}

if(grepl('hg19', species, ignore.case=T)){
	   genomeFasta = HG19_FASTA
	   genomeBWA = HG19_BWA_INDEX
	   genomeFAI = HG19_FAI
	   DB_SNP = paste0(dataDir, "/hg19/dbsnp_138.hg19.vcf")
	   FP_INT = paste0(dataDir, "/hg19/Agilent51MBExome__hg19__FP_intervals.list")
	   FP_TG  = paste0(dataDir, "/hg19/Agilent51MBExome__hg19__FP_tiling_genotypes.txt")
	   ExAC_VCF = paste0(VEP, "/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz");
	   FACETS_DB_SNP = paste0(dataDir, "/hg19/dbsnp_137.hg19__RmDupsClean__plusPseudo50__DROP_SORT.vcf");
	   MILLS_1000G = paste0(dataDir, "/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf");
	   HAPMAP = paste0(dataDir, "/hg19/hapmap_3.3.hg19.vcf");
	   OMNI_1000G = paste0(dataDir, "/hg19/1000G_omni2.5.hg19.vcf");
	   PHASE1_SNPS_1000G = paste0(dataDir, "/hg19/1000G_phase1.snps.high_confidence.hg19.vcf");
	   COSMIC = paste0(dataDir, "/hg19/CosmicCodingMuts_v67_20131024.vcf");
}

if(grepl('b37', species, ignore.case=T)){
	genomeFasta = B37_FASTA
	genomeBWA =  B37_BWA_INDEX
	genomeFAI = B37_FAI
	DB_SNP = paste0(dataDir, "/b37/dbsnp_138.b37.vcf")
	FP_INT = paste0(dataDir, "/b37/Agilent51MBExome__b37__FP_intervals.list")
	FP_TG  = paste0(dataDir, "/b37/Agilent51MBExome__b37__FP_tiling_genotypes.txt")
	ExAC_VCF = paste0(VEP, "/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz");
	FACETS_DB_SNP = paste0(dataDir, "/b37/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf");
	MILLS_1000G = paste0(dataDir, "/b37/Mills_and_1000G_gold_standard.indels.b37.vcf");
	HAPMAP = paste0(dataDir, "/b37/hapmap_3.3.b37.vcf");
	OMNI_1000G = paste0(dataDir, "/b37/1000G_omni2.5.b37.vcf");
	PHASE1_SNPS_1000G = paste0(dataDir, "/b37/1000G_phase1.snps.high_confidence.b37.vcf");
	COSMIC = paste0(dataDir, "/b37/CosmicCodingMuts_v67_b37_20131024__NDS.vcf");
	COSMIC_HOTSPOTS = paste0(dataDir, "/b37/dmp_cosmic_for_hotspots.vcf");
	vep.species = 'homo_sapiens'
}

exe.jobs = function(jobnames=NULL, logger=logger, outputs=NULL, size=2000){
	if(is.null(jobnames) & !is.vector(jobnames)){
		msg = "jobnames parameter must vector of character"
		fatal(logger, msg)
		stop(msg)
	}
	for(i in 1:length(jobnames)){
		if(is.null(outputs)){ # run job without check
			jobname = jobnames[i]
			info(logger, jobname)
			ret = system(jobname)
			info(logger, paste0('return value for the above cmd: ', ret))
		}else if(!is.null(outputs[i]) & file.exists(outputs[i]) & file.size(output[i]) > size){ # no run if output exist
			msg = paste0('Output file ', outputs[i], ' exist, job not submitted for this')
			info(logger, msg)
		}else{ # either output file not exists or is too small
			jobname = jobnames[i]
			info(logger, jobname)
			ret = system(jobname)
			info(logger, paste0('return value for the above cmd: ', ret))
		}
	}
}

bjoblist = function(){
	tmpfile = tempfile()
	system(paste0('bjobs -w > ', tmpfile), wait=T)
	if(file.exists(tmpfile) & file.size(tmpfile) > 5 ) {
		joblist = suppressWarnings(fread(tmpfile))
		system(paste0('rm -f ', tmpfile))
		return(joblist$V7)
	}else{
		system(paste0('rm -f ', tmpfile))
		return('')
	}
}

wait4jobs = function(jobnames, logger){
	tmpfile = tempfile()
	system(paste0('bjobs -w > ', tmpfile), wait=T)
	joblist = 1
	if(file.exists(tmpfile) & file.size(tmpfile) > 5 ) joblist = fread(tmpfile)
	while(joblist != 1 & length(intersect(joblist$V7, jobnames)) > 0){
		info(logger, 'waiting for the below jobs to finish')
		info(logger, paste(intersect(joblist$V7, jobnames), collapse = ' '))
		Sys.sleep(100);
	}
	system(paste0('rm -f ', tmpfile))
}

#>Parameters
isPE			= '1'
clusters		= 'LSF'
#>Versions

#>Convention
##name1fastq		= 'base_xx_R1.fastq.gz'
##name2fastq		= 'base_xx_R2.fastq.gz'
##name1mrgedfq		= 'base_R1.fq.gz'
##name2mrgedfq		= 'base_R2.fq.gz'
##name1trimmedfq		= 'base_ptrimmed_R1.fq.gz'
##name2trimmedfq		= 'base_ptrimmed_R2.fq.gz'
##sai1			= 'base_R1.sai'
##sai2			= 'base_R2.sai'
##sam			= 'base_genome.sam'
##bam			= 'base_genome'
##bamext			= 'base_genome.bam'
##bamsorted		= 'base_genome_sorted'
##bamsortedext		= 'base_genome_sorted.bam'
##bamsortedrmdup		= 'base_genome_sorted_rmdup.bam'
##htseqcount		= 'base.htseq.count'
##bigwig			= 'base_genome_sorted_rmdup_10m.bw'
##MACS14_peaks_file 	= 'base_genome_macs14_peaks.bed'
##MACS2_peaks_file	= 'base_genome_macs2_peaks.xls'
##MACS2_broad_peaks_file  = 'base_genome_macs2broad_peaks.broadPeak'
##peakannoout		= 'base_genome_macs2_peaks_anno.txt'
##mcspeaks		= 'base_genome_macs2_peaks.narrowPeak'
##mcsbroadpeaks		= 'base_genome_macs2broad_peaks.broadPeak'
##peakannout		= 'base_genome_macs2_peaks_anno.txt'
##peakannout2		= 'base_genome_macs2broad_peaks_anno.txt'
##peakannoFiltered	= 'base_genome_macs2_peaks_anno_filtered.txt'
##peakannoFiltered2	= 'base_genome_macs2broad_peaks_anno_filtered.txt'
##ngsplot_file		= 'ngsplot_base.tss.avgprof.pdf'
##ngsplot_file_heatmap	= 'ngsplot_base.tss.heatmap.pdf'
##TSSCounts		= 'base_genome_tss2kb.bed'
##tssout			= 'ngsplot.base.tss'
##genebodyout		= 'ngsplot.base.genebody'
##enhancerout		= 'ngsplot.base.enhancer'
##motifout		= 'motif'
##starbam        		= 'base.Aligned.sortedByCoord.out.bam'
##
##>Checkpoints
##mergefastq		= 'name1mrgedfq'
##trimfq			= 'name1trimmedfq'
##alignStarAlign		= 'bamsortedext'
##alignBowtie2		= 'bamsrotedext'
##alignBwaAln		= 'sai1'
##alignBwaSampe		= 'sam'
##sortSam			= 'bamsortedext'
##sam2bam			= 'bamsortedext'
##removeDup		= 'bamsortedrmdup'
##bigwig			= 'bigwig'
##indexBam		= 'bamInput_bai'
##MACS2			= 'mcspeaks'
##MACS2BroadPeak		= 'mcsbroadpeaks'
##annotatePeaksHomer	= 'peakannoout_peakannoout2'
##annotateFilter		= 'peakannoFiltered_peakannoFiltered2'
##macs14			= 'mcspeaks'
##TSSCounts		= 'TSSCounts'
##Ngsplot			= 'tssout.avgprof.pdf'
##motifHomerFinder	= 'motifout'
##htseqcount		= 'htseqcount'

checkfiles = function(files, sizes){
	if(!(all(file.exists(files))) & !(all(file.size(files) > sizes))){
		msg = paste0('files not passed check: ', paste(files, collapse=' '))
		fatal(logger, msg)
		stop(msg)
	}
}

bsub.head = function(jobname, mem, cpu, We='', cwd='', postdone=NULL, statsDir = ''){
	if(We == '') We = '1:55'
	if(cwd == '') cwd = getwd()
	if(statsDir == ''){
		if(dir.exists(paste0(resDir, '/stats'))){
			statsDir = paste0(resDir, '/stats')
		}else{
			statsDir = '.'
		}
	}
	cmd = paste0(BSUB, " -J ", jobname, " -e ", statsDir, '/', jobname, ".err -o ", statsDir, '/', jobname, ".std ")
	cmd = paste0(cmd, " -cwd ", cwd, ' -We ', We, ' -R "rusage[mem=', mem, ']" -R "rusage[iounits=0]" -n ', cpu, ' ')
	if(!missing(postdone)){
		postdone = unlist(strsplit(postdone, split=" "))
		runningjobs = bjoblist()
		for( i in 1:length(postdone)){
			if(length(intersect(runningjobs, postdone[i])) > 0){
				cmd = paste0(cmd, ' -w "post_done(', postdone[i], ')" ')
			}
		}
	}
	cmd
}

## below could import from ~/program/fun/param.r
simple_table = data.table(
			  simple_val = c('Truncating', 'Truncating', 'Inframe', 'Inframe', 
					 'Missense_Mutation', 'Truncating', 'Inframe', 
					 'Inframe', 'Inframe', 'Inframe', 'Inframe'),
			  simple_name = c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 
					  'Missense_Mutation','Nonsense_Mutation', 'Nonstop_Mutation', 
					  'Silent', 'Splice_Region', 'Splice_Site', 'Translation_Start_Site')
			  )
simple_table
setkey(simple_table, 'simple_name')
color_table  = data.table(
			 color_index = c("Missense_Mutation Oncogenic",
					 "Truncating Oncogenic",
					 "Inframe Oncogenic",
					 "Missense_Mutation Passage",
					 "Truncating Passage",
					 "Inframe Passage") ,
			 color_value = c("#bebebe", 
					 "#000000", 
					 "#993404", 
					 adjustcolor("#bebebe", alpha.f=0.6), 
					 adjustcolor("#000000", alpha.f=0.6),
					 adjustcolor("#993404", alpha.f=0.6)) 
			 )
color_table
setkey(color_table, 'color_index')

rsem.fpkm = function(rsem){
	rsem$fpkm = rsem$counts * 10^9 / rsem$length
	rsem$fpkm = sweep(rsem$fpkm, 2, colSums(rsem$counts), '/')
	rsem
}

del.files = function(files, logger){
	if(all(file.exists(files))){
		for(i in 1:nrow(files)){
			info(logger, paste('deleting file:', files[i]))
			system(paste0('rm ', files[i]))
		}
	}else{
		info(logger, paste0('Not all files can be found. Nothing are deleted!'))
	}
}

