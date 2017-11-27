library(data.table)
library(log4r)

#setwd(resDir)
cwd = getwd()

mapfile = 'Proj_07813_D_sample_mapping.txt'
groupfile = 'Proj_07813_D_sample_grouping.txt'
pairfile = 'Proj_07813_D_sample_pairing.txt'

resDir = 'r_fang' # result folder

species = 'B37' # genome
wesImpactTarget = 'AgilentExon_51MB_b37_v3' # baits
assay = 'wes'

delFiles = F # delete files after each step 

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

statsDir = paste0(resDir, '/stats')
initDir = paste0(resDir, '/initFiles')
progressDir = paste0(resDir, '/progrss')
matricsDir = paste0(resDir, '/matrics')
varDir = paste0(resDir, '/variation')
haploDir = paste0(varDir, '/haplotypecaller')
mutectDir = paste0(varDir, '/mutect')
sniperDir = paste0(varDir, '/somaticsniper')
facetsDir = paste0(varDir, '/facets')
strvarDir = paste0(varDir, '/strvar')
system(paste0('mkdir -p ', resDir)) ## r_fang
system(paste0('mkdir -p ', statsDir)) ## err std files
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

##logfile = paste0('_', resDir, '_log')
##if(file.exists(logfile)){
##	stop(paste0(logfile, 'exist, already running?'))
##}else{
##	create.logger() -> logger
##	logfile(logger) = logfile
##	level(logger) = 'INFO' 
##	#info(logger, 'message') others like warn(), error(), fatal()
##}


## this configure for the locations for common programs
## the task.conf under the working directory provide the tasks. 
## save the input and base as special words here
## _genome will sub with $genome
#>Locations
belly2			= '/home/huw/program/delly_v0.7.7_linux_x86_64bit'
exeDir			= '/ifs/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/'
targetsDir		= '/ifs/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/targets/'
facets			= '/ifs/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/facets'
vepDir		= '/ifs/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/vep'
dataDir			=  '/ifs/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/data'
ConvertQualityScore     = '/home/huw/program/ConvertQualityScore'
STAR			= '/home/huw/program/STAR-STAR_2.5.0a/bin/Linux_x86_64_static/STAR'
STAR_DIR		= '/home/huw/program/STAR-STAR_2.5.0a/bin/Linux_x86_64_static'
RSEM			= '/opt/common/CentOS_6/rsem/RSEM-1.2.25/'
KALLISTO		= '/home/huw/local/bin/kallisto '
##STAR			= '/opt/common/CentOS_6/star/STAR-STAR_2.5.0a/bin/Linux_x86_64/STAR'
samtools		= '/opt/common/CentOS_6-dev/bin/current/samtools'
sambamba		= '/home/huw/program/sambamba' 
java			= '/home/huw/program/jdk7/bin/java'
java			= '/opt/common/CentOS_6-dev/java/jdk1.8.0_31/bin/java'
picard			= '/home/huw/program/picard.jar'
htseqcountexe		= '/home/huw/local/bin/htseq-count'
NGSPLOTexe		= '/home/huw/program/ngsplot/bin/ngs.plot.r'
QSUB			= '/common/sge/bin/lx24-amd64/qsub'
BSUB			= '/common/lsf/9.1/linux2.6-glibc2.3-x86_64/bin/bsub'
trim_galore		= '/home/huw/program/trim_galore'
bowtie2                 = '/home/huw/program/bowtie2/bowtie2'
macs14                  = '/home/huw/local/bin/macs14'
macs2                  	= '/home/huw/local/bin/macs2'
homerFindMotif          = '/home/huw/program/homer/bin/findMotifsGenome.pl'
homerAnnotatePeaks	= '/home/huw/program/homer/bin/annotatePeaks.pl'
tssCounts		= '/home/huw/program/tsscounts.pl'
bwa			= '/home/huw/program/bwa-0.7.12/bwa '
bedtools		= '/home/huw/program/bedtools/bin'
bam2bigwig		= '/home/huw/program/bam2bigwig.sh'
annotationFilter	= '/home/huw/program/annotationFilter.pl'

enhancerMm10Bed		= '/ifs/work/solitlab//huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.enhancer.bed'
enhancerMm9Bed          = '/home/program/PLtest/Enhancers_mm9_CreyghtonPNAS/PutativeEnhancers_mm9_5types_SortMerge300_noPromOL_GT50bp.txt'
tssMm9Bed               = '/home/huw/program/PLtest/refFlat_mm9_TSSpm2kbSORTED.txt'
tssMm10Bed		= '/ifs/work/solitlab//huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.chr.symbol.tss2kb.bed'
genomeMm9Bowtie2	= '/ifs/work/solitlab//huw/study/db/mmu/mm9-bowtie2/mm9'
genomeMm9StarDir	= '/ifs/work/solitlab//huw/study/db/mmu/mm9_star'
genomeMm9GTF		= '/ifs/work/solitlab//huw/study/db/mmu/dna/Mus_musculus.NCBIM37.67.chr.gtf'
genomeMm10Bowtie2	= '/ifs/work/solitlab//huw/study/db/mmu/mm10-bowtie2/mm10'
genomeMm10StarDir	= '/ifs/work/solitlab//huw/study/db/mmu/mm10star'
genomeMm10GTF		= '/ifs/work/solitlab//huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.chr.gtf'
mm10chromsize		= '/ifs/work/solitlab//huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.chrom.size'

ABRA			= '/opt/common/CentOS_6/abra2/abra2-2.07'
BCFTOOLS		= '/opt/common/CentOS_6/bcftools/bcftools-1.2/bin'
BEDTOOLS		= '/opt/common/CentOS_6/bedtools/bedtools-2.22.0/bin'
BWA			= '/opt/common/CentOS_6/bwa/bwa-0.7.12/bwa'
CUTADAPT		= '/opt/common/CentOS_6/cutadapt/cutadapt-1.9.1/bin/cutadapt'
DELLY			= '/opt/common/CentOS_6/delly/delly_v0.6.1'
dRANGER			= '/opt/common/CentOS_6/dranger/dRanger_annotate_v2.0'
FACETS_SUITE		= '/home/huw/program/facets-suite-1.0.1'
FACETS_LIB		= '/opt/common/CentOS_6/facets_lib/facets-0.3.30'
FIXMULTIINDEL		= '/opt/common/CentOS_6/FixMultiInDel/FixMultiInDel-2.0.1'
GATK			= '/home/huw/program/GenomeAnalysisTK-3.8-0-ge9d806836'
JAVA			= '/opt/common/CentOS_6/java/jdk1.8.0_31/bin'
JAVA7_MUTECT            = '/opt/common/CentOS_6/java/jdk1.7.0_75/bin'
MCR			= '/opt/common/CentOS_6/matlab/v81'
MUTECT			= '/opt/common/CentOS_6/muTect/muTect-1.1.7'
PERL			= '/opt/common/CentOS_6/perl/perl-5.22.0/bin'
PICARD			= '/opt/common/CentOS_6/picard/picard-tools-1.124'
PYTHON			= '/opt/common/CentOS_6/python/python-2.7.8/bin/python'
R			= '/opt/common/CentOS_6/R/R-3.1.2/bin'
SAMTOOLS                = '/opt/common/CentOS_6/samtools/samtools-1.2'
SCALPEL			= '/opt/common/CentOS_6/scalpel/scalpel-0.2.2'
SOMATICSNIPER		= '/opt/common/CentOS_6/somaticsniper/somatic-sniper-1.0.4'
STRELKA			= '/opt/common/CentOS_6/strelka/strelka_1.0.11'
TABIX			= '/opt/common/CentOS_6/samtools/samtools-1.2/htslib-1.2.1'
VARSCAN			= '/opt/common/CentOS_6/varscan/v2.3.7'
VCF2MAF			= '/home/huw/program/vcf2maf-1.6.14'
VCFTOOLS		= '/opt/common/CentOS_6-dev/vcftools/v0.1.14/bin'
VEP                     = '/opt/common/CentOS_6/vep/v86'
VIRMID			= '/opt/common/CentOS_6/virmid/Virmid-1.1.1'
WES_FILTER		= '/opt/common/CentOS_6/wes-filter/wes-filters-1.1.3'


B37_BWA_INDEX		= '/ifs/depot/assemblies/H.sapiens/b37/index/bwa/0.7.12/b37.fasta'
B37_FASTA		= '/ifs/depot/assemblies/H.sapiens/b37/b37.fasta'
B37_FAI			= '/ifs/depot/assemblies/H.sapiens/b37/b37.fasta.fai'
B37_MM10_HYBRID_BWA_INDEX	= '/ifs/depot/assemblies/hybrids/H.sapiens_M.musculus/b37_mm10/index/bwa/0.7.12/b37_mm10.fasta'
B37_MM10_HYBRID_FASTA		= '/ifs/depot/assemblies/hybrids/H.sapiens_M.musculus/b37_mm10/b37_mm10.fasta'
B37_MM10_HYBRID_FAI	= '/ifs/depot/assemblies/hybrids/H.sapiens_M.musculus/b37_mm10/b37_mm10.fasta.fai'
HG19_BWA_INDEX		= '/ifs/depot/assemblies/H.sapiens/hg19/index/bwa/0.7.4-r385/hg19.fasta'
HG19_FASTA		= '/ifs/depot/assemblies/H.sapiens/hg19/hg19.fasta'
HG19_FAI		= '/ifs/depot/assemblies/H.sapiens/hg19/hg19.fasta.fai'
HG19_MM10_HYBRID_BWA_INDEX	= '/ifs/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/bwa/0.7.8/hg19_mm10.fasta'
HG19_MM10_HYBRID_FASTA	= '/ifs/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/hg19_mm10.fasta'
HG19_MM10_HYBRID_FAI	= '/ifs/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/hg19_mm10.fasta.fai'
MM9_BWA_INDEX		= '/ifs/depot/assemblies/M.musculus/mm9/index/bwa/0.7.4-r385/mm9.fasta'
MM9_FASTA		= '/ifs/depot/assemblies/M.musculus/mm9/mm9.fasta'
MM9_FAI			= '/ifs/depot/assemblies/M.musculus/mm9/mm9.fasta.fai'
MM10_BWA_INDEX		= '/ifs/depot/assemblies/M.musculus/mm10/index/bwa/0.7.8/mm10.fasta'
MM10_FASTA		= '/ifs/depot/assemblies/M.musculus/mm10/mm10.fasta'
MM10_FAI		= '/ifs/depot/assemblies/M.musculus/mm10/mm10.fasta.fai'
clipR1			= 'AGATCGGAAGAGCACACGTCT'
clipR2			= 'AGATCGGAAGAGCACACGTCT'
bqtrim			= '3'

## for homer annotatePeaks.pl mm10
homerMm10AnnotationFile	= '/ifs/work/solitlab//huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.chr.gtf.annotations.final.txt'

genomeHg19Bwa		= '/ifs/work/solitlab//huw/study/db/hsa/hg37-bwa/hg37 '
genomeHg19StarDir	= '/ifs/work/solitlab//huw/study/db/hsa/hg37star'
genomeHg19Bowtie2	= '/ifs/work/solitlab//huw/study/db/hsa/hg37bowtie2/hg37'
genomeHg38Bowtie2	= '/ifs/work/solitlab//huw/study/db/hsa/hg38_bowtie2/hg38'
genomeHg19GTF		= '/ifs/work/solitlab//huw/study/db/hsa/dna/Homo_sapiens.GRCh37.75.gtf'
hg19chromsize		= '/ifs/work/solitlab//huw/study/db/hsa/dna/Homo_sapiens.GRCh37.75.chrom.size'

genomeGRCh37StarDir	= '/ifs/depot/assemblies/H.sapiens/b37/index/star/2.4.1d/gencode/v18/overhang49 /ifs/work/solitlab/huw/park/study/hiseq/chipseq/template.conf'

genomeGRCh38GTF		= '/ifs/work/solitlab/huw/study/db/hsa/dna/Homo_sapiens.GRCh38.88.gtf'
genomeGRCh38Fasta	= '/ifs/work/solitlab/huw/study/db/hsa/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
genomeGRCh38StarDir	= '/ifs/work/solitlab/huw/study/db/hsa/hg38_star'
genomeGRCh38Rsem	= '/ifs/work/solitlab/huw/program/RSEM_tutorial/ref/GRCh38/GRCh38'
RSEM_REF		= '/ifs/work/solitlab/huw/program/RSEM_tutorial/ref/GRCh38/GRCh38'

genomeGRCh38Kallisto	= '/ifs/work/solitlab/huw/study/db/hsa/hg38_kallisto/kallisto_GRCh38'

genomeGRCm38Fasta	= '/ifs/work/solitlab/huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.dna.primary_assembly.chr.fa'
genomeGRCm38GTF		= '/ifs/work/solitlab/huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.chr.gtf'
genomeGRCm38Rsem	= '/ifs/work/solitlab/huw/program/RSEM_tutorial/ref/GRCm38/GRCm38'

## this for set up @ENV
path			= '/home/huw/program/homer/bin:/home/huw/local/bin:/home/huw/perl5/CPAN/bin:/home/huw/program/cufflink221/:/home/huw/program/tophat213/:/home/huw/program/bowtie2/:/home/huw/program/bin:/home/huw/program/blat:/home/huw/program/weblogo:/home/huw/program/gs/bin:/home/huw/program/ngsplot:$PATH'
JAVA_HOME 		= '/home/huw/program/jdk7/bin/java'
JAVA_HOME 		= '/opt/common/CentOS_6-dev/java/jdk1.8.0_31/'
NGSPLOT 		= '/home/huw/program/ngsplot'
TMPDIR			= '/ifs/e63data/bergerm1/Resources/TMPDIR/'

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

baits_ilist = paste0(targetsDir, '/', wesImpactTarget, '/', wesImpactTarget, '_baits.ilist')
targets_ilist = paste0(targetsDir, '/', wesImpactTarget, '/', wesImpactTarget, '_targets.ilist')
targets_bed = paste0(targetsDir, '/', wesImpactTarget, '/', wesImpactTarget, '_targets.bed')
targets5bp_ilist = paste0(targetsDir, '/', wesImpactTarget, '/', wesImpactTarget, '_targets_plus5bp.ilist')
targets5bp_bed = paste0(targetsDir, '/', wesImpactTarget, '/', wesImpactTarget, '_targets_plus5bp.bed')



if(grep('b37', species, ignore.case=T)){
	genomeFasta = B37_FASTA
	genomeBWA =  B37_BWA_INDEX
	genomeFAI = B37_FAI
	DB_SNP = paste0(dataDir, "/b37/dbsnp_138.b37.vcf")
	FP_INT = paste0(dataDir, "/b37/Agilent51MBExome__b37__FP_intervals.list")
	FP_TG  = paste0(dataDir, "/b37/Agilent51MBExome__b37__FP_tiling_genotypes.txt")
	REF_SEQ = B37_FASTA;
	REF_FAI = B37_FAI;
	BWA_INDEX = B37_BWA_INDEX;
	ExAC_VCF = paste0(VEP, "/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz");
	FACETS_DB_SNP = paste0(dataDir, "/b37/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf");
	MILLS_1000G = paste0(dataDir, "/b37/Mills_and_1000G_gold_standard.indels.b37.vcf");
	HAPMAP = paste0(dataDir, "/b37/hapmap_3.3.b37.vcf");
	OMNI_1000G = paste0(dataDir, "/b37/1000G_omni2.5.b37.vcf");
	PHASE1_SNPS_1000G = paste0(dataDir, "/b37/1000G_phase1.snps.high_confidence.b37.vcf");
	COSMIC = paste0(dataDir, "/b37/CosmicCodingMuts_v67_b37_20131024__NDS.vcf");
	COSMIC_HOTSPOTS = paste0(dataDir, "/b37/dmp_cosmic_for_hotspots.vcf");
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

deldir = function(ddir=NULL){
	if(!is.null(ddir) & dir.exists(ddir)){
		msg = paste0('deleting folder ', ddir)
		warn(logger, msg)
		system(paste0('rm -rf ', ddir))
		return(0)
	}
	return(1)
}

delfile.2 = function(dfile=NULL, logger){
	if(!is.null(dfile) & file.exists(dfile)){
		msg = paste0('deleting file ', dfile)
		warn(logger, msg)
		system(paste0('rm -f ', dfile))
		return(0)
	}
	return(1)
}

delfile = function(dfile=NULL, logger, cfile=NULL, s){
	if(!is.null(dfile) & !file.exists(dfile)){
		msg = paste0(dfile, ' not exists')
		info(logger, msg);stop(msg)
	}
	if(!is.null(cfile) & !file.exists(cfile)){
		msg = paste0(cfile, ' not exists')
		info(logger, msg);stop(msg)
	}else if(file.size(cfile) < s){
		msg = paste0(cfile, ' is too small')
		info(logger, msg);stop(msg)
	}else{
		msg = paste0('deleting file ', dfile)
		warn(logger, msg)
		system(paste0('rm -f ', dfile))
		return(0)
	}
	return(1)
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

##REF_SEQ = "B37_MM10_HYBRID_FASTA";
##REF_FAI = "B37_MM10_HYBRID_FAI";
##BWA_INDEX = "B37_MM10_HYBRID_BWA_INDEX";
##ExAC_VCF = "VEP/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz";
##DB_SNP = "Bin/data/b37/dbsnp_138.b37.excluding_sites_after_129.vcf";
##FACETS_DB_SNP = "Bin/data/b37/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf";
##MILLS_1000G = "Bin/data/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";
##HAPMAP = "Bin/data/b37/hapmap_3.3.b37.vcf";
##OMNI_1000G = "Bin/data/b37/1000G_omni2.5.b37.vcf";
##PHASE1_SNPS_1000G = "Bin/data/b37/1000G_phase1.snps.high_confidence.b37.vcf";
##COSMIC = "Bin/data/b37/CosmicCodingMuts_v67_b37_20131024__NDS.vcf";
##COSMIC_HOTSPOTS = "Bin/data/b37/dmp_cosmic_for_hotspots.vcf";

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

bsub.head2 = function(jobname, mem, cpu, We, cwd, postdone, submit, resDir = 'r_fang'){
	cmd = paste0(BSUB, " -J ", jobname, " -e ", statsDir, '/', jobname, ".err -o ", statsDir, '/', jobname, ".std ")
	cmd = paste0(cmd, " -cwd ", cwd, ' -We ', We, ' -R "rusage[mem=', mem, ']" -R "rusage[iounits=0]" -n ', cpu)
	if(!missing(postdone)){
		postdone = unlist(strsplit(postdone, split=" "))
		for( i in 1:length(postdone)){
			cmd = paste0(cmd, ' -w "post_done(', postdone[i], ')" ')
		}
	}
	if(!missing(submit)){
		submit2 = unlist(strsplit(submit, " "))
		submit3 = paste0('ls ', cwd, '/progress/', resDir, '/', submit2)
		submit4 = paste0(submit3, collapse = ' && ')
		cmd = paste0(cmd, ' -E "', submit4, '"')
	}
	cmd
}

bsub.head = function(jobname, mem, cpu, We, cwd, postdone, submit, resDir = 'r_fang'){
	cmd = paste0(BSUB, " -J ", jobname, " -e ", statsDir, '/', jobname, ".err -o ", statsDir, '/', jobname, ".std ")
	cmd = paste0(cmd, " -cwd ", cwd, ' -We ', We, ' -R "rusage[mem=', mem, ']" -R "rusage[iounits=0]" -n ', cpu)
	if(!missing(postdone)){
		postdone = unlist(strsplit(postdone, split=" "))
		runningjobs = bjoblist()
		for( i in 1:length(postdone)){
			if(length(intersect(runningjobs, postdone[i])) > 0){
				cmd = paste0(cmd, ' -w "post_done(', postdone[i], ')" ')
			}
		}
	}
	if(!missing(submit)){
		submit2 = unlist(strsplit(submit, " "))
		submit3 = paste0('ls ', cwd, '/progress/', resDir, '/', submit2)
		submit4 = paste0(submit3, collapse = ' && ')
		cmd = paste0(cmd, ' -E "', submit4, '"')
	}
	cmd
}

checkfiles = function(files, sizes){
	if(!(all(file.exists(files))) & !(all(file.size(files) > sizes))){
		msg = paste0('files not passed check: ', paste(files, collapse=' '))
		fatal(logger, msg)
		stop(msg)
	}
}
