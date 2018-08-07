fastqc 			= '/home/huw/program/FastQC/fastqc'
hsaMethylBaitIntervalOrig='/home/huw/program/truseq-methyl-capture-epic-manifest-file.bed'
hsaMethylBaitInterval	='/home/huw/program/truseq-methyl-capture-epic-manifest-file-bait.interval'
hsaMethylTargetInterval	= '/home/huw/program/truseq-methyl-capture-epic-manifest-file-plusminus5bp.interval'

options(scipen=999)
system(paste0('cp /home/huw/program/truseq-methyl-capture-epic-header.txt ', hsaMethylBaitInterval))
tmp = fread(hsaMethylBaitIntervalOrig)
tmp = tmp[order(V1),]
library(tibble)
tmp = add_column(tmp, string = "+", .after = 3)
tmp[, V1:=sub("chr", "", V1)]
fwrite(tmp, file=hsaMethylBaitInterval, quote=F, sep="\t", col.names=F, append=T)
write.table(as.data.frame(tmp), file=hsaMethylBaitInterval, quote=F, sep="\t", col.names=F, append=T, row.names=F)
tmp[, V2:=V2-5]
tmp[, V3:=V3+5]
head(tmp)
tmp[, V1 := as.character(V1)]
tmp[, V2 := as.character(V2)]
tmp[, V3 := as.character(V3)]
tmp[, V4:=paste0("chr", V1, "_", V2, "-", V3)]
system(paste0('cp /home/huw/program/truseq-methyl-capture-epic-header.txt ', hsaMethylTargetInterval))
write.table(as.data.frame(tmp), file=hsaMethylTargetInterval, quote=F, sep="\t", col.names=F, append=T, row.names=F)
rm(tmp)

bismarkGRCh37		= '/ifs/work/solitlab/huw/study/db/hsa/grch37_bismark'
bismark			= '/home/huw/program/Bismark'
methylTarget		= '/home/huw/program/truseq-methyl-capture-epic-manifest-file.bed'
SOAPHLA			= '/home/huw/program/SOAP-HLA'
bamReadCount		= '/home/huw/local/bam-readcount'
IEDBMHC1		= '/home/huw/program/mhc_i'
IEDBMHC1		= '/home/huw/program/mhc_ii'
pvacseq			= '/home/huw/program/anaconda3/bin/pvacseq'
tmpdir			= '/scratch/huw'
polysolver		= '/home/huw/program/polysolver'
ponfile			= '/home/huw/program/pon56Sample.vcf.gz'
make.trinuc		= '/home/huw/program/make_trinuc_maf.py'
mutsig			= '/opt/common/CentOS_6-dev/mutsig/cv-1.4/'
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
samtoolsdir 		= '/opt/common/CentOS_6-dev/samtools/samtools-1.3.1'
sambamba		= '/home/huw/program/sambamba' 
java			= '/home/huw/program/jdk7/bin/java'
java			= '/opt/common/CentOS_6-dev/java/jdk1.8.0_31/bin/java'
picard			= '/home/huw/program/picard.jar'
htseqcountexe		= '/home/huw/local/bin/htseq-count'
NGSPLOTexe		= '/home/huw/program/ngsplot/bin/ngs.plot.r'
NGSFilter		= '/home/huw/program/ngs-filters/'
QSUB			= '/common/sge/bin/lx24-amd64/qsub'
BSUB			= '/common/lsf/9.1/linux2.6-glibc2.3-x86_64/bin/bsub'
trim_galore		= '/home/huw/program/trim_galore'
bowtie2                 = '/home/huw/program/bowtie2/bowtie2'
bowtie2dir              = '/home/huw/program/bowtie2'
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
FACETS_LIB		= '/opt/common/CentOS_6-dev/facets_lib/facets-0.5.6'
FIXMULTIINDEL		= '/opt/common/CentOS_6/FixMultiInDel/FixMultiInDel-2.0.1'
GATK			= '/home/huw/program/GenomeAnalysisTK-3.8-0-ge9d806836'
GATK4			= '/home/huw/program/gatk-4.0.2.1/gatk'
#GATK4			= '/home/huw/program/gatk-4.0.2.1/gatk-package-4.0.2.1-local.jar'
gnomadfile			= '/home/huw/af-only-gnomad.raw.sites.b37.vcf.gz'


JAVA			= '/opt/common/CentOS_6/java/jdk1.8.0_31/bin'
JAVA7_MUTECT            = '/opt/common/CentOS_6/java/jdk1.7.0_75/bin'
MCR			= '/opt/common/CentOS_6/matlab/v81'
MUTECT			= '/opt/common/CentOS_6/muTect/muTect-1.1.7'
PERL			= '/opt/common/CentOS_6/perl/perl-5.22.0/bin'
PERL			= '/home/huw/local/perl/perl-5.26.0/bin'
PICARD			= '/home/huw/program/picard.jar'
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
VEP                     = '/opt/common/CentOS_6-dev/vep/v86'
VEPv88                  = '/opt/common/CentOS_6-dev/vep/v88/variant_effect_predictor.pl'
VEP_plugins             = '/home/huw/program/VEP_plugins'
VIRMID			= '/opt/common/CentOS_6/virmid/Virmid-1.1.1'
WES_FILTER		= '/opt/common/CentOS_6/wes-filter/wes-filters-1.1.3'
CMOBIN			= '/opt/common/CentOS_6-dev/python/python-2.7.10/bin/'


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

genomeGRCh37Fasta	= '/ifs/work/solitlab/huw/study/db/hsa/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
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
