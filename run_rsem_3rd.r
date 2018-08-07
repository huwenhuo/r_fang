# exome sequence vlad
# exome sequence
## cnv  gene specific
## facets plot gene specific t2 vs t1
## gene level by t2 vs t1
# rnaseq res.david
## res.david heatmap
## res.david gsea
## res.david dotplot
## immuno convolution analysis based on david
## rnaseq based on david
# for grant 12/21/2017
## rnaseq based on rsem, the whole pipeline
## rnaseq based on rsem
## rnaseq my count results with 11 pairs
## top 1000 variant genes, pca
## top 1000 variant genes heatmap
## top 100 bladder genes 
##immune gene target heatmap
##basal luminal classifier genes
## concencus clustering
## pamr
library(log4r)
source('~/program/fun/sync.r')

sync = function(loc1, loc2){ ## sync('exon', 'maftools')
	base.selene = '/ifs/work/solitlab/huw/solit/study/hiseq/blca_cmo_06155_2016/'
	base.mski1925 = 'mski1925:~/solit/study/hiseq/blca_cmo_06155_2016/'
	ll1  = paste0(base.selene, loc1, '/', loc2, '/')
	ll2  = paste0(base.mski1925, loc1, '/', loc2, '/')
	cmd = paste0('rsync -avur ', ll1, ' ', ll2)
	cat(cmd, "\n")
	system(cmd)
}

# exome sequence vlad
# neoantigen
vlad.neo = fread('exonseq/Results_Scc/neoantigens/SCC.AllBindersPan_wgenes.9.txt')
vlad.neo[, Dif_Rank := MT.H_Avg_Ranks - WT.H_Avg_Ranks]
vlad.neo[, pat.max.wt.score := max(WT.Score), by = list(Sample, CHROM_POS_REF_ALT)]
vlad.neo[, pat.max.mt.score := max(MT.Score), by = list(Sample, CHROM_POS_REF_ALT)]
vlad.neo.v2 = vlad.neo[!duplicated(Sample, CHROM_POS_REF_ALT),]
vlad.neo.v2
tmp = data.table::melt(vlad.neo.v2[, .(Sample, CHROM_POS_REF_ALT, pat.max.wt.score, pat.max.mt.score)])
tmp[, id := paste0(Sample, "_", CHROM_POS_REF_ALT, "_", toupper(substr(variable, 9, 10))) ]
tmp
tmp2 = data.table::dcast(tmp[, .(Sample, id, value)], Sample ~ id)
tmp2[is.na(tmp2)] = 0
rn = tmp2[,1]
tmp2 = tmp2[, -1]
tmp2[1:4,1:4]
tmp3 = as.data.frame(tmp2); tmp3
tmp3[1:4,1:4]
row.names(tmp3) = rn
tmp3 = as.matrix(tmp3)
tmp3[1:4,1:4]

pdf(file='exonseq/maftools/heatmap_neoantigen.pdf', width=10, height=12)
heatmap.2(tmp3,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')
#pheatmap(as.matrix(vlad.neo.v2.w), scale='none', cluster_col = F, cluster_row = F)
dev.off()
sync('exonseq', 'maftools')

## compare vlad with nick's
scc.maf[, tag:=paste0(Hugo_Symbol, ':', Start_Position)]
fread('exonseq/Results_Scc/variants/SCC.MAF.hg19.AutoFilt.txt') -> vlad.maf
fread('exonseq/Results_Scc/variants/SCC.MAF.hg19.merged.txt') -> vlad.maf
vlad.maf[, tag:=paste0(Hugo_Symbol, ':', Start_Position)]
vlad.maf[, Tumor_Sample_Barcode := sub(".bam", "", Tumor_Sample_Barcode)]
as.data.table(as.data.frame(table(vlad.maf$Tumor_Sample_Barcode))) -> cmp
as.data.table(as.data.frame(table(scc.maf$Tumor_Sample_Barcode))) -> tmp
cmp
merge(cmp, tmp, by.x = 'Var1', by.y='Var1') -> cmp; cmp
setnames(cmp, c('Freq.x', 'Freq.y'), c('chan', 'nick'))
cmp[, dif := nick - chan]
cmp[, id := substr(Var1, 10, 12)]
cmp[, chan.t12 := 0]
cmp[, nick.t12 := 0]
for(i in 1:nrow(cmp)){
	id = cmp[i, id]
	id1 = paste0('s_DS_bla_', id, '_T1')
	id2 = paste0('s_DS_bla_', id, '_T2')
	cmp[i, nick.t12 := length(intersect(scc.maf[Tumor_Sample_Barcode == id1,tag], scc.maf[Tumor_Sample_Barcode == id2,tag]))]
	cmp[i, vlad.t12 := length(intersect(vlad.maf[Tumor_Sample_Barcode == id1,tag], vlad.maf[Tumor_Sample_Barcode == id2,tag]))]
}
tmp = vlad.maf
colnames(vlad.maf)
vlad.maf[, isMutect := F]
vlad.maf[grep('MT', Caller) & Taf >  0.05, ]
tmp = vlad.maf[grep('MT', Caller), ]
tmp = tmp[ Taf > 0.05 & Tcov > 3,]
tmp[, Var1 := sub(".bam", "", Var1)]
tmp
cmp
setkey(tmp, 'Var1')
setkey(cmp, 'Var1')
cmp = merge(cmp, tmp)
setnames(cmp, 'Freq', 'chan.mutect')
cmp
tmp2 = as.data.table(as.data.frame(table(tmp$Tumor_Sample_Barcode)))
tmp[
cbind(cmp, tmp$Freq)

table(vlad.maf$isMutect)
vlad.maf$Caller[1:10]

fread('exonseq/Results_Scc/variants/SCC.CNT.txt') -> vlad.cnt
sum(vlad.cnt$SNV_merged)
colnames(vlad.cnt)
vlad.maf

# exome sequence
setwd('..')
library(data.table)
library(autospy)
library(WriteXLS)
library(readxl)
library(ggplot2)
library(metafolio)
library(VennDiagram)
library(maftools)
library(Hmisc)

options(width=155)

## /ifs/res/share/solit/alahmadh/Proj_07813_DF
pid = "bla_[0-9]+_"
sid = "_[T|M][0-9]?$"

maffile = "exonseq/Proj_07813_DF/r_001/post/Proj_07813_DF___SOMATIC.vep.filtered.facets.V3.maf"
maffile = "exonseq/scc_oncokb.maf"
scc.maf = fread(maffile)
unique(scc.maf$Variant_Classification)
cn = c("Hugo_Symbol", "Chromosome", "Variant_Classification", "Variant_Type", "Reference_Allele", 
       "Tumor_Seq_Allele1", 'Tumor_Seq_Allele2','Tumor_Sample_Barcode', 't_depth', 't_ref_count', 't_alt_count',
       'n_depth', 'n_ref_count', 'n_alt_count', 'HGVSp_Short')
scc.maf[Variant_Classification == 'Frame_Shift_Del', cn, with=F]
colnames(scc.maf)

## mutsig for mc3
cmd = bsub.head('mutsig', mem=180, cpu=4, We='8:26')
cmd = paste0(cmd, ' "', "/opt/common/CentOS_6-dev/mutsig/cv-1.4/run_MutSigCV.sh /opt/common/CentOS_6-dev/matlab/R2013a/v81/ ../Proj_06230/tcga_mc3_20171217_oncokb_facets_lessc1_trinuc.tsv /opt/common/CentOS_6-dev/mutsig/cv-1.4/lib/exome_full192.coverage.txt /opt/common/CentOS_6-dev/mutsig/cv-1.4/lib/gene.covariates.txt mc3_mutsig /opt/common/CentOS_6-dev/mutsig/cv-1.4/lib/mutation_type_dictionary_file.txt /opt/common/CentOS_6-dev/mutsig/cv-1.4/lib/chr_files_hg19", '"')
cmd
system(cmd)

## mutsig
system(" /opt/common/CentOS_6-dev/mutsig/cv-1.4/run_MutSigCV.sh /opt/common/CentOS_6-dev/matlab/R2013a/v81/ scc_oncokb.maf /opt/common/CentOS_6-dev/mutsig/cv-1.4/lib/exome_full192.coverage.txt /opt/common/CentOS_6-dev/mutsig/cv-1.4/lib/gene.covariates.txt myresults /opt/common/CentOS_6-dev/mutsig/cv-1.4/lib/mutation_type_dictionary_file.txt /opt/common/CentOS_6-dev/mutsig/cv-1.4/lib/chr_files_hg19")

fwrite(scc.maf, file='exonseq/scc_oncokb.maf', sep="\t", quote=F)

# missense putative driver   #008000
# missense putative passager #53d400
# truncating driver	     #000000
# truncating passager        #708090
# no alteration		     #bebebe
# duplication                #ff0000
# deletion                   #ff0000
# inframe driver	     #993404
# inframe passager           #fe9929
# fusion 		     #860009
simple_table = data.table(
			  simple_value = c('Truncating', 'Truncating', 'Inframe', 'Inframe', 
					 'Missense_Mutation', 'Truncating', 'Inframe', 
					 'Silent', 'Splice_Region', 'Inframe', 'Inframe'),
			  simple_name = c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 
					  'Missense_Mutation','Nonsense_Mutation', 'Nonstop_Mutation', 
					  'Silent', 'Splice_Region', 'Splice_Site', 'Translation_Start_Site')
			  )
simple_table[, simple_value := factor(simple_value, levels=c('Truncating', 'Missense_Mutation', 'Inframe'))]
simple_table

setkey(simple_table, 'simple_name')

color_table  = data.table(
	color_index = c("Missense_Mutation Oncogenic",
			 "Truncating Oncogenic",
			 "Inframe Oncogenic",
			 "Missense_Mutation Passage",
			 "Truncating Passage",
			 "Inframe Passage",
			 "AML",
			 "DEL",
			 "Multi_Hit"),
	color_value = c("forestgreen", 
			 "black", 
			 "brown4", 
			 adjustcolor("forestgreen", alpha.f=0.6), 
			 adjustcolor("black", alpha.f=0.6),
			 adjustcolor("brown4", alpha.f=0.6),
			 "red",
			 "blue",
			 "magenta3") )
color_table[, color_index := factor(color_index, levels=c("Missense_Mutation Oncogenic", "Truncating Oncogenic", "Inframe Oncogenic", "Missense_Mutation Passage", "Truncating Passage", "Inframe Passage", "AML", "DEL", 'Multi_Hit') )]
color_vector = color_table$color_value
names(color_vector) = color_table$color_index
color_vector

vc.key = c("Missense_Mutation Oncogenic", "Inframe Oncogenic", "Inframe Passage", "Missense_Mutation Passage", "Truncating Passage", "Truncating Oncogenic",  "AML", "DEL", 'Multi_Hit')

setkey(color_table, 'color_index')
oncogenic_table = data.table(
			     oncogenic_name = c('', 'Likely Neutral', 'Inconclusive', 'Predicted Oncogenic', 'Likely Oncogenic', 'Oncogenic'),
			     oncogenic_value = factor(c('Passage', 'Passage', 'Passage', 'Oncogenic', 'Oncogenic', 'Oncogenic'), levels=c('Oncogenic', 'Passage')))
setkey(oncogenic_table, 'oncogenic_name')
oncogenic_table

maffile
fread(maffile) -> scc.maf
## export bed files for each patient
scc.maf[, {
	pat = .BY
	poses = .SD[!duplicated(paste0(Chromosome, Start_Position, End_Position)), .(Chromosome, Start_Position, End_Position + 1)]
	poses[, name := paste0(Chromosome, ':', Start_Position, ':', V3)]
	poses[, score := 100]
	poses[, strands := '+']
	poses = poses[order(Chromosome, Start_Position),]
	txtfile = paste0('pos_', pat, '.txt')
	fwrite(poses, quote=F, sep="\t", row.names=F, file=txtfile, col.names=F)
	system(paste0('sort -k1,1 -k2,2n ', txtfile, ' > tmp'))
	system(paste0('mv tmp ', txtfile))
	}, by='patient']

bamfiles = fread("bamlist", header=F)
bamfiles[, group := substr(basename(V1), 41, 48)]
bamfiles$group

# convert bam to bed
bamfiles[, bam2bed.jobname := paste(basename(V1))]
bamfiles[, bedfile := paste0('bed/', basename(V1), '.bed')]
bamfiles[, bam2bed.cmd := bsub.head(jobname=bam2bed.jobname, We='1:11', cpu=1, mem=6, cwd=getwd()), by = 1:nrow(bamfiles)]
bamfiles[, bam2bed.cmd := paste0(bam2bed.cmd, ' " ', bb, '/bamToBed -i ', V1, ' > ', bedfile, ' "')]
bamfiles$bam2bed.cmd[1]

exe.jobs(bamfiles$bam2bed.cmd, logger)

bamfiles$bam2bed.jobname[1]

# sort bed file
bamfiles[, sortbed.jobname := paste0('sortbed.', basename(V1))]
bamfiles[, bedsorted := paste0('bed/', basename(V1), '.sorted.bed')]
bamfiles[, sortbed.cmd := bsub.head(jobname=sortbed.jobname, We='1:11', cpu=10, mem=56, cwd=getwd(), postdone=bam2bed.jobname), by=1:nrow(bamfiles)]
bamfiles[, sortbed.cmd := paste0(sortbed.cmd, ' " sort -k1,1 -k2,2n ', bedfile, ' > ', bedsorted, ' "')]
bamfiles$sortbed.cmd[1]

exe.jobs(bamfiles$sortbed.cmd, logger)

## read counts corresponding positions for each patient
bb = '/home/huw/program/bedtools2/bin'
bamfiles[, countfile := paste0('bed/', basename(V1), '.count')]
bamfiles[, {
	pat = .BY
	posfile = paste0('pos_', pat, '.txt')
	tmp = copy(.SD)
	tmp[, jobname := paste0('count.', basename(V1))]
	tmp[, cmd := bsub.head(jobname=jobname, We='1:11', cpu=1, mem=63, cwd=getwd())]
	#tmp[, cmd := paste0(cmd, ' " ', bb, '/bedtools coverage -counts -sorted -g ', genomeFasta, '.fai -a ', posfile, ' -b ', V1, ' > ', countfile,  '"')]
	tmp[, cmd := paste0(cmd, ' " ', bb, '/bedtools coverage -counts -a ', posfile, ' -b ', bedsorted, ' > ', countfile,  '"')]
	exe.jobs(tmp$cmd, logger)
	}, by='group']

bamfiles[, cov.jobname := paste0(basename(V1))]
bamfiles[, outfile := paste0('stats/', basename(V1), '.txt')]
bamfiles
bamfiles[, cmd := bsub.head(jobname=cov.jobname, We='1:11', cpu=1, mem=63, cwd=getwd()), by=1:nrow(bamfiles)]
bamfiles[, cmd := paste0(cmd, ' " ', bb, '/bedtools coverage -a a.bed -b ', V1, ' > ', outfile,  '"')]
#bamfiles[,cmd := paste0(cmd, ' " ', BEDTOOLS, '/bedtools coverage -a pos.txt -b ', V1, ' > ', outfile,  '"')]
bamfiles$cmd[1]

exe.jobs(bamfiles$cmd[1], logger)

bamfiles[, {
	pat = .BY
	posfile = paste0('pos_', pat, '.txt')
	tmp = copy(.SD)
	tmp[, jobname := paste0('count.', basename(V1))]
	tmp[, cmd := bsub.head(jobname=jobname, We='1:11', cpu=1, mem=63, cwd=getwd())]
	#tmp[, cmd := paste0(cmd, ' " ', bb, '/bedtools coverage -counts -sorted -g ', genomeFasta, '.fai -a ', posfile, ' -b ', V1, ' > ', countfile,  '"')]
	tmp[, cmd := paste0(cmd, ' " ', bb, '/bedtools coverage -counts -a ', posfile, ' -b ', bedsorted, ' > ', countfile,  '"')]
	exe.jobs(tmp$cmd, logger)
	}, by='group']

#fread(maffile) -> bak
IMPACT468 <- scan("/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv6/genelist", what = "")
scc.maf[, IMPACT_468 := Hugo_Symbol %in% IMPACT468]

## reassign Variant_Classification for simple version
scc.maf[, id := paste(Variant_Classification, oncogenic, sep=" ")]
setkey(scc.maf, 'Variant_Classification')
scc.maf = scc.maf[simple_table,] 
table(scc.maf$simple_value)
## reassign oncogenic value
setkey(scc.maf, 'oncogenic')
scc.maf = scc.maf[oncogenic_table, ]
table(scc.maf$oncogenic_value)
## assign color
scc.maf[, id := paste(simple_value, oncogenic_value, sep=' ')]
setkey(scc.maf, 'id')
scc.maf = scc.maf[color_table,] 
table(scc.maf$color_value)

scc.maf[, Variant_Classification_old := Variant_Classification]

scc.maf[, Variant_Classification := id]
scc.maf[t_alt_count > 3 & t_alt_count / t_depth >0.05,] -> scc.maf.fil
scc.maf.fil
fwrite(scc.maf.fil, file='scc_maf_fil.txt', sep="\t")

## gene level facets copy number results
cn.table = fread('exonseq/autospy/genelevel_facets_cnv.txt')
cn.table = cn.table[, .(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Type)]
fwrite(cn.table, file='cntable.txt', sep="\t", quote=F)
head(cn.table)
cn.table[Hugo_Symbol == 'BRCA',]

## clinical data
clin = data.table(Tumor_Sample_Barcode = unique(scc.maf.fil[,Tumor_Sample_Barcode]))
clin[, Cell_Type := sub(".*_", "", Tumor_Sample_Barcode)]
fwrite(clin, file='clinical.txt', sep="\t", quote=F)
clin

#read.maf(scc.maf.fil, cnTable='cntable.txt') -> scc.maf.o
#fread('scc_maf_fil.txt') -> a
#read.maf('scc_maf_fil.txt') -> scc.maf.o
read.maf('scc_maf_fil.txt', clinicalData='clinical.txt', vc_nonSyn = vc.key) -> scc.maf.o
as.character(scc.maf.o@data$Tumor_Sample_Barcode) -> tsb
tsb.t1 = tsb[grep("T1", tsb)]
tsb.t2 = tsb[grep("T2", tsb)]
scc.maf.t1.o = subsetMaf(scc.maf.o, tsb=tsb.t1, mafObj = T)
scc.maf.t2.o = subsetMaf(scc.maf.o, tsb=tsb.t2, mafObj = T)

scc.maf.fil[Hugo_Symbol == 'KMT2D', .(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Variant_Classification_old, oncogenic)]
table(scc.maf.o@data$Variant_Classification)

# most amplied known function
mut.2017 = c('TP53', 'RB1', 'RHOB', 'PIK3CA', 'KDM6A', 'TSC1', 'ELF3', 'KMT2D', 'CREBBP', 'CDKN1A', 'EP300', 'ZFP36L1', 'ARID1A', 'STAG2', 'CDKN2A', 'HRAS', 'KRAS', 'FBXW7', 'ERCC2', 'ASXL2', 'RHOA', 'KMT2A', 'FGFR3', 'NFE2L2', 'KMT2C', 'PSIP1', 'KANSL1', 'C3orf70', 'FAT1', 'SPTAN1', 'RXRA', 'ZBTB7B', 'PTEN', 'ATM', 'KLF5', 'PARD3', 'CUL1', 'NRAS', 'SF3B1', 'GNA13', 'RBM10', 'ACTB', 'MBD1', 'CASP8', 'HIST1H3B', 'TAF11', 'ERBB2', 'NUP93', 'SF1', 'ERBB3', 'METTL3', 'SPN', 'MB21D2', 'SSH3', 'USP28', 'ASXL1', 'TMCO4', 'HES1', 'ZNF773')
mut.2014 = c('TP53', 'MLL2', 'ARID1A', 'KDM6A', 'PIK3CA', 'EP300', 'CDKN1A', 'RB1', 'ERCC2', 'FGFR3', 'STAG2', 'ERBB3', 'FBXW7', 'RXRA', 'ELF3', 'NFE2L2', 'TSC1', 'KLF5', 'TXNIP', 'FOXQ1', 'CDKN2A', 'RHOB', 'FOXA1', 'PAIP1', 'BTG2', 'HRAS', 'ZFP36L', 'RHOA', 'CCND3')
amp.2013 = c('AHR', 'BCL2L1', 'CCND1', 'CCNE1', 'E2F3', 'EGFR', 'ERBB2', 'FGFR3', 'GATA3', 'KRAS', 'MDM2', 'MYCL1', 'PPARG', 'PVRL4', 'SOX4', 'TERT', 'YWHAZ', 'ZNF703') 


## plot summary
maffile = "exonseq/Proj_07813_DF/r_001/post/Proj_07813_DF___SOMATIC.vep.filtered.facets.V3.maf"
fread(maffile) -> tmp
tmp[t_alt_count > 3 & t_alt_count / t_depth >0.05,] -> tmp
read.maf(tmp) -> tmp.o
plotmafSummary(tmp.o, rmOutlier=T, addStat = 'median', dashboard = T, file='exonseq/maftools/maftools_summary_scc.pdf', width=10)
subsetMaf(tmp.o, tsb=tsb.t1, mafObj = T) -> tmp.t1.o
plotmafSummary(tmp.t1.o, rmOutlier=T, addStat = 'median', dashboard = T, file='exonseq/maftools/maftools_summary_scc_t1', width=10)
subsetMaf(tmp.o, tsb=tsb.t2, mafObj = T) -> tmp.t2.o
plotmafSummary(tmp.t2.o, rmOutlier=T, addStat = 'median', dashboard = T, file='exonseq/maftools/maftools_summary_scc_t2', width=10)
rm(tmp, tmp.o, tmp.t1.o, tmp.t2.o)

w = 12; h = 8 
mut.2017.sel = c('KMT2D', 'TP53', 'PIK3CA', 'ATM', 'ARID1A', 'EP300', 'FGFR3', 'STAG2', 'KDM6A', 'KMT2C', 
		 'RB1', 'CUL1', 'ELF3', 'ERBB3', 'KMT2A', 'NFE2L2', 'RHOA', 'ACTB', 'CREBBP', 'SPTAN1')

pdf('maftools/oncoplot_mut2017.pdf', width=w, height=h)
oncoplot(maf = scc.maf.o, genes=mut.2017.sel, color = color_vector, showTumorSampleBarcodes=F, top=20, clinicalFeatures='Cell_Type')
dev.off()

unique(scc.maf.t1.o@data$Tumor_Sample_Barcode)
w = 8; h = 8 
pdf('maftools/oncoplot_t1.pdf', width=w, height=h)
oncoplot(maf = scc.maf.t1.o, genes=mut.2017.sel, color = color_vector, showTumorSampleBarcodes=F, top=20)
dev.off()
pdf('maftools/oncoplot_t2.pdf', width=w, height=h)
oncoplot(maf = scc.maf.t2.o, genes=mut.2017.sel, color = color_vector, showTumorSampleBarcodes=F, top=20)
dev.off()
system('rsync -avur maftools/ -e ssh mski1925:/Volumes/LaCie/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/maftools/')

scc.maf[t_alt_count > 3 & t_alt_count / t_depth >0.05,] -> scc.maf.fil
tmp = read.maf(scc.maf.fil)
pdf('exonseq/maftools/lolli_kmt2d.pdf', width=6, height=3)
lollipopPlot(scc.maf.o, gene='KMT2D', cBioPortal=F, colors = color_vector)
dev.off()
sync('exonseq/maftools')

## tcga blca
tcga_sc_id = c('TCGA-BT-A0YX-01A', 'TCGA-BT-A20U-01A', 'TCGA-BT-A2LD-01A', 'TCGA-C4-A0F1-01A', 'TCGA-C4-A0F7-01A', 'TCGA-CU-A0YN-01A',
	       'TCGA-DK-A2I2-01A', 'TCGA-FD-A3B5-01A', 'TCGA-FD-A3N5-01A', 'TCGA-G2-A2ES-01A', 'TCGA-G2-A3IB-01A', 'TCGA-GC-A3I6-01A',
	       'TCGA-GD-A3OS-01A', 'TCGA-BT-A42E', 'TCGA-GU-A42Q', 'TCGA-FD-A43Y', 'TCGA-FD-A5BU', 'TCGA-K4-A4AC', 'TCGA-FD-A5BY', 'TCGA-PQ-A6FI',
	       'TCGA-GU-A766', 'TCGA-CU-A72E', 'TCGA-E7-A7XN', 'TCGA-XF-A8HE', 'TCGA-YC-A89H', 'TCGA-E7-A97P', 'TCGA-XF-A8HH', 'TCGA-ZF-A9RE',
	       'TCGA-ZF-AA4W', 'TCGA-4Z-AA80', 'TCGA-4Z-AA82', 'TCGA-4Z-AA89', 'TCGA-XF-A9SJ', 'TCGA-XF-A9T4', 'TCGA-XF-A9T8', 'TCGA-ZF-AA53',
	       'TCGA-XF-AAME', 'TCGA-XF-AAMH', 'TCGA-XF-AAMT', 'TCGA-XF-AAN2', 'TCGA-XF-AAN5', 'TCGA-ZF-A9RD', 'TCGA-ZF-A9RG', 'TCGA-BT-A20X-01A',
	       'TCGA-FD-A3B4-01A', 'TCGA-FD-A3N6-01A')

tcga_sc_id = sub("-01A", "", tcga_sc_id)
length(unique(tcga_sc_id))

tcga_basal_id = c('TCGA-DK-A3WW', 'TCGA-SY-A9G5', 'TCGA-DK-A2I4', 'TCGA-E7-A7XN', 'TCGA-XF-A9T8', 'TCGA-XF-A9T5', 'TCGA-ZF-AA4V', 'TCGA-BT-A0YX', 'TCGA-DK-AA6S',
		  'TCGA-UY-A9PB', 'TCGA-4Z-AA7W', 'TCGA-4Z-AA84', 'TCGA-BT-A3PJ', 'TCGA-XF-A8HD', 'TCGA-ZF-AA4W', 'TCGA-BT-A20J', 'TCGA-FD-A3SS', 'TCGA-XF-AAN2',
		  'TCGA-XF-A9SY', 'TCGA-DK-A3IU', 'TCGA-G2-A2ES', 'TCGA-XF-A8HE', 'TCGA-XF-A9T3', 'TCGA-DK-A3IN', 'TCGA-DK-AA6L', 'TCGA-FD-A3SN', 'TCGA-FD-A3SO',
		  'TCGA-FD-A5BX', 'TCGA-GC-A3RC', 'TCGA-XF-A9SI', 'TCGA-HQ-A5NE', 'TCGA-K4-A5RJ', 'TCGA-K4-A6FZ', 'TCGA-BT-A3PK', 'TCGA-DK-AA6R', 'TCGA-FD-A3N5',
		  'TCGA-E7-A7DV', 'TCGA-UY-A8OB', 'TCGA-XF-AAN5', 'TCGA-XF-A9T6', 'TCGA-XF-A9SJ', 'TCGA-UY-A78P', 'TCGA-G2-A2EJ', 'TCGA-DK-A2I6', 'TCGA-FD-A3B6',
		  'TCGA-G2-A2EF', 'TCGA-XF-A9SM', 'TCGA-GC-A3YS', 'TCGA-ZF-A9RF', 'TCGA-CF-A1HS', 'TCGA-G2-AA3C', 'TCGA-4Z-AA7Q', 'TCGA-CU-A3KJ', 'TCGA-GC-A6I1',
		  'TCGA-DK-AA74', 'TCGA-E5-A2PC', 'TCGA-4Z-AA86', 'TCGA-4Z-AA81', 'TCGA-ZF-AA58', 'TCGA-DK-A1AB', 'TCGA-UY-A78L', 'TCGA-E7-A97P', 'TCGA-GU-A766',
		  'TCGA-ZF-AA54', 'TCGA-K4-A5RH', 'TCGA-GU-A42Q', 'TCGA-ZF-A9RN', 'TCGA-BT-A42F', 'TCGA-GC-A3I6', 'TCGA-BT-A20O', 'TCGA-ZF-AA53', 'TCGA-C4-A0F1',
		  'TCGA-BT-A42E', 'TCGA-ZF-A9RD', 'TCGA-GU-AATQ', 'TCGA-FD-A3NA', 'TCGA-FT-A61P', 'TCGA-UY-A9PA', 'TCGA-FD-A3B3', 'TCGA-XF-AAN3', 'TCGA-ZF-AA4N',
		  'TCGA-CU-A0YN', 'TCGA-FD-A6TH', 'TCGA-BL-A5ZZ', 'TCGA-GV-A40E', 'TCGA-GC-A3WC', 'TCGA-K4-A4AC', 'TCGA-PQ-A6FI', 'TCGA-FD-A3B5', 'TCGA-FD-A5C1',
		  'TCGA-FD-A6TD', 'TCGA-LC-A66R', 'TCGA-K4-A5RI', 'TCGA-FD-A3SP', 'TCGA-GU-A764', 'TCGA-E7-A3X6', 'TCGA-FT-A3EE', 'TCGA-DK-A6B5', 'TCGA-DK-A2I2',
		  'TCGA-FD-A6TK', 'TCGA-DK-A1AE', 'TCGA-FD-A3N6', 'TCGA-FD-A5BT', 'TCGA-DK-AA6Q', 'TCGA-4Z-AA82', 'TCGA-BT-A2LD', 'TCGA-GD-A3OQ', 'TCGA-GU-A762',
		  'TCGA-CU-A72E', 'TCGA-ZF-AA56', 'TCGA-BL-A13I', 'TCGA-UY-A8OC', 'TCGA-XF-A9SX', 'TCGA-XF-AAMT', 'TCGA-BT-A20V', 'TCGA-G2-A3IB', 'TCGA-FD-A62N',
		  'TCGA-HQ-A5ND', 'TCGA-C4-A0F0', 'TCGA-GV-A3QG', 'TCGA-DK-A3IM', 'TCGA-BT-A20U', 'TCGA-BT-A20X', 'TCGA-PQ-A6FN', 'TCGA-FD-A3B7', 'TCGA-FD-A5BU',
		  'TCGA-FD-A62S', 'TCGA-DK-AA6M', 'TCGA-FD-A3B4', 'TCGA-ZF-AA5H', 'TCGA-FD-A43Y', 'TCGA-ZF-AA5N', 'TCGA-FD-A5BY', 'TCGA-DK-A3WX', 'TCGA-XF-A9T4',
		  'TCGA-XF-AAN4', 'TCGA-XF-AAMW', 'TCGA-FD-A3B8', 'TCGA-FD-A5BS', 'TCGA-XF-AAME', 'TCGA-DK-A3WY', 'TCGA-XF-AAN8');

tcga_sc_basal_id = c(tcga_sc_id, tcga_basal_id)
head(tcga_sc_basal_id)
length(unique(tcga_sc_basal_id)) ## 151

# import tcga blca data
mc3.blca.file = '../../tcga/tcga_mc3_blca_20170830_facets.tsv'
fread(mc3.blca.file) -> mc3.blca
mc3.blca[, vaf := t_alt_count / t_depth]
dim(mc3.blca)
mc3.blca=mc3.blca[t_alt_count > 3 & t_depth > 8 & vaf > 0.05,] 
mc3.blca[, bcr := substr(mc3.blca$Tumor_Sample_Barcode, 1, 12)]
length(unique(mc3.blca$bcr))
length(unique(mc3.blca$Tumor_Sample_Barcode))
dim(mc3.blca)

mc3.blca.sb = mc3.blca[bcr %in% tcga_sc_basal_id,] # squamous and basal  cell
mc3.blca.nn = mc3.blca[!(bcr %in% tcga_sc_basal_id),] # non basal and squamous cell
dim(mc3.blca.sb) # 47247
dim(mc3.blca.nn) # 86527
nrow(mc3.blca.sb) / length(unique(mc3.blca.sb$Tumor_Sample_Barcode)) # 319 per case
nrow(mc3.blca.nn) / length(unique(mc3.blca.nn$Tumor_Sample_Barcode)) # 346 per case

blca.clin = data.table(Tumor_Sample_Barcode = unique(mc3.blca$Tumor_Sample_Barcode))
blca.clin[, cellType := 'Uro']
blca.clin[, bcr := substr(Tumor_Sample_Barcode, 1, 12)]
blca.clin[bcr %in% tcga_basal_id, cellType := 'Basal']
blca.clin[bcr %in% tcga_sc_id, cellType := 'Sq']
blca.clin[bcr %in% intersect(tcga_sc_id, tcga_basal_id), cellType := 'Basal & Sq']
blca.clin
fwrite(blca.clin, file='blca.clin.txt', sep="\t", quote=F)

mc3.blca[, id := paste(Variant_Classification, oncogenic, sep=" ")]
setkey(mc3.blca, 'Variant_Classification')
mc3.blca = mc3.blca[simple_table,] 
mc3.blca
table(mc3.blca$simple_val)
## reassign oncogenic value
setkey(mc3.blca, 'oncogenic')
mc3.blca = mc3.blca[oncogenic_table, ]
table(mc3.blca$oncogenic_value)
## assign color
mc3.blca[, id := paste(simple_val, oncogenic_value, sep=' ')]
setkey(mc3.blca, 'id')
mc3.blca = mc3.blca[color_table,] 
table(mc3.blca$color_value)

mc3.blca[, Variant_Classification_old := Variant_Classification]
mc3.blca[, Variant_Classification := id]

read.maf(mc3.blca, clinicalData=blca.clin, vc_nonSyn = vc.key) -> mc3.blca.o
mc3.blca.sb.o = subsetMaf(mc3.blca.o, tsb = unique(mc3.blca.sb$Tumor_Sample_Barcode), mafObj = T) # squmous or basal
mc3.blca.nn.o = subsetMaf(mc3.blca.o, tsb = unique(mc3.blca.nn$Tumor_Sample_Barcode), mafObj = T) # non squmouse or basal
mc3.blca.sq.o = subsetMaf(mc3.blca.o, tsb = unique(mc3.blca[bcr %in% tcga_sc_id, Tumor_Sample_Barcode]), mafObj = T) # squamous
mc3.blca.sq1.o = subsetMaf(mc3.blca.o, tsb = unique(mc3.blca[bcr %in% setdiff(tcga_sc_id, tcga_basal_id), Tumor_Sample_Barcode]), mafObj = T) # squamous
mc3.blca.ba1.o = subsetMaf(mc3.blca.o, tsb = unique(mc3.blca[bcr %in% setdiff(tcga_basal_id, tcga_sc_id), Tumor_Sample_Barcode]), mafObj = T) # squamous
mc3.blca.ba.o = subsetMaf(mc3.blca.o, tsb = unique(mc3.blca[bcr %in% tcga_basal_id, Tumor_Sample_Barcode]), mafObj = T) # basal
mc3.blca.2sb.o = subsetMaf(mc3.blca.o, tsb = unique(mc3.blca[bcr %in% intersect(tcga_basal_id, tcga_sc_id), Tumor_Sample_Barcode]), mafObj = T) # samples have both basal and squamouse cell types

ll = c('TP53', 'KMT2D', 'ARID1A', 'KDM6A', 'PIK3CA', 'KMT2C', 'RB1', 'STAG2', 'EP300', 'ATM', 'FGFR3', 'SPTAN1', 'ERBB2', 'ELF3', 'CREBBP', 'FAT1', 'KMT2A', 'ERBB3', 'ERCC2', 'CDKN1A')
w = 22; h = 7 
pdf('maftools/oncoplot_mc3_blca_mut2017.pdf', width=w, height=h)
tmp = oncoplot(maf = mc3.blca.o, genes=ll, color = color_vector, showTumorSampleBarcodes=F, top=20, clinicalFeatures='cellType', fontSize=8)
dev.off()

system('rsync -avur maftools/ -e ssh mski1925:/Volumes/LaCie/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/maftools/')

w = 15; h = 7 
ll = c('TP53', 'ARID1A', 'KDM6A', 'KMT2D', 'PIK3CA', 'KMT2C', 'STAG2', 'FGFR3', 'SPTAN1', 'ELF3', 'ATM', 'RB1', 'EP300', 'ERBB2', 'ERBB3', 'CREBBP', 'FAT1', 'CDKN1A', 'TSC1', 'KMT2A')
pdf('maftools/oncoplot_mc3_blca_uro_mut2017.pdf', width=w, height=h)
oncoplot(maf = mc3.blca.nn.o, genes=ll, color = color_vector, showTumorSampleBarcodes=F, top=20, clinicalFeatures='cellType', fontSize=8)
dev.off()
system('rsync -avur maftools/ -e ssh mski1925:/Volumes/LaCie/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/maftools/')

pdf('maftools/oncoplot_mc3_blca_BaSq_mut2017.pdf', width=w, height=h)
oncoplot(maf = mc3.blca.sb.o, genes=mut.2017, color = color_vector, showTumorSampleBarcodes=F, top=20, clinicalFeatures='cellType')
dev.off()

w = 7; h = 7 
ll = c('TP53', 'KMT2D', 'KDM6A', 'KMT2C', 'PIK3CA', 'STAG2', 'FAT1', 'ATM', 'NFE2L2', 'CDKN2A', 'ERBB2', 'RB1', 'RXRA', 'ARID1A', 'FBXW7', 'FGFR3', 'KMT2A', 'CREBBP', 'ERBB3', 'ERCC2')
pdf('maftools/oncoplot_mc3_blca_Sq_mut2017.pdf', width=w, height=h)
oncoplot(maf = mc3.blca.sq.o, genes=ll, color = color_vector, showTumorSampleBarcodes=F, top=20, clinicalFeatures='cellType', fontSize=8, legendFontSize=7, annotationFontSize=7, annotationTitleFontSize=7)
dev.off()
system('rsync -avur maftools/ -e ssh mski1925:/Volumes/LaCie/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/maftools/')


## test
pdf('exonseq/maftools/oncoplot_mc3_test.pdf', width=w, height=h)
tmp = oncoplot(maf = mc3.blca.sq.o, genes=ll, color = color_vector, showTumorSampleBarcodes=T, top=20, clinicalFeatures='cellType', fontSize=8, legendFontSize=7, annotationFontSize=7, annotationTitleFontSize=7)
dev.off()
sync('exonseq', 'maftools')

blca.cc -> blca.cc.v2
colnames(blca.cc.v2) = substr(colnames(blca.cc), 1, 12)
substr(colnames(tmp@ht_list$matrix_3@matrix), 1, 12) -> tmp.cn
ov = intersect(colnames(blca.cc.v2), tmp.cn)
blca.cc.v2['ENSG00000167548', ov]

intersect(colnames(blca.cc.v2),  colnames(tmp@ht_list$matrix_3@matrix))
substr(colnames(blca.cc.v2), 1, 10)
substr(colnames(tmp@ht_list$matrix_3@matrix), 1, 10)
intersect(substr(colnames(blca.cc.v2), 1, 10),substr(colnames(tmp@ht_list$matrix_3@matrix), 1, 10))

blca.cc['ENSG00000167548', colnames(tmp@ht_list$matrix_3@matrix)]

pdf('maftools/oncoplot_mc3_blca_Ba_mut2017.pdf', width=w, height=h)
oncoplot(maf = mc3.blca.ba.o, genes=mut.2017, color = color_vector, showTumorSampleBarcodes=F, top=20, clinicalFeatures='cellType')
dev.off()
pdf('maftools/oncoplot_mc3_blca_2sb_mut2017.pdf', width=w, height=h)
oncoplot(maf = mc3.blca.2sb.o, genes=mut.2017, color = color_vector, showTumorSampleBarcodes=F, top=20, clinicalFeatures='cellType')
dev.off()
pdf('maftools/oncoplot_mc3_bblca_sq1_mut2017.pdf', width=w, height=h) ## only squmaouse
oncoplot(maf = mc3.blca.sq1.o, genes=mut.2017, color = color_vector, showTumorSampleBarcodes=F, top=20, clinicalFeatures='cellType')
dev.off()

ll = c('TP53', 'KMT2D', 'RB1', 'ARID1A', 'PIK3CA', 'EP300', 'KMT2C', 'KDM6A', 'KMT2A', 'CREBBP', 'ERBB2', 'ERCC2', 'ATM', 'ERBB3', 'FAT1', 'ELF3', 'NFE2L2', 'CDKN1A', 'PARD3', 'SPTAN1')
w = 10; h = 7 
pdf('maftools/oncoplot_mc3_bblca_ba1_mut2017.pdf', width=w, height=h) ## only squmaouse
oncoplot(maf = mc3.blca.ba1.o, genes=ll, color = color_vector, showTumorSampleBarcodes=F, top=20, clinicalFeatures='cellType', fontSize=8, legendFontSize=7, annotationFontSize=7, annotationTitleFontSize=7)
dev.off()
system('rsync -avur maftools/ -e ssh mski1925:/Volumes/LaCie/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/maftools/')

w = 10; h = 27 
pdf('maftools/oncoplot_mc3_bblca_ba1_mut2017_moregenes.pdf', width=w, height=h) ## only squmaouse
oncoplot(maf = mc3.blca.ba1.o, genes=mut.2017, color = color_vector, showTumorSampleBarcodes=F, top=20, clinicalFeatures='cellType', fontSize=8, legendFontSize=7, annotationFontSize=7, annotationTitleFontSize=7)
dev.off()
system('rsync -avur maftools/ -e ssh mski1925:/Volumes/LaCie/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/maftools/')

system("find /ifs/res/taylorlab/tcga_wes_facets/blca_2017/ -name '*hisens.out' > aa")
fread("aa", header=F) -> out
out
out$purity = sapply(out$V1, function(x){system(paste("grep Purity ", x, " |sed 's/# Purity = //' "), intern=T)})
out$ploidy = sapply(out$V1, function(x){system(paste("grep Ploidy ", x, " |sed 's/# Ploidy = //' "), intern=T)})
out[, bcr := substr(basename(V1), 1, 12)]
out[, cellType := 'Uro']
out[bcr %in% tcga_sc_id, cellType := 'Sq'] 
out[bcr %in% setdiff(tcga_basal_id, tcga_sc_id), cellType := 'Ba'] 
out[, .(cellType, ploidy), by=cellType]
out[, ploidy := as.numeric(ploidy)]

## no difference
anova(as.data.frame(out[, .(cellType, ploidy)]))
aggregate(out$ploidy, by=list(out$cellType), )
out[ploidy > 2.5,] -> tmp
x1 = aggregate(rep(1,nrow(out)), by=list(out$cellType), sum); x1
x2 = aggregate(rep(1,nrow(tmp)), by=list(tmp$cellType), sum); x2
x2$x/x1$x
out[, cellType := factor(cellType, levels=c('Uro', 'Sq', 'Ba'))]

g = ggplot(out, aes(x =cellType, y = ploidy))
g = g + geom_jitter(width=0.2)
g = g + stat_summary(fun.data = "mean_cl_boot", geom='crossbar', width=.5, mapping=aes(group='cellType', colour = "red"),size=1, show.legend = F)
g = g + ylab('Ploidy') + ggtitle("Anova p < 0.01)") + xlab('')
g = g + theme(plot.title = element_text(face="bold", size=10, hjust=.5))
ggsave(g, file=paste0('maftools/ploidy_tcga_blca.pdf'), width=3, height=4)
system('rsync -avur maftools/ -e ssh mski1925:/Volumes/LaCie/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/maftools/')

out$cellType
summary(aov(out$ploidy~ out$cellType))


## tcga rnaseq comparing uro to squamous cell, basal type cases
load("../../bcg/tcga_blca.RData")
tmp = assay(blca.data)
indi = substr(colnames(tmp), 1, 12)

design = data.frame(tsb = colnames(tmp), celltype = rep('Uro', ncol(tmp)), stringsAsFactors=F)
design$celltype[indi %in% tcga_sc_id] = 'Sq'
design$celltype[indi %in% setdiff(tcga_basal_id, tcga_sc_id)] = 'Ba'
design$celltype = factor(design$celltype, levels=c('Uro', 'Sq', 'Ba'))

ddsmat = DESeqDataSetFromMatrix(countData = tmp,
				colData = design,
				design = ~ celltype);

dds.ds <- estimateSizeFactors(ddsmat);
dds <- DESeq(dds.ds, parallel=T);

table(design$celltype)
res.sq = results(dds, contrast = c('celltype', 'Sq',  'Uro'), nrow=nrow(tmp))
res.sq = res.sq[order(res.sq$padj),]
res.sq$symbol = blca.rows[row.names(res.sq), 'external_gene_name']
head(res.sq)

res.sq[grep("gata3", res.sq$symbol, ignore.case=T),]
aggregate(tmp['ENSG00000107485',], by=list(design$celltype), mean)
aggregate(tmp['ENSG00000229647',], by=list(design$celltype), mean)

res.ba = results(dds, contrast = c('celltype', 'Ba', 'Uro'), nrow=nrow(tmp))
res.ba$symbol = blca.rows[row.names(res.ba), 'external_gene_name']
res.ba = res.ba[order(res.ba$padj),]

res.ba[grep("gata3", res.ba$symbol, ignore.case=T),]
res.ba[grep("foxa1", res.ba$symbol, ignore.case=T),]
res.ba[grep("ppar", res.ba$symbol, ignore.case=T),]
res.ba[grep("KMT2D", res.ba$symbol, ignore.case=T),]

blca.cc = counts(dds, normalized = T)
blca.cc = as.data.frame(tmp.cc)
blca.cc = tmp.cc

p3 =blca.cc[c('ENSG00000167548', 'ENSG00000107485', 'ENSG00000129514', 'ENSG00000132170'),]
p3 = as.data.table(p3)
genes = c('KMT2D', 'GATA3', 'FOXA1', 'PPARG')
p3[, gene := factor(genes, levels = genes)]
p3[1:4,1:2]
p3 = melt(p3)
p3[, bcr := substr(variable, 1, 12)]
p3[, celltype := 'Carcinoma']
p3[bcr %in% tcga_sc_id, celltype := 'Squamous']
p3[bcr %in% setdiff(tcga_basal_id, tcga_sc_id), celltype := 'Basal']
p3[, celltype := factor(celltype, levels=c('Carcinoma', 'Basal', 'Squamous'))]
p3[, log2Reads := log2(value+1)]
p3
gg = ggplot(data = p3, aes(celltype, log2Reads)) + 
	geom_jitter(width=.2) + 
      	stat_summary(fun.data = "mean_cl_boot", geom='crossbar', width=.5, mapping=aes(group='celltype', colour = "red"),size=1, show.legend = F) + 
	xlab('') + ylab('Log2 # Reads') + facet_grid(. ~ gene)
ggsave(gg, file="maftools/three_gene.pdf", width=9, height=4)
system('rsync -avur maftools/ -e ssh mski1925:/Volumes/LaCie/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/maftools/')

## methylation, see the solit/study/tcga/tcga.r folder
## see hiseq/tcga
cc = c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'HGVSp_Short', 'Variant_Classification')
scc.maf[Hugo_Symbol %in% genes, cc, with=F]

##
t1.gs = row.names(scc.maf.t1.o@oncoMatrix)[1:50]
t2.gs = row.names(scc.maf.t2.o@oncoMatrix)[1:50]
intersect(t1.gs, t2.gs)
setdiff(t2.gs, t1.gs)
setdiff(t1.gs, t2.gs)

t1.gs = row.names(scc.maf.t1.o@oncoMatrix)[1:250]
t2.gs = row.names(scc.maf.t2.o@oncoMatrix)[1:250]
t1.gs.i = t1.gs[t1.gs %in% IMPACT468]
t2.gs.i = t2.gs[t2.gs %in% IMPACT468]
intersect(t1.gs, t2.gs)
setdiff(t2.gs, t1.gs)
setdiff(t1.gs, t2.gs)

scc.maf.ip = scc.maf.fil[IMPACT_468 == T, ]
scc.maf.t1.ip = scc.maf.t1[IMPACT_468 == T, ]
scc.maf.t2.ip = scc.maf.t2[IMPACT_468 == T, ]
scc.maf.ip.o =read.maf(scc.maf.ip)
scc.maf.t1.ip.o =read.maf(scc.maf.t1.ip)
scc.maf.t2.ip.o =read.maf(scc.maf.t2.ip)
pdf('oncoplot_impact_all.pdf', width=w, height=h)
oncoplot(maf = scc.maf.ip.o, top=top)
dev.off()
pdf('oncoplot_t1_impact.pdf', width=w, height=h)
oncoplot(maf = scc.maf.t1.ip.o, top=top)
dev.off()
pdf('oncoplot_t2_impact.pdf', width=w, height=h)
oncoplot(maf = scc.maf.t2.ip.o, top=top)
dev.off()


## for autospy
## adjusted:
scc.maf[ccf_Mcopies > 0.3] -> tmp
ds202.t1 = tmp[Tumor_Sample_Barcode == 's_DS_bla_202_T1', Hugo_Symbol]
ds202.t2 = tmp[Tumor_Sample_Barcode == 's_DS_bla_202_T2', Hugo_Symbol]
ds202.m1 = tmp[Tumor_Sample_Barcode == 's_DS_bla_202_M1', Hugo_Symbol]
setdiff(ds202.m1, c(ds202.t1, ds202.t2)) # m1 unique 16
setdiff(ds202.t1, c(ds202.m1, ds202.t2)) # t1 unique 44 
setdiff(ds202.t2, c(ds202.m1, ds202.t1)) # t1 unique 44 
setdiff(intersect(ds202.m1, ds202.t2), ds202.t1) # t1 unique 44 
intersect(ds202.t2, intersect(ds202.m1, ds202.t1)) # t1 unique 44 

ds211.t1 = tmp[Tumor_Sample_Barcode == 's_DS_bla_211_T1', Hugo_Symbol]
ds211.t2 = tmp[Tumor_Sample_Barcode == 's_DS_bla_211_T2', Hugo_Symbol]
ds211.t1
ds211.t2
intersect(ds211.t1, ds211.t2) # t1 unique 44 
setdiff(ds211.t1, ds211.t2) # t1 unique 44 
setdiff(ds211.t2, ds211.t1) # t1 unique 44 




fread(maffile) -> scc.maf
IMPACT468 <- scan("/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv6/genelist", what = "")
scc.maf.fil[, IMPACT_468 := Hugo_Symbol %in% IMPACT468]
scc.maf[, Variant_Bioportal := paste0(Variant_Classification, oncogenic)]

scc.maf.fil[, patient := stringr::str_extract(Tumor_Sample_Barcode, pid)]
scc.maf.fil[, sample := stringr::str_extract(Tumor_Sample_Barcode, sid)]
dim(scc.maf.fil)


scc.maf.fil[, t_var_freq := as.numeric(t_alt_count) / t_depth]
scc.maf.fil[, tm := paste0(Hugo_Symbol, " ", as.numeric(gsub("[^\\d]+", "", HGVSp_Short, perl=T)))] # 
scc.maf.fil[, TAG := paste0('chr', Chromosome, ':', Start_Position, '-', End_Position, ':', Reference_Allele, ':', Tumor_Seq_Allele2)]
scc.maf.fil[, variant := paste0(TAG, '::', Hugo_Symbol, ':', HGVSp_Short)]

sampleID = data.table(`Sample ID` = unique(scc.maf.fil$Tumor_Sample_Barcode))
sampleID[, `Sample Class` := 'Primary']
sampleID[grep("T2", `Sample ID`), `Sample Class` := 'Tumor']
sampleID[grep("M", `Sample ID`), `Sample Class` := 'Tumor']
sampleID[, patient := stringr::str_extract(`Sample ID`, pid)]
sampleID[, sample := stringr::str_extract(`Sample ID`, sid)]
sampleID

o = getwd()
setwd('exonseq/autospy')

IMPACT410 <- scan("/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv5/genelist", what = "")
IMPACT410


## 
source('~/program/autospy/R/process_autopsy_maf.R')
source('~/program/autospy/R/stratton_plot.R')

### refilter
filter_results <- filter_maf_report(scc.maf.fil)
dim(scc.maf.fil)
scc.maf.fil <- filter_results$maf
nrow(tmp)
nrow(scc.maf.fil)

write.table(scc.maf.fil, 'filtered.maf', sep="\t", quote=F, row.names=F)

### overlap plot
make_mutation_overlap_plot(scc.maf.fil, pid = pid, log = F, out = "pre_filter")

## mutation signature
cmd = 'python ~/program/autospy/inst/mutation-signatures/main.py --seed 100 ~/program/autospy/inst/mutation-signatures/Stratton_signatures29.txt filtered.maf signature.out'
system(cmd)
fread("signature.out") -> mut.sig
fwrite(filter_results$report, file='report.txt', sep="\t", quote=F)
ggsave(plot_mutation_signatures(mut.sig, pid, sid, fraction_threshold = 0.5), 
       filename='mutation_singaure_decompoisition_combined.pdf', width=5, height=6)

### patient wise
patients = unique(scc.maf.fil$patient)
tumor_samples = unique(scc.maf.fil$Tumor_Sample_Barcode)

for(pat in patients){
	primary <- sampleID[patient == pat & `Sample Class` == 'Primary', `Sample ID`]
	primary
	pat.maf <- scc.maf.fil[patient== pat,] 
	pat.maf
	print(paste("make_variant_classification_plot for patient ", pat))
	make_variant_classification_plot(pat.maf, out = pat)
	print("make_ccf_plots")
	make_ccf_plots(pat.maf, tumor_samples)
	print("make_stratton_plots")
	make_stratton_plots(pat.maf, tumor_samples, out = pat)
	print("make_binary_tree")
	print("make_mutation_signatures")
	write.table(pat.maf, 'pat.maf', sep="\t", quote=F, row.names=F)
	cmd = paste0('python ~/program/autospy/inst/mutation-signatures/main.py --seed 100 ')
	cmd = paste0(cmd, '~/program/autospy/inst/mutation-signatures/Stratton_signatures29.txt pat.maf signature.out')
	system(cmd)
	fread("signature.out") -> mut.sig
	ggsave(plot_mutation_signatures(mut.sig, pid, sid, fraction_threshold = 0.5), 
	       filename=paste0('mutation_singaure_decompoisition_', pat, '.pdf'), width=5, height=6)
	samples = unique(pat.maf$sample)
	sample_pairs <- combn(samples, 2, simplify = F)
	dir.create('ccf_2d_plots')
	for(sample_pairs in sample_pairs){
		make_ccf_2d(pat.maf,
			    sample_pairs,
			    out = pat,
			    directory = 'ccf_2d_plots')
	}
}

for(pat in patients){
	primary <- sampleID[patient == pat & `Sample Class` == 'Primary', `Sample ID`]
	pat.maf <- scc.maf.fil[patient== pat,] 
	make_binary_tree(pat.maf, primary, hotspots = autospy::hotspots, vertical = TRUE, margins = TRUE, out = pat)
}

### maftools
scc.maf.fil.o = read.maf(scc.maf.fil)
plotmafSummary(scc.maf.fil.o, rmOutlier=T, addStat = 'median', dashboard = T, file='maftools_summary.pdf')

patients = unique(scc.maf.fil$Tumor_Sample_Barcode)
lapply(unique(patients),
       function(tsb){
	       maftools::rainfallPlot(scc.maf.fil.o, tsb = tsb, savePlot = T)
       })

scc.maf.fil$sample

intersect(unique(scc.maf.t1$Hugo_Symbol), unique(scc.maf.t2$Hugo_Symbol))
setdiff(unique(scc.maf.t1$Hugo_Symbol), unique(scc.maf.t2$Hugo_Symbol))
setdiff(unique(scc.maf.t2$Hugo_Symbol), unique(scc.maf.t1$Hugo_Symbol))
unique(scc.maf.t1$Hugo_Symbol)
scc.maf.t1$Hugo_Symbol
scc.maf.t2$Hugo_Symbol

## exonseq/autospy
## cnv
IMPACT341_targets <- suppressWarnings(fread(paste0('grep -v "^@" /ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv6/picard_targets.interval_list')))
setnames(IMPACT341_targets, c("chr", "start", "end", "strand", "name"))
IMPACT341_targets[, symbol := sub("_.*", "", name)]
IMPACT341_targets[, gene.start := min(start), by=symbol][, gene.end := max(end), by = symbol]
IMPACT341_targets = IMPACT341_targets[!duplicated(symbol),]
IMPACT341_targets

source('~/program/facets-suite/geneLevel.R')
fread('cncf_hisens.list', header=F) -> genelevel
genelevel = setNames(genelevel, 'cncfFile')
genelevel[, cncfFile := paste0('../', cncfFile)]
genelevel[, bcr := str_extract(cncfFile, "DS_bla_..._[T|M][0-9]?")]
genelevel[, outfile := paste0('genelevel/', bcr, '.genelevel')]
genelevel[, infofile := sub("cncf.txt", "out", cncfFile)]
genelevel[, rdatafile := sub("out", "Rdata", infofile)] # 

dir.create('genelevel')
head(genelevel)
genelevel[, {
	cnv = get_gene_level_calls(cncf_files = cncfFile)
	cnv[, Tumor_Sample_Barcode := bcr]
	fwrite(cnv, file=outfile, sep="\t", quote=F, row.names=F)
}, by=1:nrow(genelevel)]

cnv.list = lapply(genelevel$outfile, fread)
rbindlist(cnv.list) -> all.cnv
head(all.cnv)
dim(all.cnv)
all.cnv = all.cnv[ !is.na(FACETS_CNA) & FACETS_CNA != 0,]
all.cnv
all.cnv = all.cnv[, Variant_Classification := 'CNV']
all.cnv = all.cnv[, Variant_Type := 'AMP']
all.cnv = all.cnv[FACETS_CNA < 0, Variant_Type := 'DEL']
fwrite(all.cnv, file='genelevel_facets_cnv.txt', row.names=F, sep="\t", quote=F)
colnames(all.cnv)
head(all.cnv)

# facets ploidy and purity inforation
facet.info = data.table(bcr = genelevel$bcr)
facet.info$purity = sapply(genelevel$infofile, function(x){system(paste("grep Purity ", x, " |sed 's/# Purity = //' "), intern=T)})
facet.info$ploidy = sapply(genelevel$infofile, function(x){system(paste("grep Ploidy ", x, " |sed 's/# Ploidy = //' "), intern=T)})
facet.info[, purity := as.numeric(purity)]
facet.info[, ploidy := as.numeric(ploidy)]
facet.info[, sampletype := substr(bcr, 12, 13)]
facet.info[, patient := substr(bcr, 1, 10)]
facet.info[grep('188', bcr),]
facet.info[grep('202', bcr),]

# ploidy
facet.info
tmp = facet.info[grep("M", bcr, invert=T),]
tmp.w = dcast(tmp,  patient ~ sampletype, value.var='ploidy'); tmp.w
mbar = apply(tmp.w[,2:3], 2, median, na.rm=T) 
mbar = as.data.table(mbar)
g = ggplot(tmp, aes(x = sampletype, y = ploidy))
g = g + geom_jitter(width=0.2)
g = g + geom_segment(data=mbar, aes(x = c(.7, 1.7), y = mbar, xend = c(1.3, 2.3), yend = mbar), alpha=.5, lwd=1.5, colour='blue')
g = g + geom_segment(data=tmp.w, aes(x = 1.2, y = T1, xend = 1.8, yend = T2), alpha=0.2, colour='red')
g = g + ylab('Ploidy') + ggtitle("wilcox rank test(paired, p=0.01755)") + xlab('')
g = g + theme(plot.title = element_text(face="bold", size=10, hjust=.5))
ggsave(g, file=paste0('../maftools/ploidy.pdf'), width=3, height=4)
tmp.w
wilcox.test(tmp.w$T1, tmp.w$T2, paired=T)

gg = ggplot(data=facet.info[grep("M", bcr, invert=T),], aes(x = sampletype, y = ploidy)) + geom_point()
ggsave(file="../maftools/ploidy.pdf", width=3, height=4)
facet.info$ploidy

all.cnv[, patient := substr(Tumor_Sample_Barcode, 1, 10)]
all.cnv[, sampletype := substr(Tumor_Sample_Barcode, 12, 13)]

tmp = all.cnv[Hugo_Symbol == gene, .(sampletype,Variant_Type)];tmp
table(tmp)
chisq.test(unlist(tmp[,1]), unlist(tmp[,2]))
tmp
tmp[,1]
tmp[,2]

# cnv  gene specific
all.cnv
genes = c('GATA3', 'FOXA1', 'PPARG')
source('~/program/facets-suite/fPlots_ggplot2.R')
tmp = all.cnv[Hugo_Symbol %in%  genes[3], .(Tumor_Sample_Barcode, Hugo_Symbol, sampletype, Variant_Type, patient, FACETS_CNA)];tmp
tmp.w = dcast(tmp,  patient ~ sampletype, value.var='FACETS_CNA'); tmp.w

for(gene in genes){
	tmp = all.cnv[Hugo_Symbol %in%  gene, .(Tumor_Sample_Barcode, Hugo_Symbol, sampletype, Variant_Type, patient, FACETS_CNA)];tmp
	tmp.w = dcast(tmp,  patient ~ sampletype, value.var='FACETS_CNA'); tmp.w
	for(pat in tmp.w$patient){
		patient = paste0(pat, '_T1')
		load(genelevel[bcr == patient, rdatafile])
		plot.facets.all.output(out, fit, type='png', main=paste0(patient, '| cval: 100'), plotname=paste0('../maftools/', patient, '_', gene), gene.name=gene)
		patient = paste0(pat, '_T2')
		load(genelevel[bcr == patient, rdatafile])
		plot.facets.all.output(out, fit, type='png', main=paste0(patient, '| cval: 100'), plotname=paste0('../maftools/', patient, '_', gene), gene.name=gene)
	}
}



tmp.w[is.na(T1), T1 := 0]
tmp.w[is.na(T2), T2 := 0];tmp.w
tmp.w = tmp.w[T2 - T1 != 0,]; tmp.w
tmp.w[, dd := T2 -T1]
tmp.w
tmp.w$T2 - tmp.w$T1
melt(tmp.w) -> tmp; colnames(tmp) = c('patient', 'sampletype', 'FACETS_CNA'); tmp
mbar = apply(tmp.w[,2:3], 2, median, na.rm=T) 
mbar = as.data.table(mbar); mbar
g = ggplot(tmp, aes(x = sampletype, y = FACETS_CNA))
g = g + geom_jitter(width=0.1)
g = g + geom_segment(data=mbar, aes(x = c(.7, 1.7), y = mbar, xend = c(1.3, 2.3), yend = mbar), alpha=.5, lwd=1.5, colour='blue')
g = g + geom_segment(data=tmp.w, aes(x = 1.2, y = T1, xend = 1.8, yend = T2), alpha=0.2, colour='red')
g = g + ylab('CNV') + ggtitle(gene) + xlab('')
g = g + theme(plot.title = element_text(face="bold", size=10, hjust=.5))
ggsave(g, file=paste0('../maftools/cnv_', gene, '.pdf'), width=3, height=4)
wilcox.test(tmp.w$T1, tmp.w$T2, paired=T)

table(all.cnv[Hugo_Symbol == 'GATA3', .(sampletype, Variant_Type)])
unique(all.cnv$Tumor_Sample_Barcode)

out = export.oncoprint(maf, cnv, res.cc, res.sig)
dim(out)
fwrite(out, file='output_4_oncoprint.txt', sep="\t", quote=F)

##
library(VennDiagram)
scc.maf[grep('202', Tumor_Sample_Barcode),] -> tmp
tmp = tmp[, .(Tumor_Sample_Barcode, Hugo_Symbol)]
tmp$patient = substr(tmp$Tumor_Sample_Barcode, 3, 12)
tmp$celltype = substr(tmp$Tumor_Sample_Barcode, 14, 15)
tmp.m = tmp[grep('M1', Tumor_Sample_Barcode), Hugo_Symbol]
tmp.t1 = tmp[grep('T1', Tumor_Sample_Barcode), Hugo_Symbol]
tmp.t2 = tmp[grep('T2', Tumor_Sample_Barcode), Hugo_Symbol]
tmp.m
pdf("exonseq/autospy/DS_bla_202_ven.pdf", width=8, height=8)
grid.newpage()
draw.triple.venn(area1 = length(unique(tmp.t1)),
		 area2 = length(unique(tmp.t2)),
		 area3 = length(unique(tmp.m)),
		 n12 = length(intersect(tmp.t1, tmp.t2)), 
		 n23 = length(intersect(tmp.t2, tmp.m)), 
		 n13 = length(intersect(tmp.t1, tmp.m)), 
		 n123 = length(intersect(intersect(tmp.t1, tmp.t2), tmp.m)), 
		 category = c("T1", "T2", "M1"), scaled=T)
dev.off()


#merge(all.cnv, scc.maf.fil, by = 'Tumor_Sample_Barcode') -> scc.maf.cnv
## export for oncoprint
## http://www.cbioportal.org/oncoprinter.jsp#
## Sample Gene 
## Alteration Type
## R248W	MISSENSE
## R306fs	TRUNC
## MCN237del	INFRAME
## FUSION	FUSION
## HOMDEL	CNA
## HETLOSS	CNA
## GAIN		CNA
## AMP		CNA
## UP		EXP
## DOWN		EXP
export.oncoprint = function(maf, cnv, res.cc, res.sig){
	maf.sel = maf[, .(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Variant_Type, HGVSp_Short)]
	IMPACT468 <- scan("/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv6/genelist", what = "")
	maf.sel[, IMPACT_468 := Hugo_Symbol %in% IMPACT468]
	maf.sel = maf.sel[IMPACT_468 == T, ]
	#write.table(unique(maf.sel$Variant_Classification), '~/program/fun/oncoprint_table', quote=F, sep="\t")
	fread("~/program/fun/oncoprint_table", header = T) -> tbl
	setkey(tbl, 'Variant_Classification')
	setkey(maf.sel, 'Variant_Classification')
	maf.sel.2 = maf.sel[tbl,]
	maf.sel.2 = maf.sel.2[!is.na(Type),]
	maf.sel.2[, Alteration := sub("p.", "", HGVSp_Short)]
	maf.sel.2$Alteration
	## CNA
	## HOMDEL	CNA
	## HETLOSS	CNA
	## GAIN		CNA
	## AMP		CNA
	cnv = all.cnv
	cnv = cnv[, .(Tumor_Sample_Barcode, Hugo_Symbol, chr, Variant_Classification, Variant_Type, tcn, lcn, FACETS_CNA)]
	cnv = cnv[, Type := 'CNA']
	cnv = cnv[FACETS_CNA == 1, Alteration := 'GAIN']
	cnv = cnv[FACETS_CNA > 1, Alteration := 'AMP']
	cnv = cnv[FACETS_CNA == -1, Alteration := 'HETLOSS']
	cnv = cnv[FACETS_CNA < -1, Alteration := 'HOMDEL']
	cnv = cnv[chr == 'X' & FACETS_CNA == -1, Alteration := 'HOMDEL']
	cnv = cnv[chr == 'X' & FACETS_CNA > 0, Alteration := 'AMP']
	out = rbind(maf.sel.2[, .(Tumor_Sample_Barcode, Hugo_Symbol, Alteration, Type)],cnv[, .(Tumor_Sample_Barcode, Hugo_Symbol, Alteration, Type)])
	res.cc.sig = res.cc[row.names(res.sig),] 
	rs = apply(res.cc.sig, 1, mean)
	res.cc.sig = res.cc.sig[rs > 100,]
	dim(res.cc.sig)
	res.cc.sig.sd = apply(res.cc.sig, 1, sd)
	res.cc.sig.zscore = t(apply(res.cc.sig, 1, scale))
	colnames(res.cc.sig.zscore) = colnames(res.cc.sig)
	res.cc.sig.zscore > 1.96
	res.cc.sig.exp = res.cc.sig.zscore
	res.cc.sig.exp[res.cc.sig.zscore > 1.96] = 'UP'
	res.cc.sig.exp[res.cc.sig.zscore < -1.96] = 'DOWN'
	res.cc.sig.exp = melt(res.cc.sig.exp)
	res.cc.sig.exp = res.cc.sig.exp[res.cc.sig.exp$value %in% c('UP', 'DOWN'), ]
	colnames(res.cc.sig.exp) = c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'Alteration')
	res.cc.sig.exp$Type = 'EXP'
	res.cc.sig.exp$Hugo_Symbol == sub(".*_", "", res.cc.sig.exp$Hugo_Symbol) 
	head(res.cc.sig.exp)
	dim(res.cc.sig.exp)
	rbind(out, res.cc.sig.exp) -> out
	out[, Tumor_Sample_Barcode := sub("s_", "", out$Tumor_Sample_Barcode)]
	out
}

res$ID = row.names(res)
res.dt = as.data.table(res)
res.sel = res.dt[symbol %in% IMPACT468 & !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >1,]
res.sel
res.cc.sel = res.cc[row.names(res.cc) %in% res.sel$ID,]
row.names(res.cc.sel) = sub(".*_", "", row.names(res.cc.sel))
res.cc.sel=res.cc.sel[c(1:5, 8),]
res.cc.sel
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(res.cc.sel, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("res/heatmap_impact_sig_genes.pdf", width=10, height=3.5)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), show_colnames = T,
	      scale='none', show_rownames = T, breaks = breaks, 
	      clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
	      annotation_col =anno.col)
dev.off()


# rnaseq
library(data.table)
source("~/program/configure.R")
cwd = getwd()

## results from David
rna.res = read.csv("rnaseq/scc_patient_batchAdj_outlierExc.csv")
res.inc = read.csv("rnaseq/scc_patient_batchAdj_outlierInc.csv")
res.exc = read.csv("rnaseq/scc_patient_batchAdj_outlierExc.csv")

head(rna.exc)


## old RNA-Seq was from outside, only have bam file
## rnaseq
# convert bam file to oq fastq
bam2fq = '~/program/fun/bam2fa_oq.sh'
bamlist = fread('bamlist', header = T)
bamlist
bamlist[, bases := gsub("-", "_", samplename)]
bamlist

bamlist[, bam2fq.jobname := paste0("bam2fq.", bamlist$bases)]
bamlist[, bam2fq.cmd := paste0(BSUB, " -J ", bam2fq.jobname, " -e ", bam2fq.jobname, ".err -o ", bam2fq.jobname, ".std ")]
bamlist[, bam2fq.cmd := paste0(bam2fq.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 2 "')]
bamlist[, bam2fq.cmd := paste0(bam2fq.cmd, bam2fq, ' ', bamfile, " ", bases, ' "')]
bamlist$bam2fq.cmd[1]

file.exists(bamlist$bamfile)
for(i in 1:nrow(bamlist)){ system(bamlist$bam2fq.cmd[i])}

## old target
old.target = fread("old.target")
old.target = dcast(old.target, samplename~R12, value.var = 'fastq', fun.aggregate = function(x){paste(x, collapse = " ")})
old.target

all(sub("R1", "R2", old.target$R1) == old.target$R2)

old.target[, batch := 'old']
old.target[, bases := samplename]

old.target[, R1cat := paste0(bases, "/", bases, "_R1.fastq.gz")]
old.target[, R2cat := paste0(bases, "/", bases, "_R2.fastq.gz")]
old.target[, mvcmd1 := paste0("mv ", R1, " ", R1cat)]
old.target[, mvcmd2 := paste0("mv ", R2, " ", R2cat)]
old.target$mvcmd1

for(i in 1:nrow(old.target)){ system(old.target$mvcmd1[i])}
for(i in 1:nrow(old.target)){ system(old.target$mvcmd2[i])}

old.target[, R1cat := paste0(bases, "/", bases, "_R1.fastq.gz")]
old.target[, catr1.jobname := paste0("catr1.", old.target$bases)]
old.target[, catr1.cmd := paste0(BSUB, " -J ", catr1.jobname, " -e ", catr1.jobname, ".err -o ", catr1.jobname, ".std ")]
old.target[, catr1.cmd := paste0(catr1.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 2 "')]
old.target[, catr1.cmd := paste0(catr1.cmd, 'cat ', R1, " > ", R1cat, '"')]
old.target$catr1.cmd[1]
#for(i in 1:nrow(old.target)){ system(old.target$catr1.cmd[i])}

old.target[, R2cat := paste0(bases, "/", bases, "_R2.fastq.gz")]
old.target[, catr2.jobname := paste0("catr2.", old.target$bases)]
old.target[, catr2.cmd := paste0(BSUB, " -J ", catr2.jobname, " -e ", catr2.jobname, ".err -o ", catr2.jobname, ".std ")]
old.target[, catr2.cmd := paste0(catr2.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 2 "')]
old.target[, catr2.cmd := paste0(catr2.cmd, 'cat ', R2, " > ", R2cat, '"')]
old.target$catr2.cmd[1]
#for(i in 1:nrow(old.target)){ system(old.target$catr2.cmd[i])}

all(file.exists(old.target$R1cat))
all(file.exists(old.target$R2cat))

old.target[, star.jobname := paste0("star.", old.target$bases)]
old.target[, star.cmd := paste0(BSUB, " -J ", star.jobname, " -e ", star.jobname, ".err -o ", star.jobname, ".std ")]
old.target[, star.cmd := paste0(star.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=129]" -R "rusage[iounits=0]" -n 9 ')]
#old.target[, star.cmd := paste0(star.cmd, '-w "post_done(', catr1.jobname, ')" ')]
#old.target[, star.cmd := paste0(star.cmd, '-w "post_done(', catr2.jobname, ')" ')]
old.target[, star.cmd := paste0(star.cmd, '"', RSEM, "/rsem-calculate-expression -p 8 --paired-end --star --star-path ", STAR_DIR, " --gzipped-read-file --append-names --estimate-rspd ")]
old.target[, star.cmd := paste0(star.cmd, ' ', R1cat, ' ', R2cat, ' ', RSEM_REF, ' ', bases, '/', bases, '"')]
old.target$star.cmd[1]
old.target$bases

for(i in 1:nrow(old.target)){ system(old.target$star.cmd[i])}


## new batch of RNA-Seq
new.target = fread("new.target")
colnames(new.target)
new.target$R12
new.target = dcast(new.target, samplename~R12, value.var = 'fastq', fun.aggregate = function(x){paste(x, collapse = " ")})
new.target$bases = gsub("-", "_", new.target$samplename)
colnames(new.target)
new.target$R1
new.target$R2
new.target[1:5,1:4]

new.target[, mkdir := paste0("mkdir -p ", bases)]
for(i in 1:nrow(new.target)){ system(new.target$mkdir[i])}

new.target[, R1cat := paste0(bases, "/", bases, "_R1.fastq.gz")]
new.target[, R2cat := paste0(bases, "/", bases, "_R2.fastq.gz")]
new.target[, catr1.jobname := paste0("catr1.", new.target$bases)]
new.target[, catr1.cmd := paste0(BSUB, " -J ", catr1.jobname, " -e ", catr1.jobname, ".err -o ", catr1.jobname, ".std ")]
new.target[, catr1.cmd := paste0(catr1.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 2 "')]
new.target[, catr1.cmd := paste0(catr1.cmd, 'cat ', R1, " > ", R1cat, '"')]
new.target$catr1.cmd[1]
#for(i in 1:nrow(new.target)){ system(new.target$catr1.cmd[i])}

new.target[, catr2.jobname := paste0("catr2.", new.target$bases)]
new.target[, catr2.cmd := paste0(BSUB, " -J ", catr2.jobname, " -e ", catr2.jobname, ".err -o ", catr2.jobname, ".std ")]
new.target[, catr2.cmd := paste0(catr2.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 2 "')]
new.target[, catr2.cmd := paste0(catr2.cmd, 'cat ', R2, " > ", R2cat, '"')]
new.target$catr2.cmd[1]
#for(i in 1:nrow(new.target)){ system(new.target$catr2.cmd[i])}

all(file.exists(new.target$R1cat))
all(file.exists(new.target$R2cat))

new.target[, star.jobname := paste0("star.", new.target$bases)]
new.target[, star.cmd := paste0(BSUB, " -J ", star.jobname, " -e ", star.jobname, ".err -o ", star.jobname, ".std ")]
new.target[, star.cmd := paste0(star.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=129]" -R "rusage[iounits=0]" -n 9 ')]
new.target[, star.cmd := paste0(star.cmd, '-w "post_done(', catr1.jobname, ')" ')]
new.target[, star.cmd := paste0(star.cmd, '-w "post_done(', catr2.jobname, ')" ')]
new.target[, star.cmd := paste0(star.cmd, '"', RSEM, "/rsem-calculate-expression -p 8 --paired-end --star --star-path ", STAR_DIR, " --gzipped-read-file --append-names --estimate-rspd ")]
new.target[, star.cmd := paste0(star.cmd, ' ', R1cat, ' ', R2cat, ' ', RSEM_REF, ' ', bases, '/', bases, '"')]
new.target$star.cmd[1]

for(i in 1:nrow(new.target)){ system(new.target$star.cmd[i])}

new.target[, batch := 'new']
old.target[, batch := 'old']

## combine old and new target
ov = intersect(colnames(new.target), colnames(old.target))
target = rbind(new.target[, ov, with=F], old.target[, ov, with=F])
target
target[, rsem := paste0(bases, "/", bases, ".genes.results")]
target[, .(bases)]

target[, R1trim := paste0(bases, '/', bases, '_trim_R1.fastq')]
target[, R2trim := paste0(bases, '/', bases, '_trim_R2.fastq')]
target[, R1trimorphan := paste0(bases, '/', bases, '_trim_orphan_R1.fastq')]
target[, R2trimorphan := paste0(bases, '/', bases, '_trim_orphan_R2.fastq')]

target[, trim.jobname := paste0("trim.", bases)]
target[, trim.cmd := paste0(BSUB, " -J ", trim.jobname, " -e ", trim.jobname, ".err -o ", trim.jobname, ".std ")]
target[, trim.cmd := paste0(trim.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=9]" -R "rusage[iounits=0]" -n 2 ')]
target[, trim.cmd := paste0(trim.cmd, '" /home/huw/program/sickle/sickle pe --pe-file1 ', R1cat, ' --pe-file2 ', R2cat, ' --output-pe1 ', R1trim, ' --output-pe2 ', R2trim)]
target[, trim.cmd := paste0(trim.cmd, ' --output-single ', R1trimorphan, ' -t sanger "')]
target$trim.cmd[1]

write.table(target$trim.cmd[1], file='a')

for(i in 1:nrow(target)){ system(target$trim.cmd[i])}

target[, zipfq.jobname := paste0("zipfq.", bases)]
target[, zipfq.cmd := bsub.head(jobname = zipfq.jobname, mem=3, cpu=1, We="1:11", cwd=cwd, postdone=trim.jobname), by = 1:nrow(target)]
target[, zipfq.cmd := paste0(zipfq.cmd, ' " gzip -f ', R1trim, ' "')]
target$zipfq.cmd[2]

exe.jobs(target$zipfq.cmd[file.exists(target$R1trim)], logger)

file.exists(target$R1trim)
file.exists(target$R2trim)

target[, zipfq2.jobname := paste0("zipfq2.", bases)]
target[, zipfq2.cmd := bsub.head(jobname = zipfq2.jobname, mem=3, cpu=1, We="1:11", cwd=cwd, postdone=trim.jobname), by = 1:nrow(target)]
target[, zipfq2.cmd := paste0(zipfq2.cmd, ' " gzip -f ', R2trim, ' "')]
target$zipfq2.cmd[2]
exe.jobs(target$zipfq2.cmd, logger)

target[, R1trim.gz := paste0(bases, '/', bases, '_trim_R1.fastq.gz')]
target[, R2trim.gz := paste0(bases, '/', bases, '_trim_R2.fastq.gz')]

target[, rsem.jobname := paste0("rsem.", bases)]
target[, rsem.cmd := bsub.head(jobname = rsem.jobname, mem=90, cpu=19, cwd = cwd, We='2:50', postdone=paste(c(zipfq.jobname, zipfq2.jobname), collapse=' ')), by=1:nrow(target)]
target[, rsem.cmd := paste0(rsem.cmd, ' "', RSEM, "/rsem-calculate-expression -p 20 --paired-end --star --star-path ", STAR_DIR, " --gzipped-read-file --append-names --estimate-rspd ")]
target[, rsem.cmd := paste0(rsem.cmd, ' ', R1trim.gz, ' ', R2trim.gz, ' ', RSEM_REF, ' ', bases, '/', bases, '"')]
target$rsem.cmd[1]

exe.jobs(target$rsem.cmd, logger)

## align by star
target[, star.jobname := paste0("star.", bases)]
target[, star.cmd := bsub.head(jobname = star.jobname, mem=90, cpu=20, cwd = cwd, We='2:50'), by=1:nrow(target)]
target[, star.cmd := paste0(star.cmd, ' "', STAR, " --runThreadN 19 --genomeDir ", genomeGRCh38StarDir)]
target[, star.cmd := paste0(star.cmd, ' --sjdbOverhang 49 ')]  
target[, star.cmd := paste0(star.cmd, ' --outSAMtype SAM Unsorted --readFilesCommand zcat  --readFilesIn ', R1trim.gz, ' ', R2trim.gz)]
target[, star.cmd := paste0(star.cmd, ' --outFileNamePrefix ', bases, '/', bases, ' "')]
target$star.cmd[1]
file.size(target$R1trim.gz)
file.exists(target$R1trim.gz)
file.exists(target$R2trim.gz)

exe.jobs(target$star.cmd, logger)

## sort by Coord after star
target[, samfile := paste0(bases, '/', bases, 'Aligned.out.sam')]
target[, sort.jobname := paste0("sort.", bases)]
target[, sortedByCoordBam := paste0(bases, '/', bases, '.Aligned.sortedByCoord.out.bam')]
target[, sort.cmd := bsub.head(jobname = sort.jobname, mem=90, cpu=20, cwd = cwd, We='2:50', postdone=star.jobname), by=1:nrow(target)]
target[, sort.cmd := paste0(sort.cmd, ' "', samtools, ' view -S -b ', samfile, ' | ', samtools, ' sort -@ 20 -m 4G -o ', sortedByCoordBam, ' - "')]
target$sort.cmd[1]

exe.jobs(target$sort.cmd, logger)

## index sorted bam file 
target[, index.jobname := paste0("index.", bases)]
target[, index.cmd := bsub.head(jobname = index.jobname, mem=30, cpu=3, cwd = cwd, We='2:50', postdone=sort.jobname), by=1:nrow(target)]
target[, index.cmd := paste0(index.cmd, ' "', samtools, ' index ', sortedByCoordBam, ' "')]
target$index.cmd[1]

exe.jobs(target$index.cmd, logger)

## use htseq-count to get reads
target[, sortedByCoordBam := paste0(bases, '/', bases, '.Aligned.sortedByCoord.out.bam')]
target[, htseqcount.jobname := paste0("htseqcount.", bases)]
target[, htseqcount.output.file := paste0(bases, '/htseqcount_', bases, '.txt')]
target[, htseqcount.cmd := bsub.head(jobname = htseqcount.jobname, mem=20, cpu=1, cwd = cwd, We='2:50', postdone=index.jobname), by=1:nrow(target)]
target[, htseqcount.cmd := paste0(htseqcount.cmd, ' "', htseqcountexe, " -f bam -r pos -s no --mode=intersection-nonempty -i gene_id ", sortedByCoordBam, ' ', genomeGRCh38GTF, ' > ', htseqcount.output.file, ' "')]
target$htseqcount.cmd[1]

exe.jobs(target$htseqcount.cmd, logger)

options(warnings = -1)

## rnaseq based on rsem
library(DESeq2); 
library(pheatmap); 
library(SomaticSignatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tximport)

target$patientID = sub("_T.", "", target$samplename)
target$patientID = sub("-T.", "", target$patientID)
target$patientID = factor(target$patientID)
target$patientID
target[, rsem := paste0(bases, "/", bases, ".genes.results")]
target[!grep("190", bases),] -> target
file.exists(target$rsem)
target$rsem

rsem = tximport(target$rsem, type = "rsem")
colnames(rsem$counts) = target$bases
rsem.counts = rsem$counts
rsem.length = rsem$length
head(rsem.counts)
head(rsem.length)
colnames(rsem.counts) = target$bases
colnames(rsem.length) = target$bases
rsem.counts[1:10,1:5]
colnames(rsem.counts)
mode(rsem.counts) <- "integer"

## fpkm
head(rsem$counts)
tmp = rsem$counts
head(tmp)
tmp = sweep(tmp,2,colSums(tmp),`/`)
tmp = 10 ^ 9 * (tmp/rsem$length)
rsem$fpkm = tmp
rsem$log2fpkm = log2(tmp+1)
head(rsem$log2fpkm)
rm(tmp)

class(rsem.counts)
#rsem.counts.fil = rsem.counts[rsem.counts.rowsum > 20,] 
rsem.counts.fil = rsem.counts
dim(rsem.counts)
dim(rsem.counts.fil)
rsem.counts[1:10,1:5]

target[, group := sub(".*_", "", bases)]
condition = target$group
condition
patientID = factor(target$patientID)
patientID
design = data.frame(
		    row.names       = colnames(rsem.counts.fil),
		    condition       = condition,
		    patientID 	    = patientID,
		    libType         = rep("PE", ncol(rsem.counts.fil)));
design$patientID
design

ddsmat = DESeqDataSetFromMatrix(countData = rsem.counts.fil,
				colData = design,
				design = ~ patientID + condition);
dds.ds <- estimateSizeFactors(ddsmat);
dds.ds
dds <- DESeq(dds.ds, parallel=T);

res = results(dds)
res.cc = counts(dds, normalized = T)
res.cc = as.data.frame(res.cc)
class(res.cc)

library(future)
 
#rld <- rlog(dds, blind = FALSE)
f <- future({
	rld <- rlog(dds, blind = FALSE)
	rld
}) %plan% multiprocess
rld  <- value(f)


res = res[order(res$pvalue),]
res$symbol = sub(".*_", "", row.names(res))
res.sig = res[!is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 1 & !is.na(res$padj) & res$padj < 0.05,]
dim(res.sig)
head(res.sig)
head(res)
source("~/program/fun/write_rnk.r")
write_rnk(res, file='res/res.rnk')
write.csv(res, file='res/res.csv')

source("~/program/fun/run_gsea.R")
run_gsea('res.rnk')

anno.col = data.frame(row.names = colnames(res.cc), sampletype = target$group, patient = target$patientID)

tmp = res.cc
tmp[tmp < 1] = 1
res.cc.log2 = log2(tmp)
rm(tmp)

lls = seq(0.9, 3, by=0.2)
res.sig2 = res.sig[res.sig$baseMean > 50,]
for(ll in lls){
	res.sig.id = row.names(res.sig2)
	matrix_tmp = res.cc.log2[res.sig.id, ]
	cn = colnames(matrix_tmp)
	matrix_tmp = t(apply(matrix_tmp, 1, scale))
	colnames(matrix_tmp) = cn
	matrix_tmp[matrix_tmp > ll] = ll
	matrix_tmp[matrix_tmp < -ll] = ll 
	breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
	pdf(paste0("res/heatmap_sig_genes_", ll, ".pdf"), width=8, height=8)
	    hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), show_colnames = T,
			  scale='none', show_rownames = F, breaks = breaks, 
			  main='Different genes in paired T2 vs T1', clustering_method = 'complete', 
			  cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
			  annotation_col =anno.col)
	    dev.off()
}

library(ggplot2)
library(data.table)

res.cc.2 = as.data.table(res.cc)
res.cc.2[, symbol := sub(".*_", "", row.names(res.cc))]
res.cc.3 = melt.data.table(res.cc.2, id.vars = c('symbol'), measure.vars=colnames(res.cc))
res.cc.3[, group := substr(variable, 12, 13)]
res.cc.3[, patient := substr(variable, 1, 10)]
head(res.cc.3)
res.cc.3[ value < 1, value := 1]

for (i in 1:10){
for (i in 1:nrow(res.sig)){
	gene = res.sig$symbol[i]
	fdr = round(res.sig$padj[i],3)
	if(is.na(gene)){next}
	tmp = res.cc.3[symbol == gene, ]
	tmp[, log_reads := log2(value)]
	tmp
	tmp.w = dcast(tmp, patient ~ group, value.var='log_reads'); tmp.w
	mbar = apply(tmp.w[,2:3], 2, mean, na.rm=T) 
	mbar = as.data.table(mbar)
	g = ggplot(tmp, aes(x = group, y =log_reads))
	g = g + geom_jitter(width=0.2)
	g = g + geom_segment(data=mbar, aes(x = c(.7, 1.7), y = mbar, xend = c(1.3, 2.3), yend = mbar), alpha=.5, lwd=1.5, colour='blue')
	g = g + geom_segment(data=tmp.w, aes(x = 1.2, y = T1, xend = 1.8, yend = T2 - .2), alpha=0.2, colour='red')
	g = g + ylab('log2 # reads') + ggtitle(gene) + xlab('')
	g = g + theme(plot.title = element_text(face="bold", size=10, hjust=.5))
	ggsave(g, file=paste0('res/dotplot/dots_', gene, '_', res.sig$log2FoldChange[i], '_', res.sig$padj[i], '.pdf'), width=2.5, height=3)
}
mbar
tmp.w

## top 1000 variant genes, pca
pdf("res/rld_pca.pdf", width=8, height=8)
plotPCA(rld)
dev.off()

library(ggfortify)
row.sd = apply(res.cc.log2, 1, sd)
row.sd = row.sd[order(row.sd, decreasing=T)]
names(row.sd[1:1000]) -> sel.id
matrix_tmp = res.cc.log2[row.names(res.cc.log2) %in% sel.id, ]
pca1 = prcomp(t(matrix_tmp), scale. = TRUE)
pca1$x -> pca.x
pca.x[1:10,1:10]
summary(pca1) -> pca.sum
pca.sum.key = pca.sum$importance
rm(pca.sum, pca1)
pca.sum.key
xlab=paste0('PC1 ', round(100*pca.sum.key[2, 1], 0), '%') 
ylab=paste0('PC2 ', round(100*pca.sum.key[2, 2], 0), '%') 
xlab
ylab
t12 = sub(".*_", "", rownames(pca.x)); t12
t12 = factor(t12)
t12
batch = sub("DS_bla_", "", rownames(pca.x))
batch = sub("_.*", "", batch); batch
as.numeric(batch) -> batch
batch
batch[batch < 200] = 1
batch[batch >= 200] = 2
batch = paste0('batch', batch)
batch = factor(batch)
batch

pca.x = as.data.frame(pca.x)
pca.x$t12 = t12
pca.x$batch = batch
pca.x
rm(t12, batch)

pdf("res/ggplot_pca.pdf", width=8, height=8)
ggplot(data=pca.x, aes(x=PC1, y = PC2, shape=batch, fill=t12, color=t12)) +
	geom_hline(yintercept = 0, colour = "gray65") +
	geom_vline(xintercept = 0, colour = "gray65") +
	geom_point(size = 4) +
	ggtitle("PCA analysis of RNA-Seq") +
	xlab(xlab) +
	ylab(ylab)
dev.off()


matri_tmp
pdf("res/topvar1000_pca.pdf", width=8, height=8)
plotPCA(rld)
dev.off()

## rnaseq my count results with 11 pairs
head(rsem.counts)
rsem.counts.11pair = rsem.counts[, gsub("-", "_", incl.id)]
head(rsem.counts.11pair)

design = data.frame(
		    row.names       = colnames(rsem.counts.11pair),
		    condition       = substr(colnames(rsem.counts.11pair), 12, 13),
		    patientID 	    = substr(colnames(rsem.counts.11pair), 1, 10),
		    libType         = rep("PE", ncol(rsem.counts.11pair)));
design

ddsmat = DESeqDataSetFromMatrix(countData = rsem.counts.11pair,
				colData = design,
				design = ~ patientID + condition);

dds.ds <- estimateSizeFactors(ddsmat);
dds.ds
dds <- DESeq(dds.ds, parallel=T);

res.11pair = results(dds)
res.11pair$id = row.names(res.11pair)
res.11pair = as.data.table(res.11pair)
res.11pair = res.11pair[order(pvalue),]
res.11pair[, symbol := sub(".*_", "", id)]
res.11pair[, ensembl := sub("_.*", "", id)]
res.11pair
res.11pair.sig = res.11pair[abs(log2FoldChange) > 1 & padj < 0.05,]
res.11pair.sig
res.11pair.cc = counts(dds, normalized = T)
res.11pair.cc = as.data.frame(res.cc)
class(res.11pair.cc)

## top 1000 variant genes heatmap
row.sd = apply(res.cc.log2, 1, sd)
row.sd = row.sd[order(row.sd, decreasing=T)]
names(row.sd[1:1000]) -> sel.id
matrix_tmp = res.cc.log2[row.names(res.cc.log2) %in% sel.id, ]
matrix_tmp
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
ll = 2
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > ll] = ll
matrix_tmp[matrix_tmp < -ll] = ll 
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
pdf(paste0("res/heatmap_topvar1000_genes_", ll, ".pdf"), width=8, height=8)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), show_colnames = T,
	      scale='none', show_rownames = F, breaks = breaks, 
	      main='top varied 1000 genes', clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
	      annotation_col =anno.col)
dev.off()

## top 100 bladder genes 
gtex_bladder_top100 = read.csv("../gtex_bladder_top100.txt")
gtex_bladder_top100
gs = sub(".*_", "", row.names(res.cc.log2))
matrix_tmp = res.cc.log2[gs %in% gtex_bladder_top100$Gene.Symbol, ]
matrix_tmp
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > ll] = ll
matrix_tmp[matrix_tmp < -ll] = ll 
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
pdf(paste0("res/heatmap_bladder100_genes_", ll, ".pdf"), width=8, height=8)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), show_colnames = T,
	      scale='none', show_rownames = F, breaks = breaks, 
	      main='Different genes in paired T2 vs T1', clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
	      annotation_col =anno.col)
dev.off()

##immune gene target heatmap
# 'C10orf54' -> VSIR
immuno_target_symbol = c('PDCD1', 'CD274', 'PDCD1LG2', 'HAVCR2', 'LAG3', 'BTLA', 'VSIR','CD8A', 'CD4', 'PTPRC', 'MKI67', 'FOXP3', 'CD68', 'STAT1', 'GZMB')
immuno_target_alias = c('PD-1', 'PD-L1', 'PD-L2', 'TIM-3', 'LAG-3', 'BTLA', 'VISTA', 'CD8', 'CD4', 'CD45', 'Ki67', 'Fox-P3', 'CD68', 'pSTAT1', 'granzyme B')
immuno_target = data.frame(row.names = immuno_target_symbol, symbol = immuno_target_symbol, alias = immuno_target_alias, stringsAsFactors = F)
immuno_target
mtx = res.david.cc[res.david.sig[Gene %in% immuno_target$symbol, `entrez`], incl.id]
setkey(res.david, `Entrez ID`)
row.names(mtx) = res.david[row.names(mtx), symbol]
mtx = log2(mtx+1)
mtx
cn = colnames(mtx)
mtx = t(apply(mtx, 1, scale))
colnames(mtx) = cn
ll = 2
mtx[mtx > ll] = ll
mtx[mtx < -ll] = ll 
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
pdf(paste0("rnaseq/res/heatmap_11pair_immuno_targets_sig.pdf"), width=8, height=5)
hp = pheatmap(mtx, color=greenred(length(breaks)+1), show_colnames = T,
	      scale='none', show_rownames = T, breaks = breaks, 
	      main='immune target gene expressions between T1/T2', clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 8,
	      annotation_col =anno.col)
dev.off()
sync('rnaseq', 'res')

##basal luminal classifier genes
load("../../bcg/uc_basal_luminal_marks.RData")

ll = 2
mtx = res.david.cc[res.david.sig[Gene %in% base47, `entrez`], incl.id]
setkey(res.david, `Entrez ID`)
row.names(mtx) = res.david[row.names(mtx), symbol]
mtx = log2(mtx+1)
cn = colnames(mtx)
mtx = t(apply(mtx, 1, scale))
colnames(mtx) = cn
mtx[mtx > ll] = ll
mtx[mtx < -ll] = ll 
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
anno_rows = data.frame(row.names = row.names(mtx), Class = rep("basal", nrow(mtx)), stringsAsFactors = F)
anno_rows$Class[row.names(anno_rows) %in% base47.luminal] = "luminal"
anno_rows
anno_cols = data.frame(row.names=colnames(mtx), CellType = rep("T1", ncol(mtx)) , stringsAsFactors = F)
anno_cols$CellType[grep("T2", row.names(anno_cols))] = "T2"
anno_cols$CellType[grep("_M", row.names(anno_cols))] = "M"
anno_cols
pdf(paste0("rnaseq/res/heatmap_base47_12pair_sig.pdf"), width=6, height=2)
hp = pheatmap(mtx, color=greenred(length(breaks)+1), show_colnames = F,
	      scale='none', show_rownames = T, breaks = breaks, 
	      main='Luminal/basal gene expression (BASE47)', clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = F, fontsize_row = 6,
	      annotation_col =anno_cols, annotation_row = anno_rows)
dev.off()
sync('rnaseq/res')

mtx = res.david.cc[res.david.sig[Gene %in% mcconkey, `entrez`], incl.id]
setkey(res.david, `Entrez ID`)
row.names(mtx) = res.david[row.names(mtx), symbol]
mtx = log2(mtx+1)
cn = colnames(mtx)
mtx = t(apply(mtx, 1, scale))
colnames(mtx) = cn
mtx[mtx > ll] = ll
mtx[mtx < -ll] = ll 
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
anno_rows = data.frame(row.names = row.names(mtx), Class = rep("basal", nrow(mtx)), stringsAsFactors = F)
anno_rows$Class[row.names(anno_rows) %in% base47.luminal] = "luminal"
anno_rows
anno_cols = data.frame(row.names=colnames(mtx), CellType = rep("T1", ncol(mtx)) , stringsAsFactors = F)
anno_cols$CellType[grep("T2", row.names(anno_cols))] = "T2"
anno_cols$CellType[grep("_M", row.names(anno_cols))] = "M"
anno_cols
pdf(paste0("rnaseq/res/heatmap_mdAnderson_12pair_sig.pdf"), width=6, height=2.7)
hp = pheatmap(mtx, color=greenred(length(breaks)+1), show_colnames = F,
	      scale='none', show_rownames = T, breaks = breaks, 
	      main='Luminal/basal gene expression (BASE47)', clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = F, fontsize_row = 6,
	      annotation_col =anno_cols, annotation_row = anno_rows)
dev.off()
sync('rnaseq/res')


## concencus clustering
library(ConsensusClusterPlus)
mtx = res.david.cc[res.david.sig[, `entrez`], incl.id]
dim(mtx)
setkey(res.david, `Entrez ID`)
row.names(mtx) = res.david[row.names(mtx), symbol]
mtx = log2(mtx+1)
mtx
#cn = colnames(mtx)
#mtx = t(apply(mtx, 1, scale))
#colnames(mtx) = cn
pdf(paste0("rnaseq/res/concensusClustering_v2.pdf"), width=8, height=5)
results = ConsensusClusterPlus(mtx, maxK=5,reps=500, writeTable=T,pItem=0.8,pFeature=1, title='consensus',clusterAlg="kmdist",distance="pearson",seed=1262118388.71279)
dev.off()
sync('rnaseq', 'res')

icl = calcICL(results,title=title)
icl

##https://www.biostars.org/p/198789/
Kvec = 2:5
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
	  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
}
# The optimal K
optK = Kvec[which.min(PAC)]
optK

r5 = results[[2]]
r5
r5.order = r5$consensusTree[['order']]
r5.order
r5.mat = r5$consensusMatrix
dim(r5.mat)
r5.mat.2 = r5.mat[r5.order, r5.order]
r5.mat.2
anno.col.tmp2 = anno.col.tmp[r5.order,]
anno.col.tmp2$col = ann.colors[['samp_class']][anno.col.tmp2$samp.class]
pdf('rnaseq/res/consensus_k2.pdf', width=10, height=3)
plot(NULL, xlim=c(0, nrow(anno_col-<- tmp2)), ylim=c(0,4))
nrow(anno <- col <- tmp2) -> nn
rect(1:(nn-1), rep(1, (nn-1)), 2:nn, rep(2, (nn-1)), col=anno <- col <- tmp2$col)
dev.off()

## pamr
library(pamr)
head(res.11pair.sig)
mtx = res.11pair.cc[res.11pair.sig$id,]
tmp = cbind(row.names(mtx), mtx)
write.table(tmp, file='mtx.txt', sep="\t", quote=F)
p.d = pamr.from.excel('mtxv2.txt', 24, sample.labels=T)
p.d.train = pamr.train(p.d)
p.d.train
p.d.res = pamr.cv(p.d.train,p.d)

pdf(file='rnaseq/res/pamr_plotcv.pdf')
pamr.plotcv(p.d.res)
dev.off()
sync('rnaseq', 'res')

## Compute the confusion matrix for a particular model (threshold=4.0) 
pamr.confusion(p.d.res, threshold=4.0)

## Plot the cross-validated class probabilities by class
pdf(file='rnaseq/res/pamr_plotcvprob.pdf')
pamr.plotcvprob(p.d.res, p.d, threshold=4.0)
dev.off()
sync('rnaseq', 'res')

## Plot the class centroids
pdf(file='rnaseq/res/pamr_plotcen.pdf')
pamr.plotcen(p.d.train, p.d, threshold=4.0)
dev.off()
sync('rnaseq', 'res')

## Make a gene plot of the most significant genes
pdf(file='rnaseq/res/pamr_geneplot.pdf')
pamr.geneplot(p.d.train, p.d, threshold=3)
dev.off()
sync('rnaseq', 'res')


# Estimate false discovery rates and plot them
fdr.obj<- pamr.fdr(p.d.train, p.d)
fdr.obj

pdf(file='rnaseq/res/plotfdr.pdf')
pamr.plotfdr(fdr.obj)
dev.off()
sync('rnaseq', 'res')


## List the significant genes
pamr.listgenes(p.d.train, p.d, threshold=3)
pamr.list = as.data.table(pamr.listgenes(p.d.train, p.d, threshold=3)) # 
pamr.list
pamr.list.sig = res.david.sig[symbol %in% pamr.list$id,]
## gene expression of these in tcga blca
mtx = res.david.cc[pamr.list.sig$entrez, incl.id]
pamr.list.sig
row.names(mtx) = pamr.list.sig[row.names(mtx), symbol]
mtx = log2(mtx+1)
mtx
cn = colnames(mtx)
mtx = t(apply(mtx, 1, scale))
colnames(mtx) = cn
mtx
mtx[mtx > 2] = 2
mtx[mtx < -2] = -2
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
pdf("rnaseq/res/heatmap_pamrlist_11pair.pdf", width=6, height=3.5)
hp = pheatmap(mtx, color=greenred(length(breaks)+1), show_colnames = T,
	      scale='none', show_rownames = T, breaks = breaks, 
	      clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
	      annotation_col =anno.col, fontsize_col=6)
dev.off()
sync('rnaseq', 'res')

setkey(res.11pair, 'ensembl')
gs = res.11pair[row.names(blca.cc), symbol]
mtx = blca.cc[gs %in% pamr.list$id,]
res.11pair[row.names(mtx), symbol]
row.names(mtx) = res.11pair[row.names(mtx), symbol]
cn = colnames(mtx)
mtx = t(apply(mtx, 1, scale))
colnames(mtx) = cn
mtx[mtx > 2] = 2
mtx[mtx < -2] = -2
mtx
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
pdf("rnaseq/res/heatmap_pamrlist_tcga_blca.pdf", width=19, height=3.5)
hp = pheatmap(mtx, color=greenred(length(breaks)+1), show_colnames = F,
	      scale='none', show_rownames = T, breaks = breaks, 
	      clustering_method = 'ward', 
	      cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
	      fontsize_col=6)
dev.off()
sync('rnaseq', 'res')

## Try heterogeneity analysis, with class "BL" taken to be the normal group
p.d.train2 <- pamr.train(p.d, hetero="BL")
p.d.res2 <-  pamr.cv(p.d.train2, p.d)

## Look for better threshold scalings
p.d.scale <- pamr.adaptthresh(p.d.train)
p.d.train3 <- pamr.train(p.d, threshold.scale=p.d.scale)
p.d.res3 <-  pamr.cv(p.d.train3, p.d)
p.d.res3


pamr.menu(p.d)

## rnaseq my count results with 11 pairs
condition.new = target$group[target$batch == 'new']
condition.new
patientID.new = factor(target$patientID[target$batch == 'new'])
patientID.new
sample.new = target$bases[target$batch == 'new']
sample.new
design = data.frame(
		    row.names       = colnames(sample.new),
		    condition       = condition.new,
		    patientID 	    = patientID.new,
		    libType         = rep("PE", length(sample.new)));

counts.new = rsem.counts.fil[, target$batch == 'new']
ddsmat = DESeqDataSetFromMatrix(countData = counts.new,
				colData = design,
				design = ~ patientID + condition);
dds.ds <- estimateSizeFactors(ddsmat);
dds.ds
dds.new <- DESeq(dds.ds, parallel=T);

res.new = results(dds.new)
res.cc.new = counts(dds.new, normalized = T)
res.cc.new[res.cc.new < 1] = 1
res.cc.new.log2 = as.data.frame(log2(res.cc.new))

rld.new <- rlog(dds.new, blind = FALSE)
rld.old <- rlog(dds.old, blind = FALSE)


## rnaseq old batch
target[, group := sub(".*_", "", bases)]
condition.old = target$group[target$batch == 'old']
condition.old
patientID.old = factor(target$patientID[target$batch == 'old'])
patientID.old
sample.old = target$bases[target$batch == 'old']
sample.old
design = data.frame(
		    row.names       = colnames(sample.old),
		    condition       = condition.old,
		    patientID 	    = patientID.old,
		    libType         = rep("PE", length(sample.old)));

counts.old = rsem.counts.fil[, target$batch == 'old']
ddsmat = DESeqDataSetFromMatrix(countData = counts.old,
				colData = design,
				design = ~ patientID + condition);
dds.ds <- estimateSizeFactors(ddsmat);
dds.ds
dds.old <- DESeq(dds.ds, parallel=T);

res.old = results(dds.old)
res.cc.old = counts(dds.old, normalized = T)
res.cc.old[res.cc.old < 1] = 1
res.cc.old.log2 = as.data.frame(log2(res.cc.old))
rld <- rlog(dds.old, blind = FALSE)

## rnaseq new batch
condition.new = target$group[target$batch == 'new']
condition.new
patientID.new = factor(target$patientID[target$batch == 'new'])
patientID.new
sample.new = target$bases[target$batch == 'new']
sample.new
design = data.frame(
		    row.names       = colnames(sample.new),
		    condition       = condition.new,
		    patientID 	    = patientID.new,
		    libType         = rep("PE", length(sample.new)));

counts.new = rsem.counts.fil[, target$batch == 'new']
ddsmat = DESeqDataSetFromMatrix(countData = counts.new,
				colData = design,
				design = ~ patientID + condition);
dds.ds <- estimateSizeFactors(ddsmat);
dds.ds
dds.new <- DESeq(dds.ds, parallel=T);

res.new = results(dds.new)
res.cc.new = counts(dds.new, normalized = T)
res.cc.new[res.cc.new < 1] = 1
res.cc.new.log2 = as.data.frame(log2(res.cc.new))

rld.new <- rlog(dds.new, blind = FALSE)
rld.old <- rlog(dds.old, blind = FALSE)

## rnaseq compare between old and new batches
res.old$symbol = sub(".*_", "", row.names(res.old))
res.old.sig = res.old[!is.na(res.old$log2FoldChange) & abs(res.old$log2FoldChange) > 1 & !is.na(res.old$padj) & res.old$padj < 0.05,]
source("~/program/fun/write_rnk.r")
write_rnk(res.old, file='res/res_old.rnk')
write.csv(res.old, file='res/res_old.csv')
source("~/program/fun/run_gsea.R")
run_gsea('res_old.rnk')

res.new$symbol = sub(".*_", "", row.names(res.new))
res.new.sig = res.new[!is.na(res.new$log2FoldChange) & abs(res.new$log2FoldChange) > 1 & !is.na(res.new$padj) & res.new$padj < 0.05,]
source("~/program/fun/write_rnk.r")
write_rnk(res.new, file='res/res_new.rnk')
write.csv(res.new, file='res/res_new.csv')
source("~/program/fun/run_gsea.R")
run_gsea('res_new.rnk')

library(grid)
library(VennDiagram)
pdf("rnaseq/res/ven.pdf", width=8, height=8)
grid.newpage()
draw.triple.venn(area1 = length(unique(res.sig$symbol)),
		 area2 = length(unique(res.new.sig$symbol)),
		 area3 = length(unique(res.old.sig$symbol)),
		 n12 = length(intersect(res.sig$symbol, res.new.sig$symbol)), 
		 n23 = length(intersect(res.new.sig$symbol, res.old.sig$symbol)), 
		 n13 = length(intersect(res.sig$symbol, res.old.sig$symbol)), 
		 n123 = length(intersect(intersect(res.sig$symbol, res.new.sig$symbol), res.old.sig$symbol)), 
		 category = c("Both", "New batch", "Old batch"))
dev.off()

intersect(res.old.sig$symbol, res.new.sig$symbol)
nrow(res.old.sig)
nrow(res.new.sig)

res.old.cc = counts(dds.old, normalized = T)
res.new.cc = counts(dds.new, normalized = T)
res.old.cc[res.old.cc < 1] = 1
res.new.cc[res.new.cc < 1] = 1
res.new.cc.log2 = log2(res.new.cc)
res.old.cc.log2 = log2(res.old.cc)

source('~/program/fun/mypheatmap.R')

res.sig.id = row.names(res.old.sig)
mat = res.old.cc.log2[res.sig.id, ]
mat.sum = apply(mat, 1, sum)
mat = mat[mat.sum > 50,]
anno.col = data.frame(row.names=colnames(mat))
anno.col$cellType = 'T1'
anno.col$cellType[grep("T2", row.names(anno.col))] = 'T2'
mypheatmap(mat, filename='res/heatmap_old_batch_sig.pdf', show.rownames = F, wi=5, hi=6, anno.col = anno.col)

res.sig.id = row.names(res.new.sig)
mat = res.new.cc.log2[res.sig.id, ]
mat.sum = apply(mat, 1, sum)
mat = mat[mat.sum > 50,]
anno.col = data.frame(row.names=colnames(mat))
anno.col$cellType = 'T1'
anno.col$cellType[grep("T2", row.names(anno.col))] = 'T2'
mypheatmap(mat, filename='res/heatmap_new_batch_sig.pdf', show.rownames = F, wi=5, hi=6, anno.col = anno.col)

## immuno convolution analysis based on david
pat.sel=c(200, 202, 204, 206, 207, 208, 210, 211, 185, 186, 189, 193, 194, 195)
pat.sel=incl
pat.sel
conv = fread("rnaseq/tim/SCC.Sample_Info.ImmuneDeconvolutionTest_sheet1.csv", skip=1)
conv
conv.ssgsea = fread("rnaseq/tim/SCC.Sample_Info.ImmuneDeconvolutionTest_sheet2.csv")
conv.ssgsea
conv.cybersort = fread("rnaseq/tim/SCC.Sample_Info.ImmuneDeconvolutionTest_sheet3.csv")
conv.cybersort

conv.ssgsea[, Significance := 'p > 0.05']
conv.ssgsea[pvalue < 0.05, Significance := 'p < 0.05']
gg = ggplot(conv.ssgsea, aes(x=reorder(pop, mean_diff), y=mean_diff, color = Significance, size=-log(pvalue))) + geom_point() +
	theme(axis.text.x = element_text(angle=-90, vjust = 0.5, hjust=0)) + ggtitle('Immune cell infilteration \nby ssGSEA (T2 vs T1) ') + 
	ylab("Mean difference") + xlab("Immune cells") + coord_flip() + 
	geom_hline(yintercept=0, color=adjustcolor(2, .7))
ggsave(gg, file='rnaseq/res/immune_cell_deconv_ssGSEA.pdf', width=4, height=5)
sync('rnaseq/res')

conv.cybersort[, Significance := 'p > 0.05']
conv.cybersort[pvalue < 0.05, Significance := 'p < 0.05']
gg = ggplot(conv.cybersort, aes(x=reorder(pop, mean_diff), y=mean_diff, color = Significance, size=-log(pvalue))) + geom_point() +
	theme(axis.text.x = element_text(angle=-90, vjust = 0.5, hjust=0)) + ggtitle('Immune cell infilteration \nby CyberSort (T2 vs T1) ') + 
	ylab("Mean difference") + xlab("Immune cells") + coord_flip() + 
	geom_hline(yintercept=0, color=adjustcolor(2, .7))
ggsave(gg, file='rnaseq/res/immune_cell_deconv_cybersort.pdf', width=4, height=5)
sync('rnaseq/res')

conv
conv[, t12 := substr(Sample, 12, 13)]
conv = conv[Outliers != 'Y',]
conv[, set.sel, with=F]
colnames(conv)
conv[,.(t12, Mono)]
colnames(conv)
conv$Sample

set.sel = c('Sample', 't12', 'Macro.M0', 'Th2 cells', 'Eosinophils', 'NK CD56bright cells', 'PDL1')
conv.sel = conv[Sample %in% incl.id, c(1,62, 6:61),];conv.sel
save(conv.sel, file="t.RData")
tryCatch(data.table::melt(conv.sel[, 1:57], id.var=c('Sample', 't12'))) -> conv.sel
conv.sel[, pval := {x=.SD[order(Sample),]; wilcox.test(x[t12=='T1',value], x[t12=='T2',value], paired=T)[[3]]}, by=variable]; conv.sel
conv.sel[, patient := substr(Sample, 1, 10)]
conv.sel[, tag := paste0(variable, ' \np=', round(pval,2))]
conv.sel[pval < 0.05,]
gg = ggplot(conv.sel[pval < 0.05,], aes(x=t12, y=value)) +
	geom_jitter(width=.2, color=adjustcolor(1, .8)) +
	stat_summary(fun.data="mean_cl_boot", geom='crossbar', width=.3, color=adjustcolor(4, alpha=.7)) +
	geom_line(aes(group=patient),  color=adjustcolor(2, alpha=.3)) +
	xlab('') + ylab("Deconvolution Score") + facet_wrap(~ tag, scale="free") +
	theme(strip.text.x = element_text(size = 7))
ggsave(gg, file=paste0("rnaseq/res/deconv_selected_12pair.pdf"), width=4.5, height=5)
sync('rnaseq/res')

conv

## rnaseq based on david Kuo, Fengshen
exl = c(192, 194, 188, 191, 187, 190, 211, 206); exl
exl = c(192, 188, 191, 187, 190, 211, 206); exl

#incl = setdiff(substr(colnames(res.david.cc), 8,10), exl)
#incl = c(200, 202, 204, 207, 208, 210, 185, 186, 189, 193, 194, 195)
incl = c(200, 202, 204, 207, 208, 210, 185, 186, 189, 193, 195, 194)
incl.id = c(paste0('DS-bla-', incl, '-T1'),paste0('DS-bla-', incl, '-T2')); incl.id
incl.id
length(incl.id)

setwd('../rnaseq')
fread("scc_patient_batchAdj_outlierExc.csv") -> res.david
res.david[, `Entrez ID` := as.character(`Entrez ID`)]
res.david

load('rnaseq/SCC.Sample.2Batch.hg19KnownGene.RData')
dim(raw_cnt_mtx)
raw.mtx = raw_cnt_mtx[, incl.id]
dim(raw.mtx)
colnames(raw.mtx)
design = data.frame(row.names       = colnames(raw.mtx),
		    libType         = rep("PE", ncol(raw.mtx)));
design$celltype = substr(row.names(design), 12, 13)
design$patient = substr(row.names(design), 1, 10)
design
ddsmat = DESeqDataSetFromMatrix(countData =raw.mtx,
				colData = design,
				design = ~ patient + celltype);

dds.ds <- estimateSizeFactors(ddsmat);
dds <- DESeq(dds.ds, parallel=T);
dds

# my results
res.david.deseq = results(dds)
res.david.deseq$entrez = row.names(res.david.deseq)
as.data.table(res.david.deseq) -> res.david.deseq
# add symbol
setkey(res.david, `Entrez ID`)
setkey(res.david.deseq, `entrez`)
res.david
res.david.deseq = merge(res.david.deseq, res.david[, .(`Entrez ID`, Gene)], by.x = 'entrez', by.y = 'Entrez ID', all.x = T)
res.david.deseq = res.david.deseq[order(pvalue),]
res.david.deseq
# overlapp
res.david.deseq[padj < 0.05 & abs(log2FoldChange) >1,]
res.david[padj < 0.05 & abs(log2FoldChange) >1,]
intersect(res.david.deseq[padj < 0.05 & abs(log2FoldChange) >1, Gene],res.david[padj < 0.05 & abs(log2FoldChange) >1,Gene])

# this is for the heatmap and dotplot
dds
res.david.cc = counts(dds, normalized = T)
colnames(res.david.cc)
dim(res.david.cc)
res.david.sig=res.david[padj < 0.05 & abs(log2FoldChange) >1,]

setnames(res.david.sig, 'Entrez ID', 'entrez')

## res.david heatmap
incl
p5 = incl.id
#p5 = sample(unique(substr(incl.id, 8, 10)), 5)
#p5= c(185, 189, 195, 186, 193)
#p5 = c(paste0('DS-bla-', p5, '-T1'), paste0('DS-bla-', p5, '-T2')); p5
res.david.deseq.sig=res.david.deseq[padj < 0.01 & abs(log2FoldChange) >1,]
res.david.cc.sel = log2(res.david.cc[res.david.deseq.sig$entrez, p5]+1)
colnames(res.david.cc.sel)
substr(colnames(res.david.cc.sel), 8, 10)
anno.col = data.frame(row.names=colnames(res.david.cc.sel), celltypeid= substr(colnames(res.david.cc.sel), 12, 13))
anno.col$cellType[anno.col$celltypeid == 'T1'] = 'Urothelium'
anno.col$cellType[anno.col$celltypeid == 'T2'] = 'Basal'
anno.col$celltypeid = NULL
anno.col
matrix_tmp = res.david.cc.sel
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
ll=2
matrix_tmp[matrix_tmp > 2] = ll
matrix_tmp[matrix_tmp < -2] = -ll
matrix_tmp
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
pdf("rnaseq/res/heatmap_12pair_v3.pdf", width=3.5, height=5)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), show_colnames = F,
	      scale='none', show_rownames = F, breaks = breaks, 
	      main='Different genes in paired T2 vs T1', clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = TRUE, fontsize=6,fontsize_row = 6,
	      annotation_col =anno.col)
dev.off()
sync('rnaseq/res')

## res.david gsea
setnames(res.david, 'Gene', 'symbol')
source("~/program/fun/write_rnk.r")
write_rnk(res.david, file='res/res.rnk')

source("~/program/fun/run_gsea.R")
run_gsea('res.rnk', nperm=10000, mem=60, cpu=10)
sync('rnaseq', 'gsea')

system('find rnaseq/gsea -name "gsea_report*.xls" > gsea.sum.files')
gsea.sum.files = fread('gsea.sum.files', header=F); gsea.sum.files
gsea.sum.files = gsea.sum.files[grep('c2', V1), ][grep('cgp', V1, invert=T),]; gsea.sum.files
gsea.sum = lapply(gsea.sum.files$V1, fread)
rbindlist(gsea.sum) -> gsea.sum
dim(gsea.sum)
gsea.sum[order(`FDR q-val`),][grep('BIOCAR', NAME),][!duplicated(NAME),] [`RANK AT MAX` < 3500, ][,.(NAME, `FDR q-val`, `RANK AT MAX`)]
tmp = gsea.sum[order(`FDR q-val`),][!duplicated(NAME),][`RANK AT MAX` < 3500, ][`FDR q-val` < 0.05,][1:20,.(NAME, `FDR q-val`, `RANK AT MAX`, NES)];tmp
tmp[, NAME := factor(NAME, levels=NAME)]
gg = ggplot(tmp, aes(x=NAME, y=-log(`FDR q-val` + 1e-6)*sign(NES))) + geom_bar(stat='identity') + theme(axis.text.x=element_text(angle=90, size=6, vjust=1, hjust=1)) + xlab('') + coord_flip() + ylab('-log(FDR)')
ggsave(gg, file='rnaseq/res/gsea_sum.pdf', width=17, height=5)
sync('rnaseq/res')
tmp = fl.dt[geneset == 'Hallmark' & qval < 0.001 & abs(NES) > 2,]
tmp[qval < 0.0001, qval := 0.0001]
tmp
gg = ggplot(tmp, aes(x=reorder(NAME, NES), y=NES)) + geom_jitter(width=0.2, height=.05) +
	theme(axis.text.x = element_text(angle=-90, vjust = 0.5, hjust=0)) +
	ylab("NES") + xlab("Gene set names") + coord_flip()
ggsave(gg, file='rnaseq/res/gsea_hallmark.pdf', width=9, height=3.5)
sync('rnaseq/res')


## res.david maplot
source('~/program/fun/gg_maplot.r')
gg_maplot(res.david[order(padj),], labx = 'log2 FC (T2/T1)', genes = c('GATA3', 'FOXA1', 'PPARG'), filename='rnaseq/res/maplot.pdf')
sync('rnaseq/res')

fanconi = c('FANCA', 'FANCB', 'FANCC', 'BRCA2', 'FANCD2', 'FANCE', 
	    'FANCF', 'FANCG', 'FANCI', 'BRIP1', 'FANCL', 'FANCM', 
	    'PALB2', 'RAD51C', 'SLX4', 'ERCC4', 'RAD51', 'BRCA1', 
	    'UBE2T', 'XRCC2', 'MAD2L2', 'RFWD3') 
gg_maplot(res.david[order(padj),], labx = 'log2 FC (T2/T1)', genes = fanconi, topn=0, filename='rnaseq/res/maplot_fanconi.pdf')
sync('rnaseq/res')

gf = res.david[grep('^FOX', symbol), symbol]
gg_maplot(res.david[order(padj),], labx = 'log2 FC (T2/T1)', genes = gf, topn=0, filename='rnaseq/res/maplot_FOX.pdf')
gf = res.david[grep('^GJ', symbol), symbol]
gg_maplot(res.david[order(padj),], labx = 'log2 FC (T2/T1)', genes = gf, topn=0, filename='rnaseq/res/maplot_GJ.pdf')
gf = res.david[grep('ARHG', symbol), symbol]
gg_maplot(res.david[order(padj),], labx = 'log2 FC (T2/T1)', genes = gf, topn=0, filename='rnaseq/res/maplot_ARHG.pdf')
gf = res.david[grep('APOBEC', symbol), symbol]
gg_maplot(res.david[order(padj),], labx = 'log2 FC (T2/T1)', genes = gf, topn=0, filename='rnaseq/res/maplot_APOBEC.pdf')
sync('rnaseq/res')

## res.david dotplot
res.david

dim(res.david.cc)
head(res.david.cc)
gene.sel = c('GATA3', 'PPARG', 'FOXA1', 'KMT2D')
gene.sel.sig = res.david[symbol %in% gene.sel,]; 
gene.sel.cc = res.david.cc[gene.sel.sig$`Entrez ID`, incl.id]; gene.sel.cc
gene.sel.cc = as.data.table(gene.sel.cc)
gene.sel.cc[, symbol := gene.sel.sig$symbol]; gene.sel.cc
gene.sel.cc[, tag := paste0(gene.sel, ' p=', round(gene.sel.sig$padj, 2))]
gene.sel.cc = as.data.table(melt(gene.sel.cc, id.va=c('tag','symbol')))
gene.sel.cc[, value := log2(value)]
gene.sel.cc[, t12 := substr(variable, 12, 13)]
gene.sel.cc[, patient := substr(variable, 1, 10)]
gene.sel.cc[, group := patient]
gene.sel.cc
g = ggplot(gene.sel.cc, aes(x = t12, y = value, group = patient))
g = g + geom_jitter(width=0.2) + stat_summary(fun.data = 'mean_cl_boot', geom='crossbar', color=adjustcolor(2, alpha=0.8))
g = g + ylab('log2 # reads') + xlab('') + geom_line(position='dodge', color=adjustcolor(2, .2))
g = g + theme(plot.title = element_text(face="bold", size=10, hjust=.5)) + ggtitle("")
g = g + facet_wrap(~tag)
ggsave(g, file='rnaseq/res/dotplot_4genes_12pair.pdf', width=3.5, height=5)
sync('rnaseq/res')

for (i in 1:nrow(res.david.sig)){
	gene = res.david.sig$symbol[i]
	entrez = res.david.sig$entrez[i]
	fc = round(res.david.sig$log2FoldChange[i], 2)
	fdr = round(res.david.sig$padj[i],3)
	if(is.na(entrez)){next}
	if(is.na(gene)){next}
	tmp = res.david.cc[entrez, ]
	tmp = as.data.table(tmp, keep.rownames=T)
	tmp[, log_reads := log2(tmp)]
	tmp[, group := substr(rn, 12, 13)]
	tmp[, patient := substr(rn, 1, 10)]
	g = ggplot(tmp, aes(x = group, y = log_reads, group = group))
	g = g + geom_jitter(width=0.2) + stat_summary(fun.data = 'mean_cl_boot', geom='crossbar', color=adjustcolor(2, alpha=0.6), width=.3)
	g = g + ylab('log2 # reads') + ggtitle(paste0(gene, "\np=", fdr)) + xlab('') + geom_line(aes(group=patient), color=adjustcolor(4, .2))
	g = g + theme(plot.title = element_text(face="bold", size=10, hjust=.5)) 
	ggsave(g, file=paste0('rnaseq/res/dotplot12pair/dotplot_', gene, '_fdr', fdr, '_log2fc', fc, '.pdf'), width=2, height=2.5)
}
sync('rnaseq', 'res/dotplot/')


##TCGA dataset 
##CNV for HLA genes
##Haplodeficiencies GATA3/FOXA1/PPARG
res5[res5$symbol == 'FOXA1',]

options(jupyter.plot_mimetypes = c('image/jpeg'))
options(repr.plot.width=10, repr.plot.height=4)
par(mfrow=c(1,3))
genesymbols = c("GATA3", "FOXA1", "CTSE")
for(genesymbol in genesymbols){
    var1 = as.numeric(t1[genesymbol,])
    var2 = as.numeric(t2[genesymbol,])
    pval = t2_t1_sum_df[genesymbol, "pvalues"]
    title = paste0(genesymbol, " n=", round(pval, 2))
    plot_dot(var1, var2, genesymbol, title)
}

samples = unique(mut$Tumor_Sample_Barcode)

tmp = mut[mut$Tumor_Sample_Barcode == samples[1], ]
tmp_r = tmp$t_alt_count / tmp$t_depth
length(tmp_r)

shapiro.test(tmp_r)

mut_r = log2(mut$t_alt_count / mut$t_depth)
mut_r = mut$t_alt_count / mut$t_depth
plot(density(tmp_r))

## columns need to check sometimes
cn = c("Hugo_Symbol", "Chromosome", "Variant_Classification", "Variant_Type", "Reference_Allele", 
       "Tumor_Seq_Allele1", 'Tumor_Seq_Allele2','Tumor_Sample_Barcode', 't_depth', 't_ref_count', 't_alt_count',
       'n_depth', 'n_ref_count', 'n_alt_count')

mut = read.table("data_mutations_extended.txt", sep="\t", header=T, quote="", stringsAsFactors = F)
colnames(mut)


length(unique(mut$Tumor_Sample_Barcode))

snp = mut[mut$Variant_Type == "SNP", ]
sample_n = length(unique(snp$Tumor_Sample_Barcode))
#average coverage based on snp
# tumor cell coverage
mean(snp$t_depth)
mean(snp$n_depth)
mean(snp$t_depth + snp$n_depth)

unique(snp$Variant_Classification)


# nonsynonymous mutations
nonsyno_tag = c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Translation_Start_Site')
nonsyno = snp[snp$Variant_Classification %in% nonsyno_tag, ]
dim(nonsyno)
nonsyno_df = as.data.frame(table(nonsyno$Variant_Classification))
# synonymous mutations
syno = snp[snp$Variant_Classification == "Silent", ]
dim(syno)

nonsyno_df

tmp = table(nonsyno$Tumor_Sample_Barcode)
mean(tmp)
sd(tmp)

tmp = table(syno$Tumor_Sample_Barcode)
mean(tmp)
sd(tmp)

tmp = table(mut$Tumor_Sample_Barcode)
pat = substr(names(tmp), 8, 10)
tmp = aggregate(as.numeric(tmp), by=list(pat), mean)
tmp
mean(tmp$x)
sd(tmp$x)
range(tmp$x)
mean(tmp$x)*1000000/98754738

# see program/runback/mysql_from_UCSC.sh
# hg19 exon seq len is 98754738 
# see this for exon size discussion: http://seqanswers.com/forums/showthread.php?t=5298

3281/160

unique(mut$Variant_Classification)

unique(mut$Variant_Type)

colnames(mut)

tmp$Tumor_Sample_Barcode

require(lubridate)
days = 365*2
date = seq(as.Date("2000-01-01"), length = days, by = "day")
year = year(date)
month = month(date)
x1 = cumsum(rnorm(days, 0.05)) 
x2 = cumsum(rnorm(days, 0.05))
df1 = data.frame(date, year, month, x1, x2)
df_melt <- melt(df1, id = c("date", "year", "month"))
head(df_melt)

library(reshape2)

id = 187
tmp = mut[grep(id, mut$Tumor_Sample_Barcode), ]
tmp$VAF = tmp$t_alt_count/tmp$t_depth
#tmp = tmp[, c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode', 'VAF')]
range(tmp$VAF)


blca_nature_mut = c('TP53', 'KMT2D', 'ARID1A', 'KDM6A', 'PIK3CA', 'EP300', 'CDKN1A', 'RB1', 
                    'ERCC2', 'FGFR3', 'STAG2', 'ERBB3', 'FBXW7', 'RXRA', 'ELF3', 'NFE2L2', 
                    'TSC1', 'KLF5', 'TXNIP', 'FOXQ1', 'CDKN2A', 'RHOB', 'FOXA1', 'PAIP1', 'BTG2', 'HRAS', 'ZFP36L1', 'RHOA', 'CCND3')

blca_nature_cna = c('CDKN2A', 'E2F3', 'SOX4', 'CCND1', 'RB1', 'EGFR', 'PPARG', 'PVRL4', 'YWHAZ', 
                    'MDM2', 'ERBB2', 'CREBBP', 'NCOR1', 'YAP1', 'CCNE1', 'MYC', 'ZNF703', 'FGFR3', 'PTEN', 'MYCL', 'BCL2L1')
blca_nature=c(blca_nature_mut, blca_nature_cna)

onco_color = c(Mutation = "#26A818", Missense = "#26A818", Nonsense = "black", Splicing = "#ffaa00", 
               Frameshift = "#A05E35" , Promoter = "#2986E2", InFrame = "#F26529", Present = "darkorchid2", NotPresent = "#DCD9D3", 
               NotTested = "darkgrey", del = "red", LOH = "#D17878", homodel = "brown4", 
               CNLOH =  "deepskyblue", Amplification = "#EA2E49", Deletion = "#174D9D", Yes = "#155B6B", No = "#12C8F9", 
               Unknown = "azure1", Fusion =  "#D38C1F", Pathogenic="white")

mut_color = c(  #'Intron' = 'white',
                'Frame_Shift_Ins' = 'Frameshift',
                #'3\'UTR' = 'white',
                'Missense_Mutation' = 'Missense',
                #'3\'Flank' = 'white',
                #'Targeted_Region' = 'white',
                'In_Frame_Ins' = 'InFrame',
                #'5\'Flank' = 'Promoter',
                'Nonsense_Mutation' = 'Nonsense',
                #'Silent' = 'Mutation',
                'In_Frame_Del' = 'InFrame',
                'Frame_Shift_Del' = 'Frameshift',
                #'5\'UTR' = 'Promoter',
                'Splice_Site' = 'Splicing',
                #'RNA' = 'Mutation',
                #'IGR' = 'Mutation',
                'Nonstop_Mutation' = 'Missense_mutation',
                'Translation_Start_Site' = 'Nonsense'
)


focus_mut = c('Frame_Shift_Ins', 'Frame_Shift_Del',  'In_Frame_Ins', 'In_Frame_Del', 
'Missense_Mutation','Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site')


ids = unique(substr(unique(mut$Tumor_Sample_Barcode), 8, 10))
ids

options(repr.plot.width=9, repr.plot.height=9)
par(mfrow=c(2,2))
for(id in ids){
    tmp = mut[grep(id, mut$Tumor_Sample_Barcode), ]
    tmp$VAF = tmp$t_alt_count/tmp$t_depth
    #tmp = tmp[, c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode', 'VAF')] 
    tmp$col = onco_color[mut_color[tmp$Variant_Classification]]
    tmp = tmp[tmp$Variant_Classification %in% focus_mut, ]

    tmp2 = dcast(tmp, Hugo_Symbol+Chromosome+Start_Position+Variant_Classification+Variant_Type+col ~ Tumor_Sample_Barcode, sum, value.var = 'VAF')

    id1 = paste0('DS_bla_', id, '_T1')
    id2 = paste0('DS_bla_', id, '_T2')
    plot(tmp2[,id1], tmp2[,id2], cex=.7, xlab=paste0(id1, " VAF"), ylab=paste0(id2, " VAF"), col = tmp2$col, xlim=c(0,1), ylim=c(0,1))
    tmp3 = tmp2[tmp2$Hugo_Symbol %in% blca_nature & tmp2[,id1] != 0 & tmp2[,id2] !=0,]
    if(nrow(tmp3) > 0){
        text(tmp3[, id1], tmp3[, id2], tmp3$Hugo_Symbol, adj=c(-0.1,0), cex = .7, col = tmp3$col)
        points(rep(0.75, nrow(tmp3)), seq(from=1, by = -0.08, length.out = nrow(tmp3)), col=tmp3$col)
        text(rep(0.75, nrow(tmp3)), seq(from=1, by = -0.08, length.out = nrow(tmp3)), tmp3$Hugo_Symbol, adj=c(-0.2,.5))
        }
    }
plot.new()
legend('topright', legend=names(mut_color), col = onco_color[mut_color], cex=.7, pch=1)


options(repr.plot.width=9, repr.plot.height=9)
par(mfrow=c(2,2))
for(id in ids){
    tmp = mut[grep(id, mut$Tumor_Sample_Barcode), ]
    tmp$VAF = tmp$t_alt_count/tmp$t_depth
    #tmp = tmp[, c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode', 'VAF')] 
  
mut_sigs_file <- autospy::mutation.signatures(filtered_maf_file)
mut_sigs <- suppressWarnings(fread(mut_sigs_file))

## write report header
cat(paste0("#",
	   capture.output(print(as.data.frame(tdt(as.data.table(filter_results$column_definitions))),
				row.names = FALSE, right = FALSE, useSource = F)),
	   "\n"),
    file = "report.txt")

## ...and rest of report
write.tab(filter_results$report,
	  "report.txt",
	  append = TRUE)


  ###### PLOTS #####
  print("make_mutation_overlap_plot")
  make_mutation_overlap_plot(maf, pid = a$pid, log = F, out = "post_filter")
  make_mutation_signatures(mut_sigs,
                           pid = a$pid,
                           sid = a$sid,
                           out = "combined")

  ### optional maftools stuff
  if (requireNamespace("maftools", quietly = TRUE)) {
    library(maftools)
    maftools_maf = maftools::read.maf(maf = filtered_maf_file, 
                                      removeSilent = F, 
                                      useAll = T)

    maftools::plotmafSummary(maf = maftools_maf, 
                   rmOutlier = T, 
                   addStat = 'median', 
                   dashboard = TRUE,
                   file = "maftools_summary.pdf")
    
    ### make rainfall plots (to display kataegis)
    lapply(unique(maftools_maf@data$Tumor_Sample_Barcode), 
           function(tsb){
             maftools::rainfallPlot(maftools_maf, tsb = tsb, savePlot = T)
           })
  }
    



tumor_samples = unique(scc.maf.fil$Tumor_Sample_Barcode)
process_autopsy_maf()
setwd(o)

unique(scc.maf$Tumor_Sample_Barcode)

# rnaseq
library(data.table)
library(BiocParallel)
register(MulticoreParam(10))

source('configure_dir.R')
source("~/program/configure_rnaseq.R")
cwd = getwd()

## old RNA-Seq was from outside, only have bam file
# convert bam file to oq fastq
bam2fq = '~/program/fun/bam2fa_oq.sh'
bamlist = fread('bamlist', header = T)
bamlist
bamlist[, bases := gsub("-", "_", samplename)]
bamlist

bamlist[, bam2fq.jobname := paste0("bam2fq.", bamlist$bases)]
bamlist[, bam2fq.cmd := paste0(BSUB, " -J ", bam2fq.jobname, " -e ", bam2fq.jobname, ".err -o ", bam2fq.jobname, ".std ")]
bamlist[, bam2fq.cmd := paste0(bam2fq.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 2 "')]
bamlist[, bam2fq.cmd := paste0(bam2fq.cmd, bam2fq, ' ', bamfile, " ", bases, ' "')]
bamlist$bam2fq.cmd[1]

file.exists(bamlist$bamfile)
for(i in 1:nrow(bamlist)){ system(bamlist$bam2fq.cmd[i])}

## old target
old.target = fread("old.target")
old.target = dcast(old.target, samplename~R12, value.var = 'fastq', fun.aggregate = function(x){paste(x, collapse = " ")})
old.target

all(sub("R1", "R2", old.target$R1) == old.target$R2)

old.target[, batch := 'old']
old.target[, bases := samplename]

old.target[, R1cat := paste0(bases, "/", bases, "_R1.fastq.gz")]
old.target[, R2cat := paste0(bases, "/", bases, "_R2.fastq.gz")]
old.target[, mvcmd1 := paste0("mv ", R1, " ", R1cat)]
old.target[, mvcmd2 := paste0("mv ", R2, " ", R2cat)]
old.target$mvcmd1
for(i in 1:nrow(old.target)){ system(old.target$mvcmd1[i])}
for(i in 1:nrow(old.target)){ system(old.target$mvcmd2[i])}

old.target[, R1cat := paste0(bases, "/", bases, "_R1.fastq.gz")]
old.target[, catr1.jobname := paste0("catr1.", old.target$bases)]
old.target[, catr1.cmd := paste0(BSUB, " -J ", catr1.jobname, " -e ", catr1.jobname, ".err -o ", catr1.jobname, ".std ")]
old.target[, catr1.cmd := paste0(catr1.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 2 "')]
old.target[, catr1.cmd := paste0(catr1.cmd, 'cat ', R1, " > ", R1cat, '"')]
old.target$catr1.cmd[1]
#for(i in 1:nrow(old.target)){ system(old.target$catr1.cmd[i])}

old.target[, R2cat := paste0(bases, "/", bases, "_R2.fastq.gz")]
old.target[, catr2.jobname := paste0("catr2.", old.target$bases)]
old.target[, catr2.cmd := paste0(BSUB, " -J ", catr2.jobname, " -e ", catr2.jobname, ".err -o ", catr2.jobname, ".std ")]
old.target[, catr2.cmd := paste0(catr2.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 2 "')]
old.target[, catr2.cmd := paste0(catr2.cmd, 'cat ', R2, " > ", R2cat, '"')]
old.target$catr2.cmd[1]
#for(i in 1:nrow(old.target)){ system(old.target$catr2.cmd[i])}

all(file.exists(old.target$R1cat))
all(file.exists(old.target$R2cat))

old.target[, star.jobname := paste0("star.", old.target$bases)]
old.target[, star.cmd := paste0(BSUB, " -J ", star.jobname, " -e ", star.jobname, ".err -o ", star.jobname, ".std ")]
old.target[, star.cmd := paste0(star.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=129]" -R "rusage[iounits=0]" -n 9 ')]
#old.target[, star.cmd := paste0(star.cmd, '-w "post_done(', catr1.jobname, ')" ')]
#old.target[, star.cmd := paste0(star.cmd, '-w "post_done(', catr2.jobname, ')" ')]
old.target[, star.cmd := paste0(star.cmd, '"', RSEM, "/rsem-calculate-expression -p 8 --paired-end --star --star-path ", STAR_DIR, " --gzipped-read-file --append-names --estimate-rspd ")]
old.target[, star.cmd := paste0(star.cmd, ' ', R1cat, ' ', R2cat, ' ', RSEM_REF, ' ', bases, '/', bases, '"')]
old.target$star.cmd[1]
old.target$bases

for(i in 1:nrow(old.target)){ system(old.target$star.cmd[i])}

## new batch of RNA-Seq
new.target = fread("new.target")
colnames(new.target)
new.target$R12
new.target = dcast(new.target, samplename~R12, value.var = 'fastq', fun.aggregate = function(x){paste(x, collapse = " ")})
new.target$bases = gsub("-", "_", new.target$samplename)
colnames(new.target)
new.target$R1
new.target$R2
new.target

new.target[, mkdir := paste0("mkdir -p ", bases)]
for(i in 1:nrow(new.target)){ system(new.target$mkdir[i])}

new.target[, R1cat := paste0(bases, "/", bases, "_R1.fastq.gz")]
new.target[, R2cat := paste0(bases, "/", bases, "_R2.fastq.gz")]
new.target[, catr1.jobname := paste0("catr1.", new.target$bases)]
new.target[, catr1.cmd := paste0(BSUB, " -J ", catr1.jobname, " -e ", catr1.jobname, ".err -o ", catr1.jobname, ".std ")]
new.target[, catr1.cmd := paste0(catr1.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 2 "')]
new.target[, catr1.cmd := paste0(catr1.cmd, 'cat ', R1, " > ", R1cat, '"')]
new.target$catr1.cmd[1]
#for(i in 1:nrow(new.target)){ system(new.target$catr1.cmd[i])}

new.target[, catr2.jobname := paste0("catr2.", new.target$bases)]
new.target[, catr2.cmd := paste0(BSUB, " -J ", catr2.jobname, " -e ", catr2.jobname, ".err -o ", catr2.jobname, ".std ")]
new.target[, catr2.cmd := paste0(catr2.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 2 "')]
new.target[, catr2.cmd := paste0(catr2.cmd, 'cat ', R2, " > ", R2cat, '"')]
new.target$catr2.cmd[1]
#for(i in 1:nrow(new.target)){ system(new.target$catr2.cmd[i])}

all(file.exists(new.target$R1cat))
all(file.exists(new.target$R2cat))

new.target[, star.jobname := paste0("star.", new.target$bases)]
new.target[, star.cmd := paste0(BSUB, " -J ", star.jobname, " -e ", star.jobname, ".err -o ", star.jobname, ".std ")]
new.target[, star.cmd := paste0(star.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=129]" -R "rusage[iounits=0]" -n 9 ')]
new.target[, star.cmd := paste0(star.cmd, '-w "post_done(', catr1.jobname, ')" ')]
new.target[, star.cmd := paste0(star.cmd, '-w "post_done(', catr2.jobname, ')" ')]
new.target[, star.cmd := paste0(star.cmd, '"', RSEM, "/rsem-calculate-expression -p 8 --paired-end --star --star-path ", STAR_DIR, " --gzipped-read-file --append-names --estimate-rspd ")]
new.target[, star.cmd := paste0(star.cmd, ' ', R1cat, ' ', R2cat, ' ', RSEM_REF, ' ', bases, '/', bases, '"')]
new.target$star.cmd[1]

for(i in 1:nrow(new.target)){ system(new.target$star.cmd[i])}

new.target[, batch := 'new']

## combine old and new target
ov = intersect(colnames(new.target), colnames(old.target))
target = rbind(new.target[, ov, with=F], old.target[, ov, with=F])
target[, rsem := paste0(bases, "/", bases, ".genes.results")]

target[, R1trim := paste0(bases, '/', bases, '_trim_R1.fastq.gz')]
target[, R2trim := paste0(bases, '/', bases, '_trim_R2.fastq.gz')]
target[, R1trimorphan := paste0(bases, '/', bases, '_trim_orphan_R1.fastq.gz')]
target[, R2trimorphan := paste0(bases, '/', bases, '_trim_orphan_R2.fastq.gz')]

target[, trim.jobname := paste0("trim.", bases)]
target[, trim.cmd := paste0(BSUB, " -J ", trim.jobname, " -e ", trim.jobname, ".err -o ", trim.jobname, ".std ")]
target[, trim.cmd := paste0(trim.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=9]" -R "rusage[iounits=0]" -n 2 ')]
target[, trim.cmd := paste0(trim.cmd, '" ', java, ' -jar /home/huw/local/bin/trimmomatic.jar PE ', R1cat, ' ', R2cat, ' ', R1trim, ' ', R1trimorphan, ' ', R2trim, ' ', R2trimorphan, ' ILLUMINACLIP:/home/huw/program/trimmomatic/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 "')]
target$trim.cmd[1]

for(i in 1:nrow(target)){ system(target$trim.cmd[i])}

target[, rsem.jobname := paste0("rsem.", bases)]
target[, rsem.cmd := paste0(BSUB, " -J ", rsem.jobname, " -e ", rsem.jobname, ".err -o ", rsem.jobname, ".std ")]
target[, rsem.cmd := paste0(rsem.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=89]" -R "rusage[iounits=0]" -n 19 ')]
target[, rsem.cmd := paste0(rsem.cmd, '-w "post_done(', trim.jobname, ')" ')]
target[, rsem.cmd := paste0(rsem.cmd, '"', RSEM, "/rsem-calculate-expression -p 18 --paired-end --rsem --rsem-path ", STAR_DIR, " --gzipped-read-file --append-names --estimate-rspd ")]
target[, rsem.cmd := paste0(rsem.cmd, ' ', R1trim, ' ', R2trim, ' ', RSEM_REF, ' ', bases, '/', bases, '"')]
target$rsem.cmd[1]

for(i in 1:nrow(target)){ system(target$rsem.cmd[i])}


options(warnings = -1)

library(DESeq2); 
library(pheatmap); 
library(SomaticSignatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tximport)

target$patientID = sub("_T.", "", target$samplename)
target$patientID = sub("-T.", "", target$patientID)
target$patientID = factor(target$patientID)
target$patientID
target[, rsem := paste0(bases, "/", bases, ".genes.results")]
file.exists(target$rsem)
target$rsem
rsem = tximport(target$rsem, type = "rsem")
rsem.counts = rsem$counts
colnames(rsem.counts) = target$bases
rsem.counts[1:10,1:5]
colnames(rsem.counts)
mode(rsem.counts) <- "integer"

rsem.counts.rowsum = rowSums(rsem.counts)
class(rsem.counts)
rsem.counts.fil = rsem.counts[rsem.counts.rowsum > 20,] 
dim(rsem.counts)
dim(rsem.counts.fil)
rsem.counts[1:10,1:5]

condition = target$group
patientID = target$patientID
design = data.frame(
		    row.names       = colnames(rsem.counts.fil),
		    condition       = condition,
		    patientID 	    = patientID,
		    libType         = rep("PE", ncol(rsem.counts.fil)));
ddsmat = DESeqDataSetFromMatrix(countData = rsem.counts.fil,
				colData = design,
				design = ~ patientID + condition);
dds.ds <- estimateSizeFactors(ddsmat);
dds <- DESeq(dds.ds, parallel=T);

res = results(dds)
res.cc = counts(dds, normalized = T)

res = res[order(res$pvalue),]
res$symbol = sub(".*_", "", row.names(res))
res.sig = res[!is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 1 & !is.na(res$padj) & res$padj < 0.05,]
dim(res.sig)
head(res.sig)
head(res)
source("~/program/fun/write_rnk.r")
write_rnk(res, file='res/res.rnk')
write.csv(res, file='res/res.csv')

source("~/program/fun/run_gsea.R")
run_gsea('res.rnk')
sync('rnaseq', 'gsea')


library(pheatmap)
library(gplots)
breaks = c( seq(from = -2, to = 2, length.out = 30))
breaks = unique(breaks)

anno.col = data.frame(row.names = colnames(res.cc), sampletype = target$group, patient = target$patientID)

res.sig.id = row.names(res.sig)
matrix_tmp = res.cc[res.sig.id, ]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("res/heatmap_sig_genes.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), show_colnames = T,
	      scale='none', show_rownames = F, breaks = breaks, 
	      main='Different genes in paired T2 vs T1', clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
	      annotation_col =anno.col)
dev.off()

library(ggplot2)
res.cc.2 = as.data.table(res.cc)
res.cc.2[, symbol := sub(".*_", "", row.names(res.cc))]
res.cc.3 = melt.data.table(res.cc.2, id.vars = c('symbol'), measure.vars=colnames(res.cc))
res.cc.3[, group := substr(variable, 12, 13)]
res.cc.3[, patient := substr(variable, 1, 10)]
head(res.cc.3)
res.cc.3[ value < 1, value := 1]

for (i in 1:nrow(res.sig)){
	gene = res.sig$symbol[i]
	if(is.na(gene)){next}
	tmp = res.cc.3[symbol == gene, ]
	tmp[, log_reads := log2(value)]
	tmp
	tmp.w = dcast(tmp, patient ~ group, value.var='log_reads'); tmp.w
	g = ggplot(tmp, aes(x = group, y = log2(value)))
	g = g + geom_jitter(width=0.2)
	g = g + geom_segment(data=tmp.w, aes(x = 1.2, y = T1, xend = 1.8, yend = T2 - .2), alpha=0.2, colour='red')
	g = g + ylab('log2 # reads') + ggtitle(gene) + xlab('')
	g = g + theme(plot.title = element_text(face="bold", size=10, hjust=.5))
	ggsave(g, file=paste0('res/dots_', gene, '_', res.sig$log2FoldChange[i], '_', res.sig$padj[i], '.pdf'), width=2.5, height=3)
}

## start here
tmp = apply(rnaseq_filtered_fil, 1, sd)


tmp_fpkm = rnaseq_filtered_fil[tmp > 5, ]
dim(tmp_fpkm) 
options(repr.plot.width=6, repr.plot.height=6)
pheatmap(tmp_fpkm, scale="row", col=greenred(15), show_rownames = F, 
         clustering_method = "complete", cluster_rows = T, cluster_cols = T)

t2_t1_sub_rm_na = t1[1,, drop=F]
sss = t1[1,, drop=F]
log2fc = c()
pvalues = c()
t1_mean = c()
t2_mean = c()
nn = c()
for(i in 1:nrow(t2)){
    ss = t1[i, ] > -9.9 & t2[i, ] > -9.9
    ss_len = length(ss[ss])
    tmp = as.numeric(t2[i, ] - t1[i, ])
    tmp[!ss] = NA
    t2_t1_sub_rm_na = rbind(t2_t1_sub_rm_na, tmp)
    sss = rbind(sss, ss)
    pvalue = NA
    if(length(tmp[ss]) >= 2){
        ttest = t.test(as.numeric(tmp[ss]))
        pvalue = ttest$p.value 
    }
    pvalues = c(pvalues, pvalue)
    log2fc = c(log2fc, mean(tmp[ss]))
    t1_mean = c(t1_mean, mean(as.numeric(t1[i, ss])))
    t2_mean = c(t2_mean, mean(as.numeric(t2[i, ss])))
    nn = c(nn, length(ss[ss]))
}

t2_t1_sub_rm_na = t2_t1_sub_rm_na[-1,]
row.names(t2_t1_sub_rm_na) = row.names(t2)

sss = sss[-1,]
row.names(sss) = row.names(t2)

head(t2_t1_sub_rm_na)
head(sss)

t2_t1_sum_df_rm_na = data.frame(row.names = row.names(t2), 
                                t1_mean = t1_mean, 
                                t2_mean = t2_mean,
                                log2FC = log2fc, 
                                pvalues = pvalues,
                                nn = nn)
t2_t1_sum_df_rm_na$padj = p.adjust(t2_t1_sum_df_rm_na$pvalues, method = "BH", )


t2_t1_sum_df_rm_na$genesymbol = row.names(t2_t1_sum_df_rm_na)
t2_t1_sum_df_rm_na = t2_t1_sum_df_rm_na[order(t2_t1_sum_df_rm_na$pvalues), ]

head(t2_t1_sum_df_rm_na)

rm(t1_mean, t2_mean, log2fc, pvalues)

head(t2_t1_sum_df_rm_na)


options(repr.plot.width=12, repr.plot.height=6)
par(mfrow=c(2,4))
for(i in 1:8){
    genesymbol = t2_t1_sum_df_rm_na$genesymbol[i]
    ss = as.logical(sss[genesymbol, ])
    var1 = as.numeric(t1[genesymbol,ss])
    var2 = as.numeric(t2[genesymbol,ss])
    pval = t2_t1_sum_df_rm_na[genesymbol, "pvalues"]
    t2_t1_sum_df_rm_na[genesymbol, "nn"] -> nn
    title = paste0(genesymbol, "(p=", round(pval, 2), " n=", nn, ")")
    plot_dot(var1, var2, genesymbol, title)    
}

tmp = t2_t1_sum_df_rm_na
tmp[tmp$genesymbol != '' & !is.na(tmp$log2FC) & !is.na(tmp$pvalues), ] -> tmp;
tmp$gsea = sign(tmp$log2FC) * (- log10(tmp$pvalues))
tmp = tmp[order(tmp$gsea, tmp$pvalues),]
head(tmp)


head(t2_t1_sum_df)

source("/Volumes/LaCie/huw/program/fun.r")
source("/home/huw/program/autospy/R/stratton_plot.R")
tmp = t2_t1_sum_df[t2_t1_sum_df$pvalues < 0.05, c("log2FC", "pvalues", "symbol")]
colnames(tmp) = c("log2FoldChange", "pvalue", "symbol")
write_rnk(tmp, filename = "scc_t2_t1_sig.rnk")
## one selene:
## ~/program/runback/run_gsea_selene.sh scc_t2_t1_rm_na.rnk 1>2 > log

head(t2_t1_sum_df_rm_na)

t2_t1_sum_df_rm_na[!is.na(t2_t1_sum_df_rm_na$padj) & 
                   t2_t1_sum_df_rm_na$pvalues < 0.05 & abs(t2_t1_sum_df_rm_na$log2FC) > 1, ] -> tmp
dim(tmp)
dim(tmp[tmp$log2FC > 0 & tmp$pvalues < 0.05, ])
dim(tmp[tmp$log2FC < 0 & tmp$pvalues < 0.05, ])


head(tmp[tmp$log2FC > 0 & tmp$pvalues < 0.05, ])

head(tmp[tmp$log2FC < 0 & tmp$pvalues < 0.05, ])

dim(t2_t1_sum_df_rm_na[t2_t1_sum_df_rm_na$padj < 0.05, ])
intersect(t2_t1_sum_df_rm_na[t2_t1_sum_df_rm_na$padj < 0.05 & t2_t1_sum_df_rm_na$log2FC > 0, "genesymbol"], 
t2_t1_sum_df_rm_na[t2_t1_sum_df_rm_na$padj < 0.05 & t2_t1_sum_df_rm_na$log2FC < 0, "genesymbol"])

library(WriteXLS)
WriteXLS(t2_t1_sum_df_rm_na, ExcelFileName = "scc_t2_t1_rm_na.xls")

t2_t1_sub = t2 - t1

t2_t1_sum_df = data.frame(row.names = row.names(t2), mean = rowMeans(t2_t1_sub))

t2_t1_sum_df$pvalues = apply(t2_t1_sub, 1, function(x) {
    ttest = t.test(as.numeric(x))
    ttest$p.value})
t2_t1_sum_df$t1_mean = rowMeans(t1)
t2_t1_sum_df$t1_sd = apply(t1, 1, sd)
t2_t1_sum_df$t2_mean = rowMeans(t2)
t2_t1_sum_df$t2_sd = apply(t2, 1, sd)
t2_t1_sum_df$symbol = row.names(t2_t1_sum_df)

t2_t1_sum_df$padj = p.adjust(t2_t1_sum_df$pvalues, method = "BH", )
t2_t1_sum_df$log2FC = t2_t1_sum_df$t2_mean - t2_t1_sum_df$t1_mean


t2_t1_sum_df = t2_t1_sum_df[order(t2_t1_sum_df$pvalues, decreasing = F),]


head(t2_t1_sum_df)

col = t2_t1_sum_df$pvalues
col = rep(1, nrow(t2_t1_sum_df))
col[t2_t1_sum_df$pvalues<0.01 & abs(t2_t1_sum_df$log2FC) > 1] = 2
plot(t2_t1_sum_df$log2FC, -log(t2_t1_sum_df$pvalues), xlab="log2 Fold Change", ylab="-log(p values)", col=col, cex=.6)

head(t2_t1_sum_df)

library(WriteXLS)

t2_t1_sum_df$gsea = sign(t2_t1_sum_df$log2FC) * -log10(t2_t1_sum_df$pvalues)
t2_t1_sum_df$symbol = row.names(t2_t1_sum_df)
head(t2_t1_sum_df)
WriteXLS(t2_t1_sum_df, ExcelFileName = "scc_t2_t1.xls")
write.table(t2_t1_sum_df[!is.na(t2_t1_sum_df$gsea) ,c("symbol", "gsea")], file="scc_t2_t1.rnk", sep="\t", quote=F, col.names =F, row.names = F)

t1_ordered = t1[t2_t1_sum_df$symbol, ]
t2_ordered = t2[t2_t1_sum_df$symbol, ]


options(repr.plot.width=12, repr.plot.height=6)
par(mfrow=c(2,4))
for(i in 1:8){
    genesymbol = t2_t1_sum_df$symbol[i]
    
    var1 = as.numeric(t1[genesymbol,])
    var2 = as.numeric(t2[genesymbol,])
    
    pval = t2_t1_sum_df[genesymbol, "pvalues"]
    title = paste0(genesymbol, "(p=", round(pval, 2), ")")

    plot_dot(var1, var2, genesymbol, title)
    
}

t1['PDCD1LG2',]  #PD-L2 in UC
t2['PDCD1LG2',]  #PD-L2 in SCC

immuno_target_symbol = c('PDCD1', 'CD274', 'PDCD1LG2', 'HAVCR2', 'LAG3', 'BTLA', 'C10orf54', 'CD8A', 'CD4', 'PTPRC', 'MKI67', 'FOXP3', 'CD68', 'STAT1', 'GZMB')
immuno_target_alias = c('PD-1', 'PD-L1', 'PD-L2', 'TIM-3', 'LAG-3', 'BTLA', 'VISTA', 'CD8', 'CD4', 'CD45', 'Ki67', 'Fox-P3', 'CD68', 'pSTAT1', 'granzyme B')
immuno_target = data.frame(row.names = immuno_target_symbol, symbol = immuno_target_symbol, alias = immuno_target_alias, stringsAsFactors = F)
immuno_target

tmp = t2_t1_sum_df_rm_na[t2_t1_sum_df_rm_na$genesymbol %in% immuno_target$symbol, ]
tmp$alias = immuno_target[row.names(tmp), "alias"]


plot_dot_with_sample_names = function(var1, var2, genesymbol, title, samplenames1, samplenames2){
    xx_base = c(1, 1.3)
    xx = rep(xx_base, each = length(var1))
    yy1 = var1
    yy2 = var2
    yy = c(yy1, yy2)

    ## ylim
    y1 = floor(min(yy)); 
    y2 = ceiling(max(yy)); 
    if(y1 > 0) { y1 = 0}
    
    plot("", type="n", xlab="", xlim=c(0.8, 1.7), ylim=c(y1, y2), axes=F, ylab="log2(FPKM)", 
         main=title)
    axis(2)
    text(xx_base, y1, c("T1", "T2"), adj=c(.5, 2), xpd=T)
    xx_rand = xx + runif(n = length(xx), min = -0.05, max = 0.05)
    points(xx_rand, yy, cex=1.5)
    segments(xx_rand[1:length(var1)], yy[1:length(var1)], xx_rand[(1+length(var1)):(2*length(var1))], 
             yy[(1+length(var1)):(2*length(var1))])
    shift = c(-0.1, 0.1)
    # mean
    segments(xx_base - 0.07 + shift, c(mean(yy1), mean(yy2)), 
             xx_base + 0.07 + shift, c(mean(yy1), mean(yy2)), col=1, xpd=T, lwd=2) 
    # top sd
    segments(xx_base - 0.05 + shift, c(mean(yy1) + sd(yy1), mean(yy2) + sd(yy2)), 
             xx_base + 0.05 + shift, c(mean(yy1) + sd(yy1), mean(yy2) + sd(yy2)), col=1, xpd=T) 
    # bot sd
    segments(xx_base - 0.05 + shift, c(mean(yy1) - sd(yy1), mean(yy2) - sd(yy2)), 
             xx_base + 0.05 + shift, c(mean(yy1) - sd(yy1), mean(yy2) - sd(yy2)), col=1, xpd=T) 
    # vertical line
    segments(xx_base + shift, c(mean(yy1) - sd(yy1), mean(yy2) - sd(yy2)), 
                 xx_base + shift, c(mean(yy1) + sd(yy1), mean(yy2) + sd(yy2)), col=1, xpd=T)
    text(xx_rand[1:length(var1)], yy1, samplenames1, adj=c(0,0), xpd=T, col=2)
    #text(xx_rand[(1+length(var1)):(2*length(var1))], yy2, samplenames2, adj=c(0,0), xpd=T)
}

immuno_target
head(tmp)

## tmp from above
options(repr.plot.width=12, repr.plot.height=6)
par(mfrow=c(2,4))
for(i in 1:nrow(immuno_target)){
    genesymbol = immuno_target$symbol[i]
    genealias  = immuno_target$alias[i]
    ss = as.logical(sss[genesymbol, ])
    var1 = as.numeric(t1[genesymbol, ss])
    var2 = as.numeric(t2[genesymbol, ss])
    pval = t2_t1_sum_df_rm_na[genesymbol, "pvalues"]
    nn = length(ss[ss])
    title = paste0(genesymbol, " (p=", round(pval, 2), " n=", nn, ")")
    sn = colnames(sss)[ss]
    sn = sub("DS_bla_", "", sn)
    sn = sub("_T1", "", sn)
    plot_dot_with_sample_names(var1, var2, genealias, title, sn, sn)
}


## tmp from above
options(repr.plot.width=12, repr.plot.height=6)
par(mfrow=c(2,4))
for(i in 1:nrow(immuno_target)){
    genesymbol = tmp$genesymbol[i]
    genealias = tmp$alias[i]
    ss = as.logical(sss[genesymbol, ])
    var1 = as.numeric(t1[genesymbol, ss])
    var2 = as.numeric(t2[genesymbol, ss])
    pval = t2_t1_sum_df_rm_na[genesymbol, "pvalues"]
    nn = length(ss[ss])
    title = paste0(genesymbol, " (p=", round(pval, 2), " n=", nn, ")")
    plot_dot(var1, var2, genealias, title)
}


t12 = cbind(t1, t2)
t12[immuno_target$symbol, ]

#GBP5: Guanylate Binding Protein 5
#PI3: Peptidase Inhibitor 3
#SAA1: Serum Amyloid A1
#TGM1: Transglutaminase 1 
#GBP6: Guanylate Binding Protein Family Member 6 
#ELF3: E74 Like ETS Transcription Factor 3
#CTLA4
tmp = c("GBP5", "PI3", "SAA1", "TGM1", "GBP6", "ELF3", "CTLA4")
t2_t1_sum_df[t2_t1_sum_df$symbol %in% tmp, ]
options(repr.plot.width=12, repr.plot.height=6)
par(mfrow=c(2,4))
for(i in 1:length(tmp)){
    genesymbol = tmp[i]
    genealias = tmp[i]
    var1 = as.numeric(t1[genesymbol, ])
    var2 = as.numeric(t2[genesymbol, ])
    pval = t2_t1_sum_df[genesymbol, "pvalues"]
    title = paste0(genesymbol, " p=", round(pval, 2))
    n1 = substr(colnames(t1), 8, 10)
    n2 = substr(colnames(t2), 8, 10)
    plot_dot_with_sample_names(var1, var2, genealias, title, n1, n2)
}


tmp

options(repr.plot.width=12, repr.plot.height=6)
par(mfrow=c(2,4))

for(i in 1:length(tmp)){
    genesymbol = tmp[i]
    genealias = tmp[i]
    var1 = as.numeric(t1[genesymbol, ])
    var2 = as.numeric(t2[genesymbol, ])
    pval = t2_t1_sum_df[genesymbol, "pvalues"]
    title = paste0(genesymbol, " p=", round(pval, 2))
    n1 = substr(colnames(t1), 8, 10)
    n2 = substr(colnames(t2), 8, 10)
    plot_dot(var1, var2, genealias, title)
}

genesymbols = immuno_target$symbol
t12 = cbind(t1, t2)

anno_cols = data.frame(row.names=colnames(t12), T1_T2 = rep("T1", ncol(t12)) , stringsAsFactors = F)
anno_cols$T1_T2[grep("T2", row.names(anno_cols))] = "T2"

options(repr.plot.width=6, repr.plot.height=4)
tmp = t12[genesymbols, ]
row.names(tmp) = immuno_target[genesymbols, "alias"]
pheatmap(tmp, scale="row", col = greenred(15), show_rownames = T,  annotation_col = anno_cols, cluster_cols = F, border_color = "grey0")
pheatmap(tmp, scale="row", col = greenred(15), show_rownames = T,  annotation_col = anno_cols, cluster_cols = T, border_color = "grey0")


dim(t2_t1_sum_df)

t2_t1_sig = t2_t1_sum_df[!is.na(t2_t1_sum_df$pvalues) & t2_t1_sum_df$pvalues < 0.01,]
dim(t2_t1_sig)
fpkm = cbind(t1, t2)
fpkm_sig = fpkm[row.names(fpkm) %in% t2_t1_sig$symbol, ]
dim(fpkm_sig)

head(t2_t1_sig)

library(pheatmap)
library(gplots)

options(repr.plot.width=6, repr.plot.height=6)
pheatmap(as.matrix(fpkm_sig), scale = "row", col = greenred(10), cluster_cols = T, 
         cluster_rows = T, show_rownames = F, hclust="complete", clustering_distance_cols = "maximum",
        border_color="grey0",  )

options(repr.plot.width=6, repr.plot.height=12)
pheatmap(as.matrix(fpkm_sig[1:100,]), scale = "row", col = greenred(10), cluster_cols = F, cluster_rows = T, show_rownames = T,
        border_color="grey0", cex=.8)

head(t2_t1_sig)

t12['KLK12',]

options(repr.plot.height = 4, repr.plot.width = 6)
par(mfrow=c(1,2))
tmp = as.numeric(t12['KLK12',])
p = t2_t1_sum_df['KLK12', 'pvalues']
plot_dot(tmp[1:10], tmp[11:20], "KLK12", p)


options(repr.plot.height = 4, repr.plot.width = 6)
par(mfrow=c(1,2))
tmp = as.numeric(t12['FANCF',])
p = t2_t1_sum_df['FANCF', 'pvalues']
plot_dot(tmp[1:10], tmp[11:20], "FANCF", p)
tmp = as.numeric(t12['FANCA',])
p = t2_t1_sum_df['FANCA', 'pvalues']
plot_dot(tmp[1:10], tmp[11:20], "FANCA", p)

t2_t1_sig[grep("CD", t2_t1_sig$symbol), ]

head(t2_t1_sum_df)

as.numeric(t1['CD209',])


options(repr.plot.width=9, repr.plot.height=3)
par(mfrow=c(1,3))
genesymbols = c("CD209", "CD74", "CD68")
for(genesymbol in genesymbols){
    var1 = as.numeric(t1[genesymbol,])
    var2 = as.numeric(t2[genesymbol,])
    pval = t2_t1_sum_df[genesymbol, "pvalues"]
    title = paste0(genesymbol, " p=", round(pval, 2))
    plot_dot(var1, var2, genesymbol, title)
}

load("../bcg/uc_basal_luminal_marks.RData")
tmp = t2_t1_sum_df[t2_t1_sum_df$symbol %in% base47, ]
#tmp
genesymbols = tmp$symbol
t12 = cbind(t1, t2)

anno_rows = data.frame(row.names = tmp$symbol, Mark = rep("basal", nrow(tmp)), stringsAsFactors = F)
anno_rows$Mark[row.names(anno_rows) %in% base47.luminal] = "luminal"
anno_cols = data.frame(row.names=colnames(t12), T1_T2 = rep("T1", ncol(t12)) , stringsAsFactors = F)
anno_cols$T1_T2[grep("T2", row.names(anno_cols))] = "T2"

options(repr.plot.width=6, repr.plot.height=7)
pheatmap(t12[genesymbols,], scale="row", col = greenred(15), show_rownames = T, annotation_row = anno_rows, 
         annotation_col = anno_cols, cluster_cols = T, cluster_rows = T, clustering_method = "ward.D", 
         border_color = "grey0", main="BASE47 classifiers")

load("../bcg/uc_basal_luminal_marks.RData")
tmp = t2_t1_sum_df[t2_t1_sum_df$symbol %in% mcconkey, ]
#tmp
genesymbols = tmp$symbol
t12 = cbind(t1, t2)

anno_rows = data.frame(row.names = tmp$symbol, Mark = rep("basal", nrow(tmp)), stringsAsFactors = F)
anno_rows$Mark[row.names(anno_rows) %in% mcconkey.luminal] = "luminal"
anno_cols = data.frame(row.names=colnames(t12), T1_T2 = rep("T1", ncol(t12)) , stringsAsFactors = F)
anno_cols$T1_T2[grep("T2", row.names(anno_cols))] = "T2"

options(repr.plot.width=6, repr.plot.height=4.5)
pheatmap(t12[genesymbols,], col = greenred(15), show_rownames = T, annotation_row = anno_rows, 
         annotation_col = anno_cols, cluster_cols = F, scale="row", clustering_method = "average", main="McConkey classifiers")

t2["GATA3",]
t1["FOXA1",]

dev.off()

t1[c("KLF8", "KLF14"),]
t2[c("KLF8", "KLF14"),]

options(jupyter.plot_mimetypes = c('image/jpeg'))
options(repr.plot.width=10, repr.plot.height=15)
par(mfrow=c(5,4))
genesymbols = row.names(t1[grep("KLF", row.names(t1)),])
for(genesymbol in genesymbols){
    var1 = as.numeric(t1[genesymbol,])
    var2 = as.numeric(t2[genesymbol,])
    pval = t2_t1_sum_df[genesymbol, "pvalues"]
    title = paste0(genesymbol, " n=", round(pval, 2))
    plot_dot(var1, var2, genesymbol, title)
}

options(jupyter.plot_mimetypes = c('image/jpeg'))
options(repr.plot.width=10, repr.plot.height=4)
par(mfrow=c(1,3))
genesymbols = c("GATA3", "FOXA1", "CTSE")
for(genesymbol in genesymbols){
    var1 = as.numeric(t1[genesymbol,])
    var2 = as.numeric(t2[genesymbol,])
    pval = t2_t1_sum_df[genesymbol, "pvalues"]
    title = paste0(genesymbol, " n=", round(pval, 2))
    plot_dot(var1, var2, genesymbol, title)
}

samples = unique(mut$Tumor_Sample_Barcode)

tmp = mut[mut$Tumor_Sample_Barcode == samples[1], ]
tmp_r = tmp$t_alt_count / tmp$t_depth
length(tmp_r)

shapiro.test(tmp_r)

mut_r = log2(mut$t_alt_count / mut$t_depth)
mut_r = mut$t_alt_count / mut$t_depth
plot(density(tmp_r))

## columns need to check sometimes
cn = c("Hugo_Symbol", "Chromosome", "Variant_Classification", "Variant_Type", "Reference_Allele", 
       "Tumor_Seq_Allele1", 'Tumor_Seq_Allele2','Tumor_Sample_Barcode', 't_depth', 't_ref_count', 't_alt_count',
       'n_depth', 'n_ref_count', 'n_alt_count')

mut = read.table("data_mutations_extended.txt", sep="\t", header=T, quote="", stringsAsFactors = F)
colnames(mut)


length(unique(mut$Tumor_Sample_Barcode))

snp = mut[mut$Variant_Type == "SNP", ]
sample_n = length(unique(snp$Tumor_Sample_Barcode))
#average coverage based on snp
# tumor cell coverage
mean(snp$t_depth)
mean(snp$n_depth)
mean(snp$t_depth + snp$n_depth)

unique(snp$Variant_Classification)


# nonsynonymous mutations
nonsyno_tag = c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Translation_Start_Site')
nonsyno = snp[snp$Variant_Classification %in% nonsyno_tag, ]
dim(nonsyno)
nonsyno_df = as.data.frame(table(nonsyno$Variant_Classification))
# synonymous mutations
syno = snp[snp$Variant_Classification == "Silent", ]
dim(syno)

nonsyno_df

tmp = table(nonsyno$Tumor_Sample_Barcode)
mean(tmp)
sd(tmp)

tmp = table(syno$Tumor_Sample_Barcode)
mean(tmp)
sd(tmp)

tmp = table(mut$Tumor_Sample_Barcode)
pat = substr(names(tmp), 8, 10)
tmp = aggregate(as.numeric(tmp), by=list(pat), mean)
tmp
mean(tmp$x)
sd(tmp$x)
range(tmp$x)
mean(tmp$x)*1000000/98754738

# see program/runback/mysql_from_UCSC.sh
# hg19 exon seq len is 98754738 
# see this for exon size discussion: http://seqanswers.com/forums/showthread.php?t=5298

3281/160

unique(mut$Variant_Classification)

unique(mut$Variant_Type)

colnames(mut)

tmp$Tumor_Sample_Barcode

require(lubridate)
days = 365*2
date = seq(as.Date("2000-01-01"), length = days, by = "day")
year = year(date)
month = month(date)
x1 = cumsum(rnorm(days, 0.05)) 
x2 = cumsum(rnorm(days, 0.05))
df1 = data.frame(date, year, month, x1, x2)
df_melt <- melt(df1, id = c("date", "year", "month"))
head(df_melt)

library(reshape2)

id = 187
tmp = mut[grep(id, mut$Tumor_Sample_Barcode), ]
tmp$VAF = tmp$t_alt_count/tmp$t_depth
#tmp = tmp[, c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode', 'VAF')]
range(tmp$VAF)


blca_nature_mut = c('TP53', 'KMT2D', 'ARID1A', 'KDM6A', 'PIK3CA', 'EP300', 'CDKN1A', 'RB1', 
                    'ERCC2', 'FGFR3', 'STAG2', 'ERBB3', 'FBXW7', 'RXRA', 'ELF3', 'NFE2L2', 
                    'TSC1', 'KLF5', 'TXNIP', 'FOXQ1', 'CDKN2A', 'RHOB', 'FOXA1', 'PAIP1', 'BTG2', 'HRAS', 'ZFP36L1', 'RHOA', 'CCND3')

blca_nature_cna = c('CDKN2A', 'E2F3', 'SOX4', 'CCND1', 'RB1', 'EGFR', 'PPARG', 'PVRL4', 'YWHAZ', 
                    'MDM2', 'ERBB2', 'CREBBP', 'NCOR1', 'YAP1', 'CCNE1', 'MYC', 'ZNF703', 'FGFR3', 'PTEN', 'MYCL', 'BCL2L1')
blca_nature=c(blca_nature_mut, blca_nature_cna)

onco_color = c(Mutation = "#26A818", Missense = "#26A818", Nonsense = "black", Splicing = "#ffaa00", 
               Frameshift = "#A05E35" , Promoter = "#2986E2", InFrame = "#F26529", Present = "darkorchid2", NotPresent = "#DCD9D3", 
               NotTested = "darkgrey", del = "red", LOH = "#D17878", homodel = "brown4", 
               CNLOH =  "deepskyblue", Amplification = "#EA2E49", Deletion = "#174D9D", Yes = "#155B6B", No = "#12C8F9", 
               Unknown = "azure1", Fusion =  "#D38C1F", Pathogenic="white")

mut_color = c(  #'Intron' = 'white',
                'Frame_Shift_Ins' = 'Frameshift',
                #'3\'UTR' = 'white',
                'Missense_Mutation' = 'Missense',
                #'3\'Flank' = 'white',
                #'Targeted_Region' = 'white',
                'In_Frame_Ins' = 'InFrame',
                #'5\'Flank' = 'Promoter',
                'Nonsense_Mutation' = 'Nonsense',
                #'Silent' = 'Mutation',
                'In_Frame_Del' = 'InFrame',
                'Frame_Shift_Del' = 'Frameshift',
                #'5\'UTR' = 'Promoter',
                'Splice_Site' = 'Splicing',
                #'RNA' = 'Mutation',
                #'IGR' = 'Mutation',
                'Nonstop_Mutation' = 'Missense_mutation',
                'Translation_Start_Site' = 'Nonsense'
)


focus_mut = c('Frame_Shift_Ins', 'Frame_Shift_Del',  'In_Frame_Ins', 'In_Frame_Del', 
'Missense_Mutation','Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site')


ids = unique(substr(unique(mut$Tumor_Sample_Barcode), 8, 10))
ids

options(repr.plot.width=9, repr.plot.height=9)
par(mfrow=c(2,2))
for(id in ids){
    tmp = mut[grep(id, mut$Tumor_Sample_Barcode), ]
    tmp$VAF = tmp$t_alt_count/tmp$t_depth
    #tmp = tmp[, c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode', 'VAF')] 
    tmp$col = onco_color[mut_color[tmp$Variant_Classification]]
    tmp = tmp[tmp$Variant_Classification %in% focus_mut, ]

    tmp2 = dcast(tmp, Hugo_Symbol+Chromosome+Start_Position+Variant_Classification+Variant_Type+col ~ Tumor_Sample_Barcode, sum, value.var = 'VAF')

    id1 = paste0('DS_bla_', id, '_T1')
    id2 = paste0('DS_bla_', id, '_T2')
    plot(tmp2[,id1], tmp2[,id2], cex=.7, xlab=paste0(id1, " VAF"), ylab=paste0(id2, " VAF"), col = tmp2$col, xlim=c(0,1), ylim=c(0,1))
    tmp3 = tmp2[tmp2$Hugo_Symbol %in% blca_nature & tmp2[,id1] != 0 & tmp2[,id2] !=0,]
    if(nrow(tmp3) > 0){
        text(tmp3[, id1], tmp3[, id2], tmp3$Hugo_Symbol, adj=c(-0.1,0), cex = .7, col = tmp3$col)
        points(rep(0.75, nrow(tmp3)), seq(from=1, by = -0.08, length.out = nrow(tmp3)), col=tmp3$col)
        text(rep(0.75, nrow(tmp3)), seq(from=1, by = -0.08, length.out = nrow(tmp3)), tmp3$Hugo_Symbol, adj=c(-0.2,.5))
        }
    }
plot.new()
legend('topright', legend=names(mut_color), col = onco_color[mut_color], cex=.7, pch=1)


options(repr.plot.width=9, repr.plot.height=9)
par(mfrow=c(2,2))
for(id in ids){
    tmp = mut[grep(id, mut$Tumor_Sample_Barcode), ]
    tmp$VAF = tmp$t_alt_count/tmp$t_depth
    #tmp = tmp[, c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode', 'VAF')] 
    tmp$col = onco_color[mut_color[tmp$Variant_Classification]]
    tmp = tmp[tmp$Variant_Classification %in% focus_mut, ]

    tmp2 = dcast(tmp, Hugo_Symbol+Chromosome+Start_Position+Variant_Classification+Variant_Type+col ~ Tumor_Sample_Barcode, sum, value.var = 'VAF')

    id1 = paste0('DS_bla_', id, '_T1')
    id2 = paste0('DS_bla_', id, '_T2')
        
    tmp2 = tmp2[tmp2$Hugo_Symbol %in% blca_nature,]

    plot(tmp2[,id1], tmp2[,id2], cex=.7, xlab=paste0(id1, " VAF"), ylab=paste0(id2, " VAF"), col = tmp2$col, xlim=c(0,1), ylim=c(0,1))
    tmp3 = tmp2[tmp2$Hugo_Symbol %in% blca_nature,]
    if(nrow(tmp3) > 0){
        text(tmp3[, id1], tmp3[, id2], tmp3$Hugo_Symbol, adj=c(-0.1,0), cex = .7, col = tmp3$col)
        points(rep(0.75, nrow(tmp3)), seq(from=1, by = -0.08, length.out = nrow(tmp3)), col=tmp3$col)
        text(rep(0.75, nrow(tmp3)), seq(from=1, by = -0.08, length.out = nrow(tmp3)), tmp3$Hugo_Symbol, adj=c(-0.2,.5))
        }
    }
plot.new()
legend('topright', legend=names(mut_color), col = onco_color[mut_color], cex=.7, pch=1)


library(ggplot2)
library(data.table)

sig_file = "Proj_06155___SOMATIC_FACETS.vep.noheader.Ref_Tri.sig"
sig <- fread(sig_file)
setnames(sig, names(sig), gsub("Signature.", "", names(sig)))
sig[, APOBEC := `2` + `13`]
sig[, `2` := NULL]
sig[, `13` := NULL]

sig[, c("Tumor", "sample") := tstrsplit(`Sample Name`, "_")[3:4]]
sig[, `Sample Name` := NULL]
#sig[, Tumor := as.integer(Tumor) - 184]
sig[ , sample := factor(sample, levels = c("M", "T2", "T1"))]

MSI_sigs <- c("6", "15", "20", "26")
setnames(sig, MSI_sigs, paste0("MSI", MSI_sigs))
### consolidate other columns
other_cols <- setdiff(names(sig), c("1", "APOBEC", "Tumor", "sample", "Number of Mutations", paste0("MSI", MSI_sigs)))
sig[, others := rowSums(.SD), .SDcol = other_cols]
sig[, other_cols := NULL, with = F]

### rename columns
setnames(sig, "1", "Age")
sig = melt.data.table(sig, id.vars = c("Tumor", "sample", "Number of Mutations"))

old <- theme_set(theme_grey(base_size = 8))

ggplot(sig, aes(sample, fill = variable, weight = value)) +
    geom_bar(stat = "count") +
    coord_flip() +
    facet_grid(Tumor ~ ., space = "free_x", scales = "free_x", drop = T) +
    #scale_fill_manual(values = rainbow(31))
    scale_fill_manual("", values = c("orange",'#a6bddb','#74a9cf','#2b8cbe','#045a8d', "dark green", "grey")) +
    xlab("") +
    ylab("fraction of mutations") +
    scale_fill_brewer()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8))


save.image()

library(facets)

xx <- preProcSample(mut)
oo <- procSample(xx, cval = 300)
emcncf(oo)


cn2 = c('Verification_Status', 'Validation_Status', 'Mutation_Status', 
        'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score')
#mut[1000:1010, cn2]

colnames(mut)

RunAbsolute(seg.dat.fn = "blca_cmo_06155_2016_data_cna_hg19.seg", )

table(mut$Variant_Classification)

tmp = mut$Hugo_Symbol

mut$t12 = substr(mut$Tumor_Sample_Barcode, 12, 13)
tmp = table(mut[, c("Hugo_Symbol", "t12")])
tmp = as.data.frame(tmp)
mut_ordered_wide = dcast(tmp, Hugo_Symbol ~ t12, value.var = "Freq")
mut_ordered_wide$t2_t1 = mut_ordered_wide$T2 - mut_ordered_wide$T1
mut_ordered_wide = mut_ordered_wide[order(mut_ordered_wide$t2_t1, decreasing = T),]
tmp = mut_ordered_wide[abs(mut_ordered_wide$t2_t1) > 7, ]
dim(tmp)
head(tmp)

mut$t12 = "T1"
mut$t12[grep("T2", mut$Tumor_Sample_Barcode)] = "T2"
tmp3 = table(mut$Hugo_Symbol, mut$Variant_Classification, mut$t12)
tmp4 = as.data.frame(tmp3)
tmp4 = dcast(tmp4, Var1 + Var3 ~ Var2, value.var = "Freq")
tmp5 = apply(tmp4[, 3:20], 1, sum)


tmp6 = tmp4[tmp5>7, ]
tmp6 = tmp4[tmp4$Var1 %in% tmp6$Var1,]
tmp6$Var3 = as.character(tmp6$Var3)
tmp6t1 = tmp6[grep("T1", tmp6$Var3),]
tmp6t2 = tmp6[grep("T2", tmp6$Var3),]
tmp6t12 = tmp6t1
for(i in 1:nrow(tmp6t1)){
    tmp6t12[i, 3:20] = tmp6t2[i, 3:20] - tmp6t1[i, 3:20]
    tmp6t12$Var3[i] = "T2-T1"
}

tmp6t12 = tmp6t12[order(tmp6t12$Missense_Mutation, decreasing = T),]
tmp6t12 = tmp6t12[order(abs(tmp6t12$Missense_Mutation), abs(tmp6t12$Silent), decreasing = T),]
head(tmp6t12[,c("Var1", "Var3", "Missense_Mutation", "Silent", "Frame_Shift_Ins", "Intron")])
write.table(tmp6t12, file="mutation_T2_T1.xls", sep="\t")

head(tmp6t12[,c("Var1", "Var3", "Missense_Mutation", "Silent", "Frame_Shift_Ins", "Intron")], 20)

tmp2 = mut[mut$Hugo_Symbol %in% tmp$Hugo_Symbol,]
tmp2$t12 = "T1"
tmp2$t12[grep("T2", tmp2$Tumor_Sample_Barcode)] = "T2"
tmp3 = table(tmp2$Hugo_Symbol, tmp2$Variant_Classification, tmp2$t12)
tmp4 = as.data.frame(tmp3)
tmp4 = dcast(tmp4, Var1 + Var3 ~ Var2, value.var = "Freq")
tmp4

head(mut_ordered_wide[grep("CD", mut_ordered_wide$Hugo_Symbol),])

x = colnames(mut)
names(x) = 1:length(x)
x

table(mut$t12, mut$Variant_Type)

snp = mut[mut$Variant_Type == "SNP", ]
head(snp[,cn])

table(snp$t12, snp$Variant_Classification)

mut$snp = paste(mut$Reference_Allele, mut$Tumor_Seq_Allele2, sep="->")
mut_snp = mut[mut$Variant_Type == "SNP",]
snp_table = table(mut_snp$snp, mut_snp$Tumor_Sample_Barcode)
snp_table = as.data.frame.matrix(snp_table, )
snp_table = snp_table[, -7]
dim(snp_table)

snp_mean = apply(snp_table, 1, function(x){
    tmp = c(mean(x[seq(1, 22, 2)]), mean(x[seq(2, 22, 2)]));
})
snp_mean = t(snp_mean)
snp_mat = as.data.frame(snp_mean, stringsAsFactors = F)
colnames(snp_mat) = c("T1_mean", "T2_mean")
snp_mat

count = as.data.frame.matrix(table(mut_snp$snp, mut_snp$t12))
colnames(count) = c("T1_Count", "T2_Count")
snp_mat = cbind(snp_mat, count)
snp_mat

sum(snp_mat[,3])

sum(snp_mat[,4])

snp_mat$T1_percent = snp_mat$T1_Count / sum(snp_mat[,3])
snp_mat$T2_percent = snp_mat$T2_Count / sum(snp_mat[,4])

pvalues = apply(snp_table, 1, function(x){
    tmp = t.test(x[seq(1, 22, 2)], x[seq(2, 22, 2)]);
    tmp$p.value
})
snp_mat$pvalue = pvalues
snp_mat



x = snp_table[1,]

pairwise.t.test(as.numeric(x), rep(c(1, 2), times = 11));

pvalues_pair = apply(snp_table, 1, function(x){
    tmp = pairwise.t.test(x, rep(c(1, 2), times = 11));
    tmp$p.value
})
snp_mat$pvalue_pairwise = pvalues_pair
snp_mat

shapiro.test(mut$t_)

snp = mut[mut$Variant_Type == "SNP", ]


snp_gr = VRanges(
    seqnames = paste0("chr", snp$Chromosome),
    ranges = IRanges(start = snp$Start_Position, end = snp$End_Position),
    ref = snp$Reference_Allele,
    alt = snp$Tumor_Seq_Allele2,
    sampleNames = snp$Tumor_Sample_Barcode,
    seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19),
    study = snp$t12)
snp_gr = mutationContext(snp_gr, BSgenome.Hsapiens.UCSC.hg19)
#snp_gr = motifMatrix(snp_gr, group = "study", normalize = TRUE)
#head(snp_gr)

#snp_gr = as.data.frame(snp_gr)

#png(file="mut_snp_T1T2_signatures.pdf", width=9, height=2.5)
options(repr.plot.width = 9, repr.plot.height = 3)
plotMutationSpectrum(snp_gr, "study", colorby = "alteration")
#dev.off()


library(sciClone)
library(clonevol)
library(fishplot)
data(aml1)

see here https://www.biostars.org/p/191247/
feed the results from sciClone to clonevol (https://github.com/hdng/clonevol)
https://github.com/chrisamiller/fishplot

head(mut[, cn])

segfile = '/ifs/work/solitlab/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/Proj_07813_DF/r_001/variants/copyNumber/facets/Proj_07813_DF_facets_merge_hisens.seg'
cnv = fread(segfile)
cnv
colnames(cnv) = c('ID', 'chr', 'start', 'stop', 'num_mark', 'segment_mean')
cnv = cnv[,-5]
cnv = cnv[!is.na(segment_mean), ]
cnv

## test for single patient
# read in vaf data from three related tumors
#format is 5 column, tab delimited: 
#chr, pos, ref_reads, var_reads, vaf (variant allele fraction)
scc.maf[, vaf := 100 * t_alt_count / t_depth]
mut_clone =scc.maf[Variant_Type == "SNP", c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", "t_ref_count", "t_alt_count", 'vaf')]
colnames(mut_clone) = c("Tumor_Sample_Barcode", "chr", "pos", "ref_reads", "var_reads", 'vaf')

id = 194
samples = unique(scc.maf[grep(id, Tumor_Sample_Barcode), Tumor_Sample_Barcode]); samples
for(i in 1:length(samples)){
	tmp = scc.maf[Tumor_Sample_Barcode == samples[i] & Variant_Type == "SNP", c( "Chromosome", "Start_Position", "t_ref_count", "t_alt_count", 'vaf')]
	colnames(tmp) = c("chr", "pos", "ref_reads", "var_reads", 'vaf')
	assign(samples[i], tmp)
}
get(samples[1])

cnv.samples = unique(cnv[grep(id, ID), ID]); cnv.samples
cnv.clone = cnv[ID %in% cnv.samples,]

vafs.list = lapply(samples, get)
names(vafs.list) = samples
vafs.list

library(sciClone)
sciClone(vafs = lapply(samples, get), sampleNames = samples, verbose=F)

id1 = paste0("DS_bla_", id, "_T1")
id2 = paste0("DS_bla_", id, "_T2")
mut_t1 = mut_clone[mut_clone$Tumor_Sample_Barcode == id1, ]
mut_t2 = mut_clone[mut_clone$Tumor_Sample_Barcode == id2, ]
mut_t1 = mut_t1[, -1]
mut_t2 = mut_t2[, -1]

cnv_t1 = cnv[cnv$ID == id1, ]
cnv_t2 = cnv[cnv$ID == id2, ]

#set sample names
names = paste0(c("T1_", "T2_"), id)

sc = sciClone(vafs = list(mut_t1, mut_t2),  sampleNames = names, verbose = F)
sc_df = sc@vafs.merged       

#sc_df = sc@vafs.merged
sel_cn = c("cluster", paste0(c("T1_", "T1_", "T2_", "T2_"), id, c(".vaf", ".depth", ".vaf", ".depth")))
sc_sel = sc_df[!is.na(sc_df$cluster), sel_cn]
dim(sc_sel)
vaf.col.names <- grep(".vaf", colnames(sc_sel), value=TRUE)

sc_clonevol <- infer.clonal.models(variants=sc_sel,
            cluster.col.name="cluster",
            vaf.col.names=vaf.col.names,
            subclonal.test="bootstrap",
            subclonal.test.model="non-parametric",
            cluster.center="mean",
            num.boots=1000,
            founding.cluster=1,
            min.cluster.vaf=0.01,
            p.value.cutoff=0.01,
            alpha=0.1,
            random.seed=63108)

f = generateFishplotInputs(results=sc_clonevol)

fishes = createFishPlotObjects(f)
#plot with fishplot
#pdf('fish.pdf', width=8, height=5)
for (i in 1:length(fishes)){
    fish = layoutClones(fishes[[i]])
    fish = setCol(fish,f$clonevol.clone.colors)
    fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
             vlines=seq(1, length(names)), vlab=names, pad.left=0.5)
}
#dev <- dev.off()

### 188

## test for single patient
# read in vaf data from three related tumors
#format is 5 column, tab delimited: 
#chr, pos, ref_reads, var_reads, vaf (variant allele fraction)
mut_clone = mut[mut$Variant_Type == "SNP", c("Tumor_Sample_Barcode", "Chromosome", 
                                             "Start_Position", "t_ref_count", "t_alt_count")]
colnames(mut_clone) = c("Tumor_Sample_Barcode", "chr", "pos", "ref_reads", "var_reads")
mut_clone$vaf = 100 * mut_clone$var_reads / (mut_clone$var_reads + mut_clone$ref_reads)

id = 188
id1 = paste0("DS_bla_", id, "_T1")
id2 = paste0("DS_bla_", id, "_T2")
id3 = paste0("DS_bla_", id, "_M1")
mut_t1 = mut_clone[mut_clone$Tumor_Sample_Barcode == id1, ]
mut_t2 = mut_clone[mut_clone$Tumor_Sample_Barcode == id2, ]
mut_m1 = mut_clone[mut_clone$Tumor_Sample_Barcode == id2, ]
mut_t1 = mut_t1[, -1]
mut_t2 = mut_t2[, -1]
mut_m1 = mut_m1[, -1]

cnv_t1 = cnv[cnv$ID == id1, ]
cnv_t2 = cnv[cnv$ID == id2, ]
cnv_m1 = cnv[cnv$ID == id3, ]

#set sample names
names = paste0(c("T1_", "T2_"), id)
names = paste0(c("T1_", "T2_", "M1_"), id)

sc = sciClone(vafs = list(mut_t1, mut_t2, mut_m1),  sampleNames = names, verbose = F)
sc_df = sc@vafs.merged       

#sc_df = sc@vafs.merged
sel_cn = c("cluster", paste0(c("T1_", "T1_", "T2_", "T2_", "M1_"), id, c(".vaf", ".depth", ".vaf", ".depth", ".vaf", ".depth")))
sc_sel = sc_df[!is.na(sc_df$cluster), sel_cn]
dim(sc_sel)
vaf.col.names <- grep(".vaf", colnames(sc_sel), value=TRUE)

sc_clonevol <- infer.clonal.models(variants=sc_sel,
            cluster.col.name="cluster",
            vaf.col.names=vaf.col.names,
            subclonal.test="bootstrap",
            subclonal.test.model="non-parametric",
            cluster.center="mean",
            num.boots=1000,
            founding.cluster=1,
            min.cluster.vaf=0.01,
            p.value.cutoff=0.01,
            alpha=0.1,
            random.seed=63108)

f = generateFishplotInputs(results=sc_clonevol)

fishes = createFishPlotObjects(f)
#plot with fishplot
#pdf('fish.pdf', width=8, height=5)
for (i in 1:length(fishes)){
    fish = layoutClones(fishes[[i]])
    fish = setCol(fish,f$clonevol.clone.colors)
    fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
             vlines=seq(1, length(names)), vlab=names, pad.left=0.5)
}
#dev <- dev.off()

head(mut_clone)

## test for single patient
# read in vaf data from three related tumors
#format is 5 column, tab delimited: 
#chr, pos, ref_reads, var_reads, vaf (variant allele fraction)
mut_clone = mut[mut$Variant_Type == "SNP", c("Tumor_Sample_Barcode", "Chromosome", 
                                             "Start_Position", "t_ref_count", "t_alt_count")]
colnames(mut_clone) = c("Tumor_Sample_Barcode", "chr", "pos", "ref_reads", "var_reads")
mut_clone$vaf = 100 * mut_clone$var_reads / (mut_clone$var_reads + mut_clone$ref_reads)

id = 188
id1 = paste0("DS_bla_", id, "_T1")
id2 = paste0("DS_bla_", id, "_T2")
id3 = paste0("DS_bla_", id, "_M1")
mut_t1 = mut_clone[mut_clone$Tumor_Sample_Barcode == id1, ]
mut_t2 = mut_clone[mut_clone$Tumor_Sample_Barcode == id2, ]
mut_t3 = mut_clone[mut_clone$Tumor_Sample_Barcode == id2, ]
mut_t1 = mut_t1[, -1]
mut_t2 = mut_t2[, -1]
mut_t3 = mut_t3[, -1]

# copy number variation data
#4 columns - chr, start, stop, segment_mean   
# Values: -2 = homozygous deletion; -1 = hemizygous deletion; 0 = neutral / no change; 1 = gain; 2 = high level amplification.
# this is from facets
# gata3 chr10:8,096,667-8,117,164
# 

cnv_t1 = cnv[cnv$ID == id1, ]
cnv_t2 = cnv[cnv$ID == id2, ]
cnv_t3 = cnv[cnv$ID == id3, ]

#set sample names
names = paste0(c("T1_", "T2_"), id)
names = paste0(c("T1_", "T2_", "M1_"), id)

#read in regions to exclude (commonly LOH)
#format is 3-col bed
# regions = read.table("data/exclude.loh")

# no CNV change with 185 samples
#sc = sciClone(vafs = mut, copyNumberCalls = cnv, sampleNames = "ALL",  regionsToExclude=reg1)
sc = sciClone(vafs = list(mut_t1, mut_t2),  sampleNames = names, verbose = F)
#sc = sciClone(vafs = list(mut_t1, mut_t2),  sampleNames = names, verbose = F)
#sc = sciClone(vafs = list(mut_t1, mut_t2), copyNumberCalls = list(cnv_t1, cnv_t2), sampleNames = names, verbose = F)
sc_df = sc@vafs.merged       

#create output
writeClusterTable(sc, paste0("results/sciClone_", id))
#sc.plot1d(sc,"results/sciClone_all.pdf")


#sc_df = sc@vafs.merged
sel_cn = c("cluster", paste0(c("T1_", "T1_", "T2_", "T2_"), id, c(".vaf", ".depth", ".vaf", ".depth")))
sc_sel = sc_df[!is.na(sc_df$cluster), sel_cn]
dim(sc_sel)
vaf.col.names <- grep(".vaf", colnames(sc_sel), value=TRUE)

sc_clonevol <- infer.clonal.models(variants=sc_sel,
            cluster.col.name="cluster",
            vaf.col.names=vaf.col.names,
            subclonal.test="bootstrap",
            subclonal.test.model="non-parametric",
            cluster.center="mean",
            num.boots=1000,
            founding.cluster=1,
            min.cluster.vaf=0.01,
            p.value.cutoff=0.01,
            alpha=0.1,
            random.seed=63108)




f = generateFishplotInputs(results=sc_clonevol)

fishes = createFishPlotObjects(f)
#plot with fishplot
#pdf('fish.pdf', width=8, height=5)
for (i in 1:length(fishes)){
    fish = layoutClones(fishes[[i]])
    fish = setCol(fish,f$clonevol.clone.colors)
    fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
             vlines=seq(1, length(names)), vlab=names, pad.left=0.5)
}
#dev <- dev.off()

generateFishplotInputs <- function(results, rescale=T, samples=NULL){
  #no results, punt
  if (is.null(results$matched)){return(NULL)}

  #if specific samples aren't given, use them all
  if (is.null(samples)){
    samples = names(results$models)
  }

  #store the clonevol results in a list
  res = list(samples=samples, clonevol.clone.names=NULL, clonevol.clone.colors=NULL,
    timepoints=seq(1, length(samples)), num.models = nrow(results$matched$index),
    parents=list(), cell.fractions=list(), all=list())
  
  
  clonevol.clone.names = NULL
  clone.nums = NULL
  clonevol.clone.colors = NULL

  #create the needed inputs to fishplot
  for (i in 1:nrow(results$matched$index)){
    vv = NULL
    for (s in samples){
      v = results$models[[s]][[results$matched$index[i, s]]]
      if (rescale){v = rescale.vaf(v)}
      v = v[, c('lab', 'vaf', 'parent', 'color')]
      
      ## scale vaf and make cell.frac
      max.vaf = max(v$vaf)
      scale = 0.5/max.vaf*2*100
      v$vaf = v$vaf*scale
      v$vaf[v$vaf > 100] = 100# safeguard against rounding error making some vaf slightly > 100
      
      colnames(v) = c('clone', s , 'parent', 'color')
      v = v[!is.na(v$parent) & v$clone != '0',]
      if (is.null(vv)){vv = v}else{vv = merge(vv, v, all=T)}
    }
    for (s in samples){
      vv[is.na(vv[[s]]),s] = 0
    }
    vv = vv[order(as.integer(vv$clone)),]
    vv$parent[vv$parent == '-1'] = 0
    rownames(vv) = vv$clone
    
    ## fishplot requires clones to be named in sequential order. Do that, but
    ## store the clonevol-generated names and colors for pass-through
    if (is.null(clone.nums)){
      clone.nums = c(0, seq(1, nrow(vv)))
      names(clone.nums) = c(0, vv$clone)
      
      clonevol.clone.names = names(clone.nums)
      names(clonevol.clone.names) = as.character(clone.nums)
      res$clonevol.clone.names = clonevol.clone.names[-1]
      
      clonevol.clone.colors = c( 'white', vv$color)
      names(clonevol.clone.colors) = as.character(clone.nums)
      res$clonevol.clone.colors = clonevol.clone.colors[-1]

    }
      print(vv)
    vv$clone = clone.nums[vv$clone]
    vv$parent = clone.nums[vv$parent]
    
    par = vv$parent
    frac = vv[, samples]
    res$parents[[i]] = par
    res$cell.fractions[[i]] = as.matrix(frac)
    res$all[[i]] = vv
  }  
  return(res)
}


    vv$clone = clone.nums[vv$clone]
    vv$parent = clone.nums[vv$parent]
    

names

assign()

sc.plot2d(sc,"results/sciClone_all_2d.pdf")


IDs = unique(substr(sub("DS_bla_", "", unique(mut$Tumor_Sample_Barcode)), 1, 3))
#read in vaf data from three related tumors
#format is 5 column, tab delimited: 
#chr, pos, ref_reads, var_reads, vaf (variant allele fraction)
#for(i in 1:length(IDs)){
    i = 185
    id1 = paste0("DS_bla_", i, "_T1")
    id2 = paste0("DS_bla_", i, "_T2")
    mut_t1 = mut[mut$Tumor_Sample_Barcode == id1, )]
    mut_t2 = mut[mut$Tumor_Sample_Barcode == id2, )]

    cnv_t1 = cnv[cnv$ID == id1, ]
    cnv_t2 = cnv[cnv$ID == id2, ]

    #set sample names
    names = paste0(c("T1_", "T2_"), i)
    
    sc = sciClone(vafs=list(mut_t1, mut_t2),
             copyNumberCalls=list(cnv_t1, cnv_t2), 
             sampleNames=names)
    #         regionsToExclude=reg1)
    #create output
    writeClusterTable(sc, paste0("results/T_", i))
    sc.plot1d(sc,paste0("results/T_", i, "_plot.pdf")
#}

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4374044/

cnv_t2$segment_mean

mut_indel = mut[mut$Variant_Type %in% c("DEL", "INS"), ];
table(mut_indel$Variant_Classification)

mut_indel = mut[mut$Variant_Type %in% c("DEL", "INS"), ];
indel_table = table(mut_indel$Variant_Type, mut_indel$Tumor_Sample_Barcode)
indel_table = indel_table[,-7]
indel_table = as.data.frame.matrix(indel_table)

indel_mean = apply(indel_table, 1, function(x){
    tmp = c(mean(x[seq(1, 22, 2)]), mean(x[seq(2, 22, 2)]));
})
indel_mean

indel_mat = t(indel_mean)
indel_mat = as.data.frame(indel_mat, stringsAsFactors = F)
colnames(indel_mat) = c("T1_mean", "T2_mean")
indel_mat

pvalues = apply(indel_table, 1, function(x){
    tmp = t.test(x[seq(1, 22, 2)], x[seq(2, 22, 2)]);
    tmp$p.value
})
indel_mat$pvalue = pvalues
indel_mat


t(indel_table)

options(repr.plot.width=8)
par(mfrow=c(1,3))
xx = as.numeric(substr(colnames(indel_table), 13, 13))
plot(xx, indel_table[1,], xlim=c(0.5, 2.5), xlab="Deletion mutations in T1/T2", axes=F, ylab="numbers", ylim=c(0, 500))
axis(2); text(c(1,2), 0, c("T1", "T2"), xpd=T, adj=c(0.5, 2))
plot(xx, indel_table[2,], xlim=c(0.5, 2.5), xlab="Insertion mutations in T1/T2", axes=F, ylab="numbers")
axis(2); text(c(1,2), 0, c("T1", "T2"), xpd=T, adj=c(0.5, 2))
plot(xx, indel_table[3,], xlim=c(0.5, 2.5), xlab="SNP in T1/T2", axes=F, ylab="numbers", ylim=c(0, 3000))
axis(2); text(c(1,2), 0, c("T1", "T2"), xpd=T, adj=c(0.5, 2))

Ozone

attach(airquality)
Month <- factor(Month, labels = month.abb[5:9])
pairwise.t.test(Ozone, Month)
pairwise.t.test(Ozone, Month, p.adj = "bonf")
pairwise.t.test(Ozone, Month, pool.sd = FALSE)
detach()


as.data.frame(table(mut_indel$Variant_Type, mut_indel$Variant_Classification, mut_indel$Tumor_Sample_Barcode))


library(xlsx)

tmp1 = table(mut$Variant_Classification, mut$t12, mut$Variant_Type)
tmp1 = as.data.frame(tmp1)
write.xlsx(tmp1, file="mut_type.xlsx")

table(mut$Variant_Type)

data.frame(t(as.list(tmp2)))

tmp = row.names(t2_t1_sum_df[!is.na(t2_t1_sum_df$pvalues) & t2_t1_sum_df$pvalues < 0.05 & !is.na(t2_t1_sum_df$log2FC) & abs(t2_t1_sum_df$log2FC) > 1,])
length(intersect(tmp, row.names(blca.cc.fil2.log.z)))
nrow(blca.cc.fil2.log.z)
length(tmp)

tcga_sc_id = c('TCGA-BT-A0YX-01A', 'TCGA-BT-A20U-01A', 'TCGA-BT-A2LD-01A', 'TCGA-C4-A0F1-01A', 'TCGA-C4-A0F7-01A', 'TCGA-CU-A0YN-01A', 
  'TCGA-DK-A2I2-01A', 'TCGA-FD-A3B5-01A', 'TCGA-FD-A3N5-01A', 'TCGA-G2-A2ES-01A', 'TCGA-G2-A3IB-01A', 'TCGA-GC-A3I6-01A', 
  'TCGA-GD-A3OS-01A', 'TCGA-BT-A42E', 'TCGA-GU-A42Q', 'TCGA-FD-A43Y', 'TCGA-FD-A5BU', 'TCGA-K4-A4AC', 'TCGA-FD-A5BY', 'TCGA-PQ-A6FI', 
  'TCGA-GU-A766', 'TCGA-CU-A72E', 'TCGA-E7-A7XN', 'TCGA-XF-A8HE', 'TCGA-YC-A89H', 'TCGA-E7-A97P', 'TCGA-XF-A8HH', 'TCGA-ZF-A9RE', 
  'TCGA-ZF-AA4W', 'TCGA-4Z-AA80', 'TCGA-4Z-AA82', 'TCGA-4Z-AA89', 'TCGA-XF-A9SJ', 'TCGA-XF-A9T4', 'TCGA-XF-A9T8', 'TCGA-ZF-AA53', 
  'TCGA-XF-AAME', 'TCGA-XF-AAMH', 'TCGA-XF-AAMT', 'TCGA-XF-AAN2', 'TCGA-XF-AAN5', 'TCGA-ZF-A9RD', 'TCGA-ZF-A9RG', 'TCGA-BT-A20X-01A', 
  'TCGA-FD-A3B4-01A', 'TCGA-FD-A3N6-01A')

tcga_sc_id = sub("-01A", "", tcga_sc_id)
tcga_sc_id

tcga_basal_id = c('TCGA-DK-A3WW', 'TCGA-SY-A9G5', 'TCGA-DK-A2I4', 'TCGA-E7-A7XN', 'TCGA-XF-A9T8', 'TCGA-XF-A9T5', 'TCGA-ZF-AA4V', 'TCGA-BT-A0YX', 'TCGA-DK-AA6S', 
  'TCGA-UY-A9PB', 'TCGA-4Z-AA7W', 'TCGA-4Z-AA84', 'TCGA-BT-A3PJ', 'TCGA-XF-A8HD', 'TCGA-ZF-AA4W', 'TCGA-BT-A20J', 'TCGA-FD-A3SS', 'TCGA-XF-AAN2', 
  'TCGA-XF-A9SY', 'TCGA-DK-A3IU', 'TCGA-G2-A2ES', 'TCGA-XF-A8HE', 'TCGA-XF-A9T3', 'TCGA-DK-A3IN', 'TCGA-DK-AA6L', 'TCGA-FD-A3SN', 'TCGA-FD-A3SO', 
  'TCGA-FD-A5BX', 'TCGA-GC-A3RC', 'TCGA-XF-A9SI', 'TCGA-HQ-A5NE', 'TCGA-K4-A5RJ', 'TCGA-K4-A6FZ', 'TCGA-BT-A3PK', 'TCGA-DK-AA6R', 'TCGA-FD-A3N5', 
  'TCGA-E7-A7DV', 'TCGA-UY-A8OB', 'TCGA-XF-AAN5', 'TCGA-XF-A9T6', 'TCGA-XF-A9SJ', 'TCGA-UY-A78P', 'TCGA-G2-A2EJ', 'TCGA-DK-A2I6', 'TCGA-FD-A3B6', 
  'TCGA-G2-A2EF', 'TCGA-XF-A9SM', 'TCGA-GC-A3YS', 'TCGA-ZF-A9RF', 'TCGA-CF-A1HS', 'TCGA-G2-AA3C', 'TCGA-4Z-AA7Q', 'TCGA-CU-A3KJ', 'TCGA-GC-A6I1', 
  'TCGA-DK-AA74', 'TCGA-E5-A2PC', 'TCGA-4Z-AA86', 'TCGA-4Z-AA81', 'TCGA-ZF-AA58', 'TCGA-DK-A1AB', 'TCGA-UY-A78L', 'TCGA-E7-A97P', 'TCGA-GU-A766', 
  'TCGA-ZF-AA54', 'TCGA-K4-A5RH', 'TCGA-GU-A42Q', 'TCGA-ZF-A9RN', 'TCGA-BT-A42F', 'TCGA-GC-A3I6', 'TCGA-BT-A20O', 'TCGA-ZF-AA53', 'TCGA-C4-A0F1', 
  'TCGA-BT-A42E', 'TCGA-ZF-A9RD', 'TCGA-GU-AATQ', 'TCGA-FD-A3NA', 'TCGA-FT-A61P', 'TCGA-UY-A9PA', 'TCGA-FD-A3B3', 'TCGA-XF-AAN3', 'TCGA-ZF-AA4N', 
  'TCGA-CU-A0YN', 'TCGA-FD-A6TH', 'TCGA-BL-A5ZZ', 'TCGA-GV-A40E', 'TCGA-GC-A3WC', 'TCGA-K4-A4AC', 'TCGA-PQ-A6FI', 'TCGA-FD-A3B5', 'TCGA-FD-A5C1', 
  'TCGA-FD-A6TD', 'TCGA-LC-A66R', 'TCGA-K4-A5RI', 'TCGA-FD-A3SP', 'TCGA-GU-A764', 'TCGA-E7-A3X6', 'TCGA-FT-A3EE', 'TCGA-DK-A6B5', 'TCGA-DK-A2I2', 
  'TCGA-FD-A6TK', 'TCGA-DK-A1AE', 'TCGA-FD-A3N6', 'TCGA-FD-A5BT', 'TCGA-DK-AA6Q', 'TCGA-4Z-AA82', 'TCGA-BT-A2LD', 'TCGA-GD-A3OQ', 'TCGA-GU-A762', 
  'TCGA-CU-A72E', 'TCGA-ZF-AA56', 'TCGA-BL-A13I', 'TCGA-UY-A8OC', 'TCGA-XF-A9SX', 'TCGA-XF-AAMT', 'TCGA-BT-A20V', 'TCGA-G2-A3IB', 'TCGA-FD-A62N', 
  'TCGA-HQ-A5ND', 'TCGA-C4-A0F0', 'TCGA-GV-A3QG', 'TCGA-DK-A3IM', 'TCGA-BT-A20U', 'TCGA-BT-A20X', 'TCGA-PQ-A6FN', 'TCGA-FD-A3B7', 'TCGA-FD-A5BU', 
  'TCGA-FD-A62S', 'TCGA-DK-AA6M', 'TCGA-FD-A3B4', 'TCGA-ZF-AA5H', 'TCGA-FD-A43Y', 'TCGA-ZF-AA5N', 'TCGA-FD-A5BY', 'TCGA-DK-A3WX', 'TCGA-XF-A9T4', 
  'TCGA-XF-AAN4', 'TCGA-XF-AAMW', 'TCGA-FD-A3B8', 'TCGA-FD-A5BS', 'TCGA-XF-AAME', 'TCGA-DK-A3WY', 'TCGA-XF-AAN8');


tcga_sc_basal_id = intersect(tcga_sc_id, tcga_basal_sc_id)
tcga_sc_basal_id

save(tcga_sc_basal_id, tcga_sc_id, tcga_basal_id, file="tcga_sc_basal_id.RData")
getwd()


length(tcga_sc_id)
length(tcga_basal_id)
length(tcga_sc_basal_id)

## RNA Seq
load("../bcg/tcga_blca.RData")

tcga_blca_8788_408 = blca.cc.fil2.log.z[, grep("01A", substr(colnames(blca.cc.fil2.log.z), 14, 16))]
dim(tcga_blca_8788_408)
head(tcga_blca_8788_408)


#colnames(tcga_blca_8788_408) = substr(colnames(tcga_blca_8788_408), 1, 12)
tcga_blca_8788_408_gs = blca.rows[row.names(tcga_blca_8788_408), "external_gene_name"]
tcga_blca_8788_408[1:4,1:4]

annotation_col = data.frame(row.names = colnames(tcga_blca_8788_408), 
                                type=rep("UC", ncol(tcga_blca_8788_408)), stringsAsFactors = F)
annotation_col$type[substr(colnames(tcga_blca_8788_408), 1, 12) %in% tcga_sc_id] = "UC_w_SC"
annotation_col$type[substr(colnames(tcga_blca_8788_408), 1, 12) %in% tcga_basal_id] = "basal"
annotation_col$type[substr(colnames(tcga_blca_8788_408), 1, 12) %in% tcga_sc_basal_id] = "SC_basal"
head(annotation_col)

annotation_row = data.frame(row.names = row.names(tcga_blca_8788_408), SC_vs_UC=rep("NoChange", nrow(tcga_blca_8788_408)), stringsAsFactors = F)

up = row.names(t2_t1_sum_df[!is.na(t2_t1_sum_df$pvalues) & t2_t1_sum_df$pvalues < 0.01 & !is.na(t2_t1_sum_df$log2FC) & t2_t1_sum_df$log2FC > 1,])
annotation_row$SC_vs_UC[tcga_blca_8788_408_gs %in% up] = "up"


dn = row.names(t2_t1_sum_df[!is.na(t2_t1_sum_df$pvalues) & t2_t1_sum_df$pvalues < 0.01 & !is.na(t2_t1_sum_df$log2FC) & t2_t1_sum_df$log2FC < -1,])
annotation_row$SC_vs_UC[tcga_blca_8788_408_gs %in% dn] = "dn"

head(annotation_row)


table(annotation_col$type)

ann_colors = list(
    type = c(UC = "#cde6ec", UC_w_SC = "#ad7f75", basal = "#A020F0", SC_basal = "#072591"),
    SC_vs_UC = c(NoChange = "#7570B3", up = "#E7298A", dn = "#66A61E")
)

library(gplots)

breaks = c(
    seq(from = -3,  to = 0,  length.out = 8),
    seq(from = 0,   to = 3,  length.out = 8)
          )
breaks = unique(breaks)
col = greenred(length(breaks + 1))
breaks
col

options(repr.plot.width=15, repr.plot.height= 15)
pheatmap(tcga_blca_8788_408, color=col, breaks = breaks, scale='none', 
        main='408 patients, 8788 genes (SD >1.5, mean reads > 5)', 
        clustering_method = 'complete', cluster_rows = TRUE, cluster_cols = T, 
        show_rownames = F, show_colnames = F, annotation_col = annotation_col,
        annotation_row = annotation_row, annotation_colors = ann_colors
         )


dim(tcga_updn)

tcga_updn = tcga_blca_8788_408[annotation_row$SC_vs_UC != "NoChange",  ]
row.names(tcga_updn) = blca.rows[row.names(tcga_updn), "external_gene_name"]

## to change the ENSEMBL ID to gene symbol
annotation_row_updn = data.frame(row.names = row.names(tcga_updn), SC_vs_UC=rep("NoChange", nrow(tcga_updn)), stringsAsFactors = F)

up = t2_t1_sum_df$symbol[!is.na(t2_t1_sum_df$pvalues) & t2_t1_sum_df$pvalues < 0.01 & !is.na(t2_t1_sum_df$log2FC) & t2_t1_sum_df$log2FC > 1]
annotation_row_updn$SC_vs_UC[row.names(tcga_updn) %in% up] = "up"

dn = t2_t1_sum_df$symbol[!is.na(t2_t1_sum_df$pvalues) & t2_t1_sum_df$pvalues < 0.01 & !is.na(t2_t1_sum_df$log2FC) & t2_t1_sum_df$log2FC < -1]
annotation_row_updn$SC_vs_UC[row.names(tcga_updn) %in% dn] = "dn"

tcga_updn[tcga_updn > 3] = 3
tcga_updn[tcga_updn < -3] = -3
options(repr.plot.width=10, repr.plot.height= 9)
pheatmap(tcga_updn, color=col, breaks = breaks, scale='none', 
        main='408 patients', 
        clustering_method = 'ward.D', cluster_rows = TRUE, cluster_cols = T, 
        show_rownames = F, show_colnames = F, annotation_col = annotation_col,
        annotation_row = annotation_row_updn, annotation_colors = ann_colors
         )


tcga_updn = tcga_blca_8788_408[annotation_row$SC_vs_UC != "NoChange",  ]
row.names(tcga_updn) = blca.rows[row.names(tcga_updn), "external_gene_name"]

## to change the ENSEMBL ID to gene symbol
annotation_row_updn = data.frame(row.names = row.names(tcga_updn), SC_vs_UC=rep("NoChange", nrow(tcga_updn)), stringsAsFactors = F)

up = t2_t1_sum_df$symbol[!is.na(t2_t1_sum_df$pvalues) & t2_t1_sum_df$pvalues < 0.01 & !is.na(t2_t1_sum_df$log2FC) & t2_t1_sum_df$log2FC > 1]
annotation_row_updn$SC_vs_UC[row.names(tcga_updn) %in% up] = "up"

dn = t2_t1_sum_df$symbol[!is.na(t2_t1_sum_df$pvalues) & t2_t1_sum_df$pvalues < 0.01 & !is.na(t2_t1_sum_df$log2FC) & t2_t1_sum_df$log2FC < -1]
annotation_row_updn$SC_vs_UC[row.names(tcga_updn) %in% dn] = "dn"

tcga_updn[tcga_updn > 3] = 3
tcga_updn[tcga_updn < -3] = -3
options(repr.plot.width=15, repr.plot.height= 20)
pheatmap(tcga_updn, color=col, breaks = breaks, scale='none', 
        main='408 patients', 
        clustering_method = 'ward.D', cluster_rows = TRUE, cluster_cols = T, 
        show_rownames = T, show_colnames = F, annotation_col = annotation_col,
        annotation_row = annotation_row_updn, annotation_colors = ann_colors
         )


table(annotation_col$type)

annotation_col_updn = annotation_col[annotation_col$type %in% c('basal', 'UC_w_SC', "SC_basal"), , drop = F]

head(annotation_col_updn)

tcga_updn = tcga_blca_8788_408[annotation_row$SC_vs_UC != "NoChange",  ]
row.names(tcga_updn) = blca.rows[row.names(tcga_updn), "external_gene_name"]

## to change the ENSEMBL ID to gene symbol
annotation_row_updn = data.frame(row.names = row.names(tcga_updn), SC_vs_UC=rep("NoChange", nrow(tcga_updn)), stringsAsFactors = F)
  
up = t2_t1_sum_df$symbol[!is.na(t2_t1_sum_df$pvalues) & t2_t1_sum_df$pvalues < 0.01 & !is.na(t2_t1_sum_df$log2FC) & t2_t1_sum_df$log2FC > 1]
annotation_row_updn$SC_vs_UC[row.names(pdn) %in% up] = "up"

dn = t2_t1_sum_df$symbol[!is.na(t2_t1_sum_df$pvalues) & t2_t1_sum_df$pvalues < 0.01 & !is.na(t2_t1_sum_df$log2FC) & t2_t1_sum_df$log2FC < -1]
annotation_row_updn$SC_vs_UC[row.names(tcga_updn) %in% dn] = "dn"

annotation_col_updn = annotation_col[annotation_col$type %in% c('basal', 'UC_w_SC', "SC_basal"), , drop = F]

tcga_basal_updn = tcga_updn[,annotation_col$type %in% c('basal', 'UC_w_SC', "SC_basal")]

tcga_basal_updn[tcga_basal_updn > 3] = 3
tcga_basal_updn[tcga_basal_updn < -3] = -3

options(repr.plot.width=15, repr.plot.height= 10)
pheatmap(tcga_basal_updn, color=col, breaks = breaks, scale='none', 
        main='', 
        clustering_method = 'ward.D', cluster_rows = TRUE, cluster_cols = T, 
        show_rownames = F, show_colnames = F, annotation_col = annotation_col_updn,
        annotation_row = annotation_row_updn, annotation_colors = ann_colors
         )


t2_t1_sum_df['KRT1',]

tmp = t12['KRT1', ]
plot_dot(tmp[1:10, 11:20, genesymbol = "KRT1", pval = 0.00037])


condition = rep("UC", times = ncol(blca.assay));
condition[substr(colnames(blca.assay), 1, 12) %in% tcga_sc_id] = "UC_w_SC"

design = data.frame(
                row.names       = colnames(blca.assay),
                condition       = condition,
                libType = rep("PE", ncol(blca.assay)));

tcga_ddsmat = DESeqDataSetFromMatrix(countData = blca.assay,
                colData = design,
                design = ~ condition); 


tcga_dds_ds <- estimateSizeFactors(tcga_ddsmat);


#rld = rlog(dds)

#save.image()

tcga_dds <- DESeq(tcga_dds_ds, parallel=T);
tcga_res = results(tcga_dds, contrasts = c("condition", "UC_w_SC", "UC"), nrow = nrow(tcga_dds))


tcga_res$symbol = blca.rows[row.names(tcga_res), "external_gene_name"]
head(tcga_res)

dim(tcga_res[!is.na(tcga_res$padj) & tcga_res$padj < 0.01 & !is.na(tcga_res$log2FoldChange) & abs(tcga_res$log2FoldChange) > 1, ])

sel = !is.na(tcga_res$padj) & tcga_res$padj < 0.01 & !is.na(tcga_res$log2FoldChange) & abs(tcga_res$log2FoldChange) > 1.5
tmp = blca.assay[row.names(tcga_res[sel,]), ]
dim(tmp        )

anno_col = data.frame(row.names=colnames(blca.assay), type=rep("UC", ncol(blca.assay)), stringsAsFactors = F)
anno_col$type[substr(colnames(blca.assay), 1, 12) %in% tcga_sc_id] = "UC_w_SC"
head(anno_col)



plot_dot = function(var1, var2, genesymbol, title){
    xx_base = c(1, 1.3)
    xx = rep(xx_base, each = length(var1))
    yy1 = var1
    yy2 = var2
    yy = c(yy1, yy2)
    ##
    ## ylim
    y1 = floor(min(yy)); 
    y2 = ceiling(max(yy)); 
    if(y1 > 0) { y1 = 0}
    ##
    plot("", type="n", xlab="", xlim=c(0.8, 1.7), ylim=c(y1, y2), axes=F, ylab="log2(FPKM)", 
         main=title)
    axis(2)
    text(xx_base, y1, c("T1", "T2"), adj=c(.5, 2), xpd=T)
    xx_rand = xx + runif(n = length(xx), min = -0.05, max = 0.05)
    points(xx_rand, yy, cex=1.5)
    segments(xx_rand[1:length(var1)], yy[1:length(var1)], 
             xx_rand[(length(var1)+1):(length(var1)*2)], yy[(length(var1)+1):(length(var1)*2)])
    shift = c(-0.1, 0.1)
    # mean
    segments(xx_base - 0.07 + shift, c(mean(yy1), mean(yy2)), 
             xx_base + 0.07 + shift, c(mean(yy1), mean(yy2)), col=1, xpd=T, lwd=2) 
    # top sd
    segments(xx_base - 0.05 + shift, c(mean(yy1) + sd(yy1), mean(yy2) + sd(yy2)), 
             xx_base + 0.05 + shift, c(mean(yy1) + sd(yy1), mean(yy2) + sd(yy2)), col=1, xpd=T) 
    # bot sd
    segments(xx_base - 0.05 + shift, c(mean(yy1) - sd(yy1), mean(yy2) - sd(yy2)), 
             xx_base + 0.05 + shift, c(mean(yy1) - sd(yy1), mean(yy2) - sd(yy2)), col=1, xpd=T) 
    # vertical line
    segments(xx_base + shift, c(mean(yy1) - sd(yy1), mean(yy2) - sd(yy2)), 
                 xx_base + shift, c(mean(yy1) + sd(yy1), mean(yy2) + sd(yy2)), col=1, xpd=T)
}

ls()

##### for grant 12/21/2017
library("DESeq2")
cn = colnames(rsem.counts)
cn
ids = c(200, 202, 204, 206, 207, 208, 210, 211)
sel1 = paste0('DS_bla_', ids, '_T1'); sel1
sel2 = paste0('DS_bla_', ids, '_T2'); sel2
sel = c(sel1, sel2)
rsem.sel = rsem.counts.fil[, sel]
head(rsem.sel)

sel
condition = sub(".*_", "", sel)
patientID = sub("_..$", "", sel)
condition
patientID
design = data.frame(
		    row.names       = colnames(rsem.sel),
		    condition       = condition,
		    patientID 	    = patientID,
		    libType         = rep("PE", ncol(rsem.sel)));
design$patientID

ddsmat = DESeqDataSetFromMatrix(countData = rsem.sel,
				colData = design,
				design = ~ patientID + condition);
dds.ds <- estimateSizeFactors(ddsmat);
dds <- DESeq(dds.ds, parallel=T);

res5 = results(dds)
res5.cc = counts(dds, normalized = T)
res5.cc = as.data.frame(res5.cc)

class(res5.cc)

library(future)
 
#rld <- rlog(dds, blind = FALSE)
f <- future({
	rld <- rlog(dds, blind = FALSE)
	rld
}) %plan% multiprocess

res5.rld  <- value(f)

res5 = res5[order(res5$pvalue),]
res5$symbol = sub(".*_", "", row.names(res5))
res5.sig = res5[!is.na(res5$log2FoldChange) & abs(res5$log2FoldChange) > 0.5 & !is.na(res5$padj) & res5$padj < 0.05,]
dim(res5.sig)
head(res5.sig)
head(res5)

source("~/program/fun/write_rnk.r")
write_rnk(res5, file='res5/res5.rnk')
write.csv(res5, file='res5/res5.csv')

source("~/program/fun/run_gsea.R")
run_gsea('res5.rnk')

anno.col = data.frame(row.names = colnames(res5.cc), sampletype = condition, patient = patientID)
anno.col

tmp = res5.cc
tmp[tmp < 1] = 1
res5.cc.log2 = log2(tmp)
rm(tmp)

lls = seq(0.9, 3, by=0.2)
res5.sig2 = res5.sig[res5.sig$baseMean > 50,]
ll = 1.9
res5.sig.id = row.names(res5.sig2)
matrix_tmp = res5.cc.log2[res5.sig.id, ]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > ll] = ll
matrix_tmp[matrix_tmp < -ll] = ll 
row.names(matrix_tmp) = sub(".*_", "", row.names(matrix_tmp))
matrix_tmp
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
pdf(paste0("res5/heatmap_sig_genes_", ll, ".pdf"), width=5, height=5)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), show_colnames = T,
	      scale='none', show_rownames = F, breaks = breaks, 
	      main='Different genes in paired T2 vs T1', clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
	      annotation_col =anno.col)
dev.off()

library(ggplot2)
library(data.table)

res5.cc.2 = as.data.table(res5.cc)
res5.cc.2[, symbol := sub(".*_", "", row.names(res5.cc))]
head(res5.cc.2)
res5.cc.3 = melt.data.table(res5.cc.2, id.vars = c('symbol'), measure.vars=colnames(res5.cc))
res5.cc.3[, group := substr(variable, 12, 13)]
res5.cc.3[, patient := substr(variable, 1, 10)]
head(res5.cc.3)
res5.cc.3[ value < 1, value := 1]

as.data.frame(res5[res5$symbol %in% immuno_target$symbol, ])
as.data.frame(res.inc[res.inc$Gene %in% immuno_target$symbol, ])
as.data.frame(res.exc[res.exc$Gene %in% immuno_target$symbol, ])

for (i in 1:nrow(res5.sig)){
for (i in 1:nrow(immuno_target)){
	gene = immuno_target$symbol[i]
	fdr = res5$padj[res5$symbol == gene]
	logfc= res5$log2FoldChange[res5$symbol == gene]
	if(is.na(gene)){next}
	gene
	tmp = res5.cc.3[symbol == gene, ]
	tmp[, log_reads := log2(value)]
	tmp
	tmp.w = dcast(tmp, patient ~ group, value.var='log_reads'); tmp.w
	mbar = apply(tmp.w[,2:3], 2, median, na.rm=T) 
	mbar = as.data.table(mbar)
	g = ggplot(tmp, aes(x = group, y =log_reads))
	g = g + geom_jitter(width=0.2)
	g = g + geom_segment(data=mbar, aes(x = c(.7, 1.7), y = mbar, xend = c(1.3, 2.3), yend = mbar), alpha=.5, lwd=1.5, colour='blue')
	g = g + geom_segment(data=tmp.w, aes(x = 1.2, y = T1, xend = 1.8, yend = T2), alpha=0.2, colour='red')
	g = g + ylab('log2 # reads') + ggtitle(gene) + xlab('')
	g = g + theme(plot.title = element_text(face="bold", size=10, hjust=.5))
	ggsave(g, file=paste0('rnaseq/res5/dotplot/dots_immune_', gene, '_', logfc, '_', fdr, '.pdf'), width=2.5, height=3)
}


gene = 'FOXA1'
fdr = round(res5$padj[res5$symbol == gene],3)
fdr
if(is.na(gene)){next}
tmp = res5.cc.3[symbol == gene, ]
tmp[, log_reads := log2(value)]
tmp
tmp.w = dcast(tmp, patient ~ group, value.var='log_reads'); tmp.w
mbar = apply(tmp.w[,2:3], 2, mean, na.rm=T) 
mbar = as.data.table(mbar)
g = ggplot(tmp, aes(x = group, y =log_reads))
g = g + geom_jitter(width=0.2)
g = g + geom_segment(data=mbar, aes(x = c(.7, 1.7), y = mbar, xend = c(1.3, 2.3), yend = mbar), alpha=.5, lwd=1.5, colour='blue')
g = g + geom_segment(data=tmp.w, aes(x = 1.2, y = T1, xend = 1.8, yend = T2), alpha=0.2, colour='red')
g = g + ylab('log2 # reads') + ggtitle(gene) + xlab('')
g = g + theme(plot.title = element_text(face="bold", size=10, hjust=.5))
ggsave(g, file=paste0('res5/dotplot/dots_', gene, '_', res5.sig$log2FoldChange[i], '_', res5.sig$padj[i], '.pdf'), width=2.5, height=3)

## top 1000 variant genes, pca
pdf("res5/rld_pca.pdf", width=8, height=8)
plotPCA(res5.rld)
dev.off()

library(ggfortify)
row.sd = apply(res5.cc.log2, 1, sd)
row.sd = row.sd[order(row.sd, decreasing=T)]
names(row.sd[1:1000]) -> sel.id
matrix_tmp = res5.cc.log2[row.names(res5.cc.log2) %in% sel.id, ]
pca1 = prcomp(t(matrix_tmp), scale. = TRUE)
pca1$x -> pca.x
pca.x[1:10,1:10]
summary(pca1) -> pca.sum
pca.sum.key = pca.sum$importance
rm(pca.sum, pca1)
pca.sum.key
xlab=paste0('PC1 ', round(100*pca.sum.key[2, 1], 0), '%') 
ylab=paste0('PC2 ', round(100*pca.sum.key[2, 2], 0), '%') 
xlab
ylab
t12 = sub(".*_", "", rownames(pca.x)); t12
t12 = factor(t12)
t12
batch = sub("DS_bla_", "", rownames(pca.x))
batch = sub("_.*", "", batch); batch
as.numeric(batch) -> batch
batch
batch[batch < 200] = 1
batch[batch >= 200] = 2
batch = paste0('batch', batch)
batch = factor(batch)
batch

pca.x = as.data.frame(pca.x)
pca.x$t12 = t12
pca.x$batch = batch
pca.x
rm(t12, batch)

pdf("res5/ggplot_pca.pdf", width=8, height=8)
ggplot(data=pca.x, aes(x=PC1, y = PC2, shape=batch, fill=t12, color=t12)) +
	geom_hline(yintercept = 0, colour = "gray65") +
	geom_vline(xintercept = 0, colour = "gray65") +
	geom_point(size = 4) +
	ggtitle("PCA analysis of RNA-Seq") +
	xlab(xlab) +
	ylab(ylab)
dev.off()

pdf("res5/topvar1000_pca.pdf", width=8, height=8)
plotPCA(rld)
dev.off()


## top 1000 variant genes heatmap
row.sd = apply(res5.cc.log2, 1, sd)
row.sd = row.sd[order(row.sd, decreasing=T)]
names(row.sd[1:1000]) -> sel.id
matrix_tmp = res5.cc.log2[row.names(res5.cc.log2) %in% sel.id, ]
matrix_tmp
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
ll = 2
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > ll] = ll
matrix_tmp[matrix_tmp < -ll] = ll 
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
pdf(paste0("res5/heatmap_topvar1000_genes_", ll, ".pdf"), width=8, height=8)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), show_colnames = T,
	      scale='none', show_rownames = F, breaks = breaks, 
	      main='top varied 1000 genes', clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
	      annotation_col =anno.col)
dev.off()

## top 100 bladder genes 
gtex_bladder_top100 = read.csv("../gtex_bladder_top100.txt")
gtex_bladder_top100
gs = sub(".*_", "", row.names(res5.cc.log2))
matrix_tmp = res5.cc.log2[gs %in% gtex_bladder_top100$Gene.Symbol, ]
matrix_tmp
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > ll] = ll
matrix_tmp[matrix_tmp < -ll] = ll 
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
pdf(paste0("res5/heatmap_bladder100_genes_", ll, ".pdf"), width=8, height=8)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), show_colnames = T,
	      scale='none', show_rownames = F, breaks = breaks, 
	      main='Different genes in paired T2 vs T1', clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
	      annotation_col =anno.col)
dev.off()

##immune gene target heatmap
# 'C10orf54' -> VSIR
immuno_target_symbol = fread("immuno_gene.tsv", header=F)
immuno_target_symbol = immuno_target_symbol$V1 
immuno_target_symbol

#immuno_target_symbol = c('PDCD1', 'CD274', 'PDCD1LG2', 'HAVCR2', 'LAG3', 'BTLA', 'VSIR','CD8A', 'CD4', 'PTPRC', 'MKI67', 'FOXP3', 'CD68', 'STAT1', 'GZMB')
#immuno_target_alias = c('PD-1', 'PD-L1', 'PD-L2', 'TIM-3', 'LAG-3', 'BTLA', 'VISTA', 'CD8', 'CD4', 'CD45', 'Ki67', 'Fox-P3', 'CD68', 'pSTAT1', 'granzyme B')
#immuno_target = data.frame(row.names = immuno_target_symbol, symbol = immuno_target_symbol, alias = immuno_target_alias, stringsAsFactors = F)
#immuno_target

res5.cc.symbol = toupper(sub(".*_", "", row.names(res5.cc.log2)))
res5.cc.log2[res5.cc.symbol %in% immuno_target_symbol, ] -> res5.cc.sel
row.names(res5.cc.sel) = sub(".*_", "", row.names(res5.cc.sel))
res5.cc.sel
matrix_tmp = res5.cc.sel
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
ll = 2
matrix_tmp[matrix_tmp > ll] = ll
matrix_tmp[matrix_tmp < -ll] = ll 
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
pdf(paste0("res5/heatmap_immuno_targets.pdf"), width=7, height=6.5)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), show_colnames = T,
	      scale='none', show_rownames = T, breaks = breaks, 
	      main='immune target gene expressions between T1/T2', clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols =T, fontsize_row = 8,
	      annotation_col =anno.col)
dev.off()

##basal luminal classifier genes
load("../../bcg/uc_basal_luminal_marks.RData")

incl.id.v2 = gsub("-", "_", incl.id); incl.id.v2
gs = sub(".*_", "", row.names(res.11pair.cc))
res.11pair.cc[gs %in% mcconkey, incl.id.v2] -> res.11pair.cc.sel
row.names(res.11pair.cc.sel) = sub(".*_", "", row.names(res.11pair.cc.sel))
res.11pair.cc.sel
mtx = log2(res.11pair.cc.sel+1)
cn = colnames(mtx)
mtx = t(apply(mtx, 1, scale))
colnames(mtx) = cn
ll = 2
mtx[mtx > ll] = ll
mtx[mtx < -ll] = ll 
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks);breaks
anno_rows = data.frame(row.names = row.names(mtx), Class = rep("basal", nrow(mtx)), stringsAsFactors = F)
anno_rows$Class[row.names(anno_rows) %in% mcconkey.luminal] = "luminal"
anno_rows
anno_cols = data.frame(row.names=colnames(mtx), CellType = rep("T1", ncol(mtx)) , stringsAsFactors = F)
anno_cols$CellType[grep("T2", row.names(anno_cols))] = "T2"
anno_cols$CellType[grep("_M", row.names(anno_cols))] = "M"
pdf(paste0("rnaseq/res/heatmap_12pair_mdAnderson_colF.pdf"), width=6, height=5)
hp = pheatmap(mtx, color=greenred(length(breaks)+1), show_colnames = F,
	      scale='none', show_rownames = T, breaks = breaks, 
	      main='Luminal/basal gene expression', clustering_method = 'ward.D2', 
	      cluster_rows = TRUE, cluster_cols = F, fontsize_row = 7, fontsize_col = 7,
	      annotation_col =anno_cols, annotation_row = anno_rows)
dev.off()
sync('rnaseq/res')

## exome seq for grant
setwd('..')
library(data.table)
library(autospy)
library(WriteXLS)
library(ggplot2)
library(metafolio)
library(VennDiagram)

options(width=155)

bcr = paste0('s_', sel)
bcr

scc.maf[Tumor_Sample_Barcode %in% bcr,] -> maf.g

## /ifs/res/share/solit/alahmadh/Proj_07813_DF
pid = "bla_[0-9]+_"
sid = "_[T|M][0-9]?$"

tmp = as.data.table(table(scc.maf.fil$Tumor_Sample_Barcode))
tmp
tmp[, patient := substr(V1, 3, 12)]
tmp[, group := substr(V1, 14, 15)]
tmp = tmp[grep("M", group, invert=T), ]
tmp
t.test(tmp$N ~ tmp$group)
tmp.w = dcast(tmp, patient ~ group, value.var='N'); tmp.w
tmp.w = as.data.table(tmp.w)
tmp.w
wilcox.test(tmp.w$T1, tmp.w$T2, paired=T)
mbar = apply(tmp.w[,2:3], 2, median, na.rm=T) 
mbar = as.data.table(mbar)
g = ggplot(tmp, aes(x = group, y = N))
g = g + geom_jitter(width=0.2)
g = g + scale_y_log10()
g = g + geom_segment(data=mbar, aes(x = c(.7, 1.7), y = mbar, xend = c(1.3, 2.3), yend = mbar), alpha=.5, lwd=1.5, colour='blue')
g = g + geom_segment(data=tmp.w, aes(x = 1.2, y = T1, xend = 1.8, yend = T2), alpha=0.2, colour='red')
g = g + ylab('# of mutations') + ggtitle('') + xlab('')
g = g + theme(plot.title = element_text(face="bold", size=10, hjust=.5))
ggsave(g, file=paste0('exonseq/maftools/dotplot_mutation_N.pdf'), width=2.5, height=3)
sync('exonseq/maftools')

#fread(maffile) -> bak

IMPACT468 <- scan("/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv6/genelist", what = "")
scc.maf[, IMPACT_468 := Hugo_Symbol %in% IMPACT468]

## gene level facets copy number results
cn.table = fread('exonseq/autospy/genelevel_facets_cnv.txt')
cn.table = cn.table[, .(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Type)]
fwrite(cn.table, file='cntable.txt', sep="\t", quote=F)

## clinical data
clin = data.table(Tumor_Sample_Barcode = unique(maf.g[,Tumor_Sample_Barcode]))
clin[, CellType := sub(".*_", "", Tumor_Sample_Barcode)]
clin
fwrite(clin, file='clinical.txt', sep="\t", quote=F)

maf.g[, id := paste(Variant_Classification, oncogenic, sep=" ")]
setkey(maf.g, 'Variant_Classification')
maf.g = maf.g[simple_table,] 
table(maf.g$simple_value)
## reassign oncogenic value
setkey(maf.g, 'oncogenic')
maf.g = maf.g[oncogenic_table, ]
table(maf.g$oncogenic_value)
## assign color
maf.g[, id := paste(simple_value, oncogenic_value, sep=' ')]
setkey(maf.g, 'id')
maf.g = maf.g[color_table,] 
table(maf.g$color_value)

maf.g[, Variant_Classification_old := Variant_Classification]
maf.g[, Variant_Classification := id]
maf.g[t_alt_count > 3 & t_alt_count / t_depth >0.05,] -> maf.g.fil
maf.g.fil

fwrite(maf.g.fil, file='maf_g_fil.txt', sep="\t", quote=F)

read.maf('maf_g_fil.txt', clinicalData='clinical.txt', vc_nonSyn = vc.key) -> maf.g.o

scc.maf.fil[Hugo_Symbol == 'KMT2D', .(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Variant_Classification_old, oncogenic)]
table(scc.maf.o@data$Variant_Classification)

# most amplied known function
mut.2017 = c('TP53', 'RB1', 'RHOB', 'PIK3CA', 'KDM6A', 'TSC1', 'ELF3', 'KMT2D', 'CREBBP', 'CDKN1A', 'EP300', 'ZFP36L1', 'ARID1A', 'STAG2', 'CDKN2A', 'HRAS', 'KRAS', 'FBXW7', 'ERCC2', 'ASXL2', 'RHOA', 'KMT2A', 'FGFR3', 'NFE2L2', 'KMT2C', 'PSIP1', 'KANSL1', 'C3orf70', 'FAT1', 'SPTAN1', 'RXRA', 'ZBTB7B', 'PTEN', 'ATM', 'KLF5', 'PARD3', 'CUL1', 'NRAS', 'SF3B1', 'GNA13', 'RBM10', 'ACTB', 'MBD1', 'CASP8', 'HIST1H3B', 'TAF11', 'ERBB2', 'NUP93', 'SF1', 'ERBB3', 'METTL3', 'SPN', 'MB21D2', 'SSH3', 'USP28', 'ASXL1', 'TMCO4', 'HES1', 'ZNF773')
mut.2014 = c('TP53', 'MLL2', 'ARID1A', 'KDM6A', 'PIK3CA', 'EP300', 'CDKN1A', 'RB1', 'ERCC2', 'FGFR3', 'STAG2', 'ERBB3', 'FBXW7', 'RXRA', 'ELF3', 'NFE2L2', 'TSC1', 'KLF5', 'TXNIP', 'FOXQ1', 'CDKN2A', 'RHOB', 'FOXA1', 'PAIP1', 'BTG2', 'HRAS', 'ZFP36L', 'RHOA', 'CCND3')
amp.2013 = c('AHR', 'BCL2L1', 'CCND1', 'CCNE1', 'E2F3', 'EGFR', 'ERBB2', 'FGFR3', 'GATA3', 'KRAS', 'MDM2', 'MYCL1', 'PPARG', 'PVRL4', 'SOX4', 'TERT', 'YWHAZ', 'ZNF703') 

maf.g[t_alt_count > 3 & t_alt_count / t_depth >0.05,] -> maf.g.fil
maf.g.o = read.maf(maf.g.fil)
plotmafSummary(maf.g.o, rmOutlier=T, addStat = 'median', dashboard = T, 
	       file='exonseq/maftools5/maftools_summary_scc', width=8, height=6)

mut.2017.sel = c('KMT2D', 'TP53', 'PIK3CA', 'ATM', 'ARID1A', 'EP300', 'FGFR3', 'STAG2', 'KDM6A', 'KMT2C', 
		 'RB1', 'CUL1', 'ELF3', 'ERBB3', 'KMT2A', 'NFE2L2', 'RHOA', 'ACTB', 'CREBBP', 'SPTAN1')
anno.df = data.frame(Tumor_Sample_Barcode = unique(maf.g.fil$Tumor_Sample_Barcode))
anno.df$CellType = sub(".*_", "", anno.df$Tumor_Sample_Barcode)
anno.df

## the data set for maf summary and oncoplot is different
w = 7; h = 8 
pdf('exonseq/maftools5/oncoplot_mut2017.pdf', width=w, height=h)
oncoplot(maf = maf.g.o, genes=mut.2017.sel, color = color_vector, 
	 annotationDat = anno.df, showTumorSampleBarcodes=T, top=20, clinicalFeatures='CellType')
dev.off()

maf.g.t1
maf.g.t2
maf.g.t1 = maf.g[grep("T1", Tumor_Sample_Barcode),]
maf.g.t2 = maf.g[grep("T2", Tumor_Sample_Barcode),]
maf.g.t1.o = read.maf(maf.g.t1, vc_nonSyn = vc.key)
maf.g.t2.o = read.maf(maf.g.t2, vc_nonSyn = vc.key)

plotmafSummary(maf.g.t1.o, rmOutlier=T, addStat = 'median', dashboard = T, file='exonseq/maftools5/maftools_summary_scc_t1.pdf')
plotmafSummary(maf.g.t2.o, rmOutlier=T, addStat = 'median', dashboard = T, file='exonseq/maftools5/maftools_summary_scc_t2.pdf')

pdf('exonseq/maftools5/oncoplot_t1.pdf', width=w, height=h)
oncoplot(maf = maf.g.t1.o, genes=mut.2017.sel, color = color_vector, showTumorSampleBarcodes=T, top=20, clinicalFeatures='Cell_Type', annotationDat=clin)
dev.off()
pdf('exonseq/maftools5/oncoplot_t2.pdf', width=w, height=h)
oncoplot(maf = maf.g.t2.o, genes=mut.2017.sel, color = color_vector, showTumorSampleBarcodes=T, top=20, clinicalFeatures='Cell_Type', annotationDat=clin)
dev.off()


## for autospy
#fread(maffile) -> maf.g
IMPACT468 <- scan("/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv6/genelist", what = "")
maf.g.fil[, IMPACT_468 := Hugo_Symbol %in% IMPACT468]
maf.g[, Variant_Bioportal := paste0(Variant_Classification, oncogenic)]

maf.g.fil[, patient := stringr::str_extract(Tumor_Sample_Barcode, pid)]
maf.g.fil[, sample := stringr::str_extract(Tumor_Sample_Barcode, sid)]
dim(maf.g.fil)

maf.g.fil[, t_var_freq := as.numeric(t_alt_count) / t_depth]
maf.g.fil[, tm := paste0(Hugo_Symbol, " ", as.numeric(gsub("[^\\d]+", "", HGVSp_Short, perl=T)))] # 
maf.g.fil[, TAG := paste0('chr', Chromosome, ':', Start_Position, '-', End_Position, ':', Reference_Allele, ':', Tumor_Seq_Allele2)]
maf.g.fil[, variant := paste0(TAG, '::', Hugo_Symbol, ':', HGVSp_Short)]

sampleID = data.table(`Sample ID` = unique(maf.g.fil$Tumor_Sample_Barcode))
sampleID[, `Sample Class` := 'Primary']
sampleID[grep("T2", `Sample ID`), `Sample Class` := 'Tumor']
sampleID[grep("M", `Sample ID`), `Sample Class` := 'Tumor']
sampleID[, patient := stringr::str_extract(`Sample ID`, pid)]
sampleID[, sample := stringr::str_extract(`Sample ID`, sid)]
sampleID

o = getwd()
setwd('exonseq/autospy5')

## 
source('~/program/autospy/R/process_autopsy_maf.R')
source('~/program/autospy/R/stratton_plot.R')

### refilter
library(autospy)
filter_results <- filter_maf_report(scc.maf.fil)
filter_results
dim(scc.maf.fil)
maf.g.fil <- filter_results$maf

write.table(maf.g.fil, 'filtered.maf', sep="\t", quote=F, row.names=F)

### overlap plot
make_mutation_overlap_plot(scc.maf.fil, pid = pid, log = F, out = "pre_filter")

## mutation signature
cmd = 'python ~/program/autospy/inst/mutation-signatures/main.py --seed 100 ~/program/autospy/inst/mutation-signatures/Stratton_signatures30_bbn.txt filtered.maf signature.out'
system(cmd)
fread("signature.out") -> mut.sig
mut.sig
ggsave(plot_mutation_signatures(mut.sig, pid, sid, fraction_threshold = 0.15), 
       filename='mutation_singaure_decompoisition_combined.pdf', width=25, height=6)
sync('exonseq', 'autospy')

### patient wise
patients = unique(maf.g.fil$patient)
tumor_samples = unique(maf.g.fil$Tumor_Sample_Barcode)

for(pat in patients){
	primary <- sampleID[patient == pat & `Sample Class` == 'Primary', `Sample ID`]
	primary
	pat.maf <- scc.maf.fil[patient== pat,] 
	pat.maf
	print(paste("make_variant_classification_plot for patient ", pat))
	make_variant_classification_plot(pat.maf, out = pat)
	print("make_ccf_plots")
	make_ccf_plots(pat.maf, tumor_samples)
	print("make_stratton_plots")
	make_stratton_plots(pat.maf, tumor_samples, out = pat)
	print("make_binary_tree")
	print("make_mutation_signatures")
	write.table(pat.maf, 'pat.maf', sep="\t", quote=F, row.names=F)
	cmd = paste0('python ~/program/autospy/inst/mutation-signatures/main.py --seed 100 ')
	cmd = paste0(cmd, '~/program/autospy/inst/mutation-signatures/Stratton_signatures29.txt pat.maf signature.out')
	system(cmd)
	fread("signature.out") -> mut.sig
	ggsave(plot_mutation_signatures(mut.sig, pid, sid, fraction_threshold = 0.5), 
	       filename=paste0('mutation_singaure_decompoisition_', pat, '.pdf'), width=5, height=6)
	samples = unique(pat.maf$sample)
	sample_pairs <- combn(samples, 2, simplify = F)
	dir.create('ccf_2d_plots')
	for(sample_pairs in sample_pairs){
		make_ccf_2d(pat.maf,
			    sample_pairs,
			    out = pat,
			    directory = 'ccf_2d_plots')
	}
}

for(pat in patients){
	primary <- sampleID[patient == pat & `Sample Class` == 'Primary', `Sample ID`]
	pat.maf <- scc.maf.fil[patient== pat,] 
	make_binary_tree(pat.maf, primary, hotspots = autospy::hotspots, vertical = TRUE, margins = TRUE, out = pat)
}

### maftools
scc.maf.fil.o = read.maf(scc.maf.fil)
plotmafSummary(maf.g.o, rmOutlier=T, addStat = 'median', dashboard = T, file='exonseq/autospy5/maftools_summary_v2.pdf')

patients = unique(maf.g$Tumor_Sample_Barcode)
lapply(unique(patients),
       function(tsb){
	       maftools::rainfallPlot(maf.g.o, tsb = tsb, savePlot = T)
       })

source('~/program/facets-suite/geneLevel.R')
fread('../cncf.list', header=F) -> genelevel
genelevel[grep("hisens", V1),] -> genelevel
genelevel = setNames(genelevel, 'cncfFile')
genelevel[, cncfFile := paste0('../', cncfFile)]
genelevel[, bcr := str_extract(cncfFile, "DS_bla_..._[T|M][0-9]?")]
genelevel[, outfile := paste0('genelevel/', bcr, '.genelevel')]
dir.create('genelevel')
head(genelevel)
genelevel[, {
	cnv = get_gene_level_calls(cncf_files = cncfFile)
	cnv[, Tumor_Sample_Barcode := bcr]
	fwrite(cnv, file=outfile, sep="\t", quote=F, row.names=F)
}, by=1:nrow(genelevel)]

## left cut
cnv.list = lapply(genelevel$outfile, fread)
rbindlist(cnv.list) -> all.cnv
head(all.cnv)
dim(all.cnv)
all.cnv = all.cnv[ !is.na(FACETS_CNA) & FACETS_CNA != 0,]
all.cnv
all.cnv = all.cnv[, Variant_Classification := 'CNV']
all.cnv = all.cnv[, Variant_Type := 'AMP']
all.cnv = all.cnv[FACETS_CNA < 0, Variant_Type := 'DEL']
fwrite(all.cnv, file='genelevel_facets_cnv.txt', row.names=F, sep="\t", quote=F)
colnames(all.cnv)

out = export.oncoprint(maf, cnv, res.cc, res.sig)
dim(out)
fwrite(out, file='output_4_oncoprint.txt', sep="\t", quote=F)


## rerun facets by comparing T2 to T1 directly
genelevel
genelevel[, countFile := paste0(dirname(rdatafile), '/tmp/Proj_07813_DF_indelRealigned_recal_s_', bcr, ".dat")]
file.exists(genelevel$countFile)
genelevel[, patient := substr(bcr, 1, 10)]
genelevel[, group := sub(".*_", "", bcr)]
genelevel.w = dcast(genelevel[, .(patient, group, countFile)], patient~group)
genelevel.w 

genelevel.w[, facetsMergeOutFile  := paste0('facets/count_merged_', patient, '_T2_T1.dat.gz')]
genelevel.w[, facetsMerge.jobname := paste0('facetsMerge.', patient)]
genelevel.w[, facetsMerge.cmd := bsub.head(facetsMerge.jobname, mem=18, cpu=4, We='8:26')]
#genelevel.w[, facetsMerge.cmd := bsub.head(facetsMerge.jobname, mem=18, cpu=4, We='8:26', cwd=cwd)]
genelevel.w[, facetsMerge.cmd := paste0(facetsMerge.cmd, ' "', facets, '/mergeTN.R ', T2, ' ', T1, ' ', facetsMergeOutFile, ' "')]
genelevel.w$facetsMerge.cmd[1]
genelevel.w
file.exists(genelevel.w$T2)
genelevel.w$T2

exe.jobs(genelevel.w$facetsMerge.cmd, logger)

# below is the save as the facets/facets_RUN.sh
pcval = 200;cval = 100
genelevel.w[, facetsRunOutDir  := paste0('facets/', patient)] ; 
for(i in 1:nrow(genelevel.w)){system(paste0('mkdir -p ', genelevel.w$facetsRunOutDir[i]))}
genelevel.w[, facetsRun.jobname := paste0('facetsRun.', patient)]
genelevel.w[, facetsRun.cmd := bsub.head(facetsRun.jobname, mem=3, cpu=3, We='8:26')]
genelevel.w[, facetsRun.cmd := paste0(facetsRun.cmd, ' "', FACETS_SUITE, '/doFacets.R -D facets/', patient, ' -r ', Rlib, ' -t ', patient, ' -f ', facetsMergeOutFile, ' -G T -pc ', pcval, ' -c ', cval, ' "')]
genelevel.w$facetsRun.cmd[2]

write.table(genelevel.w$facetsRun.cmd, file='tt', quote=F, sep="\t", row.names=F, col.names=F)

exe.jobs(genelevel.w$facetsRun.cmd, logger)

## gene level by t2 vs t1
genelevel.w[, cncfFile := paste0('facets/', patient, '/', patient, '_hisens.cncf.txt')]
genelevel.w[,  outfile := paste0('facets/', patient, '_genelevel.txt')]
genelevel.w$outfile
genelevel.w$cncfFile
genelevel.w[, {
	cnv = get_gene_level_calls(cncf_files = cncfFile)
	cnv[, Tumor_Sample_Barcode := bcr]
	fwrite(cnv, file=outfile, sep="\t", quote=F, row.names=F)
}, by=1:nrow(genelevel.w)]

cnv.list = lapply(genelevel.w$outfile, fread)
rbindlist(cnv.list) -> t1t2.cnv
head(t1t2.cnv)
dim(t1t2.cnv)
t1t2.cnv = t1t2.cnv[ !is.na(FACETS_CNA) & FACETS_CNA != 0,]
t1t2.cnv
t1t2.cnv = t1t2.cnv[, Variant_Classification := 'CNV']
t1t2.cnv = t1t2.cnv[, Variant_Type := 'AMP']
t1t2.cnv = t1t2.cnv[FACETS_CNA < 0, Variant_Type := 'DEL']
fwrite(t1t2.cnv, file='genelevel_facets_t1t2_cnv.txt', row.names=F, sep="\t", quote=F)
colnames(t1t2.cnv)
head(all.cnv)

## facets plot gene specific t2 vs t1
genes = c('GATA3', 'FOXA1', 'PPARG')
source('~/program/facets-suite/fPlots_ggplot2.R')
gene = genes[1]
genelevel.w[, rdatafile := paste0('facets/', patient, '/', patient, '_hisens.Rdata')]
genelevel.w[, bcr := patient]
genelevel.w
file.exists(genelevel.w$rdatafile)
for(gene in genes){
	tmp = t1t2.cnv[Hugo_Symbol %in%  gene, .(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Type, FACETS_CNA)];tmp
	for(pat in tmp$Tumor_Sample_Barcode){
		genelevel.w[bcr == pat, rdatafile]
		load(genelevel.w[bcr == pat, rdatafile])
		plot.facets.all.output(out, fit, type='png', main=paste0(patient, '| cval: 100'), plotname=paste0('maftools/', patient, '_t1t2_', gene), gene.name=gene)
	}
}

system("rsync -avur /ifs/work/solitlab/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/facets/ mski1925:/Volumes/LaCie/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/facets/")
system("rsync -avur /ifs/work/solitlab/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/maftools/ mski1925:/Volumes/LaCie/huw/solit/study/hiseq/blca_cmo_06155_2016/exonseq/maftools/")



#merge(all.cnv, scc.maf.fil, by = 'Tumor_Sample_Barcode') -> scc.maf.cnv
## export for oncoprint
## http://www.cbioportal.org/oncoprinter.jsp#
## Sample Gene 
## Alteration Type
## R248W	MISSENSE
## R306fs	TRUNC
## MCN237del	INFRAME
## FUSION	FUSION
## HOMDEL	CNA
## HETLOSS	CNA
## GAIN		CNA
## AMP		CNA
## UP		EXP
## DOWN		EXP
export.oncoprint = function(maf, cnv, res.cc, res.sig){
	maf.sel = maf[, .(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Variant_Type, HGVSp_Short)]
	IMPACT468 <- scan("/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv6/genelist", what = "")
	maf.sel[, IMPACT_468 := Hugo_Symbol %in% IMPACT468]
	maf.sel = maf.sel[IMPACT_468 == T, ]
	#write.table(unique(maf.sel$Variant_Classification), '~/program/fun/oncoprint_table', quote=F, sep="\t")
	fread("~/program/fun/oncoprint_table", header = T) -> tbl
	setkey(tbl, 'Variant_Classification')
	setkey(maf.sel, 'Variant_Classification')
	maf.sel.2 = maf.sel[tbl,]
	maf.sel.2 = maf.sel.2[!is.na(Type),]
	maf.sel.2[, Alteration := sub("p.", "", HGVSp_Short)]
	maf.sel.2$Alteration
	## CNA
	## HOMDEL	CNA
	## HETLOSS	CNA
	## GAIN		CNA
	## AMP		CNA
	cnv = all.cnv
	cnv = cnv[, .(Tumor_Sample_Barcode, Hugo_Symbol, chr, Variant_Classification, Variant_Type, tcn, lcn, FACETS_CNA)]
	cnv = cnv[, Type := 'CNA']
	cnv = cnv[FACETS_CNA == 1, Alteration := 'GAIN']
	cnv = cnv[FACETS_CNA > 1, Alteration := 'AMP']
	cnv = cnv[FACETS_CNA == -1, Alteration := 'HETLOSS']
	cnv = cnv[FACETS_CNA < -1, Alteration := 'HOMDEL']
	cnv = cnv[chr == 'X' & FACETS_CNA == -1, Alteration := 'HOMDEL']
	cnv = cnv[chr == 'X' & FACETS_CNA > 0, Alteration := 'AMP']
	out = rbind(maf.sel.2[, .(Tumor_Sample_Barcode, Hugo_Symbol, Alteration, Type)],cnv[, .(Tumor_Sample_Barcode, Hugo_Symbol, Alteration, Type)])
	res.cc.sig = res.cc[row.names(res.sig),] 
	rs = apply(res.cc.sig, 1, mean)
	res.cc.sig = res.cc.sig[rs > 100,]
	dim(res.cc.sig)
	res.cc.sig.sd = apply(res.cc.sig, 1, sd)
	res.cc.sig.zscore = t(apply(res.cc.sig, 1, scale))
	colnames(res.cc.sig.zscore) = colnames(res.cc.sig)
	res.cc.sig.zscore > 1.96
	res.cc.sig.exp = res.cc.sig.zscore
	res.cc.sig.exp[res.cc.sig.zscore > 1.96] = 'UP'
	res.cc.sig.exp[res.cc.sig.zscore < -1.96] = 'DOWN'
	res.cc.sig.exp = melt(res.cc.sig.exp)
	res.cc.sig.exp = res.cc.sig.exp[res.cc.sig.exp$value %in% c('UP', 'DOWN'), ]
	colnames(res.cc.sig.exp) = c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'Alteration')
	res.cc.sig.exp$Type = 'EXP'
	res.cc.sig.exp$Hugo_Symbol == sub(".*_", "", res.cc.sig.exp$Hugo_Symbol) 
	head(res.cc.sig.exp)
	dim(res.cc.sig.exp)
	rbind(out, res.cc.sig.exp) -> out
	out[, Tumor_Sample_Barcode := sub("s_", "", out$Tumor_Sample_Barcode)]
	out
}

res$ID = row.names(res)
res.dt = as.data.table(res)
res.sel = res.dt[symbol %in% IMPACT468 & !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >1,]
res.sel
res.cc.sel = res.cc[row.names(res.cc) %in% res.sel$ID,]
row.names(res.cc.sel) = sub(".*_", "", row.names(res.cc.sel))
res.cc.sel=res.cc.sel[c(1:5, 8),]
res.cc.sel
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(res.cc.sel, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("res/heatmap_impact_sig_genes.pdf", width=10, height=3.5)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), show_colnames = T,
	      scale='none', show_rownames = T, breaks = breaks, 
	      clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
	      annotation_col =anno.col)
dev.off()

# rnaseq by David's data
## modify the id names of immune target
library(org.Hs.eg.db)
tmp = select(org.Hs.eg.db, keys = immuno_target$symbol, keytype='SYMBOL', column=c('ENTREZID', 'GENENAME'))
row.names(tmp) = tmp$SYMBOL
tmp
immuno_target$entrez = tmp[immuno_target$symbol, 'ENTREZID']
immuno_target

load("rnaseq/SCC.Sample.2Batch.hg19KnownGene.RData")
colnames(fpkm_cnt_mtx) = gsub("-", "_", colnames(fpkm_cnt_mtx))
ov = intersect(sel, colnames(fpkm_cnt_mtx))
grant.fpkm = fpkm_cnt_mtx[, 1:21]
grant.fpkm = fpkm_cnt_mtx[, 22:37]

for (i in 1:nrow(immuno_target)){
	gene = immuno_target$symbol[i]
	geneid = immuno_target$entrez[i]
	fdr = round(res5$padj[res5$symbol == gene],3)
	log2fc = round(res5$log2FoldChange[res5$symbol == gene], 2)
	if(is.na(gene)){next}
	tmp = grant.fpkm[geneid, ]
	#tmp[, log_reads := log2(value)]
	tmp = data.table(sampleName = names(tmp), fpkm=tmp)
	tmp[, T12 := substr(sampleName, 12, 13)]
	tmp[, patient := substr(sampleName, 1, 10)]
	tmp
	tmp.w = dcast(tmp, patient ~ T12, value.var='fpkm'); tmp.w
	mbar = apply(tmp.w[,2:3], 2, median, na.rm=T) 
	mbar = as.data.table(mbar)
	mbar
	g = ggplot(tmp, aes(x = T12, y = fpkm))
	g = g + geom_jitter(width=0.2)
	g = g + geom_segment(data=mbar, aes(x = c(.7, 1.7), y = mbar, xend = c(1.3, 2.3), yend = mbar), alpha=.5, lwd=1.5, colour='blue')
	g = g + geom_segment(data=tmp.w, aes(x = 1.2, y = T1, xend = 1.8, yend = T2), alpha=0.2, colour='red')
	g = g + ylab('FPKM') + ggtitle(gene) + xlab('')
	g = g + theme(plot.title = element_text(face="bold", size=10, hjust=.5))
	ggsave(g, file=paste0('rnaseq/res5/dotplot/dots_fpkm_batch2_', gene, '_', log2fc, '_', fdr, '.pdf'), width=2.5, height=3)
}

cmo = fread("data_RNA_Seq_v2_expression_median.txt")
cmo[, 2:ncol(cmo)] -> tmp
head(tmp)
tmp[tmp == 0] = 1
tmp = log2(tmp)
cmo[, 2:ncol(cmo)] <- tmp
head(cmo)

head(gs)
gs = sub(".*_", "", row.names(res.cc.log2))
rsem$log2fpkm = log2(rsem$fpkm+1)
head(gs)
for (i in 1:nrow(immuno_target)){
	gene = immuno_target$symbol[i]
	geneid = immuno_target$entrez[i]
	if(is.na(gene)){next}
	tmp = res.cc.log2[gs == gene, ]
	#t1 = tmp[seq(from=1, to=16, by=2)]
	#t2 = tmp[seq(from=2, to=16, by=2)]
	log2fc = round(res$log2FoldChange[gs==gene],2)
	pval = round(res$pvalue[gs==gene], 2)
	#if(nrow(tmp) < 1) next
	tmp = unlist(tmp)
	tmp = data.table(sampleName = names(tmp), Freq=tmp)
	tmp
	tmp[, T12 := substr(sampleName, 12, 13)]
	tmp[, patient := substr(sampleName, 1, 10)]
	tmp
	tmp.w = dcast(tmp, patient ~ T12, value.var='Freq'); tmp.w
	mbar = apply(tmp.w[,2:3], 2, mean, na.rm=T) 
	mbar = as.data.table(mbar)
	mbar
	g = ggplot(tmp, aes(x = T12, y = Freq))
	g = g + geom_jitter(width=0.2)
	g = g + geom_segment(data=mbar, aes(x = c(.7, 1.7), y = mbar, xend = c(1.3, 2.3), yend = mbar), alpha=.5, lwd=1.5, colour='blue')
	g = g + geom_segment(data=tmp.w, aes(x = 1.2, y = T1, xend = 1.8, yend = T2), alpha=0.2, colour='red')
	g = g + ylab('log2(#reads)') + ggtitle(paste0(gene, '(log2FC=', log2fc, ' p=', pval, ')')) + xlab('')
	g = g + theme(plot.title = element_text(face="bold", size=10, hjust=.5))
	ggsave(g, file=paste0('rnaseq/res5/dotplot/dots_res_', gene, '_', log2fc, '_', pval, '.pdf'), width=2.5, height=3)
}

scc.maf[, .(HGVSp_Short, SWISSPROT, Transcript_ID, AA_MAF, Feature, HGVSp, HGVSc, HGVS_OFFSET)]
grep("HGvs", colnames(scc.maf), ignore.case=T)
colnames(scc.maf)[c(35, 36, 37, 97)]

immuno_target
grant.fpkm[968,]
fpkm_cnt_mtx[968, 1:21]
median(fpkm_cnt_mtx[968, 1:21]) # batch 1
median(fpkm_cnt_mtx[968, 22:37]) # batch2
median(unlist(cmo[Hugo_Symbol == 'CD68',2:ncol(cmo)]))
res.cc[grep("CD68", row.names(res.cc)),17:36]
median(as.numeric(res.cc[grep("CD68", row.names(res.cc)), 17:36]))
rsem.counts.fil[grep("CD68", row.names(rsem.counts.fil)),]

head(res.cc)
res5.cc.4 = res5.cc.2
res5.cc.4[,9:16] = log2(res5.cc.2[,9:16]+1)
res5.cc.4[,1:8] = log2(res5.cc.2[,1:8]+1)
tmp = res5.cc.4[,9:16]-res5.cc.4[,1:8]
res5.cc.4$log2FC = apply(tmp, 1, mean)
res5.cc.4$pvalue = apply(tmp, 1, function(x){unlist(t.test(x[1:8])[3])})
res5.cc.4$padj = p.adjust(res5.cc.4$pvalue)
res5.cc.4 = res5.cc.4[order(res5.cc.4$pvalue), ]
head(res5.cc.4)


pvalue = apply(res5.cc.4[,1:16], 1, function(x){unlist(wilcox.test(x[1:8], x[9:16])[3])})
res5.cc.4$pvalue = pvalue
res5.cc.4 = res5.cc.4[order(res5.cc.4$pvalue), ]
head(res5.cc.4)
t.test(tmp[1,])

head(rsem.counts)
rsem.length = rsem$length
apply(rsem.counts, 2, sum)

## netMHCpan 
# export the protein sequence

coverage = 32931921
tmp = scc.maf[, .(muts = .N), by = Tumor_Sample_Barcode]
tmp[, t12 := substr(Tumor_Sample_Barcode,14,15)]
tmp[, pat := substr(Tumor_Sample_Barcode,1,13)]
tmp = dcast(tmp, pat ~ t12, value.var='muts')
tmp[, t1.perMb := 1000000*T1/coverage]
tmp[, t2.perMb := 1000000*T2/coverage]
mean(tmp$t1.perMb)
mean(tmp$t2.perMb)

