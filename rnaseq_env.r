options(width=222)
options(warnings = -1)

library(log4r)

cwd = '/ifs/work/solitlab/huw/solit/study/hiseq/Project_0706'
if(!dir.exists(cwd)) {dir.create(cwd)}
resDir = cwd

logfile = paste0(cwd, '_', resDir, '_log')
create.logger() -> logger
logfile(logger) = logfile
level(logger) = 'INFO' 

library(data.table)
library(DESeq2)
library(gplots);
library(pheatmap);
library(ggplot2)
library(tximport)

library(readxl)
library(survival)
library(RColorBrewer)
library(Cairo)

## target file
## base: sample ID
## input: folder that contain fastq files
## fastq1.file full path 
## fastq2.file full path 
## group

cfg = new.env()

cfg$resDir = 'r_fang'
cfg$inPre = 'AKT_E17K'
cfg$pre = 'AKT'
cfg$species = 'b37'
cfg$assay = 'rnaseq'

setwd(cwd)

source('~/pipeline/configure_env.R')

fastq.dir = paste0('/ifs/archive/BIC/share/solitd/', c("MOMO_0244_BCBL3EANXX", "MOMO_0247_BCBGPCANXX", "MOMO_0250_BCBN5PANXX"))
fastq.dir

system('echo "" > target')
for(dd in fastq.dir){
	system(paste0('find ', dd, ' -name *fastq.gz >> target'))
}

fastq.files = scan('target', character())
target.ori = data.table(fastq = fastq.files, 
		    samplename = sub("(DS-bla-.*)_IGO.*", "\\1", basename(fastq.files)))
target.ori[, samplename := gsub("-", "_", samplename)]
target.ori[, r12 := sub(".*_(R.)_.*", "\\1", basename(fastq.files))]
target.ori[, base := basename(fastq.files)]
fwrite(target.ori, file='target.xls', sep="\t", quote=F)
source('~/program/fun/sync.r')
sync('.')
target.ori
target = dcast(target.ori, formula = samplename ~ r12, fun.aggregate =function(x){return(paste(x, collapse=","))}, value.var = 'fastq')
target

fread('target', header=F) -> target
setnames(target, 1, 'fastq')
target[, R12 := 'R1']
target[grep("R2_00", fastq), R12 := 'R2']            
target[, samplename := sub(".*/(.*)_IGO.*", "\\1", fastq)]            
target[, group := sub("(.*)_.", "\\1", samplename)]            
target = dcast(target, samplename + group ~ R12, value.var = 'fastq')
target[, base := samplename]
target

target[, rsem.jobname := paste0(base, ".rsem")]
target[, rsem.cmd := bsub.head(rsem.jobname, mem=55, cpu = 9)]
target[, rsem.cmd := paste0(rsem.cmd,  " \"", cfg$RSEM, "/rsem-calculate-expression  --paired-end -p 8 --star --star-path ", cfg$STAR_DIR,  " --gzipped-read-file --append-names --estimate-rspd --output-genome-bam ")]
target[, rsem.cmd := paste0(rsem.cmd, fastq1.file, " ", fastq2.file, " ",  cfg$RSEM_REF, " ", input, "/", base, " \"")]
target[1, rsem.cmd]

exe.jobs(target$rsem.cmd, logger)

## bam files
target[, bamfile := paste0(input, '/', base, '.bam', )]

wait4jobs(target$rsem.jobname, logger)

## RSeQC geneBody_coverage.py
bamfiles = paste(target$bamfile, collapse=',')
rseqc.jobname = paste0(cfg$pre, ".rseqc")
rseqc.cmd = bsub.head(rseqc.jobname, mem=55, cpu = 9, postdone=paste(target$rsem.jobname, collapse=' '))
rseqc.cmd = paste0(rseqc.cmd,  " \"", cfg$anaconda2, '/bin/geneBody_coverage.py ') 
rseqc.cmd = paste0(rseqc.cmd,  " -i ", bamfiles, ' -r ', cfg$hsa.housekeeping.bed, ' -o ', cfg$rseqcDir)
rseqc.cmd = paste0(rseqc.cmd,  " \"")
rseqc.cmd

exe.jobs(target$rseqc.cmd, logger)

## fastqc
fastqc.jobname = paste0("fastqc")
fastq = paste(c(target$R1, target$R2), collapse=" ")
fastqc.cmd = bsub.head(fastqc.jobname, mem=10, cpu=20, postdone=paste(target$rsem.jobname, collapse=' '))
fastqc.cmd = paste0(fastqc.cmd, " \"", cfg$fastqc, " -o ", cfg$fastqcDir, " -j ", cfg$java, " -f fastq -t 20 ", fastq, " \"")
fastqc.cmd
         
exe.jobs(fastqc.cmd, logger)

##
target[, countfile := paste0(sampledir, '/', samplename, '.genes.results')]
file.exists(target$countfile)

## import reads 
rsem = tximport(target$countfile, type = "rsem")
rsem$counts -> rsem
mode(rsem) = 'integer'
colnames(rsem) = target$base
head(rsem)

save(rsem, file='rsem.RData')

design = data.frame(
        row.names       = colnames(target$ID),
        condition       = target$condition,
        libType         = rep("PE", nrow(target)));

ddsmat = DESeqDataSetFromMatrix(countData = rsem,
        colData = design,
        design = ~ condition);

dds.ds <- estimateSizeFactors(ddsmat);
dds <- DESeq(dds.ds, parallel=T);

resultsNames(dds)

res = results(dds, contrast = c('condition', 'cond1', 'cond2'))
res = res[order(res$pvalue),]
save(res, file='res.Rdata')

## export data
source("~/program/fun/write_rnk.r")
rnkfile = ''
write_rnk(data = x, filename=rnkfile, pvalue.cutoff=1)
source("~/program/fun/run_gsea.R")
run_gsea(rnkfile )

## gsea summary
system('find gsea -name "*gsea_report_for*.xls" > gsea.list')
scan('gsea.list', character()) -> fl
fl[grep('h.all', fl)]
fl
lapply(fl, fread) -> fl.list
names(fl.list) = fl
as.data.table(rbindlist(fl.list, idcol = TRUE)) -> fl.dt
fl.dt[grep('rt4', .id), cell := 'RT4']
fl.dt[grep('mghu3', .id), cell := 'MGHU3']
fl.dt[grep('gata3', .id), shrna := 'GATA3']
fl.dt[grep('foxa1', .id), shrna := 'FOXA1']
fl.dt[grep('_neg_', .id), updn := 'DN']
fl.dt[grep('_pos_', .id), updn := 'UP']
fl.dt[, geneset := sub(".*rnk.(.*).v5.0.symbols.gmt.G.*", "\\1", .id)]
fl.dt[grep('h.all', .id), geneset := 'Hallmark']
fl.dt[, qval := `FDR q-val`]
fl.dt[, cell := factor(cell, levels=c('RT4', 'MGHU3'))]
fl.dt[, name := sub("HALLMARK_", "", NAME)]
fl.dt[, name := gsub("_", " ", name)]
fl.dt[, cellshrna := paste0(cell, ', shRNA ', shrna)]
fl.dt

for(gs.name in unique(fl.dt$geneset)){
	tmp = fl.dt[grep(gs.name, geneset),][qval < 0.001 & abs(NES) > 2,]
	gg = ggplot(tmp, aes(x=reorder(name, NES), y=NES, group = name)) + geom <- jitter(width=0.2, height=.05) +
		theme(axis.text.x = element <- text(angle=-90, vjust = 0.5, hjust=0)) +
		ylab("NES") + xlab("Gene set names") + coord_flip()
	ggsave(gg, file=paste0('r_fang/res/gsea_', gs.name, '.pdf'), width=7, height=8)
}
sync('res')

source('~/program/fun/sync.r')
sync('gsea')

res.cc = counts(dds, normalized=T);
save(res.cc, file='counts_normalized.RData')


rld = rlog(dds, blind = F)
vsd = vst(dds, blind=F)

## sample dist
sample.dists = dist(t(assay(vsd)))
sample.dist.mtx = as.matrix(smaple.dists)
rownames(sample.dist.mtx) = paste(vsd$dex, vsd$cell, sep='-')
colnames(sample.dist.mtx) = NULL
clrs = colorRampPalette(rev(brewer.pal(9, 'Blues')))(255)
pdf(file='sample_dist.pdf')
pheatmap(sample.dist.mtx, 
	 clustering_distance_rows = sample.dists, 
	 clustering_distance_cols = sampleDists,
	 col = clrs)
dev.off()
## end sample dist

## pca plot
pdf(file='pca.pdf')
plotPCA(vsd, intgroup='condition')
dev.off()
## end pca plot

## heatmap of sig genes
res.sig = res[!is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 1 & !is.na(res$padj) & res$padj < 0.05,]
res.sig.cc = res.cc[row.names(res.cc) %in% row.names(res.sig), ]
head(res.sig.cc)

anno.col = data.frame(row.names =target$samplename,
		      sampleClass = target$condition)
anno.col
#anno.row = data.frame()

row.sel = anno.row[]
col.sel = anno.col[]

breaks = c( seq(from = -2, to = 2, length.out = 30))
breaks = unique(breaks)
clrs = greenred(length(breaks + 1))

res.sig.cc -> mtx
cn = colnames(mtx)
mtx = t(apply(mtx, 1, scale))
colnames(mtx) = cn
mtx[mtx > 2] = 2
mtx[mtx < -2] = -2
heatmapfile = paste0()
pdf(heatmapfile, width=10, height=6)
hp = pheatmap(mtx, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks, border_color = NA,
        main='', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno.row, annotation_row = anno.col,
dev.off()

## venn diagram
source("~/program/vennDia.R")
venfile = paste0()
pdf(venfile, width=6, height=6)
venndiagram(x = xx, y = yy, type= '2', labels = c("", ""))
dev.off()

## tcga RNA-seq
xx = load("../bcg/tcga_blca.RData");
xx
blca.cc.log2 = log2(blca.cc + 1)
head(blca.cc.log2)
blca.assay
blca.rows$id = paste(blca.rows$ensembl_gene_id, blca.rows$external_gene_name, sep = "_")
blca_count_symbol = blca.assay
row.names(blca_count_symbol) = blca.rows[row.names(blca_count_symbol), "id"]
blca_count_symbol[1:10, 1:2]

## pca analysis
pcafile = paste0()
pdf(pcafile)
rd = plotPCA(rld, intgroup = 'condition', returnData = T)
dev.off()

pcafile = paste0()
pdf(pcafile, width=8, height=6)
ggplot(rd, aes(PC1, PC2, label = samp_class,  color=samp_class)) + geom_point() + theme(legend.position="right") 
dev.off()

pcafile = paste0()
percentvar = round(100 * attr(rd, 'percentVar'))
pdf('res/pca_2.pdf')
ggplot(rd, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ", percentvar[1], "% variance")) +
        ylab(paste0("PC2: ", percentvar[2], "% variance")) +
        coord_fixed()
dev.off()


# PCA by using top 100 genes in bladder from https://www.gtexportal.org/home/
tmp = prcomp(t(res.cc))
tmp = tmp$x
tmp = as.data.frame(tmp)

pcafile = paste0()
pdf(pcafile, width=8, height=6)
p = ggplot(tmp, aes(PC1, PC2, color=samp_class)) + geom_point() + theme(legend.position="right"); print(p)
dev.off()

## 
load("../bcg/uc_basal_luminal_marks.RData")

setdiff(base47, c(base47.luminal, base47.basal))

colnames(tcga_shen_cc_log)
anno_col = data.frame(row.names = colnames(tcga_shen_cc_log),
        Tumor_Stage = c(blca.data$tumor_stage, rep("organoids", ncol(rsem_shen$counts))), 
        Definition = c(blca.data$definition, rep("organoids", ncol(rsem_shen$counts))),
        stringsAsFactors = F)
dim(anno_col)

anno_col$samp_class= "TCGA BLCA"
anno_col$samp_class[grep("Org", row.names(anno_col))] = "organoid"
anno_col$samp_class[grep("tumor", row.names(anno_col))] = "tumor"
tail(anno_col)

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

anno_col$SC = "NOT"
anno_col$SC[434:474] = "organoids"
anno_col$SC[substr(row.names(anno_col), 1, 12) %in% tcga_sc_id] = "SC"
anno_col$SC[substr(row.names(anno_col), 1, 12) %in% tcga_basal_id] = "Basal"
anno_col$SC[substr(row.names(anno_col), 1, 12) %in% intersect(tcga_sc_id, tcga_basal_id)] = "Basal & SC"
table(anno_col$SC)

unique(anno_col$Definition)

load("../../genelist/early117.RData")
read.table("../../genelist/early_117_subset.txt", head=F) -> early.subset
colnames(early.subset) = c('class', 'genesymbol')
early.subset

gs = row.names(tcga_shen_cc_log)
gs = sub(".*_", "", gs)
anno_row = data.frame(row.names=row.names(tcga_shen_cc_log), 
                      BASE47 = rep("NOT", nrow(tcga_shen_cc_log)), stringsAsFactors = F);
anno_row$BASE47[gs %in% base47.luminal] = "Luminal";
anno_row$BASE47[gs %in% base47.basal] = "Basal";
anno_row$McConkey = "NOT"
anno_row$McConkey[gs %in% mcconkey.luminal] = "Luminal"
anno_row$McConkey[gs %in% mcconkey.basal] = "Basal"
anno_row$Early117 = "NOT"
anno_row$Early117[gs %in% early.subset$genesymbol[early.subset$class == 'Class12']] = "Class12"
anno_row$Early117[gs %in% early.subset$genesymbol[early.subset$class == 'Class13']] = "Class13"
anno_row$Early117[gs %in% early.subset$genesymbol[early.subset$class == 'Class2']] = "Class2"
anno_row$sd = apply(tcga_shen_cc_log, MARGIN = 1, FUN = function(x) sd(x) != 0)
head(anno_row$sd)

breaks = c( seq(from = -2, to = 2, length.out = 30))
breaks = unique(breaks)
col = greenred(length(breaks + 1))


## immune cell infiltration
## ssGSEA 
immuno_sig =fread("~/program/immuno_sig.txt")
gg = unique(immuno_sig$CellType)
geneset = list(immuno_sig[CellType == gg[1], Symbol])
for(i in 2:length(gg)){
	geneset = c(geneset, list(immuno_sig[CellType == gg[i], Symbol]))
}
names(geneset) = gg
length(geneset)
names(geneset)

adaptive = c( "Bcells", "Tcells", "Thelpercells", "Tcm", "Th1cells", "Th2cells", "TFH", "Th17cells", "TReg", "CD8Tcells", "Tgd", "Cytotoxiccells")
innate = c( "NKcells", "NKCD45dimencells", "NKCD56brightcells", "DC", "iDC", "aDC", "pDC", "Eosinophils", "Macrophages", "Mastcells", "Neutrophils") 
anno_row = data.frame(row.names = names(geneset), ImmuneType = rep("Adaptive", 28), stringsAsFactors = F)
anno_row$ImmuneType[row.names(anno_row) %in% innate] = "Innate"
names(geneset)
head(geneset)

gsva_score  = gsva(res.cc, geneset, method="ssgsea")

library(limma)
design <- model.matrix(~ substr(colnames(gsva_score), 12, 13))
colnames(design) <- c("ALL", "T2T1")
head(design)
fit <- lmFit(gsva_score, design)
fit <- eBayes(fit)
gsva.res <- topTable(fit, coef="T2T1", number=Inf)
gsva.res = gsva.res[order(gsva.res$P.Value),] 
gsva.res

gsva.res = as.data.table(gsva.res, keep.rownames=T)
gsva.res[, Significance := 'Not']
gsva.res[P.Value < 0.05, Significance := 'Sig']
gsva.res

gg = ggplot(gsva.res, aes(x=reorder(rn, P.Value), y=logFC, color = Significance, size=-log(P.Value))) + geom_point() +
	theme(axis.text.x = element_text(angle=-90, vjust = 0.5, hjust=0)) + ggtitle('Immune cell infilteration \nby ssGSEA (T2 vs T1) ') +
	ylab("logFC") + xlab("Immune cells") + coord_flip() +
	geom_hline(yintercept=0, color=adjustcolor(2, .3))
ggsave(gg, file='r_fang/res/immune_cell_deconv_ssGSEA_T2_T1.pdf', width=5, height=5)
sync('r_fang/res')

## need revise
tmp = as.data.table(melt(gsva_score))
= merge(gsva.score.sel.melt, gsva.anno.sel.dt, by.x = 'Var2', by.y = 'rn')
gsva.score.sel.melt
gg = ggplot(gsva.score.sel.melt, aes(x=type, y=value)) + geom_jitter(width=.2) +
	ylab('Score in ssGSEA analyiss') + facet_wrap(~Var1, scales = 'free')  +
	stat_summary(fun.data = "mean_cl_boot", colour = "red", size = .2) +
	scale_y_continuous(trans='log2')
ggsave(gg, file='r_fang/res/dotplots_gsva_score_primary_secondary.pdf', width=12, height=12)
sync('r_fang/res')

# need revise
gs = sub(".*_", "", row.names(res.cc.log))
res.cc.log[gs %in% immuno_sig$Symbol, ] -> tmp
#my_scale(tmp) -> tmp
cn = colnames(tmp)
matrix_tmp = t(apply(tmp, 1, scale))
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
colnames(matrix_tmp) = cn
breaks = c( seq(from = -2, to = 2, length.out = 30))
remove_null_row(matrix_tmp) -> matrix_tmp
pdf("res/heatmap_immuno_sig.pdf")
hp = pheatmap(tmp, color=greenred(length(breaks)+1), scale='none', clustering_method = 'ward.D2',
	              cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 9, show_rownames = F)
dev.off()

## immune cell infiltration
## cybersort 
tmp = as.data.frame(res.cc)
gs = sub(".*_", "", row.names(tmp))
gs1 = gs[!duplicated(gs)]
tmp = tmp[!duplicated(gs),]
row.names(tmp) = gs1
rm(gs1, gs)
write.table(tmp, file='tmp.tsv', sep="\t", quote=F) ## add a tab in the head of this file!!
source('~/program/cibersort/CIBERSORT.R')
ciber = CIBERSORT(paste0(cfg$cibersort, '/LM22.txt'), 'tmp.tsv', perm=1000, QN=TRUE, absolute=T)
ciber.t = t(ciber)
ciber.t = ciber.t[1:22, ]

design <- model.matrix(~ condition)
colnames(design) <- c("ALL", "PS")
fit <- lmFit(ciber.t.ps, design)
fit <- eBayes(fit)
summary(decideTests(fit, p.value=0.05))
ciber.res <- topTable(fit, coef="PS", number=Inf)
ciber.res = as.data.table(ciber.res, keep.rownames=T)
ciber.res[, Significance := 'Not']
ciber.res[P.Value < 0.05, Significance := 'Sig']
ciber.res

save(ciber.res, ciber.t.ps, ciber.t, file='cibersort_score_primary_vs_secondary.RData')

gg = ggplot(ciber.res, aes(x=reorder(rn, P.Value), y=logFC, color = Significance, size=-log(adj.P.Val))) + geom <- point() +
	theme(axis.text.x = element <- text(angle=-90, vjust = 0.5, hjust=0)) + ggtitle('Immune cell infilteration \nby cibersort (primary vs secondary) ') +
	ylab("logFC") + xlab("Immune cells") + coord <- flip() +
	geom <- hline(yintercept=0, color=adjustcolor(2, .3))
ggsave(gg, file='res/immune_cell_deconv_primary_vs_secondary_cibersort.pdf', width=6, height=5)
sync('res')

ciber.t.ps.melt = as.data.table(melt(ciber.t.ps))
ciber.t.ps.melt = merge(ciber.t.ps.melt, gsva.anno.sel.dt, by.x = 'Var2', by.y = 'rn')
ciber.t.ps.melt
ciber.res
gg = ggplot(ciber.t.ps.melt, aes(x=type, y=value)) + geom <- jitter(width=.2) +
	ylab('Score in CIBERSORT analyiss') + facet <- wrap(~Var1, scales = 'free')  +
	stat <- summary(fun.data = "mean_cl_boot", colour = "red", size = .2) +
	scale <- y <- continuous(trans='log2')
ggsave(gg, file='res/dotplots_cibersort_score_primary_secondary.pdf', width=12, height=12)
sync('res')

cn = colnames(ciber.t)
matrix <- tmp = t(apply(ciber.t, 1, scale))
matrix <- tmp[matrix <- tmp > 2] = 2
matrix <- tmp[matrix <- tmp < -2] = -2
colnames(matrix <- tmp) = cn
breaks = c( seq(from = -2, to = 2, length.out = 30))
pdf("res/heatmap_cibersort_score.pdf", width=22, height=5)
hp = pheatmap(matrix <- tmp, color=greenred(length(breaks)+1), scale='none', clustering <- method = 'complete',
	      cluster <- rows = TRUE, cluster <- cols = TRUE, show <- colnames = F, fontsize <- row = 9, show <- rownames = T)
dev.off()
sync('res')

## heatmap of the immuno sig genes from cibersort
ciber.gs = fread(paste0(cfg$cibersort, '/LM22.txt'))
ciber.gs = ciber.gs$`Gene symbol`
ciber.gs
gs = sub(".*_", "", row.names(res.cc.log))
res.cc.log[gs %in% ciber.gs, ] -> tmp
#my <- scale(tmp) -> tmp
cn = colnames(tmp)
matrix <- tmp = t(apply(tmp, 1, scale))
matrix <- tmp[matrix <- tmp > 2] = 2
matrix <- tmp[matrix <- tmp < -2] = -2
colnames(matrix <- tmp) = cn
breaks = c( seq(from = -2, to = 2, length.out = 30))
remove <- null <- row(matrix <- tmp) -> matrix <- tmp
pdf("res/heatmap_immuno_sig_cibersort.pdf")
hp = pheatmap(tmp, color=greenred(length(breaks)+1), scale='none', clustering <- method = 'ward.D2',
	      cluster <- rows = TRUE, cluster <- cols = TRUE, show <- colnames = T, fontsize <- row = 9, show <- rownames = F)
dev.off()


#blca.cc[c(base47.basal.ensembl, base47.luminal.ensembl),] -> blca_base47
#dim(blca_base47)

#row.names(blca_base47) = blca.rows[row.names(blca_base47), "external_gene_name"]

#anno_row = data.frame(row.names = row.names(blca_base47), BASE47 = rep("luminal", nrow(blca_base47)), stringsAsFactors = F)
#anno_row$BASE47[row.names(anno_row) %in% base47.basal] = "basal"

# add annotation of bladder cluster of 129 samples
anno_col$nat_cl = "NOT"
anno_col$nat_cl[grep("MS0", row.names(anno_col))] = "organoids"
anno_col[row.names(nat2014_cluster), "nat_cl"] = nat2014_cluster[, "cluster"]
table(anno_col$nat_cl)
dim(anno_col)
head(nat2014_cluster)


colnames(anno_col)
unique(anno_col$Tumor_Stage)
unique(anno_col$Definition)
unique(anno_col$SC)
unique(anno_col$nat_cl)
colnames(anno_col)
unique(anno_col$Tumor_Stage)
unique(anno_col$Definition)
unique(anno_col$SC)
unique(anno_col$nat_cl)
colnames(anno_row)
unique(anno_row$BASE47)
unique(anno_row$McConkey)
unique(anno_row$Early117)
unique(anno_col$passage)

del = c("MS020_SuB9_9", "MS029_SuB38_2", 'xx', 'yy')
anno_col$del = 0
anno_col$del[row.names(anno_col) %in% del] = 1

#anno_color
ann_colors = list(
mark = c(basal = "#1B9E77", luminal = "#D95F02"),
#type = c(primary = "#7570B3", secondary = "#E7298A", Shen = "#02427D",  BLCA= "#66A61E"),
type = c(Shen = "#02427D",  BLCA= "#66A61E"),
BASE47 = c(NOT="#cde6ec", Luminal="#E7298A", Basal="#02427D"),
McConkey = c(NOT="#cde6ec", Luminal="#E7298A", Basal="#02427D"),
Early117 = c(NOT="#cde6ec", Class2="#E7298A", Class12="#02427D", Class13="#299d02"),
#Early117 = c(NOT="#cde6ec", Uromol="#E7298A"),
SC = c(Basal="#5351c7", NOT="#cde6ec", SC="#e81b1b", "Basal & SC" ="#f20bda", organoids="#9d7802"),
nat_cl = c(c1="#eec8c2", c2="#9190d7", c3="#f05840", c4 ="#0a07ab", NOT="#cde6ec", organoids="#9d7802"),
Tumor_Stage = c("not reported"="#cde6ec", "stage i" = "#adb6ec", "stage ii" = "#9190d7", "stage iii" = "#5c72eb", "stage iv" = "#0a07ab", organoids="#9d7802" ),
BASE47_129 = c("NOT"="#cde6ec", "C1" = "#eec8c2", "C2" = "#9190d7", "C3" = "#f05840", "C4" = "#0a07ab" ),
Definition= c('Primary solid Tumor' ="#dedae8",  'Solid Tissue Normal' = "#fdd4c1", organoids="#9d7802"),
samp_class = c("TCGA BLCA"='#cde6ec', 'organoid'=adjustcolor('#9d7802', alpha.f=.3), 'tumor'='#9d7802'),
passage = c(TCGA="#ffffff", "P14"=adjustcolor("#2a027d", alpha.f = .06*2),"P13"=adjustcolor("#2a027d", alpha.f = .06*3),
           "P12"=adjustcolor("#2a027d", alpha.f = .06*3.5),   "P11"=adjustcolor("#2a027d", alpha.f = .06*4), 
           "P10"=adjustcolor("#2a027d", alpha.f = .06*4.5),   "P9"=adjustcolor("#2a027d", alpha.f = .06*5), 
           "P8"=adjustcolor("#2a027d", alpha.f = .06*5.5),    "P7"=adjustcolor("#2a027d", alpha.f = .06*6), 
           "P6"=adjustcolor("#2a027d", alpha.f = .06*7),      "P5"=adjustcolor("#2a027d", alpha.f = .06*8), 
           "P4"=adjustcolor("#2a027d", alpha.f = .06*9),      "P3"=adjustcolor("#2a027d", alpha.f = .06*10), 
           "P2"=adjustcolor("#2a027d", alpha.f = .06*11),     "P1"=adjustcolor("#2a027d", alpha.f = .06*12), 
           "P0"=adjustcolor("#2a027d", alpha.f = .06*13),     "tumor"=adjustcolor("#2a027d", alpha.f = .06*14)
	   )
)
#BASE47_129 four clusters by BASE47 genes on 129 samples

## BASE47 classifiers 129 samples
col_sel = anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT" & anno_col$del == 0 
row_sel = anno_row$BASE47 != "NOT" & anno_row$sd
anno_col_tmp = anno_col[col_sel, c("nat_cl", "Tumor_Stage", "SC")]
anno_row_tmp = anno_row[row_sel, "BASE47", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_base47v2.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks, border_color = NA,
        main='TCGA BLCA RNA-Seq 129 samples, BASE47 genes', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()

tail(anno_col)
table(anno_col_tmp$nat_cl, anno_col_tmp$SC)

# assign the clusters of the 129 patient samples based on the BASE47 genes
anno_col$BASE47_129 = rep("NOT", nrow(anno_col))
row.names(anno_col_tmp) -> tmp
tmp = tmp[hp$tree_col$order]
anno_col$BASE47_129[row.names(anno_col) %in% tmp[1:51]] = "C1"
anno_col$BASE47_129[row.names(anno_col) %in% tmp[52:68]] = "C2"
anno_col$BASE47_129[row.names(anno_col) %in% tmp[69:94]] = "C3"
anno_col$BASE47_129[row.names(anno_col) %in% tmp[95:129]] = "C4"

## 51, 35, 25, 18
col_sel = anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT" & anno_col$del == 0
row_sel = anno_row$BASE47 != "NOT" & anno_row$sd
anno_col_tmp = anno_col[col_sel, c("nat_cl", "Tumor_Stage", "SC", "BASE47_129")]
anno_row_tmp = anno_row[row_sel, "BASE47", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_base47_cluster_v2.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks, border_color = NA,
        main='TCGA BLCA RNA-Seq 129 samples, BASE47 genes', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()


# survive curve for these nat_cl clusters
library(survival)
death = blca.data$days_to_death
death = death[anno_col$nat_cl != "NOT"]
followup = blca.data$days_to_last_follow_up
followup = followup[anno_col$nat_cl != "NOT"]
event = rep(0, length(death))
event[death > 0] = 1
time = death
time[is.na(time)] = followup[is.na(time)]
tmp = anno_col[anno_col$nat_cl != "NOT", "nat_cl"]
sf = survfit(Surv(time, event) ~ tmp)
pdf("results/surv_tcga_4cluster.pdf")
plot(sf, col=2:5, mark.time=T, cex=.5, xpd=T, xlab="days", ylab="% survived")
legend(4000, 1, col=2:5, legend=c("Basal", "Basal & SC", "NOT", "SC"), lty=1, cex=.5, bty="n", y.intersp=3)
text(4000, .5, "P value 0.187", cex=.7)
dev.off()
#data.frame(death = death, followup = followup, time=time)

sf$n
survdiff(Surv(time, event) ~ tmp)

# survive curve for these SC clusters
death = blca.data$days_to_death
death = death[anno_col$nat_cl != "NOT"]
followup = blca.data$days_to_last_follow_up
followup = followup[anno_col$nat_cl != "NOT"]
event = rep(0, length(death))
event[death > 0] = 1
time = death
time[is.na(time)] = followup[is.na(time)]
tmp = anno_col[anno_col$nat_cl != "NOT", "SC"]
tmp
tmp2 = tmp
tmp2[tmp == "C1" | tmp == "C2"] = "C1/2"
tmp2[tmp == "C3" | tmp == "C4"] = "C3/4"
sf = survfit(Surv(time, event) ~ tmp)
pdf("results/surv_tcga_pathology.pdf")
plot(sf, col=2:5, mark.time=T, cex=.5, xpd=T, xlab="days", ylab="% survived")
legend(4000, 1, col=2:5, legend=c("Basal", "Basal & SC", "NOT", "SC"), lty=1, cex=.5, bty="n", y.intersp=3)
text(4000, .5, "P value 0.0173", cex=.7)
dev.off()
sf$n
survdiff(Surv(time, event) ~ tmp)

# survive curve for these BASE47_129 clusters
death = blca.data$days_to_death
death = death[anno_col$nat_cl != "NOT"]
followup = blca.data$days_to_last_follow_up
followup = followup[anno_col$nat_cl != "NOT"]
event = rep(0, length(death))
event[death > 0] = 1
time = death
time[is.na(time)] = followup[is.na(time)]
tmp = anno_col[anno_col$nat_cl != "NOT", "BASE47_129"]
tmp2 = tmp
tmp2[tmp == "C1" | tmp == "C3"] = "C1/3"
tmp2[tmp == "C2" | tmp == "C4"] = "C2/4"
sf = survfit(Surv(time, event) ~ tmp2)
pdf("results/surv_tcga_base47_base47cluster.pdf")
plot(sf, col=2:3, mark.time=T, cex=.5, xpd=T, xlab="days", ylab="% survived")
legend(4000, 1, col=2:3, legend=c("C1/C3 n=75", "C2/C4 n=53"), lty=1, cex=.5, bty="n", y.intersp=3)
text(4000, .4, "P value 0.02", cex=.7)
dev.off()
sf$n
survdiff(Surv(time, event) ~ tmp2)

#genes from nature 2014
tmp = c("FGFR3", "KRT5", "KRT6A", "KRT14", "EGFR", "GATA3", "FOXA1", "UPK3A", "ERBB2", "ESR2", "CDH1")
col_sel = anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT"
row_sel = gs %in% tmp & anno_row$sd
anno_col_tmp = anno_col[col_sel, c("nat_cl", "Tumor_Stage", "SC")]
anno_row_tmp = anno_row[row_sel, "BASE47", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_natureSig.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=col, scale='none', show_rownames = T, breaks = breaks, border_color = NA, main='TCGA BLCA RNA-Seq 129 samples, Nature 2014 gene list', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, #annotation_row = anno_row_tmp, 
        annotation_colors = ann_colors)
dev.off()


#BASE47 genes, all samples 
col_sel = anno_col$Definition == "Primary solid Tumor" 
row_sel = anno_row$BASE47 != "NOT" & anno_row$sd 
anno_col_tmp = anno_col[col_sel, c("nat_cl", "Tumor_Stage", "SC", "BASE47_129")]
anno_row_tmp = anno_row[row_sel, "BASE47", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = 2
pdf("results/heatmap_tcga414_base47.pdf", width=10, height=7)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks, border_color = NA,
        main='TCGA BLCA RNA-Seq 414 samples, BASE47 genes', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 8,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()

#BASE47 genes, blca plus SHEN samples 
tcga_shen_cc_log[1:10, 404:414]
unique(anno_col$samp_class)
head(anno_col)
col_sel = anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT"| anno_col$Definition == "organoids"  & anno_col$del == 0
row_sel = anno_row$BASE47 != "NOT" & anno_row$sd 
anno_col_tmp = anno_col[col_sel, c("nat_cl", "Tumor_Stage", "SC", "samp_class")]
anno_row_tmp = anno_row[row_sel, "BASE47", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
ll = 2
breaks = c( seq(from = -ll, to = ll, length.out = 30));breaks = unique(breaks)
matrix_tmp[matrix_tmp > ll] = ll
matrix_tmp[matrix_tmp < -ll] = -ll
pdf("results2/heatmap_tcga414_organoids_base47.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks, border_color = NA,
        main='414 samples plus Organoid samples', clustering_method = 'ward.D2', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()

#BASE47 genes, SHEN samples 
col_sel = anno_col$Definition == "organoids" & anno_col$del == 0
row_sel = anno_row$BASE47 != "NOT" & anno_row$sd 
anno_col_tmp = anno_col[col_sel, c("samp_class", "passage")]
anno_row_tmp = anno_row[row_sel, "BASE47", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
dim(matrix_tmp)
row.names(matrix_tmp) = sub(".*_", "", row.names(matrix_tmp))
row.names(anno_row_tmp) = sub(".*_", "", row.names(anno_row_tmp))
pdf("results/heatmap_organoids_42samples_base47.pdf", width=7, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks, border_color = NA,
        main='Organoid/tumor samples, BASE47 genes', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 6, fontsize_col=6,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()
## cutree
## 1 basal, 2 luminal, 3 mix
cutree(hp$tree_col, 3) -> tree
tree[tree == 1] = 'basal'
tree[tree == 2] = 'luminal'
tree[tree == 3] = 'mix'
length(tree)

rsem_shen_counts = rsem_shen$counts
dim(rsem_shen_counts)
mode(rsem_shen_counts) <- "integer"
condition =  tree
design = data.frame(
        row.names       = colnames(rsem_shen_counts),
        condition       = condition,
        libType         = rep("PE", ncol(rsem_shen_counts)));
ddsmat = DESeqDataSetFromMatrix(countData = rsem_shen_counts,
        colData = design,
        design = ~ condition);
dds.ds <- estimateSizeFactors(ddsmat);
dds <- DESeq(dds.ds, parallel=T);
rsem_shen_cc = counts(dds, normalized=T)
rsem_shen_cc[rsem_shen_cc == 0] = 1
rsem_shen_cc_log = log2(rsem_shen_cc)

res_basal_luminal = results(dds, contrast = c('condition', 'basal', 'luminal'))
res_mix_luminal     = results(dds, contrast = c('condition', 'mix', 'luminal'))

#source("~/program/fun/post_deseq2_results.R")
vv = c('res_basal_luminal', 'res_mix_luminal')
vv_sig = paste0(vv, "_sig")
padj = 0.01; log2FC = 2; 
for(i in 1:length(vv)){
	v = vv[i]
	v_sig = vv_sig[i]
	x = get(v)
	x$symbol = sub(".*_", "", row.names(x))
	x[order(x$pvalue),] -> x
	x.sig = x[!is.na(x$padj) & x$padj < padj & !is.na(x$log2FoldChange) & abs(x$log2FoldChange) > log2FC,] 
	source("~/program/fun/write_rnk.r")
	write_rnk(data = x, filename=paste0(v, ".rnk"), pvalue.cutoff=1)
	write.table(x, file = paste0(v, ".csv"), sep="\t", quote=F)
	assign(v, x); 
	assign(v_sig, x.sig)
}

source("~/program/fun/run_gsea.R")
run_gsea("res_basal_luminal.rnk" )
run_gsea("res_mix_luminal.rnk")

source("~/program/vennDia.R")
pdf("results/ven_basal_mix_vs_lumina.pdf", width=6, height=6)
venndiagram(x = res_basal_luminal_sig$symbol, y = res_mix_luminal_sig$symbol, type= '2', labels = c("basal vs luminal", "mixed type vs luminal"))
dev.off()

dim(res_basal_luminal_sig)
dim(res_mix_luminal_sig)

## heatmap
id = union(row.names(res_basal_luminal_sig), row.names(res_mix_luminal_sig))
rsem_shen_cc_log[id, ] -> matrix_tmp
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = names(tree)
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
ann = data.frame(row.names = names(tree), type = tree)
pdf("results/heatmap_organoids_42sample_subtype.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = F, breaks = breaks, border_color = NA,
        main='Heatmap of differentially expressed genes in tumor/organoids subtypes', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
        annotation_col = ann)
dev.off()


#BASE47 genes, SHEN samples organoids, not primary
col_sel = anno_col$samp_class == "organoid" & anno_col$del == 0
col_sel
row_sel = anno_row$BASE47 != "NOT" & anno_row$sd 
anno_col_tmp = anno_col[col_sel, c("samp_class", "passage")]
anno_row_tmp = anno_row[row_sel, "BASE47", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
dim(matrix_tmp)
row.names(matrix_tmp) = sub(".*_", "", row.names(matrix_tmp))
row.names(anno_row_tmp) = sub(".*_", "", row.names(anno_row_tmp))
pdf("results/heatmap_organoidsonly_base47.pdf", width=7, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks, border_color = NA,
        main='Organoid/tumor samples, BASE47 genes', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 6, fontsize_col=6,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()
## cutree
## 1 basal, 2 luminal
cutree(hp$tree_col, 2) -> tree_orgonly
tree_orgonly[tree_orgonly == 1] = 'basal'
tree_orgonly[tree_orgonly == 2] = 'luminal'
length(tree_orgonly)
tree_orgonly

rsem_shen_counts = rsem_shen$counts
rsem_shen_counts = rsem_shen_counts[, grep("tumor", target$id, invert = T)]
dim(rsem_shen_counts)
mode(rsem_shen_counts) <- "integer"
condition =  tree
design = data.frame(
        row.names       = colnames(rsem_shen_counts),
        condition       = condition,
        libType         = rep("PE", ncol(rsem_shen_counts)));
ddsmat = DESeqDataSetFromMatrix(countData = rsem_shen_counts,
        colData = design,
        design = ~ condition );
dds.ds <- estimateSizeFactors(ddsmat);
dds <- DESeq(dds.ds, parallel=T);
rsem_shen_cc = counts(dds, normalized=T)
rsem_shen_cc[rsem_shen_cc == 0] = 1
rsem_shen_cc_log = log2(rsem_shen_cc)

res_basal_luminal_orgonly = results(dds, contrast = c('condition', 'basal', 'luminal'))

vv = c('res_basal_luminal_orgonly')
vv_sig = paste0(vv, "_sig")
padj = 0.01; log2FC = 2; 
for(i in 1:length(vv)){
	v = vv[i]
	v_sig = vv_sig[i]
	x = get(v)
	x$symbol = sub(".*_", "", row.names(x))
	x[order(x$pvalue),] -> x
	x.sig = x[!is.na(x$padj) & x$padj < padj & !is.na(x$log2FoldChange) & abs(x$log2FoldChange) > log2FC,] 
	source("~/program/fun/write_rnk.r")
	write_rnk(data = x, filename=paste0(v, ".rnk"), pvalue.cutoff=1)
	write.table(x, file = paste0(v, ".csv"), sep="\t", quote=F)
	assign(v, x); 
	assign(v_sig, x.sig)
}
dim(res_basal_luminal_orgonly_sig)

source("~/program/fun/run_gsea.R")
run_gsea("res_basal_luminal_orgonly.rnk" )

## heatmap
dim(res_basal_luminal_orgonly_sig)
table(sign(res_basal_luminal_orgonly_sig$log2FoldChange))
head(res_basal_luminal_orgonly_sig)
id = row.names(res_basal_luminal_orgonly_sig)
rsem_shen_cc_log[id, target$passage != 'p' ] -> matrix_tmp
dim(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = names(tree_orgonly)
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
ann = data.frame(row.names = names(tree_orgonly), type = tree_orgonly)
dim(ann)
pdf("results/heatmap_organoids_34sample_subtype.pdf", width=7, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = F, breaks = breaks, border_color = NA,
        main='Heatmap of differentially expressed genes in tumor/organoids subtypes', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 6,
        annotation_col = ann)
dev.off()


#BASE47 genes, blca plus SHEN samples , less header
col_sel = anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT"| anno_col$Definition == "organoids"   & anno_col$del == 0
row_sel = anno_row$BASE47 != "NOT" & anno_row$sd 
anno_col_tmp = anno_col[col_sel, c("samp_class"), drop = F]
anno_row_tmp = anno_row[row_sel, "BASE47", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results2/heatmap_tcga129_organoids_base47_lesshead.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks, border_color = NA,
        main='129 samples plus Organoid samples', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()

#BASE47 genes, blca plus SHEN samples, with col names
col_sel = anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT"| anno_col$Definition == "organoids" & anno_col$del == 0
row_sel = anno_row$BASE47 != "NOT" & anno_row$sd 
anno_col_tmp = anno_col[col_sel, c("nat_cl",  "SC", "samp_class")]
anno_row_tmp = anno_row[row_sel, "BASE47", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results2/heatmap_tcga129_organoids_base47_colnames.pdf", width=24, height=8)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks, border_color = NA,
        main='129 samples plus Organoid samples', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 6, fontsize_col = 7,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()

head(newid)
anno_col$patients = "TCGA"
anno_col$patients[434:475] = row.names(anno_col)[434:475]
anno_col$passage = "TCGA"
target$passage2 = paste0("P", target$passage)
target$passage2[target$passage2 == 'Pp'] = 'tumor'
anno_col$passage[434:475] = target$passage2
tail(anno_col)

#BASE47 genes, blca plus SHEN samples, colnames and passage 
col_sel = anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT"| anno_col$Definition == "organoids"  & anno_col$del == 0
row_sel = anno_row$BASE47 != "NOT" & anno_row$sd 
anno_col_tmp = anno_col[col_sel, c("nat_cl",  "SC", "samp_class", "passage")]
anno_row_tmp = anno_row[row_sel, "BASE47", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_organoids_base47_colname_passage.pdf", width=20, height=10)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks,
        main='129 samples plus Organoid samples', clustering_method = 'complete', border_color = NA,
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 8, fontsize_col = 8,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()

## patient wise
tail(matrix_tmp)
tmp = row.names(anno_col)
tmp
tmp = sub("MS..._", "", tmp)
tmp = sub("_.*", "", tmp)
tmp
anno_col$pat="TCGA BLCA"
anno_col$pat[grep("MS0", row.names(anno_col))]= tmp[grep("MS0", row.names(anno_col))]
anno_col

col_sel = anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT"| anno_col$Definition == "organoids"  & anno_col$del == 0
row_sel = anno_row$BASE47 != "NOT" & anno_row$sd 
anno_col_tmp = anno_col[col_sel, c( "SC", "samp_class", "passage", "patients")]
anno_row_tmp = anno_row[row_sel, "BASE47", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_organoids_base47_colname_passage_pat.pdf", width=20, height=10)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks,
        main='129 samples plus Organoid samples', clustering_method = 'complete', border_color = NA,
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 8, fontsize_col = 8,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()

tmp = cutree(hp$tree_col, 2)
tmp[tmp == 2] = "basal"
tmp[tmp == 1] = "luminal"
tmp = tmp[grep("MS0", names(tmp))]
base47_org_cluster = tmp
head(base47_org_cluster)
base47_org_cluster

#library(xlsx)
## see http://wlab.ethz.ch/cspa/#highlights
#hsa_cd = read.xlsx("../../genelist/Surface_protein.xlsx", 1)

dim(hsa_cd)

#all organoids samples
row.names(anno_col)
col_sel = anno_col$Definition == "organoids"  & anno_col$del == 0
tmp = apply(tcga_shen_cc_log, 1, sd)
row_sel = tmp > 2.5

colnames(anno_col)
anno_col_tmp = anno_col[col_sel, c("patients", "passage")]
anno_row_tmp = anno_row[row_sel, "BASE47", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
dim(tcga_shen_cc_log)
length(col_sel)
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
matrix_tmp[is.na(matrix_tmp)] = 0.001
dim(matrix_tmp)

pdf("results/heatmap_organoids_top1787.pdf", width=10, height=16)
hp = pheatmap(matrix_tmp, color=col, scale='none', show_rownames = F, breaks = breaks,
	      main='Organoid samples, top SD (2.5), 1787 genes ', clustering_method = 'complete', 
	      cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 12,
	      annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()

#surface protein genes, blca 129 
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT") 
row_sel = gs %in% hsa_cd$ENTREZ.gene.symbol
anno_col_tmp = anno_col[col_sel, c("SC", "Definition")]
anno_row_tmp = anno_row[row_sel, , drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
options(repr.plot.width=10, repr.plot.height=5)
pdf("results/heatmap_tcga129_CD.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = F, breaks = breaks,
        main='129 samples, cell surface protein', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_colors = ann_colors)
dev.off()

#surface protein genes, blca 129 plus SHEN samples
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT") | anno_col$Definition == "organoids" 
row_sel = gs %in% hsa_cd$ENTREZ.gene.symbol
anno_col_tmp = anno_col[col_sel, c("SC", "Definition", "passage")]
anno_row_tmp = anno_row[row_sel, , drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_organoids_CD.pdf", width=20, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = F, breaks = breaks,
        main='129 samples plus Organoid samples', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 6, fontsize_col = 6,
        annotation_col = anno_col_tmp, annotation_colors = ann_colors)
dev.off()

#surface protein genes, blca 129 plus SHEN samples, top SD
tmp = apply(tcga_shen_cc_log, 1, sd)
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT") | anno_col$Definition == "organoids" 
row_sel = gs %in% hsa_cd$ENTREZ.gene.symbol & tmp >2.5
anno_col_tmp = anno_col[col_sel, c("SC", "samp_class", "passage")]
anno_row_tmp = anno_row[row_sel, , drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_organoids_CD_topSD.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = F, breaks = breaks,
        main='129 samples plus Organoid samples', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_colors = ann_colors)
dev.off()

row.names(matrix_tmp)
tr = cutree(hp$tree_row,3)

## top SD of all genes
tmp = apply(tcga_shen_cc_log, 1, sd)
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT") | anno_col$Definition == "organoids" 
row_sel = tmp >2.5 
table(row_sel)
anno_col_tmp = anno_col[col_sel, c("SC", "samp_class", "passage")]
anno_row_tmp = anno_row[row_sel, , drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_organoids_topSD.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = F, breaks = breaks,
        main='129 samples plus Organoid samples', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_colors = ann_colors)
dev.off()

#surface protein genes, blca 129 plus SHEN samples, no passage information, top SD
tmp = apply(tcga_shen_cc_log, 1, sd)
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT") | anno_col$Definition == "organoids" 
row_sel = gs %in% hsa_cd$ENTREZ.gene.symbol & tmp > 2.7
anno_col_tmp = anno_col[col_sel, c("SC", "Definition", "passage")]
anno_row_tmp = anno_row[row_sel, , drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
an = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
dim(matrix_tmp)
pdf("results/heatmap_tcga129_organoids_topSD_2.pdf", width=10, height=10)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks+1)), scale='none', show_rownames = T, breaks = breaks,
        main='129 samples plus Organoid samples', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_colors = ann_colors)
dev.off()


### blca plus organoids, both base47 and mcconkey genes
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT") | anno_col$Definition == "organoids"  & anno_col$del == 0
row_sel = anno_row$BASE47 != "NOT" | anno_row$McConkey != "NOT"
anno_col_tmp = anno_col[col_sel, c("SC", "samp_class")]
anno_row_tmp = anno_row[row_sel,c("BASE47", "McConkey")]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
matrix_tmp[1:4,1:4]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results2/heatmap_tcga129_organoids_base47_mcconkey.pdf", width=10, height=9)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks+1)), scale='none', show_rownames = T, breaks = breaks,
        main='129 samples plus Organoid samples', clustering_method = 'complete', border_color = NA, fontsize_col = 7,
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F,
        annotation_col = anno_col_tmp, annotation_colors = ann_colors, annotation_row = anno_row_tmp)
dev.off()

anno_col[434:474,]

## mcconkey, organoids
col_sel = anno_col$Definition == "organoids"  & anno_col$del == 0
row_sel = anno_row$McConkey != "NOT"
anno_col_tmp = anno_col[col_sel, c("samp_class", "passage")]
anno_row_tmp = anno_row[row_sel,c("McConkey"), drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
matrix_tmp[1:4,1:4]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
row.names(matrix_tmp) = sub(".*_", "", row.names(matrix_tmp))
row.names(anno_row_tmp) = sub(".*_", "", row.names(anno_row_tmp))
pdf("results/heatmap_organoids_mcconkey.pdf", width=7, height=4.5)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks+1)), scale='none', show_rownames = T, breaks = breaks,
        main='Organoid/tumor samples McConkey genes', clustering_method = 'ward.D2', border_color = NA, fontsize_col = 6,
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_colors = ann_colors, annotation_row = anno_row_tmp)
dev.off()

## blca plus organoids, both base47 and mcconkey genes, with colnames
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT") | anno_col$Definition == "organoids" & anno_col$del == 0
row_sel = anno_row$BASE47 != "NOT" | anno_row$McConkey != "NOT"
anno_col_tmp = anno_col[col_sel, c("SC", "samp_class")]
anno_row_tmp = anno_row[row_sel,c("BASE47", "McConkey")]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
matrix_tmp[1:4,1:4]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_organoids_base47_mcconkey_colname.pdf", width=24, height=9)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks+1)), scale='none', show_rownames = T, breaks = breaks,
        main='129 samples plus Organoid samples', clustering_method = 'complete', border_color = NA, 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 6, fontsize_col = 6,
        annotation_col = anno_col_tmp, annotation_colors = ann_colors, annotation_row = anno_row_tmp)
dev.off()

#### concensus clustering
#### concensus clustering
#BASE47 McConkey genes, blca plus SHEN samples
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT") | anno_col$Definition == "organoids" & anno_col$del == 0
row_sel = anno_row$BASE47 != "NOT" | anno_row$McConkey != "NOT"
anno_col_tmp = anno_col[col_sel, c("SC", "samp_class")]
anno_row_tmp = anno_row[row_sel,c("BASE47", "McConkey")]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
matrix_tmp[1:4,1:4]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
title=tempdir()
pdf("results/consensus.pdf")
results = ConsensusClusterPlus(matrix_tmp, maxK=5,reps=50,pItem=0.8,pFeature=1, 
			       title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279)
dev.off()
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
r5_order = r5$consensusTree[['order']]
r5_mat = r5$consensusMatrix
dim(r5_mat)
r5_mat_2 = r5_mat[r5_order, r5_order]
anno_col_tmp2 = anno_col_tmp[r5_order,] 
anno_col_tmp2$col = ann_colors[['samp_class']][anno_col_tmp2$samp_class]
pdf('results/b.pdf', width=10, height=3)
plot(NULL, xlim=c(0, nrow(anno_col_tmp2)), ylim=c(0,4))
nrow(anno_col_tmp2) -> nn
rect(1:(nn-1), rep(1, (nn-1)), 2:nn, rep(2, (nn-1)), col=anno_col_tmp2$col) 
dev.off()

### consencus clustering use top sd genes
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT") | anno_col$Definition == "organoids" 
sd = apply(tcga_shen_cc_log, 1, sd)
row_sel = sd > 2.5
anno_col_tmp = anno_col[col_sel, c("SC", "samp_class")]
anno_row_tmp = anno_row[row_sel,c("BASE47", "McConkey")]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
matrix_tmp[1:4,1:4]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
title=tempdir()
pdf("results/consensus.pdf")
results = ConsensusClusterPlus(matrix_tmp, maxK=5,reps=50,pItem=0.8,pFeature=1, 
			       title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279)
dev.off()
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
r5_order = r5$consensusTree[['order']]
r5_mat = r5$consensusMatrix
dim(r5_mat)
r5_mat_2 = r5_mat[r5_order, r5_order]
anno_col_tmp2 = anno_col_tmp[r5_order,] 
anno_col_tmp2$col = ann_colors[['samp_class']][anno_col_tmp2$samp_class]
pdf('results/b.pdf', width=10, height=3)
plot(NULL, xlim=c(0, nrow(anno_col_tmp2)), ylim=c(0,4))
nrow(anno_col_tmp2) -> nn
rect(1:(nn-1), rep(1, (nn-1)), 2:nn, rep(2, (nn-1)), col=anno_col_tmp2$col) 
dev.off()
#### concensus clustering
#### concensus clustering


#early117 genes, blca 129 
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT")
row_sel = anno_row$Early117 != "NOT"
#anno_col_tmp = anno_col[col_sel, c("nat_cl", "Tumor_Stage", "SC", "Definition", "BASE47_129")]
anno_col_tmp = anno_col[col_sel, c("nat_cl", "SC", "samp_class")]
anno_row_tmp = anno_row[row_sel, "Early117", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results2/heatmap_tcga129_arly117.pdf", width=10, height=11)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks,
        main='129 samples', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()

#early117 genes, blca 129 plus SHEN samples 
anno_col
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT") | anno_col$Definition == "organoids" 
row_sel = anno_row$Early117 != "NOT"
sum(row_sel)
anno_col_tmp = anno_col[col_sel, c("nat_cl", "SC", "samp_class")]
anno_row_tmp = anno_row[row_sel, "Early117", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
matrix_tmp
cn = colnames(matrix_tmp)
cn
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results2/heatmap_tcga129_organoids_early117.pdf", width=10, height=11)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks,
        main='129 samples plus our sample', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()

#early117 genes, SHEN samples 
col_sel = (anno_col$passage != 'TCGA') 
row_sel = anno_row$Early117 != "NOT"
anno_col_tmp = anno_col[col_sel, c("passage", "samp_class")]
anno_row_tmp = anno_row[row_sel, "Early117", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
head(matrix_tmp)
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
cn
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
row.names(matrix_tmp) = sub(".*_", "", row.names(matrix_tmp))
matrix_tmp
row.names(anno_row_tmp) = sub(".*_", "", row.names(anno_row_tmp))
anno_row_tmp
pdf("results2/heatmap_organoids_early117.pdf", width=10, height=12)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks,
        main='only our samples', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()
hp

#McConkey genes, blca 129 
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT")
row_sel = anno_row$McConkey != "NOT"
anno_col_tmp = anno_col[col_sel, c("nat_cl", "SC", "samp_class")]
anno_row_tmp = anno_row[row_sel, "McConkey", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_mcconkey.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks, border_color = NA,
        main='129 samples, classifiers from McConkey, 2015', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()


#McConkey genes, blca 129  organoids
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT") | anno_col$Definition == "organoids"  & anno_col$del == 0
row_sel = anno_row$McConkey != "NOT"
anno_col_tmp = anno_col[col_sel, c("nat_cl", "SC", "samp_class")]
anno_row_tmp = anno_row[row_sel, "McConkey", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_organoids_mcconkey.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks,, border_color = NA,
        main='129 samples, classifiers from McConkey, 2015', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()


#McConkey genes, blca 129 , organoids, with more heads
col_sel = (anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT") | anno_col$Definition == "organoids"  & anno_col$del == 0
row_sel = anno_row$McConkey != "NOT"
anno_col_tmp = anno_col[col_sel, c("nat_cl", "SC", "Definition", "samp_class", "passage")]
anno_row_tmp = anno_row[row_sel, "McConkey", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = sweep(matrix_tmp, 1, apply(matrix_tmp,1,median,na.rm=T))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_organoids_mcconkey_morehead.pdf", width=10, height=6)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks,
        main='129 samples, classifiers from McConkey, 2015', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()

#McConkey genes, blca plus SHEN samples , with sample names
col_sel = anno_col$Definition == "Primary solid Tumor" & anno_col$nat_cl != "NOT"| anno_col$Definition == "organoids" & anno_col$del == 0
row_sel = anno_row$McConkey != "NOT" & anno_row$sd 
anno_col_tmp = anno_col[col_sel, c("nat_cl",  "SC", "samp_class", "passage")]
anno_row_tmp = anno_row[row_sel, "McConkey", drop=F]
matrix_tmp = tcga_shen_cc_log[row_sel, col_sel]
cn = colnames(matrix_tmp)
matrix_tmp = t(apply(matrix_tmp, 1, scale))
colnames(matrix_tmp) = cn
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
pdf("results/heatmap_tcga129_organoids_mcconkey_colnames.pdf", width=30, height=10)
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', show_rownames = T, breaks = breaks,
        main='129 samples plus Organoid samples', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 10,
        annotation_col = anno_col_tmp, annotation_row = anno_row_tmp, annotation_colors = ann_colors)
dev.off()

tmp = cutree(hp$tree_col, 2)
tmp[tmp == 2] = "basal"
tmp[tmp == 1] = "luminal"
tmp = tmp[grep("MS0", names(tmp))]
mcconkey_org_cluster = tmp
head(mcconkey_org_cluster)

cv=data.frame(mcconkey_org_cluster, base47_org_cluster)
write.table(cv, file="results/base47_mcconkey_crossvalidation.txt")


## survive curve for these 129 samples McConkey clusters
#library(survival)
#
#death = blca.data$days_to_death
#death = death[anno_col$nat_cl != "NOT" & anno_col$nat_cl != "organoids"]
#followup = blca.data$days_to_last_follow_up
#followup = followup[anno_col$nat_cl != "NOT" & anno_col$nat_cl != "organoids"]
#event = rep(0, length(death))
#event[death > 0] = 1
#
#time = death
#time[is.na(time)] = followup[is.na(time)]
#tmp = anno_col[anno_col$nat_cl != "NOT" & anno_col$nat_cl != "organoids", "McConkey_129"]
#tmp
#
#
#sf = survfit(Surv(time, event) ~ tmp)
#options(repr.plot.width=4, repr.plot.height=4)
#plot(sf, col=2:5, mark.time=T, cex=.5, xpd=T, xlab="days", ylab="% survived")
#legend(4000, 1, col=2:5, legend=c("C1 n=52", "C2 n=76"), lty=1, cex=.5, bty="n", y.intersp=3)
#text(4000, .5, "P value .3", cex=.7)
##data.frame(death = death, followup = followup, time=time)
#sf$n
#survdiff(Surv(time, event) ~ tmp)
#
#
#
### only BASE47 genes
#yesno = anno_row$BASE47 != "NOT"
#anno_row[yesno > 0, 1, drop=F] -> anno_row_sel
#tcga_shen_cc_log[yesno > 0,] -> tcga_shen_cc_log_sel
#dim(anno_row_sel)
#row.names(tcga_shen_cc_log_sel) = blca.rows[row.names(tcga_shen_cc_log_sel), "external_gene_name"]
#row.names(anno_row_sel) = blca.rows[row.names(anno_row_sel), "external_gene_name"]
#png(file="heatmap_base47_all.png", res=200, height=1200, width=1500)
#pdf("results/heatmap_tcga129_base47.pdf", width=10, height=6)
#hp = pheatmap(tcga_shen_cc_log_sel, color=col, breaks = breaks, scale='row', show_rownames = T,
#        main='TCGA w/ Shen lab RNA-Seq', clustering_method = 'complete',
#        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
#        annotation_col = anno_col, annotation_row = anno_row_sel, annotation_colors = ann_colors)
#dev.off()
#dev.copy2png(file="heatmap_base47_all.pdf", res=200)
#system("open heatmap_base47_all.pdf")
#
#
#
### only mcconkey genes
#yesno = anno_row$McConkey != "NOT"
#anno_row[yesno, 2, drop=F] -> anno_row_sel
#tcga_shen_cc_log[yesno,] -> tcga_shen_cc_log_sel
#dim(anno_row_sel)
#apply(tcga_shen_cc_log_sel, 1, sd) -> sd
#tcga_shen_cc_log_sel = tcga_shen_cc_log_sel[sd > 0.5,]
#anno_row_sel = anno_row_sel[row.names(tcga_shen_cc_log_sel),, drop=F]
#
#tcga_shen_cc_log_sel = as.matrix(tcga_shen_cc_log_sel)
#row.names(tcga_shen_cc_log_sel) = blca.rows[row.names(tcga_shen_cc_log_sel), "external_gene_name"]
#anno_row_sel = as.matrix(anno_row_sel)
#row.names(anno_row_sel) = blca.rows[row.names(anno_row_sel), "external_gene_name"]
#
#png(file="heatmap_mcconkey_all.png", res=200, width=1500, height=600)
#pdf("results/heatmap_tcga129_base47.pdf", width=10, height=6)
#hp = pheatmap(tcga_shen_cc_log_sel, color=col, breaks = breaks, scale='row', show_rownames = T,
#        main='TCGA w/ Shen lab RNA-Seq', clustering_method = 'complete',
#        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
#        annotation_col = anno_col, annotation_row = anno_row_sel, annotation_colors = ann_colors)
#dev.off()
#dev.copy2pdf(file="heatmap_mcconkey_all.pdf")
#system("open heatmap_base47_all.pdf")
#
### only early 117 genes
#yesno = anno_row$Early117!= "NOT"
#anno_row[yesno, 3, drop=F] -> anno_row_sel
#tcga_shen_cc_log[yesno,] -> tcga_shen_cc_log_sel
#dim(anno_row_sel)
#tmp = apply(tcga_shen_cc_log_sel, 1, sd)
#tcga_shen_cc_log_sel[tmp > 0.5,] -> tcga_shen_cc_log_sel
#anno_row_sel = anno_row_sel[row.names(tcga_shen_cc_log_sel),, drop=F]
#
#tcga_shen_cc_log_sel = as.matrix(tcga_shen_cc_log_sel)
#row.names(tcga_shen_cc_log_sel) = blca.rows[row.names(tcga_shen_cc_log_sel), "external_gene_name"]
#sel = !duplicated(row.names(tcga_shen_cc_log_sel))
#tcga_shen_cc_log_sel = tcga_shen_cc_log_sel[sel,]
#anno_row_sel = anno_row_sel[sel,, drop=F]
#row.names(anno_row_sel) = blca.rows[row.names(anno_row_sel), "external_gene_name"]
#
#
#png(file="heatmap_early117_all_genenames.png", width=1200, res=200, height=3200)
#pdf("results/heatmap_tcga129_base47.pdf", width=10, height=6)
#hp = pheatmap(tcga_shen_cc_log_sel, color=col, breaks = breaks, scale='row', show_rownames = T,
#        main='TCGA w/ Shen lab RNA-Seq', clustering_method = 'complete',
#        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
#        annotation_col = anno_col, annotation_row = anno_row_sel, annotation_colors = ann_colors)
#dev.off()
#pdf(file="heatmap_early117_all_genenames.pdf", width=10, height=25)
#hp = pheatmap(tcga_shen_cc_log_sel, color=col, breaks = breaks, scale='row', show_rownames = T,
#        main='TCGA w/ Shen lab RNA-Seq', clustering_method = 'complete',
#        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F, fontsize_row = 6,
#        annotation_col = anno_col, annotation_row = anno_row_sel, annotation_colors = ann_colors)
#dev.off()
#
### selected primary, secondary, and Shen samples.
#sel = anno_col$type %in% c("primary", "secondary", "Shen")
#anno_col[sel, ] -> anno_col_sel
#tcga_shen_cc_log_sel[, sel] -> tcga_shen_cc_log_sel2
#hp = pheatmap(tcga_shen_cc_log_sel2, color=col, breaks = breaks, scale='row', show_rownames = F,
#        main='TCGA w/ Shen lab RNA-Seq', clustering_method = 'complete',
#        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F,
#        annotation_col = anno_col, annotation_row = anno_row_sel, annotation_colors = ann_colors)


### gene expression 
table(blca.data$definition)

colnames(rsem_shen$counts)
condition = c(blca.data$definition, rep("Org", 28), rep("Org Primary", 8), rep("Org", 5)); ##ncol(rsem_shen$counts)))
condition
design = data.frame(
        row.names       = colnames(tcga_shen_rsem),
        condition       = condition,
        libType         = rep("PE", ncol(tcga_shen_rsem)));

ddsmat = DESeqDataSetFromMatrix(countData = tcga_shen_rsem,
        colData = design,
        design = ~ condition);

dds.ds <- estimateSizeFactors(ddsmat);

dds <- DESeq(dds.ds, parallel=T);

load("res.RData")
head(res)

res = results(dds, contrast=c("condition", "Org", "Primary solid Tumor"))
res_2 = results(dds, contrast=c("condition", "Org", "Org Primary"))

res = res[order(res$pvalue),]
res$ensembl = sub("_.*", "", row.names(res))
res$symbol = blca.rows[res$ensembl, "external_gene_name"]
head(res)

res_2 = res_2[order(res_2$pvalue),]
res_2$ensembl = sub("_.*", "", row.names(res_2))
res_2$symbol = blca.rows[res_2$ensembl, "external_gene_name"]
head(res_2)
class(res_2)

res = as.data.frame(res)
res_2 = as.data.frame(res_2)
head(res_2)

write.table(res, file="results/res_org_to_blca.txt", sep="\t")
write.table(res_2, file="results/res_org_to_primary.txt", sep="\t")

source("~/program/fun/write_rnk.r")
write_rnk(res, file="results/res_shen_to_blca.rnk")
write_rnk(res_2, file="results/res_shen_to_blca_2.rnk")
run_gsea("res_shen_to_blca_2.rnk")

res = as.data.frame(res); res = as.data.table(res)
res_2 = as.data.frame(res_2); res_2 = as.data.table(res_2)
res$symbol[1:100]
res$threshold = as.factor(res$symbol %in% res$symbol[1:100])
res$logpvalue = -log(res$padj)
res$logpvalue[is.infinite(res$logpvalue)] = 750
res$symbol2 = rep(NULL, nrow(res))
res$threshold == "TRUE" -> ss
res$symbol2[ss] = res$symbol[ss]
pdf("results/plotma_x.pdf", width=10, height=10)
ggplot(data=res, aes(y=logpvalue, x=log2FoldChange, color = threshold, label = symbol2)) +
#ggplot(data=res, aes(y=-log(padj), x=log2FoldChange, color = threshold)) +
	geom_point(alpha=0.4, size=1.4) + 
	geom_text(check_overlap=F, vjust=1, hjust = 1, show.legend = F, angle = 90, size=2, color=1) +
	#geom_hline(aes(yintercept=0), color='blue', size=1.2) + 
	ylab("-log FDR") +
	ylim(c(0,800)) +
	xlab("Log2 Fold Change") +
	theme(axis.title.x = element_text(face = "bold", size = 15),
	        axis.text.x = element_text(face = "bold", size = 12)) +
	theme(axis.title.y = element_text(face = "bold", size = 15),
		axis.text.y = element_text(face = "bold", size = 12)) +
	scale_colour_discrete(name = "Top 100") +
	theme(legend.title = element_text(face = "bold", size = 15)) +
	theme(legend.text = element_text(size = 14))  
dev.off()

res$threshold = as.factor(res$padj < 0.001 & !is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 2)
res$logpvalue = -log(res$padj)
res$logpvalue[is.infinite(res$logpvalue)] = 750
res$symbol2 = rep(NULL, nrow(res))
res$threshold == "TRUE" & res$symbol %in% intersect(res_2$symbol[1:1000], res$symbol[1:1000]) -> ss
res$symbol2[ss] = res$symbol[ss]
pdf("results/plotma_shared.pdf", width=10, height=10)
ggplot(data=res, aes(y=logpvalue, x=log2FoldChange, color = threshold, label = symbol2)) +
#ggplot(data=res, aes(y=-log(padj), x=log2FoldChange, color = threshold)) +
	geom_point(alpha=0.4, size=1.4) + 
	geom_text(check_overlap=F, vjust=1, hjust = 1, show.legend = F, angle = 90, size=2, color=1) +
	#geom_hline(aes(yintercept=0), color='blue', size=1.2) + 
	ylab("-log FDR") +
	ylim(c(0,800)) +
	xlab("Log2 Fold Change") +
	theme(axis.title.x = element_text(face = "bold", size = 15),
	        axis.text.x = element_text(face = "bold", size = 12)) +
	theme(axis.title.y = element_text(face = "bold", size = 15),
		axis.text.y = element_text(face = "bold", size = 12)) +
	scale_colour_discrete(name = "p.adjusted < 0.001") +
	theme(legend.title = element_text(face = "bold", size = 15)) +
	theme(legend.text = element_text(size = 14))  
dev.off()

res_2$threshold = as.factor(res_2$padj < 0.001 & !is.na(res_2$log2FoldChange) & abs(res_2$log2FoldChange) > 2)
res_2$logpvalue = -log(res_2$padj)
res_2$logpvalue[is.infinite(res_2$logpvalue)] = 750
res_2$symbol2 = rep("", nrow(res_2))
res_2$threshold == "TRUE" & res_2$symbol %in% intersect(res_2$symbol[1:1000], res_2$symbol[1:1000]) -> ss
ss = 1:100
res_2$symbol2[ss] = res_2$symbol[ss]
pdf("results/plotma_x2.pdf", width=10, height=10)
ggplot(data=res_2, aes(y=logpvalue, x=log2FoldChange, color = threshold, label = symbol2)) +
#ggplot(data=res, aes(y=-log(padj), x=log2FoldChange, color = threshold)) +
	geom_point(alpha=0.4, size=1.4) + 
	geom_text(check_overlap=F, vjust=1, hjust = 1, show.legend = F, angle = 90, size=2, color=1) +
	#geom_hline(aes(yintercept=0), color='blue', size=1.2) + 
	ylab("-log FDR") +
	ylim(c(0,120)) +
	xlab("Log2 Fold Change") +
	theme(axis.title.x = element_text(face = "bold", size = 15),
	        axis.text.x = element_text(face = "bold", size = 12)) +
	theme(axis.title.y = element_text(face = "bold", size = 15),
		axis.text.y = element_text(face = "bold", size = 12)) +
	scale_colour_discrete(name = "p.adjusted < 0.001") +
	theme(legend.title = element_text(face = "bold", size = 15)) +
	theme(legend.text = element_text(size = 14))  
dev.off()

plot(res$log2FoldChange, -log(res$pvalue), xlab="log2FC", ylab="-log(pvalue)")
tmp = res[!is.na(res$padj) & res$padj < 0.001 & !is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 2, ]
tmp = res[1:200,]
points(tmp$log2FoldChange, -log(tmp$pvalue), col=2)

## heatmap of top genes 
res[!is.na(res$padj) & res$padj < 0.01 & !is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 2,] -> tmp
n=100
row.names(res)[1:n] -> top500
col_sel = anno_col$Definition == "organoids"
top500 = tcga_shen_cc_log[top500,col_sel]
remove_null_row(top500) -> top500
row.names(top500) = res$symbol[1:n]
cn = colnames(top500)
top500 = t(apply(top500, 1, scale))
colnames(top500) = cn
top500[top500 > 2] = 2
top500[top500 < -2] = -2
options(repr.plot.width=8, repr.plot.height=10)
pdf("results/heatmap_orgVsPrimary.pdf", width=7, height=9)
hp = pheatmap(top500, color=col, scale='none', show_rownames = T, show_colnames = T, fontsize_row = 6,
        main='organoids vs primary RNA-Seq, top 150 genes', clustering_method = 'ward.D', fontsize_col = 6,
        cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()

## heatmap of top genes org vs tcga
res_2[!is.na(res_2$padj) & res_2$padj < 0.01 & !is.na(res_2$log2FoldChange) & abs(res_2$log2FoldChange) > 2,] -> tmp
n=100
row.names(res_2)[1:n] -> top500
top500 = tcga_shen_cc_log[top500,]
row.names(top500) = res_2$symbol[1:n]
cn = colnames(top500)
top500 = t(apply(top500, 1, scale))
colnames(top500) = cn
top500[top500 > 2] = 2
top500[top500 < -2] = -2
options(repr.plot.width=8, repr.plot.height=10)
pdf("results/heatmap_shenVsTcga.pdf", width=10, height=9)
hp = pheatmap(top500, color=col, scale='none', show_rownames = T, show_colnames = F, fontsize_row = 6,
        main='organoids vs TCGA BLCA lab RNA-Seq, top 100 genes', clustering_method = 'complete', fontsize_col = 7,
        cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()

intersect(res_2$symbol[1:200], res$symbol[1:200])
res[res$symbol == 'CD34',]
####  END

### count number per samples comparison show Shen's data has much low reads
assays(blca.data) -> tcga_cc
tcga_cc = tcga_cc[[1]]
### counttable

pdf("hist.pdf", height=4, width=4)
dim(as.numeric(tcga_cc))
hist(as.matrix(tcga_cc[,3:435]), ylim=c(0, 3*100000), xlab="reads per gene", main="TCGA BLCA RNA-Seq samples")
hist(counttable, ylim=c(0, 1.5*10000000), xlab="reads per gene", main="Shen RNA-Seq samples")
dev.off()
system("open hist.pdf")
median(as.matrix(tcga_cc[,3:435]))
median(counttable)

breaks=seq(0, 20000, by=500)
breaks = c(breaks, 2000000)
pdf("hist2.pdf")
hist(tcga_cc[,3],  xlab="reads per gene", main="TCGA BLCA RNA-Seq samples", breaks = breaks, xlim=c(0, 10000), ylim=c(0, 0.002))
hist(counttable[6:nrow(counttable),1],  xlab="reads per gene", main="TCGA BLCA RNA-Seq samples", breaks = breaks, xlim=c(0, 10000), ylim=c(0, 0.002))
dev.off()
system("open hist2.pdf")

apply(tcga_cc, 2, mean)
apply(counttable, 2, mean)
apply(tcga_cc, 2, median)
apply(counttable, 2, median)

## fpkm
read.delim("Homo_sapiens.GRCh38.83.gtf_genelength.txt", header=F,row.names=1,  sep=" ") -> gl
head(gl)
gl$len = gl$V4 -gl$V3

tcga_shen[row.names(tcga_shen) %in% row.names(gl),] -> tcga_shen_protein

gl1 = gl[row.names(tcga_shen_protein), "len"]
sumOfReads = apply(tcga_shen_protein, 2, sum)
tcga_shen_protein_fpkm = tcga_shen_protein
for(i in 1:ncol(tcga_shen)){
        tcga_shen_protein_fpkm[,i] = (tcga_shen_protein[,i]/gl1) * (1000000000 / sumOfReads[i])
}
head(tcga_shen_protein_fpkm)


## use the fpkm draw the heatmap
gs = blca.rows[row.names(tcga_shen_protein_fpkm), "external_gene_name"]

anno_row = data.frame(row.names=row.names(tcga_shen_protein_fpkm), BASE47 = rep("NOT", nrow(tcga_shen_protein_fpkm)), stringsAsFactors = F);
anno_row$BASE47[gs %in% base47.luminal] = "Luminal";
anno_row$BASE47[gs %in% base47.basal] = "Basal";

anno_row$McConkey = "NOT"
anno_row$McConkey[gs %in% mcconkey.luminal] = "Luminal"
anno_row$McConkey[gs %in% mcconkey.basal] = "Basal"

anno_row$Early117 = "NOT"
anno_row$Early117[gs %in% early_uc_117_list[["class2"]] ] = "Class2"
anno_row$Early117[gs %in% early_uc_117_list[["class12"]] ] = "Class12"
anno_row$Early117[gs %in% early_uc_117_list[["class13"]] ] = "Class13"

apply(anno_row, 1, function(x){length(x[x!="NOT"])}) -> yesno
anno_row[yesno > 0,] -> anno_row_sel   
tcga_shen_protein_fpkm[yesno > 0,] -> tcga_shen_protein_fpkm_sel
dim(anno_row_sel)

pdf("heatmap_2.pdf", width=15, height=10)
hp = pheatmap(tcga_shen_protein_fpkm_sel, color=col, breaks = breaks, scale='row', show_rownames = F,
        main='TCGA w/ Shen lab RNA-Seq', clustering_method = 'complete',
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F,
        annotation_col = anno_col, annotation_row = anno_row_sel, annotation_colors = ann_colors)
dev.off()
                                                                                                                                                                          


filenames = c('MS001.counts.txt','MS002.counts.txt', 'MS003.counts.txt', 'MS004.counts.txt',
              'MS005.counts.txt', 'MS008.counts.txt', 'MS009.counts.txt', 'MS010.counts.txt',
              'MS011.counts.txt', 'MS013.counts.txt', 'MS014.counts.txt', 'MS015.counts.txt',
              'MS016.counts.txt', 'MS018.counts.txt')


if (exists("test")) {rm(test_count)}

for (i in 1:length(filenames)) {
        sample <- read.table(as.character(filenames[i]), stringsAsFactors = F, header=T);
        if (exists("test_count")) {test_count <- merge(test_count, sample, by= "gene")} else {test_count <- sample}
        rm(sample)
}



head(test_count)

row.names(test_count) = test_count$gene
test_count = test_count[,-1]
head(test_count)

condition = rep("tumor", times =ncol(test_count));
design = data.frame(
        row.names       = colnames(test),
        condition       = condition,
        libType         = rep("PE", ncol(test)));



tt = DESeqDataSetFromMatrix(countData = test,
        colData = design,
        design = ~ 1);

tt <- estimateSizeFactors(tt);

tt <- DESeq(tt, parallel=T);

tt = counts(tt, normalized=T);


tmp = apply(test_count, 2, sum)
options(repr.plot.width=3, repr.plot.height=10)
barplot(tmp, xlab="Total reads in \ncoding genes (orignal)", col=3, names.arg = "",  ylim=c(0, 1.5e+8))

mean(as.numeric(tt['FAP', ]))
mean(as.numeric(blca_cc['FAP', ]))


# exon sequence analysis
# import exon mutations
# Tumor_Sample_Barcode is comID
# need comID for impact samples
library(autospy)
library(data.table)
library(readxl)
setwd("exon")
getwd()                                # 
#check args.tsv: if info_file is exon_sampleID_v2.xlsx, MaB30 has two primary 
#check args.tsv: if info_file is exon_sampleID.xlsx, MaB30 has two primary 
# need change the exon_patientID.xlsx when switch
## pdf1 folder combine the patient MaB30 two primary samples
## see args.tsv
## exon_sampleID.xlsx is part of exon_patientID.xlsx
## current not used here
setwd("./pdf2")

IDs = read_excel("exon_patientID_v2.xlsx", 3)
IDs = as.data.table(IDs)
IDs[exon == 'yes',] -> IDs
IDs
IDs[, patientID2 := paste0("Pat", patientID)]
IDs[, sampleID2 := paste0("Sam", sampleID)]
setkey(IDs, seqID)
IDs$newid
IDs$seqID

## reorganize maf file
maffile = "/ifs/res/share/solit/solitd/Proj_06048_P/r_001/post/Proj_06048_P___SOMATIC.vep.filtered.facets.V3.maf"
maf = data.table::fread(maffile)
maf
unique(maf$Tumor_Sample_Barcode)
colnames(maf)
maf[, seqID := Tumor_Sample_Barcode] 
maf[, VAF := t_alt_count / t_depth] 
setkey(IDs, seqID)
maf[, patientID2 := IDs[maf$seqID, 'patientID2'] ]
maf[, sampleID2 := IDs[maf$seqID, 'sampleID2'] ]
maf[,old_barcode := Tumor_Sample_Barcode]
maf[, Tumor_Sample_Barcode := IDs[maf[, old_barcode], comID] ]
maf[, .(Tumor_Sample_Barcode, old_barcode)]
maf[, FILTER2 := FILTER]
maf[FILTER2 == "PASS", FILTER := "."]
maf[FILTER2 == "RESCUE", FILTER := "."]
maf[Chromosome == "MT", Chromosome := 'M']
maf[, vaf := 100 * t_alt_count / t_depth]
maf[,TAG2 := paste0(Hugo_Symbol, ":", TAG)]
maf[FILTER == '.',] -> maf
setkey(maf, Tumor_Sample_Barcode)
setkey(IDs, comID)
maf

## filter low vaf mutations
maf$vaf
maf
maf2 = maf[, .SD[max(ccf_Mcopies) > .1], by = list(patientID2, Hugo_Symbol)]
maf2$ccf_Mcopies

table(maf$Tumor_Sample_Barcode)
table(maf2$Tumor_Sample_Barcode)

#based on AgilentExon_51MB_b37_v3_targets.bed
coverage = 32931921
tmp = maf[, .(meanVAF = median(VAF), meanMut = length(VAF)/length(unique(sampleID2)) , mutNum = length(VAF), sampleNum = length(unique(sampleID2))), by = patientID2]
tmp[, perMb := 1000000*meanMut/coverage]
sum(tmp$meanMut) / 5

tmp2 = maf2[, .(meanVAF = median(VAF), meanMut = length(VAF)/length(unique(sampleID2)) , mutNum = length(VAF), sampleNum = length(unique(sampleID2))), by = patientID2]
tmp2[, perMb := 1000000*meanMut/coverage]
tmp2
sum(tmp2$meanMut) / 5

fwrite(maf2, file="Proj_6048.maf", sep="\t", quote=F)

process_autopsy_maf()

xp = read.maf("Proj_6048.maf")
xp = xp@data
xp[,.(Tumor_Sample_Barcode)]
maf[,.(Tumor_Sample_Barcode)]
intersect(xp[,.(TAG)], maf[,.(TAG)])
intersect(xp$TAG, maf$TAG)
ii = setdiff(maf$TAG, xp$TAG)
maf[TAG %in% ii,]
read.maf
xx = maf
maf="Proj_6048.maf"



nrow(maf)
maf$newid
## maf is clean filtered, maf2 has FILTER2 == "."
nrow(maf)
nrow(maf2)
xx = table(maf$Tumor_Sample_Barcode)
sum(as.numeric(xx))
table(maf$FILTER2)

pdf("../../results/test.pdf")
plot_mutation_overlap(maf, freq=0, log = F)
plot_mutation_overlap(maf, freq=1, log = F)
dev.off()

plot_mutation_overlap <- function(maf, pid = NULL, log = F, freq=0) {
	dc <- dcast.data.table(maf,
			       TAG ~ Tumor_Sample_Barcode,
			       value.var = "Tumor_Sample_Barcode",
			       fun.aggregate = function(x) {
				       as.integer(length(x)>0)
			       }, fill = 0)
	mat <- as.matrix(dc[, -1, with = F])
	### take the pair-wise overlap of columns (Number of mutations in common)
	cm <- crossprod(mat)
	diag(cm) <- 0
	cm[upper.tri(cm)] = 0
	fill_palette <- colorRampPalette(c("white", "red"))(n = 499)
	if(freq == 1){ 
		cm2 = cm
		tb = table(maf$Tumor_Sample_Barcode)
		cm3 = apply(cm2, 1, function(x){100*x/tb})
		cm = cm3
		fill_palette <- colorRampPalette(c("white", "blue"))(n = 99)
		key.xlab = "% overlap based on bottom sample"
	} else{
		if(log == TRUE) {
			cm <- log10(cm+1)
			key.xlab = "log10( Nmutations )"
		} else {
			key.xlab = "Nmutations"
		}
	}
	if(!is.null(pid)){
		ipatient <- as.integer(factor(stringr::str_extract(rownames(cm), pid)))
		side_palette <- gg_color_hue(uniqueN(ipatient))[ipatient]
		gplots::heatmap.2(cm, trace = "none", col = fill_palette,
				  Rowv = NULL, Colv = NULL, lhei = c(2, 8), margins=c(10,10),
				  key.title = NA, key.ylab = NA, key.xlab = key.xlab,
				  RowSideColors = side_palette
				  )
	} else {
		gplots::heatmap.2(cm, trace = "none", col = fill_palette, Rowv = NULL, Colv = NULL, lhei = c(2, 8), margins=c(10,10), key.title = NA, key.ylab = NA, key.xlab = key.xlab)
	}
}

library(ggplot2)
library(data.table)
library(methods)
library(stringdist)
library(ape)
library(phytools)
library(stringr)
library(RColorBrewer)
library(readxl)
library(gplots)
library(cowplot)


## see args.tsv
## exon_sampleID.xlsx is part of exon_patientID.xlsx

## the same but more clean, and used in the comparison between impact and exon sequence
#maffile = "exon_filtered.maf"
## this is for maftools
setwd("..")
exon_maf = read.maf(maf)
exon_maf_data = exon_maf@data
nrow(exon_maf_data)
maf

## IDs for impact seq
## all tumor_sample_barcode ID are converted to comID
setwd("exon")
impact_id = fread("impact_clinical_mixed_cbe_solitd_DS_blorg_clinical_datav2.tsv")
impact_id[, seqID := `Sample ID`]
tmp = impact_id$"Collaboration Id"
tmp = unlist(lapply(strsplit(tmp, "_"), "[[", 1))
impact_id[, patientID := tmp]
impact_id[, patientID2 := paste0("Pat", patientID)]
tmp = impact_id$"Collaboration Id"
tmp = sub(".*?_", "", tmp)
impact_id[, sampleID := tmp]
impact_id[, sampleID2 := paste0("Sam", sampleID)]
impact_id[,  comID := factor(paste0(patientID2, "_", sampleID2))]
impact_id[,  Tumor_Sample_Barcode := comID]
impact_id[grep("PDXOrgP0", Tumor_Sample_Barcode), ]
impact_id[, newid := factor(tmp, levels = c('Par_slides', 'Par_Pellet', 'Par_Slides', 'Par', 'Par_FFPE', 'Par_resection_A', 'Par_resection_B1', 'OrgP0', 'PDX474L_475RP6_OrgP0', 'Org', 'Organoid', 'OrgAE', 'OrgP1', 'OrgP3', 'OrgP5', 'OrgP6', 'OrgP7', 'OrgP10', 'OrgP11', 'OrgP14', 'OrgP8', 'OrgP4', 'OrgP13', 'OrgP2', 'OrgP9', 'PDX', 'PDX474LP6', 'PDX475RP6', 'PDXOrgP0', 'PDXOrgP1', 'PDXOrgP2', 'OrgP0PDX', 'P0_PDX', 'P1PDX'))]

newid$Line..old.name = sub(" ", "", newid$Line..old.name)
intersect(impact_id$patientID, newid$Line..old.name)
setdiff(impact_id$patientID, newid$Line..old.name)


## impact results for maftools
## import impact data
setwd("exon")
getwd()
read.maf("data_mutations_extended.maf") -> impact_maf
impact_maf_data = impact_maf@data
impact_maf_data[, seqID := as.character(Tumor_Sample_Barcode)]
impact_maf_data$Tumor_Sample_Barcode
setkey(impact_id, seqID)
setkey(impact_maf_data, Tumor_Sample_Barcode)
impact_maf_data$Tumor_Sample_Barcode   # 
impact_maf_data
#impact_maf_data = impact_maf_data[impact_id] ## add comID by merge
impact_maf_data = impact_maf_data[impact_id[,.(seqID, comID, patientID, patientID2, sampleID, sampleID2, newid)]]
impact_maf_data[, .(patientID, patientID2, sampleID, sampleID2, newid, Tumor_Sample_Barcode, seqID)]
impact_maf_data[, vaf := 100 * t_alt_count / (t_ref_count + t_alt_count)]
impact_maf_data$TAG = paste0(impact_maf_data$Chromosome, ":",  impact_maf_data$Start_Position, "-",  impact_maf_data$End_Position, ":", impact_maf_data$Reference_Allele, ":", impact_maf_data$Tumor_Seq_Allele2)
## set the Tumor_Sample_Barcode the same as in the exon-seq ID
impact_maf_data[, Tumor_Sample_Barcode2 := Tumor_Sample_Barcode]
impact_maf_data$Tumor_Sample_Barcode
impact_maf_data$comID
impact_maf_data[, Tumor_Sample_Barcode := comID]
impact_maf_data_orig = impact_maf_data
impact_maf_data = impact_maf_data[vaf > 0, ]
impact_maf = read.maf(impact_maf_data)
nrow(impact_maf_data) # 1224
unique(impact_maf_data$patientID)  # 25
nrow(impact_maf_data)/length(unique(impact_maf_data$patientID))  # 48.96
impact_maf_data$comID
impact_maf_data$newid
write.table(impact_maf_data, sep="\t", file="data_mutations_extended_filtered.maf")

##\subsection{overlapped between exon and impact}
exon_maf_data$Tumor_Sample_Barcode
impact_maf_data$Tumor_Sample_Barcode
impact_maf_data$newid
exon_impact_data = exon_maf_data[Tumor_Sample_Barcode %in% impact_maf_data$Tumor_Sample_Barcode, ]
exon_impact_maf = read.maf(exon_impact_data)
nrow(exon_impact_data) # 6877

## overlapped between impact and exon
impact_exon_maf_data = impact_maf_data[Tumor_Sample_Barcode %in% exon_maf_data$Tumor_Sample_Barcode, ]
unique(impact_exon_maf_data$Hugo_Symbol) # 80 genes
nrow(impact_exon_maf_data) # 318
impact_exon_maf = read.maf(impact_exon_maf_data)
unique(impact_exon_maf@data$comID) # 22
impact_exon_maf = read.maf(impact_exon_maf_data)
#### end of importing data

## analysis pdf2
getwd()
setwd("exon/pdf2")

## start for pdf1: combine MaB30 all samples
setwd("../pdf1")

IDs1 = read_excel("exon_patientID.xlsx", 3)
IDs1 = as.data.frame(IDs1)
IDs1 = as.data.table(IDs1)
IDs1[exon == 'yes',] -> IDs1
IDs1[, patientID2 := paste0("Pat", patientID)]
IDs1[, sampleID2 := paste0("Sam", sampleID)]
setkey(IDs1, seqID)
## reorganize maf file
maffile = "/ifs/res/share/solit/solitd/Proj_06048_P/r_001/post/Proj_06048_P___SOMATIC.vep.filtered.facets.V3.maf"
#maffile = "Proj_06048_P___SOMATIC.vep.filtered.facets.V3.maf"
maf = data.table::fread(maffile)
maf[,old_barcode := Tumor_Sample_Barcode]
setkey(IDs, seqID)
maf[, Tumor_Sample_Barcode := IDs[maf[, old_barcode], comID] ]
maf[, .(Tumor_Sample_Barcode, old_barcode)]
maf[, FILTER2 := FILTER]
maf[FILTER2 == "PASS", FILTER := "."]
maf[FILTER2 == "RESCUE", FILTER := "."]
fwrite(maf, file="Proj_6048.maf", sep="\t", quote=F)

## see args.tsv
## notice: the Proj_6048.maf was used according the args.tsv
process_autopsy_maf()
## end of pdf1

## SomaticSignature
maf_tmp = maf[Variant_Type == 'SNP',]
plot_mutation_spectrum(maf_tmp, group = maf_tmp$patientID, pdffile = "../results/combined_mutation_signature.pdf", ww=7, hh=5.3)

#plot_mutation_spectrum(maf_tmp, group = maf_tmp$Tumor_Sample_Barcode, pdffile = "per_sample_somaticSignature.pdf")
patIDs = unique(maf_tmp$patientID2)

for(i in 1:length(patIDs)){
	tmp = maf_tmp[patientID2 == patIDs[i],]
	group = tmp$Tumor_Sample_Barcode
	pdffile=paste0("per_sample_", patIDs[i], "_somaticSignature.pdf")
	plot_mutation_spectrum(tmp, group = group, pdffile = paste0("per_sample_", patIDs[i], "_somaticSignature.pdf"))
}
plot_mutation_spectrum(maf_tmp, pdffile = "per_patient_somaticSignature.pdf")

plot_mutation_spectrum = function(maf_tmp, group=maf_group, pdffile, genome = BSgenome.Hsapiens.UCSC.hg19, ww, hh){
	maf_gr = VRanges(
			 seqnames = paste0("chr", maf_tmp$Chromosome),
			 ranges = IRanges(start = maf_tmp$Start_Position, end = maf_tmp$End_Position),
			 ref = maf_tmp$Reference_Allele,
			 alt = maf_tmp$Tumor_Seq_Allele2,
			 sampleNames = maf_tmp$Tumor_Sample_Barcode,
			 seqinfo = seqinfo(genome),
			 group = group)
	maf_gr = mutationContext(maf_gr, BSgenome.Hsapiens.UCSC.hg19)
	gg = plotMutationSpectrum(maf_gr, "group", colorby = "alteration")
	ggsave(gg, file=pdffile, width=ww, height=hh)
}
## see /home/huw/program/autospy/R/process_autopsy_maf.R
## plot_mutation_signature
## trying to incorporate into autospy package


## for sciClone
library(sciClone)
library(clonevol)
library(fishplot)
## reorganize the cnv data
cnvfile = '/ifs/res/share/solit/solitd/Proj_06048_P/r_001/variants/copyNumber/facets/Proj_06048_P_facets_merge_hisens.seg'
cnvfile = 'Proj_06048_P_facets_merge_hisens.seg'
cnv = read.table( cnvfile, sep="\t", header=T, stringsAsFactors = F)
colnames(cnv) = c('ID', 'chr', 'start', 'stop', 'num_mark', 'segment_mean')
#seq(-2.5, 4, by = 1) -> ct
#ct[1] = -5
#cut(cnv$segment_mean, ct, labels = c(-2, -1, 0, 1, 2, 3)) -> cnv_cut
#cnv$segment_mean = cnv_cut
#rm(cnv_cut)
cnv = cnv[,-5]
cnv = cnv[!is.na(cnv$segment_mean), ]
setkey(IDs, cnvID)
cnv$patientID2 = IDs[cnv$ID, patientID2]
cnv$sampleID2 = IDs[cnv$ID, sampleID2]
cnv$ID = IDs[cnv$ID, comID]
head(cnv)

fwrite(cnv, file="Proj_6048_cnv.seg", sep="\t", quote=F)
rm(process_autopsy_maf)
head(cnv)
## Porj_6048_cnv.seg is listed in args.tsv

## use maf directly

maf_clone = maf[Variant_Type == "SNP", c("Tumor_Sample_Barcode", "Chromosome",
						   "Start_Position", "t_ref_count", "t_alt_count", "vaf", "patientID2")]
colnames(maf_clone) = c("Tumor_Sample_Barcode", "chr", "pos", "ref_reads", "var_reads", 'vaf', "patientID2")

patIDs = unique(IDs$patientID2)
patID = patIDs[3]
var_mut = unique(maf_clone[patientID2 == patID, Tumor_Sample_Barcode]) # 
var_mut2 = paste0("mut_", var_mut) # for names of variable
for(i in 1:length(var_mut)){
	tmp = as.data.frame(maf_clone[patientID2 == patID & Tumor_Sample_Barcode == var_mut[i],])
	tmp = tmp[,-1]
	tmp = tmp[,-7]
	assign(var_mut2[i], tmp)
	rm(tmp)
}
var_cnv = unique(cnv[cnv$patientID2 == patID, 'ID']) # 
var_cnv2 = paste0("cnv_", var_cnv)
for(i in 1:length(var_cnv)){
	tmp = as.data.frame(cnv[cnv$ID == var_cnv[i] & cnv$patientID2 == patID,])
	assign(var_cnv2[i], tmp)
	rm(tmp)
}
var_mut_list = mget(var_mut2)
var_cnv_list = mget(var_cnv2)
names(var_mut_list)
names(var_cnv_list)

sc = sciClone(vafs = var_mut_list, copyNumberCalls = var_cnv_list,  sampleNames = var_mut, verbose = F)
sc = sciClone(vafs = var_mut_list, sampleNames = var_mut, verbose = F)
sc_df = sc@vafs.merged
sc_sel = sc_df[!is.na(sc_df$cluster), ] 
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
##t_alt_count / t_depth
##plot with fishplot
##pdf('fish.pdf', width=8, height=5)
for (i in 1:length(fishes)){
	fish = layoutClones(fishes[[i]])
	fish = setCol(fish,f$clonevol.clone.colors)
	fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
		 vlines=seq(1, length(names)), vlab=names, pad.left=0.5)
}
#dev <- dev.off()

## SNPs per patient
snp = maf[Variant_Type == 'SNP', ]
tmp = snp[, .(patientID2, sampleID2)]
tmp = tmp[!duplicated(tmp),]
t1 = table(tmp$patientID2)
t2 = table(snp$patientID2)
t2/t1

### pdf2 for point line plot
## same as before, dont run
setwd("pdf2")

newid_conv = data.frame(row.names=c('SamPar',  'SamOrgP1', 'SamOrgP2', 'SamOrgP4', 'SamOrgP6', 'SamOrgP8', 'SamOrgP9', 'SamOrgP10','SamPDX', 'SamPDXOrgP0', 'SamPDXOrgP2'), newid=c("tumor", "orgP1", "orgP2", "orgP4", "orgP6", "orgP8", "orgP9", "orgP10", "xen", "xen-orgp0", "xen-orgP2"))
maf[, newid := factor(newid_conv[sampleID2, "newid"], levels=c("tumor", "orgP1", "orgP2", "orgP4", "orgP6", "orgP8", "orgP9", "orgP10", "xen", "xen-orgp0", "xen-orgP2"))]
maf$newid

maf_sel = maf[IMPACT_410  == T,]
table(maf_sel$FILTER)
maf_sel[, N := nrow(.SD), by = list(patientID2, TAG2)]
maf_sel = maf_sel[N > 1,]
gg = ggplot(maf_sel, aes(x=newid, y=vaf, group = patientID2, color = patientID2)) +
	geom_point() +
	geom_line() + 
	facet_wrap( ~ paste0(SYMBOL, ":", HGVSp_Short), ncol=8) + 
	#facet_wrap( ~ TAG2, ncol=8) + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        theme(strip.text.x = element_text(size = 8, hjust=0))
ggsave(gg, file=paste0("results/all_Pat_IMPACT_dynamic", ".pdf"), width=16, height=22)

## vaf, impact genes
getwd()
tmp = maf_sel ## save for pat1 pat2 ...
tmp$newid
patIDs = unique(tmp[, patientID2])
for(i in 1:length(patIDs)){
	patID = patIDs[i]
	maf_sel2 = tmp[patientID2 == patID,]
	maf_sel2[, N := nrow(.SD), by = list(patientID2, TAG2)]
	maf_sel2 = maf_sel2[N > 1,];
	maf_sel2[, newid:= factor(newid, levels=levels(droplevels(newid)))]
	nn = length(unique(maf_sel2$TAG2))
	hh = 2; ww = 2; h_adj = 1.2; w_adj = 1
	w = ww*8 + w_adj; h = hh * ceiling(nn/8) + h_adj; 
	if(nn < 8) {w = ww * nn + w_adj; } 
	gg = ggplot(maf_sel2, aes(x=newid, y=vaf, group = patientID2, color=patientID2)) +
		geom_point() +
		geom_line() + 
		facet_wrap( ~ paste0(SYMBOL, ":", HGVSp_Short), ncol=8) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
		theme(strip.text.x = element_text(size = 8, hjust=0))
	ggsave(gg, file=paste0("results/pat", i, "_IMPACT_", patID, ".pdf"), width=w, height=h)
}

## non impact genes
maf_sel = maf[IMPACT_410  == F,]
maf_sel[, sampIDFull := factor(sampleID2, levels=c('SamPar', 'SamPDX', 'SamPDXOrgP0', 'SamPDXOrgP2', 
					       'SamOrgP1', 'SamOrgP2', 'SamOrgP4', 'SamOrgP6', 'SamOrgP8', 'SamOrgP9', 'SamOrgP10'))]
tmp = maf_sel ## save for pat1 pat2 ...
for(i in 1:length(patIDs)){
	patID = patIDs[i]
	maf_sel = tmp[patientID2 == patID,]
	maf_sel[, N := nrow(.SD), by = list(patientID2, TAG2)]
	maf_sel = maf_sel[N > 1,];
	maf_sel[, sampID := factor(sampIDFull, levels=levels(droplevels(sampIDFull)))]
	nn = length(unique(maf_sel$TAG2))
	hh = 2; ww = 2; h_adj = 1.2; w_adj = 1
	w = ww*8 + w_adj; h = hh * ceiling(nn/8) + h_adj; 
	if(nn < 8) {w = ww * nn + w_adj; } 
	gg = ggplot(maf_sel, aes(x=sampID, y=vaf, group = patientID2, color=patientID2)) +
		geom_point() +
		geom_line() + 
		facet_wrap( ~ paste0(SYMBOL, ":", HGVSp_Short), ncol=8) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
		theme(strip.text.x = element_text(size = 8, hjust=0))
	ggsave(gg, file=paste0("pat", i, "_nonIMPACT_", patID, ".pdf"), width=w, height=h, limitsize = FALSE)
}

## ccf_Mcopies, impact genes
maf_sel = maf[IMPACT_410  == T,]
maf_sel[, sampIDFull := factor(sampleID2, levels=c('SamPar', 'SamPDX', 'SamPDXOrgP0', 'SamPDXOrgP2', 
					       'SamOrgP1', 'SamOrgP2', 'SamOrgP4', 'SamOrgP6', 'SamOrgP8', 'SamOrgP9', 'SamOrgP10'))]
tmp = maf_sel ## save for pat1 pat2 ...
patIDs = unique(tmp[, patientID2])
for(i in 1:length(patIDs)){
	patID = patIDs[i]
	maf_sel = tmp[patientID2 == patID,]
	maf_sel[, N := nrow(.SD), by = list(patientID2, TAG2)]
	maf_sel = maf_sel[N > 1,];
	maf_sel[, sampID := factor(sampIDFull, levels=levels(droplevels(sampIDFull)))]
	nn = length(unique(maf_sel$TAG2))
	hh = 2; ww = 2; h_adj = 1.2; w_adj = 1
	w = ww*8 + w_adj; h = hh * ceiling(nn/8) + h_adj; 
	if(nn < 8) {w = ww * nn + w_adj; } 
	gg = ggplot(maf_sel, aes(x=sampID, y=ccf_Mcopies, group = patientID2, color=patientID2)) +
		geom_point() +
		geom_line() + 
		facet_wrap( ~ paste0(SYMBOL, ":", HGVSp_Short), ncol=8) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
		theme(strip.text.x = element_text(size = 8, hjust=0))
	ggsave(gg, file=paste0("pat", i, "_ccf_IMPACT_", patID, ".pdf"), width=w, height=h, limitsize = F)
}

## ccf non impact genes
maf_sel = maf[IMPACT_410  == F,]
maf_sel[, sampIDFull := factor(sampleID2, levels=c('SamPar', 'SamPDX', 'SamPDXOrgP0', 'SamPDXOrgP2', 
					       'SamOrgP1', 'SamOrgP2', 'SamOrgP4', 'SamOrgP6', 'SamOrgP8', 'SamOrgP9', 'SamOrgP10'))]
for(i in 1:length(patIDs)){
	patID = patIDs[i]
	maf_sel = tmp[patientID2 == patID,]
	maf_sel[, N := nrow(.SD), by = list(patientID2, TAG2)]
	maf_sel = maf_sel[N > 1,];
	maf_sel[, sampID := factor(sampIDFull, levels=levels(droplevels(sampIDFull)))]
	nn = length(unique(maf_sel$TAG2))
	hh = 2; ww = 2; h_adj = 1.2; w_adj = 1
	w = ww*8 + w_adj; h = hh * ceiling(nn/8) + h_adj; 
	if(nn < 8) {w = ww * nn + w_adj; } 
	gg = ggplot(maf_sel, aes(x=sampID, y=ccf_Mcopies, group = patientID2, color=patientID2)) +
		geom_point() +
		geom_line() + 
		facet_wrap( ~ paste0(SYMBOL, ":", HGVSp_Short), ncol=8) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
		theme(strip.text.x = element_text(size = 8, hjust=0))
	ggsave(gg, file=paste0("pat", i, "_ccf_nonIMPACT_", patID, ".pdf"), width=w, height=h, limitsize = FALSE)
}

## description summary
# 24 samples from 4 patients. 5 are parental samples, and 12 parental derived organoid samples, 4 PDX samples, 3 PDX derived samples. 
# maf: filtered
# maf2 original with FILTER2 == "." for filter
dim(maf) # 27753 165, 27753 mutations SNPs
unique(maf$Hugo_Symbol) # 6446 genes passed the filter
unique(maf$patientID2) #  "PatSuB2"  "PatJuB3"  "PatMaB30" "PatMaB33" for 4 patients
table(maf$patientID2, maf$sampleID2) 
setwd("exon")

plotmafSummary(exon_maf, width=12, file="../results/mafSummary_exon.pdf")
plotmafSummary(impact_maf, file="../results/mafSummary_impact.pdf", width=12)

pdf(file="../results/oncoplot_exon_impact.pdf", width=12)
oncoplot(maf = exon_maf, genes=c('TP53', 'RB1', 'STAG2', 'KDM6A'), removeNonMutated = TRUE) 
oncoplot(maf = impact_maf, genes=c('TP53', 'RB1', 'STAG2', 'KDM6A'), removeNonMutated = TRUE) 
dev.off()

titv = titv(exon_maf, useSyn = T, file="../results/mafTitv_exon.pdf")
titv = titv(impact_maf, useSyn = T, file="../results/mafTitv_impact.pdf")

pdf("../results/mafLolliplot_exon.pdf", width=9, height=4)
lollipopPlot(exon_maf, gene="TP53", AACol = 'HGVSp', showMutationRate = TRUE, domainLabelSize = 3, defaultYaxis = FALSE, repel = TRUE)
lollipopPlot(exon_maf, gene="KDM6A", AACol = 'HGVSp', showMutationRate = TRUE, domainLabelSize = 3, defaultYaxis = FALSE, repel = T)
dev.off()

#lollipopPlot(exon_maf, gene="STAG2", AACol = 'HGVSp', showMutationRate = TRUE, domainLabelSize = 3, defaultYaxis = FALSE, repel = T)
#lollipopPlot(exon_maf, gene="RB1", AACol = 'HGVSp', showMutationRate = TRUE, domainLabelSize = 3, defaultYaxis = FALSE, repel = T)

##\section{analyze by maftools}
## read impact seq maf
## maftools impact seq
table(impact_maf_data$exac_filter)
length(unique(impact_maf_data$Hugo_Symbol))
unique(impact_maf_data$Hugo_Symbol)
unique(impact_maf_data[,Variant_Classification])
getSampleSummary(impact_maf)
getGeneSummary(impact_maf)
getFields(impact_maf)

pdf(file="../results/impact_summary.pdf", width=10, height=6)
plotmafSummary(maf = impact_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()

pdf(file="../results/exon_summary.pdf", width=10, height=6)
plotmafSummary(maf = exon_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()

impact_anno = as.data.frame(impact_id[, .(Tumor_Sample_Barcode, patientID)])
exon_anno = as.data.frame(IDs[, .(Tumor_Sample_Barcode, patientID)])
annoColor=list(patientID=c("JuB3" = 'red', "MaB30" = 'violetred', "SuB2" = 'green', "MaB33" = 'yellow', "MaBX30" = 'blueviolet'))


pdf(width=10, height=15, file="../results/impact_oncoplot.pdf")
#oncoplot(maf = impact_maf, top = 20, removeNonMutated = TRUE, annotation=impact_anno, showTumorSampleBarcode = F, fontSize = 8)
oncoplot(maf = impact_maf)#, top = 20, removeNonMutated = TRUE, annotation=impact_anno, showTumorSampleBarcode = F, fontSize = 8)
dev.off()

impact_id$patientID
impact_maf_data$Tumor_Sample_Barcode


tmp = read.table("a")
colnames(tmp) = 'oldid'
tmp$newid = sub("_.*", "", tmp$oldid)
tmp$newid
	 
exon_maf
pdf(width=6.5, height=7, file="../results/exon_oncoplot.pdf")
oncoplot(maf = exon_maf, removeNonMutated = TRUE, annotation=exon_anno, showTumorSampleBarcode = F, fontSize = 8, annotationColor = annoColor)
dev.off()

## plot impact mutations from patients that also have exome sequence
## see import data
pdf(file="../results/impact_exon_oncoplot_top20.pdf", width=6.5, height=6)
oncoplot(maf = impact_exon_maf, top = 20, removeNonMutated = TRUE, 
	 annotation=exon_anno, showTumorSampleBarcode = T, annotationColor=annoColor,
	 fontSize = 8)
dev.off()
pdf(file="../results/impact_exon_oncoplot_top80.pdf", width=6.5, height=16)
oncoplot(maf = impact_exon_maf, top = 80, removeNonMutated = TRUE, 
	 annotation=exon_anno, showTumorSampleBarcode = T, annotationColor=annoColor,
	 fontSize = 8)
dev.off()

## 
pdf(file="../results/exon_impact_oncoplot_top20.pdf", width=6.5, height=6)
oncoplot(maf = exon_impact_maf, top = 20, removeNonMutated = TRUE, 
	 annotation=exon_anno, showTumorSampleBarcode = T, annotationColor=annoColor,
	 fontSize = 8)
dev.off()
pdf(file="../results/exon_impact_oncoplot_top80.pdf", width=6.5, height=16)
oncoplot(maf = exon_impact_maf, top = 80, removeNonMutated = TRUE, 
	 annotation=exon_anno, showTumorSampleBarcode = T, annotationColor=annoColor,
	 fontSize = 8)
dev.off()
## vaf plot
source("~/program/fun/merge_maf.r")
source("~/program/fun/gg_color_hue.r")
xy = merge_maf(impact_maf_data, exon_maf_data)
xy$x[is.na(xy$x)] = 0
xy$y[is.na(xy$y)] = 0
xy

xy[Tumor_Sample_Barcode %in% exon_maf_data$Tumor_Sample_Barcode, ] -> xy
myColors = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(myColors) = c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'Missense_Mutation', 'Splice_Site', 'Translation_Start_Site', 'Nonsense_Mutation', 'In_Frame_Ins')

pdf(file="../results/exon_impact_vaf.pdf", width=10, height=8)
gg = ggplot(data=xy, aes(x = x, y = y, label=Hugo_Symbol)) +
	geom_point(aes(col = Variant_Classification), alpha = 0.8, size=1) +
	labs(x = "VAF (IMPACT-Seq)", y = "VAF (Exome-Seq)") + 
	facet_wrap(facets = ~ Tumor_Sample_Barcode, ncol = 5) +
	geom_text(size = 1.3, hjust = 1.2, check_overlap = T, nudge_y = 1) 
gg$labels$colour = "Variant Classification"
plot(gg)
dev.off()

xy[x > 0.0001 & y > 0.0001,] -> xy
pdf(file="../results/exon_impact_vaf_shared.pdf", width=8, height=9.6)
gg = ggplot(data=xy, aes(x = x, y = y, label=Hugo_Symbol)) +
	geom_point(aes(col = Variant_Classification), alpha = 0.8, size=1) +
	labs(x = "VAF (IMPACT-Seq)", y = "VAF (Exome-Seq)") + 
	facet_wrap(facets = ~ Tumor_Sample_Barcode, ncol = 4) +
	geom_text(size = 1.3, hjust = 1.2, check_overlap = T, nudge_y = 1)  + 
	scale_color_manual(values = myColors)
gg$labels$colour = "Variant Classification"
plot(gg)
dev.off()

cor.test(xy$x, xy$y)
cor(xy$x, xy$y)
pdf(file="../results/exon_impact_vaf_shared_combined.pdf", width=1.3*5.5, height=3.7*1.3)
gg = ggplot(data=xy, aes(x = x, y = y, label=Hugo_Symbol)) +
	geom_point(aes(col = Variant_Classification), alpha = 0.8, size=0.6) +
	labs(x = "VAF (IMPACT-Seq)", y = "VAF (Exome-Seq)") + 
	geom_text(size = 1.3, hjust = 1.2, check_overlap = T, nudge_y = 1) +
	scale_color_manual(values = myColors)
gg$labels$colour = "Variant Classification"
plot(gg)
dev.off()

col.sel = c("Tumor_Sample_Barcode", "Hugo_Symbol", "vaf", 'HGVSp', 'HGVSc', 'HGVSp_Short')
impact_maf_data[Hugo_Symbol == 'NF1' & patientID == 'JuB3', col.sel, with=F]
exon_maf_data[Hugo_Symbol == 'NF1' & patientID == 'JuB3', col.sel, with=F]
impact_maf_data[Hugo_Symbol == 'ERBB2' & patientID == 'MaBe3', col.sel, with=F]
exon_maf_data[Hugo_Symbol == 'ERBB3', col.sel, with=F]
impact_maf_data[Hugo_Symbol == 'TP53', col.sel, with=F]

oncoplot(maf = exon_impact_maf, top = 80, removeNonMutated = TRUE, 
	 annotation=exon_anno, showTumorSampleBarcode = T, annotationColor=annoColor, fontSize = 8)

impact_maf_data[Hugo_Symbol == 'CTNNB1', col, with=F]
impact_maf_data_orig[Hugo_Symbol == 'CTNNB1', ]
exon_maf_data[Hugo_Symbol == 'CTNNB1', col, with=F]

# from exom sequence, 78 genes from these 4 patients are impact_410 genes
# 75 of them are also detected by impact_410 sequence except E2F3, TERT, EGFR which have low VAF 0.021 - 0.036% in exome seqeunce
unique(exon_maf_data[Tumor_Sample_Barcode %in% impact_maf_data$Tumor_Sample_Barcode & IMPACT_410 == T, 
       Hugo_Symbol])  # 78
intersect(unique(exon_maf_data[Tumor_Sample_Barcode %in% impact_maf_data$Tumor_Sample_Barcode & IMPACT_410 == T,
	Hugo_Symbol]), unique(impact_exon_maf_data[, Hugo_Symbol]))
setdiff(unique(exon_maf_data[Tumor_Sample_Barcode %in% impact_maf_data$Tumor_Sample_Barcode & IMPACT_410 == T,
	Hugo_Symbol]), unique(impact_exon_maf_data[, Hugo_Symbol]))

impact_maf_data[Hugo_Symbol == 'E2F3',]
impact_maf_data[Hugo_Symbol == 'EGFR',]
impact_maf_data[Hugo_Symbol == 'TERT',]
col = c('seqID', 'Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome', 't_alt_count', 't_ref_count', 'vaf')
exon_maf_data[Hugo_Symbol %in% c('TERT', 'E2F3', 'EGFR'), col, with=F]

nrow(exon_maf_data[t_ref_count > 1000,])
print(exon_maf_data[vaf > 80, col, with=F], topn=100)
nrow(exon_maf_data[vaf > 80, col, with=F]) # 294

pdf(file="../results/hist.pdf")
hist(exon_maf_data$t_ref_count, main="tumor ref count in exome sequence")
hist(exon_maf_data$vaf, main="tumor vaf in exome sequence")
dev.off()

exon_impact_maf_data = exon_maf_data[Tumor_Sample_Barcode %in% impact_maf_data$Tumor_Sample_Barcode & IMPACT_410 == T, ]
unique(exon_impact_maf_data$Hugo_Symbol)
nrow(impact_exon_maf_data)
impact_exon_maf = read.maf(impact_exon_maf_data)
anno = as.data.frame(IDs[, .(Tumor_Sample_Barcode, patientID)])
annoColor=list(patientID=c("JuB3" = 'red', "MaB30" = 'violetred', "SuB2" = 'green', "MaB33" = 'yellow', "MaBX30" = 'mediumvioletred'))
length(unique(impact_exon_maf@data$Hugo_Symbol)) # 80
unique(impact_exon_maf@data$comID)
pdf(file="../results/impact_exon_oncoplot.pdf", width=6.5, height=16)
oncoplot(maf = impact_exon_maf, top = 80, removeNonMutated = TRUE, 
	 annotation=anno, showTumorSampleBarcode = T, annotationColor=annoColor,
	 fontSize = 8)
dev.off()

exon_maf_data[IMPACT_410 == T, .(Hugo_Symbol, Tumor_Sample_Barcode)]
tmp = read.maf(exon_maf_data[IMPACT_410 == T,])
tmp@data[,.(Tumor_Sample_Barcode, Hugo_Symbol)]
pdf(file="../results/exon_oncoplot_impact.pdf", width=10, height=9)
oncoplot(maf = tmp, removeNonMutated = TRUE, annotation=anno, showTumorSampleBarcode = T, annotationColor=annoColor)
dev.off()

titv = titv(maf = impact_maf, plot = FALSE, useSyn = TRUE)
pdf("impact_titv.pdf")
plotTiTv(res = titv)
dev.off()

pdf("impact_lolli_kdm6a.pdf")
lolli_kdm6a= lollipopPlot(maf = impact_maf, gene = 'KDM6A', AACol = 'Protein_Change', 
			  showMutationRate = TRUE, domainLabelSize = 3, defaultYaxis = FALSE)
dev.off()

## kit.lpop = lollipopPlot(maf = laml, gene = 'KIT', AACol = 'Protein_Change', labelPos = c(416, 418), refSeqID = 'NM_000222', domainLabelSize = 3)
## kit.lpop = lollipopPlot(maf = laml, gene = 'KIT', AACol = 'Protein_Change', labelPos = c(416, 418), refSeqID = 'NM_000222', repel = TRUE, domainLabelSize = 3)
## laml.dnmt3a = lollipopPlot(maf = laml, gene = 'DNMT3A', AACol = 'Protein_Change', refSeqID = 'NM_175629', labelPos = 882, collapsePosLabel = TRUE, cBioPortal = TRUE, domainLabelSize = 3, defaultYaxis = FALSE)


## coad.rf = rainfallPlot(maf = coad, detectChangePoints = TRUE, fontSize = 12, pointSize = 0.6)

## compared to TCGA
##laml.mutload = tcgaCompare(maf = laml, cohortName = 'Example-LAML')

# /ifs/work/solitlab/huw/solit/study/hiseq/tcga/tcga_maf_summary.txt.gz
source("../tcgacompare.R")
getwd()
maffile = "exon/Proj_6048.maf"
nrow(exon_maf_data)
tcgaCompare(maf = exon_maf, cohortName = "Organoids", fn="../results/exon_tcgaCompare.pdf")
tcgaCompare

tmp = exon_maf_data
setkey(tmp, Tumor_Sample_Barcode)
merge(tmp, IDs[, .(seqID, comID, patientID)], all.x=T, by.x='Tumor_Sample_Barcode', by.y='seqID') -> tmp2
anno = unique(tmp2[, .(Tumor_Sample_Barcode, patientID)])
anno = as.data.frame(anno)
anno
annoColor=list(patientID=c("JuB3" = 'red', "MaB30" = 'violetred', "SuB2" = 'green', "MaB33" = 'yellow'))
tmp2 = read.maf(tmp2)

pdf(file="../results/exon_oncoplot.pdf", width=10, height=5)
oncoplot(maf = tmp2, top = 20, removeNonMutated = TRUE, annotation=anno, annotationColor=annoColor, annotationColor = annoColor)
dev.off()

tmp = exon_maf_data
tmp = tmp[IMPACT_410 == T,]
setkey(tmp, Tumor_Sample_Barcode)
merge(tmp, IDs[, .(seqID, comID, patientID)], all.x=T, by.x='Tumor_Sample_Barcode', by.y='seqID') -> tmp2
anno = unique(tmp2[, .(Tumor_Sample_Barcode, patientID)])
anno = as.data.frame(anno)
anno
annoColor=list(patientID=c("JuB3" = 'red', "MaB30" = 'violetred', "SuB2" = 'green', "MaB33" = 'yellow'))
tmp2 = read.maf(tmp2)

pdf(file="../results/exon_oncoplot_impact.pdf", width=10, height=5)
oncoplot(maf = tmp2, top = 10, removeNonMutated = TRUE, annotation=anno, annotationColor=annoColor)
dev.off()
tmp2


## regenerate facets results of copy number data
IDs$facetsDir = sub("_hisens", "_facets", IDs$cnvID)
IDs$facetsDir = paste0("./r_001/variants/copyNumber/facets/", IDs$facetsDir)
IDs$facetsData = paste0(IDs$facetsDir, "/", IDs$cnvID, ".Rdata")
IDs$facetsPng = paste0(IDs$facetsDir, "/", IDs$cnvID, ".CNCF.png")
getwd()
source("~/program/facets-suite-1.0.1/fPlots2.R")
IDs$facetsData
i = 16
cnv_Rdata = IDs$facetsData[i]
cnv_Rdata
(load(cnv_Rdata))
names(out)
names(fit)
tmp = out$IGV
head(tmp)
tmp[tmp$chrom==1 & tmp$seg.mean > 1,]
head(fit$ploidy)
png(paste0(IDs$facetsDir[i], "/facets_custom2.png"), width=7*72, height=14*72)
plotSampleCNCF.custom.png(fit = fit, out = out$out, jointseg = out$jointseg, main = IDs$sampleID2[i])
dev.off()
pdf(paste0(IDs$facetsDir[i], "/facets_custom2.pdf"), width=7, height=14)
source("/Volumes/LaCie/huw/program/facets-suite-1.0.1/fPlots2.R")
plotSampleCNCF.custom.pdf(fit = fit, out = out$out, jointseg = out$jointseg, main = IDs$sampleID2[i])
dev.off()

source("/Volumes/LaCie/huw/program/facets-suite-1.0.1/fPlots_ggplot2v2.R")
for(i in 1:nrow(IDs)){
	cnv_Rdata = IDs$facetsData[i]
	load(cnv_Rdata)
	#plot.facets.all.output(out, fit, type='pdf', main=IDs$sampleID2[i], plotname=paste0(IDs$facetsDir[i], "/facets.pdf"))
	plot.facets.all.output(out, fit, type='png', main=IDs$sampleID2[i], plotname=paste0(IDs$facetsDir[i], "/facets.png"))
}

for(i in 1:nrow(IDs)){
	cnv_Rdata = IDs$facetsData[i]
	load(cnv_Rdata)
	#plot.facets.all.output(out, fit, type='pdf', main=IDs$sampleID2[i], plotname=paste0(IDs$facetsDir[i], "/facets.pdf"))
	plot.facets.all.output(out, fit, type='png', main=IDs$sampleID2[i], plotname=paste0(IDs$facetsDir[i], "/facets.png"))
}

library(imager)

for(i in 1:nrow(IDs)){
	png = paste0(IDs$facetsPng[i])
	load.image(png) -> im
	a = extract_patches(im, 423, 100, 842, 168)
	save.image(a[[1]], file=paste0(IDs$facetsDir[i], "/facets_custom_1.png"), quality=1)
	a = extract_patches(im, 423, 280, 842, 168)
	save.image(a[[1]], file=paste0(IDs$facetsDir[i], "/facets_custom_2.png"), quality=1)
	a = extract_patches(im, 423, 822, 842, 168)
	save.image(a[[1]], file=paste0(IDs$facetsDir[i], "/facets_custom_5.png"), quality=1)
	a = extract_patches(im, 423, 1003, 842, 168)
	save.image(a[[1]], file=paste0(IDs$facetsDir[i], "/facets_custom_6.png"), quality=1)
}


cat_img = function(patID, ii, pic, IDs){
	var_ii = paste0("var", ii)
	par(mfrow=c(length(ii),1)) ## for test
	for(j in 1:length(ii)){
		i = ii[j]
		print(i)
		tmp = load.image(paste0(IDs$facetsDir[i], "/facets_custom_", pic, ".png")); 
		#tmp = implot(tmp, text(50,14,IDs$comID[i], adj=c(0,0)))
		plot(tmp) ## for test
		assign(value = tmp, var_ii[j])
	}
	mget(var_ii) -> tmp
	save.image(imappend(tmp, 'y'), file=paste0(patID, '_', pic, '.png'))
	system(paste0("open ", paste0(patID, '_', pic, '.png')))
}

## SuB2
patID = 'JuB2'
ii = c(1,24,2,23)
for(pic in c(1,2,5,6)){
	cat_img(patID = patID, ii, pic, IDs)
}

## JuB3
patID = 'JuB3'
ii = c(3, 14, 5, 4, 15, 13)
for(pic in c(1,2,5,6)){
	cat_img(patID = patID, ii, pic, IDs)
}

## MaB33
patID = 'MaB33'
ii = c(12, 21, 22)
for(pic in c(1,2,5,6)){
	cat_img(patID = patID, ii, pic, IDs)
}

## MaB30
patID = 'MaB30'
ii = c(6, 19, 10, 8, 18, 20)
for(pic in c(1,2,5,6)){
	cat_img(patID = patID, ii, pic, IDs)
}

## MaB30 X
patID = 'MaB30'
ii = c(6, 19, 10, 8, 18, 20)
for(pic in c(1,2,5,6)){
	cat_img(patID = patID, ii, pic, IDs)
}

maf2[Variant_Type == 'SNP', ] -> tmp
tmp
tmp[, xx := paste(Reference_Allele, Tumor_Seq_Allele2, sep=":")]
table(tmp$xx) ->  xxx
xxx
xxx/sum(xxx)


getwd()
## reRun facets
## Pipelines/FACETS/FACETSv2dev/getFacetCountsTN.sh
library(facets)
library(Cairo)

source(system.file("extRfns", "readSnpMatrixMSK.R", package='facets'))

filenames = c(
"Proj_06048_P_indelRealigned_recal_s_C_000269_T001_d_pileup", 
"Proj_06048_P_indelRealigned_recal_s_C_000269_X001_d_pileup", 
"Proj_06048_P_indelRealigned_recal_s_C_000271_T001_d_pileup", 
"Proj_06048_P_indelRealigned_recal_s_C_000271_X001_d_pileup", 
"Proj_06048_P_indelRealigned_recal_s_C_000271_X002_d_pileup", 
"Proj_06048_P_indelRealigned_recal_s_C_000273_T001_d_pileup", 
"Proj_06048_P_indelRealigned_recal_s_C_000273_T002_d_pileup", 
"Proj_06048_P_indelRealigned_recal_s_C_000273_X001_d_pileup", 
"Proj_06048_P_indelRealigned_recal_s_C_000273_X002_d_pileup", 
"Proj_06048_P_indelRealigned_recal_s_C_000273_X003_d_pileup", 
"Proj_06048_P_indelRealigned_recal_s_C_000273_X004_d_pileup", 
"Proj_06048_P_indelRealigned_recal_s_C_000274_T001_d_pileup", 
"Proj_06048_P_indelRealigned_recal_s_JuB3_P10_Pellet_6048D_pileup", 
"Proj_06048_P_indelRealigned_recal_s_JuB3_P1_Pellet_6048D_pileup", 
"Proj_06048_P_indelRealigned_recal_s_JuB3_P6_Pellet_6048D_pileup", 
"Proj_06048_P_indelRealigned_recal_s_MaB30_2_Org_P1_06048I_pileup", 
"Proj_06048_P_indelRealigned_recal_s_MaB30_2_Org_P8_06048I_pileup", 
"Proj_06048_P_indelRealigned_recal_s_MaB30_Org_P4_06048L_pileup", 
"Proj_06048_P_indelRealigned_recal_s_MaB30_Organoid_6048B_pileup", 
"Proj_06048_P_indelRealigned_recal_s_MaB30_P8_Organoid_6048F_pileup", 
"Proj_06048_P_indelRealigned_recal_s_MaB33_Organoid_6048B_pileup", 
"Proj_06048_P_indelRealigned_recal_s_MaB33_P10_Organoid_6048F_pileup", 
"Proj_06048_P_indelRealigned_recal_s_SuB2_Org_P9_06048I_pileup", 
"Proj_06048_P_indelRealigned_recal_s_SuB2_Pellet_6048D_pileup"
)

cval=300
cat("", file=paste0("log_", cval))
for( i in 1:length(filenames)){
	filename = filenames[i]
	id = sub(".*recal_", "", filename) 
	id = sub("_pileup", "", id) 
	id
	outbase = as.character(IDs[IDs$seqID == id, "comID"])
	outbase
	cat("filename: ", filename, "\nid: ", id, "\noutbase: ", outbase, "\n") 
	cat("filename: ", filename, "\nid: ", id, "\noutbase: ", outbase, "\n", file=paste0("log_", cval), append = T)
	readSnpMatrixMSK(file=paste0("./reRunFacets/", filename, ".gz"), skip=1, type='ns2') -> tmp
	xx = preProcSample(tmp)
	oo = procSample(xx, cval=cval)
	ff = emcncf(oo)
	cat("ff$ploidy: ", ff$ploidy, "\n")
	cat("ff$ploidy: ", ff$ploidy, "\n", file=paste0("log_", cval), append = T)
	save(oo, file=paste0("oo_", outbase, ".RData"))
	save(ff, file=paste0("ff_", outbase, ".RData"))
	CairoPNG(file=paste0("../results/", outbase, "_", cval, ".png"), width=900, height=570)
	plotSample(oo, ff, sname=paste0(outbase, " ploidy: ", round(ff$ploidy,2)))
	dev.off()
}

cval=300
for( i in 1:length(filenames)){
	filename = filenames[i]
	id = sub(".*recal_", "", filename) 
	id = sub("_pileup", "", id) 
	outbase = as.character(IDs[IDs$seqID == id, "comID"])
	load(file=paste0("oo_", outbase, ".RData"))
	load(file=paste0("ff_", outbase, ".RData"))
	#png(file=paste0("../results/", outbase, "_", cval, "_small.png"), width=450, height=235)
	png(file=paste0("../results/", outbase, "_", cval, "_small.png"), width=500, height=370, res=150)
	plotSample(oo, ff, sname=paste0(outbase, " ploidy: ", round(ff$ploidy,2)))
	dev.off()
}

cval=300
for( i in 1:length(filenames)){
	filename = filenames[i]
	id = sub(".*recal_", "", filename) 
	id = sub("_pileup", "", id) 
	outbase = as.character(IDs[IDs$seqID == id, "comID"])
	load(file=paste0("oo_", outbase, ".RData"))
	load(file=paste0("ff_", outbase, ".RData"))
	#png(file=paste0("../results/", outbase, "_", cval, "_small.png"), width=450, height=235)
	pdf(file=paste0("../results/", outbase, "_", cval, "_small.pdf"), width=5, height=2.8)
	plotSample(oo, ff, sname=paste0(outbase, " ploidy: ", round(ff$ploidy,2)))
	dev.off()
}

getwd()
IDs$dipLogR=0
for(i in 1:nrow(IDs)){
	outbase = IDs$comID[i]
	load(paste0("ff_", outbase, ".RData"))
	IDs[IDs$comID == IDs$comID[i], 'dipLogR'] = ff$dipLogR
}
IDs$dipLogR

IDs$parDipLogR = IDs$dipLogR
IDs[, parDipLogR := .SD[grep("SamPar", .SD$comID), dipLogR], by = patientID]
IDs$parDipLogR

## run again by using new dipLogR
cval=200
cat("", file=paste0("log_", cval))
for( i in 1:length(filenames)){
	
	i=14
	filename = filenames[i]
	id = sub(".*recal_", "", filename) 
	id = sub("_pileup", "", id) 
	outbase = as.character(IDs[IDs$seqID == id, "comID"])
	pardipLogR = as.numeric(IDs[IDs$seqID == id, "parDipLogR"])
	cat("pardipLogR: ", pardipLogR, "\n")
	cat("filename: ", filename, "\nid: ", id, "\noutbase: ", outbase, "\n") 
	cat("filename: ", filename, "\nid: ", id, "\noutbase: ", outbase, "\n", file=paste0("log_", cval), append = T)
	readSnpMatrixMSK(file=paste0("./reRunFacets/", filename, ".gz"), skip=1, type='ns2') -> tmp
	xx = preProcSample(tmp)
	oo = procSample(xx, cval=cval, dipLogR = pardipLogR)
	ff = emcncf2(oo, trace=T, min.nhet=30, maxiter=10, difcf=0.05, maxk=5, eps=1e-3)
	#ff = emcncf(oo, trace = T, min.het = 50)
	cat("ff$ploidy: ", ff$ploidy, "\n")
	cat("ff$ploidy: ", ff$ploidy, "\n", file=paste0("log_", cval), append = T)
	save(oo, file=paste0("oo_dipLogR_", outbase, ".RData"))
	save(ff, file=paste0("ff_dipLogR_", outbase, ".RData"))
	CairoPNG(file=paste0("../results/", outbase, "_", cval, "_dipLogR.png"), width=320, height=200)
	plotSample(oo, ff)
	dev.off()

}

pp =  c('SuB2', 'JuB3', 'MaB30', 'MaB33')

for(i in 1:length(filenames)){
	filename = filenames[i]
	id = sub(".*recal_", "", filename) 
	id = sub("_pileup", "", id) 
	outbase = as.character(IDs[IDs$seqID == id, "comID"])
	cat("dealing with ", outbase, "\n")
	load(file=paste0("oo_dipLogR_", outbase, ".RData"))
	load(file=paste0("ff_dipLogR_", outbase, ".RData"))
	pdf(file=paste0("../results/", outbase, "_", cval, ".pdf"), width=10, height=8)
	plotSample(oo, ff, sname=paste0(outbase, " ploidy: ", ff$ploidy))
	dev.off()
}
IDs[,.(comID, seqID)]

oncomatrix2 = aggregate(


oncomatrix2[[2]]
lapply(oncomatrix2, unlist) -> tmpp
head(tmpp)

exon_maf_data$Tumor_Sample_Barcode
oncomatrix = data.table::dcast(exon_maf_data[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode, patientID)], 
			       formula = Hugo_Symbol ~ patientID + Tumor_Sample_Barcode, 
			       fun.aggregate = function(x){
				       x = unique(as.character(x))
				       xad = x[x %in% c('Amp', 'Del')]
				       xvc = x[!x %in% c('Amp', 'Del')]
				       if(length(xvc)>0){
					       xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
				       }
				       x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
				       x = gsub(pattern = ';$', replacement = '', x = x)
				       x = gsub(pattern = '^;', replacement = '', x = x)
				       return(x)
			       } , value.var = 'Variant_Classification', fill = '')

x.n = names(oncomatrix)
x.n = sub("_.*", "", x.n)

Hugo_Symbol = as.character(oncomatrix[1, 1])
jub3 = paste(as.character(oncomatrix[1, 2:7]), collapse = ";"); jub3
mab30 = paste(as.character(oncomatrix[1, 8:13]), collapse = ";"); mab30
mab33 = paste(as.character(oncomatrix[1, 14:16]), collapse = ";"); mab33
mab302 = paste(as.character(oncomatrix[1, 17:21]), collapse = ";"); mab302
sub2 = paste(as.character(oncomatrix[1, 22:25]), collapse = ";"); sub2
oncomatrix3 = data.frame(Hugo_Symbol = Hugo_Symbol, mab30 = mab30, mab302 = mab302, mab33 = mab33, sub2 = sub2, jub3 =jub3)

for(i in 2:nrow(oncomatrix)){
	Hugo_Symbol = as.character(oncomatrix[i, 1])
	jub3 = paste(as.character(oncomatrix[i, 2:7]), collapse = ";"); jub3
	mab30 = paste(as.character(oncomatrix[i, 8:13]), collapse = ";"); mab30
	mab33 = paste(as.character(oncomatrix[i, 14:16]), collapse = ";"); mab33
	mab302 = paste(as.character(oncomatrix[i, 17:21]), collapse = ";"); mab302
	sub2 = paste(as.character(oncomatrix[i, 22:25]), collapse = ";"); sub2
	oncomatrix3 = rbind(oncomatrix3, c(Hugo_Symbol, mab30, mab302, mab33, sub2, jub3))
}
oncomatrix[i, 14:16]
oncomatrix$Hugo_Symbol

oncomatrix3[1:8,1:6]                   # 


oncomatrix2

gg = ggplot(exon_maf_data, aes(x=patientID, y = Variant_Classification, fill = Tumor_Sample_Barcode)) + 
	geom_bar(stat="count")
ggsave(gg, file="../results/test.pdf")

impact_maf_data$Tumor_Sample_Barcode

suk.id = c( "s_MaB19_P11_Organoid_6048F", "s_MaB28_Org_P11_06048I", "s_MaB30_P8_Organoid_6048F", "s_MaB30_2_Org_P14_06048T",
	   "s_MaB33_Org_P13_06048T", "s_JuB3_P10_Pellet_6048D", "s_SuB2_Org_P9_06048I", "s_SuB17_Org_P0_06048I", "s_SuB18_Org_P7_06048T",
	   "s_SuB19_Org_P5_06048T", "s_SuB22_Org_P5_06048T", "s_SuB25_Org_P5_06048T", "s_SuB27_Org_P0_06048L", "s_SuB33_Org_P2_06048T",
	   "s_SuB36_Org_P0_06048T", "s_UTUC17T_PDXP6_Org_P0_06048T")
tmp = sub("s_", "", suk.id)
tmp = sub("_.*", "", tmp)
table(tmp)
suk.sel = impact_maf_data[seqID %in% suk.id, ]
unique(suk.sel$seqID)
unique(impact_maf_data[, .(seqID, patientID, sampleID)]) -> tmp
tmp=tmp[order(tmp[,2]),]
write.table(tmp, file="x", sep="\t", quote=F)
## edit x for batch for oncoplot
tmp = read.table("x", header=T)

gene.sel = c( "FGFR3", "KDM6A", "KMT2D", "KMT2C", "TP53", 
	     "RBM10", "TSC1", "STAG2", "ARID1A", "CTNNB1", 
	     "TP63", "CUL3", "PIK3CA", "DNMT3A", "EPHA3", 
	     "FAT1", "HIST1H3C", "TERT", "MYD88", "RPTOR", "TGFBR2")
gene.sel2 = c('TERT', 'KMT2D', 'FGFR3', 'KDM6A', 'KMT2C', 'TP53', 'CDKN1A', 'RBM10', 'STAG2', 'TSC1', 'KMT2A', 'TP63', 'ARID1A', 'PIK3CA', 'CUL3', 'HIST1H3C', 'TGFBR2', 'EPHA3', 'DAXX', 'RIT1', 'RPTOR')
intersect(gene.sel, gene.sel2)
setdiff(gene.sel, gene.sel2)
unique(c(gene.sel, gene.sel2)) -> gene.sel3
write.table(gene.sel3, file="yy", quote=F, sep="\t") # 

cc = c("Missense_Mutation" = "green", "Frame_Shift_Del"='blue', "Multi_Hit" = "#9d0298", "Nonsense_Mutation"='black', "Splice_Site" = '#914b25', "Frame_Shift_Ins" = 'blue')

for(i in 1:max(tmp$batch)){
	sel = tmp[tmp$batch == i, 'seqID']
	impact_maf_data[seqID %in% sel, ] -> b1
	pdf(width=10, height=15, file=paste0("../results/impact_oncoplot_b_genelis3_", i, ".pdf"))
	oncoplot(maf = impact_maf, genes = gene.sel3)
	dev.off()
}

write.table(impact_id[,.(comID, seqID)], file="impact_id2")

maffile = "/ifs/res/share/solit/solitd/Proj_06048_P/r_001/post/Proj_06048_P___SOMATIC.vep.filtered.facets.V3.maf"
#maffile = "../../Proj_06048_P___SOMATIC.vep.filtered.facets.V3.maf"
tmp = data.table::fread(maffile)
tmp[, seqID := Tumor_Sample_Barcode]
setkey(tmp, 'Tumor_Sample_Barcode')
tmp = tmp[IDs[,.(seqID, comID, patientID, patientID2, sampleID, sampleID2, newid)]]
tmp[, FILTER2 := FILTER]
tmp[FILTER2 == "PASS", FILTER := "."]
tmp[FILTER2 == "RESCUE", FILTER := "."]
tmp[Chromosome == "MT", Chromosome := 'M']
tmp[, vaf := 100 * t_alt_count / t_depth]
tmp[,TAG2 := paste0(Hugo_Symbol, ":", TAG)]
tmp[, Tumor_Sample_Barcode := newid]
lll = c('Tumor_Sample_Barcode', 'comID', 'Hugo_Symbol', 'ccf_Mcopies', 'ccf_Mcopies_prob95', 'HGVSp_Short', 'FILTER')
tmp2 = tmp[patientID == 'JuB3', lll, with=F]
gg = c('ERBB3', 'JAK2', 'ERBB2', 'TP53', 'RB1', 'TERT', 'STAG2', 'KMT2D', 'NF1', 'ERCC4' )
tmp3 = tmp2[Hugo_Symbol %in% gg, ]
nrow(tmp3)
tmp3 = as.data.frame(tmp3)
tmp3
write.table(tmp3, file="results/a.txt")
maf[patientID == 'JuB3' & Hugo_Symbol == 'ERBB3', c(lll, 'Variant_Classification', 'Chromosome', 'Start_Position', 'RefSeq'), with=F]
maf[patientID == 'JuB3' & Hugo_Symbol == 'ERBB3', c(lll, 'Variant_Classification', 'Chromosome', 'Start_Position', 'RefSeq', 'Tumor_Seq_Allele2', 'Reference_Allele'), with=F]

colnames(tmp)
lll = c('Tumor_Sample_Barcode', 'comID', 'Hugo_Symbol', 'newid', 'vaf')
pid = c('JuB3', 'MaB33', 'SuB2', 'MaB30', 'MaB30-2', 'MaBX30')
gs = c('TSC1', 'TERT', 'FGFR3', 'KDM6A', 'STAG2', 'CTNNB1', 'TP53', 'ARID1A', 'ERBB3', 'JAK2', 'ERBB2', 'RB1', 'KMT2D', 'NF1', 'ERCC4')
tmp = impact_maf_data[patientID %in% pid & Hugo_Symbol %in% gs , ]
unique(tmp$patientID)
unique(tmp$sampleID)
tmp[, sampleID := factor(sampleID, levels=c('Par', 'Par_slides', 'Par_Pellet', 'OrgP1', 'OrgP4', 'OrgP6', 'OrgP8', 'OrgP10', 'OrgP13', 'OrgP14', 'PDX', 'PDXOrgP2'))]
tmp[, TAG2 := paste0(Hugo_Symbol, ":", HGVSp_Short)]

tert = fread('tert1.tsv')
setkey(tert, "Sample ID")
setkey(impact_id, "seqID")
ttt = tert[impact_id[, .(Tumor_Sample_Barcode, seqID,patientID, sampleID)]]
tttt = ttt[, .(`Sample ID`, Var, sampleID, patientID, `Allele Freq (T)`)]
tttt$Hugo_Symbol = 'TERT'
tttt$TAG = 'Promoter'
tttt$TAG2 = 'TERT:Promoter'
setnames(tttt, "Allele Freq (T)", "vaf")
tttt[, vaf := 100 * vaf]
tttt
tail(tmp3)

tmp2 = rbind, 'OrgP9'(tmp[,.(sampleID, patientID, Hugo_Symbol, TAG2, vaf)], tttt[, .(sampleID, patientID, Hugo_Symbol, TAG2, vaf)])
tmp2[is.na(sampleID), 1] = 'PDXOrgP0'
tmp2[is.na(sampleID), 1]

pat = 'MaB30'
tmp3 = tmp2[patientID == pat,]
tmp3
tmp3[, sampleID := factor(sampleID, levels = levels(droplevels(sampleID)))]
gg = ggplot(tmp3, aes(x = sampleID, y = vaf, color =TAG2, group =TAG2)) + 
	geom_point() + geom_line() 
ggsave(paste0("../results/impact_vaf_phylogeny_pat", pat, ".pdf"), width=7, height=3)

pat = 'SuB2'
tmp3 = tmp2[patientID == pat,]
tmp3
pa = data.table(t(c('PDXOrgP0', 'SuB2', 'TERT', 'TERT:Promoter', 77.3)))
setnames(pa, c('sampleID', 'patientID', 'Hugo_Symbol', 'TAG2', 'vaf'))
tmp4 = rbindlist(list(tmp3, pa))
tmp4[, sampleID := factor(sampleID, levels =c('Par_slides', 'Par_Pellet', 'OrgP9', 'PDX', 'PDXOrgP0'))]
#tmp4 = as.data.table(tmp4)
#tmp4
#tmp3[, sampleID := factor(sampleID, levels = levels(droplevels(sampleID)))]
gg = ggplot(tmp4, aes(x = sampleID, y = vaf, color =TAG2, group =TAG2)) + 
	geom_point() + geom_line() 
ggsave(paste0("../results/impact_vaf_phylogeny_pat", pat, ".pdf"), width=7, height=3)

pats = c('SuB2', 'JuB3', 'MaB30', 'MaBX30', 'MaB33')
for(pat in pats){
tmp3 = tmp2[patientID == pat,]
tmp3[, sampleID := factor(sampleID, levels = levels(droplevels(sampleID)))]
gg = ggplot(tmp3, aes(x = sampleID, y = vaf, color =TAG2, group =TAG2)) + 
	geom_point() + geom_line() 
ll = length(unique(tmp3$TAG2))
ggsave(paste0("../results/impact_vaf_phylogeny2_pat", pat, ".pdf"), width=7, height=3)
}

unique(maf$patientID)

impact_maf_data[patientID == 'SuB2', Hugo_Symbol]
impact_maf_data[patientID == 'MaBX30', Hugo_Symbol]

write.table(tmp[, .(sampleID, vaf, Hugo_Symbol, patientID, HGVSp_Short)], file="../results/tmp.txt")
sum(c(509, 563, 533, 1510, 1716, 1961, 1591, 1097, 1090, 146, 103, 128, 126, 95, 106, 128, 136, 123, 114, 84, 503, 503, 532))
table(maf$Variant_Classification)

## uromol
head(uromol)
rsem_shen$fpkm = rsem_shen$counts 
rsem.fpkm(rsem_shen) -> rsem_shen
dim(rsem_shen$fpkm)
dim(uromol)
uromo
read.table("UROMOL_gene_fpkm_gtf.txt", row.names=1, header = T) -> uromol
row.names(uromol) = sub("\\..*$", "", row.names(uromol))
uromol[1:10,1:10]
dim(uromol)
rsem_shen$fpkm[1:10,1:10]
ids = data.frame(row.names = sub("_.*", "", row.names(rsem_shen$fpkm)), ensg = row.names(rsem_shen$fpkm), stringsAsFactors = F) 
ov = intersect(row.names(uromol), row.names(ids))
uromol = uromol[ov, ]
row.names(uromol) = ids[ov, 'ensg']
uromol[1:10,1:10]

cbind(uromol, rsem_shen$fpkm[row.names(uromol),]) -> uro
uro[, 1:3] = NULL
uro[1:10,1:10]

## uromol pca
uro.sd = apply(uro, 1, sd)

tmp = uro[order(uro.sd, decreasing = T), ]
tmp = tmp[1:1000,]
tmp = log2(tmp+1)
head(tmp)
tmp = prcomp(t(tmp))
tmp = tmp$x
tmp = as.data.frame(tmp)
tmp[1:10,1:10]
tmp.bak = tmp

tmp = tmp.bak
tmp$samp_class= "Orgnoid"
tmp$samp_class[grep("^U", row.names(tmp))] = "UROMOL"
tmp$samp_class[grep("umor", row.names(tmp))] = "tumor"
tmp$col = 2
tmp$col[grep("Org", row.names(tmp))] = 3
tmp$col[grep("UROMOL", row.names(tmp))] = 4
tmp$col[grep("umor", row.names(tmp))] = 5
tail(row.names(tmp))

range(tmp$PC1)
range(tmp$PC2)
tmp = tmp[tmp$PC1 > 0,] 
tmp = tmp[tmp$PC2 <  2500,] 
dim(tmp)
tmp[1:10,1:10]
row.names(tmp)

tmp$samp_class
pdf("results/uromol_pca.pdf", width=8, height=6)
ggplot(tmp, aes(PC1, PC2, label = samp_class,  color=samp_class)) + geom_point() + theme(legend.position="right") 
dev.off()
system('rsync -avur  results/ mski1925:/Volumes/LaCie/huw/solit/study/hiseq/shen/results/')

rm(uro.sd, tmp)


## immune cell infiltration
## ssGSEA 
library(GSVA)
immuno_sig =fread("~/program/immuno_sig.txt")
immuno_sig
gg = unique(immuno_sig$CellType)
geneset = list(immuno_sig[CellType == gg[1], Symbol])
for(i in 2:length(gg)){
	geneset = c(geneset, list(immuno_sig[CellType == gg[i], Symbol]))
}
names(geneset) = gg
length(geneset)
names(geneset)
geneset

adaptive = c( "Bcells", "Tcells", "Thelpercells", "Tcm", "Th1cells", "Th2cells", "TFH", "Th17cells", "TReg", "CD8Tcells", "Tgd", "Cytotoxiccells")
innate = c( "NKcells", "NKCD45dimencells", "NKCD56brightcells", "DC", "iDC", "aDC", "pDC", "Eosinophils", "Macrophages", "Mastcells", "Neutrophils") 
anno_row = data.frame(row.names = names(geneset), ImmuneType = rep("Adaptive", 28), stringsAsFactors = F)
anno_row$ImmuneType[row.names(anno_row) %in% innate] = "Innate"
names(geneset)
head(geneset)

tmp = blca.cc.log2
row.names(tmp) = blca.rows[row.names(tmp), 2]
gsva.score = gsva(tmp, geneset, method="ssgsea")
head(gsva.score)
dim(gsva.score)
blca.clin

## significance test of the GSVA score by using limma
library(limma)
library(xlsx)
bcg.ps = read.xlsx("../bcg/Updated Database for TME in PvS.xlsx", 1, stringsAsFactors=F)
bcg.ps = as.data.table(bcg.ps)
colnames(bcg.ps)
head(bcg.ps)

gsva.anno = data.table(tsb =colnames(gsva.score))
gsva.anno[, sid := substr(tsb, 1, 15)]
gsva.anno = merge(gsva.anno, bcg.ps[, .(Sample.ID, Primary.MIBC.V.Secondary)], by.x='sid', by.y='Sample.ID', all.x=T)
gsva.anno

cn = colnames(gsva.score)
cn
condition = rep(1, length(cn))
condition[grep('primary', cn)] = 0
condition

design <- model.matrix(~ condition)
design
colnames(design) <- c("ALL", "PDXvsPrimary")
fit <- lmFit(gsva.score, design)
fit <- eBayes(fit)
gsva.res <- topTable(fit, coef="PDXvsPrimary", number=Inf)
gsva.res
summary(decideTests(fit, p.value=0.05))
summary(gsva.res)
gsva.res[, Significance := 'Not']
gsva.res[P.Value < 0.05, Significance := 'Sig']
gsva.res = as.data.table(gsva.res, keep.rownames=T)
gsva.res

save(gsva.score, gsva.res, file='gsva_score.RData')

gg = ggplot(gsva.res, aes(x=reorder(rn, P.Value), y=logFC, color = Significance, size=-log(adj.P.Val))) + geom_point() +
	theme(axis.text.x = element_text(angle=-90, vjust = 0.5, hjust=0)) + ggtitle('Immune cell infilteration \nby ssGSEA (PDX vs Primary) ') + 
	ylab("logFC") + xlab("Immune cells") + coord_flip() + 
	geom_hline(yintercept=0, color=adjustcolor(2, .3))
ggsave(gg, file='res/immune_cell_deconv_ssGSEA.pdf', width=5, height=5)
sync('res')

cn = colnames(gsva_score)
matrix_tmp = t(apply(gsva_score, 1, scale))
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
colnames(matrix_tmp) = cn
breaks = c( seq(from = -2, to = 2, length.out = 30))
pdf("res/heatmap_gsva.pdf")
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 9, show_rownames = T,
        annotation_row = anno_row)
dev.off()
sync('res')

gsva_hp = hp
cutree(hp$tree_col, 2) -> gsva_cluster

source("~/program/fun/topn_pheatmap.R")
source("~/program/fun/my_pheatmap.R")
source("~/program/fun/my_scale.R")
source("~/program/fun/top_rows.R")
source("~/program/fun/remove_null_row.R")
source("~/program/fun/consensus_cluster.R")

## heatmap of the immuno sig genes from ssGSEA
gs = sub(".*_", "", row.names(res.cc.log))
res.cc.log[gs %in% immuno_sig$Symbol, ] -> tmp
#my_scale(tmp) -> tmp
cn = colnames(tmp)
matrix_tmp = t(apply(tmp, 1, scale))
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
colnames(matrix_tmp) = cn
breaks = c( seq(from = -2, to = 2, length.out = 30))
remove_null_row(matrix_tmp) -> matrix_tmp
pdf("res/heatmap_immuno_sig.pdf")
hp = pheatmap(tmp, color=greenred(length(breaks)+1), scale='none', clustering_method = 'ward.D2', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 9, show_rownames = F)
dev.off()

## cibersort
tmp = as.data.frame(res.cc)
head(tmp)
gs = sub(".*_", "", row.names(tmp))
gs1 = gs[!duplicated(gs)]
tmp = tmp[!duplicated(gs),]
row.names(tmp) = gs1
head(tmp)
dim(tmp)
rm(gs1, gs)
write.table(tmp, file='tmp.tsv', sep="\t", quote=F)
## add one more column names

##library(Rserve)
##Rserve(args="--no-save")
##cmd = paste0(cfg$java, ' -Xmx5g -Xms3g -jar ', cfg$cibersort, '/CIBERSORT.jar ') 
##cmd = paste0(cmd, ' -M tmp.tsv -B ', cfg$cibersort, '/LM22.txt > cibersort.tsv')
##cmd
##exe.jobs(cmd, logger)

source('~/program/cibersort/CIBERSORT.R')
ciber = CIBERSORT(paste0(cfg$cibersort, '/LM22.txt'), 'tmp.tsv', perm=1000, QN=TRUE, absolute=T)
ciber.t = t(ciber)
ciber.t
ciber.t = ciber.t[1:22, ]
tail(ciber.t)
dim(ciber.t)
ciber.t[-1,]

library(limma)
cn = colnames(ciber.t)
condition = rep(1, length(cn))
condition[grep('primary', cn)] = 0
condition

design <- model.matrix(~ condition)
design
colnames(design) <- c("ALL", "PDXvsPrimary")
fit <- lmFit(ciber.t, design)
fit <- eBayes(fit)
summary(decideTests(fit, p.value=0.05))
ciber.res <- topTable(fit, coef="PDXvsPrimary", number=Inf)
ciber.res = as.data.table(ciber.res, keep.rownames=T)
ciber.res[, Significance := 'Not']
ciber.res[P.Value < 0.05, Significance := 'Sig']
ciber.res

save(ciber.res, ciber.t, file='cibersort_score_39_samples.RData')

gg = ggplot(ciber.res, aes(x=reorder(rn, P.Value), y=logFC, color = Significance, size=-log(adj.P.Val))) + geom_point() +
	theme(axis.text.x = element_text(angle=-90, vjust = 0.5, hjust=0)) + ggtitle('Immune cell infilteration \nby cibersort (PDX vs Primary) ') + 
	ylab("logFC") + xlab("Immune cells") + coord_flip() + 
	geom_hline(yintercept=0, color=adjustcolor(2, .3))
ggsave(gg, file='res/immune_cell_deconv_39_samples_cibersort.pdf', width=4, height=5)
sync('res')

cn = colnames(ciber.t)
matrix_tmp = t(apply(ciber.t, 1, scale))
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
colnames(matrix_tmp) = cn
breaks = c( seq(from = -2, to = 2, length.out = 30))
pdf("res/heatmap_cibersort_score.pdf")
hp = pheatmap(matrix_tmp, color=greenred(length(breaks)+1), scale='none', clustering_method = 'complete', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 9, show_rownames = T,
        annotation_row = anno_row)
dev.off()
sync('res')

## heatmap of the immuno sig genes from cibersort
ciber.gs = fread(paste0(cfg$cibersort, '/LM22.txt'))
ciber.gs = ciber.gs$`Gene symbol`
ciber.gs
gs = sub(".*_", "", row.names(res.cc.log))
res.cc.log[gs %in% ciber.gs, ] -> tmp
#my_scale(tmp) -> tmp
cn = colnames(tmp)
matrix_tmp = t(apply(tmp, 1, scale))
matrix_tmp[matrix_tmp > 2] = 2
matrix_tmp[matrix_tmp < -2] = -2
colnames(matrix_tmp) = cn
breaks = c( seq(from = -2, to = 2, length.out = 30))
remove_null_row(matrix_tmp) -> matrix_tmp
pdf("res/heatmap_immuno_sig_cibersort.pdf")
hp = pheatmap(tmp, color=greenred(length(breaks)+1), scale='none', clustering_method = 'ward.D2', 
        cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = T, fontsize_row = 9, show_rownames = F)
dev.off()
