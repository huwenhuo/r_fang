library(data.table)
library(GenomicRanges)
options(width = 400)

setwd('/ifs/depot/resources/dmp/data_20150226/mskdata/interval-lists/VERSIONS')
list.files('cv5')

scan('cv3/genelist', character()) # 341
scan('cv5/genelist', character()) # 410
scan('cv6/genelist', character()) # 468

bd = fread('cv3/picard_baits.interval_list', skip='tiled', header=F)
setnames(bd, c('chr', 'start', 'end', 'strand', 'name'))
bd[, ww := end - start + 1]
sum(bd$ww) # 1656966

bd = fread('cv5/picard_baits.interval_list', skip='tiled', header=F)
setnames(bd, c('chr', 'start', 'end', 'strand', 'name'))
bd[, ww := end - start + 1]
sum(bd$ww) # 1847924

bd = fread('cv6/picard_baits.interval_list', skip='target', header=F)
setnames(bd, c('chr', 'start', 'end', 'strand', 'name'))
bd[, ww := end - start + 1]
sum(bd$ww) # 2077695
