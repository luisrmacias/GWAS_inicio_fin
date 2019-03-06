##idenficación de variantes duplicadas en archivo .bim
args <- commandArgs(TRUE)
snp <- read.table(args)
chr_pos <- snp[,c(1,4)]
dobles <- snp[duplicated(chr_pos),2]
dobles <- dobles[which(dobles!=levels(dobles)[1])]
dobles <- droplevels(dobles)
write.table(dobles, "~/analisis_asociacion/NMO/listas/eliminacion_QC/SNP_duplicados.list", row.names=FALSE, col.names=FALSE, quote=FALSE)

