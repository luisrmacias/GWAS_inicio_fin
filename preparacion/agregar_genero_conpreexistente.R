archivos <- list.files(path="calidad/5_fases", pattern="*.fam", full.names=TRUE)
ids <- lapply(archivos, read.table)
ids <- rbind(ids[[1]], ids[[2]], ids[[3]], ids[[4]])
sexos <- read.table("cirugia_bar/imputados/fusionados/fusion_agosto2017.fam")
sexos <- rbind(sexos, read.table("calidad/cluster2/terceros_adultos.fam"))
sexos <- cbind(ids[,1:2], sexos$V5[pmatch(ids$V2, sexos$V2)])
names(sexos) <- c("FID", "IID", "SEX")
##sexos$SEX[is.na(sexos$SEX)] <- as.numeric(relevel(bioquimicos_GEA$Sexo[pmatch(sexos$IID[is.na(sexos$SEX)], bioquimicos_GEA$IID)], ref="Hombre"))
##sexos$SEX[is.na(sexos$SEX)] <- as.numeric(relevel(bioquimicos_GEA$Sexo[pmatch(sexos$IID[is.na(sexos$SEX)], bioquimicos_GEA$id_MA)], ref="Hombre"))
##faltan <- grep("GC", sexos$IID[is.na(sexos$SEX)][!sexos$IID[is.na(sexos$SEX)] %in% bioquimicos_GEA$id_MA], value=TRUE)
##sexos$SEX[grep(paste(faltan, collapse="|"), sexos$IID)] <-  as.numeric(relevel(bioquimicos_GEA$Sexo[grep(paste(gsub("GC0", "GC_0", faltan), collapse="|"), bioquimicos_GEA$IID)], ref="Hombre"))

sexos$SEX[is.na(sexos$SEX)] <- 0
write.table(sexos, "listas/integrar_sexos_cero.txt", row.names=FALSE, quote=FALSE)
