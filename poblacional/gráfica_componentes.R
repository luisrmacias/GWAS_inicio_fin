carpeta <- "ancestrales.adultos.082019.indep.eigenvec"
pca <- read.table(carpeta, header=TRUE)
ids_ceu <- scan("~/analisis_asociacion/colaboraciones/CANDELA/haplotipos/listas/grupos/fundadores_CEU.txt", what='list')
ids_yri <- scan("~/analisis_asociacion/colaboraciones/CANDELA/haplotipos/listas/grupos/fundadores_YRI.txt", what='list')
por_centro <- read.table("../grupo_asig.txt")
europeos <- na.omit(pmatch(ids_ceu, pca$IID))
yoruba <- na.omit(pmatch(ids_yri, pca$IID))
indigenas <- grep("FN|TOT|ZAP|NFM", pca$IID)
ref <- c(indigenas, europeos, yoruba)
status <- read.table("../../../cirugia_bar/imputados/fusionados/adultos.032019.QC.fam")
casos <- which(pca$IID %in% status[status[,6]==2,2])
controles <- which(pca$IID %in% status[status[,6]==1,2])
desconocidos <- which(pca$IID %in% status[status[,6]==-9,2])
medix <- grep("MZ1|CON|SOL|RED|TER", pca$IID)
qx_ob <- grep("RL|GEA", pca$IID)
controles <- grep("GC|ICN", pca$IID)        
pob <- character(dim(pca)[1])
pob[medix] <- "Medix"
pob[qx_ob] <- "Cirugía bariátrica"
pob[casos] <- "Con obesidad"
pob[controles] <- "Normopeso"
pob[desconocidos] <- "Otro estado nutricio"
pob[europeos] <- "CEU"
pob[indigenas] <- "Indígenas"
pob[yoruba] <- "YRI"
table(pob)
pob[pob==""] <- NA
pob <- factor(pob)
pca <- cbind(pob, pca)
rm(pob)
pca <- pca[pca$pob!="",]
png('mds_porestudio.png')
	plot(pca$PC1, pca$PC2, col=c("turquoise", "orange","green","black","purple","brown")[pca$pob], pch=18, xlab="First principal component", ylab="Second principal component")
        ##title(main="Distribución de componentes principales de los individuos\ndel estudio con respecto a poblaciones ancestrales")
        legend(0, 0.07, levels(pca$pob), bg="white", col=c("turquoise", "orange","green","black","purple","brown"), pch=18, pt.cex=2, y.intersp=1.15)
	dev.off()
