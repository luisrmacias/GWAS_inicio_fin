library(ggplot2)
library(broom)
carpeta <- "~/analisis_asociacion/adultos_mega/calidad/cluster3"
estudio <- "ancestrales_adultos_indep"
setwd(carpeta)
ances <- read.table(paste0(estudio,".3.Q"))
ances <- cbind(read.table(paste0(estudio, ".fam"))[c(1:2)], ances)
ids_ceu <- scan("~/CANDELA/haplotipos/ref_beagle/listas/grupos/fundadores_CEU.txt", what='list')
ids_yri <- scan("~/CANDELA/haplotipos/ref_beagle/listas/grupos/fundadores_YRI.txt", what='list')
indigenas <- grepl("FN|TOT|ZAP|NFM", ances$V2)
summary(ances[ances$V2 %in% ids_ceu,3:5])
##V3 CEU
columna_CEU <- which(apply(ances[ances$IID %in% ids_ceu,3:5], 2, median)>0.95)
summary(ances[ances$V2 %in% ids_yri,3:5])
##V2 YRI
columna_YRI <- which(apply(ances[ances$IID %in% ids_yri,3:5], 2, median)>0.95)
summary(ances[indigenas,3:5])
##V1 NAT
summary(ances[!indigenas & !grepl("NA", ances$V2),3:5])
columna_NAT <- which(apply(ances[indigenas,3:5], 2, median)>0.95)

names(ances)[1:2] <- c("FID", "IID")
names(ances)[c(columna_NAT, columna_CEU, columna_YRI)+2] <- c("NAT", "CEU", "YRI")
mestizos <- ances[!indigenas & !grepl("NA", ances$IID),]
	tst <- wilcox.test(mestizos$CEU ~ grepl("GC|ICN", mestizos$IID))
	pval = tidy(tst)$p.value
	plt1 <- ggplot(mestizos, aes(x=grepl("GC|ICN", mestizos$IID), y=NAT))+ 
	geom_boxplot()
	plt1 + geom_text(aes(x=0.7, y = 0.82, label=paste('P-value:', format.pval(pval, digits=2))))
boxplot(mestizos$CEU ~ grepl("GC|ICN", mestizos$IID), col=c("honeydew3", "indianred"), notch=TRUE, names=c("Casos", "Controles"), ylab="Proporción de componente ancestral europeo", cex.lab=1.15, cex.axis=1.15)
text(0.7, 0.83, bquote(p== .(signif(pval,4))), cex=1.15)
savePlot("prop_europeo_cc.png", type="png")
