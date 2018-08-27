library(ggplot2)
	library(broom)
	setwd("~/analisis_asociacion/cardio/ancestria")
	ances <- read.table("ancestrales_GEA.3.Q")
	ances <- cbind(read.table("ancestrales_GEA.fam")[c(1:2)], ances)
	ids_ceu <- scan("~/CANDELA/haplotipos/ref_beagle/listas/grupos/ids_CEU.txt", what='list')
	ids_yri <- scan("~/CANDELA/haplotipos/ref_beagle/listas/grupos/ids_YRI.txt", what='list')
	indigenas <- grepl("FN|TOT|ZAP|NFM", ances$V2)
	summary(ances[ances$V2 %in% ids_ceu,3:5])
	##V1 CEU
	summary(ances[ances$V2 %in% ids_yri,3:5])
	##V3 YRI
	summary(ances[indigenas,3:5])
	##V2 NAT
	summary(ances[!indigenas & !grepl("NA", ances$V2),3:5])
	names(ances) <- c("FID", "IID", "CEU", "NAT", "YRI")	
	mestizos <- ances[!indigenas & !grepl("NA", ances$NAT),]
	##tst <- wilcox.test(mestizos$CEU ~ grepl("GC", mestizos$IID))
	##pval = tidy(tst)$p.value
	##plt1 <- ggplot(mestizos, aes(x=grepl("GC", mestizos$IID), y=NAT))+ 
	##geom_boxplot()
	##plt1 + geom_text(aes(x=0.7, y = 0.82, label=paste('P-value:', format.pval(pval, digits=2))))
	boxplot(mestizos[,3] ~ grepl("GC", mestizos$V2), col=c("honeydew3", "indianred"), notch=TRUE, names=c("Casos", "Controles"), ylab="Proporción de componente ancestral europeo", cex.lab=1.15, cex.axis=1.15)
	text(2.3, 0.83, "p=0.0020", cex=1.4)
	savePlot("prop_europeo_cc.png", type="png")
