library(openxlsx)
library(car)
setwd("/export/home/lmacias/analisis_asociacion/convivencias_medix_qxob/adultos_mega")
lipid <- loadWorkbook("cirugia_bar/fenotipos/metabolomica/lipidos/AHV_NAFLD_Plasma_Lipidomics.xlsx")
lip_desglos <- read.xlsx(lipid)
muestras_elim <- c("GEA23", "GEA24", "GEA37", "RL315", "RL333", "RL335")
tmp <- loadWorkbook("cirugia_bar/fenotipos/Sujetos_hipolipemiantes.xlsx")
muestras_elim <- c(muestras_elim, read.xlsx(tmp)[,1])
muestras_elim[1:3] <- gsub("GEA", "GEAo0", muestras_elim[1:3])
rm(tmp)
lip_desglos$Sample.Type[nchar(lip_desglos$Sample.Type)==4] <- gsub("RL", "RL0", lip_desglos$Sample.Type[nchar(lip_desglos$Sample.Type)==4])
lip_desglos$Sample.Type[grep("GEA", lip_desglos$Sample.Type)] <- gsub("GEA", "GEAo0", lip_desglos$Sample.Type[grep("GEA", lip_desglos$Sample.Type)])
rownames(lip_desglos) <- lip_desglos$Sample.Type
lip_desglos <- lip_desglos[!rownames(lip_desglos) %in% muestras_elim,]
lip_desglos<- lip_desglos[-(1:2)]

names(lip_desglos) <- gsub("\\(", "_", names(lip_desglos))
names(lip_desglos) <- gsub("\\)", "", names(lip_desglos))
names(lip_desglos) <- gsub("\\:", "", names(lip_desglos))
names(lip_desglos) <- gsub("\\/", "", names(lip_desglos))
names(lip_desglos) <- gsub("\\[", "", names(lip_desglos))
names(lip_desglos) <- gsub("\\]", "", names(lip_desglos))
names(lip_desglos) <- gsub("\\-", "", names(lip_desglos))
names(lip_desglos) <- gsub('\\\\', "", names(lip_desglos))
names(lip_desglos) <- gsub("\\+", "", names(lip_desglos))

blancos <- c("Cer(d18:0/18:0)", "Cer(d18:0/24:1)", "Cer(d17:1/18:0)", "Cer(d18:2/26:0)", "Hex1Cer(d16:1/18:0)", "Hex2Cer(d18:2/24:1)", "Hex3Cer(d18:1/16:0)", "PE(18:0/18:1)", "PE(16:0/16:1)")
blancos <- gsub("\\(", "_", blancos)
blancos <- gsub("[:punct:]|\\/|\\)", "", blancos)
blancos <- gsub("PE_180", "PE_180_", blancos)
blancos <- gsub("PE_160", "PE_160_", blancos)
blancos <- lip_desglos[grep(paste(blancos, collapse="$|^"), names(lip_desglos))]

load("cirugia_bar/fenotipos/histologia_revisada.RData")
rm(biliares, pando, muestras_elim)
all(rownames(blancos) %in% morbidos$ID)
##TRUE
blancos <- cbind(morbidos[pmatch(rownames(blancos), morbidos$ID), c("IID", "GENERO", "EDAD", "IMC", "PC1", "PC2")], blancos)
boxCox(lm(Cer_d171180~GENERO + EDAD + IMC + PC1 + PC2, data=blancos))
#raíz cuarta
boxCox(lm(Cer_d182260~GENERO + EDAD + IMC + PC1 + PC2, data=blancos))
#raíz cuadrada
blancos$Cer_d180180[blancos$Cer_d180180==0] <- 1e-6
boxCox(lm(Cer_d180180~GENERO + EDAD + IMC + PC1 + PC2, data=blancos))
#raíz cuadrada
boxCox(lm(Cer_d180241~GENERO + EDAD + IMC + PC1 + PC2, data=blancos))
#logaritmo
blancos$Hex1Cer_d161180[blancos$Hex1Cer_d161180==0] <- 1e-6
boxCox(lm(Hex1Cer_d161180~GENERO + EDAD + IMC + PC1 + PC2, data=blancos))
#raíz cuadrada
boxCox(lm(Hex2Cer_d182241~GENERO + EDAD + IMC + PC1 + PC2, data=blancos))
#raíz cuadrada
boxCox(lm(Hex3Cer_d181160~GENERO + EDAD + IMC + PC1 + PC2, data=blancos))
#sin cambios
boxCox(lm(PE_160_161~GENERO + EDAD + IMC + PC1 + PC2, data=blancos))
#logaritmo
boxCox(lm(PE_180_181~GENERO + EDAD + IMC + PC1 + PC2, data=blancos))
#logaritmo

blancos[7] <- blancos[7]^0.25
blancos[c(8:9,11:12)] <- sqrt(blancos[c(8:9,11:12)])
blancos[c(10,14:15)] <- log(blancos[c(10,14:15)])

cgeno <- read.table("cirugia_bar/imputados/fusionados/adultos.082019.QC.fam")
cgeno <- cgeno[grepl("^RL|^GEA", cgeno$V2),]
blancos$IID[!(blancos$IID %in% cgeno$V2)]
##revisamos uno por uno, ninguno tiene microarreglo ADN

fenos <- data.frame(FID=0, blancos[blancos$IID %in% cgeno$V2,c(1,7:15)])
write.table(fenos, "cirugia_bar/fenotipos/para_asoc/Cer_Hex2Cer_PE_selectos.phen", quote=FALSE, row.names=FALSE)
