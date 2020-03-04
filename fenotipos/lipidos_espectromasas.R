library(openxlsx)
library(car)
setwd("/export/home/lmacias/analisis_asociacion/convivencias_medix_qxob/adultos_mega")
lipid <- loadWorkbook("cirugia_bar/fenotipos/metabolomica/lipidos/AHV_NAFLD_Plasma_Lipidomics.xlsx")
lip_sumas <- read.xlsx(lipid, sheet=2)
muestras_elim <- c("GEA23", "GEA24", "GEA37", "RL315", "RL333", "RL335")
tmp <- loadWorkbook("cirugia_bar/fenotipos/Sujetos_hipolipemiantes.xlsx")
muestras_elim <- c(muestras_elim, read.xlsx(tmp)[,1])
muestras_elim[1:3] <- gsub("GEA", "GEAo0", muestras_elim[1:3])
rm(tmp)
lip_sumas$Sample.Type[nchar(lip_sumas$Sample.Type)==4] <- gsub("RL", "RL0", lip_sumas$Sample.Type[nchar(lip_sumas$Sample.Type)==4])
lip_sumas$Sample.Type[grep("GEA", lip_sumas$Sample.Type)] <- gsub("GEA", "GEAo0", lip_sumas$Sample.Type[grep("GEA", lip_sumas$Sample.Type)])
rownames(lip_sumas) <- lip_sumas$Sample.Type
lip_sumas <- lip_sumas[!rownames(lip_sumas) %in% muestras_elim,]
lip_sumas<- lip_sumas[-(1:2)]

names(lip_sumas) <- gsub("\\(", "_", names(lip_sumas))
names(lip_sumas) <- gsub("\\)", "", names(lip_sumas))
names(lip_sumas) <- gsub("\\:", "", names(lip_sumas))
names(lip_sumas) <- gsub("\\/", "", names(lip_sumas))
names(lip_sumas) <- gsub("\\[", "", names(lip_sumas))
names(lip_sumas) <- gsub("\\]", "", names(lip_sumas))
names(lip_sumas) <- gsub("\\-", "", names(lip_sumas))
names(lip_sumas) <- gsub('\\\\', "", names(lip_sumas))
names(lip_sumas) <- gsub("\\+", "", names(lip_sumas))

blancos <- c("Total.PE", "Total.THC", "Total.dhCer")
blancos <- lip_sumas[grep(paste(blancos, collapse="$|^"), names(lip_sumas))]

load("cirugia_bar/fenotipos/histologia_revisada.RData")
rm(biliares, pando, muestras_elim)
all(rownames(blancos) %in% morbidos$ID)
##TRUE
blancos <- cbind(morbidos[pmatch(rownames(blancos), morbidos$ID), c("IID", "GENERO", "EDAD", "IMC", "PC1", "PC2")], blancos)
boxCox(lm(Total.PE~GENERO + EDAD + IMC + PC1 + PC2, data=blancos))
#raíz cuarta
boxCox(lm(Total.THC~GENERO + EDAD + IMC + PC1 + PC2, data=blancos))
#raíz cuadrada
boxCox(lm(Total.dhCer~GENERO + EDAD + IMC + PC1 + PC2, data=blancos))
#logaritmo natural

blancos[7:9] <- cbind(blancos[7]^0.25, sqrt(blancos[8]), log(blancos[9]))

cgeno <- read.table("cirugia_bar/imputados/fusionados/adultos.082019.QC.fam")
cgeno <- cgeno[grepl("^RL|^GEA", cgeno$V2),]
blancos$IID[!(blancos$IID %in% cgeno$V2)]
##revisamos uno por uno, ninguno tiene microarreglo ADN

fenos <- with(blancos[blancos$IID %in% cgeno$V2,], data.frame(FID=0, IID=IID, Total.PE, Total.THC, Total.dhCer))
write.table(fenos, "cirugia_bar/fenotipos/para_asoc/PE_THC_dhCer.phen", quote=FALSE, row.names=FALSE)
