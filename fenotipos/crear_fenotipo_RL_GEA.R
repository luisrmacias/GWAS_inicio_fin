library(MASS)
library(stringr)
library(foreign)

setwd("~/analisis_asociacion/adultos_mega")
load("cirugia_bar/fenotipos/histologia_revisada.RData")
genotipos <- "cirugia_bar/imputados/fusionados/fusion_agosto2018.QC"
ids <- read.table(paste0(genotipos, ".fam"))
ids <- ids[grepl("RL|GEA", ids$V2),]
ids$V2[!(ids$V2 %in% morbidos$IID)]
##RL_109, RL_252, RL_253, RL_258, RL_282, RL436 y RL_441 se genotiparon tanto con y sin guión bajo. Se cambia IID (se retira guión bajo).
morbidos$IID[grep("252|253|258|109|282|436|441", morbidos$IID)] <- gsub("_", "", morbidos$IID[grep("252|258|253|109|282|436|441", morbidos$IID)])

recientes <- read.spss("cirugia_bar/fenotipos/SPSS/metabolómica_2aetapa.sav", to.data.frame=TRUE)
recientes$ID <- str_trim(recientes$ID)
recientes$FECHA_NACIMIENTO <- as.Date(recientes$FECHA_NACIMIENTO, "%d/%m/%Y")
names(recientes)[!names(recientes) %in% names(morbidos)]
names(recientes)[grep("[0-9]\\.[1-5]", names(recientes))] <- gsub("\\.", "", grep("[0-9]\\.[1-5]", names(recientes), value=TRUE))
names(recientes)[grep("sexo", names(recientes))] <- "GENERO"
recientes$GENERO <- as.numeric(recientes$GENERO)*-1+3
recientes$IID <- recientes$ID
recientes <- recientes[-(212:213),]
recientes <- recientes[names(recientes) %in% names(morbidos)]
recientes <- recientes[,-grep("GENERO.1", names(recientes))]
recientes[pmatch(names(morbidos), names(recientes))]
agregar <- matrix(NA, nrow(recientes), ncol(morbidos[!names(morbidos) %in% names(recientes)]))
colnames(agregar) <- names(morbidos)[!names(morbidos) %in% names(recientes)]
recientes <- cbind(recientes, data.frame(agregar))

morbidos <- rbind(morbidos, recientes[pmatch(names(morbidos), names(recientes))])
rm(agregar, genotipos, ids, recientes)
save.image("cirugia_bar/fenotipos/histologia_revisada.RData")

grep("IID", names(morbidos))
##2
grep("EDAD", names(morbidos))
##33
grep("IMC", names(morbidos))
##40
grep("GENERO", names(morbidos))
##34

##primer archivo de covariables generales
covar <- read.table("fenotipos/GCTA/covar.cov", header=TRUE)
covar <- covar[!covar$IID %in% morbidos$IID,]
covar <- rbind(covar, with(morbidos[morbidos$IID!="",], data.frame(FID=0,  IID=IID, EDAD=EDAD, IMC=IMC)))
write.table(covar, "fenotipos/GCTA/covar.cov", row.names=FALSE, quote=FALSE)

##car_fuerarango <- c("C4OHC3DC", "C51", "C6", "C16", "C161", "C161OH", "C16OH", "C18", "C181OH", "C182", "C18OH")
##morbidos[,!names(morbidos) %in% car_fuerarango]
amino <- 183:193
carni <- 194:214
##3ARG raíz cuadrada/4CIT log/5GLY logaritmo/6ALA logaritmo/7LEU raíz cuadrada/8MET log/9PHE raíz cuadrada/10TYR raíz cuadrada/11VAL igual/12ORN logaritmo/13PRO raíz cuadrada
amino_trans <- sqrt(morbidos[amino[c(1,5,7,8,11)]])
names(amino_trans) <- paste("raiz2", names(amino_trans), sep="_")
amino_trans <- cbind(amino_trans, log(morbidos[amino[c(2:4,7,10)]]))
names(amino_trans)[6:ncol(amino_trans)] <- paste("log", names(amino_trans)[6:ncol(amino_trans)], sep="_")


raiz <- c("SA", "C0", "C5OHC4DC", "C5DCC6OH", "C5", "C102")
inverso <- c("C6DC", "C101")
loga <- c("C3", "C4", "C51", "C6", "C8", "C81", "C161", "C16OH", "C10", "C14", "C18")
raiz_inversa <- c("C2", "C4OHC3DC", "C16", "C12", "C121", "C141", "C142", "C14OH", "C181", "C182")

carni_trans <- sqrt(morbidos[,names(morbidos) %in% raiz])
names(carni_trans) <- paste("raiz2", names(carni_trans), sep="_")
carni_trans <- cbind(carni_trans, 1/morbidos[,names(morbidos) %in% inverso])
names(carni_trans)[names(carni_trans) %in% inverso] <- paste("uno_sobre", inverso, sep="_")
carni_trans <- cbind(carni_trans, log(morbidos[,names(morbidos) %in% loga]))
names(carni_trans)[names(carni_trans) %in% loga] <- paste("log", names(carni_trans)[names(carni_trans) %in% loga], sep="_")
carni_trans <- cbind(carni_trans, 1/sqrt(morbidos[,names(morbidos) %in% raiz_inversa]))
names(carni_trans)[names(carni_trans) %in% raiz_inversa] <- paste("uno_sobre_raiz2", names(carni_trans)[names(carni_trans) %in% raiz_inversa], sep="_")

write.table(data.frame(FID=0, IID=morbidos$IID[complete.cases(amino_trans)], amino_trans[complete.cases(amino_trans),]), "cirugia_bar/fenotipos/para_asoc/aminoacidos_transformados.phen", row.names=FALSE, quote=FALSE)
write.table(data.frame(FID=0, IID=morbidos$IID[complete.cases(carni_trans)], amino_trans[complete.cases(carni_trans),]), "cirugia_bar/fenotipos/para_asoc/carnitinas_transformadas.phen", row.names=FALSE, quote=FALSE)

morbidos <- cbind(morbidos, amino_trans, carni_trans)
rm(carni, amino, inverso, loga, raiz, raiz_inversa)
save.image("cirugia_bar/fenotipos/histologia_revisada")
