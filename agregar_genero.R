##creación de lista de identificadores individuales

##cut -f 1,2  fuentes/MEGA_EX_Canizlaes_U15_118_GS/PLINK/U15_118.ped > listas/eliminacion_QC/segundos_controles/ids.txt
ids <- read.table("calidad/MEGAex_cluster/agosto_2017.fam")
trios <- read.table("~/CANDELA/haplotipos/indigenas/MEGA/fusion/OMNI_MEGA.fam")
load("fenotipos/combinados/muy_combinado.RData")
##ICN$Sexo <- factor(ICN$Sexo, labels=c("Mujer", "Hombre"))


llc <- read.csv("../leucemia/listas/LLC.csv")
llc$ID[nchar(llc$ID)==6] <- gsub("_", "_0", llc$ID[nchar(llc$ID)==6])
llc$ID[nchar(llc$ID)==5] <- gsub("_", "_00", llc$ID[nchar(llc$ID)==5])
llc$SEX <- llc$SEX + 1
llc$ID <- paste0("LLC_", llc$ID)
llc$ID[nchar(llc$ID)==6] <- gsub("LLC_", "LLC_0", llc$ID[nchar(llc$ID)==6])
llc$ID[nchar(llc$ID)==5] <- gsub("LLC_", "LLC_00", llc$ID[nchar(llc$ID)==5])
##tenemos datos de GC, ICN, LL, 36 CON, 36 RED, RL, 30 SOL y 30 TER
sexos <- cbind(ids[,1:2], NA)
names(sexos) <- c("FID", "IID", "SEX")
sexos$SEX <- casosycontroles$Sexo[pmatch(ids$V2, casosycontroles$IID)]

sexos$SEX[sexos$IID %in% llc$ID] <- llc$SEX[pmatch(sexos$IID[sexos$IID %in% llc$ID], llc$ID)]
temp <- as.character(sexos$IID)

sexos$SEX[temp %in% morbidos$IID] <- morbidos$GENERO[pmatch(temp[temp %in% morbidos$IID], morbidos$IID)]
sexos$SEX[temp %in% trios$V2] <- trios$V5[pmatch(temp[temp %in% trios$V2], trios$V2)]
sexos$SEX[is.na(sexos$SEX)] <- 0
write.table(sexos, "listas/integrar_sexos_primeros.txt", row.names=FALSE, quote=FALSE)
