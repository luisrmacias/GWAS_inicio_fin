##este guión se sale del flujo de trabajo.
##tenemos listado de SNP asociados a HDL en metanálisis, buscamos frecuencias, asociación con HDL y genotipos tanto en niños como en adultos
##Del listado faltan 4 SNP y un quinto tiene otro nombre
##copiamos listado de Excel y pegamos
cd ~/analisis_asociacion/convivencias_medix_qxob
cat > conjunto_niños_adultos/listas/SNPS/arti_lipidos.txt
##cambiamos rs35980001 por rs35980001,rs397821255 en nano
##buscamos 4 faltantes
##buscamos fusionar los genotipos de niños y adultos
adultos=adultos_mega/cirugia_bar/imputados/fusionados/adultos.082019
convi=niños/imputados/fusionados/fusion_noviembre2018
salida=conjunto_niños_adultos/genotipos/niños_adultos/
plink --bfile $adultos --bmerge $convi --allow-no-sex --make-bed --out $salida/niños_adultos_072020
##no hubo necesidad de cambiar hebras (¡Wow!) 

cd conjunto_niños_adultos/
grep -wf listas/SNPS/arti_lipidos.txt genotipos/niños_adultos/niños_adultos_072020.bim > listas/SNPS/disponibles_hdl.tmp
cut -f2 listas/SNPS/disponibles_hdl.tmp | diff  - listas/SNPS/arti_lipidos.txt 
cut -f2 listas/SNPS/disponibles_hdl.tmp | diff - listas/SNPS/arti_lipidos.txt | grep rs | tr -d "> " > listas/SNPS/faltantes.tmp
rm listas/SNPS/disponibles_hdl.tmp
#encontramos los genotipos en la carpeta 15 placas
plink --extract listas/SNPS/faltantes.tmp --bfile ../adultos_mega/calidad/15placas/15placas --make-bed --out genotipos/snp_perdidos_hdl_adultos
plink --extract listas/SNPS/faltantes.tmp --bfile ../niños/calidad/MEGAex_cluster/niños --make-bed --out genotipos/snp_perdidos_hdl_niños
awk '{print $1,$2,0,$2}' genotipos/snp_perdidos_hdl_niños.fam > listas/personas/cambio_ids_perdidos_hdl.tmp
plink --bfile genotipos/snp_perdidos_hdl_niños --update-ids listas/personas/cambio_ids_perdidos_hdl.tmp --make-bed --out genotipos/snp_perdidos_hdl_niños
awk '{print $1,$2,0,$2}' genotipos/snp_perdidos_hdl_adultos.fam > listas/personas/cambio_ids_perdidos_hdl.tmp
plink --bfile genotipos/snp_perdidos_hdl_adultos --update-ids listas/personas/cambio_ids_perdidos_hdl.tmp --make-bed --out genotipos/snp_perdidos_hdl_adultos
plink --bfile ../snp_perdidos_hdl_adultos --bmerge ../snp_perdidos_hdl_niños --allow-no-sex --make-bed --out snp_perdidos_hdl
plink --bfile niños_adultos_072020 --bmerge ../snp_perdidos_hdl --allow-no-sex --make-bed --out niños_adultos_072020
plink --bfile niños_adultos_072020 --geno 0.05 --maf 0.01 --mind 0.05 --make-bed --out niños_adultos_072020.QC
cd ../..
grep -wf listas/SNPS/arti_lipidos.txt genotipos/niños_adultos/niños_adultos_072020.QC.bim > listas/SNPS/disponibles_hdl.tmp
cut -f2 listas/SNPS/disponibles_hdl.tmp | diff  - listas/SNPS/arti_lipidos.txt | grep rs | tr -d "> " > listas/SNPS/faltantes.tmp
##vuelven a faltar 5, eliminados en QC
cd genotipos/niños_adultos
plink --bfile niños_adultos_072020 --extract ../../listas/SNPS/faltantes.tmp --make-bed --out snp_sueltos
plink --bfile niños_adultos_072020.QC --bmerge snp_sueltos --allow-no-sex --make-bed --out niños_adultos_072020.QC
plink --bfile niños_adultos_072020 --freq --extract ../../listas/SNPS/arti_lipidos.txt --out asociados_hdl
##terminamos haciendo el cálculo de frecuencia y la extracción de genotipos sin QC
##pasamos de muuk a tobago e integramos información a hoja de cálculo
##en tobago
cd /home/lr/Documentos/ugepas/analisis_asociacion/agrupado/asociacion
scp -P7637 lmacias@muuk.inmegen.gob.mx:/export/home/lmacias/analisis_asociacion/convivencias_medix_qxob/conjunto_niños_adultos/genotipos/niños_adultos/asociados_hdl.frq lipidos/
##integramos en R
library(openxlsx)
hdl <- loadWorkbook("lipidos/Table1. Metaanalisis_ultimo_HDL216_Julio_2020.xlsx")
hdl <- read.xlsx(hdl) 
##verificamos que SNP en adendum que no tienen las mismas columnas están incluídos en columnas
hdl$CHR[66:72] %in% hdl$rs.ID
##todos TRUE
##quitamos esas líneas
hdl <- hdl[complete.cases(hdl$rs.ID),]

freq <- read.table("lipidos/asociados_hdl.frq", header=TRUE)
##verificamos alelos
ver_al <- data.frame(hdl$Allele.Test, A1=freq$A1[pmatch(hdl$rs.ID, freq$SNP)], hdl$Allele2, A2=freq$A2[pmatch(hdl$rs.ID, freq$SNP)])
ver_al$hdl.Allele.Test <- toupper(as.character(ver_al$hdl.Allele.Test))
ver_al$hdl.Allele2 <- toupper(as.character(ver_al$hdl.Allele2))
table(ver_al$hdl.Allele.Test==ver_al$A1 & ver_al$hdl.Allele2==ver_al$A2, exclude=NULL)
no_conciden <- which(!(ver_al$hdl.Allele.Test==ver_al$A1 & ver_al$hdl.Allele2==ver_al$A2))
##42 coincidentes, 16 no
##tengo el cerebro trabado y no se me ocurré cómo automatizar la diferenciación entre cambio de hebra y de alelos. Lo haré uno por uno.
difer <- c("hebra", "alelos", "alelos", "hebra", "alelos", "hebra", "hebra", "alelos", "alelos", "alelos", "alelos", "alelos", "alelos", "alelos", "alelos", "alelos")
rm(difer, no coinciden)
##después de muchas vueltas descubrí cómo identificar los que se tienen que cambiar
difer <- ver_al$hdl.Allele.Test==ver_al$A2 & ver_al$hdl.Allele2==ver_al$A1
hdl$Allele.Test.frequency[which(!difer)] <- freq$MAF[pmatch(hdl$rs.ID[which(!difer)], freq$SNP)]
hdl$Allele.Test.frequency[which(difer)] <- 1-freq$MAF[pmatch(hdl$rs.ID[which(difer)], freq$SNP)]
hdl$Allele.2.frequency <- 1-hdl$Allele.Test.frequency

##regresamos a muuk/bash para agregar un SNP que encontramos
##hacemos metanálisis de asociación con HDL con meta (R), nos basamos en /export/home/lmacias/analisis_asociacion/convivencias_medix_qxob/conjunto_niños_adultos/guiones/metanalisis/metanalisis_urico.R y en /export/home/lmacias/analisis_asociacion/convivencias_medix_qxob/conjunto_niños_adultos/guiones/lineal_mixto.R

meta <- read.table("lipidos/beta_se_hdl.txt", header=TRUE)
hdl$rs.ID[!hdl$rs.ID %in% rownames(meta)]
##hay que cambiar dos rs
rownames(meta) <- gsub("\\.", ",", rownames(meta))
hdl[,9] <- meta$meta_beta[pmatch(hdl$rs.ID, rownames(meta))]
hdl$SE <- meta$meta_se[pmatch(hdl$rs.ID, rownames(meta))
rownames(meta)[grep("rs35980001", rownames(meta))] <- "rs35980001"
write.xlsx(hdl, "lipidos/Table1. Metaanalisis_ultimo_HDL222_Julio_2020.xlsx")
