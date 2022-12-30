set.seed(123)
     Geno.mtx = matrix(rbinom(10000,2,0.3),1000,10)
     rownames(Geno.mtx) = egData$IID
     colnames(Geno.mtx) = paste0("rs",1:10)
     chrVec = "1"  # equivalant to chrVec = rep("1", ncol(Geno.mtx))
     outPOLMM = POLMM(objNull, Geno.mtx, chrVec)
     
     outPOLMM
     round(as.numeric(outPOLMM$pval.spa),2)

objpolmm2 <- POLMM.plink(objNull, PlinkFile, chrVec.plink = 1:22, output.file="adversos_tnf2_sinplac");  # if you only want to analyze variants in chromosomes 1-4 
     outPOLMM = read.table(outFile, header = T)
     head(outPOLMM)

