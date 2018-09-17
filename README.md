# GWAS_inicio_fin

Guiones para manejo de datos, control de calidad y asociación genómica. 

##carpeta preparación
- cambio_ids_filtrado_variantes_duplicadas.txt
Convertimos de formato de texto a binario, identificamos y filtramos variantes duplicadas, cambiamos identificadores de Illumina a rs. Depende de sitios_dup_bim.R
- agregar_genero.R o agregar_genero_conpreexistente.R
Creamos archivo para cambio de sexo en archivo .fam
- agregar_datos_muestra.txt
Incluimos dato de sexo y estado de caso o control para obesidad en archivo .fam

##carpeta QC
fusionar_grupos.txt
-unimos genotipos de los grupos/poblaciones de estudio
identificacion_eliminacion_individuos.txt
-Tras descartar SNP de baja calidad y no informativos identificamos y eliminamos individuos con tasa de genotipificación baja, sexo discordante y duplicados
evaluación_homogeneidad.txt
-Evaluamos si hay variantes cuya frecuencia alélica difiere por grupo de asignación de genotipos o por estudio a niveles sospechosos de errores en la genotipificación
filtrado_variantes_preimput.txt
-Creamos genotipos útiles para cálculo de indicadores poblacionales, para imputación usaremos genotipos con un filtro de frecuencia menos estricto

##carpeta poblacional
acomodo.txt
-Eliminamos archivos temporales
fusion_ancestrales.txt
-Unimos genotipos de población de estudio con poblaciones de referencia 
calculo_ancestria_pca.txt
calculo_semejanza_genética.txt
gráfica_componente.R
cajas_ancestría.R
-me parece los últimos 4 guiones tienen un nombre lo suficientemente claro como para requerir explicación

##carpeta imputación
conversion_vcf.txt
-pasamos genotipos de formato plink a vcf
faseo_desdevcf.txt
-Paso del faseo e imputación. Supone que ya contamos con genotipos de referencia (1000 genomas + indígenas) faseados.
filtrado_calidad_imput.txt
-removemos variantes de baja frecuencia (en este momento <0.01 >) y variantes cuya imputación es poco confiable.
reunion_faseados_resto.txt
-agregamos a genotipos faseados variantes que no se encontraron en las referencias. Estas variantes permanecen sin fasear pero para fines de asociación son útiles
controlcalidad_imputados.txt
-filtrado final que genera genotipos listos para asociación

##carpeta fenotipos
archivo_covariable_sexos_estudio.txt
-GCTA requiere que las covariables contínuas y categóricas se encuentren en archivos separados. En este caso creamos las covariables categóricas.
crear_fenotipo_RL_GEA.R
-Ejemplo de archivo de fenotipos, tal vez sea difícil generalizar (desde el punto de vista computacional) esta parte y se tenga que integrar un archivo de fenotipos dependiendo del objetivo de investigación
