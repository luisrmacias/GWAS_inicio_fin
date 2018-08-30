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
-Creamos genotipos temporales para cálculo de indicadores poblacionales

##carpeta poblacional
acomodo.txt
-Eliminamos archivos temporales
fusion_ancestrales.txt
-Unimos genotipos de población de estudio con poblaciones de referencias
calculo_ancestria_pca.txt
calculo_semejanza_genética.txt
gráfica_componente.R
cajas_ancestría.R
