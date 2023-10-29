
###### CÓDIGO DE R, DE LEYRE HUARTE, PARA LA REALIZACIÓN DE SU TRABAJO DE FINAL DE MÁSTER. 

#.libPaths("C:/Users/huart/OneDrive/Escritorio/R_Library")


### PASO 1º Instalar y cargar paquetes necesarios para la actividad

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("GOSeSim")

BiocManager::install("TCGAbiolinks")

if (!requireNamespace("BiocManager", quietly=TRUE)){install.packages("BiocManager")}
BiocManager::install("biomaRt")
BiocManager::install("ComplexHeatmap")
BiocManager::install("pathview")
install.packages("dplyr")         # Paquete para manipulación de datos
install.packages("biomaRt")       # Paquete para acceso a datos biológicos desde BioMart
install.packages("SummarizedExperiment")  # Paquete para almacenar datos omics
install.packages("AnnotationHub")  # Paquete para obtener metadatos de datos biológicos
install.packages("GEOquery")      # Paquete para descargar datos de GEO (Gene Expression Omnibus)
install.packages("survival")      # Paquete para análisis de supervivencia

if (!requireNamespace("tidyverse", quietly=TRUE)){install.packages("tidyverse")}
if (!requireNamespace("tidyverse", quietly=TRUE)){BiocManager::install("biomaRt")}
if (!requireNamespace("tidyverse", quietly=TRUE)){BiocManager::install("enrichR")}
if (!requireNamespace("tidyverse", quietly=TRUE)){BiocManager::install("ComplexHeatmap")}
if (!requireNamespace("tidyverse", quietly=TRUE)){BiocManager::install("pathview")}
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GOSemSim")
if (!requireNamespace("NOISeq", quietly = TRUE)) {install.packages("NOISeq")}


knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)
version
packageVersion("TCGAbiolinks")

# Cargar librerías necesarias:
library(GOSemSim)
library(NOISeq)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(TCGAbiolinks)
library(biomaRt)
library(AnnotationDbi)
library(survival)
library(tidyverse)


# Establecemos la ruta para nuestro directorio de trabajo
getwd()
setwd("C:/Users/huart/OneDrive/Escritorio/TFM")



### PASO 2: TCGAbiolinks: Búsqueda en la base de datos TCGA:
# Podemos buscar fácilmente datos de TCGA utilizando la función GDCquery.
# Utilizando un resumen de filtros como los utilizados en el portal TCGA, la función trabaja con los siguientes argumentos:

knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

BRCA_clinical_data <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab",
)

GDCdownload(BRCA_clinical_data)
clinical.BCRtab.all <- GDCprepare(BRCA_clinical_data)
View(clinical.BCRtab.all)

###########################################################################################################################
# En total hay 9 :

####clinical_follow_up_v4.0_brca (718 entries, 13 total columns)
clinical_follow_up_v4.0_brca <- as.data.frame(clinical.BCRtab.all$clinical_follow_up_v4.0_brca)
View(clinical_follow_up_v4.0_brca)
lista_de_variables <- names(clinical_follow_up_v4.0_brca)
print(lista_de_variables)

#clinical_follow_up_v4.0_brca %in% c("[Discrepancy]", "[Not Applicable]", "[Not Available]", "[Unknown]") |
View(clinical_follow_up_v4.0_brca)


###clinical_follow_up_v2.1_brca (525 entries, 45 total columns)
clinical_follow_up_v2.1_brca <- as.data.frame(clinical.BCRtab.all$clinical_follow_up_v2.1_brca)
View(clinical_follow_up_v2.1_brca)
lista_de_variables <- names(clinical_follow_up_v2.1_brca)
print(lista_de_variables)

clinical_follow_up_v2.1_brca <- clinical_follow_up_v2.1_brca %>%
  mutate_all(~ifelse(. %in% c("[Not Available]", "[Not Applicable]"), NA, .))
View(clinical_follow_up_v2.1_brca)

###clinical_follow_up_v1.5_brca ( 116 entries, 63 total columns)
clinical_follow_up_v1.5_brca <- as.data.frame(clinical.BCRtab.all$clinical_follow_up_v1.5_brca)
View(clinical_follow_up_v1.5_brca)

###clinical_patient_brca (1,099 entries, 112 total columns)
clinical_patient_brca <- as.data.frame(clinical.BCRtab.all$clinical_patient_brca)
View(clinical_patient_brca)

###clinical_radiation_brca (620 entries, 18 total columns)
clinical_radiation_brca <- as.data.frame(clinical.BCRtab.all$clinical_radiation_brca)
View(clinical_radiation_brca)
######################################################################################################


#################### La INFORMACIÓN que más nos interesa es la de RESPUESTA A TRATAMIENTO:####################
Respuesta_tratamiento <- as.data.frame(clinical.BCRtab.all$clinical_drug_brca)
names(Respuesta_tratamiento)
View(Respuesta_tratamiento)

# Necesito saber que medicamneto es el mas frecuente, para ver cual podria resultarme útil: ANTES DE FILTRAR
frecuencia_medic <- table(Respuesta_tratamiento$pharmaceutical_therapy_drug_name)

# Ordenar en orden descendente para encontrar los más comunes
frecuencia_medic <- sort(frecuencia_medic, decreasing = TRUE)
View(frecuencia_medic)

#A partir de aqui podemos explorar las diferentes colummnas del dataframe Respuesta_tratamiento:

## Para ver el número de pacientes que hay con cada categoría de tratamiento utilizamos:
table(Respuesta_tratamiento$pharmaceutical_therapy_type) 

# Y nuestra variable más importante que pretendemos predecir, podemos obtener las diferentes respuestas al tratamiento:
table(Respuesta_tratamiento$treatment_best_response)
table(Respuesta_tratamiento$pharmaceutical_therapy_drug_name)

########################    Crear un DATAFRAME con las columnas requeridas     ########################   

Variables_utiles_tratamiento <- data.frame(
  #bcr_patient_uuid = Respuesta_tratamiento$bcr_patient_uuid,
  #bcr_drug_uuid = Respuesta_tratamiento$bcr_drug_uuid,
  bcr_patient_barcode = Respuesta_tratamiento$bcr_patient_barcode,
  pharmaceutical_therapy_drug_name = Respuesta_tratamiento$pharmaceutical_therapy_drug_name,
  #pharmaceutical_therapy_type = Respuesta_tratamiento$pharmaceutical_therapy_type,
  treatment_best_response = Respuesta_tratamiento$treatment_best_response
)
# Ver el nuevo dataframe
View(Variables_utiles_tratamiento)

# Eliminar las filas con valores NO deseados
Variables_utiles_tratamiento_filt <- Variables_utiles_tratamiento[!(Variables_utiles_tratamiento$pharmaceutical_therapy_drug_name == "[Not Available]" |
                                                                      Variables_utiles_tratamiento$treatment_best_response %in% c("[Discrepancy]", "[Not Applicable]", "[Not Available]", "[Unknown]", "[CDE_ID:2857291]")
), ]

### Categorizamosla variable treatment_best_response

Variables_utiles_tratamiento_filt <- Variables_utiles_tratamiento_filt %>%
  mutate(Grado_Respuesta = case_when(
    treatment_best_response %in% c("Complete Response", "Partial Response") ~ "Respondedor",
    treatment_best_response %in% c("Stable Disease", "Clinical Progressive Disease") ~ "No respondedor",
    TRUE ~ NA_character_  # En caso de que haya algún valor inesperado en treatment_best_response
  ))


medicamento_mas_comun <- table(Variables_utiles_tratamiento_filt$pharmaceutical_therapy_drug_name) # Cytoxan, 
# Ordenar en orden descendente para encontrar los más comunes
medicamento_mas_comun <- sort(medicamento_mas_comun, decreasing = TRUE)
View(medicamento_mas_comun)

### Verificar si hay duplicaods
duplicados <- duplicated(Variables_utiles_tratamiento_filt$bcr_patient_barcode)


### PASO 3º: TCGAbiolinks: Búsqueda en la base de datos TCGA de los DATOS DE EXPRESIÓN 

########################    DATOS DE EXPRESIÓN SANOS #####################

Expresion_data_sano <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"
)
GDCdownload(Expresion_data_sano)
expression_data_sano <- GDCprepare(Expresion_data_sano)

print(expression_data_sano)

########################    DATOS DE EXPRESION TUMOR ########################    

Expresion_data_tumor <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type =  "Primary Tumor"
)
GDCdownload(Expresion_data_tumor)
expression_data_tumor <- GDCprepare(Expresion_data_tumor)
print(expression_data_tumor)

########################    CREAR UN DATAFRAME NUEVO CON LA INFORMACIÓN COLDATA DE LOS DATOS DE EXPRESION DE TUMOR  ########################
# TUMOR
coldata_tumor = colData(expression_data_tumor)
coldata_tumor = as.data.frame(coldata_tumor)

#SANO
coldata_sano = colData(expression_data_sano)
coldata_sano = as.data.frame(coldata_sano)

########################    CREAR UN DATAFRAME NUEVO CON LOS COUNTS BRUTOS DE LOS DATOS DE EXPRESION DE TUMOR  ########################
# TUMOR
rawcounts_tumor_matrix = assay(expression_data_tumor) 
rawcounts_tumor_df = as.data.frame(rawcounts_tumor_matrix)
# SANO
rawcounts_sano_matrix = assay(expression_data_sano)
rawcounts_sano_df = as.data.frame(rawcounts_sano_matrix)


###########################   PASO 4º: NORMALIZACIÓN

## Para realizar la NORMALIZACIÓN  TMM (Trimmed Mean of M-values) en los datos de expresión génica, necesitas obtener los recuentos brutos (raw counts) de los datos. 
##  Los recuentos brutos son el número de lecturas de secuenciación que se asignan a cada gen en cada muestra.

# El argumento "unstranded" se utiliza para acceder a los recuentos brutos sin tener en cuenta la orientación de las secuencias. 

# TUMOR
rawcounts_tumor_matrix = assay(expression_data_tumor) 
rawcounts_tumor_df = as.data.frame(rawcounts_tumor_matrix)
# SANO
rawcounts_sano_matrix = assay(expression_data_sano)
rawcounts_sano_df = as.data.frame(rawcounts_sano_matrix)

##### NORMALIZACIÓN TMM

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("NOISeq")
library(NOISeq)

datos_expresion_norm_sano<- tmm(rawcounts_sano_matrix)
datos_expresion_norm_sano <-  log2(datos_expresion_norm_sano + 1)

datos_expresion_norm_tumor<- tmm(rawcounts_tumor_matrix)
datos_expresion_norm_tumor <-  log2(datos_expresion_norm_tumor + 1)


#### ENSG to Symbol ####

#Function to convert ENSG to Symbol

####### SANO #############

convert.ENSG.Symbol<-function(genes){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
  return(G_list)
}

expr.symbol_sano <- datos_expresion_norm_sano

rownames(expr.symbol_sano) <- gsub("\\..*", "", rownames(expr.symbol_sano))

##### # converting ensenmbl gene ids to huugo gymbols using Biomart package

conversion.table <- convert.ENSG.Symbol(rownames(expr.symbol_sano))

conversion.inter <- intersect(conversion.table[-which(conversion.table$hgnc_symbol==""),]$ensembl_gene_id, rownames(expr.symbol_sano))
conversion.table2 <- conversion.table[which(conversion.table$ensembl_gene_id %in% conversion.inter),]
dups_to_remove <- conversion.table2$ensembl_gene_id[duplicated(conversion.table2$ensembl_gene_id)]
dups_to_remove
conversion.table2 <- conversion.table2[!conversion.table2$ensembl_gene_id%in%dups_to_remove,]
rownames(conversion.table2) <- conversion.table2$ensembl_gene_id

conversion.inter2 <- intersect(rownames(conversion.table2), rownames(expr.symbol_sano))
conversion.table2 <- conversion.table2[conversion.inter2,]

expr.symbol_sano <- expr.symbol_sano[rownames(conversion.table2),]
rownames(expr.symbol_sano) <- conversion.table2$hgnc_symbol

any(is.na(rownames(conversion.table2)))
any(rownames(conversion.table2)=="")
any(is.na(rownames(expr.symbol_sano)))
any(rownames(expr.symbol_sano)=="")


table(duplicated(rownames(expr.symbol_sano))) # some duplicated gene symbols
expr.symbol_sano <- expr.symbol_sano[!duplicated(rownames(expr.symbol_sano)),]


datos_expresion_norm_sano <- expr.symbol_sano


####### TUMOR #############

convert.ENSG.Symbol<-function(genes){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
  return(G_list)
}

expr.symbol_tumor <- datos_expresion_norm_tumor

rownames(expr.symbol_tumor) <- gsub("\\..*", "", rownames(expr.symbol_tumor))

##### # converting ensenmbl gene ids to huugo gymbols using Biomart package

conversion.table <- convert.ENSG.Symbol(rownames(expr.symbol_tumor))

conversion.inter <- intersect(conversion.table[-which(conversion.table$hgnc_symbol==""),]$ensembl_gene_id, rownames(expr.symbol_tumor))
conversion.table2 <- conversion.table[which(conversion.table$ensembl_gene_id %in% conversion.inter),]
dups_to_remove <- conversion.table2$ensembl_gene_id[duplicated(conversion.table2$ensembl_gene_id)]
dups_to_remove
conversion.table2 <- conversion.table2[!conversion.table2$ensembl_gene_id%in%dups_to_remove,]
rownames(conversion.table2) <- conversion.table2$ensembl_gene_id

conversion.inter2 <- intersect(rownames(conversion.table2), rownames(expr.symbol_tumor))
conversion.table2 <- conversion.table2[conversion.inter2,]

expr.symbol_tumor <- expr.symbol_tumor[rownames(conversion.table2),]
rownames(expr.symbol_tumor) <- conversion.table2$hgnc_symbol

any(is.na(rownames(conversion.table2)))
any(rownames(conversion.table2)=="")
any(is.na(rownames(expr.symbol_tumor)))
any(rownames(expr.symbol_tumor)=="")


table(duplicated(rownames(expr.symbol_tumor))) # some duplicated gene symbols
expr.symbol_tumor <- expr.symbol_tumor[!duplicated(rownames(expr.symbol_tumor)),]


datos_expresion_norm_tumor <- expr.symbol_tumor



##### UTILIZA LAS PRIMERAS 20 COLUMNAS DEBIDO AL TAMAÑO
primeras_20_columnas_sano <- datos_expresion_norm_sano[, 1:20]
primeras_20_columnas_tumor <- datos_expresion_norm_tumor[, 1:20]



#############################   PASO 5º: CÁLCULO DE M-SCORES MEDIANTE PATHMED #####################

# El primer paso consiste en calcular las puntuaciones M de los conjuntos de datos de referencia. 
# Las M-scores miden el nivel de alteración de las vías moleculares entre los pacientes y los controles sanos. Este enfoque permite comparar los resultados obtenidos con las puntuaciones M de los nuevos datos, incluso cuando los genes medidos no se solapan completamente.

install.packages("devtools")
library(devtools)
devtools::install_github("jordimartorell/pathMED")
library(pathMED)
library(FactoMineR)
library(factoextra)

?getMscoresRef
?getML
?GDCquery


#### Prepare data for pathMED ####

save.image("C:/Users/huart/OneDrive/Escritorio/TFM/TFM_environment.RData")


### Organizar la estructura de los datos y guardarlos en un archivo RData.
# 20 prim
my_refData_20prim <- list(
  dataset1 = list(
    Disease = primeras_20_columnas_tumor,
    Healthy = primeras_20_columnas_sano)
)
# total
my_refData <- list(
  dataset1 = list(
    Disease = datos_expresion_norm_tumor,
    Healthy = datos_expresion_norm_sano)
)

save(my_refData, file = "C:/Users/huart/OneDrive/Escritorio/TFM/my_refData.RData")
save(my_refData_20prim, file = "C:/Users/huart/OneDrive/Escritorio/TFM/my_refData_20prim.RData")

#### Apply pathMED reference ####

#remove.packages("pathMED")
#devtools::install_github("jordimartorell/pathMED")
library(pathMED)
library(FactoMineR)
library(factoextra)
library(dplyr)

# 20 prim
my_refMscore_20prim <- getMscoresRef(data=my_refData_20prim,
                                     genesets = "reactome",
                                     cores = 8)

my_relevantPaths_20prim <- diseasePaths(MRef=my_refMscore_20prim, Pcutoff = 0.001)

save(my_refMscore_20prim, file = "C:/Users/huart/OneDrive/Escritorio/TFM/my_refMscore_20prim.RData")
save(my_relevantPaths_20prim, file = "C:/Users/huart/OneDrive/Escritorio/TFM/my_relevantPaths_20prim.RData")

?diseasePaths
##total
my_refMscore <- getMscoresRef(data=my_refData,
                              genesets = "reactome",
                              cores = 8)

my_relevantPaths <- diseasePaths(MRef=my_refMscore)

#----------------------------------------------------------------

###### PREPARACIÓN DE LA INFORMACION CLINICA PARA CADA MEDICAMENTO Y NECESARIA PARA ENTRENAR LOS MODELOS #################

##########################################################
############ FÁRMACO CYTOXAN/CYCLOPHOSPHAMIDE  ###########
##########################################################


####### Filtramos los datos para incluir sólo a los pacientes tratados con el fármaco deseado ############

Cytoxan_pacientes <- subset(Variables_utiles_tratamiento_filt, pharmaceutical_therapy_drug_name %in% c("Cytoxan", "Cyclophosphamide", "CYTOXAN", "cyclophosphamide", "cyclophosphamid", "Cyclophospamide", "doxorubicine cyclophosphamide tamoxifen", "doxorubicine+cyclophosphamide+tamoxifen", "adrimicin+cyclophosphamide", "adriamycin+cyclophosphamide
", "cyclophosphamidum", "Cyclophosphane", "Cytoxan and Taxotere"))

### A pharmaceutical_therapy_drug_name, le asignamos el mismo valor, "CYTOXAN" y comprobamos que no existan duplicados, Y SI LOS HAY LOS ELIMINAMOOS Y POR ULTIMO ESCOGEMOS LOS 20 PACIENTES:

hay_duplicados <- any(duplicated(Cytoxan_pacientes$bcr_patient_barcode))
print(hay_duplicados)
Cytoxan_pacientes$pharmaceutical_therapy_drug_name <- "CYTOXAN"
Cytoxan_pacientes <- Cytoxan_pacientes %>% distinct(bcr_patient_barcode, .keep_all = TRUE)
Cytoxan_pacientes_20prim <- head(Cytoxan_pacientes, 20)



################ OBTENEMOS TODOS LOS BARCODE DE PACIENTES QUE HAN USADO CYTOXAN ############# 
Id_cytoxan_clin <- Cytoxan_pacientes$bcr_patient_barcode
Id_cytoxan_clin
Id_pacientes_expr <- coldata_tumor$patient
CYTOXAN_pacientes_coin <- intersect(Id_cytoxan_clin, Id_pacientes_expr)
CYTOXAN_pacientes_coin = as.data.frame(CYTOXAN_pacientes_coin) 

################ OBTENEMOS LOS 20 BARCODE DE PACIENTES QUE HAN USADO CYTOXAN ############# 

Id_cytoxan_clin_20 <- Cytoxan_pacientes_20prim$bcr_patient_barcode ####PASO UNO 
Id_cytoxan_clin_20
Id_pacientes_expr <- expression_data_tumor$patient
Id_pacientes_expr

CYTOXAN_pacientes_coincidentes <- intersect(Id_cytoxan_clin_20, Id_pacientes_expr)
CYTOXAN_pacientes_coincidentes = as.data.frame(CYTOXAN_pacientes_coincidentes) ##### los pacientes que estan tanto en clinicos, como en datos de expresión. 


###### Crear un nuevo dataframe con las columnas correspondientes a los ID de pacientes coincidentes, en coldata_tumor

expr_CYTOXAN <- coldata_tumor[coldata_tumor$patient %in% CYTOXAN_pacientes_coin$CYTOXAN_pacientes_coin, ]

expr_CYTOXAN_20prim <- coldata_tumor[coldata_tumor$patient %in% CYTOXAN_pacientes_coincidentes$CYTOXAN_pacientes_coincidentes, ]

###### 
Cytoxan_pacientes <- Cytoxan_pacientes[match(expr_CYTOXAN$patient, Cytoxan_pacientes$bcr_patient_barcode) ,]
rownames_temp <- rownames(expr_CYTOXAN)
rownames(expr_CYTOXAN) <- rownames(Cytoxan_pacientes)
rownames(Cytoxan_pacientes) <- rownames_temp

Cytoxan_pacientes_20prim <- Cytoxan_pacientes_20prim[match(expr_CYTOXAN_20prim$patient, Cytoxan_pacientes_20prim$bcr_patient_barcode) ,]
rownames_temp <- rownames(expr_CYTOXAN_20prim)
rownames(expr_CYTOXAN_20prim) <- rownames(Cytoxan_pacientes_20prim)
rownames(Cytoxan_pacientes_20prim) <- rownames_temp


########### OBTENEMOS LA MATRIZ DE EXPRESION DE LOS PACIENTES TRATADOS CON CYTOXAN ############# 

barcode <- expr_CYTOXAN$barcode
Cytoxan_expr <- datos_expresion_norm_tumor[, colnames(datos_expresion_norm_tumor) %in% barcode]

barcode_20prim <- expr_CYTOXAN_20prim$barcode
Cytoxan_expr_20prim <- datos_expresion_norm_tumor[, colnames(datos_expresion_norm_tumor) %in% barcode_20prim]

##########################################################
############ FÁRMACO ADRYAMICIN  ###########
##########################################################

####### Filtramos los datos para incluir sólo a los pacientes tratados con el fármaco deseado ############

Adryamicin_pacientes <- subset(Variables_utiles_tratamiento_filt, pharmaceutical_therapy_drug_name %in% c("Adriamycin", "ADRIAMYCIN", "adrimicin+cyclophosphamide", "adriamycin+cyclophosphamide"))

### A pharmaceutical_therapy_drug_name, le asignamos el mismo valor, "CYTOXAN" y comprobamos que no existan duplicados, Y SI LOS HAY LOS ELIMINAMOOS Y POR ULTIMO ESCOGEMOS LOS 20 PACIENTES:

hay_duplicados <- any(duplicated(Adryamicin_pacientes$bcr_patient_barcode))
print(hay_duplicados)
Adryamicin_pacientes$pharmaceutical_therapy_drug_name <- "ADRYAMICIN"
Adryamicin_pacientes <- Adryamicin_pacientes %>% distinct(bcr_patient_barcode, .keep_all = TRUE)
Adryamicin_pacientes_20prim <- head(Adryamicin_pacientes, 20)



################ TOTAL 
Id_Adryamicin_clin <- Adryamicin_pacientes$bcr_patient_barcode
Id_Adryamicin_clin
Id_pacientes_expr <- coldata_tumor$patient
ADRYAMICIN_pacientes_coin <- intersect(Id_Adryamicin_clin, Id_pacientes_expr)
ADRYAMICIN_pacientes_coin = as.data.frame(ADRYAMICIN_pacientes_coin) 

############## 2O PRIMEROS 
Id_Adryamicin_clin_20 <- Adryamicin_pacientes_20prim$bcr_patient_barcode ####PASO UNO 
Id_pacientes_expr <- expression_data_tumor$patient

########################    Encuentra los ID de pacientes que coinciden en ambos dataframes  ########################

Adryamicin_pacientes_coincidentes <- intersect(Id_Adryamicin_clin_20, Id_pacientes_expr)
Adryamicin_pacientes_coincidentes = as.data.frame(Adryamicin_pacientes_coincidentes) ##### los pacientes que estan tanto en clinicos, como en datos de expresión. 


###### Crear un nuevo dataframe con las columnas correspondientes a los ID de pacientes coincidentes, en coldata_tumor
expr_ADRYAMICIN <- coldata_tumor[coldata_tumor$patient %in% ADRYAMICIN_pacientes_coin$ADRYAMICIN_pacientes_coin, ]

expr_ADRYAMICIN_20prim <- coldata_tumor[coldata_tumor$patient %in% Adryamicin_pacientes_coincidentes$Adryamicin_pacientes_coincidentes, ]


Adryamicin_pacientes <- Adryamicin_pacientes[match(expr_ADRYAMICIN$patient, Adryamicin_pacientes$bcr_patient_barcode) ,]
rownames_temp <- rownames(expr_ADRYAMICIN)
rownames(expr_ADRYAMICIN) <- rownames(Adryamicin_pacientes)
rownames(Adryamicin_pacientes) <- rownames_temp

Adryamicin_pacientes_20prim <- Adryamicin_pacientes_20prim[match(expr_ADRYAMICIN_20prim$patient, Adryamicin_pacientes_20prim$bcr_patient_barcode) ,]
rownames_temp <- rownames(expr_ADRYAMICIN_20prim)
rownames(expr_ADRYAMICIN_20prim) <- rownames(Adryamicin_pacientes_20prim)
rownames(Adryamicin_pacientes_20prim) <- rownames_temp


##################

barcode <- expr_ADRYAMICIN$barcode
Adryamicin_expr <- datos_expresion_norm_tumor[, colnames(datos_expresion_norm_tumor) %in% barcode]

barcode_20prim <- expr_ADRYAMICIN_20prim$barcode
Adryamicin_expr_20prim <- datos_expresion_norm_tumor[, colnames(datos_expresion_norm_tumor) %in% barcode_20prim]



##########################################################
############ FÁRMACO TAXOTERE  ###########
##########################################################

####### Filtramos los datos para incluir sólo a los pacientes tratados con el fármaco deseado ############

Taxotere_pacientes <- subset(Variables_utiles_tratamiento_filt, pharmaceutical_therapy_drug_name %in% c("Taxotere", "ADRIAMYCIN",  "TAXOTERE", "taxotere", "Cytoxan and Taxotere"))

### A pharmaceutical_therapy_drug_name, le asignamos el mismo valor, "CYTOXAN" y comprobamos que no existan duplicados, Y SI LOS HAY LOS ELIMINAMOOS Y POR ULTIMO ESCOGEMOS LOS 20 PACIENTES:

hay_duplicados <- any(duplicated(Taxotere_pacientes$bcr_patient_barcode))
print(hay_duplicados)
Taxotere_pacientes$pharmaceutical_therapy_drug_name <- "TAXOTERE"
Taxotere_pacientes <- Taxotere_pacientes %>% distinct(bcr_patient_barcode, .keep_all = TRUE)
Taxotere_pacientes_20prim <- head(Taxotere_pacientes, 20)



################ TOTAL 
Id_Taxotere_clin <- Taxotere_pacientes$bcr_patient_barcode
Id_pacientes_expr <- coldata_tumor$patient
TAXOTERE_pacientes_coin <- intersect(Id_Taxotere_clin, Id_pacientes_expr)
TAXOTERE_pacientes_coin = as.data.frame(TAXOTERE_pacientes_coin) 

############## 2O PRIMEROS 
Id_Taxotere_clin_20 <- Taxotere_pacientes_20prim$bcr_patient_barcode ####PASO UNO 
Id_pacientes_expr <- expression_data_tumor$patient

########################    Encuentra los ID de pacientes que coinciden en ambos dataframes  ########################

Taxotere_pacientes_coincidentes <- intersect(Id_Taxotere_clin_20, Id_pacientes_expr)
Taxotere_pacientes_coincidentes = as.data.frame(Taxotere_pacientes_coincidentes) ##### los pacientes que estan tanto en clinicos, como en datos de expresión. 


###### Crear un nuevo dataframe con las columnas correspondientes a los ID de pacientes coincidentes, en coldata_tumor
expr_TAXOTERE <- coldata_tumor[coldata_tumor$patient %in% TAXOTERE_pacientes_coin$TAXOTERE_pacientes_coin, ]

expr_TAXOTERE_20prim <- coldata_tumor[coldata_tumor$patient %in% Taxotere_pacientes_coincidentes$Taxotere_pacientes_coincidentes, ]


Taxotere_pacientes <- Taxotere_pacientes[match(expr_TAXOTERE$patient, Taxotere_pacientes$bcr_patient_barcode) ,]
rownames_temp <- rownames(expr_TAXOTERE)
rownames(expr_TAXOTERE) <- rownames(Taxotere_pacientes)
rownames(Taxotere_pacientes) <- rownames_temp

Taxotere_pacientes_20prim <- Taxotere_pacientes_20prim[match(expr_TAXOTERE_20prim$patient, Taxotere_pacientes_20prim$bcr_patient_barcode) ,]
rownames_temp <- rownames(expr_TAXOTERE_20prim)
rownames(expr_TAXOTERE_20prim) <- rownames(Taxotere_pacientes_20prim)
rownames(Taxotere_pacientes_20prim) <- rownames_temp


##################

barcode <- expr_TAXOTERE$barcode
Taxotere_expr <- datos_expresion_norm_tumor[, colnames(datos_expresion_norm_tumor) %in% barcode]

barcode_20prim <- expr_TAXOTERE_20prim$barcode
TAxotere_expr_20prim <- datos_expresion_norm_tumor[, colnames(datos_expresion_norm_tumor) %in% barcode_20prim]

save.image("C:/Users/huart/OneDrive/Escritorio/TFM/TFM_environment.RData")


###################################################################
##################   ENTRENAMIENTO DEL MODELO #####################
###################################################################

################ CYTOXAN ###################

Cytoxan_MScores <- getMscores(genesets = my_relevantPaths,
                              Patient = Cytoxan_expr,
                              Healthy = datos_expresion_norm_sano,
                              cores = 12)  

pathways_annotated <- ann2term(Cytoxan_MScores)

Cytoxan_pacientes[Cytoxan_pacientes$Grado_Respuesta == "No respondedor", "Grado_Respuesta"] <- "NoRespondedor"

set.seed(123)
Cytoxan_model <- getML(expData=Cytoxan_MScores,
                       metadata=Cytoxan_pacientes,
                       var2predict="Grado_Respuesta",
                       models = methodsML(algorithms = c("rf", "nb", "glm", "lda", 'svmLinear'),
                                          tuneLength = 100),
                       # models = methodsML(algorithms = "all", tuneLength = 100),
                       subsamples = 1,
                       repeatsCV = 5,
                       foldsCV = 3,
                       positiveClass = "Respondedor")



#---------------------------------
################ ADRYAMICIN ###################

Adryamicin_MScores <- getMscores(genesets = my_relevantPaths,
                                 Patient = Adryamicin_expr,
                                 Healthy = datos_expresion_norm_sano,
                                 cores = 8)  


Adryamicin_pacientes[Adryamicin_pacientes$Grado_Respuesta == "No respondedor", "Grado_Respuesta"] <- "NoRespondedor"

set.seed(123)
Adryamicin_model <- getML(expData=Adryamicin_MScores,
                          metadata=Adryamicin_pacientes,
                          var2predict="Grado_Respuesta",
                          models = methodsML(algorithms = c("rf", "nb", "glm", "lda", 'svmLinear'),
                                             tuneLength = 1000),
                          subsamples = 1,
                          repeatsCV = 5,
                          foldsCV = 3,
                          positiveClass = "Respondedor"
)

save(Cytoxan_MScores, Cytoxan_model, Adryamicin_MScores, Adryamicin_model, pathways_annotated,
     file = "pathMED_results.RData")






###################################################################
##################   RESULTADOS #####################
###################################################################

################## AÑADIR COLUMNA DE ESTADIO DEL TUMOR A LOS DATAFRAMES DE CADA MEDICAMENTO ###################

################# PARA UTILIZARLOS EN LA CREACION DEL HEATMAP #################


################# CYTOXAN #############

Resultados_cytoxan <- inner_join(Cytoxan_pacientes, clinical_patient_brca, by = "bcr_patient_barcode")

row.names(Resultados_cytoxan) <- row.names(Cytoxan_pacientes)
# #IDs_repetidos <- Resultados_cytoxan %>%
#   group_by(bcr_patient_barcode) %>%
#   filter(n() > 1) %>%
#   pull(bcr_patient_barcode)

#Resultados_cytoxan <- Resultados_cytoxan %>% distinct(bcr_patient_barcode, .keep_all = TRUE)

# Esto asume que 'ID' es el nombre de tu columna de identificación y 'nueva_columna' es la columna que quieres transferir

Cytoxan_pacientes$ajcc_pathologic_tumor_stage <- Resultados_cytoxan$ajcc_pathologic_tumor_stage

Cytoxan_nuevo <- Cytoxan_pacientes


################# ADYAMICIN ##############

Resultados_adryamicin<- inner_join(Adryamicin_pacientes, clinical_patient_brca, by = "bcr_patient_barcode")
row.names(Resultados_adryamicin) <- row.names(Adryamicin_pacientes)

# Esto asume que 'ID' es el nombre de tu columna de identificación y 'nueva_columna' es la columna que quieres transferir
Adryamicin_pacientes$ajcc_pathologic_tumor_stage <- Resultados_adryamicin$ajcc_pathologic_tumor_stage


Adryamicin_nuevo <- Adryamicin_pacientes


save(Adryamicin_nuevo, Cytoxan_nuevo, file = "C:/Users/huart/OneDrive/Escritorio/TFM/clinical_cytoxan_adriamicin.RData")




############# HEATMAPS #################

# Cargar la biblioteca necesaria

install.packages("pheatmap")
install.packages("RColorBrewer")

library(pheatmap)
library(RColorBrewer)
library(dplyr)


# Para Cytoxan

# 1. Preparar los datos de anotación para las columnas basados en tus datos clínicos.


anotaciones_columnas <- data.frame(
  Respuesta = Cytoxan_nuevo$Grado_Respuesta,  # o el nombre correspondiente de la columna
  Estadio = Cytoxan_nuevo$ajcc_pathologic_tumor_stage # puedes agregar más variables clínicas aquí
)

# Convertir factores a caracteres 
anotaciones_columnas <- data.frame(lapply(anotaciones_columnas, as.character), stringsAsFactors = FALSE)



# 3. Crear el heatmap sin nombres de muestras y con las anotaciones de las columnas.
# También vamos a especificar un archivo para guardar el heatmap.

rownames(anotaciones_columnas) <- rownames(Cytoxan_nuevo)

pheatmap(
  Cytoxan_MScores,
  annotation_col = anotaciones_columnas,  
  show_colnames = FALSE, 
  fontsize_row = 13,  
  filename = "C:/Users/huart/OneDrive/Escritorio/TFM/Heatmap_cytoxan1.png",  
  width = 14,  # Ajusta el ancho del heatmap
  height = 14,
  main = "Mapa de calor de M-scores para Cytoxan", 
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  show_rownames = TRUE, 
  
)

# Para Adryamicin

# 1. Preparar los datos de anotación para las columnas basados en tus datos clínicos.

anotaciones_columnas <- data.frame(
  Respuesta = Adryamicin_nuevo$Grado_Respuesta,  # o el nombre correspondiente de la columna
  Estadio = Adryamicin_nuevo$ajcc_pathologic_tumor_stage # puedes agregar más variables clínicas aquí
)

# Convertir factores a caracteres 
anotaciones_columnas <- data.frame(lapply(anotaciones_columnas, as.character), stringsAsFactors = FALSE)



# 3. Crear el heatmap sin nombres de muestras y con las anotaciones de las columnas.
# También vamos a especificar un archivo para guardar el heatmap.

rownames(anotaciones_columnas) <- rownames(Adryamicin_nuevo)

pheatmap(
  Adryamicin_MScores,
  annotation_col = anotaciones_columnas,  
  show_colnames = FALSE,  # Esto quitará los nombres de las columnas (nombres de las muestras)
  fontsize_row = 13,  
  filename = "C:/Users/huart/OneDrive/Escritorio/TFM/Heatmap_adryamicin.png",  # Esto guardará tu heatmap como un archivo PDF
  width = 14,  # Ajusta el ancho del heatmap
  height = 14,
  main = "Mapa de calor de M-scores para Adryamicin", 
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  show_rownames = TRUE
)


############# GGPLOTS #################

# Cargar la biblioteca necesaria

install.packages("gplots")
install.packages("ggplot2")
library(gplots)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)

########## CYTOXAN #########

write.table(Cytoxan_model$stats)

# Recreando tu tabla de datos.
data <- data.frame(
  metrica = c("MCC", "AUC", "Accuracy", "Sensitivity", "Specificity", 
              "Pos_Pred_Value", "Neg_Pred_Value", "Precision", "Recall", 
              "F1", "Prevalence", "Detection_Rate", "Detection_Prevalence", 
              "Balanced_Accuracy"),
  nb = c(0.37, 0.86, 0.925, 0.947368421052632, 0.5, 0.972972972972973, 
         0.333333333333333, 0.972972972972973, 0.947368421052632, 0.96, 
         0.95, 0.9, 0.925, 0.723684210526316),
  lda = c(-0.115, 0.45, 0.75, 0.789473684210526, 0, 0.9375, 0, 0.9375, 
          0.789473684210526, 0.857142857142857, 0.95, 0.75, 0.8, 
          0.394736842105263),
  glm = c(-0.124, 0.24, 0.725, 0.763157894736842, 0, 0.935483870967742, 
          0, 0.935483870967742, 0.763157894736842, 0.840579710144928, 
          0.95, 0.725, 0.775, 0.381578947368421),
  rf = c(NA, 0.63, 0.95, 1, 0, 0.95, NA, 0.95, 1, 0.974358974358974, 
         0.95, 0.95, 1, 0.5),
  svmLinear = c(NA, 0.64, 0.95, 1, 0, 0.95, NA, 0.95, 1, 0.974358974358974, 
                0.95, 0.95, 1, 0.5)
)

# Transformamos los datos de formato ancho a largo
datos_largos <- data %>% 
  pivot_longer(
    cols = -metrica, 
    names_to = "modelo", 
    values_to = "valor"
  )

# Ahora puedes verificar tus datos transformados
print(datos_largos)

# Crear el diagrama de barras agrupado
ggplot(datos_largos, aes(fill=modelo, y=valor, x=metrica)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + # Ajustar el ancho de las barras
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotar y aumentar el tamaño de las etiquetas del eje x
        axis.text.y = element_text(size = 12),  # Aumentar el tamaño de las etiquetas del eje y
        axis.title.x = element_text(size = 14),  # Aumentar el tamaño del título del eje x
        axis.title.y = element_text(size = 14))  # Aumentar el tamaño del título del eje y


######## ADRYAMICIN ########

write.table(Adryamicin_model$stats)

# Crear un dataframe con tus datos

data_adr <- data.frame(
  metrica = c("MCC", "AUC", "Accuracy", "Sensitivity", "Specificity", 
              "Pos_Pred_Value", "Neg_Pred_Value", "Precision", "Recall", 
              "F1", "Prevalence", "Detection_Rate", "Detection_Prevalence", 
              "Balanced_Accuracy"),
  glm = c(0.267, 0.5, 0.538461538461538, 0.5, 1, 1, 0.142857142857143, 1, 0.5, 0.666666666666667, 
          0.923076923076923, 0.461538461538462, 0.461538461538462, 0.75),
  lda = c(-0.158, 0.5, 0.692307692307692, 0.75, 0, 0.9, 0, 0.9, 0.75, 0.818181818181818, 
          0.923076923076923, 0.692307692307692, 0.769230769230769, 0.375),
  nb = c(-0.192, 0.21, 0.615384615384615, 0.666666666666667, 0, 0.888888888888889, 0, 0.888888888888889, 
         0.666666666666667, 0.761904761904762, 0.923076923076923, 0.615384615384615, 0.692307692307692, 0.333333333333333),
  rf = c(NA, 0.08, 0.923076923076923, 1, 0, 0.923076923076923, NA, 0.923076923076923, 1, 0.96, 
         0.923076923076923, 0.923076923076923, 1, 0.5),
  svmLinear = c(NA, 0.17, 0.923076923076923, 1, 0, 0.923076923076923, NA, 0.923076923076923, 1, 0.96, 
                0.923076923076923, 0.923076923076923, 1, 0.5)
)


# Transformamos los datos de formato ancho a largo
datos_largos <- data_adr %>% 
  pivot_longer(
    cols = -metrica, 
    names_to = "modelo", 
    values_to = "valor"
  )

# Crear el diagrama de barras agrupado
ggplot(datos_largos, aes(fill=modelo, y=valor, x=metrica)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + # Ajustar el ancho de las barras
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotar y aumentar el tamaño de las etiquetas del eje x
        axis.text.y = element_text(size = 12),  # Aumentar el tamaño de las etiquetas del eje y
        axis.title.x = element_text(size = 14),  # Aumentar el tamaño del título del eje x
        axis.title.y = element_text(size = 14))  # Aumentar el tamaño del título del eje y







##### extraer tablas etc, para la memoria  #############

Pathways_20 <- pathways_annotated[1:20, ]

write.csv(Pathways_20, file = "Pathways.csv", row.names = FALSE)
install.packages("officer")
install.packages("flextable")
library(officer)
library(flextable)

my_flextable <- flextable::flextable(pathways_annotated)
my_flextable <- flextable::autofit(my_flextable)

# Create a new Word document
Pathways_20 <- officer::read_docx()

# Add the table to the document
Pathways_20 <- flextable::add_flextable(Pathways_20, value = my_flextable)
# Save the document to a file
output_filepath <- "C:/Users/huart/OneDrive/Escritorio/TFM/pathways_20.docx"  # specify the path and filename
officer::print(Pathways_20, target = output_filepath)


write.csv(pathways_annotated, file = "C:/Users/huart/OneDrive/Escritorio/TFM/Pathways_total.csv", row.names = FALSE)

write.csv(Cytoxan_nuevo, file = "C:/Users/huart/OneDrive/Escritorio/TFM/Cytoxan_clinico.csv", row.names = FALSE)
write.csv(Adryamicin_nuevo, file = "C:/Users/huart/OneDrive/Escritorio/TFM/Adryamicin_clinico.csv", row.names = FALSE)






