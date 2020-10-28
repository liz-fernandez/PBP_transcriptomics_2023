---
layout: page
title: Análisis transcriptómicos
subtitle: Análisis de expresión diferencial
minutes: 5
---
> ## Objetivos de aprendizaje {.objectives}
>
> *  Aprender a realizar análisis de expresión diferencial usando DESeq2

> ## Datos requeridos {.prereq}
>
> Para proceder con esta práctica, se requieren:
>
> *  El genoma `Sp_genome.fa`
> *  Las lecturas no filtradas 

Utilizaremos los archivos bam filtrados por calidad que utilizamos en la clase 
de mapeo. Los análisis se realizaran en la plataforma R.

Generemos un directorio de trabajo en su directorio de preferencia

~~~ {.bash}
$ mkdir DESEQ2
$ cd DESEQ2/
~~~

## Contando lecturas

El primer paso es contar lecturas. Tanto DESeq2 como otros programas que usan
la distribución negativa binomial utilizan conteos **crudos** (raw). Esto 
se refiere la número de lecturas que alinean a un transcrito específico sin 
sufrir ningún tipo de transformación. Como vieron en la clase pasada, esto se 
debe a que el programa estima un factor de normalización. Si usamos datos 
previamente normalizados estaremos violando los supuestos. 

Descargaremos la anotación de los de *Schizosaccharomyces pombe* genes de aquí:

~~~ {.bash}
$ wget https://liz-fernandez.github.io/PBI_transcriptomics_2020/DATA/Sp_genes.gtf
~~~

Usaremos el programa [HTSeq](https://pypi.python.org/packages/source/H/HTSeq/). Ya se encuentra 
instalado en Docker.

Este programa cuenta el número de lecturas que alinean con cada gen anotado en un archivo 
de referencia. Nosotros le proporcionaremos nuestro archivo `Sp_genes.gtf`.

~~~ {.bash}
$ htseq-count -f bam -s no -r pos Sp_ds_sorted.bam /usr/local/data/Sp_genes.gtf > Sp_ds.counts 
$ htseq-count -f bam -s no -r pos Sp_hs_sorted.bam /usr/local/data/Sp_genes.gtf > Sp_hs.counts 
$ htseq-count -f bam -s no -r pos Sp_log_sorted.bam /usr/local/data/Sp_genes.gtf > Sp_log.counts 
$ htseq-count -f bam -s no -r pos Sp_plat_sorted.bam /usr/local/data/Sp_genes.gtf > Sp_plat.counts 
~~~

Veamos uno de los archivos de cuentas, primero las primeras 10 líneas:

~~~ {.bash}
$ head Sp_ds.count
~~~

~~~ {.output}
SPAC1002.19	23
SPAC1093.06c	48
SPAC10F6.01c	362
SPAC10F6.05c	91
SPAC11D3.18c	25
SPAC11E3.06	28
SPAC11E3.14	399
SPAC11H11.04	16
SPAC1250.07	71
SPAC12B10.14c	70
~~~

Y las últimas 10 líneas:

~~~ {.bash}
$ tail Sp_ds.count
~~~

~~~ {.output}
SPCC736.05	18
SPCC736.12c	57
SPCC757.03c	1064
SPCC794.10	68
SPCP25A2.02c	49
__no_feature	0
__ambiguous	2
__too_low_aQual	5483
__not_aligned	1579
__alignment_not_unique	3
~~~

Estas líneas al final del ejemplo muestran lecturas que:

1. **No feature**	Lecturas que no alinean a una región anotada (gen)
1. **Ambiguous**	Lecturas que alinean a más de un gen 
1. **Too low alignment quality**	Lecturas con baja calidad de alineamiento
1. **Not aligned**	Lecturas no alineadas
1. **Alignment not unique**	Lecturas que alinean más de una vez a la referencia


> ## Contando librerías de manera independiente {.challenge}
>
> ¿Por qué contamos las lecturas en cada librería de manera independiente si
> ensamblamos el transcriptoma de referencia usando todas las librerías (muestras)? 
> 

Es importante verificar que el conteo de HTSeq engloba a la mayoría de las lecturas. 
Si estos porcentajes son muy bajos podrían indicar problemas con los parámetros de 
conteo (e.g. si le decimos que nuestra muestra no tiene direccionalidad (strandness)
cuando si la tiene).

Para realizar nuestros análisis de expresión debemos eliminar estas líneas. Uniremos los resultados usando el siguiente comando:

~~~ {.bash}
$ paste *counts | cut -f 1,2,4,6,8 | grep '__' -v > Sp_counts_table.txt
$ head Sp_counts_table.txt
~~~

Esto unió los distintos archivos en uno solo el cual ya podemos utilizar en R. 
Es importante verificar el orden de pegado de los archivos para no confundir una 
muestra con otra. 

~~~ {.output}
SPAC1002.19	23	44	2	230
SPAC1093.06c	48	22	34	59
SPAC10F6.01c	362	323	527	200
SPAC10F6.05c	91	149	54	203
SPAC11D3.18c	25	11	12	23
SPAC11E3.06	28	8	5	4
SPAC11E3.14	399	148	53	387
SPAC11H11.04	16	6	12	23
SPAC1250.07	71	36	29	72
SPAC12B10.14c	70	37	52	47
~~~

### Análisis de expresión diferencial

Ahora usaremos R. Puedes hacer esto directamente en Docker o descargar la tabla de cuentas a tu computadora. Cargamos el paquete `DESeq2`:

~~~ {.r}
> library(“DESeq2”)
~~~

Este paquete estima una distribución binomial negativa.

Leemos las cuentas en R usando read.delim:

~~~ {.r}
> countData <- read.delim("./Sp_counts_table.txt",header=FALSE, row.names=1)
> head(countData)
~~~

~~~ {.output}
              V2  V3  V4  V5
SPAC1002.19   23  44   2 230
SPAC1093.06c  48  22  34  59
SPAC10F6.01c 362 323 527 200
SPAC10F6.05c  91 149  54 203
SPAC11D3.18c  25  11  12  23
SPAC11E3.06   28   8   5   4
~~~ 

Y añadimos la descripción de cada una de las muestras. 
Este es el orden en el cual las "pegamos" usando el comando `paste`. 
Como ya mencionamos, en sus experimentos deberán cerciorarse de que este orden es 
el correcto.

~~~ {.r}
> condition <- factor(c("ds","hs","log","plat"))
> colData <- data.frame(row.names=colnames(countData), condition)
> head(colData)
~~~

~~~ {.output}
   condition
V2        ds
V3        hs
V4       log
V5      plat
~~~ 

Una vez ensamblados nuestros datos podemos usar la función DESeqDataSetFromMatrix para 
ponerlo en el formato requerido por DESeq2:

~~~ {.r}
> dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition )
> dds
~~~

~~~ {.output}
class: DESeqDataSet 
dim: 200 4 
exptData(0):
assays(1): counts
rownames(200): SPAC1002.19 SPAC1093.06c ... SPCC794.10 SPCP25A2.02c
rowData metadata column names(0):
colnames(4): V2 V3 V4 V5
colData names(1): condition
~~~

Filtramos genes con muy muy baja expresión, en este caso todos aquellos que no tengan al
menos una lectura. Nota: DESeq2 aplica otros filtros más estrictos por default.

~~~ {.r}
> dds <- dds[ rowSums(counts(dds)) > 1, ]
> dds
~~~

~~~ {.output}
class: DESeqDataSet 
dim: 198 4 
exptData(0):
assays(1): counts
rownames(198): SPAC1002.19 SPAC1093.06c ... SPCC794.10 SPCP25A2.02c
rowData metadata column names(0):
colnames(4): V2 V3 V4 V5
colData names(1): condition
~~~ 

La funcion DESeq realiza el análisis estándar de expresión diferencial, 
incluyendo los pasos de estimación de la distribución. 
Esta función genera una tabla con cambios de expresión en base log2, 
así como p-values y p-values ajustados por pruebas múltiples. 

~~~ {.r}
> dds <- DESeq(dds)
~~~

~~~ {.output}
estimating size factors
estimating dispersions
Error in checkForExperimentalReplicates(object, modelMatrix) :

  The design matrix has the same number of samples and coefficients to fit,
  so estimation of dispersion is not possible. Treating samples
  as replicates was deprecated in v1.20 and no longer supported since v1.22.
~~~

Este error nos indica que DESeq2 no nos permite realizar el análisis debido a la falta de réplicas; no tiene manera de saber cuanta dispersión en 
los niveles de expresión se debe a variación biológica o a diferencia entre los 
distintos tratamientos. Por lo tanto no podemos hacer expresión diferencial con estos datos. 

Realizaremos el análisis con un set de datos distinto que descargaremos así:

~~~ {.bash}
$ wget XXXXX GSE60450_Lactation-GenewiseCounts.txt
~~~
############ CHECK

Leemos las cuentas en R usando read.delim:

~~~ {.r}
> countData <- read.delim("./Sp_counts_table_complete.txt",header=FALSE, row.names=1)
> head(countData)
~~~

~~~ {.output}
                V2    V3    V4    V5    V6    V7    V8    V9   V10   V11   V12
SPAC1002.19     15    10     5     6     6     2     4     5     2     1     2
SPAC1093.06c  6257  7345  6361  5570  3585  3388  9602  9742  7702  6697  2781
SPAC10F6.01c     1     0     0     0     0     0     0     0     0     0     0
SPAC10F6.05c   380   260   360   346   362   361   356   399   414   331   275
SPAC11D3.18c    36    30   147   176   174   151     0     0     0     0     0
SPAC11E3.06  19491 18797 27019 24121 19250 15293 23490 21142 26245 25419 10869
              V13
SPAC1002.19     4
SPAC1093.06c 2802
SPAC10F6.01c    0
SPAC10F6.05c  260
SPAC11D3.18c    0
SPAC11E3.06  7484
~~~ 

c("Sp.DG","Sp.DH","Sp.DI","Sp.DJ","Sp.DK","Sp.DL","Sp.LA","Sp.LB","Sp.LC","Sp.LD","Sp.LE","Sp.LF")

Y añadimos la descripción de cada una de las muestras. 
Este es el orden en el cual las "pegamos" usando el comando `paste`. 
Como ya mencionamos, en sus experimentos deberán cerciorarse de que este orden es 
el correcto.

~~~ {.r}
> sample_info <- read.delim("./Sp_sample_info.txt")
#> condition <- factor(c("ds","hs","log","plat"))
> colData <- data.frame(row.names=colnames(countData), condition)
> head(colData)
~~~

~~~ {.output}
   condition
V2        ds
V3        hs
V4       log
V5      plat
~~~ 

Una vez ensamblados nuestros datos podemos usar la función DESeqDataSetFromMatrix para 
ponerlo en el formato requerido por DESeq2:

~~~ {.r}
> dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ CellType + Status )
> dds
~~~

~~~ {.output}
class: DESeqDataSet 
dim: 200 4 
exptData(0):
assays(1): counts
rownames(200): SPAC1002.19 SPAC1093.06c ... SPCC794.10 SPCP25A2.02c
rowData metadata column names(0):
colnames(4): V2 V3 V4 V5
colData names(1): condition
~~~

Filtramos genes con muy muy baja expresión, en este caso todos aquellos que no tengan al
menos una lectura. Nota: DESeq2 aplica otros filtros más estrictos por default.

~~~ {.r}
> dds <- dds[ rowSums(counts(dds)) > 1, ]
> dds
~~~

~~~ {.output}
class: DESeqDataSet 
dim: 198 4 
exptData(0):
assays(1): counts
rownames(198): SPAC1002.19 SPAC1093.06c ... SPCC794.10 SPCP25A2.02c
rowData metadata column names(0):
colnames(4): V2 V3 V4 V5
colData names(1): condition
~~~ 

La funcion DESeq realiza el análisis estándar de expresión diferencial, 
incluyendo los pasos de estimación de la distribución. 
Esta función genera una tabla con cambios de expresión en base log2, 
así como p-values y p-values ajustados por pruebas múltiples. 

~~~ {.r}
> dds <- DESeq(dds)
~~~

~~~ {.output}
estimating size factors
estimating dispersions
Error in checkForExperimentalReplicates(object, modelMatrix) :

  The design matrix has the same number of samples and coefficients to fit,
  so estimation of dispersion is not possible. Treating samples
  as replicates was deprecated in v1.20 and no longer supported since v1.22.
~~~

~~~ {.r}
dds
~~~

~~~ {.output}
class: DESeqDataSet 
dim: 198 4 
exptData(0):
assays(3): counts mu cooks
rownames(198): SPAC1002.19 SPAC1093.06c ... SPCC794.10 SPCP25A2.02c
rowData metadata column names(37): baseMean baseVar ... deviance maxCooks
colnames(4): V2 V3 V4 V5
colData names(2): condition sizeFactor
~~~

> ## Diseño multifactorial {.callout}
> 
> Experimentos con mas de un factor se pueden analizar usando DESeq2 al incluirse 
> las variables adicionales. Por ejemplo, si se realizaron los mismos experimentos 
> en distintos grupos o batches podemos corregir errores experimentales. 
> Cada factor añadido incrementa la sensibilidad a diferencias causadas por la condición. 
> 
> Por ejemplo, podríamos añadir los siguientes factores
>
> ~~~ {.output}
> colData(dds)
> ## DataFrame with 7 rows and 3 columns
> ## condition type sizeFactor
> ## <factor> <factor> <numeric>
> ## treated1fb treated single-read 1.512
> ## treated2fb treated paired-end 0.784
> ## treated3fb treated paired-end 0.896
> ## untreated1fb untreated single-read 1.050
> ## untreated2fb untreated single-read 1.659
> ## untreated3fb untreated paired-end 0.712
> ## untreated4fb untreated paired-end 0.784
> ~~~
>
> Y analizarlos así para comparar el tipo de librería con el tratamiento
> 
> ~~~ {.r}
> design(ddsMF) <- formula(~ type + condition)
> ddsMF <- DESeq(ddsMF)
> ~~~
>

Podemos extraer los resultados usando la función res:

~~~ {.r}
res <- results(dds)
res
~~~
~~~ {.output}
log2 fold change (MAP): condition plat vs ds 
Wald test p-value: condition plat vs ds 
DataFrame with 6 rows and 6 columns
              baseMean log2FoldChange     lfcSE       stat    pvalue      padj
             <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
SPAC1002.19   68.10299     1.37929086 1.0630096  1.2975338 0.1944476 0.9704666
SPAC1093.06c  38.84823     0.06599578 0.6115834  0.1079097 0.9140673 0.9704666
SPAC10F6.01c 368.22680    -0.80731941 0.9035895 -0.8934582 0.3716119 0.9704666
SPAC10F6.05c  84.74091     0.83338935 0.8311666  1.0026743 0.3160180 0.9704666
SPAC11D3.18c  16.79657    -0.27280407 0.7300479 -0.3736797 0.7086426 0.9704666
SPAC11E3.06   10.98535    -2.00717671 1.0185420 -1.9706371 0.0487654 0.9704666
~~~

Podemos obtener más información acerca de cada columna:

~~~ {.r}
mcols(res)$description
~~~
~~~ {.output}
[1] "mean of normalized counts for all samples"   
[2] "log2 fold change (MAP): condition plat vs ds"
[3] "standard error: condition plat vs ds"        
[4] "Wald statistic: condition plat vs ds"        
[5] "Wald test p-value: condition plat vs ds"     
[6] "BH adjusted p-values"
~~~ 

La función `summary` nos proporciona un resumen de los indicadores de nuestros
resultados.

~~~ {.r}
summary(res)
~~~
~~~ {.output}
out of 198 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 0, 0% 
LFC < 0 (down)   : 0, 0% 
outliers [1]     : 0, 0% 
low counts [2]   : 0, 0% 
(mean count < 0.7)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
~~~

¿Qué nos sugieren estos resultados? Si vemos una gráfica MA vemos que no hay 
ningún gen con expresión diferencial significativa:

~~~ {.r}
plotMA(res, main="DESeq2", ylim=c(-2,2))
~~~

![MA Plot](fig/MA_plot.jpeg)

Como pueden ver estos resultados solo muestran la comparación entre las condiciones
'plat' y 'ds'. Para extraer otros contrastes extraemos distintas parte de los 
resultados guardados en dds usando la función `results`.

Extraemos todas las combinaciones:

~~~ {.r}
res_ds_vs_hs <- results(dds,contrast=c("condition","ds","hs"))
res_ds_vs_log <- results(dds,contrast=c("condition","ds","log"))
res_ds_vs_plat <- results(dds,contrast=c("condition","ds","plat"))
res_hs_vs_log <- results(dds,contrast=c("condition","hs","log"))
res_hs_vs_plat <- results(dds,contrast=c("condition","hs","plat"))
res_log_vs_plat <- results(dds,contrast=c("condition","log","plat"))
~~~

Y los escribimos a archivos planos (Solo se muestra el ejemplo del primer archivo):

~~~ {.r}
res_ds_vs_hs <- res_ds_vs_hs[order(res_ds_vs_hs$padj),]  # Ordenando
# Datos sin filtrar
write.csv(as.data.frame(subset(res_ds_vs_hs)),file="Results_ds_vs_hs_unfiltered.csv", quote=FALSE)
# Datos filtrados
write.csv(as.data.frame(subset(res_ds_vs_hs, padj < 0.1)),file="Results_ds_vs_hs_padj0.05.csv", quote=FALSE)
~~~

### Revisando efectos de grupo (batch effects)

Un efecto de grupo es cuando hay una diferencia entre tratamientos no debido
a cambios biológicos si no a la manera en que fue procesada la muestra. 

Para verificar los datos usamos transformaciones y visualizaciones. 
Como ejemplo, transformaremos los datos de expresión usando la transformación 
rlog. Estas transformaciones permiten visualizar mejor los datos aunque, como en todos
los casos, es bueno verificar que nuestros datos cumplen con los supuestos requeridos.
En este caso se requiere que los datos no difieran mucho entre las distintas condiciones. 

####### CHECK INSTALL LIBS

~~~ {.r}
rld <- rlogTransformation(dds, blind=FALSE) 
sampleDistsRLD <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrixRLD <- as.matrix(sampleDistsRLD)
rownames(sampleDistMatrixRLD) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrixRLD) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrixRLD,
         clustering_distance_rows=sampleDistsRLD,
         clustering_distance_cols=sampleDistsRLD,
         col=colors)
~~~

![Samples Heatmap](fig/SampleHeatmap.jpeg)

Ahora imprimimos un análisis de componentes principales (PCA). 

~~~ {.r}
pdf("Sp_Samples_PCA_rlogTransformed.pdf",10,10)
plotPCA(rld)
dev.off()
~~~

![Samples PCA](fig/Sp_Samples_PCA_rlogTransformed.jpg)

### Visualizando grupos de genes

Una de las maneras más comunes de visualizar genes es usando un heatmap. Un heatmap
es un gráfico con recuadros que representan la expresión de cada gen. Es importante
mencionar que los heatmaps son útiles pero pueden sesgar nuestros resultados, por ello
es importante seleccionar los genes a visualizar por medio de bases estadísticas y 
estar concientes de los parámetros de visualización que cambian la escala de los datos. 

Visualicemos los 20 genes más expresados en nuestro datos usando el paquete "pheatmap". 

~~~ {.r}
install.packages("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds))
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
+          cluster_cols=FALSE, annotation_col=df)
~~~ 

![Heatmap](fig/Heatmap.jpeg)

Existen muchas otras opciones de análisis y visualizaciones, las cuales pueden
consultar en la [Viñeta de DESeq2](https://www.bioconductor.org/packages/3.3/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf).


