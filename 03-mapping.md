---
layout: page
title: Análisis transcriptómicos
subtitle: Alineamiento de lecturas
minutes: 5
---
> ## Objetivos de aprendizaje {.objectives}
>
> *  Alinear datos de secuenciación a una referencia (genomas o transcriptomas)
> *  Entender como interpretar datos de alineamiento de secuenciación masiva
> *  Entender los formatos de coordenadas bam y sam

> ## Datos requeridos {.prereq}
>
> Para proceder con esta práctica, se requieren los resultados de tarea de 
> ensamble de transcriptoma de la clase pasada:
>
> *  Archivos fastq de lecturas filtradas por calidad usando Trimmomatic
> *  Archivo fasta de ensamble de transcriptoma completo

Usaremos los archivos fastq que utilizamos la clase pasada, pero 
filtrados por calidad por Trimmomatic via Trinity 
resultantes de su tarea. Estos se encuentran el en directorio
`trinity_out_dir` con la extensión `P.qtrim.gz`. 

Así como el genoma de referencia de nuestro organismo, 
*Saccharomyces pombe*:

[Sp_genome.fa](datasets/genome/Sp_genome.fa) 

La lista de archivos debe ser similar a esta:

~~~ {.output}
~~~

Como ya sabemos, es buena idea primero ver que los datos están en el 
formato correcto:

~~~ {.bash}
$ for i in *qtrim.gz ; do zcat $i | head ; wait ; done 
~~~

Una vez que verificamos que las lecturas están el formato correcto,
vamos a mapear los datos al genoma utilizando Bowtie2 via TopHat. 

Pueden encontrar el manual en el siguiente [link](https://ccb.jhu.edu/software/tophat/manual.shtml).

## Mapeando los transcritos al genoma

En los casos en los que contamos con una referencia, es muy útil analizar los transcritos
en su contexto genómico. En este ejercicio usaremos el programa de mapeo 
[GMAP](http://research-pub.gene.com/gmap/).

El primer paso para realizar el mapeo es generar el índice por medio del comando:

~~~ {.bash}
$ gmap_build -d genome -D . -k 13 Sp_genome.fa
~~~
~~~ {.output}
~~~

Una vez generado el índice mapeamos los transcritos generados en la tarea al genoma:

~~~ {.bash}
$ gmap -n 0 -D . -d genome trinity_out_dir/Trinity.fasta -f samse > trinity_gmap.sam
~~~
~~~ {.output}
~~~

Si revisamos el resultado, es un archivo con coordenadas en un formato llamado
SAM (Sequence Alignment/Map format). Veamos en que consiste este formato:

~~~ {.bash}
$ head -n trinity_gmap.sam
~~~
~~~ {.output}
~~~


### El formato SAM

Supongamos que tenemos el siguiente alineamiento:

~~~ {.output}
Coor	12345678901234 5678901234567890123456789012345
ref	AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT
+r001/1	      TTAGATAAAGGATA*CTG
+r002	     aaaAGATAA*GGATA
+r003	   gcctaAGCTAA
+r004	                 ATAGCT..............TCAGC
-r003	                        ttagctTAGGC
-r001/2	                                      CAGCGGCAT
~~~

El formato SAM correspondiente será el siguiente:

~~~ {.output}
@HD	VN:1.5	SO:coordinate
@SQ	SN:ref	LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r002	0	ref	9	30	3S6M1P1I4M	*	0	0	AAAAGATAAGGATA	*
r003	0	ref	9	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:ref,29,-,6H5M,17,0;
r004	0	ref	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
r003	2064	ref	29	17	6H5M	*	0	0	TAGGC	*	SA:Z:ref,9,+,5S6M,30,1;
r001	147	ref	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
~~~

El formato SAM es un formato de texto plano que nos permite guardar datos de 
secuenciación en formato ASCII delimitado por tabulaciones. 

Está compuesto de dos secciones principales:

*  El encabezado
*  El alineamiento

La sección del **encabezado** comienza con el caracter `@` seguido por uno de lo códigos 
de dos letras que sirven para características de los alineamientos en este archivo. 
Cada línea esta delimitada por tabulaciones y, además de las líneas que comienzan con 
`@CO`, cada campo de datos tiene el formato `TAG:VALUE`, en donde `TAG` es una cadena
de dos caracteres que define el formato y contenido de `VALUE`. 

El encabezado no es indispensable pero contiene información acerca de la versión del 
archivo así como si está ordenado o no. Por ello es recomendable incluirlo.

La sección de *alineamiento* contiene la siguiente información:

1. **QNAME** Nombre de la referencia, QNAME (SAM)/Nombre de la lectura(BAM). 
Se utiliza para agrupar alineamientos que están juntos, como es el caso de alineamientos
de lecturas por pares o una lectura que aparece en alineamientos múltiples. 
2.  **FLAG** Set de información describiendo el alineamiento. Provee la siguiente información:
..* ¿Hay múltiples fragmentos?are there multiple fragments?
..* ¿Todos los fragmentos están bien alineados?are all fragments properly aligned?
..* ¿Está alineado este fragmento?is this fragment unmapped?
..* ¿No ha sido alineado el siguiente fragmento?is the next fragment unmapped?
..* ¿Es esta referencia la cadena inversa?is this query the reverse strand?
..* ¿El el siguiente fragmento la cadena reversa?is the next fragment the reverse strand?
..* ¿Es este el primer fragmento?is this the 1st fragment?
..* ¿Es este el último fragmento?is this the last fragment?
..* ¿Es este un alineamiento secundario?is this a secondary alignment?
..* ¿Esta lectura falló los filtros de calidad?did this read fail quality controls?
..* ¿Es esta lectura un duplicado por PCR o óptico?is this read a PCR or optical duplicate?
3. **RNAME** Nombre de la secuencia de referencia.
3. **POS** Posición de alineamiento izquierda (base 1).
3. **MAPQ** Calidad del alineamiento.
3. **CIGAR** cadena CIGAR.
3. **RNEXT** Nombre de referencia del par (mate) o la siguiente lectura.
3. **PNEXT** Posición del par (mate) o la siguiente lectura.
3. **TLEN** Longitud del alineamiento.
3. **SEQ** La secuencia de prueba de este alineamiento (en este caso la secuencia de la lectura).
3. **QUAL** La calidad de la lectura.
3. **TAGs** Información adicional. 

> ## Cadenas CIGAR {.callout}
>
> La secuencia alineada a la referencia puede tener bases adicionales que no están en 
> la referencia o puede no tener bases en la lectura que si están en la referencial.
> La cadena CIGAR es una cadena que codifica cada base y la caracteristica de cada una en 
> el alineamiento.
> 
> Por ejemplo, la cadena CIGAR:
> 
> ~~~ {.output}
> CIGAR: 3M1I3M1D5M
> ~~~
> 
> indica que las primera 3 bases de la lectura alinea con la referencia (3M), la siguiente base
> no existe en la referencia (1I), las siguientes 3 bases alinean con la referencia (3M), la 
> siguiente base no existe en la lectura (1D), y 5 bases más alinean con la referencia (5M). 

Como pueden ver estos archivos contienen muchísima información

La versión comprimida de los archivos tipo SAM se conoce como BAM (binary sam). 
Convirtamos el archivo SAM a BAM usando samtools:

~~~ {.bash}
$ samtools view -Sb trinity_gmap.sam > trinity_gmap.bam
~~~

No abriremos este archivo ya que, dado que esta en formato binario, es ilegible. 
Sin embargo, tenemos que realizar dos últimos pasos visualizar 
los resultados:

~~~ {.bash}
$ samtools sort trinity_gmap.bam trinity_gmap
$ samtools index trinity_gmap.bam
~~~

El primer paso ordena los resultados por sus coordenadas y el segundo crea índices
para hacer más rápida la visualización usando un navegador.


## Mapeando las lecturas filtradas al genoma

Primero generaremos un índice de bowtie2 para el genoma:

~~~ {.bash}
$ bowtie2-build Sp_genome.fa Sp_genome 
~~~ 
~~~ {.output}
~~~

Usamos tophat2 para mapear las lecturas. Este programa nos permite dividir lecturas
que atraviesan sitios de splicing 

~~~ {.bash}
$ tophat2 -I 300 -i 20 genome \
    RNASEQ_data/Sp_log.left.fq.gz,RNASEQ_data/Sp_hs.left.fq.gz,RNASEQ_data/Sp_ds.left.fq.gz,RNASEQ_data/Sp_plat.left.fq.gz \
    RNASEQ_data/Sp_log.right.fq.gz,RNASEQ_data/Sp_hs.right.fq.gz,RNASEQ_data/Sp_ds.right.fq.gz,RNASEQ_data/Sp_plat.right.fq.gz
~~~ 

Exploramos el resultado, el cuál es un archivo tipo sam. 

~~~ {.bash}
$ head X 
~~~ 

> ## Tarea - Alineando las lecturas filtradas al transcriptoma {.challenge}
>
> Hemos alineado las lecturas al genoma pero queremos alinearlas también directamente
> al transcriptoma. La tarea consiste en revisar el manual de TopHat2 y usar las opciones
> que nos permite mapear lecturas directamente a transcriptomas. 
> 
> Agreguen el archivo Trinity.fasta de referencia y el archivo de mapeo en formato BAM
> (ordenado e indizado) a su repositorio en un nuevo directorio llamado Transcriptome_Mapping.
> Lo deberán agregar antes de la clase del viernes 8 de abril.
> 
> *Pista 1:* No podrán utilizar el índice generado previamente.
> *Pista 2:* No deberán iniciar un nuevo repositorio, solo tienen que agregar el nuevo 
> directorio y su contenido. 







