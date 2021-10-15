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
> *  Primer acercamiento a los formatos de coordenadas SAM y BAM

> ## Datos requeridos {.prereq}
>
> Para proceder con esta práctica, se requieren los resultados de tarea de 
> ensamble de transcriptoma de la clase pasada:
>
> *  Archivos fastq de lecturas filtradas por calidad usando Trimmomatic via Trinity
> *  Archivo fasta de ensamble de transcriptoma completo

Usaremos los archivos fastq que utilizamos la clase pasada.

También usaremos el genoma de referencia de nuestro organismo, 
*Saccharomyces pombe*:

~~~ {.output}
$ wget https://liz-fernandez.github.io/PBI_transcriptomics/datasets/genome/Sp_genome.fa
~~~

Los datos ya se encuentran en su directorio de datos en Docker. 

Vamos a alinear las lecturas y transcritos al genoma utilizando HISAT2 

Pueden encontrar el manual en el siguiente [link](http://daehwankimlab.github.io/hisat2/manual/).

Creemos un directorio para guardar nuestros análisis:

~~~ {.bash}
$ cd /usr/local/
$ mkdir ANALYSIS
$ cd ANALYSIS
$ mkdir MAPPING
$ cd MAPPING
~~~

Copiemos los archivos de las lecturas y el genoma:

~~~ {.bash}
$ cp /usr/local/data/Sp*fq.gz .
$ cp /usr/local/data/Sp_genome.fa.gz .
~~~

Y descomprimimos el genoma:

~~~ {.bash}
$ gunzip Sp_genome.fa.gz
~~~

## Alineando las lecturas filtradas al genoma

Primero generaremos un índice de hisat2 para el genoma:

~~~ {.bash}
$ hisat2-build Sp_genome.fa Sp_genome 
~~~

~~~ {.output}
Settings:
  Output files: "Sp_genome.*.ht2"
  Line rate: 6 (line is 64 bytes)
  Lines per side: 1 (side is 64 bytes)
  Offset rate: 4 (one in 16)
  FTable chars: 10
  Strings: unpacked
  Local offset rate: 3 (one in 8)
  Local fTable chars: 6
  Local sequence length: 57344
  Local sequence overlap between two consecutive indexes: 1024
  Endianness: little
  Actual local endianness: little
  Sanity checking: disabled
  Assertions: disabled
  Random seed: 0
  Sizeofs: void*:8, int:4, long:8, size_t:8
Input files DNA, FASTA:
  Sp_genome.fa
Reading reference sizes
  Time reading reference sizes: 00:00:00
Calculating joined length
Writing header
Reserving space for joined string
Joining reference sequences
  Time to join reference sequences: 00:00:00
  Time to read SNPs and splice sites: 00:00:00
Using parameters --bmax 80836 --dcv 1024
  Doing ahead-of-time memory usage test
  Passed!  Constructing with these parameters: --bmax 80836 --dcv 1024
Constructing suffix-array element generator
Building DifferenceCoverSample
  Building sPrime
  Building sPrimeOrder
  V-Sorting samples
  V-Sorting samples time: 00:00:00
  Allocating rank array
  Ranking v-sort output
  Ranking v-sort output time: 00:00:00
  Invoking Larsson-Sadakane on ranks
  Invoking Larsson-Sadakane on ranks time: 00:00:00
  Sanity-checking and returning
Building samples
Reserving space for 12 sample suffixes
Generating random suffixes
QSorting 12 sample offsets, eliminating duplicates
QSorting sample offsets, eliminating duplicates time: 00:00:00
Multikey QSorting 12 samples
  (Using difference cover)
  Multikey QSorting samples time: 00:00:00
Calculating bucket sizes
Splitting and merging
  Splitting and merging time: 00:00:00
Split 1, merged 7; iterating...
Splitting and merging
  Splitting and merging time: 00:00:00
Split 1, merged 1; iterating...
Splitting and merging
  Splitting and merging time: 00:00:00
Avg bucket size: 71853.5 (target: 80835)
Converting suffix-array elements to index image
Allocating ftab, absorbFtab
Entering GFM loop
Getting block 1 of 6
  Reserving size (80836) for bucket 1
  Calculating Z arrays for bucket 1
  Entering block accumulator loop for bucket 1:
  bucket 1: 10%
  bucket 1: 20%
  bucket 1: 30%
  bucket 1: 40%
  bucket 1: 50%
  bucket 1: 60%
  bucket 1: 70%
  bucket 1: 80%
  bucket 1: 90%
  bucket 1: 100%
  Sorting block of length 79345 for bucket 1
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 79346 for bucket 1
Getting block 2 of 6
  Reserving size (80836) for bucket 2
  Calculating Z arrays for bucket 2
  Entering block accumulator loop for bucket 2:
  bucket 2: 10%
  bucket 2: 20%
  bucket 2: 30%
  bucket 2: 40%
  bucket 2: 50%
  bucket 2: 60%
  bucket 2: 70%
  bucket 2: 80%
  bucket 2: 90%
  bucket 2: 100%
  Sorting block of length 78903 for bucket 2
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 78904 for bucket 2
Getting block 3 of 6
  Reserving size (80836) for bucket 3
  Calculating Z arrays for bucket 3
  Entering block accumulator loop for bucket 3:
  bucket 3: 10%
  bucket 3: 20%
  bucket 3: 30%
  bucket 3: 40%
  bucket 3: 50%
  bucket 3: 60%
  bucket 3: 70%
  bucket 3: 80%
  bucket 3: 90%
  bucket 3: 100%
  Sorting block of length 73828 for bucket 3
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 73829 for bucket 3
Getting block 4 of 6
  Reserving size (80836) for bucket 4
  Calculating Z arrays for bucket 4
  Entering block accumulator loop for bucket 4:
  bucket 4: 10%
  bucket 4: 20%
  bucket 4: 30%
  bucket 4: 40%
  bucket 4: 50%
  bucket 4: 60%
  bucket 4: 70%
  bucket 4: 80%
  bucket 4: 90%
  bucket 4: 100%
  Sorting block of length 70864 for bucket 4
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 70865 for bucket 4
Getting block 5 of 6
  Reserving size (80836) for bucket 5
  Calculating Z arrays for bucket 5
  Entering block accumulator loop for bucket 5:
  bucket 5: 10%
  bucket 5: 20%
  bucket 5: 30%
  bucket 5: 40%
  bucket 5: 50%
  bucket 5: 60%
  bucket 5: 70%
  bucket 5: 80%
  bucket 5: 90%
  bucket 5: 100%
  Sorting block of length 74487 for bucket 5
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 74488 for bucket 5
Getting block 6 of 6
  Reserving size (80836) for bucket 6
  Calculating Z arrays for bucket 6
  Entering block accumulator loop for bucket 6:
  bucket 6: 10%
  bucket 6: 20%
  bucket 6: 30%
  bucket 6: 40%
  bucket 6: 50%
  bucket 6: 60%
  bucket 6: 70%
  bucket 6: 80%
  bucket 6: 90%
  bucket 6: 100%
  Sorting block of length 53694 for bucket 6
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 53695 for bucket 6
Exited GFM loop
fchr[A]: 0
fchr[C]: 135086
fchr[G]: 215577
fchr[T]: 296650
fchr[$]: 431126
Exiting GFM::buildToDisk()
Returning from initFromVector
Wrote 4340624 bytes to primary GFM file: Sp_genome.1.ht2
Wrote 107788 bytes to secondary GFM file: Sp_genome.2.ht2
Re-opening _in1 and _in2 as input streams
Returning from GFM constructor
Returning from initFromVector
Wrote 201185 bytes to primary GFM file: Sp_genome.5.ht2
Wrote 109640 bytes to secondary GFM file: Sp_genome.6.ht2
Re-opening _in5 and _in5 as input streams
Returning from HierEbwt constructor
Headers:
    len: 431126
    gbwtLen: 431127
    nodes: 431127
    sz: 107782
    gbwtSz: 107782
    lineRate: 6
    offRate: 4
    offMask: 0xfffffff0
    ftabChars: 10
    eftabLen: 0
    eftabSz: 0
    ftabLen: 1048577
    ftabSz: 4194308
    offsLen: 26946
    offsSz: 107784
    lineSz: 64
    sideSz: 64
    sideGbwtSz: 48
    sideGbwtLen: 192
    numSides: 2246
    numLines: 2246
    gbwtTotLen: 143744
    gbwtTotSz: 143744
    reverse: 0
    linearFM: Yes
Total time for call to driver() for forward index: 00:00:01
~~~

El índice se guarda en archivos con extension .ht2:
 
~~~ {.bash}
$ ls *ht2
~~~ 

~~~ {.output}
Sp_genome.1.ht2
Sp_genome.2.ht2
Sp_genome.3.ht2
Sp_genome.4.ht2
Sp_genome.5.ht2
Sp_genome.6.ht2
Sp_genome.7.ht2
Sp_genome.8.ht2
~~~

Usamos hisat2 para alinear las lecturas. Este programa nos permite dividir lecturas
que atraviesan sitios de splicing:

~~~ {.bash}
$ hisat2 --max-intronlen 300 --min-intronlen 20 -x ./Sp_genome -1 Sp_log.left.fq.gz  -2 Sp_log.right.fq.gz -S Sp_log.sam
~~~ 

~~~ {.output}
37915 reads; of these:
  37915 (100.00%) were paired; of these:
    2945 (7.77%) aligned concordantly 0 times
    34969 (92.23%) aligned concordantly exactly 1 time
    1 (0.00%) aligned concordantly >1 times
    ----
    2945 pairs aligned concordantly 0 times; of these:
      278 (9.44%) aligned discordantly 1 time
    ----
    2667 pairs aligned 0 times concordantly or discordantly; of these:
      5334 mates make up the pairs; of these:
        3354 (62.88%) aligned 0 times
        1980 (37.12%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
95.58% overall alignment rate
~~~

Exploramos el resultado, el cuál es un archivo tipo SAM. 

~~~ {.bash}
$ head Sp_log.sam
~~~ 

### El formato SAM

Veamos un ejemplo más pequeño de este formato.
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

La sección de **alineamiento** contiene la siguiente información:

1. **QNAME**	Nombre de la referencia, QNAME (SAM)/Nombre de la lectura(BAM). 
Se utiliza para agrupar alineamientos que están juntos, como es el caso de alineamientos
de lecturas por pares o una lectura que aparece en alineamientos múltiples. 
2.  **FLAG**	Set de información describiendo el alineamiento. Provee la siguiente información:
	* ¿Hay múltiples fragmentos?
	* ¿Todos los fragmentos están bien alineados?
	* ¿Está alineado este fragmento?
	* ¿No ha sido alineado el siguiente fragmento?
	* ¿Es esta referencia la cadena inversa?
	* ¿El el siguiente fragmento la cadena reversa?
	* ¿Es este el primer fragmento?
	* ¿Es este el último fragmento?
	* ¿Es este un alineamiento secundario?
	* ¿Esta lectura falló los filtros de calidad?
	* ¿Es esta lectura un duplicado por PCR o óptico?
3. **RNAME**	Nombre de la secuencia de referencia.
3. **POS**	Posición de alineamiento izquierda (base 1).
3. **MAPQ**	Calidad del alineamiento.
3. **CIGAR**	cadena CIGAR.
3. **RNEXT**	Nombre de referencia del par (mate) o la siguiente lectura.
3. **PNEXT**	Posición del par (mate) o la siguiente lectura.
3. **TLEN**	Longitud del alineamiento.
3. **SEQ**	La secuencia de prueba de este alineamiento (en este caso la secuencia de la lectura).
3. **QUAL**	La calidad de la lectura.
3. **TAGs**	Información adicional. 

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

Como pueden ver estos archivos contienen muchísima información que puede ser analizada
usando scripts que arrojen estadísticas del alineamiento. Programas como 
[Picard](http://broadinstitute.github.io/picard/) realizan
este tipo de análisis.

La versión comprimida de los archivos tipo SAM se conoce como BAM (binary sam). 
Convirtamos el archivo SAM a BAM usando samtools:

~~~ {.bash}
$ samtools view -S -b Sp_log.sam > Sp_log.bam
~~~

No abrimos los archivos BAM ya que, dado que esta en formato binario, son ilegible.
Para poder leerlos necesitamos usar samtools con el comando:

~~~ {.bash}
$ samtools view Sp_log.bam | head
~~~ 

> ## Tarea 1 - Alineando las muestras restantes al genoma {.challenge}
>
> Alineamos una muestra al genoma pero queremos alinear todas las muestras.  
> Alinea el resto de las muestras y guarda los resultados (en formato bam), ya que los necesitaras en 
> prácticas subsecuentes. 
>

> ## Tarea 2 (Opcional) - Alineando las lecturas filtradas al transcriptoma {.challenge}
>
> Hemos alineado las lecturas al genoma pero queremos alinearlas también directamente
> al transcriptoma. Revisa el manual de HISAT2 y usar las opciones
> que nos permite mapear lecturas directamente a transcriptomas. Deberás alinear cada 
> muestra de manera independiente.  Alinea el resto de las muestras y guarda los 
> resultados (en formato bam), ya que los necesitaras en prácticas subsecuentes.
> 
> * **Pista:** No podrán utilizar el índice generado previamente.
>










