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

Usaremos los archivos fastq que utilizamos la clase pasada, pero 
filtrados por calidad por Trimmomatic via Trinity 
resultantes de su tarea. Estos se encuentran el en directorio
`trinity_out_dir` con la extensión `*.qtrim.gz`. 

Con lecturas en pares, habrá casos en los que ambas lecturas sean filtradas (P), 
y casos en los que una lectura se conservó y la otra se desechó (U). Las lecturas 
filtradas se encuentran en los siguientes archivos (ejemplo):

~~~ {.output}
reads.left.fq.gz.P.qtrim.gz     
reads.left.fq.gz.U.qtrim.gz
reads.right.fq.gz.P.qtrim.gz
reads.right.fq.gz.U.qtrim.gz
~~~

También usaremos el genoma de referencia de nuestro organismo, 
*Saccharomyces pombe*:

~~~ {.output}
$ wget https://liz-fernandez.github.io/PBI_transcriptomics/datasets/genome/Sp_genome.fa
~~~

Como ya sabemos, es buena idea primero ver que los datos están en el 
formato correcto. Para este ejercicio solo utilizaremos las lecturas
que conservaron sus pares después de ser filtradas:

~~~ {.bash}
$ for i in trinity_out_dir/*.P.qtrim.gz ; do zcat $i | head ; done 
$ ln -s ./trinity_out_dir/*.P.qtrim.gz .
~~~

Una vez que verificamos que las lecturas están el formato correcto,
vamos a alinear las lecturas y transcritos al genoma utilizando Bowtie2 via TopHat. 

Pueden encontrar el manual en el siguiente [link](https://ccb.jhu.edu/software/tophat/manual.shtml).

## Mapeando las lecturas filtradas al genoma

Primero generaremos un índice de bowtie2 para el genoma:

~~~ {.bash}
$ bowtie2-build Sp_genome.fa Sp_genome 
~~~

~~~ {.output}
Settings:
  Output files: "Sp_genome.*.bt2"
  Line rate: 6 (line is 64 bytes)
  Lines per side: 1 (side is 64 bytes)
  Offset rate: 4 (one in 16)
  FTable chars: 10
  Strings: unpacked
  Max bucket size: default
  Max bucket size, sqrt multiplier: default
  Max bucket size, len divisor: 4
  Difference-cover sample period: 1024
  Endianness: little
  Actual local endianness: little
  Sanity checking: disabled
  Assertions: disabled
  Random seed: 0
  Sizeofs: void*:8, int:4, long:8, size_t:8
Input files DNA, FASTA:
  Sp_genome.fa
Building a SMALL index
Reading reference sizes
  Time reading reference sizes: 00:00:00
Calculating joined length
Writing header
Reserving space for joined string
Joining reference sequences
  Time to join reference sequences: 00:00:00
bmax according to bmaxDivN setting: 107781
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
Entering Ebwt loop
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
Exited Ebwt loop
fchr[A]: 0
fchr[C]: 135086
fchr[G]: 215577
fchr[T]: 296650
fchr[$]: 431126
Exiting Ebwt::buildToDisk()
Returning from initFromVector
Wrote 4340604 bytes to primary EBWT file: Sp_genome.1.bt2
Wrote 107788 bytes to secondary EBWT file: Sp_genome.2.bt2
Re-opening _in1 and _in2 as input streams
Returning from Ebwt constructor
Headers:
    len: 431126
    bwtLen: 431127
    sz: 107782
    bwtSz: 107782
    lineRate: 6
    offRate: 4
    offMask: 0xfffffff0
    ftabChars: 10
    eftabLen: 20
    eftabSz: 80
    ftabLen: 1048577
    ftabSz: 4194308
    offsLen: 26946
    offsSz: 107784
    lineSz: 64
    sideSz: 64
    sideBwtSz: 48
    sideBwtLen: 192
    numSides: 2246
    numLines: 2246
    ebwtTotLen: 143744
    ebwtTotSz: 143744
    color: 0
    reverse: 0
Total time for call to driver() for forward index: 00:00:00
Reading reference sizes
  Time reading reference sizes: 00:00:00
Calculating joined length
Writing header
Reserving space for joined string
Joining reference sequences
  Time to join reference sequences: 00:00:00
  Time to reverse reference sequence: 00:00:00
bmax according to bmaxDivN setting: 107781
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
Split 1, merged 6; iterating...
Splitting and merging
  Splitting and merging time: 00:00:00
Split 1, merged 0; iterating...
Splitting and merging
  Splitting and merging time: 00:00:00
Avg bucket size: 53889.9 (target: 80835)
Converting suffix-array elements to index image
Allocating ftab, absorbFtab
Entering Ebwt loop
Getting block 1 of 8
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
  Sorting block of length 41923 for bucket 1
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 41924 for bucket 1
Getting block 2 of 8
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
  Sorting block of length 79784 for bucket 2
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 79785 for bucket 2
Getting block 3 of 8
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
  Sorting block of length 45293 for bucket 3
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 45294 for bucket 3
Getting block 4 of 8
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
  Sorting block of length 39705 for bucket 4
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 39706 for bucket 4
Getting block 5 of 8
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
  Sorting block of length 51531 for bucket 5
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 51532 for bucket 5
Getting block 6 of 8
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
  Sorting block of length 49697 for bucket 6
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 49698 for bucket 6
Getting block 7 of 8
  Reserving size (80836) for bucket 7
  Calculating Z arrays for bucket 7
  Entering block accumulator loop for bucket 7:
  bucket 7: 10%
  bucket 7: 20%
  bucket 7: 30%
  bucket 7: 40%
  bucket 7: 50%
  bucket 7: 60%
  bucket 7: 70%
  bucket 7: 80%
  bucket 7: 90%
  bucket 7: 100%
  Sorting block of length 77916 for bucket 7
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 77917 for bucket 7
Getting block 8 of 8
  Reserving size (80836) for bucket 8
  Calculating Z arrays for bucket 8
  Entering block accumulator loop for bucket 8:
  bucket 8: 10%
  bucket 8: 20%
  bucket 8: 30%
  bucket 8: 40%
  bucket 8: 50%
  bucket 8: 60%
  bucket 8: 70%
  bucket 8: 80%
  bucket 8: 90%
  bucket 8: 100%
  Sorting block of length 45270 for bucket 8
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 45271 for bucket 8
Exited Ebwt loop
fchr[A]: 0
fchr[C]: 135086
fchr[G]: 215577
fchr[T]: 296650
fchr[$]: 431126
Exiting Ebwt::buildToDisk()
Returning from initFromVector
Wrote 4340604 bytes to primary EBWT file: Sp_genome.rev.1.bt2
Wrote 107788 bytes to secondary EBWT file: Sp_genome.rev.2.bt2
Re-opening _in1 and _in2 as input streams
Returning from Ebwt constructor
Headers:
    len: 431126
    bwtLen: 431127
    sz: 107782
    bwtSz: 107782
    lineRate: 6
    offRate: 4
    offMask: 0xfffffff0
    ftabChars: 10
    eftabLen: 20
    eftabSz: 80
    ftabLen: 1048577
    ftabSz: 4194308
    offsLen: 26946
    offsSz: 107784
    lineSz: 64
    sideSz: 64
    sideBwtSz: 48
    sideBwtLen: 192
    numSides: 2246
    numLines: 2246
    ebwtTotLen: 143744
    ebwtTotSz: 143744
    color: 0
    reverse: 1
Total time for backward call to driver() for mirror index: 00:00:00
~~~
 
~~~ {.bash}
$ ls *bt2
~~~ 

~~~ {.output}
Sp_genome.1.bt2        
Sp_genome.2.bt2              
Sp_genome.3.bt2       
Sp_genome.4.bt2           
Sp_genome.rev.1.bt2          
Sp_genome.rev.2.bt2
~~~

Creamos hiperlinks para los archivos pareados en nuestro directorio local:
~~~ {.bash}
$ ln -s ./trinity_out_dir/*P.qtrim.gz .
~~~

Usamos tophat2 para mapear las lecturas. Este programa nos permite dividir lecturas
que atraviesan sitios de splicing:

~~~ {.bash}
$ tophat2 -I 300 -i 20 Sp_genome \
 Sp_log.left.fq.gz.P.qtrim.gz,Sp_hs.left.fq.gz.P.qtrim.gz,Sp_ds.left.fq.gz.P.qtrim.gz,Sp_plat.left.fq.gz.P.qtrim.gz \
 Sp_log.right.fq.gz.P.qtrim.gz,Sp_hs.right.fq.gz.P.qtrim.gz,Sp_ds.right.fq.gz.P.qtrim.gz,Sp_plat.right.fq.gz.P.qtrim.gz
~~~ 

~~~ {.output}
[2017-01-24 13:53:23] Beginning TopHat run (v2.1.1)
-----------------------------------------------
[2017-01-24 13:53:23] Checking for Bowtie
		  Bowtie version:	 2.2.8.0
[2017-01-24 13:53:23] Checking for Bowtie index files (genome)..
[2017-01-24 13:53:23] Checking for reference FASTA file
[2017-01-24 13:53:23] Generating SAM header for Sp_genome
[2017-01-24 13:53:23] Preparing reads
	 left reads: min. length=25, max. length=68, 329455 kept reads (431 discarded)
	right reads: min. length=25, max. length=68, 329725 kept reads (161 discarded)
[2017-01-24 13:53:31] Mapping left_kept_reads to genome Sp_genome with Bowtie2
[2017-01-24 13:53:44] Mapping left_kept_reads_seg1 to genome Sp_genome with Bowtie2 (1/2)
[2017-01-24 13:53:44] Mapping left_kept_reads_seg2 to genome Sp_genome with Bowtie2 (2/2)
[2017-01-24 13:53:44] Mapping right_kept_reads to genome Sp_genome with Bowtie2
[2017-01-24 13:53:57] Mapping right_kept_reads_seg1 to genome Sp_genome with Bowtie2 (1/2)
[2017-01-24 13:53:57] Mapping right_kept_reads_seg2 to genome Sp_genome with Bowtie2 (2/2)
[2017-01-24 13:53:57] Searching for junctions via segment mapping
	Coverage-search algorithm is turned on, making this step very slow
	Please try running TopHat again with the option (--no-coverage-search) if this step takes too much time or memory.
[2017-01-24 13:54:02] Retrieving sequences for splices
[2017-01-24 13:54:02] Indexing splices
Building a SMALL index
[2017-01-24 13:54:02] Mapping left_kept_reads_seg1 to genome segment_juncs with Bowtie2 (1/2)
[2017-01-24 13:54:02] Mapping left_kept_reads_seg2 to genome segment_juncs with Bowtie2 (2/2)
[2017-01-24 13:54:02] Joining segment hits
[2017-01-24 13:54:03] Mapping right_kept_reads_seg1 to genome segment_juncs with Bowtie2 (1/2)
[2017-01-24 13:54:03] Mapping right_kept_reads_seg2 to genome segment_juncs with Bowtie2 (2/2)
[2017-01-24 13:54:03] Joining segment hits
[2017-01-24 13:54:04] Reporting output tracks
-----------------------------------------------
[2017-01-24 13:54:30] A summary of the alignment counts can be found in ./tophat_out/align_summary.txt
[2017-01-24 13:54:30] Run complete: 00:01:07 elapsed
~~~

Exploramos el resultado, el cuál es un archivo tipo SAM. 

~~~ {.bash}
$ samtools view tophat_out/accepted_hits.bam | head
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
Convirtamos el archivo BAM a SAM usando samtools:

~~~ {.bash}
$ samtools view tophat_out/accepted_hits.bam > samtools view tophat_out/accepted_hits.sam
~~~

No abrimos los archivos SAM ya que, dado que esta en formato binario, son ilegible. 
Sin embargo, tenemos que realizar dos últimos pasos visualizar 
los resultados:

~~~ {.bash}
$ samtools sort tophat_out/accepted_hits.bam accepted_hits
$ samtools index tophat_out/accepted_hits.bam
~~~

El primer paso ordena los resultados por sus coordenadas y el segundo crea índices
para hacer más rápida la visualización usando un navegador. 

> ## Tarea - Alineando las lecturas filtradas al transcriptoma {.challenge}
>
> Hemos alineado las lecturas al genoma pero queremos alinearlas también directamente
> al transcriptoma. La vamos a revisar el manual de TopHat2 y usar las opciones
> que nos permite mapear lecturas directamente a transcriptomas. 
> 
> * **Pista:** No podrán utilizar el índice generado previamente.
>
> ### Solución
>
> ~~~ {.bash}
> $ bowtie2-build Trinity.fasta Trinity_assembly_Sp
> $ tophat2 -I 300 -i 20 Trinity_assembly_Sp Sp_ds.left.fq.gz.P.qtrim.gz,Sp_hs.left.fq.gz.P.qtrim.gz,Sp_log.left.fq.gz.P.qtrim.gz,Sp_plat.left.fq.gz.P.qtrim.gz Sp_ds.right.fq.gz.P.qtrim.gz,Sp_hs.right.fq.gz.P.qtrim.gz,Sp_log.right.fq.gz.P.qtrim.gz,Sp_plat.right.fq.gz.P.qtrim.gz
> ~~~










