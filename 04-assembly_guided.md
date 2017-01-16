---
layout: page
title: Análisis transcriptómicos
subtitle: Ensamble de transcriptomas guiado
minutes: 5
---
> ## Objetivos de aprendizaje {.objectives}
>
> *  Aprender como funciona el ensamble de RNA-Seq *de novo*.
> *  Aprender a usar Trinity para ensamblar datos de novo.

Para empezar, instalaremos Trinity en sus computadoras:

~~~ {.bash}
$ sudo apt-get install trinityrnaseq
~~~

Y bajaremos los siguientes archivos fastq:

[Sp_ds.left.fq.gz](datasets/Sp_ds.left.fq.gz)
  
[Sp_ds.right.fq.gz](datasets/Sp_ds.right.fq.gz)
   
[Sp_hs.right.fq.gz](datasets/Sp_hs.right.fq.gz)
 
[Sp_hs.left.fq.gz](datasets/Sp_hs.left.fq.gz)
    
[Sp_plat.left.fq.gz](datasets/Sp_plat.left.fq.gz)
  
[Sp_plat.right.fq.gz](datasets/Sp_plat.right.fq.gz)
 
[Sp_log.left.fq.gz](datasets/Sp_log.left.fq.gz)
   
[Sp_log.right.fq.gz](datasets/Sp_log.right.fq.gz) 

Exploremos las lecturas en cada archivo:

~~~ {.bash}
$ for i in S*fq.gz ; do zcat $i | head ; wait ; done 
~~~

Este comando nos mostrará las primeras 10 líneas de cada archivo
de manera iterativa. Esta estructura es conocida como un `for loop`. 
Un for loop nos permite ejecutar un comando en varios archivos 
de manera secuencial. Es extremadamente útil cuando tenemos varios 
archivos.

Una vez que verificamos que los datos son correctos,
usaremos el programa Trinity para ensamblar los transcritos.

Primero empezaremos usando un comando genérico:

~~~ {.bash}
$ Trinity --seqType fq --SS_lib_type RF \
--left Sp_log.left.fq.gz Sp_hs.left.fq.gz \
--right Sp_log.right.fq.gz Sp_hs.right.fq.gz \
--CPU 2 --max_memory 1G
~~~

Este comando tardará aproximadamente 15 minutos en ensamblar un 
transcriptoma. 

Las opciones (banderas) que hemos utilizado son las siguientes:

* --seqType fq - Indicamos que estamos usando archivos tipo fastq
* --SS_lib_type RF - Indicamos que la librería fue construida usando 
lecturas en pares en la orientación RF (reverse-forward)
* --left - Lecturas del lado izquierdo (o R1)
* --right - Lecturas del lado derecho (o R2)
* --CPU 2 - Utilizar 2 CPUs
* --max_memory 1G - Indicar que la memoria RAM máxima a usar es de 1GB

Estás son las opciones más esenciales para llevar a cabo el análisis. 
Por ello es muy importante saber la orientación de la librería 
dado el protocolo que se usó para generarla. 

> ## Orientación de librerías {.callout}
>
> La orientación de las librerías es un parámetro esencial para el 
> ensamble, ya que esto dictará la orientación de los transcritos. 
> 
> Esta es una clave que les permitirá decidir que parámetro utilizar
> para indicar a Trinity la orientación de su librería. 
> 
> ![Orientación de librerías](fig/strand_specificity.jpg)
>

> ## ¿Porqué mezclamos librerías? {.challenge}
>
> En el comando usado en la parte superior podemos ver que 
> Estamos ensamblando el transcriptoma con la librería log 
> así como la librería hs. ¿Por qué no ensamblarlas de forma
> independiente? ¿Qué ventajas o desventajas crees que tendría 
> ensamblarlas juntas o separadas?

Uno de los problemas comunes con Trinity es la falta de memoria. 
Un análisis típico de Trinity requiere ~1 hora y ~1GB de RAM 
por ~1 millón de lecturas en pares (paired-end). Es por ello 
que estos análisis se efectúan en un servidor con grandes 
cantidades de memoria y se le deja correr por varios días. 

Trinity tiene muchas otras opciones las cuales podemos 
explorar escribiendo (en una terminal diferente a la que estamos
usando para nuestro análisis):

~~~ {.bash}
$ Trinity 
~~~
~~~ {.output}
###############################################################################
#
#     ______  ____   ____  ____   ____  ______  __ __
#    |      ||    \ |    ||    \ |    ||      ||  |  |
#    |      ||  D  ) |  | |  _  | |  | |      ||  |  |
#    |_|  |_||    /  |  | |  |  | |  | |_|  |_||  ~  |
#      |  |  |    \  |  | |  |  | |  |   |  |  |___, |
#      |  |  |  .  \ |  | |  |  | |  |   |  |  |     |
#      |__|  |__|\_||____||__|__||____|  |__|  |____/
#
###############################################################################
#
# Required:
#
#  --seqType <string>      :type of reads: ( fa, or fq )
#
#  --max_memory <string>      :suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting, etc)
#                            provied in Gb of RAM, ie.  '--max_memory 10G'
#
#  If paired reads:
#      --left  <string>    :left reads, one or more file names (separated by commas, no spaces)
#      --right <string>    :right reads, one or more file names (separated by commas, no spaces)
#
#  Or, if unpaired reads:
#      --single <string>   :single reads, one or more file names, comma-delimited (note, if single file contains pairs, can use flag: --run_as_paired )
#
####################################
##  Misc:  #########################
#
#  --SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.
#                                   if paired: RF or FR,
#                                   if single: F or R.   (dUTP method = RF)
#                                   See web documentation.
#
#  --CPU <int>                     :number of CPUs to use, default: 2
#  --min_contig_length <int>       :minimum assembled contig length to report
#                                   (def=200)
#
#  --long_reads <string>           :fasta file containing error-corrected or circular consensus (CCS) pac bio reads
#
#  --genome_guided_bam <string>    :genome guided mode, provide path to coordinate-sorted bam file.
#                                   (see genome-guided param section under --show_full_usage_info)
#
#  --jaccard_clip                  :option, set if you have paired reads and
#                                   you expect high gene density with UTR
#                                   overlap (use FASTQ input file format
#                                   for reads).
#                                   (note: jaccard_clip is an expensive
#                                   operation, so avoid using it unless
#                                   necessary due to finding excessive fusion
#                                   transcripts w/o it.)
#
#  --trimmomatic                   :run Trimmomatic to quality trim reads
#                                        see '--quality_trimming_params' under full usage info for tailored settings.
#                                  
#
#  --normalize_reads               :run in silico normalization of reads. Defaults to max. read coverage of 50.
#                                       see '--normalize_max_read_cov' under full usage info for tailored settings.
#     
#  --no_distributed_trinity_exec   :do not run Trinity phase 2 (assembly of partitioned reads), and stop after generating command list.
#
#
#  --output <string>               :name of directory for output (will be
#                                   created if it doesn't already exist)
#                                   default( your current working directory: "/home/usuario/trinity_out_dir" 
#                                    note: must include 'trinity' in the name as a safety precaution! )
#  
#  --full_cleanup                  :only retain the Trinity fasta file, rename as ${output_dir}.Trinity.fasta
#
#  --cite                          :show the Trinity literature citation
#
#  --verbose                       :provide additional job status info during the run.
#
#  --version                       :reports Trinity version (v2.0.6) and exits.
#
#  --show_full_usage_info          :show the many many more options available for running Trinity (expert usage).
#
#
###############################################################################
#
#  *Note, a typical Trinity command might be:
#
#        Trinity --seqType fq --max_memory 50G --left reads_1.fq  --right reads_2.fq --CPU 6
#
#
#    and for Genome-guided Trinity:
#
#        Trinity --genome_guided_bam rnaseq_alignments.csorted.bam --max_memory 50G
#                --genome_guided_max_intron 10000 --CPU 6
#
#     see: /usr/lib/trinityrnaseq/sample_data/test_Trinity_Assembly/
#          for sample data and 'runMe.sh' for example Trinity execution
#
#     For more details, visit: http://trinityrnaseq.github.io
#
###############################################################################
~~~

Una de las opciones más interesante es la normalización *in silico* 
ejecutada a través de la opción `--normalize_reads`. Esta opción es
especialmente útil para sets de datos de >300 millones de pares. 

Esta cambia la probabilidad de las lecturas dada su abundancia en 
relación a la abundancia de sus kmeros. Este tipo de normalización 
resulta en ahorros significativos de tiempo de análisis. 

Nuestro trabajo debe haber terminado, revisemos los transcritos 
generados, los cuales se encuentran en el archivo `Trinity.fasta`:

~~~ {.bash}
$ head trinity_out_dir/Trinity.fasta
~~~

Observamos que los resultados son secuencias en formato `fasta`.

~~~ {.output}
>TR1|c0_g1_i1 len=344 path=[322:0-343] [-1, 322, -2]
GAAAAAATTTATCATCGTCAAATCTTTAAGCAGTTTCTGACTAACAAAATTTTGAAAGAT
CCAAAACAAAGAAGGTTTTTTAAAATGACGGACTTGCACGATTTGTTTACGTTAGGCGAT
AACAAGACTGAGGGCACTGAGACAGGCAGCATGTTTTTGGGATCTGAACGAGTACTTCGA
AAGGATAATTCCTCAAGAAATGGCAATGAAGCTGAAGATATTCCAGCTCGTGACCGAAAA
AAGCACAAAATTCACGACAAAGGTAAAAAAGTTAACAGCTCCAAAGTGTTTGAAAAAATG
GGGATTGCATCGATGGAAAAGTATAAACCACCGCAAGAGTCAAA
>TR1|c3_g1_i1 len=229 path=[207:0-228] [-1, 207, -2]
ATTGGTGATGAAATGGGACTAGGAAAAACTATTCAAATAGTATCTTTCCTTTCGTCTTTG
CATCACTCTGGCAAATTTCAAAAGCCTGCACTTATCGTTTGTCCAGCGACTTTAATGAAG
~~~

## Analisando las estadísticas del transcriptoma ensamblado

Podemos capturar algunas estadística acerca de este ensamble
usando un programa que es parte de Trinity:

~~~ {.bash}
$ /usr/lib/trinityrnaseq/util/TrinityStats.pl trinity_out_dir/Trinity.fasta 
~~~

El cuál generará los siguientes datos: 

~~~ {.output}
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':    333
Total trinity transcripts:    338
Percent GC: 39.27

########################################
Stats based on ALL transcript contigs:
########################################

    Contig N10: 2605
    Contig N20: 2016
    Contig N30: 1596
    Contig N40: 1394
    Contig N50: 1160

    Median contig length: 496
    Average contig: 773.81
    Total assembled bases: 261549


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

    Contig N10: 2576
    Contig N20: 1917
    Contig N30: 1591
    Contig N40: 1377
    Contig N50: 1143

    Median contig length: 494
    Average contig: 761.29
    Total assembled bases: 253508
~~~

Este resumen nos indica: 

* Cuantos genes y transcritos fueron ensamblados
* El contenido de GC
* Estadísticas sobre el tamaño medio de los contigs
* Número de bases ensambladas
* Estadísticas sobre el tamaño medio de la isoforma más larga de cada gen

De este resumen, es particularmente importante el concepto de Contig N50, N40, etc. 
Por ejemplo, N50 indica el tamaño del contig medio (o 50%) cuando 
todos los contigs son ordenados por tamaño. N40 es el 40% etc. Esta medida nos
ayuda a no basarnos simplemente en el contig más largo y nos permite observar
como va incrementando la longitud en todos los transcritos. 

Finalmente, realizaremos un blast de las primeras 5 secuencias 
para identificar de que organismo provienen. Navegen hacia:

[NCBI Blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi)

> ## ¿Qué tipo de Blast debo utilizar? {.challenge}
>
> Dado que tenemos datos nucleotídicos pero queremos buscar en las
> bases de datos de proteínas de NCBI.
> ¿Qué tipo de blast debemos utilizar?

> ## Tarea {.challenge}
>
> Si se pudieron dar cuenta, los datos que ensamblamos son los datos crudos, es decir
> sin filtrar por calidad. Afortunadamente Trinity incluye una opción para usar 
> Trimmomatic de manera sencilla. 
>
> Su tarea consiste en re-ensamblar el transcriptoma esta vez utilizando todos los 
> archivos proporcionados, pero filtrando la calidad de las secuencias. 
> Añadan a su repositorio un archivo con las primeras 20 secuencias del 
> transcriptoma con el nombre "Ensamble_Final_NOMBRE_APELLIDO.fasta". 
> Una vez más, añadanlo a su repositorio en cuanto empecemos la siguiente clase.
>
> **Nota: Por favor agreguen también el reporte generado por el programa `TrinityStats.pl` 
en su transcriptoma completo con el nombre "Ensamble_Final_NOMBRE_APELLIDO_Stats.txt".**
>
> Conserven el archivo fasta completo, ya que lo necesitarán para las siguientes
> prácticas.

