---
layout: page
title: Análisis transcriptómicos
subtitle: Anotación de transcritos
minutes: 5
---
> ## Objetivos de aprendizaje {.objectives}
>
> *  Aprender a usar Trinotate y Blast2GO para anotar un transcriptoma *de novo*

> ## Datos requeridos {.prereq}
>
> Para proceder con esta práctica, se requieren resultados de las clases anteriores:
>
> *  Archivo Trinity.fasta con ensamble de secuencias *de novo*

## Anotación usando Trinotate

Trinotate es el programa de anotación generada por los programadores
autores de Trinity. 

Utiliza tres tipos de evidencia para anotar transcritos:

*  Búsqueda de homología con secuencias protéicas conocidas (BLAST+/SwissProt)
*  Identificación de dominios protéicos (HMMER/PFAM)
*  Identificación de péptidos señal (signalP)
*  Identificación de dominios transmembranales (PFAM)
*  Anotación de ontología GO (Gene Ontology) (eggNOG/GO/KEGG)

Estas anotaciones se integran en el una base datos generada usando 
el protocolo SQLite, permitiéndonos accesar estos datos rápidamente. 

Comenzaremos descargando los archivos que requerimos:

[Trinotate 3.0.0](https://github.com/Trinotate/Trinotate/archive/v3.0.0.tar.gz)

Installando SQLite:

~~~ {.bash}
sudo apt-get install sqlite3 libsqlite3-dev
~~~

Descargamos las bases de datos que necesitamos. *Los links de abajo les permitirán bajar las bases de datos pero durante la clase las copiarán de una
memoria para hacer este proceso más rápido*.

[SwissProt](https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/uniprot_sprot.pep.gz)

[Pfam-A](https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Pfam-A.hmm.gz)

También descargamos los ejecutables de los programas:

*  [signalP v4](signalp-4.1c.Linux.tar.gz)
*  [tmhmm v2](programs/tmhmm-2.0c.Linux.tar.gz)
*  [RNAMMER](programs/rnammer-1.2.src.tar.Z)

### Paso 1. Extraer marcos de lectura abierta (ORFs) de todos los transcritos.

Por el momento, en su archivo Trinity.fasta solo contiene secuencias de 
nucleótidos en formato fasta. Vamos a generar (recuperar) las secuencias de proteinas en estos 
transcritos usando el programa 
(TransDecoder)[https://transdecoder.github.io/]. Este ya está incluido en el 
paquete que acompaña a Trinity.

TransDecoder identifica secuencias con potencial codificante basándose en:

*  La longitud mínima del ORF en cada transcrito.
*  Un log-score > 0. Este es mayor cuando el péptido se codifica en el marco 1 de lectura.

Si un ORF esta contenido dentro de otro, se reportará el más largo. Sin 
embargo un mismo transcrito tiene ORFs múltiples.

Identificamos los ORFs más largos:

~~~ {.bash}
$ ~/TransDecoder/TransDecoder.LongOrfs -t trinity_out_dir/Trinity.fasta
~~~

TransDecoder genera los siguientes archivos principales:

*  transcripts.fasta.transdecoder.pep : peptide sequences for the final candidate ORFs; all shorter candidates within longer ORFs were removed.
*  transcripts.fasta.transdecoder.cds  : nucleotide sequences for coding regions of the final candidate ORFs
*  transcripts.fasta.transdecoder.gff3 : positions within the target transcripts of the final selected ORFs
*  transcripts.fasta.transdecoder.bed  : bed-formatted file describing ORF positions, best for viewing using GenomeView or IGV.

Revisemos los resultados:
~~~ {.output}
$ less transcripts.fasta.transdecoder.pep
~~~

Este archivo esta en formato fasta y contiene los ORFs más largos de cada
transcrito.

### Paso 2. Ejecutar programas de información para la anotación.

Vamos a ejecutar los distintos programas listados previamente que nos proveen 
con información para la anotación. 

Primero, creemos un directorio llamado `Anotacion` para guardar nuestros 
datos y copiemos los archivos con el transcriptoma y el genoma al mismo:

~~~ {.bash}
$ mkdir Anotacion
$ cd Anotacion
$ cp ../transcripts.fasta.transdecoder.pep . 
$ cp ../trinity_out_dir/Trinity.fasta . 
$ mv ~/Downloads/uniprot_sprot.pep.gz . 
$ mv ~/Downloads/Pfam-A-hmm.gz . 
~~~

Preparamos las bases de datos de BLAST:
~~~ {.bash}
$ gunzip uniprot_sprot.pep.gz
$ makeblastdb -in uniprot_sprot.pep -dbtype prot
~~~

Y de HMMER:
~~~ {.bash}
$ gunzip Pfam-A.hmm.gz
$ hmmpress Pfam-A.hmm
~~~

Ejecutamos entonces los programas 

~~~ {.bash}
$ blastx -query Trinity.fasta -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
$ blastp -query transdecoder.pep -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
~~~

> ## ¿Doble BLAST? {.challenge}
>
>  ¿Porqué alinear nucleótidos y proteínas? ¿No sería suficiente con hacerlo 
> con una sola de estas opciones?
> 

~~~ {.bash}
$ hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm transdecoder.pep > pfam.log
$ signalp -f short -n signalp.out transdecoder.pep
$ tmhmm --short < transdecoder.pep > tmhmm.out
$ $TRINOTATE_HOME/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer /usr/bin/software/rnammer_v1.2/rnammer
~~~

Cada uno de estos programas produce un archivo tabulado que indica en que parte
de cada transcrito se encuentran algunas de estos dominios, regiones de identidad o señales. Pueden explorarlos usando el comando `head`

### Paso 3. Agregar los resultados de la anotación a la base de datos.

Primero creamos un archivo que especifique que transcrito corresponde a que 
proteina:

~~~ {.bash}
TRINITY_HOME/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta >  Trinity.fasta.gene_trans_map
~~~

Iniciamos la base de datos importando el transcriptoma y ORFs:

~~~ {.bash}
$ Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta Trinity.fasta --transdecoder_pep transdecoder.pep
~~~ 

Las secuencias homólogas y dominios encontrados:
~~~ {.bash}
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
Trinotate Trinotate.sqlite LOAD_signalp signalp.out
~~~ 

Y finalmente generamos un reporte:
~~~ {.bash}
Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls
~~~ 

## Anotación usando Blast2GO

El uso de Blast2GO es bastante sencillo, aunque de no tener la versión PRO
puede llegar a ser muy lento. 

