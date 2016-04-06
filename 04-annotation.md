---
layout: page
title: Análisis transcriptómicos
subtitle: Anotación de transcritos
minutes: 5
---
> ## Objetivos de aprendizaje {.objectives}
>
> *  Aprender a usar Trinotate para anotar un transcriptoma *de novo*

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

Comenzamos instalando SQLite:

~~~ {.bash}
sudo apt-get install sqlite3 libsqlite3-dev
~~~

Descargamos las bases de datos y programas que necesitamos. *Los links de abajo les permitirán bajar las bases de datos pero durante la clase las copiarán de una
memoria para hacer este proceso más rápido*.

[Trinotate 3.0.0](https://github.com/Trinotate/Trinotate/archive/v3.0.0.tar.gz)

[SwissProt](https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/uniprot_sprot.pep.gz)
[Pfam-A](https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Pfam-A.hmm.gz)

También descargamos los ejecutables de los programas:

*  [signalP v4](programs/signalp-4.1c.Linux.tar.gz)
*  [tmhmm v2](programs/tmhmm-2.0c.Linux.tar.gz)

Descomprimimos estos archivos dando doble click en la carpeta `Downloads`. 
Los dejaremos en esta carpeta para los propósitos de esta practica. Editaremos 
el archivo signalP haciendo doble click en el archivo llamado `signalP`. 
Cambiaremos la línea:

~~~ {.output}
ENV{SIGNALP} = '/usr/cbs/bio/src/signalp-4.1'
~~~ 

por 

~~~ {.output}
ENV{SIGNALP} = '/home/usuario/Downloads/signalp-4.1'
~~~ 

y en el archivo tmhmm y thmhmmformat.pl 

~~~ {.output}
!/usr/local/bin/perl
~~~ 

por

~~~ {.output}
!/usr/bin/perl
~~~ 

También descargaremos y descomprimiremos el templato de la base de datos del Broad Institute. 
*Esta base de datos también esta en la memoria*. 

~~~ {.bash}
$ wget "https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Trinotate_v3.sqlite.gz" -O Trinotate.sqlite.gz
gunzip Trinotate.sqlite.gz
~~~

### Paso 1. Extraer marcos de lectura abierta (ORFs) de todos los transcritos.

Por el momento, en su archivo Trinity.fasta solo contiene secuencias de 
nucleótidos en formato fasta. Vamos a generar (predecir) las secuencias de proteinas en estos 
transcritos usando el programa 
[TransDecoder](https://transdecoder.github.io/). Este ya está incluido en el 
paquete que acompaña a Trinity.

TransDecoder identifica secuencias con potencial codificante basándose en:

*  La longitud mínima del ORF en cada transcrito.
*  Un log-score > 0. Este es mayor cuando el péptido se codifica en el marco 1 de lectura.

Si un ORF esta contenido dentro de otro, se reportará el más largo. Sin 
embargo un mismo transcrito tiene ORFs múltiples.

Identificamos los ORFs más largos:

~~~ {.bash}
$ TransDecoder.LongOrfs -t trinity_out_dir/Trinity.fasta
$ TransDecoder.Predict -t trinity_out_dir/Trinity.fasta
~~~

TransDecoder genera los siguientes archivos principales:

*  **Trinity.fasta.transdecoder.pep**: Secuencias peptídicas de los ORFs identificados.
*  **Trinity.fasta.transdecoder.mRNA**: Secuencias nucleotídicas de estos ORFs.

Revisemos los resultados:

~~~ {.bash}
$ more Trinity.fasta.transdecoder.pep
$ more Trinity.fasta.transdecoder.mRNA
~~~

Estos archivos están en formato fasta y contiene los ORFs más largos de cada
transcrito.

### Paso 2. Ejecutar programas de información para la anotación.

Vamos a ejecutar los distintos programas listados previamente que nos proveen 
con información para la anotación. 

Primero, creemos un directorio llamado `Anotacion` para guardar nuestros 
datos y copiemos los archivos con el transcriptoma y el genoma al mismo:

~~~ {.bash}
$ mkdir Anotacion
$ cd Anotacion
~~~

Copiamos los transcritos y los archivos generados por TransDecoder a este 
directorio:

~~~ {.bash}
$ cp ../Trinity.fasta.transdecoder.pep . 
$ cp ../trinity_out_dir/Trinity.fasta . 
$ mv ~/Downloads/uniprot_sprot.pep.gz . 
$ mv ~/Downloads/Pfam-A-hmm.gz . 
~~~

Movemos las bases de datos `uniprot_sprot.pep.gz` y `Pfam-A.hmm.gz` a este directorio. 

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

Ejecutamos entonces la búsqueda de homólogos por BLAST:

~~~ {.bash}
$ blastx -query Trinity.fasta -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6 &
$ blastp -query Trinity.fasta.transdecoder.pep -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6 &
~~~

> ## ¿Doble BLAST? {.challenge}
>
>  ¿Porqué alinear nucleótidos y proteínas? ¿No sería suficiente con hacerlo 
> con una sola de estas opciones?
> 

~~~ {.bash}
$ hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm Trinity.fasta.transdecoder.pep > pfam.log &
$ ~/Downloads/signalp-4.1/signalp -f short -n signalp.out Trinity.fasta.transdecoder.pep > signalp.log &
$ ~/Downloads/tmhmm-2.0c/tmhmm --short < Trinity.fasta.transdecoder.pep > tmhmm.out &
~~~

Cada uno de estos programas produce un archivo tabulado que indica en que parte
de cada transcrito se encuentran algunas de estos dominios, regiones de identidad o señales. Pueden explorarlos usando el comando `head`. 

~~~ {.bash}
$ head blastx.outfmt6 blastp.outfmt6 TrinotatePFAM.out signalp.out tmhmm.out
~~~
### Paso 3. Agregar los resultados de la anotación a la base de datos.

Primero creamos un archivo que especifique que transcrito corresponde a que 
proteina:

~~~ {.bash}
/usr/lib/trinityrnaseq/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta >  Trinity.fasta.gene_trans_map
~~~

Iniciamos la base de datos importando el transcriptoma y ORFs:

~~~ {.bash}
$ ~/Downloads/Trinotate-3.0.0/Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta ../trinity_out_dir/Trinity.fasta --transdecoder_pep Trinity.fasta.transdecoder.pep
~~~ 

Las secuencias homólogas y dominios encontrados:

~~~ {.bash}
~/Downloads/Trinotate-3.0.0/Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
~/Downloads/Trinotate-3.0.0/Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
~/Downloads/Trinotate-3.0.0/Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out 
~/Downloads/Trinotate-3.0.0/Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out 
~/Downloads/Trinotate-3.0.0/Trinotate Trinotate.sqlite LOAD_signalp signalp.out 
~~~ 

Y finalmente generamos un reporte:

~~~ {.bash}
~/Downloads/Trinotate-3.0.0/Trinotate Trinotate.sqlite report > Trinotate_annotation_report.xls
~~~ 

Este reporte tiene las siguientes columnas:

1. \#gene_id
1. transcript_id
2. sprot_Top_BLASTX_hit
3. RNAMMER
4. prot_id
5. prot_coords
6. sprot_Top_BLASTP_hit
7. pep_BLASTX
8. pep_BLASTP
9. Pfam
10. SignalP
11. TmHMM
12. eggnog
13. Kegg
14. gene_ontology_blast
15. gene_ontology_pfam
16. transcript
17. peptide

> ## Tarea - Dominios de proteínas {.challenge}
>
> Usando su reporte de Trinotate, generen una lista de los dominios de PFAM, ya sea
> en R o en excel. Generen una gráfica de barras de los 10 dominios más abundantes en la 
> anotación de su transcriptoma. Debajo de la gráfica discutan brevemente (1 párrafo de máximo 
> 10 líneas) cual es la relevancia biológica de conocer estos dominios. 
> 
> Añadan este reporte (la 
> gráfica y el párrafo) a su repositorio de análisis de transcriptomas. El formato del 
> reporte es libre así que lo pueden realizar en Rmd o incluso en Word. Si eligen el último
> formato **por favor súbanlo en PDF, no en Word**. 



