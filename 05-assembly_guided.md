---
layout: page
title: Análisis transcriptómicos
subtitle: Ensamble de transcriptomas guiado
minutes: 5
---
> ## Objetivos de aprendizaje {.objectives}
>
> *  Aprender como funciona el ensamble guiado.
> *  Aprender a usar Cufflinks para ensamblar datos de novo.

Usaremos las mismas lecturas que ya descargamos en la lección de [ensamble de novo](02-assembly_denovo.html). 

Para empezar, creemos un directorio para trabajar con nuestros datos:

~~~ {.bash}
$ mkdir CUFFLINKS
$ cd CUFFLINKS/
~~~ 

Generemos un acceso directo para el archivo del genoma:

~~~ {.bash}
$ ln -s ../MAP/Sp_genome.fa
~~~

Preparemos el archivo `Sp_genome.fa` para alinear usando TopHat2:

~~~ {.bash}
$ bowtie2-build Sp_genome.fa Sp_genome
~~~

Y alineamos las lecturas usando Top Hat

~~~ {.bash}
$ tophat2 -I 1000 -i 20 --bowtie1 --library-type fr-firststrand -o tophat.Sp_ds.dir Sp_genome Sp_ds.left.fq.gz Sp_ds.right.fq.gz
~~~

Las opciones (banderas) y archivos que hemos utilizado son las siguientes:

* -I 1000 - Tamaño máximo del intron. 
* -i 20 - Tamaño mínimo del intron. 
* --library-type fr-firststrand - Indicamos que la librería fue construida usando 
lecturas en pares en la orientación RF (reverse-forward)
* -o tophat.Sp_ds.dir Directorio en donde queremos que se guarden nuestros resultados (output)
* Sp_ds.left.fq.gz Lecturas del lado izquierdo (o R2)
* Sp_ds.right.fq.gz Lecturas del lado derecho (o R2)

Una vez finalizado, renombraremos el resultado del alineamiento (bam) para que 
concuerde con el nombre de la muestra analizada:

~~~ {.bash}
$ mv tophat.Sp_ds.dir/accepted_hits.bam tophat.Sp_ds.dir/Sp_ds.bam
~~~ 

Y crearemos un índice para que podamos visualizar este archivo en nuestra siguiente 
lección.

~~~ {.bash}
$ samtools index tophat.Sp_ds.dir/Sp_ds.bam
~~~ 

Reconstruiremos los transcritos para esta muestra usando Cufflinks:

~~~ {.bash}
$ cufflinks --overlap-radius 1 \
             --library-type fr-firststrand \
             -o cufflinks.Sp_ds.dir tophat.Sp_ds.dir/Sp_ds.bam
~~~ 

Y una vez más lo renombraremos para poder identificar de que muestra vino:

~~~ {.bash}
$ mv cufflinks.Sp_ds.dir/transcripts.gtf cufflinks.Sp_ds.dir/Sp_ds.transcripts.gtf
~~~ 
 
Esto finaliza el análisis para esta muestra, ahora deberás hacer los mismo para las 
demás muestras, por ejemplo, para la muestra Sp_log:

~~~ {.bash}
$ tophat2 -I 1000 -i 20 --library-type fr-firststrand \
           -o tophat.Sp_log.dir Sp_genome \
           Sp_log.left.fq.gz Sp_log.right.fq.gz

$ mv tophat.Sp_log.dir/accepted_hits.bam tophat.Sp_log.dir/Sp_log.bam

$ samtools index tophat.Sp_log.dir/Sp_log.bam

$ cufflinks --overlap-radius 1 \
             --library-type fr-firststrand \
             -o cufflinks.Sp_log.dir tophat.Sp_log.dir/Sp_log.bam

$ mv cufflinks.Sp_log.dir/transcripts.gtf cufflinks.Sp_log.dir/Sp_log.transcripts.gtf
~~~ 

**NOTA:** Los comandos deberán ejecutar de manera secuencial, es decir, deberás esperar
a que terminen los anteriores para continuar.

Finalmente vamos a crear un archivo *maestro* con coordenadas gtf, que utilizaremos para
visualizar nuestro transcriptoma:

~~~ {.bash}
echo cufflinks.Sp_ds.dir/Sp_ds.transcripts.gtf > assemblies.txt

echo cufflinks.Sp_hs.dir/Sp_hs.transcripts.gtf >> assemblies.txt

echo cufflinks.Sp_log.dir/Sp_log.transcripts.gtf >> assemblies.txt

echo cufflinks.Sp_plat.dir/Sp_plat.transcripts.gtf >> assemblies.txt
~~~
  
 Verifica que este archivo contiene todos los nombres de los archivos `.gtf` para que 
 se incluyan en el merge:

~~~ {.bash}
cat assemblies.txt 
~~~

~~~ {.output}
cufflinks.Sp_ds.dir/Sp_ds.transcripts.gtf
cufflinks.Sp_hs.dir/Sp_hs.transcripts.gtf
cufflinks.Sp_log.dir/Sp_log.transcripts.gtf
cufflinks.Sp_plat.dir/Sp_plat.transcripts.gtf
~~~

Entonces estamos listos para mezclar los transcriptomas usando `cuffmerge`:

~~~ {.bash}
cuffmerge -s Sp_genome.fa assemblies.txt
~~~ 

El grupo de transcritos combinados (merged) esta ahora en el archivo 'merged_asm/merged.gtf'.
Usaremos este archivo para visualizarlo en el navegador [IGV](http://software.broadinstitute.org/software/igv/).
