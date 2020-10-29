---
layout: page
title: Análisis transcriptómicos
subtitle: Ensamble de transcriptomas guiado
minutes: 5
---
> ## Objetivos de aprendizaje {.objectives}
>
> *  Aprender como funciona el ensamble guiado.
> *  Aprender a usar Stringtie para ensamblar datos de novo.

Usaremos las mismas lecturas que ya descargamos en la lección de [ensamble de novo](02-assembly_denovo.html). 
También necesitamos las lecturas que alineamos al genoma en nuestra práctica anterior.

Para empezar, creemos un directorio para trabajar con nuestros datos:

~~~ {.bash}
$ mkdir STRINGTIE
$ cd STRINGTIE/
~~~ 

Reconstruiremos los transcritos para esta muestra usando Stringtie.

Primero tenemos que ordenar los archivos bam usando samtools - haz esto para todos los archivos:

~~~ {.bash}
samtools sort ../Sp_ds.bam -o Sp_ds_sorted.bam
~~~

Y después ejecutamos stringtie, para ensamblar un transcriptoma guiado para cada una de nuestras muestras. 

~~~ {.bash}
$ stringtie Sp_ds_sorted.bam -o Sp_ds_stringtie.gtf 
~~~ 

Repite este comando para cada una de las condiciones restantes hs, log y plat. 

Entonces estamos listos para mezclar los transcriptomas usando `stringtie merge`:

~~~ {.bash}
stringtie --merge Sp*gtf -o Sp_stringtie_merged.gtf
~~~ 

El grupo de transcritos combinados (merged) esta ahora en el archivo `Sp_stringtie_merged.gtf`.
Este archivo es nuestro resultado final. Un archivo gtf es un 
formato de coordenadas genómicas jerárquico que nos permite saber dónde se ubican los genes
en una referencia y su estructura (exones, intrones, CDS, etc). Una descripción más extensa 
de este archivo se encuentra [aquí](https://www.ensembl.org/info/website/upload/gff.html).

Usaremos este archivo "merged" para realizar el análisis de expresión diferencial. 

