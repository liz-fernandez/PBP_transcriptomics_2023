---
layout: page
title: Análisis transcriptómicos
subtitle: Control de calidad de datos de secuenciación masiva
minutes: 5
---
> ## Objetivos de aprendizaje {.objectives}
>
> *   Entender el formato FastQ.
> *   Entender el concepto de calidad de una secuencia. 
> *   Aprender a utilizar software para medir la calidad de datos de secuenciación masiva.

Lo primero que hacemos al recibir nuestros datos es descargarlos usando el comando wget:

~~~ {.bash}
$ wget XX
~~~
 
Los archivos de secuenciación generalmente están en un formato comprimido llamado tar. 
Podemos descomprimir los datos usando ese mismo comando:

~~~ {.bash}
$ tar -xvf XX
~~~

La bandera `-x`

Esto genera un directorio llamado `FastQC_Short`. Entramos en ese directorio:

~~~ {.bash}
$ cd FastQC_Short
~~~
 
Revisamos su contenido:

~~~ {.bash}
$ ls
~~~
~~~ {.output}

~~~
 
Revisemos uno de los archivos:
~~~ {.bash}
head X
~~~
~~~ {.output}

~~~


http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/RNA-Seq_fastqc.html
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/small_rna_fastqc.html
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/RRBS_fastqc.html
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/pacbio_srr075104_fastqc.html
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/454_SRR073599_fastqc.html




> ## La larga historia de los sistemas de control de versiones{.callout}
>
> Los sistemas automatizados de control de versiones no son nada nuevo.
> Herramientas como RCS, CVS o Subversion existen desde el principio de los años 80s y son utilizadas por varias compañias grandes.
> Sin embargo, varios de estos se empiezan a considerar sistemas heredados debido a las diversas limitaciones en sus capacidades. 
> A los sistemas más modernos como Git y [Mercurial](http://swcarpentry.github.io/hg-novice/) se les conoce como
> *distribuidos*, ya que no cuentan con un servidor centralizado que almacena el repositorio.
> Estos sistemas modernos también incluyen poderosas herramientas e unión o "merging" que permiten que multiples contribuidores trabajen 
> en los mismos archivos de manera simultánea. 
