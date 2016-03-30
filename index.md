---
layout: page
title: Version Control with Git
---

Estamos realizando un estudio transcriptomico para entender la expresion de genes
en la levadura (Saccharomyces cerevisiae) cuando se expresa una 
mutante de la proteina Puf2 que cambia la especificidad 

http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73227

Hemos logrado aislar los transcritos del cromosoma I en dos condiciones: mutante y 
wild type (WT).

Queremos realizar un analisis completo de la expresión diferencial de estos transcritos. 
Durante esta practica realizaremos el análisis desde identificar la calidad de las secuencias
hasta medir si los cambios en expresión son significativos usando métodos estadísticos.

Diagrama de flujo de trabajo para analisis de datos de secuenciacion masiva

Esta mañana recibimos un correo de nuestro colaborador:

> Good day, 
> 
> I have great news. All samples passed the QC and where successfully sequenced in 
> our Hi-Seq 2000. You can download the data from the URL below. The md5sum files are also 
> in the tar-ball.  
> 
> XXxxxxx
> 
> Can’t wait to see if there are indeed major changes to the transcriptome in the mutant! 
> 
> Happy data wranggling!
> 
> Professor X. 

Bajamos los datos de esta dirección y nos preparamos para llevar a cabo el análisis.

> ## Conocimientos previos {.prereq}
>
> Se espera que los estudiantes tengan alguna experiencia con la terminal,
> asi con el uso del lenguaje de programación R. 

## Temas

1.  [Control de calidad de datos de secuenciación masiva](01-quality.html)
2.  [Ensamblando datos de novo](02-assembly.html)
3.  [Mapeo](03-mapping.html)
4.  [Anotación](04-annotation.html)
5.  [Explorador de genomas (genome browser)](05-genomebrowser.html)
6.  [Analisis de expresión diferencial](06-diffexpression.html)
7.  [Otras herramientas](07-others.html)

## Otros recursos

*   [Analisis de microarreglos - Inglés](microarrays.html)
*   [Referencias](reference.html)
