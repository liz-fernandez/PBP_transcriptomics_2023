---
layout: page
title: Control de versiones usando Git
---

Billy y Mandy han sido contratados para investigar si es posible enviar un 
robot a la superficie de Marte. Quieren trabajar de manera conjunta en los 
planes, pero han tenido problemas para editar sus notas simultáneamente en 
proyectos previos. Si deciden tomar turnos, cada uno pasará mucho tiempo esperando
a que el otro termine, pero si trabajan simultaneamente en su copia individual 
del reporte y mandan correos para enviarse cambios entre sí siempre hay cambios que
se pierden, se sobre escriben o se duplican.  

Un colega les sugiere utilizar [control de versiones](reference.html#version-control) 
para manejar su trabajo. El control de versiones es mejor que enviar versiones por
correo electrónico ya que:

*   Ningún cambio que haya sido confirmado (commited) al sistema de control de versiones
    se pierde. Dado que todas la versiones anteriores de los archivos están guardadas, 
    siempre es posible regresar en el tiempo y ver exactamente quien escribió que en un 
    día en particular, o que versiones de que programa fueron utilizadas para generar 
    ciertos resultados.

*   Ya que tenemos un registro de quien hizo que cambios y cuando los realizó, 
    sabemos a quien preguntar si tenemos preguntas acerca de estos cambios en 
    el futuro, y, de ser necesario, revertirlos a una versión previa, similar a 
    la función "deshacer" en un editor de texto. 

*   Cuando varias personas colaboran en un mismo proyectos, es probable que ignoren o 
    sobre escriban los cambios de otro colaborador de forma accidental. El sistema 
    de control de versiones automáticamente notifica a los usuarios cada vez que se 
    detecta un conflicto entre los cambios propuestos por dos colaboradores. 

El sistema de control de versiones no es solo útil para trabajos en grupo: 
investigadores individuales se pueden beneficiar inmensamente de sus capacidades. 
Mantener un registro de que ha cambiado, cuando y porque es extremadamente práctico
y útil en casos en los que se requieran regresar a un proyecto después de mucho tiempo
(por ejemplo después de un año o más, cuando la memoria ya no está fresca). 

El sistema de control de versiones es la bitácora de laboratorio del mundo digital:
Es lo que utilizan los profesionales para mantener un registro de lo que han hecho 
así como para colaborar con otras personas. Cada proyecto de desarrollo de 
software a gran escala depende de este sistema, y muchos programadores lo 
utilizan también para proyectos pequeños. Y no solo se utiliza para software:
libros, artículos, sets de datos pequeños y cualquier otra entidad de datos sobre 
la que se realizan cambios a través de tiempo o tiene que ser compartida puede y
debe ser almacenada utilizando un sistema de control de versiones. 

> ## Conocimientos previos {.prereq}
>
> En esta clase aprenderemos a usar Git desde la terminal de Linux/Unix. 
> Se espera que los estudiantes tengan alguna experiencia con la terminal,
> *aunque no es indispensable*.

> ## Preparativos {.getready}
>
> Tener git instalado o accesarlo a traves del cuaderno en linea Jupyter. 
> Si ya lo instalaste: Estas listo!

## Temas

1.  [Control Automatizado de Versiones](01-basics.html)
2.  [Configurando Git](02-setup.html)
3.  [Creando un repositorio](03-create.html)
4.  [Seguimiento de los cambios](04-changes.html)
5.  [Explorando el historial](05-history.html)
6.  [Ignorando Cosas](06-ignore.html)
7.  [Repositorios remotos en GitHub](07-github.html)
8.  [Colaborando](08-collab.html)
9.  [Conflictos](09-conflict.html)


[//]: # (10. [Open Science](10-open.html))
[//]: # (11. [Licencias](11-licensing.html))
[//]: # (12. [Referencias](12-citation.html))
[//]: # (13. [Albergando información](13-hosting.html))

## Otros recursos

*   [Referencias - Inglés](reference.html)
*   [Discusión - Inglés](discussion.html)

