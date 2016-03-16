---
layout: page
title: Control de versiones usando Git
subtitle: Control Automatizado de Versiones
minutes: 5
---
> ## Objetivos de aprendizaje {.objectives}
>
> *   Entender los beneficios de un sistema automatizado de control de versiones. 
> *   Entender las bases del funcionamiento de Git. 

Comenzaremos explorando como podemos utilizar el control de versiones
para mantener un registro de que cambios hizo una sola persona, así como
del cuando se realizaron estos cambios. Estos cambios pueden ser en un escrito
o en una pieza de código. 
Incluso si no realizas trabajos en colaboración con otros personas, 
usar un sistema de control de versiones es mucho mejor que terminar en una 
situación como esta: 

[![Piled Higher and Deeper by Jorge Cham, http://www.phdcomics.com](fig/phd101212s.gif)](http://www.phdcomics.com)

"Piled Higher and Deeper" by Jorge Cham, http://www.phdcomics.com

Todos hemos estado en una situación similar: parece ridículo tener multiples versiones prácticamente idénticas del mismo documento. Algunos procesadores de texto nos ayudan a lidiar con estas situaciones, como por exemplo la función "Resaltar Cambios" (Track Changes) de Microsoft Word o el historial de versiones de Google Docs. 

Los sistemas de control de versiones comienzan con una version base del documento y después guardan los cambios hechos en cada paso. Podemos pensar en estos sistemas como si fuesen una película: si la regresamos al inicio y comenzamos con el documento base podemos reproducir cada cambio y terminar con la versión final del documento. 

![Los cambios se guardan de manera secuencial](fig/play-changes.svg)

Si comenzamos a pensar en los cambios como entidades independientes al documento mismo, podemos empezar a pensar en reproducir diferentes grupos de cambios sobre el documento base para así obtener diferentes versiones del documento. Por ejemplo, dos usuarios pueden hacer cambios independientes basándose en el mismo documento. 

![Se pueden guardar versiones distintas](fig/versions.svg)

Si no hay cambios conflictivos incluso podemos tratar de reproducir dos grupos de cambios sobre el mismo documento base.

![Se pueden unir versiones multiples](fig/merge.svg)

Un sistema de control de versiones es una herramienta que mantiene un registro de cambios de forma automática y nos ayuda a controlar versiones y a unir cambios hechos a un documento. Nos permite decidir que cambios se deben incorporar a la próxima version, lo que se le conoce como un [commit](reference.html#commit). Además, guarda metadatos útiles acerca de estos cambios permanentes (commits). El historial completo de commits de un proyecto en particular así como sus metadatos son lo que constituye un repositorio o [repository](reference.html#repository). Los Repositorios pueden mantenerse sincronizados en varias computadoras facilitando así las colaboraciones entre distintos individuos. 

> ## La larga historia de los sistemas de control de versiones{.callout}
>
> Los sistemas automatizados de control de versiones no son nada nuevo.
> Herramientas como RCS, CVS o Subversion existen desde el principio de los años 80s y son utilizadas por varias compañias grandes.
> Sin embargo, varios de estos se empiezan a considerar sistemas heredados debido a las diversas limitaciones en sus capacidades. 
> A los sistemas más modernos como Git y [Mercurial](http://swcarpentry.github.io/hg-novice/) se les conoce como
> *distribuidos*, ya que no cuentan con un servidor centralizado que almacena el repositorio.
> Estos sistemas modernos también incluyen poderosas herramientas e unión o "merging" que permiten que multiples contribuidores trabajen 
> en los mismos archivos de manera simultanea. 
