---
layout: page
title: Control de versiones usando Git
subtitle: Creando un repositorio
minutes: 10
---
> ## Objetivos de aprendizaje {.objectives}
> 
> *   Create a local Git repository.

Una vez que Git esta configurado podemos comenzar a usarlo. 
Creemos un directorio que contenga nuestro trabajo y entremos 
a ese directorio:

~~~ {.bash}
$ mkdir planets
$ cd planets
~~~

Una vez dentro, digámosle a Git que cree un [repository](reference.html#repository)&mdash llamado `planets` ; un lugar en donde Git pueda almacenar versiones de nuestros archivos:

~~~ {.bash}
$ git init
~~~

Si utilizamos el comando `ls` para mostrar el contenido del directorio,
pareciera que nada a cambiado:

~~~ {.bash}
$ ls
~~~

Pero si agregamos la bandera `-a` para mostrarlo todo,
podemos ver que Git genera un directorio oculto dentro de `planets` llamado `.git`:

~~~ {.bash}
$ ls -a
~~~
~~~ {.output}
.	..	.git
~~~

Git almacena información del proyecto en el este subdirectorio especial. 
Si alguna vez lo borramos perderemos todo el historial de ese proyecto. 

Podemos revisar que todo está configurado correctamente pidiéndole a Git 
que nos indique el estado (status) actual del nuestro proyecto:

~~~ {.bash}
$ git status
~~~
~~~ {.output}
# On branch master
#
# Initial commit
#
nothing to commit (create/copy files and use "git add" to track)
~~~

> ## Lugares para crear repositorios de Git {.challenge}
>
> Mandy comienza un nuevo proyecto llamado `moons`, relacionado a su proyecto `planets`.
> A pesar de las advertencias de Billy, ejecuta la siguiente secuencia de comandos 
> para crear un repositorio de Git dentro de otro:
> 
> ~~~ {.bash}
> cd             # return to home directory
> mkdir planets  # make a new directory planets
> cd planets     # go into planets
> git init       # make the planets directory a Git repository
> mkdir moons    # make a sub-directory planets/moons
> cd moons       # go into planets/moons
> git init       # make the moons sub-directory a Git repository
> ~~~
> 
> ¿Porqué es una mala idea hacer esto?
> ¿Cómo puede Mandy "deshacer" su último `git init`?
