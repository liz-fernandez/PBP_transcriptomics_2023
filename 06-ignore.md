---
layout: page
title: Control de versiones usando Git
subtitle: Ignorando Cosas
minutes: 5
---
> ## Objetivos de aprendizaje {.objectives}
>
> *   Configurar Git para ignorar archivos específicos.
> *   Explicar porque puede ser útil ignorar cosas.

¿Qué sucede si hay archivos en el repositorio que no queremos sean rastreados
por Git? Por ejemplo, archivos de respaldo creados por nuestro editor de texto 
o archivos intermediarios creados durante el análisis de datos. 
Creemos algunos archivos de prueba:


~~~ {.bash}
$ mkdir results
$ touch a.dat b.dat c.dat results/a.out results/b.out
~~~

y veamos que nos dice Git:

~~~ {.bash}
$ git status
~~~
~~~ {.output}
# On branch master
# Untracked files:
#   (use "git add <file>..." to include in what will be committed)
#
#	a.dat
#	b.dat
#	c.dat
#	results/
nothing added to commit but untracked files present (use "git add" to track)
~~~

Poner estos archivos bajo control de versiones sería una pérdida de espacio de disco.
Aún peor, tenerlos listados puede distraernos de cambios que realmente importan. 
Digámosle a Git que los ignore. 

Hacemos esto creando un archivo en el directorio base (o root) de nuestro proyecto
llamado `.gitignore`:

~~~ {.bash}
$ nano .gitignore
$ cat .gitignore
~~~
~~~ {.output}
*.dat
results/
~~~

Estas expresiones le dicen a Git que ignore cualquier archivo cuyo nombre tiene la 
extensión `.dat` así como todo lo que se encuentra en el directorio `results`. 
(Si cualquiera de estos archivos ya estaba siendo rastreado, Git continuará 
registrando sus cambios).  

Una vez creado este archivo, 
el mensaje de salida de `git status` es mucho más limpio:


~~~ {.bash}
$ git status
~~~
~~~ {.output}
# On branch master
# Untracked files:
#   (use "git add <file>..." to include in what will be committed)
#
#	.gitignore
nothing added to commit but untracked files present (use "git add" to track)
~~~

Lo único que Git nota es el recién creado archivo `.gitignore`. 
Podrías pensar que no querríamos rastrearlo, 
pero todos aquellos con los que compartimos el repositorio probablemente 
quieran ignorar las mismas cosas que nosotros estamos ignorando. 
Agreguemos y confirmemos `.gitignore` en el repositorio:

~~~ {.bash}
$ git add .gitignore
$ git commit -m "Add the ignore file"
$ git status
~~~
~~~ {.output}
# On branch master
nothing to commit, working directory clean
~~~

Como un extra, utilizar `.gitignore` nos ayuda a evitar agregar accidentalmente al repositorio archivos que ya habíamos 
decidido ignorar:

~~~ {.bash}
$ git add a.dat
~~~
~~~ {.output}
The following paths are ignored by one of your .gitignore files:
a.dat
Use -f if you really want to add them.
fatal: no files added
~~~

Si realmente deseamos ignorar nuestras configuraciones previas, 
podemos utilizar el comando `git add -f` para forzar a Git a añadir algo. 
También podemos ver el estado de archivos ignorados si así lo deseamos:

~~~ {.bash}
$ git status --ignored
~~~
~~~ {.output}
# On branch master
# Ignored files:
#  (use "git add -f <file>..." to include in what will be committed)
#
#        a.dat
#        b.dat
#        c.dat
#        results/

nothing to commit, working directory clean
~~~

> ## Ignorando archivos anidados {.challenge}
>
> Dado un directorio con la siguiente estructura:
>
> ~~~
> results/data
> results/plots
> ~~~
>
> ¿Cómo ignorarías solo `results/plots` y no `results/data`?

> ## Incluyendo archivos específicos {.challenge}
>
> ¿Como ignorarías todos los archivos con la extensión `.data` en tu directorio base, con la excepción de
> `final.data`?
> Pista: Averigua lo que el operador `!` (signo de exclamación) hace. 

