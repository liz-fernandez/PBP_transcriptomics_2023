---
layout: page
title: Control de versiones usando Git
subtitle: Explorando el historial
minutes: 25
---
> ## Objetivos de aprendizaje {.objectives}
>
> *   Identificar y utilizar los códigos identificadores de commit generados por Git.
> *   Comparar distintas versiones de archivos rastreados. 
> *   Restaurar versiones previas de archivos.

Si queremos ver que cambios se han realizado entre distintas versiones podemos utilizar nuevamente el comando `git diff`, pero con la notación `HEAD~1`, `HEAD~2`, etc, para referirnos a commits previos:

~~~ {.bash}
$ git diff HEAD~1 mars.txt
~~~
~~~ {.output}
diff --git a/mars.txt b/mars.txt
index 315bf3a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1,2 +1,3 @@
 Cold and dry, but everything is my favorite color
 The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
~~~
~~~ {.bash}
$ git diff HEAD~2 mars.txt
~~~
~~~ {.output}
diff --git a/mars.txt b/mars.txt
index df0654a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1 +1,3 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
~~~

De este modo podemos construir una cadena de commits. Al eslabón más reciente de la cadena se le conoce como `HEAD`. Podemos nombrar commits previos utilizando la 
notación `~`, entonces `HEAD~1` ("head minus one") significa "el commit previo" mientras que `HEAD~123` regresa 123 commits antes del commit actual. 

También podemos llamar estos commits utilizando las largas series de dígitos y letras que despliega el comando `git log`. 
Estos son identificadores únicos para cada uno de los cambios, y en este caso son realmente únicos:
cada cambio en cualquier archivo en el repositorio 
en cualquier computadora tiene un identificador único 
de 40 caracteres.
Por ejemplo, a nuestro primer commit se le asigno el 
identificador
f22b25e3233b4645dabd0d81e651fe074bd8e73b. 
Intentemos lo siguiente:

~~~ {.bash}
$ git diff f22b25e3233b4645dabd0d81e651fe074bd8e73b mars.txt
~~~
~~~ {.output}
diff --git a/mars.txt b/mars.txt
index df0654a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1 +1,3 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
~~~

Esa es la respuesta correcta,That's the right answer,
Pero teclear cadenas de 40 caracteres aleatorios es ineficiente, 
por lo tanto Git nos permite utilizar solo los primeros caracteres:

~~~ {.bash}
$ git diff f22b25e mars.txt
~~~
~~~ {.output}
diff --git a/mars.txt b/mars.txt
index df0654a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1 +1,3 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
~~~


¡Muy bien! Podemos guardar los cambios a los archivos 
y ver que ha cambiado. Ahora, ¿cómo podemos restaurar 
versiones previas de archivos? Supongamos que sobre escribimos accidentalmente nuestro archivo:

~~~ {.bash}
$ nano mars.txt
$ cat mars.txt
~~~
~~~ {.output}
We will need to manufacture our own oxygen
~~~

`git status` nos dice que el archivo ha cambiado,
pero estos cambios no han sido agregados al área de stagin:

~~~ {.bash}
$ git status
~~~
~~~ {.output}
# On branch master
# Changes not staged for commit:
#   (use "git add <file>..." to update what will be committed)
#   (use "git checkout -- <file>..." to discard changes in working directory)
#
#	modified:   mars.txt
#
no changes added to commit (use "git add" and/or "git commit -a")
~~~

Podemos regresar las cosas a su estado previo 
utilizando el comando `git checkout`:

~~~ {.bash}
$ git checkout HEAD mars.txt
$ cat mars.txt
~~~
~~~ {.output}
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
~~~

Como su nombre sugiere,
`git checkout` checa (restaura) una version previa de un archivo.
En este caso,
le estamos diciendo a Git que queremos recuperar
la versión del archivo guardado en `HEAD`,
que es la última versión o commit guardado. 
Si queremos regresar a versiones anteriores, 
podemos usar el identificador del commit en lugar de `HEAD`:

~~~ {.bash}
$ git checkout f22b25e mars.txt
~~~

> ## No pierdas la cabeza (HEAD) {.callout}
> En la parte superior utilizamos
>
> ~~~ {.bash}
> $ git checkout f22b25e mars.txt
> ~~~
>
> para revertir mars.txt a su estado después del commit f22b25e.
> Si olvidas especificar `mars.txt` en ese comando, git te mostrará 
> el mensaje "You are in
> 'detached HEAD' state." En este estado no debes realizar ningún cambio.
> Puedes reparar esto re-uniendo el `HEAD` usando el comando ``git checkout master``

Es importante recordar que debemos usar el número 
de commit que identifica el estado del repositorio 
*antes* del cambio que estamos tratando de revertir. 
Un error común es utilizar el número del commit 
en el cual hicimos el cambio que estamos tratando de 
deshacer. 
En el ejemplo siguiente, tratamos de revertir el 
repositorio a su estado antes del último commit, 
que es el commit `f22b25e`:

![Git Checkout](fig/git-checkout.svg)

En resumen:

> ## Caricatura acerca de como funciona Git {.callout}
> ![http://figshare.com/articles/How_Git_works_a_cartoon/1328266](fig/git_staging.svg)

> ## Simplificando el caso común {.callout}
>
> Si lees cuidadosamente el resultado del comando `git status`,
> verás que incluye la siguiente pista:
>
> ~~~ {.bash}
> (use "git checkout -- <file>..." to discard changes in working directory)
> ~~~
>
> Como dice, As it says,
> `git checkout` sin un identificador de una versión restaura archivos al estado guardado en el `HEAD`.
> El guión doble `--` se necesita para separar los nombres de los archivos que estan siendo restaurados
> del comando mismo:
> sin estos,
> Git tratará de utilizar el nombre del archivo como el número identificador del commit.

El hecho de que los archivos pueden ser restaurados uno por uno
tiende a cambiar la manera en la que la gente organiza su trabajo. 
Si todo está en un solo documento es difícil (aunque no imposible) 
revertir cambios a una parte del documento (por ejemplo la introducción)
sin deshacer cambios en otras parte (ej: la discusión). 
Por el contrario, si la introducción y la conclusión se almacenan 
en archivos independientes cambiar a distintas versiones de las distintas
partes de forma independiente se vuelve mucho más sencillo. 


> ## Recuperando Versiones previas de un Archivo {.challenge}
>
> Jennifer ha hecho cambios a un script de Python en el cual ha estado trabajado por 
> semanas, y las modificaciones que realizó esta mañana lo "rompieron" por lo que el 
> script ya no funciona. Se ha pasado más de una hora tratando de arreglarlo sin lograrlo ... 
>
> Por suerte, ¡Jennifer ha utilizado control de versiones de su proyecto usando Git!
> ¿Cuáles de los comandos de abajo le permitirán recuperar la última versión confirmada 
> (commited) de su script de Python llamado `data_cruncher.py`?
>
> 1. 
>
>     ~~~
>     $ git checkout HEAD
>     ~~~
> 2. 
>
>     ~~~
>     $ git checkout HEAD data_cruncher.py
>     ~~~
> 3. 
>
>     ~~~
>     $ git checkout HEAD~1 data_cruncher.py
>     ~~~
> 4. 
>
>     ~~~
>     $ git checkout <unique ID of last commit> data_cruncher.py
>     ~~~
> 5. Both 2 & 4


> ## Entendiendo el historial y el flujo de datos (workflow) {.challenge}
>
> ¿Cuál es el resultado del comando `cat venus.txt` al final de esta serie de comandos?
>
> ~~~ {.bash}
> $ cd planets
> $ nano venus.txt #input the following text: Venus is beautiful and full of love
> $ git add venus.txt
> $ nano venus.txt #add the following text: Venus is too hot to be suitable as a base
> $ git commit -m "comments on Venus as an unsuitable base"
> $ git checkout HEAD venus.txt
> $ cat venus.txt #this will print the contents of venus.txt to the screen
> ~~~
>
> 1. 
> 
>     ~~~ {.output}
>     Venus is too hot to be suitable as a base
>     ~~~
>
> 2. 
> 
>     ~~~ {.output}
>     Venus is beautiful and full of love
>     ~~~
>
> 3. 
> 
>     ~~~ {.output}
>     Venus is beautiful and full of love
>     Venus is too hot to be suitable as a base
>     ~~~
>
> 4. 
> 
>     ~~~ {.output}
>     Error dado que cambiaste venus.txt sin confirmar los cambios
>     ~~~

