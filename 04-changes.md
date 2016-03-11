---
layout: page
title: Control de versiones usando Git
subtitle: Seguimiento de los cambios
minutes: 20
---
> ## Objetivos de aprendizaje {.objectives}
> 
> *   Realizar el ciclo modificar-añadir-confirmar (modify-add-commit) para un solo archivo.
> *   Explicar donde se almacena la información en cada etapa del ciclo de Git.

Creemos un archivo llamado `mars.txt` que contenga algunas notas
acerca de las posibilidades de utilizar el planeta rojo como base espacial.  
(Para el ejercicio utilizaremos el editor `nano` para editar el archivo;
puedes utilizar tu editor de texto preferido si así lo deseas. Mantén en mente 
que este no tiene que ser el `core.editor` que asignamos como editor global 
previamente.)

~~~ {.bash}
$ nano mars.txt
~~~

Escribe el siguiente texto dentro del archivo `mars.txt`:

~~~ {.output}
Cold and dry, but everything is my favorite color
~~~

`mars.txt` ahora contiene una sola línea, la cual podemos ver ejecutando:

~~~ {.bash}
$ ls
~~~
~~~ {.output}
mars.txt
~~~
~~~ {.bash}
$ cat mars.txt
~~~
~~~ {.output}
Cold and dry, but everything is my favorite color
~~~

Si revisamos de nuevo el estado del proyecto,
Git nos dice que ha encontrado el nuevo archivo:

~~~ {.bash}
$ git status
~~~
~~~ {.output}
# On branch master
#
# Initial commit
#
# Untracked files:
#   (use "git add <file>..." to include in what will be committed)
#
#	mars.txt
nothing added to commit but untracked files present (use "git add" to track)
~~~

El mensaje notificándonos de archivos no rastreados "untracked files" nos indica
que hay un archivo que Git no está rastreando:
Le podemos decir a Git que rastree este archivo utilizando `git add`:

~~~ {.bash}
$ git add mars.txt
~~~

y revisar que el cambio se efectuó correctamente:

~~~ {.bash}
$ git status
~~~
~~~ {.output}
# On branch master
#
# Initial commit
#
# Changes to be committed:
#   (use "git rm --cached <file>..." to unstage)
#
#	new file:   mars.txt
#
~~~

Git ahora sabe que debe rastrear cambios en el archivo `mars.txt`,
pero no ha registrado estos cambios como confirmados en un "commit". 
Para lograr que lo haga, tenemos que ejecutar un comando más:

~~~ {.bash}
$ git commit -m "Start notes on Mars as a base"
~~~
~~~ {.output}
[master (root-commit) f22b25e] Start notes on Mars as a base
 1 file changed, 1 insertion(+)
 create mode 100644 mars.txt
~~~

Cuando ejecutamos `git commit`,
Git toma todo lo que le hemos indicado agregar usando `git add`
y almacena una copia permanentemente dentro del directorio especial `.git`.
Esta copia permanente se llama un [commit](reference.html#commit)
(o [revision](reference.html#revision)) y su código de identificación corto 
es `f22b25e`. (Tu commit puede tener otro identificador). 

Utilizamos la bandera `-m` (de "mensaje")
para guardar un mensaje específico, corto y descriptivo que nos permita
recordar que cambios hicimos y porque. 
Si solo ejecutamos el comando `git commit` sin la opción `-m`,  
Git iniciará `nano` (o cualquier otro editor que hayamos asignado como `core.editor`)
para permitirnos escribir un mensaje más largo. 

Los [mensajes buenos en un commit][commit-messages] comienzan con un pequeño resumen 
(menos de 50 caracteres) de los cambios hechos como parte del commit. 
Si deseas agregar detalles más descriptivos puedes añadir una línea vacía
entre la línea de resumen y otros notas adicionales. 

Si después de esto ejecutamos `git status`:

~~~ {.bash}
$ git status
~~~
~~~ {.output}
# On branch master
nothing to commit, working directory clean
~~~

Git nos informa que todo está actualizado. 
Si queremos saber que cambios se han hecho recientemente,
podemos pedirle a Git que nos muestre el historial reciente
del proyecto utilizando el comando `git log`:

~~~ {.bash}
$ git log
~~~
~~~ {.output}
commit f22b25e3233b4645dabd0d81e651fe074bd8e73b
Author: Vlad Dracula <vlad@tran.sylvan.ia>
Date:   Thu Aug 22 09:51:46 2013 -0400

    Start notes on Mars as a base
~~~

`git log` lista todos los commits realizados en ese repositorio en ordern cronológico
inverso. El listado de cada commit incluye el identificador (el cual comienza
con los mismo caracteres que el identificador corto impreso por el comando 
`git commit` que ejecutamos previamente), el autor del commit, cuando fue creado, 
y el mensaje de registro (log) que le otorgamos a Git cuando creamos el commit. 

> ## ¿Dónde están mis cambios? {.callout}
>
> Si ejecutamos `ls` en este momento, aún veremos un solo archivo llamado `mars.txt`.
> Esto se debe a que Git guarda información acerca del historial de cada archivo 
> en el directorio especial llamado `.git` de modo que nuestro espacio de 
> trabajo o filesystem no se llene con otros archivos (así como para prevenir
> que accidentalmente editemos o borremos una versión previa). 

Ahora supongamos que Mandy quiere añadir más información al archivo. 
(Lo editaremos de nuevo usando `nano` seguido por el comando `cat` para mostrar
el contenido del archivo. Puedes utilizar otro editor y alternativas a `cat`
si lo prefieres). 

~~~ {.bash}
$ nano mars.txt
$ cat mars.txt
~~~
~~~ {.output}
Cold and dry, but everything is my favorite color
The two moons may be a problem for Mandy
~~~

Ahora, si ejecutamos `git status`
nos informa que un archivo que ya conoce ha sido modificado:

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

La última linea proporciona la frase clave:
"no changes added to commit".
Hemos modificado este archivo, 
pero no le hemos dicho a Git que queremos agregar estos cambios
(lo cual hacemos con el comando `git add`)
como tampoco los hemos guardado permanentemente (lo que hacemos usando `git commit`).
Hagámoslo ahora. Es un buen hábito siempre revisar nuestros cambios
antes de guardarlos. Hacemos esto utilizando `git diff`. 
Este comando nos muestra las diferencias entre el estado actual del
archivo y la última versión guardada:


~~~ {.bash}
$ git diff
~~~
~~~ {.output}
diff --git a/mars.txt b/mars.txt
index df0654a..315bf3a 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1 +1,2 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Mandy
~~~

El texto de salida es críptico ya que
es una serie de comandos para herramientas tales como editores y `patch`
indicándoles como reconstruir un archivo dado los cambios en el otro. 
Si lo dividimos en partes:

1.  La primera linea nos indica que Git esta imprimiendo un resultado similar al del
	commando `diff` al comparar la versión antigua y nueva del archivo.
2.  La segunda línea nos dice exactamente que versiones del archivo está comparando
    Git;
    `df0654a` y `315bf3a` son etiquetas únicas generadas por la computadora para 
    cada una de esas versiones.
3.  La tercera y cuarta línea muestran el nombre del archivo que cambiamos.
4.  Las líneas restantes son las más interesantes, ya que nos muestras las 
    diferencias en sí así como las líneas en las que ocurrieron estos cambios. 
    En particular, los símbolos de `+` en la primera columna muestras donde se
    han añadido líneas.

Después de revisar estos cambios ha llegado el momento de aprobarlos a través de un commit:

~~~ {.bash}
$ git commit -m "Add concerns about effects of Mars' moons on Mandy"
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

Whoops:
Git no acepta el commit ya que no utilizamos `git add` primero.
Arreglémoslo:

~~~ {.bash}
$ git add mars.txt
$ git commit -m "Add concerns about effects of Mars' moons on Mandy"
~~~
~~~ {.output}
[master 34961b1] Add concerns about effects of Mars' moons on Mandy
 1 file changed, 1 insertion(+)
~~~

Git insiste en que agreguemos archivos que hemos cambiado antes
de aprobarlos usando commit ya que puede ser que no queramos aprobar
todos los cambios al mismo tiempo.
Por ejemplo, 
supongamos que agregamos algunas referencias a nuestro trabajo de tesis. 
Tal vez queramos aprobar esos cambios y añadir estas a la bibliografía, 
pero *sin* aprobar (commit) los cambios que realizamos a la conclusión 
(los cuales aún no hemos terminado). 

Para permitir esto, 
Git tiene un área especial conocida como *staging area*
en la cual rastrea cambios que han sido añadidos al grupo de cambios o 
[change set](reference.html#change-set) actual, pero que aún no han sido 
aprobados por medio de un commit.

> ## Staging area {.callout}
> 
> Si pensamos en Git como una cámera que toma fotos de los cambios a lo largo
> de la vida de un proyecto, entonces `git add` indica *que* se va a incluir
> en la fotografía y `git commit` esta encargado de *tomar* la fotografía, así como de 
> crear un registro permanente de la misma (como un commit). 
> Si no hay nada en fase "staged" cuando ejecutamos `git commit`, 
> Git te pedirá que utilices `git commit -a` o `git commit --all`,
> lo cual sería el equivalente de !invitar a *todos* los cambios a salir en la foto!
> Sin embargo, casi siempre es mejor añadir cada cambio de manera explícita
> al área de staging, ya que puede ser que apruebes cambios que hiciste pero olvidaste.
> (Siguiendo con la analogía de las fotografías puede que incluyas a un 
> extra con maquillaje incompleto caminando en medio de la foto por 
> haber utilizado `-a`!). 
> Trata de añadir cosas manualmente o frecuentemente te encontrarás buscando 
> el comando "git undo commit". 


![El "Staging Area" de Git](fig/git-staging-area.svg)

Observemos como los cambios a un archivo se mueven de 
nuestro editor de texto a el área de staging y dentro del 
registro de cambios a largo plazo. 
Primero agreguemos otra línea al archivo:

~~~ {.bash}
$ nano mars.txt
$ cat mars.txt
~~~
~~~ {.output}
Cold and dry, but everything is my favorite color
The two moons may be a problem for Mandy
But Billy will appreciate the lack of humidity
~~~
~~~ {.bash}
$ git diff
~~~
~~~ {.output}
diff --git a/mars.txt b/mars.txt
index 315bf3a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1,2 +1,3 @@
 Cold and dry, but everything is my favorite color
 The two moons may be a problem for Mandy
+But Billy will appreciate the lack of humidity
~~~

Hasta ahora todo va bien:
añadimos una línea al final del archivo
(mostrada con un `+` en la primera columna). 
Ahora, pongamos este cambio en el área de staging 
y veamos que reporta `git diff`:

~~~ {.bash}
$ git add mars.txt
$ git diff
~~~

No hay resultado: 
en lo que a Git concierne 
no hay diferencia entre lo que se le esta pidiendo que guarde de forma
permanente y lo que esta actualmente en el directorio.
Sin embargo, 
si hacemos lo siguiente: 

~~~ {.bash}
$ git diff --staged
~~~
~~~ {.output}
diff --git a/mars.txt b/mars.txt
index 315bf3a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1,2 +1,3 @@
 Cold and dry, but everything is my favorite color
 The two moons may be a problem for Mandy
+But Billy will appreciate the lack of humidity
~~~

nos muestra un cambio entre 
el último cambio aprobado (commited) 
y el archivo que se encuentra en el área de staging. 
Guardemos nuestros cambios:

~~~ {.bash}
$ git commit -m "Discuss concerns about Mars' climate for Billy"
~~~
~~~ {.output}
[master 005937f] Discuss concerns about Mars' climate for Billy
 1 file changed, 1 insertion(+)
~~~

y revisemos una vez más el status:

~~~ {.bash}
$ git status
~~~
~~~ {.output}
# On branch master
nothing to commit, working directory clean
~~~

y revisemos el historial para ver que hemos hecho hasta el momento:

~~~ {.bash}
$ git log
~~~
~~~ {.output}
commit 005937fbe2a98fb83f0ade869025dc2636b4dad5
Author: Vlad Dracula <vlad@tran.sylvan.ia>
Date:   Thu Aug 22 10:14:07 2013 -0400

    Discuss concerns about Mars' climate for Billy

commit 34961b159c27df3b475cfe4415d94a6d1fcd064d
Author: Vlad Dracula <vlad@tran.sylvan.ia>
Date:   Thu Aug 22 10:07:21 2013 -0400

    Add concerns about effects of Mars' moons on Mandy

commit f22b25e3233b4645dabd0d81e651fe074bd8e73b
Author: Vlad Dracula <vlad@tran.sylvan.ia>
Date:   Thu Aug 22 09:51:46 2013 -0400

    Start notes on Mars as a base
~~~

Para recapitular, cada vez que queremos añadir cambios a nuestro repositorio, 
tenemos primero que añadir los archivos modificados al "staging area" 
(`git add`) y subsecuentemente aprobar los cambios para que sean agregados al 
repositorio (`git commit`):


![El Flujo de trabajo de Git Commit](fig/git-committing.svg)

> ## Eligiendo un mensaje para el commit {.challenge}
>
> Cual de los siguiente mensajes sería más apropiado para el último commit realizado en 
> el archivo `mars.txt`?
> 
> 1. 
>
>     "Cambios"
> 2. 
>
>     "Línea añadida 'But Billy will appreciate the lack of humidity' a mars.txt"
> 3. 
>
>     "Discución acerca de los efectos del clima marciano en Billy"

> ## Enviando commits a Git {.challenge}
>
> ¿Cuál o cuáles de los siguientes comandos guardarían los cambios hechos a `myfile.txt` a mi repositorio de Git local?
>
> 1. 
>
>     ~~~
>     $ git commit -m "my recent changes"
>     ~~~
> 2. 
>
>     ~~~
>     $ git init myfile.txt
>     $ git commit -m "my recent changes"
>     ~~~
> 3. 
>
>     ~~~
>     $ git add myfile.txt
>     $ git commit -m "my recent changes"
>     ~~~
> 4. 
>
>     ~~~
>     $ git commit -m myfile.txt "my recent changes"
>     ~~~

> ## `bio` Repositorio {.challenge}
>
> Crea un nuevo repositorio de Git en tu computadora llamado `bio`.
> Escribe un biografía de tres líneas en un archivo llamado `me.txt`,
> agrega y aprueba tus cambios,
> modifica una línea, agrega una cuarta línea fourth line, y muestra las diferencias
> entre su estado actualizado y su estado original.

> ## Author y Committer {.challenge}
>
> Por cada uno de los commits que haz realizado, Git ha almacenado tu nombre 
> dos veces. Se almacena tu nombre como el autor así como la persona realizando
> los commits o "commiter". Puedes revisar esto pidiéndole a Git que te muestre 
> más información acerca de tus commits más recientes: 
>
> ~~~
> $ git log --format=full
> ~~~
>
> Al realizar un commit puedes nombrar a alguien más como autor:
>
> ~~~
> $ git commit --author="Stanley Kubrick <s.kubrick@sciland.com>"
> ~~~
>
> Crea un nuevo repositorio y crea dos commits: uno sin la opción 
> `--author` y otra nombrando a un compañero tuyo como autor.
> Ejecuta `git log` y `git log --format=full`. Piensa en que manera
> esto puede ayudarte a colaborar con tus compañeros.

[commit-messages]: http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html
