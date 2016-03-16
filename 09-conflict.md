---
layout: page
title: Control de versiones usando Git
subtitle: Conflictos
minutes: 15
---
> ## Objetivos de aprendizaje {.objectives}
>
> *   Explicar que son los conflictos y cuando pueden ocurrir. 
> *   Resolver conflictos resultantes de una fusión (merge). 

Tan pronto como varias personas comienzan a trabajar en el mismo proyecto en 
paralelo, es probable que ocurran errores. Esto puede ocurrir aún con una sola persona:
si estás trabajando en una pieza de código en tu laptop y en el server del lab, puedes
hacer cambios diferentes en cada copia. El control de versiones nos permite administrar estos
[conflictos](reference.html#conflicts) por medio de sus herramientas para
[resolver](reference.html#resolve) cambios que superponen o son conflictivos. 

Para entender mejor como resolver conflictos, debemos crear uno. El archivo
`mars.txt` actualmente tiene la siguiente información en cada una de las copias de
el repositorio `planets`:

~~~ {.bash}
$ cat mars.txt
~~~
~~~ {.output}
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
~~~

Añadamos una línea a una de las copias de los compañeros de cada par:

~~~ {.bash}
$ nano mars.txt
$ cat mars.txt
~~~
~~~ {.output}
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
This line added to Wolfman's copy
~~~

y subamos este cambio a GitHub:

~~~ {.bash}
$ git add mars.txt
$ git commit -m "Adding a line in our home copy"
~~~
~~~ {.output}
[master 5ae9631] Adding a line in our home copy
 1 file changed, 1 insertion(+)
~~~
~~~ {.bash}
$ git push origin master
~~~
~~~ {.output}
Counting objects: 5, done.
Delta compression using up to 4 threads.
Compressing objects: 100% (3/3), done.
Writing objects: 100% (3/3), 352 bytes, done.
Total 3 (delta 1), reused 0 (delta 0)
To https://github.com/vlad/planets
   29aba7c..dabb4c8  master -> master
~~~

Ahora, dejemos que el otro compañero 
realice un cambio diferente en su copia
*sin* actualizar desde GitHub:

~~~ {.bash}
$ nano mars.txt
$ cat mars.txt
~~~
~~~ {.output}
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
We added a different line in the other copy
~~~

Podemos confirmar el cambio de manera local:

~~~ {.bash}
$ git add mars.txt
$ git commit -m "Adding a line in my copy"
~~~
~~~ {.output}
[master 07ebc69] Adding a line in my copy
 1 file changed, 1 insertion(+)
~~~

pero Git no nos permitirá subir los cambios a GitHub:

~~~ {.bash}
$ git push origin master
~~~
~~~ {.output}
To https://github.com/vlad/planets.git
 ! [rejected]        master -> master (non-fast-forward)
error: failed to push some refs to 'https://github.com/vlad/planets.git'
hint: Updates were rejected because the tip of your current branch is behind
hint: its remote counterpart. Merge the remote changes (e.g. 'git pull')
hint: before pushing again.
hint: See the 'Note about fast-forwards' in 'git push --help' for details.
~~~

![Conflictos en los cambios](fig/conflict.svg)

Git detecta que los cambios hechos a una copia se superponen a aquellos hechos 
en la otra y nos impide reescribir nuestro trabajo previo.
Lo que tenemos que hacer es bajar (pull) los cambios de GitHub, 
y fusionarlos ([merge](reference.html#merge)) dentro de la copia en la que estamos
trabajando actualmente y subir (push) es última versión. 
Empecemos bajando los cambios:

~~~ {.bash}
$ git pull origin master
~~~
~~~ {.output}
remote: Counting objects: 5, done.        
remote: Compressing objects: 100% (2/2), done.        
remote: Total 3 (delta 1), reused 3 (delta 1)        
Unpacking objects: 100% (3/3), done.
From https://github.com/vlad/planets
 * branch            master     -> FETCH_HEAD
Auto-merging mars.txt
CONFLICT (content): Merge conflict in mars.txt
Automatic merge failed; fix conflicts and then commit the result.
~~~

`git pull` nos indica que hay un conflicto, 
e indica el conflicto en el archivo afectado:

~~~ {.bash}
$ cat mars.txt
~~~
~~~ {.output}
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
<<<<<<< HEAD
We added a different line in the other copy
=======
This line added to Wolfman's copy
>>>>>>> dabb4c8c450e8475aee9b14b4383acc99f42af1d
~~~

Nuestro cambio&mdash;aquel en `HEAD`&mdash;es precedido por `<<<<<<<`.
Git ha insertado `=======` como un separador entre cambios conflictivos
y marcado el fin del contenido descargado de Github con el símbolo `>>>>>>>`.
(La cadena de letras y dígitos después de ese marcador
identifica el commit que acabamos de descargar). 

Ahora depende de nosotros editar este archivo y remover estos marcadores
para reconciliar los cambios. 
Podemos hacer lo que queramos: guardar los cambios hechos en el repositorio local, 
guardar los cambios hechos en el repositorio remoto, escribir algo nuevo para 
reemplazar a ambos o deshacernos del cambio completamente. 
Reemplacemos ambos de modo que el archivo quede así:

~~~ {.bash}
$ cat mars.txt
~~~
~~~ {.output}
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
We removed the conflict on this line
~~~

Para finalizar la fusión, 
añadimos `mars.txt` con los cambios hechos en la fusión (merge)
y los aprobamos (commit):

~~~ {.bash}
$ git add mars.txt
$ git status
~~~
~~~ {.output}
# On branch master
# All conflicts fixed but you are still merging.
#   (use "git commit" to conclude merge)
#
# Changes to be committed:
#
#	modified:   mars.txt
#
~~~
~~~ {.bash}
$ git commit -m "Merging changes from GitHub"
~~~
~~~ {.output}
[master 2abf2b1] Merging changes from GitHub
~~~

Ahora podemos subir nuestros cambios a GitHub:

~~~ {.bash}
$ git push origin master
~~~
~~~ {.output}
Counting objects: 10, done.
Delta compression using up to 4 threads.
Compressing objects: 100% (6/6), done.
Writing objects: 100% (6/6), 697 bytes, done.
Total 6 (delta 2), reused 0 (delta 0)
To https://github.com/vlad/planets.git
   dabb4c8..2abf2b1  master -> master
~~~

Git mantiene un record de que se unió con que, 
de modo que no tengamos que reparar cosas manualmente
cuando el colaborador que hizo el primer cambio baje los cambios (pull):

~~~ {.bash}
$ git pull origin master
~~~
~~~ {.output}
remote: Counting objects: 10, done.        
remote: Compressing objects: 100% (4/4), done.        
remote: Total 6 (delta 2), reused 6 (delta 2)        
Unpacking objects: 100% (6/6), done.
From https://github.com/vlad/planets
 * branch            master     -> FETCH_HEAD
Updating dabb4c8..2abf2b1
Fast-forward
 mars.txt | 2 +-
 1 file changed, 1 insertion(+), 1 deletion(-)
~~~

Obtenemos el siguiente archivo fusionado:

~~~ {.bash}
$ cat mars.txt 
~~~
~~~ {.output}
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
We removed the conflict on this line
~~~

No tenemos que fusionar nuevamente ya que Git sabe que alguien ya lo ha hecho.

La capacidad del control de versiones de fusionar cambios conflictivos
es otra de las razones por la que los usuarios tienden a dividir sus programas y
artículos en archivos multiples en lugar de almacenar todo en 
un solo archivo grande. 
Hay otro beneficio:
cada vez que hay conflictos repetidos en un archivo en particular, 
el sistema de control de versiones dice a sus usuarios
que deben aclarar quién es responsable de que, 
o encontrar una manera diferente de dividir el trabajo.


> ## Resolviendo conflictos creados por tí {.challenge}
>
> Clona el repositorio creado por tu instructor,
> Añádele un nuevo archivo, 
> y modifica un archivo existente (tu instructor te dirá cual). 
> Cuando te lo pida el instructor, 
> baja los cambios del repositorio para crear un conflicto
> y resuélvelo. 

> ## Conflictos en archivos no textuales {.challenge}
>
> ¿Qué hace Git cuando hay un conflicto en una imagen o en otro archivo
> no textual almacenado usando control de versiones?

> ## Una sesión de trabajo común {.challenge}
>
> Te sientas en tu computadora a trabajar en un proyecto compartido 
> rastreado en un repositorio remote de Git. Durante tu sesión de trabajo, haces lo siguiente,
> aunque no en ese orden: 
> 
> - *Realizar cambios* añadiendo el número `100` al archivo `numbers.txt`
> - *Actualizar el repositorio remoto* para que sea igual que el repositorio local
> - *Celebrar* tu éxito con una cerveza
> - *Actualizar el repositorio local* para que sea igual al repositorio remoto
> - *Añadir los cambios* al área de staging para que sean aprobados
> - *Aprobar los cambios* al repositorio local
>
> ¿En que orden se deben realizar estas acciones para minimizar las probabilidades 
> de generar conflictos? Coloca los comandos en el orden correcto en la columna *acción* de la tabla
> de abajo. Cuando tengas el orden correcto, escribe los comandos correspondientes en la
> columna *comandos*. Hemos escrito algunos pasos para ayudarte a empezar. 
>
> |orden|acción . . . . . . . . . . |comando . . . . . . . . . . |
> |-----|---------------------------|----------------------------|
> |1    |                           |                            |
> |2    |                           | `echo 100 >> numbers.txt`  |
> |3    |                           |                            |
> |4    |                           |                            |
> |5    |                           |                            |
> |6    | Celebrate!                | `AFK`                      |
