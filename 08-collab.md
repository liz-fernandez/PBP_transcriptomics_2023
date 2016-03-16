---
layout: page
title: Control de versiones usando Git
subtitle: Colaborando
minutes: 25
---
> ## Objetivos de aprendizaje {.objectives}
>
> *  Clona un repositorio remoto.
> *  Colabora subindo cambios a un repositorio común. 

Para el siguiente paso deben trabajar en parejas. 
Una persona será el "Dueño" (esta es la persona cuyo repositorio de GitHub usaremos al inicio de este ejercicio)
y la otra persona será el "Colaborador" (esta es la persona que clonará el repositorio del Dueño y 
realizará cambios al mismo). 

> ## Practicando solo {.callout}
>
> Si estás estudiando esta lección solo(a), puedes continuar abriendo 
> una segunda terminal, y cambiando a otro directorio (e.g. `/tmp`).
> Esta ventana representará a tu pareja, trabajando en otra computadora. 
> En este caso no tendrás que otorgar acceso a alguien más en GitHub, 
> ya que tu serás ambos compañeros. 

El Dueño deberá darle acceso al Colaborador. 
En Github, da click en el botón settings en el lado derecho, 
después selecciona Colaboradores, e inserta el nombre de tu compañero. 

![Añadiendo colaboradores en GitHub](fig/github-add-collaborators.png)

El Colaborador debe trabajar en este proyecto de manera local. El o ella deberá 
cambiar de directorio usando el comando `cd` para entrar un directorio que no contenga
el folder `planets`, y entonces hacer una copia del repositorio del Dueño:

~~~ {.bash}
$ git clone https://github.com/vlad/planets.git
~~~

Remplaza 'vlad' con el nombre del Dueño:

`git clone` crea una copia local fresca en tu repositorio remoto.

![After Creating Clone of Repository](fig/github-collaboration.svg)

El Colaborador puede ahora hacer cambios en su copia del respositorio:

~~~ {.bash}
$ cd planets
$ nano pluto.txt
$ cat pluto.txt
~~~
~~~ {.output}
It is so a planet!
~~~
~~~ {.bash}
$ git add pluto.txt
$ git commit -m "Some notes about Pluto"
~~~
~~~ {.output}
 1 file changed, 1 insertion(+)
 create mode 100644 pluto.txt
~~~

Y subir estos cambios a GitHub:

~~~ {.bash}
$ git push origin master
~~~
~~~ {.output}
Counting objects: 4, done.
Delta compression using up to 4 threads.
Compressing objects: 100% (2/2), done.
Writing objects: 100% (3/3), 306 bytes, done.
Total 3 (delta 0), reused 0 (delta 0)
To https://github.com/vlad/planets.git
   9272da5..29aba7c  master -> master
~~~

Nota que no tuvimos que crear un remoto llamado `origin`:
Git lo crea automaticamente, 
usando este nombre, 
cuando clonamos un repositorio. 
(Es por ellos que `origin` fue una buena elección
de nombre cuando configuramos los remotos manualmente).

Ahora podemos descargar los cambios en el repositorio original en nuestra máquina:

~~~ {.bash}
$ git pull origin master
~~~
~~~ {.output}
remote: Counting objects: 4, done.
remote: Compressing objects: 100% (2/2), done.
remote: Total 3 (delta 0), reused 3 (delta 0)
Unpacking objects: 100% (3/3), done.
From https://github.com/vlad/planets
 * branch            master     -> FETCH_HEAD
Updating 9272da5..29aba7c
Fast-forward
 pluto.txt | 1 +
 1 file changed, 1 insertion(+)
 create mode 100644 pluto.txt
~~~

> ## Revisando los cambios {.challenge}
>
> El Dueño sube (push) commits al repositorio sin informar a su Colaborador. 
> ¿Cómo puede el Colaborador averiguar que ha cambiado en la 
> línea de comandos? ¿Y en GitHub?
> 
> ## Comentando cambios en GitHub {.challenge}
>
> El Colaborados tiene algunas preguntas acerca de un cambio en una línea hecha por el Dueño y
> tiene algunas sugerencias al respecto
> 
> GitHub permite comentar en un `diff` o en un `commit`. Sobre la línea del código  
> a comentar aparece un ícono azul para abrir una ventana de comentarios. c
> 
> El Colaborador publica sus comentarios y sugerencias usando la interfaz de GitHub. 
