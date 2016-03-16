---
layout: page
title: Control de versiones usando Git
subtitle: Repositorios remotos en GitHub
minutes: 30
---
> ## Objetivos de aprendizaje {.objectives}
>
> *   Explicar que son los repositorios remotos y porque son útiles.
> *   Subir o bajar (push o pull) de un repositorio remoto.

El sistema de control de versiones realmente brilla
cuando comenzamos a colaborar con otras personas. 
Ya tenemos la mayoría de la maquinaria que requerimos para hacer esto;
la única cosa que nos falta el copiar cambios de un repositorio a otro. 

Sistemas como Git nos permiten mover trabajo entre dos repositorios. 
Sin embargo en la práctica es más sencillo utilizar uno de estos 
como un hub central, y mantenerlo en la red en vez de en la computadora
de alguien. La mayoría de los programadores utilizan servicios de almacenamiento como
[GitHub](http://github.com),
[BitBucket](http://bitbucket.org) or
[GitLab](http://gitlab.com/)
para guardar las copias maestras.
Exploraremos los pros y contras de este método en la sección final de esta lección.

Comencemos por compartir los cambios que hemos realizado a nuestro proyecto con el 
mundo. 
Inicia sesión en GitHub,
y da click en el ícono en la parte superior derecha para crear un nuevo repositorio 
llamado `planets`:

![Creating a Repository on GitHub (Step 1)](fig/github-create-repo-01.png)

Llama tu repositorio "planets" y da click en en el botón "Create Repository":

![Creating a Repository on GitHub (Step 2)](fig/github-create-repo-02.png)

Tan pronto se haya creado el repositorio, 
GitHub despliega una página con la dirección URL y alguna información acerca de como
configurar tu repositorio local:

![Creating a Repository on GitHub (Step 3)](fig/github-create-repo-03.png)

Esto ejecuta los siguientes comandos en los servidores de GitHub:

~~~ {.bash}
$ mkdir planets
$ cd planets
$ git init
~~~

Nuestro repositorio local aún contiene trabajo previo realizado en 
en el archivo `mars.txt`, pero el repositorio remoto en GitHub aún
no tiene ningún archivo:

![Freshly-Made GitHub Repository](fig/git-freshly-made-github-repo.svg)

El siguiente paso es conectar ambos repositorios.
Para esto convertimos el repositorio de GitHub en un [remote](reference.html#remote)
del repositorio local. 
La página de entrada del repositorio en GitHub incluye
el identificador que requerimos para identificarlo:

![Where to Find Repository URL on GitHub](fig/github-find-repo-string.png)

Da click en el link 'HTTPS' para cambiar el [protocol](reference.html#protocol) de SSH a HTTPS.

> ## HTTPS vs SSH {.callout}
>
> Utilizamos HTTPS ya que no requiere configuración adicional. 
> Después del curso puedes configurar acceso via SSH, el cual es un poco más seguro. 
> Sigue alguno de estos excelentes tutoriales:
> [GitHub](https://help.github.com/articles/generating-ssh-keys),
> [Atlassian/BitBucket](https://confluence.atlassian.com/display/BITBUCKET/Set+up+SSH+for+Git)
> y [GitLab](https://about.gitlab.com/2014/03/04/add-ssh-key-screencast/)
> (este último incluye videos).

![Changing the Repository URL on GitHub](fig/github-change-repo-string.png)

Copia esa dirección URL del navegador, 
ve al repositorio local `planets` y 
ejecuta este comando:

~~~ {.bash}
$ git remote add origin https://github.com/vlad/planets.git
~~~

Segurarate de utilizar la URL para tu repositorio en vez de la de Vlad:
debes utilizar tu nombre de usuario en vez de `vlad`.

Podemos verificar que el comando funcionó ejecutando `git remote -v`:

~~~ {.bash}
$ git remote -v
~~~
~~~ {.output}
origin   https://github.com/vlad/planets.git (push)
origin   https://github.com/vlad/planets.git (fetch)
~~~

El nombre `origin` es un apodo local para tu repositorio remoto:
podríamos utilizar otro nombre si lo desearamos, 
pero `origin` es el nombre más común. 

Una vez que el apodo `origin` se ha asignado, 
este comando subirá (push) los cambios realizados en nuestro repositorio local
al repositorio en GitHub:

~~~ {.bash}
$ git push origin master
~~~
~~~ {.output}
Counting objects: 9, done.
Delta compression using up to 4 threads.
Compressing objects: 100% (6/6), done.
Writing objects: 100% (9/9), 821 bytes, done.
Total 9 (delta 2), reused 0 (delta 0)
To https://github.com/vlad/planets
 * [new branch]      master -> master
Branch master set up to track remote branch master from origin.
~~~

> ## Proxy {.callout}
>
> Si la red a la que estás conectado(a) utiliza un proxy existe la posibilidad de que el
> comando anterior haya fallado y mostrado el mensaje de error "Could not resolve hostname". 
> Para resolver este problema necesitas informarle a Git acerca del proxy:
>
> ~~~ {.bash}
> $ git config --global http.proxy http://user:password@proxy.url
> $ git config --global https.proxy http://user:password@proxy.url
> ~~~
>
> Cuando te conectas a otras redes que no utilizan un proxy deberás
> decirle a Git que desactive el proxy utilizando:
>
> ~~~ {.bash}
> $ git config --global --unset http.proxy
> $ git config --global --unset https.proxy
> ~~~

> ## Administradores de contraseñas {.callout}
>
> If your operating system has a password manager configured, `git push` will
> try to use it when it needs your username and password. If you want to type
> your username and password at the terminal instead of using
> a password manager, type:
>
> ~~~ {.bash}
> $ unset SSH_ASKPASS
> ~~~
>
> You may want to add this command at the end of your `~/.bashrc` to make it the
> default behavior.

Our local and remote repositories are now in this state:

![GitHub Repository After First Push](fig/github-repo-after-first-push.svg)

> ## The '-u' Flag {.callout}
>
> You may see a `-u` option used with `git push` in some documentation.
> It is related to concepts we cover in our intermediate lesson,
> and can safely be ignored for now.

We can pull changes from the remote repository to the local one as well:

~~~ {.bash}
$ git pull origin master
~~~
~~~ {.output}
From https://github.com/vlad/planets
 * branch            master     -> FETCH_HEAD
Already up-to-date.
~~~

Pulling has no effect in this case
because the two repositories are already synchronized.
If someone else had pushed some changes to the repository on GitHub,
though,
this command would download them to our local repository.

> ## GitHub GUI {.challenge}
> 
> Browse to your `planets` repository on GitHub.
> Under the Code tab, find and click on the text that says "XX commits" (where "XX" is some number). 
> Hover over, and click on, the three buttons to the right of each commit.
> What information can you gather/explore from these buttons?
> How would you get that same information in the shell?

> ## GitHub Timestamp {.challenge}
>
> Create a remote repository on GitHub.
> Push the contents of your local repository to the remote.
> Make changes to your local repository and push these changes.
> Go to the repo you just created on Github and check the [timestamps](reference.html#timestamp) of the files.
> How does GitHub record times, and why?

> ## Push vs. commit {.challenge}
>
> In this lesson, we introduced the "git push" command.
> How is "git push" different from "git commit"?

> ## Fixing up remote settings {.challenge}
>
> It happens quite often in practice that you made a typo in the
> remote URL. This exercice is about how to fix this kind of issues.
> First start by adding a remote with an invalid URL:
>
> ~~~ {.bash}
> git remote add broken https://github.com/this/url/is/invalid
> ~~~
>
> Do you get an error when adding the remote? Can you think of a
> command that would make it obvious that your remote URL was not
> valid? Can you figure out how to fix the URL (tip: use `git remote
> -h`)? Don't forget to clean up and remove this remote once you are
> done with this exercise.
