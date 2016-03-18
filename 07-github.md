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
> Si tu sistema operativo tiene un administrador de contraseñas configurado, `git push` 
> tratará de utilizarlo cuando necesite tu nombre de usuario y contraseña. Si prefieres teclear tu 
> usuario y contraseña en la terminal en vez del administrador, teclea:
>
> ~~~ {.bash}
> $ unset SSH_ASKPASS
> ~~~
>
> Puedes añadir este comando al final del archivo `~/.bashrc` para hacerlo el 
> default.

Nuestro repositorio local y remoto están ahora en este estado:

![GitHub Repository After First Push](fig/github-repo-after-first-push.svg)

> ## La opción '-u' {.callout}
>
> En ciertos tutoriales encontrarás la opción `-u` usada en conjunto con `git push`.
> Se relaciona a conceptos que cubriremos en nuestra lección intermedia y 
> por el momento podemos ignorarla. 

También podemos bajar cambios del repositorio remoto al local:

~~~ {.bash}
$ git pull origin master
~~~
~~~ {.output}
From https://github.com/vlad/planets
 * branch            master     -> FETCH_HEAD
Already up-to-date.
~~~

Bajar en este caso no tienen ningún efecto, 
ya que ambos repositorios están sincronizados.
Si alguien más hubiese subido cambios al repositorio de GitHub, 
este comando nos permitiría bajarlos a nuestro repositorio local. 

Sin embargo se debe tener cuidado de actualizar el repositorio local si 
si otra  persona se encuentra modificando archivos de forma simultanea, ya que 
al tratar de subir archivos sin antes bajar las últimas modificaciones, el 
usuario caerá en un conflicto. 


> ## GitHub GUI {.challenge}
> 
> Navega hasta tu repositorio `planets` en GitHub.
> Bajo la pestaña "Code", encuentra y da click en el texto de dice "XX commits" (donde "XX" es un número). 
> Pon el mouse sobre el mismo y da click a los tres botones a la derecha de cada commit. 
> ¿Qué información te proporciona cada uno de estos botones?
> ¿Cómo podrías obtener la misma información en la línea de comandos?

> ## Marca del tiempo en GitHub {.challenge}
>
> Crea un repositorio remoto en GitHub.
> Sube (push) los contenidos de tu repositorio local al remoto. 
> Realiza cambios en tu repositorio local y súbelos a GitHub.
> Ve al repositorio que acabase de crear y revisa la marca de tiempo o [timestamp](reference.html#timestamp) del archivo.
> ¿Cómo registra GitHub los tiempos y por qué?

> ## Push vs. commit {.challenge}
>
> En esta lección introdujimos el comando "git push".
> ¿En qué difiere "git push" de "git commit"?

