---
layout: page
title: Control de versiones usando Git
subtitle: Configurando Git
minutes: 5
---
> ## Objetivos de aprendizaje {.objectives}
>
> * Instalar 'git' en dos plataformas: Windows y Linux
> *  Configurar `git` la primera vez que es utilizado en una computadora.
> *  Entender el significado de la bandera de configuración `--global`.

Para descargar e instalar 'git' en Linux ejecuta el siguiente comando en bash:

sudo apt-get install git-all

Existen diferentes versiones de visualizar Git en plaforma Windows, una
de ellas es emulando la consola de Linux, la podrás encontrar en el siguiente link:

https://git-scm.com/download/win 


Cuando utilizamos Git en una nueva computadora por primera vez,
tenemos que configurar algunas cosas. Así es como Mandy configura Git 
en su nueva laptop:

~~~ {.bash}
$ git config --global user.name "Mandy Grim"
$ git config --global user.email "mandy.grim@purohueso.org"
$ git config --global color.ui "auto"
~~~

(Por favor utilicen su propio nombre y correo electrónico en vez del de Mandy.)

También configura su editor de texto favorito, de acuerdo a los comandos en esta tabla:

| Editor             | Comando de configuración                    |
|:-------------------|:-------------------------------------------------|
| nano               | `$ git config --global core.editor "nano -w"`    |
| Text Wrangler      | `$ git config --global core.editor "edit -w"`    |
| Sublime Text (Mac) | `$ git config --global core.editor "subl -n -w"` |
| Sublime Text (Win, 32-bit install) | `$ git config --global core.editor "'c:/program files (x86)/sublime text 3/sublime_text.exe' -w"` |
| Sublime Text (Win, 64-bit install) | `$ git config --global core.editor "'c:/program files/sublime text 3/sublime_text.exe' -w"` |
| Notepad++ (Win)    | `$ git config --global core.editor "'c:/program files (x86)/Notepad++/notepad++.exe' -multiInst -notabbar -nosession -noPlugin"`|
| Kate (Linux)       | `$ git config --global core.editor "kate"`       |
| Gedit (Linux)      | `$ git config --global core.editor "gedit -s -w"`   |
| emacs              | `$ git config --global core.editor "emacs"`   |
| vim                | `$ git config --global core.editor "vim"`   |


Los comandos de Git se escriben `git verb`, donde `verb` es lo que queremos que git haga.
En este caso, le estamos diciendo a Git:

*   nuestro nombre y dirección de correo electrónico,
*   que coloree el texto de salida (output),
*   cual es nuestro editor de texto preferido,
*   y que queremos utilizar esta configuración globalmente (i.e. en cada proyecto)

Los cuatro comandos que escribimos anteriormente solo se tienen que ejecutar una sola
vez: la bandera `--global` le indica a Git que debe utilizar esta configuración en 
cada proyecto, con ese usuario, en esa computadora. 

Puedes verificar tu configuración en cualquier momento usando:

~~~ {.bash}
$ git config --list
~~~

Puedes cambiar la configuración tantas veces como lo desees: solo utiliza los
mismos comandos para elegir otro editor de texto o para actualizar tu dirección 
de correo electrónico.

> ## Proxy {.callout}
>
> En algunas redes necesitaras usar un
> [proxy](https://en.wikipedia.org/wiki/Proxy_server). De ser necesario también
> tendrás que informarle a Git acerca del proxy:
>
> ~~~ {.bash}
> $ git config --global http.proxy proxy-url
> $ git config --global https.proxy proxy-url
> ~~~
>
> Para desactivar el proxy utiliza:
>
> ~~~ {.bash}
> $ git config --global --unset http.proxy
> $ git config --global --unset https.proxy
> ~~~
