# Utilisation de cmake pour la compilation de MIDAS

## Compilation

La version avec cmake minimise le nombre de changements pour que la
compilation et l'installation restent le plus proche possible du système
actuellement utilisé.

Un script (cado) permettant de simplifier l'utilisation de cmake a également
été installé: il est surtout utile si un lien symbolique est utilisé pour le
répertoire de compilation. La façon plus traditionnelle est également proposée.

Les variables utilisées pour configurer les répertoires de build et
d'installation sont ${MIDAS_COMPILE_DIR_MAIN} et ${MIDAS_ABS_LEAFDIR}, tels
que définis dans le fichier de configuration.

L'utilisation est la suivante, surtout si un lien symbolique pour le
répertoire de compilation est utilisé:
```
. ./src/config-cmake.dot.sh
cado cmake # configuration
cado build -j # compilation
cado install # installation (binaires et scripts)
cado package # pour générer le package SSM
```

L'utilisation traditionnelle avec cmake est la suivante :
```
. ./src/config-cmake.dot.sh
mkdir build
cd build
cmake .. # configuration
make -j # compilation
make install # installation (binaires et scripts)
make package # pour générer le package SSM
```

## Changements apportés

Les changements sont les suivants:

- ajout du fichier MANIFEST: informations sur MIDAS et liste des librairies
  et de leur version minimale requise. Cette liste sera analysée par cmake
  pour trouver la bonne version à charger.

- ajout d'un fichier CMakeLists.txt:
  . dans le répertoire principal
  . dans le répertoire scripts/convenient_tools/dumpBmatrix
  . dans le répertoire src/modules
  . dans le répertoire src/programs
  . dans le répertoire tools

  Ces fichiers permettent d'indiquer quoi et comment compiler les sources
  On peut regarder src/modules/CMakeLists.txt, qui montre un exemple simple

- le fichier src/config.dot.sh a été copié dans config-cmake.dot.sh et modifié, notamment pour
  ne plus utiliser makedepf90, ou traduire les options habituellement
  envoyées aux outils comme s.f90. Exemple:
```
-    FOPTMIZ=4
+    FOPTMIZ="-O3 -fast-transcendentals -no-prec-div -fpic -ip -no-prec-sqrt"
```

- renommer les fichiers `*ftn90` (extension non reconnue de cmake) vers `*F90`:
```
  src/modules/clibInterfaces_mod.ftn90 -> src/modules/clibInterfaces_mod.F90
  src/modules/codePrecision_mod.ftn90 -> src/modules/codePrecision_mod.F90
  src/modules/rttovInterfaces_mod.ftn90 -> src/modules/rttovInterfaces_mod.F90
```

- ajout du répertoire `.ssm.d`, qui pourrait être utilisé pour l'installation
  de MIDAS avec la commande `make package`

- ajout du fichier `config.in`, qui permet de créer un fichier nommé
  `midas-config`, listant le compilateur, la version de MIDAS, les
  librairies utilisées, etc.

- le fichier `tools/splitobs/splitobs.c` a été modifié par Jean-Philippe,
  pour montrer l'utilisation de la librairie App (installée avec librmn)
  pour gérer les messages d'erreur

À noter que l'option de compilation `-warn errors` est enlevée, car elle fait
échouer la compilation à cause de l'option `-stand f08`, qui est maintenant
ajoutée à nos règles.

Par ailleurs, nous pourrions rediscuter du nom des programmes installés, qui
inclut la plateforme et la version de midas: ce ne serait sans doute pas
nécessaire.

Pour l'instant, les binaires sont installés directement dans le répertoire
d'installation, je propose des les installer dans un sous-répertoire `bin`,
comme dans l'installation `ssm`.

Vos commentaires sont bienvenus!
