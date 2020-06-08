Ces scripts servent à lancer les programmes de
[MIDAS](https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas)
et à manipuler les fichiers d'input ou bien d'output.

## Variables d'environnement

Les variables d'environnement suivantes sont possibles:
 * `MIDAS_RAMDISKPATH`: Cette variable définit le path vers le RAMDisk
   Par defaut, elle est égale à `${TMPDIR}`.
 * `MIDAS_SAVE_SPLITOBS`: Si cette variable est égale à `yes`, on
   copie les observations qui ont été splittés mais conservés en
   RAMDisk dans le répertoire de travail de `midas.Abs`.  Par défaut,
   cette variable est vide.  Cela est utile pour débogguer le
   splitting des observations et pour examiner ce que les programmes
   MIDAS lisent comme fichier d'observations.
 * `MIDAS_OBS_MPI_ORDERING`: Si cette variable est égale à `inverse`,
   alors on utilisera la répartition MPI des observations telle que
   cela était utilisée dans les versions `v_3.3.*` et antérieures de
   MIDAS.  Si elle est égale à `regular`, alors on utilisera la
   répartition définie dans les versions de MIDAS postérieures à
   `v_3.3`.
 * `MIDAS_MPI_BARRIER_VERBOSE`: Si cette variable est égale à `yes`,
   alors un `set -x` sera effectué dans le script `midas.mpi_barrier`.
 * `MIDAS_INTERPENSTRIALS_THREADS_PER_MEMBER' : Cette variable définit le
   nombre de threads utilisés pour traiter en parallèle chacun des
   trials d'ensemble dans le script `midas.interpEnsTrials.ksh`.  La
   parallélisation est faite sur les dates contenues dans chacun des
   trials d'ensemble.  Par défaut, cette variable est égale à 1.

## midas.launch

Le script `midas.launch` est le script principal qui lance les
binaires `midas-var`, `midas-oMinusf` qui servent dans l'assimilation
et le contrôle de qualité.  Ce script appelle les scripts suivants:
 * `midas.check_ensemble`
 * `midas.tripotenml`
 * `midas.mpirun`
 * `midas.var_mpi`
   * `midas.splitobs`
   * `midas.mpi_barrier`
 * `midas.reunir_obs`
   * `midas.reunir_obs_mpi`
 * `midas.interpEnsTrials_driver.ksh`
   * `midas.interpEnsTrials.ksh`
 * `midas.combineFstd`

### midas.check_ensemble

Ce script verifie si tous les trials d'ensemble sont présent avant de
demarrer une assimilation EnVar.  Dans le cas où on appelle le script
avec `-fallback_mode` égal à `continue` et qu'il manque des membres
d'ensemble, on applique l'[algorithme de
contigence](https://wiki.cmc.ec.gc.ca/wiki/RPN-AD/Ensemble_contingency/FullDescription).

### midas.tripotenml

Ce script sert à modifier des entrées dans un namelist fortran.  On
l'utilise pour changer la valeur de l'étiquette.

#### midas.mpirun

Ce script est un wrapper autour de `r.run_in_parallel` pour rassembler
le code qui doit être exécuté avant de lancer le MPI.

### midas.var_mpi

Ce script est lancé en MPI et est celui qui appelle vraiment les
programmes de MIDAS.  C'est ce script qui appelle `midas.splitobs` qui
fait le splitting des observations avant d'appeler les véritables
programmes MPI.

#### midas.splitobs

Ce script fait le splitting des observations avant de lancer le
programme MIDAS qui a besoin des observations déjà splittées.  On
supporte les fichiers BURP et SQLite.

#### midas.mpi_barrier

Ce script permet de resynchronier toutes les tuiles MPI après le
splitting des observations et avant de lancer le programme MPI.
Sinon, on obtient des timeout avec des erreurs `Alarm call` sur les
PPP.

### midas.reunir_obs

Ce script est le driver principal pour rassembler les observations qui
ont ete splittées par `midas.splitobs` et modifiées par le programme
MIDAS.  On appelle `midas.reunir_obs_mpi` pour chaque famille.

#### midas.reunir_obs_mpi

Ce script est lancé par `midas.reunir_obs` et fait le rassemblage des
observations pour une seule famille.
On peut l'appeler interactivement avec:
```bash
fam=al
rm -f obsfiles_${fam}.updated/{obs${fam},TABLES_REUNIR} obs${fam}
midas.reunir_obs_mpi -obsin $PWD -obsout $PWD -families2process ${fam}
```

## midas.interpEnsTrials_driver.ksh

Ce script gère l'interpolation des ensembles sur la grille horizontale
et verticale de l'analyse.  Ce script appelle
`midas.interpEnsTrials.ksh`.

### midas.interpEnsTrials.ksh

Ce script fait l'interpolation pour un membre d'ensemble en
particulier.  On y utilise amplement l'outil `pxs2pxt`.

## midas.combineFstd

Ce script rassemble plusieurs fichiers standard RPN en un seul.  Ce
script n'est pas utilisé mais on le laisse dans la librairie.
