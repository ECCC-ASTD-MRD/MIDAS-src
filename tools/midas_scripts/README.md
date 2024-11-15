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

 * `MIDAS_CONCURRENT_SPLITOBS`: Cette variable controle le nombre de
   processus de splitting lancés en parallèle dans le script
   `midas.mpi`.  Par défaut, elle n'est pas définie et ainsi aucune
   limite n'est imposée.

 * `MIDAS_OBS_MPI_ORDERING`: Si cette variable est égale à `inverse`,
   alors on utilisera la répartition MPI des observations telle que
   cela était utilisée dans les versions `v_3.3.*` et antérieures de
   MIDAS.  Si elle est égale à `regular`, alors on utilisera la
   répartition définie dans les versions de MIDAS postérieures à
   `v_3.3`.

 * `MIDAS_MPI_BARRIER_VERBOSE`: Si cette variable est égale à `yes`,
   alors un `set -x` sera effectué dans le script `midas.mpi_barrier`.

## `midas.launch`

Le script `midas.launch` est le script principal qui lance les
programmes MIDAS qui servent dans l'assimilation et le contrôle de
qualité.  Ce script appelle les scripts suivants:
 * `midas.check_ensemble`
 * `midas.tripotenml`
 * `midas.mpirun`
 * `midas.mpi`
   * `midas.splitobs`
   * `midas.mpi_barrier`
 * `midas.reunir_obs`
   * `midas.reunir_obs_mpi`

### `midas.check_ensemble`

Ce script verifie si tous les trials d'ensemble sont présents.  Dans
le cas où on appelle le script avec `-fallback_mode` égal à `continue`
et qu'il manque des membres d'ensemble, on applique l'[algorithme de
contigence](https://wiki.cmc.ec.gc.ca/wiki/RPN-AD/Ensemble_contingency/FullDescription).

### `midas.tripotenml`

Ce script sert à modifier des entrées dans un namelist fortran.  On
l'utilise pour changer la valeur de l'étiquette.

#### `midas.mpirun`

Ce script est un wrapper autour de `r.run_in_parallel` pour rassembler
le code qui doit être exécuté avant de lancer le MPI.

### `midas.mpi`

Ce script est lancé en MPI et est celui qui appelle vraiment les
programmes de MIDAS.  C'est ce script qui appelle `midas.splitobs` qui
fait le splitting en parallèle des observations avant d'appeler les
véritables programmes MPI.  On peut limiter le nombre de processus en
parallèle envoyés avec la variable `MIDAS_CONCURRENT_SPLITOBS`.

#### `midas.splitobs`

Ce script fait le splitting des observations avant de lancer le
programme MIDAS qui a besoin des observations déjà splittées.  On
supporte les fichiers BURP et SQLite.

#### `midas.mpi_barrier`

Ce script permet de resynchronier toutes les tuiles MPI après le
splitting des observations et avant de lancer le programme MPI.
Sinon, on obtient des timeout avec des erreurs `Alarm call` sur les
PPP.

### `midas.reunir_obs`

Ce script est le driver principal pour rassembler les observations qui
ont ete splittées par `midas.splitobs` et modifiées par le programme
MIDAS.  On appelle `midas.reunir_obs_mpi` pour chaque famille.

#### `midas.reunir_obs_mpi`

Ce script est lancé par `midas.reunir_obs` et fait le rassemblage des
observations pour une seule famille.
On peut l'appeler interactivement avec:
```bash
fam=al
rm -f obsfiles_${fam}.updated/{obs${fam},TABLES_REUNIR} obs${fam}
midas.reunir_obs_mpi -obsin $PWD -obsout $PWD -families2process ${fam}
```

## Running MIDAS as a stand alone program

The scripts `midas.prepare_workdir` and `midas.launch_program` can be
used to prepare a working directory and then run a MIDAS program as a
standalone execution.
  * `midas.prepare_workdir`: This script prepares a working directory
     using inputs from `cmcarc` archives or a directory for a fixed MPI
     topology.  The last point is essential because the options will be
     prepared for that MPI topology.
  * `midas.launch_program`: This script launches a job on the
     supercomputer using the working directory prepared by
     `midas.prepare_workdir`.  The job is using the MPI topology that
     has been used to prepare the observations.

There is another script, `midas.unsplitobs`, that can be used to
recombine the observations to be in their original order before they
were split with `midas.splitobs.Abs`.

All three scripts supports the option `-h` to have a description of
the available options.

### Example

Here is an example of utilization of those three scripts where we use
the test `/Tests/letkf/glb_xc40` as a reference:

```bash
cd $(git rev-parse --show-toplevel)/tools/midas_scripts

program_directory=/fs/ssm/eccc/mrd/rpn/anl/midas/4.0.2/rhel-8-icelake-64/bin
inputs=/home/sanl000/data_maestro/ppp5/UnitTests/midas/letkf/glb_xc40/0013
namelist=$(git rev-parse --show-toplevel)/maestro/suites/midas_system_tests/config/Tests/letkf/glb_xc40/nml

workdir=${HOME}/data_maestro/ppp6/tmp/midas_letkf_workdir_3
observations=${HOME}/data_maestro/ppp6/tmp/midas_letkf_observations

./midas.unsplitobs -input ${inputs}/inputs_obsfiles.ca -output ${observations}             \
                   -prefix_dir burpfiles_ -prefix_obs brp                                  \
                   -workdir ${HOME}/data_maestro/ppp6/tmp/midas_letkf_unsplit_observations \
                   -splitobs ${program_directory}/midas.splitobs.Abs

./midas.prepare_workdir -workdir ${workdir}                    \
                        -nml ${namelist}                       \
                        -ensemble ${inputs}/inputs_ensemble.ca \
                        -observations ${observations}          \
                        -constants ${inputs}/inputs.ca         \
                        -npex 36 -npey 18                      \
                        -splitobs ${program_directory}/midas.splitobs.Abs

./midas.launch_program -workdir ${workdir} -pgm ${program_directory}/midas-letkf.Abs
```

The results have been verified with the reference for the test
`/Tests/letkf/glb_xc40`.

If we would have the original unsplit files, running
`midas.unsplitobs` can be avoided.

### `midas.energyNorm`

The program `midas-energyNorm` is useful to compare two atmospheric
states.

The script `midas.energyNorm` is the tool to prepare the working
directory for that program:
```bash
midas.energyNorm [positionnels] This script launches the program 'midas-energyNorm.Abs'
 IN       -pgm [:] path to MIDAS program
 IN       -nml [:] path to namelist file
 IN       -date [:] date of the input files
 IN       -reference [:] path to reference file
 IN       -states [:] path to several files on which the energy norm will be computed
 IN       -workdir [:] path to working directory prepared by 'midas.prepare_workdir'
 IN       -supercomputer [ppp6:ppp6] supercomputer to launch the job (default: 'ppp6')
 IN       -jobname [energyNorm:energyNorm] name of the job (default: 'letkf')
 IN       -memory [:] memory request per MPI rank (default: 2.25G*${omp_num_threads})
 IN       -npex [1:1] Number of mpi components in the 'x' direction (default: 1)
 IN       -npey [1:1] Number of mpi components in the 'y' direction (default: 1)
 IN       -omp_num_threads [1:1] Number of Open MP threads per MPI rank (default: 1)
 IN       -wallclocktime [10:10] wall clock time for job in minutes (default: 10
```

We suggest to use the namelist from the test
[`/Tests/energyNorm/analmean`](maestro/suites/midas_system_tests/config/Tests/energyNorm/analmean/nml).
