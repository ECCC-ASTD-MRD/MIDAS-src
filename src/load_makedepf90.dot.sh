#! /bin/sh
#CONDA_VENV=/fs/homeu1/eccc/aq/ords/mad001/conda/envs/midas
CONDA_VENV=/home/mad001/.conda/envs/midas-dev


##  loading conda
. ssmuse-sh -x hpco/exp/mib002/anaconda2/anaconda2-5.0.1-hpcobeta2

source activate ${CONDA_VENV}

##  disable conda prompt
##  (otherwise very very long.  This is all in the mean time until 
##  `makedepf90` is globally installed)
##  a cleaner hack would be to change the prompt...
conda config --set changeps1 False
