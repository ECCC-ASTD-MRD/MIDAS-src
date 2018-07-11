
Obtain the original sphinx-fortran extension and put in the directory sphinx-fortran:

git clone https://github.com/VACUMM/sphinx-fortran.git ${HOME}/sphinx-fortran


Install sphinx-fortran extension in the midas directory (assuming midas depot is ${HOME}/midas/:

cd ${HOME}/sphinx-fortran
python setup.py install --install-lib=${HOME}/midas/scripts/sphinx/lib/python2.7

