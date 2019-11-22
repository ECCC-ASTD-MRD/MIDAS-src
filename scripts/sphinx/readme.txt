
Obtain the original sphinx-fortran extension and put in the directory sphinx-fortran:

From the centreal github repository:
git clone https://github.com/VACUMM/sphinx-fortran.git ${HOME}/sphinx-fortran

Or from our local (modified) copy:
git clone git@gitlab.science.gc.ca:mab001/sphinx-fortran.git ${HOME}/sphinx-fortran

There are two ways to connect to sphinx-fortran

1) INSTALL THE EXTENSION AS A *.egg FILE
========================================
Install sphinx-fortran extension in the midas directory (assuming midas depot is ${HOME}/midas/:

cd ${HOME}/sphinx-fortran
rm ${HOME}/midas/scripts/sphinx/lib/python2.7/*
python setup.py install --install-lib=${HOME}/midas/scripts/sphinx/lib/python2.7

2) CREATE A SOFT LINK TO THE EXTENSION
======================================
N.B.:  This is for development only because, with this method, sphinx-fortran is not
       git-archived with the midas code
cd …/sphinx/lib/python2.7
rm sphinx_fortran-1.0.1-py2.7.egg easy-install.pth
cd …/sphinx-fortran
python setup.py develop --install-dir=/home/jbl001/GIT-depots/envar/scripts/sphinx/lib/python2.7

TO BUILD THE DOCUMENTATION
==========================
Run the script to build the documentation:

cd ${HOME}/midas/scripts/sphinx/
./build_html.sh
