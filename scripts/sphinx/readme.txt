
Obtain the original sphinx-fortran extension and put in the directory sphinx-fortran:

From the centreal github repository:
git clone https://github.com/VACUMM/sphinx-fortran.git ${HOME}/sphinx-fortran

Or from our local (modified) copy:
git clone git@gitlab.science.gc.ca:mab001/sphinx-fortran.git ${HOME}/sphinx-fortran

Install sphinx-fortran extension in the midas directory (assuming midas depot is ${HOME}/midas/:

cd ${HOME}/sphinx-fortran
rm ${HOME}/midas/scripts/sphinx/lib/python2.7/*
python setup.py install --install-lib=${HOME}/midas/scripts/sphinx/lib/python2.7

Run the script to build the documentation:

cd ${HOME}/midas/scripts/sphinx/
./build_html.sh
