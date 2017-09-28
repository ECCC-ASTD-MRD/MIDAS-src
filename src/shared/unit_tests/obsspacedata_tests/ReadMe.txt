How to run unit tests on the ObsSpaceData code:
===============================================
1.  The tests should run on either Linux or the supercomputers.
2.  cd <the directory where you have checked out obsspacedata_mod.ftn90>
3.  cd unit_tests/obsspacedata_tests
4.  On supercomputer, enter:  gmake clean; gmake unit_tests
5.  On Linux,         enter:   make clean;  make unit_tests
6.  On the screen, you will be told where to find the listing.
7.  In the listing, search for "T E S T S   R U N" and it will tell you pass or fail.

PART 2
======
8.  Open testobsdriver.ftn90 in an editor.
9.  Reverse the commenting of lines 30, 31. (They start with "call add(suite_master".)
10. Repeat steps 4-7 above.
11. Undo your edit in testobsdriver.ftn90, so as to avoid unnecessary archived edits.

Any difficulties may be addressed to Jeff Blezius x4728.
