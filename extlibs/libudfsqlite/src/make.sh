#!/bin/bash


d.compile -conly -o libudfsqlite-shared.so                 -src udfsqlite.c -librmn rmn       -libsys z   -debug
d.compile -conly -o libudfsqlite.a -defines =-DSQLITE_CORE -src udfsqlite.c -librmn rmn       -libsys z m -debug

