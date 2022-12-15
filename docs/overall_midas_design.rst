Overall MIDAS design
====================

The MIDAS code is contained within a set of main programs and FORTRAN
modules. Each program "uses" a subset of the FORTRAN modules. These modules also
"use" other modules. In addition to the MIDAS programs and modules, several
external libraries are also employed to provide specific functionality
(e.g. MPI, RTTOV, SQLITE, RMNLIB). The programs and modules are located in the
subdirectories::

  /src/programs
  /src/modules

Because modules can "use" other modules, it is important to be aware of the
hierarchical relationship between the modules to avoid any circular
dependencies (i.e. ``moduleA`` "uses" ``moduleB`` and ``moduleB`` "uses"
``moduleA``). In general, this can be avoided by limiting the type of
responsibilities each module is assigned. For example, a high-level module
responsible for implementing a high-level algorithm that manipulates high-level
data objects should not also contain low-level routines that could be useful to
many other modules, since this is likely to create a circular dependency.

The MIDAS modules can be divided into several general categories, listed in the
order that they exist in the dependency hierarchy (i.e. higher-level modules can
"use" lower-level modules, but not vice versa):

* High-level functionality (e.g. ``innovation_mod``, ``minimization_mod``)
* Transformation of data objects (e.g. ``stateToColumn_mod``, ``gridVariableTransforms_mod``)
* High-level data objects (e.g. ``gridStateVector_mod``, ``columnData_mod``, ``obsSpaceData_mod``)
* Low-level data objects (e.g. ``horizontalCoord_mod``, ``verticalCoord_mod``)
* Low-level utilities (e.g. ``physicsFunctions_mod``, ``utilities_mod``)
