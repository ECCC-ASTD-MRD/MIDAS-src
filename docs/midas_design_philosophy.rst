MIDAS design philosophy
=======================

* MIDAS uses the **FORTRAN** programming language exclusively.

.. |br| raw:: html

* We aim to apply **an incremental approach** to improving the design and
  programming style of the MIDAS code. In particular, it is required that any
  new code and significantly modified existing code is consistent with the
  desired MIDAS code design and programming style. Otherwise, making arbitrary
  changes to code design and style should be avoided, except for specific gitlab
  issues focused soley on such improvements

.. |br| raw:: html

* **Object-oriented design:** To the extent possible, data (i.e. derived type
  definition and module variables) and related code (subroutines and functions)
  should be grouped together within a FORTRAN module. Where appropriate, the
  composition type of inheritance (i.e. "has a" inheritance) is used to
  construct more complex object by combining several simpler objects. For
  example, the ``gridStateVector`` object is composed of the objects defined in
  the ``horizontalCoord_mod``, ``verticalCoord_mod`` and ``timeCoord_mod``
  modules.

.. |br| raw:: html

* **Modularity and encapsulation:** In FORTRAN this means avoiding the use of
  ``public`` variables, except for very limited cases. In addition, subroutines
  and functions should be ``private`` whenever appropriate. All subroutines,
  functions and variables accessible by code outside of the module must be
  explicitly defined as ``public`` at the beginning of the module, while the
  default is set to ``private`` for all variables, subroutines and functions.

.. |br| raw:: html

* **Clear relationship between MIDAS modules:** Each FORTRAN module in MIDAS has
  a short (preferrably 3 letters) prefix associated with it to clearly identify
  all public subroutines/functions/variables. In addition, public derived types
  defined in a module are usually named with the prefix ``struct_`` followed by
  the module prefix, for example, ``struct_hco`` for the module
  ``horizontalCoord_mod``. To make obvious the relationship between a FORTRAN
  module and the rest of the MIDAS code, all ``use`` statements appear at the
  beginning of the module declaration (and not within contained
  subroutines/functions). For a similar reason, all public entities are
  explicitly declared ``public`` near the beginning of the source file just
  after the ``use`` statements.
