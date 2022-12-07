MPI in MIDAS
============

MIDAS uses Message Passing Interface (MPI) to allow programs to be parallelized
over a potentially large number of processors on the supercomputers. When
submitted for execution, a MIDAS program will run many times independently in
parallel on separate groups of processors consisting of one or more threads
each. For example, if the `var` program is submitted in a maestro task using the
following topology::

  cpu="10x54x4"

then the program will run 540 (i.e. 10x54) times in parallel, each using 4
threads. These 540 executions are each referred to as an "MPI task". We use the
``rpn_comm`` library (part of ``rmnlib``) as a wrapper around the standard MPI
library routines. The ``rpn_comm`` library allows the decomposition of the MPI
tasks in two dimensions. In the above example, the 540 MPI tasks are decomposed
into 10 MPI tasks in the "x" direction and 54 MPI tasks in the "y"
direction. This is used in various contexts in MIDAS, including to distribute
the data of the data object associated with the ``gridStateVector_mod`` module
in the two horizontal dimensions. This is known as the "Tiles" MPI distribution.

.. image:: ../../docs/images/mpi_2D_decomposition.png
  :width: 700

**Figure 1:** Example of 2D domain decomposition of a global data object
associated with the ``gridStateVector_mod`` module with an MPI topology of
``8x4x1``. MPI task (7,1) is shown with its local grid values that cover a
sub-region over the South Pacific.

However, when performing various horizontal transformations on the data
(e.g. horizontal interpolation from one grid to another) it is more practical to
have access to all of the data for a particular vertical level and variable on a
single MPI task. In this case the same data can instead be distributed over the
MPI tasks using the so-called "VarsLevs" distribution. Instead of a subset of
horizontal grid points on each MPI task, the contatenation of the vertical
levels with the different variables is decomposed over all MPI tasks
(i.e. combining the tasks in the "x" and "y" directions).

.. image:: ../../docs/images/mpi_varslevs_decomposition.png
  :width: 700

**Figure 2:** Example of the "VarsLevs" MPI decomposition that distributes
complete 2D horizontal fields over the MPI tasks. In this case there are a total
of 6 MPI tasks, such as from an MPI topology of ``3x2x1``.

In addition to the "Tiles" and "VarsLevs" MPI distributions, MIDAS can also
store the entire gridded data on a single MPI task. This global distribution is
used when writing and reading to/from files on disk. In some situations,
individual time steps can be distributed over multiple MPI tasks to allow each
of the global 3D data to be written or read to/from the disk in parallel.

It is often necessary to be able to transform data that exists in one type of
MPI distribution into another MPI distribution. For example, to perform
horizontal and vertical interpolation with data that is originally in the
"Tiles" MPI distribution is performed by:

1. Transpose the original data from "Tiles" to "VarsLevs" MPI distribution.
2. Perform horizontal interpolation on all MPI tasks in parallel for a subset of
   vertical levels/variables on each.
3. Transpose the result from "VarsLevs" back to "Tiles" MPI distribution.
4. Perform vertical interpolation on all MPI tasks in parallel for a subset of
   horizontal grid points on each.

Therefore subroutines are included in MIDAS that perform various transpositions
between different pairs of MPI distributions.
