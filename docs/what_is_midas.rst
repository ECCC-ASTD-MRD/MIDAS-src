What is MIDAS?
==============

The Modular and Integrated Data Assimilation System (MIDAS) aims to provide a
unified data assimilation framework for research and operational use at
Environment and Climate Change Canada with the objective of avoiding redundant
work within the organization and to facilitate the transition towards coupled
data assimilation of different components of the Earth system.

MIDAS evolved out of the operational FORTRAN software developed for variational
data assimilation (3D-Var and 4D-Var) for global and regional atmospheric
applications. It consists of a collection of programs and modules written
exclusively in FORTRAN for a growing number of applications related to data
assimilation, including the currently operational data assimilation, observation
quality control and observation bias correction for the atmospheric
deterministic and ensemble prediction systems (i.e. GDPS, GEPS, RDPS, HRDPS,
REPS). Other Earth system applications, including ocean and sea-ice, have also
been implemented in MIDAS. In addition, a MIDAS program for estimating the
impact of observations recently became operational at ECCC using the Hybrid
Forecast Sensitivity to Observation Impact approach

To facilitate the efficient development and continual improvement of MIDAS by
many developers in parallel, a :doc:`code design philosophy
<midas_design_philosophy>` was developed based on modularity, encapsulation and
some principles of object-oriented programming.
