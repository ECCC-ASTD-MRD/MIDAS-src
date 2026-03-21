midas_sources = run_command('sh', '-c', 'ls *.f90 *.F90', check: true).stdout().split()

midas_deps = [
  rmn_dep,
  rpncomm_dep,
  vgrid_dep,
  hpcoperf_dep,
  hdf5_dep,
  netcdf_dep,
  zlib_dep,
  curl_dep
]

if not mkl_support
  midas_deps += [blas_dep, lapack_dep]
endif

midas_lib = library('midas',
  midas_sources,
  dependencies: midas_deps,
  include_directories: [include_directories('../../')] + rttov_inc, # Accessing build root and RTTOV
  install: true
)

midas_dep = declare_dependency(
  link_with: midas_lib,
  include_directories: [include_directories('.'), include_directories('../../')] + rttov_inc,
  dependencies: midas_deps
)
