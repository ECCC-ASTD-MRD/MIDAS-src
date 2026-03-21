program_list = ['energyNorm', 'letkf', 'ensPostProcess']

foreach prog : program_list
  source = prog + '.f90'
  deps = [mpi_dep, midas_dep]
  
  if prog == 'energyNorm'
    deps += [random_dep, hdf5_dep]
  elif prog == 'letkf'
    deps += [sqlite3_dep, burp_dep, random_dep, hdf5_dep] + rttov_libs
  elif prog == 'ensPostProcess'
    deps += [sqlite3_dep, random_dep, hdf5_dep] + rttov_libs
  endif
  
  executable('midas-' + prog + '.Abs', source,
    dependencies: deps,
    install: true
  )
endforeach
