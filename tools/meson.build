splitobs_sources = ['splitobs/splitobs.c', 'splitobs/fstdlib.c']

executable('midas.splitobs.Abs',
  splitobs_sources,
  dependencies: [sqlite3_dep, burp_dep, rmn_dep],
  include_directories: include_directories('..'), # for the root build dir
  install: true
)
