project('midas', 'c', 'fortran',
  version: run_command('sh', '-c', './midas.version --abbrev=0 HEAD', check: true).stdout().strip(),
  meson_version: '>= 0.54.0',
  default_options: ['buildtype=release', 'fortran_std=f2018']
)

fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')

# Project-wide options
sqlite_fortran_support = get_option('sqlite_fortran_support')
mkl_support = get_option('mkl_support')

if fc.get_id() == 'intel' or fc.get_id() == 'intel-cl'
  mkl_support = true
endif

# Dependencies
mpi_dep = dependency('mpi', language: 'fortran')
openmp_dep = dependency('openmp', language: 'fortran')

# For the EC libraries, we'll try to find them using dependency() first, 
# then fallback to find_library if needed.
# Since I don't know if they have pkg-config files, I'll start with find_library
# for those that look like internal/EC-specific.

rmn_dep = dependency('rmn', required: false)
if not rmn_dep.found()
  rmn_dep = fc.find_library('rmn', required: true)
endif

vgrid_dep = dependency('vgrid', required: false)
if not vgrid_dep.found()
  vgrid_dep = fc.find_library('vgrid', required: true)
endif

rpncomm_dep = dependency('rpncomm', required: false)
if not rpncomm_dep.found()
  rpncomm_dep = fc.find_library('rpncomm', required: true)
endif

random_dep = dependency('random', required: false)
if not random_dep.found()
  random_dep = fc.find_library('random', required: true)
endif

sqlite3_dep = dependency('sqlite3', required: true)
burp_dep = dependency('burp', required: false)
if not burp_dep.found()
  burp_dep = fc.find_library('burp', required: true)
endif

hdf5_dep = dependency('hdf5', language: 'fortran', required: true)
netcdf_dep = dependency('netcdf', language: 'fortran', required: true)
curl_dep = dependency('libcurl', required: true)
zlib_dep = dependency('zlib', required: true)

if not mkl_support
  blas_dep = dependency('blas', required: true)
  lapack_dep = dependency('lapack', required: true)
endif

rttov_inc = []
rttov_libs = []
rttov_installdir = run_command('sh', '-c', 'echo $rttov_INSTALLDIR', check: true).stdout().strip()
if rttov_installdir != ''
  rttov_inc += include_directories(rttov_installdir + '/include')
  rttov_inc += include_directories(rttov_installdir + '/mod')
endif

rttov_components = ['rttov_coef_io', 'rttov_hdf', 'rttov_parallel', 'rttov_emis_atlas', 'rttov_other', 'rttov_mw_scatt', 'rttov_main']
foreach comp : rttov_components
  rttov_libs += fc.find_library(comp, required: true)
endforeach

hpcoperf_dep = fc.find_library('hpcoperf', required: true)

# Compiler flags
fflags = []
if fc.get_id() == 'intel' or fc.get_id() == 'intel-cl'
  if get_option('buildtype') == 'release'
    fflags += ['-O3', '-fast-transcendentals', '-no-prec-div', '-fpic', '-ip', '-no-prec-sqrt']
  elif get_option('buildtype') == 'debug'
    fflags += ['-check', 'all', '-fp-speculation=safe', '-init=snan,arrays']
  endif
  fflags += ['-qmkl', '-check', 'noarg_temp_created', '-no-wrap-margin', '-warn', 'all', '-warn', 'noexternal', '-static-intel']
elif fc.get_id() == 'intel-llvm'
  if get_option('buildtype') == 'release'
    # Default O3 is usually fine for icx/ifx
  elif get_option('buildtype') == 'debug'
    fflags += ['-check', 'all', '-fp-speculation=safe', '-init=snan,arrays']
  endif
  fflags += ['-qmkl', '-check', 'noarg_temp_created', '-no-wrap-margin', '-warn', 'all', '-warn', 'noexternal']
endif

add_project_arguments(fflags, language: 'fortran')
add_project_link_arguments('-Wl,--disable-new-dtags', language: ['c', 'fortran'])

# Configuration and Build Info
conf_data = configuration_data()
conf_data.set_quoted('GIT_VERSION', meson.project_version())
configure_file(output : 'midas_build_info.h', configuration : conf_data)

# midas-config generation
config_in = configuration_data()
config_in.set('NAME', meson.project_name())
config_in.set('MIDAS_BASE_VERSION', meson.project_version())
config_in.set('CMAKE_C_COMPILER_ID', cc.get_id())
config_in.set('CMAKE_C_COMPILER_VERSION', cc.version())
config_in.set('CMAKE_Fortran_COMPILER_ID', fc.get_id())
config_in.set('CMAKE_Fortran_COMPILER_VERSION', fc.version())

# For flags, we can use the ones we've collected
config_in.set('EC_C_FLAGS', ' '.join(get_option('c_args')))
config_in.set('EC_Fortran_FLAGS', ' '.join(fflags))
config_in.set('EC_CMAKE_DEFINITIONS', '') # Simplified

# Dependency versions
# For dependencies found via find_library, .version() might not work.
# We'll try to get it, or use 'unknown'.
rmn_version = rmn_dep.found() ? rmn_dep.version() : 'unknown'
rpncomm_version = rpncomm_dep.found() ? rpncomm_dep.version() : 'unknown'
vgrid_version = vgrid_dep.found() ? vgrid_dep.version() : 'unknown'
hdf5_version = hdf5_dep.found() ? hdf5_dep.version() : 'unknown'
sqlite3_version = sqlite3_dep.found() ? sqlite3_dep.version() : 'unknown'

config_in.set('rmn_VERSION', rmn_version)
config_in.set('rpncomm_VERSION', rpncomm_version)
config_in.set('vgrid_VERSION', vgrid_version)
config_in.set('HDF5_VERSION', hdf5_version)
config_in.set('SQLite3_VERSION', sqlite3_version)
config_in.set('rttov_VERSION', 'unknown') # Simplified

config_in.set('EC_CI_ENV', 'unknown')
config_in.set('EC_ARCH', 'unknown')

configure_file(
  input : 'config.in',
  output : 'midas-config',
  configuration : config_in,
  install_dir : get_option('bindir')
)

# Subdirectories
subdir('src/modules')
subdir('src/programs')
subdir('tools')
