project('stripy', 'c')

name = 'stripy'

add_languages('fortran')

platform = host_machine.system()

# These link arguments make each extension module self-contained, but they are
# GNU/gfortran conventions. Other Fortran toolchains (LLVM flang, Intel, the
# MSVC-based linkers) warn about, reject, or silently ignore them and provide
# their runtime libraries differently, so only emit them for a GNU (gfortran/GCC)
# compiler (see gh-116). Non-GNU builds rely on the packaging step
# (delocate/auditwheel, or conda compiler activation) to resolve the runtime.
#
# This stays correct under cross-compilation: both host_machine.system() and
# meson.get_compiler('fortran') refer to the *host* (target) machine, so a cross
# build with a gfortran cross-compiler (e.g. conda-forge linux-aarch64 or
# osx-arm64) still selects the right per-target flags.
if meson.get_compiler('fortran').get_id() == 'gcc'
  if platform == 'windows'
    # Statically link libgfortran/libquadmath (MinGW) so the .pyd is standalone.
    add_project_link_arguments('-static', language: ['fortran', 'c'])
  elif platform == 'darwin'
    add_project_link_arguments('-Wl,-rpath,@loader_path', language: ['fortran', 'c'])
  else
    add_project_link_arguments('-Wl,-rpath,$ORIGIN', language: ['fortran', 'c'])
  endif
endif

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
  ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
  check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

# Install the pure-Python package, including the shipped example Notebooks that
# stripy.documentation relies on. The excludes drop common editor/interpreter
# cruft; note meson matches exact relative paths (no globbing), so these only
# catch top-level occurrences. Clean release wheels are built by CI from a fresh
# checkout, so this is just a best-effort safety net for local (dirty) installs.
install_subdir(name,
  install_dir: py.get_install_dir() / name,
  strip_directory: true,
  exclude_files: ['.DS_Store'],
  exclude_directories: ['__pycache__', '.ipynb_checkpoints'],
)

# f2py turns each .pyf interface into a C extension source (<name>module.c) and,
# when the interface requires it, a Fortran wrapper (<name>-f2pywrappers.f). The
# `output` lists below must match exactly what the pinned f2py (numpy>=2.0)
# emits; if a future f2py changes its output set, update these lists to match or
# the build will fail with a missing/unexpected output error.
tripack_source = custom_target('_tripackmodule.c',
  input : ['src/tripack.pyf'],
  output : ['_tripackmodule.c', '_tripack-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@']
)

stripack_source = custom_target('_stripackmodule.c',
  input : ['src/stripack.pyf'],
  output : ['_stripackmodule.c', '_stripack-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@']
)

srfpack_source = custom_target('_srfpackmodule.c',
  input : ['src/srfpack.pyf'],
  output : ['_srfpackmodule.c'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@']
)

ssrfpack_source = custom_target('_ssrfpackmodule.c',
  input : ['src/ssrfpack.pyf'],
  output : ['_ssrfpackmodule.c', '_ssrfpack-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@']
)

fortran_source = custom_target('_fortranmodule.c',
  input : ['src/stripyf.pyf'],
  output : ['_fortranmodule.c', '_fortran-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@']
)


py.extension_module('_tripack',
  ['src/tripack.f90', tripack_source],
  incdir_f2py / 'fortranobject.c',
  subdir: 'stripy',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
)

py.extension_module('_stripack',
  ['src/stripack.f90', stripack_source],
  incdir_f2py / 'fortranobject.c',
  subdir: 'stripy',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
)

py.extension_module('_srfpack',
  ['src/srfpack.f', srfpack_source],
  incdir_f2py / 'fortranobject.c',
  subdir: 'stripy',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
)

py.extension_module('_ssrfpack',
  ['src/ssrfpack.f', ssrfpack_source],
  incdir_f2py / 'fortranobject.c',
  subdir: 'stripy',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
)

py.extension_module('_fortran',
  ['src/stripyf.f90', fortran_source],
  incdir_f2py / 'fortranobject.c',
  subdir: 'stripy',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
)
