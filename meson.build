project('stripy', 'c')

add_languages('fortran')

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

tripack_source = custom_target('_tripackmodule.c',
  input : ['src/tripack.pyf', 'src/tripack.f90'],  # .f so no F90 wrappers
  output : ['_tripackmodule.c', '_tripack-f2pywrappers.f', '_tripack-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_tripack', '--lower']
)

stripack_source = custom_target('_stripackmodule.c',
  input : ['src/stripack.pyf', 'src/stripack.f90'],  # .f so no F90 wrappers
  output : ['_stripackmodule.c', '_stripack-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_stripack', '--lower']
)

srfpack_source = custom_target('_srfpackmodule.c',
  input : ['src/srfpack.pyf', 'src/srfpack.f'],  # .f so no F90 wrappers
  output : ['_srfpackmodule.c', '_srfpack-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_srfpack', '--lower']
)

ssrfpack_source = custom_target('_ssrfpackmodule.c',
  input : ['src/ssrfpack.pyf', 'src/ssrfpack.f'],  # .f so no F90 wrappers
  output : ['_ssrfpackmodule.c', '_ssrfpack-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_ssrfpack', '--lower']
)

fortran_source = custom_target('_fortranmodule.c',
  input : ['src/stripyf.pyf', 'src/stripyf.f90'],  # .f so no F90 wrappers
  output : ['_fortranmodule.c', '_fortran-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_fortran', '--lower']
)


py.extension_module('_tripack',
  ['src/tripack.f90', tripack_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
)

py.extension_module('_stripack',
  ['src/stripack.f90', stripack_source],
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
)

py.extension_module('_srfpack',
  ['src/srfpack.f', srfpack_source],
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
)

py.extension_module('_ssrfpack',
  ['src/ssrfpack.f', ssrfpack_source],
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
)

py.extension_module('_fortran',
  ['src/stripyf.f90', fortran_source],
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
)