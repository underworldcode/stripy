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
  input : ['tripack.pyf', 'tripack.f90'],  # .f so no F90 wrappers
  output : ['_tripackmodule.c', '_tripack-f2pywrappers.f', '_tripack-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_tripack', '--lower']
)

py.extension_module('_tripack',
  ['tripack.f90', tripack_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
)

# py.extension_module('_tripack',
#   ['src/tripack.pyf',	'src/tripack.f90'],
#   include_directories: inc_np,
#   dependencies: py_dep,
#   install: true,
# )

# py.extension_module(
# 	'_srfpack',
#   'src/srfpack.pyf',
# 	'src/srfpack.f',
#   include_directories: inc_np,
#   dependencies: py_dep,
#   install: true,
# )

# py.extension_module(
# 	'_ssrfpack',
#   'src/ssrfpack.pyf',
# 	'src/ssrfpack.f',
#   include_directories: inc_np,
#   dependencies: py_dep,
#   install: true,
# )

# py.extension_module(
# 	'_fortran',
#   'src/stripyf.pyf',
# 	'src/stripyf.f90',
#   include_directories: inc_np,
#   dependencies: py_dep,
#   install: true,
# )
