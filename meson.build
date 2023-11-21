project('stripy', 'c')

name = 'stripy'

add_languages('fortran')

platform = host_machine.system()
if platform == 'windows'
  add_project_link_arguments('-static', language: ['fortran', 'c'])
elif platform == 'darwin'
  add_project_link_arguments('-Wl,-rpath, "@loader_path"', language: ['fortran', 'c'])
else
  add_project_link_arguments('-Wl,-rpath,"$ORIGIN"', language: ['fortran', 'c'])
endif

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

install_subdir(name, install_dir: py.get_install_dir() / name, strip_directory: true)

tripack_source = custom_target('_tripackmodule.c',
  input : ['src/tripack.pyf'],  # .f so no F90 wrappers
  output : ['_tripackmodule.c', '_tripack-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@']
)

stripack_source = custom_target('_stripackmodule.c',
  input : ['src/stripack.pyf'],  # .f so no F90 wrappers
  output : ['_stripackmodule.c', '_stripack-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@']
)

srfpack_source = custom_target('_srfpackmodule.c',
  input : ['src/srfpack.pyf'],  # .f so no F90 wrappers
  output : ['_srfpackmodule.c'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@']
)

ssrfpack_source = custom_target('_ssrfpackmodule.c',
  input : ['src/ssrfpack.pyf'],  # .f so no F90 wrappers
  output : ['_ssrfpackmodule.c', '_ssrfpack-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@']
)

fortran_source = custom_target('_fortranmodule.c',
  input : ['src/stripyf.pyf'],  # .f so no F90 wrappers
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
