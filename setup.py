## To install locally: python setup.py build && python setup.py install
## (If there are problems with installation of the documentation, it may be that
##  the egg file is out of sync and will need to be manually deleted - see error message
##  for details of the corrupted zip file. )
##
## To push a version through to pip.
##  - Make sure it installs correctly locally as above
##  - Update the version information in this file
##  - python setup.py sdist upload -r pypitest  # for the test version
##  - python setup.py sdist upload -r pypi      # for the real version
##
## (see http://peterdowns.com/posts/first-time-with-pypi.html)


from setuptools import setup, find_packages
from numpy.distutils.core import setup, Extension
from os import path
import io

## in development set version to none and ...
PYPI_VERSION = "1.1"

# Return the git revision as a string (from numpy)
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', '--short', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


if PYPI_VERSION is None:
    PYPI_VERSION = git_version()


this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# interface for Renka's algorithm 772 fortran code
ext1 = Extension(name    = 'stripy._stripack',
                 sources = ['src/stripack.pyf','src/stripack.f90'])
ext2 = Extension(name    = 'stripy._tripack',
                 sources = ['src/tripack.pyf', 'src/tripack.f90'])
ext3 = Extension(name    = 'stripy._srfpack',
                 sources = ['src/srfpack.pyf', 'src/srfpack.f'])
ext4 = Extension(name    = 'stripy._ssrfpack',
                 sources = ['src/ssrfpack.pyf', 'src/ssrfpack.f'])

if __name__ == "__main__":
    setup(name = 'stripy',
          author            = "Louis Moresi",
          author_email      = "louis.moresi@unimelb.edu.au",
          url               = "https://github.com/underworldcode/stripy",
          version           = PYPI_VERSION,
          description       = "Python interface to TRIPACK and STRIPACK fortran code for triangulation/interpolation in Cartesian coordinates and on a sphere",
          long_description  = long_description,
          long_description_content_type='text/markdown',
          ext_modules       = [ext1, ext2, ext3, ext4],
          install_requires  = ['numpy', 'scipy>=0.15.0'],
          python_requires   = '>=2.7, >=3.5',
          setup_requires    = ["pytest-runner", 'webdav'],
          tests_require     = ["pytest", 'webdav'],
          packages          = ['stripy'],
          package_data      = {'stripy': ['Notebooks/*ipynb', # Worked Examples is not currently used
                                          'Notebooks/CartesianTriangulations/*ipynb',
                                          'Notebooks/SphericalTriangulations/*ipynb',
                                          'Notebooks/Data/*'] },
          include_package_data = True,
          classifiers       = ['Programming Language :: Python :: 2',
                               'Programming Language :: Python :: 2.6',
                               'Programming Language :: Python :: 2.7',
                               'Programming Language :: Python :: 3',
                               'Programming Language :: Python :: 3.3',
                               'Programming Language :: Python :: 3.4',
                               'Programming Language :: Python :: 3.5',
                               'Programming Language :: Python :: 3.6',
                               'Programming Language :: Python :: 3.7'
                               ]
          )
