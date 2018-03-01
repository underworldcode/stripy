from numpy.distutils.core import setup, Extension

# interface for Renka's algorithm 772 fortran code
ext1 = Extension(name    = 'stripy._stripack',
                 sources = ['src/stripack.pyf','src/stripack.f90'])
ext2 = Extension(name    = 'stripy._tripack',
                 sources = ['src/tripack.pyf', 'src/tripack.f90'])
ext3 = Extension(name    = 'stripy._srfpack',
                 sources = ['src/srfpack.pyf', 'f77-src/srfpack.f'])
ext4 = Extension(name    = 'stripy._ssrfpack',
                 sources = ['src/ssrfpack.pyf', 'f77-src/ssrfpack.f'])

if __name__ == "__main__":
    setup(name = 'stripy',
          author            = "Louis Moresi",
          author_email      = "louis.moresi@unimelb.edu.au",
          url               = "https://github.com/University-of-Melbourne-Geodynamics/stripy",
          version           = "0.3.0",
          description       = "Python interface to TRIPACK and STRIPACK fortran code for triangulation/interpolation in Cartesian coordinates and on a sphere",
          ext_modules       = [ext1, ext2, ext3, ext4],
          packages          = ['stripy'],
          classifiers       = ['Programming Language :: Python :: 2',
                               'Programming Language :: Python :: 2.6',
                               'Programming Language :: Python :: 2.7',
                               'Programming Language :: Python :: 3',
                               'Programming Language :: Python :: 3.3',
                               'Programming Language :: Python :: 3.4',
                               'Programming Language :: Python :: 3.5',
                               'Programming Language :: Python :: 3.6']
          )
