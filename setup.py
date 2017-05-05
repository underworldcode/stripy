from numpy.distutils.core import setup, Extension

# interface for Renka's algorithm 772 fortran code
ext1 = Extension(name    = '_stripack',
                 sources = ['src/stripack.pyf','src/stripack.f90'])
ext2 = Extension(name    = '_tripack',
                 sources = ['src/tripack.pyf', 'src/tripack.f90'])
ext3 = Extension(name    = '_srfpack',
                 sources = ['src/srfpack.pyf', 'src/srfpack.f'])
ext4 = Extension(name    = '_ssrfpack',
                 sources = ['src/ssrfpack.pyf', 'src/ssrfpack.f'])

if __name__ == "__main__":
    setup(name = 'stripy',
          author            = "LM",
          author_email      = "louis.moresi@unimelb.edu.au",
          url               = "https://github.com/University-of-Melbourne-Geodynamics/stripy",
          version           = "0.1.1",
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
