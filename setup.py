from numpy.distutils.core  import setup, Extension

# interface for Renka's algorithm 772 fortran code
ext = Extension(name  = '_stripy',
                sources       = ['_stripy.pyf','_stripy.f90'])

if __name__ == "__main__":
    setup(name = 'stripy',
          author            = "LM",
          author_email      = "louis.moresi@unimelb.edu.au",
          url               = "",
          download_url      = "",
          version           = "0.1",
          description       = "Python interface to STRIPACK fortran code for triangulation/interpolation on a sphere and a few more things",
          ext_modules       = [ext],
          packages          = ['stripy'],
          )
