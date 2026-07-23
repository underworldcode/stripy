# Stripy changelog

## 2.3.4

 - Build against NumPy 2.x so wheels run on both NumPy 1.x and 2.x at runtime ([#114](https://github.com/underworldcode/stripy/issues/114)).
 - Replace NumPy aliases removed in NumPy 1.24 (`np.int`, `np.bool`) for NumPy 2.x compatibility.
 - Fix the release wheel workflow for `actions/upload-artifact@v4` (unique per-job artifact names).
 - Meson build robustness fixes and CI modernisation; remove the obsolete `numpy.distutils` `setup.py` build path.

## 2.1.0 

 - Jupyterbook
 - Bug fixes


## 2.0.5 (beta) 2

  - [Notebooks](https://underworldcode.github.io/stripy/2.0.5b2)
  - [API docs](https://underworldcode.github.io/stripy/2.0.5b2_api)
    - Moving the notebook examples under a jupyterbook build for the purposes of providing browseable documentation.  
    - Updating workflows to autobuild / deploy conda on geo-down-under channel 
    - Improved testing

## 2.0.0

  - Spline tensions exposed
  - Voronoi diagrams (See Notebook example #9) 

## 1.x.x

  - JOSS publication release

