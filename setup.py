from distutils.core import setup, Extension
from Cython.Build import cythonize


ext = Extension("lib",
        sources=["lib.pyx", "cilindric.cpp","gamma.cpp","ede.cpp","file_processor.cpp",
                  "spline.cpp","tov.cpp","general_cilindric.cpp"],
        language="c++")

setup(name="ObsCalc",ext_modules=cythonize(ext))
