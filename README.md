# SphOut
A Multinest-based outflow model fitter for observed spectra

Requirements:
  -numpy
  -matplotlib
  -scipy
  -astropy
  -uncertainties
  -Multinest (https://github.com/farhanferoz/MultiNest)
  -PyMultinest  (https://johannesbuchner.github.io/PyMultiNest/)
  
SphOut is made of two parts: a code that takes in various geometrical parameters of a galactic outflow and calculates its expected observed spectrum, and a code that uses Multinest to fit these models to an observed spectrum.

Last Updated: 6 August, 2019
