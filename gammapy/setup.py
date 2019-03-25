from setuptools import setup, find_packages
import setuptools





setup(name='gammapy',
	  license='Creative Commons Attribution-Noncommercial-Share Alike license',
	  author= 'David Drumond',
	  version ='0.1',
	  author_email= 'david.engminas.ufmg@gmail.com',
	  description ='Gammaypy is a computer development to calculate spatial continuity functions,\
	  				such as variograms, covariograms, madograms, etc.',
	  packages=setuptools.find_packages(),
	  url="https://github.com/DavidDrumond/gammapy-",
	  long_description = open('README.txt').read(),
	  zip_safe = False)

