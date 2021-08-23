import os
from setuptools import setup, find_packages
 
classifiers = [
  'Development Status :: 5 - Production/Stable',
  'Intended Audience :: Science/Research',
  'Intended Audience :: Developers',
  'Operating System :: Microsoft :: Windows :: Windows 10',
  'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
  'Programming Language :: Python :: 3',
  'Topic :: Scientific/Engineering :: GIS',
  'Topic :: Scientific/Engineering :: Information Analysis',
  'Topic :: Scientific/Engineering :: Mathematics',
  'Topic :: Scientific/Engineering :: Visualization',
]

with open('README.md') as f:
    long_description = f.read()

# only specify install_requires if not in RTD environment
if os.getenv("READTHEDOCS") == "True":
    INSTALL_REQUIRES = []
else:
    with open("requirements.txt") as f:
        INSTALL_REQUIRES = [line.strip() for line in f.readlines()]
        
setup(
  name='x2polygons',
  version='0.0.16',
  description='A package to find the distance between two polygons',
  long_description=long_description,
  long_description_content_type='text/markdown',  
  url='https://github.com/banbar/x2polygons',  
  author='Berk Anbaroğlu',
  author_email='banbar@hacettepe.edu.tr',
  license='GNU General Public License v2.0', 
  classifiers=classifiers,
  keywords=['GIS', 'spatial analysis'], 
  packages=["x2polygons",
            "x2polygons.tests"],
  install_requires=INSTALL_REQUIRES
)

# The package may be installed to: C:\Users\banbar\anaconda3\Lib\site-packages
