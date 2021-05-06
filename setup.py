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

setup(
  name='x2polygons',
  version='0.0.4',
  description='A package to find the distance between two polygons',
  long_description=long_description,
  long_description_content_type='text/markdown',  
  url='https://github.com/banbar/polygon2polygon-distance',  
  author='Berk Anbaroğlu',
  author_email='banbar@hacettepe.edu.tr',
  license='GNU General Public License v2.0', 
  classifiers=classifiers,
  keywords=['GIS', 'spatial analysis'], 
  packages=find_packages(),
  install_requires=['matplotlib>=3.2.2', 
                    'shapely>=1.7.1']
)
