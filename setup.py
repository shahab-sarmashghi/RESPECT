from setuptools import setup
import sys

from respect import __version__

if sys.version_info[:2] in [(3, 7), (3, 8)]:
    install_requires = ["numpy<=1.19.5", "scipy<=1.5.0", "pandas==1.0.0", "gurobipy<=10.0.2"]
else:
    sys.exit('Sorry, this version of Python is not supported. Please use Python 3.7 or 3.8')

setup(name='respect',
      python_requires='>=3.6',
      version=__version__,
      description='Estimating repeat spectra and genome length from low-coverage genome skims',
      author='Shahab Sarmashghi',
      author_email='ssarmash@ucsd.edu',
      license='BSD-3-Clause',
      url='https://github.com/shahab-sarmashghi/RESPECT',
      packages=['respect'],
      package_dir={'respect': 'respect'},
      install_requires=install_requires,
      dependency_links=['https://pypi.gurobi.com/gurobipy/'],
      include_package_data=True,
      zip_safe=True,
      provides=["respect"],
      entry_points={
            'console_scripts': ['respect=respect.__main__:main']
      },
      classifiers=["Environment :: Console",
                   "Intended Audience :: Developers",
                   "Intended Audience :: Science/Research",
                   "Operating System :: Unix",
                   "Programming Language :: Python",
                   "Programming Language :: Python :: 3",
                   "Topic :: Scientific/Engineering :: Bio-Informatics"],
      )

