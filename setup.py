import os.path
import sys
from setuptools import find_packages
from distutils.core import setup


here = os.path.abspath(os.path.dirname(__file__))
exec(open(os.path.join(here, 'pyconserve/version.py')).read())

setup(name='pyconserve',
      version=__version__,
      description='Extract conservation for given BED file',
      url='http://github.com/kcha/pyconserve',
      author='Kevin Ha',
      author_email='k.ha@mail.utoronto.ca',
      license='MIT',
      packages=find_packages(),
      install_requires=['setuptools',
                        'pandas >= 0.24',
                        'pybedtools >= 0.7.9'],
      python_requires='~=3.6',
      entry_points={
            'console_scripts': [
                'pyconserve = pyconserve.pyconserve:main',
                'summarize_conserve = pyconserve.summarize:main'
            ]
      }
)
