"""pmx: a toolkit for free-energy calculation setup/analysis and
biomolecular structure handling.

For installation type the command::
  python setup.py install
or
  pip install .
"""

from setuptools import setup, Extension


def readme():
    with open('README.md') as f:
        return f.read()

pmx = Extension('pmx/_pmx',
                libraries=['m'],
                include_dirs=['src/pmx'],
                sources=['src/pmx/Geometry.c', 'src/pmx/wrap_Geometry.c',
                         'src/pmx/init.c', 'src/pmx/Energy.c']
                )

xdrio = Extension('pmx/_xdrio',
                  libraries=['m'],
                  include_dirs=['src/xdr'],
                  sources=['src/xdr/xdrfile.c', 'src/xdr/xdrfile_trr.c',
                           'src/xdr/xdrfile_xtc.c']
                  )

setup(name='pmx',
      version='1.1.0dev',
      description='Toolkit for free-energy calculation setup/analysis '
                  'and biomolecular structure handling',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
      ],
      url='https://github.com/dseeliger/pmx/',
      author='Daniel Seeliger',
      author_email='seeliger.biosoft@gmail.de',
      license='GPL 3',
      packages=['pmx'],
      package_dir={'pmx': 'pmx'},
      include_package_data=True,
      package_data={'pmx': ['data/*', 'data/*/*', 'data/*/*/*']},
      zip_safe=False,
      ext_modules=[pmx, xdrio],
      install_requires=['numpy', 'scipy', 'matplotlib']
      )
