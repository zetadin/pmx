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

# ----------
# Extensions
# ----------
pmx = Extension('pmx._pmx',
                libraries=['m'],
                include_dirs=['pmx/extensions/pmx'],
                sources=['pmx/extensions/pmx/Geometry.c',
                         'pmx/extensions/pmx/wrap_Geometry.c',
                         'pmx/extensions/pmx/init.c',
                         'pmx/extensions/pmx/Energy.c']
                )

xdrio = Extension('pmx._xdrio',
                  libraries=['m'],
                  include_dirs=['pmx/extensions/xdr'],
                  sources=['pmx/extensions/xdr/xdrfile.c',
                           'pmx/extensions/xdr/xdrfile_trr.c',
                           'pmx/extensions/xdr/xdrfile_xtc.c']
                  )
extensions = [pmx, xdrio]

# -----
# Setup
# -----
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
      ext_modules=extensions,
      install_requires=['numpy', 'scipy', 'matplotlib']
      )
