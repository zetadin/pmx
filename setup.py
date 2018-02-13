# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2013 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

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
      description='Toolkit for free-energy calculation setup/analysis and biomolecular structure handling',
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
