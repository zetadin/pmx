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
import versioneer


def readme():
    with open('README.md') as f:
        return f.read()


# ----------
# Extensions
# ----------
pmx = Extension('pmx._pmx',
                libraries=['m'],
                include_dirs=['src/pmx'],
                sources=['src/pmx/Geometry.c',
                         'src/pmx/wrap_Geometry.c',
                         'src/pmx/init.c',
                         'src/pmx/Energy.c']
                )

xdrio = Extension('pmx._xdrio',
                  libraries=['m'],
                  include_dirs=['src/xdr'],
                  sources=['src/xdr/xdrfile.c',
                           'src/xdr/xdrfile_trr.c',
                           'src/xdr/xdrfile_xtc.c']
                  )
extensions = [pmx, xdrio]

# -----
# Setup
# -----
setup(name='pmx',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Python Toolbox structure file editing and writing simulation setup/analysis tools',
      long_description=readme(),
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Programming Language :: Python',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
                  ],
      author='Daniel Seeliger',
      author_email='seeliger.biosoft@gmail.de',
      url='https://github.com/deGrootLab/pmx',
      license='GPL 3',
      packages=['pmx'],
      include_package_data=True,
      package_data={'data': ['bbdep.pkl',
                             'bp.pkl',
                             'ffamber99sb.rtp',
                             'ffamber99sbbon.itp',
                             'ffamber99sbnb.itp',
                             'blosum62_new.mat'],
                    },
      zip_safe=False,
      ext_modules=extensions,
      python_requires=">=2.7, <3",
      install_requires=['numpy>=1.14', 'scipy>=1.1', 'matplotlib>=2.2']
      )
