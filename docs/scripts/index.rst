.. _scripts:

Scripts
=======

pmx provides a few command line scripts that can be used to setup and analyse
free energy calculations. After installing the library, ``pmx``
should be in your  ``$PATH``. You can check this with the following command. ::

  $ which pmx

The ``pmx`` help will show the scripts available. ::

   $ pmx -h
   usage: pmx [-h]

    ------------------------
    pmx command line scripts
    ------------------------

    Available commands are:
        mutate     Mutate protein or DNA/RNA
        analyse    Estimate free energy from Gromacs xvg files

    optional arguments:
        -h, --help  show this help message and exit

A description of these scripts can be found here:

.. toctree::
   :maxdepth: 1

   mutate
   analysis

In the :ref:`examples` you can find instead how the same tasks can be carried out
using the API.
