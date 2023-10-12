.. highlight:: shell

Installation
============

Prerequisites
~~~~~~~~~~~~~

Oktoberfest requires python >=3.8,<=3.11. Best practise is to use a clean conda environment (`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_).

If you provide RAW files, you need ThermoRawFileParser for conversion to mzML.
Please download the latest release from the `github repository <https://github.com/compomics/ThermoRawFileParser>`_ using the provided zip file and unpack it in the desired location.
On linux and MacOS, the default location is "/opt/compomics/". On Windows, the default location is the directory from which Oktoberfest is executed.
You can provide the location of the executable in the config file when starting an Oktoberfest run.

To make ThermoRawFileParser work on linux or MacOS, make sure mono `mono <https://www.mono-project.com/>`_ is installed.

If you want to use percolator for rescoring, make sure you install version 3.06.1 (`percolator <https://github.com/percolator/percolator/releases/tag/rel-3-06-01>`_).

Using pip (recommended)
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install oktoberfest

Docker image
~~~~~~~~~~~~

Prerequisites:
  - `make <https://www.gnu.org/software/make/>`_
  - `docker <https://www.docker.com/>`_

After cloning the repository of oktoberfest, checkout the branch you want to build the container from.
The latest stable version is always on the main branch. Then build the container using:

.. code-block:: bash

   make build

