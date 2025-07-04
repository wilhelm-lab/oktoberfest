.. highlight:: shell

Installation
============

Oktoberfest can be installed on all three major platforms (Linux, MacOS, Windows).

An installer script is provided for Linux (installers for MacOS and Windows are currently prepared), which takes care of all dependencies. For MacOS and Windows, the Docker container or manual installation can be used.

Installer script (Debian / Ubuntu only)
---------------------------------------

The installer script automatically installs dependencies and creates a new conda environment for oktoberfest. The installation takes roughly 5 minutes and follows the steps outlined in the manual installation below. Get the installer and execute the script using

.. code-block:: bash

   wget https://raw.githubusercontent.com/wilhelm-lab/oktoberfest/main/installer.sh -O install_oktoberfest.sh
   bash install_oktoberfest.sh

The installer searches for an existing anaconda / miniconda installation. If none is found, it will download and install miniconda.

Docker Image
------------

This is only recommended for experienced users and can take up to 30 minutes to install, depending on experience.

Prerequisites:
  - `make <https://www.gnu.org/software/make/>`_
  - `docker <https://www.docker.com/>`_

After cloning the repository of oktoberfest, checkout the branch you want to build the container from.
The latest stable version is always on the main branch. Then build the container using:

.. code-block:: bash

   make build


Manual installation
-------------------

This is a step-by-step guide for the manual installation of all mandatory and optional dependencies which should roughly take up to 20 minutes if all dependencies are installed.

Install Python
~~~~~~~~~~~~~~

Oktoberfest requires python >=3.9 and <=3.11. Best practice is to use a clean conda environment (`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_).
Follow the installation guide for your operating system, then create a new environment using

.. code-block:: bash

   conda create -y -n oktoberfest python==3.10


Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

There are multiple optional dependencies depending on job types. Detailed notes and installation instructions can be found below.

.. table::

    +------------------------------------------------------------------+---------------------+---------------------------------------------------------+
    | Job type                                                         | Dependency          | Notes                                                   |
    +==================================================================+=====================+=========================================================+
    | Pre-processing                                                   | ThermoRawFileParser | not required if only mzML files are provided            |
    +                                                                  +---------------------+---------------------------------------------------------+
    |                                                                  | `mono`              | required for ThermoRawFileParser to work on Linux/macOS |
    +------------------------------------------------------------------+---------------------+---------------------------------------------------------+
    | :ref:`Rescoring <jobs:a) without refinement>`                    | Percolator          |                                                         |
    +------------------------------------------------------------------+---------------------+---------------------------------------------------------+
    | :ref:`Rescoring + refinement learning <jobs:b) with refinement>` | DLomix              |                                                         |
    +------------------------------------------------------------------+---------------------+---------------------------------------------------------+

**ThermoRawFileParser**
`ThermoRawFileParser v1.4.3 <https://github.com/compomics/ThermoRawFileParser/releases/tag/v1.4.3>`_:
For conversion of RAW to mzML format. Download and unpack the zip or tar.gz file. The default locations Oktoberfest expects the executable to be at `/opt/compomics/` (Linux/MacOS) or the folder from which you want to execute Oktoberfest (Windows).
You do not need this package if you only ever provide mzML files. However, it is recommended to let Oktoberfest convert RAW files for you, to ensure the mzML files are formatted in the way Oktoberfest expects it.

**`mono`**
For ThermoRawFileParser to work on Linux, you also need to ensure `mono` is installed using

.. code-block:: bash

   sudo apt -y update && sudo apt -y install mono-devel  # Debian / Ubuntu

For MacOS, follow the instructions provided by `Mono <https://www.mono-project.com/docs/getting-started/install/mac/>`_.

**Percolator**
`Percolator v3.06.1 <https://github.com/percolator/percolator/releases/tag/rel-3-06-01>`_:
This is the tool Mokapot is based on. As it has more options and is generally more stable wrt. to FDR cutoffs and deduplication, it is recommended to use this tool instead of Mokapot.
Installable packages are provided for Linux/MacOS/Windows.

**DLomix**
`DLomix <https://github.com/wilhelm-lab/dlomix>`_ is a Python framework for deep learning in proteomics. Oktoberfest uses DLomix to refinement-learn intensity predictors on input spectra. It is listed as an optional dependency and can be installed using

.. code-block:: bash

    poetry install -E dlomix

Installing Oktoberfest
~~~~~~~~~~~~~~~~~~~~~~

Oktoberfest is listed on the Python Package Index (PyPI) and can be installed with pip. Activate your conda environment (or skip this if you use a system wide python installation) and install the package (and optionally jupyterlab) using

.. code-block:: bash

   conda activate oktoberfest
   pip install oktoberfest jupyterlab

For local prediction & refinement learning, you have to install Oktoberfest with the `dlomix` extra:

.. code-block:: bash

   conda activate oktoberfest
   pip install oktoberfest[dlomix]

Installing xiFDR
~~~~~~~~~~~~~~~~

Oktoberfest generates input files for **xiFDR**, which is required for FDR estimation of XL-MS data. It does not run xiFDR directly.

To install and run xiFDR, please follow the instructions at:

- https://www.rappsilberlab.org/software/xifdr/
- https://github.com/Rappsilber-Laboratory/xiFDR
