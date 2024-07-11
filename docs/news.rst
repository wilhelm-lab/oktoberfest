News
====

.. role:: date
    :class: date


Oktoberfest provides access to AlphaPept and MS2PIP predictions :date:`2024-05-29`
----------------------------------------------------------------------------------

| Oktoberfest 0.7.0 is published, providing access to predictions from AlphaPept and MS2PIP. This comes with a large API overhaul that allows easier access to underlying data for manual insights and analysis.
| Please note that you can now specify the instrument type in the configuration for AlphaPept predictions.
| Stay tuned for a documentation update including jupyter notebooks on how to navigate and make use of your data for manual in-depth analysis.
| This release drops python 3.8 support!

Support for timsTOF added / spectral library generation overhaul :date:`2024-01-31`
------------------------------------------------------------------------------------

| Oktoberfest 0.6.0 is out supporting the new Prosit model for timsTOF data. Check out the preprint at `BioRxiv <https://doi.org/10.1101/2023.07.17.549401>`__.
| Spectral library generation was completely overhauled, resulting in more stability, dramatic decreases in runtime and storage requirements and more customization in the config.
| Note that this release introduces changes in the configuration file that might break your previous setup! Please check the `new configuration format <./config.html>`_.

Oktoberfest now provides a low level API :date:`2023-10-03`
------------------------------------------------------------

| With our new release 0.5.0, we refactored our code base to provide a low level API to the users. Check the new `API <./API.html>`_ section for more details.
| Please note: We dropped support for percolator version 3.05. Please update to version 3.06.1. See `installation <./installation.html>`_ for details.


Oktoberfest is published! :date:`2023-09-07`
---------------------------------------------

| Alongside the new 0.4.0 release, Oktoberfest was published in Proteomics.
| Check out the publication here: `https://doi.org/10.1002/pmic.202300112 <https://doi.org/10.1002/pmic.202300112>`_.
| Accordingly, please check `How to cite <./reference.html>`_.


The first Oktoberfest release is ready! :date:`2023-06-25`
-----------------------------------------------------------

Version 0.3.0 is available on `PyPI <https://pypi.org/project/oktoberfest/>`_.

We are happy to announce that our first standalone version of Oktoberfest is finally ready.
Besides a lot of bugfixes, new features and performance enhancements over the existing online Prosit service at `ProteomicsDB <https://proteomicsdb.org/prosit>`_, the Oktoberfest codebase is now open source and available on `GitHub <https://github.com/wilhelm-lab/oktoberfest>`_.
In addition, we switched to our new online prediction service `koina <https://koina.proteomicsdb.org>`_ which allows connecting Oktoberfest to any public or self-hosted koina instance for easy prediction model integration and being independent of ProteomicsDB!