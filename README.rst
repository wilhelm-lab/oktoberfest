Oktoberfest: Rescoring, Collision Energy Calibration and Spectral Library Generation for Proteomics
===================================================================================================

|PyPI| |Python Version| |License| |Read the Docs| |CI| |Codecov| |pre-commit| |Ruff|

.. |PyPI| image:: https://img.shields.io/pypi/v/oktoberfest.svg
   :target: https://pypi.org/project/oktoberfest/
   :alt: PyPI
.. |Python Version| image:: https://img.shields.io/pypi/pyversions/oktoberfest
   :target: https://pypi.org/project/oktoberfest
   :alt: Python Version
.. |License| image:: https://img.shields.io/github/license/wilhelm-lab/oktoberfest
   :target: https://opensource.org/licenses/MIT
   :alt: License
.. |Read the Docs| image:: https://img.shields.io/readthedocs/oktoberfest/latest.svg?label=Read%20the%20Docs
   :target: https://oktoberfest.readthedocs.io/
   :alt: Read the documentation at https://oktoberfest.readthedocs.io/
.. |CI| image:: https://github.com/wilhelm-lab/oktoberfest/workflows/CI/badge.svg
   :target: https://github.com/wilhelm-lab/oktoberfest/actions?workflow=CI
   :alt: CI Status
.. |Codecov| image:: https://codecov.io/gh/wilhelm-lab/oktoberfest/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/wilhelm-lab/oktoberfest
   :alt: Codecov
.. |pre-commit| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
   :alt: pre-commit
.. |Ruff| image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
   :target: https://github.com/astral-sh/ruff
   :alt: Ruff

Oktoberfest is a python tool for collision energy calibration, rescoring search results and generating spectral libraries for proteomics research within the Prosit ecosystem. It offers an end to end pipeline that takes search results, predicts peptide properties using koina, plots summaries and quality control figures and performs FDR estimation with either mokapot or percolator.

Documentation
==============

The official Oktoberfest documentation can be found at `https://oktoberfest.readthedocs.io <https://oktoberfest.readthedocs.io>`__.

How to cite
===========

Please always cite the main publication
----

[Oktoberfest] Picciani M, Gabriel W, Giurcoiu VG et al. (2023), _Oktoberfest: Open-source spectral library generation and rescoring pipeline based on Prosit_, [Proteomics](https://doi.org/10.1002/pmic.202300112)

Please cite the Prosit model(s) you are using
----

When using Prosit

[Prosit] Gessulat S, Schmidt T, Zolg DP et al. (2019), _PROSIT: Proteome-wide prediction of peptide tandem mass spectra by deep learning_, [Nature Methods](https://doi.org/10.1038/s41592-019-0426-7)

[Prosit-HLA] Wilhelm M, Zolg DP, Graber M et al. (2021), _Deep learning boosts sensitivity of mass spectrometry-based immunopeptidomics_, [Nature Communications](https://doi.org/10.1038/s41467-021-23713-9)

When using Prosit-TMT

[Prosit-TMT] Gabriel W, The M, Zolg D et al. (2022), _Prosit-TMT: Deep Learning Boosts Identification of TMT-Labeled Peptides_, [Analytical Chemistry](https://doi.org/10.1021/acs.analchem.1c05435)

When using Prosit-timsTOF

[Prosit-timsTOF] Adams C, Gabriel W, Laukens K et al. (2024), _Fragment ion intensity prediction improves the identification rate of non-tryptic peptides in timsTOF_, [Nature Communications](https://doi.org/10.1038/s41467-024-48322-0)

When using Prosit-XL

[Prosit-XL] Kalhor M, Saylan CC, Picciani M et al. (2024), _Prosit-XL: enhanced cross-linked peptide identification by accurate fragment intensity prediction to study protein-protein interactions and protein structures_, [Nature Communications](https://doi.org/10.1038/s41467-025-61203-4)

When using Prosit-MultiFrag

[Prosit-MultiFrag] Levin N, Saylan CC, Lapin J (2026) _Integration of alternative fragmentation techniques into standard LC-MS workflows using a single deep learning model enhances proteome coverage_, [Nature Methods](https://doi.org/10.1038/s41592-026-03042-9)

Please cite when using protein grouping and quantification
----

[PickedGroupFDR] The M, Samaras P, Kuster B, Wilhelm, M. (2022), _Reanalysis of ProteomicsDB using an accurate, sensitive, and scalable false discovery rate estimation approach for protein groups_, [Molecular & Cellular Proteomics](https://doi-org.org/10.1016/j.mcpro.2022.100437)

