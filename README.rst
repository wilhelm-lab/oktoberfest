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

How to cite
===========

Please always cite the main publication
---------------------------------------

[1] Picciani M, Gabriel W, Giurcoiu VG et al. (2023),
*Oktoberfest: Open-source spectral library generation and rescoring pipeline based on Prosit*,
`Proteomics <https://doi.org/10.1002/pmic.202300112>`__.

Please cite the Prosit model(s) you are using
---------------------------------------------

**When using Prosit**

[2] Gessulat S, Schmidt T, Zolg DP et al. (2019),
*PROSIT: Proteome-wide prediction of peptide tandem mass spectra by deep learning*,
`Nature Methods <https://doi.org/10.1038/s41592-019-0426-7>`__.

[3] Wilhelm M, Zolg DP, Graber M et al. (2021),
*Deep learning boosts sensitivity of mass spectrometry-based immunopeptidomics*,
`Nature Communications <https://doi.org/10.1038/s41467-021-23713-9>`__.

**When using Prosit-TMT**

[4] Gabriel W, The M, Zolg D et al. (2022),
*Prosit-TMT: Deep Learning Boosts Identification of TMT-Labeled Peptides*,
`Analytical Chemistry <https://doi.org/10.1021/acs.analchem.1c05435>`__.

**When using Prosit-timsTOF**

[5] Adams C, Gabriel W, Laukens K et al. (2024),
*Fragment ion intensity prediction improves the identification rate of non-tryptic peptides in timsTOF*,
`Nature Communications <https://doi.org/10.1038/s41467-024-48322-0>`__.

**When using Prosit-XL**

[6] Kalhor M, Saylan C, Picciani M et al. (2025),
*Prosit-XL: enhanced cross-linked peptide identification by fragment intensity prediction to study protein interactions and structures*,
`Nature Communications <https://www.nature.com/articles/s41467-025-61203-4>`__.

**When using Prosit-MultiFrag**

[7] Levin N, Saylan CC, Lapin J et al. (2026),
*Integration of alternative fragmentation techniques into standard LC-MS workflows using a single deep learning model enhances proteome coverage*,
`Nature Methods <https://doi.org/10.1038/s41592-026-03042-9>`__.

Please cite in case you use quantification
-------------------------------------------

[8] The M, Samaras P, Kuster B, Wilhelm, M. (2022),
*Reanalysis of ProteomicsDB using an accurate, sensitive, and scalable false discovery rate estimation approach for protein groups*,
`Molecular & Cellular Proteomics <https://doi-org.org/10.1016/j.mcpro.2022.100437>`__.
