.. module:: oktoberfest

.. automodule:: oktoberfest
   :noindex:

API
===

Import Oktoberfest using

.. code-block:: python

   import oktoberfest as ok

Preprocessing: :code:`pp`
-------------------------

.. module:: oktoberfest.pp

.. currentmodule:: oktoberfest

Generating libraries

.. autosummary::
   :toctree: api/pp

   pp.digest
   pp.gen_lib
   pp.merge_spectra_and_peptides

Spectra preprocessing

.. autosummary::
   :toctree: api/pp

   pp.list_spectra
   pp.convert_spectra
   pp.load_spectra


Peptide preprocessing

.. autosummary::
   :toctree: api/pp

   pp.convert_search
   pp.load_search
   pp.filter_peptides


Predicting: :code:`pr`
----------------------

.. module:: oktoberfest.pr

.. currentmodule:: oktoberfest

All things predictions

Koina interface
~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: api/pr

   pr.grpc_predict
   pr.infer_predictions

Postprocessing koina response
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: api/pr


   pr.parse_fragment_labels


Rescoring: :code:`re`
---------------------

.. module:: oktoberfest.re

.. currentmodule:: oktoberfest

General
~~~~~~~

.. autosummary::
   :toctree: api/re

   re.generate_features
   re.merge_input
   re.rescore_with_mokapot
   re.rescore_with_percolator


Plotting: :code:`pl`
---------------------

.. module:: oktoberfest.pl

.. currentmodule:: oktoberfest

General
~~~~~~~

.. autosummary::
   :toctree: api/pl

   pl.plot_score_distribution
   pl.joint_plot
   pl.plot_gain_loss
   pl.plot_mean_sa_ce
   pl.plot_violin_sa_ce
   pl.plot_pred_rt_vs_irt
   pl.plot_all
