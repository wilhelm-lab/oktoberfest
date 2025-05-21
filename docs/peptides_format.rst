Custom in-silico digestion
==========================

While Oktoberfest can do in-silico digestion by providing a fasta file, you can also provide a list of peptides yourself, or follow the below internal format for the highest level of customization.

Providing a list of peptides and associated proteins
----------------------------------------------------

In this case, you need to have the following parameter in your config file:

.. code-block:: json

    "library_input_type": "peptides",

Oktoberfest will then create the table of peptides with associated metadata in internal format (see below) based on the configuration in the spectralLibraryOptions of your configuration file. For a list of these options, check the `configuration options <./config.html>`_.

Description of peptide list columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Oktoberfest expects a csv formatted file, where each row represent a peptide and optional mappings to proteins.

.. table::

    +-------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | Column Header     | Explanation                                                                                                                                                                                                  |
    +===================+==============================================================================================================================================================================================================+
    | peptide           | The unmodified peptide sequence. "C" will always be carbamidomethylated (fixed modification), and a TMT modification is always added to the N-term and "K" if a tag is specified in the configuration file.  |
    +-------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | proteins          | An optional list of protein ids separated by ';'. If this column is left out, or if no protein is provided, the string "unknown" will be used as a proteinID in the spectral library.                        |
    +-------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Example of peptide list
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::

    peptide,proteins
    ASPTQPIQL,
    KIEKLKVEL,
    AAAAAWEEPSSGNGTAR,Q9P258
    KDVDGAYMTK,P04264;CON__P04264
    VIGRGSYAK,P11216;P11217
    TTENIPGGAEEISEVLDSLENLMR,tr|A0A075B6G3|A0A075B6G3_HUMAN;sp|P11532|DMD_HUMAN;tr|A0A5H1ZRP8|A0A5H1ZRP8_HUMAN
    TYCDATKCFTVTE

Internal file format specification
----------------------------------

If you want to have full control, you can provide the table in internal format directly. In this case, you need to have the following parameter in your config file:

.. code-block:: json

    "library_input_type": "internal",

Oktoberfest will then read the table directly.

Description of internal file columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Oktoberfest expects a csv formatted file where each row represents a peptide with given metadata. The following table provides the file format specification.

.. table::

    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | Column Header     | Explanation                                                                                                                                                                                                                                                                         |
    +===================+=====================================================================================================================================================================================================================================================================================+
    | modified_sequence | The peptide sequence including variable modifications in unimod format (only M[UNIMOD:35] is supported). "C" will always be carbamidomethylated (fixed modification), and a TMT modification is always added to the N-term and "K" if a tag is specified in the configuration file. |
    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | collision_energy  | The collision energy to use in peptide property prediction                                                                                                                                                                                                                          |
    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | precursor_charge  | Charge state of the precursor ion                                                                                                                                                                                                                                                   |
    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | fragmentation     | Method used for fragmentation; can be "HCD" or "CID"                                                                                                                                                                                                                                |
    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | peptide_length    | An optional column containing the sequence list. Needed only when predicting intensities with AlphaPept.                                                                                                                                                                            |
    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | instrument_types  | An optional column containing the type of mass spectrometer. Only needed when predicting intensities with AlphaPept. Choose one of ["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"].                                                                                                          |
    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | proteins          | An optional list of protein ids separated by ';'                                                                                                                                                                                                                                    |
    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Example of internal file
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::

    modified_sequence,collision_energy,precursor_charge,fragmentation,peptide_length,instrument_types,proteins
    ASPTQPIQL,31,1,HCD,,,
    KIEKLKVEL,31,2,HCD,9,QE,
    AAAAAWEEPSSGNGTAR,30,3,HCD,,,Q9P258
    AAAAAWEEPSSGNGTAR,31,2,HCD,,,Q9P258
    KDVDGAYM[UNIMOD:35]TK,30,2,HCD,10,LUMOS,P04264;CON__P04264
    VIGRGSYAK,35,2,HCD,9,TIMSTOF,P11216;P11217
    TTENIPGGAEEISEVLDSLENLMR,30,1,hcd,tr|A0A075B6G3|A0A075B6G3_HUMAN;sp|P11532|DMD_HUMAN;tr|A0A5H1ZRP8|A0A5H1ZRP8_HUMAN
    TYCDATKCFTVTE,34,2,HCD,,,
