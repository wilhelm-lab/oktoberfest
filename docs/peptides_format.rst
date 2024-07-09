Custom in-silico digestion
==========================

While Oktoberfest can do in-silico digestion by providing a fasta file, you can provide a list of peptides with the required metadata directly. In this case, you need to have the following parameter in your config file:

.. code-block:: json

    "library_input_type": "peptides",

Internal file format specification
----------------------------------

Oktoberfest expects a csv formatted file where each row represents a peptide fragment with given metadata. The following table provides the file format specification.

.. table::

    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | Column Header     | Explanation                                                                                                                                                                                                                                                                         |
    +===================+=====================================================================================================================================================================================================================================================================================+
    | modified_sequence | The peptide sequence including modifications in unimod format (only M[UNIMOD:35] supported) and excluding the fixed modification C[UNIMOD:4] (Carbamidomethylation) as this modification will be added automatically. If you add C[UNIMOD:4] manually, you will get wrong results.  |
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


Example
-------

.. code-block::

    modified_sequence,collision_energy,precursor_charge,fragmentation,peptide_lengths,instrument_type,proteins
    ASPTQPIQL,31,1,HCD,,,
    KIEKLKVEL,31,2,HCD,9,QE,
    AAAAAWEEPSSGNGTAR,30,3,HCD,,,Q9P258
    AAAAAWEEPSSGNGTAR,31,2,HCD,,,Q9P258
    KDVDGAYM[UNIMOD:35]TK,30,2,HCD,10,LUMOS,P04264;CON__P04264
    VIGRGSYAK,35,2,HCD,9,TIMSTOF,P11216;P11217
    TTENIPGGAEEISEVLDSLENLMR,30,1,hcd,tr|A0A075B6G3|A0A075B6G3_HUMAN;sp|P11532|DMD_HUMAN;tr|A0A5H1ZRP8|A0A5H1ZRP8_HUMAN
    TYCDATKCFTVTE,34,2,HCD,,,
