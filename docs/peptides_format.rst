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
    | modified_sequence | The peptide sequence including modifications in unimod format (only M[UNIMOD:35] supported) and excluding the fixed modification C[UNIMOD:35] (Carbamidomethylation) as this modification will be added automatically. If you add C[UNIMOD:35] manually, you will get wrong results.|
    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | collision_energy  | The collision energy to use in peptide property prediction                                                                                                                                                                                                                          |
    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | precursor_charge  | Charge state of the precursor ion                                                                                                                                                                                                                                                   |
    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | fragmentation     | Method used for fragmentation; can be "HCD" or "CID"                                                                                                                                                                                                                                |
    +-------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


Example
-------

.. code-block::

    modified_sequence,collision_energy,precursor_charge,fragmentation
    ASPTQPIQL,35,1,HCD
    KIIDRAITSL,34,2,HCD
    KIEKLKVEL,35,2,HCD
    KINQQKLKL,34,3,HCD
    MLGNM[UNIMOD:35]NVFMAVLGIILC[UNIMOD:4]SGFLAAYFSHK,30,4,HCD
    TYC[UNIMOD:4]DATKC[UNIMOD:4]FTVTE,34,2,HCD
    VIPSIAYTEPEVAWVGLTEKEAKEK,30,2,HCD
