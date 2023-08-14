Custom search results
=====================

If the search engine you get your results from is not directly supported by Oktoberfest, the outputs can be manually transformed into the internal prosit file format. If you want to use this format, you need have the following parameter in your config file:

.. code-block:: json

   "search_results_type": "internal",

Internal file format specification
----------------------------------

Oktoberfest expects a csv formatted file where each row represents a PSM. The following provides the file format specification.


.. table::

    +-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------+
    | Column Header     | Explanation                                                                                                                                              |
    +===================+==========================================================================================================================================================+
    | RAW_FILE          | Name of the RAW file or mzml file associated with the PSM without the file extension                                                                     |
    +-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------+
    | SCAN_NUMBER       | RAW file derived sequential number of an individal scan during the mass spectrometry run                                                                 |
    +-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------+
    | MODIFIED_SEQUENCE | Peptide sequence including modifications in UNIMOD format                                                                                                |
    +-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------+
    | PRECURSOR_CHARGE  | Charge state of the precursor ion                                                                                                                        |
    +-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------+
    | SCAN_EVENT_NUMBER | Optional number of the specific scan event in relation to the acquisition cycle that depends on the used search engine and might not be present          |
    +-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------+
    | MASS              | Monoisotopic mass of the peptide including modifications                                                                                                 |
    +-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------+
    | SCORE             | Search engine derived score (highest score is considered "best", i.e. needs manual transformation if the search engine score is not following this rule) |
    +-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------+
    | REVERSE           | Whether the PSM is a decoy ("True") or a target ("False")                                                                                                |
    +-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------+
    | SEQUENCE          | Unmodified peptide sequence                                                                                                                              |
    +-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------+
    | PEPTIDE_LENGTH    | Length of the unmodified peptide sequence                                                                                                                |
    +-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------+

Example
-------

.. code-block::

    RAW_FILE,SCAN_NUMBER,MODIFIED_SEQUENCE,PRECURSOR_CHARGE,SCAN_EVENT_NUMBER,MASS,SCORE,REVERSE,SEQUENCE,PEPTIDE_LENGTH
    GN20170722_SK_HLA_G0103_R1_02,51825,AAAAPPWAC[UNIMOD:4]FAAV,2,2,1301.6227,3.9695,True,AAAAPPWACFAAV,13
    GN20170722_SK_HLA_G0103_R1_01,23561,AMILGKM[UNIMOD:35]VL,2,6,990.56059,47.574,False,AMILGKMVL,9
    GN20170722_SK_HLA_G0103_R2_01,30017,APAKPGGP,1,9,693.38097,2.9045,False,APAKPGGP,8
    GN20170722_SK_HLA_G0103_R2_02,27809,APIC[UNIMOD:4]SEAYSHC[UNIMOD:4]C[UNIMOD:4]DC[UNIMOD:4]F,6,10,1875.6685,0.0,True,APICSEAYSHCCDCF,15
    GN20170722_SK_HLA_G0103_R1_01,13509,KM[UNIMOD:35]PAC[UNIMOD:4]NIM[UNIMOD:35]L,2,11,1108.5079,57.047,False,KMPACNIML,9
    GN20170722_SK_HLA_G0103_R2_01,29784,KNHGRARW,3,2,1023.5475,2.6794,False,KNHGRARW,8
    GN20170722_SK_HLA_G0103_R1_02,26938,KNHKKSHK,3,5,1005.5832,5.4817,True,KNHKKSHK,8
    GN20170722_SK_HLA_G0103_R1_02,31215,KNHKKSHK,3,10,1005.5832,4.4157,True,KNHKKSHK,8
    GN20170722_SK_HLA_G0103_R2_02,25018,KNHKKSHK,3,4,1005.5832,6.8953,True,KNHKKSHK,8


