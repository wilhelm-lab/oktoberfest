Custom modifications
====================

By default, Oktoberfest uses only static carbamidomethylation and variable methionine oxidation.
Since adding the required UNIMOD id and monoisotopic modification mass for any modification one could think of is quite cumbersome, the user can provide them manually, via the configuration file.

Important: TMT modifications must be provided using the "tag" flag in the configuration file, they are not a custom modification!


Required information
--------------------

Oktoberfest needs to map the search engine specific output to UNIMOD format to create the internal file format, as well as the monoisotopic modification mass.
The configuration file accepts two flags, "static_mods", and "var_mods", which both expect key-value mappings of the form:

"<key>": [<UNIMOD_ID>, <mod_mass>]

The UNIMOD_ID and the corresponding monoisotopic modification mass can be retrieved from `unimod.org <https://unimod.org/>`_ (click on "Login as Guest" if prompted for credentials).
Search for your desired modification, then find the UNIMOD ID in the "Accession #" column and the modification mass in the "Monoisotopic mass" column.

The format of the key depends on the search engine. Please consult the documentation for the search engine of your choice. Below are a few examples for common modifications for the search engines that Oktoberfest supports out of the box.

This will translate <key> to "[UNIMOD:<UNIMOD_ID>]" in the internal format (ProForma standard) and use <mod_mass> to add to the peptide mass.
Important: If you want to provide n-terminal or c-terminal modifications, you have to add a "^" or "$", respectively, in front of the key since Oktoberfest automatically appends or prepends the modification identifier with a "-" (ProForma standard). This also means, that "-" is part of the key, if the search engine already uses this notation!


Example configuration
---------------------

Oktoberfest supports two configuration flags, one for static and one for variable modifications, containing the search engine specific format as a key, and the UNIMOD ID and monoisotopic modification mass as a value:

.. code-block:: json

    "static_mods": {
        "^(Lactylation (Lac))": [2114, 72.021129]
        "C": [4, 57.021464]
    },
    "var_mods": {
        "M(ox)": [35, 15.994915],
        "(Deamidation (NQ))": [7, 0.984016]
    },


The above example is a MaxQuant specific mapping. It converts any carbamidomethylated "C" to "C[UNIMOD:4]", variable oxidation of "M" into "M[UNIMOD:35]", and, since no specific aminoacid was provided any deamidation on any aminoacid to "[UNIMOD:7]". Additionally, any n-terminal lactylation is converted to "[UNIMOD:2114]-" (note the added "-").


Table of default modifications
------------------------------

.. table::
   :class: fixed-table

   +-------------------------+-------------+-----------+---------------------------+-----------+-------------------+-------------------------+
   | MaxQuant                | Sage        | MSFragger | Modification              | UNIMOD ID | Monoisotopic Mass | Internal representation |
   +=========================+=============+===========+===========================+===========+===================+=========================+
   | C                       | C[+57.0215] | C[160]    | carbamidomethylation of C | 4         | 57.0215           | C[UNIMOD:4]             |
   +-------------------------+-------------+-----------+---------------------------+-----------+-------------------+-------------------------+
   | M(ox); M(Oxidation (M)) | M[+15.9949] | M[147]    | oxidation of M            | 35        | 15.9949           | M[UNIMOD:35]            |
   +-------------------------+-------------+-----------+---------------------------+-----------+-------------------+-------------------------+
   | R(Citrullination)       | R[+0.98402] | R[157]    | citrullination of R       | 7         | 0.984016          | R[UNIMOD:7]             |
   +-------------------------+-------------+-----------+---------------------------+-----------+-------------------+-------------------------+
   | Q(Deamidation (NQ))     | Q[+0.98402] | Q[129]    | deamidation of Q          | 7         | 0.984016          | R[UNIMOD:7]             |
   +-------------------------+-------------+-----------+---------------------------+-----------+-------------------+-------------------------+
   | N(Deamidation (NQ))     | N[+0.98402] | N[115]    | deamidation of N          | 7         | 0.984016          | R[UNIMOD:7]             |
   +-------------------------+-------------+-----------+---------------------------+-----------+-------------------+-------------------------+
