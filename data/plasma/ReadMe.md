# Data origin

Bian, Y., Bayer, F. P., Chang, Y. C., Meng, C., Hoefer, S., Deng, N., ... & Kuster, B. (2021). Robust microflow LC-MS/MS for proteome analysis: 38 000 runs and counting. Analytical Chemistry, 93(8), 3686-3690.

https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/02/PXD023650/

# Running this example

First, you will need to download the raw file "Plamsa_03401_RD2_30min_Top10_86ms_RP13_R1" from the above resource and store it in this directory.

You can then run this from the base directory using

```
python tests/integration_tests/test_re_score.py
```

Note: If ThermoRawFileParser.exe is not residing in the base directory, you need to alter the thermoExe parameter in the config to point to the executable.

In case you have installed oktoberfest in a docker container, please execute

```
DATA=data/plasma make run_oktoberfest
```

from the base directory (ThermoRawFileParser.exe is installed in the base directory automatically).
