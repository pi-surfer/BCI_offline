# Brain Computer Interface Simulation Offline in Matlab

Offline simulation of a BCI loop using data collected during a 3-day experiment with 8 healthy subjects. Data has been recorded with 16-channel EEG amplifier (g.USBamp, g.Tec) at 512 Hz. Electrodes were placed accordingly to the 10-20 international layout.  Each subject participated in at least 2 recording days. In the first day, subjects performed 3 “offline” (calibration, no real feedback) and 2 “online” (with real feedback) runs. In the second day and third day, they performed 2 “online” runs each.

This project is part of the "Neurorobotics and Neurorehabilitation" class, Bioengineering Master, University of Padua.
Authors: Andrea Rachele Aparo, Virginia Furlani, Francesco Pio Monaco, Anna Pegreffi, Lorenzo Sterzi.

## Methodology

Two types of analyses were done:

  1. Grand average analyses on the whole population and on representative subjects

  2. Analyses on BMI decoding on each subject
       - Calibration phase:
         - Data processing of offline runs;
         - Feature extraction and selection of the most discriminant features;
         -  Creation of a classifier based on those features.

       - Evaluation phase:
         - Evaluation of the classifier created during the calibration phase using data from online runs;
         - Implementation of a evidence accumulation framework on the posterior probabilities.
