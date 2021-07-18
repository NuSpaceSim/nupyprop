.. _cross-section:

Neutrino Cross-Sections
=======================

This file contains neutrino-nucleon and anti neutrino-nucleon cross-section models. They are stored in lookup_tables.h5 under path. The user is able to choose from the default models included with the code package or define/import their own cross-section model. The imported/defined model needs to be binned in a 91 binned energy logspace ie. np.logspace(3,12,91) from E=10^3 GeV to 10^12 GeV. They will need to have a charged current cross-section and a neutral current cross-section value in cm^2.

Provide example table.