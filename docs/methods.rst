Methods to Prescribe Soil Moisture
==================================

pSMtype = 1
  Prescribe ``SOILLIQ`` and ``SOILICE`` at every time step (PRES_LIQ+ICE). Ignores gridcells/ levels where ``SOILLIQ`` or ``SOILICE`` is negative.

pSMtype = 2
  Prescribe ``SOILLIQ`` if soil temperature is above freezing (PRES_LIQ). Prescribes ``SOILLIQ + SOILICE`` at every time step. Ignores gridcells/ levels where ``SOILLIQ`` is negative.

pSMtype = 3
  Not tested, do not use. Prescribe ``SOILLIQ`` if soil temperature is above freezing (PRES_LIQ) and only if runoff (or potentially water in the reservoir is available). Does not remove water. Prescribes ``SOILLIQ + SOILICE`` at every time step.

pSMtype = 4
  As pSMtype = 2 but does not remove water. 


