Namelist
========

PrescribeSM has its own namelist that is called :file:`prescribe_SM_nl`.
This namelist must go in the :file:`run` directory of CESM.
The location of the :file:`run` directory is defined in :file:`$CASEROOT/env_run.xml` under the entry "RUNDIR".

.. NOTE::
   an example namelist is given with the code

Namelist Options
----------------
pSMfile : string
  path to the history file containing SOILLIQ and SOILICE
pSMtype : int
 Defines how to prescribe SM, see :doc:`Methods to Prescribe Soil Moisture</methods>`.
monthly : bool
  If ``.false.`` uses daily input. If ``.true.`` assumes monthly mean SM is provided.
interp_day : bool
  If ``.false.`` uses the daily mean, if ``.true.`` linearly interpolates between daily mean. Applies only if daily input data is used (``monthly = .false.``).
levstart : int, optional
  First level to prescribe SM, default = 1.
levstop : int, optional
  Last level to prescribe SM, default = 10.
nudge : float between 0. and 1., optional
  Nudging parameter. Not tested. Default = 1 (i.e. no nudging).
use_qdrai : bool, optional
  Only for pSMtype XX. Do not use, default = ``.true.``.
reservoir_capacity, float, optional
  Only for pSMtype XX. Do not use, default = 0.
 
.. WARNING::
   if you set monthly=.true. but have a daily SM input file it will still work (uses the first 12 days as the months)
