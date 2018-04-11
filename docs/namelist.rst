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
  path to the history file containing SOILLIQ and SOILICE. If ``one_file_per_day`` is ``.true.``, the string patterns ``%y``, ``%m``, and ``%d`` are expanded to the year, day, and month of the simulation.
one_file_per_day : bool
  If ``.false.`` uses the ``pSMfile`` as is, and iterates through the timesteps (days or months). If ``.true.`` uses one file per day and uses its first timestep.
  Default = ``.true.``.
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
use_qdrai : bool, optional
  Only for pSMtype 3. Whether subsurface runoff (qdrai) is also used for irrigation, default = ``.true.``.
reservoir_capacity, float, optional
  Only for pSMtype 3. Size of the reservoir where water can be transfered in time, default = 0.
 
.. WARNING::
   if you set ``monthly=.true.`` but have a daily SM input file it will still work (uses the first 12 days as the months)
