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
monthly : bool
  If ``.true.`` assumes monthly mean SM is provided. If ``.false.`` uses daily input.
pSMtype : int, {1, 2, 3}
 Defines how to prescribe SM:

 1. Default: prescribe SOILLIQ and SOILICE.
 2. Prescribe SM only if there is no SOILICE is present (never prescribes SOILICE). Lets the model decide if there is ICE or not at a certain gridpoint.
 3. Use runoff (QOVER) to prescribe SM. Is currently work in progress.
 
.. WARNING::
   if you set monthly=.true. but have a daily SM input file it will still work (uses the first 12 days as the months)
