Introduction
============

prescribeSM_CESM_1.2.x allows to prescribe SOILLIQ and SOILICE in CLM 4.0 for CESM 1.2.x.

The code is originally by Ruth Lorenz (with the help of Dave Lawrence) and has be adapted by me.

Features
--------

- Prescribe SOILLIQ and SOILICE (-> SM)
- Prescribe SM globally or regionally
- Prescribe daily mean or monthly mean SM
- Two ways to do this:

  - Prescribe SOILLIQ *and* SOILICE (standard)
  - Prescribe SOILLIQ only if no SOILICE is present
  
- Define used Inputfile via namelist


Procedure
---------

- Compile CESM/ CLM with the changed 
- Prepare the :doc:`Soil Moisture forcing </SMforcing>` 
- Save namelist
- Run the model 
- 