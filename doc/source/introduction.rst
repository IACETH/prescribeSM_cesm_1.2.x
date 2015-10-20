Introduction
============

``prescribeSM_CESM_1.2.x`` allows to prescribe ``SOILLIQ`` and ``SOILICE`` in CLM 4.0 for CESM 1.2.x.

The code is originally by Ruth Lorenz (with the help of Dave Lawrence) and has been adapted by me for CESM 1.2.1 and CESM 1.2.2.

.. NOTE::
      The code of ``clm_driver`` and ``BalanceCheckMod`` are identical in 1.2.1 and 1.2.2.

Features
--------

- Prescribe ``SOILLIQ`` and ``SOILICE`` (-> SM) in CLM
- Prescribe SM globally or regionally
- Prescribe daily mean or monthly mean SM (linear interpolation in between)
- Two ways to do this:

  - Prescribe ``SOILLIQ`` *and* ``SOILICE`` (standard)
  - Prescribe ``SOILLIQ`` only if no ``SOILICE`` is present
  
- Define used Inputfile via namelist

Procedure
---------

- :doc:`Compile CESM/ CLM</compile>` with the changed source files
- Prepare the :doc:`Soil Moisture forcing </SMforcing>` 
- Create the :doc:`namelist` and put it in the `run` directory
- Run the model 
