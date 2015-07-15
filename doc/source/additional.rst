Additional Info
===============

.. _source_files:

Source Files
------------
The tree changed source files are:

- :file:`prescribeSoilMoistureMod.F90`
- :file:`clm_driver.F90`
- :file:`BalanceCheckMod.F90`

.. object:: prescribeSoilMoistureMod.F90

The main file for prescribing soil moisture. Contains the following routines:

- prescribeSoilMoisture
- initPrescribeSoilMoisture
- interpSoilMoisture
- readSoilMoisture

.. object:: clm_driver.F90

The subroutine prescribeSoilMoisture is called from the main program of CLM.

.. object:: BalanceCheckMod.F90

  Prescribing soil moisture violates the water and energy balance.
  Therefore these two checks have to be turned off.

Glossary
--------
.. glossary::

   4D
      The usual file.

      time x level x latitude x longitude
      
   
   3D
      time x level x index



