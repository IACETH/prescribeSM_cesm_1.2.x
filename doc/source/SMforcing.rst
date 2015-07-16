Soil Moisture Forcing
=====================

The forcing dataset is a netCDF file with the variables ``SOILICE`` and ``SOILLIQ`` on it.
The variables are not in the usual lat/lon format (named :term:`4D`) but in an indexed format (named :term:`3D`).

The file must have 12 (365) time steps on it to prescribe SM monthly (daily) means. 

Regional Forcing
----------------

- if you set SOILLIQ or SOILICE to -1 at a gridpoint/ level, ``SOILLIQ`` and ``SOILICE`` are calculated interactively at this point!
- this can be used to prescribe SM only regionally
- or at certain depths only

.. WARNING::
   to prescribe SM only at certain time steps may not work with the current implementation (because of the temporal interpolation)



.. _3D-CLM:

Obtain 3D from CLM
------------------
CLM can directly output the needed :term:`3D` soil moisture fields.
Assume you want to output daily :term:`4D` files in ``history tape 2`` and daily :term:`3D` fields in ``history tape 3``, you have to add the following lines to :file:`user_nl_clm`::

  # daily SM, lat, lon grid:
  hist_fincl2 = 'SOILICE', 'SOILLIQ'
  # daily SM, column form
  hist_fincl3 = 'SOILLIQ','SOILICE'
  hist_nhtfrq = 0, -24, -24
  hist_mfilt = 1, 365, 365
  # conversion to lat/ lon format?
  hist_dov2xy = .true., .true., .false.


Convert 4D to 3D
----------------
If you only have your SM forcing file in :term:`4D` format it is possible to convert it to a :term:`3D` file.
A small python script that shows how to do this is given: 
:file:`clm_col_to_xy_example.py`.

It translates the lat/ lon information to index information.
This requires an :term:`3D` file with the used set up (create it as in :ref:`3D-CLM`).





