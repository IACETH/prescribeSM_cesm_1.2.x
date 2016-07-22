Soil Moisture Data Set
======================

This file defines the target values which are prescribed in CLM. It is a netCDF file and needs to have the variables ``SOILICE`` and ``SOILLIQ`` on it (even if you only prescribe ``SOILLIQ``).
The variables are not in the usual lat/ lon format (named :term:`4D`) but in an indexed format (named :term:`3D`). The data is obtained from a reference simulation of CLM. 
The file must have 365 (12) time steps on it to prescribe SM daily (monthly) data. 

.. _3D-CLM:

Obtain SM Data Set from CLM
---------------------------
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
If you only have the SM forcing file in :term:`4D` format it is possible to convert it to a :term:`3D` file.
A small python script that shows how to do this is given in the repository: 
:file:`clm_col_to_xy_example.py`.

It translates the lat/ lon information to index information.
This requires one :term:`3D` file with the used set up (create it as in :ref:`3D-CLM`).


Regional Forcing
----------------

- if you set SOILLIQ or SOILICE to -1 at a gridpoint/ level, ``SOILLIQ`` and ``SOILICE`` are calculated interactively at this point
- this can be used to prescribe SM only regionally
- or at certain depths only

.. WARNING::
   to prescribe SM only at certain time days of the year will only work if you use (1) daily data and (2) set 
