Soil Moisture Forcing
=====================

The forcing dataset is a netCDF file with the variables `SOILICE` and `SOILLIQ` on it.
The variables are not in the usual lat/lon format (named 4D [#f1]_) but in an indexed format (named 3D [#f2]_).




Obtain 3D from CLM
------------------
CLM can directly output the needed 3D soil moisture fields. 

Example namelist for output::

  # daily SM, lat, lon grid:
  hist_fincl2 = 'SOILICE', 'SOILLIQ'
  # daily SM, column form
  hist_fincl3 = 'SOILLIQ','SOILICE'
  hist_nhtfrq = 0, -24, -24
  hist_mfilt = 1, 365, 365
  hist_dov2xy = .true., .true., .false.


Convert 4D to 3D
----------------





.. rubric:: Footnotes

.. [#f1] time x level x latitude x longitude
.. [#f2] time x level x index






