prescribeSM CESM 1.2.x
======================

*prescribeSM_CESM_1.2.x* allows to prescribe SOILLIQ and SOILICE in CLM 4.0 for CESM 1.2.x.

Features
--------

- Prescribe SOILLIQ and SOILICE (-> SM)
- Prescribe SM globally or regionally
- Prescribe daily mean or monthly mean SM
- Two ways to do this:

  - Prescribe SOILLIQ *and* SOILICE (standard)
  - Prescribe SOILLIQ only if no SOILICE is present
  
- Define used Inputfile via namelist

Documentation
-------------

You can find documentation on http://prescribesm-cesm-12x.readthedocs.org/

Prerequisites
-------------
- CLM 4.0 history file with the variables SOILLIQ and SOILICE

  - Annual file with monthly averages (12 time steps on file)
  - Annual file with daily averages (365 time steps on file)

- You need the output in column form and not as a lat lon grid!

Example namelist for output::

  # daily SM, lat, lon grid:
  hist_fincl2 = 'SOILICE', 'SOILLIQ'
  # daily SM, column form
  hist_fincl3 = 'SOILLIQ','SOILICE'
  hist_nhtfrq = 0, -24, -24
  hist_mfilt = 1, 365, 365
  hist_dov2xy = .true., .true., .false.


- A small python script that shows how to get from lat, lon form to col form is given (clm_col_to_xy_example.py)

Installation
------------
Get the code with:
``git clone https://github.com/mathause/prescribeSM_cesm_1.2.x prescribeSM``
Then add the 3 source files in the SourceMod/src.clm folder and build cesm normally.

Usage
-----
- if you set SOILLIQ or SOILICE to -1 at a gridpoint, LIQ and ICE are calculated interactively at this point!
- this can be used to prescribe SM only regionally
- to prescribe SM only at certain time steps may not work with the current implementation

.. NOTE::
   CLM does not provide a sensible result if you prescribe only *one* of SOILLIQ/ SOILICE.
  
- add the namelist called *&prescribe_SM* to a file called *prescribe_SM_nl* to the run folder of the simulation (RUNDIR in env_run.xml)
- see also the example namelist provided

Namelist
^^^^^^^^
  
pSMfile : string
  path to the history file containing SOILLIQ and SOILICE
monthly : bool
  If .true. assumes monthly mean SM is provided. If false uses daily input.
pSMtype : int, {1, 2, 3}
 Defines how to prescribe SM. 1: default: prescribe SOILLIQ and SOILICE. 2: prescribe SM only if there is no SOILICE
 present (never prescribes SOILICE). Lets the model decide if there is ICE or not at a certain gridpoint.
 3: yields no physically sensible output.
 
.. WARNING::
   if you set monthly=.true. but have a daily SM input file it will still work (uses the first 12 days as the months)


Source Code
-----------

- Source Code: github.com/mathause/prescribeSM_cesm_1.2.x

License
-------

?? Probably the same as CESM ??


Info
----

Author: Mathias Hauser
Date:   June 2014

Model Version:
CESM 1.2.1
CLM 4.0
