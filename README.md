### prescribe SM in CLM/ CESM

Author: Mathias Hauser
Date:   June 2014

Model Version:
CESM 1.2.1
CLM 4.0

introduce prescribed soil moisture

- you can prescribe monthly mean or daily mean SM (linearly interpolated in between)
- you need the variables SOILLIQ and SOILICE from CLM (history file)
- you need the output in column form and not as a lat lon grid!

example namelist for output:

```
hist_fincl2 = 'SOILICE', 'SOILLIQ' -> daily SM, lat, lon grid
hist_fincl3 = 'SOILLIQ','SOILICE'  -> daily SM, column form
hist_nhtfrq = 0, -24, -24
hist_mfilt = 1, 365, 365
hist_dov2xy = .true., .true., .false.
```

- there is a small python script that shows how to get from lat, lon form to col form
- add the prescribe_SM_nl to the run directory of the model and adjust it accordingly
- Warning: if you set monthly=.true. but have a daily SM input file it will still work (uses the first 12 days as the months)
- if you set SOILLIQ or SOILICE to -1 at a gridpoint, LIQ and ICE are calculated interactively at this point!



