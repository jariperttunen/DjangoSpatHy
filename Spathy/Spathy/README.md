# Spathy
If you read this, you are working with version for python3

Spatial Hydrology of boreal catchments

This repository contains working copies of 'new' version being develoed by Samuli. Uploaded Jan 5th, 2018 files that 'work' but are not final.

Ongoing revisions to previous Spathy-versions + ToDo -list:

CanopyGrid:

* modified evaporation calculations. Both Efloor and Evap from canopy storage are now computed using P-M equation and resistance terms include dependency to wind speed. This should alter sensitivity to climatic forcing.
* Snow sublimation from canopy now according to Jules -model and Essery et al., Hedström & Pomeroy. Seems to increase variability of SWEmax across landscape (more interception of snow)
* was able to show that gsref is not free parameter but constrained by leaf Amax and WUE.

TODO: 
* Efloor restrictions; as in APES moss module?
* clean initialization routine
* Add topographical constraints to radiation

BucketGrid:

* now 2-layer bucket where toplayer acts as interception storage resembling organic layer. Interception computed as in CanopyGrid. Efloor is taken from this layer and Transpi from the lower layer.
* added top layer properties into .ini -file
* ground water returnflow now handled in 'watbal'; in case of inflow excess toplayer and pond storages are updated and remaining is routed to surface runoff

TODO:
* clean initialization routine, change state variable names etc.
* add bypass fraction through toplayer?

Topmodel:

* re-classify tail of TWI-distribution to avoind arbitrary large TWI's. Now hard-coded that highest 2.5% of distribution is set to 97.5th percentile. Give parameter as input insted.
* changed outputs of Topmodel_homogenous so that it matches Spathy_4 structure
* In old (and still in current) versions: TWI = ln[a/tan(b)] where a is flow accumulation area. To be consistent with TWI definition, this needs to be TWI = ln[(a/dx) / tan(b)], where a/dx is flow accumulation area per unit grid length. Implementation of this needs re-calibration of parameter m.

TODO:
* finish above modifications & clean code

Spathy4:

* changed the way BucketGrid and Topmodel storages are linked. The current version seems more realistic and handles returnflow better. Note that now surface runoff is produced only in BucketGrid while Topmodel gives total subsurface flow to stream network, and grid of returnflow (that is then used to update BucketGrid storages)
* Outputs are slightly changed - the changes are not yet in NetCDF -structure

TODO:
* clean inputs & initialization
* clean outputs & code
* re-calibrate and test
