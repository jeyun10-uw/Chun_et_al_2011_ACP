github-hosted UW Version of SAM, the System for Atmospheric Modeling

This repository was created by Peter Blossey to maintain a version of
SAM for use at the University of Washington with git as the version control
system.

#S12_CTLmod
Case setup files for CHUN et al 2022. 

prm: specify model parameter, initial/boundary conditions
	
File extention indicates setup for simulations.

-lc: Initial 12-hours simulation for spinup for Pristine6 case
-lcn: Control run (No plume) for Pristine6 case
-lcs: Plume run (with plume) for Pristine6 case

-lm: Initial 12-hours simulation for spinup for Middle6 case
-lmn: Control run (No plume) for Middle6 case
-lms: Plume run (with plume) for Middle6 case

-mph: Initial 12-hours simulation for spinup for Polluted6 case
-mphn: Control run (No plume) for Polluted6 case
-mphs: Plume run (with plume) for Polluted6 case

-mpm: Initial 12-hours simulation for spinup for Polluted3.5 case
-mpmn: Control run (No plume) for Polluted3.5 case
-mpms: Plume run (with plume) for Polluted3.5 case

-mpd: Initial 12-hours simulation for spinup for Polluted1.5 case
-mpdn: Control run (No plume) for Polluted1.5 case
-mpds: Plume run (with plume) for Polluted1.5 case

* IOP files for initial/boundary conditions

lc: S12_CTL_Dec2010c_Turnedv2_QuasiEquil_MixedLayerInit_q6w1.5.nc
lm: S12_CTL_Dec2010c_Turnedv2_QuasiEquil_MixedLayerInit_q6w1.3.nc
mph: S12_CTL_Dec2010c_Turnedv2_QuasiEquil_MixedLayerInit_q6w1.nc
mpm: S12_CTL_Dec2010c_Turnedv2_QuasiEquil_MixedLayerInit_q3.5w1.nc
mpd: S12_CTL_Dec2010c_Turnedv2_QuasiEquil_MixedLayerInit_q1.5w1.nc

* An example to submit a job to machine

resub.lc

* For more Information, please contanct jeyun.chun@gmail.com
