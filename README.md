# DALESsubroutines
This repository contains all the subroutines that were changed for the master's thesis: 'Data assimilation of observed cloud fields in LES model' by Mats Steerneman, supervised by Stephan de Roode.  
It also contains Python scripts used to create DALES input files. Note that these are meant more as an example and will likely need to be changed for use by others.  
In the thesis, a 3D-nudging method was added to DALES. Furthermore, the possibility of running DALES using only advection was added, as well as a way to start a regular LES run using pre-described fields.  
In this document, it is described how to use the newly added inputs to 'namoptions' to use the different methods. The changes made to the subroutines are indicated by comments in these subroutines, which are always preluded by 'MS'.

## Nudging 
*adapted subroutines: 'modnudge.f90'*

1D-nudging already existed in DALES, 3D-nudging was added. 3D-nudging is currently only applied to thl and qt. 1D-nudging requires an input file called 'nudge.inp.iexpnr', where iexpnr is the experiment number in namoptions.  
This file should contain the desired vertical mean profiles at the desired times. This is standard DALES so will not be elaborated on here.  
3D-nudging also requires an input file: 'nudge3D.inp.iexpnr'. This should be a binary Fortran file (see nudgeinputfiles.py for an example of how to create these).  
It should first contain an array that has all the points in time (in seconds) nudging should be done towards.  
Then it should contain an array of the fluctuation in thl at all these points in time and all positions.  
Then it should contain an array of the fluctuation in qt at all these points in time and all positions.  
So:  
t(tn)  
thl_acc(tn,k,j,i)  
qt_acc(tn,k,j,i)  

The 'namoptions' used in DALES nudging are (the full 'namoptions' file I used for my reference run can also be found in this repository):

``
&NAMNUDGE  
lnudge        = .false.      # Indicate whether nudging should be used in the experiment (true = use, false = do not use)  
knudgestart   = 1            # Lower height limit for nudging (index)  
knudgestop    = 79           # Upper height limit for nudging (index)  
t3Dnudgefac   = 10           # Nudging time scale for 3D nudging (s)  
t1Dnudgefac   = 10           # Nudging time scale for 1D nudging (s)  
tnudgestop    = 7200         # Time at which to deactivate nudging (s)  
tnudgestart   = 0            # Time at which to activate nudging (s)  
ntnudge3D     = 12           # Amount of time points for which desired nudge fields are given (tn in the above)  
/
``

The provided subroutine 'modnudge.f90' gives the nudging subroutine for general use. In 'modnudge_large.f90', the input file 'nudge.inp.iexpnr' became to large and was hardcoded to use two input files. This could be improved upon for later versions.

## Persistence or advection only method
*adapted subroutines: 'modstartup.f90', 'modglobal.f90', 'tstep.f90'*

Using the persistence method in DALES is straightforward. It needs an input file called 'advfield.inp.iexpnr', containing the fields of 'thl' and 'qt' to be advected. An example python script to create this file is given in 'input_Pers_or_RefCold.py'.
So, the file contains:

thl(k,j,i)

qt(k,j,i)

Additionally, one should give the desired mean horizontal wind speeds at each height in 'prof.inp.iexpnr', as in standard DALES.

The only input in 'namoptions' is given in the &PHYSICS block:

``
&PHYSICS
ladvectonly=  .false.        # Set to true to use the persistence method, to false not to use it.
/
``

## Standard LES with reference fields
*adapted subroutines: 'modstartup.f90', 'modglobal.f90'*

It is also possible to use regular LES runs, that start with the correct reference fields. For this, one must supply the initial fields of 'thl' and 'qt' in a file named: 'coldinifield.inp.iexpnr'. It is created in the same way as described for 'advfield.inp.iexpnr' above. Its only input in namoptions is given in the &RUN block:

``
&RUN
lcoldstartfiles = .false.    # Set to true to use a regular LES run with thermodynamic reference field initialization, and to false not to.
/
``

## Printing LWP and SWD(down) at the surface
*adapted subroutines: 'modtimestat.f90'*

In this thesis, the functionality was added to DALES to print LWP and SWD(down) values at the surface for various points in the simulation. This was done to check these values against reference values, to find the accuracy of an experiment's forecast. In the end, it was not used. However, for completeness, a short description is given here. The added commands in 'namoptions':
``
&NAMTIMESTAT
llwp_pos    = .true.         # Indicate whether the LWP and SWD should be printed at specified locations.
ips1        = 32             # i index for position 1
jps1        = 32             # j index for position 1
ips2        = 12             # etc.
jps2        = 46
ips3        = 53
jps3        = 22
ips4        = 5
jps4        = 18
/
``
The LWP and SWD(down) are then plotted at time intervals specified by 'dtav' (in &NAMTIMESTAT) in output files 'tmradpsn.iexpnr' (where n = 1,2,3,4 indicates the position).

This concludes the description of all 'namoptions' added to DALES in my thesis. The subroutines used are found in this repository, as well as example Python scripts for creating the input files.
