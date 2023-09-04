# DALESsubroutines
This repository contains all the subroutines that were changed for the master's thesis: 'Data assimilation of observed cloud fields in LES model' by Mats Steerneman, supervised by Stephan de Roode.
It also contains python scripts used to create DALES input files. Note that these are meant more as an example and will likely need to be changed for use by others.
In the thesis, a 3D-nudging method was added to DALES. Furthermore, the possibility of running DALES using only advection was added, as well as a way to start a regular LES run using pre-described fields.
In this document, it is described how to use the newly added inputs to 'namoptions' to use the different methods. Also, a description of all the changes in the subroutines is given.

Nudging.
1D-nudging already existed in DALES, 3D-nudging was added. 3D-nudging is currently only applied to thl and qt. 1D-nudging requires an input file called 'nudge.inp.iexpnr', where iexpnr is the experiment number in namoptions.
This file should contain the desired vertical mean profiles at the desired times. This is standard DALES so will not be elaborated on here.
3D-nudging also requires an input file: 'nudge3D.inp.iexpnr'. This should be a binary Fortran file (see nudgeinputfiles.py for an example of how to create these).
It should first contain an array which has all the points in time (in seconds) nudging should be done towards.
Then it should contain an array of the fluctuation in thl at all these points in time and all positions.
Then it should contain an array of the fluctuation in qt at all these points in time and all positions.
So:
t(tn)
thl_acc(tn,k,j,i)
qt_acc(tn,k,j,i)

The namoptions used in DALES nudging are (the full namoptions file I used for my reference run can also be found in this repository):
&NAMNUDGE
lnudge        = .false.      # Indicate whether nudging should be used in the experiment
knudgestart   = 1            # Lower height limit for nudging
knudgestop    = 79           # Upper height limit for nudging
t3Dnudgefac   = 10           # Nudging time scale for 3D nudging
t1Dnudgefac   = 10           # Nudging time scale for 1D nudging
tnudgestop    = 7200         # 
tnudgestart   = 3600
ntnudge3D     = 12
/

