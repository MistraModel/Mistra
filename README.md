            _     _             
           (_)   | |            
  _ __ ___  _ ___| |_ _ __ __ _ 
 | '_ ` _ \| / __| __| '__/ _` |
 | | | | | | \__ \ |_| | | (_| |
 |_| |_| |_|_|___/\__|_|  \__,_|
                                
                                

!
! Copyright 1996-2017 the Authors
!
! Licensed under the EUPL, Version 1.1 only (the "Licence");
!
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at:
!   https://joinup.ec.europa.eu/software/page/eupl
!
! Unless required by applicable law or agreed to in writing,
! software distributed under the Licence is distributed on an
! "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
! either express or implied.
!
! See the Licence for the specific language governing permissions
! and limitations under the Licence.


Unless specifically specified in each file and/or routine, the principal Authors are:
  ** Andreas Bott ** (a.bott@uni-bonn.de)
  ** Roland von Glasow **

Additional work to make mistra code available, and significant technical changes in the code
have been made by:
  ** Josue Bock ** (josue.bock@univ-smb.fr)

---

__Mistra__ is 

For the license see

For Credits see


HOW TO RUN MISTRA -- in a nutshell:
-----------------------------------

1. Install __KPP-Mistra__. Please refer to the `README.md` file in the
[KPP-Mistra repository](https://github.com/MistraModel/KPP-Mistra) for
instructions.

1. Install a Fortran compiler (preferably __ifort__, alternatively
__gfortran__) and the __netcdf__ development libraries.


3) Check that C shell (csh) is installed on your system, then generate the mechanism files: in ./src/mech, run make

4) Check that netcdf libraries are installed on your system (locate netcdf.inc). If not, install the package libnetcdff-dev
Edit the Makefile (in ./src [or main?] directory), set the correct path to the netcdf libraries and include file. Depending on netcdf distribution, use "-lnetcdf" or "-lnetcdf -lnetcdff".

5) Compile the code using "make".

6) In the param/param_... file, set the appropriate path to Mistra input files. It is advised to use a directory that will not be duplicated for each version the user will run, but choose a generic directory instead.
Set the appropriate path to parent output directory, to the executable, and set the namelist to be used.
Run the param_... file to run the model.

For more detailed information see the manual in the ./doc directory
