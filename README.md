            _     _             
           (_)   | |            
  _ __ ___  _ ___| |_ _ __ __ _   
 | '_ ` _ \| / __| __| '__/ _` |  
 | | | | | | \__ \ |_| | | (_| |  
 |_| |_| |_|_|___/\__|_|  \__,_|  
                                  
                                  


__Mistra__ is a one dimensional chemical model of the Marine Boundary layer.

For the license, see the `LICENSE` file.

For authorship, credits, acknowledgements, funding information, see the `CREDITS.md` file.


HOW TO RUN MISTRA -- in a nutshell:
-----------------------------------

1. Install __KPP-Mistra__. Please refer to the `README.md` file in the [KPP-Mistra repository](https://github.com/MistraModel/KPP-Mistra) for instructions.

2. Install a Fortran compiler (preferably __ifort__, alternatively __gfortran__) and the __netcdf__ development libraries.


3) Check that C shell (__csh__) is installed on your system, then generate the mechanism files: in `./src/mech/`, run `make`

4) Check that the __netcdf__ libraries are installed on your system (`locate netcdf.inc`). If not, install the package `libnetcdff-dev`

5) Edit the `Makefile` (in `./src/` [or main?] directory), set the correct path to the netcdf libraries and include file. Depending on netcdf distribution, use `-lnetcdf` or `-lnetcdf -lnetcdff`.

6) Compile the code using `make`.

7) In the `param/param_...` file, set the appropriate path to Mistra input files. It is advised to use a directory that will not be duplicated for each version the user will run, but choose a generic directory instead. Set the appropriate path to parent output directory, to the executable, and set the namelist to be used. Run the `param_...` file to run the model.

For more information, see the manual (`manual_Mistra_v9.0.pdf`) in the `doc/` directory.
