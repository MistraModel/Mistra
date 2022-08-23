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

To install Mistra the following packages are required:

1. a Fortran compiler. The default is GNU `gfortran`. Use the command `gfortran --version` to check the compiler version. Alternatively the Intel `ifort` compiler can be used.

2. the NetCDF Fortran development library. On Debian (and Debian derivatives) the package is called `libnetcdff-dev`. Depending on the distribution, it may have a different name. Use the command `locate netcdf.inc` to check that the library is installed.

3. the C shell

To set up and run the Mistra model, follow these steps:

1. Install __KPP-Mistra__. Please refer to the `README.md` file in the [KPP-Mistra repository](https://github.com/MistraModel/KPP-Mistra) for instructions.

2. Start the C shell (__csh__), then generate the mechanism files: in `./src/mech/`, run `make`.

3. Edit the `Makefile` (in `./src/` directory), set the correct path to the netcdf libraries and include file. Depending on netcdf distribution, use the flag `-lnetcdf` or -lnetcdff` or `-lnetcdf -lnetcdff`. Select the compiler if needed.

4. Compile the code using `make`.

5. In the `param/param_...` file, set the appropriate path to Mistra input files. It is advised to use a directory that will not be duplicated for each version the user will run, but choose a generic directory instead. Set the appropriate path to parent output directory, to the executable, and set the namelist to be used. Run the `param_...` file to run the model.

For more information, see the manual (`manual_Mistra_v9.0.pdf`) in the `doc/` directory.
