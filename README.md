            _     _             
           (_)   | |            
  _ __ ___  _ ___| |_ _ __ __ _ 
 | '_ ` _ \| / __| __| '__/ _` |
 | | | | | | \__ \ |_| | | (_| |
 |_| |_| |_|_|___/\__|_|  \__,_|
                                  

**Mistra** is a one dimensional chemical model of the Marine Boundary layer.

For the license, see the `LICENSE` file.

For authorship, credits, acknowledgements, funding information, see the `CREDITS.md` file.


How to install Mistra
---------------------

The following packages are required:

1. a *Fortran* compiler. The default is GNU `gfortran`. Use the command `gfortran --version` to check the compiler version. Alternatively the Intel `ifort` compiler can be used.

2. the *NetCDF Fortran* development library. On Debian (and Debian derivatives) the package is called `libnetcdff-dev`. Depending on the distribution, it may have a different name. Use the command `locate netcdf.inc` to check that the library is installed and its location.

3. the *C* shell (or the *tcsh* shell).


How to run Mistra
-----------------

1. Install **KPP-Mistra**. Refer to the `README.md` file in the [KPP-Mistra repository](https://github.com/MistraModel/KPP-Mistra) for instructions.

2. Start the **C** or the **tcsh** shell.

3. Generate the mechanism files: in `./src/mech/`, execute `make`.

4. Edit `./src/Makefile` to set the path to the NetCDF library (`netcdf.inc`, see above) and the NetCDF compilation flags. Change the Fortran compiler if needed (the default is `gfortran`).

5. Compile the code: in `./src/`, execute `make`.

6. In the `param/param_...` file, set the appropriate path to the Mistra input files. It is advised to use a directory that will not be duplicated for each version the user will run, but choose a generic directory instead. Set the appropriate path to parent output directory, to the executable, and set the namelist to be used. Run the `param_...` file to run the model.

For more information, see the manual (`manual_Mistra_v9.0.pdf`) in the `doc/` directory.
