__Mistra__ is a one dimensional chemical model of the Marine Boundary layer.

The [MistraModel](https://github.com/MistraModel) version in this repository is similar -- with very small, mostly cosmetic, differences -- to the __Mistra-UEA__ version. It is released under an open source license (GPLv2, see the `LICENSE` file), as agreed by the majority of the Mistra contributors.

For information about authorship, credits, acknowledgements, and funding, see the `CREDITS.md` file.

> [!Warning]
>
> Note that the Mistra-UEA version of the model was published (Bock et al., Geosci. Model Dev., 2022) without the knowledge or approval of most of the Mistra contributors, and under a licence which was not agreed upon by them. For more information, see the open discussion and the editor comments during the peer-review of the paper on the GMD website: https://gmd.copernicus.org/articles/15/5807/2022/gmd-15-5807-2022-discussion.html
>
> Therefore, __the use of the Mistra model in publications is *discouraged*__ until all the authorship and licensing issues are resolved. However, the code is open source and can be reused in other projects under the terms of the GPLv2 license.


How to install Mistra
---------------------

The following packages are required:

1. a *Fortran* compiler. The default is GNU `gfortran`. Use the command `gfortran --version` to check the compiler version. Alternatively the Intel `ifort` compiler can be used.

2. the *NetCDF Fortran* development library. On Debian (and Debian derivatives) the package is called `libnetcdff-dev`. Depending on the distribution, it may have a different name. Use the command `locate netcdf.inc` to check that the library is installed and its location.

3. the *C* shell (or the *tcsh* shell).

4. __optional__: the plotting software *Ferret* (https://ferret.pmel.noaa.gov/Ferret/) to use the scripts in the `scripts/` directory.


How to run Mistra
-----------------

1. Install __KPP-Mistra__. Refer to the `README.md` file in the [KPP-Mistra repository](https://github.com/MistraModel/KPP-Mistra) for instructions.

2. Start the *C* shell (or the *tcsh* shell).

3. Generate the mechanism files: in `./src/mech/`, execute `make`.

4. Edit `./src/Makefile` to set the path to the NetCDF library (`netcdf.inc`, see above) and the NetCDF compilation flags (`-lnetcdf` or -lnetcdff` or `-lnetcdf -lnetcdff`). Change the Fortran compiler if needed (the default is `gfortran`).

5. Compile the code: in `./src/`, execute `make`.

6. In the `param/param_...` file, set the appropriate path to the Mistra input files. It is advised to use a directory that will not be duplicated for each version the user will run, but choose a generic directory instead. Set the appropriate path to parent output directory, to the executable, and set the namelist to be used. Run the `param_...` file to run the model.

For more information, see the manual (`manual_Mistra_v9.0.pdf`) in the `doc/` directory.
