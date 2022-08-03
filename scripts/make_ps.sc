#!/bin/bash
# script that writes a ferret script to plot the wanted reaction rates and executes it

# name of plot
# ------------
pl_name="rxn_f"

# determine number of meta print files

metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt $pl_name.20.plt $pl_name.21.plt $pl_name.22.plt $pl_name.23.plt $pl_name.24.plt $pl_name.25.plt $pl_name.26.plt $pl_name.27.plt $pl_name.28.plt $pl_name.29.plt $pl_name.30.plt $pl_name.31.plt $pl_name.32.plt $pl_name.33.plt $pl_name.34.plt $pl_name.35.plt $pl_name.36.plt $pl_name.37.plt $pl_name.38.plt $pl_name.39.plt $pl_name.40.plt $pl_name.41.plt $pl_name.42.plt $pl_name.43.plt $pl_name.44.plt $pl_name.45.plt $pl_name.46.plt $pl_name.47.plt $pl_name.48.plt $pl_name.49.plt $pl_name.50.plt $pl_name.51.plt $pl_name.52.plt $pl_name.53.plt $pl_name.54.plt $pl_name.55.plt $pl_name.56.plt $pl_name.57.plt $pl_name.58.plt $pl_name.59.plt $pl_name.60.plt $pl_name.61.plt $pl_name.62.plt $pl_name.63.plt $pl_name.64.plt $pl_name.65.plt $pl_name.66.plt $pl_name.67.plt $pl_name.68.plt $pl_name.69.plt $pl_name.70.plt $pl_name.71.plt $pl_name.72.plt $pl_name.73.plt $pl_name.74.plt $pl_name.75.plt $pl_name.76.plt $pl_name.77.plt $pl_name.78.plt $pl_name.79.plt $pl_name.80.plt $pl_name.81.plt $pl_name.82.plt $pl_name.83.plt $pl_name.84.plt $pl_name.85.plt $pl_name.86.plt $pl_name.87.plt $pl_name.88.plt $pl_name.89.plt $pl_name.90.plt $pl_name.91.plt $pl_name.92.plt $pl_name.93.plt $pl_name.94.plt $pl_name.95.plt $pl_name.96.plt $pl_name.97.plt $pl_name.98.plt $pl_name.99.plt $pl_name.100.plt"


gksm2ps -p portrait -l cps -d cps -o $pl_name.pre.ps $metafiles

# add page numbering to ps file and delete empty pages
ps2ps $pl_name.pre.ps $pl_name.ps



# clean up ----
# temporary files
#rm -rf pgtmp*
# prelim PS file
#rm -f $pl_name.pre.ps
# meta print files
#rm -f *.plt
