#!/bin/csh
# called by make_kpp.sc
# edit master_aqueous.eqn file so that identical mechanisms for the aqueous phase are produced
#
# created by Roland von Glasow
# updated in June 2016 by Josué Bock when upgrading KPP

#
# Copyright 1996-2017 the Authors
#
# Licensed under the EUPL, Version 1.1 only (the "Licence");
#
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#   https://joinup.ec.europa.eu/software/page/eupl
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the Licence is distributed on an
# "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
# either express or implied.
#
# See the Licence for the specific language governing permissions
# and limitations under the Licence.



# number of aqueous classes to be produced
# ========================================
set nm=("1" "2" "3" "4")

# define filenames
# ================
set in_file=("master_aqueous.eqn")
set out_file=("aq1.eqn" "aq2.eqn" "aq3.eqn" "aq4.eqn")
set head_file=("aer_eqn.head" "tot_eqn12.head" "tot_eqn34.head")
set title_file=("title.1" "title.2" "title.3" "title.4")

# edit mechanism
# ==============
foreach i ($nm)
    echo " duplicating aqueous mechanism for class $i"
    set f_out=$out_file[$i]
#    echo $f_out
# this special quoting ('"$i"') allows use of shell variables in sed
    sed 's/z/'"$i"'/g' $in_file > $f_out
end

# concatenate the different classes to aer.eqn and tot.eqn
# ========================================================

# aerosol
echo " creating equation file for 'aer' mechanism"
cat $head_file[1] $title_file[1] $out_file[1] $title_file[2] $out_file[2] > aer.eqn

# total
echo " creating equation file for 'tot' mechanism"
cat $head_file[2] $title_file[1] $out_file[1] $title_file[2] $out_file[2] $head_file[3] $title_file[3] $out_file[3] $title_file[4] $out_file[4] > tot.eqn

# clean up
# ========

rm -f aq[1-4].eqn
unset nm in_file out_file head_file title_file f_out
