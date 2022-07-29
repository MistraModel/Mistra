#!/bin/csh
# change kpp output of different mechanisms
# so that they can be used in one programme
#
# created by Roland von Glasow
# updated in June 2016 by Josue Bock when upgrading KPP

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

# ============================================================================================
# Modifications :
# -------------
  # 31-Oct-2017  Josue Bock  - replaced \+ by \{1,\}
  #                          The former is GNU extension, and is not recognised by all systems
  #                          - changed rm, mv, cp into \rm, \mv, \cp to use native functions
  #                          without aliases. If aliases were defined (for instance, cp='cp -i'),
  #                          they were not overwritten by the local flag, and the user was asked
  #                          to confirm all steps.
  #
  # 29-May-2021  Josue Bock  After converting kpp.f into kpp.f90, need to make header files
  #                          (*_Global.h and *_Parameters.h) compliant with both f77 and f90, by
  #                          converting comment signs (c or C) into !. Apply sed substitution
  # 24-Jun-2021  Josue Bock  Move this to mult_ch.sc

# == End of modifications =====================================================================


# define filenames
# ================
# prefix = mechanisms names
set prefix=$1

# file_k : input files for KPP
set file_k=$1.k

# suffix = KPP output file names, necessary for MISTRA
set suffix=("_Main.f" "_Rates.f" "_Integrator.f" "_Function.f" "_Jacobian.f" "_LinearAlgebra.f" "_JacobianSP.f" "_Monitor.f")

# file_f : necessary output files from KPP concatenated in a single file (per mechanism) for MISTRA
set file_f=$1.f

# file_r : output file with reaction rates, used for budget files
set file_r=$1_Function.f

# file_p : another header from which is read the total number of reactions
#          which is needed for budget file (do loop at the end, see tail.bud)
set file_p=$1_Parameters.h

# define appendices
# =================
set appendix=$2


# ===============
# make mechanisms
# ===============
echo "creating mechanisms from master files"
./make_aq_mech.sc

# ===========
# execute kpp
# ===========
echo "executing kpp on mechanism $prefix"
kpp $file_k

# ========================================================
# concatenate kpp files into one single file per mechanism
# ========================================================
# first need to remove existing file to use >> thereafter
rm -f $file_f
foreach j ($suffix)
    cat $prefix$j >> $file_f
end

# =========================================
# change concatenated files with mult_ch.sc
# =========================================
echo "tuning KPP files for mechanism $prefix"
    ./mult_ch.sc $prefix $appendix

# =======================
# make budget subroutines
# =======================
echo "make budget subroutine for mechanism $prefix"
# step 1: copy rates expression from file_r
sed -n '/ \{6\}A([0-9]\{1,\}) = RCT([0-9]\{1,\}/p' $file_r > ! tmp0
# step 2: substitute A(...) (KPP variable) by bg(1,kl,...) (MISTRA variable)
sed 's/A(\([0-9]\{1,\}\))/bg(1,\1,kl)/g' tmp0 > ! tmp1
# step 3: substitute RCT, V and F (local KPP variable) by RCONST, VAR and FIX (global KPP variable)
sed 's/RCT(/RCONST(/g' tmp1 > ! tmp2
sed 's/V(/VAR(/g' tmp2 > ! tmp3
sed 's/F(/FIX(/g' tmp3 > ! tmp4
# step 4: concatenate start and end of budget file
cat head.bud tmp4 tail.bud > ! tmp5
# step 5: rename subroutines and include files
sed 's/KPP_ROOT/'"$prefix"'/' tmp5 > ! tmp6
sed 's/appendix/'"$appendix"'/' tmp6 > ! tmp7
# step 6: get number of reactions from KPP _Parameters.h file
set nreact=`sed -n '/NREACT =/p' $file_p | awk -F= '{printf "%d",$2}'`
#        and use it in the integrated budget calculation (see tail.bud)
sed 's/do i=1,xxxx/do i=1,'"$nreact"'/' tmp7 > ! tmp8
\mv -f tmp8 bud_$appendix.f
\rm -f tmp[0-7]


# ==========================================
# Copy all necessary files to main directory
# ==========================================
echo "copy chemistry files with option "-f" (y) ?"
set flag_yes=$<
set flag='-i'
if ($flag_yes == 'y') set flag='-f'

\cp $flag  $1.f $1*.h $1_*.dat bud_$2.f ..


# ========
# clean up
# ========
# aer and tot eqn files can now be deleted (remember: they are created by make_aq_mech.sc)
\rm -f aer.eqn tot.eqn
# output files from KPP are also deleted
\rm -f $1_*.f Makefile_*


unset nm
unset prefix suffix appendix
unset file_k file_f file_r file_p
unset nreact
unset flag_yes flag
