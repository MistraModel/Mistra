#!/bin/csh
# called by make_kpp.sc
# change *.f *.h for each input
#
# created by Roland von Glasow
# updated in June 2016 by Josue Bock when upgrading KPP

# ======================================================================= #
# Copyright 1996-2017 the Authors
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ======================================================================= #


# Modifications :
# -------------
  # 31-Oct-2017  Josue Bock  - replaced \? and \+ by \{0,1\} and \{1,\} respectively
  #                          The former are GNU extensions, and are not recognised by all systems
  #                          - corrected the first line #!/bin/csh instead of ##/bin/csh
  #                          - replaced rm, mv by \rm, \mv to use native functions without potential
  #                          user defined aliases (for instance, rm='rm -i')
  #                          - reorganised the tmp files so that they are all saved until the end
  #                          in order to facilitate debug, if needed.
  #
  # 14-May-2021  Josue Bock  add "exit from Rosenbrock" in the string list so that in case of error,
  #                          the user can know which mechanism was run
  #
  # 24-Jun-2021  Josue Bock  move the comment sign replacement (c, C => !) here instead of make_kpp.sc
  #                          and remove sed -i option ("in-place editing", to edit the input file directly)
  #                          since it was proven to be less portable (issues with MacOS)

# == End of modifications =======================================================================


# define variables
# ----------------
set file_f=$1.f
set file_h=("$1_Global.h" "$1_Parameters.h")
set file_s=$1_Sparse.h

set appendix=$2

# ==============
# work on file_f
# ==============

# delete unnecessary subroutines
# ------------------------------
set del_SR=("Update_SUN" "Update_PHOTO" "KppDecompCmplx" "KppSolveIndirect" "KppSolveCmplx" "KppSolveTR" "WCOPY" "WSCAL" "WLAMCH" "SET2ZERO" "WDOT" "WADD" "Jac_SP_Vec" "JacTR_SP_Vec")
set fold=$file_f
set inn=0
foreach i ($del_SR)
  @ inn=$inn + 1
  set fnew=tmp_f_$inn
# [CEN] accomodates either a comment, the end of word "SUBROUTINE" or the end of word "FUNCTION"
# the special quoting ('"$i"') allows use of shell variables in sed
  sed '/[CEN] '"$i"'/,/End of '"$i"'/d' $fold > ! $fnew
  set fold=$fnew
end

# delete calls to Update_SUN and Update_PHOTO
# -------------------------------------------
set del_CALL=("Update_SUN" "Update_PHOTO")
foreach i ($del_CALL)
  @ inn=$inn + 1
  set fnew=tmp_f_$inn
  sed '/CALL '"$i"'/,/)/d' $fold > ! $fnew
  set fold=$fnew
end

# rename remaining subroutines
# ----------------------------
set ch_SR=("Update_RCONST" "INTEGRATE" "RosenbrockIntegrator" "ros_ErrorNorm" "ros_FunTimeDerivative" "ros_PrepareMatrix" "ros_ErrorMsg" "Ros2" "Ros3" "Ros4" "Rodas3" "Rodas4" "DecompTemplate" "SolveTemplate" "FunTemplate" "JacTemplate" "KppDecomp" "WAXPY" "KppSolve" "WLAMCH_ADD" "Jac_SP" "exit from Rosenbrock")
# :q specify the string format/object for csh. It allows to have spaces in the searched strings
foreach i ($ch_SR:q)
  @ inn=$inn + 1
  set fnew=tmp_f_$inn
  sed 's/'"$i"'/&_'"$appendix"'/g' $fold > ! $fnew
  set fold=$fnew
end

# Special cases:
# 1) for SR Fun:
#   needs to avoid replace every occurrence of Fun*, which is a substring of "FunTemplate" and is obviously also a substring of "Function"
#   the search pattern accomodates "Fun(", "Fun (", "Fun -" and "Fun f"
@ inn=$inn + 1
set fnew=tmp_f_$inn
# " \?" matches zero or one space (backslash is needed in pattern context)
# it is a GNU extension, for portability reasons it has been replaced by \{0,1\}
#
# [f(-] matches one of the 3 character inside the brackets
# \( ... \) saves the pattern in between, which is recalled with \1
sed 's/Fun\( \{0,1\}[f(-]\)/Fun_'"$appendix"'\1/g' $fold > ! $fnew
set fold=$fnew

# 2) for SR Rosenbrock : add ( to the searched pattern
@ inn=$inn + 1
set fnew=tmp_f_$inn
sed 's/Rosenbrock(/Rosenbrock_'"$appendix"'(/g' $fold > ! $fnew
set fold=$fnew

# delete unnecessary LU_IROW data
# -------------------------------
# explanations, if future change is needed:
#   \+ means one or more. It is a GNU extension, for portability reasons it is now replaced by \{1,\}
#   \{1\} means exactly one
# searched ending string:
#   - starts by "*"
#   - optional pattern: optional space(s) followed by at least 1 number followed by comma
#   - mandatory pattern: optional space(s) followed by at least 1 number followed by 1 space and /
@ inn=$inn + 1
set fnew=tmp_f_$inn
sed '/DATA( LU_IROW(i)/,/\*\( *[0-9]\{1,\},\)*\( *[0-9]\{1,\} \/\)\{1\}/d' $fold > ! $fnew
set fold=$fnew

# rename LU_ variables, BLOCK DATA JACOBIAN_SPARSE_DATA and MONITOR_DATA
# ----------------------------------------------------------------------
set ch_VAR=("LU_ICOL" "LU_CROW" "LU_DIAG" "JACOBIAN_SPARSE_DATA" "MONITOR_DATA")
foreach i ($ch_VAR)
  @ inn=$inn + 1
  set fnew=tmp_f_$inn
  sed 's/'"$i"'/&_'"$appendix"'/g' $fold > ! $fnew
  set fold=$fnew
end
\mv -f $fnew $file_f



# ==============
# work on file_h
# ==============

# delete unnecessary variables in _Global.h
#  (these ones are only declared)
set del_VAR=("SUN" "TEMP" "RTOLS" "TSTART" "TEND" "CFACTOR" "LOOKAT" "MONITOR" "SMASS" "EQN_TAGS")
set fold=$file_h[1]
set inn=0
foreach i ($del_VAR)
  @ inn=$inn + 1
  set fnew=tmp_h1_$inn
  sed '/C '"$i"' -/,/\/ '"$i"'/d' $fold > ! $fnew
  set fold=$fnew
end

# change names of common blocks *GDATA in _Global.h
set ch_VAR=("\/GDATA" "\/CHARGDATA" "\/INTGDATA")
foreach i ($ch_VAR)
  @ inn=$inn + 1
  set fnew=tmp_h1_$inn
  sed 's/'"$i"'/&_'"$appendix"'/g' $fold > ! $fnew
  set fold=$fnew
end

# change comment signs in _Global.h to make them compliant with both F77 and F90
# NB: the escape \! in sed is necessary for csh script, in bash this would simply be !
set fold=$fnew
@ inn=$inn + 1
set fnew=tmp_h1_$inn
sed 's/^[cC]/\!/g' $fold > ! $fnew
\mv -f $fnew $file_h[1]

# delete unnecessary variables in _Parameters.h
#  (these one are defined as parameters, the searched pattern is different)
set del_VAR2=("NONZERO" "CNVAR" "NLOOKAT" "NMONITOR" "NMASS")
set fold=$file_h[2]
set inn=0
foreach i ($del_VAR2)
  @ inn=$inn + 1
  set fnew=tmp_h2_$inn
  sed '/C '"$i"' -/,/PARAMETER ( '"$i"' =/d' $fold > ! $fnew
  set fold=$fnew
end

# change comment signs in _Parameters.h to make them compliant with both F77 and F90
# NB: the escape \! in sed is necessary for csh script, in bash this would simply be !
set fold=$fnew
@ inn=$inn + 1
set fnew=tmp_h1_$inn
sed 's/^[cC]/\!/g' $fold > ! $fnew
\mv -f $fnew $file_h[2]



# ==============
# work on file_s
# ==============
# delete unnecessary variables
set del_VAR3=("LU_IROW")
set fold=$file_s
set inn=0
foreach i ($del_VAR3)
  @ inn=$inn + 1
  set fnew=tmp_s_$inn
  sed '/C '"$i"' -/,/\/ '"$i"'/d' $fold > ! $fnew
  set fold=$fnew
end

# change names of common block SDATA
@ inn=$inn + 1
set fnew=tmp_s_$inn
sed 's/SDATA/&_'"$appendix"'/g' $fold > ! $fnew
set fold=$fnew

# change names of LU_ variables
set ch_VAR2=("LU_ICOL" "LU_CROW" "LU_DIAG")
foreach i ($ch_VAR2)
  @ inn=$inn + 1
  set fnew=tmp_s_$inn
  sed 's/'"$i"'/&_'"$appendix"'/g' $fold > ! $fnew
  set fold=$fnew
end
\mv -f $fnew $file_s



# ========
# clean up
# ========
\rm -f tmp_[fhs]*
unset file_f file_h file_s appendix
unset fold fnew inn
unset del_SR ch_SR ch_VAR del_VAR del_VAR2 del_VAR3 ch_VAR2
