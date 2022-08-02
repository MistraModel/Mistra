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

# Intel Fotran Compiler (tested up to version 2021.3.0)
FC=ifort
# -- Preferred flags for production run
CFLAGS= -O2 -fpe0 -traceback
# -- Alternative flags for production run
#CFLAGS= -O2 -fpe0
#CFLAGS= -O2 -fp-model source -traceback
#CFLAGS= -O2 -r8 -fpe0
#CFLAGS= -O2 -r8 -fpe0 -traceback
#CFLAGS= -O2 -r8 -fpe0 -traceback -noalign
#CFLAGS= -O2 -r8 -fpe0 -traceback -fltconsistency
#CFLAGS= -O2 -r8 -fpe0 -traceback -fp-model source
#CFLAGS= -O2 -noalign -fp-model source -traceback # some (not all) GEOS-Chem options for standard runs
#CFLAGS= -O2 #-check all -warn declarations -warn uncalled -warn unused#-g #-warn all
#CFLAGS= -O2 -ipo #-check all -warn all
#CFLAGS= -fast # = -O3 -ipo -static (-xP)
# -- Debug flags for development stages
#CFLAGS= -O0 -r8 -fpe0 -traceback -fltconsistency -check uninit -warn all -check all
#CFLAGS= -O0 -r8 -init=snan,arrays -traceback
#CFLAGS= -O0 -debug all -g -check bounds -check arg_temp_created -fpe0 -ftrapuv # GEOS-Chem debug options
#CFLAGS= -O0 -r8 -fpe0 -traceback -warn all -check all
#CFLAGS= -check all -traceback -warn all

# GNU Fortran (tested up to version 7.5.0)
#FC=gfortran
# --  Preferred run flags
#CFLAGS= -O3 -w -fbacktrace # -fautomatic -fno-align-commons
# --  Alternative run flags
#CFLAGS= -O3 -fdefault-real-8 -fdefault-double-8 -fbacktrace
# -- debug flags
#CFLAGS= -O0 -g -fdefault-real-8 -fdefault-double-8 -fbacktrace -fcheck=bounds -ffpe-trap=invalid,zero,overflow -finit-real=snan -Wall -fcheck=all
#CFLAGS= -g -fbacktrace -ffpe-trap=zero,underflow,overflow,invalid
#CFLAGS= -fbacktrace -fdump-fortran-optimized -fdump-fortran-original -fdump-parse-tree -ffpe-trap= -ffpe-summary= -g -Wsurprising -Walign-commons -Wall


# Makefile memo
# -------------
# $@ = target name
# $< = first dependency
# $^ = all dependencies


# Inclusion of Makefile in dependencies lists:
# --------------------------------------------
# when compilations option are changed in the Makefile, all files must be recompiled.
# To avoid mistakes in this situation, the Makefile has been added in the dependency list of all modules,
# and in the dependency list of all objects which don't already include a module

mistra: precision.o constants.o global_params.o common_modules.o\
	data_surface.o\
	file_unit.o\
	str.o radinit.o  nrad.o outp.o nuc.o jrate.o activity.o \
	tot.o bud_t.o bud_s_t.o \
	aer.o bud_a.o bud_s_a.o \
	gas.o bud_g.o bud_s_g.o \
	kpp.o utils.o \
	config.o \
	mod_out_netCDF.o out_netCDF.o
	$(FC) $(CFLAGS) -o $@ $^ -lnetcdf -lnetcdff

# modules
global_params.o global_params.mod: global_params.f90 tot_Parameters.h Makefile
	rm -f global_params.mod
	$(FC) $(CFLAGS) -c $<

common_modules.o gas_common.mod: common_modules.f90 global_params.mod Makefile
	rm -f gas_common.mod
	$(FC) $(CFLAGS) -c $<

constants.o constants.mod: constants.f90 precision.mod Makefile
	rm -f constants.mod
	$(FC) $(CFLAGS) -c $<

data_surface.o data_surface.mod: data_surface.f90 precision.mod Makefile
	rm -f data_surface.mod
	$(FC) $(CFLAGS) -c $<

file_unit.o file_unit.mod: file_unit.f90 Makefile
	rm -f file_unit.mod
	$(FC) $(CFLAGS) -c $<

mod_out_netCDF.o cdf_var_gas.mod: mod_out_netCDF.f90 gas_common.mod Makefile
	rm -f cdf_var_gas.mod
	$(FC) $(CFLAGS) -c $<

precision.o precision.mod: precision.f90 Makefile
	rm -f precision.mod
	$(FC) $(CFLAGS) -c $<

# files containing modules
tot.o kpp_tot_parameters.mod kpp_tot_global.mod: tot.f tot_Parameters.h tot_Global.h tot_Sparse.h aer_mk.dat aer_km.dat \
	constants.mod gas_common.mod global_params.mod Makefile
	rm -f kpp_tot_parameters.mod kpp_tot_global.mod
	$(FC) $(CFLAGS) -c $<

aer.o kpp_aer_parameters.mod kpp_aer_global.mod: aer.f aer_Parameters.h aer_Global.h aer_Sparse.h aer_mk.dat aer_km.dat \
	constants.mod gas_common.mod global_params.mod Makefile
	rm -f kpp_aer_parameters.mod kpp_aer_global.mod
	$(FC) $(CFLAGS) -c $<

gas.o kpp_gas_parameters.mod kpp_gas_global.mod: gas.f gas_Parameters.h gas_Global.h gas_Sparse.h \
	constants.mod gas_common.mod global_params.mod Makefile
	rm -f kpp_gas_parameters.mod kpp_gas_global.mod
	$(FC) $(CFLAGS) -c $<

config.o config.mod: config.f90 file_unit.mod precision.mod Makefile
	rm -f config.mod
	$(FC) $(CFLAGS) -c $<

# other files
utils.o: utils.f90 config.mod file_unit.mod gas_common.mod global_params.mod \
	kpp_aer_global.mod \
	kpp_gas_global.mod \
	kpp_tot_global.mod \
	precision.mod \
	Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

str.o: str.f90 config.mod constants.mod file_unit.mod gas_common.mod global_params.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

radinit.o: radinit.f90 config.mod constants.mod file_unit.mod global_params.mod precision.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

nrad.o: nrad.f90 config.mod constants.mod file_unit.mod global_params.mod precision.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

outp.o: outp.f90 config.mod constants.mod file_unit.mod global_params.mod gas_common.mod precision.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

out_netCDF.o: out_netCDF.f cdf_var_gas.mod config.mod file_unit.mod global_params.mod gas_common.mod precision.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $< -I/usr/include

jrate.o: jrate.f config.mod file_unit.mod global_params.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

nuc.o: nuc.f90 config.mod constants.mod file_unit.mod global_params.mod gas_common.mod precision.mod \
	Makefile
	rm -f mod_nuc.mod
	$(FC) $(CFLAGS) -o $(@F) -c $<

kpp.o: kpp.f90 tot_Parameters.h tot_Global.h aer_Parameters.h aer_Global.h gas_Parameters.h gas_Global.h \
	kpp_tot_global.mod kpp_tot_parameters.mod \
	kpp_aer_global.mod kpp_aer_parameters.mod \
	config.mod constants.mod file_unit.mod gas_common.mod global_params.mod precision.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

activity.o: activity.f90 constants.mod global_params.mod precision.mod Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

bud_g.o: bud_g.f global_params.mod gas_Parameters.h gas_Global.h Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

bud_a.o: bud_a.f global_params.mod aer_Parameters.h aer_Global.h Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

bud_t.o: bud_t.f global_params.mod tot_Parameters.h tot_Global.h Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

bud_s_g.o: bud_s_g.f global_params.mod gas_Parameters.h gas_Global.h Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

bud_s_a.o: bud_s_a.f global_params.mod aer_Parameters.h aer_Global.h Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<

bud_s_t.o: bud_s_t.f global_params.mod tot_Parameters.h tot_Global.h Makefile
	$(FC) $(CFLAGS) -o $(@F) -c $<
