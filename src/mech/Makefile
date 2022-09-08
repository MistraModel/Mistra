# ======================================================================= #
# Copyright 2015-2017 the Authors
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

# List dependencies
DEPEN_ALL = head.bud tail.bud master.spc master_gas.eqn
DEPEN_AER = aer_km.dat aer_mk.dat master_aqueous.eqn title.1 title.2 title.3 title.4

# Script to use
SC = ./make_kpp.sc

# ---------------------------------------------------

all: gas.f aer.f tot.f

gas.f: gas.k gas.def gas.eqn drive_gas.f $(DEPEN_ALL)
	$(SC) gas g

aer.f: aer.k aer.def aer_eqn.head drive_aer.f $(DEPEN_ALL) $(DEPEN_AER)
	$(SC) aer a

tot.f: tot.k tot.def tot_eqn12.head tot_eqn34.head drive_tot.f $(DEPEN_ALL) $(DEPEN_AER)
	$(SC) tot t

clean:
	rm -f bud_*.f gas.f aer.f tot.f *.h *.map
