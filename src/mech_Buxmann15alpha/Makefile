#!/bin/csh
# List dependencies

#
# Copyright 2015-2017 Josue Bock
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

DEPEN_ALL = head.bud tail.bud master.spc master_gas.eqn
DEPEN_AER = aer_km.dat aer_mk.dat master_aqueous.eqn title.1 title.2 title.3 title.4

# Script to use
SC = ./make_kpp.sc

all: gas.f aer.f tot.f

gas.f: gas.k gas.def gas.eqn drive_gas.f $(DEPEN_ALL)
	$(SC) gas g

aer.f: aer.k aer.def aer_eqn.head drive_aer.f $(DEPEN_ALL) $(DEPEN_AER)
	$(SC) aer a

tot.f: tot.k tot.def tot_eqn12.head tot_eqn34.head drive_tot.f $(DEPEN_ALL) $(DEPEN_AER)
	$(SC) tot t

clean:
	rm -f bud_*.f gas.f aer.f tot.f *.h *.map
