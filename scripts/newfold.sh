#!/bin/bash
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

for i in plot*sc
do
  sed "s|data/Elise/Mistra_Elise_data/Elise/Mistra_Elise_|data/Elise/Mistra_Elise_|g" $i > tmp.sc
  sed "s|_2_unfinished_run||g" tmp.sc > $i
#  mv tmp.sc $i
  chmod u+x $i
done
