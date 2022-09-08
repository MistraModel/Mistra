#!/bin/csh
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

# script that writes a ferret script to plot the wanted meteorological params and executes it

# name of plot
# ------------
set pl_name = ("met")

# define number of runs and directories of netCDF files
# -----------------------------------------------------

# max number of runs is 5!
set n_runs = ("1")
set m_runs = ("1")
# set n_runs = ("2")
# set m_runs = ("1" "2")
# set n_runs = ("4")
# set m_runs = ("1" "2" "3" "4")
# set n_runs = ("5")
# set m_runs = ("1" "2" "3" "4" "5")
set runs = ("/data/Elise/Mistra_Elise_v008.4/")



# how to plot
# -----------

# 1 - evolution w/ time; 2 - vertical profile; 3 - contour plot; 4 - box model
# set pl_ty = ("1")
# set pl_ty = ("2")
set pl_ty = ("3")
# set pl_ty = ("4")

# times (h from model start)/altitudes (m); n_pl: number of heights/times to be overplotted
set n_pl   = ("2")
set m_pl   = ("1" "2")

# so far put levels and timesteps instead of height and time...
set pld = ("25" "75")
# set pld = ("48" "96")

# maximum and minimum heights if plot type = 2,3 [now in levels, later in m]
# set minheight = 2 - doesn't work!! set to 2 by hand
set maxheight = 150

# NOTE: if you chose more than 4 lines it'll be hard to tell the differences!

# unit of plot: 1 - unchanged
# --------------------------------------------------------
set pl_un = 1


# what to plot
# ------------

set n_spec = ("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11")
set species = ("u" "v" "w" "theta" "thetl" "temp" "q" "LWC" "rh" "p" "rho" )


# ===========================================================================
# HANDS OFF FROM BELOW!! NO USER SERVICABLE PARTS INSIDE!! DEVELOPPERS ONLY!!
# ===========================================================================

# make sure that each plotprog file is saved and none overwritten to be able to
# reproduce runs; this plotprog file name is put into plot

# numbers of plots per page
set pl_n = ("7")

# define "overflow" value for viewports
set pl_np = $pl_n
@ pl_np+=1

# define header
set n_hl=("1" "2" "3" "4" "5" "6" "7" "8")
set hl1=("! Description: Ferret plotprogram for Mistra netCDf output: meteorology\n")
set hl2=("! input unit is mol m-3\n")
set hl3=(" \n")
set hl4=("! reset everything\n")
set hl5=("\cancel mode verify ! like "unset echo" under unix\n")
set hl6=("cancel data/all\n")
set hl7=("cancel variable/all\n")
set hl8=("cancel symbol/all\n")

# write header
echo $hl1 $hl2 $hl3 $hl4 $hl6 $hl7 $hl8  >> pgtmp0

# write special definitions
echo "let q=1000.*xm1" >> pgtmp0
echo "let LWC=1000.*xm2" >> pgtmp0
echo "let rh=100.*feu" >> pgtmp0

# open files
foreach i ($m_runs)
    echo "use" ' "'  "$runs[$i]"'meteo.nc"' >> pgtmp0
end

# define timeaxis
# mtime as in netCDF is "local time", get "time since model start"
echo "let mmtime=mtime-mtime[l=1]" >> pgtmp0

# define viewports and page style
# 3x4
echo "define viewport /xlim=0.00,0.49 /ylim=0.65,1.00 VP1"  >> pgtmp0
echo "define viewport /xlim=0.00,0.49 /ylim=0.45,0.80 VP2"  >> pgtmp0
echo "define viewport /xlim=0.00,0.49 /ylim=0.25,0.60 VP3"  >> pgtmp0
echo "define viewport /xlim=0.00,0.49 /ylim=0.05,0.40 VP4"  >> pgtmp0
echo "define viewport /xlim=0.51,1.00 /ylim=0.65,1.00 VP5"  >> pgtmp0
echo "define viewport /xlim=0.51,1.00 /ylim=0.45,0.80 VP6"  >> pgtmp0
echo "define viewport /xlim=0.51,1.00 /ylim=0.25,0.60 VP7"  >> pgtmp0
echo "define viewport /xlim=0.51,1.00 /ylim=0.05,0.40 VP8"  >> pgtmp0

# set size of symbols
echo "ppl axlsze .16,.16"  >> pgtmp0
echo "ppl labset .16,.16,.16"  >> pgtmp0
# for legend font ppl command might have to be used instead of the normal "shade"

# output in metafile format
echo "SET MODE METAFILE:"$pl_name.1".plt" >> pgtmp0

# plot command is being composed of the different single steps:

# base plot command
if ($pl_ty == 1) then
# 1 - evolution w/ time;
 set plcmd = ("plot /vs /line /nolabel /k=")
 set plcmdo = ("plot /vs /line /nolabel /overlay /k=")
 set pll1 = ("k=")
 set lim1 = ("/vlimits=0:")
 set lim2 = (",l=@max")
 set xax = ("mmtimedd1 ,")
 set yax = (" ")
endif
if ($pl_ty == 2) then
# 2 - vertical profile
# set plcmd = ("plot /nolabel /k=$minheight:$maxheight /l=")
 set plcmd = ("plot /vs /line /nolabel /k=2:$maxheight /l=")
# set plcmdo = ("plot /nolabel /overlay /k=$minheight:$maxheight /l=")
 set plcmdo = ("plot /vs /line /nolabel /overlay /k=2:$maxheight /l=")
 set pll1 = ("l=")
 set lim1 = ("/hlimits=0:")
# set lim2 = (",k=$minheight:$maxheight@max")
 set lim2 = (",k=2:$maxheight@max")
 set xax = (" ")
 set yax = (",etadd1")
endif
if ($pl_ty == 3) then
# 3 - contour plot
# set plcmd = ("shade /nolabel /k=$minheight:")
 set plcmd = ("shade /nolabel /k=2:")
# set plcmdo = ("shade /nolabel /overlay /k=$minheight:")
 set plcmdo = ("shade /nolabel /overlay /k=2:")
 set pll1 = ("k=")
 set lim1 = (" ")
 set lim2 = (",l=@max")
 set pld=$maxheight
 set n_pl = ("1")
 set m_pl = ("1")
 set n_runs = ("1")
 set m_runs = ("1")
 set xax = (" ")
 set yax = (" ")
endif
# 4 - box model
if ($pl_ty == 4) then
# like evolution w/ time, but /k=2
 set plcmd = ("plot /nolabel /k=")
 set plcmdo = ("plot /nolabel /overlay /k=")
 set pll1 = ("k=")
 set lim1 = ("/vlimits=0:")
 set lim2 = (",l=@max")
 set pld=("2")
 set n_pl = ("1")
 set m_pl = ("1")
 set xax = (" ")
 set yax = (" ")
endif
# /title='"'$species[$k]'"'

# unit conversion: 1 - unchanged
# unchanged
if ($pl_un == 1) then
    set unit1a=1.
    set unit1b=1.
    set unit1c=1.
    set unit3a=1.
    set unit3b=1.
    set unit3c=1.
    set unit5a=1.
    set unit5b=1.
    set unit5c=1.
    set unit7a=1.
    set unit7b=1.
    set unit7c=1.
    set unit9a=1.
    set unit9b=1.
    set unit9c=1.
endif

# run info always on view port1 (VP1)
set countVP = 1
echo "set viewport VP$countVP" >> pgtmp0
# define linestyle master
set lstym=("" "/dash=(.1,.1,.1,.1)" "/dash=(.01,.1,.01,.1)" "/dash=(.3,.1,.3,.1)")
set lsty=("" "" "" "" "" "" "" "" "" "")
# define labels
set line1 =(" ")
set line2 =(" ")
set line3 =(" ")
set line4 =(" ")
set line5 =(" ")
set line6 =(" ")
set line7 =(" ")
set line8 =(" ")
set line9 =(" ")

# shorten directory names, max 32. chars
# n_runs - number of runs
# runs - complete directory names of runs
# TODO! maybe via awk and substr command - output inextra file, but how to read in from that file into var?

# one line explanation of run
if ($n_runs == 1) then
    set countl = 0
    foreach kk ($m_pl)
       @ countl+=1
       set line$countl =("$runs[1] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
    end
endif
if ($n_runs == 2) then
    set countl = 0
    foreach kk ($m_pl)
       @ countl+=1
       set line$countl =("$runs[1] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[2] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
    end
endif
if ($n_runs == 3) then
    set countl = 0
    foreach kk ($m_pl)
       @ countl+=1
       set line$countl =("$runs[1] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[2] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[3] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
    end
endif
if ($n_runs == 4) then
    set countl = 0
    foreach kk ($m_pl)
       @ countl+=1
       set line$countl =("$runs[1] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[2] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[3] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[4] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
    end
endif
if ($n_runs == 5) then
    set countl = 0
    foreach kk ($m_pl)
       @ countl+=1
       set line$countl =("$runs[1] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[2] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[3] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[4] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[5] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
    end
endif

# put unit on plot
# set unit_c=(" ")
# if ($pl_un == 1) set unit_c=("unit: mol m-3")
# if ($pl_un == 2) set unit_c=("mol mol1")
# if ($pl_un == 3) set unit_c=("unit: molec cm-3")
set unit_c=("")
set line9 = $unit_c[$pl_un]

# empty plot
echo "plot /i=0:100 /hlimits=0:100 /vlimits=0:100 /noaxis /nolabel 0" >> pgtmp0
# plot explanation
echo "label 10,110,  -1, 0, .25 @P1Meteorology" >> pgtmp0
echo "label 0, 90, -1, 0, .15 @P1$line1" >> pgtmp0
echo "plot /overlay /vs /nolabel /line=1 $lsty[1] {80,99,99}, {95,95,95}" >> pgtmp0
echo "label 0, 80, -1, 0, .15 @P1$line2" >> pgtmp0
if ($countl >= 2) echo "plot /overlay /vs /nolabel /line=2 $lsty[2] {80,99,99}, {85,85,85}" >> pgtmp0
echo "label 0, 70, -1, 0, .15 @P1$line3" >> pgtmp0
if ($countl >= 3) echo "plot /overlay /vs /nolabel /line=3 $lsty[3] {80,99,99}, {75,75,75}" >> pgtmp0
echo "label 0, 60, -1, 0, .15 @P1$line4" >> pgtmp0
if ($countl >= 4) echo "plot /overlay /vs /nolabel /line=4 $lsty[4] {80,99,99}, {65,65,65}" >> pgtmp0
echo "label 0, 50, -1, 0, .15 @P1$line5" >> pgtmp0
if ($countl >= 5) echo "plot /overlay /vs /nolabel /line=5 $lsty[5] {80,99,99}, {55,55,55}" >> pgtmp0
echo "label 0, 40, -1, 0, .15 @P1$line6" >> pgtmp0
if ($countl >= 6) echo "plot /overlay /vs /nolabel /line=6 $lsty[6] {80,99,99}, {45,45,45}" >> pgtmp0
echo "label 0, 30, -1, 0, .15 @P1$line7" >> pgtmp0
if ($countl >= 7) echo "plot /overlay /vs /nolabel /line=7 $lsty[7] {80,99,99}, {35,35,35}" >> pgtmp0
echo "label 0, 20, -1, 0, .15 @P1$line8" >> pgtmp0
if ($countl >= 8) echo "plot /overlay /vs /nolabel /line=8 $lsty[8] {80,99,99}, {25,25,25}" >> pgtmp0
echo "label 0, 10, -1, 0, .15 @P1$line9" >> pgtmp0
if ($countl >= 9) echo "plot /overlay /vs /nolabel /line=9 $lsty[9] {80,99,99}, {15,15,15}" >> pgtmp0

# counter for pages
set countPG = 1

# plot loop
foreach k ($n_spec)
# number of runs to be overplotted (= n_runs)
    @ countVP+=1
    echo "set viewport VP$countVP" >> pgtmp0
# counter, set viewport w/ counter: set viewport VP$count
# d1, d3, d5, .. have to be replaced with sed in final ferret file to give [d=1] etc

# 1 run
    if ($n_runs == 1) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
       end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       if ($pl_ty != 3) then
           echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`'  $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       else
           echo $plcmd$pld[1] $lstym[1] $species[$k]dd1/$unit1a/$unit1c >> pgtmp0
       endif
    endif
# 2 runs
    if ($n_runs == 2) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de2,$pll1$pld[$kk]$lim2]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd2/$unit3a/$unit3c$yax >> pgtmp0
    endif
# 3 runs
    if ($n_runs == 3) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de2,$pll1$pld[$kk]$lim2]/$unit2a/$unit2b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd2/$unit3a/$unit3c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd3/$unit5a/$unit5c$yax >> pgtmp0
    endif
# 4 runs
    if ($n_runs == 4) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de2,$pll1$pld[$kk]$lim2]/$unit2a/$unit2b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de4,$pll1$pld[$kk]$lim2]/$unit4a/$unit4b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd2/$unit3a/$unit3c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd3/$unit5a/$unit5c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd4/$unit7a/$unit7c$yax >> pgtmp0
    endif
# 5 runs
    if ($n_runs == 5) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de2,$pll1$pld[$kk]$lim2]/$unit2a/$unit2b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de4,$pll1$pld[$kk]$lim2]/$unit4a/$unit4b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de5,$pll1$pld[$kk]$lim2]/$unit5a/$unit5b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd2/$unit3a/$unit3c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd3/$unit5a/$unit5c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd4/$unit7a/$unit7c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd5/$unit9a/$unit9c$yax >> pgtmp0
    endif

#   add name of species
    echo "let label_y = 1.15*vmax" >> pgtmp0
    echo "if "'`'"$pl_ty eq 2 "'`'" then let label_y="'`'"1.05*eta[de1,l=1,k=$maxheight]"'`'" endif" >> pgtmp0
#    echo "if "'`'"$pl_ty eq 3 "'`'" then let label_y="'`'"1.05*eta[de1,l=1,k=$maxheight]"'`'" endif" >> pgtmp0
    echo "if "'`'"$pl_ty eq 3 "'`'" then let label_y=1.05*$maxheight endif" >> pgtmp0
    echo "label 0.," '`'"label_y"'`'",0,0,.2 @P1$species[$k]" >> pgtmp0
#   heights/times to be overplotted
    foreach kk ($m_pl)
        if ($n_runs == 1 && $kk != 1) echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       if ($n_runs == 2 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd2/$unit3a/$unit3c$yax >> pgtmp0
       endif
       if ($n_runs == 3 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd2/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd3/$unit5a/$unit5c$yax >> pgtmp0
       endif
       if ($n_runs == 4 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd2/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd3/$unit5a/$unit5c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd4/$unit7a/$unit7c$yax >> pgtmp0
       endif
       if ($n_runs == 5 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd2/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd3/$unit5a/$unit5c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd4/$unit7a/$unit7c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd5/$unit9a/$unit9c$yax >> pgtmp0
       endif
    end
    if ($countVP == $pl_np) then
#       next page
       echo "CANCEL MODE METAFILE" >> pgtmp0
       @ countPG+=1
       echo "SET MODE METAFILE:"$pl_name.$countPG".plt" >> pgtmp0
        set countVP = 1
        echo "SET WINDOW /clear" >> pgtmp0
#       put basic info again on VP1
        echo "set viewport VP$countVP" >> pgtmp0
        echo "plot /i=0:100 /hlimits=0:100 /vlimits=0:100 /noaxis /nolabel 0" >> pgtmp0
        echo "label 10,110,  -1, 0, .25 @P1Meteorology" >> pgtmp0
       echo "label 0, 90, -1, 0, .15 @P1$line1" >> pgtmp0
       echo "plot /overlay /vs /nolabel /line=1 $lsty[1] {80,99,99}, {95,95,95}" >> pgtmp0
       echo "label 0, 80, -1, 0, .15 @P1$line2" >> pgtmp0
       if ($countl >= 2) echo "plot /overlay /vs /nolabel /line=2 $lsty[2] {80,99,99}, {85,85,85}" >> pgtmp0
       echo "label 0, 70, -1, 0, .15 @P1$line3" >> pgtmp0
       if ($countl >= 3) echo "plot /overlay /vs /nolabel /line=3 $lsty[3] {80,99,99}, {75,75,75}" >> pgtmp0
       echo "label 0, 60, -1, 0, .15 @P1$line4" >> pgtmp0
       if ($countl >= 4) echo "plot /overlay /vs /nolabel /line=4 $lsty[4] {80,99,99}, {65,65,65}" >> pgtmp0
       echo "label 0, 50, -1, 0, .15 @P1$line5" >> pgtmp0
       if ($countl >= 5) echo "plot /overlay /vs /nolabel /line=5 $lsty[5] {80,99,99}, {55,55,55}" >> pgtmp0
       echo "label 0, 40, -1, 0, .15 @P1$line6" >> pgtmp0
       if ($countl >= 6) echo "plot /overlay /vs /nolabel /line=6 $lsty[6] {80,99,99}, {45,45,45}" >> pgtmp0
       echo "label 0, 30, -1, 0, .15 @P1$line7" >> pgtmp0
       if ($countl >= 7) echo "plot /overlay /vs /nolabel /line=7 $lsty[7] {80,99,99}, {35,35,35}" >> pgtmp0
       echo "label 0, 20, -1, 0, .15 @P1$line8" >> pgtmp0
       if ($countl >= 8) echo "plot /overlay /vs /nolabel /line=8 $lsty[8] {80,99,99}, {25,25,25}" >> pgtmp0
       echo "label 0, 10, -1, 0, .15 @P1$line9" >> pgtmp0
       if ($countl >= 9) echo "plot /overlay /vs /nolabel /line=9 $lsty[9] {80,99,99}, {15,15,15}" >> pgtmp0
    endif

end

echo "CANCEL MODE METAFILE" >> pgtmp0
# if total number of plots is multiple of "plots per page" there is an empty output file: ignore it!
if ($countVP == 1) @ countPG-=1

# finalize plot file
# d1, d3, d5, .. have to be replaced with sed in final ferret file to give [d=1] etc
sed 's/dd10/[d=10]/g' pgtmp0 > pgtmp0a
sed 's/dd1/[d=1]/g' pgtmp0a > pgtmp1
sed 's/dd2/[d=2]/g' pgtmp1 > pgtmp2
sed 's/dd3/[d=3]/g' pgtmp2 > pgtmp3
sed 's/dd4/[d=4]/g' pgtmp3 > pgtmp4
sed 's/dd5/[d=5]/g' pgtmp4 > pgtmp5
sed 's/dd6/[d=6]/g' pgtmp5 > pgtmp6
sed 's/dd7/[d=7]/g' pgtmp6 > pgtmp6a
sed 's/dd8/[d=8]/g' pgtmp6a > pgtmp6b
sed 's/dd9/[d=9]/g' pgtmp6b > pgtmp6c
sed 's/de10/d=10/g' pgtmp6c > pgtmp6d
sed 's/de1/d=1/g' pgtmp6d > pgtmp7
sed 's/de2/d=2/g' pgtmp7 > pgtmp8
sed 's/de3/d=3/g' pgtmp8 > pgtmp9
sed 's/de4/d=4/g' pgtmp9 > pgtmp10
sed 's/de5/d=5/g' pgtmp10 > pgtmp11
sed 's/de6/d=6/g' pgtmp11 > pgtmp12
sed 's/de7/d=7/g' pgtmp12 > pgtmp13
sed 's/de8/d=8/g' pgtmp13 > pgtmp14
sed 's/de9/d=9/g' pgtmp14 > pgtmp15

# name of plot: pl_name
mv -f pgtmp15 $pl_name.jnl

# plot

# source /soft/ferret_paths_RH9
# # ferret -batch $pl_name.ps -script $pl_name.jnl
ferret -script $pl_name.jnl

# determine number of meta print files
if ($countPG == 1) set metafiles="$pl_name.1.plt"
if ($countPG == 2) set metafiles="$pl_name.1.plt $pl_name.2.plt"
if ($countPG == 3) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt"
if ($countPG == 4) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt"
if ($countPG == 5) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt"
if ($countPG == 6) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt"
if ($countPG == 7) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt"
if ($countPG == 8) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt"
if ($countPG == 9) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt"
if ($countPG == 10) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt"
if ($countPG == 11) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt"
if ($countPG == 12) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt"
if ($countPG == 13) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt"


gksm2ps -p portrait -l cps -d cps -o $pl_name.pre.ps $metafiles

# add page numbering to ps file and delete empty pages
ps2ps $pl_name.pre.ps $pl_name.ps



# clean up ----
# temporary files
rm -rf pgtmp*
# prelim PS file
rm -f $pl_name.pre.ps
# meta print files
rm -f *.plt
