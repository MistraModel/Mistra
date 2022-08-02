#!/bin/csh
# script that writes a ferret script to plot the wanted aqueous phase params and executes it

# name of plot
# ------------
set pl_name = ($2)

# define number of runs and directories of netCDF files
# -----------------------------------------------------

# max number of runs is 5!
set n_runs = ("1")
set m_runs = ("1")
# set n_runs = ("2")
# set m_runs = ("1" "2")
# set n_runs = ("3")
# set m_runs = ("1" "2" "3")
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
set n_pl   = ("1")
# set m_pl   = ("1" "2" "3")
# set m_pl   = ("1" "2")
set m_pl   = ("1")

# so far put levels and timesteps instead of height and time...
# set pld = ("5" "25" "55" "75")
# set pld = ("1" "2" "3")
# set pld = ("40" "64")
set pld = ("25")
# set pld = ("144" "192")
# set pld = ("25" "75")
# set pld = ("50" "250" "550" "750")
# set pld = ("24" "48" "72" "96")
# set pld = ("240" "288")
# set pld = ("10" "20" "30")
# set pld = ("12" "36" "50" "74")
# set pld = ("48" "144" "200" "296")
# set pld = ("144" "296")

# maximum and minimum heights if plot type = 2,3 [now in levels, later in m]
# set minheight = 2 - doesn't work!! set to 2 by hand
# set maxheight = 150
set maxheight = 100

# NOTE: if you chose more than 4 lines it'll be hard to tell the differences!


# numbers of plots per page - better don't change this!
set pl_n = ("7")

# unit of plot: 1 - mol m-3; 2 - mol mol-1; 3 - molec cm-3; 4 - mol l-1
# ---------------------------------------------------------------------
set pl_un = 1


# what to plot
# ------------

# chose aqueous bin: 1 - 4, if cloud free run only 1 - 2
# set pl_bin = ("1")
# set pl_bin = ("2")
# set pl_bin = ("3")
# set pl_bin = ("4")
set pl_bin = $1


# max possible (incl. those commented in output):
# set n_spec = ("125")
# set species = ( "NO" "NO2" "HNO3" "NH3" "SO2" "H2SO4" "O3" "CH4" "C2H6" "C3H8" "ALKA" "ETHE" "ALKE" "AROM" "ACO2" "HCOOH" "ACTA" "HCHO" "ALD2" "H2O2" "ROOH" "CH3OOH" "HONO" "PAN" "TPAN" "KET" "CRES" "DIAL" "GLYX" "MGL" "NH4NO3" "HCl" "R3N2" "RAN2" "RAN1" "N2O5" "HNO4" "NO3" "DMS" "HOCl" "ClNO2" "ClNO3" "Cl2" "HBr" "HOBr" "BrNO2" "BrNO3" "Br2" "BrCl" "HI" "HOI" "I2O2" "INO2" "INO3" "I2" "ICl" "IBr" "CH3I" "CH2I2" "CH2ClI" "C3H7I" "DMSO" "CH3SO2" "CH3SO3" "CH3SO3H" "CO" "Cl2O2" "DMOO" "CH3S" "CH3SO" "CH3SO2H" "DMSO2" "CH2BrI" "CHBr2I" "C2H5I" "OH" "HO2" "DOM" "CH3OO" "MO2" "IO" "Cl" "Br" "CH3OH" "O2" "CO2" "Hp" "NH4p" "OHm" "CH2OHSO3m" "HSO3m" "SO32m" "SO4m" "SO42m" "HCO3m" "CO3m" "O2m" "NO2m" "NO3m" "Clm" "Cl2m" "HCOOm" "FE3p" "MN2p" "HSO4m" "Nap" "NO4m" "ClOm" "ClOHm" "Brm" "Br2m" "BrOm" "BrOHm" "BrCl2m" "Br2Clm" "CH3SO3m" "HSO5m" "SO3m" "SO5m" "Im" "IO2m" "IO3m" "ICl2m" "IBr2m" "CH3SO2m")
# species actually in output:
# set species = ("NO2" "HNO3" "NH3" "SO2" "H2SO4" "O3" "HCOOH" "HCHO" "H2O2" "CH3OOH" "HONO" "HCl" "HNO4" "NO3" "DMS" "HOCl" "Cl2" "HBr" "HOBr" "Br2" "BrCl" "HOI" "I2" "ICl" "IBr" "DMSO" "DMSO2" "OH" "HO2" "DOM" "CH3OO" "IO" "Cl" "Br" "CH3OH" "O2" "CO2" "Hp" "NH4p" "OHm" "CH2OHSO3m" "HSO3m" "SO32m" "SO4m" "SO42m" "HCO3m" "CO3m" "O2m" "NO2m" "NO3m" "Clm" "Cl2m" "HCOOm" "HSO4m" "Nap" "NO4m" "ClOm" "ClOHm" "Brm" "Br2m" "BrOm" "BrOHm" "BrCl2m" "Br2Clm" "CH3SO3m" "HSO5m" "SO3m" "SO5m" "Im" "IO2m" "IO3m" "ICl2m" "IBr2m" "CH3SO2m")
# the useful ones in useful order (w/o iodine):
#set species = ("NO2" "HNO3" "NO3m" "HNO4" "NO4m" "NO3" "HONO" "NO2m" "NH3" "NH4p" "Hp" "OHm" "O3" "O2" "O2m" "OH" "HO2" "H2O2" "CH3OO" "HCOOH" "HCOOm" "HCHO" "CH3OH" "CH3OOH" "CO2" "HCO3m" "CO3m" "SO2" "HSO3m" "SO32m" "H2SO4" "HSO4m" "SO42m" "SO3m" "SO4m" "HSO5m" "SO5m" "DMS" "DMSO" "DMSO2" "CH3SO2m" "CH3SO3m" "CH2OHSO3m" "HCl" "Clm" "HOCl" "ClOm" "Cl2" "Cl" "Cl2m" "ClOHm" "HBr" "Brm" "HOBr" "BrOm" "Br2" "BrCl" "Br" "Br2m" "BrOHm" "BrCl2m" "Br2Clm"  "Nap" "DOM" "SIV" "SVI"  "pH" "cw" "rc" "totNO3m" "Brm_Nap" "Clm_Nap")
#set n_spec = ("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "39" "40" "41" "42" "43" "44" "45" "46" "47" "48" "49" "50" "51" "52" "53" "54" "55" "56" "57" "58" "59" "60" "61" "62" "63" "64" "65" "66" "67" "68" "69" "70" "71" "72")
# the useful ones in useful order (w/ iodine):
  set species = ("NO2" "HNO3" "NO3m" "HNO4" "NO4m" "NO3" "HONO" "NO2m" "NH3" "NH4p" "Hp" "OHm" "O3" "O2" "O2m" "OH" "HO2" "H2O2" "CH3OO" "HCOOH" "HCOOm" "HCHO" "CH3OH" "CH3OOH" "CO2" "HCO3m" "CO3m" "SO2" "HSO3m" "SO32m" "H2SO4" "HSO4m" "SO42m" "SO3m" "SO4m" "HSO5m" "SO5m" "DMS" "DMSO" "DMSO2" "CH3SO2m" "CH3SO3m" "CH2OHSO3m" "HCl" "Clm" "HOCl" "ClOm" "Cl2" "Cl" "Cl2m" "ClOHm" "HBr" "Brm" "HOBr" "BrOm" "Br2" "BrCl" "Br" "Br2m" "BrOHm" "BrCl2m" "Br2Clm" "HOI" "ICl" "IBr" "I2" "IO" "Im" "IO2m" "IO3m" "ICl2m" "IBr2m" "ICl2m" "IBr2m" "Nap" "DOM" "SIV" "SVI"  "pH" "cw" "rc" "totNO3m" "Brm_Nap" "Clm_Nap")
 set n_spec = ("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "39" "40" "41" "42" "43" "44" "45" "46" "47" "48" "49" "50" "51" "52" "53" "54" "55" "56" "57" "58" "59" "60" "61" "62" "63" "64" "65" "66" "67" "68" "69" "70" "71" "72" "73" "74" "75" "76" "77" "78" "79" "80" "81" "82" "83" "84")

# set species = ("totNO3m" "Hp" "HCO3m" "HCl" "SVI" "Clm" "Brm" "Nap" "cw" "pH")
# set n_spec = ("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
# set species = ("SVI"  "pH" "cw" "totNO3m" "Hp" "HCO3m" "Clm")
# set n_spec = ("1" "2" "3" "4" "5" "6" "7")

# set species = ("NO2" "HNO3" "NH3" "SO2" "H2SO4" "O3" "HCOOH" "HCHO" "H2O2" "CH3OOH" "HONO" "HNO4" "NO3")
# set n_spec = ("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13")


# ===========================================================================
# HANDS OFF FROM BELOW!! NO USER SERVICABLE PARTS INSIDE!! DEVELOPPERS ONLY!!
# ===========================================================================

# make sure that each plotprog file is saved and none overwritten to be able to
# reproduce runs; this plotprog file name is put into plot

# define "overflow" value for viewports
set pl_np = $pl_n
@ pl_np+=1

# define header
set n_hl=("1" "2" "3" "4" "5" "6" "7" "8")
set hl1=("! Description: Ferret plotprogram for Mistra netCDf output: aqueous phase\n")
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
echo "let totNO3m=HNO3+NO3m"  >> pgtmp0
echo "let SVI=H2SO4+HSO4m+SO42m"  >> pgtmp0
echo "let SIV=SO2+HSO3m+SO32m"  >> pgtmp0
echo "let pH=(-1.)*log(hp/(cw*1000.))"  >> pgtmp0
# the factor is to account for the diff between "real" Na+ and model Na+ which is the
# sum of additional positive charges
echo "let Brm_Nap=Brm/Nap/.806" >> pgtmp0
echo "let Clm_Nap=Clm/Nap/.806" >> pgtmp0

# open files
foreach i ($m_runs)
    echo "use" ' "'  "$runs[$i]"'aq.nc"' '; use "'  "$runs[$i]"'meteo.nc"' >> pgtmp0
    echo $i
end

# define timeaxis
# mtime as in netCDF is "local time", get "time since model start"
echo "let mmtime=mtime-mtime[l=1]" >> pgtmp0

# define viewports and page style
# 3x4
echo "define viewport /xlim=0.02,0.475 /ylim=0.69,1.00 VP1"  >> pgtmp0
echo "define viewport /xlim=0.02,0.475 /ylim=0.49,0.80 VP2"  >> pgtmp0
echo "define viewport /xlim=0.02,0.475 /ylim=0.29,0.60 VP3"  >> pgtmp0
echo "define viewport /xlim=0.02,0.475 /ylim=0.09,0.40 VP4"  >> pgtmp0
echo "define viewport /xlim=0.49,0.955 /ylim=0.69,1.00 VP5"  >> pgtmp0
echo "define viewport /xlim=0.49,0.955 /ylim=0.49,0.80 VP6"  >> pgtmp0
echo "define viewport /xlim=0.49,0.955 /ylim=0.29,0.60 VP7"  >> pgtmp0
echo "define viewport /xlim=0.49,0.955 /ylim=0.09,0.40 VP8"  >> pgtmp0

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
# set plcmd = ("plot /vs /line /nolabel /i=$pl_bin /k=")
# set plcmdo = ("plot /vs /line /nolabel /i=$pl_bin /overlay /k=")
 set plcmd = ("plot /vs /line /nolabel /k=")
 set plcmdo = ("plot /vs /line /nolabel /overlay /k=")
 set pll1 = ("k=")
 set lim1 = ("/vlimits=0:")
 set lim2 = (",l=@max")
 set bin = ("i=$pl_bin")
 set bin2 = (",i=$pl_bin")
 set xax = ("mmtimedd2 ,")
 set yax = (" ")
endif
if ($pl_ty == 2) then
# 2 - vertical profile
# set plcmd = ("plot /nolabel /i=$pl_bin /k=$minheight:$maxheight /l=")
 set plcmd = ("plot /vs /line /nolabel /k=2:$maxheight /l=")
# set plcmdo = ("plot /nolabel /overlay /i=$pl_bin /k=$minheight:$maxheight /l=")
 set plcmdo = ("plot /vs /line /nolabel /overlay /k=2:$maxheight /l=")
 set pll1 = ("l=")
 set lim1 = ("/hlimits=0:")
# set lim2 = (",k=$minheight:$maxheight@max")
 set lim2 = (",k=2:$maxheight@max")
 set bin = ("i=$pl_bin")
 set bin2 = (",i=$pl_bin")
 set xax = (" ")
 set yax = (",etadd2")
endif
if ($pl_ty == 3) then
# 3 - contour plot
echo "define axis /from_data /z /name=zaxis eta[d=2,l=1,k=2:$maxheight]" >> pgtmp0
# echo "define axis /from_data /z /name=zaxis eta[d=2,l=1,k=2:$maxheight]+1635." >> pgtmp0
# set size of axis  labels (for some strange reasons, the key labels have to be defined right before plot):
echo "ppl axlsze .16,.16" >> pgtmp0
# set plcmd = ("shade /nolabel /i=$pl_bin /k=$minheight:")
# set plcmd = ("shade /nolabel /i=$pl_bin /k=2:")
 set plcmd = ("shade /set_up /nolabel /i=$pl_bin /k=2:")
# set plcmdo = ("shade /nolabel /i=$pl_bin /overlay /k=$minheight:")
# set plcmdo = ("shade /nolabel /i=$pl_bin /overlay /k=2:")
 set plcmdo = ("shade /set_up /nolabel /i=$pl_bin /overlay /k=2:")
 set pll1 = ("k=")
 set lim1 = (" ")
 set lim2 = (",l=@max")
 set pld=$maxheight
 set n_pl = ("1")
 set m_pl = ("1")
 set n_runs = ("1")
 set m_runs = ("1")
 set bin = ("i=$pl_bin")
 set bin2 = (",i=$pl_bin")
 set xax = (" ")
 set yax = (" ")
# define axes
# echo " define axis/from_data/z/name=zaxis eta[l=1,k=2:150,d=2]+1635."  >> pgtmp0
echo " define axis/from_data/z/name=zaxis eta[l=1,k=2:150,d=2]"  >> pgtmp0
# echo " define axis/from_data/t/name=taxis mtime[d=2,k=1]"  >> pgtmp0
echo " define axis/from_data/t/name=taxis mtime[d=2,k=1]"  >> pgtmp0
endif
# 4 - box model
if ($pl_ty == 4) then
# like evolution w/ time, but /k=2
 set plcmd = ("plot /nolabel /i=$pl_bin /k=")
 set plcmdo = ("plot /nolabel /i=$pl_bin /overlay /k=")
 set pll1 = ("k=")
 set lim1 = ("/vlimits=0:")
 set lim2 = (",l=@max")
 set bin = ("i=$pl_bin")
 set bin2 = (",i=$pl_bin")
 set pld=("2" "0" "0" "0")
 set n_pl = ("1")
 set m_pl = ("1")
# set xax = ("mmtimedd2 ,")
 set xax = (" ")
 set yax = (" ")

endif
# /title='"'$species[$k]'"'

# unit conversion: 1 - mol m-3; 2 - mol mol-1; 3 - molec cm-3; 4 - mol l-1
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
# mol m-3 --> mol mol-1: 1./(rho / 29.e-3)
if ($pl_un == 2) then
    set unit1a=34.4827
# isse probblem hier: nix aendere pld[1] fuer andere hoehe, kann sein skalierprobblem auch??
    set unit1b="rho[de2,$pll1$pld[1]$lim2]"
#    set unit1b=("rho[de2,$pll1$pld[1]$lim2]" "rho[de2,$pll1$pld[2]$lim2]" "rho[de2,$pll1$pld[3]$lim2]" "rho[de2,$pll1$pld[4]$lim2]")
    set unit1c="rhodd2"
    set unit3a=34.4827
    set unit3b="rho[de4,$pll1$pld[1]$lim2]"
    set unit3c="rhodd4"
    set unit5a=34.4827
    set unit5b="rho[de6,$pll1$pld[1]$lim2]"
    set unit5c="rhodd6"
    set unit7a=34.4827
    set unit7b="rho[de8,$pll1$pld[1]$lim2]"
    set unit7c="rhodd8"
    set unit9a=34.4827
    set unit9b="rho[de10,$pll1$pld[1]$lim2]"
    set unit9c="rhodd10"
endif
# mol m-3 --> molec cm-3: 1./(1.e6 / 6.022e23)
if ($pl_un == 3) then
    set unit1a=1.6606e-18
    set unit1b=1.
    set unit1c=1.
    set unit3a=1.6606e-18
    set unit3b=1.
    set unit3c=1.
    set unit5a=1.6606e-18
    set unit5b=1.
    set unit5c=1.
    set unit7a=1.6606e-18
    set unit7b=1.
    set unit7c=1.
    set unit9a=1.6606e-18
    set unit9b=1.
    set unit9c=1.
endif
# mol m-3 --> mol l-1
if ($pl_un == 4) then
    set unit1a=1.e3
    set unit1b="cwde1,$pll1$pld[1]$lim2,$bin]"
# maybe introduce unitXd and split "b" in 2 to be able to put $kk in?
    set unit1c="cwde1,$bin]"
    set unit3a=1.e3
    set unit3b="cwde3,$pll1$pld[1]$lim2,$bin]"
    set unit3c="cwde3,$bin]"
    set unit5a=1.e3
    set unit5b="cwde5,$pll1$pld[1]$lim2,$bin]"
    set unit5c="cwde5,$bin]"
    set unit7a=1.e3
    set unit7b="cwde7,$pll1$pld[1]$lim2,$bin]"
    set unit7c="cwde7,$bin]"
    set unit9a=1.e3
    set unit9b="cwde9,$pll1$pld[1]$lim2,$bin]"
    set unit9c="cwde9,$bin]"
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
# if ($pl_un == 4) set unit_c=("unit: mol l-1")
set unit_c=("mol/m3" "mol/mol" "molec/cm3" "mol/l")
set line9 = $unit_c[$pl_un]

# empty plot
echo "plot /i=0:100 /hlimits=0:100 /vlimits=0:100 /noaxis /nolabel 0" >> pgtmp0
# plot explanation
echo "label 10,110,  -1, 0, .25 @P1Aqueous phase" >> pgtmp0
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
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2,$bin]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
       end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       if ($pl_ty != 3) then
           echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`'  $xax$species[$k]"[de1,$bin]"/$unit1a/$unit1c$yax >> pgtmp0
       else
#           echo $plcmd$pld[1] $lstym[1] $species[$k]dd1/$unit1a/$unit1c >> pgtmp0
#axzx ,gz=zaxis@asn,gt=taxis@asn
# define new variable in order to be able to use "nice" axes
           echo "let spe_ax $species[$k]axzx " >> pgtmp0
           echo "let vmax=eta[l=1,d=2,k=maxheight]" >> pgtmp0
#           echo $plcmd$pld[1] $lstym[1] $species[$k]"[d=1,gz=zaxis@asn]"/$unit1a/$unit1c >> pgtmp0
           echo $plcmd$pld[1] $lstym[1] spe_ax/$unit1a/$unit1c >> pgtmp0
#         echo "ppl shakey 1,1,0.13,2,4,8" >> pgtmp0
           echo "ppl shakey 1,1,0.14,2,4,8,,,," >> pgtmp0
           echo "ppl shade" >> pgtmp0
       endif
    endif
# 2 runs
    if ($n_runs == 2) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo $kk >> /dev/null
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2,$bin]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2,$bin]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]"[de1,$bin]"/$unit1a/$unit1c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]"[de3,$bin]"/$unit3a/$unit3c$yax >> pgtmp0
    endif
# 3 runs
    if ($n_runs == 3) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2,$bin]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2,$bin]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de5,$pll1$pld[$kk]$lim2,$bin]/$unit5a/$unit5b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]"[de1,$bin]"/$unit1a/$unit1c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]"[de3,$bin]"/$unit3a/$unit3c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]"[de5,$bin]"/$unit5a/$unit5c$yax >> pgtmp0
    endif
# 4 runs
    if ($n_runs == 4) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2,$bin]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2,$bin]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de5,$pll1$pld[$kk]$lim2,$bin]/$unit5a/$unit5b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de7,$pll1$pld[$kk]$lim2,$bin]/$unit7a/$unit7b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]"[de1,$bin]"/$unit1a/$unit1c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]"[de3,$bin]"/$unit3a/$unit3c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]"[de5,$bin]"/$unit5a/$unit5c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]"[de7,$bin]"/$unit7a/$unit7c$yax >> pgtmp0
    endif
# 5 runs
    if ($n_runs == 5) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2,$bin]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2,$bin]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de5,$pll1$pld[$kk]$lim2,$bin]/$unit5a/$unit5b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de7,$pll1$pld[$kk]$lim2,$bin]/$unit7a/$unit7b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de9,$pll1$pld[$kk]$lim2,$bin]/$unit9a/$unit9b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd5/$unit5a/$unit5c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd7/$unit7a/$unit7c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd9/$unit9a/$unit9c$yax >> pgtmp0
    endif

#   add name of species
    echo "let label_y = 1.05*vmax" >> pgtmp0
    echo "if "'`'"$pl_ty eq 2 "'`'" then let label_y="'`'"1.05*eta[de2,l=1,k=$maxheight]"'`'" endif" >> pgtmp0
#    echo "if "'`'"$pl_ty eq 3 "'`'" then let label_y="'`'"1.05*eta[de2,l=1,k=$maxheight]"'`'" endif" >> pgtmp0
#    echo "if "'`'"$pl_ty eq 3 "'`'" then let label_y=1.05*$maxheight endif" >> pgtmp0
#    echo "if "'`'"$pl_ty eq 3 "'`'" then let label_y=1.02*(eta[d=2,l=1,k=$maxheight]+1635.) endif" >> pgtmp0
    echo "if "'`'"$pl_ty eq 3 "'`'" then let label_y=1.05*(eta[d=2,l=1,k=$maxheight]) endif" >> pgtmp0
    echo "label 0.," '`'"label_y"'`'",0,0,.2 @P1$species[$k]" >> pgtmp0
#   heights/times to be overplotted
    foreach kk ($m_pl)
        if ($n_runs == 1 && $kk != 1) echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de1,$bin]"/$unit1a/$unit1c$yax >> pgtmp0
       if ($n_runs == 2 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de1,$bin]"/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de3,$bin]"/$unit3a/$unit3c$yax >> pgtmp0
       endif
       if ($n_runs == 3 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de1,$bin]"/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de5,$bin]"/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de3,$bin]"/$unit5a/$unit5c$yax >> pgtmp0
       endif
       if ($n_runs == 4 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de1,$bin]"/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de3,$bin]"/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de5,$bin]"/$unit5a/$unit5c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de7,$bin]"/$unit7a/$unit7c$yax >> pgtmp0
       endif
       if ($n_runs == 5 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de1,$bin]"/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de3,$bin]"/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de5,$bin]"/$unit5a/$unit5c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de7,$bin]"/$unit7a/$unit7c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]"[de9,$bin]"/$unit9a/$unit9c$yax >> pgtmp0
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
        echo "label 10,110,  -1, 0, .25 @P1Aqueous phase" >> pgtmp0
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
sed 's/axzx/[d=1,gz=zaxis\@asn,gt=taxis\@asn]/g' pgtmp0 > pgtmp0a
sed 's/dd10/[d=10]/g' pgtmp0a > pgtmp0b
sed 's/dd1/[d=1]/g' pgtmp0b > pgtmp1
sed 's/dd2/[d=2]/g' pgtmp1 > pgtmp2
sed 's/dd3/[d=3]/g' pgtmp2 > pgtmp3
sed 's/dd4/[d=4]/g' pgtmp3 > pgtmp4
sed 's/dd5/[d=5]/g' pgtmp4 > pgtmp5
sed 's/dd6/[d=6]/g' pgtmp5 > pgtmp6
sed 's/dd7/[d=7]/g' pgtmp6 > pgtmp6a
sed 's/dd8/[d=8]/g' pgtmp6a > pgtmp6b
sed 's/dd9/[d=9]/g' pgtmp6b > pgtmp6c
sed 's/cwde/cw[de/g' pgtmp6c > pgtmp6d
sed 's/de10/d=10/g' pgtmp6d > pgtmp6e
sed 's/de1/d=1/g' pgtmp6e > pgtmp7
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
if ($countPG == 14) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt"
if ($countPG == 15) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt"


gksm2ps -p portrait -l cps -d cps -o $pl_name.pre.ps $metafiles

# add page numbering to ps file and delete empty pages
ps2ps $pl_name.pre.ps $pl_name.ps



# clean up ----
# temporary files
rm -f pgtmp*
# prelim PS file
rm -f $pl_name.pre.ps
# meta print files
rm -f *.plt
