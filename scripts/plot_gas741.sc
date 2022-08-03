#!/bin/csh
# script that writes a ferret script to plot the wanted gas phase params and executes it

# note: location of label currently set by hand for contour plots !! don't use "maxheight" to calculate it, do it "differently"...
# note: x-axis offset done by hand as model run was restart, see: "define axis/from_data/t/name=taxis.."

# name of plot
# ------------
# do NOT use "ave", "dd1", "de1", .. as part of plot title
set pl_name = ("v741_gc")

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

# 1 - evolution w/ time; 2 - vertical profile; 3 - contour plot; 4 - box model; 5 - BL average
# set pl_ty = ("1")
# set pl_ty = ("2")
set pl_ty = ("3")
# set pl_ty = ("4")
# set pl_ty = ("5")


# times (h from model start)/altitudes (m); n_pl: number of heights/times to be overplotted
set n_pl   = ("2")
set m_pl   = ("1" "2")
# set m_pl   = ("1" )

# so far put levels and timesteps instead of height and time...
# set pld = ("240" "288")
# set pld = ("4" "6" "8")
# set pld = ("25" "75")
set pld = ("25" "75")
# set pld = ("25")
# set pld = ("144" "192")

# maximum and minimum heights if plot type = 2,3,5 [now in levels, later in m]
# set minheight = 2 - doesn't work!! set to 2 by hand
# set maxheight = 150
# set maxheight = 125
# set maxheight = 90
set maxheight = 80

# start and end of run - for pl_ty=1 only, change "set lim1" below if you want whole run; this feature has
#  been introduce to be able to plot runs of different lengths
# set istart = 0
# set iend   = 346


# NOTE: if you chose more than 4 lines it'll be hard to tell the differences!

# unit of plot: 1 - mol m-3; 2 - mol mol-1; 3 - molec cm-3
# --------------------------------------------------------
set pl_un = 2


# what to plot
# ------------


# w/ iodine:
set n_spec = ("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "39" "40" "41" "42" "43" "44" "45" "46" "47" "48" "49" "50" "51" "52" "53" "54" "55" "56" "57" "58" "59" "60" "61" "62" "63" "64" "65" "66"  "67" "68" "69" "70" "71" "72" "73" "74" "75" "76" "77" "78" "79" "80" "81" "82" "83" "84" "85" "86" "87" "88" "89" "90" "91" "92" "93" "94" "95" "96" "97" "98" "99" "100" "101" "102" "103" "104" "105" "106" "107" "108" "109" "110" "111" "112 " "113" "114" "115" "116" "117" "118" "119" "120" "121" "122" "123" "124" "125" "126" "127" "128" "129" "130" "131" "132" "133" "134" "135" "136" "137")
# "138")
set species = ("NO" "NO2" "NONO2" "NOx" "N2O5" "NO3" "HNO3" "HNO4" "HONO" "PAN" "TPAN" "NOy" "NH3" "O3" "O1D" "O3P" "H2" "OH" "HO2" "H2O2" "CO2" "CO" "CH4" "C2H6" "C2H4" "ALKA" "ALKE" "CH3OH" "C2H5OH" "HCOOH" "CH3COOH" "HCHO" "ALD2" "ROOH" "HOCH2O2" "CH3CO2" "CH3O2" "C2H5O2" "EO2" "KO2" "R3O2" "RAO2" "TO2" "TCO3" "ZO2" "PO2" "CH2O2" "CH3CHO2" "PRN1" "AROM" "KET" "CRES" "DIAL" "CHOCHO" "CH3COCHO" "R3N2" "RAN2" "RAN1" "SO2" "SO3" "HOSO2" "H2SO4" "DMS" "DMSO" "DMSO2" "DMOO" "CH3S" "CH3SO" "CH3SO2" "CH3SO3" "MSIA" "MSA" "SPAN" "NHS" "SOR" "HCl" "HOCl" "ClONO" "ClNO2" "ClNO3" "Cl2" "Cl" "ClO" "OClO" "Cl2O2" "ClO3" "Cl2O3" "RCl" "ClRO2" "Clx" "Cltot" "HOClClx" "HBr" "HOBr" "BrNO2" "BrNO3" "Br2" "BrCl" "Br" "BrO" "Br2O" "RBr" "BrRO2" "Brx" "Brtot" "HOBrBrx"  "BrOBrx" "Br2Brx" "HI" "HIO3" "HOI" "INO" "INO2" "INO3" "I2" "ICl" "IBr" "Irad" "IO" "OIO" "I2O2" "I2O3" "I2O" "I2O4" "I2O5" "CH3I" "CH2I2" "CH2ClI" "C3H7I" "CH2BrI" "CHBr2I" "C2H5I" "IRO2" "Ix" "Iorg" "Itot" "XOR" )
# "NUCV")

# additional diagnostics
# "logClx" "logCltot" "ClOSO2" "BrOClO" "SO2HCl" "HClCltot" "HOClClx" "ClNO2Clx" "ClNO3Clx" "Cl2Clx" "ClClx" "ClOClx" "OClOClx" "Cl2O2Clx" "ClxCltot" "BrClClx"
#  "NO2_ClNO3"  "logBrx" "logBrtot" "BrOSO2" "BrOClO" "logSO2HBr" "HBrBrtot" "HOBrBrx" "BrNO2Brx" "BrNO3Brx" "Br2Brx" "BrClBrx" "BrBrx" "BrOBrx"  "BrxBrtot"



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
set hl1=("! Description: Ferret plotprogram for Mistra netCDf output: gas phase\n")
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
# C
echo "let HCOOH=ACO2" >> pgtmp0
echo "let C2H4=ETHE" >> pgtmp0
echo "let CH3CO2=MCO3" >> pgtmp0
echo "let CH3O2=MO2" >> pgtmp0
echo "let C2H5O2=ETO2" >> pgtmp0
echo "let CH2O2=CHO2" >> pgtmp0
echo "let CH3CHO2=CRO2" >> pgtmp0
echo "let CH3COOH=ACTA" >> pgtmp0
echo "let CHOCHO=GLYX" >> pgtmp0
echo "let CH3COCHO=MGLY" >> pgtmp0
echo "let HOCH2O2=AHO2" >> pgtmp0
# N
echo "let NOx=NO+NO2" >> pgtmp0
echo "let NONO2=NO/NO2" >> pgtmp0
echo "let logNOx=log(NOx)" >> pgtmp0
echo "let NOy=NOx+NO3+2*N2O5+HNO3+HNO4+PAN+HONO+ClNO2+ClNO3+BrNO2+BrNO3" >> pgtmp0
echo "let logNOy=log(NOy)" >> pgtmp0
# S
echo "let logH2SO4=log(H2SO4)" >> pgtmp0
echo "let logSO2=log(SO2)" >> pgtmp0
# Cl
echo "let Clx=HOCl+ClONO+ClNO2+ClNO3+2*Cl2+Cl+ClO+OClO+2*Cl2O2+ClO3+2*Cl2O3+BrCl+RCl+ClRO2"  >> pgtmp0
echo "let logClx=log(Clx)" >> pgtmp0
echo "let Cltot=HCl+HOCl+ClONO+ClNO2+ClNO3+2*Cl2+Cl+ClO+OClO+2*Cl2O2+ClO3+2*Cl2O3+BrCl+RCl+ClRO2"  >> pgtmp0
echo "let logCltot=log(Cltot)" >> pgtmp0
echo "let logHCl=log(HCl)" >> pgtmp0
echo "let HClCltot=HCl/Cltot"  >> pgtmp0
echo "let HOClClx=HOCl/Clx"  >> pgtmp0
echo "let ClNO2Clx=ClNO2/Clx"  >> pgtmp0
echo "let ClNO3Clx=ClNO3/Clx"  >> pgtmp0
echo "let ClNO2NOy=ClNO2/NOy"  >> pgtmp0
echo "let ClNO3NOy=ClNO3/NOy"  >> pgtmp0
echo "let Cl2Clx=2.*Cl2/Clx"  >> pgtmp0
echo "let BrClClx=BrCl/Clx"  >> pgtmp0
echo "let ClClx=Cl/Clx"  >> pgtmp0
echo "let ClOClx=ClO/Clx"  >> pgtmp0
echo "let OClOClx=OClO/Clx"  >> pgtmp0
echo "let Cl2O2Clx=Cl2O2/Clx"  >> pgtmp0
echo "let ClxCltot=Clx/Cltot"  >> pgtmp0
echo "let ClOSO2=ClO/SO2"  >> pgtmp0
echo "let SO2HCl=SO2/HCl"  >> pgtmp0
echo "let NO2_ClNO3=NO2/ClNO3" >> pgtmp0
# Br
echo "let Brx=HOBr+BrNO2+BrNO3+2*Br2+BrCl+Br+BrO+2*Br2O+RBr+BrRO2"  >> pgtmp0
echo "let logBrx=log(Brx)" >> pgtmp0
echo "let Brtot=HBr+HOBr+BrNO2+BrNO3+2*Br2+BrCl+Br+BrO+2*Br2O+RBr+BrRO2"  >> pgtmp0
echo "let logBrtot=log(Brtot)" >> pgtmp0
echo "let logHBr=log(HBr)" >> pgtmp0
echo "let HBrBrtot=HBr/Brtot"  >> pgtmp0
echo "let HOBrBrx=HOBr/Brx"  >> pgtmp0
echo "let BrNO2Brx=BrNO2/Brx"  >> pgtmp0
echo "let BrNO3Brx=BrNO3/Brx"  >> pgtmp0
echo "let BrNO2NOy=BrNO2/NOy"  >> pgtmp0
echo "let BrNO3NOy=BrNO3/NOy"  >> pgtmp0
echo "let Br2Brx=2.*Br2/Brx"  >> pgtmp0
echo "let BrClBrx=BrCl/Brx"  >> pgtmp0
echo "let BrBrx=Br/Brx"  >> pgtmp0
echo "let BrOBrx=BrO/Brx"  >> pgtmp0
echo "let BrxBrtot=Brx/Brtot"  >> pgtmp0
echo "let BrOClO=BrO/ClO"  >> pgtmp0
echo "let BrOSO2=BrO/SO2"  >> pgtmp0
echo "let logSO2HBr=log(SO2/HBr)"  >> pgtmp0
# I
echo "let Ix=HIO3+HOI+INO+INO2+INO3+2.*I2+ICl+IBr+Irad+IO+OIO+2.*I2O2+2*I2O3+2*I2O+2*I2O4+2*I2O5+IRO2"  >> pgtmp0
echo "let Iy=HI+HIO3+HOI+INO+INO2+INO3+2.*I2+ICl+IBr+Irad+IO+OIO+2.*I2O2+2*I2O3+2*I2O+2*I2O4+2*I2O5+IRO2"  >> pgtmp0
echo "let Iorg=CH3I+2.*CH2I2+CH2ClI+C3H7I+CH2BrI+CHBr2I+C2H5I"  >> pgtmp0
echo "let Itot=Iy+Iorg"  >> pgtmp0

# open files
foreach i ($m_runs)
    echo "use" ' "'  "$runs[$i]"'gas.nc"' '; use "'  "$runs[$i]"'meteo.nc"' >> pgtmp0
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
 set plcmd = ("plot /vs /line /nolabel /k=")
 set plcmdo = ("plot /vs /line /nolabel /overlay /k=")
 set pll1 = ("k=")
#  set lim1 = ("/hlimits=$istart:$iend /vlimits=0:")
 set lim1 = ("/vlimits=0:")
 set lim2 = (",l=@max")
 set xax = ("mmtimedd2 ,")
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
 set yax = (",etadd2")
endif
if ($pl_ty == 3) then
# 3 - contour plot
 echo "define axis /from_data /z /name=zaxis eta[d=2,l=1,k=2:$maxheight]" >> pgtmp0
 echo " define axis/from_data/t/name=taxis mtime[d=2,k=1]"  >> pgtmp0
# echo " define axis/from_data/t/name=taxis mtime[d=2,k=1]"  >> pgtmp0
# set plcmd = ("shade /nolabel /k=$minheight:")
#  set plcmd = ("shade /nolabel /k=2:")
 set plcmd = ("shade /set_up /nolabel /k=2:")
# set plcmdo = ("shade /nolabel /overlay /k=$minheight:")
# set plcmdo = ("shade /nolabel /overlay /k=2:")
 set plcmdo = ("shade /set_up /nolabel /overlay /k=2:")
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
if ($pl_ty == 5) then
# 5 - BL average
 set plcmd = ("plot /vs /line /nolabel ")
 set plcmdo = ("plot /vs /line /nolabel /overlay")
 set pll1 = ("k=")
 set n_pl   = ("1")
 set m_pl   = ("1")
 set pld = ($maxheight)
#  set lim1 = ("/hlimits=$istart:$iend /vlimits=0:")
 set lim1 = ("/vlimits=0:")
 set lim2 = (",l=@max")
 set xax = ("mmtimed2k2 ,")
 set yax = (" ")
endif
# /title='"'$species[$k]'"'

# unit conversion: 1 - mol m-3; 2 - mol mol-1; 3 - molec cm-3
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
    set unit1b="rho[de2,$pll1$pld[1]$lim2]"
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

# run info always on view port1 (VP1)
set countVP = 1
echo "set viewport VP$countVP" >> pgtmp0
# define linestyle master
set lstym=("" "/dash=(.1,.1,.1,.1)" "/dash=(.3,.1,.3,.1)" "/dash=(.01,.1,.01,.1)" "/dash=(.03,.1,.03,.1)" "" "")
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
set unit_c=("mol/m3" "mol/mol" "molec/cm3")
set line9 = $unit_c[$pl_un]

# empty plot
echo "plot /i=0:100 /hlimits=0:100 /vlimits=0:100 /noaxis /nolabel 0" >> pgtmp0
# plot explanation
echo "label 10,110,  -1, 0, .25 @P1Gas phase" >> pgtmp0
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
       if (($pl_ty == 1) || ($pl_ty == 2) || ($pl_ty == 4)) then
           echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`'  $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       endif
       if ($pl_ty == 3) then
# "old" and working version with metres as y-axis
#           echo "let vmax=eta[l=1,d=2,k=maxheight]" >> pgtmp0
#           echo $plcmd$pld[1] $lstym[1] $species[$k]"[d=1,gz=zaxis@asn]"/$unit1a/$unit1c >> pgtmp0
#           echo "ppl shakey 1,1,0.14,2,4,8,,,," >> pgtmp0
#           echo "ppl shade" >> pgtmp0
# new, "posh" version with time as x-axis and metres as y-axis:
# axzx ,gz=zaxis@asn,gt=taxis@asn
# define new variable in order to be able to use "nice" axes
           echo "let spe_ax $species[$k]axzx " >> pgtmp0
           if ($unit1c == 1.) then
              echo "let spe_ax_rho 1." >> pgtmp0
           else
              echo "let spe_ax_rho axzx2" >> pgtmp0
           endif
#           echo $plcmd$pld[1] $lstym[1] spe_ax/$unit1a/$unit1c >> pgtmp0
           echo $plcmd$pld[1] $lstym[1] spe_ax/$unit1a/spe_ax_rho  >> pgtmp0
           echo "ppl shakey 1,1,0.14,2,4,8,,,," >> pgtmp0
           echo "ppl shade" >> pgtmp0
       endif
       if ($pl_ty == 5) then
# first define a new "species" that already includes the unit conversion; then
# average this new "species"
#           echo "let spe_ax $species[$k,k=2:$maxheight@ave]axzx " >> pgtmp0
           echo "let spe_mr $species[$k]dd1 /$unit1a/$unit1c" >> pgtmp0
           echo "let spe_ax spe_mr[k2m1 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
           echo $plcmd $lstym[1] $lim1'`'"1.05*vmax"'`'  $xax spe_ax >> pgtmp0
       endif
    endif
# 2 runs
    if ($n_runs == 2) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       if ($pl_ty != 5) then
           echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[1] $lstym[2] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
       else
# run 1 prep
           echo "let spe_mr1 $species[$k]dd1 /$unit1a/$unit1c" >> pgtmp0
           echo "let spe_ax1 spe_mr1[k2m1 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax1[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# run 2 prep
           echo "let spe_mr3 $species[$k]dd3 /$unit3a/$unit3c" >> pgtmp0
           echo "let spe_ax3 spe_mr3[k2m3 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax3[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# plot
           echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
           echo $plcmd $lstym[1] $lim1'`'"1.05*vmax"'`'  $xax spe_ax1 >> pgtmp0
           echo $plcmdo $lstym[2] $xax spe_ax3 >> pgtmp0
       endif
    endif
# 3 runs
    if ($n_runs == 3) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de5,$pll1$pld[$kk]$lim2]/$unit5a/$unit5b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       if ($pl_ty != 5) then
           echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[1] $lstym[2] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[1] $lstym[3] $xax$species[$k]dd5/$unit5a/$unit5c$yax >> pgtmp0
       else
# run 1 prep
           echo "let spe_mr1 $species[$k]dd1 /$unit1a/$unit1c" >> pgtmp0
           echo "let spe_ax1 spe_mr1[k2m1 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax1[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# run 2 prep
           echo "let spe_mr3 $species[$k]dd3 /$unit3a/$unit3c" >> pgtmp0
           echo "let spe_ax3 spe_mr3[k2m3 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax3[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# run 3 prep
           echo "let spe_mr5 $species[$k]dd5 /$unit5a/$unit5c" >> pgtmp0
           echo "let spe_ax5 spe_mr5[k2m5 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax5[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# plot
           echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
           echo $plcmd $lstym[1] $lim1'`'"1.05*vmax"'`'  $xax spe_ax1 >> pgtmp0
           echo $plcmdo $lstym[2] $xax spe_ax3 >> pgtmp0
           echo $plcmdo $lstym[3] $xax spe_ax5 >> pgtmp0
       endif
    endif
# 4 runs
    if ($n_runs == 4) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de5,$pll1$pld[$kk]$lim2]/$unit5a/$unit5b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de7,$pll1$pld[$kk]$lim2]/$unit7a/$unit7b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       if ($pl_ty != 5) then
           echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[1] $lstym[2] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[1] $lstym[3] $xax$species[$k]dd5/$unit5a/$unit5c$yax >> pgtmp0
           echo $plcmdo$pld[1] $lstym[4] $xax$species[$k]dd7/$unit7a/$unit7c$yax >> pgtmp0
       else
# run 1 prep
           echo "let spe_mr1 $species[$k]dd1 /$unit1a/$unit1c" >> pgtmp0
           echo "let spe_ax1 spe_mr1[k2m1 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax1[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# run 2 prep
           echo "let spe_mr3 $species[$k]dd3 /$unit3a/$unit3c" >> pgtmp0
           echo "let spe_ax3 spe_mr3[k2m3 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax3[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# run 3 prep
           echo "let spe_mr5 $species[$k]dd5 /$unit5a/$unit5c" >> pgtmp0
           echo "let spe_ax5 spe_mr5[k2m5 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax5[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# run 4 prep
           echo "let spe_mr7 $species[$k]dd7 /$unit7a/$unit7c" >> pgtmp0
           echo "let spe_ax7 spe_mr7[k2m7 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax7[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# plot
           echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
           echo $plcmd $lstym[1] $lim1'`'"1.05*vmax"'`'  $xax spe_ax1 >> pgtmp0
           echo $plcmdo $lstym[2] $xax spe_ax3 >> pgtmp0
           echo $plcmdo $lstym[3] $xax spe_ax5 >> pgtmp0
           echo $plcmdo $lstym[4] $xax spe_ax7 >> pgtmp0
       endif
    endif
# 5 runs
    if ($n_runs == 5) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de5,$pll1$pld[$kk]$lim2]/$unit5a/$unit5b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de7,$pll1$pld[$kk]$lim2]/$unit7a/$unit7b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de9,$pll1$pld[$kk]$lim2]/$unit9a/$unit9b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       if ($pl_ty != 5) then
           echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[1] $lstym[2] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[1] $lstym[3] $xax$species[$k]dd5/$unit5a/$unit5c$yax >> pgtmp0
           echo $plcmdo$pld[1] $lstym[4] $xax$species[$k]dd7/$unit7a/$unit7c$yax >> pgtmp0
           echo $plcmdo$pld[1] $lstym[5] $xax$species[$k]dd9/$unit9a/$unit9c$yax >> pgtmp0
       else
# run 1 prep
           echo "let spe_mr1 $species[$k]dd1 /$unit1a/$unit1c" >> pgtmp0
           echo "let spe_ax1 spe_mr1[k2m1 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax1[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# run 2 prep
           echo "let spe_mr3 $species[$k]dd3 /$unit3a/$unit3c" >> pgtmp0
           echo "let spe_ax3 spe_mr3[k2m3 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax3[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# run 3 prep
           echo "let spe_mr5 $species[$k]dd5 /$unit5a/$unit5c" >> pgtmp0
           echo "let spe_ax5 spe_mr5[k2m5 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax5[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# run 4 prep
           echo "let spe_mr7 $species[$k]dd7 /$unit7a/$unit7c" >> pgtmp0
           echo "let spe_ax7 spe_mr7[k2m7 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax7[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# run 5 prep
           echo "let spe_mr9 $species[$k]dd9 /$unit9a/$unit9c" >> pgtmp0
           echo "let spe_ax9 spe_mr9[k2m9 $maxheight av] " >> pgtmp0
           echo "let vmax1=spe_ax9[l=@max]" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
# plot
           echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
           echo $plcmd $lstym[1] $lim1'`'"1.05*vmax"'`'  $xax spe_ax1 >> pgtmp0
           echo $plcmdo $lstym[2] $xax spe_ax3 >> pgtmp0
           echo $plcmdo $lstym[3] $xax spe_ax5 >> pgtmp0
           echo $plcmdo $lstym[4] $xax spe_ax7 >> pgtmp0
           echo $plcmdo $lstym[5] $xax spe_ax9 >> pgtmp0
       endif
    endif

#   add name of species
    echo "let label_y = 1.15*vmax" >> pgtmp0
    echo "if "'`'"$pl_ty eq 2 "'`'" then let label_y="'`'"1.05*eta[de2,l=1,k=$maxheight]"'`'" endif" >> pgtmp0
    echo "if "'`'"$pl_ty eq 3 "'`'" then let label_y="'`'"1.1*eta[de2,l=1,k=$maxheight]"'`'" endif" >> pgtmp0
#    echo "if "'`'"$pl_ty eq 3 "'`'" then let label_y=1.05*$maxheight endif" >> pgtmp0
#    echo "if "'`'"$pl_ty eq 3 "'`'" then let label_y=1.05*4150 endif" >> pgtmp0
    echo "label 0.," '`'"label_y"'`'",0,0,.2 @P1$species[$k]" >> pgtmp0
#   heights/times to be overplotted
    foreach kk ($m_pl)
        if ($n_runs == 1 && $kk != 1) echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       if ($n_runs == 2 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
       endif
       if ($n_runs == 3 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd5/$unit5a/$unit5c$yax >> pgtmp0
       endif
       if ($n_runs == 4 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd5/$unit5a/$unit5c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd7/$unit7a/$unit7c$yax >> pgtmp0
       endif
       if ($n_runs == 5 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd5/$unit5a/$unit5c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd7/$unit7a/$unit7c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd9/$unit9a/$unit9c$yax >> pgtmp0
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
        echo "label 10,110,  -1, 0, .25 @P1Gas phase" >> pgtmp0
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
sed 's/k2m1/d=1,k=2:/g' pgtmp0 > pgtmp0a
sed 's/k2m3/d=3,k=2:/g' pgtmp0a > pgtmp0b
sed 's/k2m5/d=5,k=2:/g' pgtmp0b > pgtmp0c
sed 's/k2m7/d=7,k=2:/g' pgtmp0c > pgtmp0d
sed 's/k2m9/d=9,k=2:/g' pgtmp0d > pgtmp0e
sed 's/d2k2/[d=2,k=2]/g' pgtmp0e > pgtmp0f
sed 's/av/@ave/g' pgtmp0f > pgtmp0g
sed 's/axzx2/rho[d=2,gz=zaxis\@asn,gt=taxis\@asn]/g' pgtmp0g > pgtmp0h
sed 's/axzx/[d=1,gz=zaxis\@asn,gt=taxis\@asn]/g' pgtmp0h > pgtmp0i
sed 's/dd10/[d=10]/g' pgtmp0i > pgtmp0j
sed 's/dd1/[d=1]/g' pgtmp0j > pgtmp1
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
if ($countPG == 14) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt"
if ($countPG == 15) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt"
if ($countPG == 16) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt"
if ($countPG == 17) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt"
if ($countPG == 18) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt"
if ($countPG == 19) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt"
if ($countPG == 20) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt $pl_name.20.plt"



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
