! ======================================================================= !
! Copyright 1996-2017 the Authors
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
! ======================================================================= !


! outp.f90: output
! box: the initial (const*) and overview (prof*) output SRs are not adjusted,
!      i.e. they produce a heap of output that's irrelevant. To avoid huge
!     output files, the size of the plou* SRs has been adjusted by using
!     n_bl and n_bln. Output starts at k=1 to save the deposited/surface
!     values as well.

! This file contains the following subroutines:
!
!     - outm: restart files, meteorological data
!     - outc: restart files, chemical data
!
! plout* routines are called if binout=true
!     - ploutj: output of photolysis rates up to level n_bln
!     - ploutc: driving subroutine of chemistry output, calls the following subroutines:
!              - ploutcg: output of gas phase chemical species
!              - ploutci: output of liquid phase ionic species (sion1)
!              - ploutcl: output of liquid phase species sl1
!              - ploutcr:
!              - ploutcgr: output of instantaneous rates
!              - ploutcgs:
!     - ploutm: output of meteorological data
!     - ploutp: output of particle distributions
!     - ploutr: output of radiation data
!     - ploutt: output of turbulence data
!
!     - constm: output of constants and parameters
!     - constc: output of chemical constants and parameters
!
!     - profm: output of meteorological profiles
!     - profc: vertical profiles of chemical species
!     - profr: output of radiation variables



subroutine outm
!
! Description:
!    profiles of meteorological data needed to restart the program
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
!          05/2021   Bugfix: add z0 in restart file           <Josue Bock>
! 1.2      08/2016   Use module for parameters                <Josue Bock>
!                    Comments / header
!                    Removed one unused argument
!
! 1.1       ?        Output every 12h with increasing names   <Roland von Glasow>
!                    Changed most common blocks
!
! 1.0       ?        Original code.                           <Andreas Bott>
!

!
! Declarations:
! Modules used:
  USE data_surface, ONLY : &
       tw, &                       ! water surface temperature
       ustern, z0,&                ! frictional velocity, roughness length
       gclu, gclt                  ! coefficients for momentum, and temperature and humidity

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunrstm

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nb, &
       nka, &
       nkt, &
       mb

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  character (len=10) :: fname
  character (len=1) :: sub

! Common blocks:
  common /cb11/ totrad (mb,n)
  real (kind=dp) :: totrad
  common /cb18/ alat,declin                ! for the SZA calculation
  real (kind=dp) :: alat,declin
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct
  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy
  common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
  real (kind=dp) :: gm, gh, sm, sh, xl
  common /cb44/ a0m,b0m(nka)
  real (kind=dp) :: a0m,b0m
  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w
  common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb), &
                ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
  real (kind=dp) :: zb, dzb, dzbw, tb, eb, ak, d, &
       ajb, ajq, ajl, ajt, ajd, ajs, ds1, ds2, ajm, reif, tau, trdep
  common /cb48/ sk,sl,dtrad(n),dtcon(n)
  real (kind=dp) :: sk, sl, dtrad, dtcon
  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb53a/ thet(n),theti(n)
  real(kind=dp) :: thet, theti
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a
  common /cb63/ fcs(nka),xmol3(nka)
  real(kind=dp) :: fcs, xmol3

! == End of declarations =======================================================

  fname='rstm .dat'
  sub='_'

  if (lday*24+lst.eq.12) sub='A'
  if (lday*24+lst.eq.24) sub='B'
  if (lday*24+lst.eq.36) sub='C'
  if (lday*24+lst.eq.48) sub='D'
  if (lday*24+lst.eq.60) sub='E'
  if (lday*24+lst.eq.72) sub='F'
  if (lday*24+lst.eq.84) sub='G'
  if (lday*24+lst.eq.96) sub='H'
  if (lday*24+lst.eq.108) sub='I'
  if (lday*24+lst.eq.120) sub='J'
  if (lday*24+lst.eq.132) sub='K'
  if (lday*24+lst.eq.144) sub='L'
  if (lday*24+lst.eq.156) sub='M'
  if (lday*24+lst.eq.168) sub='N'
  if (lday*24+lst.eq.180) sub='O'
  if (lday*24+lst.eq.192) sub='P'
  fname(5:5)=sub

 3000 continue
  open (jpfunrstm,file=fname,status='unknown',form='unformatted',err=3000)
  write (jpfunrstm) &
! double precision arrays
       atkm,atkh,b0m,dfddt,dtrad,eb,ff,fcs,feu,fsum, &
       gh,p,rho,t,talt,tb,thet,theta,theti,thetl,tke,tkep, &
       totrad,u,v,w,xl,xm1,xm1a,xm2,xmol3, &
! double precision single vars
       a0m,alat,declin,ds1,ds2,gclt,gclu,reif,sk,sl,tau, &
       trdep,tw,ustern,z0, &
! integer arrays
       nar, &
! integer single vars
       it,lcl,lct,lday,lmin,lst
  close (jpfunrstm)

end subroutine outm
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine outc
!
! Description:
!    output of chemical data needed to restart the program
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      08/2016   Use module for parameters                <Josue Bock>
!                    Comments / header
!                    Removed one unused argument
!                    Removed unnecessary and outdated output
!
! 1.1       ?        Output every 12h with increasing names   <Roland von Glasow>
!                    Changed most common blocks
!
! 1.0       ?        Original code.                           <Andreas Bott>
!

!
! Declarations:
! Modules used:
  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunrstc

  USE gas_common, ONLY: &
! Imported Array Variables with intent (in):
       s1, &
       es1, &
       s3

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j6, &
       nf, &
       n, &
       nka, &
       nkt, &
       nkc, &
       nphrxn, &
       nlev, &
       nrxn

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  character (len=10) :: fname
  character (len=1) :: sub

! Common blocks:
  common /band_rat/ photol_j(nphrxn,n)
  real (kind=dp) :: photol_j
  common /blck01/ am3(n),cm3(n)
  real (kind=dp) :: am3, cm3
  common /blck11/ rc(nkc,n)
  real (kind=dp) :: rc
  common /blck12/ cw(nkc,n),cm(nkc,n)
  real (kind=dp) :: cw, cm
  common /blck13/ conv2(nkc,n) ! conversion factor = 1/(1000*cw)
  real (kind=dp) :: conv2
  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1
  common /blck78/ sa1(j2,nka)
  real (kind=dp) :: sa1
  common /budg/ bg(2,nrxn,nlev),il(nlev)
  real (kind=dp) :: bg
  integer :: il
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct
  common /kinv_i/ kinv
  integer :: kinv
  common /kpp_crys/ xcryssulf,xcrysss,xdelisulf,xdeliss
  real (kind=dp) :: xcryssulf,xcrysss,xdelisulf,xdeliss
  common /kpp_l1/ cloudt(nkc,n)
  logical :: cloudt
  common /kpp_mol/ xgamma(j6,nkc,nf)
  real (kind=dp) :: xgamma
  common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
  real (kind=dp) :: vt, vd, vdm

! == End of declarations =======================================================

  fname='rstc .dat'
  sub='_'
  if (lday*24+lst.eq.12) sub='A'
  if (lday*24+lst.eq.24) sub='B'
  if (lday*24+lst.eq.36) sub='C'
  if (lday*24+lst.eq.48) sub='D'
  if (lday*24+lst.eq.60) sub='E'
  if (lday*24+lst.eq.72) sub='F'
  if (lday*24+lst.eq.84) sub='G'
  if (lday*24+lst.eq.96) sub='H'
  if (lday*24+lst.eq.108) sub='I'
  if (lday*24+lst.eq.120) sub='J'
  if (lday*24+lst.eq.132) sub='K'
  if (lday*24+lst.eq.144) sub='L'
  if (lday*24+lst.eq.156) sub='M'
  if (lday*24+lst.eq.168) sub='N'
  if (lday*24+lst.eq.180) sub='O'
  if (lday*24+lst.eq.192) sub='P'
  fname(5:5)=sub

3000 continue
  open (jpfunrstc,file=fname,status='unknown',form='unformatted',err=3000)
  write (jpfunrstc) &
! double precision arrays
       am3,cm,cm3,conv2,cw,es1,photol_j,rc,s1,s3,sa1, &
       sl1,sion1,vd,vdm,vt,xgamma, &
! double precision, single values
       xcryssulf,xcrysss,xdelisulf,xdeliss, &
! logicals
       cloudt, &
! integers
       il,kinv,lday,lmin,lst
  close (jpfunrstc)

end subroutine outc
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ploutj (fogtype,n_bln)
!
! Description:
!    output of photolysis rates
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      08/2016   Use module for parameters                <Josue Bock>
!                    Comments / header
!                    Removed one unused argument
!                    Cleaning
!
! 1.0       ?        Original code.                          <Roland von Glasow>
!

!
! Declarations:
! Modules used:
  USE config, ONLY : &
       coutdir

  USE file_unit, ONLY : &
       jpfunjra

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nphrxn

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype
  integer, intent(in) :: n_bln

! Local scalars:
  character (len=10) :: fname
  character (len=150) :: clpath  ! complete path to file
  integer :: j,k
  real (kind=dp) :: xday,xst,xmin

! Local arrays:
  real (kind=dp) :: i0(nphrxn,n_bln)   ! the local array where are written the photolysis rates, up to the selected level (1:n_bln)

! Common blocks:
  common /band_rat/ photol_j(nphrxn,n)
  real (kind=dp) :: photol_j

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

! == End of declarations =======================================================

  do k=1,n_bln
     do j=1,nphrxn
        i0(j,k)=photol_j(j,k)
     enddo
  enddo

  xday = real(lday,dp)
  xst  = real(lst,dp)
  xmin = real(lmin,dp)

  fname='jra .out'
  fname(4:4)=fogtype
  clpath=trim(coutdir)//trim(fname)

  open (jpfunjra, file=trim(clpath), status='old',form='unformatted', position='append')
!  write (jpfunjra) lday,lst,lmin,i0
  write (jpfunjra) xday,xst,xmin,i0
  close (jpfunjra)

end subroutine ploutj
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ploutc (fogtype,mic,n_bl,n_bl8)

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype
  logical, intent(in) :: mic
  integer, intent(in) :: n_bl,n_bl8

! == End of declarations =======================================================

  call ploutcg (fogtype,n_bl)
  if (mic) then
     call ploutci (fogtype,n_bl)
     call ploutcl (fogtype,n_bl)
  end if
  call ploutcr (fogtype,n_bl)
  call ploucgr (fogtype,n_bl8)
  call ploucgs (fogtype,n_bl)

end subroutine ploutc
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ploutcg (fogtype,n_bl)
!
! Description:
!    output of gas phase chemical species s1:
!
!    1. NO      2. NO2     3. HNO3     4. NH3     5. SO2      6. H2SO4
!    7. O3      8. CH4     9. C2H6    10. C3H8   11. ALKA    12. ETHE
!   13. ALKE   14. AROM   15. HCOOH   16. ACTA   17. HCHO    18. ALD2
!   19. H2O2   20. CH3OOH 21. HONO    22. PAN    23. TPAN    24. KET
!   25. CRES   26. DIAL   27. GLYX    28. MGLY   29. NH4NO3  30. HCl
!   31. R3N2   32. RAN2   33. RAN1    34. N2O5   35. HNO4    36. NO3
!   37. DMS    38. HOCl   39. ClNO2   40. ClNO3  41. Cl2     42. HBr
!   43. HOBr   44. BrNO2  45. BrNO3   46. Br2    47. BrCl    48. HI
!   49. HOI    50. I2O2   51. INO2    52. INO3   53. I2      54. ICl
!   55. IBr    56. CH3I   57. CH2I2   58. CH2ClI 59. C3H7I   60. DMSO
!   61. CH3SO2 62. CH3SO3 63. CH3SO3H 64. CO     65. Cl2O2   66. DMOO
!   67. CH3S   68. CH3SO  69. MSIA    70. DMSO2  71. CH2BrI  72. CHBr2I
!   73. C2H5I  74. HIO3   75. NUCV    76. SO3    77. HOSO2   78. CO2
!   79. I2O    80. I2O3   81. I2O4    82. I2O5   83. INO     84. Br2O
!   85. ClONO  86. ClO3   87. Cl2O3   88. CH3OH  89. C2H5OH  90. H2
!   91. NHS    92. RCl    93. RBr     94. XOR    95. SOR     96. SPAN
!   97. Hg     98. HgO    99. HgCl   100. HgCl2 101. HgBr   102. HgBr2
!   102 -- 130 undefined
!
! DMOO = CH3SCH2OO, MSIA = CH3S(O)OH
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      08/2016   Use module          <Josue Bock>
!                    Comments / header
!
! 1.0       ?        Original code.      <Roland von Glasow>
!

!
! Declarations:
! Modules used:
  USE config, ONLY : &
       coutdir

  USE file_unit, ONLY : &
       jpfunsg1

  USE gas_common, ONLY : &
! Imported Parameters:
       j1, &
! Imported Array Variables with intent (in):
       s1

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype
  integer, intent(in) :: n_bl

! Local scalars:
  character (len=10) :: fname
  character (len=150) :: clpath  ! complete path to file
  integer :: j,k, nmax

! Local arrays:
  real (kind=dp) :: i0(n_bl,j1)

! Common blocks:
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

! == End of declarations =======================================================

  nmax = n_bl
  if (n_bl.gt.3) nmax = n_bl-2
  do k=1,nmax
     do j=1,j1
        i0(k,j)=s1(j,k)
     enddo
  enddo
! only for 1D:
  if (n_bl.gt.3) then
     do j=1,j1
        i0(nf-1,j)=etw(lcl)
        i0(nf,j)=etw(lct)
     enddo
  endif

  fname='sg1 .out'
  fname(4:4)=fogtype
  clpath=trim(coutdir)//trim(fname)

  open (jpfunsg1, file=trim(clpath), status='old',form='unformatted', position='append')
  write (jpfunsg1) lday,lst,lmin,i0
  close (jpfunsg1)

end subroutine ploutcg
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ploutci (fogtype,n_bl)
!
! Description:
!    output of ions sion1:
!
!    1. H+       2. NH4+     3. OH-      4. CH2OHSO3-  5. HSO3-
!    6. SO3=     7. SO4-     8. SO4=     9. HCO3-     10. CO3-
!   11. O2-     12. NO2-    13. NO3-    14. Cl-       15. Cl2-
!   16. HCOO-   17. Fe3+    18. Mn2+    19. HSO4-     20. Na+ (check electroneg)
!   21. NO4-    22. ClO-    23. ClOH-   24. Br-       25. Br2-
!   26. BrO-    27. BrOH-   28. BrCl2-  29. Br2Cl-    30. CH3SO3-
!   31. HSO5-   32. SO3-    33. SO5-    34. I-        35. IO2-
!   36. IO3-    37. ICl2-   38. IBr2-   39. MS-       40. Hg+
!   41. Hg2+    42. HgOH+   43. HgCl+   44. HgCl3-    45. HgCl42-
!   46. HgBr+   47. HgBr3-  48. HgBr42- 49. Hg(SO3)22- 50. --
!
! MS- = CH3S(O)O-, CH3SO3- = CH3S(OO)O-
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      08/2016   Use module          <Josue Bock>
!                    Comments / header
!
! 1.0       ?        Original code.      <Roland von Glasow>
!

!
! Declarations:
! Modules used:
  USE config, ONLY: &
       coutdir, &
       nkc_l

  USE file_unit, ONLY : &
       jpfunion

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j6, &
       nkc, &
       nf, &
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype
  integer, intent(in) :: n_bl

  ! Local scalars:
  integer :: i, j, k, nmax
  character (len=10) :: fname
  character (len=150) :: clpath  ! complete path to file

! Local arrays:
  real (kind=dp) :: i0(n_bl,j6,nkc_l)

! Common blocks:
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1

! == End of declarations =======================================================

  nmax = n_bl
  if (n_bl.gt.3) nmax = n_bl-2
  do k=1,nmax
     do j=1,j6
        do i=1,nkc_l
           i0(k,j,i)=sion1(j,i,k)
        enddo
     enddo
  enddo
! only for 1D:
  if (n_bl.gt.3) then
     do j=1,j6
        do i=1,nkc_l
           i0(nf-1,j,i)=etw(lcl)
           i0(nf,j,i)=etw(lct)
        enddo
     enddo
  endif

  fname='ion .out'
  fname(4:4)=fogtype
  clpath=trim(coutdir)//trim(fname)

  open (jpfunion, file=trim(clpath), status='old',form='unformatted', position='append')
  write (jpfunion) lday,lst,lmin,i0
  close (jpfunion)

end subroutine ploutci
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ploutcl (fogtype,n_bl)
!
! Description:
!    output of liquid phase species sl1:
!
!    1. NO      2. NO2     3. HNO3     4. NH3     5. SO2      6. H2SO4
!    7. O3      8. CH4     9. C2H6    10. C3H8   11. ALKA    12. ETHE
!   13. ALKE   14. AROM   15. HCOOH   16. ACTA   17. HCHO    18. ALD2
!   19. H2O2   20. CH3OOH 21. HONO    22. PAN    23. TPAN    24. KET
!   25. CRES   26. DIAL   27. GLYX    28. MGLY   29. NH4NO3  30. HCl
!   31. R3N2   32. RAN2   33. RAN1    34. N2O5   35. HNO4    36. NO3
!   37. DMS    38. HOCl   39. ClNO2   40. ClNO3  41. Cl2     42. HBr
!   43. HOBr   44. BrNO2  45. BrNO3   46. Br2    47. BrCl    48. HI
!   49. HOI    50. I2O2   51. INO2    52. INO3   53. I2      54. ICl
!   55. IBr    56. CH3I   57. CH2I2   58. CH2ClI 59. C3H7I   60. DMSO
!   61. CH3SO2 62. CH3SO3 63. CH3SO3H 64. CO     65. Cl2O2   66. DMOO
!   67. CH3S   68. CH3SO  69. MSIA    70. DMSO2  71. CH2BrI  72. CHBr2I
!   73. C2H5I  74. HIO3   75. NUCV    76. SO3    77. HOSO2   78. CO2
!   79. I2O    80. I2O3   81. I2O4    82. I2O5   83. INO     84. Br2O
!   85. ClONO  86. ClO3   87. Cl2O3   88. CH3OH  89. C2H5OH  90. H2
!   91. NHS    92. RCl    93. RBr     94. XOR    95. SOR     96. SPAN
!   97. Hg     98. HgO    99. HgCl   100. HgCl2 101. HgBr   102. HgBr2
!
! DMOO = CH3SCH2OO, MSIA = CH3S(O)OH.
!
!   j2-j3+1. --     j2-j3+2. OH    j2-j3+3. HO2    j2-j3+4. DOM     j2-j3+5.  HIO2
!   j2-j3+6. CH3OO  j2-j3+7. IO    j2-j3+8. Cl     j2-j3+9. Br      j2-j3+10. --
!   j2-j3+11.O2     j2-j3+12.OIO   j2-j3+13. HgOH2 j2-j3+14. HgOHCl j2-j3+15. HgSO3
!   j2-j3+16. HgOHBr j2-j3+17. --  j2-j3+18. --    j2-j3+19. --     j2-j3+20. --
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      08/2016   Use module          <Josue Bock>
!                    Comments / header
!
! 1.0       ?        Original code.      <Roland von Glasow>
!

!
! Declarations:
! Modules used:
  USE config, ONLY: &
       coutdir, &
       nkc_l

  USE file_unit, ONLY : &
       jpfunsl1

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j6, &
       nkc, &
       nf, &
       n, &
       nka, &
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp


  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype
  integer, intent(in) :: n_bl

! Local scalars:
  character (len=10) :: fname
  character (len=150) :: clpath  ! complete path to file
  integer :: i, ia, j, jt, k, kc, nmax
  real (kind=dp) :: x0

! Local arrays:
  real (kind=dp) :: i0(n_bl,j2,nkc_l),irc(n_bl,nkc_l),icw(n_bl,nkc_l)

! Common blocks:
  common /blck11/ rc(nkc,n)
  real (kind=dp) :: rc
  common /blck12/ cw(nkc,n),cm(nkc,n)
  real (kind=dp) :: cw, cm
  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

! == End of declarations =======================================================

  nmax = n_bl
  if (n_bl.gt.3) nmax = n_bl-2
  do k=1,nmax
     do j=1,j2
        do i=1,nkc_l
           i0(k,j,i)=sl1(j,i,k)
        enddo
     enddo
     x0=0._dp
     do ia=1,nka
        do jt=1,nkt
           x0=x0+ff(jt,ia,k)*en(ia)
        enddo
     enddo
     i0(k,j2,1)=x0
     i0(k,j2,2)=fsum(k)
  enddo
! only for 1D:
  if (n_bl.gt.3) then
     do j=1,j2
        do i=1,nkc_l
           i0(nf-1,j,i)=etw(lcl)
           i0(nf,j,i)=etw(lct)
        enddo
     enddo
  endif

  fname='sl1 .out'
  fname(4:4)=fogtype
  clpath=trim(coutdir)//trim(fname)

  do kc=1,nkc_l
     do k=1,n_bl
        irc(k,kc)=rc(kc,k) ! jjb this transposition is maybe useless, just to stick to old write
        icw(k,kc)=cw(kc,k)
     enddo
  enddo

  open (jpfunsl1, file=trim(clpath), status='old',form='unformatted', position='append')
  write (jpfunsl1) lday,lst,lmin,i0,irc,icw
  close (jpfunsl1)

end subroutine ploutcl
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ploutcr (fogtype,n_bl)


! output of radical species s3:
! species 1-3 are treated as long lived
!    1. ----    2. ----    3. ---     4. OH      5. HO2
!    6. AHO2    7. MCO3    8. CH3OO   9. ETO2   10. KO2
!   11. R3O2   12. RAO2   13. TO2    14. TCO3   15. ZO2
!   16. EO2    17. PO2    18. CHO2   19. CRO2   20. PRN1
!   21. O(1D)  22. Cl     23. ClO    24. OClO   25. Br
!   26. BrO    27. I      28. IO     29. OIO    30. O3P
!   31. ClRO2  32. BrRO2  33. IRO2

! Modifications:
!     30-11-2016   J. Bock   implicit none, header
!     30-11-2016   J. Bock   correction in 1D case: indexes n_bl / n_bl-1
!                            instead of nf / nf-1
!     30-11-2016   J. Bock   loop order: priority for write (efficiency)
!     14-02-2017   J. Bock   renamed MO2 -> CH3OO
!                            MO2 was only used in the gas phase, while CH3OOlz was used in aq. phase,
!                            thus it could not be recognised by the newly developed routine searching
!                            for exchanged species.

  USE config, ONLY : &
       coutdir

  USE file_unit, ONLY : &
       jpfunsr1

  USE gas_common, ONLY: &
! Imported Parameters:
       j5, &
! Imported Array Variables with intent (in):
       s3

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype
  integer, intent(in) :: n_bl

! Local scalars:
  character (len=10) :: fname
  character (len=150) :: clpath  ! complete path to file
  integer :: nmax, j, k

! Local arrays:
  real (kind=dp) :: i0(n_bl,j5)

! Common blocks:
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

! == End of declarations =======================================================

  nmax = n_bl
  if (n_bl.gt.3) nmax = n_bl-2

  do j=1,j5
     do k=1,nmax
        i0(k,j)=s3(j,k)
     enddo
  enddo
! only for 1D:
  if (n_bl.gt.3) then
     do j=1,j5
        i0(n_bl-1,j)=etw(lcl)
        i0(n_bl,j)=etw(lct)
     enddo
  endif

  fname='sr1 .out'
  fname(4:4)=fogtype
  clpath=trim(coutdir)//trim(fname)

  open (jpfunsr1, file=trim(clpath), status='old',form='unformatted', position='append')
  write (jpfunsr1) lday,lst,lmin,i0
  close (jpfunsr1)

end subroutine ploutcr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ploucgr (fogtype,n_bl8)
!
! Description:
!    output of instantaneous rates [mol/(m3 s)]

! Modifications:
!     30-11-2016   J. Bock   USE module instead of local parameters
!     30-11-2016   J. Bock   implicit none

! Declarations:
! Modules used:
  USE config, ONLY : &
       coutdir

  USE file_unit, ONLY : &
       jpfungr

  USE global_params, ONLY : &
! Imported Parameters:
       nlev, &
       nrxn

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype
  integer, intent(in) :: n_bl8

! Local scalars:
  character (len=10) :: fname
  character (len=150) :: clpath  ! complete path to file
  integer :: j, k

! Local arrays:
  real (kind=dp) :: i0(n_bl8,nrxn)  ! the local array where are written the photolysis
                                    !   rates, up to the selected level (1:n_bl8)


! Common blocks:
  common /budg/ bg(2,nrxn,nlev),il(nlev)
  real (kind=dp) :: bg
  integer :: il

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

! == End of declarations =======================================================

  do j=1,nrxn
     do k=1,n_bl8
        i0(k,j)=bg(1,j,k)
     enddo
  enddo

  fname='gr .out'
  fname(3:3)=fogtype
  clpath=trim(coutdir)//trim(fname)

  open (jpfungr, file=trim(clpath), status='old',form='unformatted', position='append')
  write (jpfungr) lday,lst,lmin,il,i0
  close (jpfungr)

end subroutine ploucgr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ploucgs (fogtype,n_bl)

!
! Description:
!     output of gas phase rates

! Modifications:
!     30-11-2016   J. Bock   USE module instead of local parameters
!     30-11-2016   J. Bock   implicit none, header


! Declarations:
! Modules used:
  USE config, ONLY : &
       coutdir

  USE file_unit, ONLY : &
       jpfungs

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype
  integer, intent(in) :: n_bl

! Local scalars:
  character (len=10) :: fname
  character (len=150) :: clpath  ! complete path to file
  integer :: j, k

! Local arrays:
  real (kind=dp) :: i0(n_bl,122),i1(n_bl,122)

! Common blocks:
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /budgs/ bgs(2,122,n)
  real (kind=dp) :: bgs

! == End of declarations =======================================================

  do j=1,122
     do k=1,n_bl
        i0(k,j)=bgs(1,j,k)
        i1(k,j)=bgs(2,j,k)
     enddo
  enddo

  fname='gs .out'
  fname(3:3)=fogtype
  clpath=trim(coutdir)//trim(fname)

  open (jpfungs, file=trim(clpath), status='old',form='unformatted', position='append')
  write (jpfungs) lday,lst,lmin,i0,i1
  close (jpfungs)

end subroutine ploucgs
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ploutm (fogtype,n_bln)
! output of meteorological data

  USE config, ONLY : &
       coutdir

  USE file_unit, ONLY : &
       jpfunpm, jpfunpb        ! ploutm files: pm*, pb*

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nrlay, &
       nb, &
       nka, &
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype
  integer, intent(in) :: n_bln

! Local scalars:
  character (len=10) :: fname
  character (len=150) :: clpath  ! complete path to file
  integer :: j, k
! Local arrays:
  real (kind=dp) :: i0(12,n_bln)

! Common blocks:
  common /cb15/ fnseb,flgeg,hr(nrlay)
  real (kind=dp) :: fnseb, flgeg, hr

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy
  common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb), &
                ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
  real (kind=dp) :: zb, dzb, dzbw, tb, eb, ak, d, &
       ajb, ajq, ajl, ajt, ajd, ajs, ds1, ds2, ajm, reif, tau, trdep
  common /cb48/ sk,sl,dtrad(n),dtcon(n)
  real (kind=dp) :: sk, sl, dtrad, dtcon

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a

! == End of declarations =======================================================

  do k=1,n_bln
     i0(1,k)=rho(k)
     i0(2,k)=atkh(k)
     i0(3,k)=theta(k)
     i0(4,k)=t(k)
     i0(5,k)=thetl(k)
     i0(6,k)=dtcon(k)
     i0(7,k)=dtrad(k)
     i0(8,k)=fsum(k)
     i0(9,k)=xm1(k)
     i0(10,k)=xm2(k)
     i0(11,k)=feu(k)
  enddo
! only for 1D:
  if (n_bln.gt.3) then
     do j=1,12
        i0(j,n-1)=etw(lcl)
        i0(j,n)=etw(lct)
     enddo
  endif

! ploutm file 1 = pm*
  fname='pm .out'
  fname(3:3)=fogtype
  clpath=trim(coutdir)//trim(fname)
  open (jpfunpm, file=trim(clpath), status='old',form='unformatted', position='append')
  write (jpfunpm) lday,lst,lmin,i0
  close (jpfunpm)

! ploutm file 2 = pb*
  fname(2:2)='b'
  clpath=trim(coutdir)//trim(fname)
  open (jpfunpb, file=trim(clpath), status='old',form='unformatted', position='append')
  write (jpfunpb) lday,lst,lmin,tb,eb,ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2, &
       ajm,reif,tau,trdep,sl,sk,fnseb,flgeg
  close (jpfunpb)

end subroutine ploutm
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ploutp (fogtype)
! output of particle distributions


! 22-Sep-2020   Josue Bock   Remove remaining hard coded nf,n,nka, nkt parameter values

  USE config, ONLY : &
       coutdir

  USE file_unit, ONLY : &
       jpfunf1, jpfunf2, jpfunf3  ! ploutp files: f1*, f2*, f3*

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nka, &
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype

! Local scalars:
  character (len=10) :: fname
  character (len=150) :: clpath  ! complete path to file
  integer :: kk, ia, jt, k     ! loop indexes
  real (kind=dp) :: x0, x1, x2
! Local arrays:
  real (kind=dp) :: ff2(nka,nkt)

! Common blocks:
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a

! == End of declarations =======================================================

  fname='f1 .out'
  fname(3:3)=fogtype
  clpath=trim(coutdir)//trim(fname)
  open (jpfunf1, file=trim(clpath), status='old',form='unformatted', position='append')
  do kk=1,3
!     if (kk.eq.1) k=lcl-5
!     if (kk.eq.2) k=lcl
!     if (kk.eq.3) k=(lct+lcl)/2
!      k=max0(k,5)
     if (kk.eq.1) k=2
     if (kk.eq.2) k=12
     if (kk.eq.3) k=22
     do ia=1,nka
        do jt=1,nkt
           ff2(ia,jt)=ff(jt,ia,k)
        enddo
     enddo
     x0=xm2(k)
     x1=feu(k)
     x2=eta(k)
     write (jpfunf1) lday,lst,lmin,k,x0,x1,x2,ff2
  enddo
  close (jpfunf1)

  fname(2:2)='2'
  clpath=trim(coutdir)//trim(fname)
3010 continue
  open (jpfunf2, file=trim(clpath), status='old',form='unformatted', position='append',err=3010)
  do kk=1,3
!     if (kk.eq.1) k=lct-4
!     if (kk.eq.2) k=lct-2
!     if (kk.eq.3) k=lct
!     k=max0(k,5)
     if (kk.eq.1) k=32
     if (kk.eq.2) k=42
     if (kk.eq.3) k=52
     do ia=1,nka
        do jt=1,nkt
           ff2(ia,jt)=ff(jt,ia,k)
        enddo
     enddo
     x0=xm2(k)
     x1=feu(k)
     x2=eta(k)
     write (jpfunf2) lday,lst,lmin,k,x0,x1,x2,ff2
  enddo
  close (jpfunf2)

  fname(2:2)='3'
  clpath=trim(coutdir)//trim(fname)
3020 continue
  open (jpfunf3, file=trim(clpath), status='old',form='unformatted', position='append',err=3020)
  do kk=1,3
!     if (kk.eq.1) k=lct+2
!     if (kk.eq.2) k=lct+4
!     if (kk.eq.3) k=lct+6
!     k=max0(k,5)
     if (kk.eq.1) k=62
     if (kk.eq.2) k=72
     if (kk.eq.3) k=82
     do ia=1,nka
        do jt=1,nkt
           ff2(ia,jt)=ff(jt,ia,k)
        enddo
     enddo
     x0=xm2(k)
     x1=feu(k)
     x2=eta(k)
     write (jpfunf3) lday,lst,lmin,k,x0,x1,x2,ff2
  enddo
  close (jpfunf3)

end subroutine ploutp
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ploutr (fogtype,n_bln)
! output of radiation data

  USE config, ONLY : &
       coutdir

  USE file_unit, ONLY : &
       jpfunpr                ! ploutr file: pr*

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nrlay, &
       nrlev, &
       nka, &
       nkt, &
       mb, &
       mbs

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype
  integer, intent(in) :: n_bln

! Local scalars:
  character (len=10) :: fname
  character (len=150) :: clpath  ! complete path to file
  integer :: ia, ib, j, jt, k
  real (kind=dp) :: rsum1, rsum3
  real (kind=dp) :: x0, x1
! Local arrays:
  real (kind=dp) :: i0(12,n_bln)

! Common blocks:
  common /cb10/ totrad (mb,nrlay)
  real (kind=dp) :: totrad

  common /cb15/ fnseb,flgeg,hr(nrlay)
  real (kind=dp) :: fnseb, flgeg, hr

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /kurz/ fs1(nrlev),fs2(nrlev),totds(nrlev),ss(nrlev), &
                fsn(nrlev),dtdts(nrlay)
  real (kind=dp) :: fs1, fs2, totds, ss, fsn, dtdts

  common /lang/ fl1(nrlev),fl2(nrlev),fln(nrlev),dtdtl(nrlay)
  real (kind=dp) :: fl1, fl2, fln, dtdtl

! == End of declarations =======================================================

  do k=1,n_bln
     i0(1,k)=fs1(k)
     i0(2,k)=fs2(k)
     i0(3,k)=totds(k)
     i0(4,k)=dtdts(k)
     i0(5,k)=fl1(k)
     i0(6,k)=fl2(k)
     i0(7,k)=dtdtl(k)
     do ib=1,mbs
        i0(8,k)=i0(8,k)+totrad(ib,k)
     enddo
     do ib=mbs+1,mb
        i0(9,k)=i0(9,k)+totrad(ib,k)
     enddo
     rsum1=0._dp
!     rsum2=0._dp
     rsum3=0._dp
     do ia=1,nka
        do jt=2,nkt
           x0=rq(jt,ia)
           x1=(rw(jt,ia)-rw(jt-1,ia))
           rsum1=rsum1+ff(jt,ia,k)*x0*x1
!           rsum2=rsum2+ff(jt,ia,k)*x0**2*x1
           rsum3=rsum3+ff(jt,ia,k)*x0**3*x1
        enddo
     enddo
     i0(10,k)=rsum1
!     i0(11,k)=rsum2
     i0(12,k)=rsum3
     i0(11,k)=totrad(1,k)
  enddo
! only for 1D:
  if (n_bln.gt.3) then
     do j=1,12
        i0(j,n-1)=etw(lcl)
        i0(j,n)=etw(lct)
     enddo
     i0(1,n-2)=fnseb
     i0(1,n-3)=flgeg
  endif

! ploutr file = pr*
  fname='pr .out'
  fname(3:3)=fogtype
  clpath=trim(coutdir)//trim(fname)
  open (jpfunpr, file=trim(clpath), status='old',form='unformatted', position='append')
  write (jpfunpr) lday,lst,lmin,i0
  close (jpfunpr)

end subroutine ploutr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ploutt (fogtype,n_bln)
! output of turbulence data


! 22-Sep-2020   Josue Bock   Remove remaining hard coded nf,n,nka, nkt parameter values

  USE config, ONLY : &
       coutdir

  USE file_unit, ONLY : &
       jpfunpt                 ! ploutt files: pt*

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype
  character (len=150) :: clpath  ! complete path to file
  integer, intent(in) :: n_bln

! Local scalars:
  character (len=10) :: fname
  integer :: j, k
! Local arrays:
  real (kind=dp) :: i0(12,n_bln)

! Common blocks:
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy

  common /cb42a/ tkeps(n),tkepb(n),tkepd(n)
  real (kind=dp) :: tkeps, tkepb, tkepd

  common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
  real (kind=dp) :: gm, gh, sm, sh, xl

  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w

! == End of declarations =======================================================

  do k=1,n_bln
     i0(1,k)=atkm(k)
     i0(2,k)=u(k)
     i0(3,k)=v(k)
     i0(4,k)=w(k)
     i0(5,k)=tke(k)
     i0(6,k)=tkeps(k)
     i0(7,k)=tkepb(k)
     i0(8,k)=tkepd(k)
     i0(9,k)=-atkh(k)*buoy(k)
     i0(10,k)=buoy(k)
     i0(11,k)=xl(k)
  enddo
! only for 1D:
  if (n_bln.gt.3) then
     do j=1,12
        i0(j,n-1)=etw(lcl)
        i0(j,n)=etw(lct)
     enddo
  endif
  fname='pt .out'
  fname(3:3)=fogtype
  clpath=trim(coutdir)//trim(fname)
  open (jpfunpt, file=trim(clpath), status='old',form='unformatted', position='append')
  write (jpfunpt) lday,lst,lmin,i0
  close (jpfunpt)

end subroutine ploutt
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine constm

! Description :
! -----------
!   output of constants and parameters used in the current run
!    (called only once during model initialisation)

! Author :
! ------
!    Andreas Bott

! Modifications :
! -------------
!    Josue Bock  Review, add header, missing declarations, f90
!    Josue Bock  Add a test to print soil data only if isurf==1 (21-May-2021)

! == End of header =============================================================

  USE config, ONLY : &
! Imported Parameters:
       chem, mic, rst, &
       ug, vg, wmin, wmax, &
       isurf, scaleo3_m

  USE constants, ONLY : &
! Imported Parameters:
       cp, &              ! Specific heat of dry air, in J/(kg.K)
       g,  &              ! Gravitational acceleration (m/s**2)
       r0, &              ! Specific gas constant of dry air, in J/(kg.K)
       r1                 ! Specific gas constant of water vapour, in J/(kg.K)

  USE data_surface, ONLY : &
       tw, &                       ! water surface temperature
       z0                          ! roughness length

  USE file_unit, ONLY : &
       jpfunprofm

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nb, &
       nka, &
       nkt, &
       nrlay, &
       nrlev

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  real (kind=dp) :: xxsum
  integer :: ia, jt, k
! Local arrays:
  real (kind=dp) :: xsum(n)

! Common blocks:
  common /blck06/ kw(nka),ka
  integer :: kw, ka
  common /cb18/ alat,declin                ! for the SZA calculation
  real (kind=dp) :: alat,declin

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw
  common /cb44/ a0m,b0m(nka)
  real (kind=dp) :: a0m,b0m
  common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb), &
                ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
  real (kind=dp) :: zb, dzb, dzbw, tb, eb, ak, d, &
       ajb, ajq, ajl, ajt, ajd, ajs, ds1, ds2, ajm, reif, tau, trdep
  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb51/ dlgew,dlgenw,dlne
  real (kind=dp) :: dlgew,dlgenw,dlne

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

! == End of declarations =======================================================

  write (jpfunprofm,6000)
6000 format (16x,'constants and parameters of the current run' &
          ,///,6x,'numerical grid',/,6x,'eta:')
  write (jpfunprofm,6010) (eta(k),k=1,n)
6010 format (1x,13e12.5)
  write (jpfunprofm,6020)
6020 format (6x,'etw:')
  write (jpfunprofm,6010) (etw(k),k=1,n)
  write (jpfunprofm,6030)
6030 format (6x,'deta:')
  write (jpfunprofm,6010) (deta(k),k=1,n)
  write (jpfunprofm,6040)
6040 format (6x,'detw:')
  write (jpfunprofm,6010) (detw(k),k=1,n)
  write (jpfunprofm,6050)
6050 format (6x,'enw:')
  write (jpfunprofm,6010) (enw(k),k=1,nka)
  write (jpfunprofm,6060)
6060 format (6x,'en:')
  write (jpfunprofm,6010) (en(k),k=1,nka)
  write (jpfunprofm,6070)
6070 format (6x,'rn:')
  write (jpfunprofm,6010) (rn(k),k=1,nka)
  write (jpfunprofm,6090)
6090 format (6x,'ew:')
  write (jpfunprofm,6010) (ew(k),k=1,nkt)
  write (jpfunprofm,6100)
6100 format (6x,'e:')
  write (jpfunprofm,6010) (e(k),k=1,nkt)
  write (jpfunprofm,6110)
6110 format (6x,'dew:')
  write (jpfunprofm,6010) (dew(k),k=1,nkt)
  write (jpfunprofm,6120)
6120 format (6x,'rq(k,1):')
  write (jpfunprofm,6010) (rq(k,1),k=1,nkt)
  write (jpfunprofm,6130)
6130 format (6x,'rq(k,nka):')
  write (jpfunprofm,6010) (rq(k,nka),k=1,nkt)
  if (isurf == 1) then
     write (jpfunprofm,6140)
6140 format (6x,'zb:')
     write (jpfunprofm,6010) (zb(k),k=1,nb)
     write (jpfunprofm,6150)
6150 format (6x,'dzb:')
     write (jpfunprofm,6010) (dzb(k),k=1,nb)
     write (jpfunprofm,6160)
6160 format (6x,'dzbw:')
     write (jpfunprofm,6010) (dzbw(k),k=1,nb)
  end if
  write (jpfunprofm,6170)
6170 format (//,16x,'constants and parameters'/)
  write (jpfunprofm,6173)
6173 format (6x,'aerosol type: 1=urban 2=rural 3=ocean 4=background')
  write (jpfunprofm,6176) (nar(k),k=1,n)
6176 format (6x,60i2)
  write (jpfunprofm,6177)
6177 format (6x,'solution term b0m of droplet growth equation')
  write (jpfunprofm,6010) (b0m(k),k=1,nka)
  write (jpfunprofm,6180)
6180 format (6x,'r0,r1,g,cp,a0m,ug,vg,dlgew,dlgenw,dlne')
  write (jpfunprofm,6010) r0,r1,g,cp,a0m,ug,vg,dlgew,dlgenw,dlne
  write (jpfunprofm,6190)
6190 format (6x,'z0')
  write (jpfunprofm,6010) z0
  write (jpfunprofm,6195) ka,kw
6195 format (//,16x,'ka and kw for aqueous phase reactions',21i6)
  write (jpfunprofm,6200)
6200 format (//,16x,'dimensions of arrays: n,nrlay,nrlev,nb,nka,nkt')
  write (jpfunprofm,6210) n,nrlay,nrlev,nb,nka,nkt
6210 format (/,1x,6i10)
  write (jpfunprofm,6220) alat,declin,wmin,wmax,tw,scaleo3_m
6220 format (//,6x,'geogr. latitude ',f9.1,' declination ',f9.1, &
          ' large scale subsidence in m/s',2f9.5, &
          ' water temperature',f8.2,' O3 column (hv only) ',f4.0,//)
  write (jpfunprofm,6230) chem,mic,rst
6230 format (6x,'current program evaluation: ','   chem: ',l1, &
          ' mic: ',l1,'   rst: ',l1,//)
  xxsum=0._dp
  xsum(:)=0._dp
  do k=2,nf
     do ia=1,nka
        do jt=1,nkt
           xsum(k)=xsum(k)+ff(jt,ia,k)*en(ia)
        enddo
     enddo
     xsum(k)=xsum(k)*1.e+09
     xxsum=xxsum+xsum(k)*detw(k)
  enddo

  write (jpfunprofm,6240)
6240 format (/,6x,'aerosol mass in ug m**-3 in layers 2 - nf')
  write (jpfunprofm,6250) xsum
6250 format (1x,15f8.3)
  write (jpfunprofm,6260) xxsum
6260 format(6x,'total aerosol mass in ug m**-2 of layers 2 - nf',f12.3)

end subroutine constm
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine constc
!
! Description:
!    output of chemical constants and parameters used in current run
!    This subroutine is called only once during initialisation (no restart case)
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      11/2016   Cleaned + header + implicit none         <Josue Bock>
!                    Removed outdated features
!
! 1.1       ?        Added advected species in the output     <Roland von Glasow>
!
! 1.0       ?        Original code.                           <Andreas Bott>
!


! Declarations:
! Modules used:
  USE config, ONLY : &
       neula

  USE file_unit, ONLY : &
       jpfunprofc

  USE gas_common, ONLY : &
! Imported Parameters:
       j1, &
       j5, &
       nadvmax, nindadv, xadv ! for Eulerian advection

  USE global_params, ONLY : &
! Imported Parameters:
       j2

  implicit none

! Local scalars:
  integer :: j

! == End of declarations =======================================================

  write (jpfunprofc,5900) neula
  write (jpfunprofc,5910) (nindadv(j),xadv(j),j=1,nadvmax)
5900 format ('euler (=0) or lagrangean view (=1): ',i3,' the following' &
          ,' species are advected only if neula=0')
5910 format (i3,d12.3)

  write (jpfunprofc,6090)
6090 format (//,16x,'dimension of arrays: j1, j2, j5')
  write (jpfunprofc,6100) j1,j2,j5
6100 format (/,1x,5i10,//)

end subroutine constc
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine profm (dt)

! output of meteorological profiles

! Declarations:
! Modules used:

  USE config, ONLY : &
       isurf

  USE data_surface, ONLY : &
       ustern, z0                  ! frictional velocity, roughness length

  USE file_unit, ONLY : &
       jpfunprofm

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nrlay, &
       nb, &
       nka, &
       nkt, &
       mbs

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real(kind=dp), intent(in) :: dt

! Local scalars:
  character (len=10) :: srname
  integer :: ia, jt, k
  real(kind=dp) :: xxm1, xxm2, xxsum
! Local arrays:
  real(kind=dp) :: xsum(n)

! Common blocks:
  common /cb16/ u0,albedo(mbs),thk(nrlay)
  real (kind=dp) :: u0, albedo, thk

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy
  common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
  real (kind=dp) :: gm, gh, sm, sh, xl
  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w
  common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb), &
                ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
  real (kind=dp) :: zb, dzb, dzbw, tb, eb, ak, d, &
       ajb, ajq, ajl, ajt, ajd, ajs, ds1, ds2, ajm, reif, tau, trdep
  common /cb48/ sk,sl,dtrad(n),dtcon(n)
  real (kind=dp) :: sk, sl, dtrad, dtcon

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a

! == End of declarations =======================================================

  srname='          '
  write (jpfunprofm,6000) it,dt,lday,lst,lmin
6000 format (//,6x,i8,'-th. timestep dt = ',f4.1,' sec ',i2,' day ', &
          i2,' hour ',i2,' min '/)
  write (jpfunprofm,6010)
6010 format (1x,'k  height  d/dz(thetl)*100  atkm  atkh  xl', &
          '  tke  tkep  thetl   d/dt(dtrad)*3600   d/dt(dtcon)*3600')
  write (jpfunprofm,6020) (k,etw(k),(thetl(k+1)-thetl(k))/deta(k)*100._dp, &
       atkm(k),atkh(k),xl(k), &
       tke(k),tkep(k),thetl(k),dtrad(k)*3600._dp,dtcon(k)*3600._dp,k=n-1,1,-1)
6020 format (1x,i3,f7.1,3f9.3,f7.1,f9.4,e12.4,f9.3,2f12.4)
  write (jpfunprofm,6030)
6030 format (/,1x,'k  height       u         v         t          p', &
          '      theta      feu          q        m2        fsum')
  write (jpfunprofm,6040) (k,eta(k),u(k),v(k), &
       t(k)-273.15_dp,p(k)/100._dp,theta(k), &
       feu(k)*100._dp,1000._dp*xm1(k),1000._dp*xm2(k), &
       fsum(k),k=n,1,-1)
  xxm1=0._dp
  xxm2=0._dp
  do k=1,n
!     xxm1=xm1(k)*detw(k)+xxm1
     xxm1=xm1(k)*detw(k)*rho(k)*1000+xxm1
     xxm2=xm2(k)*detw(k)*1000+xxm2
! xxm1 in g/m**2 vapour content of atm., xxm2 in g/m**2 liquid water content of atm.
  enddo
6040 format (1x,i3,f10.1,9f10.3)
  write (jpfunprofm,6050) xxm1,xxm2,tau*1000._dp,reif*1000._dp,trdep*1000._dp,ds1*1000._dp,ds2*1000._dp
6050 format (/,1x,'sum(xm1*detw*rho)',f10.3,1x,'sum(xm2*detw)',f10.3, &
          1x,'dew',f10.3,5x,'rime',f10.3,5x,'particles',f10.3,5x,'aerosol', &
          f10.3,5x,'droplets',f10.3)
  write (jpfunprofm,6060) u0,sk,sl,-5.6697d-8*t(1)**4
6060 format (1x,'surface radiative fluxes'/,10x,'u0: ',f10.4, &
          10x,'solar: ',e11.4,3x,'infrared: ',e11.4,3x,'emission: ',e11.4)
  write (jpfunprofm,6070) z0,ustern,ajq,ajs,ajm,ajd
6070 format (1x,'roughness length ',e10.3,' friction velocity u*', &
          f8.3,/,1x,'surface moisture fluxes'/, &
          10x,'water vapor: ',e11.4,3x,'droplet sedimentation: ',e11.4,3x, &
          'ground moisture: ',e11.4,3x,'dew storage: ',e11.4)
  write (jpfunprofm,6080) ajb,ajl,ajt,sk+sl-5.669d-8*t(1)**4
6080 format (1x,'surface heat fluxes'/, &
          10x,'ground heat: ',e11.4,3x,'latent heat: ',e11.4,3x, &
          'sensible heat: ',e11.4,3x,'net radiation: ',e11.4)

  if (isurf == 1) then
     write (jpfunprofm,6090)
6090 format (/,1x,'temperature and volumetric moisture content in', &
          ' ground:')
     write (jpfunprofm,6100) (zb(k),k=1,nb)
     write (jpfunprofm,6110) (tb(k)-273.15_dp,k=1,nb)
     write (jpfunprofm,6120) (eb(k),k=1,nb)
6100 format (4x,'zb:',/,10f10.3,/,10f10.3)
6110 format (4x,'tb:',/,10f10.3,/,10f10.3)
6120 format (4x,'eb:',/,10f10.3,/,10f10.3)
  end if

  xxsum=0._dp
  xsum(:)=0._dp
  do k=2,nf
     do ia=1,nka
        do jt=1,nkt
           xsum(k)=xsum(k)+ff(jt,ia,k)*en(ia)
        enddo
     enddo
     xsum(k)=xsum(k)*1.e+09
     xxsum=xxsum+xsum(k)*detw(k)
  enddo
  write (jpfunprofm,6240)
6240 format (/,6x,'aerosol mass in ug m**-3 in layers 2 - nf')
  write (jpfunprofm,6250) xsum
6250 format (1x,15f8.3)
  write (jpfunprofm,6260) xxsum
6260 format(6x,'total aerosol mass in ug m**-2 of layers 2 - nf',f12.3)

  call ion_mass (srname)

! 141  format (2i3,9d12.4)
! 142  format (6x,9d12.4)

!  do k=2,nf
!     do kc=1,nkc_l
!        write (*,141) (k,kc,(dss(k,l,kc),l=1,lsp))
!        write (*,142) (svc(k,kc,kkc),kkc=1,nkc_l)
!!        write (*,142) (fss(k,kc,1),fss(k,kc,2),svc(k,kc,1),svc(k,kc,2))
!     enddo
!  enddo

end subroutine profm
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine profc (dt,mic)

! vertical profiles of chemical species

! technical comment:
!     when writing the blocks headers ('height' + 'species_names') the format
!     can include 10a12 even if less than 10 species names are written.
!     However, when writting the arrays (one line per model level), the exact
!     number of species has to be used in the format descriptor (fmt), because
!     it is repeated for each line.

  USE config, ONLY: &
       nkc_l

  USE file_unit, ONLY : &
       jpfunprofc

  USE gas_common, ONLY : &
! Imported parameter
       j1, &
       j5, &
! Imported Array Variables with intent (in):
       s1, &
       s3, &
       gas_name, &
       rad_name

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j3, &
       j6, &
       nkc, &
       nlev, &
       nrxn, &
       nf, &
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  real (kind=dp), intent(in) :: dt
  logical, intent(in) :: mic

  real (kind=dp) :: si(10,n),xfac(nf)

  character (len=17) :: fmt
  integer :: k, kc, l
  integer :: jblock, jlay, jspec, jspec_max
  integer :: nrate
  real (kind=dp) :: zx0(n) ! conversion factor

  common /blck01/ am3(n),cm3(n)
  real (kind=dp) :: am3, cm3

  common /blck11/ rc(nkc,n)
  real (kind=dp) :: rc

  common /blck12/ cw(nkc,n),cm(nkc,n)
  real (kind=dp) :: cw, cm

  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1

  common /budg/ bg(2,nrxn,nlev),il(nlev) ! jjb corrected, was written il(nrxn)
  real (kind=dp) :: bg
  integer :: il

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw
! == End of declarations =======================================================


! ==============================================================================
!  -- 1.1 -- Gas phase species: header and common conversion factor
! ==============================================================================

!     Write header
!     ------------
  write (jpfunprofc,6000) it,dt,lday,lst,lmin
6000 format (//,6x,i8,'-th. timestep dt = ',f4.1,' sec ',i2,' day ', &
          i2,' hour ',i2,' min '/)
  write (jpfunprofc,6010)
6010 format (10x,'gas phase species in ppb; at the ground total', &
          ' deposition in molecules/cm**2')


!     output in ppb
!     define conversion factor (zx0):
  do jlay=2,n
     zx0(jlay) = 1.e+09/am3(jlay)
  end do



! ==============================================================================
!  -- 1.2 -- Non radical, gas phase species
! ==============================================================================

!     Write data blocks of 10 species
!     -------------------------------
  jblock = 0
  do while (10*jblock < j1)
     jspec_max = MIN(10,j1-10*jblock) ! Number of output species in the current block (10 or less)

     ! No conversion for layer 1
     do jspec=1,jspec_max
        si(jspec,1)=s1(10*jblock+jspec,1)
     enddo

     do jlay=2,n
        do jspec=1,jspec_max
           si(jspec,jlay)=s1(10*jblock+jspec,jlay)*zx0(jlay)
        end do
     end do

     ! write gas names
     write (jpfunprofc,6015)'height',(TRIM(gas_name(10*jblock+jspec)),jspec=1,jspec_max)
6015 format (3x,a7,10a12)
     ! write gas concentrations (and deposition at the ground)
     write (fmt,'( "(f10.1,",i2,"e12.5)" )')jspec_max
     write (jpfunprofc,fmt) (eta(jlay),(si(jspec,jlay),jspec=1,jspec_max),jlay=n,1,-1)

     jblock = jblock + 1
  end do



! ==============================================================================
!  -- 1.3 -- Radical, gas phase species
! ==============================================================================

!     Write data blocks of 10 species
!     -------------------------------
  jblock = 0
  do while (10*jblock < j5)
     jspec_max = MIN(10,j5-10*jblock) ! Number of output species in the current block (10 or less)

     ! no deposition for radicals, start at layer 2
     do jlay=2,n
        do jspec=1,jspec_max
           si(jspec,jlay)=s3(10*jblock+jspec,jlay)*zx0(jlay)
        end do
     end do

     ! write radicals names
     write (jpfunprofc,6015)'height',(TRIM(rad_name(10*jblock+jspec)),jspec=1,jspec_max)
     ! write radicals concentrations
     write (fmt,'( "(f10.1,",i2,"e12.5)" )')jspec_max
     write (jpfunprofc,fmt) (eta(jlay),(si(jspec,jlay),jspec=1,jspec_max),jlay=n,2,-1)

     jblock = jblock + 1
  end do



! ==============================================================================
!  -- 2345 -- to be done
! ==============================================================================

!     output of reaction rates integrated (later: over 1 hour), converted to mol/(mol*h)
  write (jpfunprofc,6170)
  write (jpfunprofc,6180) eta(il(1)),eta(il(2)),eta(il(3))
  nrate=600
  if (nkc_l.eq.4) nrate=1120
  do l=1,nrate
     write (jpfunprofc,6190) l,bg(2,l,1)/am3(il(1)),bg(1,l,1)/am3(il(1)), &
          bg(2,l,2)/am3(il(2)),bg(1,l,2)/am3(il(2)), &
          bg(2,l,3)/am3(il(3)),bg(1,l,3)/am3(il(3))
  enddo
! 6170 format (//,'reaction rates integrated over 1 hour, converted to'
6170 format (//,'accumulated reaction rates [mol/(m^3(air) s)]')
6180 format (/,'height = ',19x,f10.2,' m',4x,f10.2,' m',4x,f10.2,' m')

! 6190 format ('reaction no. ',i4,' rate : ',4d16.4)
6190 format ('no. ',i4,' : ',6d16.8)



!  if (.not.mic.or.lct.le.1) return
  if (.not.mic) return
! liquid phase chemistry
  do kc=1,nkc_l
     write (jpfunprofc,6040) kc
6040 format (/,10x,'aqueous phase species in mole/m**3;', &
          ' at the ground total deposition in mole/m**2 for bin:',i3,//, &
          4x,'height',5x,'cw',7x,'rc',7x,'hno3',6x,'nh3',7x,'so2',5x, &
          'h2so4',7x,'o3',9x,'h2o2',7x,'fe(3)',6x,'mn(2)')
     write (jpfunprofc,6020) (eta(k),cw(kc,k),rc(kc,k),sl1(3,kc,k), &
          sl1(4,kc,k),sl1(5,kc,k),sl1(6,kc,k),sl1(7,kc,k), &
          sl1(19,kc,k),sl1(44,kc,k),sl1(45,kc,k),k=nf,1,-1)
!          sl1(19,kc,k),sl1(44,kc,k),sl1(45,kc,k),k=lct,lcl,-1)
     write (jpfunprofc,6041)
6041 format (4x,'height',5x,'OH',7x,'HO2',7x,'NO3',6x,'NO',7x,'NO2',5x, &
          'HONO',7x,'HCHO')
     write (jpfunprofc,6042) (eta(k),sl1(j2-j3+2,kc,k),sl1(j2-j3+3,kc,k), &
          sl1(36,kc,k),sl1(1,kc,k),sl1(2,kc,k),sl1(21,kc,k),sl1(17,kc,k), &
          k=nf,1,-1)
!          k=lct,lcl,-1)
6042 format (f10.1,7e10.3)
     write (jpfunprofc,6050) kc
! 6050 format (/,20x,'ion concentrations in mole/liter',//, &
6050 format (/,20x,'ion concentrations in mol/m^3 for bin:',i3,//, &
!          4x,'height',7x,'h+',7x,'nh4+',6x,'cl-',2x,'ch2ohso3-',2x, &
          4x,'height',7x,'h+',7x,'nh4+',6x,'cl-',2x,'Br-',2x, &
!          3x,'hso3-',6x,'so3=',6x,'so4-',5x,'so4=',5x,'no3-',6x,'fe')
          3x,'hso3-',6x,'so3=',6x,'so4-',5x,'so4=',5x,'no3-',6x,'tracer')
     do k=1,nf
        xfac(k)=1._dp             !mol/m^3_air       --> mol/m^3_air
     enddo

     write (jpfunprofc,6020) (eta(k),sion1(1,kc,k)*xfac(k),sion1(2,kc,k)* &
          xfac(k),sion1(14,kc,k)*xfac(k),sion1(24,kc,k)*xfac(k), &
          (sion1(l,kc,k)*xfac(k),l=5,8), &
!          sion1(17,kc,k)*xfac(k),k=lct,lcl,-1)
!          sion1(13,kc,k)*xfac(k),sion1(3,kc,k)*xfac(k),k=lct,lcl,-1)
          sion1(13,kc,k)*xfac(k),sion1(21,kc,k),k=nf,1,-1)
6020 format (f10.1,10e10.3)
  end do

  write (jpfunprofc,*) 'done with profc'

end subroutine profc
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine profr

! Description :
! -----------
!   output of radiation variables


! Interface :
! ---------

! Input :
! -----

! Output :
! ------

! Externals :
! ---------
!    none


! Author :
! ------
!    Andreas Bott


! Modifications :
! -------------
!           2016  Josue Bock  removed hard coded parameters, use module instead
!                             correction and homogeneisation of array size.
!    16-Mar-2017  Josue Bock  missing declarations, implicit none
!    18-Mar-2017  Josue Bock  improvement of output formats, to fit the data
!                             added header, comments
!    19-Oct-2017  Josue Bock  double precision => real(kind=dp)
!    13-Nov-2017  Josue Bock  removed ntypa and ntypd from /cb02/, unused

! == End of header =============================================================

! Declarations:
! ------------
! Modules used:
  USE file_unit, ONLY : &
       jpfunprofr

  USE global_params, ONLY : &
! Imported Parameters:
       nrlay, &
       nrlev, &
       mb, &
       mbs

  USE precision, ONLY : &
       dp

  implicit none

! Local scalars:
  integer :: i, ib

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real(kind=dp) :: t,p,rho,xm1,ts

  common /cb10/ totrad (mb,nrlay)
  real(kind=dp) :: totrad

  common /cb15/ fnseb,flgeg,hr(nrlay)
  real(kind=dp) :: fnseb, flgeg, hr

  common /cb16/ u0,albedo(mbs),thk(nrlay)
  real(kind=dp) :: u0, albedo, thk

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real(kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /kurz/ fs1(nrlev),fs2(nrlev),totds(nrlev),ss(nrlev), &
                fsn(nrlev),dtdts(nrlay)
  real(kind=dp) :: fs1, fs2, totds, ss, fsn, dtdts

  common /lang/ fl1(nrlev),fl2(nrlev),fln(nrlev),dtdtl(nrlay)
  real(kind=dp) :: fl1, fl2, fln, dtdtl

! == End of declarations =======================================================


  ! Time stamp
  write (jpfunprofr,6000) lday,lst,lmin,u0

  ! Solar bands (middle of layer (or half level) values)
  write (jpfunprofr,6010)
  do i=1,nrlay
     write (jpfunprofr,6011) i,(totrad(ib,i),ib=1,mbs),hr(i)
  enddo

  ! IR bands (middle of layer (or half level) values)
  write (jpfunprofr,6020)
  do i=1,nrlay
     write (jpfunprofr,6021) i,(totrad(ib,i),ib=mbs+1,mb)
  enddo

  ! P, T, and fluxes (level values)
  write (jpfunprofr,6030)
  do i=1,nrlev
     write (jpfunprofr,6031) i,p(i),t(i),fs1(i),fs2(i),ss(i),fl1(i),fl2(i)
  enddo


! Formats :
! -------
6000 format (/,'day: ',i3,10x,'hour: ',i3,10x,'minute: ',i3, &
          10x,'cosine of zenith distance: ',f8.2,/)

6010 format (/,'#layer',20x,'totrad(l,i) l=1,6 (solar)',20x,'hr')
6011 format (i4,6f10.3,es14.4)

6020 format (/,'#layer',38x,'totrad(l,i) l=7,18 (IR)')
6021 format (i4,12f10.3)

6030 format (/,'#level',4x,'pres',4x,'temp',10x,'fs1',10x,'fs2',12x, &
          'ss',11x,'fl1',11x,'fl2')
6031 format(i4,f10.1,f8.1,5f14.3)


end subroutine profr
