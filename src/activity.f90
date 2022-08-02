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

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! activity.f:  calculation of activity coefficients after Pitzer

!  called by SR activ (in kpp.f)
!              |
!              |___SR pitzer
!                    |______SR calpar
!                    |______FN gammann
!                    |        |_______FN g
!                    |        |_______FN gs
!                    |        |_______SR efunct
!                    |
!                    |______FN gammasn
!                             |_______SR efunct
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Authors:
! --------
!    < Beiping Luo >

!      + contributions from:
!         < Roland von Glasow >  adapted for Mistra
!         < Josue Bock >         partly re-written, corrected, cleaned and improved*
!
!         * : technical improvements only, scientific content neither changed nor checked.

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine pitzer (k,kc,wact)
!
! Description:

!      This only a simple version, which only take the unsymmetrical
!      factor e-theta and E-theta ' into account.

!      Anion list:
!         1 HSO4-
!         2 SO4=
!         3 NO3-
!         4 Cl-

!      Cation list:
!         1 H+
!         2 NH4+
!         3 Na+
!

! History:
! Version   Date     Comment
! -------   ----     -------
!           1990     Original code.                <Beiping Luo>
!
!           ?        Adapted for Mistra            <Roland von Glasow>
!
!         2016/2017  Multiple code improvement     <Josue Bock>
!                     to prepare release version:
!                     - general file header
!                     - all subroutines headers
!                     - removal of unused variables
!                     - all missing declarations and implicit none everywhere
!                     - lots of cleaning, properly indenting
!                     - integer vs real conflicts corrected in the whole module
!                     - use module for parameters
!                     - updated Fortran features such as intent, external
!                     - several bugfix

! Recent changes:
! ---------------
!      /  /2016  <jjb>  in case N==1, removed unused statement rh=rh*exp(...) where rh uninitialised

!    03/02/2017  <jjb>  xs array size is only 11 instead of 100, corrected in calpar as well

!    03/02/2017  <jjb>  bugfix: c1(3,:) was not defined along with b0, b1, and c0 (see after call calpar).
!                           Now set to 0. This was causing troubles when running the model with gfortran

!    03/02/2017  <jjb>  defined the dimension of all arrays using NC, NA instead of nmax (removed)

!    03/02/2017  <jjb>  use Mistra indexes for species from hardcoded "ionind" list

!    03/03/2017  <jjb>  more cleaning / consistency of variable naming with other places
!                       computing efficiency : if / elseif / endif  insead of repeated if ... endif
!                       za and zc are integers: changed declaration, use "real" to convert kind when necessary
!                       calculation of ZI, I and I2 only once, passed to SRs as arguments

!     01/11/2017 <jjb> Fortran90
!                      changed T0 => T1 to be consistent throughout the code


! Declarations:
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       j2,                  &
       j6,                  &
       nf,                  &
       n,                   &
       nkc

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
  integer,       intent(in)  :: k, kc ! Layer index, liquid bin index
  real(kind=dp), intent(out) :: wact

! Common blocks
  common /blck12/ cw(nkc,n),cm(nkc,n)                       ! LWC input
  real(kind=dp) :: cw, cm

  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)             ! ion mixing ratio
  real(kind=dp) :: sl1, sion1

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)  ! get temp from MISTRA
  real(kind=dp) :: theta, thetl, t, talt, p, rho

  common /kpp_mol/ xgamma(j6,nkc,nf)                        ! gamma output
  real(kind=dp) :: xgamma

! External functions:
  real(kind=dp), external :: gammann
  real(kind=dp), external :: gammasn   ! water activity

! Local parameters:
  integer,parameter :: NA = 4 ! number of anions
  integer,parameter :: NC = 3 ! number of cations
  real(kind=dp), parameter :: T1 = 298.15_dp

! Local scalars:
  integer :: ia, ic                 ! loop indexes
  integer :: NI                     ! loop index, index of the ion whose activity is calculated
  real(kind=dp) :: g
  real(kind=dp) :: I, I2            ! I = ionic strength
  real(kind=dp) :: rhmix1, r2mix1
  real(kind=dp) :: TK               ! temperature in layer k (in [K])
  real(kind=dp) :: xam, xap         ! sum of negative / positive charges (signed value)
  real(kind=dp) :: xmix3, xu2, xu3
  real(kind=dp) :: ZI               ! ZI = sum (mc*zc) + sum(ma*za)

! Local arrays:
  integer       :: Iflag(nc,na)
  real(kind=dp) :: b0(nc,na),B1(nc,na)  ! Pitzer coefficients
  real(kind=dp) :: c0(nc,na),C1(nc,na)
  real(kind=dp) :: MC(nc),MA(na)        ! concentration of cations and anions
  integer       :: ZC(nc),ZA(na)        ! ionic charge of cations and anions
  real(kind=dp) :: omega(nc,na)
  real(kind=dp) :: xs(11)

!  character (len=6) :: ionn(7)        ! jjb currently unused, but left commented for later use
!                                      !     could be tested along with Mistra ion names (index)
!  data ionn /'H+', 'NH4+', 'Na+', 'HSO4-','SO4--','NO3-', 'Cl-'/

  integer, parameter :: ionind(NC+NA) = (/1,2,20,19,8,13,14/) !indices of ions in sion1

  data ZC /1, 1, 1/
  data ZA /1, 2, 1, 1/

!- End of header ---------------------------------------------------------------


! -- Iflag definition, used in calpar
  Iflag(1,1)=3
  Iflag(1,2)=4
  Iflag(1,3)=1
  Iflag(1,4)=2

  Iflag(2,1)=5
  Iflag(2,2)=6
  Iflag(2,3)=7
  Iflag(2,4)=8

  Iflag(3,1)=0 ! <jjb> added this for security, even if never used.
  Iflag(3,2)=0 !        If a mistake was done in calpar, this is safer.
  Iflag(3,3)=0 !        See notes below for complete story
  Iflag(3,4)=0

  TK=t(k)            ! temperature in layer k (in [K])

! concentrations --> molality:
!  mol/m^3_air * m^3_air/m^3_solvent * 10^-3 * 1 dm^3/kg (density of water) = mol/kg_solvent
!   sion1      * cm^-1               * 1.d-3

! input molality of cations
  ! NB : cm /= 0. (already checked in the calling SR: activ, in kpp.f)
  do ic=1,nc
     mc(ic) = sion1(ionind(ic),kc,k)*1.e-3_dp/cm(kc,k)
  enddo

! input of molality of anions
  do ia=1,na
     ma(ia) = sion1(ionind(NC+ia),kc,k)*1.e-3_dp/cm(kc,k)
  enddo

! -- Calculate xap, xam, ZI and I
  I = 0._dp
  xap = 0._dp
  do ic=1,NC
     I = I + mc(ic)*real(ZC(ic)**2,dp)
     xap = xap + mc(ic)*real(zc(ic),dp)
  enddo
  xam=0.
  do IA=1,NA
     I = I + ma(ia)*real(ZA(ia)**2,dp)
     xam = xam - ma(ia)*real(za(ia),dp)
  enddo
  I=0.5_dp*I
  I2=sqrt(I)

  ZI = xap - xam

!x! charge balance
!x      nadd=0
!x!     if( abs(xam+xap).ge.1.e-10_dp) then
!x         if( abs(xam+xap).ge.1.e-5_dp) then
!x 2010       continue
!x            nadd=nadd+1
!xcd            print *,xam,xap,abs(xam+xap)
!x            if (abs(xam).gt.abs(xap)) then
!x               xc=abs(xap/xam)
!x               do kn=1,na
!x                  ma(kn)=ma(kn)*xc
!x               enddo
!x            else
!x               xc=abs(xam/xap)
!x               do kn=1,nc
!x                  mc(kn)=mc(kn)*xc
!x               enddo
!x            endif
!x            do I=1,nc
!x               xap=xap+zc(I)*mc(I)
!xcd               print *,ii, mc(I)
!x             enddo
!x             do I=1,na
!x                xam=xam-za(I)*ma(I)
!xcd                print *,ii, mA(I)
!x             enddo
!x!             print *,xam,xap,abs(xam+xap),xc
!x!             print*, ' THE SOLUTION IS NOT NEUTRAL !'
!x!             if( abs(xam+xap).ge.1.e-10_dp) goto 2010
!x             if( (1-xc).ge.1.e-4_dp) goto 2010
!x!             goto 1000
!x         endif


  call calpar(Iflag,NC,NA,TK, b0,b1,c0,c1,omega,xs) ! <jjb> note NC is used here, to get array size in calpar
                                                    !      however, NC=3 values are set below, not in calpar
                                                    !      thus, NC-1=2 is used in calpar do loop.

! get the parameters for Na-anions (NA=3)
  ! -- Na-HSO4 ! Die Temperatureabhaengigkeit fuer Na-HSO4 fehlt noch.
  b0(3,1)=.0454_dp
  b1(3,1)=.398_dp
  c0(3,1)=0._dp
  c1(3,1)=0._dp
  omega(3,1)=2._dp
  ! -- Na-SO4
  b0(3,2)=0.0261_dp + (TK-T1)* 2.36e-3_dp
  b1(3,2)=1.484_dp  + (TK-T1)* 5.63e-3_dp
  c0(3,2)=.00938_dp - (TK-T1)* .172e-3_dp
  c1(3,2)=0._dp
  omega(3,2)=2._dp
  ! -- Na-NO3
  b0(3,3)=.0068_dp + (TK-T1)* 12.66e-4_dp
  b1(3,3)=.1783_dp + (TK-T1)* 20.6e-4_dp
  c0(3,3)=-.00072_dp/2._dp - (TK-T1)* 23.16e-5_dp/2._dp
  c1(3,3)=0._dp
  omega(3,3)=2._dp
  ! -- NaCl
  b0(3,4)=0.0765_dp + (TK-T1)* 7.159e-4_dp
  b1(3,4)=.2664_dp  + (TK-T1)* 7.e-4_dp
  c0(3,4)=.00127_dp/2._dp - (TK-T1)* 10.5e-5_dp/2_dp
  c1(3,4)=0._dp
  omega(3,4)=2._dp

! -- Calculate ions activities
  do NI=1, NC+NA
     ! gammann calculate on the binary terms
     g=gammann(TK,NI,NC,NA,mC,mA,zC,zA,ZI,I,I2,b0,b1,C0,C1,omega)

     ! calculate the ternary terms of H+
     if( NI.eq.1) then
        xmix3=ma(2)*mc(2)*xs(7)+mC(2)*mA(1)*xs(9)
        rhmix1=ma(1)*ma(3)*xs(1)+xs(3)*ma(2)*ma(3)
        g=g*exp(rhmix1)

     ! calculate the ternary terms of NH4+
     else if(NI.eq.2) then
        xmix3=ma(2)*mc(1)*xs(8)+mA(1)*ma(2)*xs(7)+ma(1)*mc(1)*xs(9)
        g=g*exp((xmix3)+2*mc(1)*xs(10))


     ! calculate the ternary terms of HSO4-
     else if(NI-NC.eq. 1) then
        xu2=mc(1)*ma(3)*xs(1)
        xu2=xu2+mc(1)*ma(4)*xs(4) + ma(4)*xs(5)*2
        ! mixture of NH4,SO4,HSO4
        xmix3=ma(2)*mc(2)*xs(7)+mc(1)*MC(2)*xs(9)
        g=g*exp(xu2+xmix3)

     ! calculate the ternary terms of SO4--
     else if(NI-NC.eq. 2) then
        xu3=ma(3)*mc(1)*xs(3) + ma(3)*xs(2)*2
        xu3=xu3+ma(4)*mc(1)*xs(6)
        xmix3=ma(1)*mc(2)*xs(7) +mc(1)*mc(2)*xs(8)
        g=g*exp(xu3+xmix3)

     ! calculate the ternary terms of NO3-
     else if(NI-NC.eq. 3) then
        r2mix1=ma(1)*mc(1)*xs(1)   & !H,HSO4,NO3
              + 2._dp*ma(2)*xs(2)  & !SO4,NO3
              + ma(2)*mc(1)*xs(3)  & ! H,SO4,NO3
              + mc(2)*ma(2)*xs(11)   ! NH4,SO4,NO3
        g=g*exp(r2mix1)
     endif

     xgamma(ionind(NI),kc,k)=g

  end do

! -- Calculate water activity      
  wact=GAMMASN(TK,NC,NA,MC,MA,ZC,ZA,ZI,I,I2,B0,B1,C0,C1,omega)


end subroutine pitzer

!------------------------------------------------------------------

subroutine calpar(Iflag,NC,NA,T,b0,b1,C0,C1,omega,xs)

! 03/02/2017 <jjb> corrected array sizes:
!     - bb(27) instead of bb(50)
!     - b2(21) instead of b2(50)
!     - b3(42) instead of b3(50)
!     - b5(31) instead of b5(50) ! notice unused indexes
!     - b6(13) instead of b6(50)
!     - b7(13) instead of b7(50)
!     - b8(9)  instead of b8(50)
!     - b9(17) instead of b9(50) ! unused, thus commented
!
! 03/02/2017 <jjb> aligned data blocks for clarity and ease futur check
!                  removed duplicated definitions of dt inside the if cases
!                  declared all variables, implicit none
!
! 03/02/2017 <jjb> defined the dimension of all arrays using NC, NA instead of a larger nmax (removed)
!                  changed the main do-loop bound, from NC to NC-1 (previously, calpar was called with
!                  a "fake" NC = NC-1). Now, called with the "real" NC value.
!                  Data for NC = 3 are in SR Pitzer.
!
! 03/03/2017 <jjb> xs(1:3) and xs(11) are defined twice the same (H-NO3 and NH4-Cl cases)
!                  commented the second, no need to overwrite with same values
!                  Also removed the initialisation (data xs /11*0./) in pitzer header: useless
!
! 04/03/2017 <jjb> introduced dt2, dt3, dt4 to avoid multiple ** calculations
!
! 01/11/2017 <jjb> Fortran90
!                  changed T0 => T1 to be consistent throughout the code

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  integer,       intent(in) :: NA, NC
  real(kind=dp), intent(in) :: T
! Array arguments with intent(in):
  integer,       intent(in) :: Iflag(NC,NA)

! Array arguments with intent(out):
  real(kind=dp), intent(out) :: b0(nc,na), b1(nc,na)
  real(kind=dp), intent(out) :: c0(nc,na), c1(nc,na)

  real(kind=dp), intent(out) :: omega(nc,na)
  real(kind=dp), intent(out) :: xs(11)

! Local parameters:
  real(kind=dp), parameter :: T1 = 298.15_dp

! Local scalars:
  real(kind=dp) :: dt, dt2, dt3, dt4
  integer :: i, j

! Local arrays:
  real(kind=dp) :: bb(27)
  real(kind=dp) :: b2(21)
  real(kind=dp) :: b3(42)
  real(kind=dp) :: b5(31)
  real(kind=dp) :: b6(13)
  real(kind=dp) :: b7(13)
  real(kind=dp) :: b8(9)
!  real(kind=dp) :: b9(17) ! <jjb> unused

!- End of header ---------------------------------------------------------------


! ----------- data for H-NO3 -------------------------------
  data (bb(I),I=1,27)/ &
    3.895835e-3_dp, -1.55571e-2_dp,  1.703729e-2_dp, -5.6173712e-3_dp, &
    5.732047e-3_dp,  0.91622_dp,     0.613523_dp,    -0.68489_dp,      &
    0.3038_dp,      -0.32888_dp,     7.6086113e-7_dp, 7.2714678e-5_dp, &
   -1.0037e-4_dp,    3.475e-5_dp,   -3.62927e-5_dp,   5.380465e-2_dp,  &
   -2.2163e-2_dp,   -1.0166e-2_dp,   6.5423e-3_dp,   -8.80248e-3_dp,   &
    0.907342_dp,    -6.78428e-4_dp,  9.576e-4_dp,     0._dp,           &
    0._dp,           7.769e-3_dp,   -5.819e-4_dp /

! ---------- data for H-Cl ---------------------------------
  data (b2(I),i=1,21) / &
    0.23378_dp,     -7.21238e-2_dp,  -1.7335667e-2_dp, 5.760665e-3_dp, &
   -8.29279e-3_dp,   0.2897_dp,       7.575434e-2_dp, -1.1474e-3_dp,   &
    0.38038_dp,     -0.309442_dp,    -2.794885e-3_dp,  2.309349e-4_dp, &
    9.322982e-4_dp, -2.398e-4_dp,     2.85959e-4_dp,  -0.21154_dp,     &
    0.101481_dp,     5.945618e-2_dp, -0.107864_dp,     8.81749e-2_dp,  &
    1.9916_dp /

! ----------- data for H-HSO4, H-SO4 -----------------------
  data (b3(I),i=1,42) /                0.148843_dp,    -7.769e-2_dp,     &
    2.8062e-2_dp,    4.7903e-4_dp,     7.25e-4_dp,      0.17843_dp,      &
    0.678_dp,        8.7381e-2_dp,    -0.57881_dp,      7.58e-2_dp,      &
   -9.878e-4_dp,     5.447651e-4_dp,  -2.58798e-4_dp,   1.8466527e-5_dp, &
    1.23457e-5_dp,   0.37138_dp,      -9.24874e-2_dp,  -9.21372e-3_dp,   &
   -1.065158e-2_dp,  5.4987733e-2_dp,  0.2726312_dp,   -1.34824e-3_dp,   &
   -0.24711_dp,      1.25978e-2_dp,    0.11919_dp,      0.7397_dp,       &
   -3.01755_dp,     -4.5305_dp,       -3.1072_dp,      -0.8555842_dp,    &
    9.2223e-4_dp,   -4.1694532e-3_dp,  7.141266e-3_dp,  2.32984e-3_dp,   &
   -6.98191e-4_dp,  -2.242_dp,         0.71925_dp,      2.52_dp,         &
   -0.7391_dp,      -1.548503_dp,      1.5452_dp,       2._dp /

! ----------- data for NH4-HSO4 ----------------------------
! -- without chan's data
  data (b5(I),I=1,31)/ &
   -8.746e-4_dp,  -2.3125_dp,    -9.56785e-6_dp, 2.58238_dp,    2.38_dp,      &
   -3.1314e-4_dp,  1.6896e-2_dp, -0.7351_dp,     0.6883_dp,     1.813e-3_dp,  &   ! indexes 8 & 9 never used
   -0.1012515_dp, -2.66e-2_dp,   -2.86617e-3_dp, 0.22925_dp,    0.438188_dp,  &
    2.522e-4_dp,  -2.90117e-5_dp, 0.9014_dp,     0.41774_dp,   -1035.9_dp,    &   ! indexes 20 - 23 never used
    0.0_dp,       -299.69_dp,     0.0_dp,       -4.9687e-4_dp,  0.0_dp,       &   !
    1.21485e-2_dp, 0.0_dp,       -1.0334e-3_dp,  0.0_dp,        8.48374e-2_dp,&
    0.0_dp /

! -- with chan's data
!  data (b5(I),I=1,31)/ &
!   -7.8224e-3_dp,  -1.722_dp,      7.5882e-5_dp,  0.962_dp,      1.8914_dp,     &
!   -4.285e-4_dp,    3.99e-4_dp,   -0.76151_dp,    0.53133_dp,    1.16726e-3_dp, &  ! indexes 8 & 9 never used
!   -9.45855e-2_dp, -2.007e-2_dp,   6.0482e-3_dp, -0.1392_dp,     2.9544_dp,     &
!    2.45133e-4_dp, -5.492e-5_dp,   0.8358_dp,    -0.3791_dp,    -980.81_dp,     &  ! indexes 20 - 23 never used
!    0._dp,          265.7_dp,      0._dp,        -1.7568e-3_dp, -7.72e-4_dp,    &  !
!    3.437e-4_dp,   -1.6456e-2_dp, -1.2034e-3_dp, -2.2943e-3_dp,  7.9281e-2_dp,  &
!    3.772e-2_dp /

! ----------- data for NH4-SO4 -----------------------------
  data (b6(I),I=1,13) /           -1.2058223e-2_dp, 1.1043_dp,     4.79018e-5_dp, &
    2.14346e-2_dp, 0.58_dp,       -2.9146e-2_dp,    1.9631e-4_dp,  1.1378_dp,     &
    0.9283_dp,     1.28548e-4_dp,  1.684e-5_dp,     2.6267e-2_dp, -2.6e-4_dp /       ! wt=1


!  data (b6(I),I=1,13)/           -1.2058223e-2_dp,  1.1043_dp,    4.79018e-5_dp, &
!    2.14346e-2_dp,  0.58_dp,     -0.1188_dp,        8.5e-2_dp,    2.10514_dp,    &
!    0.5942_dp,      7.888e-4_dp, -5.503e-4_dp,      5.815e-2_dp, -3.766e-2_dp /     ! wt=5d-5


! ----------- data for NH4-NO3 -----------------------------
  data (b7(I),I=1,13) /-2.3275e-2_dp, 0.15_dp, 1.1634e-4_dp,  1.62e-3_dp, &
    0.43_dp,  & ! 0.107264_dp,   0.221_dp,    -3.8442e-4_dp, -1.3872e-2_dp, 4*0._dp/
    8.78e-2_dp,   0.2753645_dp, -3.349e-4_dp, -1.093e-2_dp,               &
   -4.769e-2_dp,  0.1776_dp,     1.25e-4_dp,   6.9751e-3_dp /

!   2.15438e-2_dp, 0.67073_dp, 2*0._dp, -1.3662e-2_dp, 3.4747e-2_dp, 2*0._dp/


! ----------- data for NH4-Cl ------------------------------
  data (b8(I),I=1,9) /       -6.333e-4_dp,  -3.99546e-4_dp, 0.3155_dp, &
    0.1414_dp, -3.837e-5_dp,  1.08331e-4_dp, 5.2436e-2_dp,  1.6827e-2_dp, 1.19_dp /

!! ----------- data for Na-SO4 ------------------------------              ! <jjb> unused
!  data (b9(I),I=1,17) /           1.63e-3_dp,    2.092e-3_dp,   3.484156e-2_dp, &
!   -1.057e-2_dp,   1.0775_dp,     0.9_dp,        0.8206_dp,    -9.6425e-2_dp,   &
!    2.7492e-3_dp, -9.2838e-4_dp, -6.8268e-4_dp,  4.8126e-5_dp,  7.182e-2_dp,    &
!    0.7586_dp,    -0.30291_dp,    0.2311_dp,     1.7_dp/

  xs(4)=0._dp
  xs(5)=0._dp
  xs(6)=0._dp

! -- Common definition
  dt=(T-T1)/100._dp
  dt2 = dt*dt
  dt3 = dt*dt2
  dt4 = dt2*dt2

  do I=1,NC-1   ! Note the use of NC-1 here. Case NC=3: data are already in Pitzer subroutine
     do J=1,NA
        ! -- Iflag == 1 : H-NO3
        if ( Iflag(I,J) .eq. 1) then
           B0(i,j)=bb(1) +dt*bb(2) +dt2*bb(3) +dt3*bb(4) +dt4*bb(5)
           B1(i,j)=bb(6) +dt*bb(7) +dt2*bb(8) +dt3*bb(9) +dt4*bb(10)
           c0(i,j)=bb(11)+dt*bb(12)+dt2*bb(13)+dt3*bb(14)+dt4*bb(15)
           c1(i,j)=bb(16)+dt*bb(17)+dt2*bb(18)+dt3*bb(19)+dt4*bb(20)

           omega(i,j)=bb(21)

           xs(1)=bb(22)+bb(23)*dt
           xs(2)=bb(24)+bb(25)*dt
           xs(3)=bb(26)+bb(27)*dt
           xs(11)=4.75458e-4_dp - 4.0577e-3_dp*dt
                 !4.633e-4_dp - 4.093e-3_dp*dt

        ! -- Iflag == 2 : H-Cl
        else if (Iflag(I,J) .eq. 2) then
           B0(i,j)=b2(1) +dt*b2(2) +dt2*b2(3) +dt3*b2(4) +dt4*b2(5)
           B1(i,j)=b2(6) +dt*b2(7) +dt2*b2(8) +dt3*b2(9) +dt4*b2(10)
           c0(i,j)=b2(11)+dt*b2(12)+dt2*b2(13)+dt3*b2(14)+dt4*b2(15)
           c1(i,j)=b2(16)+dt*b2(17)+dt2*b2(18)+dt3*b2(19)+dt4*b2(20)

           omega(i,j)=b2(21)


        ! -- Iflag == 3 : H-HSO4  ( and Iflag == 4 : H-SO4 )
        else if (Iflag(I,j).eq.3) then ! NB: Iflag(i,j) == 3 when i=1, j=1.
                                       !    this case also fills arrays(1,2) which is Iflag == 4
           B0(i,j)=b3(1) +dt*b3(2) +dt2*b3(3) +dt3*b3(4) +dt4*b3(5)
           B1(i,j)=b3(6) +dt*b3(7) +dt2*b3(8) +dt3*b3(9) +dt4*b3(10)
           c0(i,j)=b3(11)+dt*b3(12)+dt2*b3(13)+dt3*b3(14)+dt4*b3(15)
           c1(i,j)=b3(16)+dt*b3(17)+dt2*b3(18)+dt3*b3(19)+dt4*b3(20)

           B0(i,j+1)=b3(21)+dt*b3(22)+dt2*b3(23)+dt3*b3(24)+dt4*b3(25)
           B1(i,j+1)=b3(26)+dt*b3(27)+dt2*b3(28)+dt3*b3(29)+dt4*b3(30)
           c0(i,j+1)=b3(31)+dt*b3(32)+dt2*b3(33)+dt3*b3(34)+dt4*b3(35)
           c1(i,j+1)=b3(36)+dt*b3(37)+dt2*b3(38)+dt3*b3(39)+dt4*b3(40)

           omega(I,J)=b3(41)
           omega(I,J+1)=b3(42)

        ! -- Iflag == 5 : NH4-HSO4
        else if ( Iflag(I,J) .eq. 5) then
           B0(I,J)=b5(1)+b5(11+1)*dt+b5(11+2)*dt2
           B1(I,J)=b5(2)+b5(11+3)*dt+b5(11+4)*dt2
           c0(I,J)=b5(3)+b5(11+5)*dt+b5(11+6)*dt2
           c1(I,J)=b5(4)+b5(11+7)*dt+b5(11+8)*dt2

           omega(I,J)=b5(5)

           xs(7) =b5(6) +b5(11+13)*dt+dt2*b5(11+14)
           xs(8) =b5(7) +b5(11+15)*dt+dt2*b5(11+16)
           xs(9) =b5(10)+b5(11+17)*dt+dt2*b5(11+18)
           xs(10)=b5(11)+b5(11+19)*dt+dt2*b5(11+20)

        ! -- Iflag == 6 : NH4-SO4
        else if ( Iflag(I,J) .eq. 6) then
           b0(I,J)=b6(1)+b6(6)*dt+b6(7)*dt2
           b1(I,J)=b6(2)+b6(8)*dt+b6(9)*dt2
           C0(I,J)=b6(3)+b6(10)*dt+b6(11)*dt2
           C1(I,J)=b6(4)+b6(12)*dt+b6(13)*dt2

           omega(I,J)=b6(5)

        ! -- Iflag == 7 : NH4-NO3
        else  if (Iflag(I,j).eq.7) then
           b0(I,J)=b7(1)+b7(6)*dt+b7(10)*dt2
           b1(I,J)=b7(2)+b7(7)*dt+b7(11)*dt2
           c0(I,J)=b7(3)+b7(8)*dt+b7(12)*dt2
           c1(I,J)=b7(4)+b7(9)*dt+b7(13)*dt2

           omega(I,j)=b7(5)

           !xs(1)=bb(22)+bb(23)*dt            ! <jjb> already defined (identical) in H-NO3 case
           !xs(2)=bb(24)+bb(25)*dt
           !xs(3)=bb(26)+bb(27)*dt
           !xs(11)=4.75458e-4_dp - 4.0577e-3_dp * dt
           !      !4.633e-4_dp - 4.093e-3_dp * dt

        ! -- Iflag == 8 : NH4-Cl
        else  if (Iflag(I,j).eq.8) then
           b0(I,J)=b8(1)+b8(2)*dt
           b1(I,J)=b8(3)+b8(4)*dt
           c0(I,J)=b8(5)+b8(6)*dt
           c1(I,J)=b8(7)+b8(8)*dt

           omega(I,j)=b8(9)

!        else  IF (Iflag(I,J).eq.9) then                       ! <jjb> unused. See SR pitzer
!           B0(i,j)=b9( 1)+dt*b9( 2)+dt2*b9( 3)+dt3*b9( 4)
!           B1(i,j)=b9( 5)+dt*b9( 6)+dt2*b9( 7)+dt3*b9( 8)
!           c0(i,j)=b9( 9)+dt*b9(10)+dt2*b9(11)+dt3*b9(12)
!           c1(i,j)=b9(13)+dt*b9(14)+dt2*b9(15)+dt3*b9(16)
!
!           omega(I,j)=b9(17)

        endif
     enddo
  enddo

end subroutine calpar

!-------------------------------------------------------------

function gammann(T,N,NC,NA,mC,mA,ZC,ZA,ZI,I,I2,B0,B1,C0,C1,omega)
!
!      This only a simple version, which only take the unsymmetrical
!      factor e-theta and E-theta ' into account.
!

! 25/01/2017 <jjb> cleaned, properly indented. Later: rewritten without labels, goto
!
! 03/02/2017 <jjb> changed nested do loops order for efficiency (4 places)
!                  za and zc are integers: changed declaration, use "real" to convert kind when necessary
!                  imported ZI from calling subroutine, instead of recalculating it
!
! 03/03/2017 <jjb> added a test to avoid division by zero: ABS(XO4) > TINY(0.)
!                  this case should already be handled in calling subroutine activ, but left here for security
!
! 01/11/2017 <jjb> Fortran90
!                  turned alpha into a parameter, and 273.15 into the parameter T0


  USE precision, ONLY : &
! Imported Parameters:
       dp, tiny_dp

  implicit none

! Function declaration
  real(kind=dp) :: gammann

! Subroutine arguments
! Scalar arguments with intent(in):
  real(kind=dp), intent(in) :: T      ! temperature, in K
  integer,       intent(in) :: N      ! the index of the ion in the C+A list
  integer,       intent(in) :: NA, NC ! number of anion and cation species
! Array arguments with intent(in):
  real(kind=dp), intent(in) :: mC(nc), mA(na)      ! molality of cations and anions
  integer,       intent(in) :: ZC(nc), ZA(na)      ! ionic charge of cations and anions
  real(kind=dp), intent(in) :: ZI                  ! ZI = sum (zc*mc) + sum(za*ma)
  real(kind=dp), intent(in) :: I, I2               ! I = sum (mc*zc**2) + sum(ma*za**2) and I2=sqrt(I)
  real(kind=dp), intent(in) :: B0(nc,na),B1(nc,na) ! Pitzer coefficients
  real(kind=dp), intent(in) :: C0(nc,na),C1(nc,na)
  real(kind=dp), intent(in) :: omega(nc,na)

! External functions:
  real(kind=dp), external :: g, gs

! Local parameters:
  real(kind=dp), parameter :: alpha = 2._dp
  real(kind=dp), parameter :: T0 = 273.15_dp

! Local scalars:
  real(kind=dp) :: Aphi
  real(kind=dp) :: a1, a2, a3
  real(kind=dp) :: E, ED
  real(kind=dp) :: F1, F2, F3, F4, F
  real(kind=dp) :: gam
  real(kind=dp) :: gg, ggs
  real(kind=dp) :: omega1
  real(kind=dp) :: x, xo, xo4, xhx, xhxs
  integer :: icharg, jcharg                 ! arguments passed to SR efunc
  integer :: IA,IA1,IA2, IC,IC1,IC2         ! loop indexes
  integer :: N1
  integer :: Z1, Z2

! Local arrays:
  real(kind=dp) :: B(nc,na),Bs(nc,na)
  real(kind=dp) :: C(nc,na),Cs(nc,na)

!- End of header ---------------------------------------------------------------


! --  Calculate B, Bs, C, Cs
  x=I2*alpha
  gg=g(x)
  ggs=gs(x)

  do IA = 1, Na
     do IC = 1, NC
        B(IC,IA) = b0(IC,Ia) + gg* b1(IC,Ia)
        BS(IC,IA) = ggs*b1(IC,Ia)/I

        omega1=OMEGA(IC,IA)
        xo=omega1*I2
        xo4 = xo*xo*xo*xo

        if(abs(XO4) > tiny_dp) then
           !xhx=1._dp/xo**4*(6._dp-exp(-xo)*(6._dp+6._dp*xo+3._dp*xo**2+xo**3))
           xhx=1._dp/xo4*(6._dp-exp(-xo)*(6._dp+6._dp*xo+3._dp*xo**2+xo**3))
        else
           print*,'Warning: in SR gammann, xo4 too small'
           xhx=0._dp
        endif
        xhxs=exp(-xo)/2._dp - 2._dp*xhx

        C(IC,Ia)  = (C0(IC,Ia) + 4._dp * C1(IC,Ia) * xhx)
        Cs(IC,Ia) =  C1(IC,Ia) / I * xhxs
     enddo
  enddo

  Aphi=.377_dp + 4.684e-4_dp * (T-T0) + 3.74e-6_dp*(T-T0)**2

  F1= -Aphi*(I2/(1._dp + 1.2_dp * I2) + 2._dp / 1.2_dp*log(1._dp + 1.2_dp*I2))

  F2=0._dp
  do IA = 1, NA
     do IC = 1, NC
        F2=F2+mc(IC)*Ma(Ia)*(BS(IC,Ia)+2._dp*ZI*CS(IC,IA) )
     enddo
  enddo

  ICHARG=1
  JCHARG=2
  call EFUNC(ICHARG,JCHARG,Aphi,I,E,ED)

  F3= 0.
  do IC1 = 1, NC
     do IC2 = IC1+1, NC
        z1=ZC(IC1)
        z2=ZC(IC2)
        if(Z1.ne.Z2) F3=F3+ED*MC(IC1)*MC(IC2)
     enddo
  enddo

  F4=0.
  do IA1 = 1, NA
     do IA2 = IA1+1, NA
        Z1=ZA(IA1)
        Z2=ZA(IA2)
        if(Z1.ne.Z2) F4=F4+ED*MA(IA1)*MA(IA2)
     enddo
  enddo

  F=F1+F2+F3+F4

  ! -- Cations
  if (N.le.NC) then
     a1=real(ZC(N)**2,dp) *F

     a2=0.
     do Ia=1,NA
        a2=a2+MA(IA)*( 2*B(N,IA)+ZI*C(N,IA) )
     enddo

     a3=0.
     do IA=1,NA
        do IC=1,NC
           a3=a3+MA(IA)*MC(IC)*C(IC,IA)
        enddo
     enddo
     a3=real(ZC(N),dp)*a3

     F3=0.
     do IC = 1, NC
        Z1=ZC(N)
        Z2=ZC(IC)
        if(Z1.ne.Z2) F3=F3+E*MC(IC)
     enddo

     gam=a1+a2+a3 +F3

  ! -- Anions
  else
     N1=N-NC

     a1=real(ZA(N1)**2,dp) *F

     a2=0.
     do IC=1,NC
        a2=a2+MC(IC)*( 2._dp*B(IC,N1)+ZI*C(IC,N1) )
     enddo

     a3=0.
     do IA=1,NA
        do IC=1,NC
           a3=a3+MA(IA)*MC(IC)*C(IC,IA)
        enddo
     enddo
     a3=real(ZA(N1),dp)*a3

     F4=0.
     do IA = 1, NA
        Z1=ZA(N1)
        Z2=ZA(IA)
        if(Z1.ne.Z2) F4=F4+E*MA(IA)
     enddo

     gam=a1+a2+a3+F4

  endif

  gammann= exp(gam)

end function gammann

!-------------------------------------------------------------

function g(x)

!    01/2017 <jjb> implicit none
! 01/11/2017 <jjb> Fortran90

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none
  real(kind=dp) :: g,x

  g=2._dp*(1._dp-(1._dp+x)*exp(-x))/x**2
end function g

!-------------------------------------------------------------

function gs(x)

!    01/2017 <jjb> implicit none
! 01/11/2017 <jjb> Fortran90

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none
  real(kind=dp) :: gs, x

  gs=2._dp*(-1._dp+(1._dp+x+x**2/2._dp)*exp(-x) )/x**2
end function gs

!-------------------------------------------------------------

subroutine EFUNC(ICHARG,JCHARG,APHI,XI,E,ED)

! 03/02/2017 <jjb> declared all variables + organised header
!                  removed internal variables XE, XED which were identical to E, ED
!                  renamed A -> APHI to be consistent with calling name
!                  in calling SRs, renamed J1 -> ICHARG and J2 -> JCHARG
!
! 01/11/2017 <jjb> Fortran90

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  integer,       intent(IN) :: ICHARG, JCHARG
  real(kind=dp), intent(IN) :: APHI, XI
! Scalar arguments with intent(out):
  real(kind=dp), intent(OUT) :: E, ED

! Local scalars:
  integer :: I
  real(kind=dp) :: DUM
! Local arrays:
  real(kind=dp) :: X(3)
  real(kind=dp) :: J0(3), J1(3)
!- End of header ---------------------------------------------------------------

  if((ICHARG.eq.JCHARG) .or. (XI .le. 1.e-30_dp) )then
     E=0._dp
     ED=0._dp
  else
     X(1)=6._dp*ICHARG*JCHARG*sqrt(XI)*APHI
     X(2)=6._dp*ICHARG*ICHARG*sqrt(XI)*APHI
     X(3)=6._dp*JCHARG*JCHARG*sqrt(XI)*APHI

     do I=1,3
        DUM=-1.2e-2_dp*X(I)**.528_dp
        J0(I)=X(I)/(4._dp + 4.581_dp * X(I)**(-.7238_dp)*exp(DUM))
        J1(I)=(4._dp + 4.581_dp * X(I)**(-.7238_dp)*exp(DUM)*(1.7238_dp-DUM*.528_dp)) &
             /(4._dp + 4.581_dp * X(I)**(-.7238_dp)*exp(DUM))**2
     enddo

     E =ICHARG*JCHARG/4._dp/XI*(J0(1)-.5_dp*J0(2)-.5_dp*J0(3))
     ED=ICHARG*JCHARG/8._dp/XI**2*(X(1)*J1(1)-.5_dp*X(2)*J1(2)-.5_dp*X(3)*J1(3))-E/XI
  endif

end subroutine EFUNC

!------------------------------------------------------------------

function GAMMASN(T,NC,NA,MC,MA,ZC,ZA,ZI,I,I2,B0,B1,C0,C1,OMEGA)

!      GAMMASN === WATER ACTIVITY = PH2O/PH2O(PURE WATER)
!      ALL ARGUMENTS ARE INPUTS
!      C1 == 0.
!      OMEGA =1.


! 03/02/2017 <jjb> changed nested do loops order for efficiency (2 places)
!                  removed unreferenced variable xmh0
!                  za and zc are integers: changed declaration, use "real" to convert kind when necessary
!                  imported ZI, I and I2 from calling subroutine, instead of recalculating it

! 01/11/2017 <jjb> Fortran90
!                  turned alpha into a parameter, and 273.15 into the parameter T0
!                  use module for M_wat instead of 18./1000.

  USE constants, ONLY : &
! Imported Parameters:
       M_wat               ! Water molar mass [kg/mol]

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Function declaration:
  real(kind=dp) :: gammasn

! Subroutine arguments
! Scalar arguments with intent(in):
  real(kind=dp), intent(in) :: T      ! temperature, in K
  integer,       intent(in) :: NA, NC ! number of anion and cation species
! Array arguments with intent(in):
  real(kind=dp), intent(in) :: MC(nc), MA(na)      ! molality of cations and anions
  integer,       intent(in) :: ZC(nc), ZA(na)      ! ionic charge of cations and anions
  real(kind=dp), intent(in) :: ZI                  ! ZI = sum (mc*zc) + sum(ma*za)
  real(kind=dp), intent(in) :: I, I2               ! I = sum (mc*zc**2) + sum(ma*za**2) and I2=sqrt(I)
  real(kind=dp), intent(in) :: B0(nc,na),B1(nc,na) ! Pitzer coefficients
  real(kind=dp), intent(in) :: C0(nc,na),C1(nc,na)
  real(kind=dp), intent(in) :: OMEGA(nc,na)


! Local parameters:
  real(kind=dp), parameter :: ALPHA = 2._dp
  real(kind=dp), parameter :: T0 = 273.15_dp

! Local scalars:
  real(kind=dp) :: APHI, AS
  real(kind=dp) :: E, ED
  real(kind=dp) :: F3, F4, FPHI1
  real(kind=dp) :: OMEGA1
  real(kind=dp) :: PP, PHI, PHIX
  real(kind=dp) :: X, XO, XMI, XS
  integer :: icharg, jcharg            ! arguments passed to SR efunc
  integer :: IA,IA1,IA2, IC,IC1,IC2    ! loop indexes
  integer :: Z1, Z2

! Local arrays:
  real(kind=dp) :: BPHI(NC,NA),C(NC,NA)

!- End of header ---------------------------------------------------------------


! --  CALCULATE BPHI, C
  X=I2*ALPHA
  do IA = 1, NA
     do IC = 1, NC
        BPHI(IC,IA) = B0(IC,IA) + exp(-X)*B1(IC,IA)
        OMEGA1 = OMEGA(IC,IA)                       ! see comment in fct header
        XO = OMEGA1*I2
        C(IC,IA) = C0(IC,IA) + C1(IC,IA)*exp(-XO)   ! see comment in fct header
     enddo
  enddo

  APHI=.377_dp + 4.684e-4_dp*(T-T0)+3.74e-6_dp*(T-T0)**2

! --  CALCULATE SUM MI
  XMI=0._dp
  do IC=1,NC
     XMI=XMI+MC(IC)
  enddo
  do IA=1,NA
     XMI=XMI+MA(IA)
  enddo

  FPHI1 = -APHI*I**1.5_dp / (1._dp+1.2_dp*I2) ! 1.5=3./2.
  XS = 0._dp
  do IA=1,NA
     do IC=1,NC
        XS = XS + MA(IA)*MC(IC)*(ZI*C(IC,IA)+BPHI(IC,IA))
     enddo
  enddo

  ICHARG=1
  JCHARG=2
  call EFUNC(ICHARG,JCHARG,APHI,I,E,ED)

  PP=E+I*ED

  F3 = 0._dp
  do IC1 = 1, NC
     do IC2 = IC1+1, NC
        Z1=ZC(IC1)
        Z2=ZC(IC2)
        if(Z1.ne.Z2) F3 = F3 + PP*MC(IC1)*MC(IC2)
     enddo
  enddo

  F4 = 0._dp
  do IA1 = 1, NA
     do IA2 = IA1+1, NA
        Z1 = ZA(IA1)
        Z2 = ZA(IA2)
        if(Z1.ne.Z2) F4 = F4 + PP*MA(IA1)*MA(IA2)
     enddo
  enddo

  PHIX =  FPHI1 + XS  + F3 + F4
  PHI = 1._dp + PHIX*2._dp/XMI
  AS = (-PHI) * M_wat * XMI
  GAMMASN = exp(AS)

end function GAMMASN
