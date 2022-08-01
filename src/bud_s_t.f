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


      subroutine bud_s_tot (dtg,k)
!
! Description:
!    chemical reactions budget
!
! Authors:
! --------
!     Roland von Glasow
!

! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      08/2016   Use module for parameters                <Josue Bock>
!                    Comments / header / implicit none
!                    Removed three unused argument
!                    Cleaning
!
! 1.0       ?        Original code.                          <Roland von Glasow>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:
      USE global_params, ONLY :
! Imported Parameters:
     &     n

      implicit none

! Include statements:
      include 'tot_Parameters.h' !additional common blocks and other definitions
      include 'tot_Global.h' !additional common blocks and other definitions

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision dtg
      integer k

! Local scalars:
      integer i

! Common blocks:
      common /budgs/ bgs(2,122,n)
      double precision bgs

!- End of header ---------------------------------------------------------------


c rates (mol/(m^3*s)

      bgs(1,  1,k)=RCONST( 67)*C(ind_SO2)*C(ind_OH)
      bgs(1,  2,k)=RCONST(346)*C(ind_H2SO4)-RCONST(347)*C(ind_SO4l1)
      bgs(1,  3,k)=RCONST(595)*C(ind_H2SO4)-RCONST(596)*C(ind_SO4l2)
      bgs(1,  4,k)=RCONST(195)*C(ind_HSO3ml1)*C(ind_O3l1)
      bgs(1,  5,k)=RCONST(196)*C(ind_SO32ml1)*C(ind_O3l1)
      bgs(1,  6,k)=RCONST(197)*C(ind_HSO3ml1)*C(ind_OHl1)
      bgs(1,  7,k)=RCONST(198)*C(ind_SO32ml1)*C(ind_OHl1)
      bgs(1,  8,k)=RCONST(199)*C(ind_HSO3ml1)*C(ind_HO2l1)
      bgs(1,  9,k)=RCONST(200)*C(ind_HSO3ml1)*C(ind_O2ml1)
      bgs(1, 10,k)=RCONST(201)*C(ind_HSO3ml1)*C(ind_H2O2l1)
      bgs(1, 11,k)=RCONST(205)*C(ind_HSO3ml1)*C(ind_HNO4l1)
      bgs(1, 12,k)=RCONST(206)*C(ind_HSO3ml1)*C(ind_CH3OOHl1)
      bgs(1, 13,k)=RCONST(207)*C(ind_SO32ml1)*C(ind_CH3OOHl1)
      bgs(1, 14,k)=RCONST(262)*C(ind_HOCll1)*C(ind_SO32ml1)
      bgs(1, 15,k)=RCONST(263)*C(ind_HOCll1)*C(ind_HSO3ml1)
      bgs(1, 16,k)=RCONST(289)*C(ind_BrOml1)*C(ind_SO32ml1)
      bgs(1, 17,k)=RCONST(293)*C(ind_HOBrl1)*C(ind_SO32ml1)
      bgs(1, 18,k)=RCONST(294)*C(ind_HOBrl1)*C(ind_HSO3ml1)
      bgs(1, 19,k)=RCONST(344)*C(ind_SO2)
      bgs(1, 20,k)=RCONST(345)*C(ind_SO2l1)
      bgs(1, 21,k)=RCONST(444)*C(ind_HSO3ml2)*C(ind_O3l2)
      bgs(1, 22,k)=RCONST(445)*C(ind_SO32ml2)*C(ind_O3l2)
      bgs(1, 23,k)=RCONST(446)*C(ind_HSO3ml2)*C(ind_OHl2)
      bgs(1, 24,k)=RCONST(447)*C(ind_SO32ml2)*C(ind_OHl2)
      bgs(1, 25,k)=RCONST(448)*C(ind_HSO3ml2)*C(ind_HO2l2)
      bgs(1, 26,k)=RCONST(449)*C(ind_HSO3ml2)*C(ind_O2ml2)
      bgs(1, 27,k)=RCONST(450)*C(ind_HSO3ml2)*C(ind_H2O2l2)
      bgs(1, 28,k)=RCONST(454)*C(ind_HSO3ml2)*C(ind_HNO4l2)
      bgs(1, 29,k)=RCONST(455)*C(ind_HSO3ml2)*C(ind_CH3OOHl2)
      bgs(1, 30,k)=RCONST(456)*C(ind_SO32ml2)*C(ind_CH3OOHl2)
      bgs(1, 31,k)=RCONST(511)*C(ind_HOCll2)*C(ind_SO32ml2)
      bgs(1, 32,k)=RCONST(512)*C(ind_HOCll2)*C(ind_HSO3ml2)
      bgs(1, 33,k)=RCONST(538)*C(ind_BrOml2)*C(ind_SO32ml2)
      bgs(1, 34,k)=RCONST(542)*C(ind_HOBrl2)*C(ind_SO32ml2)
      bgs(1, 35,k)=RCONST(543)*C(ind_HOBrl2)*C(ind_HSO3ml2)
      bgs(1, 36,k)=RCONST(593)*C(ind_SO2)
      bgs(1, 37,k)=RCONST(594)*C(ind_SO2l2)
      bgs(1, 38,k)=RCONST(862)*C(ind_H2SO4)-RCONST(863)*C(ind_SO4l3)
      bgs(1, 39,k)=RCONST(1111)*C(ind_H2SO4)-RCONST(1112)*C(ind_SO4l4)
      bgs(1, 40,k)=RCONST(711)*C(ind_HSO3ml3)*C(ind_O3l3)
      bgs(1, 41,k)=RCONST(712)*C(ind_SO32ml3)*C(ind_O3l3)
      bgs(1, 42,k)=RCONST(713)*C(ind_HSO3ml3)*C(ind_OHl3)
      bgs(1, 43,k)=RCONST(714)*C(ind_SO32ml3)*C(ind_OHl3)
      bgs(1, 44,k)=RCONST(715)*C(ind_HSO3ml3)*C(ind_HO2l3)
      bgs(1, 45,k)=RCONST(716)*C(ind_HSO3ml3)*C(ind_O2ml3)
      bgs(1, 46,k)=RCONST(717)*C(ind_HSO3ml3)*C(ind_H2O2l3)
      bgs(1, 47,k)=RCONST(721)*C(ind_HSO3ml3)*C(ind_HNO4l3)
      bgs(1, 48,k)=RCONST(722)*C(ind_HSO3ml3)*C(ind_CH3OOHl3)
      bgs(1, 49,k)=RCONST(723)*C(ind_SO32ml3)*C(ind_CH3OOHl3)
      bgs(1, 50,k)=RCONST(778)*C(ind_HOCll3)*C(ind_SO32ml3)
      bgs(1, 51,k)=RCONST(779)*C(ind_HOCll3)*C(ind_HSO3ml3)
      bgs(1, 52,k)=RCONST(805)*C(ind_BrOml3)*C(ind_SO32ml3)
      bgs(1, 53,k)=RCONST(809)*C(ind_HOBrl3)*C(ind_SO32ml3)
      bgs(1, 54,k)=RCONST(810)*C(ind_HOBrl3)*C(ind_HSO3ml3)
      bgs(1, 55,k)=RCONST(860)*C(ind_SO2)
      bgs(1, 56,k)=RCONST(861)*C(ind_SO2l3)
      bgs(1, 57,k)=RCONST(960)*C(ind_HSO3ml4)*C(ind_O3l4)
      bgs(1, 58,k)=RCONST(961)*C(ind_SO32ml4)*C(ind_O3l4)
      bgs(1, 59,k)=RCONST(962)*C(ind_HSO3ml4)*C(ind_OHl4)
      bgs(1, 60,k)=RCONST(963)*C(ind_SO32ml4)*C(ind_OHl4)
      bgs(1, 61,k)=RCONST(964)*C(ind_HSO3ml4)*C(ind_HO2l4)
      bgs(1, 62,k)=RCONST(965)*C(ind_HSO3ml4)*C(ind_O2ml4)
      bgs(1, 63,k)=RCONST(966)*C(ind_HSO3ml4)*C(ind_H2O2l4)
      bgs(1, 64,k)=RCONST(970)*C(ind_HSO3ml4)*C(ind_HNO4l4)
      bgs(1, 65,k)=RCONST(971)*C(ind_HSO3ml4)*C(ind_CH3OOHl4)
      bgs(1, 66,k)=RCONST(972)*C(ind_SO32ml4)*C(ind_CH3OOHl4)
      bgs(1, 67,k)=RCONST(1026)*C(ind_HOCll4)*C(ind_SO32ml4)
      bgs(1, 68,k)=RCONST(1028)*C(ind_HOCll4)*C(ind_HSO3ml4)
      bgs(1, 69,k)=RCONST(1054)*C(ind_BrOml4)*C(ind_SO32ml4)
      bgs(1, 70,k)=RCONST(1058)*C(ind_HOBrl4)*C(ind_SO32ml4)
      bgs(1, 71,k)=RCONST(1059)*C(ind_HOBrl4)*C(ind_HSO3ml4)
      bgs(1, 72,k)=RCONST(1109)*C(ind_SO2)
      bgs(1, 73,k)=RCONST(1110)*C(ind_SO2l4)


c DMS related rxns
      bgs(1, 75,k)=RCONST( 68)*C(ind_DMS)*C(ind_OH)
      bgs(1, 76,k)=RCONST( 69)*C(ind_DMS)*C(ind_OH)
      bgs(1, 77,k)=RCONST( 70)*C(ind_DMS)*C(ind_NO3)
      bgs(1, 78,k)=RCONST( 71)*C(ind_DMS)*C(ind_Cl)
      bgs(1, 79,k)=RCONST( 72)*C(ind_DMS)*C(ind_Br)
      bgs(1, 80,k)=RCONST( 73)*C(ind_DMS)*C(ind_BrO)
      bgs(1, 81,k)=RCONST(227)*C(ind_DMSl1)*C(ind_O3l1)
      bgs(1, 82,k)=RCONST(228)*C(ind_DMSl1)*C(ind_OHl1)
      bgs(1, 83,k)=RCONST(229)*C(ind_DMSOl1)*C(ind_OHl1)
      bgs(1, 84,k)=RCONST(230)*C(ind_CH3SO2ml1)*C(ind_OHl1)
      bgs(1, 85,k)=RCONST(476)*C(ind_DMSl2)*C(ind_O3l2)
      bgs(1, 86,k)=RCONST(477)*C(ind_DMSl2)*C(ind_OHl2)
      bgs(1, 87,k)=RCONST(478)*C(ind_DMSOl2)*C(ind_OHl2)
      bgs(1, 88,k)=RCONST(479)*C(ind_CH3SO2ml2)*C(ind_OHl2)
      bgs(1, 89,k)=RCONST(743)*C(ind_DMSl3)*C(ind_O3l3)
      bgs(1, 90,k)=RCONST(744)*C(ind_DMSl3)*C(ind_OHl3)
      bgs(1, 91,k)=RCONST(745)*C(ind_DMSOl3)*C(ind_OHl3)
      bgs(1, 92,k)=RCONST(746)*C(ind_CH3SO2ml3)*C(ind_OHl3)
      bgs(1, 93,k)=RCONST(992)*C(ind_DMSl4)*C(ind_O3l4)
      bgs(1, 94,k)=RCONST(993)*C(ind_DMSl4)*C(ind_OHl4)
      bgs(1, 95,k)=RCONST(994)*C(ind_DMSOl4)*C(ind_OHl4)
      bgs(1, 96,k)=RCONST(995)*C(ind_CH3SO2ml4)*C(ind_OHl4)
      bgs(1, 97,k)=
     &     RCONST(356)*C(ind_CH3SO3H)-RCONST(357)*C(ind_CH3SO3ml1)
      bgs(1, 98,k)=
     &     RCONST(605)*C(ind_CH3SO3H)-RCONST(606)*C(ind_CH3SO3ml2)
      bgs(1, 99,k)=
     &     RCONST(872)*C(ind_CH3SO3H)-RCONST(873)*C(ind_CH3SO3ml3)
      bgs(1,100,k)=
     &     RCONST(1121)*C(ind_CH3SO3H)-RCONST(1122)*C(ind_CH3SO3ml4)
      bgs(1,101,k)=RCONST( 78)*C(ind_CH3SO)*C(ind_NO2)
      bgs(1,102,k)=RCONST( 80)*C(ind_CH3SO2)
      bgs(1,103,k)=RCONST( 83)*C(ind_CH3SO3)*C(ind_HO2)
      bgs(1,104,k)=RCONST( 87)*C(ind_CH3SO2H)*C(ind_O3)
      bgs(1,105,k)=RCONST(257)*C(ind_Cl2ml1)*C(ind_DMSl1)
      bgs(1,106,k)=RCONST(284)*C(ind_Br2ml1)*C(ind_DMSl1)
      bgs(1,107,k)=RCONST(506)*C(ind_Cl2ml2)*C(ind_DMSl2)
      bgs(1,108,k)=RCONST(533)*C(ind_Br2ml2)*C(ind_DMSl2)
      bgs(1,109,k)=RCONST(773)*C(ind_Cl2ml3)*C(ind_DMSl3)
      bgs(1,110,k)=RCONST(800)*C(ind_Br2ml3)*C(ind_DMSl3)
      bgs(1,111,k)=RCONST(1022)*C(ind_Cl2ml4)*C(ind_DMSl4)
      bgs(1,112,k)=RCONST(1049)*C(ind_Br2ml4)*C(ind_DMSl4)
      bgs(1,113,k)=RCONST(231)*C(ind_CH3SO3ml1)*C(ind_OHl1)
      bgs(1,114,k)=RCONST(480)*C(ind_CH3SO3ml2)*C(ind_OHl2)
      bgs(1,115,k)=RCONST(747)*C(ind_CH3SO3ml3)*C(ind_OHl3)
      bgs(1,116,k)=RCONST(996)*C(ind_CH3SO3ml4)*C(ind_OHl4)
      bgs(1,117,k)=
     &     RCONST(354)*C(ind_CH3SO2H)-RCONST(355)*C(ind_CH3SO2ml1)
      bgs(1,118,k)=
     &     RCONST(603)*C(ind_CH3SO2H)-RCONST(604)*C(ind_CH3SO2ml2)
      bgs(1,119,k)=
     &     RCONST(870)*C(ind_CH3SO2H)-RCONST(871)*C(ind_CH3SO2ml3)
      bgs(1,120,k)=
     &     RCONST(1119)*C(ind_CH3SO2H)-RCONST(1120)*C(ind_CH3SO2ml4)
      bgs(1,121,k)=RCONST( 86)*C(ind_CH3SO2H)*C(ind_OH)
      bgs(1,122,k)=RCONST( 84)*C(ind_CH3SO3)
c accumulated rates (mlc/cm^3)
      do i=1+0,122
         bgs(2,i,k)=bgs(2,i,k)+dtg*bgs(1,i,k)
      enddo
      return
      end
