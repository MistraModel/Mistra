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


      subroutine bud_s_gas (dtg,k)
!
! Description:
!    gas phase chemical reactions budget
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
      include 'gas_Parameters.h' !additional common blocks and other definitions
      include 'gas_Global.h' !additional common blocks and other definitions

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


! rates (mol/(m^3*s)

      bgs(1,  1,k)=RCONST( 67)*C(ind_SO2)*C(ind_OH)
! DMS related rxns
      bgs(1, 75,k)=RCONST( 68)*C(ind_DMS)*C(ind_OH)
      bgs(1, 76,k)=RCONST( 69)*C(ind_DMS)*C(ind_OH)
      bgs(1, 77,k)=RCONST( 70)*C(ind_DMS)*C(ind_NO3)
      bgs(1, 78,k)=RCONST( 71)*C(ind_DMS)*C(ind_Cl)
      bgs(1, 79,k)=RCONST( 72)*C(ind_DMS)*C(ind_Br)
      bgs(1, 80,k)=RCONST( 73)*C(ind_DMS)*C(ind_BrO)
      bgs(1,101,k)=RCONST( 78)*C(ind_CH3SO)*C(ind_NO2)
      bgs(1,102,k)=RCONST( 80)*C(ind_CH3SO2)
      bgs(1,103,k)=RCONST( 83)*C(ind_CH3SO3)*C(ind_HO2)
      bgs(1,104,k)=RCONST( 87)*C(ind_CH3SO2H)*C(ind_O3)
      bgs(1,121,k)=RCONST( 86)*C(ind_CH3SO2H)*C(ind_OH)
      bgs(1,122,k)=RCONST( 84)*C(ind_CH3SO3)

! accumulated rates (mlc/cm^3)
      do i=1,1
         bgs(2,i,k)=bgs(2,i,k)+dtg*bgs(1,i,k)
      enddo
      do i=75,80
         bgs(2,i,k)=bgs(2,i,k)+dtg*bgs(1,i,k)
      enddo
      do i=101,104
         bgs(2,i,k)=bgs(2,i,k)+dtg*bgs(1,i,k)
      enddo
      do i=121,122
         bgs(2,i,k)=bgs(2,i,k)+dtg*bgs(1,i,k)
      enddo
      return
      end subroutine bud_s_gas
