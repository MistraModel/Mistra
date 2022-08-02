! ======================================================================= !
! Copyright 2015-2017 the Authors
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
MODULE cdf_var_gas

! contains the old common /cdf_var_gas/
! data exchanged between open / write / close gas subroutines

! Author:
! -------
!     Josue Bock

! History:
! --------
!     20/11/2016 -- Josue Bock -- First version

USE gas_common, ONLY : &
&     j1, &
&     j5

IMPLICIT NONE
SAVE

INTEGER :: idgasfile, igascount
INTEGER, ALLOCATABLE :: idvar_gas(:)

END MODULE cdf_var_gas
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
