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


module precision

! Define kind of variables used throughout the model


! Author:
! -------
!     Josue Bock


! Modifications :
! -------------
  ! 28-Oct-2017  Josue Bock  introduce tiny_dp to replace == 0. tests by < tiny_dp

implicit none

save

integer, parameter :: dp = kind(1.d0)
real(kind=dp), parameter :: tiny_dp = tiny(0._dp)

end module precision
