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


module data_surface

! Description :
! -----------
  ! Global declaration of model variables for surface and canopy module

! Author :
! ------
  ! Werner Schneider (original version)
  ! Josue Bock (rewritten for this version of Mistra)

! Declarations:
! ------------
! Modules used:
  use precision, only : &
       dp                    ! double precision kind

  implicit none

  public ! confirm this property with an explicit statement
  save   ! not mandatory if called from main program


! Prandtl layer, Clarke functions, frictional velocity
!-----------------------------------------------------
  real (kind=dp) :: &
       xzpdl(18), & ! tabled values of zpdl for interpolation (claf)        (1)
       xzpdz0(7), & ! tabled values of zpdz0 for interpolation (claf)       (1)
       fu(18,7), &  ! tabled values of clarke function for momentum         (1)
       ft(18,7), &  ! tabled values of clarke function for temp., moisture  (1)
       gclu, &      ! Clarke function for momentum                          (1)
       gclt, &      ! Clarke function for temperature, humidity             (1)
       ustern, &    ! frictional velocity                                 (m/s)
       z0           ! roughness length of the surface                       (m)

! surface temperature for water surface
!--------------------------------------
  real (kind=dp) :: &
       tw           ! water surface temperature                          (K)


  !-------------------------------------
  !  data_soil
  !-------------------------------------
! soil constants for sandy loam
  real (kind=dp), parameter :: aks = 3.41e-5_dp    ! hydraulic conductivity for saturated soil     (m/s)
  ! (renamed hcs in the latest version)
  real (kind=dp), parameter :: anu0 = 43.415524_dp ! reference for thermal conductivity
  real (kind=dp), parameter :: bs = 4.9_dp         ! exponent (calculation of moisture potential)    (1)
  ! (renamed b in the latest version)
  real (kind=dp), parameter :: bs0 = 2.128043_dp   ! reference for exponent bs
  real (kind=dp), parameter :: ebc = .0742724_dp   ! reference value for soil moisture
  real (kind=dp), parameter :: ebs = .435_dp       ! volumetric porosity of the soil           (m^3/m^3)
  real (kind=dp), parameter :: psis = -.218_dp     ! moist potential for saturated soil              (m)

end module data_surface
