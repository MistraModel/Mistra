!
! Copyright 2015-2017 the Authors
!
! Licensed under the EUPL, Version 1.1 only (the "Licence");
!
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at:
!   https://joinup.ec.europa.eu/software/page/eupl
!
! Unless required by applicable law or agreed to in writing,
! software distributed under the Licence is distributed on an
! "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
! either express or implied.
!
! See the Licence for the specific language governing permissions
! and limitations under the Licence.


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
