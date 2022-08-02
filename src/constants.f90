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


module constants

! Declare and initialise all physical and chemical constants used in Mistra


! Author:
! -------
!     Josue Bock

USE precision, only: &
  &   dp                   ! kind double precision real

implicit none

save

! Avogadro constant
real(kind=dp), parameter :: Avogadro = 6.022140857e+23_dp        ! [1/mol]

! Calorie:
!     amount of energy required to warm one gram of air-free water from 14.5 to 15.5 C at standard atmospheric pressure
real(kind=dp), parameter :: cal15 = 4.1855_dp                    ! [J]

! Conversion factor 1:
!     conv1 = Avogadro / 10^6
!     multiply by conv1 to convert cm^3(air)/mlc --> m^3(air)/mol
real(kind=dp), parameter :: conv1 = Avogadro * 1.e-6_dp          ! [m3/cm3/mol]

! Air molar mass
real(kind=dp), parameter :: M_air = 28.96546e-3_dp               ! [kg/mol]

! Water molar mass
real(kind=dp), parameter :: M_wat = 18.01528e-3_dp               ! [kg/mol]

! Pi
real(kind=dp), parameter :: pi = 3.1415926535897932_dp

! Rad(1 deg)
real(kind=dp), parameter :: rad = pi / 180._dp

! Molar gas constant
!     ref: Gavioso et al. (2015), Metrologia 52 (S274-S304) doi: 10.1088/0026-1394/52/5/S274
real(kind=dp), parameter :: gas_const = 8.3144743_dp             ! [J/K/mol]

! Specific gas constant of dry air
!     r0 = R / M_air ~= 287.05
real(kind=dp), parameter :: r0 = gas_const / M_air               ! [J/(kg.K)]

! Specific gas constant of water vapour
!     r1 = R / M_water ~= 461.51
real(kind=dp), parameter :: r1 = gas_const / M_wat               ! [J/(kg.K)]

! Water density
!     This could be improved by a parameterisation as a function of T and P
real(kind=dp), parameter :: rhow = 1000.0_dp                     ! [kg/m3]

! Aerosol density
real(kind=dp), parameter :: rho3 = 2000.0_dp                     ! [kg/m3]

! Gravitational acceleration
real(kind=dp), parameter :: g = 9.80665_dp                       ! [m/s2]

! Specific heat capacity of dry air at constant pressure
!     ref: Seinfeld & Pandis, 2nd Ed., Table A.7 p. 1178
real(kind=dp), parameter :: cp = 1005.0_dp                       ! [J/kg/K]

! von Karman constant
real(kind=dp), parameter :: kappa = 0.4_dp                       ! [1]



end module constants
