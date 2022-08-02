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


module global_params

! Description :
! -----------
!   Declaration of global parameters of Mistra


! Author :
! ------  
!   Josue Bock


! Modifications :
! -------------
  ! 12-Aug-2016  Josue Bock  First version (introduced in my version u1)
  ! 31-Oct-2017  Josue Bock  Fortran90, public, added comments

  implicit none

  include 'tot_Parameters.h'

  public
  save
! param_met
  ! number of layers, constant height
  integer, parameter :: nf = 100

  ! number of layers, constant height + exponentially increasing height
  integer, parameter :: n = nf+50

  integer, parameter :: nm = n-1

! param_mic
  ! number of dry radius classes in the microphysical 2D grid
  integer, parameter :: nka = 70

  ! number of total radius classes in the microphysical 2D grid
  integer, parameter :: nkt = 70

! param_rad
  ! number of layers including the extra ones for radiation calculation
  !   see SR initr (in radinit.f) for height definition
  !   note that the infinitesimally thin layer is ignored, thus only
  !   (n-1) layers are taken from the main model to the radiative code
  integer, parameter :: nrlay = (n -1) + 11

  integer, parameter :: nrlev = nrlay + 1

  ! Total number of wavelength bands in the radiative code
  integer, parameter :: mb = 18

  ! Wavelength bands for the solar region
  integer, parameter :: mbs = 6

  ! Wavelength bands for the infrared region
  integer, parameter :: mbir = 12

  ! Number of aerosol types whose radiative properties (qabs, qext and asym)
  !   have been tabulated with Mie calculations. See SR intrad.
  !   1=urban, 2=rural, 3=ocean
  integer, parameter :: jptaerrad = 3

  ! Number of tabulated classes for the parameters related to liquid water
  !   in the radiation code (see SR water)
  integer, parameter :: ncw = 8

! parameters soil grid
  integer, parameter :: nb = 20

  integer, parameter :: nbm = nb-1

! param_che
  ! Number of non-radical gases
  !integer, parameter :: j1 = 96
  integer, parameter :: j1_fake = 96

  ! Number of radical gas species
  integer, parameter :: j3 = 25

  integer, parameter :: j2 = j1_fake + j3

  integer, parameter :: j6 = 55

  ! Number of chemical bins
  integer, parameter :: nkc = 4

  ! Number of photolysis reactions
  !  (if updated here, *.def files for KPP must also be updated for the same variable)
  integer, parameter :: nphrxn = 47

  ! Number of levels to compute budget
  integer, parameter :: nlev = 15

  ! Number of chemical reactions in the tot mechanism
  !   (improvement: get this value from KPP header files)
  !integer, parameter :: nrxn = 1627
  integer, parameter :: nrxn = NREACT

  ! Maximum level to compute prognostic chemistry in the aerosols
  integer, parameter :: nmax_chem_aer = nf

end module global_params
