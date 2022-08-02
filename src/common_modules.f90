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


! In this file are all modules which replace old common block

module gas_common

! Description :
! -----------
!    gas_common holds the variables related to gas phase species

! Author :
! ------
!    Josue Bock

! Modifications :
! -------------
!  26-Sep-2020   Josue Bock   Fortran 90 conversion (straightforward)

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  save

! GASES (non radical)
  ! Numbers of gas species
  integer :: j1                             ! Number of non radical gases actually used
  integer :: j1_br, j1_cl, j1_iod, j1_halo  ! Number of brominated / chlorinated / ...  gases

  ! Indexes
  integer, allocatable :: ind_gas (:)        ! User defined index for non radical gas species
  integer, allocatable :: ind_gas_rev (:)    ! revert old hard coded indexes. Should disappear soon...
  integer, allocatable :: ind_gas_br  (:,:)  ! Same for halogenated species:
  integer, allocatable :: ind_gas_cl  (:,:)  !   - first column holds stoechiometric nb
!  integer, allocatable :: ind_gas_iod (:,:)  !   - second column holds the (compressed) index

  ! Logical
  logical, allocatable :: gas_is_halo (:)    ! True is halogenated (Br, Cl, I) gas. Used to initialise differently these gases.

  ! Names
  character(len=12),  allocatable :: gas_name (:)
  character(len=100), allocatable :: gas_name_long (:)

  ! Data:
  !   - molar mass of non radical gas [kg/mol]
  real (kind=dp), allocatable :: gas_mass (:)
  !   - concentration of non radical gas [mol/m3]
  real (kind=dp), allocatable :: s1 (:,:)
  !   - user defined initial concentrations at ground and top level [ppb]
  real (kind=dp), allocatable :: s1_init_grd(:)
  real (kind=dp), allocatable :: s1_init_top(:)
  !   - ground emission of the gas []
  real (kind=dp), allocatable :: es1 (:)
  real (kind=dp), allocatable :: vg (:)

! Conversion tables, gas species
  ! KPP to Mistra: in Mistra order
  integer, allocatable :: gas_k2m_g (:)
  integer, allocatable :: gas_k2m_a (:)
  integer, allocatable :: gas_k2m_t (:)
  ! Mistra to KPP: in KPP order
  integer, allocatable :: gas_m2k_g (:,:)
  integer, allocatable :: gas_m2k_a (:,:)
  integer, allocatable :: gas_m2k_t (:,:)

  ! Advection (Eulerian view)
  integer :: nadvmax
  integer, allocatable :: nindadv (:)
  real (kind=dp), allocatable :: xadv (:)


! RADICALS (gas phase)
  ! Numbers of radical species
  integer :: j5                            ! Number of radicals actually used
  integer :: j5_br, j5_cl, j5_iod, j5_halo ! Number of brominated / chlorinated / ...  radicals

  ! Indexes
  integer, allocatable :: ind_rad (:)        ! User defined index for radical species
  integer, allocatable :: ind_rad_br  (:,:)  ! Same for halogenated species:
  integer, allocatable :: ind_rad_cl  (:,:)  !   - first column holds stoechiometric nb
!  integer, allocatable :: ind_rad_iod (:,:)  !   - second column holds the (compressed) index

  ! Logical
!  logical, allocatable :: rad_is_halo (:)

  ! Names
  character(len=12),  allocatable :: rad_name (:)
  character(len=100), allocatable :: rad_name_long (:)

  ! Data:
  !   - molar mass of radical [kg/mol]
  real (kind=dp), allocatable :: rad_mass (:)
  !   - concentration of radical [mol/m3]
  real (kind=dp), allocatable :: s3 (:,:)

! Conversion tables, raedical species
  ! KPP to Mistra: in Mistra order
  integer, allocatable :: rad_k2m_g (:)
  integer, allocatable :: rad_k2m_a (:)
  integer, allocatable :: rad_k2m_t (:)
  ! Mistra to KPP: in KPP order
  integer, allocatable :: rad_m2k_g (:,:)
  integer, allocatable :: rad_m2k_a (:,:)
  integer, allocatable :: rad_m2k_t (:,:)


! FIXED species which are exchanged
  integer, parameter :: j4 = 1
  real (kind=dp) :: fix_mass(j4)
  data fix_mass /32.d-3/

  real (kind=dp), allocatable :: vmean_init(:)
  real (kind=dp), allocatable :: vmean(:,:)

end module gas_common
