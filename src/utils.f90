!
! Copyright 2015-2017 Josue Bock
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


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine mk_interface

! Description :
! -----------
!   Prepare the interface between Mistra and KPP

! Author
! ------
!     Josue Bock


USE config, ONLY : &
     cgaslistfile, &
     cradlistfile, &
     coutdir

USE file_unit, ONLY : &
     jpfunint

USE gas_common, ONLY :    &
     j1,                  &
     gas_name,            &
     gas_m2k_g,gas_k2m_g, &
     gas_m2k_a,gas_k2m_a, &
     gas_m2k_t,gas_k2m_t, &
     ind_gas,ind_gas_rev, &

     j5,                  &
     rad_name,            &
     rad_m2k_g,rad_k2m_g, &
     rad_m2k_a,rad_k2m_a, &
     rad_m2k_t,rad_k2m_t

!USE global_params, ONLY : &
!     j2,                  &
!     j3

USE kpp_gas_Global, ONLY : &
     nspec_g=>NSPEC,       &
     spc_names_g=>SPC_NAMES

USE kpp_aer_Global, ONLY : &
     nspec_a=>NSPEC,       &
     spc_names_a=>SPC_NAMES

USE kpp_tot_Global, ONLY : &
     nspec_t=>NSPEC,       &
     spc_names_t=>SPC_NAMES

implicit none


! =======================================================
! -- 1. -- Output file to write initialisation information
! =======================================================
open(jpfunint,file=trim(coutdir)//'interface.out')

call global_parameters_check

! =======================================================
! -- 2.1 -- Non radical gas species
! =======================================================
call read_mistra_gas_data (trim(cgaslistfile))

! Allocate KPP to Mistra tables (in Mistra order) and match indexes
allocate ( gas_k2m_g(j1), gas_k2m_a(j1), gas_k2m_t(j1) )

call match_mk_indexes (j1,gas_name,nspec_g,spc_names_g,gas_k2m_g)
call match_mk_indexes (j1,gas_name,nspec_a,spc_names_a,gas_k2m_a)
call match_mk_indexes (j1,gas_name,nspec_t,spc_names_t,gas_k2m_t)

! Allocate Mistra to KPP tables (in KPP order) and fill-out these tables by sorting k2m tables
allocate (gas_m2k_g(2,j1),gas_m2k_a(2,j1),gas_m2k_t(2,j1))

call sort_2D_array (j1,gas_k2m_g,gas_m2k_g)
call sort_2D_array (j1,gas_k2m_a,gas_m2k_a)
call sort_2D_array (j1,gas_k2m_t,gas_m2k_t)  

! Build temporary revert indexes array for hard coded indexes
allocate ( ind_gas_rev(ind_gas(j1)) )
call revert_index (j1,ind_gas(j1),ind_gas,ind_gas_rev)

write(jpfunint,*)"Exchange lists for non radical gas species"
write(jpfunint,*)"For gas mechanism"
write(jpfunint,*)"list in Mistra order"
write(jpfunint,*)gas_k2m_g(:)
write(jpfunint,*)"list in KPP order"
write(jpfunint,*)gas_m2k_g(:,:)
write(jpfunint,*)"For aer mechanism"
write(jpfunint,*)"list in Mistra order"
write(jpfunint,*)gas_k2m_a(:)
write(jpfunint,*)"list in KPP order"
write(jpfunint,*)gas_m2k_a(:,:)
write(jpfunint,*)"For tot mechanism"
write(jpfunint,*)"list in Mistra order"
write(jpfunint,*)gas_k2m_t(:)
write(jpfunint,*)"list in KPP order"
write(jpfunint,*)gas_m2k_t(:,:)

! Define j2
!j2 = j1 + j3


! =======================================================
! -- 2.2 -- Radical gas species
! =======================================================
call read_mistra_rad_data (trim(cradlistfile))

! Allocate KPP to Mistra tables (in Mistra order) and match indexes
allocate ( rad_k2m_g(j5), rad_k2m_a(j5), rad_k2m_t(j5) )

call match_mk_indexes (j5,rad_name,nspec_g,spc_names_g,rad_k2m_g)
call match_mk_indexes (j5,rad_name,nspec_a,spc_names_a,rad_k2m_a)
call match_mk_indexes (j5,rad_name,nspec_t,spc_names_t,rad_k2m_t)

! Allocate Mistra to KPP tables (in KPP order) and fill-out these tables by sorting k2m tables
allocate (rad_m2k_g(2,j5),rad_m2k_a(2,j5),rad_m2k_t(2,j5))

call sort_2D_array (j5,rad_k2m_g,rad_m2k_g)
call sort_2D_array (j5,rad_k2m_a,rad_m2k_a)
call sort_2D_array (j5,rad_k2m_t,rad_m2k_t)

write(jpfunint,*)"Exchange lists for radical gas species"
write(jpfunint,*)"For gas mechanism"
write(jpfunint,*)"list in Mistra order"
write(jpfunint,*)rad_k2m_g(:)
write(jpfunint,*)"list in KPP order"
write(jpfunint,*)rad_m2k_g(:,:)
write(jpfunint,*)"For aer mechanism"
write(jpfunint,*)"list in Mistra order"
write(jpfunint,*)rad_k2m_a(:)
write(jpfunint,*)"list in KPP order"
write(jpfunint,*)rad_m2k_a(:,:)
write(jpfunint,*)"For tot mechanism"
write(jpfunint,*)"list in Mistra order"
write(jpfunint,*)rad_k2m_t(:)
write(jpfunint,*)"list in KPP order"
write(jpfunint,*)rad_m2k_t(:,:)


close(jpfunint)

end subroutine mk_interface
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine global_parameters_check

! Description :
! -----------
!   Check the consistency of a few global parameters, as compared to KPP parameters

! Author
! ------
!     Josue Bock


USE config, ONLY : &
! Imported Routines:
     abortM

USE file_unit, ONLY : &
     jpfunerr, jpfunint

USE kpp_tot_Global, ONLY : &
     nreact_t=>NREACT

USE global_params, ONLY : &
     nrxn

implicit none

if(nrxn.lt.nreact_t) then
   write(jpfunerr,*)"Error in global_params.f90"
   write(jpfunerr,*)"  nrxn < NREACT (tot mechanism)"
   write(jpfunerr,*)"  nrxn should be set to: ",nreact_t
   call abortM('Error in SR global_parameters_check')
else if(nrxn.gt.nreact_t) then
   write(jpfunint,*)"Warning: in global_params.f90"
   write(jpfunint,*)"  nrxn > NREACT (tot mechanism)"
   write(jpfunint,*)"  nrxn can be set to: ",nreact_t
end if
end subroutine global_parameters_check
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine read_mistra_gas_data (file_name)

!     Imports and checks user defined gas data file,
!     and creates conversion tables between Mistra and KPP indexes
!

! Author
! ------
!     Josue Bock

USE config, ONLY : &
! Imported Parameters:
     cmechdir,     &
     halo,         &
     iod,          &
! Imported Routines:
     abortM

USE file_unit, ONLY : &
! Imported Parameters:
     jpfunerr, jpfungas, jpfunint

USE global_params, ONLY : &
     n

USE gas_common, ONLY : &
! Imported Scalar Variables with intent (out):
     j1,               & ! number of gases actually used
     j1_br,            & ! number of brominated species
     j1_cl,            & ! number of chlorinated species
     j1_halo,          & ! number of halogenated (Br and/or Cl) species
     j1_iod,           & ! number of iodinated species

! Imported Array Variables with intent (out):
     ind_gas,          &
     ind_gas_br,       & ! includes the number of atoms (first dimension), for mass check
     ind_gas_cl,       & ! includes the number of atoms (first dimension), for mass check
!     ind_gas_halo,     & ! includes only the index of halogenated species (Br, Cl) for switch purpose
!     ind_gas_iod,      & ! includes only the index of iodinated species for switch purpose
     gas_is_halo,      &
     gas_name,         &
     gas_name_long,    &
     gas_mass,         &
     s1,               &
     s1_init_grd,      & ! initial gas concentrations, ground level (in ppb)
     s1_init_top,      & ! initial gas concentrations, top level (in ppb)
     es1,              &
     vg

USE kpp_gas_Global, ONLY : &
     nspec_g=>NSPEC

USE precision, ONLY : &
! Imported Parameters:
     dp

implicit none

! Local parameters:
character (len=*), intent(in) :: file_name

! Local scalars:
integer :: nb_gas      ! The total number of gas species actually used
integer :: ind_tmp     ! User defined index
integer :: igj         ! ind_gas(j)
integer :: i, i_br, i_cl, i_iod, j, k
integer :: ind_br_start, ind_cl_start ! final write section
integer :: nb_br, nb_cl               ! final write section

character (len=300) :: line          ! Line read in files
character (len=12)  :: name_tmp      ! Gas name read in file
character (len=100) :: name_long_tmp ! Long name of rad read in file
character (len=3) halo_txt, iod_txt  ! final write section

real (kind=dp) :: mass_tmp ! Gas mass read in file

! NB: the user enters ground and top level mixing ratio in [nmol/mol]=[ppb]
!     there are converted to [mol/m3] later during the initialisation (see SR initc)
real (kind=dp) :: conc_ground_tmp ! Mixing ratio at ground level [ppb]
real (kind=dp) :: conc_top_tmp    ! Mixing ratio at top level [ppb]

real (kind=dp) :: emission_tmp    ! Emission rates of gas phase species [molecules/(cm2*s)]

! Local arrays:
character (len=12),  allocatable :: gas_name_tmp(:)      ! Gas name read in file
character (len=100), allocatable :: gas_name_long_tmp(:) ! Long gas name read in file

integer, allocatable :: ind_gas_tmp(:)
integer, allocatable :: ind_gas_br_tmp(:,:)
integer, allocatable :: ind_gas_cl_tmp(:,:)

logical, allocatable :: gas_is_halo_tmp(:)
logical, allocatable :: gas_is_iod_tmp(:)

real (kind=dp), allocatable :: gas_mass_tmp(:) ! Gas mass read in file
real (kind=dp), allocatable :: s1_init_grd_tmp(:)
real (kind=dp), allocatable :: s1_init_top_tmp(:)
real (kind=dp), allocatable :: es1_tmp(:)

! External function:
integer, external :: get_atom_nb

!- End of header ---------------------------------------------------------------


! ==============================================================================
! -- 1 --  Allocate temporary arrays using species number from KPP files
! ==============================================================================
allocate (ind_gas_tmp(nspec_g))
allocate (ind_gas_br_tmp(2,nspec_g))
allocate (ind_gas_cl_tmp(2,nspec_g))

allocate (gas_is_halo_tmp(nspec_g))
allocate (gas_is_iod_tmp(nspec_g))

allocate (gas_name_tmp(nspec_g))
allocate (gas_name_long_tmp(nspec_g))
allocate (gas_mass_tmp(nspec_g))
allocate (s1_init_grd_tmp(nspec_g))
allocate (s1_init_top_tmp(nspec_g))
allocate (es1_tmp(nspec_g))


! ==============================================================================
! -- 2 --  Read the file to get all data
! ==============================================================================

! Initialisation, local
nb_gas = 0
! Initialisation, global
j1 = 0       ! Total number of non radical gases actually used
j1_br = 0    ! Total number of brominated non radical gases
j1_cl = 0    ! Total number of chlorinated non radical gases
j1_iod = 0   ! Total number of iodinated non radical gases
j1_halo = 0  ! Total number of halogenated non radical gases


open(unit=jpfungas,file=trim(cmechdir)//file_name,status='old',err=100)


do
   ! Read a line as a string
   read(jpfungas,'(a300)',err=101,end=102)line

   ! Eliminates the leading and trailing blanks
   line = trim(adjustl(line))
   ! Eliminates the leading tabs (ASCII code 9) and embedded blanks
   do while (line(1:1).eq.achar(9))
      line = line(2:len_trim(line))
      line = trim(adjustl(line))
   end do

   ! Ignore blank lines, or lines starting with an exclamation mark
   if(len_trim(line).eq.0 .or. line(1:1).eq.'!') cycle

   ! Eliminates the comment at end of line, if present
   if(scan(line,'!').gt.0) line = line(1:scan(line,'!')-1)

   ! Converts the string into the relevant gas species data
   read(line,*,err=101)ind_tmp,name_tmp,mass_tmp, &
                       conc_ground_tmp,conc_top_tmp, &
                       emission_tmp, &
                       name_long_tmp

! ....................................................
! Skip if halo = .false. and name contains Cl, Br or I
!   or if iod  = .false. and name contains I
   i_br = index(name_tmp,'Br')
   i_cl = index(name_tmp,'Cl')
   i_iod = index(name_tmp,'I')

   if(.not.(halo) .and. i_br.gt.0) cycle
   if(.not.(halo) .and. i_cl.gt.0) cycle
   if(.not.(halo) .and. i_iod.gt.0) cycle
   if(.not.(iod)  .and. i_iod.gt.0) cycle
! ....................................................

! --------------------------------------------
! At this stage, the specie has to be recorded
   nb_gas = nb_gas+1
! --------------------------------------------

! Check that user data are correct:
   ! All indexes must be greater than zero
   if(ind_tmp.le.0) then
      write(jpfunerr,*)"Error in "//file_name
      write(jpfunerr,*)"  All indexes must be strictly positive"
      write(jpfunerr,*)"  See after index ",ind_gas_tmp(max(1,nb_gas-1))
      write(jpfunerr,*)"  Index read is: ",ind_tmp
      close(jpfungas)
      call abortM("Error 1 in SR read_mistra_gas_data")
   end if

   ! Molar masses must be greater than zero
   if(mass_tmp.le.0.) then
      write(jpfunerr,*)"Error in "//file_name
      write(jpfunerr,*)"  All molar masses must be strictly positive"
      write(jpfunerr,*)"  See after index ",ind_gas_tmp(max(1,nb_gas-1))
      write(jpfunerr,*)"  Index read is: ",ind_tmp
      close(jpfungas)
      call abortM("Error 2 in SR read_mistra_gas_data")
   end if

   ! Concentrations must be greater equal to zero
   if(conc_ground_tmp.lt.0. .or. conc_top_tmp.lt.0.) then
      write(jpfunerr,*)"Error in "//file_name
      write(jpfunerr,*)"  All concentrations must be positive"
      write(jpfunerr,*)"  See after index ",ind_gas_tmp(max(1,nb_gas-1))
      write(jpfunerr,*)"  Index read is: ",ind_tmp
      close(jpfungas)
      call abortM("Error 3 in SR read_mistra_gas_data")
   end if

   ! All emissions must be greater equal to zero
   if(emission_tmp.lt.0.) then
      write(jpfunerr,*)"Error in "//file_name
      write(jpfunerr,*)"  All emissions must be positive"
      write(jpfunerr,*)"  See after index ",ind_gas_tmp(max(1,nb_gas-1))
      write(jpfunerr,*)"  Index read is: ",ind_tmp
      close(jpfungas)
      call abortM("Error 4 in SR read_mistra_gas_data")
   end if


   if(nb_gas.gt.1) then
      ! Check that indexes are ordered
      if(ind_tmp.le.ind_gas_tmp(nb_gas-1)) then
         write(jpfunerr,*)"Error in "//file_name//", sort out gas indexes"
         write(jpfunerr,*)"  Gas nb: ",ind_tmp
         close(jpfungas)
         call abortM("Error 5 in SR read_mistra_gas_data")
      end if

      ! Check that names are all different
      do i=1,nb_gas-1
         if(name_tmp.eq.gas_name_tmp(i)) then
            write(jpfunerr,*)"Error in "//file_name//", two identical names"
            write(jpfunerr,*)"  Gas numbers: ",ind_gas_tmp(i),ind_tmp
            write(jpfunerr,*)"  Gas name: ",name_tmp
            close(jpfungas)
            call abortM("Error 6 in SR read_mistra_gas_data")
         end if
      end do
   end if
! All checks are done


! Replace "_" by " " in long names
   i = index(name_long_tmp,"_")
   do while(i.gt.0)
      name_long_tmp(i:i) = " "
      i = index(name_long_tmp,"_")
   end do


! Store values read in final arrays (with user defined index)
   ind_gas_tmp(nb_gas) = ind_tmp
   gas_name_tmp(nb_gas) = trim(name_tmp)
   gas_name_long_tmp(nb_gas) = trim(name_long_tmp)
   gas_mass_tmp(nb_gas) = mass_tmp
   s1_init_grd_tmp(nb_gas) = conc_ground_tmp ! ground concentrations are layer 2
   s1_init_top_tmp(nb_gas) = conc_top_tmp
   es1_tmp(nb_gas) = emission_tmp

! Deal with halogenated species
   if(i_br.gt.0 .or. i_cl.gt.0) then
      ! Feed halo list with gas index only
      j1_halo = j1_halo+1
!      ind_gas_halo_tmp(j1_halo) = nb_gas
      gas_is_halo_tmp(nb_gas) = .TRUE.

      ! Feed Br list, with atoms numbers (for mass check) and gas index
      if(i_br.gt.0) then
         j1_br = j1_br+1
         ind_gas_br_tmp(1,j1_br) = get_atom_nb(name_tmp,'Br')
         ind_gas_br_tmp(2,j1_br) = nb_gas
      end if
      ! Feed Cl list, with atoms numbers (for mass check) and gas index
      if(i_cl.gt.0) then
         j1_cl = j1_cl+1
         ind_gas_cl_tmp(1,j1_cl) = get_atom_nb(name_tmp,'Cl')
         ind_gas_cl_tmp(2,j1_cl) = nb_gas
      end if
   else
      gas_is_halo_tmp(nb_gas) = .FALSE.
   end if

! Deal with iodinated specied
   if(i_iod.gt.0) then
      j1_iod = j1_iod+1
!      ind_gas_iod_tmp(j1_iod) = ind_gas(nb_gas)
      gas_is_halo_tmp(nb_gas) = .TRUE.            ! note that I ==> halo ...
      gas_is_iod_tmp(nb_gas) = .TRUE.
   else
      gas_is_iod_tmp(nb_gas) = .FALSE.
      ! gas_is_halo_tmp(ind_gas(nb_gas)) = ?      ! ... but not(I) =/=> not(halo)
   end if

end do ! read all lines

! Missing file
 100  write(jpfunerr,*)"Error: the file '"//file_name//"' is missing in mech directory"
      call abortM("Error 100 in SR read_mistra_gas_data")

! Error during read
 101  write(jpfunerr,*)"Error when reading "//file_name
      write(jpfunerr,*)"  Check if numbers are correctly written,"
      write(jpfunerr,*)"  and if line breaks are unix-like."
      write(jpfunerr,*)line
      close(jpfungas)
      call abortM("Error 101 in SR read_mistra_gas_data")

! End of file
 102  close(jpfungas)

j1 = nb_gas

! ==============================================================================
! -- 3.1 --  Allocate global arrays
! ==============================================================================

! Indexes
allocate ( ind_gas(j1) )
! allocate ( ind_gas_rev(ind_gas_tmp(nb_gas)) ) ! utilitary array to revert old hard coded indexes
allocate ( ind_gas_br(2,j1_br) )
allocate ( ind_gas_cl(2,j1_cl) )

! Logical
allocate ( gas_is_halo(j1) )    ! Used to initialise differently these gases.
                                      !   see SR initc in kpp.f
! Names
allocate ( gas_name(j1) )
allocate ( gas_name_long(j1) )

! Data
allocate ( gas_mass(j1) )

! Other arrays related to non radical gases
allocate ( s1(j1,n) )
allocate ( s1_init_grd(j1) )
allocate ( s1_init_top(j1) )
allocate ( es1(j1) )
allocate ( vg(j1) )


! ==============================================================================
! -- 3.2 --  Reshape tmp arrays into global arrays
! ==============================================================================

ind_gas(:) = ind_gas_tmp(:j1)
ind_gas_br(:,:) = ind_gas_br_tmp(:,:j1_br)
ind_gas_cl(:,:) = ind_gas_cl_tmp(:,:j1_cl)
gas_name(:) = gas_name_tmp(:j1)
gas_name_long(:) = gas_name_long_tmp(:j1)
gas_is_halo(:) = gas_is_halo_tmp(:j1)
gas_mass(:) = gas_mass_tmp(:j1)
s1_init_grd(:) = s1_init_grd_tmp(:j1)
s1_init_top(:) = s1_init_top_tmp(:j1)
es1(:) = es1_tmp(:j1)


! ==============================================================================
! -- 3.3 --  Write non radical gas species list, to check that retrieval went well
! ==============================================================================

write(jpfunint,*)"Information: the number of gases read in ",file_name," is ",j1
write(jpfunint,*)"j1_br= ",j1_br," j1_cl= ",j1_cl," j1_iod= ",j1_iod
write(jpfunint,*)"j1_halo= ",j1_halo
write(jpfunint,*)"Index  Name      halo? iod ? n_Br n_Cl mass"
do j=1,j1
   igj = ind_gas(j)
   ind_br_start = 1
   ind_cl_start = 1
   nb_br = 0
   nb_cl = 0
   if(gas_is_halo(j)) then
      halo_txt = 'yes'
      do k=ind_br_start,j1_br
         if(ind_gas_br(2,k).eq.j) then
            nb_br = ind_gas_br(1,k)
            ind_br_start = ind_br_start + 1
            exit
         else if(ind_gas_br(2,k).gt.j) then
            nb_br = 0
            exit
         end if
      end do
      do k=ind_cl_start,j1_cl
         if(ind_gas_cl(2,k).eq.j) then
            nb_cl = ind_gas_cl(1,k)
            ind_cl_start = ind_cl_start + 1
            exit
         else if(ind_gas_cl(2,k).gt.j) then
            nb_cl = 0
            exit
         end if
      end do
   else
      halo_txt = 'no'
   end if
   if(gas_is_iod_tmp(j)) then
      iod_txt = 'yes'
   else
      iod_txt = 'no'
   end if
         
   write(jpfunint,110)j,igj,gas_name(j),halo_txt,iod_txt,nb_br,nb_cl,gas_mass(j)
110 format(2i4,1x,a12,1x,2(a4,1x),2(i4,1x),f6.3)
end do

! ==============================================================================
! -- 3.4 --  Deallocate temporary arrays
! ==============================================================================
deallocate (ind_gas_tmp)
deallocate (ind_gas_br_tmp)
deallocate (ind_gas_cl_tmp)

deallocate (gas_is_halo_tmp)
deallocate (gas_is_iod_tmp)

deallocate (gas_name_tmp)
deallocate (gas_name_long_tmp)
deallocate (gas_mass_tmp)
deallocate (s1_init_grd_tmp)
deallocate (s1_init_top_tmp)
deallocate (es1_tmp)

end subroutine read_mistra_gas_data
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine read_mistra_rad_data (file_name)

!     Imports and checks user defined radical (gas phase) data file,
!     and creates conversion tables between Mistra and KPP indexes
!

! Author
! ------
!     Josue Bock

USE config, ONLY : &
! Imported Parameters:
     cmechdir,     &
     halo,         &
     iod,          &
! Imported Routines:
     abortM

USE file_unit, ONLY : &
! Imported Parameters:
     jpfunerr, jpfunint, jpfunrad

USE global_params, ONLY : &
     n

USE gas_common, ONLY : &
! Imported Scalar Variables with intent (out):
     j5,               &
     j5_br,            & ! number of brominated species
     j5_cl,            & ! number of chlorinated species
     j5_iod,           & ! number of iodinated species
     j5_halo,          & ! number of halogenated species
     
! Imported Array Variables with intent (out):
     ind_rad,          &
     ind_rad_br,       & ! includes the number of atoms (first dimension), for mass check
     ind_rad_cl,       & ! includes the number of atoms (first dimension), for mass check
!     ind_rad_halo,     & ! includes only the index of halogenated species (Br, Cl) for switch purpose
!     ind_rad_iod,      & ! includes only the index of iodinated species for switch purpose
!     rad_is_halo,      &
     rad_name,         &
     rad_name_long,    &
     rad_mass,         &
     s3

USE kpp_gas_Global, ONLY : &
     nspec_g=>NSPEC

USE precision, ONLY : &
! Imported Parameters:
     dp

implicit none

! Local parameters:
character (len=*), intent(in) :: file_name

! Local scalars:
integer :: nb_rad      ! The total number of rad species actually used
integer :: ind_tmp     ! User defined index
integer :: irj         ! ind_rad(j)
integer :: i, i_br, i_cl, i_iod, j, k
integer :: ind_br_start, ind_cl_start ! final write section
integer :: nb_br, nb_cl               ! final write section

character (len=300) :: line          ! Line read in files
character (len=12)  :: name_tmp      ! Rad name read in file
character (len=100) :: name_long_tmp ! Long name of rad read in file
character (len=3) halo_txt, iod_txt  ! final write section

real (kind=dp) :: mass_tmp ! Rad mass read in file

! Local arrays:
character (len=12),  allocatable :: rad_name_tmp(:)
character (len=100), allocatable :: rad_name_long_tmp(:)

integer, allocatable :: ind_rad_tmp(:)
integer, allocatable :: ind_rad_br_tmp(:,:)
integer, allocatable :: ind_rad_cl_tmp(:,:)
! integer, allocatable :: ind_rad_halo_tmp(:)

real (kind=dp), allocatable :: rad_mass_tmp(:)

logical, allocatable :: rad_is_halo_tmp(:)
logical, allocatable :: rad_is_iod_tmp(:)


! External function:
integer, external :: get_atom_nb

!- End of header ---------------------------------------------------------------


! ==============================================================================
! -- 1 --  Allocate temporary arrays using species number from KPP files
! ==============================================================================
allocate (ind_rad_tmp(nspec_g))
! allocate (ind_rad_halo_tmp(nspec_g))
allocate (ind_rad_br_tmp(2,nspec_g))
allocate (ind_rad_cl_tmp(2,nspec_g))

allocate (rad_is_halo_tmp(nspec_g))
allocate (rad_is_iod_tmp(nspec_g))

allocate (rad_name_tmp(nspec_g))
allocate (rad_name_long_tmp(nspec_g))

allocate (rad_mass_tmp(nspec_g))


! ==============================================================================
! -- 2 --  Read the file to get all data
! ==============================================================================

! Initialisation, local
nb_rad = 0
! Initialisation, global
j5_br = 0    ! Total number of brominated radicals
j5_cl = 0    ! Total number of chlorinated radicals
j5_iod = 0   ! Total number of iodinated radicals
j5_halo = 0  ! Total number of halogenated radicals


open(unit=jpfunrad,file=trim(cmechdir)//file_name,status='old',err=100)

do
   ! Read a line as a string
   read(jpfunrad,'(a300)',err=101,end=102)line

   ! Eliminates the leading and trailing blanks
   line = trim(adjustl(line))
   ! Eliminates the leading tabs (ASCII code 9) and embedded blanks
   do while (line(1:1).eq.achar(9))
      line = line(2:len_trim(line))
      line = trim(adjustl(line))
   end do

   ! Ignore blank lines, or lines starting with an exclamation mark
   if(len_trim(line).eq.0 .or. line(1:1).eq.'!') cycle

   ! Eliminates the comment at end of line, if present
   if(scan(line,'!').gt.0) line = line(1:scan(line,'!')-1)

   ! Converts the string into the relevant radical gas species data
   read(line,*,err=101)ind_tmp,name_tmp,mass_tmp,name_long_tmp

! ....................................................
! Skip if halo = .false. and name contains Cl, Br or I
!   or if iod  = .false. and name contains I
   i_br = index(name_tmp,'Br')
   i_cl = index(name_tmp,'Cl')
   i_iod = index(name_tmp,'I')

   if(.not.(halo) .and. i_br.gt.0) cycle
   if(.not.(halo) .and. i_cl.gt.0) cycle
   if(.not.(halo) .and. i_iod.gt.0) cycle
   if(.not.(iod)  .and. i_iod.gt.0) cycle
! ....................................................

! --------------------------------------------
! At this stage, the specie has to be recorded
   nb_rad = nb_rad+1
! --------------------------------------------

! Check that user data are correct:
   ! All indexes must be greater than zero
   if(ind_tmp.le.0) then
      write(jpfunerr,*)"Error in "//file_name
      write(jpfunerr,*)"  All indexes must be strictly positive"
      write(jpfunerr,*)"  See after index ",ind_rad_tmp(max(1,nb_rad-1))
      write(jpfunerr,*)"  Index read is: ",ind_tmp
      close(jpfunrad)
      call abortM("Error 1 in SR read_mistra_rad_data")
   end if

   if(nb_rad.gt.1) then
      ! Check that indexes are ordered
      if(ind_tmp.le.ind_rad_tmp(nb_rad-1)) then
         write(jpfunerr,*)"Error in "//file_name//", indexes have to be sorted out"
         write(jpfunerr,*)"  Rad nb: ",ind_tmp
         close(jpfunrad)
         call abortM("Error 2 in SR read_mistra_rad_data")
      end if

      ! Check that names are all different
      do i=1,nb_rad-1
         if(name_tmp.eq.rad_name_tmp(i)) then
            write(jpfunerr,*)"Error in "//file_name//", two identical names"
            write(jpfunerr,*)"  Radical numbers: ",ind_rad_tmp(i),ind_tmp
            write(jpfunerr,*)"  Radical name: ",name_tmp
            close(jpfunrad)
            call abortM("Error 3 in SR read_mistra_rad_data")
         end if
      end do
   end if
! All checks are done


! Replace "_" by " " in long names
   i = index(name_long_tmp,"_")
   do while(i.gt.0)
      name_long_tmp(i:i) = " "
      i = index(name_long_tmp,"_")
   end do

! Store values read in final arrays (with user defined index)
   ind_rad_tmp(nb_rad) = ind_tmp
   rad_name_tmp(nb_rad) = trim(name_tmp)
   rad_name_long_tmp(nb_rad) = trim(name_long_tmp)
   rad_mass_tmp(nb_rad) = mass_tmp

! Deal with halogenated species
   if(i_br.gt.0 .or. i_cl.gt.0) then
      ! Feed halo list with gas index only
      j5_halo = j5_halo+1
!      ind_rad_halo_tmp(j5_halo) = nb_rad
      rad_is_halo_tmp(nb_rad) = .TRUE.

      ! Feed Br list, with atoms numbers (for mass check) and gas index
      if(i_br.gt.0) then
         j5_br = j5_br+1
         ind_rad_br_tmp(1,j5_br) = get_atom_nb(name_tmp,'Br')
         ind_rad_br_tmp(2,j5_br) = nb_rad
      end if
      ! Feed Cl list, with atoms numbers (for mass check) and gas index
      if(i_cl.gt.0) then
         j5_cl = j5_cl+1
         ind_rad_cl_tmp(1,j5_cl) = get_atom_nb(name_tmp,'Cl')
         ind_rad_cl_tmp(2,j5_cl) = nb_rad
      end if
   else
      rad_is_halo_tmp(nb_rad) = .FALSE.
   end if

! Deal with iodinated specied
   if(i_iod.gt.0) then
      j5_iod = j5_iod+1
      rad_is_iod_tmp(nb_rad) = .TRUE.
      rad_is_halo_tmp(nb_rad) = .TRUE.   ! note that I ==> halo ...
   else
      rad_is_iod_tmp(nb_rad) = .FALSE.
      ! rad_is_halo_tmp(nb_rad) = ?      ! ... but not(I) =/=> not(halo)
   end if

end do ! read all lines

! Missing file
 100  write(jpfunerr,*)"Error: the file '"//file_name//"' is missing in mech directory"
      call abortM("Error 100 in SR read_mistra_rad_data")

! Error during read
 101  write(jpfunerr,*)"Error when reading "//file_name
      write(jpfunerr,*)"  Check if numbers are correctly written,"
      write(jpfunerr,*)"  and if line breaks are unix-like."
      write(jpfunerr,*)line
      close(jpfunrad)
      call abortM("Error 101 in SR read_mistra_rad_data")

! End of file
 102  close(jpfunrad)

j5 = nb_rad


! ==============================================================================
! -- 3.1 --  Allocate global arrays
! ==============================================================================

! Indexes
allocate ( ind_rad(j5) )
allocate ( ind_rad_br(2,j5_br) )
allocate ( ind_rad_cl(2,j5_cl) )

! Logical
! allocate ( rad_is_halo(j5) )       ! currently not used

! Names
allocate ( rad_name(j5) )
allocate ( rad_name_long(j5) )

! Data
allocate ( rad_mass(j5) )

! Other arrays related to radicals
allocate ( s3(j5,n) )


! ==============================================================================
! -- 3.2 --  Reshape tmp arrays into global arrays
! ==============================================================================

ind_rad(:) = ind_rad_tmp(:j5)
ind_rad_br(:,:) = ind_rad_br_tmp(:,:j5_br)
ind_rad_cl(:,:) = ind_rad_cl_tmp(:,:j5_cl)
rad_name(:) = rad_name_tmp(:j5)
rad_name_long(:) = rad_name_long_tmp(:j5)
rad_mass(:) = rad_mass_tmp(:j5)

! ==============================================================================
! -- 3.3 --  Write radical gas species list, to check that retrieval went well
! ==============================================================================
write(jpfunint,*)"Information: the number of radicals read in ",trim(file_name)," is: j5 = ",j5
write(jpfunint,*)"j5_br= ",j5_br," j5_cl= ",j5_cl," j5_iod= ",j5_iod
write(jpfunint,*)"j5_halo = ",j5_halo
write(jpfunint,*)"Index  Name      halo ? iod ? n_Br n_Cl mass"
do j=1,j5
   irj = ind_rad(j)
   ind_br_start = 1
   ind_cl_start = 1
   nb_br = 0
   nb_cl = 0
   if(rad_is_halo_tmp(j)) then
      halo_txt = 'yes'
      do k=ind_br_start,j5_br
         if(ind_rad_br(2,k).eq.j) then
            nb_br = ind_rad_br(1,k)
            ind_br_start = ind_br_start + 1
            exit
         else if(ind_rad_br(2,k).gt.j) then
            nb_br = 0
            exit
         end if
      end do
      do k=ind_cl_start,j5_cl
         if(ind_rad_cl(2,k).eq.j) then
            nb_cl = ind_rad_cl(1,k)
            ind_cl_start = ind_cl_start + 1
            exit
         else if(ind_rad_cl(2,k).gt.j) then
            nb_cl = 0
            exit
         end if
      end do
   else
      halo_txt = 'no'
   end if
   if(rad_is_iod_tmp(j)) then
      iod_txt = 'yes'
   else
      iod_txt = 'no'
   end if

   write(jpfunint,110)j,irj,rad_name(j),halo_txt,iod_txt,nb_br,nb_cl,rad_mass(j)
110 format(2i4,1x,a12,1x,2(a4,1x),2(i4,1x),f6.3)
end do

! ==============================================================================
! -- 3.4 --  Deallocate temporary arrays
! ==============================================================================
deallocate (ind_rad_tmp)
deallocate (rad_name_tmp)
deallocate (rad_name_long_tmp)
deallocate (rad_mass_tmp)
deallocate (rad_is_halo_tmp)
deallocate (rad_is_iod_tmp)
! deallocate (ind_rad_halo_tmp)
deallocate (ind_rad_br_tmp)
deallocate (ind_rad_cl_tmp)

end subroutine read_mistra_rad_data
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine match_mk_indexes &
                (nm,m_name, &
                 nk,k_name, &
                 k2m_table)

! Important note:
!     New Mistra arrays are now allocated the size of actually used species.
!     Thus, m_name(1:nm) is full, and used indexes are 1, 2, 3, ..., nm
!     
!     k2m conversion table (in Mistra order) is thus a 1-D dimension array,
!       which contains the equivalent KPP indexes

! Author
! ------
!     Josue Bock


USE config, ONLY : &
! Imported Routines:
     abortM

USE file_unit, ONLY : &
! Imported Parameters:
     jpfunerr

implicit none

! Scalar arguments with intent(in):
integer, intent(in) :: nm ! Size for Mistra arrays
integer, intent(in) :: nk ! Size for KPP names array

! Array arguments with intent(in):
character(len=*), intent(in) :: m_name (nm)
character(len=*), intent(in) :: k_name (nk)

! Array arguments with intent(out):
integer, intent(out) :: k2m_table (nm)

! Local scalars:
integer :: i,j

!- End of header ---------------------------------------------------------------

do i=1,nm
   do j=1,nk
      if(m_name(i).eq.k_name(j)) then
         k2m_table(i) = j
         exit
      end if
      if(j.eq.nk .and. m_name(i).ne.k_name(j)) then
         write(jpfunerr,*)"Mistra species not found in KPP names list:"
         write(jpfunerr,*)"Species name: ",m_name(i)," index ",i
         call abortM(' Error in SR match_mk_indexes_2')
      end if
   end do
end do
end subroutine match_mk_indexes
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine sort_2D_array (n,array_in,array_out)
!
! This subroutine sorts array_in into the first column of array_out, and writes
!    in the second column or array_out the original index
!
! The overall aim is to make efficient conversion tables between Mistra and KPP
!    so that read / write are done in the proper order, to reduce computing time

! Author
! ------
!     Josue Bock


implicit none

integer, intent(in)  :: n               ! array size along the second dimension
integer, intent(in)  :: array_in (n)
integer, intent(out) :: array_out (2,n)

integer :: i,j
integer :: col2sort (n)

col2sort = array_in

! Neither max indexes (of column 1 or column 2) are necessaryly related to n
! Thus, it is easier to fill-out arrays from n to 1, and overwrite the maxloc value with 0
! rather than filling it from 1 to n and overwritting the minloc value with an unknown value
! (even if that would be easy to do)
do i=n,1,-1
   j = maxloc(col2sort,1)
   array_out(1,i) = col2sort(j)
   array_out(2,i) = j
   col2sort(j) = 0
end do

end subroutine sort_2D_array
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine revert_index (nc,nnc,array_comp,array_non_comp)
!
! This subroutine is used temporarily to keep hard coded indexes
!    and convert into compressed index
!
! Compressed indexes:
!  1 | 2 | 4 | 9 | 10 | 11 | 18 | ...
! Output:
!  1 | 2 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 4 | 5 | 0 | ...

! Author
! ------
!     Josue Bock


implicit none

integer, intent(in) :: nc
integer, intent(in) :: nnc
integer, intent(in) :: array_comp (nc)
integer, intent(out) :: array_non_comp (nnc)

integer :: ic,inc

ic = 1
do inc=1,nnc
   if(array_comp(ic) == inc) then
      array_non_comp(inc) = ic
      ic = ic + 1
      cycle
   else
      array_non_comp(inc) = 0
   end if
end do

end subroutine revert_index
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function get_atom_nb (spec_name,atom_name)
!
! Description:
!     This utilitary function returns the number of atoms in a molecule, after analysing the name.
!     Example: find the number of carbon atoms in C3H7Cl (which can be passed either as C3H7Cl, CH3CHClCH3)
!
!     The function first checks that the atom name is present.
!     Then it progressively reduces the species name, counts the numbers following the atoms,
!       and skip the atom name in the case of single letter atom name, followed by a lower case letter in the species name.
!     Last, this function produces a warning in case a number 0 or 1 is found following the atom name.
!     Thus, it is not conceived for very large molecules containing more than 9 atoms.
!
! Method:
!     note about ASCII code for characters: 48 - 57  <=> 0-9
!                                           97 - 122 <=> a-z
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.0      08/2016   Original code, developped to improve Mistra       <Josue Bock>
!                    (in order to remove hard coded species indexes)
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:

USE file_unit, ONLY : &
! Imported Parameters:
     jpfunint

implicit none

! Function declaration:
integer :: get_atom_nb

! Function arguments
! Scalar arguments with intent(in):
character (len=*), intent(in) :: spec_name
character (len=*), intent(in) :: atom_name

! Local scalars:
character (len=len_trim(spec_name)) :: s_name
character (len=len_trim(atom_name)) :: a_name

integer :: i      ! position of the atom name in the species name (if present)
integer :: ascii  ! ASCII code of the character following the atom name (if present)
integer :: lta    ! length of the atom name

!- End of header ---------------------------------------------------------------

! Initialize
get_atom_nb = 0
s_name = trim(spec_name)
a_name = trim(atom_name)
lta = len_trim(atom_name)

i = index(s_name,a_name)

!if(i.eq.0) then ! atom not found. Return with value 0
if(i.ne.0) then
   do while (len_trim(s_name)-lta.ge.0 .and. i.gt.0)
      if(len_trim(s_name).eq.lta) then ! The (remaining) species name is exactly the length of the atom name
         get_atom_nb = get_atom_nb + 1
         exit ! exit do while and return
      end if
      s_name = s_name(i+lta:len_trim(s_name))
      ascii = iachar(s_name(1:1))
      if(ascii.eq.48 .or. ascii.eq.49) then
         write(jpfunint,*)"Warning in function get_atom_nb:"
         write(jpfunint,*)"  numerical character 0 or 1 following atom"
         write(jpfunint,*)"  atom: ",atom_name," species: ",spec_name
         get_atom_nb = get_atom_nb + 1 ! suppose only one atom in this (unexpected) case
      else if(ascii.ge.50 .and. ascii.le.57) then
         get_atom_nb = get_atom_nb + ascii -48
       ! Next case:
       ! avoid single letter atom name ('C' for instance) followed by a lower case letter ('l') => not the searched atom
      else if(lta.gt.1 .or. lta.eq.1 .and. (ascii.lt.97 .or. ascii.gt.122)) then
         get_atom_nb = get_atom_nb + 1
      end if
      i = index(s_name,a_name)
   end do
end if
end function get_atom_nb
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
