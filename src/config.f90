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


module config

! Description :
! -----------
!   General configuration of the model
!   Declaration of all user-defined parameters
!   Contains utilitary subroutines to read/write these parameters

! Interface :
! ---------

! Input :
! -----

! Output :
! ------

! Externals :
! ---------

! Method :
! ------

! Author :
! ------
!   Josue Bock

! Modifications :
! -------------
  !   12-Mar-2017  Josue Bock  First version based on the initial parameter reading (file istart) which was done in MAIN program
  !
  !   03-Nov-2017  Josue Bock  Defined namelist /mistra_cfg/ that contains the former "istart" file parameters
!-----------------------------------------------------------------------------------------------------------------------------------



! Declarations:
! ------------
! Modules used:

use data_surface, only : &
     tw,&
     z0
  
use precision, only : &
  dp                    ! double precision kind


implicit none


public
save

logical :: &
  rst,     & ! rst      : restart or initialization of program run
  netCDF,  & ! netCDF   : output in netCDF format
  binout,  & ! binout   : output in binary format
  mic,     & ! mic      : microphysics included
  chem,    & ! chem     : chemistry included
  halo,    & ! halo     : halogen chemistry on/off
  iod,     & ! iod      : iodine chemistry on/off
  box,     & ! box      : box model run
  BL_box,  & ! BL_box   : box only, average init cond over BL and/or mean of J-values over BL
  chamber, & ! chamber  : chamber model run
  nuc,     & ! nuc      : nucleation on/off
  Napari,  & ! Napari   : nuc only, Napari = ternary H2SO4-H2O-NH3 nucleation
  Lovejoy, & ! Lovejoy  : nuc only, Lovejoy = homogeneous OIO nucleation
  ltwcst,  & ! ltwcst   : constant tw (if not, ntwopt must be set)
  lpmona,  & ! lpmona   : activate Monahan scheme for particle emission (see SR aer_source)
  lpsmith    ! lpsmith  : activate Smith scheme for particle emission (see SR aer_source)

integer :: &
  jpPartDistSet, &      ! Whole set of particle distributions (for one or several aerosol types)
  iaertyp, & ! iaertyp  : type of aerosol; 1=urban, 2=rural, 3=ocean, 4=background
  ifeed,   & ! ifeed    : retroaction over microphysics, and/or chemistry. See manual.
  isurf,   & ! isurf    : type of surface, (0) for water or snow, (1) for bare soil
  lstmax,  & ! lstmax   : integration time in hours
  neula,   & ! neula    : eulerian (0) or lagrangian (1) view
  nlevbox, & ! nlevbox  : box only, level to be used for init cond of box if  BL_box=false
  nkc_l,   & ! nkc_l    : number of output classes for aq. chem.
  ntwopt,  & ! ntwopt   : option for tw varying with time, see SR surf0
  jpOutPart2dOpt !      : option for netCDF output layers of 2D particle spectrum

integer :: &
     nday,        & ! starting time day
     nmonth,      & ! starting time month
     nyear,       & ! starting time year
     nhour          ! starting time hour
real (kind=dp) :: &
     alon,        & ! longitude (in degree, -180 ; 180)
     alat           ! latitude (in degree)

integer :: &
     nuvProfOpt,  & ! Option for the profile of geostrophic wind components, 0 = cst except layers 1-3, 3=linearly decreasing
     nwProfOpt      ! Option for the profile of subsidence, 1=BTZ96, 2=default, 3=linearly decreasing (Bott2000)
real (kind=dp) :: &
     rhMaxBL,     & ! Maximum relative humidity in the boundary layer (model initialisation)
     rhMaxFT,     & ! Maximum relative humidity above inversion = in the free troposphere
     rp0,         & ! Surface pressure [Pa]
     ug,          & ! geostrophic wind, x-direction [m/s]
     vg,          & ! geostrophic wind, y-direction [m/s]
     wmin,        & ! large scale subsidence (min) [m/s]
     wmax,        & ! large scale subsidence (max) [m/s]
     xm1w,        & ! specific humidity below inversion layer (kg/kg)
     xm1i,        & ! specific humidity above inversion layer (kg/kg)
     zinv,        & ! initial inversion height [m]
     dtinv          ! inversion strength = temperature drop at inversion [K]

real (KIND=dp) :: &
  detamin,        & ! atmospheric grid: height of constant layers [m]
  etaw1,          & ! atmospheric grid: top of the grid [m]
  rnw0,           & ! microphysics grid: min radius of dry aerosol [um]
  rnw1,           & ! microphysics grid: max radius of dry aerosol [um]
  rw0,            & ! microphysics grid: min radius of particle [um]
  rw1,            & ! microphysics grid: max radius of particle [um]
  rhsurf,         & ! rhsurf   : relative humidity at the surface, forced at each timestep (see SR surf0)
  scaleo3_m,      & ! scaleo3_m: total O3 in DU (for photolysis only)
  z_box             ! z_box    : height of MBL (if box run)

! Surface settings
integer :: jpAlbedoOpt ! albedo of the surface (set related configuration in radinit.f90)

! Special runs switchs
logical :: lpBuxmann15alph ! Switch on some special configuration to reproduce Buxmann et al 2015 (alpha case)
logical :: lpBuys13_0D ! Switch on some special configuration of Buys et al 2013 (0D case)
logical :: lpJoyce14bc ! Switch on some special configuration of Joyce et al 2014 (base case)

character (len=100) :: cnmlfile

character (len=100) :: cinpdir      ! input directory: general data files for Mistra
character (len=109) :: cinpdir_phot ! input directory for photolysis data files
character (len=100) :: coutdir      ! output directory
character (len=100) :: cmechdir     ! mechanism directory
character (len=100) :: cgaslistfile ! user file holding the list of gas phase species
character (len=100) :: cradlistfile ! user file holding the list of gas phase radical species


namelist /mistra_cfg/ &
     rst,             &
     lstmax,          &
     netcdf,          &
     binout,          &
     jpOutPart2dOpt,  &
! model grids
     detamin, etaw1, rnw0, rnw1, rw0, rw1, &
! timing and geography
     nday, nmonth, nyear, nhour, alon, alat,  &
! meteorological data
     rp0, zinv, dtinv, xm1w, xm1i, rhMaxBL, rhMaxFT, ug, vg, wmin, wmax, nuvProfOpt, nwProfOpt, &
! Surface setings
     isurf, tw, ltwcst, ntwopt, rhsurf, z0, jpAlbedoOpt, &
     mic,             &
     jpPartDistSet, iaertyp,      &
! Chemistry setings
     chem, halo, iod, nkc_l,      &
     cgaslistfile, cradlistfile,  &
     lpmona, lpsmith, &
     neula,           &
! Box settings
     box, bl_box, nlevbox, z_box, &
! Chamber settings
     chamber, &
! Nucleation settings
     nuc, ifeed, Napari, Lovejoy, &
     scaleo3_m,       &
! Special configuration
     lpBuxmann15alph, lpBuys13_0D, lpJoyce14bc

contains


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine read_config

! Description :
! -----------
!   Read configuration file


! Interface :
! ---------
!   SR read_config is called by main program during initialisation

! Input :
! -----
!   configuration file 'config.txt'

! Output :
! ------

! Externals :
! ---------
!   none

! Method :
! ------

! Author :
! ------
!   Josue Bock

! Modifications :
! -------------
  !   12-Mar-2017  Josue Bock  First version based on the initial parameter reading (file istart)
  !                            which was done in MAIN program
  !
  !      Nov-2017  Josue Bock  Rewritten (almost) from scratch, with a different structure (env variables, ...)
!--------------------------------------------------------------------------------------------------------------

  use file_unit, only : &
! Imported Parameters:
       jpfunnam,        &
       jpfunerr, jpfunout, &
       jpfuncfgout

  use precision, only : &
! Imported Parameters:
       dp

  implicit none

  character (len=3) :: clstat
  integer :: istat


! =======================================================
! -- 1. -- get I/O directories from environment variables
! =======================================================
  call getenv ('INPDIR',cinpdir)
  if (trim(cinpdir) == '') then
     call abortM ('Error during initialisation: no input directory specified as environment variable INPDIR')
  end if
  cinpdir_phot = trim(cinpdir)//'photolys/'

  call getenv ('OUTDIR',coutdir)
  if (trim(coutdir) == '') then
     call abortM ('Error during initialisation: no output directory specified as environment variable OUTDIR')
  end if

  call getenv ('MECHDIR',cmechdir)
  if (trim(cmechdir) == '') then
     call abortM ('Error during initialisation: no mechanism directory specified as environment variable MECHDIR')
  end if

! ===============================================================
! -- 2. -- Mistra_cfg namelist: default values, and read namelist
! ===============================================================
! Default values
rst = .false.
lstmax = 1
netCDF = .false.
binout = .false.
jpOutPart2dOpt = 0

! timing and geography
 nday = 01
 nmonth = 07
 nyear = 2021
 nhour = 0
 alon = 0._dp
 alat = 0._dp

! model grids
detamin = 10._dp
etaw1 = 2000._dp
rnw0 = 0.005_dp
rnw1 =  15._dp
rw0  = 0.005_dp
rw1  = 150._dp

! meteorological data
rp0 = 101325._dp
xm1w = 8.5e-3_dp
xm1i = 4.0e-3_dp
rhMaxBL = 1._dp
rhMaxFT = 1._dp
zinv = 700._dp
dtinv = 6._dp
ug = 6._dp
vg = 6._dp
nuvProfOpt = 0
nwProfOpt = 2
wmin = 0._dp
wmax = -0.006_dp

isurf = 0
tw = 293._dp
ltwcst = .true.
ntwopt = 1
rhsurf = 1._dp
z0 = 0.01_dp
jpAlbedoOpt = 0
mic = .false.
jpPartDistSet = 0
iaertyp = 3
! Chemistry setings
chem = .true.
halo = .true.
iod = .true.
nkc_l = 4
cgaslistfile='gas_species.csv'
cradlistfile='gas_radical_species.csv'
lpmona = .true.
lpsmith = .false.
neula = 1
! Box settings
box = .false.
bl_box = .false.
nlevbox = 2
z_box = 700._dp
! Chamber settings
chamber = .false.
! Nucleation settings
nuc = .false.
ifeed = 0
Napari = .true.
Lovejoy = .true.
! Ozone column (only for photolysis rates)
scaleo3_m = 300._dp

! Special configuration
lpBuxmann15alph = .false.
lpBuys13_0D = .false.
lpJoyce14bc = .false.

call getenv ('NAMELIST',cnmlfile)

if (trim(cnmlfile) /= '') then
   open (UNIT=jpfunnam, FILE=cnmlfile, STATUS='old', FORM='formatted', IOSTAT=istat)
   if (istat /= 0) call abortM ('Error in SR read_config: cannot open namelist file: '//cnmlfile)

   read (UNIT=jpfunnam, NML=mistra_cfg, IOSTAT=istat)
   if (istat /= 0) call abortM ('Error in SR read_config: cannot read namelist mistra_cfg in file: '//cnmlfile)

   close (UNIT=jpfunnam)
else
   write(jpfunout,'(a)') 'Warning: no namelist specified, only hardcoded default settings will be used'
end if

! ======================================================
! -- 3. -- Perform some checks over the resulting values
! ======================================================
if (nuc.and..not.chem) then
   nuc = .false.
   write(jpfunout,*) 'Warning: nuc has been set to false since chemistry is off'
end if
if (.not.nuc) then
   ifeed = 0
   ! do not display warning message in this case
end if
if (nuc.and.(.not.Napari).and.(.not.Lovejoy)) then
   write(jpfunerr,*) 'Error: Napari or Lovejoy must be true is nuc is used'
   call abortM ('Stopped by SR read_config')
end if

if (lpmona .and. lpsmith) then
   write (jpfunerr,*) 'Error in namelist settings: choose either lpmona or lpsmith for aer emission scheme'
   call abortM ('Stopped by SR read_config')
end if

if (chamber .and. box) then
   write (jpfunerr,*) 'Error in namelist settings: choose either chamber or box'
   call abortM ('Stopped by SR read_config')
end if
if (chamber .and. neula==0) then
   write (jpfunerr,*) 'Error in namelist settings: neula=0 is not possible in chamber mode'
   call abortM ('Stopped by SR read_config')
end if

if (BL_box.and..not.box) then
   BL_box = .false.
   write(jpfunout,*) 'Warning: BL_box has been set to false since box is off'
end if

if (chamber .and. jpPartDistSet.ne.4) then
   jpPartDistSet = 4
   write(jpfunout,*) 'Warning: jpPartDistSet has been set to 4, only possible choice for chamber version'
end if
if (chamber) then
   alat = 45
   write(jpfunout,*) 'Warning: alat has been set to 45 for chamber version'
end if

! =====================================================
! -- 4. -- Export current configuration in file cfg.out
! =====================================================

  if (rst) then
     clstat = 'old'
  else
     clstat = 'new'
  end if
  open (unit=jpfuncfgout, FILE=trim(coutdir)//'cfg.out', STATUS=clstat, FORM='formatted', POSITION='append', IOSTAT=istat)
  if (istat /= 0) call abortM ('Error in SR read_config: cannot open cfg.out file in dir: '//coutdir)

  write (jpfuncfgout,'(a)') 'Mistra configuration:'
  write (jpfuncfgout,'(a)') '  input directory:',cinpdir
  write (jpfuncfgout,'(a)') '  output directory:',coutdir
  write (jpfuncfgout,'(a)') '  mechanism directory:',cmechdir
  write (jpfuncfgout,'(a)') '  mistra_cfg namelist:'

  write (unit=jpfuncfgout, NML=mistra_cfg, IOSTAT=istat)
  if (istat /= 0) call abortM ('Error in SR read_config: cannot write current mistra_cfg namelist values in cfg.out file.')

  close (unit=jpfuncfgout)

end subroutine read_config
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine abortM (cderrmessage)

  ! Description :
  ! -----------
  !    abortM is called after and error arise somewhere in the code.
  !    The message is written in stderr file, and the program is stopped
  !
  !    Calling abort (Fortran intrinsic function) did not work well with the
  !    current calling param file (probably because of a returned non-zero
  !    STATUS value). This might be imporved in the future.

  ! Author :
  ! ------
  !    Josue Bock


  use file_unit, only : &
       jpfunerr
  implicit none
  character (len=*), intent(in) :: cderrmessage
  write (jpfunerr,'(a)') cderrmessage

  if (netcdf) then
     write (jpfunerr,'(a)') 'Trying to close netCDF files'
     call close_netcdf(mic,chem,box,nuc)
  end if

  stop '  --> stopped by SR abort'

end subroutine abortM
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module config
