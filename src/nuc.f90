!
! Copyright 1996-2017 the Authors
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
!                  *************************************
!                  * nucleation subroutines for MISTRA *
!                  *************************************
!
! Authors :
! -------
!     Susanne Pechtl (b. Marquart)
!     Roland von Glasow
!
!
! See detailed explanations in the header of SR appnucl
! See MISTRA manual for full references and details.
!
! This file contains the following subroutines and functions:
!     - appnucl2
!     - appnucl
!     - dmean
!     - ternucl
!     - oionucl
!     - J_nuc   (function)
!     - nucout1
!     - nucout2
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module mod_nuc

! Purpose:
! --------
!     - declare and set the number of nucleating vapours nvap
!     - declare variables used throughout the nucleation scheme


! Author :
! ------
!    Josue Bock


! Modifications :
! -------------
  ! 30-11-2016   Josue Bock   first version of this module
  ! 04-11-2017   Josue Bock   Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  use precision, only : &
! Imported Parameters:
       dp

  implicit none
  save

! Local parameterss:
  integer, parameter :: nvap = 1 ! Number of condensible vapors, please adjust!

! Local scalars:
  integer :: ind_H2SO4, ind_NH3, ind_OIO

! Local arrays:
  integer :: ical(nvap)
  integer :: ivap(nvap)
  real (kind=dp) :: m_vap(nvap), concsat(nvap)
  character (len=12) :: nuc_name(nvap)

! == End of declarations =======================================================

end module mod_nuc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine nuc_init (Napari,Lovejoy,iod)

! Purpose:
! --------
!     This SR checks consistency between several options, and gets indexes of
!     relevant species. The user has to declare the names of the species used in
!     the nucleation code.


! Author :
! ------
!    Josue Bock


! Modifications :
! -------------
  ! 30-11-2016   Josue Bock   first version of this SR
  ! 04-11-2017   Josue Bock   Fortran90
  ! 18-02-2019   Josue Bock   update error messages: use file unit jpfunerr for long messages,
  !                              and call abortM to stop the program


! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY: &
! Imported Routines:
       abortM

  USE file_unit, ONLY: &
! Imported Parameters:
       jpfunerr, jpfunout

  USE gas_common, ONLY: &
! Imported Parameters:
       j1,              &
       j5,              &
! Imported Array Variables with intent (in):
       gas_name,        &
       gas_mass,        &
       rad_name,        &
       rad_mass

  USE mod_nuc, ONLY: &
! Imported Parameters:
       nvap,         &
! Imported Scalar Variables with intent (out):
       ind_H2SO4,    &
       ind_NH3,      &
       ind_OIO,      &
! Imported Array Variables with intent (out):
       ical,         &
       ivap,         &
       m_vap,        &
       concsat,      &
       nuc_name

  use precision, only : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  logical, intent(in) :: Napari, Lovejoy, iod

! Local scalars:
  integer :: jvap, jspec

! == End of declarations =======================================================

! ==============================================================================
!        USER DEFINED LIST OF GASES HAS TO BE DEFINED BELOW
! ==============================================================================

!   --- Prescribe "initial" concentration(s) (in ppt) in gas_species.csv ---
!   --- also, NUCV (i.e. hypothetical unreactive condensible vapor, gas 75) can be used ---
!              ivap(i): condensible vapors, please adjust kind and number!
!              ical(i): radical(s3) =1, no radical(s1) =0, please adjust!
!   --- saturation vapor density in [molec/cm3] (=0. for non-volatile vapors), please adjust!
!       note: ideally, non-zero concsat should be similar in magnitude to avoid errors!

  jvap = 0

! OIO
  jvap = jvap+1
  nuc_name(jvap) = 'OIO'
  ical(jvap) = 1          ! radical
  concsat(jvap) = 0._dp

!! I2O2
!  jvap = jvap+1
!  nuc_name(jvap) = 'I2O2'
!  ical(jvap) = 0          ! non radical
!  concsat(jvap) = 0._dp

!! IO
!  jvap = jvap+1
!  nuc_name(jvap) = 'IO'
!  ical(jvap) = 1          ! radical
!  concsat(jvap) = 0._dp

!! HOI
!  jvap = jvap+1
!  nuc_name(jvap) = 'HOI'
!  ical(jvap) = 0          ! non radical
!  concsat(jvap) = 0._dp

! H2SO4
!  jvap = jvap+1
!  nuc_name(jvap) = 'H2SO4'
!  ical(jvap) = 0          ! non radical
!  concsat(jvap) = 0._dp

  if (jvap /= nvap) then
     write (jpfunerr,'(a)')'Error in SR nuc_init:'
     write (jpfunerr,'(a)')'  the number of condensible vapours (nvap = ',nvap,')'
     write (jpfunerr,'(a)')'  differs from the number of species names (jvap = ',jvap,').'
     write (jpfunerr,'(a)')'  Please adjust.'
     call abortM ('Error 1 SR nuc_init')
  end if

! ==============================================================================
!   Retrieve the species names from gas (radical or non radical) name lists
! ==============================================================================
  write (jpfunout,*)'Proceed to nuc initialisation: retrieve the species indexes'

  do jvap = 1,nvap
     jspec = 1 ! initialise index

     write (jpfunout,*)'Search for ',nuc_name(jvap)

     if (ical(jvap) == 0) then
        do while ( trim(gas_name(jspec)) /= trim(nuc_name(jvap)) )
           jspec = jspec+1
           if (jspec > j1) then
              write (jpfunerr,'(a)')'Error in SR nuc_init:'
              write (jpfunerr,'(a)')'  ',trim(nuc_name(jvap)),' not found in the list of gas names'
              call abortM ('Error 2 SR nuc_init')
           end if
        end do
        ivap(jvap) = jspec
        m_vap(jvap) = gas_mass(jspec)
        write (jpfunout,1010)trim(nuc_name(jvap)),' found in gas list with index ',jspec,', mass = ',gas_mass(jspec)

     else if (ical(jvap) == 1) then
        do while ( trim(rad_name(jspec)) /= trim(nuc_name(jvap)) )
           jspec = jspec+1
           if (jspec > j5) then
              write (jpfunerr,'(a)')'Error in SR nuc_init:'
              write (jpfunerr,'(a)')'  ',trim(nuc_name(jvap)),' not found in the list of radicals names'
              call abortM ('Error 3 SR nuc_init')
           end if
        end do
        ivap(jvap) = jspec
        m_vap(jvap) = rad_mass(jspec)
        write (jpfunout,1010)trim(nuc_name(jvap)),' found in rad list with index ',jspec,', mass = ',rad_mass(jspec)

     else
        write (jpfunerr,'(a)')'Error in SR nuc_init:'
        write (jpfunerr,'(a)')'  ical must be 0 or 1, ical = ',ical(jvap),jvap
        call abortM ('Error 4 SR nuc_init')
     end if

  end do

! ==============================================================================
!   Retrieve specific indexes depending on the options
! ==============================================================================
  if (Napari) then ! Search for H2SO4 and NH3 indexes
     jspec = 1              ! initialise index
     do while ( trim(gas_name(jspec)) /= 'H2SO4' )
        jspec = jspec+1
        if (jspec > j1) then
           write (jpfunerr,'(a)')'Error in SR nuc_init:'
           write (jpfunerr,'(a)')'  H2SO4 not found in the list of gas names'
           call abortM ('Error 5 SR nuc_init')
        end if
     end do
     ind_H2SO4 = jspec
     write (jpfunout,*)'Napari: H2SO4 found with index ',jspec

     jspec = 1              ! initialise index
     do while ( trim(gas_name(jspec)) /= 'NH3' )
        jspec = jspec+1
        if (jspec > j1) then
           write (jpfunerr,'(a)')'Error in SR nuc_init:'
           write (jpfunerr,'(a)')'  NH3 not found in the list of gas names'
           call abortM ('Error 6 SR nuc_init')
        end if
     end do
     ind_NH3 = jspec
     write (jpfunout,*)'Napari: NH3 found with index ',jspec
  else
     ind_H2SO4 = 0
     ind_NH3   = 0
  end if

  if (Lovejoy) then ! Search for OIO index
     if (.not.iod) then
        write (jpfunerr,'(a)')'Error in SR nuc_init: inconsistency in the nucleation module'
        write (jpfunerr,'(a)')'  Lovejoy scheme require OIO, thus iod must be true'
        call abortM ('Error 7 SR nuc_init')
     end if

     jspec = 1              ! initialise index
     do while ( trim(rad_name(jspec)) /= 'OIO' )
        jspec = jspec+1
        if (jspec > j5) then
           write (jpfunerr,'(a)')'Error in SR nuc_init:'
           write (jpfunerr,'(a)')'  OIO not found in the list of gas names'
           call abortM ('Error 8 SR nuc_init')
        end if
     end do
     ind_OIO = jspec
     write (jpfunout,*)'Lovejoy: OIO found with index ',jspec
  else
     ind_OIO = 0
  end if

!  print*,'ind_H2SO4, ind_NH3, ind_OIO',ind_H2SO4, ind_NH3, ind_OIO
!  write(*,1011)(jvap,ivap(jvap),ical(jvap),m_vap(jvap),jvap=1,nvap)
1010 format(a,a,i4,a,f6.3)
! 1011 format(3i4,f6.3)


end subroutine nuc_init
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine appnucl2 (dt,both)
!
! Description :
! -----------
!     help routine for apparent nucleation if apparent nucleation
!     has to be called twice because of two different real nucleation mechanisms
!


! Authors :
! -------
!    Susanne Pechtl (b. Marquart) (original code in Mistra v7.4.0)
!    Roland von Glasow


! Modifications :
! -------------
  !    Jul-2016   Josue Bock   Header, cleaning (unused CBs or variables, commented tests)
  !                            Use module for parameters
  !
  ! 04-Nov-2017   Josue Bock   Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       n

  use precision, only : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: dt
  logical, intent(in) :: both

! Local scalars:
  integer :: jz

! Local arrays:
  real (kind=dp) :: xn_app1(n), dnucv1(n), grorate1(n), grorate2(n), concnuc1(n), concnuc2(n)

! Common blocks:
  common /nuclapp/ xn_app(n), xn_apacc(n), xv_apacc(n),bn_ges(n), &
                   bd_mean(n), dnucv(n), grorate(n), concnuc(n)
  real (kind=dp) :: xn_app, xn_apacc, xv_apacc, bn_ges, &
                    bd_mean, dnucv, grorate, concnuc

! == End of declarations =======================================================


  call appnucl (dt,.true.,.false.,both)

  do jz = 2, n-1
     xn_app1(jz)  = xn_app(jz)
     dnucv1(jz)   = dnucv(jz)
     grorate1(jz) = grorate(jz)
     concnuc1(jz) = concnuc(jz)
  enddo

  call appnucl (dt,.false.,.true.,both)

  do jz = 2, n-1
     xn_app(jz)   = xn_app(jz) + xn_app1(jz)
     dnucv(jz)    = dnucv(jz)  + dnucv1(jz)
     grorate2(jz) = grorate(jz)
     grorate(jz)  = (grorate1(jz) + grorate2(jz)) / 2._dp
     concnuc2(jz) = concnuc(jz)
     if (concnuc2(jz).ge.concnuc1(jz)) then
        if ((xn_app(jz)-xn_app1(jz)).gt.0.01_dp) then
           concnuc(jz) = concnuc2(jz)*xn_app(jz)/(xn_app(jz)-xn_app1(jz))
        else
           concnuc(jz) = concnuc2(jz)
        endif
     else
        concnuc(jz) = concnuc1(jz) * xn_app(jz)/xn_app1(jz)
     endif
  enddo

end subroutine appnucl2
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine appnucl (dt,Napari,Lovejoy,both)
!
! Description :
! -----------
!     Calculation of the "apparent" nucleation rate of non-valotile vapors
!     after Kerminen and Kulmala (2002) @149.
!     With amendments for semivolatile vapors after Kerminen et al. (2004) @429

! Method :
! ------
!     assumptions of Kerminen and Kulmala (2002):
!     (1) only important sink for nuclei is coagulation to
!     larger pre-existing particles (i.e., no self-coagulation!)
!     (2) nuclei grow by condensation at a constant rate
!     (3) pre-existing particle population remains unchanged during nuclei growth
!     main subsequent restrictions:
!     total nuclei number concentration must remain sufficiently low
!     to prevent effective self-coagulation (< 10^5-10^6 nuclei/cm3)

!     modifications of Kerminen et al. (2004):
!     First, nuclei grow at a constant rate by condensation of non-volatile vapors.
!     After the nucleus has reached a critical size, it grows by non-volatile
!     and semi-volatile vapors (i.e. two different growth rates).
!
!     switches that must be set in str.f:
!       Napari:  optional coupling with ternary H2SO4-H2O-NH3 nucleation
!       Lovejoy: optional coupling with homogeneous OIO nucleation
!       ifeed:   optional feedback of nucleation with microphysics/chemistry
!              ifeed = 0 !no feedback to background particles
!              ifeed = 1 !with feedback to background particles
!              ifeed = 2 !feedback to background particles, but particles do not grow
!                         out of the 1st microphysical bin and chemistry does not "see"
!                         the 1st microphysical bin (i.e. no feedback nucl.->chemistry)
!
!     if chosen (Napari = .true.):
!     "Real" ternary H2SO4-H2O-NH3 nucleation rate and cluster size can be calculated
!     after Napari et al. (2002) @116
!     validity ranges: temp = 240 - 300 K
!                      RH = 0.05 - 0.95
!                      NH3 = 0.1 - 100 ppt
!                      H2SO4 = 1e4 - 1e9 molec./cm3
!                      J_real = 1e-5 - 1e6 /(cm3s)
!     if chosen (Lovejoy = .true.):
!     "Real" homogeneous OIO nucleation rate from data of Lovejoy/Burkholder
!
!     Output: J_real und d_nucini (= Input fuer APPNUCL)
!
!     to adjust in nuc.f: nvap, ivap, ical, m_vap, conc, concsat, (fcs), (xnue),
!                           d_nucini, J_real, uppest nucleation level (k-loop!)
!
!


! Authors :
! -------
!    Roland von Glasow (original code in Mistra v7.3.2)
!    Susanne Pechtl (b. Marquart)


! Modifications :
! -------------
  !    Jul-2016   Josue Bock   Header, cleaning
  !                            Use module for parameters
  !                            Explicit declaration of Knnuc and Kncrit
  !                              (problem with implicit typing, names starting with "K")
  !
  !    Aug-2016   Josue Bock   Remaining missing declarations (with comments) and implicit none
  !
  !    Nov-2016   Josue Bock   BUGFIX: transfer gas ==> liq was nested in the if (ical(jv))
  !                                    ... else construct, thus H2SO4 was not transferred
  !                                    properly to the liquid phase. Only OIO was ok.
  !
  ! 04-Nov-2017   Josue Bock   Fortran90
  !                            changed several test (... .ne. 0.) => (... .gt. 0.)
  !                            changed "over" (real) => llnovflw
  !
  ! 28-Apr-2021   Josue Bock   use rho3 (in module constants) to define ro_nuc
  !                            turn hard-coded "2000" into conversion factors zconv_r2d / d2r

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
       coutdir,      &
       ifeed              ! ifeed = 0/1/2: feedback with background particles no/yes/partly


  USE constants, ONLY : &
       conv1,           & ! multiply by conv1 to get cm^3(air)/mlc --> m^3(air)/mol
       pi,              &
       r1,              &
       ro_nuc => rho3,  & ! nuclei density [kg/m3] (same as assumed aesosol density rho3)
       rhow

  USE file_unit, ONLY : &
       jpfunnuc0

  USE gas_common, ONLY: &
! Imported Parameters:
       j1,              &
       j5,              &
! Imported Array Variables with intent (inout):
       s1,              &
       s3

  USE global_params, ONLY : &
! Imported Parameters:
       j2,                  &
       j3,                  &
       j6,                  &
       nkc,                 &
       nka,                 &
       nkt,                 &
       n

  USE mod_nuc, ONLY: &
! Imported Parameters:
       nvap,         &
       ical,         &
       ivap,         &
       nuc_name,     &
       m_vap,        &    ! Molar weight of nucleation vapor (e.g. NUCV=OIO) [kg/mol]
       concsat            ! saturation vapor density of condensible vapor [/cm3]
                          !  (measure for saturation vapor pressure over particle surface)

  use precision, only : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: dt
  logical, intent(in) :: Napari        ! Napari=true: J_real and d_nucini are calculated from ternary nucleation
  logical, intent(in) :: Lovejoy       ! Lovejoy=true: J_real and d_nucini are calculated from OIO nucleation
  logical, intent(in) :: both          ! true if both Napari and Lovejoy are true (needed)

! Local parameters:
  ! optimisation: define parameters that will be computed only once
  real (kind=dp), parameter :: zrho_frac = ro_nuc / rhow      ! = rho3 / rhow
  real (kind=dp), parameter :: zconv_r2d = 2000._dp           ! convert radius, in micro meter, to diameter, in nm
  real (kind=dp), parameter :: zconv_d2r = 1._dp / zconv_r2d  ! convert diameter, in nm, to radius, in micro meter

!     --- Local Variables ---
  integer :: icount,jts
  integer :: ia,jt,jtt                   ! loop indexes, 2D particle grid
  integer :: jv                          ! loop index, number of condensible vapour
  integer :: jz                          ! loop index, vertical (bottom-up)
  integer :: isum, mult

  real (kind=dp) :: alphaa !accomodation coefficient of condensing vapor
  real (kind=dp) :: betanuc, betacrit
  real (kind=dp) :: temp, press !Temperature [K] and Pressure [Pa]
  real (kind=dp) :: rh     ! relative humidity in layer k
  real (kind=dp) :: zdp(nkt)  !Particle diameter [nm]
  real (kind=dp) :: zdpmin    !Particle diameter of smallest dry size bin [nm]
  real (kind=dp) :: zdpmint   !Particle diameter of according total size bin [nm]
  real (kind=dp) :: Np(nkt)  !Particle number [/cm3]
  real (kind=dp) :: Nges     !total particle number [/cm3] (calculated)
  real (kind=dp) :: lambda   !Mean free path [m]
  real (kind=dp) :: v_mean(nvap) !Mean molecular speed of vapor [m/s]
  real (kind=dp) :: d_nucini !initial particle diameter [nm]; from subroutine ternucl
  real (kind=dp) :: J_real   !real nucleation rate [/(cm3*s)]; from subroutine ternucl
  real (kind=dp) :: d_mean   !Number mean diameter [nm]; (calculated)
  real (kind=dp) :: conc(nvap)    !concentration of condensible vapor [/cm3]
  real (kind=dp) :: conc_nuc
  real (kind=dp) :: concsemi      !sum of all semi-volatile vapors
  real (kind=dp) :: concsatmean   !mean saturation vapor density for semi-volatile vapors
  real (kind=dp) :: Scrit !Ratio of supersaturation of semi-volatile vapors
  real (kind=dp) :: dcrit !Critical nuclei diameter from which particle grows also
                          !  by condensation of semivolatile vapors [nm]
  real (kind=dp) :: Kn(nkt), beta(nkt), cs,gr,grs,gamma,eta,etages !calculated
  real (kind=dp) :: J_app    !apparent nucleation rate [/(cm3*s)] (main output)
  real (kind=dp) :: Knnuc, Kncrit
  real (kind=dp) :: factor
  real (kind=dp) :: a0, a0mn, b0, b0mn
  real (kind=dp) :: deltax, facul, gg, rmin, summ
  real (kind=dp) :: rg !equilibrium total aerosol radius
!     --- necessary for MISTRA mass bilances: ---
  real (kind=dp) :: rel_vap(nvap) !rel. contribution of vapor to nuclei growth
  real (kind=dp) :: m_vapmean     !mean molar weight of vapor
  real (kind=dp) :: s1old(j1,n), s3old(j5,n)   ![mol/m3]
  real (kind=dp) :: soldsum(n) ![mol/m3]

  logical :: llnovflw ! no risk of overflow

! Common blocks:
  common /backpart/ partd(nkt,n), partN(nkt,n), partsa(n)
  real (kind=dp) :: partd, partN, partsa

  common /blck01/ am3(n),cm3(n)
  real (kind=dp) :: am3, cm3

  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real (kind=dp) :: theta, thetl, t, talt, p, rho

  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real (kind=dp) :: xm1, xm2, feu, dfddt, xm1a

  common /nuclapp/ xn_app(n), xn_apacc(n), xv_apacc(n),bn_ges(n), &
                   bd_mean(n), dnucv(n), grorate(n), concnuc(n)
  real (kind=dp) :: xn_app, xn_apacc, xv_apacc, bn_ges, &
                    bd_mean, dnucv, grorate, concnuc

  common /nucl2/ Jn, dc      ! input from real nucleation
  real (kind=dp) :: Jn         ! nucleation rate [1/(cm3 s)]
  real (kind=dp) :: dc         ! diameter of cluster [nm]

  common /nucl3/ Jnio, dcio  ! input from real nucleation
  real (kind=dp) :: Jnio       ! nucleation rate [1/(cm3 s)]
  real (kind=dp) :: dcio       ! diameter of cluster [nm]

  common /nucsum/ ssum(n)
  real (kind=dp) :: ssum

! External function:
  real (kind=dp), external :: rgl

! == End of declarations =======================================================


  open(unit=jpfunnuc0,file=trim(coutdir)//'nuc.out',status='unknown',form='formatted')

  zdpmin = rn(1) * zconv_r2d ! rn is the radius in micro m, dp is the diameter in nm

!     --- alphaa: for values for some species see alpha(jz,spec) in subroutine st_coeff_t ---
!         but: in the current version only 1 alphaa can be used for all species together
  alphaa = 1._dp         !mass absorption coefficient of condensing vapor


  do jz=2,n-1  !Loop over vertical levels (does not necessariliy have to extend to n-1!)

     temp = t(jz)
     press = p(jz)
     rh = feu(jz)
     if (rh.ge.1._dp) rh = 0.999_dp !for equilibrium growth (see function rgl)

     ! Initialise concentration from non radical (s1) or radical (s3) gas concentration
     do jv=1,nvap
        if(ical(jv)==0) then
           conc(jv)=s1(ivap(jv),jz) * conv1
        else
           conc(jv)=s3(ivap(jv),jz) * conv1
        end if
     end do

!     --- 1D size distribution of background particles wrt total droplet diameter ---
     write (jpfunnuc0,1070)
     partsa(jz) = 0._dp
     do jt = 1,nkt
        zdp(jt) = zconv_r2d * rq(jt,1) !first dry aerosol class defines size bins for 1D distribution
        Np(jt)  = 0._dp
        do ia = 1, nka
           if (rn(ia).gt.rw(jt,1)) goto 2001
           do jtt = 1,nkt
              if (rq(jtt,ia).le.rw(jt,1)) then
                 if (jt.gt.1) then
                    if (rq(jtt,ia).gt.rw(jt-1,1)) Np(jt) = Np(jt) + ff(jtt,ia,jz)

                 !else if ((jt.eq.1).and.(rq(jtt+1,ia).gt.rw(jt,1))) then     ! jjb index out of bounds if jtt = nkt
                 else if ((jt.eq.1).and.(jtt.lt.nkt).and. &                   ! jjb quick fix, CHECK!
     &                    (rq(jtt+1,ia).gt.rw(jt,1))) then
                    Np(jt) = Np(jt) + ff(jtt,ia,jz)
                 endif
              else
                 goto 2002
              endif
           enddo
2002       continue
        enddo
2001    continue
        write (jpfunnuc0,1080) jz, jt, zdp(jt), Np(jt)
!         -- background particle surface area --
        partsa(jz) = partsa(jz) + Np(jt) * zdp(jt)**2 *3.1416* 1.d-6 !in um2/cm3
!         -- 1D particle distribution --
        partN(jt,jz) = Np(jt)
        partd(jt,jz) = zdp(jt)
     enddo  !jt
1070 format ('jz  jt  zdp(jt)   Np(jt)')
1080 format (2I4,2F16.5)

!     --- lambda = freep(jz) in kpp.f ---
     lambda = 2.28e-5_dp * temp / press

!     --- Initialize ---
     cs = 0._dp
     gr = 0._dp
     grs = 0._dp

     do jv = 1,nvap
!     --- vmean = sqrt((8*R_gas*T/(M*pi)) ---
!     --- sqrt(8*R_gas/pi)=4.60138 ---
        v_mean(jv) = sqrt(temp/m_vap(jv)) * 4.60138_dp
     enddo

! Note: both nucleation routines are called any time (for output reasons),
!       but nucleation is "physically seen" only if respective switches are .true.
!       I.e., if a switch is .false., nuclei are not seen by the model and
!       concentrations of nucleating species do not change.

!     --- ternary H2SO4-H2O-NH3 nucleation routine: ---
     if ((both).and.(.not.Napari).and.(Lovejoy)) then
!       do not call ternucl
     else
        call ternucl(dt,jz,Napari)
     endif

!     --- homogeneous OIO nucleation routine: ---
     if ((both).and.(Napari).and.(.not.Lovejoy)) then
!       do not call oionucl
     else
        call oionucl(dt,jz,Lovejoy)
     endif

     if (Napari) then
        d_nucini = dc
        J_real = Jn
     else if (Lovejoy) then
        d_nucini = dcio
        J_real = Jnio
     else !prescribe values for clusters
        d_nucini = 1.0_dp       ![nm]           please adjust!
        J_real = 1000._dp       ![/(cm3*s)]     please adjust!
     endif

!     --- semi-volatile vapors ---
!     --- concsemi can be used as switch for the existence of semi-volatile vapors ---
     concsemi = 0._dp !sum of all semi-volatile vapors
     concsatmean = 0._dp !mean saturation vapor density for semi-volatile vapors
     icount = 0
     do jv=1, nvap !mean values for semi-volatile vapors
        if (concsat(jv).gt.0._dp) then
           concsemi = concsemi + conc(jv)
           concsatmean = concsatmean + concsat(jv)
           icount = icount + 1
        endif
     enddo
     if (icount .ge. 1) concsatmean = concsatmean / real(icount, dp)
!     --- supersauration ratio and critical diameter for semi-volatile vapors after @429 ---
     if (concsatmean.gt.0._dp)  then
        Scrit = concsemi/concsatmean
        if (Scrit.le.1) then !no growth for subsaturation
           do jv=1, nvap
              if (concsat(jv).gt.0._dp) conc(jv)=0._dp
           enddo
           dcrit = zdpmin
        else
           dcrit = ( 6.49_dp - 0.01556_dp*temp + 0.039_dp*log(Scrit) ) / &
                    ( 1._dp - 0.002_dp*temp + 0.174_dp*log(Scrit) )
           if (dcrit.le.d_nucini) concsemi = 0._dp !semi-volatiles are defined as non-volatile
        endif
     endif

     call dmean(nkt,zdp,Np,d_mean,Nges)

     do jt = 1,nkt
!     --- Knudsen number ---
        Kn(jt) = 2._dp * 1.e9_dp * lambda / zdp(jt)
!     --- transition correction for condensational mass flux (Fuchs & Sutugin, 1971) ---
        beta(jt) = (1._dp + Kn(jt)) / (1._dp + 0.377_dp * Kn(jt) + 1.33_dp * Kn(jt) * (1._dp + Kn(jt)) / alphaa)
!     --- Condensation sink (condensation of vapor on pre-existing particles) ---
        cs = cs + 0.5_dp * zdp(jt) * 1.e-7_dp * beta(jt) * Np(jt) ![1/cm2]
     enddo

!     --- growth rate of nuclei: change of diameter with time (eq.(20) of @149) ---
!     --- gr: growth rate for all volatile vapors;
!     --- grs: growth rate for all semi-volatile vapors
     Knnuc = 2._dp * 1.e9_dp * lambda / d_nucini
     betanuc = (1._dp + Knnuc) / (1._dp + 0.377_dp * Knnuc + 1.33_dp * Knnuc * (1._dp + Knnuc) / alphaa)
     if (concsemi.gt.0._dp) then !semi-volatiles
        Kncrit = 2._dp * 1.e9_dp * lambda / dcrit
        betacrit = (1._dp + Kncrit) / (1._dp + 0.377_dp * Kncrit + 1.33_dp * Kncrit * (1._dp + Kncrit) / alphaa)
     endif
     do jv = 1,nvap
        if ((concsat(jv).gt.0._dp).and.(concsemi.gt.0._dp)) then !semi-volatiles
           grs = grs + v_mean(jv) * m_vap(jv) * (conc(jv)-concsat(jv))
        else !non-volatiles
           gr = gr + v_mean(jv) * m_vap(jv) * conc(jv)
        end if
     enddo
     gr= gr * 7969.45_dp * lambda * betanuc / d_nucini / ro_nuc   ![nm/h], non-volatiles
     if (concsemi.gt.0._dp) &
          grs= grs * 7969.45_dp * lambda * betacrit / dcrit / ro_nuc   ![nm/h], semi-volatiles

!     ---  relative contribution of vapor i to new aerosol mass ---
!          (calculated according to definition of gr and grs)
!          factor determines mass fraction of semi-volatile vapors
     m_vapmean = 0._dp
     if (concsemi.gt.0._dp) then
        factor = grs * (zdpmin**3-dcrit**3) / ( grs * (zdpmin**3-dcrit**3) + gr * (zdpmin**3-d_nucini**3) )
     else
        factor = 0._dp
     endif
     do  jv = 1,nvap
        if ((concsat(jv).gt.0._dp).and.(concsemi.gt.0._dp)) then !semi-volatiles
           if (grs .gt. 1.e-2_dp) then
              rel_vap(jv) = (v_mean(jv) * m_vap(jv)* (conc(jv)-concsat(jv))) &
                         / (grs *ro_nuc *dcrit /betacrit /lambda /7969.45_dp) * factor
           else
              rel_vap(jv) = 1._dp / real(nvap,dp)
           endif
        else
           if (gr .gt. 1.e-2_dp) then
              rel_vap(jv) = (v_mean(jv) * m_vap(jv)* conc(jv)) &
                     / (gr *ro_nuc *d_nucini /betanuc /lambda /7969.45_dp) * (1._dp - factor)
           else
              rel_vap(jv) = 1._dp / real(nvap,dp)
           endif
        endif
        m_vapmean = m_vapmean + rel_vap(jv) * m_vap(jv)
     enddo

!     --- calculate equilibrium with ambient water vapor (RH) ---
!     parameters a0m, b0m of koehler curve of subroutine subkon:
!     sr(ia,jt)=exp(a0m/(r(ia,jt)*t)-b0m(ia)*en(ia)/ew(jt)))
!     a0m see p. 141 pruppacher and klett a0m=2 sigma/(r1*t*rhow*a)
!     152200= 2 sigma*10**6 with sigma is surface tension = 76.1*10**-3
!     see pruppacher and klett p. 104
!     a0mn = a0m=152200./(r1*rhow)
!     b0mn = b0m=fcs*xnue*xmol2/xmol3: fcs(ia) fraction of soluble species
!     xnue number of ions; xmol2 (xmol3) mol masses of water (aerosol)
     a0mn = 152200._dp / (r1 * ro_nuc) ! jjb mistake here? Should this be rhow instead of ro_nuc?
     b0mn = 1._dp * 1._dp * 0.018_dp/m_vapmean  !fcs=1 and xnue=1 are assumed -> please adjust!
     a0 = a0mn/temp
!       ---  b0=b0m*rho3/rhow; rho3=2000; rhow=1000
     b0 = b0mn*zrho_frac ! jjb: not sure zrho_frac is correct: a0mn is computed with ro_nuc instead of rhow, shouldn't
                         !      that be the same here? I.e. rho3/ro_nuc?
     rmin = zdpmin * zconv_d2r
     rg = rgl(rmin,a0,b0,rh) !equilibrium total aerosol radius
     do jt = 1,nkt !find correct droplet size bin jts to smallest dry bin
        if (rw(jt,1) .ge. rg) then
           jts = jt
           goto 2003
        endif
     enddo
2003 continue
     zdpmint = rq(jts,1) * zconv_r2d
!     --- Update growth rate (including equilibrium with ambient RH) after Kerminen (pers. comm.)---
     gr = gr * zdpmint/zdpmin
     grs = grs * zdpmint/zdpmin

!     --- semi-empirical proportionality factor
!     (describing coagulation with larger pre-existing particles)
     gamma = 2300._dp *(d_nucini/1.d0)**0.2_dp *(zdpmint/3._dp)**0.075_dp &
          * (d_mean/150._dp)**0.048_dp * (ro_nuc/1000._dp)**(-0.33_dp) &
          * (temp/293._dp)**(-0.75_dp) ![nm2cm2/h]

!     --- Apparent nucleation rate for above calculated size bin of total diameter ---
     if ((gr .gt. 1.e-2_dp).and.(J_real .gt. 0.01_dp)) then
        eta = gamma * cs / gr   ![nm]
        etages = gamma * cs / (gr + grs)  ![nm]
        if (concsemi.gt.0._dp) then
           J_app = J_real * exp (etages/zdpmint + (eta-etages)/dcrit - eta/d_nucini) ![molecules/(cm3*s)]
        else
           J_app = J_real * exp (eta/zdpmint - eta/d_nucini)
        endif
!         --- diagnose total nuclei concentration < zdpmint ---
!             in order to check the validity range (19) in @149
!             conc_nunc should by < 1.e6 nuclei/cm3
!             (exactly valid only if no semi-volatile vapors are present)
        SUMM = 0._dp
        do isum = 1,20
           facul = 1._dp
           mult = isum
           do while (mult.GE.1)  !Calculation of k!
              facul = facul * mult
              mult = mult - 1
           enddo !while
           SUMM = SUMM + ( eta**(isum+1) * (d_nucini**(-isum) - zdpmint**(-isum)) ) / isum/facul
        enddo
! avoid floating overflow (exp(700) =(approx) 1.d308 which is highest number on manolito):
! (this is output only)
        llnovflw=.true.
        if (eta/zdpmint.gt.700._dp) llnovflw=.false.
        if (eta/d_nucini.gt.700._dp) llnovflw=.false.
        if (llnovflw) &
             GG = zdpmint*exp(eta/zdpmint) - d_nucini*exp(eta/d_nucini) &
                + eta*log(zdpmint/d_nucini) + SUMM
        conc_nuc = J_app/gr * exp (-eta/zdpmint) * GG * 3600._dp  ![nuclei/cm3]
     else
        eta      = 0._dp
        etages   = 0._dp
        J_app    = 0._dp
        conc_nuc = 0._dp
     endif

!     --- Update of variables ---
     xn_app(jz) = 0._dp
     bn_ges(jz) = Nges
     bd_mean(jz) = d_mean
     soldsum(jz) = 0._dp
     concnuc(jz) = 0._dp

     do jv = 1,nvap
        if (ical(jv).eq.0) then      ! non radical gas phase species
           s1old(ivap(jv),jz) = s1(ivap(jv),jz)
           soldsum(jz) = soldsum(jz) + s1old(ivap(jv),jz)
        else if (ical(jv).eq.1) then ! radical gas phase species
           s3old(ivap(jv),jz) = s3(ivap(jv),jz)
           soldsum(jz) = soldsum(jz) + s3old(ivap(jv),jz)
        endif
     enddo !jv

     ssum(jz) = soldsum(jz)

     if (J_app .gt. 0.1_dp) then
        xn_app(jz) = J_app  !nucleation rate of new "apparent" particles [/(cm3*s)]
        concnuc(jz) = conc_nuc ! 1D (real) nuclei concentration for output
        if (ifeed.ne.0) ff(jts,1,jz) = ff(jts,1,jz) + xn_app(jz)*dt
        deltax = xn_app(jz)*dt * pi/6._dp * (zdpmin**3-d_nucini**3) &
      &           * ro_nuc / m_vapmean * 1.e-21_dp    ![mol/m3(air)] !additional dry material
        ssum(jz) = 0._dp
        do jv = 1,nvap
           if (ical(jv).eq.0) then
              s1(ivap(jv),jz) = s1(ivap(jv),jz) - deltax*rel_vap(jv)
              if (s1(ivap(jv),jz) .lt. 0._dp) s1(ivap(jv),jz) = 0._dp
              ssum(jz) = ssum(jz) + s1(ivap(jv),jz)
           else
              s3(ivap(jv),jz) = s3(ivap(jv),jz) - deltax*rel_vap(jv)
              if (s3(ivap(jv),jz) .lt. 0._dp) s3(ivap(jv),jz) = 0._dp
              ssum(jz) = ssum(jz) + s3(ivap(jv),jz) !ssum = sum of nucleating vapors
           end if
! OIO goes into liquid phase (for mass conserving reasons only)
! note: other condensing vapors than OIO and H2SO4 are NOT conserved!
!    OIO goes to unreactive OIO in the liquid phase
           if (trim(nuc_name(jv)).eq.'OIO') &
                sl1(j2-j3+12,1,jz) = sl1(j2-j3+12,1,jz) + s3old(ivap(jv),jz) - s3(ivap(jv),jz)
!    H2SO4 goes to H2SO4 in the liquid phase (which will further equilibrate with HSO4- and SO4=)
           if (trim(nuc_name(jv)).eq.'H2SO4') &
                sl1(6,1,jz) = sl1(6,1,jz) + s1old(ivap(jv),jz) - s1(ivap(jv),jz)
        enddo !jv

        xn_apacc(jz) = xn_apacc(jz) + xn_app(jz)*dt                        ! accumulated number [/cm3]
        xv_apacc(jz) = xv_apacc(jz) + pi/6._dp * zdpmint**3 *xn_app(jz)*dt ! acculumated volume [nm3/cm3]
        write (jpfunnuc0,109)
        write (jpfunnuc0,110) jz,jts,xn_apacc(jz),xn_app(jz),ssum(jz)
     endif
     !dnucv(jz) = soldsum(jz) - ssum(jz)                      ! change in nucleating vapors (sum) [mol/m3]
     dnucv(jz) = (soldsum(jz) - ssum(jz)) / am3(jz) *1.e12_dp  ! change in nucleating vapors (sum) [ppt]
     grorate(jz) = gr                                       ! 1D growth rate for output (for non-volatiles)

!     --- Output (Test) ---
     write(jpfunnuc0,1003)
     do jv=1,nvap
        if (j_real.gt.0._dp) then
           write(jpfunnuc0,1004) jz,d_nucini,zdpmin,zdpmint,d_mean,Nges,temp, &
                press,rh,ro_nuc,m_vap(jv),conc(jv),cs,gr,eta,gamma, &
                J_real,J_app,(J_app/J_real)
        else
           write(jpfunnuc0,1004) jz,d_nucini,zdpmin,zdpmint,d_mean,Nges,temp, &
                press,rh,ro_nuc,m_vap(jv),conc(jv),cs,gr,eta,gamma, &
                J_real,J_app
        endif
     enddo

  enddo !Loop over vertical levels (jz)

1003 format ('jz,d_nucini(nm),zdpmin(nm),zdpmint(nm),d_mean(nm),Nges, ' &
      // 'temp(K),press(Pa),rh,ro_nuc(kg/m3),m_vap(jv)(kg/mol), '       &
      // 'conc(jv)(/cm3),cs(/cm2),gr(nm/h),eta(nm),gamma(nm2cm2/h), '   &
      // 'J_real(/cm3s),J_app(/cm3s),(J_app/J_real)')
 109  format ('jz,jts,xn_apacc(jz),xn_app(jz),ssum(jz)')
 1004 format (I4,10F14.3,E12.5,7F14.5)
 110  format (2i4,9d16.8)

  close (jpfunnuc0)

end subroutine appnucl
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine dmean (n,pdp,Np,d_mean,Nges)
!
! Description :
! -----------
!     --- This routine calculates the number mean diameter d_mean ---
!         and the total particle number Np_ges of the background particles
!


! Authors :
! -------
!    Roland von Glasow (original code in Mistra v7.3.2)
!    Susanne Pechtl (b. Marquart)


! Modifications :
! -------------
  !    Jul-2016   Josue Bock   Header, implicit none
  !
  ! 04-Nov-2017   Josue Bock   Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  use precision, only : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  integer, intent(in) :: n
! Array arguments with intent(in):
  real (kind=dp), intent(in) :: pdp(n), Np(n)
! Scalar arguments with intent(out):
  real (kind=dp), intent(out) :: d_mean, Nges
! Local scalars:
  integer :: i

! == End of declarations =======================================================

  d_mean = 0._dp
  Nges = 0._dp
  do i=1,n
     d_mean = d_mean + pdp(i)*Np(i)
     Nges = Nges + Np(i)
  enddo
  if (Nges.gt.0._dp) then
     d_mean = d_mean / Nges
  else
     d_mean = 0._dp
  endif

end subroutine dmean
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine ternucl (dt,kz,Napari)
!
! Description :
! -----------
!    calculation of nucleation rate and cluster size
!    using the parameterization of Napari et al., JGR, 107, 4381,
!    doi:10.1029/2002JD002132, 2002
!    validity range: h2so4: 1.d4 - 1.d9 molec/cm3
!                    nh3: 0.1 - 100 ppt
!                    J_n: 1.d-5 - 1.d6 /(cm3 s)


! Authors :
! -------
!    Roland von Glasow (original code in Mistra v7.1.1)
!    Susanne Pechtl (b. Marquart)


! Modifications :
! -------------
  !    Jul-2016   Josue Bock   Header, implicit none
  !                            Use modules for constant and parameters
  !
  ! 02-Mar-2017   Josue Bock   BUGFIX: missing initialisation of nn and nh
  !
  ! 04-Nov-2017   Josue Bock   Fortran90

! == End of header =============================================================

  USE constants, ONLY : &
       conv1 ! multiply by conv1 to convert cm^3(air)/mlc --> m^3(air)/mol

  USE gas_common, ONLY: &
       s1

  USE file_unit, ONLY : &
       jpfunnuc0

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE mod_nuc, ONLY: &
! Imported Parameters:
       ind_H2SO4, &
       ind_NH3

  use precision, only : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: dt
  integer, intent(in) :: kz
  logical, intent(in) :: Napari

! Local scalars:
  real (kind=dp) :: &
       c_h2so4,     &
       c_nh3,       &
       nh,          &   ! number of H2SO4 molecules in critical cluster
       nn,          &   ! number of NH3 molecules in critical cluster
!       n_tot,nt,    &   ! total number of molecules in critical cluster
       rc,          &   ! radius of cluster [nm]
       rh,          &
       Temp

! Internal functions:
  real (kind=dp) :: &
       a,b,         &   ! arguments
       n_h2so4,     &   ! number of H2SO4 molecules in critical cluster
       n_nh3,       &   ! number of NH3 molecules in critical cluster
       r_clust          ! radius of cluster [nm]

! External function:
  real (kind=dp), external :: &
       J_nuc            ! nucleation rate [1/(cm3 s)]

! Common blocks:
  common /blck01/ am3(n),cm3(n)
  real (kind=dp) :: am3, cm3

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real (kind=dp) :: theta, thetl, t, talt, p, rho

  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real (kind=dp) :: xm1, xm2, feu, dfddt, xm1a

  common /nucl/ xn_new(n), xn_acc(n), xv_acc(n), dh2so4(n), dnh3(n)
  real (kind=dp) :: xn_new, xn_acc, xv_acc, dh2so4, dnh3

  common /nucl2/ Jn, dc   !SM: Jn and dc are used in subroutine appnucl
  real (kind=dp) :: &
       Jn,          &     ! nucleation rate [1/(cm3 s)]
       dc                 ! diameter of cluster [nm]
!- End of header ---------------------------------------------------------------

! functions for diagnostics of critical cluster; a = J_nuc, b=Temp
  n_h2so4 (a,b) = 38.1645_dp   + 0.774106_dp  *log(a)   + 2.98879e-3_dp*log(a)**2 &
               - 0.357605_dp*b - 3.66358e-3_dp*log(a)*b + 8.553e-4_dp  *b**2

  n_nh3 (a,b)   = 26.8982_dp   + 0.682905_dp  *log(a)   + 3.57521e-3_dp*log(a)**2 &
               - 0.265748_dp*b - 3.41895e-3_dp*log(a)*b + 6.73454e-4_dp*b**2

!  n_tot (a,b)   = 79.3484_dp   + 1.7384_dp    *log(a)   + 7.11403e-3_dp*log(a)**2 &
!               - 0.744993_dp*b - 8.20608e-3_dp*log(a)*b + 1.7855e-3_dp *b**2

  r_clust (a,b) = 0.141027_dp   - 1.22625e-3_dp*log(a)   - 7.82211e-6_dp*log(a)**2 &
              - 1.56727e-3_dp*b - 3.076e-5_dp  *log(a)*b + 1.08375e-5_dp*b**2

! convert units
  rh      = feu(kz)                            ! relative humidity [0..1]
  c_nh3   = s1(ind_NH3,kz) / am3(kz) *1.e12_dp ! NH3 in [pmol/mol]
  c_nh3   = min(100._dp,c_nh3)                 ! NH3 <= 100 ppt
  c_h2so4 = s1(ind_H2SO4,kz) * conv1           ! H2SO4 in [molec/cm3]
  !c_h2so4 = c_h2so4 * 1000._dp                ! artificial amplification for test
  write (jpfunnuc0,120)
  write (jpfunnuc0,110) kz,rh,c_nh3,c_h2so4
  Temp    = t(kz)                             ! temperature [K]
  Jn = 0._dp
  nn = 0._dp
  nh = 0._dp
  xn_new(kz) = 0._dp
  rc = 1._dp !initialize
  dc = 2._dp * rc
  if (c_h2so4.gt.1.e4_dp) then
! calculate the nucleation rate
     Jn=min(1.e6_dp,J_nuc(rh,c_nh3,c_h2so4,Temp))     ! nucleation rate [1/(cm3 s)]
     if (Jn .ge. 0.01_dp) then
        nh=n_h2so4(Jn,Temp) ! number of H2SO4 molecules in critical cluster
        nn=n_nh3(Jn,Temp)   ! number of NH3 molecules in critical cluster
        !nt=n_tot(Jn,Temp)   ! total number of molecules in critical cluster
        rc=r_clust(Jn,Temp) ! radius of cluster [nm]
        dc = 2._dp * rc     ! diameter of cluster [nm]
        xn_new(kz) = Jn      ! nucleation rate of new clusters [/(cm3*s)]
        if (Napari) then
           s1(ind_NH3,kz) = s1(ind_NH3,kz) - xn_new(kz) * dt * nn / conv1
           if (s1(ind_NH3,kz) .lt. 0._dp) s1(ind_NH3,kz) = 0._dp

           s1(ind_H2SO4,kz) = s1(ind_H2SO4,kz) - xn_new(kz) * dt * nh / conv1
           if (s1(ind_H2SO4,kz) .lt. 0._dp) s1(ind_H2SO4,kz) = 0._dp
        endif
        xn_acc(kz) = xn_acc(kz) + xn_new(kz)*dt                  ! accumulated number of new clusters [/cm3]
        xv_acc(kz) = xv_acc(kz) + 4.18879_dp*rc**3*xn_new(kz)*dt ! acculumated volume of new clusters [nm3/cm3]
        write (jpfunnuc0,109)
        write (jpfunnuc0,110) kz,dc,nn,nh,xn_acc(kz),xn_new(kz)
     endif
  endif
!  dnh3(kz)   = (c_nh3*am3(kz)/1.e12_dp) - s1(ind_NH3,kz)  ![mol/m3]
!  dh2so4(kz) = (c_h2so4/conv1) - s1(ind_H2SO4,kz)   ![mol/m3]
!  dnh3(kz)   = xn_new(kz)*dt*nn/conv1           ![mol/m3]
!  dh2so4(kz) = xn_new(kz)*dt*nh/conv1         ![mol/m3]
  dnh3(kz)   = xn_new(kz)*dt*nn/conv1 / am3(kz) *1.e12_dp ![pmol/mol]
  dh2so4(kz) = xn_new(kz)*dt*nh                     ![molec/cm3]

! save xnew to some output array
 109  format ('kz,dc,nn,nh,xn_acc(kz),xn_new(kz)')
 110  format (i4,9d16.8)
 120  format ('kz,rh,c_nh3[ppt],c_h2so4[/cm3]')

end subroutine ternucl

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine oionucl (dt,kz,Lovejoy)
!
! Description :
! -----------
!    calculation of nucleation rate from function
!    inferred from model simulation of homogeneous
!    OIO nucleation in the lab by Lovejoy/Burkholder (pers. comm.)
!
!    function is valid for nuclei with d=2.0nm and density=2g/cm3
!    (i.e. 1 nuclei has 34 OIO molecules)
!    validity range: 0.4 ppt < OIO < 22 ppt, extrapolation seems ok
!
!    it is assumed that only OIO contributes to nucleation
!


! Authors :
! -------
!    Susanne Pechtl (b. Marquart) (original code in Mistra v7.4.0)
!    Roland von Glasow


! Modifications :
! -------------
  !    Jul-2016   Josue Bock   Header
  !                            Use module for constant and parameters
  !                            Declarations (with comments) and implicit none
  !                            BUGFIX: renamed internal function J_nuc => J_nuc2
  !                                    which was conflicting with external fct J_nuc
  !
  ! 04-Nov-2017   Josue Bock   Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE constants, ONLY : &
       conv1 ! multiply by conv1 to get cm^3(air)/mlc --> m^3(air)/mol

  USE file_unit, ONLY : &
       jpfunnuc0

  USE gas_common, ONLY: &
       s3

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE mod_nuc, ONLY: &
       ind_OIO

  use precision, only : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: dt
  integer, intent(in) :: kz
  logical, intent(in) :: Lovejoy

! Local parameters:
  real (kind=dp), parameter :: n_oio = 34._dp   ! number of OIO molecules in critical cluster
  real (kind=dp), parameter :: rc = 1._dp       ! radius of cluster [nm]

! Local scalar:
  real (kind=dp) :: c_oio          ! OIO in [pmol/mol]
  real (kind=dp) :: Temp           ! temperature [K]

! Internal function:
  real (kind=dp) :: J_nuc2    ! (statement function)
  real (kind=dp) :: a, b      ! internal function arguments

! Common blocks:
  common /blck01/ am3(n),cm3(n)
  real (kind=dp) :: am3, cm3

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real (kind=dp) :: theta, thetl, t, talt, p, rho

  common /nuclio/ xn_newio(n), xn_accio(n), xv_accio(n), doio(n)
  real (kind=dp) :: xn_newio, xn_accio, xv_accio, doio

  common /nucl3/ Jnio, dcio
  real (kind=dp) :: Jnio         ! nucleation rate [1/(cm3 s)]
  real (kind=dp) :: dcio         ! diameter of cluster [nm]

! Internal function for diagnostics of nucleation rate; a = oio[ppt], b=temp[K]
  J_nuc2(a,b) = a**(0.030657_dp * b - 4.4471_dp) * exp(-0.30947_dp * b + 81.097_dp)

! == End of declarations =======================================================

! convert units
  c_oio   = s3(ind_OIO,kz) / am3(kz) *1.e12_dp
  write (jpfunnuc0,120)
  write (jpfunnuc0,110) kz,c_oio

  dcio = 2._dp * rc
  Temp = t(kz)
  Jnio = 0._dp
  xn_newio(kz) = 0._dp
  if (c_oio .gt. 0.01_dp) then
! calculate the nucleation rate
     Jnio=min(1.e4_dp,J_nuc2(c_oio,Temp))
     if (Jnio .ge. 0.01_dp) then
        xn_newio(kz) = Jnio    !nucleation rate of new clusters [/(cm3*s)]
        if (Lovejoy) then
           s3(ind_OIO,kz) = s3(ind_OIO,kz)-xn_newio(kz)*dt*n_oio/conv1
           if (s3(ind_OIO,kz) .lt. 0._dp) s3(ind_OIO,kz)=0._dp
        endif
        xn_accio(kz) = xn_accio(kz) + xn_newio(kz)*dt !accumulated number of new clusters [/cm3]
        xv_accio(kz) = xv_accio(kz) + 4.18879_dp*rc**3*xn_newio(kz)*dt ! acculumated volume of new clusters [nm3/cm3]
        write (jpfunnuc0,109)
        write (jpfunnuc0,110) kz,dcio,n_oio,xn_accio(kz),xn_newio(kz)
     endif
  endif

!  doio(kz) = (c_oio*am3(kz)/1.d12) - s3(ind_OIO,kz)          ![mol/m3]
!  doio(kz) = xn_newio(kz)*dt*n_oio/conv1                     ![mol/m3]
  doio(kz) = xn_newio(kz)*dt*n_oio/conv1 / am3(kz) *1.e12_dp ![pmol/mol]

! save xnew to some output array
 109  format ('kz,dcio,n_oio,xn_accio(kz),xn_newio(kz)')
 110  format (i4,9d16.8)
 120  format ('kz,c_oio[ppt]')

end subroutine oionucl

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function J_nuc (rH,nh3,h2so4,Temp)
!
! Description :
! -----------
!    nucleation rate
!    checked and compared to figs in #2158/@116: okay

! the parameterization is valid until for nucleation rates up to 10^6 1/(cm3 s), here
! it's used also for higher rates, so it's clearly an upper limit - magnitude of error
! not checked yet

! ln c  = ln [H2SO4], []=concentration
! ln S  = ln [NH3],   []=mixing ratio
! ln rH = ln rH
!


! Authors :
! -------
!    Susanne Pechtl (b. Marquart) (original code in Mistra v7.3.2)
!    Roland von Glasow


! Modifications :
! -------------
  !    Jul-2016   Josue Bock   Header, implicit none
  !                            Use module for parameters
  !
  ! 04-Nov-2017   Josue Bock   Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  use precision, only : &
! Imported Parameters:
       dp

  implicit none

! Function declaration
  real (kind=dp) :: J_nuc

! Function arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: rh,nh3,h2so4,Temp

! Local parameters:
  real (kind=dp), parameter :: fpd(4,20) = reshape( &
    (/-0.355297_dp,    -3.38448e+1_dp,  0.34536_dp,    -8.24007e-4_dp, &
       3.13735_dp,     -0.772861_dp,    5.61204e-3_dp, -9.74576e-6_dp, &
       1.90359e+1_dp,  -0.170957_dp,    4.79808e-4_dp, -4.14699e-7_dp, &
       1.07605_dp,      1.48932_dp,    -7.96052e-3_dp,  7.61229e-6_dp, &
       6.0916_dp,      -1.25378_dp,     9.39836e-3_dp, -1.74927e-5_dp, &
       0.31176_dp,      1.64009_dp,    -3.43852e-3_dp, -1.09753e-5_dp, &
      -2.00735e-2_dp,  -0.752115_dp,    5.25813e-3_dp, -8.98038e-6_dp, &
       0.165536_dp,     3.26623_dp,    -4.89703e-2_dp,  1.46967e-4_dp, &
       6.52645_dp,     -0.258002_dp,    1.43456e-3_dp, -2.02036e-6_dp, &
       3.68024_dp,     -0.204098_dp,    1.06259e-3_dp, -1.26560e-6_dp, &
      -6.6514e-2_dp,   -7.82382_dp,     1.22938e-2_dp,  6.18554e-5_dp, &
       0.65874_dp,      0.190542_dp,   -1.65718e-3_dp,  3.41744e-6_dp, &
       5.99321e-2_dp,   5.96475_dp,    -3.62432e-2_dp,  4.93337e-5_dp, &
      -0.732731_dp,    -1.84179e-2_dp,  1.47186e-4_dp, -2.37711e-7_dp, &
       0.728429_dp,     3.64736_dp,    -2.7422e-2_dp,   4.93478e-5_dp, &
       4.13016e+1_dp,  -0.35752_dp,     9.04383e-4_dp, -5.73788e-7_dp, &
      -0.160336_dp,     8.89881e-3_dp, -5.39514e-5_dp,  8.39522e-8_dp, &
       8.57868_dp,     -0.112358_dp,    4.72626e-4_dp, -6.48365e-7_dp, &
       5.301767e-2_dp, -1.98815_dp,     1.57827e-2_dp, -2.93564e-5_dp, &
      -2.32736_dp,      2.34646e-2_dp, -7.6519e-5_dp,   8.0459e-8_dp/), &
       shape(fpd))

! Local scalars:
  real (kind=dp) :: lnc, lnS, lnrH

! Internal function:
  real (kind=dp) :: fpol    ! (statement function)
  real (kind=dp) :: T
  integer :: n
  fpol(n,T) = fpd(1,n) + fpd(2,n)*T + fpd(3,n)*T**2 + fpd(4,n)*T**3

! == End of declarations =======================================================

  lnc  = log(h2so4)
  lnS  = log(nh3)
  lnrH = log(rH)

  J_nuc = exp(-84.7551_dp + fpol(1,Temp)/lnc + fpol(2,Temp)*lnc + &
       fpol(3,Temp)*lnc**2 + fpol(4,Temp)*lnS +                   &
       fpol(5,Temp)*lnS**2 + fpol(6,Temp)*rH +                    &
       fpol(7,Temp)*lnrH + fpol(8,Temp)*lnS/lnc +                 &
       fpol(9,Temp)*lnS*lnc + fpol(10,Temp)*rH*lnc +              &
       fpol(11,Temp)*rH/lnc + fpol(12,Temp)*rH*lnS +              &
       fpol(13,Temp)*lnrH/lnc + fpol(14,Temp)*lnrH*lnS +          &
       fpol(15,Temp)*lnS**2/lnc + fpol(16,Temp)*lnc*lnS**2 +      &
       fpol(17,Temp)*lnc**2*lnS +  fpol(18,Temp)*rH*lnS**2 +      &
       fpol(19,Temp)*rH*lnS/lnc + fpol(20,Temp)*lnc**2*lnS**2)

end function J_nuc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine nucout1
!
! Description:
! -----------
!    output nucleation: .asc output for gnu-plotting


! Authors :
! -------
!    Susanne Pechtl (b. Marquart) (original code in Mistra v7.3.2)
!    Roland von Glasow


! Modifications :
! -------------
  !    Aug-2016   Josue Bock   Header, implicit none
  !                            Use module for parameters
  !
  ! 04-Nov-2017   Josue Bock   Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
       coutdir

  USE constants, ONLY : &
       conv1 ! multiply by conv1 to get cm^3(air)/mlc --> m^3(air)/mol

  USE gas_common, ONLY: &
       s1,              & ! non radical gas concentration
       s3                 ! radical gas concentration

  USE file_unit, ONLY : &
       jpfunnuc1, jpfunnuc2, &
       jpfunnuc3, jpfunnuc4

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE mod_nuc, ONLY: &
       ind_H2SO4, &
       ind_NH3, &
       ind_OIO

  use precision, only : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  character (len=150) :: clpath  ! complete path to file

! Common blocks:
  common /blck01/ am3(n),cm3(n)
  real (kind=dp) :: am3,cm3

  common /nucl/ xn_new(n), xn_acc(n), xv_acc(n), dh2so4(n), dnh3(n)
  real (kind=dp) :: xn_new, xn_acc, xv_acc, dh2so4, dnh3

  common /nuclio/ xn_newio(n), xn_accio(n), xv_accio(n), doio(n)
  real (kind=dp) :: xn_newio, xn_accio, xv_accio, doio

  common /nuclapp/ xn_app(n), xn_apacc(n), xv_apacc(n),bn_ges(n), &
                   bd_mean(n), dnucv(n), grorate(n), concnuc(n)
  real (kind=dp) :: xn_app, xn_apacc, xv_apacc, bn_ges, &
                    bd_mean, dnucv, grorate, concnuc

  common /nucsum/ ssum(n)
  real (kind=dp) :: ssum

! == End of declarations =======================================================

  clpath=trim(coutdir)//'nuc+part.asc'
  open (unit=jpfunnuc1, file=clpath, status='unknown')
  write (jpfunnuc1,106) xn_app(2),xn_apacc(2),bd_mean(2),bn_ges(2),grorate(2),concnuc(2)
  close (jpfunnuc1)

  clpath=trim(coutdir)//'ternuc.asc'
  open (unit=jpfunnuc2, file=clpath, status='unknown')
  write (jpfunnuc2,104) xn_new(2),xn_acc(2),xn_newio(2),xn_accio(2)
  close (jpfunnuc2)

  clpath=trim(coutdir)//'nh3_h2so4.asc'
  open (unit=jpfunnuc3, file=clpath, status='unknown')
  write (jpfunnuc3,106) (s1(ind_NH3,2)/am3(2)*1.e12_dp),(s1(ind_H2SO4,2)*conv1), &
                 dnh3(2),dh2so4(2), &
                 (s3(ind_OIO,2)/am3(2)*1.e12_dp),doio(2)
  close (jpfunnuc3)

  clpath=trim(coutdir)//'NUCV.asc'
  open (unit=jpfunnuc4, file=clpath, status='unknown')
  write (jpfunnuc4,102) (ssum(2)/am3(2)*1.e12_dp),dnucv(2)
  close (jpfunnuc4)

 102  format (2d16.8)
 104  format (4d16.8)
 106  format (6d16.8)

end subroutine nucout1

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine nucout2
!
! Description:
! -----------
!    output nucleation


! Authors :
! -------
!    Susanne Pechtl (b. Marquart) (original code in Mistra v7.3.2)
!    Roland von Glasow


! Modifications :
! -------------
  !    Aug-2016   Josue Bock   Header, implicit none
  !                            Use module for parameters
  !
  ! 04-Nov-2017   Josue Bock   Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
       coutdir

  USE file_unit, ONLY : &
       jpfunnuc5

  USE global_params, ONLY : &
! Imported Parameters:
       n

  use precision, only : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  character (len=150) :: clpath  ! complete path to file
  integer :: jz

! Common blocks:
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday,lst,lmin,it,lcl,lct

  common /nucl/ xn_new(n), xn_acc(n), xv_acc(n), dh2so4(n), dnh3(n)
  real (kind=dp) :: xn_new, xn_acc, xv_acc, dh2so4, dnh3

  common /nuclapp/ xn_app(n), xn_apacc(n), xv_apacc(n),bn_ges(n), &
                   bd_mean(n), dnucv(n), grorate(n), concnuc(n)
  real (kind=dp) :: xn_app, xn_apacc, xv_apacc, bn_ges, &
                    bd_mean, dnucv, grorate, concnuc

  common /nuclio/ xn_newio(n), xn_accio(n), xv_accio(n), doio(n)
  real (kind=dp) :: xn_newio, xn_accio, xv_accio, doio

! == End of declarations =======================================================

  clpath=trim(coutdir)//'nucout2.out'
  open (unit=jpfunnuc5, file=clpath, status='unknown')

  write (jpfunnuc5,104) lday,lst,lmin
  do jz=1,n
     write (jpfunnuc5,105) jz,xn_new(jz),xn_acc(jz),xv_acc(jz)
  enddo
  write (jpfunnuc5,104) lday,lst,lmin
  do jz=1,n
     write (jpfunnuc5,105) jz,xn_newio(jz),xn_accio(jz),xv_accio(jz)
  enddo
  write (jpfunnuc5,104) lday,lst,lmin
  do jz=1,n
     write (jpfunnuc5,105) jz,xn_app(jz),xn_apacc(jz),xv_apacc(jz)
  enddo

  close (jpfunnuc5)

 104  format (3i4)
 105  format (i4,3d16.8)

end subroutine nucout2
