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
! nrad.f : radiation
! ------------------
! calculation of radiative fluxes for 18 spectral bands using the
! cumulative probability method of Fu (1991)
! ---------------------------------------------------------------
!
! contains the following subroutines and functions:
!     - nstrahl
!     - frr
!     - water
!     - gascon
!     - qopcon
!     - planck
!     - plkavg    (function)
!     - plancktab             (unused, commented)
!     - fst4      (function)  (unused, commented)
!     - gase
!     - qks
!     - qki
!     - qkio3
!     - qopo3s
!     - qoph2o
!     - qopch4
!     - qopno2
!     - qopo3i
!     - qophc
!     - tau
!     - kurzw
!     - langw
!     - jeanfr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine nstrahl
!
! Description :
! -----------
!    Main program of the radiation module:
!    Calls subroutines for calculation of radiative transfer equation for 18
!    spectral bands (6 solar and 12 ir) and associated cumulative probabilities.
!    Finally integrates results over the spectral intervals
!
!    Changes were done in the old part of the program (PIFM1) to upgrade it to 18 spectral
!    bands. This includes reading in more information in the routine 'initr'
!    with qext, qabs, and asymm. The upgrade then follows on to other variables
!    which are calculated from these three, namely bea baa and ga. There are also
!    changes in the arrays which store the radiation data as a function of wavelength.
!


! Author :
! ------
!    developed at the University of Mainz by W.-G. Panhans and others


! Main variables:
!     as  is the surface albedo for the 6 solar spectral regions
!     ee  is the emissivity of the ground
!
!  all radiative fluxes in W m**-2
!     f1  is the diffuse upward flux
!     fs1 is the diffuse upward flux (short wavelenght <=> solar bands)
!     fl1 is the diffuse upward flux (long wavelenght <=> IR bands)
!     f1f is the diffuse upward flux, cloud free case (short or long wavelength depending on context)
!     f1w is the diffuse upward flux, cloudy case (short or long wavelength depending on context)
!
!     f2  is the diffuse downward flux
!     fs2 is the diffuse downward flux (short wavelenght <=> solar bands)
!     fl2 is the diffuse downward flux (long wavelenght <=> IR bands)
!     f2f is the diffuse downward flux, cloud free case (short or long wavelength depending on context)
!     f2w is the diffuse downward flux, cloudy case (short or long wavelength depending on context)
!
!     fnseb is the radiative net flux in the solar spectrum
!
!     hr(i) is the radiative heating rate in Kelvin per second
!
!     s0 is solar constant


! Modifications :
! -------------
  !                       Original code, PIFM1
  !                       Upgraded from 6 to 18 spectral bands: PIFM2
  !
  !           Roland      Introduced specific outputs that are further used in the
  !           von Glasow  photolysis module: taer_s, taer_a, ga_pl
  !
  ! Jul-2015  Josue Bock  First check of the code using Forcheck (Mistra v7.3.3)
  !                       minor cleaning
  ! Mar-2016  Josue Bock  Transferred to Mistra v7.4.1, further cleaning
  ! Jul-2016  Josue Bock  BUGFIX: berayl was defined twice:
  !                       - in SR initr (4 values only, read from file initr.dat: likely from PIFM1)
  !                       - and here with a "data" instruction.
  !                       But the data instruction was not correct Fortran:
  !                       "[223 W] initialization of named COMMON should be in BLOCKDATA" (Forcheck W)
  !                       thus the first four values already initialised were not overwritten!
  !
  !                       Plus cosmetic: header, comments and cleaning (namely /cb07/, which
  !                       was defined in this routine, but used nowhere else)
  !
  ! Oct-2016  Josue Bock  Major work, cleaning, debuging, ...
  !                       Checked consistency of array size throughout the code
  !                       Reindexed some arrays (totrad, ...) to improve efficiency
  !                       Cosmetic: rearranged the routines order to facilitate reading
  !
  ! Mar-2017  Josue Bock  Reindexed all arrays from top to bottom. This version was not in the
  !                       first GitHub release, and has been done again in a more readable way.
  !
  ! Oct-2017  Josue Bock  Converted to f90, further cleaning/improvements in the code
  !
  ! Nov-2017  Josue Bock  Removed ntypa and ntypd from /cb02/, unused
  !                       Replaced hardcoded value of total solar energy (1340.0 W/m2) by
  !                         s0tot (which is slightly different: 1339.94507 W/m2)

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       mb,                  &
       mbs,                 &
       mbir,                &
       nrlay,               &
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Local parameterss:
  real (kind=dp), parameter :: u0min = 1.0e-2_dp ! minimum solar angle to compute solar bands

! Local scalars:
  integer :: ib, ig, ibanf
  integer :: jl      ! Loop index for Legendre coefficients
  integer :: jz, jzp ! Loop indexes in the vertical dimension (top-down)
  real (kind=dp) :: hk
  real (kind=dp) :: s0
  real (kind=dp) :: zbsca ! scattering coefficient for aerosols
  real (kind=dp) :: zdopr
  real (kind=dp) :: zfni, zfnip
  real (kind=dp) :: zfuq1, zfuq2
  real (kind=dp) :: zx0

! Local arrays:
  integer :: kg(mb)                   ! number of cumulative probabilities
  data kg / 10, 8, 12, 7, 12, 5, 2, 3, 4, 4, 3, 5, 2, 10, 12, 7, 7, 8 /
  real (kind=dp) :: dlam(7,nrlev,mb)   ! array with solar/ir rad. fluxes
  real (kind=dp) :: zgaer(nrlay)       ! needed for Legendre-coefficients
  real (kind=dp) :: sss(nrlev)

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

  common /cb10/ totrad (mb,nrlay)
  real (kind=dp) :: totrad

  common /cb15/ fnseb,flgeg,hr(nrlay)
  real (kind=dp) :: fnseb, flgeg, hr

  common /cb16/ u0,albedo(mbs),thk(nrlay)
  real (kind=dp) :: u0, albedo, thk

  common /cb19/ berayl(6),bea(mb,nrlay),baa(mb,nrlay),ga(mb,nrlay)
  real (kind=dp) :: berayl, bea, baa, ga

  common /extra/ waer(nrlay),taer(nrlay),plaer(2,nrlay)
  real (kind=dp) :: waer, taer, plaer

  common /extra_2/ taer_s(nrlay),taer_a(nrlay),ga_pl(nrlay) ! for photolysis calculation
  real (kind=dp) :: taer_s, taer_a, ga_pl

  common /kurz/ fs1(nrlev),fs2(nrlev),totds(nrlev),ss(nrlev), &
                fsn(nrlev),dtdts(nrlay)
  real (kind=dp) :: fs1, fs2, totds, ss, fsn, dtdts          ! integrated radiation fluxes

  common /lang/ fl1(nrlev),fl2(nrlev),fln(nrlev),dtdtl(nrlay)           ! integrated radiation fluxes
  real (kind=dp) :: fl1, fl2, fln, dtdtl

  common /leck2/ sf(nrlev),sw(nrlev),ssf(nrlev),ssw(nrlev), &             ! radiation fluxes
                 f2f(nrlev),f2w(nrlev),f1f(nrlev),f1w(nrlev)
  real (kind=dp) :: sf, sw, ssf, ssw, f2f, f2w, f1f, f1w

  common /planci/ pib(nrlev),pibs                              ! black-body radiation
  real (kind=dp) :: pib, pibs

  common /ray/ dtaur(nrlay), plr(2,nrlay)                         ! rayleigh
  real (kind=dp) :: dtaur, plr

  common /sol/ s0b(mbs), s0tot
  real (kind=dp) :: s0b, s0tot

  common /tmp2/ as(mbs),ee(mbir)                            ! albedo, and emissivity of the surface
  real (kind=dp) :: as, ee

  common /umcon/ umco2,umch4,umn2o                          ! gas concentrations
  real (kind=dp) :: umco2,umch4,umn2o

! == End of declarations =======================================================

  s0 = 1355.3_dp ! jjb this is the sum of the 4 former wavelength bands (1018.3+230.1+39.2+67.7)
                 !      see paper from Loughlin et al, (1997) QJRMS vol. 123, pp. 1985-2007, table 1
                 !     SHOULD THIS BE UPDATED? see also zfuq1 below: 1340.0
! s0 = 1360.3_dp
! s0 = 1368.0_dp ! latest AB version

! concentrations of trace gases ! jjb improve by reading in namelist
  umco2 = 330._dp
  umch4 = 1.6_dp
  umn2o = 0.28_dp

!
! ----------------------------------------------------------------------
! cloud fractions: frr includes the cloud fraction section taken from
! the new PIFM code.
!
  call frr

! Initialisation or arrays
  fs1(:) = 0._dp
  fs2(:) = 0._dp
  ss(:)  = 0._dp
  sss(:) = 0._dp
  totds(:) = 0._dp
  fsn(:) = 0._dp
  fl1(:) = 0._dp
  fl2(:) = 0._dp
  fln(:) = 0._dp
  totrad(:,:) = 0._dp
  dlam(:,:,:) = 0._dp

! at night: no solar radiation (ibanf=7)
  if ( u0 <= u0min ) then
     ibanf = mbs + 1
  else
     ibanf = 1
  endif

  zdopr = 2._dp * rho(nrlev)

! loop for 18 spectral bands
!---------------------------
  do ib = ibanf,mb

! rayleigh scattering
     if (ib <= mbs) then
        do jz=1,nrlay
           dtaur(jz) = berayl(ib) * thk(jz) * (rho(jz)+rho(jz+1))/zdopr
           plr(1,jz) = 0.0_dp
           plr(2,jz) = 0.5_dp
        enddo
     else
        dtaur(:) = 0.0_dp
        plr(:,:) = 0.0_dp
     endif

! ......................................................................
! aerosols and clouds particles with explicit optical depth
! (taer = optical depth)
!
!     Output of this section:
!     taer, waer and plaer, for use in SR tau.
!
!     In the solar case (band 1 only), also calculate variables for the
!     photolysis model: taer_s, taer_a and ga_pl

     if (ib <= mbs) then
        do jz = 1, nrlay
           ! -- taer
           taer(jz) = bea(ib,jz) * thk(jz)
           zbsca = bea(ib,jz) - baa(ib,jz)
           if (ib == 1) then
              taer_s(jz) = zbsca * thk(jz)            ! photolysis
              taer_a(jz) = baa(ib,jz) * thk(jz)       ! photolysis
           endif

           ! -- plaer
           if (zbsca+(dtaur(jz)/thk(jz)) >= 1.0e-20_dp) then
              zgaer(jz) = ga(ib,jz) * zbsca / (zbsca + (dtaur(jz) / thk(jz)) )
              do jl = 1, 2
                 plaer(jl,jz) = real(2*jl+1,dp) * zgaer(jz)**jl
              enddo
              if(ib == 1) ga_pl(jz) = ga(ib,jz)       ! photolysis
           else
              plaer(:,jz) = 0.0_dp
              if(ib == 1) ga_pl(jz) = 0.0_dp          ! photolysis
           endif

           ! -- waer
           if (bea(ib,jz) > 1.e-20_dp) then
              waer(jz) = 1.0_dp - (baa(ib,jz)/bea(ib,jz))
           else
              waer(jz) = 0.0_dp
           endif
        enddo
     else
        do jz=1,nrlay
           ! -- taer
           taer(jz) = bea(ib,jz) * thk(jz)
           ! -- plaer
           do jl=1,2
              plaer(jl,jz)=real(2*jl+1,dp)*ga(ib,jz)**jl
           enddo
           ! -- waer
           if (bea(ib,jz) > 1.e-20_dp) then
              waer(jz)=1.0_dp - (baa(ib,jz)/bea(ib,jz))
           else
              waer(jz)=0.0_dp
           endif
        enddo
     endif
! ......................................................................



     ! water droplets
     call water(ib)

     ! continuum absorption of water vapour
     call gascon(ib)

     ! black body radiation
         ! planck: explicit calculation
         ! plancktab: calculation from tabulated values by interpolation
     call planck(ib)
     ! call plancktab(ib)

! gas absorption, loop for cumulative probabilities
!--------------------------------------------------
     do ig=1,kg(ib)
        call gase(ib,ig,hk)

        ! total optical depth:
        call tau

        if ( ib <= mbs ) then ! solar bands

! solution of short-wave radiative transfer equation
!---------------------------------------------------

! radiation fluxes:
! ss    is direct solar downward
! sss   is direct solar downward + Delta-Eddington peak
! fs1   is total diffusive solar upward
! fs2   is diffusive solar downward

           call kurzw(ib,u0)
           call jeanfr(ib)

           ss(:)  =  ss(:) + ( sf(:) +  sw(:))*hk
           sss(:) = sss(:) + (ssf(:) + ssw(:))*hk
           fs1(:) = fs1(:) + (f1f(:) + f1w(:))*hk
           fs2(:) = fs2(:) + (f2f(:) + f2w(:))*hk
           dlam(1,:,ib) = dlam(1,:,ib) + ( sf(:) +  sw(:))*hk
           dlam(2,:,ib) = dlam(2,:,ib) + (ssf(:) + ssw(:))*hk
           dlam(3,:,ib) = dlam(3,:,ib) + (f1f(:) + f1w(:))*hk
           dlam(4,:,ib) = dlam(4,:,ib) + (f2f(:) + f2w(:))*hk

        else ! IR bands

! solution of long-wave radiative transfer equation
!--------------------------------------------------

! radiation fluxes:
! fl1   is longwave upward
! fl2   is longwave downward
!
! dlam contains the different contributions from solar diffuse and ir
! in each of the first index (see above). The last index is the
! wavelength bin identifier and the middle index is the height.
! Note that dlam is only filled over (1:4,:,1:mbs) and (5:7,:,mbs+1:mb)

           call langw(ib)
           call jeanfr(ib)

           fl1(:) = fl1(:) + (pib(:) - f1f(:) - f1w(:))*hk
           fl2(:) = fl2(:) + (pib(:) - f2f(:) - f2w(:))*hk
           dlam(5,:,ib) = dlam(5,:,ib) + (pib(:) - f1f(:) - f1w(:))*hk
           dlam(6,:,ib) = dlam(6,:,ib) + (pib(:) - f2f(:) - f2w(:))*hk
           dlam(7,:,ib) = dlam(7,:,ib) + pib(:)*hk

        endif
     enddo ! ig
  enddo ! ib

!---------------------------------------------------------------------------

! flux correction for given solar constant s0 and emissivity ee

  zfuq1 = s0 / s0tot
  zfuq2 = pibs * 0.03_dp * ee(mbir)

  ! Solar wavelength bands corrections
  if ( u0 > u0min ) then
     ss(:)  =  ss(:)*zfuq1
     sss(:) = sss(:)*zfuq1
     fs1(:) = fs1(:)*zfuq1
     fs2(:) = fs2(:)*zfuq1
     dlam(:4,:,:mbs) = dlam(:4,:,:mbs)*zfuq1
  end if

! dummy variables for additional output:
! totds is total solar downward
! fsn   is total solar net radiation flux
! fln   is total ir net radiation flux
! f1 and f2 are total upward / downward radiation fluxes
! fn    is total net radiation flux

  ! Delta-Eddington peak (sss-ss) is added to fs2
  if ( u0 > u0min ) then
     totds(:) = sss(:) + fs2(:)
     fs2(:) = totds(:) - ss(:)
     fsn(:) = fs1(:) - totds(:)
  end if
  fl1(:) = fl1(:) + zfuq2
  dlam(5,:,18) = dlam(5,:,18) + zfuq2 ! jjb mysterious ! fl2?
  fln(:) = fl1(:) - fl2(:)

  ! surface radiation fluxes
  flgeg = fl2(nrlev)
  fnseb = fs2(nrlev) + ss(nrlev) - fs1(nrlev)

! calculation of heating rates {hr},
! dtdts and dtdtl are the solar / ir heating rates (currently not used)
  do jz = 1, nrlay
     jzp = jz+1
     if (jz == 1) zfni = fl1(1) - fl2(1) + fs1(1) - ss(1) - fs2(1)
     zfnip     = fl1(jzp) - fl2(jzp) + fs1(jzp) - ss(jzp) - fs2(jzp)
     zx0       = (thk(jz) * (rho(jz) + rho(jzp)) * 502.5_dp)
     dtdts(jz) = (fs1(jzp)-ss(jzp)-fs2(jzp) - fs1(jz)+ss(jz)+fs2(jz)) / zx0
     dtdtl(jz) = (fl1(jzp)-fl2(jzp) - fl1(jz)+fl2(jz)) / zx0
     hr(jz)    = (zfnip-zfni) / zx0
     zfni      = zfnip
  enddo

! ------------------------------------------------------------------------------
  do jz=1,nrlay
     jzp = jz+1
     do ib=1,mbs
        totrad(ib,jz)=((dlam(2,jz,ib)+dlam(2,jzp,ib))/(2.0_dp*u0))+ &
                      (dlam(3,jz,ib)+dlam(3,jzp,ib)+ &
                       dlam(4,jz,ib)+dlam(4,jzp,ib))
     enddo
     do ib=mbs+1,mb
!       if(jz == 1) then                                        ! jjb problem here, same as below
        totrad(ib,jz)=-(dlam(7,jz,ib)+dlam(7,jzp,ib))*2._dp  & !     Should there be a different
                      +dlam(6,jz,ib)+dlam(6,jzp,ib) &          !     calculation for the first
                      +dlam(5,jz,ib)+dlam(5,jzp,ib)            !     (or last = ground) level?
!       else
!       totrad(ib,jz)=-(dlam(7,jz,ib)+dlam(7,jz+1,ib))*2._dp
!                     +dlam(6,jz,ib)+dlam(6,jz+1,ib)
!                     +dlam(5,jz,ib)+dlam(5,jz+1,ib)
!       endif
     enddo
  enddo

end subroutine nstrahl

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
! the following routines are cut and paste from the new PIFM code with the
! exception of several common blocks which have been altered!!
! the routines tau and water are the only two exceptions. Water is only used
! if load0 is called otherwise it returns a zero optical depth. The new version
! of PIFM includes a routine to combine all the optical depths etc...that has been
! re-written for this program and is called tau.
!

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine frr
!
! Description :
! -----------
!    Calculation of random overlapping cloud coverage.
!    --------------------------------------------------------------------------
!    bb: cloud free part, cc: cloudy part

! References :
! ----------
!    Geleyn and Hollingsworth, Betrei. Phys. Atmosph. Vol. 52 No. 1, February 1979
!        see page 7, formulas 16. Hereafter GH.
!
!    also cited by Zdunkowski et al., Betrei. Phys. Atmosph. Vol. 55 No. 3, August 1982
!        pp. 215-238 (see page 219: coefficients b_1,i to b_4,i) Note that these authors
!        oriented the vertical indexes from bottom to top


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY", reference
  !                       Declarations and implicit none
  !
  ! Oct-2016  Josue Bock  BUGFIX 1: frac was indexed in the wrong line (ground -> top)
  !                       BUGFIX 2a: the general definitions are valid and have to be applied
  !                                  to b1 and b3 for the bottom layer,
  !                                  and to b2 and b4 for the top layer.
  !                       BUGFIX 2b: the initialisation bb(1,1)=1.-frac(1) was thus wrong (but
  !                                  probably without effect since frac(1) was 0.)
  !                       Rewritten with if...else... structure. The way it is done avoids
  !                       to have equalities in if tests.
  !
  ! Oct-2017  Josue Bock  Converted to f90, further cleaning/improvements in the code
  !
  ! May-2021  Josue Bock  Double checked, commented, the <= and >= tests have been changed
  !                       so that equality cases lead to the simplest calculation

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlay

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  integer :: jz, jzm, jzp

! Common blocks:
  common /cb09/ frac(nrlay),rew(nrlay),rho2w(nrlay)
  real (kind=dp) :: frac, rew, rho2w

  common /part/ cc(4,nrlay),bb(4,nrlay)
  real (kind=dp) :: cc, bb

! == End of declarations =======================================================

  do jz = 1,nrlay
     jzm = jz-1
     jzp = jz+1

     ! GH equation 16a: b(1,j) = (1 - max(C(j),C(j-1))) / (1 - C(j-1))
     ! GH equation 16c: b(3,j) = min(C(j),C(j-1)) / C(j-1)
     !   where C = frac
     !   all indetermination cases = 1
     if (jz==1)then ! uppermost model layer
        bb(1,1) = 1._dp
        bb(3,1) = 1._dp
     else
        ! dependance of cloud coverage in the layer above
        if (frac(jzm) > 0._dp) then
           if (frac(jzm) < 1._dp) then
              if (frac(jzm) < frac(jz)) then
                 bb(1,jz) = (1._dp-frac(jz)) / (1._dp-frac(jzm))
                 bb(3,jz) = 1._dp
              else ! 0. < frac(jz) <= frac(jzm) < 1.
                 bb(1,jz) = 1._dp
                 bb(3,jz) = frac(jz) / frac(jzm)
              end if
           else    ! frac(jzm) = 1
              bb(1,jz) = 1._dp
              bb(3,jz) = frac(jz)
           end if
        else       ! frac(jzm) = 0
           bb(1,jz) = 1._dp - frac(jz)
           bb(3,jz) = 1._dp
        end if
     end if

     ! GH equation 16b: b(2,j) = (1 - max(C(j),C(j+1))) / (1 - C(j+1))
     ! GH equation 16d: b(4,j) = min(C(j),C(j+1)) / C(j+1)
     !   where C = frac
     !   all indetermination cases = 1
     if (jz==nrlay) then ! surface
        bb(2,nrlay) = 1._dp
        bb(4,nrlay) = 1._dp
     else
        ! dependance of cloud coverage in the underneath layer
        if (frac(jzp) > 0._dp) then
           if (frac(jzp) < 1._dp) then
              if (frac(jz) > frac(jzp)) then
                 bb(2,jz) = (1._dp-frac(jz)) / (1._dp-frac(jzp))
                 bb(4,jz) = 1._dp
              else ! 0. < frac(jz) <= frac(jzp) < 1.
                 bb(2,jz) = 1._dp
                 bb(4,jz) = frac(jz) / frac(jzp)
              end if
           else    ! frac(jzp) = 1
              bb(2,jz) = 1._dp
              bb(4,jz) = frac(jz)
           end if
        else       ! frac(jzp) = 0
           bb(2,jz) = 1._dp-frac(jz)
           bb(4,jz) = 1._dp
        end if
     end if

  end do

  ! cloudy part: cc = 1 - bb
  cc(:,:) = 1._dp - bb(:,:)

end subroutine frr

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine water(ib)
!
! Description :
! -----------
!   Calculation of optical depth {t2w}, single scattering albedo {w2w} and
!   Legendre coefficients {pl2w} of the phase-function for Mie-scattering
!   of cloud droplets. (eqs. 4.25 - 4.27 diss. Fu (1991))
!-------------------------------------------------------------------------
!   rho2w and reff=rew are the partial density of cloud water and the effective
!   radius of the cloud droplets. Effective radii lower or higher than the
!   boundary values of 4.18 um and 31.18 um are set to the corresponding
!   boundary values. Clouds with rho2w less than 10**-5 kg m**-3 are not
!   calculated.
!-------------------------------------------------------------------------
!

! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  !                       pi no longer hardcoded (removed in a later version)
  !
  ! Oct-2016  Josue Bock  Correction of size of array r(nkt) instead of r(150)
  !                       Missing initialisation for rew -- probably no consequence
  !                         since rew was defined from nfog0 to nfog1 (lcl to lct),
  !                         and apart from this range rho2w < threshold=1e-5
  !                       BUGFIX: rho2w(i) was used instead of rhow in r(jt) calculation
  !
  ! Oct-2017  Josue Bock  Reactivated this routine by removing the inactivation switch
  !                       icld, and correction of several bugs:
  !                       BUGFIX 1: r(jt) was 0.01 times too small (error in units conversion)
  !                       BUGFIX 2: rew(jz) must be shifted by -1 as compared to ff, to account
  !                                 that layer 2 (not 1) in Mistra corresponds to layer nrlay
  !                                 in the radiative code.
  !
  !                       Moved the definition of rew(jz) in SR load1 (with the calculation of
  !                         r(nkt) in SR grid, since this doesn't change with time)
  !
  !                       Conversion of the code in Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlay,               &
       ncw,                 &
       mb,                  &
       mbs

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  integer, intent(in) :: ib    ! current wavelength band

! Local scalars:
  integer :: jl         ! Loop index for Legendre coefficients
  integer :: jz         ! Loop index in the vertical dimension (top-down)
  integer :: k          ! Lower index of tabulated ret values such that ret(k) < rew=reff <= ret(k+1)
  real (kind=dp) :: gg  ! g2wt interpolation

! Common blocks:
  common /cb09/ frac(nrlay),rew(nrlay),rho2w(nrlay)
  real (kind=dp) :: frac, rew, rho2w

  common /cb16/ u0,albedo(mbs),thk(nrlay)
  real (kind=dp) :: u0, albedo, thk

  common /h2o/ t2w(nrlay), w2w(nrlay), pl2w(2,nrlay)
  real (kind=dp) :: t2w, w2w, pl2w

  common /was1/ ret(ncw),r2wt(ncw),b2wt(ncw,mb),w2wt(ncw,mb), g2wt(ncw,mb) ! 't' = tabulated values
  real (kind=dp) :: ret, r2wt, b2wt, w2wt, g2wt

! == End of declarations =======================================================

  do jz=1,nrlay

! no clouds
     if ( rho2w(jz) < 1.0e-5_dp ) then
        t2w(jz)=0.0_dp
        w2w(jz)=0.0_dp
        pl2w(:,jz)=0.0_dp
     else

! interpolation for given reff and rho2w from tabulated values
! {ret} and {r2wt}

        ! lower limit 4.18 um
        if ( rew(jz) <= ret(1) ) then
           t2w(jz) = thk(jz) * rho2w(jz) * b2wt(1,ib) / r2wt(1)
           w2w(jz) = w2wt(1,ib)
           do jl=1,2
              pl2w(jl,jz) = real(2*jl+1,dp) * g2wt(1,ib)**jl
           end do

        ! upper limit 31.18 um
        elseif ( rew(jz) >= ret(ncw) ) then
           t2w(jz) = thk(jz) * rho2w(jz) * b2wt(ncw,ib) / r2wt(ncw)
           w2w(jz) = w2wt(ncw,ib)
           do jl=1,2
              pl2w(jl,jz) = real(2*jl+1,dp) * g2wt(ncw,ib)**jl
           end do

        ! linear interpolation for values between the limits
        else
           k = 1
           do while (rew(jz) > ret(k+1))
              k = k + 1
           end do
           t2w(jz) = thk(jz) * rho2w(jz) * ( b2wt(k,ib)/r2wt(k) + &
                    ( b2wt(k+1,ib)/r2wt(k+1) - b2wt(k,ib)/r2wt(k) ) / &
                    ( 1.0_dp/ret(k+1) - 1.0_dp/ret(k) ) * &
                    ( 1.0_dp/rew(jz)  - 1.0_dp/ret(k) ) )
           w2w(jz) = w2wt(k,ib) + ( w2wt(k+1,ib) - w2wt(k,ib) ) / &
                    ( ret(k+1) - ret(k) ) * ( rew(jz) - ret(k) )
           gg = g2wt(k,ib) + ( g2wt(k+1,ib) - g2wt(k,ib) ) / &
                    ( ret(k+1) - ret(k) ) * ( rew(jz) - ret(k) )
           do jl=1,2
              pl2w(jl,jz) = real(2*jl+1,dp) * gg**jl
           end do
        endif
     endif
  end do

end subroutine water

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine gascon ( ib )
!
! Description :
! -----------
!    Continuum absorption of water vapour, bands 11 - 17

! References :
! ----------
!    For partial reference, see Loughlin et al. 1997, QJRMS vol 123 (543),
!       pp. 1985-2007 doi: 10.1002/qj.49712354311


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  ! Oct-2016  Josue Bock  implicit none
  ! Oct-2017  Josue Bock  Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlay

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  integer, intent(in) :: ib

! Local arrays:
  real (kind=dp), parameter :: vv(7) = &    ! central wavenumbers of the spectral bands in in cm**-1
       (/ 1175.0_dp, 1040.0_dp, 890.0_dp, 735.0_dp, &
           605.0_dp,  470.0_dp, 340.0_dp /)

! Common blocks:
  common /con/ tgcon(nrlay)
  real (kind=dp) :: tgcon

! == End of declarations =======================================================

  if ( ib >= 11 .and. ib <= 17 ) then
     call qopcon ( vv(ib-10) )
  else
     tgcon(:)=0.0_dp
  endif

end subroutine gascon

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine qopcon ( vv )
!
! Description :
! -----------
!    e-Type- and continuum absorption after eq (A.19) of
!    dissertation Fu. (1991)
!
!    vv are the central wavenumbers of the spectral bands in in cm**-1.
!    All other variables in MKS:
!    p1  =  partial pressure of water vapour in Pa
!    ff  =  continuum/ e-Type- absorption coefficients in m**2 kg**-1
!    delta-tau  =  ff * xm1 *rho * delta-z  =  - ff * xm1 * delta-p / gE,
!    ff * xm1 arithmetical mean values
!


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  ! Oct-2016  Josue Bock  implicit none, removed arithmetic do loops
  ! Oct-2017  Josue Bock  Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlay,               &
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: vv

! Local scalars:
  integer :: jz
  real (kind=dp), parameter :: &
       r = 0.00002_dp,         &
       x = 418._dp,            &
       y = 557780._dp,         &
       zz = 0.00787_dp
  real (kind=dp) :: s,w

! Local arrays:
  real (kind=dp) :: ff(nrlev), p1(nrlev)

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

  common /con/ tgcon(nrlay)
  real (kind=dp) :: tgcon

! == End of declarations =======================================================

  s = ( x + y * exp ( - zz * vv ) ) / 101325._dp

  do jz=1,nrlev
     p1(jz) = p(jz) * xm1(jz) / ( 0.622_dp + 0.378_dp * xm1(jz) )
     w = exp ( 1800.0_dp / t(jz) - 6.08108_dp )
     ff(jz) = s * ( p1(jz)/100._dp + r * p(jz) ) * w
  enddo

  do jz=1,nrlay
     tgcon(jz) = ( ff(jz) * xm1(jz) + ff(jz+1) * xm1(jz+1) )* &
                 ( p(jz+1) - p(jz) ) * 0.00509892_dp
  enddo

end subroutine qopcon

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine planck ( ib )
!
! Description :
! -----------
!   Calculation of over {ib} integrated Planck-function for every model layer,
!   and of the radiation temperature of the surface.


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  !                       implicit none
  !                       pi imported from constants module
  !
  ! Oct-2017  Josue Bock  Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE constants, ONLY : &
! Imported Parameters:
       pi

  USE global_params, ONLY : &
! Imported Parameters:
       nrlev,               &
       mbir,                &
       mbs

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  integer, intent(in) :: ib

! Local scalars:
  integer :: jz

! Local arrays:
  real (kind=dp), parameter :: wvl(mbir+1) = &
       (/ 2200._dp, 1900._dp, 1700._dp, 1400._dp, 1250._dp, 1100._dp, &
           980._dp,  800._dp,  670._dp,  540._dp,  400._dp,  280._dp, 0._dp /)

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

  common /planci/ pib(nrlev),pibs
  real (kind=dp) :: pib, pibs

! External function:
  real (kind=dp), external :: plkavg

! == End of declarations =======================================================

  if(ib > mbs) then
     do jz=1,nrlev
        pib(jz)=pi*plkavg(wvl(ib-5),wvl(ib-6),t(jz))
     enddo
     pibs=pi*plkavg(wvl(ib-5),wvl(ib-6),ts)
  end if

end subroutine planck

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
function plkavg ( WNUMLO, WNUMHI, xT )
!
! Description :
! -----------
!   Calculation of black body radiation for given temperature {xt}
!   For low wavenumbers the calculation is done with power series, for high
!   wavenumbers with exponential series. Critical wavenumber is 'vcut'.


! Modifications :
! -------------
  ! 16-Jul-2015  Josue Bock  First Forcheck correction, real*4 => real*8 (inconsistency)
  !    Jul-2016  Josue Bock  Header, comments and cleaning
  !
  ! 15-Mar-2017  Josue Bock  removed unused part of code (situation that never happen,
  !                            due to the calling structure, with tabulated values for
  !                            wnumlo and wnumhi)
  !                          use module constants for pi (instead of defining it as asin(1.0) )
  !                          missing declarations and implicit none
  !                          explicit conversion I -> R using real( )
  !                          updated labeled do loops
  !
  ! 28-Oct-2017  Josue Bock  Fortran90
  !                          replaced final test == 0. by < tiny_dp

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! External subroutine
       abortM

  USE constants, ONLY : &
! Imported Parameters:
       pi

  USE file_unit, ONLY :     &
! Imported Parameters:
       jpfunerr

  USE precision, ONLY :     &
! Imported Parameters:
       dp,                  &
       tiny_dp

  implicit none

! Function declaration:
  real (kind=dp) :: plkavg

! Function arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: &
       wnumlo, wnumhi,          &  ! wavenumbers
       xt                          ! temperature

! Local parameters:
  real (kind=dp), parameter :: &   ! coefficients for power series
       a1 = 1._dp/3._dp,       &
       a2 =-1._dp/8._dp,       &
       a3 = 1._dp/60._dp,      &
       a4 =-1._dp/5040._dp,    &
       a5 = 1._dp/272160._dp,  &
       a6 =-1._dp/13305600._dp

  real (kind=dp), parameter :: &
       c2 = 1.438786_dp,       &   ! Planck radiation constant (K*m*10^2)
       conc = 15._dp / pi**4,  &
       sigma = 5.67032E-8_dp,  &
       sigdpi = sigma / pi,    &
       vcut = 1.5_dp               ! critical wave number
  real (kind=dp), parameter :: &
       vcp(7) = (/10.25_dp, 5.7_dp, 3.9_dp, 2.9_dp, 2.3_dp, 1.9_dp, 0.0_dp/) ! wavenumbers for exp series

! Local scalars:
  integer :: ismallv              ! number of "small" cases (v<vcut ==> power series)
  integer :: jj, jm               ! loop indexes
  integer :: mmax
  real (kind=dp) :: ex, exm, mv   ! exponential series
  real (kind=dp) :: vsq           ! power series, v squared

! Local arrays:
  real (kind=dp) :: d(2), p(2), v(2)

! == End of declarations =======================================================

! Check input temperature
  if( xt < 0.0_dp ) then
     write(jpfunerr,*)'Error in SR PLKAVG -- negative temperature'
     call abortM ('Stopped by SR plkavg (radiative code)')
  end if

  if ( xt < 1.e-4_dp ) then
     plkavg = 0.0_dp
     return
  endif

! planck
  v(1) = c2 * wnumlo / xt
  v(2) = c2 * wnumhi / xt


  ismallv = 0
  ! lower (1) and upper (2) wavenumber
  do  jj = 1, 2

     if ( v(jj) < vcut ) then
! ** use power series
        ismallv = ismallv + 1
        vsq = v(jj)**2
        p(jj) = conc * vsq * v(jj) * ( a1 + v(jj) * ( a2 + v(jj) * &
                ( a3 + vsq * ( a4 + vsq * ( a5 + vsq * a6 ) ) ) ) )

     else
! ** use exponential series
        mmax = 1
        ! ** find upper limit of series
        do while ( v(jj) < vcp(mmax) )
           mmax = mmax + 1
        end do
        ex = exp ( - v(jj) )
        exm = 1.0_dp
        d(jj) = 0.0_dp

        do jm = 1, mmax
           mv = real(jm,dp) * v(jj)
           exm = ex * exm
           d(jj) = d(jj) + exm * ( 6._dp + mv * ( 6._dp + mv * ( 3._dp + mv ) ) ) &
                / real(jm**4,dp)
        end do
        d(jj) = conc * d(jj)
     end if
  end do

! black-body radiation in current spectral band:
! difference between integral of upper and lower boundary
  if ( ismallv == 2 ) then
     ! ** wnumlo and wnumhi both small
     plkavg = p(2) - p(1)
  else if ( ismallv == 1 ) then
     ! ** wnumlo small, wnumhi large
     plkavg = 1._dp - p(1) - d(2)
  else
     ! ** wnumlo and wnumhi both large
     plkavg = d(1) - d(2)
  end if
  plkavg = sigdpi * xt**4 * plkavg

! Final test: warn if zero
  if ( plkavg < tiny_dp ) then
     print*,'plkavg--returns zero; possible underflow'
  end if

end function plkavg

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
! jjb 24/07/2016 unused, thus commented
!
!!$subroutine plancktab ( ib )
!!$!
!!$! Description :
!!$! -----------
!!$!   Calculation of Planck-function integrated over {ib},
!!$!   from tabled values, for every model layer.
!!$
!!$
!!$! Modifications :
!!$! -------------
!!$  !    Jul-2016  Josue Bock  Header including "USE ... ONLY"
!!$  ! 28-Oct-2017  Josue Bock  Fortran90, cleaned and formatted the same as planck
!!$
!!$! == End of header =============================================================
!!$
!!$! Declarations :
!!$! ------------
!!$! Modules used:
!!$
!!$  USE global_params, ONLY : &
!!$! Imported Parameters:
!!$       nrlay,               &
!!$       nrlev,               &
!!$       mbs
!!$
!!$  USE precision, ONLY :     &
!!$! Imported Parameters:
!!$       dp
!!$
!!$  implicit none
!!$
!!$! Subroutine arguments
!!$! Scalar arguments with intent(in):
!!$  integer, intent(in) :: ib
!!$
!!$! Local scalars:
!!$  integer :: ibir ! reduced index for IR only (1 - 12)
!!$  integer :: jz   ! loop index, top to bottom
!!$
!!$! Common blocks:
!!$  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
!!$  real (kind=dp) :: t, p, rho, xm1, ts
!!$
!!$  common /planci/ pib(nrlev),pibs
!!$  real (kind=dp) :: pib, pibs
!!$
!!$! External function:
!!$  real (kind=dp), external :: fst4
!!$
!!$! == End of declarations =======================================================
!!$
!!$  if(ib > mbs) then
!!$     ibir=ib-mbs
!!$     do jz=1,nrlev
!!$        pib(jz) = fst4(ibir,t(jz))
!!$     end do
!!$     pibs = fst4(ibir,ts)
!!$  end if
!!$
!!$end subroutine plancktab

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
! jjb 25/07/2016: used only by SR plancktab, whose call is commented
!                  -> not used anymore, thus commented out
!
!!$function fst4 (ibir,t)
!!$!
!!$! Description :
!!$! -----------
!!$!   Interpolation of black body radiation from tabulated values {pibtab}
!!$!   for 35 temperatures {ttab}.
!!$!   Data in common block plancd.
!!$
!!$
!!$! Modifications :
!!$! -------------
!!$  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
!!$  ! Oct-2017  Josue Bock  Fortran90
!!$  !                       added test to detect temperature out of range
!!$  !                       interpolation rewritten, removed arithmetic if
!!$  !                       reindexed pibtab(mbir,35) => pibtab(35,mbir)
!!$
!!$! == End of header =============================================================
!!$
!!$! Declarations :
!!$! ------------
!!$! Modules used:
!!$
!!$  USE global_params, ONLY : &
!!$! Imported Parameters:
!!$       mbir
!!$
!!$  USE precision, ONLY :     &
!!$! Imported Parameters:
!!$       dp
!!$
!!$  implicit none
!!$
!!$! Function declaration:
!!$  real (kind=dp) :: fst4
!!$
!!$! Function arguments
!!$! Scalar arguments with intent(in):
!!$  integer, intent(in) :: ibir
!!$  real(kind=dp), intent(in) :: t
!!$
!!$! Local scalars:
!!$  integer :: itl, ith ! indexes in table, low, high, such that ttab(itl) <= t <= ttab(ith)
!!$
!!$! Common blocks:
!!$  common /plancd/ ttab(35),pibtab(35,mbir)
!!$  real(kind=dp) :: ttab, pibtab
!!$
!!$! == End of declarations =======================================================
!!$
!!$  if ( t<ttab(1) .or. t>ttab(35) ) then
!!$     write(0,*) 'Error in FN fst4: temperature out of range',t
!!$     stop 'Stopped by FN fst4 (radiative code)'
!!$  else
!!$     ith = 2
!!$     do while ( t > ttab(ith) )
!!$        ith = ith + 1
!!$     end do
!!$     itl = ith - 1
!!$     fst4 = pibtab(itl,ibir) + (pibtab(ith,ibir)-pibtab(itl,ibir)) * (t-ttab(itl)) * 0.2_dp
!!$  end if
!!$
!!$end function fst4

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine gase ( ib, ig, hk )
!
! Description :
! -----------
!    Calculation of gas absorption for different gases and 18 spectral bands.
!    The subroutines qks (solar) and qki (ir) calculate the absorption
!    coefficients 'fkg' from tabulated values.
!    tg(nrlay) are the optical depths due to nongray gaseous absorption, in
!    nrlay layers for a given band ib and cumulative probability ig. They are
!    calculated from fkg in the subroutines qop{'gas formula'}.


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header
  !                       Explicit array index (start:end) in several calls
  !                       Removal of labeled do-loops, cleaning
  !
  ! Oct-2017  Josue Bock  Fortran90, plus select case instead of "computed" goto
  !                       matrix operations instead of vertical loops
  !
  ! Nov-2017  Josue Bock  Replaced hardcoded values of solar energies by the values
  !                         read from data file, and carried in CB /sol/

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! External subroutine
       abortM

  USE file_unit, ONLY :     &
! Imported Parameters:
       jpfunerr

  USE global_params, ONLY : &
! Imported Parameters:
       mbs,                 &
       nrlay,               &
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  integer, intent(in) :: ib, ig

! Scalar arguments with intent(out):
  real (kind=dp), intent(out) :: hk

! Local scalars:
  integer :: jz ! loop indexes
  real (kind=dp) :: fk

! Local arrays:
  real (kind=dp) :: fkg(nrlev), fkga(nrlev), fkgb(nrlev) ! absorption coefficients
  real (kind=dp) :: pq(nrlev)                            ! specific humidity for bands 14 and 15
  real (kind=dp) :: tg1(nrlay), tg2(nrlay), tg3(nrlay)   ! optical depths

! Common blocks:
  common /band1/ hk1(10), fk1o3(10)
  real (kind=dp) :: hk1, fk1o3
  common /band2/ hk2(8), c2h2o(3,11,8)
  real (kind=dp) :: hk2, c2h2o
  common /band3/ hk3(12), c3h2o(3,11,12)
  real (kind=dp) :: hk3, c3h2o
  common /band4/ hk4(7), c4h2o(3,11,7)
  real (kind=dp) :: hk4, c4h2o
  common /band5/ hk5(12), c5h2o(3,11,12)
  real (kind=dp) :: hk5, c5h2o
  common /band6/ hk6(5), c6h2o(3,11,5)
  real (kind=dp) :: hk6, c6h2o
  common /band7/ hk7(2), c7h2o(3,19,2)
  real (kind=dp) :: hk7, c7h2o
  common /band8/ hk8(3), c8h2o(3,19,3)
  real (kind=dp) :: hk8, c8h2o
  common /band9/ hk9(4), c9h2o(3,19,4)
  real (kind=dp) :: hk9, c9h2o
  common /band10/ hk10(4),c10h2o(3,19,4),c10ch4(3,19),c10n2o(3,19)
  real (kind=dp) :: hk10, c10h2o, c10ch4, c10n2o
  common /band11/ hk11(3),c11h2o(3,19,3),c11ch4(3,19),c11n2o(3,19)
  real (kind=dp) :: hk11, c11h2o, c11ch4, c11n2o
  common /band12/ hk12(5), c12o3(3,19,5), c12h2o(3,19)
  real (kind=dp) :: hk12, c12o3, c12h2o
  common /band13/ hk13(2), c13h2o(3,19,2)
  real (kind=dp) :: hk13, c13h2o
  common /band14/ hk14(10), c14hca(3,19,10), c14hcb(3,19,10)
  real (kind=dp) :: hk14, c14hca, c14hcb
  common /band15/ hk15(12), c15hca(3,19,12), c15hcb(3,19,12)
  real (kind=dp) :: hk15, c15hca, c15hcb
  common /band16/ hk16(7), c16h2o(3,19,7)
  real (kind=dp) :: hk16, c16h2o
  common /band17/ hk17(7), c17h2o(3,19,7)
  real (kind=dp) :: hk17, c17h2o
  common /band18/ hk18(8), c18h2o(3,19,8)
  real (kind=dp) :: hk18, c18h2o

  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

  common /gas/ tg(nrlay)
  real (kind=dp) :: tg

  common /sol/ s0b(mbs), s0tot
  real (kind=dp) :: s0b, s0tot

  common /umcon/ umco2,umch4,umn2o
  real (kind=dp) :: umco2,umch4,umn2o

! == End of declarations =======================================================

  select case (ib)

  case (1)
! Band 1:
! In this band ( 50000 - 14500 cm**-1 ), we consider the nongray gaseous absorption of O3.
! solar energy contained in the band: s0b(1) = 619.618 W m**-2
     fk = fk1o3(ig)
     call qopo3s ( fk, tg )
     hk = s0b(1) * hk1(ig)

  case (2)
! Band 2:
! In this band ( 14500 - 7700 cm**-1 ), we consider the nongray gaseous absorption of H2O.
! solar energy contained in the band: s0b(2) = 484.295 W m**-2
     call qks ( c2h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg )
     hk = s0b(2) * hk2(ig)

  case (3)
! Band 3:
! In this band ( 7700 - 5250 cm**-1 ), we consider the nongray gaseous absorption of H2O.
! solar energy contained in the band: s0b(3) = 149.845 W m**-2
     call qks ( c3h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg )
     hk = s0b(3) * hk3(ig)

  case (4)
! Band 4:
! In this band ( 5250 - 4000 cm**-1 ), we consider the nongray gaseous absorption of H2O.
! solar energy contained in the band: s0b(4) = 48.7302 W m**-2
     call qks ( c4h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg )
     hk = s0b(4) * hk4(ig)

  case (5)
! Band 5:
! In this band ( 4000 - 2850 cm**-1 ), we consider the nongray gaseous absorption of H2O.
! solar energy contained in the band: s0b(5) = 31.6576 W m**-2
     call qks ( c5h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg )
     hk = s0b(5) * hk5(ig)

  case (6)
! Band 6:
! In this band ( 2850 - 2500 cm**-1 ), we consider the nongray gaseous absorption of H2O.
! solar energy contained in the band: s0b(6) = 5.79927 W m**-2
     call qks ( c6h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg )
     hk = s0b(6) * hk6(ig)

  case (7)
! Band 7:
! In this band ( 2200 - 1900 cm**-1 ), we consider the nongray gaseous absorption of H2O.
     call qki ( c7h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg )
     hk = hk7(ig)

  case (8)
! Band 8:
! In this band ( 1900 - 1700 cm**-1 ), we consider the nongray gaseous absorption of H2O.
     call qki ( c8h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg )
     hk = hk8(ig)

  case (9)
! Band 9:
! In this band ( 1700 - 1400 cm**-1 ), we consider the nongray gaseous absorption of H2O.
     call qki ( c9h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg )
     hk = hk9(ig)

  case (10)
! Band 10:
! In this band ( 1400 - 1250 cm**-1 ), we consider the overlapping
! absorption of H2O, CH4, and N2O by approach one of Fu(1991).
     call qki ( c10h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg1 )
     call qki ( c10ch4, fkg )
     call qopch4 ( fkg, tg2 )
     call qki ( c10n2o, fkg )
     call qopn2o ( fkg, tg3 )
     tg(:) = tg1(:) + tg2(:)/1.6_dp*umch4 + tg3(:)/0.28_dp*umn2o
     hk = hk10(ig)

  case (11)
! Band 11:
! In this band ( 1250 - 1100 cm**-1 ), we consider the overlapping
! absorption of H2O, CH4, and N2O by approach one of Fu(1991).
     call qki ( c11h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg1 )
     call qki ( c11ch4, fkg )
     call qopch4 ( fkg, tg2 )
     call qki ( c11n2o, fkg )
     call qopn2o ( fkg, tg3 )
     tg(:) = tg1(:) + tg2(:)/1.6_dp*umch4 + tg3(:)/0.28_dp*umn2o
     hk = hk11(ig)

  case (12)
! Band 12:
! In this band ( 1100 - 980 cm**-1 ), we consider the overlapping
! absorption of H2O and O3 by approach one of Fu(1991).
     call qkio3 ( c12o3(:,:,ig), fkg )
     call qopo3i ( fkg, tg1 )
     call qki ( c12h2o, fkg )
     call qoph2o ( fkg, tg2 )
     tg(:) = tg1(:) + tg2(:)
     hk = hk12(ig)

  case (13)
! Band 13:
! In this band ( 980 - 800 cm**-1 ), we consider the nongray gaseous absorption of H2O.
     call qki ( c13h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg )
     hk = hk13(ig)

  case (14)
! Band 14
! In this band ( 800 - 670 cm**-1), we consider the overlapping
! absorption of H2O and CO2 by approach two of Fu(1991).
     do jz=1, nrlev
        if ( p(jz) >= 6310._dp ) then
           pq(jz)=xm1(jz)
        else
           pq(jz)=0.0_dp
        endif
     end do
     call qki ( c14hca(:,:,ig), fkga )
     call qki ( c14hcb(:,:,ig), fkgb )
     fkg(:) = fkga(:)/330.0_dp*umco2 + pq(:) * fkgb(:)
     call qophc ( fkg, tg)
     hk = hk14(ig)

  case (15)
! Band 15:
! In this band ( 670 - 540 cm**-1), we consider the overlapping
! absorption of H2O and CO2 by approach two of Fu(1991).
     do jz=1, nrlev
        if ( p(jz) >= 6310._dp ) then
           pq(jz)=xm1(jz)
        else
           pq(jz)=0.0_dp
        endif
     end do
     call qki ( c15hca(:,:,ig), fkga )
     call qki ( c15hcb(:,:,ig), fkgb )
     fkg(:) = fkga(:)/330.0_dp*umco2 + pq(:) * fkgb(:)
     call qophc ( fkg, tg)
     hk = hk15(ig)

  case (16)
! Band 16:
! In this band ( 540 - 400 cm**-1 ), we consider the nongray
! gaseous absorption of H2O.
     call qki ( c16h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg )
     hk = hk16(ig)

  case (17)
! Band 17:
! In this band ( 400 - 280 cm**-1 ), we consider the nongray gaseous absorption of H2O.
     call qki ( c17h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg )
     hk = hk17(ig)

  case (18)
! Band 18:
! In this band ( 280 - 000 cm**-1 ), we consider the nongray gaseous absorption of H2O.
     call qki ( c18h2o(:,:,ig), fkg )
     call qoph2o ( fkg, tg )
     hk = hk18(ig)

  case default
     write(jpfunerr,*)'Error in SR gase: wrong spectral band index'
     call abortM ('Stopped by SR gase')

  end select

end subroutine gase

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine qks ( coefks, fkg )
!
! Description :
! -----------
!    Calculation of the absorption coefficients for solar spectral bands.
!
! fkg(nrlev) are the gaseous absorption coefficients in units of (cm-atm)**-1
! for a given cumulative probability in nrlev layers.
! coefks(3,11) are the coefficients to calculate the absorption coefficient
! at the temperature t for the 11 pressures by
!         ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
! and the absorption coefficient at conditions other than those eleven
! pressures is interpolated linearly with pressure (Fu, 1991).


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  !                       Declarations and implicit none
  !
  ! Oct-2017  Josue Bock  Fortran90
  !                       minor improvements of code

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Array arguments with intent(in):
  real (kind=dp), intent(in) :: coefks(3,11) ! tabulated coefficients
! Array arguments with intent(out):
  real (kind=dp), intent(out) :: fkg(nrlev)  ! gaseous absorption coefficients

! Local parameters:
  real (kind=dp), parameter :: ztref = 245.0_dp

! Local scalars:
  integer :: ipl, iph                ! index low, high such that stanp(ipl) < p <stanp(iph)
  integer :: jz                      ! vertical loop index (top-down)
  real (kind=dp) :: x1, x2           ! b and c of the above formula
  real (kind=dp) :: ztfact, ztfact2  ! temperature factor, squared

! Local arrays:
  real (kind=dp), parameter :: stanp(11) = &
       (/ 1000._dp,  1580._dp,  2510._dp,  3980._dp,  6310._dp, &
         10000._dp, 15800._dp, 25100._dp, 39800._dp, 63100._dp, &
        100000._dp /)

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

! == End of declarations =======================================================

  do jz=1, nrlev
     ztfact  = t(jz) - ztref
     ztfact2 = ztfact * ztfact

     if ( p(jz) <= stanp(1) ) then
        ! pressure lower than 10 hPa: interpolate between 0. and 10 hPa
        x1 = exp( coefks(1,1) + coefks(2,1)*ztfact + coefks(3,1)*ztfact2 )
        fkg(jz) = x1 * p(jz) / stanp(1)

     else if ( p(jz) >= stanp(11) ) then
        ! pressure higher than 1000 hPa: extrapolate the slope of the latest two
        ! tabulated values
        x1 = exp( coefks(1,10) + coefks(2,10)*ztfact + coefks(3,10)*ztfact2 )
        x2 = exp( coefks(1,11) + coefks(2,11)*ztfact + coefks(3,11)*ztfact2 )
        fkg(jz) = x1 + ( x2 - x1 ) / ( stanp(11) - stanp(10) ) &
                 * ( p(jz) - stanp(10) )

     else
        ! other cases: linear interpolation between stanp(ipl) and stanp(iph)
        iph = 2
        do while ( p(jz) > stanp(iph) )
           iph = iph + 1
        end do
        ipl = iph - 1

        x1 = exp( coefks(1,ipl) + coefks(2,ipl)*ztfact + coefks(3,ipl)*ztfact2 )
        x2 = exp( coefks(1,iph) + coefks(2,iph)*ztfact + coefks(3,iph)*ztfact2 )
        fkg(jz) = x1 + ( x2 - x1 ) / ( stanp(iph) - stanp(ipl) ) &
                 * ( p(jz) - stanp(ipl) )
     endif
  end do

end subroutine qks

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine qki ( coefki, fkg )
!
! Description :
! -----------
!    Calculation of the absorption coefficients for ir spectral bands.
!
! fkg(nrlev) are the gaseous absorption coefficients in units of (cm-atm)**-1
! for a given cumulative probability in nrlev layers.
! coefki(3,19) are the coefficients to calculate the absorption coefficient
! at the temperature t for the 19 pressures by
!         ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
! and the absorption coefficient at  conditions  other  than  those 19
! pressures is interpolated linearly with pressure (Fu, 1991).


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  !                       Declarations and implicit none
  !
  ! Oct-2017  Josue Bock  Fortran90
  !                       minor improvements of code

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Array arguments with intent(in):
  real(kind=dp), intent(in) :: coefki(3,19) ! tabulated coefficients
! Array arguments with intent(out):
  real(kind=dp), intent(out) :: fkg(nrlev)  ! gaseous absorption coefficients

! Local parameters:
  real (kind=dp), parameter :: ztref = 245.0_dp

! Local scalars:
  integer :: ipl, iph                ! index low, high such that stanp(ipl) < p <stanp(iph)
  integer :: jz                      ! vertical loop index (top-down)
  real (kind=dp) :: x1, x2           ! b and c of the above formula
  real (kind=dp) :: ztfact, ztfact2  ! temperature factor, squared

! Local arrays:
  real(kind=dp), parameter :: stanp(19) = &
       (/ 25.1_dp,   39.8_dp,   63.1_dp,   100._dp,   158._dp, &
          251._dp,   398._dp,   631._dp,  1000._dp,  1580._dp, &
         2510._dp,  3980._dp,  6310._dp, 10000._dp, 15800._dp, &
        25100._dp, 39800._dp, 63100._dp,100000._dp /)

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

! == End of declarations =======================================================

  do jz=1, nrlev
     ztfact  = t(jz) - ztref
     ztfact2 = ztfact * ztfact

     if ( p(jz) <= stanp(1) ) then
        ! pressure lower than 25.1 Pa: linear interpolation between 0. and the first value
        x1 = exp( coefki(1,1) + coefki(2,1)*ztfact + coefki(3,1)*ztfact2 )
        fkg(jz) = x1 * p(jz) / stanp(1)

     elseif ( p(jz) >= stanp(19) ) then
        ! pressure higher than 1000 hPa: extrapolate the slope of the latest two
        ! tabulated values
        x1 = exp( coefki(1,18) + coefki(2,18)*ztfact + coefki(3,18)*ztfact2 )
        x2 = exp( coefki(1,19) + coefki(2,19)*ztfact + coefki(3,19)*ztfact2 )
        fkg(jz) = x1 + ( x2 - x1 ) / ( stanp(19) - stanp(18) ) &
                 * ( p(jz) - stanp(18) )

     else
        ! other cases: linear interpolation between stanp(ipl) and stanp(iph)
        iph = 2
        do while ( p(jz) > stanp(iph) )
           iph = iph + 1
        end do
        ipl = iph - 1

        x1 = exp( coefki(1,ipl) + coefki(2,ipl)*ztfact + coefki(3,ipl)*ztfact2 )
        x2 = exp( coefki(1,iph) + coefki(2,iph)*ztfact + coefki(3,iph)*ztfact2 )
        fkg(jz) = x1 + ( x2 - x1 ) / ( stanp(iph) - stanp(ipl) ) &
                 * ( p(jz) - stanp(ipl) )
     endif
  end do

end subroutine qki

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine qkio3 ( coefki, fkg )
!
! Description :
! -----------
!    Calculation of the absorption coefficients for ozone in spectral band 12
!
! fkg(nrlev) are the gaseous absorption coefficients in units of (cm-atm)**-1
! for a given cumulative probability in nrlev layers.
! coefki(3,19) are the coefficients to calculate the absorption coefficient
! at the temperature t for the 19 pressures by
!         ln k = a + b * ( t - 250 ) + c * ( t - 250 ) ** 2
! and the absorption coefficient at  conditions  other  than  those 19
! pressures is interpolated linearly with pressure (Fu, 1991).


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  !                       Declarations and implicit none
  !
  ! Oct-2017  Josue Bock  Fortran90
  !                       minor improvements of code

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Array arguments with intent(in):
  real (kind=dp), intent(in) :: coefki(3,19)
! Array arguments with intent(out):
  real (kind=dp), intent(out) :: fkg(nrlev)

! Local parameters:
  real (kind=dp), parameter :: ztref = 250.0_dp

! Local scalars:
  integer :: ipl, iph                ! index low, high such that stanp(ipl) < p <stanp(iph)
  integer :: jz                      ! vertical loop index (top-down)
  real (kind=dp) :: x1, x2           ! b and c of the above formula
  real (kind=dp) :: ztfact, ztfact2  ! temperature factor, squared

! Local arrays:
  real(kind=dp), parameter :: stanp(19) = &
       (/ 25.1_dp,   39.8_dp,   63.1_dp,   100._dp,   158._dp, &
          251._dp,   398._dp,   631._dp,  1000._dp,  1580._dp, &
         2510._dp,  3980._dp,  6310._dp, 10000._dp, 15800._dp, &
        25100._dp, 39800._dp, 63100._dp,100000._dp /)

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

! == End of declarations =======================================================

  do jz=1, nrlev
     ztfact  = t(jz) - ztref
     ztfact2 = ztfact * ztfact

     if ( p(jz) <= stanp(1) ) then
        ! pressure lower than 25.1 Pa: linear interpolation between 0. and the first value
        x1 = exp( coefki(1,1) + coefki(2,1)*ztfact + coefki(3,1)*ztfact2 )
        fkg(jz) = x1 * p(jz) / stanp(1)

     elseif ( p(jz) >= stanp(19) ) then
        ! pressure higher than 1000 hPa: extrapolate the slope of the latest two
        ! tabulated values
        x1 = exp( coefki(1,18) + coefki(2,18)*ztfact + coefki(3,18)*ztfact2 )
        x2 = exp( coefki(1,19) + coefki(2,19)*ztfact + coefki(3,19)*ztfact2 )
        fkg(jz) = x1 + ( x2 - x1 ) / ( stanp(19) - stanp(18) ) &
                 * ( p(jz) - stanp(18) )

     else
        ! other cases: linear interpolation between stanp(ipl) and stanp(iph)
        iph = 2
        do while ( p(jz) > stanp(iph) )
           iph = iph + 1
        end do
        ipl = iph - 1

        x1 = exp( coefki(1,ipl) + coefki(2,ipl)*ztfact + coefki(3,ipl)*ztfact2 )
        x2 = exp( coefki(1,iph) + coefki(2,iph)*ztfact + coefki(3,iph)*ztfact2 )
        fkg(jz) = x1 + ( x2 - x1 ) / ( stanp(iph) - stanp(ipl) ) &
                 * ( p(jz) - stanp(ipl) )
     endif
  end do

end subroutine qkio3

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine qopo3s ( fk, tg )
!
! Description :
! -----------
!     Calculation of ozone absorption in band 1 ( 50000 - 14500 cm**-1 )


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  !                       Declarations and implicit none
  !
  ! Oct-2017  Josue Bock  Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlay,               &
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real(kind=dp), intent(in) :: fk
! Array arguments with intent(out):
  real(kind=dp), intent(out) :: tg(nrlay)

! Local scalars:
  real(kind=dp) :: fq
  integer :: jz

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

  common /ozon/ qmo3(nrlev)
  real(kind=dp) :: qmo3

! == End of declarations =======================================================

  fq = 2.3808_dp * fk
  do jz = 1,nrlay
     tg(jz) = ( qmo3(jz) + qmo3(jz+1) ) * ( p(jz+1) - p(jz) ) * fq
  end do

end subroutine qopo3s

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine qoph2o ( fkg, tg )
!
! Description :
! -----------
!    Calculation of optical depth of water vapour absorption in bands
!    2 ( 14500 - 7700 cm**-1 ), 3  ( 7700 - 5250 cm**-1 ),
!    4  ( 5250 - 4000 cm**-1 ), 5  ( 4000 - 2850 cm**-1 ),
!    6  ( 2850 - 2500 cm**-1 ), 7  ( 2200 - 1900 cm**-1 ),
!    8  ( 1900 - 1700 cm**-1 ), 9  ( 1700 - 1400 cm**-1 ),
!    10 ( 1400 - 1250 cm**-1 ), 11 ( 1250 - 1100 cm**-1 ),
!    12 ( 1100 -  980 cm**-1 ), 13 (  980 -  800 cm**-1 ),
!    16 (  540 -  400 cm**-1 ), 17 (  400 -  280 cm**-1 ),
!    and 18 ( 280 - 0 cm**-1 ).


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  ! Oct-2016  Josue Bock  Declarations and implicit none
  !
  ! Oct-2017  Josue Bock  Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlay,               &
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real(kind=dp), intent(in) :: fkg(nrlev)
! Array arguments with intent(out):
  real(kind=dp), intent(out) :: tg(nrlay)

! Local scalars:
  integer :: jz

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

! == End of declarations =======================================================

  do jz = 1,nrlay
     tg(jz) = ( fkg(jz) * xm1(jz) + fkg(jz+1) * xm1(jz+1) ) &
            * ( p(jz+1) - p(jz) ) * 6.349205_dp
  end do

end subroutine qoph2o

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine qopch4 ( fkg, tg )
!
! Description :
! -----------
!   Calculation of optical depth of methane (CH4) absorption in bands
!   10 ( 1400 - 1250 cm**-1 ) and 11 ( 1250 - 1100 cm**-1 )


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  ! Oct-2016  Josue Bock  Declarations and implicit none
  !
  ! Oct-2017  Josue Bock  Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlay,               &
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real(kind=dp), intent(in) :: fkg(nrlev)
! Array arguments with intent(out):
  real(kind=dp), intent(out) :: tg(nrlay)

! Local scalars:
  integer :: jz

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

! == End of declarations =======================================================

  do jz = 1,nrlay
     tg(jz) = ( fkg(jz)+fkg(jz+1) ) * ( p(jz+1)-p(jz) ) * 6.3119e-6_dp
  end do

end subroutine qopch4

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine qopn2o ( fkg, tg )
!
! Description :
! -----------
!    Calculation of optical depth of nitrogene oxide (N2O) absorption in bands
!    10 ( 1400 - 1250 cm**-1 ) and 11 ( 1250 - 1100 cm**-1 )


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  !                       Declarations and implicit none
  !
  ! Oct-2017  Josue Bock  Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlay,               &
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real(kind=dp), intent(in) :: fkg(nrlev)
! Array arguments with intent(out):
  real(kind=dp), intent(out) :: tg(nrlay)

! Local scalars:
  integer :: jz

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

! == End of declarations =======================================================

  do jz = 1,nrlay
     tg(jz) = ( fkg(jz)+fkg(jz+1) ) * ( p(jz+1)-p(jz) ) * 1.10459e-6_dp
  end do

end subroutine qopn2o
!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine qopo3i ( fkg, tg )
!
! Description :
! -----------
!    Calculation of optical depth of ozone absorption in band
!    12 ( 1100 - 980 cm**-1 ).


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  !                       Declarations and implicit none
  !
  ! Oct-2017  Josue Bock  Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlay,               &
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Array arguments with intent(in):
  real(kind=dp), intent(in) :: fkg(nrlev)
! Array arguments with intent(out):
  real(kind=dp), intent(out) :: tg(nrlay)

! Local scalars:
  integer :: jz

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

  common /ozon/ qmo3(nrlev)
  real(kind=dp) :: qmo3

! == End of declarations =======================================================

  do jz = 1,nrlay
     tg(jz) = ( fkg(jz) * qmo3(jz) + fkg(jz+1) * qmo3(jz+1) ) &
            * ( p(jz+1) - p(jz) ) * 2.3808_dp
  end do

end subroutine qopo3i

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine qophc ( fkg, tg )
!
! Description :
! -----------
!    Calculation of optical depth for overlapping absorption of water vapour and
!    carbon dioxide (CO2) in bands
!    14 ( 800 - 670 cm**-1) and 15 ( 670 - 540 cm**-1).
!    See page 86 of Fu (1991).
!

! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  !                       Declarations and implicit none
  !
  ! Oct-2017  Josue Bock  Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlay,               &
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Array arguments with intent(in):
  real(kind=dp), intent(in) :: fkg(nrlev)
! Array arguments with intent(out):
  real(kind=dp) , intent(out) :: tg(nrlay)

! Local scalars:
  integer :: jz

! Common blocks:
  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

! == End of declarations =======================================================

  do jz = 1,nrlay
     tg(jz) = ( fkg(jz) + fkg(jz+1) ) * ( p(jz+1) - p(jz) ) * 0.005_dp
  end do

end subroutine qophc

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine tau
!
! Description :
! -----------
!   Calculation of total optical depth {dtau}, single scattering albedo {om}
!   and Legendre coefficients {pl} for current spectral band {ib}
!   and cumulative probability {ig}.
!   The first index (second for {pl}) stands for cloud free (1) and
!   cloudy (2) parts.


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  !
  ! Oct-2016  Josue Bock  Removed unused argument (ig)
  !                       Removed NO2 in the first wavelength band
  !                          ('tauno2' was set to 0 anyway)
  !                          (see Loughlin et al., QJRMS 1997: the absorption of NO2 is not in PIFM2)
  !
  ! Oct-2017  Josue Bock  Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       nrlay

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  integer :: jl                  ! Loop index for Legendre coefficients
  integer :: jz                  ! Loop index in the vertical dimension (top-down)
  real(kind=dp) :: zsum1, zsum2
  real(kind=dp) :: zx1_na
  real(kind=dp) :: zf

! Common blocks:
  common /con/ tgcon(nrlay)
  real(kind=dp) :: tgcon

  common /extra/ waer(nrlay),taer(nrlay),plaer(2,nrlay)
  real(kind=dp) :: waer, taer, plaer

  common /gas/ tg(nrlay)
  real(kind=dp) :: tg

  common /h2o/ t2w(nrlay), w2w(nrlay), pl2w(2,nrlay)
  real(kind=dp) :: t2w, w2w, pl2w

  common /opohne/ dtau(2,nrlay),om(2,nrlay),pl(2,2,nrlay)
  real(kind=dp) :: dtau, om, pl

  common /ray/ dtaur(nrlay), plr(2,nrlay)
  real(kind=dp) :: dtaur, plr

! == End of declarations =======================================================

! total optical depth {dtau}
  do jz=1,nrlay
     dtau(1,jz)=dtaur(jz)+taer(jz)+tgcon(jz)+tg(jz)
     dtau(2,jz)=dtau(1,jz)+t2w(jz)
  enddo

  do jz=1,nrlay
     zx1_na= taer(jz)*waer(jz)
     zsum1=dtaur(jz) + zx1_na
     zsum2=zsum1 + t2w(jz)*w2w(jz)
! total single scattering albedo {om}
     if(dtau(1,jz) > 1.0e-20_dp) then
        om(1,jz)=zsum1/dtau(1,jz)
        om(2,jz)=zsum2/dtau(2,jz)
     else
        om(1,jz)=0.0_dp
        om(2,jz)=0.0_dp
     endif
! Legendre-coefficients for total scattering {pl}
     do jl=1,2
        zf=dtaur(jz)*plr(jl,jz) + zx1_na*plaer(jl,jz)
        if(zsum1 >= 1.e-20_dp) then
           pl(jl,1,jz)=zf/zsum1
           pl(jl,2,jz)=(zf+t2w(jz)*w2w(jz)*pl2w(jl,jz))/zsum2
        else
           pl(jl,1,jz)=0.0_dp
           pl(jl,2,jz)=0.0_dp
        endif
     enddo
  enddo

end subroutine tau

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine kurzw(ib,u0)
!
! Description :
! -----------
!   Solution of the radiative transfer equation for current spectral band
!   (solar -> ib = 1 ... 6) and cumulative probability {ig}.
!
!   The direct downward solar fluxes are calculated twice:
!   -- sf and sw are the direct downward solar fluxes without the forward
!      scattering peak of the Delta-Eddington method
!   -- ssf and ssw are the direct downward solar fluxes with Delta-Eddington
!   ('f' and 'w' stand for cloud free and cloudy parts)
!   The arrays f2f, f2w, f1f, f1w are the diffuse downward (2) and upward (1)
!   solar fluxes.
!
!   Mostly based on the following reference (and references therein):
!   Zdunkowski et al. 1982, Contrib. Atmos. Phys. Vol. 55 (3), pp. 215-238


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  ! Nov-2016  Josue Bock  ua ... wd changed from 1D aray (nrlay) to scalars
  !                       missing declarations, implicit none
  !                       added u0red to avoid changing a global variable
  !                       further bugfix related to u0red: previously, u0
  !                          was changed within this SR, but not properly
  !                          reset for all other cases. It has to be used
  !                          ONLY in the resonance case!
  !
  ! Oct-2017  Josue Bock  Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! External subroutine
       abortM

  USE file_unit, ONLY :     &
! Imported Parameters:
       jpfunerr

  USE global_params, ONLY : &
! Imported Parameters:
       mbs,                 &
       mbir,                &
       nrlay,               &
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Array arguments with intent(in):
  integer, intent(in) :: ib                         ! current spectral band
  real(kind=dp), intent(in) :: u0                ! cosine (solar zenith angle)

! Local parameters:
  real(kind=dp), parameter :: u=2._dp            ! diffusivity factor = reciprocal of the mean effective cos(SZA)
  real(kind=dp), parameter :: delu0=0.001_dp     ! resonance case correction
  real(kind=dp), parameter :: reson=0.1e-6_dp    ! boundary value for resonance case
  real(kind=dp), parameter :: absfr=0.001_dp     ! boundary value for no absorption case
  real(kind=dp), parameter :: strfr=0.03_dp      ! boundary value for no scattering case
  real(kind=dp), parameter :: p1ray=0.1_dp       ! boundary value for ray. scattering

! Local scalars:
  real(kind=dp) :: ua, ub, uc, ud                   ! flux modification for clouds
  real(kind=dp) :: va, vb, vc, vd
  real(kind=dp) :: wa, wb, wc, wd
  real(kind=dp) :: ak                               ! ratio of absorption and extinction
  real(kind=dp) :: alph1, alph2, alph3, alph4       ! coefficients 'alpha'
  real(kind=dp) :: u0kw                             !
  real(kind=dp) :: u0red                            ! check of resonance case
  real(kind=dp) :: f, emf, emfkw                    ! asymmetry factor of phase function
  real(kind=dp) :: b0, bu0                          ! backscattering coefficient
  real(kind=dp) :: dtu0, dtu, eps2, eps, omf, emomf ! equation factors
  real(kind=dp) :: u02, ueps2, emu, e, e2, m, m2
  real(kind=dp) :: e2m2, ouf, te, u0a1, u0a2
  real(kind=dp) :: gam1, gam2, g1a1, da

  integer :: jc       ! loop index, cloud free (1) or cloudy (2)
  integer :: jz, jzp  ! loop indexes, vertical (top-down)

! Local arrays:
  real(kind=dp) :: a1(2,nrlay),a2(2,nrlay),a3(2,nrlay),a6(2,nrlay)   ! matrix coefficients (local)

! Common blocks:
  common /leck1/ a4(2,nrlay),a5(2,nrlay)                             ! matrix coefficients
  real(kind=dp) :: a4, a5

  common /leck2/ sf(nrlev),sw(nrlev),ssf(nrlev),ssw(nrlev), &        ! radiation fluxes
                 f2f(nrlev),f2w(nrlev),f1f(nrlev),f1w(nrlev)
  real(kind=dp) :: sf, sw, ssf, ssw, f2f, f2w, f1f, f1w

  common /opohne/ dtau(2,nrlay),om(2,nrlay),pl(2,2,nrlay)            ! optical variables
  real(kind=dp) :: dtau, om, pl

  common /part/ cc(4,nrlay),bb(4,nrlay)                              ! cloudiness (continuity factors)
  real(kind=dp) :: cc, bb

  common /tmp2/ as(mbs),ee(mbir)                                     ! albedo and emissivity (unused here)
  real(kind=dp) :: as, ee

! == End of declarations =======================================================


! Definition of multiply used variables
  u0kw  = 1._dp / u0

  do jz = 1,nrlay

     ! Loop for cloud free (index 1) and cloudy (index 2) parts
     !---------------------------------------------------------
     do jc = 1,2

! Coefficients a for matrix solution of the radiative transfer equation.
!-----------------------------------------------------------------------
!    the coefficients a1 without Delta-Eddington are set to a6( *, *).


! Case 1: no extinction
        if(dtau(jc,jz) <= 1.e-7_dp) then
           a1(jc,jz) = 1._dp
           a2(jc,jz) = 0._dp
           a3(jc,jz) = 0._dp
           a4(jc,jz) = 1._dp
           a5(jc,jz) = 0._dp
           a6(jc,jz) = 1._dp

        else
           ! General definitions for all other cases
           dtu0      = dtau(jc,jz) * u0kw
           a6(jc,jz) = exp(-min(dtu0,7.5e1_dp))
           dtu       = dtau(jc,jz) * u

! Case 2: no scattering (om < strfr)
           if(om(jc,jz) < strfr) then
              a1(jc,jz) = a6(jc,jz)
              a2(jc,jz) = 0._dp
              a3(jc,jz) = 0._dp
              a4(jc,jz) = exp(-min(dtu,7.5e1_dp))
              a5(jc,jz) = 0._dp

           else
              ! General definitions for all other cases
              !
              ! backscattering coefficients: beta-0 (here b0), beta(u0) (here bu0)
              ! In the case of pure rayleigh scattering (p1 < p1ray) they
              ! are set to 0.5.
              ! p1 is the first Legendre-coefficient of the scattering function.
              ! f is the asymmetry factor of Delta-Eddington
              !
              ! Caution: The jump from b0 and bu0 to 0.5 in the case of p1 < p1ray
              ! is not satisfying. It is an interim solution.                        ! jjb NOK
              ak    = 1._dp - om(jc,jz)
              b0    = 0.5_dp
              bu0   = 0.5_dp
              f     = pl(2,jc,jz) / 5._dp
              emf   = 1._dp - f
              emfkw = 1._dp / emf
              if(pl(1,jc,jz) >= p1ray) then
                 b0  = (3.0_dp - pl(1,jc,jz)) / 8._dp
                 bu0 = 0.5_dp - u0/4._dp * (pl(1,jc,jz) - 3._dp*f) * emfkw
              end if

! Case 3: no absorption = pure scattering (ak < absfr, om = 1 )
              if(ak < absfr) then
                 alph1 = u*b0
                 alph3 = bu0 ! jjb NOK 1-f missing
                 alph4 = 1._dp - alph3
                 gam1  = alph3-alph1*u0*emfkw ! jjb NOK
                 a1(jc,jz) = exp(-min(dtu0*emf,7.5e1_dp))
                 a4(jc,jz) = 1._dp / (1._dp + alph1*dtau(jc,jz))
                 a2(jc,jz) = a4(jc,jz) * (1._dp - gam1 * (1._dp - a1(jc,jz))) - a1(jc,jz) ! jjb NOK
                 a3(jc,jz) = 1._dp - a1(jc,jz) - a2(jc,jz) ! jjb NOK
                 a5(jc,jz) = 1._dp - a4(jc,jz)

! Case 4: absorption and scattering (normal case)
              else
                 alph2 = u * b0 * om(jc,jz)
                 alph1 = u * ak + alph2
                 alph3 = bu0 * om(jc,jz) ! jjb NOK 1-f missing
                 alph4 = om(jc,jz) - alph3 ! jjb NOK 1-f missing
                 eps2  = alph1**2 - alph2**2
                 eps   = sqrt(eps2)
                 omf   = om(jc,jz) * f
                 emomf = 1._dp - omf

                 u02   = u0**2
                 ueps2 = u02 * eps2
                 emu   = emomf**2 - ueps2

                 ! First initialisation of u0red, used in resonance case
                 u0red = u0

! Resonance case if the sun angle u0 is in resonance
! with the reference angle u:
! -> u0 is reduced by delu0
                 do while(abs(emu) <= reson)
                    u0red = u0red - delu0
                    if(u0red <= 0.) then
                       write(jpfunerr,*)'SR kurzw: Error in resonance-case correction:' &
                               //' The sun is moved back below the horizon.' &
                               //' Thus it is nighttime and shortwave fluxes are zero'
                       call abortM ('Error in subroutine kurzw.')
                    end if

                    u02   = u0red**2
                    ueps2 = u02*eps2
                    emu   = emomf**2 - ueps2
                 end do

                 a1(jc,jz) = exp(-min(dtu0*emomf, 7.5e1_dp))
                 e         = exp(-min(dtau(jc,jz)*eps, 7.5e1_dp))
                 m         = alph2 / (alph1+eps)
                 e2        = e**2
                 m2        = m**2
                 e2m2      = e2 * m2
                 ouf       = 1._dp / (1._dp - e2m2)
                 a4(jc,jz) = e * (1._dp - m2) * ouf
                 a5(jc,jz) = m * (1._dp - e2) * ouf
                 te        = emf / emu
                 u0a1      = u0red * alph1
                 u0a2      = u0red * alph2
                 gam1      =  (alph3*(emomf-u0a1) - u0a2*alph4)*te
                 gam2      = -(alph4*(emomf+u0a1) + u0a2*alph3)*te
                 g1a1      = gam1*a1(jc,jz)
                 da        = a1(jc,jz) - a4(jc,jz)
                 a2(jc,jz) =  gam2*da - a5(jc,jz)*g1a1
                 a3(jc,jz) = -gam2*a5(jc,jz) - a4(jc,jz)*g1a1+gam1
              end if
           end if
        end if
     end do
  end do

!------------------------------------------------------------------------------

! Calculation of the parallel fluxes and the right hand side of the
! diffuse system A * F  =  R0.

  sf(1)  = u0
  sw(1)  = 0._dp
  ssf(1) = sf(1)
  ssw(1) = 0._dp
  ua     = bb(1,1) * ssf(1)
  ub     = ssf(1) - ua
  f2f(1) = 0._dp
  f2w(1) = 0._dp
  ! cloud free part
  ssf(2) = a1(1,1)*ua
  f2f(2) = a2(1,1)*ua
  f1f(1) = a3(1,1)*ua
  sf(2)  = a6(1,1)*ua
  ! cloudy part
  ssw(2) = a1(2,1)*ub
  f2w(2) = a2(2,1)*ub
  f1w(1) = a3(2,1)*ub
  sw(2)  = a6(2,1)*ub
  do jz = 2,nrlay
     jzp = jz+1

     ua = bb(1,jz) * ssf(jz)
     ub = ssf(jz) - ua
     uc = bb(1,jz) * sf(jz)
     ud = sf(jz) - uc
     va = cc(3,jz) * ssw(jz)
     vb = ssw(jz) - va
     vc = cc(3,jz) * sw(jz)
     vd = sw(jz) - vc
     wa = ua + va
     wb = ub + vb
     wc = uc + vc
     wd = ud + vd
     ! cloud free part
     ssf(jzp) = a1(1,jz) * wa
     f2f(jzp) = a2(1,jz) * wa
     f1f(jz)  = a3(1,jz) * wa
     sf(jzp)  = a6(1,jz) * wc
     ! cloudy part
     ssw(jzp) = a1(2,jz) * wb
     f2w(jzp) = a2(2,jz) * wb
     f1w(jz)  = a3(2,jz) * wb
     sw(jzp)  = a6(2,jz) * wd
  end do

! Last two elements of the vector R0
!
  f1f(nrlev) = as(ib) * ssf(nrlev)
  f1w(nrlev) = as(ib) * ssw(nrlev)

end subroutine kurzw

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine langw(ib)
!
! Description :
! -----------
!   Solution of the radiative transfer equation for current spectral band
!   (ir -> ib = 7 ... 18) and cumulative probability {ig}.
!
!   The arrays f2f, f2w, f1f, f1w are the diffuse downward (2) and upward (1)
!   ir fluxes. ('f' and 'w' stand for cloud free and cloudy parts)


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  ! Nov-2016  Josue Bock  missing declarations, implicit none
  !                       added an undetermination case for a6 calculation
  !                       reorganised a(4:6,...) ==> a4(...) - a6(...)
  ! Oct-2017  Josue Bock  Fortran90,
  !                       replaced a test /= 0. by abs(...) >= tiny_dp

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       mbs,                 &
       mbir,                &
       nrlay,               &
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp,                  &
       tiny_dp

  implicit none

! Subroutine arguments
! Array arguments with intent(in):
  integer, intent(in) :: ib                         ! current spectral band

! Local parameters:
  real(kind=dp), parameter :: u=1.66_dp             ! diffusivity factor = reciprocal of the mean effective cos(SZA)

! Local scalars:
  real(kind=dp) :: agdb, ak, alph1, alph2, at
  real(kind=dp) :: b0, db, dtu
  real(kind=dp) :: e, eps, epstau, eq
  real(kind=dp) :: ha, hb
  real(kind=dp) :: rm, rmq, rn
  integer :: i, ip, ii0, l                             ! loop indexes

! Local arrays
  real(kind=dp) :: a6(2,nrlay)                         ! matrix coefficients (local)

! Common blocks:
  common /cb09/ frac(nrlay),rew(nrlay),rho2w(nrlay)
  real (kind=dp) :: frac, rew, rho2w

  common /leck1/ a4(2,nrlay),a5(2,nrlay)                       ! matrix coefficients
  real(kind=dp) :: a4, a5

  common /leck2/ sf(nrlev),sw(nrlev),ssf(nrlev),ssw(nrlev), &  ! radiation fluxes
                 f2f(nrlev),f2w(nrlev),f1f(nrlev),f1w(nrlev)
  real(kind=dp) :: sf, sw, ssf, ssw, f2f, f2w, f1f, f1w

  common /opohne/ dtau(2,nrlay),om(2,nrlay),pl(2,2,nrlay)      ! optical variables
  real(kind=dp) :: dtau, om, pl

  common /part/ cc(4,nrlay),bb(4,nrlay)                        ! cloudiness (continuity factors)
  real(kind=dp) :: cc, bb

  common /planci/ pib(nrlev),pibs                              ! black body radiation
  real(kind=dp) :: pib, pibs

  common /tmp2/ as(mbs),ee(mbir)                               ! albedo (unused here) and emissivity
  real(kind=dp) :: as, ee

! == End of declarations =======================================================


  do i=1,nrlay               ! Index i goes from top to bottom

     ! Loop for cloud free (index 1) and cloudy parts (index 2)
     !---------------------------------------------------------
     do l=1,2

!
! Case 1: no extinction
!
        if(dtau(l,i) <= 1.e-7_dp) then
           a4(l,i) = 1._dp
           a5(l,i) = 0._dp
           a6(l,i) = 1._dp

!
! Case 2: no scattering (om < 1.0e-7)
!
        else if(om(l,i) <= 1.e-7_dp) then
           dtu = dtau(l,i) * u
           a4(l,i) = exp(-dtu)
           a5(l,i) = 0._dp
           a6(l,i) = (1._dp - a4(l,i)) / dtu

!
! Case 3: no absorption (ak < 1.0e-7)
!
        else
           ak = 1._dp - om(l,i)
           b0 = (3._dp - pl(1,l,i)) / 8._dp
           alph1 = u * (1._dp - (1._dp - b0) * om(l,i))
           alph2 = u * b0 * om(l,i)

           if(ak <= 1.e-7_dp) then
              at = alph1 * dtau(l,i)
              a4(l,i) = 1._dp / (1._dp + at)
              a5(l,i) = a4(l,i) * at
              a6(l,i) = 0._dp

!
! Case 4: absorption and scattering (normal case)
!
           else
              eps = sqrt(alph1**2 - alph2**2)
              epstau = eps * dtau(l,i)
              e = 0._dp
              if(epstau < 87._dp) e=exp(-epstau)
              rm  = alph2 / (alph1+eps)
              eq  = e**2
              rmq = rm**2
              rn  = 1._dp - eq*rmq
              a4(l,i) = e * (1._dp - rmq) / rn
              a5(l,i) = rm * (1._dp - eq) / rn
              if (abs(alph1+alph2) >= tiny_dp) then ! jjb see [Z82] p. 219: indetermination case
                 a6(l,i) = (1._dp - a4(l,i) - a5(l,i)) / &
                      ((alph1+alph2) * dtau(l,i))
              else
                 a6(l,i) = 1.0_dp
              end if
           end if
        end if
     end do
  end do

!------------------------------------------------------------------------------

! Right hand side of the equation system  A * (pi*B - F)  =  LK(a6,pi*B)
! set to the arrays f1f, f1w, f2f, f2w, where the solution vectors
! can be found after the matrix solution.

!
! uppermost layer
  f2f(1) = pib(1)
  f2w(1) = 0._dp
!
! regular layers
  do i = 1,nrlay
     ip = i+1
     db = pib(i)-pib(ip)
     f1f(i)  = (1._dp-frac(i))*a6(1,i)*db
     f1w(i)  = frac(i)*a6(2,i)*db
     f2f(ip) = -f1f(i)
     f2w(ip) = -f1w(i)
  end do
!
! surface
  ii0 = ib-mbs
  agdb = ee(ii0) * (pib(nrlev)-pibs) + (1._dp - ee(ii0)) * &
      (1._dp - ee(ii0)) * (pib(nrlev) - pib(nrlay)) * a6(1,nrlay) * &
      (1._dp - frac(nrlay))
  f1w(nrlev) = agdb * frac(nrlay)
  f1f(nrlev) = agdb - f1w(nrlev)
!
! upper boundary condition for the right hand side of the third to sixth
! equation
  ha = bb(1,1)*f2f(1)
  hb = f2f(1)-ha
  f2f(2) = f2f(2)+a4(1,1)*ha
  f1f(1) = f1f(1)+a5(1,1)*ha
  f2w(2) = f2w(2)+a4(2,1)*hb
  f1w(1) = f1w(1)+a5(2,1)*hb

end subroutine langw

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
subroutine jeanfr (ib)
!
! Description :
! -----------
!   Solution of the linear equation system A(diffus) * F  =  R0.
!   The right hand side is calculated in KURZW or LANGW.


! Modifications :
! -------------
  ! Jul-2016  Josue Bock  Header including "USE ... ONLY"
  !                       Removal of labeled do-loops
  ! Nov-2016  Josue Bock  Changed td from 2-D to 1-D
  ! Oct-2017  Josue Bock  Fortran90

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       mbs,                 &
       mbir,                &
       nrlay,               &
       nrlev

  USE precision, ONLY :     &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Array arguments with intent(in):
  integer, intent(in) :: ib

! Local scalars:
  real(kind=dp) :: ae
  real(kind=dp) :: fa
  real(kind=dp) :: ga,gb,gc,gd,ge,gf
  real(kind=dp) :: ha,hb,hc,hd
  real(kind=dp) :: tds1,tds2,tds3,tus1
  integer :: jz,jzm,jzp                       ! Loop indexes

! Local arrays:
  real(kind=dp) :: tu(9,nrlay),td(7)

! Common blocks:
  common /leck1/ a4(2,nrlay),a5(2,nrlay)                          ! matrix coefficients
  real(kind=dp) :: a4, a5

  common /leck2/ sf(nrlev),sw(nrlev),ssf(nrlev),ssw(nrlev), &             ! radiation fluxes
                 f2f(nrlev),f2w(nrlev),f1f(nrlev),f1w(nrlev)
  real(kind=dp) :: sf, sw, ssf, ssw, f2f, f2w, f1f, f1w

  common /part/ cc(4,nrlay),bb(4,nrlay)
  real(kind=dp) :: cc, bb

  common /tmp2/ as(mbs),ee(mbir)
  real(kind=dp) :: as, ee

! == End of declarations =======================================================

!   Diffusive System A * F = R0
!
! The right side is provided by subroutines KURZW or LANGW
!
! equation 3. to 6.: save of matrix elements above the main diagonal
! in arrays -tu(k,jz)

  tu(1,1)=0._dp
  tu(2,1)=a4(1,1)*bb(2,1)
  tu(3,1)=a4(1,1)*cc(4,1)
  tu(4,1)=a4(2,1)*cc(2,1)
  tu(5,1)=a4(2,1)*bb(4,1)
  tu(6,1)=a5(1,1)*bb(2,1)
  tu(7,1)=a5(1,1)*cc(4,1)
  tu(8,1)=a5(2,1)*cc(2,1)
  tu(9,1)=a5(2,1)*bb(4,1)

! Blocks of four equations:
! Elimination of matrix elements below the main diagonal, matrix elements
! above the main diagonal are saved in -tu(k, jz).
! Right hand side saved in f2f, f2w, f1f, f1w.
  do jz=2,nrlay
     jzm=jz-1
     jzp=jz+1

     ga=bb(1,jz)*tu(6,jzm)
     gb=tu(6,jzm)-ga
     gc=cc(3,jz)*tu(8,jzm)
     gd=tu(8,jzm)-gc
     ha=ga+gc
     hc=gb+gd
     ga=bb(1,jz)*tu(7,jzm)
     gb=tu(7,jzm)-ga
     gc=cc(3,jz)*tu(9,jzm)
     gd=tu(9,jzm)-gc
     hb=ga+gc
     hd=gb+gd
     ga=bb(1,jz)*f2f(jz)
     ge=f2f(jz)-ga
     gc=cc(3,jz)*f2w(jz)
     gf=f2w(jz)-gc
     gb=ga+gc
     gd=ge+gf
     td(1)=1._dp / (1._dp - a5(1,jz)*ha)
     f1f(jz)=td(1)*(f1f(jz)+a5(1,jz)*gb)
     tu(1,jz)=td(1)*a5(1,jz)*hb
     fa=td(1)*a4(1,jz)
     tu(2,jz)=fa*bb(2,jz)
     tu(3,jz)=fa*cc(4,jz)
     td(2)=a5(2,jz)*hc
     td(3)=1._dp / (1._dp - a5(2,jz)*hd-td(2)*tu(1,jz))
     f1w(jz)=td(3)*(f1w(jz)+a5(2,jz)*gd+td(2)*f1f(jz))
     td(4)=a4(1,jz)*ha
     td(5)=a4(1,jz)*hb+td(4)*tu(1,jz)
     f2f(jzp)=f2f(jzp)+a4(1,jz)*gb+td(4)*f1f(jz)+td(5)*f1w(jz)
     tu(4,jz)=td(3)*(a4(2,jz)*cc(2,jz)+td(2)*tu(2,jz))
     tu(5,jz)=td(3)*(a4(2,jz)*bb(4,jz)+td(2)*tu(3,jz))
     tu(6,jz)=a5(1,jz)*bb(2,jz)+td(4)*tu(2,jz)+td(5)*tu(4,jz)
     tu(7,jz)=a5(1,jz)*cc(4,jz)+td(4)*tu(3,jz)+td(5)*tu(5,jz)
     td(6)=a4(2,jz)*hc
     td(7)=a4(2,jz)*hd+td(6)*tu(1,jz)
     f2w(jzp)=f2w(jzp)+a4(2,jz)*gd+td(6)*f1f(jz)+td(7)*f1w(jz)
     tu(8,jz)=a5(2,jz)*cc(2,jz)+td(6)*tu(2,jz)+td(7)*tu(4,jz)
     tu(9,jz)=a5(2,jz)*bb(4,jz)+td(6)*tu(3,jz)+td(7)*tu(5,jz)
  end do

! Elimination and backsubstitution of the last two equations.
! Calculation with albedo (solar) and 1-emissivity (ir)
  if(ib <= mbs) then
     ae = as(ib)
  else
     ae = 1._dp - ee(ib-mbs)
  endif
  tds1 = 1._dp / (1._dp - ae*tu(6,nrlay))
  f1f(nrlev) = tds1*(f1f(nrlev)+ae*f2f(nrlev))
  tus1 = tds1*ae*tu(7,nrlay)
  tds2 = ae*tu(8,nrlay)
  tds3 = 1._dp / (1._dp - ae*tu(9,nrlay)-tds2*tus1)
  f1w(nrlev) = tds3*(f1w(nrlev)+ae*f2w(nrlev)+tds2*f1f(nrlev))
  f1f(nrlev) = f1f(nrlev)+tus1*f1w(nrlev)

! Now we have a upper-triangle-matrix with elements -tu(k, i) or 0 or 1.
! Resubstitution with results on the arrays f2f, f2w, f1f, f1w
  do jz=nrlay,1,-1
     jzp=jz+1
     f2w(jzp) = f2w(jzp) + tu(8,jz)*f1f(jzp)+tu(9,jz)*f1w(jzp)
     f2f(jzp) = f2f(jzp) + tu(6,jz)*f1f(jzp)+tu(7,jz)*f1w(jzp)
     f1w(jz)  = f1w(jz)  + tu(4,jz)*f1f(jzp)+tu(5,jz)*f1w(jzp)
     f1f(jz)  = f1f(jz)  + tu(2,jz)*f1f(jzp)+tu(3,jz)*f1w(jzp)+tu(1,jz)*f1w(jz)
  end do

end subroutine jeanfr
