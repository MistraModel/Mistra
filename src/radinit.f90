! ======================================================================= !
! Copyright 1996-2017 the Authors
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


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! radinit.f: radiation initialisation and start (SR str)
!
! contains the following subroutines:
!     - radiation
!     - intrad
!     - ipdata
!     - initr
!     - load1
!     - rotate_in
!     - rotate_out
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine radiation (ldinit)
!
! Description:
! -----------
  !    interface between Mistra and the radiative code:
  !      read tabulated data during model initialisation,
  !      then regrid some variables (T, P, ...) and call the radiative code
  !      to calculate radiative fluxes and heating rates
  !
  !    This routine results from the merge of former SR initstr, which was called during initialisation,
  !      and SR str, which was called at the end of the initialisation and during time integration.
  !      The history of modification from both SRs has been merged below.
  !    Note that str is a shortcut for "strahlung", which means "radiation" in German.


! Author :
! ------
  !    Andreas Bott and others (original SRs initstr and str)
  !    Josue Bock              (current version)


! Modifications :
! -------------
  !     ?        Roland      Change to add 'mic' case, and call load1 every 2 minutes only
  !              von Glasow
  !
  !    Jul-2016  Josue Bock  <SR initstr> Removal of labeled do loops.
  !                          <SR str> Common block /cb20/ was missing for 'clouds'
  !                          <SR str> Comments / header
  !                          Use module for parameters
  !                          All explicit declarations and implicit none
  !
  !    Oct-2016  Josue Bock  <SR str> 'clouds' removed after discussion with Andeas Bott
  !
  ! 10-Nov-2017  Josue Bock  <SR initstr> Removal of cb55, which is useless
  !                          <SR initstr> Calling nstrahl twice during initialisation has
  !                                         absolutely no effect, since the only apparent
  !                                         reason was to fill cb55 data, which were not
  !                                         actually used. This simplifies much this subroutine.
  !
  ! 16-Nov-2017  Josue Bock  <SR str> The mic=false case is now directly handled in SR load1.
  !                                     This avoids to duplicate the part of this SR which is
  !                                     common to both cases.
  !                          Then merged initstr and str into this "new" subroutine, renamed to
  !                            highlight the change.

! == End of header =============================================================

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  logical, intent(in) :: ldinit

! == End of declarations =======================================================

  if (ldinit) then
  ! initialisation of radiation code:
  ! --------------------------------
     ! interpolation of radiation coefficients for prognostic aerosols
     call intrad
     ! input data for the radiative code
     call ipdata
     ! define the vertical grid for the radiation code, and load variables (T, P, ...) from the main program
     call initr
  else
     ! load variables (T, P, ...) from the main program
     call load1
     ! if (lmin/2*2.eq.lmin) call load1
     ! difference is nearly unplottable (cloudless case) but saves 20 % CPU time !
  end if

  ! rotate the grid and the variables so that they are all indexed top-down
  call rotate_in (ldinit)

  ! call the radiative code
  call nstrahl

  ! rotate the output variables
  call rotate_out

  ! print the output of the first calculation
  if (ldinit) call profr

end subroutine radiation
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine intrad

! Description :
! -----------
!    Linear interpolation of radiation parameters qabs0, qext0, asym0
!    which are tabulated for droplet radii xw0 (in microns)
!    and percentage part of pure water xa0 on the refractive index m
!      with: m = ma + (mw -ma) *xa0,
!      mw: refractive index pure water, ma: refractive index pure aerosol
!    on the radiation parameters qabs, qext and asym
!    for actual droplet radii xw1 (in microns) and percentage part of
!    pure water xa1 on the refractive index.
!
!    More information about the calculations of the tabulated values qabs0, qext0, asym0
!    using Mie theorie can be found in the following paper:
!    Bott, Sievers & Zdunkowski, J. Atmos. Sci., 47 (18), 2153-2166, 1990
!    Note that in this paper, alpha = 1-xa


! Interface :
! ---------
!    SR intrad is called during initialisation

! Input :
! -----
!    - data files (3 aerosol types, short (visible) and long (IR) bands)
!    - rn (dry particle radius) and rq (total particle radius) from /cb50/

! Output :
! ------
!    - interpolated values qabs, qext, asym in /cb49/

! Externals :
! ---------
!    none


! Method :
! ------
!    Find index iw0 (2 <= iw0 <= nw0) such that xw0(iw0-1) < xw1 <= xw0(iw0)
!    Then define an interpolation factor dx (proportionality factor, 0 <= dx <= 1)
!      such that the interpolated value is dx*f(iw0) + (1-dx)*f(iw0-1) where f is one of the
!      radiation parameter qabs, qext, asym.

!    If xw1 < xw0(1) or xw1 > xw0(nw0), then the boundary value (xw0(1) or xw0(nw0)) will be
!       used without interpolation nor extrapolation (dx=0. or dx=1., respectively). Warning
!       messages are displayed in these specific cases.

!    Similarly, find index ia0 (2 <= ia0 <= na0) such that xa0(ia0-1) < xw1 <= xa0(ia0), and define an
!       interpolation factor dy.
!    Note the different treatment of special cases:
!       - if xa1 < 0, an error message arise and the program is stopped
!       - if xa1 > 1, another error message arise and the program is stopped


! Author :
! ------
!    Andreas Bott


! Modifications :
! -------------
!     ?         Roland       Introduction of special cases to avoid out-of-bounds errors, and deal
!               von Glasow      with particle radius smaller / larger than tabulated values
!
!    Jul-2016   Josue Bock   Header with comments
!                            Use module for parameters
!                            Explicit declaration of all variables, implicit none
!                            Reindexed qabs/qext/asym for efficiency (innermost index is leftmost)
!                            Removed labeled do loops
!                            Rewritten, with more generic cases, error/warning messages, comments and cleaned code structure
!
!    Oct-2016   Josue Bock   Further commenting after personnal communication with A. Bott
!
! 10-Mar-2017   Josue Bock   Improvement of j case (now accept xa1=0, error if xa<0)
!                            Cosmetic reorganisation, new header organisation
!                            Set nkaer as a new global parameter (later replaced by jptaerrad)
!
! 04-Apr-2017   Josue Bock   Fortran90 conversion
!
! 05-Nov-2017   Josue Bock   Removed module directories, replaced by config
!
! 09-Nov-2017   Josue Bock   Final cleaning after second F90 conversion (GitHub version),
!                              rewritten interpolation, coding standards for variable names
!
! 18-Feb-2019   Josue Bock   Stop the program using abortM
!
! 22-Sep-2019   Josue Bock   Print the warning messages only if jaer==1. Remove format statements (coding standards)

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! Imported Parameters:
       cinpdir,      &
! Imported Routines:
       abortM

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunerr,      &
       jpfunout,      &
       jpfunaerrad

  USE global_params, ONLY : &
! Imported Parameters:
       jptaerrad,           &
       nka,                 &
       nkt,                 &
       mb,                  &
       mbs

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local parameters:
  character (len=*), parameter :: cpfmt1='(3e16.8)'                   ! input format to read tabulated files
  character (len=*), parameter :: cpfmt2='(a,2(i3,x,a),2(f9.4,x,a))'  ! output format for warning messages
  integer, parameter :: na0 = 11  ! Number of tabulated values, volume fraction of water
  integer, parameter :: nw0 = 40  ! Number of tabulated values, total aerosol radius

  real (kind=dp), parameter :: xa0(na0) = &  ! volume fraction of pure water (tabulated values) []
       (/0.0_dp, 0.2_dp, 0.4_dp, 0.6_dp, 0.7_dp, 0.8_dp, 0.85_dp, 0.9_dp, 0.95_dp, 0.975_dp, 1.0_dp/)
  real (kind=dp), parameter :: xw0(nw0) = &  ! total radius of the scattering spheres (tabulated values) [micro m]
       (/0.01_dp, 0.0125_dp, 0.015_dp, 0.02_dp, 0.025_dp, 0.03_dp, 0.04_dp, 0.05_dp, 0.06_dp, 0.08_dp, &
         0.1_dp,  0.125_dp,  0.15_dp,  0.2_dp,  0.25_dp,  0.3_dp,  0.4_dp,  0.5_dp,  0.6_dp,  0.8_dp,  &
         1.0_dp,  1.25_dp,   1.5_dp,   2.0_dp,  2.5_dp,   3.0_dp,  4.0_dp,  5.0_dp,  6.0_dp,  8.0_dp,  &
        10.0_dp, 12.5_dp,   15.0_dp,  20.0_dp, 25.0_dp,  30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp, 80.0_dp/)


! Local scalars:
  character (len=len_trim(cinpdir)+11) :: fname
  real (kind=dp) :: dx,dy,dxdy,dxdy1,dx1dy,dx1dy1 ! Interpolation coefficients
  real (kind=dp) :: xa1                           ! Actual percentage of water computed for each particle class (2D spectrum)
  real (kind=dp) :: xw1                           ! Actual total aerosol radius [micro m] for each particle class (2D spectrum)
  integer :: ifun                                 ! File unit number
  integer :: ia0,iw0                              ! upper indexes of tabulated values for interpolation
  integer :: ja0,jw0                              ! loop indexes, tabulated values
  integer :: jaer                                 ! loop index for aerosol type
  integer :: jb                                   ! loop index for spectral band number
  integer :: jka,jkt                              ! loop indexes, 2D microphysical grid

! Local arrays:
  real (kind=dp) :: qabs0(mb,nw0,na0) ! tabulated values of absorption coefficient
  real (kind=dp) :: qext0(mb,nw0,na0) ! tabulated values of extinction coefficient
  real (kind=dp) :: asym0(mb,nw0,na0) ! tabulated values of asymmetry factor


! Common blocks:
  common /cb49/ qabs(mb,nkt,nka,jptaerrad),qext(mb,nkt,nka,jptaerrad),asym(mb,nkt,nka,jptaerrad)
  real (kind=dp) :: qabs,qext,asym

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)                      ! only rn (dry radius), rq (total radius) are used here
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

! == End of declarations =======================================================


! Open input data files
! ---------------------

! first unit number to read
  ifun=jpfunaerrad

! input qabs, qext and asym for radiation parameters of particles
!   *kw.dat: short wave data
!   *lw.dat: long wave data
  fname=TRIM(cinpdir)//"urbankw.dat"
  open (unit=ifun, file=fname, status='old')
  fname=TRIM(cinpdir)//"urbanlw.dat"
  open (unit=ifun+1, file=fname, status='old')
  fname=TRIM(cinpdir)//"ruralkw.dat"
  open (unit=ifun+2, file=fname, status='old')
  fname=TRIM(cinpdir)//"rurallw.dat"
  open (unit=ifun+3, file=fname, status='old')
  fname=TRIM(cinpdir)//"ozeankw.dat"
  open (unit=ifun+4, file=fname, status='old')
  fname=TRIM(cinpdir)//"ozeanlw.dat"
  open (unit=ifun+5, file=fname, status='old')


  do jaer=1,jptaerrad  ! jaer=1 urban, jaer=2 rural, jaer=3 ocean

! Read input data file: tabulated values
! --------------------------------------
     do ja0=1,na0
        do jw0=1,nw0
           do jb=1,mbs
              read(ifun,cpfmt1) qabs0(jb,jw0,ja0),qext0(jb,jw0,ja0),asym0(jb,jw0,ja0)
           enddo
        enddo
     enddo
     close (ifun)
     ifun=ifun+1
     do ja0=1,na0
        do jw0=1,nw0
           do jb=mbs+1,mb
              read(ifun,cpfmt1) qabs0(jb,jw0,ja0),qext0(jb,jw0,ja0),asym0(jb,jw0,ja0)
           enddo
        enddo
     enddo
     close (ifun)
     ifun=ifun+1

! Loop over the aerosol class
! ---------------------------
     do jka=1,nka
        do jkt=1,nkt
           ! Values of total radius (xw1) and volume mixing ratio of water (xa1) for each aerosol class
           ! ------------------------------------------------------------------------------------------
           xw1=rq(jkt,jka)
           xa1=1._dp-(rn(jka)/rq(jkt,jka))**3

           ! Locate xw1 in the tabulated xw0 grid, and calculate interpolation factor dx
           ! ---------------------------------------------------------------------------
           if (xw1 < xw0(1)) then
              ! case rq < xw0(1), warn user, set iw0 and dx values as if xw1=xw0(1)
              iw0 = 2
              dx  = 0._dp
              if (jaer==1) then ! Print this message only once
                 write(jpfunout,cpfmt2)'Warning: in SR intrad, rq < xw0(1), jka =',jka,'jkt =',jkt,'rq =', &
                      & xw1,'the value tabulated for xw0(1) =',xw0(1),'will be used.'
              end if
           else if (xw1 > xw0(nw0)) then
              ! case rq > xw0(nw0), warn user, set iw0 and dx values as if xw1=xw0(nw0)
              iw0 = nw0
              dx  = 1._dp
              if (jaer==1) then ! Print this message only once
                 write(jpfunout,cpfmt2)'Warning: in SR intrad, rq > xw0(nw0), jka =',jka,'jkt =',jkt,'rq =', &
                      & xw1,'the value tabulated for xw0(nw0) =',xw0(nw0),'will be used.'
              end if
           else
              ! general case: xw0(iw-1) < rq <= xw0(iw)
              iw0 = 2
              do while (xw1 > xw0(iw0))
                 iw0 = iw0 + 1
              end do
              dx=(xw1-xw0(iw0-1))/(xw0(iw0)-xw0(iw0-1))
           end if

           ! Locate xa1 in the tabulated xa0 grid, and calculate interpolation factor dy
           ! ---------------------------------------------------------------------------
           ! First test to detect potential model inconsistency
           !  (this case should never happend: if the particle grid has been properly defined,
           !   rq(jkt,jka) >= rn(jka) thus xa1 >= 0.)
           if (xa1 < xa0(1)) then
              write(jpfunerr,*)'Error in SR intrad: xa1 < xa0(1)=0.',jka,jkt,xa1
              call abortM(' Error 1 in SR intrad')
           else if (xa1 > xa0(na0)) then
              write(jpfunerr,*)'Error in SR intrad: xa1 > xa0(na0)=1.',jka,jkt,xa1
              call abortM(' Error 2 in SR intrad')
           else
              ! general case: xa0(ia0-1) < xa1 <= xa0(ia0)
              ia0 = 2
              do while (xa1 > xa0(ia0))
                 ia0 = ia0 + 1
              end do
              dy=(xa1-xa0(ia0-1))/(xa0(ia0)-xa0(ia0-1))
           end if


           ! Interpolation factors
           ! ---------------------
           dxdy   = dx*dy
           dxdy1  = dx*(1._dp-dy)
           dx1dy  = (1._dp-dx)*dy
           dx1dy1 = (1._dp-dx)*(1._dp-dy)

           ! Interpolated values of qabs, qext and asym
           ! ------------------------------------------
           do jb=1,mb
              qabs(jb,jkt,jka,jaer) = dxdy   * qabs0(jb,iw0,ia0)   &
                                    + dxdy1  * qabs0(jb,iw0,ia0-1) &
                                    + dx1dy  * qabs0(jb,iw0-1,ia0)   &
                                    + dx1dy1 * qabs0(jb,iw0-1,ia0-1)

              qext(jb,jkt,jka,jaer) = dxdy   * qext0(jb,iw0,ia0)   &
                                    + dxdy1  * qext0(jb,iw0,ia0-1) &
                                    + dx1dy  * qext0(jb,iw0-1,ia0)   &
                                    + dx1dy1 * qext0(jb,iw0-1,ia0-1)

              asym(jb,jkt,jka,jaer) = dxdy   * asym0(jb,iw0,ia0)   &
                                    + dxdy1  * asym0(jb,iw0,ia0-1) &
                                    + dx1dy  * asym0(jb,iw0-1,ia0)   &
                                    + dx1dy1 * asym0(jb,iw0-1,ia0-1)
           enddo

        end do ! jkt loop
     end do ! jka loop

  end do ! jaer loop

end subroutine intrad
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine ipdata
!
! Description :
! -----------
!    input data for radiation code:
!    reads from file 'pifm2_yymmdd.dat'
!    writes in several common blocks
!    descriptions of variables can be found at the end of file 'pifm2_yymmdd.dat'
!    (part of these description written in German)


! Modifications :
! -------------
  !    Jul-2016  Josue Bock  Header, comments
  !                          Use modules for parameters
  !
  !    Oct-2016  Josue Bock  More cleaning, after checking the variables actually used
  !                            (several old variables from PIFM1, apparently)
  !
  ! 28-Oct-2017  Josue Bock  Reindexed pibtab (thus changed input file) so that it is read
  !                            in optimum Fortran order (leftmost is innermost in loops)
  !                          Replaced dimension 12 by mbir
  !
  ! 15-Nov-2017  Josue Bock  Merged input data from pifm2*.dat and initr*.dat
  !                          Cleaned some unused data from PIFM1 (mostly in cb56)
  !                          Adopted the same data file structure as WS did in his Mistra-F90 version
  !                            ==> all comment lines are from his version
  !                          Each data block is now separated by a comment line.
  !                          Also added a former common block /sol/ to remove hardcoded values from
  !                            the radiative code.
  !                          Final cleaning after Fortran90 conversion
  !
  ! 22-Sep-2019  Josue Bock  Remove format statements (coding standards)

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! Imported Parameters:
       cinpdir,      &
! Imported Routines:
       abortM

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfundatarad

  USE global_params, ONLY : &
! Imported Parameters:
       mb,                  &
       mbs,                 &
       mbir,                &
       ncw,                 &
       nrlay

  USE precision, ONLY : &
       dp

  implicit none

! Local parameters:
  character (len=*), parameter :: cpfmt1='(a)'       ! format for separation lines
  character (len=*), parameter :: cpfmt2='(8e16.8)'  ! format for tabulated data

! Local scalars:
  character (len=len_trim(cinpdir)+16) :: clfname
  integer :: istat

! Common blocks:
  common /aeag/ feux(8),seanew(8,mb,4),saanew(8,mb,4),ganew(8,mb,4)
  real (kind=dp) :: feux, seanew, saanew, ganew
  common /band1/ hk1(10),fk1o3(10)
  real (kind=dp) :: hk1, fk1o3
  common /band2/ hk2(8),c2h2o(3,11,8)
  real (kind=dp) :: hk2, c2h2o
  common /band3/ hk3(12),c3h2o(3,11,12)
  real (kind=dp) :: hk3, c3h2o
  common /band4/ hk4(7),c4h2o(3,11,7)
  real (kind=dp) :: hk4, c4h2o
  common /band5/ hk5(12),c5h2o(3,11,12)
  real (kind=dp) :: hk5, c5h2o
  common /band6/ hk6(5),c6h2o(3,11,5)
  real (kind=dp) :: hk6, c6h2o
  common /band7/ hk7(2),c7h2o(3,19,2)
  real (kind=dp) :: hk7, c7h2o
  common /band8/ hk8(3),c8h2o(3,19,3)
  real (kind=dp) :: hk8, c8h2o
  common /band9/ hk9(4),c9h2o(3,19,4)
  real (kind=dp) :: hk9, c9h2o
  common /band10/ hk10(4),c10h2o(3,19,4),c10ch4(3,19),c10n2o(3,19)
  real (kind=dp) :: hk10, c10h2o, c10ch4, c10n2o
  common /band11/ hk11(3),c11h2o(3,19,3),c11ch4(3,19),c11n2o(3,19)
  real (kind=dp) :: hk11, c11h2o, c11ch4, c11n2o
  common /band12/ hk12(5),c12o3(3,19,5),c12h2o(3,19)
  real (kind=dp) :: hk12, c12o3, c12h2o
  common /band13/ hk13(2),c13h2o(3,19,2)
  real (kind=dp) :: hk13, c13h2o
  common /band14/ hk14(10),c14hca(3,19,10),c14hcb(3,19,10)
  real (kind=dp) :: hk14, c14hca, c14hcb
  common /band15/ hk15(12),c15hca(3,19,12),c15hcb(3,19,12)
  real (kind=dp) :: hk15, c15hca, c15hcb
  common /band16/ hk16(7),c16h2o(3,19,7)
  real (kind=dp) :: hk16, c16h2o
  common /band17/ hk17(7),c17h2o(3,19,7)
  real (kind=dp) :: hk17, c17h2o
  common /band18/ hk18(8),c18h2o(3,19,8)
  real (kind=dp) :: hk18, c18h2o

  common /cb19/ berayl(mbs),bea(mb,nrlay),baa(mb,nrlay),ga(mb,nrlay)
  real (kind=dp) :: berayl, bea, baa, ga

  common /cb56/ o3un(52)
  real (kind=dp) :: o3un

  common /sol/ s0b(mbs), s0tot
  real (kind=dp) :: s0b, s0tot

  common /plancd/ ttab(35),pibtab(35,mbir) ! jjb was used in FN fst4 (called by SR plancktab), both unused now
  real (kind=dp) :: ttab, pibtab           !     left for backward compatibility, might be deleted in the future

  common /was1/ ret(ncw),r2wt(ncw),b2wt(ncw,mb),w2wt(ncw,mb),g2wt(ncw,mb)
  real (kind=dp) :: ret, r2wt, b2wt, w2wt, g2wt

! == End of declarations =======================================================

  clfname=TRIM(cinpdir)//'pifm2_171115.dat'

  open(unit=jpfundatarad, file=clfname, status='old', iostat=istat)
  if (istat /= 0) call abortM ('Error in SR ipdata: cannot open radiation data file: '//clfname)

! planck table                 ! could be deleted if fst4 and plancktab are deleted
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) ttab
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) pibtab
! effective droplet radii and radiation coefficients for water
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) ret
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) r2wt
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) b2wt
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) w2wt
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) g2wt
! table of reference relative humidities and tabulated radiation coefficients for aerosols
! last index of seanew, saanew, ganew = aerosol type: 1=rural, 2=urban, 3=ocean, 4=background
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) feux
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) seanew
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) saanew
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) ganew
! solar energy in the 6 SW spectral bands
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) s0b
! radiation coefficients for gas absorption for every spectral band
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk1
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) fk1o3
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk2
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c2h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk3
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c3h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk4
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c4h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk5
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c5h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk6
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c6h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk7
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c7h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk8
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c8h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk9
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c9h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk10
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c10h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c10ch4
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c10n2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk11
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c11h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c11ch4
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c11n2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk12
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c12o3
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c12h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk13
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c13h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk14
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c14hca
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c14hcb
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk15
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c15hca
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c15hcb
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk16
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c16h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk17
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c17h2o
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) hk18
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) c18h2o
! unreduced ozone amount from Craig table
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) o3un
! coefficients for rayleigh scattering
  read(jpfundatarad,cpfmt1)
  read(jpfundatarad,cpfmt2) berayl

  close(jpfundatarad)

end subroutine ipdata
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine initr
!
! Description :
! -----------
!    standard atmosphere between etw(n) and 50000 meters.
!    input of constant radiation parameters


! Modifications :
! -------------
  !     ?        Roland      Calculation of qmo3 for layers 1-80 (bug fix)
  !              von Glasow
  !
  !    Jul-2016  Josue Bock  Removal of labeled do loops.
  !                          Use modules for parameters
  !                          BUGFIX: data file updated with berayl(6)
  !                                 (previously, an old version berayl(4) was read here, then
  !                                 berayl was redefined with 'data' in strahl). But for some
  !                                 unknown reason, only the last two values were updated, while
  !                                 the first four outdated values were NOT overwritten)
  !
  !    Oct-2016  Josue Bock  corrected several array size
  !                          initialisation of qmo3(nrlev) was missing
  !                          initialisation of rnaer(nrlay) was missing
  !                          removed rnaer from /cb02/, turned into a local array
  !
  ! 10-Mar-2017  Josue Bock  corrected rf=0.08 for level 30000 (was 0.8)
  !
  ! 14-Nov-2017  Josue Bock  removed ntypdx(nrlay) (droplet type), initialised =4 but unused
  !                          removed ntypax(nrlay) (aerosol type) from cb01, local declaration
  !
  ! 17-Nov-2017  Josue Bock  Recent changes: cb02 is now cb01
  !                                          added CB /sol/, and s0tot calculation in this SR
  !                                          rewritten interpolations
  !                                          reorganised this SR, uniformisation of some variable names
  !                          BUGFIX: the qmo3x calculation was wrong (it rotated the indexes twice, and
  !                                  used the pressure array in the wrong line
  ! 18-Nov-2017              BUGFIX: missing initialisation of beax, baax and gax for layers 1:n-1
  !
  ! 18-Feb-2019  Josue Bock  Stop the program using abortM

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! Imported Routines:
       abortM,       &
! Imported Parameters:
       jpAlbedoOpt

  USE constants, ONLY : &
! Imported Parameters:
       g, &             ! Gravitational acceleration (m/s**2)
       r0               ! Specific gas constant of dry air, in J/(kg.K)

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunerr

  USE global_params, ONLY : &
! Imported Parameters:
       n,                   &
       nrlay,               &
       nrlev,               &
       mb,                  &
       mbs,                 &
       mbir

  USE precision, ONLY : &
! Imported parameters:
       dp

  implicit none

! Local scalars:
  integer :: jb, jz, k, na
  integer :: io3_inf, io3_sup                  ! indexes for O3 interpolation
  integer :: irh_inf, irh_sup                  ! indexes for relative humidity interpolation
  real (kind=dp) :: dz                         ! Height increment of extra layers up to the tropopause [m]
  real (kind=dp) :: dzo3, zo3_inf, zo3_sup     ! heights for O3 interpolation [m]
  real (kind=dp) :: gamma                      ! Vertical temperature gradients [K/m]
  real (kind=dp) :: rf                         ! relative humidity
  real (kind=dp) :: drh, omdrh, xdrh, xomdrh   ! coefficients for aerosol coefficients interpolation
  real (kind=dp) :: xnaer                      ! total number concentration of aerosol particles [m**-3]

! Local arrays:
  integer :: ntypax(nrlay)
  real (kind=dp) :: eta_o3(nrlev), u_o3(nrlay) ! intermediate variables for O3 interpolation
  real (kind=dp) :: rnaer(nrlay)               ! total number concentration of aerosol particles [cm**-3]

! External function:
  real (kind=dp), external :: p21              ! saturation water vapour pressure [Pa]

! Common blocks:
  common /aeag/ feux(8),seanew(8,mb,4),saanew(8,mb,4),ganew(8,mb,4)
  real (kind=dp) :: feux, seanew, saanew, ganew

  common /cb01/ tx(nrlev),px(nrlev),rhox(nrlev),xm1x(nrlev),      &
                rho2wx(nrlay),fracx(nrlay),zx(nrlev),thkx(nrlay), &
                beax(mb,nrlay),baax(mb,nrlay),gax(mb,nrlay),      &
                qmo3x(nrlev),rewx(nrlay), tsx
  real (kind=dp) :: tx,px,rhox,xm1x,rho2wx,fracx,zx,thkx, &
                    beax, baax, gax, qmo3x, rewx, tsx

  common /cb16/ u0,albedo(mbs),thk(nrlay)
  real (kind=dp) :: u0, albedo, thk

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb56/ o3un(52)
  real (kind=dp) :: o3un

  common /sol/ s0b(mbs), s0tot
  real (kind=dp) :: s0b, s0tot

  common /tmp2/ as(mbs),ee(mbir)
  real (kind=dp) :: as, ee

! == End of declarations =======================================================

! albedo of the ground for the six solar wavelength regions of the radiation code
  select case (jpAlbedoOpt)
  case (0)
     albedo(:)=0.05_dp
  case (1)
     albedo(:)=0.8_dp       ! snow
  case default
     write (jpfunerr,*) "Wrong choice in namelist for jpAlbedoOpt, must be 0 or 1"
     call abortM ("Error in SR initr")
  end select

  as = albedo
! emissivity of the ground
  ee(:)=1.0_dp

! total solar energy = sum of solar energies for each spectral bands
  s0tot = sum(s0b(:))


! radiation model level
!----------------------
  ! Use meteorological grid for the lowest n-1 layers (note that etw(1) = 0.)
  ! Then add 7 equidistant layers up to 11 km (tropopause height in standard atmosphere)
  ! Then add 4 layers at 20, 30, 40, and 50 km height, and the last one at 100 km height
  ! Note that the same vertical structure is used in the photolysis code (except the uppermost level,
  ! whose temperature and pressure are calculated differently in jrate.f (see SR read_data)).

  ! First check consistency: nlev = n + 11 in the current settings
  if(nrlev /= n+11) then
     write(jpfunerr,*) 'Error: nlev /= n+11, change nlev in global_params,'
     write(jpfunerr,*) '  or change the extra layers settings in SR initr'
     call abortM(' Error 1 in SR initr')
  end if

  ! Main model grid, define radiation level <=> "wall" values
  zx(1:n-1) = etw(1:n-1)

  ! Even if etw(n-1) should be much lower than 11 km, check anyway
  if(zx(n-1) >= 11000._dp) then
     write(jpfunerr,*) 'Error: etw(n-1) >= 11 km, check the vertical grid'
     call abortM(' Error 2 in SR initr')
  end if
  ! dz: height increment up to the tropopause [m]
  dz = (11000._dp - zx(n-1)) / 7._dp

  do k=n,n+6
     zx(k)=zx(k-1)+dz
  end do
  zx(n+7)  =  20000._dp
  zx(n+8)  =  30000._dp
  zx(n+9)  =  40000._dp
  zx(n+10) =  50000._dp
  zx(n+11) = 100000._dp

  ! Layer thicknesses calculated during initialisation
  do k=1,nrlay
     thkx(k) = zx(k+1) - zx(k)
  end do


! Initialise aerosol coefficients baax, beax and gax to zero
!   (they will be overwritten in SR load1 from index 1 to n, only if mic=true)
! -------------------------------------------
  baax(:,:n-1) = 0._dp
  beax(:,:n-1) = 0._dp
  gax (:,:n-1) = 0._dp

! Initialise liquid water variables: fracx, rho2wx, and rewx: no clouds
!   (they will be overwritten in SR load1 from index 1 to nf, only if mic=true)
! -------------------------------------------
  fracx(:)  = 0._dp
  rewx(:)   = 0._dp
  rho2wx(:) = 0._dp


! load actual meteorological profiles in the model domain z(1)-z(n-1)
! -------------------------------------------------------------------
  call load1

! define meteorological profile from n-1 to nrlay / nrlev
! -------------------------------------------------------

  ! temperature gradient and relative humidity for the next layers:
  gamma=0.0065_dp
  rf=.3_dp
  do k=n,n+6
     tx(k)=tx(k-1)-gamma*thkx(k-1)
     px(k)=px(k-1)*(tx(k)/tx(k-1))**(g/(r0*gamma))
     xm1x(k)=0.62198_dp*rf/(px(k)/p21(tx(k))-0.37802_dp*rf)
     rhox(k)=px(k)/(r0*tx(k)*(1._dp+.608_dp*xm1x(k)))
     rnaer(k)=100._dp
  enddo

  k=n+7
  ! relative humidity for the next layer:
  rf=.02_dp
  tx(k)=tx(k-1)
  px(k)=px(k-1)*dexp(-g*(zx(k)-zx(k-1))/(r0*tx(k)))
  xm1x(k)=0.62198_dp*rf/(px(k)/p21(tx(k))-0.37802_dp*rf)
  rhox(k)=px(k)/(r0*tx(k)*(1._dp+.608_dp*xm1x(k)))
  rnaer(k)=0._dp

  k=n+8
  ! temperature gradient and relative humidity for the next layer:
  gamma=-0.001_dp
  rf=.005_dp
  tx(k)=tx(k-1)-gamma*thkx(k-1)
  px(k)=px(k-1)*(tx(k)/tx(k-1))**(g/(r0*gamma))
  xm1x(k)=0.62198_dp*rf/(px(k)/p21(tx(k))-0.37802_dp*rf)
  rhox(k)=px(k)/(r0*tx(k)*(1._dp+.608_dp*xm1x(k)))
  rnaer(k)=0._dp

  k=n+9
  ! temperature gradient and relative humidity for the next layer:
  gamma=-0.0026_dp
  rf=.00005_dp
  tx(k)=tx(k-1)-gamma*thkx(k-1)
  px(k)=px(k-1)*(tx(k)/tx(k-1))**(g/(r0*gamma))
  xm1x(k)=0.62198_dp*rf/(px(k)/p21(tx(k))-0.37802_dp*rf)
  rhox(k)=px(k)/(r0*tx(k)*(1._dp+.608_dp*xm1x(k)))
  rnaer(k)=0._dp

  k=n+10
  ! temperature gradient and relative humidity for the next layer:
  gamma=-0.0018_dp
  rf=.000002_dp
  tx(k)=tx(k-1)-gamma*thkx(k-1)
  px(k)=px(k-1)*(tx(k)/tx(k-1))**(g/(r0*gamma))
  xm1x(k)=0.62198*rf/(px(k)/p21(tx(k))-0.37802_dp*rf)
  rhox(k)=px(k)/(r0*tx(k)*(1._dp+.608_dp*xm1x(k)))
  rnaer(k)=0._dp

! fictitious level at infinity
  tx(nrlev)=210._dp
  px(nrlev)=0._dp
  xm1x(nrlev)=0._dp
  rhox(nrlev)=0._dp

! aerosol type: 1 rural 2 urban 3 maritme 4 tropospheric
  ntypax(:)=4

! ozone concentration profile
!----------------------------
! interpolate unreduced ozone amounts from craig table
! save preliminarily the interpolated total amount in array etao3
  ! Notes about interpolation:
  !   - while the tabulated o3un has 52 values, the latest one (=0) is not used here
  !   - tabulated values from the surface (zo3_inf=0) to 51 km height, every km
  !   - zx is a monotonic function, thus io3_inf is not reset to 1 in the jz do-loop.
  !
  ! In the radiative code, qmo3 is used twice, after being  multiplied by 2.3808*(delta_P);
  !   this is thus identical to u_o3 (ozone path length).
  io3_inf = 1
  zo3_sup = real(io3_inf,dp) * 1000._dp
  do jz=1,nrlev
     do while (zx(jz) > zo3_sup .and. io3_inf < 51)
        io3_inf = io3_inf + 1
        zo3_sup = real(io3_inf,dp) * 1000._dp
     end do
     if (io3_inf < 51) then
        io3_sup = io3_inf + 1
        zo3_inf = zo3_sup - 1000._dp
        dzo3 = (zx(jz)-zo3_inf) / 1000._dp
        eta_o3(jz) = o3un(io3_inf) + (o3un(io3_sup)-o3un(io3_inf)) * dzo3
     else
        ! this case should happend only for the latest level at 100 km
        eta_o3(jz) = 0._dp
     end if
  end do

  do jz=1,nrlay
     ! path lengths for every layer
     u_o3(jz)=(eta_o3(jz)-eta_o3(jz+1))*0.01_dp
     qmo3x(jz)=u_o3(jz)/(2.3808_dp*(px(jz)-px(jz+1)))
  end do
  qmo3x(nrlev)=0._dp


! diagnostic aerosol properties above n
!--------------------------------------
  do jz=n,nrlay
     na=ntypax(jz)
     if (na.gt.0.and.rnaer(jz).gt.0._dp) then
        rf=xm1x(jz)*px(jz)/(p21(tx(jz))*(.62198_dp+.37802_dp*xm1x(jz)))

        irh_sup = 2
        do while (rf > feux(irh_sup) .and. irh_sup < 8)
           irh_sup = irh_sup + 1
        end do
        irh_inf = irh_sup - 1

        drh = (rf-feux(irh_inf)) / (feux(irh_sup)-feux(irh_inf))
        omdrh = 1 - drh

        xnaer = rnaer(jz) * 1.e6_dp
        xdrh   = xnaer * drh
        xomdrh = xnaer * omdrh
        do jb=1,mb
           beax(jb,jz) = xomdrh * seanew(irh_inf,jb,na) + xdrh * seanew(irh_sup,jb,na)
           baax(jb,jz) = xomdrh * saanew(irh_inf,jb,na) + xdrh * saanew(irh_sup,jb,na)
           gax (jb,jz) = omdrh  * ganew (irh_inf,jb,na) +  drh * ganew (irh_sup,jb,na)
        enddo
     else
        beax(:,jz) = 0._dp
        baax(:,jz) = 0._dp
        gax (:,jz) = 0._dp
     endif
  enddo

end subroutine initr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine load1
!
! Description :
! -----------
!    loading actual variables and actual optical parameters


! Method :
! ------
  ! tx, px, xm1x and rhox are defined for half-level (= layer) in Mistra.
  ! tx, px and xm1x are linearly interpolated to get level values for the radiative code
  ! rhox is then recalculated using the interpolated values.
  !
  ! At the ground level, px(1)=p(1), but tx(1)=t(2) and xm1x(1)=xm1(2).
  ! Indeed, while the pressure varies continuously, the surface temperature (t(1)) can not be used as
  ! the air temperature at the surface. Thus, it is expected that using t(2) is less wrong than t(1).


! Author :
! ------
  !    Andreas Bott and others


! Modifications :
! -------------
  !    Jul-2016  Josue Bock  Removal of labeled do loops.
  !                          Header
  !                          Use module for parameters
  !                          All explicit declarations and implicit none
  !
  !    Oct-2016  Josue Bock  Removed /nox/ tauno2(nrlay) which was set to 0
  !                            (see SR tau in nrad.f90 for explanations)
  !
  ! 13-Oct-2016  Josue Bock  Added CB /neu/ icld and declared icld=0. This was missing for two routines
  !
  ! 22-Oct-2017  Josue Bock  icld further removed along with an important bugfix in SR water
  !
  ! 18-Nov-2017  Josue Bock  introduced a test in the final calculation of gax (if beax-baax > 0)
  !                            previously, this case was handled by adding a small value (gax = gax / (beax-baax+1e-15),
  !                            which might lead to unexpectedly high values (just a guess). It seems safer to do a proper test.
  !                          BUGFIX: beax, baax and gax must be initialised to 0. for all spectral bands. If not, at night, the
  !                                  last solar values calculated during the day were still used.

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! Imported Parameters:
       mic

  USE constants, ONLY : &
! Imported Parameters:
       pi,              &
       r0                 ! Specific gas constant of dry air, in J/(kg.K)

  USE global_params, ONLY : &
! Imported Parameters:
       jptaerrad,           &
       n,                   &
       nf,                  &
       nrlay,               &
       nrlev,               &
       nka,                 &
       nkt,                 &
       mb,                  &
       mbs

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  logical :: llswitch
  integer :: ib0                                  ! first spectral band index (day => 1, night => 7)
  integer :: jb                                   ! loop index for spectral band number
  integer :: jka,jkt                              ! loop indexes, 2D microphysical grid
  integer :: jz, jzm                              ! loop indexes, vertical
  integer :: ka, ka0
  real (kind=dp) :: x0                            ! interpolation factor
  real (kind=dp) :: zeit, horang, rlat, rdec, u00, ru0
  real (kind=dp) :: znum,zdenom,zfix   ! calculation of effective drop radius

! Common blocks:
  common /cb01/ tx(nrlev),px(nrlev),rhox(nrlev),xm1x(nrlev),      &
                rho2wx(nrlay),fracx(nrlay),zx(nrlev),thkx(nrlay), &
                beax(mb,nrlay),baax(mb,nrlay),gax(mb,nrlay),      &
                qmo3x(nrlev),rewx(nrlay), tsx
  real (kind=dp) :: tx,px,rhox,xm1x,rho2wx,fracx,zx,thkx, &
                    beax, baax, gax, qmo3x, rewx, tsx

  common /cb08/ re1(nkt), re2(nkt), re3(nkt)
  real (kind=dp) :: re1, re2, re3

  common /cb16/ u0,albedo(mbs),thk(nrlay)
  real (kind=dp) :: u0, albedo, thk

  common /cb18/ alat,declin                ! for the SZA calculation
  real (kind=dp) :: alat,declin

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb49/ qabs(mb,nkt,nka,jptaerrad),qext(mb,nkt,nka,jptaerrad),asym(mb,nkt,nka,jptaerrad)
  real (kind=dp) :: qabs,qext,asym

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff,fsum
  integer :: nar

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real (kind=dp) :: theta, thetl, t, talt, p, rho

  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real (kind=dp) :: xm1, xm2, feu, dfddt, xm1a

! == End of declarations =======================================================

  tsx     = t(1)
  tx(1)   = t(2)
  px(1)   = p(1)
  xm1x(1) = xm1(2)
  rhox(1) = px(1)/(r0*tx(1)*(1._dp+.608_dp*xm1x(1)))
  do jz = 2,n-1
     x0      = 0.5_dp*detw(jz)/deta(jz)
     tx(jz)   = t(jz)+(t(jz+1)-t(jz))*x0
     px(jz)   = p(jz)+(p(jz+1)-p(jz))*x0
     xm1x(jz) = xm1(jz)+(xm1(jz+1)-xm1(jz))*x0
     rhox(jz) = px(jz)/(r0*tx(jz)*(1._dp+.608_dp*xm1x(jz)))
  enddo

! calculate u0 from geogr. latitude, declination and hourangle
! make correction because of spherical surface of the earth
  zeit=lst*3600._dp+lmin*60._dp
  horang=7.272205e-05_dp*zeit-pi
! pi/180=1.745329e-02
  rlat=alat*1.745329e-02_dp
  rdec=declin*1.745329e-02_dp
  u00=cos(rdec)*cos(rlat)*cos(horang)+sin(rdec)*sin(rlat)
  ru0=6371._dp*u00
  u0=8._dp/(dsqrt(ru0**2+102000._dp)-ru0)

  if (mic) then

     ! Wavelenght bands: 1-6 = shortwave (sun), 7-18 = longwave (earth)
     ! 24h per day: earth longwave radiation taken into account
     ! Depending on solar zenith angle: sun shortwave radiation also accounted for
     if (u0.gt.1.e-02_dp) then
        ib0 = 1
     else
        ib0 = 7
     end if

! =====================================================================================================
! -- 1 -- aerosol coefficients: baax, beax, gax
! =====================================================================================================
! the quantities needed in the radiation code are given by
! the sum over jka and jkt of: qabs(jka,jkt,l)*pi*r**2*ff(jkt,jka,jz)
! (same formula for qext) asym has weighting factor qsca*f(jz)

     ! Initialisation for all spectral bands
     ! (it is mandatory to set all these parameters to zero, even if lu0=7. If not, the previous, non-zero
     !  values, would still exist).
     baax(:,1:n-1) = 0._dp
     beax(:,1:n-1) = 0._dp
     gax (:,1:n-1) = 0._dp

     do jz=1,n-1
        jzm = jz+1
        ! index jz is related to radiative code variables
        ! index jzm ("m" as "mistra") is shifted by +1 to account for the first, infinitesimally thin layer, ignored here

        ka0=nar(jzm)

        do jka=1,nka
           ka=ka0
           if (rn(jka).lt.0.5_dp.and.ka0.eq.3) ka=2 ! ka=1 urban, ka=2 rural, ka=3 ocean
           do jkt=1,nkt
              x0=pi*1.e-6_dp*rq(jkt,jka)**2*ff(jkt,jka,jzm)
              do jb = ib0,mb
                 baax(jb,jz) = baax(jb,jz) + qabs(jb,jkt,jka,ka)*x0
                 beax(jb,jz) = beax(jb,jz) + qext(jb,jkt,jka,ka)*x0
                 gax (jb,jz) = gax (jb,jz) + asym(jb,jkt,jka,ka)*x0*(qext(jb,jkt,jka,ka)-qabs(jb,jkt,jka,ka))
              end do
           end do
        end do

        do jb=ib0,mb
           if (beax(jb,jz)-baax(jb, jz) > 0._dp) then
              gax(jb,jz) = gax(jb,jz) / (beax(jb,jz)-baax(jb,jz))
           else
              gax(jb,jz) = 0._dp
           end if
        end do

     end do ! jz loop over vertical layers


! =====================================================================================================
! -- 2 -- liquid water related variables: rho2wx, fracx, and rewx
! =====================================================================================================
! note that this loop goes up to nf-1 only, not n-1
! indeed, the cloud top is Mistra cannot be higher than nf
!
! Also note that the 2D microphysical grid is read in an unusual order here (jkt first, then jka)
! This is because a given jkt class has the same amount of water (expressed as an equivalent radius,
! surface, or volume in variables re1, re2 and re3) for all jka classes.

     llswitch = .false. ! the next section has to be switched off for mistra

     if (llswitch) then
     do jz=1,nf-1
        jzm = jz + 1

! rho2wx is used in SR water, both as a switch (if rho2 < 1.e-5 : no absorption in liquid water),
!   and in the calculations if they are done. NB: frac could be used as well as a switch, since it
!   is set = 1 between lcl and lct, and lcl and lct are defined according to the LWC.
!   For the same reason, it is not necessary to define it above nf-1.
! It is a layer (=half level) value, thus no interpolation is required, just an index shift to account
!   for Mistra layer 1 which is the infinitesimally thin layer, while the first "true" layer is #2.
        rho2wx(jz) = xm2(jzm)

! define frac, for use in SR frr and langw,
!     and rew, for use in SR water
        if (jzm.ge.lcl .and. jzm.le.lct) then
           fracx(jz) = 1._dp

           znum   = 0._dp
           zdenom = 0._dp
           do jkt = 1, nkt
              zfix = 0._dp
              do jka = 1, nka
                 zfix = zfix + ff(jkt,jka,jzm)
              end do
              znum   = znum   + re3(jkt)*zfix
              zdenom = zdenom + re2(jkt)*zfix
           end do
           rewx(jz) = znum/zdenom

        else
           fracx(jz) = 0._dp
           rewx(jz)  = 0._dp
        end if

     end do

     else ! llswitch
        rho2wx(1:nf-1) = 0._dp
        fracx(1:nf-1) = 0._dp
        rewx(1:nf-1)  = 0._dp
     end if ! llswitch

  end if ! mic = .true.

end subroutine load1
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine rotate_in (ldinit)

! Description :
! -----------
!    Rotate arrays filled from main program, that are indexed from ground to top (bottom-up: _bu),
!    towards radiative code input arrays, indexed from top to ground (top-down: _td)


! Interface :
! ---------
!    SR rotate_in is called by SR initstr during initialisation, and by SR str during the run

! Input :
! -----
!    - /cb01/ contains all "_bu" arrays that need to be rotated to be used in the radiative code

! Output :
! ------
!    - /cb02/
!    - /cb09/    liquid water variables
!    - /cb16/    thk(nrlay)
!    - /cb19/    aerosol variables
!    - /height/  zx (nrlev)

! Externals :
! ---------
!    none


! Method :
! ------
!    All extra layer variables (T, P, ...) are constant.
!    Thus, the arrays are rotated over all layers/levels during initialisation only.
!    Else, only the meteorological layers will be rotated.



! Author :
! ------
!    Josue Bock


! Modifications :
! -------------
  !    02-Apr-2017  Josue Bock   First version of this routine
  !    13-Nov-2017  Josue Bock   removed ntypa and ntypd from /cb02/, unused
  !    14-Nov-2017  Josue Bock   further removed ntypd from /cb01/
  !    19-Nov-2017  Josue Bock   reorganised CBs: cb09 for water related variables
  !                              efficiency: qmo3 rotated only during initialisation, and aerosol/water related
  !                                          variables rotated only if mic (and during initialisation)

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! Imported Parameters:
       mic

  USE global_params, ONLY : &
! Imported Parameters:
       n,                   &
       nrlay,               &
       nrlev,               &
       mb, mbs

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  logical, intent(in) :: ldinit          ! called during initialisation, or not.

! Local scalars:
  integer :: jlbu, jltd                  ! running indexes (bu = bottom-up, td = top-down)
  integer :: nmin

! Common blocks:
  common /cb01/ tx(nrlev),px(nrlev),rhox(nrlev),xm1x(nrlev),         &
                rho2wx(nrlay),fracx(nrlay),zx_bu(nrlev),thkx(nrlay), &
                beax(mb,nrlay),baax(mb,nrlay),gax(mb,nrlay),         &
                qmo3x(nrlev),rewx(nrlay), tsx
  real (kind=dp) :: tx,px,rhox,xm1x,rho2wx,fracx,zx_bu,thkx, &
                    beax, baax, gax, qmo3x, rewx, tsx

  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),ts
  real (kind=dp) :: t, p, rho, xm1, ts

  common /cb09/ frac(nrlay),rew(nrlay),rho2w(nrlay)
  real (kind=dp) :: frac, rew, rho2w

  common /cb16/ u0,albedo(mbs),thk(nrlay)
  real (kind=dp) :: u0, albedo, thk

  common /cb19/ berayl(6),bea(mb,nrlay),baa(mb,nrlay),ga(mb,nrlay)
  real (kind=dp) :: berayl, bea, baa, ga

  common /height/ zx(nrlev)        ! Level height [m] indexed from ground to top
  real (kind=dp) :: zx

  common /ozon/ qmo3(nrlev)
  real (kind=dp) :: qmo3

! == End of declarations =======================================================


! Vertical grid arrays and ozone profile: rotate only during initialisation
  if(ldinit) then
     ! thk(nrlay)
     do jltd = 1,nrlay
        jlbu = nrlay - jltd + 1

        thk(jltd) = thkx(jlbu)
     end do
     ! zx(nrlev)
     ! qmo3x(nrlev)
     do jltd = 1,nrlev
        jlbu = nrlev - jltd + 1

        zx(jltd) = zx_bu(jlbu)
        qmo3(jltd) = qmo3x(jlbu)
     end do
  end if


! Define max index that will be rotated
  if(ldinit) then
     nmin = 1
  else
     nmin = nrlev - n + 1 ! fits nrlev size arrays. For nrlev size arrays, use nmin-1 instead
  end if

! nrlev size variables, rotated any time
  do jltd = nmin, nrlev
     jlbu = nrlev - jltd + 1

     t(jltd) = tx(jlbu)
     p(jltd) = px(jlbu)
     rho(jltd) = rhox(jlbu)
     xm1(jltd) = xm1x(jlbu)
     !(jltd) = x(jlbu)
  end do

! nrlay size variables, which are all related to aerosols or water.
! Then need to be rotated only if mic=true, and during initialisation
  if (ldinit .or. mic) then
     do jltd = MAX(nmin-1,1), nrlay
        jlbu = nrlay - jltd + 1

        rho2w(jltd) = rho2wx(jlbu)
        rew(jltd) = rewx(jlbu)
        frac(jltd) = fracx(jlbu)
        bea(:,jltd) = beax(:,jlbu)
        baa(:,jltd) = baax(:,jlbu)
        ga(:,jltd) = gax(:,jlbu)

     end do
  end if

  ts = tsx

end subroutine rotate_in
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine rotate_out

! Description :
! -----------
  !    Rotate arrays outputed by the radiative code, that are indexed from top to the ground
  !    (top-down: _td) to the corresponding variables used in the main program, that are indexed
  !    from the ground to the top (bottom-up: _bu).


! Interface :
! ---------
  !    SR rotate_out is called after the radiative code

! Input :
! -----
  !    - /cb10/ totrad(mb,nrlay)
  !    - /cb15/ fnseb,flgeg,hr(nrlay)

! Output :
! ------
  !    - /cb11/ totrad(mb,n)
  !    - /cb48/ sk,sl,dtrad(n)

! Externals :
! ---------
  !    none


! Method :
! ------
  !    Mistra layer #1 is surface
  !    Thus, fill dtrad and totrad_bu from indexes #2 to #n
  !    Index #2 in Mistra must match index #nrlay in the radiative code


! Author :
! ------
  !    Josue Bock


! Modifications :
! -------------
  !    13-Nov-2017  Josue Bock   First version of this routine

! == End of header =============================================================


! Declarations:
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       n,                   &
       nrlay,               &
       mb

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  integer :: jzbu, jztd                 ! running indexes (bu = bottom-up, td = top-down)

! Common blocks:
  common /cb10/ totrad_td (mb,nrlay)
  real (kind=dp) :: totrad_td

  common /cb11/ totrad_bu (mb,n)
  real (kind=dp) :: totrad_bu

  common /cb15/ fnseb,flgeg,hr(nrlay)
  real (kind=dp) :: fnseb, flgeg, hr

  common /cb48/ sk,sl,dtrad(n),dtcon(n)
  real (kind=dp) :: sk, sl, dtrad, dtcon

! == End of declarations =======================================================


  sk = fnseb
  sl = flgeg

  dtrad(1)       = 0._dp
  totrad_bu(:,1) = 0._dp

  do jzbu = 2,n
     jztd = nrlay-jzbu+2

     dtrad(jzbu)=hr(jztd)
     totrad_bu(:,jzbu)=totrad_td(:,jztd)
  end do

end subroutine rotate_out
