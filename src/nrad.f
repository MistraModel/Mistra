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
! Description:
!    Main program of the radiation module:
!    Calls subroutines for calculation of radiative transfer equation for 18
!    spectral bands (6 solar and 12 ir) and associated cumulative probabilities.
!    Finally integrates results over the spectral intervals
!
!    Changes were done in the old part of the program to upgrade it to 18 spectral
!    bands. This includes reading in more information in the routine 'initr'
!    with qext, qabs, and asymm. The upgrade then follows on to other variables
!    which are calculated from these three, namely bea baa and ga. There are also
!    changes in the arrays which store the radiation data as a function of wavelength.
!
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
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2.2       10/2016  Major work, cleaning, debuging                 <Josue Bock>
!                    Consistency of array size throughout the code
!
!           07/2016  Header, comments and cleaning                  <Josue Bock>
!
! 2.1       ?        Original code, PIFM2                           <Andreas Bott>
!
! 1.0       ?        Original code, PIFM1                           <Andreas Bott>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     mb,
     &     mbs,
     &     mbir,
     &     nrlay,
     &     nrlev

      implicit none

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

      common /cb11/ totrad (mb,nrlay)
      double precision totrad

      common /cb15/ fnseb,flgeg,hr(nrlay)
      double precision fnseb, flgeg, hr

      common /cb16/ u0,albedo(mbs),thk(nrlay)
      double precision u0, albedo, thk

      common /cb19/ berayl(6),bea(mb,nrlay),baa(mb,nrlay),ga(mb,nrlay)
      double precision berayl, bea, baa, ga

      common /extra/ waer(nrlay),taer(nrlay),plaer(2,nrlay)
      double precision waer, taer, plaer

      common /extra_2/ taer_s(nrlay),taer_a(nrlay),ga_pl(nrlay) ! for photolysis calculation
      double precision taer_s, taer_a, ga_pl

      common /kurz/ fs1(nrlev),fs2(nrlev),totds(nrlev),ss(nrlev),
     &              fsn(nrlev),dtdts(nrlay)
      double precision fs1, fs2, totds, ss, fsn, dtdts          ! integrated radiation fluxes

      common /lang/ fl1(nrlev),fl2(nrlev),fln(nrlev),dtdtl(nrlay)           ! integrated radiation fluxes
      double precision fl1, fl2, fln, dtdtl

      common /leck2/ sf(nrlev),sw(nrlev),ssf(nrlev),ssw(nrlev),             ! radiation fluxes
     &     f2f(nrlev),f2w(nrlev),f1f(nrlev),f1w(nrlev)
      double precision sf, sw, ssf, ssw, f2f, f2w, f1f, f1w

      common /neu/ icld
      integer icld

      common /planci/ pib(nrlev),pibs                              ! black-body radiation
      double precision pib, pibs

      common /ray/ dtaur(nrlay), plr(2,nrlay)                         ! rayleigh
      double precision dtaur, plr

      common /tmp2/ as(mbs),ee(mbir)                            ! albedo, and emissivity of the surface
      double precision as, ee

      common /umcon/ umco2,umch4,umn2o                          ! gas concentrations
      double precision umco2,umch4,umn2o

! Local scalars:
      integer i, ib, ig, ip, ix, iz, ibanf, j, l
      double precision hk
      double precision s0
      double precision zbsca
      double precision zdopr
      double precision zfni, zfnip
      double precision zfuq1, zfuq2
      double precision zx0

! Local arrays:
      double precision dlam(7,nrlev,mb)   ! array with solar/ir rad. fluxes
      double precision zgaer(nrlay)       ! needed for Legendre-coefficients
      double precision sss(nrlev)
      integer kg(mb)                   ! number of cumulative probabilities
      data kg / 10, 8, 12, 7, 12, 5,
     &     2, 3, 4, 4, 3, 5, 2, 10, 12, 7, 7, 8 /
!- End of header ---------------------------------------------------------------

      s0=1355.3 ! jjb this is the sum of the 4 former wavelength bands (1018.3+230.1+39.2+67.7)
                !      see paper from Loughlin et al, (1997) QJRMS vol. 123, pp. 1985-2007, table 1
                !     SHOULD THIS BE UPDATED? see also zfuq1 below: 1340.0
c1360.3
! concentrations of trace gases ! jjb improve by reading in the actual values if chemistry is activated
      umco2=330.
      umch4=1.6
      umn2o=0.28

!
! ----------------------------------------------------------------------
! cloud fractions: frr includes the cloud fraction section taken from
! the new PIFM code.
!
      call frr(icld)

! Initialisation or arrays
      fs1(:) = 0.
      fs2(:) = 0.
      ss(:)  = 0.
      sss(:) = 0.
      totds(:) = 0.
      fsn(:) = 0.
      fl1(:) = 0.
      fl2(:) = 0.
      fln(:) = 0.
      totrad(:,:) = 0.
      dlam(:,:,:) = 0.

! at night: no solar radiation (ibanf=7)
      if ( u0 <= 1.0e-2 ) then
         ibanf=mbs + 1
      else
         ibanf=1
      endif

      zdopr=2.*rho(1)


! loop for 18 spectral bands
!---------------------------
      do ib=ibanf,mb

! rayleigh scattering

         if (ib <= mbs) then
            do i=1,nrlay ! Index i goes from top to bottom
               j=nrlev-i    ! Index j goes from bottom to top
               dtaur(i)=berayl(ib)*thk(j)*(rho(j+1)+rho(j))/zdopr
               plr(1,i)=0.
               plr(2,i)=0.5
            enddo
         else
            do i=1,nrlay
               dtaur(i)=0.0
               plr(1,i)=0.0
               plr(2,i)=0.0
            enddo
         endif

! ......................................................................
! aerosols (taer = optical depth)
!
!     First define local variable zgaer (2 cases, solar or IR), which is
!     then used to calculate the output of this section: 
!     taer, waer and plaer, for use in SR tau.
!
!     In the solar case (band 1 only), also calculate variables for the
!     photolysis model: taer_s, taer_a and ga_pl

         if (ib <= mbs) then
            do iz=1,nrlay        ! Index iz goes from top to bottom
               j=nrlev-iz          ! Index j goes from bottom to top
               zbsca=bea(ib,j)-baa(ib,j)
               taer(iz)=bea(ib,j)*thk(j)
               if (ib == 1) then
                  taer_s(iz)=zbsca*thk(j)            ! photolysis
                  taer_a(iz)=baa(ib,j)*thk(j)       ! photolysis
               endif
               if (zbsca+(dtaur(iz)/thk(j)) >= 1.0e-20) then
                  zgaer(iz)=ga(ib,j)*zbsca/(zbsca+(dtaur(iz)/thk(j)))
                  if(ib == 1) ga_pl(iz)=ga(ib,j)     ! photolysis
               else
                  zgaer(iz)=0.0
                  if(ib == 1) ga_pl(iz)=0.0          ! photolysis
               endif
               if (bea(ib,j) > 1.e-20) then
                  waer(iz)=1.0-(baa(ib,j)/bea(ib,j))
               else
                  waer(iz)=0.0
               endif
               do l=1,2
                  plaer(l,iz)=dble(2*l+1)*zgaer(iz)**l
               enddo
            enddo
         else
            do iz=1,nrlay        ! Index iz goes from top to bottom
               j=nrlev-iz          ! Index j goes from bottom to top
               zbsca=bea(ib,j)-baa(ib,j)
               taer(iz)=bea(ib,j)*thk(j)
               zgaer(iz)=ga(ib,j)
               if (bea(ib,j) > 1.e-20) then
                  waer(iz)=1.0-(baa(ib,j)/bea(ib,j))
               else
                  waer(iz)=0.0
               endif
               do l=1,2
                  plaer(l,iz)=dble(2*l+1)*zgaer(iz)**l
               enddo
            enddo
         endif
! ......................................................................



         ! water droplets
         call water(ib,icld)

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
            call tau(ib)

            if ( ib <= mbs ) then

! solution of short-wave radiative transfer equation
!---------------------------------------------------

! radiation fluxes:
! ss    is direct solar downward
! sss   is direct solar downward + Delta-Eddington peak
! fs2   is diffusive solar downward
! fs1   is total diffusive solar upward
               call kurzw(ib,u0)
               call jeanfr(ib)
               do i=1, nrlev  ! Index i goes from ground to top
                  j=nrlev-i+1 ! Index j goes from top to bottom
                  ss(i)=ss(i)+(sf(j)+sw(j))*hk
                  sss(i)=sss(i)+(ssf(j)+ssw(j))*hk
                  fs1(i)=fs1(i)+(f1f(j)+f1w(j))*hk
                  fs2(i)=fs2(i)+(f2f(j)+f2w(j))*hk
                  dlam(1,i,ib)=dlam(1,i,ib)+(sf(j)+sw(j))*hk
                  dlam(2,i,ib)=dlam(2,i,ib)+(ssf(j)+ssw(j))*hk
                  dlam(3,i,ib)=dlam(3,i,ib)+(f1f(j)+f1w(j))*hk
                  dlam(4,i,ib)=dlam(4,i,ib)+(f2f(j)+f2w(j))*hk
               enddo

            else

! solution of long-wave radiative transfer equation
!--------------------------------------------------

! radiation fluxes:
! fl2   is longwave downward
! fl1   is longwave upward
!
! dlam contains the different contributions from solar diffuse and ir
! in each of the middle index (see above). the first index is the
! wavelength bin identifier and the last is the height.

               call langw(ib)
               call jeanfr(ib)
               do i=1, nrlev  ! Index i goes from ground to top
                  j=nrlev-i+1 ! Index j goes from top to bottom
                  fl1(i)=fl1(i)+(pib(j)-f1f(j)-f1w(j))*hk
                  fl2(i)=fl2(i)+(pib(j)-f2f(j)-f2w(j))*hk
                  dlam(5,i,ib)=dlam(5,i,ib)+(pib(j)-f1f(j)-f1w(j))*hk
                  dlam(6,i,ib)=dlam(6,i,ib)+(pib(j)-f2f(j)-f2w(j))*hk
                  dlam(7,i,ib)=dlam(7,i,ib)+pib(j)*hk
               enddo
            endif
         enddo ! ig
      enddo ! ib

!---------------------------------------------------------------------------

! flux correction for given solar constant s0 and emissivity ee

      zfuq1=s0/ 1340.0
      zfuq2=pibs * 0.03 *  ee(mbir)

      ! Solar wavelength bands corrections
      do i=1,nrlev ! This loop goes from ground to top
         ss(i)=ss(i)*zfuq1
         sss(i)=sss(i)*zfuq1
         fs1(i)=fs1(i)*zfuq1
         fs2(i)=fs2(i)*zfuq1
      enddo
      do ib=1,mbs ! optimised loop order: innermost is leftmost
         do i=1,nrlev
            do ix=1,4
               dlam(ix,i,ib)=dlam(ix,i,ib)*zfuq1
            enddo
         enddo
      enddo

! dummy variables for additional output:
! totds is total solar downward
! fsn   is total solar net radiation flux
! fln   is total ir net radiation flux
! f1 and f2 are total upward / downward radiation fluxes
! fn    is total net radiation flux

      do i=1,nrlev
         ! Delta-Eddington peak (sss-ss) is added to fs2
         totds(i)=sss(i)+fs2(i)
         fs2(i)=totds(i)-ss(i)
         fsn(i)=fs1(i)-totds(i)
         fl1(i)=fl1(i)+zfuq2
         dlam(5,i,18)=dlam(5,i,18)+zfuq2 ! jjb mysterious ! fl2?
         fln(i)=fl1(i)-fl2(i)
      enddo

      ! surface radiation fluxes
      flgeg=fl2(1)
      fnseb=fs2(1)+ss(1)-fs1(1)

! calculation of heating rates {hr},
! dtdts and dtdtl are the solar / ir heating rates (currently not used)
      do i=1,nrlay ! This loop goes from ground to top
         ip=i+1
         if (i == 1) zfni=fl1(1)-fl2(1)+fs1(1)-ss(1)-fs2(1)
         zfnip=fl1(ip)-fl2(ip)+fs1(ip)-ss(ip)-fs2(ip)
         zx0=(thk(i)*(rho(ip)+rho(i))*502.5)
         dtdts(i)=(fs1(i)-ss(i)-fs2(i)-fs1(ip)+ss(ip)+fs2(ip))/zx0
         dtdtl(i)=(fl1(i)-fl2(i)-fl1(ip)+fl2(ip))/zx0
         hr(i)=(zfni-zfnip)/zx0
         zfni=zfnip
      enddo

! ------------------------------------------------------------------------------
      do i=1,nrlay
         do ib=1,mbs
            totrad(ib,i)=((dlam(2,i,ib)+dlam(2,i+1,ib))/(2.0*u0))+
     &                   (dlam(3,i,ib)+dlam(3,i+1,ib)+dlam(4,i,ib)
     &                   +dlam(4,i+1,ib))
         enddo
      enddo
      do i=1,nrlay
         do ib=mbs+1,mb
!            if(i == 1) then                                 ! jjb problem here, same as below
               totrad(ib,i)=-(dlam(7,i,ib)+dlam(7,i+1,ib))*2
     &                      +dlam(5,i,ib)+dlam(6,i,ib)
     &                      +dlam(5,i+1,ib)+dlam(6,i+1,ib)
!            else
!               totrad(ib,i)=-(dlam(7,i,ib)+dlam(7,i+1,ib))*2
!     &                      +dlam(5,i,ib)+dlam(6,i,ib)
!     &                      +dlam(5,i+1,ib)+dlam(6,i+1,ib)
!            endif
         enddo
      enddo

! ss    is down sol
! fs2   is down diff
! totds is total downward
! fs1   is up diff
! fsn   is total net
!
! fl2   is longwave down
! fl1   is longwave up
!
! dlam now contains the different contributions from solar diffuse and ir
! in each of the middle index (see above). the first index is the
! wavelength bin identifier and the last is the height.
!
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

      subroutine frr(icld)
!
! Description:
!    Calculation of random overlapping cloud coverage.
! -----------------------------------------------------------------------------
!    bb: cloud free part, cc: cloudy part
!
!    See Zdunkowski et al., Betrei. Phys. Atmosph. Vol. 55 No. 3, August 1982
!        pp. 215-238 (see page 219: coefficients b_1,i to b_4,i).
!    And reference therein.
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      10/2016   Rewritten with if...else... structure  <Josue Bock>
!
! 1.1      07/2016   Header including "USE ... ONLY"        <Josue Bock>
!                    Declarations and implicit none
!
! 1.0       ?        Original code.                         <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      integer, intent(in) :: icld

! Local scalars:
      integer i, ip, j, k

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay), ! only: frac
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /part/ cc(4,nrlay),bb(4,nrlay)
      double precision cc, bb

!- End of header ---------------------------------------------------------------

      if(icld == 0) then
         do i=1,nrlay
            ip=i+1
            j=nrlay-i+1
            if(ip >= lcl.and.ip <= lct) then ! lcl and lct are defined in the vertical grid including
               frac(j)=1.0                   ! the first infinitesimally thin layer (ground) which is
            else                             ! disregarded when dealing with the radiative code
               frac(j)=0.0
            endif
         end do
      endif

      do j=1,nrlay
         i=j-1
         k=j+1

         if(j==1)then
            bb(1,1)=1.
            bb(3,1)=1.
         else
            if(frac(i) > 0.) then
               if(frac(i) < 1.) then
                  if(frac(i) < frac(j)) then
                     bb(3,j)=1.
                     bb(1,j)=(1.-frac(j))/(1.-frac(i))
                  else
                     bb(3,j)=frac(j)/frac(i)
                     bb(1,j)=1.
                  end if
               else             ! frac(i) = 1
                  bb(3,j)=frac(j)
                  bb(1,j)=1.
               end if
            else                ! frac(i) = 0
               bb(3,j)=1.
               bb(1,j)=1.-frac(j)
            end if
         end if

         if(j==nrlay)then
            bb(2,nrlay)=1.
            bb(4,nrlay)=1.
         else
            if(frac(k) > 0.) then
               if(frac(k) < 1.) then
                  if(frac(j) > frac(k)) then
                     bb(4,j)=1.
                     bb(2,j)=(1.-frac(j))/(1.-frac(k))
                  else
                     bb(4,j)=frac(j)/frac(k)
                     bb(2,j)=1.
                  end if
               else             ! frac(k) = 1
                  bb(4,j)=frac(j)
                  bb(2,j)=1.
               end if
            else                ! frac(k) = 0
               bb(4,j)=1.
               bb(2,j)=1.-frac(j)
            end if
         end if

      end do

      do i=1,nrlay
         do j=1,4
            cc(j,i)=1.-bb(j,i)
         end do
      enddo

      end subroutine frr

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine water(ib,icld)
!
! Description:
! *********************************************************************
!   Calculation of optical depth {t2w}, single scattering albedo {w2w} and
!   Legendre coefficients {pl2w} of the phase-function for Mie-scattering
!   of cloud droplets. (eqs. 4.25 - 4.27 diss. Fu (1991))
!------------------------------------------------------------------------------
!   rho2 and reff are the partial density of cloud water and the effective
!   radius of the cloud droplets. Effective radii lower or higher than the
!   boundary values of 4.18 um and 31.18 um are set to the corresponding
!   boundary values. Clouds with rho2 less than 10**-5 kg m**-3 are not 
!   calculated.
! *********************************************************************
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1       10/2016  Correction of size of array r(nkt)                                <Josue Bock>
!                    Avoided a division by 0 in w2w and gg definitions, by slighly changing
!                        the tests ( <=  ret(1), then  >  ret(k))
!                        This also avoids that the calculations after "30 continue" could be
!                        done twice (case rew(i) == ret(k)).
!                    Missing initialisation for 
!
!           07/2016  Header including "USE ... ONLY"                                   <Josue Bock>
!                    pi no longer hardcoded
!
! 1.0       ?        Original code.                                                    <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE constants, ONLY :
! Imported Parameters:
     &     pi,
     &     rhow => water_density ! [kg/m**3]

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nrlay,
     &     nrlev,
     &     nka,
     &     nkt,
     &     ncw,
     &     mb,
     &     mbs

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      integer, intent(in) :: ib    ! current wavelength band
      integer, intent(in) :: icld  ! clouds (/=0) or not (=0)

! Local scalars:
      integer i, j            ! loop indexes (vertical)
      integer ia, jt          ! 2D aerosol grid indexes
      integer jl, k
      double precision gg                 ! interpolation
      double precision znum,zdenom,zfix   ! calculation of effective drop radius

! Local arrays:
      double precision dz(nrlay) ! layer thickness, from bottom to top
      double precision r(nkt)
      double precision rew(nrlay), rho2w(nrlay)

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

      common /cb16/ u0,albedo(mbs),thk(nrlay)
      double precision u0, albedo, thk

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &     e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      double precision ff, fsum
      integer nar

      common /h2o/ t2w(nrlay), w2w(nrlay), pl2w(2,nrlay)
      double precision t2w, w2w, pl2w

      common /was1/ ret(ncw),r2wt(ncw),b2wt(ncw,mb),w2wt(ncw,mb), ! 't' = tabulated values
     &     g2wt(ncw,mb)
      double precision ret, r2wt, b2wt, w2wt, g2wt

!- End of header ---------------------------------------------------------------

      if(icld == 0) then
         t2w(:)=0.0
         w2w(:)=0.0
         pl2w(:,:)=0.0
         return
      else
         ! Initialisations:
         ! ----------------
         ! Define local name for layer thicknesses
         do i=1,nrlay
            dz(i)=thk(i)
         enddo
         do i=1,nrlay
            rho2w(i)=rho2(i)
         enddo
         rew(:) = 0.d0
! -------------------------------------------------------
         do i=lcl,lct
            znum=0.0
            zdenom=0.0
            do jt=1,nkt
               r(jt)=0.0001*(3*e(jt)/(4*pi*
     &              rhow))**0.333333333
               zfix=0.0
               do ia=1,nka
                  zfix=zfix+ff(jt,ia,i)
               enddo
               znum =znum+(r(jt)**3)*zfix
               zdenom=zdenom+(r(jt)**2)*zfix
            enddo
            rew(i)=znum/zdenom
         enddo
! -------------------------------------------------------
         do i=1,nrlay ! here, i goes frop the top to the ground
            j=nrlev-i   ! here, j goes frop the ground to the top

! no clouds
            if ( rho2w(j) < 1.0e-5 ) then
               t2w(i)=0.0
               w2w(i)=0.0
               do jl=1,2
                  pl2w(jl,i)=0.0
               end do
            else

! interpolation for given reff and rho2 from tabulated values
! {ret} and {r2wt}
! lower limit 4.18 um
               if ( rew(j) <= ret(1) ) then
                  t2w(i)=dz(i) * rho2w(j)*b2wt(1,ib)/r2wt(1)
                  w2w(i)=w2wt(1,ib)
                  do jl=1,2
                     pl2w(jl,i)=dble(2*jl+1)*g2wt(1,ib)**jl
                  end do

! upper limit 31.18 um
               elseif ( rew(j) >= ret(ncw) ) then
                  t2w(i)=dz(i) * rho2w(j)*b2wt(ncw,ib)/r2wt(ncw)
                  w2w(i)=w2wt(ncw,ib)
                  do jl=1,2
                     pl2w(jl,i)=dble(2*jl+1)*g2wt(ncw,ib)**jl
                  end do

! linear interpolation for values between the limits
               else
                  k = 1
                  do while (rew(j) > ret(k+1))
                     k = k + 1
                  end do

                  t2w(i)=dz(i) * rho2w(j) * ( b2wt(k,ib)
     &                 / r2wt(k) +( b2wt(k+1,ib) / r2wt(k+1) -
     &                 b2wt(k,ib) / r2wt(k) ) /
     &                 ( 1.0 / ret(k+1) - 1.0 / ret(k) )
     &                 * ( 1.0 / rew(j)
     &                 - 1.0 / ret(k) ) )
                  w2w(i)=w2wt(k,ib) + ( w2wt(k+1,ib) - w2wt(k,ib) ) /
     &                 ( ret(k+1) - ret(k) ) * ( rew(j) - ret(k) )
                  gg=g2wt(k,ib) + ( g2wt(k+1,ib) - g2wt(k,ib) ) /
     &                 ( ret(k+1) - ret(k) ) * ( rew(j) - ret(k) )
                  do jl=1,2
                     pl2w(jl,i)=dble(2*jl+1)*gg**jl
                  end do
               endif
            endif
         end do
      endif

      end subroutine water

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine gascon ( ib )
!
! Description: continuum absorption of water vapour, bands 11 - 17
!
!              for partial reference, see Loughlin et al. 1997, QJRMS vol 123 (543), pp. 1985-2007
!              doi: 10.1002/qj.49712354311
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      10/2016   implicit none
!          07/2016   Header including "USE ... ONLY"    <Josue Bock>
!
! 1.0       ?        Original code.                     <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     mb,
     &     nrlay

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      integer, intent(in) :: ib

! Local scalars:
      integer i

! Local arrays:
      double precision vv(mb)
      data vv / 10*0.0, 1175.0, 1040.0, 890.0, 735.0,
     &     605.0, 470.0, 340.0, 0.0 /

! Common blocks:
      common /con/ tgcon(nrlay)
      double precision tgcon
!- End of header ---------------------------------------------------------------

      if ( ib >= 11 .and. ib <= 17 ) then
         call qopcon ( vv(ib) )
      else
         do i=1, nrlay
            tgcon(i)=0.0
         enddo
      endif

      end subroutine gascon

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine qopcon ( vv )
!
! Description:
! **************************************************************
!   e-Type- and continuum absorption after eq (A.19) of
!   dissertation Fu. (1991)
!   vv are the central wavenumbers of the spectral bands in in cm**-1.
!   All other variables in MKS: 
!   p1  =  partial pressure of water vapour in Pa
!   ff  =  continuum/ e-Type- absorption coefficients in m**2 kg**-1
!   delta-tau  =  ff * xm1 *rho * delta-z  =  - ff * xm1 * delta-p / gE,
!   ff * xm1 arithmetical mean values
! **************************************************************
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      10/2016   Removed arithmetic do loops       <Josue Bock>
!                    Implicit none
!          07/2016   Header including "USE ... ONLY"
!
! 1.0       ?        Original code.                    <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision, intent(in) :: vv

! Local scalars:
      integer i,j
      double precision r,x,y,zz
      double precision s,w

! Local arrays:
      double precision ff(nrlev), p1(nrlev)

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

      common /con/ tgcon(nrlay)
      double precision tgcon
!- End of header ---------------------------------------------------------------

      x=418.
      y=557780.
      zz=0.00787
      r=0.00002
      s=( x + y * exp ( - zz * vv ) ) / 101325.
      do j=1,nrlev                                            ! jjb should this be nrlay instead ???? CHECK !
         i=nrlev+1-j
         p1(i)=p(j) * xm1(j) / ( 0.622 + 0.378 * xm1(j) )
         w=exp ( 1800.0 / t(j) - 6.08108 )
         ff(i)=s * ( p1(i)/100. + r * p(j) ) * w
      enddo
      do i=1,nrlay
         j=nrlev+1-i
         tgcon(i)=( ff(i) * xm1(j) + ff(i+1) * xm1(j-1) )*
     &            ( p(j-1) - p(j) ) * 0.00509892
      enddo

      end subroutine qopcon

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine planck(ib)
!
! Description:
! ***************************************************************
!   Calculation of over {ib} integrated Planck-function for every model layer,
!   and of the radiation temperature of the surface.
! ***************************************************************
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header including "USE ... ONLY"   <Josue Bock>
!                    implicit none
!                    pi defined through a module
!
! 1.0       ?        Original code.                    <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE constants, ONLY :
! Imported Parameters:
     & pi

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev,
     &     mbs

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      integer, intent(in) :: ib

! Local scalars:
      integer i

! Local arrays:
      double precision wvl(13)
      data wvl / 2200., 1900., 1700., 1400., 1250., 1100., 980.,
     &              800., 670., 540., 400., 280., 0. /

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

      common /planci/ pib(nrlev),pibs
      double precision pib, pibs

! External function:
      double precision plkavg
!- End of header ---------------------------------------------------------------

      if(ib <= mbs) return
      do i=1,nrlev
         pib(i)=pi*plkavg(wvl(ib-5),wvl(ib-6),t(nrlev-i+1))
      enddo
      pibs=pi*plkavg(wvl(ib-5),wvl(ib-6),ts)

      end subroutine planck

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      double precision function plkavg( WNUMLO, WNUMHI, xT )
!
! Description:
! ***************************************************************
!   Calculation of black body radiation for given temperature {xt}
!   For low wavenumbers the calculation is done with power series, for high
!   wavenumbers with exponential series. Critical wavenumber is 'vcut'. 
! ***************************************************************
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header, comments and cleaning  <Josue Bock>
!
! 1.0       ?        Original code.                 <Andreas Bott>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
      implicit double precision(a-h,o-z)

      double precision     xT, WNUMLO, WNUMHI
      PARAMETER  ( A1=1./3., A2=-1./8., A3=1./60., A4=-1./5040.,
     &             A5=1./272160., A6=-1./13305600. )
      INTEGER  SMALLV

      REAL*8     C2, CONC, D(2), EPSIL, EX, MV, P(2), SIGMA, SIGDPI,
     &         V(2), VCUT, VCP(7), VSQ

      SAVE     CONC, VMAX, EPSIL, SIGDPI
      DATA     C2 / 1.438786 /,  SIGMA / 5.67032E-8 /,
     &         VCUT / 1.5 /, VCP / 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /
      F(X)=X**3 / ( EXP(X) - 1. )
!
!
      PI=2. * ASIN( 1.0 )
!      VMAX=DLOG( D1MACH(2) )
!      EPSIL=D1MACH(3)
      VMAX=DLOG( huge(1d0) ) ! jjb this should probably be tested: shouldn't that be a log10? If so, could be replaced by range(1d0)
      EPSIL=epsilon(1d0)

      SIGDPI=SIGMA / PI
      CONC=15. / PI**4
!
      IF( xT < 0.0 .OR. WNUMHI <= WNUMLO .OR. WNUMLO < 0.)
     &    print *,'PLKAVG--TEMPERATURE OR WAVENUMS. WRONG'
!
      IF ( xT < 1.E-4 )  THEN
         PLKAVG=0.0
         RETURN
      ENDIF
!
      V(1)=C2 * WNUMLO / xT
      V(2)=C2 * WNUMHI / xT
      IF ( V(1) > EPSIL .AND. V(2) < VMAX .AND.
     &     (WNUMHI-WNUMLO)/WNUMHI < 1.E-2 )  THEN
!
!                          ** WAVENUMBERS ARE VERY CLOSE.  GET INTEGRAL
!                          ** BY ITERATING SIMPSON RULE TO CONVERGENCE.
         HH=V(2) - V(1)
         OLDVAL=0.0
         xx1=V(1)
         xx2=V(2)
         VAL0=F( xx1 ) + F( xx2 )
!
         DO  2  N=1, 10
            DEL=HH / (2*N)
            VAL=VAL0
            DO  1  K=1, 2*N-1
               VAL=VAL + 2*(1+MOD(K,2)) * F( V(1) + K*DEL )
    1       CONTINUE
            VAL=DEL/3. * VAL
            IF ( ABS( (VAL-OLDVAL)/VAL ) <= 1.E-6 )  GO TO 3
            OLDVAL=VAL
    2    CONTINUE
         print *,'PLKAVG--SIMPSON RULE DIDNT CONVERGE'
    3    PLKAVG=SIGDPI * xT**4 * CONC * VAL
         RETURN
      END IF
      SMALLV=0
      DO  50  I=1, 2
         IF( V(I) < VCUT )  THEN
! ** USE POWER SERIES
            SMALLV=SMALLV + 1
            VSQ=V(I)**2
            P(I)= CONC * VSQ * V(I) * ( A1 + V(I) * ( A2 + V(I) *
     &           ( A3 + VSQ * ( A4 + VSQ * ( A5 + VSQ*A6 ) ) ) ) )
         ELSE
! ** USE EXPONENTIAL SERIES
            MMAX=0
! ** FIND UPPER LIMIT OF SERIES
 20         MMAX=MMAX + 1
            IF ( V(I) < VCP( MMAX ) )  GO TO 20
            EX=EXP( - V(I) )
            EXM=1.0
            D(I)=0.0
!
            DO  30  M=1, MMAX
               MV=M * V(I)
               EXM=EX * EXM
               D(I)=D(I) +
     &              EXM * ( 6. + MV*( 6. + MV*( 3. + MV ) ) ) / M**4
 30         CONTINUE
!
            D(I)=CONC * D(I)
         END IF
!
 50   CONTINUE
!
      IF ( SMALLV == 2 ) THEN
! ** WNUMLO AND WNUMHI BOTH SMALL
         PLKAVG=P(2) - P(1)
      ELSE IF ( SMALLV == 1 ) THEN
! ** WNUMLO SMALL, WNUMHI LARGE
         PLKAVG=1. - P(1) - D(2)
      ELSE
! ** WNUMLO AND WNUMHI BOTH LARGE
         PLKAVG=D(1) - D(2)
      END IF
      PLKAVG=SIGDPI * xT**4 * PLKAVG
      IF( PLKAVG == 0.0 )
     &     print *,'PLKAVG--RETURNS ZERO; POSSIBLE UNDERFLOW'

      END function plkavg

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
! jjb 24/07/2016 unused, thus commented
!
c$$$      subroutine plancktab(ib)
c$$$!
c$$$! Description:
c$$$! ***************************************************************
c$$$!   Berechnung der ueber den Bereich {ib} integrierten Planck-
c$$$!   Funktion auf allen Rechenflaechen und der Erdoberflaechen-
c$$$!   temperatur (kann ungleich der Temperatur der Luft am Boden
c$$$!   sein)
c$$$! ***************************************************************
c$$$!
c$$$
c$$$!
c$$$! History:
c$$$! Version   Date     Comment
c$$$! -------   ----     -------
c$$$! 1.1      07/2016   Header including "USE ... ONLY"  <Josue Bock>
c$$$!
c$$$! 1.0       ?        Original code.                   <unknown>
c$$$!
c$$$! Code Description:
c$$$!   Language:          Fortran 77 (with Fortran 90 features)
c$$$!
c$$$! Declarations:
c$$$! Modules used:
c$$$
c$$$      USE global_params, ONLY :
c$$$! Imported Parameters:
c$$$     &     nrlay,
c$$$     &     nrlev,
c$$$     &     mbs
c$$$
c$$$      implicit double precision(a-h,o-z)
c$$$
c$$$! Common blocks:
c$$$      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
c$$$     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
c$$$      double precision t,p,rho,xm1,rho2,frac,ts
c$$$      integer ntypa,ntypd
c$$$
c$$$      common /planci/ pib(nrlev),pibs
c$$$      double precision pib, pibs
c$$$!- End of header ---------------------------------------------------------------
c$$$
c$$$      do 1 i=1,nrlev
c$$$      pib(i)=0.
c$$$ 1         continue
c$$$      pibs=0.
c$$$      if(ib <= mbs) return
c$$$      ibir=ib-mbs
c$$$      do 2 i=1,nrlev
c$$$      pib(i)=fst4(ibir,t(nrlev+1-i)) ! function fst4 has to be uncommented if used
c$$$ 2         continue
c$$$      pibs=fst4(ibir,ts)
c$$$
c$$$      end subroutine plancktab

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
! jjb 25/07/2016: used only by SR plancktab, whose call is commented
!                  -> not used anymore, thus commented out
!
c$$$      double precision function fst4(ib,t)
c$$$!
c$$$! Description:
c$$$! ***************************************************************
c$$$!   Interpolation of black body radiation from tabled values {pibtab}
c$$$!   for 35 temperatures {ttab}.
c$$$!   Daten in Block Data PLANC1.
c$$$! ***************************************************************
c$$$!
c$$$
c$$$!
c$$$! History:
c$$$! Version   Date     Comment
c$$$! -------   ----     -------
c$$$! 1.1      07/2016   Header including "USE ... ONLY"   <Josue Bock>
c$$$!
c$$$! 1.0       ?        Original code.                     <unknown>
c$$$!
c$$$! Code Description:
c$$$!   Language:          Fortran 77 (with Fortran 90 features)
c$$$!
c$$$! Declarations:
c$$$
c$$$      implicit double precision(a-h,o-z)
c$$$
c$$$! Common blocks:
c$$$      common /plancd/ ttab(35),pibtab(12,35)
c$$$!- End of header ---------------------------------------------------------------
c$$$      do 1 k=1,35
c$$$         kk=k
c$$$         if(t-ttab(k)) 3,2,1
c$$$ 1    continue
c$$$ 2    fst4=pibtab(ib,kk)
c$$$      return
c$$$ 3    continue
c$$$      kkm=kk-1
c$$$      fst4=pibtab(ib,kkm)+(pibtab(ib,kk)-pibtab(ib,kkm))
c$$$     &     *(t-ttab(kkm))*0.2
c$$$
c$$$      end function fst4

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine gase ( ib, ig, hk )
!
! Description:
!    Calculation of gas absorption for different gases and 18 spectral bands.
!    The subroutines qks (solar) and qki (ir) calculate the absorption
!    coefficients 'fkg' from tabled values.
!    tg(nrlay) are the optical depths due to nongray gaseous absorption, in
!    nrlay layers for a given band ib and cumulative probability ig. They are 
!    calculated from fkg in the subroutines qop{'gas formula'}.

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header                                              <Josue Bock>
!                    Explicit array index (start:end) in several calls
!                    Removal of labeled do-loops, cleaning
!
! 1.0       ?        Original code.                                      <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      integer, intent(in) :: ib, ig
! Scalar arguments with intent(out):
      double precision, intent(out) :: hk

! Local scalars:
      integer i,j ! loop indexes
      double precision fk

! Common blocks:
      common /band1/ hk1(10), fk1o3(10)
      double precision hk1, fk1o3
      common /band2/ hk2(8), c2h2o(3,11,8)
      double precision hk2, c2h2o
      common /band3/ hk3(12), c3h2o(3,11,12)
      double precision hk3, c3h2o
      common /band4/ hk4(7), c4h2o(3,11,7)
      double precision hk4, c4h2o
      common /band5/ hk5(12), c5h2o(3,11,12)
      double precision hk5, c5h2o
      common /band6/ hk6(5), c6h2o(3,11,5)
      double precision hk6, c6h2o
      common /band7/ hk7(2), c7h2o(3,19,2)
      double precision hk7, c7h2o
      common /band8/ hk8(3), c8h2o(3,19,3)
      double precision hk8, c8h2o
      common /band9/ hk9(4), c9h2o(3,19,4)
      double precision hk9, c9h2o
      common /band10/ hk10(4),c10h2o(3,19,4),c10ch4(3,19),c10n2o(3,19)
      double precision hk10, c10h2o, c10ch4, c10n2o
      common /band11/ hk11(3),c11h2o(3,19,3),c11ch4(3,19),c11n2o(3,19)
      double precision hk11, c11h2o, c11ch4, c11n2o
      common /band12/ hk12(5), c12o3(3,19,5), c12h2o(3,19)
      double precision hk12, c12o3, c12h2o
      common /band13/ hk13(2), c13h2o(3,19,2)
      double precision hk13, c13h2o
      common /band14/ hk14(10), c14hca(3,19,10), c14hcb(3,19,10)
      double precision hk14, c14hca, c14hcb
      common /band15/ hk15(12), c15hca(3,19,12), c15hcb(3,19,12)
      double precision hk15, c15hca, c15hcb
      common /band16/ hk16(7), c16h2o(3,19,7)
      double precision hk16, c16h2o
      common /band17/ hk17(7), c17h2o(3,19,7)
      double precision hk17, c17h2o
      common /band18/ hk18(8), c18h2o(3,19,8)
      double precision hk18, c18h2o

      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

      common /gas/ tg(nrlay)
      double precision tg

      common /umcon/ umco2,umch4,umn2o
      double precision umco2,umch4,umn2o

! Local arrays:
      double precision fkg(nrlev), fkga(nrlev), fkgb(nrlev), pq(nrlev)
      double precision tg1(nrlay), tg2(nrlay), tg3(nrlay)
!- End of header ---------------------------------------------------------------

      goto ( 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 ) ib
      stop 'SR gase'

! Band 1:
! In this band ( 50000 - 14500 cm**-1 ), we consider the nongray
! gaseous absorption of O3.    619.618 is the solar energy contained in
! the band in units of Wm**-2.
 1    fk = fk1o3(ig)
      call qopo3s ( fk, tg )
      hk=619.618 * hk1(ig)
      goto 20

! Band 2:
! In this band ( 14500 - 7700 cm**-1 ), we consider the nongray
! gaseous absorption of H2O.  484.295 is the solar energy contained in
! the band in units of Wm**-2.
 2    call qks ( c2h2o(1:3,1:11,ig), fkg )
      call qoph2o ( fkg, tg )
      hk=484.295 * hk2(ig)
      goto 20

! Band 3:
! In this band ( 7700 - 5250 cm**-1 ), we consider the nongray
! gaseous absorption of H2O. 149.845 is the solar energy contained in
! the band in units of Wm**-2.
 3    call qks ( c3h2o(1:3,1:11,ig), fkg )
      call qoph2o ( fkg, tg )
      hk=149.845 * hk3(ig)
      goto 20

! Band 4:
! In this band ( 5250 - 4000 cm**-1 ), we consider the nongray
! gaseous absorption of H2O. 48.7302 is the solar energy contained in
! the band in units of Wm**-2.
 4    call qks ( c4h2o(1:3,1:11,ig), fkg )
      call qoph2o ( fkg, tg )
      hk=48.7302 * hk4(ig)
      goto 20

! Band 5:
! In this band ( 4000 - 2850 cm**-1 ), we consider the nongray
! gaseous absorption of H2O. 31.6576 is the solar energy contained in
! the band in units of Wm**-2.
 5    call qks ( c5h2o(1:3,1:11,ig), fkg )
      call qoph2o ( fkg, tg )
      hk=31.6576 * hk5(ig)
      goto 20

! Band 6:
! In this band ( 2850 - 2500 cm**-1 ), we consider the nongray
! gaseous absorption of H2O. 5.79927 is the solar energy contained in
! the band in units of Wm**-2.
 6    call qks ( c6h2o(1:3,1:11,ig), fkg )
      call qoph2o ( fkg, tg )
      hk=5.79927 * hk6(ig)
      goto 20

! Band 7:
! In this band ( 2200 - 1900 cm**-1 ), we consider the nongray
! gaseous absorption of H2O.
 7    call qki ( c7h2o(1:3,1:19,ig), fkg )
      call qoph2o ( fkg, tg )
      hk=hk7(ig)
      goto 20

! Band 8:
! In this band ( 1900 - 1700 cm**-1 ), we consider the nongray
! gaseous absorption of H2O.
 8    call qki ( c8h2o(1:3,1:19,ig), fkg )
      call qoph2o ( fkg, tg )
      hk=hk8(ig)
      goto 20

! Band 9:
! In this band ( 1700 - 1400 cm**-1 ), we consider the nongray
! gaseous absorption of H2O.
 9    call qki ( c9h2o(1:3,1:19,ig), fkg )
      call qoph2o ( fkg, tg )
      hk=hk9(ig)
      goto 20

! Band 10:
! In this band ( 1400 - 1250 cm**-1 ), we consider the overlapping
! absorption of H2O, CH4, and N2O by approach one of Fu(1991).
 10   call qki ( c10h2o(1:3,1:19,ig), fkg )
      call qoph2o ( fkg, tg1 )
      call qki ( c10ch4, fkg )
      call qopch4 ( fkg, tg2 )
      call qki ( c10n2o, fkg )
      call qopn2o ( fkg, tg3 )
      do i=1,nrlay
         tg(i)=tg1(i) + tg2(i)/1.6*umch4 + tg3(i)/0.28*umn2o
      end do
      hk=hk10(ig)
      goto 20

! Band 11:
! In this band ( 1250 - 1100 cm**-1 ), we consider the overlapping
! absorption of H2O, CH4, and N2O by approach one of Fu(1991).
 11   call qki ( c11h2o(1:3,1:19,ig), fkg )
      call qoph2o ( fkg, tg1 )
      call qki ( c11ch4, fkg )
      call qopch4 ( fkg, tg2 )
      call qki ( c11n2o, fkg )
      call qopn2o ( fkg, tg3 )
      do i=1,nrlay
         tg(i)=tg1(i) + tg2(i)/1.6*umch4 + tg3(i)/0.28*umn2o
      end do
      hk=hk11(ig)
      goto 20

! Band 12:
! In this band ( 1100 - 980 cm**-1 ), we consider the overlapping
! absorption of H2O and O3 by approach one of Fu(1991).
 12   call qkio3 ( c12o3(1:3,1:19,ig), fkg )
      call qopo3i ( fkg, tg1 )
      call qki ( c12h2o, fkg )
      call qoph2o ( fkg, tg2 )
      do i=1,nrlay
         tg(i)=tg1(i) + tg2(i)
      end do
      hk=hk12(ig)
      goto 20

! Band 13:
! In this band ( 980 - 800 cm**-1 ), we consider the nongray
! gaseous absorption of H2O.
 13   call qki ( c13h2o(1:3,1:19,ig), fkg )
      call qoph2o ( fkg, tg )
      hk=hk13(ig)
      goto 20

! Band 14
! In this band ( 800 - 670 cm**-1), we consider the overlapping
! absorption of H2O and CO2 by approach two of Fu(1991).
 14   do j=1, nrlev
         i=nrlev+1-j
         if ( p(j) >= 6310. ) then
            pq(i)=xm1(j)
         else
            pq(i)=0.0
         endif
      end do
      call qki ( c14hca(1:3,1:19,ig), fkga )
      call qki ( c14hcb(1:3,1:19,ig), fkgb )
      do i=1, nrlev
         fkg(i)=fkga(i)/330.0*umco2 + pq(i) * fkgb(i)
      end do
      call qophc ( fkg, tg)
      hk=hk14(ig)
      goto 20

! Band 15:
! In this band ( 670 - 540 cm**-1), we consider the overlapping
! absorption of H2O and CO2 by approach two of Fu(1991).
 15   do j=1, nrlev
         i=nrlev+1-j
         if ( p(j) >= 6310. ) then
            pq(i)=xm1(j)
         else
            pq(i)=0.0
         endif
      end do
      call qki ( c15hca(1:3,1:19,ig), fkga )
      call qki ( c15hcb(1:3,1:19,ig), fkgb )
      do i=1, nrlev
         fkg(i)=fkga(i)/330.0*umco2 + pq(i) * fkgb(i)
      end do
      call qophc ( fkg, tg)
      hk=hk15(ig)
      goto 20

! Band 16:
! In this band ( 540 - 400 cm**-1 ), we consider the nongray
! gaseous absorption of H2O.
 16   call qki ( c16h2o(1:3,1:19,ig), fkg )
      call qoph2o ( fkg, tg )
      hk=hk16(ig)
      goto 20

! Band 17:
! In this band ( 400 - 280 cm**-1 ), we consider the nongray
! gaseous absorption of H2O.
 17   call qki ( c17h2o(1:3,1:19,ig), fkg )
      call qoph2o ( fkg, tg )
      hk=hk17(ig)
      goto 20

! Band 18:
! In this band ( 280 - 000 cm**-1 ), we consider the nongray
! gaseous absorption of H2O.
 18   call qki ( c18h2o(1:3,1:19,ig), fkg )
      call qoph2o ( fkg, tg )
      hk=hk18(ig)

 20   continue

      end subroutine gase

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine qks ( coefks, fkg )
!
! Description:
!     Calculation of the absorption coefficients for solar spectral bands.
! 
! fkg(nrlev) are the gaseous absorption coefficients in units of (cm-atm)**-1
! for a given cumulative probability in nrlev layers. coefks(3,11)
! are the coefficients to calculate the absorption coefficient at the
! temperature t for the 11 pressures by
!         ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
! and the absorption coefficient at conditions other than those eleven
! pressures is interpolated linearly with pressure (Fu, 1991).
!   Obwohl der Druck dem U.P in MKS gegeben wird, haben die fkg die
!   originaere Fu-Einheit. stanp ist ebenfalls in Pascal
! *********************************************************************
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header including "USE ... ONLY"   <Josue Bock>
!                    Declarations and implicit none
!
! 1.0       ?        Original code.                    <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Array arguments with intent(in):
      double precision, intent(in) :: coefks(3,11)
! Array arguments with intent(out):
      double precision, intent(out) :: fkg(nrlev)

! Local scalars:
      integer i, i1, j
      double precision x1, x2, y1
! Local arrays:
      double precision stanp(11)
      data stanp / 1000., 1580., 2510., 3980., 6310., 10000.,
     &     15800.,25100.,39800.,63100.,100000. /

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

!- End of header ---------------------------------------------------------------

      do 5 j=1, nrlev
        i1=1
        i=nrlev+1-j

! absorption coefficient for boundary values: pressure lower than 
! 10 hPa or higher than 1000 hPa
        if ( p(j) < stanp(1) ) then
           x1=exp( coefks(1,1) + coefks(2,1) * ( t(j) - 245.0 )
     &          + coefks(3,1) * ( t(j) - 245.0 ) ** 2 )
           fkg(i)=x1 * p(j) / stanp(1)

        elseif ( p(j) >= stanp(11) ) then
           y1=( t(j) - 245.0 ) * ( t(j) - 245.0 )
           x1=exp( coefks(1,10) + coefks(2,10) * ( t(j) - 245.0 )
     &          + coefks(3,10) * y1 )
           x2=exp( coefks(1,11) + coefks(2,11) * ( t(j) - 245.0 )
     &          + coefks(3,11) * y1 )
           fkg(i)=x1 + ( x2 - x1 ) / ( stanp(11) - stanp(10) )
     &          * ( p(j) - stanp(10) )

! linear interpolation of absorption coefficients coefks between reference
! values of pressure in array stanp.
        else
 30        continue
           if ( p(j) >= stanp(i1) ) goto 20
           y1=( t(j) - 245.0 ) * ( t(j) - 245.0 )
           x1=exp( coefks(1,i1-1) + coefks(2,i1-1) * (t(j)-245.0)
     &          + coefks(3,i1-1) * y1 )
           x2=exp( coefks(1,i1) + coefks(2,i1) * ( t(j) - 245.0 )
     &          + coefks(3,i1) * y1 )
           fkg(i)=x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )
     &          * ( p(j) - stanp(i1-1) )
           goto 5
 20        i1=i1 + 1
           goto 30
        endif
 5    continue

      end subroutine qks

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine qki ( coefki, fkg )
!
! Description:
!    Calculation of the absorption coefficients for ir spectral bands.
!
! fkg(nrlev) are the gaseous absorption coefficients in units of (cm-atm)**-1
! for a given cumulative probability in nrlev layers. coefki(3,19)
! are the coefficients to calculate the absorption coefficient at the
! temperature t for the 19 pressures by
!         ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
! and the absorption coefficient at  conditions  other  than  those 19
! pressures is interpolated linearly with pressure (Fu, 1991).
!   Fuer die Einheiten gilt dasselbe wie in U.P. qks
! *********************************************************************
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header including "USE ... ONLY"    <Josue Bock>
!                    Declarations and implicit none
!
! 1.0       ?        Original code.                     <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Array arguments with intent(in):
      double precision, intent(in) :: coefki(3,19)
! Array arguments with intent(out):
      double precision, intent(out) :: fkg(nrlev)

! Local scalars:
      integer i, i1, j
      double precision x1, x2, y1
! Local arrays:
      double precision stanp(19)
      data stanp / 25.1,  39.8,  63.1,  100.,  158.,  251.,
     &     398.,  631., 1000., 1580., 2510., 3980.,
     &     6310.,10000.,15800.,25100.,39800.,63100.,
     &     100000. /

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

!- End of header ---------------------------------------------------------------

      do 5 j=1, nrlev
         i1=1
         i=nrlev+1-j

! absorption coefficient for boundary values: pressure lower than
! 0.25 hPa or higher than 1000 hPa
         if ( p(j) < stanp(1) ) then
            x1=exp( coefki(1,1) + coefki(2,1) * ( t(j) - 245.0 )
     &           + coefki(3,1) * ( t(j) - 245.0 ) ** 2 )
            fkg(i)=x1 * p(j) / stanp(1)
         elseif ( p(j) >= stanp(19) ) then
            y1=( t(j) - 245.0 ) * ( t(j) - 245.0 )
            x1=exp( coefki(1,18) + coefki(2,18) * ( t(j) - 245.0 )
     &           + coefki(3,18) * y1 )
            x2=exp( coefki(1,19) + coefki(2,19) * ( t(j) - 245.0 )
     &           + coefki(3,19) * y1 )
            fkg(i)=x1 + ( x2 - x1 ) / ( stanp(19) - stanp(18) )
     &           * ( p(j) - stanp(18) )

! linear interpolation of absorption coefficients coefks between reference
! values of pressure in array stanp.
         else
 30         continue
            if ( p(j) >= stanp(i1) ) goto 20
            y1=( t(j) - 245.0 ) * ( t(j) - 245.0 )
            x1=exp( coefki(1,i1-1) + coefki(2,i1-1) * (t(j)-245.0)
     &           + coefki(3,i1-1) * y1 )
            x2=exp( coefki(1,i1) + coefki(2,i1) * ( t(j) - 245.0 )
     &           + coefki(3,i1) * y1 )
            fkg(i)=x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )
     &           * ( p(j) - stanp(i1-1) )
            goto 5
 20         i1=i1 + 1
            goto 30
         endif
 5    continue

      end subroutine qki

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine qkio3 ( coefki, fkg )
!
! Description:
!    Calculation of the absorption coefficients for ozone in spectral band 12
!
! fkg(nrlev) are the gaseous absorption coefficients in units of (cm-atm)**-1
! for a given cumulative probability in nrlev layers. coefki(3,19)
! are the coefficients to calculate the absorption coefficient at the
! temperature t for the 19 pressures by
!         ln k = a + b * ( t - 250 ) + c * ( t - 250 ) ** 2
! and the absorption coefficient at  conditions  other  than  those 19
! pressures is interpolated linearly with pressure (Fu, 1991).
!   Fuer die Einheiten gilt dasselbe wie in U.P. qks
! *********************************************************************
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header including "USE ... ONLY"   <Josue Bock>
!                    Declarations and implicit none
!
! 1.0       ?        Original code.                    <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Array arguments with intent(in):
      double precision, intent(in) :: coefki(3,19)
! Array arguments with intent(out):
      double precision, intent(out) :: fkg(nrlev)

! Local arrays:
      integer i, i1, j
      double precision x1, x2, y1
      double precision stanp(19)
      data stanp / 25.1,  39.8,  63.1,  100.,  158.,  251.,
     &     398.,  631., 1000., 1580., 2510., 3980.,
     &     6310.,10000.,15800.,25100.,39800.,63100.,
     &     100000. /

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

!- End of header ---------------------------------------------------------------

      do 5 j=1, nrlev
      i1=1
      i=nrlev+1-j

! absorption coefficient for boundary values: pressure lower than
! 0.25 hPa or higher than 1000 hPa
      if ( p(j) < stanp(1) ) then
         x1=exp( coefki(1,1) + coefki(2,1) * ( t(j) - 250.0 )
     &        + coefki(3,1) * ( t(j) - 250.0 ) ** 2 )
         fkg(i)=x1 * p(j) / stanp(1)
      elseif ( p(j) >= stanp(19) ) then
         y1=( t(j) - 250.0 ) * ( t(j) - 250.0 )
         x1=exp( coefki(1,18) + coefki(2,18) * ( t(j) - 250.0 )
     &        + coefki(3,18) * y1 )
         x2=exp( coefki(1,19) + coefki(2,19) * ( t(j) - 250.0 )
     &        + coefki(3,19) * y1 )
         fkg(i)=x1 + ( x2 - x1 ) / ( stanp(19) - stanp(18) )
     &        * ( p(j) - stanp(18) )

! linear interpolation of absorption coefficients coefks between reference
! values of pressure in array stanp.
      else
 30      continue
         if ( p(j) >= stanp(i1) ) goto 20
         y1=( t(j) - 250.0 ) * ( t(j) - 250.0 )
         x1=exp( coefki(1,i1-1) + coefki(2,i1-1) * (t(j)-250.0)
     &        + coefki(3,i1-1) * y1 )
         x2=exp( coefki(1,i1) + coefki(2,i1) * ( t(j) - 250.0 )
     &        + coefki(3,i1) * y1 )
         fkg(i)=x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )
     &        * ( p(j) - stanp(i1-1) )
         goto 5
 20      i1=i1 + 1
         goto 30
      endif
 5    continue

      end subroutine qkio3

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine qopo3s ( fk, tg )
!
! Description:
!     Calculation of ozone absorption in band 1 ( 50000 - 14500 cm**-1 )
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header including "USE ... ONLY"  <Josue Bock>
!                    Declarations and implicit none
!
! 1.0       ?        Original code.                   <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision, intent(in) :: fk
! Array arguments with intent(out):
      double precision, intent(out) :: tg(nrlay)

! Local scalars:
      double precision fq
      integer i,j

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

      common /ozon/ qmo3(nrlev)
      double precision qmo3

!- End of header ---------------------------------------------------------------

      fq=2.3808 * fk
      do i=1,nrlay
         j=nrlev+1-i
         tg(i)=( qmo3(j) + qmo3(j-1) ) * ( p(j-1) - p(j) ) * fq
      end do

      end subroutine qopo3s

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine qoph2o ( fkg, tg )
!
! Description:
!    Calculation of optical depth of water vapour absorption in bands
!    2 ( 14500 - 7700 cm**-1 ), 3  ( 7700 - 5250 cm**-1 ),
!    4  ( 5250 - 4000 cm**-1 ), 5  ( 4000 - 2850 cm**-1 ),
!    6  ( 2850 - 2500 cm**-1 ), 7  ( 2200 - 1900 cm**-1 ),
!    8  ( 1900 - 1700 cm**-1 ), 9  ( 1700 - 1400 cm**-1 ),
!    10 ( 1400 - 1250 cm**-1 ), 11 ( 1250 - 1100 cm**-1 ),
!    12 ( 1100 -  980 cm**-1 ), 13 (  980 -  800 cm**-1 ),
!    16 (  540 -  400 cm**-1 ), 17 (  400 -  280 cm**-1 ),
!    and 18 ( 280 - 0 cm**-1 ).
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header including "USE ... ONLY"   <Josue Bock>
!          10/2016   Declarations and implicit none
!
! 1.0       ?        Original code.                    <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision, intent(in) :: fkg(nrlev)
! Array arguments with intent(out):
      double precision, intent(out) :: tg(nrlay)

! Local scalars:
      integer i,j

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

!- End of header ---------------------------------------------------------------

      do i=1,nrlay
         j=nrlev-i+1
         tg(i)=( fkg(i) * xm1(j) + fkg(i+1) * xm1(j-1) )
     &        * ( p(j-1) - p(j) ) * 6.349205
      end do

      end subroutine qoph2o

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine qopch4 ( fkg, tg )
!
! Description:
!   Calculation of optical depth of methane (CH4) absorption in bands
!   10 ( 1400 - 1250 cm**-1 ) and 11 ( 1250 - 1100 cm**-1 )
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      10/2016   Declarations and implicit none   <Josue Bock>
!          07/2016   Header including "USE ... ONLY"
!
! 1.0       ?        Original code.                    <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision, intent(in) :: fkg(nrlev)
! Array arguments with intent(out):
      double precision, intent(out) :: tg(nrlay)

! Local scalars:
      integer i,j

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

!- End of header ---------------------------------------------------------------

      do i=1,nrlay
         j=nrlev+1-i
         tg(i)=( fkg(i)+fkg(i+1) ) * ( p(j-1)-p(j) ) * 6.3119e-6
      end do

      end subroutine qopch4

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine qopn2o ( fkg, tg )
!
! Description:
!    Calculation of optical depth of nitrogene oxide (N2O) absorption in bands
!    10 ( 1400 - 1250 cm**-1 ) and 11 ( 1250 - 1100 cm**-1 )
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      10/2016   Declarations and implicit none   <Josue Bock>
!          07/2016   Header including "USE ... ONLY"
!
! 1.0       ?        Original code.                    <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision, intent(in) :: fkg(nrlev)
! Array arguments with intent(out):
      double precision, intent(out) :: tg(nrlay)

! Local scalars:
      integer i,j

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

!- End of header ---------------------------------------------------------------

      do i=1,nrlay
         j=nrlev+1-i
         tg(i)=( fkg(i)+fkg(i+1) ) * ( p(j-1)-p(j) ) * 1.10459e-6
      end do

      end subroutine qopn2o
!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine qopo3i ( fkg, tg )
!
! Description:
!    Calculation of optical depth of ozone absorption in band
!    12 ( 1100 - 980 cm**-1 ).
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header including "USE ... ONLY"  <Josue Bock>
!                    Declarations and implicit none
!
! 1.0       ?        Original code.                   <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Array arguments with intent(in):
      double precision, intent(in) :: fkg(nrlev)
! Array arguments with intent(out):
      double precision, intent(out) :: tg(nrlay)

! Local scalars:
      integer i,j

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

      common /ozon/ qmo3(nrlev)
      double precision qmo3

!- End of header ---------------------------------------------------------------

      do i=1,nrlay
         j=nrlev+1-i
         tg(i)=( fkg(i) * qmo3(j) + fkg(i+1) * qmo3(j-1) )
     &        * ( p(j-1) - p(j) ) * 2.3808
      end do

      end subroutine qopo3i

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine qophc ( fkg, tg )
!
! Description:
!    Calculation of optical depth for overlapping absorption of water vapour and
!    carbon dioxide (CO2) in bands
!    14 ( 800 - 670 cm**-1) and 15 ( 670 - 540 cm**-1).
!    See page 86 of Fu (1991).
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header including "USE ... ONLY"  <Josue Bock>
!                    Declarations and implicit none
!
! 1.0       ?        Original code.                   <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Array arguments with intent(in):
      double precision, intent(in) :: fkg(nrlev)
! Array arguments with intent(out):
      double precision , intent(out) ::tg(nrlay)

! Local scalars:
      integer i,j

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

!- End of header ---------------------------------------------------------------

      do i=1,nrlay
         j=nrlev+1-i
         tg(i)=( fkg(i) + fkg(i+1) ) * ( p(j-1) - p(j) ) * 0.005
      end do
! See page 86 of Fu (1991).

      end subroutine qophc

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine tau(ib)
!
! Description:
!   Calculation of total optical depth {dtau}, single scattering albedo {om}
!   and Legendre coefficients {pl} for current spectral band {ib}
!   and cumulative probability {ig}.
!   The first index (second for {pl}) stands for cloud free (1) and
!   cloudy (2) parts.
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1       10/2016  Removed unused argument (ig)
!                    Removed NO2 in the first wavelength band
!                       ('tauno2' was set to 0 anyway)
!                       (see Loughlin et al., QJRMS 1997: the absorption of NO2 is not in PIFM2)

!           07/2016  Header including "USE ... ONLY"   <Josue Bock>
!
! 1.0       ?        Original code.                    <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      implicit none

! Subroutine arguments
! Array arguments with intent(in):
      integer, intent(in) :: ib

! Local scalars:
      integer iz, l
      double precision zsum1, zsum2
      double precision zx1_na
      double precision zf

! Common blocks:
      common /con/ tgcon(nrlay)
      double precision tgcon

      common /extra/ waer(nrlay),taer(nrlay),plaer(2,nrlay)
      double precision waer, taer, plaer

      common /gas/ tg(nrlay)
      double precision tg

      common /h2o/ t2w(nrlay), w2w(nrlay), pl2w(2,nrlay)
      double precision t2w, w2w, pl2w

      common /opohne/ dtau(2,nrlay),om(2,nrlay),pl(2,2,nrlay)
      double precision dtau, om, pl

      common /ray/ dtaur(nrlay), plr(2,nrlay)
      double precision dtaur, plr
!- End of header ---------------------------------------------------------------

! total optical depth {dtau}
      do iz=1,nrlay ! Index iz goes from top to bottom
         dtau(1,iz)=dtaur(iz)+taer(iz)+tgcon(iz)+tg(iz)
         dtau(2,iz)=dtau(1,iz)+t2w(iz)
      enddo

! total single scattering albedo {om}
      do iz=1,nrlay              ! Index iz goes from top to bottom
         zx1_na= taer(iz)*waer(iz)
         zsum1=dtaur(iz) + zx1_na
         zsum2=zsum1 + t2w(iz)*w2w(iz)
         if(dtau(1,iz) > 1.0e-20) then
            om(1,iz)=zsum1/dtau(1,iz)
            om(2,iz)=zsum2/dtau(2,iz)
         else
            om(1,iz)=0.0
            om(2,iz)=0.0
         endif
! Legendre-coefficients for total scattering {pl}
         do l=1,2
            zf=dtaur(iz)*plr(l,iz) + zx1_na*plaer(l,iz)
            if(zsum1 >= 1.d-20) then
               pl(l,1,iz)=zf/zsum1
               pl(l,2,iz)=(zf+t2w(iz)*w2w(iz)*pl2w(l,iz))/zsum2
            else
               pl(l,1,iz)=0.0
               pl(l,2,iz)=0.0
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
! Description:
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


!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header including "USE ... ONLY"                  <Josue Bock>
!          11/2016   ua ... wd changed from 1D aray (nrlay)  to scalars
!          11/2016   missing declarations, implicit none
!          11/2016   added u0red to avoid changing a global variable
!          11/2016   further bugfix related to u0red: previously, u0
!                     was changed within this SR, but not properly
!                     reset for all other cases. It has to be used
!                     ONLY in the resonance case!
!
! 1.0       ?        Original code.                                   <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     mbs,
     &     mbir,
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Array arguments with intent(in):
      integer, intent(in) :: ib                         ! current spectral band
      double precision, intent(in) :: u0                ! cosine (solar zenith angle)

! Local parameters:
      double precision, parameter :: u=2.d0             ! diffusivity factor = reciprocal of the mean effective cos(SZA)
      double precision, parameter :: delu0=0.001d0      ! resonance case correction
      double precision, parameter :: reson=0.1d-6       ! boundary value for resonance case
      double precision, parameter :: absfr=0.001d0      ! boundary value for no absorption case
      double precision, parameter :: strfr=0.03d0       ! boundary value for no scattering case
      double precision, parameter :: p1ray=0.1d0        ! boundary value for ray. scattering

! Local scalars:
      double precision ua, ub, uc, ud                   ! flux modification for clouds
      double precision va, vb, vc, vd
      double precision wa, wb, wc, wd
      double precision ak                               ! ratio of absorption and extinction
      double precision alph1, alph2, alph3, alph4       ! coefficients 'alpha'
      double precision u0kw                             !
      double precision u0red                            ! check of resonance case
      double precision f, emf, emfkw                    ! asymmetry factor of phase function
      double precision b0, bu0                          ! backscattering coefficient
      double precision dtu0, dtu, eps2, eps, omf, emomf ! equation factors
      double precision u02, ueps2, emu, e, e2, m, m2
      double precision e2m2, ouf, te, u0a1, u0a2
      double precision gam1, gam2, g1a1, da

      integer i, ip, l                                  ! loop indexes

! Local arrays:
      double precision a1(2,nrlay),a2(2,nrlay),a3(2,nrlay),a6(2,nrlay)     ! matrix coefficients (local)

! Common blocks:
      common /leck1/ a4(2,nrlay),a5(2,nrlay)                          ! matrix coefficients
      double precision a4, a5

      common /leck2/ sf(nrlev),sw(nrlev),ssf(nrlev),ssw(nrlev),             ! radiation fluxes
     &     f2f(nrlev),f2w(nrlev),f1f(nrlev),f1w(nrlev)
      double precision sf, sw, ssf, ssw, f2f, f2w, f1f, f1w

      common /opohne/ dtau(2,nrlay),om(2,nrlay),pl(2,2,nrlay)            ! optical variables
      double precision dtau, om, pl

      common /part/ cc(4,nrlay),bb(4,nrlay)                           ! cloudiness (continuity factors)
      double precision cc, bb

      common /tmp2/ as(mbs),ee(mbir)                            ! albedo and emissivity (unused here)
      double precision as, ee

!- End of header ---------------------------------------------------------------


! Definition of multiply used variables
      u0kw  = 1./u0

! Coefficients a for matrix solution of the radiative transfer equation.
!-----------------------------------------------------------------------
!    the coefficients a1 without Delta-Eddington are set to a6( *, *).

      do i=1,nrlay               ! Index i goes from top to bottom

         ! Loop for cloud free (index 1) and cloudy (index 2) parts
         !---------------------------------------------------------
         do l=1,2

! Case 1: no extinction
            if(dtau(l,i) <= 1.d-7) then
               a1(l,i) = 1.
               a2(l,i) = 0.
               a3(l,i) = 0.
               a4(l,i) = 1.
               a5(l,i) = 0.
               a6(l,i) = 1.

            else
               ! General definitions for all other cases
               dtu0    = dtau(l,i)*u0kw
               a6(l,i) = exp(-min(dtu0,7.5d1))
               dtu     = dtau(l,i)*u

! Case 2: no scattering (om < strfr)
               if(om(l,i) < strfr) then
                  a1(l,i)=a6(l,i)
                  a2(l,i)=0.
                  a3(l,i)=0.
                  a4(l,i)=exp(-min(dtu,7.5d1))
                  a5(l,i)=0.

               else
                  ! General definitions for all other cases
                  !
                  ! backscattering coefficients: beta-0 (here b0), beta(u0) (here bu0)
                  ! In the case of pure rayleigh scattering (p1 < p1ray) they
                  ! are set to 0.5.
                  ! p1 is the first Legendre-coefficient of the scattering function.
                  ! f is the asymmetry factor of Delta-Eddington
                  !
                  !   Anmerkung: Der Sprung von b0 und bu0 auf 0.5 im Fall p1 <
                  !   p1ray ist sicher nicht zufriedenstellend. Ist ein Provisorium. ! jjb NOK
                  ak    = 1.-om(l,i)
                  b0    = 0.5
                  bu0   = 0.5
                  f     = pl(2,l,i)/5.
                  emf   = 1.-f
                  emfkw = 1./emf
                  if(pl(1,l,i) >= p1ray) then
                     b0=(3.0-pl(1,l,i))/8.
                     bu0=0.5-u0/4.*(pl(1,l,i)-3.*f)*emfkw
                  end if

! Case 3: no absorption = pure scattering (ak < absfr, om = 1 )
                  if(ak < absfr) then
                     alph1=u*b0
                     alph3=bu0 ! jjb NOK 1-f missing
                     alph4=1.-alph3
                     gam1=alph3-alph1*u0*emfkw ! jjb NOK
                     a1(l,i)=exp(-min(dtu0*emf,7.5d1))
                     a4(l,i)=1./(1.+alph1*dtau(l,i))
                     a2(l,i)=a4(l,i)*(1.-gam1*(1.-a1(l,i)))-a1(l,i) ! jjb NOK
                     a3(l,i)=1.-a1(l,i)-a2(l,i) ! jjb NOK
                     a5(l,i)=1.-a4(l,i)

! Case 4: absorption and scattering (normal case)
                  else
                     alph2=u*b0*om(l,i)
                     alph1=u*ak+alph2
                     alph3=bu0*om(l,i) ! jjb NOK 1-f missing
                     alph4=om(l,i)-alph3 ! jjb NOK 1-f missing
                     eps2=alph1**2-alph2**2
                     eps=sqrt(eps2)
                     omf=om(l,i)*f
                     emomf=1.-omf

                     u02=u0**2
                     ueps2=u02*eps2
                     emu=emomf**2-ueps2

                     ! First initialisation of u0red, used in resonance case
                     u0red = u0

! Resonance case if the sun angle u0 is in resonance
! with the reference angle u:
! -> u0 is reduced by delu0
                     do while(abs(emu) <= reson)
                        u0red=u0red-delu0
                        if(u0red <= 0.) then
                           write(89,5)          ! jjb improve this by using globaly defined file unit
 5                         format(/1x,'Error in resonance-case'
     & //' correction: The sun is moved back below the horizon.'
     & //' Thus it is nighttime and shortwave fluxes are zero')
                           stop 'Error in subroutine kurzw.'
                        end if

                        u02=u0red**2
                        ueps2=u02*eps2
                        emu=emomf**2-ueps2
                     end do

                     a1(l,i) = exp(-min(dtu0*emomf,7.5d1))
                     e       = exp(-min(dtau(l,i)*eps,7.5d1))
                     m       = alph2/(alph1+eps)
                     e2      = e**2
                     m2      = m**2
                     e2m2    = e2*m2
                     ouf     = 1./(1.-e2m2)
                     a4(l,i) = e*(1.-m2)*ouf
                     a5(l,i) = m*(1.-e2)*ouf
                     te      = emf/emu
                     u0a1    = u0red*alph1
                     u0a2    = u0red*alph2
                     gam1    =  (alph3*(emomf-u0a1)-u0a2*alph4)*te
                     gam2    = -(alph4*(emomf+u0a1)+u0a2*alph3)*te
                     g1a1    = gam1*a1(l,i)
                     da      = a1(l,i)-a4(l,i)
                     a2(l,i) =  gam2*da-a5(l,i)*g1a1
                     a3(l,i) = -gam2*a5(l,i)-a4(l,i)*g1a1+gam1
                  end if
               end if
            end if
         end do
      end do

!------------------------------------------------------------------------------

! Calculation of the parallel fluxes and the right hand side of the
! diffuse system A * F  =  R0.
!
      sf(1)  = u0
      sw(1)  = 0.d0
      ssf(1) = sf(1)
      ssw(1) = 0.d0
      ua     = bb(1,1) * ssf(1)
      ub     = ssf(1) - ua
      f2f(1) = 0.d0
      f2w(1) = 0.d0
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
      do i=2,nrlay
         ip=i+1

         ua = bb(1,i) * ssf(i)
         ub = ssf(i) - ua
         uc = bb(1,i) * sf(i)
         ud = sf(i) - uc
         va = cc(3,i) * ssw(i)
         vb = ssw(i) - va
         vc = cc(3,i) * sw(i)
         vd = sw(i) - vc
         wa = ua + va
         wb = ub + vb
         wc = uc + vc
         wd = ud + vd
         ! cloud free part
         ssf(ip) = a1(1,i) * wa
         f2f(ip) = a2(1,i) * wa
         f1f(i)  = a3(1,i) * wa
         sf(ip)  = a6(1,i) * wc
         ! cloudy part
         ssw(ip) = a1(2,i) * wb
         f2w(ip) = a2(2,i) * wb
         f1w(i)  = a3(2,i) * wb
         sw(ip)  = a6(2,i) * wd
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
! Description:
!   Solution of the radiative transfer equation for current spectral band
!   (ir -> ib = 7 ... 18) and cumulative probability {ig}.
!
!   The arrays f2f, f2w, f1f, f1w are the diffuse downward (2) and upward (1)
!   ir fluxes. ('f' and 'w' stand for cloud free and cloudy parts)
!
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header including "USE ... ONLY"          <Josue Bock>
!          11/2016   missing declarations and implicit none
!
! 1.0       ?        Original code.                           <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     mbs,
     &     mbir,
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Array arguments with intent(in):
      integer, intent(in) :: ib                         ! current spectral band

! Local parameters:
      double precision, parameter :: u=1.66d0           ! diffusivity factor = reciprocal of the mean effective cos(SZA)

! Local scalars:
      double precision agdb, ak, alph1, alph2, at
      double precision b0, db, dtu
      double precision e, eps, epstau, eq
      double precision ha, hb
      double precision rm, rmq, rn
      integer i, ip, ii0, l                             ! loop indexes

! Local arrays
      double precision a6(2,nrlay)                         ! matrix coefficients (local)

! Common blocks:
      common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2(nrlay),
     &              frac(nrlay),ts,ntypa(nrlay),ntypd(nrlay)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

      common /leck1/ a4(2,nrlay),a5(2,nrlay)                          ! matrix coefficients
      double precision a4, a5

      common /leck2/ sf(nrlev),sw(nrlev),ssf(nrlev),ssw(nrlev),             ! radiation fluxes
     &     f2f(nrlev),f2w(nrlev),f1f(nrlev),f1w(nrlev)
      double precision sf, sw, ssf, ssw, f2f, f2w, f1f, f1w

      common /opohne/ dtau(2,nrlay),om(2,nrlay),pl(2,2,nrlay)            ! optical variables
      double precision dtau, om, pl

      common /part/ cc(4,nrlay),bb(4,nrlay)                           ! cloudiness (continuity factors)
      double precision cc, bb

      common /planci/ pib(nrlev),pibs                              ! black body radiation
      double precision pib, pibs

      common /tmp2/ as(mbs),ee(mbir)                            ! albedo (unused here) and emissivity
      double precision as, ee

!- End of header ---------------------------------------------------------------


      do i=1,nrlay               ! Index i goes from top to bottom

         ! Loop for cloud free (index 1) and cloudy parts (index 2)
         !---------------------------------------------------------
         do l=1,2

!
! Case 1: no extinction
!
            if(dtau(l,i) <= 1.d-7) then
               a4(l,i)=1.
               a5(l,i)=0.
               a6(l,i)=1.

!
! Case 2: no scattering (om < 1.0d-7)
!
            else if(om(l,i) <= 1.d-7) then
               dtu=dtau(l,i)*u
               a4(l,i)=exp(-dtu)
               a5(l,i)=0.
               a6(l,i)=(1.-a4(l,i))/dtu

!
! Case 3: no absorption (ak < 1.0d-7)
!
            else
               ak=1.-om(l,i)
               b0=(3.-pl(1,l,i))/8.
               alph1=u*(1.-(1.-b0)*om(l,i))
               alph2=u*b0*om(l,i)

               if(ak <= 1.d-7) then
                  at=alph1*dtau(l,i)
                  a4(l,i)=1./(1.+at)
                  a5(l,i)=a4(l,i)*at
                  a6(l,i)=0.

!
! Case 4: absorption and scattering (normal case)
!
               else
                  eps=sqrt(alph1**2-alph2**2)
                  epstau=eps*dtau(l,i)
                  e=0.
                  if(epstau < 87.) e=exp(-epstau)
                  rm=alph2/(alph1+eps)
                  eq=e**2
                  rmq=rm**2
                  rn=1.-eq*rmq
                  a4(l,i)=e*(1.-rmq)/rn
                  a5(l,i)=rm*(1.-eq)/rn
                  if(alph1+alph2 /= 0.d0) then ! jjb see [Z82] p. 219: indetermination case
                     a6(l,i)=(1.-a4(l,i)-a5(l,i)) /
     &                        ((alph1+alph2)*dtau(l,i))
                  else
                     a6(l,i)=1.d0
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
      f2f(1)=pib(1)
      f2w(1)=0.
!
! surface
      ii0=ib-mbs
      agdb=ee(ii0)*(pib(nrlev)-pibs)+(1.-ee(ii0))*
     & (1.-ee(ii0))*(pib(nrlev)-pib(nrlay))*a6(1,nrlay)*
     & (1.-frac(nrlay))
      f1w(nrlev)=agdb*frac(nrlay)
      f1f(nrlev)=agdb-f1w(nrlev)
!
! remaining layers
      do i=1,nrlay
         ip=i+1
         db=pib(i)-pib(ip)
         f1f(i)=(1.-frac(i))*a6(1,i)*db
         f1w(i)=frac(i)*a6(2,i)*db
         f2f(ip)=-f1f(i)
         f2w(ip)=-f1w(i)
      end do
!
! upper boundary condition for the right hand side of the third to sixth
! equation
      ha=bb(1,1)*f2f(1)
      hb=f2f(1)-ha
      f2f(2)=f2f(2)+a4(1,1)*ha
      f1f(1)=f1f(1)+a5(1,1)*ha
      f2w(2)=f2w(2)+a4(2,1)*hb
      f1w(1)=f1w(1)+a5(2,1)*hb

      end subroutine langw

!
! ---------------------------------------------------------------------
! *********************************************************************
! ---------------------------------------------------------------------
!
      subroutine jeanfr(ib)
!
! Description:
!   Solution of the linear equation system A(diffus) * F  =  R0.
!   The right hand side is calculated in KURZW or LANGW.
!
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header including "USE ... ONLY"  <Josue Bock>
!                    Removal of labeled do-loops
!
! 1.0       ?        Original code.                   <unknown>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     mbs,
     &     mbir,
     &     nrlay,
     &     nrlev

      implicit none

! Subroutine arguments
! Array arguments with intent(in):
      integer, intent(in) :: ib

! Local scalars:
      double precision ae
      double precision fa
      double precision ga,gb,gc,gd,ge,gf
      double precision ha,hb,hc,hd
      double precision tds1,tds2,tds3,tus1
      integer i,im,ip                       ! Loop indexes

! Local arrays:
      double precision tu(9,nrlay),td(7)

! Common blocks:
      common /leck1/ a4(2,nrlay),a5(2,nrlay)                          ! matrix coefficients
      double precision a4, a5

      common /leck2/ sf(nrlev),sw(nrlev),ssf(nrlev),ssw(nrlev),             ! radiation fluxes
     &     f2f(nrlev),f2w(nrlev),f1f(nrlev),f1w(nrlev)
      double precision sf, sw, ssf, ssw, f2f, f2w, f1f, f1w

      common /part/ cc(4,nrlay),bb(4,nrlay)
      double precision cc, bb

      common /tmp2/ as(mbs),ee(mbir)
      double precision as, ee
!- End of header ---------------------------------------------------------------

!   Diffusive System A * F = R0
!
! The right side is provided by subroutines KURZW or LANGW
!
! equation 3. to 6.: save of matix elements above the main diagonal
! in arrays -tu(k,i)
!
      tu(1,1)=0.d0
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
! above the main diagonal are saved in -tu(k, i).
! Right hand side saved in f2f, f2w, f1f, f1w.
!
      do i=2,nrlay
         im=i-1
         ip=i+1

         ga=bb(1,i)*tu(6,im)
         gb=tu(6,im)-ga
         gc=cc(3,i)*tu(8,im)
         gd=tu(8,im)-gc
         ha=ga+gc
         hc=gb+gd
         ga=bb(1,i)*tu(7,im)
         gb=tu(7,im)-ga
         gc=cc(3,i)*tu(9,im)
         gd=tu(9,im)-gc
         hb=ga+gc
         hd=gb+gd
         ga=bb(1,i)*f2f(i)
         ge=f2f(i)-ga
         gc=cc(3,i)*f2w(i)
         gf=f2w(i)-gc
         gb=ga+gc
         gd=ge+gf
         td(1)=1.d0/(1.d0-a5(1,i)*ha)
         f1f(i)=td(1)*(f1f(i)+a5(1,i)*gb)
         tu(1,i)=td(1)*a5(1,i)*hb
         fa=td(1)*a4(1,i)
         tu(2,i)=fa*bb(2,i)
         tu(3,i)=fa*cc(4,i)
         td(2)=a5(2,i)*hc
         td(3)=1.d0/(1.d0-a5(2,i)*hd-td(2)*tu(1,i))
         f1w(i)=td(3)*(f1w(i)+a5(2,i)*gd+td(2)*f1f(i))
         td(4)=a4(1,i)*ha
         td(5)=a4(1,i)*hb+td(4)*tu(1,i)
         f2f(ip)=f2f(ip)+a4(1,i)*gb+td(4)*f1f(i)+td(5)*f1w(i)
         tu(4,i)=td(3)*(a4(2,i)*cc(2,i)+td(2)*tu(2,i))
         tu(5,i)=td(3)*(a4(2,i)*bb(4,i)+td(2)*tu(3,i))
         tu(6,i)=a5(1,i)*bb(2,i)+td(4)*tu(2,i)+td(5)*tu(4,i)
         tu(7,i)=a5(1,i)*cc(4,i)+td(4)*tu(3,i)+td(5)*tu(5,i)
         td(6)=a4(2,i)*hc
         td(7)=a4(2,i)*hd+td(6)*tu(1,i)
         f2w(ip)=f2w(ip)+a4(2,i)*gd+td(6)*f1f(i)+td(7)*f1w(i)
         tu(8,i)=a5(2,i)*cc(2,i)+td(6)*tu(2,i)+td(7)*tu(4,i)
         tu(9,i)=a5(2,i)*bb(4,i)+td(6)*tu(3,i)+td(7)*tu(5,i)
      end do

! Elimination and backsubstitution of the last two equations.
! Calculation with albedo (solar) and 1-emissivity (ir)
!
      if(ib <= mbs) then
         ae=as(ib)
      else
         ae=1.-ee(ib-mbs)
      endif
      tds1=1./(1.-ae*tu(6,nrlay))
      f1f(nrlev)=tds1*(f1f(nrlev)+ae*f2f(nrlev))
      tus1=tds1*ae*tu(7,nrlay)
      tds2=ae*tu(8,nrlay)
      tds3=1./(1.-ae*tu(9,nrlay)-tds2*tus1)
      f1w(nrlev)=tds3*(f1w(nrlev)+ae*f2w(nrlev)+tds2*f1f(nrlev))
      f1f(nrlev)=f1f(nrlev)+tus1*f1w(nrlev)

! Now we have a upper-triangle-matrix with elements -tu(k, i) or 0 or 1.
! Resubstitution with results on the arrays f2f, f2w, f1f, f1w
!
      do i=nrlay,1,-1
         ip=i+1
         f2w(ip)=f2w(ip)+tu(8,i)*f1f(ip)+tu(9,i)*f1w(ip)
         f2f(ip)=f2f(ip)+tu(6,i)*f1f(ip)+tu(7,i)*f1w(ip)
         f1w(i)=f1w(i)+tu(4,i)*f1f(ip)+tu(5,i)*f1w(ip)
         f1f(i)=f1f(i)+tu(2,i)*f1f(ip)+tu(3,i)*f1w(ip)+tu(1,i)*f1w(i)
      end do

      end subroutine jeanfr
