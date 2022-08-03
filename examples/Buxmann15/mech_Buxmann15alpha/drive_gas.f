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

      module kpp_KPP_ROOT_Parameters
      include 'KPP_ROOT_Parameters.h' ! KPP parameters
      end module kpp_KPP_ROOT_Parameters

      module kpp_KPP_ROOT_Global
      include 'KPP_ROOT_Parameters.h' ! KPP parameters
      include 'KPP_ROOT_Global.h'     ! KPP common blocs and additional user common blocks and other definitions
      end module kpp_KPP_ROOT_Global



      subroutine KPP_ROOT_drive
     &     (tkpp,dt_ch,k,yhal,yiod,yhet1,yhet2,air,h2o,xph_rat)


      USE constants, ONLY :
! Imported Parameters:
     &     xconv1=>conv1 ! multiply by conv1 to get cm^3(air)/mlc --> m^3(air)/mol

      USE gas_common, ONLY :
     &     j1, j5,
     &     s1, s3,
     &     gas_k2m_g, gas_m2k_g,
     &     rad_k2m_g, rad_m2k_g

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     n,
     &     nkc,
     &     nlev,
     &     nrxn

      implicit none

      include 'KPP_ROOT_Parameters.h' ! KPP parameters
      include 'KPP_ROOT_Global.h'     ! KPP common blocs and additional user common blocks and other definitions
! Subroutine arguments
      double precision tkpp,dt_ch,yhal,yiod,yhet1,yhet2,air,h2o

      double precision xph_rat(nphrxn)
      integer k

! Local scalar:
      integer kl ! do loop index to search if k belongs to il(nlev) list
      integer j

! Common blocks:
      common /blck12/ cwd(nkc,n),cm(nkc,n)
      double precision cwd, cm

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      double precision sl1, sion1

      common /budg/ bg(2,nrxn,nlev),il(nlev)
      double precision bg ! reaction rates (bg(1,:,:): instantaneous, bg(2,:,:): cumulative)
      integer il          ! indexes of the selected levels for reaction rates output

      common /kpp_dryg/ xkmtd(NSPEC,2,n),henry(NSPEC,n),xeq(NSPEC,n)
      double precision xkmtd, henry, xeq
!     common /kpp_dryp/ rcd(n,2),cwd(n,2)


! parameters for /kpp_rate_g/
      xhal=yhal
      xiod=yiod
      xhet1=yhet1
      xhet2=yhet2
      conv1=xconv1
      ph_rat=xph_rat ! jjb

c the following data is needed only for heterogeneous reactions on dry aerosol
      ycwd(:)=cwd(:2,k)
      yxkmtd(:,:)=xkmtd(:,:,k)
      yhenry(:)  =henry(:,k)
      yxeq(:)    =xeq(:,k)

c concentrations are handed over HERE (and not in seperate SRs) because the
c parameter (ind_XXX) are different for each KPP block


! Transfer Mistra concentration arrays towards KPP arrays
      do j=1,j1
         C(gas_m2k_g(1,j)) = s1(gas_m2k_g(2,j),k)
      end do

      do j=1,j5
         C(rad_m2k_g(1,j)) = s3(rad_m2k_g(2,j),k)
      end do

c PRN2,PRPN,OZID are products, that don't react further, so no transport is needed
c maybe they are interesting as output ?? #



! #DEFFIX
      FIX(indf_O2)  = 0.21*air
      FIX(indf_H2O) = h2o
      FIX(indf_N2) = 0.79*air

! define "heterogeneous species"
      C(ind_HNO3l1) = max(0.d0,sl1(3,1,k))
      C(ind_NH3l1)  = max(0.d0,sl1(4,1,k))
      C(ind_SO4l1)  = max(0.d0,sl1(6,1,k))
      C(ind_HNO3l2) = max(0.d0,sl1(3,2,k))
      C(ind_NH3l2)  = max(0.d0,sl1(4,2,k))
      C(ind_SO4l2)  = max(0.d0,sl1(6,2,k))

!         C(ind_OHml1) =sion1(3,1,k)
!         C(ind_NO3ml1) =sion1(13,1,k)
!         C(ind_CLml1) =sion1(14,1,k)
!         C(ind_Brml1) = sion1(24,1,k)

!         C(ind_OHml2) =sion1(3,2,k)
!         C(ind_NO3ml2) =sion1(13,2,k)
!         C(ind_CLml2) =sion1(14,2,k)
!         C(ind_Brml2) = sion1(24,2,k)


! integrate

      dt=dt_ch
      call Update_RCONST ()
      call INTEGRATE (tkpp, tkpp+dt_ch )

!     call bud_g (h2o,co,air,dt_ch,k) ! jjb h2o, co, air unused
!     call bud_g (dt_ch,k)

! Call budget subroutine only for selected levels
      do kl=1,nlev
         if(k.eq.il(kl)) then
            call bud_KPP_ROOT (dt_ch,kl)
            exit
         end if
      end do

! Call specific budget subroutine for all levels
      call bud_s_KPP_ROOT (dt_ch,k)

! hand-over concentrations: KPP --> MISTRA
! #DEFVAR
      do j=1,j1
         s1(j,k) = C(gas_k2m_g(j))
      end do
! #DEFRAD
      do j=1,j5
         s3(j,k) = C(rad_k2m_g(j))
      end do

! define "heterogeneous species"
      sl1(3,1,k) = max(0.d0,C(ind_HNO3l1))
      sl1(4,1,k) = max(0.d0,C(ind_NH3l1))
      sl1(6,1,k) = max(0.d0,C(ind_SO4l1))
      sl1(3,2,k) = max(0.d0,C(ind_HNO3l2))
      sl1(4,2,k) = max(0.d0,C(ind_NH3l2))
      sl1(6,2,k) = max(0.d0,C(ind_SO4l2))

!         sion1(3,1,k) =C(ind_OHml1)
!         sion1(13,1,k) =C(ind_NO3ml1)
!         sion1(14,1,k) =C(ind_CLml1)
!         sion1(24,1,k) = C(ind_Brml1)

!         sion1(3,2,k) =C(ind_OHml2)
!         sion1(13,2,k) =C(ind_NO3ml2)
!         sion1(14,2,k) =C(ind_CLml2)
!         sion1(24,2,k) = C(ind_Brml2)

      end subroutine KPP_ROOT_drive
