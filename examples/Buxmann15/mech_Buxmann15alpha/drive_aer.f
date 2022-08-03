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
     &     (tkpp,dt_ch,k,xcvv1,xcvv2,yhal,yiod,
     &      yliq1,yliq2,yhet1,yhet2,air,h2o,xph_rat)

      USE constants, ONLY :
! Imported Parameters:
     &     xconv1=>conv1 ! multiply by conv1 to get cm^3(air)/mlc --> m^3(air)/mol

      USE gas_common, ONLY :
     &     j1, j5,
     &     s1,s3,
     &     gas_k2m_a, gas_m2k_a,
     &     rad_k2m_a, rad_m2k_a

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j3,
     &     j6,
     &     n,
     &     nf,
     &     nkc,
     &     nlev,
     &     nrxn

      implicit none

      include 'KPP_ROOT_Parameters.h' ! KPP parameters
      include 'KPP_ROOT_Global.h'     ! KPP common blocs and additional user common blocks and other definitions

! Subroutine arguments
      double precision tkpp,dt_ch,xcvv1,xcvv2,yhal,yiod,
     &      yliq1,yliq2,yhet1,yhet2,air,h2o

      double precision xph_rat(nphrxn)
      integer k

! Local scalar:
      integer kl ! do loop index to search if k belongs to il(nlev) list
      integer j

! Common blocks:
      common /blck12/ cw(nkc,n),cm(nkc,n)
      double precision cw, cm

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      double precision sl1, sion1

      common /budg/ bg(2,nrxn,nlev),il(nlev)
      double precision bg ! reaction rates (bg(1,:,:): instantaneous, bg(2,:,:): cumulative)
      integer il          ! indexes of the selected levels for reaction rates output

      common /kpp_laer/ henry(NSPEC,nf),xkmt(NSPEC,nkc,nf),
     &     xkef(NSPEC,nkc,nf),xkeb(NSPEC,nkc,nf)
      double precision henry, xkmt,xkef,xkeb

      common /kpp_drya/ xkmtd(NSPEC,2,nf),xeq(NSPEC,nf)
      double precision xkmtd,xeq

!      common /k_surf/ xkmt_OHClm(nf,nkc) ! jjb gamma_surf now commented, but keep this!

c parameters for /kpp_rate_a/
      cvv1=xcvv1
      cvv2=xcvv2
      xliq1=yliq1
      xliq2=yliq2
      xhet1=yhet1
      xhet2=yhet2
      xhal=yhal
      xiod=yiod
      conv1=xconv1
      ph_rat=xph_rat ! jjb

c the following data is needed only for heterogeneous reactions on dry aerosol
      ycwd(:)=cw(:2,k)
      yxkmtd(:,:)=xkmtd(:,:,k)
      yxeq(:)    =xeq(:,k)

c liquid phase rates
      if (xliq1.eq.1..or.xliq2.eq.1.) then ! jjb this test is probably useless
         yhenry(:)=henry(:,k)
      else
         yhenry(:)=0.d0
      endif
      if (xliq1.eq.1.) then
         ycw(1)=cw(1,k)
         yxkmt(:,1)=xkmt(:,1,k)
         ykef(:,1)=xkef(:,1,k)
         ykeb(:,1)=xkeb(:,1,k)
!         ykmt_OHClm(1) = xkmt_OHClm(k,1) ! jjb gamma_surf now commented, but keep this!
      else
         ycw(1)=0.
         yxkmt(:,1)=0.
         ykef(:,1)=0.
         ykeb(:,1)=0.
!         ykmt_OHClm(1) = 0. ! jjb gamma_surf now commented, but keep this!
      endif
      if (xliq2.eq.1.) then
         ycw(2)=cw(2,k)
         yxkmt(:,2)=xkmt(:,2,k)
         ykef(:,2)=xkef(:,2,k)
         ykeb(:,2)=xkeb(:,2,k)
!         ykmt_OHClm(2) = xkmt_OHClm(k,2) ! jjb gamma_surf now commented, but keep this!
      else
         ycw(2)=0.
            yxkmt(:,2)=0.
            ykef(:,2)=0.
            ykeb(:,2)=0.
!            ykmt_OHClm(2) = 0. ! jjb gamma_surf now commented, but keep this!
      endif


c concentrations are handed over HERE (and not in separate SRs) because the
c parameter (I_XXX) are different for each KPP block



c PRN2,PRPN,OZID are products, that don't react further, so no transport is needed
c maybe they are interesting as output ?? #

c include C(ind_)=s1/3(k,)
      do j=1,j1
         C(gas_m2k_a(1,j)) = s1(gas_m2k_a(2,j),k)
      end do

      do j=1,j5
         C(rad_m2k_a(1,j)) = s3(rad_m2k_a(2,j),k)
      end do

C#DEFFIX
      FIX(indf_O2)  = 0.21*air
      FIX(indf_H2O) = h2o
      FIX(indf_N2) = 0.79*air

c liquid phase
      sl1(:,:,k)=max(0.d0,sl1(:,:,k)) ! eliminate negative values
      sion1(:,:,k)=max(0.d0,sion1(:,:,k)) ! eliminate negative values

      if (cvv1.gt.0) then
         FIX(indf_H2Ol1)=55.55/cvv1
      else
         FIX(indf_H2Ol1)=0.
      endif
      if (cvv2.gt.0) then
         FIX(indf_H2Ol2)=55.55/cvv2
      else
         FIX(indf_H2Ol2)=0.
      endif

c aerosol
c include C(ind_)=sl1/sion1(k,,1/2)
      include 'aer_mk.dat'


c integrate

      dt=dt_ch
      call Update_RCONST ()
      call INTEGRATE (tkpp,tkpp+dt_ch)
!        call bud_a (h2o,co,air,dt_ch,k) ! jjb h2o, co, air unused
!         call bud_a (dt_ch,k)

! Call budget subroutine only for selected levels
      do kl=1,nlev
         if(k.eq.il(kl)) then
            call bud_KPP_ROOT (dt_ch,kl)
            exit
         end if
      end do

! Call specific budget subroutine for all levels
      call bud_s_KPP_ROOT (dt_ch,k)

c hand-over concentrations: KPP --> MISTRA
c#DEFVAR
      do j=1,j1
         s1(j,k) = C(gas_k2m_a(j))
      end do
C#DEFRAD
      do j=1,j5
         s3(j,k) = C(rad_k2m_a(j))
      end do

c aerosol
c include sl1/sion1(k,,1/2)=C(ind_)
      include 'aer_km.dat'

      end subroutine KPP_ROOT_drive
