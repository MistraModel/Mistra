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
     &     (tkpp,dt_ch,k,xcvv1,xcvv2,xcvv3,xcvv4,yhal,yiod,
     &      yliq1,yliq2,yliq3,yliq4,yhet1,yhet2,air,h2o,xph_rat)


      USE constants, ONLY :
! Imported Parameters:
     &     xconv1=>conv1 ! multiply by conv1 to get cm^3(air)/mlc --> m^3(air)/mol

      USE gas_common, ONLY :
     &     j1, j5,
     &     s1, s3,
     &     gas_k2m_t, gas_m2k_t,
     &     rad_k2m_t, rad_m2k_t

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
      double precision tkpp,dt_ch,xcvv1,xcvv2,xcvv3,xcvv4,yhal,yiod,
     &      yliq1,yliq2,yliq3,yliq4,yhet1,yhet2,air,h2o

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

      common /kpp_ltot/ henry(NSPEC,nf),xkmt(NSPEC,nkc,nf),
     &     xkef(NSPEC,nkc,nf),xkeb(NSPEC,nkc,nf)
      double precision henry, xkmt,xkef,xkeb
!      common /k_surf/ xkmt_OHClm(nf,nkc) ! jjb gamma_surf now commented, but keep this!

      common /kpp_dryt/ xkmtd(NSPEC,2,nf),xeq(NSPEC,nf) ! jjb copied from drive_aer.f
      double precision xkmtd,xeq
!      common /kpp_dryp/ rcd(n,2),cwd(n,2) ! jjb copied from drive_aer.f

c parameters for /kpp_rate_t/
      cvv1=xcvv1
      cvv2=xcvv2
      cvv3=xcvv3
      cvv4=xcvv4
      xliq1=yliq1
      xliq2=yliq2
      xliq3=yliq3
      xliq4=yliq4
      xhet1=yhet1
      xhet2=yhet2
      xhal=yhal
      xiod=yiod
      conv1=xconv1
      ph_rat=xph_rat ! jjb

c the following data is needed only for heterogeneous reactions on dry aerosol ! jjb copied from drive_aer.f
      ycwd(:)=cw(:2,k)
      yxkmtd(:,:)=xkmtd(:,:,k)
      yxeq(:)    =xeq(:,k)

c liquid phase rates
 ! jjb this test is probably useless
      if (xliq1.eq.1..or.xliq2.eq.1..or.xliq3.eq.1..or.xliq4.eq.1.) then
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
      if (xliq3.eq.1.) then
         ycw(3)=cw(3,k)
         yxkmt(:,3)=xkmt(:,3,k)
         ykef(:,3)=xkef(:,3,k)
         ykeb(:,3)=xkeb(:,3,k)
!         ykmt_OHClm(3) = xkmt_OHClm(k,3) ! jjb gamma_surf now commented, but keep this!
      else
         ycw(3)=0.
         yxkmt(:,3)=0.
         ykef(:,3)=0.
         ykeb(:,3)=0.
!         ykmt_OHClm(3) = 0. ! jjb gamma_surf now commented, but keep this!
      endif
      if (xliq4.eq.1.) then
         ycw(4)=cw(4,k)
         yxkmt(:,4)=xkmt(:,4,k)
         ykef(:,4)=xkef(:,4,k)
         ykeb(:,4)=xkeb(:,4,k)
!         ykmt_OHClm(4) = xkmt_OHClm(k,4) ! jjb gamma_surf now commented, but keep this!
      else
         ycw(4)=0.
         yxkmt(:,4)=0.
         ykef(:,4)=0.
         ykeb(:,4)=0.
!         ykmt_OHClm(4) = 0. ! jjb gamma_surf now commented, but keep this!
      endif


c concentrations are handed over HERE (and not in seperate SRs) because the 
c parameter (I_XXX) are different for each KPP block



c PRN2,PRPN,OZID are products, that don't react further, so no transport is needed
c maybe they are interesting as output ?? #

! Transfer Mistra concentration arrays towards KPP arrays
      do j=1,j1
         C(gas_m2k_t(1,j)) = s1(gas_m2k_t(2,j),k)
      end do

      do j=1,j5
         C(rad_m2k_t(1,j)) = s3(rad_m2k_t(2,j),k)
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
         if (cvv3.gt.0) then 
            FIX(indf_H2Ol3)=55.55/cvv3
         else
            FIX(indf_H2Ol3)=0.
         endif
         if (cvv4.gt.0) then 
            FIX(indf_H2Ol4)=55.55/cvv4
         else
            FIX(indf_H2Ol4)=0.
         endif

c aerosol
c include C(ind_)=sl1/sion1(,1/2,k)
         include 'aer_mk.dat'

c liquid phase
c _l3: small droplets
         C(ind_NOl3) = sl1(1,3,k)
         C(ind_NO2l3) = sl1(2,3,k)
         C(ind_HNO3l3) = sl1(3,3,k)
         C(ind_NH3l3) = sl1(4,3,k)
         C(ind_SO2l3) = sl1(5,3,k)
         C(ind_SO4l3) = sl1(6,3,k)
         C(ind_O3l3) = sl1(7,3,k)
c         C(ind_CH4l3) = sl1(8,3,k)
c         C(ind_C2H6l3) = sl1(9,3,k)
c         C(ind_C3H8l3) = sl1(10,3,k)
c         C(ind_ALKAl3) = sl1(11,3,k)
c         C(ind_ETHEl3) = sl1(12,3,k)
c         C(ind_ALKEl3) = sl1(13,3,k)
c         C(ind_AROMl3) = sl1(14,3,k)
         C(ind_HCOOHl3) = sl1(15,3,k)  ! HCOOH is ACO2 in gas phase
         C(ind_ACTAl3) = sl1(16,3,k)
         C(ind_HCHOl3) = sl1(17,3,k)
c         C(ind_ALD2l3) = sl1(18,3,k)
         C(ind_H2O2l3) = sl1(19,3,k)
         C(ind_CH3OOHl3) = sl1(20,3,k)  ! CH3OOH is ROOH in gas phase
         C(ind_HONOl3) = sl1(21,3,k)
c         C(ind_PANl3) = sl1(22,3,k)
c         C(ind_TPANl3) = sl1(23,3,k)
c         C(ind_KETl3) = sl1(24,3,k)
c         C(ind_CRESl3) = sl1(25,3,k)
c         C(ind_DIALl3) = sl1(26,3,k)
c         C(ind_GLYXl3) = sl1(27,3,k)
c         C(ind_MGLYl3) = sl1(28,3,k)
c         C(ind_NH4NO3l3) = sl1(29,3,k)
         C(ind_HCll3) = sl1(30,3,k)
c         C(ind_R3N2l3) = sl1(31,3,k)
c         C(ind_RAN2l3) = sl1(32,3,k)
c         C(ind_RAN1l3) = sl1(33,3,k)
c         C(ind_N2O5l3) = sl1(34,3,k)
         C(ind_HNO4l3) = sl1(35,3,k)
         C(ind_NO3l3) = sl1(36,3,k)
         C(ind_DMSl3) = sl1(37,3,k)
         C(ind_HOCll3) = sl1(38,3,k)
c         C(ind_ClNO2l3) = sl1(39,3,k)
c         C(ind_ClNO3l3) = sl1(40,3,k)
         C(ind_Cl2l3) = sl1(41,3,k)
         C(ind_HBrl3) = sl1(42,3,k)
         C(ind_HOBrl3) = sl1(43,3,k)
c         C(ind_BrNO2l3) = sl1(44,3,k)
c         C(ind_BrNO3l3) = sl1(45,3,k)
         C(ind_Br2l3) = sl1(46,3,k)
         C(ind_BrCll3) = sl1(47,3,k)
cc         C(ind_HIl3) = sl1(48,3,k)
         C(ind_HOIl3) = sl1(49,3,k)
cc         C(ind_I2O2l3) = sl1(50,3,k)
cc         C(ind_INO2l3) = sl1(51,3,k)
cc         C(ind_INO3l3) = sl1(52,3,k)
         C(ind_I2l3) = sl1(53,3,k)
         C(ind_ICll3) = sl1(54,3,k)
         C(ind_IBrl3) = sl1(55,3,k)
cc         C(ind_CH3Il3) = sl1(56,3,k)
cc         C(ind_CH2I2l3) = sl1(57,3,k)
cc         C(ind_CH2ClIl3) = sl1(58,3,k)
cc         C(ind_C3H7Il3) = sl1(59,3,k)
         C(ind_DMSOl3) = sl1(60,3,k)
c         C(ind_CH3SO2l3) = sl1(61,3,k)
c         C(ind_CH3SO3l3) = sl1(62,3,k)
c         C(ind_CH3SO3Hl3) = sl1(63,3,k)
c         C(ind_COl3) = sl1(64,3,k)
c         C(ind_Cl2O2l3) = sl1(65,3,k)
c         C(ind_DMOOl3) = sl1(66,3,k)
c         C(ind_CH3Sl3) = sl1(67,3,k)
c         C(ind_CH3SOl3) = sl1(68,3,k)
c         C(ind_CH3SO2Hl3) = sl1(69,3,k)
         C(ind_DMSO2l3) = sl1(70,3,k)
cc         C(ind_CH2BrIl3) = sl1(71,3,k)  
cc         C(ind_CHBr2Il3) = sl1(72,3,k)  
cc         C(ind_C2H5Il3) = sl1(73,3,k)  
!         C(ind_HIO3l3) = sl1(74,3,k)
c         C(ind_NUCVl3) = sl1(75,3,k)
c         C(ind_SO3l3) = sl1(76,3,k)
c         C(ind_HOSO2l3) = sl1(77,3,k)
         C(ind_CO2l3) = sl1(78,3,k)
c         C(ind_I2Ol3) = sl1(79,3,k)
c         C(ind_I2O3l3) = sl1(80,3,k)
c         C(ind_I2O4l3) = sl1(81,3,k)
c         C(ind_I2O5l3) = sl1(82,3,k)
c         C(ind_INOl3) = sl1(83,3,k)
c         C(ind_Br2Ol3) = sl1(84,3,k)
c         C(ind_ClONOl3) = sl1(85,3,k)
c         C(ind_ClO3l3) = sl1(86,3,k)
c         C(ind_Cl2O3l3) = sl1(87,3,k)
         C(ind_CH3OHl3) = sl1(88,3,k)
         C(ind_C2H5OHl3) = sl1(89,3,k)
c         C(ind_H2l3) = sl1(90,3,k)
c         C(ind_NHSl3) = sl1(91,3,k)
c         C(ind_RCll3) = sl1(92,3,k)
c         C(ind_RBrl3) = sl1(93,3,k)
         C(ind_XORl3) = sl1(94,3,k)
         C(ind_SORl3) = sl1(95,3,k)
c         C(ind_SPANl3) = sl1(96,3,k)
c         C(ind_Hgl3) = sl1(97,3,k)
c         C(ind_HgOl3) = sl1(98,3,k)
c         C(ind_HgCll3) = sl1(99,3,k)
c         C(ind_HgCl2l3) = sl1(100,3,k)
c         C(ind_HgBrl3) = sl1(101,3,k)
c         C(ind_HgBr2l3) = sl1(102,3,k)

         C(ind_OHl3) = sl1(j2-j3+2,3,k)
         C(ind_HO2l3) = sl1(j2-j3+3,3,k)
         C(ind_DOMl3) = sl1(j2-j3+4,3,k)
         C(ind_HIO2l3) = sl1(j2-j3+5,3,k)
         C(ind_CH3OOl3) = sl1(j2-j3+6,3,k)  ! CH3OO is MO2 in gas phase
         C(ind_IOl3) = sl1(j2-j3+7,3,k)
         C(ind_Cll3) = sl1(j2-j3+8,3,k)
         C(ind_Brl3) = sl1(j2-j3+9,3,k)
c         C(ind_) = sl1(j2-j3+10,3,k)
         C(ind_O2l3) = sl1(j2-j3+11,3,k)
c         C(ind_OIOl3) = sl1(j2-j3+12,3,k)
c         C(ind_HgOH2l3) = sl1(j2-j3+13,3,k) 
c         C(ind_HgOHCll3) = sl1(j2-j3+14,3,k)
c         C(ind_HgSO3l3) = sl1(j2-j3+15,3,k)
c         C(ind_HgOHBrl3) = sl1(j2-j3+16,3,k)

c ions
         C(ind_Hpl3) =sion1(1,3,k)  
         C(ind_NH4pl3) =sion1(2,3,k)  
         C(ind_OHml3) =sion1(3,3,k)  
         C(ind_CH2OHSO3ml3) =sion1(4,3,k)  
         C(ind_HSO3ml3) =sion1(5,3,k)  
         C(ind_SO32ml3) =sion1(6,3,k)  
         C(ind_SO4ml3) =sion1(7,3,k) 
         C(ind_SO42ml3) =sion1(8,3,k)  
         C(ind_HCO3ml3) =sion1(9,3,k)  
         C(ind_CO3ml3) =sion1(10,3,k)  
         C(ind_O2ml3) =sion1(11,3,k)  
         C(ind_NO2ml3) =sion1(12,3,k)  
         C(ind_NO3ml3) =sion1(13,3,k)  
         C(ind_Clml3) =sion1(14,3,k)  
         C(ind_Cl2ml3) =sion1(15,3,k)  
         C(ind_HCOOml3) =sion1(16,3,k)  
c         C(ind_FE3pl3) =sion1(17,3,k)  
c         C(ind_MN2pl3) =sion1(18,3,k)  
         C(ind_HSO4ml3) =sion1(19,3,k)
c         C(ind_Napl3) =sion1(20,3,k)  
         C(ind_NO4ml3) = sion1(21,3,k)
         C(ind_ClOml3) = sion1(22,3,k)
         C(ind_ClOHml3) = sion1(23,3,k)
         C(ind_Brml3) = sion1(24,3,k)
         C(ind_Br2ml3) = sion1(25,3,k)
         C(ind_BrOml3) = sion1(26,3,k)
         C(ind_BrOHml3) = sion1(27,3,k)
         C(ind_BrCl2ml3) = sion1(28,3,k)
         C(ind_Br2Clml3) = sion1(29,3,k)
         C(ind_CH3SO3ml3) = sion1(30,3,k)
         C(ind_HSO5ml3) = sion1(31,3,k)
         C(ind_SO3ml3) = sion1(32,3,k)
         C(ind_SO5ml3) = sion1(33,3,k)
         C(ind_Iml3) = sion1(34,3,k)
         C(ind_IO2ml3) = sion1(35,3,k)
         C(ind_IO3ml3) = sion1(36,3,k)
         C(ind_ICl2ml3) = sion1(37,3,k)
         C(ind_IBr2ml3) = sion1(38,3,k)
         C(ind_CH3SO2ml3) = sion1(39,3,k)
c         C(ind_Hgpl3) = sion1(40,3,k)
c         C(ind_Hg2pl3) = sion1(41,3,k)
c         C(ind_HgOHpl3) = sion1(42,3,k)
c         C(ind_HgClpl3) = sion1(43,3,k)
c         C(ind_HgCl3ml3) = sion1(44,3,k)
c         C(ind_HgCl42ml3) = sion1(45,3,k)
c         C(ind_HgBrpl3) = sion1(46,3,k)
c         C(ind_HgBr3ml3) = sion1(47,3,k)
c         C(ind_HgBr42ml3) = sion1(48,3,k)
c         C(ind_HgSO322ml3) = sion1(49,3,k)
c         C(ind_l3) = sion1(50,3,k)

c liquid phase
c _l4: large droplets
         C(ind_NOl4) = sl1(1,4,k)
         C(ind_NO2l4) = sl1(2,4,k)
         C(ind_HNO3l4) = sl1(3,4,k)
         C(ind_NH3l4) = sl1(4,4,k)
         C(ind_SO2l4) = sl1(5,4,k)
         C(ind_SO4l4) = sl1(6,4,k)
         C(ind_O3l4) = sl1(7,4,k)
c         C(ind_CH4l4) = sl1(8,4,k)
c         C(ind_C2H6l4) = sl1(9,4,k)
c         C(ind_C3H8l4) = sl1(10,4,k)
c         C(ind_ALKAl4) = sl1(11,4,k)
c         C(ind_ETHEl4) = sl1(12,4,k)
c         C(ind_ALKEl4) = sl1(13,4,k)
c         C(ind_AROMl4) = sl1(14,4,k)
         C(ind_HCOOHl4) = sl1(15,4,k)  ! HCOOH is ACO2 in gas phase
         C(ind_ACTAl4) = sl1(16,4,k)
         C(ind_HCHOl4) = sl1(17,4,k)
c         C(ind_ALD2l4) = sl1(18,4,k)
         C(ind_H2O2l4) = sl1(19,4,k)
         C(ind_CH3OOHl4) = sl1(20,4,k)  ! CH3OOH is ROOH in gas phase
         C(ind_HONOl4) = sl1(21,4,k)
c         C(ind_PANl4) = sl1(22,4,k)
c         C(ind_TPANl4) = sl1(23,4,k)
c         C(ind_KETl4) = sl1(24,4,k)
c         C(ind_CRESl4) = sl1(25,4,k)
c         C(ind_DIALl4) = sl1(26,4,k)
c         C(ind_GLYXl4) = sl1(27,4,k)
c         C(ind_MGLYl4) = sl1(28,4,k)
c         C(ind_NH4NO3l4) = sl1(29,4,k)
         C(ind_HCll4) = sl1(30,4,k)
c         C(ind_R3N2l4) = sl1(31,4,k)
c         C(ind_RAN2l4) = sl1(32,4,k)
c         C(ind_RAN1l4) = sl1(33,4,k)
c         C(ind_N2O5l4) = sl1(34,4,k)
         C(ind_HNO4l4) = sl1(35,4,k)
         C(ind_NO3l4) = sl1(36,4,k)
         C(ind_DMSl4) = sl1(37,4,k)
         C(ind_HOCll4) = sl1(38,4,k)
c         C(ind_ClNO2l4) = sl1(39,4,k)
c         C(ind_ClNO3l4) = sl1(40,4,k)
         C(ind_Cl2l4) = sl1(41,4,k)
         C(ind_HBrl4) = sl1(42,4,k)
         C(ind_HOBrl4) = sl1(43,4,k)
c         C(ind_BrNO2l4) = sl1(44,4,k)
c         C(ind_BrNO3l4) = sl1(45,4,k)
         C(ind_Br2l4) = sl1(46,4,k)
         C(ind_BrCll4) = sl1(47,4,k)
cc         C(ind_HIl4) = sl1(48,4,k)
         C(ind_HOIl4) = sl1(49,4,k)
cc         C(ind_I2O2l4) = sl1(50,4,k)
cc         C(ind_INO2l4) = sl1(51,4,k)
cc         C(ind_INO3l4) = sl1(52,4,k)
         C(ind_I2l4) = sl1(53,4,k)
         C(ind_ICll4) = sl1(54,4,k)
         C(ind_IBrl4) = sl1(55,4,k)
cc         C(ind_CH3Il4) = sl1(56,4,k)
cc         C(ind_CH2I2l4) = sl1(57,4,k)
cc         C(ind_CH2ClIl4) = sl1(58,4,k)
cc         C(ind_C3H7Il4) = sl1(59,4,k)
         C(ind_DMSOl4) = sl1(60,4,k)
c         C(ind_CH3SO2l4) = sl1(61,4,k)
c         C(ind_CH3SO3l4) = sl1(62,4,k)
c         C(ind_CH3SO3Hl4) = sl1(63,4,k)
c         C(ind_COl4) = sl1(64,4,k)
c         C(ind_Cl2O2l4) = sl1(65,4,k)
c         C(ind_DMOOl4) = sl1(66,4,k)
c         C(ind_CH3Sl4) = sl1(67,4,k)
c         C(ind_CH3SOl4) = sl1(68,4,k)
c         C(ind_CH3SO2Hl4) = sl1(69,4,k)
         C(ind_DMSO2l4) = sl1(70,4,k)
cc         C(ind_CH2BrIl4) = sl1(71,4,k)  
cc         C(ind_CHBr2Il4) = sl1(72,4,k)  
cc         C(ind_C2H5Il4) = sl1(73,4,k) 
!         C(ind_HIO3l4) = sl1(74,4,k)
c         C(ind_NUCVl4) = sl1(75,4,k)
c         C(ind_SO3l4) = sl1(76,4,k)
c         C(ind_HOSO2l4) = sl1(77,4,k)
         C(ind_CO2l4) = sl1(78,4,k)
c         C(ind_I2Ol4) = sl1(79,4,k)
c         C(ind_I2O3l4) = sl1(80,4,k)
c         C(ind_I2O4l4) = sl1(81,4,k)
c         C(ind_I2O5l4) = sl1(82,4,k)
c         C(ind_INOl4) = sl1(83,4,k)
c         C(ind_Br2Ol4) = sl1(84,4,k)
c         C(ind_ClONOl4) = sl1(85,4,k)
c         C(ind_ClO3l4) = sl1(86,4,k)
c         C(ind_Cl2O3l4) = sl1(87,4,k)
         C(ind_CH3OHl4) = sl1(88,4,k)
         C(ind_C2H5OHl4) = sl1(89,4,k)
c         C(ind_H2l4) = sl1(90,4,k)
c         C(ind_NHSl4) = sl1(91,4,k)
c         C(ind_RCll4) = sl1(92,4,k)
c         C(ind_RBrl4) = sl1(93,4,k)
         C(ind_XORl4) = sl1(94,4,k)
         C(ind_SORl4) = sl1(95,4,k)
c         C(ind_SPANl4) = sl1(96,4,k)
c         C(ind_Hgl4) = sl1(97,4,k)
c         C(ind_HgOl4) = sl1(98,4,k)
c         C(ind_HgCll4) = sl1(99,4,k)
c         C(ind_HgCl2l4) = sl1(100,4,k)
c         C(ind_HgBrl4) = sl1(101,4,k)
c         C(ind_HgBr2l4) = sl1(102,4,k)

         C(ind_OHl4) = sl1(j2-j3+2,4,k)
         C(ind_HO2l4) = sl1(j2-j3+3,4,k)
         C(ind_DOMl4) = sl1(j2-j3+4,4,k)
         C(ind_HIO2l4) = sl1(j2-j3+5,4,k)
         C(ind_CH3OOl4) = sl1(j2-j3+6,4,k)  ! CH3OO is MO2 in gas phase
         C(ind_IOl4) = sl1(j2-j3+7,4,k)
         C(ind_Cll4) = sl1(j2-j3+8,4,k)
         C(ind_Brl4) = sl1(j2-j3+9,4,k)
c         C(ind_) = sl1(j2-j3+10,4,k)
         C(ind_O2l4) = sl1(j2-j3+11,4,k)
c         C(ind_OIOl4) = sl1(j2-j3+12,4,k)
c         C(ind_HgOH2l4) = sl1(j2-j3+13,4,k) 
c         C(ind_HgOHCll4) = sl1(j2-j3+14,4,k)
c         C(ind_HgSO3l4) = sl1(j2-j3+15,4,k)
c         C(ind_HgOHBrl4) = sl1(j2-j3+16,4,k)

c ions
         C(ind_Hpl4) =sion1(1,4,k)  
         C(ind_NH4pl4) =sion1(2,4,k)  
         C(ind_OHml4) =sion1(3,4,k)  
         C(ind_CH2OHSO3ml4) =sion1(4,4,k)  
         C(ind_HSO3ml4) =sion1(5,4,k)  
         C(ind_SO32ml4) =sion1(6,4,k)  
         C(ind_SO4ml4) =sion1(7,4,k) 
         C(ind_SO42ml4) =sion1(8,4,k)  
         C(ind_HCO3ml4) =sion1(9,4,k)  
         C(ind_CO3ml4) =sion1(10,4,k)  
         C(ind_O2ml4) =sion1(11,4,k)  
         C(ind_NO2ml4) =sion1(12,4,k)  
         C(ind_NO3ml4) =sion1(13,4,k)  
         C(ind_Clml4) =sion1(14,4,k)  
         C(ind_Cl2ml4) =sion1(15,4,k)  
         C(ind_HCOOml4) =sion1(16,4,k)  
c         C(ind_FE3pl4) =sion1(17,4,k)  
c         C(ind_MN2pl4) =sion1(18,4,k)  
         C(ind_HSO4ml4) =sion1(19,4,k)
c         C(ind_Napl4) =sion1(20,4,k)  
         C(ind_NO4ml4) = sion1(21,4,k)
         C(ind_ClOml4) = sion1(22,4,k)
         C(ind_ClOHml4) = sion1(23,4,k)
         C(ind_Brml4) = sion1(24,4,k)
         C(ind_Br2ml4) = sion1(25,4,k)
         C(ind_BrOml4) = sion1(26,4,k)
         C(ind_BrOHml4) = sion1(27,4,k)
         C(ind_BrCl2ml4) = sion1(28,4,k)
         C(ind_Br2Clml4) = sion1(29,4,k)
         C(ind_CH3SO3ml4) = sion1(30,4,k)
         C(ind_HSO5ml4) = sion1(31,4,k)
         C(ind_SO3ml4) = sion1(32,4,k)
         C(ind_SO5ml4) = sion1(33,4,k)
         C(ind_Iml4) = sion1(34,4,k)
         C(ind_IO2ml4) = sion1(35,4,k)
         C(ind_IO3ml4) = sion1(36,4,k)
         C(ind_ICl2ml4) = sion1(37,4,k)
         C(ind_IBr2ml4) = sion1(38,4,k)
         C(ind_CH3SO2ml4) = sion1(39,4,k)
c         C(ind_Hgpl4) = sion1(40,4,k)
c         C(ind_Hg2pl4) = sion1(41,4,k)
c         C(ind_HgOHpl4) = sion1(42,4,k)
c         C(ind_HgClpl4) = sion1(43,4,k)
c         C(ind_HgCl3ml4) = sion1(44,4,k)
c         C(ind_HgCl42ml4) = sion1(45,4,k)
c         C(ind_HgBrpl4) = sion1(46,4,k)
c         C(ind_HgBr3ml4) = sion1(47,4,k)
c         C(ind_HgBr42ml4) = sion1(48,4,k)
c         C(ind_HgSO322ml4) = sion1(49,4,k)
c         C(ind_l4) = sion1(50,4,k)

c integrate

         dt=dt_ch
         call Update_RCONST ()
         call INTEGRATE (tkpp,tkpp+dt_ch)
!        call bud_t (h2o,co,air,dt_ch,k) ! jjb h2o, co,air  unused
!         call bud_t (dt_ch,k)

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
         s1(j,k) = C(gas_k2m_t(j))
      end do
C#DEFRAD
      do j=1,j5
         s3(j,k) = C(rad_k2m_t(j))
      end do

c aerosol
c include sl1/sion1(,1/2,k)=C(ind_)
      include 'aer_km.dat'



c liquid phase
c _l3: small droplets
         sl1(1,3,k) = C(ind_NOl3)
         sl1(2,3,k) = C(ind_NO2l3)
         sl1(3,3,k) = C(ind_HNO3l3)
         sl1(4,3,k) = C(ind_NH3l3)
         sl1(5,3,k) = C(ind_SO2l3)
         sl1(6,3,k) = C(ind_SO4l3)
         sl1(7,3,k) = C(ind_O3l3)
c         sl1(8,3,k) = C(ind_CH4l3)
c         sl1(9,3,k) = C(ind_C2H6l3)
c         sl1(10,3,k) = C(ind_C3H8l3)
c         sl1(11,3,k) = C(ind_ALKAl3)
c         sl1(12,3,k) = C(ind_ETHEl3)
c         sl1(13,3,k) = C(ind_ALKEl3)
c         sl1(14,3,k) = C(ind_AROMl3)
         sl1(15,3,k) = C(ind_HCOOHl3)  ! HCOOH is ACO2 in gas phase
         sl1(16,3,k) = C(ind_ACTAl3)
         sl1(17,3,k) = C(ind_HCHOl3)
c         sl1(18,3,k) = C(ind_ALD2l3)
         sl1(19,3,k) = C(ind_H2O2l3)
         sl1(20,3,k) = C(ind_CH3OOHl3)  ! CH3OOH is ROOH in gas phase
         sl1(21,3,k) = C(ind_HONOl3)
c         sl1(22,3,k) = C(ind_PANl3)
c         sl1(23,3,k) = C(ind_TPANl3)
c         sl1(24,3,k) = C(ind_KETl3)
c         sl1(25,3,k) = C(ind_CRESl3)
c         sl1(26,3,k) = C(ind_DIALl3)
c         sl1(27,3,k) = C(ind_GLYXl3)
c         sl1(28,3,k) = C(ind_MGLYl3)
c         sl1(29,3,k) = C(ind_NH4NO3l3)
         sl1(30,3,k) = C(ind_HCll3)
c         sl1(31,3,k) = C(ind_R3N2l3)
c         sl1(32,3,k) = C(ind_RAN2l3)
c         sl1(33,3,k) = C(ind_RAN1l3)
c         sl1(34,3,k) = C(ind_N2O5l3)
         sl1(35,3,k) = C(ind_HNO4l3)
         sl1(36,3,k) = C(ind_NO3l3)
         sl1(37,3,k) = C(ind_DMSl3)
         sl1(38,3,k) = C(ind_HOCll3)
c         sl1(39,3,k) = C(ind_ClNO2l3)
c         sl1(40,3,k) = C(ind_ClNO3l3)
         sl1(41,3,k) = C(ind_Cl2l3)
         sl1(42,3,k) = C(ind_HBrl3)
         sl1(43,3,k) = C(ind_HOBrl3)
c         sl1(44,3,k) = C(ind_BrNO2l3)
c         sl1(45,3,k) = C(ind_BrNO3l3)
         sl1(46,3,k) = C(ind_Br2l3)
         sl1(47,3,k) = C(ind_BrCll3)
cc         sl1(48,3,k) = C(ind_HIl3)
         sl1(49,3,k) = C(ind_HOIl3)
cc         sl1(50,3,k) = C(ind_I2O2l3)
cc         sl1(51,3,k) = C(ind_INO2l3)
cc         sl1(52,3,k) = C(ind_INO3l3)
         sl1(53,3,k) = C(ind_I2l3)
         sl1(54,3,k) = C(ind_ICll3)
         sl1(55,3,k) = C(ind_IBrl3)
cc         sl1(56,3,k)= C(ind_CH3Il3)
cc         sl1(57,3,k)= C(ind_CH2I2l3)
cc         sl1(58,3,k)= C(ind_CH2ClIl3)
cc         sl1(59,3,k)= C(ind_C3H7Il3)
         sl1(60,3,k) = C(ind_DMSOl3)
c         sl1(61,3,k)= C(ind_CH3SO2l3)
c         sl1(62,3,k)= C(ind_CH3SO3l3)
c         sl1(63,3,k)= C(ind_CH3SO3Hl3)
c         sl1(64,3,k) = C(ind_COl3) 
c         sl1(65,3,k) = C(ind_Cl2O2l3)  
c         sl1(66,3,k) = C(ind_DMOOl3)  
c         sl1(67,3,k) = C(ind_CH3Sl3)  
c         sl1(68,3,k) = C(ind_CH3SOl3) 
c         sl1(69,3,k) = C(ind_CH3SO2Hl3)  
         sl1(70,3,k) = C(ind_DMSO2l3)
cc         sl1(71,3,k) = C(ind_CH2BrIl3)
cc         sl1(72,3,k) = C(ind_CHBr2Il3)
cc         sl1(73,3,k) = C(ind_C2H5Il3)
!         sl1(74,3,k) = C(ind_HIO3l3)
c         sl1(75,3,k) = C(ind_NUCVl3)
c         sl1(76,3,k) = C(ind_SO3l3)
c         sl1(77,3,k) = C(ind_HOSO2l3)
         sl1(78,3,k) = C(ind_CO2l3)
c         sl1(79,3,k) = C(ind_I2Ol3)
c         sl1(80,3,k) = C(ind_I2O3l3)
c         sl1(81,3,k) = C(ind_I2O4l3)
c         sl1(82,3,k) = C(ind_I2O5l3)
c         sl1(83,3,k) = C(ind_INOl3)
c         sl1(84,3,k) = C(ind_Br2Ol3)
c         sl1(85,3,k) = C(ind_ClONOl3)
c         sl1(86,3,k) = C(ind_ClO3l3)
c         sl1(87,3,k) = C(ind_Cl2O3l3)
         sl1(88,3,k) = C(ind_CH3OHl3)
         sl1(89,3,k) = C(ind_C2H5OHl3)
c         sl1(90,3,k) = C(ind_H2l3)
c         sl1(91,3,k) = C(ind_NHSl3)
c         sl1(92,3,k) = C(ind_RCll3)
c         sl1(93,3,k) = C(ind_RBrl3)
         sl1(94,3,k) = C(ind_XORl3)
         sl1(95,3,k) = C(ind_SORl3)
c         sl1(96,3,k) = C(ind_SPANl3)
c         sl1(97,3,k) = C(ind_Hgl3)  
c         sl1(98,3,k) = C(ind_HgOl3)  
c         sl1(99,3,k) = C(ind_HgCll3) 
c         sl1(100,3,k) = C(ind_HgCl2l3)  
c         sl1(101,3,k) = C(ind_HgBrl3)  
c         sl1(102,3,k) = C(ind_HgBr2l3) 

         sl1(j2-j3+2,3,k) = C(ind_OHl3)
         sl1(j2-j3+3,3,k) = C(ind_HO2l3)
         sl1(j2-j3+4,3,k) = C(ind_DOMl3)
         sl1(j2-j3+5,3,k) = C(ind_HIO2l3)
         sl1(j2-j3+6,3,k) = C(ind_CH3OOl3)  ! CH3OO is MO2 in gas phase
         sl1(j2-j3+7,3,k) = C(ind_IOl3)
         sl1(j2-j3+8,3,k) = C(ind_Cll3)
         sl1(j2-j3+9,3,k) = C(ind_Brl3)
c         sl1(j2-j3+10,3,k) = C(ind_)
         sl1(j2-j3+11,3,k) = C(ind_O2l3)
c         sl1(j2-j3+12,3,k) = C(ind_OIOl3)
c         sl1(j2-j3+13,3,k) = C(ind_HgOH2l3)   
c         sl1(j2-j3+14,3,k) = C(ind_HgOHCll3)  
c         sl1(j2-j3+15,3,k) = C(ind_HgSO3l3)  
c         sl1(j2-j3+16,3,k) = C(ind_HgOHBrl3) 

c ions
         sion1(1,3,k) =C(ind_Hpl3) 
         sion1(2,3,k) =C(ind_NH4pl3) 
         sion1(3,3,k) =C(ind_OHml3) 
         sion1(4,3,k) =C(ind_CH2OHSO3ml3) 
         sion1(5,3,k) =C(ind_HSO3ml3) 
         sion1(6,3,k) =C(ind_SO32ml3) 
         sion1(7,3,k) =C(ind_SO4ml3) 
         sion1(8,3,k) =C(ind_SO42ml3) 
         sion1(9,3,k) =C(ind_HCO3ml3) 
         sion1(10,3,k) =C(ind_CO3ml3) 
         sion1(11,3,k) =C(ind_O2ml3) 
         sion1(12,3,k) =C(ind_NO2ml3) 
         sion1(13,3,k) =C(ind_NO3ml3) 
         sion1(14,3,k) =C(ind_Clml3) 
         sion1(15,3,k) =C(ind_Cl2ml3) 
         sion1(16,3,k) =C(ind_HCOOml3) 
c         sion1(17,3,k) =C(ind_FE3pl3) 
c         sion1(18,3,k) =C(ind_MN2pl3)
         sion1(19,3,k) =C(ind_HSO4ml3) 
c         sion1(20,3,k) =C(ind_Napl3) 
         sion1(21,3,k) = C(ind_NO4ml3)
         sion1(22,3,k) = C(ind_ClOml3)
         sion1(23,3,k) = C(ind_ClOHml3)
         sion1(24,3,k) = C(ind_Brml3)
         sion1(25,3,k) = C(ind_Br2ml3)
         sion1(26,3,k) = C(ind_BrOml3)
         sion1(27,3,k) = C(ind_BrOHml3)
         sion1(28,3,k) = C(ind_BrCl2ml3)
         sion1(29,3,k) = C(ind_Br2Clml3)
         sion1(30,3,k) = C(ind_CH3SO3ml3)
         sion1(31,3,k) = C(ind_HSO5ml3)
         sion1(32,3,k) = C(ind_SO3ml3)
         sion1(33,3,k) = C(ind_SO5ml3)
         sion1(34,3,k) = C(ind_Iml3)
         sion1(35,3,k) = C(ind_IO2ml3)
         sion1(36,3,k) = C(ind_IO3ml3)
         sion1(37,3,k) = C(ind_ICl2ml3)
         sion1(38,3,k) = C(ind_IBr2ml3)
         sion1(39,3,k) = C(ind_CH3SO2ml3)
c         sion1(40,3,k) = C(ind_Hgpl3) 
c         sion1(41,3,k) = C(ind_Hg2pl3) 
c         sion1(42,3,k) = C(ind_HgOHpl3) 
c         sion1(43,3,k) = C(ind_HgClpl3) 
c         sion1(44,3,k) = C(ind_HgCl3ml3) 
c         sion1(45,3,k) = C(ind_HgCl42ml3) 
c         sion1(46,3,k) = C(ind_HgBrpl3) 
c         sion1(47,3,k) = C(ind_HgBr3ml3) 
c         sion1(48,3,k) = C(ind_HgBr42ml3) 
c         sion1(49,3,k) = C(ind_HgSO322ml3) 
c         sion1(50,3,k) = C(ind_l3) 

c liquid phase
c _l4: large droplets
         sl1(1,4,k) = C(ind_NOl4)
         sl1(2,4,k) = C(ind_NO2l4)
         sl1(3,4,k) = C(ind_HNO3l4)
         sl1(4,4,k) = C(ind_NH3l4)
         sl1(5,4,k) = C(ind_SO2l4)
         sl1(6,4,k) = C(ind_SO4l4)
         sl1(7,4,k) = C(ind_O3l4)
c         sl1(8,4,k) = C(ind_CH4l4)
c         sl1(9,4,k) = C(ind_C2H6l4)
c         sl1(10,4,k) = C(ind_C3H8l4)
c         sl1(11,4,k) = C(ind_ALKAl4)
c         sl1(12,4,k) = C(ind_ETHEl4)
c         sl1(13,4,k) = C(ind_ALKEl4)
c         sl1(14,4,k) = C(ind_AROMl4)
         sl1(15,4,k) = C(ind_HCOOHl4)  ! HCOOH is ACO2 in gas phase
         sl1(16,4,k) = C(ind_ACTAl4)
         sl1(17,4,k) = C(ind_HCHOl4)
c         sl1(18,4,k) = C(ind_ALD2l4)
         sl1(19,4,k) = C(ind_H2O2l4)
         sl1(20,4,k) = C(ind_CH3OOHl4)  ! CH3OOH is ROOH in gas phase
         sl1(21,4,k) = C(ind_HONOl4)
c         sl1(22,4,k) = C(ind_PANl4)
c         sl1(23,4,k) = C(ind_TPANl4)
c         sl1(24,4,k) = C(ind_KETl4)
c         sl1(25,4,k) = C(ind_CRESl4)
c         sl1(26,4,k) = C(ind_DIALl4)
c         sl1(27,4,k) = C(ind_GLYXl4)
c         sl1(28,4,k) = C(ind_MGLYl4)
c         sl1(29,4,k) = C(ind_NH4NO3l4)
         sl1(30,4,k) = C(ind_HCll4)
c         sl1(31,4,k) = C(ind_R3N2l4)
c         sl1(32,4,k) = C(ind_RAN2l4)
c         sl1(33,4,k) = C(ind_RAN1l4)
c         sl1(34,4,k) = C(ind_N2O5l4)
         sl1(35,4,k) = C(ind_HNO4l4)
         sl1(36,4,k) = C(ind_NO3l4)
         sl1(37,4,k) = C(ind_DMSl4)
         sl1(38,4,k) = C(ind_HOCll4)
c         sl1(39,4,k) = C(ind_ClNO2l4)
c         sl1(40,4,k) = C(ind_ClNO3l4)
         sl1(41,4,k) = C(ind_Cl2l4)
         sl1(42,4,k) = C(ind_HBrl4)
         sl1(43,4,k) = C(ind_HOBrl4)
c         sl1(44,4,k) = C(ind_BrNO2l4)
c         sl1(45,4,k) = C(ind_BrNO3l4)
         sl1(46,4,k) = C(ind_Br2l4)
         sl1(47,4,k) = C(ind_BrCll4)
cc         sl1(48,4,k) = C(ind_HIl4)
         sl1(49,4,k) = C(ind_HOIl4)
cc         sl1(50,4,k) = C(ind_I2O2l4)
cc         sl1(51,4,k) = C(ind_INO2l4)
cc         sl1(52,4,k) = C(ind_INO3l4)
         sl1(53,4,k) = C(ind_I2l4)
         sl1(54,4,k) = C(ind_ICll4)
         sl1(55,4,k) = C(ind_IBrl4)
cc         sl1(56,4,k)= C(ind_CH3Il4)
cc         sl1(57,4,k)= C(ind_CH2I2l4)
cc         sl1(58,4,k)= C(ind_CH2ClIl4)
cc         sl1(59,4,k)= C(ind_C3H7Il4)
         sl1(60,4,k) = C(ind_DMSOl4)
c         sl1(61,4,k)= C(ind_CH3SO2l4)
c         sl1(62,4,k)= C(ind_CH3SO3l4)
c         sl1(63,4,k)= C(ind_CH3SO3Hl4)
c         sl1(64,4,k) = C(ind_COl4) 
c         sl1(65,4,k) = C(ind_Cl2O2l4)  
c         sl1(66,4,k) = C(ind_DMOOl4)  
c         sl1(67,4,k) = C(ind_CH3Sl4)  
c         sl1(68,4,k) = C(ind_CH3SOl4) 
c         sl1(69,4,k) = C(ind_CH3SO2Hl4)  
         sl1(70,4,k) = C(ind_DMSO2l4)
cc         sl1(71,4,k) = C(ind_CH2BrIl4)
cc         sl1(72,4,k) = C(ind_CHBr2Il4)
cc         sl1(73,4,k) = C(ind_C2H5Il4)
!         sl1(74,4,k) = C(ind_HIO3l4)
c         sl1(75,4,k) = C(ind_NUCVl4)
c         sl1(76,4,k) = C(ind_SO3l4)
c         sl1(77,4,k) = C(ind_HOSO2l4)
         sl1(78,4,k) = C(ind_CO2l4)
c         sl1(79,4,k) = C(ind_I2Ol4)
c         sl1(80,4,k) = C(ind_I2O3l4)
c         sl1(81,4,k) = C(ind_I2O4l4)
c         sl1(82,4,k) = C(ind_I2O5l4)
c         sl1(83,4,k) = C(ind_INOl4)
c         sl1(84,4,k) = C(ind_Br2Ol4)
c         sl1(85,4,k) = C(ind_ClONOl4)
c         sl1(86,4,k) = C(ind_ClO3l4)
c         sl1(87,4,k) = C(ind_Cl2O3l4)
         sl1(88,4,k) = C(ind_CH3OHl4)
         sl1(89,4,k) = C(ind_C2H5OHl4)
c         sl1(90,4,k) = C(ind_H2l4)
c         sl1(91,4,k) = C(ind_NHSl4)
c         sl1(92,4,k) = C(ind_RCll4)
c         sl1(93,4,k) = C(ind_RBrl4)
         sl1(94,4,k) = C(ind_XORl4)
         sl1(95,4,k) = C(ind_SORl4)
c         sl1(96,4,k) = C(ind_SPANl4)
c         sl1(97,4,k) = C(ind_Hgl4)  
c         sl1(98,4,k) = C(ind_HgOl4)  
c         sl1(99,4,k) = C(ind_HgCll4) 
c         sl1(100,4,k) = C(ind_HgCl2l4)  
c         sl1(101,4,k) = C(ind_HgBrl4)  
c         sl1(102,4,k) = C(ind_HgBr2l4) 

         sl1(j2-j3+2,4,k) = C(ind_OHl4)
         sl1(j2-j3+3,4,k) = C(ind_HO2l4)
         sl1(j2-j3+4,4,k) = C(ind_DOMl4)
         sl1(j2-j3+5,4,k) = C(ind_HIO2l4)
         sl1(j2-j3+6,4,k) = C(ind_CH3OOl4)  ! CH3OO is MO2 in gas phase
         sl1(j2-j3+7,4,k) = C(ind_IOl4)
         sl1(j2-j3+8,4,k) = C(ind_Cll4)
         sl1(j2-j3+9,4,k) = C(ind_Brl4)
c         sl1(j2-j3+10,4,k) = C(ind_)
         sl1(j2-j3+11,4,k) = C(ind_O2l4)
c         sl1(j2-j3+12,4,k) = C(ind_OIOl4)
c         sl1(j2-j3+13,4,k) = C(ind_HgOH2l4)   
c         sl1(j2-j3+14,4,k) = C(ind_HgOHCll4)  
c         sl1(j2-j3+15,4,k) = C(ind_HgSO3l4)  
c         sl1(j2-j3+16,4,k) = C(ind_HgOHBrl4) 

c ions
         sion1(1,4,k) =C(ind_Hpl4) 
         sion1(2,4,k) =C(ind_NH4pl4) 
         sion1(3,4,k) =C(ind_OHml4) 
         sion1(4,4,k) =C(ind_CH2OHSO3ml4) 
         sion1(5,4,k) =C(ind_HSO3ml4) 
         sion1(6,4,k) =C(ind_SO32ml4) 
         sion1(7,4,k) =C(ind_SO4ml4) 
         sion1(8,4,k) =C(ind_SO42ml4) 
         sion1(9,4,k) =C(ind_HCO3ml4) 
         sion1(10,4,k) =C(ind_CO3ml4) 
         sion1(11,4,k) =C(ind_O2ml4) 
         sion1(12,4,k) =C(ind_NO2ml4) 
         sion1(13,4,k) =C(ind_NO3ml4) 
         sion1(14,4,k) =C(ind_Clml4) 
         sion1(15,4,k) =C(ind_Cl2ml4) 
         sion1(16,4,k) =C(ind_HCOOml4) 
c         sion1(17,4,k) =C(ind_FE3pl4) 
c         sion1(18,4,k) =C(ind_MN2pl4)
         sion1(19,4,k) =C(ind_HSO4ml4) 
c         sion1(20,4,k) =C(ind_Napl4) 
         sion1(21,4,k) = C(ind_NO4ml4)
         sion1(22,4,k) = C(ind_ClOml4)
         sion1(23,4,k) = C(ind_ClOHml4)
         sion1(24,4,k) = C(ind_Brml4)
         sion1(25,4,k) = C(ind_Br2ml4)
         sion1(26,4,k) = C(ind_BrOml4)
         sion1(27,4,k) = C(ind_BrOHml4)
         sion1(28,4,k) = C(ind_BrCl2ml4)
         sion1(29,4,k) = C(ind_Br2Clml4)
         sion1(30,4,k) = C(ind_CH3SO3ml4)
         sion1(31,4,k) = C(ind_HSO5ml4)
         sion1(32,4,k) = C(ind_SO3ml4)
         sion1(33,4,k) = C(ind_SO5ml4)
         sion1(34,4,k) = C(ind_Iml4)
         sion1(35,4,k) = C(ind_IO2ml4)
         sion1(36,4,k) = C(ind_IO3ml4)
         sion1(37,4,k) = C(ind_ICl2ml4)
         sion1(38,4,k) = C(ind_IBr2ml4)
         sion1(39,4,k) = C(ind_CH3SO2ml4)
c         sion1(40,4,k) = C(ind_Hgpl4) 
c         sion1(41,4,k) = C(ind_Hg2pl4) 
c         sion1(42,4,k) = C(ind_HgOHpl4) 
c         sion1(43,4,k) = C(ind_HgClpl4) 
c         sion1(44,4,k) = C(ind_HgCl3ml4) 
c         sion1(45,4,k) = C(ind_HgCl42ml4) 
c         sion1(46,4,k) = C(ind_HgBrpl4) 
c         sion1(47,4,k) = C(ind_HgBr3ml4) 
c         sion1(48,4,k) = C(ind_HgBr42ml4) 
c         sion1(49,4,k) = C(ind_HgSO322ml4) 
c         sion1(50,4,k) = C(ind_l4) 


      end subroutine KPP_ROOT_drive

