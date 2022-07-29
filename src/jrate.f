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


! Description of jrate.f (as of 15/08/2016)
! =========================================
!
!  This file contains all SR of Jochen Landgraf's BAND MODEL to
! calculate the photolysis rates. SR photol is the former MAIN PROGRAM
! Two SR are needed to link MISTRA and BAND:
!    - SR read_data
!    - SR optic_data

! Major changes done on this file:
! ================================
!     - removal of an index used for 2D models (MJ = latitude), unused in Mistra
!       This implied to re-index most arrays, to remove one dimension. And to remove all do J=1,MJ loops
!     - major cleaning all over this file: unused variables, arguments, parameters and common blocks
!     - initialisation subroutine to read data only once


! I/O files
! =========
!
!   input files are: - flux.dat        \
!                    - sig0900.dat      |- read in SR cross_init
!                    - cheb_coeff.dat  /
!                    - lookt0900.dat   --- read by SR lookup
!
!  output files are: - prof.out ---------- written by SR atm_out    (switch SW4)
!                    - check_four.out ---- written by SR check_four (switch SW13)
!                    - f4st.out ---------- written by SR flux_out   (switch SW12)

! Flow chart
! ==========
!
! SR photol_initialize
!   |------SR cross_init
!   |------SR lookup
!
! SR photol
!   |------SR read_data
!   |        |---------SR atm_out
!   |------SR column
!   |------SR cross_atm
!   |        |---------SR sr_o2_km
!   |        |           |--------FN chebev
!   |        |---------SR sr_no_af
!   |------SR optic_data
!   |------SR four_intf
!   |        |---------SR check_four
!   |        |---------SR QTFS
!   |                    |--------SR adjust
!   |                    |--------SR QCCFE
!   |                               |------SR coefft0
!   |                               |------SR coefft1
!   |                               |------SR coefft2
!   |                               |------SR coefft4
!   |                               |------SR coefftl
!   |                               |------SR QCFEL
!   |------SR flux_out
!   |------SR photo_cal


*******************************************************************
      subroutine photol_initialize

      ! jjb this SR has been set up so that data files are read only once during model initialization
      !     previously, cross_init was called before cross_atm
      !                     lookup was called before photo_cal
      implicit none

      CALL CROSS_INIT
      CALL LOOKUP

      end subroutine photol_initialize
*******************************************************************


*******************************************************************
      Subroutine photol

*****************************************************
*  RUN BAND MODEL TO CACULATE PHOTOLYSIS RATES      *
*****************************************************

! History
!   13/08/2016 jjb final changes and comments before release
!              resulting from cleaning / debuging work between 11/2015 and 08/2016
!
!   done in this SR:
!       - USE global_params to get grid parameters from Mistra
!       - unused variables commented, but not deleted in case of (unlikely) future use
!       - removal of unused COMMON
!       - update of subroutine calls after arguments removal
!       - explicit definition of SW* switches (were assumed as real, when implicit DP was used)
!       - reorganisation of SRs band, lookup and photo_cal:
!            the first removed, the last 2 directly called by this SR
!            see more comments below
!       - implicit none
!

      USE global_params, ONLY :
     &     n_m=>n,
     &     nphrxn,
     &     nrlay

      IMPLICIT NONE

      INTEGER MAXLAY,MAXWAV,NMOM,NW
      PARAMETER(
     $     MAXLAY=nrlay,         !maximal number of layers
c top layer in  MISTRA is already "infinity" in Jochens code "infinity" is layer "0"
     $     MAXWAV=176,          !maximal number of wavelength intervals
     $     NMOM =4,             !number of moments of phase functions
     $     NW=7)                !number of wavelengths intervals

      DOUBLE PRECISION
     $     TEMP(0:MAXLAY),   !temperature at model levels [K]
     $     PRESS(0:MAXLAY),  !pressure at model levels [hPa]
     $     RELO3(0:MAXLAY),  !O3 volumn mixing ratio [ppm]
     $     ALBEDO(MAXWAV),   !shortwave albedo
     $     U0,               !cosine of solar zenith angle
     $     SCALEO3           !total vert. ozone column for scaling

!      INTEGER
!     $     ITYPE(MAXLAY)         !aerosol type

      DOUBLE PRECISION
     $     V2(0:MAXLAY),    !O2 column density [part./cm^2]
     $     DV2(MAXLAY),     !diff. O2 column density [part./cm^2]
     $     V2S(0:MAXLAY),   !slant O2 column density [part./cm^2]
     $     DV2S(MAXLAY),    !diff. slant O2 column density [part./cm^2]
     $     V3(0:MAXLAY),    !O3 column density [part./cm^2]
     $     DV3(MAXLAY),     !diff. O3 column density [part./cm^2]
     $     V3S(0:MAXLAY),   !slant O3 column density [part./cm^2]
     $     DV3S(MAXLAY),    !diff. slant O2 column density [part./cm^2]
     $     DENS(0:MAXLAY)   !density of air molecules [part./cm^3]

c     actinic fluxes

      DOUBLE PRECISION
     $     FACT(0:MAXLAY,NW)  !act.   flux

      DOUBLE PRECISION
     $     CST_O3(0:MAXLAY,MAXWAV),      ! jjb 16/08/2016: reindexed for CPU efficiency
     $     CST_O2(MAXWAV,0:MAXLAY)

      DOUBLE PRECISION
     $     QYNO3n(0:MAXLAY),
     $     QYO1D(MAXWAV,0:MAXLAY)

      DOUBLE PRECISION
     $     TAUS_CLR(NW,MAXLAY),  !scattering optical depth for clear sky
     $     TAUA_CLR(NW,MAXLAY),  !absorption optical depth for clear sky
     $     TAUS_AER(NW,MAXLAY),  !scattering optical depth for aerosol
     $     TAUA_AER(NW,MAXLAY),  !absorption optical depth for aerosol
     $     TAUS_CLD(NW),         !scattering optical depth for cloudy sky
     $     TAUA_CLD(NW)          !absorption optical depth for cloudy sky

      DOUBLE PRECISION
     $     PRAY(NW,MAXLAY,0:NMOM), !Rayleigh phase fct. moments
     $     PCLD(NW,0:NMOM),        !phase fct. moments for clouds
     $     PAER(NW,MAXLAY,0:NMOM)  !phase fct. moments for aerosols

      DOUBLE PRECISION
     $     ALB(NW),                !shortwave albedo
     $     FLX(NW)                 !extraterrestic flux

c     photolysis rates

      DOUBLE PRECISION
     $     RJO2(0:MAXLAY),       RJO3P(0:MAXLAY),
     $     RJH2O2(0:MAXLAY),
     $     RJHNO3(0:MAXLAY),     RJN2O(0:MAXLAY),
     $     RJNO2(0:MAXLAY),      RJN2O5(0:MAXLAY),
     $     RJCOH2(0:MAXLAY),     RJCHOH(0:MAXLAY),
     $     RJClONO2(0:MAXLAY),
     $                              RJCFC11(0:MAXLAY),
     $     RJCFC12(0:MAXLAY),
     $                              RJO1D(0:MAXLAY),
     $     RJHOCl(0:MAXLAY),     RJCH3OOH(0:MAXLAY),
     $     RJHNO4(0:MAXLAY),     RJNO2O(0:MAXLAY),
     $     RJNOO2(0:MAXLAY),
     $     RJCl2O2(0:MAXLAY),    RJBrNO3(0:MAXLAY),
     $                              RJNO3n(0:MAXLAY),
     $                              RJBrCl_noT(0:MAXLAY),
     $     RJClNO2(0:MAXLAY),    RJBrNO2(0:MAXLAY),
     $     RJBr2(0:MAXLAY),      RJIO(0:MAXLAY),
     $     RJINO3(0:MAXLAY),     RJCH3I(0:MAXLAY),
     $     RJI2(0:MAXLAY),       RJICl(0:MAXLAY),
     $     RJIBr(0:MAXLAY),      RJC3H7I(0:MAXLAY),
     $     RJCH2ClI(0:MAXLAY),   RJCH2I2(0:MAXLAY),
     $     RJINO2(0:MAXLAY),     RJBrO_noT(0:MAXLAY),
     $     RJOClO_noT(0:MAXLAY), RJCl2_noT(0:MAXLAY),
     $     RJHOI_jen91(0:MAXLAY),RJHOBr(0:MAXLAY),
     $     RJHONO(0:MAXLAY),     RJNO2m(0:MAXLAY),
     $     RJdumm23(0:MAXLAY),
     $     RJdumm24(0:MAXLAY),   RJdumm25(0:MAXLAY),
     $     RJdumm26(0:MAXLAY)

      DOUBLE PRECISION
     $     H_O2(0:MAXLAY),
     $     H_O3(0:MAXLAY),
     $     H_O4(0:MAXLAY)

      INTEGER
     $     NWS(NW),   !specification of interval: 1 < NWS(L) < MAXWAV
     $     NBAN(7)    !specification of interval for band calculation

      DATA NBAN/15,  43,  49,  56,  67,  80,  122/

      COMMON/WL/WAVE(MAXWAV),  !wavelength in the middle of the interval [cm]
     $          DWAVE(MAXWAV)  !width of the wavelength intervals [cm]
      DOUBLE PRECISION WAVE,DWAVE

      INTEGER K,L ! indexes of do loops

c jrates for MISTRA
      common /band_rat/ photol_j(nphrxn,n_m)
      double precision photol_j

      INTEGER SW2, SW4, SW7, SW12, SW13 ! jjb added 18/07/2015 for correct use below. Such integer definition was probably in Swich.dat (commented below) but was missing here.
c----------------------------------------------------------------------

       SW2 = 1     !four stream
       SW4 = 0     !output of atmospheric profiles after importing from MISTRA
       SW7 = 1     !photolysis rates by band calculation
       SW12= 0     !output of fluxes
       SW13= 0     !output of check_four


      IF (SW7.eq.1) THEN
         DO L = 1,NW
            NWS(L) = NBAN(L)
         ENDDO
      ELSE
         DO L = 1,NW
            NWS(L) = L
         ENDDO
      ENDIF

      WRITE(*,'(8I4)')NWS

!      CALL READ_DATA(TEMP,  PRESS,   RELO3,
!     $               ITYPE, ALBEDO,  U0,    SCALEO3,
!     $               SW4)
      CALL READ_DATA(TEMP,  PRESS,   RELO3,
     $                      ALBEDO,  U0,    SCALEO3,
     $               SW4)

      CALL COLUMN(
     $            PRESS,  TEMP,     RELO3,     U0,
     $          SCALEO3,
     $               V2,   V2S,       DV2,   DV2S,
     $               V3,   V3S,       DV3,   DV3S,
     $             DENS)


      IF (SW7 .eq. 1) THEN
         WRITE(*,*)'band model, wavelength range: 178 - 690 nm'
      ELSE
         WRITE(*,'(A16,2F8.2,A5)')'wavelength range',
     $        (WAVE(NWS(L))*1.E7,l=1,NW,NW-1),' [nm]'
      ENDIF

      CALL CROSS_ATM(
     $                         V2S,      TEMP,
     $         CST_O3,
     $         CST_O2,
     $         QYNO3n,        QYO1D)

      CALL OPTIC_DATA(
     $         NWS,
     $         DV2,      DV3,       V2S,      U0,
     $         CST_O2,   CST_O3,    ALBEDO,
     $         TAUS_CLR, TAUA_CLR,  TAUS_AER, TAUA_AER,
     $         TAUS_CLD, TAUA_CLD,
     $         PRAY,     PCLD,      PAER,
     $         FLX,      ALB)
c end of initialization

      IF (SW2 .eq. 1) THEN
         CALL FOUR_INTF(
     $        TAUS_CLR, TAUA_CLR, TAUS_AER,  TAUA_AER,
     $        TAUS_CLD, TAUA_CLD,
     $        PRAY,     PAER,     PCLD,
     $        ALB,      FLX,       U0,
     $        FACT,     SW13)

         IF (SW12 .eq. 1) THEN
            CALL FLUX_OUT('f4st',FACT,FLX, NWS)
         ENDIF

      ENDIF


      IF (SW7 .eq. 1) THEN
        WRITE(*,*)'band'

      IF (NW.ne.7) THEN
         WRITE(*,*)'*************','ERROR IN band.f','*************'
         RETURN
      ENDIF

      CALL PHOTO_CAL(
     $           V2S,       V3S,
     $           U0,        TEMP,     RELO3,
     $           FACT,      PRESS,
     $           RJO3P,     RJO1D,    RJNO2,    RJHNO3,
     $           RJCOH2,    RJCHOH,   RJN2O5,   RJHNO4,
     $           RJNO2O,    RJNOO2,   RJH2O2,   RJCH3OOH,
     $           RJO2,      RJCFC11,  RJCFC12,  RJN2O,
     $           RJClONO2,  RJBrNO3,  RJCl2O2,  RJHOCl,
     $           RJBrCl_noT,RJClNO2,  RJBrNO2,  RJBr2,
     $           RJIO,      RJINO3,   RJCH3I,   RJI2,
     $           RJICl,     RJIBr,    RJC3H7I,  RJCH2ClI,
     $           RJCH2I2,   RJINO2,   RJBrO_noT,RJOClO_noT,
     $           RJCl2_noT, RJHOI_jen91,RJHOBr, RJHONO,
     $           RJNO2m,    RJNO3n,   RJdumm23, RJdumm24,
     $           RJdumm25,  RJdumm26,
     $           H_O2,      H_O3,     H_O4,   QYNO3n)

      ENDIF

c jrates for MISTRA
      do k=MAXLAY,MAXLAY-n_m+1,-1
         photol_j( 1,MAXLAY-k+1)=RJNO2(k)
         photol_j( 2,MAXLAY-k+1)=RJNOO2(k)
         photol_j( 3,MAXLAY-k+1)=RJO1D(k)
         photol_j( 4,MAXLAY-k+1)=RJHONO(k)
         photol_j( 5,MAXLAY-k+1)=RJHNO3(k)
         photol_j( 6,MAXLAY-k+1)=RJH2O2(k)
         photol_j( 7,MAXLAY-k+1)=RJHNO4(k)*2./3. ! 1.channel
         photol_j( 8,MAXLAY-k+1)=RJCHOH(k)
         photol_j( 9,MAXLAY-k+1)=RJCOH2(k)
         photol_j(10,MAXLAY-k+1)=RJNO2O(k)
         photol_j(11,MAXLAY-k+1)=RJHNO4(k)*1./3. ! 2.channel
         photol_j(12,MAXLAY-k+1)=RJN2O5(k)
         photol_j(13,MAXLAY-k+1)=RJHOCl(k)
         photol_j(14,MAXLAY-k+1)=RJClONO2(k)
         photol_j(15,MAXLAY-k+1)=RJBrNO3(k)
         photol_j(16,MAXLAY-k+1)=RJCl2O2(k)
         photol_j(17,MAXLAY-k+1)=RJCH3OOH(k)
         photol_j(18,MAXLAY-k+1)=RJClNO2(k)
         photol_j(19,MAXLAY-k+1)=RJCl2_noT(k)
         photol_j(20,MAXLAY-k+1)=RJHOBr(k)
         photol_j(21,MAXLAY-k+1)=RJBrNO2(k)
         photol_j(22,MAXLAY-k+1)=RJBr2(k)
         photol_j(23,MAXLAY-k+1)=RJBrCl_noT(k)
         photol_j(24,MAXLAY-k+1)=RJBrO_noT(k)
         photol_j(25,MAXLAY-k+1)=RJIO(k)
         photol_j(26,MAXLAY-k+1)=RJHOI_jen91(k)
         photol_j(27,MAXLAY-k+1)=RJI2(k)
         photol_j(28,MAXLAY-k+1)=RJICl(k)
         photol_j(29,MAXLAY-k+1)=RJIBr(k)
         photol_j(30,MAXLAY-k+1)=RJINO3(k)
         photol_j(31,MAXLAY-k+1)=RJCH3I(k)
         photol_j(32,MAXLAY-k+1)=RJC3H7I(k)
         photol_j(33,MAXLAY-k+1)=RJCH2ClI(k)
         photol_j(34,MAXLAY-k+1)=RJCH2I2(k)
         photol_j(35,MAXLAY-k+1)=RJOClO_noT(k)
         photol_j(36,MAXLAY-k+1)=9.*photol_j(16,MAXLAY-k+1) ! I2O2
         photol_j(37,MAXLAY-k+1)=RJINO2(k)
         photol_j(38,MAXLAY-k+1)=RJNO2m(k)
         photol_j(39,MAXLAY-k+1)=RJNO3n(k)
         photol_j(40,MAXLAY-k+1)=photol_j(35,MAXLAY-k+1) !OIO
         photol_j(41,MAXLAY-k+1)=RJdumm24(k)
         photol_j(42,MAXLAY-k+1)=RJdumm25(k)
         photol_j(43,MAXLAY-k+1)=RJdumm26(k)
         photol_j(44,MAXLAY-k+1)=1./17.*photol_j(34,MAXLAY-k+1) !CH2BrI
c         photol_j(45,MAXLAY-k+1)= ??? !CHBr2I
         photol_j(46,MAXLAY-k+1)=photol_j(31,MAXLAY-k+1) !C2H5I
         photol_j(47,MAXLAY-k+1)=RJO3P(k)
      enddo

      END Subroutine photol

*******************************************************************


*******************************************************************
!     SUBROUTINE READ_DATA(TEMP,    PRESS,  RELO3,
!    $                     ITYPE,   ALBEDO,  U0,      SCALEO3,
!    $                     SW4)
      SUBROUTINE READ_DATA(TEMP,    PRESS,  RELO3,
     $                              ALBEDO,  U0,      SCALEO3,
     $                     SW4)

! Description:
!   changed to get model init from MISTRA

c   this subroutine links MISTRA with Jochen Landgraf's code
c   it transfers the profiles from MISTRA to the relevant quantities
c   in Jochen's code. Lower case variables are from MISTRA and are passed in common
c   blocks. Upper case variables are from Jochen's code and are passed explicitly.


! 13/08/2016 jjb
!     changes done in this SR:
!        - removal of 2 unused arguments
!        - USE global_params for vertical grid parameters
!        - removal of 2 unused parameters
!        - a few unused variables commented after checking consistensistency all over this file
!        - update of call after removing arguments for atm_out
!        - little cleaning
!        - add of a switch (SW4) to decide whether output file should be written or not
!        - implicit none
!
! 21/10/2016 jjb
!     several bug fix related to array max size.
!        - variables in cb02 (t_m, p_m, ... )
!        - qmo3, z_mi
!        - also used thk_m instead of recalculating layer thicknesses (see below)
!
! 23/10/2016 jjb
!      removed Z and DZ from argument list, removed DZ declaration (based on thk_m : but unused)
!      removed FRAC, which was initialised to 0, but used nowhere in jrate
!
! 13/11/2017 jjb
!      removed ntypa and ntypd from /cb02/, defined but unused


      USE config, ONLY :
     &     scaleo3_m

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay,
     &     nrlev

      IMPLICIT NONE

      INTEGER MAXLAY,MAXWAV,SW4
      PARAMETER(MAXLAY=nrlay, MAXWAV=176)

C     OUTPUTS

      DOUBLE PRECISION
     $     TEMP(0:MAXLAY),   !temperature at model levels [K]
     $     PRESS(0:MAXLAY),  !pressure at model levels [hPa]
     $     RELO3(0:MAXLAY),  !O3 volume mixing ratio [1]
     $     ALBEDO(MAXWAV),   !shortwave albedo [1]
     $     U0,               !cosine of solar zenith angle
     $     SCALEO3,          !total vert. ozone column for scaling
     $     Z(MAXLAY)         !hight of levels [m]                      ! jjb for output purpose only
!     INTEGER
!    $     ITYPE(MAXLAY)     !aerosol type

c commom blocks from MISTRA:
      common /cb02/ t_m(nrlev),p_m(nrlev),rho_m(nrlev),xm1_m(nrlev),ts_m
      double precision t_m,p_m,rho_m,xm1_m,ts_m

      common /cb16/ u0_m,albedo_m(6),thk_m(nrlay)
      double precision u0_m, albedo_m, thk_m

      common /ozon/ qmo3_m(nrlev)
      double precision qmo3_m

      common /height/ z_mi(nrlev)
      double precision z_mi

! indexes for do loops
      INTEGER K,L

C----------------------------------------------------------------------

c      SCALEO3 = 300.
      SCALEO3 = scaleo3_m
c read in from istart, used ONLY here (i.e. NOT for calculations with
c PIFM for heating rates etc.)

c----------------------------------------------------------------------
c     SOLAR ZENITH ANGLE
c----------------------------------------------------------------------

      U0 = u0_m
      write (*,*) 'u0(band)= ',U0

c-----------------------------------------------------------------
c     Fix the model atmosphere here
c-----------------------------------------------------------------


C----------------------------------------------------------------------
C change the grid refer to altitudes
c----------------------------------------------------------------------

      DO K=1,MAXLAY
         Z(K)      = z_mi(K+1)
         TEMP(K)   = t_m(K+1)
         RELO3(K)  = qmo3_m(K+1)
         PRESS(K)  = p_m(K+1)/100. !Pa --> hPa
      ENDDO

c     PRESSURE, TEMPERATURE AND MIXING RATIOS AT VIRTUAL LEVEL 0
         PRESS(0)  = 0.37 * PRESS(1)  !0.37 = 1/e
c        0.63 = 1 -1/e , linear extrapolation in pressure
         TEMP(0)   = (TEMP(2)-TEMP(1))/(PRESS(2)-PRESS(1)) *
     $                 (-0.63) * PRESS(1) + TEMP(1)
         RELO3(0)  = (RELO3(2)-RELO3(1))/(PRESS(2)-PRESS(1)) *
     $                 (-0.63) * PRESS(1) + RELO3(1)


C----------------------------------------------------------------------
C     CLOUD AND AEROSOL DATA, SURFACE ALBEDO
C----------------------------------------------------------------------

C     SURFACE ALBEDO
      DO L=1,MAXWAV
         ALBEDO(L)=albedo_m(1)
      ENDDO

C     LIQUID WATER CONTENT FOR CLOUDS

C     CLOUD LIQUID WATER PATH

c     AEROSOL TYPES
c all aerosol data needed are tau, g see SR optic_data


C----------------------------------------------------------------------
C     CHECK DATA
C----------------------------------------------------------------------

!     IF(SW4.EQ.1) CALL ATM_OUT(ALBEDO,PRESS,TEMP,Z,RELO3,ITYPE)
      IF(SW4.EQ.1) CALL ATM_OUT(ALBEDO,PRESS,TEMP,Z,RELO3)

      END SUBROUTINE READ_DATA
*******************************************************************


********************************************************************
!     subroutine atm_out(ALBEDO,PRESS,TEMP,Z,RELO3,ITYPE)
      subroutine atm_out(ALBEDO,PRESS,TEMP,Z,RELO3)

! This SR writes atmospheric profiles after importing them from Mistra in the previous SR (read_data)
! In SR photol, this can be switched on/off using SW4

! 13/08/2016 jjb
!     changes done in this SR:
!        - removal of 4 unused arguments
!        - USE global_params for vertical grid parameters
!        - removal of 2 unused parameters
!        - a few unused variables commented, and write commented
!        - file unit defined as a parameter for easier change
!        - implicit none
!
! 23/10/2016 jjb
!     removed FRAC, which was initialised to 0 in READ_DATA, and used nowhere
!
! 05/06/2021 jjb
!     commented out ITYPE output, this variable is not linked to actual mistra aerosol type in jrate


      USE file_unit, ONLY :
! Imported Parameters:
     &     fun=>jpfunprofj   ! file unit number

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      IMPLICIT NONE

      INTEGER MAXLAY,MAXWAV
      PARAMETER(MAXLAY=nrlay, MAXWAV=176)

      DOUBLE PRECISION
     $     ALBEDO(MAXWAV),   !shortwave albedo [1]
     $     PRESS(0:MAXLAY),  !pressure at model levels [hPa]
     $     TEMP(0:MAXLAY),   !temperature at model levels [K]
     $     Z(MAXLAY),        !hight of levels [m]
     $     RELO3(0:MAXLAY)   !O3 volumn mixing ratio [1]

!      INTEGER
!     $     ITYPE(MAXLAY)     !aerosol type

      INTEGER K


 11   FORMAT(5(1P,E9.2))
 12   FORMAT(1P,6E11.4)
 13   FORMAT(1P,6(E9.2,2X))

      OPEN(FUN,FILE='prof.out',status='unknown')

      WRITE(FUN,*)'SURFACE ALBEDO [1]'
      WRITE(FUN,11)ALBEDO(1)

      WRITE(FUN,*)'SURFACE PRESSURE [mbar]'
      WRITE(FUN,12)PRESS(MAXLAY)

      WRITE(FUN,*)'TEMPERATURE PROFILE [K]'
      WRITE(FUN,13)(TEMP(K),K=1,MAXLAY)

      WRITE(FUN,*)'ALTIUDE GRID [m]'
      WRITE(FUN,13) (Z(K),K=1,MAXLAY)

      WRITE(FUN,*)'MIXING RATIO OF OZON [1]'
      WRITE(FUN,13)(RELO3(K), K=1,MAXLAY)

!      WRITE(FUN,*)'AEROSOL TYPE '//
!     $          '(1=rural,2=maritime,3=urban,4=free troposphere)'
!      WRITE(FUN,'(10I2)')(ITYPE(K),K=1,MAXLAY)

      CLOSE(FUN)

      END subroutine atm_out
**************************************************************************


********************************************************************
      SUBROUTINE COLUMN(
     $               PRESS,    TEMP,   RELO3,     U0,
     $             SCALEO3,
     $                  V2,     V2S,     DV2,   DV2S,
     $                  V3,     V3S,     DV3,   DV3S,
     $                DENS)

********************************************************************
*  calculate column densities of O2 and O3                         *
********************************************************************


! 13/08/2016 jjb
!     changes done in this SR:
!        - USE global_params for vertical grid parameters
!        - removal of 2 unused parameters
!        - little cleaning
!        - implicit none
!
! 13/05/2021 jjb
!        - USE constants


      USE constants, ONLY :
! Imported Parameters:
     &     AIR_M => M_air,
     &     AVOGA => Avogadro,
     &     GRAV  => g

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      IMPLICIT NONE

      INTEGER MAXLAY
      PARAMETER(MAXLAY=nrlay)

C     input:

      DOUBLE PRECISION
     $     TEMP(0:MAXLAY),  !temperature at model levels [K]
     $     RELO3(0:MAXLAY), !ozone mixing ratio at model levels [1]
     $     PRESS(0:MAXLAY), !pressure at model levels [hPa]
     $     U0,             !cosine of SZA
     $     SCALEO3             !total vert. ozone column for scaling

c      output:

       DOUBLE PRECISION
     $     V2(0:MAXLAY),    !O2 column density [part./cm^2]
     $     DV2(MAXLAY),     !diff. O2 column density [part./cm^2]
     $     V2S(0:MAXLAY),   !slant O2 column density [part./cm^2]
     $     DV2S(MAXLAY),    !diff. slant O2 column density [part./cm^2]
     $     V3(0:MAXLAY),    !O3 column density [part./cm^2]
     $     DV3(MAXLAY),     !diff. O3 column density [part./cm^2]
     $     V3S(0:MAXLAY),   !slant O3 column density [part./cm^2]
     $     DV3S(MAXLAY),    !diff. slant O2 column density [part./cm^2]
     $     DENS(0:MAXLAY)   !density of air molecules [part./cm^3]

c      internal variables

       DOUBLE PRECISION
     $     SECA            !secants od SZA
C--------------------------------------------------------------------------
      INTEGER K

      !DOUBLE PRECISION GRAV,BOLTZ,RELO2,AVOGA,AIR_M,SP,CONST,CONSTANT
      DOUBLE PRECISION BOLTZ,RELO2,SP,CONST,CONSTANT


      !GRAV=9.81
      BOLTZ=1.381D-23
      RELO2 = 0.2095D0
      !AVOGA = 6.022E23
      !AIR_M = 28.97D-3

      SP=AVOGA/(AIR_M*GRAV)*1.D-2               ![part./cm^2 * 1/hPa]
      CONST=SP * RELO2                          ![part./cm^2 * 1/hPa]

      CONSTANT=3.767D-20    !change units: part./cm^2 -> DU


C     O2 AND O3 COLUMNS AND SLANT COLUMNS

      IF (U0.gt.0.) THEN
         SECA = 1.D0 / U0
      ELSE
         SECA = 0.D0
      ENDIF

      V2(0)=CONST*PRESS(0)
      V3(0)=0.7D0 * SP * PRESS(0) * RELO3(0)
      V2S(0)=SECA*V2(0)
      V3S(0)=SECA*V3(0)
      DENS(0)=100.D0 * PRESS(0)/(BOLTZ*TEMP(0))*1.D-6

      DO K=1,MAXLAY
         V2(K)  = CONST*PRESS(K)
         V3(K)  = V3(K-1) + SP*(PRESS(K)-PRESS(K-1))*
     1            0.5D0 *( RELO3(K) + RELO3(K-1) )
         V2S(K) = SECA*V2(K)
         V3S(K) = SECA*V3(K)
         DENS(K)=100.D0 * PRESS(K)/(BOLTZ*TEMP(K))*1.D-6
      ENDDO

c     SCALING OF O3 COLUMN

      DO K=1,MAXLAY
         V3(K)  = V3(K) * SCALEO3/(V3(MAXLAY)*CONSTANT*1.D3)
         V3S(K) = SECA*V3(K)
      ENDDO

C     DIFFERENTIAL O2 AND O3 COLUMNS AND SLANT COLUMNS

      DV2(1)  = V2(1)
      DV3(1)  = V3(1)
      DV2S(1) = SECA*DV2(1)
      DV3S(1) = SECA*DV3(1)
      DO K=2,MAXLAY
         DV2(K)  = V2(K)-V2(K-1)
         DV3(K)  = V3(K)-V3(K-1)
         DV2S(K) = SECA*DV2(K)
         DV3S(K) = SECA*DV3(K)
      ENDDO

      WRITE(*,'(A10,1P,E12.4)')'O3 column:',V3(MAXLAY)*CONSTANT*1.D3
      WRITE(*,'(A10,1P,E12.4)')'O3 slantc:',V3S(MAXLAY)*CONSTANT*1.D3

      END SUBROUTINE COLUMN
**************************************************************************


**************************************************************************
      Subroutine CROSS_INIT
**********************************************************************
* INITIALIZATION OF ABSORPTION CROSS SECTIONS AND QUANTUM YIELDS     *
**********************************************************************


! 13/08/2016 jjb
!     changes done in this SR:
!        - USE global_params for vertical grid parameters and include data directory
!        - removal of 2 unused parameters
!        - unused COMMON commented, and corresponding read commented as well
!        - a few unused variables commented, and write commented
!        - COMMONS defined as data moved in a block data after this SR
!        - added mandatory repeat spec. for x descriptor (last two READ)
!
! 17/08/2016 jjb
!     further cleaning of unused, commented parts
!
! 05/11/2017 jjb: replaced directories by config

      USE config, ONLY : cinpdir_phot ! input directory for photolysis data files

      USE file_unit, ONLY :
! Imported Parameters:
     &     jpfunflx,
     &     jpfunsig,
     &     jpfuncheb

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      USE precision, ONLY :
! Imported Parameters:
     &     dp

      IMPLICIT NONE

      INTEGER MAXLAY,MAXWAV
      PARAMETER(MAXLAY=nrlay, MAXWAV=176)

      CHARACTER (LEN=14) :: CLFMT    ! format for read
      CHARACTER (LEN=130) :: CLPATH  ! path to input files
      INTEGER :: I, L                ! running indexes
      REAL (KIND=DP) :: sk, WL, X

      COMMON/WL/WAVE(MAXWAV),  !wavelength in the middle of the interval [cm]
     $          DWAVE(MAXWAV)  !width of the wavelength intervals [cm]
      COMMON/FL/FLUX(MAXWAV)   !extraterrestic flux per interval
c                               [photons/(cm^2 s)]
      REAL (KIND=DP) :: FLUX, WAVE, DWAVE

      COMMON/RAY_J/CS_RAY(MAXWAV)!Rayleigh scattering cross section [cm^2/part.]
      REAL (KIND=DP) :: CS_RAY

      COMMON/CROSS_SEC/        !cross sections [cm^2/part.]
     $     CS_H2O(MAXWAV),    CS_HNO3(MAXWAV),    CS_HNO4(MAXWAV),
     $     CS_SO2(MAXWAV),    CS_HCl(MAXWAV),     CS_HOCl(MAXWAV),
     $     CS_BrNO3(MAXWAV),  CS_CF3Cl(MAXWAV),   CS_CCl3F(MAXWAV),
     $     CS_CCl4(MAXWAV),   CS_CCl2O(MAXWAV),   CS_F115(MAXWAV),
     $     CS_F114(MAXWAV),   CS_F113(MAXWAV),    CS_CF2O(MAXWAV),
     $     CS_CClFO(MAXWAV),  CS_O2(MAXWAV),      CS_CH3OH(MAXWAV),
     $     CS_H2O2(MAXWAV),   CS_F22(MAXWAV),     CS_F13B1(MAXWAV),
     $     CS_F12B1(MAXWAV),  CS_CH3Br(MAXWAV),   CS_CCl2F2(MAXWAV),
     $     CS_CH3OOH(MAXWAV), CS_Cl2(MAXWAV),     CS_CHBr3(MAXWAV),
     $     CS_Cl2O2(MAXWAV),  CS_N2O5(MAXWAV),    CS_O4(MAXWAV),
     $     CS_NO3n(MAXWAV),   CS_O3H2O(MAXWAV),   CS_HOI_Jen91(MAXWAV),
     $     CS_HOCH2OOH(MAXWAV),CS_HOBr_JPL(MAXWAV),CS_HOBr(MAXWAV),
     $     CS_BrCl(MAXWAV),   CS_BrCl_noT(MAXWAV),   CS_ClNO2(MAXWAV),
     $     CS_BrNO2(MAXWAV),  CS_Br2(MAXWAV),     CS_IO(MAXWAV),
     $     CS_INO3(MAXWAV),   CS_CH3I(MAXWAV),    CS_I2(MAXWAV),
     $     CS_ICl(MAXWAV),    CS_IBr(MAXWAV),     CS_C3H7I(MAXWAV),
     $     CS_CH2ClI(MAXWAV), CS_CH2I2(MAXWAV),   CS_INO2(MAXWAV),
     $     CS_BrO_noT(MAXWAV), CS_OClO_noT(MAXWAV), CS_Cl2_noT(MAXWAV),
     $     CS_HONO(MAXWAV),
     $     CS_NO2m(MAXWAV),   CS_dumm23(MAXWAV),
     $     CS_dumm24(MAXWAV),  CS_dumm25(MAXWAV),   CS_dumm26(MAXWAV)
      REAL (KIND=DP) ::
     $     CS_H2O,    CS_HNO3,    CS_HNO4,
     $     CS_SO2,    CS_HCl,     CS_HOCl,
     $     CS_BrNO3,  CS_CF3Cl,   CS_CCl3F,
     $     CS_CCl4,   CS_CCl2O,   CS_F115,
     $     CS_F114,   CS_F113,    CS_CF2O,
     $     CS_CClFO,  CS_O2,      CS_CH3OH,
     $     CS_H2O2,   CS_F22,     CS_F13B1,
     $     CS_F12B1,  CS_CH3Br,   CS_CCl2F2,
     $     CS_CH3OOH, CS_Cl2,     CS_CHBr3,
     $     CS_Cl2O2,  CS_N2O5,    CS_O4,
     $     CS_NO3n,   CS_O3H2O,   CS_HOI_Jen91,
     $     CS_HOCH2OOH,CS_HOBr_JPL,CS_HOBr,
     $     CS_BrCl,   CS_BrCl_noT,   CS_ClNO2,
     $     CS_BrNO2,  CS_Br2,     CS_IO,
     $     CS_INO3,   CS_CH3I,    CS_I2,
     $     CS_ICl,    CS_IBr,     CS_C3H7I,
     $     CS_CH2ClI, CS_CH2I2,   CS_INO2,
     $     CS_BrO_noT, CS_OClO_noT, CS_Cl2_noT,
     $     CS_HONO,
     $     CS_NO2m,   CS_dumm23,
     $     CS_dumm24,  CS_dumm25,   CS_dumm26

c     cross sections CS_X at temperature T_X (X=specie)
      COMMON/CROSS_SEC_T/
     $     CS_O3(MAXWAV,3),   CS_NO3(MAXWAV,2),   CS_NO2(MAXWAV,2),
     $     CS_OCS(MAXWAV,2),  CS_ClONO2(MAXWAV,3),CS_CH3CCl3(MAXWAV,3),
     $     CS_CO2(MAXWAV,3),  COEFF_HNO3(MAXWAV), CS_HOI(MAXWAV,3),
     $     CS_CH2O(MAXWAV,2), CS_CH3Cl(MAXWAV,3),
     $     T_O3(3),           T_NO3(2),           T_NO2(2),
     $     T_OCS(2),          T_ClONO2(3),        T_CH3CCl3(3),
     $     T_CO2(3),          T_HOI(3),           T_CH2O(2),
     $     T_CH3Cl(3)
      REAL (KIND=DP) ::
     $     CS_O3, CS_NO3, CS_NO2, CS_OCS, CS_ClONO2, CS_CH3CCl3, CS_CO2,
     $     COEFF_HNO3, CS_HOI, CS_CH2O, CS_CH3Cl,
     $     T_O3, T_NO3, T_NO2, T_OCS, T_ClONO2, T_CH3CCl3, T_CO2,
     $     T_HOI, T_CH2O, T_CH3Cl

C     coefficients for Koppers and Murtagh  parameterization of CS_O2 in
C     Schuhmann-Runge band

      COMMON/CHEB_COEFF/CHEB_COEFF_A(20,13),
     $                  CHEB_COEFF_B(20,13)
      REAL (KIND=DP) :: CHEB_COEFF_A, CHEB_COEFF_B



C--------------------------------------------------------------------
C     wavelengths in the middle of the intervals
C--------------------------------------------------------------------
C     see: C. Bruehl, P.J. Crutzen, Scenarios of possible changes in ....
C          Climate Dynamics (1988)2: 173-203

      DO 1 L=1,13
        WAVE(L)=1.D0 /(56250.D0 -500.D0 *L)
    1 CONTINUE
       DO 2 L=14,45
       WAVE(L)=1.D0 /(49750.D0 -(L-13)*500.D0 )
    2 CONTINUE
      DO 3 L=46,68
        WAVE(L)=(266.D0 +(L-13))*1.D-7
    3 CONTINUE
      DO 4 L=69,71
        WAVE(L)=(320.5D0 +2.D0 *(L-68))*1.D-7
    4 CONTINUE
      DO 5 L=72,176
        WAVE(L)=(325.D0 +5.D0 *(L-71))*1.D-7
    5 CONTINUE

      DO L=2,MAXWAV-1
         DWAVE(L)=0.5D0*(WAVE(L+1)-WAVE(L-1))
      ENDDO

      DWAVE(1)      = DWAVE(2)
      DWAVE(MAXLAY) = DWAVE(MAXLAY-1)

C---------------------------------------------------------------------
C     Rayleigh scattering cross sections
C     emp. formula of Nicolet   plan. space sci.32, 1467f (1984)
C---------------------------------------------------------------------

      DO  L=1,176
        WL=WAVE(L)*1.E4
        X=0.389D0 *WL+0.09426D0 /WL-0.3228D0
        CS_RAY(L)=4.02D-28/WL**(4.D0 +X)
      ENDDO

      CLPATH = trim(cinpdir_phot)//'flux.dat'
      OPEN(UNIT=jpfunflx, FILE=trim(CLPATH),STATUS='OLD')
      READ(jpfunflx,10) FLUX

      CLOSE(jpfunflx)

 10   FORMAT(1P,7E10.2)

      CLPATH = trim(cinpdir_phot)//'sig0900.dat'
      OPEN(UNIT=jpfunsig, FILE=trim(CLPATH), STATUS='OLD')

      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_H2O
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_HNO3
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_HNO4
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_SO2
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_HCl
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_HOCl
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_BrNO3
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CF3Cl
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CCl3F
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CCl4
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CCl2O
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_F115
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_F114
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_F113
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CF2O
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CClFO
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_O2
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CH3OH
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_H2O2
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_F22
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_F13B1
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_F12B1
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CH3Br
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CCl2F2
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CH3OOH
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_Cl2
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CHBr3
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_Cl2O2
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_N2O5
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_O4
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_NO3n
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_O3H2O
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_HOI_Jen91
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_HOCH2OOH
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_HOBr_JPL
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_HOBr
      READ(jpfunsig,*)
c      READ(jpfunsig,101)CS_BrCl
c      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_BrCl_noT
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_ClNO2
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_BrNO2
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_Br2
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_IO
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_INO3
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CH3I
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_I2
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_ICl
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_IBr
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_C3H7I
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CH2ClI
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_CH2I2
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_INO2
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_BrO_noT
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_OClO_noT
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_Cl2_noT
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_HONO
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_NO2m
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_dumm23
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_dumm24
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_dumm25
      READ(jpfunsig,*)
      READ(jpfunsig,101)CS_dumm26
      READ(jpfunsig,*)
      READ(jpfunsig,201)T_O3
      DO I=1,3
         READ(jpfunsig,101)(CS_O3(L,I),L=1,MAXWAV)
      ENDDO
      READ(jpfunsig,*)
      READ(jpfunsig,201)T_NO3
      DO I=1,2
         READ(jpfunsig,101)(CS_NO3(L,I),L=1,MAXWAV)
      ENDDO
      READ(jpfunsig,*)
      READ(jpfunsig,201)T_NO2
      DO I=1,2
         READ(jpfunsig,101)(CS_NO2(L,I),L=1,MAXWAV)
      ENDDO
      READ(jpfunsig,*)
      READ(jpfunsig,201)T_OCS
      DO I=1,2
         READ(jpfunsig,101)(CS_OCS(L,I),L=1,MAXWAV)
      ENDDO
      READ(jpfunsig,*)
      READ(jpfunsig,201)T_ClONO2
      DO I=1,3
         READ(jpfunsig,101)(CS_ClONO2(L,I),L=1,MAXWAV)
      ENDDO
      READ(jpfunsig,*)
      READ(jpfunsig,201)T_CH3CCl3
      DO I=1,3
         READ(jpfunsig,101)(CS_CH3CCl3(L,I),L=1,MAXWAV)
      ENDDO
      READ(jpfunsig,*)
      READ(jpfunsig,201)T_CO2
      DO I =1,3
         READ(jpfunsig,101)(CS_CO2(L,I),L=1,MAXWAV)
      ENDDO
      READ(jpfunsig,*)
      READ(jpfunsig,101)COEFF_HNO3
      READ(jpfunsig,*)
      READ(jpfunsig,201)T_HOI
      DO I =1,3
         READ(jpfunsig,101)(CS_HOI(L,I),L=1,MAXWAV)
      ENDDO
      READ(jpfunsig,*)
      READ(jpfunsig,201)T_CH2O
      DO I =1,2
         READ(jpfunsig,101)(CS_CH2O(L,I),L=1,MAXWAV)
      ENDDO
      READ(jpfunsig,*)
      READ(jpfunsig,201)T_CH3Cl
      DO I=1,3
         READ(jpfunsig,101)(CS_CH3Cl(L,I),L=1,MAXWAV)
      ENDDO

 101  FORMAT(7E10.2)
 201  FORMAT(4F5.0)

      CLOSE(jpfunsig)


C     CEBESHEV COEFFICIENT A AND B FOR KOPPER MURTG. PARAMETERIZATION

      CLPATH = trim(cinpdir_phot)//'cheb_coeff.dat'
      CLFMT = '(17(E23.6,1X))'
      OPEN(jpfuncheb,file=trim(CLPATH),status='old')

      READ(jpfuncheb,*)
      READ(jpfuncheb,*)
      DO I=1,20
         READ(jpfuncheb,CLFMT)sk,sk,(CHEB_COEFF_A(I,L),L=1,13),sk,sk ! sk = skip
      ENDDO

      READ(jpfuncheb,*)
      READ(jpfuncheb,*)
      DO I=1,20
         READ(jpfuncheb,CLFMT)sk,sk,(CHEB_COEFF_B(I,L),L=1,13),sk,sk ! sk = skip
      ENDDO
      CLOSE(jpfuncheb)

      END SUBROUTINE CROSS_INIT
*****************************************************************


*****************************************************************
      block data cross_init_data

! 13/08/2016 jjb
!     the data here come from the abose SR "cross_init"
!     Definition of COMMON variable as data in a SR in incorrect
!     The correct way to proceed is to define data in block data

      implicit none

      integer i, l

C     coefficients for H2O2 cross section (JPL 94)
      COMMON/C_H2O2/ A_H2O2(0:7),B_H2O2(0:4)
      DOUBLE PRECISION A_H2O2, B_H2O2

      DATA A_H2O2 /
     $      6.4761000D+04, -9.2170972D+02,    4.5356490D+00,
     $     -4.4589016D-03, -4.0351010D-05,    1.6878206D-07,
     $     -2.6520140D-10,  1.5534675D-13/
      DATA B_H2O2 /
     $      6.8123000D+03, -5.1351000D+01,    1.1522000D-01,
     $     -3.0493000D-05, -1.0924000D-07/


C     coefficients for N2O cross section (JPL 94)
      COMMON/C_N2O/  A_N2O(0:4),  B_N2O(0:3)
      DOUBLE PRECISION A_N2O, B_N2O

      DATA A_N2O /
     $      6.821023D+01,  -4.071805D+00,      4.301146D-02,
     $     -1.777846D-04,   2.520672D-07/
      DATA B_N2O /
     $      1.234014D+02,  -2.116255D+00,      1.111572D-02,
     $     -1.881058D-05/

C     coefficients for NO cross section (Allen, Frederic)
c$$$      COMMON/C_NO/ B_NO(2,5),  A_NO(2,9)
c$$$      DOUBLE PRECISION B_NO, A_NO
c$$$
c$$$      DATA A_NO/
c$$$     $     -1.790868D+1, -1.654245D+1, -1.924701D-1, +5.836899D-1,
c$$$     $     -7.217717D-2, +3.449436D-1,  5.648282D-2,  1.700653D-1,
c$$$     $      4.569175D-2, -3.324717D-2,  8.353572D-3, -4.952424D-2,
c$$$     $      0.0D0,        1.579306D-2,  0.0D0,        1.835462D-2,
c$$$     $      0.0D0,        3.368125D-3/
c$$$
c$$$      DATA B_NO/
c$$$     $      7.836832D+3,  1.297581D+4, -1.549880D+3, -2.582981D+3,
c$$$     $      1.148342D+2,  1.927709D+2, -3.777754,    -6.393008,
c$$$     $      4.655696D-2,  7.949835D-2/

C     coefficients for O(1D) QUANTUM YIELD (MICHELSON)
      COMMON/C_O1D/ A_O1D(19),  B_O1D(19)
      DOUBLE PRECISION A_O1D, B_O1D

      DATA A_O1D/1.01D0, 1.01D0, 1.05D0, 1.15D0, 1.39D0, 1.90D0, 2.93D0,
     1           4.87D0, 8.21D0, 13.3D0, 17.6D0, 20.4D0, 18.0D0, 21.8D0,
     2           18.1D0, 17.2D0, 7.99D0, 12.9D0, 11.25D0 /
      DATA B_O1D/3.933D0, 11.51D0, 33.09D0, 79.39D0, 159.9D0, 272.5D0,
     1           407.9D0, 551.4D0, 682.3D0, 791.6D0, 851.3D0, 903.8D0,
     2           900.3D0, 948.4D0, 891.1D0, 1066.D0, 969.4D0, 1191.5D0,
     3           1293.5D0 /

c     coefficients for optical depth of Schumann-Runge band above TOA
      COMMON/C_O2_TOP/ CT_TOP(4,13)
      DOUBLE PRECISION CT_TOP

      DATA ((CT_TOP(I,L),I=1,4),L=1,13) /
     $     -2.5488D+02,  1.5900D+01, -3.4078D-01,  2.5083D-03,
     $     -5.8222D+02,  3.5825D+01, -7.4328D-01,  5.2068D-03,
     $     -5.8239D+02,  3.5637D+01, -7.3537D-01,  5.1210D-03,
     $     -5.6359D+02,  3.4235D+01, -7.0220D-01,  4.8652D-03,
     $     -5.5623D+02,  3.3538D+01, -6.8358D-01,  4.7115D-03,
     $     -6.4776D+02,  3.8519D+01, -7.7292D-01,  5.2339D-03,
     $     -5.7035D+02,  3.3504D+01, -6.6617D-01,  4.4825D-03,
     $     -5.7514D+02,  3.3451D+01, -6.5964D-01,  4.4075D-03,
     $     -9.3045D+02,  5.3921D+01, -1.0505D+00,  6.8803D-03,
     $     -8.9272D+02,  5.1460D+01, -1.0005D+00,  6.5579D-03,
     $     -7.1078D+02,  4.0599D+01, -7.8842D-01,  5.1978D-03,
     $     -1.4366D+02,  6.1527D+00, -9.5919D-02,  5.8395D-04,
     $     -1.1535D+02,  4.5631D+00, -6.6966D-02,  4.1305D-04/

      end block data cross_init_data
*****************************************************************


*****************************************************************
      SUBROUTINE CROSS_ATM(
     $                         V2S,      TEMP,
     $         CST_O3,
     $         CST_O2,
     $         QYNO3n,        QYO1D)

******************************************************************
* CALCULATION OF CROSS SECTION DEPENDING ON ATMOSPHERICAL        *
* CONDITIONS, LIKE TEMERATURE AND O2-COLUMN                      *
******************************************************************


! 13/08/2016 jjb
!     changes done in this SR:
!        - USE global_params for vertical grid parameters
!        - removal of 2 unused parameters
! 16/08/2016 jjb
!     further changes in this SR:
!        - most CST_* deleted, as unused, except CST_O3 and CST_O2
!        - removal of all relevant parts (common blocks, internal functions, calculations)
!        - re-indexing of CST_O3 for CPU efficiency
! 25/06/2021 jjb: all missing declarations and implicit none


      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      USE precision, ONLY :
! Imported Parameters:
     &     dp

      IMPLICIT NONE

      INTEGER MAXLAY,MAXWAV
      PARAMETER(MAXLAY=nrlay, MAXWAV=176)

C     inputs

      DOUBLE PRECISION
     $     V2S(0:MAXLAY),       !slant O2 column density [part./cm^2]
     $     TEMP(0:MAXLAY)       !temperature at model levels [K]

C     outputs

      DOUBLE PRECISION  !temp. dependent cross section for profile TEMP
     $     CST_O3(0:MAXLAY,MAXWAV),
     $     CST_O2(MAXWAV,0:MAXLAY)

      DOUBLE PRECISION !temp. dependent quantum yield for profile TEMP
     $     QYNO3n(0:MAXLAY),            QYO1D(MAXWAV,0:MAXLAY)

c     cross section initzialization

      COMMON/CROSS_SEC/        !cross sections [cm^2/part.]
     $     CS_H2O(MAXWAV),    CS_HNO3(MAXWAV),    CS_HNO4(MAXWAV),
     $     CS_SO2(MAXWAV),    CS_HCl(MAXWAV),     CS_HOCl(MAXWAV),
     $     CS_BrNO3(MAXWAV),  CS_CF3Cl(MAXWAV),   CS_CCl3F(MAXWAV),
     $     CS_CCl4(MAXWAV),   CS_CCl2O(MAXWAV),   CS_F115(MAXWAV),
     $     CS_F114(MAXWAV),   CS_F113(MAXWAV),    CS_CF2O(MAXWAV),
     $     CS_CClFO(MAXWAV),  CS_O2(MAXWAV),      CS_CH3OH(MAXWAV),
     $     CS_H2O2(MAXWAV),   CS_F22(MAXWAV),     CS_F13B1(MAXWAV),
     $     CS_F12B1(MAXWAV),  CS_CH3Br(MAXWAV),   CS_CCl2F2(MAXWAV),
     $     CS_CH3OOH(MAXWAV), CS_Cl2(MAXWAV),     CS_CHBr3(MAXWAV),
     $     CS_Cl2O2(MAXWAV),  CS_N2O5(MAXWAV),    CS_O4(MAXWAV),
     $     CS_NO3n(MAXWAV),   CS_O3H2O(MAXWAV),   CS_HOI_Jen91(MAXWAV),
     $     CS_HOCH2OOH(MAXWAV),CS_HOBr_JPL(MAXWAV),CS_HOBr(MAXWAV),
     $     CS_BrCl(MAXWAV),   CS_BrCl_noT(MAXWAV),   CS_ClNO2(MAXWAV),
     $     CS_BrNO2(MAXWAV),  CS_Br2(MAXWAV),     CS_IO(MAXWAV),
     $     CS_INO3(MAXWAV),   CS_CH3I(MAXWAV),    CS_I2(MAXWAV),
     $     CS_ICl(MAXWAV),    CS_IBr(MAXWAV),     CS_C3H7I(MAXWAV),
     $     CS_CH2ClI(MAXWAV), CS_CH2I2(MAXWAV),   CS_INO2(MAXWAV),
     $     CS_BrO_noT(MAXWAV), CS_OClO_noT(MAXWAV), CS_Cl2_noT(MAXWAV),
     $     CS_HONO(MAXWAV),
     $     CS_NO2m(MAXWAV),   CS_dumm23(MAXWAV),
     $     CS_dumm24(MAXWAV),  CS_dumm25(MAXWAV),   CS_dumm26(MAXWAV)
      REAL (KIND=DP) ::
     $     CS_H2O,    CS_HNO3,    CS_HNO4,
     $     CS_SO2,    CS_HCl,     CS_HOCl,
     $     CS_BrNO3,  CS_CF3Cl,   CS_CCl3F,
     $     CS_CCl4,   CS_CCl2O,   CS_F115,
     $     CS_F114,   CS_F113,    CS_CF2O,
     $     CS_CClFO,  CS_O2,      CS_CH3OH,
     $     CS_H2O2,   CS_F22,     CS_F13B1,
     $     CS_F12B1,  CS_CH3Br,   CS_CCl2F2,
     $     CS_CH3OOH, CS_Cl2,     CS_CHBr3,
     $     CS_Cl2O2,  CS_N2O5,    CS_O4,
     $     CS_NO3n,   CS_O3H2O,   CS_HOI_Jen91,
     $     CS_HOCH2OOH,CS_HOBr_JPL,CS_HOBr,
     $     CS_BrCl,   CS_BrCl_noT,   CS_ClNO2,
     $     CS_BrNO2,  CS_Br2,     CS_IO,
     $     CS_INO3,   CS_CH3I,    CS_I2,
     $     CS_ICl,    CS_IBr,     CS_C3H7I,
     $     CS_CH2ClI, CS_CH2I2,   CS_INO2,
     $     CS_BrO_noT, CS_OClO_noT, CS_Cl2_noT,
     $     CS_HONO,
     $     CS_NO2m,   CS_dumm23,
     $     CS_dumm24,  CS_dumm25,   CS_dumm26

c     cross sections CS_X at temperature T_X (X=specie)
      COMMON/CROSS_SEC_T/
     $     CS_O3(MAXWAV,3),   CS_NO3(MAXWAV,2),   CS_NO2(MAXWAV,2),
     $     CS_OCS(MAXWAV,2),  CS_ClONO2(MAXWAV,3),CS_CH3CCl3(MAXWAV,3),
     $     CS_CO2(MAXWAV,3),  COEFF_HNO3(MAXWAV), CS_HOI(MAXWAV,3),
     $     CS_CH2O(MAXWAV,2), CS_CH3Cl(MAXWAV,3),
     $     T_O3(3),           T_NO3(2),           T_NO2(2),
     $     T_OCS(2),          T_ClONO2(3),        T_CH3CCl3(3),
     $     T_CO2(3),          T_HOI(3),           T_CH2O(2),
     $     T_CH3Cl(3)
      REAL (KIND=DP) ::
     $     CS_O3, CS_NO3, CS_NO2, CS_OCS, CS_ClONO2, CS_CH3CCl3, CS_CO2,
     $     COEFF_HNO3, CS_HOI, CS_CH2O, CS_CH3Cl,
     $     T_O3, T_NO3, T_NO2, T_OCS, T_ClONO2, T_CH3CCl3, T_CO2,
     $     T_HOI, T_CH2O, T_CH3Cl

      COMMON/WL/WAVE(MAXWAV), !wavelength in the middle of the interval [cm]
     $          DWAVE(MAXWAV)  !width of the wavelength intervals [cm]

C     coefficients for O(1D) QUANTUM YIELD (MICHELSON)

      COMMON/C_O1D/ A_O1D(19),  B_O1D(19)

c     local arrays

      DOUBLE PRECISION
     $     SRO2(13,0:MAXLAY)

c     local scalars
      INTEGER K, L
      REAL (KIND=DP) :: A_O1D, B_O1D, C1_O3, C2_O3, C3_O3,
     $     DWAVE, TEMPER, WAVE

C--------------------------------------------------------------------
C temperature dependent cross section
C--------------------------------------------------------------------

      DO  L = 1,MAXWAV
         C1_O3    = CS_O3(L,1)
         C2_O3    = (CS_O3(L,2) - CS_O3(L,1))/(T_O3(2) - T_O3(1))
         C3_O3    = ((CS_O3(L,3) - CS_O3(L,2))/
     $               (T_O3(3) - T_O3(2)) - C2_O3 ) /
     $               (T_O3(3) - T_O3(1))

         DO K  = 0,MAXLAY
            CST_O3(K,L)= DMAX1(0.D0,
     $           ((TEMP(K)-T_O3(2))*C3_O3 +
     $           C2_O3)*(TEMP(K)-T_O3(1)) +
     $           C1_O3)
         ENDDO
      ENDDO


C-------------------------------------------------------------------
C O2 cross section in the Schumann-Runge bands
C-------------------------------------------------------------------

!      SEC = 1./U0

C     O2 cross section of Schumann-Runge-band

C     Allen Frederic Parameterization

c      CALL SR_O2_AF(MAXLAY,MJ,SEC,V2,PRESS,TEMP,SRO2)

C     Koopers Murtagh Parameterization

      CALL SR_O2_KM(MAXLAY,V2S,TEMP,SRO2)

      DO K  = 0,MAXLAY
c        Schumann-Runge band
         DO L=1,13
            CST_O2(L,K)=SRO2(L,K)
         ENDDO
c        Herzberg band
         DO L=14,MAXWAV
            CST_O2(L,K)=CS_O2(L)
         ENDDO
      ENDDO


C-----------------------------------------------------------------
C temperature dependent quantum yield of NO3- and O1D
C-----------------------------------------------------------------
      DO K  = 0,MAXLAY

C        NO3- QUANTUM YIELD IN DROPLETS (ZELLNER ET AL., JAC 10, 411, 1990)

         QYNO3n(K)=1.7D-2*DEXP(1800.D0 *(1.D0 / 298.D0 - 1./TEMP(K)))

c        QUANTUM YIELDS OF O3-> O1D AS GIVEN BY MICHELSEN ET AL.
C        GEOPHYS. RES. LETT. 21, 2227-2230, 1994

C        L=1 ,38  <->  179.4 - 268.5 nm
C        L=39,51  <->  272.1 - 304.0 nm
C        L=52,70  <->  305.0 - 324.5 nm
C        L=71,176 <->  326.5 - 850.0 nm

         TEMPER=TEMP(K)
         IF(TEMPER.LT.185) TEMPER=185
         IF(TEMPER.GT.320) TEMPER=320
         DO L=1,38
           QYO1D(L,K)=0.87D0
         ENDDO
         DO L=39,51
           QYO1D(L,K)=1.98D0 - 301.D0/(WAVE(L)*1.D7)
         ENDDO
         DO L=52,70
           QYO1D(L,K)=A_O1D(L-51)*EXP(-1.439D0 * B_O1D(L-51)/TEMPER)
         ENDDO
         DO L=71,176
           QYO1D(L,K)=0.
         ENDDO
      ENDDO

      END SUBROUTINE CROSS_ATM
*******************************************************************


*******************************************************************
c$$$  jjb 2015-11-21 : commented, as unused. Was formerly called in
c$$$    SR Cross_atm, now replaced by call SR_02_KM
c$$$
c$$$
c$$$      SUBROUTINE SR_O2_AF(MAXLAY,MJ,SEC,V2,PRESS,TEMP,SRO2)
c$$$
c$$$*******************************************************************
c$$$* SCHUMANN-RUNGE PARAMETERIZATION FOR EFFECTIVE O2 CROSS SECTION  *
c$$$* ALLEN M, FREDERICK JE, J. ATMOS. SCI. 39, 2066 FF               *
c$$$* SPECTRAL RANGE: 179.4 - 201.0 NM,  WMO(1985) intervals          *
c$$$*******************************************************************
c$$$
c$$$      IMPLICIT NONE
c$$$C     INPUTS
c$$$
c$$$      INTEGER MAXLAY,MJ
c$$$
c$$$      DOUBLE PRECISION
c$$$     $     SEC,             ! 1./COS(POLAR ANGLE)
c$$$     $     V2(0:MAXLAY),     ! O2 COLUMN
c$$$     $     PRESS(0:MAXLAY),  ! PRESSURE LEVELS
c$$$     $     TEMP(0:MAXLAY)    ! TEMPERATURE PROFILE
c$$$
c$$$C     OUTPUTS
c$$$
c$$$      DOUBLE PRECISION
c$$$     $     SRO2(13,0:MAXLAY) ! EFFECTIVE O2 CROSS SECTION AT LEVEL INTERFACE
c$$$
c$$$C     HERZBERG CONTINUUM BACKGROUND
c$$$
c$$$C     COEFFICIENTS AK(L,I), BK(L,I)
c$$$
c$$$      COMMON/C_O2/ AK(13,9),  BK(13,5)
c$$$      DOUBLE PRECISION AK, BK
c$$$
c$$$      INTEGER I, J, K, L
c$$$      DOUBLE PRECISION ALP, ALV2, CLN, SF, SLN, XL
c$$$
c$$$C----------------------------------------------------------------
c$$$
c$$$c      WRITE(*,*)' ALLEN FREDERIC'
c$$$
c$$$C     O2 ABSORPTION CROSS SECTIONS
c$$$
c$$$      DO J  = 1,MJ
c$$$      DO K=0,MAXLAY
c$$$
c$$$         ALP  = DMIN1(2.D0,DLOG10(PRESS(K)))
c$$$         ALV2 = DLOG10(V2(K))
c$$$
c$$$         DO L=1,13
c$$$
c$$$            IF (L.EQ.13) THEN
c$$$               XL=TEMP(K)
c$$$            ELSE
c$$$               XL=ALP
c$$$            ENDIF
c$$$
c$$$
c$$$            SLN=0.0
c$$$            CLN=0.0
c$$$
c$$$            DO I=1,9
c$$$               SLN=SLN + AK(L,I)*XL**(I-1)
c$$$            ENDDO
c$$$            DO I=1,5
c$$$               CLN=CLN + BK(L,I)*ALV2**(I-1)
c$$$            ENDDO
c$$$
c$$$            SF=SEC**(-10.**CLN)
c$$$            SRO2(L,K)=SF*10.**SLN
c$$$
c$$$C           CORRECTION FOR OVERESTIMATED BACKGROUND CROSS SECTION BY 40 %
c$$$
c$$$            IF (L .GE. 11) THEN
c$$$             SRO2(L,K)=(1.-0.4) * SRO2(L,K)
c$$$            ENDIF
c$$$         ENDDO
c$$$
c$$$      ENDDO
c$$$      ENDDO
c$$$
c$$$      END SUBROUTINE SR_O2_AF
***********************************************************************


***********************************************************************
      SUBROUTINE SR_O2_KM(MAXLAY,V2S,TEMP,SRO2)

****************************************************************
* SCHUMANN-RUNGE PARAMETERIZATION                              *
* Koppers GAA, Murtagh DP, Ann. Geophysicae 14, 68-79          *
****************************************************************

      IMPLICIT NONE

C     INPUTS

      INTEGER MAXLAY

      DOUBLE PRECISION
     $     V2S(0:MAXLAY),       !SLANT O2 COLUMN
     $     TEMP(0:MAXLAY)       !TEMPERATURE PROFILE

C     OUTPUTS

      DOUBLE PRECISION
     $     SRO2(13,0:MAXLAY)    !EFFECTIVE O2 CROSS SECTION


C     CEBESHEV COEFFICIENT A AND B

      COMMON/CHEB_COEFF/CHEB_COEFF_A(20,13),
     $                  CHEB_COEFF_B(20,13)
      DOUBLE PRECISION CHEB_COEFF_A, CHEB_COEFF_B

      DOUBLE PRECISION
     $                 COEFF_A(20),
     $                 COEFF_B(20)

C     EXTERNAL FUNCTIONS

      DOUBLE PRECISION, EXTERNAL :: CHEBEV

C     LOCAL SCALARS
      INTEGER I, K, L
      DOUBLE PRECISION A, B, DL_O2

c-----------------------------------------------------------------

C     CALCULATION OF COEFFICIENTS A AND B AND EFFECTIVE O2 CROSS SECTION

      DO L=1,13

         DO I=1,20
            COEFF_A(I) = CHEB_COEFF_A(I,L)
            COEFF_B(I) = CHEB_COEFF_B(I,L)
         ENDDO

         DO K=0,MAXLAY

            IF (V2S(K).ge.EXP(38.)) THEN
               DL_O2 =MIN(56.D0,DLOG(V2S(K)))

               A= CHEBEV(38.D0,56.D0,COEFF_A,20,DL_O2)
               B= CHEBEV(38.D0,56.D0,COEFF_B,20,DL_O2)

               SRO2(L,K)= EXP(A * (TEMP(K)-220.D0) + B)
            ELSE
               SRO2(L,K) = 0.D0
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE SR_O2_KM
*****************************************************************


*****************************************************************
      DOUBLE PRECISION FUNCTION CHEBEV(A,B,C,M,X)

*****************************************************************
*     CHEBESHEV POLYNOMIAL EVALUATION: All arguments are input. *
*     C(1:M) is an array of Chebyshev coefficients. The Cheb-   *
*     yshev ploynomial sum_{k=1}^{M} C_k*T_{k-1}(y)-c_1/2 is    *
*     evaluated at a point y=(x-(a+b)/2)/(a-b)/2.               *
*     see: Numerical recipes in fortran chapter 5.8 P.184 ff    *
*****************************************************************
      IMPLICIT NONE

C     INPUTS

      INTEGER M

      DOUBLE PRECISION  A,B,C(M),X

C     INTERNALS

      INTEGER J

      DOUBLE PRECISION D, DD, SV, Y, Y2
C----------------------------------------------------------------

      IF ( (X-A)*(X-B).GT.0.) print*,'X NOT IN RANGE IN CHEBEV !!!'

      D=0.D0
      DD=0.D0

      Y  = (2.D0*X-A-B)/(B-A)            !CHANGE OF VARIABLE TO RANGE [-1,1]
      Y2 = 2.D0*Y
      DO J = M,2,-1                   !CLENSHAW'S RECURRENCE
         SV = D
         D  = Y2*D-DD+C(J)
         DD = SV
      ENDDO

      CHEBEV = Y*D-DD + 0.5D0 * C(1)      !LAST STEP IS DIFFERENT

      END FUNCTION CHEBEV
*****************************************************************


****************************************************************
      SUBROUTINE OPTIC_DATA(
     $     NWS,
     $     DV2,      DV3,       V2S,      U0,
     $     CST_O2,   CST_O3,    ALBEDO,
     $     TAUS_CLR, TAUA_CLR,  TAUS_AER, TAUA_AER,
     $     TAUS_CLD, TAUA_CLD,
     $     PRAY,     PCLD,      PAER,
     $     FLX,      ALB)

! Description:
c   this subroutine links MISTRA with Jochen Landgraf's code.
c   it transfers the optical data from MISTRA to the relevant quantities
c   in Jochen's code. Lower case variables are from MISTRA and are passed in common
c   blocks. Upper case variables are from Jochen's code and are passed explicitly.

! 13/08/2016 jjb
!     changes done in this SR:
!        - removal of 8 unused arguments
!        - USE global_params for vertical grid parameters
!        - several unused variables commented
!        - /RAY_J/ CS_RAY instead of CRRAY (old name ?) for homogeneity
!
! 17/08/2016 & 23/08/2016 jjb
!     further changes:
!        - major cleaning of all unused, comented parts except potentially useful comments
!
! 17/05/2021 jjb
!     add missing initialisation of TOTABS and TOTSCA, but commented since unused
!     missing declarations then implicit none

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      IMPLICIT NONE

      INTEGER I,K,L,MAXLAY,MAXWAV,NMOM,NW,NL
      PARAMETER(MAXLAY=nrlay, MAXWAV=176, NMOM =4, NW=7)

      INTEGER
     $     NWS(NW)   !specification of interval: 1 < NWS(L) < MAXWAV

      DOUBLE PRECISION
     $     V2S(0:MAXLAY),       !O2 column density [part./cm^2]
     $     DV2(MAXLAY),         !diff. O2 column density [part./cm^2]
     $     DV3(MAXLAY),         !diff. O3 column density [part./cm^2]
     $     U0                  !cosine of solar zenith angle

      DOUBLE PRECISION
c                                  !per layer [part./cm^2]
     $     ALBEDO(MAXWAV)       !shortwave albedo

      DOUBLE PRECISION
     $     TAUS_CLR(NW,MAXLAY),  !scattering optical depth for clear sky
     $     TAUA_CLR(NW,MAXLAY),  !absorption optical depth for clear sky
     $     TAUS_AER(NW,MAXLAY),  !scattering optical depth for aerosol
     $     TAUA_AER(NW,MAXLAY),  !absorption optical depth for aerosol
     $     TAUS_CLD(NW),  !scattering optical depth for cloudy sky
     $     TAUA_CLD(NW)   !absorption optical depth for cloudy sky

      DOUBLE PRECISION
     $     PRAY(NW,MAXLAY,0:NMOM), !Rayleigh phase fct. moments
     $     PCLD(NW,0:NMOM),        !phase fct. moments for clouds
     $     PAER(NW,MAXLAY,0:NMOM)  !phase fct. moments for aerosols

      DOUBLE PRECISION
     $     ALB(NW),                      !shortwave albedo
     $     FLX(NW)                          !extraterrestic flux

      DOUBLE PRECISION
     $     CST_O3(0:MAXLAY,MAXWAV),     CST_O2(MAXWAV,0:MAXLAY)

c     coefficients for optical depth of Schumann-Runge band above TOA
C     asume fixed temperature of T=220K. Total error in JO2 and FINT of the
c     fit is less than 4% for 45.5 .le. DLOG(V2S) .le. 54.0.

      COMMON/C_O2_TOP/ CT_TOP(4,13)
      DOUBLE PRECISION CT_TOP

      COMMON/RAY_J/ CS_RAY(MAXWAV)
      DOUBLE PRECISION CS_RAY
      COMMON/FL/   FLUX(MAXWAV)   !extraterrestic flux per interval
c                                  [photons/(cm^2 s)]
      DOUBLE PRECISION FLUX

      DOUBLE PRECISION
     $     DLV2S,
     $     TA_O2,
     $     TA_O3

!      DOUBLE PRECISION
!     $     Ti_SCA(NW,MAXLAY),  !integrated absorption optical depth
!     $     Ti_ABS(NW,MAXLAY),  !integrated absorption optical depth
!     $     TOTABS(NW),         !total absorption optical depth
!     $     TOTSCA(NW)          !total absorption optical depth


c MISTRA common block
      common /extra_2/ taer_s(nrlay),taer_a(nrlay),ga_pl(nrlay)
      double precision taer_s, taer_a, ga_pl
C-------------------------------------------------------------------
C     internal functions

      DOUBLE PRECISION P3
      DOUBLE PRECISION C0, C1, C2, C3, X
      P3(C0,C1,C2,C3,X) = C0 + (C1 + (C2+C3*X)*X)*X

C-------------------------------------------------------------------

      IF (1.D0/U0 .ge. 0.D0) THEN

      !! jjb initialisation was missing
      !TOTABS(:) = 0.D0
      !TOTSCA(:) = 0.D0

      DO NL = 1,NW
         L = NWS(NL)

         DO  K=1,MAXLAY

c           OPTICAL DEPTHS

            TA_O2 = 0.5D0 * (CST_O2(L,K-1) + CST_O2(L,K)) * DV2(K)

            IF (K.eq.1 .AND. L.le.13) THEN
               DLV2S = DLOG(V2S(1))
               TA_O2 = U0*EXP(P3(CT_TOP(1,L),CT_TOP(2,L),CT_TOP(3,L),
     $                              CT_TOP(4,L),DLV2S))
            ENDIF

c crray method is the correct one
            TA_O3 = 0.5D0*(CST_O3(K-1,L) + CST_O3(K,L)) * DV3(K)

            TAUA_CLR(NL,K) = TA_O2 + TA_O3

c           in the Schumann-Runge band only absorption

            IF(L.le.13) THEN
               TAUS_CLR(NL,K) = 0.D0
            ELSE
               TAUS_CLR(NL,K) = CS_RAY(L)*1./0.21*DV2(K)
            ENDIF

            TAUA_AER(NL,K) = taer_a(K)
            TAUS_AER(NL,K) = taer_s(K)

            TAUA_CLD(NL) = 0.D0     !no clouds
            TAUS_CLD(NL) = 0.D0

C           PHASE FUNCTIONS

            PRAY(NL,K,0) = 1.D0
            PCLD(NL,0) = 1.D0
            PAER(NL,K,0) = 1.D0

            DO I = 1,NMOM
               PRAY(NL,K,I) = 0.D0
               PCLD(NL,I) = 0.D0   !no clouds
               PAER(NL,K,I) = ga_pl(K)**I    !Henyey-Greenstein Fct
            ENDDO


            PRAY(NL,K,2) = 0.1D0

            ! jjb commented out, unused
            !TOTABS(NL) = TOTABS(NL) + TAUA_CLR(NL,K)
            !TOTSCA(NL) = TOTSCA(NL) + TAUS_CLR(NL,K)

            !IF(K.eq.1)THEN
            !   Ti_SCA(NL,K)=TAUS_CLR(NL,K)
            !   Ti_ABS(NL,K)=TAUA_CLR(NL,K)
            !ELSE
            !   Ti_SCA(NL,K)=Ti_SCA(NL,K-1)+TAUS_CLR(NL,K)
            !   Ti_ABS(NL,K)=Ti_ABS(NL,K-1)+TAUA_CLR(NL,K)
            !ENDIF

         ENDDO
      ENDDO
      ENDIF


      DO NL = 1,NW
         FLX(NL) = FLUX(NWS(NL))
         ALB(NL) = ALBEDO(NWS(NL))
      ENDDO

      END SUBROUTINE OPTIC_DATA
***************************************************************************


***************************************************************************
      SUBROUTINE FOUR_INTF(
     $     TAUS_CLR, TAUA_CLR, TAUS_AER,  TAUA_AER,
     $     TAUS_CLD, TAUA_CLD,
     $     PRAY,     PAER,     PCLD,
     $     ALB,      FLX,       U0,
     $     FACT,     SW13)

C****************************************************************
C  run analytic FOUR-stream code to calculate actinic fluxes    C
C****************************************************************

! 13/08/2016 jjb comments before release
!     changes done in this SR:
!        - removal of one unused argument
!        - USE global_params for vertical grid parameters
!        - removal of unused parameters
!        - commented unused variables


      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      IMPLICIT NONE

      INTEGER MAXLAY,NMOM,NW,SW13
      PARAMETER(MAXLAY=nrlay, NMOM =4,NW=7)

c     inputs

      DOUBLE PRECISION
     $     TAUS_CLR(NW,MAXLAY),    !scattering optical depth for clear sky
     $     TAUA_CLR(NW,MAXLAY),    !absorption optical depth for clear sky
     $     TAUS_AER(NW,MAXLAY),    !scattering optical depth for aerosol
     $     TAUA_AER(NW,MAXLAY),    !absorption optical depth for aerosol
     $     TAUS_CLD(NW),    !scattering optical depth for cloudy sky
     $     TAUA_CLD(NW)     !absorption optical depth for cloudy sky

      DOUBLE PRECISION
     $     PRAY(NW,MAXLAY,0:NMOM), !Rayleigh phase fct. moments
     $     PCLD(NW,0:NMOM),        !phase fct. moments for clouds
     $     PAER(NW,MAXLAY,0:NMOM)  !phase fct. moments for aerosols

      DOUBLE PRECISION
     $     ALB(NW),                !shortwave albedo
     $     FLX(NW),                   !extraterrestic flux
     $     U0                     !cosine of solar zenith angle

      DOUBLE PRECISION
     $     FACT(MAXLAY+1,NW)      !actinic flux

      DOUBLE PRECISION
     $     TC(MAXLAY),               !ext. opt. depth for each layer
     $     TT(MAXLAY),               !ext. opt. depth (counted from TOA)
     $     WC(MAXLAY),               !single scat. albedo (not delta-scaled)
     $     WW1(MAXLAY),              !legendre
     $     WW2(MAXLAY),              !...phase
     $     WW3(MAXLAY),              !......funct.
     $     WW4(MAXLAY)               !.........exp. coeff. (not scaled)

      DOUBLE PRECISION
     $     FSD(MAXLAY+1),            !parallel solar flux (delta scaled)
     $     FFD(MAXLAY+1),            !downward diffuse flux
     $     FFU(MAXLAY+1),            !upward diffuse flux
     $     UAV(MAXLAY+1)             !actinic flux

      INTEGER I, K, L, NRFL, NP
      DOUBLE PRECISION AS, FF0, PI, TAUSCAT, TAUTOT, UU0, XX

C-------------------------------------------------------------------
      PI=2.D0*DASIN(1.0D0)

      DO 100 L  = 1,NW            !wavel. loop

         IF (U0.ge.0.) THEN
         AS=ALB(L)
         UU0=U0
         FF0=FLX(L)/PI

         DO K=1,MAXLAY
            TAUSCAT  = TAUS_CLR(L,K) + TAUS_CLD(L) +
     $                 TAUS_AER(L,K)
            TAUTOT   = TAUA_CLR(L,K) + TAUA_CLD(L) +
     $                 TAUA_AER(L,K) + TAUSCAT
            TC(K) = TAUTOT
            IF(TAUTOT .LT. 1.D-20) THEN
               WC(K)=1.D0
            ELSE
               WC(K)=TAUSCAT/TAUTOT
            ENDIF

C  setup of phase function expansion coefficients

            IF (TAUSCAT .LT. 1.D-20 ) THEN
               WW1(K) = 0.D0
               WW2(K) = 0.D0
               WW3(K) = 0.D0
               WW4(K) = 0.D0
            ELSE
            I = 1
            XX = DBLE(2*I+1)
               WW1(K) = ( XX*PCLD(L,I)*TAUS_CLD(L) +
     $                    XX*PRAY(L,K,I)*TAUS_CLR(L,K) +
     $                    XX*PAER(L,K,I)*TAUS_AER(L,K) ) /
     $                    TAUSCAT
            I = 2
            XX = DBLE(2*I+1)
               WW2(K) = ( XX*PCLD(L,I)*TAUS_CLD(L) +
     $                    XX*PRAY(L,K,I)*TAUS_CLR(L,K) +
     $                    XX*PAER(L,K,I)*TAUS_AER(L,K) ) /
     $                    TAUSCAT
            I = 3
            XX = DBLE(2*I+1)
               WW3(K) = ( XX*PCLD(L,I)*TAUS_CLD(L) +
     $                    XX*PRAY(L,K,I)*TAUS_CLR(L,K) +
     $                    XX*PAER(L,K,I)*TAUS_AER(L,K) ) /
     $                    TAUSCAT
            I = 4
            XX = DBLE(2*I+1)
               WW4(K) = ( XX*PCLD(L,I)*TAUS_CLD(L) +
     $                    XX*PRAY(L,K,I)*TAUS_CLR(L,K) +
     $                    XX*PAER(L,K,I)*TAUS_AER(L,K) ) /
     $                    TAUSCAT
            ENDIF

         ENDDO

c  setup total extinction optical depth (counted from TOA)

         TT(1)=TC(1)
         DO K=2,MAXLAY
            TT(K)=TT(K-1)+TC(K)
         ENDDO

         NRFL = MAXLAY
         NP = MAXLAY+1

         DO K=1,NP
            FFU(K)=0.D0
            FFD(K)=0.D0
            FSD(K)=0.D0
            UAV(K)=0.D0
         ENDDO

         IF(SW13.EQ.1) CALL CHECK_FOUR(
     $            TT,     WC,
     $            WW1,    WW2,   WW3,   WW4,
     $            UU0,    AS,    FF0)


         CALL QFTS(NP, NRFL,  AS, UU0, FF0,
     $        WW1, WW2, WW3, WW4,  WC,  TT,
     $        FFU, FFD, FSD, UAV)

c     actinic flux

         DO K=1,MAXLAY+1
            FACT(K,L) = MAX(0.d0 , 4.D0*PI*UAV(K))
         ENDDO
      ELSE

         DO K=1,MAXLAY+1
            FACT(K,L) = 0.0D0
         ENDDO

      ENDIF
 100  CONTINUE

! jjb: it the following block is uncommented, make sure that the unit number is not already open
c      OPEN(22,FILE='result/fact_four.dat',STATUS='UNKNOWN')

c      DO J=1,MJ
c      DO L = 1,NW
c         WRITE(22,'(1P,6E12.4)')(FACT(K,L),K=1,MAXLAY+1)
c      ENDDO
c      ENDDO
c      CLOSE(22)

      END SUBROUTINE FOUR_INTF
*************************************************************


*************************************************************
      SUBROUTINE CHECK_FOUR(
     $            TT,     WC,
     $            WW1,    WW2,   WW3,   WW4,
     $            U0,     ALB,   FLX)

! 13/08/2016 jjb comments before release
!     changes done in this SR:
!        - USE global_params for vertical grid parameters
!        - little cleaning


      USE file_unit, ONLY :
! Imported Parameters:
     &     jpfunchk4

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      USE precision, ONLY :
! Imported Parameters:
     &     dp

      IMPLICIT NONE

      INTEGER MAXLAY
      PARAMETER(MAXLAY=nrlay)

      REAL (KIND=DP), INTENT(IN) ::
     $     TT(MAXLAY),
     $     WC(MAXLAY),
     $     WW1(MAXLAY),
     $     WW2(MAXLAY),
     $     WW3(MAXLAY),
     $     WW4(MAXLAY)
      REAL (KIND=DP), INTENT(IN) :: U0, ALB, FLX

      REAL (KIND=DP), PARAMETER :: PI = 2._dp * ASIN(1._dp)

      INTEGER :: K

C-----------------------------------------------------------
      OPEN(jpfunchk4,FILE='check_four.out',STATUS='UNKNOWN')

      WRITE(jpfunchk4,*)'total extinction optical depth'
      WRITE(jpfunchk4,'(1P,6E12.4)')(TT(K),K=1,MAXLAY)

      WRITE(jpfunchk4,*)'single scattering albedo (not delta scaled)'
      WRITE(jpfunchk4,'(1P,6E12.4)')(WC(K),K=1,MAXLAY)

      WRITE(jpfunchk4,*)'phase fn. exp. coeff  I=1 (unscaled)'
      WRITE(jpfunchk4,'(1P,6E12.4)')(WW1(K),K=1,MAXLAY)

      WRITE(jpfunchk4,*)'phase fn. exp. coeff  I=2 (unscaled)'
      WRITE(jpfunchk4,'(1P,6E12.4)')(WW2(K),K=1,MAXLAY)

      WRITE(jpfunchk4,*)'phase fn. exp. coeff  I=3 (unscaled)'
      WRITE(jpfunchk4,'(1P,6E12.4)')(WW3(K),K=1,MAXLAY)

      WRITE(jpfunchk4,*)'phase fn. exp. coeff  I=4 (unscaled)'
      WRITE(jpfunchk4,'(1P,6E12.4)')(WW4(K),K=1,MAXLAY)

      WRITE(jpfunchk4,*)'solar cosine'
      WRITE(jpfunchk4,'(1P,6E12.4)') U0

      WRITE(jpfunchk4,*)'surface albedo'
      WRITE(jpfunchk4,'(1P,6E12.4)') ALB

      WRITE(jpfunchk4,*)'extraterrestrial flux'
      WRITE(jpfunchk4,'(1P,6E12.4)') FLX*PI

      CLOSE(jpfunchk4)

      END SUBROUTINE CHECK_FOUR
*************************************************************


************************************************************************
*                 ANALYTIC FOUR STREAM METHOD                          *
*                 to calculate actinic fluxes                          *
************************************************************************
*     This version is not suitable for the calculation of radiative    *
*     transfer through layers with partial cloud covers!               *
************************************************************************


c **********************************************************************
c The delta-four-stream approximation for nonhomogeneous atmospheres
c in the solar wavelengths (Fu, 1991). The input parameters are ndfs,
c mdfs, and ndfs4 through 'para.file', as, u0, f0 for solar through
c arguments of  'qfts', and
c ww1(ndfs), ww2(ndfs), ww3(ndfs), ww4(ndfs), ww(ndfs), and tt(ndfs)
c through common statement 'dfsin'.
c **********************************************************************
      subroutine qfts(np,  nrfl,   as,  u0,   f0,
     I               ww1,  ww2,  ww3,  ww4,  ww,  tt,
     O               ffu,  ffd,  fsd,  uav)

! 13/08/2016 jjb: final changes before release
!        - USE global_params for vertical grid parameters
!        - removal of unused parameters and common blocks
!        - cleaning, more DP declaration
! 25/06/2021 jjb: all missing declarations and implicit none


      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      implicit none

      INTEGER MAXLAY,ndfs,mdfs
      PARAMETER (MAXLAY=nrlay)
      parameter (ndfs  = MAXLAY,
     $           mdfs  = MAXLAY + 1)

c input
      double precision ww1(ndfs), ww2(ndfs), ww3(ndfs), ww4(ndfs),
     1          ww(ndfs), tt(ndfs)
      double precision w1(ndfs), w2(ndfs), w3(ndfs), w(ndfs),
     1                 t(ndfs), u0a(ndfs), f0a(ndfs)
      double precision fk1(ndfs), fk2(ndfs), a4(4,4,ndfs),
     1                 z4(4,ndfs), g4(4,ndfs)
      double precision as
c output
      double precision ffu(mdfs), ffd(mdfs), fsd(mdfs), uav(mdfs)
c local
      double precision x(4), fi(4)
      integer i, ii, jj, k, n, np, nrfl, m
      double precision asbs, f0, fourpi, fw1, fw2, fw3, fw4,
     $     pi, u0, y, y1

      n = nrfl
      m = np
      asbs = as
      call adjust(nrfl, ww1,  ww2,  ww3,  ww4,   ww,   tt,
     O            w1,    w2,   w3,    w,    t)
      do 5 i = 1, n
         u0a(i) = u0
         f0a(i) = f0
 5    continue
      call qccfe(nrfl, asbs,  w1,  w2,  w3,   w,   t, u0a, f0a,
     O           fk1,   fk2,  a4,  z4,  g4)
!      pi = 3.14159265   ! jjb changed:
      PI=2.*DASIN(1.0D0) !     same definition of pi as in two other places in jrate.f
      fw1 = 0.6638961D0
      fw2 = 2.4776962D0
      fw3 = u0 * pi * f0
      fourpi = 4.0D0 * pi
      do 10 i = 1, m
         if ( i .eq. 1 ) then
            x(1) = 1.0D0
            x(2) = 1.0D0
            x(3) = exp ( - fk1(1) * t(1) )
            x(4) = exp ( - fk2(1) * t(1) )
            k = 1
            y = 1.0D0
         elseif ( i .eq. 2 ) then
            x(1) = exp ( - fk2(1) * t(1) )
            x(2) = exp ( - fk1(1) * t(1) )
            x(3) = 1.0D0
            x(4) = 1.0D0
            k = 1
            y = exp ( - t(1) / u0 )
         else
            k = i - 1
            y1 = t(k) - t(k-1)
            x(1) = exp ( - fk2(k) * y1 )
            x(2) = exp ( - fk1(k) * y1 )
            x(3) = 1.0D0
            x(4) = 1.0D0
            y = exp ( - t(k) / u0 )
         endif
         do 37 jj = 1, 4
            fi(jj) = z4(jj,k) * y
 37      continue
         do 40 ii = 1, 4
            fw4 = g4(ii,k) * x(ii)
            do 45 jj = 1, 4
               fi(jj) = fi(jj) + a4(jj,ii,k) * fw4
 45         continue
 40      continue
         ffu(i)= fw1 * fi(2) + fw2 * fi(1)
         ffd(i)= fw1 * fi(3) + fw2 * fi(4)
         fsd(i)= fw3 * y
         uav(i)= 0.25D0*(fi(1)+fi(2)+fi(3)+fi(4))+y*f0*pi/fourpi
 10   continue

      end subroutine qfts
*******************************************************************


c **********************************************************************
c In this subroutine, we incorporate a delta-function adjustment to
c account for the  forward  diffraction  peak in the context of the
c four-stream approximation ( Liou, Fu and Ackerman, 1988 ). w1(n),
c w2(n), w3(n), w(n), and t(n) are the adjusted parameters.
c **********************************************************************
      subroutine adjust(nrfl, ww1,  ww2,  ww3,  ww4,
     I                  ww,    tt,
     O                  w1,    w2,   w3,    w,    t)

! 13/08/2016 jjb: final changes before release
!        - USE global_params for vertical grid parameters
!        - removal of unused parameters
!        - cleaning
! 20/08/2016 jjb: missing declarations and implicit none


      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      implicit none

      INTEGER MAXLAY,ndfs
      PARAMETER (MAXLAY=nrlay)
      parameter (ndfs  = MAXLAY)

c input
      double precision ww1(ndfs), ww2(ndfs), ww3(ndfs), ww4(ndfs),
     1          ww(ndfs), tt(ndfs)
      integer nrfl
c output
      double precision w1(ndfs), w2(ndfs), w3(ndfs), w(ndfs),
     1          t(ndfs)
c local
      double precision dtt(ndfs), dt(ndfs), f, fw, tt0

      integer i,n
      n = nrfl
      tt0 = 0.0D0
      do 10 i = 1, n
         f = ww4(i) / 9.0D0
         fw = 1.0D0 - f * ww(i)
         w1(i) = ( ww1(i) - 3.0D0 * f ) / ( 1.0D0 - f )
         w2(i) = ( ww2(i) - 5.0D0 * f ) / ( 1.0D0 - f )
         w3(i) = ( ww3(i) - 7.0D0 * f ) / ( 1.0D0 - f )
         w(i) = ( 1.0D0 - f ) * ww(i) / fw
         dtt(i) = tt(i) - tt0
         tt0 = tt(i)
         dt(i) = dtt(i) * fw
 10   continue
      t(1) = dt(1)
      do 20 i = 2, n
         t(i) = dt(i) + t(i-1)
 20   continue

      end subroutine adjust
*******************************************************************


c **********************************************************************
c In the solar band  asbs is the surface albedo, while in the infrared
c band asbs is  blackbody intensity emitted at the surface temperature
c times surface emissivity.  In this subroutine, the delta-four-stream
c is applied to nonhomogeneous atmospheres. See comments in subroutine
c 'qcfel' for array AB(13,4*n).
c **********************************************************************
      subroutine qccfe(nrfl, asbs,  w1,  w2,  w3,   w,   t,  u0,  f0,
     O                 fk1,   fk2,  a4,  z4,  g4)

! 13/08/2016 jjb: final changes before release
!        - USE global_params for vertical grid parameters
!        - removal of unused parameters and common blocks,
!        - updating calls after changes of argument lists
!        - cleaning, more DP declaration


      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      implicit none

      INTEGER MAXLAY,ndfs,ndfs4
      PARAMETER (MAXLAY=nrlay)
      parameter (ndfs  = MAXLAY,
     $           ndfs4 = 4 * MAXLAY)

      integer nrfl
      double precision asbs

      double precision w1(ndfs), w2(ndfs), w3(ndfs), w(ndfs),
     1          t(ndfs), u0(ndfs), f0(ndfs)

      double precision fk1(ndfs), fk2(ndfs), a4(4,4,ndfs),
     1          z4(4,ndfs), g4(4,ndfs)

      double precision aa(4,4,2), zz(4,2), a1(4,4), z1(4)

      double precision b(4,3)
      double precision af2(2,2,2), d(4)
      double precision z(4)
      double precision ab(13,ndfs4), bx(ndfs4), xx(ndfs4)
c local
      double precision fu(4,4), wu(4)
      double precision b1, c1, f0n, fk1t, fk2t, fw1, fw2, t0n, t1n,
     $     u0n, v1, v2, v3, w1n, w2n, w3n, wn
      integer i, i1, i2, i3, i8, j, j1, j2, j3, k, kf,
     $     m1, m2, m18, m28, n, n4

      n = nrfl
      n4 = n*4
      do 333 i = 1, n4
         do 333 j = 1, 13
            ab(j,i) = 0.0D0
 333     continue
      wn = w(1)
      w1n = w1(1)
      w2n = w2(1)
      w3n = w3(1)
      t0n = 0.0
      t1n = t(1)
      u0n = u0(1)
      f0n = f0(1)
      if ( wn .ge. 0.99999999999D0 ) then
         wn = 0.99999999999D0
      endif
      if ( wn .le. 1.0d-12 ) then
         call coefft0(t0n,  t1n,
     O                z1,  fk1t, fk2t,   a1,   zz,   aa)
         fk1(1) = fk1t
         fk2(1) = fk2t
      else
         call coeff1(wn,  w1n,  w2n,  w3n,  u0n,
     O               b)
         call coeff2(u0n,   b,
     O               af2,   d)
         call coeff4(u0n, af2,   d,
     O               b1,   c1,   z)
         call coeffl(t0n, t1n, u0n, f0n,
     I               b,   af2,   b1,   c1,    z,
     O               z1, fk1t, fk2t,   a1,   zz,   aa)
         fk1(1) = fk1t
         fk2(1) = fk2t
      endif
      do 10 i = 1, 4
         z4(i,1) = z1(i)
         do 10 j = 1, 4
            a4(i,j,1) = a1(i,j)
 10      continue
      do 20 i = 1, 2
         bx(i) = - zz(i+2,1)
         i8 = i + 8
         do 20 j = 1, 4
            ab(i8-j,j) = aa(i+2,j,1)
 20      continue
      do 30 i = 1, 4
         wu(i) = zz(i,2)
         do 30 j = 1, 4
            fu(i,j) = aa(i,j,2)
 30      continue
      do 40 k = 2, n
         wn = w(k)
         w1n = w1(k)
         w2n = w2(k)
         w3n = w3(k)
         t0n = t(k-1)
         t1n = t(k)
         u0n = u0(k)
         f0n = f0(k)
         if ( wn .ge. 0.99999999999D0 ) then
            wn = 0.99999999999D0
         endif
         if ( wn .le. 1.0d-12 ) then
            call coefft0(t0n,  t1n,
     O                   z1,  fk1t, fk2t,   a1,   zz,   aa)
            fk1(k) = fk1t
            fk2(k) = fk2t
         else
            call coeff1(wn,  w1n,  w2n,  w3n,  u0n,
     O                  b)
            call coeff2(u0n,   b,
     O                  af2,   d)
            call coeff4(u0n, af2,  d,
     O                  b1,   c1,  z)
            call coeffl(t0n, t1n, u0n, f0n,
     I                  b,   af2,   b1,   c1,    z,
     O                  z1, fk1t, fk2t,   a1,   zz,   aa)
            fk1(k) = fk1t
            fk2(k) = fk2t
         endif
         do 50 i = 1, 4
            z4(i,k) = z1(i)
            do 50 j = 1, 4
               a4(i,j,k) = a1(i,j)
 50         continue
         kf = k + k + k + k
         i1 = kf - 5
         i2 = i1 + 3
         j1 = kf - 7
         j2 = j1 + 3
         i3 = 0
         do 55 i = i1, i2
            i3 = i3 + 1
            bx(i) = - wu(i3) + zz(i3,1)
            j3 = 0
            i8 = i + 8
            do 60 j = j1, j2
               j3 = j3 + 1
               ab(i8-j,j) = fu(i3,j3)
 60         continue
            j3 = 0
            do 65 j = j2 + 1, j2 + 4
               j3 = j3 + 1
               ab(i8-j,j) = - aa(i3,j3,1)
 65         continue
 55      continue
         do 70 i = 1, 4
            wu(i) = zz(i,2)
            do 70 j = 1, 4
               fu(i,j) = aa(i,j,2)
 70         continue
 40   continue
      v1 = 0.2113247D0 * asbs
      v2 = 0.7886753D0 * asbs
      v3 = asbs * u0(1) * f0(1) * exp ( - t(n) / u0(1) )
      m1 = n4 - 1
      m2 = n4
      m18 = m1 + 8
      m28 = m2 + 8
      fw1 = v1 * wu(3)
      fw2 = v2 * wu(4)
      bx(m1) = - ( wu(1) - fw1 - fw2 - v3 )
      bx(m2) = - ( wu(2) - fw1 - fw2 - v3 )
      do 80 j = 1, 4
         j1 = n4 - 4 + j
         fw1 = v1 * fu(3,j)
         fw2 = v2 * fu(4,j)
         ab(m18-j1,j1) = fu(1,j) - fw1 - fw2
         ab(m28-j1,j1) = fu(2,j) - fw1 - fw2
 80   continue

      call qcfel(nrfl,  ab,   bx,   xx)
      do 90 k = 1, n
         j = k + k + k + k - 4
         do 90 i = 1, 4
            j = j + 1
            g4(i,k) = xx(j)
 90      continue

      end subroutine qccfe

c **********************************************************************
c Double-Gauss quadratures and weights (Sykes, 1951).
c **********************************************************************
      block data gaus2
      implicit none
      common /point/ u(4)
      double precision u
      data u / -0.7886752D0, -0.2113247D0, 0.2113247D0, 0.7886752D0 /
      end

c *********************************************************************
c p0, p1, p2 and p3 are Legendre polynomials for l = 1, 2, 3.
c *********************************************************************
c      function p0 ( x )
c      p0 = 1.0
c      return
c      end
c      function p1 ( x )
c      p1 = x
c      return
c      end
c      function p2 ( x )
c      p2 = 1.5 * x * x - 0.5
c      return
c      end
c      function p3 ( x )
c      p3 = ( 2.5 * x * x - 1.5 ) * x
c      return
c      end

c **********************************************************************
c p0d(4), p1d(4), p2d(4), and p3d(4) are Legendre polynomials p0(x),
c p1(x), p2(x), and p3(x) when x = u(1), u(2), u(3), and u(4).
c **********************************************************************
      block data legend
      implicit none
      common /legen/ p0d(4), p1d(4), p2d(4), p3d(4)
      double precision p0d, p1d, p2d, p3d
      data p0d /  .100000D+01,  .100000D+01,  .100000D+01, .100000D+01 /
      data p1d / -.788675D+00, -.211325D+00,  .211325D+00, .788675D+00 /
      data p2d /  .433013D+00, -.433013D+00, -.433013D+00, .433013D+00 /
      data p3d / -.433940D-01,  .293394D+00, -.293394D+00, .433940D-01 /
      end

c *********************************************************************
c p11d(4,4), p22d(4,4), and p33d(4,4) are defined as 0.5*p1d(i)*p1d(j),
c 0.5*p2d(i)*p2d(j), and 0.5*p3d(i)*p3d(j), respectively.
c *********************************************************************
      block data legenf
      implicit none
      common /legen1/ p11d(4,4), p22d(4,4), p33d(4,4)
      double precision p11d, p22d, p33d
      data p11d / .311004D+00, .833334D-01,-.833334D-01,-.311004D+00,
     1            .833334D-01, .223291D-01,-.223291D-01,-.833334D-01,
     1           -.833334D-01,-.223291D-01, .223291D-01, .833334D-01,
     1           -.311004D+00,-.833334D-01, .833334D-01, .311004D+00 /
      data p22d / .937501D-01,-.937501D-01,-.937501D-01, .937501D-01,
     1           -.937501D-01, .937501D-01, .937501D-01,-.937501D-01,
     1           -.937501D-01, .937501D-01, .937501D-01,-.937501D-01,
     1            .937501D-01,-.937501D-01,-.937501D-01, .937501D-01 /
      data p33d / .941520D-03,-.636577D-02, .636577D-02,-.941520D-03,
     1           -.636577D-02, .430400D-01,-.430400D-01, .636577D-02,
     1            .636577D-02,-.430400D-01, .430400D-01,-.636577D-02,
     1           -.941520D-03, .636577D-02,-.636577D-02, .941520D-03 /
      end

c **********************************************************************
c coefficient calculations for four first-order differential equations.
c **********************************************************************
      subroutine coeff1(w,   w1,   w2,   w3,   u0,  b)

! 13/08/2016 jjb: final changes before release
!                 removal of unused parameters and common blocks
!                 cleaning

      implicit none

      common /point/ u(4)
      double precision u
      common /legen/ p0d(4), p1d(4), p2d(4), p3d(4)
      double precision p0d, p1d, p2d, p3d
      common /legen1/ p11d(4,4), p22d(4,4), p33d(4,4)
      double precision p11d, p22d, p33d

c input
      double precision, intent(in) :: w, w1, w2, w3, u0

c output
      double precision b(4,3)

c local
      integer i, j
      double precision c(4,5)
      double precision fq, q1, q2, q3, fw, w0w, w1w, w2w, w3w, x

c-----------------------------------------------------------------------
      x = 0.5D0 * w
      w0w = x
      w1w = x * w1
      w2w = x * w2
      w3w = x * w3
      fw = u0 * u0
      q1 = - w1w * u0
      q2 = w2w * ( 1.5D0 * fw - 0.5D0 )
      q3 = - w3w * ( 2.5D0 * fw - 1.5D0 ) * u0
      fq = 0.5D0 * w0w
      do 10 i = 3, 4
         do 20 j = 1, 4
            c(i,j) = fq + w1w * p11d(i,j) +
     1               w2w * p22d(i,j) + w3w * p33d(i,j)
            if ( i .eq. j ) then
               c(i,j) = ( c(i,j) - 1.0D0 ) / u(i)
            else
               c(i,j) = c(i,j) / u(i)
            endif
 20      continue
 10   continue
      do 30 i = 1, 4
         c(i,5) = w0w + q1 * p1d(i) +
     1            q2 * p2d(i) + q3 * p3d(i)
         c(i,5) = c(i,5) / u(i)
30       continue
      b(1,1) = c(4,4) - c(4,1)
      b(1,2) = c(4,4) + c(4,1)
      b(2,1) = c(4,3) - c(4,2)
      b(2,2) = c(4,3) + c(4,2)
      b(3,1) = c(3,4) - c(3,1)
      b(3,2) = c(3,4) + c(3,1)
      b(4,1) = c(3,3) - c(3,2)
      b(4,2) = c(3,3) + c(3,2)
      b(1,3) = c(4,5) - c(1,5)
      b(2,3) = c(3,5) - c(2,5)
      b(3,3) = c(3,5) + c(2,5)
      b(4,3) = c(4,5) + c(1,5)

      end subroutine coeff1

c **********************************************************************
c coefficient calculations for second order differential equations.
c **********************************************************************
      subroutine coeff2(u0,  b,
     O                   a,   d)
      implicit none
c input
      double precision, intent(in) :: u0, b(4,3)
c output
      double precision, intent(out) :: a(2,2,2), d(4)
      double precision fw1, fw2, fw3, fw4
      fw1 = b(1,1) * b(1,2)
      fw2 = b(2,1) * b(3,2)
      fw3 = b(3,1) * b(2,2)
      fw4 = b(4,1) * b(4,2)
      a(2,2,1) = fw1 + fw2
      a(2,1,1) = b(1,1) * b(2,2) + b(2,1) * b(4,2)
      a(1,2,1) = b(3,1) * b(1,2) + b(4,1) * b(3,2)
      a(1,1,1) = fw3 + fw4
      a(2,2,2) = fw1 + fw3
      a(2,1,2) = b(1,2) * b(2,1) + b(2,2) * b(4,1)
      a(1,2,2) = b(3,2) * b(1,1) + b(4,2) * b(3,1)
      a(1,1,2) = fw2 + fw4
      d(1) = b(3,2) * b(4,3) + b(4,2) * b(3,3) + b(2,3) / u0
      d(2) = b(1,2) * b(4,3) + b(2,2) * b(3,3) + b(1,3) / u0
      d(3) = b(3,1) * b(1,3) + b(4,1) * b(2,3) + b(3,3) / u0
      d(4) = b(1,1) * b(1,3) + b(2,1) * b(2,3) + b(4,3) / u0

      end subroutine coeff2

c **********************************************************************
c coefficient calculations for fourth-order differential equations.
c **********************************************************************
      subroutine coeff4(u0,   a,    d,
     O                  b1,  c1,    z)
      implicit none
c input
      double precision, intent(in) :: u0, a(2,2,2), d(4)
c output
      double precision, intent(out) :: b1, c1, z(4)
      double precision x
      x = u0 * u0
      b1 = a(2,2,1) + a(1,1,1)
      c1 = a(2,1,1) * a(1,2,1) - a(1,1,1) * a(2,2,1)
      z(1) = a(2,1,1) * d(3) + d(4) / x - a(1,1,1) * d(4)
      z(2) = a(1,2,1) * d(4) - a(2,2,1) *d(3) + d(3) / x
      z(3) = a(2,1,2) * d(1) + d(2) / x - a(1,1,2) * d(2)
      z(4) = a(1,2,2) * d(2) - a(2,2,2) * d(1) + d(1) / x
      end subroutine coeff4

c **********************************************************************
c fk1 and fk2 are the eigenvalues.
c **********************************************************************

      subroutine coeffl(t0,  t1,  u0,  f0,
     I                  b,    a,  b1,  c1,   z,
     O                  z1, fk1, fk2,  a1,  zz,  aa)

! 13/08/2016 jjb: final changes before release
!                 removal of one unused argument,
!                 removal of unused parameters
!                 cleaning.

      implicit none

c input
      double precision b(4,3)
      double precision a(2,2,2)
      double precision z(4)
      double precision b1, c1, f0, t0, t1, u0

c output
      double precision aa(4,4,2), zz(4,2), a1(4,4), z1(4)
      double precision fk1, fk2

c local
      integer i
      double precision a2, b2, dt, fq0, fq1, fw, fw1, fw2, x, y, zx
c--------------------------------------------------------------
      dt = t1 - t0
      x = sqrt ( b1 * b1 + 4.0D0 * c1 )
      fk1 = sqrt ( ( b1 + x ) * 0.5D0 )
      fk2 = sqrt ( ( b1 - x ) * 0.5D0 )
      fw = u0 * u0
      x = 1.0D0 / ( fw * fw ) - b1 / fw - c1
      fw = 0.5D0 * f0 / x
      z(1) = fw * z(1)
      z(2) = fw * z(2)
      z(3) = fw * z(3)
      z(4) = fw * z(4)
      z1(1) = 0.5D0 * ( z(1) + z(3) )
      z1(2) = 0.5D0 * ( z(2) + z(4) )
      z1(3) = 0.5D0 * ( z(2) - z(4) )
      z1(4) = 0.5D0 * ( z(1) - z(3) )
      a2 = ( fk1 * fk1 - a(2,2,1) ) / a(2,1,1)
      b2 = ( fk2 * fk2 - a(2,2,1) ) / a(2,1,1)
      x = b(1,1) * b(4,1) - b(3,1) * b(2,1)
      fw1 = fk1 / x
      fw2 = fk2 / x
      y = fw2 * ( b2 * b(2,1) - b(4,1) )
      zx = fw1 * ( a2 * b(2,1) - b(4,1) )
      a1(1,1) = 0.5D0 * ( 1 - y )
      a1(1,2) = 0.5D0 * ( 1 - zx )
      a1(1,3) = 0.5D0 * ( 1 + zx )
      a1(1,4) = 0.5D0 * ( 1 + y )
      y = fw2 * ( b(3,1) - b2 * b(1,1) )
      zx = fw1 * ( b(3,1) - a2 * b(1,1) )
      a1(2,1) = 0.5D0 * ( b2 - y )
      a1(2,2) = 0.5D0 * ( a2 - zx )
      a1(2,3) = 0.5D0 * ( a2 + zx )
      a1(2,4) = 0.5D0 * ( b2 + y )
      a1(3,1) = a1(2,4)
      a1(3,2) = a1(2,3)
      a1(3,3) = a1(2,2)
      a1(3,4) = a1(2,1)
      a1(4,1) = a1(1,4)
      a1(4,2) = a1(1,3)
      a1(4,3) = a1(1,2)
      a1(4,4) = a1(1,1)
      fq0 = exp ( - t0 / u0 )
      fq1 = exp ( - t1 / u0 )
      x = exp ( - fk1 * dt )
      y = exp ( - fk2 * dt )
      do 40 i = 1, 4
         zz(i,1) = z1(i) * fq0
         zz(i,2) = z1(i) * fq1
         aa(i,1,1) = a1(i,1)
         aa(i,2,1) = a1(i,2)
         aa(i,3,1) = a1(i,3) * x
         aa(i,4,1) = a1(i,4) * y
         aa(i,3,2) = a1(i,3)
         aa(i,4,2) = a1(i,4)
         aa(i,1,2) = a1(i,1) * y
         aa(i,2,2) = a1(i,2) * x
 40   continue
      end subroutine coeffl


c **********************************************************************
c In the limits of no scattering ( Fu, 1991 ), fk1 = 1.0 / u(3) and
c fk2 = 1.0 / u(4).
c **********************************************************************
      subroutine coefft0(t0,   t1,
     O                   z1,  fk1,  fk2,   a1,   zz,   aa)

! 13/08/2016 jjb: final changes before release
!                 removal of 2 unused argument,
!                 removal of unused parameters and common block,
!                 cleaning
! 25/06/2021 jjb: all missing declarations and implicit none

      implicit none

! input
      double precision t0, t1
! output
      double precision z1(4), fk1, fk2, a1(4,4), zz(4,2), aa(4,4,2)
! Local scalars:
      integer i, j, k
      double precision dt, x, y

      fk1 = 4.7320545D0
      fk2 = 1.2679491D0
!     y = exp ( - ( t1 - t0 ) / u0 ) ! jjb internal function not used
!     fw = 0.5 * f0                  ! jjb unreferenced
      do 10 i = 1, 4
         z1(i) = 0.0D0
         zz(i,1) = 0.0D0
         zz(i,2) = 0.0D0
         do 10 j = 1, 4
            a1(i,j) = 0.0D0
            do 10 k = 1, 2
               aa(i,j,k) = 0.0D0
 10   continue
      do 20 i = 1, 4
         j = 5 - i
         a1(i,j) = 1.0D0
 20   continue
      dt = t1 - t0
      x = exp ( - fk1 * dt )
      y = exp ( - fk2 * dt )
      aa(1,4,1) = y
      aa(2,3,1) = x
      aa(3,2,1) = 1.0D0
      aa(4,1,1) = 1.0D0
      aa(1,4,2) = 1.0D0
      aa(2,3,2) = 1.0D0
      aa(3,2,2) = x
      aa(4,1,2) = y
      end subroutine coefft0


c **********************************************************************
      subroutine qcfel(nrfl,  ab,    b,    x)
c **********************************************************************
c 1. `qcfel' is the abbreviation of ` qiu constants for each layer'.
c 2. The inhomogeneous atmosphere is divided into n adjacent homogeneous
c    layers where the  single scattering properties are constant in each
c    layer and allowed to vary from one to another. Delta-four-stream is
c    employed for each homogeneous layer. The boundary conditions at the
c    top and bottom of the atmosphere,  together with  continuity condi-
c    tions  at  layer interfaces lead to a system of algebraic equations
c    from which 4*n unknown constants in the problom can be solved.
c 3. This subroutine is used for solving the 4*n unknowns of A *X = B by
c    considering the fact that the coefficient matrix is a sparse matrix
c    with the precise pattern in this special problom.
c 4. The method is not different in principle from the general scheme of
c    Gaussian elimination with backsubstitution, but carefully optimized
c    so as to minimize arithmetic operations.  Partial  pivoting is used
c    to quarantee  method's numerical stability,  which will  not change
c    the basic pattern of sparsity of the matrix.
c 5. Scaling special problems so as to make  its nonzero matrix elements
c    have comparable magnitudes, which will ameliorate the stability.
c 6. a, b and x present A, B and X in A*X=B, respectively. and n4=4*n.
c 7. AB(13,4*n) is the matrix A in band storage, in rows 3 to 13; rows 1
c    and 2 and other unset elements should be set to zero on entry.
c 8. The jth column of A is stored in the jth column of the array AB  as
c    follows:
c            AB(8+i-j,j) = A(i,j) for max(1,j-5) <= i <= min(4*n,j+5).
c    Reversedly, we have
c            A(ii+jj-8,jj) = AB(ii,jj).
c **********************************************************************

! 13/08/2016 jjb: final changes before release
!        - USE global_params for vertical grid parameters
!        - removal of unused parameters
!        - cleaning
! 25/06/2021 jjb: all missing declarations and implicit none

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      implicit none

      INTEGER MAXLAY,ndfs4
      PARAMETER (MAXLAY=nrlay)
      parameter (ndfs4 = 4 * MAXLAY)

c input/output
      double precision ab(13,ndfs4), b(ndfs4), x(ndfs4)
c local
      integer i, i0, i0f, i0m1, ifq, im1, j, k, k44, l,
     $     m, m1, m18, m1f, m2, m28, m3, m38, m4, m48,
     $     n, n1, n2, n3, n4, n44, nrfl
      double precision p, t, xx, yy

      n = nrfl
      n4 = n*4
      do 5 k = 1, n - 1
         k44 = 4 * k - 4
         do 3 l= 1, 4
            m1 = k44 + l
            p = 0.0D0
            do 10 i = 8, 14 - l
               if ( abs ( ab(i,m1) ) .gt. abs ( p ) ) then
                  p = ab(i,m1)
                  i0 = i
               endif
 10         continue
            i0m1 = i0 + m1
            m18 = m1 + 8
            if ( i0 .eq. 8 ) goto 20
            do 15 j = m1, m1 + 8 - l
               i0f = i0m1 - j
               m1f = m18 - j
               t = ab(i0f,j)
               ab(i0f,j) = ab(m1f,j)
               ab(m1f,j) = t
 15         continue
            i0f = i0m1 - 8
            t = b(i0f)
            b(i0f) = b(m1)
            b(m1) = t
 20         continue
            yy = ab(8,m1)
            ab(8,m1) = 1.0D0
            do 25 j = m1 + 1, m1 + 8 - l
               m1f = m18 - j
               ab(m1f,j) = ab(m1f,j) / yy
 25         continue
            b(m1) = b(m1) / yy
            do 30 i = 9, 14 - l
               xx = ab(i,m1)
               ab(i,m1) = 0.0D0
               im1 = i + m1
               do 35 j = m1 + 1, m1 + 8 - l
                  ifq = im1 - j
                  m1f = m18 - j
                  ab(ifq,j) = ab(ifq,j) - ab(m1f,j) * xx
 35            continue
               ifq = im1 - 8
               b(ifq) = b(ifq) - b(m1) * xx
 30         continue
 3       continue
 5    continue
      n44 = n4 - 4
      do 40 l = 1, 3
         m1 = n44 + l
         p = 0.0D0
         do 45 i = 8, 12 - l
            if ( abs ( ab(i,m1) ) .gt. abs ( p ) ) then
               p = ab(i,m1)
               i0 = i
            endif
 45      continue
         i0m1 = i0 + m1
         m18 = m1 + 8
         if( i0 .eq. 8 ) goto 55
         do 50 j = m1, m1 + 4 - l
            i0f = i0m1 - j
            m1f = m18 - j
            t = ab(i0f,j)
            ab(i0f,j) = ab(m1f,j)
            ab(m1f,j) = t
 50      continue
         i0f = i0m1 - 8
         t = b(i0f)
         b(i0f) = b(m1)
         b(m1) = t
 55      continue
         yy = ab(8,m1)
         ab(8,m1) = 1.0D0
         do 60 j = m1 + 1, m1 + 4 - l
            m1f = m18 - j
            ab(m1f,j) = ab(m1f,j) / yy
 60      continue
         b(m1) = b(m1) / yy
         do 65 i = 9, 12 - l
            xx = ab(i,m1)
            ab(i,m1) = 0.0D0
            im1 = i + m1
            do 70 j = m1 + 1, m1 + 4 - l
               ifq = im1 - j
               m1f = m18 - j
               ab(ifq,j) = ab(ifq,j) - ab(m1f,j) * xx
 70         continue
            ifq = im1 - 8
            b(ifq) = b(ifq) - b(m1) * xx
 65      continue
 40   continue
      yy = ab(8,n4)
      ab(8,n4) = 1.0D0
      b(n4) = b(n4) / yy
      n3 = n4 - 1
      n2 = n3 - 1
      n1 = n2 - 1
      x(n4) = b(n4)
      x(n3) = b(n3) - ab(7,n4) * x(n4)
      x(n2) = b(n2) - ab(7,n3) * x(n3) - ab(6,n4) * x(n4)
      x(n1) = b(n1) - ab(7,n2) * x(n2) - ab(6,n3) * x(n3) -
     1     ab(5,n4) * x(n4)
      do 80 k = 1, n - 1
         m4 = 4 * ( n - k )
         m3 = m4 - 1
         m2 = m3 - 1
         m1 = m2 - 1
         m48 = m4 + 8
         m38 = m3 + 8
         m28 = m2 + 8
         m18 = m1 + 8
         x(m4) = b(m4)
         do 85 m = m4 + 1, m4 + 4
            x(m4) = x(m4) - ab(m48-m,m) * x(m)
 85      continue
         x(m3) = b(m3)
         do 90 m = m3 + 1, m3 + 5
            x(m3) = x(m3) - ab(m38-m,m) * x(m)
 90      continue
         x(m2) = b(m2)
         do 95 m = m2 + 1, m2 + 6
            x(m2) = x(m2) - ab(m28-m,m) * x(m)
 95      continue
         x(m1) = b(m1)
         do 100 m = m1 + 1, m1 + 7
            x(m1) = x(m1) - ab(m18-m,m) * x(m)
 100     continue
 80   continue
      end

*******************************************************************
      SUBROUTINE FLUX_OUT(FILEN, FACT, FLX, NWS)

! 13/08/2016 jjb
!     changes done in this SR:
!        - removal of 3 unused arguments
!        - USE global_params for vertical grid parameters
!        - removal of unused parameters
!        - a few unused variables commented
!        - correction of FLX dimension
!        - little cleaning
!        - implicit none

      USE file_unit, ONLY :
! Imported Parameters:
     &     jpfunf4st

      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      IMPLICIT NONE

      INTEGER MAXLAY,MAXWAV,NW
      PARAMETER (MAXLAY=nrlay, MAXWAV=176, NW=7)

      DOUBLE PRECISION
     $      FACT(0:MAXLAY,NW), !act. flux from DISORT
     $      FLX(NW)               !extraterrestrial flux   ! jjb corrected 2015-11-19 see Forcheck list nb 149 & notebook p. 119

      DOUBLE PRECISION TP(0:MAXLAY)

      COMMON/WL/WAVE(MAXWAV),  !wavelength in the middle of the interval [cm]
     $          DWAVE(MAXWAV)  !width of the wavelength intervals [cm]
      DOUBLE PRECISION WAVE,DWAVE

      CHARACTER (LEN=4) FILEN

      INTEGER
     $     NWS(NW)             !specification of interval: 1 < NWS(L) < MAXWAV

      INTEGER K,L ! indexes for do loops

c-------------------------------------------------------------------


      OPEN(jpfunf4st,FILE=FILEN//'.out',STATUS='UNKNOWN')

      DO L  = 1,NW
         DO K=0,MAXLAY
            IF (FACT(K,L).lt.1.D-25) THEN
               TP(K)=0D0
            ELSE
               TP(K)=FACT(K,L)/(DWAVE(NWS(L))*1.D7) !-> [photons/(cm^2 s nm)]
            ENDIF
         ENDDO
         WRITE(jpfunf4st,'(1p,6E12.4)')TP
      ENDDO
      WRITE(jpfunf4st,'(1p,6E12.4)')(FLX(L)/(DWAVE(NWS(L))*1.D7),L=1,NW)

      CLOSE(jpfunf4st)

      END SUBROUTINE FLUX_OUT
********************************************************************


********************************************************************
      SUBROUTINE LOOKUP
********************************************************************
* initialization of lookup table for photolysis and heating rate   *
* calculation with PHOTO_CAL                                      *
********************************************************************

! 13/08/2016 jjb: final changes before release
!                 removal of unused parameters
!                 use of module "directories"
!                 removal of unused common blck ("haloe")
!                 removal of unused data
!                 initialisation on named COMMON moved in block data
!                 cleaning
! 05/11/2017 jjb: replaced directories by config

      USE config, ONLY : cinpdir_phot     ! input directory for photolysis data files
      USE file_unit, ONLY : nlt=>jpfunnlt ! unit to open and read lookup table ('lookt0900.dat')
      USE precision, ONLY : dp

      IMPLICIT NONE

! Local scalar:
      INTEGER i, j, k

*--------------------------------------------------------------------*
*     LOOK-UP TABLE                                                  *
*--------------------------------------------------------------------*
      COMMON/F_O2 /
     $              CS_O2(58,3,2),     FS_O2(58,2,2),
     $              CS_O3(58,2,2),     FS_O3(58,2,2),
     $              CS_N2O(58,2,2),    FS_N2O(58,2,2),
     $              CS_CFC11(58,3,2),  FS_CFC11(58,2,2),
     $              CS_CFC12(58,3,2),  FS_CFC12(58,2,2),
     $              CS_H2O2(58,2,2),   FS_H2O2(58,2,2),
     $              CS_HNO3(58,2,2),   FS_HNO3(58,2,2),
     $              CS_HNO4(58,2,2),   FS_HNO4(58,2,2),
     $              CS_ClONO2(58,2,2), FS_ClONO2(58,1,2),
     $              CS_BrNO3(58,2,2),
     $              CS_Cl2O2(58,2,2),
     $              CS_HOCl(58,2,2),
     $              CS_N2O5(58,2,2),   FS_N2O5(58,1,2),
     $              CS_HO2(58,3,2),    FS_HO2(58,2,2),
     $              CS_HO3(58,2,2),    FS_HO3(58,2,2),
     $              CS_BrCl_noT(58,3,2),
     $              CS_ClNO2(58,3,2),
     $              CS_BrNO2(58,3,2),
     $              CS_Br2(58,2,2),
     $              CS_CH3I(58,2,2),
     $              CS_NO3n(58,3,2),   FS_NO3n(58,2,2),
     $              CS_dumm23(58,2,2),
     $              CS_dumm24(58,2,2),
     $              CS_dumm25(58,2,2),
     $              CS_dumm26(58,2,2)
      REAL (KIND=DP) ::
     $     CS_O2, FS_O2, CS_O3, FS_O3, CS_N2O, FS_N2O,
     $     CS_CFC11, FS_CFC11, CS_CFC12, FS_CFC12,
     $     CS_H2O2, FS_H2O2, CS_HNO3, FS_HNO3, CS_HNO4, FS_HNO4,
     $     CS_ClONO2, FS_ClONO2, CS_BrNO3, CS_Cl2O2, CS_HOCl,
     $     CS_N2O5, FS_N2O5, CS_HO2, FS_HO2, CS_HO3, FS_HO3,
     $     CS_BrCl_noT, CS_ClNO2, CS_BrNO2, CS_Br2, CS_CH3I,
     $     CS_NO3n, FS_NO3n, CS_dumm23, CS_dumm24, CS_dumm25, CS_dumm26
      
      COMMON/LOOK/
     $          TAUA1(55,3),     TAUB1(55,3),
     $          A1_O3(55,3),     B1_O3(55,3),
     $          A1_O2(55,2),     B1_O2(55,2),
     $          A1_H2O2(55,2),   B1_H2O2(55,2),
     $          A1_HNO3(55,2),   B1_HNO3(55,2),
     $          A1_HNO4(55,2),   B1_HNO4(55,2),
     $          A1_N2O5(55,2),   B1_N2O5(55,2),
     $          A1_CH3OOH(55,2), B1_CH3OOH(55,2),
     $          A1_N2O(55,3),    B1_N2O(55,3),
     $          A1_CFC11(55,2),  B1_CFC11(55,2),
     $          A1_CFC12(55,3),  B1_CFC12(55,3),
     $          A1_ClONO2(55,2), B1_ClONO2(55,2),
     $          A1_BrNO3(55,2),  B1_BrNO3(55,2),
     $          A1_Cl2O2(55,3),  B1_Cl2O2(55,3),
     $          A1_HOCl(55,2),   B1_HOCl(55,2),
     $          A1_H_O2(55,2),   B1_H_O2(55,2),
     $          A1_H_O3(55,2),   B1_H_O3(55,2),
     $          A1_BrCl_noT(55,2),  B1_BrCl_noT(55,2),
     $          A1_ClNO2(55,2),  B1_ClNO2(55,2),
     $          A1_BrNO2(55,2),  B1_BrNO2(55,2),
     $          A1_Br2(55,2),    B1_Br2(55,2),
     $          A1_CH3I(55,2),   B1_CH3I(55,2),
     $          A1_ICl(55,2),    B1_ICl(55,2),
     $          A1_IBr(55,2),    B1_IBr(55,2),
     $          A1_C3H7I(55,2),  B1_C3H7I(55,2),
     $          A1_CH2ClI(55,2), B1_CH2ClI(55,2),
     $          A1_CH2I2(55,2),  B1_CH2I2(55,2),
     $          A1_INO2(55,2),   B1_INO2(55,2),
     $          A1_Cl2_noT(55,2),  B1_Cl2_noT(55,2),
     $          A1_NO3n(55,3),   B1_NO3n(55,3),
     $          A1_dumm23(55,2),  B1_dumm23(55,2),
     $          A1_dumm24(55,2),  B1_dumm24(55,2),
     $          A1_dumm25(55,2),  B1_dumm25(55,2),
     $          A1_dumm26(55,2),  B1_dumm26(55,2),
     $          TAUA2(55),       TAUB2(55),
     $          A2_TO3(55),      B2_TO3(55),
     $          A2_O1D(55),      B2_O1D(55),
     $          A2_O3P(55),      B2_O3P(55),
     $          A2_H2O2(55),     B2_H2O2(55),
     $          A2_HNO3(55),     B2_HNO3(55),
     $          A2_HNO4(55),     B2_HNO4(55),
     $          A2_N2O5(55),     B2_N2O5(55),
     $          A2_CH3OOH(55),   B2_CH3OOH(55),
     $          A2_NO2(55),      B2_NO2(55),
     $          A2_ClONO2(55),   B2_ClONO2(55),
     $          A2_BrNO3(55),    B2_BrNO3(55),
     $          A2_Cl2O2(55),    B2_Cl2O2(55),
     $          A2_HOCl(55),     B2_HOCl(55),
     $          A2_H_O3(55),     B2_H_O3(55),
     $          A2_H_NO2(55),    B2_H_NO2(55), ! jjb H_NO2 is no longer calculated, but these values are still read
     $          A2_BrCl_noT(55),    B2_BrCl_noT(55),
     $          A2_ClNO2(55),    B2_ClNO2(55),
     $          A2_BrNO2(55),    B2_BrNO2(55),
     $          A2_Br2(55),      B2_Br2(55),
     $          A2_INO3(55),     B2_INO3(55),
     $          A2_CH3I(55),     B2_CH3I(55),
     $          A2_ICl(55),      B2_ICl(55),
     $          A2_IBr(55),      B2_IBr(55),
     $          A2_C3H7I(55),    B2_C3H7I(55),
     $          A2_CH2ClI(55),   B2_CH2ClI(55),
     $          A2_CH2I2(55),    B2_CH2I2(55),
     $          A2_INO2(55),     B2_INO2(55),
     $          A2_OClO_noT(55), B2_OClO_noT(55),
     $          A2_Cl2_noT(55),  B2_Cl2_noT(55),
     $          A2_HOBr(55),     B2_HOBr(55),
     $          A2_NO3n(55),     B2_NO3n(55),
     $          A2_dumm23(55),    B2_dumm23(55),
     $          A2_dumm24(55),    B2_dumm24(55),
     $          A2_dumm25(55),    B2_dumm25(55),
     $          A2_dumm26(55),    B2_dumm26(55),
     $          TAUA3(55),       TAUB3(55),
     $          A3_TO3(55),      B3_TO3(55),
     $          A3_O1D(55),      B3_O1D(55),
     $          A3_O3P(55),      B3_O3P(55),
     $          A3_H2O2(55),     B3_H2O2(55),
     $          A3_HNO3(55),     B3_HNO3(55),
     $          A3_HNO4(55),     B3_HNO4(55),
     $          A3_N2O5(55),     B3_N2O5(55),
     $          A3_CH3OOH(55),   B3_CH3OOH(55),
     $          A3_NO2(55),      B3_NO2(55),
     $          A3_COH2(55),     B3_COH2(55),
     $          A3_CHOH(55),     B3_CHOH(55),
     $          A3_ClONO2(55),   B3_ClONO2(55),
     $          A3_BrNO3(55),    B3_BrNO3(55),
     $          A3_Cl2O2(55),    B3_Cl2O2(55),
     $          A3_HOCl(55),     B3_HOCl(55),
     $          A3_H_O3(55),     B3_H_O3(55),
     &          A3_H_NO2(55),    B3_H_NO2(55), ! jjb H_NO2 is no longer calculated, but these values are still read
     $          A3_BrCl_noT(55), B3_BrCl_noT(55),
     $          A3_ClNO2(55),    B3_ClNO2(55),
     $          A3_BrNO2(55),    B3_BrNO2(55),
     $          A3_Br2(55),      B3_Br2(55),
     $          A3_INO3(55),     B3_INO3(55),
     $          A3_CH3I(55),     B3_CH3I(55),
     $          A3_ICl(55),      B3_ICl(55),
     $          A3_IBr(55),      B3_IBr(55),
     $          A3_C3H7I(55),    B3_C3H7I(55),
     $          A3_CH2ClI(55),   B3_CH2ClI(55),
     $          A3_CH2I2(55),    B3_CH2I2(55),
     $          A3_INO2(55),     B3_INO2(55),
     $          A3_BrO_noT(55),  B3_BrO_noT(55),
     $          A3_OClO_noT(55), B3_OClO_noT(55),
     $          A3_Cl2_noT(55),  B3_Cl2_noT(55),
     $          A3_HOI_jen91(55),B3_HOI_jen91(55),
     $          A3_HOBr(55),     B3_HOBr(55),
     $          A3_NO2m(55),     B3_NO2m(55),
     $          A3_NO3n(55),     B3_NO3n(55),
     $          A3_dumm23(55),   B3_dumm23(55),
     $          A3_dumm24(55),   B3_dumm24(55),
     $          A3_dumm25(55),   B3_dumm25(55),
     $          A3_dumm26(55),   B3_dumm26(55)
      REAL (KIND=DP) ::
     $     TAUA1, TAUB1, A1_O3, B1_O3, A1_O2, B1_O2, A1_H2O2, B1_H2O2,
     $     A1_HNO3, B1_HNO3, A1_HNO4, B1_HNO4, A1_N2O5, B1_N2O5,
     $     A1_CH3OOH, B1_CH3OOH, A1_N2O, B1_N2O, A1_CFC11, B1_CFC11,
     $     A1_CFC12, B1_CFC12, A1_ClONO2, B1_ClONO2, A1_BrNO3, B1_BrNO3,
     $     A1_Cl2O2, B1_Cl2O2, A1_HOCl, B1_HOCl, A1_H_O2, B1_H_O2,
     $     A1_H_O3, B1_H_O3, A1_BrCl_noT, B1_BrCl_noT,
     $     A1_ClNO2, B1_ClNO2, A1_BrNO2, B1_BrNO2, A1_Br2, B1_Br2,
     $     A1_CH3I, B1_CH3I, A1_ICl, B1_ICl, A1_IBr, B1_IBr,
     $     A1_C3H7I, B1_C3H7I, A1_CH2ClI, B1_CH2ClI, A1_CH2I2, B1_CH2I2,
     $     A1_INO2, B1_INO2, A1_Cl2_noT, B1_Cl2_noT, A1_NO3n, B1_NO3n,
     $     A1_dumm23, B1_dumm23, A1_dumm24, B1_dumm24,
     $     A1_dumm25, B1_dumm25, A1_dumm26, B1_dumm26,
     $     TAUA2, TAUB2,
     $     A2_TO3, B2_TO3, A2_O1D, B2_O1D, A2_O3P, B2_O3P,
     $     A2_H2O2, B2_H2O2, A2_HNO3, B2_HNO3, A2_HNO4, B2_HNO4,
     $     A2_N2O5, B2_N2O5, A2_CH3OOH, B2_CH3OOH, A2_NO2, B2_NO2,
     $     A2_ClONO2, B2_ClONO2, A2_BrNO3, B2_BrNO3, A2_Cl2O2, B2_Cl2O2,
     $     A2_HOCl, B2_HOCl, A2_H_O3, B2_H_O3, A2_H_NO2, B2_H_NO2,
     $     A2_BrCl_noT, B2_BrCl_noT, A2_ClNO2, B2_ClNO2,
     $     A2_BrNO2, B2_BrNO2, A2_Br2, B2_Br2, A2_INO3, B2_INO3,
     $     A2_CH3I, B2_CH3I, A2_ICl, B2_ICl, A2_IBr, B2_IBr,
     $     A2_C3H7I, B2_C3H7I, A2_CH2ClI, B2_CH2ClI, A2_CH2I2, B2_CH2I2,
     $     A2_INO2, B2_INO2, A2_OClO_noT, B2_OClO_noT,
     $     A2_Cl2_noT, B2_Cl2_noT, A2_HOBr, B2_HOBr, A2_NO3n, B2_NO3n,
     $     A2_dumm23, B2_dumm23, A2_dumm24, B2_dumm24,
     $     A2_dumm25, B2_dumm25, A2_dumm26, B2_dumm26,
     $     TAUA3, TAUB3,
     $     A3_TO3, B3_TO3, A3_O1D, B3_O1D, A3_O3P, B3_O3P,
     $     A3_H2O2, B3_H2O2, A3_HNO3, B3_HNO3, A3_HNO4, B3_HNO4,
     $     A3_N2O5, B3_N2O5, A3_CH3OOH, B3_CH3OOH, A3_NO2, B3_NO2,
     $     A3_COH2, B3_COH2, A3_CHOH, B3_CHOH, A3_ClONO2, B3_ClONO2,
     $     A3_BrNO3, B3_BrNO3, A3_Cl2O2, B3_Cl2O2, A3_HOCl, B3_HOCl,
     $     A3_H_O3, B3_H_O3, A3_H_NO2, B3_H_NO2,
     $     A3_BrCl_noT, B3_BrCl_noT, A3_ClNO2, B3_ClNO2,
     $     A3_BrNO2, B3_BrNO2, A3_Br2, B3_Br2, A3_INO3, B3_INO3,
     $     A3_CH3I, B3_CH3I, A3_ICl, B3_ICl, A3_IBr, B3_IBr,
     $     A3_C3H7I, B3_C3H7I, A3_CH2ClI, B3_CH2ClI, A3_CH2I2, B3_CH2I2,
     $     A3_INO2, B3_INO2, A3_BrO_noT, B3_BrO_noT,
     $     A3_OClO_noT, B3_OClO_noT, A3_Cl2_noT, B3_Cl2_noT,
     $     A3_HOI_jen91, B3_HOI_jen91,
     $     A3_HOBr, B3_HOBr, A3_NO2m, B3_NO2m, A3_NO3n, B3_NO3n,
     $     A3_dumm23, B3_dumm23, A3_dumm24, B3_dumm24,
     $     A3_dumm25, B3_dumm25, A3_dumm26, B3_dumm26

      COMMON/C_POLY/
     $   C4_TAU(4),     C5_TAU(3),                    C7_TAU(2),
     $   C4_T_O3(4),    C5_T_O3(3),                   C7_T_O3(2),
     $   C4_O1D(4),     C5_O1D(4),
     $   C4_O3P(4),     C5_O3P(3),     C6_O3P(2),     C7_O3P(2),
     $   C4_H2O2(3),    C5_H2O2(3),    C6_H2O2(2),
     $   C4_HNO3(3),    C5_HNO3(3),    C6_HNO3(2),
     $   C4_HNO4(3),    C5_HNO4(3),
     $   C4_N2O5(3),    C5_N2O5(3),    C6_N2O5(2),
     $   C4_CH3OOH(3),  C5_CH3OOH(3),  C6_CH3OOH(2),
     $   C4_NO2(3),     C5_NO2(3),     C6_NO2(2),
     $   C4_COH2(4),    C5_COH2(3),    C6_COH2(2),
     $   C4_CHOH(4),    C5_CHOH(3),    C6_CHOH(2),
     $                                 C6_NO2O(2),    C7_NO2O(2),
     $                                                C7_NOO2(2),
     $   C4_ClONO2(3),  C5_ClONO2(3),  C6_ClONO2(1),
     $   C4_BrNO3(2),   C5_BrNO3(2),   C6_BrNO3(1),   C7_BrNO3(2),
     $   C4_Cl2O2(2),   C5_Cl2O2(2),   C6_Cl2O2(1),
     $   C4_HOCl(2),    C5_HOCl(2),    C6_HOCl(1),
     $   C4_H_O3(4),    C5_H_O3(3),    C6_H_O3(2),    C7_H_O3(2),
     $                                                C7_H_O2(2),
     $   C4_H_NO2(3),   C5_H_NO2(2),   C6_H_NO2(2),   C7_H_NO2(2), ! jjb H_NO2 is no longer calculated, but these values are still read
     $                                 C6_H_O4(2),    C7_H_O4(2),
     $   C4_BrCl_noT(3),C5_BrCl_noT(3),C6_BrCl_noT(1),C7_BrCl_noT(2),
     $   C4_ClNO2(2),   C5_ClNO2(2),   C6_ClNO2(1),
     $   C4_BrNO2(3),   C5_BrNO2(2),   C6_BrNO2(1),   C7_BrNO2(2),
     $   C4_Br2(3),     C5_Br2(3),     C6_Br2(1),     C7_Br2(2),
     $                                 C6_IO(1),      C7_IO(2),
     $   C4_INO3(3),    C5_INO3(2),    C6_INO3(1),    C7_INO3(2),
     $   C4_CH3I(3),    C5_CH3I(3),    C6_CH3I(1),
     $                                 C6_I2(1),      C7_I2(2),
     $   C4_ICl(2),     C5_ICl(2),     C6_ICl(1),     C7_ICl(2),
     $   C4_IBr(3),     C5_IBr(3),     C6_IBr(1),     C7_IBr(2),
     $   C4_C3H7I(3),   C5_C3H7I(3),   C6_C3H7I(1),
     $   C4_CH2ClI(3),  C5_CH2ClI(3),  C6_CH2ClI(1),
     $   C4_CH2I2(2),   C5_CH2I2(2),   C6_CH2I2(1),
     $   C4_INO2(3),    C5_INO2(2),    C6_INO2(1),
     $   C4_BrO_noT(3),  C5_BrO_noT(2),  C6_BrO_noT(1),
     $   C4_OClO_noT(3), C5_OClO_noT(3), C6_OClO_noT(1), C7_OClO_noT(2),
     $   C4_Cl2_noT(3), C5_Cl2_noT(2), C6_Cl2_noT(1), C7_Cl2_noT(2),
     $   C4_HOI_jen91(3),C5_HOI_jen91(3),C6_HOI_jen91(1),
     $   C7_HOI_jen91(2),
     $   C4_HOBr(3),   C5_HOBr(2),   C6_HOBr(1),   C7_HOBr(2),
     $   C4_HONO(3),   C5_HONO(3),   C6_HONO(1),
     $   C4_NO2m(2),   C5_NO2m(2),   C6_NO2m(1),
     $   C4_NO3n(4),   C5_NO3n(4),   C6_NO3n(2),
     $   C4_dumm23(2),   C5_dumm23(2),   C6_dumm23(1),   C7_dumm23(2),
     $   C4_dumm24(2),   C5_dumm24(2),   C6_dumm24(1),   C7_dumm24(2),
     $   C4_dumm25(2),   C5_dumm25(2),   C6_dumm25(1),   C7_dumm25(2),
     $   C4_dumm26(2),   C5_dumm26(2),   C6_dumm26(1),   C7_dumm26(2)
      REAL (KIND=DP)
     $   C4_TAU,     C5_TAU,                 C7_TAU,
     $   C4_T_O3,    C5_T_O3,                C7_T_O3,
     $   C4_O1D,     C5_O1D,
     $   C4_O3P,     C5_O3P,     C6_O3P,     C7_O3P,
     $   C4_H2O2,    C5_H2O2,    C6_H2O2,
     $   C4_HNO3,    C5_HNO3,    C6_HNO3,
     $   C4_HNO4,    C5_HNO4,
     $   C4_N2O5,    C5_N2O5,    C6_N2O5,
     $   C4_CH3OOH,  C5_CH3OOH,  C6_CH3OOH,
     $   C4_NO2,     C5_NO2,     C6_NO2,
     $   C4_COH2,    C5_COH2,    C6_COH2,
     $   C4_CHOH,    C5_CHOH,    C6_CHOH,
     $                           C6_NO2O,    C7_NO2O,
     $                                       C7_NOO2,
     $   C4_ClONO2,  C5_ClONO2,  C6_ClONO2,
     $   C4_BrNO3,   C5_BrNO3,   C6_BrNO3,   C7_BrNO3,
     $   C4_Cl2O2,   C5_Cl2O2,   C6_Cl2O2,
     $   C4_HOCl,    C5_HOCl,    C6_HOCl,
     $   C4_H_O3,    C5_H_O3,    C6_H_O3,    C7_H_O3,
     $                                       C7_H_O2,
     $   C4_H_NO2,   C5_H_NO2,   C6_H_NO2,   C7_H_NO2,
     $                           C6_H_O4,    C7_H_O4,
     $   C4_BrCl_noT,C5_BrCl_noT,C6_BrCl_noT,C7_BrCl_noT,
     $   C4_ClNO2,   C5_ClNO2,   C6_ClNO2,
     $   C4_BrNO2,   C5_BrNO2,   C6_BrNO2,   C7_BrNO2,
     $   C4_Br2,     C5_Br2,     C6_Br2,     C7_Br2,
     $                           C6_IO,      C7_IO,
     $   C4_INO3,    C5_INO3,    C6_INO3,    C7_INO3,
     $   C4_CH3I,    C5_CH3I,    C6_CH3I,
     $                           C6_I2,      C7_I2,
     $   C4_ICl,     C5_ICl,     C6_ICl,     C7_ICl,
     $   C4_IBr,     C5_IBr,     C6_IBr,     C7_IBr,
     $   C4_C3H7I,   C5_C3H7I,   C6_C3H7I,
     $   C4_CH2ClI,  C5_CH2ClI,  C6_CH2ClI,
     $   C4_CH2I2,   C5_CH2I2,   C6_CH2I2,
     $   C4_INO2,    C5_INO2,    C6_INO2,
     $   C4_BrO_noT, C5_BrO_noT, C6_BrO_noT,
     $   C4_OClO_noT,C5_OClO_noT,C6_OClO_noT,C7_OClO_noT,
     $   C4_Cl2_noT, C5_Cl2_noT, C6_Cl2_noT, C7_Cl2_noT,
     $   C4_HOI_jen91,C5_HOI_jen91,C6_HOI_jen91,C7_HOI_jen91,
     $   C4_HOBr,    C5_HOBr,    C6_HOBr,   C7_HOBr,
     $   C4_HONO,    C5_HONO,    C6_HONO,
     $   C4_NO2m,    C5_NO2m,    C6_NO2m,
     $   C4_NO3n,    C5_NO3n,    C6_NO3n,
     $   C4_dumm23,  C5_dumm23,  C6_dumm23,  C7_dumm23,
     $   C4_dumm24,  C5_dumm24,  C6_dumm24,  C7_dumm24,
     $   C4_dumm25,  C5_dumm25,  C6_dumm25,  C7_dumm25,
     $   C4_dumm26,  C5_dumm26,  C6_dumm26,  C7_dumm26

c-------------------------------------------------------------------


      OPEN(NLT, FILE=trim(cinpdir_phot)//'lookt0900.dat', STATUS='OLD')

      READ(NLT,*)
      DO J=1,3
      DO I=1,2
        READ(NLT,'(6E12.4)')(CS_O2(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_N2O(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,3
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_CFC11(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,3
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_CFC12(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_H2O2(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_HNO3(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_HNO4(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_O3(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_ClONO2(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_BrNO3(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_HOCl(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_N2O5(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,3
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_HO2(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_HO3(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,3
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_BrCl_noT(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,3
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_ClNO2(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,3
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_BrNO2(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_Br2(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_CH3I(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,3
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_NO3n(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_dumm23(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_dumm24(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_dumm25(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(CS_dumm26(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_O2(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_N2O(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_CFC11(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_CFC12(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_H2O2(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_HNO3(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_HNO4(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_O3(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,1
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_ClONO2(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,1
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_N2O5(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_NO3n(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_HO2(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      DO J=1,2
      DO I=1,2
         READ(NLT,'(6E12.4)')(FS_HO3(K,J,I),K=1,58)
      ENDDO
      ENDDO
      READ(NLT,*)
      READ(NLT,*)
      DO I=1,3
         READ(NLT,'(6E12.4)')(TAUA1(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(TAUB1(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,3
         READ(NLT,'(6E12.4)')(A1_O3(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_O3(K,I),K=1,55)
      ENDDO

      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_O2(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_O2(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_H2O2(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_H2O2(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_HNO3(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_HNO3(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_HNO4(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_HNO4(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_N2O5(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_N2O5(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_CH3OOH(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_CH3OOH(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,3
         READ(NLT,'(6E12.4)')(A1_N2O(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_N2O(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_CFC11(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_CFC11(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,3
         READ(NLT,'(6E12.4)')(A1_CFC12(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_CFC12(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_ClONO2(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_ClONO2(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_BrNO3(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_BrNO3(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,3
         READ(NLT,'(6E12.4)')(A1_Cl2O2(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_Cl2O2(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_HOCl(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_HOCl(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_H_O2(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_H_O2(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_H_O3(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_H_O3(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_BrCl_noT(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_BrCl_noT(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_ClNO2(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_ClNO2(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_BrNO2(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_BrNO2(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_Br2(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_Br2(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_CH3I(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_CH3I(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_ICl(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_ICl(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_IBr(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_IBr(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_C3H7I(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_C3H7I(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_CH2ClI(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_CH2ClI(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_CH2I2(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_CH2I2(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_INO2(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_INO2(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_Cl2_noT(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_Cl2_noT(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,3
         READ(NLT,'(6E12.4)')(A1_NO3n(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_NO3n(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_dumm23(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_dumm23(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_dumm24(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_dumm24(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_dumm25(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_dumm25(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      DO I=1,2
         READ(NLT,'(6E12.4)')(A1_dumm26(K,I),K=1,55)
         READ(NLT,'(6E12.4)')(B1_dumm26(K,I),K=1,55)
      ENDDO
      READ(NLT,*)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(TAUA2(K),K=1,55)
      READ(NLT,'(6E12.4)')(TAUB2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_TO3(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_TO3(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_O1D(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_O1D(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_O3P(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_O3P(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_H2O2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_H2O2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_HNO3(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_HNO3(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_HNO4(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_HNO4(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_N2O5(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_N2O5(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_CH3OOH(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_CH3OOH(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_NO2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_NO2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_ClONO2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_ClONO2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_BrNO3(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_BrNO3(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_Cl2O2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_Cl2O2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_HOCl(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_HOCl(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_H_O3(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_H_O3(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_H_NO2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_H_NO2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_BrCl_noT(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_BrCl_noT(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_ClNO2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_ClNO2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_BrNO2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_BrNO2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_Br2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_Br2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_INO3(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_INO3(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_CH3I(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_CH3I(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_ICl(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_ICl(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_IBr(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_IBr(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_C3H7I(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_C3H7I(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_CH2ClI(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_CH2ClI(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_CH2I2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_CH2I2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_INO2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_INO2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_OClO_noT(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_OClO_noT(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_Cl2_noT(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_Cl2_noT(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_HOBr(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_HOBr(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_NO3n(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_NO3n(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_dumm23(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_dumm23(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_dumm24(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_dumm24(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_dumm25(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_dumm25(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A2_dumm26(K),K=1,55)
      READ(NLT,'(6E12.4)')(B2_dumm26(K),K=1,55)
      READ(NLT,*)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(TAUA3(K),K=1,55)
      READ(NLT,'(6E12.4)')(TAUB3(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_TO3(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_TO3(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_O1D(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_O1D(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_O3P(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_O3P(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_H2O2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_H2O2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_HNO3(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_HNO3(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_HNO4(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_HNO4(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_N2O5(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_N2O5(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_CH3OOH(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_CH3OOH(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_NO2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_NO2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_COH2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_COH2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_CHOH(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_CHOH(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_ClONO2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_ClONO2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_BrNO3(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_BrNO3(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_Cl2O2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_Cl2O2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_HOCl(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_HOCl(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_H_O3(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_H_O3(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_H_NO2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_H_NO2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_BrCl_noT(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_BrCl_noT(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_ClNO2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_ClNO2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_BrNO2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_BrNO2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_Br2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_Br2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_INO3(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_INO3(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_CH3I(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_CH3I(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_ICl(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_ICl(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_IBr(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_IBr(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_C3H7I(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_C3H7I(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_CH2ClI(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_CH2ClI(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_CH2I2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_CH2I2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_INO2(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_INO2(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_BrO_noT(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_BrO_noT(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_OClO_noT(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_OClO_noT(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_Cl2_noT(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_Cl2_noT(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_HOI_jen91(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_HOI_jen91(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_HOBr(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_HOBr(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_NO2m(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_NO2m(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_NO3n(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_NO3n(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_dumm23(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_dumm23(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_dumm24(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_dumm24(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_dumm25(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_dumm25(K),K=1,55)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')(A3_dumm26(K),K=1,55)
      READ(NLT,'(6E12.4)')(B3_dumm26(K),K=1,55)
      READ(NLT,*)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_TAU
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_T_O3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_O1D
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_O3P
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_H2O2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_HNO3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_HNO4
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_N2O5
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_CH3OOH
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_NO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_COH2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_CHOH
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_ClONO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_BrNO3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_Cl2O2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_HOCl
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_H_O3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_H_NO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_BrCl_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_ClNO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_BrNO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_Br2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_INO3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_CH3I
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_ICl
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_IBr
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_C3H7I
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_CH2ClI
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_CH2I2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_INO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_BrO_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_OClO_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_Cl2_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_HOI_jen91
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_HOBr
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_HONO
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_NO2m
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_NO3n
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_dumm23
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_dumm24
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_dumm25
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C4_dumm26
      READ(NLT,*)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_TAU
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_T_O3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_O1D
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_O3P
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_H2O2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_HNO3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_HNO4
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_N2O5
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_CH3OOH
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_NO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_COH2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_CHOH
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_ClONO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_BrNO3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_Cl2O2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_HOCl
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_H_O3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_H_NO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_BrCl_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_ClNO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_BrNO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_Br2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_INO3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_CH3I
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_ICl
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_IBr
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_C3H7I
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_CH2ClI
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_CH2I2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_INO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_BrO_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_OClO_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_Cl2_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_HOI_jen91
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_HOBr
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_HONO
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_NO2m
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_NO3n
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_dumm23
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_dumm24
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_dumm25
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C5_dumm26
      READ(NLT,*)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_O3P
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_H2O2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_HNO3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_N2O5
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_CH3OOH
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_NO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_COH2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_CHOH
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_NO2O
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_ClONO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_BrNO3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_Cl2O2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_HOCl
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_H_O3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_H_NO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_H_O4
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_BrCl_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_ClNO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_BrNO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_Br2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_IO
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_INO3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_CH3I
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_I2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_ICl
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_IBr
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_C3H7I
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_CH2ClI
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_CH2I2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_INO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_BrO_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_OClO_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_Cl2_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_HOI_jen91
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_HOBr
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_HONO
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_NO2m
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_NO3n
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_dumm23
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_dumm24
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_dumm25
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C6_dumm26
      READ(NLT,*)
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_TAU
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_T_O3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_O3P
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_BrNO3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_NO2O
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_NOO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_H_O3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_H_O2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_H_NO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_H_O4
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_BrCl_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_BrNO2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_Br2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_IO
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_INO3
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_I2
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_ICl
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_IBr
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_OClO_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_Cl2_noT
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_HOI_jen91
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_HOBr
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_dumm23
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_dumm24
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_dumm25
      READ(NLT,*)
      READ(NLT,'(6E12.4)')C7_dumm26
      close(NLT)

      END SUBROUTINE LOOKUP
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      block data lookup_data

! 13/08/2016 jjb final changes before release
!    moved data initialisation for common blocks from the SR above (LOOKUP)
!    this is the proper way of initialisation data for common blocks, and has already lead to bugs
!    in other parts of the program when done differently

      IMPLICIT NONE

      COMMON/FIL/ IFIL(115)
      INTEGER IFIL

      DATA IFIL /
     $      1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     $     11,  12,  13,  14,  15,  16,  17,  18,  19,  20,
     $     21,  22,  23,  24,  25,  26,  27,  28,  29,  30,
     $     31,  32,  33,  34,  35,  36,  37,  38,  39,  40,
     $     41,  41,  41,  41,  41,  42,  42,  42,  42,  42,
     $     43,  43,  43,  43,  43,  44,  44,  44,  44,  44,
     $     45,  45,  45,  45,  45,  46,  46,  46,  46,  46,
     $     47,  47,  47,  47,  47,  48,  48,  48,  48,  48,
     $     49,  49,  49,  49,  49,  50,  50,  50,  50,  50,
     $     51,  51,  51,  51,  51,  52,  52,  52,  52,  52,
     $     53,  53,  53,  53,  53,  54,  54,  54,  54,  54,
     $     55,  55,  55,  55,  55/

      COMMON /T_COEFF/   TJ_O1D(4,7),   TJ_O3P(3,7),  TJ_NO2(2,7),
     $    TJ_HNO3(4,7),  TJ_N2O5(4,7),  TJ_H2O2(3,7), TJ_NOO2(2,7),
     $    TJ_NO2O(2,7),  TJ_N2O(3,7),   TJ_CFC12(4,7),TJ_CFC11(3,7),
     $    TJ_ClONO2(3,7),TJ_CHOH(2,7),  TJ_COH2(2,7), TJ_NO3n(3,6),
     $    TH_O3(2,7),    TH_O4(2,7),   TH_NO2(2,7)
      DOUBLE PRECISION TJ_O1D, TJ_O3P, TJ_NO2, TJ_HNO3, TJ_N2O5,TJ_H2O2,
     $     TJ_NOO2, TJ_NO2O, TJ_N2O,  TJ_CFC12, TJ_CFC11, TJ_ClONO2,
     $     TJ_CHOH, TJ_COH2, TJ_NO3n, TH_O3,    TH_O4,    TH_NO2

      DATA TJ_O1D /1.00D0,  0.00D0,  0.00D0,  0.00D0,
     $             1.00D0,  0.08D0,  0.00D0,  0.00D0,
     $             1.00D0, -0.03D0,  0.00D0,  0.00D0,
     $             1.00D0,  0.88D0,  0.10D0, -0.58D0,
     $             1.00D0,  5.14D0,  7.92D0,  0.91D0,
     $             0.00D0,  0.00D0,  0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,  0.00D0,  0.00D0/

      DATA TJ_HNO3/1.00D0,  0.42D0,  0.07D0,  0.00D0,
     $             1.00D0,  0.51D0,  0.16D0,  0.00D0,
     $             1.00D0,  0.80D0,  0.30D0,  0.00D0,
     $             1.00D0,  1.05D0,  0.56D0,  0.20D0,
     $             1.00D0,  1.62D0,  1.40D0,  0.82D0,
     $             1.00D0,  2.72D0,  3.88D0,  3.50D0,
     $             0.00D0,  0.00D0,  0.00D0,  0.00D0/

      DATA TJ_N2O5/1.00D0,  0.00D0,  0.00D0,  0.00D0,
     $             1.00D0,  0.47D0, -0.21D0,  0.00D0,
     $             1.00D0,  1.79D0, -0.14D0, -0.49D0,
     $             1.00D0,  2.31D0,  0.42D0, -0.99D0,
     $             1.00D0,  3.29D0,  2.22D0, -1.37D0,
     $             1.00D0,  5.13D0,  8.01D0,  2.30D0,
     $             0.00D0,  0.00D0,  0.00D0,  0.00D0/


      DATA TJ_CFC12/1.0D0,  2.10D0, 2.31D0,  1.68D0,
     $             0.00D0,  0.00D0, 0.00D0,  0.00D0,
     $             0.00D0,  0.00D0, 0.00D0,  0.00D0,
     $             0.00D0,  0.00D0, 0.00D0,  0.00D0,
     $             0.00D0,  0.00D0, 0.00D0,  0.00D0,
     $             0.00D0,  0.00D0, 0.00D0,  0.00D0,
     $             0.00D0,  0.00D0, 0.00D0,  0.00D0/

      DATA TJ_CFC11/1.0D0,  0.52D0, 0.13D0,
     $             0.00D0,  0.00D0, 0.00D0,
     $             0.00D0,  0.00D0, 0.00D0,
     $             0.00D0,  0.00D0, 0.00D0,
     $             0.00D0,  0.00D0, 0.00D0,
     $             0.00D0,  0.00D0, 0.00D0,
     $             0.00D0,  0.00D0, 0.00D0/

      DATA TJ_ClONO2/1.0D0, 0.04D0,-0.53D0,
     $             1.00D0,  0.67D0, 0.23D0,
     $             1.00D0,  0.89D0, 0.63D0,
     $             1.00D0,  1.06D0, 0.85D0,
     $             1.00D0,  1.12D0, 1.73D0,
     $             1.00D0,  0.68D0, 0.32D0,
     $             0.00D0,  0.00D0, 0.00D0/

      DATA TJ_O3P /1.00D0,  0.00D0,  0.00D0,
     $             1.00D0,  0.07D0,  0.00D0,
     $             1.00D0, -0.40D0,  0.11D0,
     $             1.00D0, -0.93D0, -0.11D0,
     $             1.00D0, -0.36D0, -0.55D0,
     $             1.00D0,  0.00D0,  0.00D0,
     $             1.00D0,  0.00D0,  0.00D0/

      DATA TJ_H2O2/1.00D0,  0.00D0,  0.00D0,
     $             1.00D0,  0.21D0,  0.36D0,
     $             1.00D0,  0.32D0,  0.53D0,
     $             1.00D0,  0.42D0,  0.68D0,
     $             1.00D0,  0.63D0,  1.03D0,
     $             1.00D0,  0.78D0,  1.28D0,
     $             0.00D0,  0.00D0,  0.00D0/

      DATA TJ_N2O/ 1.00D0,  1.33D0,  0.91D0,
     $             0.00D0,  0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,  0.00D0/

      DATA TJ_CHOH/1.00D0,  0.00D0,
     $             1.00D0,  0.00D0,
     $             1.00D0, -0.17D0,
     $             1.00D0,  0.00D0,
     $             1.00D0, -0.25D0,
     $             1.00D0, -0.37D0,
     $             1.00D0,  0.00D0/

      DATA TJ_COH2/1.00D0,  0.00D0,
     $             1.00D0,  0.00D0,
     $             1.00D0, -0.17D0,
     $             1.00D0,  0.00D0,
     $             1.00D0, -0.31D0,
     $             1.00D0, -0.32D0,
     $             1.00D0,  0.00D0/

      DATA TJ_NO2 /1.00D0,  0.00D0,
     $             1.00D0,  0.00D0,
     $             1.00D0,  0.00D0,
     $             1.00D0,  0.00D0,
     $             1.00D0,  0.00D0,
     $             1.00D0,  0.00D0,
     $             1.00D0,  0.00D0/

      DATA TJ_NO2O/0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,
     $             1.00D0,  0.00D0,
     $             1.00D0, -0.22D0/

      DATA TJ_NOO2/0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,
     $             0.00D0,  0.00D0,
     $             1.00D0, -0.32D0/

      DATA TJ_NO3n/1.00D0, -0.01D0, -0.01D0,
     $             1.00D0,  0.00D0, -0.03D0,
     $             1.00D0,  0.00D0,  0.00D0,
     $             1.00D0,  0.00D0, -0.00D0,
     $             1.00D0,  0.00D0, -0.00D0,
     $             1.00D0,  0.00D0, -0.00D0/


      DATA TH_O3/  1.00D0,  0.00D0,
     $             1.00D0,  0.08D0,
     $             1.00D0, -0.04D0,
     $             1.00D0,  0.00D0,
     $             1.00D0,  0.00D0,
     $             1.00D0,  0.00D0,
     $             1.00D0,  0.00D0/

      DATA TH_NO2/ 1.00D0,  0.00D0,
     $             1.00D0, -0.01D0,
     $             1.00D0, -0.05D0,
     $             1.00D0, -0.04D0,
     $             1.00D0, -0.07D0,
     $             1.00D0, -0.05D0,
     $             1.00D0,  0.00D0/

      end block data lookup_data
***************************************************************************


***************************************************************************
       SUBROUTINE PHOTO_CAL(
     $        V2S,      V3S,
     $        U0,       TEMP,     RELO3,
     $        FACT,     PRESS,                        ! jjb PRESS added (needed in SR photo_cal)
     $        RJ_O3P,   RJ_O1D,   RJ_NO2,  RJ_HNO3,
     $        RJ_COH2,  RJ_CHOH,  RJ_N2O5, RJ_HNO4,
     $        RJ_NO2O,  RJ_NOO2,  RJ_H2O2, RJ_CH3OOH,
     $        RJ_O2,    RJ_CFC11, RJ_CFC12,RJ_N2O,
     $        RJ_ClONO2,RJ_BrNO3, RJ_Cl2O2, RJ_HOCl,
     $        RJ_BrCl_noT, RJ_ClNO2, RJ_BrNO2, RJ_Br2,
     $        RJ_IO,    RJ_INO3,  RJ_CH3I, RJ_I2,
     $        RJ_ICl,   RJ_IBr,   RJ_C3H7I, RJ_CH2ClI,
     $        RJ_CH2I2, RJ_INO2,  RJ_BrO_noT, RJ_OClO_noT,
     $        RJ_Cl2_noT, RJ_HOI_jen91, RJ_HOBr, RJ_HONO,
     $        RJ_NO2m,  RJ_NO3n,  RJ_dumm23, RJ_dumm24,
     $        RJ_dumm25, RJ_dumm26,
     $        H_O2,     H_O3,     H_O4,    QYNO3n )


! 13/08/2016 jjb comments before release
!     changes done in this SR:
!        - removal of unused argument
!        - USE global_params for vertical grid parameters
!        - removal of 3 unused parameters
!        - removal of unused COMMON (HALOE, Fo, TAU_TOP)
!        - commented unused variables


      USE global_params, ONLY :
! Imported Parameters:
     &     nrlay

      USE precision, ONLY :
! Imported Parameters:
     &     dp

      IMPLICIT NONE

      INTEGER MAXLAY,NW
      PARAMETER(MAXLAY=nrlay, NW=7)
*--------------------------------------------------------------------*
*     LOOK-UP TABLE                                                  *
*--------------------------------------------------------------------*
      COMMON/FIL/ IFIL(115)
      INTEGER IFIL

      COMMON/F_O2 /
     $              CS_O2(58,3,2),     FS_O2(58,2,2),
     $              CS_O3(58,2,2),     FS_O3(58,2,2),
     $              CS_N2O(58,2,2),    FS_N2O(58,2,2),
     $              CS_CFC11(58,3,2),  FS_CFC11(58,2,2),
     $              CS_CFC12(58,3,2),  FS_CFC12(58,2,2),
     $              CS_H2O2(58,2,2),   FS_H2O2(58,2,2),
     $              CS_HNO3(58,2,2),   FS_HNO3(58,2,2),
     $              CS_HNO4(58,2,2),   FS_HNO4(58,2,2),
     $              CS_ClONO2(58,2,2), FS_ClONO2(58,1,2),
     $              CS_BrNO3(58,2,2),
     $              CS_Cl2O2(58,2,2),
     $              CS_HOCl(58,2,2),
     $              CS_N2O5(58,2,2),   FS_N2O5(58,1,2),
     $              CS_HO2(58,3,2),    FS_HO2(58,2,2),
     $              CS_HO3(58,2,2),    FS_HO3(58,2,2),
     $              CS_BrCl_noT(58,3,2),
     $              CS_ClNO2(58,3,2),
     $              CS_BrNO2(58,3,2),
     $              CS_Br2(58,2,2),
     $              CS_CH3I(58,2,2),
     $              CS_NO3n(58,3,2),   FS_NO3n(58,2,2),
     $              CS_dumm23(58,2,2),
     $              CS_dumm24(58,2,2),
     $              CS_dumm25(58,2,2),
     $              CS_dumm26(58,2,2)
      REAL (KIND=DP) ::
     $     CS_O2, FS_O2, CS_O3, FS_O3, CS_N2O, FS_N2O,
     $     CS_CFC11, FS_CFC11, CS_CFC12, FS_CFC12,
     $     CS_H2O2, FS_H2O2, CS_HNO3, FS_HNO3, CS_HNO4, FS_HNO4,
     $     CS_ClONO2, FS_ClONO2, CS_BrNO3, CS_Cl2O2, CS_HOCl,
     $     CS_N2O5, FS_N2O5, CS_HO2, FS_HO2, CS_HO3, FS_HO3,
     $     CS_BrCl_noT, CS_ClNO2, CS_BrNO2, CS_Br2, CS_CH3I,
     $     CS_NO3n, FS_NO3n, CS_dumm23, CS_dumm24, CS_dumm25, CS_dumm26

      COMMON/LOOK/
     $          TAUA1(55,3),     TAUB1(55,3),
     $          A1_O3(55,3),     B1_O3(55,3),
     $          A1_O2(55,2),     B1_O2(55,2),
     $          A1_H2O2(55,2),   B1_H2O2(55,2),
     $          A1_HNO3(55,2),   B1_HNO3(55,2),
     $          A1_HNO4(55,2),   B1_HNO4(55,2),
     $          A1_N2O5(55,2),   B1_N2O5(55,2),
     $          A1_CH3OOH(55,2), B1_CH3OOH(55,2),
     $          A1_N2O(55,3),    B1_N2O(55,3),
     $          A1_CFC11(55,2),  B1_CFC11(55,2),
     $          A1_CFC12(55,3),  B1_CFC12(55,3),
     $          A1_ClONO2(55,2), B1_ClONO2(55,2),
     $          A1_BrNO3(55,2),  B1_BrNO3(55,2),
     $          A1_Cl2O2(55,3),  B1_Cl2O2(55,3),
     $          A1_HOCl(55,2),   B1_HOCl(55,2),
     $          A1_H_O2(55,2),   B1_H_O2(55,2),
     $          A1_H_O3(55,2),   B1_H_O3(55,2),
     $          A1_BrCl_noT(55,2),  B1_BrCl_noT(55,2),
     $          A1_ClNO2(55,2),  B1_ClNO2(55,2),
     $          A1_BrNO2(55,2),  B1_BrNO2(55,2),
     $          A1_Br2(55,2),    B1_Br2(55,2),
     $          A1_CH3I(55,2),   B1_CH3I(55,2),
     $          A1_ICl(55,2),    B1_ICl(55,2),
     $          A1_IBr(55,2),    B1_IBr(55,2),
     $          A1_C3H7I(55,2),  B1_C3H7I(55,2),
     $          A1_CH2ClI(55,2), B1_CH2ClI(55,2),
     $          A1_CH2I2(55,2),  B1_CH2I2(55,2),
     $          A1_INO2(55,2),   B1_INO2(55,2),
     $          A1_Cl2_noT(55,2),B1_Cl2_noT(55,2),
     $          A1_NO3n(55,3),   B1_NO3n(55,3),
     $          A1_dumm23(55,2), B1_dumm23(55,2),
     $          A1_dumm24(55,2), B1_dumm24(55,2),
     $          A1_dumm25(55,2), B1_dumm25(55,2),
     $          A1_dumm26(55,2), B1_dumm26(55,2),
     $          TAUA2(55),       TAUB2(55),
     $          A2_TO3(55),      B2_TO3(55),
     $          A2_O1D(55),      B2_O1D(55),
     $          A2_O3P(55),      B2_O3P(55),
     $          A2_H2O2(55),     B2_H2O2(55),
     $          A2_HNO3(55),     B2_HNO3(55),
     $          A2_HNO4(55),     B2_HNO4(55),
     $          A2_N2O5(55),     B2_N2O5(55),
     $          A2_CH3OOH(55),   B2_CH3OOH(55),
     $          A2_NO2(55),      B2_NO2(55),
     $          A2_ClONO2(55),   B2_ClONO2(55),
     $          A2_BrNO3(55),    B2_BrNO3(55),
     $          A2_Cl2O2(55),    B2_Cl2O2(55),
     $          A2_HOCl(55),     B2_HOCl(55),
     $          A2_H_O3(55),     B2_H_O3(55),
     $          A2_H_NO2(55),    B2_H_NO2(55),
     $          A2_BrCl_noT(55), B2_BrCl_noT(55),
     $          A2_ClNO2(55),    B2_ClNO2(55),
     $          A2_BrNO2(55),    B2_BrNO2(55),
     $          A2_Br2(55),      B2_Br2(55),
     $          A2_INO3(55),     B2_INO3(55),
     $          A2_CH3I(55),     B2_CH3I(55),
     $          A2_ICl(55),      B2_ICl(55),
     $          A2_IBr(55),      B2_IBr(55),
     $          A2_C3H7I(55),    B2_C3H7I(55),
     $          A2_CH2ClI(55),   B2_CH2ClI(55),
     $          A2_CH2I2(55),    B2_CH2I2(55),
     $          A2_INO2(55),     B2_INO2(55),
     $          A2_OClO_noT(55), B2_OClO_noT(55),
     $          A2_Cl2_noT(55),  B2_Cl2_noT(55),
     $          A2_HOBr(55),     B2_HOBr(55),
     $          A2_NO3n(55),     B2_NO3n(55),
     $          A2_dumm23(55),   B2_dumm23(55),
     $          A2_dumm24(55),   B2_dumm24(55),
     $          A2_dumm25(55),   B2_dumm25(55),
     $          A2_dumm26(55),   B2_dumm26(55),
     $          TAUA3(55),       TAUB3(55),
     $          A3_TO3(55),      B3_TO3(55),
     $          A3_O1D(55),      B3_O1D(55),
     $          A3_O3P(55),      B3_O3P(55),
     $          A3_H2O2(55),     B3_H2O2(55),
     $          A3_HNO3(55),     B3_HNO3(55),
     $          A3_HNO4(55),     B3_HNO4(55),
     $          A3_N2O5(55),     B3_N2O5(55),
     $          A3_CH3OOH(55),   B3_CH3OOH(55),
     $          A3_NO2(55),      B3_NO2(55),
     $          A3_COH2(55),     B3_COH2(55),
     $          A3_CHOH(55),     B3_CHOH(55),
     $          A3_ClONO2(55),   B3_ClONO2(55),
     $          A3_BrNO3(55),    B3_BrNO3(55),
     $          A3_Cl2O2(55),    B3_Cl2O2(55),
     $          A3_HOCl(55),     B3_HOCl(55),
     $          A3_H_O3(55),     B3_H_O3(55),
     $          A3_H_NO2(55),    B3_H_NO2(55),
     $          A3_BrCl_noT(55), B3_BrCl_noT(55),
     $          A3_ClNO2(55),    B3_ClNO2(55),
     $          A3_BrNO2(55),    B3_BrNO2(55),
     $          A3_Br2(55),      B3_Br2(55),
     $          A3_INO3(55),     B3_INO3(55),
     $          A3_CH3I(55),     B3_CH3I(55),
     $          A3_ICl(55),      B3_ICl(55),
     $          A3_IBr(55),      B3_IBr(55),
     $          A3_C3H7I(55),    B3_C3H7I(55),
     $          A3_CH2ClI(55),   B3_CH2ClI(55),
     $          A3_CH2I2(55),    B3_CH2I2(55),
     $          A3_INO2(55),     B3_INO2(55),
     $          A3_BrO_noT(55),  B3_BrO_noT(55),
     $          A3_OClO_noT(55), B3_OClO_noT(55),
     $          A3_Cl2_noT(55),  B3_Cl2_noT(55),
     $          A3_HOI_jen91(55),B3_HOI_jen91(55),
     $          A3_HOBr(55),     B3_HOBr(55),
     $          A3_NO2m(55),     B3_NO2m(55),
     $          A3_NO3n(55),     B3_NO3n(55),
     $          A3_dumm23(55),   B3_dumm23(55),
     $          A3_dumm24(55),   B3_dumm24(55),
     $          A3_dumm25(55),   B3_dumm25(55),
     $          A3_dumm26(55),    B3_dumm26(55)
      REAL (KIND=DP) ::
     $     TAUA1, TAUB1, A1_O3, B1_O3, A1_O2, B1_O2, A1_H2O2, B1_H2O2,
     $     A1_HNO3, B1_HNO3, A1_HNO4, B1_HNO4, A1_N2O5, B1_N2O5,
     $     A1_CH3OOH, B1_CH3OOH, A1_N2O, B1_N2O, A1_CFC11, B1_CFC11,
     $     A1_CFC12, B1_CFC12, A1_ClONO2, B1_ClONO2, A1_BrNO3, B1_BrNO3,
     $     A1_Cl2O2, B1_Cl2O2, A1_HOCl, B1_HOCl, A1_H_O2, B1_H_O2,
     $     A1_H_O3, B1_H_O3, A1_BrCl_noT, B1_BrCl_noT,
     $     A1_ClNO2, B1_ClNO2, A1_BrNO2, B1_BrNO2, A1_Br2, B1_Br2,
     $     A1_CH3I, B1_CH3I, A1_ICl, B1_ICl, A1_IBr, B1_IBr,
     $     A1_C3H7I, B1_C3H7I, A1_CH2ClI, B1_CH2ClI, A1_CH2I2, B1_CH2I2,
     $     A1_INO2, B1_INO2, A1_Cl2_noT, B1_Cl2_noT, A1_NO3n, B1_NO3n,
     $     A1_dumm23, B1_dumm23, A1_dumm24, B1_dumm24,
     $     A1_dumm25, B1_dumm25, A1_dumm26, B1_dumm26,
     $     TAUA2, TAUB2,
     $     A2_TO3, B2_TO3, A2_O1D, B2_O1D, A2_O3P, B2_O3P,
     $     A2_H2O2, B2_H2O2, A2_HNO3, B2_HNO3, A2_HNO4, B2_HNO4,
     $     A2_N2O5, B2_N2O5, A2_CH3OOH, B2_CH3OOH, A2_NO2, B2_NO2,
     $     A2_ClONO2, B2_ClONO2, A2_BrNO3, B2_BrNO3, A2_Cl2O2, B2_Cl2O2,
     $     A2_HOCl, B2_HOCl, A2_H_O3, B2_H_O3, A2_H_NO2, B2_H_NO2,
     $     A2_BrCl_noT, B2_BrCl_noT, A2_ClNO2, B2_ClNO2,
     $     A2_BrNO2, B2_BrNO2, A2_Br2, B2_Br2, A2_INO3, B2_INO3,
     $     A2_CH3I, B2_CH3I, A2_ICl, B2_ICl, A2_IBr, B2_IBr,
     $     A2_C3H7I, B2_C3H7I, A2_CH2ClI, B2_CH2ClI, A2_CH2I2, B2_CH2I2,
     $     A2_INO2, B2_INO2, A2_OClO_noT, B2_OClO_noT,
     $     A2_Cl2_noT, B2_Cl2_noT, A2_HOBr, B2_HOBr, A2_NO3n, B2_NO3n,
     $     A2_dumm23, B2_dumm23, A2_dumm24, B2_dumm24,
     $     A2_dumm25, B2_dumm25, A2_dumm26, B2_dumm26,
     $     TAUA3, TAUB3,
     $     A3_TO3, B3_TO3, A3_O1D, B3_O1D, A3_O3P, B3_O3P,
     $     A3_H2O2, B3_H2O2, A3_HNO3, B3_HNO3, A3_HNO4, B3_HNO4,
     $     A3_N2O5, B3_N2O5, A3_CH3OOH, B3_CH3OOH, A3_NO2, B3_NO2,
     $     A3_COH2, B3_COH2, A3_CHOH, B3_CHOH, A3_ClONO2, B3_ClONO2,
     $     A3_BrNO3, B3_BrNO3, A3_Cl2O2, B3_Cl2O2, A3_HOCl, B3_HOCl,
     $     A3_H_O3, B3_H_O3, A3_H_NO2, B3_H_NO2,
     $     A3_BrCl_noT, B3_BrCl_noT, A3_ClNO2, B3_ClNO2,
     $     A3_BrNO2, B3_BrNO2, A3_Br2, B3_Br2, A3_INO3, B3_INO3,
     $     A3_CH3I, B3_CH3I, A3_ICl, B3_ICl, A3_IBr, B3_IBr,
     $     A3_C3H7I, B3_C3H7I, A3_CH2ClI, B3_CH2ClI, A3_CH2I2, B3_CH2I2,
     $     A3_INO2, B3_INO2, A3_BrO_noT, B3_BrO_noT,
     $     A3_OClO_noT, B3_OClO_noT, A3_Cl2_noT, B3_Cl2_noT,
     $     A3_HOI_jen91, B3_HOI_jen91,
     $     A3_HOBr, B3_HOBr, A3_NO2m, B3_NO2m, A3_NO3n, B3_NO3n,
     $     A3_dumm23, B3_dumm23, A3_dumm24, B3_dumm24,
     $     A3_dumm25, B3_dumm25, A3_dumm26, B3_dumm26

      COMMON/C_POLY/
     $   C4_TAU(4),     C5_TAU(3),                    C7_TAU(2),
     $   C4_T_O3(4),    C5_T_O3(3),                   C7_T_O3(2),
     $   C4_O1D(4),     C5_O1D(4),
     $   C4_O3P(4),     C5_O3P(3),     C6_O3P(2),     C7_O3P(2),
     $   C4_H2O2(3),    C5_H2O2(3),    C6_H2O2(2),
     $   C4_HNO3(3),    C5_HNO3(3),    C6_HNO3(2),
     $   C4_HNO4(3),    C5_HNO4(3),
     $   C4_N2O5(3),    C5_N2O5(3),    C6_N2O5(2),
     $   C4_CH3OOH(3),  C5_CH3OOH(3),  C6_CH3OOH(2),
     $   C4_NO2(3),     C5_NO2(3),     C6_NO2(2),
     $   C4_COH2(4),    C5_COH2(3),    C6_COH2(2),
     $   C4_CHOH(4),    C5_CHOH(3),    C6_CHOH(2),
     $                                 C6_NO2O(2),    C7_NO2O(2),
     $                                                C7_NOO2(2),
     $   C4_ClONO2(3),  C5_ClONO2(3),  C6_ClONO2(1),
     $   C4_BrNO3(2),   C5_BrNO3(2),   C6_BrNO3(1),   C7_BrNO3(2),
     $   C4_Cl2O2(2),   C5_Cl2O2(2),   C6_Cl2O2(1),
     $   C4_HOCl(2),    C5_HOCl(2),    C6_HOCl(1),
     $   C4_H_O3(4),    C5_H_O3(3),    C6_H_O3(2),    C7_H_O3(2),
     $                                                C7_H_O2(2),
     $   C4_H_NO2(3),   C5_H_NO2(2),   C6_H_NO2(2),   C7_H_NO2(2),
     $                                 C6_H_O4(2),    C7_H_O4(2),
     $   C4_BrCl_noT(3),C5_BrCl_noT(3),C6_BrCl_noT(1),C7_BrCl_noT(2),
     $   C4_ClNO2(2),   C5_ClNO2(2),   C6_ClNO2(1),
     $   C4_BrNO2(3),   C5_BrNO2(2),   C6_BrNO2(1),   C7_BrNO2(2),
     $   C4_Br2(3),     C5_Br2(3),     C6_Br2(1),     C7_Br2(2),
     $                                 C6_IO(1),      C7_IO(2),
     $   C4_INO3(3),    C5_INO3(2),    C6_INO3(1),    C7_INO3(2),
     $   C4_CH3I(3),    C5_CH3I(3),    C6_CH3I(1),
     $                                 C6_I2(1),      C7_I2(2),
     $   C4_ICl(2),     C5_ICl(2),     C6_ICl(1),     C7_ICl(2),
     $   C4_IBr(3),     C5_IBr(3),     C6_IBr(1),     C7_IBr(2),
     $   C4_C3H7I(3),   C5_C3H7I(3),   C6_C3H7I(1),
     $   C4_CH2ClI(3),  C5_CH2ClI(3),  C6_CH2ClI(1),
     $   C4_CH2I2(2),   C5_CH2I2(2),   C6_CH2I2(1),
     $   C4_INO2(3),    C5_INO2(2),    C6_INO2(1),
     $   C4_BrO_noT(3), C5_BrO_noT(2), C6_BrO_noT(1),
     $   C4_OClO_noT(3), C5_OClO_noT(3), C6_OClO_noT(1), C7_OClO_noT(2),
     $   C4_Cl2_noT(3), C5_Cl2_noT(2), C6_Cl2_noT(1), C7_Cl2_noT(2),
     $   C4_HOI_jen91(3),C5_HOI_jen91(3),C6_HOI_jen91(1),
     $   C7_HOI_jen91(2),
     $   C4_HOBr(3),   C5_HOBr(2),   C6_HOBr(1),   C7_HOBr(2),
     $   C4_HONO(3),   C5_HONO(3),   C6_HONO(1),
     $   C4_NO2m(2),   C5_NO2m(2),   C6_NO2m(1),
     $   C4_NO3n(4),   C5_NO3n(4),   C6_NO3n(2),
     $   C4_dumm23(2),   C5_dumm23(2),   C6_dumm23(1),   C7_dumm23(2),
     $   C4_dumm24(2),   C5_dumm24(2),   C6_dumm24(1),   C7_dumm24(2),
     $   C4_dumm25(2),   C5_dumm25(2),   C6_dumm25(1),   C7_dumm25(2),
     $   C4_dumm26(2),   C5_dumm26(2),   C6_dumm26(1),   C7_dumm26(2)
      REAL (KIND=DP)
     $   C4_TAU,     C5_TAU,                 C7_TAU,
     $   C4_T_O3,    C5_T_O3,                C7_T_O3,
     $   C4_O1D,     C5_O1D,
     $   C4_O3P,     C5_O3P,     C6_O3P,     C7_O3P,
     $   C4_H2O2,    C5_H2O2,    C6_H2O2,
     $   C4_HNO3,    C5_HNO3,    C6_HNO3,
     $   C4_HNO4,    C5_HNO4,
     $   C4_N2O5,    C5_N2O5,    C6_N2O5,
     $   C4_CH3OOH,  C5_CH3OOH,  C6_CH3OOH,
     $   C4_NO2,     C5_NO2,     C6_NO2,
     $   C4_COH2,    C5_COH2,    C6_COH2,
     $   C4_CHOH,    C5_CHOH,    C6_CHOH,
     $                           C6_NO2O,    C7_NO2O,
     $                                       C7_NOO2,
     $   C4_ClONO2,  C5_ClONO2,  C6_ClONO2,
     $   C4_BrNO3,   C5_BrNO3,   C6_BrNO3,   C7_BrNO3,
     $   C4_Cl2O2,   C5_Cl2O2,   C6_Cl2O2,
     $   C4_HOCl,    C5_HOCl,    C6_HOCl,
     $   C4_H_O3,    C5_H_O3,    C6_H_O3,    C7_H_O3,
     $                                       C7_H_O2,
     $   C4_H_NO2,   C5_H_NO2,   C6_H_NO2,   C7_H_NO2,
     $                           C6_H_O4,    C7_H_O4,
     $   C4_BrCl_noT,C5_BrCl_noT,C6_BrCl_noT,C7_BrCl_noT,
     $   C4_ClNO2,   C5_ClNO2,   C6_ClNO2,
     $   C4_BrNO2,   C5_BrNO2,   C6_BrNO2,   C7_BrNO2,
     $   C4_Br2,     C5_Br2,     C6_Br2,     C7_Br2,
     $                           C6_IO,      C7_IO,
     $   C4_INO3,    C5_INO3,    C6_INO3,    C7_INO3,
     $   C4_CH3I,    C5_CH3I,    C6_CH3I,
     $                           C6_I2,      C7_I2,
     $   C4_ICl,     C5_ICl,     C6_ICl,     C7_ICl,
     $   C4_IBr,     C5_IBr,     C6_IBr,     C7_IBr,
     $   C4_C3H7I,   C5_C3H7I,   C6_C3H7I,
     $   C4_CH2ClI,  C5_CH2ClI,  C6_CH2ClI,
     $   C4_CH2I2,   C5_CH2I2,   C6_CH2I2,
     $   C4_INO2,    C5_INO2,    C6_INO2,
     $   C4_BrO_noT, C5_BrO_noT, C6_BrO_noT,
     $   C4_OClO_noT,C5_OClO_noT,C6_OClO_noT,C7_OClO_noT,
     $   C4_Cl2_noT, C5_Cl2_noT, C6_Cl2_noT, C7_Cl2_noT,
     $   C4_HOI_jen91,C5_HOI_jen91,C6_HOI_jen91,C7_HOI_jen91,
     $   C4_HOBr,    C5_HOBr,    C6_HOBr,   C7_HOBr,
     $   C4_HONO,    C5_HONO,    C6_HONO,
     $   C4_NO2m,    C5_NO2m,    C6_NO2m,
     $   C4_NO3n,    C5_NO3n,    C6_NO3n,
     $   C4_dumm23,  C5_dumm23,  C6_dumm23,  C7_dumm23,
     $   C4_dumm24,  C5_dumm24,  C6_dumm24,  C7_dumm24,
     $   C4_dumm25,  C5_dumm25,  C6_dumm25,  C7_dumm25,
     $   C4_dumm26,  C5_dumm26,  C6_dumm26,  C7_dumm26


      COMMON /T_COEFF/   TJ_O1D(4,7),   TJ_O3P(3,7),  TJ_NO2(2,7),
     $    TJ_HNO3(4,7),  TJ_N2O5(4,7),  TJ_H2O2(3,7), TJ_NOO2(2,7),
     $    TJ_NO2O(2,7),  TJ_N2O(3,7),   TJ_CFC12(4,7),TJ_CFC11(3,7),
     $    TJ_ClONO2(3,7),TJ_CHOH(2,7),  TJ_COH2(2,7), TJ_NO3n(3,6),
     $    TH_O3(2,7),    TH_O4(2,7),   TH_NO2(2,7)
      DOUBLE PRECISION TJ_O1D, TJ_O3P, TJ_NO2, TJ_HNO3, TJ_N2O5,TJ_H2O2,
     $     TJ_NOO2, TJ_NO2O, TJ_N2O,  TJ_CFC12, TJ_CFC11, TJ_ClONO2,
     $     TJ_CHOH, TJ_COH2, TJ_NO3n, TH_O3,    TH_O4,    TH_NO2

C     INPUT

      DOUBLE PRECISION
     $     TEMP(0:MAXLAY),   !temperature at model levels [K]
     $     PRESS(0:MAXLAY),  !pressure at model levels [hPa]
     $     RELO3(0:MAXLAY),  !O3 volumn mixing ratio [ppm]
     $     V2S(0:MAXLAY),    !slant O2 column density [part./cm^2]
     $     V3S(0:MAXLAY),    !slant O3 column density [part./cm^2]
     $     U0,               !cosine of solar zenith angle              ! jjb missing definition in this subroutine
     $     FACT(0:MAXLAY,NW) !act. flux from DISORT


C     output: J-values and heating rates

      DOUBLE PRECISION            !photolysis rates [1/s]
     $     RJ_O3P(0:MAXLAY),   RJ_O1D(0:MAXLAY),
     $     RJ_NO2(0:MAXLAY),   RJ_HNO3(0:MAXLAY),
     $     RJ_COH2(0:MAXLAY),  RJ_CHOH(0:MAXLAY),
     $     RJ_N2O5(0:MAXLAY),  RJ_HNO4(0:MAXLAY),
     $     RJ_NO2O(0:MAXLAY),  RJ_NOO2(0:MAXLAY),
     $     RJ_H2O2(0:MAXLAY),  RJ_CH3OOH(0:MAXLAY),
     $     RJ_O2(0:MAXLAY),    RJ_CFC11(0:MAXLAY),
     $     RJ_CFC12(0:MAXLAY), RJ_N2O(0:MAXLAY),
     $     RJ_ClONO2(0:MAXLAY),RJ_BrNO3(0:MAXLAY),
     $     RJ_Cl2O2(0:MAXLAY), RJ_HOCl(0:MAXLAY),
     $     RJ_BrCl_noT(0:MAXLAY), RJ_ClNO2(0:MAXLAY),
     $     RJ_BrNO2(0:MAXLAY), RJ_Br2(0:MAXLAY),
     $     RJ_IO(0:MAXLAY),    RJ_INO3(0:MAXLAY),
     $     RJ_CH3I(0:MAXLAY),  RJ_I2(0:MAXLAY),
     $     RJ_ICl(0:MAXLAY),   RJ_IBr(0:MAXLAY),
     $     RJ_C3H7I(0:MAXLAY), RJ_CH2ClI(0:MAXLAY),
     $     RJ_CH2I2(0:MAXLAY), RJ_INO2(0:MAXLAY),
     $     RJ_BrO_noT(0:MAXLAY), RJ_OClO_noT(0:MAXLAY),
     $     RJ_Cl2_noT(0:MAXLAY), RJ_HOI_jen91(0:MAXLAY),
     $     RJ_HOBr(0:MAXLAY),  RJ_HONO(0:MAXLAY),
     $     RJ_NO2m(0:MAXLAY),  RJ_NO3n(0:MAXLAY),
     $     RJ_dumm23(0:MAXLAY), RJ_dumm24(0:MAXLAY),
     $     RJ_dumm25(0:MAXLAY), RJ_dumm26(0:MAXLAY)

      DOUBLE PRECISION            !heating rates [K/s]
     $     H_O2(0:MAXLAY),       H_O3(0:MAXLAY),
     $     H_O4(0:MAXLAY)

C     internal variables

      DOUBLE PRECISION
     $     FINT(0:7),          !integrated actinic fluxes per interval
     $     V2S2,            !V2S  in the allowed range of int. 1
     $     V3S1,            !V3S  in the allowed range of int. 0-2
     $     V3S2,            !V3S  in the allowed range of int. 3-7
     $     DLV2,            !LOG(V2)
     $     V2S_m,           !slant O2 column [m/cm^2]
     $     V3_DU1,          !slant O3 column {DU] in int 0-2
     $     V3_DU2,          !slant O3 column {DU] in int 3-7
     $     TAU_0(MAXLAY),
     $     TAU_1,
     $     TAU_2,
     $     TAU_3,
     $     TAU_4,           !effective optical depths
     $     TAU_5,
     $     TAU_7

      INTEGER
     $     I0,
     $     I1,              !indices for lookup table
     $     I2,
     $     I3

      DOUBLE PRECISION
     $     QYNO3n(0:MAXLAY)    ! quantum yield for NO3-/NO3n


c---------------------------------------------------------------------
C     O2 AND O3 CROSS SECTIONS FOR T=250 K

      DOUBLE PRECISION CRS_O3(7)

      DATA  CRS_O3 /3.6448D-19,  1.8004D-18,  2.7700D-19,
     $              1.0500D-19,  2.6000D-20,  0.0000D+00,
     $              4.5500D-21/

      DOUBLE PRECISION CRS_O2(7)

      DATA  CRS_O2 /7.6300D-24,  0.0000D+00,  0.0000D+00,
     $              0.0000D+00,  0.0000D+00,  0.0000D+00,
     $              0.0000D+00/

C     effective cross sections

      DOUBLE PRECISION
!     $           SIG_O2(0:7,0:MAXLAY), SIG_O3(0:7,0:MAXLAY), ! jjb used up to index 1 only. SIG_03 definition up to index 7 still present, but unused
     $           SIG_O2(0:1,0:MAXLAY), SIG_O3(0:1,0:MAXLAY), ! jjb indexes reduced for efficiency
     $           SIG_CFC11(0:7),          SIG_CFC12(0:7),
     $           SIG_H2O2(0:7),           SIG_HNO3(0:7),
     $           SIG_HNO4(0:7),           SIG_N2O(0:7),
     $           SIG_CH3OOH(0:7),         SIG_H_O3(0:7),
     $           SIG_N2O5(0:7),           SIG_NO2(0:7),
     $           SIG_O1D(0:7),            SIG_O3P(0:7),
     $           SIG_COH2(0:7),           SIG_CHOH(0:7),
     $           SIG_NO2O(0:7),           SIG_NOO2(0:7),
     $           SIG_ClONO2(0:7),         SIG_BrNO3(0:7),
     $           SIG_Cl2O2(0:7),          SIG_HOCl(0:7),
     $           SIG_H_O2(0:7),!           SIG_H_NO2(0:7), ! jjb SIG_H_NO2 is no longer calculated
     $           SIG_H_O4(0:7),           SIG_BrCl_noT(0:7),
     $           SIG_ClNO2(0:7),          SIG_BrNO2(0:7),
     $           SIG_Br2(0:7),            SIG_IO(0:7),
     $           SIG_INO3(0:7),           SIG_CH3I(0:7),
     $           SIG_I2(0:7),             SIG_ICl(0:7),
     $           SIG_IBr(0:7),            SIG_C3H7I(0:7),
     $           SIG_CH2ClI(0:7),         SIG_CH2I2(0:7),
     $           SIG_INO2(0:7),           SIG_BrO_noT(0:7),
     $           SIG_OClO_noT(0:7),       SIG_Cl2_noT(0:7),
     $           SIG_HOI_jen91(0:7),      SIG_HOBr(0:7),
     $           SIG_HONO(0:7),           SIG_NO2m(0:7),
     $           SIG_NO3n(0:7),!           SIG_dumm23(0:7), ! jjb unreferenced
     $           SIG_dumm24(0:7),          SIG_dumm25(0:7),
     $           SIG_dumm26(0:7)

c     polynomial coeff. to calculate TAU_O above 0.1 hPa
      DOUBLE PRECISION
     $     CT(3)

      DATA CT /0.168306D0, 2.236551D0, 42.78577D0 /

      DOUBLE PRECISION
     $     F0(7)         ! integrated extraterrestric flux /
c                               extraterrestric flux (lambda(i))

      DATA F0 /7.3268D+01,  5.0944D+00,  1.6440D+01,  8.0198D+00,
     $         2.8000D+01,  1.8638D+01,  4.8184D+01/

C     change the units (Part./cm^2 -> Dobson Units)
C     Boltzmann constant K=1.38D-23 [J/K], normal conditions T0=273 [K],
C     P0=1000 [mbar], so CONST=K*T0/P0=3.767e-20 cm^3

      REAL (KIND=DP), PARAMETER :: CONSTANT = 3.767E-20_dp
      REAL (KIND=DP), PARAMETER :: BOLTZ = 1.381D-23
      REAL (KIND=DP), PARAMETER :: RELO2 = 0.2095D0

C     arrays to calculate TAU_0

      DOUBLE PRECISION
     $     P1, P2, P3, C0, C1, C2, C3, X, ! internal functions
     $     SIGI_O2,
     $     SIGI_O3,
     $     TEMP1(MAXLAY),
     $     TEMP2(MAXLAY)

      REAL (KIND=DP)
     $     DENS, DLV2I, DTAU_0,
     $     P2_O1D, P3_O1D, P4_O1D, P5_O1D,
     $     P2_O3P, P3_O3P, P4_O3P, P5_O3P, P6_O3P, P7_O3P,
     $     SO2, SO3, TAU3_LIM, TAU4_LIM,
     $     V3DU_3L, V3DU_4L, V3S_DU, V3_3L, V3_DU

      INTEGER
     $     I, II, K

c-------------------------------------------------------------------

**********************************************************************
C     internal functions

      P1(C0,C1,X)      = C0 + C1*X
      P2(C0,C1,C2,X)   = C0 + (C1 + C2*X)*X
      P3(C0,C1,C2,C3,X)= C0 + (C1 + (C2+C3*X)*X)*X

**********************************************************************


c      Do first the calculation of the optical depth TAU_0

C==================================================================
C         SCALING OF TEMPERATURE VARIABLE
C==================================================================
      DO K=1,MAXLAY
         TEMP1(K) = (TEMP(K)-240.D0)/240.D0 !for interval 0
         TEMP2(K) = (TEMP(K)-250.D0)/250.D0 !for interval 1-7
      ENDDO

c     ii.) first layer
      IF (U0.gt.0.) THEN

         DLV2I   = MIN(56.D0,DLOG(V2S(1)))
         V3S_DU = MIN(300.D0, V3S(1)* CONSTANT * 1.D+3)
         II      = MIN(58,INT(AINT((DLV2I-44.5D0)/0.2D0) + 1.00001D0))

         SIGI_O2=  P2(CS_O2(II,1,1)*DLV2I + CS_O2(II,1,2),
     $                   CS_O2(II,2,1)*DLV2I + CS_O2(II,2,2),
     $                   CS_O2(II,3,1)*DLV2I + CS_O2(II,3,2),
     $                   V3S_DU) *
     $                P2(1.D0,FS_O2(II,1,1)*DLV2I+FS_O2(II,1,2),
     $                      FS_O2(II,2,1)*DLV2I+FS_O2(II,2,2),
     $                      TEMP1(1))
         SIGI_O3=  P1(CS_O3(II,1,1)*DLV2I + CS_O3(II,1,2),
     $                   CS_O3(II,2,1)*DLV2I + CS_O3(II,2,2),
     $                   V3S_DU) *
     $                 P2(1.D0,FS_O3(II,1,1)*DLV2I+FS_O3(II,1,2),
     $                       FS_O3(II,2,1)*DLV2I+FS_O3(II,2,2),
     $                       TEMP1(1))

         TAU_0(1) = P2(CT(1),CT(2),CT(3),(DLV2I-47.D0)/47.D0)
      ENDIF

c     iii.) layers with pressure < 100 hPa

      DO K=2,MAXLAY
         IF (U0.gt.0.) THEN

               DLV2I   = MIN(56.D0,DLOG(V2S(K)))
               V3S_DU = MIN(300.D0, V3S(K)* CONSTANT * 1.D+3)
               II  = MIN(58,INT(AINT((DLV2I-44.5D0)/0.2D0) + 1.00001D0))

               SO2 =  P2(CS_O2(II,1,1)*DLV2I + CS_O2(II,1,2),
     $                   CS_O2(II,2,1)*DLV2I + CS_O2(II,2,2),
     $                   CS_O2(II,3,1)*DLV2I + CS_O2(II,3,2),
     $                   V3S_DU) *
     $                P2(1.D0,FS_O2(II,1,1)*DLV2I+FS_O2(II,1,2),
     $                      FS_O2(II,2,1)*DLV2I+FS_O2(II,2,2),
     $                      TEMP1(K))
               SO3 =  P1(CS_O3(II,1,1)*DLV2I + CS_O3(II,1,2),
     $                   CS_O3(II,2,1)*DLV2I + CS_O3(II,2,2),
     $                   V3S_DU) *
     $              P2(1.D0,FS_O3(II,1,1)*DLV2I+FS_O3(II,1,2),
     $                    FS_O3(II,2,1)*DLV2I+FS_O3(II,2,2),
     $                    TEMP1(K))

               DTAU_0  = 0.5D0*(SIGI_O2 + SO2) *
     $                       (V2S(K)-V2S(K-1)) +
     $                   0.5D0*(SIGI_O3 + SO3) *
     $                       (V3S(K)-V3S(K-1))

               TAU_0(K) = TAU_0(K-1) + DTAU_0

               SIGI_O2 = SO2
               SIGI_O3 = SO3

         ENDIF
      ENDDO


c------------------------------------------------------

c      BOUNDARIES FOR INTERVAL 3 and 4

      V3DU_3L = 1000.D0

      I  = MIN(INT((V3DU_3L-5D0)/25D0) +1,115)

      TAU3_LIM = P1(TAUB3(I),TAUA3(I),V3DU_3L)
      V3_3L =  V3DU_3L * 1.D-3 / CONSTANT

      V3DU_4L = 2500.D0
      TAU4_LIM = P2(C4_TAU(1),C4_TAU(2),C4_TAU(3),V3DU_4L)

c------------------------------------------------------

      DO 1111 K = 1,MAXLAY                ! altitude loop

      IF (U0.gt.0.D0) THEN

C        chance of units [V2S_m]=meter, [V3_DU]=DU

         V2S_m = V2S(K) * CONSTANT * 1.D-2
         V3_DU    = V3S(K) * CONSTANT * 1.D+3
         DLV2  = LOG(V2S(K))

C=================================================================
c         allowed ranges for columns
C=================================================================
c        interval: 0
         IF (DLV2.ge.56.D0) THEN
            DLV2  = 56.D0
         ENDIF

c        interval: 1
         IF (V2S_m.ge.500.D0) THEN
            V2S_m = 500.D0
            V2S2  = 1.3602D+24
         ELSE
            V2S2  = V2S(K)
         ENDIF

c        interval: 0 - 2
         IF (V3_DU.ge.300.D0) THEN
            V3_DU1= 300.D0
            V3S1  = 8.161D+18
         ELSE
            V3_DU1= V3_DU
            V3S1  = V3S(K)
         ENDIF

c        interval: 3 - 7
         IF (V3_DU.ge.3000.D0) THEN
            V3_DU2= 3000.D0
            V3S2  = 8.161D+19
         ELSE
            V3_DU2= V3_DU
            V3S2  = V3S(K)
         ENDIF

C==================================================================
c         Indices for lookup table
C==================================================================

         I0  = MIN(58,INT(AINT((DLV2-44.5D0)/0.2D0) + 1.00001D0))
         I1  = MIN(INT((V3_DU-0.5D0)/2.5D0) +1,115)
         I1  = IFIL(I1)
         I2  = MIN(INT((V3_DU-0.5D0)/2.5D0) +1,115)
         I2  = IFIL(I2)
         I3  = MIN(INT((V3_DU-5D0)/25D0) +1,115)
         I3  = IFIL(I3)
      ENDIF

      IF (U0.gt.0.D0) THEN

C==================================================================
C CALCULATION OF EFFECTIVE CROSS SECTIONS FOR OPTICAL DEPTH
C==================================================================
C-------------------
C O2 - O + O
C-------------------
         SIG_O2(0,K)=P2(CS_O2(I0,1,1)*DLV2 + CS_O2(I0,1,2),
     $                    CS_O2(I0,2,1)*DLV2 + CS_O2(I0,2,2),
     $                    CS_O2(I0,3,1)*DLV2 + CS_O2(I0,3,2),
     $                    V3_DU1) *
     $                 P2(1.D0,FS_O2(I0,1,1)*DLV2+
     $                         FS_O2(I0,1,2),
     $                         FS_O2(I0,2,1)*DLV2+
     $                         FS_O2(I0,2,2),
     $                         TEMP1(K))
        SIG_O2(1,K)=P1(B1_O2(I1,1),A1_O2(I1,1),V3_DU1) +
     $                 P1(B1_O2(I1,2),A1_O2(I1,2),V3_DU1) *
     $                    V2S_m
C--------------------
C O3
C--------------------
         SIG_O3(0,K)=P1(CS_O3(I0,1,1)*DLV2 + CS_O3(I0,1,2),
     $                    CS_O3(I0,2,1)*DLV2 + CS_O3(I0,2,2),
     $                    V3_DU1) *
     $                 P2(1.D0,
     $                    FS_O3(I0,1,1)*DLV2+FS_O3(I0,1,2),
     $                    FS_O3(I0,2,1)*DLV2+FS_O3(I0,2,2),
     $                    TEMP1(K))
         SIG_O3(1,K)=P1(B1_O3(I1,1),A1_O3(I1,1),V3_DU1) +
     $                 P1(B1_O3(I1,2),A1_O3(I1,2),V3_DU1) *
     $                    V2S_m +
     $                 P1(B1_O3(I1,3),A1_O3(I1,3),V3_DU1)*
     $                    V2S_m**2


      ENDIF

      IF (U0.gt.0.D0) THEN

c           calculation of optical depths above model atmosphere,
c           yields only in the range of the HALOE data

            TAU_1 = P1(TAUB1(I1,1),TAUA1(I1,1),V3_DU1) +
     $                 P1(TAUB1(I1,2),TAUA1(I1,2),V3_DU1) *
     $                    V2S_m +
     $                 P1(TAUB1(I1,3),TAUA1(I1,3),V3_DU1) *
     $                    V2S_m**2
            TAU_2 = P1(TAUB2(I2),TAUA2(I2),V3_DU1)
            TAU_3 = P1(TAUB3(I3),TAUA3(I3),V3_DU2)
            TAU_4 = P3(C4_TAU(1),C4_TAU(2),C4_TAU(3),C4_TAU(4),
     $                    V3_DU2)
            TAU_5 = P2(C5_TAU(1),C5_TAU(2),C5_TAU(3),V3_DU2)
            TAU_7 = P1(C7_TAU(1),C7_TAU(2),V3_DU2)

c        here 8.7990D+12 is the integrated flux over the SR bands at TOA

         FINT(0) = EXP(-TAU_0(K))* 8.7990D+12
         FINT(1) = EXP(-TAU_1 + CRS_O3(1)*V3S1 +
     $                    CRS_O2(1)*V2S2)*F0(1)*FACT(K,1)
         FINT(2) = EXP(-TAU_2 + CRS_O3(2)*V3S1)*F0(2)*
     $               FACT(K,2)
         IF (V3_DU2 .le. V3DU_3L) THEN
            FINT(3) = EXP(-TAU_3 + CRS_O3(3)*V3S2)*F0(3)*
     $                  FACT(K,3)
         ELSE
            FINT(3) = EXP(-TAU3_LIM + CRS_O3(3)*V3_3L)*F0(3)*
     $                  FACT(K,3)
            FINT(3) = 0.
         ENDIF
         IF (V3_DU2 .le. V3DU_4L) THEN
            FINT(4) = EXP(-TAU_4 + CRS_O3(4)*V3S2)*F0(4)*
     $                  FACT(K,4)
         ELSE
            FINT(4) = EXP(-TAU4_LIM + CRS_O3(4)*V3S2)*F0(4)*
     $                  FACT(K,4)
            FINT(4) = 0.
         ENDIF
         FINT(5) = EXP(-TAU_5 + CRS_O3(5)*V3S2)*F0(5)*
     $               FACT(K,5)
         FINT(6) = F0(6)*FACT(K,6)
         FINT(7) = EXP(-TAU_7 + CRS_O3(7)*V3S2)*F0(7)*
     $               FACT(K,7)

      ENDIF

      IF (U0.gt.0.D0) THEN

C==================================================================
C CALCULATION OF EFFECTIVE CROSS SECTIONS AND PHOTOLYSIS AND
C HEATING RATE CALCUTLATION
C==================================================================
          RJ_O2(K) = SIG_O2(0,K) * FINT(0) +
     $                    SIG_O2(1,K) * FINT(1)

C--------------------
C O3 -> O(1D) + O2
C--------------------
         P2_O1D     = P1(B2_O1D(I2),A2_O1D(I2),V3_DU1)
         SIG_O1D(2) = P3(TJ_O1D(1,2),TJ_O1D(2,2),TJ_O1D(3,2),
     $                   TJ_O1D(4,2),TEMP2(K)) * P2_O1D
         P3_O1D     = P1(B3_O1D(I3),A3_O1D(I3),V3_DU2)
         SIG_O1D(3) = P3(TJ_O1D(1,3),TJ_O1D(2,3),TJ_O1D(3,3),
     $                   TJ_O1D(4,3),TEMP2(K)) * P3_O1D
         P4_O1D     = P3(C4_O1D(1),C4_O1D(2),C4_O1D(3),C4_O1D(4),
     $                   V3_DU2)
         SIG_O1D(4) = P3(TJ_O1D(1,4),TJ_O1D(2,4),TJ_O1D(3,4),
     $                   TJ_O1D(4,4),TEMP2(K)) * P4_O1D
         P5_O1D     = P3(C5_O1D(1),C5_O1D(2),C5_O1D(3),C5_O1D(4),
     $                   V3_DU2)
         SIG_O1D(5) = P3(TJ_O1D(1,5),TJ_O1D(2,5),TJ_O1D(3,5),
     $                   TJ_O1D(4,5),TEMP2(K)) * P5_O1D

C         0.87 = fixed quantum yield for O3 -> O2 + O(1D)

          RJ_O1D(K) = 0.87D0*(SIG_O3(0,K) * FINT(0) +
     $                          SIG_O3(1,K) * FINT(1))+
     $                          SIG_O1D(2)          * FINT(2) +
     $                          SIG_O1D(3)          * FINT(3) +
     $                          SIG_O1D(4)          * FINT(4) +
     $                          SIG_O1D(5)          * FINT(5)

C--------------------
C O3 -> O(3P) + O2
C--------------------
          P2_O3P     = P1(B2_O3P(I2),A2_O3P(I2),V3_DU1)
          SIG_O3P(2) = P2(TJ_O3P(1,2),TJ_O3P(2,2),TJ_O3P(3,2),
     $                    TEMP2(K))* P2_O3P
          P3_O3P     = P1(B3_O3P(I3),A3_O3P(I3),V3_DU2)
          SIG_O3P(3) = P2(TJ_O3P(1,3),TJ_O3P(2,3),TJ_O3P(3,3),
     $                    TEMP2(K))* P3_O3P
          P4_O3P     = P3(C4_O3P(1),C4_O3P(2),C4_O3P(3),C4_O3P(4),
     $                 V3_DU2)
          SIG_O3P(4) = P2(TJ_O3P(1,4),TJ_O3P(2,4),TJ_O3P(3,4),
     $                    TEMP2(K)) * P4_O3P
          P5_O3P     = P2(C5_O3P(1),C5_O3P(2),C5_O3P(3),V3_DU2)
          SIG_O3P(5) = P2(TJ_O3P(1,5),TJ_O3P(2,5),TJ_O3P(3,5),
     $                    TEMP2(K)) * P5_O3P
          P6_O3P     = P1(C6_O3P(1),C6_O3P(2),V3_DU2)
          SIG_O3P(6) = P2(TJ_O3P(1,6),TJ_O3P(2,6),TJ_O3P(3,6),
     $                    TEMP2(K)) * P6_O3P
          P7_O3P     = P1(C7_O3P(1),C7_O3P(2),V3_DU2)
          SIG_O3P(7) = P2(TJ_O3P(1,7),TJ_O3P(2,7),TJ_O3P(3,7),
     $                    TEMP2(K)) * P7_O3P

          RJ_O3P(K) = 0.13D0*(SIG_O3(0,K) * FINT(0) +
     $                          SIG_O3(1,K) * FINT(1))+
     $                          SIG_O3P(2)          * FINT(2) +
     $                          SIG_O3P(3)          * FINT(3) +
     $                          SIG_O3P(4)          * FINT(4) +
     $                          SIG_O3P(5)          * FINT(5) +
     $                          SIG_O3P(6)          * FINT(6) +
     $                          SIG_O3P(7)          * FINT(7)
C--------------------
C N2O
C--------------------
          SIG_N2O(0)=P1(CS_N2O(I0,1,1)*DLV2 + CS_N2O(I0,1,2),
     $                  CS_N2O(I0,2,1)*DLV2 + CS_N2O(I0,2,2),
     $                  V3_DU1) *
     $               P2(1.D0,
     $                  FS_N2O(I0,1,1)*DLV2+FS_N2O(I0,1,2),
     $                  FS_N2O(I0,2,1)*DLV2+FS_N2O(I0,2,2),
     $                  TEMP1(K))
          SIG_N2O(1)=P2(TJ_N2O(1,1),TJ_N2O(2,1),TJ_N2O(3,1),
     $                  TEMP2(K))  *
     $              (P1(B1_N2O(I1,1),A1_N2O(I1,1),V3_DU1) +
     $               P1(B1_N2O(I1,2),A1_N2O(I1,2),V3_DU1) *
     $               V2S_m +
     $               P1(B1_N2O(I1,3),A1_N2O(I1,3),V3_DU1) *
     $               V2S_m**2)

          RJ_N2O(K) =    SIG_N2O(0)      * FINT(0) +
     $                     SIG_N2O(1)      * FINT(1)

C--------------------
C NO2
C--------------------
          SIG_NO2(2)=P1(TJ_NO2(1,2),TJ_NO2(2,2),TEMP2(K)) *
     $               P1(B2_NO2(I2),A2_NO2(I2),V3_DU1)
          SIG_NO2(3)=P1(TJ_NO2(1,3),TJ_NO2(2,3),TEMP2(K)) *
     $               P1(B3_NO2(I3),A3_NO2(I3),V3_DU2)
          SIG_NO2(4)=P1(TJ_NO2(1,4),TJ_NO2(2,4),TEMP2(K)) *
     $               P2(C4_NO2(1),C4_NO2(2),C4_NO2(3),V3_DU2)
          SIG_NO2(5)=P1(TJ_NO2(1,5),TJ_NO2(2,5),TEMP2(K)) *
     $               P2(C5_NO2(1),C5_NO2(2),C5_NO2(3),V3_DU2)
          SIG_NO2(6)=P1(TJ_NO2(1,5),TJ_NO2(2,5),TEMP2(K)) *
     $               P1(C6_NO2(1),C6_NO2(2),V3_DU2)

          RJ_NO2(K) = SIG_NO2(2)      * FINT(2) +
     $                     SIG_NO2(3)      * FINT(3) +
     $                     SIG_NO2(4)      * FINT(4) +
     $                     SIG_NO2(5)      * FINT(5) +
     $                     SIG_NO2(6)      * FINT(6)

C--------------------
C CFC-11
C--------------------
          SIG_CFC11(0)=P2(CS_CFC11(I0,1,1)*DLV2 +
     $                    CS_CFC11(I0,1,2),
     $                    CS_CFC11(I0,2,1)*DLV2 +
     $                    CS_CFC11(I0,2,2),
     $                    CS_CFC11(I0,3,1)*DLV2 +
     $                    CS_CFC11(I0,3,2),
     $                    V3_DU1) *
     $                 P2(1.D0,FS_CFC11(I0,1,1)*DLV2+
     $                       FS_CFC11(I0,1,2),
     $                       FS_CFC11(I0,2,1)*DLV2+
     $                       FS_CFC11(I0,2,2),TEMP1(K))
          SIG_CFC11(1)=P2(TJ_CFC11(1,1),TJ_CFC11(2,1),TJ_CFC11(3,1),
     $                    TEMP2(K)) *
     $                 P1(B1_CFC11(I1,1),A1_CFC11(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_CFC11(I1,2),A1_CFC11(I1,2),
     $                    V3_DU1)*V2S_m

          RJ_CFC11(K)=   SIG_CFC11(0)    * FINT(0) +
     $                     SIG_CFC11(1)    * FINT(1)
C--------------------
C CFC-12
C--------------------
          SIG_CFC12(0)=P2(CS_CFC12(I0,1,1)*DLV2 +
     $                    CS_CFC12(I0,1,2),
     $                    CS_CFC12(I0,2,1)*DLV2 +
     $                    CS_CFC12(I0,2,2),
     $                    CS_CFC12(I0,3,1)*DLV2 +
     $                    CS_CFC12(I0,3,2),
     $                    V3_DU1) *
     $                 P2(1.D0,FS_CFC12(I0,1,1)*DLV2 +
     $                       FS_CFC12(I0,1,2),
     $                       FS_CFC12(I0,2,1)*DLV2 +
     $                       FS_CFC12(I0,2,2),
     $                       TEMP1(K))
          SIG_CFC12(1)=P3(TJ_CFC12(1,1),TJ_CFC12(2,1),TJ_CFC12(3,1),
     $                    TJ_CFC12(4,1),TEMP2(K)) *
     $                (P1(B1_CFC12(I1,1),A1_CFC12(I1,1),
     $                    V3_DU1)+
     $                 P1(B1_CFC12(I1,2),A1_CFC12(I1,2),
     $                    V3_DU1)* V2S_m +
     $                 P1(B1_CFC12(I1,3),A1_CFC12(I1,3),
     $                    V3_DU1)* V2S_m**2)

          RJ_CFC12(K)=SIG_CFC12(0)    * FINT(0) +
     $                  SIG_CFC12(1)    * FINT(1)
      ENDIF

      IF (U0.gt.0.D0) THEN

C--------------------
C H2O2
C--------------------
          SIG_H2O2(0)=P1(CS_H2O2(I0,1,1)*DLV2 +
     $                   CS_H2O2(I0,1,2),
     $                   CS_H2O2(I0,2,1)*DLV2 +
     $                   CS_H2O2(I0,2,2),V3_DU1) *
     $                P2(1.D0,FS_H2O2(I0,1,1)*DLV2 +
     $                      FS_H2O2(I0,1,2),
     $                      FS_H2O2(I0,2,1)*DLV2 +
     $                      FS_H2O2(I0,2,2),TEMP1(K))
          SIG_H2O2(1)=P2(TJ_H2O2(1,1),TJ_H2O2(2,1),TJ_H2O2(3,1),
     $                   TEMP2(K))*
     $               (P1(B1_H2O2(I1,1),A1_H2O2(I1,1),V3_DU1)+
     $                P1(B1_H2O2(I1,2),A1_H2O2(I1,2),V3_DU1)*
     $                V2S_m)
          SIG_H2O2(2)=P2(TJ_H2O2(1,2),TJ_H2O2(2,2),TJ_H2O2(3,2),
     $                   TEMP2(K))*
     $                P1(B2_H2O2(I2),A2_H2O2(I2),V3_DU1)
          SIG_H2O2(3)=P2(TJ_H2O2(1,3),TJ_H2O2(2,3),TJ_H2O2(3,3),
     $                   TEMP2(K))*
     $                P1(B3_H2O2(I3),A3_H2O2(I3),V3_DU2)
          SIG_H2O2(4)=P2(TJ_H2O2(1,4),TJ_H2O2(2,4),TJ_H2O2(3,4),
     $                   TEMP2(K))*
     $                P2(C4_H2O2(1),C4_H2O2(2),C4_H2O2(3),V3_DU2)
          SIG_H2O2(5)=P2(TJ_H2O2(1,5),TJ_H2O2(2,5),TJ_H2O2(3,5),
     $                   TEMP2(K))*
     $                P2(C5_H2O2(1),C5_H2O2(2),C5_H2O2(3),V3_DU2)
          SIG_H2O2(6)=P2(TJ_H2O2(1,6),TJ_H2O2(2,6),TJ_H2O2(3,6),
     $                   TEMP2(K))*
     $                P1(C6_H2O2(1),C6_H2O2(2),V3_DU2)

          RJ_H2O2(K)= SIG_H2O2(0)     * FINT(0) +
     $                     SIG_H2O2(1)     * FINT(1) +
     $                     SIG_H2O2(2)     * FINT(2) +
     $                     SIG_H2O2(3)     * FINT(3) +
     $                     SIG_H2O2(4)     * FINT(4) +
     $                     SIG_H2O2(5)     * FINT(5) +
     $                     SIG_H2O2(6)     * FINT(6)
C--------------------
C HNO3
C--------------------
          SIG_HNO3(0)=P1(CS_HNO3(I0,1,1)*DLV2 +
     $                   CS_HNO3(I0,1,2),
     $                   CS_HNO3(I0,2,1)*DLV2 +
     $                   CS_HNO3(I0,2,2),V3_DU1) *
     $                P2(1.D0,FS_HNO3(I0,1,1)*DLV2 +
     $                      FS_HNO3(I0,1,2),
     $                      FS_HNO3(I0,2,1)*DLV2 +
     $                      FS_HNO3(I0,2,2),TEMP1(K))
          SIG_HNO3(1)=(P1(B1_HNO3(I1,1),A1_HNO3(I1,1),
     $                    V3_DU1)+
     $                P1(B1_HNO3(I1,2),A1_HNO3(I1,2),
     $                   V3_DU1)*V2S_m)*
     $                P3(TJ_HNO3(1,1),TJ_HNO3(2,1),TJ_HNO3(3,1),
     $                   TJ_HNO3(4,1),TEMP2(K))
          SIG_HNO3(2)=P3(TJ_HNO3(1,2),TJ_HNO3(2,2),TJ_HNO3(3,2),
     $                   TJ_HNO3(4,2),TEMP2(K))*
     $                P1(B2_HNO3(I2),A2_HNO3(I2),V3_DU1)
          SIG_HNO3(3)=P3(TJ_HNO3(1,3),TJ_HNO3(2,3),TJ_HNO3(3,3),
     $                   TJ_HNO3(4,3),TEMP2(K))*
     $                P1(B3_HNO3(I3),A3_HNO3(I3),V3_DU2)
          SIG_HNO3(4)=P3(TJ_HNO3(1,4),TJ_HNO3(2,4),TJ_HNO3(3,4),
     $                   TJ_HNO3(4,4),TEMP2(K))*
     $                P2(C4_HNO3(1),C4_HNO3(2),C4_HNO3(3),V3_DU2)
          SIG_HNO3(5)=P3(TJ_HNO3(1,5),TJ_HNO3(2,5),TJ_HNO3(3,5),
     $                   TJ_HNO3(4,5),TEMP2(K))*
     $                P2(C5_HNO3(1),C5_HNO3(2),C5_HNO3(3),V3_DU2)
          SIG_HNO3(6)=P3(TJ_HNO3(1,6),TJ_HNO3(2,6),TJ_HNO3(3,6),
     $                   TJ_HNO3(4,6),TEMP2(K))*
     $                P1(C6_HNO3(1),C6_HNO3(2),V3_DU2)

          RJ_HNO3(K)= SIG_HNO3(0)     * FINT(0) +
     $                     SIG_HNO3(1)     * FINT(1) +
     $                     SIG_HNO3(2)     * FINT(2) +
     $                     SIG_HNO3(3)     * FINT(3) +
     $                     SIG_HNO3(4)     * FINT(4) +
     $                     SIG_HNO3(5)     * FINT(5) +
     $                     SIG_HNO3(6)     * FINT(6)
C--------------------
C HNO4
C--------------------
          SIG_HNO4(0)=P1(CS_HNO4(I0,1,1)*DLV2 +
     $                   CS_HNO4(I0,1,2),
     $                   CS_HNO4(I0,2,1)*DLV2 +
     $                   CS_HNO4(I0,2,2),V3_DU1) *
     $                P2(1.D0,FS_HNO4(I0,1,1)*DLV2+
     $                      FS_HNO4(I0,1,2),
     $                      FS_HNO4(I0,2,1)*DLV2+
     $                      FS_HNO4(I0,2,2),TEMP1(K))
          SIG_HNO4(1)=P1(B1_HNO4(I1,1),A1_HNO4(I1,1),
     $                   V3_DU1) +
     $                P1(B1_HNO4(I1,2),A1_HNO4(I1,2),
     $                   V3_DU1) *V2S_m
          SIG_HNO4(2)=P1(B2_HNO4(I2), A2_HNO4(I2),V3_DU1)
          SIG_HNO4(3)=P1(B3_HNO4(I3),A3_HNO4(I3),V3_DU2)
          SIG_HNO4(4)=P2(C4_HNO4(1),C4_HNO4(2),C4_HNO4(3),V3_DU2)
          SIG_HNO4(5)=P2(C5_HNO4(1),C5_HNO4(2),C5_HNO4(3),V3_DU2)

          RJ_HNO4(K)= SIG_HNO4(0)     * FINT(0) +
     $                     SIG_HNO4(1)     * FINT(1) +
     $                     SIG_HNO4(2)     * FINT(2) +
     $                     SIG_HNO4(3)     * FINT(3) +
     $                     SIG_HNO4(4)     * FINT(4) +
     $                     SIG_HNO4(5)     * FINT(5)
C--------------------
C N2O5
C--------------------
          SIG_N2O5(0)=P1(CS_N2O5(I0,1,1)*DLV2+
     $                   CS_N2O5(I0,1,2),
     $                   CS_N2O5(I0,2,1)*DLV2+
     $                   CS_N2O5(I0,2,2),V3_DU1) *
     $                P1(1.D0,FS_N2O5(I0,1,1)*DLV2+
     $                        FS_N2O5(I0,1,2),TEMP1(K))
          SIG_N2O5(1)=P3(TJ_N2O5(1,1),TJ_N2O5(2,1),TJ_N2O5(3,1),
     $                   TJ_N2O5(4,1),TEMP2(K))*
     $               (P1(B1_N2O5(I1,1),A1_N2O5(I1,1),
     $                   V3_DU1) +
     $                P1(B1_N2O5(I1,2),A1_N2O5(I1,2),
     $                   V3_DU1) *V2S_m)
          SIG_N2O5(2)=P3(TJ_N2O5(1,2),TJ_N2O5(2,2),TJ_N2O5(3,2),
     $                   TJ_N2O5(4,2),TEMP2(K))*
     $                P1(B2_N2O5(I2),A2_N2O5(I2),V3_DU1)
          SIG_N2O5(3)=P3(TJ_N2O5(1,3),TJ_N2O5(2,3),TJ_N2O5(3,3),
     $                   TJ_N2O5(4,3),TEMP2(K))*
     $                P1(B3_N2O5(I3),A3_N2O5(I3),V3_DU2)
          SIG_N2O5(4)=P3(TJ_N2O5(1,4),TJ_N2O5(2,4),TJ_N2O5(3,4),
     $                   TJ_N2O5(4,4),TEMP2(K))*
     $                P2(C4_N2O5(1),C4_N2O5(2),C4_N2O5(3),V3_DU2)
          SIG_N2O5(5)=P3(TJ_N2O5(1,5),TJ_N2O5(2,5),TJ_N2O5(3,5),
     $                   TJ_N2O5(4,5),TEMP2(K))*
     $                P2(C5_N2O5(1),C5_N2O5(2),C5_N2O5(3),V3_DU2)
          SIG_N2O5(6)=P3(TJ_N2O5(1,6),TJ_N2O5(2,6),TJ_N2O5(3,6),
     $                   TJ_N2O5(4,6),TEMP2(K))*
     $                P1(C6_N2O5(1),C6_N2O5(2),V3_DU2)

          RJ_N2O5(K)=    SIG_N2O5(0)     * FINT(0) +
     $                     SIG_N2O5(1)     * FINT(1) +
     $                     SIG_N2O5(2)     * FINT(2) +
     $                     SIG_N2O5(3)     * FINT(3) +
     $                     SIG_N2O5(4)     * FINT(4) +
     $                     SIG_N2O5(5)     * FINT(5) +
     $                     SIG_N2O5(6)     * FINT(6)
C--------------------
C ClONO2
C--------------------
          SIG_ClONO2(0)=P1(CS_ClONO2(I0,1,1)*DLV2+
     $                     CS_ClONO2(I0,1,2),
     $                     CS_ClONO2(I0,2,1)*DLV2+
     $                     CS_ClONO2(I0,2,2),V3_DU1) *
     $                  P1(1.D0,FS_ClONO2(I0,1,1)*DLV2+
     $                     FS_ClONO2(I0,1,2),TEMP1(K))
          SIG_ClONO2(1)=P1(B1_ClONO2(I1,1),A1_ClONO2(I1,1),
     $                     V3_DU1) +
     $                  P1(B1_ClONO2(I1,2),A1_ClONO2(I1,2),
     $                     V3_DU1) * V2S_m
          SIG_ClONO2(2)=P2(TJ_ClONO2(1,2),TJ_ClONO2(2,2),
     $                     TJ_ClONO2(3,2),TEMP2(K))*
     $                  P1(B2_ClONO2(I2),A2_ClONO2(I2),
     $                     V3_DU1)
          SIG_ClONO2(3)=P2(TJ_ClONO2(1,3),TJ_ClONO2(2,3),
     $                     TJ_ClONO2(3,3),TEMP2(K))*
     $                  P1(B3_ClONO2(I3),A3_ClONO2(I3),
     $                     V3_DU2)
          SIG_ClONO2(4)=P2(TJ_ClONO2(1,4),TJ_ClONO2(2,4),
     $                     TJ_ClONO2(3,4),TEMP2(K))*
     $                  P2(C4_ClONO2(1),C4_ClONO2(2),C4_ClONO2(3),
     $                     V3_DU2)
          SIG_ClONO2(5)=P2(TJ_ClONO2(1,5),TJ_ClONO2(2,5),
     $                     TJ_ClONO2(3,5),TEMP2(K))*
     $                  P2(C5_ClONO2(1),C5_ClONO2(2),C5_ClONO2(3),
     $                     V3_DU2)
          SIG_ClONO2(6)=P2(TJ_ClONO2(1,6),TJ_ClONO2(2,6),
     $                     TJ_ClONO2(3,6),TEMP2(K))*C6_ClONO2(1)

          RJ_ClONO2(K)=SIG_ClONO2(0)    * FINT(0) +
     $                      SIG_ClONO2(1)    * FINT(1) +
     $                      SIG_ClONO2(2)    * FINT(2) +
     $                      SIG_ClONO2(3)    * FINT(3) +
     $                      SIG_ClONO2(4)    * FINT(4) +
     $                      SIG_ClONO2(5)    * FINT(5) +
     $                      SIG_ClONO2(6)    * FINT(6)
C--------------------
C BrONO2
C--------------------
          SIG_BrNO3(0)=P1(CS_BrNO3(I0,1,1)*DLV2 +
     $                    CS_BrNO3(I0,1,2),
     $                    CS_BrNO3(I0,2,1)*DLV2 +
     $                    CS_BrNO3(I0,2,2),V3_DU1)
          SIG_BrNO3(1)=P1(B1_BrNO3(I1,1),A1_BrNO3(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_BrNO3(I1,2),A1_BrNO3(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_BrNO3(2)=P1(B2_BrNO3(I2), A2_BrNO3(I2),V3_DU1)
          SIG_BrNO3(3)=P1(B3_BrNO3(I3), A3_BrNO3(I3),V3_DU2)
          SIG_BrNO3(4)=P1(C4_BrNO3(1),C4_BrNO3(2),V3_DU2)
          SIG_BrNO3(5)=P1(C5_BrNO3(1),C5_BrNO3(2),V3_DU2)
          SIG_BrNO3(6)=C6_BrNO3(1)
          SIG_BrNO3(7)=P1(C7_BrNO3(1),C7_BrNO3(2),V3_DU2)

          RJ_BrNO3(K)=   SIG_BrNO3(0)    * FINT(0) +
     $                     SIG_BrNO3(1)    * FINT(1) +
     $                     SIG_BrNO3(2)    * FINT(2) +
     $                     SIG_BrNO3(3)    * FINT(3) +
     $                     SIG_BrNO3(4)    * FINT(4) +
     $                     SIG_BrNO3(5)    * FINT(5) +
     $                     SIG_BrNO3(6)    * FINT(6) +
     $                     SIG_BrNO3(7)    * FINT(7)

C--------------------
C HOCl
C--------------------
          SIG_HOCl(0)=P1(CS_HOCl(I0,1,1)*DLV2 +
     $                   CS_HOCl(I0,1,2),
     $                   CS_HOCl(I0,2,1)*DLV2 +
     $                   CS_HOCl(I0,2,2),V3_DU1)
          SIG_HOCl(1)=P1(B1_HOCl(I1,1),A1_HOCl(I1,1),
     $                   V3_DU1) +
     $                P1(B1_HOCl(I1,2),A1_HOCl(I1,2),
     $                   V3_DU1) *V2S_m
          SIG_HOCl(2)=P1(B2_HOCl(I2),A2_HOCl(I2),V3_DU1)
          SIG_HOCl(3)=P1(B3_HOCl(I3),A3_HOCl(I3),V3_DU2)
          SIG_HOCl(4)=P1(C4_HOCl(1), C4_HOCl(2), V3_DU2)
          SIG_HOCl(5)=P1(C5_HOCl(1), C5_HOCl(2), V3_DU2)
          SIG_HOCl(6)=C6_HOCl(1)

          RJ_HOCl(K)= SIG_HOCl(0)     * FINT(0) +
     $                     SIG_HOCl(1)     * FINT(1) +
     $                     SIG_HOCl(2)     * FINT(2) +
     $                     SIG_HOCl(3)     * FINT(3) +
     $                     SIG_HOCl(4)     * FINT(4) +
     $                     SIG_HOCl(5)     * FINT(5) +
     $                     SIG_HOCl(6)     * FINT(6)
C--------------------
C CH2O -> CO + H2
C--------------------
          SIG_COH2(3)=P1(TJ_COH2(1,3),TJ_COH2(2,3),TEMP2(K))*
     $                P1(B3_COH2(I3),A3_COH2(I3),V3_DU2)
          SIG_COH2(4)=P3(C4_COH2(1),C4_COH2(2),C4_COH2(3),
     $                   C4_COH2(4),V3_DU2)
          SIG_COH2(5)=P1(TJ_COH2(1,5),TJ_COH2(2,5),TEMP2(K))*
     $                P2(C5_COH2(1),C5_COH2(2),C5_COH2(3),V3_DU2)
          SIG_COH2(6)=P1(TJ_COH2(1,6),TJ_COH2(2,6),TEMP2(K))*
     $                P1(C6_COH2(1),C6_COH2(2),V3_DU2)

          RJ_COH2(K)=
     $                     SIG_COH2(3)     * FINT(3) +
     $                     SIG_COH2(4)     * FINT(4) +
     $                     SIG_COH2(5)     * FINT(5) +
     $                     SIG_COH2(6)     * FINT(6)
C--------------------
C CH2O -> COH + H
C--------------------

          SIG_CHOH(3)=P1(TJ_CHOH(1,3),TJ_CHOH(2,3),TEMP2(K))*
     $                P1(B3_CHOH(I3),A3_CHOH(I3),V3_DU2)
          SIG_CHOH(4)=P3(C4_CHOH(1),C4_CHOH(2),C4_CHOH(3),
     $                   C4_CHOH(4),V3_DU2)
          SIG_CHOH(5)=P1(TJ_CHOH(1,5),TJ_CHOH(2,5),TEMP2(K))*
     $                P2(C5_CHOH(1),C5_CHOH(2),C5_CHOH(3),V3_DU2)
          SIG_CHOH(6)=P1(TJ_CHOH(1,6),TJ_CHOH(2,6),TEMP2(K))*
     $                P1(C6_CHOH(1),C6_CHOH(2),V3_DU2)

          RJ_CHOH(K)=
     $                     SIG_CHOH(3)     * FINT(3) +
     $                     SIG_CHOH(4)     * FINT(4) +
     $                     SIG_CHOH(5)     * FINT(5) +
     $                     SIG_CHOH(6)     * FINT(6)
C--------------------
C CH3OOH
C--------------------
          SIG_CH3OOH(1)=P1(B1_CH3OOH(I1,1),A1_CH3OOH(I1,1),
     $                     V3_DU1) +
     $                  P1(B1_CH3OOH(I1,2),A1_CH3OOH(I1,2),
     $                     V3_DU1) * V2S_m
          SIG_CH3OOH(2)=P1(B2_CH3OOH(I2),A2_CH3OOH(I2),
     $                     V3_DU1)
          SIG_CH3OOH(3)=P1(B3_CH3OOH(I3),A3_CH3OOH(I3),
     $                     V3_DU2)
          SIG_CH3OOH(4)=P2(C4_CH3OOH(1),C4_CH3OOH(2),C4_CH3OOH(3),
     $                     V3_DU2)
          SIG_CH3OOH(5)=P2(C5_CH3OOH(1),C5_CH3OOH(2),C5_CH3OOH(3),
     $                     V3_DU2)
          SIG_CH3OOH(6)=P1(C6_CH3OOH(1),C6_CH3OOH(2),V3_DU2)

          RJ_CH3OOH(K)=   SIG_CH3OOH(1) * FINT(1) +
     $                      SIG_CH3OOH(2) * FINT(2) +
     $                      SIG_CH3OOH(3) * FINT(3) +
     $                      SIG_CH3OOH(4) * FINT(4) +
     $                      SIG_CH3OOH(5) * FINT(5) +
     $                      SIG_CH3OOH(6) * FINT(6)
C--------------------
C Cl2O2
C--------------------
          SIG_Cl2O2(1)=P1(B1_Cl2O2(I1,1),A1_Cl2O2(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_Cl2O2(I1,2),A1_Cl2O2(I1,2),
     $                    V3_DU1) *V2S_m +
     $                 P1(B1_Cl2O2(I1,3),A1_Cl2O2(I1,3),
     $                    V3_DU1) *V2S_m**2
          SIG_Cl2O2(2)=P1(B2_Cl2O2(I2), A2_Cl2O2(I2),
     $                    V3_DU1)
          SIG_Cl2O2(3)=P1(B3_Cl2O2(I3), A3_Cl2O2(I3),
     $                    V3_DU2)
          SIG_Cl2O2(4)=P1(C4_Cl2O2(1),C4_Cl2O2(2),V3_DU2)
          SIG_Cl2O2(5)=P1(C5_Cl2O2(1),C5_Cl2O2(2),V3_DU2)
          SIG_Cl2O2(6)=C6_Cl2O2(1)

          RJ_Cl2O2(K)=   SIG_Cl2O2(1)    * FINT(1) +
     $                     SIG_Cl2O2(2)    * FINT(2) +
     $                     SIG_Cl2O2(3)    * FINT(3) +
     $                     SIG_Cl2O2(4)    * FINT(4) +
     $                     SIG_Cl2O2(5)    * FINT(5) +
     $                     SIG_Cl2O2(6)    * FINT(6)
C--------------------
C NO3 -> NO2 + O
C--------------------
          SIG_NO2O(6)=P1(TJ_NO2O(1,6),TJ_NO2O(2,6),TEMP2(K))*
     $                P1(C6_NO2O(1),C6_NO2O(2),V3_DU2)
          SIG_NO2O(7)=P1(TJ_NO2O(1,7),TJ_NO2O(2,7),TEMP2(K))*
     $                P1(C7_NO2O(1),C7_NO2O(2),V3_DU2)

          RJ_NO2O(K)=    SIG_NO2O(6)     * FINT(6) +
     $                     SIG_NO2O(7)     * FINT(7)
C--------------------
C NO3 -> NO  + O2
C--------------------
          SIG_NOO2(7)=P1(TJ_NOO2(1,7),TJ_NOO2(2,7),TEMP2(K)) *
     $                P1(C7_NOO2(1),C7_NOO2(2),V3_DU2)

          RJ_NOO2(K)= SIG_NOO2(7)     * FINT(7)
C--------------------
C O2 HEATING
C--------------------
          SIG_H_O2(0)=P2(CS_HO2(I0,1,1)*DLV2+CS_HO2(I0,1,2),
     $                   CS_HO2(I0,2,1)*DLV2+CS_HO2(I0,2,2),
     $                   CS_HO2(I0,3,1)*DLV2+CS_HO2(I0,3,2),
     $                   V3_DU1) *
     $                P2(1.D0,FS_HO2(I0,1,1)*DLV2+
     $                      FS_HO2(I0,1,2),
     $                      FS_HO2(I0,2,1)*DLV2+
     $                      FS_HO2(I0,2,2),TEMP1(K))
          SIG_H_O2(1)=P1(B1_H_O2(I1,1),A1_H_O2(I1,1),V3_DU1)+
     $                P1(B1_H_O2(I1,2),A1_H_O2(I1,2),V3_DU1)*
     $                V2S_m

          SIG_H_O2(7)=P1(C7_H_O2(1),C7_H_O2(2),V3_DU2)

          H_O2(K)  =(   SIG_H_O2(0)     * FINT(0) +
     $                    SIG_H_O2(1)     * FINT(1) +
     $                    SIG_H_O2(7)     * FINT(7) ) * RELO2
C--------------------
C O3 HEATING
C--------------------
          SIG_H_O3(0)=P1(CS_HO3(I0,1,1)*DLV2+ CS_HO3(I0,1,2),
     $                   CS_HO3(I0,2,1)*DLV2+ CS_HO3(I0,2,2),
     $                   V3_DU1) *
     $                P2(1.D0,FS_HO3(I0,1,1)*DLV2+
     $                      FS_HO3(I0,1,2),
     $                      FS_HO3(I0,2,1)*DLV2+
     $                      FS_HO3(I0,2,2),TEMP1(K))
          SIG_H_O3(1)=P1(B1_H_O3(I1,1),A1_H_O3(I1,1),V3_DU1)+
     $                P1(B1_H_O3(I1,2),A1_H_O3(I1,2),V3_DU1)*
     $                  V2S_m
          SIG_H_O3(2)=P1(TH_O3(1,2),TH_O3(2,2),TEMP2(K))*
     $                P1(B2_H_O3(I2),A2_H_O3(I2),V3_DU1)
          SIG_H_O3(3)=P1(TH_O3(1,3),TH_O3(2,3),TEMP2(K))*
     $                P1(B3_H_O3(I3),A3_H_O3(I3),V3_DU2)
          SIG_H_O3(4)=P1(TH_O3(1,4),TH_O3(2,4),TEMP2(K))*
     $                P3(C4_H_O3(1),C4_H_O3(2),C4_H_O3(3),C4_H_O3(4),
     $                V3_DU2)
          SIG_H_O3(5)=P1(TH_O3(1,5),TH_O3(2,5),TEMP2(K)) *
     $                P2(C5_H_O3(1),C5_H_O3(2),C5_H_O3(3),V3_DU2)
          SIG_H_O3(6)=P1(TH_O3(1,6),TH_O3(2,6),TEMP2(K)) *
     $                P1(C6_H_O3(1),C6_H_O3(2),V3_DU2)
          SIG_H_O3(7)=P1(TH_O3(1,7),TH_O3(2,7),TEMP2(K))*
     $                P1(C7_H_O3(1),C7_H_O3(2),V3_DU2)

          H_O3(K)   = (SIG_H_O3(0)    * FINT(0) +
     $                      SIG_H_O3(1)    * FINT(1) +
     $                      SIG_H_O3(2)    * FINT(2) +
     $                      SIG_H_O3(3)    * FINT(3) +
     $                      SIG_H_O3(4)    * FINT(4) +
     $                      SIG_H_O3(5)    * FINT(5) +
     $                      SIG_H_O3(6)    * FINT(6) +
     $                      SIG_H_O3(7)    * FINT(7) ) *
     $                      RELO3(K)*1.E6


C--------------------
C O2-O2 HEATING
C--------------------
c         DENS = density of moelcules [1/cm^3]
c         1.D-46: scaling factor of SIG_H_O4
C         1.D-6 because of: 1/m^3 -> 1/cm^3

          DENS = PRESS(K)/(BOLTZ*TEMP(K))*1.D-46*1.D-6

          SIG_H_O4(6)=P1(C6_H_O4(1),C6_H_O4(2),V3_DU2)
          SIG_H_O4(7)=P1(C7_H_O4(1),C7_H_O4(2),V3_DU2)

          H_O4(K)= (SIG_H_O4(6) * FINT(6) +
     $                   SIG_H_O4(7) * FINT(7) ) * RELO2**2 *DENS

C--------------------
C BrCl_noT
C--------------------
          SIG_BrCl_noT(0)=P2(CS_BrCl_noT(I0,1,1)*DLV2 +
     $                    CS_BrCl_noT(I0,1,2),
     $                       CS_BrCl_noT(I0,2,1)*DLV2 +
     $                    CS_BrCl_noT(I0,2,2),
     $                       CS_BrCl_noT(I0,3,1)*DLV2 +
     $                    CS_BrCl_noT(I0,3,2),V3_DU1)
          SIG_BrCl_noT(1)=P1(B1_BrCl_noT(I1,1),A1_BrCl_noT(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_BrCl_noT(I1,2),A1_BrCl_noT(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_BrCl_noT(2)=P1(B2_BrCl_noT(I2), A2_BrCl_noT(I2),
     $         V3_DU1)
          SIG_BrCl_noT(3)=P1(B3_BrCl_noT(I3), A3_BrCl_noT(I3),
     $         V3_DU2)
          SIG_BrCl_noT(4)=P2(C4_BrCl_noT(1),C4_BrCl_noT(2),C4_BrCl_noT
     $         (3),V3_DU2)
          SIG_BrCl_noT(5)=P2(C5_BrCl_noT(1),C5_BrCl_noT(2),C5_BrCl_noT
     $         (3),V3_DU2)
          SIG_BrCl_noT(6)=C6_BrCl_noT(1)
          SIG_BrCl_noT(7)=P1(C7_BrCl_noT(1),C7_BrCl_noT(2),V3_DU2)

          RJ_BrCl_noT(K)=   SIG_BrCl_noT(0)    * FINT(0) +
     $                     SIG_BrCl_noT(1)    * FINT(1) +
     $                     SIG_BrCl_noT(2)    * FINT(2) +
     $                     SIG_BrCl_noT(3)    * FINT(3) +
     $                     SIG_BrCl_noT(4)    * FINT(4) +
     $                     SIG_BrCl_noT(5)    * FINT(5) +
     $                     SIG_BrCl_noT(6)    * FINT(6) +
     $                     SIG_BrCl_noT(7)    * FINT(7)
C--------------------
C  ClNO2
C--------------------
          SIG_ClNO2(0)=P2(CS_ClNO2(I0,1,1)*DLV2 +
     $                    CS_ClNO2(I0,1,2),
     $                    CS_ClNO2(I0,2,1)*DLV2 +
     $                    CS_ClNO2(I0,2,2),
     $                    CS_ClNO2(I0,3,1)*DLV2 +
     $                    CS_ClNO2(I0,3,2),V3_DU1)
          SIG_ClNO2(1)=P1(B1_ClNO2(I1,1),A1_ClNO2(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_ClNO2(I1,2),A1_ClNO2(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_ClNO2(2)=P1(B2_ClNO2(I2), A2_ClNO2(I2),V3_DU1)
          SIG_ClNO2(3)=P1(B3_ClNO2(I3), A3_ClNO2(I3),V3_DU2)
          SIG_ClNO2(4)=P1(C4_ClNO2(1),C4_ClNO2(2),V3_DU2)
          SIG_ClNO2(5)=P1(C5_ClNO2(1),C5_ClNO2(2),V3_DU2)
          SIG_ClNO2(6)=C6_ClNO2(1)

          RJ_ClNO2(K)=   SIG_ClNO2(0)    * FINT(0) +
     $                     SIG_ClNO2(1)    * FINT(1) +
     $                     SIG_ClNO2(2)    * FINT(2) +
     $                     SIG_ClNO2(3)    * FINT(3) +
     $                     SIG_ClNO2(4)    * FINT(4) +
     $                     SIG_ClNO2(5)    * FINT(5) +
     $                     SIG_ClNO2(6)    * FINT(6)


C--------------------
C  BrNO2
C--------------------
          SIG_BrNO2(0)=P2(CS_BrNO2(I0,1,1)*DLV2 +
     $                    CS_BrNO2(I0,1,2),
     $                    CS_BrNO2(I0,2,1)*DLV2 +
     $                    CS_BrNO2(I0,2,2),
     $                    CS_BrNO2(I0,3,1)*DLV2 +
     $                    CS_BrNO2(I0,3,2),V3_DU1)
          SIG_BrNO2(1)=P1(B1_BrNO2(I1,1),A1_BrNO2(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_BrNO2(I1,2),A1_BrNO2(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_BrNO2(2)=P1(B2_BrNO2(I2), A2_BrNO2(I2),V3_DU1)
          SIG_BrNO2(3)=P1(B3_BrNO2(I3), A3_BrNO2(I3),V3_DU2)
          SIG_BrNO2(4)=P2(C4_BrNO2(1),C4_BrNO2(2),C4_BrNO2(3),V3_DU2)
          SIG_BrNO2(5)=P1(C5_BrNO2(1),C5_BrNO2(2),V3_DU2)
          SIG_BrNO2(6)=C6_BrNO2(1)
          SIG_BrNO2(7)=P1(C7_BrNO2(1),C7_BrNO2(2),V3_DU2)

          RJ_BrNO2(K)=   SIG_BrNO2(0)    * FINT(0) +
     $                     SIG_BrNO2(1)    * FINT(1) +
     $                     SIG_BrNO2(2)    * FINT(2) +
     $                     SIG_BrNO2(3)    * FINT(3) +
     $                     SIG_BrNO2(4)    * FINT(4) +
     $                     SIG_BrNO2(5)    * FINT(5) +
     $                     SIG_BrNO2(6)    * FINT(6) +
     $                     SIG_BrNO2(7)    * FINT(7)

C--------------------
C  Br2
C--------------------
          SIG_Br2(0)=P1(CS_Br2(I0,1,1)*DLV2 +
     $                    CS_Br2(I0,1,2),
     $                    CS_Br2(I0,2,1)*DLV2 +
     $                    CS_Br2(I0,2,2),V3_DU1)
          SIG_Br2(1)=P1(B1_Br2(I1,1),A1_Br2(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_Br2(I1,2),A1_Br2(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_Br2(2)=P1(B2_Br2(I2), A2_Br2(I2),V3_DU1)
          SIG_Br2(3)=P1(B3_Br2(I3), A3_Br2(I3),V3_DU2)
          SIG_Br2(4)=P2(C4_Br2(1),C4_Br2(2),C4_Br2(3),V3_DU2)
          SIG_Br2(5)=P2(C5_Br2(1),C5_Br2(2),C5_Br2(3),V3_DU2)
          SIG_Br2(6)=C6_Br2(1)
          SIG_Br2(7)=P1(C7_Br2(1),C7_Br2(2),V3_DU2)

          RJ_Br2(K)=   SIG_Br2(0)    * FINT(0) +
     $                     SIG_Br2(1)    * FINT(1) +
     $                     SIG_Br2(2)    * FINT(2) +
     $                     SIG_Br2(3)    * FINT(3) +
     $                     SIG_Br2(4)    * FINT(4) +
     $                     SIG_Br2(5)    * FINT(5) +
     $                     SIG_Br2(6)    * FINT(6) +
     $                     SIG_Br2(7)    * FINT(7)

C--------------------
C  IO
C--------------------
          SIG_IO(6)=C6_IO(1)
          SIG_IO(7)=P1(C7_IO(1),C7_IO(2),V3_DU2)

          RJ_IO(K)=      SIG_IO(6)    * FINT(6) +
     $                     SIG_IO(7)    * FINT(7)

C--------------------
C  INO3
C--------------------
          SIG_INO3(2)=P1(B2_INO3(I2), A2_INO3(I2),V3_DU1)
          SIG_INO3(3)=P1(B3_INO3(I3), A3_INO3(I3),V3_DU2)
          SIG_INO3(4)=P1(C4_INO3(1),C4_INO3(2),V3_DU2)
          SIG_INO3(5)=P1(C5_INO3(1),C5_INO3(2),V3_DU2)
          SIG_INO3(6)=C6_INO3(1)
          SIG_INO3(7)=P1(C7_INO3(1),C7_INO3(2),V3_DU2)

          RJ_INO3(K)=    SIG_INO3(2)    * FINT(2) +
     $                     SIG_INO3(3)    * FINT(3) +
     $                     SIG_INO3(4)    * FINT(4) +
     $                     SIG_INO3(5)    * FINT(5) +
     $                     SIG_INO3(6)    * FINT(6) +
     $                     SIG_INO3(7)    * FINT(7)

C--------------------
C  CH3I
C--------------------
          SIG_CH3I(0)=P1(CS_CH3I(I0,1,1)*DLV2 +
     $                    CS_CH3I(I0,1,2),
     $                    CS_CH3I(I0,2,1)*DLV2 +
     $                    CS_CH3I(I0,2,2),V3_DU1)
          SIG_CH3I(1)=P1(B1_CH3I(I1,1),A1_CH3I(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_CH3I(I1,2),A1_CH3I(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_CH3I(2)=P1(B2_CH3I(I2), A2_CH3I(I2),V3_DU1)
          SIG_CH3I(3)=P1(B3_CH3I(I3), A3_CH3I(I3),V3_DU2)
          SIG_CH3I(4)=P2(C4_CH3I(1),C4_CH3I(2),C4_CH3I(3),V3_DU2)
          SIG_CH3I(5)=P2(C5_CH3I(1),C5_CH3I(2),C5_CH3I(3),V3_DU2)
          SIG_CH3I(6)=C6_CH3I(1)

          RJ_CH3I(K)=   SIG_CH3I(0)    * FINT(0) +
     $                     SIG_CH3I(1)    * FINT(1) +
     $                     SIG_CH3I(2)    * FINT(2) +
     $                     SIG_CH3I(3)    * FINT(3) +
     $                     SIG_CH3I(4)    * FINT(4) +
     $                     SIG_CH3I(5)    * FINT(5) +
     $                     SIG_CH3I(6)    * FINT(6)


C--------------------
C  I2
C--------------------

          SIG_I2(6)=C6_I2(1)
          SIG_I2(7)=P1(C7_I2(1),C7_I2(2),V3_DU2)

          RJ_I2(K)=      SIG_I2(6)    * FINT(6) +
     $                     SIG_I2(7)    * FINT(7)

C--------------------
C  ICl
C--------------------
          SIG_ICl(1)=P1(B1_ICl(I1,1),A1_ICl(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_ICl(I1,2),A1_ICl(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_ICl(2)=P1(B2_ICl(I2), A2_ICl(I2),V3_DU1)
          SIG_ICl(3)=P1(B3_ICl(I3), A3_ICl(I3),V3_DU2)
          SIG_ICl(4)=P1(C4_ICl(1),C4_ICl(2),V3_DU2)
          SIG_ICl(5)=P1(C5_ICl(1),C5_ICl(2),V3_DU2)
          SIG_ICl(6)=C6_ICl(1)
          SIG_ICl(7)=P1(C7_ICl(1),C7_ICl(2),V3_DU2)

          RJ_ICl(K)=     SIG_ICl(1)    * FINT(1) +
     $                     SIG_ICl(2)    * FINT(2) +
     $                     SIG_ICl(3)    * FINT(3) +
     $                     SIG_ICl(4)    * FINT(4) +
     $                     SIG_ICl(5)    * FINT(5) +
     $                     SIG_ICl(6)    * FINT(6) +
     $                     SIG_ICl(7)    * FINT(7)

C--------------------
C  IBr
C--------------------

          SIG_IBr(1)=P1(B1_IBr(I1,1),A1_IBr(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_IBr(I1,2),A1_IBr(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_IBr(2)=P1(B2_IBr(I2), A2_IBr(I2),V3_DU1)
          SIG_IBr(3)=P1(B3_IBr(I3), A3_IBr(I3),V3_DU2)
          SIG_IBr(4)=P2(C4_IBr(1),C4_IBr(2),C4_IBr(3),V3_DU2)
          SIG_IBr(5)=P2(C5_IBr(1),C5_IBr(2),C5_IBr(3),V3_DU2)
          SIG_IBr(6)=C6_IBr(1)
          SIG_IBr(7)=P1(C7_IBr(1),C7_IBr(2),V3_DU2)

          RJ_IBr(K)=     SIG_IBr(1)    * FINT(1) +
     $                     SIG_IBr(2)    * FINT(2) +
     $                     SIG_IBr(3)    * FINT(3) +
     $                     SIG_IBr(4)    * FINT(4) +
     $                     SIG_IBr(5)    * FINT(5) +
     $                     SIG_IBr(6)    * FINT(6) +
     $                     SIG_IBr(7)    * FINT(7)

C--------------------
C  C3H7I
C--------------------

          SIG_C3H7I(1)=P1(B1_C3H7I(I1,1),A1_C3H7I(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_C3H7I(I1,2),A1_C3H7I(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_C3H7I(2)=P1(B2_C3H7I(I2), A2_C3H7I(I2),V3_DU1)
          SIG_C3H7I(3)=P1(B3_C3H7I(I3), A3_C3H7I(I3),V3_DU2)
          SIG_C3H7I(4)=P2(C4_C3H7I(1),C4_C3H7I(2),C4_C3H7I(3),V3_DU2)
          SIG_C3H7I(5)=P2(C5_C3H7I(1),C5_C3H7I(2),C5_C3H7I(3),V3_DU2)
          SIG_C3H7I(6)=C6_C3H7I(1)

          RJ_C3H7I(K)=   SIG_C3H7I(1)    * FINT(1) +
     $                     SIG_C3H7I(2)    * FINT(2) +
     $                     SIG_C3H7I(3)    * FINT(3) +
     $                     SIG_C3H7I(4)    * FINT(4) +
     $                     SIG_C3H7I(5)    * FINT(5) +
     $                     SIG_C3H7I(6)    * FINT(6)

C--------------------
C  CH2ClI
C--------------------

          SIG_CH2ClI(1)=P1(B1_CH2ClI(I1,1),A1_CH2ClI(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_CH2ClI(I1,2),A1_CH2ClI(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_CH2ClI(2)=P1(B2_CH2ClI(I2), A2_CH2ClI(I2),V3_DU1)
          SIG_CH2ClI(3)=P1(B3_CH2ClI(I3), A3_CH2ClI(I3),V3_DU2)
          SIG_CH2ClI(4)=P2(C4_CH2ClI(1),C4_CH2ClI(2),C4_CH2ClI(3)
     $         ,V3_DU2)
          SIG_CH2ClI(5)=P2(C5_CH2ClI(1),C5_CH2ClI(2),C5_CH2ClI(3)
     $         ,V3_DU2)
          SIG_CH2ClI(6)=C6_CH2ClI(1)


          RJ_CH2ClI(K)=  SIG_CH2ClI(1)    * FINT(1) +
     $                     SIG_CH2ClI(2)    * FINT(2) +
     $                     SIG_CH2ClI(3)    * FINT(3) +
     $                     SIG_CH2ClI(4)    * FINT(4) +
     $                     SIG_CH2ClI(5)    * FINT(5) +
     $                     SIG_CH2ClI(6)    * FINT(6)


C--------------------
C  CH2I2
C--------------------

          SIG_CH2I2(1)=P1(B1_CH2I2(I1,1),A1_CH2I2(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_CH2I2(I1,2),A1_CH2I2(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_CH2I2(2)=P1(B2_CH2I2(I2), A2_CH2I2(I2),V3_DU1)
          SIG_CH2I2(3)=P1(B3_CH2I2(I3), A3_CH2I2(I3),V3_DU2)
          SIG_CH2I2(4)=P1(C4_CH2I2(1),C4_CH2I2(2),V3_DU2)
          SIG_CH2I2(5)=P1(C5_CH2I2(1),C5_CH2I2(2),V3_DU2)
          SIG_CH2I2(6)=C6_CH2I2(1)


          RJ_CH2I2(K)=   SIG_CH2I2(1)    * FINT(1) +
     $                     SIG_CH2I2(2)    * FINT(2) +
     $                     SIG_CH2I2(3)    * FINT(3) +
     $                     SIG_CH2I2(4)    * FINT(4) +
     $                     SIG_CH2I2(5)    * FINT(5) +
     $                     SIG_CH2I2(6)    * FINT(6)


C--------------------
C  INO2
C--------------------

          SIG_INO2(1)=P1(B1_INO2(I1,1),A1_INO2(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_INO2(I1,2),A1_INO2(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_INO2(2)=P1(B2_INO2(I2), A2_INO2(I2),V3_DU1)
          SIG_INO2(3)=P1(B3_INO2(I3), A3_INO2(I3),V3_DU2)
          SIG_INO2(4)=P2(C4_INO2(1),C4_INO2(2),C4_INO2(3),V3_DU2)
          SIG_INO2(5)=P1(C5_INO2(1),C5_INO2(2),V3_DU2)
          SIG_INO2(6)=C6_INO2(1)


          RJ_INO2(K)=    SIG_INO2(1)    * FINT(1) +
     $                     SIG_INO2(2)    * FINT(2) +
     $                     SIG_INO2(3)    * FINT(3) +
     $                     SIG_INO2(4)    * FINT(4) +
     $                     SIG_INO2(5)    * FINT(5) +
     $                     SIG_INO2(6)    * FINT(6)


C--------------------
C  BrO_noT
C--------------------

          SIG_BrO_noT(3)=P1(B3_BrO_noT(I3), A3_BrO_noT(I3),
     $         V3_DU2)
          SIG_BrO_noT(4)=P2(C4_BrO_noT(1),C4_BrO_noT(2),C4_BrO_noT(3),
     $         V3_DU2)
          SIG_BrO_noT(5)=P1(C5_BrO_noT(1),C5_BrO_noT(2),V3_DU2)
          SIG_BrO_noT(6)=C6_BrO_noT(1)


          RJ_BrO_noT(K)= SIG_BrO_noT(3)    * FINT(3) +
     $                     SIG_BrO_noT(4)    * FINT(4) +
     $                     SIG_BrO_noT(5)    * FINT(5) +
     $                     SIG_BrO_noT(6)    * FINT(6)


C--------------------
C  OClO_noT
C--------------------

          SIG_OClO_noT(2)=P1(B2_OClO_noT(I2), A2_OClO_noT(I2),
     $         V3_DU1)
          SIG_OClO_noT(3)=P1(B3_OClO_noT(I3), A3_OClO_noT(I3),
     $         V3_DU2)
          SIG_OClO_noT(4)=P2(C4_OClO_noT(1),C4_OClO_noT(2),C4_OClO_noT
     $         (3),V3_DU2)
          SIG_OClO_noT(5)=P2(C5_OClO_noT(1),C5_OClO_noT(2),C5_OClO_noT
     $         (3),V3_DU2)
          SIG_OClO_noT(6)=C6_OClO_noT(1)
          SIG_OClO_noT(7)=P1(C7_OClO_noT(1),C7_OClO_noT(2),V3_DU2)

          RJ_OClO_noT(K)= SIG_OClO_noT(2)    * FINT(2) +
     $                     SIG_OClO_noT(3)    * FINT(3) +
     $                     SIG_OClO_noT(4)    * FINT(4) +
     $                     SIG_OClO_noT(5)    * FINT(5) +
     $                     SIG_OClO_noT(6)    * FINT(6) +
     $                     SIG_OClO_noT(7)    * FINT(7)


C--------------------
C  Cl2_noT
C--------------------

          SIG_Cl2_noT(1)=P1(B1_Cl2_noT(I1,1),A1_Cl2_noT(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_Cl2_noT(I1,2),A1_Cl2_noT(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_Cl2_noT(2)=P1(B2_Cl2_noT(I2), A2_Cl2_noT(I2),
     $         V3_DU1)
          SIG_Cl2_noT(3)=P1(B3_Cl2_noT(I3), A3_Cl2_noT(I3),
     $         V3_DU2)
          SIG_Cl2_noT(4)=P2(C4_Cl2_noT(1),C4_Cl2_noT(2),C4_Cl2_noT(3),
     $         V3_DU2)
          SIG_Cl2_noT(5)=P1(C5_Cl2_noT(1),C5_Cl2_noT(2),V3_DU2)
          SIG_Cl2_noT(6)=C6_Cl2_noT(1)
          SIG_Cl2_noT(7)=P1(C7_Cl2_noT(1),C7_Cl2_noT(2),V3_DU2)

          RJ_Cl2_noT(K)= SIG_Cl2_noT(1)    * FINT(1) +
     $                     SIG_Cl2_noT(2)    * FINT(2) +
     $                     SIG_Cl2_noT(3)    * FINT(3) +
     $                     SIG_Cl2_noT(4)    * FINT(4) +
     $                     SIG_Cl2_noT(5)    * FINT(5) +
     $                     SIG_Cl2_noT(6)    * FINT(6) +
     $                     SIG_Cl2_noT(7)    * FINT(7)

C--------------------
C  HOI_jen91
C--------------------

          SIG_HOI_jen91(3)=P1(B3_HOI_jen91(I3), A3_HOI_jen91(I3),
     $         V3_DU2)
          SIG_HOI_jen91(4)=P2(C4_HOI_jen91(1),C4_HOI_jen91(2),
     $         C4_HOI_jen91(3),V3_DU2)
          SIG_HOI_jen91(5)=P2(C5_HOI_jen91(1),C5_HOI_jen91(2),
     $         C5_HOI_jen91(3),V3_DU2)
          SIG_HOI_jen91(6)=C6_HOI_jen91(1)
          SIG_HOI_jen91(7)=P1(C7_HOI_jen91(1),C7_HOI_jen91(2),V3_DU2)

          RJ_HOI_jen91(K)=   SIG_HOI_jen91(3)    * FINT(3) +
     $                     SIG_HOI_jen91(4)    * FINT(4) +
     $                     SIG_HOI_jen91(5)    * FINT(5) +
     $                     SIG_HOI_jen91(6)    * FINT(6) +
     $                     SIG_HOI_jen91(7)    * FINT(7)

C--------------------
C   HOBr
C--------------------
          SIG_HOBr(2)=P1(B2_HOBr(I2), A2_HOBr(I2),V3_DU1)
          SIG_HOBr(3)=P1(B3_HOBr(I3), A3_HOBr(I3),V3_DU2)
          SIG_HOBr(4)=P2(C4_HOBr(1),C4_HOBr(2),C4_HOBr(3),V3_DU2)
          SIG_HOBr(5)=P1(C5_HOBr(1),C5_HOBr(2),V3_DU2)
          SIG_HOBr(6)=C6_HOBr(1)
          SIG_HOBr(7)=P1(C7_HOBr(1),C7_HOBr(2),V3_DU2)

          RJ_HOBr(K)=  SIG_HOBr(2)    * FINT(2) +
     $                     SIG_HOBr(3)    * FINT(3) +
     $                     SIG_HOBr(4)    * FINT(4) +
     $                     SIG_HOBr(5)    * FINT(5) +
     $                     SIG_HOBr(6)    * FINT(6) +
     $                     SIG_HOBr(7)    * FINT(7)

C--------------------
C  HONO
C--------------------
          SIG_HONO(4)=P2(C4_HONO(1),C4_HONO(2),C4_HONO(3),V3_DU2)
          SIG_HONO(5)=P2(C5_HONO(1),C5_HONO(2),C5_HONO(3),V3_DU2)
          SIG_HONO(6)=C6_HONO(1)

          RJ_HONO(K)=
     $                     SIG_HONO(4)    * FINT(4) +
     $                     SIG_HONO(5)    * FINT(5) +
     $                     SIG_HONO(6)    * FINT(6)


C--------------------
C  NO2m
C--------------------
          SIG_NO2m(3)=P1(B3_NO2m(I3), A3_NO2m(I3),V3_DU2)
          SIG_NO2m(4)=P1(C4_NO2m(1),C4_NO2m(2),V3_DU2)
          SIG_NO2m(5)=P1(C5_NO2m(1),C5_NO2m(2),V3_DU2)
          SIG_NO2m(6)=C6_NO2m(1)

          RJ_NO2m(K)=
     $                     SIG_NO2m(3)    * FINT(3) +
     $                     SIG_NO2m(4)    * FINT(4) +
     $                     SIG_NO2m(5)    * FINT(5) +
     $                     SIG_NO2m(6)    * FINT(6)

C--------------------
C NO3n (NO3-)
C--------------------
          SIG_NO3n(0)=P2(CS_NO3n(I0,1,1)*DLV2+ CS_NO3n(I0,1,2),
     $                   CS_NO3n(I0,2,1)*DLV2+ CS_NO3n(I0,2,2),
     $                   CS_NO3n(I0,3,1)*DLV2+ CS_NO3n(I0,3,2),
     $                   V3_DU1)
          SIG_NO3n(1)=P1(B1_NO3n(I1,1),A1_NO3n(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_NO3n(I1,2),A1_NO3n(I1,2),
     $                    V3_DU1) *V2S_m+
     $                 P1(B1_NO3n(I1,3),A1_NO3n(I1,3),
     $                    V3_DU1) *V2S_m**2
          SIG_NO3n(2)=P1(B2_NO3n(I2), A2_NO3n(I2),V3_DU1)
          SIG_NO3n(3)=P1(B3_NO3n(I3), A3_NO3n(I3),V3_DU2)
          SIG_NO3n(4)=P3(C4_NO3n(1),C4_NO3n(2),C4_NO3n(3),C4_NO3n(4),
     $                V3_DU2)
          SIG_NO3n(5)=P3(C5_NO3n(1),C5_NO3n(2),C5_NO3n(3),C5_NO3n(4),
     $                V3_DU2)
          SIG_NO3n(6)= P1(C6_NO3n(1),C6_NO3n(2),V3_DU2)


          RJ_NO3n(K)=  ( SIG_NO3n(0)    * FINT(0) +
     $                     SIG_NO3n(1)    * FINT(1) +
     $                     SIG_NO3n(2)    * FINT(2) +
     $                     SIG_NO3n(3)    * FINT(3) +
     $                     SIG_NO3n(4)    * FINT(4) +
     $                     SIG_NO3n(5)    * FINT(5) +
     $                     SIG_NO3n(6)    * FINT(6) )*
     &                     QYNO3n(K)


C--------------------
C 2nd version for NO3-
C--------------------
          SIG_NO3n(0)=P2(CS_NO3n(I0,1,1)*DLV2+ CS_NO3n(I0,1,2),
     $                   CS_NO3n(I0,2,1)*DLV2+ CS_NO3n(I0,2,2),
     $                   CS_NO3n(I0,3,1)*DLV2+ CS_NO3n(I0,3,2),
     $                   V3_DU1) *
     $               P2(1.D0,
     $                  FS_NO3n(I0,1,1)*DLV2+FS_NO3n(I0,1,2),
     $                  FS_NO3n(I0,2,1)*DLV2+FS_NO3n(I0,2,2),
     $                  TEMP1(K))
         SIG_NO3n(1)=P2(TJ_NO3n(1,1),TJ_NO3n(2,1),TJ_NO3n(3,1),
     $                   TEMP2(K))*
     $                (P1(B1_NO3n(I1,1),A1_NO3n(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_NO3n(I1,2),A1_NO3n(I1,2),
     $                    V3_DU1) *V2S_m+
     $                 P1(B1_NO3n(I1,3),A1_NO3n(I1,3),
     $                    V3_DU1) *V2S_m**2)
          SIG_NO3n(2)=P2(TJ_NO3n(1,2),TJ_NO3n(2,2),TJ_NO3n(3,2),
     $                   TEMP2(K))*
     $                P1(B2_NO3n(I2), A2_NO3n(I2),V3_DU1)
          SIG_NO3n(3)=P2(TJ_NO3n(1,3),TJ_NO3n(2,3),TJ_NO3n(3,3),
     $                   TEMP2(K))*
     $                P1(B3_NO3n(I3), A3_NO3n(I3),V3_DU2)
          SIG_NO3n(4)=P2(TJ_NO3n(1,4),TJ_NO3n(2,4),TJ_NO3n(3,4),
     $                   TEMP2(K))*
     $                P3(C4_NO3n(1),C4_NO3n(2),C4_NO3n(3),C4_NO3n(4),
     $                V3_DU2)
          SIG_NO3n(5)=P2(TJ_NO3n(1,5),TJ_NO3n(2,5),TJ_NO3n(3,5),
     $                   TEMP2(K))*
     $                P3(C5_NO3n(1),C5_NO3n(2),C5_NO3n(3),C5_NO3n(4),
     $                V3_DU2)
          SIG_NO3n(6)=P2(TJ_NO3n(1,6),TJ_NO3n(2,6),TJ_NO3n(3,6),
     $                   TEMP2(K))*
     $                P1(C6_NO3n(1),C6_NO3n(2),V3_DU2)

c          RJ_NO3n(K)=   SIG_NO3n(0)    * FINT(0) +
          RJ_dumm23(K)=   SIG_NO3n(0)    * FINT(0) +
     $                     SIG_NO3n(1)    * FINT(1) +
     $                     SIG_NO3n(2)    * FINT(2) +
     $                     SIG_NO3n(3)    * FINT(3) +
     $                     SIG_NO3n(4)    * FINT(4) +
     $                     SIG_NO3n(5)    * FINT(5) +
     $                     SIG_NO3n(6)    * FINT(6)



C--------------------
C  ---- testeinbau: dumm24
C--------------------
          SIG_dumm24(0)=P1(CS_dumm24(I0,1,1)*DLV2 +
     $                    CS_dumm24(I0,1,2),
     $                    CS_dumm24(I0,2,1)*DLV2 +
     $                    CS_dumm24(I0,2,2),V3_DU1)
          SIG_dumm24(1)=P1(B1_dumm24(I1,1),A1_dumm24(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_dumm24(I1,2),A1_dumm24(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_dumm24(2)=P1(B2_dumm24(I2), A2_dumm24(I2),V3_DU1)
          SIG_dumm24(3)=P1(B3_dumm24(I3), A3_dumm24(I3),V3_DU2)
          SIG_dumm24(4)=P1(C4_dumm24(1),C4_dumm24(2),V3_DU2)
          SIG_dumm24(5)=P1(C5_dumm24(1),C5_dumm24(2),V3_DU2)
          SIG_dumm24(6)=C6_dumm24(1)
          SIG_dumm24(7)=P1(C7_dumm24(1),C7_dumm24(2),V3_DU2)

          RJ_dumm24(K)=   SIG_dumm24(0)    * FINT(0) +
     $                     SIG_dumm24(1)    * FINT(1) +
     $                     SIG_dumm24(2)    * FINT(2) +
     $                     SIG_dumm24(3)    * FINT(3) +
     $                     SIG_dumm24(4)    * FINT(4) +
     $                     SIG_dumm24(5)    * FINT(5) +
     $                     SIG_dumm24(6)    * FINT(6) +
     $                     SIG_dumm24(7)    * FINT(7)

C--------------------
C  ---- testeinbau: dumm25
C--------------------
          SIG_dumm25(0)=P1(CS_dumm25(I0,1,1)*DLV2 +
     $                    CS_dumm25(I0,1,2),
     $                    CS_dumm25(I0,2,1)*DLV2 +
     $                    CS_dumm25(I0,2,2),V3_DU1)
          SIG_dumm25(1)=P1(B1_dumm25(I1,1),A1_dumm25(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_dumm25(I1,2),A1_dumm25(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_dumm25(2)=P1(B2_dumm25(I2), A2_dumm25(I2),V3_DU1)
          SIG_dumm25(3)=P1(B3_dumm25(I3), A3_dumm25(I3),V3_DU2)
          SIG_dumm25(4)=P1(C4_dumm25(1),C4_dumm25(2),V3_DU2)
          SIG_dumm25(5)=P1(C5_dumm25(1),C5_dumm25(2),V3_DU2)
          SIG_dumm25(6)=C6_dumm25(1)
          SIG_dumm25(7)=P1(C7_dumm25(1),C7_dumm25(2),V3_DU2)

          RJ_dumm25(K)=   SIG_dumm25(0)    * FINT(0) +
     $                     SIG_dumm25(1)    * FINT(1) +
     $                     SIG_dumm25(2)    * FINT(2) +
     $                     SIG_dumm25(3)    * FINT(3) +
     $                     SIG_dumm25(4)    * FINT(4) +
     $                     SIG_dumm25(5)    * FINT(5) +
     $                     SIG_dumm25(6)    * FINT(6) +
     $                     SIG_dumm25(7)    * FINT(7)

C--------------------
C  ---- testeinbau: dumm26
C--------------------
          SIG_dumm26(0)=P1(CS_dumm26(I0,1,1)*DLV2 +
     $                    CS_dumm26(I0,1,2),
     $                    CS_dumm26(I0,2,1)*DLV2 +
     $                    CS_dumm26(I0,2,2),V3_DU1)
          SIG_dumm26(1)=P1(B1_dumm26(I1,1),A1_dumm26(I1,1),
     $                    V3_DU1) +
     $                 P1(B1_dumm26(I1,2),A1_dumm26(I1,2),
     $                    V3_DU1) *V2S_m
          SIG_dumm26(2)=P1(B2_dumm26(I2), A2_dumm26(I2),V3_DU1)
          SIG_dumm26(3)=P1(B3_dumm26(I3), A3_dumm26(I3),V3_DU2)
          SIG_dumm26(4)=P1(C4_dumm26(1),C4_dumm26(2),V3_DU2)
          SIG_dumm26(5)=P1(C5_dumm26(1),C5_dumm26(2),V3_DU2)
          SIG_dumm26(6)=C6_dumm26(1)
          SIG_dumm26(7)=P1(C7_dumm26(1),C7_dumm26(2),V3_DU2)

          RJ_dumm26(K)=   SIG_dumm26(0)    * FINT(0) +
     $                     SIG_dumm26(1)    * FINT(1) +
     $                     SIG_dumm26(2)    * FINT(2) +
     $                     SIG_dumm26(3)    * FINT(3) +
     $                     SIG_dumm26(4)    * FINT(4) +
     $                     SIG_dumm26(5)    * FINT(5) +
     $                     SIG_dumm26(6)    * FINT(6) +
     $                     SIG_dumm26(7)    * FINT(7)



      ENDIF

1111  continue

      END SUBROUTINE PHOTO_CAL
********************************************************************
