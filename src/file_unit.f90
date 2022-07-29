module file_unit

! Description :
! -----------
!   Declare and define the unit of each I/O file used in Mistra


! Author :
! ------
!   Josue Bock


  ! JPFUNERR: unit number for error messages
  ! JPFUNOUT: unit number for standard output

  ! JPFUNNAM: unit number for namelist
  ! JPFUNCFGOUT: unit number for config output

  ! JPFUNAERRAD:  unit number for aerosol (particles) radiation parameters
  ! JPFUNDATARAD: unit number for radiative code data

  implicit none

  public
  save

  integer, parameter :: jpfunerr = 0
  integer, parameter :: jpfunout = 6

  integer, parameter :: jpfunnam = 1
  integer, parameter :: jpfuncfgout = 2

  integer, parameter :: jpfuneul = 3
  integer, parameter :: jpfunpi = 4       ! initial output of atmospheric params for plotting
  integer, parameter :: jpfungas = 7
  integer, parameter :: jpfunrad = 8
  integer, parameter :: jpfunint = 9
  integer, parameter :: jpfunfi = 10      ! initial output of aerosol params for plotting

! prof* files
  integer, parameter :: jpfunprofm = 11
  integer, parameter :: jpfunprofc = 12
  integer, parameter :: jpfunprofr = 13

! out_mass
  integer, parameter :: jpfunom = 14

! restart files
  integer, parameter :: jpfunrstm = 15
  integer, parameter :: jpfunrstc = 16

! ion balance files (aerosol/droplets)
  integer, parameter :: jpfuniba = 17
  integer, parameter :: jpfunibd = 18
! final output of aerosol size distribution
  integer, parameter :: jpfunae = 19

! ploutm files: pm*, pb*
  integer, parameter :: jpfunpm = 20
  integer, parameter :: jpfunpb = 21
! ploutp files: f1*, f2*, f3*
  integer, parameter :: jpfunf1 = 22
  integer, parameter :: jpfunf2 = 23
  integer, parameter :: jpfunf3 = 24
! ploutr file: pr*
  integer, parameter :: jpfunpr = 25
! ploutt file: pt*
  integer, parameter :: jpfunpt = 26
! ploutc files:
  integer, parameter :: jpfunsg1 = 27
  integer, parameter :: jpfunion = 28
  integer, parameter :: jpfunsl1 = 29
  integer, parameter :: jpfunsr1 = 30
  integer, parameter :: jpfungr = 31
  integer, parameter :: jpfungs = 32
! ploutj file:
  integer, parameter :: jpfunjra = 33

! nuc files:
  integer, parameter :: jpfunnuc0 = 34
  integer, parameter :: jpfunnuc1 = 35
  integer, parameter :: jpfunnuc2 = 36
  integer, parameter :: jpfunnuc3 = 37
  integer, parameter :: jpfunnuc4 = 38
  integer, parameter :: jpfunnuc5 = 39

  integer, parameter :: jpfunpph = 40
  integer, parameter :: jpfungam = 41
! mass balance -- deposition
  integer, parameter :: jpfunmass = 42

! jrate I/O files
  integer, parameter :: jpfunflx = 43
  integer, parameter :: jpfunsig = 44
  integer, parameter :: jpfuncheb = 45
  integer, parameter :: jpfunnlt = 46
  integer, parameter :: jpfunprofj = 47
  integer, parameter :: jpfunchk4 = 48
  integer, parameter :: jpfunf4st = 49

  integer, parameter :: jpfunclarke = 50

  integer, parameter :: jpfunaerrad  = 51 ! units used: 51-56
  integer, parameter :: jpfundatarad = 57

  integer, parameter :: jpfunJchamb = 58

end module file_unit
