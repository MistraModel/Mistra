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



! netCDF-output of mistra
c written by Astrid Kerkweg 16. Juli 2002
c modified by Roland von Glasow, April-Aug 2004
c supplemented by Susanne Marquart, Sep 2004
! jjb general modif (2015-2021):
!    - cleaning unused variables
!    - size idvar_jrat(50) corrected (was 49)
!    - implicit none everywhere
!    - headers and declaration blocks

      subroutine open_netcdf (n_bln,chem,mic,halo,iod,box,chamber,nuc)

      implicit none

      integer, intent(in) :: n_bln
      logical, intent(in) :: chem, mic, halo, iod, box, chamber, nuc

      logical :: true

! ================================================================
!                          15/12/2016
! jjb a plot_cst should be added
!    example, plus translation to f90 can be found in pafog
!
!    Note that in pafog, the type "plot_file" contains a logical,
!        plot_flag, which is maybe useless (hard coded)
!    In the same way, the files names should be read in a used file
!        (namelist)
!
!    Investigate what NF90_SYNC does exactly
!
! ================================================================

c open netCDF-files
      true=.true. ! to be able to compare hal vs. nohal
      call open_met (n_bln)                         ! thermodynamics
      if (mic.and..not.box)   call open_mic         ! microphysics
c      if (chem)  call open_chem_gas(n_bln,halo,iod) ! gas phase
c      if (chem)  call open_chem_aq(n_bln,halo,iod) ! aqueous phase
!     if (chem)  call open_chem_gas(n_bln,true,iod,nuc) ! gas phase ! jjb halo(=.true.), iod & nud are no longer used
      if (chem)  call open_chem_gas(n_bln)              ! gas phase ! jjb removed
      if (chem)  call open_chem_aq(n_bln,true,iod,nuc) ! aqueous phase
      if (chem)  call open_jrate (n_bln, chamber)   ! photolysis rates
      if (chem)  call open_rxn                      ! reaction rates
      if (nuc)   call open_nuc                      ! nucleation
      call open_grid ! writes information on grid that is not f(t)
c      call write_grid ! writes information on grid that is not f(t)
c      this call is now in main program AFTER call to SR grid

      end subroutine open_netcdf

c
c----------------------------------------------------------------
c

      subroutine write_netcdf (n_bln,chem,mic,halo,iod,box,chamber,nuc)

      implicit none

      integer, intent(in) :: n_bln
      logical, intent(in) :: chem, mic, halo, iod, box, chamber, nuc

      logical :: true

c write netCDF-files
      true=.true.
      call write_met (n_bln)                          ! thermodynamics

      if (mic.and..not.box) call write_mic            ! microphysics
c      if (chem)  call write_chem_aq (n_bln,halo,iod) ! aqueous phase

      if (chem)  call write_chem_gas (n_bln)             ! gas phase
      if (chem)  call write_chem_aq (n_bln,true,iod,nuc) ! aqueous phase
      if (chem)  call write_jrate (n_bln,chamber)     ! photolysis rates
      if (chem)  call write_rxn                       ! reaction rates
      if (nuc)   call write_nuc                       ! nucleation

      end subroutine write_netcdf
c
c----------------------------------------------------------------
c

      subroutine close_netcdf (mic,chem,box,nuc)

      implicit none

      logical, intent(in) :: chem, mic, box, nuc

c close netCDF-files
      call close_met
      if (mic.and..not.box)   call close_mic
      if (chem)  call close_chem_gas
      if (chem)  call close_chem_aq
      if (chem)  call close_jrate
      if (chem)  call close_rxn
      if (nuc)   call close_nuc

      end subroutine close_netcdf

c
c----------------------------------------------------------------
c

      subroutine open_met (n_bln)
c initialize plot file for meteorological output

      implicit none

! Include statements:
      include 'netcdf.inc'

      integer, intent(in) :: n_bln
! Local parameters:
      integer, parameter :: x=1,y=1,noz=1
      character (len=*), parameter :: fname = "meteo.nc"
      integer :: k, jddim1(4)
      integer :: id_n, id_noz, id_x, id_y
      common /cdf_var/ id_rec,idvar(39),idfile,icount,jddim(4)
      integer :: id_rec, idvar, idfile, icount, jddim

      icount=0

      k=nf_create(fname,nf_clobber,idfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,nf_global,'title',11,'meteorology')
      if (k.ne.nf_noerr) call ehandle(k,fname)
c dimension
      k=nf_def_dim(idfile,'n',n_bln,id_n)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idfile,'x',x,id_x)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idfile,'y',y,id_y)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idfile,'noz',noz,id_noz)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idfile,'rec',nf_unlimited,id_rec)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c variables
c time variables
      jddim(1)=id_x
      jddim(2)=id_y
      jddim(3)=id_noz
      jddim(4)=id_rec

      k=nf_def_var(idfile,'lday',nf_int,4,jddim,idvar(1))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(1),'long_name',3,'day')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idfile,'lst',nf_int,4,jddim,idvar(2))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(2),'long_name',4,'hour')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idfile,'lmin',nf_int,4,jddim,idvar(3))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(3),'long_name',6,'minute')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c wind field
      jddim1(1)=id_x
      jddim1(2)=id_y
      jddim1(3)=id_n
      jddim1(4)=id_rec
! jjb '4,jddim1' can be replaced by the necessary jddim1(i)
!    it might be necessary to create another 2D array (jddim2 for instance)
!    containing ONLY the necessary dimensions (rec & n, for instance)
      k=nf_def_var(idfile,'u',nf_float,4,jddim1,idvar(4))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(4),'long_name',22,
     & 'wind speed x-direction')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(4),'units',7,
     & 'm s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'v',nf_float,4,jddim1,idvar(5))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(5),'long_name',22,
     & 'wind speed y-direction')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(5),'units',7,
     & 'm s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'w',nf_float,4,jddim1,idvar(6))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(6),'long_name',22,
     & 'wind speed z-direction')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(6),'units',7,
     & 'm s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c temperature variables
      k=nf_def_var(idfile,'theta',nf_float,4,jddim1,idvar(7))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(7),'long_name',21,
     & 'potential temperature')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(7),'units',1,
     & 'K')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'thetl',nf_float,4,jddim1,idvar(8))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(8),'long_name',35,
     & 'liquid water potential temperature')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(8),'units',1,
     & 'K')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'temp',nf_float,4,jddim1,idvar(9))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(9),'long_name',11,
     & 'temperature')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(9),'units',1,
     & 'K')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c humidity variables
      k=nf_def_var(idfile,'xm1',nf_float,4,jddim1,idvar(10))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(10),'long_name',17,
     & 'specific humidity')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(10),'units',9,
     & 'kg kg-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'xm2',nf_float,4,jddim1,idvar(11))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(11),'long_name',20,
     & 'liquid water content')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(11),'units',9,
     & 'kg m-3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'feu',nf_float,4,jddim1,idvar(12))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(12),'long_name',17,
     & 'relative humidity')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(12),'units',1,
     & '1 ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c other variables
      k=nf_def_var(idfile,'p',nf_float,4,jddim1,idvar(13))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(13),'long_name',12,
     & 'air pressure')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(13),'units',3,
     & 'Pa')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'rho',nf_float,4,jddim1,idvar(14))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(14),'long_name',11,
     & 'air density')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(14),'units',8,
     & 'kg m-3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c these are necessary for easy plotting with Ferret:
      k=nf_def_var(idfile,'eta',nf_float,4,jddim1,idvar(15))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(15),'long_name',17,
     & 'height (cell mid)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(15),'units',3,'m')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'etw',nf_float,4,jddim1,idvar(16))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(16),'long_name',17,
     & 'height (cell lim)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(16),'units',3,'m')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'mtime',nf_float,4,jddim1,idvar(17))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(17),'long_name',10,'model time')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(17),'units',3,'h')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c heating rates
      k=nf_def_var(idfile,'dtrad',nf_float,4,jddim1,idvar(18))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(18),'long_name',22,
     &     'heating rate radiation')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(18),'units',5, 'K s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'dtcon',nf_float,4,jddim1,idvar(19))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(19),'long_name',25,
     &     'heating rate condensation')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(19),'units',5, 'K s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'dtdts',nf_float,4,jddim1,idvar(20))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(20),'long_name',15,
     &     'heating rate SW')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(20),'units',5,'K s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'dtdtl',nf_float,4,jddim1,idvar(21))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(21),'long_name',15,
     &     'heating rate LW')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(21),'units',5,'K s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c turbulence parameters
      k=nf_def_var(idfile,'atke',nf_float,4,jddim1,idvar(22))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(22),'long_name',17,
     &     'exch. coeff. TKE')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(22),'units',6, 'm2 s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'atkh',nf_float,4,jddim1,idvar(23))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(23),'long_name',17,
     &     'exch. coeff. heat')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(23),'units',6, 'm2 s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'atkm',nf_float,4,jddim1,idvar(24))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(24),'long_name',21,
     &     'exch. coeff. momentum')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(24),'units',6, 'm2 s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'tke',nf_float,4,jddim1,idvar(25))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(25),'long_name',24,
     &     'turbulent kinetic energy')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(25),'units',6, 'm2 s-2')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'tkep',nf_float,4,jddim1,idvar(26))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(26),'long_name',8,'TKE prod')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(26),'units',6, 'm2 s-3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'xl',nf_float,4,jddim1,idvar(27))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(27),'long_name',13,'mixing length')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(27),'units',1, 'm')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c flux divergences
      k=nf_def_var(idfile,'fd_u',nf_float,4,jddim1,idvar(28))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(28),'long_name',15,
     &     'd(K_m du/dz)/dz')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(28),'units',6,'m s-2')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'fd_v',nf_float,4,jddim1,idvar(29))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(29),'long_name',15,
     &     'd(K_m dv/dz)/dz')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(29),'units',6,'m s-2')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'fd_q',nf_float,4,jddim1,idvar(30))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(30),'long_name',15,
     &     'd(K_h dq/dz)/dz')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(30),'units',11,'kg kg-1 s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'fd_theta',nf_float,4,jddim1,idvar(31))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(31),'long_name',15,
     &     'd(K_h dT/dz)/dz')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(31),'units',5, 'K s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'fd_tke',nf_float,4,jddim1,idvar(32))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(32),'long_name',17,
     &     'd(K_h dTKE/dz)/dz')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(32),'units',7, 'm2 s-3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c particle numbers
      k=nf_def_var(idfile,'fsum_0',nf_float,4,jddim1,idvar(33))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(33),'long_name',11,
     &     'fsum: total')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(33),'units',4, 'cm-3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'fsum_1',nf_float,4,jddim1,idvar(34))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(34),'long_name',11,
     &     'fsum: bin 1')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(34),'units',4, 'cm-3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'fsum_2',nf_float,4,jddim1,idvar(35))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(35),'long_name',11,
     &     'fsum: bin 2')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(35),'units',4, 'cm-3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'fsum_3',nf_float,4,jddim1,idvar(36))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(36),'long_name',11,
     &     'fsum: bin 3')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(36),'units',4, 'cm-3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'fsum_4',nf_float,4,jddim1,idvar(37))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(37),'long_name',11,
     &     'fsum: bin 4')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(37),'units',4, 'cm-3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c lower and upper limit of cloud
      k=nf_def_var(idfile,'lcl',nf_int,4,jddim,idvar(38))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(38),'long_name',3,'lcl')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'lct',nf_int,4,jddim,idvar(39))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar(39),'long_name',3,'lct')
      if (k.ne.nf_noerr) call ehandle(k,fname)




c end define mode
      k=nf_enddef(idfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine open_met

c
c----------------------------------------------------------------
c


      subroutine open_grid
c initialize plot file for grid output
c ferret complains about variable ordering but the way I set it up "i", "j", "k" refer to the same parameters (nka, nkt, n) in grid and f etc

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nka,
     &     nkt

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = "grid.nc"
      integer, parameter :: x=1,y=1
      integer :: k, jddim1(4)
      integer :: id_n, id_nka, id_nkt, id_x, id_y
      common /cdf_var_grid/ id_rec,idvar_g(9),idfile,icount,jddim(4)
      integer :: id_rec, idvar_g, idfile, icount, jddim
      icount=0

      k=nf_create(fname,nf_clobber,idfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,nf_global,'title',4,'grid')
      if (k.ne.nf_noerr) call ehandle(k,fname)
c dimension
      k=nf_def_dim(idfile,'nka',nka,id_nka)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idfile,'nkt',nkt,id_nkt)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idfile,'n',n,id_n)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idfile,'x',x,id_x)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idfile,'y',y,id_y)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idfile,'rec',nf_unlimited,id_rec)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c variables
      jddim1(1)=id_x
      jddim1(2)=id_y
      jddim1(3)=id_n
      jddim1(4)=id_rec

c these are necessary for easy plotting with Ferret:
      k=nf_def_var(idfile,'eta',nf_float,4,jddim1,idvar_g(1))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(1),'long_name',17,
     & 'height (cell mid)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(1),'units',3,'m')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'etw',nf_float,4,jddim1,idvar_g(2))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(2),'long_name',17,
     & 'height (cell lim)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(2),'units',3,'m')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      jddim1(1)=id_nka
      jddim1(2)=id_x
      jddim1(3)=id_x

      k=nf_def_var(idfile,'enw',nf_float,4,jddim1,idvar_g(3))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(3),'long_name',28,
     & 'aerosol mass (middle of bin)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(3),'units',3,'mg')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'en',nf_float,4,jddim1,idvar_g(4))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(4),'long_name',26,
     & 'aerosol mass (max of bin)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(4),'units',3,'mg')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      jddim1(1)=id_x
      jddim1(2)=id_nkt
      jddim1(3)=id_x

      k=nf_def_var(idfile,'ew',nf_float,4,jddim1,idvar_g(5))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(5),'long_name',26,
     & 'water mass (middle of bin)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(5),'units',3,'mg')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'e',nf_float,4,jddim1,idvar_g(6))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(6),'long_name',24,
     & 'water mass (max of bin)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(6),'units',3,'mg')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      jddim1(1)=id_nka
      jddim1(2)=id_x
      jddim1(3)=id_x

      k=nf_def_var(idfile,'rn',nf_float,4,jddim1,idvar_g(7))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(7),'long_name',12,
     & 'radius (dry)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(7),'units',3,'um')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      jddim1(1)=id_nka
      jddim1(2)=id_nkt
      jddim1(3)=id_x

      k=nf_def_var(idfile,'rq',nf_float,4,jddim1,idvar_g(8))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(8),'long_name',30,
     & 'radius (total, middle of bin)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(8),'units',3,'um')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idfile,'rw',nf_float,4,jddim1,idvar_g(9))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(9),'long_name',27,
     & 'radius (total, max of bin)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idfile,idvar_g(9),'units',3,'um')
      if (k.ne.nf_noerr) call ehandle(k,fname)



c end define mode
      k=nf_enddef(idfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine open_grid


c
c----------------------------------------------------------------
c

      subroutine open_mic
c open netCDF file for microphysics


      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     nka,
     &     nkt

      USE precision, ONLY :
     &     dp

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = 'mic.nc'
      integer, parameter :: x=1, x2=2, y=1, noz=1
      !integer, parameter :: nat=nkt
      integer, parameter :: n1d=nka-1
! Local scalars:
      integer :: id_nf,id_n10
      integer :: id_nka,id_nkt  !,id_nat
      integer :: id_n1d,id_nka2
      integer :: id_noz,id_x,id_x2,id_y
      integer :: k
      integer :: n10
! Local arrays:
      integer :: jddim1(4)

! Common blocks:
      common /cdf_var_mic/ id_mic_rec,idvar_mic(6),idmicfile,
     & imiccount,jddim_mic(4)
      integer :: id_mic_rec, idvar_mic, idmicfile, imiccount, jddim_mic
      common /oneDsj/ rpw(nka), part1D(nka-1,nf)
      real (kind=dp) :: rpw, part1D
!- End of header ---------------------------------------------------------------

      imiccount=0
      n10=nf/10
      k=nf_create(fname,nf_clobber,idmicfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idmicfile,nf_global,'title',12,
     &  'microphysics')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idmicfile,nf_global,'title',16,
     & 'particle spectra')
      if (k.ne.nf_noerr) call ehandle(k,fname)
c dimensions
      k=nf_def_dim(idmicfile,'nka',nka,id_nka)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idmicfile,'nka2',nka,id_nka2)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idmicfile,'nkt',nkt,id_nkt)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      !k=nf_def_dim(idmicfile,'nat',nat,id_nat)
      !if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idmicfile,'n1d',n1d,id_n1d)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idmicfile,'n',n10,id_n10)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idmicfile,'nf',nf,id_nf)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idmicfile,'x',x,id_x)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idmicfile,'x2',x2,id_x2)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idmicfile,'y',y,id_y)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idmicfile,'noz',noz,id_noz)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idmicfile,'rec',nf_unlimited,id_mic_rec)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      jddim_mic(1)=id_x
      jddim_mic(2)=id_y
      jddim_mic(3)=id_noz
      jddim_mic(4)=id_mic_rec

c variables
c time variables
      k=nf_def_var(idmicfile,'lday',nf_int,4,jddim_mic,idvar_mic(1))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idmicfile,idvar_mic(1),'long_name',3,'day')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idmicfile,'lst',nf_int,4,jddim_mic,idvar_mic(2))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idmicfile,idvar_mic(2),'long_name',4,'hour')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idmicfile,'lmin',nf_int,4,jddim_mic,idvar_mic(3))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idmicfile,idvar_mic(3),'long_name',6,'minute')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      jddim1(1)=id_nka
      jddim1(2)=id_nkt
      jddim1(3)=id_n10
      jddim1(4)=id_mic_rec

      k=nf_def_var(idmicfile,'ff',nf_float,4,jddim1,idvar_mic(4))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idmicfile,idvar_mic(4),'long_name',17,
     & 'particle spectrum')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idmicfile,idvar_mic(4),'units',9,
     & 'part cm-3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c constant variables
      jddim1(1)=id_nka2 ! this must be passed as an array, even if dimension is one
      k=nf_def_var(idmicfile,'rpw',nf_float,1,jddim1(1:1),idvar_mic(5))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idmicfile,idvar_mic(5),'long_name',14,
     & '1D radius wall')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idmicfile,idvar_mic(5),'units',2,
     & 'um')
      if (k.ne.nf_noerr) call ehandle(k,fname)


      jddim1(1)=id_n1d
      jddim1(2)=id_nf
      jddim1(3)=id_mic_rec

      k=nf_def_var(idmicfile,'part1D',nf_float,3,jddim1(1:3),
     &     idvar_mic(6))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idmicfile,idvar_mic(6),'long_name',23,
     & '1D particle spectrum: N')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idmicfile,idvar_mic(6),'units',9,
     & 'part cm-3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

!      jddim1(1)=id_x

!      k=nf_def_var(idmicfile,'partr',nf_float,4,jddim1,idvar_mic(6))
!      if (k.ne.nf_noerr) call ehandle(k,fname)
!      k=nf_put_att_text(idmicfile,idvar_mic(6),'long_name',9,
!     & '1D radius')
!      if (k.ne.nf_noerr) call ehandle(k,fname)
!      k=nf_put_att_text(idmicfile,idvar_mic(6),'units',2,
!     & 'um')
!      if (k.ne.nf_noerr) call ehandle(k,fname)

c end define mode
      k=nf_enddef(idmicfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c write constant variable
      k=nf_put_var_double(idmicfile, idvar_mic(5),rpw)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine open_mic


c
c----------------------------------------------------------------
c

!     subroutine open_chem_gas(n_bln,halo,iod,nuc) ! jjb halo(=.true.), iod & nud are no longer used
      subroutine open_chem_gas (n_bln)              ! jjb removed
c open netCDF file for gas phase chemistry output

      USE gas_common, ONLY :
! Imported Parameters:
     &     j1,
     &     j5,
! Imported Array Variables with intent (in):
     &     gas_name,      rad_name,
     &     gas_name_long, rad_name_long

      USE cdf_var_gas, ONLY :
     &     idgasfile,
     &     igascount,
     &     idvar_gas

      implicit none

! Include statements:
      include 'netcdf.inc'

! Subroutine arguments
! Scalar arguments with intent(in):
      integer, intent(in) :: n_bln

! Local parameters:
      character (len=*), parameter :: fname = 'gas.nc'
      integer, parameter :: x=1, y=1, noz=1
      integer, parameter :: j0 = 3              ! number of time variables
! Local scalars:
      integer :: i0                             ! index offset
      integer :: id_n, id_noz, id_x, id_y
      integer :: idgas_rec
      integer :: j, k
! Local arrays:
      integer :: jddim_gas(4), jddim1(4)

!- End of header ---------------------------------------------------------------


      allocate ( idvar_gas(j0+j1+j5) )

      igascount=0

      k=nf_create(fname,nf_clobber,idgasfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idgasfile,nf_global,'title',9,'gas phase')
      if (k.ne.nf_noerr) call ehandle(k,fname)
c dimensions
      k=nf_def_dim(idgasfile,'n',n_bln,id_n)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idgasfile,'x',x,id_x)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idgasfile,'y',y,id_y)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idgasfile,'noz',noz,id_noz)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idgasfile,'rec',nf_unlimited,idgas_rec)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c variables
c time variables
      jddim_gas(1)=id_x
      jddim_gas(2)=id_y
      jddim_gas(3)=id_noz
      jddim_gas(4)=idgas_rec

      k=nf_def_var(idgasfile,'lday',nf_int,4,jddim_gas,idvar_gas(1))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idgasfile,idvar_gas(1),'long_name',3,'day')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idgasfile,'lst',nf_int,4,jddim_gas,idvar_gas(2))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idgasfile,idvar_gas(2),'long_name',4,'hour')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idgasfile,'lmin',nf_int,4,jddim_gas,idvar_gas(3))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idgasfile,idvar_gas(3),'long_name',6,'minute')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c chemical species

c long lived gas phase species

      jddim1(1)=id_x
      jddim1(2)=id_y
      jddim1(3)=id_n
      jddim1(4)=idgas_rec

c offset due to time variables
      i0 = j0

      do j = 1,j1
         k=nf_def_var(idgasfile,gas_name(j),nf_float,4,jddim1,
     &                idvar_gas(i0+j))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idgasfile,idvar_gas(i0+j),'long_name',
     &                len_trim(gas_name_long(j)),trim(gas_name_long(j)))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idgasfile,idvar_gas(i0+j),'units',
     &                     7,'mol m-3')
         if (k.ne.nf_noerr) call ehandle(k,fname)
      end do


c radicals
c offset: i0+j1
      i0 = i0 + j1

      do j = 1,j5
         k=nf_def_var(idgasfile,rad_name(j),nf_float,4,jddim1,
     &                idvar_gas(i0+j))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idgasfile,idvar_gas(i0+j),'long_name',
     &                len_trim(rad_name_long(j)),trim(rad_name_long(j)))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idgasfile,idvar_gas(i0+j),'units',
     &                     7,'mol m-3')
         if (k.ne.nf_noerr) call ehandle(k,fname)
      end do

c end define mode
      k=nf_enddef(idgasfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine open_chem_gas


c
c----------------------------------------------------------------
c

      subroutine open_chem_aq(n_bln,halo,iod,nuc)
c open netCDF file for aqueous phase chemistry output

      USE config, ONLY:
     &     nkc_l

      USE global_params, ONLY :
! Imported Parameters:
     &     j1=>j1_fake,
     &     j2,
     &     j3,
     &     j6

      implicit none

! Include statements:
      include 'netcdf.inc'

      logical, intent(in) :: iod,halo,nuc
      integer, intent(in) :: n_bln
! Local parameters:
      character (len=*), parameter :: fname = "aq.nc"
      integer, parameter :: x=1,y=1,noz=1
! Local scalars:
      integer :: i0, i1, i2, id_n, id_nkc_l, id_noz,id_x, id_y, k
      integer :: jddim1(4)
      common /cdf_var_aq/ idaq_rec,idvar_aq(j2+j6+7),idaqfile,
     &     iliqcount,jddim_aq(4)
      integer :: idaq_rec,idvar_aq,idaqfile,iliqcount,jddim_aq

      iliqcount=0

      k=nf_create(fname,nf_clobber,idaqfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,nf_global,'title',26,
     &                  'aqueous phase')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c dimensions

      k=nf_def_dim(idaqfile,'n',n_bln,id_n)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idaqfile,'x',x,id_x)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idaqfile,'y',y,id_y)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idaqfile,'noz',noz,id_noz)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idaqfile,'nkc_l',nkc_l,id_nkc_l)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idaqfile,'rec',nf_unlimited,idaq_rec)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c variables
c time variables

      jddim_aq(1)=id_x
      jddim_aq(2)=id_y
      jddim_aq(3)=id_noz
      jddim_aq(4)=idaq_rec

      k=nf_def_var(idaqfile,'lday',nf_int,4,jddim_aq,idvar_aq(1))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(1),'long_name',3,'day')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idaqfile,'lst',nf_int,4,jddim_aq,idvar_aq(2))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(2),'long_name',4,'hour')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idaqfile,'lmin',nf_int,4,jddim_aq,idvar_aq(3))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(3),'long_name',6,'minute')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c chemical species

c aqueous phase species

      jddim1(1)=id_nkc_l
      jddim1(2)=id_y
      jddim1(3)=id_n
      jddim1(4)=idaq_rec
c offset due to time variables
      i0 = 3

      k=nf_def_var(idaqfile,'NO',nf_float,4,jddim1,idvar_aq(i0+1))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+1),'long_name',18,
     & 'NO (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+1),'units',16,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'NO2',nf_float,4,jddim1,idvar_aq(i0+2))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+2),'long_name',19,
     & 'NO2 (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+2),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'HNO3',nf_float,4,jddim1,idvar_aq(i0+3))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+3),'long_name',20,
     & 'HNO3 (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+3),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'NH3',nf_float,4,jddim1,idvar_aq(i0+4))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+4),'long_name',19,
     & 'NH3 (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+4),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'SO2',nf_float,4,jddim1,idvar_aq(i0+5))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+5),'long_name',19,
     & 'SO2 (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+5),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'H2SO4',nf_float,4,jddim1,idvar_aq(i0+6))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+6),'long_name',21,
     & 'H2SO4 (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+6),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'O3',nf_float,4,jddim1,idvar_aq(i0+7))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+7),'long_name',18,
     & 'O3 (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+7),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'CH4',nf_float,4,jddim1,idvar_aq(i0+8))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+8),'long_name',28,
c     & 'CH4 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+8),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'C2H6',nf_float,4,jddim1,idvar_aq(i0+9))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+9),'long_name',19,
c     & 'C2H6 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+9),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'C3H8',nf_float,4,jddim1,idvar_aq(i0+10))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+10),'long_name',19,
c     & 'C3H8 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+10),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'ALKA',nf_float,4,jddim1,idvar_aq(i0+11))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+11),'long_name',26,
c     & 'Alkane > C3 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+101),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'ETHE',nf_float,4,jddim1,idvar_aq(i0+12))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+12),'long_name',20,
c     & 'Ethene (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+12),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'ALKE',nf_float,4,jddim1,idvar_aq(i0+13))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+13),'long_name',26,
c     & 'Alkene > C2 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+13),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'AROM',nf_float,4,jddim1,idvar_aq(i0+14))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+14),'long_name',28,
c     & 'Alkylbenzenes (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+14),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
      k=nf_def_var(idaqfile,'HCOOH',nf_float,4,jddim1,idvar_aq(i0+15))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+15),'long_name',27,
     & 'Formic acid (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+15),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'ACTA',nf_float,4,jddim1,idvar_aq(i0+16))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+16),'long_name',26,
     & 'Acetic acid (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+16),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'HCHO',nf_float,4,jddim1,idvar_aq(i0+17))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+17),'long_name',28,
     & 'Formaldehyde (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+17),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'ALD2',nf_float,4,jddim1,idvar_aq(i0+18))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+18),'long_name',29,
c     & 'aldehydes > C1 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+18),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'H2O2',nf_float,4,jddim1,idvar_aq(i0+19))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+19),'long_name',33,
     & 'Hydrogen Peroxide (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+19),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'CH3OOH',nf_float,4,jddim1,
     &     idvar_aq(i0+20))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+20),'long_name',36,
     & 'Alkyl hydroperoxides (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+20),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'HONO',nf_float,4,jddim1,idvar_aq(i0+21))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+21),'long_name',28,
     & 'Nitrous acid (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+21),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'PAN',nf_float,4,jddim1,idvar_aq(i0+22))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+22),'long_name',18,
c     & 'PAN (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+22),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'TPAN',nf_float,4,jddim1,idvar_aq(i0+23))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+23),'long_name',19,
c     & 'TPAN (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+23),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'KET',nf_float,4,jddim1,idvar_aq(i0+24))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+24),'long_name',29,
c     & 'Lumped ketones (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+24),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'CRES',nf_float,4,jddim1,idvar_aq(i0+25))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+25),'long_name',21,
c     & 'Cresol (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+25),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'DIAL',nf_float,4,jddim1,idvar_aq(i0+26))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+26),'long_name',37,
c     & 'Unsaturated dicarbonyl (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+26),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'GLYX',nf_float,4,jddim1,idvar_aq(i0+27))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+27),'long_name',22,
c     & 'Glyoxal (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+27),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'MGLY',nf_float,4,jddim1,idvar_aq(i0+28))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+28),'long_name',28,
c     & 'Methylglyoxal (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+28),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'NH4NO3',nf_float,4,jddim1,idvar_aq(i0+29))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+29),'long_name',31,
c     & 'Ammonium nitrate (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+29),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'R3N2',nf_float,4,jddim1,idvar_aq(i0+31))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+31),'long_name',28,
c     & 'Propane RONO2 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+31),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'RAN2',nf_float,4,jddim1,idvar_aq(i0+32))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+32),'long_name',32,
c     & 'alkane > C3 RONO2 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+32),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'RAN1',nf_float,4,jddim1,idvar_aq(i0+33))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+33),'long_name',31,
c     & 'alkane > C3 RNO2 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+33),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'N2O5',nf_float,4,jddim1,idvar_aq(i0+34))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+34),'long_name',35,
c     & 'Dinitrogen pentoxide (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+34),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'HNO4',nf_float,4,jddim1,idvar_aq(i0+35))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+35),'long_name',33,
     & 'Peroxynitric acid (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+35),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'NO3',nf_float,4,jddim1,idvar_aq(i0+36))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+36),'long_name',27,
     & 'NO3 radical (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+36),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'DMS',nf_float,4,jddim1,idvar_aq(i0+37))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+37),'long_name',30,
     & 'Dimethylsulfid (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+37),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      if (halo) then
         k=nf_def_var(idaqfile,'HCl',nf_float,4,jddim1,idvar_aq(i0+30))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+30),'long_name',33,
     &        'Hydrogen chloride (aqueous phase)')
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+30),'units',14,
     &        'mol m-3 (air) ')
         if (k.ne.nf_noerr) call ehandle(k,fname)

         k=nf_def_var(idaqfile,'HOCl',nf_float,4,jddim1,idvar_aq(i0+38))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+38),'long_name',20,
     &        'HOCl (aqueous phase)')
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+38),'units',14,
     &        'mol m-3 (air) ')
         if (k.ne.nf_noerr) call ehandle(k,fname)

c         k=nf_def_var(idaqfile,'ClNO2',nf_float,4,jddim1,idvar_aq(i0+39))
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+39),'long_name',14,
c     &        'ClNO2 (aqueous phase)')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+39),'units',14,
c     &        'mol m-3 (air) ')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c
c         k=nf_def_var(idaqfile,'ClNO3',nf_float,4,jddim1,idvar_aq(i0+40))
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+40),'long_name',20,
c     &        'ClNO3 (aqueous phase)')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+40),'units',14,
c     &        'mol m-3 (air) ')
c         if (k.ne.nf_noerr) call ehandle(k,fname)

         k=nf_def_var(idaqfile,'Cl2',nf_float,4,jddim1,idvar_aq(i0+41))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+41),'long_name',19,
     &        'Cl2 (aqueous phase)')
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+41),'units',14,
     &        'mol m-3 (air) ')
         if (k.ne.nf_noerr) call ehandle(k,fname)

         k=nf_def_var(idaqfile,'HBr',nf_float,4,jddim1,idvar_aq(i0+42))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+42),'long_name',19,
     &        'HBr (aqueous phase)')
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+42),'units',14,
     &        'mol m-3 (air) ')
         if (k.ne.nf_noerr) call ehandle(k,fname)

         k=nf_def_var(idaqfile,'HOBr',nf_float,4,jddim1,idvar_aq(i0+43))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+43),'long_name',20,
     &        'HOBr (aqueous phase)')
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+43),'units',14,
     &        'mol m-3 (air) ')
         if (k.ne.nf_noerr) call ehandle(k,fname)

c         k=nf_def_var(idaqfile,'BrNO2',nf_float,4,jddim1,idvar_aq(i0+44))
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+44),'long_name',20,
c     &        'BrNO2 (aqueous phase)')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+44),'units',14,
c     &        'mol m-3 (air) ')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c
c         k=nf_def_var(idaqfile,'BrNO3',nf_float,4,jddim1,idvar_aq(i0+45))
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+45),'long_name',20,
c     &        'BrNO3 (aqueous phase)')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+45),'units',14,
c     &        'mol m-3 (air) ')
c         if (k.ne.nf_noerr) call ehandle(k,fname)

         k=nf_def_var(idaqfile,'Br2',nf_float,4,jddim1,idvar_aq(i0+46))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+46),'long_name',19,
     &        'Br2 (aqueous phase)')
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+46),'units',14,
     &        'mol m-3 (air) ')
         if (k.ne.nf_noerr) call ehandle(k,fname)

         k=nf_def_var(idaqfile,'BrCl',nf_float,4,jddim1,idvar_aq(i0+47))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+47),'long_name',20,
     &     'BrCl (aqueous phase)')
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+47),'units',14,
     &        'mol m-3 (air) ')
         if (k.ne.nf_noerr) call ehandle(k,fname)

c         k=nf_def_var(idaqfile,'Cl2O2',nf_float,4,jddim1,
c     &        idvar_aq(i0+65))
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+65),'long_name',19,
c     &     'Cl2O2 (aqueous phase)')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+65),'units',14,
c     &        'mol m-3 (air) ')
c         if (k.ne.nf_noerr) call ehandle(k,fname)

c         k=nf_def_var(idaqfile,'Br2O',nf_float,4,jddim1,
c     &        idvar_aq(i0+84))
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+84),'long_name',20,
c     &     'Br2O (aqueous phase)')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+84),'units',14,
c     &        'mol m-3 (air) ')
c         if (k.ne.nf_noerr) call ehandle(k,fname)

c         k=nf_def_var(idaqfile,'ClONO',nf_float,4,jddim1,
c     &        idvar_aq(i0+85))
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+85),'long_name',21,
c     &     'ClONO (aqueous phase)')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+85),'units',14,
c     &        'mol m-3 (air) ')
c         if (k.ne.nf_noerr) call ehandle(k,fname)

c         k=nf_def_var(idaqfile,'ClO3',nf_float,4,jddim1,
c     &        idvar_aq(i0+86))
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+86),'long_name',20,
c     &     'ClO3 (aqueous phase)')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+86),'units',14,
c     &        'mol m-3 (air) ')
c         if (k.ne.nf_noerr) call ehandle(k,fname)

c         k=nf_def_var(idaqfile,'Cl2O3',nf_float,4,jddim1,
c     &        idvar_aq(i0+87))
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+87),'long_name',21,
c     &     'Cl2O3 (aqueous phase)')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+87),'units',14,
c     &        'mol m-3 (air) ')
c         if (k.ne.nf_noerr) call ehandle(k,fname)

c         k=nf_def_var(idaqfile,'RCl',nf_float,4,jddim1,
c     &        idvar_aq(i0+92))
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+92),'long_name',19,
c     &     'RCl (aqueous phase)')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+92),'units',14,
c     &        'mol m-3 (air) ')
c         if (k.ne.nf_noerr) call ehandle(k,fname)

c         k=nf_def_var(idaqfile,'RBr',nf_float,4,jddim1,
c     &        idvar_aq(i0+93))
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+93),'long_name',19,
c     &     'RBr (aqueous phase)')
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c         k=nf_put_att_text(idaqfile,idvar_aq(i0+93),'units',14,
c     &        'mol m-3 (air) ')
c         if (k.ne.nf_noerr) call ehandle(k,fname)

         k=nf_def_var(idaqfile,'XOR',nf_float,4,jddim1,
     &        idvar_aq(i0+94))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+94),'long_name',19,
     &     'XOR (aqueous phase)')
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i0+94),'units',14,
     &        'mol m-3 (air) ')
         if (k.ne.nf_noerr) call ehandle(k,fname)

         if (iod) then
c            k=nf_def_var(idaqfile,'HI',nf_float,4,jddim1,idvar_aq(i0+48))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+48),'long_name',17,
c     &           'HI (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+48),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)

            k=nf_def_var(idaqfile,'HOI',nf_float,4,jddim1,
     &       idvar_aq(i0+49))
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i0+49),'long_name',19,
     &           'HOI (aqueous phase)')
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i0+49),'units',14,
     &           'mol m-3 (air) ')
            if (k.ne.nf_noerr) call ehandle(k,fname)

c            k=nf_def_var(idaqfile,'I2O2',nf_float,4,jddim1,
c     &        idvar_aq(i0+50))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+50),'long_name',19,
c     &           'I2O2 (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+50),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c
c            k=nf_def_var(idaqfile,'INO2',nf_float,4,jddim1,
c     &        idvar_aq(i0+51))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+51),'long_name',19,
c     &           'INO2 (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+51),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c
c            k=nf_def_var(idaqfile,'INO3',nf_float,4,jddim1,
c     &        idvar_aq(i0+52))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+52),'long_name',19,
c     &           'INO3 (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+52),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)

            k=nf_def_var(idaqfile,'I2',nf_float,4,jddim1,
     &       idvar_aq(i0+53))
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i0+53),'long_name',18,
     &           'I2 (aqueous phase)')
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i0+53),'units',14,
     &           'mol m-3 (air) ')
            if (k.ne.nf_noerr) call ehandle(k,fname)

            k=nf_def_var(idaqfile,'ICl',nf_float,4,jddim1,
     &        idvar_aq(i0+54))
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i0+54),'long_name',19,
     &           'ICl (aqueous phase)')
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i0+54),'units',14,
     &           'mol m-3 (air) ')
            if (k.ne.nf_noerr) call ehandle(k,fname)

            k=nf_def_var(idaqfile,'IBr',nf_float,4,jddim1,
     &           idvar_aq(i0+55))
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i0+55),'long_name',19,
     &           'IBr (aqueous phase)')
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i0+55),'units',14,
     &           'mol m-3 (air) ')
            if (k.ne.nf_noerr) call ehandle(k,fname)

c            k=nf_def_var(idaqfile,'CH3I',nf_float,4,jddim1,
c     &        idvar_aq(i0+56))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+56),'long_name',19,
c     &           'CH3I (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+56),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c
c            k=nf_def_var(idaqfile,'CH2I2',nf_float,4,jddim1,
c     &       idvar_aq(i0+57))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+57),'long_name',20,
c     &           'CH2I2 (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+57),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c
c            k=nf_def_var(idaqfile,'CH2ClI',nf_float,4,jddim1,
c     &       idvar_aq(i0+58))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+58),'long_name',23,
c     &           'CH2ClI (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+58),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c
c            k=nf_def_var(idaqfile,'C3H7I',nf_float,4,jddim1,
c     &         idvar_aq(i0+59))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+59),'long_name',20,
c     &           'C3H7I (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+59),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c
c            k=nf_def_var(idaqfile,'CH2BrI',nf_float,4,jddim1,
c     &         idvar_aq(i0+71))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+71),'long_name',20,
c     &           'CH2BrI (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+71),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c
c            k=nf_def_var(idaqfile,'CHBr2I',nf_float,4,jddim1,
c     &         idvar_aq(i0+72))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+72),'long_name',20,
c     &           'CHBr2I (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+72),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c
c            k=nf_def_var(idaqfile,'C2H5I',nf_float,4,jddim1,
c     &         idvar_aq(i0+73))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+73),'long_name',20,
c     &           'C2H5I (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+73),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c
c            k=nf_def_var(idaqfile,'HIO3',nf_float,4,jddim1,idvar_aq
c     &           (i0+74))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+74),'long_name',20,
c     &           'HIO3 (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+74),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c
c            k=nf_def_var(idaqfile,'I2O',nf_float,4,jddim1,idvar_aq
c     &           (i0+79))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+79),'long_name',19,
c     &           'I2O (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+79),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)

c            k=nf_def_var(idaqfile,'I2O3',nf_float,4,jddim1,idvar_aq
c     &           (i0+80))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+80),'long_name',20,
c     &           'I2O3 (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+80),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)

c            k=nf_def_var(idaqfile,'I2O4',nf_float,4,jddim1,idvar_aq
c     &           (i0+81))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+81),'long_name',20,
c     &           'I2O4 (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+81),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)

c            k=nf_def_var(idaqfile,'I2O5',nf_float,4,jddim1,idvar_aq
c     &           (i0+82))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+82),'long_name',20,
c     &           'I2O5 (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+82),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)

c            k=nf_def_var(idaqfile,'INO',nf_float,4,jddim1,idvar_aq
c     &           (i0+83))
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+83),'long_name',19,
c     &           'INO (aqueous phase)')
c            if (k.ne.nf_noerr) call ehandle(k,fname)
c            k=nf_put_att_text(idaqfile,idvar_aq(i0+83),'units',14,
c     &           'mol m-3 (air) ')
c            if (k.ne.nf_noerr) call ehandle(k,fname)

         endif  ! iod
      endif     ! halo

      k=nf_def_var(idaqfile,'DMSO',nf_float,4,jddim1,idvar_aq(i0+60))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+60),'long_name',20,
     & 'DMSO (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+60),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'CH3SO2',nf_float,4,jddim1,idvar_aq(i0+61))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+61),'long_name',21,
c     & 'CH3SO2 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+61),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'CH3SO3',nf_float,4,jddim1,idvar_aq(i0+62))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+62),'long_name',21,
c     & 'CH3SO3 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+62),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'MSA',nf_float,4,jddim1,idvar_aq(i0+63))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+63),'long_name',22,
c     & 'CH3SO3H (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+63),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'CO',nf_float,4,jddim1,idvar_aq(i0+64))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+64),'long_name',21,
c     & 'CO (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+64),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'DMOO',nf_float,4,jddim1,idvar_aq(i0+66))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+66),'long_name',22,
c     & 'CH3SCH3OO (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+66),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'CH3S',nf_float,4,jddim1,idvar_aq(i0+67))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+67),'long_name',22,
c     & 'CH3S (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+67),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'CH3SO',nf_float,4,jddim1,idvar_aq(i0+68))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+68),'long_name',22,
c     & 'CH3SO (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+68),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'MSIA',nf_float,4,jddim1,idvar_aq(i0+69))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+69),'long_name',22,
c     & 'CH3S(O)OH (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+69),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'DMSO2',nf_float,4,jddim1,idvar_aq(i0+70))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+70),'long_name',29,
     & 'CH3S(O)(O)CH3 (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+70),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

! jjb v733 uncommented, v741 commented
c      if (nuc) then
c        k=nf_def_var(idaqfile,'NUCV',nf_float,4,jddim1,idvar_aq(i0+75))
c        if (k.ne.nf_noerr) call ehandle(k,fname)
c        k=nf_put_att_text(idaqfile,idvar_aq(i0+75),'long_name',20,
c     &       'NUCV (aqueous phase)')
c        if (k.ne.nf_noerr) call ehandle(k,fname)
c        k=nf_put_att_text(idaqfile,idvar_aq(i0+75),'units',14,
c     &       'mol m-3 (air) ')
c        if (k.ne.nf_noerr) call ehandle(k,fname)
c      endif

c      k=nf_def_var(idaqfile,'SO3',nf_float,4,jddim1,idvar_aq(i0+76))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+76),'long_name',19,
c     & 'SO3 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+76),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'HOSO2',nf_float,4,jddim1,idvar_aq(i0+77))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+77),'long_name',21,
c     & 'HOSO2 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+77),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'CO2',nf_float,4,jddim1,idvar_aq(i0+78))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+78),'long_name',19,
     & 'CO2 (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+78),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'CH3OH',nf_float,4,jddim1,idvar_aq(i0+88))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+88),'long_name',24,
     & 'Methanol (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+88),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'C2H5OH',nf_float,4,jddim1,idvar_aq(i0+89))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+89),'long_name',23,
     & 'Ethanol (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+89),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'H2',nf_float,4,jddim1,idvar_aq(i0+90))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+90),'long_name',18,
c     & 'H2 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+90),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'NHS',nf_float,4,jddim1,idvar_aq(i0+91))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+91),'long_name',19,
c     & 'NHS (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+91),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'SOR',nf_float,4,jddim1,idvar_aq(i0+95))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+95),'long_name',19,
     & 'SOR (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i0+95),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'SPAN',nf_float,4,jddim1,idvar_aq(i0+96))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+96),'long_name',20,
c     & 'SPAN (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+96),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'Hg',nf_float,4,jddim1,idvar_aq(i0+97))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+97),'long_name',21,
c     & 'Hg(O) (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+97),'units',7,
c     & 'mol m-3')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgO',nf_float,4,jddim1,idvar_aq(i0+98))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+98),'long_name',21,
c     & 'HgO (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+98),'units',7,
c     & 'mol m-3')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgCl',nf_float,4,jddim1,idvar_aq(i0+99))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+99),'long_name',22,
c     & 'HgCl (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+99),'units',7,
c     & 'mol m-3')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgCl2',nf_float,4,jddim1,idvar_aq(i0+100))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+100),'long_name',23,
c     & 'HgCl2 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+100),'units',7,
c     & 'mol m-3')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgBr',nf_float,4,jddim1,idvar_aq(i0+101))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+101),'long_name',22,
c     & 'HgBr (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+101),'units',7,
c     & 'mol m-3')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgBr2',nf_float,4,jddim1,idvar_aq(i0+102))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+102),'long_name',23,
c     & 'HgBr2 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i0+102),'units',7,
c     & 'mol m-3')
c      if (k.ne.nf_noerr) call ehandle(k,fname)



c radicals and "other", the "j3-species"
c offset: i0+j1 = i1 ; note that there are diffs to gas phase!

      i1 = i0+j1

c index 1 is empty!

      k=nf_def_var(idaqfile,'OH',nf_float,4,jddim1,idvar_aq(i1+2))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i1+2),'long_name',32,
     & 'Hydroxyl Radical (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i1+2),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'HO2',nf_float,4,jddim1,idvar_aq(i1+3))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i1+3),'long_name',36,
     & 'Hydroperoxyl radical (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i1+3),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'DOM',nf_float,4,jddim1,idvar_aq(i1+4))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i1+4),'long_name',33,
     & 'Dis. Org. Matter  (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i1+4),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'CH3OO',nf_float,4,jddim1,idvar_aq(i1+6))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i1+6),'long_name',21,
     & 'CH3O2 (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i1+6),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      if (halo) then

         if (iod) then

            k=nf_def_var(idaqfile,'HIO2',nf_float,4,jddim1,
     &           idvar_aq(i1+5))
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i1+5),'long_name',20,
     &           'HIO2 (aqueous phase)')
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i1+5),'units',14,
     &           'mol m-3 (air) ')
            if (k.ne.nf_noerr) call ehandle(k,fname)

            k=nf_def_var(idaqfile,'IO',nf_float,4,jddim1,
     &           idvar_aq(i1+7))
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i1+7),'long_name',17,
     &           'IO (aqueous phase)')
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i1+7),'units',14,
     &           'mol m-3 (air) ')
            if (k.ne.nf_noerr) call ehandle(k,fname)

            k=nf_def_var(idaqfile,'OIO',nf_float,4,jddim1,idvar_aq
     &           (i1+12))
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i1+12),'long_name',19,
     &           'OIO (aqueous phase)')
            if (k.ne.nf_noerr) call ehandle(k,fname)
            k=nf_put_att_text(idaqfile,idvar_aq(i1+12),'units',14,
     &           'mol m-3 (air) ')
            if (k.ne.nf_noerr) call ehandle(k,fname)

         endif  ! iod

         k=nf_def_var(idaqfile,'Cl',nf_float,4,jddim1,idvar_aq(i1+8))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i1+8),'long_name',18,
     &        'Cl (aqueous phase)')
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i1+8),'units',14,
     &        'mol m-3 (air) ')
         if (k.ne.nf_noerr) call ehandle(k,fname)

         k=nf_def_var(idaqfile,'Br',nf_float,4,jddim1,idvar_aq(i1+9))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i1+9),'long_name',18,
     &        'Br (aqueous phase)')
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idaqfile,idvar_aq(i1+9),'units',14,
     &        'mol m-3 (air) ')
         if (k.ne.nf_noerr) call ehandle(k,fname)

      endif   ! halo

      k=nf_def_var(idaqfile,'O2',nf_float,4,jddim1,idvar_aq(i1+11))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i1+11),'long_name',18,
     & 'O2 (aqueous phase)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i1+11),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'HgOH2',nf_float,4,jddim1,idvar_aq(i1+13))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i1+13),'long_name',21,
c     & 'HgOH2 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i1+13),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgOHCl',nf_float,4,jddim1,idvar_aq(i1+14))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i1+14),'long_name',22,
c     & 'HgOHCl (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i1+14),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgSO3',nf_float,4,jddim1,idvar_aq(i1+15))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i1+15),'long_name',21,
c     & 'HgSO3 (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i1+15),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgOHBr',nf_float,4,jddim1,idvar_aq(i1+16))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i1+16),'long_name',18,
c     & 'HgOHBr (aqueous phase)')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i1+16),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)


c ions

c rvg: it's NOT possible to replace "p" and "m" with"+" and "-"

c offset due to time variables and non-ionic species:
      i2=j1+j3+i0


      k=nf_def_var(idaqfile,'Hp',nf_float,4,jddim1,idvar_aq(i2+1))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+1),'long_name',3,
     & 'H+ ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+1),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'NH4p',nf_float,4,jddim1,idvar_aq(i2+2))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+2),'long_name',5,
     & 'NH4+ ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+2),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'OHm',nf_float,4,jddim1,idvar_aq(i2+3))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+3),'long_name',4,
     & 'OH- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+3),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'CH2OHSO3m',nf_float,4,jddim1,
     & idvar_aq(i2+4))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+4),'long_name',10,
     & 'CH2OHSO3- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+4),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'HSO3m',nf_float,4,jddim1,idvar_aq(i2+5))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+5),'long_name',6,
     & 'HSO3- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+5),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'SO32m',nf_float,4,jddim1,idvar_aq(i2+6))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+6),'long_name',5,
     & 'SO3= ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+6),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'SO4m',nf_float,4,jddim1,idvar_aq(i2+7))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+7),'long_name',5,
     & 'SO4- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+7),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'SO42m',nf_float,4,jddim1,idvar_aq(i2+8))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+8),'long_name',5,
     & 'SO4= ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+8),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'HCO3m',nf_float,4,jddim1,idvar_aq(i2+9))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+9),'long_name',6,
     & 'HCO3- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+9),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'CO3m',nf_float,4,jddim1,idvar_aq(i2+10))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+10),'long_name',5,
     & 'CO3- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+10),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'O2m',nf_float,4,jddim1,idvar_aq(i2+11))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+11),'long_name',4,
     & 'O2- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+11),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'NO2m',nf_float,4,jddim1,idvar_aq(i2+12))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+12),'long_name',5,
     & 'NO2- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+12),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'NO3m',nf_float,4,jddim1,idvar_aq(i2+13))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+13),'long_name',5,
     & 'NO3- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+13),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'Clm',nf_float,4,jddim1,idvar_aq(i2+14))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+14),'long_name',4,
     & 'Cl- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+14),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'Cl2m',nf_float,4,jddim1,idvar_aq(i2+15))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+15),'long_name',5,
     & 'Cl2- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+15),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'HCOOm',nf_float,4,jddim1,idvar_aq(i2+16))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+16),'long_name',6,
     & 'HCOO- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+16),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'Fe3p',nf_float,4,jddim1,idvar_aq(i2+17))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+17),'long_name',5,
c     & 'Fe3+ ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+17),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'Mn2p',nf_float,4,jddim1,idvar_aq(i2+18))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+18),'long_name',5,
c     & 'Mn2+ ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+18),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'HSO4m',nf_float,4,jddim1,idvar_aq(i2+19))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+19),'long_name',6,
     & 'HSO4- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+19),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'Nap',nf_float,4,jddim1,idvar_aq(i2+20))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+20),'long_name',38,
     & 'Na+ (sum of all additional pos. ions')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+20),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'NO4m',nf_float,4,jddim1,idvar_aq(i2+21))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+21),'long_name',4,
     & 'NO4- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+21),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      if (halo) then

      k=nf_def_var(idaqfile,'ClOm',nf_float,4,jddim1,idvar_aq(i2+22))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+22),'long_name',5,
     & 'ClO- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+22),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'ClOHm',nf_float,4,jddim1,idvar_aq(i2+23))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+23),'long_name',6,
     & 'ClOH- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+23),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'Brm',nf_float,4,jddim1,idvar_aq(i2+24))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+24),'long_name',4,
     & 'Br- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+24),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'Br2m',nf_float,4,jddim1,idvar_aq(i2+25))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+25),'long_name',5,
     & 'Br2- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+25),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'BrOm',nf_float,4,jddim1,idvar_aq(i2+26))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+26),'long_name',5,
     & 'BrO- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+26),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'BrOHm',nf_float,4,jddim1,idvar_aq(i2+27))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+27),'long_name',6,
     & 'BrOH- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+27),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'BrCl2m',nf_float,4,jddim1,
     &     idvar_aq(i2+28))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+28),'long_name',7,
     & 'BrCl2- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+28),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'Br2Clm',nf_float,4,jddim1,
     &     idvar_aq(i2+29))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+29),'long_name',7,
     & 'Br2Cl- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+29),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      endif  ! halo

      k=nf_def_var(idaqfile,'CH3SO3m',nf_float,4,jddim1,
     &     idvar_aq(i2+30))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+30),'long_name',8,
     & 'CH3SO3- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+30),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'HSO5m',nf_float,4,jddim1,idvar_aq(i2+31))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+31),'long_name',6,
     & 'HSO5- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+31),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'SO3m',nf_float,4,jddim1,idvar_aq(i2+32))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+32),'long_name',5,
     & 'SO3- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+32),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'SO5m',nf_float,4,jddim1,idvar_aq(i2+33))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+33),'long_name',5,
     & 'SO5- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+33),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      if (iod) then

      k=nf_def_var(idaqfile,'Im',nf_float,4,jddim1,idvar_aq(i2+34))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+34),'long_name',3,
     & 'I- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+34),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'IO2m',nf_float,4,jddim1,idvar_aq(i2+35))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+35),'long_name',5,
     & 'IO2- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+35),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'IO3m',nf_float,4,jddim1,idvar_aq(i2+36))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+36),'long_name',5,
     & 'IO3- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+36),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'ICl2m',nf_float,4,jddim1,idvar_aq(i2+37))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+37),'long_name',6,
     & 'ICl2- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+37),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idaqfile,'IBr2m',nf_float,4,jddim1,idvar_aq(i2+38))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+38),'long_name',6,
     & 'IBr2- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+38),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      endif  ! iod

      k=nf_def_var(idaqfile,'CH3SO2m',nf_float,4,jddim1,
     &     idvar_aq(i2+39))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+39),'long_name',7,
     & 'CH3S(O)O- ')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+39),'units',14,
     & 'mol m-3 (air) ')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'Hgp',nf_float,4,jddim1,idvar_aq(i2+40))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+40),'long_name',3,
c     & 'Hg+')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+40),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'Hg2p',nf_float,4,jddim1,idvar_aq(i2+41))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+41),'long_name',4,
c     & 'Hg2+')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+41),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgOHp',nf_float,4,jddim1,idvar_aq(i2+42))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+42),'long_name',5,
c     & 'HgOH+')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+42),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgClp',nf_float,4,jddim1,idvar_aq(i2+43))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+43),'long_name',5,
c     & 'HgCl+')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+43),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgCl3m',nf_float,4,jddim1,idvar_aq(i2+44))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+44),'long_name',6,
c     & 'HgCl3-')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+44),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgCl42m',nf_float,4,jddim1,idvar_aq(i2+45))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+45),'long_name',7,
c     & 'HgCl42-')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+45),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgBrp',nf_float,4,jddim1,idvar_aq(i2+46))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+46),'long_name',5,
c     & 'HgBr+')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+46),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgBr3m',nf_float,4,jddim1,idvar_aq(i2+47))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+47),'long_name',6,
c     & 'HgBr3-')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+47),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgBr42m',nf_float,4,jddim1,idvar_aq(i2+48))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+48),'long_name',7,
c     & 'HgBr42-')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+48),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c      k=nf_def_var(idaqfile,'HgSO322m',nf_float,4,jddim1,
c     &     idvar_aq(i2+49))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+49),'long_name',10,
c     & 'Hg(SO3)22-')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+49),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

c      k=nf_def_var(idaqfile,'',nf_float,4,jddim1,idvar_aq(i2+50))
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+50),'long_name',5,
c     & ' ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)
c      k=nf_put_att_text(idaqfile,idvar_aq(i2+50),'units',14,
c     & 'mol m-3 (air) ')
c      if (k.ne.nf_noerr) call ehandle(k,fname)

c LWC per aqueous bin
      k=nf_def_var(idaqfile,'cw',nf_float,4,jddim1,idvar_aq(i2+41))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+41),'long_name',20,
     & 'liquid water content')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+41),'units',15,
     & 'm3(aq) m-3(air)')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c radius
      k=nf_def_var(idaqfile,'rc',nf_float,4,jddim1,idvar_aq(i2+42))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+42),'long_name',15,
     & 'particle radius')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idaqfile,idvar_aq(i2+42),'units',1,
     & 'm')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c end define mode
      k=nf_enddef(idaqfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine open_chem_aq

c
c----------------------------------------------------------------
c

      subroutine open_jrate (n_bln, chamber)

! jjb work done
!     - use module (and solve inconsistency in ph_rates numbers
!     - implicit dp was missing --> implicit none with missing declarations

      USE global_params, ONLY :
     &     n_jrates=>nphrxn

      implicit none

! Include statements:
      include 'netcdf.inc'

! Subroutine arguments
! Scalar arguments with intent(in):
      integer, intent(in) :: n_bln
      logical, intent(in) :: chamber

! Local parameters:
      character (len=*), parameter :: fname = "jrate.nc"
      integer, parameter :: x=1, y=1, noz=1
! Local scalars:
      integer :: id_n, id_noz, id_x, id_y, ispec
      integer :: k
! Local arrays:
      integer :: jddim1(4)
      character (len=8) :: j_name(n_jrates)
      data j_name/
     & 'J_NO2','J_NOO2','J_O1D','J_HONO','J_HNO3','J_H2O2','J_HNO4',
     & 'J_HCHO','J_COH2','J_NO2O','J_HNO4_2','J_N2O5','J_HOCl',
     & 'J_ClONO2','J_BrNO3','J_Cl2O2','J_CH3OOH','J_ClNO2','J_Cl2',
     & 'J_HOBr','J_BrNO2','J_Br2','J_BrCl','J_BrO','J_IO','J_HOI','J_I2'
     & ,'J_ICl','J_IBr','J_INO3','J_CH3I','J_C3H7I','J_CH2ClI','J_CH2I2'
     & ,'J_OClO','J_I2O2','J_INO2','J_NO2m','J_NO3n','J_OIO',        ! jjb NO3m ???
!    & 'J_dumm3','J_dumm4','J_dumm5','J_CH2BrI','J_CHBr2I','J_C2H5I'/
     & 'J_dumm3','J_dumm4','J_dumm5','J_CH2BrI','J_CHBr2I','J_C2H5I',
     & 'J_O3P'/ ! jjb

! Common blocks:
      common /cdf_var_jrat/ idjrat_rec,idvar_jrat(n_jrates+3),
     &   idjratfile,ijratcount,jddim_jrat(4)
      integer idjrat_rec, idvar_jrat, idjratfile, ijratcount, jddim_jrat


      ijratcount=0

      k=nf_create(fname,nf_clobber,idjratfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idjratfile,nf_global,'title',16,
     &                  'photolysis rates')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c dimensions

      k=nf_def_dim(idjratfile,'n',n_bln,id_n)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idjratfile,'x',x,id_x)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idjratfile,'y',y,id_y)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idjratfile,'noz',noz,id_noz)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idjratfile,'rec',nf_unlimited,idjrat_rec)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      jddim_jrat(1)=id_x
      jddim_jrat(2)=id_y
      jddim_jrat(3)=id_noz
      jddim_jrat(4)=idjrat_rec

c variables
c time variables
      k=nf_def_var(idjratfile,'lday',nf_int,4,jddim_jrat,idvar_jrat(1))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idjratfile,idvar_jrat(1),'long_name',3,'day')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idjratfile,'lst',nf_int,4,jddim_jrat,idvar_jrat(2))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idjratfile,idvar_jrat(2),'long_name',4,'hour')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idjratfile,'lmin',nf_int,4,jddim_jrat,idvar_jrat(3))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idjratfile,idvar_jrat(3),'long_name',6,'minute')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c photolysis jrates

      jddim1(1)=id_x
      jddim1(2)=id_y
      if (.not.chamber) then
         jddim1(3)=id_n
      else
         jddim1(3)=id_noz
      end if
      jddim1(4)=idjrat_rec

      do ispec=1,n_jrates
         k=nf_def_var(idjratfile,j_name(ispec),nf_float,4,jddim1,
     &        idvar_jrat(ispec+3))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idjratfile,idvar_jrat(ispec+3),'long_name',8,
     &     j_name(ispec))
         if (k.ne.nf_noerr) call ehandle(k,fname)
         k=nf_put_att_text(idjratfile,idvar_jrat(ispec+3),'units',3,
     &        '1/s')
         if (k.ne.nf_noerr) call ehandle(k,fname)
      enddo

c end define mode
       k=nf_enddef(idjratfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine open_jrate

c
c----------------------------------------------------------------
c

      subroutine open_rxn
c open netcdf file for reaction rates

      USE global_params, ONLY :
! Imported Parameters:
     &     nlev,
     &     nrxn

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = "rxnrate.nc"
      integer, parameter :: x=1, y=1, noz=1
! Local scalars:
      integer :: id_nrxn, id_n, id_noz, id_x, id_y
      integer :: k
! Local arrays:
      integer :: jddim1(4)
! Common blocks:
      common /cdf_var_rxn/ idrxn_rec,idvar_rxn(4),idrxnfile,
     &   irxncount,jddim_rxn(4)
      integer :: idrxn_rec,idvar_rxn,idrxnfile,irxncount,jddim_rxn


      irxncount=0

      k=nf_create(fname,nf_clobber,idrxnfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idrxnfile,nf_global,'title',16,
     &                  'reaction rates')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c dimensions

      k=nf_def_dim(idrxnfile,'nrxn',nrxn,id_nrxn)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idrxnfile,'nlev',nlev,id_n)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idrxnfile,'x',x,id_x)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idrxnfile,'y',y,id_y)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idrxnfile,'noz',noz,id_noz)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idrxnfile,'rec',nf_unlimited,idrxn_rec)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      jddim_rxn(1)=id_x
      jddim_rxn(2)=id_y
      jddim_rxn(3)=id_noz
      jddim_rxn(4)=idrxn_rec

c variables
c time variables
      k=nf_def_var(idrxnfile,'lday',nf_int,4,jddim_rxn,idvar_rxn(1))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idrxnfile,idvar_rxn(1),'long_name',3,'day')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idrxnfile,'lst',nf_int,4,jddim_rxn,idvar_rxn(2))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idrxnfile,idvar_rxn(2),'long_name',4,'hour')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idrxnfile,'lmin',nf_int,4,jddim_rxn,idvar_rxn(3))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idrxnfile,idvar_rxn(3),'long_name',6,'minute')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c reaction rates

      jddim1(1)=id_x
      jddim1(2)=id_nrxn
      jddim1(3)=id_n
      jddim1(4)=idrxn_rec

      k=nf_def_var(idrxnfile,'Rconst',nf_float,4,jddim1,
     &        idvar_rxn(4))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idrxnfile,idvar_rxn(4),'long_name',14,
     &     'reaction rates')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idrxnfile,idvar_rxn(4),'units',11,'mol m-3 s-1')
      if (k.ne.nf_noerr) call ehandle(k,fname)

c end define mode
       k=nf_enddef(idrxnfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine open_rxn

c
c----------------------------------------------------------------

      subroutine open_nuc
! initialize plot file for nucleation output

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nkt

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = "nuc.nc"
      integer, parameter :: x=1, y=1, noz=1
! Local scalars:
      integer :: id_n, id_nf, id_nkt, id_noz, id_x, id_y, k
! Local arrays:
      integer :: jddim1(4)
! Common blocks:
      common /cdf_var_nuc/ idnuc_rec,idvar_nuc(23),idnucfile,
     &   inuccount,jddim_nuc(4)
      integer :: idnuc_rec,idvar_nuc,idnucfile,inuccount,jddim_nuc

      inuccount=0

      k=nf_create(fname,nf_clobber,idnucfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,nf_global,'title',19,
     &                  'nucleation')
      if (k.ne.nf_noerr) call ehandle(k,fname)
c dimension
      k=nf_def_dim(idnucfile,'x',x,id_x)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idnucfile,'y',y,id_y)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idnucfile,'noz',noz,id_noz)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idnucfile,'n',n,id_n)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idnucfile,'nf',nf,id_nf)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idnucfile,'nkt',nkt,id_nkt)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_dim(idnucfile,'rec',nf_unlimited,idnuc_rec)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c variables
c time variables

      jddim_nuc(1)=id_x
      jddim_nuc(2)=id_y
      jddim_nuc(3)=id_noz
      jddim_nuc(4)=idnuc_rec

      k=nf_def_var(idnucfile,'lday',nf_int,4,jddim_nuc,idvar_nuc(1))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(1),'long_name',3,'day')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idnucfile,'lst',nf_int,4,jddim_nuc,idvar_nuc(2))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(2),'long_name',4,'hour')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_def_var(idnucfile,'lmin',nf_int,4,jddim_nuc,idvar_nuc(3))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(3),'long_name',6,'minute')
      if (k.ne.nf_noerr) call ehandle(k,fname)
c
c nucleation parameters

      jddim1(1)=id_x
      jddim1(2)=id_y
      jddim1(3)=id_n
      jddim1(4)=idnuc_rec

      k=nf_def_var(idnucfile,'xn_new',nf_float,4,jddim1,idvar_nuc(4))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(4),'long_name',23,
     & 'ternary nucleation rate')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(4),'units',13,
     & 'molec/(cm3*s)')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'xn_acc',nf_float,4,jddim1,idvar_nuc(5))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(5),'long_name',30,
     & 'accum. ternary nucleation rate')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(5),'units',9,
     & 'molec/cm3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'xv_acc',nf_float,4,jddim1,idvar_nuc(6))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(6),'long_name',31,
     & 'accum. ternary nucleated volume')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(6),'units',7,
     & 'nm3/cm3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'xn_newio',nf_float,4,jddim1,idvar_nuc(7))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(7),'long_name',23,
     & 'hom.OIO nucleation rate')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(7),'units',13,
     & 'molec/(cm3*s)')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'xn_accio',nf_float,4,jddim1,idvar_nuc(8))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(8),'long_name',30,
     & 'accum. hom.OIO nucleation rate')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(8),'units',9,
     & 'molec/cm3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'xv_accio',nf_float,4,jddim1,idvar_nuc(9))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(9),'long_name',31,
     & 'accum. hom.OIO nucleated volume')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(9),'units',7,
     & 'nm3/cm3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'xn_app',nf_float,4,jddim1,idvar_nuc(10))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(10),'long_name',24,
     & 'apparent nucleation rate')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(10),'units',13,
     & 'molec/(cm3*s)')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'xn_apacc',nf_float,4,jddim1,idvar_nuc(11))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(11),'long_name',31,
     & 'accum. apparent nucleation rate')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(11),'units',9,
     & 'molec/cm3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'xv_apacc',nf_float,4,jddim1,idvar_nuc(12))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(12),'long_name',32,
     & 'accum. apparent nucleated volume')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(12),'units',7,
     & 'nm3/cm3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'bn_ges',nf_float,4,jddim1,idvar_nuc(13))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(13),'long_name',26,
     & 'background particle number')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(13),'units',4,
     & '/cm3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'bd_mean',nf_float,4,jddim1,idvar_nuc(14))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(14),'long_name',33,
     & 'mean background particle diameter')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(14),'units',2,
     & 'nm')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'dnh3',nf_float,4,jddim1,idvar_nuc(15))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(15),'long_name',24,
     & 'delta NH3 (before-after)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(15),'units',3,
     & 'ppt')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'dh2so4',nf_float,4,jddim1,idvar_nuc(16))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(16),'long_name',26,
     & 'delta h2so4 (before-after)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(16),'units',9,
     & 'molec/cm3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'doio',nf_float,4,jddim1,idvar_nuc(17))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(17),'long_name',24,
     & 'delta OIO (before-after)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(17),'units',3,
     & 'ppt')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'dnucv',nf_float,4,jddim1,idvar_nuc(18))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(18),'long_name',25,
     & 'delta NUCV (before-after)')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(18),'units',3,
     & 'ppt')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'grorate',nf_float,4,jddim1,idvar_nuc(19))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(19),'long_name',12,
     & 'growth rate')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(19),'units',4,
     & 'nm/h')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'concnuc',nf_float,4,jddim1,idvar_nuc(20))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(20),'long_name',25,
     & 'real nuclei concentration')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(20),'units',10,
     & 'nuclei/cm3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'partsa',nf_float,4,jddim1,idvar_nuc(21))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(21),'long_name',21,
     & 'particle surface area')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(21),'units',7,
     & 'um2/cm3')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      jddim1(1)=id_x
      jddim1(2)=id_nkt
      jddim1(3)=id_n
      jddim1(4)=idnuc_rec

      k=nf_def_var(idnucfile,'partd',nf_float,4,jddim1,idvar_nuc(22))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(22),'long_name',25,
     & 'background part. diameter')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(22),'units',2,
     & 'nm')
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_def_var(idnucfile,'partNu',nf_float,4,jddim1,idvar_nuc(23))
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(23),'long_name',25,
     & 'background part. number')
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_put_att_text(idnucfile,idvar_nuc(23),'units',5,
     & '1/cm3')
      if (k.ne.nf_noerr) call ehandle(k,fname)


c end define mode
      k=nf_enddef(idnucfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine open_nuc


c
c----------------------------------------------------------------
c

      subroutine write_grid

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nka,
     &     nkt

      USE precision, ONLY :
! Imported Parameters:
     &     dp

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = "grid.nc"
! Local scalars:
      integer :: ia, jt, k
! Local arrays:
      integer :: idimcount(4), idimstart(4)
      real (kind=dp) :: field(1,1,n),field2(nka,1,1),field3(1,nkt,1),
     &     field4(nka,nkt,1)

! Common blocks:
      common /cdf_var_grid/ id_rec,idvar_g(9),idfile,icount,jddim(4)
      integer :: id_rec, idvar_g, idfile, icount, jddim
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      real (kind=dp) :: detw, deta, eta, etw

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      real (kind=dp) :: enw, ew, rn, rw, en, e, dew, rq

! == End of declarations =======================================================

      icount=icount+1

      idimcount(1)=1
      idimcount(2)=1
      idimcount(3)=n
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=icount

      field(1,1,:)=eta(1:n)
      k=nf_put_vara_double(idfile,idvar_g(1),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=etw(1:n)
      k=nf_put_vara_double(idfile,idvar_g(2),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      idimcount(1)=nka
      idimcount(3)=1

      field2(:,1,1)=enw(1:nka)
      k=nf_put_vara_double(idfile,idvar_g(3),idimstart,idimcount,field2)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field2(:,1,1)=en(1:nka)
      k=nf_put_vara_double(idfile,idvar_g(4),idimstart,idimcount,field2)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      idimcount(1)=1
      idimcount(2)=nkt

      field3(1,:,1)=ew(1:nkt)
      k=nf_put_vara_double(idfile,idvar_g(5),idimstart,idimcount,field3)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field3(1,:,1)=e(1:nkt)
      k=nf_put_vara_double(idfile,idvar_g(6),idimstart,idimcount,field3)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      idimcount(1)=nka
      idimcount(2)=1

      field2(:,1,1)=rn(1:nka)
      k=nf_put_vara_double(idfile,idvar_g(7),idimstart,idimcount,field2)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      idimcount(1)=nka
      idimcount(2)=nkt

      do jt=1,nkt
         do ia=1,nka
            field4(ia,jt,1)=rq(jt,ia)
         end do
      end do
      k=nf_put_vara_double(idfile,idvar_g(8),idimstart,idimcount,field4)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      do jt=1,nkt
         do ia=1,nka
            field4(ia,jt,1)=rw(jt,ia)
         end do
      end do
      k=nf_put_vara_double(idfile,idvar_g(9),idimstart,idimcount,field4)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_sync(idfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)


c close file
      k=nf_close(idfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine write_grid



c
c----------------------------------------------------------------
c

      subroutine write_met (n_bln)
c output of meteorological variables

! jjb work done
!     - corrected character length for fname to be consistent with SR ehandle
!     - replaced idexes in CB lang, kurz for homogeneity

! 21-Sep-2020   Josue Bock   Minor bugfix: correct the third dimension of local array 'field': n_bln instead of n

! 19-May-2021   Josue Bock   Bugfix: correct lcl/lct netCDF function, nf_put_var1_int replaced by nf_put_vara_int
!                            Missing declarations and implicit none

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nrlay,
     &     nrlev,
     &     nka,
     &     nkt

      USE precision, ONLY :
! Imported Parameters:
     &     dp

      implicit none

! Include statements:
      include 'netcdf.inc'

      integer, intent(in) :: n_bln
! Local parameters:
      character (len=*), parameter :: fname = "meteo.nc"
! Local scalars:
      integer :: k, kk
      real (kind=dp) :: ddz
! Local arrays:
      integer :: ifield(1,1,1), idimcount(4), idimstart(4)
      real (kind=dp) :: field(1,1,n_bln),blowitup(n)
      real (kind=dp) :: fd_u(n),fd_v(n),fd_q(n),fd_t(n),fd_tke(n)

! Common blocks
      common /cdf_var/ id_rec,idvar(39),idfile,icount,jddim(4)
      integer :: id_rec, idvar, idfile, icount, jddim
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) ::  time
      integer :: lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      real (kind=dp) :: detw, deta, eta, etw

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy

      common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
      real (kind=dp) :: gm, gh, sm, sh, xl

      common /cb45/ u(n),v(n),w(n)
      real (kind=dp) :: u, v, w
      common /cb48/ sk,sl,dtrad(n),dtcon(n)
      real (kind=dp) ::  sk, sl, dtrad, dtcon

!     common /cb52/ ff(nkt,nka,n),fsum(n,0:nkc),nar(n) ! jjb wrong
      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)       ! jjb corrected, but mess up below, see comments
      real (kind=dp) :: ff, fsum
      integer :: nar

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      real(kind=dp) :: theta, thetl, t, talt, p, rho

      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
      real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a
      common /kurz/ fs1(nrlev),fs2(nrlev),totds(nrlev),ss(nrlev),
     &              fsn(nrlev),dtdts(nrlay)
      real (kind=dp) ::  fs1, fs2, totds, ss, fsn, dtdts

      common /lang/ fl1(nrlev),fl2(nrlev),fln(nrlev),dtdtl(nrlay)
      real (kind=dp) ::  fl1, fl2, fln, dtdtl

! == End of declarations =======================================================

      icount=icount+1

      idimcount(1)=1
      idimcount(2)=1
      idimcount(3)=1
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=icount

c time variables
      ifield(1,1,1)=lday
      k=nf_put_vara_int(idfile,idvar(1),idimstart,idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      ifield(1,1,1)=lst
      k=nf_put_vara_int(idfile,idvar(2),idimstart,idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      ifield(1,1,1)=lmin
      k=nf_put_vara_int(idfile,idvar(3),idimstart,idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      idimcount(1)=1
      idimcount(2)=1
      idimcount(3)=n_bln
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=icount
c wind field
      field(1,1,:)=u(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(4),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=v(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(5),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=w(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(6),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c temperature variables
      field(1,1,:)=theta(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(7),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=thetl(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(8),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=t(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(9),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c humidity variables
      field(1,1,:)=xm1(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(10),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=xm2(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(11),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=feu(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(12),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c other variables
      field(1,1,:)=p(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(13),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=rho(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(14),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c the next 3 variables are necessary to be f(t)/f(z) for comfortable plotting
c in Ferret: they have to have same dimensions as vars to be plotted
      field(1,1,:)=eta(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(15),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=etw(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(16),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      do kk=1,n_bln
         blowitup(kk)=lday*24.+lst+lmin/60.
      enddo
      field(1,1,:)=blowitup(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(17),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c heating rates
      field(1,1,:)=dtrad(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(18),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=dtcon(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(19),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=dtdts(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(20),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=dtdtl(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(21),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c turbulence
      field(1,1,:)=atke(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(22),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=atkh(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(23),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=atkm(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(24),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=tke(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(25),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=tkep(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(26),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=xl(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(27),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c flux divergences
c     calculate fd's
      do k=1,n-1
         ddz       = deta(k)*detw(k)
         fd_u(k)   = atkm(k)*(u(k+1)-u(k))/ddz
         fd_v(k)   = atkm(k)*(v(k+1)-v(k))/ddz
         fd_q(k)   = atkh(k)*(xm1(k+1)-xm1(k))/ddz
         fd_T(k)   = atkh(k)*(theta(k+1)-theta(k))/ddz
         fd_tke(k) = atke(k)*(tke(k+1)-tke(k))/ddz
      enddo
c     write fd's
      field(1,1,:)=fd_u(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(28),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=fd_v(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(29),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=fd_q(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(30),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=fd_T(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(31),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=fd_tke(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(32),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c particle sum
! jjb original version, which uses a wrong definition of fsum
! in any other places in the model, fsum is a one dimension array (fsum(n))
! see SR open_met for the explanation about the index 0:nkc
! the index 0 is the grand-total, and the indexes 1-4 for bins 1-4
! this could be done in the future (check Roland's todo list)
! anyway, at the moment, only fsum(n) is outputed (once) in the following change:

c$$$      field(1,1,:)=fsum(1:n_bln,0)
c$$$      k=nf_put_vara_double(idfile,idvar(33),idimstart,idimcount,field)
c$$$      if (k.ne.nf_noerr) call ehandle(k,fname)
      field(1,1,:)=fsum(1:n_bln)
      k=nf_put_vara_double(idfile,idvar(33),idimstart,idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
c$$$      field(1,1,:)=fsum(1:n_bln,1)
c$$$      k=nf_put_vara_double(idfile,idvar(34),idimstart,idimcount,field)
c$$$      if (k.ne.nf_noerr) call ehandle(k,fname)
c$$$      field(1,1,:)=fsum(1:n_bln,2)
c$$$      k=nf_put_vara_double(idfile,idvar(35),idimstart,idimcount,field)
c$$$      if (k.ne.nf_noerr) call ehandle(k,fname)
c$$$      field(1,1,:)=fsum(1:n_bln,3)
c$$$      k=nf_put_vara_double(idfile,idvar(36),idimstart,idimcount,field)
c$$$      if (k.ne.nf_noerr) call ehandle(k,fname)
c$$$      field(1,1,:)=fsum(1:n_bln,4)
c$$$      k=nf_put_vara_double(idfile,idvar(37),idimstart,idimcount,field)
c$$$      if (k.ne.nf_noerr) call ehandle(k,fname)

c lower and upper limit of cloud
      idimcount(3)=1
      ifield(1,1,1)=lcl
      k=nf_put_vara_int(idfile,idvar(38),idimstart,idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      ifield(1,1,1)=lct
      k=nf_put_vara_int(idfile,idvar(39),idimstart,idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      k=nf_sync(idfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine write_met

c
c----------------------------------------------------------------
c

      subroutine write_mic
! output of microphysics

! jjb work done
!     - corrected character length for fname to be consistent with SR ehandle
!     - missing declarations and implicit none
!     - introduce 'jpOutPart2dOpt' for easier tuning of data request

      USE config, ONLY :
! Imported Parameters:
     &     jpOutPart2dOpt,
! Imported Routines:
     &     abortM

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nka,
     &     nkt

      USE precision, ONLY :
! Imported Parameters:
     &     dp

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = "mic.nc"
! Local scalars:
      integer :: ia, ik, ind, jt, k
! Local arrays:
      integer ifield(1,1,1), idimcount(4), idimstart(4), indlist(nf/10)
      real (kind=dp) :: field(nka,nkt,nf/10)!, field2(2,nkt,nf),
      !&   field3(1,nkt,nf)

! Common blocks
      common /cdf_var_mic/ id_mic_rec,idvar_mic(6),idmicfile,
     & imiccount,jddim_mic(4)
      integer :: id_mic_rec, idvar_mic, idmicfile, imiccount, jddim_mic

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) :: time
      integer :: lday, lst, lmin, it, lcl, lct

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      real (kind=dp) :: ff, fsum
      integer :: nar

      common /oneDs/ partN(n,nkt,2),partr(n,nkt),drp(nkt)
      real (kind=dp) :: partN, partr, drp
      common /oneDsj/ rpw(nka), part1D(nka-1,nf)
      real (kind=dp) :: rpw, part1D

! == End of declarations =======================================================

      imiccount=imiccount+1

      idimcount(1)=1
      idimcount(2)=1
      idimcount(3)=1
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=imiccount

c time variables
      ifield(1,1,1)=lday
      k=nf_put_vara_int(idmicfile,idvar_mic(1),idimstart,
     &  idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      ifield(1,1,1)=lst
      k=nf_put_vara_int(idmicfile,idvar_mic(2),idimstart,
     & idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      ifield(1,1,1)=lmin
      k=nf_put_vara_int(idmicfile,idvar_mic(3),idimstart,
     &  idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      idimcount(1)=nka
      idimcount(2)=nkt
      idimcount(3)=nf/10

      select case (jpOutPart2dOpt)
      case (0)
         do ik=1,nf,10
            ind=ik/10 +1
            do ia=1,nka
               do jt=1,nkt
                  field(ia,jt,ind)=ff(jt,ia,ik)
               enddo
            enddo
         enddo

      case (1)
         data indlist /54,55,56,61,62,63,78,79,80,81/
         do ik=1,10
            ind=indlist(ik)
            do ia=1,nka
               do jt=1,nkt
                  field(ia,jt,ik)=ff(jt,ia,ind)
               enddo
            enddo
         enddo

      case default
         call abortM ('Error in SR write_mic: wrong choice for'//
     &        'jpOutPart2dOpt')
      end select

      k=nf_put_vara_double(idmicfile,idvar_mic(4),idimstart,
     & idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_sync(idmicfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      ! initialise to avoid being trapped in debuging option checking uninitialised arrays
!      field2(:,:,1) = 0._dp
!      field3(:,:,1) = 0._dp
!      do k=2,nf
!         do jt=1,nkt
!            field2(1,jt,k)=partN(k,jt,1)
!            field2(2,jt,k)=partN(k,jt,2)
!            field3(1,jt,k)=partr(k,jt)
!         enddo
!      enddo

!      idimcount(1)=2
!      idimcount(2)=nkt
!      idimcount(3)=nf

!      k=nf_put_vara_double(idmicfile,idvar_mic(5),idimstart,
!     & idimcount,field2)
!      if (k.ne.nf_noerr) call ehandle(k,fname)
!      k=nf_sync(idmicfile)
!      if (k.ne.nf_noerr) call ehandle(k,fname)

      idimcount(1)=nka-1
      idimcount(2)=nf
      idimcount(3)=1
      idimstart(3)=imiccount

       k=nf_put_vara_double(idmicfile,idvar_mic(6),idimstart(1:3),
     & idimcount(1:3),part1D)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_sync(idmicfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine write_mic


c
c----------------------------------------------------------------
c

      subroutine write_chem_gas (n_bln)

! jjb work done:
!     - character size consistent with ehandle subroutine
!     - unused declarations / variables / parameters removed
!     - use module

! 21-Sep-2020   Josue Bock   Minor bugfix: correct the third dimension of local array 'field': n_bln instead of n

      USE gas_common, ONLY :
! Imported Parameters:
     &     j1,
     &     j5,
! Imported Array Variables with intent (in):
     &     s1,
     &     s3

      USE cdf_var_gas, ONLY :
     &     idgasfile,
     &     igascount,
     &     idvar_gas

      USE precision, ONLY :
! Imported Parameters:
     &     dp

      implicit none

! Include statements:
      include 'netcdf.inc'

! Subroutine arguments
! Scalar arguments with intent(in):
      integer, intent(in) :: n_bln

! Local parameters:
      character (len=*), parameter :: fname = 'gas.nc'
      integer, parameter :: j0 = 3              ! number of time variables

! Local scalars:
      integer :: i0                             ! index offset
      integer :: jspec, k

! Local arrays:
      integer :: idimcount(4), idimstart(4), ifield(1,1,1)
      real (kind=dp) :: field(1,1,n_bln)

! Common blocks:
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) ::  time
      integer :: lday, lst, lmin, it, lcl, lct

!- End of header ---------------------------------------------------------------

      igascount=igascount+1

      idimcount(1)=1
      idimcount(2)=1
      idimcount(3)=1
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=igascount

c time variables
      ifield(1,1,1)=lday
      k=nf_put_vara_int(idgasfile,idvar_gas(1),idimstart,
     &      idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      ifield(1,1,1)=lst
      k=nf_put_vara_int(idgasfile,idvar_gas(2),idimstart,
     &  idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      ifield(1,1,1)=lmin
      k=nf_put_vara_int(idgasfile,idvar_gas(3),idimstart,
     &  idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c chemical species

      idimcount(1)=1
      idimcount(2)=1
      idimcount(3)=n_bln
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=igascount
c offset due to time variables
      i0 = j0

c write all defined species
      do jspec=1,j1
         field(1,1,:)=s1(jspec,1:n_bln)
         k=nf_put_vara_double(idgasfile,idvar_gas(i0+jspec),idimstart,
     &                        idimcount,field)
         if (k.ne.nf_noerr) call ehandle(k,fname)
      enddo

c radicals
! increase offset
      i0 = i0 + j1

      do jspec = 1,j5
         field(1,1,:)=s3(jspec,1:n_bln)
         k=nf_put_vara_double(idgasfile,idvar_gas(i0+jspec),idimstart,
     &                        idimcount,field)
         if (k.ne.nf_noerr) call ehandle(k,fname)
      enddo

      k=nf_sync(idgasfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine write_chem_gas


c
c----------------------------------------------------------------
c
      subroutine write_chem_aq (n_bln,halo,iod,nuc)

! 21-Sep-2020   Josue Bock   Minor bugfix: correct the third dimension of local array 'field': n_bln instead of n

      USE config, ONLY:
     &     nkc_l

      USE global_params, ONLY :
! Imported Parameters:
     &     j1=>j1_fake,
     &     j2,
     &     j6,
     &     n,
     &     nkc

      USE precision, ONLY :
! Imported Parameters:
     &     dp

      implicit none

! Include statements:
      include 'netcdf.inc'

      integer, intent(in) :: n_bln
      logical, intent(in) :: halo, iod, nuc

! Local parameters:
      character (len=*), parameter :: fname = "aq.nc"
      integer, parameter :: nliq=27, nhalo=10, niod=7
      integer, parameter :: nion=23, nionh=9,  nioni=5

! Local scalars:
      integer :: i0, i2, ispec, jspec
      integer :: ij, k, kc

      integer :: idimcount(4),idimstart(4),
     &        ifield(1,1,1),mliq(nliq),mhalo(nhalo),miod(niod),
     &        mion(nion),mionh(nionh),mioni(nioni)
      real (kind=dp) :: field(nkc_l,1,n_bln)

! Common blocks
      common /cdf_var_aq/ idaq_rec,idvar_aq(j2+j6+7),idaqfile,
     &  iliqcount,jddim_aq(4)
      integer :: idaq_rec, idvar_aq, idaqfile, iliqcount, jddim_aq

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) :: time
      integer :: lday, lst, lmin, it, lcl, lct

      common /blck11/ rc(nkc,n)
      real (kind=dp) :: rc
      common /blck12/ cw(nkc,n),cm(nkc,n)
      real (kind=dp) :: cw, cm
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      real (kind=dp) :: sl1, sion1

c add mercury/Hg

c this lists which indices are defined and can be output - this HAS to be the same as in
c SR open_chem_aq
c "normal" aqueous phase species plus "normal" uncharged additional species
      data mliq /1,2,3,4,5,6,7,15,16,17,19,20,21,35,36,37,60,70,78,88,
     &     89,95,0,0,0,0,0/
c no arithmetic expressions allowed within data statement!!

      mliq(23) = j1+2
      mliq(24) = j1+3
      mliq(25) = j1+4
      mliq(26) = j1+6
      mliq(27) = j1+11

c normal not defined/used: 1,8,9,10,11,12,13,14,16,18,22,23,24,25,26,27,28,29,31,32,33,34,37,61,62,63,64,66,67,68,69,j1+4,j1+5
c Cl+Br species and radicals! HCl is on purpose not defined as halogen

      data mhalo /30,38,41,42,43,46,47,94,0,0/
      mhalo(9) = j1+8
      mhalo(10) = j1+9

c Cl+Br not defined/used: 39,40,44,45,65
c I species and radicals

      data miod /49,53,54,55,0,0,0/
      miod(5) = j1+5
      miod(6) = j1+7
      miod(7) = j1+12
c I not defined/used: 48,50,51,52,56,57,58,59,71,72,73
c ions
      data mion /1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,19,20,21,30,31,
     &     32,33,39/
c ions not used:17,18,40 !Cl- is on purpose not defined as halogen
c Cl and Br ions
      data mionh /15,22,23,24,25,26,27,28,29/
c I ions
      data mioni /34,35,36,37,38/


      iliqcount=iliqcount+1

      idimcount(1)=1
      idimcount(2)=1
      idimcount(3)=1
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=iliqcount

c time variables
      ifield(1,1,1)=lday
c      k=nf_put_var1_int(idfile,idvar(1),icount,lday)
      k=nf_put_vara_int(idaqfile,idvar_aq(1),idimstart,
     &      idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      ifield(1,1,1)=lst
c      k=nf_put_var1_int(idfile,idvar(2),icount,lst)
      k=nf_put_vara_int(idaqfile,idvar_aq(2),idimstart,
     &  idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      ifield(1,1,1)=lmin
c      k=nf_put_var1_int(idfile,idvar(3),icount,lmin)
      k=nf_put_vara_int(idaqfile,idvar_aq(3),idimstart,
     &  idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c chemical species

      idimcount(1)=nkc_l
      idimcount(2)=1
      idimcount(3)=n_bln
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=iliqcount

c here loop for "normal" species
      i0=3
      do ispec=1,nliq
         jspec=mliq(ispec)
         do k=1,n_bln
            do ij=1,nkc_l
               field(ij,1,k)=sl1(jspec,ij,k)
            enddo
         end do
         k=nf_put_vara_double(idaqfile,idvar_aq(jspec+i0),idimstart,
     &        idimcount,field)
         if (k.ne.nf_noerr) call ehandle(k,fname)
      enddo

      if (halo) then
         do ispec=1,nhalo
            jspec=mhalo(ispec)
            do k=1,n_bln
               do ij=1,nkc_l
                  field(ij,1,k)=sl1(jspec,ij,k)
               enddo
            end do
            k=nf_put_vara_double(idaqfile,idvar_aq(jspec+i0),
     &           idimstart,idimcount,field)
            if (k.ne.nf_noerr) call ehandle(k,fname)
         enddo

         if (iod) then
            do ispec=1,niod
               jspec=miod(ispec)
               do k=1,n_bln
                  do ij=1,nkc_l
                     field(ij,1,k)=sl1(jspec,ij,k)
                  enddo
               end do
               k=nf_put_vara_double(idaqfile,idvar_aq(jspec+i0),
     &              idimstart,idimcount,field)
               if (k.ne.nf_noerr) call ehandle(k,fname)
            enddo
         endif
      endif

! jjb uncommented in v733, commented in v741. Should be deleted anyway?
c      if (nuc) then !NUCV
c         jspec=75
c         do k=1,n_bln
c            do ij=1,nkc_l
c               field(ij,1,k)=sl1(jspec,ij,k)
c            end do
c         enddo
c         k=nf_put_vara_double(idaqfile,idvar_aq(jspec+i0),
c     &        idimstart,idimcount,field)
c         if (k.ne.nf_noerr) call ehandle(k,fname)
c      endif

c start of ions
      i2=i0+j2
      do ispec=1,nion
         jspec=mion(ispec)
         do k=1,n_bln
            do ij=1,nkc_l
               field(ij,1,k)=sion1(jspec,ij,k)
            enddo
         end do
         k=nf_put_vara_double(idaqfile,idvar_aq(jspec+i2),idimstart,
     &        idimcount,field)
         if (k.ne.nf_noerr) call ehandle(k,fname)
      enddo

      if (halo) then
         do ispec=1,nionh
            jspec=mionh(ispec)
            do k=1,n_bln
               do ij=1,nkc_l
                  field(ij,1,k)=sion1(jspec,ij,k)
               enddo
            end do
            k=nf_put_vara_double(idaqfile,idvar_aq(jspec+i2),
     &           idimstart,idimcount,field)
            if (k.ne.nf_noerr) call ehandle(k,fname)
         enddo
         if (iod) then
            do ispec=1,nioni
               jspec=mioni(ispec)
               do k=1,n_bln
                  do ij=1,nkc_l
                     field(ij,1,k)=sion1(jspec,ij,k)
                  enddo
               end do
               k=nf_put_vara_double(idaqfile,idvar_aq(jspec+i2),
     &              idimstart,idimcount,field)
               if (k.ne.nf_noerr) call ehandle(k,fname)
            enddo
         endif
      endif

c LWC
      do k=1,n_bln
         do kc=1,nkc_l
            field(kc,1,k)=cw(kc,k)
         end do
      end do
      k=nf_put_vara_double(idaqfile,idvar_aq(i2+41),idimstart,
     &     idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
c radius
      do k=1,n_bln
         do kc=1,nkc_l
            field(kc,1,k)=rc(kc,k)
         end do
      end do
      k=nf_put_vara_double(idaqfile,idvar_aq(i2+42),idimstart,
     &     idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)



      k=nf_sync(idaqfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine write_chem_aq

c
c----------------------------------------------------------------
c

      subroutine write_jrate (n_bln,chamber)

! 21-Sep-2020   Josue Bock   Minor bugfix: correct the third dimension of local array 'field': n_bln instead of n

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nphrxn

      USE precision, ONLY :
! Imported Parameters:
     &     dp

      implicit none

! Include statements:
      include 'netcdf.inc'

      integer, intent(in) :: n_bln
      logical, intent(in) :: chamber
! Local parameters:
      character (len=*), parameter :: fname = "jrate.nc"
! Local scalars:
      integer :: ispec, k
! Local arrays:
      integer :: idimcount(4),idimstart(4),ifield(1,1,1)
      real (kind=dp) :: field(1,1,n_bln)

! Common blocks:
      common /cdf_var_jrat/ idjrat_rec,idvar_jrat(50),idjratfile,
     &     ijratcount,jddim_jrat(4)
      integer :: idjrat_rec,idvar_jrat,idjratfile,ijratcount,jddim_jrat
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) :: time
      integer :: lday, lst, lmin, it, lcl, lct

      common /band_rat/ photol_j(nphrxn,n)
      real (kind=dp) :: photol_j
      common /chamber_ph_r/ ph_rat_chamber(nphrxn)
      real (kind=dp) :: ph_rat_chamber

      ijratcount=ijratcount+1

      idimcount(1)=1
      idimcount(2)=1
      idimcount(3)=1
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=ijratcount

c time variables
      ifield(1,1,1)=lday
      k=nf_put_vara_int(idjratfile,idvar_jrat(1),idimstart,
     &      idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      ifield(1,1,1)=lst
      k=nf_put_vara_int(idjratfile,idvar_jrat(2),idimstart,
     &  idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      ifield(1,1,1)=lmin
      k=nf_put_vara_int(idjratfile,idvar_jrat(3),idimstart,
     &  idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c chemical species

      idimcount(1)=1
      idimcount(2)=1
      if (.not.chamber) then
         idimcount(3)=n_bln
      else
         idimcount(3)=1
      end if
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=ijratcount

      if (.not.chamber) then
         do ispec=1,nphrxn
            field(1,1,:)=photol_j(ispec,1:n_bln)
            k=nf_put_vara_double(idjratfile,idvar_jrat(ispec+3),
     &           idimstart,idimcount,field)
            if (k.ne.nf_noerr) call ehandle(k,fname)
         enddo
      else
         do ispec=1,nphrxn
            field(1,1,1)=ph_rat_chamber(ispec)
            k=nf_put_vara_double(idjratfile,idvar_jrat(ispec+3),
     &           idimstart,idimcount,field(1:1,1:1,1:1))
            if (k.ne.nf_noerr) call ehandle(k,fname)
         enddo
      end if

      k=nf_sync(idjratfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine write_jrate
c
c----------------------------------------------------------------
c


      subroutine write_rxn

      USE global_params, ONLY :
! Imported Parameters:
     &     nlev,
     &     nrxn

      USE precision, ONLY :
! Imported Parameters:
     &     dp

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = "rxnrate.nc"
! Local scalars:
      integer :: ilev, irxn, k
! Local arrays:
      integer :: ifield(1,1,1),idimcount(4),idimstart(4)
      real (kind=dp) :: field(1,nrxn,nlev)
! Common blocks:
      common /cdf_var_rxn/ idrxn_rec,idvar_rxn(4),idrxnfile,
     &   irxncount,jddim_rxn(4)
      integer :: idrxn_rec,idvar_rxn,idrxnfile,irxncount,jddim_rxn
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) :: time
      integer :: lday, lst, lmin, it, lcl, lct

      common /budg/ bg(2,nrxn,nlev),il(nlev)
      real (kind=dp) :: bg
      integer :: il

      irxncount=irxncount+1

      idimcount(1)=1
      idimcount(2)=1
      idimcount(3)=1
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=irxncount

c time variables
      ifield(1,1,1)=lday
      k=nf_put_vara_int(idrxnfile,idvar_rxn(1),idimstart,
     &      idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      ifield(1,1,1)=lst
      k=nf_put_vara_int(idrxnfile,idvar_rxn(2),idimstart,
     &  idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      ifield(1,1,1)=lmin
      k=nf_put_vara_int(idrxnfile,idvar_rxn(3),idimstart,
     &  idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c reaction rates, written into one variable as FERRET has a max number of 2000 variables

      idimcount(1)=1
      idimcount(2)=nrxn
      idimcount(3)=nlev

      do ilev=1,nlev
         do irxn=1,nrxn
            field(1,irxn,ilev)=bg(1,irxn,ilev)
         enddo
      enddo

      k=nf_put_vara_double(idrxnfile,idvar_rxn(4),idimstart,
     &     idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      k=nf_sync(idrxnfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine write_rxn

c
c----------------------------------------------------------------
c
c
      subroutine write_nuc
c output of nucleation parameters

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nkt

      USE precision, ONLY :
! Imported Parameters:
     &     dp

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = "nuc.nc"
! Local scalars:
      integer :: k
! Local arrays:
      integer :: ifield(1,1,1), idimcount(4), idimstart(4)
      real (kind=dp) :: field(1,1,n), field2(1,nkt,n)

! Common blocks:
      common /cdf_var_nuc/ idnuc_rec,idvar_nuc(23),idnucfile,
     &   inuccount,jddim_nuc(4)
      integer :: idnuc_rec,idvar_nuc,idnucfile,inuccount,jddim_nuc

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      real (kind=dp) :: time
      integer :: lday, lst, lmin, it, lcl, lct

      common /nucl/ xn_new(n), xn_acc(n), xv_acc(n), dh2so4(n), dnh3(n)
      real (kind=dp) :: xn_new, xn_acc, xv_acc, dh2so4, dnh3

      common /nuclapp/ xn_app(n), xn_apacc(n), xv_apacc(n),bn_ges(n),
     &                 bd_mean(n), dnucv(n), grorate(n), concnuc(n)
      real (kind=dp) :: xn_app, xn_apacc, xv_apacc, bn_ges,
     &                 bd_mean, dnucv, grorate, concnuc

      common /nuclio/ xn_newio(n), xn_accio(n), xv_accio(n), doio(n)
      real (kind=dp) :: xn_newio, xn_accio, xv_accio, doio

      common /backpart/ partd(nkt,n), partNu(nkt,n), partsa(n)
      real (kind=dp) :: partd, partNu, partsa

! == End of declarations =======================================================

      inuccount=inuccount+1

      idimcount(1)=1
      idimcount(2)=1
      idimcount(3)=1
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=inuccount

c time variables
      ifield(1,1,1)=lday
      k=nf_put_vara_int(idnucfile,idvar_nuc(1),idimstart,
     &       idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      ifield(1,1,1)=lst
      k=nf_put_vara_int(idnucfile,idvar_nuc(2),idimstart,
     &        idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)
      ifield(1,1,1)=lmin
      k=nf_put_vara_int(idnucfile,idvar_nuc(3),idimstart,
     &        idimcount,ifield)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c nucleation parameters

      idimcount(1)=1
      idimcount(2)=1
      idimcount(3)=n
      idimcount(4)=1

      idimstart(1)=1
      idimstart(2)=1
      idimstart(3)=1
      idimstart(4)=inuccount
c
      field(1,1,:)=xn_new(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(4),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=xn_acc(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(5),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=xv_acc(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(6),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=xn_newio(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(7),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=xn_accio(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(8),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=xv_accio(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(9),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=xn_app(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(10),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=xn_apacc(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(11),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=xv_apacc(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(12),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=bn_ges(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(13),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=bd_mean(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(14),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=dnh3(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(15),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=dh2so4(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(16),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=doio(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(17),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=dnucv(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(18),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=grorate(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(19),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=concnuc(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(20),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field(1,1,:)=partsa(1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(21),idimstart,
     &  idimcount,field)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      idimcount(1)=1
      idimcount(2)=nkt
      idimcount(3)=n
      idimcount(4)=1

      field2(1,1:nkt,1:n)=partd(1:nkt,1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(22),idimstart,
     &  idimcount,field2)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      field2(1,1:nkt,1:n)=partNu(1:nkt,1:n)
      k=nf_put_vara_double(idnucfile,idvar_nuc(23),idimstart,
     &  idimcount,field2)
      if (k.ne.nf_noerr) call ehandle(k,fname)

c

      k=nf_sync(idnucfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine write_nuc

c
c----------------------------------------------------------------
c
      subroutine close_met

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = "met.nc"
! Local scalars:
      integer :: k
      common /cdf_var/ id_rec,idvar(39),idfile,icount,jddim(4)
      integer :: id_rec,idvar,idfile,icount,jddim

      k=nf_close(idfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine close_met
c
c----------------------------------------------------------------
c

      subroutine close_mic

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = "mic.nc"
! Local scalars:
      integer :: k
      common /cdf_var_mic/ id_mic_rec,idvar_mic(6),idmicfile,
     & imiccount,jddim_mic(4)
      integer :: id_mic_rec, idvar_mic, idmicfile, imiccount, jddim_mic

      k=nf_close(idmicfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine close_mic

c
c----------------------------------------------------------------
c

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine close_chem_gas

      USE cdf_var_gas, ONLY :
     &     idgasfile,
     &     idvar_gas

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = 'gas.nc'
! Local scalars:
      integer :: k

!- End of header ---------------------------------------------------------------

      deallocate ( idvar_gas )

      k=nf_close(idgasfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine close_chem_gas
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine close_chem_aq

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = "aq.nc"
! Local scalars:
      integer :: k
      common /cdf_var_aq/ idaq_rec,idvar_aq(j2+j6+7),idaqfile,
     &  iliqcount,jddim_aq(4)
      integer :: idaq_rec,idvar_aq,idaqfile,iliqcount,jddim_aq

      k=nf_close(idaqfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine close_chem_aq
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine close_jrate

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname  = "jrate.nc"
! Local scalars:
      integer :: k
      common /cdf_var_jrat/ idjrat_rec,idvar_jrat(50),idjratfile,
     & ijratcount,jddim_jrat(4)
      integer :: idjrat_rec,idvar_jrat,idjratfile,ijratcount,jddim_jrat

      k=nf_close(idjratfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine close_jrate

c----------------------------------------------------------------
c

      subroutine close_rxn

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = "rxnrate.nc"
! Local scalars:
      integer :: k
      common /cdf_var_rxn/ idrxn_rec,idvar_rxn(4),idrxnfile,
     &   irxncount,jddim_rxn(4)
      integer :: idrxn_rec,idvar_rxn,idrxnfile,irxncount,jddim_rxn

      k=nf_close(idrxnfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine close_rxn

c
c----------------------------------------------------------------
c
      subroutine close_nuc

      implicit none

! Include statements:
      include 'netcdf.inc'

! Local parameters:
      character (len=*), parameter :: fname = "nuc.nc"
! Local scalars:
      integer :: k
      common /cdf_var_nuc/ idnuc_rec,idvar_nuc(23),idnucfile,
     &     inuccount,jddim_nuc(4)
      integer :: idnuc_rec,idvar_nuc,idnucfile,inuccount,jddim_nuc

      k=nf_close(idnucfile)
      if (k.ne.nf_noerr) call ehandle(k,fname)

      end subroutine close_nuc

c
c----------------------------------------------------------------
c

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ehandle (nb_err,file)

! the routine handle_err makes use of the function nf_STRERROR:
!
!     The function nf_STRERROR returns a static reference to an error
!     message string corresponding to an integer netCDF error status or
!     to a system error number, presumably returned by a previous call
!     to some other netCDF function. The list of netCDF error status
!     codes is available in the appropriate include file for each
!     language binding.

! jjb: calling abortM here instead of stop would lead to a potential recursive call

      USE file_unit, ONLY :
     &     jpfunerr ! standard error files

      implicit none

! Include statements:
      include 'netcdf.inc'

      integer, intent(in) :: nb_err
      character (len=*), intent(in) :: file

      write(jpfunerr,*) 'netCDF I/O error with file: ',file
      write(jpfunerr,*) 'error number: ',nb_err
      write(jpfunerr,*) nf_strerror(nb_err)
      stop 'Stopped by SR ehandle'

      end subroutine ehandle
