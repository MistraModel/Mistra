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


! MISTRA with chemistry
!
! str.f : main, meteo-init, microphysics, turbulence



! Author:
! ------
  !    Andreas Bott


! box = .true. is option to run model in one level only
!     so far without dynamics, microphysics and radiation.
!     Only the photolysis rates are updated and T, rh, ..
!     are changed with sinus function.
!     Init is the same as in 1D run, to make sure that all
!     variables are defined and to avoid too much differences
!     between the box and the 1D version.

!
! INITIALISATION
!
! ... no restart
!       |
!       |_____SR equil (case 0)
!       |      |_____FN rgl

! TIME LOOP
! ...   |_____SR difm
!              |_____SR atk1
!      ... if (chem)
!       |_____SR difc
!       |
!     if (mic)
!       |_____SR difp
!       |_____SR kon
!       |      |    _SR equil (case 1)   ! if rH < 70%
!       |      |   /  |_____FN rgl
!       |      |__/
!       |      |  \
!       |      |   \_SR subkon           ! if rH >= 70%
!       |      |      |_____SR advec
!       |      |
!       |     ... if (chem)
!       |      |_____SR konc
!       |
!       |_____SR sedp
!       |_____SR equil (case 2: levels nf+1 -> n)
!
!

! == End of header =============================================================

program mistra

  USE config, ONLY : &
! External subroutine
       read_config,  &
       abortM,       &
! Config switches
       binout,       &
       BL_box,       &
       box,          &
       chamber,      &
       nlevbox,      &
       z_box,        &
       chem,         &
       halo,         &
       iod,          &
       isurf,        &
       mic,          &
       netCDF,       &
       nuc, Napari, Lovejoy, &
       rst,          &
       lstmax,       &
       lpBuys13_0D,  &
       lpJoyce14bc

  USE file_unit, ONLY : &
       jpfunae, jpfunout

  USE global_params, ONLY : &
! Imported Parameters:
       nf,                  &
       n,                   &
       nm,                  &
       nrlay,               &
       nka,                 &
       nkt,                 &
       nphrxn,              &
       mbs

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  interface
     subroutine equil (ncase,kk)
       integer,           intent(in)  :: ncase
       integer, optional, intent(in)  :: kk    ! model level for calculation
     end subroutine equil
  end interface

  character (len=1), parameter :: fogtype = 'a'  ! kept for historical reasons, appendix to many file names
  real (kind=dp), parameter :: dt = 60._dp       ! timestep for slowly varying processes
  real (kind=dp), parameter :: dd = 10._dp       ! fractional timestep for faster processes
  character (len=10) :: fname
  integer :: ilmin, it0, itmax, ia, ij, jt, k
  integer :: n_bl, n_bln, n_bl8, nz_box
  logical :: llboxch, llnucboth
  logical :: llinit, llcallphotol, llsetjrates0
  real (kind=dp) :: atmax, tkemax, xm2max, xra
  real (kind=dp) :: box_switch, u0min
  real (kind=dp) :: aer(n,nka)

! Common blocks:
  common /cb16/ u0,albedo(mbs),thk(nrlay)
  real (kind=dp) :: u0, albedo, thk

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb48/ sk,sl,dtrad(n),dtcon(n)
  real (kind=dp) :: sk, sl, dtrad, dtcon

  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a

  common /band_rat/ photol_j(nphrxn,n)
  real (kind=dp) :: photol_j


  ! initialisation switch
  llinit = .true.

  ! Read namelist
  call read_config

  if (box) then
     write(jpfunout,*)'box model run'
  elseif (chamber) then
     write(jpfunout,*)'chamber model run'
  else
     write(jpfunout,*)'1-D model run'
  endif
  llboxch = box.or.chamber

! it's important to keep this n_bl = 2 for box runs as loops are designed that way
! (especially output)
  if (box) then
     n_bl  = 2
     n_bln = 2
     n_bl8 = 1
  else if (chamber) then
     n_bl  = 2 ! 3
     n_bln = 2 ! 3
     n_bl8 = 1
     xra = 0._dp ! initialise since partdep will not be called
  else
     n_bl  = nf
     n_bln = n
     n_bl8 = 15
  endif
! numerical gridpoints
  call grid
! open input/output files
  call openm (fogtype)
  if (chem) call openc (fogtype)

  ! Initialise the Mistra-KPP interface to handle the species
  if (chem) call mk_interface

  ! Initialise the nucleation module: search for relevant species in user list
  if (chem .and. nuc) then
     call nuc_init(Napari,Lovejoy,iod)
     llnucboth = Napari.and.Lovejoy
  end if

! netCDF output
  if (netCDF) call open_netcdf(n_bln,chem,mic,halo,iod,box,chamber,nuc)

  if (box.or.chamber) call get_n_box (z_box,nz_box)
  call write_grid ! writes information on grid that is not f(t)

  if (.not.rst) then

! Continue the initialisation, no-restart case
! --------------------------------------------
! initial meteorological and chemical input
     call initm (fogtype,rst)
     if (chem) call initc(n_bl)

! number of iterations
     it0=0
     itmax=60*lstmax
! initial exchange coefficients and turbulent kinetic energy
     call atk0

! initial position of humidified aerosols on Köhler curve
     call equil(0)

  else

! Continue the initialisation, restart case
! -----------------------------------------
! read meteorological and chemical input from output of previous run
     call startm (fogtype)
! init some microphysical data that's not in SR startm
     call initm (fogtype,rst)
!+      it0=it    ! use when time stamp from restart run is to be preserved
     it0=0
     it=0
     if (chem) call startc (fogtype)

! number of iterations
     itmax=it0+60*lstmax

! jjb under development, commented on purpose
!! allocate arrays and initialise vmean
!     if (chem) call v_mean_init
!     if (chem) call v_mean (t(:nmax_chem_aer))

  end if
! Continue the initialisation, both cases
! ---------------------------------------

! initialization of radiation code, and first calculation
  call radiation (llinit)

! 1D distribution of particles at model start
  call oneD_dist_jjb

! output of meteorological and chemical constants of current run
  call constm
  if (chem) call constc
! output of initial vertical profiles
  call profm (dt)
  if (chem) call profc (dt,mic)

! initial photolysis rates
!  if (chem) call photol
  if (chem) then
     call photol_initialize
     call photol
  end if

! initial output for plotting
  if (binout) then
     call ploutm (fogtype,n_bln)
     if (mic.and..not.llboxch) call ploutp (fogtype)
     call ploutr (fogtype,n_bln)
     call ploutt (fogtype,n_bln)
     if (chem) then
        call ploutc (fogtype,mic,n_bl,n_bl8)
        call ploutj (fogtype,n_bln)
     end if
  endif

  if (chem) call out_mass
  if (netCDF) call write_netcdf(n_bln,chem,mic,halo,iod,box,chamber,nuc)

  time = dt * real(it0,dp)
! local time: day (lday), hours (lst), minutes (lmin)
  fname='tim .out'
  fname(4:4)=fogtype
  open (99, file=fname,status='unknown',err=2005)
  atmax=0._dp
  write (99,6000) lday,lst,lmin,atmax
  close (99)
2005 continue

  if (box) call box_init (nlevbox,nz_box,n_bl,BL_box)
  if (box) box_switch=1._dp
  if (chamber) call chamb_init (n_bl)
  if (chamber) box_switch=0._dp

  ! initialisation switch
  llinit = .false.

  ! Define minimum solar angle for photolysis rates calculation
  if (lpBuys13_0D) then
     u0min = 1.75e-2_dp
  else
     u0min = 3.48e-2_dp
  end if

  write(jpfunout,*)'end initialisation str.f'

! ====================integration in time=====================
! outer time loop: minutes
  do it=it0+1,itmax
     if (lct.gt.nf) stop 'lct.gt.nf'
!         time=time+dt

     lmin = lmin + 1
     if (lmin.eq.60) then
        lmin = lmin - 60
        lst  = lst + 1
        if (lst.eq.24) then
           lst  = 0
           lday = lday + 1
        endif
     endif

! dry dep velocities
     if (.not.chamber) call partdep (xra)

! inner time loop: 10 sec
     do ij=1,6
        time=time+dd
! --------1D only start------------
! skip dynamics, microphysics and radiation for box and chamber model run
        if (.not.llboxch) then

! if w-field variable in time call wfield
!   WARNING: check this routine before using, it is not the same as latest version used by Bott
!         call wfield

! turbulent exchange of thermodynamic variables
           call difm (dd)
           if (chem) call difc (dd)   ! turbulent exchange of chemical species

! microphysics
           if (mic) then
              call difp (dd)          ! turbulent exchange of particles
! condensation/evaporation, update of chemical concentrations
              call kon (dd,chem)
! gravitational settling of particles
              call sedp (dd)
! put aerosol into equilibrium with current rel hum for k>nf
              call equil (2)

           else
! put aerosol into equilibrium with current rel hum
!  jjb: this seems strange, why only run for k=n_bl?
              call equil (1,n_bl)
           endif

! radiative heating
           do k=2,nm
              t(k) = t(k) + dtrad(k) * dd
           enddo

! Surface routines: flux balances at the earth's surface in SR surf*
           select case (isurf)
! water surface
           case (0)
              call surf0 (dd)
           case (1)
! temperature and humidity within the soil
              call soil (dd)
              call surf1 (dd)
           case default
              call abortM ('wrong choice for isurf, must be 0 (water) or 1 (soil), change namelist')
           end select

! dry deposition and emission of chemical species
           if (chem) then
              call sedc (dd)
! wet deposition of chemical species
              call sedl (dd)
! chemical reactions
              call stem_kpp (dd,xra,z_box,n_bl,box,chamber,nuc)
              if (nuc) then
                 if (llnucboth) then
                    call appnucl2 (dd,llnucboth)
                 else
                    call appnucl (dd,Napari,Lovejoy,llnucboth)
                 endif
!                   --- de-comment if .asc output for gnu-plotting is desired ---
!                    call nucout1
              endif
           endif
        else                ! if .not.box .not.chamber
! --------1D only end------------
! --------box model version only start -------------------
           if (box) then
! put aerosol into equilibrium with current rel hum
!               if (.not.mic) call equil (1,n_bl)
! call u0, T, rh, J-values  .. update
           !call box_update(box_switch,ij,nlevbox,nz_box,n_bl,chem,halo,iod,BL_box) ! jjb 3 unused arguments
              call box_update(box_switch,ij,nlevbox,nz_box,n_bl,BL_box)
! gas phase emissions and deposition
              call sedc_box (dd,z_box,n_bl)
! particle and aqueous phase deposition
              call box_partdep (dd,z_box,n_bl)
! aerosol emission and chemical reactions
              if (chem) then
                 call stem_kpp (dd,xra,z_box,n_bl,box,chamber,nuc)
              endif

           else if (chamber) then
! --------chamber model version only start -------------------
! call u0, T, rh, J-values  .. update
              call chamb_update (n_bl,ij)
! aerosol emission and chemical reactions
              if (chem) then
                 call stem_kpp (dd,xra,z_box,n_bl,box,chamber,nuc)
              endif
! --------chamber model version only end -------------------
           end if

        endif               ! if .not.box .not.chamber
! --------box model version only end -------------------
     enddo                  ! ij-loop : end of fractional timestep loop

! radiative fluxes and heating rates
     if (.not.llboxch) call radiation (llinit)
! new photolysis rates
! for chamber mode photolysis rates are read by SR photol_chamber
! from file chamber.dat (see SR chamb_init and chamb_update)
     if (chem.and..not.chamber) then

        ! Condition to call photolysis code (may depend on configuration)
        if (lpJoyce14bc) then
           if (u0.gt.1.0d-2) then
              llcallphotol = .true.
              llsetjrates0 = .false.
           else
              llcallphotol = .false.
              llsetjrates0 = .true.
           end if
        else
           if (u0.gt.u0min .and. lmin/2*2.eq.lmin) then
              llcallphotol = .true.
              llsetjrates0 = .false.
           else
              llcallphotol = .false.
              if (u0.gt.u0min) then ! in this case, keep using the previously calculated rates
                 llsetjrates0 = .false.
              else
                 llsetjrates0 = .true.
              end if
           end if
        end if
        ! Call photolysis code
        if (llcallphotol) then
           call photol
           if (box.and.BL_box) call ave_j (nz_box,n_bl)
        else if (llsetjrates0) then
           photol_j(:,:) = 0._dp
        endif
     endif

! output of meteorological and chemical variables ----------------------
     ilmin=15
     if (chamber) ilmin=1 !output every minute
     if (lmin/ilmin*ilmin.eq.lmin) then
! calc 1D size distribution for output
        !call oneD_dist_new
        call oneD_dist_jjb
! binary output
        if (binout) then
           call ploutm (fogtype,n_bln)
           if (lmin/30*30.eq.lmin.and.mic.and..not.llboxch) call ploutp (fogtype)
           call ploutr (fogtype,n_bln)
           call ploutt (fogtype,n_bln)
           if (chem) call ploutc (fogtype,mic,n_bl,n_bl8)
!         if (chem.and.lmin/60*60.eq.lmin) call ploutj(fogtype,n_bln)
        endif

! netCDF output
        if (netCDF) call write_netcdf(n_bln,chem,mic,halo,iod,box,chamber,nuc)
! output of data from nucleation
        if (chem.and.nuc) call nucout2
! output from mass balance
        if (chem) call out_mass
     endif

! hourly output of profiles in ascii files
     if (lmin/60*60.eq.lmin) then
        call profm (dt)
!       call profr
        if (chem) call profc (dt,mic)
     endif
! output for restart option
!     comment these calls to save disk space for production runs:
     if (lst/12*12.eq.lst.and..not.box.and.lmin.eq.0) then
        call outm
        if (chem) call outc
     endif

! output of "tima.out"
     atmax  = 0._dp
     tkemax = 0._dp
     xm2max = 0._dp
     do k=lcl,nf
        atmax  = max(atmax, atkh(k))
        tkemax = max(tkemax, tke(k))
        xm2max = max(xm2max, xm2(k) * 1000._dp / rho(k))
     enddo
     open (99, file=fname,status='unknown',err=1000)
     write (99,6010) lday,lst,lmin,tkemax,atmax,xm2max,eta(lcl),eta(lct)
     write (jpfunout,6010) lday,lst,lmin,tkemax,atmax,xm2max,eta(lcl),eta(lct)
6000 format (' time: ',i2,':',i2,':',i2,3x,' iteration: ',f10.3,3x,'cloudy region: ',f7.1,' - ',f7.1)
6010 format (1x,i2,':',i2,':',i2,3f10.3,3x,'cloudy region: ',f7.1,' - ',f7.1)
     close (99)
1000 continue

  end do
! =========================end of time integration=====================


! final output of restart files
  call outm
  if (chem) call outc
! final output of aerosol size distribution
  do k=1,n
     do ia=1,nka
        aer(k,ia)=0._dp
        do jt=1,nkt
           aer(k,ia)=aer(k,ia)+ff(jt,ia,k)
        enddo
     enddo
  enddo
  fname='ae .out'
  fname(3:3)=fogtype
  open (jpfunae, file=fname,status='unknown',form='unformatted')
  write (jpfunae) aer
  close (jpfunae)

  if (netCDF) call close_netcdf(mic,chem,box,nuc)

  stop 'main program'
end program mistra

!
!-----------------------------------------------------------------------
!

subroutine openm (fogtype)
! input/output files


! Author:
! ------
  !    Bott and RvG?


! Modifications :
! -------------
  !

! == End of header =============================================================


  USE config, ONLY : &
       binout, &
       cinpdir,&
       coutdir, &
       rst, mic

  USE file_unit, ONLY : &
       jpfunclarke, &
       jpfunprofm, jpfunprofr, &
       jpfunpm, jpfunpb, &        ! ploutm files: pm*, pb*
       jpfunpr, &                 ! ploutr file: pr*
       jpfunpt, &                 ! ploutt file: pt*
       jpfunf1, jpfunf2, jpfunf3  ! ploutp files: f1*, f2*, f3*

  USE data_surface, ONLY : fu, ft, xzpdl, xzpdz0

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  character (len=1), intent(in) :: fogtype

! Local scalars:
  character (len=10) :: fname    ! I/O files names
  character (len=150) :: clpath  ! complete path to file
  integer :: i,k                 ! implied do loops indexes

!- End of header ---------------------------------------------------------------

! input Clarke-tables
  open (jpfunclarke, file=trim(cinpdir)//'clarke.dat', status='old')
  read (jpfunclarke,5000) ((fu(i,k),i=1,9),(fu(i,k),i=10,18),k=1,7)
  read (jpfunclarke,5000) ((ft(i,k),i=1,9),(ft(i,k),i=10,18),k=1,7)
  read (jpfunclarke,5010) (xzpdl(i),i=1,9),(xzpdl(i),i=10,18)
  read (jpfunclarke,5020) (xzpdz0(i),i=1,7)
  close (jpfunclarke)

5000 format (9f8.4)
5010 format (9f6.2)
5020 format (7f5.0)

! output vertical profiles of meteorological variables
  fname='profm .out'
  fname(6:6)=fogtype
  open (jpfunprofm, file=fname,status='unknown')

! Create output files by opening them once, only if no restart
  if (.not.rst) then
     ! plout* output
     if (binout) then
        fname='pm .out'
        fname(3:3)=fogtype
        clpath=trim(coutdir)//trim(fname)
        open (jpfunpm, file=trim(clpath), status='new',form='unformatted')
        close (jpfunpm)

        fname(2:2)='b'
        clpath=trim(coutdir)//trim(fname)
        open (jpfunpb, file=trim(clpath), status='new',form='unformatted')
        !close (jpfunpb) ! closed by initm

        fname(2:2)='r'
        clpath=trim(coutdir)//trim(fname)
        open (jpfunpr, file=trim(clpath), status='new',form='unformatted')
        close (jpfunpr)

        fname(2:2)='t'
        clpath=trim(coutdir)//trim(fname)
        open (jpfunpt, file=trim(clpath), status='new',form='unformatted')
        close (jpfunpt)

        if (mic) then
           fname='f1 .out'
           fname(3:3)=fogtype
           clpath=trim(coutdir)//trim(fname)
           open (jpfunf1, file=trim(clpath), status='new',form='unformatted')
           close (jpfunf1)

           fname(2:2)='2'
           fname(3:3)=fogtype
           clpath=trim(coutdir)//trim(fname)
           open (jpfunf2, file=trim(clpath), status='new',form='unformatted')
           close (jpfunf2)

           fname(2:2)='3'
           fname(3:3)=fogtype
           clpath=trim(coutdir)//trim(fname)
           open (jpfunf3, file=trim(clpath), status='new',form='unformatted')
           close (jpfunf3)
        end if

     end if
  end if

! output radiative fluxes and heating rates
  fname='profr .out'
  fname(6:6)=fogtype
  open (jpfunprofr, file=fname,status='unknown')

end subroutine openm

!
!-------------------------------------------------------------
!

subroutine openc (fogtype)
! input/output files of chemical species


! Author:
! ------
  !    RvG?


! Modifications :
! -------------
  !

! == End of header =============================================================
  USE config, ONLY : &
       binout, &
       coutdir, &
       rst

  USE file_unit, ONLY : &
       jpfunprofc, jpfunmass, &
       jpfunsg1, jpfunion, jpfunsl1, jpfunsr1, jpfungr, jpfungs, &
       jpfunjra
  character (len=10) :: fname
  character (len=1)  :: fogtype
  character (len=150) :: clpath  ! complete path to file

! chemical concentrations for initialization
!      this is now done in SR mk_interface (see utils.f90)

! vertical profiles of chemical species
  fname='profc .out'
  fname(6:6)=fogtype
  open (jpfunprofc,file=fname,status='unknown')

! all plotfiles for chemical species
! Create output files by opening them once, only if no restart
  if (.not.rst) then
     ! plout* output
     if (binout) then
        fname='sg1 .out'
        fname(4:4)=fogtype
        clpath=trim(coutdir)//trim(fname)
        open (jpfunsg1, file=trim(clpath), status='new',form='unformatted')
        !close (jpfunsg1) !do not close, used in kpp without opening
        fname='ion .out'
        fname(4:4)=fogtype
        clpath=trim(coutdir)//trim(fname)
        open (jpfunion, file=trim(clpath), status='new',form='unformatted')
        close (jpfunion)
        fname='sl1 .out'
        fname(4:4)=fogtype
        clpath=trim(coutdir)//trim(fname)
        open (jpfunsl1, file=trim(clpath), status='new',form='unformatted')
        close (jpfunsl1)
        fname='sr1 .out'
        fname(4:4)=fogtype
        clpath=trim(coutdir)//trim(fname)
        open (jpfunsr1, file=trim(clpath), status='new',form='unformatted')
        !close (jpfunsr1) !do not close, used in kpp without opening
        fname='gr .out'
        fname(3:3)=fogtype
        clpath=trim(coutdir)//trim(fname)
        open (jpfungr, file=trim(clpath), status='new',form='unformatted')
        close (jpfungr)
        fname='gs .out'
        fname(3:3)=fogtype
        clpath=trim(coutdir)//trim(fname)
        open (jpfungs, file=trim(clpath), status='new',form='unformatted')
        close (jpfungs)
        fname='jra .out'
        fname(4:4)=fogtype
        clpath=trim(coutdir)//trim(fname)
        open (jpfunjra, file=trim(clpath), status='new',form='unformatted')
        close (jpfunjra)
     end if
  end if

  clpath=trim(coutdir)//'mass.out'
  open (jpfunmass,file=trim(clpath),status='unknown')
  write (jpfunmass,101)
  write (jpfunmass,102)
  write (jpfunmass,103)
  close (jpfunmass)
 101  format ('output of molecule burden/deposit/source; unit is', &
     & ' [mol/m2]')
 102  format ('to get balance: divide last output by first; to get', &
     & ' emitted salt mass (in [g/m2])')
 103  format ('multiply xnass with 68.108 (=23 g(Na)/mol(Na) / 0.3377', &
     & ' g(Na)/g(seasalt))')

end subroutine openc

!
!-------------------------------------------------------------
!

subroutine initm (fogtype,rst) !change also SR surf0 !_aerosol_nosub


! Author:
! ------
  !    Andreas Bott, RvG


! Modifications :
! -------------
  ! 12-May-2021  Josue Bock  Initialisation of buoy, allows smoother model spinup (see SR atk1:
  !                            filtering 80% old + 20% new values, but old has to be initialised though)
  !                          Introduce nwProfOpt for different subsidence profiles.
  !                            1- original BTZ96 paper. 2- current parameterisation

! == End of header =============================================================


  USE config, ONLY : &
! Imported Routines:
       abortM,&
! Imported Scalar Variables with intent (in):
       nyear, nmonth, nday, nhour, &
       zalat=>alat, alon, & ! mind that alat is already used in cb16, import alat from config as zalat
       jpPartDistSet, iaertyp, &
       rp0, zinv, dtinv, xm1w, xm1i, rhMaxBL, rhMaxFT, &
       ug, vg, nuvProfOpt, wmin, wmax, nwProfOpt, &
       isurf, binout, coutdir, &
       chamber, &
       lpBuys13_0D, lpJoyce14bc


  USE constants, ONLY : &
! Imported Parameters:
       cp,           &
       g,           &
       pi, rad,          &
       r0,            &      ! Specific gas constant of dry air, in J/(kg.K)
       r1,            &      ! Specific gas constant of water vapour, in J/(kg.K)
       rhow                  ! Water density [kg/m**3]

  USE data_surface, ONLY : &
       ebs, &                    ! saturation moisture potential [cm^3 cm^-3]
       tw, &                       ! water surface temperature
       ustern, z0,&                ! frictional velocity, roughness length
       gclu, gclt                  ! coefficients for momentum, and temperature and humidity

  USE file_unit, ONLY : &
       !jpfunae,            & ! initial aerosol distribution from previous run (if used)
       jpfunerr, jpfunout, & ! standard error/output files
       jpfunpb,            & ! ploutm files: pb*
       jpfunfi, jpfunpi      ! initial output for plotting

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nb, &
       nka, &
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none
! initial profiles of meteorological variables

  character (len=1), intent(in) :: fogtype
  logical, intent(in) :: rst
! External function:
  real (kind=dp), external :: p21              ! saturation water vapour pressure [Pa]

  integer, parameter :: jpdaypermonth(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  real (kind=dp), parameter :: gamma = 0.0098_dp ! = g / cp = dry adiabatic lapse rate (K/m)
  real (kind=dp), parameter :: xmol2 = 18._dp

  character (len=10) :: fname   ! file name
  character (len=110) :: clpath ! path to file
  integer :: jm ! running indexes
  integer :: idayjul, itotyear, istort, immort ! julian day, total day per year, local hr, local min
  integer :: ia, jt, k, k0, ka
  real (kind=dp) :: cc, ctq, cu, dd, deltat
  real (kind=dp) :: poben, punten, rdec, rk, tkorr, x0, xm21s, xnue, zgamma
  real (kind=dp) :: vbt, zp, zpdl, zpdz0

! Local arrays
  real(kind=dp) :: wn(4,3),wr(4,3),ws(4,3),sr(nka,nkt)
!  dimension aer(nf,nka)
!  dimension fnorm(n)

! Common blocks:
  common /cb18/ alat,declin                ! for the SZA calculation
  real(kind=dp) :: alat,declin

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real(kind=dp) :: detw, deta, eta, etw

  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy

  common /cb44/ a0m,b0m(nka)
  real(kind=dp) :: a0m,b0m

  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w
  common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb), &
                ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
  real (kind=dp) :: zb, dzb, dzbw, tb, eb, ak, d, &
       ajb, ajq, ajl, ajt, ajd, ajs, ds1, ds2, ajm, reif, tau, trdep
  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb51/ dlgew,dlgenw,dlne
  real (kind=dp) :: dlgew, dlgenw, dlne

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb53a/ thet(n),theti(n)
  real(kind=dp) :: thet, theti
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a
  common /cb63/ fcs(nka),xmol3(nka)
  real(kind=dp) :: fcs, xmol3
  common /ff_0/ ff_0(nka)
  real(kind=dp) :: ff_0
  common /kinv_i/ kinv
  integer :: kinv

! Statement functions:
  real (kind=dp) :: dfdlogr, dfdlogr2, rr
  real (kind=dp) :: dfdlogrch, dfdlogrch2 ! chamber version

! aerosol distribution; f=dfdlogr*dlogr=dfdlogr*dlgenw/3
  dfdlogr(rr,ka)=wn(ka,1)*exp(-ws(ka,1)*log10(rr/wr(ka,1))**2)+ &
                 wn(ka,2)*exp(-ws(ka,2)*log10(rr/wr(ka,2))**2)+ &
                 wn(ka,3)*exp(-ws(ka,3)*log10(rr/wr(ka,3))**2)
  dfdlogr2(rr,ka)=wn(ka,1)*exp(-ws(ka,1)*log10(rr/wr(ka,1))**2)+ &
                  wn(ka,2)*exp(-ws(ka,2)*log10(rr/wr(ka,2))**2)
! after Jaenicke/Sander/Kim:
!      dfdlogr(rr,ka)=2.8d2/(0.1106*sqrt(2*pi))*exp(-log10(rr/8.8d-2) &
!     & **2/(2*0.1106**2))+ &
!     &               6.6d-1/(0.1906*sqrt(2*pi))*exp(-log10(rr/1.7d0) &
!     & **2/(2*0.1906**2))
! see below: call adjust_f

  dfdlogrch(rr)=wn(2,3)/(ws(2,3)*sqrt(2*pi))*exp(-log10(rr/wr(2,3))**2/(2*ws(2,3)**2))+ &
                wn(2,2)/(ws(2,2)*sqrt(2*pi))*exp(-log10(rr/wr(2,2))**2/(2*ws(2,2)**2))
  dfdlogrch2 = 1.0e-100_dp

! == End of declarations =======================================================

! constants for aerosol distributions
! 3 modes j=1,2,3
! 4 aerosol types i=iaertyp: 1=urban; 2=rural; 3=ocean; 4=background
  select case (jpPartDistSet)
  case (0)
! Distribution after jaenicke (1988)
     ! ((wn(i,j),i=1,4),j=1,3) &
     wn = reshape( &
          (/1.6169e+05_dp, 1.1791e+04_dp,  80.76_dp, 79.788_dp, &
                 664.9_dp,     105.29_dp, 126.52_dp, 94.138_dp, &
            4.3091e+04_dp, 2.9846e+03_dp, 3.0827_dp, 0.0596_dp/), &
            shape=(/4,3/) )
     ! ((wr(i,j),i=1,4),j=1,3) &
     wr = reshape( &
          (/6.51e-03_dp, 7.39e-03_dp, 3.9e-03_dp, 3.6e-03_dp, &
            7.14e-03_dp,    .0269_dp,    .133_dp,    .127_dp, &
               .0248_dp,    .0419_dp,     .29_dp,    .259_dp/), &
            shape=(/4,3/) )
     ! ((ws(i,j),i=1,4),j=1,3) &
     ws = reshape( &
          (/8.3299_dp, 9.8765_dp, 1.1583_dp, 1.2019_dp, &
            1.1273_dp, 1.6116_dp, 11.338_dp, 7.8114_dp, &
            4.4026_dp, 7.0665_dp, 3.1885_dp, 2.7682_dp/), &
            shape=(/4,3/) )

  case (1)
! constants for aerosol distributions after jaenicke (1988)
! except constants for maritime aerosol distribution
! after Hoppel et al. 1990 JGR 95, pp. 3659-3686
     wn = reshape( &
          (/1.6169e+05_dp, 1.1791e+04_dp, 159.576_dp, 79.788_dp, &
                 664.9_dp,     105.29_dp, 427.438_dp, 94.138_dp, &
            4.3091e+04_dp, 2.9846e+03_dp,   5.322_dp, 0.0596_dp/), &
            shape=(/4,3/) )
     wr = reshape( &
          (/6.51e-03_dp, 7.39e-03_dp, 0.027_dp, 3.6e-03_dp, &
            7.14e-03_dp,    .0269_dp,  .105_dp,    .127_dp, &
               .0248_dp,    .0419_dp,   .12_dp,    .259_dp/), &
            shape=(/4,3/) )
     ws = reshape( &
          (/8.3299_dp, 9.8765_dp,    8._dp, 1.2019_dp, &
            1.1273_dp, 1.6116_dp, 39.86_dp, 7.8114_dp, &
            4.4026_dp, 7.0665_dp, 2.469_dp, 2.7682_dp/), &
            shape=(/4,3/) )

  case (2)
! maritime size distribution after hoppel et al. 1994, jgr. 14,443
     wr(3,1) = 0.02_dp
     wr(3,2) = 0.05_dp
     wr(3,3) = 0.15_dp
     wn(3,1) = 110._dp
     wn(3,2) = 72._dp
     wn(3,3) = 7._dp
     ws(3,1) = 0.14_dp
     ws(3,2) = 0.16_dp
     ws(3,3) = 0.18_dp
     x0 = sqrt(2._dp * pi)
     do k=1,3
        wn(3,k) = wn(3,k)/(x0*ws(3,k))
        ws(3,k) = 1._dp / (2._dp * ws(3,k)**2)
     enddo

     if (iaertyp .ne. 3) then
        write(jpfunerr,*) "Inconsistency in namelist settings:"
        write(jpfunerr,*) "  jpPartDistSet=2 can only be used with maritime aerosol (iaertyp=3)"
        call abortM ("Error in initm")
     end if

  case (3)
! polar size distribution after Jaenicke, 1988, #164
     wr(3,1) = 6.89e-2_dp
     wr(3,2) = 3.75e-1_dp
     wr(3,3) = 4.29_dp
     wn(3,1) = 21.7_dp
     wn(3,2) = 0.186_dp
     wn(3,3) = 3.04e-4_dp
     ws(3,1) = 0.245_dp
     ws(3,2) = 0.300_dp
     ws(3,3) = 0.291_dp
     x0 = sqrt(2._dp * pi)
     do k=1,3
        wn(3,k) = wn(3,k)/(x0*ws(3,k))
        ws(3,k) = 1._dp / (2._dp * ws(3,k)**2)
     enddo

     if (iaertyp .ne. 3) then
        write(jpfunerr,*) "Inconsistency in namelist settings:"
        write(jpfunerr,*) "  jpPartDistSet=3 can only be used with maritime aerosol (iaertyp=3)"
        call abortM ("Error in initm")
     end if

  case (4)
! special chamber case
     wn(1:2,1:3) = reshape( &
          (/   0.0_dp, 0.0_dp, 0.0_dp,        &
          -0.17e+2_dp, 0.0_dp, 0.53e+02_dp/), &
          shape=(/2,3/) )
     wr(1:2,1:3) = reshape( &
          (/ 0.0_dp, 0.0_dp, 0.0_dp,     &
             1.4_dp, 0.0_dp, 0.357_dp/), &
          shape=(/2,3/) )
     ws(1:2,1:3) = reshape( &
          (/ 0.0_dp, 0.0_dp, 0.0_dp,        &
          -0.125_dp, 0.0_dp, 0.126_dp/), &
          shape=(/2,3/) )

  case default
     write(jpfunerr,*) "Wrong choice for jpPartDistSet"
     call abortM ("Error in initm")
  end select

! initialisation of solar time: nyear, nmonth, nday, nhour are read from namelist
! -----------------------------
  if (.not.rst) then
! day of year = julian day
     idayjul = 0
     do jm=1,nmonth-1
        idayjul = idayjul + jpdaypermonth(jm)
     end do
     ! leap year?
     if (mod(nyear,4).eq.0 .and. nmonth.ge.3) then
        idayjul = idayjul + 1
        itotyear = 366
     else
        itotyear = 365
     end if
     idayjul = idayjul + nday

! starting time of calculations
     lday = 0
     lmin = 0
     lst = nhour
     if (chamber) lst = 12

! equation of time:
     ! gamma is the day angle, also called fractional year, in radians
     zgamma  = 2._dp * pi * real(idayjul - 1,dp) / real(itotyear,dp)
     ! equation of time (in hour)
     deltat = 24._dp / (2._dp * pi) * (0.0000075_dp &
          + 0.001868_dp * cos(zgamma) - 0.032077_dp * sin(zgamma) &
          - 0.014615_dp * cos(2._dp * zgamma) - 0.040849_dp * sin(2._dp * zgamma))
! time correction
     tkorr  = (4._dp * alon / 60._dp) + deltat
     ! Local time, for output only
     istort = lst + floor(tkorr)
     immort = nint((tkorr + real(lst - istort,dp)) * 60._dp)

! declination of the sun
     rdec   = 0.006918_dp - 0.399912_dp * cos(zgamma) + 0.070257_dp * sin(zgamma) &
          - 0.006758_dp * cos(2._dp * zgamma) + 0.000907_dp * sin(2._dp * zgamma)
     ! in cb16, declination is in degree, convert here
     declin = rdec / rad
     if (chamber) declin = 18._dp

! Fill alat in cb16 with alat=zalat from config
     alat = zalat

! output
     write (jpfunout,6000)
     write (jpfunout,6300)
     write (jpfunout,6301) nday, nmonth, nyear
     write (jpfunout,6302) int(idayjul)
     write (jpfunout,6303) tkorr
     write (jpfunout,6304) istort,immort
     write (jpfunout,6305) declin
     write (jpfunout,6000)

6000 format ('---------------------------------------------------------------&
             &-----------------')
6300 format ('Initialization of local time:')
6302 format ('Julian day: ', i3)
6301 format ('Date: ',i2.2,'.',i2.2,'.',i4.4)
6303 format ('Time correction: ',f5.2,' hours')
6304 format ('Start at local time (incl. time correction): ',i3.2,':',i2.2)
6305 format ('Declination of the sun: ',f6.2,' deg')

  end if

! initialisation of temperature and humidity profile
!---------------------------------------------------

  if (.not.rst) then
! initial inversion height zinv, find the corresponding layer
     kinv = 1
     do k=2,nf
        if (zinv.gt.eta(k) .and. zinv.le.eta(k+1)) then
           kinv = k
           exit
        end if
     enddo
     if (kinv.eq.1) then
        write (jpfunerr,*) 'Error in SR initm: zinv = ',zinv
        write (jpfunerr,*) 'It must be set lower than the uppermost prognostic layer'
        write (jpfunerr,*) 'whose height is: ',eta(nf)
        call abortM ('Error in SR initm')
     end if

! initial temperature profile
! dry adiabatic profile below inversion height
     t(1)=tw
      do k=2,kinv
         t(k)=t(1)-gamma*eta(k)
      enddo
! stable profile above inversion height
      x0=t(kinv)+dtinv
      do k=kinv+1,n
         t(k)=x0-0.006_dp*(eta(k)-eta(kinv))
      enddo

! large scale hydrostatic pressure
      ! jjb: in this formulation, the ground pressure (rp0) is taken as the bottom boundary of
      !      layer #1, whose thickness (detw(1)) is NOT zero in this context (see SR grid:
      !      detw(1)=detamin). This ensures the continuity for the calculation of diffusional
      !      processes (see for instance difm: xc coefficients). This leads to a minor shift in
      !      the pressure vertical profile.
      poben = rp0
      cc = g / (2._dp * r0)
      do k=1,n
         punten = poben
         dd     = detw(k) * cc / t(k)
         poben  = punten * (1._dp - dd) / (1._dp + dd)
         p(k)   = 0.5_dp * (poben + punten)
      enddo

! initial profiles of humidity and wind field
! xm1: = specific humidity in kg/kg
! xm2: = liquid water content in kg/m**3
      do k=1,n
         talt(k) = t(k)
         thet(k) = (p(1)/p(k))**0.286_dp
         theti(k) = 1._dp / thet(k)
         theta(k) = t(k)*thet(k)
! r0/r1=287.05/461.51=0.62198; 1-r0/r1=0.37802;
         xm21s = 0.62198_dp * p21(t(k)) / (p(k) - 0.37802_dp * p21(t(k)))
         if (k <= kinv) then
            xm1(k) = min(xm1w, rhMaxBL*xm21s)
         else
            xm1(k) = min(xm1i, rhMaxFT*xm21s)
         end if
         xm1a(k) = xm1(k)
         xm2(k)  = 0._dp
         feu(k) = xm1(k)*p(k)/((0.62198_dp + 0.37802_dp * xm1(k))*p21(t(k)))
! r1/r0-1=0.61
         rho(k) = p(k)/(r0*(t(k)*(1._dp + 0.61_dp * xm1(k))))
         thetl(k) = theta(k)*(1._dp + 0.61_dp * xm1(k))
         dfddt(k) = 0._dp
         buoy(k) = -1e-4_dp

! profile of geostrophic wind
         select case (nuvProfOpt)
         case (0)
            u(k) = ug
            v(k) = vg
         case (3) ! Bott 2020 (p.160): decreasing linearly to zero at the surface
            if (k <= kinv) then
               u(k) = ug/zinv * eta(k)
               v(k) = vg/zinv * eta(k)
            else
               u(k) = ug
               v(k) = vg
            end if
         case default
            call abortM ('Wrong option for nuvProfOpt, default = 0, other possible option = 3')
         end select

! profile of large scale subsidence
         select case (nwProfOpt)
         ! Hyperbolic form: Bott et al 1996
         case (1)
            w(k) = 0.5_dp * wmax * (tanh((eta(k)-500._dp) / 250._dp) + 1._dp)

         case (2)
            w(k) = eta(k)/1000._dp * 0.5_dp * (wmin+wmax)
         ! Bott 2020: wmax above kinv, decreasing linearly to zero (wmin) at the surface
         case (3)
            if (k <= kinv) then
               w(k) = (wmax-wmin)/zinv * eta(k) + wmin
            else
               w(k) = wmax
            end if
         case default
            call abortM ('Wrong option for nwProfOpt, choose 1, 2 or 3')
         end select

         if (k <= kinv) then
            tke(k) = 0.05_dp
         else
            tke(k) = 1.e-5_dp
         end if
      enddo
      if (nuvProfOpt == 0) then
         u(1) = 0._dp
         v(1) = 0._dp
         u(2) = 0.25_dp * ug
         v(2) = 0.25_dp * vg
         u(3) = 0.75_dp * ug
         v(3) = 0.75_dp * vg
      end if
      do k=n,1,-1
         w(k) = w(k) - w(1)
      enddo

! initial cloud layers
     lcl=1
     lct=1

  end if

! no restart only for the following block
  if (.not.rst) then

     ff(:,:,:) = 0._dp
     fsum(:)   = 0._dp
     do k=1,n
        nar(k) = iaertyp

        ! Define a scaling factor to tune the original distributions
        ! General case
        if (.not.lpJoyce14bc) then
           x0 = 1._dp
           if (iaertyp.lt.3.and.k.gt.nf) x0 = 0.2_dp
           !if (iaertyp.lt.3.and.k.gt.kinv) x0 = 0.2_dp
        ! special case for Joyce et al 2014 study
        else if (lpJoyce14bc) then
           x0 = 1.e-4_dp
        end if

        do ia=1,nka
           if (.not.chamber) then
              ff(1,ia,k)=dfdlogr(rn(ia),nar(k))*dlgenw/3.*x0
              ! special case: save the distribution to replenish (save only once, k==1)
              if (lpJoyce14bc.and.k==1) ff_0(ia)=dfdlogr(rn(ia),nar(k))*dlgenw/3.
              if (k.gt.kinv) ff(1,ia,k)=dfdlogr2(rn(ia),nar(k))*dlgenw/3.*x0
           else
              ff(1,ia,k)=dfdlogrch(rn(ia))*dlgenw/3.*x0
              if (k.gt.kinv) ff(1,ia,k)=dfdlogrch2*dlgenw/3.*x0
           end if
           fsum(k) = fsum(k) + ff(1,ia,k)
        enddo
!        write (199,*)"k,fsum",k,fsum(k)
!        fnorm(k)=fsum(k)
      enddo

! read initial aerosol distribution from previous run
!#      fname='ae .out'
!#      fname(3:3)=fogtype
!#      open (jpfunae, file=fname,status='old',form='unformatted')
!#      read (jpfunae) aer
!#      close (jpfunae)
!#      do k=1,nf
!#         fsum(k)=0.
!#         do ia=1,nka
!#            ff(1,ia,k)=aer(k,ia)
!#            fsum(k)=fsum(k)+ff(1,ia,k)
!#         enddo
!#      enddo
! normalization of aerosol size distribution
!      do k=2,nf
!         x0=fnorm(k)/fsum(k)
!         fsum(k)=0.
!         do ia=1,nka
!            ff(1,ia,k)=ff(1,ia,k)*x0
!            fsum(k)=fsum(k)+ff(1,ia,k)
!         enddo
!      enddo

   end if ! no restart

! parameters a0m, b0m of koehler curve of subroutine subkon:
! sr(ia,jt)=exp(a0m/(r(ia,jt)*t)-b0m(ia)*en(ia)/ew(jt)))
! a0m see p. 141 pruppacher and klett a0m=2 sigma/(r1*t*rhow*a)
! 152200= 2 sigma*10**6 with sigma is surface tension = 76.1*10**-3
! see pruppacher and klett p. 104
  if (.not.rst) then
     a0m = 152200._dp / (r1*rhow)
  !else restart case: a0m is read in startm
  end if

! b0m=fcs*xnue*xmol2/xmol3: fcs(ia) fraction of soluble species
! xnue number of ions; xmol2 (xmol3) mol masses of water (aerosol)
  k0=nar(2) !=iaertyp
  do ia=1,nka
! aerosol types: 1=urban 2=rural 3=ocean 4=tropospheric
     select case (k0)

! Aerosol type 1 = urban
     case (1)
! NH4NO3 mole mass 80; (NH4)2SO4 mole mass 132
! soluble part of urban aerosol: 2 mole NH4NO3 and 1 mole (NH4)2SO4

      ! general case
        if (.not.lpJoyce14bc) then
           if (rn(ia).le.1._dp) then
              fcs(ia) = 0.4_dp - rn(ia) * (0.4_dp - 0.1_dp)
           else
              fcs(ia) = 0.1_dp
           end if
           xnue = (3._dp + 2._dp * 2._dp) / 3._dp
           xmol3(ia) = (132._dp + 80._dp * 2._dp) / 3._dp

      ! special case for Joyce et al 2014 study
      ! making urban aerosol H2SO4 with remaining mass DOM and gas uptake
        else if (lpJoyce14bc) then
           if (rn(ia).le.1._dp) then
              fcs(ia) = 0.9_dp - rn(ia) * (0.9_dp - 0.5_dp)
           else
              fcs(ia) = 0.1_dp
           end if
         ! "average" number of dissoc. ions, H+,SO4=,DOM, Cl- (4)
           xnue = (3._dp ) / 4._dp
         ! "average" MW of ion = MW components/#
         ! sulfacid(98.08)+octene(122.21)+Cl- / 4 diss. com
           xmol3(ia) = (98.08_dp + 122.21_dp + 35.45_dp) / 4._dp
        end if


! Aerosol type 2 = rural
     case (2)
! soluble part of rural aerosol: pure (NH4)2SO4
        if (rn(ia).le.1._dp) then
           !fcs(ia) = 0.5_dp - rn(ia) * (0.5_dp - 0.1_dp)
           fcs(ia) = 0.9_dp - rn(ia) * (0.9_dp - 0.5_dp)
        else
           !fcs(ia) = 0.1_dp
           fcs(ia) = 0.5_dp
        end if
        xnue = 3._dp
        xmol3(ia) = 132._dp

! Aerosol type 3 = maritime
     case (3)
!c soluble part of ocean aerosol: small pure (NH4)2SO4; large pure NaCl
! soluble part of ocean aerosol: pure (NH4)2SO4;
        fcs(ia) = 1._dp
!        xnue      = 3._dp
!        xmol3(ia) = 132._dp
! 32% (NH4)2SO4, 64% NH4HSO4, 4% NH4NO3
        xnue      = 0.32_dp*3. + 0.64_dp*2. + 0.04_dp*2.
        xmol3(ia) = 0.32_dp*132. + 0.64_dp*115. + 0.04_dp*80.
! large are NaCl
        if (rn(ia).ge.0.5_dp) then
           xnue      = 2._dp     !no change in microphysics due to halogen chemistry
           xmol3(ia) = 58.4_dp
        end if

        if (lpBuys13_0D) then
           fcs(ia)   = 0._dp ! jjb: there was a mistake in Buys settings, fcs was commented out, thus not set
                             !      use 0 instead, even if this is wrong, just to reproduce the same settings
           xnue      = 2._dp
           xmol3(ia) = 58.4_dp
        end if

     end select

     if (.not.rst) then
        b0m(ia)=fcs(ia)*xnue*xmol2/xmol3(ia)
        !else restart case: b0m is read in startm
     end if
  end do

  if (.not.rst) then

! in radiation code: background aerosol = rural aerosol
     do k=2,n
        if (nar(k).eq.4) nar(k)=2
     enddo


!      call adjust_f !only if Kim aerosol is used!!


! initial calculation of Clarke-functions and frictional velocity ustern
     vbt   = sqrt(u(2)*u(2)+v(2)*v(2))
     zp    = deta(1) + z0
     zpdz0 = log(zp/z0)
     zpdl  = g * (theta(2) - t(1)) * zp / (theta(2) * vbt)
     call claf (zpdl,zpdz0,cu,ctq)
     ustern = max(0.01_dp, vbt/cu)
     gclu   = cu
     gclt   = ctq

     ajs   = 0._dp
     ds1   = 0._dp
     ds2   = 0._dp
     trdep = 0._dp
     tau   = 0._dp
     reif  = 0._dp

! temperature and volumetric moisture within soil
     if(isurf == 1) then
        x0 = 0.5_dp * ebs
!        if (iaertyp.eq.1) x0=x0*.9
        do k=1,nb
           tb(k) = 285.0_dp
           eb(k) = x0
           if (zb(k).lt.0.1_dp) tb(k)=(t(1)*(0.1_dp - zb(k)) + 285._dp*zb(k)) / 0.1_dp
        enddo
        if (binout) then
           write (jpfunpb) zb
           close (jpfunpb)
        end if
     end if
  end if ! no restart only

! initial output for plotting: both restart and no-restart cases
  fname='pi .out'
  fname(3:3)=fogtype
  clpath = trim(coutdir)//fname
  open (jpfunpi, file=trim(clpath),status='unknown')
  write (jpfunpi,6001) (eta(k),etw(k),rho(k),p(k),w(k),k=1,n)
 6001 format (5e16.8)
  close (jpfunpi)

  do jt=1,nkt
     do ia=1,nka
        rk=rw(jt,ia)
        sr(ia,jt)=max(0.1_dp, exp(a0m/(rk*t(2)) - b0m(ia)*en(ia)/ew(jt)))
     enddo
  enddo
  fname='fi .out'
  fname(3:3)=fogtype
  clpath = trim(coutdir)//fname
  open (jpfunfi, file=clpath,status='unknown')
  write (jpfunfi,6010) rn,en,rq,e,sr
  close (jpfunfi)
 6010 format (5e16.8)

end subroutine initm

!
!-------------------------------------------------------------
!

subroutine grid
!
! Description:
! -----------
  ! numerical grid for the atmosphere, soil, aerosols and water droplets


! Author:
! ------
  !    Andreas Bott
  !    Roland von Glasow (aerosol + cloud chemistry part)


! Modifications :
! -------------
  ! Josue Bock : minor bugfix in dlgenw and dlgew calculations: divide by nka and nkt
  !              instead of nka-1 nkt-1, respectively (to be consistent with enw(1) and
  !              ew(1) expressions: lower bound of first class, higher bound of last class)

! == End of header =============================================================


  USE config, ONLY : &
! External subroutine
       abortM, &
       chamber, &
       isurf, &
       detamin, etaw1,       & ! atmospheric grid: constant height, and top of grid [m]
       rnw0, rnw1, rw0, rw1    ! microphysics grid: min/max dry aerosol, min/max particle radius [um]

  USE constants, ONLY : &
! Imported Parameters:
       pi, &
       rho3, &               ! Aerosol density [kg/m**3]
       rhow                  ! Water density [kg/m**3]

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunerr, jpfunout

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nb, &
       nka, &
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local parameters:
  real (kind=dp), parameter :: etaw1_max = 2500._dp ! max authorised height for prognostic grid
  real (kind=dp), parameter :: dzbw0 = 0.001_dp, zbw1 = 1._dp ! minimum-maximum values of soil grid
! Local scalars:
  integer :: ia, ij, iij, jt, j, k, nf1, ndistrp
  real (kind=dp) :: ax, dlgzbw
  real (kind=dp) :: enwmin, enwmax
  real (kind=dp) :: ewmin, ewmax
  real (kind=dp) :: rqmin, rqmax, rfact
  real (kind=dp) :: x0, x1, x2, x3
  real (kind=dp) :: xfac
  real (kind=dp) :: zbw0, zbw, zradthres
! Local arrays:
  real (kind=dp) :: rq1D(nkt*nka),rqcp(nkt,nka)

! Common blocks:
  common /blck06/ kw(nka),ka
  integer :: kw, ka
  common /cb08/ re1(nkt), re2(nkt), re3(nkt)
  real (kind=dp) :: re1, re2, re3
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw
  common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb), &
       &              ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
  real (kind=dp) :: zb, dzb, dzbw, tb, eb, ak, d, &
       ajb, ajq, ajl, ajt, ajd, ajs, ds1, ds2, ajm, reif, tau, trdep
  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
       &              e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq
  common /cb51/ dlgew,dlgenw,dlne
  real (kind=dp) :: dlgew, dlgenw, dlne
  common /oneDsj/ rpw(nka), part1D(nka-1,nf)
  real (kind=dp) :: rpw, part1D

! == End of declarations =======================================================

! atmospheric vertical grid
! -------------------------

! equidistant grid between earth's surface and eta(nf)
  etw(1) = 0._dp
  do k=2,nf
     etw(k) = real(k-1,dp) * detamin
  end do

! logarithmically equidistant grid above eta(nf)
! maximum prognostic height: etaw1
  nf1 = nf +1

! check input
  ! Case 1: detamin too high as compared to etaw1, or too many layers requested between nf and n
  if (etaw1 .lt. etw(nf) + (n-nf)*detamin) then
     write (jpfunerr,*) 'Error: impossible to get n-nf layers with increasing height'
     write (jpfunerr,*) '  decrease detamin, or increase etaw1, or change the layer numbers'
     call abortM ('Error in SR grid')
  end if
  ! Case 2: etaw1 set too high
  if (etaw1 .gt. etaw1_max) then
     write (jpfunout,*) 'Warning: uppermost prognostic level input too high'
     write (jpfunout,*) 'The value has been reset: etw1 = 2500 m'
     etaw1 = etaw1_max
  end if

  j  = 0
  x0 = detamin
  x1 = etaw1
  do while (x1 .gt. etaw1-etw(nf))
     x0 = x0 + detamin
     j = j + 1
     x3 = detamin / x0 + 1._dp
     etw(nf1) = x0
     do k=nf+2,n
        etw(k) = etw(k-1) * x3
     enddo
     x1 = etw(n) - etw(nf1)
     if (j.gt.10000) call abortM ('Error in atmospheric grid')
  end do

! fill the etw grid
  x0 = nf * detamin - etw(nf1)
  do k=nf1,n
     etw(k) = etw(k) + x0
  enddo
! deduce detw, eta and deta from etw
  detw(1) = detamin ! jjb: this is necessary for the boundary continuity in the diffusion scheme (and other?), not a mistake
  eta(1)  = 0._dp
  do k=2,n
     detw(k)   = etw(k) - etw(k-1)
     eta(k)    = 0.5_dp * (etw(k) + etw(k-1))
     deta(k-1) = eta(k) - eta(k-1)
  enddo
  deta(n) = (1._dp + x3) * 0.5_dp * etw(n) - eta(n)

! grid within the soil
! --------------------
  if (isurf == 1) then
     zbw0 = 0._dp
     x2 = 0._dp
     do while (x2.lt.dzbw0)
        zbw0   = zbw0 + 0.0001_dp
        dlgzbw = log10(zbw1/zbw0)/real(nb,dp)
        x3     = 10._dp**dlgzbw
        zbw    = zbw0 * x3
        x2     = zbw-zbw0
     end do

     zb(1)   = zbw
     dzbw(1) = zbw-zbw0
     do k=2,nb
        zbw0     = zbw
        zbw      = zbw0 * x3
        zb(k)    = 0.5_dp * (zbw + zbw0)
        dzbw(k)  = zbw - zbw0
        dzb(k-1) = zb(k) - zb(k-1)
     enddo
     dzb(nb) = (1._dp + x3) * 0.5 * zbw - zb(nb)
     ! substract zb(1) to the whole grid so that the first layer has 0-depth
     x0 = zb(1)
     do k=1,nb
        zb(k)=zb(k)-x0
     enddo
  end if

! aerosol grid
! ------------
  x0=1._dp/3._dp
  x1=4._dp*x0*pi*rhow
  x2=4._dp*x0*pi*rho3

! aerosol mass en(i) in mg: prescribed lowest and largest aerosol class
  enwmin = x2 * rnw0**3 * 1.e-12_dp
  enwmax = x2 * rnw1**3 * 1.e-12_dp
! logarithmic equidistant mass grid en(i), enw(i)
  dlgenw = log10(enwmax/enwmin) / real(nka,dp)
  x3 = 10._dp**dlgenw

  enw(1) = enwmin * x3
  en(1)  = 0.5_dp * (enw(1) + enwmin)
  rn(1)  = (en(1) / x2)**x0 * 1.e4_dp
  do ia=2,nka
     enw(ia) = enw(ia-1) * x3
     en(ia)  = 0.5_dp * (enw(ia) + enw(ia-1))
     rn(ia)  = (en(ia)/x2)**x0 * 1.e4_dp
  enddo

! water mass ew in mg  (lowest and largest class)
  ewmin = x1 * rw0**3 * 1.e-12_dp
  ewmax = x1 * rw1**3 * 1.e-12_dp
  dlgew = log10(ewmax/ewmin) / real(nkt,dp) ! jjb nkt, not nkt-1
  ax = 10._dp**dlgew
! ax: growth factor for consecutive masses
! ln(10)=2.3025851  ln(x)=ln(10)*log(x)
  dlne = log(10._dp) * dlgew

  ew(1)  = ewmin * ax
  e(1)   = 0.5_dp * (ew(1)+ewmin)
  dew(1) = ew(1) - ewmin
  do jt=2,nkt
     ew(jt)  = ew(jt-1) * ax
     e(jt)   = 0.5_dp * (ew(jt) + ew(jt-1))
     dew(jt) = ew(jt) - ew(jt-1)
  enddo

  ! equivalent radius of (pure) water droplet (in [m]), for effective radius calculation used in the radiative code
  do jt=1,nkt
     re1(jt) = (e(jt) * 1.e-6_dp / x1)**x0
     re2(jt) = re1(jt)*re1(jt)
     re3(jt) = re1(jt)*re2(jt)
  enddo
  do ia=1,nka
     do jt=1,nkt
        rq(jt,ia) = (e(jt)  * 1.e-6_dp / x1 + (rn(ia)*1.e-6_dp)**3)**x0 * 1.e6_dp
        rw(jt,ia) = (ew(jt) * 1.e-6_dp / x1 + (rn(ia)*1.e-6_dp)**3)**x0 * 1.e6_dp
     enddo
  enddo

!
! Set 1d radius for oneD_dist routine
!

! jjb explanations: setting 1-D radius bounds to project the 2D particle spectrum onto a 1D one for output
!                   is not straightforward.
!                   If the 1D bounds are too small, there will be bounds without particle
!                   leading to ugly output
!                   If 1D bounds are too large, the output will be poorly defined.
!                   Also, by construction, the space between 2D bins is not equal depending on the position
!                   on the 2D spectrum, with much smaller bins for small (a, r) and much coarser
!                   for large (a,r).
!                   Below are various tests

  ndistrp = 7 ! <-- after several tests, satisfying results with this one
  select case (ndistrp)

  ! as many particle bins from rq(:,:) array in each rp(:) array (i.e. 1/70 th percentile spliting)
  case (1)
  ! local copy of rq that will be altered
  rqcp(:,:) = rq(:,:)
  rqmax = maxval(rqcp(:,:))
  ! Sort rq values in increasing order into a 1D array containing all values
  do ij=1,nka*nkt
     rqmin = minval(rqcp(:,:))
     rq1D(ij) = rqmin
     ! remove the min value by setting twice the max
     where(rqcp == rqmin) rqcp = 2._dp * rqmax
  end do
  ! Create rpw by picking one every nka value in rq1D array
  do ij=1,nka+11
     iij = ij + (ij - 1) * (nka-10)
     !iij = ij + (ij - 1) * (nka+10)
     iij = min(iij,nka*nkt)
     rpw(ij) = rq1D(iij)
  end do

  ! linear increase of rp(:)
  case (2)
     rqmin = minval(rq(:,:))
     rqmax = maxval(rq(:,:))
     do ij=1,nka
        rpw(ij) = rqmin + (ij-1) * (rqmax-rqmin)/(nka-1)
     end do

  ! exp increase of rp(:) (thus linear on a log scale)
  case (3)
     rqmin = minval(rq(:,:))
     rqmax = maxval(rq(:,:))
     rfact = 10**(log10(rqmax/rqmin)/(nka-1-10))
     rpw(1) = rqmin
     do ij=2,nka
        rpw(ij)=rpw(ij-1)*rfact
     enddo

  ! same as 1 but with coarser steps for small rq: 1/140th percentile
  case (4)
     ! local copy of rq that will be altered
     rqcp(:,:) = rq(:,:)
     rqmax = maxval(rqcp(:,:))
     ! Sort rq values in increasing order into a 1D array containing all values
     do ij=1,nka*nkt
        rqmin = minval(rqcp(:,:))
        rq1D(ij) = rqmin
        ! remove the min value by setting twice the max
        where(rqcp == rqmin) rqcp = 2._dp * rqmax
     end do
     ! Create rpw by picking one every nka value in rq1D array
     ij=1
     iij=1
     do while (rq1D(iij) <= rq(1,nka) * 4.)
        rpw(ij) = rq1D(iij)
        ij = ij+1
        iij = iij+140
     end do
     print*,'end large steps',ij-1
     iij = iij-70
     do while (iij <= nka*nkt)
        rpw(ij) = rq1D(iij)
        ij = ij+1
        iij = iij+70
     end do
     print*,'end small steps',ij-1
     do while (ij < nka+1)
        rpw(ij) = rpw(ij-1)*1.001
        ij = ij+1
     end do

  ! exp increase of rp(:) (thus linear on a log scale) with 2 distinct slopes
  case (5)
     rqmin = minval(rq(:,:))
     rqmax = maxval(rq(:,:))
     rfact = 10**(log10(rqmax/rqmin)/(nka-1))
     rpw(1) = rqmin
     ij = 1
     ! apply a constant increase factor rfact**1.5 for small radius
     do while (rpw(ij)*rfact*rfact <= rq(1,nka) * 4.)
        ij = ij + 1
        rpw(ij) = rpw(ij-1)*rfact**1.5
     end do
     ! apply a constant increase factor rfact for large radius
     do while (rpw(ij)*rfact <= rq(nkt,nka))
        ij = ij + 1
        rpw(ij) = rpw(ij-1)*rfact
     enddo
     ! fill the remaining indexes with slightly increasing radius, to avoid zeros
     do while (ij<nka)
        ij = ij + 1
        rpw(ij) = rpw(ij-1)*1.001
     end do

  ! another idea: use "diagonal" rw values
  case (6)
     do ia=1,nka
        rpw(ia) = rw(ia,ia)
     end do

  ! another idea (modified): use "diagonal" rq values
  case (62)
     do ia=1,nka
        rpw(ia) = rq(ia,ia)
     end do

  ! similar to 6 but coarser up to rw(1,nka) (1 over 2 bins)
  case (7)
     rpw(1) = rw(1,1)**2 / rw(3,3)
     ij = 2
     ! one every 2 indexes up to rw(1,nka)
     iij = 2*(ij-1)-1
     do while (rw(iij,iij) <= rw(1,nka))
        rpw(ij) = rw(iij,iij)
        ij = ij+1
        iij = 2*(ij-1)-1
     end do
     iij = iij+1
     ! one every index after
     ia = 0
     do while (iij+ia <= nka)
        rpw(ij+ia) = rw(iij+ia, iij+ia)
        ia = ia + 1
     end do
     ! fill the remaining indexes with slightly increasing values
     do while (ij+ia <= nka)
        rpw(ij+ia) = rpw(ij+ia-1)*1.001
        ia = ia + 1
     end do

  ! same as 7 but even coarser: 1 over 3 bins up to rw(1,nka)
  case (72)
     rpw(1) = rw(1,1)**2 / rw(4,4)
     ij = 2
     ! one every 2 indexes up to rw(1,nka)
     iij = 3*(ij-1)-2
     do while (rw(iij,iij) <= rw(1,nka))
        rpw(ij) = rw(iij,iij)
        ij = ij+1
        iij = 3*(ij-1)-2
     end do
     iij = iij+1
     ! one every index after
     ia = 0
     do while (iij+ia <= nka)
        rpw(ij+ia) = rw(iij+ia, iij+ia)
        ia = ia + 1
     end do
     ! fill the remaining indexes with slightly increasing values
     do while (ij+ia <= nka)
        rpw(ij+ia) = rpw(ij+ia-1)*1.001
        ia = ia + 1
     end do


  end select
  ! Initialise part1D
  part1D(:,:) = 0._dp

! aerosol + cloud chemistry: partition into small and large aerosol  : ka
!                       and small and large water content (aer-drop) : kw(nka)
  ka=-1

  ! Radius threshold between small and large aerosol
  if (.not. chamber) then
     zradthres = 0.5_dp
  else
     zradthres = 0.1_dp
  end if

  do ia=1,nka
     if (rn(ia).gt.zradthres .and. ka.lt.0) ka=ia-1
     kw(ia)=-1
  enddo
  if (ka.lt.0) ka=nka
! xfac: dilution as volume ratio: V_dry*x = V_water (r_dry*(x)^(1./3.)=r_water)
  xfac=10._dp ! volume ratio is 1000
  do ia=1,nka
     do jt=1,nkt
        if ((e(jt)*1.e-6_dp/x1)**x0 * 1.e6_dp .gt.xfac*rn(ia).and.kw(ia).lt.0) kw(ia)=jt-1
     enddo
     if (kw(ia).lt.0) kw(ia)=nkt
  enddo
! end aerosol + cloud chemistry

end subroutine grid

!
!-------------------------------------------------------------
!

subroutine startm (fogtype)
! vertical profiles of meteorological data if the program is restarted


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  ! Josue Bock: add some data in rst file to avoid recomputing missing ones here

! == End of header =============================================================


  USE data_surface, ONLY : &
       tw, &                       ! water surface temperature
       ustern, z0,&                ! frictional velocity, roughness length
       gclu, gclt                  ! coefficients for momentum, and temperature and humidity

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunout, jpfunrstm

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nb, &
       nka, &
       nkt, &
       mb

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  character (len=1), intent(in) :: fogtype
  character (len=10) :: fname

! Common blocks:
  common /cb11/ totrad (mb,n)
  real (kind=dp) :: totrad
  common /cb18/ alat,declin                ! for the SZA calculation
  real (kind=dp) :: alat,declin
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct
  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy
  common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
  real (kind=dp) :: gm, gh, sm, sh, xl
  common /cb44/ a0m,b0m(nka)
  real (kind=dp) :: a0m,b0m
  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w
  common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb), &
                ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
  real (kind=dp) :: zb, dzb, dzbw, tb, eb, ak, d, &
       ajb, ajq, ajl, ajt, ajd, ajs, ds1, ds2, ajm, reif, tau, trdep
  common /cb48/ sk,sl,dtrad(n),dtcon(n)
  real (kind=dp) :: sk, sl, dtrad, dtcon
  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb53a/ thet(n),theti(n)
  real(kind=dp) :: thet, theti
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a
  common /cb63/ fcs(nka),xmol3(nka)
  real(kind=dp) :: fcs, xmol3

! == End of declarations =======================================================

  fname='rstm .dat'
  fname(5:5)=fogtype

  open (jpfunrstm,file=fname,status='old',form='unformatted')
  read (jpfunrstm)  &
! double precision arrays
       atkm,atkh,b0m,dfddt,dtrad,eb,ff,fcs,feu,fsum, &
       gh,p,rho,t,talt,tb,thet,theta,theti,thetl,tke,tkep, &
       totrad,u,v,w,xl,xm1,xm1a,xm2,xmol3, &
! double precision single vars
       a0m,alat,declin,ds1,ds2,gclt,gclu,reif,sk,sl,tau, &
       trdep,tw,ustern,z0, &
! integer arrays
       nar, &
! integer single vars
       it,lcl,lct,lday,lmin,lst
  close (jpfunrstm)

  write(jpfunout,*)"restart file for meteo read, filename: ",fname
  write(jpfunout,*)lday,lst,lmin

end subroutine startm

!
!-------------------------------------------------------------
!

subroutine startc (fogtype)
! profiles of chemical data if the program is restarted


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  !

! == End of header =============================================================


  USE config, ONLY : &
       binout, coutdir

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunout, jpfunrstc,&
       jpfunsg1, jpfunsr1

  USE gas_common, ONLY : &
       s1, &
       es1, &
       s3

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j6, &
       nf, &
       n, &
       nka, &
       nkt, &
       nkc, &
       nlev, &
       nrxn, &
       nphrxn

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  character (len=1), intent(in) :: fogtype
  character (len=10) :: fname
  character (len=150) :: clpath  ! complete path to file
  integer :: k
  real (kind=dp) :: is4(n,3)

! Common blocks:
  common /band_rat/ photol_j(nphrxn,n)
  real (kind=dp) :: photol_j
  common /blck01/ am3(n),cm3(n)
  real (kind=dp) :: am3, cm3
  common /blck11/ rc(nkc,n)
  real (kind=dp) :: rc
  common /blck12/ cw(nkc,n),cm(nkc,n)
  real (kind=dp) :: cw, cm
  common /blck13/ conv2(nkc,n) ! conversion factor = 1/(1000*cw)
  real (kind=dp) :: conv2
  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1
  common /blck78/ sa1(j2,nka)
  real (kind=dp) :: sa1
  common /budg/ bg(2,nrxn,nlev),il(nlev)
  real (kind=dp) :: bg
  integer :: il
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct
  common /kinv_i/ kinv
  integer :: kinv
  common /kpp_crys/ xcryssulf,xcrysss,xdelisulf,xdeliss
  real (kind=dp) :: xcryssulf,xcrysss,xdelisulf,xdeliss
  common /kpp_l1/ cloudt(nkc,n)
  logical :: cloudt
  common /kpp_mol/ xgamma(j6,nkc,nf)
  real (kind=dp) :: xgamma
  common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
  real (kind=dp) :: vt, vd, vdm

! == End of declarations =======================================================

  fname='rstc .dat'
  fname(5:5)=fogtype
  open (jpfunrstc,file=fname,status='unknown',form='unformatted')
! alpha, henry, vmean, xkmt are not read in because they depend on the
! chemical mechanism used (aer, tot) and are calculated every time step
! anyways

  read (jpfunrstc) &
! double precision arrays
       am3,cm,cm3,conv2,cw,es1,photol_j,rc,s1,s3,sa1, &
       sl1,sion1,vd,vdm,vt,xgamma, &
! double precision, single values
       xcryssulf,xcrysss,xdelisulf,xdeliss, &
! logicals
       cloudt, &
! integers
       il,kinv,lday,lmin,lst
  close (jpfunrstc)

! initial output for plotting - only needed for binary output
  if (binout) then
     do k=1,n
        is4(k,1)=am3(k) ! stay consistent with plot routine array size!
        is4(k,2)=0._dp
        is4(k,3)=0._dp
     enddo

     fname='sg1 .out'
     fname(4:4)=fogtype
     clpath=trim(coutdir)//trim(fname)
     open (jpfunsg1, file=trim(clpath), status='old',form='unformatted')
     write (jpfunsg1) is4
     close (jpfunsg1)

     fname='sr1 .out'
     fname(4:4)=fogtype
     clpath=trim(coutdir)//trim(fname)
     open (jpfunsr1, file=trim(clpath), status='old',form='unformatted')
     write (jpfunsr1) is4
     close (jpfunsr1)
  end if

  write(jpfunout,*)"restart file for chemistry read, filename: ",fname
  write(jpfunout,*)lday,lst,lmin
  write(jpfunout,*)'deliquescense',xdelisulf,xdeliss
  write(jpfunout,*)'crystallization',xcryssulf,xcrysss

! get alpha's
  call st_coeff_a
  call st_coeff_t

end subroutine startc

!
!-------------------------------------------------------------
!

function rgl (r_dry,a,b,feu)
!
! Description:
! -----------
  ! Equilibrium radius of aerosol particle at given relative humidity

!
! Method:
! -----------
  ! Newton iterations: solve f(x) = 0
  !   noting that f(x) ~= f(x0) + f'(x0) * (x - x0)
  !   leads to x_k+1 = x_k - f(x_k) / f'(x_k)
  !
  ! Here f(x) = (x**3 - 1) * (x*ln(rH) - A/r_dry) + Bx
  !           = x**4 ln(rH) - x ln(rH) - (x**3 - 1)*A/R_dry + Bx
  ! its derivative
  !     f'(x) = 4x**3 ln(rH) - ln(rH) - 3x**2 *A/R_dry + B
  !           = (4x**3 - 1) * ln(rH)  - 3x**2 *A/R_dry + B
  !
  !  x is the ratio between particle radius and its dry radius


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  ! 29-Apr-2021  Josue Bock  Review and merge with latest version provided by A. Bott
  !                          Remove statement function f(x) and fstr(x)
  !                          Add methods section in header

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunout

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Function arguments
  real (kind=dp), intent(in)  :: r_dry ! dry particle radius
  real (kind=dp), intent(in)  :: a     ! scaling radius for surface effect
  real (kind=dp), intent(in)  :: b     ! factor accounting for the solution effect
  real (kind=dp), intent(in)  :: feu   ! relative humidity

! Local scalars:
  integer :: ij      ! running indices
  real (kind=dp) :: rgl            ! function result = particle radius
  real (kind=dp) :: alpha
  real (kind=dp) :: falt, fstralt  ! function to solve, and its derivative
  real (kind=dp) :: xalt, xneu     ! xold, xnew or x_k, x_k+1
  real (kind=dp) :: zlogf

! == End of declarations =======================================================

  if (feu.ge.1._dp) then
     write (jpfunout,*)'  FN rgl: initial value of relative humidity exceeding one'
     rgl = r_dry
     return
  endif

  zlogf = log(feu)
  alpha = a / r_dry
! initial value
  xalt = exp(feu)
! Newton iteration
  do ij=1,100
     falt    = (xalt**3 - 1._dp)*(xalt*zlogf-alpha) + b*xalt
     fstralt = (4._dp*xalt**3 - 1._dp)*zlogf - 3._dp*xalt**2 *alpha + b
     xneu = xalt - falt/fstralt
     if (abs(xneu-xalt) .lt. 1.e-7_dp * xalt) exit
     xalt = xneu
  enddo

  rgl = r_dry*xneu

end function rgl

!
!-------------------------------------------------------------
!

subroutine sedp (dt)
!
! Description:
! -----------
  ! gravitational settling of particles with terminal velocity w in m/s
  ! For further details on the determination of w see function vterm


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  !
  ! jjb bugfix ds1 ds2, they were the wrong way
  ! jjb note the different definition of psi: ff here, while ff*detw in latest Bott version of the code
  ! jjb minor bugfix in c(2) definition, deta(k) was used instead of deta(2). No consequence since they are identical

! == End of header =============================================================

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nb, &
       nka, &
       nkt, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Scalar arguments with intent(in):
  real(kind=dp), intent(in) :: dt

! External function
  real(kind=dp), external :: vterm

! Local scalars:
  integer :: ia, jt, k                  ! running indices
  real (kind=dp) :: dt0, dtmax          ! time splitting
  real (kind=dp) :: ww                  ! sedimentation velocity
  real (kind=dp) :: x0, x1, x2, x3 !,x4 ! equation factors
  real (kind=dp) :: xsum                ! particles of one class
! Local arrays:
  real (kind=dp) :: c(nf), psi(nf)      ! Courant number and variable to be advected (formerly in cb58)

! Common blocks:
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb), &
                ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
  real (kind=dp) :: zb, dzb, dzbw, tb, eb, ak, d, &
       ajb, ajq, ajl, ajt, ajd, ajs, ds1, ds2, ajm, reif, tau, trdep
  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /blck06/ kw(nka),ka
  integer :: kw, ka
  common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
  real (kind=dp) :: vt, vd, vdm

! == End of declarations =======================================================

  ajs   = 0._dp
  c(nf) = 0._dp
  x3    =-deta(2)

  do ia=1,nka
     do jt=1,nkt
!        x4=rq(jt,ia)
!        ww=-1.25d-4*x4*x4*(1.+8.6d-02/x4)
        ww   = -1.*vterm(rq(jt,ia)*1.d-6,t(nf),p(nf)) !"first guess" for determination of ww
        dt0  = dt
        xsum = 0._dp
        do k=2,nf
           psi(k) = ff(jt,ia,k) * detw(k)
           xsum   = xsum + psi(k)
        enddo

! restrict sedimentation to classes with significant particle content
        if (xsum.gt.1.e-6_dp) then
           x0 = 0._dp

! time step control
! multiple calls of advection scheme for large courant numbers
           do while (dt0 .gt. 0.1_dp)
!              dtmax=min(dt0,x3/(ww+w(nf)))
! see also SR difp: subsidence treated now consistently (i.e. like for other
! tracers): w df/dz instead of d(wf)/dz
              dtmax=min(dt0,x3/(ww))

              do k=2,nf
!                 c(k)=dtmax/deta(k)*(ww+w(k))
!                 c(k)=dtmax/deta(k)*(-1.*vterm(rq(jt,ia)*1.d-6, t(k), p(k))+w(k))
                 c(k)=dtmax/deta(k)*(-1.*vterm(rq(jt,ia)*1.d-6, t(k), p(k)))
              enddo
! particle dry deposition velocity in lowest model layer:
              c(2)=min(c(2),dtmax/deta(2)*vd(jt,ia)*(-1.))
              c(1)   = c(2)
              dt0    = dt0 - dtmax
              x1     = psi(2)
              psi(1) = x1

! vertical advection of ff(ia, jt)
! simplified upstream advection for small particles
              if (rq(jt,ia) .lt. 1._dp) then
                 call advsed0(c, psi)
              else
                 call advsed1(c, psi)
              endif
              x0 = x0 + psi(1) - x1
           end do

! new values of ff
           do k=2,nf-1
              ff(jt,ia,k) = psi(k) / detw(k)
           enddo
           ff(jt,ia,nf) = ff(jt,ia,nf-1)
        endif

! Droplet sedimentation has been evaluated at ground
! trdep :  cumulative deposition [kg liquid water/m^2]
!          being equivalent to [mm] precipitation
! ajs   :  sedimentation rate [kg liquid water/m^2/sec]
!          being equivalent to precip. rate of [mm/sec]
! update total liquid water [kg/m^3]
        x2    = x0 * e(jt) * detw(2)
        ajs   = ajs + x2 / dt
        trdep = trdep + x2

! division of sedimentation into small aerosol and large droplets
! for output
        if (jt.le.kw(ia)) then
           ds1=ds1+x2
        else
           ds2=ds2+x2
        endif
     enddo
  enddo

end subroutine sedp

!
!-------------------------------------------------------------
!

subroutine sedc (dt)
! dry deposition and emission of gaseous species

! jjb work done = implicit none, missing declarations, little cleaning, modules including constants

  USE config, ONLY : &
! Imported Parameters:
       lpJoyce14bc

  USE constants, ONLY : &
! Imported Parameters:
       Avogadro

  USE gas_common, ONLY: &
! Imported Parameters:
       j1, &
! Imported Array Variables with intent (in):
       es1, &
       ind_gas_rev, &
! Imported Array Variables with intent (inout):
       s1, &
       vg

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) ::  dt

! Local scalars:
  integer :: j
  real (kind=dp) :: s12old
  real (kind=dp) :: x4, w

! Common blocks:
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

! == End of declarations =======================================================


! constant dry deposition velocities in m/sec, not for all species defined
! old data set:
!      data vg/0.10e-2,0.10e-2,0.10e-1,0.1e-1,0.8e-02,
!c     &      0.01e-2,0.50e-2,7*0.10e-2,2*0.0,0.30e-2,
!     &      0.01e-2,0.50e-3,7*0.10e-2,2*0.0,0.30e-2,
!     &      3*0.2e-2,0.50e-2,2*0.20e-02,2*0.10e-2,
!c     &      3*0.20e-2,0.10e-1,0.,3*0.01e-2,3*0.10e-1,
!c     &      4*0.0/
!     &      3*0.20e-2,0.10e-1,0.2e-1,3*0.01e-2,3*0.10e-1,
!     &      0.1e-2,0.2e-3,3*0.0,0.2e-1,0.2e-2,4*0.,
!     &      6*0.1e-1,2*0.0,0.1e-1,9*0.0/

! Laurens data set:
!x      data vg/0.0,0.10e-4,0.050e-1,0.27e-2,0.5e-02,
!x     &      1.0e-2,0.40e-3,7*0.0,0.5e-2,0.0,0.30e-2,
!x     &      0.2e-2,0.4e-2,0.2e-2,0.0,2*0.2e-02,2*0.10e-2,
!x     &      3*0.20e-2,0.10e-1,0.2e-1,3*0.01e-2,0.10e-1,0.0,0.0,
!x     &      0.0,0.2e-2,3*0.0,0.2e-1,0.2e-2,4*0.,
!x     &      5*0.1e-1,4*0.0,3*0.0,1.e-2,2*0.0,1.e-2,2*0.0/

!      vg(4)=0.27e-2 ! NH3 old value, that fitted "nicely" in model    !=0. ! emission is net flux or
!      vg(34)=vg(30) ! N2O5=HCl
!      vg(37)=0.     ! DMS           emission is net flux
!      vg(38)=vg(30) ! HOCl = HCl
!      vg(43)=vg(30) ! HOBr = HCl
!      vg(50)=vg(49) ! I2O2=HOI
!      vg(51)=vg(49) ! INO2=HOI
!      vg(56)=0.     ! CH3I          emission is net flux
!      vg(57)=0.     ! CH2I2         emission is net flux
!      vg(58)=0.     ! CH2ClI        emission is net flux
!      vg(59)=0.     ! C3H7I         emission is net flux
!      vg(63)=vg(30) ! CH3SO3H = HCl
!      vg(71)=0.     ! CH2BrI         emission is net flux
!      vg(72)=0.     ! CHBr2I         emission is net flux
!      vg(73)=0.     ! C2H5I          emission is net flux

  if (.not.lpJoyce14bc) then
     if(ind_gas_rev(4) /= 0) &
          vg(ind_gas_rev(4))=0.27e-2_dp              ! NH3 old value, that fitted "nicely" in model    !=0. ! emission is net flux or
  end if

  ! N2O5
  ! special case for Joyce et al 2014 study
  if (lpJoyce14bc) then
     if(ind_gas_rev(34) /= 0) &
          vg(ind_gas_rev(34))=5.9e-3_dp            ! N2O5, median field value, (Huff et al, 2011)  PJ
  ! general case
  else
     if(ind_gas_rev(34) /= 0 .and. ind_gas_rev(30) /= 0) &
          vg(ind_gas_rev(34))=vg(ind_gas_rev(30))  ! N2O5=HCl
  end if

  if(ind_gas_rev(37) /= 0) &
       vg(ind_gas_rev(37))=0._dp                  ! DMS           emission is net flux
  if(ind_gas_rev(38) /= 0 .and. ind_gas_rev(30) /= 0) &
       vg(ind_gas_rev(38))=vg(ind_gas_rev(30)) ! HOCl = HCl
  if(ind_gas_rev(43) /= 0 .and. ind_gas_rev(30) /= 0) &
       vg(ind_gas_rev(43))=vg(ind_gas_rev(30)) ! HOBr = HCl
  if(ind_gas_rev(50) /= 0 .and. ind_gas_rev(49) /= 0) &
       vg(ind_gas_rev(50))=vg(ind_gas_rev(49)) ! I2O2=HOI
  if(ind_gas_rev(51) /= 0 .and. ind_gas_rev(49) /= 0) &
       vg(ind_gas_rev(51))=vg(ind_gas_rev(49)) ! INO2=HOI
  if(ind_gas_rev(56) /= 0) &
       vg(ind_gas_rev(56))=0._dp                  ! CH3I          emission is net flux
  if(ind_gas_rev(57) /= 0) &
       vg(ind_gas_rev(57))=0._dp                  ! CH2I2         emission is net flux
  if(ind_gas_rev(58) /= 0) &
       vg(ind_gas_rev(58))=0._dp                  ! CH2ClI        emission is net flux
  if(ind_gas_rev(59) /= 0) &
       vg(ind_gas_rev(59))=0._dp                  ! C3H7I         emission is net flux
  if(ind_gas_rev(63) /= 0 .and. ind_gas_rev(30) /= 0) &
       vg(ind_gas_rev(63))=vg(ind_gas_rev(30)) ! CH3SO3H = HCl
  if(ind_gas_rev(71) /= 0) &
       vg(ind_gas_rev(71))=0._dp                  ! CH2BrI         emission is net flux
  if(ind_gas_rev(72) /= 0) &
       vg(ind_gas_rev(72))=0._dp                  ! CHBr2I         emission is net flux
  if(ind_gas_rev(73) /= 0) &
       vg(ind_gas_rev(73))=0._dp                  ! C2H5I          emission is net flux

  if (lst/4*4.eq.lst.and.lmin.eq.1) then
     print *,lday,lst,lmin
     print *,' dry deposition velocities'
     do j=1,j1
        print *,j,vg(j)
     enddo
  endif


!      x3=deta(2)
! emission rates: 80% of total emission during day 20% during night
! 1.6e-2=1.6/100: factor 100 due to detw(2) with units in cm
!      x0=float(lmin)/60.
! emission rates of chemical species variable in time
!      x4=1.6e-02*dt
!      if (lst.ge.19.or.lst.lt.6) x4=.4e-02*dt
!      if (lst.eq.6) x4=(1.6e-02*x0+0.4e-02*(1.-x0))*dt
!      if (lst.eq.18) x4=(0.4e-02*x0+1.6e-02*(1.-x0))*dt
  x4=1._dp
!      dt0=dt
  do j=1,j1
     w=vg(j)
     if (w.ge.1.e-5_dp) then
!            psi2=s1(j,2)*detw(2)
!            psi3=s1(j,3)*detw(3)
!            dtmax=dmin1(dt0,x3/w)
!            cl=dtmax*w/x3
!            a0=(25.*psi2-psi3)/24.
!            a1=(psi3-psi2)/16.
!            a2=(psi3-psi2)/48.
!            x1=1.-2.*cl
!            x2=x1*x1
!            xxx2=a0*cl-a1*(1.0-x2)+a2*(1.0-x1*x2)
!            fm=dmax1(0.d0,xxx2)
!            flux1=fm*psi2/dmax1(fm+1.e-15,psi2)
!            s1(j,2)=(psi2-flux1)/detw(2)
!c total dry deposition of gas phase phase species in mol/m^2        : depfac=1.
!            depfac=1.
!c            s1(j,1)=s1(j,1)+flux1*100.
!            s1(j,1)=s1(j,1)+flux1*depfac
        s12old=s1(j,2)
        s1(j,2)=s1(j,2)*exp(-dt/deta(2)*vg(j))
        s1(j,1)=s1(j,1)+(s12old-s1(j,2))*deta(2)
     endif
! es1: emission rates in molec./cm**2/s, s1 in mol/m**3
     s1(j,2)=s1(j,2)+es1(j)*x4*dt*1.e+4_dp/(detw(2)*Avogadro)
  enddo

! Scenario emission for special model study
  if (lpJoyce14bc) then
!  surface gas emission  "gasPJ"
!  if(lst.le.2) then
!  if(lst.ge.2.AND.lst.le.3) then
     if(lday.eq.0.AND.lst.ge.2.AND.lst.le.3) then
!     es1(1)=1.0d12       ! NO
        es1(ind_gas_rev(1))=1.4e12_dp        ! NO (b25)
!     es1(4)=4.2e10       ! NH3 (base23)
!     es1(4)=3.4e11       ! NH3
!     es1(4)=5.2d10       ! NH3 (3.7% of NO, based on inventory)
        es1(ind_gas_rev(4))=6.7e10_dp        ! NH3 (4.76% of NO, Yokelson 97)
!     es1(4)=4.2d11       ! NH3
!     es1(4)=6.3d11       ! NH3 (b25)
!     es1(4)=1.0d12       ! NH3
!     es1(5)=1.5d11       ! SO2 (Dick,2003)
        es1(ind_gas_rev(5))=3.0e11_dp        ! SO2 (DEC, 2008)
     else
        es1(ind_gas_rev(1))=0._dp
        es1(ind_gas_rev(4))=0._dp
        es1(ind_gas_rev(5))=0._dp
     endif
  end if

end subroutine sedc

!
!----------------------------------------------------------------
!

subroutine sedl (dt)
! new aqueous phase concentrations due to
! gravitational settling of droplets


! Author:
! ------
  !    Andreas Bott, Roland von Glasow


! Modifications :
! -------------
  ! jjb minor bugfix in c(2) definition, deta(k) was used instead of deta(2). No consequence since they are identical
  !

! == End of header =============================================================

  USE config, ONLY : &
       nkc_l

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j6, &
       nf, &
       n, &
       nka, &
       nkt, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Scalar arguments with intent(in):
  real(kind=dp), intent(in) :: dt

! External function
  real(kind=dp), external :: vterm

! Local scalars:
  integer :: k, kc, l
  real (kind=dp) :: dt0, dtmax
  real (kind=dp) :: xfac, xxx, xxxt, x0, x1, x4
! Local arrays:
  real (kind=dp) :: c(nf), psi(nf) ! Courant number and variable to be advected (formerly in cb58)
  real (kind=dp) :: cc(nf)

! Common blocks:
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /blck11/ rc(nkc,n)
  real (kind=dp) :: rc
  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
  real (kind=dp) :: vt, vd, vdm

! == End of declarations =======================================================

  c(nf)=0._dp
! changes have to be made BOTH here and below for ions
! rc in m
  xfac=1.e6_dp
  do kc=1,nkc_l
     do k=2,nf
        xxx=0.01_dp ! jjb apparently, in um
        x4=max(xxx,xfac*rc(kc,k)) ! jjb here as well
! subsidence see SR difl
!        cc(k)=(-1.25e-4*x4*x4*(1.+8.6e-02/x4))/deta(k)
        cc(k)=(-1._dp*vterm(x4*1.e-6_dp,t(k),p(k)))/deta(k) ! jjb: in vterm, radius in m
! mass weighted terminal velocity
        cc(k)=min(cc(k),-1._dp*vt(kc,k)/deta(k))
     enddo
! particle dry deposition velocity in lowest model layer:
     cc(2)=min(cc(2),-1._dp/deta(2)*vdm(kc))
     do l=1,j2
        do k=2,nf
           psi(k)=sl1(l,kc,k) * detw(k)
        enddo
        dt0  = dt
        x0   = 0._dp
        xxxt =-.999_dp / cc(2) ! jjb difference with sedp: cc(2) is used instead of ww (computed with t(nf) and p(nf) ) and here -.999 instead of -deta(2)

! time step control
! multiple calls of advection scheme for large courant numbers
        do while (dt0 .gt. 0.1_dp)
           dtmax=min(dt0,xxxt)
           dt0=dt0-dtmax
           do k=2,nf
              c(k)=cc(k)*dtmax
           enddo
           ! jjb maybe missing here: c(2) = min(c(2),dtmax/deta(2)*vdm(kc)*(-1.))
           c(1)   = c(2)
           x1     = psi(2)
           psi(1) = x1
           call advsed1(c,psi)
           x0 = x0 + psi(1) - x1
        end do

! new values of sl1
        do k=2,nf-1
           sl1(l,kc,k) = psi(k) / detw(k)
        enddo
! wet deposition to the ground in mole/m**2
        sl1(l,kc,1) = sl1(l,kc,1) + x0 * deta(2)
     enddo
  enddo

! dito for ions
  c(nf)=0._dp
  do kc=1,nkc_l
     do k=2,nf
        xxx=0.01_dp
        x4=max(xxx,xfac*rc(kc,k))
! subsidence see SR difl
!        cc(k)=(-1.25e-4*x4*x4*(1.+8.6e-02/x4))/deta(k)
        cc(k)=(-1._dp*vterm(x4*1.e-6_dp,t(k),p(k)))/deta(k)
! mass weighted terminal velocity
        cc(k)=min(cc(k),-1._dp*vt(kc,k)/deta(k))
     enddo
! particle dry deposition velocity in lowest model layer:
     cc(2)=min(cc(2),-1._dp/deta(2)*vdm(kc))
     do l=1,j6
        do k=2,nf
           psi(k)=sion1(l,kc,k) * detw(k)
        enddo
        dt0  = dt
        x0   = 0._dp
        xxxt = -.999_dp / cc(2)

! time step control
! multiple calls of advection scheme for large courant numbers
        do while (dt0 .gt. 0.1_dp)
           dtmax=min(dt0,xxxt)
           dt0=dt0-dtmax
           do k=2,nf
              c(k)=cc(k)*dtmax
           enddo
           c(1)   = c(2)
           x1     = psi(2)
           psi(1) = x1
           call advsed1(c,psi)
           x0 = x0 + psi(1) - x1
        enddo

! new values of sion1
        do k=2,nf-1
           sion1(l,kc,k) = psi(k) / detw(k)
        enddo
! wet deposition to the ground in mole/m**2
        sion1(l,kc,1) = sion1(l,kc,1) + x0 * deta(2)
     enddo
  enddo !ions

end subroutine sedl

!
!------------------------------------------------------------
!

      function vterm(a,t,p)

! terminal velocity for droplets
! after Stokes with Cunningham correction for small droplets (regime 1)
! and after Beard for large droplets (r > 10 micron, regime 2)
! all formulas after Pruppacher and Klett Chapter 10.

      USE constants, ONLY : &
! Imported Parameters:
     &     g, &
     &     r0, &                   ! Specific gas constant of dry air, in J/(kg.K)
     &     rhow                  ! Water density [kg/m**3]

  USE precision, ONLY : &
! Imported Parameters:
       dp

      implicit none

      real (kind=dp) :: vterm

      real (kind=dp), intent(in) :: a ! radius       in [m]
      real (kind=dp), intent(in) :: t ! temperature  in [K]
      real (kind=dp), intent(in) :: p ! pressure     in [Pa]

      ! Polynomial coefficients for Beard approximation
      ! Pruppacher & Klett, p. 417, equation (10-145)
      real (kind=dp), parameter :: b0=-.318657d+1
      real (kind=dp), parameter :: b1= .992696d+0
      real (kind=dp), parameter :: b2=-.153193d-2
      real (kind=dp), parameter :: b3=-.987059d-3
      real (kind=dp), parameter :: b4=-.578878d-3
      real (kind=dp), parameter :: b5=+.855176d-4
      real (kind=dp), parameter :: b6=-.327815d-5

      real (kind=dp), parameter :: c1 = 2.d0 * g / 9.d0       ! constant factor in equation (10-138)
      real (kind=dp), parameter :: c2 = 1.26d0                ! constant in equation (10-139)
      real (kind=dp), parameter :: P0 = 101325                ! Standard pressure    [Pa]
      real (kind=dp), parameter :: T0 = 293.15_dp             ! Standard temperature [K]
      real (kind=dp), parameter :: lambda0 = 6.6d-8           ! mean free path STP   [m]
      real (kind=dp), parameter :: c3 = c2 * lambda0 * P0/T0  ! constant factor in equation (10-139) and (10-140)
      real (kind=dp), parameter :: c4 = 32.d0 * g / 3.d0      ! constant factor in equation (10-142)

      real (kind=dp) :: best  ! Davies or Best number, see equation (10-142) [-]
      real (kind=dp) :: x     ! ln(best)                                     [-]
      real (kind=dp) :: rho_a ! air density                                  [kg/m3]
      real (kind=dp) :: eta   ! dynamic viscosity                            [kg/(m.s)]
      real (kind=dp) :: y     ! Reynolds number = exp(y) in Regime 2         [-]

      rho_a=p/(r0*t)
      eta=3.7957d-06+4.9d-08*t

      if (a.le.1.d-5) then
         ! Regime 1: see P & K pp. 415-417, equations (10-138) to (10-140)
         vterm=c1*a*a*(rhow-rho_a)/eta*(1._dp+c3*t/(a*p))
         ! 2.17... = 2*g/9, see (10-138)
         ! 3.0849d-5 = 1.26 * lamda_0 * P_0 / T_0 with a mistake over T_0: 273.15 instead of 293.15 in P & K
      else
         ! Regime 2
         best=c4*a**3*(rhow-rho_a)*rho_a/(eta*eta) ! equation (10-142) with 104.60267 = 32*g/3
         x=log(best)
         ! evaluation of BEARD-polynomial with Horner-scheme
         y=b6*x+b5
         y=y*x+b4
         y=y*x+b3
         y=y*x+b2
         y=y*x+b1
         y=y*x+b0
         vterm=eta*exp(y)/(2.*rho_a*a)
      end if

      end function vterm

!
!----------------------------------------------------------------
!

subroutine wfield
!
! Description:
! -----------
  ! Optional subroutine: calculation of subsidence if chosen to be time dependent


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  ! 29-Apr-2021  Josue Bock  Review: warning, this is not the same in latest Mistra version from A. Bott

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! Imported Parameters:
       wmin, wmax                  ! large scale subsidence

  USE constants, ONLY : &
! Imported Parameters:
       pi

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  integer :: k           ! running indices
  real (kind=dp) :: u0   ! solar angle
  real (kind=dp) :: zeit ! time [s]

! Common blocks:
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)               ! eta: level height
  real (kind=dp) :: detw, deta, eta, etw

  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w                                 ! w: subsidence

! == End of declarations =======================================================

  zeit = lst*3600._dp + lmin*60._dp
  u0 = cos(2._dp*pi*zeit/86400._dp)
  do k=1,n
     w(k)=eta(k)/1000._dp*0.5_dp*((wmax+wmin)+(wmin-wmax)*u0)
  enddo
  do k=n,1,-1
     w(k)=w(k)-w(1)
  enddo

end subroutine wfield

!
!-------------------------------------------------------------------
!

subroutine difm (dt)
!
! Description:
! -----------
  ! fully implicit procedure for the solution of the diffusion equations
  ! after Roache, 1972: Computational fluid dynamics, Appendix A.
  ! all quantities are similarly defined with a(k) --> xa(k), etc.
  ! except xd(k) which has another meaning than d(k) of Roache's program.
  ! dirichlet conditions at the surface and at the top of the atmosphere


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  !

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:
  USE config, ONLY : &
! Imported Parameters:
       ug, vg                      ! geostrophic wind

  USE constants, ONLY : &
! Imported Parameters:
       r0               ! Specific gas constant of dry air, in J/(kg.K)

  USE data_surface, ONLY : &
       ustern                      ! frictional velocity

  USE global_params, ONLY : &
! Imported Parameters:
       n,                   &
       nm

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
  real (kind=dp), intent(in)  :: dt          ! fractional timestep

! Local parameters:
  real (kind=dp), parameter :: fcor=1.e-4_dp ! mean coriolis parameter (1/s)

! Local scalars:
  integer :: k, kp      ! running indices
  real (kind=dp) :: fdt
  real (kind=dp), external :: p21

! Local arrays:
  real (kind=dp) :: c(n)
  real (kind=dp) :: xa(nm),xb(nm),xc(nm),xd(nm),xe(nm),xf(nm),oldu(nm)

! Common blocks:
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy

  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb53a/ thet(n),theti(n)
  real(kind=dp) :: thet, theti
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a

! == End of declarations =======================================================

! update of prognostic variables
  tke(1)=max(1.e-06_dp,3.2537_dp*ustern**2)
  do k=1,n
     rho(k)=p(k)/(r0*(t(k)*(1._dp+0.61_dp*xm1(k))))
     theta(k)=t(k)*thet(k)
     tke(k)=max(1.e-05_dp,tke(k)+tkep(k)*dt)
     c(k)=w(k)*dt/deta(k)
  enddo

! calculation of turbulent exchange coefficients
  call atk1

! -----------------------------------------------------------------------------

! solution of the diffusive equations
!------------------------------------

! turbulent exchange with 'k_m' (atkm):
  xa(1)=atkm(1)*dt/(detw(1)*deta(1))
  xe(1)=0._dp
  do k=2,nm
     xa(k)=atkm(k)*dt/(detw(k)*deta(k))
     xc(k)=xa(k-1)*detw(k-1)/detw(k)
     xb(k)=1._dp+xa(k)+xc(k)
     xd(k)=xb(k)-xc(k)*xe(k-1)
     xe(k)=xa(k)/xd(k)
  enddo

! u-component of horizontal wind field
  fdt=fcor*dt
  xf(1)=0._dp
  do k=2,nm
     oldu(k)=u(k) ! save u for v calculation below
     xf(k)=(u(k)+fdt*(v(k)-vg)+xc(k)*xf(k-1))/xd(k)
  enddo
  do k=nm,2,-1
     u(k)=xe(k)*u(k+1)+xf(k)
  enddo

! v-component of the horizontal wind field
  do k=2,nm
     xf(k)=(v(k)-fdt*(oldu(k)-ug)+xc(k)*xf(k-1))/xd(k)
  enddo
  do k=nm,2,-1
     v(k)=xe(k)*v(k+1)+xf(k)
  enddo

! turbulent kinetic energy with 'k_e' (atke):
  xa(1)=atke(1)*dt/(detw(1)*deta(1))
  !xe(1)=0._dp
  do k=2,nm
     xa(k)=atke(k)*dt/(detw(k)*deta(k))
     xc(k)=xa(k-1)*detw(k-1)/detw(k)
     xb(k)=1._dp+xa(k)+xc(k)
     xd(k)=xb(k)-xc(k)*xe(k-1)
     xe(k)=xa(k)/xd(k)
  enddo
  xf(1)=tke(1)
  do k=2,nm
     xf(k)=(tke(k)+xc(k)*xf(k-1))/xd(k)
  enddo
  do k=nm,2,-1
     tke(k)=xe(k)*tke(k+1)+xf(k)
  enddo

! turbulent exchange with 'k_h' (atkh):
  xa(1)=atkh(1)*dt/(detw(1)*deta(1))
  !xe(1)=0._dp
  do k=2,nm
     xa(k)=atkh(k)*dt/(detw(k)*deta(k))
     xc(k)=xa(k-1)*detw(k-1)/detw(k)
     xb(k)=1._dp+xa(k)+xc(k)
     xd(k)=xb(k)-xc(k)*xe(k-1)
     xe(k)=xa(k)/xd(k)
  enddo

! specific humidity
  xf(1)=xm1(1)
  do k=2,nm
     xf(k)=(xm1(k)+xc(k)*xf(k-1))/xd(k)
  enddo
  do k=nm,2,-1
     xm1(k)=xe(k)*xm1(k+1)+xf(k)
  enddo
! new temperature obtained from potential temperature
  xf(1)=theta(1)
  do k=2,nm
     xf(k)=(theta(k)+xc(k)*xf(k-1))/xd(k)
  enddo
  do k=nm,2,-1
     theta(k)=xe(k)*theta(k+1)+xf(k)
  enddo

! large scale subsidence
  do k=2,nm
     kp=k+1
     theta(k)=theta(k)-c(k)*(theta(kp)-theta(k))
     u(k)=u(k)-c(k)*(u(kp)-u(k))
     v(k)=v(k)-c(k)*(v(kp)-v(k))
     xm1(k)=xm1(k)-c(k)*(xm1(kp)-xm1(k))
     tke(k)=tke(k)-0.5_dp*(c(k)+c(k+1))*(tke(kp)-tke(k))
  enddo
  do k=2,nm
     t(k)=theta(k)*theti(k)
     feu(k)=xm1(k)*p(k)/((0.62198_dp+0.37802_dp*xm1(k))*p21(t(k)))
  enddo

end subroutine difm

!
!----------------------------------------------------------------------
!

subroutine difp (dt)
!
! Description:
! -----------
  ! Turbulent diffusion of particles.
  ! fully implicit procedure for the solution of the turbulent transport
  ! of aerosols and cloud droplets. For further details see subroutine difm
  !
  ! for diffusion mixing ratio is needed --> factor 1./am3(k,1)
  ! (#/cm^3 --> #/mol)


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  ! 07-May-2021  Josue Bock  several changes in the units conversion, am3 undefined in restart
  !                          There was actually a bug if chem = .false.
  !                          First redefined a local am3 = rho/M_air, with a further 1e6 factor
  !                            in the conversion (thus acm3 = am3 * 1e6)
  !                          After tests, using rho instead of acm3 gives exactly the same result
  !                          But am3(k-1)/am3(k) has also been introduced in xc coefficient
  !                          It seems wrong, thus removed
  !
  ! 08-May-2021  Josue Bock  Turn xf(nm) into xf(nkt,nka,n) for matrix operation.
  !                          On a test case, it was 10 times faster this way

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:
  USE global_params, ONLY : &
! Imported Parameters:
       n,                   &
       nm,                  &
       nka,                 &
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
  real (kind=dp), intent(in)  :: dt          ! fractional timestep

! Local scalars:
  integer :: k, kp, ia, jt      ! running indices

! Local arrays:
  real (kind=dp) :: c(n)
  real (kind=dp) :: xa(nm),xb(nm),xc(nm),xd(nm),xe(nm)
  real (kind=dp) :: xf(nkt,nka,nm)

! Common blocks:
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy

  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w
  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho

! == End of declarations =======================================================

! solution of the diffusive equation for particles
!-------------------------------------------------

! mass specific particle conversion
! constant particle distribution in layer nm+1 is used to update ff(nm)
  do k=2,nm+1
     ff(:,:,k) = ff(:,:,k) / rho(k)
  end do

! turbulent exchange of aerosol and droplets with 'k_h' (atkh)
  xa(1) = atkh(1) * dt / (detw(1) * deta(1))
  xe(1) = 0._dp
!  do k=2,nf
  do k=2,nm
     xa(k) = atkh(k) * dt / (detw(k) * deta(k))
     xc(k) = xa(k-1) * detw(k-1) / detw(k)
     xb(k) = 1._dp + xa(k) + xc(k)
     xd(k) = xb(k) - xc(k) * xe(k-1)
     xe(k) = xa(k) / xd(k)
     c(k)  = w(k) * dt / deta(k)
  enddo

  xf(:,:,1) = ff(:,:,2)
  do k=2,nm
     xf(:,:,k) = (ff(:,:,k) + xc(k) * xf(:,:,k-1)) / xd(k)
  end do
  do k=nm,2,-1
     ff(:,:,k) = xe(k) * ff(:,:,k+1) + xf(:,:,k)
  enddo
! conversion back to 1/cm^3
  do k=2,nm+1
     ff(:,:,k) = ff(:,:,k) * rho(k)
  end do

! large scale subsidence
  do k=2,nm
     kp=k+1
     ff(:,:,k) = ff(:,:,k) - c(k) * (ff(:,:,kp) - ff(:,:,k))
  enddo

! update of fsum
!  do k=2,nf
  do k=2,n
     fsum(k) = 0._dp
     do ia=1,nka
        do jt=1,nkt
           fsum(k) = fsum(k) + ff(jt,ia,k)
        enddo
     enddo
  enddo

end subroutine difp

!
!-----------------------------------------------------------------
!

subroutine difc (dt)
! fully implicit procedure for the solution of the turbulent transport
! of chemical species
! for further details see subroutine difm


! Author:
! ------
  !    Andreas Bott, Roland von Glasow


! Modifications :
! -------------
  !
! jjb work done: removal of unused arguments
!     missing declarations and implicit none

! == End of header =============================================================

  USE config, ONLY : &
       nkc_l

  USE gas_common, ONLY : &
! Imported Parameters:
       j1,               &
       j5,               &
! Imported Array Variables with intent (inout):
       s1,               &
       s3

  USE global_params, ONLY : &
! Imported Parameters:
       j2,                  &
       j6,                  &
       n,                   &
       nm,                  &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: dt

! Local scalars:
  integer :: k, kc, kp, j
! Local arrays:
  real (kind=dp) :: c(n)
  real (kind=dp) :: xa(nm),xb(nm),xc(nm),xd(nm),xe(nm),xf(nm)

! Common blocks:
  common /blck01/ am3(n),cm3(n)
  real (kind=dp) :: am3, cm3

  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy

  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w

! == End of declarations =======================================================


! calculation of exchange coefficients
  xa(1)=atkh(1)*dt/(detw(1)*deta(1))
  xe(1)=0._dp
  do k=2,nm
     xa(k)=atkh(k)*dt/(detw(k)*deta(k))
     xc(k)=xa(k-1)*detw(k-1)/detw(k)
     xb(k)=1._dp+xa(k)+xc(k)
     xd(k)=xb(k)-xc(k)*xe(k-1)
     xe(k)=xa(k)/xd(k)
     c(k)=w(k)*dt/deta(k) ! large scale subsidence
  enddo

! gase phase species
! s1, s3, sl1, sion1 in mol/m^3
! for diffusion mixing ratio is needed --> factor 1./am3(k)
! (mol/m^3 --> mol/mol)
  do j=1,j1
     xf(1)=s1(j,2)/am3(2)
     do k=2,nm
        xf(k)=(s1(j,k)/am3(k)+xc(k)*xf(k-1))/xd(k)
     enddo
     do k=nm,2,-1
        s1(j,k)=(xe(k)*s1(j,k+1)/am3(k+1)+xf(k))*am3(k)
     enddo
!    large scale subsidence
     do k=2,nm
        s1(j,k)=s1(j,k)-c(k)*(s1(j,k+1)-s1(j,k))
     enddo
  enddo

! dito for radicals
  do j=1,j5
     xf(1)=s3(j,2)/am3(2)
     do k=2,nm
        xf(k)=(s3(j,k)/am3(k)+xc(k)*xf(k-1))/xd(k)
     enddo
     do k=nm,2,-1
        s3(j,k)=(xe(k)*s3(j,k+1)/am3(k+1)+xf(k))*am3(k)
     enddo
!    large scale subsidence
     do k=2,nm
        s3(j,k)=s3(j,k)-c(k)*(s3(j,k+1)-s3(j,k))
     enddo
  enddo

! aqueous phase species
!      if (lct.lt.2) return
!      if (ndt.lt.2) return

  do kc=1,nkc_l
     do j=1,j2
        xf(1)=sl1(j,kc,2)/am3(2)
!        do k=2,nf
        do k=2,nm
!            do k=lcl,lct
!            xf(ndb-1)=sl1(j,kc,ndb)
!            xf(ndt+1)=sl1(j,kc,ndt)
!            do k=ndb,ndt
           xf(k)=(sl1(j,kc,k)/am3(k)+xc(k)*xf(k-1))/xd(k)
        enddo
!        do k=nf,2,-1
        do k=nm,2,-1
!            do k=lct,lcl,-1
!            do k=ndt,ndb,-1
           sl1(j,kc,k)=(xe(k)*sl1(j,kc,k+1)/am3(k+1)+xf(k))*am3(k)
        enddo
!       large scale subsidence
!            do k=2,nf
        do k=2,nm
           kp=k+1
           sl1(j,kc,k)=sl1(j,kc,k)-c(k)*(sl1(j,kc,kp)-sl1(j,kc,k))
        enddo
     enddo
  enddo
! aqueous phase ion species
  do kc=1,nkc_l
     do j=1,j6
        xf(1)=sion1(j,kc,2)/am3(2)
!        do k=2,nf
        do k=2,nm
!            do k=lcl,lct
!            xf(ndb-1)=sion1(j,ndb,kc)
!            xf(ndt+1)=sion1(j,ndt,kc)
!            do k=ndb,ndt
           xf(k)=(sion1(j,kc,k)/am3(k)+xc(k)*xf(k-1))/xd(k)
        enddo
!        do k=nf,2,-1
        do k=nm,2,-1
!            do k=lct,lcl,-1
!            do k=ndt,ndb,-1
           sion1(j,kc,k)=(xe(k)*sion1(j,kc,k+1)/am3(k+1)+xf(k))*am3(k)
        enddo
!       large scale subsidence
!        do k=2,nf
        do k=2,nm
           kp=k+1
           sion1(j,kc,k)=sion1(j,kc,k)-c(k)*(sion1(j,kc,kp)-sion1(j,kc,k))
        enddo
     enddo
  enddo

end subroutine difc

!
!--------------------------------------------------------------------
!

subroutine atk0
! calculation of exchange coefficients, mixing length etc at model start


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  !

! == End of header =============================================================

  USE config, ONLY : &
! Imported Parameters:
       ug, vg                      ! geostrophic wind

  USE constants, ONLY : &
! Imported Parameters:
       g, kappa

  USE data_surface, ONLY : &
       ustern, z0, &               ! frictional velocity, roughness length
       gclu, gclt                  ! coefficients for momentum, and temperature and humidity

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  integer :: k
  real (kind=dp) :: x0, x1, x2
  real (kind=dp) :: st, vh, zz

! Common blocks:
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw
  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy
  common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
  real (kind=dp) :: gm, gh, sm, sh, xl
  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho

! == End of declarations =======================================================

! initialisation of mixing length xl
  xl(1) = 0._dp
  x1 = (ug + vg) * 2.7_dp
  do k=2,n
     x2    = kappa * etw(k)
     xl(k) = x2 * x1 / (x2 + x1)
     xl(k) = min(xl(k), deta(k))
  enddo

! initial calculation of exchange coefficients
  atkm(1) = 0.5_dp * eta(2) * ustern / gclu
  atkh(1) = 0.5_dp * eta(2) * ustern / gclt
  do k=2,n-1
     vh = ((u(k+1)-u(k))**2 + (v(k+1)-v(k))**2) / deta(k)**2
     zz = etw(k) + z0
     x0 = (0.4_dp * zz / (1._dp + 0.4_dp * zz / xl(k)))**2
     st = g * (theta(k+1)-theta(k)) / (deta(k)*theta(k))
     if (st.le.0._dp) then
        ! unstable case and neutral case
        atkm(k) = x0 * sqrt(vh - 11._dp * st)
        if ((vh-3._dp*st).eq.0._dp) then
           atkh(k) = atkm(k)
        else
           atkh(k) = 1.35_dp * atkm(k) * (vh - 5.5_dp * st) / (vh - 3._dp * st)
        endif
     else
        ! stable case
        atkm(k) = x0 * vh / sqrt(vh + 6._dp * st)
        atkh(k) = 1.35_dp * atkm(k) * vh / (vh + 6._dp * st)
     endif
     atkm(k) = max(1.e-3_dp, atkm(k))
     atkh(k) = max(1.e-3_dp, atkh(k))
  enddo
  atkm(n) = 0._dp
  atkh(n) = 0._dp

end subroutine atk0

!
!-------------------------------------------------------------
!

subroutine atk1
!
! Description:
! -----------
  ! Calculation of the tke production (tkep) and of the turbulent
  ! exchange coefficients for heat (atkh), moisture (atkm), and
  ! kinetic energy (atke) after 2.5 level model of Mellor and Yamada
  ! (Reviews of Geophysics and Space Physics (4), 1982, pp 851-875).



! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  !

! == End of header =============================================================

  USE constants, ONLY : &
! Imported Parameters:
       g, &                        ! Gravitational acceleration
       kappa

  USE data_surface, ONLY : &
       ustern, &                   ! frictional velocity
       gclu, gclt                  ! coefficients for momentum, and temperature and humidity

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunout

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nm

  USE precision, ONLY : &
! Imported Parameters:
       dp

! jjb work done
!     - removal of unused parameters, remaining one now in modules

! jjb questions
!     - why lct.le.lcl+2, not +1, or not lct.ne.lcl ?
!     - why hardcoded values in several loops: k=10,n (why 10, why n not nm), then lct-4,lct+4

  implicit none

! External function:
  real (kind=dp), external :: p21              ! saturation water vapour pressure [Pa]

! Local parameters
  real (kind=dp), parameter :: ppFBO = 0.8_dp          ! Filter for buoyancy, fraction of old value
  real (kind=dp), parameter :: ppFBN = 1._dp - ppFBO   ! Filter for buoyancy, fraction of new value
  real (kind=dp), parameter :: ppFShO = 0.8_dp         ! Filter for Sh, fraction of old value
  real (kind=dp), parameter :: ppFShN = 1._dp - ppFShO ! Filter for Sh, fraction of new value
  real (kind=dp), parameter :: ppFSmO = 0.8_dp         ! Filter for Sh, fraction of old value
  real (kind=dp), parameter :: ppFSmN = 1._dp - ppFSmO ! Filter for Sh, fraction of new value
  real (kind=dp), parameter :: ppFxlO = 0.95_dp        ! Filter for xl, fraction of old value
  real (kind=dp), parameter :: ppFxlN = 1._dp - ppFxlO ! Filter for xl, fraction of new value
  real (kind=dp), parameter :: ppGhMin = -0.6_dp     ! Gh min (see BTZ96 eq 7, modified; MY Fig3)
  real (kind=dp), parameter :: ppGhMax = +0.03_dp    ! Gh max (see BTZ96 eq 7 and text p639)
  ! Mellor & Yamada coefficients A1, A2, B1, B2 and C1 (eq. 45 p.858) as used in eq. 34-35
  real (kind=dp), parameter :: ppMY_A1 = 0.92_dp
  real (kind=dp), parameter :: ppMY_B1 = 16.6_dp
  real (kind=dp), parameter :: ppMY_A2 = 0.74_dp
  real (kind=dp), parameter :: ppMY_B2 = 10.1_dp
  real (kind=dp), parameter :: ppMY_C1 = 0.08_dp
  ! BTZ96 coefficients a1-a9 (eq. 9) and epsilon (eq. 6)
  real (kind=dp), parameter :: ppa1 = ppMY_A2
  real (kind=dp), parameter :: ppa2 = -9_dp * ppMY_A1 * ppMY_A2**2
  real (kind=dp), parameter :: ppa3 = 18_dp * ppMY_A1**2 * ppMY_A2 * ppMY_C1
  real (kind=dp), parameter :: ppa4 = ppMY_A1 * (1_dp - 3_dp * ppMY_C1)
  real (kind=dp), parameter :: ppa5 = 3_dp * ppMY_A1 * ppMY_A2 * &
       (3_dp * ppMY_A2 + ppMY_B2 * (3_dp * ppMY_C1 - 1_dp) + 12_dp * ppMY_A1 * ppMY_C1)
  real (kind=dp), parameter :: ppa6 = -3_dp * ppMY_A2 * (7_dp * ppMY_A1 + ppMY_B2)
  real (kind=dp), parameter :: ppa7 = 27_dp * ppMY_A1 * ppMY_A2**2 * (4_dp * ppMY_A1 + ppMY_B2)
  real (kind=dp), parameter :: ppa8 = 6_dp * ppMY_A1**2
  real (kind=dp), parameter :: ppa9 = 18_dp * ppMY_A1**2 * ppMY_A2 * (3_dp * ppMY_A2 - ppMY_B2)
  real (kind=dp), parameter :: ppeps = 1_dp / ppMY_B1

! Local scalars
  integer :: k
  real (kind=dp) :: betal, betat, betaw
  real (kind=dp) :: esat, ql, qsat, qslt
  real (kind=dp) :: ghn, gmn, shn, smn
  real (kind=dp) :: x0, x1, x2, x4, xa, xb, zinv
! Local arrays
  real (kind=dp) :: es(n),xlo(n),dmw(n),dthetl(n),xmw(n)

! Common blocks:
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy

  common /cb42a/ tkeps(n),tkepb(n),tkepd(n)
  real (kind=dp) :: tkeps, tkepb, tkepd

  common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
  real (kind=dp) :: gm, gh, sm, sh, xl

  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb53a/ thet(n),theti(n)
  real(kind=dp) :: thet, theti
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a
  common /kinv_i/ kinv
  integer :: kinv

! == End of declarations =======================================================

! calculation of the buoyancy term Pb (eq. 25b)
!----------------------------------------------

! cloud free situation
  if (lct.le.lcl+2) then
     do k=2,nm
        x0 = ((1_dp + 0.61_dp * xm1(k)) * (theta(k+1) - theta(k)) &
             + 0.61_dp * theta(k) * (xm1(k+1) - xm1(k))) / deta(k)
        sm(k) = x0
        sh(k) = x0
        buoy(k) = ppFBO * buoy(k) + ppFBN * x0                  !filter
        thetl(k) = (1_dp + 0.61_dp*xm1(k)) * theta(k)           !for output only
     enddo
! calculate inversion level
     do k=10,n
        kinv=k
        if (buoy(k).gt.1.d-5) exit
     enddo
     write(jpfunout,*)'kinv : ',kinv

! cloudy situation
  else

! liquid water potential temperature
! l21/cp=2.4774d06/1005=2465.1
! l21/r0=2.4774d06/287.05=8630.6
! l21/r1=2.4774d06/461.51=5368.
! thet(k)=theta(k)/t(k); theti(k)=t(k)/theta(k)
     do k=1,nm
        thetl(k) = theta(k) - 2465.1_dp * thet(k) * xm2(k) / rho(k)
        xmw(k)   = xm1(k) + xm2(k) / rho(k)
     enddo
     thetl(n) = thetl(nm) + 1._dp
     xmw(n)   = xm1(n) + xm2(n) / rho(n)

! buoyancy term with eq. 72
     do k=2,nm
        dthetl(k) = (thetl(k+1) - thetl(k)) / deta(k)
        dmw(k) = (xmw(k+1) - xmw(k)) / deta(k)
        x0 = (1._dp + 0.61_dp * xmw(k)) * dthetl(k) + 0.61_dp * thetl(k) * dmw(k)
        sh(k) = ppFShO * sh(k) + ppFShN * x0                        !filter
     enddo
     do k=2,lct-1
        ql   = xm2(k) / rho(k)
        esat = p21(t(k))
        qsat = 0.62198_dp * esat / (p(k) - 0.37802_dp * esat)
        qslt = 5368._dp * qsat / (t(k) * t(k))
! terms a and b of eq. 70
        xa = 1._dp / (1._dp + 2465.1_dp * qslt)
        xb = xa * theti(k) * qslt
! terms betat, betal, betaw of eq.72
        betat = (1._dp + 0.61_dp * xm1(k) - ql)
        betaw = 0.61_dp * (thetl(k) + 2465.1_dp * thet(k) * ql)
        betal = (1._dp + 0.61_dp * xmw(k) - 3.22_dp * ql) * 2465.1_dp * thet(k) - 1.61 * thetl(k)
        x0    = (betat - xb * betal) * dthetl(k) + (betaw + xa * betal) * dmw(k)
        sm(k) = ppFSmO * sm(k) + ppFSmN * x0                         !filter

        x1      = 60._dp * (min(feu(k),1._dp) - 1._dp)                 !calc. of alpha (Bott 1997 eq. 6 and text)
        betal   = betal * exp(x1)
        x0      = (betat - xb * betal) * dthetl(k) + (betaw + xa * betal) * dmw(k)
        buoy(k) = ppFBO * buoy(k) + ppFBN * x0                     !filter
     enddo
     do k=lct,nm
        buoy(k) = sh(k)
     enddo
! calculate inversion level
     kinv=lct+5
     do k=lct-4,lct+4
        if (buoy(k).gt.1.e-5_dp) then
           kinv = min(kinv,k)
        endif
     enddo
     kinv = kinv-1
     write(jpfunout,*)'cloud kinv ',kinv
  endif

! -----------------------------------------------------------------------------

! calculation of the mixing length xl (eq. 50)
!---------------------------------------------

  do k=1,n
     xlo(k) = xl(k)
     es(k) = sqrt(2._dp * tke(k))
  enddo
  x0 = 0._dp
  x1 = 0._dp
  do k=2,kinv-1
     x2 = es(k) * deta(k)
     x0 = x0 + x2 * etw(k)
     x1 = x1 + x2
  enddo
  x2 = x0 / x1
  zinv = etw(kinv)
  x4 = 0.1_dp - detw(kinv) / x2
  xl(1) = 0._dp
  do k=2,kinv-1
     x0 = kappa * etw(k)
     x1 = max(detw(k), x2 * (0.1_dp - x4 * exp((etw(k) - zinv) / 15._dp)))
     xl(k) = x0 * x1 / (x0 + x1)
  enddo
  do k=kinv,n
     x0 = kappa * etw(k)
     x1 = detw(k)
     xl(k) = x0 * x1 / (x0 + x1)
  enddo
  do k=2,nm
     xl(k) = ppFxlO * xlo(k) + ppFxlN * xl(k)                    !filter (xl(1) and xl(n) do not change)
  enddo

! -----------------------------------------------------------------------------

! calculation of stability functions gh and gm (eq. 33a, b or 74a)
!-----------------------------------------------------------------
  gh(1) = 0._dp
  gm(1) = 0._dp
  do k=2,nm
     x1 = xl(k) * xl(k) / (es(k) * es(k))
     ghn = -g * x1 / theta(k) * buoy(k)
     gh(k) = max(ppGhMin, min(ghn, ppGhMax))
     gmn = x1 * ((u(k+1) - u(k))**2 + (v(k+1) - v(k))**2) / (deta(k) * deta(k)) ! BTZ96 eq. (7)
     gm(k) = min(gmn, 25._dp * (ppGhMax - gh(k)))
  enddo

! exchange coefficients
! turbulent kinetic energy production/dissipation
! first atmospheric level with frictional velocity ustern
  atkh(1) = 0.5_dp * eta(2) * ustern / gclt
  atkm(1) = 0.5_dp * eta(2) * ustern / gclu
  atke(1) = atkm(1)

! stability functions Sm, Sh (solution of eqs. 34, 35)
  do k=2,nm
     ghn = gh(k)
     gmn = gm(k)
     ! x0 is a10 in BTZ96 eq. 9
     x0 = 1._dp / (1._dp + (ppa6 + ppa7 * ghn) * ghn + (ppa8 + ppa9 * ghn) * gmn)
     shn = (ppa1 + ppa2 * ghn + ppa3 * gmn) * x0
     smn = (ppa4 + ppa5 * ghn) * x0

! TKE production: Ps = shear production
!                 Pb = buoyant production
!                 Pd = dissipation (epsilon from eq. 25c)
     x1 = es(k)**3 / xl(k)
     tkeps(k) = x1 * smn * gmn
     tkepb(k) = x1 * shn * ghn
     tkepd(k) = -x1 * ppeps
     tkep(k) = tkeps(k) + tkepb(k) + tkepd(k)

! exchange coefficients for heat (atkh)
!                           moisture (atkm)
!                           kinetic energy (atke)
     x2 = es(k) * xl(k)
     atkh(k) = x2 * shn
     atkm(k) = x2 * smn
     atke(k) = min(atkm(k), x2 * 0.2_dp)
  enddo
  do k=1,nm
     atke(k) = 0.5_dp * (atke(k) + atke(k+1))
  enddo

end subroutine atk1

!
!-------------------------------------------------------------
!

subroutine soil (dt)
! fully implicit procedure for the solution of the diffusive heat
! and moisture transport within the soil.
! for further details see subroutine difm.

  ! jjb: this routine has not been maintained, an updated version exists, please contact
  !      the authors for more information.
  !      Here, parameters (ebs, psis, ...) have been set for sandy loam only.
  !      The updated version available contains data for 12 types of soil


! Authors:
! ------
  !     (to be clarified) J. Siebert, A. Bott, R. von Glasow, other?


  USE constants, ONLY : &
! Imported Parameters:
       rhow                  ! Water density [kg/m**3]

  USE data_surface, ONLY : &
       aks, anu0, bs, bs0, ebc, ebs, psis

  USE global_params, ONLY : &
! Imported Parameters:
       nb, &
       nbm

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: dt

  real (kind=dp), parameter :: rhoc  = 1.34e6_dp  ! volumetric heat capacity for pure soil (sandy loam)  (J/m^3/K)
  real (kind=dp), parameter :: rhocw = 4.186e6_dp ! density * specific heat capacity, water [J m^-3 K^-1]
! Local scalars:
  integer :: k
  real (kind=dp) :: akb, x0, x1, x2, x3
! Local arrays:
  real (kind=dp) :: xa(nbm),xb(nbm),xc(nbm),xd(nbm),xe(nbm),xf(nbm)

! Common blocks:
  common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb), &
                ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
  real (kind=dp) :: zb, dzb, dzbw, tb, eb, ak, d, &
       ajb, ajq, ajl, ajt, ajd, ajs, ds1, ds2, ajm, reif, tau, trdep

! == End of declarations =======================================================
! soil temperature
  xe(1) = 0._dp
  x0    = max(eb(1),ebc)
  akb   = anu0 * x0**bs0 / ((1._dp - ebs) * rhoc + eb(1) * rhocw)
  xa(1) = akb * dt / (dzbw(1) * dzb(1))
  do k=2,nbm
     x0    = max(eb(k),ebc)
     akb   = anu0 * x0**bs0 / ((1._dp - ebs) * rhoc + eb(k) * rhocw)
     xa(k) = akb * dt / (dzbw(k) * dzb(k))
     xc(k) = xa(k-1) * dzbw(k-1) / dzbw(k)
     xb(k) = 1._dp + xa(k) + xc(k)
     xd(k) = xb(k) - xc(k) * xe(k-1)
     xe(k) = xa(k) / xd(k)
  end do

  xf(1)=tb(1)
  do k=2,nbm
     xf(k)=(tb(k)+xc(k)*xf(k-1))/xd(k)
  end do
  do k=nbm,2,-1
     tb(k)=xe(k)*tb(k+1)+xf(k)
  end do

! volumetric moisture content
  x0 = 2._dp * bs + 3._dp
  x1 = bs + 2._dp
  x2 = -bs * aks * psis / ebs
  do k=2,nbm
     x3    = (eb(k)+dzbw(k)*(eb(k+1)-eb(k)) / (2._dp*dzb(k))) / ebs
     ak(k) = aks * x3**x0
     d(k)  = x2 * x3**x1
  end do
  ak(1) = 0._dp
  d(1)  = 0._dp
  if (abs(eb(2)-eb(1)).gt.1.e-5_dp) then
     d(1) = ajm*dzb(1) / (rhow*(eb(2)-eb(1)))
  end if
  xa(1) = d(1)*dt / (dzbw(1)*dzb(1))
  do k=2,nbm
     xa(k) = d(k) * dt / (dzbw(k) * dzb(k))
     xc(k) = xa(k-1) * dzbw(k-1) / dzbw(k)
     xb(k) = 1._dp + xa(k) + xc(k)
     xd(k) = xb(k) - xc(k) * xe(k-1)
     xe(k) = xa(k) / xd(k)
  end do
  xf(1) = eb(1)
  do k=2,nbm
     xf(k) = (eb(k) + dt/dzbw(k) * (ak(k-1) - ak(k)) + xc(k)*xf(k-1)) / xd(k)
  end do
  do k=nbm,2,-1
     eb(k) = xe(k) * eb(k+1) + xf(k)
  end do

end subroutine soil

!
!-------------------------------------------------------------
!

subroutine surf0 (dt)
!
! Description:
! -----------
  ! lower boundary condition for water surface
  ! constant temperature and saturation specific humidity


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  ! 06-May-2021  Josue Bock  Introduce ltwcst, ntwopt and rhsurf from namelist instead of hard-coded

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
       ltwcst,       &
       ntwopt,       &
       rhsurf,       &
! Imported Routines:
       abortM

  USE constants, ONLY : &
! Imported Parameters:
       g

  USE data_surface, ONLY : &
       tw, &                       ! water surface temperature
       ustern, z0,&                ! frictional velocity, roughness length
       gclu, gclt                  ! coefficients for momentum, and temperature and humidity

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: dt

! External function:
  real (kind=dp), external :: p21              ! saturation water vapour pressure [Pa]

! Local scalars:
  real (kind=dp) :: uu, vv, vqr, vbt                      ! wind velocities
  real (kind=dp) :: xnvl, zp, zpdz0, zpdl, cu, ctq        ! Clarke functions
  real (kind=dp) :: zp21 ! p21(tw)

! Common blocks:
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw
  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a

! == End of declarations =======================================================

! tw varying with time
  if (.not.ltwcst) then
     select case (ntwopt)
     ! different options used so far. Add other functions of time if needed
     case (1)
        tw = tw - 5.787e-6_dp * dt
     case (2)
        tw = tw - 6.94444e-6_dp * dt
     case default
        call abortM ('Error in SR surf0: wrong choice for ntwopt, change in namelist')
     end select
  endif

! temperature and specific humidity at the surface
  t(1)   = tw
  zp21   = p21(tw)
  xm1(1) = rhsurf * 0.62198_dp * zp21 / (p(1) - 0.37802_dp * zp21)

! frictional velocity
  uu  = u(2)
  vv  = v(2)
  vqr = uu*uu+vv*vv
  vbt = sqrt(vqr)

  zp    = 0.5_dp * eta(2) + z0
  zpdz0 = log(zp/z0)
  xnvl  = g * (theta(2) - tw) * 2._dp / (theta(2) + tw)
  zpdl  = zp * xnvl / vqr

  call claf (zpdl,zpdz0,cu,ctq)

  ustern = max(0.01_dp, vbt/cu)
! roughness length:
! charnock's relation: z0=0.015*ustern**2/g
  z0   = 0.015_dp * ustern * ustern / g
  gclu = cu
  gclt = ctq

end subroutine surf0

!
!-------------------------------------------------------------
!

subroutine surf1 (dt)

! calculation of surface temperature and volumetric moisture content
! by means of balance of fluxes at the surface following Pielke,
! mesoscale meteorological modelling, chapter 11.

  USE constants, ONLY : &
! Imported Parameters:
       cp, &                   ! Specific heat of dry air, in J/(kg.K)
       g,  &                   ! Gravitational acceleration (m/s**2)
       r1, &                   ! Specific gas constant of water vapour, in J/(kg.K)
       rhow                  ! Water density [kg/m**3]

  USE data_surface, ONLY : &
       aks, anu0, bs, bs0, ebc, ebs, psis, &
       ustern, z0,&                ! frictional velocity, roughness length
       gclu, gclt                  ! coefficients for momentum, and temperature and humidity

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunprofm

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nb

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  real (kind=dp), intent(in) :: dt

  real (kind=dp), external :: p21  ! saturation vapour pressure
  real (kind=dp), external :: p31  ! saturation vapour pressure (ice)
  real (kind=dp), external :: xl21 ! latent heat of vaporisation = f(temperature)

! Local parameters
  real (kind=dp), parameter :: al31 = 2.835e+6_dp
  real (kind=dp), parameter :: sigma = 5.6697e-8_dp  ! Stefan-Boltzmann-constant
  real (kind=dp), parameter :: t0 = 273.15_dp

  ! Local scalars:
  integer :: i, iaa
  logical :: l1
  real (kind=dp) :: anu, bs3
  real (kind=dp) :: ctq, cu ! output of claf function
  real (kind=dp) :: ddew, det, djbde, djbdt, djmde, djqde, djqdt, djtdt
  real (kind=dp) :: eb1, eb10
  real (kind=dp) :: f1e, f1t, f2e, f2t, fqs, fts
  real (kind=dp) :: ps, psi1, psi2, qq2, qs, qst, rak1, rrho
  real (kind=dp) :: ts, ts0, tst, uwr
  real (kind=dp) :: uu, vv, vqr, vbt ! surf variables for claf function
  real (kind=dp) :: x0, x1, xa, xm21s, xnvl, zp, zpdl, zpdz0

! Local arrays:
  real (kind=dp) :: ebb(40),tss(40),ftss(40),fqss(40)

! Common blocks:
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw
  common /cb45/ u(n),v(n),w(n)
  real (kind=dp) :: u, v, w
  common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb), &
                ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
  real (kind=dp) :: zb, dzb, dzbw, tb, eb, ak, d, &
       ajb, ajq, ajl, ajt, ajd, ajs, ds1, ds2, ajm, reif, tau, trdep
  common /cb48/ sk,sl,dtrad(n),dtcon(n)
  real (kind=dp) :: sk, sl, dtrad, dtcon
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a

! statement function
  real (kind=dp) :: cm, pp
  cm(pp)=0.62198_dp*pp/(ps-0.37802_dp*pp)

! == End of declarations =======================================================
  rrho=rho(1)
  uu=u(2)
  vv=v(2)
  vqr=max(uu*uu+vv*vv,1.d-12)
  vbt=sqrt(vqr)
!  gamma=g/cp
  bs3=2._dp*bs+3._dp
!  bs2=bs+2._dp
  psi2=psis*(ebs/eb(2))**bs
  qq2=xm1(2)
  ts=t(1)
  ps=p(1)
  eb1=eb(1)
  l1=.false.

  zp    = deta(1) + z0
  zpdz0 = log(zp/z0)
  xnvl  = g * (theta(2) - ts) * 2._dp / (theta(2) + ts)
  zpdl  = zp * xnvl / vqr
  call claf (zpdl,zpdz0,cu,ctq)
  ustern=max(0.01_dp,vbt/cu)

! specific humidity at the surface
  if (ts.ge.t0) then
     xm21s = cm(p21(ts))
  else
     xm21s = cm(p31(ts))
  endif
  psi1=psis*(ebs/eb1)**bs
  qs=xm21s*exp(g*psi1/(r1*ts))
! tstern,qstern
  tst=(theta(2)-ts*(1._dp+0.608_dp*qs))/ctq
  qst=(qq2-qs)/ctq
! heatflux from ground
  x1=max(eb1,ebc)
  anu=anu0*x1**bs0
  ajb=anu*(tb(2)-ts)/dzb(1)
! microturbulent flux of water vapor
  ajq=rrho*ustern*qst
! latent microturbulent enthalpy flux
! ajs is water flux due to droplet sedimentation;
  if (ts.lt.t0) then
     ajl=al31*ajq-(al31-xl21(ts))*ajs
  else
     ajl=xl21(ts)*ajq
  endif
! sensible microturbulent enthalpy flux
  ajt=rrho*cp*ustern*tst
! ground moisture flux
!  rak1=rhow*aks*(eb1/ebs)**bs3
  xa=0.5_dp
  rak1=rhow*aks*((xa*eb1+(1._dp-xa)*eb(2))/ebs)**bs3
  ajm=rak1*((psi2-psi1)/dzb(1)-1._dp)
!  debs=-bs*aks*psis/ebs
!  ajm=(debs*((eb1+eb(2))/(2._dp*ebs))**bs2*(eb(2)-eb1)/dzb(1)- &
!      aks*((eb1+eb(2))/(2._dp*ebs))**bs3)*rhow
! dew flux
  ajd=0._dp
  x0=ajq+ajm+ajs
  if (eb1.ge.ebs) then
     ajd=-x0
     ddew=tau/dt
     if (x0.lt.0._dp) ajd=min(-x0,ddew)
     ddew=ddew-ajd
  end if
! flux balances
  fts=sl+sk+ajb+ajl+ajt-sigma*ts**4
  fqs=x0+ajd
  if (l1) goto 2010

! 2-dimensional newton-raphson iteration for ts and eb(1)
  do iaa=1,20
     ebb(iaa)=eb1
     ftss(iaa)=fts
     fqss(iaa)=fqs
     tss(iaa)=ts
!     fts0=fts
!     fqs0=fqs
     ts0=ts
     eb10=eb1
! calculation of derivatives d(fluxes)/d(eta) and d(fluxes)/d(T)
! djbde=d(ajb)/d(eta); djbdt=d(ajb)/d(T) etc.
     djbde=0._dp
     if (eb1.gt.ebc) djbde=ajb*bs0/eb1
     djbdt=-anu/dzb(1)
     djqde=rrho*ustern*qs*g*bs*psi1/(ctq*r1*ts*eb1)
     x0=p21(ts)
     djqdt=rrho*ustern*qs/ctq*(g*psi1/(r1*ts*ts)+ &
           x0*4027.163_dp/((x0-.37802_dp*ps)*(ts-38.33_dp)**2))
     djtdt=-rrho*cp*ustern/ctq
!     djmde=(ajm*bs3+rak1/dzb(1)*psi1*bs)/eb1
     djmde=rak1/dzb(1)*psi1*bs/eb1
!     djmde=(debs*(eb1/ebs)**bs2/dzb(1)*((eb(2)-eb1)*bs2/eb1-1._dp)-
!    & bs3*aks/eb1*(eb1/ebs)**bs3)*rhow
!     djmde=min(djmde,0._dp)
! coefficients for solution of linear equation system
     x0=xl21(ts)
     f1e=djbde+x0*djqde
     f1t=djbdt-2335.5_dp*ajq+x0*djqdt+djtdt-4._dp*sigma*ts*ts*ts
     f2e=djqde+djmde
     f2t=djqdt
     det=f1e*f2t-f1t*f2e
!     if (abs(det).lt.1.e-10_dp) x0=sqrt(1.d0-2.d0)
     if (abs(det).lt.1.e-10_dp) stop 'SR surf1'
! new values of ts and eb1
     ts=ts+(fts*f2e-fqs*f1e)/det
     eb1=eb1+(fqs*f1t-fts*f2t)/det
     eb1=min(eb1,ebs)
     eb1=max(eb1,ebs/15._dp)
     if (ddew.gt.0._dp) eb1=ebs
     if (ts.gt.300._dp .or. ts.lt.250._dp) ts = ts0 - 0.01_dp
! new surface fluxes
     if (ts.ge.t0) then
        xm21s = cm(p21(ts))
     else
        xm21s = cm(p31(ts))
     endif
     psi1=psis*(ebs/eb1)**bs
     qs=xm21s*exp(g*psi1/(r1*ts))
     qst=(qq2-qs)/ctq
     tst=(theta(2)-ts*(1._dp + 0.608_dp * qs))/ctq
     x1=max(eb(1),ebc)
     anu=anu0*x1**bs0
     ajb=anu*(tb(2)-ts)/dzb(1)
     ajq=rrho*ustern*qst
     if (ts.lt.t0) then
        ajl=al31*ajq-(al31-xl21(ts))*ajs
     else
        ajl=xl21(ts)*ajq
     endif
     ajt=rrho*cp*ustern*tst
     ajm=rak1*((psi2-psi1)/dzb(1)-1._dp)
     ajd=0._dp
     x0=ajq+ajm+ajs
     if (eb1.ge.ebs) then
        ajd=-x0
        ddew=tau/dt
        if (x0.lt.0._dp) ajd=min(-x0,ddew)
        ddew=ddew-ajd
     end if
! flux balances
     fts=sl+sk+ajb+ajl+ajt-sigma*ts**4
     fqs=x0+ajd
! convergence criteria
     if(dabs(ts-ts0).le.1.d-2.and.dabs(eb1-eb10).le.1.d-3)go to 2030
     if (dabs(fts).le.1.d-1.and.dabs(fqs).le..1*dabs(ajq)) goto 2030
  end do

  write (jpfunprofm,6000) eb1,ts,fts,fqs
 6000 format (10x,'no convergence of ts- and eb1-iteration:'/ &
     & 'eb1',f16.4,'ts',f16.4,'fts',f16.4,'fqs',f16.4)
  write (jpfunprofm,6010) (ebb(i),tss(i),ftss(i),fqss(i),i=1,20)
 6010 format (10x,3f16.4,e16.4)

2030 if (tau.gt.0..and.ts.lt.t0) l1=.true.
  if (ts.gt.t0.and.reif.gt.0) l1=.true.
  if (.not.l1) go to 2010
  ts=t0
! change of dew or rime by evaporation/condensation
 2010 continue
  if (ts.ge.t0) tau=tau-ajd*dt
  if (ts.lt.t0) reif=reif-ajd*dt
  if (.not.l1) goto 2040
! coexistence of dew and rime
  uwr=min(dt*fts/3.35d+5,reif)
  uwr=max(uwr,-tau)
  tau=tau+uwr
  reif=reif-uwr
2040 tau=max(0._dp,tau)
  reif=max(0._dp,reif)
! storage of surface values
  t(1)=ts
  tb(1)=ts
  xm1(1)=qs
  feu(1)=qs/xm21s
  eb(1)=eb1

  xnvl = g * (theta(2) - ts) * 2._dp / (theta(2) + ts)
  zpdl = zp * xnvl / vqr
  call claf (zpdl,zpdz0,cu,ctq)
  ustern = max(0.01_dp,vbt/cu)
  gclu   = cu
  gclt   = ctq

end subroutine surf1

!
!-------------------------------------------------------------
!

function p31(t)
! water vapour pressure over ice

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  real (kind=dp), intent(in) :: t
  real (kind=dp) :: p31
  real (kind=dp), parameter :: t1=273.16_dp
  real (kind=dp) :: xlog10

  xlog10 = -9.09685_dp * (t1/t - 1._dp) &
           -3.56654_dp * log10(t1/t) &
           +0.87682_dp * (1._dp - t/t1) + 0.78614_dp

  p31=100._dp*(10._dp**xlog10)

end function p31

!
!-------------------------------------------------------------
!

subroutine claf (zpdl,zpdz0,u,tq)
!
! Description:
!    interpolation of Clarke functions by means of tabulated values
!    u: clarke function for momentum; tq: for temperature, humidity etc.


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  !

!
! History:
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      01/2017   Removed labels and goto's                     <Josue Bock>
!                    Use Fortran generic functions min and max
! 1.1      08/2016   Use Andreas Bott's updated version from str   <Josue Bock>
!                    Header
!                    Declaration of all variables and imput none
!
! 1.0       ?        Original code used in Mistra v741             <Andreas Bott>
!
!
! == End of header =============================================================

! Declarations:
! ------------
! Modules used:

  USE data_surface, ONLY : &
       fu, ft, xzpdl, xzpdz0 ! Clarke tabled data

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: zpdl, zpdz0
! Scalar arguments with intent(out):
  real (kind=dp), intent(out) :: u, tq

! Local scalars:
  integer :: i,nl,nz
  real (kind=dp) :: dx, dy
  real (kind=dp) :: zpdla, zpdz0a

! == End of declarations =======================================================

! xzpdl tabled values range from -5.5 to 3.0. Here, zpdla is forced to be within this range
  zpdla = max(zpdl,-5.5_dp)
  zpdla = min(zpdla,3._dp)
! xzpdz0 tabled values range from 3. to 17. zpdz0a is forced to be lower equal to the max value.
  zpdz0a = min(zpdz0,17._dp)

! Find nl and nz: indexes of the tabled values greater or equal to the actual values
  do i=2,18                     ! note that zpdla > xzpdl(1), thus nl >= 2
     nl=i                       !  thus no out-of-bounds problem with index nl-1 in dx (see below)
     if (zpdla.lt.xzpdl(i)) exit
  enddo

  do i=1,7
     nz=i
     if (zpdz0a.lt.xzpdz0(i)) exit
  enddo

! dx: proportionality factor
  dx = (zpdla - xzpdl(nl-1)) / (xzpdl(nl) - xzpdl(nl-1))

  if (nz.eq.1) then
     dy = zpdz0a / xzpdz0(1)
     u = (fu(nl,1) * dx + fu(nl-1,1) * (1._dp - dx)) * dy
     if (zpdl.ge.0._dp) then
        tq = u / 1.35_dp
     else
        tq = (ft(nl,1) * dx + ft(nl-1,1) * (1._dp - dx)) * dy / 1.35_dp
     endif
  else
     dy = (zpdz0a - xzpdz0(nz-1)) / (xzpdz0(nz) - xzpdz0(nz-1))
     u = fu(nl-1,nz-1) &
          + (fu(nl,nz-1) - fu(nl-1,nz-1)) * dx &
          + (fu(nl-1,nz) - fu(nl-1,nz-1)) * dy &
          + (fu(nl,nz) - fu(nl-1,nz) + fu(nl-1,nz-1) - fu(nl,nz-1)) * dx * dy
     if (zpdl.ge.0._dp) then
        tq = u / 1.35_dp
     else
        tq = (ft(nl-1,nz-1) &
             + (ft(nl,nz-1) - ft(nl-1,nz-1)) * dx &
             + (ft(nl-1,nz) - ft(nl-1,nz-1)) * dy &
             + (ft(nl,nz) - ft(nl-1,nz) + ft(nl-1,nz-1) - ft(nl,nz-1)) * dx * dy) / 1.35_dp
     endif
  end if

end subroutine claf

!
!-------------------------------------------------------------
!

subroutine kon (dt,chem)
!
! Description:
! -----------
  ! Driving routine for the calculation  of the diffusional droplet growth
  ! by condensation. Formulation for explicit cloud microphysics with
  ! two-dimensional droplet and aerosol distribution


! Author:
! ------
  !    Andreas Bott
  !    Roland von Glasow (chemistry part)


! Modifications :
! -------------
  !     ?        Josue Bock  corrected initialisation, nkc_l replaced by nkc
  !

  !     jjb work done
  !      removed one unused parameter; use module instead
  !      removed k0=k, used only in equil argument. For the sake of clarity

  ! 01-May-2021  Josue Bock  remove potot, computed but unused
  !                          corrected the array size in blck07 and blck08: nf+1 instead of n
  !                            (note that indexes 2:nf only are used in SR konc)
  !                          potential bugfix: compute "chemical" arrays whatever rH (<=70% or >70%)

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

!      USE config, ONLY :
!     &     nkc_l

  USE constants, ONLY : &
! Imported Parameters:
       pi

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nka, &
       nkt, &
       nkc, &
       mb

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Interface (required for SR equil, which has an optional argument)
  interface
     subroutine equil (ncase,kk)
       integer,           intent(in)  :: ncase
       integer, optional, intent(in)  :: kk    ! model level for calculation
     end subroutine equil
  end interface

! Subroutine arguments
  logical,       intent(in) :: chem
  real (kind=dp),intent(in) :: dt

! Local parameters:
  ! optimisation: define parameters that will be computed only once
  real (kind=dp), parameter :: z4pi3 = 4._dp * pi / 3._dp

! Local scalars:
  integer :: ia, ib, jt, k, kr
  real (kind=dp), external :: p21
  real (kind=dp) :: dfdt, feualt, pp, to, tn, xm1o, xm1n ! arguments of SR subkon
  real (kind=dp) :: ffk(nkt,nka), totr(mb)               ! arguments of SR subkon
  !real (kind=dp) :: potot(nkc,n)
  ! JJB temporary, check if my bugfix is justified or not
  logical :: lfeu, lcheck
  ! end JJB temproray

! Common blocks:
  common /cb11/ totrad (mb,n)
  real (kind=dp) :: totrad

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb48/ sk,sl,dtrad(n),dtcon(n)
  real (kind=dp) :: sk, sl, dtrad, dtcon

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a

  common /blck06/ kw(nka),ka
  integer :: kw, ka
  common /blck07/ part_o_a(nka,nf+1),part_o_d(nka,nf+1), &              ! old aerosol and droplet
                  part_n_a(nka,nf+1),part_n_d(nka,nf+1),pntot(nkc,nf+1) ! new aerosol and droplet, total per chem bin
  real (kind=dp) :: part_o_a, part_o_d, part_n_a, part_n_d, pntot
  common /blck08/ vol1_a(nka,nf+1),vol1_d(nka,nf+1),vol2(nkc,nf+1)
  real (kind=dp) :: vol1_a, vol1_d, vol2

! == End of declarations =======================================================

! Initialise arrays for chemistry
  if (chem) then
     ! (nka,nf+1) arrays
     vol1_a(:,:)=0._dp
     vol1_d(:,:)=0._dp
     part_o_a(:,:)=0._dp
     part_o_d(:,:)=0._dp
     part_n_a(:,:)=0._dp
     part_n_d(:,:)=0._dp
     ! (nkc,nf+1) arrays
     vol2(:,:)=0._dp
     !potot(:,:)=0._dp
     pntot(:,:)=0._dp
     ! JJB temproray
     lfeu = .false.
     lcheck = .false.
     ! end JJB temproray
  end if

  kloop: do k=2,nf+1
     dtcon(k) = 0._dp
     ! set local variables
     tn   = t(k)
     xm1n = xm1(k)
     pp   = p(k)
     ffk(:,:) = ff(:,:,k) ! jjb: moved from below (rH>70% case) for chem calculations

     ! Compute particle/volume for chemistry before
     !   (this was done only in the rH>0.7 case, now done whatever the case)
     if (chem) then

! aerosol + cloud chemistry:
! vol2 and vol1 old liquid volume in class (1-nkc) and row/class (1-nkc,1-nka)
! (um^3/cm^3), used in SR konc to shift moles from aerosol to drop
! and vice versa.
! part_o and part_n old and new part. conc. in row/class (1-nkc,1-nka) (cm^-3)
! take care if non-soluble parts are in aerosol

        ! small (dry) aerosol: aerosol (1) and droplet (3)
        do ia=1,ka
           ! aerosol: jt=1,kw(ia)
           vol1_a(ia,k)   = sum( ffk(1:kw(ia),ia) * z4pi3 * rq(1:kw(ia),ia)**3 )
           part_o_a(ia,k) = sum( ffk(1:kw(ia),ia) )
           ! droplets: jt=kw(ia)+1,nkt
           vol1_d(ia,k)   = sum( ffk(kw(ia)+1:nkt,ia) * z4pi3 * rq(kw(ia)+1:nkt,ia)**3 )
           part_o_d(ia,k) = sum( ffk(kw(ia)+1:nkt,ia) )
        end do
        vol2(1,k) = sum(vol1_a(1:ka,k))
        vol2(3,k) = sum(vol1_d(1:ka,k))
        !potot(1,k) = sum(part_o_a(1:ka,k))
        !potot(3,k) = sum(part_o_d(1:ka,k))

        ! large (dry) aerosol: aerosol (2) and droplet (4)
        do ia=ka+1,nka
           vol1_a(ia,k)   = sum( ffk(1:kw(ia),ia) * z4pi3 * rq(1:kw(ia),ia)**3)
           part_o_a(ia,k) = sum( ffk(1:kw(ia),ia) )
           vol1_d(ia,k)   = sum( ffk(kw(ia)+1:nkt,ia) * z4pi3 * rq(kw(ia)+1:nkt,ia)**3)
           part_o_d(ia,k) = sum( ffk(kw(ia)+1:nkt,ia) )
        end do
        vol2(2,k) = sum(vol1_a(ka+1:nka,k))
        vol2(4,k) = sum(vol1_d(ka+1:nka,k))
        !potot(2,k) = sum(part_o_a(ka+1:nka,k))
        !potot(4,k) = sum(part_o_d(ka+1:nka,k))
     end if

! dry case: update of humidified aerosol with Koehler curve
!----------------------------------------------------------
     if (feu(k).lt.0.7_dp) then
        feu(k)=xm1n*pp/((0.62198_dp + 0.37802_dp * xm1n)*p21(tn))
        call equil (1,k)
        if (chem) then
           ! equil update directly ff array; for chemistry calculations below, update ffk
           ffk(:,:) = ff(:,:,k)
           !! JJB temproray
           !lfeu = .true.
           !! end JJB temproray
        end if


! moist case, calculate condensational droplet growth
! ---------------------------------------------------
     else
        !! JJB temproray
        !lfeu = .false.
        !! end JJB temproray

! set input values for condensation calculation
        dfdt   = dfddt(k)
        to     = talt(k)
        xm1o   = xm1a(k)
        feualt = feu(k)
        ! radiation effect
        do ib=1,mb
           totr(ib)=totrad(ib,k)
        enddo
        kr=nar(k)
        ! jjb: moved above for chem calculations
        !do ia=1,nka
        !    do jt=1,nkt
        !       ffk(jt,ia)=ff(jt,ia,k)
        !    enddo
        ! enddo

! condensation and evaporation of droplets
! to: temperature after condensation of last timestep;
! tn: temperature before condensation of current timestep;
! tn = to + diffusion, radiation etc.
! same with xm1o and xm1n

        call subkon (dt, ffk, totr, dfdt, feualt, pp, to, tn, xm1o, xm1n, kr)

! new values
        t(k)     = to
        talt(k)  = to
        xm1(k)   = xm1o
        xm1a(k)  = xm1o
        feu(k)   = xm1o*pp/((0.62198_dp + 0.37802_dp * xm1o)*p21(to))
        dfddt(k) = (feu(k)-feualt)/dt
        xm2(k) = 0._dp
        do ia=1,nka
           do jt=1,nkt
              xm2(k) = xm2(k)+ffk(jt,ia)*e(jt)
              ff(jt,ia,k) = ffk(jt,ia)
           enddo
        enddo
        dtcon(k)=(to-tn)/dt
     end if

     if (chem) then
        ! small (dry) aerosol: aerosol (1) and droplet (3)
        do ia=1,ka
           part_n_a(ia,k)=sum(ffk(1:kw(ia),ia))
           part_n_d(ia,k)=sum(ffk(kw(ia)+1:nkt,ia))
           !! JJB temporary check
           !if (lfeu) then
           !   if (abs(part_n_a(ia,k)-part_o_a(ia,k)) > 1d-12*part_o_a(ia,k)) then
           !      print*,'JJB SR kon: bugfix justified, change of particle(a1) in equil case (rH<=70%)',ia
           !      print*,part_n_a(ia,k),part_o_a(ia,k)
           !      print*,k,feu(k)
           !      lcheck=.true.
           !   end if
           !   if (abs(part_n_d(ia,k)-part_o_d(ia,k)) > 1d-12*part_o_d(ia,k)) then
           !      print*,'JJB SR kon: bugfix justified, change of particle(d3) in equil case (rH<=70%)',ia
           !      print*,part_n_d(ia,k),part_o_d(ia,k)
           !      print*,k,feu(k)
           !      lcheck=.true.
           !   end if
           !end if
           !! end JJB
        end do
        pntot(1,k)=sum(part_n_a(1:ka,k))
        pntot(3,k)=sum(part_n_d(1:ka,k))

        ! large (dry) aerosol: aerosol (2) and droplet (4)
        do ia=ka+1,nka
           part_n_a(ia,k)=sum(ffk(1:kw(ia),ia))
           part_n_d(ia,k)=sum(ffk(kw(ia)+1:nkt,ia))
           !! JJB temporary check
           !if (lfeu) then
           !   if (abs(part_n_a(ia,k)-part_o_a(ia,k)) > 1d-12*part_o_a(ia,k)) then
           !      print*,'JJB SR kon: bugfix justified, change of particle(a2) in equil case (rH<=70%)',ia
           !      print*,part_n_a(ia,k),part_o_a(ia,k)
           !      lcheck=.true.
           !   end if
           !   if (abs(part_n_d(ia,k)-part_o_d(ia,k)) > 1d-12*part_o_d(ia,k)) then
           !      print*,'JJB SR kon: bugfix justified, change of particle(d4) in equil case (rH<=70%)',ia
           !      print*,part_n_d(ia,k),part_o_d(ia,k)
           !      lcheck=.true.
           !   end if
           !end if
           !! end JJB
        enddo
        pntot(2,k)=sum(part_n_a(ka+1:nka,k))
        pntot(4,k)=sum(part_n_d(ka+1:nka,k))
     endif

  end do kloop

! cloudy region: lowest (lcl) and highest (lct) cloud layer
  do k=nf+1,1,-1
     lct=k
     if (xm2(k).gt.1.e-5_dp) exit
  enddo
  do k=1,lct
     lcl=k
     if (xm2(k).gt.1.e-5_dp) exit
  enddo

  !! JJB temproray
  !if (chem.and.lcheck) then
  !   if (lcheck) stop 'special case SR kon: bugfix was justified, please remove stop in SR kon and proceed'
  !end if
  !! end JJB temproray

! update chemical species
  if (chem) call konc



end subroutine kon

!
!-------------------------------------------------------------
!

subroutine equil (ncase,kk)
!
! Description:
! -----------
  ! Equilibrium values for radius rg(i) and water mass eg(i)
  ! of humidified aerosol particles at given relative humidity
  ! New distribution of the particles on their equilibrium positions


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  !     ?        Roland      Add the ncase switch, expecially for the box case
  !              von Glasow?
  !
  ! 27-Apr-2021  Josue Bock  Review and merge with latest version provided by A. Bott
  !                          Add the case construct and consistency checks (arguments)
  !                          The "optional" attribute require interface when this routine is called
  !
  ! 28-Apr-2021  Josue Bock  Create case 0 for initialisation, and remove SR vgleich which was mostly a duplicate

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! Imported Routines:
       abortM

  USE constants, ONLY : &
! Imported Parameters:
       pi,              &
       rho3,            &       ! aerosol density
       rhow                     ! water density

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunerr

  USE global_params, ONLY : &
! Imported Parameters:
       nf,                  &   ! model levels
       n,                   &
       nka,                 &   ! droplet and aerosol classes
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
  integer,           intent(in)  :: ncase
  integer, optional, intent(in)  :: kk    ! model level for calculation

! Local parameters:
  ! optimisation: define parameters that will be computed only once
  real (kind=dp), parameter :: zrho_frac = rho3 / rhow
  real (kind=dp), parameter :: z4pi3 = 4.e-09_dp * pi / 3._dp
  real (kind=dp), parameter :: zfeumax = 0.99999_dp

! Local scalars:
  integer :: kmin, kmax
  integer :: k, ia, jt      ! running indices
  real (kind=dp) :: a0, b0
  real (kind=dp), external :: rgl

! Local arrays:
  real (kind=dp) :: rg(nka),eg(nka)

! Common blocks:
  common /cb44/ a0m,b0m(nka)         ! a0m, b0m: Koehler curve
  real (kind=dp) :: a0m,b0m

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), & ! e, ew, rn: aerosol / water grid
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)                    ! ff: particles
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)      ! t: temperature
  real(kind=dp) :: theta, thetl, t, talt, p, rho                ! rho: density
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)           ! feu: humidity
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a                   ! xm2: liquid water content

! == End of declarations =======================================================

! get equilibrium distribution for
!     case 0: all layers 2-n during initialisation
!     case 1: active layer (1D: k<=nf (called from SR kon), 0D: active layer)
!     case 2: all layers above nf

  select case (ncase)
  case (0)
     if (present(kk)) then
        write(jpfunerr,*)"Warning, inconsistency when calling SR equil:"
        write(jpfunerr,*)"  case 0 means equil is computed for levels 2 to n,"
        write(jpfunerr,*)"  the optional argument #2 is thus unused"
     end if
     kmin=2
     kmax=n
     feu(kmin:kmax)=min(feu(kmin:kmax),zfeumax)

  case (1)
     if (.not.present(kk)) then
        write(jpfunerr,*)"Error, inconsistency when calling SR equil:"
        write(jpfunerr,*)"  case 1 require optional argument #2 (level where the"
        write(jpfunerr,*)"  equilibrium has to be computed)"
        call abortM('Error in SR equil')
     end if
     kmin=kk
     kmax=kk

  case (2)
     if (present(kk)) then
        write(jpfunerr,*)"Warning, inconsistency when calling SR equil:"
        write(jpfunerr,*)"  case 2 means equil is computed for levels nf+1 to n,"
        write(jpfunerr,*)"  the optional argument #2 is thus unused"
     end if
     kmin=nf+1
     kmax=n

  case default
     write(jpfunerr,*)"Error, SR equil, ncase must be 1 or 2"
     call abortM('Error in SR equil')
  end select


  kloop: do k=kmin,kmax

! collect aerosol in lowest water class
!    (not needed during initialisation, all particles are already in the first water class (jt=1))
     if (ncase.gt.0) then
        do ia=1,nka
           ff(1,ia,k) = sum(ff(:,ia,k))
           ff(2:nkt,ia,k) = 0._dp
        enddo
     end if

! calculate equilibrium radius with Newton iteration (rgl)
     a0 = a0m / t(k)
     do ia=1,nka
        b0 = b0m(ia) * zrho_frac
! b0=b0m*rho3/rhow; rho3=2000; rhow=1000
        rg(ia) = rgl(rn(ia), a0, b0, feu(k))
        eg(ia) = z4pi3 * (rg(ia)**3 - rn(ia)**3)
     enddo

! new particle distribution
     do ia=1,nka
        jt=1
        do while (eg(ia).gt.ew(jt))
           jt = jt+1
        end do
        if (jt.ne.1) then
           ff(jt,ia,k) = ff(1,ia,k)
           ff(1,ia,k)  = 0._dp
        endif
     enddo

! update of total liquid water content
     xm2(k) = 0._dp
     do ia=1,nka
        do jt=1,nkt
           xm2(k) = xm2(k) + ff(jt,ia,k) * e(jt)
        enddo
     enddo

  enddo kloop

end subroutine equil

!
!-------------------------------------------------------------
!

subroutine subkon (dt, ffk, totr, dfdt, feualt, pp, to, tn, xm1o, xm1n, kr)
!
! Description:
! -----------
  !  Calculation of the diffusional droplet growth by condensation.
  !  Formulation for explicit cloud microphysics with two-dimensional droplet
  !  and aerosol distribution.
  !  All formulas and constants after Pruppacher and Klett Chapter 13.
  !  Droplet growth equation after Davies, J. Atmos. Sci., 1987


! Author:
! ------
  !    Andreas Bott


!
! Variables:
! -----------
  ! xl21 latent heat of evaporation; p21t water vapor pressure at satur.
  ! xl mean free path;
  ! deltav water vapor jump; deltat temperature jump
  ! xdvs,xkas diffusivity, conductivity corrected for small droplets
  ! xdv0, xka0 factors used for calculating xdvs and xkas.
  ! all terms in MKSA units only droplet mass and radius in mg and microns
  ! a0[microns]=2*sigma/(r1*rhow*t)*1.e6; sigma=76.1e-3=surface tension
  ! b0 solution term of koehler equation: b0=xnue*xmol2/xmol3*xa*ma/mw
  ! xnue: number of ions; xmol2 (xol3) mol masses of pure water (aerosols)
  ! xa volume fraction of soluble to unsoluble part of aerosol
  ! mw (ma) masses of water (aerosol nucleus)
  ! p21(tr)=p21(t)*exp(a0/r-b0*ma/mw)


! Modifications :
! -------------
  ! jjb work done:
  !     - reindexed all (nka,nkt) arrays for computing efficiency
  !     - defined water diffusivity as an external function
  !     - bugfix: residue (res) was not initialised. Ok if compiler sets all variables=0

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

  USE constants, ONLY : &
! Imported Parameters:
       cp, &                   ! Specific heat of dry air, in J/(kg.K)
       pi, &
       r0, &                   ! Specific gas constant of dry air, in J/(kg.K)
       r1, &                   ! Specific gas constant of water vapour, in J/(kg.K)
       rhow                    ! Water density [kg/m**3]

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunout

  USE global_params, ONLY : &
! Imported Parameters:
       jptaerrad, &
       nka, &
       nkt, &
       mb

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
  integer,       intent(in)    :: kr
  real (kind=dp),intent(in)    :: dt, totr(mb), dfdt, feualt, pp, tn, xm1n ! tn, xm1n: before current condensation calculation
  real (kind=dp),intent(inout) :: ffk(nkt,nka), to, xm1o                   ! to, xm1o: after current condensation calculation

! Local scalars
  integer :: ia, ib0, ib, itk, jt, jtp, kr0

  real (kind=dp) :: xldcp, xka, xdv, xl, &
       xdvs, xkas, &
       deltav, deltat, &
       rho, rho21, rho21s,&
       a0, xdv0, xka0, &
       de0, dep, de0p, &
       rk, x1, aa0, aa, &
       rad, p1, &
       feuneu, fquer, fqa, &
       dwsum, dmsum, dtsum, &
       resold, res, dres

  real (kind=dp), external :: diff_wat_vap
  real (kind=dp), external :: therm_conduct_air
  real (kind=dp), external :: p21  ! saturation vapour pressure
  real (kind=dp), external :: xl21 ! latent heat of vaporisation = f(temperature)
  real (kind=dp) :: zxl21          ! latent heat of vaporisation (local value for t)

! Local arrays
  real (kind=dp) :: cd(nkt,nka),cr(nkt,nka),sr(nkt,nka),falt(nkt,nka),c(nkt)
  real (kind=dp) :: psi(nkt),u(nkt)

! Common blocks:
  common /cb44/ a0m,b0m(nka)
  real (kind=dp) :: a0m,b0m

  common /cb49/ qabs(18,nkt,nka,jptaerrad), & ! only qabs is used here
                qext(18,nkt,nka,jptaerrad), &
                asym(18,nkt,nka,jptaerrad)
  real (kind=dp) :: qabs,qext,asym

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb51/ dlgew,dlgenw,dlne
  real (kind=dp) :: dlgew, dlgenw, dlne

! == End of declarations =======================================================

  zxl21  = xl21(to)
  xldcp  = zxl21/cp
  xka    = therm_conduct_air(to)
  xdv    = diff_wat_vap(to,pp)
  xl     = 24.483_dp *to/pp        ! P&K eq. (10-140) but T0=273.15, mistake?
  deltav = 1.3_dp * xl
  deltat = 2.7_dp * xl
  rho    = pp/(r0*to*(1._dp + 0.61_dp * xm1o))
  rho21  = p21(to)/(r1*to)
  rho21s = (zxl21/(r1*to)-1._dp)*rho21/to
  a0     = a0m/to
  xdv0   = xdv*sqrt(2._dp*pi/(r1*to)) / 3.6e-08_dp             ! factor for Dv* calculation, P&K (13-14)
  xka0   = xka*sqrt(2._dp*pi/(r0*to)) / (7.e-07_dp * rho * cp) ! factor for Ka* calculation, P&K (13-20)
  kr0    = kr

  if (totr(1).lt.1._dp) then
     ib0=7                  ! solar bands excluded, IR only
  else
     ib0=1
  end if

  do ia=1,nka
     do jt=1,nkt
        jtp  = min(jt+1,nkt)
        de0  = dew(jt)
        dep  = dew(jtp)
        de0p = de0 + dep

        rk = rw(jt,ia)
        sr(jt,ia) = max(0.1_dp, exp(a0 / rk - b0m(ia) * en(ia) / ew(jt)))
        xdvs = xdv / (rk / (rk + deltav) + xdv0 / rk)
        xkas = xka / (rk / (rk + deltat) + xka0 / rk)
        x1 = rhow * (zxl21 + xkas / (xdvs * rho21s * sr(jt,ia)))
        cd(jt,ia) = 3.e12_dp * rho21 * xkas / (x1 * rk * rk * rho21s * sr(jt,ia))
        if (kr0.eq.3.and.rn(ia).lt.0.5_dp) kr0=2 ! jjb need investigation, maybe a bug: once changed, will not change again to its original value
        rad=0._dp
        do ib=ib0,mb
           rad = rad + totr(ib) * (qabs(ib,jt,ia,kr0)  * de0 + &
                                   qabs(ib,jtp,ia,kr0) * dep) / de0p
        enddo
        cr(jt,ia) = rad * 7.5e5_dp / (rk * x1) - rhow * 4190._dp * (tn-to)/(dt*x1)
     enddo
  enddo
  falt(:,:) = ffk(:,:)

  feuneu = feualt + dfdt * dt
  if (feualt.lt.0.95_dp) then
     feuneu=xm1n*pp/(p21(tn)*(.62198_dp + .37802_dp * xm1n))
  end if
  fquer=0.5_dp * (feuneu+feualt)
  ! Initialisation
  res = 0._dp ! residue
  aa0 = 1._dp / dt

  ! condensation iteration
  do itk=1,10
     dwsum=0._dp
     do ia=1,nka
        do jt=1,nkt
           psi(jt)=falt(jt,ia)
           c(jt)=(cd(jt,ia)*(fquer-sr(jt,ia))-cr(jt,ia))/dlne
        enddo

        u(1) = max(0._dp,c(1))
        do jt=2,nkt-1
           u(jt) = 0.5_dp*(c(jt)+abs(c(jt))+c(jt-1)-abs(c(jt-1)))
        enddo
        u(nkt) = min(0._dp,c(nkt-1))

        call advec (dt,u,psi)

        do jt=1,nkt
           ffk(jt,ia)=psi(jt)
           dwsum=dwsum+(psi(jt)-falt(jt,ia))*e(jt)
        enddo
     enddo
     dmsum = dwsum/rho
     dtsum = xldcp*dmsum
     xm1o  = xm1n-dmsum
     to    = tn+dtsum
     p1    = xm1o*pp/(0.62198_dp + 0.37802_dp * xm1o)
     feuneu = p1/p21(to)
     resold = res
     res    = feuneu + feualt - 2._dp * fquer
     if (abs(res).lt.1.e-6_dp) return
! calculation of fquer by Newton-Interpolation
     dres = res-resold
     aa = aa0
     if (itk.gt.1.and.abs(dres).gt.1.e-8_dp) then
        aa=(fqa-fquer)/dres
     end if
     fqa = fquer
     fquer=fquer+aa*res
  enddo

  write (jpfunout,*)'warning SR subkon: no convergence of condensation iteration'

end subroutine subkon

!
!-------------------------------------------------------------
!

function diff_wat_vap(temperature,pressure)

  ! Diffusivity of water vapour in air
  ! for temperatures between -40 and +40 Celsius

  ! Dv = 0.211d-4 *(T/T0)**1.94 *(P0/P)

  ! where T0=273.15 K and P0=101325 Pa

  ! see Pruppacher and Klett, Microphysics of clouds and precipitations
  ! equation (13-3) p. 503, here expressed in MKS units


! Author
! ------
!     Josue Bock

! History
! -------
!     12-01-2017

! Declarations :
! ------------
! Modules used:

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunout

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  real (kind=dp) :: diff_wat_vap

  real (kind=dp), intent(in) :: temperature     ! in [K]
  real (kind=dp), intent(in) :: pressure        ! in [Pa]

  real (kind=dp), parameter :: cst=0.211e-4_dp  ! in [m2/s]
  real (kind=dp), parameter :: exponent=1.94_dp
  real (kind=dp), parameter :: T0=273.15_dp     ! in [K]
  real (kind=dp), parameter :: P0=101325_dp     ! in [Pa]

  real (kind=dp), parameter :: cst2=cst*P0/(T0**exponent)

! == End of declarations =======================================================

  diff_wat_vap = cst2 * temperature**exponent / pressure

  if (temperature < T0-40._dp .or. temperature > T0+40._dp) then
     write(jpfunout,*)'Warning, in FN diff_wat_vap, the parameterisation'
     write(jpfunout,*)'  is valid between -40 and +40 Celsius'
     write(jpfunout,*)'  The temperature is: ',temperature,' K.'
     write(jpfunout,*)'  The parameterisation has nonetheless been used'
     write(jpfunout,*)'  The calculated diffusivity is: ',diff_wat_vap,' m2/s'
  end if

end function diff_wat_vap

!
!-------------------------------------------------------------
!

function therm_conduct_air(temperature)

  ! Thermal conductivity of air

  ! ka = 1.d-3 * (4.39 + 0.071*T)

  ! where T is in K, and ka is in J/(m.s.K)
  ! see Seinfeld and Pandis 2nd Ed., equation (17.71) p.786


! Author
! ------
!     Josue Bock

! History
! -------
!     13-01-2017

! Declarations :
! ------------
! Modules used:

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  real (kind=dp) :: therm_conduct_air

  real (kind=dp), intent(in) :: temperature  ! in [K]

  real (kind=dp), parameter :: cst1 = 4.39e-3_dp
  real (kind=dp), parameter :: cst2 = 7.1e-5_dp

! == End of declarations =======================================================

  therm_conduct_air = cst1 + cst2*temperature


end function therm_conduct_air

!
!-------------------------------------------------------------
!

subroutine advec (dt,u,y)

! advection of particles for condensational growth

! one of many implementations of Bott's advection scheme
! specific for 2D particle spectrum


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  !

! == End of header =============================================================

! jjb work done
!    - removed arithmetic if for polynomial order (1, 2 or 4) at the end
!    - rewritten using up to date features: do, do while, exit and cycle
!    - also rewritten the initial search for min/max indexes by reading/writting arrays in increasing indexes, for computing efficiency
!    - removed archaic forms of Fortran intrinsic functions
!    - passed y and u as arguments instead of common block (formerly cb59)

  USE config, ONLY : &
! Imported Routines:
       abortM

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunerr

  USE global_params, ONLY : &
! Imported Parameters:
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  real (kind=dp), intent(in) :: dt        ! fractional timestep
  real (kind=dp), intent(in) :: u(nkt)    ! advection velocity
  real (kind=dp), intent(inout) :: y(nkt) ! advected quantity

  real (kind=dp), parameter :: ymin=1.e-32_dp ! 1.e-8_dp in Bott

  integer :: i0, i1    ! lowest and highest bin indexes where y(i) >= ymin
  integer :: i         ! loop index
  integer :: k, k1, k2 ! integer jump (k) and its previous values when iterated
  integer :: k_low, k_high

  real (kind=dp) :: a0, a1, a2, a3, a4
  real (kind=dp) :: al, al2, al3
  real (kind=dp) :: c0
  real (kind=dp) :: dt0, dt1
  real (kind=dp) :: x0       ! floating position, in the nkt-grid, between k-1 and k+1
  real (kind=dp) :: x1       ! polynomial advection of fractional part
  real (kind=dp) :: z(nkt)

! == End of declarations =======================================================

! set input values to z, initialize output y
  z(:) = y(:)
  y(:) = 0._dp

! find the lowest (=i0) particle class that contains enough (>=ymin) particles
  i0 = 1
  do while (z(i0) .lt. ymin)
     ! return if not enough particles in all nkt classes
     if (i0 .eq. nkt) return
     i0 = i0 + 1
  end do

! find highest (=i1) particle class that contains enough (>=ymin) particles
  i1 = nkt
  do while (z(i1) .lt. ymin)
     i1 = i1 - 1
  end do


! do growth/advection calculation only in between lowest and highest bin with sufficient particles
  iloop: do i=i0,i1

     ! if not enough particles in the current bin, skip
     if (z(i).lt.ymin) cycle iloop

     k2  = 0
     dt0 = 0  ! this one does not need to be initialised, just prevent forcheck to complain [313 I]
     dt1 = dt
     k   = i

     if(abs(u(k))>0._dp) then
        dt0=min(1._dp/(abs(u(k))),dt1)
     else ! u(k) = 0.d0
        y(k)=y(k)+z(i)
        cycle iloop
     end if
     x0=real(k,dp)+u(k)*dt0         ! from dt0 definition above, x0 will lie in the [k-1 ; k+1] interval
     dt1=dt1-dt0
     k1 = k

! calculate integer jump of droplet class
     do while (dt1 .gt. 1.e-7_dp)
        if (u(k) < 0._dp) then
           k = k-1             ! evaporation
        else                       ! u(k) > 0, the case u(k)==0 has already be excluded
           k = k+1             ! condensation
        end if

        ! u change of sign between two adjacent class
        if (k==k2) then
           ! jjb maybe something needed here. Amongst both options (k2=k or k1) should always k2 be the one chosen?
           y(k)=y(k)+z(i)
           cycle iloop ! cycle external i loop
        end if

        ! save the last two values of k
        k2 = k1
        k1 = k

        if(abs(u(k))>0._dp) then
           dt0=min(1._dp/(abs(u(k))),dt1)
        else ! u(k) = 0.d0
           y(k)=y(k)+z(i)
           cycle iloop
        end if
        x0=real(k,dp)+u(k)*dt0
        dt1=dt1-dt0

     end do

     k_low  = int(floor(x0))
     k_high = k_low + 1

     c0 = x0 - real(k_low,dp) ! = x0 - floor(x0)

     ! Check these final indexes
     !  (k_high is actually used only if c0 > 0.)
     if (k_low .lt. 1 .or. (k_high .gt. nkt .and. c0 > 0._dp) ) then
        write(jpfunerr,*)i,k_low,k_high,x0,c0,dt
        write(jpfunerr,*)u(max(1,k_low):min(nkt,k_high))
        call abortM('SR advec: error with k_high or k_low')
     end if

     ! General case: polynomial advection for fractional part of courant number
     if (c0 > 0._dp) then

        if (i==1 .or. i==nkt) then
           ! at first and last grid point (i=1,nkt): first order polynomial
           x1 = c0 * z(i)
        else if (i==2 .or. i==nkt-1) then
           ! at second and second last grid point (i=2,nkt-1): second order polynomial
           al  = 1._dp - 2._dp * c0
           al2 = al * al
           a0  = (26._dp * z(i) - z(i+1) - z(i-1)) / 24._dp
           a1  = (z(i+1) - z(i-1)) / 16._dp
           a2  = (z(i+1) + z(i-1) - 2._dp * z(i)) / 48._dp
           x1  = min(z(i), a0*c0 + a1*(1._dp-al2) + a2*(1._dp-al2*al))
        else
           ! i=3,nkt-2: fourth order polynomial
           al  = 1._dp - 2._dp * c0
           al2 = al  * al
           al3 = al2 * al
           a0  = (9._dp*(z(i+2)+z(i-2))-116._dp*(z(i+1)+z(i-1))+2134._dp*z(i)) / 1920._dp
           a1  = (-5._dp*(z(i+2)-z(i-2))+34._dp*(z(i+1)-z(i-1))) / 384._dp
           a2  = (-z(i+2)+12._dp*(z(i+1)+z(i-1))-22._dp*z(i)-z(i-2)) / 384._dp
           a3  = (z(i+2)-2._dp*(z(i+1)-z(i-1))-z(i-2)) / 768._dp
           a4  = (z(i+2)-4._dp*(z(i+1)+z(i-1))+6._dp*z(i)+z(i-2)) / 3840._dp
           x1  = min(z(i),a0*c0+a1*(1._dp-al2)+a2*(1.d0-al3) &
                           +a3*(1._dp-al2*al2)+a4*(1._dp-al2*al3))
        end if

        x1 = max(0._dp,x1)

        y(k_low)  = y(k_low) + z(i) - x1
        y(k_high) = y(k_high)+ x1

     else ! c0 == 0.d0
        ! This case will happen in two situations:
        !   - if u(k)*dt0 is exactly = +/- 1.000d0 (this really happen sometimes)
        !   - if u(k)*dt0 is so small that it has been numerically rounded to 0.
        !     (very unlikely, would imply that u(k) ~ tiny(0.d0) and dt1 ~ [1.1e-7 -- 1.e-1] )

        ! In this case, computing polynomial expression is not necessary, it would
        ! lead to x1 = 0
        y(k_low)=y(k_low)+z(i)

     end if

  end do iloop

end subroutine advec

!
!-------------------------------------------------------------
!

subroutine advsed0 (c, y)

! Vertical advection for sedimentation, upstream procedure


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  !

! == End of header =============================================================

!     jjb removed declaration of 6 unused variables
! changed internal parameter n=nf, misleading since in most places n=nf+50
! (obviously, led to a bug when including erroneously "n" from global_params

! jjb cleaning of unused varaibles 15/12/16
! jjb declaration of all variables, implicit none 15/12/16
! jjb fortran generic functions min and max 15/12/16
! jjb checked identical to AB str code

  ! 03-May-2021  Josue Bock  Merge with latest Mistra version from A.Bott:
  !                           - remove cb58, use arguments instead

  USE global_params, ONLY : &
! Imported Parameters:
       nf

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
  real (kind=dp), intent(in)    :: c(nf)
  real (kind=dp), intent(inout) :: y(nf)

! Local scalars:
  integer :: i                      ! running index
! Local arrays:
  real (kind=dp) :: fm(nf),fp(nf)   ! advection fluxes

! == End of declarations =======================================================

  do i=1,nf-1
     fm(i) = -min(0._dp, c(i)) * y(i+1)
     fp(i) =  max(0._dp, c(i)) * y(i)
  end do
  do i=2,nf-1
     y(i) = y(i) - fm(i-1) + fp(i-1) + fm(i) - fp(i)
  end do

end subroutine advsed0

!
!-------------------------------------------------------------
!

subroutine advsed1 (c, y)

! Vertical advection of quantity psi with the positive definite advection


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  !

! == End of header =============================================================

! area preserving flux form; Bott (1989): Monthly Weather Review.
! fourth order monotone version.
! y(i) is transport quantity, input and output.
! boundary conditions are y(1)=const, y(n)=const.
! c(i) is Courant number satisfying the CFL criterion, input.
! fm(i), fp(i) are fluxes for u(i)<0 and u(i)>0, respectively.
! a0, a1, a2, a3, a4 are coefficients of polynomials in gridbox i.
! At i=1 and i=n first order polynomial,
! at i=2 and i=n-1 second order polynomial,
! at 3<=i<=n-2 fourth order polynomial.
! w(i) are weighting factors.
! the numerical grid is equidistant.
! the procedure is one dimensional.
! for multidimensional applications time splitting has to be used.
! the quantities c(i), fm(i), fp(i)  are given at the right
! boundary of grid cell i.
! Thus, fm(i) is flux from gridbox i+1 into gridbox i for c(i)<0,
! fp(i) is flux from gridbox i into gridbox i+1 for c(i)>0.

  USE global_params, ONLY : &
! Imported Parameters:
       nf

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
  real (kind=dp), intent(in)    :: c(nf)
  real (kind=dp), intent(inout) :: y(nf)

! Local scalars:
  integer :: i                      ! running index
  real (kind=dp) :: cl, clm
  real (kind=dp) :: fmim
  real (kind=dp) :: w, x1, x2, x3
  real (kind=dp) :: ymin, ymax
! Local arrays:
  real (kind=dp) :: a0(2:nf-1),a1(2:nf-1),a2(2:nf-1),a3(2:nf-1),a4(2:nf-1)
  real (kind=dp) :: fm(nf-1)

! == End of declarations =======================================================

  a0(2)=(26._dp * y(2) - y(3) - y(1)) / 24._dp
  a1(2)=(y(3) - y(1)) / 16._dp
  a2(2)=(y(3) + y(1) - 2._dp * y(2)) / 48._dp
  a3(2)=0._dp
  a4(2)=0._dp
  do i=3,nf-2
     a0(i)=(9._dp*(y(i+2)+y(i-2))-116._dp*(y(i+1)+y(i-1)) &
              +2134._dp*y(i))/1920._dp
     a1(i)=(-5._dp*(y(i+2)-y(i-2))+34._dp*(y(i+1)-y(i-1)))/384._dp
     a2(i)=(-y(i+2)+12._dp*(y(i+1)+y(i-1))-22._dp*y(i)-y(i-2))/384._dp
     a3(i)=(y(i+2)-2._dp*(y(i+1)-y(i-1))-y(i-2))/768._dp
     a4(i)=(y(i+2)-4._dp*(y(i+1)+y(i-1))+6._dp*y(i)+y(i-2))/3840._dp
  enddo
  a0(nf-1)=(26._dp*y(nf-1)-y(nf)-y(nf-2))/24._dp
  a1(nf-1)=(y(nf)-y(nf-2))/16._dp
  a2(nf-1)=(y(nf)+y(nf-2)-2._dp*y(nf-1))/48._dp
  a3(nf-1)=0._dp
  a4(nf-1)=0._dp
  cl=-c(nf-1)
  fm(nf-1)=min(y(nf),cl*(y(nf)-(1._dp-cl)*(y(nf)-y(nf-1))*0.5_dp))
  clm=cl
  do i=nf-1,2,-1
     cl=clm
     clm=-c(i-1)
     x1 = 1._dp - 2._dp * cl
     x2 = x1 * x1
     x3 = x1 * x2
     ymin=min(y(i),y(i+1))
     ymax=max(y(i),y(i+1))
     fmim=max(0.d0,a0(i)*cl-a1(i)*(1._dp-x2)+a2(i)*(1._dp-x3) &
                   -a3(i)*(1._dp-x1*x3)+a4(i)*(1._dp-x2*x3))
     fmim=min(fmim,y(i)-ymin+fm(i))
     fmim=max(fmim,y(i)-ymax+fm(i))
     fmim=max(0._dp,fmim-(cl-clm)*y(i))
     w=y(i)/max(fmim+1.e-15_dp,y(i))
     fm(i-1)=fmim*w
  enddo

! new values of y
  y(1) = y(1) + fm(1)
  do i=2,nf-1
     y(i) = y(i) - fm(i-1) + fm(i)
  enddo
  y(nf) = y(nf) - fm(nf-1)

end subroutine advsed1

!
!-------------------------------------------------------------
!

subroutine advseda (c, y)
! Vertical advection of quantity y with the positive definite advection
! scheme after Bott (1989).


! Author:
! ------
  !    Andreas Bott


! Modifications :
! -------------
  ! jjb: currently unused, see with ABott if advseda should replace advsed1

! == End of header =============================================================



  USE global_params, ONLY : &
! Imported Parameters:
       nf

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
  real (kind=dp), intent(in)    :: c(nf)
  real (kind=dp), intent(inout) :: y(nf)

  integer :: i
  real (kind=dp) :: a0, a1, a2, a3, a4
  real (kind=dp) :: cl, fm, w
  real (kind=dp) :: x1, x2, x3
  real (kind=dp) :: flux(nf)

! == End of declarations =======================================================

! initialisation
  flux(:) = 0._dp

! flux of y in the first level
  a0 = (26._dp*y(2)-y(3)-y(1))/24._dp
  a1 = (y(3)-y(1))/16._dp
  a2 = (y(3) + y(1) - 2._dp * y(2)) / 48._dp
  cl = -c(1)
  x1 = 1._dp - 2._dp * cl
  x2 = x1*x1
  fm = max(0._dp,a0*cl-a1*(1._dp-x2)+a2*(1._dp-x1*x2))
  w  = y(2)/max(fm+1.d-15,a0+2._dp*a2)
  flux(1) = fm*w

! loop for flux of y in the fog levels
  do i=3,nf-2
     a0 = (9._dp*(y(i+2)+y(i-2))-116._dp*(y(i+1)+y(i-1))+2134._dp*y(i))/1920._dp
     a1 = (-5._dp*(y(i+2)-y(i-2))+34._dp*(y(i+1)-y(i-1)))/384._dp
     a2 = (-y(i+2)+12._dp*(y(i+1)+y(i-1))-22._dp*y(i)-y(i-2))/384._dp
     a3 = (y(i+2)-2._dp*(y(i+1)-y(i-1))-y(i-2))/768._dp
     a4 = (y(i+2)-4._dp*(y(i+1)+y(i-1))+6._dp*y(i)+y(i-2))/3840._dp
     cl = -c(i-1)
     x1 = 1._dp - 2._dp * cl
     x2 = x1 * x1
     x3 = x1 * x2
     fm = max(0._dp, a0*cl-a1*(1._dp-x2)+a2*(1._dp-x3)-a3*(1._dp-x1*x3) &
                    +a4*(1._dp-x2*x3))
     w = y(i) / max(fm+1.d-15, a0+2._dp*(a2+a4))
     flux(i-1) = fm * w
  end do

! flux of y in the second highest level
  a0 = (26._dp*y(nf-1)-y(nf)-y(nf-2))/24._dp
  a1 = (y(nf)-y(nf-2))/16._dp
  a2 = (y(nf)+y(nf-2)-2._dp*y(nf-1))/48._dp
  cl = -c(nf-2)
  x1 = 1._dp - 2._dp*cl
  x2 = x1*x1
  fm = max(0._dp,a0*cl-a1*(1._dp-x2)+a2*(1._dp-x1*x2))
  w=y(nf-1)/ max(fm+1.e-15_dp,a0+2._dp*a2)
  flux(nf-2) = fm*w

! flux of y in the highest level
  cl=-c(nf-1)
  flux(nf-1)=dmin1(y(nf),cl*(y(nf)-(1._dp-cl)*(y(nf)-y(nf-1))*0.5_dp))

! new values of y
  y(1)=y(1)+flux(1)
  do i=2,nf-1
     y(i)=y(i)-flux(i-1)+flux(i)
  end do
  y(nf)=y(nf)-flux(nf-1)

end subroutine advseda


!
!------------------------------------------------------------------------
!

subroutine stem_kpp (dd,xra,z_box,n_bl,box,chamber,nuc)
! chemical reactions
! aerosol mass change due to chemical reactions
!
! NOTE: the aerosol processing as implemented here is mass conserving
!       but shifts too many particles into the smallest bins, therefore
!       the particle spectra are a bit odd
!       this has to be improved in future versions


! Author:
! ------
  !    Roland von Glasow


! Modifications :
! -------------
  ! jjb bugfix: added nuc in argument list, needed in two if tests. Also added ifeed (from /nucfeed/, then from config)

! == End of header =============================================================

  USE config, ONLY : &
       ifeed,        &
       nkc_l

  USE constants, ONLY : &
! Imported Parameters:
       pi

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunout

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j6, &
       nf, &
       n, &
       nka, &
       nkt, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  real (kind=dp), intent(in) :: dd, xra, z_box
  integer, intent(in) :: n_bl
  logical, intent(in) :: box, chamber, nuc

  integer, parameter :: lsp=9
  real (kind=dp), parameter :: fpi = 4._dp / 3._dp * pi

  integer :: ia, ial, iau, iend, iia, iinkr, istart, ix, jt, jtl, jtu, k, kc, kkc, l, ll
  integer :: nfrom, nto, nmin, nmax
  integer tix,tixp
  logical :: llboxch
  real (kind=dp) :: c0, den, vol_ch, x0, x1, xch, xfact

! Common blocks:
  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /blck06/ kw(nka),ka
  integer :: kw, ka
  common /blck12/ cw(nkc,n),cm(nkc,n)
  real (kind=dp) :: cw, cm
  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1

  integer lj2(lsp)
  real (kind=dp) sion1o(lsp,nkc,n) ! sion1o old ion conc. [mole m**-3]
  real (kind=dp) dsion1(lsp,nkc,n) ! dsion1 change in ion conc. [mole/part.]
      ! jjb fs is not used over all its dimensions (only in an old line commented)
      ! maybe worth to reduce its dimensions after double check, to save cpu time
  real (kind=dp) fs(nka,nkc,n)
  real (kind=dp) sap(nkc,n) ! sap total number of aerosols [cm**-3]
  real (kind=dp) smp(nkc,n) ! smp total aerosol mass [mg cm**-3]
  real (kind=dp) vc(nkc,nkc,n)
  data lj2/1,2,8,9,13,14,19,20,30/

! == End of declarations =======================================================

  nmin = 2
  nmax = nf
  if (box) then
     nmin = n_bl
     nmax = n_bl
  else if (chamber) then
     nmax = n_bl
  endif

  if (.not.chamber) call aer_source (box,dd,z_box,n_bl)

  llboxch = box.or.chamber
  call liq_parm (xra,llboxch,n_bl)

!*************************** no aerosol processing *****************
!   call kpp_driver (box,dd,n_bl)
!   return
!*************************** no aerosol processing *****************

! sion1 ion conc. [mole m**-3], sion1o old ion conc. [mole m**-3],
! dsion1 change in ion conc. [mole/part.]
! sap total number of aerosols [cm**-3], smp total aerosol mass [mg cm**-3]
! fs total aerosol mass integrated over all liquid classes [mg cm**-3]
! lj2 species that define the aerosol mass

!      if (.not.box) then

! initialize variables
  vc(:,:,:) = 0._dp
  sap(:,:) = 0._dp
  smp(:,:) = 0._dp

! loop over aqueous chemistry layers to get values of sion1, fs, smp, sap
! before chemistry integration (at t=t_0)
  do k=nmin,nmax
     do kc=1,nkc_l
        if (cm(kc,k).eq.0._dp) goto 900
!       define upper and lower limits of ia loop
        if (kc.eq.1.or.kc.eq.3) then
!           if ((nuc).and.(ifeed.eq.1)) then ! jjb corrected (?) below CHECK WITH Susanne PECHTL
!              ial=1
!           else
!              ial=1+1
!           endif
           if ((nuc).and.(ifeed.eq.2)) then
              ial=2
           else
              ial=1
           endif
           iau=ka
        else
           ial=ka+1
           iau=nka-1
        endif
!       loop over all ia that are in the current kc bin
        do ia=ial,iau
           fs(ia,kc,k) = 0._dp
!          define upper and lower limits of jt loop
           if (kc.eq.1.or.kc.eq.2) then
              jtl = 1
              jtu = kw(ia)
           else
              jtl = kw(ia)+1
              jtu = nkt
           endif
!          loop over all jt that are in the current kc bin
           do jt=jtl,jtu
              fs(ia,kc,k) = fs(ia,kc,k) + ff(jt,ia,k) * en(ia)
              sap(kc,k) = sap(kc,k) + ff(jt,ia,k)
           enddo
           smp(kc,k) = smp(kc,k) + fs(ia,kc,k)
        enddo
        do l=1,lsp
           ll=lj2(l)
           sion1o(l,kc,k)=sion1(ll,kc,k)
        enddo
900     continue
     enddo   ! kc
  enddo      ! k

!      endif ! .not.box
! chemistry SR

  call kpp_driver (box,dd,n_bl)

!      if (box) return

! redistribution of particles along aerosol grid due to modified aerosol mass
  do k=nmin,nmax
     do kc=1,nkc_l
        if (cm(kc,k).eq.0._dp) goto 1000
        if (sap(kc,k).gt.1.e-6_dp) then
!          change in mass determining chemical species
           do l=1,lsp
              ll=lj2(l)
              dsion1(l,kc,k)=(sion1(ll,kc,k)-sion1o(l,kc,k))*1.e-06/sap(kc,k)
           enddo
! den: new aerosol mass in mg/particle due to chemical reactions
! mole masses for l=1,2,8,9,13,14,19,20,30: H+=1g/mole, NH4=18g/mole, SO4(2-)=96g/mole,
! HCO3-=61g/mole, NO3=62g/mole, Cl-=35.5g/mole, HSO4-=97g/mole, Na+=23g/mole, CH3SO3-=95g/mole
! HCO3-=61g/mole --> 44 g/mole as water remains in particle when CO2 degasses due to acidification
           den=(dsion1(1,kc,k)*1.+dsion1(2,kc,k)*18.+dsion1(3,kc,k) &
!                *96.+dsion1(4,kc,k)*61.+dsion1(5,kc,k)*62.+ &
                *96.+dsion1(4,kc,k)*44.+dsion1(5,kc,k)*62.+ &
                dsion1(6,kc,k)*35.5+dsion1(7,kc,k)*97.+ &
                dsion1(8,kc,k)*23.+dsion1(9,kc,k)*95.)*1000.

!          define upper and lower limits of ia loop
           if (kc.eq.1.or.kc.eq.3) then
              if ((nuc).and.(ifeed.eq.2)) then
                 ial=2
              else
                 ial=1
              endif
              iau=ka
           else
              ial=ka+1
              iau=nka-1
           endif

!          loop over all ia that are in the current kc bin
!          c0: courant number for redistribution
! if growing reverse loop order to avoid increasing the mass of some
! particles twice
           istart=ial
           iend=iau
           iinkr=1
           if (den.ge.0._dp) then
              istart=iau
              iend=ial
              iinkr=-1
           endif
!           do ia=ial,iau
           do ia=istart,iend,iinkr
              if (den.gt.0._dp) then
                 x0=en(ia)+den*en(ia)/smp(kc,k)*sap(kc,k)
              else
!                 x0=en(ia)+den*fs(ia,kc,k)/smp(kc,k)*sap(kc,k)
                 x0=en(ia)+den*en(ia)/smp(kc,k)*sap(kc,k)
!                 x0=dmax1(en(ia)+den,0.d0)
                 if (x0.le.0._dp) write(jpfunout,*) k,kc,ia,'aerosol growth'
              endif
              do iia=1,nka-1
                 if (en(iia).le.x0.and.en(iia+1).gt.x0) then
                    ix=iia
                    c0=(en(iia+1)-x0)/(en(iia+1)-en(iia))
!                    c0=(log10(en(iia+1))-log10(x0))/dlgenw
                    go to 2000
                 endif
              enddo
              if (en(1).gt.x0) then
                 ix=1
                 c0=1._dp
              else
                 ix=nka-1
                 c0=0._dp
              endif
2000          continue
!             define upper and lower limits of jt loop
              if (kc.eq.1.or.kc.eq.2) then
                 jtl=1
                 jtu=kw(ia)
              else
                 jtl=kw(ia)+1
                 jtu=nkt
              endif
!             loop over all jt that are in the current kc bin
              do jt=jtl,jtu
                 if (ff(jt,ia,k).gt.0._dp) then
                    x1=ff(jt,ia,k)
                    ff(jt,ia,k)=0._dp
                    ff(jt,ix,k)=ff(jt,ix,k)+x1*c0
                    ff(jt,ix+1,k)=ff(jt,ix+1,k) + x1*(1._dp - c0)
! find "targetbin" for ix and ix+1:
! ix
                    if (ix.gt.ka) then
                       if (jt.gt.kw(ix)) then
                          tix=4
                       else
                          tix=2
                       endif
                    else
                       if (jt.gt.kw(ix)) then
                          tix=3
                       else
                          tix=1
                       endif
                    endif
! ix+1
                    if (ix+1.gt.ka) then
                       if (jt.gt.kw(ix+1)) then
                          tixp=4
                       else
                          tixp=2
                       endif
                    else
                       if (jt.gt.kw(ix+1)) then
                          tixp=3
                       else
                          tixp=1
                       endif
                    endif
! store change of volume for ix --> ix and ia --> ix+1
                    if (tix.ne.kc)  vc(tix,kc,k) =vc(tix,kc,k) +x1*c0*fpi*rq(jt,ia)**3
                    if (tixp.ne.kc) vc(tixp,kc,k)=vc(tixp,kc,k)+x1*(1._dp-c0)*fpi*rq(jt,ia)**3
                 endif
              enddo   ! jt loop
           enddo      ! ia loop
        endif

1000    continue
     enddo  ! kc loop
  enddo     ! k  loop

! move chemical species if transport above chemistry bins took place
  do k=nmin,nmax
     do kc=1,nkc_l
        do kkc=1,nkc_l
           if (kkc.eq.kc) goto 1001
! "kc" (from) bin and "kkc" (to) bin, i.e. bin that loses moles and bin that gains moles
           if (vc(kkc,kc,k).eq.0._dp) goto 1001
           if (cw(kc,k).gt.0._dp) then
              nfrom=kc
              nto=kkc
! cw in m^3(aq)/m^3(air), vc in um^3/cm^3: 10^-12
              vol_ch=vc(kkc,kc,k)*1.d-12
              xfact=0._dp
!            if (vol_ch.gt.1.e-4*cw(k,nfrom)) then
              xfact=1._dp-(cw(nfrom,k)-vol_ch)/cw(nfrom,k)
              do l=1,j2
                 xch=sl1(l,nfrom,k)*xfact
                 sl1(l,nfrom,k)=sl1(l,nfrom,k)-xch
                 sl1(l,nto,k)  =sl1(l,nto,k) +xch
              enddo
              do l=1,j6
                 xch=sion1(l,nfrom,k)*xfact
                 sion1(l,nfrom,k)=sion1(l,nfrom,k)-xch
                 sion1(l,nto,k)  =sion1(l,nto,k) +xch
              enddo
           end if
!            endif  ! if vol_ch > 10^-4* cw
1001       continue
!         enddo  ! ic loop
        enddo               ! kkc loop
! 1002       continue ! jjb statement label unreferenced
     enddo                  !kc loop
  enddo     ! k loop

end subroutine stem_kpp

!
!----------------------------------------------------------------------
!

subroutine adjust_f
! adjustment of initial aerosol size distribution for specific scenarios


! Author:
! ------
  !    RvG?


! Modifications :
! -------------
  !

! == End of header =============================================================

  USE constants, ONLY : &
!! Imported Parameters:
!       pi,              &
       rho3,            &       ! aerosol density
       rhow                     ! water density

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nka, &
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local parameters:
  ! optimisation: define parameters that will be computed only once
  real (kind=dp), parameter :: zrho_frac = rho3 / rhow
  !real (kind=dp), parameter :: z4pi3 = 4.e-09_dp * pi / 3._dp

  integer :: ia, k, kl
  real (kind=dp) :: a0, b0, x0, rg !,eg
  real (kind=dp), external :: rgl
  real (kind=dp) :: f_inter(nka)

! Common blocks:
  common /cb44/ a0m,b0m(nka)
  real (kind=dp) :: a0m,b0m
  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq
  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho

! == End of declarations =======================================================

  x0 = 1._dp

  do k=2,n
     do ia=1,nka
        f_inter(ia) = ff(1,ia,k)
     enddo

     fsum(k)=0._dp
     do ia=nka,1,-1
        if (rn(ia).gt.0.5_dp) x0 = 0.1_dp
! init spectrum:aerosols in equilibrium with rH=76.2 % (below: last argument in rgl)
! "dry" them !(this is not really exact..)
! equilibrium radius at 76 %
        a0 = a0m / t(k)
        b0 = b0m(ia) * zrho_frac
! b0=b0m*rho3/rhow; rho3=2000; rhow=1000
        rg = rgl(rn(ia),a0,b0,.762_dp)
!        eg = z4pi3 * (rg**3 - rn(ia)**3)
        do kl=1,nka
           if (rg.le.rn(kl)) then
!           if (eg.le.enw(kl)) then
              ff(1,ia,k) = x0 * f_inter(kl)
              exit
           endif
        enddo
        fsum(k) = fsum(k) + ff(1,ia,k)
     enddo

  enddo

end subroutine adjust_f

!----------------------------------------------------------------------

subroutine partdep (ra)
! calculate particle dry deposition velocity after Seinfeld and Pandis,
! 1998, p.958ff


! Author:
! ------
  !    RvG?


! Modifications :
! -------------
  !

! == End of header =============================================================

  USE constants, ONLY : &
! Imported Parameters:
       g, &                        ! Gravitational acceleration (m/s**2)
       pi,&
       kappa

  USE data_surface, ONLY : &
       ustern, z0                  ! frictional velocity, roughness length

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nka, &
       nkt, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  real (kind=dp), intent(out) :: ra

  real (kind=dp), external :: vterm

  integer :: ia,jt, k, kc
  real (kind=dp) :: cc, phi, rb, rx, Sc, St, vs
  real (kind=dp) :: xd, xeta, xk, xlam, xnu, z
  real (kind=dp) :: xx1(nkc)

! Common blocks:
  common /blck06/ kw(nka),ka
  integer :: kw, ka
  common /blck12/ cw(nkc,n),cm(nkc,n)
  real (kind=dp) :: cw, cm
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /kinv_i/ kinv
  integer :: kinv
  common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
  real (kind=dp) :: vt, vd, vdm
!      common /kpp_kg/ vol2(nkc,n),vol1(n,nkc,nka),part_o &
!     &     (n,nkc,nka),part_n(n,nkc,nka),pntot(nkc,n),kw(nka),ka
!      real (kind=dp) :: ra,rb,vs,vd,z,xD,sc,st,rx,Cc

! == End of declarations =======================================================

  call monin (phi)

! particle dry deposition velocity:v_d=1/(ra + rb + ra rb v_s)+ v_s
!    ra=1/(kappa ustar) (ln (z/z0) +Phi) ;where z=height of surface (constant flux) layer
!                                         Phi takes stratification into account
!    rb=1/(ustar(Sc^(-2/3) + 10^(-3/St))) ;Sc=nu/D  St=v_s*ustar^2/(g nu)

  xk=1.38066d-23 !Boltzmann number
  z=0.1*eta(kinv) !surface layer height: 10% of BL (Stull), insensitive parameter
  ra=1./(kappa*ustern)*(log(z/z0)+phi)  !ra:aerodynamic resistance; kappa=0.4

  k=2 ! only in lowest model layer
  xeta=1.8325e-5*(416.16/(t(k)+120._dp))*((t(k)/296.16)**1.5) !dynamic viscosity of air, Jacobsen p. 92
  xnu=xeta/rho(k)  !kinematic viscosity of air
!  write (110,*) 'nu eta',xnu,xeta

! set xx1(kc)=0.
  do kc=1,nkc
     xx1(kc) = 0.d0
     vdm(kc) = 0.d0 ! jjb added this initialisation. If not, old values are still used when cw goes to 0 from t to t+1
  enddo
! free path length
  xlam = 2.28e-5 * t(k) / p(k)

  do ia=1,nka
     do jt=1,nkt
        rx=rq(jt,ia)*1.e-6
        vs=vterm(rx,t(k),p(k)) ! Stokes fall velocity incl Cc
        Cc=1.+xlam/rx*(1.257_dp+.4*exp(-1.1*rx/xlam)) ! Cunningham slip flow corr.
        xD=xk*t(k)*Cc/(6*pi*xeta*rx) ! aerosol diffusivity
        Sc=xnu/xD !Schmidt number
        St=vs*ustern**2/(g*xnu) !Stokes number
        rb=1./(ustern*(sc**(-2./3.)+10**(-3./st))) !quasi laminar resistance
        vd(jt,ia)=1./(ra+rb+ra*rb*vs)+vs !deposition velocity
        !        write (110,20) ia,jt,rq(jt,ia),vs*100.,vd*100.,100./(ra+rb+ra*rb*vs)

! calculate mass weighted mean dry deposition velocities
! loop over the nkc different chemical bins
        do kc=1,nkc
           if (cw(kc,k).eq.0.) goto 1001
! define kc limits  ---
           if (kc.eq.1.and.(ia.gt.ka.or.jt.gt.kw(ia))) goto 1001
           if (kc.eq.2.and.(ia.lt.(ka+1).or.jt.gt.kw(ia))) goto 1001
           if (kc.eq.3.and.(ia.gt.ka.or.jt.lt.(kw(ia)+1))) goto 1001
           if (kc.eq.4.and.(ia.lt.(ka+1).or.jt.lt.(kw(ia)+1))) goto 1001
! LWC weighted deposition velocity
           xx1(kc)=xx1(kc)+rx*rx*rx*vd(jt,ia)*ff(jt,ia,k)*1.e6
! deposition velocity:
!!              vdm(kc)=4.*3.1415927/(3.*cw(kc,k))*xx1(kc)
!              vdm(kc)=4.*pi/(3.*cw(kc,k))*xx1(kc)
1001       continue
        enddo               ! kc
     enddo
  enddo

  do kc=1,nkc
     if (cw(kc,k).gt.0.d0) &
     &              vdm(kc)=4.*pi/(3.*cw(kc,k))*xx1(kc) ! don't make this calculation nka * nkt * nkc times!
  end do

!      do kc=1,nkc
!         rx=rc(2,kc)
!         if (rx.gt.0.) vs=vterm(rx,t(2),p(2))
!         write (110,21) kc,vdm(kc),vs,rc(2,kc)
!      enddo
! 20   format (1p,2i3,5d16.8)
! 21   format (1p,i3,3d16.8)

end subroutine partdep

!
!-------------------------------------------------------
!

subroutine monin (phi)
! calculate the Monin-Obukhov length after Seinfeld and Pandis, 1998, p.862


! Author:
! ------
  !    RvG?


! Modifications :
! -------------
  !

! == End of header =============================================================

! jjb work done
!     - integer exponents (**2 instead of **2.)
!     - removed archaic forms of intrinsic functions (dlog, datan)

  USE constants, ONLY : &
! Imported Parameters:
       cp,&              ! Specific heat of dry air, in J/(kg.K)
       g, &                        ! Gravitational acceleration (m/s**2)
       kappa                       ! von Karman constant

  USE data_surface, ONLY : &
       ustern, z0                  ! frictional velocity, roughness length

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  real (kind=dp), intent(out) :: phi ! see S & P 1st Ed, p. 963, equation (19.14)

  integer :: k
  real (kind=dp) :: dtdz, q3, xeta, xeta0, xmo, z, zeta, zeta0

! Common blocks:
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
  real (kind=dp) :: atke, atkh, atkm, tke, tkep, buoy

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /kinv_i/ kinv
  integer :: kinv

! == End of declarations =======================================================

! check inversion height - it is diagnosed in SR atk1 but might be zero after restart
  if (kinv.eq.0)  then
     kinv = 70
     print *,'SR monin: kinv = 0., set to kinv=70'
  endif

! reference height = 10 % of BL height
  z=0.1*eta(kinv)
  do k=1,kinv
     if (eta(k).ge.z) goto 100
  enddo
100 continue

! L=-rho c_p T_0 Ustar^3/(kappa g \bar(q_3))
! with \bar(q_3)=rho c_p \bar(w'theta')=rho c_p (-1.) atkh d theta/d z

! dtheat'/dz in height k
  if (k.eq.1) then
     k=2
     print *,'SR monin: index out of bounds (1)'
  endif
  if (k.eq.n) then
     k=n-1
     print *,'SR monin: index out of bounds (2)'
  endif
  dtdz=((theta(k+1)-theta(k))/deta(k)+(theta(k)-theta(k-1))/deta(k-1))/2.
  q3=rho(k)*cp*(-1.)*atkh(k)*dtdz

  xmo=-1.*rho(k)*cp*t(1)*ustern**3/(kappa*g*q3) ! Seinfeld 2, p. 747, (16.70)

! effect on ra
  zeta=z/xmo
  zeta0=z0/xmo
! |L|>10^5  (neutral) => phi=0.
  if (abs(xmo) > 1.d5) then
     phi=0.
  else
! stable
     if (xmo.gt.0.) then
        phi=4.7*(zeta-zeta0)
! unstable
     else if (xmo.lt.0) then
        xeta0=(1._dp - 15.*zeta0)**0.25
        xeta=(1._dp - 15.*zeta)**0.25
        phi=log( (xeta0**2+1.)*(xeta0+1.)**2 /((xeta**2+1.)*(xeta+1.)**2) )  &
             +2.*(atan(xeta)-atan(xeta0))
     else            ! jjb: note that the case xmo=0 is not explained in S & P
        print*,'Warning, in SR monin, xmo=0',xmo
     endif
  endif



  print *,'L, Ri',xmo,z/xmo
  print *,'Phi= ',phi

end subroutine monin

!
!----------------------------------------------------------------------
!


subroutine ion_mass (srname)
! calculation of ion mass for ion balance checks

  ! This routine is called by SR profm (in outp.f90)

! Author:
! ------
  !    RvG?


! Modifications :
! -------------
  !

! == End of header =============================================================

  USE config, ONLY : &
! Imported Parameters:
       nkc_l

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunout

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j6, &
       nf, &
       n, &
       nka, &
       nkt, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  character (len=10), intent(in) :: srname
  integer :: ia, j, jt, k
  real (kind=dp) :: xHp, xNHp, xNap
  real (kind=dp) :: xSOm, xHCOm, xNOm, xClm, xHSOm, xCHSO
  real (kind=dp) :: xxsum, xsumi
  real (kind=dp) :: xsum(2:nf)

! Common blocks:
  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw
  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq
  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

! == End of declarations =======================================================

  xHp   = 0._dp
  xNHp  = 0._dp
  xSOm  = 0._dp
  xHCOm = 0._dp
  xNOm  = 0._dp
  xClm  = 0._dp
  xHSOm = 0._dp
  xNap  = 0._dp
  xCHSO = 0._dp
  xxsum = 0._dp

  do k=2,nf
     do j=1,nkc_l
        xHp   = xHp   + sion1(1,j,k) *detw(k)*1._dp  *1.e6_dp
        xNHp  = xNHp  + sion1(2,j,k) *detw(k)*19._dp *1.e6_dp
        xSOm  = xSOm  + sion1(8,j,k) *detw(k)*96._dp *1.e6_dp
        xHCOm = xHCOm + sion1(9,j,k) *detw(k)*61._dp *1.e6_dp
        xNOm  = xNOm  + sion1(13,j,k)*detw(k)*62._dp *1.e6_dp
        xClm  = xClm  + sion1(14,j,k)*detw(k)*35.5_dp*1.e6_dp
        xHSOm = xHSOm + sion1(19,j,k)*detw(k)*97._dp *1.e6_dp
        xNap  = xNap  + sion1(20,j,k)*detw(k)*23._dp *1.e6_dp
        xCHSO = xCHSO + sion1(30,j,k)*detw(k)*95._dp *1.e6_dp
     enddo

     xsum(k)=0._dp
     do ia=1,nka
        do jt=1,nkt
           xsum(k)=xsum(k)+ff(jt,ia,k)*en(ia)
        enddo
     enddo
     xsum(k)=xsum(k)*1.e+09_dp
     xxsum=xxsum+xsum(k)*detw(k)
  enddo

  xsumi=xHp+xNHp+xSOm+xHCOm+xNOm+xClm+xHSOm+xNap+xCHSO

  write (jpfunout,21) srname
  write (jpfunout,22) xxsum,xsumi,xNHp,xSOm,xHCOm,xNOm,xClm,xHSOm,xNap,xCHSO

21 format (a10)
22 format (10d14.6)

end subroutine ion_mass


!
!-------------------------------------------------------
!

subroutine box_init (nlevbox,nz_box,n_bl,BL_box)
!     initialisation for box models runs


! Author:
! ------
  !    RvG?


! Modifications :
! -------------
  !

! == End of header =============================================================

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  integer, intent(in) :: nlevbox, nz_box, n_bl
  logical, intent(in) :: BL_box

  real(kind=dp), external :: p21

  integer :: k
  real(kind=dp) :: t0, tsum, xm10, xmsum

! Common blocks:
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a

!      common /boxdat/ t0, xm10 ! this CB was fed here, but used nowhere else
  common /kinv_i/ kinv
  integer :: kinv

! == End of declarations =======================================================

! initialize kinv (needed in SR kpp_driver)
  kinv=nf
! important for restart only: call init_konc only to start with "fresh" aerosol
!      call init_konc

  if (.not.BL_box) then
! init vals from certain level:
     t0   = t(nlevbox)
     xm10 = xm1(nlevbox)
     print *,' level used for box run: ',nlevbox
  else
!  arithmetic average over box
     tsum  = 0.
     xmsum = 0.
     do k=2,nz_box
        tsum  = tsum + t(k)
        xmsum = xmsum + xm1(k)
     enddo
     t0   = tsum/(nz_box-1)
     xm10 = xmsum/(nz_box-1)
  endif
  t(n_bl)  = t0
  xm1(n_bl)= xm10
  feu(n_bl)=xm1(n_bl)*p(n_bl)/((0.62198+0.37802*xm1(n_bl))*p21(t(n_bl)))
  print *,"box temp and hum: ",t(n_bl),xm1(n_bl),feu(n_bl)

! for strange reasons it didn't work to init the whole Mistra column and
! not producing a crash or strange results, therefore this 1D column
! init is repeated every hour in SR box_update

! if smogchamber run: adjust roughness length z0
!     z0 = to be determined; maybe scale with the decline of the particle population
!     in the chamber; make sure that this value is not overwritten later during the
!     run

! also make sure for smogchamber that deposition occurs also to the sidewalls of the chamber

end subroutine box_init

!
!-------------------------------------------------------
!

!      subroutine box_update (box_switch,ij,nlevbox,nz_box,n_bl,chem, ! jjb unused arguments CHEM, HALO, IOD
subroutine box_update (box_switch,ij,nlevbox,nz_box,n_bl, &
!           halo,iod,BL_box)
          BL_box)

! update of astronomical, aerosol and meteorological properties for box model runs


! Author:
! ------
  !    RvG?


! Modifications :
! -------------
  !

! == End of header =============================================================

  USE constants, ONLY : &
! Imported Parameters:
       pi

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nrlay, &
       nkc, &
       mbs

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  interface
     subroutine equil (ncase,kk)
       integer,           intent(in)  :: ncase
       integer, optional, intent(in)  :: kk    ! model level for calculation
     end subroutine equil
  end interface

!      logical chem,halo,iod,fa_lse,BL_box ! jjb unused arguments removed: chem, halo, iod
  logical fa_lse,BL_box
  real (kind=dp) :: box_switch
  integer :: ij, nlevbox, n_bl, nz_box

  real (kind=dp), external :: p21

  integer :: k
  real (kind=dp) :: xph3, xph4
  real (kind=dp) :: horang, rlat, rdec, ru0, u00, zeit

! Common blocks:
  common /cb16/ u0,albedo(mbs),thk(nrlay)
  real (kind=dp) :: u0, albedo, thk

  common /cb18/ alat,declin                ! for the SZA calculation
  real (kind=dp) :: alat,declin

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a

  common /blck12/ cw(nkc,n),cm(nkc,n)
  real (kind=dp) :: cw, cm

!     dimension rc(nf,nkc),freep(nf) ! jjb rc now in blck11
  real (kind=dp) :: freep(n)

! == End of declarations =======================================================

!      nmin = n_bl ! jjb variable unreferenced
!      nmax = n_bl ! jjb variable unreferenced

  fa_lse = .false. ! jjb was missing, thus undefined when calling SRs fast_k_mt_* below

! calculate u0
  rlat=alat*1.745329e-02
  rdec=declin*1.745329e-02
  zeit=lst*3600.+real(lmin-1,dp)*60.
! greater intervals and variable dtg is known only in the old chemical module
  horang=7.272205e-05*zeit-pi
  u00=cos(rdec)*cos(rlat)*cos(horang)+sin(rdec)*sin(rlat)
  ru0=6371.*u00
  u0=8./(sqrt(ru0**2+102000.)-ru0)

! new photolysis rates are calculated from main program!

! get data for whole column for cw, rc, xkmt for box runs with averaging
! over BL:

! if T, rH, .. vary in time update every hour or at "reasonable times" if not,
! just call this at the beginning of the run

!      if (lmin.eq.1.and.ij.eq.1) then
  if (box_switch.eq.1..and.lmin.eq.1.and.ij.eq.1) then
!    free path length (lambda=freep):
     do k=2,n
        freep(k)=2.28e-5 * t(k) / p(k)
     enddo

     call cw_rc (nf)

     call v_mean_a (t,nf)
     call henry_a (t,nf)
     call fast_k_mt_a(freep,fa_lse,nf)
     call equil_co_a (t,nf)
     call activ (fa_lse,nf)
     call dry_cw_rc (nf)
     call dry_rates_g (t,freep,nf)
     call dry_rates_a (freep,nf)
     xph3=0.
     xph4=0.
     if (cm(3,n_bl).gt.0.) xph3 = 1.
     if (cm(4,n_bl).gt.0.) xph4 = 1.
     if (xph3.eq.1..or.xph4.eq.1.) then
        call v_mean_t (t,nf)
        call henry_t (t,nf)
        call fast_k_mt_t(freep,fa_lse,nf)
        call equil_co_t (t,nf)
     endif

     if (BL_box) then
!           average parameters over depth of BL if BL_box=.true.
        call ave_parms (n_bl,nz_box)
        call ave_aer (n_bl,nz_box)
        if (xph3.eq.1..or.xph4.eq.1.) call ave_tot (n_bl,nz_box)
     else
        call set_box_gas (nlevbox,n_bl)
        call set_box_lev_a (nlevbox,n_bl)
        feu(n_bl)=xm1(n_bl)*p(n_bl)/((0.62198+0.37802*xm1(n_bl))*p21(t(n_bl)))
        call equil (1,n_bl)
        if (xph3.eq.1..or.xph4.eq.1.) call set_box_lev_t (nlevbox,n_bl)
     endif
!         call print_vals (nlevbox,n_bl)
     box_switch = 0.
  endif


! the following parameters are used at the beginning of SR kpp_driver, so they
! have an effect on the chemistry

! prescribe a variation in T, rh as function of u0
!     put SINUS here, t0 xm10
!      t(n_bl)    = xsinus * t0        ! temp [K]
!      xm1(n_bl)  = xsinus * xm10      ! spec humidity [kg/kg]
!      feu(n_bl)=xm1(n_bl)*p(n_bl)/((0.62198+0.37802*xm1(n_bl))*p21(t(n_bl)))

!----------tests from run bug_fix_n: str.f SR surf0 (see also SST.f in bug_fix_n)
!      common /tw_0/ tw0
!c      tw=tw-5.787d-6*dt
!      tw0=tw0-6.94444d-6*dt
!c assume diurnal variation in SST:
!      th      = mod(time/3600.,24.) ! time in [h] (0..24)
!!     pi05    = 0.5 * 3.1425927     ! 1/2 pi
!      pi05    = 0.5 * pi            ! 1/2 pi
!      tmax    = 21.                 ! shift to get maximum of SST at 15:00
!      sst_amp = 1.                  ! amplitude of SST diu var
!      sst_var = sst_amp * sin(pi05*(tmax - th)/6.)
!      tw      = tw0 + sst_var
!      t(1)=tw
!----------
! change only for different LWC
!      conv2(n_bl,1-4)    ! conversion factor [1/(1000 m3/m3)]
!      cm(k,1-4)          ! LWC in particle class [m3/m3]
!      cloud(n_bl,1-4)    ! cloud present in layer [true/false] - should be false everywhere
! don't have to be changed
!      cm3(n_bl,1)        ! air dens. [mlc/cm3]
!      am3(n_bl,1)        ! air dens. [mol/m3]
!      am3(n_bl,2)        ! [CO] [mol/m3]
!      rho(n_bl)          ! density [kg/m3]
!      p(n_bl)            ! pressure [Pa]

end subroutine box_update


!
!---------------------------------------------------------------
!

subroutine sedc_box (dt,z_box,n_bl)
! dry deposition and emission of gaseous species for box runs


! Author:
! ------
  !    Roland von Glasow, duplicate of SR sedc with tuning for box case


! Modifications :
! -------------
  !
! jjb work done = implicit none, missing declarations, little cleaning, modules including constants

! == End of header =============================================================


  USE constants, ONLY : &
! Imported Parameters:
       Avogadro

  USE gas_common, ONLY : &
! Imported Parameters:
       j1, &
! Imported Array Variables with intent (in):
       es1, &
       ind_gas_rev, &
! Imported Array Variables with intent (inout):
       s1, &
       vg

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
  real (kind=dp), intent(in) :: dt, z_box
  integer, intent(in) :: n_bl

! Local scalars:
  integer :: j
  real (kind=dp) :: s12old

! Common blocks
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

! == End of declarations =======================================================



!      vg(4)=0.27e-2 ! NH3 old value, that fitted "nicely" in model    !=0. ! emission is net flux or
!      vg(34)=vg(30) ! N2O5=HCl
!      vg(37)=0.     ! DMS           emission is net flux
!      vg(38)=vg(30) ! HOCl = HCl
!      vg(43)=vg(30) ! HOBr = HCl
!      vg(50)=vg(49) ! I2O2=HOI
!      vg(51)=vg(49) ! INO2=HOI
!      vg(56)=0.     ! CH3I          emission is net flux
!      vg(57)=0.     ! CH2I2         emission is net flux
!      vg(58)=0.     ! CH2ClI        emission is net flux
!      vg(59)=0.     ! C3H7I         emission is net flux
!      vg(63)=vg(30) ! CH3SO3H = HCl
!      vg(71)=0.     ! CH2BrI         emission is net flux
!      vg(72)=0.     ! CHBr2I         emission is net flux
!      vg(73)=0.     ! C2H5I          emission is net flux

  if(ind_gas_rev(4) /= 0) &
      vg(ind_gas_rev(4))=0.27e-2_dp              ! NH3 old value, that fitted "nicely" in model    !=0. ! emission is net flux or
  if(ind_gas_rev(34) /= 0 .and. ind_gas_rev(30) /= 0) &
       vg(ind_gas_rev(34))=vg(ind_gas_rev(30)) ! N2O5=HCl
  if(ind_gas_rev(37) /= 0) &
       vg(ind_gas_rev(37))=0._dp                  ! DMS           emission is net flux
  if(ind_gas_rev(38) /= 0 .and. ind_gas_rev(30) /= 0) &
       vg(ind_gas_rev(38))=vg(ind_gas_rev(30)) ! HOCl = HCl
  if(ind_gas_rev(43) /= 0 .and. ind_gas_rev(30) /= 0) &
       vg(ind_gas_rev(43))=vg(ind_gas_rev(30)) ! HOBr = HCl
  if(ind_gas_rev(50) /= 0 .and. ind_gas_rev(49) /= 0) &
       vg(ind_gas_rev(50))=vg(ind_gas_rev(49)) ! I2O2=HOI
  if(ind_gas_rev(51) /= 0 .and. ind_gas_rev(49) /= 0) &
       vg(ind_gas_rev(51))=vg(ind_gas_rev(49)) ! INO2=HOI
  if(ind_gas_rev(56) /= 0) &
       vg(ind_gas_rev(56))=0._dp                  ! CH3I          emission is net flux
  if(ind_gas_rev(57) /= 0) &
       vg(ind_gas_rev(57))=0._dp                  ! CH2I2         emission is net flux
  if(ind_gas_rev(58) /= 0) &
       vg(ind_gas_rev(58))=0._dp                  ! CH2ClI        emission is net flux
  if(ind_gas_rev(59) /= 0) &
       vg(ind_gas_rev(59))=0._dp                  ! C3H7I         emission is net flux
  if(ind_gas_rev(63) /= 0 .and. ind_gas_rev(30) /= 0) &
       vg(ind_gas_rev(63))=vg(ind_gas_rev(30)) ! CH3SO3H = HCl
  if(ind_gas_rev(71) /= 0) &
       vg(ind_gas_rev(71))=0._dp                  ! CH2BrI         emission is net flux
  if(ind_gas_rev(72) /= 0) &
       vg(ind_gas_rev(72))=0._dp                  ! CHBr2I         emission is net flux
  if(ind_gas_rev(73) /= 0) &
       vg(ind_gas_rev(73))=0._dp                  ! C2H5I          emission is net flux


  if (lst/4*4.eq.lst.and.lmin.eq.1) then
     print *,lday,lst,lmin
     print *,' dry deposition velocities'
     do j=1,j1
        print *,j,vg(j)
     enddo
  endif

  do j=1,j1
! deposition, vg in m/s
     if (vg(j).ge.1.e-5_dp) then
        s12old=s1(j,n_bl)
        s1(j,n_bl)=s1(j,n_bl)*exp(-dt/z_box*vg(j))
        s1(j,1)=s1(j,1)+(s12old-s1(j,n_bl))*z_box
     endif

! emission
! es1: emission rates in molec./cm**2/s, s1 in mol/m**3
     s1(j,n_bl)=s1(j,n_bl)+es1(j)*dt*1.e+4_dp/(z_box*Avogadro)
  enddo

end subroutine sedc_box


!
!----------------------------------------------------------------
!

subroutine box_partdep (dt, z_box, n_bl)


! Author:
! ------
  !    RvG?


! Modifications :
! -------------
  !

! == End of header =============================================================

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j6, &
       nf, &
       n, &
       nka, &
       nkt, &
       nkc

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  real (kind=dp), intent(in) :: dt, z_box
  integer, intent(in) :: n_bl

  integer :: ia, jt, kc, l
  real (kind=dp) :: ff_old, s_old, x_depterm

! dry deposition of particles and aqueous constituents in box
! Common blocks:
  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1
  common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
  real (kind=dp) :: vt, vd, vdm

! == End of declarations =======================================================

! calculation of deposition velocity is done in SR partdep; for smog chamber runs
! the roughness length z0 has to be adjusted in SRs box_init

! apply v_dep to particles
  do ia = 1, nka
     do jt = 1, nkt
!        x = x0 * exp(-dt/z_box*vdep(l))
        ff_old = ff(jt,ia,n_bl)
        ff(jt,ia,n_bl) =  ff_old * exp(-dt/z_box*vd(jt,ia))
! add deposited numbers
        ff(jt,ia,1) = ff(jt,ia,1) + (ff_old - ff(jt,ia,n_bl))*z_box
     enddo
  enddo
! apply v_dep to non-ionic aqueous constituents
  do kc = 1, nkc
     x_depterm = exp(-dt/z_box*vdm(kc))
     do l = 1, j2
        s_old = sl1(l,kc,n_bl)
        sl1(l,kc,n_bl) = s_old * x_depterm
! add deposited numbers
        sl1(l,kc,1) = sl1(l,kc,1) + (s_old - sl1(l,kc,n_bl))*z_box
     enddo
  enddo
! apply v_dep to ionic aqueous constituents
  do kc = 1, nkc
     x_depterm = exp(-dt/z_box*vdm(kc))
     do l = 1, j6
        s_old = sion1(l,kc,n_bl)
        sion1(l,kc,n_bl) = s_old * x_depterm
! add deposited numbers
        sion1(l,kc,1) = sion1(l,kc,1) + (s_old - sion1(l,kc,n_bl))*z_box
     enddo
  enddo

end subroutine box_partdep


!
!---------------------------------------------------------------------
!
!
!----------------------------------------------------------------------
!

subroutine out_mass
! subroutine to print aerosol and ion mass

  USE config, ONLY : &
       coutdir,      &
       nkc_l

  USE global_params, ONLY : &
! Imported Parameters:
       j2, &
       j6, &
       n, &
       nka, &
       nkt, &
       nkc

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunom

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local parameters
  character (len=*), parameter :: clfname = "out_mass.out"
  integer, parameter :: lsp = 9                               ! number of tracked ions
  ! Tracked ions are: H+, NH4+, SO42-, HCO3-, NO3-, Cl-, HSO4-, Na+, CH3SO3-
  integer, parameter :: lj2(lsp) = (/1,2,8,9,13,14,19,20,30/) ! indexes of tracked ions
  real (kind=dp), parameter :: xmm(lsp) = &                   ! molar mass of mass-determining ions in g/mole
       (/1._dp, 18._dp, 96._dp, 61._dp, 62._dp, 35.5_dp, 97._dp, 23._dp, 95._dp/)

! Local scalars
  integer :: ia, jt, k, kc, l, ll
  real (kind=dp) :: xxsum
! Local arrays
  real (kind=dp) :: xsum(n), xionmass(nkc,n), xion(n)

! Common blocks:
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
  real (kind=dp) :: sl1, sion1

! == End of declarations =======================================================

  open (unit=jpfunom, file=trim(coutdir)//clfname, status='unknown', form='formatted', position='append')
  write (jpfunom,*) ' '
  write (jpfunom,*) 'output:', lday,lst,lmin

! Initialisations
  xxsum         = 0._dp
  xsum(1)       = 0._dp
  xionmass(:,:) = 0._dp ! must initialise all nkc, if nkc_l < nkc (but used until nkc anyways)

! output of aerosol mass-------------
  do k=2,n
     xsum(k)=0._dp
     do ia=1,nka
        do jt=1,nkt
           xsum(k) = xsum(k) + ff(jt,ia,k) * en(ia)
        enddo
     enddo
     xsum(k) = xsum(k) * 1.e+09_dp
     xxsum = xxsum + xsum(k) * detw(k)
  enddo

  write (jpfunom,6240)
  write (jpfunom,6250) xsum
  write (jpfunom,6260) xxsum

! output of ion mass-------------
  do k=1,n
     do kc=1,nkc_l
! calculate the mass
        xionmass(kc,k)=0._dp
        do l=1,lsp
           ll = lj2(l)
           xionmass(kc,k) = xionmass(kc,k) + sion1(ll,kc,k) * xmm(l)
        enddo
     enddo
     xion(k) = (xionmass(1,k) + xionmass(2,k) + xionmass(3,k) + xionmass(4,k)) * 1.e6_dp
  enddo

  write (jpfunom,6280)
  write (jpfunom,6250) xion

 6240 format (/,6x,'aerosol mass in ug m**-3 in layers 2 - nf')
 6250 format (1x,15f8.3)
 6260 format (6x,'total aerosol mass in ug m**-2 of layers 2 - nf',f12.3)
 6280 format (/,6x,'ion mass in ug m**-3 in layers 2 - nf')

  close (jpfunom)

end subroutine out_mass

!
!-----------------------------------------------------------------------------
!

subroutine get_n_box (z_box,nz_box)
! get grid level for box height

  USE global_params, ONLY : &
! Imported Parameters:
       n

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  real (kind=dp), intent(in) :: z_box
  integer, intent(out) :: nz_box
  integer :: k
! Common blocks:
  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

! == End of declarations =======================================================


  nz_box=0
  do k=1,n
     if (etw(k).ge.z_box) then
        nz_box=k
        exit
     endif
  enddo

  if (nz_box.eq.0) then
     nz_box=70
     print *,"WARNING!! box height too high, set to ",etw(nz_box)
  endif

  print *,"box height set to ",etw(nz_box)
  print *,"box height corresponds to level ",nz_box

end subroutine get_n_box

!
!-------------------------------------------------------------------------
!

subroutine ffssuumm
! calculation of particle number

!      jjb cleaning
!          do end do without label

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nka, &
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  integer :: ia, jt, k
  real (kind=dp) :: fsum1, fsum2
! Common blocks:
  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /blck06/ kw(nka),ka
  integer :: kw, ka

! == End of declarations =======================================================


  do k=2,n
     fsum1=0.
     fsum2=0.
! small aerosol
     do ia=1,ka
        do jt=1,kw(ia)
           fsum1=fsum1+ff(jt,ia,k)
        enddo
     enddo
! large aerosol
     do ia=ka+1,nka
        do jt=1,kw(ia)
           fsum2=fsum2+ff(jt,ia,k)
        enddo
     enddo
     print *,k,fsum1,fsum2,fsum(k)
  end do

end subroutine ffssuumm

!
!----------------------------------------------------------------
!

subroutine oneD_dist
!  calculate 1D size distribution of 2D particles dist.

  USE constants, ONLY : &
! Imported Parameters:
       pi

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nka, &
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  integer :: ia, ij, jt, jtt, k
  real (kind=dp) :: fpi

  real (kind=dp) :: rp(nka+nkt)  !particle radius [um]
  real (kind=dp) :: Np(nka+nkt)  !particle number [part cm-3]
  real (kind=dp) :: Ap(nka+nkt)  !particle surface [um2 cm-3]
  real (kind=dp) :: Vp(nka+nkt)  !particle volume [um3 cm-3]

! Common blocks:
  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /oneDs_0/partN(n,nka+nkt),partA(n,nka+nkt),partV(n,nka+nkt), &
       partr(n,nka+nkt),drp(nka+nkt),nrp
  real (kind=dp) :: partN, partA, partV, partr, drp
  integer :: nrp

! == End of declarations =======================================================

!     set up radius range problem: if rq(nkt,1) is used to map all other
!     radii on, every now and then 2 rq(jt,i) bins will fall into the
!     same rq(1:nkt,1) bin and produce artifial spikes in the size
!     distribution. The reason is that the radius resolution for small
!     ia is very fine (see plot of 2D particle grid). To avoid this use
!     rq(1,ia) until ia=nka and then rq(1:nkt,nka); the dimension of rq
!     is nka+nkt; to avoid having drp=0. the


  do ia=1,nka
     rp(ia)=rq(1,ia)
  enddo
  do jt=1,nkt
     jtt=jt
     if (rq(jt,nka).gt.rp(nka)*1.05) goto 1021
  enddo
1021 continue
  do jt=1,nkt-jtt
     rp(nka+jt)=rq(jt+jtt,nka)
  enddo
  do ij=1,nka+nkt-1
     drp(ij)=rp(ij+1)-rp(ij)
!     print *,ij,rp(ij),drp(ij)
  enddo
  nrp=nka+nkt-jtt
  drp(nrp)=drp(nrp-1)

  do k=2,n
     do ij = 1,nrp
        Np(ij) = 0.
        Ap(ij) = 0.
        Vp(ij) = 0.
        fpi=4.*pi
        do ia = 1, nka
           if (rn(ia).gt.rp(ij)) goto 2001
           do jt = 1,nkt
              if (rq(jt,ia).le.rp(ij)) then
                 if (ij.gt.1) then
                    if (rq(jt,ia).gt.rp(ij-1)) then
                       Np(ij) = Np(ij) + ff(jt,ia,k)
                       Ap(ij) = Ap(ij) + ff(jt,ia,k)*fpi*rp(ij)*rp(ij)
                       Vp(ij) = Vp(ij) + ff(jt,ia,k)*fpi*rp(ij)*rp(ij)*rp(ij)/3.
                    endif
                 else if ((ij.eq.1).and.(rq(jt+1,ia).gt.rp(ij))) then           ! jjb BUG here, jt+1 leads to out of bounds index when jt = nkt
!                write (*,100) ij,ia,jt,rq(jt,1),rp(ij),rq(jt,ia),ff(jt,ia,k)
                    Np(ij) = Np(ij) + ff(jt,ia,k)
                    Ap(ij) = Ap(ij) + ff(jt,ia,k)*fpi*rp(ij)*rp(ij)
                    Vp(ij) = Vp(ij) + ff(jt,ia,k)*fpi*rp(ij)*rp(ij)*rp(ij)/3.
                 endif
              else
                 goto 2002
              endif
           enddo                 !jt
2002       continue
        enddo                    !ia
2001    continue

! to plot dN/dr, dA/dr, dV/dr:
        partN(k,ij) = Np(ij)/drp(ij)
        partA(k,ij) = Ap(ij)/drp(ij)
        partV(k,ij) = Vp(ij)/drp(ij)
        partr(k,ij) = rp(ij)

     enddo                     !ij
  enddo                     ! k

! 100  format(3i4,4d16.8)

end subroutine oneD_dist


!
!----------------------------------------------------------------
!

subroutine oneD_dist_new
!  calculate 1D size distribution of 2D particles dist.

! jjb rewritten, but needs improvements.
!     idea: map all rq onto a specific which would range between rq(1,1) and rq(nkt,nka) with XXX values
!     (XXX could be equal to any of the already used parameters, nka, nkt, nka+nkt for finer resolution, or whatever)
!     currently, many particles are likely to be in the last rp class

  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nka, &
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  integer :: ia, ij, jt, k
  real (kind=dp) :: rmin, rmax, rfact, xnsum
  real (kind=dp) :: Np(nkt)  !particle number [part cm-3]
  real (kind=dp) :: rp(0:nkt)  !particle radius [um]
  real (kind=dp) :: xlogdrp(nkt) ! jjb removed from /oneDs/

! Common blocks:
  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  !     common /oneDs/  partN(n,nkt,2),partr(n,nkt),drp(nkt),xlogdrp(nkt) ! jjb last variable shouldn't be in this CB
  ! jjb : reindexing would be good, and decrease n to nf of increase the loop below (and in out_netcdf)
  common /oneDs/  partN(n,nkt,2),partr(n,nkt),drp(nkt)              ! jjb removed and declared below
  real (kind=dp) :: partN, partr, drp

! == End of declarations =======================================================

! Ap, Vp can easily be calculated in ferret, so reduce output file size

! ignore radius range problem as outlined in SR oneD_dist_old and use rq(1:nkt,1)
! to map all particles onto; this might lead to adding a few bins of 2D spectrum
! into one on 1D but I think this is just a "cosmetic" problem which in the end
! might be overcome by smoothing the plot. Benefit: everything is easier and there
! are no "empty" bins as when using the old radius 2D --> 1D mapping nka+nkt

! define 1D grid
  rmin = rw(1,1)
  rmax = rw(nkt,nka)
  rfact = 10**(log10(rmax/rmin)/(nkt-1))
  rp(0) = rmin/rfact
  do jt=1,nkt
     rp(jt)=rp(jt-1)*rfact
  enddo

! calculate width of each bin in 1D grid
!      drp(1)=rp(1)
!      xlogdrp(1)=log10(rp(1))
  do ij=1,nkt
     drp(ij)=rp(ij)-rp(ij-1)
! unit: "implicit" division of rp by 1um to get a unit-less property to be able to use log
     xlogdrp(ij)=log10(rp(ij))-log10(rp(ij-1))
!         print *,ij,rp(ij),drp(ij)
!         print *,drp(ij),xlogdrp(ij),log10(drp(ij))
  enddo

! map 2D spectrum on 1D spectrum
  do k=2,nf
     Np(:)=0.
     do ij = 1,nkt-1
        do ia = 1, nka
           if (rn(ia).gt.rp(ij+1)) goto 2001 ! save time
           do jt = 1,nkt
              if (rq(jt,ia).lt.rp(ij+1)) then
                 if (rq(jt,ia).ge.rp(ij)) then
                    Np(ij) = Np(ij) + ff(jt,ia,k)
                 endif
              else
                 if (ij == nkt-1) then ! special case
                    Np(nkt) = Np(nkt) + ff(jt,ia,k)
                 else
                    goto 2002 ! save time
                 end if
              endif
           enddo           !jt
2002       continue
        enddo              !ia
2001    continue
     end do                !ij
     Np(2) = Np(2)+ff(1,1,k) ! jjb correction to be tested
     do ij = 1,nkt
! to plot dN/dr, dA/dr, dV/dr:
        partN(k,ij,1) = Np(ij)/drp(ij)      ! #/um/cm3
! to plot dN/dlogr, dA/dlogr, dV/dlogr:
        partN(k,ij,2) = Np(ij)/xlogdrp(ij)  ! #/cm3
        partr(k,ij) = rp(ij)
!         print *,partN(k,ij),partA(k,ij),partV(k,ij),partr(k,ij)
     enddo                    !ij

! check if bins were "missed":
     xnsum=0.
     do jt=1,nkt
        xnsum=xnsum+Np(jt)
     enddo
     if (xnsum.lt..99*fsum(k)) print *,'N too small',k,fsum(k),xnsum
     if (xnsum.gt.1.01*fsum(k)) print *,'N too big',k,fsum(k),xnsum


  enddo                     ! k

! 100  format(3i4,4d16.8)

end subroutine oneD_dist_new


!
!----------------------------------------------------------------
!

subroutine oneD_dist_jjb
!  calculate 1D size distribution of 2D particles dist.


  USE global_params, ONLY : &
! Imported Parameters:
       nf, &
       n, &
       nka, &
       nkt

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  integer :: ia, ij, jt, k
  real (kind=dp) :: xnsum
  real (kind=dp) :: ff1D(nka-1)

! Common blocks:
  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  double precision enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff, fsum
  integer :: nar

  common /oneDsj/ rpw(nka), part1D(nka-1,nf)
  real (kind=dp) :: rpw, part1D

! == End of declarations =======================================================


! map 2D spectrum on 1D spectrum
  do k=2,nf
     ff1D(:) = 0._dp
     do ia = 1, nka
        ij=1 ! initialisation within ia loop is not a mistake:
             !  for a given ia, rq is strictly increasing when jt increases
             !  thus once the first ij has been updated for a given (ia,1), the next ij
             !  will be at least equal to the previous ij
        do jt = 1,nkt
           do while (rq(jt,ia) > rpw(ij+1))
              ij = ij+1
           end do
           ff1D(ij) = ff1D(ij) + ff(jt,ia,k)
        enddo           !jt
     enddo              !ia

     part1D(:,k) = ff1D(:)

! check if bins were "missed":
     xnsum = sum(ff1D(:))
     if (xnsum.lt..99*fsum(k)) print *,'N too small',k,fsum(k),xnsum
     if (xnsum.gt.1.01*fsum(k)) print *,'N too big',k,fsum(k),xnsum

  enddo                  ! k

end subroutine oneD_dist_jjb


!
!--------------------------------------------------------------------------
!

! Compute the latent heat of vaporisation as a function of temperature
function xl21(temperature)


  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Function arguments
  ! Scalar arguments with intent(in):
  real(kind=dp), intent(in) :: temperature         ! in [K]

! Local parameters:
  real(kind=dp)             :: ppA = 3138708._dp   ! in J/kg
  real(kind=dp)             :: ppB = -2339.4_dp    ! in J/kg/K
! Local scalars:
  real(kind=dp)             :: xl21                ! in J/kg

! == End of declarations =======================================================

  xl21 = ppA + ppB*temperature

end function xl21


!
!--------------------------------------------------------------------------
!

! Compute the saturation water vapour pressure as a function of temperature
!  after Magnus formula
function p21(ttt)


  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Function arguments
  ! Scalar arguments with intent(in):
  real(kind=dp), intent(in) :: ttt         ! temperature, in [K]

! Local scalars:
  real(kind=dp)             :: p21         ! saturation vapour pressure, in [Pa]

! == End of declarations =======================================================


  p21 = 610.7_dp * exp(17.15_dp*(ttt-273.15_dp)/(ttt-38.33_dp))

end function p21
!
!--------------------------------------------------------------------------
!

subroutine chamb_init (n_bl)
! initialisation for chamber models runs

  USE config, ONLY : &
! Imported Parameters:
       cinpdir_phot

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunout, jpfunJchamb

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nf, &
       nphrxn

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  interface
     subroutine equil (ncase,kk)
       integer,           intent(in)  :: ncase
       integer, optional, intent(in)  :: kk    ! model level for calculation
     end subroutine equil
  end interface

  integer, intent(in) :: n_bl
  character (len=110) :: clpath
  integer :: k
  logical :: lights
  real (kind=dp) :: zp21
  real (kind=dp) :: freep(n)
  real (kind=dp), external :: p21

! Common blocks:
  common /boxdat/ t0,rh0
  real(kind=dp) :: t0,rh0
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n)
  real(kind=dp) :: xm1, xm2, feu, dfddt, xm1a
  common /chamber_ph_r/ ph_rat_chamber(nphrxn)
  real (kind=dp) :: ph_rat_chamber
  common /kinv_i/ kinv
  integer :: kinv


! initialize kinv (needed in SR kpp_driver)
  kinv=nf
! important for restart only: call init_konc only to start with "fresh" aerosol
!      call init_konc

! get chamber initial conditions from file
  clpath=trim(cinpdir_phot)//'chamber.dat'
  open (jpfunJchamb,file=trim(clpath),status='old')
  read (jpfunJchamb,5101) t0
  read (jpfunJchamb,5101) rh0
  close(jpfunJchamb)
 5101 format (f6.0)

  t(n_bl)  = t0
  feu(n_bl)= rh0*1.0d-2
  zp21 = p21(t(n_bl))
  xm1(n_bl)=(0.62198*feu(n_bl)*zp21)/(p(n_bl)-0.37802*feu(n_bl)*zp21)

  write(jpfunout,*)"box temp and hum: ",t(n_bl),xm1(n_bl),feu(n_bl)

! photolysis rates
! start chamber run with lights off (typical experiment, otherwise set lights=TRUE)
! model run start at midday as set in SR initm
  lights = .false.
  if (lights) then
     call photol_chamber
  else
     ph_rat_chamber(:) = 0.0_dp
  endif

!! comment this to turn off aerosol chemistry in the chamber
! free path length (lambda=freep):
  do k=2,n
     freep(k)=2.28e-5 * t(k) / p(k)
  enddo

  call cw_rc (nf)

  call v_mean_a  (t,nf)
  call henry_a (t,nf)
  call fast_k_mt_a(freep,.false.,nf)
  call equil_co_a (t,nf)
  call activ (.false.,nf)
  call dry_cw_rc (nf)
  call dry_rates_g (t,freep,n)
  call dry_rates_a (freep,nf)

  call equil (1,n_bl)

!     if smogchamber run: adjust roughness length z0 z0 = to be
!     determined; maybe scale with the decline of the particle
!     population in the chamber; make sure that this value is not
!     overwritten later during the run

!     also make sure for smogchamber that deposition occurs also to the
!     sidewalls of the chamber


end subroutine chamb_init

!
!-------------------------------------------------------
!

subroutine chamb_update (n_bl,ij)
! turn on and off photolysis rates. if lights are not OFF at the beginning of
! the experiment change SR chamb_init

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nf, &
       nphrxn

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  interface
     subroutine equil (ncase,kk)
       integer,           intent(in)  :: ncase
       integer, optional, intent(in)  :: kk    ! model level for calculation
     end subroutine equil
  end interface

  integer, intent(in) :: ij, n_bl
  integer :: k
  logical :: lights
  real (kind=dp) :: on_hr, on_min, off_hr, off_min, on_time, off_time
  real (kind=dp) :: freep(n)
  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct
  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real(kind=dp) :: theta, thetl, t, talt, p, rho
  common /chamber_ph_r/ ph_rat_chamber(nphrxn)
  real (kind=dp) :: ph_rat_chamber



! model time in chamber mode starts at 12:00 (set in SR initm)
! time since beginning of experiment when lights are turned ON
  on_hr = 0._dp
  on_min = 15._dp
  on_time = (on_hr*3600.)+(on_min*60.)

! time since beginning of experiment when lights are turned OFF
  off_hr = 2._dp
  off_min = 0._dp
  off_time = (off_hr*3600.)+(off_min*60.)

  lights=.false.
  if (time.ge.on_time) then
     lights=.true.
  endif
  if (time.ge.off_time) then
     lights=.false.
  endif

! if lights are ON calculate photolysis rates, otherwise set photolysis rates to zero
  if (lights) then
     call photol_chamber
  else
     ph_rat_chamber(:) = 0.0_dp
  endif

!! comment this to turn off aerosol chemistry in the chamber
  if (lmin.eq.1.and.ij.eq.1) then
!    free path length (lambda=freep):
     do k=2,n
        freep(k)=2.28e-5 * t(k) / p(k)
     enddo

     call cw_rc (nf)

     call v_mean_a  (t,nf)
     call henry_a (t,nf)
     call fast_k_mt_a(freep,.false.,nf)
     call equil_co_a (t,nf)
     call activ (.false.,nf)
     call dry_cw_rc (n)
     call dry_rates_g (t,freep,n)
     call dry_rates_a (freep,nf)

     call equil (1,n_bl)
  endif

end subroutine chamb_update
