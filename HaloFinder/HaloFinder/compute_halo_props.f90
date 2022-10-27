module compute_halo_props

  use input_output
  use halo_defs
#ifdef STARS
  use utils
#endif

  public

contains

!//////////////////////////////////////////////////////////////////////////
!**************************************************************************
  subroutine init()

    implicit none

    write_resim_masses = .false.

    ! initialize gravitationalsoftening
    call initgsoft()

    ! initialize cosmological and technical parameters of the simulation
    call init_cosmo()

    return

  end subroutine init

!*************************************************************************
  subroutine initgsoft

  ! subroutine to initialize the required arrays for the gravitational
  ! field smoothing interpolation in routine SOFTGRAV. This is only
  ! required for a cubic spline kernel; interpolation is performed in
  ! distance.

    implicit none

    integer(kind=4) :: i
    real(kind=8)    :: deldrg,xw,xw2,xw3,xw4,tiny,one,two

    tiny = 1.d-19
    one  = 1.0
    two  = 2.0

    if (gravsoft .eq. 'harmonic') return

    deldrg      = 2./ninterp
    phsmooth(0) = 7.*sqrt(tiny)/5.

    do i=1,1+ninterp
       xw  = i*deldrg
       xw2 = xw*xw
       xw3 = xw2*xw
       xw4 = xw2*xw2
       if (xw .le. one) then
          phsmooth(i) = -2.*xw3*(one/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.
       else
          phsmooth(i) = -one/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-xw3/30.)
       endif
       if (xw .ge. two) then
          phsmooth(i) = one
       endif
    enddo

    return

  end subroutine initgsoft

!*************************************************************************
  subroutine init_cosmo()

    ! This routine reads in the input_HaloMaker.dat file which contains the cosmological
  ! and technical parameters of the N-Body simulation to analyze.

    use neiKDtree
    use fof
    implicit none

    integer(kind=4)      :: i
    ! cannot do otherwise than setting all the strings to a value larger than that of a line
    ! in the input file and trim them whenever it is needed
    character(len=200)   :: line,name,value


    ! set default parameter values
    ! ----------------------------
    ai             = 1.0         ! Initial (beginning of the simulation) expansion factor
    omega_f        = 0.3333      !   mass density at final timestep
    omega_lambda_f = 0.6667      !   lambda at final timestep
    af             = 36.587      !   expansion factor of the final timestep
    Lf             = 150.0       !   final length of box in physical Mpc
    H_f            = 66.667      !   Hubble constant in km/s/Mpc, at final timestep
    b_init         = 0.2         !   linking length friend-of-friend parameter @ z=0
    Nmembers       = 20          !   minimum number of particles for a halo
    nsteps         = 1           !   number of files to analyse (listed in 'files.dat')
    method         = "FOF"
#ifdef ANG_MOM_OF_R
    agor_file      = "ang_mom_of_r"
#endif

    write(errunit,*) '> Values of input parameters:  '
    write(errunit,*) '> ---------------------------- '
    write(errunit,*) ''

    write(errunit,*) '> Looking for input_HaloMaker.dat in directory: ',trim(data_dir)
    call open_exist(20,'input_HaloMaker.dat')
    do
       read (20,'(a)',end=2) line
       i = scan(line,'=')
       if (i == 0 .or. line(1:1) == '#') cycle
       name  = trim(adjustl(line(:i-1)))
       value = trim(adjustl(line(i+1:)))
       ! check for a comment at end of line !
       i     = scan(value,'!')
       if (i /= 0) value = trim(adjustl(value(:i-1)))
       if(verbose) write(errunit,'(1x,a1,a15,a3,a10)') '>',trim(name),' : ',trim(value)
       select case (trim(name))
       case ('omega_0' , 'Omega_0' , 'omega_f')
          read(value,*) omega_f
       case ('omega_l' , 'lambda_0' , 'lambda_f')
          read(value,*) omega_lambda_f
       case ('af' , 'afinal' , 'a_f')
          read(value,*) af
       case ('Lf' , 'lf' , 'lbox')
          read(value,*) Lf
       case ('H_f', 'H_0', 'H')
          read(value,*) H_f
       case('FlagPeriod')
          read(value,*) FlagPeriod
       case ('n', 'N', 'npart')
          read(value,*) nMembers
       case('cdm')
          read(value,*) cdm
       case ('method' )
          write(method,'(a3)') trim(value)
       case ('b')
          read(value,*) b_init
       case ('nvoisins')
          read(value,*) nvoisins
       case ('nhop')
          read(value,*) nhop
       case ('rhot')
          read(value,*) rho_threshold
       case ('fudge')
          read(value,*) fudge
       case ('fudgepsilon')
          read(value,*) fudgepsilon
       case ('alphap')
          read(value,*) alphap
       case ('verbose')
          read(value,*) verbose
       case ('megaverbose')
          read(value,*) megaverbose
       case ('DPMMC')
          read(value,*) DPMMC
       case ('SC')
          read(value,*) SC
       case ('dcell_min')
          read(value,*) dcell_min
       case ('eps_SC')
          read(value,*) eps_SC
       case ('MS')
          read(value,*) MS
       case ('dms_min')
          read(value,*) dms_min
       case('nsteps','nsteps_do')
          read(value,*) nsteps
       case('dump_stars')
          read(value,*) dump_stars
       case('proper_time')
          read(value,*) proper_time
       ! TODO: case dump_dm, or rename to dump_parts, etc
       case('nchem')
          read(value,*) nchem
#ifdef ANG_MOM_OF_R
       case ('agor_file')
          write(agor_file,'(a,a)') trim(data_dir),trim(value)
#endif
       case default
          write(errunit,*) 'dont recognise parameter: ',trim(name)
       end select
    end do
2   close (20)

    ! Some consistency check
    if (DPMMC .and. SC) then
       write(errunit,*) '> DPMMC and SC both activated: please select only one.'
       stop
    end if
    ! Some consistency check
    if (DPMMC .and. cdm) then
       write(errunit,*) '> DPMMC does not work with CDM def for the centre'
       stop
    end if

    ! initial size of the box in physical Mpc (NB: f index stands for final quantities)
    Lboxp = Lf*(ai/af)

    if (((omega_f+omega_lambda_f) /= 1.0) .and. (omega_lambda_f /= 0.0)) then
       write(errunit,*) '> lambda + non flat Universe not implemented yet'
       stop
    endif

    ! In the most general of cases:
    !     af/aexp         = 1+z
    !     omega(z)        = (H_f/H(z))**2 * (1+z)**3  * omega_f
    !     omega_lambda(z) = omega_lambda_f * (H_f/H(z))**2
    !     omega_c(z)      = omega_c_f * (1+z)**2 * (H_f/H(z))**2
    !     H(z)**2         = H_f**2*( omega_f*(1+z)**3 + omega_c_f*(1+z)**2 + omega_lambda_f)

    omega_c_f = 1. - omega_f - omega_lambda_f
    H_i       = H_f * sqrt( omega_f*(af/ai)**3 + omega_c_f*(af/ai)**2 + omega_lambda_f)
    ! rho_crit = 2.78782 h^2  (critical density in units of 10**11 M_sol/Mpc^3)
    mboxp     = 2.78782*(Lf**3)*(H_f/100.)**2*omega_f

    write(errunit,*)
    write(errunit,*) '> Initial/Final values of parameters:  '
    write(errunit,*) '> -----------------------------------  '
    write(errunit,*)
    write(errunit,*) '> redshift                         : ',af/ai-1.
    write(errunit,*) '> box size (Mpc)                   : ',Lboxp
    write(errunit,*) '> Hubble parameter  (km/s/Mpc)     : ',H_i
    write(errunit,*) '> box mass (10^11 Msol)            : ',mboxp
    write(errunit,*)

    return

  end subroutine init_cosmo

!*************************************************************************
  subroutine new_step()

  ! This is the main subroutine: it builds halos from the particle simulation and
  ! computes their properties ...

    implicit none

    integer(i8b)     :: indexp,i
    integer(kind=4)  :: ierr
    integer(kind=4)  :: found,n_halo_contam,n_subs_contam
    real(kind=8)                         :: read_time_ini,read_time_end
    real(kind=8)                         :: t0,t1
    logical          :: printdatacheckhalo !put to true if bug after make_linked_list
    integer(kind=4)  :: ih
    ! leo
    real(kind=8)     :: fhalo
    ! end leo
#ifdef ANG_MOM_OF_R
    character(200)   :: filename
#endif
#ifdef Test_FOF
    character(200)   :: filelisteparts
#endif
    integer(i8b)     :: npartcontam

    write(errunit,'(1x,a16,i5)') '> Timestep  --->',numero_step
    write(errunit,*) '> -------------------'
    write(errunit,*) ''
#ifdef STARS
    write(errunit,*) 'npart (new_step)=',npart
#endif

    call cpu_time(read_time_ini)

    ! read N-body info for this new step
    call read_data

    ! determine the age of the universe and current values of cosmological parameters
    call det_age()

    ! first compute the virial overdensity (with respect to the average density of the universe)
    ! predicted by the top-hat collapse model at this redshift
    call virial()

    ! Leo:
    ! use delta crit from Bryan&Norman98 fit to set rho_threshold
#ifdef BN98
    write(errunit,*) '> using B&N98 fit to set rho_threshold...'
    call fit_deltac
    fhalo = 1./3.    ! rho_halo_edge = fhalo * rho_halo_mean, fhalo ~ 1/3 assuming an isothermal profile
    rho_threshold = fhalo * delta_crit_z / omega_z
    ! also, overwrite rho_mean and vir_overdens by rho_crit_z and delta_crit_z
    ! used in det_vir_props to correct Mvir such that  mvir = vir_overdens*rho_mean*4./3.*pi*avir*bvir*cvir
    ! if BN98 then mvir = delta_crit_z*rho_crit_z*4./3.*pi*avir*bvir*cvir
    rho_mean = rho_crit_z
    vir_overdens = delta_crit_z
    write(errunit,*) '> rho_threshold, rho_crit, delta_crit :',rho_threshold,rho_crit_z,delta_crit_z
#endif
    write(errunit,*) '> rho_threshold               : ',rho_threshold
    write(errunit,*) '> rho_mean                    : ',rho_mean
    write(errunit,*) '> vir_overdens                : ',vir_overdens
    write(errunit,*)
    ! end leo

    if(method.ne."FOF") then
       allocate(liste_parts(1:nbodies),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) 'Cannot allocate liste_parts'
       endif
    end if
    call make_halos()

#ifdef Test_FOF
    write(filelisteparts,'(a,i3.3)') 'liste_parts_',numstep
    open(unit=35,form="unformatted",status="unknown",file=filelisteparts)
    write(35) nbodies,nb_of_halos
    write(35) liste_parts(1:nbodies)
    close(35)
    ! reset nb of halos not to construct halos
    nb_of_halos = 0
#endif
    ! if there are no halos go to the next timestep
    if (nb_of_halos == 0) then
       write(errunit,*) 'no halos deallocating'
       deallocate(pos,vel)
       if(allocated(density)) deallocate(density)
       if(allocated(mass)) deallocate(mass)
#ifdef STARS
       if(allocated(age_st))deallocate(age_st)
       if(allocated(met_st))deallocate(met_st)
       if(allocated(pf_st))deallocate(pf_st) ! RS
       if(allocated(pz_st))deallocate(pz_st) ! RS
       if(allocated(chem_st))deallocate(chem_st)
       call ct_clear_cosmo()
#endif
       if(allocated(family)) deallocate(family)
       deallocate(liste_parts)
       return
    end if

    allocate(first_part(0:(nb_of_halos+nb_of_subhalos)),stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) 'Cannot allocate first_part'
    endif

    allocate(nb_of_parts(0:(nb_of_halos+nb_of_subhalos)),stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) 'Cannot allocate nb_of_parts'
    endif

    ! make a linked list of the particles so that each member of a halo points to the next
    ! until there are no more members (last particles points to -1)
    allocate(linked_list(0:nbodies+1), stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) 'Cannot allocate linked_list'
    endif
    call make_linked_list()

    ! deallocate liste_parts bc it has been replaced by linked_list
    deallocate(liste_parts)

    ! in case of resimulation (or individual particle masses) count how many halos are contaminated
    ! (i.e. contain "low res" particles)
    ! NB: for numerical reasons, it can happen that mass is sligthly larger
    !     than massp. Thus, for detecting LR particles, the criterion
    !     if (mass(indexp) > massp) ...  may not work (it can get ALL parts, also HR parts!).
    !     Thus, to be conservative, use: if (mass(indexp) > massp * (1+1e-5))) ...
    if (allocated(mass)) then
       n_halo_contam = 0
       n_subs_contam = 0
       do i = 1,nb_of_halos+nb_of_subhalos
          indexp = first_part(i)
          found  = 0
          npartcontam = 0
          do while (indexp /= -1)
             if ((mass(indexp) > massp* 1.00001) .and. (family(indexp) == FAM_DM) ) then
                if (found == 0) then
                   if(fsub) then
                      if (level(i) == 1) then
                         n_halo_contam = n_halo_contam + 1
                      else
                         n_subs_contam = n_subs_contam + 1
                      endif
                   else
                      n_halo_contam = n_halo_contam + 1
                   end if
                end if
                found = 1
                npartcontam=npartcontam+1
             endif
             indexp = linked_list(indexp)
          enddo
          if(found==1)write(333,'(3(i8,2x))')i,npartcontam,nb_of_parts(i)
       enddo
       write(errunit,*) '> # of halos, # of CONTAMINATED halos :',nb_of_halos,n_halo_contam
       write(errunit,*) '> # of subhalos, # of CONTAMINATED subhalos :',nb_of_subhalos,n_subs_contam
       !open(222,file='ncontam_halos.dat',status='unknown',position='append')
       !write(222,'(5(i6,2x))') numero_step,nb_of_halos,n_halo_contam,nb_of_subhalos,n_subs_contam
       !close(222)
       close(333)
       !YDstuff
    endif

    ! allocation and initialization of the halo list
    allocate(liste_halos(0:(nb_of_halos+nb_of_subhalos)),stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) 'Cannot allocate liste_halos'
       stop
    endif
    call init_halos

! TODO compute contam fraction
#ifdef CONTAM
    ! flag contaminated halos or sub-halos
    do i = 1,nb_of_halos+nb_of_subhalos
       liste_halos(i)%contaminated = 0
       indexp = first_part(i)
       found  = 0
       do while (indexp /= -1 .and. found == 0)
          if ((mass(indexp) > massp* 1.00001) .and. (family(indexp) == FAM_DM)) then
             liste_halos(i)%contaminated = 1
             found = 1
          end if
          indexp = linked_list(indexp)
       end do
    end do
    ! contaminate main halos which contain a contaminated halo ...
    do i = 1, nb_of_halos
       ih = i
       do while (ih > 0)
          if (liste_halos(ih)%contaminated ==1) liste_halos(i)%contaminated = 1
          ih = liste_halos(ih)%nextsub
       end do
    end do
    ! contaminate sub-halos which belong to a contaminated main halo ...
    do i = 1,nb_of_halos
       if (liste_halos(i)%contaminated==1) then
          ih = i
          do while (ih > 0)
             liste_halos(ih)%contaminated = 1
             ih = liste_halos(ih)%nextsub
          end do
       end if
    end do
#endif

    ! until now we were using code units for positions and velocities
    ! this routine changes that to physical (non comoving) coordinates for positions
    ! and peculiar (no Hubble flow) velocities in km/s
    ! The masses are changed from code units into 10^11 M_sun as well.
    call change_units()


#ifdef ANG_MOM_OF_R
    write(filename,'(a,a,i3.3)') trim(agor_file),'.',numstep
    open(unit=agor_unit,file=filename,status='unknown',form='unformatted')
    write(agor_unit) nb_of_halos,nb_of_subhalos
    write(agor_unit) nshells
#endif

    printdatacheckhalo = .false.
    write(errunit,*) '> Computing halo properties'

    ! leo: openmp: use of SCHEDULE to distribute packets ?
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) &
    !$OMP PRIVATE(i)
    do i = 1,nb_of_halos + nb_of_subhalos
       if(printdatacheckhalo) then
          write(errunit,*) '> halo:', i,'nb_of_parts',nb_of_parts(i)
          call cpu_time(t0)
       end if
       ! determine mass of halo
       call det_mass(liste_halos(i))
       if(printdatacheckhalo) write(errunit,*) '> mass:', liste_halos(i)%m
       ! compute center of halo there as position of "most dense" particle
       ! and give it the velocity of the true center of mass of halo
       call det_center(liste_halos(i))
       if(printdatacheckhalo) write(errunit,*) '> center:',liste_halos(i)%p
       ! compute angular momentum of halos
       call compute_ang_mom(liste_halos(i))
       if(printdatacheckhalo) write(errunit,*) '> angular momentum:',liste_halos(i)%L
       ! compute r = max(distance of halo parts to center of halo)
       call r_halos(liste_halos(i))
       if(printdatacheckhalo) write(errunit,*) '> radius:',liste_halos(i)%r
       ! compute energies and virial properties of the halos depending on density profile
       ! (so this profile is also computed in the routine)
       call det_vir(liste_halos(i))
       if(printdatacheckhalo) write(errunit,*) '> mvir,rvir:',liste_halos(i)%datas%mvir,liste_halos(i)%datas%mvir
       ! compute dimensionless spin parameter of halos
       call compute_spin_parameter(liste_halos(i))
       if(printdatacheckhalo) write(errunit,*) '> spin:',liste_halos(i)%spin
       ! Compute 3D profiles
       call compute_3d_profile(liste_halos(i))
#ifdef STARS
!!$       ! compute Bulge properties of galaxies
!!$       call compute_bulge_prop(liste_halos(i))
!!$       if(printdatacheckhalo) write(errunit,*) '> m_bulge:',liste_halos(i)%m_bulge
       ! compute stellar surface density profiles
       call compute_stellar_profile(liste_halos(i))
       ! compute mass-weighted stellar metallicity
       call compute_stellar_metallicity(liste_halos(i))
       ! compute mass-weighted stellar age and SFR
       call compute_stellar_age(liste_halos(i))
       ! compute kinematical properties of galaxies
       call compute_stellar_kinematics(liste_halos(i))
#endif

       if(printdatacheckhalo) then
          call cpu_time(t1)
          write(errunit,*) '> halo computation took:',int(t1- t0) ,'s'
          write(errunit,*)
       end if
    end do
    !$OMP END PARALLEL DO
    ! end leo

    ! Need to loop a second time for substructure-based computations
    ! This is because of the OpenMP parallelization
    ! Maybe something like ATOMIC would be better? but also, defeats the purpose.
    write(errunit,*) '> Second loop for inclusive properties'
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) &
    !$OMP PRIVATE(i)
    do i = 1,nb_of_halos + nb_of_subhalos
       ! determine mass and number of particles including substructures
       call det_total_mass(liste_halos(i))
       ! Something for the contaminated fraction:
       ! 3 definitions: P_self, P_host, P_tot
       ! - P_self is the fraction of LR particles
       ! - P_tot is the inclusive fraction (including subhalos)
       ! - P_host is the (max) contam fraction of the host
       ! 1/ compute the mass/nb fraction of LR DM particles in halo, P_self
       ! 2/ compute the same for halos including subhalos, P_tot
       ! 3/ propagate down, P_host: if a l1 halo has contam fraction P1
       !    then the sub will have contam fraction max(P1, p_i)
       ! This mean that:
       ! - a pure host with contaminated subs will be contaminated (P_self=0, P_host=0, P_tot>0)
       ! - a pure sub with a contaminted host will be contaminated (P_self=0, P_host>0, P_tot=0)
    end do
    !$OMP END PARALLEL DO


#ifdef ANG_MOM_OF_R
    close(agor_unit)
#endif

    call write_tree_brick

    if(numero_step.eq.1) then
       open(unit=123,form='formatted',status='unknown',file='info_run.tmp')
       write(123,*) ' massp    :',massp*1.e11
       write(123,*) ' mass min :', real(nMembers)*massp*1.e11
       write(123,*) ' Method   : ', method
       write(123,*) ' cdm      :',cdm
       write(123,*) 'step,aexp,redshift,age_univ,nb_of_halos,nb_of_subhalos'
    end if
    write(123,'(1x,i3,3(1x,E16.4),2(1x,i10))') &
         numstep,aexp,af/aexp - 1.,age_univ,nb_of_halos,nb_of_subhalos
    if(numero_step.eq.nsteps) then
       close(123)
    end if

    deallocate(liste_halos)
    deallocate(nb_of_parts,first_part,linked_list)
    deallocate(pos,vel)
    if(allocated(mass)) deallocate(mass)
#ifdef STARS
    if(allocated(age_st))deallocate(age_st)
    if(allocated(met_st))deallocate(met_st)
    if(allocated(pf_st))deallocate(pf_st)  ! RS - pristie fraction
    if(allocated(pz_st))deallocate(pz_st)  ! RS - primordial Z
    if(allocated(chem_st))deallocate(chem_st)
    call ct_clear_cosmo()
#endif

    if(allocated(family)) deallocate(family)
    if(.not.cdm) deallocate(density)

    call cpu_time(read_time_end)

    write(errunit,*) '> time_step computations took : ',nint(read_time_end - read_time_ini),' seconds'
    write(errunit,*)

    return

  end subroutine new_step

!***********************************************************************
  subroutine make_halos()

  ! subroutine which builds the halos from particle data using fof or adaptahop

    use fof
    use neiKDtree
    implicit none

    real(kind=8)    :: read_time_ini,read_time_end
#ifdef STARS
    integer(i8b)             :: nstar
#endif

    write(errunit,*) '> In routine make_halos '
    write(errunit,*) '> ----------------------'

    write(errunit,*)
    fPeriod    = real(FlagPeriod,8)
    if(FlagPeriod.eq.1) then
       write(errunit,*) '> WARNING: Assuming PERIODIC boundary conditions --> make sure this is correct'
       periodic = .true.
    else
       write(errunit,*) '> WARNING: Assuming NON PERIODIC boundary conditions --> make sure this is correct'
       periodic = .false.
    end if

    if(numero_step.eq.1) then
       if(cdm) then
          write(errunit,*) '> Center of haloes and subhaloes are defined as the particle the closest to the cdm'
       else
          write(errunit,*) '> Center of haloes and subhaloes are defined as the particle with the highest density'
       end if
       select case(method)
       case("FOF")
          write(errunit,*) '> HaloMaker is using Friend Of Friend algorithm'
          call fof_init
          fsub = .false.
       case("HOP")
          write(errunit,*) '> HaloMaker is using Adaptahop in order to'
          write(errunit,*) '> Detect halos, subhaloes will not be selected'
          call init_adaptahop
          fsub = .false.
       case("DPM")
          write(errunit,*) '> HaloMaker is using Adaptahop in order to'
          write(errunit,*) '> Detect halos, and subhaloes with the Density Profile Method'
          call init_adaptahop
          fsub = .true.
       case("MSM")
          write(errunit,*) '> HaloMaker is using Adaptahop in order to'
          write(errunit,*) '> Detect halos, and subhaloes with the Most massive Subhalo Method'
          call init_adaptahop
          fsub = .true.
       case("BHM")
          write(errunit,*) '> HaloMaker is using Adaptahop in order to'
          write(errunit,*) '> Detect halos, and subhaloes with the Branch History Method'
          call init_adaptahop
          fsub = .true.
       case default
          write(errunit,*) '> Selection method: ',method,' is not included'
          write(errunit,*) '> Please check input file:, input_HaloMaker.dat'
       end select
    end if

#ifdef Test_Gd
    write(errunit,*) '> end of test for read_gadget'
    stop
#endif
    call cpu_time(read_time_ini)
    if(method .eq.'FOF') then
       ! compute density when most dense particle is choosen as center
       call fof_main
       if(.not.cdm.and.nb_of_halos.gt.0) call compute_density_for_fof
       nb_of_subhalos = 0
    else
       call compute_adaptahop
       ! no nead to keep information about density when choosing particle closest to the cdm as center
       if(cdm) deallocate(density)
    end if
    call cpu_time(read_time_end)

    write(errunit,'(a,i3,a,i8)') ' > Number of halos with more than     ', nMembers,' particles:',nb_of_halos
    write(errunit,'(a,i3,a,i8)') ' > Number of sub-halos with more than ', nMembers,' particles:',nb_of_subhalos
    write(errunit,*) '> time_step computations took : ',nint(read_time_end - read_time_ini),' seconds'

    return

  end subroutine make_halos

!***********************************************************************
  subroutine make_linked_list()

  ! Subroutine builds a linked list of parts for each halo which
  ! contains all its particles.

    implicit none

    integer(i8b)                :: i,index1,index2
    integer(kind=4)             :: ierr
    integer(i8b),allocatable :: current_ptr(:)

    allocate(current_ptr(0:(nb_of_halos+nb_of_subhalos)),stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) '> cannot allocate current_ptr in make_linked_list'
    endif

    ! initialization of linked list
    first_part(0:(nb_of_halos+nb_of_subhalos))  = -1
    nb_of_parts(0:(nb_of_halos+nb_of_subhalos)) =  0
    current_ptr(0:(nb_of_halos+nb_of_subhalos)) = -1
    linked_list(0:nbodies+1)                    = -1

    ! make linked list: a few (necessary) explanations ....
    !   1/ index1 (liste_parts(i)) is the number of the halo to which belongs particle i (i is in [1..nbodies])
    !   2/ first_part(j) is the number of the first particle of halo j  (j is in [1..nhalo])
    !   3/ current_ptr(index1) is the number of the latest particle found in halo number index1
    ! to sum up, first_part(i) contains the number of the first particle of halo i,
    ! linked_list(first_part(i)) the number of the second particle of halo i, etc ... until
    ! the last particle which points to number -1

    do i=1,nbodies
       index1 = liste_parts(i)
       if(index1.gt.(nb_of_halos+nb_of_subhalos)) stop 'error in liste_parts'
       if (first_part(index1) == -1) then
          first_part(index1)  = i
          nb_of_parts(index1) = 1
          current_ptr(index1) = i
       else
          index2              = current_ptr(index1)
          linked_list(index2) = i
          current_ptr(index1) = i
          nb_of_parts(index1) = nb_of_parts(index1)+1
       endif
    end do

    ! close linked list
    do i=0,(nb_of_halos + nb_of_subhalos)
       if (current_ptr(i) == -1) cycle
       index2              = current_ptr(i)
       linked_list(index2) = -1
    end do

    deallocate(current_ptr)

    return

  end subroutine make_linked_list

!*************************************************************************
  subroutine init_halos

    implicit none
    integer(kind=4) :: ihalo,ihtmp,imother

    do ihalo = 1, nb_of_halos + nb_of_subhalos

       call clear_halo(liste_halos(ihalo))
       liste_halos(ihalo)%my_number   = ihalo
       liste_halos(ihalo)%my_timestep = numstep
       if(fsub) then
          if(level(ihalo).eq.1) then
             if(first_daughter(ihalo).gt.0) then
                liste_halos(ihalo)%nextsub = first_daughter(ihalo)
             end if
             liste_halos(ihalo)%hosthalo   = ihalo
          else
             liste_halos(ihalo)%level    = level(ihalo)
             liste_halos(ihalo)%hostsub = mother(ihalo)
             imother = liste_halos(ihalo)%hostsub
             liste_halos(imother)%nbsub = liste_halos(imother)%nbsub + 1
             if(first_daughter(ihalo).gt.0) then
                liste_halos(ihalo)%nextsub = first_daughter(ihalo)
             else if(first_sister(ihalo).gt.0) then
                liste_halos(ihalo)%nextsub = first_sister(ihalo)
             else
                ihtmp = ihalo
                do while((first_sister(ihtmp).le.0).and.(level(ihtmp).gt.1))
                   ihtmp = mother(ihtmp)
                end do
                if(level(ihtmp).gt.1) then
                   ihtmp = first_sister(ihtmp)
                   liste_halos(ihalo)%nextsub = ihtmp
                end if
             end if
             imother = liste_halos(ihalo)%hostsub
             imother = liste_halos(imother)%hosthalo
             if(liste_halos(imother)%level.ne.1) stop 'wrong id for halo host'
             liste_halos(ihalo)%hosthalo = imother
          end if
       end if
    end do

    if(fsub) then
       deallocate(mother,first_daughter,first_sister,level)
    end if
    return

  end subroutine init_halos

!*************************************************************************
  subroutine compute_density_for_fof

    use neikdtree
    implicit none

    omegaL   = omega_lambda_f
    omega0   = omega_f
    aexp_max = af
    hubble   = H_f*1d-2
    boxsize2 = Lf
    xlong    = boxsize2
    ylong    = xlong
    zlong    = xlong
    xlongs2  = xlong*0.5d0
    ylongs2  = ylong*0.5d0
    zlongs2  = zlong*0.5d0

    pos_renorm   = xlong
    pos_shift(1) = 0.0d0
    pos_shift(2) = 0.0d0
    pos_shift(3) = 0.0d0
    mass_in_kg   = (xlong*ylong*zlong/dble(npart))*mega_parsec**3 &
         &       * omega0*critical_density*hubble**2
    Hub_pt       = 100.*hubble * sqrt(omega0*(aexp_max/aexp)**3   &
         &       + (1-omegaL-omega0)*(aexp_max/aexp)**2 + omegaL)

    nvoisins = nMembers
    nhop     = nMembers
    write(errunit,*) '> Computing density using adatahop routines'
    call change_pos
    ! action "neighbors
    call create_tree_structure
    call compute_mean_density_and_np
    call change_pos_back
    deallocate(iparneigh)

    return

  end subroutine compute_density_for_fof

!*************************************************************************
  subroutine det_mass(h)

    ! adds up masses of particles to get total mass of the halo.

    implicit none
    integer(i8b)       :: indexp,npch
    real(kind=8)       :: masshalo, mstar, mdm
    integer(i8b)       :: nstar, ndm
    type(halo)         :: h
#ifdef CONTAM
    real(kind=8)       :: mcontam
    integer(i8b)       :: ncontam
#endif

    masshalo      = 0d0
    mstar         = 0d0
    mdm           = 0d0
    npch          = 0
    nstar         = 0
    ndm           = 0
    indexp        = first_part(h%my_number)
#ifdef CONTAM
    ncontam = 0
    mcontam = 0d0
#endif
    do while(indexp /= -1)
       if (allocated(mass)) then
          masshalo = masshalo + real(mass(indexp),8)
#ifdef ALLPARTS
          if (family(indexp) == FAM_STAR) then
             mstar = mstar + real(mass(indexp),8)
             nstar = nstar + 1
          else
             mdm = mdm + real(mass(indexp),8)
             ndm = ndm + 1
          end if
#endif
       else
          masshalo = masshalo + real(massp,8)
       endif
#ifdef CONTAM
       if ((mass(indexp) > massp* 1.00001) .and. family(indexp) == FAM_DM) then
          mcontam = mcontam + real(mass(indexp),8)
          ncontam = ncontam + 1
       end if
#endif
       npch = npch + 1
       indexp = linked_list(indexp)
    end do
    h%m = real(masshalo,8)  ! in 10^11 M_sun
#ifdef ALLPARTS
    h%m_star = mstar
    h%n_star = nstar
    h%m_DM   = mdm
    h%n_DM   = ndm
#endif
#ifdef CONTAM
    h%m_contam = mcontam
    h%n_contam = ncontam
#endif

    if(npch.ne.nb_of_parts(h%my_number)) then
       write(errunit,*) '> Fatal error in det_mass for', h%my_number
       write(errunit,*) 'nb_of_parts, npch:',h%my_number,npch
       stop
    end if

    return

  end subroutine det_mass

!***********************************************************************
  subroutine det_total_mass(h)

    ! adds up masses of structures and substructures

    implicit none
    integer(kind=4)    :: ih
    type(halo)         :: h, hsub

    ih = h%my_number
    do while (ih > 0)
       hsub = liste_halos(ih)
       ! Do not include "sister" structures...
       if ((hsub%my_number .ne. h%my_number) .and. (hsub%level .le. h%level)) exit

       h%mtot = h%mtot + hsub%m
       h%ntot = h%ntot + nb_of_parts(ih)
#ifdef ALLPARTS
       ! Masses
       h%mtot_star = h%mtot_star + hsub%m_star
       h%mtot_DM = h%mtot_DM + hsub%m_DM
       ! Particle numbers
       h%ntot_star = h%ntot_star + hsub%n_star
       h%ntot_DM = h%ntot_DM + hsub%n_DM
#endif
#ifdef CONTAM
       h%mtot_contam = h%mtot_contam + hsub%m_contam
       h%ntot_contam = h%ntot_contam + hsub%n_contam
#endif

       ih = hsub%nextsub
    end do

    return

  end subroutine det_total_mass

!***********************************************************************
  subroutine compute_ang_mom(h)

  ! compute angular momentum of all halos

    implicit none

    integer(i8b)    :: indexp
    real(kind=8)    :: lx,ly,lz
    real(kind=8)    :: lxstar, lystar, lzstar, lxdm, lydm, lzdm
    type (halo)     :: h
    type (vector)   :: dr,p
#ifdef ALLPARTS
    type (vector)   :: drtype,ptype
#endif

    ! we compute r * m * v, where r & v are pos and vel of halo particles relative to center of halo
    ! (particle closest to center of mass or most dense particle)

    indexp = first_part(h%my_number)
    lx =0d0 ; ly = 0d0 ; lz = 0d0
    lxstar =0d0 ; lystar = 0d0 ; lzstar = 0d0
    lxdm =0d0 ; lydm = 0d0 ; lzdm = 0d0
    do while(indexp /= -1)

       dr%x   = pos(indexp,1) - h%p%x
       dr%y   = pos(indexp,2) - h%p%y
       dr%z   = pos(indexp,3) - h%p%z
       call correct_for_periodicity(dr)
#ifdef ALLPARTS
       if (family(indexp) == FAM_STAR) then
          drtype%x   = pos(indexp,1) - h%p_star%x
          drtype%y   = pos(indexp,2) - h%p_star%y
          drtype%z   = pos(indexp,3) - h%p_star%z
       else
          drtype%x   = pos(indexp,1) - h%p_DM%x
          drtype%y   = pos(indexp,2) - h%p_DM%y
          drtype%z   = pos(indexp,3) - h%p_DM%z
       end if
       call correct_for_periodicity(drtype)
#endif


       if (allocated(mass)) then
          p%x = mass(indexp)*(vel(indexp,1)-h%v%x)
          p%y = mass(indexp)*(vel(indexp,2)-h%v%y)
          p%z = mass(indexp)*(vel(indexp,3)-h%v%z)
       else
          p%x = massp*(vel(indexp,1)-h%v%x)
          p%y = massp*(vel(indexp,2)-h%v%y)
          p%z = massp*(vel(indexp,3)-h%v%z)
       endif
#ifdef ALLPARTS
       if (family(indexp) == FAM_STAR) then
          ptype%x = mass(indexp)*(vel(indexp,1)-h%v_star%x)
          ptype%y = mass(indexp)*(vel(indexp,2)-h%v_star%y)
          ptype%z = mass(indexp)*(vel(indexp,3)-h%v_star%z)
       else
          ptype%x = mass(indexp)*(vel(indexp,1)-h%v_DM%x)
          ptype%y = mass(indexp)*(vel(indexp,2)-h%v_DM%y)
          ptype%z = mass(indexp)*(vel(indexp,3)-h%v_DM%z)
       end if
#endif

       lx  = lx + real(dr%y*p%z - dr%z*p%y,8)   ! in 10**11 Msun * km/s * Mpc
       ly  = ly + real(dr%z*p%x - dr%x*p%z,8)
       lz  = lz + real(dr%x*p%y - dr%y*p%x,8)
#ifdef ALLPARTS
          if (family(indexp) == FAM_STAR) then
             lxstar  = lxstar + real(drtype%y*ptype%z - drtype%z*ptype%y,8)
             lystar  = lystar + real(drtype%z*ptype%x - drtype%x*ptype%z,8)
             lzstar  = lzstar + real(drtype%x*ptype%y - drtype%y*ptype%x,8)
          else
             lxdm  = lxdm + real(drtype%y*ptype%z - drtype%z*ptype%y,8)
             lydm  = lydm + real(drtype%z*ptype%x - drtype%x*ptype%z,8)
             lzdm  = lzdm + real(drtype%x*ptype%y - drtype%y*ptype%x,8)
          end if
#endif
       indexp = linked_list(indexp)

    end do

    h%L%x = real(lx,8)
    h%L%y = real(ly,8)
    h%L%z = real(lz,8)
#ifdef ALLPARTS
    h%L_star%x = real(lxstar,8)
    h%L_star%y = real(lystar,8)
    h%L_star%z = real(lzstar,8)
    h%L_dm%x = real(lxdm,8)
    h%L_dm%y = real(lydm,8)
    h%L_dm%z = real(lzdm,8)
#endif
    return

  end subroutine compute_ang_mom

!***********************************************************************
  subroutine r_halos(h)

  ! compute distance of the most remote particle (with respect to center of halo, which
  ! is either center of mass or most bound particle)

    implicit none

    integer(i8b)    :: indexp
    real(kind=8)    :: dr2max,dr2, dr2star, dr2dm
    type (vector)   :: dr
    type (halo)     :: h

    dr2max  = 0.d0
    dr2star = 0.d0
    dr2dm   = 0.d0
    indexp = first_part(h%my_number)

    do while(indexp /= -1)

       dr%x = pos(indexp,1) - h%p%x
       dr%y = pos(indexp,2) - h%p%y
       dr%z = pos(indexp,3) - h%p%z
       call correct_for_periodicity(dr)
       dr2  = (dr%x*dr%x + dr%y*dr%y + dr%z*dr%z)
       if (dr2 > dr2max) then
          dr2max         = dr2
       endif

#ifdef ALLPARTS
       if (family(indexp) == FAM_STAR) then
          dr%x = pos(indexp,1) - h%p_star%x
          dr%y = pos(indexp,2) - h%p_star%y
          dr%z = pos(indexp,3) - h%p_star%z
       else
          dr%x = pos(indexp,1) - h%p_DM%x
          dr%y = pos(indexp,2) - h%p_DM%y
          dr%z = pos(indexp,3) - h%p_DM%z
       end if
       call correct_for_periodicity(dr)
       dr2 = (dr%x*dr%x + dr%y*dr%y + dr%z*dr%z)

       if ((family(indexp) == FAM_STAR) .and. (dr2 > dr2star)) then
          dr2star        = dr2
       endif
       if ((family(indexp) == FAM_DM) .and. (dr2 > dr2dm)) then
          dr2dm          = dr2
       endif
#endif

       indexp=linked_list(indexp)

    end do

    h%r = sqrt(dr2max)
#ifdef ALLPARTS
    h%r_star = sqrt(dr2star)
    h%r_DM = sqrt(dr2dm)
#endif
    return

  end subroutine r_halos

!***********************************************************************
  subroutine correct_for_periodicity(dr)

  ! subroutine corrects for the fact that if you have periodic boundary conditions,
  ! then groups of particles that sit on the edge of the box can have part of their
  ! particles on one side and part on the other. So we have to take out a box-length
  ! when measuring the distances between group members if needed.

    implicit none

    type (vector) :: dr

    if (FlagPeriod == 0) return  !--> NO PERIODIC BCs

    if (dr%x > + Lbox_pt2) dr%x = dr%x - Lbox_pt
    if (dr%x < - Lbox_pt2) dr%x = dr%x + Lbox_pt

    if (dr%y > + Lbox_pt2) dr%y = dr%y - Lbox_pt
    if (dr%y < - Lbox_pt2) dr%y = dr%y + Lbox_pt

    if (dr%z > + Lbox_pt2) dr%z = dr%z - Lbox_pt
    if (dr%z < - Lbox_pt2) dr%z = dr%z + Lbox_pt

    return

  end subroutine correct_for_periodicity

!***********************************************************************
  subroutine det_center(h)

  ! compute position of center of mass of halo, and its velocity.

    implicit none

    type (halo)        :: h
    integer(i8b)       :: indexp, icenter,ifirst
    real(kind=8)       :: maxdens, distmin
    real(kind=8)       :: pcx,pcy,pcz,vcx,vcy,vcz,vmean,v2mean,sigma2,vnorm
    real(kind=8)       :: distr
#ifdef ALLPARTS
    real(kind=8)       :: pcxstar,pcystar,pczstar,pcxDM,pcyDM,pczDM
    real(kind=8)       :: vcxstar,vcystar,vczstar,vcxDM,vcyDM,vczDM
    real(kind=8)       :: vmeanstar, vmeanDM
    type(vector)       :: pcstar, pcDM
    integer(i8b)       :: icenterDM, icenterstar
    real(kind=8)       :: maxdensDM, maxdensstar
    real(kind=8)       :: distminDM, distminstar
#endif
    real(kind=8)       :: sigma2_star, sigma2_dm
    type(vector)       :: dr,pc
    integer(i8b)       :: i,nmax

    icenter = -1
#ifdef ALLPARTS
    icenterDM = -1
    icenterstar = -1
#endif

    if (cdm) then

       ! compute cdm
       pcx   = 0d0 ; pcy   = 0d0 ; pcz = 0d0
#ifdef ALLPARTS
       pcxstar   = 0d0 ; pcystar   = 0d0 ; pczstar = 0d0
       pcxDM   = 0d0 ; pcyDM   = 0d0 ; pczDM = 0d0
#endif
       ifirst = first_part(h%my_number)
       indexp = ifirst
       do while (indexp /= -1)
          dr%x = pos(indexp,1) - pos(ifirst,1)
          dr%y = pos(indexp,2) - pos(ifirst,2)
          dr%z = pos(indexp,3) - pos(ifirst,3)
          call correct_for_periodicity(dr)
          if (allocated(mass)) then
             pcx = pcx + real(mass(indexp)*dr%x,8)
             pcy = pcy + real(mass(indexp)*dr%y,8)
             pcz = pcz + real(mass(indexp)*dr%z,8)
#ifdef ALLPARTS
             if (family(indexp)==FAM_STAR) then
                pcxstar = pcxstar + real(mass(indexp)*dr%x,8)
                pcystar = pcystar + real(mass(indexp)*dr%y,8)
                pczstar = pczstar + real(mass(indexp)*dr%z,8)
             else
                pcxDM = pcxDM + real(mass(indexp)*dr%x,8)
                pcyDM = pcyDM + real(mass(indexp)*dr%y,8)
                pczDM = pczDM + real(mass(indexp)*dr%z,8)
             endif
#endif
          else
             pcx = pcx + real(massp*dr%x,8)
             pcy = pcy + real(massp*dr%y,8)
             pcz = pcz + real(massp*dr%z,8)
          end if
          indexp = linked_list(indexp)
       end do
       pcx  = pcx / real(h%m,8) + real(pos(ifirst,1),8)
       pcy  = pcy / real(h%m,8) + real(pos(ifirst,2),8)
       pcz  = pcz / real(h%m,8) + real(pos(ifirst,3),8)
       pc%x = real(pcx,8)
       pc%y = real(pcy,8)
       pc%z = real(pcz,8)
       call correct_for_periodicity(pc)
#ifdef ALLPARTS
       if (h%n_star >= nMembers) then
          pcstar%x = real(pcxstar,8) / real(h%m_star,8) + real(pos(ifirst,1),8)
          pcstar%y = real(pcystar,8) / real(h%m_star,8) + real(pos(ifirst,2),8)
          pcstar%z = real(pczstar,8) / real(h%m_star,8) + real(pos(ifirst,3),8)
          call correct_for_periodicity(pcstar)
       end if
       if (h%n_DM >= nMembers) then
          pcDM%x = real(pcxDM,8) / real(h%m_DM,8) + real(pos(ifirst,1),8)
          pcDM%y = real(pcyDM,8) / real(h%m_DM,8) + real(pos(ifirst,2),8)
          pcDM%z = real(pczDM,8) / real(h%m_DM,8) + real(pos(ifirst,3),8)
          call correct_for_periodicity(pcDM)
       end if
#endif
       ! search particule closest to the cdm
       indexp  = ifirst
       distmin = Lbox_pt
#ifdef ALLPARTS
       distminDM = Lbox_pt
       distminstar = Lbox_pt
#endif
       do while (indexp /= -1)
          dr%x = pos(indexp,1) - pc%x
          dr%y = pos(indexp,2) - pc%y
          dr%z = pos(indexp,3) - pc%z
          call correct_for_periodicity(dr)
          distr = sqrt(dr%x**2+dr%y**2+dr%z**2)
          if (distr.lt.distmin) then
             icenter = indexp
             distmin = distr
          end if
#ifdef ALLPARTS
          if (family(indexp) == FAM_DM) then
             dr%x = pos(indexp,1) - pcDM%x
             dr%y = pos(indexp,2) - pcDM%y
             dr%z = pos(indexp,3) - pcDM%z
             call correct_for_periodicity(dr)
             distr = sqrt(dr%x**2+dr%y**2+dr%z**2)
             if (distr.lt.distminDM) then
                icenterDM = indexp
                distminDM = distr
             end if
          end if
          if (family(indexp) == FAM_STAR) then
             dr%x = pos(indexp,1) - pcstar%x
             dr%y = pos(indexp,2) - pcstar%y
             dr%z = pos(indexp,3) - pcstar%z
             call correct_for_periodicity(dr)
             distr = sqrt(dr%x**2+dr%y**2+dr%z**2)
             if (distr.lt.distminstar) then
                icenterstar = indexp
                distminstar = distr
             end if
          endif
#endif
          indexp = linked_list(indexp)
       end do

    else
       maxdens = 0.d0
#ifdef ALLPARTS
       maxdensDM = 0d0
       maxdensstar = 0d0
#endif
       indexp  = first_part(h%my_number)
       do while (indexp /= -1)
          if (density(indexp).gt.maxdens) then
             maxdens = density(indexp)
             icenter = indexp
          end if
#ifdef ALLPARTS
          if ((density(indexp).gt.maxdensstar) .and. (family(indexp)==FAM_STAR)) then
             maxdensstar = density(indexp)
             icenterstar = indexp
          end if
          if ((density(indexp).gt.maxdensDM) .and. (family(indexp)==FAM_DM)) then
             maxdensDM = density(indexp)
             icenterDM = indexp
          end if
#endif
          indexp = linked_list(indexp)
       enddo

    end if


    if (icenter.lt.0) then
       write(errunit,*) '> Could not find a center for halo: ',h%my_number,icenter
       write(errunit,*) '  h%m,massp,h%m/massp             : ',h%m,massp,h%m/massp
       write(errunit,*) '  Lbox_pt,distmin                 : ',Lbox_pt,distmin
       write(errunit,*) '  pcx,pcy,pcz                  : ',pcx,pcy,pcz
       write(errunit,*) '  periodicity flag                : ',FlagPeriod
       write(errunit,*) '> Check routine det_center'
       stop
    end if

    h%p%x  = pos(icenter,1)
    h%p%y  = pos(icenter,2)
    h%p%z  = pos(icenter,3)
#ifdef ALLPARTS
    ! Default centre is identical for DM, stars and total
    h%p_star%x=h%p%x
    h%p_star%y=h%p%y
    h%p_star%z=h%p%z
    h%p_DM%x=h%p%x
    h%p_DM%y=h%p%y
    h%p_DM%z=h%p%z
    if (h%n_DM >= nMembers) then
       h%p_DM%x  = pos(icenterDM,1)
       h%p_DM%y  = pos(icenterDM,2)
       h%p_DM%z  = pos(icenterDM,3)
    endif
    if (h%n_star >= nMembers) then
       h%p_star%x  = pos(icenterstar,1)
       h%p_star%y  = pos(icenterstar,2)
       h%p_star%z  = pos(icenterstar,3)
    end if
#endif

    ! velocity of center is set equal velocity of center of mass:

    indexp = first_part(h%my_number)
    vcx = 0d0 ; vcy = 0d0 ; vcz =0d0
    v2mean= 0d0
#ifdef ALLPARTS
       pcxstar   = 0d0 ; pcystar   = 0d0 ; pczstar = 0d0
       pcxDM   = 0d0 ; pcyDM   = 0d0 ; pczDM = 0d0
       vcxDM   = 0d0 ; vcyDM   = 0d0 ; vczDM = 0d0
       vcxstar   = 0d0 ; vcystar   = 0d0 ; vczstar = 0d0
       vmeanstar = 0d0; vmeanDM = 0d0;
#endif
    i=0
    do while (indexp /= -1)
       i=i+1
       if (allocated(mass)) then
          vcx = vcx + real(mass(indexp)*vel(indexp,1),8)
          vcy = vcy + real(mass(indexp)*vel(indexp,2),8)
          vcz = vcz + real(mass(indexp)*vel(indexp,3),8)
#ifdef ALLPARTS
          if (family(indexp)==FAM_STAR) then
             vcxstar = vcxstar + real(mass(indexp)*vel(indexp,1),8)
             vcystar = vcystar + real(mass(indexp)*vel(indexp,2),8)
             vczstar = vczstar + real(mass(indexp)*vel(indexp,3),8)
          else
             vcxDM = vcxDM + real(mass(indexp)*vel(indexp,1),8)
             vcyDM = vcyDM + real(mass(indexp)*vel(indexp,2),8)
             vczDM = vczDM + real(mass(indexp)*vel(indexp,3),8)
          end if
#endif
       else
          vcx = vcx + real(massp*vel(indexp,1),8)
          vcy = vcy + real(massp*vel(indexp,2),8)
          vcz = vcz + real(massp*vel(indexp,3),8)
       endif
       indexp   = linked_list(indexp)
    end do
    nmax=i

    h%v%x = real(vcx,8)/h%m
    h%v%y = real(vcy,8)/h%m
    h%v%z = real(vcz,8)/h%m
    vmean=sqrt(vcx**2+vcy**2+vcz**2)/h%m

#ifdef ALLPARTS
    ! Default centre is identical for DM, stars and total
    h%v_star%x=h%v%x
    h%v_star%y=h%v%y
    h%v_star%z=h%v%z
    h%v_DM%x=h%v%x
    h%v_DM%y=h%v%y
    h%v_DM%z=h%v%z
    if (h%n_DM >= nMembers) then
       h%v_DM%x = real(vcxDM,8)/h%m_DM
       h%v_DM%y = real(vcyDM,8)/h%m_DM
       h%v_DM%z = real(vczDM,8)/h%m_DM
       vmeanDM  = sqrt(vcxDM**2+vcyDM**2+vczDM**2)/h%m_DM
    endif
    if (h%n_star >= nMembers) then
       h%v_star%x = real(vcxstar,8)/h%m_star
       h%v_star%y = real(vcystar,8)/h%m_star
       h%v_star%z = real(vczstar,8)/h%m_star
       vmeanstar  = sqrt(vcxstar**2+vcystar**2+vczstar**2)/h%m_star
    end if
#endif


    indexp = first_part(h%my_number)
    sigma2= 0d0
    i=0
    sigma2_star = 0d0
    sigma2_dm   = 0d0
    do while (indexp /= -1)
       i=i+1
       vnorm = real(sqrt(vel(indexp,1)**2+vel(indexp,2)**2+vel(indexp,3)**2),8)
       if (allocated(mass)) then
          sigma2 = sigma2 + mass(indexp)*(vnorm-vmean)**2
#ifdef ALLPARTS
          if (family(indexp) == FAM_STAR) then
             sigma2_star = sigma2_star + mass(indexp)*(vnorm-vmeanstar)**2
          else
             sigma2_dm = sigma2_dm + mass(indexp)*(vnorm-vmeanDM)**2
          end if
#endif
       else
          sigma2 = sigma2 + massp       *(vnorm-vmean)**2
       endif
       indexp   = linked_list(indexp)
    end do

    ! Normalize to the mass
    h%sigma = real(sqrt(sigma2/h%m),8)
#ifdef ALLPARTS
    if (h%n_star>=nMembers) h%sigma_star = real(sqrt(sigma2_star/h%m_star),8)
    if (h%n_DM >= nMembers) h%sigma_DM = real(sqrt(sigma2_dm/h%m_dm),8)
#endif


    return

  end subroutine det_center

!***********************************************************************
  subroutine compute_3d_profile(h)

    ! compute 3D profile of particles (and maybe stars and DM)

    implicit none

    type (halo)        :: h
    integer(i8b)       :: indexp
    real(kind=8)       :: rx,ry,rz,rr
    integer(i8b)       :: j, j50, j90
    real(kind=8)       :: dx, dxstar, dxDM, m_in
    real(kind=8),dimension(1:h%nbin):: m_tot,dv, mstar_tot, mDM_tot, dvstar, dvDM

    dx = h%r/dble(h%nbin)
#ifdef ALLPARTS
    dxstar = min(h%r_star, h%datas%rvir)/dble(h%nbin)
    dxDM   = h%r_DM/dble(h%nbin)
    mDM_tot=0d0
    mstar_tot =0d0
#endif

    indexp = first_part(h%my_number)
    m_tot=0d0
    do while (indexp /= -1)
       rx=pos(indexp,1)-h%p%x
       ry=pos(indexp,2)-h%p%y
       rz=pos(indexp,3)-h%p%z
       rr=sqrt(rx**2d0+ry**2d0+rz**2d0)

       j=int(rr/dx)+1
       if (j .le. h%nbin)then
          if (allocated(mass)) then
             m_tot(j)=m_tot(j)+mass(indexp)
          else
             m_tot(j)=m_tot(j)+massp
          endif
       endif
#ifdef ALLPARTS
       if ((family(indexp) == FAM_STAR) .and. (h%n_star >= nMembers)) then
          rx=pos(indexp,1)-h%p_star%x
          ry=pos(indexp,2)-h%p_star%y
          rz=pos(indexp,3)-h%p_star%z
          rr=sqrt(rx**2d0+ry**2d0+rz**2d0)
          j=int(rr/dxstar)+1
          if (j .le. h%nbin) then
             mstar_tot(j) = mstar_tot(j) + mass(indexp)
          endif
       else if ((family(indexp) == FAM_DM) .and. (h%n_DM >= nMembers)) then
          rx=pos(indexp,1)-h%p_DM%x
          ry=pos(indexp,2)-h%p_DM%y
          rz=pos(indexp,3)-h%p_DM%z
          rr=sqrt(rx**2d0+ry**2d0+rz**2d0)
          j=int(rr/dxDM)+1
          if (j .le. h%nbin)then
             mDM_tot(j) = mDM_tot(j) + mass(indexp)
          endif
       endif
#endif
       indexp   = linked_list(indexp)
    end do

    do j=1,h%nbin
       h%rr3D(j)=(j-0.5)*dx
#ifdef ALLPARTS
       if (h%n_star >= nMembers) h%rr3D_star(j) = (j-0.5)*dxstar
       if (h%n_DM >= nMembers)   h%rr3D_DM(j) = (j-0.5)*dxDM
#endif
    enddo


    ! Compute R50, R90
    ! TODO: could be done without binning, but would require to reorder the particles
    m_in = 0d0
    j50 = -1
    j90 = -1
    do j=1,h%nbin
       m_in = m_in + m_tot(j)
       if ((m_in > sum(m_tot)*0.5d0) .and. j50<0) j50 = j
       if ((m_in > sum(m_tot)*0.9d0) .and. j90<0) j90 = j
    end do
    if(j50<0) j50 = h%nbin
    if(j90<0) j90 = h%nbin
    h%r50 = (j50 - 0.5d0)*dx
    h%r90 = (j90 - 0.5d0)*dx
#ifdef ALLPARTS
    ! Same for stars
    if (h%n_star >= nMembers) then
       m_in = 0d0
       j50 = -1
       j90 = -1
       do j=1,h%nbin
          m_in = m_in + mstar_tot(j)
          if ((m_in > sum(mstar_tot)*0.5d0) .and. (j50<0)) j50 = j
          if ((m_in > sum(mstar_tot)*0.9d0) .and. (j90<0)) j90 = j
       end do
       if(j50<0) j50 = h%nbin
       if(j90<0) j90 = h%nbin
       h%r50_star = (j50 - 0.5d0)*dxstar
       h%r90_star = (j90 - 0.5d0)*dxstar
    else
       h%r50_star = 0d0
       h%r90_star = 0d0
    end if
    ! And for DM
    if (h%n_DM >= nMembers) then
       m_in = 0d0
       j50 = -1
       j90 = -1
       do j=1,h%nbin
          m_in = m_in + mDM_tot(j)
          if ((m_in > sum(mDM_tot)*0.5d0) .and. j50<0) j50 = j
          if ((m_in > sum(mDM_tot)*0.9d0) .and. j90<0) j90 = j
       end do
       if(j50<0) j50 = h%nbin
       if(j90<0) j90 = h%nbin
       h%r50_DM = (j50 - 0.5d0)*dxDM
       h%r90_DM = (j90 - 0.5d0)*dxDM
    else
       h%r50_DM = 0d0
       h%r90_DM = 0d0
    end if
#endif



    ! Compute densities
    do j=1,h%nbin
       ! volume of the shell is 4/3 pi * (r_{i+1}^3 - r_i^3)
       dv(j)=4.d0/3.d0*acos(-1d0)*( (h%rr3D(j)+0.5d0*dx)**3d0 - (h%rr3D(j)-0.5d0*dx)**3d0 )
#ifdef ALLPARTS
       dvstar(j)=4.d0/3.d0*acos(-1d0)*( (h%rr3D_star(j)+0.5d0*dx)**3d0 - (h%rr3D_star(j)-0.5d0*dx)**3d0 )
       dvDM(j)=4.d0/3.d0*acos(-1d0)*( (h%rr3D_DM(j)+0.5d0*dx)**3d0 - (h%rr3D_DM(j)-0.5d0*dx)**3d0 )
#endif
    enddo

    h%rr3D=h%rr3D*1d3                 ! kpc
    h%rho3D=m_tot*1d11/(dv*1d6**3d0)  ! Msun/pc**3
#ifdef ALLPARTS
    if (h%n_star >= nMembers) then
       h%rr3D_star  = h%rr3D_star*1d3
       h%rho3D_star = mstar_tot*1d11/(dvstar*1d6**3d0)
    else
       h%rr3D_star (1:h%nbin) = 0d0
       h%rho3D_star(1:h%nbin) = 0d0
    end if
    if (h%n_DM >= nMembers) then
       h%rr3D_DM = h%rr3D_DM*1d3
       h%rho3D_DM = mDM_tot*1d11/(dvDM*1d6**3d0)
    else
       h%rr3D_DM (1:h%nbin) = 0d0
       h%rho3D_DM(1:h%nbin) = 0d0
    end if
#endif

    return

  end subroutine compute_3d_profile


#ifdef STARS
!***********************************************************************
  subroutine compute_bulge_prop(h)

  ! compute position of center of mass of halo, and its velocity.

    implicit none

    type (halo)        :: h
    integer(i8b)       :: indexp
    real(kind=8)       :: vcz,vcr,vct,vmean,sigma2,vnorm
    real(kind=8)       :: Lx,Ly,Lz,Ltot,mbulge,dxx,dyy,dzz,dzcyl
    real(kind=8)       :: planecen_x,planecen_y,planecen_z,rx,ry,rz,modr
    real(kind=8)       :: thetax,thetay,thetaz
    integer(i8b)       :: i
    real(kind=8)       :: xc, yc, zc

#ifdef ALLPARTS
    Ltot=sqrt(h%L_star%x**2+h%L_star%y**2+h%L_star%z**2)
    Lx=h%L_star%x/Ltot
    Ly=h%L_star%y/Ltot
    Lz=h%L_star%z/Ltot
    xc = h%p_star%x
    yc = h%p_star%y
    zc = h%p_star%z
    ! Assume that the mean bulge velocity is the mean galaxy velocity
    vmean=sqrt(h%v_star%x**2+h%v_star%y**2+h%v_star%z**2)
#else
    Ltot=sqrt(h%L%x**2+h%L%y**2+h%L%z**2)
    Lx=h%L%x/Ltot
    Ly=h%L%y/Ltot
    Lz=h%L%z/Ltot
    xc = h%p%x
    yc = h%p%y
    zc = h%p%z
    ! Assume that the mean bulge velocity is the mean galaxy velocity
    vmean=sqrt(h%v%x**2+h%v%y**2+h%v%z**2)
#endif


    indexp = first_part(h%my_number)
    mbulge= 0d0
    sigma2= 0d0
    i=0
    do while (indexp /= -1)
       i=i+1
       if (family(indexp) == FAM_STAR) then
          dxx=pos(indexp,1)-xc
          dyy=pos(indexp,2)-yc
          dzz=pos(indexp,3)-zc
          dzcyl=dxx*Lx+dyy*Ly+dzz*Lz
          ! Compute centre coord. of the plane where the particle lies
          ! and perp. to the disc plane
          planecen_x=xc+dzcyl*Lx
          planecen_y=yc+dzcyl*Ly
          planecen_z=zc+dzcyl*Lz
          ! Compute unit radial vector from this centre to the particle
          rx=pos(indexp,1)-planecen_x
          ry=pos(indexp,2)-planecen_y
          rz=pos(indexp,3)-planecen_z
          modr=sqrt(rx**2d0+ry**2d0+rz**2d0)

          ! Check if the vector is defined
          if (modr .ne. 0d0)then
             rx=rx/modr
             ry=ry/modr
             rz=rz/modr
          endif

          ! Compute unit vector perp. to this vector and in the same plane
          thetax=Ly*rz-Lz*ry
          thetay=Lz*rx-Lx*rz
          thetaz=Lx*ry-Ly*rx

          ! Compute velocity components in this cylindrical coord.
          vcz = vel(indexp,1)*Lx     + vel(indexp,2)*Ly     + vel(indexp,3)*Lz
          vcr = vel(indexp,1)*rx     + vel(indexp,2)*ry     + vel(indexp,3)*rz
          vct = vel(indexp,1)*thetax + vel(indexp,2)*thetay + vel(indexp,3)*thetaz

          ! Check if the particle belongs to the bulge
          !if(vct .lt. vcr .and. vct .lt. vcz) then
          if(abs(vct) .lt. 1d0*sqrt(vcr**2+vcz**2)) then
             vnorm = real(sqrt(vcz**2+vcr**2+vct**2),8)
             if (allocated(mass)) then
                mbulge=mbulge+mass(indexp)
                ! Compute the velocity disperson in the bulge
                sigma2 = sigma2 + mass(indexp)*(vnorm-vmean)**2
             else
                mbulge=mbulge+massp
                ! Compute the velocity disperson in the bulge
                sigma2 = sigma2 + massp       *(vnorm-vmean)**2
             endif
          endif
       end if
       indexp   = linked_list(indexp)
    end do

    h%sigma_bulge = real(sqrt(sigma2/mbulge),8)
    h%m_bulge     = real(mbulge,8)

    return

  end subroutine compute_bulge_prop

!***********************************************************************
  subroutine compute_stellar_profile(h)

  ! compute position of center of mass of halo, and its velocity.

    implicit none

    type (halo)        :: h
    integer(i8b)       :: indexp
    real(kind=8)       :: xc, yc, zc
    real(kind=8)       :: Lx,Ly,Lz,Ltot,dxx,dyy,dzz,dzcyl
    real(kind=8)       :: planecen_x,planecen_y,planecen_z,rx,ry,rz,modr
    integer(i8b)       :: i,j,jeff
    real(kind=8)       :: rmax,dx,m_in
    real(kind=8),dimension(1:h%nbin):: m_tot,dv

#ifdef ALLPARTS
    ! Set default values for too small galaxies
    if ((h%n_star < nMembers) .or. (h%r_star == 0d0)) then
       h%datas%Reff = 0d0
       h%rr(1:h%nbin)=0d0
       h%rho(1:h%nbin)=0d0
       return
    end if

    Ltot=sqrt(h%L_star%x**2+h%L_star%y**2+h%L_star%z**2)
    Lx=h%L_star%x/Ltot
    Ly=h%L_star%y/Ltot
    Lz=h%L_star%z/Ltot
    xc = h%p_star%x
    yc = h%p_star%y
    zc = h%p_star%z
    rmax = min(h%r_star, h%datas%rvir)  ! instead of 3d0*1d-3
#else
    Ltot=sqrt(h%L%x**2+h%L%y**2+h%L%z**2)
    Lx=h%L%x/Ltot
    Ly=h%L%y/Ltot
    Lz=h%L%z/Ltot
    xc = h%p%x
    yc = h%p%y
    zc = h%p%z
    rmax = h%r  ! instead of 3d0*1d-3
#endif

    dx=rmax/dble(h%nbin)
    indexp = first_part(h%my_number)
    i=0
    m_tot=0d0
    do while (indexp /= -1)
       i=i+1
       if (family(indexp) == FAM_STAR) then
          dxx=pos(indexp,1)-xc
          dyy=pos(indexp,2)-yc
          dzz=pos(indexp,3)-zc
          dzcyl=dxx*Lx+dyy*Ly+dzz*Lz
          ! Compute centre coord. of the plane where the particle lies
          ! and perp. to the disc plane
          planecen_x=xc+dzcyl*Lx
          planecen_y=yc+dzcyl*Ly
          planecen_z=zc+dzcyl*Lz
          ! Compute unit radial vector from this centre to the particle
          rx=pos(indexp,1)-planecen_x
          ry=pos(indexp,2)-planecen_y
          rz=pos(indexp,3)-planecen_z
          modr=sqrt(rx**2d0+ry**2d0+rz**2d0)

          j=int(modr/dx)+1
          if (j .le. h%nbin)then
             if (allocated(mass)) then
                m_tot(j)=m_tot(j)+mass(indexp)
             else
                m_tot(j)=m_tot(j)+massp
             endif
          endif
       end if
       indexp   = linked_list(indexp)
    end do

    do j=1,h%nbin
       h%rr(j)=(j-0.5)*dx
    enddo
    do j=1,h%nbin
       dv(j)=acos(-1d0)*( (h%rr(j)+0.5d0*dx)**2d0 - (h%rr(j)-0.5d0*dx)**2d0 )
    enddo
    h%rr=h%rr*1d3
    h%rho=m_tot*1d11/(dv*1d6**2d0)

    m_in = 0d0
    jeff = -1
    do j=1,h%nbin
       m_in = m_in + m_tot(j)
       if ((m_in > sum(m_tot)*0.5d0) .and. jeff<0) jeff = j
    end do
    if (jeff<0) jeff=h%nbin
    h%datas%Reff = (jeff - 0.5d0)*dx


    return

  end subroutine compute_stellar_profile

!***********************************************************************
  subroutine compute_stellar_kinematics(h)

    ! compute bulge/disc decomposition, velocity dispersion, etc

    implicit none

    type (halo)              :: h
    integer(i8b)             :: indexp
    real(kind=8)             :: xc, yc, zc
    real(kind=8)             :: vcz,vcr,vct,vbulge,vnorm
    real(kind=8)             :: Lx,Ly,Lz,Ltot,mbulge,dxx,dyy,dzz,dzcyl
    real(kind=8)             :: planecen_x,planecen_y,planecen_z,rx,ry,rz,modr
    real(kind=8)             :: thetax,thetay,thetaz
    integer(i8b)             :: i,nb_of_stars
    real(kind=8)             :: mtot,vrmean,vzmean,vtmean
    real(kind=8)             :: dvrmean,dvzmean,dvtmean
    real(kind=8)             :: vr2mean,vz2mean,vt2mean
    real(kind=8)             :: sigma2
    real(kind=8)             :: dvr2mean,dvz2mean,dvt2mean
    real(kind=8)             :: vzmean_disc,vrmean_disc,vtmean_disc
    real(kind=8)             :: dvzmean_disc,dvrmean_disc,dvtmean_disc
    real(kind=8)             :: dvz2mean_disc,dvr2mean_disc,dvt2mean_disc
    real(kind=8),allocatable :: vczp(:),vcrp(:),vctp(:),mtmp(:)

#ifdef ALLPARTS
    Ltot=sqrt(h%L_star%x**2+h%L_star%y**2+h%L_star%z**2)
    Lx=h%L_star%x/Ltot
    Ly=h%L_star%y/Ltot
    Lz=h%L_star%z/Ltot
    xc = h%p_star%x
    yc = h%p_star%y
    zc = h%p_star%z
    nb_of_stars = h%n_star

    if (nb_of_stars < nMembers) then
       ! If not enough stars, just bypass this
       h%sigma_bulge = 0d0
       h%m_bulge     = 0d0
       h%sigma1D     = 0d0
       h%Vsigma      = 0d0
       h%sigma1D_disc= 0d0
       h%Vsigma_disc = 0d0
       return
    end if
#else
    Ltot=sqrt(h%L%x**2+h%L%y**2+h%L%z**2)
    Lx=h%L%x/Ltot
    Ly=h%L%y/Ltot
    Lz=h%L%z/Ltot
    xc = h%p%x
    yc = h%p%y
    zc = h%p%z
    nb_of_stars = nb_of_parts(h%my_number)
#endif

    ! vtmean: average tangential velocity
    ! vzmean: average vertical velocity
    ! vrmean: average radial velocity
    ! dvtmean: average tangential sigma
    ! dvzmean: average vertical sigma
    ! dvrmean: average radial sigma
    ! sigma1D: sqrt((dvzmean**2+dvrmean**2+dvtmean)/3d0)
    ! V/sigma: vtmean/sigma1D

    allocate(vczp(1:nb_of_stars),vcrp(1:nb_of_stars),vctp(1:nb_of_stars),mtmp(1:nb_of_stars))

    indexp = first_part(h%my_number)
    i=0
    mtot=0d0
    mbulge=0d0
    vbulge=0d0
    vrmean=0d0;vr2mean=0d0
    vzmean=0d0;vz2mean=0d0
    vtmean=0d0;vt2mean=0d0
    vrmean_disc=0d0;vzmean_disc=0d0;vtmean_disc=0d0
    do while (indexp /= -1)
       if (family(indexp) == FAM_STAR) then
          i=i+1
          dxx=pos(indexp,1)-xc
          dyy=pos(indexp,2)-yc
          dzz=pos(indexp,3)-zc
          dzcyl=dxx*Lx+dyy*Ly+dzz*Lz
          ! Compute centre coord. of the plane where the particle lies
          ! and perp. to the disc plane
          planecen_x=xc+dzcyl*Lx
          planecen_y=yc+dzcyl*Ly
          planecen_z=zc+dzcyl*Lz
          ! Compute unit radial vector from this centre to the particle
          rx=pos(indexp,1)-planecen_x
          ry=pos(indexp,2)-planecen_y
          rz=pos(indexp,3)-planecen_z
          modr=sqrt(rx**2d0+ry**2d0+rz**2d0)

          ! Check if the vector is defined
          if (modr .ne. 0d0)then
             rx=rx/modr
             ry=ry/modr
             rz=rz/modr
          endif

          ! Compute unit vector perp. to this vector and in the same plane
          thetax=Ly*rz-Lz*ry
          thetay=Lz*rx-Lx*rz
          thetaz=Lx*ry-Ly*rx

          ! Compute velocity components in this cylindrical coord.
          vcz = vel(indexp,1)*Lx     + vel(indexp,2)*Ly     + vel(indexp,3)*Lz
          vcr = vel(indexp,1)*rx     + vel(indexp,2)*ry     + vel(indexp,3)*rz
          vct = vel(indexp,1)*thetax + vel(indexp,2)*thetay + vel(indexp,3)*thetaz
          vczp(i)=vcz
          vcrp(i)=vcr
          vctp(i)=vct
          mtmp(i)=mass(indexp)
          vrmean  =vrmean  +vcr     *mass(indexp)
          vr2mean =vr2mean +vcr*vcr *mass(indexp)
          vzmean  =vzmean  +vcz     *mass(indexp)
          vz2mean =vz2mean +vcz*vcz *mass(indexp)
          vtmean  =vtmean  +vct     *mass(indexp)
          vt2mean =vt2mean +vct*vct *mass(indexp)
          mtot    =mtot    +         mass(indexp)
          if (abs(vct) .lt. sqrt(vcr**2+vcz**2)) then
             ! Star is in the (kinematical) bulge
             vnorm = real(sqrt(vcz**2+vcr**2+vct**2),8)
             mbulge = mbulge + mass(indexp)
             vbulge = vbulge + mass(indexp)*vnorm
          else
             ! Star is in the disc
             vrmean_disc = vrmean_disc + vcr*mass(indexp)
             vzmean_disc = vzmean_disc + vcz*mass(indexp)
             vtmean_disc = vtmean_disc + vct*mass(indexp)
          end if
       end if
       indexp   = linked_list(indexp)
    end do
    vrmean=vrmean/mtot;vr2mean=vr2mean/mtot
    vzmean=vzmean/mtot;vz2mean=vz2mean/mtot
    vtmean=vtmean/mtot;vt2mean=vt2mean/mtot
    ! Same for disc, if it exists
    if (mtot>mbulge) then
       vrmean_disc = vrmean_disc/(mtot-mbulge)
       vzmean_disc = vzmean_disc/(mtot-mbulge)
       vtmean_disc = vtmean_disc/(mtot-mbulge)
    end if

    dvr2mean=0d0;dvz2mean=0d0;dvt2mean=0d0
    dvr2mean_disc=0d0;dvz2mean_disc=0d0;dvt2mean_disc=0d0
    do i=1,nb_of_stars
       dvr2mean=dvr2mean+(vcrp(i)-vrmean)**2*mtmp(i)
       dvz2mean=dvz2mean+(vczp(i)-vzmean)**2*mtmp(i)
       dvt2mean=dvt2mean+(vctp(i)-vtmean)**2*mtmp(i)
       if (abs(vctp(i)) .lt. sqrt(vcrp(i)**2+vczp(i)**2)) then ! Same for disc
          dvr2mean_disc=dvr2mean_disc+(vcrp(i)-vrmean_disc)**2*mtmp(i)
          dvz2mean_disc=dvz2mean_disc+(vczp(i)-vzmean_disc)**2*mtmp(i)
          dvt2mean_disc=dvt2mean_disc+(vctp(i)-vtmean_disc)**2*mtmp(i)
       end if
    enddo
    dvrmean=sqrt(dvr2mean/mtot)
    dvzmean=sqrt(dvz2mean/mtot)
    dvtmean=sqrt(dvt2mean/mtot)
    ! Same for disc, if it exists
    if (mtot>mbulge) then
       dvrmean_disc = sqrt(dvr2mean_disc/(mtot-mbulge))
       dvzmean_disc = sqrt(dvz2mean_disc/(mtot-mbulge))
       dvtmean_disc = sqrt(dvt2mean_disc/(mtot-mbulge))
    end if

    ! Compute bulge 3D velocity dispersion
    sigma2 = 0d0
    do i=1,nb_of_stars
       if (abs(vctp(i)) .lt. sqrt(vcrp(i)**2+vczp(i)**2)) then
          vnorm = real(sqrt(vczp(i)**2+vcrp(i)**2+vctp(i)**2),8)
          sigma2 = sigma2 + mtmp(i)*(vnorm-vbulge)**2
       end if
    end do


    h%sigma_bulge = real(sqrt(sigma2/mbulge),8)
    h%m_bulge     = real(mbulge,8)
    h%sigma1D     = sqrt((dvzmean**2+dvrmean**2+dvtmean**2)/3d0)
    h%Vsigma      = vtmean/h%sigma1D
    h%sigma1D_disc= sqrt((dvzmean_disc**2+dvrmean_disc**2+dvtmean_disc**2)/3d0)
    h%Vsigma_disc = vtmean_disc/h%sigma1D_disc



    deallocate(vcrp, vczp, vctp, mtmp)

    return

  end subroutine compute_stellar_kinematics

!***********************************************************************
  subroutine compute_stellar_metallicity(h)

  ! compute mass-weighted metallicity of the stars

    implicit none

    type (halo)        :: h
    integer(i8b)       :: indexp
    real(kind=8)       :: met, mstar

    met = 0d0
    mstar = 0d0

    indexp = first_part(h%my_number)
    do while (indexp /= -1)
       if (family(indexp) == FAM_STAR) then
          met = met + real(met_st(indexp),8)*real(mass(indexp),8)
          mstar = mstar + real(mass(indexp),8)
       end if
       indexp   = linked_list(indexp)
    end do

    if (mstar>0) then
       h%metallicity = met / mstar
    else
       h%metallicity = 0d0
    end if

    return

  end subroutine compute_stellar_metallicity

!***********************************************************************
  subroutine compute_stellar_age(h)

    ! compute mass-weighted age of the stars

    implicit none

    type (halo)        :: h
    integer(i8b)       :: indexp
    real(kind=8)       :: age, mstar, agepart
    real(kind=8)       :: sfr10, sfr100, sfr1000

    age     = 0d0
    mstar   = 0d0
    sfr10   = 0d0
    sfr100  = 0d0
    sfr1000 = 0d0

    indexp = first_part(h%my_number)
    do while (indexp /= -1)
       if (family(indexp) == FAM_STAR) then
          if (proper_time) then
             agepart = lbtime - ct_proptime2time(age_st(indexp),H_f)
          else
             agepart = lbtime - ct_conftime2time(age_st(indexp))
          end if
          age = age + real(agepart,8)*real(mass(indexp),8)
          mstar = mstar + real(mass(indexp),8)
          if (agepart <= 10d6) sfr10 = sfr10 + real(mass(indexp),8)
          if (agepart <= 100d6) sfr100 = sfr100 + real(mass(indexp),8)
          if (agepart <= 1000d6) sfr1000 = sfr1000 + real(mass(indexp),8)
       end if
       indexp   = linked_list(indexp)
    end do

    if (mstar>0) then
       h%age = age / mstar
       h%sfr10 = sfr10 / 10d6
       h%sfr100 = sfr100 / 100d6
       h%sfr1000 = sfr1000 / 1000d6
    else
       h%age = -1d0
       h%sfr10 = 0d0
       h%sfr100 = 0d0
       h%sfr1000 = 0d0
    end if

    return

  end subroutine compute_stellar_age
#endif
!***********************************************************************
  function interact(i,j)

    implicit none

    integer(i8b)    :: i,j,ifirst
    real(kind=8)    :: dist2ij
    real(kind=8)    :: interact,rinveff,r3inveff
    real(kind=8)    :: epstmp,massp2,lbox2
    type (vector)   :: dr

    save ifirst,epstmp,massp2
    data ifirst /1/

    if (ifirst == 1) then
       ! epstmp is mean interparticular distance / 20.0 (i.e. smoothing length)
       ! true for tree code but RAMSES and ENZO ?
       epstmp = (massp/mboxp)**(1./3.) / 20.0
       massp2 = massp*massp
       ifirst = 0
    end if

    Lbox2   = Lbox_pt**2
    dr%x    = pos(j,1) - pos(i,1)
    dr%y    = pos(j,2) - pos(i,2)
    dr%z    = pos(j,3) - pos(i,3)
    call correct_for_periodicity(dr)
    dist2ij = (dr%x**2) + (dr%y**2) + (dr%z**2)
    dist2ij = dist2ij / Lbox2

    if (allocated(mass)) then
       if (allocated(epsvect)) then
          call softgrav(epsvect(i),epsvect(j),dist2ij,rinveff,r3inveff)
          interact = -mass(i) * mass(j) * rinveff
       else
          ! do not correct for softening --> have to change that
          interact =-mass(i)*mass(j)/sqrt(dist2ij)
       endif
    else
       call softgrav(epstmp,epstmp,dist2ij,rinveff,r3inveff)
       interact = -massp2 * rinveff
    endif

    return

  end function interact

!***********************************************************************
  subroutine softgrav(epsp,epsi,drdotdr,rinveff,r3inveff)

  ! subroutine to compute the effective distance between particles of
  ! smoothing lengths, epsp and epsi, given their real distance**2,
  ! drdotdr, in order to get the smoothed values of the potential and
  ! acceleration phi and acc (in GRAVSUM). Calculations are for an
  ! harmonic smoothing or a cubic spline smoothing kernel.
  ! For the spline smoothing, phsmooth and acsmooth must have
  ! been initialized by initgsoft.

    implicit none

    integer(kind=4) :: smindex
    real(kind=8)    :: epsp,epsi,drdotdr,sdrdotdr,rinveff,r3inveff,drdeldrg
    real(kind=8)    :: drsm,phsm,epseff,epsinv,dr,tiny,one

    tiny = 1.d-19
    one  = 1.0

    if (gravsoft == 'harmonic') then

       drdotdr  = drdotdr+tiny
       rinveff  = 1.0/sqrt(drdotdr)
       r3inveff = rinveff*rinveff*rinveff
       epseff   = 0.5*(epsp+epsi)
       if (drdotdr .lt. epseff*epseff) then
          epsinv   = 1.0/epseff
          r3inveff = epsinv*epsinv*epsinv
          rinveff  = 1.5*epsinv-0.5*drdotdr*r3inveff
       endif

    else if (gravsoft == 'cubsplin') then

       dr       = epsp+epsi
       drdotdr  = drdotdr+tiny*0.25*dr**2
       sdrdotdr = sqrt(drdotdr)
       ! rinveff  = 1.0/sdrdotdr
       ! we never use r3inveff at the mo, disable it.
       ! r3inveff = rinveff*rinveff*rinveff
       drdeldrg = sdrdotdr*ninterp/dr
       smindex  = drdeldrg
       ! if (ninterp .lt. smindex) smindex = ninterp
       ! if (one .lt. drdeldrg-smindex) then
       !    drsm = one
       ! else
       !    drsm = drdeldrg-smindex
       ! endif
       ! phsm = (1.-drsm)*phsmooth(smindex) + drsm*phsmooth(1+smindex)
       if (smindex > ninterp) then
          phsm = phsmooth(1+ninterp)
       else
          if (one < drdeldrg-smindex) then
             phsm = phsmooth(1+smindex)
          else
             drsm = drdeldrg-smindex
             phsm = (1.-drsm)*phsmooth(smindex)+drsm*phsmooth(1+smindex)
          endif
       end if
       ! rinveff = phsm*rinveff
       rinveff = phsm/sdrdotdr
       ! NB: r3inveff is used to compute the acceleration and necessitates
       !     the definition of accsmooth which is not available here.
       !     For the treecode, the relation should be:
       !     r3inveff=accsmooth*r3inveff (to be checked)

    endif

    return

  end subroutine softgrav

!***********************************************************************
  subroutine tab_props_inside(h,nr,tabm2,tabk2,tabp2,v,amax,bmax,cmax,fam)

  ! returns the cumulative mass contained in concentric ellipsoids centered on the center of the
  ! halo (cdm or mbp)

    implicit none

    integer(i8b)           :: nr,indexp,i,louped_parts
    integer(kind=4)        :: num_h, i_ell
    real(kind=8)           :: amax,bmax,cmax
    real(kind=8)           :: v(3,3)
    real(kind=8)           :: tabm2(0:nr-1),tabk2(0:nr-1),tabp2(0:nr-1) !! double
    real(kind=8)           :: srm,srk                                   !! double
    real(kind=8)           :: rmax,dra,drb,drc
    real(kind=8)           :: r_ell,v2,rf
    real(kind=8),parameter :: epsilon = 1.d-2
    type (vector)          :: posp,vt
    type (halo)            :: h
    integer(kind=1)        :: fam
    logical                :: sel=.false.
    real(kind=8)           :: xc, yc, zc, vxc, vyc, vzc, sha, shb, shc


    ! rescale to get ellipsoid  concentric to principal ellipsoid
    ! which contains all the particles of the halo
    rmax = 0.0
    indexp = first_part(h%my_number)
#ifdef ALLPARTS
    if (fam == FAM_STAR) then
       xc = h%p_star%x
       yc = h%p_star%y
       zc = h%p_star%z
       vxc = h%v_star%x
       vyc = h%v_star%y
       vzc = h%v_star%z
       sha = h%sh_star%a
       shb = h%sh_star%b
       shc = h%sh_star%c
    elseif (fam == FAM_DM) then
       xc = h%p_DM%x
       yc = h%p_DM%y
       zc = h%p_DM%z
       vxc = h%v_DM%x
       vyc = h%v_DM%y
       vzc = h%v_DM%z
       sha = h%sh_DM%a
       shb = h%sh_DM%b
       shc = h%sh_DM%c
    else
#endif
       xc = h%p%x
       yc = h%p%y
       zc = h%p%z
       vxc = h%v%x
       vyc = h%v%y
       vzc = h%v%z
       sha = h%sh%a
       shb = h%sh%b
       shc = h%sh%c
#ifdef ALLPARTS
    end if
#endif

    do while(indexp.gt.0)
       if (fam == FAM_UNDEF) then
          ! No family defined
          sel = .true.
       else
          ! Select only particles with the right type
          sel = (family(indexp) == fam)
       end if
       if (sel) then
          posp%x = pos(indexp,1) - xc
          posp%y = pos(indexp,2) - yc
          posp%z = pos(indexp,3) - zc
          call correct_for_periodicity(posp)
          ! project position vector along the principal ellipsoid axis
          dra    = posp%x*v(1,1)+posp%y*v(2,1)+posp%z*v(3,1)
          drb    = posp%x*v(1,2)+posp%y*v(2,2)+posp%z*v(3,2)
          drc    = posp%x*v(1,3)+posp%y*v(2,3)+posp%z*v(3,3)
          r_ell  = sqrt((dra / sha)**2 + (drb / shb)**2 + (drc / shc)**2)
          rmax   = max(rmax,r_ell)
       end if
       indexp = linked_list(indexp)
    end do

    amax = rmax * sha * (1.0 + epsilon)
    bmax = rmax * shb * (1.0 + epsilon)
    cmax = rmax * shc * (1.0 + epsilon)

    ! initialize loop quantities
    tabm2        = 0d0
    tabk2        = 0d0
    louped_parts = 0
    num_h        = h%my_number
    indexp       = first_part(num_h)

    do while (indexp /= -1)

       if (fam == FAM_UNDEF) then
          ! No family defined
          sel = .true.
       else
          ! Select only particles with the right type
          sel = (family(indexp) == fam)
       end if
       if (sel) then
          posp%x = pos(indexp,1) - xc
          posp%y = pos(indexp,2) - yc
          posp%z = pos(indexp,3) - zc
          call correct_for_periodicity(posp)
          ! compute velocities in the halo frame adding in the Hubble flow
          vt%x   = vel(indexp,1) - vxc + posp%x * Hub_pt
          vt%y   = vel(indexp,2) - vyc + posp%y * Hub_pt
          vt%z   = vel(indexp,3) - vzc + posp%z * Hub_pt
          v2     = vt%x**2 + vt%y**2 + vt%z**2
          ! project position vector along the principal ellipsoid axis
          dra    = posp%x*v(1,1)+posp%y*v(2,1)+posp%z*v(3,1)
          drb    = posp%x*v(1,2)+posp%y*v(2,2)+posp%z*v(3,2)
          drc    = posp%x*v(1,3)+posp%y*v(2,3)+posp%z*v(3,3)
          ! biggest ellipsoid is divided in nr concentric ellipsoid shells: we
          ! calculate below the ellipsoid bin in which each particle falls and fill up the
          ! mass and energy tables accordingly
          ! NB: if by chance the most distant particle from the center is ON the shortest
          !     axis, then r_ell is equal to 1-epsilon (one of the terms below is one and the others 0)
          !     otherwise r_ell is between 0 and 1-epsilon and so we just multiply it by nr to
          !     find the ellipsoid shell containing the particle.
          r_ell  = sqrt((dra / amax)**2 + (drb / bmax)**2 + (drc / cmax)**2)
          i_ell  = int(r_ell*nr)
          if (i_ell .lt. nr) then
             if (allocated(mass)) then
                tabm2(i_ell) = tabm2(i_ell)+real(mass(indexp),8)
                tabk2(i_ell) = tabk2(i_ell)+0.5*real(mass(indexp),8)*v2
             else
                tabm2(i_ell) = tabm2(i_ell)+real(massp,8)
                tabk2(i_ell) = tabk2(i_ell)+0.5*real(massp,8)*v2
             endif
          else
             louped_parts = louped_parts + 1
          endif
       end if
       indexp = linked_list(indexp)

    end do

    if (louped_parts .gt. 0) then
       write(errunit,*) ''
       write(errunit,*) '> Problem in tab_props_inside : missed ',louped_parts,' particles'
       write(errunit,*) ''
       stop
    end if

    srm = tabm2(0)
    srk = tabk2(0)
    do i = 1,nr-1
       srm      = srm+tabm2(i)
       srk      = srk+tabk2(i)
       tabm2(i) = srm
       tabk2(i) = srk
       ! approximation based on appendix B of paper GALICS 1:
       ! better than 10-15 % accuracy on average
       tabp2(i) = -0.3 * gravconst * tabm2(i)**2 * rf(sha**2,shb**2,shc**2)
    end do
    ! correct potential energy estimate for small halos which are calculated by direct summation
    if (h%ep /= tabp2(nr-1)) tabp2 = tabp2/tabp2(nr-1)*h%ep

    return

  end subroutine tab_props_inside

!***********************************************************************
#ifdef ANG_MOM_OF_R
  subroutine det_vir_props(h,v,amax,bmax,cmax)
#else
  subroutine det_vir_props(h,v)
#endif

    ! computes the virial properties (radius, mass) of a halo

    implicit none

    integer(kind=4)           :: i,ii,ivmax
    ! ttab = 1000 bins for virial radius precision better than 1% of halo size
    integer(kind=i8b),parameter :: ttab = 1000
    real(kind=8)              :: rvir,mvir,kvir,pvir
    real(kind=8)              :: v(3,3)
    real(kind=8)              :: amax,bmax,cmax,avir,bvir,cvir
    real(kind=8)              :: tab_mass(0:ttab-1),tab_ekin(0:ttab-1),tab_epot(0:ttab-1)  !! double
    real(kind=8)              :: virth,virth_old,volmin
    real(kind=8)              :: vmax,rvmax,m200
    type (halo)               :: h

    ! compute properties inside ttab concentric principal ellipsoids centered on center of halo
    call tab_props_inside(h,ttab,tab_mass,tab_ekin,tab_epot,v,amax,bmax,cmax, FAM_UNDEF)

    ! find the outermost ellipsoid bin where virial theorem is either satisfied better than 20 %
    ! or satisfied best ...
    mvir      = tab_mass(ttab-1)
    kvir      = tab_ekin(ttab-1)
    pvir      = tab_epot(ttab-1)
    ! initialize rvir to be the geometric average of the axis radii of the outermost ellipsoid shell
    ! which contains at least one particle
    do i = ttab-1,1,-1
       if (tab_mass(i-1) < tab_mass(i)) exit
    enddo
    avir      = real(i,8)/real(ttab-1,8)*amax
    bvir      = real(i,8)/real(ttab-1,8)*bmax
    cvir      = real(i,8)/real(ttab-1,8)*cmax
    rvir      = (avir*bvir*cvir)**(1./3.)
    ! assume initial departure from virialization is 100 %
    virth_old = 1.0
    virth     = 1.0
    do i = ttab-1,0,-1
       ! if region is unbound, it cannot be virialized in the same time
       if (tab_ekin(i)+tab_epot(i) >= 0.0) cycle
       ! region is bound so compute relative virialization |2*K+P|/|K+P|
       virth = abs((2.0*tab_ekin(i)+tab_epot(i))/(tab_ekin(i)+tab_epot(i)))
       ! if region is better virialized then update virial quantities
       if (virth < virth_old) then
          mvir      = tab_mass(i)
          ! take the min here bc initialization throws away all the empty outer shells
          avir      = min(avir,real(i,8)/real(ttab-1,8)*amax)
          bvir      = min(bvir,real(i,8)/real(ttab-1,8)*bmax)
          cvir      = min(cvir,real(i,8)/real(ttab-1,8)*cmax)
          rvir      = (avir*bvir*cvir)**(1./3.)
          kvir      = tab_ekin(i)
          pvir      = tab_epot(i)
          virth_old = virth
          ! if virial theorem holds with better than 20 % accuracy, exit do loop
          if (virth <= 0.20) exit
       endif
    end do
    !write(*,*)'Mvir=',mvir

    ! for small halos it may happen that virial theorem is not enforced to within 15 %
    ! bc the halo is not fully virialized yet or that it is valid by fluke (right combination
    ! of potential and kinetic energy) ... so ....
    ! 1/ in the latter case, we further check that the halo density is high enough
    ! (order of what predicted by the spherical top-hat model, vir_overdens)
    ! for us to believe in the measurement of the virial theorem.
    ! 2/ in the former case, as an inner region measurement of the virialization would be too noisy for
    ! lack of mass resolution (not enough particles) we only use the overdensity criterion.
    ! NB: this criterion is similar to NFW but with vir_overdens * average density and NOT 200. * critical density
    !     which does not makes sense when Omega_matter(z) /= 1 (not critical) and we take ellipsoids
    !     NOT spheres bc they yield more accurate volumes and average halo densities
    ! average density of the universe at current timestep (in 10^11 M_sun / Mpc^3)

    ! volume of the smallest concentric ellipsoid
    volmin   = 4./3.*pi*(amax/real(ttab-1,8))*(bmax/real(ttab-1,8))*(cmax/real(ttab-1,8))
    if (virth > 0.20 .or. mvir < vir_overdens*rho_mean*4./3.*pi*avir*bvir*cvir) then
       do ii = ttab-1,1,-1
          ! assume that the virial mass and radii are obtained when the density inside the ellipsoid
          ! is greater than vir_overdens * rho_mean AND there is at least one particle inside the outermost
          ! ellipsoid shell
          mvir = vir_overdens * rho_mean * volmin * real(ii,8)**3
          !mvir = 200d0 * 3d0*Hub_pt**2/8d0/acos(-1d0)/gravconst * volmin * real(ii,4)**3
          if (tab_mass(ii) >= mvir .and. tab_mass(ii-1) < tab_mass(ttab-1)) exit
       enddo
       mvir   = tab_mass(ii)
       kvir   = tab_ekin(ii)
       pvir   = tab_epot(ii)
       avir   = real(ii,8)/real(ttab-1,8)*amax
       bvir   = real(ii,8)/real(ttab-1,8)*bmax
       cvir   = real(ii,8)/real(ttab-1,8)*cmax
       rvir   = (avir*bvir*cvir)**(1./3.)
    endif
    !write(*,*)ii,ttab-1
    !write(*,*)'rhoc=',3d0*Hub_pt**2/8d0/acos(-1d0)/gravconst,'1d11 Msun/Mpc^3'
    !write(*,*)tab_mass(ttab-1)
    !write(*,*)'M200=',mvir
    !write(*,*)'r200=',rvir

    ! check if virialization conditions were met --> if not set relevant quantities to zero: this is a
    ! non-virialized halo ....
    if (mvir > 0.0 .and. rvir > 0.0) then
       ! it may happen (although it is very rare bc vector linking center of halo to most distant
       ! particle has to be roughly parallel to minor axis of the halo in such a case) that the virial
       ! radius estimated from the geometric mean of the 3 principal axis is slightly larger than the
       ! distance of the most distant particle to the center of the halo (what we call r)
       ! when this happens we set the virial radius to be r
       h%datas%rvir         = min(rvir,h%r)    ! in Mpc
       h%datas%mvir         = mvir             ! in 10^11 M_sun
       ! circular velocity at r_vir in km/s
       h%datas%cvel         = sqrt(gravconst*h%datas%mvir/h%datas%rvir)
       ! temperature at r_vir in K
       h%datas%tvir         = 35.9*gravconst*h%datas%mvir/h%datas%rvir
       ! compute halo density profile within the virialized region
       call compute_halo_profile(h)

       ! Non-inclusive R200 and M200
       ! First radius (from outside) for which 200*rhom*volume
       ! This gives R200
       do ii = ttab-1,1,-1
          m200 = 200.d0 * rho_mean * volmin * real(ii,8)**3
          if (tab_mass(ii) >= m200 .and. tab_mass(ii-1) < tab_mass(ttab-1)) exit
       enddo
       h%m200 = tab_mass(ii)
       h%r200  = real(ii,8)/real(ttab-1,8) * (amax*bmax*cmax)**(1d0/3d0)


       ! We now estimate the radius at which sqrt(GM(<r)/r) is maximal: this would be Rmax and Vmax
       !!! TEST: only DM
#ifdef ALLPARTS
       if (h%n_DM >= nMembers) then
          ! compute properties inside ttab concentric principal ellipsoids centered on center of halo
          call tab_props_inside(h,ttab,tab_mass,tab_ekin,tab_epot,v,amax,bmax,cmax, FAM_DM)
#endif
          vmax = 0d0
          ivmax = 0
          do ii = ttab-1,1,-1
             rvmax = real(ii,8)/real(ttab-1,8) * (amax*bmax*cmax)**(1d0/3d0)
             ! vmax to be compared to sqrt(gravconst*tab_mass(ii)/rvmax)
             if ((vmax .lt. sqrt(gravconst*tab_mass(ii)/rvmax)) .and. (rvmax.le.h%datas%rvir)) then
                vmax = sqrt(gravconst*tab_mass(ii)/rvmax)
                ivmax = ii
             end if
          enddo
          rvmax   = real(ivmax,8)/real(ttab-1,8) * (amax*bmax*cmax)**(1d0/3d0)
          h%rvmax = rvmax
          h%vmax  = vmax

          ! NFW concentration from Prada+2012
          ! vmax/cvel ~ sqrt(.216*c/f(c))
          call compute_concentration(h)
#ifdef ALLPARTS
       else
          h%rvmax = 0d0
          h%vmax  = 0d0
          h%cNFW  = 0d0
       end if
#endif

    else
       write(*,*) 'halo bugged',h%my_number,mvir,rvir
       write(*,*)h%p%x,h%p%y,h%p%z
       stop
    endif

    return

  end subroutine det_vir_props

!***********************************************************************
  subroutine compute_halo_profile(h)

    implicit none

    type (halo)     :: h

    if (profile == 'TSIS') then
       ! for the singular isothermal sphere the profile is defined @ rvir for it is singular at r=0.
       h%halo_profile%rho_0 = h%datas%mvir / (4.0 * pi * h%datas%rvir**3)
       h%halo_profile%r_c   = h%datas%rvir
    else
       write(errunit,*) 'Other profiles than TSIS not yet fully implemented'
       stop
    endif

    return

  end subroutine compute_halo_profile

!***********************************************************************
  subroutine compute_concentration(h)

    implicit none

    type (halo)     :: h
    integer(kind=4) :: i, iNFW
    real(kind=8)    :: v
    real(kind=8),parameter :: cmin=2.2d0, cmax=1000,dlogc=log10(cmax/cmin)/1000d0
    real(kind=8),parameter :: logc(1:1000) = (/ (log10(cmin) + (i-1)*dlogc, i=1,1000) /)
    real(kind=8),parameter :: cc(1:1000) = 10d0**logc
    real(kind=8),parameter :: cfc(1:1000) = cc/(log(1.d0+cc) - cc/(1.d0+cc))

    v = max(h%datas%cvel, h%vmax)
    call locate(cfc,1000,(v/h%datas%cvel)**2d0 / 0.216d0 ,iNFW)
    h%cNFW = cc(iNFW)

    return

  end subroutine compute_concentration

!***********************************************************************
  subroutine change_units

    ! subroutine which goes from positions in code units to physical (non comoving)
    ! Mpc and from velocities in code units to peculiar (no Hubble flow) velocities
    ! in km/s. masses are also changed from code units to 10^11 M_sun

    implicit none

    pos   = pos * Lbox_pt
    if (type == 'SN') vel = vel*Hub_pt*Lbox_pt
    massp = massp * mboxp
    if (allocated(mass)) mass  = mass * mboxp
    if(verbose) then
       write(errunit,*) '> xmin,xmax:',minval(pos(:,1)),maxval(pos(:,1))
       write(errunit,*) '> ymin,ymax:',minval(pos(:,2)),maxval(pos(:,2))
       write(errunit,*) '> zmin,zmax:',minval(pos(:,3)),maxval(pos(:,3))
       write(errunit,*) '> particle mass:', massp*1.e11
    end if

    return

  end subroutine change_units

!***********************************************************************
  function temps(x,omm,oml)

    implicit none

    real(kind=8) :: temps,x,omm,oml

    temps = 1./sqrt(omm*x**5+(1.0-omm-oml)*x**4+oml*x**2)

    return

  end function temps

!***********************************************************************
  function temps_turn_around(x,omm,oml)

    implicit none

    real(kind=8) :: temps,temps_turn_around,x,omm,oml

    temps = omm*x**5+(-omm-oml)*x**4+oml*x**2
    if (temps > 0.0) then
       temps_turn_around = 1./sqrt(temps)
    else
       temps = 0.0
    endif

    return

  end function temps_turn_around

!***********************************************************************
  function dtda_turn_around(x,omm,oml)

    ! find dt/da given a, from the Friedmann Equation.
    ! here, time "t" is understood to be in units of the inverse Hubble constant
    ! (i.e. "t" = H0*t)
    ! define the curvature so that turnaround occurs for a = 1.0
    ! Note that omega_matter+omega_lambda+omega_curvature sum to ZERO, not one, in that case
    ! definitions for parameters are as in Peebles 1993, eqn (5.53)

    implicit none

    real(kind=8) :: x,omm,oml,temp,dtda_turn_around

    temp = omm*x + oml/x**2 - (omm+oml)
    if (temp .gt. 0.0) then
       dtda_turn_around = sqrt(1./temp)
    else
       dtda_turn_around = 0.0
    end if

    return

  end function dtda_turn_around

!***********************************************************************
  subroutine det_vir(h)

  ! determine virial properties of the halos, energies and profiles

    implicit none

    type (halo)     :: h
    real(kind=8)    :: v(3,3)
    real(kind=8)    :: x0,y0,z0,r0
#ifdef ANG_MOM_OF_R
    real(kind=8)    :: amax,bmax,cmax
#endif
    integer(kind=4) :: nx

    if(method.ne."FOF".and.DPMMC) then
       x0=dble(h%p%x);y0=dble(h%p%y);z0=dble(h%p%z);r0=dble(h%r)
       nx=9 ! The initial mesh have a size of 2 times the maximum radius (~virial radius)
       call det_halo_center_multiscale(h,x0,y0,z0,r0,nx,FAM_UNDEF)
#ifdef ALLPARTS
       if (h%n_star >= nMembers) then
          call det_halo_center_multiscale(h,h%p_star%x,h%p_star%y,h%p_star%z,h%r_star,nx,FAM_STAR)
!!$       else
!!$          h%p_star%x=h%p%x
!!$          h%p_star%y=h%p%y
!!$          h%p_star%z=h%p%z
       end if
       if (h%n_DM >= nMembers) then
          call det_halo_center_multiscale(h,h%p_DM%x,h%p_DM%y,h%p_DM%z,h%r_DM,nx,FAM_DM)
!!$       else
!!$          h%p_DM%x=h%p%x
!!$          h%p_DM%y=h%p%y
!!$          h%p_DM%z=h%p%z
       end if
#endif
    endif
    if(method.ne."FOF".and.SC) then
       x0=dble(h%p%x);y0=dble(h%p%y);z0=dble(h%p%z);r0=dble(h%r)*0.2d0
       call det_halo_center_sphere(h,x0,y0,z0,r0,FAM_UNDEF)
#ifdef ALLPARTS
       if (h%n_star >= nMembers) then
          call det_halo_center_sphere(h,h%p_star%x,h%p_star%y,h%p_star%z,h%r_star*0.2d0,FAM_STAR)
       else
          h%p_star%x=h%p%x
          h%p_star%y=h%p%y
          h%p_star%z=h%p%z
       end if
       if (h%n_DM >= nMembers) then
          ! print*,'-DM-', h%my_number, h%n_star, h%n_DM, nmembers
          call det_halo_center_sphere(h,h%p_DM%x,h%p_DM%y,h%p_DM%z,h%r_DM*0.2d0,FAM_DM)
       else
          h%p_DM%x=h%p%x
          h%p_DM%y=h%p%y
          h%p_DM%z=h%p%z
       end if
#endif
!!$       write(*,*)h%p%x,h%p%y,h%p%z
    endif

    ! compute principal axis of halo
    call det_main_axis(h,v)

    ! compute halo energies if necessary i.e. in the case where the center of the halo is the
    ! center of mass bc if it is the most bound particle then energies are already available
    call det_halo_energies(h)

    ! compute virial properties based on conditions necessary for the virial theorem to apply
#ifdef ANG_MOM_OF_R
    call det_vir_props(h,v,amax,bmax,cmax)
    call det_ang_momentum_per_shell(h,amax,bmax,cmax,v)
#else
     call det_vir_props(h,v)
#endif
    return

  end subroutine det_vir

!***********************************************************************
  subroutine det_age()

  ! subroutine which determines the current age of the Universe (Gyr) and
  ! values of current cosmological parameters

    implicit none

    external midpnt
    real(kind=8) :: age0,age1,somme0,somme1,omm,oml

    save age0

    omm  = omega_f
    oml  = omega_lambda_f
    ! age0 is the age of the universe @ the beginning of the simulation (in Gyr)
    call qromo(temps,omm,oml,(af/ai),real(10001.,8),somme0,midpnt)
    age0 = 977.78*somme0/H_f

    ! age1 is the time between the beginning of the simulation and the current timestep (in Gyr)
    if (aexp /= ai) then
       call qromo(temps,omm,oml,(af/aexp),(af/ai),somme1,midpnt)
    else
       somme1 = 0.0
    endif
    age1 = 977.78*somme1/H_f

    ! finally the age of the universe is the sum
    age_univ = age0+age1

    write(errunit,*)
    write(errunit,*) '> Current values of parameters: '
    write(errunit,*) '> ----------------------------- '
    write(errunit,*)
    write(errunit,*) '> Redshift                    : ',af/aexp-1.
    write(errunit,*) '> Age of the Universe (Gyr)   : ',age_univ

    Hub_pt   = H_f * sqrt(omega_f*(af/aexp)**3 +  omega_c_f*(af/aexp)**2 + omega_lambda_f)
    Lbox_pt  = Lboxp*(aexp/ai)
    Lbox_pt2 = Lbox_pt / 2.0

    write(errunit,*) '> Hubble Parameter  (km/s/Mpc): ',Hub_pt
    write(errunit,*) '> Box Length (Mpc)            : ',Lbox_pt
#ifdef STARS
    ! Initialize conformal time things
    call ct_init_cosmo(omega_f,omega_lambda_f,omega_c_f,H_f)
    ! lbtime of -10 Gyr means "10 Gyr ago"
    ! so lbtime(now) - lbtime(z) is the time between z and now
    ! it is positive if z is before now.
    lbtime = ct_aexp2time(aexp)
    write(errunit,*) '> Lookback time (Gyr)         : ',lbtime/1d9
#endif
    write(errunit,*)

    return

  end subroutine det_age

!***********************************************************************
  subroutine virial()

    ! compute the overdensity factor for virialization in a tophat collapse model at a given redshift
    ! for a given cosmological model

    implicit none

    real(kind=8) :: a,b,eta,omega_maxexp,age,reduce,cubic

    ! convert age of universe from Gyr back into inverse Hubble parameter units
    age       = age_univ/977.78*H_f
    ! compute the overdensity needed to reach maximum expansion by half the age of the universe
    omega_maxexp = omega_f
    call collapse(age/2.0,omega_maxexp,1.d-6)
    ! calculate how far an object collapses to virial equilibrium
    eta = 2.0*omega_lambda_f/omega_maxexp*(af/aexp)**3
    if (eta == 0.0) then
       reduce = 0.5
    else
       a      = 2.0*eta
       b      = -(2.0+eta)
       reduce = cubic(a,0.0,b,1.0)
    end if

    vir_overdens = omega_maxexp/omega_f/reduce**3*(aexp/af)**3
    rho_mean     = mboxp/Lbox_pt**3

    return

  end subroutine virial

!***********************************************************************
  subroutine collapse(age0,omm,acc)
    ! this subroutine performs a brute-force search using bracketing to find the value of the cosmic curvature
    ! that gives turnaround at the specified expansion parameter. The expansion parameter is defined to
    ! be 1 at turnaround.
    !
    ! note that for a constant cosmological constant, the age at a given
    ! expansion factor decreases monotonically with omegam (before turnaround;
    ! after, there is more than one age at a given expansion factor).
    !
    ! third argument is the desired fractional accurracy.

    implicit none

    real(kind=8)            :: age0,omm,acc,oml
    real(kind=8)            :: age,omax,omin,age_f
    real(kind=8), parameter :: omax0=1.d7,omin0=1.d0
    ! IMPORTANT NOTE: omax0 corresponds to a perturbation turning around at z=140 in a LCDM standard cosmology
    !                 and needs to be increased if you analyze outputs before this redshift ...

    external midpnt

    age  = -1.0  ! impossible value
    omax = omax0
    omin = omin0
    oml  = omega_lambda_f
    do while (abs(age-age0) .gt. acc*age0 .and. omax-omin .gt. acc*omm)
       omm = 0.5*(omax+omin)
       call qromo(temps_turn_around,omm,oml,real(101,8),real(10001,8),age_f,midpnt)
       call qromo(temps_turn_around,omm,oml,real(1,8),real(101,8),age,midpnt)
       if (age+age_f .gt. age0) then
          omin = omm
       else
          omax = omm
       end if
    end do

    if (omax .eq. omax0 .or. omin .eq. omin0) then
       write (errunit,*) 'WARNING: presumed bounds for omega are inadequate in collapse.'
       write (errunit,*) 'WARNING: omax,omax0,omin,omin0=',omax,omax0,omin,omin0
    end if

    return

  end subroutine collapse

!***********************************************************************
  subroutine det_halo_energies(h)

    implicit none

    integer(kind=i8b)           :: indexp,indexpp,np
    integer(kind=4),parameter :: full_PE = 1000 ! below this number of parts, we calculate full potential energy
    real(kind=8)              :: v2,rf          ! rf is elliptic integral function from numrec
    real(kind=8)              :: ped,ked,ped_tmp ! need hi precision for potential and ke energy sum.
!!$#ifdef ALLPARTS
!!$    real(kind=8)              :: ped_star, ked_star, ped_DM, ked_DM
!!$#endif
    logical(kind=4)           :: count_pairs
    type (halo)               :: h
    type (vector)             :: vt,dr

    np          = nb_of_parts(h%my_number)
    count_pairs = (np < full_PE)
    ! get potential energy
    if (.not.count_pairs) then
       ! formula B1 of appendix of paper GalICS 1 :
       ! EP = -3/5 G M^2 Rf,
       ! with Rf = 0.5 * Int_0^infty dt/sqrt((t+x)(t+y)(t+z)) = 0.5 * rf_numrec
       ! NB : rf_numrec returns RF in inverse Mpc (because x is in Mpc^2)
       h%ep = (-0.3) * gravconst * h%m**2 * rf(h%sh%a**2,h%sh%b**2,h%sh%c**2)
!!$#ifdef ALLPARTS
!!$       h%ep_star = (-0.3) * gravconst * h%m_star**2 * rf(h%sh_star%a**2,h%sh_star%b**2,h%sh_star%c**2)
!!$       h%ep_DM = (-0.3) * gravconst * h%m_DM**2 * rf(h%sh_DM%a**2,h%sh_DM%b**2,h%sh_DM%c**2)
!!$#endif
    else
       indexp = first_part(h%my_number)
       ped    = 0d0
!!$#ifdef ALLPARTS
!!$       ped_star = 0d0
!!$       ped_DM   = 0d0
!!$#endif
       do while (indexp /= -1)
          indexpp = linked_list(indexp) ! only count pairs once
          do while (indexpp /= -1)
             ped_tmp = real(interact(indexp,indexpp),8)
             ped     = ped + ped_tmp
!!$#ifdef ALLPARTS
!!$          if (family(indexp) == FAM_STAR) then
!!$             ped_star = ped_star + ped_tmp
!!$          else
!!$             ped_DM = ped_DM + ped_tmp
!!$          end if
!!$#endif
             indexpp = linked_list(indexpp)
          end do
          indexp = linked_list(indexp)
       end do
       h%ep = real(ped,8) * gravconst / Lbox_pt
!!$#ifdef ALLPARTS
!!$       h%ep_star = real(ped_star,8) * gravconst / Lbox_pt
!!$       h%ep_DM   = real(ped_DM,8) * gravconst / Lbox_pt
!!$#endif
    end if

    ! get kinetic energy (in center-of-halo frame)
    indexp = first_part(h%my_number)
    ked    = 0d0
!!$#ifdef ALLPARTS
!!$    ked_star = 0d0
!!$    ked_DM   = 0d0
!!$#endif
    do while (indexp /= -1)
       vt%x = vel(indexp,1) - h%v%x
       vt%y = vel(indexp,2) - h%v%y
       vt%z = vel(indexp,3) - h%v%z
       dr%x = pos(indexp,1) - h%p%x
       dr%y = pos(indexp,2) - h%p%y
       dr%z = pos(indexp,3) - h%p%z
       call correct_for_periodicity(dr)
       ! add Hubble flow
       vt%x = vt%x + dr%x * Hub_pt
       vt%y = vt%y + dr%y * Hub_pt
       vt%z = vt%z + dr%z * Hub_pt
       v2   = vt%x**2 + vt%y**2 + vt%z**2
       if (allocated(mass)) then
          ked = ked + real(mass(indexp)*v2,8)
!!$#ifdef ALLPARTS
!!$          if (family(indexp) == FAM_STAR) then
!!$             ked_star = ked_star + real(mass(indexp)*v2,8)
!!$          else
!!$             ked_DM = ked_DM + real(mass(indexp)*v2,8)
!!$          end if
!!$#endif
       else
          ked = ked + real(massp*v2,8)
       endif
       indexp  = linked_list(indexp)
    end do
    h%ek = 0.5*real(ked,8)
!!$#ifdef ALLPARTS
!!$    h%ek_star = 0.5*real(ked_star,8)
!!$    h%ek_DM = 0.5*real(ked_DM,8)
!!$#endif

    ! get total energy
    h%et = h%ek + h%ep
!!$#ifdef ALLPARTS
!!$    ! Should this be ek_star+ep, and ek_DM + ep?
!!$    h%et_star = h%ek_star + h%ep_star
!!$    h%et_DM = h%ek_DM + h%ep_DM
!!$#endif

    return

  end subroutine det_halo_energies

!***********************************************************************
recursive subroutine det_halo_center_multiscale(h,x0,y0,z0,r0,nx,fam)

    implicit none

    integer(i8b)              :: indexp, itarget
    type (halo)               :: h
    integer(kind=4)           :: ii,jj,kk,imax,jmax,kmax,nx,nxnew
    real(kind=8)              :: xmin,ymin,zmin,deltax,mass_max
    real(kind=8)              :: x0,y0,z0,r0,xc,yc,zc
    real(kind=8),dimension(:,:,:),allocatable:: mass_grid
    integer(kind=1)    :: fam
    logical            :: sel=.false.

    allocate(mass_grid(1:nx,1:nx,1:nx))
    mass_grid=0.0d0

    xmin=x0-r0
    ymin=y0-r0
    zmin=z0-r0
    deltax=2d0*r0/dble(nx)

    ! Assign mass to a uniform mesh with NGP
    indexp = first_part(h%my_number)
    do while (indexp /= -1)
       if (fam == FAM_UNDEF) then
          ! No family defined
          sel = .true.
       else
          ! Select only particles with the right type
          sel = (family(indexp) == fam)
       end if
       if (sel) then
          ii= int( (pos(indexp,1) - xmin)/deltax )+1
          jj= int( (pos(indexp,2) - ymin)/deltax )+1
          kk= int( (pos(indexp,3) - zmin)/deltax )+1
          if(ii>0.and.ii<=nx.and.jj>0.and.jj<=nx.and.kk>0.and.kk<=nx)then
             if (allocated(mass)) then
                mass_grid(ii,jj,kk)=mass_grid(ii,jj,kk)+real(mass(indexp),8)
             else
                mass_grid(ii,jj,kk)=mass_grid(ii,jj,kk)+real(massp,8)
             endif
          endif
       end if
       indexp  = linked_list(indexp)
    enddo

    ! Search for the cell containing the maximum mass
    mass_max=0.0d0
    do ii=1,nx
    do jj=1,nx
    do kk=1,nx
       if(mass_grid(ii,jj,kk)>mass_max)then
          imax=ii
          jmax=jj
          kmax=kk
          mass_max=mass_grid(ii,jj,kk)
       endif
    enddo
    enddo
    enddo

    deallocate(mass_grid)


    !if(deltax > 0.02*h%r)then
    !if(deltax > 2.0d0*dcell_min)then

    nxnew=3
    if(deltax > dble(nxnew)*dcell_min)then
       xc=xmin+(dble(imax)-0.5d0)*deltax
       yc=ymin+(dble(jmax)-0.5d0)*deltax
       zc=zmin+(dble(kmax)-0.5d0)*deltax
       call det_halo_center_multiscale(h,xc,yc,zc,deltax,nxnew,fam)
    else

       ! Find the particle with the maximum density within
       ! the cell with the maximum mass
       mass_max=0.0d0
       indexp = first_part(h%my_number)
       do while (indexp /= -1)
          if (fam == FAM_UNDEF) then
             ! No family defined
             sel = .true.
          else
             ! Select only particles with the right type
             sel = (fam == family(indexp))
          end if
          if (sel) then
             ii= int( (pos(indexp,1) - xmin)/deltax )+1
             if(ii==imax)then
                jj= int( (pos(indexp,2) - ymin)/deltax )+1
                if(jj==jmax)then
                   kk= int( (pos(indexp,3) - zmin)/deltax )+1
                   if(kk==kmax)then

                      if(density(indexp)>mass_max)then
                         mass_max=density(indexp)
                         itarget=indexp
                      endif

                   endif
                endif
             endif
          endif
          indexp  = linked_list(indexp)
       enddo

       ! Assign the new halo center
       if (fam == FAM_UNDEF) then
          h%p%x=pos(itarget,1)
          h%p%y=pos(itarget,2)
          h%p%z=pos(itarget,3)
#ifdef ALLPARTS
       else
          if (fam == FAM_STAR) then
             h%p_star%x=pos(itarget,1)
             h%p_star%y=pos(itarget,2)
             h%p_star%z=pos(itarget,3)
          else
             h%p_DM%x=pos(itarget,1)
             h%p_DM%y=pos(itarget,2)
             h%p_DM%z=pos(itarget,3)
          end if
#endif
       end if
    endif

    return

  end subroutine det_halo_center_multiscale
!***********************************************************************

!***********************************************************************
recursive subroutine det_halo_center_sphere(h,x0,y0,z0,r0,fam)

    implicit none

    integer(i8b)       :: indexp,ifirst,itarget
    type (halo)        :: h
    integer(kind=4)    :: nxnew
    real(kind=8)       :: x0,y0,z0,r0,r02,xc,yc,zc,mtot
    real(kind=8)       :: distmin
    real(kind=8)       :: pcx,pcy,pcz,dr2,deltax,dmsmove
    type(vector)       :: dr,pc
    integer(kind=1)    :: fam
    logical            :: sel=.false.

    deltax=2d0*r0
    r02=r0*r0

    ! compute cdm
    pcx   = 0d0 ; pcy   = 0d0 ; pcz = 0d0 ; mtot=0d0
    ifirst = first_part(h%my_number)
    indexp = ifirst
    do while (indexp /= -1)
       if (fam == FAM_UNDEF) then
          ! No family defined
          sel = .true.
       else
          ! Select only particles with the right type
          sel = (family(indexp) == fam)
       end if
       if (sel) then
          dr%x = pos(indexp,1) - x0
          dr%y = pos(indexp,2) - y0
          dr%z = pos(indexp,3) - z0
          call correct_for_periodicity(dr)
          dr2=    dr%x*dr%x
          if(dr2.le.r02)then
             dr2=dr2+dr%y*dr%y
             if(dr2.le.r02)then
                dr2=dr2+dr%z*dr%z
                if(dr2.le.r02)then
                   if (allocated(mass)) then
                      pcx = pcx + real(mass(indexp)*dr%x,8)
                      pcy = pcy + real(mass(indexp)*dr%y,8)
                      pcz = pcz + real(mass(indexp)*dr%z,8)
                      mtot= mtot+ real(mass(indexp),8)
                   else
                      pcx = pcx + real(massp*dr%x,8)
                      pcy = pcy + real(massp*dr%y,8)
                      pcz = pcz + real(massp*dr%z,8)
                      mtot= mtot+ real(massp,8)
                   end if
                endif
             end if
          endif
       endif
       indexp = linked_list(indexp)
    end do
    if(mtot > 0)then
       xc  = pcx / real(mtot,8) + real(x0,8)
       yc  = pcy / real(mtot,8) + real(y0,8)
       zc  = pcz / real(mtot,8) + real(z0,8)
    else
       xc = x0
       yc = y0
       zc = z0
    endif

    dmsmove=sqrt( (xc-x0)**2 +(yc-y0)**2 +(zc-z0)**2 )

    nxnew=4 ! This means 2*r0 > 4*d_cell_min
    if(deltax > dble(nxnew)*dcell_min .and. mtot>0)then
       call det_halo_center_sphere(h,xc,yc,zc,(1d0-eps_SC)*r0,fam)
    else if(MS .and. dmsmove>dms_min)then
       write(*,*)h%my_number,'MS:',dmsmove,dms_min
       call det_halo_center_sphere(h,xc,yc,zc,dcell_min,fam)
    else

       pc%x = real(xc,8)
       pc%y = real(yc,8)
       pc%z = real(zc,8)
       call correct_for_periodicity(pc)
       ! search particule closest to the cdm
       distmin = r02
       itarget = -1
       do while (itarget == -1)
          indexp  = ifirst
          do while (indexp /= -1)
             if (fam == FAM_UNDEF) then
                ! No family defined
                sel = .true.
             else
                ! Select only particles with the right type
                sel = (fam == family(indexp))
             end if
             dr%x = pos(indexp,1) - pc%x
             dr%y = pos(indexp,2) - pc%y
             dr%z = pos(indexp,3) - pc%z
             call correct_for_periodicity(dr)
             dr2=dr%x**2+dr%y**2+dr%z**2
             if (dr2.lt.distmin .and. sel) then
                itarget = indexp
                distmin = dr2
             end if
             indexp = linked_list(indexp)
          end do
          distmin=distmin*2d0**2
          if(distmin > 1.0) then
             stop
          end if
       enddo
       ! Assign the new halo center
       if (fam == FAM_UNDEF) then
          h%p%x=pos(itarget,1)
          h%p%y=pos(itarget,2)
          h%p%z=pos(itarget,3)
#ifdef ALLPARTS
       else
          if (fam == FAM_STAR) then
             h%p_star%x=pos(itarget,1)
             h%p_star%y=pos(itarget,2)
             h%p_star%z=pos(itarget,3)
          else
             h%p_DM%x=pos(itarget,1)
             h%p_DM%y=pos(itarget,2)
             h%p_DM%z=pos(itarget,3)
          end if
#endif
       end if
    endif

    return

  end subroutine det_halo_center_sphere
!***********************************************************************

  subroutine compute_spin_parameter(h)

    implicit none

    real(kind=8)                :: hl,spin
    type (halo)                 :: h

    hl                  = h%L%x**2 + h%L%y**2 + h%L%z**2
    hl                  = sqrt(hl)
    spin                = hl * sqrt(abs(h%et)) / h%m**2.5
    spin                = spin / gravconst
    h%spin              = spin
!!$#ifdef ALLPARTS
!!$    ! First, do the DM
!!$    hl                  = h%L_DM%x**2 + h%L_DM%y**2 + h%L_DM%z**2
!!$    hl                  = sqrt(hl)
!!$    spin                = hl * sqrt(abs(h%et_DM)) / h%m_DM**2.5
!!$    spin                = spin / gravconst
!!$    h%spin_DM           = spin
!!$    ! Then, the star
!!$    hl                  = h%L_star%x**2 + h%L_star%y**2 + h%L_star%z**2
!!$    hl                  = sqrt(hl)
!!$    spin                = hl * sqrt(abs(h%et_star)) / h%m_star**2.5
!!$    spin                = spin / gravconst
!!$    h%spin_star         = spin
!!$#endif

    return

  end subroutine compute_spin_parameter

!***********************************************************************
  subroutine det_inertial_tensor(h,mat,fam)

  ! Compute inertial tensor with respect to center of halo (either cdm or mbp)

    implicit none

    integer(kind=4) :: num_h
    integer(i8b)    :: indexp
    integer(kind=1) :: fam
    logical         :: sel=.false.
    real(kind=8)    :: mat(1:3,1:3)
    real(kind=8)    :: md(1:3,1:3)
    type (vector)   :: dr
    type (halo)     :: h

    num_h  = h%my_number
    md     = 0d0
    indexp = first_part(num_h)

    ! select all particles if no family is specified
    do while (indexp /= -1)
       if (fam == FAM_UNDEF) then
          ! No family defined
          sel = .true.
       else
          ! Select only particles with the right type
          sel = (family(indexp) == fam)
       end if

       dr%x=pos(indexp,1)-h%p%x
       dr%y=pos(indexp,2)-h%p%y
       dr%z=pos(indexp,3)-h%p%z

       call correct_for_periodicity(dr)

       if (sel) then
          if (allocated(mass)) then
             md(1,1) = md(1,1) + real(mass(indexp)*dr%x*dr%x,8)
             md(1,2) = md(1,2) + real(mass(indexp)*dr%x*dr%y,8)
             md(1,3) = md(1,3) + real(mass(indexp)*dr%x*dr%z,8)
             md(2,1) = md(2,1) + real(mass(indexp)*dr%x*dr%y,8)
             md(2,2) = md(2,2) + real(mass(indexp)*dr%y*dr%y,8)
             md(2,3) = md(2,3) + real(mass(indexp)*dr%y*dr%z,8)
             md(3,1) = md(3,1) + real(mass(indexp)*dr%x*dr%z,8)
             md(3,2) = md(3,2) + real(mass(indexp)*dr%y*dr%z,8)
             md(3,3) = md(3,3) + real(mass(indexp)*dr%z*dr%z,8)
          else
             md(1,1) = md(1,1) + real(massp*dr%x*dr%x,8)
             md(1,2) = md(1,2) + real(massp*dr%x*dr%y,8)
             md(1,3) = md(1,3) + real(massp*dr%x*dr%z,8)
             md(2,1) = md(2,1) + real(massp*dr%x*dr%y,8)
             md(2,2) = md(2,2) + real(massp*dr%y*dr%y,8)
             md(2,3) = md(2,3) + real(massp*dr%y*dr%z,8)
             md(3,1) = md(3,1) + real(massp*dr%x*dr%z,8)
             md(3,2) = md(3,2) + real(massp*dr%y*dr%z,8)
             md(3,3) = md(3,3) + real(massp*dr%z*dr%z,8)
          endif
       end if
       indexp = linked_list(indexp)

    end do

!!$    mat(1,1) = real(md(1,1),8)
!!$    mat(1,2) = real(md(1,2),8)
!!$    mat(1,3) = real(md(1,3),8)
!!$    mat(2,1) = real(md(2,1),8)
!!$    mat(2,2) = real(md(2,2),8)
!!$    mat(2,3) = real(md(2,3),8)
!!$    mat(3,1) = real(md(3,1),8)
!!$    mat(3,2) = real(md(3,2),8)
!!$    mat(3,3) = real(md(3,3),8)

    mat(1,1) = md(1,1)
    mat(1,2) = md(1,2)
    mat(1,3) = md(1,3)
    mat(2,1) = md(2,1)
    mat(2,2) = md(2,2)
    mat(2,3) = md(2,3)
    mat(3,1) = md(3,1)
    mat(3,2) = md(3,2)
    mat(3,3) = md(3,3)

    return

  end subroutine det_inertial_tensor

!***********************************************************************
  subroutine det_main_axis(h,v)

  ! determine the principal axis of the halo (h%sh%a,b,c)

    implicit none

    integer(kind=4) :: nrot
    real(kind=8)    :: mat(1:3,1:3)
    real(kind=8)    :: d(3),v(3,3)
    type (halo)     :: h


    call det_inertial_tensor(h,mat,FAM_UNDEF)

    call jacobi(mat,3,d,v,nrot)

    d      = sqrt(d/h%m)
    h%sh%a = d(1)
    h%sh%b = d(2)
    h%sh%c = d(3)


#ifdef ALLPARTS
    ! Do the same for stars
    call det_inertial_tensor(h,mat,FAM_STAR)
    call jacobi(mat,3,d,v,nrot)
    d      = sqrt(d/h%m_star)
    h%sh_star%a = d(1)
    h%sh_star%b = d(2)
    h%sh_star%c = d(3)

    ! And for the DM
    call det_inertial_tensor(h,mat,FAM_DM)
    call jacobi(mat,3,d,v,nrot)
    d      = sqrt(d/h%m_DM)
    h%sh_DM%a = d(1)
    h%sh_DM%b = d(2)
    h%sh_DM%c = d(3)
#endif

    return

  end subroutine det_main_axis

!***********************************************************************
  subroutine fit_deltac

    ! compute delta_crit(z) using the fit form from Bryan&Norman(1998)

    implicit none

    real(kind=4)    :: z,E2z,x,rho_crit

    z = af/aexp-1.
    E2z = omega_f * (1.+z)**3. + omega_lambda_f
    omega_z = omega_f * (1.+z)**3. / E2z
    x = omega_z - 1.

    delta_crit_z = (18*pi**2 + 82.*x - 39*x**2)

    rho_crit = 2.78782*(H_f/100.)**2  ! (critical density at z=0 in units of 10**11 M_sol/Mpc^3)
    rho_crit_z = rho_crit * E2z

    return

  end subroutine fit_deltac

!***********************************************************************
#ifdef ANG_MOM_OF_R
  subroutine det_ang_momentum_per_shell(h,amax,bmax,cmax,v)

    implicit none

    integer(i8b) :: indexp
    type(halo)      :: h
    type(vector)    :: dr, dv
    real(kind=8)    :: amax,bmax,cmax ! computed in tab_props_inside
    real(kind=8)    :: v(3,3) ! computed in tab_props_inside
    real(kind=8)    :: dra, drb, drc, r_ell
    integer(kind=4) :: i_ell
    type(vector)    :: L(nshells)  ! ang mom of a shell (in 10**11 Msun * km/s * Mpc)
    real(kind=8)    :: m(nshells)  ! mass of a shell (in 10**11 Msun)
    integer(kind=4) :: i

    do i = 1,nshells
       L(i)%x = 0.0
       L(i)%y = 0.0
       L(i)%z = 0.0
       m(i)   = 0.0
    end do
    indexp = first_part(h%my_number)
    do while(indexp /= -1)
       ! particle positions relative to center of halo
       dr%x   = pos(indexp,1) - h%p%x
       dr%y   = pos(indexp,2) - h%p%y
       dr%z   = pos(indexp,3) - h%p%z
       call correct_for_periodicity(dr)
       ! convert dr into ellipsoid coords
       dra    = dr%x*v(1,1)+dr%y*v(2,1)+dr%z*v(3,1)
       drb    = dr%x*v(1,2)+dr%y*v(2,2)+dr%z*v(3,2)
       drc    = dr%x*v(1,3)+dr%y*v(2,3)+dr%z*v(3,3)
       ! index of shell containing particle
       r_ell  = sqrt((dra / amax)**2 + (drb / bmax)**2 + (drc / cmax)**2)
       i_ell  = int(r_ell*nshells)
       if (i_ell > nshells) then
          write(errunit,'(a)') '> Problem in get_ang_momentum_per_shell : i_ell > nshells '
          write(errunit,*) i_ell,nshells
          stop
       end if
       ! velocity relative to halo velocity (cdm)
       dv%x   = vel(indexp,1)-h%v%x
       dv%y   = vel(indexp,2)-h%v%y
       dv%z   = vel(indexp,3)-h%v%z
       ! update mass and angular momentum of shell
#ifdef BIG_RUN
       m(i_ell)   = m(i_ell) + real(massp,8)
       L(i_ell)%x = L(i_ell)%x + massp*(dr%y*dv%z - dr%z*dv%y) ! in 10**11 Msun * km/s * Mpc
       L(i_ell)%y = L(i_ell)%y + massp*(dr%z*dv%x - dr%x*dv%z)
       L(i_ell)%z = L(i_ell)%z + massp*(dr%x*dv%y - dr%y*dv%x)
#else
       m(i_ell)   = m(i_ell) + real(mass(indexp),8)
       L(i_ell)%x = L(i_ell)%x + mass(indexp)*(dr%y*dv%z - dr%z*dv%y)
       L(i_ell)%y = L(i_ell)%y + mass(indexp)*(dr%z*dv%x - dr%x*dv%z)
       L(i_ell)%z = L(i_ell)%z + mass(indexp)*(dr%x*dv%y - dr%y*dv%x)
#endif
       indexp = linked_list(indexp)
    end do

    write(agor_unit) h%my_number, h%p%x, h%p%y, h%p%z, h%v%x, h%v%y, h%v%z, &
         h%datas%rvir, h%datas%mvir, h%m, h%r, h%spin,                      &
         amax, bmax, cmax, v(1,1:3), v(2,1:3), v(3,1:3),                    &
         m(1:nshells), L(1:nshells)%x, L(1:nshells)%y, L(1:nshells)%z

    return

  end subroutine det_ang_momentum_per_shell
#endif
!***********************************************************************
!///////////////////////////////////////////////////////////////////////

end module compute_halo_props
