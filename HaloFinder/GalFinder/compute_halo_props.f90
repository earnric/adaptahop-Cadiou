module compute_halo_props

  use input_output
  use halo_defs

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

    ! This routine reads in the input_GalaxyMaker.dat file which contains the cosmological 
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

    integer(kind=4)                      :: indexp,ierr,i
    integer(kind=4)                      :: found,n_halo_contam,n_subs_contam
    real(kind=8)                         :: read_time_ini,read_time_end
    real(kind=8)                         :: t0,t1
    logical                              :: printdatacheckhalo !put to true if bug after make_linked_list 
    integer(kind=4)  :: ih
    ! leo            
    real(kind=4)     :: fhalo
    ! end leo        
#ifdef ANG_MOM_OF_R
    character(200)                       :: filename
#endif
#ifdef Test_FOF
    character(200)                       :: filelisteparts
#endif

    write(errunit,'(1x,a16,i5)') '> Timestep  --->',numero_step
    write(errunit,*) '> -------------------'
    write(errunit,*) ''
    write(errunit,*) 'npart (new_step)=',npart

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
       if(allocated(age_st))deallocate(age_st)
       if(allocated(met_st))deallocate(met_st)
       if(allocated(pf_st))deallocate(pf_st)  ! RS - Deallocate pristine frac array
       if(allocated(pz_st))deallocate(pz_st)  ! RS - Deallocate primordial Z array
       if(allocated(chem_st))deallocate(chem_st)
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
          do while (indexp /= -1 .and. found == 0) 
             if (mass(indexp) > massp* 1.00001 ) then
                if(fsub) then
                   if (level(i) == 1) then 
                      n_halo_contam = n_halo_contam + 1
                   else
                      n_subs_contam = n_subs_contam + 1
                   endif
                else
                      n_halo_contam = n_halo_contam + 1
                end if
                found = 1
             endif
             indexp = linked_list(indexp)
          enddo
       enddo
       write(errunit,*) '> # of halos, # of CONTAMINATED halos :',nb_of_halos,n_halo_contam
       write(errunit,*) '> # of subhalos, # of CONTAMINATED subhalos :',nb_of_subhalos,n_subs_contam
       !open(222,file='ncontam_halos.dat',status='unknown',position='append')
       !write(222,'(5(i6,2x))') numero_step,nb_of_halos,n_halo_contam,nb_of_subhalos,n_subs_contam
       !close(222)
    endif
     
    ! until now we were using code units for positions and velocities
    ! this routine changes that to physical (non comoving) coordinates for positions 
    ! and peculiar (no Hubble flow) velocities in km/s
    ! The masses are changed from code units into 10^11 M_sun as well.
    call change_units()

    ! allocation and initialization of the halo list
    allocate(liste_halos(0:(nb_of_halos+nb_of_subhalos)),stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) 'Cannot allocate liste_halos'
       stop
    endif

    call init_halos

#ifdef ANG_MOM_OF_R
    write(filename,'(a,a,i3.3)') trim(agor_file),'.',numstep
    open(unit=agor_unit,file=filename,status='unknown',form='unformatted')
    write(agor_unit) nb_of_halos,nb_of_subhalos
    write(agor_unit) nshells
#endif

    printdatacheckhalo = .false.
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
       ! compute Bulge properties of galaxies
       call compute_bulge_prop(liste_halos(i))
       if(printdatacheckhalo) write(errunit,*) '> m_bulge:',liste_halos(i)%m_bulge
       ! compute stellar surface density profiles
       call compute_stellar_profile(liste_halos(i))
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
       
       if(printdatacheckhalo) then
          call cpu_time(t1)
          write(errunit,*) '> halo computation took:',int(t1- t0) ,'s'
          write(errunit,*)
       end if
    end do
    !$OMP END PARALLEL DO
    ! end leo

#ifdef ANG_MOM_OF_R
    close(agor_unit)
#endif

    call write_tree_brick

    deallocate(liste_halos)
    deallocate(nb_of_parts,first_part,linked_list)
    deallocate(pos,vel)
    if(allocated(mass)) deallocate(mass)
    if(allocated(age_st))deallocate(age_st)
    if(allocated(met_st))deallocate(met_st)
    if(allocated(pf_st))deallocate(pf_st)  ! RS - Deallocate prist fraction
    if(allocated(pz_st))deallocate(pz_st)  ! RS - Deallocate primordial Z
    if(allocated(chem_st))deallocate(chem_st)
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
    integer(kind=4)             :: nstar

    write(errunit,*) '> In routine make_halos '
    write(errunit,*) '> ----------------------'
    write(errunit,*) 'npart=',nstar
    
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
          write(errunit,*) '> GalaxyMaker is using Friend Of Friend algorithm'
          call fof_init
          fsub = .false.
       case("HOP")
          write(errunit,*) '> GalaxyMaker is using Adaptahop in order to' 
          write(errunit,*) '> Detect halos, subhaloes will not be selected'    
          call init_adaptahop
          fsub = .false.
       case("DPM")
          write(errunit,*) '> GalaxyMaker is using Adaptahop in order to' 
          write(errunit,*) '> Detect halos, and subhaloes with the Density Profile Method'    
          call init_adaptahop
          fsub = .true.
       case("MSM")
          write(errunit,*) '> GalaxyMaker is using Adaptahop in order to' 
          write(errunit,*) '> Detect halos, and subhaloes with the Most massive Subhalo Method'
          call init_adaptahop
          fsub = .true.
       case("BHM")
          write(errunit,*) '> GalaxyMaker is using Adaptahop in order to' 
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

    integer(kind=4)             :: i,index1,index2,ierr
    integer(kind=4),allocatable :: current_ptr(:)  

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
    integer(kind=4)    :: indexp,npch
    real(kind=8)       :: masshalo
    type(halo)         :: h
    
    masshalo      = 0d0
    npch          = 0
    indexp        = first_part(h%my_number)
    do while(indexp /= -1)
       if (allocated(mass)) then 
          masshalo = masshalo + real(mass(indexp),8)
       else
          masshalo = masshalo + real(massp,8) 
       endif
       npch = npch + 1
       indexp = linked_list(indexp)       
    end do
    h%m = real(masshalo,8)  ! in 10^11 M_sun

    if(npch.ne.nb_of_parts(h%my_number)) then
       write(errunit,*) '> Fatal error in det_mass for', h%my_number
       write(errunit,*) 'nb_of_parts, npch:',h%my_number,npch
       stop
    end if

    return

  end subroutine det_mass

!***********************************************************************
  subroutine compute_ang_mom(h)

  ! compute angular momentum of all halos

    implicit none

    integer(kind=4) :: indexp
    real(kind=8)    :: lx,ly,lz
    type (halo)     :: h
    type (vector)   :: dr,p

    ! we compute r * m * v, where r & v are pos and vel of halo particles relative to center of halo
    ! (particle closest to center of mass or most dense particle)

    indexp = first_part(h%my_number)
    lx =0d0 ; ly = 0d0 ; lz = 0d0
    do while(indexp /= -1)
       
       dr%x   = pos(indexp,1) - h%p%x
       dr%y   = pos(indexp,2) - h%p%y
       dr%z   = pos(indexp,3) - h%p%z
       
       call correct_for_periodicity(dr)
       
       if (allocated(mass)) then
          p%x = mass(indexp)*(vel(indexp,1)-h%v%x)
          p%y = mass(indexp)*(vel(indexp,2)-h%v%y)
          p%z = mass(indexp)*(vel(indexp,3)-h%v%z)
       else
          p%x = massp*(vel(indexp,1)-h%v%x)
          p%y = massp*(vel(indexp,2)-h%v%y)
          p%z = massp*(vel(indexp,3)-h%v%z)
       endif

       lx  = lx + real(dr%y*p%z - dr%z*p%y,8)   ! in 10**11 Msun * km/s * Mpc
       ly  = ly + real(dr%z*p%x - dr%x*p%z,8)
       lz  = lz + real(dr%x*p%y - dr%y*p%x,8)        
       
       indexp = linked_list(indexp)   
       
    end do

    h%L%x = real(lx,8)
    h%L%y = real(ly,8)
    h%L%z = real(lz,8)

    return

  end subroutine compute_ang_mom

!***********************************************************************
  subroutine r_halos(h)

  ! compute distance of the most remote particle (with respect to center of halo, which
  ! is either center of mass or most bound particle)

    implicit none

    integer(kind=4) :: indexp
    real(kind=8)    :: dr2max,dr2
    type (vector)   :: dr
    type (halo)     :: h
 
    dr2max  = 0.0
    indexp = first_part(h%my_number)

    do while(indexp /= -1)
       
       dr%x = pos(indexp,1) - h%p%x
       dr%y = pos(indexp,2) - h%p%y
       dr%z = pos(indexp,3) - h%p%z         
       
       call correct_for_periodicity(dr)
       
       dr2    = (dr%x*dr%x + dr%y*dr%y + dr%z*dr%z)
       
       if (dr2 > dr2max) then
          dr2max         = dr2
       endif
       
       indexp=linked_list(indexp)   
       
    end do
    
    h%r = sqrt(dr2max)

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
    integer(kind=4)    :: indexp, icenter,ifirst 
    real(kind=8)       :: maxdens, distmin
    real(kind=8)       :: pcx,pcy,pcz,vcx,vcy,vcz,vmean,v2mean,sigma2,vnorm
    type(vector)       :: dr,pc

    real(kind=8)       ::aaa,bbb,ccc,ddd,eee,half_radius,mhalf
    integer(kind=4)    ::i,j,nmax
    real(kind=8),dimension(:),allocatable::drr,mm,vxx,vyy,vzz

    icenter = -1

    if (cdm) then

       ! compute cdm
       pcx   = 0d0 ; pcy   = 0d0 ; pcz = 0d0
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
       ! search particule closest to the cdm
       indexp  = ifirst
       distmin = Lbox_pt
       do while (indexp /= -1)
          dr%x = pos(indexp,1) - pc%x
          dr%y = pos(indexp,2) - pc%y
          dr%z = pos(indexp,3) - pc%z
          call correct_for_periodicity(dr)
          if (sqrt(dr%x**2+dr%y**2+dr%z**2).lt.distmin) then
             icenter = indexp
             distmin = sqrt(dr%x**2+dr%y**2+dr%z**2)
          end if
          indexp = linked_list(indexp)
       end do

    else
    
       maxdens = 0.0
       indexp  = first_part(h%my_number)
       do while (indexp /= -1)
          if (density(indexp).gt.maxdens) then
             maxdens = density(indexp)
             icenter = indexp
          end if
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

    ! velocity of center is set equal velocity of center of mass:

    indexp = first_part(h%my_number)
    vcx = 0d0 ; vcy = 0d0 ; vcz =0d0
    v2mean= 0d0
    i=0
    do while (indexp /= -1)
       i=i+1
       if (allocated(mass)) then
          vcx = vcx + real(mass(indexp)*vel(indexp,1),8)
          vcy = vcy + real(mass(indexp)*vel(indexp,2),8)
          vcz = vcz + real(mass(indexp)*vel(indexp,3),8)
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

    indexp = first_part(h%my_number)
    sigma2= 0d0
    i=0
    do while (indexp /= -1)
       i=i+1
       vnorm = real(sqrt(vel(indexp,1)**2+vel(indexp,2)**2+vel(indexp,3)**2),8)
       if (allocated(mass)) then
          sigma2 = sigma2 + mass(indexp)*(vnorm-vmean)**2
       else
          sigma2 = sigma2 + massp       *(vnorm-vmean)**2
       endif
       indexp   = linked_list(indexp)
    end do
    h%sigma = real(sqrt(sigma2/h%m),8)

!!$
!!$    allocate(drr(1:nmax),mm(1:nmax),vxx(1:nmax),vyy(1:nmax),vzz(1:nmax))
!!$    indexp = first_part(h%my_number)
!!$    i=0
!!$    do while (indexp /= -1)
!!$       i=i+1
!!$       drr(i)=real(sqrt((pos(indexp,1)-h%p%x)**2+(pos(indexp,2)-h%p%y)**2 &
!!$            + (pos(indexp,3)-h%p%z)**2),8)
!!$       vxx(i)=real(vel(indexp,1),8)
!!$       vyy(i)=real(vel(indexp,2),8)
!!$       vzz(i)=real(vel(indexp,3),8)
!!$       if (allocated(mass)) then
!!$          mm(i)=real(mass(indexp),8)
!!$       else
!!$          mm(i)=real(massp,8)
!!$       endif
!!$       indexp   = linked_list(indexp)
!!$    end do
!!$
!!$    ! Sort radius of particles
!!$    do j=2,nmax
!!$       aaa=drr(j)
!!$       bbb=vxx(j)
!!$       ccc=vyy(j)
!!$       ddd=vzz(j)
!!$       eee=mm (j)
!!$       i=j-1
!!$       do while(drr(i).gt.aaa.and.i .gt.0)
!!$          drr(i+1)=drr(i)
!!$          vxx(i+1)=vxx(i)
!!$          vyy(i+1)=vyy(i)
!!$          vzz(i+1)=vzz(i)
!!$          mm (i+1)=mm (i)
!!$          i=i-1
!!$       enddo
!!$       drr(i+1)=aaa
!!$       vxx(i+1)=bbb
!!$       vyy(i+1)=ccc
!!$       vzz(i+1)=ddd
!!$       mm (i+1)=eee
!!$    enddo
!!$
!!$    ! Compute the radius of half mass and the velocity dispersion within
!!$    i=0
!!$    mhalf = 0d0
!!$    sigma2= 0d0
!!$    do while (mhalf.le.h%m/2d0)
!!$       i=i+1
!!$       mhalf=mhalf+mm(i)
!!$       half_radius=drr(i)
!!$       vnorm = real(sqrt(vel(indexp,1)**2+vel(indexp,2)**2+vel(indexp,3)**2),8)
!!$       sigma2 = sigma2 + mm(i)*(vnorm-vmean)**2
!!$    end do
!!$    write(*,*)mhalf,h%m,i
!!$    h%sigma = real(sqrt(sigma2/mhalf),8)
!!$    deallocate(drr,mm,vxx,vyy,vzz)

    return

  end subroutine det_center
  
!***********************************************************************                
  subroutine compute_bulge_prop(h)

  ! compute position of center of mass of halo, and its velocity.

    implicit none

    type (halo)        :: h
    integer(kind=4)    :: indexp, ifirst 
    real(kind=8)       :: vcz,vcr,vct,vmean,sigma2,vnorm
    real(kind=8)       :: Lx,Ly,Lz,Ltot,mbulge,dxx,dyy,dzz,dzcyl
    real(kind=8)       :: planecen_x,planecen_y,planecen_z,rx,ry,rz,modr
    real(kind=8)       :: thetax,thetay,thetaz
    integer(kind=4)    :: i,j

    Ltot=sqrt(h%L%x**2+h%L%y**2+h%L%z**2)
    Lx=h%L%x/Ltot
    Ly=h%L%y/Ltot
    Lz=h%L%z/Ltot

    ! Assume that the mean bulge velocity is the mean galaxy velocity
    vmean=sqrt(h%v%x**2+h%v%y**2+h%v%z**2)

    indexp = first_part(h%my_number)
    mbulge= 0d0
    sigma2= 0d0
    i=0
    do while (indexp /= -1)
       i=i+1
       dxx=pos(indexp,1)-h%p%x
       dyy=pos(indexp,2)-h%p%y
       dzz=pos(indexp,3)-h%p%z
       dzcyl=dxx*Lx+dyy*Ly+dzz*Lz
       ! Compute centre coord. of the plane where the particle lies
       ! and perp. to the disc plane 
       planecen_x=h%p%x+dzcyl*Lx
       planecen_y=h%p%y+dzcyl*Ly
       planecen_z=h%p%z+dzcyl*Lz
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
    integer(kind=4)    :: indexp, ifirst 
    real(kind=8)       :: vcz,vcr,vct,vmean,sigma2,vnorm
    real(kind=8)       :: Lx,Ly,Lz,Ltot,mbulge,dxx,dyy,dzz,dzcyl
    real(kind=8)       :: planecen_x,planecen_y,planecen_z,rx,ry,rz,modr
    real(kind=8)       :: thetax,thetay,thetaz
    integer(kind=4)    :: i,j
    real(kind=8)       :: rmax,dx
    real(kind=8),dimension(1:h%nbin):: m_tot,dv

    rmax=3d0*1d-3
    dx=rmax/dble(h%nbin)
    Ltot=sqrt(h%L%x**2+h%L%y**2+h%L%z**2)
    Lx=h%L%x/Ltot
    Ly=h%L%y/Ltot
    Lz=h%L%z/Ltot

    indexp = first_part(h%my_number)
    i=0
    m_tot=0d0
    do while (indexp /= -1)
       i=i+1
       dxx=pos(indexp,1)-h%p%x
       dyy=pos(indexp,2)-h%p%y
       dzz=pos(indexp,3)-h%p%z
       dzcyl=dxx*Lx+dyy*Ly+dzz*Lz
       ! Compute centre coord. of the plane where the particle lies
       ! and perp. to the disc plane 
       planecen_x=h%p%x+dzcyl*Lx
       planecen_y=h%p%y+dzcyl*Ly
       planecen_z=h%p%z+dzcyl*Lz
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

    return

  end subroutine compute_stellar_profile
!***********************************************************************
  function interact(i,j)

    implicit none

    integer(kind=4) :: i,j,ifirst
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
  subroutine tab_props_inside(h,nr,tabm2,tabk2,tabp2,v,amax,bmax,cmax)

  ! returns the cumulative mass contained in concentric ellipsoids centered on the center of the
  ! halo (cdm or mbp)

    implicit none

    integer(kind=4)        :: nr,num_h,indexp,i,i_ell,louped_parts  
    real(kind=8)           :: amax,bmax,cmax
    real(kind=8)           :: v(3,3)
    real(kind=8)           :: tabm2(0:nr-1),tabk2(0:nr-1),tabp2(0:nr-1) !! double
    real(kind=8)           :: srm,srk                                   !! double
    real(kind=8)           :: rmax,dra,drb,drc
    real(kind=8)           :: r_ell,v2,rf
    real(kind=8),parameter :: epsilon = 1.d-2
    type (vector)          :: posp,vt
    type (halo)            :: h

    ! rescale to get ellipsoid  concentric to principal ellipsoid
    ! which contains all the particles of the halo
    rmax = 0.0
    indexp = first_part(h%my_number)
    do while(indexp.gt.0)
       posp%x = pos(indexp,1) - h%p%x
       posp%y = pos(indexp,2) - h%p%y
       posp%z = pos(indexp,3) - h%p%z
       call correct_for_periodicity(posp)
        ! project position vector along the principal ellipsoid axis
       dra    = posp%x*v(1,1)+posp%y*v(2,1)+posp%z*v(3,1)
       drb    = posp%x*v(1,2)+posp%y*v(2,2)+posp%z*v(3,2)
       drc    = posp%x*v(1,3)+posp%y*v(2,3)+posp%z*v(3,3)
       r_ell  = sqrt((dra / h%sh%a)**2 + (drb / h%sh%b)**2 + (drc / h%sh%c)**2)
       rmax   = max(rmax,r_ell)
       indexp = linked_list(indexp)
    end do
    
    amax = rmax * h%sh%a * (1.0 + epsilon)
    bmax = rmax * h%sh%b * (1.0 + epsilon)
    cmax = rmax * h%sh%c * (1.0 + epsilon)
!!$    if(h%my_number==5708)then 
!!$       write(*,*)'buggy stuff'
!!$       write(*,*)'rmax=',rmax
!!$       write(*,*)'sha,b,c=',h%sh%a,h%sh%b,h%sh%c
!!$    endif

    ! initialize loop quantities
    tabm2        = 0d0
    tabk2        = 0d0
    louped_parts = 0
    num_h        = h%my_number
    indexp       = first_part(num_h)

    do while (indexp /= -1)

       posp%x = pos(indexp,1) - h%p%x
       posp%y = pos(indexp,2) - h%p%y
       posp%z = pos(indexp,3) - h%p%z
       call correct_for_periodicity(posp)
       ! compute velocities in the halo frame adding in the Hubble flow
       vt%x   = vel(indexp,1) - h%v%x + posp%x * Hub_pt
       vt%y   = vel(indexp,2) - h%v%y + posp%y * Hub_pt
       vt%z   = vel(indexp,3) - h%v%z + posp%z * Hub_pt
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
!!$             if(i_ell==-2147483648)then 
!!$                write(*,*)'buggy stuff'
!!$                write(*,*)'id=',h%my_number
!!$                write(*,*)'i_ell=',i_ell
!!$                write(*,*)'nr=',nr
!!$                write(*,*)'r_ell=',r_ell
!!$                write(*,*)'dra,b,c=',dra,drb,drc
!!$                write(*,*)'a,b,cmax=',amax,bmax,cmax
!!$             endif
             tabm2(i_ell) = tabm2(i_ell)+real(mass(indexp),8)
             tabk2(i_ell) = tabk2(i_ell)+0.5*real(mass(indexp),8)*v2
          else
             tabm2(i_ell) = tabm2(i_ell)+real(massp,8)
             tabk2(i_ell) = tabk2(i_ell)+0.5*real(massp,8)*v2
          endif
       else
          louped_parts = louped_parts + 1
       endif

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
       tabp2(i) = -0.3 * gravconst * tabm2(i)**2 * rf(h%sh%a**2,h%sh%b**2,h%sh%c**2)
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

    integer(kind=4)           :: i,ii
    ! ttab = 1000 bins for virial radius precision better than 1% of halo size 
    integer(kind=4),parameter :: ttab = 1000 
    real(kind=8)              :: rvir,mvir,kvir,pvir
    real(kind=8)              :: v(3,3)
    real(kind=8)              :: amax,bmax,cmax,avir,bvir,cvir
    real(kind=8)              :: tab_mass(0:ttab-1),tab_ekin(0:ttab-1),tab_epot(0:ttab-1)  !! double
    real(kind=8)              :: virth,virth_old,volmin
    type (halo)               :: h

    ! compute properties inside ttab concentric principal ellipsoids centered on center of halo
    call tab_props_inside(h,ttab,tab_mass,tab_ekin,tab_epot,v,amax,bmax,cmax)

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


        !compute the half mass radius : Reff
        do ii = ttab-1,1,-1
           !total mass is tab_mass(ttab-1)
          if (tab_mass(ii) < tab_mass(ttab-1)/2) exit
       enddo
       !those are not the virial parameters, but 
       !as the variable is not used afterwards and
       !is defined, we use it
       avir   = real(ii,8)/real(ttab-1,8)*amax
       bvir   = real(ii,8)/real(ttab-1,8)*bmax
       cvir   = real(ii,8)/real(ttab-1,8)*cmax
       rvir   = (avir**2+bvir**2+cvir**2)**(1./2.)
       h%datas%Reff         = rvir    ! in Mpc
       if (rvir > h%r) then
           write(*,*) 'halo bugged',h%my_number,mvir,rvir
           write(*,*) 'The half mass radius (ellipsoid)'
           write(*,*) 'is larger than the farthest particle'
         !   stop
       end if
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

    if(method.ne."FOF".and.SC) then
       x0=dble(h%p%x);y0=dble(h%p%y);z0=dble(h%p%z);r0=dble(h%r)
       call det_halo_center_sphere(h,x0,y0,z0,r0)
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

    integer(kind=4)           :: indexp,indexpp,np
    integer(kind=4),parameter :: full_PE = 1000 ! below this number of parts, we calculate full potential energy 
    real(kind=8)              :: v2,rf          ! rf is elliptic integral function from numrec
    real(kind=8)              :: ped,ked        ! need hi precision for potential and ke energy sum.   
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
    else
       indexp = first_part(h%my_number)
       ped    = 0d0
       do while (indexp /= -1)              
          indexpp = linked_list(indexp) ! only count pairs once
          do while (indexpp /= -1)
             ped     = ped + real(interact(indexp,indexpp),8)
             indexpp = linked_list(indexpp)
          end do
          indexp = linked_list(indexp)
       end do
       h%ep = real(ped,8) * gravconst / Lbox_pt
    end if

    ! get kinetic energy (in center-of-halo frame)
    indexp = first_part(h%my_number)
    ked    = 0d0
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
       else
          ked = ked + real(massp*v2,8)
       endif
       indexp  = linked_list(indexp)
    end do
    h%ek = 0.5*real(ked,8)
    ! get total energy 
    h%et = h%ek + h%ep
    
    return
    
  end subroutine det_halo_energies

!***********************************************************************
recursive subroutine det_halo_center_sphere(h,x0,y0,z0,r0)

    implicit none

    integer(kind=4)    :: indexp,ifirst,itarget
    type (halo)        :: h
    integer            :: nxnew
    real(kind=8)       :: x0,y0,z0,r0,r02,xc,yc,zc,mtot
    real(kind=8)       :: distmin
    real(kind=8)       :: pcx,pcy,pcz,dr2,dmsmove
    type(vector)       :: dr,pc

    r02=r0*r0

    ! compute cdm
    pcx   = 0d0 ; pcy   = 0d0 ; pcz = 0d0 ; mtot=0d0
    ifirst = first_part(h%my_number)
    indexp = ifirst
    do while (indexp /= -1) 
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
    
!!$    if(h%my_number == 94)then 
!!$       write(*,'(A,4f10.5)')'***',xc,yc,zc,r0
!!$    endif
    
    if(r0 > dcell_min .and. mtot>0)then
!!$       if(h%my_number == 2)write(*,*)'SC:',r0,dcell_min,dmsmove
       call det_halo_center_sphere(h,xc,yc,zc,(1d0-eps_SC)*r0)
    else if(MS .and. dmsmove>dms_min)then
       write(*,*)h%my_number,'MS:',dmsmove,dms_min
       call det_halo_center_sphere(h,xc,yc,zc,dcell_min)       
    else
!!$    if(r0 > dcell_min .and. mtot>0)then
!!$       call det_halo_center_sphere(h,xc,yc,zc,(1d0-eps_SC)*r0)
!!$    else

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
             dr%x = pos(indexp,1) - pc%x
             dr%y = pos(indexp,2) - pc%y
             dr%z = pos(indexp,3) - pc%z
             call correct_for_periodicity(dr)
             dr2=dr%x**2+dr%y**2+dr%z**2
             if (dr2.lt.distmin) then
                itarget = indexp
                distmin = dr2
             end if
             indexp = linked_list(indexp)
          end do
          distmin=distmin*2d0**2
          if(distmin > 1.0)stop
       enddo
       ! Assign the new halo center
       h%p%x=pos(itarget,1)
       h%p%y=pos(itarget,2)
       h%p%z=pos(itarget,3)
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

    return

  end subroutine compute_spin_parameter

!***********************************************************************
  subroutine det_inertial_tensor(h,mat)

  ! Compute inertial tensor with respect to center of halo (either cdm or mbp)

    implicit none

    integer(kind=4) :: num_h,indexp
    real(kind=8)    :: mat(1:3,1:3)
    real(kind=8)    :: md(1:3,1:3)
    type (vector)   :: dr
    type (halo)     :: h

    num_h  = h%my_number
    md     = 0d0
    indexp = first_part(num_h)   

    do while (indexp /= -1)

       dr%x=pos(indexp,1)-h%p%x
       dr%y=pos(indexp,2)-h%p%y
       dr%z=pos(indexp,3)-h%p%z

       call correct_for_periodicity(dr)

!!$       if(h%my_number==5708)then 
!!$          write(*,*)dr%x,dr%y,dr%z ! buggy stuff
!!$       endif

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


    call det_inertial_tensor(h,mat)

    call jacobi(mat,3,d,v,nrot)

    d      = sqrt(d/h%m)
    h%sh%a = d(1)
    h%sh%b = d(2)
    h%sh%c = d(3)

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
    
    integer(kind=4) :: indexp
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
