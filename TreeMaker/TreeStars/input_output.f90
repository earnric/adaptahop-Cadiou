module input_output

  use tree_defs

  public

contains

!//////////////////////////////////////////////////////////////////////////
!**************************************************************************
  subroutine get_max_n_parts

    implicit none
    integer(kind=4) :: i,number_of_halos,number_of_subhalos,number_of_particles
    real(kind=8)    :: dummy
    character(80)   :: file,file2
    
    write(file,'(a,a)') trim(data_dir),'./input_TreeMaker.dat'

    open(unit=12,file=file,status='old')
    read(12,*) nsteps, n_tree_files

    max_groups = 0
    nbodies    = 0
    do i = 1,nsteps    
       read(12,*) file2
       write(file,'(a,a)') trim(data_dir),trim(file2)
       open(unit=13,file=file,form='unformatted',status='old')
       read(13) number_of_particles
       read(13) dummy
       read(13) dummy
       read(13) dummy
       read(13) dummy
       read(13) number_of_halos, number_of_subhalos
       close(13)
       ! the number of star particles varies between outputs so set it to the max
       nbodies = max(nbodies,number_of_particles)
       if (max_groups .lt. (number_of_halos + number_of_subhalos )) then
          max_groups = number_of_halos + number_of_subhalos
       end if
    end do
    
    close(12)

    return

  end subroutine get_max_n_parts

!***********************************************************************
  subroutine open_exist(unit,name)

    implicit none

    integer(kind=4) :: unit,lnblnk
    character(*)    :: name
    character(80)   :: file

    write(file,'(a,a)') data_dir(1:lnblnk(data_dir)),name(1:lnblnk(name))
    open(unit=unit,file=file,status='unknown')

    return

  end subroutine open_exist

!**************************************************************************
  subroutine open_unk(unit,name)

    implicit none

    integer(kind=4) :: unit
    character(*)    :: name
    character(80)   :: file

    write(file,'(a,a)') trim(data_dir),trim(name)
    open(unit=unit,file=file,status='unknown')

    return

  end subroutine open_unk

!**************************************************************************
  subroutine open_append(unit,name)

    implicit none

    integer(kind=4) :: unit
    character(*)    :: name    
    character(80)   ::  file

    write(file,'(a,a)') trim(data_dir),trim(name)
    if (numero_step == 1) then 
       open(unit=unit,file=file,status='unknown')
    else
       open(unit=unit,file=file,status='old',position='append')
    end if

    return

  end subroutine open_append

!**************************************************************************
  subroutine write_ts_data

    implicit none

    integer(kind=4) :: i,j
    character(200)   :: file

    i = numero_step
 
    write(file,'(a,a,i3.3)') trim(data_dir),'./halos_results.',i
    open(unit=44,file=file,status='unknown',form='unformatted')
    write(44) i,nb_of_halos(i),nb_of_subhalos(i),aexp(i),nsteps
    ! all units are physical (sizes & pos in Mpc, masses in 10^11 M_sun, speeds in km/s)
    write(44)(liste_halos(j)%level,liste_halos(j)%hosthalo,liste_halos(j)%hostsub,        &
         & liste_halos(j)%nbsub,liste_halos(j)%nextsub,                                   &
         & liste_halos(j)%p%x,liste_halos(j)%p%y, liste_halos(j)%p%z,                     &
         & liste_halos(j)%v%x,liste_halos(j)%v%y,liste_halos(j)%v%z,                      &
         & liste_halos(j)%m,liste_halos(j)%r,liste_halos(j)%spin,                         &
         & liste_halos(j)%sh%a,liste_halos(j)%sh%b,liste_halos(j)%sh%c,                   & 
         & liste_halos(j)%et,liste_halos(j)%ek,liste_halos(j)%ep,                         & 
         & liste_halos(j)%L%x,liste_halos(j)%L%y,liste_halos(j)%L%z,                      &
         & liste_halos(j)%datas%rvir,liste_halos(j)%datas%mvir,                           &
         & liste_halos(j)%datas%tvir,liste_halos(j)%datas%cvel,                           &
         & liste_halos(j)%datas%Reff,                                                     &
         & liste_halos(j)%halo_profile%rho_0,liste_halos(j)%halo_profile%r_c,             &
         & liste_halos(j)%macc,                                                           &
         & j = 1,nb_of_halos(i)+nb_of_subhalos(i))
    close(44)
#ifndef STAR
#ifndef BIG_RUN
    write(file,'(a,a,i3.3)') trim(data_dir),'halos_contam.',i
    open(unit=44,file=file,status='unknown',form='unformatted')
    write(44)(liste_halos(j)%ncont,j=1,nb_of_halos(i)+nb_of_subhalos(i))
    close(44)
#endif 
#endif
    
    return

  end subroutine write_ts_data

!**************************************************************************

  subroutine read_ts_data(i)

    implicit none

    integer(kind=4) :: i,j,k
    real(kind=8)    :: a
    character(200)  :: file
    
    write(file,'(a,a,i3.3)') trim(data_dir),'./halos_results.',i
    open(unit=44,file=file,status='unknown',form='unformatted')
    read(44) i,k,k,a,k
    ! all units are physical (sizes & pos in Mpc, masses in 10^11 M_sun, speeds in km/s)
    read(44)(liste_halos(j)%level,liste_halos(j)%hosthalo,liste_halos(j)%hostsub,         &
         & liste_halos(j)%nbsub,liste_halos(j)%nextsub,                                   &
         & liste_halos(j)%p%x,liste_halos(j)%p%y, liste_halos(j)%p%z,                     &
         & liste_halos(j)%v%x,liste_halos(j)%v%y,liste_halos(j)%v%z,                      &
         & liste_halos(j)%m,liste_halos(j)%r,liste_halos(j)%spin,                         &
         & liste_halos(j)%sh%a,liste_halos(j)%sh%b,liste_halos(j)%sh%c,                   & 
         & liste_halos(j)%et,liste_halos(j)%ek,liste_halos(j)%ep,                         & 
         & liste_halos(j)%L%x,liste_halos(j)%L%y,liste_halos(j)%L%z,                      &
         & liste_halos(j)%datas%rvir,liste_halos(j)%datas%mvir,                           &
         & liste_halos(j)%datas%tvir,liste_halos(j)%datas%cvel,                           &
         & liste_halos(j)%datas%Reff,                                                     &
         & liste_halos(j)%halo_profile%rho_0,liste_halos(j)%halo_profile%r_c,             &
         & liste_halos(j)%macc, &
         & j = 1,nb_of_halos(i)+nb_of_subhalos(i))
    close(44)
#ifndef STAR
#ifndef BIG_RUN
    write(file,'(a,a,i3.3)') trim(data_dir),'halos_contam.',i
    open(unit=44,file=file,status='unknown',form='unformatted')
    read(44)(liste_halos(j)%ncont,j=1,nb_of_halos(i)+nb_of_subhalos(i))
    close(44)
#endif 
#endif    

    return

  end subroutine read_ts_data

!***********************************************************************
  subroutine read_tree_brick
 
  ! read brick files (analysis of individual snapshots of an N-body simulation)
  ! generated with part_to_halo.F

    implicit none

    integer(kind=4)             :: i,j,unitfile,nparts_inhalo,nstep_in,ndummy,imm,nmm
    real(kind=8)                :: masstest
    character(80)               :: filename,file,file2
    integer(kind=4),allocatable :: indices(:)
    real(kind=8), allocatable   :: masstmp(:)

    unitfile = 44

    if (numero_step == 1) then  
       write(file,'(a,a)') trim(data_dir),'./input_TreeMaker.dat'
       open(unit=12,file=file,status='old')
       read(12,*) nstep_in
!!$#ifndef STAR
!!$#ifndef BIG_RUN
!!$       write(file,'(a,a)') trim(data_dir), 'resim_masses.dat'
!!$       open(45,file=file,form='unformatted',status='old')    
!!$       read(45) nbodies
!!$       allocate(mass(nbodies))
!!$       read(45) mass
!!$       close(45)       
!!$#endif
!!$#endif
    endif

    if (numero_step == 1) then  
!!$#ifdef STAR
       mass = 0.0
       write(file,'(a,a)') trim(data_dir), 'resim_masses.dat'
       open(45,file=file,form='unformatted',status='old')    
       read(45) ndummy
       allocate(mass(nbodies),masstmp(ndummy))
       read(45) masstmp
       masstest = minval(masstmp)
       mass = masstest
       close(45)       
       deallocate(masstmp)
!!$#endif
    endif

    read(12,*) file2
    write(filename,'(a,a)') trim(data_dir),trim(file2)
    write(errunit,*) '> reading tree brick file :  ',trim(filename)

    open(unit=unitfile,file=filename,form='unformatted',status='old')

    read(unitfile) nbodies
    read(unitfile) massp
    read(unitfile) aexp(numero_step)
    read(unitfile) omega_t(numero_step)
    read(unitfile) age_univ(numero_step)
    read(unitfile) nb_of_halos(numero_step), nb_of_subhalos(numero_step)
    nmm = 0 ; imm=0
    do i=1,nb_of_halos(numero_step) + nb_of_subhalos(numero_step)
       read(unitfile) nparts_inhalo
       allocate(indices(nparts_inhalo))
       read(unitfile) indices
       do j = 1,nparts_inhalo
          liste_parts(indices(j)) = i
       end do
       deallocate(indices)
       call read_halo(liste_halos(i),unitfile)
       if (nmm .lt. nparts_inhalo) then 
          imm = liste_halos(i)%my_number
          nmm = nparts_inhalo
       endif
    enddo

    print*,'> most massive halo:',imm,' has',nmm,' particles'


    close(unitfile)

    if (numero_step == nsteps) then  
       close(12)
    endif

    return

  end subroutine read_tree_brick

!***********************************************************************
  subroutine read_halo(h,unitfile)

    ! reads halo properties computed by HaloMaker

    implicit none

    integer(kind=4) :: unitfile
    type (halo)     :: h

    read(unitfile) h%my_number
    read(unitfile) h%my_timestep
    ! change the time step number to match the number of branching (time resolution) 
    ! you decided your tree is going to have 
    h%my_timestep = numero_step
    read(unitfile) h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub
    ! Old bug if structure tree does not exist, make, shouldn't occur with the latest version of HaloMaker2.0
    if(h%hosthalo.le.0) h%hosthalo = h%my_number
    if(h%level   .le.0) h%level    = 1
    read(unitfile) h%m
    read(unitfile) h%p%x,h%p%y,h%p%z
    read(unitfile) h%v%x,h%v%y,h%v%z
    read(unitfile) h%L%x,h%L%y,h%L%z 
    read(unitfile) h%r, h%sh%a, h%sh%b, h%sh%c
    read(unitfile) h%ek,h%ep,h%et
    read(unitfile) h%spin
    read(unitfile) h%sigma,h%sigma_bulge,h%m_bulge
    read(unitfile) h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel,h%datas%Reff
    read(unitfile) h%halo_profile%rho_0,h%halo_profile%r_c
    read(unitfile) h%nbin
    read(unitfile) h%rr
    read(unitfile) h%rho

    return

  end subroutine read_halo

!***********************************************************************
  subroutine write_halo(h,st,ih,unitfile,new_halo_number,BushID)

    ! write all the data necessary to build galaxy merger tree and 
    ! compute galaxy properties with GalaxyMaker

    implicit none

    integer(kind=4)             :: unitfile,k,ierr,new_halo_number,st,ih
    type (halo)                 :: h
    integer(kind=4),allocatable :: tabint(:)
    real(kind=8),allocatable    :: tabreal(:)
    integer(kind=4)             :: BushID

    if (new_halo_number == 0) then
       if(h%my_number.ne.ih) h%my_number = ih
       write(unitfile) h%my_number
    else
       write(unitfile) new_halo_number
    endif
    write(unitfile) BushID
    write(unitfile) st
    write(unitfile) h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub
    write(unitfile) h%m
    write(unitfile) h%macc
    write(unitfile) h%p%x,h%p%y,h%p%z
    write(unitfile) h%v%x,h%v%y,h%v%z
    write(unitfile) h%L%x,h%L%y,h%L%z 
    write(unitfile) h%r, h%sh%a, h%sh%b, h%sh%c
    write(unitfile) h%ek,h%ep,h%et
    write(unitfile) h%spin

    write(unitfile) tree(st,ih)%my_fathers%nb_fathers
    if (tree(st,ih)%my_fathers%nb_fathers /= 0) then
       allocate(tabint(tree(st,ih)%my_fathers%nb_fathers),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) '> not enough memory for tabint fathers allocation'
          stop
       endif
       allocate(tabreal(tree(st,ih)%my_fathers%nb_fathers),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) '> not enough memory for tabreal fathers allocation '
          stop
       endif
       do k = 1,tree(st,ih)%my_fathers%nb_fathers
          tabint(k) = tree(st,ih)%my_fathers%list_fathers(k)
       end do
       write(unitfile) tabint(1:tree(st,ih)%my_fathers%nb_fathers)
       do k = 1,tree(st,ih)%my_fathers%nb_fathers
          tabreal(k) = tree(st,ih)%my_fathers%mass_fathers(k)
       end do
       write(unitfile) tabreal(1:tree(st,ih)%my_fathers%nb_fathers)
       deallocate(tabint,tabreal)
    endif

    write(unitfile) tree(st,ih)%my_sons%nb_sons
    if (tree(st,ih)%my_sons%nb_sons /= 0) then
       allocate(tabint(tree(st,ih)%my_sons%nb_sons),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) '> not enough memory for tabint sons allocation '
          stop
       endif
       do k = 1,tree(st,ih)%my_sons%nb_sons
          tabint(k) = tree(st,ih)%my_sons%list_sons(k)
       end do
       write(unitfile) tabint(1:tree(st,ih)%my_sons%nb_sons)     
       write(unitfile) tree(st,ih)%my_sons%main_son
       deallocate(tabint)
    endif

    write(unitfile) h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel,h%datas%Reff
    write(unitfile) h%halo_profile%rho_0,h%halo_profile%r_c
#ifndef BIG_RUN
    write(unitfile) h%ncont
#endif

    return

  end subroutine write_halo

!**************************************************************************

  subroutine write_single_tree

    ! in fact, write all the trees in a single tree.dat file 

    implicit none

    integer(kind=4) :: j,st,unitfile,lnblnk
    character(80)   :: treefile

    write(treefile,'(a,a)') data_dir(1:lnblnk(data_dir)),'tree.dat'
    unitfile = 43

    write(errunit,*) '> Writing trees in old format in file:',trim(treefile)

    open(unit=unitfile,file=treefile,status='unknown',form='unformatted') 

    write(unitfile) nsteps
    write(unitfile) nb_of_halos(1:nsteps), nb_of_subhalos(1:nsteps)
    write(unitfile) aexp(1:nsteps)
    print*,aexp(1:nsteps)
    write(unitfile) omega_t(1:nsteps)
    write(unitfile) age_univ(1:nsteps)

    ! no renumbering of halos
    do st = 1,nsteps
       call read_ts_data(st)
       do j = 1,nb_of_halos(st)+nb_of_subhalos(st)
          call write_halo(liste_halos(j),st,j,unitfile,0,0)
       end do
    end do

    close(unitfile)

    return

  end subroutine write_single_tree

!***********************************************************************
  subroutine write_bushes
    
    ! JB -
    ! write bushes to a number of files trying to balance the number of halos per file.
    ! NB: n_tree_files is interpreted as a number of desired files.
    
    implicit none
    
    integer(kind=4)             :: ts, ih, NBushes, bush, nwish, n, ifile
    integer(kind=4),allocatable :: NHInBush(:),BushToFile(:)
    integer(kind=4),allocatable :: NInFile(:,:,:)
    integer(kind=4)             :: wcnt

    ! Count bushes and halos per bushes 
    NBushes  = maxval(tree(:,:)%BushID)
    allocate(NHInBush(NBushes))
    NHInBush = 0
    do ts = 1,nsteps
       do ih = 1, nb_of_halos(ts) + nb_of_subhalos(ts)
          bush = tree(ts,ih)%BushID
          if (bush > 0) NHInBush(bush) = NHInBush(bush) + 1
       end do
    end do
    
    ! try to balance 
    allocate(BushToFile(NBushes))
    nwish = int(real(sum(NHInBush)) / real(n_tree_files))
    n     = 0
    ifile = 1
    do bush = 1,NBushes
       BushToFile(bush) = ifile
       n = n + NHInBush(bush)
       if (n >= nwish) then 
          ifile = ifile + 1
          n = 0
       end if
    end do
    ! reduce number of files to what is needed.
    n_tree_files = ifile

    ! compute nb of (sub)halos per file per ts
    allocate(NInFile(n_tree_files,nsteps,2))
    NInFile(:,:,:) = 0
    do ts = 1,nsteps
       do ih = 1, nb_of_halos(ts)
          bush  = tree(ts,ih)%BushID
          if (bush > 0) then 
             ifile = BushToFile(bush)
             NInFile(ifile,ts,1) = NInFile(ifile,ts,1) + 1
          end if
       end do
       do ih = nb_of_halos(ts)+1, nb_of_halos(ts) + nb_of_subhalos(ts)
          bush  = tree(ts,ih)%BushID
          if (bush > 0) then 
             ifile = BushToFile(bush)
             NInFile(ifile,ts,2) = NInFile(ifile,ts,2) + 1
          end if
       end do
    end do

    ! do dump the stuff 
    wcnt = 0
    call open_all_tree_files(NInFile)
    do ts = 1,nsteps
       call read_ts_data(ts)
       do ih = 1, nb_of_halos(ts) + nb_of_subhalos(ts)
          bush = tree(ts,ih)%BushID
          if (bush > 0) then 
             ifile = BushToFile(bush) + dumptree_unit
             if (BushToFile(bush) == 1) write(132,*) 'halo ',ih,' bushID ',bush
             call write_halo(liste_halos(ih),ts,ih,ifile,0,bush)
             wcnt = wcnt + 1
          end if
       end do
    end do
    call close_all_tree_files

    if (wcnt .ne. sum(NHInBush)) then 
       write(errunit,*) 'ooops, missed some halos'
       stop
    end if

    deallocate(NHInBush,BushToFile,NInFile)

    return

  end subroutine write_bushes

!***********************************************************************
  subroutine open_all_tree_files(NInFile)
    
    implicit none
    
    integer(kind=4)        :: unit,ifile
    integer(kind=4)        :: NInFile(n_tree_files,nsteps,2)
    character(200)         :: file
    
    do ifile = 1, n_tree_files
       unit = dumptree_unit + ifile
       write(file,'(a,a,i3.3,a,i3.3)') trim(data_dir),'./tree_file_',nsteps,'.',ifile
       open(unit=unit,file=file,status="unknown",form="unformatted")
       ! write header 
       write(unit) nsteps    
       write(unit) NInFile(ifile,1:nsteps,1),NInFile(ifile,1:nsteps,2)
       write(unit) aexp(1:nsteps)
       write(unit) omega_t(1:nsteps)
       write(unit) age_univ(1:nsteps)
    end do
    
    return
    
  end subroutine open_all_tree_files
  
!***********************************************************************
  subroutine close_all_tree_files
    
    implicit none
    
    integer(kind=4) :: unit,ifile
    
    do ifile = 1, n_tree_files
       unit = dumptree_unit + ifile
       close(unit)
    end do
    
    return
    
  end subroutine close_all_tree_files

!***********************************************************************
  
  subroutine write_one_bush(bushid)
    
    implicit none
    
    integer(kind=4) :: bushid
    integer(kind=4) :: unit,ih,ts
    integer(kind=4),allocatable :: NInFile(:,:)
    character(200)              :: file

    ! count nb of halos in bush with bushid
    allocate(NInFile(nsteps,2))
    NInFile = 0
    do ts = 1,nsteps
       do ih = 1, nb_of_halos(ts)
          if (tree(ts,ih)%BushID == bushid) NInFile(ts,1) = NInFile(ts,1) + 1
       end do
       do ih = nb_of_halos(ts)+1, nb_of_halos(ts) + nb_of_subhalos(ts)
          if (tree(ts,ih)%BushID == bushid) NInFile(ts,2) = NInFile(ts,2) + 1
       end do
    end do
    
    ! write selected bush to a single file
    unit = dumptree_unit 
    write(file,'(a,a,i3.3,a,i8.8)') trim(data_dir),'./tree_file_',nsteps,'.',bushid
    open(unit=unit,file=file,status="unknown",form="unformatted")
    ! write header 
    write(unit) nsteps    
    write(unit) NInFile(1:nsteps,1),NInFile(1:nsteps,2)
    write(unit) aexp(1:nsteps)
    write(unit) omega_t(1:nsteps)
    write(unit) age_univ(1:nsteps)
    ! write each halo
    do ts = 1,nsteps
       call read_ts_data(ts)
       do ih = 1, nb_of_halos(ts) + nb_of_subhalos(ts)
          if (tree(ts,ih)%BushID == bushid) call write_halo(liste_halos(ih),ts,ih,unit,0,bushid)
       end do
    end do
    close(unit)

    return

  end subroutine write_one_bush

!***********************************************************************
!//////////////////////////////////////////////////////////////////////////

end module input_output



