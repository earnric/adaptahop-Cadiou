module input_output

  use tree_defs

  public

contains

!//////////////////////////////////////////////////////////////////////////
!**************************************************************************
  subroutine get_max_n_parts

    implicit none
    integer(kind=4) :: ierr
    integer(kind=4) :: i,number_of_halos,number_of_subhalos,number_of_particles,idid
    real(kind=4)    :: dummy
    character(200)   :: file,file2
    
    ! we now take into account that some tree_brick specified 
    ! in the input_TreeMaker file might not exist for some reason 
    ! or other, we count the steps where we really have tree_bricks

    write(file,'(a,a)') trim(data_dir),'input_TreeMaker.dat'

    open(unit=12,file=file,status='old')
    read(12,*) nsteps_do, n_tree_files

    max_groups = 0
    nbodies    = 0
    nsteps     = nsteps_do
    idid       = 0
    do i = 1,nsteps    
       read(12,*) file2
       write(file,'(a)') trim(file2)
       if(file(1:1) .ne."/") then
          write(file,'(a,a)') trim(data_dir),trim(file2) 
       end if
       open(unit=13,file=file,form='unformatted',status='old',iostat=ierr)
       if(ierr.eq.0) then
          read(13) number_of_particles
          read(13) dummy
          read(13) dummy
          read(13) dummy
          read(13) dummy
          read(13) number_of_halos, number_of_subhalos
          close(13)
          if (idid .gt. 0 .and. nbodies /= number_of_particles) then 
#ifdef BIG_RUN
             write(errunit,*) '> pb in get_max_n_parts: nb of particles cannot be changed between time steps'
             stop
#else
             ! leo: it could happen for zoom simulations where haloes have been identified only in a sub-volume
             write(errunit,*) '> WARNING in get_max_n_parts: nb of particles has changed between time steps'
#endif
          end if
          nbodies = number_of_particles
          if (max_groups .lt. (number_of_halos + number_of_subhalos )) then
             max_groups = number_of_halos + number_of_subhalos
          end if
          idid = idid + 1
       else

          write(errunit,*) '> missing tree_brick:',trim(file)
          ! tree_brick missing have to do without
          nsteps = nsteps - 1
       end if
    end do
    
    close(12)
    
    if(nsteps.le.1) then
       write(errunit,*) '> not enough tree_brick to analyse'
       stop
    end if

    return

  end subroutine get_max_n_parts

!***********************************************************************
  subroutine open_exist(unit,name)

    implicit none

    integer(kind=4) :: unit,lnblnk
    character(*)    :: name
    character(200)   :: file

    write(file,'(a,a)') data_dir(1:lnblnk(data_dir)),name(1:lnblnk(name))
    open(unit=unit,file=file,status='unknown')

    return

  end subroutine open_exist

!**************************************************************************
  subroutine open_unk(unit,name)

    implicit none

    integer(kind=4) :: unit
    character(*)    :: name
    character(200)   :: file

    write(file,'(a,a)') trim(data_dir),trim(name)
    open(unit=unit,file=file,status='unknown')

    return

  end subroutine open_unk

!**************************************************************************
  subroutine open_append(unit,name)

    implicit none

    integer(kind=4) :: unit
    character(*)    :: name    
    character(200)   ::  file

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
 
    write(file,'(a,a,i3.3)') trim(data_dir),'halos_results.',i
    open(unit=44,file=file,status='unknown',form='unformatted')
    write(44) i,nb_of_halos(i),nb_of_subhalos(i),aexp(i),nsteps
    ! all units are physical (sizes & pos in Mpc, masses in 10^11 M_sun, speeds in km/s)
#ifndef SIMPL
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
         & liste_halos(j)%halo_profile%rho_0,liste_halos(j)%halo_profile%r_c,             &
         & liste_halos(j)%macc,                                                           &
         & j = 1,nb_of_halos(i)+nb_of_subhalos(i))
#endif
#ifdef SIMPL
   write(44)(liste_halos(j)%p%x,liste_halos(j)%p%y, liste_halos(j)%p%z,                   &
         & liste_halos(j)%v%x,liste_halos(j)%v%y,liste_halos(j)%v%z,                      &
         & liste_halos(j)%m,liste_halos(j)%r,liste_halos(j)%spin,                         &
         & liste_halos(j)%sh%a,liste_halos(j)%sh%b,liste_halos(j)%sh%c,                   & 
         & liste_halos(j)%et,liste_halos(j)%ek,liste_halos(j)%ep,                         & 
         & liste_halos(j)%L%x,liste_halos(j)%L%y,liste_halos(j)%L%z,                      &
         & liste_halos(j)%datas%rvir,liste_halos(j)%datas%mvir,                           &
         & liste_halos(j)%datas%tvir,liste_halos(j)%datas%cvel,                           &
         & liste_halos(j)%halo_profile%rho_0,liste_halos(j)%halo_profile%r_c,             &
         & liste_halos(j)%macc,                                                           &
         & j = 1,nb_of_halos(i)+nb_of_subhalos(i))  
#endif
    close(44)
#ifndef BIG_RUN
    write(file,'(a,a,i3.3)') trim(data_dir),'halos_contam.',i
    open(unit=44,file=file,status='unknown',form='unformatted')
    write(44)(liste_halos(j)%ncont,j=1,nb_of_halos(i)+nb_of_subhalos(i))
    close(44)
#endif 
    
    return

  end subroutine write_ts_data

!**************************************************************************

  subroutine read_ts_data(i)

    implicit none

    integer(kind=4) :: i,j,k
    real(kind=4)    :: a
    character(200)  :: file
    
    write(file,'(a,a,i3.3)') trim(data_dir),'halos_results.',i
    open(unit=44,file=file,status='unknown',form='unformatted')
    read(44) i,k,k,a,k
    ! all units are physical (sizes & pos in Mpc, masses in 10^11 M_sun, speeds in km/s)
#ifndef SIMPL
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
         & liste_halos(j)%halo_profile%rho_0,liste_halos(j)%halo_profile%r_c,             &
         & liste_halos(j)%macc, &
         & j = 1,nb_of_halos(i)+nb_of_subhalos(i))
#endif
#ifdef SIMPL
    read(44)(liste_halos(j)%p%x,liste_halos(j)%p%y, liste_halos(j)%p%z,                   &
         & liste_halos(j)%v%x,liste_halos(j)%v%y,liste_halos(j)%v%z,                      &
         & liste_halos(j)%m,liste_halos(j)%r,liste_halos(j)%spin,                         &
         & liste_halos(j)%sh%a,liste_halos(j)%sh%b,liste_halos(j)%sh%c,                   & 
         & liste_halos(j)%et,liste_halos(j)%ek,liste_halos(j)%ep,                         & 
         & liste_halos(j)%L%x,liste_halos(j)%L%y,liste_halos(j)%L%z,                      &
         & liste_halos(j)%datas%rvir,liste_halos(j)%datas%mvir,                           &
         & liste_halos(j)%datas%tvir,liste_halos(j)%datas%cvel,                           &
         & liste_halos(j)%halo_profile%rho_0,liste_halos(j)%halo_profile%r_c,             &
         & liste_halos(j)%macc, &
         & j = 1,nb_of_halos(i)+nb_of_subhalos(i))
#endif
    close(44)
#ifndef BIG_RUN
    write(file,'(a,a,i3.3)') trim(data_dir),'halos_contam.',i
    open(unit=44,file=file,status='unknown',form='unformatted')
    read(44)(liste_halos(j)%ncont,j=1,nb_of_halos(i)+nb_of_subhalos(i))
    close(44)
#endif 
    
    return

  end subroutine read_ts_data

!***********************************************************************
  subroutine read_tree_brick(ierr)
 
  ! read brick files (analysis of individual snapshots of an N-body simulation)
  ! generated with part_to_halo.F

    implicit none
    integer(kind=4)             :: ierr
    integer(kind=4)             :: i,j,unitfile,nparts_inhalo,nstep_in,ndum
    character(200)               :: filename,file,file2
    integer(kind=4),allocatable :: indices(:)
#ifndef BIG_RUN
    integer(kind=4)             :: nbodies_full_sim
#endif

    unitfile = 44

    if (st_do == 1) then  
       write(file,'(a,a)') trim(data_dir),'input_TreeMaker.dat'
       open(unit=12,file=file,status='old')
       read(12,*) nstep_in,ndum
#ifndef BIG_RUN
       write(file,'(a,a)') trim(data_dir), 'resim_masses.dat'
       open(45,file=file,form='unformatted',status='old')    
       read(45) nbodies
       allocate(mass(nbodies))
       read(45) mass
       close(45)       
#endif
    endif

    read(12,*) file2
    write(filename,'(a)') trim(file2)
    if(filename(1:1) .ne."/") then
       write(filename,'(a,a)') trim(data_dir),trim(file2)
    end if

    open(unit=unitfile,file=filename,form='unformatted',status='old',iostat=ierr)
    if(ierr.eq.0) then
       write(errunit,*) '> reading tree brick file :  ',trim(filename)
       ! if file exist read and save data in new step
       numero_step = numero_step + 1
       if(numero_step.gt.nsteps) then
          write(errunit,*) '> Error in read_tree_brick'
          write(errunit,*) '> numero_step,nsteps:',numero_step,nsteps
          stop
       end if
#ifndef BIG_RUN
       read(unitfile) nbodies_full_sim
#else
       read(unitfile) nbodies
#endif
       read(unitfile) massp
       read(unitfile) aexp(numero_step)
       read(unitfile) omega_t(numero_step)
       read(unitfile) age_univ(numero_step)
       read(unitfile) nb_of_halos(numero_step), nb_of_subhalos(numero_step)

       do i=1,nb_of_halos(numero_step) + nb_of_subhalos(numero_step)
          read(unitfile) nparts_inhalo
          allocate(indices(nparts_inhalo))
          read(unitfile) indices
          do j = 1,nparts_inhalo
             liste_parts(indices(j)) = i
          end do
          deallocate(indices)
#ifndef SIMPL
          call read_halo(liste_halos(i),unitfile)
#endif
#ifdef SIMPL
          call read_halo(liste_halos(i),tree(numero_step,i),unitfile)
#endif
       enddo
       close(unitfile)
    else
       write(errunit,*) '> missing tree_brick:',trim(filename)
    endif
    if (st_do == nsteps_do) then  
       close(12)
    endif

    return

  end subroutine read_tree_brick

!***********************************************************************
#ifndef SIMPL
  subroutine read_halo(h,unitfile)
#endif
#ifdef SIMPL
  subroutine read_halo(h,t,unitfile)
#endif

    ! reads halo properties computed by HaloMaker

    implicit none

    integer(kind=4)   :: unitfile
    type (halo)       :: h
#ifdef SIMPL
    type(tree_struct) :: t
    integer(kind=4)   :: nbsub
#endif
    read(unitfile) h%my_number
    read(unitfile) h%my_timestep
    ! change the time step number to match the number of branching (time resolution) 
    ! you decided your tree is going to have 
    h%my_timestep = numero_step
#ifndef SIMPL
    read(unitfile) h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub
    ! Old bug if structure tree does not exist, make, shouldn't occur with the latest version of HaloMaker2.0
    if(h%hosthalo.le.0) h%hosthalo = h%my_number
    if(h%level   .le.0) h%level    = 1
#endif
#ifdef SIMPL
    read(unitfile) t%level,t%hosthalo,t%hostsub,nbsub,t%nextsub
    ! Old bug if structure tree does not exist, make, shouldn't occur with the latest version of HaloMaker2.0
    if(t%hosthalo.le.0) t%hosthalo = h%my_number
    if(t%level   .le.0) t%level    = 1
#endif
    read(unitfile) h%m
    read(unitfile) h%p%x,h%p%y,h%p%z
    read(unitfile) h%v%x,h%v%y,h%v%z
    read(unitfile) h%L%x,h%L%y,h%L%z 
    read(unitfile) h%r, h%sh%a, h%sh%b, h%sh%c
    read(unitfile) h%ek,h%ep,h%et
    read(unitfile) h%spin
    read(unitfile) h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel
    read(unitfile) h%halo_profile%rho_0,h%halo_profile%r_c
#ifdef CONTAM
    read(unitfile) h%contamination
#endif

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
    real(kind=4),allocatable    :: tabreal(:)
    integer(kind=4)             :: BushID

    if (new_halo_number == 0) then
       if(h%my_number.ne.ih) h%my_number = ih
       write(unitfile) h%my_number
    else
       write(unitfile) new_halo_number
    endif
    write(unitfile) BushID
    write(unitfile) st
#ifndef SIMPL
    write(unitfile) h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub
#endif
    write(unitfile) h%m
    write(unitfile) h%macc
    write(unitfile) h%p%x,h%p%y,h%p%z
    write(unitfile) h%v%x,h%v%y,h%v%z
    write(unitfile) h%L%x,h%L%y,h%L%z 
    write(unitfile) h%r, h%sh%a, h%sh%b, h%sh%c
    write(unitfile) h%ek,h%ep,h%et
    write(unitfile) h%spin
#ifndef SIMPL
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
       !! JB : a simpler solution would be : 
       !!write(unitfile) (tree(st,ih)%my_fathers%list_fathers(k),k=1,tree(st,ih)%my_fathers%nb_fathers)
       do k = 1,tree(st,ih)%my_fathers%nb_fathers
          tabreal(k) = real(tree(st,ih)%my_fathers%mass_fathers(k),4)
       end do
       write(unitfile) tabreal(1:tree(st,ih)%my_fathers%nb_fathers)
       deallocate(tabint,tabreal)
    endif
#endif

#ifndef SIMPL  
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
#endif

    write(unitfile) h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel
    write(unitfile) h%halo_profile%rho_0,h%halo_profile%r_c
#ifndef BIG_RUN
    write(unitfile) h%ncont
#endif

#ifdef SIMPL    
    write(unitfile) tree(st,ih)%level,tree(st,ih)%hosthalo,tree(st,ih)%hostsub,tree(st,ih)%nextsub
    write(unitfile) tree(st,ih)%my_dads%nb_dads
    if (tree(st,ih)%my_dads%nb_dads /= 0) then
       allocate(tabint(tree(st,ih)%my_dads%nb_dads),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) '> not enough memory for tabint dads allocation'
          stop
       endif
       allocate(tabreal(tree(st,ih)%my_dads%nb_dads),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) '> not enough memory for tabreal dads allocation '
          stop
       endif
       do k = 1,tree(st,ih)%my_dads%nb_dads
          tabint(k)  = tree(st,ih)%my_dads%list_dads(k)
          tabreal(k) = tree(st,ih)%my_dads%mass_dads(k)/h%m
          if(tabreal(k).gt.(1.+1e-5)) then
             write(errunit,*) '> Error in write_halo for halo:',st,ih
             write(errunit,*) '> idad,mass_dad:',tabint(k),tabreal(k)
             write(errunit,*) '> halo mass, mass transmitted'
             write(errunit,'(a,2(1x,E13.5))') ' >',&
                  liste_halos(ih)%m,tree(st,ih)%my_dads%mass_dads(k)
             write(errunit,*) '> mass fraction obtained from progenitor is greeter than 1'
             stop
          else if(tabreal(k).gt.1.) then
             tabreal(k) = 1.
          end if
       end do
       write(unitfile) tabint(1:tree(st,ih)%my_dads%nb_dads)
       write(unitfile) tabreal(1:tree(st,ih)%my_dads%nb_dads)
       deallocate(tabint,tabreal)
    endif
    write(unitfile) tree(st,ih)%my_sons%main_son
    write(unitfile) tree(st,ih)%frag
#endif  

    return

  end subroutine write_halo

  !**************************************************************************
  subroutine write_single_tree

    ! in fact, write all the trees in a single tree.dat file 

    implicit none

    integer(kind=4) :: j,st,unitfile,lnblnk
    character(200)   :: treefile

    write(treefile,'(a,a)') data_dir(1:lnblnk(data_dir)),'tree.dat'
    unitfile = 43

    write(errunit,*) '> Writing trees in old format in file:',trim(treefile)

    open(unit=unitfile,file=treefile,status='unknown',form='unformatted') 

    write(unitfile) nsteps
    write(unitfile) nb_of_halos(1:nsteps), nb_of_subhalos(1:nsteps)
    write(unitfile) aexp(1:nsteps)
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
    integer(kind=4),allocatable :: NOutFile(:,:),FragInFile(:,:,:)
    integer(kind=4)             :: wcnt
    character(200)              :: file_log,file_list,file
    
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
    if(n_tree_files.gt.1) then 
       nwish = int(real(sum(NHInBush)) / real(n_tree_files))
    else
       ! we might want only one tree_file so make sure there is only one
       nwish = sum(NHInBush)
    end if
    n     = 0
    ifile = 1
    do bush = 1,NBushes
       BushToFile(bush) = ifile
       n = n + NHInBush(bush)
       if (n > nwish) then 
          ifile = ifile + 1
          n = 0
       end if
    end do
    ! reduce number of files to what is needed.
    n_tree_files = ifile

    ! compute nb of (sub)halos per file per ts
    allocate(NInFile(n_tree_files,nsteps,2))
    allocate(FragInFile(n_tree_files,nsteps,2))
    allocate(NOutFile(nsteps,3))
    NInFile(:,:,:)    = 0
    FragInFile(:,:,:) = 0
    NOutFile(:,:)     = 0
    do ts = 1,nsteps
       do ih = 1, nb_of_halos(ts)
          bush  = tree(ts,ih)%BushID
          if (bush > 0) then 
             ifile = BushToFile(bush)
             NInFile(ifile,ts,1) = NInFile(ifile,ts,1) + 1
#ifdef SIMPL             
             if(tree(ts,ih)%frag.eq.1)  FragInFile(ifile,ts,1) = FragInFile(ifile,ts,1) + 1
#endif
          else
             NOutFile(ts,1) = NOutFile(ts,1) + 1
          end if
       end do
       do ih = nb_of_halos(ts)+1, nb_of_halos(ts) + nb_of_subhalos(ts)
          bush  = tree(ts,ih)%BushID
          if (bush > 0) then 
             ifile = BushToFile(bush)
             NInFile(ifile,ts,2) = NInFile(ifile,ts,2) + 1
#ifdef SIMPL             
             if(tree(ts,ih)%frag.eq.1)  FragInFile(ifile,ts,2) = FragInFile(ifile,ts,2) + 1
          else
             if(tree(ts,ih)%level.gt.1) then
                NOutFile(ts,2) = NOutFile(ts,2) + 1
             else
                NOutFile(ts,3) = NOutFile(ts,3) + 1
             end if
          end if
#else
          else
             NOutFile(ts,2) = NOutFile(ts,2) + 1
          endif
#endif
       end do
    end do

    ! do dump the stuff 
    wcnt = 0
    call open_all_tree_files(NInFile)
    do ts = 1,nsteps
       call read_ts_data(ts)
#ifdef SIMPL
#ifdef R_DADLESS_SUB
       call recover_masses(ts)
#else
#ifdef COL_TREE
       call recover_masses(ts)
#endif
#endif   
#endif
       do ih = 1, nb_of_halos(ts) + nb_of_subhalos(ts)
          bush = tree(ts,ih)%BushID
          if (bush > 0) then 
             ifile = BushToFile(bush) + dumptree_unit
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

    !! DT: lets write log file to know what is outputted and what is not.
    write(file_log,*) trim(data_dir),"tree_files.log" 
    open(unit=42,form='formatted',status='unknown',file=file_log)
    write(42,'(a20,i3)') ' > nb of tree_files:',n_tree_files
    write(42,'(a20,i3)') ' > nb of time steps:',nsteps
    do ifile = 1, n_tree_files
       write(42,*)
       write(42,'(a,i3.3,a1,i3.3,a1)') &
            ' > file ''tree_file_',nsteps,'.',ifile,''''
       write(42,*) '> st, nb of halos,nb of subhalos, nb of halo frag, nb of subhalo frag'
       do ts = 1,nsteps
          write(42,'(a5,i3,4(1x,i10))') ' >   ',ts,NInFile(ifile,ts,1:2),FragInFile(ifile,ts,1:2)
       end do
       write(42,'(a8,4(1x,i10))') ' > total',sum(NInFile(ifile,1:nsteps,1)),sum(NInFile(ifile,1:nsteps,2)),&
            sum(FragInFile(ifile,1:nsteps,1)),sum(FragInFile(ifile,1:nsteps,2))
    end do
    write(42,*)
    write(42,*) '> nb of elements not outputed in any tree_file_'
    write(42,*) '> st, nb of halos,nb of subhalos,nb of removed subhalos'
    do ts = 1,nsteps
       write(42,'(a5,i3,3(1x,i10))') ' >   ',ts,NOutFile(ts,1:3)
    end do
    write(42,'(a8,3(1x,i10))') ' > total',sum(NOutFile(1:nsteps,1)),sum(NOutFile(1:nsteps,2)),sum(NOutFile(1:nsteps,3))
    close(42)
    
    !! DT: lets write a list of tree_files 
    write(file_list,*) trim(data_dir),"tree_files.list" 
    open(unit=43,form='formatted',status='unknown',file=file_list)
    write(43,*) n_tree_files
    do ifile = 1,n_tree_files
       write(file,'(a,i3.3,a,i3.3)')'tree_file_',nsteps,'.',ifile
       write(43,*) trim(file)
    end do
    close(43)

    deallocate(FragInFile,NOutFile)
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
       write(file,'(a,a,i3.3,a,i3.3)') trim(data_dir),'tree_file_',nsteps,'.',ifile
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
    write(file,'(a,a,i3.3,a,i8.8)') trim(data_dir),'tree_file_',nsteps,'.',bushid
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
#ifdef SIMPL
#ifdef R_DADLESS_SUB
       call recover_masses(ts)
#endif
#ifdef COL_TREE
       call recover_masses(ts)
#endif   
#endif
       do ih = 1, nb_of_halos(ts) + nb_of_subhalos(ts)
          if (tree(ts,ih)%BushID == bushid) call write_halo(liste_halos(ih),ts,ih,unit,0,bushid)
       end do
    end do
    close(unit)

    return

  end subroutine write_one_bush

#ifdef SIMPL
  ! This routine is only neaded when #ifdef R_DADLESS_SUB or #ifdef COL_TREE
  !***********************************************************************
  subroutine recover_masses(st)
    
    implicit none
    integer(kind=4) :: st
    integer(kind=4) :: i,ihost
    integer(kind=4) :: stbug,ihbug,n_removed
    

    if(nb_of_subhalos(st).le.0) return

    stbug = -1
    ihbug = -1

    if(st.eq.stbug.and.ihbug.gt.0) then
       write(errunit,*) '> halo:',stbug,ihbug
       write(errunit,*) '> before recover masses:',liste_halos(ihbug)%m
    end if
    n_removed = 0

    do i = nb_of_halos(st)+nb_of_subhalos(st),nb_of_halos(st)+1,-1
       if(tree(st,i)%level.le.0) then
          ihost = tree(st,i)%hostsub
          if(ihost.le.0) then 
             write(errunit,*) '> fatal error in recover_masses of:',st,i
             write(errunit,*) '> tree(st,i)%hostsub:',tree(st,i)%hostsub
             stop
          end if
          if(tree(st,ihost)%level.le.0) then
             write(errunit,*) '> fatal error in recover_masses of:',st,i
             write(errunit,*) '> ihost=tree(st,i)%hostsub:',tree(st,i)%hostsub
             write(errunit,*) '> tree(st,ihost)%level:',tree(st,ihost)%level
             stop
          end if
          if(st.eq.stbug.and.ihost.eq.ihbug) then
             write(errunit,'(1x,a,i3,1x,i6,E14.5)') '> adding mass:',st,i,liste_halos(i)%m*1e11
             write(errunit,'(1x,a,E14.5)') '> old mass:',liste_halos(ihost)%m*1.e11
          end if
          liste_halos(ihost)%m = liste_halos(ihost)%m + liste_halos(i)%m
          if(st.eq.stbug.and.ihost.eq.ihbug) then
             write(errunit,'(1x,a,E14.5)') '> new mass:',liste_halos(ihost)%m*1.e11
          end if
          if(st.eq.stbug) n_removed = n_removed + 1
       end if
    end do
    if(st.eq.stbug) then
       if(ihbug .gt. 0) write(errunit,'(1x,a,E14.5)') '> after recover masses:',liste_halos(ihbug)%m
       write(errunit,*) '> nb of subhalos removed:',n_removed
    end if
    
    return
    
  end subroutine recover_masses
#endif

!***********************************************************************
!//////////////////////////////////////////////////////////////////////////

end module input_output



