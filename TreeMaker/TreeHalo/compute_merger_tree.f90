module compute_tree

  use tree_defs
  use input_output

  public

contains 

!*************************************************************************
  subroutine init

    implicit none

    integer(kind=4) :: ierr

    ! determine maximum number of groups (over all timesteps) and number of particles
    call get_max_n_parts
    
    ! allocation of particle lists
    allocate(liste_parts(nbodies),ex_liste_parts(nbodies),linked_list(0:nbodies+1),stat=ierr) 
    if (ierr /= 0) then
       write(errunit,*) '> not enough memory for particle lists allocation'
       stop
    endif
    liste_parts = 0 ; ex_liste_parts = 0 ; linked_list = 0

    allocate(MassToAccrete(nbodies)) ! will be initialized later in new_step (when numero_step==1)

    ! allocate and initialise liste_halos which contains all halos in two consecutive timesteps (1 is earlier than 2...)
    allocate(liste_halos(max_groups),stat=ierr)  
    if (ierr /= 0) then
       write(errunit,*) '> not enough memory for liste_halos allocation'
       stop
    endif
    call init_liste_halos

    ! allocate and initialise tree structure for all halos at all timesteps
    allocate(tree(nsteps,max_groups),stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) '> not enough memory for tree allocation'
       stop
    endif
    call init_tree 

    ! allocation of global variables containing data of brick files 
    allocate(nb_of_halos(nsteps),nb_of_subhalos(nsteps),aexp(nsteps),age_univ(nsteps),omega_t(nsteps),stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) '> not enough memory for brick file global variables'
       stop
    endif
    nb_of_halos = 0 ; nb_of_subhalos= 0; aexp = 0.0 ; age_univ = 0.0 ; omega_t = 0.0

    write(errunit,*) 
    write(errunit,*) '> Technical parameter values : '
    write(errunit,*) '> ---------------------------  '
    write(errunit,*) 
    write(errunit,*) '> number of time steps                           : ',nsteps
    write(errunit,*) '> maximum number of groups (over all time steps) : ',max_groups
    write(errunit,*) '> number of particles in the original N-body sim : ',nbodies
    write(errunit,*) 

    return

  end subroutine init

!*************************************************************************
  subroutine init_tree
    
    implicit none
    
    integer(kind=4) :: ts,hno 
    
    do ts = 1,nsteps 
       do hno = 1,max_groups
#ifdef SIMPL
          ! halo structure tree
          tree(ts,hno)%hosthalo                 = 0
          tree(ts,hno)%hostsub                  = 0
          tree(ts,hno)%nextsub                  = -1
          tree(ts,hno)%level                    = 0
          tree(ts,hno)%frag                     = 0
#endif
#ifndef SIMPL
          ! halo structure tree
          tree(ts,hno)%host                     = 0
          tree(ts,hno)%nextsub                  = -1
#endif
          ! father type
          tree(ts,hno)%my_fathers%nb_fathers    = 0
          nullify(tree(ts,hno)%my_fathers%list_fathers)
          nullify(tree(ts,hno)%my_fathers%mass_fathers)
          ! son type
          tree(ts,hno)%my_sons%nb_sons          = 0
          nullify(tree(ts,hno)%my_sons%list_sons)
#ifdef SIMPL
          tree(ts,hno)%my_sons%main_son         = -1
          ! dads type
          tree(ts,hno)%my_dads%nb_dads          = 0
          nullify(tree(ts,hno)%my_dads%list_dads)       
          nullify(tree(ts,hno)%my_dads%mass_dads)       
#endif
          ! BushID
          tree(ts,hno)%BushID                   = -1
       end do
    end do
    
    return 
    
  end subroutine init_tree

!*************************************************************************
  subroutine init_liste_halos
    
    implicit none
      
    integer(kind=4) :: hno 
    
    do hno = 1,max_groups
       ! baryon type
       liste_halos(hno)%datas%rvir               = 0.0
       liste_halos(hno)%datas%mvir               = 0.0  
       liste_halos(hno)%datas%tvir               = 0.0
       liste_halos(hno)%datas%cvel               = 0.0
       ! shape type
       liste_halos(hno)%sh%a                     = 0.0
       liste_halos(hno)%sh%b                     = 0.0
       liste_halos(hno)%sh%c                     = 0.0
       ! vector type
       liste_halos(hno)%p%x                      = 0.0
       liste_halos(hno)%p%y                      = 0.0
       liste_halos(hno)%p%z                      = 0.0
       liste_halos(hno)%v%x                      = 0.0 
       liste_halos(hno)%v%y                      = 0.0 
       liste_halos(hno)%v%z                      = 0.0
       liste_halos(hno)%L%x                      = 0.0 
       liste_halos(hno)%L%y                      = 0.0 
       liste_halos(hno)%L%z                      = 0.0
       ! profile type
       liste_halos(hno)%halo_profile%rho_0       = 0.0
       liste_halos(hno)%halo_profile%r_c         = 0.0  
       ! standard type
       liste_halos(hno)%my_number                = 0
       liste_halos(hno)%my_timestep              = 0
       liste_halos(hno)%m                        = 0.0
       liste_halos(hno)%r                        = 0.0
       liste_halos(hno)%spin                     = 0.0
       liste_halos(hno)%ek                       = 0.0
       liste_halos(hno)%ep                       = 0.0
       liste_halos(hno)%et                       = 0.0
       liste_halos(hno)%macc                     = 0.0d0
#ifndef SIMPL
       ! structure tree type
       liste_halos(hno)%level                    = 0
       liste_halos(hno)%hosthalo                 = 0
       liste_halos(hno)%hostsub                  = 0
       liste_halos(hno)%nbsub                    = 0
       liste_halos(hno)%nextsub                  = -1
#endif
#ifndef BIG_RUN
       liste_halos(hno)%ncont                    = 0
#endif
    end do
         
    return
    
  end subroutine init_liste_halos
    
  !***********************************************************************
  subroutine new_step

    implicit none

    integer(kind=4)        :: i,ierr
#ifndef BIG_RUN
    integer(kind=4)        :: good_halos,good_subs
#endif
    real(kind=4),parameter :: rel_prec = 1e-5 
    real(kind=4)           :: read_time_ini,read_time_end
    real(kind=8)           :: macc

    write(errunit,*) 
    write(errunit,*) '> time step number ---> ',st_do 

    call cpu_time(read_time_ini)
#ifdef SIMPL    
    call init_liste_halos
#endif    
    ! read a brick file created with part_to_halo.F
    call read_tree_brick(ierr)
    if(ierr.ne.0) then
       ! if file doesn't exit move to the next step
       return
    end if
#ifndef BIG_RUN
    ! make nlr the number of low-resolution particles 
    if (numero_step == 1) then
       massp = minval(mass)
       nlr   = 0
       do i = 1,nbodies
          if (mass(i) > massp*(1.+rel_prec)) then 
             nlr = nlr + 1
          end if
       end do
       write(errunit,*) '> number of hi-resolution particles  : ',nbodies-nlr
    end if
#endif

    ! initialize massToAccrete of each particle at first timestep.
    if (numero_step == 1) then 
#ifdef BIG_RUN 
       MassToAccrete(:) = real(massp,8)
#else
       do i = 1,nbodies
          MassToAccrete(i) = real(mass(i),8)
       end do
#endif
    end if

    allocate(first_part(0:(nb_of_halos(numero_step)+nb_of_subhalos(numero_step))),stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) '> not enough memory for first_part allocation'
       stop
    endif

    call make_linked_list
    write(errunit,*) '> number of halos         : ',nb_of_halos(numero_step)
    write(errunit,*) '> number of subhalos      : ',nb_of_subhalos(numero_step)
    
    ! compute mass of newly accreted particles onto each halo
    do i = 1, nb_of_halos(numero_step) + nb_of_subhalos(numero_step)
       macc = 0.d0
       ierr = first_part(i)
       do while (ierr /= -1) 
          macc = macc + MassToAccrete(ierr)
          MassToAccrete(ierr) = 0.0d0
          ierr  = linked_list(ierr)
       end do
       liste_halos(i)%macc = macc
    end do
    ! output sum (over haloes) of accreted mass during current timestep.
    write(444,'(i3,1x,e14.6)') numero_step, sum(MassToAccrete)

!!$#ifndef BIG_RUN    
!!$    good_halos = 0
!!$    good_subs  = 0 
!!$    do i = 1,nb_of_halos(numero_step)+ nb_of_subhalos(numero_step)
!!$       liste_halos(i)%ncont = 0
!!$       ierr                 = first_part(i)
!!$       do while (ierr /= -1)
!!$          if (mass(ierr) > massp*(1.+rel_prec)) liste_halos(i)%ncont = liste_halos(i)%ncont +1
!!$          ierr = linked_list(ierr)
!!$       end do
!!$       if (liste_halos(i)%ncont == 0) then 
!!$          if (liste_halos(i)%level == 1) then 
!!$             good_halos = good_halos + 1
!!$          else
!!$             good_subs  = good_subs + 1
!!$          endif
!!$       endif
!!$    end do
!!$    write(errunit,*) '> number of uncontaminated halos, subhalos : ',good_halos,good_subs
!!$#endif 

    if (numero_step == 1) then
       do i = 1,nb_of_halos(numero_step)+ nb_of_subhalos(numero_step)
          ! initialize first halos so that they have one father --> the background
          tree(numero_step,i)%my_fathers%nb_fathers = 1
          allocate(tree(numero_step,i)%my_fathers%list_fathers(1),stat=ierr)
          if (ierr /= 0) then
             write(errunit,*) '> not enough memory to allocate list_fathers in new_step'
             stop
          endif
          ! father number is 0 --> background
          tree(numero_step,i)%my_fathers%list_fathers(1) = 0
          allocate(tree(numero_step,i)%my_fathers%mass_fathers(1),stat=ierr)
          if (ierr /= 0) then
             write(errunit,*) '> not enough memory to allocate mass_fathers in new_step'
             stop
          endif
#ifndef SIMPL
          ! percentage of mass given to the halo by its father is 100.
          tree(numero_step,i)%my_fathers%mass_fathers(1) = 100d0
#endif
#ifdef SIMPL
          ! percentage of mass given to the halo by its father is 100.
          tree(numero_step,i)%my_fathers%mass_fathers(1) = real(liste_halos(i)%m,8)
#endif
       end do
    else
       call det_ts_fathers
       call det_ts_sons
    end if
    deallocate(first_part)

    if (numero_step == nsteps) then
       deallocate(liste_parts,ex_liste_parts,linked_list)
    else
       ! copy list of halo appartenance of particles into another variable to 
       ! be able to keep track of their whereabouts at the next time step
       ex_liste_parts = liste_parts
       ! re-initialize liste of parts :
       liste_parts = 0
    end if

    call write_ts_data
#ifndef SIMPL
    call update_tree
#endif
    call cpu_time(read_time_end)

    write(errunit,*) '> time step took          : ',nint(read_time_end-read_time_ini),' seconds'

    return

  end subroutine new_step

!*************************************************************************
  subroutine make_linked_list

  ! creates a linked list of particles (identical as make_linked_list in DM_FOF.F)

    implicit none

    integer(kind=4)             :: i,index1,index2,ierr
    integer(kind=4),allocatable :: current_ptr(:)

    allocate(current_ptr(0:(nb_of_halos(numero_step)+nb_of_subhalos(numero_step))),stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) '> not enough memory for current_ptr allocation'
    endif

    first_part  = -1
    current_ptr = -1
    linked_list = -1
    
    do i = 1,nbodies
       index1 = liste_parts(i)
       if (first_part(index1) ==  -1) then
          first_part(index1)  = i
          current_ptr(index1) = i
       else
          index2              = current_ptr(index1)
          linked_list(index2) = i
          current_ptr(index1) = i
       endif
    end do

    do i = 0,nb_of_halos(numero_step)
       if (current_ptr(i) == -1) cycle
       index2              = current_ptr(i)
       linked_list(index2) = -1
    end do

    deallocate(current_ptr)

    return

  end subroutine make_linked_list


!***********************************************************************
  subroutine det_ts_fathers

  ! determines the progenitors of all halos at a given time step, and the 
  ! mass they get from each progenitor

    implicit none

    integer(kind=4)    :: i,indexp
    real(kind=4)       :: mtot_fath
    type(halo),pointer :: h
    type(tree_struct), pointer :: htr
    integer(kind=4)    :: j
    
    do i = 1,nb_of_halos(numero_step)  + nb_of_subhalos(numero_step)
       htr    => tree(numero_step,i)
       h      => liste_halos(i)
       indexp = first_part(i)
       do while (indexp /= -1)
          ! determine if particles of halo liste_halos(numero_step,i) belong (at the previous timestep) to a halo h_prog which 
          ! already is in the list of progenitors of liste_halos(numero_step,i) or not:
          !    -- if it is in this list, we increase the mass contribution of h_prog to the mass of liste_halos(numero_step,i)
          !    -- if it is not, we add h_prog to this list  
          call new_father(htr,indexp)
          ! next particle
          indexp = linked_list(indexp)
       end do
       if(htr%my_fathers%nb_fathers.gt.0) then
          mtot_fath = 0.
          
          do j = 1,htr%my_fathers%nb_fathers
             mtot_fath = mtot_fath + real(htr%my_fathers%mass_fathers(j),4)
             if(htr%my_fathers%mass_fathers(j).gt.h%m*(1.+1e-5)) then
                write(errunit,'(1x,a,i2,2(1x,E13.5))') '> ERROR,j, mass_father(j),h%m:',j, htr%my_fathers%mass_fathers(j)*1e11,h%m*1e11
             end if
          end do
          if(mtot_fath.gt.h%m*(1.+1e-5)) then
             write(errunit,*) '> Error for halo:',numero_step,i
             write(errunit,'(1x,a,2(1x,E14.6))') '> sum(mass_fathers),h%m:',mtot_fath*1e11,h%m*1e11
             !stop
          end if
          
       end if
#ifndef SIMPL 
       ! compute the percentage of the total mass of the halo contributed by each 
       ! one of its progenitors (background is included)   
       do j = 1,htr%my_fathers%nb_fathers
          htr%my_fathers%mass_fathers(j) = htr%my_fathers%mass_fathers(j) / real(h%m,8) * 100.0
       end do
#endif

    end do

    return

  end subroutine det_ts_fathers



  !******************************************************************************************
  subroutine new_father(htr,indexp)

  ! determine where particle indexp of halo h belongs (at the previous timestep) 
  ! and update progenitor list accordingly 

    implicit none

    integer(kind=4)             :: indexp,i
    type(tree_struct), pointer  :: htr
    integer(kind=4)             :: numero_ancien_halo
    integer(kind=4)             :: ierr,new_fathers_nb
    integer(kind=4),allocatable :: tabtemp(:)
    real(kind=8),allocatable    :: tabtempr(:)

    ! number of halo to which particle indice belonged at previous time step
    numero_ancien_halo = ex_liste_parts(indexp)

    ! NB: the first time we go in this loop htr%my_fathers%nb_fathers = 0
    do i = 1,htr%my_fathers%nb_fathers
       ! if progenitor is already in the list increase its contribution and exit the routine
       if (htr%my_fathers%list_fathers(i) == numero_ancien_halo) then
#ifndef BIG_RUN
          htr%my_fathers%mass_fathers(i) = htr%my_fathers%mass_fathers(i) + real(mass(indexp),8)
#endif
#ifdef BIG_RUN
          htr%my_fathers%mass_fathers(i) = htr%my_fathers%mass_fathers(i) + real(massp,8) 
#endif
          return
       endif
    end do

    ! if progenitor has not been found then add a new one to the list 
    new_fathers_nb             = htr%my_fathers%nb_fathers + 1
    htr%my_fathers%nb_fathers  = new_fathers_nb

    if (associated(htr%my_fathers%list_fathers)) then

       ! temporary copy of progenitor list
       allocate(tabtemp(new_fathers_nb),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) '> not enough memory for tabtemp allocation in new_fathers 1'
          stop
       endif
       tabtemp(1:new_fathers_nb - 1) = htr%my_fathers%list_fathers(1:new_fathers_nb - 1)
       tabtemp(new_fathers_nb)       = numero_ancien_halo
       deallocate(htr%my_fathers%list_fathers)

       ! create new list of dimension = old_dimension + 1
       allocate(htr%my_fathers%list_fathers(new_fathers_nb),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) '> not enough memory for list_fathers allocation in new_father 1'
          stop
       endif
       htr%my_fathers%list_fathers(1:new_fathers_nb) = tabtemp(1:new_fathers_nb)
       deallocate(tabtemp)

       ! idem with other properties
       allocate(tabtempr(new_fathers_nb),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) '> not enough memory for tabtempr allocation in new_fathers 1'
          stop
       endif
       tabtempr(1:new_fathers_nb - 1) = htr%my_fathers%mass_fathers(1:new_fathers_nb - 1)
#ifndef BIG_RUN
       tabtempr(new_fathers_nb)     = mass(indexp)
#endif
#ifdef BIG_RUN
       tabtempr(new_fathers_nb)     = massp
#endif
       deallocate(htr%my_fathers%mass_fathers)

       allocate(htr%my_fathers%mass_fathers(new_fathers_nb),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) '> not enough memory for mass_fathers allocation in new_father 1'
          stop
       endif
       htr%my_fathers%mass_fathers(1:new_fathers_nb) = tabtempr(1:new_fathers_nb)
       deallocate(tabtempr)

    else

       allocate(htr%my_fathers%list_fathers(new_fathers_nb),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) '> not enough memory for list_fathers allocation in new_father 2'
          stop
       endif
       htr%my_fathers%list_fathers(new_fathers_nb) = numero_ancien_halo

       allocate(htr%my_fathers%mass_fathers(new_fathers_nb),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) '> not enough memory for mass_fathers allocation in new_father 2'
          stop
       endif
#ifndef BIG_RUN
       htr%my_fathers%mass_fathers(new_fathers_nb) = real(mass(indexp),8)
#endif
#ifdef BIG_RUN
       htr%my_fathers%mass_fathers(new_fathers_nb) = real(massp,8)
#endif

    endif

    return

  end subroutine new_father

  !***********************************************************************
  subroutine det_ts_sons
  
  !  determine halo sons 

    implicit none

    integer(kind=4)             :: st,j,k,l,next_st,indice
    logical(kind=4)             :: found_it
    integer(kind=4)             :: ierr
    integer(kind=4),allocatable :: tabtemp(:)
    real(kind=8)                :: my_massfrac,mainson_massfrac
    integer(kind=4)             :: idmainson 

    st = numero_step -1 
    next_st = st + 1

    
    do j = 1,nb_of_halos(next_st)+nb_of_subhalos(next_st)
       
       do k = 1,tree(next_st,j)%my_fathers%nb_fathers
          
          ! progenitor of son is background so discard it 
          if (tree(next_st,j)%my_fathers%list_fathers(k) == 0) then
             cycle
          endif
          
          indice   = tree(next_st,j)%my_fathers%list_fathers(k)
          ! determine if son is already in list of sons
          found_it = .false.
          ! NB: if we arrive here for the first time liste_halos(st,indice)%my_sons%nb_sons = 0
          do l = 1,tree(st,indice)%my_sons%nb_sons
             if (tree(st,indice)%my_sons%list_sons(l) == j) then
                found_it = .true.
                exit
             endif
          end do
          
          ! if son not already in list of sons then add it to list of sons
          if (.not. found_it) then
             
             if (associated(tree(st,indice)%my_sons%list_sons)) then              
                
                ! temporary copy
                allocate(tabtemp(tree(st,indice)%my_sons%nb_sons+1),stat=ierr)
                if (ierr /= 0) then
                   write(errunit,*) '> not enough memory for tabtemp allocation'
                   stop
                endif
                tabtemp(1:tree(st,indice)%my_sons%nb_sons)   = &
                     & tree(st,indice)%my_sons%list_sons(1:tree(st,indice)%my_sons%nb_sons)
                tabtemp(tree(st,indice)%my_sons%nb_sons+1) = j
                deallocate(tree(st,indice)%my_sons%list_sons)
                
                ! new list of sons of dimension = old_dimension + 1 
                allocate(tree(st,indice)%my_sons%list_sons(tree(st,indice)%my_sons%nb_sons+1),stat=ierr)
                if (ierr /= 0) then
                   write(errunit,*) '> not enough memory for list_sons allocation in det_sons 1'
                   stop
                endif
                tree(st,indice)%my_sons%list_sons(1:tree(st,indice)%my_sons%nb_sons+1) = &
                     & tabtemp(1:tree(st,indice)%my_sons%nb_sons+1)
                deallocate(tabtemp)
                
                tree(st,indice)%my_sons%nb_sons = tree(st,indice)%my_sons%nb_sons + 1           
                my_massfrac = &
                    &tree(next_st,j)%my_fathers%mass_fathers(k)*liste_halos(j)%m

                idmainson = tree(st,indice)%my_sons%main_son
                do l = 1, tree(next_st,idmainson)%my_fathers%nb_fathers
                    if (tree(next_st,idmainson)%my_fathers%list_fathers(l) == indice) then
                        mainson_massfrac = &
                            &tree(next_st,idmainson)%my_fathers%mass_fathers(l)*liste_halos(idmainson)%m
                    end if
                end do

                if (my_massfrac > mainson_massfrac) then
                    tree(st,indice)%my_sons%main_son = j
                end if
                
             else
                
                allocate(tree(st,indice)%my_sons%list_sons(1),stat=ierr)
                if (ierr /= 0) then
                   write(errunit,*) '> not enough memory for list_sons allocation in det_sons 2'
                   stop
                endif
                tree(st,indice)%my_sons%list_sons(1) = j
                tree(st,indice)%my_sons%nb_sons      = 1
                tree(st,indice)%my_sons%main_son      = j
                
             endif
             
          endif
          
       end do
       
    end do

    return

  end subroutine det_ts_sons

#ifndef SIMPL
!*************************************************************************
  subroutine update_tree
    
    implicit none
    integer(kind=4) :: hno
    
    do hno = 1, nb_of_halos(numero_step) + nb_of_subhalos(numero_step)
       tree(numero_step,hno)%host     = liste_halos(hno)%hostsub  !! replace hosthalo with hostsub and name it host
       tree(numero_step,hno)%nextsub  = liste_halos(hno)%nextsub
    end do
    
    call init_liste_halos

    return

  end subroutine update_tree
#endif

!*************************************************************************
  subroutine finish_all
    
    implicit none
    
    deallocate(liste_halos,nb_of_halos,nb_of_subhalos,aexp,age_univ,omega_t)
    deallocate(tree)
    
  end subroutine finish_all

!*************************************************************************

  subroutine define_bushes

    ! JB -
    ! assign a BushID to all halos and subhalos
    ! -> a bushID is a unique ID shared by all haloes connected in anyway at anytime (descendent/progenitor and/or halo/sub-halo). 

    implicit none

    integer(kind=4) :: ih,ts
    integer(kind=4) :: BushID
    real(kind=4)    :: t_start,t_end
#ifndef SIMPL
    integer(kind=4) :: st
    integer(kind=4) :: ip,ifa,iso
#endif
    ! initialise BushID
    BushID = 0
    
    call CPU_time(t_start)
    write(errunit,*)
    write(errunit,*) '> computing BushIDs'
#ifndef SIMPL
    ! loop on all haloes, starting at z=0
    do st = nsteps,1,-1
       do ih = 1, nb_of_halos(st)
          ! only consider halos which are not already part of a bush
          if (tree(st,ih)%BushID == -1) then  
             BushID                 = BushID + 1
             ts                     = st
             call get_my_tree(ts,ih,BushID)
          end if
       end do
    end do
#else
    ! only consider halo that reach the last timestep
    do ih = 1, nb_of_halos(nsteps)
       if (tree(nsteps,ih)%BushID == -1) then  
          BushID                 = BushID + 1
          ts                     = nsteps
          call get_my_tree(ts,ih,BushID)
       end if
    end do
#endif
    
    write(errunit,*) '> check BushIDs'

    ! check that all halos (and subhalos) have been assigned a correct BushID
    do ts = 1, nsteps 
       do ih = 1, nb_of_halos(ts) + nb_of_subhalos(ts)
#ifdef SIMPL
          call check_bushes(ts,ih)
#endif
          
#ifndef SIMPL
          ! check that my dads have my BushID
          if (ts > 1) then 
             
             do ip = 1,tree(ts,ih)%my_fathers%nb_fathers
                ifa = tree(ts,ih)%my_fathers%list_fathers(ip)
                if (ifa > 0) then ! not background ... 
                   if (tree(ts-1,ifa)%BushID .ne. tree(ts,ih)%BushID) then 
                      print*,'bye bye father ... ',tree(ts-1,ifa)%BushID,tree(ts,ih)%BushID
                      stop
                   end if
                end if
             end do

          end if
          
          ! check that my sons have my BushID
          if (ts < nsteps) then 

             do ip = 1,tree(ts,ih)%my_sons%nb_sons
                iso = tree(ts,ih)%my_sons%list_sons(ip)
                if (iso > 0) then ! not background ...
                   if (tree(ts+1,iso)%BushID .ne. tree(ts,ih)%BushID) then 
                      print*,'bye bye sons ... ',tree(ts+1,iso)%BushID,tree(ts,ih)%BushID
                      stop
                   end if
                end if
             end do

          end if

          ! check myself
          if (tree(ts,ih)%BushID <= 0) then 
             write(errunit,*) 'no :',ih,tree(ts,ih)%host,tree(ts,ih)%nextsub,nb_of_halos(ts),nb_of_subhalos(ts)
             stop
          end if
#endif

       end do
    end do
    call CPU_time(t_end)
    write(errunit,*) '> operation took          :',int(t_end-t_start),' seconds'
    
    return

  end subroutine define_bushes

#ifndef SIMPL 
!*************************************************************************

  recursive subroutine get_my_tree(ts,ih,BushID)
    
    ! JB -

    implicit none
    
    integer(kind=4) :: ts,ih,ip,ifa,iso,sub
    integer(kind=4) :: BushID
    
    if (tree(ts,ih)%BushID == -1) then 

       tree(ts,ih)%BushID = BushID

       ! go down
       do ip = 1,tree(ts,ih)%my_fathers%nb_fathers
          ifa = tree(ts,ih)%my_fathers%list_fathers(ip)
          if (ifa > 0) then ! not background ... 
             if (tree(ts-1,ifa)%BushID == -1) then 
                call get_my_tree(ts-1,ifa,BushID)
             end if
          end if
       end do

       ! go up 
       do ip = 1,tree(ts,ih)%my_sons%nb_sons
          iso = tree(ts,ih)%my_sons%list_sons(ip)
          if (iso > 0) then ! not background ...
             if (tree(ts+1,iso)%BushID == -1) then 
                call get_my_tree(ts+1,iso,BushID)
             end if
          end if
       end do

       ! go sideways
       ! NB: host == hostsub from HaloMaker -> 
       ! - if (and only if) host = 0 then the thing is a main halo.
       ! - if (and only if) host > 0 then the thing is a sub with ilevel > 1. 
       !   In that case, the thing with index tree(ts,ih)%host contains the thing ih.
       ! NB: nextsub spans all levels of substructure.
       if (tree(ts,ih)%host > 0) then ! this is a sub. 
          sub = tree(ts,ih)%host
       else                           ! this is a main halo : start loop on its subs
          sub = tree(ts,ih)%nextsub
       end if
       do while (sub > 0) 
          if (tree(ts,sub)%BushID == -1) then 
             call get_my_tree(ts,sub,BushID)
          end if
          sub = tree(ts,sub)%nextsub
       end do

    end if
       
    return
    
  end subroutine get_my_tree

#endif

#ifdef SIMPL       
!*************************************************************************
  recursive subroutine get_my_tree(ts,ih,BushID)
    
    ! JB - but DT was here

    implicit none
    
    integer(kind=4) :: ts,ih,ip,ifa,iso,sub
    integer(kind=4) :: BushID
    
    if (tree(ts,ih)%BushID.gt.0) return

    if (tree(ts,ih)%level.le.0) then
       write(errunit,*) '> Error in get_my_tree for halo:',ts,ih
       write(errunit,*) '> tree(ts,ih)%level,tree(ts,ih)%hostsub:',tree(ts,ih)%level,tree(ts,ih)%hostsub
       stop
    end if
    
    tree(ts,ih)%BushID = BushID

    ! go down
    do ip = 1,tree(ts,ih)%my_dads%nb_dads
       ifa = tree(ts,ih)%my_dads%list_dads(ip)
       if (ifa > 0) then ! not background ... 
          if (tree(ts-1,ifa)%BushID == -1) then 
             call get_my_tree(ts-1,ifa,BushID)
          end if
       end if
    end do
    
    iso = tree(ts,ih)%my_sons%main_son
    ! go up 
    if (iso > 0) then ! not background ...
       if (tree(ts+1,iso)%BushID == -1) then 
          call get_my_tree(ts+1,iso,BushID)
       end if
    end if
    
    ! go sideways
    ! - if level == 1 this is the main halo
    ! - if level >  1 this is a subhalo
    ! - if level < it was a subhalo and was removed
    if (tree(ts,ih)%level .gt. 1) then ! this is a sub. 
       sub = tree(ts,ih)%hosthalo
    else                           ! this is a main halo : start loop on its subs
       sub = tree(ts,ih)%nextsub
    end if
    ! starting to look sideway from the main halo
    ! NB: nextsub spans all levels of substructure.
    do while (sub > 0) 
       if (tree(ts,sub)%BushID == -1) then 
          call get_my_tree(ts,sub,BushID)
       end if
       sub = tree(ts,sub)%nextsub
    end do
    
    return
    
  end subroutine get_my_tree

  !*************************************************************************
  subroutine check_bushes(st,ih)
    
    implicit none
    
    integer(kind=4) :: ih,st
    integer(kind=4) :: BushID
    integer(kind=4) :: ip,ifa,iso,iho,isb
    
    if(tree(st,ih)%BushID.le.0) return
    
    BushID = tree(st,ih)%BushID
    
    if(BushID.le.0) then
       write(errunit,*) '> Error for halo:',st,ih,tree(st,ih)%level
       write(errunit,*) '> BushID:',BushID
       stop
    end if
    
    if(st.gt.1) then
       ! check that my dads have my BushID
       do ip = 1,tree(st,ih)%my_dads%nb_dads
          ifa = tree(st,ih)%my_dads%list_dads(ip)
          if(ifa.le.0) then
             write(errunit,*) '> Error for halo:',st,ih,tree(st,ih)%level
             write(errunit,*) '> dad id:',ifa
             stop
          end if
          if(tree(st-1,ifa)%BushID.ne.BushID) then
             write(errunit,*) '> Error for halo:',st,ih,tree(st,ih)%level
             write(errunit,*) '> dad:',ifa,tree(st-1,ifa)%level
             write(errunit,*) '> halo BushID, dad BushID:',BushID,tree(st-1,ifa)%BushID
             stop
          end if
       end do
    end if
    
    if(st.lt.nsteps) then
       ! check that main_son have my BushID
       iso = tree(st,ih)%my_sons%main_son
       if (iso .gt. 0) then ! not background ...
          if(tree(st+1,iso)%BushID.ne.BushID) then
             write(errunit,*) '> Error for halo:',st,ih,tree(st,ih)%level
             write(errunit,*) '> main son:',iso,tree(st+1,iso)%level
             write(errunit,*) '> halo BushID, main son BushID:',BushID,tree(st+1,iso)%BushID
             stop
          end if
       end if
    end if
    
    iho = tree(st,ih)%hostsub
    if(iho.gt.0) then ! not background ...
       if(tree(st,iho)%BushID.ne.BushID) then
          write(errunit,*) '> Error for halo:',st,ih,tree(st,ih)%level
          write(errunit,*) '> host:', iho,tree(st,iho)%level
          write(errunit,*) '> halo BushID, host BushID:',BushID,tree(st,iho)%BushID
          stop
       end if
    end if
    
    isb = tree(st,ih)%nextsub
    if(isb.gt.0) then ! not background ...
       if(tree(st,isb)%BushID.ne.BushID) then
          write(errunit,*) '> Error for halo:',st,ih,tree(st,ih)%level
          write(errunit,*) '> next sub:',isb,tree(st,isb)%level
          write(errunit,*) '> halo BushID, next sub BushID:',BushID,tree(st,isb)%BushID
          stop
       end if
    end if
  
    return
  
  end subroutine check_bushes

#endif
  
!*************************************************************************
!/////////////////////////////////////////////////////////////////////////

end module compute_tree





