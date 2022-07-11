module simpl_merger_tree

  use tree_defs
  
  !***********************************************************************
  ! debugging parameter
  !***********************************************************************
  ! set value above one to print out info about 
  ! stbug... is a step number and ihbug... is a halo number 
  ! how halo stbugson, ihbugson main_son is computed
  ! how halo stbugdad, ihbugdad dads are computed
  integer(kind=4),parameter :: stbugson = -1
  integer(kind=4),parameter :: ihbugson = -1 
  integer(kind=4),parameter :: stbugdad = -1 
  integer(kind=4),parameter :: ihbugdad = -1 
  !***********************************************************************
  
contains

  !***********************************************************************
  subroutine simplify_merger_tree
    
    implicit none
    integer(kind=4) :: st,old_st,nxt_st,i
    real(kind=4)    :: t_start,t_end

    call CPU_time(t_start)
    write(errunit,*) 
    write(errunit,*) '> tree simplification in progress'
    if(stbugson.gt.0.and.ihbugson.gt.0) then
       write(errunit,*) '> printing debugging son information'
       write(errunit,*) '> for halo:',stbugson,ihbugson
       write(errunit,*)
    end if
    if(stbugdad.gt.0.and.ihbugdad.gt.0) then
       write(errunit,*) '> printing debugging dad information'
       write(errunit,*) '> for halo:',stbugdad,ihbugdad
       write(errunit,*)
    end if

    do st = 1, nsteps
       old_st = st - 1
       nxt_st = st + 1
       
       do i = 1, nb_of_halos(st) + nb_of_subhalos(st)
          ! compute_main son
          
          if(st.eq.stbugson.and.i.eq.ihbugson) then
             write(errunit,*) '> calling compute_son for:',st,i
          end if
          if(nxt_st.le.nsteps) call compute_son(nxt_st,st,i)
          ! compute mainprog, and dads list
          if(st.eq.stbugdad.and.i.eq.ihbugdad) then
             write(errunit,*) '> calling compute_dads for:',st,i
          end if
          if(old_st.ge.1) call compute_dads(st,old_st,i)
       end do
       
#ifdef S_DAD_FOR_SUB
       ! correction of the merger tree for substructures
       if(nb_of_subhalos(st).gt.0.and.old_st.ge.1) then
          do i = nb_of_halos(st) + 1, nb_of_halos(st) + nb_of_subhalos(st) 
             if(tree(st,i)%my_dads%nb_dads.le.0) then
                tree(st,i)%frag = 1
                if(st.eq.stbugdad.and.i.eq.ihbugdad) then
                   write(errunit,*)
                   write(errunit,*) '> no progenitor for:',st,i,tree(st,i)%level
                   write(errunit,*) '> searching for new progenitors', st,old_st,i
                end if
                call search_dad_for_subhalo(st,old_st,i)
             end if
          end do
       end if
#endif
    end do
    
#ifdef R_DADLESS_SUB
    if(maxval(nb_of_subhalos).gt.0) then
       write(errunit,*)
       write(errunit,*) '> removing subhalos without progenitors'
       do st = 1, nsteps
          call remove_subs_without_dad(st)
       end do
    end if
#endif

#ifdef COL_TREE
    if(maxval(nb_of_subhalos).gt.0) then
       write(errunit,*)
       write(errunit,*) '> removing all subhalos and keep fly-bys'
       do st =nsteps,1,-1
          call collapse_tree(st)
       end do
    end if
#endif
  
    write(errunit,*)
    write(errunit,*) '> detect fragments and deallocate old merger tree'
    do st = 1,nsteps
       do i = 1,nb_of_halos(st) + nb_of_subhalos(st)
          call detect_fragment(st,i)
          call deallocate_old_tree(st,i)
       end do
    end do
    call CPU_time(t_end)
    write(errunit,*) '> tree simplification took:',int(t_end-t_start),' seconds'
    
    return
    
    return
    
  end subroutine simplify_merger_tree

  !***********************************************************************
  subroutine compute_son(nxt_st,st,i)

    implicit none
    integer(kind=4)      :: nxt_st,st,i,ifath,idfath,ison,idson
    real(kind=4)         :: mass,mfath
    type(son),pointer    :: sons
    type(father),pointer :: fathson


    if(st.eq.stbugson.and.i.eq.ihbugson) then
       write(errunit,*) '> in compute_son:',nxt_st,st,i
    end if
    
    sons  => tree(st,i)%my_sons

    if(sons%nb_sons.le.0) then
       sons%main_son = -1
       if(st.eq.stbugson.and.i.eq.ihbugson) then
          write(errunit,*) '> no sons, no main son'
       end if
       return
    else
       if(sons%main_son.gt.0) then
          if(tree(nxt_st,sons%main_son)%frag.gt.0.and.tree(nxt_st,sons%main_son)%level.gt.0) then
             if(st.eq.stbugson.and.i.eq.ihbugson) then
                write(errunit,*) '> main_son,frag,level'
                write(errunit,*) '>',sons%main_son,tree(nxt_st,sons%main_son)%frag,tree(nxt_st,sons%main_son)%level
                write(errunit,*) '> main son was corrected with search_dads_for_sub'
             end if
             return
          end if
       end if
    end if
   
    
    !reinitialise main_son
    sons%main_son = -1
 
    mass  = 0.
    if(st.eq.stbugson.and.i.eq.ihbugson) then
       write(errunit,*) '> sons%nb_sons:',sons%nb_sons
       write(errunit,*) '> sons,ifath,mfath'
    end if
    do ison = 1, sons%nb_sons
       idson = sons%list_sons(ison)
       if(idson.le.0) then
          write(errunit,*) '> Error in compute_son for halo:',st,i,tree(st,i)%level
          write(errunit,*) '> idson is nil'
          stop
       end if
       fathson => tree(nxt_st,idson)%my_fathers
       ifath   = 0
       idfath  = -1
       do while(idfath.ne.i.and.ifath.lt.fathson%nb_fathers) 
          ifath  = ifath + 1
          idfath = fathson%list_fathers(ifath)
       end do
       if(idfath.ne.i) then
          write(errunit,*) '> Error in compute_son for halo:',st,i,tree(st,i)%level
          write(errunit,*) '> halo is not among its sons'' fathers'
          write(errunit,*) '> for its son:',nxt_st,idson
          write(errunit,*) '> fathson%list_fathers'
          do ifath = 1, fathson%nb_fathers
             write(errunit,*) '>',fathson%list_fathers(ifath)
          end do
          stop
       end if
       mfath     = real(fathson%mass_fathers(ifath),4)
       if(mfath.gt.mass) then
          mass          = mfath
          sons%main_son = idson
       end if
       if(st.eq.stbugson.and.i.eq.ihbugson) &
            write(errunit,*) '>',idson,ifath,mfath
       
    end do
    if(st.eq.stbugson.and.i.eq.ihbugson) &
         write(errunit,*) '> main_son:',sons%main_son

    return

  end subroutine compute_son
  
  !***********************************************************************
  subroutine compute_dads(st,old_st,i)
    
    implicit none
    integer(kind=4)      :: st,old_st,i
    integer(kind=4)      :: ndads,ifath,idfath
    integer(kind=4)      :: ison,idson
    type(father),pointer :: fathers
    type(dad),pointer    :: dads
    

    dads => tree(st,i)%my_dads

    ! reinitialise
    if(associated(dads%list_dads)) deallocate(dads%list_dads)
    if(associated(dads%mass_dads)) deallocate(dads%mass_dads)
    dads%nb_dads = 0

    if(st.eq.stbugdad.and.i.eq.ihbugdad) then
       write(errunit,*) '> in compute_dads for:',st,old_st,i
    end if

    if(tree(st,i)%my_fathers%nb_fathers.le.0) then
       if(st.eq.stbugdad.and.i.eq.ihbugdad) then
          write(errunit,*) '> no fathers: no dads'
       end if
       return
    else
       fathers => tree(st,i)%my_fathers
       call ord_fathers(fathers)
    end if
    
    ! counting dads
    if(st.eq.stbugdad.and.i.eq.ihbugdad) then
       write(errunit,*) '> fath,mfath,son'
       write(errunit,*) '> fathers%nb_fathers:',fathers%nb_fathers
    end if
    ndads       = 0
    do ifath = 1, fathers%nb_fathers
       idfath = fathers%list_fathers(ifath)
       if(idfath.gt.0) then
          if(st.eq.stbugdad.and.i.eq.ihbugdad) &
               write(errunit,*) '>',old_st,fathers%list_fathers(ifath),real(fathers%mass_fathers(ifath),4),&
               tree(old_st,idfath)%my_sons%main_son
          if(tree(old_st,idfath)%my_sons%main_son.eq.i) then
             ndads = ndads + 1
          end if
          ! lets check that i is among idfath sons
          ison  = 0
          idson = -1
          do while(idson.ne.i.and.ison.lt.tree(old_st,idfath)%my_sons%nb_sons)
             ison  = ison + 1
             idson = tree(old_st,idfath)%my_sons%list_sons(ison)
          end do
          if(idson.ne.i) then
             write(errunit,*) '> Error in compute_dads for ',st,i,tree(st,i)%level
             write(errunit,*) '> halo is not among its fathers'' sons'
             write(errunit,*) '> for its father:',old_st,idfath
             write(errunit,*) '> tree(old_st,idfath)%list_sons'
             do ison = 1,tree(old_st,idfath)%my_sons%nb_sons
                write(errunit,*) '>',tree(old_st,idfath)%my_sons%list_sons(ison)
             end do
             stop
          end if
       else if(st.eq.stbugdad.and.i.eq.ihbugdad) then
          write(errunit,*) '>',fathers%list_fathers(ifath),real(fathers%mass_fathers(ifath),4),0
       end if
    end do
    dads%nb_dads = ndads
    if(st.eq.stbugdad.and.i.eq.ihbugdad) &
         write(errunit,*) '> dads%nb_dads:',dads%nb_dads

    if(ndads.le.0) return

    ! compute list_dads and mass_dads
    ndads = 0
    allocate(dads%list_dads(dads%nb_dads))
    allocate(dads%mass_dads(dads%nb_dads))
    if(st.eq.stbugdad.and.i.eq.ihbugdad) &
         write(errunit,*) '> dads%list_dads,dads%mass_dads'
    do ifath = 1, fathers%nb_fathers
       idfath = fathers%list_fathers(ifath)
       if(idfath.gt.0) then
          if(tree(old_st,idfath)%my_sons%main_son.eq.i) then
             ndads = ndads + 1
             if(ndads.gt.dads%nb_dads) then
                write(errunit,*) '> Error in compute_dads for halo:',st,i
                write(errunit,*) '> dads%nb_dads,ndads:',dads%nb_dads,ndads
                write(errunit,*) '> Wrong number of dads'
                stop 
             end if
             dads%list_dads(ndads) = idfath
             dads%mass_dads(ndads) = real(fathers%mass_fathers(ifath),4)
             if(st.eq.stbugdad.and.i.eq.ihbugdad) &
                  write(errunit,*) '>', dads%list_dads(ndads),dads%mass_dads(ndads)
             ! check that mass_dads is ordered properly
             if(ndads.gt.1) then
                if(dads%mass_dads(ndads).gt.dads%mass_dads(ndads-1)) then
                   write(errunit,*) '> Error in compute_dads for halo:',st,i
                   write(errunit,*) '> id_dad,id_son,mass_dad'
                   do ndads = 1,dads%nb_dads
                      write(errunit,'(a,2(1x,i10),1x,F8.2)') ' >',dads%list_dads(ndads),&
                           tree(st-1,dads%list_dads(ndads))%my_sons%main_son,dads%mass_dads(ndads)
                   end do
                   stop
                end if
             end if
          end if
       end if
    end do
    if(ndads.ne.dads%nb_dads) then
       write(errunit,*) '> Error in compute_dads for halo:',st,i
       write(errunit,*) '> dads%nb_dads,ndads:',dads%nb_dads,ndads
       write(errunit,*) '> Wrong number of dads'
       stop 
    end if
    
    return

  end subroutine compute_dads

  !******************************************************************************************
  subroutine ord_fathers(fath)

    implicit none
    type(father)      :: fath
    integer(kind=4)   :: j,jj,ord_j,ierr
    integer(kind=4),allocatable :: tabtemp(:),ord(:)
    real(kind=8),allocatable    :: tabtempr(:)

    ! might as well sort fathers per decreasing mass
    if(fath%nb_fathers .le.1) return
    
    if(.not.associated(fath%list_fathers)) then
       write(errunit,*) '> Error in ord_father'
       write(errunit,*) '> list_fathers is not associated'
       stop
    end if
    if(.not.associated(fath%mass_fathers)) then
       write(errunit,*) '> Error in ord_father'
       write(errunit,*) '> mass_fathers is not associated'
       stop
    end if
    ! temporary copy of progenitor list
    allocate(tabtemp(fath%nb_fathers),tabtempr(fath%nb_fathers),stat=ierr)
    if (ierr /= 0) then
       write(errunit,*) '> not enough memory for tabtemp tabtempr allocation in ord_fathers'
       stop
    endif
    tabtemp(1:fath%nb_fathers)  = fath%list_fathers(1:fath%nb_fathers)
    tabtempr(1:fath%nb_fathers) = fath%mass_fathers(1:fath%nb_fathers)
    if(fath%nb_fathers.gt.2) then
       allocate(ord(fath%nb_fathers),stat=ierr)
       if (ierr /= 0) then
          write(errunit,*) '> not enough memory for ord allocation in ord_fathers'
          stop
       endif
       ord(1:fath%nb_fathers)      = -1            
       ! sort fathers decreasing order with mass
       call rd_sort(fath%nb_fathers,tabtempr,ord)
       if(minval(ord).le.0.or.maxval(ord).gt.fath%nb_fathers) then
          
          write(errunit,*) '> Error in ord_father'
          write(errunit,*) '> fath%nb_fathers:',fath%nb_fathers
          write(errunit,*) '> jj,ord(jj),tabtemp(jj),tabtempr(jj)'
          do jj = 1,fath%nb_fathers
             write(errunit,'(a,2(1x,i3),1x,i10,1x,E14.5)') &
                  ' >', jj,ord(jj),tabtemp(jj),tabtempr(jj)*1e11
          end do
       end if
       do j = 1,fath%nb_fathers
          fath%list_fathers(j) = tabtemp(ord(j))
          fath%mass_fathers(j) = tabtempr(ord(j))
          if(j.gt.1) then
             if(fath%mass_fathers(j-1).lt.fath%mass_fathers(j)) then
                write(errunit,*) '> Error in ord_fathers'
                write(errunit,*) '> fathers have not been ordered properly'
                write(errunit,*) '> jj,ord(jj),tabtemp(ord(jj),tabtempr(ord(jj))'
                do jj = 1,fath%nb_fathers
                   write(errunit,'(a,2(1x,i3),1x,i10,1x,F10.5)') &
                        ' >', jj,ord(jj),tabtemp(ord(jj)),tabtempr(ord(jj))
                end do
                stop
             end if
          end if
       end do
       deallocate(ord)
    else if(fath%nb_fathers.eq.2) then
       if(tabtempr(2).gt.tabtempr(1)) then
          ! two fathers sorted into ascending order change that
          do j = 1,2
             ord_j                = fath%nb_fathers - j + 1
             fath%list_fathers(j) = tabtemp(ord_j)
             fath%mass_fathers(j) = tabtempr(ord_j)
          end do
       end if
    else
       write(errunit,*) '> Error in ord_fathers'
       write(errunit,*) '> fath%nb_fathers:',fath%nb_fathers
       write(errunit,*) '> should have returned'
       stop
    end if
    deallocate(tabtemp,tabtempr)
      
    return

  end subroutine ord_fathers


#ifdef S_DAD_FOR_SUB
  !***********************************************************************
  recursive subroutine search_dad_for_subhalo(st,old_st,i)

    implicit none
    integer(kind=4)      :: st,old_st,i
    integer(kind=4)      :: ifath,idfath,fathson,prog,iprog
    real(kind=4)         :: massson,mprogmin
    type(father),pointer :: fathers

    if(tree(st,i)%my_fathers%nb_fathers.le.0) return
    
    fathers => tree(st,i)%my_fathers
      
    ! search among list_fathers a progenitor that can be the main progenitor of the subhalo 
    mprogmin = real(sum(fathers%mass_fathers),4)/2.*(1.-1.e-4)
    ifath    = 1
    prog     = -1
    iprog    = -1
    do while(ifath.gt.0) 
       idfath = fathers%list_fathers(ifath)
       if(idfath .gt. 0.and.fathers%mass_fathers(ifath).gt.mprogmin) then
          fathson = tree(old_st,idfath)%my_sons%main_son
          massson = real(fathers%mass_fathers(ifath),4)
          if(fathson.le.0) then
             write(errunit,*) '> Error in search_dad_for_subhalo for halo:',old_st,idfath
             write(errunit,*) '> should have a son.'
             stop
          end if
          if(tree(st,fathson)%my_dads%nb_dads.le.0) then
             write(errunit,*) '> Error in search_dad_for_subhalo for halo:',st,fathson
             write(errunit,*) '> should have at least one dad.'
             stop
          end if
          if(tree(st,fathson)%my_dads%list_dads(1).ne.idfath) then
             ! idfath is not the main prog of it son, we can make the change
             iprog = ifath
          else if(massson.gt.tree(st,fathson)%my_dads%mass_dads(1)  &
               .and.tree(st,fathson)%level.ge.tree(st,i)%level) then
             ! idfath is the main prog of its son, but it gave less mass to its son than to this subhalo
             ! this father is du to a precedent change
             iprog = ifath
          end if
       end if
       ifath = ifath + 1
       if(ifath.gt.fathers%nb_fathers) ifath = -1
    end do

    ! if no candidate found return
    if(iprog.le.0) return

    ! candidate found get data
    prog    = fathers%list_fathers(iprog)
    fathson = tree(old_st,prog)%my_sons%main_son
    ! make the change
    if(old_st.eq.stbugson.and.prog.eq.ihbugson) then
       write(errunit,*) '> in search_dad_for_subhalo for',st,i
       write(errunit,*) '> for halo:',old_st,prog
       write(errunit,*) '> changing son from',fathson,'to',i
    end if
    tree(old_st,prog)%my_sons%main_son = i
    if(st.eq.stbugdad.and.i.eq.ihbugdad) then
       write(errunit,*) '> in search_dad_for_subhalo for',st,i
       write(errunit,*) '> calling compute_dads for',st,old_st,i
    end if
    call compute_dads(st,old_st,i)
    if(st.eq.stbugdad.and.fathson.eq.ihbugdad) then
       write(errunit,*) '> in search_dad_for_subhalo for',st,i
       write(errunit,*) '> calling compute_dads for',st,old_st,fathson
    end if
    call compute_dads(st,old_st,fathson)
   
    if(tree(st,fathson)%my_dads%nb_dads.le.0.and.tree(st,fathson)%level.gt.1) then
       ! new modif we change the routine to make it recursive
       ! if we undid a precedent modification of the merger tree 
       ! then subhalo fathson as no longer any progenitor
       if(st.eq.stbugson.and.prog.eq.ihbugson) then
          write(errunit,*) '> in search_dad_for_subhalo for',st,i
          write(errunit,*) '> no longer any progenitor for:',fathson
          write(errunit,*) '> searching for new progenitors'
       end if
       if(st.eq.stbugdad.and.fathson.eq.ihbugdad) then
          write(errunit,*) '> in search_dad_for_subhalo for',st,i
          write(errunit,*) '> no longer any progenitor for:',fathson
          write(errunit,*) '> searching for new progenitors',st,old_st,fathson
       end if
       call search_dad_for_subhalo(st,old_st,fathson)
    end if
  
    return
    
  end subroutine search_dad_for_subhalo

#endif

#ifdef R_DADLESS_SUB
!***********************************************************************
  subroutine remove_subs_without_dad(st)
    
    implicit none
    integer(kind=4) :: st
    integer(kind=4) :: i,host

    if(nb_of_subhalos(st) .le. 0) return
    
    do i = nb_of_halos(st)+1, nb_of_halos(st)+nb_of_subhalos(st)
       ! subhalo is always checked after its host
       if(tree(st,i)%my_dads%nb_dads.le.0) then
          ! remove subhalo without dad
          ! search unremoved host
          host = tree(st,i)%hostsub
          if(host.gt.0) then
             do while(tree(st,host)%level.gt.1.and.tree(st,host)%my_dads%nb_dads.le.0)
                host = tree(st,host)%hostsub
                if(host.le.0) stop '> host shouldn''t be nil'
             end do
          else
             stop '> host is le 0'
          end if
          
          call mod_structure_tree(st,i,host)
          call mod_merger_tree(st,i,host)
          
       end if
    end do

    return

  end subroutine remove_subs_without_dad

#endif

#ifdef COL_TREE
  !***********************************************************************
  subroutine collapse_tree(st)

    implicit none
    integer(kind=4) :: st,i,host,ison
    logical         :: remove,ifhost

    if(nb_of_subhalos(st).le.0) return

    ! removing subhalos from merger trees.
    do i = nb_of_halos(st)+1,nb_of_halos(st)+nb_of_subhalos(st)
        if(tree(st,i)%level.le.1) then
          write(errunit,*) '> Error in collapse tree for :',st,i
          write(errunit,*) '> subhalo level tree(st,i)%level:',tree(st,i)%level
          stop
       end if
       remove = .false.
       if(st.lt.nsteps) then
          ison = tree(st,i)%my_sons%main_son
          if(ison.gt.0) then
             if(tree(st+1,ison)%level.le.0) then
                write(errunit,*) '> Error in collapse tree for :',st,i
                write(errunit,*) '> main son, level:',ison,tree(st+1,ison)%level
                ! son has been removed, itshouldn't occur after mod merger tree on ison
                remove = .true.
             else if(tree(st+1,ison)%my_dads%list_dads(1).ne.i) then
                ! isub is secondary progenitor of its son
                ! should occur most of the time
                remove = .true.
                !else
                ! this is a fly-by we keep the subhalo
                !remove =.false.
             end if
          else
             remove = .true.
          end if
       else
          ! last step remove all subhalos
          remove = .true.
       end if

       if(remove) then
          ! searching a halo or a unremoved subhalo where we can put it
          ! as host is tested before its subhalo it should be hostsub
          host   = tree(st,i)%hostsub
          if(host.le.0) then
             write(errunit,*) '> Error in collapse tree for :',st,i
             write(errunit,*) '> tree(st,i)%hostsub:',tree(st,i)%hostsub
             stop
          end if
          ifhost = .false.
          do while(.not.ifhost)
             if(tree(st,host)%level.eq.1) then
                ifhost = .true.
             else
                ison = tree(st,host)%my_sons%main_son
                if(ison.gt.0) then
                   ! subhalo has a decendant this is a flyby don't remove
                   ifhost = .true.
                else
                   ifhost = .true.
                end if
             end if
             if(.not.ifhost) then
                if(tree(st,host)%hostsub.le.0) then
                   write(errunit,*) '> Error in colapse tree for subhalo:',st,i
                   write(errunit,*) '> st,host,level',st,host,tree(st,host)%level
                   write(errunit,*) '> host%hostsub,host%hosthal',tree(st,host)%hostsub,tree(st,host)%hosthalo
                   stop
                end if
                host = tree(st,host)%hostsub
             end if
          end do
          call mod_structure_tree(st,i,host)
          call mod_merger_tree(st,i,host)
       end if
    end do
    
    return

  end subroutine collapse_tree

#endif


  !***********************************************************************
  ! The 6 next routines are done to remove a subhalo from the merger tree
  ! This should be call only #ifdef R_DADLESS_SUB #ifdef COL_TREE
  !***********************************************************************
  subroutine mod_structure_tree(st,isub,ihost)

    implicit none
    integer(kind=4) :: st,isub,ihost
    integer(kind=4) :: oldsub,isubtmp
    
    ! modify structure tree
    ! remove subhalo from linked list
    
    oldsub = ihost 
    do while(tree(st,oldsub)%nextsub.ne.isub)
       oldsub = tree(st,oldsub)%nextsub
       if(oldsub.le.0) stop '> subhalo not found in structure tree'
    end do
    tree(st,oldsub)%nextsub = tree(st,isub)%nextsub

    ! change removed subhalo subhaloes into host
    ! now that subhalo i is removed from the structure change data for its subhalos
   
    isubtmp = tree(st,ihost)%nextsub
    do while(isubtmp.gt.0)
       do while(isubtmp.gt.0)
          if(tree(st,isubtmp)%level.le.tree(st,ihost)%level) then
             isubtmp = -1
          else
             if(tree(st,isubtmp)%hostsub.eq.isub) tree(st,isubtmp)%hostsub = ihost
             tree(st,isubtmp)%level = tree(st,tree(st,isubtmp)%hostsub)%level + 1
             isubtmp = tree(st,isubtmp)%nextsub
          end if
       end do
    end do
    ! put subhalo i level to -1 to show it was removed change %host to host to show in which halo it was removed
    tree(st,isub)%level    = -1
    tree(st,isub)%hostsub  = ihost

    return

  end subroutine mod_structure_tree

  !***********************************************************************
  subroutine mod_merger_tree(st,isub,ihost)
  
    implicit none
    integer(kind=4) :: st,isub,ihost
    integer(kind=4) :: ifath,idfath,ison,idson,old_son,new_son 
    type(father),pointer :: fathsub,fathson
    type(son),pointer    :: sonsub

    if(tree(st,isub)%level.gt.0) then
       write(errunit,*) '> Error in mod_merger_tree for',st,isub
       write(errunit,*) '> level',tree(st,isub)%level
       write(errunit,*) '> subhalo has not been removed from the halo structure tree'
       stop
    end if

    fathsub => tree(st,isub)%my_fathers
    if(fathsub%nb_fathers.gt.0.and.st.gt.1) then
       ! give isub fathers to ihost
       call merge_fathers_list(st,isub,ihost)
       ! We have delt with the non simplified merger tree 
       ! now we nead to correct the simplified merger tree
       
       ! reset mass_fathers to 0. so that isub cannot be the main son of its fathers
       fathsub%mass_fathers(1:fathsub%nb_fathers) = 0.
       
       do ifath = 1,fathsub%nb_fathers
          ! what have changed with merge_fathers_list
          ! for each fathers of isub we either have a new progenitor for ihost 
          ! or the corresponding mass_fathers for ihost has changed
          ! meaning that the main_son of any progenitor of isub might be different
          idfath = fathsub%list_fathers(ifath)
          if(idfath.gt.0) then
             old_son = tree(st-1,idfath)%my_sons%main_son
             if(st-1.eq.stbugson.and.idfath.eq.ihbugson) then
                write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
                write(errunit,*) '> calling compute_son for:',st-1,idfath
             end if
             call compute_son(st,st-1,idfath)
             new_son = tree(st-1,idfath)%my_sons%main_son
             if(st-1.eq.stbugson.and.idfath.eq.ihbugson) then
                write(errunit,*) '> old_son,new_son:',old_son,new_son
             end if
             if(new_son.eq.isub) then 
                write(errunit,*) '> Error in mod_merger_tree for:',st,isub,ihost
                write(errunit,*) '> st-1,idfath,old_son,new_son:', st-1,idfath,old_son,new_son
                write(errunit,*) '> new_son is isub shouldn''t happen'
                stop
             end if
             if(old_son.ne.new_son) then
                if(old_son.gt.0.and.old_son.ne.isub) then
                   if(st.eq.stbugdad.and.old_son.eq.ihbugdad) then
                      write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
                      write(errunit,*) '> calling compute_dads for:',st,old_son
                   end if
                   call compute_dads(st,st-1,old_son)
                end if
                if(new_son.gt.0) then
                   if(st.eq.stbugdad.and.new_son.eq.ihbugdad) then
                      write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
                      write(errunit,*) '> calling compute_dads for:',st,new_son
                   end if
                   call compute_dads(st,st-1,new_son)
                end if
             end if
          end if
       end do
       ! ok now that we redone the simplification reset isub fathers
       fathsub%nb_fathers = 0
       deallocate(fathsub%list_fathers,fathsub%mass_fathers)
       ! recompute dads of isub
       if(st.eq.stbugdad.and.isub.eq.ihbugdad) then
          write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
          write(errunit,*) '> calling compute_dads for:',st,isub
       end if
       call compute_dads(st,st-1,isub)
       if(tree(st,isub)%my_dads%nb_dads.gt.0) then
          write(errunit,*) '> Error in mod_merger_tree for:',st,isub,ihost
          write(errunit,*) '> tree(st,isub)%my_dads%nb_dads:',tree(st,isub)%my_dads%nb_dads
          write(errunit,*) '> isub should have no dad anymore'
          stop
       end if
       ! recompute dads for ihost
       if(st.eq.stbugdad.and.ihost.eq.ihbugdad) then
          write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
          write(errunit,*) '> calling compute_dads for:',st,ihost
       end if
       call compute_dads(st,st-1,ihost)
    end if

    
    sonsub => tree(st,isub)%my_sons
    if(sonsub%nb_sons.gt.0.and.st.lt.nsteps) then
       ! give isub sons to host
       call merge_sons_list(st,isub,ihost)
       ! We have delt with the non simplified merger tree 
       ! now we nead to correct the simplified merger tree
       do ison = 1,sonsub%nb_sons
          ! what have changed with merge_sons_list
          ! for each sons of isub isub is no longer in their list_fathers 
          ! and mass_father may have changed
          ! so the main_son of each father of isub migh be a different one 
          ! the order of progenitor may have changed as well if not the progenitors themselves 
          idson = sonsub%list_sons(ison)
          fathson => tree(st+1,idson)%my_fathers
          if(fathson%nb_fathers.gt.0) then
             do ifath = 1, fathson%nb_fathers
                idfath = fathson%list_fathers(ifath)
                if(idfath.gt.0.and.idfath.ne.ihost) then
                   old_son = tree(st,idfath)%my_sons%main_son
                   if(st.eq.stbugson.and.idfath.eq.ihbugson) then
                      write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
                      write(errunit,*) '> calling compute_son for:',st,idfath
                   end if
                   call compute_son(st+1,st,idfath)
                   new_son = tree(st,idfath)%my_sons%main_son
                   if(st.eq.stbugson.and.idfath.eq.ihbugson) then
                      write(errunit,*) '> old_son,new_son:',old_son,new_son
                   end if
                   if(old_son.gt.0.and.old_son.ne.new_son) then 
                      ! main son has changed
                      if(st+1.eq.stbugdad.and.old_son.eq.ihbugdad) then
                         write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
                         write(errunit,*) '> calling compute_dads for:',st+1,old_son
                      end if
                      call compute_dads(st+1,st,old_son)
                   end if
                   if(new_son.gt.0) then
                      ! order of dads might differ
                      if(st+1.eq.stbugdad.and.new_son.eq.ihbugdad) then
                         write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
                         write(errunit,*) '> calling compute_dads for:',st+1,new_son
                      end if
                      call compute_dads(st+1,st,new_son)
                   end if
                end if
             end do
          end if
       end do
       ! remove all sons from isub
       sonsub%nb_sons = 0
       deallocate(sonsub%list_sons)
       ! recompute son of isub
       old_son = sonsub%main_son
       if(st.eq.stbugson.and.isub.eq.ihbugson) then
          write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
          write(errunit,*) '> calling compute_son for:',st,isub
       end if
       call compute_son(st+1,st,isub)
       if(st.eq.stbugson.and.isub.eq.ihbugson) then
          new_son = sonsub%main_son
          write(errunit,*) '> old_son,new_son:',old_son,new_son
       end if

       if(sonsub%main_son.gt.0) then
          write(errunit,*) '> Error in mod_merger_tree for:',st,isub,ihost
          write(errunit,*) '> sonsub%main_son:',sonsub%main_son
          write(errunit,*) '> isub should have no main_son anymore'
          stop
       end if
       if(old_son.gt.0) then
          ! recompute isub former main_son dads 
          if(st+1.eq.stbugdad.and.old_son.eq.ihbugdad) then
             write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
             write(errunit,*) '> calling compute_dads for:',st+1,old_son
          end if
          call compute_dads(st+1,st,old_son)
       end if
       ! recompute son of ihost
       old_son = tree(st,ihost)%my_sons%main_son
       if(st.eq.stbugson.and.ihost.eq.ihbugson) then
          write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
          write(errunit,*) '> calling compute_son for:',st,ihost
       end if
       call compute_son(st+1,st,ihost)
       new_son = tree(st,ihost)%my_sons%main_son
       if(st.eq.stbugson.and.ihost.eq.ihbugson) then
          write(errunit,*) '> old_son,new_son:',old_son,new_son
       end if

       if(old_son.gt.0.and.old_son.ne.new_son) then 
          ! main son has changed
          if((st+1.eq.stbugdad.and.old_son.eq.ihbugdad).or.(st.eq.stbugson.and.ihost.eq.ihbugson)) then
             write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
             write(errunit,*) '> calling compute_dads for:',st+1,old_son
          end if
          call compute_dads(st+1,st,old_son)
       end if
       if(new_son.gt.0) then
          ! order of dads might differ
          if((st+1.eq.stbugdad.and.new_son.eq.ihbugdad).or.(st.eq.stbugson.and.ihost.eq.ihbugson)) then
             write(errunit,*) '> In mod_merger_tree for:',st,isub,ihost
             write(errunit,*) '> calling compute_dads for:',st+1,new_son
          end if
          call compute_dads(st+1,st,new_son)
       end if
    end if
    
    return
    
  end subroutine mod_merger_tree

  !***********************************************************************
  subroutine merge_fathers_list(st,ih1,ih2) 
    
    implicit none
    integer(kind=4)             :: st,ih1,ih2
    type(father),pointer        :: fath1,fath2
    integer(kind=4)             :: ifath1,ifath2,idfath1,idfath2
    integer(kind=4)             :: nfath,nfath_plus
    integer(kind=4),allocatable :: listfathtmp(:)
    real(kind=8), allocatable   :: massfathtmp(:)
    logical                     :: search
    
    fath1  => tree(st,ih1)%my_fathers
    fath2  => tree(st,ih2)%my_fathers

    if(fath1%nb_fathers.le.0) return

    nfath      = fath2%nb_fathers
    nfath_plus = 0
    
    if(st.eq.stbugdad.and.(ih1.eq.ihbugdad.or.ih2.eq.ihbugdad)) then
       write(errunit,*) '> in merge_fathers_list',st,ih1,ih2
       write(errunit,*) '> fath1%list_fathers',fath1%nb_fathers
       do ifath1  = 1, fath1%nb_fathers
          write(errunit,*) '>',fath1%list_fathers(ifath1),real(fath1%mass_fathers(ifath1),4)
       end do
       if(fath2%nb_fathers.gt.0) then
          write(errunit,*) '> fath2%list_fathers',fath2%nb_fathers
          do ifath2  = 1, fath2%nb_fathers
             write(errunit,*) '>',fath2%list_fathers(ifath2),real(fath2%mass_fathers(ifath2),4)
          end do
       end if
    end if

    do ifath1  = 1, fath1%nb_fathers
       idfath1 = fath1%list_fathers(ifath1)
       ! search for father idfath1 in fath2 list
       ifath2  = 1
       idfath2 = -1
       search = .true.
       do while(search)
          if(ifath2.le.fath2%nb_fathers) then
             idfath2 = fath2%list_fathers(ifath2)
             if(idfath2.eq.idfath1) then
                ! found it cancel search
                search = .false.
             elseif(ifath2.lt.fath2%nb_fathers)then
                ! look further if possible
                ifath2 = ifath2 + 1
             else
                ! end of fath2 list end search
                search = .false.
             end if
          end if
       end do
       if(idfath2.ne.idfath1) then
          ! fath2 should contain one more progenitor
          nfath_plus = nfath_plus + 1
          nfath      = nfath + 1
       end if
    end do

    if(nfath.ne. fath2%nb_fathers+nfath_plus) then
       write(errunit,*) '> Error in merge_fathers_list init check for:',st,ih1,ih2
       write(errunit,*) '> nfath,fath2%nb_fathers+nfath_plus:',nfath,fath2%nb_fathers+nfath_plus
       stop
    end if

    if(nfath_plus.gt.0) then
       if(fath2%nb_fathers.gt.0) then
          ! create a copy and reallocate fath2%list_fathers and fath2%mass_fathers
          allocate(listfathtmp(1:fath2%nb_fathers),massfathtmp(1:fath2%nb_fathers))
          listfathtmp(1:fath2%nb_fathers) = fath2%list_fathers(1:fath2%nb_fathers)
          massfathtmp(1:fath2%nb_fathers) = fath2%mass_fathers(1:fath2%nb_fathers)
          deallocate(fath2%list_fathers,fath2%mass_fathers)
       end if
       allocate(fath2%list_fathers(nfath),fath2%mass_fathers(nfath))
       fath2%list_fathers(1:nfath) = -1
       fath2%mass_fathers(1:nfath) = 0.
       if(fath2%nb_fathers.gt.0) then
          fath2%list_fathers(1:fath2%nb_fathers) = listfathtmp(1:fath2%nb_fathers)
          fath2%mass_fathers(1:fath2%nb_fathers) = massfathtmp(1:fath2%nb_fathers)
       end if
    end if
    
    nfath_plus = 0
    do ifath1  = 1, fath1%nb_fathers
       idfath1 = fath1%list_fathers(ifath1)
       ! search for father idfath1 in fath2 list
       ifath2  = 1
       idfath2 = -1
       search = .true.
       do while(search)
          if(ifath2.le.fath2%nb_fathers) then
             idfath2 = fath2%list_fathers(ifath2)
             if(idfath2.eq.idfath1) then
                ! found it cancel search
                search = .false.
             elseif(ifath2.lt.fath2%nb_fathers) then
                ! look further if possible
                ifath2 = ifath2 + 1
             else
                ! end of fath2 list end search
                ifath2  = -1
                idfath2 = -1
                search  = .false.
             end if
          end if
       end do
       if(idfath2.ne.idfath1) then
          ! fath2 should contain one more progenitor
          nfath_plus = nfath_plus + 1
          ifath2     = fath2%nb_fathers + nfath_plus
          ! add progenitor to list
          if(ifath2.le.0.or.ifath2.gt.nfath) then
             write(errunit,*) '> Error in merger_fathers_list for:',st,ih1,ih2
             write(errunit,*) '> ifath2,nfath:',ifath2,nfath
             stop
          end if
          fath2%list_fathers(ifath2) = idfath1 
       end if
       ! add masses
       fath2%mass_fathers(ifath2) = fath2%mass_fathers(ifath2) + fath1%mass_fathers(ifath1)
    end do
    if(fath2%nb_fathers + nfath_plus.eq.nfath) then
       fath2%nb_fathers = nfath
    else
        write(errunit,*) '> Error in merger_fathers_list,final check for:',st,ih1,ih2
        write(errunit,*) '> fath2%nb_fathers + nfath_plus,nfath:',fath2%nb_fathers + nfath_plus,nfath
        stop
    end if
    
    if(st.eq.stbugdad.and.(ih1.eq.ihbugdad.or.ih2.eq.ihbugdad)) then
       if(fath2%nb_fathers.gt.0) then
          write(errunit,*) '> fath2%list_fathers',fath2%nb_fathers
          do ifath2  = 1, fath2%nb_fathers
             write(errunit,*) '>',fath2%list_fathers(ifath2),fath2%mass_fathers(ifath2)
          end do
       end if
    end if

    do ifath1  = 1, fath1%nb_fathers
       idfath1 = fath1%list_fathers(ifath1)
       if(idfath1.gt.0) then
          ! now deal with ih1 fathers to replace ih1 by ih2 in their son lists
          call replace_son(st-1,idfath1,ih1,ih2)
       end if
    end do
    
    return
    
  end subroutine merge_fathers_list

  !***********************************************************************
  subroutine merge_sons_list(st,ih1,ih2) 
    
    implicit none
    integer(kind=4)             :: st,ih1,ih2
    type(son),pointer           :: son1,son2
    integer(kind=4)             :: ison1,ison2,idson1,idson2
    integer(kind=4)             :: nson,nson_plus
    integer(kind=4),allocatable :: listsontmp(:)
    logical                     :: search

    son1   => tree(st,ih1)%my_sons
    son2   => tree(st,ih2)%my_sons
    
    if(son1%nb_sons.le.0) return

    if(st.eq.stbugson.and.(ih1.eq.ihbugson.or.ih2.eq.ihbugson)) then
       write(errunit,*) '> in merge_sons_list',st,ih1,ih2
       write(errunit,*) '> son1%list_sons',son1%nb_sons
       do ison1  = 1, son1%nb_sons
          write(errunit,*) '>',son1%list_sons(ison1)
       end do
       if(son2%nb_sons.gt.0) then
          write(errunit,*) '> son2%list_sons',son2%nb_sons
          do ison2  = 1, son2%nb_sons
             write(errunit,*) '>',son2%list_sons(ison2)
          end do
       end if
    end if

    nson      = son2%nb_sons
    nson_plus = 0
    
    if(st.eq.stbugson.and.(ih1.eq.ihbugson.or.ih2.eq.ihbugson)) then
       write(errunit,*) '> idson1,idson2,nson_plus'
    end if
    do ison1  = 1, son1%nb_sons
       idson1 = son1%list_sons(ison1)
       ! search for son idson1 in son2 list
       ison2  = 1
       idson2 = -1
       search = .true.
       do while(search)
          if(ison2.le.son2%nb_sons) then
             idson2 = son2%list_sons(ison2)
             if(idson2.eq.idson1) then
                ! found it cancel search
                search = .false.
             elseif(ison2.lt.son2%nb_sons) then
                ! look further if possible
                ison2 = ison2 + 1
             else
                ! end of son2 list end search
                search = .false.
             end if
          end if
       end do

       if(idson2.ne.idson1) then
          ! son2 should contain one more progenitor
          nson_plus = nson_plus + 1
          nson      = nson + 1
       end if
       if(st.eq.stbugson.and.(ih1.eq.ihbugson.or.ih2.eq.ihbugson)) then
          write(errunit,*) '>',idson1,idson2,nson_plus
       end if
    end do

    if(st.eq.stbugson.and.(ih1.eq.ihbugson.or.ih2.eq.ihbugson)) then
       write(errunit,*) '> nb of son added:',nson_plus
       write(errunit,*) '> total nb of sons:',nson
    end if

    if(nson.ne. son2%nb_sons+nson_plus) then
       write(errunit,*) '> Error in merge_sons_list init check'
       write(errunit,*) '> nson,son2%nb_sons+nson_plus:',nson,son2%nb_sons+nson_plus
       stop
    end if

    if(nson_plus.gt.0) then
       if(son2%nb_sons.gt.0) then
          ! create a copy and reallocate son2%list_sons
          allocate(listsontmp(1:son2%nb_sons))
          listsontmp(1:son2%nb_sons) = son2%list_sons(1:son2%nb_sons)
          deallocate(son2%list_sons)
       end if
       allocate(son2%list_sons(nson))
       son2%list_sons(1:nson) = -1
       if(son2%nb_sons.gt.0) then
          son2%list_sons(1:son2%nb_sons) = listsontmp(1:son2%nb_sons)
       end if
    end if
    
    nson_plus = 0
    do ison1  = 1, son1%nb_sons
       idson1 = son1%list_sons(ison1)
       ! search for son idson1 in son2 list
       ison2  = 1
       idson2 = -1
       search = .true.
       do while(search)
          if(ison2.le.son2%nb_sons) then
             idson2 = son2%list_sons(ison2)
             if(idson2.eq.idson1) then
                ! found it cancel search
                search = .false.
             elseif(ison2.lt.son2%nb_sons) then
                ! look further if possible
                ison2 = ison2 + 1
             else
                ! end of son2 list end search
                ison2  = -1
                idson2 = -1
                search = .false.
             end if
          end if
       end do
       if(idson2.ne.idson1) then
          ! son2 should contain one more progenitor
          nson_plus = nson_plus + 1
          ison2     = son2%nb_sons + nson_plus
          ! add progenitor to list
          if(ison2.le.0.or.ison2.gt.nson) then
             write(errunit,*) '> Error in merger_sons_list'
             write(errunit,*) '> ison2,nson:',ison2,nson
             stop
          end if
          son2%list_sons(ison2) = idson1 
       end if
    end do

    if(son2%nb_sons + nson_plus.eq.nson) then
       son2%nb_sons = nson
    else
        write(errunit,*) '> Error in merger_sons_list,final check'
        write(errunit,*) '> son2%nb_sons + nson_plus,nson:',son2%nb_sons + nson_plus,nson
        stop
    end if
       
    if(st.eq.stbugson.and.(ih1.eq.ihbugson.or.ih2.eq.ihbugson)) then
       if(son2%nb_sons.gt.0) then
          write(errunit,*) '> son2%list_sons',son2%nb_sons
          do ison2  = 1, son2%nb_sons
             write(errunit,*) '>',son2%list_sons(ison2)
          end do
       end if
    end if


    do ison1  = 1, son1%nb_sons
       idson1 = son1%list_sons(ison1)
       if(idson1.gt.0) then
          ! now deal with ih1 sons to replace ih1 by ih2 in their father lists
          call replace_fath(st+1,idson1,ih1,ih2)
       end if
    end do

    
    return
    
  end subroutine merge_sons_list

  !***********************************************************************
  subroutine replace_fath(st,idson,ih1,ih2)
  
    implicit none
    integer(kind=4)             :: st,idson,ih1,ih2
    type(father),pointer        :: fathson
    integer(kind=4)             :: ifath,idfath,ifath1,ifath2,ifathtmp
    integer(kind=4),allocatable :: listfathtmp(:)
    real(kind=8), allocatable   :: massfathtmp(:)
    
    fathson => tree(st,idson)%my_fathers
    if(fathson%nb_fathers.le.0) return
    
    ifath1 = -1
    ifath2 = -1
    ! search for ih1 and ih2
    do ifath = 1,fathson%nb_fathers
       idfath = fathson%list_fathers(ifath)
       if(idfath.eq.ih1) ifath1 = ifath
       if(idfath.eq.ih2) ifath2 = ifath
    end do
    
    if(ifath1.le.0) then
       write(errunit,*) '> Error in replace_fath for:',st,idson,ih1,ih2
       write(errunit,*) '> could not find ih1 among fathers'
       write(errunit,*) '> list_fathers,mass_fathers'
       do ifath = 1,fathson%nb_fathers
          write(errunit,*) '>',fathson%list_fathers(ifath),fathson%mass_fathers(ifath)
       end do
       stop
    end if
    
    if(ifath2.le.0) then
       ! replace ih1 by ih2 for slot ifath1
       fathson%list_fathers(ifath1) = ih2
    else
       ! copy, realocate and merge ih1 and ih2 in the fathers list
       allocate(listfathtmp(fathson%nb_fathers),massfathtmp(fathson%nb_fathers))
       listfathtmp(1:fathson%nb_fathers) = fathson%list_fathers(1:fathson%nb_fathers)
       massfathtmp(1:fathson%nb_fathers) = fathson%mass_fathers(1:fathson%nb_fathers)
       deallocate(fathson%list_fathers,fathson%mass_fathers)
       allocate(fathson%list_fathers(fathson%nb_fathers-1),fathson%mass_fathers(fathson%nb_fathers-1))
       ifathtmp                               = 0
       fathson%list_fathers(1:fathson%nb_fathers-1) = -1
       fathson%mass_fathers(1:fathson%nb_fathers-1) = 0.
       do ifath = 1,fathson%nb_fathers
          if(ifath.ne.ifath1) then
             ifathtmp = ifathtmp + 1
             if(ifathtmp.ge.fathson%nb_fathers) then
                write(errunit,*) '> Error in replace_fath for:',st,idson,ih1,ih2
                write(errunit,*) '> ifathtmp,fathson%nb_fathers-1:',ifathtmp,fathson%nb_fathers-1
                stop
             end if
             fathson%list_fathers(ifathtmp) = listfathtmp(ifath)
             fathson%mass_fathers(ifathtmp) = massfathtmp(ifath)
             if(ifath.eq.ifath2) then
                if(ifathtmp.lt.ifath) then
                   ifath2  = ifathtmp
                elseif(ifathtmp.gt.ifath) then
                   write(errunit,*) '> Error in replace_fath for:',st,idson,ih1,ih2
                   write(errunit,*) '> ifathtmp,ifath:',ifathtmp,ifath
                   write(errunit,*) '> ifathtmp.gt.ifath, this should not happen'
                   stop
                end if
             end if
          end if
       end do
       ! adding masses
       fathson%mass_fathers(ifath2) = fathson%mass_fathers(ifath2) + massfathtmp(ifath1)
       deallocate(listfathtmp,massfathtmp)
       if(ifathtmp.ne.fathson%nb_fathers-1) then
          write(errunit,*) '> Error in replace_fath for:',st,idson,ih1,ih2
          write(errunit,*) '> ifathtmp,fathson%nb_fathers-1:',ifathtmp,fathson%nb_fathers-1
          stop
       end if
       fathson%nb_fathers = fathson%nb_fathers - 1
    end if
    
    return
    
  end subroutine replace_fath

  !***********************************************************************
  subroutine replace_son(st,idfath,ih1,ih2)
  
    implicit none

    integer(kind=4)             :: st,idfath,ih1,ih2
    type(son),pointer           :: sonfath
    integer(kind=4)             :: ison,idson,ison1,ison2,isontmp
    integer(kind=4),allocatable :: listsontmp(:)
    
    sonfath => tree(st,idfath)%my_sons
    if(sonfath%nb_sons.le.0) return
    
    ison1 = -1
    ison2 = -1
    ! search for ih1 and ih2
    do ison = 1,sonfath%nb_sons
       idson = sonfath%list_sons(ison)
       if(idson.eq.ih1) ison1 = ison
       if(idson.eq.ih2) ison2 = ison
    end do
    
    if(ison1.le.0) then
       write(errunit,*) '> Error in replace_son for:',st,idson,ih1,ih2
       write(errunit,*) '> could not find ih1 among sons'
       write(errunit,*) '> list_sons'
       do ison = 1,sonfath%nb_sons
          write(errunit,*) '>',sonfath%list_sons(ison)
       end do
       stop
    end if
    
    if(ison2.le.0) then
       ! replace ih1 by ih2 for slot ison1
       sonfath%list_sons(ison1) = ih2
    else
       ! copy, realocate and merge ih1 and ih2 in the sons list
       allocate(listsontmp(sonfath%nb_sons))
       listsontmp(1:sonfath%nb_sons) = sonfath%list_sons(1:sonfath%nb_sons)
       deallocate(sonfath%list_sons)
       allocate(sonfath%list_sons(sonfath%nb_sons-1))
       isontmp                        = 0
       sonfath%list_sons(1:sonfath%nb_sons-1) = -1
       do ison = 1,sonfath%nb_sons
          if(ison.ne.ison1) then
             isontmp = isontmp + 1
             if(isontmp.ge.sonfath%nb_sons) then
                write(errunit,*) '> Error in replace_son for:',st,idson,ih1,ih2
                write(errunit,*) '> isontmp,sonfath%nb_sons-1:',isontmp,sonfath%nb_sons-1
                stop
             end if
             sonfath%list_sons(isontmp) = listsontmp(ison)
          end if
       end do
       deallocate(listsontmp)
       if(isontmp.ne.sonfath%nb_sons-1) then
          write(errunit,*) '> Error in replace_son for:',st,idson,ih1,ih2
          write(errunit,*) '> isontmp,sonfath%nb_sons-1:',isontmp,sonfath%nb_sons-1
          stop
       end if
       sonfath%nb_sons = sonfath%nb_sons - 1
    end if
    
    return

  end subroutine replace_son

  !***********************************************************************
  subroutine detect_fragment(st,i)
    
    ! I'm following the definition form GalaxyMaker 

    implicit none
    integer(kind=4)      :: st,i
    integer(kind=4)      :: idad,iddad,ifath,idfath
    real(kind=8)         :: mfrag,mtot
    type(dad),pointer    :: dads
    type(father),pointer :: fathers
    
    
    if(st.eq.1) then
       ! there is no fragments at the initial step
       tree(st,i)%frag = 0
       return
    end if
    
    if(tree(st,i)%level.le.0) then
       ! a removed subhalo is defined as a fragment
       tree(st,i)%frag = 1
       if(tree(st,i)%my_dads%nb_dads.gt.0) then
          write(errunit,*) '> Error in detect fragment for halo',st,i
          write(errunit,*) '> level:',tree(st,i)%level
          write(errunit,*) '> tree(st,i)%my_dads%nb_dads:',tree(st,i)%my_dads%nb_dads
          do idad = 1,tree(st,i)%my_dads%nb_dads
             iddad = tree(st,i)%my_dads%list_dads(idad)
             write(errunit,*) '>',idad, tree(st,i)%my_dads%mass_dads(idad),tree(st-1,iddad)%level
          end do
          stop
       end if
       if(tree(st,i)%my_sons%main_son.gt.0) then
          write(errunit,*) '> Error in detect fragment for halo',st,i
          write(errunit,*) '> level:',tree(st,i)%level
          write(errunit,*) '> main_son:',tree(st,i)%my_sons%main_son
          stop
       end if
       return
    end if

    if(tree(st,i)%my_dads%nb_dads.gt.0) then
       ! halo has at least one progenitor
       dads  => tree(st,i)%my_dads
       tree(st,i)%frag = 1
       do idad = 1,dads%nb_dads
          ! if one of its progenitors is not a fragment the halo is not a fragment
          iddad = dads%list_dads(idad)
          if(tree(st-1,iddad)%frag.eq.0) tree(st,i)%frag = 0
          tree(st-1,iddad)%my_sons%nb_sons =  tree(st-1,iddad)%my_sons%nb_sons + 1
          if(tree(st-1,iddad)%my_sons%nb_sons.gt.1) then
             write(errunit,*) '> Error in detect frag for halo:',st-1,iddad
             write(errunit,*) '> Halo has more than one son'
             write(errunit,*) '> tree(st-1,iddad)%my_sons%main_son,i',tree(st-1,iddad)%my_sons%main_son,i
             stop
          end if
       end do
    else
       if(tree(st,i)%level.gt.1) then
          ! a subhalo with no progenitor is a fragment
          tree(st,i)%frag = 1
          return
       end if
       ! halo has no progenitors, compute the mass ratio obtained from other halos
       mfrag = 0.0d0 
       mtot  = 0.0d0
       fathers => tree(st,i)%my_fathers
       if(fathers%nb_fathers.le.0) then
          write(errunit,*) '> Error in detect_fragment for halo:',st,i
          write(errunit,*) '> faths%nb_fathers:',fathers%nb_fathers
          stop
       end if
       do ifath = 1,fathers%nb_fathers
          idfath = fathers%list_fathers(ifath)
          if(idfath.gt.0) then
             mfrag = mfrag + fathers%mass_fathers(ifath)
          end if
          mtot = mtot + fathers%mass_fathers(ifath)
       end do
       if(mfrag/mtot.gt.0.5d0) then
          ! more than half the mass acounted for comes from others halos, this is a fragment
          tree(st,i)%frag = 1
       else
          tree(st,i)%frag = 0
       end if
    end if
    
    return

  end subroutine detect_fragment

  !***********************************************************************
  subroutine deallocate_old_tree(st,i)
    
    implicit none
    integer(kind=4) :: st,i
    type(father),pointer :: faths
    type(son),pointer    :: sons 


    faths => tree(st,i)%my_fathers
    sons  => tree(st,i)%my_sons
    if(faths%nb_fathers.gt.0) then
       deallocate(faths%list_fathers)
       deallocate(faths%mass_fathers)
       faths%nb_fathers = 0
    end if
    if(sons%nb_sons.gt.0) then
       deallocate(sons%list_sons)
       sons%nb_sons = 0
    end if
    
    return
    
  end subroutine deallocate_old_tree

end module simpl_merger_tree

!*****************************************************************************************************************
subroutine indexx(n,arr,indx)
  
  ! returns indices indx of an array arr such that arr(indx) is in ascending order (chap 8) 
  
  implicit none 
  
  integer(kind=4)            :: n,indx(n)
  real(kind=8)               :: arr(n)
  integer(kind=4), parameter :: m = 7, nstack = 50
  integer(kind=4)            :: i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
  real(kind=8)               :: a
  
  do j=1,n
     indx(j) = j
  enddo
  
  jstack = 0
  l      = 1
  ir     = n
1 if (ir-l < m) then
     
     do j=l+1,ir
        indxt = indx(j)
        a     = arr(indxt)
        do i=j-1,1,-1
           if (arr(indx(i)) <= a) goto 2
           indx(i+1) = indx(i)
        enddo
        i         = 0
2       indx(i+1) = indxt
     enddo
     if (jstack == 0) return
     ir     = istack(jstack)
     l      = istack(jstack-1)
     jstack = jstack-2
     
  else
     
     k         = (l+ir)/2
     itemp     = indx(k)
     indx(k)   = indx(l+1)
     indx(l+1) = itemp
     if (arr(indx(l+1)) > arr(indx(ir))) then
        itemp     = indx(l+1)
        indx(l+1) = indx(ir)
        indx(ir)  = itemp
     endif
     if (arr(indx(l)) > arr(indx(ir))) then
        itemp    = indx(l)
        indx(l)  = indx(ir)
        indx(ir) = itemp
     endif
     if (arr(indx(l+1)) > arr(indx(l))) then
        itemp     = indx(l+1)
        indx(l+1) = indx(l)
        indx(l)   = itemp
     endif
     i     = l+1
     j     = ir
     indxt = indx(l)
     a     = arr(indxt)
3    continue
     i     = i+1
     if (arr(indx(i)) < a) goto 3
4    continue
     j     = j-1
     if (arr(indx(j)) > a) goto 4
     if (j < i) goto 5
     itemp   = indx(i)
     indx(i) = indx(j)
     indx(j) = itemp
     goto 3
5    indx(l) = indx(j)
     indx(j) = indxt
     jstack  = jstack+2
     if (jstack > nstack) then
        write(*,*) 'nstack too small in indexx'
        read *
     endif
     if (ir-i+1 >= j-l) then
        istack(jstack)   = ir
        istack(jstack-1) = i
        ir               = j-1
     else
        istack(jstack)   = j-1
        istack(jstack-1) = l
        l                = i
     endif
     
  endif
  
  goto 1
  
end subroutine indexx

!*****************************************************************************************************************
subroutine rd_sort(n,arr,indx)

  ! as indexx doesn't seem to work for me, might be stupid I'm writing my own routine 
  
  implicit none
  integer(kind=4) :: n,indx(n)
  real(kind=8)    :: arr(n)
  real(kind=4)    :: sarr(n),amin,amax
  integer(kind=4) :: j,jmin,jmax,jlocmin,jlocmax,i1,i2,i3,i4,n_left

  do j = 1, n
     sarr(j)     = arr(j)
     indx(j)     = j
  end do
 
  n_left = n
  jmin   = 1
  jmax   = n
!  print*,"arr",arr
  do while(n_left .gt. 1) 
     ! sorting by both end
     amax    = minval(sarr)
     amin    = maxval(sarr)
     jlocmin = jmax
     jlocmax = jmin
     do j = jmin,jmax
        if(sarr(j).gt.amax) then
           amax    = sarr(j)
           jlocmax = j
        end if
        if(sarr(j).lt.amin) then
           amin    = sarr(j)
           jlocmin = j
        end if
     end do
     if(jlocmax.ne.jmin) then
        i1 = indx(jmin)
        i2 = indx(jlocmax)
        ! switch jmin with jlocmax
        indx(jmin)    = i2
        sarr(jmin)    = arr(i2)
        indx(jlocmax) = i1
        sarr(jlocmax) = arr(i1)
        if(jmin.eq.jlocmin) then
           jlocmin = jlocmax
        end if
     end if
     if(jlocmin.ne.jmax) then
        i3 = indx(jlocmin)
        i4 = indx(jmax)
        ! switch 
        indx(jmax)    = i3
        sarr(jmax)    = arr(i3)
        indx(jlocmin) = i4
        sarr(jlocmin) = arr(i4)
     end if
     jmin   = jmin + 1
     jmax   = jmax - 1
     n_left = n_left - 2
  end do

  return 

end subroutine rd_sort
