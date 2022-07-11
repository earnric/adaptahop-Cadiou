module jeje
  
#ifdef SIMPL
  ! JB - 2010 - bonus for SIMPL option.
  ! 1/ compute all IDs that we'll need for database or GalaxyMaker purposes
  ! 2/ bufferize outputs (write_bushes is a massive time killer, it seems).
  ! NB: requires SIMPL option from D.Tweed (i.e. requires a halo to have at most one descendent).

  use compute_tree

  public
  
  ! IDs (below) are unique in a simulation
  type jeje_treeID
     ! IDs of each halo
     integer(kind=8) :: BushID       ! unique ID of the bush the halo is in 
     integer(kind=8) :: TreeID       ! unique ID of the tree the halo is in
     integer(kind=8) :: HaloID       ! unique ID of the halo, defined recursively along its tree
     integer(kind=8) :: haloNum      ! == my_number   -> allows match with HaloMaker outputs (or halos_results.xxx files)
     integer(kind=8) :: haloTimestep ! == my_timestep (kind=8 is an overkill but simplifies outputs)
     ! Vertical connections (progs/descendents)
     integer(kind=8) :: FirstProgenitorID  
     integer(kind=4) :: firstProg ! an index, not an id
     integer(kind=8) :: NextProgenitorID
     integer(kind=4) :: nextProg
     integer(kind=8) :: DescendentID
     integer(kind=4) :: descendent
     integer(kind=8) :: LastProgenitorID
     ! Horizontal connections (main halos / subs)
     integer(kind=8) :: HostHaloID
     integer(kind=8) :: HostSubID
     integer(kind=8) :: NextSubID
  end type jeje_treeID
  integer(kind=4)             :: nHalosInBuffer        ! size of the buffer will be fixed by maximum nb of halos in a file at any timestep... 
  integer(kind=4),parameter   :: nIDsPerHalo    = 12   ! nb of fields in jeje_treeID type
  integer(kind=8),allocatable :: bufferIDs(:,:,:)

  type simple_tree
     integer(kind=4)         :: son      ! index on only son in tree(), or -1 if no son
     integer(kind=4)         :: ndads    ! number of dads
     integer(kind=4),pointer :: dads(:)  ! indexes of dads in tree() 
     integer(kind=4)         :: halonum  ! my_number in the output of HaloMaker.
  end type simple_tree

  ! indexes which are usefull to navigate between progs/descendents and halos/sub-halos in GALICS
  ! These indexes are valid within one treeFile (and hence defined just before output)
  type jejeIndexes
     integer(kind=4) :: me
     integer(kind=4) :: descendent
     integer(kind=4) :: ndads
     integer(kind=4) :: firstProg
     integer(kind=4) :: nextProg
     integer(kind=4) :: hosthalo
     integer(kind=4) :: hostsub
     integer(kind=4) :: nextsub
     integer(kind=4) :: level
  end type jejeIndexes
  integer(kind=4),parameter   :: nIndexesPerHalo    = 9   ! nb of fields in jejeIndexes type
  integer(kind=4),allocatable :: bufferIndexes(:,:,:)

  integer(kind=4),parameter   :: nPropsPerHalo = 21  ! pos(3),vel(3),mfof,mvir,tvir,rvir,cvel,macc,rfof,spin,fragFlag, L(3), ep, ek, et
  real(kind=4),allocatable    :: propsbuffer(:,:,:)

  ! container for trees so each timestep contains more or less halos ... 
  type tsTree
     type(simple_tree),pointer :: stree(:)
     type(jeje_treeID),pointer :: jtree(:)
     type(jejeIndexes),pointer :: itree(:)
  end type tsTree
  type(tsTree),allocatable :: ttree(:)

  integer(kind=4) :: progcnt   ! global counter used to define halo IDs in recursive routine below.

contains
!*****************************************************************************************************************

  subroutine buffered_outputs

    implicit none 
    
    integer(kind=4)             :: ts,ih,i
    integer(kind=4)             :: ifile
    integer(kind=4),allocatable :: nInFile(:),renum(:,:),ibuff(:)
    real(kind=4)    :: t_start,t_end
    integer(kind=4) :: istart,iend,cntrate
    
    call CPU_time(t_start)
    call system_clock(count_rate=cntrate)
    call system_clock(count=istart)
    write(errunit,*) '> Beginning buffered outputs'

    ! Count how many halos each file will contain
    ! make each tree point to its destination file according to some most stupid criterion (BushID modulo n_tree_files)
    ! NB: it is important that a full bush (many trees) is written in a single tree file to make sure all sub-halos are in the same file as their hosts... 
    ! NB : It is only now that we can perform "renumbering", ie. convert ID links to indexes valid within each treefile
    ! -> we fill up the iTree 
    allocate(nInFile(0:n_tree_files),renum(nsteps,0:n_tree_files)) ! elt 0 for halos with BushID = -1 (i.e. trees not rooted at 0) -> will not be written out... 
    nInFile = 0
    renum   = 0
    do ts = 1,nsteps
       do ih = 1,nb_of_halos(ts)+nb_of_subhalos(ts) 
          ifile = mod(ttree(ts)%jtree(ih)%BushID,n_tree_files) + 1 ! -> number going from 1 to n_tree_Files (or actually 0 if bushID = -1)
          !! if (tree(ts,ih)%frag == 1 .and. ifile > 0) print*,'fragged halo ... '
          nInFile(ifile)  = nInFile(ifile) + 1
          renum(ts,ifile) = renum(ts,ifile) + 1     ! basically number of halos per timesteps per file... 
          ttree(ts)%itree(ih)%me = renum(ts,ifile)  ! halo number in output file, at a given ts.
       end do
    end do
    write(errunit,*) 'number of discarded halos :',nInFile(0)

    ! define buffer size (nHalosInBuffer) as the max number of halos per timestep in any file.
    nHalosInBuffer = maxval(renum(:,1:n_tree_files))
    write(errunit,*) '> SIZE OF BUFFER === ',nHalosInBuffer

    ! finish (most of) renumbering
    do ts = 1,nsteps
       do ih = 1,nb_of_halos(ts)+nb_of_subhalos(ts)
          ! descendent
          i = ttree(ts)%jtree(ih)%descendent
          if (i > 0) ttree(ts)%itree(ih)%descendent = ttree(ts+1)%itree(i)%me
          ! first progenitor
          i = ttree(ts)%jtree(ih)%firstProg
          if (i > 0) ttree(ts)%itree(ih)%firstprog  = ttree(ts-1)%itree(i)%me
          ! next progenitor
          i = ttree(ts)%jtree(ih)%nextProg
          if (i > 0) ttree(ts)%itree(ih)%nextprog   = ttree(ts)%itree(i)%me
          ! host halo 
          i = tree(ts,ih)%hosthalo
          if (i > 0) ttree(ts)%itree(ih)%hosthalo   = ttree(ts)%itree(i)%me
          ! host sub
          i = tree(ts,ih)%hostsub
          if (i > 0) ttree(ts)%itree(ih)%hostsub    = ttree(ts)%itree(i)%me
          ! next sub 
          i = tree(ts,ih)%nextsub
          if (i > 0) ttree(ts)%itree(ih)%nextsub    = ttree(ts)%itree(i)%me
       end do
    end do

    ! output timestep properties
    call write_tsno_file(renum)

    ! output things via a buffer 
    allocate(bufferIDs(nIDsPerHalo,nHalosInBuffer,n_tree_files))
    allocate(bufferIndexes(nIndexesPerHalo,nHalosInBuffer,n_tree_files))
    write(errunit,*) '> Creating output tree files'
    do ifile = 1,n_tree_files
       call open_tree_file(ifile,renum(:,ifile)) 
    end do
    write(errunit,*) '> Filling up buffer and dumping ... '
    allocate(ibuff(n_tree_files))
    bufferIDs     = -1
    bufferIndexes = -1
    do ts = 1,nsteps
       ibuff = 0
       do ih = 1,nb_of_halos(ts)+nb_of_subhalos(ts) 
          ifile = mod(ttree(ts)%jtree(ih)%BushID,n_tree_files) + 1
          if (ifile > 0) then 
             ibuff(ifile) = ibuff(ifile) + 1
             bufferIDs(:,ibuff(ifile),ifile) = (/ttree(ts)%jtree(ih)%BushID,ttree(ts)%jtree(ih)%TreeID,ttree(ts)%jtree(ih)%HaloID,&
                  & ttree(ts)%jtree(ih)%HaloNum,ttree(ts)%jtree(ih)%HaloTimestep,ttree(ts)%jtree(ih)%FirstProgenitorID,&
                  & ttree(ts)%jtree(ih)%NextProgenitorID,ttree(ts)%jtree(ih)%DescendentID,ttree(ts)%jtree(ih)%LastProgenitorID,& 
                  & ttree(ts)%jtree(ih)%HostHaloID,ttree(ts)%jtree(ih)%HostSubID,ttree(ts)%jtree(ih)%NextSubID/)
             bufferIndexes(:,ibuff(ifile),ifile) = (/ ttree(ts)%itree(ih)%me, ttree(ts)%itree(ih)%descendent, &
                  & ttree(ts)%itree(ih)%ndads, ttree(ts)%itree(ih)%firstProg, ttree(ts)%itree(ih)%nextProg, & 
                  & ttree(ts)%itree(ih)%hosthalo, ttree(ts)%itree(ih)%hostsub, ttree(ts)%itree(ih)%nextsub,tree(ts,ih)%level/)
          end if
       end do
       ! quick check : 
       do ifile = 1,n_tree_files
          if (ibuff(ifile) /= renum(ts,ifile)) then 
             write(errunit,*) 'ibuff /= renum ... ',ibuff(ifile),renum(ts,ifile),ts,ifile
             stop
          end if
       end do
       call dump_buffer(ibuff)
    end do
    do i = 1,n_tree_files
       close(dumptree_unit+i)
    end do
    deallocate(bufferIDs,bufferIndexes)

    call CPU_time(t_end)
    call system_clock(count=iend)
    write(errunit,*) '> Writing out tree files took :',int(t_end-t_start),' seconds'
    write(errunit,*) '>                 (wall clock : ',(iend-istart)/cntrate

    ! We now need to break down halos_results files into small sub-files which match the tree files ... 
    write(errunit,*) '> Creating output props files '
    do ifile = 1,n_tree_files
       call open_props_file(ifile,renum(:,ifile))
    end do
    write(errunit,*) '> Filling up buffer and dumping ... '
    allocate(propsBuffer(nPropsPerHalo,nHalosInBuffer,n_tree_Files))
    ibuff       = 0
    propsbuffer = 0.0
    do ts = 1, nsteps
       call read_ts_data(ts)
       do ih = 1,nb_of_halos(ts)+nb_of_subhalos(ts)
          ifile =  mod(ttree(ts)%jtree(ih)%BushID,n_tree_files) + 1
          if (ifile > 0) then 
             ibuff(ifile) = ibuff(ifile) + 1
             propsbuffer(:,ibuff(ifile),ifile) = (/liste_halos(ih)%p%x,liste_halos(ih)%p%y, liste_halos(ih)%p%z, &
                  & liste_halos(ih)%v%x,liste_halos(ih)%v%y,liste_halos(ih)%v%z, & 
                  & liste_halos(ih)%m,liste_halos(ih)%r,liste_halos(ih)%spin,    & 
                  & liste_halos(ih)%datas%rvir,liste_halos(ih)%datas%mvir,liste_halos(ih)%datas%tvir, & 
                  & liste_halos(ih)%datas%cvel,real(liste_halos(ih)%macc,4),real(tree(ts,ih)%frag,4), & 
                  & liste_halos(ih)%L%x, liste_halos(ih)%L%y, liste_halos(ih)%L%z,& 
                  & liste_halos(ih)%ep, liste_halos(ih)%ek, liste_halos(ih)%et/)
          end if
       end do
       ! quick check
       do ifile = 1,n_tree_files
          if (ibuff(ifile) /= renum(ts,ifile)) then 
             write(errunit,*) 'ibuff /= renum ... '
             stop
          end if
       end do
       call dump_props_buffer(ibuff)
    end do
    do i = 1,n_tree_files
       close(dumptree_unit+i)
    end do
    deallocate(nInFile,renum,ibuff,propsbuffer)

    call CPU_time(t_start)
    call system_clock(count=istart)
    write(errunit,*) '> Writing out tree props took :',int(t_start-t_end),' seconds'
    write(errunit,*) '>                 (wall clock : ',(istart-iend)/cntrate

    return

    contains 
      
      subroutine write_tsno_file(renum)
        implicit none 
        integer(kind=4) :: ifile,renum(nsteps,0:n_tree_files),unit  
        character(512)  :: file
        do ifile = 1,n_tree_files
           write(file,'(a,a,i3.3,a,i3.3)') trim(data_dir),'tstep_file_',nsteps,'.',ifile
           unit = dumptree_unit + ifile
           open(unit=unit,file=file,status="unknown",form="unformatted")
           write(unit) nsteps
           write(unit) renum(1:nsteps,ifile)   ! nb of halos per timestep in each file
           write(unit) aexp(1:nsteps)          ! exp factor
           write(unit) age_univ(1:nsteps)      ! ... 
           close(unit)
        end do
      end subroutine write_tsno_file

      subroutine open_tree_file(ifile,nhperts)
        implicit none
        integer(kind=4) :: ifile,unit,nhperts(nsteps)
        character(512)  :: file
        write(file,'(a,a,i3.3,a,i3.3)') trim(data_dir),'tree_file_',nsteps,'.',ifile
        unit = dumptree_unit + ifile
        open(unit=unit,file=file,status="unknown",form="unformatted")
        write(unit) nsteps,nIDsPerHalo,nIndexesPerHalo
        write(unit) nhperts(1:nsteps)
      end subroutine open_tree_file

      subroutine dump_buffer(ibuff)
        implicit none
        integer(kind=4) :: ibuff(n_tree_files),i,j,ifile,unit
        do ifile = 1,n_tree_files
           unit = dumptree_unit + ifile
           if (ibuff(ifile) > 0) then 
              write(unit) ((bufferIDs(i,j,ifile),i=1,nIDsPerHalo),j=1,ibuff(ifile))
              write(unit) ((bufferIndexes(i,j,ifile),i=1,nIndexesPerHalo),j=1,ibuff(ifile))
           end if
        end do
        ibuff         = 0
        bufferIDs     = -1
        bufferIndexes = -1
        return
      end subroutine dump_buffer

      subroutine open_props_file(ifile,nhperts)
        implicit none
        integer(kind=4) :: ifile,unit,nhperts(nsteps)
        character(512)  :: file
        write(file,'(a,a,i3.3,a,i3.3)') trim(data_dir),'props_',nsteps,'.',ifile
        unit = dumptree_unit + ifile
        open(unit=unit,file=file,status="unknown",form="unformatted")
        write(unit) nsteps,nPropsPerHalo
        write(unit) nhperts(1:nsteps)
      end subroutine open_props_file

      subroutine dump_props_buffer(ibuff)
        implicit none
        integer(kind=4) :: ibuff(n_tree_files),i,j,ifile,unit
        do ifile = 1,n_tree_files
           unit = dumptree_unit + ifile
           if (ibuff(ifile) > 0) then 
              write(unit) ((propsBuffer(i,j,ifile),i=1,nPropsPerHalo),j=1,ibuff(ifile))
           end if
        end do
        ibuff       = 0
        propsBuffer = 0.0
        return
      end subroutine dump_props_buffer

  end subroutine buffered_outputs

!*****************************************************************************************************************
  
  subroutine init_tsTree

    implicit none 
    
    integer(kind=4) :: ts,n,ih,nd,j

    ! allocate and fill arrays with existing information from tree(:,:)
    allocate(ttree(nsteps))
    do ts = 1,nsteps
       n = nb_of_halos(ts) + nb_of_subhalos(ts)
       allocate(ttree(ts)%stree(n))
       allocate(ttree(ts)%jtree(n))
       allocate(ttree(ts)%itree(n))
       do ih = 1,n
          ! jTree
          ttree(ts)%jtree(ih)%BushID            = tree(ts,ih)%BushID
          ttree(ts)%jtree(ih)%TreeID            = -1
          ttree(ts)%jtree(ih)%HaloID            = -1
          ttree(ts)%jtree(ih)%FirstProgenitorID = -1
          ttree(ts)%jtree(ih)%firstProg         = -1
          ttree(ts)%jtree(ih)%NextProgenitorID  = -1
          ttree(ts)%jtree(ih)%nextProg          = -1 
          ttree(ts)%jtree(ih)%DescendentID      = -1
          ttree(ts)%jtree(ih)%descendent        = -1 
          ttree(ts)%jtree(ih)%HostHaloID        = -1
          ttree(ts)%jtree(ih)%HostSubID         = -1
          ttree(ts)%jtree(ih)%NextSubID         = -1
          ttree(ts)%jtree(ih)%haloNum           = ih
          ttree(ts)%jtree(ih)%haloTimestep      = ts
          ! sTree
          ttree(ts)%stree(ih)%halonum = ih
          ttree(ts)%stree(ih)%son     = tree(ts,ih)%my_sons%main_son
          nd = tree(ts,ih)%my_dads%nb_dads
          ttree(ts)%stree(ih)%ndads   = nd
          allocate(ttree(ts)%stree(ih)%dads(nd))
          do j = 1,nd 
             ttree(ts)%stree(ih)%dads(j) = tree(ts,ih)%my_dads%list_dads(j)
          end do
          ! iTree
          ttree(ts)%itree(ih)%me         = -1
          ttree(ts)%itree(ih)%descendent = -1
          ttree(ts)%itree(ih)%ndads      = nd
          ttree(ts)%itree(ih)%firstProg  = -1 
          ttree(ts)%itree(ih)%nextProg   = -1
          ttree(ts)%itree(ih)%hosthalo   = -1
          ttree(ts)%itree(ih)%hostsub    = -1
          ttree(ts)%itree(ih)%nextsub    = -1
       end do
    end do

    return

  end subroutine init_tsTree

!*****************************************************************************************************************

  subroutine define_IDs

    ! define all IDs in jtree

    implicit none

    integer(kind=4)           :: ts,ih,nh,hostsub,hosthalo,nextsub
    integer(kind=8)           :: TreeID
    integer(kind=8),parameter :: largeInt = 1000000  ! max number of halos per tree
    integer(kind=8)           :: offset

    ! start from z = 0 and recurse down each tree
    ! assign a unique TreeID to each tree, and define depth-first-ordered HaloIDs within each tree
    TreeID = 0
    do ts = nsteps,1,-1 
       nh = nb_of_halos(ts) + nb_of_subhalos(ts)
       do ih = 1,nh
          if (ttree(ts)%jtree(ih)%TreeID == -1) then  ! if halo has not been scanned yet get its tree
             progcnt = 0                ! initialize haloid counted for current (simple) tree
             call get_my_halo_stree(treeID,ts,ih)
             ! jeje 
             ttree(ts)%jtree(ih)%LastProgenitorID = progcnt
             ! end jeje 
             TreeID = TreeID + 1        ! whole tree has been done -> increment treeID 
          end if
       end do
    end do
    
    ! add offsets so that each ID is unique within the whole simulation
    ! NB: don't add offset to -1 values, as -1 means something ...
    do ts = 1,nsteps
       nh = nb_of_halos(ts) + nb_of_subhalos(ts)
       do ih = 1,nh
          offset = ttree(ts)%jtree(ih)%TreeID * largeInt
          if (ttree(ts)%jtree(ih)%HaloID /= -1) then 
             ttree(ts)%jtree(ih)%HaloID = ttree(ts)%jtree(ih)%HaloID + offset
          else
             write(errunit,*) '> ERROR: HaloID undefined'
             stop
          end if
          if (ttree(ts)%jtree(ih)%FirstProgenitorID /= -1) & 
               ttree(ts)%jtree(ih)%FirstProgenitorID = ttree(ts)%jtree(ih)%FirstProgenitorID + offset
          if (ttree(ts)%jtree(ih)%LastProgenitorID /= -1) & 
               ttree(ts)%jtree(ih)%LastProgenitorID = ttree(ts)%jtree(ih)%LastProgenitorID + offset
          if (ttree(ts)%jtree(ih)%NextProgenitorID /= -1) & 
               ttree(ts)%jtree(ih)%NextProgenitorID = ttree(ts)%jtree(ih)%NextProgenitorID + offset
          if (ttree(ts)%jtree(ih)%DescendentID /= -1) & 
               ttree(ts)%jtree(ih)%DescendentID = ttree(ts)%jtree(ih)%DescendentID + offset
       end do
    end do

    ! now that all HaloIDs are defined, define horizontal connections 
    do ts = 1,nsteps
       nh = nb_of_halos(ts) + nb_of_subhalos(ts)
       do ih = 1,nh
          hostsub  = tree(ts,ih)%hostsub
          if (hostsub > 0) then 
             ttree(ts)%jtree(ih)%HostSubID  = ttree(ts)%jtree(hostsub)%HaloID
          else if (hostsub == 0) then  ! ih is main halo
             ttree(ts)%jtree(ih)%HostSubID  = ttree(ts)%jtree(ih)%HaloID
          end if
          hosthalo = tree(ts,ih)%hosthalo
          ttree(ts)%jtree(ih)%HostHaloID = ttree(ts)%jtree(hosthalo)%HaloID
          nextsub  = tree(ts,ih)%nextsub
          if (nextsub > 0) ttree(ts)%jtree(ih)%NextSubID  = ttree(ts)%jtree(nextsub)%HaloID
       end do
    end do

    return

  end subroutine define_IDs

!*****************************************************************************************************************
  
  recursive subroutine get_my_halo_stree(treeID,ts,ihalo)

    ! scan the tree in a depth-first fashion to assign relative HaloIDs
    ! Define FirstProg, NextProg, and LastProg on the fly

    implicit none
    
    integer(kind=8) :: treeID
    integer(kind=4) :: ts,ihalo
    integer(kind=4) :: ip,phno,tsminus1
    integer(kind=4) :: prevprog

    prevprog = -1
    ttree(ts)%jtree(ihalo)%HaloID = progcnt
    ttree(ts)%jtree(ihalo)%TreeID = treeID
    do ip = 1, ttree(ts)%stree(ihalo)%ndads   ! loop on progenitors
       tsminus1 = ts -1
       progcnt  = progcnt + 1
       phno     = ttree(ts)%stree(ihalo)%dads(ip)
       ttree(tsminus1)%jtree(phno)%DescendentID = ttree(ts)%jtree(ihalo)%HaloID
       ttree(tsminus1)%jtree(phno)%descendent   = ihalo
       if (ip == 1) then  ! flag first progenitor
          ttree(ts)%jtree(ihalo)%FirstProgenitorID = progcnt
          ttree(ts)%jtree(ihalo)%firstProg         = phno
          prevprog = phno
       else               ! define next progenitor
          ttree(tsminus1)%jtree(prevprog)%NextProgenitorID = progcnt
          ttree(tsminus1)%jtree(prevprog)%nextProg         = phno
          prevprog = phno
       end if
       call get_my_halo_stree(treeID,tsminus1,phno)   ! recurse 
       ttree(tsminus1)%jtree(phno)%LastProgenitorID = progcnt  
    end do


    return

  end subroutine get_my_halo_stree

!*****************************************************************************************************************
#endif
end module jeje


