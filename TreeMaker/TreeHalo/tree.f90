program tree

  implicit none 
  
  type tree_indexes
     integer(kind=4) :: me
     integer(kind=4) :: descendent
     integer(kind=4) :: ndads
     integer(kind=4) :: firstProg
     integer(kind=4) :: nextProg
     integer(kind=4) :: hosthalo
     integer(kind=4) :: hostsub
     integer(kind=4) :: nextsub
     integer(kind=4) :: level
  end type tree_indexes

  ! container
  type timestep_type 
     type(tree_indexes),allocatable :: tree(:)
  end type timestep_type
  type(timestep_type),allocatable :: tsno(:)

  integer(kind=4)             :: nsteps
  integer(kind=4),allocatable :: nb_of_halos(:)

!!$  character(512)              :: treefile = '/data2/blaizot/512-100h-1Mpc-W3_new/TreeMaker/tree_file_092.115'
!!$  character(512)              :: propfile = '/data2/blaizot/512-100h-1Mpc-W3_new/TreeMaker/props_092.115'
!!$  character(512)              :: treefile = '/data2/blaizot/1024-100h-1Mpc-W3/TreeMaker/tree_file_095.083'
!!$  character(512)              :: propfile = '/data2/blaizot/1024-100h-1Mpc-W3/TreeMaker/props_095.083'
  integer(kind=4),parameter   :: errunit = 6

  integer(kind=4),parameter   :: treeunit = 12
  integer(kind=4)             :: nIDsPerHalo,nIndexesPerHalo,nPropsPerHalo
  integer(kind=4),allocatable :: indexes(:,:)
  real(kind=4),allocatable    :: props(:,:)


  integer(kind=4) :: dad,ts,i,j
  integer(kind=4),parameter :: maxNprog = 20
  integer(kind=4) :: counts(maxNProg)
  character(100)  :: fmt 
  
  integer(kind=4) :: ifile 
  character(512)  :: fichier,fichier2
  
  ! loop on tree files 
!!$  open(unit=133,file='halo_props_512.dat',status='unknown',form='formatted')
!!$  open(unit=133,file='halo_props_1024.dat',status='unknown',form='formatted')
  do ifile = 1,120
     write(fichier,'(a,i3.3)') '/data2/blaizot/512-100h-1Mpc-W3_new/TreeMaker/props_092.',ifile
     write(fichier2,'(a,i3.3)') '/data2/blaizot/512-100h-1Mpc-W3_new/TreeMaker/tree_file_092.',ifile
!!$     write(fichier,'(a,i3.3)') '/data2/blaizot/1024-100h-1Mpc-W3/TreeMaker/props_094.',ifile
!!$     write(fichier2,'(a,i3.3)') '/data2/blaizot/1024-100h-1Mpc-W3/TreeMaker/tree_file_094.',ifile



     ! READ HIERARCHICAL INDEXES AND ALLOCATE TSNO
     open(unit=treeunit,file=fichier2,status='old',form='unformatted')
     read(treeunit) nsteps,nIDsPerHalo,nIndexesPerHalo
     !write(errunit,*) '> NSTEPS = ',nsteps
     allocate(nb_of_halos(nsteps),tsno(nsteps))
     read(treeunit) nb_of_halos(1:nsteps)
     do ts = 1,nsteps
        allocate(tsno(ts)%tree(nb_of_halos(ts)))
     end do
     allocate(indexes(nIndexesPerHalo,maxval(nb_of_halos)))
     do ts = 1,nsteps
        if (nb_of_halos(ts) == 0) cycle
        read(treeunit) ! skip IDs... 
        read(treeunit) ((indexes(i,j),i=1,nIndexesPerHalo),j=1,nb_of_halos(ts))
        do i = 1,nb_of_halos(ts) 
           tsno(ts)%tree(i)%me         = indexes(1,i)
           tsno(ts)%tree(i)%descendent = indexes(2,i)
           tsno(ts)%tree(i)%ndads      = indexes(3,i)
           tsno(ts)%tree(i)%firstProg  = indexes(4,i)
           tsno(ts)%tree(i)%nextProg   = indexes(5,i)
           tsno(ts)%tree(i)%hosthalo   = indexes(6,i)
           tsno(ts)%tree(i)%hostsub    = indexes(7,i)
           tsno(ts)%tree(i)%nextsub    = indexes(8,i)
           tsno(ts)%tree(i)%level      = indexes(9,i)
        end do
     end do
     deallocate(indexes,nb_of_halos)
     close(treeunit)

     ! READ HALO PROPERTIES
!!$  open(unit=treeunit,file=propfile,status='old',form='unformatted')
!!$  read(treeunit) nsteps,nPropsPerHalo
!!$  read(treeunit) ! skip nb_of_halos
!!$  allocate(props(nPropsPerHalo,maxval(nb_of_halos)))
!!$  do ts = 1,nsteps
!!$     if (nb_of_halos(ts) == 0) cycle
!!$     read(treeunit) ((props(i,j),i=1,nPropsPerHalo),j=1,nb_of_halos(ts))
!!$  end do
!!$  deallocate(props)
!!$  close(treeunit)


     open(unit=treeunit,file=fichier,status='old',form='unformatted')
     read(treeunit) nsteps,nPropsPerHalo
     allocate(nb_of_halos(nsteps))
     read(treeunit) nb_of_halos(1:nsteps)
     allocate(props(nPropsPerHalo,maxval(nb_of_halos)))
     do ts = 1,nsteps
        if (nb_of_halos(ts) == 0) cycle
        read(treeunit) ((props(i,j),i=1,nPropsPerHalo),j=1,nb_of_halos(ts))
     end do
     do j = 1,nb_of_halos(nsteps)
        write(133,'(3(e14.6,1x),i3)') props(10,j),props(11,j),props(13,j),tsno(nsteps)%tree(j)%level
     end do

     deallocate(props,nb_of_halos)
     close(treeunit)

!!$  propsbuffer(:,ibuff(ifile),ifile) = (/liste_halos(ih)%p%x,liste_halos(ih)%p%y, liste_halos(ih)%p%z, &
!!$       & liste_halos(ih)%v%x,liste_halos(ih)%v%y,liste_halos(ih)%v%z, & 
!!$       & liste_halos(ih)%m,liste_halos(ih)%r,liste_halos(ih)%spin,    & 
!!$       & liste_halos(ih)%datas%rvir,liste_halos(ih)%datas%mvir,liste_halos(ih)%datas%tvir, & 
!!$       & liste_halos(ih)%datas%cvel,real(liste_halos(ih)%macc,4)/)
  
     do ts = 1,nsteps 
        deallocate(tsno(ts)%tree)
     end do
     deallocate(tsno)
  end do
  close(133)

!!$
!!$  ! Check that links make sense 
!!$  do ts = 1,nsteps
!!$     do i = 1,nb_of_halos(ts)
!!$        ! am i the son of my dads ? 
!!$        dad = tsno(ts)%tree(i)%firstProg
!!$        do while (dad > 0)
!!$           if (tsno(ts-1)%tree(dad)%descendent /= i) then 
!!$              write(errunit,*) 'There is something rotten in TreeMaker ... '
!!$           end if
!!$           dad = tsno(ts-1)%tree(dad)%nextProg
!!$        end do
!!$     end do
!!$  end do
!!$  
!!$  ! count halos
!!$  write(fmt,'(a,i,a)') '(',maxNProg + 1,'(i8,1x))'
!!$  do ts = 1,nsteps
!!$     counts = 0
!!$     if (ts == 36) cycle
!!$     do i = 1,nb_of_halos(ts)
!!$        if (tsno(ts)%tree(i)%ndads < maxNProg) then 
!!$           counts(tsno(ts)%tree(i)%ndads+1) = counts(tsno(ts)%tree(i)%ndads+1) + 1
!!$        end if
!!$     end do
!!$     write(1024,fmt) ts,counts(:)
!!$  end do

end program tree



