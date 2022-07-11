program read_outputs

  implicit none 

  ! global parameters ... 
  character(80)                   :: data_dir
  integer(kind=4)                 :: errunit = 6
 
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
  integer(kind=4)   :: nHalosInBuffer
  integer(kind=4)   :: nIDsPerHalo      ! nb of fields in jeje_treeID type
  integer(kind=8),allocatable :: bufferIDs(:,:,:)

  ! indexes which are usefull to navigate between progs/descendents and halos/sub-halos in GALICS
  ! These indexes are valid within one treeFile (and hence defined just before output)
  type jejeIndexes
     integer(kind=4) :: me
     integer(kind=4) :: descendent
     integer(kind=4) :: firstProg
     integer(kind=4) :: nextProg
     integer(kind=4) :: hosthalo
     integer(kind=4) :: hostsub
     integer(kind=4) :: nextsub
     integer(kind=4) :: level
  end type jejeIndexes
  integer(kind=4)   :: nIndexesPerHalo  ! nb of fields in jejeIndexes type
  integer(kind=4),allocatable :: bufferIndexes(:,:,:)

  integer(kind=4)   :: nPropsPerHalo ! pos(3),vel(3),mfof,mvir,tvir,rvir,cvel,macc,rfof,spin
  real(kind=4),allocatable    :: propsbuffer(:,:,:)

  type props_type
     real(kind=4) :: pos(3)
     real(kind=4) :: vel(3)
     real(kind=4) :: mfof
     real(kind=4) :: mvir
     real(kind=4) :: tvir
     real(kind=4) :: rvir
     real(kind=4) :: cvel
     real(kind=4) :: macc
     real(kind=4) :: rfof
     real(kind=4) :: spin
  end type props_type

  ! global timestep container 
  type tsno_type
     type(jejeIndexes),pointer :: itree(:)
     type(jeje_treeID),pointer :: jtree(:)
     type(props_type),pointer  :: props(:)
  end type tsno_type
  type(tsno_type),allocatable  :: tsno(:)

  integer(kind=4) :: nsteps, i_tree_file
  character(512)  :: file  

  
  


contains

  subroutine read_a_tree_file(treefile) 
    
    ! includes allocation of tree structure... 

    implicit none

    integer(kind=4) :: nIDsPerHAlo,nIndexesPerHalo,nPropsPerHalo,nHalosInBuffer,nsteps,ok
    integer(kind=4) :: nbuff,ts,i,j,n1,n2
    integer(kind=8),allocatable :: IDbuffer(:,:)
    integer(kind=4),allocatable :: IndexBuffer(:,:)
    real(kind=4),allocatable    :: propsBuffer(:,:)

    open(unit=12,file=treefile,form='unformatted',status='old')
    read(12) nsteps,nIDsPerHalo,nIndexesPerHalo,nHalosInBuffer,nb_of_halos(1:nsteps)

    allocate(tsno(nsteps))
    do ts = 1,nsteps
       allocate(tsno(ts)%itree(nb_of_halos(ts)))
       allocate(tsno(ts)%jtree(nb_of_halos(ts)))
       allocate(tsno(ts)%props(nb_of_halos(ts)))
    end do
    
    ! allocate buffers for reading
    allocate(IDbuffer(nIDsPerHalo,nHalosInBuffer))
    allocate(IndexBuffer(nIndexesPerHalo,nHalosInBuffer))

    ok = 1
    n1 = 1
    n2 = 1
    do while (ok == 1) 
       read(12,iostat=ok) nbuff
       read(12) ((IDbuffer(i,j),i=1,nIDsPerHalo),j=1,nbuff)
       read(12) ((Indexbuffer(i,j),i=1,nIndexesPerHalo),j=1,nbuff)
       
    end do
    close(12)

  end subroutine read_a_tree_file

!*****************************************************************************************************************

  subroutine init_tsTree

    implicit none 
    
    integer(kind=4) :: ts,n,ih,nd,j

    ! allocate and fill arrays with existing information from tree(:,:)
    allocate(tsno(nsteps))
    do ts = 1,nsteps
       n = nb_of_halos(ts) + nb_of_subhalos(ts)
       allocate(tsno(ts)%stree(n))
       allocate(tsno(ts)%jtree(n))
       allocate(tsno(ts)%itree(n))
       do ih = 1,n
          ! jTree
          tsno(ts)%jtree(ih)%BushID            = tree(ts,ih)%BushID
          tsno(ts)%jtree(ih)%TreeID            = -1
          tsno(ts)%jtree(ih)%HaloID            = -1
          tsno(ts)%jtree(ih)%FirstProgenitorID = -1
          tsno(ts)%jtree(ih)%firstProg         = -1
          tsno(ts)%jtree(ih)%NextProgenitorID  = -1
          tsno(ts)%jtree(ih)%nextProg          = -1 
          tsno(ts)%jtree(ih)%DescendentID      = -1
          tsno(ts)%jtree(ih)%descendent        = -1 
          tsno(ts)%jtree(ih)%HostHaloID        = -1
          tsno(ts)%jtree(ih)%HostSubID         = -1
          tsno(ts)%jtree(ih)%NextSubID         = -1
          tsno(ts)%jtree(ih)%haloNum           = ih
          tsno(ts)%jtree(ih)%haloTimestep      = ts
          ! props
          tsno(ts)%props(ih) = 
          
          ! iTree
          tsno(ts)%itree(ih)%me         = -1
          tsno(ts)%itree(ih)%descendent = -1
          tsno(ts)%itree(ih)%firstProg  = -1 
          tsno(ts)%itree(ih)%nextProg   = -1
          tsno(ts)%itree(ih)%hosthalo   = -1
          tsno(ts)%itree(ih)%hostsub    = -1
          tsno(ts)%itree(ih)%nextsub    = -1
       end do
    end do

    return

  end subroutine init_tsTree

!*****************************************************************************************************************



end program read_outputs
  



