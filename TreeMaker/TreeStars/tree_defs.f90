module tree_defs

  public

  !======================================================================
  ! useful types :
  !======================================================================

  type vector
     real(kind=8)            :: x,y,z
  end type vector

  type father
     real(kind=8),pointer    :: mass_fathers(:)  ! percentage of father mass that goes into son halo
     integer(kind=4),pointer :: list_fathers(:)
     integer(kind=4)         :: nb_fathers
  end type father        
                             
  type shape                 
     real(kind=8)            :: a,b,c
  end type shape             
                             
  type baryon                
     real(kind=8)            :: rvir,mvir,tvir,cvel,Reff
  end type baryon

  type son
     integer(kind=4),pointer :: list_sons(:)
     integer(kind=4)         :: nb_sons
     integer(kind=4)         :: main_son
  end type son               
                             
  type hprofile              
     real(kind=8)            :: rho_0,r_c
  end type hprofile          
                             
  type tree_struct
     integer(kind=4) :: host        ! actually "hostsub" form HaloMaker -> points to host sub-structure or has -1 value for main haloes
     integer(kind=4) :: nextsub     ! 
     type(father)    :: my_fathers
     type(son)       :: my_sons
     integer(kind=4) :: BushID
  end type tree_struct

  type halo                  
     type (baryon)           :: datas
     type (shape)            :: sh
     type (vector)           :: p
     type (vector)           :: v
     type (vector)           :: L
     type (hprofile)         :: halo_profile 
     integer(kind=4)         :: my_number
     integer(kind=4)         :: my_timestep
     integer(kind=4)         :: level, hosthalo, hostsub,nbsub,nextsub ! data for structure tree
#ifndef BIG_RUN              
     integer(kind=4)         :: ncont ! number of contaminating low-resolution particles  
#endif                       
     real(kind=8)            :: m
     real(kind=8)            :: r
     real(kind=8)            :: spin
     real(kind=8)            :: sigma
     real(kind=8)            :: sigma_bulge
     real(kind=8)            :: m_bulge
     real(kind=8)            :: ek,ep,et
     real(kind=8)            :: macc
     integer(kind=4)         :: nbin=100
     real(kind=8),dimension(1:100) :: rr,rho
  end type halo

  !======================================================================

  !======================================================================
  ! Definitions specific to input/output
  !======================================================================
  character(80)                   :: data_dir
  integer(kind=4)                 :: errunit = 0
  !======================================================================
  
  !======================================================================
  ! parameters necessary to the brick file analysis
  !======================================================================
  integer(kind=4)                 :: max_groups
  integer(kind=4)                 :: nsteps,nbodies
  !======================================================================

  !======================================================================
  ! Global variables 
  !======================================================================
  real(kind=8)                           :: massp
#ifndef BIG_RUN
  real(kind=8),allocatable               :: mass(:)    ! mass of each particle
  integer(kind=4)                        :: nlr        ! number of low-resolution particles
#endif

  real(kind=8),allocatable               :: MassToAccrete(:) ! mass that a particle still transports (0 when it has been accreted)

  integer(kind=4),allocatable            :: liste_parts(:),ex_liste_parts(:),first_part(:),linked_list(:)
  type(halo),allocatable,target          :: liste_halos(:) ! list of all halos at all timesteps
  integer(kind=4),allocatable            :: nb_of_halos(:), nb_of_subhalos(:)
  real(kind=8),allocatable               :: aexp(:),age_univ(:),omega_t(:)
  integer(kind=4)                        :: numero_step
  integer(kind=4)                        :: num_halo_final
  integer(kind=4)                        :: n_tree_files  ! approximative nb of tree files outputed (exact if = 1)
  type(tree_struct),allocatable,target   :: tree(:,:)
  !======================================================================

  integer(kind=4),parameter              :: dumptree_unit = 500


end module tree_defs


  




