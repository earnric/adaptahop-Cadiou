module tree_defs

  public

  !======================================================================
  ! useful types :
  !======================================================================

  type vector
     real(kind=4)            :: x,y,z
  end type vector
                             
  type shape                 
     real(kind=4)            :: a,b,c
  end type shape             
                             
  type baryon                
     real(kind=4)            :: rvir,mvir,tvir,cvel
  end type baryon

  type father
     real(kind=8),pointer    :: mass_fathers(:)  ! #ifndef SIMPL percentage of father mass that goes into son halo  #ifdfSIMPL father mass that goes into son halo
     integer(kind=4),pointer :: list_fathers(:)
     integer(kind=4)         :: nb_fathers
  end type father        
  
#ifdef SIMPL
  type dad
     real(kind=4),pointer    :: mass_dads(:)
     integer(kind=4),pointer :: list_dads(:)
     integer(kind=4)         :: nb_dads
  end type dad
#endif

  type son
     integer(kind=4),pointer :: list_sons(:)
     integer(kind=4)         :: nb_sons
     integer(kind=4)         :: main_son    ! main son of the halo in the simplified merger tree
  end type son               
                             
  type hprofile              
     real(kind=4)            :: rho_0,r_c
  end type hprofile          
                             
  type tree_struct
!!$#ifdef SIMPL
!!$     real(kind=4)    :: m           ! copy mass from liste_halos
!!$#endif
#ifndef SIMPL
     integer(kind=4) :: host        ! actually "hostsub" from HaloMaker -> points to host sub-structure or has -1 value for main haloes
     integer(kind=4) :: nextsub     ! linked list of subhalo starting from the host halo (level 1) to the last halo within it
#endif
#ifdef SIMPL
     integer(kind=4) :: level
     integer(kind=4) :: hosthalo
     integer(kind=4) :: hostsub
     integer(kind=4) :: nextsub
     integer(kind=4) :: frag
#endif
     type(father)    :: my_fathers
#ifdef SIMPL
     type(dad)       :: my_dads
#endif
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
#ifndef SIMPL
     integer(kind=4)         :: level, hosthalo, hostsub,nbsub,nextsub ! data for structure tree
#endif
#ifndef BIG_RUN              
     integer(kind=4)         :: ncont ! number of contaminating low-resolution particles  
#endif                       
     real(kind=4)            :: m
     real(kind=4)            :: r
     real(kind=4)            :: spin
     real(kind=4)            :: ek,ep,et
     real(kind=8)            :: macc
#ifdef CONTAM
     integer(kind=4)         :: contamination
#endif
  end type halo

  !======================================================================

  !======================================================================
  ! Definitions specific to input/output
  !======================================================================
  character(200)                   :: data_dir
  integer(kind=4)                 :: errunit = 6
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
  real(kind=4)                           :: massp
#ifndef BIG_RUN
  real(kind=4),allocatable               :: mass(:)    ! mass of each particle
  integer(kind=4)                        :: nlr        ! number of low-resolution particles
#endif

  real(kind=8),allocatable               :: MassToAccrete(:) ! mass that a particle still transports (0 when it has been accreted)

  integer(kind=4),allocatable            :: liste_parts(:),ex_liste_parts(:),first_part(:),linked_list(:)
  type(halo),allocatable,target          :: liste_halos(:) ! list of all halos at all timesteps
  integer(kind=4),allocatable            :: nb_of_halos(:), nb_of_subhalos(:)
  real(kind=4),allocatable               :: aexp(:),age_univ(:),omega_t(:)
  integer(kind=4)                        :: numero_step,st_do,nsteps_do
  integer(kind=4)                        :: num_halo_final
  integer(kind=4)                        :: n_tree_files  ! approximative nb of tree files outputed (exact if = 1)
  type(tree_struct),allocatable,target   :: tree(:,:)
  !======================================================================

  integer(kind=4),parameter              :: dumptree_unit = 500


end module tree_defs


  




