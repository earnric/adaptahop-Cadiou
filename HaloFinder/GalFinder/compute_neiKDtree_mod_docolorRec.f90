!=======================================================================
!                           COMPUTE_NEIKDTREE
!=======================================================================
! Author : S. Colombi
!          Institut d'Astrophysique de Paris
!          98 bis bd Arago, F-75014, Paris, France 
!          colombi@iap.fr
!
! This program has multiple usage and can do three things
! (1) Read an input particle distribution file and compute the mean
!     square distance between each particle and its nearest neighbourgs
! (2) Read an input particle distribution file and compute the SPH
!     density associated to each particle + the list of its nearest
!     neighbourgs
! (3) Read an input particle distribution and a neighbourgs file which
!     is an output of step 2), and output the tree of the structures in
!     structures.
!
! If steps (1) and (2) are probably correct, step (3) is still under
! test and construction. Step (1) was rather extensively tested, while
! step (2) needs further tests although I am quite confident it should
! be correct.
!
!=======================================================================
!                           COMPUTE_NEIKDTREE_MOD
!=======================================================================
! Modication of compute_neiKDtree to fit inside GalaxyMaker
! Modification : D. Tweed
!
! PARAMETERS IN THE CONFIG FILE 
! +++++++++++++++++++++++++++++
! 
! A example of config file is given in compute_neiKDtree.config
!
!
! verbose     : normal verbose mode (.true. or .false.)
! megaverbose : full verbose mode (.true. or .false.)
! filein      : name of the input particles/velocities file ('myfilename')
! Ntype       : format of the file (integer number)
!               0 : PM simple precision unformatted
!               1 : GADGET simple precision unformatted, cosmological simulation
!                   of dark matter only with periodic boundaries
!               2 : RAMSES unformatted (particles) simulation
!                   MPI multiple output
!               3 : CONE unformatted 
!               4 : Ninin treecode unformatted
! remove_degenerate : if this set to .true., check for particles at the same 
!               position and add a random displacement at approximately 
!               floating representation accuracy in order to avoid infinite
!               number of KD tree cells creation. This parameter is global since
!               if a random displacement has been applied once for computing
!               SPH density, the same one must be applied again for subsequent 
!               calculations in order to be self-consistent.
!               Setting remove_degenerate=.false. will of course speep-up
!               the program at the cost of a risk of crash due to infinite
!               KD tree cell creation. This is however safe in most cases since
!               KD tree creation can now deal with particles at the same position.
!               It is therefore recommended to set remove_degenerate to .false.
!               first before trying .true.
! action      : 'distances' for computing mean square distance between each
!                           particle and its nearest neighbourgs
!               'neighbors' for computing the list of the nearest neighbourgs of
!                           each particle and its SPH density
!               'adaptahop' for computing the tree of structures and substructures
!                           in the simulation 
!               'posttreat' for computing the physical properties of dynamically
!                           selected haloes and subhaloes
!
! if (Ntype=1) : GADGET format
! --------------
!
! nmpigadget  : number of CPU used to perform the simulation (in order to read the 
!               appropriate number of files). If there is no multiple file 
!               extension (only one file), set nmpigadget=-1.
!
! if (Ntype=2) : RAMSES unformatted (dark matter particles)
! --------------
!
! mtot        : total mass in the full simulation box, of the population considered,
!               in internal RAMSES units.
!               mtot=-1 works if full simulation box is analysed.
! Ltot        : size of the full simulation box, in internal RAMSES units.
!               Ltot=-1 gives default value, valid for cosmological simulations
!               only.
! 
! if (Ntype=3) : CONE unformatted (dark matter particles) 
! ------------------
!
! boxsize     : comoving size of the simulation box in Mpc
!
! if (Ntype=4) : Ninin treecode unformatted 
! --------------
!
! boxsize2    : comoving size of the simulation box, in Mpc
! hubble      : value of H0/100 in km/s/Mpc
! omega0      : value of cosmological density parameter
! omegaL      : value of cosmological constant
! aexp_max    : value of expansion factor at present time
!
! if (action='distances') :
! -------------------------
!
! nvoisdis    : number of neighbors considered for each particle
! filedis     : name of the output file for the mean square distances 
!               ('mydistancename')
!               (this is a binary unformatted or formatted file which will be 
!               described more in detail in the near future) 
! formatted   : .true. if the output file is in ASCII, .false. if the output
!               file is in binary unformatted double float.
! ncpu        : number of virtual threads for a parallel calculation with openMP
!               This integer number should be a multiple of the real number of 
!               threads used during the run. The larger it will be, the better
!               the load balancing will be. 
! velocities  : .true. or .false. : if set to .true., include velocities 
!               treatment in the calculation and the output of the results
!
! if (action='neighbors') :
! -------------------------
!
! nvoisnei    : number of neighbors considered for computing the SPH density
!               (integer number). Typically this number should vary between
!               20 and 100 (template value 64).
! nhop        : number of stored nearest neighbourgs in the output file
!              (integer number).
!              nhop must smaller smaller or equal to nvoisnei. Typically nhop
!              should be of order 10-30 (template value 16).
! fileneinei  : name of the output file for the SPH density and the list of
!              nhop nearest neighbors for each particle ('myfileneighbors')
!              (this is a binary unformatted file which will be described more 
!              in detail in the near future)
!
! if (action='adaptahop') :
! -------------------------
!
! rho_threshold : density threshold. Particles with SPH density below this
!              density threshold are not selected. This thresholding will
!              define a number of connected regions, each of which corresponding
!              to a structure. For each of these structures, we aim to 
!              build a tree of substructures. Typically, one is interested
!              in finding the substructures of highly nonlinear objects. A value
!              often taken for rho_threshold is 80, which corresponds roughly
!              to friend-of-friend algorithm parameter b=2.
! nmembthresh: threshold on the number of particles that a structure or a
!              a substructure (above some density threshold rhot) must contain 
!              for being considered as significant.
!              The choice of nmembthresh is related to effects of N-body relaxation
!              (treecode) or over-softening of the forces (AMR or PM codes) or
!              SPH smoothing (SPH smoothing over a number N of particles tends
!              to ``erase'' structures with less than N particles).
!              Basically, due to either of these effects, structures or 
!              substructures with a number of particles smaller than nmembthresh
!              are not considered. A typical choice of nmembthresh is 64.
! fudgepsilon: This parameter can be seen as an Eulerian version of the thresholding
!              controlled by nmembthresh. It defines the size of the smallest structures
!              that can exist physically and is related to force softening.
!              Basically, if epsilon is the softening parameter, substructures with
!              typical radius smaller than epsilon are just coincidences and should
!              not be considered. The typical radius of a structure is given by the
!              mean square distance of particles belonging to it with respect to its
!              center of gravity. 
!              The criterion for selecting a structure is thus
!              
!              radius_{substructure} > epsilon,
!             
!              where epsilon=fudgepsilon*L/Npart^{1/3}. Fudgepsilon is thus expressed
!              in units of mean interparticle separation.
!              - For a treecode, fudgepsilon is typically of the order of 1/20;
!              - For a PM code, fudgepsilon is typically of the order of N_g/Npart^{1/3}
!              where N_g is the resolution of the grid used for the force calculations,
!              i.e. of the order of 1 or 0.5.
!              - For a quasi-Lagrangian code such as RAMSES, putting constrains from fudgespilon 
!              is not really useful. Basically it corresponds to N_g/Npart^{1/3} where N_g
!              would be the equivalent of the PM grid probed by the smallest AMR cells.              
! alpha      : a criterion can be applied to decide weither a substructure
!              is selected or not. Here the choice of alpha dictates by how
!              much the maximum local density in this substructure should be larger
!              than the local average density in this substructure:
!
!              rho_max_{substructure} >= alpha*< rho >_{substructure}.
!
!              Basically, the choice of alpha dictates the ``peakyness'' of a 
!              substructure. For instance, a substructure might contain itself
!              5 local maxima. For this substructure to be dynamically significant
!              we want the largest of the local maxima to be at least alpha
!              times the mean local density, i.e. the average density of all the
!              particles within the substructure. Dynamically bounded substructures
!              are expected to be very peaky. A typical choice of alpha would be
!              alpha=0 or a few unites, e.g. alpha=4.
! fudge      : a criterion can be applied to decide wither a substructure is 
!              statistically significant in terms of Poisson noise. Indeed,
!              even if Poisson noise is considerably reduced by SPH smoothing
!              it is still present to some extent, and some substructures might
!              be simply due to local Poisson fluctuations. 
!              Here the choice of fudge dictates by how many ``sigma's'' 
!              a structure must be compared to a random Poisson fluctuation.
!              The criterion applied is
!
!              < rho >_{substructure} > rhot*[1+fudge/sqrt(N)]
!
!              where N is the number of particles in the substructure and rhot
!              the density threshold corresponding to this substructure (in 
!              other worlds, rhot is the minimum density of particles whithin
!              the substructure). 
!              This criterion can be understood as follows : if a given substructure
!              contains N particles, the corresponding random poisson fluctuations are
!              of the order of sqrt(N). So typically the uncertainty on the density
!              estimate of a structure containing N particles is of the order of
!              sqrt(N) (this of course neglects the fact that SPH softening reduces
!              considerably the effects of poisson fluctuations). For the
!              fluctuation associated to the substructure to be significant, we
!              impose that (< rho >_{substructure} - rhot)/rhot > fudge*sigma with
!              sigma=1/sqrt(N)=sqrt(< (M-<M>)^2 >)/<M> where M is a Poisson process
!              of average N.
!              A typical value is fudge=4.
!
!              IMPORTANT NOTE : we should always have fudge > 0. Small values of
!              fudge will slow down the program, while large values of fudge will
!              make the hierarchical decomposition in terms of substructures less
!              accurate. 
!              The reason for that is that despite the fact we know all the saddle points
!              connecting all the substructures, it was impossible to figure out a simple
!              way of sorting these saddle points in order to construct automatically
!              the tree of structures and substructures. As a result, we operate 
!              iteratively by increasing the local threshold density as follows:
!              [rhot NEXT] = [rhot OLD]*[1+fudge/sqrt(N)] where N is the number of particles
!              in the substructure. Then we see weither saddle points within this
!              substructure are below the new value of rhot: if this happens, it
!              means that at the density level [rhot NEXT], the substructure is composed
!              of disconnected subsubstructures and the threshold value of connection
!              between these subsubstructures is approximated by [rhot NEXT] (the real value
!              should be between [rhot NEXT] and [rhot OLD]).
! filenode   : output file for the tree of structures and substructures. See below
!              for explanation of its format
! simu_unitsnei : .true. or .false : to specifie if the nodes positions and radii in the
!              file filenode are specified in the same units as in the simulation input
!              file or in Mpc
! filepartnodenei : ouput file for particule node number. To each particle, a
!              integer is associated, which is the node number the deepest possible in the tree
!              given in filenode. This way, at any level in the tree, it will be possible to
!              find recursively all the particles belonging to the tree. If the node number 
!              is zero, it means the SPH density of particle is below the threshold rho_threshold 
!              and is therefore associated to no structure. 
! formattedpartnodenei : .true. for having filepartnodenei in ascii format, .false. for
!              having filepartnodenei in binary format. The format is in both
!              case, the number of particles npart in the first line, followed by id(:), 
!              where id(:) is an array of npart integers.
! fileneihop : input file with SPH density and list of nearest neighbors of each
!              particle, obtained by running comptu_neiKDtree with action='neighbors'.
!
! if (action='posttreat') :
! -------------------------
! This option is still in development and test phase
!
!=======================================================================
!
! HISTORY
! +++++++
!
! 02/05/02/SC/IAP : first operational clean version with namelist and
!                   three options for the calculations.
!                   + action='distances' tested by comparison with the
!                   results of a different algorithm based on link list
!                   approach instead of KDtree (compute_neighbourgs2.f)
!                   The agreement between the 2 softwares is good within
!                   float calculation accuracy (5.10^{-6} relative 
!                   maximum differences). The differences can be 
!                   explained by a slightly different approach used
!                   in the 2 programs to do floating calculations.
!                   + action='neighbors' is not fully tested but should be
!                   okay. 
!                   + action='adaptahop' has been seen to work fine on
!                   a few data sets, but is not by any mean extensively
!                   tested. (No comments on the parameters of the namelist
!                   to avoid the users using this option)
!                   --> Version 0.0
! 04/05/02/SC/IAP : Full comments on the namelist parameters corresponding
!                   to action='adaptahop'. A few comments are added to
!                   the program. A new input particle file format is 
!                   added (RAMSES). 
!                   --> Version 0.1
! 04/17/02/SC&RT/IAP : Clean up the RAMSES format. Speed up the SPH 
!                   smoothing procedure
!                   --> Version 0.2
! 12/02/02/SC/IAP : Parallelize action='distances' and action='neighbors'
!                   Add RAMSES dark matter format with periodic boundaries
!                   (Ntype=-2), add CONE format (Ntype=3)
!                   --> Version 0.3
! 27/05/03/SC&RT/IAP&STR : Add RAMSES MPI dark matter + stars format with
!                   periodic boundaries (still under work, Ntype=-20) and
!                   a subroutine to remove degeneracy in particle positions
!                   Add GADGET multiple file MPI format.
! 04/06/03/SC/IAP : Dynamical selection of haloes and subhaloes :
!                   + new action='posttreat'
! 26/10/03/SC/IAP : add Ninin treecode format Ntype=4
! 16/10/06/SC/IAP : Change the treecode algorithm to be able to deal with
!                   2 or more particles at the same position. It is not
!                   yet completely safe in terms of the calculation of
!                   the SPH density (if all the particles are at the
!                   same position one gets 0/0).
!                   Add a few comments on file formats.
!                   --> Version 0.8
! 18/10/06/SC&RT/IAP : remove obsolete formats and improve RAMSES MPI
!                   format (now swapped from -20 to 2).
! 28/09/06/DT/CRAL: module compute_neiKDtree_mod.f90 compatibible with GalaxyMaker2.0
! 24/09/07/DT/CRAL: choice of 3 selection method using flag method
!                   HOP: Adaptahop is only used to detect haloes
!                   DPM: Subhaloes are selected through the density profile method, 
!                   accurate if studiing one step only
!                   MHM: Subhaloes are selected through the merger history method, 
!                   use this method if you need to obtain an accurate merger tree 
!                   containing subhaloes
!=======================================================================

module neiKDtree
  
  use halo_defs
  integer::ncall         !YDdebug
  integer::icolor_select !YDdebug
  contains

!=======================================================================
subroutine compute_adaptahop 
!=======================================================================

  implicit none
  
  call change_pos
  ! action neighbors
  call create_tree_structure
  call compute_mean_density_and_np
  ! action adaptahop
  call find_local_maxima
  call create_group_tree
  call change_pos_back
  ! check that we have halos
  call count_halos
  if(nb_of_halos.gt.0) then
     ! reinit halo and subhalo count
     nb_of_halos    = 0
     nb_of_subhalos = 0
     ! node structure tree -> halo structure tree: 3 methods avalable
     select case(method)
     case("HOP")
        call select_halos
     case("DPM")
        call select_with_DP_method
     case("MSM")
        call select_with_MS_method
     case("BHM")
        call select_with_BH_method
     case default
        write(errunit,*) '> could not recognise selection method:',method
     end select
  else
     deallocate(node)
  end if

  return

end subroutine compute_adaptahop

!=======================================================================
subroutine list_parameters
!=======================================================================

  implicit none

  write(errunit,*) '============================================================'
  write(errunit,*) 'Compute_neiKDtree was called with the following attributes :'
  write(errunit,*)     'verbose         :',verbose
  write(errunit,*)     'megaverbose     :',megaverbose
  write(errunit,*)     'boxsize2        :',boxsize2
  write(errunit,*)     'hubble          :',hubble
  write(errunit,*)     'omega0          :',omega0
  write(errunit,*)     'omegaL          :',omegaL
  write(errunit,*)     'aexp_max        :',aexp_max
  write(errunit,*)     'action          : neighbors'
  write(errunit,*)     'nvoisins        :',nvoisins
  write(errunit,*)     'nhop            :',nhop
  write(errunit,*)     'action          : adaptahop'
  write(errunit,*)     'rho_threshold   :',rho_threshold
  write(errunit,*)     'nmembthresh     :',nmembthresh
  write(errunit,*)     'alphap          :',alphap
  write(errunit,*)     'fudge           :',fudge
  write(errunit,*)     'fudgepsilon     :',fudgepsilon
  write(errunit,*)     'Selection method:',' ',method
  write(errunit,*) '============================================================'
end subroutine list_parameters

!=======================================================================
subroutine init_adaptahop
!=======================================================================

  implicit none

  omegaL   = omega_lambda_f
  omega0   = omega_f
  aexp_max = af
  hubble   = H_f*1e-2
  boxsize2 = Lf
  xlong    = boxsize2
  ylong    = xlong
  zlong    = xlong
  xlongs2  = xlong*0.5d0
  ylongs2  = ylong*0.5d0
  zlongs2  = zlong*0.5d0

  pos_renorm  =xlong
  pos_shift(1)=0.0d0
  pos_shift(2)=0.0d0
  pos_shift(3)=0.0d0
  write(*,*)'npart=',npart,dble(npart)
  mass_in_kg=(xlong*ylong*zlong/dble(npart))*mega_parsec**3 &
       &             *omega0*critical_density*hubble**2 
  Hub_pt = 100.*hubble * sqrt(omega0*(aexp_max/aexp)**3 + (1-omegaL-omega0)*(aexp_max/aexp)**2 + omegaL) 
  
  nmembthresh = nMembers  
  
  call list_parameters

  return

end subroutine init_adaptahop

!=======================================================================
subroutine change_pos
!=======================================================================

  implicit none

  npart    = nbodies
  epsilon  = fudgepsilon*xlong/dble(npart)**(1.d0/3.d0)
  pos      = pos * boxsize2

  return

end subroutine change_pos

!=======================================================================
subroutine change_pos_back
!=======================================================================

  implicit none
  
  pos = pos / boxsize2

end subroutine change_pos_back

!=======================================================================
subroutine count_halos
!=======================================================================

  implicit none
  integer(kind=4) :: inode
  
  nb_of_halos = 0
  do inode = 1,nnodes
     if(node(inode)%level.eq.1) nb_of_halos = nb_of_halos + 1
  end do

  return
  
end subroutine count_halos

!=======================================================================
subroutine select_halos
!=======================================================================

  implicit none
  integer(kind=4) :: inode, ihalo, ipar
  integer(kind=4) :: node_to_halo(nnodes)

  node_to_halo = -1

  ! counting number of halos
  nb_of_halos = 0
  do inode = 1, nnodes
     if(node(inode)%mother.le.0) nb_of_halos = nb_of_halos + 1
     ihalo = inode
     do while(node(ihalo)%mother .ne.0)
        if(node_to_halo(ihalo).le.0) then
           ihalo = node(ihalo)%mother
        else
           ihalo = node_to_halo(ihalo)
        end if
     end do
     node_to_halo(inode) = ihalo 
  end do
  write(errunit,*) 'number of nodes :',nnodes
  write(errunit,*) 'number of haloes:', nb_of_halos

  do ipar = 1, npart
     inode = liste_parts(ipar)
     if(inode.gt.0) then
        ihalo = node_to_halo(inode)
        if((ihalo.le.nb_of_halos).and.(ihalo.gt.0)) then
           liste_parts(ipar) = ihalo
        else
           stop 'ihalo .gt. nb_of_halos or nil'
        end if
     end if
  end do

  deallocate(node)
  
  return

end subroutine select_halos

!=======================================================================
subroutine select_with_DP_method
!=======================================================================

  implicit none
  integer(kind=4) :: ihalo,inode,isub,ip
  integer(kind=4) :: inodetmp, imother, istruct 
  integer(kind=4) :: nb_halos, nb_sub
  integer(kind=4) :: mostdenssub(nnodes), node_to_struct(nnodes)
  real(kind=8)    :: maxdenssub(nnodes)

  if(verbose) write(errunit,*) 'Using DP method'

  maxdenssub     = 0.0
  mostdenssub    = -1
  node_to_struct = -1

  ! detect for each node its most dens child
  do inode = 1, nnodes
     if(node(inode)%firstchild.gt.0) then
        inodetmp = node(inode)%firstchild
        do while(inodetmp.gt.0)
           if(node(inodetmp)%densmax.gt.maxdenssub(inode)) then
              maxdenssub(inode)     = node(inodetmp)%densmax
              mostdenssub(inode)    = inodetmp
           end if
           inodetmp = node(inodetmp)%sister
        end do
     end if
  end do
 
  ! Second run to write node_to_struct array and count the number of substructures
  nb_sub   = 0
  nb_halos = 0
  istruct  = 0
 
  do inode = 1, nnodes
     if(node(inode)%mother.le.0) then
        nb_halos = nb_halos + 1
        istruct  = istruct  + 1
        node_to_struct(inode) = istruct
     else
        if(node_to_struct(inode).le.0) then
           imother  = node(inode)%mother
           inodetmp = node(imother)%firstchild
           do while(inodetmp.gt.0)
              if(inodetmp.eq.mostdenssub(imother)) then
                 node_to_struct(inodetmp) = node_to_struct(imother)
              else
                 nb_sub = nb_sub + 1
                 istruct = istruct + 1
                 node_to_struct(inodetmp) = istruct
              end if
              inodetmp = node(inodetmp)%sister
           end do
        end if
     end if
  end do

  ! Check that the structure tree is ordered as it should be 
  do inode = 1, nnodes
     if(inode.le.nb_halos) then
        if(.not.((node(inode)%level.eq.1).and.(node_to_struct(inode).eq.inode))) stop 'error in node order' 
     else
        if((node(inode)%level.eq.1)) stop 'error shouldn''t be any halo here' 
     end if
  end do

  nstruct = nb_halos + nb_sub
  if(nstruct.ne.istruct) then
     write(errunit,*) 'Error in subroutines count,',nstruct, istruct
     stop
  end if
  
  write(errunit,*)
  write(errunit,*) '> number of nodes               :', nnodes
  write(errunit,*) '> number of halos               :', nb_halos
  write(errunit,*) '> number of substructures       :', nb_sub
  write(errunit,*) '> number of node removed        :', nnodes - nstruct
  write(errunit,*)

  if(verbose) write(errunit,*) '> Cleaning liste_parts'
  ! Cleaning liste_parts
  do ip = 1, nbodies
     inode  = liste_parts(ip)
     if(inode.gt.0) then
        if(node_to_struct(inode).le.0) then
           write(errunit,*) ip, inode,  node_to_struct(inode)
           stop 'error in node_to_struct'
        end if
        liste_parts(ip) = node_to_struct(inode)
     end if
  end do

  if(verbose) write(errunit,*) '> Creating new structure tree'
  ! creating new structure tree
  allocate(mother(nstruct),first_sister(nstruct),first_daughter(nstruct),level(nstruct))
  mother         = -1
  first_sister   = -1
  first_daughter = -1
  level          = 0
  ihalo          = -1 

  do inode = 1, nnodes
     write(errunit,*) inode, node_to_struct(inode)
     istruct = node_to_struct(inode)
     if(istruct.le.0) stop 'index nil for istrut'
     if(mother(istruct).lt.0) then
        if(node(inode)%mother.le.0) then
           mother(istruct) = 0
           if(ihalo.gt.0) first_sister(ihalo) = istruct
           ihalo = istruct
           level(istruct) = 1
        else
           imother          = node_to_struct(node(inode)%mother)
           mother(istruct)  = imother
           level(istruct)   = level(imother) + 1
           if(first_daughter(imother).le.0) then
              first_daughter(imother) = istruct
           else
              isub = first_daughter(imother)
              do while(first_sister(isub).gt.0)
                 isub = first_sister(isub)
              end do
              first_sister(isub) = istruct
           end if
        end if
     end if
  end do

  if(megaverbose) then
     ! For test we shall output the structure tree in file struct_tree.dat'
     open(unit=12,form='formatted',status='unknown',file='struct_tree.dat')
     do istruct = 1, nstruct
        if(mother(istruct).le.0) then
           write(12,*) '---------'
           write(12,'(a6,1x,i6,2x,a12,1x,i6,2x,a7,1x,i6)')                & 
                ' halo:',istruct, 'first_child:',first_daughter(istruct), &
                'sister:', first_sister(istruct) 
           isub = first_daughter(istruct)
           do while(isub.gt.0)
              write(12,'(a5,1x,i6,2x,a7,1x,i6,2x,a12,1x,i6,2x,a7,1x,i6)') &
                   ' sub:',isub, 'mother:', mother(isub), 'first_child:', &
                   first_daughter(isub), 'sister:', first_sister(isub)
              isub = first_sister(isub)
           end do
        else
           isub = first_daughter(istruct)
           if(isub.ne.0) then
              write(12,'(a5,1x,i6,2x,a7,1x,i6,2x,a12,1x,i6,2x,a7,1x,i6)')      &
                   ' sub:',istruct, 'mother:', mother(istruct), 'first_child:',&
                   first_daughter(istruct), 'sister:', first_sister(istruct)
              do while(isub.gt.0)
                 write(12,'(a5,1x,i6,2x,a7,1x,i6,2x,a12,1x,i6,2x,a7,1x,i6)')   &
                      ' sub:',isub, 'mother:', mother(isub), 'first_child:',   &
                      first_daughter(isub), 'sister:', first_sister(isub)
                 isub = first_sister(isub)
              end do
           end if
        end if
     end do
     close(12)
  end if

  nb_of_halos    = nb_halos
  nb_of_subhalos = nb_sub

  deallocate(node)
  
  return

end subroutine select_with_DP_method

!=======================================================================
subroutine select_with_MS_method
!=======================================================================

  implicit none
  integer(kind=4) :: ihalo,inode,isub,ip
  integer(kind=4) :: imother,istruct 
  integer(kind=4) :: nb_halos,nb_sub
  integer(kind=4) :: mostmasssub(nnodes),node_to_struct(nnodes)
  real(kind=8)    :: maxmasssub,mass_acc
!  real(kind=8)    : masstest
  integer(kind=4), allocatable :: npartcheck(:)

  mass_acc = 1.0e-2

  if(verbose) write(errunit,*) 'Using MS method'

  ! First recompute node mass accordingly to their particle list
  do inode = 1, nnodes
     isub = node(inode)%firstchild
     do while(isub.gt.0)
        node(inode)%mass = node(inode)%mass - node(isub)%mass
        node(inode)%truemass = node(inode)%truemass - node(isub)%truemass
        isub = node(isub)%sister
     end do
     if(node(inode)%mass.le.0.or.node(inode)%truemass.le.0.0) then
        write(errunit,*) ' Error in computing node',inode,'mass'
        write(errunit,*) ' mass, truemass:', node(inode)%mass, node(inode)%truemass
        stop
     end if
  end do
  if(verbose) then
     ! check that the new masses are correct
     allocate(npartcheck(0:nnodes))
     npartcheck = 0
     do ip = 1, nbodies
        if(liste_parts(ip).lt.0) stop 'liste_parts is smaller than 0'
        npartcheck(liste_parts(ip)) = npartcheck(liste_parts(ip)) + 1
     end do
     if(sum(npartcheck).ne.nbodies) stop 'Error in particles count'
     do inode = 1, nnodes
        if(node(inode)%mass.ne.npartcheck(inode)) then
           write(errunit,*) 'Error in node particle count, for node',inode
           write(errunit,*) 'it first subnode is:',node(inode)%firstchild
           if(node(inode)%firstchild.gt.0) write(errunit,*) 'it has:',node(node(inode)%firstchild)%nsisters,'subnodes' 
           stop
        end if
     end do
     deallocate(npartcheck)
  end if
  
  mostmasssub    = -1
  do inode = nnodes, 1, -1
     ! Search inode subnodes for the most massive one
     isub = node(inode)%firstchild
     if(isub.gt.0) then
        ! init all on the first subnode
        maxmasssub         = node(isub)%truemass
        mostmasssub(inode) = isub
        isub               = node(isub)%sister
     end if
     do while(isub.gt.0)
        !masstest = (node(isub)%truemass / maxmasssub) - 1.0
        !if((masstest.gt.mass_acc)) then.or. &
        !     ((abs(masstest).le.mass_acc).and.node(isub)%densmax.gt.node(mostmasssub(inode))%densmax)) then
        if(node(isub)%truemass.gt.maxmasssub) then
           maxmasssub         = node(isub)%truemass
           mostmasssub(inode) = isub
        end if
        isub               = node(isub)%sister
     end do
     ! add mostmasssub mass to inode mass and inode's densmax = mostmasssub's densmax
     if (mostmasssub(inode) .gt. 0) then 
        node(inode)%mass     = node(inode)%mass + node(mostmasssub(inode))%mass
        node(inode)%truemass = node(inode)%truemass + node(mostmasssub(inode))%truemass
        node(inode)%densmax  = node(mostmasssub(inode))%densmax
     endif
  end do

  ! Second run to write node_to_struct array and count the number of substructures
  nb_sub         = 0
  nb_halos       = 0
  istruct        = 0
  node_to_struct = -1
 
  do inode = 1, nnodes
     if(node_to_struct(inode).gt.0) stop 'node_to_struct is greater than 0'
     if(node(inode)%mother.le.0) then
        nb_halos              = nb_halos + 1
        istruct               = istruct  + 1
        node_to_struct(inode) = istruct
     else
        imother  = node(inode)%mother
        if(inode.eq.mostmasssub(imother)) then
           if(node_to_struct(imother).le.0) stop 'node_to_struct not defined for imother'
           node_to_struct(inode) = node_to_struct(imother)
        else
           nb_sub                = nb_sub + 1
           istruct               = istruct + 1
           node_to_struct(inode) = istruct
        end if
     end if
  end do
  
  ! check halos and subhalos number 
  nstruct = nb_halos + nb_sub
  if(nstruct.ne.istruct) then
     write(errunit,*) 'Error in subroutines count,',nstruct, istruct
     stop
  end if

  if(verbose) then
     ! Check that the structure tree is ordered as it should be 
     do inode = 1, nnodes
        if(inode.le.nb_halos) then
           if(.not.((node(inode)%level.eq.1).and.(node_to_struct(inode).eq.inode))) stop 'error in node order' 
        else
           if((node(inode)%level.eq.1)) stop 'error shouldn''t be any halo here' 
        end if
     end do
     write(errunit,*)
     write(errunit,*) '> number of nodes               :', nnodes
     write(errunit,*) '> number of halos               :', nb_halos
     write(errunit,*) '> number of substructures       :', nb_sub
     write(errunit,*) '> number of node removed        :', nnodes - nstruct
     write(errunit,*)   
     write(errunit,*) '> Cleaning liste_parts'
     allocate(npartcheck(nstruct))
     npartcheck = 0
  end if

  ! Cleaning liste_parts
  do ip = 1, nbodies
     inode  = liste_parts(ip)
     if(inode.gt.0) then
        if(node_to_struct(inode).le.0) then
           write(errunit,*) ip, inode,  node_to_struct(inode)
           stop 'error in node_to_struct'
        end if
        liste_parts(ip) = node_to_struct(inode)
        if(verbose) npartcheck(node_to_struct(inode)) = npartcheck(node_to_struct(inode)) + 1 
     end if
  end do
  if(verbose) then
     ! Check node(inode)%mass it should now correspond to npartcheck count
     do inode = 1, nnodes
        if(node_to_struct(inode).le.0) stop 'node_to_struct is nil'
        imother = node(inode)%mother
        if((imother.le.0).or.(imother.gt.0.and.(node_to_struct(imother).ne.node_to_struct(inode)))) then
           if(node(inode)%mass.ne.npartcheck(node_to_struct(inode))) then
              write(errunit,*) 'Wrong nb of part in struct: ', node_to_struct(inode)
              write(errunit,*) 'inode,node(inode)%mass,istruct,npartcheck(istruct)',inode,&
                   node(inode)%mass,node_to_struct(inode),npartcheck(node_to_struct(inode))
              stop
           end if
        end if
     end do
     deallocate(npartcheck)
  end if

  if(verbose) write(errunit,*) '> Creating new structure tree'
  ! creating new structure tree
  allocate(mother(nstruct), first_sister(nstruct), first_daughter(nstruct), level(nstruct)) 
  mother         = -1
  first_sister   = -1
  first_daughter = -1
  level          = 0
  ihalo          = -1 

  do inode = 1, nnodes
     istruct = node_to_struct(inode)
     if(istruct.le.0) stop 'index nil for istrut'
     if(mother(istruct).lt.0) then
        if(node(inode)%mother.le.0) then
           mother(istruct) = 0
           if(ihalo.gt.0) first_sister(ihalo) = istruct
           ihalo = istruct
           level(istruct) = 1
        else
           imother          = node_to_struct(node(inode)%mother)
           mother(istruct)  = imother
           level(istruct)   = level(imother) + 1
           if(first_daughter(imother).le.0) then
              first_daughter(imother) = istruct
           else
              isub = first_daughter(imother)
              do while(first_sister(isub).gt.0)
                 isub = first_sister(isub)
              end do
              first_sister(isub) = istruct
           end if
        end if
     end if
  end do

  if(megaverbose) then
     ! For test we shall output the structure tree in file struct_tree.dat'
     open(unit=12,form='formatted',status='unknown',file='struct_tree.dat')
     do istruct = 1, nstruct
        if(mother(istruct).le.0) then
           write(12,*) '---------'
           write(12,'(a6,1x,i6,2x,a12,1x,i6,2x,a7,1x,i6)') ' halo:',istruct,&
                'first_child:',first_daughter(istruct), 'sister:', first_sister(istruct) 
           isub = first_daughter(istruct)
           do while(isub.gt.0)
              write(12,'(a5,1x,i6,2x,a7,1x,i6,2x,a12,1x,i6,2x,a7,1x,i6)') &
                   ' sub:',isub, 'mother:', mother(isub), 'first_child:', &
                   first_daughter(isub),'sister:', first_sister(isub)
              isub = first_sister(isub)
           end do
        else
           isub = first_daughter(istruct)
           if(isub.ne.0) then
              write(12,'(a5,1x,i6,2x,a7,1x,i6,2x,a12,1x,i6,2x,a7,1x,i6)')      &
                   ' sub:',istruct, 'mother:', mother(istruct), 'first_child:',&
                   first_daughter(istruct), 'sister:', first_sister(istruct)
              do while(isub.gt.0)
                 write(12,'(a5,1x,i6,2x,a7,1x,i6,2x,a12,1x,i6,2x,a7,1x,i6)')   &
                      ' sub:',isub, 'mother:', mother(isub), 'first_child:',   &
                      first_daughter(isub), 'sister:', first_sister(isub)
                 isub = first_sister(isub)
              end do
           end if
        end if
     end do
     close(12)
  end if

  nb_of_halos    = nb_halos
  nb_of_subhalos = nb_sub

  deallocate(node)
  
  return

end subroutine select_with_MS_method

!=======================================================================
subroutine select_with_BH_method
!=======================================================================

  use input_output
  implicit none

  integer(kind=4) :: inode,imother,isub,ihalo,istruct,isubmassmax
  integer(kind=4) :: ip
  integer(kind=4) :: n_removed, node_to_struct(nnodes)
  integer(kind=4) :: subprog,prog,rsubprog,n_extrasub
  real(kind=8)    :: mass_acc,masssubmax
!  real(kind=8)    :: masstest
  integer(kind=4), allocatable :: npartcheck(:)
  character(len=20) :: filelast
 
  mass_acc = 1.0e-2

  if(numero_step == 1) then
     allocate(ex_liste_parts(nbodies))
     if(numstep.eq.1) then
        if(verbose) write(errunit,*) '> using MS method for the first step'
        call select_with_MS_method
        if(verbose) write(errunit,*) '> done saving liste_part and nb_of_struct, exlevel'
        ex_liste_parts   = liste_parts
        ex_nb_of_structs = nstruct
        allocate(ex_level(ex_nb_of_structs))
        ex_level         = level
        if(verbose) write(errunit,*) '> done returning'
        return
     else
        write(filelast,'(a11,i3.3)') 'tree_bricks',numstep - 1
        call read_last_brick(filelast)
     end if
  end if

  if(ex_nb_of_structs .le.0) then
     write(errunit,*) '> no structures at the previous step switching to MS method'
     call select_with_MS_method
     ex_liste_parts   = liste_parts
     ex_nb_of_structs = nstruct
     allocate(ex_level(ex_nb_of_structs))
     ex_level         = level
     return
  end if

  if(verbose) then
     write(errunit,*) 
     write(errunit,*) '> Selecting subhalos using previous step results'
  end if

  ! First recompute node mass accordingly to their particle list
  do inode = 1, nnodes
     isub = node(inode)%firstchild
     do while(isub.gt.0)
        node(inode)%mass     = node(inode)%mass - node(isub)%mass
        node(inode)%truemass = node(inode)%truemass - node(isub)%truemass
        isub                 = node(isub)%sister
     end do
     if(node(inode)%mass.le.0.or.node(inode)%truemass.le.0.0) then
        write(errunit,*) ' Error in computing node',inode,'mass'
        write(errunit,*) ' mass, truemass:', node(inode)%mass, node(inode)%truemass
        stop
     end if
  end do
  if(verbose) then
     ! check that the new masses are correct
     allocate(npartcheck(0:nnodes))
     npartcheck = 0
     do ip = 1, nbodies
        if(liste_parts(ip).lt.0) stop 'liste_parts is smaller than 0'
        npartcheck(liste_parts(ip)) = npartcheck(liste_parts(ip)) + 1
     end do
     if(sum(npartcheck).ne.nbodies) stop 'Error in particles count'
     do inode = 1, nnodes
        if(node(inode)%mass.ne.npartcheck(inode)) then
           write(errunit,*) 'Error in node particle count, for node',inode
           write(errunit,*) 'it first subnode is:',node(inode)%firstchild
           if(node(inode)%firstchild.gt.0) write(errunit,*) 'it has:',node(node(inode)%firstchild)%nsisters,'subnodes' 
           stop
        end if
     end do
     deallocate(npartcheck)
  end if

  ! create a linked list of particles belonging to nodes
  allocate(first_part(0:nnodes),linked_list(nbodies),npfather(0:ex_nb_of_structs))
  call make_linked_node_list

  if(verbose) write(errunit,*) '> Searching for each node a subnode to remove'
  allocate(removesub(nnodes))
  removesub  = -1
  n_extrasub = 0
  do inode = nnodes, 1, -1
     isub = node(inode)%firstchild
     if(isub.gt.0) then
        if(megaverbose) write(errunit,*) '> node:', inode,'first subnode:',isub, 'nsub:', node(isub)%nsisters
        ! first search most massive sub
        isubmassmax = isub
        masssubmax  = node(isubmassmax)%truemass
        isub        = node(isub)%sister
        do while(isub.gt.0)
           if(node(isub)%truemass.gt.masssubmax) then
              isubmassmax = isub
              masssubmax  = node(isubmassmax)%truemass
           end if
           isub = node(isub)%sister
        end do
        if(megaverbose) write(errunit,*) '> most massive subnode:', isubmassmax, masssubmax

        ! search for the subnode to remove
        if(removesub(inode).gt.0) stop 'one sub already removed'
        isub = node(inode)%firstchild
        do while(isub.gt.0)
           if((node(inode)%truemass + node(isub)%truemass).gt.masssubmax) then
              if(megaverbose) write(errunit,*)'> removable subnode:',isub
              ! if isub was removed inode mass would be greater than any of its subhaloes mass
              ! compute main prog of isub
              call main_prog(isub,subprog)
              if(megaverbose) write(errunit,*) '> subnode mainprog:',subprog
              ! compute main prog of isub + inode
              call main_prog_node_sub(inode,isub,prog)
              if(megaverbose) write(errunit,*) '> subnode + node mainprog:',prog
              if(subprog.eq.prog) then
                 if(removesub(inode).le.0) then
                    if(megaverbose) write(errunit,*) '> remove subnode',isub,' for now'
                    removesub(inode) = isub
                    rsubprog         = subprog
                 else
                    if(megaverbose) write(errunit,*) '> one sub has been removed:',removesub(inode),'prog:',rsubprog
                    if(rsubprog.eq.prog) then
                       if(megaverbose) write(errunit,*) '> choose with MS method'
                       ! use MS method
                       !masstest = (node(isub)%truemass / node(removesub(inode))%truemass) - 1.0
                       !if((masstest.gt.mass_acc).or. &
                       !     ((abs(masstest).le.mass_acc).and.node(isub)%densmax.gt.node(removesub(inode))%densmax)) then
                       if(node(isub)%truemass.gt.node(removesub(inode))%truemass) then
                          removesub(inode) = isub
                          rsubprog         = subprog
                          if(megaverbose) write(errunit,*) '> remove subnode',isub,' instead'
                       end if
                    else
                       if(megaverbose) then
                          write(errunit,*) '> Using progs level criteria'
                          write(errunit,*) '> prog:',ex_level(prog),'rsubprog:',ex_level(rsubprog)
                       end if
                       if(ex_level(prog).lt.ex_level(rsubprog)) then
                          removesub(inode) = isub
                          rsubprog         = subprog
                          if(megaverbose) write(errunit,*) '> remove subnode',isub,' instead'
                       else if(ex_level(prog).eq.ex_level(rsubprog)) then
                          if(megaverbose) write(errunit,*) '> same level, use MS method'
                          ! apply MS for now but might better use the most massive of the two progenitors
                          ! masstest = (node(isub)%truemass / node(removesub(inode))%truemass) - 1.0
                          !if((masstest.gt.mass_acc).or.                              &
                          !    ((abs(masstest).le.mass_acc).and.                      & 
                          !    (node(isub)%densmax.gt.node(removesub(inode))%densmax))) then
                          if(node(isub)%truemass.gt.node(removesub(inode))%truemass) then
                             removesub(inode) = isub
                             rsubprog         = subprog
                             if(megaverbose) write(errunit,*) '> remove subnode',isub,' instead'
                          end if
                       end if
                    end if
                 end if
              end if
           end if
           isub = node(isub)%sister
        end do
        
        if(removesub(inode).le.0) then
           if(node(inode)%firstchild.gt.0) n_extrasub = n_extrasub + 1
           if(megaverbose) write(errunit,*) '> no subnode removed'
           ! none sub has been removed, check that inode is more massive than any of its subhalos if not remove most massive subhalo
           if(node(inode)%truemass.le.node(isubmassmax)%truemass) then
              if(megaverbose) write(errunit,*) '> have to remove something, remove most massive sub'
              removesub(inode) = isubmassmax
              if(megaverbose) write(errunit,*) '> remove subnode',isubmassmax
           end if
        end if
        if(removesub(inode).gt.0) then
           ! add removedsub mass to inode change densmax
           node(inode)%mass = node(inode)%mass + node(removesub(inode))%mass 
           node(inode)%truemass = node(inode)%truemass + node(removesub(inode))%truemass
           node(inode)%densmax  = node(removesub(inode))%densmax
        end if
        if(megaverbose) write(errunit,*)
     end if
  end do
  deallocate(first_part,linked_list,npfather)

  node_to_struct = -1
  nstruct         = 0
  n_removed       = 0 
  nb_of_halos     = 0
  nb_of_subhalos  = 0
  ! write node_to_struct array
  do inode = 1, nnodes
     if(node_to_struct(inode).gt.0) stop 'node_to_struct is greater than 0'
     imother = node(inode)%mother
     if(imother.eq.0) then
        if(nb_of_subhalos.gt.0) stop 'All subhalos should have been delt with'
        nb_of_halos            = nb_of_halos + 1
        nstruct               = nstruct +1
        node_to_struct(inode) = nstruct
     else
        if(node_to_struct(imother).le.0) stop 'node_to_struct doesn''t exist for imother'
        if(inode.eq.removesub(imother)) then
           node_to_struct(inode) = node_to_struct(imother)
           n_removed             = n_removed + 1
        else
           nb_of_subhalos        = nb_of_subhalos + 1
           nstruct               = nstruct +1
           node_to_struct(inode) = nstruct
        end if
     end if
     if(node_to_struct(inode).le.0) stop 'node_to_struct is smaller than 0'
  end do

  if(verbose) then
     write(errunit,*) '> number of nodes            :',nnodes
     write(errunit,*) '> number of structs          :',nstruct
     write(errunit,*) '> number of haloes           :',nb_of_halos
     write(errunit,*) '> number of subhaloes        :',nb_of_subhalos
     write(errunit,*) '> number nodes removed       :',n_removed
     write(errunit,*) '> number of extra subhaloes  :',n_extrasub
     write(465,*) numstep,n_extrasub
  end if

  if(nstruct.ne.(nb_of_halos+nb_of_subhalos)) stop 'error in halo and subhalo count'

 ! Cleaning liste_parts
  if(verbose) then
     allocate(npartcheck(nstruct))
     npartcheck = 0
  end if
  do ip = 1, nbodies
     inode  = liste_parts(ip)
     if(inode.gt.0) then
        if(node_to_struct(inode).le.0) then
           write(errunit,*) ip, inode,  node_to_struct(inode)
           stop 'error in node_to_struct'
        end if
        liste_parts(ip) = node_to_struct(inode)
        if(verbose) npartcheck(node_to_struct(inode)) = npartcheck(node_to_struct(inode)) + 1 
     end if
  end do
  if(verbose) then
     ! Check node(inode)%mass it should now correspond to npartcheck count
     do inode = 1, nnodes
        if(node_to_struct(inode).le.0) stop 'node_to_struct is nil'
        imother = node(inode)%mother
        if((imother.le.0).or.(imother.gt.0.and.(node_to_struct(imother).ne.node_to_struct(inode)))) then
           if(node(inode)%mass.ne.npartcheck(node_to_struct(inode))) then
              write(errunit,*) 'Wrong nb of part in struct: ', node_to_struct(inode)
              write(errunit,*) 'inode,node(inode)%mass,istruct,npartcheck(istruct)',inode, &
              node(inode)%mass,node_to_struct(inode),npartcheck(node_to_struct(inode))
              stop
           end if
        end if
     end do
     deallocate(npartcheck)
  end if

  if(verbose) write(errunit,*) '> Creating new structure tree'
  ! creating new structure tree
  allocate(mother(nstruct), first_sister(nstruct), first_daughter(nstruct), level(nstruct)) 
  mother         = -1
  first_sister   = -1
  first_daughter = -1
  level          = 0
  ihalo          = -1 

  do inode = 1, nnodes
     istruct = node_to_struct(inode)
     if(istruct.le.0) stop 'index nil for istrut'
     if(mother(istruct).lt.0) then
        if(node(inode)%mother.le.0) then
           mother(istruct) = 0
           if(ihalo.gt.0) first_sister(ihalo) = istruct
           ihalo          = istruct
           level(istruct) = 1
        else
           imother          = node_to_struct(node(inode)%mother)
           mother(istruct)  = imother
           level(istruct)   = level(imother) + 1
           if(first_daughter(imother).le.0) then
              first_daughter(imother) = istruct
           else
              isub = first_daughter(imother)
              do while(first_sister(isub).gt.0)
                 isub = first_sister(isub)
              end do
              first_sister(isub) = istruct
           end if
        end if
     end if
  end do

  if(megaverbose) then
     ! For test we shall output the structure tree in file struct_tree.dat'
     open(unit=12,form='formatted',status='unknown',file='struct_tree_BHM.dat')
     do istruct = 1, nstruct
        if(mother(istruct).le.0) then
           write(12,*) '---------'
           write(12,'(a6,1x,i6,2x,a12,1x,i6,2x,a7,1x,i6)') ' halo:',istruct, &
                'first_child:',first_daughter(istruct), 'sister:', first_sister(istruct) 
           isub = first_daughter(istruct)
           do while(isub.gt.0)
              write(12,'(a5,1x,i6,2x,a7,1x,i6,2x,a12,1x,i6,2x,a7,1x,i6)') &
                   ' sub:',isub, 'mother:', mother(isub), 'first_child:', &
                   first_daughter(isub), 'sister:', first_sister(isub)
              isub = first_sister(isub)
           end do
        else
           isub = first_daughter(istruct)
           if(isub.ne.0) then
              write(12,'(a5,1x,i6,2x,a7,1x,i6,2x,a12,1x,i6,2x,a7,1x,i6)')      &
                   ' sub:',istruct, 'mother:', mother(istruct), 'first_child:',&
                   first_daughter(istruct), 'sister:', first_sister(istruct)
              do while(isub.gt.0)
                 write(12,'(a5,1x,i6,2x,a7,1x,i6,2x,a12,1x,i6,2x,a7,1x,i6)')    &
                      ' sub:',isub, 'mother:', mother(isub), 'first_child:',    &
                      first_daughter(isub), 'sister:', first_sister(isub)
                 isub = first_sister(isub)
              end do
           end if
        end if
     end do
     close(12)
  end if

  deallocate(removesub)
  deallocate(node)

  if(numero_step == nsteps) then
     deallocate(ex_liste_parts)
     deallocate(ex_level)
  else
     ex_liste_parts   = liste_parts 
     ex_nb_of_structs = nstruct
     deallocate(ex_level)
     allocate(ex_level(ex_nb_of_structs))
     ex_level(1:ex_nb_of_structs) = level(1:ex_nb_of_structs)
  end if
  
  return

end subroutine select_with_BH_method

!=======================================================================
subroutine main_prog(inode,mainprog)
!=======================================================================

  implicit none
  integer(kind=4) :: inode, mainprog
  integer(kind=4) :: isubremoved,inodetmp,ip,npprog
 

  npfather(0:ex_nb_of_structs) = 0
  mainprog = -1 
  ip = first_part(inode)
  do while(ip.gt.0)
     if(ex_liste_parts(ip).gt.ex_nb_of_structs) stop 'Error in ex_liste_parts'
     npfather(ex_liste_parts(ip)) = npfather(ex_liste_parts(ip)) + 1
     ip = linked_list(ip)
  end do

  isubremoved = removesub(inode)
  do while(isubremoved.gt.0)
     ip = first_part(isubremoved)
     do while(ip.gt.0)
        if(ex_liste_parts(ip).gt.ex_nb_of_structs) stop 'Error in ex_liste_parts'
        npfather(ex_liste_parts(ip)) = npfather(ex_liste_parts(ip)) + 1
        ip = linked_list(ip)
     end do
     isubremoved = removesub(isubremoved)
  end do

  npprog = 0
  do inodetmp = 0, ex_nb_of_structs
     if(npfather(inodetmp).gt.npprog) then
        npprog   = npfather(inodetmp)
        mainprog = inodetmp
     end if
  end do

  return

end subroutine main_prog

!=======================================================================
subroutine main_prog_node_sub(inode,isub,mainprog)
!=======================================================================

  implicit none
  integer(kind=4) :: inode,isub, mainprog
  integer(kind=4) :: isubremoved,inodetmp,ip,npprog
 

  npfather(0:ex_nb_of_structs) = 0
  mainprog = -1 
  ip = first_part(inode)
  do while(ip.gt.0)
     if(ex_liste_parts(ip).gt.ex_nb_of_structs) stop 'Error in ex_liste_parts'
     npfather(ex_liste_parts(ip)) = npfather(ex_liste_parts(ip)) + 1
     ip = linked_list(ip)
  end do

  ip = first_part(isub)
  do while(ip.gt.0)
     if(ex_liste_parts(ip).gt.ex_nb_of_structs) stop 'Error in ex_liste_parts'
     npfather(ex_liste_parts(ip)) = npfather(ex_liste_parts(ip)) + 1
     ip = linked_list(ip)
  end do

  isubremoved = removesub(isub)
  do while(isubremoved.gt.0)
     ip = first_part(isubremoved)
     do while(ip.gt.0)
        if(ex_liste_parts(ip).gt.ex_nb_of_structs) stop 'Error in ex_liste_parts'
        npfather(ex_liste_parts(ip)) = npfather(ex_liste_parts(ip)) + 1
        ip = linked_list(ip)
     end do
     isubremoved = removesub(isubremoved)
  end do

  npprog = 0
  do inodetmp = 0, ex_nb_of_structs
     if(npfather(inodetmp).gt.npprog) then
        npprog   = npfather(inodetmp)
        mainprog = inodetmp
     end if
  end do

  return

end subroutine main_prog_node_sub

!=======================================================================
subroutine make_linked_node_list
!=======================================================================

  implicit none
  integer(kind=4) :: index1, index2,i
  integer(kind=4), allocatable :: current_ptr(:)

  allocate(current_ptr(0:nnodes))
  current_ptr  = -1
  first_part  = -1
  linked_list = -1

  do i = 1, nbodies
     index1 = liste_parts(i)
     if(first_part(index1).le.0) then
        first_part(index1) = i
        if(allocated(nb_of_parts)) nb_of_parts(index1) = 1
     else
        index2 = current_ptr(index1)
        linked_list(index2) = i
        if(allocated(nb_of_parts)) nb_of_parts(index1) = nb_of_parts(index1) + 1 
     end if
     current_ptr(index1) = i
  end do

  do i = 0,nnodes
     if(current_ptr(i) ==-1) cycle
     index2              = current_ptr(i)
     linked_list(index2) = -1 
  end do
  
  deallocate(current_ptr)
  return

end subroutine make_linked_node_list

!=======================================================================
subroutine compute_mean_density_and_np
!=======================================================================

  implicit none

  integer(kind=4)                     :: ipar
  real(kind=8), dimension(0:nvoisins) :: dist2
  integer, dimension(nvoisins)        :: iparnei
  real(kind=8)                        :: densav

  if (verbose) write(errunit,*) 'Compute mean density for each particle...'

  allocate(iparneigh(nhop,npart))
  allocate(density(npart))

  if (verbose) write(errunit,*) 'First find nearest particles'

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ipar,dist2,iparnei)
  do ipar=1,npart
     call find_nearest_parts(ipar,dist2,iparnei)
     call compute_density(ipar,dist2,iparnei)
     iparneigh(1:nhop,ipar)=iparnei(1:nhop)
  enddo
!$OMP END PARALLEL DO

! Check for average density
  if (verbose) then
     densav=0.d0
     do ipar=1,npart
        densav=(densav*dble(ipar-1)+dble(density(ipar)))/dble(ipar)
     enddo
     write(errunit,*) 'Average density :',densav
  endif
  
  deallocate(mass_cell)
  deallocate(size_cell)
  deallocate(pos_cell)
  deallocate(sister)
  deallocate(firstchild)

end subroutine compute_mean_density_and_np

!=======================================================================
subroutine find_local_maxima
!=======================================================================

  implicit none

  integer(kind=4)             :: ipar,idist,iparid,iparsel,igroup,nmembmax,nmembtot
  integer(kind=4),allocatable :: nmemb(:)
  real(kind=8)                :: denstest

  if (verbose) write(errunit,*) 'Now Find local maxima...'

  allocate(igrouppart(npart))

  idpart=0
  ngroups=0
  do ipar=1,npart
     denstest=density(ipar)
     if (denstest.gt.rho_threshold) then
        iparsel=ipar
        do idist=1,nhop
           iparid=iparneigh(idist,ipar)
           if (density(iparid).gt.denstest) then
              iparsel=iparid
              denstest=density(iparid)
           elseif (density(iparid).eq.denstest) then
              iparsel=min(iparsel,iparid)
              if (verbose) &
 &               write(errunit,*) 'WARNING : equal densities in find_local_maxima.'
           endif
        enddo
        if (iparsel.eq.ipar) then
           ngroups=ngroups+1
           idpart(ipar)=-ngroups
        else
           idpart(ipar)=iparsel
        endif
     endif
  enddo
 
  if (verbose) write(errunit,*) 'Number of local maxima found :',ngroups

! Now Link the particles associated to the same maximum
  if (verbose) write(errunit,*) 'Create link list...'

  allocate(densityg(ngroups))
  allocate(nmemb(ngroups))
  allocate(firstpart(ngroups))

  do ipar=1,npart
     if (density(ipar).gt.rho_threshold) then
        iparid=idpart(ipar)
        if (iparid.lt.0) densityg(-iparid)=density(ipar)
        do while (iparid.gt.0)
           iparid=idpart(iparid)
        enddo
        igrouppart(ipar)=-iparid
     else
        igrouppart(ipar)=0
     endif
  enddo

  nmemb(1:ngroups)=0
  firstpart(1:ngroups)=0
  do ipar=1,npart
     igroup=igrouppart(ipar)
     if (igroup.gt.0) then
        idpart(ipar)=firstpart(igroup)   
        firstpart(igroup)=ipar
        nmemb(igroup)=nmemb(igroup)+1
     endif
  enddo

  nmembmax=0
  nmembtot=0
  do igroup=1,ngroups
     nmembmax=max(nmembmax,nmemb(igroup))
     nmembtot=nmembtot+nmemb(igroup)
  enddo

  if (verbose) then
     write(errunit,*) 'Number of particles of the largest group :',nmembmax
     write(errunit,*) 'Total number of particles in groups ',nmembtot
  endif
  deallocate(nmemb)
  
end subroutine find_local_maxima

!=======================================================================
subroutine create_group_tree
!=======================================================================

  implicit none

  integer(kind=4) :: inode,mass_loc,masstmp,igroup,igr1,igr2,igroupref
  real(kind=8)    :: rhot,posg(3),posref(3),rsquare,densmoy,truemass,truemasstmp

  if (verbose) write(errunit,*) 'Create the tree of structures of structures'

! End of the branches of the tree
  call compute_saddle_list

  if (verbose) write(errunit,*) 'Build the hierarchical tree'

  nnodesmax=2*ngroups
! Allocations
  allocate(node(0:nnodesmax))
  allocate(idgroup(ngroups))
  allocate(color(ngroups))
  allocate(igroupid(ngroups))
  allocate(idgroup_tmp(ngroups))

! Initializations
  liste_parts = 0

! Iterative loop to build the rest of the tree
  inode              = 0
  nnodes             = 0
  rhot               = rho_threshold
  node(inode)%mother = 0
  mass_loc           = 0
  truemass           = 0.d0
  igroupref          = 0
  do igroup=1,ngroups
     call treat_particles(igroup,rhot,posg,masstmp,igroupref,posref, &
&                          rsquare,densmoy,truemasstmp)
     mass_loc = mass_loc+masstmp
     truemass = truemass+truemasstmp
  enddo

  node(inode)%mass          = mass_loc
  node(inode)%truemass      = truemass
  node(inode)%radius        = 0.d0
  node(inode)%density       = 0.d0
  node(inode)%position(1:3) = 0.d0
  node(inode)%densmax       = maxval(densityg(1:ngroups))
  node(inode)%rho_saddle    = 0.
  node(inode)%level         = 0
  node(inode)%nsisters      = 0
  node(inode)%sister        = 0
  igr1 = 1
  igr2 = ngroups
  do igroup=1,ngroups
     idgroup(igroup)=igroup
     igroupid(igroup)=igroup
  enddo
  call create_nodes(rhot,inode,igr1,igr2)

  deallocate(idgroup)
  deallocate(color)
  deallocate(igroupid)
  deallocate(idgroup_tmp)
  deallocate(idpart)
  deallocate(group)
  deallocate(densityg)
  deallocate(firstpart)
 
end subroutine create_group_tree

!=======================================================================
recursive subroutine create_nodes(rhot,inode,igr1,igr2)
!=======================================================================

  implicit none
  integer(kind=4)              :: inode,igr1,igr2
  real(kind=8)                 :: rhot
  integer(kind=4)              :: igroup,icolor,igr,igr_eff
  integer(kind=4)              :: inc_color_tot,inode1,mass_loc,masstmp,igroupref
  integer(kind=4)              :: inodeout,igr1out,igr2out,isisters,nsisters
  integer(kind=4)              :: mass_comp,icolor_ref
  real(kind=8)                 :: posg(3),posgtmp(3),posref(3),rhotout,rsquaretmp,rsquareg
  real(kind=8)                 :: densmoyg,densmoytmp
  real(kind=8)                 :: densmoy_comp_max,truemass,truemasstmp
  real(kind=8)                 :: posfin(3)
  real(kind=8)                 :: densmaxgroup
  integer(kind=4), allocatable :: igrpos(:),igrinc(:)
  integer(kind=4), allocatable :: igrposnew(:)
  integer(kind=4), allocatable :: massg(:)
  real(kind=8),  allocatable   :: truemassg(:)
  real(kind=8),  allocatable   :: densmaxg(:)
  real(kind=8),  allocatable   :: densmoy_comp_maxg(:)
  integer(kind=4), allocatable :: mass_compg(:)
  real(kind=8),  allocatable   :: posgg(:,:)
  real(kind=8),  allocatable   :: rsquare(:)
  real(kind=8),  allocatable   :: densmoy(:)
  logical, allocatable         :: ifok(:)  

  color(igr1:igr2)=0
! Percolate the groups
  icolor_select=0
  do igr=igr1,igr2
     igroup=idgroup(igr)
     if (color(igr).eq.0) then
        icolor_select=icolor_select+1
!!$        call do_colorize(icolor_select,igroup,igr,rhot)
        ncall=0                                               !YDdebug
        call do_colorize(igroup,igr,rhot) !YDdebug
        write(errunit,'(A,3I8,e10.2,I8)')'End of do_colorize',icolor_select,igroup,igr,rhot,ncall
     endif
  enddo

! We select only groups where we are sure of having at least one
! particle above the threshold density rhot
! Then sort them to gather them on the list
  allocate(igrpos(0:icolor_select))
  allocate(igrinc(1:icolor_select))
  igrpos(0)=igr1-1
  igrpos(1:icolor_select)=0
  do igr=igr1,igr2
     icolor=color(igr)
     igroup=idgroup(igr)
     if (densityg(igroup).gt.rhot) &
 &           igrpos(icolor)=igrpos(icolor)+1
  enddo
  do icolor=1,icolor_select
     igrpos(icolor)=igrpos(icolor-1)+igrpos(icolor)
  enddo
  if (igrpos(icolor_select)-igr1+1.eq.0) then
     write(errunit,*) 'ERROR in create_nodes :'
     write(errunit,*) 'All subgroups are below the threshold.'
!     STOP
  endif      

  igrinc(1:icolor_select)=0
  do igr=igr1,igr2
     icolor=color(igr)
     igroup=idgroup(igr)
     if (densityg(igroup).gt.rhot) then
        igrinc(icolor)=igrinc(icolor)+1
        igr_eff=igrinc(icolor)+igrpos(icolor-1)
        idgroup_tmp(igr_eff)=igroup
        igroupid(igroup)=igr_eff
     endif
  enddo
  igr2=igrpos(icolor_select)
  idgroup(igr1:igr2)=idgroup_tmp(igr1:igr2)

  inc_color_tot=0
  do icolor=1,icolor_select
     if (igrinc(icolor).gt.0) inc_color_tot=inc_color_tot+1
  enddo 
  allocate(igrposnew(0:inc_color_tot))
  igrposnew(0)=igrpos(0)
  inc_color_tot=0
  do icolor=1,icolor_select
     if (igrinc(icolor).gt.0) then
        inc_color_tot=inc_color_tot+1
        igrposnew(inc_color_tot)=igrpos(icolor)
     endif
  enddo
  deallocate(igrpos)
  deallocate(igrinc)

  isisters=0
  allocate(posgg(1:3,1:inc_color_tot))
  allocate(massg(1:inc_color_tot))
  allocate(truemassg(1:inc_color_tot))
  allocate(densmaxg(1:inc_color_tot))
  allocate(densmoy_comp_maxg(1:inc_color_tot))
  allocate(mass_compg(1:inc_color_tot))
  allocate(rsquare(1:inc_color_tot))
  allocate(densmoy(1:inc_color_tot))
  allocate(ifok(1:inc_color_tot))
  ifok(1:inc_color_tot)=.false.

  do icolor=1,inc_color_tot
     posg(1:3)=0.d0
     mass_loc=0
     truemass=0.d0
     rsquareg=0.d0
     densmoyg=0.d0
     igr1=igrposnew(icolor-1)+1
     igr2=igrposnew(icolor)
     densmaxgroup=-1.
     mass_comp=0
     densmoy_comp_max=-1.
     igroupref=0
     do igr=igr1,igr2
        igroup=idgroup(igr)
        densmaxgroup=max(densmaxgroup,densityg(igroup))
        call treat_particles(igroup,rhot,posgtmp,masstmp, &
 &                           igroupref,posref,rsquaretmp, &
 &                           densmoytmp,truemasstmp)
        posg(1)=posg(1)+posgtmp(1)
        posg(2)=posg(2)+posgtmp(2)
        posg(3)=posg(3)+posgtmp(3)
        rsquareg=rsquareg+rsquaretmp
        mass_loc=mass_loc+masstmp
        truemass=truemass+truemasstmp
        densmoyg=densmoyg+densmoytmp
        densmoytmp=densmoytmp/dble(masstmp)
        mass_comp=max(mass_comp,masstmp)
        if (masstmp > 0) then
           densmoy_comp_max=max(densmoy_comp_max, &
 &             densmoytmp/(1.d0+fudge/sqrt(dble(masstmp))))
        endif
     enddo
     massg(icolor)=mass_loc
     truemassg(icolor)=truemass
     posgg(1:3,icolor)=posg(1:3)
     densmaxg(icolor)=densmaxgroup
     densmoy_comp_maxg(icolor)=densmoy_comp_max
     mass_compg(icolor)=mass_comp
     rsquare(icolor)=sqrt(abs( &
 &             (truemass*rsquareg- &
 &             (posg(1)**2+posg(2)**2+posg(3)**2) )/ &
 &              truemass**2 ))
     densmoy(icolor)=densmoyg/dble(mass_loc)

     ifok(icolor)=mass_loc.ge.nmembthresh.and. &
 &      (densmoy(icolor).gt.rhot*(1.d0+fudge/sqrt(dble(mass_loc))).or. &
 &       densmoy_comp_maxg(icolor).gt.rhot).and. &
 &      densmaxg(icolor).ge.alphap*densmoy(icolor).and. &
 &      rsquare(icolor).ge.epsilon

     if (ifok(icolor)) then
        isisters=isisters+1
        icolor_ref=icolor
     endif
  enddo
  nsisters=isisters
  if (nsisters.gt.1) then
     isisters=0
     inode1=nnodes+1
     do icolor=1,inc_color_tot
        if (ifok(icolor)) then
           isisters=isisters+1
           nnodes=nnodes+1
           if (nnodes.gt.nnodesmax) then
              write(errunit,*) 'ERROR in create_nodes :'
              write(errunit,*) 'nnodes > nnodes max'
              STOP
           endif
           if (mod(nnodes,max(nnodesmax/10000,1)).eq.0.and.megaverbose) then
              write(errunit,*) 'nnodes=',nnodes
           endif
           node(nnodes)%mother=inode
           node(nnodes)%densmax=densmaxg(icolor)
           if (isisters.gt.1) then
              node(nnodes)%sister=nnodes-1
           else
              node(nnodes)%sister=0
           endif
           node(nnodes)%nsisters=nsisters
           node(nnodes)%mass=massg(icolor)
           node(nnodes)%truemass=truemassg(icolor)
           if (mass_loc.eq.0) then
              write(errunit,*) 'ERROR in create_nodes :'
              write(errunit,*) 'NULL mass for nnodes=',nnodes
              STOP
           endif
           posfin(1:3)=real(posgg(1:3,icolor)/truemassg(icolor))
           node(nnodes)%radius=real(rsquare(icolor))    
           node(nnodes)%density=real(densmoy(icolor))
           if (posfin(1).ge.xlongs2) then
              posfin(1)=posfin(1)-xlong
           elseif (posfin(1).lt.-xlongs2) then
              posfin(1)=posfin(1)+xlong
           endif
           if (posfin(2).ge.ylongs2) then
              posfin(2)=posfin(2)-ylong
           elseif (posfin(2).lt.-ylongs2) then
              posfin(2)=posfin(2)+ylong
           endif
           if (posfin(3).ge.zlongs2) then
              posfin(3)=posfin(3)-zlong
           elseif (posfin(3).lt.-zlongs2) then
              posfin(3)=posfin(3)+zlong
           endif
           node(nnodes)%position(1:3)=posfin(1:3)
           node(nnodes)%rho_saddle=rhot
           node(nnodes)%level=node(inode)%level+1
           if (megaverbose.and.node(nnodes)%mass.ge.nmembthresh) then
              write(errunit,*) '*****************************************'
              write(errunit,*) 'new node :',nnodes
              write(errunit,*) 'level    :',node(nnodes)%level
              write(errunit,*) 'nsisters :',node(nnodes)%nsisters
              write(errunit,*) 'mass     :',node(nnodes)%mass
              write(errunit,*) 'true mass:',node(nnodes)%truemass
              write(errunit,*) 'radius   :',node(nnodes)%radius
              write(errunit,*) 'position :',node(nnodes)%position
              write(errunit,*) 'rho_saddl:',node(nnodes)%rho_saddle
              write(errunit,*) 'rhomax   :',node(nnodes)%densmax
              write(errunit,*) '*****************************************'
           endif
        endif
     enddo
     node(inode)%firstchild=nnodes
     inodeout=inode1
     do icolor=1,inc_color_tot
        if (ifok(icolor)) then
           igr1out=igrposnew(icolor-1)+1
           igr2out=igrposnew(icolor)
           do igr=igr1out,igr2out
              call paint_particles(idgroup(igr),inodeout,rhot)
           enddo
           rhotout=rhot*(1.d0+fudge/sqrt(dble(mass_compg(icolor))))
           if (igr2out.ne.igr1out) then
              call create_nodes(rhotout,inodeout,igr1out,igr2out)
           else
              node(inodeout)%firstchild=0
           endif
           inodeout=inodeout+1
        endif
     enddo
  elseif (nsisters.eq.1) then
     inodeout=inode
     rhotout=rhot*(1.d0+fudge/sqrt(dble(mass_compg(icolor_ref))))
     igr1out=igrposnew(0)+1
     igr2out=igrposnew(inc_color_tot)
     if (igr2out.ne.igr1out) then
        call create_nodes(rhotout,inodeout,igr1out,igr2out)
     else
        node(inode)%firstchild=0
     endif
  else
     node(inode)%firstchild=0
  endif
  deallocate(igrposnew)
  deallocate(posgg)
  deallocate(massg)
  deallocate(truemassg)
  deallocate(densmaxg)
  deallocate(densmoy_comp_maxg)
  deallocate(densmoy)
  deallocate(mass_compg)
  deallocate(rsquare)
  deallocate(ifok)

end subroutine create_nodes

!=======================================================================
subroutine paint_particles(igroup,inode,rhot)
!=======================================================================

  implicit none

  integer(kind=4) :: igroup,inode
  real(kind=8)    :: rhot
  integer(kind=4) :: ipar

  ipar=firstpart(igroup)
  do while (ipar.gt.0)
     if (density(ipar).gt.rhot) then
        liste_parts(ipar)=inode
     endif
     ipar=idpart(ipar)
  enddo
end subroutine paint_particles

!=======================================================================
subroutine treat_particles(igroup,rhot,posg,imass,igroupref,posref, &
&                          rsquare,densmoy,truemass)
!=======================================================================

  implicit none

  real(kind=8)    :: rhot
  real(kind=8)    :: posg(3),posref(3)
  real(kind=8)    :: posdiffx,posdiffy,posdiffz,rsquare,densmoy,truemass,xmasspart
  real(kind=8)    :: densmax,densmin
  integer(kind=4) :: imass,ipar,iparold,igroup,igroupref
  logical         :: first_good

  imass=0
  truemass=0.d0
  rsquare=0.d0
  densmoy=0.d0
  posg(1:3)=0.d0
  ipar=firstpart(igroup)
  first_good=.false.
  do while (ipar.gt.0)
     if (density(ipar).gt.rhot) then
        if (.not.first_good) then
           if (igroupref.eq.0) then
              posref(1:3)=dble(pos(ipar,1:3))
              igroupref=igroup
           endif
           first_good=.true.
           firstpart(igroup)=ipar
           densmin=density(ipar)
           densmax=densmin
        else
           idpart(iparold)=ipar
        endif
        iparold=ipar
        imass=imass+1
        if(allocated(mass)) then
           xmasspart=mass(ipar)
        else
           xmasspart=massp
        end if
        truemass=truemass+xmasspart
        posdiffx=dble(pos(ipar,1))-posref(1)
        posdiffy=dble(pos(ipar,2))-posref(2)
        posdiffz=dble(pos(ipar,3))-posref(3)
        if (posdiffx.ge.xlongs2) then
           posdiffx=posdiffx-xlong
        elseif (posdiffx.lt.-xlongs2) then
           posdiffx=posdiffx+xlong
        endif
        if (posdiffy.ge.ylongs2) then
           posdiffy=posdiffy-ylong
        elseif (posdiffy.lt.-ylongs2) then
           posdiffy=posdiffy+ylong
        endif
        if (posdiffz.ge.zlongs2) then
           posdiffz=posdiffz-zlong
        elseif (posdiffz.lt.-zlongs2) then
           posdiffz=posdiffz+zlong
        endif
        posdiffx=posdiffx+posref(1)
        posdiffy=posdiffy+posref(2)
        posdiffz=posdiffz+posref(3)
        posg(1)=posg(1)+posdiffx*xmasspart
        posg(2)=posg(2)+posdiffy*xmasspart
        posg(3)=posg(3)+posdiffz*xmasspart
        rsquare=rsquare+xmasspart*(posdiffx**2+posdiffy**2+posdiffz**2)
        densmoy=densmoy+dble(density(ipar))
        densmax=max(densmax,density(ipar))
        densmin=min(densmin,density(ipar))
     endif
     ipar=idpart(ipar)
  enddo
  if (.not.first_good) firstpart(igroup)=0

  if (densmin.le.rhot.or.densmax.ne.densityg(igroup)) then
     write(errunit,*) 'ERROR in treat_particles'
     write(errunit,*) 'igroup, densmax, rhot=',igroup,densityg(igroup),rhot
     write(errunit,*) 'denslow, denshigh    =',densmin,densmax
     STOP
  endif

end subroutine treat_particles

!=======================================================================
!!$recursive subroutine do_colorize(icolor_select,igroup,igr,rhot)
recursive subroutine do_colorize(igroup,igr,rhot) !YDdebug
!=======================================================================

  implicit none
  integer(kind=4) :: igroup,igr
  integer(kind=4) :: ineig,igroup2,igr2,neig
  real(kind=8)  :: rhot  

  ncall=ncall+1
!!$  write(errunit,'(A,3I8,e10.2,I8)')'do_colorize',icolor_select,igroup,igr,rhot,ncall
  color(igr)=icolor_select
  neig=group(igroup)%nhnei
  do ineig=1,group(igroup)%nhnei
     if (group(igroup)%rho_saddle_gr(ineig).gt.rhot) then
! We connect this group to its neighbourg
        igroup2=group(igroup)%isad_gr(ineig)
        igr2=igroupid(igroup2)
        if (color(igr2).eq.0) then
!!$           call do_colorize(icolor_select,igroup2,igr2,rhot)
           call do_colorize(igroup2,igr2,rhot) !YDdebug
        elseif (color(igr2).ne.icolor_select) then
           write(errunit,*) 'ERROR in do_colorize : color(igr2) <> icolor_select'
           write(errunit,*) 'The connections are not symmetric.'
           STOP
        endif
     else
! We do not need this saddle anymore (and use the fact that saddles 
! are ranked in decreasing order)
        neig=neig-1
     endif
  enddo
  group(igroup)%nhnei=neig

end subroutine do_colorize

!=======================================================================
subroutine compute_saddle_list
!=======================================================================
! Compute the lowest density threshold below which each group is 
! connected to an other one
!=======================================================================

  implicit none
  integer(kind=4) :: ipar1,ipar2,igroup2,ineig,idist,igroup1,ineig2
  integer(kind=4) :: neig,ineigal,in1,in2,idestroy,icon_count
  integer(kind=4) :: i
  real(kind=8)  :: density1,density2,rho_sad12
  logical :: exist
  logical, allocatable :: touch(:)
  integer, allocatable :: listg(:)
  real(kind=8),  allocatable :: rho_sad(:)  
  integer, allocatable :: isad(:)
  integer, allocatable :: indx(:)

  if (verbose) write(errunit,*) 'Fill the end of the branches of the group tree'

! Allocate the array of nodes
  allocate(group(ngroups))

  if (verbose) write(errunit,*) 'First count the number of neighbourgs'// &
&                         ' of each elementary group...'

  allocate(touch(ngroups))
  allocate(listg(ngroups))

  touch(1:ngroups)=.false.


! First count the number of neighbourgs for each group to allocate
! arrays isad_in,isad_out,rho_saddle_gr
  do igroup1=1,ngroups
     ineig=0
     ipar1=firstpart(igroup1)
! Loop on all the members of the group
     do while (ipar1.gt.0)
        do idist=1,nhop
           ipar2=iparneigh(idist,ipar1)
           igroup2=igrouppart(ipar2)
! we test that we are in a group (i.e. that density(ipar) >= rho_threshold)
! and that this group is different from the one we are sitting on
           if (igroup2.gt.0.and.igroup2.ne.igroup1) then
              if (.not.touch(igroup2)) then
                 ineig=ineig+1
                 touch(igroup2)=.true.
                 listg(ineig)=igroup2
              endif
           endif
        enddo
! Next member
        ipar1=idpart(ipar1)
     enddo
! Reinitialize touch
     do in1=1,ineig
        igroup2=listg(in1)
        touch(igroup2)=.false.
     enddo
! Allocate the nodes 
     group(igroup1)%nhnei=ineig
     ineigal=max(ineig,1)
     allocate(group(igroup1)%isad_gr(1:ineigal))
     allocate(group(igroup1)%rho_saddle_gr(1:ineigal))
  enddo 

  if (verbose) write(errunit,*) 'Compute lists of neighbourgs and saddle points...'


! arrays isad_in,isad_out,rho_saddle_gr
  do igroup1=1,ngroups
! No calculation necessary if no neighbourg
     neig=group(igroup1)%nhnei
     if (neig.gt.0) then
        ineig=0
        ipar1=firstpart(igroup1)
        allocate(rho_sad(1:neig))
! Loop on all the members of the group
        do while (ipar1.gt.0)
           density1=density(ipar1)
           do idist=1,nhop
              ipar2=iparneigh(idist,ipar1)
              igroup2=igrouppart(ipar2)
! we test that we are in a group (i.e. that density(ipar) >= rho_threshold)
! and that this group is different from the one we are sitting on
              if (igroup2.gt.0.and.igroup2.ne.igroup1) then
                 density2=density(ipar2)
                 if (.not.touch(igroup2)) then
                    ineig=ineig+1
                    touch(igroup2)=.true.
                    listg(igroup2)=ineig
                    rho_sad12=min(density1,density2)
                    rho_sad(ineig)=rho_sad12
                    group(igroup1)%isad_gr(ineig)=igroup2
                 else
                    ineig2=listg(igroup2)
                    rho_sad12=min(density1,density2)
                    rho_sad(ineig2)=max(rho_sad(ineig2),rho_sad12)
                 endif
              endif
           enddo
! Next member
           ipar1=idpart(ipar1)
        enddo
        if (ineig.ne.neig) then
! Consistency checking
           write(errunit,*) 'ERROR in compute_saddle_list :'
           write(errunit,*) 'The number of neighbourgs does not match.'
           write(errunit,*) 'ineig, neig =',ineig,neig
           STOP
        endif
        group(igroup1)%rho_saddle_gr(1:ineig)=rho_sad(1:ineig)
        deallocate(rho_sad)
! Reinitialize touch
        do in1=1,ineig
           igroup2=group(igroup1)%isad_gr(in1)
           touch(igroup2)=.false.
        enddo
! No neighbourg
     endif
  enddo 

  deallocate(touch)
  deallocate(listg)

  if (verbose) write(errunit,*) 'Establish symmetry in connections...'

! Total number of connections count
  icon_count=0

! Destroy the connections between 2 groups which are not symmetric
! This might be rather slow and might be discarded later
  idestroy=0
  do igroup1=1,ngroups
     if (group(igroup1)%nhnei.gt.0) then
        do in1=1,group(igroup1)%nhnei
           exist=.false.
           igroup2=group(igroup1)%isad_gr(in1)
           if (igroup2.gt.0) then
              do in2=1,group(igroup2)%nhnei
                 if (group(igroup2)%isad_gr(in2).eq.igroup1) then
                    exist=.true.
                    rho_sad12=min(group(igroup2)%rho_saddle_gr(in2), &
 &                                group(igroup1)%rho_saddle_gr(in1))
                    group(igroup2)%rho_saddle_gr(in2)=rho_sad12
                    group(igroup1)%rho_saddle_gr(in1)=rho_sad12
                 endif
              enddo
           endif
           if (.not.exist) then
              group(igroup1)%isad_gr(in1)=0
              idestroy=idestroy+1
           else
              icon_count=icon_count+1
           endif
        enddo
     endif
  enddo

  if (verbose) write(errunit,*) 'Number of connections removed :',idestroy
  if (verbose) write(errunit,*) 'Total number of connections remaining :',icon_count

  if (verbose) then
     write(errunit,*) 'Rebuild groups with undesired connections removed...'
  endif


! Rebuild the group list correspondingly with the connections removed
! And sort the list of saddle points
  do igroup1=1,ngroups
     neig=group(igroup1)%nhnei
     if (neig.gt.0) then
        allocate(rho_sad(neig))
        allocate(isad(neig))
        ineig=0
        do in1=1,neig
           igroup2=group(igroup1)%isad_gr(in1)
           if (igroup2.gt.0) then
              ineig=ineig+1
              rho_sad(ineig)=group(igroup1)%rho_saddle_gr(in1)
              isad(ineig)=igroup2
           endif
        enddo
        deallocate(group(igroup1)%isad_gr)
        deallocate(group(igroup1)%rho_saddle_gr)
        ineigal=max(ineig,1)
        allocate(group(igroup1)%isad_gr(ineigal))
        allocate(group(igroup1)%rho_saddle_gr(ineigal))
        group(igroup1)%nhnei=ineig
        if (ineig.gt.0) then
! sort the saddle points by decreasing order
           allocate(indx(ineig))
           call indexx(ineig,rho_sad,indx)
           do i=1,ineig
              ineig2=indx(i)
              group(igroup1)%isad_gr(ineig-i+1)=isad(ineig2)
              group(igroup1)%rho_saddle_gr(ineig-i+1)=rho_sad(ineig2)
           enddo
           deallocate(indx)
        endif
        deallocate(rho_sad)
        deallocate(isad)
     endif
  enddo

  deallocate(iparneigh)
  deallocate(igrouppart)

end subroutine compute_saddle_list

!=======================================================================
subroutine compute_density(ipar,dist2,iparnei)
!=======================================================================

  implicit none
  real(kind=8)          :: dist2(0:nvoisins)
  integer(kind=4)       :: iparnei(nvoisins)
  real(kind=8)          :: r,unsr,contrib
  real(kind=8),external :: spline
  integer(kind=4)       :: idist,ipar

  r=sqrt(dist2(nvoisins))*0.5
  unsr=1./r
  contrib=0.
  do idist=1,nvoisins-1
     if(allocated(mass)) then
        contrib=contrib+mass(iparnei(idist))*spline(sqrt(dist2(idist))*unsr)
     else
        contrib=contrib+massp*spline(sqrt(dist2(idist))*unsr)
     end if
  enddo
! Add the contribution of the particle itself and normalize properly
! to get a density with average unity (if computed on an uniform grid)
! note that this assumes that the total mass in the box is normalized to 1.
  if(allocated(mass)) then
     density(ipar)=(xlong*ylong*zlong)*(contrib+mass(ipar)) &
          &              /(pi*r**3)
  else
     density(ipar)=(xlong*ylong*zlong)*(contrib+massp) &
          &              /(pi*r**3)
  end if

end subroutine compute_density

!=======================================================================
subroutine find_nearest_parts(ipar,dist2,iparnei)
!=======================================================================

  implicit none

  integer(kind=4) :: ipar,idist,icell_identity,inccellpart
  real(kind=8)    :: dist2(0:nvoisins)
  integer(kind=4) :: iparnei(nvoisins)
  real(kind=8)    :: poshere(1:3)

  poshere(1:3)=pos(ipar,1:3)
  dist2(0)=0.
  do idist=1,nvoisins
     dist2(idist)=bignum
  enddo
  icell_identity =1
  inccellpart    =0
  call walk_tree(icell_identity,poshere,dist2,ipar,inccellpart,iparnei)

end subroutine find_nearest_parts

!=======================================================================
recursive subroutine walk_tree(icellidin,poshere,dist2,    &
 &                   iparid,inccellpart,iparnei)
!=======================================================================
  implicit none


  integer(kind=4) :: icellidin,icell_identity,iparid,inccellpart,ic,iparcell
  real(kind=8)    :: poshere(3),dist2(0:nvoisins)
  real(kind=8)    :: dx,dy,dz,distance2,sc
  integer(kind=4) :: idist,inc
  integer(kind=4) :: icellid_out
  real(kind=8)    :: discell2(0:8)
  integer(kind=4) :: iparnei(nvoisins)
  integer(kind=4) :: icid(8)

  integer(kind=4) :: i,npart_pos_this_node
  real(kind=8)    :: distance2p

  icell_identity=firstchild(icellidin)
  inc=1
  discell2(0)=0
  discell2(1:8)=1.e30
  do while (icell_identity.ne.0)
     sc=size_cell(icell_identity)
     dx=abs(pos_cell(1,icell_identity)-poshere(1))
     dx=max(0.,min(dx,real(xlong,8)-dx)-sc)
     dy=abs(pos_cell(2,icell_identity)-poshere(2))
     dy=max(0.,min(dy,real(ylong,8)-dy)-sc)
     dz=abs(pos_cell(3,icell_identity)-poshere(3))
     dz=max(0.,min(dz,real(zlong,8)-dz)-sc)
     distance2=dx**2+dy**2+dz**2
     if (distance2.lt.dist2(nvoisins)) then
        idist=inc-1
        do while (discell2(idist).gt.distance2)
           discell2(idist+1)=discell2(idist)
           icid(idist+1)=icid(idist)
           idist=idist-1
        enddo
        discell2(idist+1)=distance2
        icid(idist+1)=icell_identity
        inc=inc+1
     endif
     icell_identity=sister(icell_identity)
  enddo
  inccellpart=inccellpart+inc-1
  do ic=1,inc-1
     icellid_out=icid(ic)
     if (firstchild(icellid_out) < 0) then
        if (discell2(ic).lt.dist2(nvoisins)) then
           npart_pos_this_node=-firstchild(icellid_out)-1
           do i=npart_pos_this_node+1,npart_pos_this_node+mass_cell(icellid_out)
              iparcell=idpart(i)
              dx=abs(pos(iparcell,1)-poshere(1))
              dx=max(0.,min(dx,real(xlong,8)-dx))
              dy=abs(pos(iparcell,2)-poshere(2))
              dy=max(0.,min(dy,real(ylong,8)-dy))
              dz=abs(pos(iparcell,3)-poshere(3))
              dz=max(0.,min(dz,real(zlong,8)-dz))
              distance2p=dx**2+dy**2+dz**2
              if (distance2p .lt. dist2(nvoisins)) then 
                 if (iparcell.ne.iparid) then
                    idist=nvoisins-1
                    do while (dist2(idist).gt.distance2p)
                       dist2(idist+1)=dist2(idist)
                       iparnei(idist+1)=iparnei(idist)
                       idist=idist-1
                    enddo
                    dist2(idist+1)=distance2p
                    iparnei(idist+1)=iparcell
                 endif
              endif
           enddo
        endif              
     elseif (discell2(ic).lt.dist2(nvoisins)) then
        call walk_tree(icellid_out,poshere,dist2,iparid,inccellpart,iparnei)
     endif
  enddo
end subroutine walk_tree

!=======================================================================
subroutine create_tree_structure
!=======================================================================

  implicit none
  integer(kind=4) :: nlevel,inccell,idmother,ipar
  integer(kind=4) :: npart_this_node,npart_pos_this_node
  integer(kind=4) :: ncell
  real(kind=8)    :: pos_this_node(3)

  if (verbose) write(errunit,*) 'Create tree structure...'

  ! we modified to put 2*npart-1 instead of 2*npart so that AdaptaHOP can work on a 1024^3, 2*(1024^3)-1 is still an integer(kind=4), 2*(1024^3) is not 
  ncellmx=2*npart -1
  ncellbuffer=max(nint(0.1*npart),ncellbuffermin)
  allocate(idpart(npart))
  allocate(idpart_tmp(npart))
  allocate(mass_cell(ncellmx))
  allocate(size_cell(ncellmx))
  allocate(pos_cell(3,ncellmx))
  allocate(sister(1:ncellmx))
  allocate(firstchild(1:ncellmx))
  do ipar=1,npart
     idpart(ipar)=ipar
  enddo
  nlevel=0
  inccell=0
  idmother=0
  pos_this_node(1:3)=0.
  npart_this_node=npart
  npart_pos_this_node=0
  idpart_tmp(1:npart)=0
  pos_cell(1:3,1:ncellmx)=0.
  size_cell(1:ncellmx)=0.
  mass_cell(1:ncellmx)=0
  sister(1:ncellmx)=0
  firstchild(1:ncellmx)=0
  sizeroot=real(max(xlong,ylong,zlong))
    
  call create_KDtree(nlevel,pos_this_node, &
 &                   npart_this_node,npart_pos_this_node,inccell, &
 &                   idmother)
  ncell=inccell

  if (verbose) write(errunit,*) 'total number of cells =',ncell

  deallocate(idpart_tmp)

end subroutine create_tree_structure


!=======================================================================
recursive subroutine create_KDtree(nlevel,pos_this_node,npart_this_node, &
 &                npart_pos_this_node,inccell,idmother)
!=======================================================================
!  nlevel : level of the node in the octree. Level zero corresponds to 
!           the full box
!  pos_this_node : position of the center of this node
!  npart  : total number of particles 
!  idpart : array of dimension npart containing the id of each
!           particle. It is sorted such that neighboring particles in
!           this array belong to the same cell node.
!  idpart_tmp : temporary array of same size used as a buffer to sort
!           idpart.
!  npart_this_node : number of particles in the considered node
!  npart_pos_this_node : first position in idpart of the particles 
!           belonging to this node
!  pos :    array of dimension 3.npart giving the positions of 
!           each particle belonging to the halo
!  inccell : cell id number for the newly created structured grid site
!  pos_cell : array of dimension 3.npart (at most) corresponding
!           to the positions of each node of the structured grid
!  mass_cell : array of dimension npart (at most) giving the 
!           number of particles in each node of the structured grid
!  size_cell : array of dimension npart (at most) giving half the
!           size of the cube forming each node of the structured grid
!  sizeroot : size of the root cell (nlevel=0)
!  idmother : id of the mother cell
!  sister   : sister of a cell (at the same level, with the same mother)
!  firstchild : first child of a cell (then the next ones are found 
!             with the array sister). If it is a cell containing only
!           one particle, it gives the id of the particle.
!  ncellmx : maximum number of cells
!  megaverbose : detailed verbose mode
!=======================================================================

   implicit none

   integer(kind=4)           :: nlevel,npart_pos_this_node,npart_this_node
   real(kind=8)              :: pos_ref(3,0:7)
   real(kind=8)              :: pos_this_node(3)   
   integer(kind=4)           :: ipar,icid,j,inccell,nlevel_out
   integer(kind=4), external :: icellid
   integer(kind=4)           :: npart_pos_this_node_out,npart_this_node_out
   integer(kind=4)           :: incsubcell(0:7),nsubcell(0:7)
   real(kind=8)              :: xtest(3),pos_this_node_out(3)
   integer(kind=4)           :: idmother,idmother_out

   integer(kind=8)           :: ncellmx_old
   integer, allocatable      :: mass_cell_tmp(:),sister_tmp(:),firstchild_tmp(:)
   real(kind=8), allocatable :: size_cell_tmp(:),pos_cell_tmp(:,:)

!  pos_ref : an array used to find positions of the 8 subcells in this
!           node.
   data  pos_ref /-1.,-1.,-1., &
 &                 1.,-1.,-1., &
 &                -1., 1.,-1., &
 &                 1., 1.,-1., &
 &                -1.,-1., 1., &
 &                 1.,-1., 1., &
 &                -1., 1., 1., &
 &                 1., 1., 1.  / 


   if (npart_this_node.gt.0) then
      inccell=inccell+1
      if (mod(inccell,1000000).eq.0.and.megaverbose) write(errunit,*) 'inccell=',inccell
      if (inccell.gt.ncellmx) then
         ncellmx_old=ncellmx
         ncellmx=ncellmx+ncellbuffer
         if (megaverbose) write(errunit,*) &
 &          'ncellmx is too small. Increase it and reallocate arrays accordingly'
         allocate(mass_cell_tmp(ncellmx_old))
         mass_cell_tmp(1:ncellmx_old)=mass_cell(1:ncellmx_old)
         deallocate(mass_cell)
         allocate(mass_cell(ncellmx))
         mass_cell(1:ncellmx_old)=mass_cell_tmp(1:ncellmx_old)
         deallocate(mass_cell_tmp)
         allocate(sister_tmp(ncellmx_old))
         sister_tmp(1:ncellmx_old)=sister(1:ncellmx_old)
         deallocate(sister)
         allocate(sister(ncellmx))
         sister(1:ncellmx_old)=sister_tmp(1:ncellmx_old)
         deallocate(sister_tmp)
         allocate(firstchild_tmp(ncellmx_old))
         firstchild_tmp(1:ncellmx_old)=firstchild(1:ncellmx_old)
         deallocate(firstchild)
         allocate(firstchild(ncellmx))
         firstchild(1:ncellmx_old)=firstchild_tmp(1:ncellmx_old)
         firstchild(ncellmx_old:ncellmx)=0
         deallocate(firstchild_tmp)
         allocate(size_cell_tmp(ncellmx_old))
         size_cell_tmp(1:ncellmx_old)=size_cell(1:ncellmx_old)
         deallocate(size_cell)
         allocate(size_cell(ncellmx))
         size_cell(1:ncellmx_old)=size_cell_tmp(1:ncellmx_old)
         deallocate(size_cell_tmp)
         allocate(pos_cell_tmp(3,ncellmx_old))
         pos_cell_tmp(1:3,1:ncellmx_old)=pos_cell(1:3,1:ncellmx_old)
         deallocate(pos_cell)
         allocate(pos_cell(3,ncellmx))
         pos_cell(1:3,1:ncellmx_old)=pos_cell_tmp(1:3,1:ncellmx_old)
         deallocate(pos_cell_tmp)
      endif
      pos_cell(1:3,inccell)=pos_this_node(1:3)
      mass_cell(inccell)=npart_this_node
      size_cell(inccell)=2.**(-nlevel)*sizeroot*0.5
      if (idmother.gt.0) then
         sister(inccell)=firstchild(idmother)
         firstchild(idmother)=inccell
      endif
!     If there is only one particle in the node or we have reach
!     maximum level of refinement, we are done
      if ((npart_this_node <= npartpercell).or.(nlevel.eq.nlevelmax)) then
         firstchild(inccell)=-(npart_pos_this_node+1)
         return
      endif
   else
      return
   endif


!  Count the number of particles in each subcell of this node
   incsubcell(0:7)=0
   do ipar=npart_pos_this_node+1,npart_pos_this_node+npart_this_node
      xtest(1:3)=pos(idpart(ipar),1:3)-pos_this_node(1:3)
      icid=icellid(xtest)
      incsubcell(icid)=incsubcell(icid)+1
   enddo

!  Create the array of positions of the first particle of the lists
!  of particles belonging to each subnode
   nsubcell(0)=0
   do j=1,7
      nsubcell(j)=nsubcell(j-1)+incsubcell(j-1)
   enddo

!  Sort the array of ids (idpart) to gather the particles belonging
!  to the same subnode. Put the result in idpart_tmp.
   incsubcell(0:7)=0
   do ipar=npart_pos_this_node+1,npart_pos_this_node+npart_this_node
      xtest(1:3)=pos(idpart(ipar),1:3)-pos_this_node(1:3)
      icid=icellid(xtest)
      incsubcell(icid)=incsubcell(icid)+1
      idpart_tmp(incsubcell(icid)+nsubcell(icid)+npart_pos_this_node)=idpart(ipar)
   enddo
   
!  Put back the sorted ids in idpart
   do ipar=npart_pos_this_node+1,npart_pos_this_node+npart_this_node
      idpart(ipar)=idpart_tmp(ipar)
   enddo

!  Call again the routine for the 8 subnodes:
!  Compute positions of subnodes, new level of refinement, 
!  positions in the array idpart corresponding to the subnodes,
!  and call for the treatment recursively.
   nlevel_out=nlevel+1
   idmother_out=inccell
   do j=0,7
      pos_this_node_out(1:3)=pos_this_node(1:3)+ &
 &                           sizeroot*pos_ref(1:3,j)*2.**(-nlevel-2)
      npart_pos_this_node_out=npart_pos_this_node+nsubcell(j)
      npart_this_node_out=incsubcell(j)
      call create_KDtree(nlevel_out,pos_this_node_out,npart_this_node_out, &
 &                npart_pos_this_node_out,inccell,idmother_out)
   enddo
end subroutine create_KDtree

!=======================================================================
subroutine remove_degenerate_particles
!=======================================================================

  implicit none
  real(kind=8), parameter :: accurac=1.e-6

  integer, allocatable, dimension(:) :: idgene
  real(kind=8),allocatable, dimension(:) :: tmp
  real(kind=8),allocatable, dimension(:,:) :: possav
  logical,allocatable,dimension(:) :: move
  real(kind=8) :: xref,tolerance,tolerance2,phi,costeta,teta,sinteta
  real(kind=8) :: ran2
  logical :: doneiter,done
  integer(kind=4) :: i,niter,j,ipar,jpar,idd,jdd,incdege
  integer(kind=4) :: idum

  if (megaverbose) write(errunit,*) 'Move randomly particles at exactly the same position'

  idum=-111
  tolerance=max(xlongs2,ylongs2,zlongs2)*accurac
  allocate(idgene(npart),tmp(npart),move(npart))
  tolerance2=2.0*tolerance
  allocate(possav(3,npart))
  do i=1,npart
     possav(1:3,i)=pos(i,1:3)
  enddo

  niter=1
  doneiter=.false.
  do while (.not.doneiter) 
     if (megaverbose) write(errunit,*) 'Iteration :',niter
     niter=niter+1
     doneiter=.false.
     do i=1,npart
        tmp(i)=pos(i,1)
        move(i)=.false.
     enddo
     call indexx(npart,tmp,idgene)

     done=.false.
     i=1
     do while (.not.done)
        xref=tmp(idgene(i))
        j=i+1
        if (j > npart) then
           j=1
           done=.true.
        endif
        do while (abs(xref-tmp(idgene(j))) < tolerance)
           j=j+1
           if (j > npart) then
              j=1
              done=.true.
           endif
        enddo
        do ipar=i,j-1
           idd=idgene(ipar)
           do jpar=ipar+1,j-1
              jdd=idgene(jpar)
              if (abs(pos(jdd,2)-pos(idd,2)) < tolerance .and. &
 &                abs(pos(jdd,3)-pos(idd,3)) < tolerance) then
                 move(idgene(jpar))=.true.
              endif
           enddo
        enddo
        i=j
     enddo
     incdege=0
     do i=1,npart
        if (move(i)) then
           incdege=incdege+1
        endif
     enddo
     if (megaverbose) write(errunit,*) &
 &                 'Found the following number of degeneracies :',incdege
     do i=1,npart
        if (move(i)) then
           phi=2.*pi*ran2(idum)
           costeta=1.0-2.0*ran2(idum)
           teta=acos(costeta)           
           sinteta=sin(teta)
           pos(i,1)=possav(1,i)+tolerance2*cos(phi)*sinteta
           pos(i,2)=possav(2,i)+tolerance2*sin(phi)*sinteta
           pos(i,3)=possav(3,i)+tolerance2*costeta
        endif
     enddo
     
     if (incdege==0) doneiter=.true.
  enddo

  deallocate(possav)
  deallocate(tmp,idgene,move)
  
end subroutine remove_degenerate_particles

!=======================================================================
subroutine convtoasc(number,sstring)
!=======================================================================
! To convert an integer(kind=4) smaller than 999999 to a 6 characters string
!=======================================================================
  implicit none
  integer(kind=4) :: number, istring, num, nums10, i
  character*6 :: sstring
  character*10,parameter :: nstring='0123456789'
      
  num=1000000
  nums10=num/10
  do i=1,6
     istring=1+mod(number,num)/nums10
     sstring(i:i)=nstring(istring:istring)
     num=num/10
     nums10=nums10/10
  enddo
end subroutine convtoasc

!=======================================================================
end module neiKDtree


!=======================================================================
function spline(x)
!=======================================================================
  implicit none
  real(kind=8) :: spline,x

  if (x.le.1.) then
     spline=1.-1.5*x**2+0.75*x**3
  elseif (x.le.2.) then
     spline=0.25*(2.-x)**3
  else
     spline=0.
  endif

end function spline

!=======================================================================
function HsurH0(z,omega0,omegaL,omegaR)
!=======================================================================
  implicit none
  real(kind=8) :: z,omega0,omegaL,omegaR,HsurH0

  HsurH0=sqrt(Omega0*(1.d0+z)**3+OmegaR*(1.d0+z)**2+OmegaL)
end function HsurH0

!=======================================================================
function icellid(xtest)
!=======================================================================
!  Compute cell id corresponding to the signs of coordinates of xtest
!  as follows :
!  (-,-,-) : 0
!  (+,-,-) : 1
!  (-,+,-) : 2
!  (+,+,-) : 3
!  (-,-,+) : 4
!  (+,-,+) : 5
!  (-,+,+) : 6
!  (+,+,+) : 7
!  For self-consistency, the array pos_ref should be defined exactly 
!  with the same conventions
!=======================================================================
   implicit none
   integer(kind=4) :: icellid,j,icellid3d(3)
   real(kind=8) :: xtest(3)
   do j=1,3
      if (xtest(j).ge.0) then
         icellid3d(j)=1
      else
         icellid3d(j)=0
      endif
   enddo
   icellid=icellid3d(1)+2*icellid3d(2)+4*icellid3d(3)
end function icellid
