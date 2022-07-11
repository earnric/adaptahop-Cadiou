module fof

  use halo_defs

  ! defs for creating fof groups
  ! octtree
  type cell
     real(kind=8)    :: size                    ! half lengh of the cell
     real(kind=8)    :: pos(3)                  ! position of the center of the cell
     integer(kind=4) :: sister
     integer(kind=4) :: istart                  ! name of the cell next to this one with the same mother
     integer(kind=4) :: npart                   ! number of particule per cell
     integer(kind=4) :: mother
     integer(kind=4) :: first_child 
  end type cell

  type(cell), allocatable :: kd_cell(:)
  type(cell)              :: celltmp(8)

  integer(kind=4) :: ncell, nleaf, ncellmax, npcellmax
  real(kind=8)    :: pos_sh(3,8), p0(6,3), sizemax, sizemin
  integer(kind=4), allocatable :: cellmother(:),cellsister(:),cellchild(:)
  real(kind=8),    allocatable :: cellpos(:,:), cellsize(:)
  integer(kind=4), allocatable :: iorder(:), iordertmp(:), idcellpart(:)

  ! fof group
  real(kind=8)    ::fEps, fEps2
  integer(kind=4), allocatable :: lstart(:),lend(:),ifof_order(:),idleaf(:)
  integer(kind=4) :: nfofmax
  integer(kind=4),allocatable :: fofstart(:), fofend(:)
  integer(kind=1),allocatable :: fofsorted(:) 

contains
  
!***********************************************************************************
! Subroutines for fof
!***********************************************************************************

!***********************************************************************************
subroutine fof_init
!***********************************************************************************

  implicit none

  ! current mean interparticular distance (total mass and length box in code units are 1. so that 
  ! massp is really 1./Nbodies here, provided all mass are equal, otherwise massp should be the mass 
  ! of the highest resolution particles in the multi-mass case) multiplied by linking length parameter
  ! b_init so that all particles separated by a distance smaller than fEps will be identified as part of 
  ! the same structure by the FOF algorithm ... 


  real(kind=8) :: mass_eq

  if(massp.le.0.0) stop
! if simulation has been renormalized everything is fine, otherwise 
! we have to do it here so that fEps is truly a fraction of the mean 
! interparticular distance (i.e. mass_eq = 1/Nbodies) which is not the 
! case if omega_baryons is not zero!
  if (allocated(mass)) then 
     mass_eq = sum(real(mass,8))
     mass_eq = real(massp,8)/mass_eq
  else
     mass_eq = 1d0/real(nbodies,8)
  endif
  fEps      = b_init*real(mass_eq,8)**(1./3.)
  fEps2     = fEps*fEps
  sizemin   = fEps/sqrt(3.0)
  npcellmax = 50
  ncellmax  = ceiling(real(nbodies)/5.0)
  ! We can safely assume that if we only fof groups on nMember particules only one 1 particle in 5 would be in such a halo have 
  nfofmax   = ceiling(real(nbodies)/real(5*nMembers))
  sizemax     = 0.5
  
  pos_sh(1,1) = - 1
  pos_sh(2,1) = - 1
  pos_sh(3,1) = - 1

  pos_sh(1,2) = + 1
  pos_sh(2,2) = - 1
  pos_sh(3,2) = - 1
  
  pos_sh(1,3) = + 1
  pos_sh(2,3) = + 1
  pos_sh(3,3) = - 1

  pos_sh(1,4) = + 1
  pos_sh(2,4) = + 1
  pos_sh(3,4) = + 1

  pos_sh(1,5) = + 1
  pos_sh(2,5) = - 1
  pos_sh(3,5) = + 1

  pos_sh(1,6) = - 1
  pos_sh(2,6) = - 1
  pos_sh(3,6) = + 1

  pos_sh(1,7) = - 1
  pos_sh(2,7) = + 1
  pos_sh(3,7) = + 1

  pos_sh(1,8) = - 1
  pos_sh(2,8) = + 1
  pos_sh(3,8) = - 1

  write(errunit,*) '============================================================'
  write(errunit,*) ' FOF was called with the following attributes : '
  write(errunit,*) ' Verbose                                      : ',verbose
  write(errunit,*) ' Linking length parameter @ z=0               : ', b_init
  write(errunit,*) ' Linking length fEps (simu units)             : ', fEps
  write(errunit,*) ' Box size (simu units)                        : ', sizemax
  if(verbose) then
     write(errunit,*) ' Number max of cells                          : ', ncellmax 
     write(errunit,*) ' Minimal size of cell                         : ', sizemin
     write(errunit,*) ' Or number max of particles                   : ', npcellmax
     write(errunit,*) ' Maximum number of haloes                     : ',nfofmax
  end if
  write(errunit,*) ' Minimum number of parts per halo             : ', nMembers         
  write(errunit,*) '============================================================'

  return

end subroutine fof_init

!***********************************************************************************
subroutine fof_main
!***********************************************************************************
! here it is the main subroutine
! all position must be contained in a box of size 1.0, and center (0,0,0) 

  implicit none

  if(verbose) write(errunit,*) 
  call create_tree
  call make_fof_groups
  if(verbose) write(errunit,*)
  return

end subroutine fof_main

!***********************************************************************************
subroutine create_tree
!***********************************************************************************
! This routine build an octtree and stops when the last cell is smaller then sizemin
!***********************************************************************************

  implicit none
  
  integer(kind=4)              :: icell,ipar,nbodiesch,ileaf,imoth,icellnew,i,npmaxtmp  
  integer(kind=4), allocatable :: idcell(:)

  if(verbose) write(errunit,*) '> Creating tree'

  allocate(iorder(nbodies), iordertmp(nbodies), idcellpart(nbodies))
  do ipar = 1, nbodies
     iorder(ipar) = ipar
  end do

  allocate(kd_cell(ncellmax))

  ! Create root node  
  kd_cell(1)%pos(1:3)  = 0
  kd_cell(1)%size      = sizemax
  kd_cell(1)%sister    = 1
  kd_cell(1)%npart     = nbodies
  kd_cell(1)%istart    = 1
  kd_cell(1)%mother    = 0
  ncell                = 1

  ! Creating tree from root node
  icell = 0
  do while(icell.ne.1)
     if(icell==0) icell=1
 
     do while((kd_cell(icell)%npart.gt.npcellmax).and.(kd_cell(icell)%size.gt.sizemin))
        ! cutting cell in 8 cells until it is smaller than sizemin
        call create_cells(icell)
        if( kd_cell(icell)%first_child.ne.0) icell = kd_cell(icell)%first_child
     end do

     ! Search for a cell with more than npcell particules,or greater than sizemin
     do while((kd_cell(icell)%npart.le.npcellmax).or.(kd_cell(icell)%size.le.sizemin))
        icell = kd_cell(icell)%sister
     end do
 
     if(icell==0) stop 'Error occured in routine create tree, icell = 0'
  end do
  deallocate(iordertmp,idcellpart)

  if(verbose) write(errunit,*) '> Number of cells created  :', ncell

  nleaf     = 0
  nbodiesch = 0
  do icell = 1, ncell
     if(kd_cell(icell)%first_child.le.0) then
        if((kd_cell(icell)%size.gt.sizemin).and.(kd_cell(icell)%npart.gt.npcellmax)) stop 'error in leaf size'
        nleaf     = nleaf + 1
        nbodiesch = nbodiesch + kd_cell(icell)%npart
     end if
  end do

  if(nbodiesch.ne.nbodies) then
     write(errunit,*) nbodies, nbodiesch
     stop 'all particles are not contained in leafs'
  end if
  
  if(verbose) write(errunit,*) '> Number of leaves created :', nleaf

  allocate(idcell(0:ncell))
  idcell    = -1
  idcell(0) = 0
  imoth     = nleaf + 1
  idcell(1) = imoth
  icell     = 2   
  ileaf     = 0
  do while(icell.gt.1)
     if(idcell(icell).gt.0) stop 'error in cell reordering'
     if(kd_cell(icell)%first_child.gt.0) then
        imoth         = imoth + 1
        idcell(icell) = imoth
        icell         = kd_cell(icell)%first_child
     else
        ileaf         = ileaf + 1
        idcell(icell) = ileaf
        icell         = kd_cell(icell)%sister
     end if
  end do

  allocate(lstart(nleaf),lend(nleaf))
  allocate(cellpos(ncell,3),cellsize(ncell),cellmother(ncell),cellsister(ncell),cellchild(ncell))

  do icell = 1, ncell
     if(idcell(icell).le.0) then
        write(*,*) 'idcell is nil here', icell
        stop
     end if
     icellnew = idcell(icell) 
     cellmother(icellnew)      = idcell(kd_cell(icell)%mother)
     cellsister(icellnew)      = idcell(kd_cell(icell)%sister)
     cellchild(icellnew)       = idcell(kd_cell(icell)%first_child)
     cellpos(icellnew,1:3)     = kd_cell(icell)%pos(1:3)
     cellsize(icellnew)        = kd_cell(icell)%size
     if(cellchild(icellnew).le.0) then
        if(icellnew.gt.nleaf) stop 'Error for lstart and lend'
        lstart(icellnew)          = kd_cell(icell)%istart
        lend(icellnew)            = kd_cell(icell)%npart + kd_cell(icell)%istart - 1
     end if
  end do
  
  deallocate(kd_cell,idcell)
  
  if(verbose) then
     write(errunit,*) '> Number min of particles in leaf:', minval(lend(1:nleaf) - lstart(1:nleaf)) + 1
     write(errunit,*) '> Number max of particles in leaf:', maxval(lend(1:nleaf) - lstart(1:nleaf)) + 1
  end if
  
  ! modify Kdtree to stop when ast sister is reached
  cellsister(nleaf+1) = 0
  do icell = 1,ncell
     if(cellsister(icell).gt.0) then
        if(cellmother(cellsister(icell)).ne.cellmother(icell)) cellsister(icell) = 0
     end if
  end do

  allocate(idleaf(nbodies))
  idleaf   = 0
  npmaxtmp = 0
  do ileaf = 1,nleaf
     do i = lstart(ileaf), lend(ileaf)
        ipar         = iorder(i)
        idleaf(ipar) = ileaf
     end do
     if(cellsize(ileaf).ge.sizemin) npmaxtmp = max(npmaxtmp,(lend(ileaf) - lstart(ileaf) + 1))
  end do
  if(verbose) write(errunit,*) '> Number max of parts in big leaf:', npmaxtmp

  return
 
end subroutine create_tree

!***********************************************************************************
subroutine create_cells(icell)
!***********************************************************************************

  implicit none

  integer(kind=4):: icell

  integer(kind=4):: istart, iend, pos_inverse, ncreated
  integer(kind=4):: i, j, ip,itmp, icelltmp
  real(kind=8)   :: d(3), size

  istart = kd_cell(icell)%istart
  iend   = istart + kd_cell(icell)%npart - 1
  size   = kd_cell(icell)%size / 2
  if(istart ==0) stop 'error istart = 0'
  
  do icelltmp = 1, 8
     celltmp(icelltmp)%size   = size
     do j = 1, 3
        celltmp(icelltmp)%pos(j) = kd_cell(icell)%pos(j) + pos_sh(j,icelltmp) * size
     end do
     celltmp(icelltmp)%istart = 0
     celltmp(icelltmp)%npart  = 0
  end do

  do i = 0, kd_cell(icell)%npart -1
     ip = iorder( i + istart)
     do j = 1,3
        d(j) = pos(ip,j) - kd_cell(icell)%pos(j)
     end do
      idcellpart(ip) = pos_inverse(d)
  end do

  ncreated = 0
  itmp = kd_cell(icell)%istart - 1
  do icelltmp = 1,8
     do i = 0, kd_cell(icell)%npart - 1
        ip = iorder(istart + i)
        if(idcellpart(ip) == icelltmp) then
           itmp = itmp + 1
           iordertmp(itmp) = iorder(i+istart)
           if(celltmp(icelltmp)%npart==0) then
              celltmp(icelltmp)%istart = itmp 
              ncreated = ncreated + 1
           end if
           celltmp(icelltmp)%npart = celltmp(icelltmp)%npart + 1
        end if
     end do
  end do
  if(itmp.ne.iend) stop 'Error in cells creation'
 
  iorder(istart:iend) = iordertmp(istart:iend)
 
  if(ncreated.gt.1) then
     kd_cell(icell)%first_child = ncell + 1
     do icelltmp = 1, 8
        if(celltmp(icelltmp)%npart.ne.0) then
           ncell                   = ncell + 1
           if(ncell>ncellmax) then
              write(errunit, *) ncell, ncellmax
              stop 'To much cells to create'
           end if
           kd_cell(ncell)        = celltmp(icelltmp)
           kd_cell(ncell)%mother = icell
           kd_cell(ncell)%sister = ncell + 1
        end if
     end do
     kd_cell(ncell)%sister = kd_cell(icell)%sister
   else
      ! no new cell are created change pos and size of cell please
      kd_cell(icell)%size = size
      do icelltmp = 1, 8
         if(celltmp(icelltmp)%npart.ne.0) kd_cell(icell)%pos(1:3) = celltmp(icelltmp)%pos(1:3)
      end do
  end if

  return

end subroutine create_cells

!***********************************************************************************
subroutine make_fof_groups
!***********************************************************************************

  implicit none
  integer(kind=4) :: ip,i,ilast_index,ifirst_index,ig_fof,nfof_groups,iprint,fprint

  if(verbose) write(errunit,*) '> Making fof groups'
 
  allocate(fofstart(nfofmax),fofend(nfofmax))
  fofstart     = -1
  fofend       = -1
  allocate(fofsorted(nbodies),ifof_order(nbodies))
  fofsorted    = 0
  ifof_order   = -1

  ilast_index  = 0
  ifirst_index = 0
  nfof_groups  = 0
  iprint       = 10
  fprint       = ceiling(real(nbodies)/10.)
  if(verbose) write(errunit,'(a,i3,a)') ' > done ',0,' %'
  do i = 1,nbodies
     if(verbose) then
        if(mod(i,fprint) .eq. 0) then
           write(errunit,'(a,i3,a)') ' > done ',iprint,' %'
           iprint = iprint+10
        end if
     end if
     ip = iorder(i)   ! treat particles as ordered in leaves
     if(fofsorted(ip).eq.0) then    
        call new_fof_group(ip,ifirst_index,ilast_index)
        if(ilast_index-ifirst_index +1 .ge.nMembers) then
           nfof_groups = nfof_groups + 1
           if(nfof_groups.gt.nfofmax) stop 'Found too many halos' 
           fofstart(nfof_groups) = ifirst_index
           fofend(nfof_groups)   = ilast_index
        end if
     end if
  end do
  if(verbose) write(errunit,'(a,i3,a)') ' > done ',100,' %'

  if(ilast_index.ne.nbodies) then
     stop 'all particles have not been sorted into fof groups'
  end if

  deallocate(iorder,idleaf)
  deallocate(cellmother,cellsister,cellchild,cellpos,cellsize)
  deallocate(lstart,lend)
  deallocate(fofsorted)

  if(verbose) write(errunit,*) '> Number of initial fof groups:',nfof_groups
  
  allocate(liste_parts(nbodies)) ! allocate liste_parts at the latest possible moment
  liste_parts(1:nbodies) = 0
  nb_of_halos            = 0

  do ig_fof = 1, nfof_groups
     nb_of_halos = nb_of_halos + 1
     do i = fofstart(ig_fof),fofend(ig_fof)
        ip = ifof_order(i)
        if(liste_parts(ip).gt.0) stop 'fatal error for liste_parts'
        liste_parts(ip) = nb_of_halos
     end do
  end do
  deallocate(fofstart,fofend)
  deallocate(ifof_order)

  if(verbose) write(errunit,*) '> Number of final fof groups:',nb_of_halos

  return

end subroutine make_fof_groups

!***********************************************************************************
subroutine new_fof_group(ig_fof,ifirst_index,ilast_index)
!***********************************************************************************

  implicit none
  integer(kind=4) :: ig_fof,ifirst_index,ilast_index  
  integer(kind=4) :: ifof_index
  integer(kind=4) :: ip,icell,ileaf
  logical(kind=4) :: fsearch

  ifirst_index = ilast_index + 1
  ifof_index   = ilast_index

  call add_particle_to_fof_group(ig_fof,idleaf(ig_fof),ilast_index)

  do while(ifof_index.ne.ilast_index)
     ifof_index   = ifof_index + 1
     ip           = ifof_order(ifof_index)
     ileaf        = idleaf(ip)
     call test_part(ip,ileaf,fsearch)
     if(fsearch) then
        icell = nleaf + 1
        do while(icell.gt.0)
           call find_next_leaf(ip,icell)
           if(icell.gt.0) then
              ! search in the leaf for particles to add in fof group
              call search_leaf(ip,ileaf,icell,ilast_index)
              ! find an other cell to search
              do while(cellsister(icell).le.0.and.cellmother(icell).gt.0)
                 icell = cellmother(icell)
              end do
              if(cellsister(icell).gt.0) then
                 icell = cellsister(icell)
              else
                 ! search ended
                 icell = -1
              end if
           end if
        end do
     else
        ! search own leaf only
        call search_leaf(ip,ileaf,ileaf,ilast_index)
     end if
  end do

  return

end subroutine new_fof_group

!***********************************************************************************
subroutine find_next_leaf(ip,icell)
!***********************************************************************************

  implicit none
  integer(kind=4) :: ip,icell
  logical         :: fopen

  ! test to see if this cell should be opened
  call test_cell(ip,icell,fopen)

  ! search until we find a leaf we can iopen
  do while(cellchild(icell).gt.0.or..not.fopen)

     ! can open cell search a leaf to open
     ! go up the tree
     do while(fopen.and.cellchild(icell).gt.0) 
        icell = cellchild(icell)
        call test_cell(ip,icell,fopen)
     end do

     ! cannot open cell search an other cell to open 
     do while(.not.fopen)
        ! go sideway in the tree
        do while(.not.fopen.and.cellsister(icell).gt.0)
           icell = cellsister(icell)
           call test_cell(ip,icell,fopen)
        end do
        
        if(.not.fopen) then
           ! go down the tree
           do while(.not.fopen)
              if(cellmother(icell).le.0) then
                 ! search is ended no leaf is to be found
                 icell = -1
                 return
              else
                 ! test whether we should open mother cell
                 if(cellsister(icell).le.0) then
                    ! all sister tested don't bother
                    fopen = .false.
                 else   
                    call test_cell(ip,cellmother(icell),fopen)
                 end if
              end if
              if(.not.fopen) icell = cellmother(icell)
           end do
           
           ! go one step sideway in the tree icell has already been tested
           if(cellsister(icell).gt.0) then
              icell = cellsister(icell)
              call test_cell(ip,icell,fopen)
           else
              write(errunit,*) '> side+1 failed :',icell,cellchild(icell),cellsister(icell),cellmother(icell),fopen
              stop ' > No cellsister here'
           end if
        end if

     end do

  end do

  return

end subroutine find_next_leaf

!***********************************************************************************
subroutine test_cell(ip,icell,fopen)
!***********************************************************************************

  implicit none

  integer(kind=4) :: ip, icell
  logical         :: fopen
  integer(kind=4) :: ic
  real(kind=8)    :: d

  fopen = .true.
  do ic = 1,3
     d = abs(cellpos(icell,ic) - pos(ip,ic))
     if ((fPeriod(ic).ne.0.0).and.(d.gt.0.5)) d = 1.0 - d
     if (d.gt.cellsize(icell)) then
        d = d - cellsize(icell)
        if(d.gt.fEps) then
           fopen = .false.
           return
        end if
     end if
  end do

  return

end subroutine test_cell

!***********************************************************************************
subroutine test_part(ip,ileaf,fsearch)
!***********************************************************************************
! test part in its own leaf to check whether it can have neighbors in other leaves
  implicit none

  integer(kind=4) :: ip, ileaf
  logical         :: fsearch
  integer(kind=4) :: ic
  real(kind=8)    :: d

  fsearch = .false.
  do ic = 1,3
     d = abs(cellpos(ileaf,ic) - pos(ip,ic))
     d = cellsize(ileaf) - d
     if(d.lt.0.) stop ' > particule outside cell'
     if(d.le.fEps) then
        fsearch = .true.
        return
     end if
  end do

  return

end subroutine test_part

!***********************************************************************************
subroutine search_leaf(ip,ileaf,ileaflink,ilast_index)
!***********************************************************************************

  implicit none
  integer(kind=4) :: ip,ileaf,ileaflink,ilast_index
  integer(kind=4) :: i,ic,iplink
  real(kind=8)    :: d(3), dist2
  logical         :: quicker

  quicker = (cellsize(ileaflink).lt.sizemin)
 
  if(ileaflink.eq.ileaf) then
     if(quicker) return
  end if
  if(quicker) then
     i      = lstart(ileaflink)
     iplink = iorder(i)
     if(fofsorted(iplink).eq.1) return     
  end if
  do i = lstart(ileaflink),lend(ileaflink)
     iplink = iorder(i)
     if(fofsorted(iplink).eq.0) then
        dist2 = 0.
        do ic = 1,3
           d(ic) = abs(pos(ip,ic) - pos(iplink,ic))
           if ((fPeriod(ic).ne.0.0).and.(d(ic).gt.0.5)) d(ic) = 1.0 - d(ic)
           dist2 = dist2 + d(ic)**2
        end do
        if (dist2.le.fEps2) then
           call add_particle_to_fof_group(iplink,ileaflink,ilast_index)
           if(quicker) return
        end if
     end if
  end do
  
  return

end subroutine search_leaf

!***********************************************************************************
subroutine add_particle_to_fof_group(iplink,ileaflink,ilast_index)
!***********************************************************************************
  
  implicit none
  integer(kind=4) :: iplink,ileaflink,ilast_index
  integer(kind=4) :: i,ipl

  if(fofsorted(iplink).eq.1) stop 'particle already belongs to a FOF group'

  ! test leaf
  if(cellsize(ileaflink).lt.sizemin) then
     ! leaf is small enough add all leaf
     do i = lstart(ileaflink),lend(ileaflink) 
        ipl = iorder(i)
        fofsorted(ipl)          = 1 
        ilast_index             = ilast_index + 1
        ifof_order(ilast_index) = ipl
     end do
  else
     fofsorted(iplink)       = 1 
     ilast_index             = ilast_index + 1
     ifof_order(ilast_index) = iplink
  end if

  return

end subroutine add_particle_to_fof_group

!***********************************************************************************

end module fof


!***********************************************************************************
function pos_inverse(r)
!***********************************************************************************
! 1:  - - -
! 2:  + - -
! 3:  + + -
! 4:  + + +
! 5:  + - +
! 6:  - - +
! 7:  - + +
! 8:  - + - 
!***********************************************************************************
 
  implicit none

  real(kind=8)    :: r(3)
  integer(kind=4) :: pos_inverse

  if (r(1).lt.0.and.r(2).lt.0.and.r(3).lt.0) pos_inverse = 1
  if (r(1).ge.0.and.r(2).lt.0.and.r(3).lt.0) pos_inverse = 2
  if (r(1).ge.0.and.r(2).ge.0.and.r(3).lt.0) pos_inverse = 3 
  if (r(1).ge.0.and.r(2).ge.0.and.r(3).ge.0) pos_inverse = 4
  if (r(1).ge.0.and.r(2).lt.0.and.r(3).ge.0) pos_inverse = 5
  if (r(1).lt.0.and.r(2).lt.0.and.r(3).ge.0) pos_inverse = 6
  if (r(1).lt.0.and.r(2).ge.0.and.r(3).ge.0) pos_inverse = 7
  if (r(1).lt.0.and.r(2).ge.0.and.r(3).lt.0) pos_inverse = 8

end function pos_inverse
