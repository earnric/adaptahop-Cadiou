module subbox

  !JB: module to deal with splitting the simulation into pieces for low-memory analysis ... 

  use utils

  public 

  logical(kind=4)             :: subboxFoF
  real(kind=4),parameter      :: rel_prec = 1.e-6
  real(kind=4),parameter      :: maxnpart = 50000000    ! max nb of particles a sub-box may contain
  integer(kind=4)             :: last_part_in_subbox    
  integer(kind=4),allocatable :: particleID(:)          ! original id of the particles
  logical(kind=4),allocatable :: partInPad(:)           ! true if a particle is in pad volume

  ! sub-box boundaries in box unit (i.e. with pos from 0 to 1)
  real(kind=4)                :: xmin_subbox,xmax_subbox,ymin_subbox,ymax_subbox,zmin_subbox,zmax_subbox 
  ! sub-box boundaries in physical Mpc (_not_ comoving)
  real(kind=4)                :: xmin_subbox_pt,xmax_subbox_pt,ymin_subbox_pt,ymax_subbox_pt,zmin_subbox_pt,zmax_subbox_pt
  
  ! sub-box lengths (box units)
  real(kind=4)                :: subbox_xc,subbox_yc,subbox_zc     ! center of the sub-box
  real(kind=4)                :: subbox_xl,subbox_yl,subbox_zl     ! half side of the sub-box (with its pad)
  real(kind=4)                :: subbox_xl2,subbox_yl2,subbox_zl2  ! half side of s-b with pad minus danger zone...
  
  real(kind=4)                :: padLength,padlength_pt
  real(kind=4)                :: dangerPadLength
  logical(kind=4),allocatable :: haloHasPadPart(:)                 ! true for halos which contain a pad-particle
  integer(kind=4)             :: subboxnum                         ! number of current subbox within simulation
  integer(kind=4),parameter   :: padunit = 20
  integer(kind=4)             :: nbodies_full_sim                  ! total number of particles in the simulation

  ! re-numbering array
  integer(kind=4),allocatable :: old2new(:)
 

  ! The code assumes a total mass normalized to one. In case of a sub-box, we have only part of this mass, though
  ! we keep the same normalization -> subbox_mass_frac 
  real(kind=4)                :: subbox_mass_frac
  
contains
!***********************************************************************

  subroutine select_subbox_particles_gd(n,x,vx,id)

    implicit none
    
    integer(kind=4)             :: n
    real(kind=4)                :: x(3,n),vx(3,n)
    integer(kind=4)             :: id(n)
    integer(kind=4)             :: i,np,j,ierr
    integer(kind=4),allocatable :: list(:)
    logical(kind=4),allocatable :: padlist(:)
    real(kind=4)                :: dx,dy,dz
       
    allocate(list(n),padlist(n)) 
    list    = -1
    padlist = .false.
    np   = 0
    do i=1,n 
       dx = x(1,i) - subbox_xc
       call correct_for_periodicity_code_units(dx)
       if (abs(dx) <= subbox_xl) then 
          dy = x(2,i) - subbox_yc
          call correct_for_periodicity_code_units(dy)
          if (abs(dy) <= subbox_yl) then 
             dz = x(3,i) - subbox_zc
             call correct_for_periodicity_code_units(dz)
             if (abs(dz) <= subbox_zl) then 
                ! particle is in the padded sub-box
                np = np + 1
                list(np) = i
                ! check if it's in danger zone ... 
                if (abs(dx) >= subbox_xl2 .or. abs(dy) >= subbox_yl2 .or. abs(dz) >= subbox_zl2) then 
                   padlist(np)   = .true.
                else
                   padlist(np)   = .false.
                end if
             end if
          end if
       end if
    end do
    
    ! if there aint no particle from current input file within the sub-box, return 
    if (np == 0) then
       deallocate(list,padlist)
       return
    end if

    if (last_part_in_subbox + 1 + np <= maxnpart) then 
       do i = 1,np
          j = last_part_in_subbox + i 
          pos(j,1)      = x(1,list(i))
          pos(j,2)      = x(2,list(i))
          pos(j,3)      = x(3,list(i))
          vel(j,1)      = vx(1,list(i))
          vel(j,2)      = vx(2,list(i))
          vel(j,3)      = vx(3,list(i))
          particleID(j) = id(list(i))
          partInPad(j)  = padlist(i)
       end do
       last_part_in_subbox = j
    else
       write(errunit,*) '> ERROR: more particles in subbox than allowed -> increase MAXNPART',last_part_in_subbox+1+np
       stop
    end if

    deallocate(list,padlist)

    return

  end subroutine select_subbox_particles_gd

!***********************************************************************
  subroutine select_subbox_particles_ra(n,x,vx,id,m)
    ! leo
    ! also check if particle is a debris...

    ! exact same as above, but x and v have dim (n,3) instead  of (3,n) ... 
    implicit none
    
    integer(kind=4)             :: n
    real(kind=4)                :: x(n,3),vx(n,3)
    integer(kind=4)             :: id(n)
    real(kind=4),optional       :: m(n)
    integer(kind=4)             :: i,np,j,ierr
    integer(kind=4),allocatable :: list(:)
    logical(kind=4),allocatable :: padlist(:)
    real(kind=4)                :: dx,dy,dz
       
    allocate(list(n),padlist(n)) 
    list    = -1
    padlist = .false.
    np   = 0
    do i=1,n 
       dx = x(i,1) - subbox_xc
       call correct_for_periodicity_code_units(dx)
       if (abs(dx) <= subbox_xl) then 
          dy = x(i,2) - subbox_yc
          call correct_for_periodicity_code_units(dy)
          if (abs(dy) <= subbox_yl) then 
             dz = x(i,3) - subbox_zc
             call correct_for_periodicity_code_units(dz)
             if (abs(dz) <= subbox_zl) then 
                ! particle is in the padded sub-box
                ! now check if debris...
                if (id(i).gt.0) then
                   np = np + 1
                   list(np) = i
                   ! check if it's in danger zone ... 
                   if (abs(dx) >= subbox_xl2 .or. abs(dy) >= subbox_yl2 .or. abs(dz) >= subbox_zl2) then 
                      padlist(np)   = .true.
                   else
                      padlist(np)   = .false.
                   end if
                endif
             end if
          end if
       end if
    end do
    
    ! if there aint no particle from current input file within the sub-box, return 
    if (np == 0) then
       deallocate(list,padlist)
       return
    end if

    if (last_part_in_subbox + 1 + np <= maxnpart) then 
       do i = 1,np
          j = last_part_in_subbox + i 
          pos(j,1)      = x(list(i),1)
          pos(j,2)      = x(list(i),2)
          pos(j,3)      = x(list(i),3)
          vel(j,1)      = vx(list(i),1)
          vel(j,2)      = vx(list(i),2)
          vel(j,3)      = vx(list(i),3)
          particleID(j) = id(list(i))
          partInPad(j)  = padlist(i)
          if (present(m)) mass(j) = m(list(i))
       end do
       last_part_in_subbox = j
    else
       write(errunit,*) '> ERROR: more particles in subbox than allowed -> increase MAXNPART',last_part_in_subbox+1+np
       stop
    end if

    deallocate(list,padlist)

    return

  end subroutine select_subbox_particles_ra

!***********************************************************************
  subroutine select_subbox_particles_ra3(n,x,vx,id,t,m)

    ! exact same as above, but x and v have dim (n,3) instead  of (3,n) ... 
    implicit none
    
    integer(kind=4)             :: n
    real(kind=8)                :: x(n,3),vx(n,3),t(n)
    integer(kind=4)             :: id(n)
    real(kind=8),optional       :: m(n)
    integer(kind=4)             :: i,np,j,ierr
    integer(kind=4),allocatable :: list(:)
    logical(kind=4),allocatable :: padlist(:)
    real(kind=4)                :: dx,dy,dz
       
    allocate(list(n),padlist(n)) 
    list    = -1
    padlist = .false.
    np   = 0
    do i=1,n 
       dx = x(i,1) - subbox_xc
       call correct_for_periodicity_code_units(dx)
       if (abs(dx) <= subbox_xl) then 
          dy = x(i,2) - subbox_yc
          call correct_for_periodicity_code_units(dy)
          if (abs(dy) <= subbox_yl) then 
             dz = x(i,3) - subbox_zc
             call correct_for_periodicity_code_units(dz)
             if (abs(dz) <= subbox_zl) then 
                ! particle is in the padded sub-box
#ifndef STARS                
                if (id(i) > 0 .and. t(i) == 0.0) then ! select DM particles.
#else 
                if (id(i) > 0 .and. t(i) .ne. 0.0) then ! select star particles
#endif
                   np = np + 1
                   list(np) = i
                   ! check if it's in danger zone ... 
                   if (abs(dx) >= subbox_xl2 .or. abs(dy) >= subbox_yl2 .or. abs(dz) >= subbox_zl2) then 
                      padlist(np)   = .true.
                   else
                      padlist(np)   = .false.
                   end if
                endif
             end if
          end if
       end if
    end do
    
    ! if there aint no particle from current input file within the sub-box, return 
    if (np == 0) then
       deallocate(list,padlist)
       return
    end if

    if (last_part_in_subbox + 1 + np <= maxnpart) then 
       do i = 1,np
          j = last_part_in_subbox + i 
          pos(j,1)      = x(list(i),1)
          pos(j,2)      = x(list(i),2)
          pos(j,3)      = x(list(i),3)
          vel(j,1)      = vx(list(i),1)
          vel(j,2)      = vx(list(i),2)
          vel(j,3)      = vx(list(i),3)
          particleID(j) = id(list(i))
          partInPad(j)  = padlist(i)
          if (present(m)) mass(j) = m(list(i))
       end do
       last_part_in_subbox = j
    else
       write(errunit,*) '> ERROR: more particles in subbox than allowed -> increase MAXNPART',last_part_in_subbox+1+np
       stop
    end if

    deallocate(list,padlist)

    return

  end subroutine select_subbox_particles_ra3

!***********************************************************************

  subroutine define_padded_subbox
    
    ! define the quantities used to check whether a halo or particle is in the box. 
    implicit none 
    
    real(kind=8) :: mass_eq, xmin,xmax,ymin,ymax,zmin,zmax
    
    write(errunit,*) '> Using PadLength       = ',padLength
    write(errunit,*) '>   and dangerPadLength = ',dangerPadLength

    ! define subbox in terms of center (xc,...) and half-size (xl,...)
    xmin = xmin_subbox - padLength
    xmax = xmax_subbox + padLength
    subbox_xc   = 0.5 * (xmax + xmin)
    subbox_xl   = 0.5 * (xmax - xmin)
    ymin = ymin_subbox - padLength
    ymax = ymax_subbox + padLength
    subbox_yc   = 0.5 * (ymax + ymin) 
    subbox_yl   = 0.5 * (ymax - ymin)
    zmin = zmin_subbox - padLength
    zmax = zmax_subbox + padLength
    subbox_zc   = 0.5 * (zmax + zmin)
    subbox_zl   = 0.5 * (zmax - zmin)
    
    ! define danger zone
    subbox_xl2  = subbox_xl - dangerPadLength
    subbox_yl2  = subbox_yl - dangerPadLength
    subbox_zl2  = subbox_zl - dangerPadLength

    last_part_in_subbox = 0
    
    write(errunit,*) '------------ SUB-BOX -------------'
    write(errunit,*) '> Padded subbox is centered on :'
    write(errunit,*) '(xc,yc,zc) = ',subbox_xc,subbox_yc,subbox_zc
    write(errunit,*) '> Padded subbox has extension :'
    write(errunit,*) '(dx,dy,dz) = ',subbox_xl,subbox_yl,subbox_zl
    write(errunit,*) '----------------------------------'

    return
    
  end subroutine define_padded_subbox

!***********************************************************************

  subroutine convert_subbox_units

    ! define boundaries in Mpc too 
    implicit none
    
    ! convert to physical Mpc for later use ... 
    xmin_subbox_pt = xmin_subbox * Lboxp*(aexp/af)
    xmax_subbox_pt = xmax_subbox * Lboxp*(aexp/af)
    ymin_subbox_pt = ymin_subbox * Lboxp*(aexp/af)
    ymax_subbox_pt = ymax_subbox * Lboxp*(aexp/af)
    zmin_subbox_pt = zmin_subbox * Lboxp*(aexp/af)
    zmax_subbox_pt = zmax_subbox * Lboxp*(aexp/af)
    padlength_pt   = padlength   * lboxp*(aexp/af)
    
    return

  end subroutine convert_subbox_units

!***********************************************************************
  
  subroutine check_subbox_params

    ! check that user defined proper boundaries for sub-box and set subboxFoF to 
    ! true or false .

    implicit none
    
    if (xmin_subbox < -0.5) then 
       write(errunit,'(a)') 'xmin_subbox should be in code units (-0.5 < x < 0.5)',xmin_subbox
       stop
    end if
    if (ymin_subbox < -0.5) then 
       write(errunit,'(a)') 'ymin_subbox should be in code units (-0.5 < y < 0.5)',ymin_subbox
       stop
    end if
    if (zmin_subbox < -0.5) then 
       write(errunit,'(a)') 'zmin_subbox should be in code units (-0.5 < z < 0.5)',zmin_subbox
       stop
    end if
    if (xmax_subbox > 0.5) then 
       write(errunit,'(a)') 'xmax_subbox should be in code units (-0.5 < x < 0.5)',xmax_subbox
       stop
    end if
    if (ymax_subbox > 0.5) then 
       write(errunit,'(a)') 'ymax_subbox should be in code units (-0.5 < y < 0.5)',ymax_subbox
       stop
    end if
    if (zmax_subbox > 0.5) then 
       write(errunit,'(a)') 'zmax_subbox should be in code units (-0.5 < z < 0.5)',zmax_subbox
       stop
    end if

    ! if the whole box is asked for (either min/max not provided or explicitely set to lbox)
    ! set subboxFoF flag to false. Else set it to true to use sub-boxes.
    if (abs(xmin_subbox + 0.5) < rel_prec .and. abs(xmax_subbox - 0.5) < rel_prec .and. & 
         & abs(ymin_subbox + 0.5) < rel_prec .and. abs(ymax_subbox - 0.5) < rel_prec .and. & 
         & abs(zmin_subbox + 0.5) < rel_prec .and. abs(zmax_subbox - 0.5) < rel_prec) then 
       write(errunit,*) '> Normal run -> whole box done in one go'
       call subbox_defaults
    else
       write(errunit,*) '> Split run -> only the following sub-box will be done:'
       write(errunit,*) 'xmin,xmax = ',xmin_subbox,xmax_subbox
       write(errunit,*) 'ymin,ymax = ',ymin_subbox,ymax_subbox
       write(errunit,*) 'zmin,zmax = ',zmin_subbox,zmax_subbox
       subboxFoF = .true.
    end if
    
    return
    
  end subroutine check_subbox_params

!***********************************************************************

  subroutine subbox_defaults
  
    ! set sub-box parameters to default values before read from file
    ! Default is use the whole box.
    implicit none 
    
    xmin_subbox = -0.5
    xmax_subbox = 0.5
    ymin_subbox = -0.5
    ymax_subbox = 0.5
    zmin_subbox = -0.5
    zmax_subbox = 0.5
    padlength   = 0.0
    subboxFoF   = .false.

    return
    
  end subroutine subbox_defaults

!***********************************************************************
 
  function halo_is_in_the_box(h) 
    
    ! return true if a halo's center is in the desired sub-box, false otherwise. 
    ! Check is done in non-comoving Mpc coords (after call to change_units)
    ! NB: we have to take into account periodic boundary conditions

    implicit none
    
    type(halo)      :: h
    type(vector)    :: p
    logical(kind=4) :: halo_is_in_the_box
    
    p%x = h%p%x
    p%y = h%p%y
    p%z = h%p%z
    call correct_for_periodicity(p)
    if (p%x <= xmax_subbox_pt .and. p%x > xmin_subbox_pt .and. &
         & p%y <= ymax_subbox_pt .and. p%y > ymin_subbox_pt .and. &
         & p%z <= zmax_subbox_pt .and. p%z > zmin_subbox_pt) then 
       halo_is_in_the_box = .true.
    else
       halo_is_in_the_box = .false.
    end if
    
    return
    
  end function halo_is_in_the_box

!***********************************************************************

  subroutine read_gadget_subbox(name_file)

    ! simplified version -> DOES NOT WORK FOR RESIMULATIONS 
    ! (i.e. assumes BIG_RUN option is on)

    implicit none

    integer(kind=4)                      :: nbig,stop_read,nf,istat
    character(len=*)                     :: name_file
    character(len=3)                     :: fnumber
    character(len=len_trim(name_file)+4) :: Vfilename
    type header
       integer(kind=4) :: Vnpart(6),Vnpart_tot(6)
       real(kind=8)    :: massarr(6),npartTotalHighWord(6)
       real(kind=8)    :: time,redshift,Boxsize,Omega0,OmegaLambda,HubbleParam
       integer(kind=4) :: flag_sfr,flag_cooling,flag_feedback,flag_stellarage,flag_metals,flag_entropy_instead_u
       integer(kind=4) :: num_files
    end type header
    type(header)       :: io_header
    integer(kind=4)                      :: VN,nbaryon,id_part,idp,itype
    integer(kind=4),allocatable          :: Vid(:)
    real(kind=4),allocatable             :: Gd_pos(:,:),Gd_vel(:,:)    
    real(kind=4)                         :: Lbox_pt_local
    ! softening lengths (comoving, and maximum physical in kpc/h) used in the simulation:
    real(kind=4),parameter               :: HRepscom = 10.0, HRepsmaxph = 5.0
    real(kind=4),parameter               :: LRepscom = 10.0, LRepsmaxph = 5.0
    integer(kind=4)                      :: i,ierr
    real(kind=4)                         :: Omega0,OmegaBaryon,mass_fac,mass_box 
    character*200                        :: paramname,value,line,fileparam
    logical(kind=4),allocatable          :: padp(:)
    integer(kind=4)                      :: i1,i2,n12
    
#ifndef BIG_RUN
    write(errunit,*) '> STOP : subbox option only works with BIG_RUN'
    stop
#endif

    ! first build file name as path/snapshot_xxx.n, where:
    ! "path/snapshot_xxx" are in "name_file" variable, and 
    ! "n" is built by: //fnumber(verify(fnumber,' '):3)
    stop_read = 0
    nf        = -1
    nbig      = 1000
    nbodies   = 0
    nbaryon   = 0

    read_files: do while (nf <= nbig .and. stop_read == 0)
       istat = 1
       if(nf.lt.0) then
          Vfilename = trim(name_file)
          open(1,file=Vfilename,status='old',form='unformatted',iostat=istat)
          if(istat == 0) then
             stop_read = 1
             write(errunit,*) '> read_gadget... only one file to read'
          else
             write(errunit,*) '> read_gadget... several files to read'
          end if
       end if
       if(istat /= 0) then
          nf        = nf+1
          write(fnumber,'(i3)') nf  !--> transform the integer nf in a character
          Vfilename = trim(name_file)// '.'//fnumber(verify(fnumber,' '):3)
          print*,trim(Vfilename)
          open(1,file=Vfilename,status='old',form='unformatted',iostat=istat)
          if (istat /= 0.and.nf.le.1) stop ' > ERROR: read_gadget: input file not there'
          if(istat /=0) stop_read = 1
       endif
       if (istat /= 0) exit read_files
       write(errunit,*) '> Read file: ',trim(Vfilename)
       read (1,iostat=istat)                  &
            io_header%Vnpart,                 &  
            io_header%massarr,                &  
            io_header%time,                   &
            io_header%redshift,               &
            io_header%flag_sfr,               &
            io_header%flag_feedback,          &
            io_header%Vnpart_tot,             &
            io_header%flag_cooling,           &
            io_header%num_files,              &
            io_header%Boxsize,                &
            io_header%Omega0,                 &
            io_header%OmegaLambda,            &
            io_header%HubbleParam,            &
            io_header%flag_stellarage,        &
            io_header%flag_metals,            &
            io_header%npartTotalhighWord,     &
            io_header%flag_entropy_instead_u
       if (istat /= 0) stop ' > ERROR: read_gadget... could not read header.'
       if(nbodies.le.0) then ! if nbodies is not allocated we're reading the first file
          aexp  = io_header%time*af  ! because SN format assumes a(t_in)=1
          nbodies = io_header%Vnpart_tot(2)
          nbodies_full_sim = nbodies
          nbaryon = io_header%Vnpart_tot(1) + io_header%Vnpart_tot(5)
          if(io_header%Vnpart_tot(3).gt.0.or.io_header%Vnpart_tot(4).gt.0.or.io_header%Vnpart_tot(6).gt.0) then
             if(io_header%Vnpart(2).gt.0.and.io_header%massarr(2).eq.0.) then
                write(errunit,*)
                write(errunit,*) '> ERROR: read_gadget... mass is 0. for DM particles (type2)'
                write(errunit,*) '> Code was compiled with -DBIG_RUN option'
                write(errunit,*) '> This shouldn''t be the case'
                stop
             end if
          endif
          write(errunit,*) '> Number of DM particles            :',nbodies
          write(errunit,*) '> Number of baryons                 :',nbaryon
          
          nbodies = maxnpart
          allocate(pos(nbodies,3))
          allocate(vel(nbodies,3))
          allocate(ParticleID(nbodies),PartInPad(nbodies))
          pos        = 0.0
          vel        = 0.0
          ParticleID = -1
          PartInPad  = .false.
       end if
       
       VN    = sum(io_header%Vnpart)       
       write(errunit,*) '> Total nb of particles in file     :',sum(io_header%Vnpart)
       allocate(Gd_pos(3,VN),Gd_vel(3,VN),Vid(VN))
       read(1) Gd_pos
       read(1) Gd_vel
       read(1) Vid
       close(1)
       ! select particles (among DM particles)
       i1  = io_header%Vnpart(1)+1
       i2  = io_header%Vnpart(1)+io_header%Vnpart(2) 
       n12 = i2 - i1 + 1
       ! convert positions to box units
       gd_pos = gd_pos * io_header%time/(io_header%Boxsize*aexp) - 0.5

       ! jeje 
       print*,'VN,n12',vn,n12
       ! end jeje

       call select_subbox_particles_gd(n12,Gd_pos,gd_vel,vid)
       deallocate(Gd_pos, Gd_vel, Vid)
       
    enddo read_files
!********************************************************
    print*,'after selection:',minval(pos),maxval(pos),last_part_in_subbox
    ! re-size arrays depending on how many particles were selected
    nbodies = last_part_in_subbox

    allocate(gd_pos(nbodies,3))
    gd_pos(:,:) = pos(1:nbodies,:)
    deallocate(pos)
    allocate(pos(nbodies,3))
    pos = gd_pos

    gd_pos(:,:) = vel(1:nbodies,:)
    deallocate(vel)
    allocate(vel(nbodies,3))
    vel = gd_pos
    deallocate(gd_pos)

    allocate(vid(nbodies))
    vid = ParticleID(1:nbodies)
    deallocate(ParticleID)
    allocate(ParticleID(nbodies))
    ParticleID = vid
    deallocate(vid)

    allocate(padp(nbodies))
    padp = PartInPad(1:nbodies)
    deallocate(PartInPad)
    allocate(PartInPad(nbodies))
    PartInPad = padp
    deallocate(padp)

    print*,'after selection:',minval(pos),maxval(pos)
    ! jeje : test ... 
!!$    ierr = 0
!!$    do i=1,nbodies
!!$       if (pos(i,1) == 0.0 .and. pos(i,2) == 0 .and. pos(i,3) == 0.0) then 
!!$          !print*,i,pos(i,:)
!!$          ierr = ierr + 1
!!$       end if
!!$    end do
!!$    print*,'in pos:',ierr
    ! end jeje


    ! have to recalculate mboxp because Gadget particles are in physical units 
    mass_box = io_header%massarr(2)*10./H_f*real(nbodies_full_sim,4) ! in 10^11 M_sun
    massp    = 1./real(nbodies_full_sim,4) !in units of mboxp !io_header%massarr(2)*10./H_f/mass_box  
    
    ! define mass of subbox: 
    subbox_mass_frac = real(nbodies)/real(nbodies_full_sim)

    ! Change units of pos, vel, masses: 
    ! .. Transform comoving positions in h^-1 Kpc to physical positions in Mpc, and 
    !    normalise to physical box length in Mpc to have them in [-0.5,0.5]:
    !    Gd_pos * h^-1 * 10^-3 *io_header%time/Lbox(t), where Lbox(t)=Lboxp*aexp
    !    pos  = pos*io_header%time/(10.*H_f*Lboxp*aexp) - 0.5
    !    We can as well normalize with io_header%Boxsize = 10.*H_f*Lboxp
    ! .. Transform Gadget velocities, in SIMPLE_peculiar_velocies in units
    !    of Hubble parameter accross the box  
    !  where:
    !     Gd_vel= sqrt(io_header%time)*x_dot in km/s [io_header%time == Gadget expansion param]
    !  thus: 
    !     pec_vel=[sqrt(io_header%time)*Gd_vel/H(t)/L(t)]
    ! NB: pos is already in correct units here ... 

    Lbox_pt_local = Lboxp*(aexp/ai)
    vel           = vel *sqrt(io_header%time)

    if (nbaryon > 0) then 
       write(errunit,*)
       write(errunit,*) '> Checking the "parameters-usedvalues" file: '
       i         = index(name_file,'/snapshot_')
       fileparam = trim(name_file(1:i))//'parameters-usedvalues'
       open(unit=5,status='old',form='formatted',file=fileparam,iostat=ierr)
       if(ierr.ne.0) then
          write(errunit,*) '> Please copy or link the "parameters-usedvalues" file in directory:'
          write(errunit,*) fileparam(1:i)
          stop
       else
          Omega0       = -1.
          OmegaBaryon  = -1.0
          do
             read(5,'(a)',end=51) line
             i = scan(line,'     ')
             if(i.eq.0.or.line(1:1) .eq. '#') cycle
             paramname = trim(adjustl(line(:i-1)))
             value     = trim(adjustl(line(i+1:)))
             select case(trim(paramname))
             case('Omega0')
                read(value,*) Omega0
             case('OmegaBaryon')
                read(value,*) OmegaBaryon
             end select
          end do
51        close(5)
          write(errunit,*) '> paremeters I read form the "parameter-usedvalues" file'
          write(errunit,*) '> Omega0       : ',Omega0
          write(errunit,*) '> OmegaBaryon  : ',OmegaBaryon
          if (Omega0 .le. 0.) stop ' > read_gadget... Coudn''t read Omega0 value'
          if (OmegaBaryon .le. 0.) stop ' > read_gadget... Coudn''t read OmegaBaryon value'
          write(errunit,*)
       end if
       mass_fac = 1./(1.0-OmegaBaryon/Omega0)
       mass_box = mass_box * mass_fac
#ifndef RENORM
       massp = massp/mass_fac
#endif
    endif

    write(errunit,*) '> In Gadget'
    write(errunit,*) '> box mass (M_sun)            :',mass_box*1.e11
    write(errunit,*) '> particle mass (in M_sun)    :',massp*mass_box*1.e11         
#ifdef RENORM
    if(nbaryon > 0) then
       write(errunit,*) '> after renormalisation'
       write(errunit,*) '> box mass (M_sun)         : ',mass_box*1e11
       write(errunit,*) '> particle mass (M_sun)    :', massp*mass_box*1e11
    end if
#endif
    write(errunit,*)
    write(errunit,*) '> In HaloMaker, after rescaling'
    write(errunit,*) '> box mass (M_sun)            :',mboxp*1.e11
    write(errunit,*) '> particle mass (in M_sun)    :',massp*mboxp*1.e11         
#ifdef RENORM
    if(nbaryon > 0) then
       write(errunit,*) '> after renormalisation'
       write(errunit,*) '> box mass (M_sun)         : ',mboxp*1e11
       write(errunit,*) '> particle mass (M_sun)    :', massp*mboxp*1e11
    end if
#endif

    return

  end subroutine read_gadget_subbox

  !***********************************************************************                
  subroutine read_ramses_subbox(repository)
    ! -> same as read_ramses, but selects particles in a sub-box... 
    
    implicit none

    character(len=*)            :: repository
    integer(kind=4)             :: ndim,npart,idim,icpu,ipos,ncpu,i,ipar
    integer(kind=4)             :: ncpu2,npart2,ndim2
    integer(kind=4)             :: nx,ny,nz,nlevelmax,ngridmax,nstep_coarse
    integer(kind=4),allocatable :: idp(:)
    real(kind=4)                :: boxlen,tco,aexp_ram,hexp
    real(kind=4)                :: omega_m,omega_l,omega_k,omega_b
    real(kind=4)                :: scale_l,scale_d,scale_t
    real(kind=8)                :: mtot,massres
    real(kind=4),allocatable    :: tmpp(:,:),tmpv(:,:),tmpm(:)
    character*200               :: nomfich
    character*5                 :: nchar,ncharcpu
    logical                     :: ok
    ! jeje 
    logical(kind=4),allocatable :: padp(:)
    real(kind=4)    :: simminmass
    ! end jeje


#ifndef BIG_RUN
    write(errunit,*) '> STOP : subbox option only works with BIG_RUN'
    stop
#endif

    ! NB: repository is directory containing output files
    ! e.g. /horizon1/teyssier/ramses_simu/boxlen100_n256/output_00001/

    ! read cosmological params in header of amr file
    ipos    = index(repository,'output_')
    nchar   = repository(ipos+7:ipos+13)
    nomfich = trim(repository)//'amr_'//trim(nchar)//'.out00001'
    inquire(file=nomfich,exist=ok)
    if (.not. ok) then
       write(errunit,*)'File '//trim(nomfich)//' not found'
       stop
    else
       open(unit=10,file=nomfich,status='old',form='unformatted')
       read(10) ncpu
       read(10) ndim
       read(10) nx,ny,nz
       read(10) nlevelmax
       read(10) ngridmax
       read(10) nstep_coarse    
       read(10) boxlen
       ! temps conforme tau, expansion factor, da/dtau
       read(10) tco,aexp_ram,hexp
       read(10) omega_m,omega_l,omega_k,omega_b
       ! to get units cgs multiply by these scale factors
       read(10) scale_l,scale_d,scale_t
       close(10)
       !write(errunit,991) ncpu,ndim,nstep_coarse
       !write(errunit,993) nlevelmax,ngridmax
       !write(errunit,997) boxlen*scale_l
       !write(errunit,994) t,aexp,hexp
       !write(errunit,995) omega_m,omega_l,omega_k,omega_b
    end if
    ! use approximate comv from cm to Mpc to match Romain's conversion... 
    Lboxp          = boxlen*scale_l/3.08e24/aexp_ram ! converts cgs to Mpc comoving
    !write(errunit,*) 'af,hf,lboxp,ai,aexp',af,h_f,lboxp,ai,aexp_ram
    aexp           = aexp_ram*af  
    omega_f        = omega_m+omega_b
    omega_lambda_f = omega_l
    omega_c_f      = omega_k
990 format(' Enter output number:',i6)
991 format(' ncpu=',i6,' ndim=',i1,' nstep=',i6)
992 format(' nx=',i3,' ny=',i3,' nz=',i3)
993 format(' nlevelmax=',i3,' ngridmax=',i8)
994 format(' t=',1pe10.3,' aexp=',1pe10.3,' hexp=',1pe10.3)
995 format(' omega_m=',F6.3,' omega_l=',F6.3,' omega_k=',F6.3,' omega_b=',F6.3)
997 format(' boxlen=',1pe10.3,' h-1 Mpc')

    ! now read the particle data files
    nomfich = trim(repository)//'/part_'//trim(nchar)//'.out00001'
    inquire(file=nomfich,exist=ok) ! verify input file
    if ( .not. ok ) then
       write(errunit,*) trim(nomfich)//' not found.'
       stop
    endif
    
    open(unit=1,file=nomfich,status='old',form='unformatted')
    read(1) ncpu
    read(1) ndim
    close(1)
    
    ! get total number of particles in the simulation
    npart = 0
    do icpu = 1,ncpu
       call title(icpu,ncharcpu)
       nomfich = trim(repository)//'/part_'//trim(nchar)//'.out'//trim(ncharcpu)
       !write(errunit,*)'> Reading file '//trim(nomfich)
       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1) ncpu2
       read(1) ndim2
       read(1) npart2
       close(1)
       npart = npart+npart2
    end do
    
    ! jeje 
    !!$ nbodies = npart
    nbodies_full_sim = npart
    ! end jeje
    write(errunit,*)'> Found ',npart,' particles'
    write(errunit,*)'> Reading positions and masses...'

    ! jeje 
!!$    allocate(pos(1:npart,1:ndim))
!!$    allocate(vel(1:npart,1:ndim))
!!$    allocate(mass(1:npart))
    nbodies = maxnpart
    allocate(pos(1:nbodies,1:ndim))
    pos = 0.0
    allocate(vel(1:nbodies,1:ndim))
    pos = 0.0
    ! we don't need mass 'cause we aint dealing with re-simulations ... 
!!$    allocate(mass(1:nbodies))
!!$    mass = 0.0
    allocate(ParticleID(nbodies),PartInPad(nbodies))
    ParticleID = -1
    PartInPad  = .false.
    !end jeje 

    do icpu = 1,ncpu
       call title(icpu,ncharcpu)
       nomfich = trim(repository)//'/part_'//trim(nchar)//'.out'//trim(ncharcpu)
       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1) ncpu2
       read(1) ndim2
       read(1) npart2
       allocate(tmpp(1:npart2,1:ndim),tmpv(1:npart2,1:ndim),tmpm(1:npart2),idp(1:npart2))
       ! read all particle positions
       do idim = 1,ndim
          read(1) tmpp(1:npart2,idim)
       end do
       ! read all particle velocities
       do idim = 1,ndim
          read(1) tmpv(1:npart2,idim)
       end do
       ! read all particle masses
       read(1) tmpm(1:npart2)
       ! jeje
       if (icpu == 1) then 
          simminmass = minval(tmpm)
       end if
       ! end jeje 
       ! read all particle ids
       read(1) idp(1:npart2)
       close(1)
       
       ! jeje 
       ! 1/ convert positions before selection -> box units ranging from -0.5 to 0.5 
       tmpp = tmpp - 0.5
       ! 2/ select particles into pos,vel, and ParticleID arrays
       call select_subbox_particles_ra(npart2,tmpp,tmpv,idp)
       ! 3/ comment stuff below...
!!$       ! now sort DM particles in ascending id order
!!$       do ipar=1,npart2 
!!$          ! put all positions between -0.5 and 0.5
!!$          pos(idp(ipar),1:ndim) = tmpp(ipar,1:ndim) - 0.5
!!$          ! convert code units to km/s 
!!$          vel(idp(ipar),1:ndim) = tmpv(ipar,1:ndim)*scale_l/scale_t*1e-5
!!$          mass(idp(ipar))       = tmpm(ipar)
!!$       end do
       ! end jeje 
       deallocate(tmpp,tmpv,tmpm,idp)
    end do

    ! jeje 
    ! after selection, resize arrays, declare nbodies, and change velocities (as was done on the read before)
    nbodies = last_part_in_subbox
    write(errunit,*) '> SUBBOX: nb of selected particles = ',nbodies
    allocate(tmpp(nbodies,3))  ! screw ndim /= 3 !! 
    tmpp = pos(1:nbodies,:)
    deallocate(pos)
    allocate(pos(nbodies,3))
    pos  = tmpp
    tmpp = vel(1:nbodies,:)
    deallocate(vel)
    allocate(vel(nbodies,3))
    vel = tmpp * scale_l/scale_t*1e-5  ! unit conversion : -> km/s (i guess)
    deallocate(tmpp)
    allocate(idp(nbodies))
    idp = ParticleID(1:nbodies)
    deallocate(ParticleID)
    allocate(ParticleID(nbodies))
    ParticleID = idp
    deallocate(idp)
    allocate(padp(nbodies))
    padp = PartInPad(1:nbodies)
    deallocate(PartInPad)
    allocate(PartInPad(nbodies))
    PartInPad = padp
    deallocate(padp)
    ! end jeje

    ! jeje : all that mass business ... I have to define massp without reading the whole snapshot ... 
    ! ... and with a somewhat confused understanding of the RENORM option ... 
#ifdef RENORM 
    massp = 1.0 / nbodies_full_sim ! ta daaah ! 
    print*,'+++++++++++++++++ massp from jeje +++++++++++++++++++> ',massp,real(1.0d0/dble(nbodies_full_sim),4)
    print*,'+++++++++++++++++ minval(mass)    +++++++++++++++++++> ',simminmass
#else 
    massp = simminmass 
#endif


!!$    mtot = 0.0d0
!!$    do i = 1,npart
!!$       mtot = mtot+real(mass(i),8)
!!$    enddo
!!$    ! that is for the dark matter so let's add baryons now if there are any 
!!$    ! and renormalization flag is on !!
!!$    massres = minval(mass)*mboxp*1d11
!!$    massp   = minval(mass)
!!$    write(errunit,*) '> particle mass (in M_sun)               = ',massres
!!$#ifdef RENORM
!!$    massres = minval(mass)*mboxp*1d11/mtot
!!$    massp   = minval(mass)/real(mtot,4)
!!$    write(errunit,*) '> particle mass (in M_sun) after renorm  = ',massres
!!$#endif
!!$#ifdef BIG_RUN
!!$    deallocate(mass)
!!$#endif
    !end jeje ... 
    
  end subroutine read_ramses_subbox

  !***********************************************************************                
  subroutine read_ramses_star_subbox(repository)
    ! -> same as read_ramses, but selects particles in a sub-box... 
    ! Leo -> for star particles...

    implicit none

    character(len=*)            :: repository
    integer(kind=4)             :: ndim,npart,idim,icpu,ipos,ncpu,i,ipar
    integer(kind=4)             :: ncpu2,npart2,ndim2
    integer(kind=4)             :: nx,ny,nz,nlevelmax,ngridmax,nstep_coarse
    integer(kind=4),allocatable :: idp(:)
    real(kind=4)                :: boxlen,tco,aexp_ram,hexp
    real(kind=4)                :: omega_m,omega_l,omega_k,omega_b
    real(kind=4)                :: scale_l,scale_d,scale_t
    real(kind=8)                :: mtot,massres
    real(kind=4),allocatable    :: tmpp(:,:),tmpv(:,:),tmpm(:)
    character*200               :: nomfich
    character*5                 :: nchar,ncharcpu
    logical                     :: ok
    ! jeje 
    logical(kind=4),allocatable :: padp(:)
    !real(kind=4)    :: simminmass
    ! end jeje
    real(kind=8)                :: aaa,bbb


#ifndef BIG_RUN
    write(errunit,*) '> STOP : subbox option only works with BIG_RUN'
    stop
#endif

    ! leo debug
    write(errunit,*) '> you are in read_ramses_star_subbox'

    ! NB: repository is directory containing output files
    ! e.g. /horizon1/teyssier/ramses_simu/boxlen100_n256/output_00001/

    ! read cosmological params in header of amr file
    ipos    = index(repository,'output_')
    nchar   = repository(ipos+7:ipos+13)
    nomfich = trim(repository)//'amr_'//trim(nchar)//'.out00001'
    inquire(file=nomfich,exist=ok)
    if (.not. ok) then
       write(errunit,*)'File '//trim(nomfich)//' not found'
       stop
    else
       open(unit=10,file=nomfich,status='old',form='unformatted')
       read(10) ncpu
       read(10) ndim
       read(10) nx,ny,nz
       read(10) nlevelmax
       read(10) ngridmax
       read(10) nstep_coarse    
       read(10) boxlen
       ! temps conforme tau, expansion factor, da/dtau
       read(10) tco,aexp_ram,hexp
       read(10) omega_m,omega_l,omega_k,omega_b
       ! to get units cgs multiply by these scale factors
       read(10) scale_l,scale_d,scale_t
       close(10)
       !write(errunit,991) ncpu,ndim,nstep_coarse
       !write(errunit,993) nlevelmax,ngridmax
       !write(errunit,997) boxlen*scale_l
       !write(errunit,994) t,aexp,hexp
       !write(errunit,995) omega_m,omega_l,omega_k,omega_b
    end if
    ! use approximate comv from cm to Mpc to match Romain's conversion... 
    Lboxp          = boxlen*scale_l/3.08e24/aexp_ram ! converts cgs to Mpc comoving
    !write(errunit,*) 'af,hf,lboxp,ai,aexp',af,h_f,lboxp,ai,aexp_ram
    aexp           = aexp_ram*af  
    omega_f        = omega_m+omega_b
    omega_lambda_f = omega_l
    omega_c_f      = omega_k
990 format(' Enter output number:',i6)
991 format(' ncpu=',i6,' ndim=',i1,' nstep=',i6)
992 format(' nx=',i3,' ny=',i3,' nz=',i3)
993 format(' nlevelmax=',i3,' ngridmax=',i8)
994 format(' t=',1pe10.3,' aexp=',1pe10.3,' hexp=',1pe10.3)
995 format(' omega_m=',F6.3,' omega_l=',F6.3,' omega_k=',F6.3,' omega_b=',F6.3)
997 format(' boxlen=',1pe10.3,' h-1 Mpc')


!     print*,hexp
!     print*,omega_m
!     print*,omega_b
!     print*,omega_m/omega_b
!     !!print*,omega_b/omega_m
!     print*,scale_l
!     print*,scale_d
!     aaa = (scale_d * scale_l) * (scale_l/1.99d33) * scale_l
!     bbb = Lboxp**3 * omega_b * 2.78d11 * (H_f/100.)**2
!     print*,aaa,bbb,aaa/bbb
!     stop

    ! now read the star particle data files
    nomfich = trim(repository)//'/star_'//trim(nchar)//'.out00001'
    inquire(file=nomfich,exist=ok) ! verify input file
    if ( .not. ok ) then
       write(errunit,*) trim(nomfich)//' not found.'
       stop
    endif
    
    open(unit=1,file=nomfich,status='old',form='unformatted')
    read(1) ncpu
    read(1) ndim
    close(1)
    
    ! get total number of star particles in the simulation
    npart = 0
    do icpu = 1,ncpu
       call title(icpu,ncharcpu)
       nomfich = trim(repository)//'/star_'//trim(nchar)//'.out'//trim(ncharcpu)
       !write(errunit,*)'> Reading file '//trim(nomfich)
       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1) ncpu2
       read(1) ndim2
       read(1) npart2
       close(1)
       npart = npart+npart2
    end do
    
    ! jeje 
    !!$ nbodies = npart
    nbodies_full_sim = npart
    ! end jeje
    write(errunit,*)'> Found ',npart,' particles'
    write(errunit,*)'> Reading positions and masses...'

    ! jeje 
!!$    allocate(pos(1:npart,1:ndim))
!!$    allocate(vel(1:npart,1:ndim))
!!$    allocate(mass(1:npart))
    nbodies = maxnpart
    allocate(pos(1:nbodies,1:ndim))
    pos = 0.0
    allocate(vel(1:nbodies,1:ndim))
    pos = 0.0
    ! we don't need mass 'cause we aint dealing with re-simulations ... 
    ! Leo: yes we need mass particles since star particles have different mass...
    allocate(mass(1:nbodies))
    mass = 0.0
    allocate(ParticleID(nbodies),PartInPad(nbodies))
    ParticleID = -1
    PartInPad  = .false.
    !end jeje 

    do icpu = 1,ncpu
       call title(icpu,ncharcpu)
       nomfich = trim(repository)//'/star_'//trim(nchar)//'.out'//trim(ncharcpu)
       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1) ncpu2
       read(1) ndim2
       read(1) npart2
       allocate(tmpp(1:npart2,1:ndim),tmpv(1:npart2,1:ndim),tmpm(1:npart2),idp(1:npart2))
       ! read all particle positions
       do idim = 1,ndim
          read(1) tmpp(1:npart2,idim)
       end do
       ! read all particle velocities
       do idim = 1,ndim
          read(1) tmpv(1:npart2,idim)
       end do
       ! read all particle masses
       read(1) tmpm(1:npart2)
       ! Leo -> we don't need simminmass
!        ! jeje
!        if (icpu == 1) then 
!           simminmass = minval(tmpm)
!        end if
!        ! end jeje 
       ! read all particle ids
       read(1) idp(1:npart2)
       close(1)
       
       ! jeje 
       ! 1/ convert positions before selection -> box units ranging from -0.5 to 0.5 
       tmpp = tmpp - 0.5
       ! 2/ select particles into pos,vel, and ParticleID arrays
       !    leo: and select true stars and exclude debris...
       call select_subbox_particles_ra(npart2,tmpp,tmpv,idp,m=tmpm)
       ! 3/ comment stuff below...
!!$       ! now sort DM particles in ascending id order
!!$       do ipar=1,npart2 
!!$          ! put all positions between -0.5 and 0.5
!!$          pos(idp(ipar),1:ndim) = tmpp(ipar,1:ndim) - 0.5
!!$          ! convert code units to km/s 
!!$          vel(idp(ipar),1:ndim) = tmpv(ipar,1:ndim)*scale_l/scale_t*1e-5
!!$          mass(idp(ipar))       = tmpm(ipar)
!!$       end do
       ! end jeje 
       deallocate(tmpp,tmpv,tmpm,idp)
    end do

    ! jeje 
    ! after selection, resize arrays, declare nbodies, and change velocities (as was done on the read before)
    nbodies = last_part_in_subbox
    write(errunit,*) '> SUBBOX: nb of selected particles = ',nbodies
    allocate(tmpp(nbodies,3))  ! screw ndim /= 3 !! 
    tmpp = pos(1:nbodies,:)
    deallocate(pos)
    allocate(pos(nbodies,3))
    pos  = tmpp
    tmpp = vel(1:nbodies,:)
    deallocate(vel)
    allocate(vel(nbodies,3))
    vel = tmpp * scale_l/scale_t*1e-5  ! unit conversion : -> km/s (i guess)
    deallocate(tmpp)
    allocate(tmpm(nbodies))
    tmpm = mass(1:nbodies)
    deallocate(mass)
    allocate(mass(nbodies))
    mass  = tmpm
    deallocate(tmpm)
    allocate(idp(nbodies))
    idp = ParticleID(1:nbodies)
    deallocate(ParticleID)
    allocate(ParticleID(nbodies))
    ParticleID = idp
    deallocate(idp)
    allocate(padp(nbodies))
    padp = PartInPad(1:nbodies)
    deallocate(PartInPad)
    allocate(PartInPad(nbodies))
    PartInPad = padp
    deallocate(padp)
    ! end jeje


    ! leo : mass stuff -> convert mass in unit of per total baryonic mass
!     msuncgs = 
!     rho_c   = 
!     scale_m     = scale_d * scale_l^3  ! convert code unit to cgs
!     scalefactor = scale_m/msuncgs      ! convert cgs to Msun
!     scalefactor = scalefactor / (Lboxp^3. * omega_b * rhoc_c)
!     mass *= scalefactor
    mass = mass * omega_m/omega_b

    massp = SUM(mass)/dble(nbodies)
    print*,'+++++++++++++++++++> ',massp,SUM(mass),nbodies
    print*,'paso dealloc '


    ! jeje : all that mass business ... I have to define massp without reading the whole snapshot ... 
    ! ... and with a somewhat confused understanding of the RENORM option ... 
! #ifdef RENORM 
!     massp = 1.0 / nbodies_full_sim ! ta daaah ! 
!     print*,'+++++++++++++++++ massp from jeje +++++++++++++++++++> ',massp,real(1.0d0/dble(nbodies_full_sim),4)
!     print*,'+++++++++++++++++ minval(mass)    +++++++++++++++++++> ',simminmass
! #else 
!     massp = simminmass 
! #endif
    ! Leo : 
#ifdef RENORM
    write(errunit,*) '> STOP : GalaxyMaker does not work with RENORM'
    stop
#endif
    ! end leo...

!!$    mtot = 0.0d0
!!$    do i = 1,npart
!!$       mtot = mtot+real(mass(i),8)
!!$    enddo
!!$    ! that is for the dark matter so let's add baryons now if there are any 
!!$    ! and renormalization flag is on !!
!!$    massres = minval(mass)*mboxp*1d11
!!$    massp   = minval(mass)
!!$    write(errunit,*) '> particle mass (in M_sun)               = ',massres
!!$#ifdef RENORM
!!$    massres = minval(mass)*mboxp*1d11/mtot
!!$    massp   = minval(mass)/real(mtot,4)
!!$    write(errunit,*) '> particle mass (in M_sun) after renorm  = ',massres
!!$#endif
!!$#ifdef BIG_RUN
!!$    deallocate(mass)
!!$#endif
    !end jeje ... 

    ! leo debug
    print*,' minval(mass) = ',minval(mass)
    print*,' maxval(mass) = ',maxval(mass)
    print*,' minval(ParticleID) = ',minval(ParticleID)
    print*,' maxval(ParticleID) = ',maxval(ParticleID)
    ! print*,' minval(x) = ',minval(pos(:,1))
!     print*,' maxval(x) = ',maxval(pos(:,1))
!     print*,' minval(y) = ',minval(pos(:,2))
!     print*,' maxval(y) = ',maxval(pos(:,2))
!     print*,' minval(z) = ',minval(pos(:,3))
!     print*,' maxval(z) = ',maxval(pos(:,3))

    ! leo : test : dump star particles...
    open(unit=233,file='catstars.dat',form='unformatted',status='unknown')
    write(233) nbodies
    write(233) pos(:,1)
    write(233) pos(:,2)
    write(233) pos(:,3)
    write(233) ParticleID
    !write(233) vel(:,1)
    !write(233) vel(:,2)
    !write(233) vel(:,3)
    close(233)

    
    deallocate(mass)


  end subroutine read_ramses_star_subbox

  !***********************************************************************                

  subroutine read_ramses_new_stars_subbox(repository)
    
    implicit none
    
    character(len=*)            :: repository

    print*,'subroutine READ_RAMSES_NEW_STARS_SUBBOX is yet to be written...'
    stop
    
    return

  end subroutine read_ramses_new_stars_subbox

  !***********************************************************************                

  subroutine read_ramses_new_subbox(repository)

    implicit none

    character(len=*)            :: repository
    integer(kind=4)             :: ndim,npart,idim,icpu,ipos,ncpu,i,ipar
    integer(kind=4)             :: ncpu2,npart2,ndim2,idum,nout,nsink,nstar
    integer(kind=4)             :: nx,ny,nz,nlevelmax,ngridmax,nstep_coarse
    integer(kind=4),allocatable :: levelp(:)
    integer(kind=4),allocatable :: idp(:)
    real(kind=8)                :: boxlen,tco,aexp_ram,hexp
    real(kind=8)                :: omega_m,omega_l,omega_k,omega_b
    real(kind=8)                :: scale_l,scale_d,scale_t,dummy
    real(kind=8)                :: mtot,massres
    real(kind=8),allocatable    :: dumout(:),tmpp(:,:),tmpv(:,:),tmpm(:),tmpt(:)
    character(len=200)          :: line,name,value,nomfich
    character(len=5)            :: nchar,ncharcpu
    logical                     :: ok
    integer(kind=4) :: cnt,cntneg
    real(kind=8),parameter :: gramm_to_1011Msun = 1.d0 / 1.9891d+33 / 1.d11
    ! jeje 
    logical(kind=4),allocatable :: padp(:)

    ! read cosmological params in header of amr file
    ipos    = index(repository,'output_')
    nchar   = repository(ipos+7:ipos+13)
    nomfich = trim(repository)//'amr_'//trim(nchar)//'.out00001'
    inquire(file=nomfich,exist=ok)
    if (.not. ok) then
       write(errunit,*)'File '//trim(nomfich)//' not found'
       stop
    else
       open(unit=10,file=nomfich,status='old',form='unformatted')
       read(10) ncpu
       read(10) ndim
       read(10) nx,ny,nz
       read(10) nlevelmax
       read(10) ngridmax
       read(10) idum
       read(10) idum
       read(10) boxlen
       read(10) nout,idum,idum
       allocate(dumout(nout))
       read(10) dumout
       read(10) dumout
       deallocate(dumout)
       read(10) tco
       allocate(dumout(nlevelmax))
       read(10) dumout
       read(10) dumout
       deallocate(dumout)
       read(10) idum,nstep_coarse    
       read(10) dummy
       read(10) omega_m,omega_l,omega_k,omega_b,dummy
       ! expansion factor, da/dtau
       read(10) aexp_ram,hexp
       close(10)
    end if

    nomfich = trim(repository)//'info_'//trim(nchar)//'.txt'
    inquire(file=nomfich,exist=ok)
    if (.not. ok) then
       write(errunit,*)'File '//trim(nomfich)//' not found'
       stop
    else
       open(unit=10,file=nomfich,status='old',form='formatted')
       do
          read (10,'(a)',end=2) line
          i = scan(line,'=')
          if (i == 0 .or. line(1:1) == '#') cycle
          name  = trim(adjustl(line(:i-1)))
          value = trim(adjustl(line(i+1:)))
          ! check for a comment at end of line !
          i     = scan(value,'!')
          if (i /= 0) value = trim(adjustl(value(:i-1)))
          select case (trim(name))
          case ('unit_l' , 'scale_l')
             read(value,*) scale_l
          case ('unit_d' , 'scale_d')
             read(value,*) scale_d
          case ('unit_t' , 'scale_t')
             read(value,*) scale_t
          end select
       end do
2      close(10)
    endif

    Lboxp          = boxlen*scale_l/3.08e24/aexp_ram ! converts cgs to Mpc comoving
    aexp           = aexp_ram*af  
    omega_f        = omega_m
    omega_lambda_f = omega_l
    omega_c_f      = omega_k

    ! now read the particle data files
    nomfich = trim(repository)//'/part_'//trim(nchar)//'.out00001'
    inquire(file=nomfich,exist=ok) ! verify input file
    if ( .not. ok ) then
       write(errunit,*) trim(nomfich)//' not found.'
       stop
    endif
    open(unit=1,file=nomfich,status='old',form='unformatted')
    read(1) ncpu
    read(1) ndim
    close(1)
    
    npart = 0
    do icpu = 1,ncpu
       call title(icpu,ncharcpu)
       nomfich = trim(repository)//'/part_'//trim(nchar)//'.out'//trim(ncharcpu)
       !write(errunit,*)'> Reading file '//trim(nomfich)
       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1) ncpu2
       read(1) ndim2
       read(1) npart2
       read(1) idum
       read(1) nstar
       close(1)
       npart = npart+npart2
    end do
    npart   = npart - nstar
    ! jeje 
    !!$ nbodies = npart
    nbodies_full_sim = npart
    ! end jeje

    write(errunit,*)'> Found ',npart,' particles'
    write(errunit,*)'> Reading positions, velocities and masses...'
    
    ! jeje 
!!$    allocate(pos(1:npart,1:ndim))
!!$    allocate(vel(1:npart,1:ndim))
!!$    allocate(mass(1:npart))
    nbodies = maxnpart
    allocate(pos(1:nbodies,1:ndim))
    pos = 0.0
    allocate(vel(1:nbodies,1:ndim))
    pos = 0.0
#ifndef BIG_RUN
    allocate(mass(1:nbodies))
    mass = 0.0
#endif
    ! jeje 
    allocate(ParticleID(npart),PartInPad(nbodies))
    particleID = -1
    PartInPad  = .false.
    ! end jeje 

    cnt = 0
    cntneg = 0
    do icpu = 1,ncpu
       call title(icpu,ncharcpu)
       nomfich = trim(repository)//'/part_'//trim(nchar)//'.out'//trim(ncharcpu)

       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1) ncpu2
       read(1) ndim2
       read(1) npart2
       read(1) idum
       read(1) nstar
       read(1) dummy   
       read(1) dummy
       read(1) nsink

       allocate(tmpp(1:npart2,1:ndim),tmpv(1:npart2,1:ndim),tmpm(1:npart2),tmpt(1:npart2),idp(1:npart2),levelp(1:npart2))
       allocate(dumout(1:npart2))
       ! read all particle positions
       do idim = 1,ndim
          read(1) tmpp(1:npart2,idim)
       end do

       ! read all particle velocities
       do idim = 1,ndim
          read(1) tmpv(1:npart2,idim)
       end do

       ! read all particle masses
       read(1) tmpm(1:npart2)

       ! read all particle ids
       read(1) idp(1:npart2)

       read(1) levelp(1:npart2)

       ! read all particle creation times if necessary
       tmpt = 0.0
       if (nstar .gt. 0) then 
          read(1) tmpt(1:npart2)
          read(1) dumout(1:npart2)
       endif
       
       close(1)

       ! jeje 
       ! 1/ convert positions before selection -> box units ranging from -0.5 to 0.5 
       tmpp = tmpp - 0.5
       ! 2/ select particles into pos,vel, and ParticleID arrays
       call select_subbox_particles_ra3(npart2,tmpp,tmpv,idp,tmpt,m=tmpm)
!!$       ! now sort DM particles in ascending id order and get rid of stars !!!
!!$       do ipar=1,npart2 
!!$          if (idp(ipar) < 0) cntneg = cntneg + 1
!!$          if (tmpt(ipar) == 0.0) then 
!!$             cnt = cnt + 1
!!$             ! put all positions between -0.5 and 0.5
!!$             pos(idp(ipar),1:ndim) = real(tmpp(ipar,1:ndim),4) - 0.5
!!$             ! convert code units to km/s 
!!$             vel(idp(ipar),1:ndim) = real(tmpv(ipar,1:ndim),4)*scale_l/scale_t*1e-5
!!$             mass(idp(ipar))       = real(tmpm(ipar),4)
!!$             particleID(idp(ipar)) = idp(ipar)
!!$          end if
!!$       end do
       deallocate(tmpp,tmpv,tmpm,tmpt,idp,levelp)
       deallocate(dumout)
    end do
   
    ! jeje 
    ! after selection, resize arrays, declare nbodies, and change velocities (as was done on the read before)
    nbodies = last_part_in_subbox
    write(errunit,*) '> SUBBOX: nb of selected particles = ',nbodies
    allocate(tmpp(nbodies,3))  ! screw ndim /= 3 !! 
    tmpp = pos(1:nbodies,:)
    deallocate(pos)
    allocate(pos(nbodies,3))
    pos  = tmpp
    tmpp = vel(1:nbodies,:)
    deallocate(vel)
    allocate(vel(nbodies,3))
    vel = tmpp * scale_l/scale_t*1e-5  ! unit conversion : -> km/s (i guess)
    deallocate(tmpp)
    allocate(idp(nbodies))
    idp = ParticleID(1:nbodies)
    deallocate(ParticleID)
    allocate(ParticleID(nbodies))
    ParticleID = idp
    deallocate(idp)
    allocate(padp(nbodies))
    padp = PartInPad(1:nbodies)
    deallocate(PartInPad)
    allocate(PartInPad(nbodies))
    PartInPad = padp
    deallocate(padp)

    !! Leo: reshape mass array (useful later for contamination check)
    allocate(tmpm(nbodies))
    tmpm = mass(1:nbodies)
    deallocate(mass)
    allocate(mass(nbodies))
    mass = tmpm
    deallocate(tmpm)

!!$ 
!!$    ! re-allocate arrays with new nbodies = cnt
!!$    ! ('cause perhaps there might be debris particles around, not counted as stars ...) 
!!$    if (nbodies /= cnt) then 
!!$       nbodies = cnt
!!$       npart   = nbodies
!!$       allocate(tmpp(1:nbodies,1:ndim),tmpv(1:nbodies,1:ndim),tmpm(1:nbodies),idp(1:nbodies))
!!$       tmpp = pos(1:cnt,:)
!!$       tmpv = vel(1:cnt,:)
!!$       tmpm = mass(1:cnt)
!!$       idp  = particleID(1:cnt)
!!$       deallocate(pos,vel,mass,particleID)
!!$       allocate(pos(1:nbodies,1:ndim),vel(1:nbodies,1:ndim),mass(1:nbodies),ParticleID(1:nbodies))
!!$       pos        = tmpp
!!$       vel        = tmpv
!!$       mass       = tmpm
!!$       particleID = idp
!!$       deallocate(tmpp,tmpv,tmpm,idp)
!!$    end if
!!$    print*,nbodies,cntneg
    ! end jeje 



    ! convert masses to box units (box mass is defined in compute_halo_props:init_cosmo
    mass  = mass * (scale_d * scale_l**3 * gramm_to_1011Msun)  ! DM particle mass in 10^11 Msun
    mass  = mass / mboxp ! where mboxp is total mass of the box in 10^11Msun
    massp = minval(mass)
    write(errunit,*) '> min particle mass (in 10d11 Msun)               = ',massp*mboxp
    write(errunit,*) '> max particle mass (in 10d11 Msun)               = ',maxval(mass) * mboxp

!!$    mtot = 0.0d0
!!$    do i = 1,npart
!!$       mtot = mtot + real(mass(i),8)
!!$    enddo
!!$    ! that is for the dark matter so let's add baryons now if there are any 
!!$    ! and renormalization flag is on !!
!!$    massres = minval(mass)*mboxp*1d11
!!$    massp   = minval(mass)
!!$    write(errunit,*) '> particle mass (in M_sun)               = ',massres
!!$#ifdef RENORM
!!$    massres = minval(mass)*mboxp*1d11/mtot
!!$    massp   = minval(mass)/real(mtot,4)
!!$    write(errunit,*) '> particle mass (in M_sun) after renorm  = ',massres
!!$#endif
#ifdef RENORM
    write(errunit,*) '> RENORM to be re-implemented later ... '
    stop
#endif

#ifdef BIG_RUN
    if (minval(mass) /= maxval(mass)) then 
       write(errunit,*) 'we have particles with different masses ... '
       write(errunit,*) minval(mass), maxval(mass)
       stop
    end if
    deallocate(mass)
#endif
        
    return

  end subroutine read_ramses_new_subbox

  !***********************************************************************                
end module subbox

