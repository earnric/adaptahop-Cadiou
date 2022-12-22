module input_output

  use halo_defs

  public

contains

!///////////////////////////////////////////////////////////////////////
!***********************************************************************
  subroutine read_data

    ! This routine read the output of N-body simulations (particles positions and speeds, 
    ! cosmological and technical parameters)

    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ! WARNING: this routine just reads the data and converts positions       !
    !          and velocities from CODE units to these units                 !
    !          -- positions are between -0.5 and 0.5                         !
    !          -- velocities are in km/s                                     !
    !          -- total box mass is 1.0                                      !
    !             for simulation with hydro (only with -DRENORM) flag        !
    !          -- initial (beg of simulation) expansion factor is ai=1.0     ! 
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    implicit none

    character(len=200) :: name_of_file

    write(errunit,*)
    write(errunit,*) '> In read_data: timestep  ---> ',numero_step

    if (numero_step == 1) then  
       write(errunit,*) '> data_dir: ''',trim(data_dir),''''
       ! contains the number of snapshots to analyze and their names, type and number (see below)
       call open_exist(12,'inputfiles_HaloMaker.dat')
    endif

    ! then read name of snapshot, its type (pm, p3m, SN, Nzo, Gd), num of procs used and number of snapshot
    read(12,*) name_of_file,type,nbPes,numstep
    write(file_num,'(i3.3)') numstep
    !!write(name_of_file,'(a,a)') trim(data_dir),trim(name_of_file)

    ! Note 1: old treecode SNAP format has to be converted [using SNAP_to_SIMPLE (on T3E)] 
    !         into new treecode SIMPLE (SN) format.
    ! Note 2: of the five format (pm, p3m, SN, Nzo, Gd) listed above, only  SN, Nzo and Gd 
    !         are fully tested so the code stops for pm and p3m

    if (numero_step == nsteps) close(12)

    if (type == 'Ra3') then

       call read_ramses_new(name_of_file)

       ! Computation of omega_t = omega_matter(t)
       !
       !                            omega_f*(1+z)^3
       ! omega(z)   = ------------------------------------------------------------------
       !              omega_f*(1+z)^3+(1-omega_f-omega_lambda_f)*(1+z)^2+omega_lambda_f
       !
       !
       !                              omega_lambda_0
       ! omega_L(z) = ----------------------------------------------------------------
       !              omega_f*(1+z)^3+(1-omega_f-omega_lambda_f)*(1+z)^2+omega_lambda_f

       omega_t  = omega_f*(af/aexp)**3
       omega_t  = omega_t/(omega_t+(1.-omega_f-omega_lambda_f)*(af/aexp)**2+omega_lambda_f)

     else if (type == 'Gd') then

       call read_gadget(name_of_file)

       ! Computation of omega_t = omega_matter(t)
       !
       !                            omega_f*(1+z)^3
       ! omega(z)   = ------------------------------------------------------------------
       !              omega_f*(1+z)^3+(1-omega_f-omega_lambda_f)*(1+z)^2+omega_lambda_f
       !
       !
       !                              omega_lambda_0
       ! omega_L(z) = ----------------------------------------------------------------
       !              omega_f*(1+z)^3+(1-omega_f-omega_lambda_f)*(1+z)^2+omega_lambda_f

       omega_t  = omega_f*(af/aexp)**3
       omega_t  = omega_t/(omega_t+(1.-omega_f-omega_lambda_f)*(af/aexp)**2+omega_lambda_f)
 
    else
       write(errunit,*) '> Don''t know the snapshot format: ''',type,''''
       stop
    endif

    write(errunit,*) '> min max position (in box units)   :',minval(pos),maxval(pos)
    write(errunit,*) '> min max velocities (in km/s)      :',minval(vel),maxval(vel)
    write(errunit,*) '> Reading done.'
    write(errunit,*) '> aexp = ',aexp

    return

  end subroutine read_data

!***********************************************************************
  subroutine open_exist(unit,name)

    implicit none

    integer(kind=4)                                  :: unit,lnblnk
    character(len=*)                                 :: name
    character(len=len_trim(name)+len_trim(data_dir)) :: file

    write(file,'(a,a)') data_dir(1:lnblnk(data_dir)),name(1:lnblnk(name))
    open(unit=unit,file=file,status='unknown')

    return

  end subroutine open_exist

!***********************************************************************
  subroutine open_append(unit,name)

    implicit none

    integer(kind=4)                                  :: unit,lnblnk
    character(len=*)                                 :: name    
    character(len=len_trim(name)+len_trim(data_dir)) :: file

    write(file,'(a,a)') data_dir(1:lnblnk(data_dir)),name(1:lnblnk(name))
    if (numero_step == 1) then 
       open(unit=unit,file=file,status='unknown')
    else
       open(unit=unit,file=file,status='old',position='append')
    end if

    return

  end subroutine open_append

!***********************************************************************
  subroutine read_gadget(name_file)


    ! Useful info for reading Gadget2
    ! Particles in gadget are sorted into 6 type
    ! The number of particles in each type is found in io_header%Vnpart_tot(1:6) and their masses in io_header%massarr(1:6)
    ! > When output is splitted the number of particles per file is io_header%Vnpart(1:6)
    ! > When output is splitted the number of particles is io_header%Vnpart(1:6) = io_header%Vnpart_tot(1:6) 

    ! > In most cases
    ! particles type 1 Gas
    ! particles type 2 Dark Matter
    ! particles type 5 Star
    ! nbodies = io_header%Vnpart_tot(2) 
    ! nbaryon = io_header%Vnpart_tot(1) + io_header%Vnpart_tot(5) 
    ! id of DM particles are between nbaryon - 1 and nbaryon + nbodies - 1

    ! > For simulation -DBIG_RUN option
    ! io_header%massarr(2) > 0., Gd_mass is not written in the file
    ! > For resimulation 
    ! io_header%massarr(2) = 0., Gd_mass is written in the file
   
    ! > If resimulation dark matter only ( Gadget 1 format )
    ! particles type 2 high resolution Dark Matter particles
    ! particles type 3 low  resolution Dark Matter particles
    ! particles type 4 low  resolution Dark Matter particles
    ! particles type 5 low  resolution Dark Matter particles
    ! nbodies = io_header%Vnpart_tot(2) + io_header%Vnpart_tot(3) + io_header%Vnpart_tot(4) + io_header%Vnpart_tot(5) 
    ! nhr     = io_header%Vnpart_tot(2)
    ! nlr     = io_header%Vnpart_tot(2) + io_header%Vnpart_tot(3) + io_header%Vnpart_tot(4) + io_header%Vnpart_tot(5) 
    ! nbaryon = io_header%Vnpart_tot(2) + io_header%Vnpart_tot(1) + io_header%Vnpart_tot(6) 
    ! All values of io_header%massarr are 0, all other values of io_header%Vnpart_tot should be 0
    ! id of DM particles are between 1 and nbodies 


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
!       character*60    :: fill
    end type header
    type(header)       :: io_header
    integer(kind=4)                      :: VN,nbaryon,id_part,idp,itype
    integer(kind=4),allocatable          :: Vid(:)
    real(kind=8),allocatable             :: Gd_pos(:,:),Gd_vel(:,:)    
#ifndef BIG_RUN
    real(kind=8),allocatable             :: Gd_mass(:)
    integer(kind=4)                      :: VN_lr,VN_total,nlr
    real(kind=8)                         :: HRsoftlen, LRsoftlen    
#endif
    real(kind=8)                         :: Lbox_pt_local
    ! softening lengths (comoving, and maximum physical in kpc/h) used in the simulation:
    real(kind=8),parameter               :: HRepscom = 10.0, HRepsmaxph = 5.0
    real(kind=8),parameter               :: LRepscom = 10.0, LRepsmaxph = 5.0
    integer(kind=4)                      :: i,ierr
    real(kind=8)                         :: Omega0,OmegaBaryon,mass_fac 
    character*200                        :: paramname,value,line,fileparam
    logical(kind=4)                      :: gadget1
#ifdef Test_Gd
    integer(kind=4)                      :: typeoffset,ngood,nbad,nrest
    integer(kind=4),allocatable          :: typeindex(:)
#ifndef BIG_RUN
    integer(kind=4)                      :: nvlr
    real(kind=8)                         :: massmin,massmax
    real(kind=8), allocatable            :: massloc(:)
#endif
#endif

    ! first build file name as path/snapshot_xxx.n, where:
    ! "path/snapshot_xxx" are in "name_file" variable, and 
    ! "n" is built by: //fnumber(verify(fnumber,' '):3)
    stop_read = 0
    nf        = -1
    nbig      = 1000
    nbodies   = 0
    nbaryon   = 0
    gadget1   = .false.

!********************************************************
    read_files: do while (nf <= nbig .and. stop_read == 0)
       istat = 1
       ! serch the file to read
       if(nf.lt.0) then
          ! try to open the single file if it exist
          Vfilename = trim(name_file)
          open(1,file=Vfilename,status='old',form='unformatted',iostat=istat)
          if(istat == 0) then
             stop_read = 1
             write(errunit,*) '> read_gadget... only one file to read'
          else
             ! .. if something went wrong (iostat/=0): 
             write(errunit,*) '> read_gadget... several files to read'
          end if
       end if
       if(istat /= 0) then
          nf        = nf+1
          write(fnumber,'(i3)') nf  !--> transform the integer nf in a character
          Vfilename = trim(name_file)// '.'//fnumber(verify(fnumber,' '):3)
          ! .. try to open file # (nf-1)
          ! possibly only several snap per step exist, thus .fnumber exists 
          ! (the file name is snapshot_xxx.fnumber instead of snapshot_xxx), so 
          ! try to read again, wit .fnumber: starting at zero if .1 isn't found either we suppose no file exist
          open(1,file=Vfilename,status='old',form='unformatted',iostat=istat)
          ! .. if something went wrong again ==> stop here:
          if (istat /= 0.and.nf.le.1) stop ' > ERROR: read_gadget: input file not there'
          if(istat /=0) stop_read = 1
       endif
       if (istat /= 0) exit read_files

       write(errunit,*) '> Read file: ',trim(Vfilename)

       ! .. if everything ok, continue:
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
       ! .. if file exists but is empty, stop reading
       if (istat /= 0) stop ' > ERROR: read_gadget... could not read header.'

       if(nbodies.le.0) then
          ! if nbodies is not allocated we're reading the first file
          
          write(errunit,*)
          write(errunit,*) '> Useful data from header:'
          write(errunit,*) '> io_header%time (aexp/af)          :',real(io_header%time,8)
          write(errunit,*) '> io_header%redshift                :',real(io_header%redshift,8)
          write(errunit,*) '> io_header%Boxsize                 :',real(io_header%Boxsize,8)
          write(errunit,*) '> io_header%Omega0                  :',real(io_header%Omega0,8)
          write(errunit,*) '> io_header%OmegaLambda             :',real(io_header%OmegaLambda,8)
          write(errunit,*) '> io_header%HubbleParam             :',real(io_header%HubbleParam,8)
          write(errunit,*)

#ifdef Test_Gd
          do itype = 1,6
             write(errunit,'(1x,a,i1,1x,a1,1x,i10)') '> Tot nb of particles type ',itype,':',io_header%Vnpart_tot(itype)
             write(errunit,'(1x,a,i1,1x,a1,1x,e10.3)') '> mass for particles type  ',itype,':',real(io_header%massarr(itype))
          end do
          write(errunit,*)
#endif
          nbodies = io_header%Vnpart_tot(2) 
          nbaryon = io_header%Vnpart_tot(1) + io_header%Vnpart_tot(5)
          if(io_header%Vnpart_tot(3).gt.0.or.io_header%Vnpart_tot(4).gt.0.or.io_header%Vnpart_tot(6).gt.0) then
             gadget1 = .true.
#ifdef BIG_RUN       
             write(errunit,*)
             write(errunit,*) '> WARNING: the code was compiled with -DBIG_RUN option'
             write(errunit,*) '> Some particles have not been accounted for.'
             do itype = 1,6
                write(errunit,'(1x,a,i1,1x,a1,1x,i10)') '> Tot nb of particles type ',itype,':',io_header%Vnpart_tot(itype)
                write(errunit,'(1x,a,i1,1x,a1,1x,e10.3)') '> mass for particles type  ',itype,':',real(io_header%massarr(itype))
             end do
             write(errunit,*) '> Make sure this is correct'
             write(errunit,*)
             if(io_header%Vnpart(2).gt.0.and.io_header%massarr(2).eq.0.) then
                write(errunit,*)
                write(errunit,*) '> ERROR: read_gadget... mass is 0. for DM particles (type2)'
                write(errunit,*) '> Code was compiled with -DBIG_RUN option'
                write(errunit,*) '> This shouldn''t be the case'
                stop
             end if
#else
             do itype = 1,6
                write(errunit,'(1x,a,i1,1x,a1,1x,i10)') '> Tot nb of particles type ',itype,':',io_header%Vnpart_tot(itype)
                write(errunit,'(1x,a,i1,1x,a1,1x,e10.3)') '> mass for particles type  ',itype,':',real(io_header%massarr(itype))
             end do
             write(errunit,*)
             write(errunit,*) '> WARNING: interpreting type 2 particles as high res DM particles '
             nhr = io_header%Vnpart_tot(2)
             if(io_header%Vnpart(1).gt.0) then
                write(errunit,*) '> WARNING: interpreting type 3, type 4 particles as low res DM particles'
                nlr = io_header%Vnpart_tot(3) + io_header%Vnpart_tot(4)
                write(errunit,*) '> WARNING: interpreting type 1, type 5 particles as baryonnic matter'
                write(errunit,*) 
             else
                write(errunit,*) '> WARNING: interpreting type 3, type 4,type 5 particles as low res DM particles'
                nlr = io_header%Vnpart_tot(3) + io_header%Vnpart_tot(4) + io_header%Vnpart_tot(5)
             endif
             if(io_header%Vnpart_tot(6).gt.0) write(errunit,*) '> WARNING: don''t know how to interpret type 6 particles'
             write(errunit,*) '> Rerun after compiling with -DTest_Gd option for testing'
             write(errunit,*) '> Nb of high resolution particles   :',nhr 
             write(errunit,*) '> Nb of low resolution particles    :',nlr
             nbodies = nhr + nlr
             nbaryon = io_header%Vnpart_tot(1) + io_header%Vnpart_tot(6)
#endif
          endif
          write(errunit,*) '> Number of DM particles            :',nbodies
          write(errunit,*) '> Number of baryons                 :',nbaryon
            
          allocate(pos(nbodies,3))
          allocate(vel(nbodies,3))
          pos     = 0.0
          vel     = 0.0
#ifndef BIG_RUN
          allocate(mass(nbodies), epsvect(nbodies))
          mass    = 0.0
          epsvect = 0.0
#endif
       else
          gadget1 = .false.
       end if

#ifdef Test_Gd
       if(nf.ge.0) then
          ! For now we don't know to what kind of particles type 3,4 and 6 do correspond to
          write(errunit,*)
          do itype = 1,6
             write(errunit,'(1x,a,i1,1x,a1,i10)') '> Nb of particles type ',itype,':',io_header%Vnpart(itype)
          end do
          write(errunit,*)
       end if
#endif

       VN    = sum(io_header%Vnpart)       
       write(errunit,*) '> Total nb of particles in file     :',sum(io_header%Vnpart)

       allocate(Gd_pos(3,VN),Gd_vel(3,VN),Vid(VN))
       read(1) Gd_pos
       read(1) Gd_vel
       read(1) Vid
#ifndef BIG_RUN
       allocate(Gd_mass(VN))
       if(gadget1.and.io_header%massarr(2).gt.0..and.io_header%Vnpart(1).eq.0) then
           read(1,iostat=istat) Gd_mass(io_header%Vnpart(2)+1:) !old version Hr masses were not written it seems.
       else
          read(1,iostat=istat) Gd_mass
       end if
       if(istat /= 0) stop ' > ERROR: read_gadget... masses are not written in this file.'
#endif
       close(1)

#ifdef Test_Gd
       ! We thing that DM part id are written with an offset corresponding to nbaryons
       write(errunit,'(1x,a,6(i10,1x))')   '> io_header%Vnpart :',io_header%Vnpart
       write(errunit,'(1x,a,6(E10.3,1x))') '> io_header%massarr:',io_header%massarr
       write(errunit,*)
       typeoffset = 0
       do itype = 1, 6
          write(errunit,'(1x,a,i1,a,i10)') '> type: ',itype,' io_header%Vnpart ',io_header%Vnpart(itype)
          if(io_header%Vnpart(itype).gt.0) then
             write(errunit,*) '> between:', typeoffset +1,'and', typeoffset+io_header%Vnpart(itype)
             ngood = 0
             nbad  = 0
             nrest = 0
             if(nbaryon .gt. 0) then
                do id_part = typeoffset+1, typeoffset+io_header%Vnpart(itype)
                   
                   if(Vid(id_part).ge.0.and.Vid(id_part).lt.nbaryon) then
                      ngood = ngood + 1
                   else if(Vid(id_part).lt.nbaryon + nbodies) then
                      nbad  = nbad + 1
                   else
                      nrest = nrest + 1
                   end if
                end do
                write(errunit,*) '> id between 0 and nbaryon-1               :',ngood
                write(errunit,*) '> id between nbaryon and nbaryon+nbodies-1 :',nbad
                write(errunit,*) '> id below 0 or above nbaryon+nbodies-1    :',nrest   
             else
                do id_part = typeoffset+1, typeoffset+io_header%Vnpart(itype)
                   if(Vid(id_part).gt.0.and.Vid(id_part).le.nbodies) then
                      ngood = ngood + 1
                   else if(Vid(id_part).eq.0) then
                      nbad  = nbad + 1
                   else
                      nrest = nrest + 1
                   end if
                end do
                write(errunit,*) '> id between 1 and nbodies                 :',ngood
                write(errunit,*) '> id between eq 0                          :',nbad
                write(errunit,*) '> id below 0 or above nbodies              :',nrest                   
             end if
             allocate(typeindex(io_header%Vnpart(itype)))
             typeindex(1:io_header%Vnpart(itype)) = Vid(typeoffset+1:typeoffset + io_header%Vnpart(itype))
             write(errunit,*) '> main,max Vid',minval(typeindex),maxval(typeindex)
             deallocate(typeindex)
#ifndef BIG_RUN
             allocate(massloc(io_header%Vnpart(itype)))
             massloc(1:io_header%Vnpart(itype)) = Gd_mass(typeoffset+1:typeoffset + io_header%Vnpart(itype))
             write(errunit,*) '> main,max Gd_mass',minval(massloc),maxval(massloc)
             deallocate(massloc)
#endif
          end if
          typeoffset = typeoffset + io_header%Vnpart(itype)
          write(errunit,*)
       end do
#endif

       ! if SPH simulation, you have gas particles (number in io_header%Vnpart(1))
       ! which can turn into star particles (in io_header%Vnpart(5)) so we have to skip
       ! these to read the DM particle data only
           
       do id_part = io_header%Vnpart(1)+1,io_header%Vnpart(1)+io_header%Vnpart(2) ! loop only on DM particles 
          idp   = Vid(id_part) - nbaryon
          if (idp == 0) idp = nbodies ! offset from C to FORTRAN numbering
          if (idp < 1 .or. idp > nbodies ) then
             write(errunit,*) '> idp,nbodies',idp,nbodies
             stop ' > ERROR: read_gadget... bad particle id'
          end if
          pos(idp,:) = Gd_pos(:,id_part)     
          vel(idp,:) = Gd_vel(:,id_part)
#ifndef BIG_RUN
          mass(idp)  = Gd_mass(id_part)
#endif
       end do

#ifndef BIG_RUN
       if(gadget1) then
          !loop for low rest of particles
          do id_part = io_header%Vnpart(1)+io_header%Vnpart(2) +1,&
               io_header%Vnpart(1)+io_header%Vnpart(2) +io_header%Vnpart(3)+io_header%Vnpart(4)
             idp   = Vid(id_part) - nbaryon
             if (idp == 0) idp = nbodies ! offset from C to FORTRAN numbering
             if (idp < 1 .or. idp > nbodies ) then
                write(errunit,*) '> idp,nbodies',idp,nbodies
                stop ' > ERROR: read_gadget... bad particle id'
             end if
             pos(idp,:) = Gd_pos(:,id_part)     
             vel(idp,:) = Gd_vel(:,id_part)
             mass(idp)  = Gd_mass(id_part)
          end do
          if(io_header%Vnpart_tot(1).eq.0) then
             do id_part =   io_header%Vnpart(1)+io_header%Vnpart(2) +io_header%Vnpart(3)+io_header%Vnpart(4)+1,&
                  io_header%Vnpart(1)+io_header%Vnpart(2) +io_header%Vnpart(3)+io_header%Vnpart(4)+io_header%Vnpart(5)
                idp   = Vid(id_part) - nbaryon
                if (idp == 0) idp = nbodies ! offset from C to FORTRAN numbering
                if (idp < 1 .or. idp > nbodies ) then
                   write(errunit,*) '> idp,nbodies',idp,nbodies
                   stop ' > ERROR: read_gadget... bad particle id'
                end if

                pos(idp,:) = Gd_pos(:,id_part)     
                vel(idp,:) = Gd_vel(:,id_part)
                mass(idp)  = Gd_mass(id_part)
             end do
          end if
       end if
       deallocate(Gd_mass)
#endif
       deallocate(Gd_pos, Gd_vel, Vid)
       
    enddo read_files
!********************************************************

    aexp  = io_header%time*af  ! because SN format assumes a(t_in)=1
!    write(errunit,*) '  SN_aexp, z                                   : ',aexp,-1.+af/aexp 
#ifdef BIG_RUN
    ! have to recalculate mboxp because Gadget particles are in physical units 
    mboxp = io_header%massarr(2)*10./H_f*real(nbodies,8) ! in 10^11 M_sun
    massp = io_header%massarr(2)*10./H_f/mboxp  !in units of mboxp
#else
    ! count low res and high res DM particles
    massp = minval(mass)
    nhr     = 0
    nlr     = 0
    do idp = 1, nbodies
       if(mass(idp).eq.massp) then
          nhr = nhr + 1
       else if(mass(idp).gt.massp) then
          nlr = nlr + 1
       end if
    end do
    write(errunit,*) '> Nb of high resolution particles   :',nhr 
    write(errunit,*) '> Nb of low resolution particles    :',nlr
    
    ! .. check that reading went ok for resimulation
    VN_total = nhr+nlr
    Vn_lr    = nlr
    if (VN_total /= nbodies) then
       stop ' > read_gadget... ERROR: nlr + nhr /= nbodies !' 
    endif
    mboxp = real(sum(real(mass,8)),8)*10./H_f ! in 10^11 M_sun
    massp = massp*10./H_f/mboxp  !in units of mboxp
#endif
    
    write(errunit,*) '> particle mass (in M_sun)          :', massp*mboxp*1e11
 

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

    pos           = pos*io_header%time/(io_header%Boxsize*aexp) - 0.5
    Lbox_pt_local = Lboxp*(aexp/ai)
    vel           = vel *sqrt(io_header%time)
#ifndef BIG_RUN
    ! Physical softening length at the current time, used in the original simulation
    HRsoftlen  = min(real(io_header%time,8)*HRepscom,HRepsmaxph)
    LRsoftlen  = min(real(io_header%time,8)*LRepscom,LRepsmaxph)
    ! get epsvect = softening length in units of "fraction of Lbox(t)", as needed by subroutine softgrav:
    epsvect(1:nhr)  = HRsoftlen/(10.*H_f*Lbox_pt_local)
    epsvect(nhr+1:) = LRsoftlen/(10.*H_f*Lbox_pt_local)
    !write(errunit,*) ' '
    !write(errunit,*) 'softening length of HR parts in units of Lbox(t)=', epsvect(1:1)
    !write(errunit,*) 'softening length of LR parts in units of Lbox(t)=', epsvect(nhr+1:nhr+1)
    ! .. Transform masses from 10^10 Msol/h --> to units of mboxp
    mass          = mass * 10.0 / H_f / mboxp !in units of box total mass
    write(errunit,*) '> HR particle mass (in M_sun)       :', minval(mass)*mboxp*1e11
    ! .. min and max mass:
    !write(errunit,*) 
    minlrmrat = minval(mass(nhr+1:))/massp
    !write(errunit,*) '  min-max MASS of LR particles [in units of mhr]  : ', &
    !                    minlrmrat,maxval(mass(nhr+1:))/massp
#endif  

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
          write(errunit,*) '> paremeters I neaded form the "parameter-usedvalues" file'
          write(errunit,*) '> Omega0       : ',Omega0
          write(errunit,*) '> OmegaBaryon  : ',OmegaBaryon
          if (Omega0 .le. 0.) stop ' > read_gadget... Coudn''t read Omega0 value'
          if (OmegaBaryon .le. 0.) stop ' > read_gadget... Coudn''t read OmegaBaryon value'
          write(errunit,*)
       end if
       mass_fac = 1./(1.0-OmegaBaryon/Omega0)
       mboxp    = mboxp * mass_fac
#ifndef RENORM
       massp = massp/mass_fac
#ifndef BIG_RUN
       mass  = mass/mass_fac
#endif
#endif
#ifdef RENORM
       write(errunit,*) '> mass after renormalization (M_sun):', massp*mboxp*1e11
#endif
    endif

    return

  end subroutine read_gadget


!***********************************************************************
subroutine read_ramses_new(repository)

  ! This routine reads star particles dumped in the RAMSES 3.0 new i/o format.

    implicit none

    character(len=*)            :: repository
!!$    integer(kind=8)             :: npart
    integer(kind=4)             :: ndim,idim,icpu,ipos,ncpu,i,ipar,ichem
    integer(kind=4)             :: ncpu2,npart2,ndim2,idum,nout,nsink,nstar
    integer(kind=4)             :: nx,ny,nz,nlevelmax,ngridmax,nstep_coarse
    integer(kind=4),allocatable :: idp(:),llp(:)
    real(kind=8)                :: boxlen,tco,aexp_ram,hexp
    real(kind=8)                :: omega_m,omega_l,omega_k,omega_b
    real(kind=8)                :: scale_l,scale_d,scale_t,dummy
    real(kind=8)                :: mtot,massres
    real(kind=8),allocatable    :: dumout(:),tmpp(:,:),tmpv(:,:),tmpm(:),tmpt(:),tmpz(:),tmpchem(:,:)
    real(kind=8),allocatable    :: tmppf(:),tmppz(:) ! RS - space for prist fraction & primord Z
    ! cannot do otherwise than setting all the strings to a value larger than that of a line 
    ! in the input file and trim them whenever it is needed
    character(len=200)          :: line,name,value,nomfich
    character(len=5)            :: nchar,ncharcpu,nchargrp
    logical                     :: ok
    ! Particle types
    integer(1),parameter :: FAM_DM=1, FAM_STAR=2, FAM_CLOUD=3, FAM_DEBRIS=4, FAM_OTHER=5, FAM_UNDEF=127
    integer(1),parameter :: FAM_TRACER_GAS=0
    integer(1),parameter :: FAM_TRACER_DM=-1, FAM_TRACER_STAR=-2, FAM_TRACER_CLOUD=-3, FAM_TRACER_DEBRIS=-4, FAM_TRACER_OTHER=-5
    integer(1),allocatable :: famp(:), tagp(:)
    
    ! NB: repository is directory containing output files
    ! e.g. /horizon1/teyssier/ramses_simu/boxlen100_n256/output_00001/

    ! read cosmological params in header of amr file
    ipos    = index(repository,'output_')
    nchar   = repository(ipos+7:ipos+13)
    if (iogroupsizerep>0) then
       nomfich = trim(repository)//'group_00001/amr_'//trim(nchar)//'.out00001'
    else
       nomfich = trim(repository)//'amr_'//trim(nchar)//'.out00001'
    end if
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
       read(10) dummy
       close(10)
       !write(errunit,991) ncpu,ndim,nstep_coarse
       !write(errunit,993) nlevelmax,ngridmax
       !write(errunit,997) boxlen*scale_l
       !write(errunit,994) t,aexp,hexp
       !write(errunit,995) omega_m,omega_l,omega_k,omega_b
    end if
991 format(' ncpu=',i6,' ndim=',i1,' nstep=',i6)
993 format(' nlevelmax=',i3,' ngridmax=',i8)
994 format(' t=',1pe10.3,' aexp=',1pe10.3,' hexp=',1pe10.3)
995 format(' omega_m=',F6.3,' omega_l=',F6.3,' omega_k=',F6.3,' omega_b=',F6.3)
997 format(' boxlen=',1pe10.3,' h-1 Mpc')

    if (iogroupsizerep>0) then
       nomfich = trim(repository)//'group_00001/info_'//trim(nchar)//'.txt'
    else
       nomfich = trim(repository)//'info_'//trim(nchar)//'.txt'
    end if
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
          !       write(errunit,*) trim(name),' : ',trim(value)
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

    !write(errunit,*) 'scale_l,scale_d,scale_t: ',scale_l,scale_d,scale_t
    Lboxp          = boxlen*scale_l/3.08e24/aexp_ram ! converts cgs to Mpc comoving
    !write(errunit,*) 'af,hf,lboxp,ai,aexp: ',af,h_f,lboxp,ai,aexp_ram
    aexp           = aexp_ram*af  
    omega_f        = omega_m
    omega_lambda_f = omega_l
    omega_c_f      = omega_k

    ! now read the particle data files
    if (iogroupsizerep>0) then
       nomfich = trim(repository)//'/group_00001/part_'//trim(nchar)//'.out00001'
    else
       nomfich = trim(repository)//'/part_'//trim(nchar)//'.out00001'
    end if
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
       if(iogroupsizerep>0) call title(((icpu-1)/iogroupsizerep)+1, nchargrp)

       if (iogroupsizerep>0) then
          nomfich = trim(repository)//'/group_'//TRIM(nchargrp)//'/part_'//trim(nchar)//'.out'//trim(ncharcpu)
       else
          nomfich = trim(repository)//'/part_'//trim(nchar)//'.out'//trim(ncharcpu)
       end if

       !write(errunit,*)'> Reading file '//trim(nomfich)
       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1) ncpu2
       read(1) ndim2
       read(1) npart2
       read(1) idum
       read(1) nstar
       read(1) dummy   
       read(1) dummy
       read(1) nsink
       close(1)
       npart = npart+npart2
       !write(errunit,*) 'nstar',icpu,npart,npart2,nstar
    end do
    
    npart   = npart - nstar - nsink*2109
    nbodies = nstar
    write(errunit,*)'> Found ',npart,' DM   particles'
    write(errunit,*)'> Found ',nstar,' star particles'
    write(errunit,*)'> Reading positions, velocities and masses...'
    
    allocate(pos(1:nstar,1:ndim))
    allocate(vel(1:nstar,1:ndim))
    allocate(mass(1:nstar))
    allocate(age_st(1:nstar))
    allocate(famp(1:npart), tagp(1:npart))
#ifdef METALS
    allocate(met_st(1:nstar))
    allocate(pf_st(1:nstar)) ! RS - allocate space for prist fraction of star particle
    allocate(pz_st(1:nstar)) ! RS - allocate space for primordial Z of star particle
    if(nchem.gt.0) allocate(chem_st(1:nstar,1:nchem))
#endif

    do icpu = 1,ncpu
       call title(icpu,ncharcpu)
       if(iogroupsizerep>0) call title(((icpu-1)/iogroupsizerep)+1, nchargrp)
       if (iogroupsizerep>0) then
          nomfich = trim(repository)//'/group_'//TRIM(nchargrp)//'/part_'//trim(nchar)//'.out'//trim(ncharcpu)
       else
          nomfich = trim(repository)//'/part_'//trim(nchar)//'.out'//trim(ncharcpu)
       end if
       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1) ncpu2
       read(1) ndim2
       read(1) npart2
       read(1) idum
       read(1) nstar
       read(1) dummy
       read(1) dummy
       read(1) nsink
       !write(errunit,*)npart2,nstar,nsink
       allocate(tmpp(1:npart2,1:ndim),tmpv(1:npart2,1:ndim),tmpm(1:npart2) &
            & ,tmpt(1:npart2),idp(1:npart2),llp(1:npart2))
#ifdef METALS
       allocate(tmpz(1:npart2))
       allocate(tmppf(1:npart2))  ! RS - prist fract of star particle
       allocate(tmppz(1:npart2))  ! RS - primod Z of star particle
       if(nchem.gt.0)allocate(tmpchem(1:npart2,1:nchem))
#endif          
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
       ! read grid level of particles
       read(1) llp(1:npart2)

       ! particle types and tags
       read(1) famp(1:npart2)
       read(1) tagp(1:npart2)


       ! read all particle creation times if necessary
       tmpt = 0.0
       ! read info relative to star particles
       if ((nstar .gt. 0) .or. (nsink .gt. 0)) then 
          read(1) tmpt(1:npart2)
#ifdef METALS
          read(1) tmpz(1:npart2)
          read(1) tmppf(1:npart2)  ! RS - Prist frac follows Z
          read(1) tmppz(1:npart2)  ! RS - Primordial Z follows PF
          do ichem=1,nchem
             read(1) tmpchem(1:npart2,ichem)
          enddo
#endif          
       endif
       close(1)

       !write(errunit,*) 'min max, id',minval(idp),maxval(idp)
       ! now get rid of non-star particles
       !!do ipar=1,npart2
       !!   if (idp(ipar) .gt. 0.0 .and. tmpt(ipar) .ne. 0.0)idp(ipar)=-idp(ipar)
       !!enddo
       do ipar=1,npart2 
          !!! for star particles
          ! if (idp(ipar) .lt. 0.0 .and. tmpt(ipar) .ne. 0.0) then 
          if (famp(ipar) == FAM_STAR) then
             ! put all positions between -0.5 and 0.5
             pos(abs(idp(ipar)),1:ndim) = real(tmpp(ipar,1:ndim),8) - 0.5
             ! convert code units to km/s 
             vel(abs(idp(ipar)),1:ndim) = real(tmpv(ipar,1:ndim),8)*scale_l/scale_t*1e-5
             mass(abs(idp(ipar)))       = real(tmpm(ipar),8)
             if(dump_stars)then
                age_st(abs(idp(ipar)))  = real(tmpt(ipar),8)
#ifdef METALS
                met_st(abs(idp(ipar)))  = real(tmpz(ipar),8)
                pf_st(abs(idp(ipar)))  = real(tmppf(ipar),8) ! RS - adding prist fraction for stars
                pz_st(abs(idp(ipar)))  = real(tmppz(ipar),8) ! RS - adding primordial Z for stars
                do ichem=1,nchem
                   chem_st(abs(idp(ipar)),ichem)  = real(tmpchem(ipar,ichem),8)
                enddo
#endif
             endif
          end if
       end do
       deallocate(tmpp,tmpv,tmpm,tmpt,idp,llp)
#ifdef METALS
       deallocate(tmpz,tmppf,tmppz)
       if(nchem.gt.0)deallocate(tmpchem)
#endif          
    end do
    
    mtot = 0.0d0
    do i = 1,nstar
       mtot = mtot + real(mass(i),8)
    enddo
    ! that is for the dark matter so let's add baryons now if there are any 
    ! and renormalization flag is on !!
    massres = minval(mass)*mboxp*1d11
    massp   = minval(mass)
    write(errunit,*) '> particle mass (in M_sun)               = ',massres
#ifdef RENORM
    massres = minval(mass)*mboxp*1d11/mtot
    massp   = minval(mass)/real(mtot,8)
    write(errunit,*) '> particle mass (in M_sun) after renorm  = ',massres
#endif
#ifdef BIG_RUN
    deallocate(mass)
#endif

  end subroutine read_ramses_new


!***********************************************************************
  subroutine title(n,nchar)

    implicit none

    integer(kind=4) :: n
    character*5     :: nchar
    character*1     :: nchar1
    character*2     :: nchar2
    character*3     :: nchar3
    character*4     :: nchar4
    character*5     :: nchar5


    if(n.ge.10000)then
       write(nchar5,'(i5)') n
       nchar = nchar5
    elseif(n.ge.1000)then
       write(nchar4,'(i4)') n
       nchar = '0'//nchar4
    elseif(n.ge.100)then
       write(nchar3,'(i3)') n
       nchar = '00'//nchar3
    elseif(n.ge.10)then
       write(nchar2,'(i2)') n
       nchar = '000'//nchar2
    else
       write(nchar1,'(i1)') n
       nchar = '0000'//nchar1
    endif

  end subroutine title

!***********************************************************************
  subroutine title2(n,nchar)

    implicit none

    integer(kind=4) :: n
    character*7     :: nchar
    character*1     :: nchar1
    character*2     :: nchar2
    character*3     :: nchar3
    character*4     :: nchar4
    character*5     :: nchar5
    character*6     :: nchar6
    character*7     :: nchar7


    if(n.ge.1000000)then
       write(nchar7,'(i7)') n
       nchar = nchar7
    elseif(n.ge.100000)then
       write(nchar6,'(i6)') n
       nchar = '0'//nchar6
    elseif(n.ge.10000)then
       write(nchar5,'(i5)') n
       nchar = '00'//nchar5
    elseif(n.ge.1000)then
       write(nchar4,'(i4)') n
       nchar = '000'//nchar4
    elseif(n.ge.100)then
       write(nchar3,'(i3)') n
       nchar = '0000'//nchar3
    elseif(n.ge.10)then
       write(nchar2,'(i2)') n
       nchar = '00000'//nchar2
    else
       write(nchar1,'(i1)') n
       nchar = '000000'//nchar1
    endif

  end subroutine title2
  

!***********************************************************************
  subroutine write_tree_brick

  ! This subroutine writes the information relevant to building a halo 
  ! merging tree (using the build_tree program) i.e. for each halo:
  !   1/ the list of all the particles it contains (this enables us --- as  
  !      particle numbers are time independent --- to follow the halo history) 
  !   2/ its properties which are independent of merging history (mass ...)
 
    implicit none

    integer(kind=4)                                         :: i,unitfile,start,j,idim,ndim=3,ichem
    character(LEN=5)                                        :: nchar
    character(LEN=7)                                        :: ncharg
    character(LEN=300)                                      :: nomfich
#ifndef BIG_RUN
    character(len=len_trim(data_dir)+16)                    :: file
#endif
    character(len=len_trim(data_dir)+len_trim(file_num)+11) :: filename
    integer(kind=4) ,allocatable                            :: members(:)
    real(kind=8) ,allocatable                            :: mass_memb(:),age_memb(:),met_memb(:),mdump(:)
    real(kind=8) ,allocatable                            :: pf_memb(:), pz_memb(:) ! RS - prist frac & primord Z
    real(kind=8) ,allocatable                            :: pos_memb(:,:),vel_memb(:,:),chem_memb(:,:)
    logical                                              :: done
    
    if(dump_stars)then
       nchar   = '00'//TRIM(file_num)
       call system('mkdir GAL_'//TRIM(nchar))
    endif

    done = .false.
#ifndef BIG_RUN
    if (write_resim_masses) then 
       write(file,'(a,a)') trim(data_dir),'resim_masses.dat'
       unitfile = 44
       open(unitfile,file=file,form='unformatted',status='unknown')
       write(unitfile) nbodies
       write(unitfile) mass
       close(unitfile)     
       write_resim_masses = .false.
    end if
#endif

    if(.not.fsub) then
       write(filename,'(a,a,a3)') trim(data_dir),'tree_brick_',file_num
    else
       write(filename,'(a,a,a3)') trim(data_dir),'tree_bricks',file_num
    end if
    unitfile = 44
    write(errunit,*)
    write(errunit,*) '> Output data to build halo merger tree to: ',filename
    open(unit=unitfile,file=filename,form='unformatted',status='unknown')
    write(unitfile) nbodies
    write(unitfile) massp
    write(unitfile) aexp
    write(unitfile) omega_t
    write(unitfile) age_univ
    write(unitfile) nb_of_halos, nb_of_subhalos
    do i=1,nb_of_halos + nb_of_subhalos
       ! write list of particles in each halo
       allocate(members(1:nb_of_parts(i)))
       if(dump_stars)then
          allocate(mass_memb(1:nb_of_parts(i)),pos_memb(1:nb_of_parts(i),1:3),vel_memb(1:nb_of_parts(i),1:3) &
               &,age_memb(1:nb_of_parts(i)),mdump(1:nb_of_parts(i)))
#ifdef METALS
          allocate(met_memb(1:nb_of_parts(i)))
          allocate(pf_memb(1:nb_of_parts(i)))  ! RS - space for prist fraction
          allocate(pz_memb(1:nb_of_parts(i)))  ! RS - space for primordial Z
          if(nchem.gt.0)allocate(chem_memb(1:nb_of_parts(i),1:nchem))
#endif
       endif
       start = first_part(i)
       do j =1,nb_of_parts(i)            
          members(j) = start
          if(dump_stars)then
             mass_memb (j)=mass (start)
             pos_memb(j,1)=pos(start,1)
             pos_memb(j,2)=pos(start,2)
             pos_memb(j,3)=pos(start,3)
             vel_memb(j,1)=vel(start,1)
             vel_memb(j,2)=vel(start,2)
             vel_memb(j,3)=vel(start,3)
             age_memb  (j)=age_st(start)
#ifdef METALS
             met_memb(j)=met_st(start)
             pf_memb(j)=pf_st(start)  ! RS - Adding prist fraction
             pz_memb(j)=pz_st(start)  ! RS - Adding primord Z
             do ichem=1,nchem
                chem_memb(j,ichem)=chem_st(start,ichem)
             enddo
#endif
          endif
          start = linked_list(start)
       end do
       write(unitfile) nb_of_parts(i)
       write(unitfile) members

       if(dump_stars)then
           call title2(liste_halos(i)%my_number,ncharg)
           nomfich='GAL_'//TRIM(nchar)//'/gal_stars_'//TRIM(ncharg)
           open(unit=9,file=nomfich,form='unformatted')
           write(9)liste_halos(i)%my_number
           write(9)liste_halos(i)%level
           write(9)dble(liste_halos(i)%m)
           write(9)dble(liste_halos(i)%p%x),dble(liste_halos(i)%p%y),dble(liste_halos(i)%p%z)
           write(9)dble(liste_halos(i)%v%x),dble(liste_halos(i)%v%y),dble(liste_halos(i)%v%z)
           write(9)dble(liste_halos(i)%L%x),dble(liste_halos(i)%L%y),dble(liste_halos(i)%L%z)
           write(9)nb_of_parts(i)
           do idim=1,ndim
              mdump(1:nb_of_parts(i))=pos_memb(1:nb_of_parts(i),idim)
              write(9)dble(mdump)
           enddo
           do idim=1,ndim
              mdump(1:nb_of_parts(i))=vel_memb(1:nb_of_parts(i),idim)
              write(9)dble(mdump)
           enddo
           write(9)dble(mass_memb)
           write(9)members
           write(9)dble(age_memb)
#ifdef METALS
           write(9)dble(met_memb)
         !   write(9)dble(pf_memb)  ! RS - write the prist fraction
         !   write(9)dble(pz_memb)  ! RS - write the primordial Z
           do ichem=1,nchem
              mdump(1:nb_of_parts(i))=chem_memb(1:nb_of_parts(i),ichem)
              write(9)dble(mdump)
           enddo
#endif
           close(9)
           deallocate(mass_memb,pos_memb,vel_memb,age_memb,mdump)
#ifdef METALS
           deallocate(met_memb, pf_memb, pz_memb) ! RS - deallocate prist fraction & primord Z
           if(nchem.gt.0)deallocate(chem_memb)
#endif
       endif

       deallocate(members)
       ! write each halo properties
       call write_halo(liste_halos(i),unitfile)
    enddo

    close(unitfile)

    return

  end subroutine write_tree_brick

!***********************************************************************
  subroutine write_halo(h,unitfile)

    implicit none

    integer(kind=4) :: unitfile
    type (halo)     :: h

    ! Masses (h%m,h%datas%mvir) are in units of 10^11 Msol, and 
    ! Lengths (h%p%x,h%p%y,h%p%z,h%r,h%datas%rvir) are in units of Mpc
    ! Velocities (h%v%x,h%v%y,h%v%z,h%datas%cvel) are in km/s
    ! Energies (h%ek,h%ep,h%et) are in
    ! Temperatures (h%datas%tvir) are in K
    ! Angular Momentum (h%L%x,h%L%y,h%L%z) are in
    ! Other quantities are dimensionless (h%my_number,h%my_timestep,h%spin)  

    write(unitfile) h%my_number
    write(unitfile) h%my_timestep 
    write(unitfile) h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub
    write(unitfile) h%m
    write(unitfile) h%p%x,h%p%y,h%p%z
    write(unitfile) h%v%x,h%v%y,h%v%z
    write(unitfile) h%L%x,h%L%y,h%L%z 
    write(unitfile) h%r, h%sh%a, h%sh%b, h%sh%c
    write(unitfile) h%ek,h%ep,h%et
    write(unitfile) h%spin
    write(unitfile) h%sigma,h%sigma_bulge,h%m_bulge
    write(unitfile)h%datas%rvir,h%datas%mvir,h%datas%tvir,&
        &h%datas%cvel,h%datas%Reff
    write(unitfile) h%halo_profile%rho_0,h%halo_profile%r_c
    ! Stellar surface density profiles
    write(unitfile) h%nbin
    write(unitfile) h%rr
    write(unitfile) h%rho

    return

  end subroutine write_halo

  !***********************************************************************
  subroutine read_last_brick(filelast) 

    implicit none
    integer(kind=4) :: idummy,nh_old,nsub_old,i,npstruct,j
    real(kind=8)    :: rdummy
    integer(kind=4),allocatable :: members(:) 

    character(len=14) ::filelast
    character(len=(len_trim(data_dir)+len_trim(filelast)))   :: filename

    ex_liste_parts = 0

    write(filename,'(a,a)') trim(data_dir),trim(filelast)
    write(errunit,*) 'Reading tree_brick from previous running:',trim(filename)
    
    open(unit=30,file=filename,form='unformatted',status='unknown')
    read(30) idummy                                         ! nbodies
    if(verbose) write(errunit,*) 'nbodies  :',idummy
    read(30) rdummy                                         ! massp
    if(verbose) write(errunit,*) 'massp    :',rdummy
    read(30) rdummy                                         ! aexp
    if(verbose) write(errunit,*) 'aexp     :',rdummy
    read(30) rdummy                                         ! omega_t
    if(verbose) write(errunit,*) 'omega_t  :',rdummy
    read(30) rdummy                                         ! age_univ
    if(verbose) write(errunit,*) 'age_univ :',rdummy
    read(30) nh_old,nsub_old                                ! nb_of_halos, nb_of_subhalos
    ex_nb_of_structs = nh_old + nsub_old
    allocate(ex_level(ex_nb_of_structs)) !! 
    do i=1,ex_nb_of_structs
       read(30)  npstruct ! nb_of_parts(i)
       !ex_nb_of_parts(i) = npstruct
         ! read list of particles in each halo
       allocate(members(npstruct))
       read(30) members
       do j =1,npstruct           
          ex_liste_parts(members(j)) = i
       end do
       deallocate(members)
       ! write each halo properties
       call read_halo(30,i)         
    enddo

    close(30)
    write(errunit,*)

    return

  end subroutine read_last_brick

  !***********************************************************************
  subroutine read_halo(unitfile,ih) 

    implicit none
    integer(kind=4) :: idummy,unitfile,ih
    real(kind=8)    :: rdummy

    read(unitfile) idummy                               !h%my_number
    read(unitfile) idummy                               !h%my_timestep 
    read(unitfile) ex_level(ih),idummy,idummy,idummy,idummy   !h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub
    read(unitfile) rdummy                               !h%m
    read(unitfile) rdummy,rdummy,rdummy                 !h%p%x,h%p%y,h%p%z
    read(unitfile) rdummy,rdummy,rdummy                 !h%v%x,h%v%y,h%v%z
    read(unitfile) rdummy,rdummy,rdummy                 !h%L%x,h%L%y,h%L%z 
    read(unitfile) rdummy,rdummy,rdummy,rdummy          !h%r, h%sh%a, h%sh%b, h%sh%c
    read(unitfile) rdummy,rdummy,rdummy                 !h%ek,h%ep,h%et
    read(unitfile) rdummy                               !h%spin
    read(unitfile) rdummy,rdummy,rdummy,rdummy          !h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel
    read(unitfile) rdummy,rdummy                        !h%halo_profile%rho_0,h%halo_profile%r_c

    return

  end subroutine read_halo

!***********************************************************************
!///////////////////////////////////////////////////////////////////////

end module input_output

