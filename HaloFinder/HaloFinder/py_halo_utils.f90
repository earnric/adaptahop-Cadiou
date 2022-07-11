module py_halo_utils

public

contains

  
  subroutine get_nb_halos(file,nhalo,nsubhalo)
    
    implicit none 

    character(1000),intent(in)  :: file
    integer(kind=4),intent(out) :: nhalo
    integer(kind=4),intent(out) :: nsubhalo

    open(unit=11,file=file,status='old',form='unformatted',action='read')
    read(11) !nbodies_full_sim
    read(11) !massp
    read(11) !aexp
    read(11) !omega_t
    read(11) !age_univ
    read(11) nhalo,nsubhalo
    close(11)

    return

  end subroutine get_nb_halos

  subroutine read_all_halos_with_contam(file,halos,n)

    implicit none 

    character(1000),intent(in) :: file
    integer(kind=4),intent(in) :: n
    real(kind=4),intent(out)   :: halos(n,33)
    
    integer(kind=4) :: i,nh,ns
    integer(kind=4) :: np,num,ts,lev,hosth,hosts,nbsub,nextsub,contam
    real(kind=4)    :: mass,x,y,z,vx,vy,vz,lx,ly,lz,r,a,b,c,ek,ep,et,spin,rvir,mvir,tvir,cvel,rho_0,r_c

    open(unit=11,file=file,status='old',form='unformatted')
    read(11) !nbodies_full_sim
    read(11) !massp
    read(11) !aexp
    read(11) !omega_t
    read(11) !age_univ
    read(11) nh,ns
    if (n > nh+ns) stop
    do i = 1,n
       read(11) np
       read(11) ! skip part id's
       read(11) num
       read(11) ts
       read(11) lev,hosth,hosts,nbsub,nextsub
       read(11) mass
       read(11) x,y,z
       read(11) vx,vy,vz
       read(11) lx,ly,lz
       read(11) r,a,b,c
       read(11) ek,ep,et
       read(11) spin
       read(11) rvir,mvir,tvir,cvel
       read(11) rho_0,r_c
       read(11) contam
       halos(i,:) = (/real(np,4),real(num,4),real(ts,4),real(lev,4),real(hosth,4),real(hosts,4),real(nbsub,4),&
            & real(nextsub,4),mass,x,y,z,vx,vy,vz,lx,ly,lz,r,a,b,c,&
            & ek,ep,et,spin,rvir,mvir,tvir,cvel,rho_0,r_c,real(contam,4)/)
    end do
    close(11)

    return
    
  end subroutine read_all_halos_with_contam
  
  subroutine read_all_halos(file,halos,n)

    implicit none 

    character(1000),intent(in) :: file
    integer(kind=4),intent(in) :: n
    real(kind=4),intent(out)   :: halos(n,32)
    
    integer(kind=4) :: i,nh,ns
    integer(kind=4) :: np,num,ts,lev,hosth,hosts,nbsub,nextsub,contam
    real(kind=4)    :: mass,x,y,z,vx,vy,vz,lx,ly,lz,r,a,b,c,ek,ep,et,spin,rvir,mvir,tvir,cvel,rho_0,r_c

    open(unit=11,file=file,status='old',form='unformatted')
    read(11) !nbodies_full_sim
    read(11) !massp
    read(11) !aexp
    read(11) !omega_t
    read(11) !age_univ
    read(11) nh,ns
    if (n > nh+ns) stop
    do i = 1,n
       read(11) np
       read(11) ! skip part id's
       read(11) num
       read(11) ts
       read(11) lev,hosth,hosts,nbsub,nextsub
       read(11) mass
       read(11) x,y,z
       read(11) vx,vy,vz
       read(11) lx,ly,lz
       read(11) r,a,b,c
       read(11) ek,ep,et
       read(11) spin
       read(11) rvir,mvir,tvir,cvel
       read(11) rho_0,r_c
       halos(i,:) = (/real(np,4),real(num,4),real(ts,4),real(lev,4),real(hosth,4),real(hosts,4),real(nbsub,4),&
            & real(nextsub,4),mass,x,y,z,vx,vy,vz,lx,ly,lz,r,a,b,c,&
            & ek,ep,et,spin,rvir,mvir,tvir,cvel,rho_0,r_c/)
    end do
    close(11)

    return

  end subroutine read_all_halos

  subroutine read_all_galaxies(file,halos,n)

    ! Galaxy maker version of read_all_halos

    implicit none 

    character(1000),intent(in) :: file
    integer(kind=4),intent(in) :: n
    real(kind=4),intent(out)   :: halos(n,32)
    
    integer(kind=4) :: i,nh,ns
    integer(kind=4) :: np,num,ts,lev,hosth,hosts,nbsub,nextsub,contam
    real(kind=4)    :: mass,x,y,z,vx,vy,vz,lx,ly,lz,r,a,b,c,ek,ep,et,spin,rvir,mvir,tvir,cvel,rho_0,r_c
    !real(kind=4)    :: sm(:), sa(:), sZ(:) ! stellar particle mass, age and metallicity

    open(unit=11,file=file,status='old',form='unformatted')
    read(11) !nbodies_full_sim
    read(11) !massp
    read(11) !aexp
    read(11) !omega_t
    read(11) !age_univ
    read(11) nh,ns
    if (n > nh+ns) stop
    do i = 1,n
       read(11) np
       read(11) ! skip stellar particle mass
       read(11) ! skip stellar particle age
       read(11) ! skip stellar particle metallicity
       read(11) ! skip part id's
       read(11) num
       read(11) ts
       read(11) lev,hosth,hosts,nbsub,nextsub
       read(11) mass
       read(11) x,y,z
       read(11) vx,vy,vz
       read(11) lx,ly,lz
       read(11) r,a,b,c
       read(11) ek,ep,et
       read(11) spin
       read(11) rvir,mvir,tvir,cvel
       read(11) rho_0,r_c
       halos(i,:) = (/real(np,4),real(num,4),real(ts,4),real(lev,4),real(hosth,4),real(hosts,4),real(nbsub,4),&
            & real(nextsub,4),mass,x,y,z,vx,vy,vz,lx,ly,lz,r,a,b,c,&
            & ek,ep,et,spin,rvir,mvir,tvir,cvel,rho_0,r_c/)
    end do
    close(11)

    return

  end subroutine read_all_galaxies


end module py_halo_utils
