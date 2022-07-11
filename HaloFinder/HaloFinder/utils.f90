module utils

  ! include general routines which have to be moved here because used by subbox and compute_halo_props

  use halo_defs

  ! Some variable defined for cosmological setup
  ! conformal time things
  integer(kind=4),parameter             :: n_frw = 1000
  real(KIND=8),dimension(:),allocatable :: aexp_frw,hexp_frw,tau_frw,t_frw


  public

contains

!***********************************************************************

!!$  subroutine correct_for_periodicity(dr)
!!$
!!$  ! subroutine corrects for the fact that if you have periodic boundary conditions,
!!$  ! then groups of particles that sit on the edge of the box can have part of their
!!$  ! particles on one side and part on the other. So we have to take out a box-length
!!$  ! when measuring the distances between group members if needed.
!!$
!!$    implicit none
!!$
!!$    type (vector) :: dr
!!$
!!$    if (FlagPeriod == 0) return  !--> NO PERIODIC BCs
!!$
!!$    if (dr%x > + Lbox_pt2) dr%x = dr%x - Lbox_pt
!!$    if (dr%x <= - Lbox_pt2) dr%x = dr%x + Lbox_pt
!!$
!!$    if (dr%y > + Lbox_pt2) dr%y = dr%y - Lbox_pt
!!$    if (dr%y <= - Lbox_pt2) dr%y = dr%y + Lbox_pt
!!$
!!$    if (dr%z > + Lbox_pt2) dr%z = dr%z - Lbox_pt
!!$    if (dr%z <= - Lbox_pt2) dr%z = dr%z + Lbox_pt
!!$
!!$    return
!!$
!!$  end subroutine correct_for_periodicity

!***********************************************************************

  subroutine correct_for_periodicity_code_units(x)

    ! same as correct_for_periodicity but argument is only one coord in code units (i.e. from -0.5 to 0.5)

    implicit none

    real(kind=4) :: x

    if (x >  0.5d0) then
       x = x - 1.0d0
    end if
    if (x <= -0.5d0) then
       x = x + 1.0d0
    end if

    return

  end subroutine correct_for_periodicity_code_units

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

  ! conformal time utils :
  ! Most of these are extracted from RASCAS, but date back earlier!
    function ct_conftime2time(tau)

    ! return look-back time in yr

    implicit none

    real(kind=8),intent(in) :: tau
    real(kind=8)            :: ct_conftime2time
    integer(kind=4)         :: i


    ! locate bracketing conf. times
    i = 1
    do while(tau_frw(i) > tau .and. i < n_frw)
       i = i + 1
    end do
    ! Interploate time
    ct_conftime2time = t_frw(i) * (tau-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &
         & t_frw(i-1)        * (tau-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))

    return

  end function ct_conftime2time

  function ct_proptime2time(tau,h0)
    ! return look-back time in yr

    implicit none
    real(kind=8),intent(in) :: tau,h0
    real(kind=8)            :: ct_proptime2time

    ct_proptime2time = tau / (h0 / 3.08d19) / (365.25*24.*3600.)
    return
  end function ct_proptime2time


  function ct_aexp2time(a)

    ! return look-back time in yr

    implicit none

    real(kind=8),intent(in) :: a
    real(kind=8)            :: ct_aexp2time
    integer(kind=4)         :: i

    ! find bracketting aexp's
     i = 1
     do while(aexp_frw(i)>a.and.i<n_frw)
        i = i + 1
     end do
     ! Interploate time
     ct_aexp2time = t_frw(i) * (a-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
          & t_frw(i-1)    * (a-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))

    return

  end function ct_aexp2time


  subroutine ct_init_cosmo(omega_m,omega_l,omega_k,h0)

    ! h0 is in km/s/Mpc

    implicit none
    real(kind=8),intent(in) :: omega_m,omega_l,omega_k,h0
    real(kind=8)            :: time_tot

    allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
    allocate(tau_frw(0:n_frw),t_frw(0:n_frw))
    call ct_friedman(omega_m,omega_l,omega_k,1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)
    ! convert time to yr
    t_frw = t_frw / (h0 / 3.08d19) / (365.25*24.*3600.)

    return

  end subroutine ct_init_cosmo


  subroutine ct_clear_cosmo

    implicit none

    deallocate(aexp_frw,hexp_frw,tau_frw,t_frw)

    return

  end subroutine ct_clear_cosmo


  subroutine ct_friedman(O_mat_0,O_vac_0,O_k_0,alpha_prec,axp_min, &
       & axp_out,hexp_out,tau_out,t_out,ntable,age_tot)

    implicit none
    integer::ntable
    real(kind=8)::O_mat_0, O_vac_0, O_k_0
    real(kind=8)::alpha_prec,axp_min,age_tot
    real(kind=8),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
    ! ######################################################!
    ! This subroutine assumes that axp = 1 at z = 0 (today) !
    ! and that t and tau = 0 at z = 0 (today).              !
    ! axp is the expansion factor, hexp the Hubble constant !
    ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
    ! time, and t the look-back time, both in unit of 1/H0. !
    ! alpha_prec is the required accuracy and axp_min is the     !
    ! starting expansion factor of the look-up table.       !
    ! ntable is the required size of the look-up table.     !
    ! ######################################################!
    real(kind=8)::axp_tau, axp_t
    real(kind=8)::axp_tau_pre, axp_t_pre
    real(kind=8)::dtau,dt
    real(kind=8)::tau,t
    integer::nstep,nout,nskip

    !  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
    !     write(*,*)'Error: non-physical cosmological constants'
    !     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
    !     write(*,*)'The sum must be equal to 1.0, but '
    !     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
    !     stop
    !  end if

    axp_tau = 1.0D0
    axp_t = 1.0D0
    tau = 0.0D0
    t = 0.0D0
    nstep = 0

    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) )

       nstep = nstep + 1
       dtau = alpha_prec * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha_prec * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt

    end do

    age_tot=-t

    nskip=nstep/ntable

    axp_t = 1.d0
    t = 0.d0
    axp_tau = 1.d0
    tau = 0.d0
    nstep = 0
    nout=0
    t_out(nout)=t
    tau_out(nout)=tau
    axp_out(nout)=axp_tau
    hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) )

       nstep = nstep + 1
       dtau = alpha_prec * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha_prec * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt

       if(mod(nstep,nskip)==0)then
          nout=nout+1
          t_out(nout)=t
          tau_out(nout)=tau
          axp_out(nout)=axp_tau
          hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
       end if

    end do
    t_out(ntable)=t
    tau_out(ntable)=tau
    axp_out(ntable)=axp_tau
    hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  contains
    function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
      implicit none
      real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
      dadtau = axp_tau*axp_tau*axp_tau *  &
           &   ( O_mat_0 + &
           &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
           &     O_k_0   * axp_tau )
      dadtau = sqrt(dadtau)
      return
    end function dadtau

    function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
      implicit none
      real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
      dadt   = (1.0D0/axp_t)* &
           &   ( O_mat_0 + &
           &     O_vac_0 * axp_t*axp_t*axp_t + &
           &     O_k_0   * axp_t )
      dadt = sqrt(dadt)
      return
    end function dadt

  end subroutine ct_friedman

end module utils



