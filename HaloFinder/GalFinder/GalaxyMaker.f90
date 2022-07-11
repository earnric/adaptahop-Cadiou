program GalaxyMaker

  ! this program reads in the particle data from N-Body simulation snapshots and 
  ! puts them into halos (using a friend-of-friend algorithm). It then goes on to
  ! compute various properties for these halos (mass, spin, energy ...)

  use compute_halo_props
  use halo_defs

  implicit none

  write(errunit,*)
  write(errunit,*) '_______________________________________________________________________'
  write(errunit,*)
  write(errunit,*) '              GalaxyMaker                                                '
  write(errunit,*) '              ---------                                                ' 
  write(errunit,*)
  write(errunit,*)
  write(errunit,*) ' first version    :        S. Ninin                  (1999)            '
  write(errunit,*) ' modif.           :        J. Devriendt              (1999-2002)       '
  write(errunit,*) ' modif.           :        B. Lanzoni & S. Hatton    (2000-2001)       '
  write(errunit,*) ' galics v1.0      :        J. Blaizot & J. Devriendt (2002)            '
  write(errunit,*) ' horizon v1.0     :        J. Devriendt & D. Tweed   (2006)            '
  write(errunit,*) ' horizon v2.0     :        J. Devriendt & D. Tweed   (2007)            '
  write(errunit,*) ' horizon v2.0.2   :        J. Devriendt, D. Tweed & J.Blaizot (2008)   '
  write(errunit,*) ' parallel version :        J. Blaizot (2010)                           '
  write(errunit,*) 
  write(errunit,*) '_______________________________________________________________________'  
  write(errunit,*) 

  ! get directory where to input/output the data
  call getarg(1,data_dir)
  ! initialize cosmological and technical parameters of the N_Body simulation 
  call init

  ! loop over snapshots
  do numero_step=1,nsteps
     call new_step
     write(errunit,*)numero_step,'npart(step)=',npart
  end do

  write(errunit,*)
  write(errunit,*) '> Bricks now available to build halo merger tree with TreeMaker'
  write(errunit,*)
  write(errunit,*) '_______________________________________________________________________'  
  write(errunit,*)
  write(errunit,*) '      End of GalaxyMaker                                        '
  write(errunit,*)
  write(errunit,*) '_______________________________________________________________________'  


  stop

end program GalaxyMaker







