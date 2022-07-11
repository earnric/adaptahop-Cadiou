program TreeMaker

  ! this program reads in halo/particle data obtained after running the HaloMaker code 
  ! on N-Body simulation snapshots and builds a halo merger tree from these data.

  use compute_tree
  use input_output

  implicit none
  
  write(errunit,*)
  write(errunit,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
  write(errunit,*)
  write(errunit,*) '      TREE_MAKER                                                       '
  write(errunit,*) '      ----------                                                       ' 
  write(errunit,*)
  write(errunit,*)
  write(errunit,*) ' first version    :        S. Ninin                  (1999)            '
  write(errunit,*) ' modif.           :        J. Devriendt              (1999-2002)       '
  write(errunit,*) ' modif.           :        B. Lanzoni & S. Hatton    (2000-2001)       '
  write(errunit,*) ' galics v1.0      :        J. Blaizot & J. Devriendt (2002)            '
  write(errunit,*) ' horizon v2.0     :        J. Devriendt & D.Tweed    (2006)            '
  write(errunit,*) ' horizon v2.1     :        J. Blaizot & J. Devriendt (2006)            '
  write(errunit,*) ' horizon v2.2     :        J. Devriendt & D.Tweed    (2007)            '
  write(errunit,*) '                  option BUSHID : J.Blaizot (2008)                     ' 
  write(errunit,*) '                  option COMPUTE_ACCRETION : J.B. (2009)               '
  write(errunit,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'  
  write(errunit,*) 

  call getarg(1,data_dir)       ! get directory where to input/output the data
  
  call init                     ! initialize cosmological and technical parameters

  open(unit=444,file='accretionHistory.dat',status='unknown',form='formatted')
  do numero_step = 1,nsteps     ! loop over snapshots to determine halo fathers and sons
     call new_step
  end do
  close(444)

  ! define BushID for all haloes in the simulation
  call define_bushes
  
  ! output the trees
  if (n_tree_files <= 1) then 
     call write_single_tree
  else
     call write_bushes
     !call write_one_bush(126685)
  end if

  call finish_all               ! properly deallocate all global variables before exiting the program
 
  write(errunit,*)
  write(errunit,*) '> DM halo tree now available to build galaxy tree with GalaxyMaker'
  write(errunit,*)
  write(errunit,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
  write(errunit,*)
  write(errunit,*) '      END OF TREE_MAKER                                                '
  write(errunit,*)
  write(errunit,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'  
  write(errunit,*) 

  stop

end program TreeMaker






















