program TreeMaker

  ! this program reads in halo/particle data obtained after running the HaloMaker code 
  ! on N-Body simulation snapshots and builds a halo merger tree from these data.

  use compute_tree
  use input_output
#ifdef SIMPL
  use simpl_merger_tree
  use jeje
#endif

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
  write(errunit,*) ' mod BUSHID       :        J. Blaizot                (2008)            ' 
  write(errunit,*) ' mod accretion    :        J. Blaizot                (2009)            '
  write(errunit,*) ' option SIMPL     :        D. Tweed                  (2010)            '   
  write(errunit,*)
  write(errunit,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'  
  write(errunit,*) 

  call getarg(1,data_dir)       ! get directory where to input/output the data
  
  call init                     ! initialize cosmological and technical parameters
  
  open(unit=444,file='accretionHistory.dat',status='unknown',form='formatted')
  numero_step = 0
  do st_do    = 1,nsteps_do     ! loop over snapshots to determine halo fathers and sons
     call new_step
  end do
  close(444)

#ifdef SIMPL
  call simplify_merger_tree
#endif

  ! define BushID for all haloes in the simulation
  call define_bushes
  
  ! output the trees jeje's way
#ifdef SIMPL 
  call init_TsTree
  call define_IDs
  call buffered_outputs
#else
  if (n_tree_files < 1) then ! We can get one tree_file 
     call write_single_tree
  else
     call write_bushes
     !call write_one_bush(126685)
  end if
#endif

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






















