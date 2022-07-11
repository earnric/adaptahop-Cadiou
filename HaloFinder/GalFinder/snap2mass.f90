program snap2mass

  use input_output
  use compute_halo_props

  implicit none

  call getarg(1,data_dir) ! get directory where to input/output the data

  call init

  numero_step=1

  call read_data

  open(unit=233,file='resim_masses.dat',form='unformatted',status='unknown')
  write(233) nbodies
  write(233) mass
  close(233)
  
  stop

end program snap2mass


