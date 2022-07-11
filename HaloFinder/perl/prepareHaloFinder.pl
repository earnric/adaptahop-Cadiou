#!/usr/bin/perl
use strict;
use warnings;


# define snapshot directory (dir containing all output_xxxxx directories)
my $snapDir = "/home/pfister/data/TDE/TDE4";
my $NumberOfSnap = 51;
my $TypeofRun =  ""; #can be "_BR" or ""

# define (and create if needed) global output directory (which will host many sub-directories...)
my $runDir = "$snapDir"."/Halos/";
my $runDirGalaxies = "$snapDir"."/Galaxies/";
if (!-e $runDir) {system("mkdir $runDir");}
if (!-e $runDirGalaxies) {system("mkdir $runDirGalaxies");}
system("rm HaloMaker*");

# define list of snapshots to process
my @snapList = (1);
for (my $k = 0; $k<$NumberOfSnap; $k++) {$snapList[$k] = $k+1;}

# copy code over to run directory, as a documentation ... 
if (!-e $runDir."/code") {system("mkdir $runDir"."/code");}
if (!-e $runDirGalaxies."/code") {system("mkdir $runDirGalaxies"."/code");}
system("cp -r ../\* $runDirGalaxies"."/code/.");
system("cp -r ../\* $runDir"."/code/.");

# define executable file as the above copied executable (so we can link with no risk). 
my $execFile  = "$runDir"."/code/HaloFinder/HaloFinder".$TypeofRun;
my $execFileGalaxies  = "$runDirGalaxies"."/code/GalFinder/GalFinder";
my $cmdFile   = "$runDir" . "cmdFile.sh";    # file which will contain instructions to run HaloFinder
my $cmdFileGalaxies   = "$runDirGalaxies" . "cmdFile.sh";    # file which will contain instructions to run HaloFinder


# LOOP ON SNAPSHOTS 
foreach (@snapList)
  {
    # convert snapshot number into i5.5 format (e.g. 00001)
    my $i   = $_;
    my $num = "$i";
    if ($i < 10000) {$num = "0"."$num";}
    if ($i < 1000)  {$num = "0"."$num";}
    if ($i < 100)   {$num = "0"."$num";}
    if ($i < 10)    {$num = "0"."$num";}

    # define location of simulation files (snapshots)
    my $fileRoot  = "$snapDir" . "/Outputs/output_$num/";

    # create working directory for snapshot
    my $dir = "$runDir"."$i";
    my $cmd = "mkdir $dir";
    system($cmd);
    
    # write input_HaloMaker.dat file
    write_input_hm($dir, 0);

    # link executable and stuff to running dir.
    $cmd = "ln -s $execFile $dir"; 
    system($cmd);

    # write inputfiles_HaloMaker.dat file
    my $file = "$dir"."/inputfiles_HaloMaker.dat";
    open(JFILE,">$file");
    print JFILE "\'$fileRoot \'  Ra3   1  $num \n";
    close(JFILE);
  }

foreach (@snapList)
  {
    # convert snapshot number into i5.5 format (e.g. 00001)
    my $i   = $_;
    my $num = "$i";
    if ($i < 10000) {$num = "0"."$num";}
    if ($i < 1000)  {$num = "0"."$num";}
    if ($i < 100)   {$num = "0"."$num";}
    if ($i < 10)    {$num = "0"."$num";}

    # define location of simulation files (snapshots)
    my $fileRoot  = "$snapDir" . "/Outputs/output_$num/";

    # create working directory for snapshot
    my $dirGalaxies = "$runDirGalaxies"."$i";
    my $cmd = "mkdir $dirGalaxies";
    system($cmd);
    
    # write input_HaloMaker.dat file
    #write_input_hm($dir, true);
    write_input_hm($dirGalaxies, 1);

    # link executable and stuff to running dir.
    $cmd = "ln -s $execFileGalaxies $dirGalaxies"; 
    system($cmd);

    # write inputfiles_HaloMaker.dat file
    my $file = "$dirGalaxies"."/inputfiles_HaloMaker.dat";
    open(JFILE,">$file");
    print JFILE "\'$fileRoot \'  Ra3   1  $num \n";
    close(JFILE);
  }

# write job file 
open(JFILE,">Halo.job");
print JFILE "#!/bin/sh \n";
print JFILE "#PBS -S /bin/sh \n";
print JFILE "#PBS -j oe \n";
print JFILE "#PBS -N HaloMaker\n";
print JFILE "#PBS -l nodes=b6:ppn=5,walltime=48:00:00 \n";
print JFILE " \n";
print JFILE "module (){\n";
print JFILE "    eval \$(\/usr\/bin\/modulecmd bash \$\*)\n";
print JFILE "\}\n";
print JFILE "module load openmpi/1.8.8-ifort-15.0 \n";
print JFILE "export PATH=\"\/home\/pfister\/anaconda2\/bin:\$PATH\" \n";
print JFILE " \n";
print JFILE "cd ".$runDir."\n";
print JFILE "j=PBS_ARRAYID\n";
print JFILE "cd \$j\n";
print JFILE "./HaloFinder".$TypeofRun." 2> log.out\n";
print JFILE " if [ \$j = 1 ]\n";
print JFILE " then\n";
print JFILE " ../code/HaloFinder/snap2mass > snap2mass.log\n";
print JFILE " cp ~/Codes/HaloFinder/BH.py ".$snapDir."\n";
print JFILE " cp ~/Codes/HaloFinder/loadSnap.py ".$snapDir."\n";
#print JFILE " cd ../..\n";
#print JFILE " python BH.py\n";
print JFILE " fi\n";
print JFILE " \n";
print JFILE "cd ".$runDirGalaxies."\n";
print JFILE "cd \$j\n";
print JFILE "./GalFinder".$TypeofRun." 2> log.out\n";
print JFILE " if [ \$j = $NumberOfSnap ]\n";
print JFILE " then\n";
print JFILE " ../code/GalFinder/snap2mass > snap2mass.log\n";
print JFILE " mv resim_masses.dat ../1\n";
print JFILE " fi\n";
print JFILE "exit 0\n";
print JFILE " \n";
close(JFILE);

# write job file 
open(JFILE,">Finder.sh");
print JFILE "for i in \{1..".$NumberOfSnap."\}\n";
print JFILE "do\n";
print JFILE "     sed -i s/HaloMaker/HaloMaker_\$i/g Halo.job\n";
print JFILE "     sed -i s/j=PBS_ARRAYID/j=\$i/g Halo.job\n";
print JFILE "     qsub Halo.job\n";
print JFILE "     sed -i s/HaloMaker_\$i/HaloMaker/g Halo.job\n";
print JFILE "     sed -i s/j=\$i/j=PBS_ARRAYID/g Halo.job\n";
print JFILE "done\n";
print JFILE " \n";
close(JFILE);


sub write_input_hm {
    my $filename  = "$_[0]"."/input_HaloMaker.dat";
    open(FILE,">$filename");
    # MN PARAMETERS 
    print FILE "lbox            = 59.0493080328 ! Size of the box at z=0 (in Mpc! ==> 479 Mpc/h=684.286 Mpc) \n";
    #print FILE "lbox            = 142 ! Size of the box at z=0 (in Mpc! ==> 479 Mpc/h=684.286 Mpc) \n";
    #print FILE "lbox            = 20 ! Size of the box at z=0 (in Mpc! ==> 479 Mpc/h=684.286 Mpc) \n";
    print FILE "npart           = 50           ! Minimum number of particles per halo (ex:10) \n";
    print FILE "af              = 1.0        	! Today expansion factor == 1+zi \n";
    print FILE "H_f             = 67.7399978638 ! Hubble constant (in Km/s/Mpc) at z=0 \n";
    print FILE "omega_f         = 0.308899998665! Omega matter at z=0 (ex: 0.3) \n";
    print FILE "lambda_f        = 0.691100001335! Omega lambda at z=0 (ex: 0.7) \n"; 
    print FILE "SC              = .true.        ! Find the center using shrinking sphere \n"; 
    print FILE "eps_SC              = 0.1        ! Fraction by which the radius is reduced \n"; 
    print FILE "dcell_min       = 0.001         ! Minium radius of the SC in Mpc\n";
    print FILE "fudgepsilon     = 0.0001          ! parameter for adaptahop (usually 0.05) YOHAN 1/2**(lmax(redshift)-3)\n";   
    print FILE "FlagPeriod      = 1           	! 1 for periodic boundary conditions 0 else \n\n";
    # Zoom options : 
    #print FILE "xmin = -0.2 \n";
    #print FILE "xmax =  0.2 \n";
    #print FILE "ymin = -0.2 \n";
    #print FILE "ymax =  0.2 \n";
    #print FILE "zmin = -0.2 \n";
    #print FILE "zmax =  0.2 \n\n";
    # Back to MN PARAMETERS
    print FILE "method          = MSM        	! choose between FOF, HOP, DPM, MHM (view ReadMe for details). \n";
    print FILE "cdm             = .false.    	! false if not fof \n";
    print FILE "b               = 0.2        	! fof parameter b (usually b=0.2) \n";
    print FILE "nsteps          = 1         	! Number of time steps to analyse \n";

    # adaptaHOP parameters
    print FILE "nvoisins        = 20        	! parameter for adaptahop (usually 20)\n";
    print FILE "nhop            = 20         	! parameter for adaptahop (usually 20)\n";
    print FILE "rhot            = 80.       	! parameter for adaptahop (80 coreespond to b = 0.2)\n";
    print FILE "fudge           = 4.          	! parameter for adaptahop (usually 4.) \n";
    print FILE "alphap          = 1          	! parameter for adaptahop (usually 1.)\n";
    print FILE "verbose         = .false.     	! verbose parameter for both halo finder\n";
    print FILE "megaverbose     = .false.     	! parameter for adaptahop\n";
    if ($_[1]) {print FILE "dump_stars     = .true. !Dump stellar properties \n";}
    close(FILE);

}


