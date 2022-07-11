#!/usr/bin/env python2                                                                                                                                                                                               
# -*- coding: utf-8 -*-                                                                                                                                                                                              
import os, sys, shutil

"""
What we want to do
------------------

 1. Create a Trees directory
 2. Link all the tree_bricks files
 3. In the Trees directory, copy the "code" directory
 4. In the Trees directory, create the right input_TreeMaker.dat file

Some configuration is needed:
 * The RATS directory shoud be specified, either as a environment variable RATS
   or directly in the configuration dictionary "config". The config dict is
   always preferred.
 * The folder in which the simulation lies is specified in the config dict.

JB: modif to have a sorted list of files ... 

"""


def main(config):
    import glob


    # First, create the Trees directory if it does not exist
    treedir = config['folder'] + '/Trees'
    os.system('rm -rf '+treedir)

    try: 
        os.makedirs(treedir)
    except OSError:
        if not os.path.isdir(treedir):
            raise

    # Then, move to the Trees directory
    os.chdir(treedir)

    # Then, link all the tree_bricks files to the Trees directory
    # FIXME: is this really necessary?
    #JB- 
    #fname = config['folder']+'/Halos/*/tree_bricks*'
    #tb_list = glob.glob(fname)
    tb_list = []
    for i in range(5000):
        fname = config['folder']+'/Halos/%i/tree_bricks%.3i'%(i,i)
        if os.path.isfile(fname):
            tb_list.append(fname)

    for tree_brick in tb_list:
        tbname = os.path.basename(tree_brick)
        os.symlink(tree_brick, tbname)

    # Then, get the "code" directory (only the f90 part)
    ratsdir = config.get('rats', os.getenv('RATS'))

    try:
        os.path.isdir(ratsdir+'/TreeMaker/TreeHalo/')
    except TypeError:
        print('ERROR: rats directory is not ', ratsdir)
        print('\nExiting...\n')
        sys.exit('RATS directory not found')
    else:
        codedir = ratsdir+'/TreeMaker/TreeHalo/'

    # Copy the f90 directory to "code"
    shutil.copytree(codedir, treedir+'/code/')

    # Link the TreeMaker binary to the current directory
    # TODO: add the TreeMaker_BR option for big runs
    if config.get('bigrun', True):
        os.symlink(treedir+'/code/TreeMaker_BR', 'TreeMaker_BR')
    else:
        os.symlink(treedir+'/code/TreeMaker', 'TreeMaker')
        os.symlink(config['folder']+'/Halos/1/resim_masses.dat', 'resim_masses.dat')

    # Now create the input_TreeMaker.dat file
    ifile = open('input_TreeMaker.dat', 'w')
    ifile.write(str(len(tb_list)) + '\t1\n')
    for tree_brick in tb_list:
        ifile.write('\''+tree_brick+'\'\n')
    ifile.close()

    #create a script to run treemaker
    os.chdir(ratsdir+'/TreeMaker/py')
    os.system('rm TreeMaker.sh')
    ifile = open('TreeMaker.sh', 'w')
    ifile.write('cd '+treedir+'\n')
    if config.get('bigrun',True):
        ifile.write('./TreeMaker_BR\n')
    else:
        ifile.write('./TreeMaker \n')
    ifile.close()



    # First, create the Trees directory if it does not exist
    treedir = config['folder'] + '/TreeStars'
    os.system('rm -rf '+treedir)

    try: 
        os.makedirs(treedir)
    except OSError:
        if not os.path.isdir(treedir):
            raise

    # Then, move to the Trees directory
    os.chdir(treedir)

    # Then, link all the tree_bricks files to the Trees directory
    # FIXME: is this really necessary?
    #JB- 
    #fname = config['folder']+'/Halos/*/tree_bricks*'
    #tb_list = glob.glob(fname)
    tb_list = []
    for i in range(5000):
        fname = config['folder']+'/Galaxies/%i/tree_bricks%.3i'%(i,i)
        if os.path.isfile(fname):
            tb_list.append(fname)

    for tree_brick in tb_list:
        tbname = os.path.basename(tree_brick)
        os.symlink(tree_brick, tbname)

    # Then, get the "code" directory (only the f90 part)
    ratsdir = config.get('rats', os.getenv('RATS'))

    try:
        os.path.isdir(ratsdir+'/TreeMaker/TreeStars/')
    except TypeError:
        print('ERROR: rats directory is not ', ratsdir)
        print('\nExiting...\n')
        sys.exit('RATS directory not found')
    else:
        codedir = ratsdir+'/TreeMaker/TreeStars/'

    # Copy the f90 directory to "code"
    shutil.copytree(codedir, treedir+'/code/')

    # Link the TreeMaker binary to the current directory
    # TODO: add the TreeMaker_BR option for big runs
    os.symlink(treedir+'/code/TreeMaker', 'TreeMaker')
    os.symlink(config['folder']+'/Galaxies/1/resim_masses.dat', 'resim_masses.dat')

    # Now create the input_TreeMaker.dat file
    ifile = open('input_TreeMaker.dat', 'w')
    ifile.write(str(len(tb_list)) + '\t1\n')
    for tree_brick in tb_list:
        ifile.write('\''+tree_brick+'\'\n')
    ifile.close()

    #create a script to run treemaker
    os.chdir(ratsdir+'/TreeMaker/py')
    ifile = open('TreeMaker.sh', 'a')
    ifile.write('cd '+treedir+'\n')
    ifile.write('./TreeMaker\n')
    ifile.close()

    print("Everything is done")
    
    return

if __name__ == '__main__':

    config = dict()
    config['folder'] = "/home/pfister/data/TDE/TDE4"
    config['bigrun'] = False
    config['rats'] = "/home/pfister/Codes"
    main(config)



