#!/user/bin/envpython3
"""
Author:GustavoTamasco
Scripttorun:

Thescripttakes:

"""

#importstatements
from sys import argv
import os.path
import subprocess
import os
from os import listdir
from os.path import isfile,join
import shutil
import re

#functions and classes

def create_genomes_dir(path):
    '''
    This function creates a new directory in the local machine were the fna files will
    be stored
    :param path: path to the new dir
    :return: None
    '''
    if os.path.exists(path):
        print("Dir {} was already created".format(path))
    else:
        os.makedirs(path)
        print("{} was created".format(path))

def move_files(path, new_path, files):
    '''
    This function will move the fna files to the created folder

    :param bof:
    :param new_dir:
    :return: none
    '''
    for file in files:
        if ".tsv" in file:
            current_path = "{}/{}".format(path, file)
            noval_path = "{}/{}".format(new_path, file)
            shutil.move(current_path, noval_path)
   


def run_vfdb(genome):
    '''
    This function runs parsnp
    :param dir_path: is the location in the start of the section
    :param genomes_dir: path to the fna files
    :return:
    '''
    out_name = genome.split(".")[0]
    out_name = "{}_vfdb.tsv".format(out_name)
    
    if os.path.exists(out_name):
        print("VFDB was already executed for {}".format(genome))
    else:
        cmd_abricate = "abricate --db vfdb {} > {}".format(genome, out_name)
        exit_message = subprocess.check_call(cmd_abricate, shell = True)
        print("Exit status: {0}".format(exit_message))
        print("{0} Parsnp was executed".format(out_name))

    return(out_name)


def run_argannot(genome):
    '''
    This function runs parsnp
    :param dir_path: is the location in the start of the section
    :param genomes_dir: path to the fna files
    :return:
    '''
    out_name = genome.split(".")[0]
    out_name = "{}_argannot.tsv".format(out_name)
    if os.path.exists(out_name):
        print("Argannot was already executed for {}".format(genome))
    else:
        cmd_abricate = "abricate --db argannot {} > {}".format(genome, out_name)
        exit_message = subprocess.check_call(cmd_abricate, shell = True)
        print("Exit status: {0}".format(exit_message))
        print("{0} Parsnp was executed".format(out_name))

    return(out_name)

def run_card(genome):
    '''
    This function runs parsnp
    :param dir_path: is the location in the start of the section
    :param genomes_dir: path to the fna files
    :return:
    '''
    out_name = genome.split(".")[0]
    out_name = "{}_card.tsv".format(out_name)
  
    
    if os.path.exists(out_name):
        print("CARD was already executed for {}".format(genome))
    else:
        cmd_abricate = "abricate --db card {} > {}".format(genome, out_name)
        exit_message = subprocess.check_call(cmd_abricate, shell = True)
        print("Exit status: {0}".format(exit_message))
        print("{0} Parsnp was executed".format(out_name))

    return(out_name)

def run_resfinder(genome):
    '''
    This function runs parsnp
    :param dir_path: is the location in the start of the section
    :param genomes_dir: path to the fna files
    :return:
    '''
    out_name = genome.split(".")[0]
    out_name = "{}_resfinder.tsv".format(out_name)
    
    
    if os.path.exists(out_name):
        print("Resfinder was already executed for {}".format(genome))
    else:
        cmd_abricate = "abricate --db resfinder {} > {}".format(genome, out_name)
        exit_message = subprocess.check_call(cmd_abricate, shell = True)
        print("Exit status: {0}".format(exit_message))
        print("{0} Parsnp was executed".format(out_name))
    return(out_name)

def list_files(path):
    '''
    This function lists all file in a the path
    :param path: path of a dir
    :return: all files path in a list
    '''
    files = [file for file in listdir(path) if isfile(join(path,file))]
    return (files)

def main():
    """Main code of the script"""

    #Get the files
    
    dirpath=os.getcwd()
    files = list_files(dirpath)
    genomes = [file for file in files if ".fna" in file]
    for genome in genomes:
        genome_index = genome.split(".")[0]
        if os.path.exists(genome_index):
            print("Genome {} was already evaluated".format(genome_index))
        else:
            genome_path = "{}/{}".format(dirpath,genome.split(".")[0])
            create_genomes_dir(genome_path)
            agr_file = run_argannot(genome)
            vfdb_file = run_vfdb(genome)
            card_file = run_card(genome)
            resfinder_file = run_resfinder(genome)
            genome_files = list_files(dirpath)
            genome_files = [x for x in genome_files if genome_index in x]
            move_files(dirpath, genome_path, genome_files)

       
        
   
    #os.chdir(path_to_all_info)
    #directories = list_directories(path_to_all_info)





#main
if __name__=='__main__':
    main()
