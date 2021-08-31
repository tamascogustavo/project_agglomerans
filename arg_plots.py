#!/user/bin/envpython3
"""
Author:GustavoTamasco
Scripttorun:

"""

#importstatements
from sys import argv
import os.path
import subprocess
import os
from os import listdir
from os.path import isfile,join
import shutil
from collections import Counter
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt





#fnctionsandclasses
def list_directories(path):
    '''
    This function lists all file in a the path
    :param path:path of a dir
    :return: allfiles path in a list
    '''
    all_dirs = listdir(path)
    return all_dirs

def list_files_simple(path):
    '''
    This function lists all file in a the path
    :param path: path of a dir
    :return: all files path in a list
    '''
    files = []
    files_ = [file for file in listdir(path) if isfile(join(path,file))]
    for f in files_:
        if ".tsv" in f:
            files.append(f)
    return (files)


def list_files(all_paths, path,organism):
    '''
    This function lists all file in a the path

    :param path:path of adir
    :return:all files path in a list

    '''
    complete_path = "{}/{}/ABRICATE".format(path,organism)
    files = list_directories(complete_path)
    for file in files:
        if "_argannot.tsv" in file:
            file_path = "{}/{}".format(complete_path, file)
            all_paths.append(file_path)
    return all_paths

def print_status(files):
    '''
    This function prints info about the files selected

    :param files: List of fna files
    :return: None
    '''

    message = "A total of {} file containing _argannot.tsv in their name were found." \
              "Next step, building a directory of these files".format(len(files))
    print(message)
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

def move_file(file_path, new_dir):
    '''
    This function will move the fna files to the created folder

    :param bof:
    :param new_dir:
    :return: none
    '''
    out_name = file_path.strip().split("/")[-1]
    if os.path.exists(out_name):
        print("{} already in the dir".format(out_name))
    else:
        shutil.copy(file_path, new_dir)
        print("{} was moved to the dir".format(file_path))

def list_files_new_source(path):
    '''
    This function lists all file in a the path
    :param path: path of a dir
    :return: all files path in a list
    '''
    files = [file for file in listdir(path) if isfile(join(path,file))]
    return (files)

def parse_genes(arg_db):
    all_genes = []
    for org, genes in arg_db.items():
        all_genes.append(genes)
    flat_list = [item for sublist in all_genes for item in sublist]
    final_list = set(flat_list)
    return(final_list)


def get_all_genes(metadata):

    all_genes = []
    for k, v in metadata.items():
        for gene in v:
            all_genes.append(gene)
    return all_genes

def build_class_df(df_info, classes, metadata):
    '''
    Creates a dictionary of organism and their genes

    df_info = if a dict that is updated
    classes = a list of all possible genes across all genomes
    metadata = dict of k= organism and v= the present genes
    '''
    for organism, genes in metadata.items():
        #organism_classes = sorted(([s[s.find("(") + 1:s.find(")")] for s in genes]))
        abundance = Counter(genes)
        array = get_array(abundance, classes)
        df_info[organism] = array
    return df_info

        
  

def get_array(genes_dict, classes):
    '''
    Build the individual array for each organism 

    If the organism has a gene, it receives the count 
    If not receives 0 
    '''
    array = []
    for x in classes:
        if x in genes_dict.keys():
            array.append(genes_dict[x])
        else:
            value = 0
            array.append(value)
    return (array)

def parse_arg(path, meta):
    '''
     Parse arg info from the file and build a dictionary
    '''

    arg_genes = []
    arg_file = str([file for file in list_files_simple(path) if "_argannot.tsv" in file])
    file_path = "{}/{}".format(path,arg_file[2:-2])
    with open(file_path) as arg:
        for line in arg:
            if not line.startswith("#FILE"):
                info = line.strip().split()[4]
                arg_genes.append(info)
                name = file_path.strip().split("/")[-1][0:-13]
                meta[name] = arg_genes
    return(meta)


def parse_card(path, meta):

    '''
     Parse card info from the file and build a dictionary
    '''


    arg_genes = []
    arg_file = str([file for file in list_files_simple(path) if "_card.tsv" in file])
    file_path = "{}/{}".format(path,arg_file[2:-2])
    with open(file_path) as arg:
        for line in arg:
            if not line.startswith("#FILE"):
                ident = float(line.strip().split()[9])
                if ident >= 60:
                    info = line.strip().split()[4]
                    arg_genes.append(info)
                    name = file_path.strip().split("/")[-1][0:-9]
                    meta[name] = arg_genes
    return(meta)

def parse_vir(path, meta):
    '''
     Parse vir info from the file and build a dictionary
    '''


    arg_genes = []
    arg_file = str([file for file in list_files_simple(path) if "_vfdb.tsv" in file])
    file_path = "{}/{}".format(path,arg_file[2:-2])
    with open(file_path) as arg:
        for line in arg:
            if not line.startswith("#FILE"):
                ident = float(line.strip().split()[9])
                if ident >= 60:
                    info = line.strip().split()[4]
                    arg_genes.append(info)
                    name = file_path.strip().split("/")[-1][0:-9]
                    meta[name] = arg_genes
    return(meta)
   

def main():

    """Main code of the script"""

    #Get thefiles
    
    dirpath=os.getcwd()
    file_path = "{}/agr_vir".format(dirpath)
    directories = list_directories(file_path)
    '''
    ARG analysis

    '''

    #Get all agr info
    metadata_arg = {}

    metadata_card = {}

    metadata_vir = {}

    for organism in directories:
        if ".DS_Store" not in organism:
            organism_path = "{}/{}".format(file_path, organism)
            parse_arg(organism_path, metadata_arg)
            parse_card(organism_path, metadata_card)
            parse_vir(organism_path, metadata_vir)

    #Get all possible genes

    all_genes = parse_genes(metadata_arg)
    all_genes = sorted(all_genes)

    #Build dataframe

    df_info = {}

    df_major_classes = build_class_df(df_info, all_genes, metadata_arg)
    df = pd.DataFrame.from_dict(df_major_classes, orient='index', columns=['(Bla)AmpH', '(Bla)Penicillin_Binding_Protein_Ecoli', '(Flq)OqxBgb', '(Phe)CatB4'])
    #df.to_csv("pantoea_arg.csv")
    #scale
    sns.set(font_scale=0.5)
    
    ##Plots
    sns.set_context("talk", font_scale=0.4)
    bar_plot = df.plot.barh(stacked=True, colormap='vlag')
    plt.legend(bbox_to_anchor=(0.6, 0.95), loc='upper left', borderaxespad=0)
    
    #plt.show()
    sns.set_color_codes()
    cluster_plot = sns.clustermap(data=df, method="single", cmap="vlag")
    cluster_plot.fig.subplots_adjust(right=0.7)
    cluster_plot.ax_cbar.set_position((0.8, .2, .03, .4))
    
    #plt.show()
    

    
    bar_plot.figure.savefig("arg_barplot.pdf", bbox_inches='tight')
    cluster_plot.savefig("arg_cluster.pdf", bbox_inches='tight')




    '''
    Card analysis

    '''
    
    all_card_genes = parse_genes(metadata_card)
    all_card_genes = sorted(all_card_genes)

    
    #Build dataframe

    df__card_info = {}

    df_major_card_classes = build_class_df(df__card_info, all_card_genes, metadata_card)
    df_card = pd.DataFrame.from_dict(df_major_card_classes, orient='index', columns=['CRP', 'Escherichia_coli_emrE', 'H-NS', 'Streptomyces_cinnamoneus_EF-Tu_mutants_conferring_resistance_to_elfamycin', 'YojI', 'acrB', 'aminocoumarin_resistant_alaS', 'aminocoumarin_resistant_cysB', 'arnA', 'bacA', 'baeR', 'cpxA', 'cpxR', 'emrB', 'emrR', 'mdtB', 'mdtC', 'mdtD', 'mexB', 'mexN', 'mfd', 'msbA', 'oqxB', 'sdiA', 'smeB', 'tolC', 'vgaC'])
    df_card.to_csv("pantoea_card.csv")   
    #scale
    sns.set(font_scale=0.65)
    
    ##Plots
    sns.set_context("talk", font_scale=0.4)
    bar_plot_card = df_card.plot.barh(stacked=True, colormap='vlag')
    plt.legend(bbox_to_anchor=(1.01, 0.95), loc='upper left', borderaxespad=0)
    
    sns.set_color_codes()
    cluster_plot_card = sns.clustermap(data=df_card, method="single", cmap="vlag")
    cluster_plot_card.fig.subplots_adjust(right=0.9)
    cluster_plot_card.ax_cbar.set_position((1.01, .2, .03, .4))

    #plt.show()

    ## Save all plots

    bar_plot_card.figure.savefig("card_barplot.pdf", bbox_inches='tight')
    cluster_plot_card.savefig("card_cluster.pdf", bbox_inches='tight')

    

    '''
    Vir analysis

    '''
    
    all_vir_genes = parse_genes(metadata_vir)
    all_vir_genes = sorted(all_vir_genes)

    
    #Build dataframe

    df__vir_info = {}

    df_major_vir_classes = build_class_df(df__vir_info, all_vir_genes, metadata_vir)
    df_vir = pd.DataFrame.from_dict(df_major_vir_classes, orient='index', columns=['astA', 'cheA', 'cheD', 'cheW', 'cheY', 'cheZ', 'clpV1', 'entB', 'flgD', 'flgE', 'flhA', 'flhC', 'flhD', 'fliA', 'fliC', 'fliG', 'fliI', 'fliM', 'fliN', 'fliP', 'icmF1/tssM1', 'luxS', 'ompA', 'pilL', 'pscR', 'tsr', 'tssH-5/clpV', 'vgrG1a', 'vgrG1b'])
    #scale
    #df_vir.to_csv("pantoea_vir.csv")
    sns.set(font_scale=0.65)
    
    ##Plots
    sns.set_context("talk", font_scale=0.4)
    bar_plot_vir = df_vir.plot.barh(stacked=True, colormap='vlag')
    plt.legend(bbox_to_anchor=(1.01, 0.95), loc='upper left', borderaxespad=0)
    
    sns.set_color_codes()
    cluster_plot_vir = sns.clustermap(data=df_vir, method="single", cmap="vlag")
    cluster_plot_vir.fig.subplots_adjust(right=0.68)
    cluster_plot_vir.ax_cbar.set_position((0.8, .2, .03, .4))

    
    ## Save all plots

    bar_plot_vir.figure.savefig("vir_barplo.pdf", bbox_inches='tight')
    cluster_plot_vir.savefig("vir_cluster.pdf", bbox_inches='tight')


        




#main
if __name__=='__main__':
    main()
