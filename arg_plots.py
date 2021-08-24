#!/user/bin/envpython3
"""
Author:GustavoTamasco
Scripttorun:

Thescripttakes:

Runthecodeby:ex_p5.py<reference_fasta_file><related_fasta_file>
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
    for organism, genes in metadata.items():
        #organism_classes = sorted(([s[s.find("(") + 1:s.find(")")] for s in genes]))
        abundance = Counter(genes)
        array = get_array(abundance, classes)
        df_info[organism] = array
    return df_info

        
  

def get_array(genes_dict, classes):
    array = []
    for x in classes:
        if x in genes_dict.keys():
            array.append(genes_dict[x])
        else:
            value = 0
            array.append(value)
    return (array)

def parse_arg(path, meta):

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
            

def main():

    """Main code of the script"""

    #Get thefiles
    
    dirpath=os.getcwd()
    file_path = "{}/agr_vir".format(dirpath)
    directories = list_directories(file_path)
    

    #Get all agr info
    metadata_arg = {}
    for organism in directories:
        if ".DS_Store" not in organism:
            organism_path = "{}/{}".format(file_path, organism)
            parse_arg(organism_path, metadata_arg)

    #Get all possible genes

    all_genes = parse_genes(metadata_arg)
    all_genes = sorted(all_genes)

    #Build dataframe

    df_info = {}

    df_major_classes = build_class_df(df_info, all_genes, metadata_arg)
    df = pd.DataFrame.from_dict(df_major_classes, orient='index', columns=['(Bla)AmpH', '(Bla)Penicillin_Binding_Protein_Ecoli', '(Flq)OqxBgb', '(Phe)CatB4'])
    #scale
    sns.set(font_scale=0.65)
    
    ##Plots

    bar_plot = df.plot.barh(stacked=True, colormap='vlag')
    
    #plt.show()
    sns.set_color_codes()
    sns.clustermap(data=df, method="single", cmap="vlag")
    plt.show()

    #cmap="RdYlBu_r"


    """
    fig = plt.figure(figsize=(12,8))
    ax = plt.axes(projection='3d')

    sctt = ax.scatter3D(X_pca_3[:, 0], X_pca_3[:, 1], X_pca_3[:, 2], c = df, s = 50, alpha=0.6)

    plt.title("3D scatterplot", pad=15)
    ax.set_xlabel("First PCA")
    ax.set_ylabel("Second PCA")
    ax.set_zlabel("Third PCA")

    plt.savefig("3D_scatterplot.pdf")
   

    #full_plot.savefig("plasmid_arg.pdf", bbox_inches='tight')
    #not_full.savefig("plasmid_arg_scale1.pdf", bbox_inches='tight')

    """

        


    """
    '''Plasmids'''
    plasmid_path = "{}{}".format(path_to_all_info,directories[1])
    os.chdir(plasmid_path)
    plasmid_dir = list_directories(plasmid_path)
    for organism in plasmid_dir:
        arg_files = list_files(all_file_path,plasmid_path,organism)
    print_status(arg_files)


    '''Building a dir of fna files'''
    plasmid_arg_path = "{}/plasmid_argannot_files".format(dirpath)
    create_genomes_dir(plasmid_arg_path)
    os.chdir(plasmid_arg_path)
    for file in arg_files:
        move_file(file, plasmid_arg_path)

    '''Building metadata'''
    #All genes to each genome
    plasmid_agr_metadata = {}
    plasmid_arg = list_files_simple(plasmid_arg_path)
    for f in plasmid_arg:
        with open(f) as plasmid_arg_info:
            parse_genes(f, plasmid_arg_info, plasmid_agr_metadata)
    #All genes that occured
    all_genes = sorted(set(get_all_genes(plasmid_agr_metadata)))

    #All arg classes
    all_arg_classes = list(filter(None,sorted(set([s[s.find("(")+1:s.find(")")] for s in all_genes]))))
    print(all_arg_classes)

    '''Build dataframe for the classes plot'''
    df_info = {}

    df_major_classes = build_class_df(df_info, all_arg_classes, plasmid_agr_metadata)
    df = pd.DataFrame.from_dict(df_major_classes, orient='index', columns=['AGly', 'AGly/Flqn', 'Bla', 'Flq', 'MLS', 'Phe', 'PheCmlB', 'Rif', 'Sul', 'Tet', 'Tmt'])
    #df = df.transpose()
    #df.to_csv('arg_genes.csv', sep='\t', encoding='utf-8')
    sns.set(font_scale=0.65)
    # Need both
    not_full = sns.clustermap(df, label='small', cmap="vlag", standard_scale=1, linewidths=0)
    full_plot = sns.clustermap(df, label='small', cmap="vlag", linewidths=0)
    # plt.title('Antibiotic resistance genes across 34 organism', fontsize=15)
    # sns.set(font_scale=1)
    plt.show()
    full_plot.savefig("plasmid_arg.pdf", bbox_inches='tight')
    not_full.savefig("plasmid_arg_scale1.pdf", bbox_inches='tight')


    """



#main
if __name__=='__main__':
    main()
