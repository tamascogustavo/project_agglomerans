#!/user/bin/envpython3
"""
Author:GustavoTamasco


This script is based on the information found at https://www.metagenomics.wiki/tools/fastq/ncbi-ftp-genome-download
Scripttorun:

Thescripttakes:

"""

#importstatements
import os
import pandas as pd
import subprocess

#functions and classes

def parse_refseqs(info):
    '''
    This function will get the refseqs from a csv file

    retun: a list of all refseqs

    '''
    all_info = []
    for line in info:
        newline = line.strip().split(";")[1]
        all_info.append(newline)
    return all_info

def generate_metadata_assembly(out_file):
    '''
    This function will the inital metadata

    Return a file with all the the refsummary
    
    '''
    if os.path.exists(out_file):
        print("{} was already generated".format(out_file))
    else:
        cmd = "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/{}}".format(out_file)
        exit_message = subprocess.check_call(cmd, shell=True)
        print("Exit status: {0}".format(exit_message))
        print("{0} Refseq summary was generated".format(out_name))

def get_download_links(lib, target, links):
    '''

    This function will produce the links for download
    '''
    cmd = "grep -E '{}' {} | cut -f 20".format(target, lib)
    link = subprocess.check_output(cmd, shell=True, text=True)
    final_link = "{}/{}_genomic.fna.gz".format(link.strip("\n"), link.strip('\n').split("/")[-1])
    #print(final_link)

    links.append(final_link)
    
def download(links):
    '''
    This function will download all the genome files

    '''
    for link in links:
        out_name = link.split("/")[-1]
        if os.path.exists(out_name):
            print("{} already downloaded".format(out_name))
        else:
            cmd = "wget {}".format(link)
            subprocess.check_call(cmd, shell= True)
            print("{} was generated".format(out_name))
    

def main():
    """Main code of the script"""

    #Get the ids
    dirpath = os.getcwd()
    file = 'metadata.csv'
    assembly_lib = "assembly_summary_refseq.txt"
    with open(file) as metadata:
        refs = parse_refseqs(metadata)
        refs = list(filter(None, refs))
        refs = list(filter(lambda a: a != 'RefSeq', refs))

    # Get a list of assemblys info 

    generate_metadata_assembly(assembly_lib)

    #Get FTP download link
    download_links = []
    for ref in refs:
        get_download_links(assembly_lib, ref, download_links)

    #Generate the full link 
    formated_download_link =  list(map(lambda s: s.strip(), download_links))
    download(formated_download_link)

        

#main
if __name__=='__main__':
    main()
