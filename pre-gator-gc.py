#!/bin/env python

import time
import sys
import os
import glob
import subprocess
import tempfile
import argparse
import numpy as np
import multiprocessing

from Bio import SeqIO
from datetime import datetime
from typing import List, Set, Dict

## Define constants
GENBANK_EXTENSIONS = ['*.gbk', '*.gbff', '*.gb']
MODULAR_DOMAINS_HMMDB = 'flat/modular_domains.hmmdb'
VERSION = 'v1.0.0'
DESCRIPTION = """
                                                                                                                                                                     
     -\ ---\--\ -------\ ----\ ---\--\ ---\ --\ ----\ ----\--------\ /---                         
   /--/ ---/--/ -------/ ----/ ---/--/ ---/ --/ ----/ ----/--------/ \----\                                                                    
  ________   _____  __________________  __________         _________________                                                                                
 /  _____/  /  _  \ \__    ___/_____  \ \______   \       /  _____/\_   ___ \                                                                     
/   \  ___ /  /_\  \  |    |   /   |   \ |       _/ _____    \  ___/    \  \/                     
\    \_\  \    |    \ |    |  /    |    \|    |   \/_____/    \_\  \     \_____                   
 \______  /____|__  / |____|  \_______  /|____|_  /       \______  /\_______  /                  
        \/        \/                  \/        \/               \/         \/                                                            
    -----\--\ -------\ ----\ ---\--\ ---\ --\ ------\ ----\---------\ /----                                                               
    \----/--/ -------/ ----/ ---/--/ ---/ --/ ------/ ----/---------/ \----/                                                                       

GATOR-GC: Genomic Assessment Tool for Orthologous Regions and Gene Clusters                                                                               
Developer: José D. D. Cediel-Becerra
Code reviewers: Valérie de Crécy-Lagard and Marc G. Chevrette                                                                               
Afiliation: Microbiology & Cell Science Deparment, University of Florida                                                                               
Please contact me at jcedielbecerra@ufl.edu if you have any questions                                                                                       
Version: """+VERSION

np.random.seed(53000)
stime = time.time()

def parse_arguments():
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--genomes_dir', type=str, help='Directory name containing the Genbanks (*.gbff/*.gbk/*.gb) genomes.', required=True)
    parser.add_argument('--proteins', type=str, help='Name for the precompiled protein database (.faa) of genomes_dir.', required=True)
    parser.add_argument('--dmnd_database', type=str, help='Name for the precompiled diamond database (.dnmd).', required=True)
    parser.add_argument('--e_value', type=float, default=1e-4, help='E-value threshold  wanted for hmmsearch (default: 1e-4).', required=False)    
    parser.add_argument('--modular_domtblout', type=str, help='Name for the precomputed hmmsearch domain table for modular domains', required=True)
    parser.add_argument('--threads', type=int, default= int(multiprocessing.cpu_count()), help='CPUs wanted for hmmsearch (default: all available).', required=False)
    parser.add_argument('--out', type=str, help='Output directory name that will contain the proteins, the dmnd_database, and the modular domtblout table.', required=True)
    return parser.parse_args()

def print_datetime():
    return "[" + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "]"

def create_directory(directory_name: str) -> None:
    ## What function does:                                                                                                                                
     # Create the directory to save the PRE-GATOR-GC output                                                                                                     
    ## Arguments:
     # directory_name (str): Name for the directory
    ## Returns:
     # None: This function Just create the directory
    if '--' in directory_name:
        sys.stderr.write("ERROR: Directory name "+directory_name+" contains an invalid character for GATOR-GC. Do not use '--'\n")
        sys.exit(1)
    try:
        os.mkdir(directory_name)
        print("[1]" + print_datetime(), f"The {directory_name} Directory was Created")
    except FileExistsError:
        sys.stderr.write("ERROR: Directory "+directory_name+" Already Exists"+"\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write("ERROR: Failed to Make the "+directory_name+" Directory"+"\n")
        sys.stderr.write(str(e))
        sys.exit(1)
        
def replace_dashes_filenames_genomes_dir(genomes_dir):
    genbank_files = [geno_file for geno_ext in GENBANK_EXTENSIONS for geno_file in glob.glob(os.path.join(genomes_dir, geno_ext))]
    for filepath in genbank_files:
        filename = os.path.basename(filepath)
        if '--' in filename:
            new_filename = filename.replace('--', '_')
            new_filepath = os.path.join(genomes_dir, new_filename)
            os.rename(filepath, new_filepath)
            print(f"Renamed file to avoid the '--' string: '{filepath}' to '{new_filepath}'")
        
def dbfaa_from_gb_dir(genomes_dir: str, db_faa: str) -> None:
    ## What function does:
     # Parse out protein sequences from GenBank files and writes them to a FASTA file.
    ## Arguments:
     # genomes_dir (str): paths containing GenBank files (*.gbk or *.gbff).
     # db_faa (str): output file name for the generated FASTA file containing protein sequences that will become dmnd db.
    ## Returns:
     #  None  
    genbank_files = [geno_file for geno_ext in GENBANK_EXTENSIONS for geno_file in glob.glob(os.path.join(genomes_dir, geno_ext))]
    with open(db_faa, 'w') as out_fh:
        for file_path in genbank_files:
            if '--' in file_path:
                sys.stderr.write(f"Error: Genbank File {file_path} contains an invalid character for GATOR-GC. Do not use '--'\n")
                sys.exit(1)
            genome = os.path.basename(file_path)
            with open(file_path, 'r') as in_fh:
                for rec in SeqIO.parse(in_fh, 'genbank'):
                    for feat in rec.features:
                        if feat.type == 'CDS':
                            name, seq = None, None
                            if 'locus_tag' in feat.qualifiers:
                                name = "|-|".join([genome, feat.qualifiers['locus_tag'][0], str(int(feat.location.start)), str(int(feat.location.end)), rec.id])
                            if 'translation' in feat.qualifiers:
                                seq = feat.qualifiers['translation'][0]
                            if name is not None and seq is not None:
                                out_fh.write('>'+name+"\n"+seq+"\n")

def create_diamond_database(db_faa: str, database_name: str) -> None:
    ## What function does
     # Create a dmnd db using the provided FASTA file.
    ## Arguments:
     # db_faa (str): Path to the input FASTA file containing protein sequences.
     # database_name (str): Name of the dmnd db to be created.
    ## Returns:
     # None: This function does not return anything. It simply creates the dmnd db file.
    command = ' '.join(['diamond',
                        'makedb',
                        '--in', db_faa,
                        '--quiet',
                        '-d', database_name,
                        '--threads', str(args.threads)])
    subprocess.run([command], shell=True)
    
def run_hmmsearch(hmmdb: str, query: str, out: str, cpu: int, e_value: float) -> None:
    ## What function does:
     # Run the hmmsearch pipeline to search for modular domain hits
    ## Arguments:
     # hmmdb (str): A string for the HMMs modular domains (nrps, pks) database name
     # query (str): A string for the query protein file name
     # out (str): A string for the modular domtblout table name
     # cpu (int): A integer for the number of CPUs wanted to run hmmsearch
     # e_value (float): A float for the E-value threshold wanted to run hmmsearch
    ## Returns:
     #  None: This function does not return anything. It simply creates the modular domtblout table
    command = ' '.join(['hmmsearch',
                        '--domtblout', out,
                        '--cpu', str(cpu),
                        '--domE', str(e_value),
                        '--noali',
                        '-o /dev/null',
                        hmmdb,
                        query])
    subprocess.run([command], shell=True)

#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
############################# BEGIN MAIN BLOCK ##################################

## print description
print(DESCRIPTION)

## parse arguments
args = parse_arguments()

## make output directory 
create_directory(args.out)

## replace dashes in filenames
replace_dashes_filenames_genomes_dir(args.genomes_dir)
                      
## making protein database 
print("[2]" + print_datetime(), 'Generating the Protein Database')
dbfaa_from_gb_dir(args.genomes_dir, f'{args.out}/{args.proteins}')

## making diamond database 
print("[3]" + print_datetime(), 'Generating the Diamond  Database')
create_diamond_database(f'{args.out}/{args.proteins}', f'{args.out}/{args.dmnd_database}') 

## running hmmsearch
print("[4]" + print_datetime(), 'Making the Modular Domtblout Table for the Protein Database') 
run_hmmsearch(MODULAR_DOMAINS_HMMDB, f'{args.out}/{args.proteins}', f'{args.out}/{args.modular_domtblout}', args.threads, args.e_value)

## elapsed time 
etime = time.time()
ftime = round((etime - stime) / 60, 2)
if ftime < 1:
    ftime *= 60
    time_unit = "seconds"
elif ftime < 60:
    time_unit = "minutes"
else:
    ftime /= 60
    time_unit = "hours"
print("[5]" + print_datetime(), f"Elapsed time: {ftime:.1f} {time_unit}")
