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

from Bio.Seq import Seq
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
Please contact Jose at jcedielbecerra@ufl.edu if you have any issues                                                                                       
Version: """+VERSION

np.random.seed(53000)
stime = time.time()

def parse_arguments():
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawTextHelpFormatter)
    input_group = parser.add_argument_group("Input Options")
    hmmer_group = parser.add_argument_group("HMMER Options")
    out_group = parser.add_argument_group("Output Options")
    input_group.add_argument('--genomes_dir', type=str, nargs='+', help='Directory(ies) name containing the Genbanks (*.gbff/*.gbk/*.gb) genomes.', required=True)
    hmmer_group.add_argument('--e_value', type=float, default=1e-4, help='E-value threshold  wanted for hmmsearch (default: 1e-4).', required=False)    
    hmmer_group.add_argument('--threads', type=int, default= int(multiprocessing.cpu_count()), help='CPUs wanted for hmmsearch (default: all available).', required=False)
    out_group.add_argument('--out', type=str, help='Output directory name that will contain the proteins, the dmnd_database, and the modular domtblout table.', required=True)
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

def get_list_genomes(genomes_dir):
    gfiles = []
    for pattern in genomes_dir:
        dir_pattern, file_pattern = os.path.split(pattern)
        if not file_pattern:
            for geno_ext in GENBANK_EXTENSIONS:
                dirs = glob.glob(dir_pattern)
                for dirg in dirs:
                    gfiles.extend(glob.glob(os.path.join(dirg, geno_ext)))
        else:
            dirs = glob.glob(dir_pattern)
            for dirg in dirs:
                gfiles.extend(glob.glob(os.path.join(dirg, file_pattern)))
    return gfiles

def replace_dashes_filenames_genomes_dir(gfiles):
    new_gfiles = []
    for file_path in gfiles:
        directory = os.path.dirname(file_path)
        filename = os.path.basename(file_path)
        if '--' in filename:
            new_filename = filename.replace('--', '_')
            new_filepath = os.path.join(directory, new_filename)
            os.rename(file_path, new_filepath)
            print(f"Renamed file to avoid the '--' string: '{file_path}' to '{new_filepath}'")
        else:
            new_filepath = file_path    
        new_gfiles.append(new_filepath)
    return new_gfiles

def dbfaa_from_gb_dir(genbank_files: List, db_faa: str) -> None:
    ## What function does:
     # Parse out protein sequences from GenBank files and writes them to a FASTA file.
    ## Arguments:
     # genomes_dir (str): paths containing GenBank files (*.gbk or *.gbff).
     # db_faa (str): output file name for the generated FASTA file containing protein sequences that will become dmnd db.
    ## Returns:
     #  None  
    with open(db_faa, 'w') as out_fh:
        for file_path in genbank_files:
            genome = os.path.basename(file_path)
            with open(file_path, 'r') as in_fh:
                for rec in SeqIO.parse(in_fh, 'genbank'):
                    for feat in rec.features:
                        if feat.type == 'CDS':
                            name, seq = None, None
                            if 'locus_tag' in feat.qualifiers:
                                name = "|-|".join([feat.qualifiers['locus_tag'][0]+"|_|" + genome, str(int(feat.location.start)), str(int(feat.location.end)), rec.id])
                            elif 'protein_id' in feat.qualifiers:
                                name = "|-|".join([feat.qualifiers['protein_id'][0]+"|_|" + genome, str(int(feat.location.start)), str(int(feat.location.end)), rec.id])
                            try:
                                if 'translation' in feat.qualifiers:
                                    seq = feat.qualifiers['translation'][0]
                                else:
                                    if feat.location.strand == 1:
                                        seq = str(Seq(rec.seq[feat.location.start:feat.location.end]).translate())
                                    else:
                                        seq = str(Seq(rec.seq[feat.location.start:feat.location.end].reverse_complement()).translate())
                                    print(f"Warning: {feat.location} does not have sequence, but it was fixed")
                            except Exception as e:
                                print(f"Warning: Skipping {feat.location} due to missing translation. Error: {e}")
                                seq = None
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

## getting list of genomes

genomes_list = get_list_genomes(args.genomes_dir)

## replace dashes in filenames
new_genomes_list = replace_dashes_filenames_genomes_dir(genomes_list)
                      
## making protein database 
print("[2]" + print_datetime(), 'Generating the Protein Database')
dbfaa_from_gb_dir(new_genomes_list, f'{args.out}/{args.out}.faa')

## making diamond database 
print("[3]" + print_datetime(), 'Generating the Diamond  Database')
create_diamond_database(f'{args.out}/{args.out}.faa', f'{args.out}/{args.out}.dmnd')

## running hmmsearch
print("[4]" + print_datetime(), 'Making the Modular Domtblout Table for the Protein Database') 
run_hmmsearch(MODULAR_DOMAINS_HMMDB, f'{args.out}/{args.out}.faa', f'{args.out}/{args.out}.domtblout', args.threads, args.e_value)

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