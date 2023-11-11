#!/bin/env python

import time
import sys
import csv
import os
import glob
import pandas as pd
import numpy as np
import subprocess
import tempfile
import argparse
import seaborn as sns
import multiprocessing
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from Bio import SeqIO
from scipy.stats import norm
from operator import itemgetter
from datetime import datetime
from pygenomeviz import GenomeViz
from typing import List, Set, Dict
from matplotlib.patches import Patch 

## Constants
NRPS = set(['NRPS-C', 'NRPS-A'])
PKS = set(['PKS-KS', 'PKS-AT'])
MODULE_TYPE = {'Condensation_LCL': 'NRPS-C',
               'Condensation_DCL': 'NRPS-C',
               'Condensation_Dual': 'NRPS-C',
               'Condensation_Starter': 'NRPS-C',
               'Cglyc': 'NRPS-C',
               'AMP-binding': 'NRPS-A',
               'A-OX': 'NRPS-A',
               'Enediyne-KS': 'PKS-KS',
               'Modular-KS': 'PKS-KS',
               'Trans-AT-KS': 'PKS-KS',
               'Hybrid-KS': 'PKS-KS',
               'PKS_KS': 'PKS-KS',
               'Iterative-KS': 'PKS-KS',
               'PKS_AT': 'PKS-AT'}
MODULAR_DOMAINS_HMMDB = 'flat/modular_domains.hmmdb'
VERSION = 'v0.8.5'
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

GATOR-GC: Genomic Assesment Tool for Orthologous Regions and Gene Clusters
Authors: José D.D Cediel-Becerra, Valérie de Crécy-Lagard, Marc G. Chevrette
Afiliation: Microbiology & Cell Science Deparment, University of Florida
Please contact me at jcedielbecerra@ufl.edu if you have any questions
Version:"""+VERSION

np.random.seed(53000)
plt.rcParams['font.family'] = 'Arial'

stime = time.time()

def parse_arguments():
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--required', type=str, help='Query protein fasta file containing required proteins.', required=True)
    parser.add_argument('--optional', type=str, help='Query protein fasta file containing optional proteins.', required=False)
    parser.add_argument('--genomes_dir', type=str, help='Directory name containing the Genbanks (*.gbff/*.gbk/*.gb) genomes.', required=True)
    parser.add_argument('--proteins', type=str, help='Precompiled protein database (.faa) of genomes_dir to search against.', required=True)
    parser.add_argument('--dmnd_database', type=str, help='Precompiled diamond database (.dnmd) to search against.', required=True)
    parser.add_argument('--modular_domtblout', type=str, help='Precomputed hmmsearch domain table for modular domains against the precompiled protein database.', required=True)
    parser.add_argument('--threads', type=int, default= int(multiprocessing.cpu_count()), help='CPUs wanted for diamond search and hmmsearch (default: all CPUs available).', required=False)
    parser.add_argument('--query_cover', type=int, default=70,help='Protein percent query cover for diamond search (default: 70).', required=False)
    parser.add_argument('--identity', type=int, default=35, help='Protein percent identity for diamond search (default: 35).', required=False)
    parser.add_argument('--e_value', type=float, default=1e-4, help='E-value threshold  wanted for hmmsearch (default: 1e-4).', required=False)
    parser.add_argument('--intergenic_distance', type=int, default=86, help='Distance in kilobases between required genes (default: 86 kb)', required=False)
    parser.add_argument('--window_extension', type=int, default=10, help='Extension in kilobases from the start and end positions of the windows (default: 10 kb)', required=False)
    parser.add_argument('--out', type=str, help='Output directory name that will have GATOR-GC results', required=True)
    return parser.parse_args()

def print_datetime():
    return "[" + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "]"

def create_directory(directory_name: str) -> None:
    ## What function does:
     # Create the directory to save the GATOR-GC output
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

def init_modular(protein_dict: Dict) -> Dict:
    ## What function does:
     # Create a dictionary structure for query proteins and make subkeys for boolean decisions (is_nrps, is_pks)
    ## Arguments:
     # protein_dict (Dict): Dictionary containing as keys the headers for required and optional proteins and subkeys for boolean decision (is_nrps, is_pks)  
    ## Returns:
     # protein_dict (Dict): Same protein_dict dictinoary structure
    for header in protein_dict:
        protein_dict[header]['is_nrps'] = False
        protein_dict[header]['is_pks'] = False
    return protein_dict

def parse_query_faa(query_faa: str, query_type: str, protein_dict: Dict) -> Dict:
    ## What function does:
     # It adds two more subkeys in the protein_dict dictionary (i.e., query_type and seq)
    ## Arguments:
     # query_faa (str): A string for the query protein file name
     # query_type (str): A string for the query protein type (i.e., required and optional)
     # protein_dict (str): A Dictionary containing as keys the headers for required and optional proteins and subkeys for boolean decision (is_nrps, is_pks)
    ## Returns:
     # init_modular(protein_dict) (def): Call the init_modular function to create the protein_dict
    with open(query_faa, 'r') as ifh:
        header, sequence = None, ''
        for line in ifh.read().splitlines():
            if line.startswith('>'):
                if header is not None:
                    if header not in protein_dict:
                        protein_dict[header] = {}
                        protein_dict[header]['query_type'] = query_type
                    protein_dict[header]['seq'] = sequence
                header, sequence = line[1:], ''
            else:
                sequence += line
        if header not in protein_dict:
            protein_dict[header] = {}
            protein_dict[header]['query_type'] = query_type
        protein_dict[header]['seq'] = sequence
    return init_modular(protein_dict)

def run_hmmsearch(hmmdb: str, query: str, out: str, cpu: int, e_value: float) -> None:
    ## What function does:
     # Run the hmmsearch pipeline to search for modular domain hits
    ## Arguments:
     # hmmdb (str): A string for the HMMs modular domains (nrps, pks) database name
     # query (str): A string for the query protein file name
     # out (str): A string for the modular domtblout table name
     # cpu (int): A integer for the number of CPUs wanted to run hmmsearch
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

def parse_modular_domtblout(domtbl: str, protein_dict: Dict) -> Dict:
    ## What function does:
     # Parse the modular domtblout table and categorize the domain hits according to the MODULE_TYPE constant
    ## Arguments:
     # domtbl (str): A string for the modular domtblout table name  
     # protein_dict (Dict): A dictionary containing as keys the headers for required and optional proteins, subkeys for boolean decision (is_nrps, is_pks) and query_type and seq. 
    ## Returns:
     # modular_domain_hit (Dict): A dictionary containing the MODULE_TYPE for each hit found
    modular_domain_hit = {}
    with open(domtbl, 'r') as dfh:
        for line in dfh.read().splitlines():
            if not line.startswith('#'):
                l = line.split()
                query_protein, domain_hit = l[0], l[3]
                if query_protein not in modular_domain_hit:
                    modular_domain_hit[query_protein] = set()
                modular_domain_hit[query_protein].add(MODULE_TYPE[domain_hit])
    return modular_domain_hit

def parse_modular_domain_query_hits(modular_domain_hit: Dict, protein_dict: Dict) -> Dict:
    ## What function does:
     # Make a boolean decision (True or False) for the modular_domain_hit depending on the MODULE_TYPE categorization.
    ## Arguments:
     # modular_domain_hit (Dict): A dictionary containing the MODULE_TYPE categorization for each hit found in the modular domtblout table 
     # protein_dict (Dict): Same protein_dict dictionary structure than before with is_nrps and is_pks with False values as default
    ## Returns:
     # protein_dict (Dict): Same protein_dict dictionary structure than before but making the boolean decision for is_nrps and is_pks       
    for query_protein in modular_domain_hit:
        if NRPS.issubset(modular_domain_hit[query_protein]):
            protein_dict[query_protein]['is_nrps'] = True
        if PKS.issubset(modular_domain_hit[query_protein]):
            protein_dict[query_protein]['is_pks'] = True
    return protein_dict

def parse_modular_domain_genome_hits(modular_domain_hit: Dict, protein_dict: Dict) -> Dict:
    ## What function does:
     # Adds the modular domain hits found in the precompiled protein database in the protein_dict dictionary under the subkey 'hit'
    ## Arguments:
     # modular_domain_hit (Dict): A dictionary containing the MODULE_TYPE categorization for each hit found in the modular domtblout table    
     # protein_dict (Dict): Same protein_dict dictionary structure than before but with the is_nrps and is_pks boolean decision performed
    ## Returns:
     # protein_dict (Dict): Same protein_dict dictionary structure but with a new subkey 'hit', which will be feeded with the domain hit 
    modular_hit_set = set() 
    for genome_protein in modular_domain_hit:
        if NRPS.issubset(modular_domain_hit[genome_protein]) or PKS.issubset(modular_domain_hit[genome_protein]):
            modular_hit_set.add(genome_protein)
            modular_hit_dict = {'hit': modular_hit_set}
    for query_protein in protein_dict: 
        if protein_dict[query_protein]['is_nrps'] or protein_dict[query_protein]['is_pks']:
            protein_dict[query_protein]['hit'] = modular_hit_dict['hit']
    return protein_dict

def run_diamond(query: str, dmnddb: str, dmnd_out: str, threads: int, query_cover: float, identity: float, max_matches: int = int(1e50)) -> None:
    ## What function does
     # Run the diamond blastp search
    ## Arguments:
     # query (str): A string for the query file name containing the precompiled protein database
     # dmnddb (str): A string for the precomputed diamond database
     # dmnd_out (str): A string with the desire name for the diamond table output
     # threads (int): Integer for the number of CPUs wanted to run diamond (default: maximum available in the machine).
     # query_cover (float): A float for the protein percent query cover wanted for diamond search (default: 70)
     # identity (float): A float for the protein percent identity wanted for diamond search (default: 35)
     # max_matches (int): A integer for the maximum number of matches to report per query (default: int(1e50))
    ## Returns:
     # None: This function does not return anything explicitly, as it writes the output to the specified output file.   
    command = ' '.join(['diamond',
                        'blastp',
                        '-q', query,
                        '-d', dmnddb,
                        '-o', dmnd_out,
                        '--quiet',
                        '--threads', str(threads),
                        '--query-cover', str(query_cover),
                        '--id', str(identity),
                        '-k', str(max_matches)
                        ])
    subprocess.run([command], shell=True)

def parse_diamond_search(dmnd: str, protein_dict: Dict) -> None:
    ## What function does:
     # Adds the diamond hits in the protein_dict dictionary under the subkey 'hit'
    ## Arguments:
     # dmnd (str): A string for the diamond table output name
     # protein_dict (Dict): Same protein_dict dictionary structure than before 
    ## Returns:
     # protein_dict (Dict): Same protein_dict dictionary structure but with a new subkey 'hit', which will be feeded with the diamond hits
    with open(dmnd, 'r') as mfh:
        for line in mfh:
            columns = line.strip().split('\t')
            query_protein, dmnd_hit = columns[0], columns[1]
            if query_protein not in protein_dict:
                protein_dict[query_protein] = {'hit': set()}
            if 'hit' not in protein_dict[query_protein]:
                protein_dict[query_protein]['hit'] = set()
            protein_dict[query_protein]['hit'].add(dmnd_hit)
    return protein_dict

def grouping_user_req_opt_proteins_by_contig(protein_dict: Dict, query_type) -> [Dict, Set]:
    ## What function does:
     # Groups required and optional user proteins by the same contig id and make a new key 'window' to assign na values
    ## Arguments:
    # protein_dict (Dict): Same protein_dict dictionary structure than before
    # query_type (str): A string for the query_type (req, opt)
    ## Returns:
     # hits_by_contig (Dict): A dictionary containing all the required and/or optional hits grouped by the same contig id
     # all_required_proteins (Set): A set containing all the required user proteins 
    hits_by_contig = {}
    all_required_proteins = set()
    for query_proteins in protein_dict:
        if protein_dict[query_proteins]['query_type'] == query_type:
            all_required_proteins.add(query_proteins)
            for hit in protein_dict[query_proteins]['hit']:
                locus, start_position, end_position, genome_contig_id = (hit.split('|-|')[1:4] + [hit.split('|-|')[0] + '|-|' + hit.split('|-|')[-1]])
                if genome_contig_id not in hits_by_contig:
                    hits_by_contig[genome_contig_id] = []
                hits_by_contig[genome_contig_id].append({'query_protein': query_proteins,
                                                         'query_type': protein_dict[query_proteins]['query_type'],
                                                         'start_position': int(start_position),
                                                         'end_position': int(end_position),
                                                         'locus': locus,
                                                         'genome_contig_id': genome_contig_id,
                                                         'window': 'na'})  
    return hits_by_contig, all_required_proteins

def checking_all_required_proteins(init_window: List, all_required_proteins: Set, final_window: List[Dict]) -> None:
    ## What function does:
     # Checks if all elements in `all_required_proteins` are present in the init_window. If yes, it appends the init_window to the final_window.
    ## Arguments:
     # intit_window (List): A temporary list used to hold dictionaries that belong to the same contig.
     # A new window is started when the ‘genome_contig_id’ value changes or when intergenic distance is larger than threshold. The init_window list is cleared each time a new group starts
     # all_required_proteins (Set): A set of required user proteins that are required to be present in each init_window.
     # final_window (List): A list of dictionaries containing the hits for each window
    if all_required_proteins == {data['query_protein'] for data in init_window}:
        final_window.append(init_window)
        
def checking_distance_bw_loci(req_hits_by_contig: Dict[str, List[Dict[str, int]]], all_required_proteins: Set) -> List[Dict]:
    ## What function does:
     # Define the windows based on the same contig, intergenic distance between loci and the presence of all required user proteins
    ## Arguments:
     # req_hits_by_contig (Dict[str, List[Dict[str, int]]]): A dictionary containing lists of dictionaries with required hits grouped by the same contig id
     # all_required_proteins (Set): A set containing all the required user proteins
    ## Returns:
     # final_window (List[Dict]): A list of dictionaries containing the windows that passed all conditions of contig id, intergenic distance and presence of all required user proteins
    init_window = []
    final_window = []
    has_all_req = []
    for genome_contig_id in req_hits_by_contig:
        if all_required_proteins == {req['query_protein'] for req in req_hits_by_contig[genome_contig_id]}:
            req_hits_by_contig[genome_contig_id].sort(key=itemgetter('start_position'))
            for required_dict in req_hits_by_contig[genome_contig_id]:
                has_all_req.append(required_dict)
    for required_hits in has_all_req:
        if init_window:
            is_different_contig = required_hits['genome_contig_id'] != init_window[0]['genome_contig_id']
            exceeds_intergenic_distance_cutoff = int(required_hits['start_position']) - int(init_window[-1]['end_position']) > (args.intergenic_distance * 1000)
            if is_different_contig or exceeds_intergenic_distance_cutoff:
                checking_all_required_proteins(init_window, all_required_proteins, final_window)
                init_window = []
        init_window.append(required_hits)
    checking_all_required_proteins(init_window, all_required_proteins, final_window)
    print("[5.1]" + print_datetime(), 'Number of Windows Found:', len(final_window))
    if len(final_window) == 0:
        sys.exit(1)
    return final_window

def generating_windows_genbanks(final_window: List[Dict], output_directory: str) -> None:
    ## What function does:
     # Generates the Genbank files for the windows found
    ## Arguments:
     # final_window (List[Dict]): A list of dictionaries containing the hits that belongs for each window found
     # output_directory (str): A sring containing the directory folder name that will used to save the Genbank files
    ## Returns:
     # None
    window_count = 1
    for windows in final_window:
        start_positions, end_positions = [], []
        for window in windows:
            parsing_gci, parsing_start, parsing_end  = window['genome_contig_id'].split('|-|'), window['start_position'], window['end_position']
            genome, contig = parsing_gci[0], parsing_gci[1]
            start_positions.append(parsing_start)
            end_positions.append(parsing_end)
        min_start, max_end = (min(start_positions) - (args.window_extension * 1000)), (max(end_positions) + (args.window_extension * 1000))
        if min_start < 0:
            min_start = 0
        with open(args.genomes_dir + genome, 'r') as gf: 
            for rec in SeqIO.parse(gf, 'genbank'):
                if rec.id == contig:
                    record_to_write = rec
                    if max_end > len(record_to_write):
                        max_end = len(record_to_write)
                    sub_seq = record_to_write[min_start:max_end]
                    sub_seq.annotations["molecule_type"] = "DNA"
        os.makedirs(output_directory, exist_ok=True)
        basename, extension = os.path.splitext(os.path.basename(genome))
        file_name = f"window_{window_count}--{basename}{extension}"
        output_file = os.path.join(output_directory, file_name)
        with open(output_file, 'w') as outf:
            SeqIO.write(sub_seq, outf, 'genbank')
        window_count += 1

def making_windows_and_optional_table_hits(genbank_dir: str, req_hits_by_contig: Dict[str, List[Dict[str, int]]], opt_hits_by_contig: Dict[str, List[Dict[str, int]]], windows_tsv: str) -> List:
    ## What function does:
     # Export a TSV file containig all required and optional hits organized by windows
    ## Arguments:
     # genbank_dir (str): A string for the directory folder name containing the Genabank window files
     # req_hits_by_contig (Dict[str, List[Dict[str, int]]]): A dictionary containing a list of dictionaries for all required hits
     # opt_hits_by_contig (Dict[str, List[Dict[str, int]]]): A dictionary containing a list of dictionaries for all optional hits 
    ## Returns:
     # None
    genbank_extensions = ['*.gbk', '*.gbff', '*.gb']
    genbank_files = [geno_file for geno_ext in genbank_extensions for geno_file in glob.glob(os.path.join(genbank_dir, geno_ext))]
    final_window_extension = []
    for file_path in genbank_files:
        window_genome = os.path.basename(file_path).split('--')[0]
        num_window = window_genome.split('_')[1]
        for rec in SeqIO.parse(file_path, 'genbank'):
            for feat in rec.features:
                if feat.type == 'CDS':
                    final_window_extension.append({'locus': feat.qualifiers['locus_tag'][0],
                                                   'window': num_window})
    flattened_req_hits_by_contig = [item for sublist in req_hits_by_contig.values() for item in sublist]
    flattened_opt_hits_by_contig = [item for sublist in opt_hits_by_contig.values() for item in sublist]
    req_opt_hits_list = flattened_req_hits_by_contig + flattened_opt_hits_by_contig
    for req_opt_hits in req_opt_hits_list:
        for window_req_opt_hits in final_window_extension:
            if req_opt_hits['locus'] == window_req_opt_hits['locus']:
                req_opt_hits['window'] = window_req_opt_hits['window']
    for hits in req_opt_hits_list:
        try:
            hits['window'] = int(hits['window'])
        except ValueError:
            hits['window'] = float('inf')
    window_opt_hits_by_gci = {}
    for hits in req_opt_hits_list:
        gci = hits['genome_contig_id']
        if gci not in window_opt_hits_by_gci:
            window_opt_hits_by_gci[gci] = []
        window_opt_hits_by_gci[gci].append(hits)
        window_opt_hits_by_gci[gci].sort(key=itemgetter('start_position'))
    sorted_by_window = sorted(window_opt_hits_by_gci.values(), key=lambda dictionaries: min(hits['window'] for hits in dictionaries))
    for original_na in sorted_by_window:
        for hits in original_na:
            hits['window'] = 'window_' + str(hits['window']) if hits['window'] != float('inf') else 'na'
    with open(windows_tsv, 'w', newline='') as f_output:
        window_opt_hits_tsv = csv.DictWriter(f_output, fieldnames=req_opt_hits_list[0].keys(), delimiter='\t')
        window_opt_hits_tsv.writeheader()
        for dictionaries in sorted_by_window:
            window_opt_hits_tsv.writerows(dictionaries)
    return flattened_req_hits_by_contig, flattened_opt_hits_by_contig

def dbfaa_from_gb_dir(genbank_dir: str, db_faa: str) -> None:
    ## What function does:
     # Parse out protein sequences from GenBank files and writes them to a FASTA file.
    ## Arguments:
     # genomes_dir (str): paths containing GenBank files (*.gbk or *.gbff).
     # db_faa (str): output file name for the generated FASTA file containing protein sequences that will become dmnd db.
    ## Returns:
     #  None
    genbank_extensions = ['*.gbk', '*.gbff', '*.gb']
    genbank_files = [geno_file for geno_ext in genbank_extensions for geno_file in glob.glob(os.path.join(genbank_dir, geno_ext))]
    with open(db_faa, 'w') as out_fh:
        for file_path in genbank_files:
            genome = os.path.basename(file_path)
            with open(file_path, 'r') as in_fh:
                for rec in SeqIO.parse(in_fh, 'genbank'):
                    for feat in rec.features:
                        if feat.type == 'CDS':
                            name, seq = None, None
                            if 'locus_tag' in feat.qualifiers:
                                name = "|-|".join([genome, feat.qualifiers['locus_tag'][0]])
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
    
def pa_tables_from_dmnd_windows_faa(input_file: str, output_folder: str) -> None:
    ## What function does
     # Make the gene-level presence-absence (pa) dataframes for each window 
    ## Arguments:
     # input_file (str): A string for the dmnd output for all vs all windows
     # output_folder (str): A string for the directory folder to save the presence-absence dataframes 
    ## Returns:
     # None: It just generates the pa dataframes
    dmnd_out = pd.read_csv(input_file, delimiter="\t", header=None)
    pa_dfs = {}
    for query in dmnd_out[0].unique():
        subset = dmnd_out[dmnd_out[0] == query]
        genomes = subset[1].apply(lambda x: x.split("|-|")[0]).unique()
        pa_df = pd.DataFrame(0, index=genomes, columns=subset[0].apply(lambda x: x.split("|-|")[1]).unique())
        for i in range(len(subset)):
            genome = subset.iloc[i, 1].split("|-|")[0]
            locus = subset.iloc[i, 0].split("|-|")[1]
            pa_df.loc[genome, locus] = 1
        if query.split("|-|")[0] in pa_dfs:
            pa_dfs[query.split("|-|")[0]].append(pa_df)
        else:
            pa_dfs[query.split("|-|")[0]] = [pa_df]
    for query, pa_df_list in pa_dfs.items():
        query_no_ext = os.path.splitext(query)[0]
        concatenated_pa_df = pd.concat(pa_df_list, axis=1)
        concatenated_pa_df.fillna(0, inplace=True)
        os.makedirs(output_folder, exist_ok=True)
        output_file = os.path.join(output_folder, f"{query_no_ext}.csv")
        concatenated_pa_df.to_csv(output_file, index=True)

def calculating_gator_focal_scores(pa_tables_directory: str, windows_table: str, directory_output: str) -> None:
    ## What function does
     # Calculate the Gator Focal Scores for all the Genbank Windows
    ## Arguments:
     # pa_tables_directory (str): A string for the directory name containing the presence-absence tables
    ## Returns:
     # None: It just generates dataframes containing the Gator Focal Scores.
    for pa_tables_path in glob.glob(os.path.join(pa_tables_directory, '*.csv')):
        pa_tables = pd.read_csv(pa_tables_path, index_col=[0])
        with open(windows_table, 'r') as fh:
            rows = fh.readlines()
            windows_tsv = rows[1:]
            loci = [rows.split()[4] for rows in windows_tsv]
            processed_loci = set()
            max_normal_distribution = None
            for locus in loci:
                if locus in pa_tables.columns and locus not in processed_loci:
                    loci_pa_tables = pa_tables.select_dtypes(include=np.number).columns
                    loci_index = len(loci_pa_tables) + 1
                    num_loci_array = np.arange(1, loci_index)
                    position_loci_in_window = loci_pa_tables.get_loc(locus) + 1
                    stan_dev = (loci_index * 15) / 200
                    probability_function = norm.pdf(num_loci_array, loc=position_loci_in_window, scale=stan_dev)
                    normal_distribution = probability_function / np.max(probability_function)
                    if max_normal_distribution is None:
                        max_normal_distribution = np.zeros_like(normal_distribution)
                    max_normal_distribution = np.maximum(max_normal_distribution, normal_distribution)
                    processed_loci.add(locus)
        loci_pa_tables = pa_tables.select_dtypes(include=np.number).columns
        pa_tables[loci_pa_tables] = pa_tables[loci_pa_tables].multiply(max_normal_distribution)                               
        sum_res = pa_tables[loci_pa_tables].sum(axis=1)
        basename, extension = os.path.splitext(os.path.basename(pa_tables_path))
        pa_tables['Gator Focal Window'] = f"{basename}"
        pa_tables['Gator Focal Sum'] = sum_res
        pa_tables['Gator Focal Score'] = sum_res / sum_res.max()
        pa_tables = pa_tables.sort_values(by='Gator Focal Score', ascending=False)
        os.makedirs(directory_output, exist_ok=True)
        output_filepath = os.path.join(directory_output, f"{basename}_GFS.csv")
        pa_tables.to_csv(output_filepath, index=True)

def concatenating_gator_focal_scores(directory_input: str, directory_output: str) -> None:
    ## What function does
     # Concatenate all the Gator Focal Scores
    ## Arguments:
     # directory_input (str): directory name containing the Gator Focal Scores
     # directory_output (str): directory name wanted to save the concatenated dataframe with the Gator Focal Scores
    ## Returns:
     # None: It just generates the concatenate file
    gfs_dfs = glob.glob(directory_input + '/*.csv')
    extracted_data = []
    for gfs in gfs_dfs:
        all_gfs = pd.read_csv(gfs)
        extracted_data.append(all_gfs.iloc[:, [0, -3, -2, -1]])
    concatenated_data = pd.concat(extracted_data, ignore_index=True)
    concatenated_data.columns = ['Genbank Windows', 'Gator Focal Window', 'Gator Focal Sum', 'Gator Focal Score']
    os.makedirs(directory_output, exist_ok=True)
    output_filepath = os.path.join(directory_output, "concatenated_GFS.csv")
    concatenated_data.to_csv(output_filepath, index=False)

def generating_clustermap(concatenated_GFSs: str, output_dir: str) -> None:
    ## What function does
     # Generates a clustermap for the GFSs
    ## Arguments:
     # concatenated_GFSs (str): A string containing the path for the dataframe containing the GFSs matrix
     # output_dir (str): A string containing the path for the directory wanted to save the clustermap
    ## Returns:
     # None: It just generates the clustermap
    gcfs = pd.read_csv(concatenated_GFSs, delimiter=",")
    gcfs['Genbank Windows'] = gcfs['Genbank Windows'].str.split('-').str[0]
    gcfs['Gator Focal Window'] = gcfs['Gator Focal Window'].str.split('-').str[0]
    matrix_gcfs = gcfs.pivot(index='Genbank Windows', columns='Gator Focal Window', values='Gator Focal Score')
    matrix_gcfs = matrix_gcfs.fillna(0)
    fig, ax = plt.subplots(figsize=(13, 10))
    gcfs_clustermap = sns.clustermap(matrix_gcfs, linewidths=0.1, linecolor='black', annot=False, cmap='viridis', cbar_kws={'label': 'Gator Focal Score', 'shrink': 1}, xticklabels=True, yticklabels=True, tree_kws={'linewidths':1.2})
    gcfs_clustermap.savefig(output_dir, dpi=1200)

def calculating_gator_conservation_percentages(pa_tables_directory: str) -> Dict:
    percentages = {}
    for pa_tables_path in glob.glob(os.path.join(pa_tables_directory, '*.csv')):
        window_genbank= os.path.splitext(os.path.basename(pa_tables_path))[0]
        pa_tables = pd.read_csv(pa_tables_path)
        loci = pa_tables.columns[1:]
        percentages[window_genbank] = {}
        for locus in loci:
            percentage = (pa_tables[locus].sum() / len(pa_tables[locus]))
            percentages[window_genbank][locus] = percentage
    return percentages
    
def opacity_to_hex(opacity: int) -> int:
    ## What function does
     # Convert the opacity values in hexadecimal string format
    ## Arguments:
     # opacity (int): values from 0 to 1
    ## returns:
     # A hexadecimal string format for the opacity (alpha)
    return hex(int(opacity * 255))[2:].zfill(2)

def gator_conservation_plot(genbank_dir: str, flattened_req_hits_by_contig: List, flattened_opt_hits_by_contig: List, percentages: Dict, directory_output: str) -> None:
    ## What function does
     # Generates a gator conservation figure for each window
    ## Arguments:
     # genbank_dir (str): A string containing the directory path for the windows Genbanks
     # flattened_req_hits_by_contig (List): A list containing all the required hits grouped by the same contig
     # flattened_opt_hits_by_contig (List): A list containing all the optional hits grouped by the same contig
     # directory_output (str): A string containing the directory path wanted to save the gator conservation figures
    ## Returns:
     # None: It just generates the gator conservation figures
    genbank_extensions = ['*.gbk', '*.gbff', '*.gb']
    genbank_files = [geno_file for geno_ext in genbank_extensions for geno_file in glob.glob(os.path.join(genbank_dir, geno_ext))]
    req_loci = [req_hits['locus'] for req_hits in flattened_req_hits_by_contig]
    opt_loci = [opt_hits['locus'] for opt_hits in flattened_opt_hits_by_contig]
    for file_path in genbank_files:
        loci_list = []
        window_genome = os.path.splitext(os.path.basename(file_path))[0]
        window = window_genome.split('--')[0]
        for rec in SeqIO.parse(file_path, 'genbank'):
            for feat in rec.features:
                if feat.type == 'CDS':
                    locus = feat.qualifiers['locus_tag'][0]
                    record_length = len(rec)
                    start = int(feat.location.start)
                    end = int(feat.location.end)
                    strand = feat.location.strand
                    annotation = feat.qualifiers['product'][0]
                    loci_list.append((start, end, strand, locus, annotation))
        gv = GenomeViz(
            fig_width = 20,
            fig_track_height = 0.5,
            align_type = "left",
            feature_track_ratio =  1.0,
            link_track_ratio = 1.0,
            tick_track_ratio = 1.0,
            track_spines = False,
            tick_style = "bar",
            plot_size_thr = 0,
            tick_labelsize = 15,
        )
        track = gv.add_feature_track(window, record_length, labelmargin=0.03, linecolor="#333333", linewidth=2)
        for index, cds_tuple in enumerate(loci_list, 1):
            start, end, strand, locus, annotation = cds_tuple
            opacity = 0
            for inner in percentages.values():
                if locus in inner:
                    opacity = inner[locus]
                    break
            alpha = opacity_to_hex(opacity)
            if locus in req_loci:
                color = '#7570B3' + alpha
                label = annotation + ' (req)'
            elif locus in opt_loci:
                color = '#C87137' + alpha
                label = annotation + ' (opt)'
            else:
                color = '#008423' + alpha
                label = ''
            track.add_feature(start, end, strand, facecolor="white")
            track.add_feature(start, end, strand, label=label, facecolor=color, labelsize=10, linewidth=1.5, labelvpos="top")
            track.set_sublabel(text=f"{round(record_length/1000, 2)} Kb", ymargin=1.5)
        os.makedirs(directory_output, exist_ok=True)
        output_filepath = os.path.join(directory_output, f"{window_genome}.svg")
        fig = gv.plotfig()
        gv.set_colorbar(fig, bar_colors=['#7570B3ff', '#C87137ff', '#008423ff'], alpha=1, vmin=0, vmax=100, bar_height=1.4, bar_label="Conservation", bar_labelsize=13)
        fig.savefig(output_filepath, dpi=1200)
        
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
############################# BEGIN MAIN BLOCK ##################################  

## print description
print(DESCRIPTION)

## Parse arguments
args = parse_arguments()

## Make output directory
create_directory(args.out)

## Initialize query protein dictionary
q = {}
q = parse_query_faa(args.required, 'req', q)

if args.optional is not None:
    q = parse_query_faa(args.optional, 'opt', q)

## Create combined required and optional queries fasta
all_merged_queries = tempfile.NamedTemporaryFile(dir=args.out, prefix='all_merged_queries_', suffix='.faa', mode='w', delete=False)
for header in q:
    all_merged_queries.write('>'+header+"\n"+q[header]['seq']+"\n")
all_merged_queries.close()

## hmmsearch against PKS/NRPS modular domains and parse results
with tempfile.NamedTemporaryFile(dir=args.out, prefix='all_merged_queries_', suffix='.domtbl', mode='w', delete=False) as query_modular_domtbl:
    run_hmmsearch(MODULAR_DOMAINS_HMMDB, all_merged_queries.name, query_modular_domtbl.name, args.threads, args.e_value)
    modular_domain_hit = parse_modular_domtblout(query_modular_domtbl.name, q)
    print("[2]" + print_datetime(), 'Screening for Modular Domains in Required and Optional Proteins')
    parse_modular_domain_query_hits(modular_domain_hit, q)

## Create nonmodular fasta (to be subsequent diamond query)
nonmodular_queries = tempfile.NamedTemporaryFile(dir=args.out, prefix='non_modular_queries_', suffix='.faa', mode='w', delete=False)
for header in q:
    if not q[header]['is_nrps'] and not q[header]['is_pks']:
        nonmodular_queries.write('>'+header+"\n"+q[header]['seq']+"\n")
nonmodular_queries.close()

## parse genome modular search
modular_domain_hit = parse_modular_domtblout(args.modular_domtblout, q)
print("[3]" + print_datetime(), 'Parsing  Modular Domains from the Genomes Database')
q = parse_modular_domain_genome_hits(modular_domain_hit, q)

## diamond against nonmodular fasta
with tempfile.NamedTemporaryFile(dir=args.out, prefix='dmnd_out_', suffix='.txt', mode='w', delete=False) as dmnd_out:
    print("[4]" + print_datetime(), 'Running Diamond')
    run_diamond(nonmodular_queries.name, args.dmnd_database, dmnd_out.name, args.threads, args.query_cover, args.identity) 

## parse diamond search
parse_diamond_search(dmnd_out.name, q)

## grouping query proteins hits by contig
req_hits_by_contig, all_required_proteins = grouping_user_req_opt_proteins_by_contig(q, 'req')
opt_hits_by_contig, all_optional_proteins = grouping_user_req_opt_proteins_by_contig(q, 'opt')

## checking distance between loci and presence of all required user proteins
print("[5]" + print_datetime(), 'Checking Presence for All Required User Proteins and Evaluating Intergenic Distance Betweeen Loci')
final_window = checking_distance_bw_loci(req_hits_by_contig, all_required_proteins)

## generating Genbank files
print("[5.2]" + print_datetime(), 'Generating Genbank Files for the Windows')
generating_windows_genbanks(final_window, f"{args.out}/windows_genbanks")

## making table for windows and optional hits                                                                                                                                            
flattened_req_hits_by_contig, flattened_opt_hits_by_contig = making_windows_and_optional_table_hits(f"{args.out}/windows_genbanks", req_hits_by_contig, opt_hits_by_contig, f'{args.out}/windows_genbanks/windows.tsv')

## making protein database for the windows
with tempfile.NamedTemporaryFile(dir=args.out, prefix='allvall_proteins_', suffix='.faa', mode='w', delete=False) as allvall_proteins:
    print("[6]" + print_datetime(), 'Generating the Protein Database from the Windows Genbanks')
    dbfaa_from_gb_dir(f"{args.out}/windows_genbanks", allvall_proteins.name)

## making diamond database for the windows
with tempfile.NamedTemporaryFile(dir=f"{args.out}/windows_genbanks", prefix='dmnd_db_', suffix='.dmnd', mode='w', delete=False) as dmnd_db:
    print("[7]" + print_datetime(), 'Creating Diamond Database')
    create_diamond_database(allvall_proteins.name, dmnd_db.name)

## running diamond for all vs all clusters
with tempfile.NamedTemporaryFile(dir=f"{args.out}/windows_genbanks", prefix='dmnd_out_allvall_', suffix='.txt', mode='w', delete=False) as dmnd_out_allvall:
    print("[8]" + print_datetime(), 'Running Diamond against the Windows Database')
    run_diamond(allvall_proteins.name, allvall_proteins.name, dmnd_out_allvall.name, args.threads, args.query_cover, args.identity)

## making gene-level presence-absence tables
print("[9]" + print_datetime(), 'Making Gene-Level Presence-Absence Tables')
pa_tables_from_dmnd_windows_faa(dmnd_out_allvall.name, f"{args.out}/presence_absence")

## calculating gator scores for all vs all windows
print("[10]" + print_datetime(), 'Calculating Gator Focal Scores (GFS) for All vs All Windows')
calculating_gator_focal_scores(f"{args.out}/presence_absence", f"{args.out}/windows_genbanks/windows.tsv", f"{args.out}/gator_scores")

## concatening all gator forcal scores
print("[11]" + print_datetime(), 'Concatening All GFS')
concatenating_gator_focal_scores(f"{args.out}/gator_scores", f"{args.out}/concatenated_scores")

## making the clustermap for the gator focal scores concatenate
print("[12]" + print_datetime(), 'Generating the Clustermap for the GFS concatenate')
generating_clustermap(f"{args.out}/concatenated_scores/concatenated_GFS.csv", f"{args.out}/concatenated_scores/clustermap_GFSs.svg")

## calculating locus conservation percentage
percentages = calculating_gator_conservation_percentages(f"{args.out}/presence_absence")

## making the gator conservation figures                                                                                                                     
print("[13]" + print_datetime(), 'Generating the Gator Conservation Figures')
gator_conservation_plot(f"{args.out}/windows_genbanks/", flattened_req_hits_by_contig, flattened_opt_hits_by_contig, percentages, f"{args.out}/gator_conservation_plots")

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

print("[14]" + print_datetime(), f"Elapsed time: {ftime:.1f} {time_unit}")
