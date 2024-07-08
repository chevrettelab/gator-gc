# GATOR-GC: Genomic Assessment Tool for Orthologous Regions and Gene Clusters
<img src="images/gator_gc_logo.jpg">

## Overview
<img src="images/gator_gc_alogorithm.png">

GATOR-GC is a user-friendly algorithm designed for targeted exploration of BGC and genomic islands diversity. It focuses on key biosynthetic enzymes and offers flexibility in defining the taxonomic scope of the analysis. Unlike methods relying on arbitrary cutoffs, GATOR-GC establishes BGC boundaries based on evolutionary principles and implements an enzyme-aware scoring system for assessing BGC-BGC distances, moving beyond a binary presence-absence framework. This approach enhances the tool's capability to identify and prioritize novelty, effectively mapping biosynthetic diversity into distinct groups

## Authors

- Code, design, and development: José D. D. Cediel-Becerra,  Marc G. Chevrette
- Afiliation: Microbiology & Cell Science Deparment, University of Florida
- Please contact José at jcedielbecerra@ufl.edu if there are any issues

## Features

- **Targeted Search:** Conduct targeted searches for essential key enzymes and optional tailor enzymes within Biosynthetic Gene Clusters (BGCs) and genomic islands, streamlining the discovery process
- **Modular Domain Screening:** Automatically screen user-provided protein files to identify critical modular domains, such as Non-Ribosomal Peptide Synthetases (NRPSs) and Polyketide Synthases (PKSs), using state-of-the-art HMM profiles from antiSMASH for unparalleled precision
- **Customizable Parameters:**  Customize search parameters to include required and optional proteins, ensure complete assembly with all necessary proteins, and define specific distances between required proteins. This customization enhances the specificity of GATOR window identification, offering tailored analysis to meet research needs
- **Dereplication for Gator Windows:** Remove identical Gator Windows to avoid redundant calculations in subsequent steps, while keeping track of which windows are duplicates
- **GATOR Focal Scores:** Employ a novel enzyme-aware scoring system to accurately compare GATOR windows against targeted focal windows. This approach ensures precise evaluation of genomic contexts and enzyme functionalities
- **GATOR Conservation:**  Generate dynamic gene cluster diagrams that visually differentiate between required and optional proteins using color coding. Transparency levels indicate the gene's presence within GATOR windows, providing a clear visual representation of gene conservation
- **GATOR Neighborhoods:** Visualize each GATOR window's genomic neighborhood with organized tracks based on GATOR focal scores. Homology between genes is intuitively illustrated with gray bars, facilitating easy understanding of genetic relationships and conservation

## Installation

Installation can be performed via conda and should take ~5 minutes

```bash
# 1. clone the repository:
git clone https://github.com/chevrettelab/gator-gc.git
cd gator-gc/

# 2. create conda environment using yaml file and activate it. Use mamba instead of conda for faster installation:
   # with conda:
   conda env create -f gator-gc_env.yml
   conda activate gator-gc

   # or with mamba:
   mamba env create -f gator-gc_env.yml
   conda activate gator-gc	

# 3. install the python package
pip install .
```

## Usage

Before utilizing gator-gc, it's necessary to execute pre-gator-gc to obtain the the diamond database, and the modular domain table output.

```
optional arguments:
  -h, --help  show this help message and exit
  -v          Enable verbose output (default: False).

Input Options:
  -g  [ ...]  Directories containing GenBank files (*.gbff/*.gbk/*.gb). You can specify multiple directories separated by spaces. Directories can be specified with or without wildcards.

HMMER Options:
  -e          E-value threshold for HMMER hmmsearch (default: 1e-4).
  -t          Number of CPU threads to use for hmmsearch (default: all available threads).

Output Options:
  -o          Directory where the gator databases (protein,  DIAMOND, and modular domtblout databases) will be saved.
``` 

Now we can run gator-gc to identify the gator windows. 

```
optional arguments:
  -h, --help  show this help message and exit
  -v          Enable verbose output. (Default: False)

Input Options:
  -rq         Path to the query protein FASTA file containing required proteins.
  -op         Path to the query protein FASTA file containing optional proteins.
  -g  [ ...]  Directory containing the Genbank files (*.gbff/*.gbk/*.gb).You can specify multiple directories separated by spaces. Directories can be specified with or without wildcards.
  -d          Directory containing the PRE-GATOR-GC databases (.dmnd and .domtblout files). (Required)

Diamond Options:
  -t          Number of CPUs to use for diamond search and hmmsearch. (Default: all available CPUs)
  -qc         Minimum percent query cover for diamond search. (Default: 70)
  -idt        Minimum percent identity for diamond search. (Default: 35)

HMMER Options:
  -e          E-value threshold for hmmsearch. (Default: 1e-4)

GATOR-GC Options:
  -rd         Maximum distance in kilobases between required genes to define a gator window. (Default: 86 kb)
  -we         Extension in kilobases from the start and end positions of the gator windows. (Default: 10 kb)

Output Options:
  -o          Directory to save GATOR-GC results.
  -nc         Disable creation of GATOR conservation figures.
  -nn         Disable creation of GATOR neighborhoods figures.
```

## Example

Let's explore an example involving the enzymes responsible for the production of prodigiosin and prodigiosin-like compounds, categorizing them as required and optional. Subsequently, we will define the taxonomic scope of the search, focusing  on a couple genomes from *Streptomyces*, *Serratia*, *Pseudoalteromonas*, and *Hallela*.
The inital step involves using pre-gator-gc:

```
pre-gator-gc --genomes_dir example/genomes/ --out output_name_gator_databases -v

```
after running pre-gator-gc on this example and using the verbose flag, we will see:

```
                                                                                                                                                                     
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
Afiliation: Microbiology & Cell Science Deparment, University of Florida                                                                              
Please contact José at jcedielbecerra@ufl.edu if you have any issues                                                                                       
Version: v0.9.0
[1] - 2024-07-08 14:51:12,501 - INFO - The output_name_gator_databases directory was created successfully.
[2] - 2024-07-08 14:51:12,502 - INFO - Total genome files found: 10
[3] - 2024-07-08 14:51:45,238 - INFO - Successfully created the gator protein database to output_name_gator_databases/output_name_gator_databases.faa
[4] - 2024-07-08 14:51:45,840 - INFO - Successfully created the gator DIAMOND database to output_name_gator_databases/output_name_gator_databases.dmnd
[5] - 2024-07-08 14:51:59,790 - INFO - Successfully created the gator domtblout database to output_name_gator_databases/output_name_gator_databases.domtblout
[6] - 2024-07-08 14:51:59,790 - INFO - Execution time: 47.3 seconds
```

With this foundation in place, we can proceed to use gator-gc to search gene clusters (a.k.a GATOR windows) harboring the wanted genes. We can do it running these arguments:

```
gator-gc -rq example/proteins/req.faa -op example/proteins/opt.faa -g example/genomes/ -d output_name_gator_databases/ -v -o results

```
after running gator-gc on this example and using the verbose flag, we will see:

```
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
Afiliation: Microbiology & Cell Science Deparment, University of Florida
Please contact José at jcedielbecerra@ufl.edu if you have any issues
Version:v0.9.0
[1] - 2024-07-08 14:36:16,037 - INFO - The prodigiosin directory was created successfully.
[2] - 2024-07-08 14:36:16,038 - INFO - Total genome files found: 10
[3] - 2024-07-08 14:36:17,247 - INFO - hmmsearch completed successfully. Output written to /mnt/clab1/GATOR-GC/20230922-gator_dev/gator-gc/prodigiosin/all_merged_queries_muqt3vf6.domtbl
[4] - 2024-07-08 14:36:17,248 - INFO - Parsing modular domtblout completed successfully. Found 1 proteins with domain hits.
[5] - 2024-07-08 14:36:17,248 - INFO - Protein dictionary  with boolean values set for is_nrps and is_pks completed successfully.
[6] - 2024-07-08 14:36:17,249 - INFO - Found Diamond database file: example/output_name_gator_databases/output_name_gator_databases.dmnd
[7] - 2024-07-08 14:36:17,249 - INFO - Found Domtblout database file: example/output_name_gator_databases/output_name_gator_databases.domtblout
[8] - 2024-07-08 14:36:17,268 - INFO - Parsing modular domtblout completed successfully. Found 582 proteins with domain hits.
[9] - 2024-07-08 14:36:17,268 - INFO - No need to integrate Hybrids, NRPS, and PKS hits into protein_dict.
[10] - 2024-07-08 14:36:17,269 - INFO - Running Diamond with /mnt/clab1/GATOR-GC/20230922-gator_dev/gator-gc/prodigiosin/non_modular_queries_8rw38rm_.faa against example/output_name_gator_databases/output_name_gator_databases.dmnd
[11] - 2024-07-08 14:36:20,284 - INFO - Successfully parsed diamond hits from /mnt/clab1/GATOR-GC/20230922-gator_dev/gator-gc/prodigiosin/dmnd_out_9b0sy8th.txt and integrated to protein_dict
[12] - 2024-07-08 14:36:20,285 - INFO - Number of Windows Found: 6
[13] - 2024-07-08 14:36:20,286 - INFO - Processing gator windows, writing GenBank files for each window, FASTA database, and tracks metadata for visualization and analysis
[14] - 2024-07-08 14:36:33,720 - INFO - Successfully wrote gator table with metadata to prodigiosin/windows_genbanks/gator_hits.tsv
[15] - 2024-07-08 14:36:33,876 - INFO - Starting deduplication process for gator windows.
[16] - 2024-07-08 14:36:46,211 - INFO - All windows have been processed and deduplicated.
[17] - 2024-07-08 14:36:46,212 - INFO - Number of unique gator windows after deduplication: 5
[18] - 2024-07-08 14:36:46,212 - INFO - Unique computations written to prodigiosin/deduplication_data/unq_comp.tsv
[19] - 2024-07-08 14:36:48,703 - INFO - Gene level presence absence tables for deduplicated gator windows created sucessfully in prodigiosin/presence_absence
[20] - 2024-07-08 14:36:48,713 - INFO - Gator focal scores tables for deduplicated gator windows created sucessfully in prodigiosin/gator_scores
[21] - 2024-07-08 14:36:48,727 - INFO - Concatenated of gator focal scores across deduplicated gator windows created sucessfully in prodigiosin/concatenated_scores/concatenated_gfs.csv
[22] - 2024-07-08 14:36:49,723 - INFO - Clustered heatmap with gator focal scores saved to prodigiosin/concatenated_scores/clustermap_gfs.svg
[23] - 2024-07-08 14:36:55,328 - INFO - All gator conservation figures have been generated successfully.
[24] - 2024-07-08 14:37:19,032 - INFO - All gator neighborhood figures have been generated successfully.
[25] - 2024-07-08 14:37:19,033 - INFO - Execution time: 1.05 minutes
```

Six putative regions containing the required genes (pigC, pigG, and PigI) were identified, with a specified maximum distance between them of <86 kb (this default value can be adjusted using --rd. Additionally, a window extension of 10 Kb was applied (this default value can be modified using -we).

The output folder "results/" for gator-gc includes the following subfolders:

1. windows_genbanks: This folder contains six GenBank files, each generated for a specific window.
2. presence_absence: Within this folder, six gene-level presence-absence tables are present for each gator window in CSV format. The tables are named based on the window number and the corresponding GenBank filename. In these tables, the columns represent the locus_tags in the gator focal window, and the rows correspond to the gator windows.
3. gator_scores: This folder houses six tables with a structure similar to the presence-absence tables. However, instead of binary numbers, normal distributions were applied to each required and optional protein. The highest distribution value is set to 1, and these values were multiplied by the presence-absence values. The resulting sums for each gator window were then normalized based on the maximum sum score (gator focal window), generating the gator focal scores. The tables are named based on the window number and GenBank filenames.
4. concatenated_scores: This folder contains a "concatenated_GFS.csv" file, which consolidates all the gator focal scores. Additionally, a "clustermap_GFSs" file displays the distribution of gator windows based on the gator focal scores using a heatmap and a dendrogram. This figure is a high-quality vectorized representation.
<img src="example/prodigiosin_results/concatenated_scores/clustermap_GFSs.svg">

5. gator_conservation_plots: This folder includes six high-quality vectorized figures, each corresponding to a gator focal window. These figures display the gene organization and cluster size. An example of this figure is shown below:

<img src="example/prodigiosin_results/gator_conservation_plots/window_1--GCF_017298755.1_ASM1729875v1_genomic.svg">

In these figures, genes are color-coded according to the user's input, with required genes represented in purple, optional genes in orange, and genes not present in the protein queries files in green. The transparency of each gene's color is determined by the conservation of homologous genes found in the gator windows. Additionally, labels in the figures include the user's input headers categorized as required and optional, along with annotations present in the GenBank. Genes with red edges indicate contig edges.

6. gator_neighborhoods_plots: There are six high-quality vectorized figures showcasing genome neighborhoods. Each file corresponds to a gator focal window (top track), with the remaining gator windows sorted based on the gator focal scores. This arrangement positions the second top track as the most similar gator window, and the last track at the bottom represents the most dissimilar gator window. The genomic organization of the gator windows is flipped based on the first required gene in the gator focal window, ensuring that the genomic organization aligns with that of the gator focal window. Homology rails are displayed, depicting the diamond protein alignment position hits in the genes. These figures have the following appearance:

<img src="example/prodigiosin_results/gator_neighborhoods_plots/window_1--GCF_017298755.1_ASM1729875v1_genomic_neighboorhoods.svg">
