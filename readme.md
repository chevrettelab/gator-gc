# GATOR-GC: Genomic Assessment Tool for Orthologous Regions and Gene Clusters
<img src="images/gator_gc_logo.jpg">

## Overview
<img src="images/gator_gc_alogorithm.png">
GATOR-GC is a user-friendly algorithm designed for targeted exploration of BGC and genomic islands diversity. It focuses on key biosynthetic enzymes and offers flexibility in defining the taxonomic scope of the analysis. Unlike methods relying on arbitrary cutoffs, GATOR-GC establishes BGC boundaries based on evolutionary principles and implements an enzyme-aware scoring system for assessing BGC-BGC distances, moving beyond a binary presence-absence framework. This approach enhances the tool's capability to identify and prioritize novelty, effectively mapping biosynthetic diversity into distinct groups

## Authors

- Developer: José D. D. Cediel-Becerra
- Code reviewers: Valérie de Crécy-Lagard and Marc G. Chevrette
- Afiliation: Microbiology & Cell Science Deparment, University of Florida
- Please contact José at jcedielbecerra@ufl.edu if there are any issues

## Features

- **Targeted Search:** Conduct targeted searches for essential key enzymes and optional tailor enzymes within Biosynthetic Gene Clusters (BGCs) and genomic islands, streamlining the discovery process
- **Modular Domain Screening:** Automatically screen user-provided protein files to identify critical modular domains, such as Non-Ribosomal Peptide Synthetases (NRPSs) and Polyketide Synthases (PKSs), using state-of-the-art HMM profiles from antiSMASH for unparalleled precision
- **Customizable Parameters:**  Customize search parameters to include required and optional proteins, ensure complete assembly with all necessary proteins, and define specific distances between required proteins. This customization enhances the specificity of GATOR window identification, offering tailored analysis to meet research needs
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
  -h, --help            show this help message and exit

Input Options:
  --genomes_dir GENOMES_DIR [GENOMES_DIR ...]
                        Directory(ies) name(s) containing the Genbanks (*.gbff/*.gbk/*.gb) genomes.

HMMER Options:
  --e_value E_VALUE     E-value threshold  wanted for hmmsearch (default: 1e-4).
  --threads THREADS     CPUs wanted for hmmsearch (default: all available).

Output Options:
  --out OUT             Output directory name that will contain the proteins, the dmnd_database, and the modular domtblout table.
``` 

Now we can run gator-gc to identify the gator windows. 

```
optional arguments:
  -h, --help            show this help message and exit

Input Options:
  --required REQUIRED   Query protein fasta file containing required proteins.
  --optional OPTIONAL   Query protein fasta file containing optional proteins.
  --genomes_dir GENOMES_DIR [GENOMES_DIR ...]
                        Directory containing the Genbanks (*.gbff/*.gbk/*.gb) genomes.
  --gator_databases GATOR_DATABASES
                        Folder containing the pre-gator-gc databases (.faa and .domtblout)

Diamond Options:
  --threads THREADS     CPUs wanted for diamond search and hmmsearch (default: all CPUs available).
  --query_cover QUERY_COVER
                        Protein percent query cover for diamond search (default: 70).
  --identity IDENTITY   Protein percent identity for diamond search (default: 35).

HMMER Options:
  --e_value E_VALUE     E-value threshold  wanted for hmmsearch (default: 1e-4).

GATOR-GC Options:
  --required_distance REQUIRED_DISTANCE
                        Distance in kilobases between required genes (default: 86 kb)
  --window_extension WINDOW_EXTENSION
                        Extension in kilobases from the start and end positions of the windows (default: 10 kb)

Output Options:
  --out OUT             Output directory name that will have GATOR-GC results
  --no_conservation_figs
                        Gator conservation figures will not be created
  --no_neighborhoods_figs
                        Gator neighborhoods figures will not be created

```

## Example

Let's explore an example involving the enzymes responsible for the production of prodigiosin and prodigiosin-like compounds, categorizing them as required and optional. Subsequently, we will define the taxonomic scope of the search, focusing  on a couple genomes from *Streptomyces*, *Serratia*, *Pseudoalteromonas*, and *Hallela*.
The inital step involves using pre-gator-gc:

```
pre-gator-gc --genomes_dir example/genomes/ --out output_name_gator_databases
```

With this foundation in place, we can proceed to use gator-gc to search gene clusters (a.k.a GATOR windows) harboring the wanted genes. We can do it running this arguments:

```
gator-gc --required example/proteins/req.faa --optional example/proteins/opt.faa --genomes_dir example/genomes/ --gator_databases output_name_gator_databases --out results

```
after running gator-gc on this example, we will see:

```
[1][2024-02-16 10:26:54] The results Directory was Created
[2][2024-02-16 10:26:55] Screening for Modular Domains in Required and Optional Proteins
[3][2024-02-16 10:26:55] Parsing  Modular Domains from the Genomes Database
[4][2024-02-16 10:26:55] Running Diamond
[5][2024-02-16 10:26:56] Checking Presence for All Required User Proteins and Evaluating Intergenic Distance Betweeen Loci
[5.1][2024-02-16 10:26:56] Number of Windows Found: 6
[5.2][2024-02-16 10:26:56] Generating Genbank Files for the Windows
[6][2024-02-16 10:27:00] Generating the Protein Database from the Windows Genbanks
[7][2024-02-16 10:27:01] Creating Diamond Database
[8][2024-02-16 10:27:01] Running Diamond against the Windows Database
[9][2024-02-16 10:27:01] Making Gene-Level Presence-Absence Tables
[10][2024-02-16 10:27:01] Calculating Gator Focal Scores (GFS) for All vs All Windows
[11][2024-02-16 10:27:01] Concatening All GFS
[12][2024-02-16 10:27:01] Generating the Clustermap for the GFS concatenate
[13][2024-02-16 10:27:02] Generating the Gator Conservation Figures
[14][2024-02-16 10:27:03] Generating the Windows Neighborhoods Figures
[15][2024-02-16 10:27:13] Elapsed time: 18.6 seconds
```
Six putative regions were found containing the required genes (pigC, pigG and PigI) with a required distance between them < 86 kb (this is default, it can be changed using --required_distance), and with a window extension of 10 Kb (this is also default, it can be changed using --window_extension).

The output folder results/ for gator-gc have these subfolders:
1. windows_genbanks/ have six GenBank files were generated for each window.
2. presence_absence/ have six gene-level presence absence tables for each gator window in CSV format. These tables are named based on the window numer and the GenBank filename. In these tables the columns are the locus_tags in the gator focal window, and the rows are the gator widows
3. gator_scores/ having six tables with the same structure than the presence_absence tables, but instead of binary numbers, normal distributions were applied for each required and optional proteins (highest distribution value, equals to 1) and these values were multiplied for the presence absence values, and summed for each gator window. These summ scores were normalized based on the maximum summ socre (gator focal window) and this generates the gator focal scores. Gator windows are sorted in the tables based on the gator focal scores. The gator scores tables are names based on the window number and GenBank filenames.
4. concatenated_scores/ have a concatenated_GFS.csv file containing all the gator focal scores, and a clustermap_GFSs showing the distribution of the gator windows based on the gator focal scores using a heatmpa and a dendogram. This is a high-quality vectorized figure.

<img src="example/prodigiosin_results/concatenated_scores/clustermap_GFSs.svg">

5. gator_conservation_plots/ contains six high-quality vectorized figures for each gator focal window displaying the gene organization and the cluster size. This figure looks like this:

<img src="example/prodigiosin_results/gator_conservation_plots/window_1--GCF_017298755.1_ASM1729875v1_genomic.svg">

where genes are color coded based on the user's input (required as purple, optional as orange) and genes no present in the proteins queries files as green. Note that the gene's color have transparency based on the conservation of homologous genes in the gator windows found. Also, the labels contains the user's headers inputs categorized as required and optional, as well as the annotation present in the GenBank. When genes have red edges means that it is contig edge.

6. gator_neighborhoods_plots/ have six high-quality vectorized figures displaying genome neighborhoods, where each file corresponds for a gator focal window (top track), and the rest gator windows are sorted based on the gator focal scores, meaning that the second top track is the most similar gator window and the last track at the bottom the most disimilar gator window. The genomic organization of the gator windows were flipped based on the first required gene in the gator focal window to ensure that the genomic organization for the gator windows is the same than the gator focal window. Homology rails are displayed based on the diamond protein alignment position hits in the genes. These figures looks like this:

<img src="example/prodigiosin_results/gator_neighborhoods_plots/window_1--GCF_017298755.1_ASM1729875v1_genomic_neighboorhoods.svg">
