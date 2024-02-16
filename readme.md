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

## Instalation

Installation can be performed via conda and should take ~5 minutes

```bash
# 1. clone the repository:
git clone https://github.com/chevrettelab/gator-gc.git
cd gator-gc/

# 2. create conda environment using yaml file and activate it. Use mamba instead of conda for faster installation:
conda env create -f gator-gc_env.yml or mamba env create -f gator-gc_env.yml
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
