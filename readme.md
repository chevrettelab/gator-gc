# GATOR-GC: Genomic Assessment Tool for Orthologous Regions and Gene Clusters
<img src="images/gator_gc_logo.jpg">

## Overview

GATOR-GC is a user-friendly algorithm designed for targeted exploration of BGC and genomic islands diversity. It focuses on key biosynthetic enzymes and offers flexibility in defining the taxonomic scope of the analysis. Unlike methods relying on arbitrary cutoffs, GATOR-GC establishes BGC boundaries based on evolutionary principles and implements an enzyme-aware scoring system for assessing BGC-BGC distances, moving beyond a binary presence-absence framework. This approach enhances the tool's capability to identify and prioritize novelty, effectively mapping biosynthetic diversity into distinct groups

## Authors

- Developer: José D. D. Cediel-Becerra
- Code reviewers: Valérie de Crécy-Lagard and Marc G. Chevrette
- Afiliation: Microbiology & Cell Science Deparment, University of Florida
- Please contact me at jcedielbecerra@ufl.edu if you have any questions

## Features

- **Targeted Search:** Efficiently search for specific key (required) and tailoring (optional) enzymes within BGCs and/or genes in genomic islands
- **Screening for modular domains:** User input proteins file will be screened to identify modular domains (NRPSs and PKSs)
- **Modular domains are treated differently:** NRPSs and PKSs will be searched using the HMMs profiles used by antiSMASH
- **Customizable Parameters:** Fine-tune search parameters (i.e., required/optional proteins, assembly with ALL required proteins, user-defined intergenic distance between required proteins ) to enhance specificity to define the GATOR windows
- **GATOR Scores:** Novel score to compare GATOR windows with specific focal windows implementing an enzyme-aware scoring system 
- **GATOR Conservation:** Gene cluster figure where the color of the genes are based on required and optional proteins, and their transparency depends on the  presence of the genes in the GATOR windows 
- **GATOR Neighborhoods:** Genomic neighborhoods containing all the GATOR windows, which are sorted based on the GATOR scores. Each GATOR window will have a neighborhood figure. Homology between genes are ilustrated by a gray bar. 

## Usage
Before using GATOR-GC, we need run PRE-GATOR-GC in order to get the proteins, dmnd database and the modular domtblout.

```
pre-gator-gc --genomes_dir /my_genomes --proteins name_for_proteins.faa --dmnd_database name_for_dmnd_database.dmnd --modular_domtblout name_for_modular_domtblout.txt --e_value 1e-4 --out my_output_folder
optional arguments:
  -h, --help            show this help message and exit
  --genomes_dir GENOMES_DIR
                        Directory name containing the Genbanks (*.gbff/*.gbk/*.gb) genomes.
  --proteins PROTEINS   Name for the precompiled protein database (.faa) of genomes_dir.
  --dmnd_database DMND_DATABASE
                        Name for the precompiled diamond database (.dnmd).
  --e_value E_VALUE     E-value threshold  wanted for hmmsearch (default: 1e-4).
  --modular_domtblout MODULAR_DOMTBLOUT
                        Name for the precomputed hmmsearch domain table for modular domains
  --threads THREADS     CPUs wanted for hmmsearch (default: all available).
  --out OUT             Output directory name that will contain the proteins, the dmnd_database, and the modular domtblout table.
``` 

```
gator-gc --required req.faa --optional opt.faa --genomes_dir my_genomes/ --proteins my_proteins.faa --dmnd_database my_dmnd_database.dmnd --modular_domtblout my_modular_domtblout.txt --intergenic_distance 86 --window_extension 10 --out my_folder_output

optional arguments:
  -h, --help            show this help message and exit
  --required REQUIRED   Query protein fasta file containing required proteins.
  --optional OPTIONAL   Query protein fasta file containing optional proteins.
  --genomes_dir GENOMES_DIR
                        Directory name containing the Genbanks (*.gbff/*.gbk/*.gb) genomes.
  --proteins PROTEINS   Precompiled protein database (.faa) of genomes_dir to search against.
  --dmnd_database DMND_DATABASE
                        Precompiled diamond database (.dnmd) to search against.
  --modular_domtblout MODULAR_DOMTBLOUT
                        Precomputed hmmsearch domain table for modular domains against the precompiled protein database.
  --threads THREADS     CPUs wanted for diamond search and hmmsearch (default: all CPUs available).
  --query_cover QUERY_COVER
                        Protein percent query cover for diamond search (default: 70).
  --identity IDENTITY   Protein percent identity for diamond search (default: 35).
  --e_value E_VALUE     E-value threshold  wanted for hmmsearch (default: 1e-4).
  --intergenic_distance INTERGENIC_DISTANCE
                        Distance in kilobases between required genes (default: 86 kb)
  --window_extension WINDOW_EXTENSION
                        Extension in kilobases from the start and end positions of the windows (default: 10 kb)
  --out OUT             Output directory name that will have GATOR-GC results
```

## Example 

```
python pre-gator-gc.py --genomes_dir example/genomes/ --proteins proteins.faa --dmnd_database database.dmnd --modular_domtblout modular.domtblout --threads 4 --out example/pre_gator_data
```
