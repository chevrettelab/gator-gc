# GATOR-GC: Genomic Assessment Tool for Orthologous Regions and Gene Clusters
<img src="images/gator_gc_logo.jpg">

## Overview

GATOR-GC is a user-friendly algorithm designed for targeted exploration of BGC and genomic islands diversity. It focuses on key biosynthetic enzymes and offers flexibility in defining the taxonomic scope of the analysis. Unlike methods relying on arbitrary cutoffs, GATOR-GC establishes BGC boundaries based on evolutionary principles and implements an enzyme-aware scoring system for assessing BGC-BGC distances, moving beyond a binary presence-absence framework. This approach enhances the tool's capability to identify and prioritize novelty, effectively mapping biosynthetic diversity into distinct groups

## Features

- **Targeted Search:** Efficiently search for specific key (required) and tailoring (optional) enzymes within BGCs and/or genes in genomic islands
- **Screening for modular domains:** User input proteins file will be screened to identify modular domains (NRPSs and PKSs)
- **Modular domains are treated differently:** NRPSs and PKSs will be searched using the HMMs profiles used by antiSMASH
- **Customizable Parameters:** Fine-tune search parameters (i.e., required/optional proteins, assembly with ALL required proteins, user-defined intergenic distance between required proteins ) to enhance specificity to define the GATOR windows
- **GATOR Scores:** Novel score to compare GATOR windows with specific focal windows implementing an enzyme-aware scoring system 
- **GATOR Conservation:** Gene cluster figure where the color of the genes are based on required and optional proteins, and their transparency depends on the  presence of the genes in the GATOR windows 
- **GATOR Neighborhoods:** Genomic neighborhoods containing all the GATOR windows, which are sorted based on the GATOR scores. Each GATOR window will have a neighborhood figure. Homology between genes are ilustrated by a gray bar. 

## Authors

- Developer: José D. D. Cediel-Becerra
- Code reviewers: Valérie de Crécy-Lagard and Marc G. Chevrette
- Afiliation: Microbiology & Cell Science Deparment, University of Florida
- Please contact me at jcedielbecerra@ufl.edu if you have any questions
