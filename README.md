# EPS-Operon-Mining-and-Visualization
#This repository contains custom Python scripts developed for the identification, extraction, and comparative visualization of exopolysaccharide (EPS) biosynthesis gene clusters from bacterial genomes. This workflow bridges the gap between raw functional annotations (e.g., dbCAN, NCBI RefSeq)
## Scripts Included

### 1. `find_eps_clusters.py`
A parsing script that identifies physical EPS gene clusters (operons) within a genome. 
* **Input:** A target `.gff` genomic annotation file and a list of identified Glycosyltransferases (e.g., from dbCAN).
* **Function:** Scans for co-localized structural genes (GTs) and essential EPS translocation/regulatory markers (e.g., flippases, polymerases, *wzx/wzy* components) within a user-defined base-pair window (default: 15,000 bp).
* **Output:** Terminal readout of candidate operon coordinates and functional contents.

### 2. `extract_operons.py`
A surgical extraction tool for isolating specific gene clusters while preserving genetic context.
* **Input:** A master RefSeq GenBank flat file (`.gbff`).
* **Function:** Uses Biopython to slice out targeted operon coordinates with a customizable flanking buffer (e.g., 2,000 bp) to preserve upstream/downstream regulatory context.
* **Output:** Standardized `.gbk` files ready for comparative alignment.

## Prerequisites
To run these scripts, you will need a standard bioinformatics environment (like Ubuntu) with the following installed:
* Python 3.x
* [Biopython](https://biopython.org/) (`pip install biopython`)
* [clinker](https://github.com/gamcil/clinker) (`pip install clinker`) - *Used for the downstream visualization.*

## Usage
1. Parse the genome to find your clusters:
   ```bash
   python find_eps_clusters.py eps_builders_GT.txt combined_genomic.gff
   
2.Update the coordinates in extract_operons.py based on your results, then run:

   python extract_operons.py
   
3.Generate the comparative visual alignment:

   clinker candidate1.gbk candidate2.gbk -p final_eps_plot.html
