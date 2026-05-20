#!/usr/bin/env python3
import sys
from Bio import SeqIO

def main():
    gbff_file = "genomic.gbff"
    
    print(f"[*] Loading {gbff_file} (this may take a few seconds)...")
    try:
        records = list(SeqIO.parse(gbff_file, "genbank"))
    except FileNotFoundError:
        print(f"Error: Could not find '{gbff_file}' in this directory.")
        sys.exit(1)
        
    print(f"[*] Found {len(records)} sequence(s) in the file.")

    chrom = None
    # Print the IDs so we can see what Biopython is actually reading
    for rec in records:
        print(f"    -> Scanning sequence ID: {rec.id} (Length: {len(rec.seq)} bp)")
        if "CP025781" in rec.id or "NZ_CP025781" in rec.id:
            chrom = rec
            break
            
    # Fallback: If name matching fails, just grab the biggest contig
    if not chrom:
        print("\n[*] Exact ID match failed. Grabbing the largest sequence instead...")
        chrom = max(records, key=lambda r: len(r.seq))
        print(f"[*] Selected {chrom.id} as the main chromosome.")

    print("\n[*] Slicing Operons...")
    # Slice Candidate 1
    c1_start, c1_end = 2514365 - 2000, 2542732 + 2000 #change the location co-ordinates according to the output of find_eps_clusters.py
    cand1 = chrom[c1_start:c1_end]
    cand1.id = "Candidate_1_EpsFG"
    cand1.name = "EpsFG"
    cand1.description = "Classic wzx/wzy dependent EPS Operon"
    SeqIO.write(cand1, "candidate1.gbk", "genbank")
    print(f"[+] Candidate 1 extracted: {c1_start} to {c1_end}")

    # Slice Candidate 2
    c2_start, c2_end = 1326542 - 2000, 1353353 + 2000 #change the location co-ordinates according to the output of find_eps_clusters.py
    cand2 = chrom[c2_start:c2_end]
    cand2.id = "Candidate_2_Flippase"
    cand2.name = "Flippase"
    cand2.description = "Complex capsular/heteropolymer Operon"
    SeqIO.write(cand2, "candidate2.gbk", "genbank")
    print(f"[+] Candidate 2 extracted: {c2_start} to {c2_end}")
    
    print("\nExtraction complete! You can now run clinker.")

if __name__ == "__main__":
    main()
