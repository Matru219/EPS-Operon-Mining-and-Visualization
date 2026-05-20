#!/usr/bin/env python3
import sys
import re

def main():
    if len(sys.argv) != 3:
        print("Usage: python find_eps_clusters.py <eps_builders_GT.txt> <genomic.gff>")
        sys.exit(1)

    gt_file = sys.argv[1]
    gff_file = sys.argv[2]
    cluster_distance_bp = 15000  # Max distance between genes to be considered the same operon

    # 1. Extract protein IDs from the dbCAN GT output
    gt_ids = set()
    with open(gt_file, 'r') as f:
        for line in f:
            if line.strip():
                # dbCAN usually outputs ID in the first column
                gt_id = line.split()[0].replace('.hmm', '') 
                gt_ids.add(gt_id)

    # 2. EPS specific keywords to hunt for in the GFF
    eps_keywords = ["wzx", "wzy", "wzz", "flippase", "polymerase", "chain length", "abc transporter", "tyrosine-protein kinase", "epsA", "epsB", "epsC"]
    
    hits = []
    
    # 3. Parse the GFF3 file
    print("[*] Parsing GFF and locating pathway components...")
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS": continue
            
            contig, start, end, attributes = parts[0], int(parts[3]), int(parts[4]), parts[8]
            
            # Check if this CDS is one of our high-confidence GTs
            is_gt = any(gt in attributes for gt in gt_ids)
            
            # Check if this CDS is a structural/transport EPS component
            is_eps_component = any(kw in attributes.lower() for kw in eps_keywords)
            
            if is_gt or is_eps_component:
                hits.append({
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "desc": attributes,
                    "is_gt": is_gt,
                    "is_eps": is_eps_component
                })

    # 4. Group into physical clusters
    print(f"[*] Scanning for clusters (Genes within {cluster_distance_bp} bp)...")
    clusters = []
    current_cluster = []

    for hit in hits:
        if not current_cluster:
            current_cluster.append(hit)
            continue
        
        last_hit = current_cluster[-1]
        # If on the same contig and within the distance threshold, add to cluster
        if hit["contig"] == last_hit["contig"] and (hit["start"] - last_hit["end"]) <= cluster_distance_bp:
            current_cluster.append(hit)
        else:
            clusters.append(current_cluster)
            current_cluster = [hit]
    
    if current_cluster:
        clusters.append(current_cluster)

    # 5. Filter and print valid operons (Must have at least 1 GT and 1 EPS Component)
    print("\n" + "="*50)
    print("      POTENTIAL EPS BIOSYNTHESIS OPERONS")
    print("="*50)
    
    found_operon = False
    for i, cluster in enumerate(clusters):
        has_gt = any(gene["is_gt"] for gene in cluster)
        has_eps = any(gene["is_eps"] for gene in cluster)
        
        if has_gt and has_eps:
            found_operon = True
            contig = cluster[0]["contig"]
            start = cluster[0]["start"]
            end = cluster[-1]["end"]
            print(f"\n>> Candidate Cluster on {contig} (Positions: {start} - {end})")
            
            for gene in cluster:
                gene_type = "[GT]" if gene["is_gt"] else "[Membrane/Reg]"
                # Clean up the description for readability
                product = re.search(r'product=([^;]+)', gene["desc"])
                prod_name = product.group(1) if product else gene["desc"][:50]
                print(f"   {gene['start']}-{gene['end']} : {gene_type} {prod_name}")

    if not found_operon:
        print("\nNo complete clusters found. The pathway may be split across contigs or highly novel.")

if __name__ == "__main__":
    main()
