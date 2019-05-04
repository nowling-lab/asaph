import argparse
import random

import numpy as np

format_header = "##fileformat=VCFv4.1"
header_left = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"

def generate_vcf(flname, positions, diploids):
    n_sites = len(positions)
    n_individuals = len(diploids)

    positions = list(positions)
    positions.sort()

    header = [header_left]
    for i in range(n_individuals):
        header.append("indiv_" + str(i+1))
    
    header = "\t".join(header)
            
    yield format_header
    yield header

    for site_idx, pos in enumerate(positions):
        cols = ["1", str(pos), ".", "A", "T", "0", "PASS", "AC=30;AF=0.357;AN=84;DP=804;PercentNBaseSolid=0.0000;set=AGC", "GT"]

        for indiv_idx in range(n_individuals):
            gt = diploids[indiv_idx][1].get(pos, 0)
            
            if gt == 0:
                cols.append("0/0")
            elif gt == 2:
                cols.append("1/1")
            else:
                cols.append("0/1")

        yield "\t".join(cols)

def write_vcf(flname, positions, diploids):
    stream = generate_vcf(flname, positions, diploids)
    with open(flname, "w") as fl:
        for ln in stream:
            fl.write(ln)
            fl.write("\n")

def read_snps(flname):
    chromosomes = []
    recording = False
    all_snp_positions = set()
    with open(flname) as fl:
        for ln in fl:
            if not recording and ln.startswith("<DATA>"):
                recording = True
                start = len("<DATA>")
                snp_positions = set(map(int, ln[start:].strip().split()))
                # invertFREGENE uses 0 to indicate the end of a chromosome
                # also seems to use an EOL
                snp_positions.remove(0)
                chromosomes.append((len(chromosomes), snp_positions))
                all_snp_positions.update(snp_positions)
            elif recording and ln.startswith("</DATA>"):
                recording = False
            elif recording:
                snp_positions = set(map(int, ln[5:].strip().split()))
                snp_positions.remove(0)
                chromosomes.append((len(chromosomes), snp_positions))
                all_snp_positions.update(snp_positions)

    return all_snp_positions, chromosomes

def read_karyotypes(flname):
    karyotypes = dict()
    with open(flname) as fl:
        for ln in fl:
            cols = ln.strip().split()
            karyotypes[int(cols[0])] = int(cols[1])
    
    return karyotypes

def form_diploids(chromosomes, karyotypes):
    random.shuffle(chromosomes)

    diploids = []
    diploid_karyotypes = []
    for i in range(0, len(chromosomes), 2):
        idx1, chrom1 = chromosomes[i]
        idx2, chrom2 = chromosomes[i+1]
        k = karyotypes[idx1] + karyotypes[idx2]

        genotypes = dict()
        for pos in chrom1:
            genotypes[pos] = 1

        for pos in chrom2:
            genotypes[pos] = genotypes.get(pos, 0) + 1
        
        diploids.append((k, genotypes))

    return diploids

def write_pops(basename, diploids):
    groups = { 0 : [],
               1 : [],
               2 : [] }

    for i, (kt, _) in enumerate(diploids):
        groups[kt].append("indiv_" + str(i+1))

    print(groups)
        
    with open(basename + ".pops", "w") as fl:
        for group_id, pop in groups.items():
            if group_id == 0:
                group_name = "Homo. Std."
            elif group_id == 1:
                group_name = "Hetero."
            else:
                group_name = "Homo. Inv."
            fl.write(group_name)

            for ident in pop:
                fl.write(",")
                fl.write(ident)
            fl.write("\n")

    for group_id, pop in groups.items():
        if group_id == 0:
            group_name = "homo_std"
        elif group_id == 1:
            group_name = "hetero"
        else:
            group_name = "homo_inv"

        with open(basename + "_" + group_name + ".ids", "w") as fl:
            for ident in pop:
                fl.write(ident)
                fl.write("\n")
    

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--sim-input",
                        type=str,
                        required=True)

    parser.add_argument("--karyotype-input",
                        type=str,
                        required=True)
    
    parser.add_argument("--output-base",
                        type=str,
                        required=True)
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    all_snp_positions, chromosomes = read_snps(args.sim_input)

    karyotypes = read_karyotypes(args.karyotype_input)

    print(len(all_snp_positions), len(chromosomes), len(karyotypes))

    diploids = form_diploids(chromosomes, karyotypes)
    print(len(diploids))
    for i in range(5):
        print(diploids[i][0])

    print(all_snp_positions)
    print()
    print(diploids[0])


    write_vcf(args.output_base + ".vcf",
              all_snp_positions,
              diploids)

    write_pops(args.output_base,
               diploids)
