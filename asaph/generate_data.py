import argparse
import random

header_left = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"

def snp_generator(n_individuals, n_snps):
    for i in xrange(n_snps):
        yield [(random.randint(0, 1),
                random.randint(0, 1))
               for j in xrange(n_individuals)]

def generate_lines(n_individuals, n_snps):
    header = [header_left]
    for i in xrange(n_individuals):
        header.append(str(i))
    header = "\t".join(header)

    yield header

    for i, snps in enumerate(snp_generator(n_individuals, n_snps)):
        cols = ["1", str(i), ".", "A", "T", "0", "PASS", "AC=30;AF=0.357;AN=84;DP=804;PercentNBaseSolid=0.0000;set=AGC", "GT"]
        for allele1, allele2 in snps:
            cols.append(str(allele1) + "/" + str(allele2))

        yield "\t".join(cols)

def vcf_writer(flname, stream):
    with open(flname, "w") as fl:
        for ln in stream:
            fl.write(ln)
            fl.write("\n")

def pops_writer(flname, n_individuals):
    pops = { "population1" : [],
             "population2" : [] }
    
    for i in xrange(n_individuals):
        pop = random.sample(pops.keys(), 1)[0]
        pops[pop].append(str(i))

    with open(flname, "w") as fl:
        for key, value in pops.iteritems():
            fl.write(key)
            fl.write("\t")
            fl.write(",".join(value))
            fl.write("\n")
        
def parseargs():
    parser = argparse.ArgumentParser()

    parser.add_argument("--output-vcf",
                        type=str,
                        required=True)

    parser.add_argument("--output-populations",
                        type=str,
                        required=True)

    parser.add_argument("--individuals",
                        type=int,
                        required=True)

    parser.add_argument("--snps",
                        type=int,
                        required=True)

    return parser.parse_args()
    

if __name__ == "__main__":
    args = parseargs()

    vcf_writer(args.output_vcf,
               generate_lines(args.individuals, args.snps))

    pops_writer(args.output_populations,
                args.individuals)
    
