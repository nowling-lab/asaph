import argparse

import numpy as np

def read_data(flname, chromosome_id):
    with open(flname) as fl:
        data = []
        for ln in fl:
            cols = ln.split()
            try:
                molecule = cols[0]
                if molecule == chromosome_id:
                    data.append(cols)
            except:
                pass
    return data
            
def parseargs():
    args = argparse.ArgumentParser()

    args.add_argument("--input",
                      type=str,
                      required=True)

    args.add_argument("--output",
                      type=str,
                      required=True)

    args.add_argument("--select-id",
                      type=str,
                      required=True)

    args.add_argument("--output-id",
                      type=str,
                      required=True)

    return args.parse_args()

if __name__ == "__main__":
    args = parseargs()

    data = read_data(args.input,
                     args.select_id)

    with open(args.output, "w") as fl:
        for rec in data:
            rec[0] = args.output_id
            fl.write("\t".join(rec))
            fl.write("\n")
