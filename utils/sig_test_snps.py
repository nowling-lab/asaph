import argparse

import numpy as np

def read_data(flname):
    with open(flname) as fl:
        data = []
        for ln in fl:
            cols = ln.split()
            try:
                v = float(cols[2])
                if not np.isnan(v):
                    data.append((v, ln))
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

    args.add_argument("--significance",
                      type=float)

    return args.parse_args()

if __name__ == "__main__":
    args = parseargs()

    data = read_data(args.input)

    if args.significance is not None:
        data = filter(lambda p: p[0] <= args.significance, data)

    data.sort(key=lambda p: p[0])
    
    with open(args.output, "w") as fl:
        for _, ln in data:
            fl.write(ln)
