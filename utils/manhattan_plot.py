"""
Copyright 2018 Ronald J. Nowling

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import argparse
from collections import defaultdict
import gzip

import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

import numpy as np

import seaborn as sns

SAMPLE_START_COLUMN = 9

def read_snp_probs(flname):
    with open(flname) as fl:
        for ln in fl:
            cols = ln.strip().split()

            pos = int(cols[1])
            pvalue = float(cols[2])
            
            yield pos, pvalue
            

def parseargs():
    parser = argparse.ArgumentParser()
            
    parser.add_argument("--input-tsv",
                        required=True,
                        type=str)

    parser.add_argument("--plot-fl",
                        required=True,
                        type=str)

    parser.add_argument("--y-limit",
                        type=float)

    parser.add_argument("--sig-line",
                        type=float)
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    xs = []
    ys = []
    for pos, pvalue in read_snp_probs(args.input_tsv):
        xs.append(pos)
        ys.append(-np.log10(pvalue))

    plt.scatter(xs,
                ys,
                marker=".",
                color="c")

    if args.sig_line:
        plt.plot([0, max(xs)],
                 [-np.log10(args.sig_line),
                  -np.log10(args.sig_line)],
                 "r-")

    plt.xlabel("Position (bp)", fontsize=16)
    plt.ylabel("p-value (-log10)", fontsize=16)

    if args.y_limit:
        plt.ylim([0, args.y_limit])
    else:
        plt.ylim([0.0, max(ys)])

    plt.xlim([0.0, max(xs)])

    plt.savefig(args.plot_fl,
                DPI=300)
                

        


