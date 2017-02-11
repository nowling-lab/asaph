"""
Copyright 2017 Ronald J. Nowling

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
import csv
import os

import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

def read_rankings(flname):
    rankings = []
    with open(flname) as fl:
        reader = csv.reader(fl, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        for rec in reader:
            rankings.append(tuple(rec[:2]))

    return rankings

def analyze_rankings(universe_size, ranking_flname_1, ranking_flname_2):
    rankings_1 = read_rankings(ranking_flname_1)
    rankings_2 = read_rankings(ranking_flname_2)

    thresholds = [0.01, 0.05, 0.1, 0.25, 0.5]
    min_len = min(len(rankings_1), len(rankings_2))
    if universe_size == None:
        universe_size = min_len
    for t in thresholds:
        cutoff = min(int(t * universe_size) + 1, min_len)
        intersection = set(rankings_1[:cutoff]).intersection(rankings_2[:cutoff])
        percentage = 100.0 * float(len(intersection)) / float(cutoff)
        print t, percentage

def parseargs():
    parser = argparse.ArgumentParser()

    parser.add_argument("--universe-size",
                        type=int)

    parser.add_argument("--rankings-fl-1",
                        type=str,
                        required=True)

    parser.add_argument("--rankings-fl-2",
                        type=str,
                        required=True)

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    analyze_rankings(universe_size,
                     args.rankings_fl_1,
                     args.rankings_fl_2)

