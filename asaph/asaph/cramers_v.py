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

from collections import defaultdict
from itertools import product
from math import sqrt

def count_occurrences(seq):
    counts = defaultdict(int)
    for item in seq:
        counts[item] += 1

    return counts
        

def cramers_v(seq1, seq2):
    """
    Computes the assocation between two nominal variables using Cramer's V.

    When used to analyze variants, one sequence should be the population labels.
    The other should encode the variants as categories (e.g., "aa", "tt", "at").

    >>> cramers_v([0, 1, 0, 1], [1, 2, 1, 2])
    1.0
    """
    set1 = set(seq1)
    set2 = set(seq2)

    if len(set1) == 1 and len(set2) == 1:
        return 1.0
    elif len(set1) == 0 or len(set2) == 0:
        return 0.0
    elif len(set1) == 1 or len(set2) == 1:
        return 0.0

    pairCounts = count_occurrences(zip(seq1, seq2))
    counts1 = count_occurrences(seq1)
    counts2 = count_occurrences(seq2)

    n_obs = float(len(seq1))

    chi2 = 0.0
    for value1, value2 in product(set1, set2):
        nij = pairCounts[(value1, value2)]
        ni = counts1[value1]
        nj = counts2[value2]

        b = ni * nj / n_obs
        c = (nij - b) * (nij - b) / b

        chi2 += c

    min_dim = float(min(len(set1) - 1, len(set2) - 1))

    v = sqrt(chi2 / n_obs / min_dim)

    return v
