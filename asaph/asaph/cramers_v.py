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

import numpy as np
from scipy.stats import chisquare

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
    if len(seq1) != len(seq2):
        raise Exception("Two sequences have different lengths")

    seq1 = map(int, seq1)
    seq2 = map(int, seq2)

    set1 = set(seq1)
    set2 = set(seq2)

    if len(set1) == 1 and len(set2) == 1:
        return 1.0
    elif len(set1) == 0 or len(set2) == 0:
        return 0.0

    # fudge factor
    observed = np.ones((max(set1) + 1, max(set2) + 1))
    count1_freq = np.zeros(max(set1) + 1)
    count2_freq = np.zeros(max(set2) + 1)

    for i, j in zip(seq1, seq2):
        observed[i, j] += 1.0
        count1_freq[i] += 1.0
        count2_freq[j] += 1.0

    n_obs = len(seq1)
    count1_freq /= n_obs
    count2_freq /= n_obs

    # fudge factor
    expected = np.ones((max(set1) + 1, max(set2) + 1))
    for i in xrange(max(set1) + 1):
        for j in xrange(max(set2) + 1):
            exp = count1_freq[i] * count2_freq[j] * n_obs
            expected[i, j] += int(np.around(exp))

    length = observed.shape[0] * observed.shape[1]
    observed = observed.reshape(length)
    expected = expected.reshape(length)
    chi2, _ = chisquare(observed, expected)
    phi2 = chi2 / n_obs
    min_dim = float(min(len(set1) - 1, len(set2) - 1))

    v = sqrt(phi2 / min_dim)
    return v
