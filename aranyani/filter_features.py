"""
Copyright 2015 Ronald J. Nowling

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

import sys
from ioutils import *

def filter_features(feature_labels, feature_matrix):
    n_samples, n_features = feature_matrix.shape

    kept_features = []
    kept_columns = []
    for i in xrange(n_features):
        nzs = np.count_nonzero(feature_matrix[:, i])
        # remove columns that are all 0 or all 1
        if nzs != 0 and nzs != n_samples:
            kept_features.append(feature_labels[i])
            kept_columns.append(i)

    print "Removed %s out of %s features" % (n_features - len(kept_columns),
                                             n_features)
    print "%s features remaining" % len(kept_columns)

    return kept_features, feature_matrix[:, kept_columns]

if __name__ == "__main__":
    basename = sys.argv[1]

    labels, selected_indices, selected_individuals = select_groups(basename, False)
    feature_labels = read_features(basename)
    feature_matrix = open_feature_matrix(basename)
    
    feature_matrix = feature_matrix[selected_indices, :]

    feature_labels, feature_matrix = filter_features(feature_labels, feature_matrix)
    n_individuals, n_features = feature_matrix.shape

    filtered_feature_matrix = create_feature_matrix(basename, n_individuals, \
                                                    n_features, filtered=True)
    filtered_feature_matrix[:, :] = feature_matrix
    filtered_feature_matrix.flush()
    del filtered_feature_matrix

    write_features(basename, feature_labels, filtered=True)
    write_individuals(basename, selected_individuals, filtered=True)
    


