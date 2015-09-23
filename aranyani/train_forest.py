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
from sklearn.ensemble import RandomForestClassifier
from ioutils import *

if __name__ == "__main__":
    basename = sys.argv[1]
    n_trees = int(sys.argv[2])

    labels, selected_indices, selected_individuals = select_groups(basename, True)
    feature_labels = read_features(basename, True)
    feature_matrix = open_feature_matrix(basename, True)
    
    rf1 = RandomForestClassifier(n_estimators=n_trees)
    rf1.fit(feature_matrix, labels)

    rf2 = RandomForestClassifier(n_estimators=n_trees)
    rf2.fit(feature_matrix, labels)

    serialize_models(basename, rf1, rf2)

    


