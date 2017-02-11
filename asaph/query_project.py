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
import os

from asaph.newioutils import read_features

def query_project(workdir):
    if not os.path.exists(workdir):
        raise Exception, "workdir '%s' does not exist." % workdir

    features = read_features(workdir)

    print "Number of positions:", len(features.snp_feature_map)
    print "Number of samples:", len(features.sample_labels)
    print "Number of features:", features.feature_matrix.shape[1]

def parseargs():
    parser = argparse.ArgumentParser()

    parser.add_argument("--workdir",
                        type=str,
                        required=True)

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    query_project(args.workdir)
