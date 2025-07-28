"""
This module defines data structures and named tuples for representing project summaries,
including feature counts, sample information, and dimensionality reduction parameters
for population genetics analysis workflows.

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

from collections import namedtuple

ProjectSummary = namedtuple("ProjectSummary",
                            ["n_features",
                             "n_samples",
                             "feature_type",
                             "sampling_method",
                             "sample_names",
                             "explained_variance_ratios"])
