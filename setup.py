"""
Copyright 2020 Ronald J. Nowling

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

from setuptools import setup

setup(name="asaph",
      version=2.0,
      description="SNP analysis",
      author="Ronald J. Nowling",
      author_email="rnowling@gmail.com",
      license="Apache License, Version 2.0",
      zip_safe=False,
      packages=["asaph"],
      python_requires=">=3.6",
      install_requires = ["numpy>=0.19.1", "scipy>=0.19.1", "matplotlib", "seaborn", "sklearn", "joblib", "pandas", "mmh3"],
      scripts=["bin/asaph_pca", "bin/asaph_query",
               "bin/asaph_localize", "bin/asaph_genotype", "bin/asaph_generate_data"])
