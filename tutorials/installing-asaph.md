# Installing Asaph

The [Asaph](https://github.com/rnowling/asaph) software implements methods for analyzing large inversions using SNP data.  Asaph and its methods have been described in our 2018 ACM-BCB paper [Detecting Chromosomal Inversions from Dense SNPs by Combining PCA and Association Tests](https://dl.acm.org/citation.cfm?id=3233571) and 2020 PLOS ONE paper [Detecting inversions with PCA in the presence of population structure](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0240429). In this tutorial, we describe how to install Asaph.

## Getting Asaph
Download Asaph by cloning the git repository:

```bash
$ git clone https://github.com/rnowling/asaph.git
```

## Installing Asaph
You can install Asaph using Python's setup tools:

```bash
$ cd asaph
$ python3 setup.py install
```

Asaph requires Python 3 and the [numpy and scipy](http://www.numpy.org/), [matplotlib](http://matplotlib.org/), [scikit-learn](http://scikit-learn.org/stable/), and [seaborn](https://seaborn.pydata.org/index.html) libraries.  If, for some reason, the installation of Asaph failed, you can manually update these libraries and try again:

```bash
$ pip3 install -U numpy scipy sklearn matplotlib seaborn
```

If you would prefer to use a virtual environment instead of installing system wide, you can:

```bash
$ cd asaph
$ python3 -m venv venv
$ source venv/bin/activate
# python3 setup.py install
```

Asaph will install the scripts located under `bin/` system wide.

## What Next?
Now that Asaph is installed, check out some of our other [tutorials](README.md).
