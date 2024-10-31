# Installing Asaph

Asaph requires Python 3.6 and above and the [numpy and scipy](http://www.numpy.org/), [matplotlib](http://matplotlib.org/), [scikit-learn](http://scikit-learn.org/stable/), [seaborn](https://seaborn.pydata.org/index.html), and [pandas](https://pandas.pydata.org/) libraries.

## Getting Asaph
Download Asaph by cloning the git repository:

```bash
$ git clone https://github.com/rnowling/asaph.git
```

## Option 1. Install in Virtual Environment
If you would prefer to use a virtual environment instead of installing system wide, you can:

```bash
$ cd asaph
$ python3 -m venv venv
$ source venv/bin/activate
$ pip3 install setuptools
$ python3 setup.py install
```

Asaph will install the scripts located under `venv/bin` system wide.  The binaries will be on your path when the virtual environment is activated.

## Option 2. Install System-Wide
You can install Asaph using Python's setup tools:

```bash
$ cd asaph
$ sudo pip3 install setuptools
$ sudo python3 setup.py install
```

Asaph will install the scripts located under `bin/` system wide.

## Debugging Your installation

If, for some reason, the installation of Asaph failed, install these libraries and try again:

```bash
$ pip3 install -U numpy scipy scikit-learn matplotlib seaborn pandas
```


## What Next?
Now that Asaph is installed, check out the [next tutorial](pca.md) on how to prepare, import, and perform PCA on SNP data.
