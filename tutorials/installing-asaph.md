# Installing Asaph

Asaph requires Python 3.7 and above and the [numpy and scipy](http://www.numpy.org/), [matplotlib](http://matplotlib.org/), [scikit-learn](http://scikit-learn.org/stable/), [seaborn](https://seaborn.pydata.org/index.html), and [pandas](https://pandas.pydata.org/) libraries.

## Getting Asaph
Download Asaph by cloning the git repository:

```bash
$ git clone https://github.com/rnowling/asaph.git
```

## Recommended: Install in Virtual Environment
Using a virtual environment is the recommended approach to avoid conflicts with system packages:

```bash
$ cd asaph
$ python3 -m venv venv
$ source venv/bin/activate
$ pip install -e .
```

The `-e` flag installs Asaph in "editable" mode, which is useful for development. The scripts will be available on your path when the virtual environment is activated.

## Alternative: Install System-Wide
You can install Asaph system-wide using pip:

```bash
$ cd asaph
$ pip install .
```

For system-wide installation, you may need to use `sudo` depending on your system configuration.

## Legacy Installation (Deprecated)
The old setup.py approach still works but is deprecated:

```bash
$ cd asaph
$ python3 setup.py install
```

This method will show deprecation warnings and should be avoided in favor of the modern pip-based installation.

## Debugging Your Installation

If the installation fails, it's usually due to missing dependencies. The modern pip installation should handle dependencies automatically, but you can manually install them if needed:

```bash
$ pip install -U numpy scipy scikit-learn matplotlib seaborn pandas mmh3 joblib
```

If you're still having issues, try updating pip itself:

```bash
$ pip install --upgrade pip
```


## What Next?
Now that Asaph is installed, check out the [next tutorial](pca.md) on how to prepare, import, and perform PCA on SNP data.
