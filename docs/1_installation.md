# Installation 


----------


This section will help you set up the necessary software and packages required to run SigProfilerAssignment.

----------


## Prerequisites ##

- [Python][1] version >= 3.9
- Downloaded reference genomes using [SigProfilerMatrixGenerator][2] (only if mutation calling files are used as input)
- Other dependencies and necessary packages are downloaded during the installation

## Installation ##

SigProfilerAssignment can be executed on any Windows/MacOS/Unix system. First follow the [SigProfilerMatrixGenerator][2] guide for installing `Python` and `pip`. Next, follow the download instructions for the latest stable release or the current GitHub version.

### Installation with `pip` ###

Install last `SigProfilerAssignment` PyPI version using `pip`:
```
$ pip install SigProfilerAssignment
```

### Install specific GitHub Release ###

First, download the [zip file][3] or clone the GitHub repository by:
```
$ git clone https://github.com/SigProfilerSuite/SigProfilerAssignment.git
```

Next, enter the downloaded directory and install the package by unzipping the contents of SigProfilerAssignment-master or the zip file of a corresponding branch:
```
$ cd SigProfilerAssignment
$ pip install .
```

## Download Reference Genome ##

In case you want to use SigProfilerAssignment with mutation calling files as input, you first need to download the appropriate reference genome. Current reference genomes supported include GRCh37, GRCh38, mm9, mm10, rn6 and rn7. To install the reference genome/s, you need to use [SigProfilerMatrixGenerator][2].

The last PyPI [SigProfilerMatrixGenerator][2] version is installed with SigProfilerAssignment by default. You can also install a specific version following the instructions in [SigProfilerMatrixGenerator Wiki][2].

Once [SigProfilerMatrixGenerator][2], install your desired reference genome from the command line/terminal as follows. 

### Installation from command line ###

```
$ SigProfilerMatrixGenerator install GRCh37
```
### Installation from Python terminal ###

``` python
$ python
>> from SigProfilerMatrixGenerator import install as genInstall
>> genInstall.install('GRCh37', rsync=False, bash=True)
```

In case you prefer to install a reference genome that you have saved locally, you can do the following:
``` python
$ python
>> from SigProfilerMatrixGenerator import install as genInstall
>> genInstall.install('GRCh37', offline_files_path='path/to/directory/containing/GRCh37.tar.gz')
```

  [1]: https://www.python.org/downloads
  [2]: https://sigprofilersuite.github.io/SigProfilerMatrixGenerator/Installation-Python.html
  [3]: https://github.com/SigProfilerSuite/SigProfilerAssignment/releases

