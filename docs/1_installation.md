# Installation 


----------


This section will help you set up the necessary software and packages required to run SigProfilerAssignment.

@[toc](Sections)

----------


## Prerequisites ##
- [Python][1] version >= 3.4.0
- Downloaded reference genomes using [SigProfilerMatrixGenerator][2] (only if mutation calling files are used as input)
- Other dependencies and necessary packages are downloaded during the installation

## Download Reference Genome ##
In case you want to use SigProfilerAssignment with mutation calling files as input, you first need to download the appropriate reference genome. Current reference genomes supported include GRCh37, GRCh38, mm9, mm10, and rn6. To install the reference genome/s, you need to use [SigProfilerMatrixGenerator][2].

First, install the python package using pip:
```
$ pip install SigProfilerMatrixGenerator
```
Install your desired reference genome from the command line/terminal as follows:
```
$ python
>> from SigProfilerMatrixGenerator import install as genInstall
>> genInstall.install('GRCh37', rsync=False, bash=True)
```
In case you prefer to install a reference genome that you have saved locally, you can do the following:
```
$ python
>> from SigProfilerMatrixGenerator import install as genInstall
>> genInstall.install('GRCh37', offline_files_path='path/to/directory/containing/GRCh37.tar.gz')
```

## Mac/Unix ##
First, follow the [SigProfilerMatrixGenerator][3] Mac/Unix guide for installing `Python` and `pip`. Next, follow the download instructions for the latest stable release or the current GitHub version.

### Mac/Unix Stable Release ###
Install `SigProfilerAssignment` using `pip`:
```
$ pip install SigProfilerAssignment
```

### Mac/Unix GitHub Release ###
First, download the [zip file][4] or clone the GitHub repository by:
```
$ git clone https://github.com/SigProfilerSuite/SigProfilerAssignment.git
```

Next, enter the downloaded directory and install the package by unzipping the contents of SigProfilerAssignment-master or the zip file of a corresponding branch:
```
$ cd SigProfilerAssignment
$ pip install .
```
After SigProfilerAssignment successfully installs, the environment is set up and ready to use.

## Windows ##
First, follow the [SigProfilerMatrixGenerator][3] Windows guide for installing `Python` and `pip`. Next, follow the download instructions for the latest stable release or the current GitHub version. 


### Windows Stable Release ###
Install `SigProfilerAssignment` using `pip`:
```
$ pip install SigProfilerAssignment
```

### Windows GitHub Release ###
First, download the [zip file][4] or clone the GitHub repository by:
```
$ git clone https://github.com/SigProfilerSuite/SigProfilerAssignment.git
```

Next, enter the downloaded directory and install the package by unzipping the contents of SigProfilerAssignment-master or the zip file of a corresponding branch:
```
$ cd SigProfilerAssignment
$ pip install .
```
After SigProfilerAssignment successfully installs, the environment is set up and ready to use.


  [1]: https://www.python.org/downloads
  [2]: https://osf.io/s93d5/wiki/1.%20Installation%20-%20Python/
  [3]: https://osf.io/s93d5/wiki/1.%20Installation%20-%20Python/
  [4]: https://github.com/SigProfilerSuite/SigProfilerAssignment/releases

