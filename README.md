[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://osf.io/mz79v/wiki/home/) 
[![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![CI](https://github.com/SigProfilerSuite/SigProfilerAssignment/actions/workflows/ci.yml/badge.svg)](https://github.com/SigProfilerSuite/SigProfilerAssignment/actions/workflows/ci.yml)

<img src="SigProfilerAssignment/figures/SigProfilerAssignment.png" alt="drawing" width="1000"/>

# SigProfilerAssignment
SigProfilerAssignment enables assignment of previously known mutational signatures to individual samples and individual somatic mutations. The tool refits different types of reference mutational signatures, including [COSMIC signatures](https://cancer.sanger.ac.uk/signatures/), as well as custom signature databases. SigProfilerAssignment makes use of [SigProfilerMatrixGenerator](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator) and [SigProfilerPlotting](https://github.com/SigProfilerSuite/SigProfilerPlotting), seamlessly integrating with other in [SigProfilerSuite](https://github.com/SigProfilerSuite).

## Documentation
Detailed documentation can be found at https://sigprofilersuite.github.io/SigProfilerAssignment.

## Quick Start Guide
###  Installation

Install the current stable PyPi version of SigProfilerAssignment:
```
$ pip install SigProfilerAssignment
```

If mutation calling files (MAF, VCF, or simple text files) are used as input, please install your desired reference genome as follows (available reference genomes are: GRCh37, GRCh38, mm9, mm10, rn6, and rn7):
```python
$ python
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh37')
```


### Running

Assignment of known mutational signatures to individual samples is performed using the `cosmic_fit` function. Input samples are provided using the `samples` parameter in the form of mutation calling files (VCFs, MAFs, or simple text files), segmentation files, or mutational matrices. COSMIC mutational signatures v3.5 are used as the default reference signatures, although previous COSMIC versions and custom signature databases are also supported using the `cosmic_version` and `signature_database` parameters. Results will be found in the folder specified in the `output` parameter.

```python
from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(samples, output, input_type="matrix", context_type="96", genome_build="GRCh37")
```

You can also run SigProfilerAssignment `cosmic_fit` function from command line: 

``` bash
$ SigProfilerAssignment cosmic_fit samples output --input_type "matrix" --context_type "96" --genome_build "GRCh37"

```

## Reference

Díaz-Gay M, Vangara R, Barnes M, Wang X, Islam SMA, Vermes I, Duke S, Narasimman NB, Yang T, Jiang Z, Moody S, Senkin S, Brennan P, Stratton MR, Alexandrov LB. Assigning mutational signatures to individual samples and individual somatic mutations with SigProfilerAssignment. *Bioinformatics*. 2023;39(12):btad756. [https://doi.org/10.1093/bioinformatics/btad756](https://doi.org/10.1093/bioinformatics/btad756)

## Contact

For questions, support requests, or bug reports, please contact the SigProfilerSuite team via GitHub [issues](https://github.com/SigProfilerSuite/SigProfilerAssignment/issues) or by email at [contact@sigprofilersuite.org](mailto:contact@sigprofilersuite.org).

