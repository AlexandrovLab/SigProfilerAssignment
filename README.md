[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://osf.io/t6j7u/wiki/home/) 
[![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![Build Status](https://api.travis-ci.com/AlexandrovLab/SigProfilerAssignment.svg?branch=master)](https://app.travis-ci.com/AlexandrovLab/SigProfilerAssignment)

# SigProfilerAssignment

SigProfilerAssignment is a new mutational attribution and decomposition tool that performs the following functions:
-   Attributing a known set of mutational signatures to an individual sample or multiple samples.
-   Decomposing de novo signatures to COSMIC signature database.
-   Attributing COSMIC database or a custom signature database to given samples.

The tool identifies the activity of each signature in the sample and assigns the probability for each signature to cause a specific mutation type in the sample. The tool makes use of SigProfilerMatrixGenerator, SigProfilerExtractor and SigProfilerPlotting.


## Installs
for installing dependencies in new conda environment

```
$ pip install -r requirements.txt
```

Installing this package : git clone this repo or download the zip file.
Unzip the contents of SigProfilerExtractor-master.zip or the zip file of a corresponding branch.

```bash
$ cd SigProfilerAssignment-master
$ pip install .
```


Decomposes the De Novo Signatures into COSMIC Signatures and assigns COSMIC signatures into samples

<!-- 
```python
spa_analyze(  samples,  output, signatures=None, signature_database=None,decompose_fit= True,denovo_refit=True,cosmic_fit=True, nnls_add_penalty=0.05, 
              nnls_remove_penalty=0.01, initial_remove_penalty=0.05, de_novo_fit_penalty=0.02, 
              genome_build="GRCh37",  make_decomposition_plots=True, collapse_to_SBS96=True,connected_sigs=True, verbose=False): 
```  -->
### Decompose Fit
```python
from SigProfilerAssignment import Analyzer as Analyze
Analyze.decompose_fit(samples,  output, signatures=None, signature_database=None,genome_build="GRCh37",  make_decomposition_plots=True, collapse_to_SBS96=True,connected_sigs=True, verbose=False)
```
### *De Novo* Fit
```python
from SigProfilerAssignment import Analyzer as Analyze
Analyze.denovo_fit(samples,  output, signatures=None, signature_database=None,genome_build="GRCh37",  make_decomposition_plots=True, collapse_to_SBS96=True,connected_sigs=True, verbose=False)
```
### Cosmic Fit
```python
from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(samples,  output, signatures=None, signature_database=None,genome_build="GRCh37",  make_decomposition_plots=True, collapse_to_SBS96=True,connected_sigs=True, verbose=False)
```
## Parameters
| Parameter | Variable Type | Parameter Description |
| --------------------- | -------- |-------- |
| **signatures** | String | Path to a  tab delimited file that contains the signaure table where the rows are mutation types and colunms are signature IDs. |
| **activities** | String | Path to a tab delimilted file that contains the activity table where the rows are sample IDs and colunms are signature IDs. |
| **samples** | String | Path to a tab delimilted file that contains the activity table where the rows are mutation types and colunms are sample IDs. |
| **output** | String | Path to the output folder. |
| **genome_build** | String | The genome type. Example: "GRCh37", "GRCh38", "mm9", "mm10". The default value is "GRCh37" |
| **verbose** | Boolean | Prints statements. Default value is False.  |
        

#### SPA analysis Example


```python
#import modules
import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze

#set directories and paths to signatures and samples
dir_inp = spa.__path__[0]+'/data/Examples/'
signatures = dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt"
activities=dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Activities/SBS96_S3_NMF_Activities.txt"
samples=dir_inp+"Input_scenario_8/Samples.txt"
output="output_example/"
sigs= "COSMIC_v3_SBS_GRCh37_noSBS84-85.txt" #Custom Signature Database

#Analysis of SP Assignment 
Analyze.cosmic_fit( samples, output, signatures=None,signature_database=sigs,genome_build="GRCh37", verbose=False)

```
## <a name="copyright"></a> Copyright
This software and its documentation are copyright 2022 as a part of the SigProfiler project. The SigProfilerAssignment framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## <a name="contact"></a> Contact Information
Please address any queries or bug reports to Raviteja Vangara at rvangara@health.ucsd.edu
