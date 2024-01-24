[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://osf.io/mz79v/wiki/home/) 
[![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![Build Status](https://api.travis-ci.com/AlexandrovLab/SigProfilerAssignment.svg?branch=main)](https://app.travis-ci.com/AlexandrovLab/SigProfilerAssignment)

<img src="SigProfilerAssignment/figures/SigProfilerAssignment.png" alt="drawing" width="1000"/>

# SigProfilerAssignment
SigProfilerAssignment enables assignment of previously known mutational signatures to individual samples and individual somatic mutations. The tool refits different types of reference mutational signatures, including [COSMIC signatures](https://cancer.sanger.ac.uk/signatures/), as well as custom signature databases. Refitting of known mutational signatures is a numerical optimization approach that not only identifies the set of operative mutational signatures in a particular sample, but also quantifies the number of mutations assigned to each signature found in that sample. SigProfilerAssignment makes use of [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator) and [SigProfilerPlotting](https://github.com/AlexandrovLab/SigProfilerPlotting), seamlessly integrating with other [SigProfiler tools](https://cancer.sanger.ac.uk/signatures/tools/).

For users that prefer working in an R environment, a wrapper package is provided and can be found and installed from: https://github.com/AlexandrovLab/SigProfilerAssignmentR. Detailed documentation can be found at: https://osf.io/mz79v/wiki/home/.


## Table of contents
- [Installation](#installation)
- [Running](#running)
  - [Main Parameters](#parameters)
  - [Signature Subgroups](#subgroups)
- [Examples](#examples)
- [_De novo_ extraction of mutational signatures downstream analysis](#denovo)
- [Citation](#citation)
- [Copyright](#copyright)
- [Contact Information](#contact)

## <a name="installation"></a> Installation

Install the current stable PyPi version of SigProfilerAssignment:
```
$ pip install SigProfilerAssignment
```

If mutation calling files (MAF, VCF, or simple text files) are used as input, please install your desired reference genome as follows (available reference genomes are: GRCh37, GRCh38, mm9, mm10, and rn6):
```python
$ python
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh37')
```
## <a name="running"></a> Running

Assignment of known mutational signatures to individual samples is performed using the `cosmic_fit` function. Input samples are provided using the `samples` parameter in the form of mutation calling files (VCFs, MAFs, or simple text files), segmentation files or mutational matrices. COSMIC mutational signatures v3.4 are used as the default reference signatures, although previous COSMIC versions and custom signature databases are also supported using the `cosmic_version` and `signature_database` parameters. Results will be found in the folder specified in the `output` parameter.

```python
from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(samples, output, input_type="matrix", context_type="96",
                   collapse_to_SBS96=True, cosmic_version=3.4, exome=False,
                   genome_build="GRCh37", signature_database=None,
                   exclude_signature_subgroups=None, export_probabilities=False,
                   export_probabilities_per_mutation=False, make_plots=False,
                   sample_reconstruction_plots=False, verbose=False)
```

<!-- (nnls_add_penalty=0.05, nnls_remove_penalty=0.01, initial_remove_penalty=0.05, connected_sigs=True) -->

### <a name="parameters"></a> Main Parameters

| Parameter | Variable Type | Parameter Description |
| ------ | ----------- | ----------- |
| samples | String | Path to the input somatic mutations file (if using segmentation file/mutational matrix) or input folder (mutation calling file/s). |
| output | String | Path to the output folder. |
| input_type | String | Three accepted input types:<ul><li> "vcf": if using mutation calling file/s (VCF, MAF, simple text file) as input</li><li>"seg:TYPE": if using a segmentation file as input. Please check the required format at https://github.com/AlexandrovLab/SigProfilerMatrixGenerator#copy-number-matrix-generation. The accepted callers for TYPE are the following {"ASCAT", "ASCAT_NGS", "SEQUENZA", "ABSOLUTE", "BATTENBERG", "FACETS", "PURPLE", "TCGA"}. For example:"seg:BATTENBERG"</li><li>"matrix": if using a mutational matrix as input</li></ul>The default value is "matrix". |
| context_type | String | Required context type if `input_type` is "vcf". `context_type` takes which context type of the input data is considered for assignment. Valid options include "96", "288", "1536", "DINUC", and "ID". The default value is "96". |
| cosmic_version | Float | Defines the version of the COSMIC reference signatures. Takes a positive float among 1, 2, 3, 3.1, 3.2, 3.3, and 3.4. The default value is 3.4. |
| exome | Boolean | Defines if the exome renormalized COSMIC signatures will be used. The default value is False. |
| genome_build | String | The reference genome build, used for select the appropriate version of the COSMIC reference signatures, as well as processing the mutation calling file/s. Supported genomes include "GRCh37", "GRCh38", "mm9", "mm10" and "rn6". The default value is "GRCh37". If the selected genome is not in the supported list, the default genome will be used. |
| signature_database | String | Path to the input set of known mutational signatures (only in case that COSMIC reference signatures are not used), a tab delimited file that contains the signature matrix where the rows are mutation types and columns are signature IDs. |
| exclude_signature_subgroups | List | Removes the signatures corresponding to specific subtypes to improve refitting (only available when using default COSMIC reference signatures). The usage is explained below. The default value is None, which corresponds to use all COSMIC signatures. |
| export_probabilities | Boolean | Defines if the probability matrix per mutational context for all samples is created. The default value is True. |
| export_probabilities_per_mutation | Boolean | Defines if the probability matrices per mutation for all samples are created. Only available when `input_type` is "vcf". The default value is False. |
| make_plots | Boolean | Toggle on and off for making and saving plots. The default value is True. |
| sample_reconstruction_plots | String | Select the output format for sample reconstruction plots. Valid inputs are {'pdf', 'png', 'both', None}. The default value is None. |
| verbose | Boolean | Prints detailed statements. The default value is False. |



### <a name="subgroups"></a> Signature Subgroups

When using COSMIC reference signatures, some subgroups of signatures can be removed to improve the refitting analysis. To use this feature, the `exclude_signature_subgroups` parameter should be added, following the sintax below:

```python
exclude_signature_subgroups = ['MMR_deficiency_signatures',
                               'POL_deficiency_signatures',
                               'HR_deficiency_signatures' ,
                               'BER_deficiency_signatures',
                               'Chemotherapy_signatures',
                               'Immunosuppressants_signatures'
                               'Treatment_signatures'
                               'APOBEC_signatures',
                               'Tobacco_signatures',
                               'UV_signatures',
                               'AA_signatures',
                               'Colibactin_signatures',
                               'Artifact_signatures',
                               'Lymphoid_signatures']
```

The full list of signature subgroups is included in the following table:

|Signature subgroup |           SBS signatures excluded | DBS signatures excluded | ID signatures excluded |
| ----------- | ----------- | ----------- | ----------- |
|MMR_deficiency_signatures|     6, 14, 15, 20, 21, 26, 44|      7, 10|  7|
|POL_deficiency_signatures|     10a, 10b, 10c, 10d, 28|         3|      -|
|HR_deficiency_signatures|      3|                              13|      6|
|BER_deficiency_signatures|     30, 36|                         -|      -|
|Chemotherapy_signatures|       11, 25, 31, 35, 86, 87, 90, 99|     5|      -|
|Immunosuppressants_signatures| 32|                             -|      -|
|Treatment_signatures|          11, 25, 31, 32, 35, 86, 87, 90, 99| 5|      -|
|APOBEC_signatures|             2, 13|                          -|      -|
|Tobacco_signatures |           4, 29, 92|                      2|      3|
|UV_signatures|                 7a, 7b, 7c, 7d, 38|             1|      13|
|AA_signatures|                 22a, 22b|                             20|      23|
|Colibactin_signatures|         88|                             -|      18|
|Artifact_signatures|           27, 43, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 95|14|-|
|Lymphoid_signatures|           9, 84, 85|                      -|      -|

        

## <a name="examples"></a> Examples

### Using mutation calling files (VCFs) as input

```python
import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze

Analyze.cosmic_fit(samples=spa.__path__[0]+"/data/tests/vcf_input", 
                   output="example_vcf",
                   input_type="vcf",
                   context_type="96",
                   genome_build="GRCh37",
                   cosmic_version=3.4)
```


### Using a multi-sample segmentation file as input

```python
import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze

Analyze.cosmic_fit(samples=spa.__path__[0]+"/data/tests/cnv_input/all.breast.ascat.summary.sample.tsv", 
                   output="example_sf",
                   input_type="seg:ASCAT_NGS",
                   cosmic_version=3.4,
                   collapse_to_SBS96=False)
```

### Using a mutational matrix as input

```python
import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze

Analyze.cosmic_fit(samples=spa.__path__[0]+"/data/tests/txt_input/sample_matrix_SBS.txt", 
                   output="example_mm",
                   input_type="matrix",
                   genome_build="GRCh37",
                   cosmic_version=3.4)
```

## <a name="denovo"></a> _De novo_ extraction of mutational signatures downstream analysis
Additional functionalities for downstream analysis of _de novo_ extraction of mutational signatures are also available as part of SigProfilerAssignment, including assignment of _de novo_ extracted mutational signatures and decomposition of _de novo_ signatures using a known set of signatures. More information can be found on the wiki page at https://osf.io/mz79v/wiki/5.%20Advanced%20mode/.

## <a name="citation"></a> Citation

Díaz-Gay, M., Vangara, R., Barnes, M., ... & Alexandrov, L. B. (2023). Assigning mutational signatures to individual samples and individual somatic mutations with SigProfilerAssignment, bioRxiv, 2023-07. doi: https://doi.org/10.1101/2023.07.10.548264

## <a name="copyright"></a> Copyright
This software and its documentation are copyright 2022 as a part of the SigProfiler project. The SigProfilerAssignment framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## <a name="contact"></a> Contact Information
Please address any queries or bug reports to Raviteja Vangara at rvangara@health.ucsd.edu or Marcos Díaz-Gay at mdiazgay@health.ucsd.edu.
