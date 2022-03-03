# SigProfilerSingleSamplePro (SSPro)

SigProfilerSingleSamplePro is a new mutational attribution and decomposition tool the performs the following functions:
-   Attributing a known set of mutational signatures to an individual sample.
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
$ cd SigProfilerSingleSamplePro-master
$ pip install .
```


Decomposes the De Novo Signatures into COSMIC Signatures and assigns COSMIC signatures into samples

```python
sspro_analyze(  samples,  output, signatures=None, signature_database=None,decompose_fit= True,denovo_refit=True,cosmic_fit=True, nnls_add_penalty=0.05, 
              nnls_remove_penalty=0.01, initial_remove_penalty=0.05, de_novo_fit_penalty=0.02, 
              genome_build="GRCh37",  make_decomposition_plots=True, collapse_to_SBS96=True,connected_sigs=True, verbose=False): 
``` 
<!--> decompose(signatures, activities, samples,  output, signature_database=None, nnls_add_penalty=0.05, nnls_remove_penalty=0.01, initial_remove_penalty=0.05, de_novo_fit_penalty=0.02, genome_build="GRCh37", refit_denovo_signatures=True, make_decomposition_plots=True, connected_sigs=True, verbose=False) -->
| Parameter | Variable Type | Parameter Description |
| --------------------- | -------- |-------- |
| **signatures** | String | Path to a  tab delimited file that contains the signaure table where the rows are mutation types and colunms are signature IDs. |
| **activities** | String | Path to a tab delimilted file that contains the activity table where the rows are sample IDs and colunms are signature IDs. |
| **samples** | String | Path to a tab delimilted file that contains the activity table where the rows are mutation types and colunms are sample IDs. |
| **output** | String | Path to the output folder. |
| **decompose_fit** | Boolean | Perform decomposition of denovo signatures into COSMIC. Default value is True.  |
| **denovo_refit** | Boolean | Refit the given denovo signatures to given samples. Default value is True.  |
| **cosmic_fit** | Boolean | Assignment of mutational profiles from samples to COSMIC or a custom signature database. Default value is True.  |
| **de_novo_fit_penalty** | Float | Takes any positive float. Default is 0.02. Defines the weak (remove) thresh-hold cutoff to be assigned denovo signatures to a sample. Optional parameter. |
| **nnls_add_penalty** | Float | Takes any positive float. Default is 0.01. Defines the weak (remove) thresh-hold cutoff to be assigned COSMIC signatures to a sample. Optional parameter. |
| **nnls_remove_penalty** | Float | Takes any positive float. Default is 0.01. Defines the weak (remove) thresh-hold cutoff to be assigned COSMIC signatures to a sample. Optional parameter. |
| **initial_remove_penalty** | Float | Takes any positive float. Default is 0.05. Defines the initial weak (remove) thresh-hold cutoff to be COSMIC assigned signatures to a sample. Optional parameter. |
| **genome_build** | String | The genome type. Example: "GRCh37", "GRCh38", "mm9", "mm10". The default value is "GRCh37" |
| **verbose** | Boolean | Prints statements. Default value is False.  |
        

#### SSPro analysis Example


```python
import SigProfilerSingleSamplePro as sspro
from SigProfilerSingleSamplePro import decomposition as decomp
dir_inp = sspro.__path__[0]+'/data/Examples/'
signatures = dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt"
activities=dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Activities/SBS96_S3_NMF_Activities.txt"
samples=dir_inp+"/Input_scenario_8/Samples.txt"
output="output_example/"

decomp.sspro_analyze( samples, output, signatures=signatures,genome_build="GRCh37", verbose=False,decompose_fit= True,denovo_refit=True,cosmic_fit=True)
```
## <a name="copyright"></a> Copyright
This software and its documentation are copyright 2022 as a part of the SigProfiler project. The SigProfilerSingleSamplePro framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## <a name="contact"></a> Contact Information
Please address any queries or bug reports to Raviteja Vangara at rvangara@health.ucsd.edu
