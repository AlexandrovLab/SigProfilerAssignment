<h1> Quick Start Example </h1>

----------

This section provides an example for users to quickly get started with using the SigProfilerAssignment tool. The following example will use somatic mutational data from breast cancer samples from [Nik-Zainal et al. 2012 Cell][1], and will showcase how to use SigProfilerAssignment with different types of files containing the input somatic mutations, including variant calling files (VCFs) and mutational matrices.


@[toc](Sections)

----------

## Prerequisites ##
This tutorial requires that you have completed all steps in the [installation guide][2], specifically:

 - Installed SigProfilerAssignment
 - Downloaded **GRCh37** reference genome using SigProfilerMatrixGenerator


## Downloading Input Example Data ##
This example uses somatic mutational data from a breast cancer genome. Download the example dataset `BRCA.zip` at the following location or use the command line:

    ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerAssignment/Example_data/
    
If using the command line, then enter the following command in bash on OS X or Unix systems:

    $ wget ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerAssignment/Example_data/BRCA.zip
    
Once `BRCA.zip` has been downloaded, unzip the file. The unzipped `BRCA` folder contains `BRCA.txt` and another folder `BRCA_vcf`. The file `BRCA.txt` is a mutational matrix defined using SBS-96 classification (created by [SigProfilerMatrixGenerator][3]) and `BRCA_vcf` contains the corresponding VCF file associated to the sample.


## Running SigProfilerAssignment (VCF) ##
You will be assigning reference mutational signatures from [COSMIC][4] v3.3 to the breast cancer sample in the subfolder `BRCA_vcf` used as input for this example.

First, start a Python interactive shell and import the SigProfilerAssignment library.

``` python
$ python
>>> from SigProfilerAssignment import Analyzer as Analyze
```

Next, assign reference COSMIC signatures by running the following command. **Note**: Update `"path/to/BRCA_vcf"` with the actual path to the `BRCA_vcf` folder.

``` python
Analyze.cosmic_fit(samples="path/to/BRCA_vcf", 
                   output="output_vcf",
                   input_type="vcf",
                   context_type="96",
                   genome_build="GRCh37")
```

After SigProfilerAssignment has finished running, an output directory name `output_vcf` will be created. This directory will contain the output files and is located in the directory where the Python instance was started. To learn more about the output produced by SigProfilerAssignment, please refer to the [Using the Tool - Output][5] section.

## Running SigProfilerAssignment (Mutational matrix) ##
You will be assigning reference mutational signatures from [COSMIC][6] v3.3 to the mutational matrix defined using the SBS-96 classification named `BRCA.txt` input for this example.

First, start a Python interactive shell and import the SigProfilerAssignment library.

``` python
$ python
>>> from SigProfilerAssignment import Analyzer as Analyze
```

Next, assign reference COSMIC signatures by running the following command. **Note**: Update `"path/to/BRCA.txt"` with the actual path to the `BRCA.txt` file.

``` python
Analyze.cosmic_fit(samples="path/to/BRCA.txt", 
                   output="output_mm",
                   input_type="matrix")
```

After SigProfilerAssignment has finished running, an output directory name `output_mm` will be created. This directory will contain the output files and is located in the directory where the Python instance was started. To lear more about the output produced by SigProfilerAssignment, please refer to the [Using the Tool - Output][5] section.

## Running SigProfilerAssignment (Multi-sample segmentation) ##
You will be assigning reference mutational signatures from [COSMIC][8] v3.3 to the multi-sample segmentation file obtained from one of the copy number calling tools named `all.breast.ascat.summary.sample.tsv` input for this example.

First, start a Python interactive shell and import the SigProfilerAssignment library.

``` python
$ python
>>> from SigProfilerAssignment import Analyzer as Analyze
```

Next, assign reference COSMIC signatures by running the following command. **Note**: Update `"path/to/all.breast.ascat.summary.sample.tsv"` with the actual path to the `all.breast.ascat.summary.sample.tsv` file.

``` python
Analyze.cosmic_fit(samples="path/to/all.breast.ascat.summary.sample.tsv", 
                   output="example_sf",
                   input_type="seg:ASCAT_NGS",
                   cosmic_version=3.3,
                   collapse_to_SBS96=False)
```

After SigProfilerAssignment has finished running, an output directory name `example_sf` will be created. This directory will contain the output files and is located in the directory where the Python instance was started. To lear more about the output produced by SigProfilerAssignment, please refer to the [Using the Tool - Output][5] section.

## Additional Information ##
In the above examples, the other non specified parameters are passed in with their default values. All of the function arguments and their types are explained in detail in the [Using the Tool - Input section][7]. To learn more about the files that were produced, you can refer to [Using the Tool - Output][5].
  


  [1]: https://doi.org/10.1016/j.cell.2012.04.024
  [2]: https://osf.io/mz79v/wiki/1.Installation/
  [3]: https://osf.io/s93d5/wiki/home
  [4]: https://cancer.sanger.ac.uk/signatures/
  [5]: https://osf.io/mz79v/wiki/4.%20Using%20the%20Tool%20-%20Output/
  [6]: https://cancer.sanger.ac.uk/signatures/
  [7]: https://osf.io/mz79v/wiki/3.Using%20the%20Tool%20-%20Input/
  [8]: https://cancer.sanger.ac.uk/signatures/
