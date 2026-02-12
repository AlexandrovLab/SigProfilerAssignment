Advanced Mode Examples
=====================


----------


This section provides examples for advanced users to get started with using SigProfilerAssignment for *de novo* extraction of mutational signautres downstream analysis, including:

1. Assignment of *de novo* extracted mutational signatures using `denovo_fit`.
2. Decomposition of *de novo* extracted mutational signatures using a known set of signatures (reference [COSMIC][5] signatures or customized signature databases) by `decompose_fit`.

----------



## `denovo_fit` ##
Attributes the somatic mutations of a given sample/s to a set of input *de novo* signatures. 

Two input files are required. First, a file containing the input somatic mutations, in any of the formats specified in the [Using the Tool - Input][1] section. Also, a matrix containing the *de novo* extracted mutational signatures, commonly derived from [SigProfilerExtractor][2]. 

**Note**: A reference genome build should also be specified if a mutation calling file is used as input.


```
from SigProfilerAssignment import Analyzer as Analyze
Analyze.denovo_fit(samples="path/to/input/mutations/file",
                   output="path/to/output/folder",
                   input_type="desired/input/mutation/file/format",
                   signatures="path/to/input/denovo/signatures/file",
                   genome_build="GRCh37")
```

## `decompose_fit` ##
Decomposes a set of *de novo* extracted mutational signatures into a known set of signatures (reference COSMIC signatures or customized signature databases) and assigns these known signatures into a given sample/s.

Two input files are required. First, a file containing the input somatic mutations, in any of the formats specified in the [Using the Tool - Input][3] section. Also, a matrix containing the *de novo* extracted mutational signatures, commonly derived from [SigProfilerExtractor][4]. 

An optional third input file is needed in case a custom reference signature database is used. By default, reference [COSMIC][5] signatures v3.3 are used for decomposing the set of *de novo* extracted signatures.

**Note**: A reference genome build should also be specified if a mutation calling file is used as input.

```
from SigProfilerAssignment import Analyzer as Analyze
Analyze.decompose_fit(samples="path/to/input/mutations/file", 
                      output="path/to/output/folder",
                      input_type="desired/input/mutation/file/format",
                      signatures="path/to/input/denovo/signatures/file",
                      [signature_database="path/to/optional/reference/signatures/database/file",]
                      genome_build="GRCh37")
```


  [1]: https://osf.io/mz79v/wiki/3.Using%20the%20Tool%20-%20Input/
  [2]: https://osf.io/t6j7u/wiki/home/
  [3]: https://osf.io/mz79v/wiki/3.Using%20the%20Tool%20-%20Input/
  [4]: https://osf.io/t6j7u/wiki/home/
  [5]: https://cancer.sanger.ac.uk/signatures/
