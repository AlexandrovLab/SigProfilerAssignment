# SigProfilerAssignment

![Logo](assets/images/SigProfilerAssignment.png)

----------

**SigProfilerAssignment** is a Python framework that assigns and decomposes mutational signatures in individual samples and individual somatic mutations. 

Refitting of  mutational signatures is a numerical optimization approach that not only identifies the set of operative mutational signatures in a particular sample, but also quantifies the number of mutations attributed to each signature detected in that sample. In addition to refitting reference signatures, the tool can decompose _de novo_ extracted signatures into reference signatures, facilitating biological interpretation and comparison with established catalogs. 

It supports multiple mutation contexts and genome builds, provides confidence estimation and signature activity thresholds, and generates detailed visual and tabular outputs to aid downstream analyses. The tool refits different types of known reference mutational signatures, including COSMIC [SBS][1], [DBS][2], [ID][3], and [CN][4] signatures, as well as custom signature databases. 

**SigProfilerAssignment** makes use of [SigProfilerMatrixGenerator][5] and [SigProfilerPlotting][6], enabling seamless integration with other tools in the SigProfiler suite.

The SigProfilerAssignment library is available on [GitHub](https://github.com/SigProfilerSuite/SigProfilerAssignment) and [PyPI](https://pypi.org/project/SigProfilerAssignment)


----------

### Citation

Díaz-Gay M, Vangara R, Barnes M, Wang X, Islam SMA, Vermes I, Duke S, Narasimman NB, Yang T, Jiang Z, Moody S, Senkin S, Brennan P, Stratton MR, Alexandrov LB. Assigning mutational signatures to individual samples and individual somatic mutations with SigProfilerAssignment. *Bioinformatics*. 2023;39(12):btad756. [https://doi.org/10.1093/bioinformatics/btad756](https://doi.org/10.1093/bioinformatics/btad756)

### License

This software and its documentation are part of the SigProfiler project and are copyrighted © 2022. The SigProfilerAssignment framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

### Contact

For questions, support requests, or bug reports, please contact the SigProfilerSuite team via GitHub [issues](https://github.com/SigProfilerSuite/SigProfilerAssignment/issues) or by email at [contact@sigprofilersuite.org](mailto:contact@sigprofilersuite.org).

  [1]: https://cancer.sanger.ac.uk/cosmic/signatures/SBS/
  [2]: https://cancer.sanger.ac.uk/signatures/dbs/
  [3]: https://cancer.sanger.ac.uk/signatures/id/
  [4]: https://cancer.sanger.ac.uk/signatures/cn/
  [5]: https://osf.io/s93d5/
  [6]: https://osf.io/2aj6t/
  [7]: https://cancer.sanger.ac.uk/signatures/tools/
  [8]: https://github.com/SigProfilerSuite/SigProfilerAssignment/
  [9]: https://github.com/SigProfilerSuite/SigProfilerAssignmentR/
