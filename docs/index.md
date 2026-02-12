SigProfilerAssignment
=====================


----------


SigProfilerAssignment is a [python](https://www.python.org/) framework that assigns previously known mutational signatures to individual samples and individual somatic mutations. The tool refits different types of reference mutational signatures, including COSMIC [SBS][1], [DBS][2], [ID][3], and [CN][4] signatures, as well as custom signature databases. Refitting of known mutational signatures is a numerical optimization approach that not only identifies the set of operative mutational signatures in a particular sample, but also quantifies the number of mutations assigned to each signature found in that sample. SigProfilerAssignment makes use of [SigProfilerMatrixGenerator][5] and [SigProfilerPlotting][6], seamlessly integrating with other [SigProfiler tools][7].

The SigProfilerAssignment library can be found on GitHub [here][8]. For users that prefer working in an R environment, we provide an R wrapper (SigProfilerAssignmentR) that can be found on GitHub [here][9].


----------

### Citation
Díaz-Gay, M., Vangara, R., Barnes, M., ... & Alexandrov, L. B. (2023). Assigning mutational signatures to individual samples and individual somatic mutations with SigProfilerAssignment, bioRxiv, 2023-07. doi: https://doi.org/10.1101/2023.07.10.548264

### License
This software and its documentation are copyright 2022 as a part of the SigProfiler project. The SigProfilerAssignment framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

### Contact
Please address any queries or bug reports to Raviteja Vangara at rvangara@health.ucsd.edu or Marcos Díaz-Gay at mdiazgay@health.ucsd.edu.


  [1]: https://cancer.sanger.ac.uk/cosmic/signatures/SBS/
  [2]: https://cancer.sanger.ac.uk/signatures/dbs/
  [3]: https://cancer.sanger.ac.uk/signatures/id/
  [4]: https://cancer.sanger.ac.uk/signatures/cn/
  [5]: https://osf.io/s93d5/
  [6]: https://osf.io/2aj6t/
  [7]: https://cancer.sanger.ac.uk/signatures/tools/
  [8]: https://github.com/SigProfilerSuite/SigProfilerAssignment/
  [9]: https://github.com/SigProfilerSuite/SigProfilerAssignmentR/
