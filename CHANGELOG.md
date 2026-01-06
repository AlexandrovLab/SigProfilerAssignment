# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.1] - 2026-01-06
### Added
- Added missing COSMIC v3.5 signature files: exome versions for all supported genome builds, mm39 genome build signatures, and rn7 genome build signatures.
- Updated COSMIC v3.5 signature files to match official COSMIC download (includes all signatures: CN26, ID24/ID25, SV11/SV12).

### Fixed
- Fixed missing COSMIC v3.5 signature files that were not included in the initial v3.5 release.

## [1.1.0] - 2026-01-05
### Added
- Support for COSMIC v3.5 mutational signatures as the default reference signature version.
- Added COSMIC v3.5 signature files for all supported genome builds (GRCh37, GRCh38, mm9, mm10, rn6).

### Changed
- Updated default `cosmic_version` parameter from 3.4 to 3.5 across all functions and CLI.
- Updated `Tobacco_signatures` subgroup to include new signatures SBS100 and SBS109.
- Updated README.md to reflect COSMIC v3.5 as the default version and include 3.5 in valid version options.

## [1.0.4] - 2026-01-05
### Added
- Support for `decompose_fit` in higher mutation contexts (288 and 1536) without collapsing to 96 when `collapse_to_SBS96=False`.
- Support for custom signature names in signature decomposition workflows.

### Changed
- Improved context type handling in decomposition plots to correctly pass the mutation context type.

### Fixed
- Fixed variable reusage error in nested for loop that could cause incorrect behavior.
- Fixed missing `rn7` and `mm39` reference signatures in package distribution (now included in MANIFEST.in).
- Updated README.md to include correct `exclude_signature_subgroups` parameter documentation.

## [1.0.3] - 2025-10-31
### Added
- Added `add_background_signatures` parameter (default: `True`) to control whether background signatures SBS1 and SBS5 are automatically added during signature assignment.
- Parameter available in `spa_analyze()`, `signature_decomposition()`, and all Analyzer functions (`decompose_fit()`, `denovo_fit()`, `cosmic_fit()`).
- CLI parameter `--add_background_signatures` added to match Python API functionality.
- When set to `False`, background signatures are not forced but may still be detected naturally if present in samples.

### Changed
- Background signature assignment logic now respects the `add_background_signatures` parameter instead of always forcing SBS1/SBS5 inclusion.

## [1.0.2] - 2025-10-28

### Added
- Implemented a CI/CD pipeline with Travis CI to automate the building and publishing of Docker images to Docker Hub.

## [1.0.1] - 2025-10-20

### Fixed
- Refactored the handling of the `cpu` parameter. The logic has been updated to pass the CPU count directly via the function's parameters, instead of through an internal `devopts` dictionary, to reduce complexity and improve code clarity.

### Added
- Added the `cpu` parameter to the metadata output log.

## [1.0.0] - 2025-10-07

### Added
- Parallel execution of SPA runs for improved performance on multi-core systems.
- `cpu` parameter to control the number of processor cores used during assignment. Defaults to `-1` to use all available cores.

### Documentation
- Updated README to document the new `cpu` parameter.

## [0.2.6] - 2025-09-17

### Added
- Support for assigning mutational signatures using the `rn7` and `mm39` genome builds.

### Security
- Updated dependency requirement to `pypdf>=6.0.0` (previous versions contained a security vulnerability).

## [0.2.5] - 2025-08-07

### Added
- Support for generating decomposition plots for custom signature sets.

### Changed
- Refactored the test script to be more readable and maintainable.

### Fixed
- Removed redundant plotting of COSMIC signatures in the decomposition plots.

## [0.2.4] - 2025-08-07

### Added
- Added rn7 and mm39 reference signatures.

## [0.2.3] - 2025-05-07

### Changed
- Modified rounding exposure values to handle low mutation counts.

## [0.2.2] - 2025-05-02

### Changed
- Replaced PDF-to-PNG conversion backend from `PyMuPDF` to `pdf2image` for compatibility with Conda and improved portability.
- Added new CLI parameter: `--sample_reconstruction_plots` with options `'none'` (default), `'pdf'`, `'png'`, and `'both'`.
- Updated `spa_analyze` and CLI dispatch logic to support format-based sample reconstruction plot output.
- Default behavior now skips sample reconstruction plots unless explicitly requested.
- Removed `fitz` dependency; added system requirement note for `poppler` in `setup.py` and README.

### Added
- Added a pyproject.toml file to the repository for better project management and configuration.

## [0.2.1] - 2025-04-29

### Fixed
- CLI from returning non-zero exit code when --help flag is passed.

### Changed
- Update CI/CD pipelines installation of reference genome to include timeout to prevent long waits during installation.

## [0.2.0] - 2025-02-11

### Changed
- Updated dependencies: Now requires **Pandas >= 2.0.0**, **NumPy >= 2.0.0**, and **Python >= 3.9**.
- Dropped support for **Python 3.8**

## [0.1.9] - 2024-11-12

### Changed
- Replaced `PdfMerger` with `PdfWriter` due to deprecation in `pypdf >= 5.0.0`.

### Fixed
- Addressed deprecation issues with `PdfMerger`, ensuring compatibility with recent `pypdf` versions.

## [0.1.8] - 2024-08-20

### Added
- Added a Dockerfile to the repository for containerization. Documentation on how to use the Dockerfile needs to be added to the README.

### Changed
- Removed unnecessary imports from `setup.py` to clean up the codebase.

## [0.1.7] - 2024-06-21

### Added
- Added CLI boolean handling to improve command-line interface usability.
- Added new pytests for CLI to ensure correct handling of booleans.
- Updated the README to reflect the new pytest additions.
- Updated CI/CD pipelines to accommodate pytest changes.

### Changed
- Improved input_type value check mechanism to use .lower() before checking against values 'vcf' or 'matrix'.
- Updated dependency from PyPDF2 to pypdf to increase compatibility and resolve installation issues on bioconda.
