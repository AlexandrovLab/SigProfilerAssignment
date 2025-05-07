# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
