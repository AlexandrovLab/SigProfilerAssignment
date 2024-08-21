# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
