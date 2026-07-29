"""
Unit tests for cosmic_fit function.
Run with: pytest tests/test_cosmic_fit.py -v
"""

import os
import pandas as pd
import pytest

from SigProfilerAssignment import decompose_subroutines as sub

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
DIR_INPUT = os.path.join(os.path.dirname(TESTS_DIR), "SigProfilerAssignment", "data", "tests")
TXT_INPUT_DIR = os.path.join(DIR_INPUT, "txt_input")

MATRICES = {"SBS": os.path.join(TXT_INPUT_DIR, "sample_matrix_SBS.txt"),
            "DBS": os.path.join(TXT_INPUT_DIR, "sample_matrix_DBS.txt"),
            "ID": os.path.join(TXT_INPUT_DIR, "sample_matrix_ID.txt"),
            "SV": os.path.join(TXT_INPUT_DIR, "sample_matrix_SV32.txt"),
            "CNV": os.path.join(TXT_INPUT_DIR, "sample_matrix_CNV48.txt")}
GENOMES = ["GRCh37", "GRCh38", "mm9", "mm10", "mm39", "rn6", "rn7"]
COSMIC_VERSIONS = [1, 2, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6]
EXOMES = [False, True]

class TestCosmic:
    @pytest.mark.parametrize("mtype", MATRICES.keys())
    @pytest.mark.parametrize("genome_build", GENOMES)
    @pytest.mark.parametrize("cosmic_version", COSMIC_VERSIONS)
    @pytest.mark.parametrize("exome", EXOMES)
    def test_file_exists(self, mtype, genome_build, cosmic_version, exome):

        genomes = pd.read_csv(
            MATRICES[mtype],
            sep="\t",
            index_col=0,
        )

        signature_database = sub.getProcessAvg(
            genomes,
            genome_build=genome_build,
            cosmic_version=cosmic_version,
            exome=exome,
        )[0]

        assert isinstance(signature_database, pd.DataFrame)