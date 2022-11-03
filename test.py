#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12/2/2021
@author: rvangara
"""

import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
import os
DIR_INPUT = spa.__path__[0]+'/data/tests/'
OUTPUT = "output_example/"
TXT_SIGNATURES = os.path.join(DIR_INPUT, "txt_input/SBS96_S3_Signatures.txt")
TXT_SAMPLES = os.path.join(DIR_INPUT, "txt_input/sample_matrix.txt")
exclude_signature_subgroups = None

def cosmic_fit_test():
    Analyze.cosmic_fit(TXT_SAMPLES,
        OUTPUT,
        signatures=None,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        collapse_to_SBS96=True,
        exclude_signature_subgroups=exclude_signature_subgroups)

def denovo_fit_test():
    Analyze.denovo_fit(TXT_SAMPLES,
        OUTPUT,
        signatures=TXT_SIGNATURES,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False)

def decompose_fit_test():
    Analyze.decompose_fit(TXT_SAMPLES,
        OUTPUT,
        signatures=TXT_SIGNATURES,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        new_signature_thresh_hold=0.8,
        exclude_signature_subgroups=exclude_signature_subgroups)

if __name__ == '__main__':
    cosmic_fit_test()
    denovo_fit_test()
    decompose_fit_test()