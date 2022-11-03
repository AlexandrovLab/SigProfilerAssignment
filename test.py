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
OUTPUT_TXT = "output_txt_example/"
OUTPUT_VCF = "output_vcf_example/"
SIGNATURES_TXT = os.path.join(DIR_INPUT, "txt_input/SBS96_S3_Signatures.txt")
SAMPLES_TXT = os.path.join(DIR_INPUT, "txt_input/sample_matrix.txt")
SAMPLES_VCF = os.path.join(DIR_INPUT, "vcf_input/")
exclude_signature_subgroups = None

def cosmic_fit_txt_test():
    Analyze.cosmic_fit(SAMPLES_TXT,
        OUTPUT_TXT,
        signatures=None,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        collapse_to_SBS96=True,
        exclude_signature_subgroups=exclude_signature_subgroups)

def denovo_fit_txt_test():
    Analyze.denovo_fit(SAMPLES_TXT,
        OUTPUT_TXT,
        signatures=SIGNATURES_TXT,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False)

def decompose_fit_txt_test():
    Analyze.decompose_fit(SAMPLES_TXT,
        OUTPUT_TXT,
        signatures=SIGNATURES_TXT,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        new_signature_thresh_hold=0.8,
        exclude_signature_subgroups=exclude_signature_subgroups)

def cosmic_fit_vcf_test():
    Analyze.cosmic_fit(SAMPLES_VCF,
        OUTPUT_VCF,
        input_type="vcf",
        signatures=None,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        collapse_to_SBS96=True,
        exclude_signature_subgroups=exclude_signature_subgroups,
        make_plots=True)

def decompose_fit_vcf_test():
    Analyze.decompose_fit(SAMPLES_VCF,
        OUTPUT_VCF,
        input_type="vcf",
        signatures=SIGNATURES_TXT,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        new_signature_thresh_hold=0.8,
        exclude_signature_subgroups=exclude_signature_subgroups,
        make_plots=True)

def denovo_fit_vcf_test():
    Analyze.denovo_fit(SAMPLES_VCF,
        OUTPUT_VCF,
        signatures=SIGNATURES_TXT,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        make_plots=True,
        input_type="vcf")

if __name__ == '__main__':
    cosmic_fit_txt_test()
    denovo_fit_txt_test()
    decompose_fit_txt_test()
    cosmic_fit_vcf_test()
    decompose_fit_vcf_test()
    denovo_fit_vcf_test()