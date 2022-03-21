#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12/2/2021
@author: rvangara
"""

import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
dir_inp = spa.__path__[0]+'/data/Examples/'
def main():
    #data = sig.importdata("text")
    signatures = dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt"
    activities=dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Activities/SBS96_S3_NMF_Activities.txt"
    samples=dir_inp+"Input_scenario_8/Samples.txt"
    output="output_example/"
    sigs= "COSMIC_v3_SBS_GRCh37_noSBS84-85.txt"
    #decomp.spa_analyze( samples, output, signatures=signatures,signature_database=sigs,genome_build="GRCh37", verbose=False,decompose_fit= True,denovo_refit=True,cosmic_fit=True)
    #Analyze.decompose_fit( samples, output, signatures=signatures,signature_database=sigs,genome_build="GRCh37", verbose=False)
    Analyze.denovo_fit( samples, output, signatures=signatures,signature_database=sigs,genome_build="GRCh37", verbose=False)
    #Analyze.cosmic_fit( samples, output, signatures=None,signature_database=sigs,genome_build="GRCh37", verbose=False)



if __name__ == '__main__':
    main()