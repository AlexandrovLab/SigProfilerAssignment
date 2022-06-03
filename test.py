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
    #signatures=dir_inp+"SBS288_7_Signatures/Signatures/SBS288_S7_Signatures.txt"
    # samples=dir_inp+"SBS288_7_Signatures/Samples.txt"
    signatures=dir_inp+"SBS1536_5_Signatures/Signatures/SBS1536_S5_Signatures.txt"
    samples=dir_inp+"SBS1536_5_Signatures/Samples.txt"
    output="output_example/"
    # sigs= "COSMIC_v3_SBS_GRCh37_noSBS84-85.txt"


    signature_subgroups = [#'remove_MMR_deficiency_signatures',
                            #'remove_POL_deficiency_signatures',
                            #'remove_HR_deficiency_signatures' ,
                            #'remove_BER_deficiency_signatures',
                            'remove_Chemotherapy_signatures', 
                            'remove_APOBEC_signatures', 
                            #'remove_Tobacco_signatures', 
                            #'remove_UV_signatures', 
                            #'remove_AA_signatures',
                            'remove_Colibactin_signatures', 
                            'remove_Artifact_signatures', 
                            'remove_Lymphoid_signatures']

    signature_subgroups = None

    Analyze.decompose_fit( samples, 
                           output, 
                           signatures=signatures,
                           signature_database=None,
                           genome_build="GRCh37", 
                           verbose=False,
                           new_signature_thresh_hold=0.8,
                           signature_subgroups=signature_subgroups)

    Analyze.denovo_fit( samples,
                        output, 
                        signatures=signatures,
                        signature_database=None,
                        genome_build="GRCh37", 
                        verbose=False)

    Analyze.cosmic_fit( samples, 
                        output, 
                        signatures=None,
                        signature_database=None,
                        genome_build="GRCh37", 
                        verbose=False,
                        collapse_to_SBS96=True,
                        signature_subgroups=signature_subgroups,make_plots=False)

if __name__ == '__main__':
    main()  