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

    signatures  = dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt"
    samples     = dir_inp+"Input_scenario_8/Samples.txt"
    # samples = spa.__path__[0]+'/data/vcftest/' If input is a directory of  vcf files. 
    output="output_example/"
    sigs= "COSMIC_v3_SBS_GRCh37_noSBS84-85.txt"
    
    
    # vcf_opts={'project_name': 'test_sample','vcf_context': '288' } # Uncomment this  If  vcf files are provided as input. 

    # signature_subgroups = ['remove_MMR_deficiency_signatures',
    #                         'remove_POL_deficiency_signatures',
    #                         'remove_HR_deficiency_signatures' ,
    #                         'remove_BER_deficiency_signatures',
    #                         'remove_Chemotherapy_signatures', 
    #                         'remove_APOBEC_signatures', 
    #                         'remove_Tobacco_signatures', 
    #                         'remove_UV_signatures', 
    #                         'remove_AA_signatures',
    #                         'remove_Colibactin_signatures', 
    #                         'remove_Artifact_signatures', 
    #                         'remove_Lymphoid_signatures']

    signature_subgroups = None

    Analyze.decompose_fit( samples, 
                           output, 
                           signatures=signatures,
                           signature_database=None,
                           genome_build="GRCh37", 
                           verbose=False,
                           new_signature_thresh_hold=0.8,
                           signature_subgroups=signature_subgroups,
                        #    vcf_opts=vcf_opts
                           )

    Analyze.denovo_fit( samples,
                        output, 
                        signatures=signatures,
                        signature_database=None,
                        genome_build="GRCh37", 
                        verbose=False,
                        # vcf_opts=vcf_opts
                        )

    Analyze.cosmic_fit( samples, 
                        output, 
                        signatures=None,
                        signature_database=None,
                        genome_build="GRCh37", 
                        verbose=False,
                        collapse_to_SBS96=True,
                        signature_subgroups=signature_subgroups,
                        # vcf_opts=vcf_opts
                       )

if __name__ == '__main__':
    main()  