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
    #signatures="SBS96_S3_Signatures.txt"
    signatures="SBS288_7_Signatures/Signatures/SBS288_S7_Signatures.txt"
    #signatures = dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt"
    #activities=dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Activities/SBS96_S3_NMF_Activities.txt"
    #samples=dir_inp+"Input_scenario_8/Samples.txt"
    samples="SBS288_7_Signatures/Samples.txt"
    output="output_example_288_2/"
    #sigs= "COSMIC_v3_SBS_GRCh37_noSBS84-85.txt"
    sigs= "COSMIC_v3_duplicate.txt"
    # #decomp.spa_analyze( samples, output, signatures=signatures,signature_database=sigs,genome_build="GRCh37", verbose=False,decompose_fit= True,denovo_refit=True,cosmic_fit=True)
    signature_subgroups = {'permit_MMR_deficiency_signatures' :False,
                      'permit_POL_deficiency_signatures' :False,
                      'permit_HR_deficiency_signatures' :False,
                      'permit_BER_deficiency_signatures' :False,
                      'permit_Chemotherapy_signatures' :True,
                      'permit_APOBEC_signatures' :True,
                      'permit_Tobacco_signatures' :True,
                      'permit_UV_signatures' :True,
                      'permit_AA_signatures' :True,
                      'permit_Colibactin_signatures' :True,
                      'permit_Artifact_signatures' :True,
                      'permit_Lymphoid_signatures}' :True}


    Analyze.decompose_fit( samples, 
                           output, 
                           signatures=signatures,
                           signature_database=sigs,
                           genome_build="GRCh37", 
                           verbose=False,
                           new_signature_thresh_hold=0.8,
                           signature_subgroups=signature_subgroups)
    # Analyze.denovo_fit( samples,
    #                     output, 
    #                     signatures=signatures,
    #                     signature_database=sigs,
    #                     genome_build="GRCh37", 
    #                     verbose=False)
    Analyze.cosmic_fit( samples, 
                        output, 
                        signatures=None,
                        signature_database=sigs,
                        genome_build="GRCh37", 
                        verbose=False,
                        collapse_to_SBS96=False,
                        signature_subgroups=signature_subgroups)

# permit_MMR_deficiency_signatures
# permit_POL_deficiency_signatures
# permit_HR_deficiency_signatures
# permit_BER_deficiency_signatures
# permit_Chemotherapy_signatures
# permit_APOBEC_signatures
# permit_Tobacco_signatures
# permit_UV_signatures
# permit_AA_signatures
# permit_Colibactin_signatures
# permit_Artifact_signatures
# permit_Lymphoid_signatures


if __name__ == '__main__':
    main()  