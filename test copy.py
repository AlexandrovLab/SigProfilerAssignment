#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12/2/2021
@author: rvangara
"""

import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
dir_inp = spa.__path__[0]+'/data/Examples/'
dir_maria= "/Users/rvangara/Documents/collabs/MariaZ/"
def main():
    signatures="SBS288_7_Signatures/Signatures/SBS288_S7_Signatures.txt"
    samples="SBS288_7_Signatures/Samples.txt"
    output="output_test_plots/"
    sigs= "COSMIC_v3_duplicate.txt"

    # #signatures="SBS288_7_Signatures/Signatures/SBS288_S7_Signatures.txt"
    # samples= dir_maria+"Assignment.SBS96.all" #SBS288_7_Signatures/Samples.txt"
    # output=dir_maria+"output_Maria_mm10_v3.2_cosmic_signatures_no17_no18_no36/"
    # sigs= dir_maria+"mm10_v3.2_cosmic_signatures_no17_no18_no36.txt"
    # #sigs= "COSMIC_v3_duplicate.txt"

    # signature_subgroups = [#'remove_MMR_deficiency_signatures',
    #                         #'remove_POL_deficiency_signatures',
    #                         #'remove_HR_deficiency_signatures' ,
    #                         #'remove_BER_deficiency_signatures',
    #                         'remove_Chemotherapy_signatures', 
    #                         'remove_APOBEC_signatures', 
    #                         #'remove_Tobacco_signatures', 
    #                         #'remove_UV_signatures', 
    #                         #'remove_AA_signatures',
    #                         'remove_Colibactin_signatures', 
    #                         'remove_Artifact_signatures', 
    #                         'remove_Lymphoid_signatures']

    # signature_subgroups = {'permit_MMR_deficiency_signatures' :False}
                            # 'permit_POL_deficiency_signatures' :False,
                            # 'permit_HR_deficiency_signatures' :False,
                            # 'permit_BER_deficiency_signatures' :False,
                            # 'permit_Chemotherapy_signatures' :True,
                            # 'permit_APOBEC_signatures' :True,
                            # 'permit_Tobacco_signatures' :True,
                            # 'permit_UV_signatures' :True,
                            # 'permit_AA_signatures' :True,
                            # 'permit_Colibactin_signatures' :True,
                            # 'permit_Artifact_signatures' :True,
                            # 'permit_Lymphoid_signatures}' :True}
    signature_subgroups =None

    # Analyze.decompose_fit( samples, 
    #                        output, 
    #                        signatures=signatures,
    #                        signature_database=sigs,
    #                        genome_build="GRCh37", 
    #                        verbose=False,
    #                        new_signature_thresh_hold=0.8,
    #                        signature_subgroups=signature_subgroups)

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
                        #genome_build="GRCh37", 
                        verbose=False,
                        collapse_to_SBS96=False,
                        signature_subgroups=signature_subgroups,make_plots=False)

if __name__ == '__main__':
    main()  