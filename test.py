#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12/2/2021
@author: rvangara
"""
import SigProfilerSingleSamplePro as sspro
from SigProfilerSingleSamplePro import decomposition as decomp
dir_inp = sspro.__path__[0]+'/data/Examples/'
def main():
    #data = sig.importdata("text")
    signatures = dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt"
    activities=dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Activities/SBS96_S3_NMF_Activities.txt"
    samples=dir_inp+"/Input_scenario_8/Samples.txt"
    output="output_decompose_test/"

    decomp.Assign(samples,output, genome_build="GRCh37", verbose=False)

if __name__ == '__main__':
    main()
