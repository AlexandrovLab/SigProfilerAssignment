#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 12:21:06 2021

@author: rvangara
"""

from SigProfilerSingleSamplePro import decompose_sub_routines as sub
import numpy as np
import pandas as pd
#import SigProfilerExtractor as cosmic
import os,sys

def Assign( samples,  output, signatures=None, signature_database=None, nnls_add_penalty=0.05, 
              nnls_remove_penalty=0.01, initial_remove_penalty=0.05, 
              genome_build="GRCh37", collapse_to_SBS96=True,connected_sigs=True, verbose=False):

    
    """
    Decomposes the De Novo Signatures into COSMIC Signatures and assigns COSMIC signatures into samples.
    
    Parameters: 
        
        signatures: A string. Path to a  tab delimited file that contains the signaure table where the rows are mutation types and colunms are signature IDs. 
        activities: A string. Path to a tab delimilted file that contains the activity table where the rows are sample IDs and colunms are signature IDs.
        samples: A string. Path to a tab delimilted file that contains the activity table where the rows are mutation types and colunms are sample IDs.
        output: A string. Path to the output folder.
        genome_build = A string. The genome type. Example: "GRCh37", "GRCh38", "mm9", "mm10". The default value is "GRCh37"
        verbose = Boolean. Prints statements. Default value is False. 
        
    Values:
        The files below will be generated in the output folder. 
        
        Cluster_of_Samples.txt
        comparison_with_global_ID_signatures.csv
        Decomposed_Solution_Activities.txt
        Decomposed_Solution_Samples_stats.txt
        Decomposed_Solution_Signatures.txt
        decomposition_logfile.txt
        dendogram.pdf
        Mutation_Probabilities.txt
        Signature_assaignment_logfile.txt
        Signature_plot[MutatutionContext]_plots_Decomposed_Solution.pdf
        
    Example:
        >>>from SigProfilerExtractor import decomposition as decomp
        >>>signatures = "path/to/dDe_Novo_Solution_Signatures.txt"
        >>>activities="path/to/De_Novo_Solution_Activities.txt"
        >>>samples="path/to/Samples.txt"
        >>>output="name or path/to/output.txt"
        decomp.decompose(signatures, activities, samples, output, genome_build="GRCh37", verbose=False)

    """
    
    #processAvg = pd.read_csv(signatures, sep = "\t", index_col=0) 
    #originalProcessAvg=processAvg
    #exposureAvg = pd.read_csv(activities, sep = "\t", index_col = 0)
    genomes = pd.read_csv(samples, sep = "\t", index_col = 0)
    mutation_type = str(genomes.shape[0])
    #m=mutation_type
    index = genomes.index
    colnames = genomes.columns
    #listOfSignatures = processAvg.columns
    #Get COSMIC SIGNATURES
    #import pdb
    if signatures == None:
        processAvg = sub.getProcessAvg(genomes, genome_build, "3.2")
        processAvg = processAvg.set_index('Type').rename_axis('MutationType')
    else:
        #pdb.set_trace()
        try:
            processAvg = pd.read_csv(signatures,sep='\t', index_col=0)
        except:
            sys.exit("Something is wrong with input signatures, Pass a text file of signatures in the format of COSMIC sig database")
    #pdb.set_trace()

    ##

    #creating list of mutational type to sync with the vcf type input
    if mutation_type == "78":
        mutation_context = "DBS78"
    elif mutation_type == "83":
        mutation_context = "ID83"
    elif mutation_type=="48":
        mutation_context = "CNV48"
        print("Mutation Type is: CNV")
    else:
        mutation_context = "SBS"+mutation_type
        
    #processAvg = np.array(processAvg)    
    #signature_names = sub.make_letter_ids(idlenth = processAvg.shape[1], mtype = mutation_context)
    #exposureAvg.columns=signature_names   
    # create the folder for the final solution/ De Novo Solution
  
    try:
        if not os.path.exists(output):
            os.makedirs(output)
    except: 
        print ("The {} folder could not be created".format("output"))
    
    # make the texts for signature plotting
    layer_directory1 = output+"/Assignment_Solution"
    try:
        if not os.path.exists(layer_directory1):
            os.makedirs(layer_directory1)
    except: 
        print ("The {} folder could not be created".format("Assignment_Solution"))
    import pdb
    # pdb.set_trace() 
    # if processAvg.shape[0]==1536 and collapse_to_SBS96==True: #collapse the 1596 context into 96 only for the deocmposition 
    #     processAvg = pd.DataFrame(processAvg, index=index)
    #     processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
    #     genomes = genomes.groupby(genomes.index.str[1:8]).sum()
    #     index = genomes.index
    #     processAvg = np.array(processAvg)
    
    # if processAvg.shape[0]==288 and collapse_to_SBS96==True: #collapse the 288 context into 96 only for the deocmposition 
    #     processAvg = pd.DataFrame(processAvg, index=index)
    #     processAvg = processAvg.groupby(processAvg.index.str[2:9]).sum()
    #     genomes = pd.DataFrame(genomes, index=index)
    #     genomes = genomes.groupby(genomes.index.str[2:9]).sum()
    #     index = genomes.index
    #     processAvg = np.array(processAvg)
    if genomes.shape[0]==1536 and collapse_to_SBS96==True: #collapse the 1596 context into 96 only for the deocmposition 
        #processAvg = pd.DataFrame(processAvg, index=index)
        #processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
        genomes = genomes.groupby(genomes.index.str[1:8]).sum()
        index = genomes.index
        #processAvg = np.array(processAvg)
    
    if genomes.shape[0]==288 and collapse_to_SBS96==True: #collapse the 288 context into 96 only for the deocmposition 
        processAvg = pd.DataFrame(processAvg, index=index)
        processAvg = processAvg.groupby(processAvg.index.str[2:9]).sum()
        genomes = pd.DataFrame(genomes, index=index)
        genomes = genomes.groupby(genomes.index.str[2:9]).sum()
        index = genomes.index
        #processAvg = np.array(processAvg)
            
    
    index = genomes.index
    colnames = genomes.columns
    
    
    allsigids = processAvg.columns.to_list()
    processAvg = processAvg.values
    
    attribution={}
    for i in allsigids:
        attribution[i]= [i] 
    #only for SBS96
    if mutation_type == "96" or mutation_type=="288" or mutation_type=="1536":        
        background_sigs = sub.get_indeces(list(allsigids), ['SBS1', 'SBS5'])
        # add connected signatures   
        #different_signatures = ss.add_connected_sigs(different_signatures, list(signames))
    #for other contexts
    else:
        background_sigs = []
    exposureAvg_dummy = pd.DataFrame(np.random.rand(processAvg.shape[1],genomes.shape[1]),index=allsigids,columns=colnames.to_list()).transpose().rename_axis('Samples')
  
    result = sub.make_final_solution(processAvg, genomes, allsigids, layer_directory1, mutation_type, index, colnames, 
                            cosmic_sigs=False, attribution = attribution, denovo_exposureAvg  = exposureAvg_dummy ,  
                            background_sigs=background_sigs, verbose=verbose, genome_build=genome_build, 
                            add_penalty=nnls_add_penalty, remove_penalty=nnls_remove_penalty, 
                            initial_remove_penalty=initial_remove_penalty,connected_sigs=connected_sigs,
                            collapse_to_SBS96=collapse_to_SBS96,
                            refit_denovo_signatures=True)
   
    return result


