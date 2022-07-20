#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 19 12:21:06 2019

@author: mishugeb
"""

#from SigProfilerExtractor import subroutines as sub

from cmath import cos

# from torch import sign
from SigProfilerAssignment import decompose_sub_routines as sub
import numpy as np
import pandas as pd

#import SigProfilerExtractor as cosmic
import os,sys


def spa_analyze(  samples,  output, signatures=None, signature_database=None,decompose_fit_option= True,denovo_refit_option=True,cosmic_fit_option=True, nnls_add_penalty=0.05, 
              nnls_remove_penalty=0.01, initial_remove_penalty=0.05, de_novo_fit_penalty=0.02, 
              genome_build="GRCh37", cosmic_version=3.3, make_plots=True, collapse_to_SBS96=True,connected_sigs=True, verbose=False,devopts=None,new_signature_thresh_hold=0.8,signature_subgroups=None):

    
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
    if devopts == None:
        layer_directory1 = output+"/De_Novo_Solution"
        layer_directory2 = output+"/Decompose_Solution"
        layer_directory3 = output+"/Assignment_Solution"
    else:
        layer_directory1 = devopts['denovo_outpath']
        layer_directory2 = devopts['decompose_outpath']
        layer_directory3 = devopts['Assignment_outpath']


    if (denovo_refit_option == True or decompose_fit_option ==True) and signatures is None:
         raise Exception("If denovo_refit or decompose_fit is True, signatures cannot be empty")
    
    try:
        genomes = pd.read_csv(samples, sep = "\t", index_col = 0)
    except:
        genomes = samples
        genomes = pd.DataFrame(genomes)
        
    # if signatures is None:
    #     processAvg = sub.getProcessAvg(genomes, genome_build=genome_build, cosmic_version=cosmic_version)[0]
    #     processAvg = processAvg.rename_axis('MutationType')
    #     #processAvg = processAvg.set_index('Type').rename_axis('MutationType')   
    # else:
    #     try:
    #         processAvg = pd.read_csv(signatures,sep='\t', index_col=0)
    #     except:
    #         try:
    #             processAvg=signatures
    #         except:
    #             sys.exit("Error in formatting of input signatures, Pass a text file of signatures in the format of COSMIC sig database")
    
    default_subgroups_dict= {'remove_MMR_deficiency_signatures' :False,
                      'remove_POL_deficiency_signatures' :False,
                      'remove_HR_deficiency_signatures' :False,
                      'remove_BER_deficiency_signatures' :False,
                      'remove_Chemotherapy_signatures' :False,
                      'remove_Immunosuppressants_signatures' :False,
                      'remove_Iatrogenic_signatures' :False,
                      'remove_APOBEC_signatures' :False,
                      'remove_Tobacco_signatures' :False,
                      'remove_UV_signatures' :False,
                      'remove_AA_signatures' :False,
                      'remove_Colibactin_signatures' :False,
                      'remove_Artifact_signatures' :False,
                      'remove_Lymphoid_signatures' :False}
                      
    default_subgroups_siglists= {'remove_MMR_deficiency_signatures' :['6', '14', '15', '20', '21', '26', '44'],
                      'remove_POL_deficiency_signatures' :['10a', '10b', '10c', '10d', '28'],
                      'remove_HR_deficiency_signatures' :['3'],
                      'remove_BER_deficiency_signatures' :['30','36'],
                      'remove_Chemotherapy_signatures' :['11','25','31','35','86','87','90'],
                      'remove_Immunosuppressants_signatures' :['32'],
                      'remove_Iatrogenic_signatures' :['11','25','31','32','35','86','87','90'],
                      'remove_APOBEC_signatures' :['2','13'],
                      'remove_Tobacco_signatures' :['4','29','92'],
                      'remove_UV_signatures' :['7a','7b','7c','7d','38'],
                      'remove_AA_signatures' :['22'],
                      'remove_Colibactin_signatures' :['88'],
                      'remove_Artifact_signatures' :['27','43','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60'],
                      'remove_Lymphoid_signatures' :['9','84','85']}
    
    
    signature_subgroups_dict = default_subgroups_dict.copy()
    if signature_subgroups == None:
        pass
    else:
        if type(signature_subgroups) is not list:
            sys.exit("signature_subgroups input should be a list of appropriate flags, please refer to documentation.")
        else:

            for key in default_subgroups_dict:
                if key in signature_subgroups:
                    signature_subgroups_dict[key]=True

    sig_exclusion_list=[]
    if signature_subgroups == None:
        sig_exclusion_list=[]
    else:
        for key in signature_subgroups_dict:
            if signature_subgroups_dict[key]:
                sig_exclusion_list.append(default_subgroups_siglists[key])
    
    sig_exclusion_list = [item for sublist in sig_exclusion_list for item in sublist]
    
    mutation_type = str(genomes.shape[0])
    m=mutation_type

    try:
        if not os.path.exists(output):
            os.makedirs(output)
    except:
        print ("The {} folder could not be created".format("output"))
  

                                                                #################
                                                                # Denovo refiting #
                                                                #################
    
    if denovo_refit_option == True:
        try:
            processAvg = pd.read_csv(signatures,sep='\t', index_col=0)
        except:
            try:
                processAvg=signatures
            except:
                sys.exit("Error in formatting of input signatures, Pass a text file of signatures in the format of denovo signatures")
          
        if devopts == None:
            listOfSignatures = processAvg.columns
            index = genomes.index
            colnames = genomes.columns
        else:
            listOfSignatures=devopts['listOfSignatures']
            index=devopts['index']
            colnames=devopts['colnames']
            genomes = genomes.set_index(index)
            genomes.columns = colnames
            #genomes = genomes.rename_axis("Mutation Types", axis="columns")
       
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
        try:
            allsigids = processAvg.columns.to_list()
        except:
            allsigids = list(listOfSignatures)
        processAvg = np.array(processAvg)    
        signature_names = sub.make_letter_ids(idlenth = processAvg.shape[1], mtype = mutation_context)
       
        # create the folder for the final solution/ De Novo Solution
        #pdb.set_trace()
        exposureAvg_dummy = pd.DataFrame(np.random.rand(processAvg.shape[1],genomes.shape[1]),index=listOfSignatures,columns=colnames.to_list()).transpose().rename_axis('Samples')
        exposureAvg = exposureAvg_dummy
        exposureAvg.columns=signature_names   

        

        refit_denovo_signatures= True
        init_rem_denovo = 0.0
        
        # make the texts for signature plotting
        #layer_directory1 = output+"/De_Novo_Solution"
        try:
            if not os.path.exists(layer_directory1):
                os.makedirs(layer_directory1)
        except: 
            print ("The {} folder could not be created".format("De_Novo_Solution"))
    
        listOfSignatures = sub.make_letter_ids(idlenth = processAvg.shape[1], mtype=mutation_context)
        genomes = pd.DataFrame(genomes)
        denovo_exposureAvg = np.array(exposureAvg.T)
        print("\n Denovo Fitting .....")
        #exposureAvg = sub.make_final_solution(processAvg, genomes, listOfSignatures, layer_directory1, mutation_type, index,\
                   # colnames,denovo_exposureAvg  = denovo_exposureAvg, add_penalty=nnls_add_penalty, remove_penalty=nnls_remove_penalty, initial_remove_penalty=initial_remove_penalty, de_novo_fit_penalty=de_novo_fit_penalty, connected_sigs=connected_sigs, refit_denovo_signatures=refit_denovo_signatures)    
                   # 
        #######
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

        if devopts == None:
            exposureAvg = sub.make_final_solution(processAvg, genomes, allsigids, layer_directory1, mutation_type, index, colnames, 
                                    cosmic_sigs=True, attribution = attribution, denovo_exposureAvg  = exposureAvg_dummy , 
                                    background_sigs=background_sigs, verbose=verbose, genome_build=genome_build, 
                                    add_penalty=nnls_add_penalty, remove_penalty=nnls_remove_penalty,
                                    initial_remove_penalty=init_rem_denovo,connected_sigs=connected_sigs,refit_denovo_signatures=False,make_plots=make_plots)

        else:
            signature_stabilities=devopts['signature_stabilities']
            signature_total_mutations=devopts['signature_total_mutations']
            signature_stats = devopts['signature_stats']
            sequence=devopts['sequence']
            processSTE=devopts['processSTE']
            sequence =devopts['sequence']

            exposureAvg = sub.make_final_solution(processAvg, genomes, allsigids, layer_directory1, mutation_type, index, colnames, 
                                    cosmic_sigs=True, attribution = attribution, denovo_exposureAvg  = exposureAvg_dummy , sequence=sequence, 
                                    background_sigs=background_sigs, verbose=verbose, genome_build=genome_build, signature_total_mutations = signature_total_mutations,
                                    add_penalty=nnls_add_penalty, remove_penalty=nnls_remove_penalty, process_std_error = processSTE, signature_stabilities = signature_stabilities,
                                    initial_remove_penalty=init_rem_denovo,connected_sigs=connected_sigs,refit_denovo_signatures=True)


                                                                #################
                                                                # Decomposition       
                                                                #################
    if decompose_fit_option ==True:
        try:
            processAvg = pd.read_csv(signatures,sep='\t', index_col=0)
        except:
            try:
                processAvg=signatures
            except:
                sys.exit("Error in formatting of input signatures, Pass a text file of signatures in the format of denovo signatures")

        if devopts == None:
            listOfSignatures = processAvg.columns
            index = genomes.index
            colnames = genomes.columns
        else:
            listOfSignatures=devopts['listOfSignatures']
            index=devopts['index']
            colnames=devopts['colnames']
            genomes = genomes.set_index(index)
            genomes.columns = colnames
            #genomes = genomes.rename_axis("Mutation Types", axis="columns")

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
        try:
            allsigids = processAvg.columns.to_list()
        except:
            allsigids = list(listOfSignatures)
        processAvg = np.array(processAvg)    
        signature_names = sub.make_letter_ids(idlenth = processAvg.shape[1], mtype = mutation_context)

        exposureAvg_dummy = pd.DataFrame(np.random.rand(processAvg.shape[1],genomes.shape[1]),index=listOfSignatures,columns=colnames.to_list()).transpose().rename_axis('Samples')
        exposureAvg = exposureAvg_dummy
        exposureAvg.columns=signature_names   


        #############################
        #layer_directory2 = output+"/Decompose_Solution"
        if isinstance(processAvg, pd.DataFrame):
            pass
        else:
            originalProcessAvg=pd.DataFrame(processAvg,index=index,columns=listOfSignatures)
        try:
            if not os.path.exists(layer_directory2):
                os.makedirs(layer_directory2)
        except: 
            print ("The {} folder could not be created".format("Decomposed_Solution"))

        if processAvg.shape[0]==1536 and collapse_to_SBS96==True: #collapse the 1596 context into 96 only for the deocmposition 
            processAvg = pd.DataFrame(processAvg, index=index)
            processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
            genomes = genomes.groupby(genomes.index.str[1:8]).sum()
            index = genomes.index
            processAvg = np.array(processAvg)
        
        if processAvg.shape[0]==288 and collapse_to_SBS96==True: #collapse the 288 context into 96 only for the deocmposition 
            processAvg = pd.DataFrame(processAvg, index=index)
            processAvg = processAvg.groupby(processAvg.index.str[2:9]).sum()
            genomes = pd.DataFrame(genomes, index=index)
            genomes = genomes.groupby(genomes.index.str[2:9]).sum()
            index = genomes.index
            processAvg = np.array(processAvg)
                
        print("\n Decomposing De Novo Signatures  .....")
        final_signatures = sub.signature_decomposition(processAvg, mutation_type, layer_directory2, genome_build=genome_build,signature_database=signature_database, mutation_context=mutation_context, add_penalty=0.05, connected_sigs=connected_sigs,remove_penalty=0.01, make_decomposition_plots=make_plots, originalProcessAvg=originalProcessAvg,new_signature_thresh_hold=new_signature_thresh_hold,sig_exclusion_list=sig_exclusion_list)    
        #final_signatures = sub.signature_decomposition(processAvg, m, layer_directory2, genome_build=genome_build)
        # extract the global signatures and new signatures from the final_signatures dictionary
        globalsigs = final_signatures["globalsigs"]
        globalsigs = np.array(globalsigs)
        newsigs = final_signatures["newsigs"]
        processAvg = np.hstack([globalsigs, newsigs])  
        allsigids = final_signatures["globalsigids"]+final_signatures["newsigids"]
        attribution = final_signatures["dictionary"]
        background_sigs= final_signatures["background_sigs"]
        index = genomes.index
        colnames = genomes.columns

        print("\n Assigning decomposed signature")
        result = sub.make_final_solution(processAvg, genomes, allsigids, layer_directory2, mutation_type, index, colnames, 
                                cosmic_sigs=True, attribution = attribution, denovo_exposureAvg  = exposureAvg ,  
                                background_sigs=background_sigs, verbose=verbose, genome_build=genome_build, 
                                add_penalty=nnls_add_penalty, remove_penalty=nnls_remove_penalty, 
                                initial_remove_penalty=initial_remove_penalty,connected_sigs=connected_sigs,
                                collapse_to_SBS96=collapse_to_SBS96,
                                refit_denovo_signatures=False,make_plots=make_plots)


                                                                #################
                                                                # Cosmic Fitting       
                                                                #################
        
    if cosmic_fit_option ==True:
        try:
            if not os.path.exists(layer_directory3):
                os.makedirs(layer_directory3)
        except: 
            print ("The {} folder could not be created".format("Assignment_Solution"))
        # if signatures == None:
        #     processAvg = sub.getProcessAvg(genomes, genome_build, "3.2")
        #     processAvg = processAvg.set_index('Type').rename_axis('MutationType')
        # else:
        #     try:
        #         processAvg = pd.read_csv(signatures,sep='\t', index_col=0)
        #     except:
        #         sys.exit("Something is wrong with the format of input signatures, Pass a text file of signatures in the format of COSMIC sig database")
        index = genomes.index
        colnames = genomes.columns
        if (genomes.sum()==0).sum() >0:
            print("Removing samples with zero TMB ...... ")
            genomes=genomes.loc[:, (genomes != 0).any(axis=0)]
            colnames = genomes.columns
        
        if genomes.shape[0]==1536 and collapse_to_SBS96==True: #collapse the 1596 context into 96 only for the deocmposition 
            # processAvg = pd.DataFrame(processAvg, index=index)
            # processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
            genomes = pd.DataFrame(genomes, index=index)
            genomes = genomes.groupby(genomes.index.str[1:8]).sum()
            index = genomes.index
            #processAvg = np.array(processAvg)
        
        if genomes.shape[0]==288 and collapse_to_SBS96==True: #collapse the 288 context into 96 only for the deocmposition 
            # processAvg = pd.DataFrame(processAvg, index=index)
            # processAvg = processAvg.groupby(processAvg.index.str[2:9]).sum()
            genomes = pd.DataFrame(genomes, index=index)
            genomes = genomes.groupby(genomes.index.str[2:9]).sum()
            index = genomes.index


        if signature_database==None:
            processAvg = sub.getProcessAvg(genomes, genome_build=genome_build, cosmic_version=cosmic_version)[0]
            #processAvg = processAvg.set_index('Type').rename_axis('MutationType')
        else:
            try:
                processAvg = pd.read_csv(signature_database,sep='\t', index_col=0)
            except:
                sys.exit("Something is wrong with the format of signature database, Pass a text file of signatures in the format of COSMIC sig database")
        


        if processAvg.shape[0]==1536 and collapse_to_SBS96==True: #collapse the 1596 context into 96 only for the deocmposition 
            processAvg = pd.DataFrame(processAvg, index=index)
            processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
            genomes = pd.DataFrame(genomes, index=index)
            genomes = genomes.groupby(genomes.index.str[1:8]).sum()
            index = genomes.index
            #processAvg = np.array(processAvg)
        
        if processAvg.shape[0]==288 and collapse_to_SBS96==True: #collapse the 288 context into 96 only for the deocmposition 
            processAvg = pd.DataFrame(processAvg, index=index)
            processAvg = processAvg.groupby(processAvg.index.str[2:9]).sum()
            genomes = pd.DataFrame(genomes, index=index)
            genomes = genomes.groupby(genomes.index.str[2:9]).sum()
            index = genomes.index
            #processAvg = np.array(processAvg)

        
        

        #processAvg is sigdatabase: remove sigs corresponding to exclusion rules.
        sig_exclusion_list= ['SBS'+items for items in sig_exclusion_list]
        if sig_exclusion_list:
            print("The following signatures are excluded: "+" ".join(str(item) for item in sig_exclusion_list))
        # # 
        processAvg.drop(sig_exclusion_list, axis=1, inplace=True,errors='ignore')

        # import pdb
        # pdb.set_trace()        
        #processAvg= originalProcessAvg
        #index = genomes.index
        #colnames = genomes.columns
        allsigids = processAvg.columns.to_list()
        processAvg = processAvg.values
        attribution={}
        for i in allsigids:
            attribution[i]= [i] 
        #only for SBS96
        # pdb.set_trace()
        if mutation_type == "96" or mutation_type=="288" or mutation_type=="1536":        
            background_sigs = sub.get_indeces(list(allsigids), ['SBS1', 'SBS5'])
            # add connected signatures   
            #different_signatures = ss.add_connected_sigs(different_signatures, list(signames))
        #for other contexts
        else:
            background_sigs = []
        exposureAvg_dummy = pd.DataFrame(np.random.rand(processAvg.shape[1],genomes.shape[1]),index=allsigids,columns=colnames.to_list()).transpose().rename_axis('Samples')
        print("Assigning COSMIC sigs or Signature Database ...... ")
       
        if processAvg.shape[0] != 96:
            if genomes.shape[0] == processAvg.shape[0] and collapse_to_SBS96 ==True:
                sys.exit("Signatures Database and Samples are of same context type and is not equal to 96. please rerun by setting the flag \"collapse_to_SBS96 = False \"")
        
        # import pdb
        # pdb.set_trace()


        sub.make_final_solution(processAvg, genomes, allsigids, layer_directory3, mutation_type, index, colnames, 
                            cosmic_sigs=True, attribution = attribution, denovo_exposureAvg  = exposureAvg_dummy ,  
                            background_sigs=background_sigs, verbose=verbose, genome_build=genome_build, 
                            add_penalty=nnls_add_penalty, remove_penalty=nnls_remove_penalty, 
                            initial_remove_penalty=initial_remove_penalty,connected_sigs=connected_sigs,
                            collapse_to_SBS96=collapse_to_SBS96,refit_denovo_signatures=False, make_plots =make_plots)
   
  


