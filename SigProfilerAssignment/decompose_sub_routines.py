#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 16 2021

@author: rvangara
"""
import string 
import numpy as np
import os,sys
#from matplot,pdblib.backends.backend_pdf import PdfPages
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.backends.backend_pdf import PdfPages
#from sklearn import metrics
#import time
#import multiprocessing
#from multiprocessing import current_process
# from functools import partial
# from numpy import linalg as LA
import sigProfilerPlotting as plot
from SigProfilerExtractor import PlotDecomposition as sp
from SigProfilerExtractor import plotActivity as plot_ac
from SigProfilerExtractor import tmbplot as tmb
import string 
import PyPDF2
import scipy
#import SigProfilerAssignment as sspro
from PyPDF2 import PdfFileMerger
import SigProfilerAssignment as spa
from SigProfilerAssignment import single_sample as ss
from scipy.spatial.distance import correlation as cor
from alive_progress import alive_bar

def getProcessAvg(samples, genome_build="GRCh37", cosmic_version=3.2,signature_database=None,connected_sigs = True):
    paths = spa.__path__[0]
    
    if samples.shape[0]==96:
        sigDatabase = pd.read_csv(paths+"/data/Reference_Signatures/"+genome_build+"/COSMIC_v"+str(cosmic_version)+"_SBS_"+genome_build+".txt", sep="\t", index_col=0)
        signames = sigDatabase.columns   
        
    elif samples.shape[0]==288:
        sigDatabase = pd.read_csv(paths+"/data/Reference_Signatures/GRCh37/COSMIC_v"+str(cosmic_version)+"_SBS"+str(samples.shape[0])+"_GRCh37.txt", sep="\t", index_col=0)
        signames = sigDatabase.columns
        
    elif samples.shape[0]==1536:
        sigDatabase = pd.read_csv(paths+"/data/Reference_Signatures/"+"GRCh37"+"/COSMIC_v"+str(cosmic_version)+"_SBS"+str(samples.shape[0])+"_GRCh37.txt", sep="\t", index_col=0)
        signames = sigDatabase.columns
    
    elif samples.shape[0]==78:
        sigDatabase = pd.read_csv(paths+"/data/Reference_Signatures/"+genome_build+"/COSMIC_v"+str(cosmic_version)+"_DBS_"+genome_build+".txt", sep="\t", index_col=0)
        signames = sigDatabase.columns
        connected_sigs=False
        
    elif samples.shape[0]==83:
        sigDatabase = pd.read_csv(paths+"/data/Reference_Signatures/GRCh37/COSMIC_v"+str(cosmic_version)+"_ID_GRCh37.txt", sep="\t", index_col=0)
        signames = sigDatabase.columns
        connected_sigs=False
        
    elif samples.shape[0]==48:
        if cosmic_version < 3.3:
            print("The selected cosmic version is"+str(cosmic_version)+". CN signatures are available only for version 3.3. So, the cosmic version is reset to v3.3.")
            cosmic_version=3.3
        sigDatabase = pd.read_csv(paths+"/data/Reference_Signatures/GRCh37/COSMIC_v"+str(cosmic_version)+"_CN_GRCh37.txt", sep="\t",index_col=0)
        signames = sigDatabase.columns
        connected_sigs=False
    else:
        sigDatabase = pd.DataFrame(samples)
        sigDatabase.columns=sigDatabase.columns.astype(str)
        sigDatabase.index=sigDatabase.index.astype(str)
        signames=sigDatabase.columns
        connected_sigs=False
    return sigDatabase,signames,connected_sigs
    
    if signature_database != None:#pd.core.frame.DataFrame:
        print("################## USING CUSTOM SIGNATURE DATBASE ##################")
        signature_database= pd.read_csv(signature_database,sep="\t",index_col=0)
        if samples.shape[0]==signature_database.shape[0]:
            sigDatabase=signature_database
            signames = sigDatabase.columns 
            #make_decomposition_plots=False
            del signature_database
        else:
            sys.exit("The Signatures and the custom signature database have different context types.")
    sigDatabases = sigDatabase.reset_index()
    return sigDatabases


def signature_plotting_text(value, text, Type):
    name = text + ": "
    name_list =[]
    total=np.sum(np.array(value))
    for i in value:
        
        if Type=="integer":  
            i = int(i)
            p=round(i/total*100,1)
            i = format(i, ',d')   
            tail = str(i)+'/'+str(p)+'%'
            name_list.append(name+tail)
        elif Type=="float":
            i = round(i,2)
            name_list.append(name + str(i))
    return(name_list)

def make_letter_ids(idlenth = 10, mtype = "SBS96"):

    listOfSignatures = []
    letters = list(string.ascii_uppercase)
    letters.extend([i+b for i in letters for b in letters])
    letters = letters[0:idlenth]
    
    for j,l in zip(range(idlenth),letters):
        listOfSignatures.append(mtype+l)
    listOfSignatures = np.array(listOfSignatures)
    return listOfSignatures

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

def get_indeces(a, b):
    
    """ 
    Extracts the indices multiple items in a list.
    
    Parameters:
        a: list. where we want to get the index of the items.
        b: list. the items we want to get index of. 
    #example: 
    x = ['SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS8', 'SBS13', 'SBS40']
    y = ['SBS1',  'SBS5']
    get_indeces(x, y)
    #result
    >>> [1,3]
    """

    indeces = []
    for i in b:
        try: 
            idx = a.index(i)
            indeces.append(idx)
        except: 
            next

    return indeces

def get_items_from_index(x,y):
    """ decipher the values of items in a list from their indices.
    """
    z = []
    for i in y:
        try:
            z.append(x[i])
        except:
            pass
    return z

def signature_decomposition(signatures, mtype, directory, genome_build="GRCh37", cosmic_version=3.2,signature_database=None, add_penalty=0.05, remove_penalty=0.01, mutation_context=None, connected_sigs=True, make_decomposition_plots=True, originalProcessAvg=None,new_signature_thresh_hold=0.8,sig_exclusion_list=[]):

    originalProcessAvg = originalProcessAvg.reset_index()
    if not os.path.exists(directory+"/Solution_Stats"):
        os.makedirs(directory+"/Solution_Stats")
    # open the log file for signature decomposition 
    lognote = open(directory+"/Solution_Stats/Cosmic_"+mutation_context+"_Decomposition_Log.txt", "w") 
    lognote.write("############################ Signature Decomposition Details ################################\n\n\n")
    lognote.write("Context Type: {}\n".format(mtype))
    lognote.write("Genome Build: {}\n".format(genome_build))
    
    # paths = spa.__path__[0]
    
    # if signatures.shape[0]==96:
    #     sigDatabase = pd.read_csv(paths+"/data/Reference_Signatures/"+genome_build+"/COSMIC_v"+str(cosmic_version)+"_SBS_"+genome_build+".txt", sep="\t", index_col=0)
    #     signames = sigDatabase.columns   
        
    # elif signatures.shape[0]==288:
    #     sigDatabase = pd.read_csv(paths+"/data/Reference_Signatures/GRCh37/COSMIC_v"+str(3.2)+"_SBS"+str(signatures.shape[0])+"_GRCh37.txt", sep="\t", index_col=0)
    #     signames = sigDatabase.columns
        
    # elif signatures.shape[0]==1536:
    #     sigDatabase = pd.read_csv(paths+"/data/Reference_Signatures/"+"GRCh37"+"/COSMIC_v"+str(3.2)+"_SBS"+str(signatures.shape[0])+"_GRCh37.txt", sep="\t", index_col=0)
    #     signames = sigDatabase.columns
    
    # elif signatures.shape[0]==78:
    #     sigDatabase = pd.read_csv(paths+"/data/Reference_Signatures/"+"GRCh37"+"/COSMIC_v"+str(cosmic_version)+"_DBS_"+"GRCh37"+".txt", sep="\t", index_col=0)
    #     signames = sigDatabase.columns
    #     connected_sigs=False
        
    # elif signatures.shape[0]==83:
    #     sigDatabase = pd.read_csv(paths+"/data/Reference_Signatures/GRCh37/COSMIC_v"+str(cosmic_version)+"_ID_GRCh37.txt", sep="\t", index_col=0)
    #     signames = sigDatabase.columns
    #     connected_sigs=False
        
    # elif signatures.shape[0]==48:
    #     sigDatabase = pd.read_csv(paths+"/data/CNV_signatures.txt", sep="\t",index_col=0)
    #     signames = sigDatabase.columns
    #     connected_sigs=False
    # else:
    #     sigDatabase = pd.DataFrame(signatures)
    #     sigDatabase.columns=sigDatabase.columns.astype(str)
    #     sigDatabase.index=sigDatabase.index.astype(str)
    #     signames=sigDatabase.columns
    #     connected_sigs=False


    if signature_database==None:
        sigDatabase,signames,connected_sigs = getProcessAvg(signatures, genome_build, cosmic_version=cosmic_version)
        #processAvg = processAvg.set_index('Type').rename_axis('MutationType')
    else:
        try:
            sigDatabase = pd.read_csv(signature_database,sep='\t', index_col=0)
            #indx = sigDatabase.index()
            if sigDatabase.shape[0]==1536: #collapse the 1596 context into 96 only for the deocmposition 
                sigDatabase = sigDatabase.groupby(sigDatabase.index.str[1:8]).sum()
                
            elif sigDatabase.shape[0]==288 : #collapse the 288 context into 96 only for the deocmposition 
                #sigDatabase = pd.DataFrame(processAvg, index=index)
                sigDatabase = sigDatabase.groupby(sigDatabase.index.str[2:9]).sum()
            
            if sigDatabase.shape[0]== 78 or sigDatabase.shape[0]== 83 or sigDatabase.shape[0]== 48:
                connected_sigs=False
            lognote.write("##### Using a custom signature database for decomposition #####")
        except:
            sys.exit("Wrong format of signature database for decompose_fit, Please pass a text file of signatures in the format of COSMIC sig database")
    
    sig_exclusion_list= ['SBS'+items for items in sig_exclusion_list]
    lognote.write("The following signatures are excluded: "+" ".join(str(item) for item in sig_exclusion_list))
    sigDatabase.drop(sig_exclusion_list, axis=1, inplace=True,errors='ignore')
    signames=sigDatabase.columns
    
    # if type(signature_database)==pd.core.frame.DataFrame:
        
    #     if signatures.shape[0]==signature_database.shape[0]:
    #         sigDatabase=signature_database
    #         signames = sigDatabase.columns 
    #         #make_decomposition_plots=False
    #         del signature_database    
    sigDatabases = sigDatabase.reset_index()
    letters = list(string.ascii_uppercase)
    letters.extend([i+b for i in letters for b in letters])
    letters = letters[0:signatures.shape[1]]
    
    # replace the probability data of the process matrix with the number of mutation
    for i in range(signatures.shape[1]):
        signatures[:, i] =  signatures[:, i]*5000      #(np.sum(exposureAvg[i, :]))
    
    sigDatabase = np.array(sigDatabase)
    allsignatures = np.array([])
    newsig = list() # will create the letter based id of newsignatures
    newsigmatrixidx = list() # will create the original id of newsignature to help to record the original matrix
    fh = open(directory+"/De_Novo_map_to_COSMIC_"+mutation_context+".csv", "w")
    fh.write("De novo extracted, Global NMF Signatures, L1 Error %, L2 Error %, KL Divergence, Cosine Similarity, Correlation\n")
    fh.close()
    dictionary = {}
    # bgsigs = (np.where(np.isin(signames.tolist(),['SBS1','SBS5']))[0]).tolist()

        #only for SBS96
    if mtype == "96" or mtype=="288" or mtype=="1536":        
        bgsigs = get_indeces(list(signames), ['SBS1', 'SBS5'])
    else:
        bgsigs = []
    
    # import pdb
    # pdb.set_trace()

    # get the names of denovo signatures
    denovo_signature_names = make_letter_ids(signatures.shape[1], mtype=mutation_context)
    #lognote.write("\n********** Starting Signature Decomposition **********\n\n")
    activity_percentages=[]
    merger = PdfFileMerger()

    
    for i, j in zip(range(signatures.shape[1]), denovo_signature_names):
        
        # Only for context SBS96
        if signatures.shape[0]==96:
            lognote = open(directory+"/Solution_Stats/Cosmic_"+mutation_context+"_Decomposition_Log.txt", "a")  
            lognote.write("\n\n\n\n######################## Decomposing "+j+" ########################\n"  )
            lognote.close()
            if genome_build=="mm9" or genome_build=="mm10":
                check_rule_negatives = [1,16]
                check_rule_penalty=1.50
            else:
                check_rule_negatives = []
                check_rule_penalty=1.0
            
            # import pdb
            # pdb.set_trace()
            
            _, exposures,L2dist,similarity, kldiv, correlation, cosine_similarity_with_four_signatures = ss.add_remove_signatures(sigDatabase, 
                                                                                                         signatures[:,i], 
                                                                                                         metric="l2", 
                                                                                                         solver="nnls", 
                                                                                                         background_sigs = bgsigs,#[0,4], 
                                                                                                         permanent_sigs = bgsigs,#[0,4], 
                                                                                                         candidate_sigs="all", 
                                                                                                         allsigids = signames, 
                                                                                                         add_penalty = add_penalty, 
                                                                                                         remove_penalty = remove_penalty,
                                                                                                         check_rule_negatives = check_rule_negatives, 
                                                                                                         checkrule_penalty = check_rule_penalty, 
                                                                                                         directory = directory+"/Solution_Stats/Cosmic_"+mutation_context+"_Decomposition_Log.txt", 
                                                                                                         connected_sigs=connected_sigs,
                                                                                                         verbose=False)
    
        else:
            lognote = open(directory+"/Solution_Stats/Cosmic_"+mutation_context+"_Decomposition_Log.txt", "a")  
            lognote.write("\n\n\n\n######################## Decomposing "+j+" ########################\n"  )
            lognote.close()
            
            _, exposures,L2dist,similarity, kldiv, correlation, cosine_similarity_with_four_signatures = ss.add_remove_signatures(sigDatabase, 
                                                                                                         signatures[:,i], 
                                                                                                         metric="l2", 
                                                                                                         solver="nnls", 
                                                                                                         background_sigs = [], 
                                                                                                         candidate_sigs="all", 
                                                                                                         add_penalty = add_penalty, 
                                                                                                         remove_penalty = remove_penalty,
                                                                                                         check_rule_negatives = [], 
                                                                                                         checkrule_penalty = [], 
                                                                                                         directory = directory+"/Solution_Stats/Cosmic_"+mutation_context+"_Decomposition_Log.txt", 
                                                                                                         connected_sigs=connected_sigs,
                                                                                                         verbose=False)
        L1dist = np.linalg.norm(signatures[:,i]-np.dot(sigDatabase,exposures) , ord=1)/np.linalg.norm(signatures[:,i], ord=1)
        exposure_percentages = exposures[np.nonzero(exposures)]/np.sum(exposures[np.nonzero(exposures)])*100
        listofinformation = list("0"*len(np.nonzero(exposures)[0])*3)
        
        count =0
        decomposed_signatures = []
        contribution_percentages = []
        
        for j in np.nonzero(exposures)[0]:
            listofinformation[count*3] = signames[j]
            listofinformation[count*3+1] = round(exposure_percentages[count],2)
            contribution_percentages.append(round(exposure_percentages[count],2))
            listofinformation[count*3+2]="%"
            decomposed_signatures.append(signames[j])
            count+=1
        ListToTumple = tuple([mtype, letters[i]]+listofinformation+[L1dist*100]+[L2dist*100]+[kldiv]+[similarity]+[correlation])
        activity_percentages.append(contribution_percentages)
        
        weights=[]
        basis_names=[]
        nonzero_exposures=exposures[np.nonzero(exposures)]
        denovo_name=mutation_context+letters[i]
        for info in range(0, len(listofinformation), 3):
            #print(info)
            sigName=listofinformation[info]
            sigWeigt=str(listofinformation[info+1])+"%"
            weights.append(sigWeigt)
            basis_names.append(sigName)
        
        denovo_signames=[]
        for letter in letters:
            denovo_signames.append(mutation_context+letter)
       
        
        sigDatabases_DF=sigDatabases
        
        if mtype=="1536":
            mtype_par="1536"
        elif mtype=="288":
            mtype_par="288"
        elif mtype=="96":
            mtype_par="96"
        elif mtype=="DINUC" or mtype=="78":
            mtype_par="78"
        elif mtype=="INDEL" or mtype=="83":
            mtype_par="83"
        elif mtype=="CNV" or mtype=="48":
            mtype_par="48"
        else:
            mtype_par="none"
        try:
            # import pdb
            # pdb.set_trace() 
            if mtype_par!="none" and make_decomposition_plots==True:
                # Get the names of the columns for each dataframe
                denovo_col_names = originalProcessAvg.columns
                cosmic_col_names = sigDatabases_DF.columns
                # Get the name for the MutationTypes column
                cosmic_mut_types_col = cosmic_col_names[0]
                denovo_mut_types_col =  denovo_col_names[0]
                # create lists of implemented columns
                basis_cols = basis_names.copy()
                basis_cols.insert(0,cosmic_mut_types_col)
                denovo_cols=[denovo_mut_types_col, denovo_name]
                byte_plot = sp.run_PlotDecomposition(originalProcessAvg[denovo_cols], denovo_name, sigDatabases_DF[basis_cols], basis_names, weights, nonzero_exposures/5000, directory, "test", mtype_par)
                merger.append(byte_plot)
                with alive_bar(1, ctrl_c=False,bar='blocks', title=f'Decompositon Plot:{denovo_name}') as bar:
                    bar()
                #print("Decompositon Plot made for {}".format(denovo_name))
        except:
            print("The context-" + str(mtype_par) + " decomposition plots pages were not able to be generated.")
        
        strings ="Signature %s-%s,"+" Signature %s (%0.2f%s) &"*(len(np.nonzero(exposures)[0])-1)+" Signature %s (%0.2f%s), %0.2f,  %0.2f, %0.3f, %0.2f, %0.2f\n" 
        #new_signature_thresh_hold = 0.8
        if  similarity>new_signature_thresh_hold and cosine_similarity_with_four_signatures > new_signature_thresh_hold: ########### minimum signtatures and cosine similarity needs to be fitted to become a unique signature 
            allsignatures = np.append(allsignatures, np.nonzero(exposures))
            fh = open(directory+"/De_Novo_map_to_COSMIC_"+mutation_context+".csv", "a")
            fh.write(strings%(ListToTumple))
            fh.close()
            
            dictionary.update({"{}".format(mutation_context+letters[i]):decomposed_signatures}) 
            
        else:
            newsig.append(mutation_context+letters[i])
            newsigmatrixidx.append(i)
            fh = open(directory+"/De_Novo_map_to_COSMIC_"+mutation_context+".csv", "a")
            fh.write("Signature {}-{}, Signature {}-{}, {}, {}, {}, {}, {}\n".format(mtype, letters[i], mtype, letters[i], 0, 0, 0, 1, 1))
            fh.close()
            dictionary.update({"{}".format(mutation_context+letters[i]):["{}".format(mutation_context+letters[i])]}) 
            #dictionary.update({letters[i]:"Signature {}-{}, Signature {}-{}, {}\n".format(mtype, letters[i], mtype, letters[i], 1 )}) 
    
    try:
        if make_decomposition_plots and mtype_par != 'none':
            # Write out the decomposition plots   
            contexts = {'96':'SBS96', '288':'SBS288', '1536':'SBS1536', '78':'DBS78', '83':'ID83', "48":"CNV"}
            merger.write(directory+"/"+contexts[mtype_par]+"_Decomposition_Plots.pdf")
    except:
        print("The context-" + str(mtype_par) + " decomposition pages were not able to be merged.")
    
    different_signatures = np.unique(allsignatures)
    different_signatures=different_signatures.astype(int)
    if mtype == "96" or mtype=="288" or mtype=="1536":
        different_signatures = list(set().union(different_signatures, bgsigs))
        different_signatures.sort()    
      
    
    #get the name of the signatures
    try:
        detected_signatures = signames[different_signatures]
        globalsigmats= sigDatabases.loc[:,list(detected_signatures)]
    except:
        detected_signatures=[None]
        globalsigmats=None
    
    newsigsmats=signatures[:,newsigmatrixidx]
        
    #only for SBS96
    if mtype == "96" or mtype=="288" or mtype=="1536":        
        background_sigs = get_indeces(list(detected_signatures), ['SBS1', 'SBS5'])
        # add connected signatures   
        different_signatures = ss.add_connected_sigs(different_signatures, list(signames))
    #for other contexts
    else:
        background_sigs = []
        
    # close the lognote
    lognote.close()

    # import pdb
    # pdb.set_trace()    
    #return values
    return {"globalsigids": list(detected_signatures), "newsigids": newsig, "globalsigs":globalsigmats, "newsigs":newsigsmats/5000, "dictionary": dictionary, 
            "background_sigs": background_sigs, "activity_percentages": activity_percentages} 


############################################################################################################
######################################## MAKE THE FINAL FOLDER ##############################################
#############################################################################################################
def make_final_solution(processAvg, allgenomes, allsigids, layer_directory, m, index, allcolnames, process_std_error = "none", signature_stabilities = " ", \
                        signature_total_mutations= " ", signature_stats = "none",  cosmic_sigs=False, attribution= 0, denovo_exposureAvg  = "none", add_penalty=0.05, \
                        remove_penalty=0.01, initial_remove_penalty=0.05, de_novo_fit_penalty=0.02, background_sigs=0, genome_build="GRCh37", sequence="genome", export_probabilities=True, \
                        refit_denovo_signatures=True, collapse_to_SBS96=True, connected_sigs=True, pcawg_rule=False, verbose=False,make_plots = True):

    if processAvg.shape[0]==allgenomes.shape[0] and processAvg.shape[0] != 96:
        collapse_to_SBS96=False


    # Get the type of solution from the last part of the layer_directory name
    solution_type = layer_directory.split("/")[-1]
    solution_prefix = solution_type.split("_")
    solution_prefix = "_".join(solution_prefix[0:2])
    if refit_denovo_signatures==True:
            solution_prefix_refit=solution_prefix+"_refit"
    
    if not os.path.exists(layer_directory+"/Signatures"):
        os.makedirs(layer_directory+"/Signatures")
    if not os.path.exists(layer_directory+"/Activities"):
        os.makedirs(layer_directory+"/Activities")
    if not os.path.exists(layer_directory+"/Solution_Stats"):
        os.makedirs(layer_directory+"/Solution_Stats")
          
    
    # Create the lognote file
    if refit_denovo_signatures==True:
        lognote = open(layer_directory+"/Solution_Stats/"+solution_prefix_refit+"_Signature_Assignment_log.txt", "w")
    else:
        lognote = open(layer_directory+"/Solution_Stats/"+solution_prefix+"_Signature_Assignment_log.txt", "w")
    lognote.write("************************ Stepwise Description of Signature Assignment to Samples ************************")
    lognote.close()
    
    
    #Get the type of Signatures
    if m == 83 or m=="83":
        signature_type = "INDEL83"
        connected_sigs=False
    elif m==78 or m=="78":
        signature_type = "DINUC78"
        connected_sigs=False
    else:
        signature_type = "SBS"+str(m)
    
    
    allgenomes = np.array(allgenomes)
    if (m=="96" or m=="1536" or m=="288") and (genome_build=="mm9" or genome_build=="mm10") and (collapse_to_SBS96==True):
        check_rule_negatives = [1,16]
        check_rule_penalty=1.50
    else:
        check_rule_negatives = []
        check_rule_penalty=1.0
    exposureAvg = np.zeros([processAvg.shape[1], allgenomes.shape[1]] )  
    if cosmic_sigs==True:
        denovo_exposureAvg = denovo_exposureAvg.T
        with alive_bar(allgenomes.shape[1]) as bar:
        #print("\n")
            for r in range(allgenomes.shape[1]):
                #print("Analyzing Sample => " , str(r+1))
                bar()
                if verbose==True:
                    print("\n\n\n\n\n                                        ################ Sample "+str(r+1)+ " #################")
                        
                # Record information to lognote
                lognote = open(layer_directory+"/Solution_Stats/"+solution_prefix+"_Signature_Assignment_log.txt", "a")
                lognote.write("\n\n\n\n\n                    ################ Sample "+str(r+1)+ " #################\n") 
                
                sample_exposure = np.array(denovo_exposureAvg.iloc[:,r])
                
                init_sig_idx = np.nonzero(sample_exposure)[0]
                init_sigs = denovo_exposureAvg.index[init_sig_idx]
                
                
                init_decomposed_sigs = []
                for de_novo_sig in init_sigs:
                    
                    init_decomposed_sigs = union(init_decomposed_sigs, list(attribution[de_novo_sig]))
                    
                #print(init_decomposed_sigs) 
                init_decomposed_sigs_idx = get_indeces(allsigids, init_decomposed_sigs)
                init_decomposed_sigs_idx.sort()
                init_decomposed_sigs_idx = list(set().union(init_decomposed_sigs_idx, background_sigs))
                #print(init_decomposed_sigs_idx)
                
                # get the indices of the background sigs in the initial signatures
                background_sig_idx = get_indeces(init_decomposed_sigs_idx, background_sigs)
                
                
                fit_signatures = processAvg[:,init_decomposed_sigs_idx]
                #fit signatures
                newExposure, newSimilarity = ss.fit_signatures(fit_signatures, allgenomes[:,r])
                
                
                #create the exposureAvg vector
                #print(init_decomposed_sigs_idx)
                #print(newExposure)
                for nonzero_idx, nozero_exp in zip(init_decomposed_sigs_idx, newExposure):
                    exposureAvg[nonzero_idx, r] = nozero_exp
                
                
                if pcawg_rule==True:
                    maxmutation=np.sum(allgenomes[:,r])
                    exposureAvg[:, r], remove_distance, _ = ss.remove_all_single_signatures(processAvg, exposureAvg[:, r], allgenomes[:,r], metric="l2", verbose = False, cutoff=0.02)
                    # get the maximum value of the new Exposure
                    maxcoef = max(list(exposureAvg[:, r]))
                    idxmaxcoef = list(exposureAvg[:, r]).index(maxcoef)
                
                    exposureAvg[:, r] = np.round(exposureAvg[:, r])
                
                    # We may need to tweak the maximum value of the new exposure to keep the total number of mutation equal to the original mutations in a genome
                    if np.sum(exposureAvg[:, r])!=maxmutation:
                        exposureAvg[:, r][idxmaxcoef] = round(exposureAvg[:, r][idxmaxcoef])+maxmutation-sum(exposureAvg[:, r])
                    #print(exposureAvg[:, r]) 
                    #print("\n")
                        
                else:
                    if verbose==True:
                        print("############################# Initial Composition #################################### ") 
                        print(pd.DataFrame(exposureAvg[:, r],  index=allsigids).T)   
                        print("L2%: ", newSimilarity)  
                        
                    lognote.write("############################# Initial Composition ####################################\n")
                    exposures = pd.DataFrame(exposureAvg[:, r],  index=allsigids).T
                    lognote.write("{}\n".format(exposures.iloc[:,exposures.to_numpy().nonzero()[1]])) 
                    lognote.write("L2 Error %: {}\nCosine Similarity: {}\n".format(round(newSimilarity,2), round(cos_sim(allgenomes[:,r], np.dot(processAvg, exposureAvg[:, r] )),2)))
                    #remove signatures 
                    exposureAvg[:,r],L2dist,cosine_sim = ss.remove_all_single_signatures(processAvg, exposureAvg[:, r], allgenomes[:,r], metric="l2", \
                            solver = "nnls", cutoff=initial_remove_penalty, background_sigs= [], verbose=False)
                    if verbose==True:
                        print("############################## Composition After Initial Remove ############################### ")
                        print(pd.DataFrame(exposureAvg[:, r],  index=allsigids).T)  
                        print("L2%: ", L2dist)
                    lognote.write("############################## Composition After Initial Remove ###############################\n")
                    exposures = pd.DataFrame(exposureAvg[:, r],  index=allsigids).T
                    lognote.write("{}\n".format(exposures.iloc[:,exposures.to_numpy().nonzero()[1]])) 
                    lognote.write("L2 Error %: {}\nCosine Similarity: {}\n".format(round(L2dist,2), round(cosine_sim,2)))
                    lognote.write("\n############################## Performing Add-Remove Step ##############################\n")
                    #Close the Lognote file
                    lognote.close()
                    
                    init_add_sig_idx = list(set().union(list(np.nonzero(exposureAvg[:, r])[0]), background_sigs))
                    #print(init_add_sig_idx)
                    
                    
                    #get the background_sig_idx for the add_remove function only for the decomposed solution:
                    if background_sigs != 0:  # in the decomposed solution only 
                        background_sig_idx = get_indeces(allsigids, ["SBS1", "SBS5"])
                
                    # if the there is no other signatures to be added on top the existing signatures
                    try:
                        _, exposureAvg[:, r],L2dist,similarity, kldiv, correlation, cosine_similarity_with_four_signatures = ss.add_remove_signatures(processAvg, 
                                                                                                        allgenomes[:,r], 
                                                                                                        metric="l2", 
                                                                                                        solver="nnls", 
                                                                                                        background_sigs = init_add_sig_idx, 
                                                                                                        permanent_sigs = background_sig_idx, 
                                                                                                        candidate_sigs="all", 
                                                                                                        allsigids = allsigids, 
                                                                                                        add_penalty = add_penalty, 
                                                                                                        remove_penalty=remove_penalty,
                                                                                                        check_rule_negatives = check_rule_negatives, 
                                                                                                        checkrule_penalty = check_rule_penalty, 
                                                                                                        connected_sigs=connected_sigs,
                                                                                                        directory = layer_directory+"/Solution_Stats/"+solution_prefix+"_Signature_Assignment_log.txt", 
                                                                                                        verbose=False)
                        
                        if verbose==True:
                            print("####################################### Composition After Add-Remove #######################################\n") 
                            print(exposureAvg[:, r])
                            print("L2%: ", L2dist)
                        # Recond the information in the log file
                        lognote = open(layer_directory+"/Solution_Stats/"+solution_prefix+"_Signature_Assignment_log.txt", "a")
                        lognote.write("####################################### Composition After Add-Remove #######################################\n")
                        exposures = pd.DataFrame(exposureAvg[:, r],  index=allsigids).T
                        lognote.write("{}\n".format(exposures.iloc[:,exposures.to_numpy().nonzero()[1]])) 
                        lognote.write("L2 Error %: {}\nCosine Similarity: {}\n".format(round(L2dist,2), round(similarity,2)))
                        lognote.close()
                    except:
                        pass
            
    else:   
        # when refilt de_novo_signatures 
        refit_denovo_signatures_old = False
        if refit_denovo_signatures_old==True:
            exposureAvg=denovo_exposureAvg
            for g in range(allgenomes.shape[1]):
                print("Analyzing Sample => " , str(g+1))
                
                # Record information to lognote
                lognote = open(layer_directory+"/Solution_Stats/"+solution_prefix_refit+"_Signature_Assignment_log.txt", "a")
                lognote.write("\n\n\n\n\n                    ################ Sample "+str(g+1)+ " #################\n")

                lognote.write("############################# Initial Composition ####################################\n")
                exposures = pd.DataFrame(exposureAvg[:, g],  index=allsigids).T
                lognote.write("{}\n".format(exposures.iloc[:,exposures.to_numpy().nonzero()[1]])) 
                
                #remove signatures 
                exposureAvg[:,g],L2dist,cosine_sim = ss.remove_all_single_signatures(processAvg, exposureAvg[:, g], allgenomes[:,g], metric="l2", \
                           solver = "nnls", cutoff=de_novo_fit_penalty, background_sigs= [], verbose=False)
                if verbose==True:
                    print("############################## Composition After Remove ############################### ")
                    print(pd.DataFrame(exposureAvg[:, g],  index=allsigids).T)  
                    print("L2%: ", L2dist)
                lognote.write("############################## Composition After  Remove ###############################\n")
                exposures = pd.DataFrame(exposureAvg[:, g],  index=allsigids).T
                lognote.write("{}\n".format(exposures.iloc[:,exposures.to_numpy().nonzero()[1]])) 
                lognote.write("L2 Error %: {}\nCosine Similarity: {}\n".format(round(L2dist,2), round(cosine_sim,2)))
                lognote.close()
                
        # when use the exposures from the initial NMF
        else:
            exposureAvg=denovo_exposureAvg

    processAvg= pd.DataFrame(processAvg.astype(float))
    processes = processAvg.set_index(index)
    processes.columns = allsigids
    processes = processes.rename_axis("MutationsType", axis="columns")
    processes.to_csv(layer_directory+"/Signatures"+"/"+solution_prefix+"_"+"Signatures.txt", "\t", float_format='%.8f',index_label=[processes.columns.name]) 
    exposureAvg = pd.DataFrame(exposureAvg.astype(int))
    allsigids = np.array(allsigids)
    exposures = exposureAvg.set_index(allsigids)
    exposures.columns = allcolnames
    exposures = exposures.T
    exposures = exposures.rename_axis("Samples", axis="columns")
    if refit_denovo_signatures==True:
        exposures.to_csv(layer_directory+"/Activities"+"/"+solution_prefix+"_"+"Activities_refit.txt", "\t", index_label=[exposures.columns.name]) 
    else:
        exposures.to_csv(layer_directory+"/Activities"+"/"+solution_prefix+"_"+"Activities.txt", "\t", index_label=[exposures.columns.name]) 

    #plt tmb
    tmb_exposures = pd.melt(exposures)
    if make_plots ==True:
        if refit_denovo_signatures==True:
            tmb.plotTMB(tmb_exposures, scale=sequence, Yrange="adapt", output= layer_directory+"/Activities"+"/"+solution_prefix+"_"+"TMB_plot_refit.pdf")
        else:
            tmb.plotTMB(tmb_exposures, scale=sequence, Yrange="adapt", output= layer_directory+"/Activities"+"/"+solution_prefix+"_"+"TMB_plot.pdf")
        del tmb_exposures
        
    #plot activities
    if make_plots ==True:
        if refit_denovo_signatures==True:
            plot_ac.plotActivity(layer_directory+"/Activities"+"/"+solution_prefix+"_"+"Activities_refit.txt", output_file = layer_directory+"/Activities/"+solution_prefix+"_"+"Activity_Plots_refit.pdf", bin_size = 50, log = False)
        else:
            plot_ac.plotActivity(layer_directory+"/Activities"+"/"+solution_prefix+"_"+"Activities.txt", output_file = layer_directory+"/Activities/"+solution_prefix+"_"+"Activity_Plots.pdf", bin_size = 50, log = False)
    
    # Calcutlate the similarity matrices
    est_genomes = np.dot(processAvg, exposureAvg)
    all_similarities, cosine_similarities = calculate_similarities(allgenomes, est_genomes, allcolnames)
    all_similarities.iloc[:,[3,5]] = all_similarities.iloc[:,[3,5]].astype(str) + '%'
    
    if refit_denovo_signatures==True:
        all_similarities.to_csv(layer_directory+"/Solution_Stats/"+solution_prefix+"_Samples_Stats_refit.txt", sep="\t")
    else:
        all_similarities.to_csv(layer_directory+"/Solution_Stats/"+solution_prefix+"_Samples_Stats.txt", sep="\t")
    
    #if cosmic_sigs==False:
    if refit_denovo_signatures ==True:
        try:
            process_std_error= pd.DataFrame(process_std_error)
            processSTE = process_std_error.set_index(index)
            processSTE.columns = allsigids
            processSTE = processSTE.rename_axis("MutationType", axis="columns")
            processSTE.to_csv(layer_directory+"/Signatures"+"/"+solution_prefix+"_"+"Signatures_SEM_Error.txt", "\t", float_format='%.2E', index_label=[processes.columns.name]) 
        except:
            pass
    #if cosmic_sigs==False:
    if refit_denovo_signatures ==True:
        try: 
            signature_stats = signature_stats.set_index(allsigids)
            signature_stats = signature_stats.rename_axis("Signatures", axis="columns")
            signature_stats.to_csv(layer_directory+"/Solution_Stats"+"/"+solution_prefix+"_"+"Signatures_Stats.txt", "\t", index_label=[exposures.columns.name]) 
            signature_total_mutations = np.sum(exposureAvg, axis =1).astype(int)
            signature_total_mutations = signature_plotting_text(signature_total_mutations, "Sig. Mutations", "integer")
        except:
            pass
    else: #when it works with the decomposed solution
        signature_total_mutations = np.sum(exposureAvg, axis =1).astype(int)
        signature_total_mutations = signature_plotting_text(signature_total_mutations, "Sig. Mutations", "integer")
        if (m == "1536" or m=="288") and collapse_to_SBS96==True: # collapse the 1536 to 96
            m = "96"  
    if make_plots == True:
    ########################################### PLOT THE SIGNATURES ################################################
        if m=="DINUC" or m=="78":
            plot.plotDBS(layer_directory+"/Signatures/"+solution_prefix+"_"+"Signatures.txt", layer_directory+"/Signatures"+"/" , solution_prefix, "78", True, custom_text_upper= signature_stabilities, custom_text_middle = signature_total_mutations )        
        elif m=="INDEL" or m=="83":
            plot.plotID(layer_directory+"/Signatures/"+solution_prefix+"_"+"Signatures.txt", layer_directory+"/Signatures"+"/" , solution_prefix, "94", True, custom_text_upper= signature_stabilities, custom_text_middle = signature_total_mutations )
        elif m=="CNV" or m=="48":
            plot.plotCNV(layer_directory+"/Signatures/"+solution_prefix+"_"+"Signatures.txt", layer_directory+"/Signatures"+"/" , solution_prefix, "pdf", percentage=True, aggregate=False)
        elif m=="SV" or m=="32":
            plot.plotSV(layer_directory+"/Signatures/"+solution_prefix+"_"+"Signatures.txt", layer_directory+"/Signatures"+"/" , solution_prefix, "pdf", percentage=True, aggregate=False)
        elif (m=="96" or m=="288" or m=="384" or m=="1536") and collapse_to_SBS96==True:
            plot.plotSBS(layer_directory+"/Signatures/"+solution_prefix+"_"+"Signatures.txt", layer_directory+"/Signatures"+"/", solution_prefix, m, True, custom_text_upper= signature_stabilities, custom_text_middle = signature_total_mutations )
        elif m=="96":
            plot.plotSBS(layer_directory+"/Signatures/"+solution_prefix+"_"+"Signatures.txt", layer_directory+"/Signatures"+"/", solution_prefix, m, True, custom_text_upper= signature_stabilities, custom_text_middle = signature_total_mutations )
        elif m=="288":
            plot.plotSBS(layer_directory+"/Signatures/"+solution_prefix+"_"+"Signatures.txt", layer_directory+"/Signatures"+"/", solution_prefix, m, True, custom_text_upper= signature_stabilities, custom_text_middle = signature_total_mutations )
        elif m=="1536":
            plot.plotSBS(layer_directory+"/Signatures/"+solution_prefix+"_"+"Signatures.txt", layer_directory+"/Signatures"+"/", solution_prefix, m, True, custom_text_upper= signature_stabilities, custom_text_middle = signature_total_mutations )
        else:
            custom_signatures_plot(processes, layer_directory+"/Signatures")
    
    probability = probabilities(processAvg, exposureAvg, index, allsigids, allcolnames)
    probability=probability.set_index("Sample Names" )
    
    if cosmic_sigs==False:
        
        if refit_denovo_signatures==True:
            probability.to_csv(layer_directory+"/Activities"+"/"+"De_Novo_Mutation_Probabilities_refit.txt", "\t") 
        else:
            probability.to_csv(layer_directory+"/Activities"+"/"+"De_Novo_Mutation_Probabilities.txt", "\t") 
    if cosmic_sigs==True:
        probability.to_csv(layer_directory+"/Activities"+"/"+"Decomposed_Mutation_Probabilities.txt", "\t") 
    

    return exposures
################################################################### FUNCTION ONE ###################################################################
#function to calculate multiple similarities/distances
def calculate_similarities(genomes, est_genomes, sample_names=False):
    from numpy import inf
    
    if  sample_names is False:
        sample_names = ["None"]*genomes.shape[1]
        
    cosine_similarity_list = []
    kl_divergence_list = []
    correlation_list=[]
    l1_norm_list = []
    l2_norm_list = []
    total_mutations_list = []
    relative_l1_list = []
    relative_l2_list = []
    
    for i in range(genomes.shape[1]):
        p_i = genomes[:,i]
        q_i = est_genomes[:, i]
        cosine_similarity_list.append(round(cos_sim(p_i,q_i ),3))
        kl_divergence_list.append(round(scipy.stats.entropy(p_i,q_i),5))
        correlation_list.append(round(scipy.stats.pearsonr(p_i,q_i)[0],3))
        l1_norm_list.append(round(np.linalg.norm(p_i-q_i , ord=1),3))
        relative_l1_list.append(round((l1_norm_list[-1]/np.linalg.norm(p_i, ord=1))*100,3))
        l2_norm_list.append(round(np.linalg.norm(p_i-q_i , ord=2),3))
        relative_l2_list.append(round((l2_norm_list[-1]/np.linalg.norm(p_i, ord=2))*100,3))
        total_mutations_list.append(np.sum(p_i))

    kl_divergence_list = np.array(kl_divergence_list)
    kl_divergence_list[kl_divergence_list == inf] =1000
    similarities_dataframe = pd.DataFrame({"Sample Names": sample_names, \
                                           "Total Mutations":total_mutations_list, \
                                           "Cosine Similarity": cosine_similarity_list, \
                                           "L1 Norm": l1_norm_list, \
                                           "L1_Norm_%":relative_l1_list, \
                                           "L2 Norm": l2_norm_list, \
                                           "L2_Norm_%": relative_l2_list, \
                                           "KL Divergence": kl_divergence_list, \
                                           "Correlation": correlation_list })
    similarities_dataframe = similarities_dataframe.set_index("Sample Names")
    return [similarities_dataframe, cosine_similarity_list]

############################################################### FUNCTIONS TO CALCULATE DISTANCES BETWEEN VECTORS ##################################################
################################################################### FUNCTION ONE ###################################################################
#function to calculate the cosine similarity
def cos_sim(a, b):
      
    
    """Takes 2 vectors a, b and returns the cosine similarity according 
    to the definition of the dot product
    
    Dependencies: 
    *Requires numpy library. 
    *Does not require any custom function (constructed by me)
    
    Required by:
    * pairwise_cluster_raw
    	"""
    if np.sum(a)==0 or np.sum(b) == 0:
        return 0.0      
    dot_product = np.dot(a, b)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    return dot_product / (norm_a * norm_b)

def cor_sim(a, b):
      
    
    """Takes 2 vectors a, b and returns the corrilation similarity according 
    to the definition of the dot product
    
    Dependencies: 
    *Requires numpy library. 
    *Does not require any custom function (constructed by me)
    
    Required by:
    * pairwise_cluster_raw
    	"""
    if np.sum(a)==0 or np.sum(b) == 0:
        return 0.0      
    corr =1-cor(a, b)
    return corr

################################################### Generation of probabilities for each processes given to A mutation type ############################################
def probabilities(W, H, index, allsigids, allcolnames):  
    
    # setting up the indices 
    rows = index
    cols = allcolnames
    sigs = allsigids
    
    W = np.array(W)
    H= np.array(H)
    # rebuild the original matrix from the estimated W and H 
    genomes = np.dot(W,H)
    
    
    result = 0
    for i in range(H.shape[1]): #here H.shape is the number of sample
        
        M = genomes[:,i][np.newaxis]
        probs = W*H[:,i]/M.T        
        probs = pd.DataFrame(probs)
        probs.columns = sigs
        col1 = [cols[i]]*len(rows)
        probs.insert(loc=0, column='Sample Names', value=col1)
        probs.insert(loc=1, column='MutationTypes', value = rows)
        if i!=0:
            result = pd.concat([result, probs], axis=0)
        else:
            result = probs
    
        
    return result

def custom_signatures_plot(signatures, output):
    with PdfPages(output+'/Custom_Signature_Plots.pdf') as pdf:
        plt.figure(figsize=(10, 3))
        plt.bar(list(range(1,1+len(signatures.iloc[:,0]))),signatures.iloc[:,0])
        plt.title('Custom Signature {}'.format(0+1))
        plt.xticks([])
        plt.xlabel("Mutation Types")
        plt.ylabel("Probabilities")
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
        for i in range(1,signatures.shape[1]):
            # if LaTeX is not installed or error caught, change to `usetex=False`
            plt.rc('text', usetex=False)
            plt.figure(figsize=(10, 3))
            plt.bar(list(range(1, 1+len(signatures.iloc[:,i]))),signatures.iloc[:,i])
            plt.title('Custom Signature {}'.format(i+1))
            plt.xticks([])
            plt.xlabel("Mutation Types")
            plt.ylabel("Probabilities")
            pdf.attach_note("signature plots")  
            pdf.savefig()
            plt.close()

def merge_pdf(input_folder, output_file):
    pdf2merge = []
    for filename in os.listdir(input_folder):
        #print(filename)
        if filename.endswith('.pdf'):
            pdf2merge.append(filename)
            
    pdf2merge.sort()
    pdfWriter = PyPDF2.PdfFileWriter()
    for filename in pdf2merge:
        pdfFileObj = open(input_folder+"/"+filename,'rb')
        pdfReader = PyPDF2.PdfFileReader(pdfFileObj)
        for pageNum in range(pdfReader.numPages):
            pageObj = pdfReader.getPage(pageNum)
            pdfWriter.addPage(pageObj)
            
    pdfOutput = open(output_file+'.pdf', 'wb')
    pdfWriter.write(pdfOutput)
    #Outputting the PDF
    pdfOutput.close()
