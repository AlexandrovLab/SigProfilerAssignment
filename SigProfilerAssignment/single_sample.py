#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 10:55:51 2019

@author: mishugeb
"""
import numpy as np
from numpy import linalg as LA
from scipy.optimize import nnls
from scipy.optimize import minimize
from SigProfilerAssignment import decompose_subroutines as sub
import scipy.stats
import copy
import os


#############################################################################################################
#################################### Functions For Single Sample Algorithms #############################

################################################################### FUNCTION ONE ###################################################################
# function to calculate the cosine similarity
def cos_sim(a, b):

    """Takes 2 vectors a, b and returns the cosine similarity according
    to the definition of the dot product

    Dependencies:
    *Requires numpy library.
    *Does not require any custom function (constructed by me)

    Required by:
    * pairwise_cluster_raw
    """
    if np.sum(a) == 0 or np.sum(b) == 0:
        return 0.0
    dot_product = np.dot(a, b)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    return dot_product / (norm_a * norm_b)


def re_distribute_signature_activity(newExposure, genome):
    act = np.asarray(newExposure)
    print("maxmutation before rounding:", genome)
    maxmutation = round(np.sum(genome))
    print("maxmutation after rounding:", maxmutation)
    act2 = (act / act.sum()) * maxmutation
    return np.asarray(act2)


# def weighted_threshold(normalised_solution, genome, threshold):
# normalised_solution[normalised_solution < threshold] = 0
# threshold_solution = re_distribute_signature_activity(normalised_solution, genome)
# return threshold_solution


def weighted_threshold_split(
    normalised_solution, genome, threshold, ALLSIGIDSGLOBAL=[]
):
    # import pdb; pdb.set_trace();
    # print("this is the normalised solution:"+str(normalised_solution))
    # print("this is the normalised solution len2:"+str(normalised_solution.shape[0]))
    # print(type(normalised_solution))
    clock_threshold = -1  # make sure clock like signatures will not be removed
    if len(ALLSIGIDSGLOBAL) > 0:
        signature_ids = ALLSIGIDSGLOBAL
    else:
        signature_ids = sub.make_letter_ids(
            idlenth=normalised_solution.shape[0], mtype="Signature "
        )
    current_signatures = signature_ids

    # Apply the new threshold to the first and fifth values
    normalised_solution[current_signatures.index("SBS1")] = (
        0
        if normalised_solution[current_signatures.index("SBS1")] < clock_threshold
        else normalised_solution[current_signatures.index("SBS1")]
    )
    normalised_solution[current_signatures.index("SBS5")] = (
        0
        if normalised_solution[current_signatures.index("SBS5")] < clock_threshold
        else normalised_solution[current_signatures.index("SBS5")]
    )

    # Get indices that is not SBS1 and SBS5
    indices_except_sbs1_and_sbs5 = [
        i for i, sig in enumerate(current_signatures) if sig not in ["SBS1", "SBS5"]
    ]

    # Apply the original threshold to the rest of the values
    # Create a mask array for indices_except_sbs1_and_sbs5
    mask_indices = np.array(indices_except_sbs1_and_sbs5)
    mask = normalised_solution[mask_indices] < threshold
    normalised_solution[mask_indices[mask]] = 0

    threshold_solution = re_distribute_signature_activity(normalised_solution, genome)
    # print("this is the normalised solution after threshold:"+str(threshold_solution))
    return threshold_solution


def weighted_threshold_add_remove(finalactivities, W, M, threshold, ALLSIGIDSGLOBAL=[]):
    clock_threshold = -1  # make sure clock like signatures will not be removed
    # print("finalactivities before threshold:"+str(finalactivities))

    if len(ALLSIGIDSGLOBAL) > 0:
        signature_ids = ALLSIGIDSGLOBAL
    else:
        signature_ids = sub.make_letter_ids(
            idlenth=finalactivities.shape[0], mtype="Signature "
        )

    current_signatures = signature_ids
    # import pdb; pdb.set_trace();
    # Apply the new threshold to the first and fifth values
    finalactivities[current_signatures.index("SBS1")] = (
        0
        if finalactivities[current_signatures.index("SBS1")] < clock_threshold
        else finalactivities[current_signatures.index("SBS1")]
    )
    finalactivities[current_signatures.index("SBS5")] = (
        0
        if finalactivities[current_signatures.index("SBS5")] < clock_threshold
        else finalactivities[current_signatures.index("SBS5")]
    )

    # Get indices that is not SBS1 and SBS5
    indices_except_sbs1_and_sbs5 = [
        i for i, sig in enumerate(current_signatures) if sig not in ["SBS1", "SBS5"]
    ]

    # Apply the original threshold to the rest of the values
    # Create a mask array for indices_except_sbs1_and_sbs5
    mask_indices = np.array(indices_except_sbs1_and_sbs5)
    mask = finalactivities[mask_indices] < threshold
    finalactivities[mask_indices[mask]] = 0

    if 0 in finalactivities:
        finalactivities = re_distribute_signature_activity(finalactivities, M)
        S = np.dot(W, finalactivities)
        newFinalactivities = list(finalactivities)
        distance_after_threshold = np.linalg.norm(M - S, ord=2) / np.linalg.norm(
            M, ord=2
        )
        newFinalactivities = np.array(newFinalactivities)
        cosine_similarity = cos_sim(M, S)
        kldiv = round(scipy.stats.entropy(M, S), 3)
        correlation, _ = scipy.stats.pearsonr(M, S)
        correlation = round(correlation, 2)
        print("finalactivities after threshold:" + str(newFinalactivities))
        return (
            newFinalactivities,
            distance_after_threshold,
            cosine_similarity,
            kldiv,
            correlation,
        )

    else:
        return None


def add_connected_sigs(background_sigs, allsigids):
    connected_sigs = [
        ["SBS2", "SBS13"],
        ["SBS7a", "SBS7b", "SBS7c", "SBS7d"],
        ["SBS10a", "SBS10b"],
        ["SBS17a", "SBS17b"],
    ]

    backround_sig_names = sub.get_items_from_index(allsigids, background_sigs)
    connect_the_sigs = []
    for i in range(len(connected_sigs)):
        if len(set(connected_sigs[i]).intersection(set(backround_sig_names))) > 0:
            connect_the_sigs = connect_the_sigs + connected_sigs[i]

    backround_sig_names = list(set(backround_sig_names).union(set(connect_the_sigs)))
    background_sigs = sub.get_indeces(allsigids, backround_sig_names)
    background_sigs.sort()
    return background_sigs


def parameterized_objective2_custom(x, signatures, samples):
    rec = np.dot(signatures, x)
    try:
        y = LA.norm(samples - rec[:, np.newaxis])
    except:
        y = LA.norm(samples - rec)
    return y


def constraints1(x, samples):
    sumOfSamples = np.sum(samples, axis=0)
    # print(sumOfSamples)
    # print(x)
    result = sumOfSamples - (np.sum(x))
    # print (result)
    return result


def create_bounds(idxOfZeros, samples, numOfSignatures):
    total = np.sum(samples)
    b = (0.0, float(total))
    lst = [b] * numOfSignatures

    for i in idxOfZeros:
        lst[i] = (0.0, 0.0)

    return lst


# required for removing signatures with background
def get_changed_background_sig_idx(exposures, background_sigs):

    background_sigs_values = sub.get_items_from_index(exposures, background_sigs)
    temp_exposures = exposures[:]
    temp_exposures[:] = (value for value in temp_exposures if value != 0)

    # remove the background signatures with zero values
    background_sigs[:] = (value for value in background_sigs_values if value != 0)

    # get the new indices of the background signatures
    background_sigs = sub.get_indeces(temp_exposures, background_sigs_values)

    return background_sigs


# Fit signatures
def fit_signatures(W, genome, metric="l2"):
    # import pdb; pdb.set_trace();
    genome = np.array(genome)
    maxmutation = round(np.sum(genome))

    # initialize the guess
    x0 = np.random.rand(W.shape[1], 1) * maxmutation
    x0 = x0 / np.sum(x0) * maxmutation

    # the optimization step
    # set the bounds and constraints
    # bnds = create_bounds([], genome, W.shape[1])
    # cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genome]}

    # sol = scipy.optimize.minimize(parameterized_objective2_custom, x0, args=(W, genome),  bounds=bnds, constraints =cons1, tol=1e-30)
    # solution = sol.x

    # ### using NNLS algorithm
    # reg = nnls(W,genome)
    # weights = reg[0]
    # normalised_weights = weights/sum(weights)
    # solution = normalised_weights*sum(genome)
    # fit_threshold_solution = solution
    # est_genome = np.dot(W, fit_threshold_solution)

    ### using NNLS algorithm
    reg = nnls(W, genome)
    weights = reg[0]
    est_genome = np.dot(W, weights)
    normalised_weights = weights / sum(weights)
    solution = normalised_weights * sum(genome)

    # print(W1)
    # convert the newExposure vector into list type structure
    newExposure = list(solution)
    # print("This is the newExposure after fit_signature threshold:"+str(newExposure))

    # get the maximum value of the new Exposure
    maxcoef = max(newExposure)
    idxmaxcoef = newExposure.index(maxcoef)

    newExposure = np.round(newExposure)

    # We may need to tweak the maximum value of the new exposure to keep the total number of mutation equal to the original mutations in a genome
    if np.sum(newExposure) != maxmutation:
        newExposure[idxmaxcoef] = (
            round(newExposure[idxmaxcoef]) + maxmutation - sum(newExposure)
        )

    if metric == "cosine":
        newSimilarity = cos_sim(genome, est_genome)
    else:
        newSimilarity = np.linalg.norm(genome - est_genome, ord=2) / np.linalg.norm(
            genome, ord=2
        )

    return (newExposure, newSimilarity)


# to do the calculation, we need: W, samples(one column), H for the specific sample, tolerence
# Fit signatures with multiprocessing
def fit_signatures_pool(total_genome, W, index, metric="l2", threshold=-1):
    genome = total_genome[:, index]
    genome = np.array(genome)
    maxmutation = round(np.sum(genome))

    # initialize the guess
    x0 = np.random.rand(W.shape[1], 1) * maxmutation
    x0 = x0 / np.sum(x0) * maxmutation

    # the optimization step
    # set the bounds and constraints
    # bnds = create_bounds([], genome, W.shape[1])
    # cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genome]}

    # sol = scipy.optimize.minimize(parameterized_objective2_custom, x0, args=(W, genome),  bounds=bnds, constraints =cons1, tol=1e-30)
    # solution = sol.x

    ### using NNLS algorithm
    reg = nnls(W, genome)
    weights = reg[0]
    # est_genome = np.dot(W, weights)
    # normalised_weights = weights/sum(weights)
    # solution = normalised_weights*sum(genome)
    normalised_weights = weights / sum(weights)
    solution = normalised_weights * sum(genome)
    fit_threshold_solution = weighted_threshold(solution, genome, threshold)
    est_genome = np.dot(W, fit_threshold_solution)

    # print(W1)
    # convert the newExposure vector into list type structure
    newExposure = list(solution)

    # get the maximum value of the new Exposure
    maxcoef = max(newExposure)
    idxmaxcoef = newExposure.index(maxcoef)

    newExposure = np.round(newExposure)

    # We may need to tweak the maximum value of the new exposure to keep the total number of mutation equal to the original mutations in a genome
    if np.sum(newExposure) != maxmutation:
        newExposure[idxmaxcoef] = (
            round(newExposure[idxmaxcoef]) + maxmutation - sum(newExposure)
        )

    if metric == "cosine":
        newSimilarity = cos_sim(genome, est_genome)
    else:
        newSimilarity = np.linalg.norm(genome - est_genome, ord=2) / np.linalg.norm(
            genome, ord=2
        )

    return (newExposure, newSimilarity)


def add_signatures(
    W,
    genome,
    cutoff=0.05,
    presentSignatures=[],
    toBeAdded="all",
    metric="l2",
    solver="nnls",
    check_rule_negatives=[],
    check_rule_penalty=1.0,
    verbose=False,
):

    # This function takes an array of signature and a single genome as input, returns a dictionray of cosine similarity, exposures and presence
    # of signatures according to the indices of the original signature array

    if toBeAdded == "all":
        notToBeAdded = []
    else:
        # print(len(list(range(W.shape[1]))))
        # print(len(toBeAdded))
        notToBeAdded = list(
            set(list(range(W.shape[1]))) - set(toBeAdded)
        )  # get the indices of the signatures to be excluded

        # print(len(notToBeAdded))
    originalSimilarity = (
        np.inf
    )  # it can be also written as oldsimilarity. wanted to put a big number
    maxmutation = round(np.sum(genome))
    init_listed_idx = presentSignatures
    # print("The init_listed_idx in add_loop:"+str(init_listed_idx))

    init_nonlisted_idx = list(
        set(list(range(W.shape[1]))) - set(init_listed_idx) - set(notToBeAdded)
    )

    finalRecord = [
        ["similarity place-holder"],
        ["newExposure place-holder"],
        ["signatures place-holder"],
    ]  # for recording the cosine difference, similarity, the new exposure and the index of the best signauture

    # get the initial original similarity if some signatures are already given to be present
    if len(init_listed_idx) != 0:
        loop_liststed_idx = init_listed_idx
        loop_liststed_idx.sort()
        W1 = W[:, loop_liststed_idx]

        # initialize the guess
        x0 = np.random.rand(W1.shape[1], 1) * maxmutation
        x0 = x0 / np.sum(x0) * maxmutation

        # set the bounds and constraints
        bnds = create_bounds([], genome, W1.shape[1])
        cons1 = {"type": "eq", "fun": constraints1, "args": [genome]}

        # the optimization step
        if solver == "slsqp":
            sol = minimize(
                parameterized_objective2_custom,
                x0,
                args=(W1, genome),
                bounds=bnds,
                constraints=cons1,
                tol=1e-30,
            )
            newExposure = list(sol.x)
            est_genome = np.dot(W1, np.array(newExposure))
        if solver == "nnls":
            ### using NNLS algorithm
            reg = nnls(W1, genome[:, 0])
            weights = reg[0]
            est_genome = np.dot(W1, weights)
            normalised_weights = weights / sum(weights)
            solution = normalised_weights * sum(np.sum(genome, axis=0))
            newExposure = list(solution)
        # if solver == "nnls":
        # ### using NNLS algorithm
        #     reg = nnls(W1, genome[:,0])
        #     weights = reg[0]
        #     normalised_weights = weights/sum(weights)
        #     solution = normalised_weights*sum(np.sum(genome, axis=0))
        #     est_genome = np.dot(W1, solution)
        #     newExposure = list(solution)
        # print("add_signatures_newExposure"+str(newExposure))

        # print(W1)
        # convert the newExposure vector into list type structure

        # newExposure = list(solution) #for nnls

        # get the maximum value of the new Exposure
        maxcoef = max(newExposure)
        idxmaxcoef = newExposure.index(maxcoef)

        newExposure = np.round(newExposure)

        # print("Value in add_loop_1 after nnls:" + str(newExposure))

        # We may need to tweak the maximum value of the new exposure to keep the total number of mutation equal to the original mutations in a genome
        if np.sum(newExposure) != maxmutation:
            newExposure[idxmaxcoef] = (
                round(newExposure[idxmaxcoef]) + maxmutation - sum(newExposure)
            )

        # print("newExposure after maxmutation:" + str(newExposure))

        if metric == "cosine":
            originalSimilarity = 1 - cos_sim(genome[:, 0], est_genome)
        elif metric == "l2":
            originalSimilarity = np.linalg.norm(
                genome[:, 0] - est_genome, ord=2
            ) / np.linalg.norm(genome[:, 0], ord=2)

        finalRecord = [originalSimilarity, newExposure, init_listed_idx]

    while True:

        bestDifference = -np.inf  #### wanted to put a big initila number
        bestSimilarity = np.inf  #### actually bestdistance
        loopRecord = [
            ["newExposure place-holder"],
            ["signatures place-holder"],
            ["best loop signature place-holder"],
        ]

        for sig in init_nonlisted_idx:

            if len(init_listed_idx) != 0:
                loop_liststed_idx = init_listed_idx + [sig]
                loop_liststed_idx.sort()
                # print(loop_liststed_idx)
                W1 = W[:, loop_liststed_idx]
                # print (W1.shape)
                # initialize the guess
                x0 = np.random.rand(W1.shape[1], 1) * maxmutation
                x0 = x0 / np.sum(x0) * maxmutation

                # set the bounds and constraints
                bnds = create_bounds([], genome, W1.shape[1])
                cons1 = {"type": "eq", "fun": constraints1, "args": [genome]}
            # for the first time addintion
            else:
                W1 = W[:, sig][:, np.newaxis]
                # print (W1.shape)
                # initialize the guess
                x0 = np.ones((1, 1)) * maxmutation

                # set the bounds and constraints
                bnds = create_bounds([], genome, 1)
                cons1 = {"type": "eq", "fun": constraints1, "args": [genome]}

            # print("Value in add_loop_2 before nnls:" + str(newExposure))

            if solver == "slsqp":
                # the optimization step
                sol = minimize(
                    parameterized_objective2_custom,
                    x0,
                    args=(W1, genome),
                    bounds=bnds,
                    constraints=cons1,
                    tol=1e-30,
                )
                newExposure = list(sol.x)
                # compute the estimated genome
                est_genome = np.dot(W1, np.array(newExposure))
            if solver == "nnls":
                reg = nnls(W1, genome[:, 0])
                weights = reg[0]
                est_genome = np.dot(W1, weights)
                normalised_weights = weights / sum(weights)
                solution = normalised_weights * sum(genome)
                newExposure = list(solution)  # for nnls
            # if solver == "nnls":
            #     reg = nnls(W1, genome[:,0])
            #     weights = reg[0]
            #     normalised_weights = weights/sum(weights)
            #     solution = normalised_weights*sum(genome)
            #     est_genome = np.dot(W1, solution)
            #     newExposure = list(solution)
            # convert the newExposure vector into list type structure
            # print("Value in add_loop_2 after nnls:" + str(newExposure))

            # get the maximum value of the new Exposure
            maxcoef = max(newExposure)
            idxmaxcoef = newExposure.index(maxcoef)
            newExposure = np.round(newExposure)

            # We may need to tweak the maximum value of the new exposure to keep the total number of mutation equal to the original mutations in a genome
            if np.sum(newExposure) != maxmutation:
                newExposure[idxmaxcoef] = (
                    round(newExposure[idxmaxcoef]) + maxmutation - sum(newExposure)
                )

            if metric == "cosine":
                newSimilarity = 1 - cos_sim(genome[:, 0], est_genome)
            elif metric == "l2":
                newSimilarity = np.linalg.norm(
                    genome[:, 0] - est_genome, ord=2
                ) / np.linalg.norm(genome[:, 0], ord=2)

            if sig in check_rule_negatives:

                newSimilarity = newSimilarity * check_rule_penalty

            difference = originalSimilarity - newSimilarity

            # record the best values so far that creates the best difference
            if difference > bestDifference:
                bestDifference = difference
                bestSimilarity = newSimilarity
                loopRecord = [
                    newExposure,
                    W1,
                    sig,
                ]  # recording the cosine difference, the new exposure and the index of the best signauture

        if originalSimilarity - bestSimilarity > cutoff:

            originalSimilarity = bestSimilarity
            init_listed_idx.append(
                loopRecord[2]
            )  # add the qualified signature in the signature list
            init_nonlisted_idx.remove(
                loopRecord[2]
            )  # remover the qualified signature from the to be added list
            init_listed_idx.sort()
            # print(originalSimilarity)
            finalRecord = [
                originalSimilarity,
                loopRecord[0],
                init_listed_idx,
                loopRecord[1],
                genome,
            ]
            if verbose == True:
                print("Current signaure indices: ", init_listed_idx)
                if metric == "cosine":
                    print("Cosine similarity: ", 1 - originalSimilarity)
                elif metric == "l2":
                    print("L2 norm: ", originalSimilarity)
                print("\n")
            if len(init_nonlisted_idx) != 0:

                continue
            else:
                break
        else:

            break

    # print("Good so far")
    # print(finalRecord)

    dictExposure = {
        "similarity": finalRecord[0],
        "exposures": finalRecord[1],
        "signatures": finalRecord[2],
    }
    addExposure = np.zeros([W.shape[1]])
    addExposure[dictExposure["signatures"]] = dictExposure["exposures"]
    # print(cos_sim(genome[:,0], np.dot(W,addExposure)))
    if metric == "cosine":
        return (
            addExposure,
            cos_sim(genome[:, 0], np.dot(W, addExposure)),
            cos_sim(genome[:, 0], np.dot(W, addExposure)),
        )  # first value is the exprosure, second value is the similarity (based on the similarity matic), third value is the cosine similarity
    elif metric == "l2":
        return (
            addExposure,
            originalSimilarity,
            cos_sim(genome[:, 0], np.dot(W, addExposure)),
        )  # first value is the exprosure, second value is the similarity (based on the similarity matic), third value is the cosine similarity


def remove_all_single_signatures(
    W,
    H,
    genomes,
    metric="l2",
    solver="nnls",
    cutoff=0.05,
    background_sigs=[],
    verbose=False,
):
    # make the empty list of the successfull combinations

    # print("This is the initial H in remove_loop:"+str(H))
    signature_ids = sub.make_letter_ids(idlenth=W.shape[1], mtype="Signature ")

    current_signatures = signature_ids
    # print(current_signatures)
    successList = [0, [], 0]
    background_sig = copy.deepcopy(background_sigs)
    # print("This is the background_sig in remove:",background_sig)

    # setting the base_similarity
    base_H_index = list(np.nonzero(H)[0])
    # print("This is the base_H_index:"+str(base_H_index))

    reg = nnls(W[:, base_H_index], genomes)
    base_weights = reg[0]
    # print("This is the base_weight after nnls:"+str(base_weights))

    if metric == "cosine":
        # get the cos_similarity with sample for the oringinal W and H[:,i]
        originalSimilarity = 1 - cos_sim(
            genomes, np.dot(W[:, base_H_index], base_weights)
        )
        if verbose == True:
            print("originalSimilarity", 1 - originalSimilarity)
    elif metric == "l2":

        originalSimilarity = np.linalg.norm(
            genomes - np.dot(W[:, base_H_index], base_weights), ord=2
        ) / np.linalg.norm(genomes, ord=2)
        if verbose == True:
            print("originalSimilarity", originalSimilarity)
    # make the original exposures of specific sample round
    oldExposures = np.round(H)

    # set the flag for the while loop
    if len(oldExposures[np.nonzero(oldExposures)]) > 1:
        Flag = True
    else:
        Flag = False
        if metric == "cosine":
            return (
                oldExposures,
                1 - originalSimilarity,
                cos_sim(genomes, np.dot(W, H)),
            )  # first value is the exprosure, second value is the similarity (based on the similarity matic), third value is the cosine similarity
        else:
            return (
                oldExposures,
                originalSimilarity,
                cos_sim(genomes, np.dot(W, H)),
            )  # first value is the exprosure, second value is the similarity (based on the similarity matic), third value is the cosine similarity
    # The while loop starts here
    current_signatures = sub.get_items_from_index(
        signature_ids, np.nonzero(oldExposures)[0]
    )

    while Flag:
        # get the list of the indices those already have zero values
        if len(successList[1]) == 0:
            initialZerosIdx = list(np.where(oldExposures == 0)[0])
            # get the indices to be selected
            selectableIdx = list(np.where(oldExposures > 0)[0])
        elif len(successList[1]) > 1:
            initialZerosIdx = list(np.where(successList[1] == 0)[0])
            # get the indices to be selected
            selectableIdx = list(np.where(successList[1] > 0)[0])
        else:
            print("iteration is completed")
            # break

        # get the total mutations for the given sample
        maxmutation = round(np.sum(genomes))

        # new signature matrix omiting the column for zero
        # Winit = np.delete(W, initialZerosIdx, 1)
        Winit = W[:, selectableIdx]

        # set the initial cos_similarity or other similarity distance
        record = [np.inf, [], 1]  #
        # get the number of current nonzeros
        l = Winit.shape[1]

        background_sig = get_changed_background_sig_idx(
            list(oldExposures), background_sig
        )

        for i in range(l):

            if i in background_sig:
                continue

            loopSelection = list(range(l))
            del loopSelection[i]

            W1 = Winit[:, loopSelection]

            # initialize the guess
            x0 = np.random.rand(l - 1, 1) * maxmutation
            x0 = x0 / np.sum(x0) * maxmutation

            # set the bounds and constraints
            bnds = create_bounds([], genomes, W1.shape[1])
            cons1 = {"type": "eq", "fun": constraints1, "args": [genomes]}

            # the optimization step
            if solver == "slsqp":
                sol = minimize(
                    parameterized_objective2_custom,
                    x0,
                    args=(W1, genomes),
                    bounds=bnds,
                    constraints=cons1,
                    tol=1e-15,
                )

                # print (sol.success)
                # print (sol.x)

                # convert the newExposure vector into list type structure
                newExposure = list(sol.x)
            if solver == "nnls":
                ### using NNLS algorithm
                reg = nnls(W1, genomes)
                weights = reg[0]
                newSample = np.dot(W1, weights)
                normalised_weights = weights / sum(weights)
                solution = normalised_weights * sum(genomes)
                newExposure = list(solution)
            # if solver == "nnls":
            #     ### using NNLS algorithm
            #     reg = nnls(W1,genomes)
            #     weights = reg[0]
            #     normalised_weights = weights/sum(weights)
            #     solution = normalised_weights*sum(genomes)
            #     newSample = np.dot(W1, solution)
            #     newExposure = list(solution)
            #     # print("newExposure remove_all_sig:" + str(newExposure))

            # insert the loopZeros in its actual position
            newExposure.insert(i, 0)

            # insert zeros in the required position the newExposure matrix
            initialZerosIdx.sort()

            for zeros in initialZerosIdx:
                newExposure.insert(zeros, 0)

            # get the maximum value the new Exposure
            # maxcoef = max(newExposure)
            # idxmaxcoef = newExposure.index(maxcoef)

            # newExposure = np.round(newExposure)

            # if np.sum(newExposure)!=maxmutation:
            # newExposure[idxmaxcoef] = round(newExposure[idxmaxcoef])+maxmutation-sum(newExposure)

            newExposure = np.array(newExposure)

            if verbose == True:
                # print(newExposure)
                print("\nRemoving {}".format(current_signatures[i]))
            if metric == "cosine":
                newSimilarity = 1 - cos_sim(genomes, newSample)
                if verbose == True:
                    print("newSimilarity", 1 - newSimilarity)
            elif metric == "l2":
                newSimilarity = np.linalg.norm(
                    genomes - newSample, ord=2
                ) / np.linalg.norm(genomes, ord=2)
                if verbose == True:
                    print("newSimilarity", newSimilarity)
            difference = newSimilarity - originalSimilarity
            if difference < 0:
                difference = cutoff + 1e-100
            if verbose == True:
                print("difference", difference)
            if difference < record[0]:
                record = [difference, newExposure, newSimilarity]

                # print("difference", difference)
            if verbose == True:
                print("--------------------------------------")

        if verbose == True:
            print("\n############################\n############################")
            # print("Selected Exposure")
            # print(record[1])
            selected_sigs = sub.get_items_from_index(
                signature_ids, np.nonzero(record[1])[0]
            )
            dropped_sig = list(set(current_signatures) - set(selected_sigs))
            print("Dropped Signature: {}".format(dropped_sig))
            current_signatures = selected_sigs
            print("New Similarity", record[2])
            print("Similarity Difference: {}".format(record[2] - originalSimilarity))
            print("Current Signatures: {}".format(selected_sigs))
            print("****************************\n****************************\n\n")
            # print ("This loop's selection is {}".format(record))

        # print("This is the cutoff:"+str(record[0]))

        if record[0] > cutoff:
            Flag = False
        elif len(record[1][np.nonzero(record[1])]) == 1:
            successList = record
            Flag = False
        else:
            successList = record
            background_sig = get_changed_background_sig_idx(
                list(record[1]), background_sig
            )
            originalSimilarity = record[2]
        # print("background_sig at the end of remove", background_sig)
        # print("The loop selection is {}".format(successList))

        # print (Flag)
        # print ("\n\n")

    # print ("The final selection is {}".format(successList))

    if len(successList[1]) == 0:
        successList = [0.0, oldExposures, originalSimilarity]

    if verbose == True:
        print("\n")
        print("Final Exposure")
        print(successList[1])
        if metric == "cosine":
            print("Final Similarity: ", 1 - successList[2])
        if metric == "l2":
            print("Final Similarity: ", successList[2])
    H = successList[1]
    # print(H)
    return (
        H,
        successList[2],
        cos_sim(genomes, np.dot(W, H)),
    )  # first value is the exprosure, second value is the similarity (based on the similarity matic), third value is the cosine similarity


def initial_remove_all_single_signatures(
    W,
    H,
    genomes,
    metric="l2",
    solver="nnls",
    cutoff=0.05,
    background_sigs=[],
    verbose=False,
    threshold=30,
    ALLSIGIDSGLOBAL=[],
):
    # make the empty list of the successfull combinations
    # import pdb; pdb.set_trace();
    # print("threshold in initial remove all signatures:"+str(W))
    # print("This is the W in initial remove:"+str(W.shape))
    # print("This is the initial H in remove_loop:"+str(H))
    if len(ALLSIGIDSGLOBAL) > 0:
        signature_ids = ALLSIGIDSGLOBAL
    else:
        signature_ids = sub.make_letter_ids(idlenth=W.shape[1], mtype="Signature ")
    current_signatures = signature_ids
    # print(current_signatures)
    successList = [0, [], 0]
    background_sig = copy.deepcopy(background_sigs)
    # print("background sigs for initial removal:"+str(background_sigs))

    # setting the base_similarity
    base_H_index = list(np.nonzero(H)[0])
    reg = nnls(W[:, base_H_index], genomes)
    # print("reg[0] for 1st regular nnls:"+str(reg[0]))

    base_weights = reg[0]

    if metric == "cosine":
        # get the cos_similarity with sample for the oringinal W and H[:,i]
        originalSimilarity = 1 - cos_sim(
            genomes, np.dot(W[:, base_H_index], base_weights)
        )
        if verbose == True:
            print("originalSimilarity", 1 - originalSimilarity)
    elif metric == "l2":
        # norm_solution = (base_weights/sum(base_weights))*sum(genomes)
        # originalSimilarity = np.linalg.norm(genomes-np.dot(W[:,base_H_index], norm_solution) , ord=2)/np.linalg.norm(genomes, ord=2)
        originalSimilarity = np.linalg.norm(
            genomes - np.dot(W[:, base_H_index], base_weights), ord=2
        ) / np.linalg.norm(genomes, ord=2)
        if verbose == True:
            print("originalSimilarity", originalSimilarity)
    # make the original exposures of specific sample round
    oldExposures = np.round(H)

    # import pdb; pdb.set_trace();
    # set the flag for the while loop
    if len(oldExposures[np.nonzero(oldExposures)]) > 1:
        Flag = True
    else:
        Flag = False
        if metric == "cosine":
            return (
                oldExposures,
                1 - originalSimilarity,
                cos_sim(genomes, np.dot(W, H)),
            )  # first value is the exprosure, second value is the similarity (based on the similarity matic), third value is the cosine similarity
        else:
            return (
                oldExposures,
                originalSimilarity,
                cos_sim(genomes, np.dot(W, H)),
            )  # first value is the exprosure, second value is the similarity (based on the similarity matic), third value is the cosine similarity
    # The while loop starts here
    current_signatures = sub.get_items_from_index(
        signature_ids, np.nonzero(oldExposures)[0]
    )

    while Flag:
        # print("initial Flag:" + str(Flag))
        # get the list of the indices those already have zero values
        if len(successList[1]) == 0:
            initialZerosIdx = list(np.where(oldExposures == 0)[0])
            # get the indices to be selected
            selectableIdx = list(np.where(oldExposures > 0)[0])
        elif len(successList[1]) > 1:
            initialZerosIdx = list(np.where(successList[1] == 0)[0])
            # get the indices to be selected
            selectableIdx = list(np.where(successList[1] > 0)[0])
        else:
            print("iteration is completed")
            # break

        # get the total mutations for the given sample
        maxmutation = round(np.sum(genomes))

        # new signature matrix omiting the column for zero
        # Winit = np.delete(W, initialZerosIdx, 1)
        Winit = W[:, selectableIdx]

        # set the initial cos_similarity or other similarity distance
        record = [np.inf, [], 1]  #
        # get the number of current nonzeros
        l = Winit.shape[1]

        # print("current non-zeros:"+str(l))

        background_sig = get_changed_background_sig_idx(
            list(oldExposures), background_sig
        )
        # print("2 nd background sigs for initial removal:"+str(background_sigs))

        for i in range(l):

            if i in background_sig:
                continue

            loopSelection = list(range(l))
            del loopSelection[i]

            W1 = Winit[:, loopSelection]

            # initialize the guess
            x0 = np.random.rand(l - 1, 1) * maxmutation
            x0 = x0 / np.sum(x0) * maxmutation

            # set the bounds and constraints
            bnds = create_bounds([], genomes, W1.shape[1])
            cons1 = {"type": "eq", "fun": constraints1, "args": [genomes]}

            # the optimization step
            if solver == "slsqp":
                sol = minimize(
                    parameterized_objective2_custom,
                    x0,
                    args=(W1, genomes),
                    bounds=bnds,
                    constraints=cons1,
                    tol=1e-15,
                )

                # print (sol.success)
                # print (sol.x)

                # convert the newExposure vector into list type structure
                newExposure = list(sol.x)
            if solver == "nnls":
                ### using NNLS algorithm
                reg = nnls(W1, genomes)
                weights = reg[0]
                newSample = np.dot(W1, weights)
                normalised_weights = weights / sum(weights)
                solution = normalised_weights * sum(genomes)
                newExposure = list(solution)
            # if solver == "nnls":
            #     ### using NNLS algorithm
            #     # import pdb;pdb.set_trace()
            #     reg = nnls(W1,genomes)
            #     weights = reg[0]
            #     # newSample = np.dot(W1, weights)
            #     # normalised_weights = weights/sum(weights)
            #     # solution = normalised_weights*sum(genomes)
            #     # newExposure = list(solution)
            #     normalised_weights = weights/sum(weights)
            #     solution = normalised_weights*sum(genomes)
            #     # threshold_solution = weighted_threshold(solution, genomes, threshold=-1)
            #     threshold_solution = solution
            #     newSample = np.dot(W1, threshold_solution)
            #     newExposure = list(threshold_solution)

            # insert the loopZeros in its actual position
            newExposure.insert(i, 0)

            # insert zeros in the required position the newExposure matrix
            initialZerosIdx.sort()

            for zeros in initialZerosIdx:
                newExposure.insert(zeros, 0)

            # get the maximum value the new Exposure
            # maxcoef = max(newExposure)
            # idxmaxcoef = newExposure.index(maxcoef)

            # newExposure = np.round(newExposure)

            # if np.sum(newExposure)!=maxmutation:
            # newExposure[idxmaxcoef] = round(newExposure[idxmaxcoef])+maxmutation-sum(newExposure)

            newExposure = np.array(newExposure)

            # print("newExposure after adding zeros:" + str(newExposure))

            if verbose == True:
                # print(newExposure)
                print("\nRemoving {}".format(current_signatures[i]))
                print(
                    "Atribution", sub.non_zero_values_dict(newExposure, signature_ids)
                )
            if metric == "cosine":
                newSimilarity = 1 - cos_sim(genomes, newSample)
                if verbose == True:
                    print("newSimilarity", 1 - newSimilarity)
            elif metric == "l2":
                newSimilarity = np.linalg.norm(
                    genomes - newSample, ord=2
                ) / np.linalg.norm(genomes, ord=2)
                if verbose == True:
                    print("newSimilarity", newSimilarity)
            difference = newSimilarity - originalSimilarity
            if difference < 0:
                difference = cutoff + 1e-100
            if verbose == True:
                print("difference", difference)
            if difference < record[0]:
                # print("WOW!!difference<record[0]")
                record = [difference, newExposure, newSimilarity]

                # print("difference", difference)
            if verbose == True:
                print("--------------------------------------")
            # print("FOR LOOP" , l)

        # print("While loop ", Flag)
        # pdb.set_trace()
        if verbose == True:
            print("\n############################\n############################")
            # print("Selected Exposure")
            # print(record[1])
            selected_sigs = sub.get_items_from_index(
                signature_ids, np.nonzero(record[1])[0]
            )
            dropped_sig = list(set(current_signatures) - set(selected_sigs))
            print("Dropped Signature: {}".format(dropped_sig))
            current_signatures = selected_sigs
            print("New Similarity", record[2])
            print("Similarity Difference: {}".format(record[2] - originalSimilarity))
            print("Current Signatures: {}".format(selected_sigs))
            print("****************************\n****************************\n\n")

        if record[0] > cutoff:
            # print("!!!!! record[0]>cutoff !!!!!")
            # print(record[0])
            Flag = False
        elif len(record[1][np.nonzero(record[1])]) == 1:
            # print("!!! len(record[1][np.nonzero(record[1])])==1 !!!")
            successList = record
            # print("successList/record in loop:"+str(successList))
            Flag = False
        else:
            successList = record
            # print("successList/record in loop:"+str(successList))
            background_sig = get_changed_background_sig_idx(
                list(record[1]), background_sig
            )
            originalSimilarity = record[2]

        # print("The loop selection is {}".format(successList))

        # print ("for i in range loop Flag:" + str(Flag))
        # print ("\n\n")

    # print("successList final:"+str(successList))
    if len(successList[1]) == 0:
        successList = [0.0, oldExposures, originalSimilarity]

    if verbose == True:
        print("\n")
        print("Final Exposure")
        print(successList[1])
        if metric == "cosine":
            print("Final Similarity: ", 1 - successList[2])
        if metric == "l2":
            print("Final Similarity: ", successList[2])
    H = successList[1]
    # print("newExposure outside the loop:" + str(newExposure))
    # print("H:"+str(H))
    # print('SBS10d',dict_sigvals['SBS10d'],'--','SBS18',dict_sigvals['SBS18'])
    # print("H:"+ str(H))
    # print("successList:"+str(successList))
    # print("H outside loop before threshold:" + str(H))

    H_thresh = weighted_threshold_split(H, genomes, threshold, ALLSIGIDSGLOBAL)
    # print("This is the H in initial remove:"+str(H.shape))
    # print("H_thresh:" + str(H_thresh))
    # import pdb;pdb.set_trace()
    final_l2 = np.linalg.norm(genomes - np.dot(W, H_thresh), ord=2) / np.linalg.norm(
        genomes, ord=2
    )
    return (
        H_thresh,
        final_l2,
        cos_sim(genomes, np.dot(W, H_thresh)),
    )  # first value is the exposure, second value is the similarity (based on the similarity matic), third value is the cosine similarity


def remove_all_single_signatures_pool(indices, W, exposures, totoalgenomes):
    i = indices
    H = exposures[:, i]
    genomes = totoalgenomes[:, i]
    # make the empty list of the successfull combinations
    successList = [0, [], 0]
    # get the cos_similarity with sample for the oringinal W and H[:,i]
    originalSimilarity = cos_sim(genomes, np.dot(W, H))
    # make the original exposures of specific sample round
    oldExposures = np.round(H)

    # set the flag for the while loop
    if len(oldExposures[np.nonzero(oldExposures)]) > 1:
        Flag = True
    else:
        Flag = False
        return oldExposures
    # The while loop starts here
    while Flag:

        # get the list of the indices those already have zero values
        if len(successList[1]) == 0:
            initialZerosIdx = list(np.where(oldExposures == 0)[0])
            # get the indices to be selected
            selectableIdx = list(np.where(oldExposures > 0)[0])
        elif len(successList[1]) > 1:
            initialZerosIdx = list(np.where(successList[1] == 0)[0])
            # get the indices to be selected
            selectableIdx = list(np.where(successList[1] > 0)[0])
        else:
            print("iteration is completed")
            # break

        # get the total mutations for the given sample
        maxmutation = round(np.sum(genomes))

        # new signature matrix omiting the column for zero
        # Winit = np.delete(W, initialZerosIdx, 1)
        Winit = W[:, selectableIdx]

        # set the initial cos_similarity
        record = [0.11, []]
        # get the number of current nonzeros
        l = Winit.shape[1]

        for i in range(l):
            # print(i)
            loopSelection = list(range(l))
            del loopSelection[i]
            # print (loopSelection)
            W1 = Winit[:, loopSelection]

            # initialize the guess
            x0 = np.random.rand(l - 1, 1) * maxmutation
            x0 = x0 / np.sum(x0) * maxmutation

            # set the bounds and constraints
            bnds = create_bounds([], genomes, W1.shape[1])
            cons1 = {"type": "eq", "fun": constraints1, "args": [genomes]}

            # the optimization step
            sol = minimize(
                parameterized_objective2_custom,
                x0,
                args=(W1, genomes),
                bounds=bnds,
                constraints=cons1,
                tol=1e-15,
            )

            # print (sol.success)
            # print (sol.x)

            # convert the newExposure vector into list type structure
            newExposure = list(sol.x)

            # insert the loopZeros in its actual position
            newExposure.insert(i, 0)

            # insert zeros in the required position the newExposure matrix
            initialZerosIdx.sort()
            for zeros in initialZerosIdx:
                newExposure.insert(zeros, 0)

            # get the maximum value the new Exposure
            maxcoef = max(newExposure)
            idxmaxcoef = newExposure.index(maxcoef)

            newExposure = np.round(newExposure)

            if np.sum(newExposure) != maxmutation:
                newExposure[idxmaxcoef] = (
                    round(newExposure[idxmaxcoef]) + maxmutation - sum(newExposure)
                )

            newSample = np.dot(W, newExposure)
            newSimilarity = cos_sim(genomes, newSample)

            difference = originalSimilarity - newSimilarity
            # print(originalSimilarity)
            # print(newSample)
            # print(newExposure)
            # print(newSimilarity)

            # print(difference)
            # print (newExposure)
            # print (np.round(H))
            # print ("\n\n")

            if difference < record[0]:
                record = [difference, newExposure, newSimilarity]

        # print ("This loop's selection is {}".format(record))

        if record[0] > 0.01:
            Flag = False
        elif len(record[1][np.nonzero(record[1])]) == 1:
            successList = record
            Flag = False
        else:
            successList = record
        # print("The loop selection is {}".format(successList))

        # print (Flag)
        # print ("\n\n")

    # print ("The final selection is {}".format(successList))

    if len(successList[1]) == 0:
        successList = [0.0, oldExposures, originalSimilarity]
    # print ("one sample completed")
    return successList[1]


def add_remove_signatures(
    W,
    sample,
    metric="l2",
    solver="nnls",
    background_sigs=[],
    permanent_sigs=[],
    candidate_sigs="all",
    add_penalty=0.05,
    remove_penalty=0.01,
    check_rule_negatives=[],
    checkrule_penalty=1.00,
    allsigids=False,
    directory=os.getcwd() + "/Signature_assaignment_logfile.txt",
    connected_sigs=True,
    verbose=False,
    threshold=0,
    ALLSIGIDSGLOBAL=[],
):

    lognote = open(directory, "a")
    print("threshold in add_remove_signatures:" + str(threshold))
    maxmutation = np.sum(np.array(sample))
    # print("######## OBSERVE FROM HERE #######")

    always_background = copy.deepcopy(permanent_sigs)
    M = sample
    # print("M=sample="+str(M))
    if candidate_sigs == "all":
        candidate_sigs = list(range(W.shape[1]))
    # first add each signatures
    original_distance = np.inf  # a big number
    layer = 0

    # check the cosine_similarity with 4 signatures (highest)
    cosine_similarity_with_four_signatures = (
        1.0  # a value that will allow the signature not be reported as novel signature
    )

    # set the signature's name
    if type(allsigids) != bool:
        allsigids = sub.get_items_from_index(allsigids, candidate_sigs)
    else:
        allsigids = candidate_sigs
    # print("allsigs before while loop:",allsigids)
    # print("backsigs before while loop:",background_sigs)
    while True:
        if verbose:
            print("\n\n\n\n!!!!!!!!!!!!!!!!!!STARTING LAYER: ", layer)
        lognote.write(
            "\n\n!!!!!!!!!!!!!!!!!!!!!!!!! LAYER: {} !!!!!!!!!!!!!!!!!!!!!!!!!\n".format(
                layer
            )
        )
        layer = layer + 1
        layer_original_distance = np.inf  # a big number
        sigsToBeAdded = list(set(candidate_sigs) - set(background_sigs))
        # print("sigsToBeAdded inside while loop:",sigsToBeAdded)
        # set the signature's name
        if type(allsigids) != bool:
            allsigidsToBeAdded = sub.get_items_from_index(allsigids, sigsToBeAdded)
        else:
            allsigidsToBeAdded = sigsToBeAdded

        for i, j in zip(sigsToBeAdded, allsigidsToBeAdded):
            loop_sig = [i]

            background_sigs = list(set(background_sigs).union(set(always_background)))
            if connected_sigs == True:
                background_sigs = add_connected_sigs(background_sigs, list(allsigids))
            # print("background_sigs inside for loop:",background_sigs)
            add_exposures, add_distance, _ = add_signatures(
                W,
                M[:, np.newaxis],
                presentSignatures=copy.deepcopy(background_sigs),
                toBeAdded=loop_sig,
                metric="l2",
                verbose=False,
                check_rule_negatives=check_rule_negatives,
                check_rule_penalty=checkrule_penalty,
                cutoff=add_penalty,
            )
            # import pdb;pdb.set_trace()
            # print("After add loop: ", np.nonzero(add_exposures)[0]," errror : ",add_distance)
            if verbose:
                print(
                    "\n\n\n################## Add Index {} ########################".format(
                        j
                    )
                )
                print(np.nonzero(add_exposures)[0])
                print(sigsToBeAdded)
                print(add_distance)

            remove_exposures, remove_distance, _ = remove_all_single_signatures(
                W,
                add_exposures,
                M,
                background_sigs=always_background,
                metric="l2",
                verbose=False,
                cutoff=remove_penalty,
            )  # initial removal
            # print("always_background inside for loop:",always_background)
            # import pdb;pdb.set_trace()
            # print("After remove loop: ",np.nonzero(remove_exposures)[0]," errror : ",remove_distance)
            if verbose:
                print("\n################## Remove ########################")
                print(np.nonzero(remove_exposures)[0])
                print(remove_distance)
            # check if there is any change between the add and remove signatures and assign the distance and exposure accoridingly
            if (
                np.nonzero(add_exposures)[0].all()
                == np.nonzero(remove_exposures)[0].all()
            ) and np.nonzero(add_exposures)[0].shape == np.nonzero(remove_exposures)[
                0
            ].shape:
                distance = add_distance
                exposures = add_exposures
            else:
                distance = remove_distance
                exposures = remove_exposures

            if distance < layer_original_distance:
                selected_signatures = np.nonzero(exposures)[0]
                layer_original_distance = distance
                activities = exposures
        if verbose:
            print(
                "\n\n#################### After Add-Remove #################################"
            )
            print(selected_signatures)
            print(layer_original_distance)

        lognote.write(
            "Best Signature Composition {}\n".format(
                sub.get_items_from_index(allsigids, selected_signatures)
            )
        )
        lognote.write("L2 Error % {}\n".format(round(layer_original_distance, 2)))
        lognote.write(
            "Cosine Similarity {}\n".format(round(cos_sim(M, np.dot(W, activities)), 2))
        )
        if layer_original_distance < original_distance:  # original_distance:inf
            # print("layer_original_distance when layer_original_distance < original_distance:"+str(layer_original_distance))
            # print("original_distance when layer_original_distance < original_distance:"+str(original_distance))
            original_distance = layer_original_distance
            background_sigs = list(selected_signatures)
            finalactivities = activities
            # print("Final activities when layer_original_distance < original_distance:"+str(finalactivities))
            if len(background_sigs) == 4:
                cosine_similarity_with_four_signatures = cos_sim(
                    M, np.dot(W, finalactivities)
                )
                # print(cosine_similarity_with_four_signatures)
        else:
            # print("DANG!!! layer_original_distance > original_distance")
            # print("layer_original_distance when layer_original_distance > original_distance:"+str(layer_original_distance))
            # print("original_distance when layer_original_distance > original_distance:"+str(original_distance))
            N = np.dot(W, finalactivities)
            # print("The " + str(i)+"th loop for value N:" + str(N))
            cosine_similarity = cos_sim(M, N)
            kldiv = round(scipy.stats.entropy(M, N), 3)
            correlation, _ = scipy.stats.pearsonr(M, N)
            correlation = round(correlation, 2)
            break

    # get the maximum value of the new Exposure
    maxcoef = np.max(finalactivities)
    idxmaxcoef = list(finalactivities).index(maxcoef)

    finalactivities = np.round(finalactivities)
    # print("The finalactivities after loop but before maxmutation:"+str(finalactivities))

    # We may need to tweak the maximum value of the new exposure to keep the total number of mutation equal to the original mutations in a genome
    if np.sum(finalactivities) != maxmutation:
        finalactivities[idxmaxcoef] = (
            round(finalactivities[idxmaxcoef]) + maxmutation - np.sum(finalactivities)
        )

    # print("The finalactivities after loop and after maxmutation:"+str(type(finalactivities)))
    # import pdb; pdb.set_trace();
    (
        finalactivities,
        original_distance,
        cosine_similarity,
        kldiv,
        correlation,
    ) = weighted_threshold_add_remove(finalactivities, W, M, threshold, ALLSIGIDSGLOBAL)

    finalactivities = np.round(finalactivities)

    # We may need to tweak the maximum value of the new exposure to keep the total number of mutation equal to the original mutations in a genome
    if np.sum(finalactivities) != maxmutation:
        finalactivities[idxmaxcoef] = (
            round(finalactivities[idxmaxcoef]) + maxmutation - np.sum(finalactivities)
        )

    # print("The finalactivities after nnls:"+str(type(finalactivities)))

    if verbose:
        print("\n########################## Final ###########################")
        print(background_sigs)
        print(original_distance)
        print(finalactivities)
    lognote.write(
        "\n#################### Final Composition #################################\n"
    )
    lognote.write(
        "{}\n".format(sub.get_items_from_index(allsigids, selected_signatures))
    )
    lognote.write("L2 Error % {}\n".format(round(original_distance, 2)))
    lognote.write(
        "Cosine Similarity {}\n".format(
            round(cos_sim(M, np.dot(W, finalactivities)), 2)
        )
    )
    # close lognote
    lognote.close()
    # newExposure, newSimilarity = fit_signatures(W[:,list(background_sigs)], M, metric="l2")
    # print(newExposure, newSimilarity)
    return (
        background_sigs,
        finalactivities,
        original_distance,
        cosine_similarity,
        kldiv,
        correlation,
        cosine_similarity_with_four_signatures,
    )
