#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 19 12:21:06 2019

@author: mishugeb
"""

# from SigProfilerExtractor import subroutines as sub

from cmath import cos
import datetime
import platform

# from torch import sign
from SigProfilerAssignment import decompose_subroutines as sub
from SigProfilerAssignment.DecompositionPlots import PlotDecomposition as plot_decomp
import SigProfilerAssignment

import numpy as np
import pandas as pd
import SigProfilerMatrixGenerator
from SigProfilerMatrixGenerator.scripts import (
    SigProfilerMatrixGeneratorFunc as datadump,
)
from SigProfilerMatrixGenerator.scripts import CNVMatrixGenerator as scna
from sigProfilerPlotting import sigProfilerPlotting as sigPlot
import sigProfilerPlotting
import os, sys
from PyPDF2 import PdfMerger
import fitz
import time


def convert_PDF_to_PNG(input_file_name, output_directory, page_names):
    pdf_doc = fitz.open(input_file_name)
    zoom = 3
    magnify = fitz.Matrix(zoom, zoom)

    if pdf_doc.page_count != len(page_names):
        raise ValueError(
            "Error: The number of samples and number of plots do not match."
        )
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for sample_name, page in zip(page_names, pdf_doc):
        pix = page.get_pixmap(matrix=magnify)
        out_file_name = os.path.join(output_directory, sample_name + ".png")
        pix.save(out_file_name)


# Create sample reconstruction plots
def generate_sample_reconstruction(
    cosmic_sigs,
    samples_input,
    activities,
    output_dir,
    recon_output_types,
    execution_parameters,
):
    # recon_output_types does the following:
    #       pdf     = pdf generation
    #       both    = pdf generation + png conversion
    #       png     = pdf generation + png conversion + pdf removal
    project = "test_run"
    mtype = "96"

    final_pdf = PdfMerger()
    samples = samples_input.copy(deep=True)
    samples.reset_index(inplace=True)
    for sample_name in samples.columns[1:]:
        # basis names and their corresponding weights
        subset = activities[activities["Samples"] == sample_name]
        subset = subset.loc[:, (subset != 0).any(axis=0)]
        basis_names = (
            subset[subset["Samples"].str.contains(sample_name)].columns[1:].tolist()
        )
        recon_tmb = subset.sum(axis=1)
        weights = []
        for i in range(len(basis_names)):
            weights.append(
                str(
                    float(
                        "{0:.6f}".format(
                            100 * int(subset[basis_names[i]]) / int(recon_tmb)
                        )
                    )
                )
                + "%"
            )

        names_copy = basis_names.copy()
        names_copy.insert(0, cosmic_sigs.columns[0])
        result = plot_decomp.run_PlotSSDecomposition(
            samples[[samples.columns[0], sample_name]],
            sample_name,
            cosmic_sigs[names_copy],
            basis_names,
            weights,
            output_dir,
            project,
            mtype,
            genome_build=execution_parameters["reference_genome"],
            cosmic_version=str(execution_parameters["cosmic_version"]),
            exome=execution_parameters["exome"],
        )
        final_pdf.append(result)

    pdf_output_path = os.path.join(
        output_dir, "Reconstructed_Sample_Plots_" + str(mtype) + ".pdf"
    )
    web_png_path = os.path.join(output_dir, "WebPNGs")

    # write out pdf
    final_pdf.write(pdf_output_path)

    # write out pngs if png (then remove the pdf), or keep both
    if recon_output_types.lower() == "png" or recon_output_types.lower() == "both":
        convert_PDF_to_PNG(pdf_output_path, web_png_path, samples.columns[1:])
        if recon_output_types.lower() == "png":
            # remove the PDF
            os.remove(pdf_output_path)

    return pdf_output_path


def record_parameters(sysdata, execution_parameters, start_time):
    sysdata.write("\n--------------EXECUTION PARAMETERS--------------\n")
    sysdata.write("INPUT DATA\n")
    sysdata.write("\tinput_type: {}\n".format(execution_parameters["input_type"]))
    sysdata.write("\toutput: {}\n".format(execution_parameters["output"]))
    if isinstance(execution_parameters["samples"], str):
        sysdata.write("\tsamples: {}\n".format(execution_parameters["samples"]))
    else:
        sysdata.write("\tsamples: {}\n".format(type(execution_parameters["samples"])))
    sysdata.write(
        "\treference_genome: {}\n".format(execution_parameters["reference_genome"])
    )
    sysdata.write("\tcontext_types: {}\n".format(execution_parameters["context_type"]))
    sysdata.write("\texome: {}\n".format(execution_parameters["exome"]))

    sysdata.write("COSMIC MATCH\n")
    sysdata.write(
        "\tcosmic_version: {}\n".format(execution_parameters["cosmic_version"])
    )
    sysdata.write(
        "\tnnls_add_penalty: {}\n".format(execution_parameters["nnls_add_penalty"])
    )
    sysdata.write(
        "\tnnls_remove_penalty: {}\n".format(
            execution_parameters["nnls_remove_penalty"]
        )
    )
    sysdata.write(
        "\tinitial_remove_penalty: {}\n".format(
            execution_parameters["initial_remove_penalty"]
        )
    )
    sysdata.write(
        "\tde_novo_fit_penalty: {}\n".format(
            execution_parameters["de_novo_fit_penalty"]
        )
    )
    sysdata.write(
        "\texport_probabilities: {}\n".format(
            execution_parameters["export_probabilities"]
        )
    )
    sysdata.write(
        "\tcollapse_to_SBS96: {}\n".format(execution_parameters["collapse_to_SBS96"])
    )
    sysdata.write(
        "\tdenovo_refit_option: {}\n".format(
            execution_parameters["denovo_refit_option"]
        )
    )
    sysdata.write(
        "\tdecompose_fit_option: {}\n".format(
            execution_parameters["decompose_fit_option"]
        )
    )
    sysdata.write(
        "\tcosmic_fit_option: {}\n".format(execution_parameters["cosmic_fit_option"])
    )
    sysdata.write("\n-------Analysis Progress------- \n")
    sysdata.write("[{}] Analysis started: \n".format(str(start_time).split(".")[0]))


def spa_analyze(
    samples,
    output,
    input_type="matrix",
    context_type="96",
    signatures=None,
    signature_database=None,
    decompose_fit_option=True,
    denovo_refit_option=True,
    cosmic_fit_option=True,
    nnls_add_penalty=0.05,
    nnls_remove_penalty=0.01,
    initial_remove_penalty=0.05,
    de_novo_fit_penalty=0.02,
    genome_build="GRCh37",
    cosmic_version=3.4,
    make_plots=True,
    collapse_to_SBS96=True,
    connected_sigs=True,
    verbose=False,
    devopts=None,
    new_signature_thresh_hold=0.8,
    exclude_signature_subgroups=None,
    exome=False,
    export_probabilities=True,
    export_probabilities_per_mutation=False,
    sample_reconstruction_plots=None,
    make_metadata=True,
):
    """
    Decomposes the De Novo Signatures into COSMIC Signatures and assigns COSMIC signatures into samples.

    Parameters:

        signatures: A string. Path to a  tab delimited file that contains the signaure table where the rows are mutation types and colunms are signature IDs.
        activities: A string. Path to a tab delimilted file that contains the activity table where the rows are sample IDs and colunms are signature IDs.
        samples: A string. Path to a tab delimilted file that contains the activity table where the rows are mutation types and colunms are sample IDs.
        output: A string. Path to the output folder.
        genome_build = A string. The reference genome build. List of supported genomes: "GRCh37", "GRCh38", "mm9", "mm10" and "rn6". The default value is "GRCh37". If the selected genome is not in the supported list, the default genome will be used.
        verbose = Boolean. Prints statements. Default value is False.
        exome = Boolean. Defines if the exome renormalized signatures will be used. The default value is False.

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
        layer_directory1 = output + "/De_Novo_Solution"
        layer_directory2 = output + "/Decompose_Solution"
        layer_directory3 = output + "/Assignment_Solution"
    else:
        layer_directory1 = devopts["denovo_outpath"]
        layer_directory2 = devopts["decompose_outpath"]
        layer_directory3 = devopts["Assignment_outpath"]

    if (
        denovo_refit_option == True or decompose_fit_option == True
    ) and signatures is None:
        raise Exception(
            "If denovo_refit or decompose_fit is True, signatures cannot be empty"
        )

    if input_type == "vcf":
        project_name = "Input_vcffiles"
        vcf_context = context_type
        data = datadump.SigProfilerMatrixGeneratorFunc(
            project_name,
            genome_build,
            samples,
            exome=exome,
            bed_file=None,
            chrom_based=False,
            plot=False,
            gs=False,
        )
        genomes = data[vcf_context]

    elif input_type.split(":")[0].lower() == "seg":
        cnv_file_type = input_type.split(":")[1].upper()
        context_type = "CNV48"
        genomes = scna.generateCNVMatrix(
            cnv_file_type, samples, cnv_file_type, os.path.dirname(samples) + "/"
        )
        genomes = genomes.set_index("MutationType")

    elif input_type == "matrix":
        try:
            genomes = pd.read_csv(samples, sep="\t", index_col=0)
        except:
            genomes = samples
    else:
        sys.exit("Invalid input_type specified")

    mutation_type = str(genomes.shape[0])
    m = mutation_type

    # Re-indexing the input matrix file by using process_input function from SigProfilePlotting
    genomes = sigPlot.process_input(genomes, m)

    m_for_subgroups = ""
    if m == "96" or m == "288" or m == "1536":
        m_for_subgroups = "SBS"
    if m == "78":
        m_for_subgroups = "DBS"
    if m == "83":
        m_for_subgroups = "ID"

    default_subgroups_dict = {
        "MMR_deficiency_signatures": False,
        "POL_deficiency_signatures": False,
        "HR_deficiency_signatures": False,
        "BER_deficiency_signatures": False,
        "Chemotherapy_signatures": False,
        "Immunosuppressants_signatures": False,
        "Treatment_signatures": False,
        "APOBEC_signatures": False,
        "Tobacco_signatures": False,
        "UV_signatures": False,
        "AA_signatures": False,
        "Colibactin_signatures": False,
        "Artifact_signatures": False,
        "Lymphoid_signatures": False,
    }

    default_subgroups_siglists = {
        "MMR_deficiency_signatures": {
            "SBS": ["6", "14", "15", "20", "21", "26", "44"],
            "DBS": ["7", "10"],
            "ID": ["7"],
        },
        "POL_deficiency_signatures": {
            "SBS": ["10a", "10b", "10c", "10d", "28"],
            "DBS": ["3"],
            "ID": [],
        },
        "HR_deficiency_signatures": {"SBS": ["3"], "DBS": ["13"], "ID": ["6"]},
        "BER_deficiency_signatures": {"SBS": ["30", "36"], "DBS": [], "ID": []},
        "Chemotherapy_signatures": {
            "SBS": ["11", "25", "31", "35", "86", "87", "90", "99"],
            "DBS": ["5"],
            "ID": [],
        },
        "Immunosuppressants_signatures": {"SBS": ["32"], "DBS": [], "ID": []},
        "Treatment_signatures": {
            "SBS": ["11", "25", "31", "32", "35", "86", "87", "90", "99"],
            "DBS": ["5"],
            "ID": [],
        },
        "APOBEC_signatures": {"SBS": ["2", "13"], "DBS": [], "ID": []},
        "Tobacco_signatures": {"SBS": ["4", "29", "92"], "DBS": ["2"], "ID": ["3"]},
        "UV_signatures": {
            "SBS": ["7a", "7b", "7c", "7d", "38"],
            "DBS": ["1"],
            "ID": ["13"],
        },
        "AA_signatures": {"SBS": ["22", "22a", "22b"], "DBS": ["20"], "ID": ["23"]},
        "Colibactin_signatures": {"SBS": ["88"], "DBS": [], "ID": ["18"]},
        "Artifact_signatures": {
            "SBS": [
                "27",
                "43",
                "45",
                "46",
                "47",
                "48",
                "49",
                "50",
                "51",
                "52",
                "53",
                "54",
                "55",
                "56",
                "57",
                "58",
                "59",
                "60",
                "95",
            ],
            "DBS": ["14"],
            "ID": [],
        },
        "Lymphoid_signatures": {"SBS": ["9", "84", "85"], "DBS": [], "ID": []},
    }

    signature_subgroups_dict = default_subgroups_dict.copy()
    if exclude_signature_subgroups == None:
        pass
    else:
        if type(exclude_signature_subgroups) is not list:
            sys.exit(
                "exclude_signature_subgroups input should be a list of appropriate flags, please refer to documentation."
            )
        else:
            for key in default_subgroups_dict:
                if key in exclude_signature_subgroups:
                    signature_subgroups_dict[key] = True

    sig_exclusion_list = []
    if exclude_signature_subgroups == None:
        sig_exclusion_list = []
    else:
        for key in signature_subgroups_dict:
            if signature_subgroups_dict[key]:
                sig_exclusion_list.append(
                    default_subgroups_siglists[key][m_for_subgroups]
                )

    sig_exclusion_list = [item for sublist in sig_exclusion_list for item in sublist]

    try:
        if not os.path.exists(output):
            os.makedirs(output)
    except:
        print("The {} folder could not be created".format("output"))

    start_time = datetime.datetime.now()  # start time

    execution_parameters = {
        "input_type": input_type,
        "context_type": context_type,
        "signatures": signatures,
        "signature_database": signature_database,
        "denovo_refit_option": denovo_refit_option,
        "decompose_fit_option": decompose_fit_option,
        "cosmic_fit_option": cosmic_fit_option,
        "output": output,
        "samples": samples,
        "reference_genome": genome_build,
        "cosmic_version": cosmic_version,
        "context_type": context_type,
        "exome": exome,
        "nnls_add_penalty": nnls_add_penalty,
        "nnls_remove_penalty": nnls_remove_penalty,
        "initial_remove_penalty": initial_remove_penalty,
        "de_novo_fit_penalty": de_novo_fit_penalty,
        "collapse_to_SBS96": collapse_to_SBS96,
        "export_probabilities": export_probabilities,
        "make_plots": make_plots,
    }

    if make_metadata:
        # create JOB_METADATA_SPA
        sysdata = open(os.path.join(output, "JOB_METADATA_SPA.txt"), "w")
        sysdata.write("THIS FILE CONTAINS THE METADATA ABOUT SYSTEM AND RUNTIME\n\n\n")
        sysdata.write("-------System Info-------\n")
        sysdata.write(
            "Operating System Name: "
            + platform.uname()[0]
            + "\n"
            + "Nodename: "
            + platform.uname()[1]
            + "\n"
            + "Release: "
            + platform.uname()[2]
            + "\n"
            + "Version: "
            + platform.uname()[3]
            + "\n"
        )
        sysdata.write("\n-------Python and Package Versions------- \n")
        sysdata.write(
            "Python Version: "
            + str(platform.sys.version_info.major)
            + "."
            + str(platform.sys.version_info.minor)
            + "."
            + str(platform.sys.version_info.micro)
            + "\n"
        )
        sysdata.write(
            "SigProfilerPlotting Version: " + sigProfilerPlotting.__version__ + "\n"
        )
        sysdata.write(
            "SigProfilerMatrixGenerator Version: "
            + SigProfilerMatrixGenerator.__version__
            + "\n"
        )
        sysdata.write(
            "SigProfilerAssignment Version: " + SigProfilerAssignment.__version__ + "\n"
        )
        sysdata.write("Pandas version: " + pd.__version__ + "\n")
        sysdata.write("Numpy version: " + np.__version__ + "\n")
        record_parameters(sysdata, execution_parameters, start_time)
        sysdata.close()
    # Add sequence parameter to control the tmbplot y-axis scale
    if exome == True:
        sequence = "exome"
    else:
        sequence = "genome"

        #################
        # Denovo refiting #
        #################

    if denovo_refit_option == True:
        try:
            processAvg = pd.read_csv(signatures, sep="\t", index_col=0)
        except:
            try:
                processAvg = signatures
            except:
                sys.exit(
                    "Error in formatting of input signatures, Pass a text file of signatures in the format of denovo signatures"
                )

        if devopts == None:
            listOfSignatures = processAvg.columns
            index = genomes.index
            colnames = genomes.columns
        else:
            listOfSignatures = devopts["listOfSignatures"]
            index = devopts["index"]
            colnames = devopts["colnames"]
            genomes = genomes.set_index(index)
            genomes.columns = colnames
            # genomes = genomes.rename_axis("Mutation Types", axis="columns")

        # creating list of mutational type to sync with the vcf type input
        if mutation_type == "78":
            mutation_context = "DBS78"
        elif mutation_type == "83":
            mutation_context = "ID83"
        elif mutation_type == "48":
            mutation_context = "CNV48"
        elif mutation_type == "32":
            mutation_context = "SV32"
        else:
            mutation_context = "SBS" + mutation_type
        try:
            allsigids = processAvg.columns.to_list()
        except:
            allsigids = list(listOfSignatures)
        processAvg = np.array(processAvg)
        signature_names = sub.make_letter_ids(
            idlenth=processAvg.shape[1], mtype=mutation_context
        )

        # create the folder for the final solution/ De Novo Solution
        exposureAvg_dummy = (
            pd.DataFrame(
                np.random.rand(processAvg.shape[1], genomes.shape[1]),
                index=listOfSignatures,
                columns=colnames.to_list(),
            )
            .transpose()
            .rename_axis("Samples")
        )
        exposureAvg = exposureAvg_dummy
        exposureAvg.columns = signature_names

        refit_denovo_signatures = True
        init_rem_denovo = 0.0

        # make the texts for signature plotting
        # layer_directory1 = output+"/De_Novo_Solution"
        try:
            if not os.path.exists(layer_directory1):
                os.makedirs(layer_directory1)
        except:
            print("The {} folder could not be created".format("De_Novo_Solution"))

        listOfSignatures = sub.make_letter_ids(
            idlenth=processAvg.shape[1], mtype=mutation_context
        )
        genomes = pd.DataFrame(genomes)
        denovo_exposureAvg = np.array(exposureAvg.T)
        print("\n De Novo Fitting .....")
        current_time_start = datetime.datetime.now()
        if make_metadata:
            with open(os.path.join(output, "JOB_METADATA_SPA.txt"), "a") as sysdata:
                sysdata.write("\n De Novo Fitting .....")

        # exposureAvg = sub.make_final_solution(processAvg, genomes, listOfSignatures, layer_directory1, mutation_type, index,\
        # colnames,denovo_exposureAvg  = denovo_exposureAvg, add_penalty=nnls_add_penalty, remove_penalty=nnls_remove_penalty, initial_remove_penalty=initial_remove_penalty, de_novo_fit_penalty=de_novo_fit_penalty, connected_sigs=connected_sigs, refit_denovo_signatures=refit_denovo_signatures)
        #
        #######
        attribution = {}
        for i in allsigids:
            attribution[i] = [i]
        # only for SBS96
        if mutation_type == "96" or mutation_type == "288" or mutation_type == "1536":
            background_sigs = sub.get_indeces(list(allsigids), ["SBS1", "SBS5"])
            # add connected signatures
            # different_signatures = ss.add_connected_sigs(different_signatures, list(signames))
        # for other contexts
        else:
            background_sigs = []
        exposureAvg_dummy = (
            pd.DataFrame(
                np.random.rand(processAvg.shape[1], genomes.shape[1]),
                index=allsigids,
                columns=colnames.to_list(),
            )
            .transpose()
            .rename_axis("Samples")
        )

        if devopts == None:
            exposureAvg = sub.make_final_solution(
                processAvg,
                genomes,
                allsigids,
                layer_directory1,
                mutation_type,
                index,
                colnames,
                cosmic_sigs=True,
                attribution=attribution,
                denovo_exposureAvg=exposureAvg_dummy,
                sequence=sequence,
                background_sigs=background_sigs,
                verbose=verbose,
                genome_build=genome_build,
                add_penalty=nnls_add_penalty,
                remove_penalty=nnls_remove_penalty,
                initial_remove_penalty=init_rem_denovo,
                connected_sigs=connected_sigs,
                refit_denovo_signatures=False,
                make_plots=make_plots,
                export_probabilities=export_probabilities,
                export_probabilities_per_mutation=export_probabilities_per_mutation,
                samples=samples,
                input_type=input_type,
                denovo_refit_option=denovo_refit_option,
                exome=exome,
            )

        else:
            signature_stabilities = devopts["signature_stabilities"]
            signature_total_mutations = devopts["signature_total_mutations"]
            signature_stats = devopts["signature_stats"]
            sequence = devopts["sequence"]
            processSTE = devopts["processSTE"]

            exposureAvg = sub.make_final_solution(
                processAvg,
                genomes,
                allsigids,
                layer_directory1,
                mutation_type,
                index,
                colnames,
                cosmic_sigs=True,
                attribution=attribution,
                denovo_exposureAvg=exposureAvg_dummy,
                sequence=sequence,
                background_sigs=background_sigs,
                verbose=verbose,
                genome_build=genome_build,
                signature_total_mutations=signature_total_mutations,
                add_penalty=nnls_add_penalty,
                remove_penalty=nnls_remove_penalty,
                process_std_error=processSTE,
                signature_stabilities=signature_stabilities,
                initial_remove_penalty=init_rem_denovo,
                connected_sigs=connected_sigs,
                refit_denovo_signatures=True,
                export_probabilities=export_probabilities,
                export_probabilities_per_mutation=export_probabilities_per_mutation,
                samples=samples,
                input_type=input_type,
                denovo_refit_option=denovo_refit_option,
                exome=exome,
            )

        if make_metadata:
            with open(os.path.join(output, "JOB_METADATA_SPA.txt"), "a") as sysdata:
                current_time_end = datetime.datetime.now()
                sysdata.write(
                    f"\n Finished Denovo fitting! \nExecution time:{str(current_time_end-current_time_start)}\n"
                )

                #################
                # Decomposition
                #################
    if decompose_fit_option == True:
        try:
            processAvg = pd.read_csv(signatures, sep="\t", index_col=0)
        except:
            try:
                processAvg = signatures
            except:
                sys.exit(
                    "Error in formatting of input signatures, Pass a text file of signatures in the format of denovo signatures"
                )

        if devopts == None:
            listOfSignatures = processAvg.columns
            index = genomes.index
            colnames = genomes.columns
            make_decomposition_plots = make_plots
        else:
            listOfSignatures = devopts["listOfSignatures"]
            index = devopts["index"]
            colnames = devopts["colnames"]
            genomes = genomes.set_index(index)
            genomes.columns = colnames
            make_decomposition_plots = devopts["make_decomposition_plots"]
            # genomes = genomes.rename_axis("Mutation Types", axis="columns")

        # creating list of mutational type to sync with the vcf type input
        if mutation_type == "78":
            mutation_context = "DBS78"
        elif mutation_type == "83":
            mutation_context = "ID83"
        elif mutation_type == "48":
            mutation_context = "CNV48"
        elif mutation_type == "32":
            mutation_context = "SV32"
        else:
            mutation_context = "SBS" + mutation_type
        try:
            allsigids = processAvg.columns.to_list()
        except:
            allsigids = list(listOfSignatures)
        processAvg = np.array(processAvg)
        signature_names = sub.make_letter_ids(
            idlenth=processAvg.shape[1], mtype=mutation_context
        )

        exposureAvg_dummy = (
            pd.DataFrame(
                np.random.rand(processAvg.shape[1], genomes.shape[1]),
                index=listOfSignatures,
                columns=colnames.to_list(),
            )
            .transpose()
            .rename_axis("Samples")
        )
        exposureAvg = exposureAvg_dummy
        exposureAvg.columns = signature_names

        #############################
        # layer_directory2 = output+"/Decompose_Solution"
        if isinstance(processAvg, pd.DataFrame):
            pass
        else:
            originalProcessAvg = pd.DataFrame(
                processAvg, index=index, columns=listOfSignatures
            )
        try:
            if not os.path.exists(layer_directory2):
                os.makedirs(layer_directory2)
        except:
            print("The {} folder could not be created".format("Decomposed_Solution"))

        if (
            processAvg.shape[0] == 1536 and collapse_to_SBS96 == True
        ):  # collapse the 1596 context into 96 only for the deocmposition
            processAvg = pd.DataFrame(processAvg, index=index)
            processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
            genomes = genomes.groupby(genomes.index.str[1:8]).sum()
            index = genomes.index
            processAvg = np.array(processAvg)

        if (
            processAvg.shape[0] == 288 and collapse_to_SBS96 == True
        ):  # collapse the 288 context into 96 only for the deocmposition
            processAvg = pd.DataFrame(processAvg, index=index)
            processAvg = processAvg.groupby(processAvg.index.str[2:9]).sum()
            genomes = pd.DataFrame(genomes, index=index)
            genomes = genomes.groupby(genomes.index.str[2:9]).sum()
            index = genomes.index
            processAvg = np.array(processAvg)

        print("\n Decomposing De Novo Signatures  .....")
        current_time_start = datetime.datetime.now()
        if make_metadata:
            with open(os.path.join(output, "JOB_METADATA_SPA.txt"), "a") as sysdata:
                sysdata.write("\n Decomposing De Novo Signatures  .....")
        final_signatures = sub.signature_decomposition(
            processAvg,
            mutation_type,
            layer_directory2,
            genome_build=genome_build,
            cosmic_version=cosmic_version,
            signature_database=signature_database,
            mutation_context=mutation_context,
            add_penalty=0.05,
            connected_sigs=connected_sigs,
            remove_penalty=0.01,
            make_decomposition_plots=make_decomposition_plots,
            originalProcessAvg=originalProcessAvg,
            new_signature_thresh_hold=new_signature_thresh_hold,
            sig_exclusion_list=sig_exclusion_list,
            exome=exome,
            m_for_subgroups=m_for_subgroups,
        )
        # final_signatures = sub.signature_decomposition(processAvg, m, layer_directory2, genome_build=genome_build)
        # extract the global signatures and new signatures from the final_signatures dictionary
        globalsigs = final_signatures["globalsigs"]
        globalsigs = np.array(globalsigs)
        newsigs = final_signatures["newsigs"]
        processAvg = np.hstack([globalsigs, newsigs])
        allsigids = final_signatures["globalsigids"] + final_signatures["newsigids"]
        attribution = final_signatures["dictionary"]
        background_sigs = final_signatures["background_sigs"]
        index = genomes.index
        colnames = genomes.columns

        print("\n Assigning decomposed signature")
        if make_metadata:
            with open(os.path.join(output, "JOB_METADATA_SPA.txt"), "a") as sysdata:
                sysdata.write(
                    "\n [{}] Assigning decomposed signature \n".format(
                        str(start_time).split(".")[0]
                    )
                )

        result = sub.make_final_solution(
            processAvg,
            genomes,
            allsigids,
            layer_directory2,
            mutation_type,
            index,
            colnames,
            cosmic_sigs=True,
            attribution=attribution,
            denovo_exposureAvg=exposureAvg,
            sequence=sequence,
            background_sigs=background_sigs,
            verbose=verbose,
            genome_build=genome_build,
            add_penalty=nnls_add_penalty,
            remove_penalty=nnls_remove_penalty,
            initial_remove_penalty=initial_remove_penalty,
            connected_sigs=connected_sigs,
            collapse_to_SBS96=collapse_to_SBS96,
            refit_denovo_signatures=False,
            make_plots=make_plots,
            export_probabilities=export_probabilities,
            export_probabilities_per_mutation=export_probabilities_per_mutation,
            samples=samples,
            input_type=input_type,
            denovo_refit_option=denovo_refit_option,
            exome=exome,
        )

        if make_metadata:
            with open(os.path.join(output, "JOB_METADATA_SPA.txt"), "a") as sysdata:
                current_time_end = datetime.datetime.now()
                sysdata.write(
                    f"\n Finished Decompose fitting! \nExecution time:{str(current_time_end-current_time_start)}\n"
                )
                #################
                # Cosmic Fitting
                #################

    if cosmic_fit_option == True:
        try:
            if not os.path.exists(layer_directory3):
                os.makedirs(layer_directory3)
        except:
            print("The {} folder could not be created".format("Assignment_Solution"))
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
        if (genomes.sum() == 0).sum() > 0:
            print("Removing samples with zero TMB ...... ")
            genomes = genomes.loc[:, (genomes != 0).any(axis=0)]
            colnames = genomes.columns

        if (
            genomes.shape[0] == 1536 and collapse_to_SBS96 == True
        ):  # collapse the 1596 context into 96 only for the deocmposition
            # processAvg = pd.DataFrame(processAvg, index=index)
            # processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
            genomes = pd.DataFrame(genomes, index=index)
            genomes = genomes.groupby(genomes.index.str[1:8]).sum()
            index = genomes.index
            # processAvg = np.array(processAvg)

        if (
            genomes.shape[0] == 288 and collapse_to_SBS96 == True
        ):  # collapse the 288 context into 96 only for the deocmposition
            # processAvg = pd.DataFrame(processAvg, index=index)
            # processAvg = processAvg.groupby(processAvg.index.str[2:9]).sum()
            genomes = pd.DataFrame(genomes, index=index)
            genomes = genomes.groupby(genomes.index.str[2:9]).sum()
            index = genomes.index

        if signature_database == None:
            processAvg = sub.getProcessAvg(
                genomes,
                genome_build=genome_build,
                cosmic_version=cosmic_version,
                exome=exome,
            )[0]
            # for sample reconstruction plots
            cosmic_sig_ref = processAvg.copy(deep=True)
            cosmic_sig_ref.reset_index(inplace=True)
        else:
            try:
                processAvg = pd.read_csv(signature_database, sep="\t", index_col=0)
            except:
                sys.exit(
                    "Something is wrong with the format of signature database, Pass a text file of signatures in the format of COSMIC sig database"
                )

        # collapse the 1596 context into 96 only for the decomposition
        if processAvg.shape[0] == 1536 and collapse_to_SBS96 == True:
            processAvg = pd.DataFrame(processAvg, index=index)
            processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
            genomes = pd.DataFrame(genomes, index=index)
            genomes = genomes.groupby(genomes.index.str[1:8]).sum()
            index = genomes.index

        if (
            processAvg.shape[0] == 288 and collapse_to_SBS96 == True
        ):  # collapse the 288 context into 96 only for the deocmposition
            processAvg = pd.DataFrame(processAvg, index=index)
            processAvg = processAvg.groupby(processAvg.index.str[2:9]).sum()
            genomes = pd.DataFrame(genomes, index=index)
            genomes = genomes.groupby(genomes.index.str[2:9]).sum()
            index = genomes.index

        # processAvg is sigdatabase: remove sigs corresponding to exclusion rules.
        sig_exclusion_list = [m_for_subgroups + items for items in sig_exclusion_list]
        if sig_exclusion_list:
            print(
                "The following signatures are excluded: "
                + " ".join(str(item) for item in sig_exclusion_list)
            )
        # #
        processAvg.drop(sig_exclusion_list, axis=1, inplace=True, errors="ignore")

        # processAvg= originalProcessAvg
        # index = genomes.index
        # colnames = genomes.columns
        allsigids = processAvg.columns.to_list()
        processAvg = processAvg.values
        attribution = {}
        for i in allsigids:
            attribution[i] = [i]
        # only for SBS96
        if mutation_type == "96" or mutation_type == "288" or mutation_type == "1536":
            background_sigs = sub.get_indeces(list(allsigids), ["SBS1", "SBS5"])
            # add connected signatures
            # different_signatures = ss.add_connected_sigs(different_signatures, list(signames))
        # for other contexts
        else:
            background_sigs = []
        exposureAvg_dummy = (
            pd.DataFrame(
                np.random.rand(processAvg.shape[1], genomes.shape[1]),
                index=allsigids,
                columns=colnames.to_list(),
            )
            .transpose()
            .rename_axis("Samples")
        )
        print("Assigning COSMIC sigs or Signature Database ...... ")
        current_time_start = datetime.datetime.now()
        if make_metadata:
            with open(os.path.join(output, "JOB_METADATA_SPA.txt"), "a") as sysdata:
                sysdata.write("\n Assigning COSMIC sigs or Signature Database ...... ")
        if processAvg.shape[0] != 96:
            if genomes.shape[0] == processAvg.shape[0] and collapse_to_SBS96 == True:
                sys.exit(
                    'Signatures Database and Samples are of same context type and is not equal to 96. please rerun by setting the flag "collapse_to_SBS96 = False "'
                )

        sub.make_final_solution(
            processAvg,
            genomes,
            allsigids,
            layer_directory3,
            mutation_type,
            index,
            colnames,
            cosmic_sigs=True,
            attribution=attribution,
            denovo_exposureAvg=exposureAvg_dummy,
            sequence=sequence,
            background_sigs=background_sigs,
            verbose=verbose,
            genome_build=genome_build,
            add_penalty=nnls_add_penalty,
            remove_penalty=nnls_remove_penalty,
            initial_remove_penalty=initial_remove_penalty,
            connected_sigs=connected_sigs,
            collapse_to_SBS96=collapse_to_SBS96,
            refit_denovo_signatures=False,
            make_plots=make_plots,
            export_probabilities=export_probabilities,
            export_probabilities_per_mutation=export_probabilities_per_mutation,
            samples=samples,
            input_type=input_type,
            denovo_refit_option=denovo_refit_option,
            exome=exome,
        )
        if make_metadata:
            with open(os.path.join(output, "JOB_METADATA_SPA.txt"), "a") as sysdata:
                current_time_end = datetime.datetime.now()
                sysdata.write(
                    f"\n Finished Cosmic fitting! \nExecution time:{str(current_time_end-current_time_start)}\n"
                )

    recon_output_types = ["png", "pdf", "both"]
    # Generate sample reconstruction plots
    if (
        sample_reconstruction_plots in recon_output_types
        and mutation_type == "96"
        and signature_database is None
    ):
        ss_recon_odir = os.path.join(
            layer_directory3, "Activities", "SampleReconstruction"
        )
        if not os.path.exists(ss_recon_odir):
            os.makedirs(ss_recon_odir)
        acts_path = os.path.join(
            layer_directory3, "Activities", "Assignment_Solution_Activities.txt"
        )
        cosmic_activities = pd.read_csv(acts_path, sep="\t")
        # Create sample reconstruction plots
        generate_sample_reconstruction(
            cosmic_sig_ref,
            genomes,
            cosmic_activities,
            ss_recon_odir,
            sample_reconstruction_plots,
            execution_parameters,
        )

    # Complete JOB_METADATA_SPA
    if make_metadata:
        sysdata = open(os.path.join(output, "JOB_METADATA_SPA.txt"), "a")
        end_time = datetime.datetime.now()
        sysdata.write("\n[{}] Analysis ended: \n".format(str(end_time).split(".")[0]))
        sysdata.write("\n-------Job Status------- \n")
        sysdata.write(
            "Assignment of mutational signatures completed successfully! \nTotal execution time: "
            + str(end_time - start_time).split(".")[0]
            + " \nResults can be found in: "
            + " "
            + output
            + " "
            + " folder"
        )
        sysdata.close()
        print(
            "\n\n \nYour Job Is Successfully Completed! Thank You For Using SigProfilerAssignment.\n "
        )
