#!/usr/bin/env python3

"""
Created: March 10th, 2020
@author: Mark Barnes

PlotDecomposition works with matrix formats SigProfiler SBS-96, SBS-1536, DBS-78,
and ID-83. This program is intended to take two matrices.

(1) Sample matrix - A SigProfiler formatted SBS-96, SBS-1536, DBS-78, or ID-83
matrix.
(2) Basis matrix - A SigProfiler formatted SBS-96, SBS-1536, DBS-78, or ID-83
matrix that is the decomposition of (1).

When running the function 'run_PlotDecomposition' a plot of the decomposition will
be generated and saved to the output folder. Refer to the function below to learn
more about the parameters required to generate the decomposition plot.
"""

import os
import pandas as pd
import numpy as np
import scipy.stats
import sigProfilerPlotting as sigPlt
import shutil
import SigProfilerAssignment
import SigProfilerAssignment.DecompositionPlots
from SigProfilerAssignment.DecompositionPlots import SigProfilerPlottingMatrix as mPlt
from SigProfilerAssignment.DecompositionPlots import PlotDecomposition_SBS96 as spd_96
from SigProfilerAssignment.DecompositionPlots import PlotDecomposition_SBS288 as spd_288
from SigProfilerAssignment.DecompositionPlots import (
    PlotDecomposition_SBS1536 as spd_1536,
)
from SigProfilerAssignment.DecompositionPlots import PlotDecomposition_DBS78 as spd_78
from SigProfilerAssignment.DecompositionPlots import PlotDecomposition_ID83 as spd_83
from SigProfilerAssignment.DecompositionPlots import PlotDecomposition_CNV48 as cnv_48
from SigProfilerAssignment import decompose_subroutines as sub

# imports for working with plots in memory
import io
from PIL import Image
from reportlab.lib.utils import ImageReader
import json
import base64

# Global Variables
SBS_CONTEXTS = ["6", "24", "96", "288", "384", "1536", "6144"]
DBS_CONTEXTS = ["78", "186", "1248", "2976"]
ID_CONTEXTS = ["28", "83", "415"]
CNV_CONTEXTS = ["48"]
MTYPE_OPTIONS = [
    "6",
    "24",
    "96",
    "384",
    "1536",
    "6144",
    "28",
    "83",
    "415",
    "78",
    "186",
    "1248",
    "2976",
]
DECOMPOSITION_PATH = SigProfilerAssignment.DecompositionPlots.__path__[0]
REFSIG_PATH = os.path.join(
    SigProfilerAssignment.__path__[0], "data/Reference_Signatures"
)
TEMPLATE_PATH = os.path.join(DECOMPOSITION_PATH, "CosmicTemplates")


# Remove templates so that they can be rebuilt
def remove_cosmic_templates():
    if not os.path.exists(TEMPLATE_PATH):
        print("No files are installed at: ", TEMPLATE_PATH)
    try:
        shutil.rmtree(TEMPLATE_PATH)
    except OSError as e:
        print("Error: %s : %s" % (TEMPLATE_PATH, e.strerror))


# Create a set of serialized JSON reference signature plots for fast loading
def install_cosmic_plots(
    context_type="96", genome_build="GRCh37", cosmic_version="3.4", exome=False
):
    if not os.path.exists(TEMPLATE_PATH):
        os.mkdir(TEMPLATE_PATH)

    # determine if context is from SBS, ID, DBS, or CNV
    context_type_str = ""
    if context_type in SBS_CONTEXTS:
        context_type_str = "SBS"
        # the mtype of the cosmic reference signatures to be plotted
        cosmic_mtype = "96"
    elif context_type in DBS_CONTEXTS:
        context_type_str = "DBS"
        cosmic_mtype = "78"
    elif context_type in ID_CONTEXTS:
        context_type_str = "ID"
        cosmic_mtype = "83"
    elif context_type in CNV_CONTEXTS:
        context_type_str = "CNV"
        cosmic_mtype = "48"
    else:
        raise ValueError("ERROR: context", context_type, "not in context lists.")

    # if exome is true, append _exome to end of file name
    exome_str = ""
    if exome:
        exome_str = "_exome"

    cosmic_file_name = (
        "COSMIC_v"
        + str(cosmic_version)
        + "_"
        + context_type_str
        + "_"
        + genome_build
        + exome_str
        + ".txt"
    )
    json_file_name = (
        "COSMIC_v"
        + str(cosmic_version)
        + "_"
        + context_type_str
        + "_"
        + genome_build
        + exome_str
        + ".json"
    )

    # ID signatures exome=False, genome_build=GRCh37
    if context_type in ID_CONTEXTS:
        cosmic_file_name = "COSMIC_v" + str(cosmic_version) + "_ID_GRCh37.txt"
        json_file_name = "COSMIC_v" + str(cosmic_version) + "_ID_GRCh37.json"
        genome_build = "GRCh37"
        exome_str = ""

    # Load cosmic plots if they exist
    filename = os.path.join(TEMPLATE_PATH, json_file_name)
    if os.path.exists(filename):
        cosmic_buff_bytes = json.load(open(filename))
        cosmic_buff_plots = {}
        # Read from JSON, decode, and convert to bytesIO
        for tmp_plots in cosmic_buff_bytes.keys():
            cosmic_buff_plots[tmp_plots] = io.BytesIO(
                base64.b64decode(cosmic_buff_bytes[tmp_plots])
            )
        return cosmic_buff_plots

    # Generate cosmic plots if they were not found
    else:
        cosmic_file_path = os.path.join(REFSIG_PATH, genome_build, cosmic_file_name)
        json_file_path = os.path.join(TEMPLATE_PATH, json_file_name)
        print(
            "Generating plots for",
            "COSMIC_v"
            + str(cosmic_version)
            + "_"
            + context_type_str
            + "_"
            + genome_build
            + exome_str,
            "now...",
        )

        # Create the respective plots
        if context_type_str == "SBS":
            cosmic_buff_plots = sigPlt.plotSBS(
                cosmic_file_path,
                "buffer",
                "buffer",
                cosmic_mtype,
                percentage=True,
                savefig_format="PIL_Image",
            )
        elif context_type_str == "DBS":
            cosmic_mtx = pd.read_csv(cosmic_file_path, sep="\t")
            cosmic_buff_plots = sigPlt.plotDBS(
                cosmic_mtx,
                "buffer",
                "buffer",
                cosmic_mtype,
                percentage=True,
                savefig_format="PIL_Image",
            )
        elif context_type_str == "ID":
            cosmic_buff_plots = sigPlt.plotID(
                cosmic_file_path,
                "buffer",
                "buffer",
                cosmic_mtype,
                percentage=True,
                savefig_format="PIL_Image",
            )

        # Process the plots to be stored in JSON file
        cosmic_img_dict = {}
        for tmp_plot in cosmic_buff_plots.keys():
            # Convert to bytesIO and encode to base64
            plot_bytes = io.BytesIO()
            seek_start = cosmic_buff_plots[tmp_plot].seek(0)
            cosmic_buff_plots[tmp_plot].save(plot_bytes, "png")
            encoded = base64.b64encode(plot_bytes.getvalue())
            cosmic_img_dict[tmp_plot] = encoded.decode("ascii")
            plot_bytes.close()

        # JSON output processing
        json_object = json.dumps(cosmic_img_dict)
        with open(json_file_path, "w") as outfile:
            outfile.write(json_object)

        print(
            "Plots for",
            "COSMIC_v"
            + str(cosmic_version)
            + "_"
            + context_type_str
            + "_"
            + genome_build
            + exome_str,
            "have been successfully installed.",
        )
        return cosmic_buff_plots


# Convert array of Images to ImageReader array for ReportLab
def convert_to_imgReaderDict(image_dict):
    imgReader_dict = dict()
    for name in image_dict.keys():
        imgReader_dict[name] = ImageReader(image_dict[name])
    return imgReader_dict


def calculate_similarities(denovo, denovo_name, est_denovo):
    from numpy import inf

    # If matrix is 1536 context, then collapse it to 96 format
    if denovo.shape[0] == 1536:
        index = denovo.iloc[:, 0]
        denovo_tmp = pd.DataFrame(denovo, index=index)
        denovo_tmp = denovo.groupby(denovo_tmp.index.str[1:8]).sum()
        denovo = pd.DataFrame(denovo_tmp)
        denovo = denovo.reset_index()
    elif denovo.shape[0] == 288:
        index = denovo.iloc[:, 0]
        denovo_tmp = pd.DataFrame(denovo, index=index)
        denovo_tmp = denovo.groupby(denovo_tmp.index.str[2:9]).sum()
        denovo = pd.DataFrame(denovo_tmp)
        denovo = denovo.reset_index()

    sample_names = [denovo_name]

    if sample_names is False:
        sample_names = ["None"] * denovo.shape[1]

    cosine_similarity_list = []
    cosine_distance_list = []
    correlation_list = []
    correlation_distance_list = []
    kl_divergence_list = []
    l1_norm_list = []
    l2_norm_list = []
    relative_l1_list = []
    relative_l2_list = []

    p_i = denovo[denovo_name]
    q_i = est_denovo

    cosine_similarity_list.append(round(sub.cos_sim(p_i, q_i), 3))
    cosine_distance_list.append(round(scipy.spatial.distance.cosine(p_i, q_i), 3))
    correlation_list.append(round(scipy.stats.pearsonr(p_i, q_i)[0], 3))
    correlation_distance_list.append(round(1 - scipy.stats.pearsonr(p_i, q_i)[0], 3))
    kl_divergence_list.append(round(scipy.stats.entropy(p_i, q_i), 4))
    l1_norm_list.append(round(np.linalg.norm(p_i - q_i, ord=1), 2))
    relative_l1_list.append(
        round((l1_norm_list[-1] / np.linalg.norm(p_i, ord=1)) * 100, 3)
    )
    l2_norm_list.append(round(np.linalg.norm(p_i - q_i, ord=2), 2))
    relative_l2_list.append(
        round((l2_norm_list[-1] / np.linalg.norm(p_i, ord=2)) * 100, 3)
    )
    kl_divergence_list = np.array(kl_divergence_list)
    kl_divergence_list[kl_divergence_list == inf] = 1000

    similarities_dataframe = pd.DataFrame(
        {
            "Sample Names": sample_names,
            "Cosine Similarity": cosine_similarity_list,
            "Cosine Distance": cosine_distance_list,
            "Correlation Distance": correlation_distance_list,
            "Correlation Coefficient": correlation_list,
            "L1 Norm": l1_norm_list,
            "L1 Norm %": relative_l1_list,
            "L2 Norm": l2_norm_list,
            "L2 Norm %": relative_l2_list,
            "KL Divergence": kl_divergence_list,
        }
    )
    similarities_dataframe = similarities_dataframe.set_index("Sample Names")

    return similarities_dataframe


def genSBS_pngs(denovo_mtx, basis_mtx, output_path, project, mtype, ss_decomp=False):
    denovo_plots = dict()
    basis_plots = dict()
    if mtype == "1536" or mtype == "288":
        denovo_plots = mPlt.plotSBS(
            denovo_mtx, output_path, project, mtype, percentage=True
        )
        if basis_mtx is not None:
            basis_plots = sigPlt.plotSBS(
                basis_mtx,
                output_path,
                project,
                "96",
                percentage=True,
                savefig_format="PIL_Image",
            )
    elif mtype == "96":
        denovo_plots = sigPlt.plotSBS(
            denovo_mtx,
            output_path,
            project,
            mtype,
            percentage=(not ss_decomp),
            savefig_format="PIL_Image",
        )
        if basis_mtx is not None:
            basis_plots = sigPlt.plotSBS(
                basis_mtx,
                output_path,
                project,
                mtype,
                percentage=True,
                savefig_format="PIL_Image",
            )
    return denovo_plots, basis_plots


def genDBS_pngs(denovo_mtx, basis_mtx, output_path, project, mtype):
    denovo_plots = dict()
    basis_plots = dict()
    denovo_plots = sigPlt.plotDBS(
        denovo_mtx,
        output_path,
        project,
        mtype,
        percentage=True,
        savefig_format="PIL_Image",
    )
    if basis_mtx is not None:
        basis_plots = sigPlt.plotDBS(
            basis_mtx,
            output_path,
            project,
            mtype,
            percentage=True,
            savefig_format="PIL_Image",
        )
    return denovo_plots, basis_plots


def genID_pngs(denovo_mtx, basis_mtx, output_path, project, mtype):
    denovo_plots = dict()
    basis_plots = dict()
    denovo_plots = sigPlt.plotID(
        denovo_mtx,
        output_path,
        project,
        mtype,
        percentage=True,
        savefig_format="PIL_Image",
    )
    if basis_mtx is not None:
        basis_plots = sigPlt.plotID(
            basis_mtx,
            output_path,
            project,
            mtype,
            percentage=True,
            savefig_format="PIL_Image",
        )
    return denovo_plots, basis_plots


def genCNV_pngs(denovo_mtx, basis_mtx, output_path, project, mtype):
    denovo_plots = dict()
    basis_plots = dict()
    denovo_plots = sigPlt.plotCNV(
        denovo_mtx,
        output_path,
        project,
        percentage=True,
        aggregate=False,
        read_from_file=False,
        savefig_format="PIL_Image",
    )

    if basis_mtx is not None:
        basis_plots = sigPlt.plotCNV(
            basis_mtx,
            output_path,
            project,
            percentage=True,
            aggregate=False,
            read_from_file=False,
            savefig_format="PIL_Image",
        )
    return denovo_plots, basis_plots


# signames, weights
def gen_sub_plots(denovo_mtx, basis_mtx, output_path, project, mtype, ss_decomp):
    # Make output directory
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if mtype in SBS_CONTEXTS:
        denovo_plots, basis_plots = genSBS_pngs(
            denovo_mtx, basis_mtx, output_path, project, mtype, ss_decomp
        )
        return denovo_plots, basis_plots

    elif mtype in DBS_CONTEXTS:
        denovo_plots, basis_plots = genDBS_pngs(
            denovo_mtx, basis_mtx, output_path, project, mtype
        )
        return denovo_plots, basis_plots

    elif mtype in ID_CONTEXTS:
        denovo_plots, basis_plots = genID_pngs(
            denovo_mtx, basis_mtx, output_path, project, mtype
        )
        return denovo_plots, basis_plots

    elif mtype in CNV_CONTEXTS:
        denovo_plots, basis_plots = genCNV_pngs(
            denovo_mtx, basis_mtx, output_path, project, mtype
        )
        return denovo_plots, basis_plots

    else:
        print("ERROR: mtype is " + mtype + " and is not yet supported.")


# generate the plot for the reconstruction
def gen_reconstructed_png_percent(
    denovo_name, basis_mtx, basis_names, weights, output_path, project, mtype
):
    reconstruction_plot = dict()
    mut_col = basis_mtx.iloc[:, 0]

    recon_plot = basis_mtx[basis_names[0]] * float(weights[0].strip("%")) / 100

    for i in range(1, len(weights)):
        recon_plot = recon_plot + basis_mtx[basis_names[i]] * (
            float(weights[i].strip("%")) / 100
        )

    recon_plot = pd.Series(recon_plot, name=denovo_name)
    reconstruction_mtx = pd.concat([mut_col, recon_plot], axis=1)
    if mtype in SBS_CONTEXTS:
        if mtype == "1536" or mtype == "288":
            reconstruction_plot = sigPlt.plotSBS(
                reconstruction_mtx,
                output_path,
                "reconstruction_" + project,
                "96",
                percentage=True,
                savefig_format="PIL_Image",
            )
        else:
            reconstruction_plot = sigPlt.plotSBS(
                reconstruction_mtx,
                output_path,
                "reconstruction_" + project,
                mtype,
                percentage=True,
                savefig_format="PIL_Image",
            )
    elif mtype in DBS_CONTEXTS:
        reconstruction_plot = sigPlt.plotDBS(
            reconstruction_mtx,
            output_path,
            "reconstruction_" + project,
            mtype,
            percentage=True,
            savefig_format="PIL_Image",
        )
    elif mtype in ID_CONTEXTS:
        reconstruction_plot = sigPlt.plotID(
            reconstruction_mtx,
            output_path,
            "reconstruction_" + project,
            mtype,
            percentage=True,
            savefig_format="PIL_Image",
        )
    elif mtype in CNV_CONTEXTS:
        reconstruction_plot = sigPlt.plotCNV(
            reconstruction_mtx,
            output_path,
            "reconstruction_" + project,
            percentage=True,
            aggregate=False,
            read_from_file=False,
            savefig_format="PIL_Image",
        )
    else:
        print("ERROR: mtype is " + mtype + " and is not yet supported.")

    return reconstruction_mtx, reconstruction_plot


# generate the plot for the reconstruction
def gen_reconstructed_png_numerical(
    denovo_mtx,
    denovo_name,
    basis_mtx,
    basis_names,
    weights,
    output_path,
    project,
    mtype,
):
    sample_tmb = denovo_mtx[denovo_name].sum()
    reconstruction_plot = dict()
    mut_col = basis_mtx.iloc[:, 0]

    recon_plot = (
        basis_mtx[basis_names[0]] * float(weights[0].strip("%")) / 100 * sample_tmb
    )
    for i in range(1, len(weights)):
        recon_plot = (
            recon_plot
            + basis_mtx[basis_names[i]]
            * (float(weights[i].strip("%")) / 100)
            * sample_tmb
        )

    recon_plot = pd.Series(recon_plot, name=denovo_name)
    reconstruction_mtx = pd.concat([mut_col, recon_plot], axis=1)
    reconstruction_mtx = reconstruction_mtx.round(0)
    if mtype in SBS_CONTEXTS:
        if mtype == "1536" or mtype == "288":
            reconstruction_plot = sigPlt.plotSBS(
                reconstruction_mtx,
                output_path,
                "reconstruction_" + project,
                "96",
                percentage=False,
            )
        else:
            reconstruction_plot = sigPlt.plotSBS(
                reconstruction_mtx,
                output_path,
                "reconstruction_" + project,
                mtype,
                percentage=False,
                savefig_format="PIL_Image",
            )
    elif mtype in DBS_CONTEXTS:
        reconstruction_plot = sigPlt.plotDBS(
            reconstruction_mtx,
            output_path,
            "reconstruction_" + project,
            mtype,
            percentage=False,
            savefig_format="PIL_Image",
        )
    elif mtype in ID_CONTEXTS:
        reconstruction_plot = sigPlt.plotID(
            reconstruction_mtx,
            output_path,
            "reconstruction_" + project,
            mtype,
            percentage=False,
            savefig_format="PIL_Image",
        )
    elif mtype in CNV_CONTEXTS:
        reconstruction_plot = sigPlt.plotCNV(
            reconstruction_mtx,
            output_path,
            "reconstruction_" + project,
            percentage=True,
            aggregate=False,
            read_from_file=False,
            savefig_format="PIL_Image",
        )
    else:
        print("ERROR: mtype is " + mtype + " and is not yet supported.")

    return reconstruction_mtx, reconstruction_plot


def gen_decomposition(
    denovo_name,
    basis_names,
    weights,
    output_path,
    project,
    mtype,
    denovo_plots_dict,
    basis_plots_dict,
    reconstruction_plot_dict,
    reconstruction=False,
    statistics=None,
    cosmic_version=None,
    custom_text=None,
):
    """
    Generate the correct plot based on mtype.

    Parameters:
    ----------
    denovo_name: 				(String) 			Name of denovo signature
    basis_names: 				(List of Strings) 	Names of basis signatures
    weights:					(List of Strings) 	Percentile contribution for each basis signature
    output_path: 				(String) 			Path to existing output directory
    project: 					(String) 			Project name appended to file names
    mtype: 						(String) 			The context 'mtype_options' has valid values
    denovo_plots_dict			(Dictionary)		Signatures are keys, ByteIO plots are values
    basis_plots_dict			(Dictionary)		Signatures are keys, ByteIO plots are values
    reconstruction_plot_dict	(Dictionary)		Signatures are keys, ByteIO plots are values
    reconstruction: 			(Boolean) 			True to generate plot w/ reconstruction
    statistics: 				(Pandas Dataframe) 	Output from calculate_similarities()
    """

    if mtype == "6":
        print("Need to add support for SBS6 Decomposition")
    elif mtype == "24":
        print("Need to add support for SBS24 Decomposition")
    elif mtype == "96":
        byte_plot = spd_96.gen_decomposition(
            denovo_name,
            basis_names,
            weights,
            output_path,
            project,
            denovo_plots_dict,
            basis_plots_dict,
            reconstruction_plot_dict,
            reconstruction,
            statistics,
            cosmic_version,
            custom_text,
        )
        return byte_plot
    elif mtype == "288":
        byte_plot = spd_288.gen_decomposition(
            denovo_name,
            basis_names,
            weights,
            output_path,
            project,
            denovo_plots_dict,
            basis_plots_dict,
            reconstruction_plot_dict,
            reconstruction,
            statistics,
            cosmic_version,
            custom_text,
        )
        return byte_plot
    elif mtype == "384":
        print("Need to add support for SBS24 Decomposition")
    elif mtype == "1536":
        byte_plot = spd_1536.gen_decomposition(
            denovo_name,
            basis_names,
            weights,
            output_path,
            project,
            denovo_plots_dict,
            basis_plots_dict,
            reconstruction_plot_dict,
            reconstruction,
            statistics,
            cosmic_version,
            custom_text,
        )
        return byte_plot
    elif mtype == "6144":
        print("Need to add support for SBS6144 Decomposition")
    elif mtype == "28":
        print("Need to add support for ID28 Decomposition")
    elif mtype == "83":
        byte_plot = spd_83.gen_decomposition(
            denovo_name,
            basis_names,
            weights,
            output_path,
            project,
            denovo_plots_dict,
            basis_plots_dict,
            reconstruction_plot_dict,
            reconstruction,
            statistics,
            cosmic_version,
            custom_text,
        )
        return byte_plot
    elif mtype == "415":
        print("Need to add support for ID415 Decomposition")
    elif mtype == "78":
        byte_plot = spd_78.gen_decomposition(
            denovo_name,
            basis_names,
            weights,
            output_path,
            project,
            denovo_plots_dict,
            basis_plots_dict,
            reconstruction_plot_dict,
            reconstruction,
            statistics,
            cosmic_version,
            custom_text,
        )
        return byte_plot
    elif mtype == "186":
        print("Need to add support for DBS186 Decomposition")
    elif mtype == "1248":
        print("Need to add support for DBS1248 Decomposition")
    elif mtype == "2976":
        print("Need to add support for DBS2976 Decomposition")
    elif mtype == "48":
        byte_plot = cnv_48.gen_decomposition(
            denovo_name,
            basis_names,
            weights,
            output_path,
            project,
            denovo_plots_dict,
            basis_plots_dict,
            reconstruction_plot_dict,
            reconstruction,
            statistics,
            cosmic_version,
            custom_text,
        )
        return byte_plot


def run_PlotDecomposition(
    denovo_mtx,
    denovo_name,
    basis_mtx,
    basis_names,
    weights,
    nonzero_exposures,
    output_path,
    project,
    mtype,
    cosmic_version="3.4",
    genome_build="GRCh37",
    exome=False,
    custom_text=None,
):
    """
    Generates a decomposition plot of the denovo_mtx using the basis_mtx.

    Parameters:
    ----------

    denovo_mtx: Pandas Dataframe. This format represents the catalog of mutations seperated by tab.

    denovo_name: String. The name of the one sample in denovo_mtx to decompose.

    basis_mtx: Pandas Dataframe. This format represents the catalog of mutations seperated by tab.

    basis_names: List of Strings. The names of the samples in denovo_mtx that
    the denovo_name sample from denovo_mtx is decomposed into.
    ie. basis_names=["SBS1", "SBS5", "SBS15", "SBS20"]

    weights: List of Strings. The percentile weight corresponding to each basis
    in basis_names. Refer to example function call below for more detail.
    ie. weights=["11.58%", "42.38%", "16.46%", "29.58%"]

    output_path: String. Path to where to store the output.

    project: String. This string is appended to the file name output.

    mtype: String. The context of the data. Valid values include: "96", "1536","78", and "83".

    cosmic_version: String. The version of signatures being used.

    custom_text: String. A custom message displayed on decomposition plot.

    Returns:
    -------
    None.
    """
    # Create the denovo plots and load basis plots
    if mtype != "48":
        denovo_plots_dict = gen_sub_plots(
            denovo_mtx, None, output_path, project, mtype, ss_decomp=False
        )
        denovo_plots_dict = denovo_plots_dict[0]
    else:
        # cnv basis plots need to be generated and not loaded
        denovo_plots_dict, basis_plots_dict = gen_sub_plots(
            denovo_mtx, basis_mtx, output_path, project, mtype, ss_decomp=False
        )
    # Create the matrix and plot for the reconstructed matrix
    reconstructed_mtx, reconstruction_plot_dict = gen_reconstructed_png_percent(
        denovo_name, basis_mtx, basis_names, weights, output_path, project, mtype
    )
    # Create a subset matrix with the present signatures
    present_sigs = np.array(basis_mtx[basis_names])
    reconstructed_mtx = np.dot(present_sigs, nonzero_exposures)
    # Convert dictionary of bytes to dictionary of images
    denovo_plots_dict = convert_to_imgReaderDict(denovo_plots_dict)
    # Load in the COSMIC plots
    if mtype != "48":
        basis_plots_dict = install_cosmic_plots(
            context_type=mtype,
            genome_build=genome_build,
            cosmic_version=cosmic_version,
            exome=exome,
        )
        basis_plots_dict = {key: basis_plots_dict[key] for key in basis_names}
    basis_plots_dict = convert_to_imgReaderDict(basis_plots_dict)
    # Generate the reconstruction plot
    reconstruction_plot_dict = convert_to_imgReaderDict(reconstruction_plot_dict)
    # Get the reconstruction statistics
    statistics = calculate_similarities(denovo_mtx, denovo_name, reconstructed_mtx)
    # Return the decomposition plot as a byte array
    byte_plot = gen_decomposition(
        denovo_name,
        basis_names,
        weights,
        output_path,
        project,
        mtype,
        denovo_plots_dict,
        basis_plots_dict,
        reconstruction_plot_dict,
        reconstruction=True,
        statistics=statistics,
        cosmic_version=cosmic_version,
        custom_text=custom_text,
    )
    # Clear the plotting memory
    sigPlt.clear_plotting_memory()

    return byte_plot


# context="96", genome_build="GRCh37", cosmic_version="3.4", exome=False
def run_PlotSSDecomposition(
    denovo_mtx,
    denovo_name,
    basis_mtx,
    basis_names,
    weights,
    output_path,
    project,
    context_type,
    genome_build="GRCh37",
    cosmic_version="3.4",
    custom_text=None,
    exome=False,
):
    """
    Generates a reconstruction of a sample given a set of signatures.

    Parameters:
    ----------

    samples: Pandas Dataframe. Samples and their channel counts.
                    df_allgenomes =  pd.DataFrame(allgenomes, columns=allcolnames)

    denovo_mtx: Pandas Dataframe. This format represents the catalog of mutations seperated by tab.

    denovo_name: String. The name of the one sample in denovo_mtx to decompose.

    basis_mtx: Pandas Dataframe. This format represents the catalog of mutations seperated by tab.

    basis_names: List of Strings. The names of the samples in denovo_mtx that
    the denovo_name sample from denovo_mtx is decomposed into.
    ie. basis_names=["SBS1", "SBS5", "SBS15", "SBS20"]

    weights: List of Strings. The percentile weight corresponding to each basis
    in basis_names. Refer to example function call below for more detail.
    ie. weights=["11.58%", "42.38%", "16.46%", "29.58%"]

    output_path: String. Path to where to store the output.

    project: String. This string is appended to the file name output.

    context_type: String. The context of the data. Valid values include: "96", "1536","78", and "83".

    genome_build: String. The genome being used.

    cosmic_version: String. The version of signatures being used.

    custom_text: String. A custom message displayed on decomposition plot.

    exome: Boolean. True if using exome COSMIC signatures, and False if not.

    Returns:
    -------
    None.
    """

    # Create the denovo plots
    denovo_plots_dict = gen_sub_plots(
        denovo_mtx, None, output_path, project, context_type, ss_decomp=True
    )
    denovo_plots_dict = denovo_plots_dict[0]
    # Load in the COSMIC plots
    basis_plots_dict = install_cosmic_plots(
        context_type=context_type,
        genome_build=genome_build,
        cosmic_version=cosmic_version,
        exome=exome,
    )

    # Create reconstructed matrix and plot
    reconstructed_mtx, reconstruction_plot_dict = gen_reconstructed_png_numerical(
        denovo_mtx,
        denovo_name,
        basis_mtx,
        basis_names,
        weights,
        output_path,
        project,
        context_type,
    )

    denovo_plots_dict = convert_to_imgReaderDict(denovo_plots_dict)
    # subset basis_plots_dict to only the plots used
    basis_plots_dict = {key: basis_plots_dict[key] for key in basis_names}
    basis_plots_dict = convert_to_imgReaderDict(basis_plots_dict)

    reconstruction_plot_dict = convert_to_imgReaderDict(reconstruction_plot_dict)

    statistics = calculate_similarities(
        denovo_mtx, denovo_name, reconstructed_mtx[denovo_name]
    )
    byte_plot = gen_decomposition(
        denovo_name,
        basis_names,
        weights,
        output_path,
        project,
        context_type,
        denovo_plots_dict,
        basis_plots_dict,
        reconstruction_plot_dict,
        reconstruction=True,
        statistics=statistics,
        cosmic_version=cosmic_version,
        custom_text=custom_text,
    )
    # Clear the plotting memory
    sigPlt.clear_plotting_memory()

    return byte_plot
