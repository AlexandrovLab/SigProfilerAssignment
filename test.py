#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12/2/2021
@author: rvangara
"""

import os
import time
import numpy as np
import pandas as pd
from pypdf import PdfWriter, PdfReader
import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerAssignment.DecompositionPlots import PlotDecomposition as sp

# --- Constants and Configuration ---

# Base directory for test data
DIR_INPUT = os.path.join(spa.__path__[0], "data", "tests")
TXT_INPUT_DIR = os.path.join(DIR_INPUT, "txt_input")
VCF_INPUT_DIR = os.path.join(DIR_INPUT, "vcf_input")
SAMPLES_VCF = VCF_INPUT_DIR

# Output directories
OUTPUT_MATRIX_DIR = "TestOutput/output_matrix_example/"
OUTPUT_VCF_DIR = "TestOutput/output_vcf_example/"
OUTPUT_PLOT_DIR = "TestOutput/Results/"

# Common analysis parameters
COMMON_PARAMS = {
    "genome_build": "GRCh37",
    "verbose": False,
}

# Configuration for different analysis types (SBS, DBS, etc.)
ANALYSIS_CONFIGS = {
    "SBS": {
        "signatures": os.path.join(TXT_INPUT_DIR, "SBS96_S3_Signatures.txt"),
        "samples": os.path.join(TXT_INPUT_DIR, "sample_matrix_SBS.txt"),
        "custom_db": os.path.join(TXT_INPUT_DIR, "custom-SBS-list.txt"),
    },
    "DBS": {
        "signatures": os.path.join(TXT_INPUT_DIR, "DBS78_S3_Signatures.txt"),
        "samples": os.path.join(TXT_INPUT_DIR, "sample_matrix_DBS.txt"),
    },
    "ID": {
        "signatures": os.path.join(TXT_INPUT_DIR, "ID83_S3_Signatures.txt"),
        "samples": os.path.join(TXT_INPUT_DIR, "sample_matrix_ID.txt"),
    },
    "SV": {
        "signatures": os.path.join(TXT_INPUT_DIR, "SV32_S3_Signatures.txt"),
        "samples": os.path.join(TXT_INPUT_DIR, "sample_matrix_SV32.txt"),
    },
    "CNV": {
        "signatures": os.path.join(TXT_INPUT_DIR, "CNV48_S3_Signatures.txt"),
        "samples": os.path.join(TXT_INPUT_DIR, "sample_matrix_CNV48.txt"),
    },
}

# Configuration for decomposition plot generation
PLOT_CONFIGS = {
    "SBS96": {
        "mtype": "96",
        "denovo_file": "De_Novo_Solution_Signatures_SBS96.txt",
        "basis_file": "Decomposed_Solution_Signatures_SBS96.txt",
        "denovo_name": "SBS96A",
        "basis_names": ["SBS1", "SBS3", "SBS5", "SBS13", "SBS50", "SBS2"],
        "weights": ["0.94%", "48.72%", "28.44%", "8.42%", "13.48%", "0%"],
    },
    "SBS1536": {
        "mtype": "1536",
        "denovo_file": "De_Novo_Solution_Signatures_SBS1536.txt",
        "basis_file": "Decomposed_Solution_Signatures_SBS1536.txt",
        "denovo_name": "SBS1536A",
        "basis_names": ["SBS1", "SBS2", "SBS5", "SBS13", "SBS15", "SBS18"],
        "weights": ["0.94%", "48.72%", "28.44%", "8.42%", "13.48%", "0%"],
    },
    "SBS288": {
        "mtype": "288",
        "denovo_file": "De_Novo_Solution_Signatures_SBS288.txt",
        "basis_file": "COSMIC_SBS288_Signatures.txt",
        "denovo_name": "SBS288A",
        "basis_names": ["SBS1", "SBS2", "SBS4", "SBS5", "SBS8", "SBS90"],
        "weights": ["0.94%", "48.72%", "28.44%", "8.42%", "13.48%", "0%"],
    },
    "ID83": {
        "mtype": "83",
        "denovo_file": "De_Novo_Solution_Signatures_INDEL.txt",
        "basis_file": "COSMIC_ID83_Signatures.txt",
        "denovo_name": "ID83A",
        "basis_names": ["ID1", "ID2", "ID4", "ID5", "ID8", "ID9"],
        "weights": ["0.94%", "48.72%", "28.44%", "8.42%", "13.48%", "0%"],
    },
    "DBS78": {
        "mtype": "78",
        "denovo_file": "De_Novo_Solution_Signatures_DINUC.txt",
        "basis_file": "COSMIC_DBS78_Signatures.txt",
        "denovo_name": "DBS78A",
        "basis_names": ["DBS1", "DBS2", "DBS4", "DBS5", "DBS8", "DBS9"],
        "weights": ["0.94%", "48.72%", "28.44%", "8.42%", "13.48%", "0%"],
    },
    "CNV48": {
        "mtype": "48",
        "denovo_file": "CNV48_De-Novo_Signatures.txt",
        "basis_file": "COSMIC_CNV48_Signatures.txt",
        "denovo_name": "CNV48A",
        "basis_names": ["CN1", "CN2", "CN9", "CN20", "CN3", "CN4"],
        "weights": ["0.94%", "48.72%", "28.44%", "8.42%", "13.48%", "0%"],
    },
    "SV32": {
        "mtype": "32",
        "denovo_file": "SV32_De-Novo_Signatures.txt",
        "basis_file": "COSMIC_SV32_Signatures.txt",
        "denovo_name": "SV32A",
        "basis_names": ["SV1", "SV2", "SV3", "SV4", "SV5", "SV6", "SV7", "SV9"],
        "weights": ["10%", "10%", "10%", "10%", "10%", "10%", "10%", "20%"],
    },
}


# --- Reusable Functions ---

def run_analysis_tests(mutation_type, config):
    """
    Runs a suite of analysis tests (cosmic, denovo, decompose) for a given
    mutation type.

    Args:
        mutation_type (str): The mutation type (e.g., 'SBS', 'DBS').
        config (dict): A dictionary containing paths to sample and signature files.
    """
    print(f"Running {mutation_type} matrix input tests...")
    samples = config["samples"]
    signatures = config["signatures"]
    
    # Run cosmic_fit
    Analyze.cosmic_fit(
        samples,
        os.path.join(OUTPUT_MATRIX_DIR, f"cosmic_fit_{mutation_type}_matrix_test"),
        signatures=None,
        signature_database=None,
        collapse_to_SBS96=mutation_type == "SBS",
        sample_reconstruction_plots=True,
        cosmic_version=3.4,
        exclude_signature_subgroups=None,
        **COMMON_PARAMS
    )

    # Run denovo_fit (does not take exclude_signature_subgroups)
    Analyze.denovo_fit(
        samples,
        os.path.join(OUTPUT_MATRIX_DIR, f"denovo_fit_{mutation_type}_matrix_test"),
        signatures=signatures,
        **COMMON_PARAMS
    )

    # Run decompose_fit
    Analyze.decompose_fit(
        samples,
        os.path.join(OUTPUT_MATRIX_DIR, f"decompose_fit_{mutation_type}_matrix_test"),
        signatures=signatures,
        new_signature_thresh_hold=0.8,
        exclude_signature_subgroups=None,
        **COMMON_PARAMS
    )

def run_vcf_tests():
    """Runs the analysis tests on VCF input files."""
    print("Running VCF input tests...")
    sbs_signatures = ANALYSIS_CONFIGS["SBS"]["signatures"]

    # cosmic_fit for VCF
    Analyze.cosmic_fit(
        SAMPLES_VCF,
        OUTPUT_VCF_DIR,
        input_type="vcf",
        collapse_to_SBS96=True,
        make_plots=True,
        exclude_signature_subgroups=None,
        **COMMON_PARAMS
    )

    # decompose_fit for VCF
    Analyze.decompose_fit(
        SAMPLES_VCF,
        OUTPUT_VCF_DIR,
        input_type="vcf",
        signatures=sbs_signatures,
        new_signature_thresh_hold=0.8,
        make_plots=True,
        exclude_signature_subgroups=None,
        **COMMON_PARAMS
    )

    # denovo_fit for VCF (does not take exclude_signature_subgroups)
    Analyze.denovo_fit(
        SAMPLES_VCF,
        OUTPUT_VCF_DIR,
        input_type="vcf",
        signatures=sbs_signatures,
        make_plots=True,
        **COMMON_PARAMS
    )

def generate_decomposition_plots(plot_name, config):
    """
    Generates and saves a PDF of decomposition plots for a given configuration.

    Args:
        plot_name (str): The name for the plot set (e.g., 'SBS96').
        config (dict): A dictionary with plot-specific parameters.
    """
    start_time = time.time()
    print(f"Generating {plot_name} decomposition plots...")
    np.random.seed(1234567)
    
    plot_input_dir = os.path.join(spa.__path__[0], "DecompositionPlots", "ExampleSample")
    denovo_mtx = pd.read_csv(os.path.join(plot_input_dir, config["denovo_file"]), sep="\t")
    basis_mtx = pd.read_csv(os.path.join(plot_input_dir, config["basis_file"]), sep="\t")
    
    merger = PdfWriter()
    
    # Loop to generate plots with a decreasing number of basis signatures
    for i in range(len(config["basis_names"]), 0, -1):
        basis_names = config["basis_names"][:i]
        weights = config["weights"][:i]
        
        denovo_cols = ["MutationType", config["denovo_name"]]
        basis_cols = ["MutationType"] + basis_names
        
        nonzero_exposures = np.random.uniform(size=len(basis_names))

        result_pdf = sp.run_PlotDecomposition(
            denovo_mtx[denovo_cols],
            config["denovo_name"],
            basis_mtx[basis_cols],
            basis_names,
            weights,
            nonzero_exposures,
            OUTPUT_PLOT_DIR,
            project="test_run",
            mtype=config["mtype"],
        )
        
        # Add the generated plot to the merged PDF
        result_pdf.seek(0)
        reader = PdfReader(result_pdf)
        for page in reader.pages:
            merger.add_page(page)
            
    # Write the final merged PDF
    output_filename = os.path.join(OUTPUT_PLOT_DIR, f"Result_Decomposition_Plots_{plot_name}.pdf")
    merger.write(output_filename)
    merger.close()
    
    end_time = time.time()
    print(f"Completed {plot_name} plots in {end_time - start_time:.2f} seconds.")

def run_custom_db_test():
    """
    Validates decomposition with a custom signature database.
    """
    print("Running custom database test...")
    sbs_config = ANALYSIS_CONFIGS["SBS"]
    Analyze.decompose_fit(
        samples=sbs_config["samples"],
        output=os.path.join(OUTPUT_MATRIX_DIR, "decompose_fit_custom_db_test"),
        signatures=sbs_config["signatures"],
        signature_database=sbs_config["custom_db"],
        new_signature_thresh_hold=0.8,
        make_plots=True,
        exclude_signature_subgroups=None,
        **COMMON_PARAMS
    )


# --- Main Execution ---

def main():
    """
    Main function to run all tests.
    """
    # Create output directories if they don't exist
    for dir_path in [OUTPUT_MATRIX_DIR, OUTPUT_VCF_DIR, OUTPUT_PLOT_DIR]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    # Run matrix-based analysis tests
    run_custom_db_test()
    for m_type, config in ANALYSIS_CONFIGS.items():
        run_analysis_tests(m_type, config)
    
    # Run VCF-based analysis tests
    run_vcf_tests()

    # Generate all decomposition plots
    for plot_name, config in PLOT_CONFIGS.items():
        generate_decomposition_plots(plot_name, config)

if __name__ == "__main__":
    main()
