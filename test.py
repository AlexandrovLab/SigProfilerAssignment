#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12/2/2021
@author: rvangara
"""

import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerAssignment.DecompositionPlots import PlotDecomposition as sp
from PyPDF2 import PdfMerger
import numpy as np
import pandas as pd
import time
import os

DIR_INPUT = spa.__path__[0] + "/data/tests/"
OUTPUT_MATRIX = "TestOutput/output_matrix_example/"
OUTPUT_VCF = "TestOutput/output_vcf_example/"
SIGNATURES_MATRIX_SBS = os.path.join(DIR_INPUT, "txt_input/SBS96_S3_Signatures.txt")
SIGNATURES_MATRIX_DBS = os.path.join(DIR_INPUT, "txt_input/DBS78_S3_Signatures.txt")
SIGNATURES_MATRIX_ID = os.path.join(DIR_INPUT, "txt_input/ID83_S3_Signatures.txt")
SIGNATURES_MATRIX_CNV = os.path.join(DIR_INPUT, "txt_input/CNV48_S3_Signatures.txt")
SIGNATURES_MATRIX_SV = os.path.join(DIR_INPUT, "txt_input/SV32_S3_Signatures.txt")
SAMPLE_MATRIX_SBS = os.path.join(DIR_INPUT, "txt_input/sample_matrix_SBS.txt")
SAMPLE_MATRIX_DBS = os.path.join(DIR_INPUT, "txt_input/sample_matrix_DBS.txt")
SAMPLE_MATRIX_ID = os.path.join(DIR_INPUT, "txt_input/sample_matrix_ID.txt")
SAMPLE_MATRIX_CNV = os.path.join(DIR_INPUT, "txt_input/sample_matrix_CNV48.txt")
SAMPLE_MATRIX_SV = os.path.join(DIR_INPUT, "txt_input/sample_matrix_SV32.txt")
SAMPLES_VCF = os.path.join(DIR_INPUT, "vcf_input/")

exclude_signature_subgroups = None


########### SBS Tests ###########
def cosmic_fit_SBS_matrix_test():
    Analyze.cosmic_fit(
        SAMPLE_MATRIX_SBS,
        os.path.join(OUTPUT_MATRIX, "cosmic_fit_SBS_matrix_test"),
        signatures=None,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        collapse_to_SBS96=True,
        exclude_signature_subgroups=exclude_signature_subgroups,
        sample_reconstruction_plots="both",
        cosmic_version=3.4,
    )


def denovo_fit_SBS_matrix_test():
    Analyze.denovo_fit(
        SAMPLE_MATRIX_SBS,
        os.path.join(OUTPUT_MATRIX, "denovo_fit_SBS_matrix_test"),
        signatures=SIGNATURES_MATRIX_SBS,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
    )


def decompose_fit_SBS_matrix_test():
    Analyze.decompose_fit(
        SAMPLE_MATRIX_SBS,
        os.path.join(OUTPUT_MATRIX, "decompose_fit_SBS_matrix_test"),
        signatures=SIGNATURES_MATRIX_SBS,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        new_signature_thresh_hold=0.8,
        exclude_signature_subgroups=exclude_signature_subgroups,
    )


########### DBS Tests ###########
def cosmic_fit_DBS_matrix_test():
    Analyze.cosmic_fit(
        SAMPLE_MATRIX_DBS,
        os.path.join(OUTPUT_MATRIX, "cosmic_fit_DBS_matrix_test"),
        signatures=None,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        collapse_to_SBS96=False,
        exclude_signature_subgroups=exclude_signature_subgroups,
        sample_reconstruction_plots=True,
        cosmic_version=3.4,
    )


def denovo_fit_DBS_matrix_test():
    Analyze.denovo_fit(
        SAMPLE_MATRIX_DBS,
        os.path.join(OUTPUT_MATRIX, "denovo_fit_DBS_matrix_test"),
        signatures=SIGNATURES_MATRIX_DBS,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
    )


def decompose_fit_DBS_matrix_test():
    Analyze.decompose_fit(
        SAMPLE_MATRIX_DBS,
        os.path.join(OUTPUT_MATRIX, "decompose_fit_DBS_matrix_test"),
        signatures=SIGNATURES_MATRIX_DBS,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        new_signature_thresh_hold=0.8,
        exclude_signature_subgroups=exclude_signature_subgroups,
    )


########### ID Tests ###########
def decompose_fit_ID_matrix_test():
    Analyze.decompose_fit(
        SAMPLE_MATRIX_ID,
        os.path.join(OUTPUT_MATRIX, "decompose_fit_ID_matrix_test"),
        signatures=SIGNATURES_MATRIX_ID,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        new_signature_thresh_hold=0.8,
        exclude_signature_subgroups=exclude_signature_subgroups,
    )


def cosmic_fit_ID_matrix_test():
    Analyze.cosmic_fit(
        SAMPLE_MATRIX_ID,
        os.path.join(OUTPUT_MATRIX, "cosmic_fit_ID_matrix_test"),
        signatures=None,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        collapse_to_SBS96=False,
        exclude_signature_subgroups=exclude_signature_subgroups,
        sample_reconstruction_plots=True,
        cosmic_version=3.4,
    )


def denovo_fit_ID_matrix_test():
    Analyze.denovo_fit(
        SAMPLE_MATRIX_ID,
        os.path.join(OUTPUT_MATRIX, "denovo_fit_ID_matrix_test"),
        signatures=SIGNATURES_MATRIX_ID,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
    )


########### SV Tests ###########
def decompose_fit_SV_matrix_test():
    Analyze.decompose_fit(
        SAMPLE_MATRIX_SV,
        os.path.join(OUTPUT_MATRIX, "decompose_fit_SV_matrix_test"),
        signatures=SIGNATURES_MATRIX_SV,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        new_signature_thresh_hold=0.8,
        exclude_signature_subgroups=exclude_signature_subgroups,
    )


def cosmic_fit_SV_matrix_test():
    Analyze.cosmic_fit(
        SAMPLE_MATRIX_SV,
        os.path.join(OUTPUT_MATRIX, "cosmic_fit_SV_matrix_test"),
        signatures=None,
        signature_database=None,
        genome_build="GRCh38",
        verbose=False,
        collapse_to_SBS96=False,
        exclude_signature_subgroups=exclude_signature_subgroups,
        sample_reconstruction_plots=True,
        cosmic_version=3.4,
    )


def denovo_fit_SV_matrix_test():
    Analyze.denovo_fit(
        SAMPLE_MATRIX_SV,
        os.path.join(OUTPUT_MATRIX, "denovo_fit_SV_matrix_test"),
        signatures=SIGNATURES_MATRIX_SV,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
    )


########### CNV Tests ###########
def decompose_fit_CNV_matrix_test():
    Analyze.decompose_fit(
        SAMPLE_MATRIX_CNV,
        os.path.join(OUTPUT_MATRIX, "decompose_fit_CNV_matrix_test"),
        signatures=SIGNATURES_MATRIX_CNV,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        new_signature_thresh_hold=0.8,
        exclude_signature_subgroups=exclude_signature_subgroups,
    )


def cosmic_fit_CNV_matrix_test():
    Analyze.cosmic_fit(
        SAMPLE_MATRIX_CNV,
        os.path.join(OUTPUT_MATRIX, "cosmic_fit_CNV_matrix_test"),
        signatures=None,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        collapse_to_SBS96=False,
        exclude_signature_subgroups=exclude_signature_subgroups,
        sample_reconstruction_plots=True,
        cosmic_version=3.4,
    )


def denovo_fit_CNV_matrix_test():
    Analyze.denovo_fit(
        SAMPLE_MATRIX_CNV,
        os.path.join(OUTPUT_MATRIX, "denovo_fit_CNV_matrix_test"),
        signatures=SIGNATURES_MATRIX_CNV,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
    )


########### VCF Tests ###########
def cosmic_fit_vcf_test():
    Analyze.cosmic_fit(
        SAMPLES_VCF,
        OUTPUT_VCF,
        input_type="vcf",
        signatures=None,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        collapse_to_SBS96=True,
        exclude_signature_subgroups=exclude_signature_subgroups,
        make_plots=True,
    )


def decompose_fit_vcf_test():
    Analyze.decompose_fit(
        SAMPLES_VCF,
        OUTPUT_VCF,
        input_type="vcf",
        signatures=SIGNATURES_MATRIX_SBS,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        new_signature_thresh_hold=0.8,
        exclude_signature_subgroups=exclude_signature_subgroups,
        make_plots=True,
    )


def denovo_fit_vcf_test():
    Analyze.denovo_fit(
        SAMPLES_VCF,
        OUTPUT_VCF,
        signatures=SIGNATURES_MATRIX_SBS,
        signature_database=None,
        genome_build="GRCh37",
        verbose=False,
        make_plots=True,
        input_type="vcf",
    )


################################################
#           Decomposition Plot Tests           #
################################################
def gen_SBS96():
    np.random.seed(1234567)
    s = time.time()
    merger = PdfMerger()
    file1 = "SigProfilerAssignment/DecompositionPlots/ExampleSample/De_Novo_Solution_Signatures_SBS96.txt"
    file2 = "SigProfilerAssignment/DecompositionPlots/ExampleSample/Decomposed_Solution_Signatures_SBS96.txt"
    denovo_mtx = pd.read_csv(file1, sep="\t")
    basis_mtx = pd.read_csv(file2, sep="\t")
    output_path = "TestOutput/Results/"
    project = "test_run"
    mtype = "96"

    denovo_name = "SBS96A"
    basis_names = ["SBS1", "SBS3", "SBS5", "SBS13", "SBS50", "SBS2"]
    weights = ["0.94%", "48.72%", "28.44%", "8.42%", "13.48%", "0%"]
    denovo_cols = ["MutationType", "SBS96A"]
    basis_cols = basis_names.copy()
    basis_cols.insert(0, "MutationType")
    nonzero_exposures = np.random.uniform(size=len(basis_names))
    result = sp.run_PlotDecomposition(
        denovo_mtx[denovo_cols],
        denovo_name,
        basis_mtx[basis_cols],
        basis_names,
        weights,
        nonzero_exposures,
        output_path,
        project,
        mtype,
    )
    merger.append(result)

    for ind in range(5, 0, -1):
        basis_names = basis_names[:ind]
        weights = weights[:ind]
        denovo_cols = ["MutationType", "SBS96A"]
        basis_cols = basis_names.copy()
        basis_cols.insert(0, "MutationType")
        nonzero_exposures = np.random.uniform(size=len(basis_names))
        result = sp.run_PlotDecomposition(
            denovo_mtx[denovo_cols],
            denovo_name,
            basis_mtx[basis_cols],
            basis_names,
            weights,
            nonzero_exposures,
            output_path,
            project,
            mtype,
        )
        merger.append(result)

    merger.write(os.path.join(output_path, "Result_Decomposition_Plots_SBS96.pdf"))
    return time.time() - s


def gen_SBS1536():
    np.random.seed(1234567)
    s = time.time()
    merger = PdfMerger()
    file1 = "SigProfilerAssignment/DecompositionPlots/ExampleSample/De_Novo_Solution_Signatures_SBS1536.txt"
    file2 = "SigProfilerAssignment/DecompositionPlots/ExampleSample/Decomposed_Solution_Signatures_SBS1536.txt"
    denovo_mtx = pd.read_csv(file1, sep="\t")
    basis_mtx = pd.read_csv(file2, sep="\t")
    output_path = "TestOutput/Results/"
    project = "test_run"
    mtype = "1536"

    denovo_name = "SBS1536A"
    basis_names = ["SBS1", "SBS2", "SBS5", "SBS13", "SBS15", "SBS18"]
    weights = ["0.94%", "48.72%", "28.44%", "8.42%", "13.48%", "0%"]
    denovo_cols = ["MutationType", "SBS1536A"]
    basis_cols = basis_names.copy()
    basis_cols.insert(0, "MutationType")
    nonzero_exposures = np.random.uniform(size=len(basis_names))
    result = sp.run_PlotDecomposition(
        denovo_mtx[denovo_cols],
        denovo_name,
        basis_mtx[basis_cols],
        basis_names,
        weights,
        nonzero_exposures,
        output_path,
        project,
        mtype,
    )
    # sp.run_PlotDecomposition(denovo_mtx, basis_names, weights, output_path, project, mtype, True, statistics, "COSMICv3-GRCh37", "This is where a custom message would go.")
    merger.append(result)

    for ind in range(5, 0, -1):
        basis_names = basis_names[:ind]
        weights = weights[:ind]
        denovo_cols = ["MutationType", "SBS1536A"]
        basis_cols = basis_names.copy()
        basis_cols.insert(0, "MutationType")
        nonzero_exposures = np.random.uniform(size=len(basis_names))
        result = sp.run_PlotDecomposition(
            denovo_mtx[denovo_cols],
            denovo_name,
            basis_mtx[basis_cols],
            basis_names,
            weights,
            nonzero_exposures,
            output_path,
            project,
            mtype,
        )
        merger.append(result)

    merger.write(os.path.join(output_path, "Result_Decomposition_Plots_SBS1536.pdf"))
    return time.time() - s


def gen_SBS288():
    np.random.seed(1234567)
    s = time.time()
    merger = PdfMerger()
    file1 = "SigProfilerAssignment/DecompositionPlots/ExampleSample/De_Novo_Solution_Signatures_SBS288.txt"
    file2 = "SigProfilerAssignment/DecompositionPlots/ExampleSample/COSMIC_SBS288_Signatures.txt"
    denovo_mtx = pd.read_csv(file1, sep="\t")
    basis_mtx = pd.read_csv(file2, sep="\t")
    output_path = "TestOutput/Results/"
    project = "test_run"
    mtype = "288"

    denovo_name = "SBS288A"
    basis_names = ["SBS1", "SBS2", "SBS4", "SBS5", "SBS8", "SBS90"]
    weights = ["0.94%", "48.72%", "28.44%", "8.42%", "13.48%", "0%"]
    denovo_cols = ["MutationType", "SBS288A"]
    basis_cols = basis_names.copy()
    basis_cols.insert(0, "MutationType")
    nonzero_exposures = np.random.uniform(size=len(basis_names))
    result = sp.run_PlotDecomposition(
        denovo_mtx[denovo_cols],
        denovo_name,
        basis_mtx[basis_cols],
        basis_names,
        weights,
        nonzero_exposures,
        output_path,
        project,
        mtype,
    )
    merger.append(result)

    for ind in range(5, 0, -1):
        basis_names = basis_names[:ind]
        weights = weights[:ind]
        denovo_cols = ["MutationType", "SBS288A"]
        basis_cols = basis_names.copy()
        basis_cols.insert(0, "MutationType")
        nonzero_exposures = np.random.uniform(size=len(basis_names))
        result = sp.run_PlotDecomposition(
            denovo_mtx[denovo_cols],
            denovo_name,
            basis_mtx[basis_cols],
            basis_names,
            weights,
            nonzero_exposures,
            output_path,
            project,
            mtype,
        )
        merger.append(result)

    merger.write(os.path.join(output_path, "Result_Decomposition_Plots_SBS288.pdf"))
    return time.time() - s


def gen_ID83():
    np.random.seed(1234567)
    s = time.time()
    merger = PdfMerger()
    file1 = "SigProfilerAssignment/DecompositionPlots/ExampleSample/De_Novo_Solution_Signatures_INDEL.txt"
    file2 = "SigProfilerAssignment/DecompositionPlots/ExampleSample/COSMIC_ID83_Signatures.txt"
    denovo_mtx = pd.read_csv(file1, sep="\t")
    basis_mtx = pd.read_csv(file2, sep="\t")
    output_path = "TestOutput/Results/"
    project = "test_run"
    mtype = "83"

    denovo_name = "ID83A"
    basis_names = ["ID1", "ID2", "ID4", "ID5", "ID8", "ID9"]
    weights = ["0.94%", "48.72%", "28.44%", "8.42%", "13.48%", "0%"]
    denovo_cols = ["MutationType", "ID83A"]
    basis_cols = basis_names.copy()
    basis_cols.insert(0, "MutationType")
    nonzero_exposures = np.random.uniform(size=len(basis_names))
    result = sp.run_PlotDecomposition(
        denovo_mtx[denovo_cols],
        denovo_name,
        basis_mtx[basis_cols],
        basis_names,
        weights,
        nonzero_exposures,
        output_path,
        project,
        mtype,
    )
    merger.append(result)

    for ind in range(5, 0, -1):
        basis_names = basis_names[:ind]
        weights = weights[:ind]
        denovo_cols = ["MutationType", "ID83A"]
        basis_cols = basis_names.copy()
        basis_cols.insert(0, "MutationType")
        nonzero_exposures = np.random.uniform(size=len(basis_names))
        result = sp.run_PlotDecomposition(
            denovo_mtx[denovo_cols],
            denovo_name,
            basis_mtx[basis_cols],
            basis_names,
            weights,
            nonzero_exposures,
            output_path,
            project,
            mtype,
        )
        merger.append(result)

    merger.write(os.path.join(output_path, "Result_Decomposition_Plots_ID83.pdf"))
    return time.time() - s


def gen_DBS78():
    np.random.seed(1234567)
    s = time.time()
    merger = PdfMerger()
    file1 = "SigProfilerAssignment/DecompositionPlots/ExampleSample/De_Novo_Solution_Signatures_DINUC.txt"
    file2 = "SigProfilerAssignment/DecompositionPlots/ExampleSample/COSMIC_DBS78_Signatures.txt"
    denovo_mtx = pd.read_csv(file1, sep="\t")
    basis_mtx = pd.read_csv(file2, sep="\t")
    output_path = "TestOutput/Results/"
    project = "test_run"
    mtype = "78"

    denovo_name = "DBS78A"
    basis_names = ["DBS1", "DBS2", "DBS4", "DBS5", "DBS8", "DBS9"]
    weights = ["0.94%", "48.72%", "28.44%", "8.42%", "13.48%", "0%"]
    denovo_cols = ["MutationType", "DBS78A"]
    basis_cols = basis_names.copy()
    basis_cols.insert(0, "MutationType")
    nonzero_exposures = np.random.uniform(size=len(basis_names))
    result = sp.run_PlotDecomposition(
        denovo_mtx[denovo_cols],
        denovo_name,
        basis_mtx[basis_cols],
        basis_names,
        weights,
        nonzero_exposures,
        output_path,
        project,
        mtype,
    )
    merger.append(result)

    for ind in range(5, 0, -1):
        basis_names = basis_names[:ind]
        weights = weights[:ind]
        denovo_cols = ["MutationType", "DBS78A"]
        basis_cols = basis_names.copy()
        basis_cols.insert(0, "MutationType")
        nonzero_exposures = np.random.uniform(size=len(basis_names))
        result = sp.run_PlotDecomposition(
            denovo_mtx[denovo_cols],
            denovo_name,
            basis_mtx[basis_cols],
            basis_names,
            weights,
            nonzero_exposures,
            output_path,
            project,
            mtype,
        )
        merger.append(result)

    merger.write(os.path.join(output_path, "Result_Decomposition_Plots_DBS78.pdf"))
    return time.time() - s


def gen_CNV48():
    np.random.seed(1234567)
    s = time.time()
    merger = PdfMerger()
    file1 = "SigProfilerAssignment/DecompositionPlots/ExampleSample/CNV48_De-Novo_Signatures.txt"
    file2 = "SigProfilerAssignment/DecompositionPlots/ExampleSample/COSMIC_CNV48_Signatures.txt"
    denovo_mtx = pd.read_csv(file1, sep="\t")
    basis_mtx = pd.read_csv(file2, sep="\t")
    output_path = "TestOutput/Results/"
    project = "test_run"
    mtype = "48"

    denovo_name = "CNV48A"
    basis_names = ["CN1", "CN2", "CN9", "CN20", "CNV48B", "CNV48D"]
    weights = ["0.94%", "48.72%", "28.44%", "8.42%", "13.48%", "0%"]
    denovo_cols = ["MutationType", "CNV48A"]
    basis_cols = basis_names.copy()
    basis_cols.insert(0, "MutationType")
    nonzero_exposures = np.random.uniform(size=len(basis_names))
    result = sp.run_PlotDecomposition(
        denovo_mtx[denovo_cols],
        denovo_name,
        basis_mtx[basis_cols],
        basis_names,
        weights,
        nonzero_exposures,
        output_path,
        project,
        mtype,
    )
    # sp.run_PlotDecomposition(denovo_mtx, basis_names, weights, output_path, project, mtype, True, statistics, "COSMICv3-GRCh37", "This is where a custom message would go.")
    merger.append(result)

    for ind in range(5, 0, -1):
        basis_names = basis_names[:ind]
        weights = weights[:ind]
        denovo_cols = ["MutationType", "CNV48A"]
        basis_cols = basis_names.copy()
        basis_cols.insert(0, "MutationType")
        nonzero_exposures = np.random.uniform(size=len(basis_names))
        result = sp.run_PlotDecomposition(
            denovo_mtx[denovo_cols],
            denovo_name,
            basis_mtx[basis_cols],
            basis_names,
            weights,
            nonzero_exposures,
            output_path,
            project,
            mtype,
        )
        merger.append(result)

    merger.write(os.path.join(output_path, "Result_Decomposition_Plots_CNV48.pdf"))
    return time.time() - s


if __name__ == "__main__":
    print("Running SBS matrix input tests...")
    cosmic_fit_SBS_matrix_test()
    decompose_fit_SBS_matrix_test()
    denovo_fit_SBS_matrix_test()

    print("Running DBS matrix input tests...")
    cosmic_fit_DBS_matrix_test()
    decompose_fit_DBS_matrix_test()
    denovo_fit_DBS_matrix_test()

    print("Running ID matrix input tests...")
    cosmic_fit_ID_matrix_test()
    decompose_fit_ID_matrix_test()
    denovo_fit_ID_matrix_test()

    print("Running VCF input tests...")
    cosmic_fit_vcf_test()
    decompose_fit_vcf_test()
    denovo_fit_vcf_test()

    print("Running SV matrix input tests...")
    cosmic_fit_SV_matrix_test()
    decompose_fit_SV_matrix_test()
    denovo_fit_SV_matrix_test()

    print("Running CNV matrix input tests...")
    cosmic_fit_CNV_matrix_test()
    decompose_fit_CNV_matrix_test()
    denovo_fit_CNV_matrix_test()

    print("Decomposition Plot tests...")
    time_96 = gen_SBS96()
    print("Completed SBS96 plots in", time_96, "seconds.")
    time_288 = gen_SBS288()
    print("Completed SBS288 plots in", time_288, "seconds.")
    time_1536 = gen_SBS1536()
    print("Completed SBS1536 plots in", time_1536, "seconds.")
    time_78 = gen_DBS78()
    print("Completed DBS78 plots in", time_78, "seconds.")
    time_83 = gen_ID83()
    print("Completed ID83 plots in", time_83, "seconds.")
    time_48 = gen_CNV48()
    print("Completed CNV48 plots in", time_48, "seconds.")
