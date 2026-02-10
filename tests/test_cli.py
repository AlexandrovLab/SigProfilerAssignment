import pytest
import argparse
from SigProfilerAssignment.controllers.cli_controller import (
    parse_arguments_common,
    str2bool,
    str2list
)


def test_default_values():
    args = parse_arguments_common(
        ["dummy_sample", "dummy_output"], "Test default values"
    )
    assert args.make_plots == True
    assert args.collapse_to_SBS96 == True
    assert args.connected_sigs == True
    assert args.verbose == False
    assert args.export_probabilities == True
    assert args.export_probabilities_per_mutation == False
    assert args.exome == False


def test_argument_parsing():
    args = parse_arguments_common(
        [
            "dummy_sample",
            "dummy_output",
            "--genome_build",
            "GRCh38",
            "--cosmic_version",
            "3.4",
            "--make_plots",
            "False",
            "--collapse_to_SBS96",
            "False",
            "--connected_sigs",
            "False",
            "--verbose",
            "True",
            "--export_probabilities",
            "False",
            "--export_probabilities_per_mutation",
            "True",
            "--exome",
            "True",
            "--sample_reconstruction_plots",
            "png",
            "--exclude_signature_subgroups",
            "MMR_deficiency_signatures,POL_deficiency_signatures,HR_deficiency_signatures,BER_deficiency_signatures,Chemotherapy_signatures,Immunosuppressants_signatures,Treatment_signatures,APOBEC_signatures,Tobacco_signatures,UV_signatures,AA_signatures,Colibactin_signatures,Artifact_signatures,Lymphoid_signatures"
        ],
        "Test argument parsing",
    )
    assert args.samples == "dummy_sample"
    assert args.output == "dummy_output"
    assert args.genome_build == "GRCh38"
    assert args.cosmic_version == 3.4
    assert args.make_plots == False
    assert args.collapse_to_SBS96 == False
    assert args.connected_sigs == False
    assert args.verbose == True
    assert args.export_probabilities == False
    assert args.export_probabilities_per_mutation == True
    assert args.exome == True
    assert args.sample_reconstruction_plots == "png"
    assert args.exclude_signature_subgroups == ['MMR_deficiency_signatures', 'POL_deficiency_signatures', 'HR_deficiency_signatures' , 'BER_deficiency_signatures', 'Chemotherapy_signatures', 'Immunosuppressants_signatures', 'Treatment_signatures', 'APOBEC_signatures', 'Tobacco_signatures', 'UV_signatures', 'AA_signatures', 'Colibactin_signatures', 'Artifact_signatures', 'Lymphoid_signatures']


def test_boolean_conversion():
    assert str2bool("yes") == True
    assert str2bool("true") == True
    assert str2bool("t") == True
    assert str2bool("y") == True
    assert str2bool("1") == True
    assert str2bool("no") == False
    assert str2bool("false") == False
    assert str2bool("f") == False
    assert str2bool("n") == False
    assert str2bool("0") == False
    with pytest.raises(argparse.ArgumentTypeError):
        str2bool("maybe")

def test_str2list():
    assert str2list("arg1,arg2,arg3") == ["arg1", "arg2","arg3"]
    assert str2list("arg_unique") == ["arg_unique"]
    assert str2list("wrong.sepparator") == ["wrong.sepparator"]

if __name__ == "__main__":
    pytest.main()
