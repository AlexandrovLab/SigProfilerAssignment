import argparse
from typing import List
from SigProfilerAssignment import decomposition as decomp


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def parse_arguments_common(args: List[str], description: str) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "samples",
        help="Path to the input somatic mutations file (if using segmentation file/mutational matrix) or input folder (mutation calling file/s).",
    )
    parser.add_argument("output", help="Path to the output folder.")
    parser.add_argument(
        "--signatures", help="Path to the signatures file.", default=None
    )
    parser.add_argument(
        "--signature_database",
        help="Path to the signature database file.",
        default=None,
    )
    parser.add_argument(
        "--nnls_add_penalty", type=float, default=0.05, help="NNLS add penalty."
    )
    parser.add_argument(
        "--nnls_remove_penalty", type=float, default=0.01, help="NNLS remove penalty."
    )
    parser.add_argument(
        "--initial_remove_penalty",
        type=float,
        default=0.05,
        help="Initial remove penalty.",
    )
    parser.add_argument(
        "--genome_build",
        default="GRCh37",
        help="The reference genome build (default: GRCh37). Supported genomes: {GRCh37, GRCh38, mm9, mm10, rn6}.",
    )
    parser.add_argument(
        "--cosmic_version",
        type=float,
        default=3.4,
        help="COSMIC version (default: 3.4). Valid options: {1, 2, 3, 3.1, 3.2, 3.3, 3.4}.",
    )
    parser.add_argument(
        "--make_plots",
        type=str2bool,
        nargs="?",
        const=True,
        default=True,
        help="Generate plots (default: True).",
    )
    parser.add_argument(
        "--collapse_to_SBS96",
        type=str2bool,
        nargs="?",
        const=True,
        default=True,
        help="Collapse to SBS96 (default: True).",
    )
    parser.add_argument(
        "--connected_sigs",
        type=str2bool,
        nargs="?",
        const=True,
        default=True,
        help="Connected signatures (default: True).",
    )
    parser.add_argument(
        "--verbose",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Verbose output (default: False).",
    )
    parser.add_argument(
        "--new_signature_thresh_hold",
        type=float,
        default=0.8,
        help="New signature threshold (default: 0.8).",
    )
    parser.add_argument(
        "--exclude_signature_subgroups",
        default=None,
        help="Remove specific signature subgroups.",
    )
    parser.add_argument(
        "--exome",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Use exome renormalized COSMIC signatures (default: False).",
    )
    parser.add_argument(
        "--input_type",
        default="matrix",
        help="""Three accepted input types: "vcf", "seg:TYPE", "matrix". Default: "matrix". The accepted callers for TYPE are the following {"ASCAT", "ASCAT_NGS", "SEQUENZA", "ABSOLUTE", "BATTENBERG", "FACETS", "PURPLE", "TCGA"}.""",
    )
    parser.add_argument(
        "--context_type",
        default="96",
        help="""Required context type if `input_type` is "vcf". Valid options: "96", "288", "1536", "DINUC", "ID". Default: "96".""",
    )
    parser.add_argument(
        "--export_probabilities",
        type=str2bool,
        nargs="?",
        const=True,
        default=True,
        help="Export probabilities matrix per mutational context (default: True).",
    )
    parser.add_argument(
        "--export_probabilities_per_mutation",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Export probabilities matrices per mutation (default: False).",
    )
    parser.add_argument(
        "--volume",
        default=None,
        help="User specified directory for saving/loading template files. Note: The environment variable SIGPROFILERASSIGNMENT_VOLUME takes precedence over this parameter.",
    )

    return parser.parse_args(args)


class CliController:
    def dispatch_decompose_fit(self, user_args: List[str]) -> None:
        parsed_args = parse_arguments_common(
            user_args, "Perform decomposition fitting on the input samples."
        )
        decomp.spa_analyze(
            samples=parsed_args.samples,
            output=parsed_args.output,
            signatures=parsed_args.signatures,
            signature_database=parsed_args.signature_database,
            nnls_add_penalty=parsed_args.nnls_add_penalty,
            nnls_remove_penalty=parsed_args.nnls_remove_penalty,
            initial_remove_penalty=parsed_args.initial_remove_penalty,
            genome_build=parsed_args.genome_build,
            cosmic_version=parsed_args.cosmic_version,
            make_plots=parsed_args.make_plots,
            collapse_to_SBS96=parsed_args.collapse_to_SBS96,
            connected_sigs=parsed_args.connected_sigs,
            verbose=parsed_args.verbose,
            decompose_fit_option=True,
            denovo_refit_option=False,
            cosmic_fit_option=False,
            new_signature_thresh_hold=parsed_args.new_signature_thresh_hold,
            exclude_signature_subgroups=parsed_args.exclude_signature_subgroups,
            exome=parsed_args.exome,
            input_type=parsed_args.input_type,
            context_type=parsed_args.context_type,
            export_probabilities=parsed_args.export_probabilities,
            export_probabilities_per_mutation=parsed_args.export_probabilities_per_mutation,
        )

    def dispatch_denovo_fit(self, user_args: List[str]) -> None:
        parsed_args = parse_arguments_common(
            user_args, "Perform de novo fitting on the input samples."
        )
        decomp.spa_analyze(
            samples=parsed_args.samples,
            output=parsed_args.output,
            signatures=parsed_args.signatures,
            signature_database=parsed_args.signature_database,
            nnls_add_penalty=parsed_args.nnls_add_penalty,
            nnls_remove_penalty=parsed_args.nnls_remove_penalty,
            initial_remove_penalty=parsed_args.initial_remove_penalty,
            genome_build=parsed_args.genome_build,
            cosmic_version=parsed_args.cosmic_version,
            new_signature_thresh_hold=parsed_args.new_signature_thresh_hold,
            make_plots=parsed_args.make_plots,
            collapse_to_SBS96=parsed_args.collapse_to_SBS96,
            connected_sigs=parsed_args.connected_sigs,
            verbose=parsed_args.verbose,
            decompose_fit_option=False,
            denovo_refit_option=True,
            cosmic_fit_option=False,
            exome=parsed_args.exome,
            input_type=parsed_args.input_type,
            context_type=parsed_args.context_type,
            export_probabilities=parsed_args.export_probabilities,
            export_probabilities_per_mutation=parsed_args.export_probabilities_per_mutation,
        )

    def dispatch_cosmic_fit(self, user_args: List[str]) -> None:
        parsed_args = parse_arguments_common(
            user_args, "Perform COSMIC fitting on the input samples."
        )
        decomp.spa_analyze(
            samples=parsed_args.samples,
            output=parsed_args.output,
            signatures=parsed_args.signatures,
            signature_database=parsed_args.signature_database,
            nnls_add_penalty=parsed_args.nnls_add_penalty,
            nnls_remove_penalty=parsed_args.nnls_remove_penalty,
            initial_remove_penalty=parsed_args.initial_remove_penalty,
            genome_build=parsed_args.genome_build,
            cosmic_version=parsed_args.cosmic_version,
            make_plots=parsed_args.make_plots,
            collapse_to_SBS96=parsed_args.collapse_to_SBS96,
            connected_sigs=parsed_args.connected_sigs,
            verbose=parsed_args.verbose,
            decompose_fit_option=False,
            denovo_refit_option=False,
            cosmic_fit_option=True,
            exclude_signature_subgroups=parsed_args.exclude_signature_subgroups,
            exome=parsed_args.exome,
            input_type=parsed_args.input_type,
            context_type=parsed_args.context_type,
            export_probabilities=parsed_args.export_probabilities,
            export_probabilities_per_mutation=parsed_args.export_probabilities_per_mutation,
            sample_reconstruction_plots=False,
        )
