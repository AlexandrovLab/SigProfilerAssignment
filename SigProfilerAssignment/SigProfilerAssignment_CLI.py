#!/usr/bin/env python3

import sys
from SigProfilerAssignment.controllers import cli_controller


def main_function():
    commands = {
        "decompose_fit": "Perform decomposition fitting on the input samples.",
        "denovo_fit": "Perform de novo fitting on the input samples.",
        "cosmic_fit": "Perform COSMIC fitting on the input samples.",
    }

    if len(sys.argv) < 2 or sys.argv[1].lower() not in commands:
        print_usage(commands)
        sys.exit(1)

    command = sys.argv[1].lower()
    args = sys.argv[2:]

    controller = cli_controller.CliController()

    if command == "decompose_fit":
        controller.dispatch_decompose_fit(args)
    elif command == "denovo_fit":
        controller.dispatch_denovo_fit(args)
    elif command == "cosmic_fit":
        controller.dispatch_cosmic_fit(args)


def print_usage(commands):
    """Prints the usage message."""
    print("Usage: SigProfilerAssignment <command> [<args>]\n")
    print("Commands:")
    for cmd, desc in commands.items():
        print(f"  {cmd}: {desc}")


if __name__ == "__main__":
    main_function()
