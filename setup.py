from setuptools import setup
import shutil
import os

# remove the dist folder first if exists
if os.path.exists("dist"):
    shutil.rmtree("dist")

VERSION = "0.1.9"


def write_version_py(filename="SigProfilerAssignment/version.py"):
    # Copied from numpy setup.py
    cnt = """
# THIS FILE IS GENERATED FROM SigProfilerAssignment SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
Update = 'v0.1.9: Replace PdfMerger with PdfReader and PdfWriter for pypdf 5.0.0'

    
    """
    fh = open(filename, "w")
    fh.write(
        cnt
        % {
            "version": VERSION,
        }
    )
    fh.close()


with open("README.md") as f:
    long_description = f.read()

requirements = [
    "scipy>=1.6.3",
    "numpy>=1.21.2,<2.0.0",
    "pandas>=1.2.4,<2.0.0",
    "SigProfilerMatrixGenerator>=1.2.28",
    "sigProfilerPlotting>=1.3.24",
    "reportlab>=3.5.42",
    "pypdf>=5.0.0",
    "alive_progress>=2.4.1",
    "PyMuPDF>=1.21.0",  # required for package "fitz"
]

write_version_py()
setup(
    name="SigProfilerAssignment",
    version=VERSION,
    description="Mutational signatures attribution and decomposition tool",
    long_description=long_description,
    long_description_content_type="text/markdown",  # This is important!
    url="https://github.com/AlexandrovLab/SigProfilerAssignment.git",
    author="Raviteja Vangara",
    author_email="rvangara@health.ucsd.edu",
    license="UCSD",
    packages=["SigProfilerAssignment"],
    install_requires=requirements,
    extras_require={
        "tests": [
            "pytest",
        ],
    },
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "SigProfilerAssignment=SigProfilerAssignment.SigProfilerAssignment_CLI:main_function",
        ],
    },
    zip_safe=False,
)
