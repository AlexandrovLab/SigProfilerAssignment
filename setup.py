from setuptools import setup
import shutil
import os

#remove the dist folder first if exists
if os.path.exists("dist"):
    shutil.rmtree("dist")

VERSION = '0.0.21'


def write_version_py(filename='SigProfilerAssignment/version.py'):
    # Copied from numpy setup.py
    cnt = """
# THIS FILE IS GENERATED FROM SigProfilerAssignment SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
Update = 'Resolve reference signature default value issue for decomposition plots.'

    
    """
    fh = open(filename, 'w')
    fh.write(cnt % {'version': VERSION,})
    fh.close()

with open('README.md') as f:
	long_description = f.read()

requirements=[
          'scipy>=1.6.3',
          'numpy>=1.21.2',
          'pandas>=1.2.4',
          'SigProfilerMatrixGenerator>=1.2.13', 
          'sigProfilerPlotting>=1.3.6', 
          'pillow>=9.1.1',
          'statsmodels>=0.9.0',
          'scikit-learn>=0.24.2',
          'psutil>=5.6.1',
          'reportlab>=3.5.42',
          'PyPDF2>=1.26.0',
          'alive_progress>=2.4.1',
          'PyPDf2>=1.28.4',
          'pdf2image>=1.16.0',
          'PyMuPDF>=1.21.0'
           ]

write_version_py() 
setup(name='SigProfilerAssignment',
      version=VERSION,
      description='Mutational signatures attribution and decomposition tool',
      long_description=long_description,
      long_description_content_type='text/markdown',  # This is important!	
      url="https://github.com/AlexandrovLab/SigProfilerAssignment.git",
      author='Raviteja Vangara',
      author_email='rvangara@health.ucsd.edu',
      license='UCSD',
      packages=['SigProfilerAssignment'],
      install_requires=requirements,
      include_package_data=True,      
      zip_safe=False)
