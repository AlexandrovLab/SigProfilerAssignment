from setuptools import setup
import shutil
import os

#remove the dist folder first if exists
if os.path.exists("dist"):
    shutil.rmtree("dist")

VERSION = '0.0.17'

with open('README.md') as f:
	long_description = f.read()

def write_version_py(filename='SigProfilerAssignment/version.py'):
    # Copied from numpy setup.py
    cnt = """
# THIS FILE IS GENERATED FROM SigProfilerAssignment SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
Update = 'Fix TMB plot y-axis when using exome parameter'
    
    """
    fh = open(filename, 'w')
    fh.write(cnt % {'version': VERSION,})
    fh.close()
requirements=[
          'scipy>=1.6.3',
          'numpy>=1.21.2',
          'pandas>=1.2.4', 
          'SigProfilerExtractor>=1.1.14',
          'SigProfilerMatrixGenerator>=1.2.12', 
          'sigProfilerPlotting>=1.3.2', 
          'pillow',
          'statsmodels>=0.9.0',
          'scikit-learn>=0.24.2',
          'psutil>=5.6.1',
          'reportlab>=3.5.42',
          'PyPDF2>=1.26.0',
          'alive_progress'
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
