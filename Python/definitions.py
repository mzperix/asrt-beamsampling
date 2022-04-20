import os

ROOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
PYTHON_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT_DIR,'Data')
CONFIG_DIR = DATA_DIR
OUTPUT_DIR = os.path.join(PYTHON_DIR,'Output')
CONDOR_DIR = os.path.join(ROOT_DIR,'Condor')
TEMPLATE_DIR = os.path.join(DATA_DIR, 'templates')
