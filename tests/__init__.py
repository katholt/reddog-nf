import logging
import os
import pathlib
import sys


logging.getLogger().setLevel(51)
tests_directory = pathlib.Path(__file__).parent
# Must set add ./bin/ directory to PATH and PYTHONPATH to access scripts and enable imports
sys.path.append(os.path.abspath('bin/'))
os.environ['PATH'] += f':{os.path.abspath("bin/")}'
