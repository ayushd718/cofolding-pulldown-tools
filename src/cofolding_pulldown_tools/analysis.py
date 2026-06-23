import os
import csv
from pathlib import Path
import multiprocessing
import subprocess
from pandaprot import PandaProt
from ipsae import calculate_ipsae

def process_directory(dir):
    # generate the ipSAE and protein contact information for a single directory
    pass

def generate_csv(output_dir: os.PathLike):
    # go through all of the output dirs within the output_dir and use the process_directory helper function to calculate ipSAE and contact info, then combine them into a single csv or csv-like file to be used for downstream analysis by the user
    ## first generate ipSAEs, protein-protein interactions/contacts
    pass