#!/usr/bin/env python3

# Simple script to convert smiles strings to 3D sdf structures using RDkit. To use the script activate the rdkit (my-rdkit-env) environment with conda.
# Please note that the order of compounds in the sdf file won't be the same as in the .smi file.
# Dropped ligands will be collected at Dropped_ligands.txt
# To install openbabel: conda install -c conda-forge openbabel
# To install dimorphite_dl: python3 pip install dimorphite_dl and copy it to the proper place e.g. anaconda3/envs/my-rdkit-env/lib/python3.9/site-packages/
# To install multiprocessing: python3 pip install multiprocessing

import time
from itertools import (takewhile,repeat)
import multiprocessing as mp
import os
import argparse
import subprocess
import re
from os import path
import pandas as pd


# Getting the proper inputs

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_name',required=True,help='Input smi file, separated by space.')
parser.add_argument('-f', '--folder',required=True,help='Folder of the sub script, / is required at the end')
parser.add_argument('-o', '--output_name',required=True,help='Name of output file with scores')
parser.add_argument('-grid', '--grid_name',required=True,help='Name of the grid file (filename.zip) for Glide')
parser.add_argument('-nc', '--n_cpu',required=True,help='Number of CPUs to use per subjob')
parser.add_argument('-sj', '--n_subjobs',required=True,help='Number of subjobs')
parser.add_argument('-prot', '--protonation',required=False,default='no',choices=['yes', 'no'],help='Find the major protonation form at pH: 7.4 (default: no)')
parser.add_argument('-loc', '--localopt',required=False,default='no',choices=['yes', 'no'],help='Perform local optimization on generated 3D structure (default: no)')
parser.add_argument('-d', '--delete',required=False,default='no',choices=['yes', 'no'],help='Delete ligand pdbqt, dlg, xml and Batch files? (default: no)')
#parser.add_argument('-v', '--verbosity',required=False,default='scoresonly',choices=['scoresonly', 'verbose'],help='Set verbose if autodock log file is required (default: scoresonly)')
io_args = parser.parse_args()

# Function to calculate the number of entries in the file (the number of lines)
def rawincount(filename):
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )
    f.close()
    
# Function to divide the 3D confgen and docking into subjobs
def multiproc(count):        
       
       # Calculate the block sizes (number of smiles in one particular docking process) and send it to the subscript

       fb_file = "Ligand_file_" + str(count) + ".sdf"       
       num = rawincount(io_args.input_name)
       size = num//int(io_args.n_subjobs)
       countline = count * size
       
       # The last job contains the rest of the smiles if the compound number divided by the subjobs is not an integer
       
       if num % int(io_args.n_subjobs) != 0 and count == (int(io_args.n_subjobs)-1):
        countline2 = num  
       else:
        countline2 = (count + 1) * size
       
       # Run the ligand preparation python script for ligands between from countline to countline2 - 1
       
       cmdline_pre = "python " + str(io_args.folder) + "Glide_WorkFlow_smiles2score_sub.py -i " + str(io_args.input_name) + " -c " + str(count) + " -start " + str(countline) + " -end " + str(countline2) + " -n 1 -sj 1 -prot " + str(io_args.protonation) + " -loc " + str(io_args.localopt)	# Here the number of CPUs and the number of subjobs are set to 1, as Glide can handle ligand structures in a single sdf file. The subscripts are run in parallel in n copies defined by the number of CPUs in the main script. 
       
       print("\n\n-------------------------\nThe following command is running for ligand preparation #" + str(count) + ":\n" + str(cmdline_pre) + "\n-------------------------\n\n")	# For debugging purposes
       output_pre = subprocess.check_output(cmdline_pre, shell=True, universal_newlines=True)      

       # Preparing the input file of the docking
       
       glideinput = "Glide_docking_" + str(count) + ".in"
       glideoutput = "Glide_docking_" + str(count) + ".csv"
       if os.path.exists(glideinput):
        os.remove(glideinput) 
       gl_in = open(glideinput, 'a')
       
       text = "GRIDFILE   " + str(io_args.grid_name) + "\nLIGANDFILE   " + fb_file + "\nPRECISION   HTVS"
       gl_in.write(text)
       gl_in.close()

       # Performing the docking
       
       cmdline_glide = "\"${SCHRODINGER}/glide\" " + str(glideinput) + " -OVERWRITE -adjust -HOST localhost:1 -TMPLAUNCHDIR -WAIT"       
       print("\n\n-------------------------\nThe following command is running for Glide docking #" + str(count) + ":\n" + str(cmdline_glide) + "\n-------------------------\n\n")	# For debugging purposes
       out = subprocess.check_output(cmdline_glide, shell=True, universal_newlines=True)
       
       # Extracting docking scores to readable format
       
       usecols = ["title", "r_i_docking_score"]
       docking_data = pd.read_csv(glideoutput, usecols=usecols)
       score_file = "Glide_docking_scores_" + str(count) + ".txt"       
       docking_data.to_csv(score_file, sep=" ", header=False, index=False) 
       
       # File deletion if rquired
       
       if io_args.delete == "yes":
         ligandfilesdf = "Ligand_file_" + str(count) + ".sdf"
         preface = "Glide_docking_" + str(count)
         inputfile = str(preface) + ".in"         
         subjoblogfile = str(preface) + "_subjobs.log"
         posefile = str(preface) + "_subjob_poses.zip"
         tarfile = str(preface) + "_subjobs.tar.gz"
         skipfile = str(preface) + "_skip.csv"
         logfile = str(preface) + ".log"
         pvfile = str(preface) + "_pv.maegz"
         csvfile = str(preface) + ".csv" 
         if os.path.exists(ligandfilesdf):
          os.remove(ligandfilesdf) 
         if os.path.exists(inputfile):
          os.remove(inputfile) 
         if os.path.exists(subjoblogfile):
          os.remove(subjoblogfile)                         
         if os.path.exists(posefile):
          os.remove(posefile) 
         if os.path.exists(tarfile):
          os.remove(tarfile)        
         if os.path.exists(skipfile):
          os.remove(skipfile) 
         if os.path.exists(logfile):
          os.remove(logfile) 
         if os.path.exists(pvfile):
          os.remove(pvfile) 
         if os.path.exists(csvfile):
          os.remove(csvfile) 
          
       print("Job #" + str(count) + " is completed")       

    
def main():        
 if __name__ == '__main__':
    
    print("Docking workflow has been started.")
    
    # Creates a new output file which will store the docking results and a txt file containing the dropped ligands
    with open (io_args.output_name, "w") as scores: 
      scores.write("ZINC_ID r_i_docking_score\n")
   
    counter = range(int(io_args.n_subjobs))
    pool_size = cpus
    pool = mp.Pool(processes=pool_size,)
    collectedResults = pool.map(multiproc, counter)
    pool.close()

    # Extracting the docking scores into one file and cleaning up

    for i in range(int(io_args.n_subjobs)):
     score_file = "Glide_docking_scores_" + str(i) + ".txt"
     commandline = "cat " + str(score_file) + " >> " + str(io_args.output_name)
     #print(commandline)
     out = subprocess.check_output(commandline, shell=True, universal_newlines=True)
     os.remove(score_file)
    

### Script starts here ### 

# Preliminary calculations
 
start_time = time.time()
cpus = int(io_args.n_cpu)
main() 

# Write out timing results

print("All jobs have been completed")        
print("%.0f GPUs and %.0f CPUs Parallel - run time: %.2f seconds" % (cpus, cpus, time.time() - start_time))
