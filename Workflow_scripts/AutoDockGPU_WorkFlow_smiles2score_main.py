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


# Getting the proper inputs

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_name',required=True,help='Input smi file, separated by spaces')
parser.add_argument('-f', '--folder',required=True,help='Folder of the sub script')
parser.add_argument('-o', '--output_name',required=True,help='Name of output file with scores')
parser.add_argument('-map', '--map_name',required=True,help='Name of the maps.fld file for AutoDockGPU')
parser.add_argument('-gpudev', '--gpudevices',required=True,help='comma separated GPU devices to use')
parser.add_argument('-nc', '--n_cpu',required=True,help='Number of CPUs to use per subjob')
parser.add_argument('-sj', '--n_subjobs',required=True,help='Number of subjobs')
parser.add_argument('-prot', '--protonation',required=False,default='no',choices=['yes', 'no'],help='Find the major protonation form at pH: 7.4 (default: no)')
parser.add_argument('-loc', '--localopt',required=False,default='no',choices=['yes', 'no'],help='Perform local optimization on generated 3D structure (default: no)')
parser.add_argument('-d', '--delete',required=False,default='no',choices=['yes', 'no'],help='Delete ligand pdbqt, dlg, xml and Batch files? (default: no)')
parser.add_argument('-v', '--verbosity',required=False,default='scoresonly',choices=['scoresonly', 'verbose'],help='Set verbose if autodock log file is required (default: scoresonly)')
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

       fb_file = "Batch_py" + str(count) + ".txt"       
       num = rawincount(io_args.input_name)
       size = num//int(io_args.n_subjobs)
       countline = count * size
       
       # The last job contains the rest of the smiles if the compound number divided by the subjobs is not an integer
       
       if num % int(io_args.n_subjobs) != 0 and count == (int(io_args.n_subjobs)-1):
        countline2 = num  
       else:
        countline2 = (count + 1) * size
       
       # Run the ligand preparation python script for ligands between from countline to countline2 - 1
       
       cmdline_pre = "python " + str(io_args.folder) + "AutoDockGPU_WorkFlow_smiles2score_sub.py -i " + str(io_args.input_name) + " -c " + str(count) + " -map " + str(io_args.map_name) + " -start " + str(countline) + " -end " + str(countline2) + " -n " + str(io_args.n_cpu) + " -sj " + str(io_args.n_cpu) + " -prot " + str(io_args.protonation) + " -loc " + str(io_args.localopt)
       print("\n\n-------------------------\nThe following command is running for ligand preparation #" + str(count) + ":\n" + str(cmdline_pre) + "\n-------------------------\n\n")	# For debugging purposes
       output_pre = subprocess.check_output(cmdline_pre, shell=True, universal_newlines=True)      
       
       # Docking and docking score extraction
       try:
         with open (io_args.output_name, "a") as score:

            # GPU device selection
            
            device = count % devs
            devnum = device + 1
            command = "nvidia-smi | grep 'autodock' | awk '{print $2}'"
            devoutput = subprocess.check_output(command, shell=True, universal_newlines=True)
            new = devoutput.split("\n")
            new.pop()
            #print(devoutput)
            #print(new)
            #print(device)
            try:
             if str(device) in new:
              #print("It is used")
              for i in devlist:
               if str(i) not in new:
                #print(str(i) + " should be used")
                devnum = int(i) + 1
                break
             #else:
              #print("Not used") 
            except:
             print("Something's wrong")                 
            #print(devnum)
            
            # Scores only extraction
            
            if io_args.verbosity == "scoresonly":
             cmdline = "autodock_gpu_64wi -B " + fb_file + " -devnum " + str(devnum) + " | grep 'samples,\|Ligand file\|evaluations. Best' | grep -v '('"
             print("\n\n-------------------------\nThe following command is running for ligand docking job #" + str(count) + ":\n" + "autodock_gpu_64wi -B " + fb_file + " -devnum " + str(devnum) + "\n-------------------------\n\n")	# For debugging purposes
             output = subprocess.check_output(cmdline, shell=True, universal_newlines=True) 
             
             # Reformatting the output lines
             
             output = re.sub(r"pdbqt\s*[0-9]* samples, best energy\s*", " ", output)
             output = re.sub(r"pdbqt\s*[0-9]* evaluations. Best energy\s*", " ", output)
             output = re.sub(r"[ ]+Ligand file: \./tmp", "", output)
             output = re.sub(" kcal/mol.", "", output)
             output = re.sub("\. ", " ", output)

            # Verbose extraction
                                         
            elif io_args.verbosity == "verbose":
             #devnum = devlist[count % devs]
             cmdline = "autodock_gpu_64wi -B " + fb_file + " -devnum " + str(devnum) + " > DockLog" + str(count)
             print("\n\n-------------------------\nThe following command is running for ligand docking job #" + str(count) + ":\n" + "autodock_gpu_64wi -B " + fb_file + " -devnum " + str(devnum) + "\n-------------------------\n\n")	# For debugging purposes
             output = subprocess.check_output(cmdline, shell=True, universal_newlines=True) 
             
             # Reformatting the output lines
             
             cmdline = "grep 'samples,\|Ligand file\|evaluations. Best' DockLog" + str(count) + "| grep -v '('"
             output = subprocess.check_output(cmdline, shell=True, universal_newlines=True) 
             output = re.sub(r"pdbqt\s*[0-9]* samples, best energy\s*", " ", output)
             output = re.sub(r"pdbqt\s*[0-9]* evaluations. Best energy\s*", " ", output)
             output = re.sub(r"[ ]+Ligand file: \./tmp", "", output)
             output = re.sub(" kcal/mol.", "", output)
             output = re.sub("\. ", " ", output)                     

            score.write(str(output))
       
       except:
        print("Job number # " + str(count) + " resulted in zero poses") 


       # File deletion after the docking process is completed (if it is required)
       
       dropped = ""
       ID_file = "Lig_IDs_" + str(count) + ".txt"
       if str(io_args.delete) == "yes":
        with open(ID_file, "r") as file:
         file_list = file.readlines()
         for ID in file_list:
             ID=str(ID).replace('\n','')
          

             try:
              name = "tmp"+str(ID)+".pdbqt"
              xml = "Ligand_"+str(ID)+".xml"
              dlg = "Ligand_"+str(ID)+".dlg"
              os.remove(name)
              os.remove(xml)
              os.remove(dlg)
             except:
              #print("Ligand "+str(ID)+" has not been docked")
              dropped += "\n" + str(ID) 	# Collecting dropped ligands in case of file deletion is requested
        os.remove("Batch_py" + str(count) + ".txt")
        os.remove(ID_file)
       
       # Collecting dropped ligands if file deletion is not requested
       
       else:
        #dropped += "\n"
        with open(ID_file, "r") as file:
         file_list = file.readlines()
        for ID in file_list:
         name = "Ligand_"+str(ID)+".xml"
         if path.exists(name) == False:
          #print("Ligand "+str(ID)+" has not been docked")
          dropped += '\n' + str(ID)

       with open ("Dropped_ligands.txt", "a") as drop:
        dropped=str(dropped).replace('\n\n','\n')       
        drop.write(dropped)

       print("Job #" + str(count) + " is completed")
    
def main():        
 if __name__ == '__main__':
    
    print("Docking workflow has been started.")
    
    # Creates a new output file which will store the docking results and a tyt file containing the dropped ligands
    with open (io_args.output_name, "w") as scores: 
      scores.write("ZINC_ID r_i_docking_score\n")
    with open ("Dropped_ligands.txt", "w") as drop:
        drop.write("ID")
    
    counter = range(int(io_args.n_subjobs))
    pool_size = gpus
    pool = mp.Pool(processes=pool_size,)
    collectedResults = pool.map(multiproc, counter)
    pool.close()
    

### Script starts here ### 

# Preliminary calculations

devlist = str(io_args.gpudevices).split(",")	# For GPU device selection
devs = len(devlist)
  
start_time = time.time()
cpus = int(io_args.n_cpu) * int(devs)
gpus = int(devs)
if int(devs) > int(io_args.n_subjobs):
 print("ERROR: Number of subjobs is smaller than the number of GPUs to use")
elif int(devs) < 1:
 print("ERROR: invalid number of GPUs")
else:
 main() 

# Write out timing results

print("All jobs have been completed")        
print("%.0f GPUs and %.0f CPUs Parallel - run time: %.2f seconds" % (gpus, cpus, time.time() - start_time))
