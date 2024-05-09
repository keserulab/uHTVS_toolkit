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
from openbabel import openbabel, pybel
import re
from dimorphite_dl import DimorphiteDL
from os import path


# Getting the proper inputs

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_name',required=True,help='Input smi file, separated by space.')
parser.add_argument('-c', '--count',required=True,help='Input smi file, separated by space.')
parser.add_argument('-map', '--map_name',required=True,help='Name of the maps.fld file for AutoDockGPU')
parser.add_argument('-start', '--start',required=True,help='Start #')
parser.add_argument('-end', '--end',required=True,help='End #')
parser.add_argument('-n', '--n_cpu',required=True,help='Number of CPUs to use')
parser.add_argument('-sj', '--n_subjobs',required=True,help='Number of subjobs')
parser.add_argument('-prot', '--protonation',required=False,default='no',choices=['yes', 'no'],help='Find the major protonation form at pH: 7.4 (default: no)')
parser.add_argument('-loc', '--localopt',required=False,default='no',choices=['yes', 'no'],help='Perform local optimization on generated 3D structure (default: no)')
io_args = parser.parse_args()



# Function to return the nth line from the smile
def get_lines(fp, line_numbers):
    return (x for i, x in enumerate(fp) if i in line_numbers)


# Function to divide the 3D confgen and docking into subjobs
def multiproc(count):  

       # Calculating the starting and ending positions
       
       start = int(io_args.start) + subchunk_size * count	# The line number to start with is extracted from the other script
       if int(chunk_size) % int(io_args.n_subjobs) != 0 and count == (int(io_args.n_subjobs)-1):
        end = int(io_args.start) + int(chunk_size)	# The last line is based on the number of subjobs
       else:
        end = int(io_args.start) + subchunk_size * (count + 1)	# The last line is based on the number of subjobs
       fb = open(fb_file, 'a')
       IDfile = open(ID_file, 'a')              
              
       # Extracting smiles stings based on the line number in the .smi input file
              
       for ite in range(start, end): 	# Iterating through the line numbers between the limits calculated based on the total number of entries and the number of subjobs
        fp = open(io_args.input_name, 'r')
        lines = get_lines(fp, [ite])
        for line in lines:
         print(line)
         tr = line.strip().split(" ")
         if len(tr) != 2:
           print("Invalid number of columns for line: " + str(line))
         else:           
           smiles = tr[0]
           ID = tr[1]
           
           # Finding the correct protonation state of ligand with dimorphite_dl
           
           if str(io_args.protonation) == "yes":
            dimorphite_dl = DimorphiteDL(
              min_ph=7.4,
              max_ph=7.4,
              max_variants=10,
              label_states=False,
              pka_precision=1.0
            )             
            psmiles=dimorphite_dl.protonate(smiles)
            smiles = psmiles[0]
           
           # Converting the ligands into pdbqt and writing the batch files' entries
                
           try:
            name = "tmp"+str(ID)+".pdbqt"            
            fb.write("\n" + "./"+ name)
            fb.write("\n" + "Ligand_"+str(ID))
            IDfile.write(str(ID) + "\n")
            
            # Smiles conversion from smi to 3D pdbqt
            
            obConversion = openbabel.OBConversion()
            obConversion.SetInAndOutFormats("smi", "sdf")     
            mol = openbabel.OBMol()
            obConversion.ReadString(mol, smiles)
            pybelmol = pybel.Molecule(mol)                                    
            pybelmol.addh()	# Adding hydrogens
            pybelmol.make3D()	# Making 3D structure
            
            # Perform local optimization if required
            
            if str(io_args.localopt) == "yes":
             pybelmol.localopt()
            pybelmol.write("pdbqt", name, overwrite=True)
                         
           except:
            print("Error during conversion of ligand: "+str(ID))
            
        fp.close()	# Closes the input file                              
       fb.close()	# Closes the batch file
       IDfile.close()	# Closes the ID file
       

def main():        
 if __name__ == '__main__':
       
    # Run the multiple processes
    
    counter = range(int(io_args.n_subjobs))
    pool_size = cpus
    pool = mp.Pool(processes=pool_size,)
    collectedResults = pool.map(multiproc, counter)
    pool.close()
    
 
### Script starts here ### 
  
# Creating autodock batch file and a list file to delete ligand files later (numbered with the corresponding count number)   
    
fb_file = "Batch_py" + str(io_args.count) + ".txt"
ID_file = "Lig_IDs_" + str(io_args.count) + ".txt" 
with open (fb_file, "w") as batch:
 batch.write(io_args.map_name)
with open (ID_file, "w") as IDs:
 IDs.write("")

# Preliminary calculations

chunk_size = int(io_args.end) - int(io_args.start)	# Calculate how many ligands must be converted in this subjob
subchunk_size = int(chunk_size) // int(io_args.n_subjobs)	# Calculate how many ligands must be converted in one process
start_time = time.time()
cpus = int(io_args.n_cpu)
 
if int(io_args.n_cpu) > int(io_args.n_subjobs):
 print("ERROR: Number of subjobs is smaller than the number of cpus to use")
elif int(io_args.n_cpu) < 1:
 print("ERROR: invalid number of cpus")
else:
 main()
       
# Write out timing results

print("All jobs have been completed")        
print("%.0f procs Parallel - run time: %.2f seconds" % (cpus, time.time() - start_time))    
