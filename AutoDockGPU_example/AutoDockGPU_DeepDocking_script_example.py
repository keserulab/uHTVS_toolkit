#!/usr/bin/env python3

import time
import os
import subprocess
import pandas as pd

# Before running the script:
# The map.fld (grid) file for the protein has to be created by AutoGrid
# The smiles file containing the database has to be split and the morgan fingerprints have to be generated
# Use conda and activate the rdkit environment (default my-rdkit-env)



#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# The following section of inputs has to be set manually. The first series of inputs sets the parameters of the DeepDocking workflow, the second series of inputs handles the AutoDockGPU workflow

# The required variables for DD scripts

dd_loc = "../DeepDocking_scripts/"  								                        	# Location of the deepdocking python scripts
project_name = "AutoDockGPU_example"                                                        				# Name of the project folder (this has to be made manually)
morgan_loc = "../Example_files/Ligands/morgan_fingerprints/"                  									# Location of the morgan fingerprints
smiles_loc = "../Example_files/Ligands/split_smiles/"                  									# Location of the splitted smiles files (they are in .txt format)
num_cpus = 8                                                                                   		# Number of cpus to be used for the scripts:
sampled_num = 1000                                                                          		# Number of sampled molecules for each iterations (for the first iteration, this will be multiplied by 3)
project_folder_loc = "../"                                       						# Location of the folder that contains the folder named 'project_name'
project_folder_loc2 = ".."                                       						# Location of the folder that contains the folder named 'project_name'
mol_percent_first = 0.2                                                                        		# Percentage of the molecules labeled as hits at the first iterations (usually 0.02)
mol_percent_last = 0.02                                                                       		# Percentage of the molecules labeled as hits at the last iterations (usually 0.002)
num_mol_exported = 5000                                                                        		# Number of molecules to be exported at the end of the workflow

# The required variables for AutoDock Workflow:

autodock_workflow_script_loc = "../Workflow_scripts/"                             						# Location of the autodock workflow script
ad_grid = "../Example_files/AutoDockGPU/3HA8_empty_OK.maps.fld"                    									# Location of the AutoDock maps.fld file, including file name
gpus = '0,1,2,3'                                                                                    		# GPU-s to be used. Comma seperated GPU devices to use (use nvidia-smi for help)(e.g. 0,1,2)
cpus = 4                                                                                       		# Number of CPUs-to use per subjob, for OpenBabel (total CPU usage at a time is: number of gpus * this number) 
sj = 4                                                                                         		# Number of subjobs
proton = 'yes'                                                                                 		# Find the major protonation form at pH: 7.4? yes/no
conf = 'yes'                                                                                    		# Perform local optimization on generated 3D structure? yes/no
remove = 'yes'                                                                                  		# Delete ligand pdbqt, dlg, xml and Batch files? yes/no
verbosity = 'scoresonly'                                                                        			# Verbose or scoresonly. Set verbose if autodock log file is required. verbose/scoresonly
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# The following values are default values for DD according to the authors, but can be changed here if necessary:

num_it = 3                                                                         			# Number of iterations in integer (default is 11), now set to 3 for the example 
ttime = '00-04:00'                                                                   			# Time of set training (in minutes) (default is 00-04:00)
num_hyp = 12                                                                        			# Number of hyperparameters (default is 12)
rec_val = 0.9                                                                       			# Recall value (default is 0.9)
starting_it = 1                                                                     			# In case of an error that stops the workflow, the starting iteration can be set to another number (default is 1)

start_time = time.time()

# Commands:

txt1 = "python {}phase_1/molecular_file_count_updated.py -pt {} -it {} -cdd {} -t_pos {} -t_samp {}"
txt2 = "python {}phase_1/sampling.py -pt {} -fp {} -it {} -dd {} -t_pos {} -tr_sz {} -vl_sz {}"
txt3 = "python {}phase_1/sanity_check.py -pt {} -fp {} -it {}"
txt4 = "python {}phase_1/Extracting_morgan.py -pt {} -fp {} -it {} -md {} -t_pos {}"
txt5 = "python {}phase_1/Extracting_smiles.py -pt {} -fp {} -it {} -smd {} -t_pos {}"
txt6_test = "python {}AutoDockGPU_WorkFlow_smiles2score_main.py -i {}{}/iteration_{}/smile/test_smiles_final_updated.smi -f {} -o {}{}/iteration_{}/testing_labels_pre.txt -map {} -gpudev {} -nc {} -sj {} -prot {} -loc {} -d {} -v {}"
txt6_train = "python {}AutoDockGPU_WorkFlow_smiles2score_main.py -i {}{}/iteration_{}/smile/train_smiles_final_updated.smi -f {} -o {}{}/iteration_{}/training_labels_pre.txt -map {} -gpudev {} -nc {} -sj {} -prot {} -loc {} -d {} -v {}"
txt6_valid = "python {}AutoDockGPU_WorkFlow_smiles2score_main.py -i {}{}/iteration_{}/smile/valid_smiles_final_updated.smi -f {} -o {}{}/iteration_{}/validation_labels_pre.txt -map {} -gpudev {} -nc {} -sj {} -prot {} -loc {} -d {} -v {}"
txt7 = "python {}phase_2-3/simple_job_models.py -n_it {} -mdd {} -time {} -file_path {}{}/ -nhp {} -titr {} -n_mol {} --percent_first_mols {} -ct {} --percent_last_mols {}"
txt8 = "python -u {}phase_2-3/hyperparameter_result_evaluation.py -n_it {} --data_path {}{}/ -mdd {}"
txt9 = "python {}phase_2-3/simple_job_predictions.py --protein {} -file_path {}{} -n_it {} -mdd {}"
txt10 = "python {}final_phase/final_extraction.py -smile_dir {} -prediction_dir {}{}/iteration_{}/morgan_1024_predictions/ -processors {} -mols_to_dock {}"

# Variable to set from where to continue
cont = 0


print("DeepDocking workflow has been started with a total of {} iterations".format(num_it))

# Copying required DD scripts to the current working directory

cwd = os.getcwd()
copying1 = "cp {}phase_2-3/progressive_docking.py {}/"
out = subprocess.check_output(copying1.format(dd_loc, cwd), shell=True, universal_newlines=True)
copying2 = "cp -R {}phase_2-3/ML/ {}/"
out = subprocess.check_output(copying2.format(dd_loc, cwd), shell=True, universal_newlines=True)
copying3 = "cp {}phase_2-3/activation_script.sh {}/"
out = subprocess.check_output(copying3.format(dd_loc, cwd), shell=True, universal_newlines=True)
copying4 = "cp {}phase_2-3/Prediction_morgan_1024.py {}/"
out = subprocess.check_output(copying4.format(dd_loc, cwd), shell=True, universal_newlines=True)

# Starting the iterations

for i in range(starting_it, num_it+1):
    os.system('clear')
    print("Now doing the {}. iteration out of a total of {}".format(i, num_it))

    # Phase 1
    if cont == 0:
      if i == 1:  #for first iteration
          first_it_sampled_num = 3*sampled_num
          out = subprocess.check_output(txt1.format(dd_loc, project_name, i, morgan_loc, num_cpus, first_it_sampled_num), shell=True, universal_newlines=True)
          out = subprocess.check_output(txt2.format(dd_loc, project_name, project_folder_loc, i, morgan_loc, num_cpus, sampled_num, sampled_num), shell=True, universal_newlines=True)
      else:       #for other iterations
          new_cdd_loc = "{}{}/iteration_{}/morgan_1024_predictions/".format(project_folder_loc, project_name, i-1)
          out = subprocess.check_output(txt1.format(dd_loc, project_name, i, new_cdd_loc, num_cpus, sampled_num), shell=True, universal_newlines=True)
          out = subprocess.check_output(txt2.format(dd_loc, project_name, project_folder_loc, i, new_cdd_loc, num_cpus, sampled_num, sampled_num), shell=True, universal_newlines=True)
      out = subprocess.check_output(txt3.format(dd_loc, project_name, project_folder_loc, i), shell=True, universal_newlines=True)
      out = subprocess.check_output(txt4.format(dd_loc, project_name, project_folder_loc2, i, morgan_loc, num_cpus), shell=True, universal_newlines=True)
      out = subprocess.check_output(txt5.format(dd_loc, project_name, project_folder_loc, i, smiles_loc, num_cpus), shell=True, universal_newlines=True)
  
      # Changing tab to space in the .smi files
  
      with open('{}{}/iteration_{}/smile/test_smiles_final_updated.smi'.format(project_folder_loc, project_name, i), "r") as ofile1:
          data1 = ofile1.read()
          data1 = data1.replace('\t',' ')
      with open('{}{}/iteration_{}/smile/test_smiles_final_updated.smi'.format(project_folder_loc, project_name, i), "w") as ofile1:
          ofile1.write(data1)
      ofile1.close()
  
      with open('{}{}/iteration_{}/smile/train_smiles_final_updated.smi'.format(project_folder_loc, project_name, i), "r") as ofile2:
          data2 = ofile2.read()
          data2 = data2.replace('\t',' ')
      with open('{}{}/iteration_{}/smile/train_smiles_final_updated.smi'.format(project_folder_loc, project_name, i), "w") as ofile2:
          ofile2.write(data2)
      ofile2.close()
  
      with open('{}{}/iteration_{}/smile/valid_smiles_final_updated.smi'.format(project_folder_loc, project_name, i), "r") as ofile3:
          data3 = ofile3.read()
          data3 = data3.replace('\t',' ')
      with open('{}{}/iteration_{}/smile/valid_smiles_final_updated.smi'.format(project_folder_loc, project_name, i), "w") as ofile3:
          ofile3.write(data3)
      ofile3.close()
      
      # Usage of the AutoDockGPU workflow

      out = subprocess.check_output(txt6_test.format(autodock_workflow_script_loc, project_folder_loc, project_name,  i, autodock_workflow_script_loc, project_folder_loc, project_name, i, ad_grid, gpus, cpus, sj, proton, conf, remove, verbosity), shell=True, universal_newlines=True)
      out = subprocess.check_output(txt6_train.format(autodock_workflow_script_loc, project_folder_loc, project_name, i, autodock_workflow_script_loc, project_folder_loc, project_name, i, ad_grid, gpus, cpus, sj, proton, conf, remove, verbosity), shell=True, universal_newlines=True)
      out = subprocess.check_output(txt6_valid.format(autodock_workflow_script_loc, project_folder_loc, project_name, i, autodock_workflow_script_loc, project_folder_loc, project_name, i, ad_grid, gpus, cpus, sj, proton, conf, remove, verbosity), shell=True, universal_newlines=True)
      
      # Reformatting the label .txt files
  
      col1 = []
      col2 = []
      a = 0

      rewrite_test = open("{}{}/iteration_{}/testing_labels_pre.txt".format(project_folder_loc, project_name, i),"r")
      new_label_test = open("{}{}/iteration_{}/testing_labels.txt".format(project_folder_loc, project_name, i),"w")
  
      for row in rewrite_test: 
          if row == "ZINC_ID r_i_docking_score\n":
              pass
          else:
              a += 1
              rowsplit1 = row.replace("\n", "")
              rowsplit2 = rowsplit1.split(" ")
              rowsplit2.reverse()
              col1.append(float(rowsplit2[0]))
              col2.append(rowsplit2[1])
      df = pd.DataFrame({'col1':col1, 'col2':col2})
      dfsorted = df.sort_values(by='col1')
      dfsorted = dfsorted.reset_index(drop=True)
      new_label_test.write("r_i_docking_score,ZINC_ID\n") 
      for row2 in range(a):
          row_list = dfsorted.loc[row2, :].values.flatten().tolist()
          newrow = str(row_list[0]) + "," + str(row_list[1]) + "\n"
          new_label_test.write(newrow) 
      col1.clear()
      col2.clear()
      a = 0
  
      rewrite_test.close()
      os.remove("{}{}/iteration_{}/testing_labels_pre.txt".format(project_folder_loc, project_name, i))
      new_label_test.close()
  
      rewrite_train = open("{}{}/iteration_{}/training_labels_pre.txt".format(project_folder_loc, project_name, i),"r")
      new_label_train = open("{}{}/iteration_{}/training_labels.txt".format(project_folder_loc, project_name, i),"w")
  
      for row in rewrite_train:
          if row == "ZINC_ID r_i_docking_score\n":
              pass
          else:
              a += 1
              rowsplit1 = row.replace("\n", "")
              rowsplit2 = rowsplit1.split(" ")
              rowsplit2.reverse()
              col1.append(float(rowsplit2[0]))
              col2.append(rowsplit2[1])
      df = pd.DataFrame({'col1':col1, 'col2':col2})
      dfsorted = df.sort_values(by='col1')
      dfsorted = dfsorted.reset_index(drop=True)
      new_label_train.write("r_i_docking_score,ZINC_ID\n")
      for row2 in range(a):
          row_list = dfsorted.loc[row2, :].values.flatten().tolist()
          newrow = str(row_list[0]) + "," + str(row_list[1]) + "\n"
          new_label_train.write(newrow)
      col1.clear()
      col2.clear()
      a = 0
  
      rewrite_train.close()
      os.remove("{}{}/iteration_{}/training_labels_pre.txt".format(project_folder_loc, project_name, i))
      new_label_train.close()
  
      rewrite_valid = open("{}{}/iteration_{}/validation_labels_pre.txt".format(project_folder_loc, project_name, i),"r")
      new_label_valid = open("{}{}/iteration_{}/validation_labels.txt".format(project_folder_loc, project_name, i),"w")
  
      for row in rewrite_valid:
          if row == "ZINC_ID r_i_docking_score\n":
              pass
          else:
              a += 1
              rowsplit1 = row.replace("\n", "")
              rowsplit2 = rowsplit1.split(" ")
              rowsplit2.reverse()
              col1.append(float(rowsplit2[0]))
              col2.append(rowsplit2[1])
      df = pd.DataFrame({'col1':col1, 'col2':col2})
      dfsorted = df.sort_values(by='col1')
      dfsorted = dfsorted.reset_index(drop=True)
      new_label_valid.write("r_i_docking_score,ZINC_ID\n")
      for row2 in range(a):
          row_list = dfsorted.loc[row2, :].values.flatten().tolist()
          newrow = str(row_list[0]) + "," + str(row_list[1]) + "\n"
          new_label_valid.write(newrow)
      col1.clear()
      col2.clear()
      a = 0
  
      rewrite_valid.close()
      os.remove("{}{}/iteration_{}/validation_labels_pre.txt".format(project_folder_loc, project_name, i))
      new_label_valid.close()
      
      # Phase 2
      
      out = subprocess.check_output(txt7.format(dd_loc, i, morgan_loc, ttime, project_folder_loc, project_name, num_hyp, num_it, sampled_num, mol_percent_first, rec_val, mol_percent_last), shell=True, universal_newlines=True)
      simple_job_directory = "{}{}/iteration_{}/simple_job/".format(project_folder_loc, project_name, i)
      for j in os.listdir(simple_job_directory):
          simple_job_filepath = os.path.join(simple_job_directory, j)
          if os.path.isfile(simple_job_filepath):
              out = subprocess.check_output("chmod a+x {}".format(simple_job_filepath), shell=True, universal_newlines=True)
              out = subprocess.check_output(simple_job_filepath, shell=True, universal_newlines=True)
  
      # Phase 3
              
      out = subprocess.check_output(txt8.format(dd_loc, i, project_folder_loc, project_name, morgan_loc), shell=True, universal_newlines=True)
      out = subprocess.check_output(txt9.format(dd_loc, project_name, project_folder_loc, project_name, i, morgan_loc), shell=True, universal_newlines=True)
      simple_job_predictions_directory = "{}{}/iteration_{}/simple_job_predictions/".format(project_folder_loc, project_name, i)
      for k in os.listdir(simple_job_predictions_directory):
          simple_job_predictions_filepath = os.path.join(simple_job_predictions_directory, k)
          if os.path.isfile(simple_job_predictions_filepath):
              out = subprocess.check_output("chmod a+x {}".format(simple_job_predictions_filepath), shell=True, universal_newlines=True)
              out = subprocess.check_output(simple_job_predictions_filepath, shell=True, universal_newlines=True)
  
      # Last iteration: exporting predicted molecules
              
      if i == num_it:
          #out = subprocess.check_output("mkdir {}{}/predicted_molecules".format(project_folder_loc, project_name), shell=True, universal_newlines=True)
          #out = subprocess.check_output("cd {}{}/predicted_molecules/".format(project_folder_loc, project_name), shell=True, universal_newlines=True)
          out = subprocess.check_output(txt10.format(dd_loc, smiles_loc, project_folder_loc, project_name, i, num_cpus, num_mol_exported), shell=True, universal_newlines=True)
          os.system('clear')
          print("DeepDocking workflow finished. Output .csv files (id_score.csv and smiles.csv) for {} molecules are found in the current folder .".format(num_mol_exported))
          print("The id_score.csv file contains the molecule name, and a score from 0-1 corresponding to how good the molecule is according to the machine learning model. Larger number means better molecule.")
          print("The smiles.csv file contains the SMILES string for these molecules, and their name.")
          
print("The DeepDocking workflow has been completed")        
print("Run time: %.2f seconds" % (time.time() - start_time)) 
