# Python scripts to combine Glide HTVS and AutoDockGPU docking with Deep-Docking

The scripts found in the Workflow_scripts folder are useful python scripts that combine Glide and AutoDockGPU based docking with Deep-Docking, developed by Gentile et al. in a fully automated manner. (Gentile F, Agrawal V, Hsing M, Ton A-T, Ban F, Norinder U, et al. *Deep Docking: A Deep Learning Platform for Augmentation of Structure Based Drug Discovery.* ACS Cent Sci 2020:acscentsci.0c00229.) Deep docking (DD) is a deep learning-based tool developed to accelerate docking-based virtual screening. Using a docking program of choice, one can screen extensive chemical libraries like ZINC15 (containing > 1.3 billion molecules) 50 times faster than typical docking. The scripts required to perform the Deep-Docking wokflow are found into the developers github repository at https://github.com/jamesgleave/DD_protocol. Please copy these files into the DeepDocking_scripts folder or into the folder of your choice and copy the path in the Glide_DeepDocking_script.py or AutoDockGPU_DeepDocking_script.py python scripts' respective sections.

## Requirements
* conda
* rdkit
* tensorflow >= 1.14.0 (1.15 GPU version recommended. If you are using cuda11, please use [nvidia-tensorflow](https://developer.nvidia.com/blog/accelerating-tensorflow-on-a100-gpus/))
* pandas
* numpy
* keras
* matplotlib
* scikit-learn
* Glide or AutoDockGPU
* openbabel
* DimorphiteDL
* multiprocessing

## Content of the folders

Workflow_scripts contains the actual scripts that can be customized for the specific docking studies. If you intend to apply Glide please use the Glide_*.py files, for AutoDockGPU based docking please use the AutoDockGPU_*.py files.

Glide_example and AutoDockGPU_example folders contain previously set example workflows that carry out Deep-Docking on the compound library found in the Example_files/Ligands/Example.smi file. To perform these workflows run the Glide_DeepDocking_script.py or AutodockGPU_DeepDocking_script.py files in the corresponding folders.

Example_files contains the files required to run the example workflows. Glide and AutoDockGPU subfolders contain the grid files for the docking, while Ligands contains the smiles file with the ligand to be docked, split_smiles and morgan_fingerprints are created with the scripts described at the DeepDocking README. These subfolders contain the smiles and morgan fingerprints of the molecules to be used during the DeepDocking workflow.

DeepDocking_scripts is an empty folder, please download here the files of the DeepDocking protocol from the developers github repository. (https://github.com/jamesgleave/DD_protocol)

## Documentation of the DeepDocking, smiles2score_main and smiles2score_sub scripts

The *_DeepDocking_script.py files are input files for an automated DeepDocking workflow using Glide HTVS or AutoDockGPU as docking engines. The first 11 variables (dd_loc, project_name, morgan_loc, smiles_loc, num_cpus, sampled_num, project_folder_loc, project_folder_loc2, mol_percent_first, mol_percent_last, num_mol_exported) set the parameters of the Deep-Docking workflow. The next section controls the variables of the docking engines. (Glide: glide_workflow_script_loc, glide_grid, cpus, sj, proton, conf, remove | AutoDockGPU: autodock_workflow_script_loc, ad_grid, gpus, cpus, sj, proton, conf, remove, verbosity). The final section of the variables sets further parameters of Deep-Docking. (num_it, ttime, num_hyp, rec_val, starting_it). The number of iterations can be set by num_it. In case of an error the workflow can be restarted from the last iteration by setting starting_it to the serial number of the actual iteration. Project_folder_loc has a second instance due to the different input variable reading of the Deep-Docking scripts.

The *_smiles2score_main.py script handles the docking of the molecules found in the input smiles files created by the Deep-Docking scripts (training, test, validation sets). The script divides the set of molecules into smaller groups to enable parallelization by CPUs (Glide) or GPUs (AutoDockGPU). The given subsets of molecules are forwarded to the sub script, which creates 3D structures based on the smile strings of the molecules. Finally, the molecules are docked using the previously set grid file. Docking scores and ligand IDs are written into the *_labels.txt files. The Deep-Docking workflow uses these files to refine the models. (Note: both type of docking can be modified. Glide docking can be modified at "Performing the docking" section, while AutoDockGPU docking can be modified at "Score only extraction" and "Verbose extraction" sections)

The *_smiles2score_sub.py script creates 3D structures of the molecules described by smiles strings. The script enables the estimation of protonation states by DimorphiteDL (pH is set to 7.4 by default) and also enables local optimization on the 3D coordinates using RdKit.

## uHTVS toolkit usage

### 1. SMILES split and morgan fingerprint generation:

Before running the *_DeepDocking_script.py files the starting SMILES file - containing the complete set of molecules to be screened - must be splitted into evenly populated files. According to the developers of Deep-Docking this step is recommended "to facilitate other steps such as random sampling and inference, and place these files into a new folder." Next, the morgan fingerprints of the molecules distributed in the split SMILES files must be generated using the Morgan_fing.py script available from the developers of Deep-Docking.

#### Splitting:

```bash
split -d -l [number of molecules in a single split file] smiles.smi smile_all_ --additional-suffix=.txt
```

According to the developers: "Ideally the number of split files should be equal to the number of CPUs used for random sampling (phase 1, see below), but always larger than the number of GPUs used for inference (phase 3, see below)."

#### Morgan fingerprint generation

Morgan fingerprints can be then generated in the correct DD format using the `Morgan_fing.py` script (located in the `utilities` folder of the scripts at https://github.com/jamesgleave/DD_protocol):

```bash
python Morgan_fing.py -sfp path_smile_folder -fp path_to_morgan_folder -fn name_morgan_folder -tp num_cpus
```

which will create all the fingerprints and place them in `path_to_morgan_folder/name_morgan_folder`.

### 2. Setting up the Deep-Docking workfow file (parameters and variables)

### 3. Activate the RdKit env of conda

Activate your conda env that contain the requirements. (e.g. conda activate my-rdkit-env)

### 4. Running the *_DeepDocking_script.py script.

### 5. Final writeout

The final results including the IDs and VHL scores of the final molecules are written out into the directory containing the *_DeepDocking_script.py script as smiles.csv and id_score.csv

## Tips and troubleshoot:

Use the complete path descriptions in the *_DeepDocking_script.py scripts (e.g. use /home/user/uHTVS_toolkit/example/ instead of ./ or ../)

If the number of virtual hits reaches 0 during one of the iterations set the mol_percent_first and mol_percent_last variables a larger number between 0 and 1.

If the first compounds of the splitted SMILES files are omitted during morgan fingerpint creation delete the line "ref.readline()" in the Morgan_fing.py script.

In case of file not found errors we suggest to use the os.getcwd() method instead of using external variables in the Deep-Docking python scripts.

If the workflow has to be restarted, the library containing the current iteration (e.g. iteration_4) should be deleted before running the *_DeepDocking_script.py containg the starting_it=(actual iteration e.g. 4) line.
